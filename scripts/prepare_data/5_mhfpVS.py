#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import glob
import heapq
import logging
from pathlib import Path
import re
import sys
import time

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
SRC_DIR = REPO_ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


DATASETS_DIR = REPO_ROOT / "datasets"
DEFAULT_QUERY_CSV = DATASETS_DIR / "145-ChEMBL34_top-sellers_23_FINAL.csv"
DEFAULT_PARQUET_GLOB = str(DATASETS_DIR / "ZINC_ChEMBL34_hi_conf_1M_MQN_MHFP512_chunk_*.parquet")
DEFAULT_OUTPUT_DIR = DATASETS_DIR / "VS_MHFP"
DEFAULT_LOG_DIR = DATASETS_DIR / "LOGS"
FP_LENGTH = 512


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run MHFP virtual screening for each query and write per-query result CSVs.",
    )
    parser.add_argument(
        "--query-csv",
        type=Path,
        default=DEFAULT_QUERY_CSV,
        help="Query CSV with SMILES, ID, SUM_BIN, MQN, MHFP, DRUG_NAME.",
    )
    parser.add_argument(
        "--parquet-glob",
        default=DEFAULT_PARQUET_GLOB,
        help="Glob matching parquet chunks produced by 3b_csv2mhfp-parquet_good-MHFP.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for per-query nearest-neighbour CSV outputs.",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help="Directory for logs and timing CSV output.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=1000,
        help="Number of closest matches to keep per query.",
    )
    return parser.parse_args()


def setup_logging(log_dir: Path) -> Path:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"5_mhfpVS_{FP_LENGTH}.log"
    if log_file.exists():
        backup = log_file.with_suffix(log_file.suffix + ".OLD")
        backup.write_text(log_file.read_text())
        log_file.unlink()

    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
    )
    return log_file


def log_and_print(message: str) -> None:
    print(message)
    logging.info(message)


def norm_name(value: str) -> str:
    text = value.strip().upper()
    text = re.sub(r"[^A-Z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "UNKNOWN"


def load_queries(input_csv: Path) -> pd.DataFrame:
    required = ["SMILES", "ID", "SUM_BIN", "DRUG_NAME"]
    query_df = pd.read_csv(input_csv, usecols=required)
    return query_df


def get_parquet_files(parquet_glob: str) -> list[Path]:
    return sorted(Path(path) for path in glob.glob(parquet_glob))


def calculate_mhfp(smiles: str) -> np.ndarray | None:
    from mqn_search.mhfp_512 import calc_mhfp

    fp = calc_mhfp(smiles)
    if fp is None:
        return None
    fp_list = list(fp)
    if len(fp_list) != FP_LENGTH:
        return None
    return np.asarray(fp_list, dtype=np.uint32)


def scan_table_matches(
    qsig: np.ndarray,
    table,
    closest_matches: list[tuple[float, int, dict[str, object]]],
    num_closest: int,
    unique_counter: int,
) -> int:
    mhfp_column = table["MHFP"].combine_chunks()
    values = mhfp_column.values
    offsets = mhfp_column.offsets
    vals = values.to_numpy(zero_copy_only=False)
    offs = offsets.to_numpy()
    lens = offs[1:] - offs[:-1]
    if len(lens) == 0:
        return unique_counter
    if not (lens == FP_LENGTH).all():
        raise RuntimeError(f"Unexpected MHFP length(s): {np.unique(lens)}")

    fps = vals.reshape(len(lens), FP_LENGTH)
    sims = (fps == qsig).mean(axis=1)

    for j, sim in enumerate(sims):
        row_data = table.slice(j, 1).to_pydict()
        item = (float(sim), unique_counter, row_data)
        unique_counter += 1
        if len(closest_matches) < num_closest:
            heapq.heappush(closest_matches, item)
        elif item[0] > closest_matches[0][0]:
            heapq.heapreplace(closest_matches, item)

    return unique_counter


def find_closest_matches(
    qsig: np.ndarray,
    parquet_files: list[Path],
    num_closest: int,
) -> tuple[list[tuple[float, dict[str, object]]], float]:
    import pyarrow.parquet as pq

    closest_matches: list[tuple[float, int, dict[str, object]]] = []
    unique_counter = 0
    search_start = time.perf_counter()

    for file in parquet_files:
        with pq.ParquetFile(file) as parquet_file:
            for row_group_idx in range(parquet_file.num_row_groups):
                table = parquet_file.read_row_group(row_group_idx)
                unique_counter = scan_table_matches(
                    qsig=qsig,
                    table=table,
                    closest_matches=closest_matches,
                    num_closest=num_closest,
                    unique_counter=unique_counter,
                )

    elapsed = time.perf_counter() - search_start
    matches = [(sim, row_data) for sim, _, row_data in closest_matches]
    matches.sort(key=lambda item: item[0], reverse=True)
    return matches, elapsed


def save_results(
    drug_name: str,
    bin_value: int,
    query_id: str,
    closest_matches: list[tuple[float, dict[str, object]]],
    output_dir: Path,
) -> Path:
    def flatten_cell(value: object) -> object:
        if isinstance(value, list) and len(value) == 1:
            return value[0]
        return value

    safe_drug = norm_name(drug_name)
    output_file = output_dir / f"{safe_drug}_BIN-{int(bin_value):04d}_{query_id}_MHFP-{FP_LENGTH}.csv"

    with output_file.open("w", newline="") as file:
        writer = csv.writer(file)
        if closest_matches:
            first_row = closest_matches[0][1]
            header = ["similarity"] + list(first_row.keys())
            writer.writerow(header)
            for sim, row_data in closest_matches:
                writer.writerow([sim] + [flatten_cell(v) for v in row_data.values()])
        else:
            writer.writerow(["similarity"])

    return output_file


def init_timings_csv(path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    with path.open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "drug_name",
                "query_id",
                "sum_bin",
                "compute_mhfp_s",
                "search_s",
                "write_s",
                "total_query_s",
                "hits_written",
                "output_file",
            ]
        )


def append_timing(
    path: Path,
    drug_name: str,
    query_id: str,
    sum_bin: int,
    compute_mhfp_s: float,
    search_s: float,
    write_s: float,
    total_query_s: float,
    hits_written: int,
    output_file: Path,
) -> None:
    with path.open("a", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                drug_name,
                query_id,
                sum_bin,
                f"{compute_mhfp_s:.6f}",
                f"{search_s:.6f}",
                f"{write_s:.6f}",
                f"{total_query_s:.6f}",
                hits_written,
                str(output_file),
            ]
        )


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_file = setup_logging(args.log_dir)
    timings_csv = args.log_dir / "5_mhfpVS_timings.csv"
    init_timings_csv(timings_csv)

    from mqn_search.mhfp_512 import check_env

    if not check_env():
        print("The environment does not match the reference MHFP. Refusing to continue.")
        return 1

    query_df = load_queries(args.query_csv)
    parquet_files = get_parquet_files(args.parquet_glob)
    if not parquet_files:
        log_and_print(f"No parquet files found for glob: {args.parquet_glob}")
        return 1

    log_and_print(f"Query file loaded: {args.query_csv} ({query_df.shape[0]} rows)")
    log_and_print(f"Dataset found: {len(parquet_files)} parquet file(s)")
    log_and_print(f"Timing log: {timings_csv}")
    log_and_print(f"Text log: {log_file}")

    overall_start = time.perf_counter()
    for _, row in query_df.iterrows():
        smiles = row["SMILES"]
        query_id = str(row["ID"])
        sum_bin = int(row["SUM_BIN"])
        drug_name = str(row["DRUG_NAME"])
        query_start = time.perf_counter()

        t0 = time.perf_counter()
        qsig = calculate_mhfp(smiles)
        compute_mhfp_s = time.perf_counter() - t0
        if qsig is None:
            log_and_print(f"Skipping query_id={query_id} drug={drug_name}: MHFP failed")
            continue

        log_and_print(f"Searching query_id={query_id} drug={drug_name} bin={sum_bin}")
        closest_matches, search_s = find_closest_matches(qsig, parquet_files, args.top_k)

        t1 = time.perf_counter()
        output_file = save_results(
            drug_name=drug_name,
            bin_value=sum_bin,
            query_id=query_id,
            closest_matches=closest_matches,
            output_dir=args.output_dir,
        )
        write_s = time.perf_counter() - t1
        total_query_s = time.perf_counter() - query_start

        append_timing(
            path=timings_csv,
            drug_name=drug_name,
            query_id=query_id,
            sum_bin=sum_bin,
            compute_mhfp_s=compute_mhfp_s,
            search_s=search_s,
            write_s=write_s,
            total_query_s=total_query_s,
            hits_written=len(closest_matches),
            output_file=output_file,
        )
        log_and_print(
            f"Finished query_id={query_id} drug={drug_name} "
            f"search_s={search_s:.3f} total_s={total_query_s:.3f} hits={len(closest_matches)}"
        )

    overall_s = time.perf_counter() - overall_start
    log_and_print(f"Completed MHFP VS in {overall_s:.2f} seconds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
