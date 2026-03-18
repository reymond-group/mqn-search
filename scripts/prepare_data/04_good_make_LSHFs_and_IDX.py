#!/usr/bin/env python3

from __future__ import annotations

import argparse
import concurrent.futures as cf
import glob
import logging
import os
from pathlib import Path
import re
import sys

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
SRC_DIR = REPO_ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


DATASETS_DIR = REPO_ROOT / "datasets"
DEFAULT_INPUT_DIR = DATASETS_DIR / "BINS"
DEFAULT_LSHF_DIR = DATASETS_DIR / "LSHF"
DEFAULT_IDX_DIR = DATASETS_DIR / "IDX"
DEFAULT_LOG_DIR = DATASETS_DIR / "LOGS"

FOREST_D = 512
FOREST_L = 32
CHUNK_SIZE = 300_000
POOL_WORKERS = None

LSH_PREFIX = "100M_ZINC_ChEMBL_BIN-"
BIN_PATTERN = re.compile(r"BIN-(\d{4})")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build per-BIN LSH forests and IDX parquet files from BIN CSV inputs.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Directory containing per-BIN CSV files from 2b_make_bins.py.",
    )
    parser.add_argument(
        "--lshf-dir",
        type=Path,
        default=DEFAULT_LSHF_DIR,
        help="Directory for .lsh outputs.",
    )
    parser.add_argument(
        "--idx-dir",
        type=Path,
        default=DEFAULT_IDX_DIR,
        help="Directory for index parquet outputs.",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help="Directory for logs.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=CHUNK_SIZE,
        help="CSV chunk size used while scanning each BIN file.",
    )
    parser.add_argument(
        "--pool-workers",
        type=int,
        default=POOL_WORKERS,
        help="Worker count for fingerprint computation. Default uses os.cpu_count().",
    )
    return parser.parse_args()


def setup_logging(log_dir: Path) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=log_dir / "04_good_make_LSHFs_and_IDX.log",
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def log_and_print(message: str) -> None:
    print(message)
    logging.info(message)


def bin_from_basename(basename: str) -> str:
    match = BIN_PATTERN.search(basename)
    if not match:
        raise ValueError(f"Cannot parse BIN from '{basename}'")
    return match.group(1)


def lsh_out_path_for(basename: str, lshf_dir: Path) -> Path:
    return lshf_dir / f"{basename}-{FOREST_D}-{FOREST_L}.lsh"


def idx_out_path_for(bin4: str, idx_dir: Path) -> Path:
    return idx_dir / f"index_smiles_BIN-{bin4}.parquet"


def calc_fp(smiles: str) -> list[int] | None:
    from mqn_search.mhfp_512 import calc_mhfp

    try:
        fp = calc_mhfp(smiles)
        if fp is None:
            return None
        fp_list = list(fp)
        if len(fp_list) != FOREST_D:
            return None
        return fp_list
    except Exception:
        return None


def process_csv_file(csv_file: Path, lshf_dir: Path, idx_dir: Path, chunk_size: int, pool_workers: int | None) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq
    import tmap as tm
    from tmap import VectorUint

    basename = csv_file.stem
    bin4 = bin_from_basename(basename)
    out_lsh = lsh_out_path_for(basename, lshf_dir)
    out_idx = idx_out_path_for(bin4, idx_dir)

    if out_lsh.exists() and out_idx.exists():
        log_and_print(f"Skipping {csv_file.name}: outputs already exist.")
        return

    all_smiles: list[str] = []
    all_fps: list[list[int]] = []

    chunk_no = 0
    for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
        chunk_no += 1
        if "SMILES" not in chunk.columns:
            raise ValueError(f"'SMILES' column not found in chunk {chunk_no} of {csv_file}")

        chunk = chunk.dropna(subset=["SMILES"])
        smiles_list = chunk["SMILES"].astype(str).tolist()
        log_and_print(f"{basename} | chunk {chunk_no}: {len(smiles_list)} SMILES")

        with cf.ProcessPoolExecutor(max_workers=pool_workers) as executor:
            results = list(executor.map(calc_fp, smiles_list, chunksize=100))

        kept = 0
        for smiles, fp in zip(smiles_list, results):
            if fp is None:
                continue
            all_smiles.append(smiles)
            all_fps.append(fp)
            kept += 1

        log_and_print(f"{basename} | chunk {chunk_no}: kept {kept} valid fingerprints; total {len(all_fps)}")

    if not all_fps:
        log_and_print(f"Skipping {csv_file.name}: no valid fingerprints found.")
        return

    log_and_print(f"{basename} | building LSHForest on {len(all_fps)} fingerprints")
    tm_fps = [VectorUint(fp) for fp in all_fps]
    lsh_forest = tm.LSHForest(d=FOREST_D, l=FOREST_L)
    lsh_forest.batch_add(tm_fps)
    lsh_forest.index()
    lsh_forest.store(str(out_lsh))

    idx_col = list(range(len(all_smiles)))
    table = pa.Table.from_arrays(
        [
            pa.array(idx_col, type=pa.int32()),
            pa.array(all_smiles, type=pa.string()),
            pa.array(all_fps, type=pa.list_(pa.uint32())),
        ],
        names=["idx", "smiles", "fp"],
    )
    pq.write_table(table, out_idx, compression="zstd", use_dictionary=False)
    log_and_print(f"{basename} | saved {out_lsh.name} and {out_idx.name}")


def main() -> int:
    args = parse_args()
    setup_logging(args.log_dir)
    args.lshf_dir.mkdir(parents=True, exist_ok=True)
    args.idx_dir.mkdir(parents=True, exist_ok=True)

    from mqn_search.mhfp_512 import check_env

    if not check_env():
        print("The environment does not match the reference MHFP. Refusing to continue.")
        return 1

    csv_files = sorted(Path(p) for p in glob.glob(str(args.input_dir / "*.csv")))
    if not csv_files:
        log_and_print(f"No CSV files found in {args.input_dir}")
        return 1

    log_and_print(f"Found {len(csv_files)} BIN CSV file(s) in {args.input_dir}")
    for csv_file in csv_files:
        try:
            process_csv_file(csv_file, args.lshf_dir, args.idx_dir, args.chunk_size, args.pool_workers)
        except Exception as exc:
            logging.exception("Failed processing %s: %s", csv_file, exc)
            print(f"Failed processing {csv_file}: {exc}")

    log_and_print("Finished building LSHF and IDX outputs.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
