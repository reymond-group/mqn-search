#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[1]

    parser = argparse.ArgumentParser(
        description="Split the merged MQN CSV into per-BIN CSV files for downstream indexing.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=repo_root / "datasets" / "ZINC_ChEMBL34_hi_conf_1M_145-drugs_MQN_BIN.csv",
        help="Input CSV produced by 3_remove_drugs.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "datasets" / "BINS",
        help="Directory where per-BIN CSV files will be written.",
    )
    return parser.parse_args()


def validate_columns(df: pd.DataFrame) -> None:
    required = {"SMILES", "ID", "BIN"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {sorted(missing)}")


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    validate_columns(df)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    df = df.loc[:, ["SMILES", "ID", "BIN"]].copy()

    written = 0
    for bin_value, group in df.groupby("BIN", sort=True):
        file_name = f"100M_ZINC_ChEMBL_BIN-{int(bin_value):04d}.csv"
        group.to_csv(args.output_dir / file_name, index=False)
        written += 1

    print(f"Read {len(df)} rows from {args.input}")
    print(f"Wrote {written} BIN file(s) under {args.output_dir}")


if __name__ == "__main__":
    main()
