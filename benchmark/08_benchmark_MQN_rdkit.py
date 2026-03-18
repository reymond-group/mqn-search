#!/usr/bin/env python3
"""
Compute MQN fingerprints for a molecule file and write them to CSV.

Example:
    python compute_mqn.py \
        --input compounds.smi \
        --output mqn.csv \
        --log process.log \
        --format smi \
        --chunk 100000
"""

import argparse
import csv
import logging
import time
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def get_supplier(path: str, fmt: str):
    if fmt == "smi":
        # Assumes one SMILES string per line, optional name after a space
        return Chem.SmilesMolSupplier(path, delimiter=" ", titleLine=False, nameColumn=-1)
    if fmt == "sdf":
        return Chem.SDMolSupplier(path)
    raise ValueError(f"Unsupported input format: {fmt}")


def configure_logging(log_file: Path):
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )


def write_header(writer):
    header = ["index"] + [f"MQN{i + 1}" for i in range(42)]
    writer.writerow(header)


def main():
    parser = argparse.ArgumentParser(description="Compute RDKit MQN fingerprints.")
    parser.add_argument("--input", required=True, help="Path to input .smi or .sdf file.")
    parser.add_argument("--output", required=True, help="Destination CSV file.")
    parser.add_argument("--log", required=True, help="Log file for timing and progress.")
    parser.add_argument(
        "--format",
        choices=["smi", "sdf"],
        default="smi",
        help="Input file format. Default is smi.",
    )
    parser.add_argument(
        "--chunk",
        type=int,
        default=100000,
        help="Number of rows to buffer before each disk write. Adjust for your I/O.",
    )
    args = parser.parse_args()

    configure_logging(Path(args.log))
    start = time.perf_counter()
    logging.info("Computation started using RDKit MQN")

    supplier = get_supplier(args.input, args.format)

    processed = 0
    with open(args.output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        write_header(writer)

        buffer = []
        for idx, mol in enumerate(supplier):
            if mol is None:
                continue  # RDKit returns None for parse failures
            fp = rdMolDescriptors.MQNs_(mol)
            buffer.append([idx] + list(fp))
            processed += 1

            if len(buffer) >= args.chunk:
                writer.writerows(buffer)
                buffer.clear()
                logging.info(f"Processed {processed} molecules")

        if buffer:
            writer.writerows(buffer)

    elapsed = time.perf_counter() - start
    logging.info(f"Completed {processed} molecules in {elapsed:.2f} seconds")


if __name__ == "__main__":
    main()
