#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import logging
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import sys
import time

import pyarrow as pa
import pyarrow.parquet as pq

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
SRC_DIR = REPO_ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from mqn_search.mhfp_512 import calc_mhfp, check_env


### To be used with the correct environment made from troubleshooting/ubern.yml
if not check_env():
    print(f'The test of environment did not match the reference mhfp. Refusing to continue. Please check your environment.')
    sys.exit()
else:   
    print('The Python envorionment is consistent. Will contiue...')

### adds MHFP column
### SMILES,ID,MQN,BIN
### Adds an additional column MHFP and saves in parquet files by chunks

chunk_size = 1_000_00 #0 # odjust for the size of library
fp_length = 512  # for output filenames only


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add MHFP fingerprints to the merged MQN CSV and write chunked parquet files.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=REPO_ROOT / "datasets" / "ZINC_ChEMBL34_hi_conf_1M_145-drugs_MQN_BIN.csv",
        help="Input CSV produced by 3_remove_drugs.py.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=REPO_ROOT / "datasets" / f"ZINC_ChEMBL34_hi_conf_1M_MQN_MHFP{fp_length}",
        help="Prefix for parquet chunk files and logs.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=chunk_size,
        help="Number of CSV rows to process per parquet chunk.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=os.cpu_count() or 1,
        help="Number of worker processes to use. Default uses os.cpu_count().",
    )
    return parser.parse_args()

def setup_logging(output_prefix):
    log_file = str(output_prefix) + '_new_MHFP-' + str(fp_length) + '_LAST_PART.log'
    # If log file already exists, rename it to .OLD
    if os.path.exists(log_file):
        old_log_file = log_file + '.OLD'
        shutil.copy(log_file, old_log_file)
        os.remove(log_file)
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(message)s')

def log_and_print(message):
    print(message)
    logging.info(message)

# Function to process a chunk of data
def process_chunk(smiles_list):
    log_and_print(f"Processing chunk. {len(smiles_list)} SMILES received.")
    # Create a new MHFPEncoder instance for each chunk
    #### we don't use mhfp lib mhfp_encoder = MHFPEncoder(fp_length)
    # Apply the function to the specified column
    ### mhfp_list = [calc_mhfp(smiles, mhfp_encoder) for smiles in smiles_list]
    mhfp_list = [calc_mhfp(smiles) for smiles in smiles_list] #[tmap.VectorUint]
    #return pa.array(mhfp_list)
        # Convert tmap.VectorUint explicitly to list[int] here
    mhfp_list_list = [list(fp) if fp is not None else [] for fp in mhfp_list]
    return pa.array(mhfp_list_list)

# Function to read, process, and save a chunk (pyarrow reading of chunks is buggy)
def read_process_save_chunk(start, chunk_size, input_file, output_file):
    rows = []
    
    # Open CSV file and read the desired chunk
    with open(input_file, 'r') as f:
        reader = csv.reader(f)
        # Extract header
        header = next(reader)
        if 'SMILES' not in header:
            raise ValueError("The CSV file does not contain a 'SMILES' column.")
        # Skip rows up to the start
        for _ in range(start):
            next(reader, None)
        # Read the chunk
        for _ in range(chunk_size):
            try:
                rows.append(next(reader))
            except StopIteration:
                break

    if not rows:
        log_and_print(f"No rows to process for start {start}. Skipping this chunk.")
        return
    
    # Convert rows to PyArrow table
    try:
        # Transpose rows into columns using zip, then create a table
        columns = list(zip(*rows))
        original_table = pa.table(columns, names=header)
    except Exception as e:
        log_and_print(f"Error creating table from rows for start {start}. Error: {e}")
        return

    # Extract SMILES column for processing
    smiles_list = original_table['SMILES'].to_pylist()

    # Process the chunk to calculate MHFP
    mhfp_array = process_chunk(smiles_list)

    # Append the new column to the original table
    new_table = original_table.append_column('MHFP', mhfp_array)
    
    # Save the processed chunk to Parquet
    pq.write_table(new_table, output_file, compression='gzip')

def worker(start, chunk_size, input_file, output_prefix):
    output_file = f"{output_prefix}_chunk_{start}.parquet"
    
    if os.path.exists(output_file):
        log_and_print(f"Output file {output_file} already exists. Skipping.")
        return
    
    try:
        log_and_print(f"Worker started for start row {start}")
        read_process_save_chunk(start, chunk_size, input_file, output_file)
        log_and_print(f"File {output_file} written.")
    except Exception as e:
        log_and_print(f"Error in worker for start row {start}. Error: {e}")

# Main function to manage multiprocessing
def main(input_file, output_prefix, chunk_size, workers):
    # we want to follow what is going on
    setup_logging(output_prefix)
    start_time = time.time()
    
    # Get total number of rows in the CSV file
    with open(input_file, 'r') as f:
        total_rows = sum(1 for _ in f) - 1  # Subtract 1 for header row
        elapsed_time = time.time() - start_time
        log_and_print(f"Full dataset analyzed in {elapsed_time:.2f} seconds: Will process {total_rows} SMILES.")

    start_time = time.time()

    # Set up multiprocessing
    cpu_count = max(1, int(workers))
    pool = mp.Pool(cpu_count)
    log_and_print(f"Using {cpu_count} cores.")
    starts = range(0, total_rows, chunk_size)  # Header is read separetely, these are rows after header
    print(f"these are the starts: {list(starts)}")

    # Run the pool of workers
    pool.starmap(worker, [(start, chunk_size, input_file, output_prefix) for start in starts])
    pool.close()
    pool.join()

if __name__ == "__main__":
    args = parse_args()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    main(str(args.input), str(args.output_prefix), args.chunk_size, args.workers)
