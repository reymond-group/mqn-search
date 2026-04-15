## MQN Search

MQN search, data-preparation, and benchmarking code for a ZINC/ChEMBL search pipeline built around:

- MQN-based binning
- MHFP fingerprints
- per-bin `tmap` LSH forests
- exact-vs-approximate evaluation scripts

## Repository Layout

- [scripts/prepare_data](scripts/prepare_data): sequential preprocessing and index-building pipeline
- [benchmark](benchmark): benchmark and evaluation scripts
- [src/mqn_search](src/mqn_search): importable Python code
- [src_cpp/mqn_search](src_cpp/mqn_search): C++ MQN encoder sources
- [datasets](datasets): local input data only; derived outputs are generated, not committed

## Environment

The repository is packaged via [pyproject.toml](pyproject.toml) and the working environment is described in [environment.yml](environment.yml).

Typical setup:

```bash
conda env create -f environment.yml
conda activate mqn-search
```

Minimal validation:

```bash
pytest -q tests/test_imports.py tests/test_mhfp.py
python scripts/prepare_data/3b_csv2mhfp-parquet_good-MHFP.py --help
```

If the MHFP tests pass and `3b` starts cleanly, the environment is typically good enough for the active Python pipeline.

## Build Requirement

The preprocessing pipeline expects a compiled MQN executable at `build/bin/mqn_new`.

This repository contains the C++ source under [src_cpp/mqn_search](src_cpp/mqn_search), but it does not currently ship a checked-in build system such as `CMakeLists.txt` or a `Makefile`.

The binary has been built successfully with:

```bash
mkdir -p build/bin
cd src_cpp/mqn_search
g++ -std=c++17 -Wall -Wextra -O2 *.cpp -o ../../build/bin/mqn_new
```

The first prep script will fail early if the binary is missing.

## Data

The tracked `datasets/` contents are gzipped source inputs. See [datasets/README.md](datasets/README.md).

Generated outputs such as `BINS/`, `LSHF/`, `IDX/`, `VS_MHFP/`, parquet chunks, and logs are intended to be created locally and kept out of version control.

## Preparation Pipeline

The active prep sequence is:

1. [1_add_MQN2file.sh](scripts/prepare_data/1_add_MQN2file.sh)
   Adds MQN output from the compiled C++ binary to the SMILES input.
2. [2_extract_BIN.py](scripts/prepare_data/2_extract_BIN.py)
   Extracts the leading MQN sum into `BIN`.
3. [3_remove_drugs.py](scripts/prepare_data/3_remove_drugs.py)
   Replaces overlapping entries with the curated ChEMBL drug rows.
4. [2b_make_bins.py](scripts/prepare_data/2b_make_bins.py)
   Splits the merged CSV into per-bin CSV files.
5. [3b_csv2mhfp-parquet_good-MHFP.py](scripts/prepare_data/3b_csv2mhfp-parquet_good-MHFP.py)
   Adds MHFP to the merged CSV and writes parquet chunks.
6. [04_good_make_LSHFs_and_IDX.py](scripts/prepare_data/04_good_make_LSHFs_and_IDX.py)
   Builds one LSH forest and one IDX parquet per bin from `datasets/BINS/`.
7. [5_mhfpVS.py](scripts/prepare_data/5_mhfpVS.py)
   Runs exact MHFP virtual screening over the parquet chunks to create truth files for benchmarking.

See [scripts/prepare_data/README.md](scripts/prepare_data/README.md) for brief per-script notes.

## Benchmarks

Current benchmark entry points:

- [08_benchmark_MQN_rdkit.py](benchmark/08_benchmark_MQN_rdkit.py): RDKit MQN generation benchmark
- [09_benchmark_MQN_binary.sh](benchmark/09_benchmark_MQN_binary.sh): external C++ MQN binary benchmark
- [16_benchmark_lshforest_robust.py](benchmark/16_benchmark_lshforest_robust.py): approximate search benchmark against exact MHFP truth files

See [benchmark/README.md](benchmark/README.md) for details.
