## Prepare Data

This directory contains the sequential data-preparation pipeline.

Before running the pipeline:

1. unpack the gzipped source files in [datasets](/scratch/work/projects/AbbVie/mqn/github/repos/mqn-search/datasets)
2. activate the `mqn-search` conda environment
3. ensure `build/bin/mqn_new` exists and is executable

Typical sequence:

```bash
gunzip datasets/ZINC_ChEMBL34_hi_conf_1M.smi.gz
gunzip datasets/145-ChEMBL34_top-sellers_23_FINAL.csv.gz

conda activate mqn-search

bash scripts/prepare_data/1_add_MQN2file.sh
python scripts/prepare_data/2_extract_BIN.py
python scripts/prepare_data/3_remove_drugs.py
python scripts/prepare_data/2b_make_bins.py
python scripts/prepare_data/3b_csv2mhfp-parquet_good-MHFP.py
python scripts/prepare_data/04_good_make_LSHFs_and_IDX.py
python scripts/prepare_data/5_mhfpVS.py
```

Active scripts:

1. `1_add_MQN2file.sh`
   Runs the compiled MQN binary on the raw SMILES input and writes a CSV with `SMILES`, `ID`, and `MQN`.
2. `2_extract_BIN.py`
   Extracts the MQN sum into a `BIN` column.
3. `3_remove_drugs.py`
   Replaces overlapping rows with the curated query-drug entries.
4. `2b_make_bins.py`
   Splits the merged CSV into per-bin CSV files under `datasets/BINS/`.
5. `3b_csv2mhfp-parquet_good-MHFP.py`
   Recomputes MHFP for the merged CSV and writes parquet chunks.
6. `04_good_make_LSHFs_and_IDX.py`
   Builds one LSH forest and one IDX parquet per bin.
7. `5_mhfpVS.py`
   Performs exact MHFP virtual screening to generate truth files in `datasets/VS_MHFP/`.

The scripts are intended to be run locally against unpacked files in `datasets/`.

Expected generated outputs:

- `datasets/ZINC_ChEMBL34_hi_conf_1M_MQN.csv`
- `datasets/ZINC_ChEMBL34_hi_conf_1M_MQN_BIN.csv`
- `datasets/ZINC_ChEMBL34_hi_conf_1M_145-drugs_MQN_BIN.csv`
- `datasets/BINS/`
- `datasets/ZINC_ChEMBL34_hi_conf_1M_MQN_MHFP512_chunk_*.parquet`
- `datasets/LSHF/`
- `datasets/IDX/`
- `datasets/VS_MHFP/`
- `datasets/LOGS/`
