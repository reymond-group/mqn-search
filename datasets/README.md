## Datasets

This directory only tracks small, gzipped source inputs needed to reproduce the local pipeline.

Current tracked inputs:

- `ZINC_ChEMBL34_hi_conf_1M.smi.gz`
- `145-ChEMBL34_top-sellers_23_FINAL.csv.gz`

Before running the preparation pipeline, unpack the files locally:

```bash
gunzip datasets/ZINC_ChEMBL34_hi_conf_1M.smi.gz
gunzip datasets/145-ChEMBL34_top-sellers_23_FINAL.csv.gz
```

Generated outputs such as:

- `BINS/`
- `LSHF/`
- `IDX/`
- `VS_MHFP/`
- parquet chunks
- logs

are local working artifacts and should not be committed.
