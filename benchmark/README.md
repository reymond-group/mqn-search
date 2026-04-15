## Benchmark

This directory contains evaluation and timing scripts for the MQN/MHFP search workflow.

Files:

- `08_benchmark_MQN_rdkit.py`
  Computes RDKit MQNs for a supplied molecule file and logs runtime.
- `09_benchmark_MQN_binary.sh`
  Benchmarks the external C++ MQN binary used in the prep pipeline.
- `16_benchmark_lshforest_robust.py`
  Evaluates approximate per-bin LSH retrieval against exact MHFP truth files from `datasets/VS_MHFP/`.

The main benchmark entry point is `16_benchmark_lshforest_robust.py`.

By default it expects:

- query drug metadata in `datasets/145-ChEMBL34_top-sellers_23_FINAL.csv`
- per-bin LSH forests in `datasets/LSHF/`
- per-bin IDX parquet files in `datasets/IDX/`
- exact truth files in `datasets/VS_MHFP/`

Benchmark outputs are written under `benchmark_out/`.

Minimal benchmark run:

```bash
conda activate mqn-search
python benchmark/16_benchmark_lshforest_robust.py
```

Useful variants:

```bash
python benchmark/16_benchmark_lshforest_robust.py --only-drug APIXABAN
python benchmark/16_benchmark_lshforest_robust.py --resume
python benchmark/16_benchmark_lshforest_robust.py --span 1 --per-bin-candidates-k 5000
```

The `09_benchmark_MQN_binary.sh` helper benchmarks the external C++ MQN binary directly:

```bash
bash benchmark/09_benchmark_MQN_binary.sh
```
