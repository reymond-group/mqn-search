## mqn_search Python Package

This package contains the importable Python code used by the preparation and benchmark scripts.

Current modules:

- [__init__.py](__init__.py): package entry point
- [mhfp_512.py](mhfp_512.py): MHFP calculation and environment consistency check

The active scripts import [mhfp_512.py](mhfp_512.py) to:

- compute MHFP fingerprints from SMILES
- verify that the Python environment produces the expected reference MHFP

Minimal smoke test from the repo root:

```bash
pytest -q tests/test_imports.py tests/test_mhfp.py
```
