#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LSHForest benchmark.

This version is hardened for long-running batch benchmarks:
- coherent CLI
- resolved output paths printed at startup
- resumable via status.csv
- per-drug error breadcrumbs
- optional per-bin debug logging
- robust single-worker mode using subprocess-per-bin isolation
- retry per bin up to N times before skipping that bin
- exact rerank timings split into FP loading vs scoring/heap maintenance
"""

from __future__ import annotations

import argparse
import ast
import csv
import dataclasses
import json
import os
import pathlib as _p
import re
import subprocess
import sys
import time
from collections import OrderedDict
from typing import Any, Dict, List, Optional, Sequence, Tuple

import concurrent.futures as cf
import heapq

import numpy as np
import pyarrow as pa
import pyarrow.compute as pac
import pyarrow.parquet as pq
from rdkit import Chem
import tmap as tm
from tmap import VectorUint

REPO_ROOT = _p.Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from mqn_search.mhfp_512 import calc_mhfp as _calc_mhfp

try:
    from mqn_search.mhfp_512 import check_env as _check_env
except Exception:
    _check_env = None


def _env_ok() -> bool:
    if callable(_check_env):
        try:
            return bool(_check_env())
        except Exception:
            return False
    return True


@dataclasses.dataclass(frozen=True)
class Paths:
    drugs145_csv: _p.Path = REPO_ROOT / "datasets" / "145-ChEMBL34_top-sellers_23_FINAL.csv"
    vs_drugs_dir: _p.Path = REPO_ROOT / "datasets" / "VS_MHFP"
    lsh_root: _p.Path = REPO_ROOT / "datasets" / "LSHF"
    mqn_bin: _p.Path = REPO_ROOT / "build" / "bin" / "mqn_new"
    idx_dir: _p.Path = REPO_ROOT / "datasets" / "IDX"
    manifests_root: _p.Path = REPO_ROOT / "datasets" / "manifests"

    out_dir: _p.Path = REPO_ROOT / "benchmark_out" / "default"
    per_drug_dir: _p.Path = REPO_ROOT / "benchmark_out" / "default" / "per_drug"
    approx_ids_csv: _p.Path = REPO_ROOT / "benchmark_out" / "default" / "approx_ids.csv"
    timings_csv: _p.Path = REPO_ROOT / "benchmark_out" / "default" / "timings.csv"
    metrics_csv: _p.Path = REPO_ROOT / "benchmark_out" / "default" / "metrics.csv"
    queries_csv: _p.Path = REPO_ROOT / "benchmark_out" / "default" / "queries.csv"
    status_csv: _p.Path = REPO_ROOT / "benchmark_out" / "default" / "status.csv"


@dataclasses.dataclass(frozen=True)
class RunCfg:
    k: int = 50
    per_bin_candidates_k: int = 5000
    kc: int = -1
    span: int = 0
    lsh_workers: int = max(1, (os.cpu_count() or 4) // 2)
    idx_cache_bins: int = 0
    rerank_per_bin_topk: int = 250
    bin_retry_max: int = 3


LSH_PREFIX = "100M_ZINC_ChEMBL_BIN-"
LSH_SUFFIX = "-512-32.lsh"

_IDX_ARROW_CACHE: "OrderedDict[int, pa.Table]" = OrderedDict()
_IDX_CACHE_MAX_BINS: int = 0
_WORKER_LSH_ROOT: Optional[_p.Path] = None
_WORKER_CACHE: "OrderedDict[int, tm.LSHForest]" = OrderedDict()
_WORKER_CACHE_BYTES: int = 0
_WORKER_CACHE_MAX_BYTES: int = 0


def pad4(n: int) -> str:
    return f"{n:04d}"


def ensure_dirs(p: Paths) -> None:
    p.out_dir.mkdir(parents=True, exist_ok=True)
    p.per_drug_dir.mkdir(parents=True, exist_ok=True)


def norm_drug_name(s: str) -> str:
    s = s.strip().upper()
    s = re.sub(r"[^A-Z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def _write_headers_if_missing(path: _p.Path, header: List[str]) -> None:
    if not path.exists() or path.stat().st_size == 0:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", newline="") as f:
            csv.writer(f).writerow(header)


def _write_error_file(path: _p.Path, msg: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(msg.rstrip() + "\n")


def _append_status(
    status_csv: _p.Path,
    drug_name: str,
    sum_bin: int,
    span: int,
    per_bin: int,
    mode: str,
    status: str,
    message: str = "",
) -> None:
    with status_csv.open("a", newline="") as f:
        csv.writer(f).writerow([
            drug_name,
            sum_bin,
            span,
            per_bin,
            mode,
            int(time.time()),
            status,
            message,
        ])


def _load_completed_drugs(
    status_csv: _p.Path,
    span: int,
    per_bin: int,
    mode: str,
) -> set[str]:
    done: set[str] = set()
    if not status_csv.exists():
        return done
    with status_csv.open() as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            try:
                if (
                    int(row.get("span", "-999")) == int(span)
                    and int(row.get("per_bin", "-999")) == int(per_bin)
                    and str(row.get("mode", "")) == str(mode)
                    and str(row.get("status", "")) == "DONE"
                ):
                    dn = str(row.get("drug_name", "")).strip()
                    if dn:
                        done.add(dn)
            except Exception:
                continue
    return done


def lsh_path_for_bin(bin_id: int, root: _p.Path) -> _p.Path:
    return root / f"{LSH_PREFIX}{pad4(bin_id)}{LSH_SUFFIX}"


def mqn_bins_for(sum_bin: int, span: int = 0) -> List[int]:
    if span <= 0:
        return [sum_bin]
    out = [sum_bin]
    for d in range(1, span + 1):
        lo = sum_bin - d
        hi = sum_bin + d
        if lo >= 0:
            out.append(lo)
        out.append(hi)
    seen = set()
    centered: List[int] = []
    for x in out:
        if x >= 0 and x not in seen:
            seen.add(x)
            centered.append(x)
    return centered


_MANIFEST_RE = re.compile(rf"^{re.escape(LSH_PREFIX)}(\d{{4}}){re.escape(LSH_SUFFIX)}$")


def _manifest_path(paths: Paths, drug_name: str, span: int) -> _p.Path:
    dnorm = norm_drug_name(drug_name).lower()
    return paths.manifests_root / f"span-{span}" / f"manifest_{dnorm}_span-{span}.txt"


def bins_from_manifest(manifest_file: _p.Path) -> List[int]:
    bins: List[int] = []
    if not manifest_file.exists():
        return bins
    with manifest_file.open() as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            m = _MANIFEST_RE.match(s)
            if not m:
                continue
            bins.append(int(m.group(1)))
    seen = set()
    out: List[int] = []
    for b in bins:
        if b not in seen:
            seen.add(b)
            out.append(b)
    return out


def filter_existing_lsh_bins(bins: Sequence[int], lsh_root: _p.Path) -> List[int]:
    ok: List[int] = []
    for b in bins:
        p = lsh_path_for_bin(b, lsh_root)
        if p.exists():
            ok.append(b)
        else:
            print(f"[WARN] drop bin {b}: missing {p}")
    return ok


def bins_for_drug(
    paths: Paths,
    drug_name: str,
    sum_bin: int,
    span: int,
    require_manifest: bool = False,
) -> Tuple[List[int], str]:
    mf = _manifest_path(paths, drug_name, span)
    mbins = bins_from_manifest(mf)
    if mbins:
        return filter_existing_lsh_bins(mbins, paths.lsh_root), "manifest"
    if require_manifest:
        return [], "manifest_missing"
    return filter_existing_lsh_bins(mqn_bins_for(sum_bin, span), paths.lsh_root), "computed"


def mqn_sum_stdio(smiles: str, binary_path: _p.Path) -> Optional[int]:
    try:
        proc = subprocess.run(
            [str(binary_path)],
            input=smiles,
            text=True,
            capture_output=True,
            check=True,
        )
        out = proc.stdout.strip().splitlines()[-1]
        lst = ast.literal_eval(out)
        if not (isinstance(lst, list) and len(lst) == 43):
            return None
        return int(lst[0])
    except Exception:
        return None


def make_query_vector_from_smiles(smiles: str) -> VectorUint:
    fp = list(_calc_mhfp(smiles))
    if len(fp) != 512:
        raise ValueError(f"calc_mhfp produced length={len(fp)} (expected 512)")
    return VectorUint([int(x) for x in fp])


def _load_idx_arrow(bin_id: int, idx_dir: _p.Path) -> pa.Table:
    global _IDX_ARROW_CACHE
    if _IDX_CACHE_MAX_BINS > 0:
        tbl = _IDX_ARROW_CACHE.get(bin_id)
        if tbl is not None:
            _IDX_ARROW_CACHE.move_to_end(bin_id)
            return tbl
    p_parq = idx_dir / f"index_smiles_BIN-{pad4(bin_id)}.parquet"
    if not p_parq.exists():
        raise FileNotFoundError(f"No Parquet IDX for bin {pad4(bin_id)} at {p_parq}")
    tbl = pq.read_table(p_parq, columns=["smiles", "fp"], memory_map=True)
    if _IDX_CACHE_MAX_BINS > 0:
        _IDX_ARROW_CACHE[bin_id] = tbl
        _IDX_ARROW_CACHE.move_to_end(bin_id)
        while len(_IDX_ARROW_CACHE) > _IDX_CACHE_MAX_BINS:
            _, old_tbl = _IDX_ARROW_CACHE.popitem(last=False)
            del old_tbl
    return tbl


def _bin_file_size_bytes(bin_id: int, root: _p.Path) -> int:
    try:
        return lsh_path_for_bin(bin_id, root).stat().st_size
    except Exception:
        return 0


def _init_lsh_worker(lsh_root_str: str, cache_max_gb: float) -> None:
    global _WORKER_LSH_ROOT, _WORKER_CACHE, _WORKER_CACHE_BYTES, _WORKER_CACHE_MAX_BYTES
    _WORKER_LSH_ROOT = _p.Path(lsh_root_str)
    _WORKER_CACHE = OrderedDict()
    _WORKER_CACHE_BYTES = 0
    _WORKER_CACHE_MAX_BYTES = int(cache_max_gb * (1024 ** 3))


def _query_bin_worker_stateless(args: Tuple[int, str, List[int], int, int]) -> Dict[str, Any]:
    b, lsh_root_str, qlist, per_bin_k, kc = args
    t0 = time.perf_counter()
    try:
        lsh_path = str(_p.Path(lsh_root_str) / f"{LSH_PREFIX}{pad4(b)}{LSH_SUFFIX}")
        lf = tm.LSHForest()
        t_load_start = time.perf_counter()
        lf.restore(lsh_path)
        if not lf.is_clean():
            lf.index()
        load_s = time.perf_counter() - t_load_start
        qvec = VectorUint([int(x) for x in qlist])
        t_query_start = time.perf_counter()
        try:
            ids = lf.query(qvec, per_bin_k, kc)
            path = "query(k,kc)"
        except TypeError:
            ids = lf.query(qvec, per_bin_k)
            path = "query(k)"
        query_s = time.perf_counter() - t_query_start
        return {
            "bin": b,
            "ids": [int(i) for i in ids],
            "load_s": float(load_s),
            "query_s": float(query_s),
            "total_s": float(time.perf_counter() - t0),
            "path": path,
            "cache_hit": False,
        }
    except Exception as e:
        return {
            "bin": b,
            "ids": [],
            "load_s": float("nan"),
            "query_s": float("nan"),
            "total_s": float(time.perf_counter() - t0),
            "path": f"error:{e}",
            "cache_hit": False,
        }


def _query_bin_worker_cached(args: Tuple[int, List[int], int, int]) -> Dict[str, Any]:
    global _WORKER_LSH_ROOT, _WORKER_CACHE, _WORKER_CACHE_BYTES, _WORKER_CACHE_MAX_BYTES
    b, qlist, per_bin_k, kc = args
    if _WORKER_LSH_ROOT is None:
        raise RuntimeError("Worker not initialized. Use initializer=_init_lsh_worker.")
    t0 = time.perf_counter()
    hit = False
    load_s = 0.0
    lf = _WORKER_CACHE.get(b)
    if lf is None:
        need = _bin_file_size_bytes(b, _WORKER_LSH_ROOT)
        while _WORKER_CACHE and (_WORKER_CACHE_BYTES + need > _WORKER_CACHE_MAX_BYTES):
            old_b, old_lf = _WORKER_CACHE.popitem(last=False)
            _WORKER_CACHE_BYTES -= _bin_file_size_bytes(old_b, _WORKER_LSH_ROOT)
            del old_lf
        lf = tm.LSHForest()
        t_load = time.perf_counter()
        lf.restore(str(lsh_path_for_bin(b, _WORKER_LSH_ROOT)))
        if not lf.is_clean():
            lf.index()
        load_s = time.perf_counter() - t_load
        _WORKER_CACHE[b] = lf
        _WORKER_CACHE.move_to_end(b)
        _WORKER_CACHE_BYTES += need
    else:
        hit = True
        _WORKER_CACHE.move_to_end(b)
    qvec = VectorUint([int(x) for x in qlist])
    t_q = time.perf_counter()
    try:
        ids = lf.query(qvec, per_bin_k, kc)
        path = "query(k,kc)"
    except TypeError:
        ids = lf.query(qvec, per_bin_k)
        path = "query(k)"
    query_s = time.perf_counter() - t_q
    return {
        "bin": b,
        "ids": [int(i) for i in ids],
        "load_s": float(load_s),
        "query_s": float(query_s),
        "total_s": float(time.perf_counter() - t0),
        "path": path,
        "cache_hit": bool(hit),
    }


def lsh_query_bins_parallel(
    ex: cf.ProcessPoolExecutor,
    mode: str,
    qvec: VectorUint,
    bins: Sequence[int],
    lsh_root: _p.Path,
    per_bin_k: int,
    kc: int,
) -> Tuple[List[Tuple[int, int]], Dict[str, Any]]:
    if not bins:
        return [], {"load_s": 0.0, "query_s": 0.0, "total_s": 0.0, "bins_ok": 0, "bins_err": 0, "bins": 0, "cache_hits": 0}
    qlist = [int(x) for x in qvec]
    if mode == "cold":
        tasks = [(b, str(lsh_root), qlist, per_bin_k, kc) for b in bins]
        worker_fn = _query_bin_worker_stateless
    elif mode == "cached":
        tasks = [(b, qlist, per_bin_k, kc) for b in bins]
        worker_fn = _query_bin_worker_cached
    else:
        raise ValueError(f"Unknown mode: {mode}")
    out: List[Tuple[int, int]] = []
    load_sum = 0.0
    query_sum = 0.0
    total_sum = 0.0
    bins_ok = 0
    bins_err = 0
    cache_hits = 0
    for res in ex.map(worker_fn, tasks, chunksize=1):
        b = int(res.get("bin"))
        ids = res.get("ids") or []
        load_s = float(res.get("load_s", float("nan")))
        query_s = float(res.get("query_s", float("nan")))
        total_s = float(res.get("total_s", float("nan")))
        cache_hits += int(bool(res.get("cache_hit", False)))
        if ids:
            bins_ok += 1
            out.extend((b, int(nid)) for nid in ids)
        else:
            bins_err += 1
        if not np.isnan(load_s):
            load_sum += load_s
        if not np.isnan(query_s):
            query_sum += query_s
        if not np.isnan(total_s):
            total_sum += total_s
    return out, {
        "load_s": load_sum,
        "query_s": query_sum,
        "total_s": total_sum,
        "bins_ok": bins_ok,
        "bins_err": bins_err,
        "bins": len(bins),
        "cache_hits": cache_hits,
    }


def _bin_subprocess_script() -> str:
    return r'''
import json, sys, time, pathlib as _p
import tmap as tm
from tmap import VectorUint

LSH_PREFIX = "100M_ZINC_ChEMBL_BIN-"
LSH_SUFFIX = "-512-32.lsh"

def pad4(n):
    return f"{n:04d}"

payload = json.loads(sys.stdin.read())
b = int(payload["bin"])
lsh_root = _p.Path(payload["lsh_root"])
qlist = [int(x) for x in payload["qlist"]]
per_bin_k = int(payload["per_bin_k"])
kc = int(payload["kc"])

lsh_path = str(lsh_root / f"{LSH_PREFIX}{pad4(b)}{LSH_SUFFIX}")
lf = tm.LSHForest()
t0 = time.perf_counter()
t_load = time.perf_counter()
lf.restore(lsh_path)
if not lf.is_clean():
    lf.index()
load_s = time.perf_counter() - t_load
qvec = VectorUint(qlist)
t_q = time.perf_counter()
try:
    ids = lf.query(qvec, per_bin_k, kc)
    path = "query(k,kc)"
except TypeError:
    ids = lf.query(qvec, per_bin_k)
    path = "query(k)"
query_s = time.perf_counter() - t_q
print(json.dumps({
    "bin": b,
    "ids": [int(i) for i in ids],
    "load_s": float(load_s),
    "query_s": float(query_s),
    "total_s": float(time.perf_counter() - t0),
    "path": path,
    "cache_hit": False,
}))
'''


def _query_bin_via_subprocess(
    bin_id: int,
    lsh_root: _p.Path,
    qlist: List[int],
    per_bin_k: int,
    kc: int,
) -> Dict[str, Any]:
    payload = {
        "bin": int(bin_id),
        "lsh_root": str(lsh_root),
        "qlist": [int(x) for x in qlist],
        "per_bin_k": int(per_bin_k),
        "kc": int(kc),
    }
    proc = subprocess.run(
        [sys.executable, "-c", _bin_subprocess_script()],
        input=json.dumps(payload),
        text=True,
        capture_output=True,
    )
    stdout_lines = [x for x in proc.stdout.splitlines() if x.strip()]
    if proc.returncode != 0:
        raise RuntimeError(
            f"worker produced no JSON for bin={bin_id}; rc={proc.returncode}; "
            f"stderr={proc.stderr.strip()}; stdout_tail={stdout_lines[-3:]}"
        )
    for line in reversed(stdout_lines):
        s = line.strip()
        if s.startswith("{") and s.endswith("}"):
            return json.loads(s)
    raise RuntimeError(
        f"worker produced no JSON for bin={bin_id}; rc={proc.returncode}; "
        f"stderr={proc.stderr.strip()}; stdout_tail={stdout_lines[-3:]}"
    )


def lsh_query_bins_serial(
    mode: str,
    qvec: VectorUint,
    bins: Sequence[int],
    lsh_root: _p.Path,
    per_bin_k: int,
    kc: int,
    cache_max_gb_per_worker: float = 20.0,
    debug_bins: bool = False,
    bin_retry_max: int = 3,
) -> Tuple[List[Tuple[int, int]], Dict[str, Any]]:
    if not bins:
        return [], {
            "load_s": 0.0, "query_s": 0.0, "total_s": 0.0,
            "bins_ok": 0, "bins_err": 0, "bins": 0, "cache_hits": 0,
            "bin_failures": []
        }

    qlist = [int(x) for x in qvec]
    out: List[Tuple[int, int]] = []
    load_sum = 0.0
    query_sum = 0.0
    total_sum = 0.0
    bins_ok = 0
    bins_err = 0
    cache_hits = 0
    bin_failures: List[int] = []

    if mode == "cached":
        _init_lsh_worker(str(lsh_root), float(cache_max_gb_per_worker))

    for b in bins:
        if debug_bins:
            print(f"[BIN] start bin={b} path={lsh_path_for_bin(b, lsh_root)}")

        last_err = None
        res = None
        for attempt in range(1, max(1, int(bin_retry_max)) + 1):
            try:
                if mode == "cold":
                    res = _query_bin_via_subprocess(b, lsh_root, qlist, per_bin_k, kc)
                else:
                    res = _query_bin_worker_cached((b, qlist, per_bin_k, kc))
                break
            except Exception as e:
                last_err = e
                print(f"[WARN] bin={b} attempt={attempt}/{bin_retry_max} failed: {e}")
                time.sleep(0.2)

        if res is None:
            bins_err += 1
            bin_failures.append(int(b))
            print(f"[ERROR] bin={b} skipped after {bin_retry_max} attempts: {last_err}")
            continue

        ids = res.get("ids") or []
        load_s = float(res.get("load_s", float("nan")))
        query_s = float(res.get("query_s", float("nan")))
        total_s = float(res.get("total_s", float("nan")))
        cache_hits += int(bool(res.get("cache_hit", False)))

        if ids:
            bins_ok += 1
            out.extend((b, int(nid)) for nid in ids)
        else:
            bins_err += 1
            bin_failures.append(int(b))
            if debug_bins:
                print(f"[BIN-WARN] bin={b} returned no ids path={res.get('path', '')}")

        if not np.isnan(load_s):
            load_sum += load_s
        if not np.isnan(query_s):
            query_sum += query_s
        if not np.isnan(total_s):
            total_sum += total_s

        if debug_bins:
            print(f"[BIN] done bin={b} n_ids={len(ids)}")

    return out, {
        "load_s": load_sum,
        "query_s": query_sum,
        "total_s": total_sum,
        "bins_ok": bins_ok,
        "bins_err": bins_err,
        "bins": len(bins),
        "cache_hits": cache_hits,
        "bin_failures": bin_failures,
    }


def canonicalize(smi: str) -> Optional[str]:
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception:
        return None


def _unwrap_csv_scalar(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        if not value:
            return ""
        return str(value[0]).strip()
    text = str(value).strip()
    if text.startswith("[") and text.endswith("]"):
        try:
            parsed = ast.literal_eval(text)
            if isinstance(parsed, list):
                if not parsed:
                    return ""
                return str(parsed[0]).strip()
        except Exception:
            pass
    return text


def truth_file_for(drug_name: str, sum_bin: int, vs_dir: _p.Path) -> Optional[_p.Path]:
    drug_norm = norm_drug_name(drug_name)
    candidates = [
        vs_dir / f"{drug_norm}-{sum_bin}_MHFP-512.csv",
        vs_dir / f"{drug_norm}_BIN-{sum_bin:04d}_MHFP-512.csv",
    ]
    for path in candidates:
        if path.exists():
            return path

    patterns = [
        f"{drug_norm}_BIN-{sum_bin:04d}_*_MHFP-512.csv",
        f"{drug_norm}*_BIN-{sum_bin:04d}_*_MHFP-512.csv",
        f"{drug_norm}*-{sum_bin}_MHFP-512.csv",
    ]
    matches: list[_p.Path] = []
    for pattern in patterns:
        matches.extend(sorted(vs_dir.glob(pattern)))
    return matches[0] if matches else None


def load_truth_smiles_topk(
    drug_name: str,
    sum_bin: int,
    vs_dir: _p.Path,
    k: int = 1000,
    drop_self: bool = False,
) -> List[str]:
    p = truth_file_for(drug_name, sum_bin, vs_dir)
    if p is None:
        print(f"[WARN] truth file not found for {drug_name} bin={sum_bin}")
        return []
    rows: List[Tuple[float, str]] = []
    with p.open() as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            try:
                if row.get("Distance") not in (None, ""):
                    score = -float(row["Distance"])
                else:
                    score = float(row.get("similarity", "-inf"))
            except ValueError:
                continue
            smi = _unwrap_csv_scalar(row.get("SMILES") or row.get("smiles"))
            if smi:
                rows.append((score, smi))
    rows.sort(key=lambda x: x[0], reverse=True)
    if drop_self and rows and rows[0][0] >= 1.0:
        rows = rows[1:]
    return [s for _, s in rows[:k]]


def _canon_unique(smiles: Sequence[str], limit: Optional[int] = None) -> List[str]:
    out: List[str] = []
    seen: set[str] = set()
    for s in smiles:
        c = canonicalize(s)
        if not c or c in seen:
            continue
        seen.add(c)
        out.append(c)
        if limit is not None and len(out) >= limit:
            break
    return out


def recall_at_k_strict(approx_smiles: Sequence[str], true_smiles: Sequence[str], k: int) -> Tuple[float, int]:
    ct_k = set(_canon_unique(true_smiles, limit=k))
    eff_k = min(k, len(ct_k))
    if eff_k == 0:
        return 0.0, 0
    ca = set(_canon_unique(approx_smiles))
    return len(ca & ct_k) / float(eff_k), eff_k


def rerank_hybrid_take_heap_timed(
    query_smiles: str,
    pairs: List[Tuple[int, int]],
    idx_dir: _p.Path,
    topk: int,
    per_bin_topk: int,
    chunk_size: int = 20000,
) -> Tuple[List[Tuple[int, int, str, float]], Dict[str, float]]:
    qsig = np.asarray(list(_calc_mhfp(query_smiles)), dtype=np.uint32)
    by_bin: Dict[int, List[int]] = {}
    for b, i in pairs:
        by_bin.setdefault(b, []).append(i)

    global_heap: List[Tuple[float, int, int, str]] = []
    load_fp_s = 0.0
    score_heap_s = 0.0

    for b, ids in by_bin.items():
        seen = set()
        ids = [i for i in ids if (i not in seen) and (not seen.add(i))]
        if not ids:
            continue

        p = idx_dir / f"index_smiles_BIN-{b:04d}.parquet"
        if not p.exists():
            raise FileNotFoundError(f"Missing IDX parquet for bin={b}: {p}")

        t_load = time.perf_counter()
        tbl = pq.read_table(p, columns=["smiles", "fp"], memory_map=True)
        n_rows = tbl.num_rows
        load_fp_s += time.perf_counter() - t_load
        if n_rows == 0:
            continue

        ids = [i for i in ids if 0 <= i < n_rows]
        if not ids:
            continue

        if chunk_size <= 0:
            chunk_size = len(ids)

        bin_heap: List[Tuple[float, int, int, str]] = []

        for start in range(0, len(ids), chunk_size):
            chunk_ids = ids[start:start + chunk_size]

            t_load = time.perf_counter()
            taken = pac.take(tbl, pa.array(chunk_ids, type=pa.int32()))
            smiles = taken["smiles"].to_pylist()
            fp_list_arr = taken["fp"].combine_chunks()
            values = fp_list_arr.values
            offsets = fp_list_arr.offsets
            vals = values.to_numpy(zero_copy_only=False)
            offs = offsets.to_numpy()
            lens = offs[1:] - offs[:-1]
            if not (lens == 512).all():
                raise RuntimeError(f"Unexpected fp length(s) in bin {b}: {np.unique(lens)}")
            F = vals.reshape(len(chunk_ids), 512)
            load_fp_s += time.perf_counter() - t_load

            t_score = time.perf_counter()
            sims = (F == qsig).mean(axis=1)
            for (idv, smi, sim) in zip(chunk_ids, smiles, sims):
                item = (float(sim), b, int(idv), smi)
                if len(bin_heap) < per_bin_topk:
                    heapq.heappush(bin_heap, item)
                else:
                    if item[0] > bin_heap[0][0]:
                        heapq.heapreplace(bin_heap, item)
            score_heap_s += time.perf_counter() - t_score

        for item in bin_heap:
            if len(global_heap) < topk:
                heapq.heappush(global_heap, item)
            else:
                if item[0] > global_heap[0][0]:
                    heapq.heapreplace(global_heap, item)

    global_heap.sort(reverse=True)
    out = [(b, idv, smi, sim) for (sim, b, idv, smi) in global_heap]
    stats = {
        "load_fp_s": load_fp_s,
        "score_heap_s": score_heap_s,
        "total_s": load_fp_s + score_heap_s,
    }
    return out, stats


def load_drugs145(csv_path: _p.Path, mqn_bin: _p.Path) -> List[Tuple[str, str, int]]:
    out: List[Tuple[str, str, int]] = []
    with csv_path.open() as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            smi = (row.get("SMILES") or "").strip()
            drug = (row.get("DRUG_NAME") or "").strip()
            if not smi or not drug:
                continue
            sum_bin_raw = (row.get("SUM_BIN") or "").strip()
            if sum_bin_raw:
                try:
                    sbin = int(float(sum_bin_raw))
                except Exception:
                    sbin = None
            else:
                sbin = mqn_sum_stdio(smi, mqn_bin)
            if sbin is None:
                print(f"[WARN] MQN sum failed; skipping drug={drug}")
                continue
            out.append((drug, smi, sbin))
    return out


def run_benchmark(
    paths: Paths,
    cfg: RunCfg,
    mode: str,
    only_drug: str = "",
    cache_max_gb_per_worker: float = 0.0,
    require_manifest: bool = False,
    resume: bool = False,
    debug_bins: bool = False,
) -> None:
    global _IDX_CACHE_MAX_BINS
    _IDX_CACHE_MAX_BINS = int(cfg.idx_cache_bins)
    ensure_dirs(paths)

    print(f"[CWD] {os.getcwd()}")
    print(f"[OUT] out_dir={paths.out_dir.resolve()}")
    print(f"[OUT] per_drug_dir={paths.per_drug_dir.resolve()}")

    if not _env_ok():
        print("The test of environment did not match the reference mhfp. Refusing to continue.")
        sys.exit(2)
    print("The Python environment is consistent. Will continue...")

    _write_headers_if_missing(paths.approx_ids_csv, ["drug_name", "rank", "bin", "nn_id", "bins_used"])
    _write_headers_if_missing(paths.timings_csv, ["drug_name", "stage", "seconds"])
    _write_headers_if_missing(paths.metrics_csv, ["drug_name", "recall_at_k", "k", "per_bin_candidates_k", "kc", "bins_used"])
    _write_headers_if_missing(
        paths.queries_csv,
        ["drug_name", "smiles", "mqn_sum_bin", "bins_used", "bins_src", "k", "per_bin_candidates_k", "kc", "span", "encoder", "mode", "rerank_per_bin_topk", "bin_retry_max"],
    )
    _write_headers_if_missing(
        paths.status_csv,
        ["drug_name", "sum_bin", "span", "per_bin", "mode", "unix_time", "status", "message"],
    )

    t0 = time.perf_counter()
    queries = load_drugs145(paths.drugs145_csv, paths.mqn_bin)
    with paths.timings_csv.open("a", newline="") as f_tm:
        csv.writer(f_tm).writerow(["__all__", "load_drugs145", f"{time.perf_counter() - t0:.6f}"])

    if only_drug:
        only_norm = norm_drug_name(only_drug)
        queries = [(d, s, b) for (d, s, b) in queries if norm_drug_name(d) == only_norm]
        if not queries:
            print(f"[WARN] No drugs matched only_drug={only_drug}")
            return

    if resume:
        done = _load_completed_drugs(paths.status_csv, cfg.span, cfg.per_bin_candidates_k, mode)
        if done:
            before = len(queries)
            queries = [(d, s, b) for (d, s, b) in queries if d not in done]
            print(f"[RESUME] skipping {before - len(queries)} completed drugs for span={cfg.span} per_bin={cfg.per_bin_candidates_k} mode={mode}")
        if not queries:
            print("[RESUME] nothing left to do")
            return

    ex = None
    if cfg.lsh_workers > 1:
        if mode == "cached":
            if cache_max_gb_per_worker <= 0.0:
                cache_max_gb_per_worker = 20.0
            ex = cf.ProcessPoolExecutor(
                max_workers=cfg.lsh_workers,
                initializer=_init_lsh_worker,
                initargs=(str(paths.lsh_root), float(cache_max_gb_per_worker)),
            )
        elif mode == "cold":
            ex = cf.ProcessPoolExecutor(max_workers=cfg.lsh_workers)
        else:
            raise ValueError(f"Unknown mode: {mode}")
    else:
        if mode not in {"cold", "cached"}:
            raise ValueError(f"Unknown mode: {mode}")

    try:
        for idx, (drug, smi, sbin) in enumerate(queries, 1):
            drug_t0 = time.perf_counter()
            bins, bins_src = bins_for_drug(paths, drug, sbin, cfg.span, require_manifest=require_manifest)
            bins_str = ";".join(str(b) for b in bins)

            print(
                f"[{idx}/{len(queries)}] "
                f"Drug={drug} | MQN_bin={sbin} | bins={bins_str or '<EMPTY>'} | bins_src={bins_src} | "
                f"mode={mode} | k={cfg.k} per_bin={cfg.per_bin_candidates_k} kc={cfg.kc}"
            )

            with paths.queries_csv.open("a", newline="") as f_q:
                csv.writer(f_q).writerow([
                    drug, smi, sbin, bins_str, bins_src, cfg.k, cfg.per_bin_candidates_k,
                    cfg.kc, cfg.span, "sacha_mhfp512", mode, cfg.rerank_per_bin_topk, cfg.bin_retry_max,
                ])

            reranked_path = paths.per_drug_dir / f"{norm_drug_name(drug)}-{sbin}_approx_top{cfg.k}_reranked.csv"
            err_path = reranked_path.with_suffix(".ERROR.txt")
            _append_status(paths.status_csv, drug, sbin, cfg.span, cfg.per_bin_candidates_k, mode, "STARTED", bins_str)

            t = time.perf_counter()
            try:
                qvec = make_query_vector_from_smiles(smi)
            except Exception as e:
                msg = f"MHFP encoding failed for drug={drug} bin={sbin}: {e}"
                print(f"[ERROR] {msg}")
                _write_error_file(err_path, msg)
                _append_status(paths.status_csv, drug, sbin, cfg.span, cfg.per_bin_candidates_k, mode, "FAILED", msg)
                with paths.timings_csv.open("a", newline="") as f_tm:
                    csv.writer(f_tm).writerow([drug, "mhfp_sacha_encoding_error", "nan"])
                continue
            with paths.timings_csv.open("a", newline="") as f_tm:
                csv.writer(f_tm).writerow([drug, "mhfp_sacha_encoding", f"{time.perf_counter() - t:.6f}"])

            t = time.perf_counter()
            try:
                if cfg.lsh_workers == 1:
                    pairs, lsh_stats = lsh_query_bins_serial(
                        mode=mode,
                        qvec=qvec,
                        bins=bins,
                        lsh_root=paths.lsh_root,
                        per_bin_k=cfg.per_bin_candidates_k,
                        kc=cfg.kc,
                        cache_max_gb_per_worker=cache_max_gb_per_worker,
                        debug_bins=debug_bins,
                        bin_retry_max=cfg.bin_retry_max,
                    )
                else:
                    pairs, lsh_stats = lsh_query_bins_parallel(
                        ex=ex,
                        mode=mode,
                        qvec=qvec,
                        bins=bins,
                        lsh_root=paths.lsh_root,
                        per_bin_k=cfg.per_bin_candidates_k,
                        kc=cfg.kc,
                    )
            except cf.process.BrokenProcessPool as e:
                msg = (
                    f"LSH worker pool failed for drug={drug} bin={sbin} bins={bins_str} "
                    f"mode={mode} workers={cfg.lsh_workers}: {e}"
                )
                print(f"[ERROR] {msg}")
                _write_error_file(err_path, msg)
                _append_status(paths.status_csv, drug, sbin, cfg.span, cfg.per_bin_candidates_k, mode, "FAILED", msg)
                with paths.timings_csv.open("a", newline="") as f_tm:
                    csv.writer(f_tm).writerow([drug, "lsh_broken_process_pool_error", "nan"])
                continue
            wall = time.perf_counter() - t

            with paths.timings_csv.open("a", newline="") as f_tm:
                w = csv.writer(f_tm)
                w.writerow([drug, "lsh_wall_s", f"{wall:.6f}"])
                w.writerow([drug, "lsh_load_sum_s", f"{float(lsh_stats['load_s']):.6f}"])
                w.writerow([drug, "lsh_query_sum_s", f"{float(lsh_stats['query_s']):.6f}"])
                w.writerow([drug, "lsh_total_sum_s", f"{float(lsh_stats['total_s']):.6f}"])
                w.writerow([drug, "lsh_bins_ok", str(int(lsh_stats['bins_ok']))])
                w.writerow([drug, "lsh_bins_err", str(int(lsh_stats['bins_err']))])
                w.writerow([drug, "lsh_cache_hits", str(int(lsh_stats.get('cache_hits', 0)))])
                w.writerow([drug, "lsh_bin_failures_count", str(len(lsh_stats.get('bin_failures', [])))])

            print(
                f"  -> got {len(pairs)} candidates from {len(bins)} bins "
                f"(~{cfg.per_bin_candidates_k} each), bins_ok={lsh_stats['bins_ok']} bins_err={lsh_stats['bins_err']} "
                f"cache_hits={lsh_stats.get('cache_hits', 0)} rerank_per_bin_topk={cfg.rerank_per_bin_topk}"
            )

            with paths.approx_ids_csv.open("a", newline="") as f_nn:
                w = csv.writer(f_nn)
                for r, (b, nid) in enumerate(pairs[: min(cfg.k, 200)], 1):
                    w.writerow([drug, r, b, str(nid), bins_str])

            try:
                final_top, rr_stats = rerank_hybrid_take_heap_timed(
                    smi,
                    pairs,
                    paths.idx_dir,
                    topk=cfg.k,
                    per_bin_topk=cfg.rerank_per_bin_topk,
                    chunk_size=20000,
                )
            except Exception as e:
                msg = f"rerank failed for drug={drug} bin={sbin} pairs={len(pairs)} idx_dir={paths.idx_dir}: {e}"
                print(f"[ERROR] {msg}")
                _write_error_file(err_path, msg)
                _append_status(paths.status_csv, drug, sbin, cfg.span, cfg.per_bin_candidates_k, mode, "FAILED", msg)
                with paths.timings_csv.open("a", newline="") as f_tm:
                    csv.writer(f_tm).writerow([drug, "rerank_error", "nan"])
                continue

            with paths.timings_csv.open("a", newline="") as f_tm:
                w = csv.writer(f_tm)
                w.writerow([drug, "rerank_load_fp_s", f"{rr_stats['load_fp_s']:.6f}"])
                w.writerow([drug, "rerank_score_heap_s", f"{rr_stats['score_heap_s']:.6f}"])
                w.writerow([drug, "rerank_total_s", f"{rr_stats['total_s']:.6f}"])

            with reranked_path.open("w", newline="") as f_rr:
                w = csv.writer(f_rr)
                w.writerow(["drug_name", "rank", "bin", "nn_id", "smiles", "mhfp_minhash_jaccard"])
                for r, (b, nid, smi_nn, sim) in enumerate(final_top, 1):
                    w.writerow([drug, r, b, nid, smi_nn, f"{sim:.6f}"])

            if err_path.exists():
                try:
                    err_path.unlink()
                except Exception:
                    pass

            print(f"[WRITE] {reranked_path.resolve()} n_final={len(final_top)}")

            truth_smiles = load_truth_smiles_topk(drug, sbin, paths.vs_drugs_dir, k=cfg.k, drop_self=False)
            approx_smiles = [smi_nn for (_, _, smi_nn, _) in final_top]
            r_at_k, eff_k = recall_at_k_strict(approx_smiles, truth_smiles, cfg.k)

            with paths.metrics_csv.open("a", newline="") as f_mt:
                csv.writer(f_mt).writerow([
                    drug,
                    f"{r_at_k:.6f}",
                    cfg.k,
                    f"{cfg.per_bin_candidates_k}x{len(bins)}",
                    cfg.kc,
                    bins_str + f" | eff_k={eff_k} | mode={mode} | bin_failures={';'.join(map(str, lsh_stats.get('bin_failures', [])))}",
                ])

            _append_status(
                paths.status_csv,
                drug,
                sbin,
                cfg.span,
                cfg.per_bin_candidates_k,
                mode,
                "DONE",
                f"wall_s={time.perf_counter() - drug_t0:.6f}",
            )

    finally:
        if ex is not None:
            ex.shutdown(wait=True)

    print(
        "done. wrote:\n"
        f"- {paths.per_drug_dir}/*_reranked.csv\n"
        f"- {paths.approx_ids_csv}\n"
        f"- {paths.timings_csv}\n"
        f"- {paths.metrics_csv}\n"
        f"- {paths.queries_csv}\n"
        f"- {paths.status_csv}"
    )


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="LSHForest benchmark with robust per-bin retry and split rerank timings.")
    ap.add_argument("--k", type=int, default=50)
    ap.add_argument("--per-bin", type=int, default=5000, dest="per_bin_candidates_k")
    ap.add_argument("--span", type=int, default=0)
    ap.add_argument("--kc", type=int, default=-1)
    ap.add_argument("--lsh-workers", type=int, default=max(1, (os.cpu_count() or 4) // 2))
    ap.add_argument("--mode", choices=["cold", "cached"], default="cold")
    ap.add_argument("--cache-max-gb", type=float, default=20.0)
    ap.add_argument("--idx-cache-bins", type=int, default=0)
    ap.add_argument("--rerank-per-bin-topk", type=int, default=250)
    ap.add_argument("--bin-retry-max", type=int, default=3)
    ap.add_argument("--out-tag", type=str, default="")
    ap.add_argument("--only-drug", type=str, default="")
    ap.add_argument("--require-manifest", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--debug-bins", action="store_true")
    return ap.parse_args()


def build_paths_with_tag(base: Paths, out_tag: str) -> Paths:
    if not out_tag:
        return base
    root = REPO_ROOT / "benchmark_out" / out_tag
    return dataclasses.replace(
        base,
        out_dir=root,
        per_drug_dir=root / "per_drug",
        approx_ids_csv=root / "approx_ids.csv",
        timings_csv=root / "timings.csv",
        metrics_csv=root / "metrics.csv",
        queries_csv=root / "queries.csv",
        status_csv=root / "status.csv",
    )


def main() -> None:
    args = parse_args()
    P = build_paths_with_tag(Paths(), args.out_tag)
    C = RunCfg(
        k=int(args.k),
        per_bin_candidates_k=int(args.per_bin_candidates_k),
        kc=int(args.kc),
        span=int(args.span),
        lsh_workers=int(args.lsh_workers),
        idx_cache_bins=int(args.idx_cache_bins),
        rerank_per_bin_topk=int(args.rerank_per_bin_topk),
        bin_retry_max=int(args.bin_retry_max),
    )
    print(
        f"[CFG] k={C.k} per_bin={C.per_bin_candidates_k} kc={C.kc} span={C.span} "
        f"lsh_workers={C.lsh_workers} mode={args.mode} cache_max_gb={args.cache_max_gb} "
        f"idx_cache_bins={C.idx_cache_bins} require_manifest={bool(args.require_manifest)} "
        f"resume={bool(args.resume)} debug_bins={bool(args.debug_bins)} "
        f"rerank_per_bin_topk={C.rerank_per_bin_topk} bin_retry_max={C.bin_retry_max}"
    )
    run_benchmark(
        paths=P,
        cfg=C,
        mode=str(args.mode),
        only_drug=str(args.only_drug),
        cache_max_gb_per_worker=float(args.cache_max_gb),
        require_manifest=bool(args.require_manifest),
        resume=bool(args.resume),
        debug_bins=bool(args.debug_bins),
    )


if __name__ == "__main__":
    main()
