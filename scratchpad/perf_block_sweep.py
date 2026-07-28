"""Sweep `_SOLVE_BLOCK_BYTES` end-to-end at genome scale, and prove each setting bit-identical.

The block size cannot reach the arithmetic (every reduction is within a row), so this is purely a cache
measurement — but it is one that must be taken at the real shapes: at 3.4 k nodes the 1-D path is a single
2 MB array and every setting looks the same.
"""
from __future__ import annotations

import dataclasses
import importlib
import os
import pickle
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

import rigel.calibration.simplex_logodds as slo  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
CF = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna/_calib_cache")

d = pickle.load(open(CF / "LBX0190.pkl", "rb"))
payload, ra, sm = d["payload"], d["region_arrays"], d["strand_model"]
gpmf, rpmf = np.asarray(d["gdna_fl_pmf"]), np.asarray(d["rna_fl_pmf"])
cc = dataclasses.replace(PipelineConfig().calibration, calib_refit_iters=1)

SIZES = [int(x) for x in (sys.argv[1:] or ["1073741824", "16777216", "4194304", "1048576", "262144", "65536"])]
ref = None
for nbytes in SIZES:
    slo._SOLVE_BLOCK_BYTES = nbytes
    t0 = time.perf_counter()
    res = calmod.calibrate(payload, ra, sm, gpmf, rpmf, cc)
    dt = time.perf_counter() - t0
    cur = np.concatenate([res.mass_gdna_contained, res.mass_rna_contained,
                          res.mass_gdna_left, res.mass_rna_left,
                          res.mass_gdna_right, res.mass_rna_right]).astype(np.float64)
    if ref is None:
        ref, tag = cur, "reference"
    else:
        tag = "BIT-IDENTICAL" if np.array_equal(cur.view(np.int64), ref.view(np.int64)) else "⛔ DIFFERS"
    print(f"_SOLVE_BLOCK_BYTES={nbytes:>11} ({nbytes / 2**20:8.2f} MB)   wall={dt:7.2f} s   {tag}", flush=True)
