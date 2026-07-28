"""Genome-scale profile of the calibration solver.

Runs `calibrate` on a cached REAL cfRNA payload (human index v7, ~10^6 nodes) and on the 10 Mb synthetic,
so each hotspot cluster's SCALING can be read off directly rather than extrapolated.

    OMP_NUM_THREADS=1 python scratchpad/perf_scale.py --input cfrna:LBX0190 --refit 1 --profile
    OMP_NUM_THREADS=1 python scratchpad/perf_scale.py --input toy --refit 1 --profile
"""
from __future__ import annotations

import argparse
import cProfile
import dataclasses
import os
import pickle
import pstats
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")

import importlib  # noqa: E402

from rigel.config import PipelineConfig  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

CF = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna/_calib_cache")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
TOY_COND = "gdna_gdna300_ss_0.50_nrna_present_capture_off"


def load(spec: str):
    if spec.startswith("cfrna:"):
        d = pickle.load(open(CF / f"{spec.split(':', 1)[1]}.pkl", "rb"))
        return (
            d["payload"],
            d["region_arrays"],
            d["strand_model"],
            np.asarray(d["gdna_fl_pmf"]),
            np.asarray(d["rna_fl_pmf"]),
        )
    from selfsolve_diag import _scan_and_truth

    from rigel.index import TranscriptIndex
    from rigel.calibration.region_arrays import RegionArrays

    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    cond = spec.split(":", 1)[1] if ":" in spec else TOY_COND
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    return (
        inp["payload"],
        ra,
        inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
    )


ap = argparse.ArgumentParser()
ap.add_argument("--input", default="toy")
ap.add_argument("--refit", type=int, default=1)
ap.add_argument("--profile", action="store_true")
ap.add_argument("--repeat", type=int, default=1)
ap.add_argument("--top", type=int, default=25)
a = ap.parse_args()

payload, ra, sm, gpmf, rpmf = load(a.input)
cc = dataclasses.replace(PipelineConfig().calibration, calib_refit_iters=a.refit)
n_reg = len(ra.region_size_bp)
print(f"input={a.input} refit={a.refit}  regions={n_reg}  (nodes ≈ 2·regions+1)")


def run():
    dbg: dict = {}
    calmod.calibrate(payload, ra, sm, gpmf, rpmf, cc, _debug=dbg)
    return dbg


t0 = time.perf_counter()
run()
print(f"warm/first call: {time.perf_counter() - t0:.3f} s")
ts = []
for _ in range(a.repeat):
    t0 = time.perf_counter()
    run()
    ts.append(time.perf_counter() - t0)
if ts:
    print(f"wall: min {min(ts):.3f} s  median {sorted(ts)[len(ts) // 2]:.3f} s  over {len(ts)}")

if a.profile:
    pr = cProfile.Profile()
    pr.enable()
    run()
    pr.disable()
    st = pstats.Stats(pr)
    st.sort_stats("tottime").print_stats(a.top)
