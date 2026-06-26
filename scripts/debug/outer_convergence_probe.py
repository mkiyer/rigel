"""Measure the NESTED-loop convergence: how many OUTER (var~mean fixed-point) iterations are needed.

Each outer iteration: inner-converge the belief at fixed var~mean, then refit BOTH var~mean curves on the
converged belief. node_sweep returns the per-outer-iteration max-Δf_g of the converged belief (`outer_deltas`).
This runs calibrate() with the outer cap cranked up + a tiny outer tol and prints that trajectory per condition.

    OMP_NUM_THREADS=1 python scripts/debug/outer_convergence_probe.py [cond ...]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import dataclasses
import importlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

DEFAULT = [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_none_ss_0.50_nrna_none_capture_on",
]


def run(cond, max_outer=8):
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    cfg = dataclasses.replace(
        CalibrationConfig(), sweep_max_outer=max_outer, sweep_outer_convergence_delta=1e-9
    )
    calmod = importlib.import_module("rigel.calibration.calibrate")
    orig = calmod.node_sweep
    cap = {}

    def wrapped(*a, **k):
        belief, deltas = orig(*a, **k)
        cap["deltas"] = deltas
        return belief, deltas

    calmod.node_sweep = wrapped
    try:
        cal = calibrate(
            payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
            gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=cfg,
        )
    finally:
        calmod.node_sweep = orig
    return cap.get("deltas", []), cal


def main():
    print("NESTED-loop OUTER convergence (per-outer max-Δf_g of the inner-converged belief; cap=8, tol=1e-9)")
    for c in (sys.argv[1:] or DEFAULT):
        try:
            d, cal = run(c)
            traj = "  ".join(f"{x:.3f}" for x in d)
            print(f"\n=== {c} ===")
            print(f"  outer iterations run = {len(d) + 1}   gdna_density_global={cal.gdna_density_global:.4g}")
            print(f"  per-outer Δ (iter k vs k-1): {traj}")
        except Exception as e:  # noqa: BLE001
            import traceback
            print(f"\n=== {c} === FAILED: {e}")
            traceback.print_exc()


if __name__ == "__main__":
    main()
