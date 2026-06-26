"""Outer-loop BULK convergence: mean-Δf_g + #moving nodes between outer iterations (not just max).

The outer max-Δf_g stays high for gDNA-present conditions, but that may be a few knife-edge flip-floppers.
Run calibrate with sweep_max_outer = k for k=1..N (deterministic ⇒ run-k == state after outer iter k of one
run), capture the full converged belief, and report the bulk distribution: mean-Δ, #nodes>.05, #nodes>.20.

    OMP_NUM_THREADS=1 python scripts/debug/outer_meandelta_probe.py [cond] [--n N]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

DEFAULT = "gdna_gdna300_ss_0.50_nrna_none_capture_on"


def belief_at(ra, blob, k):
    cfg = dataclasses.replace(CalibrationConfig(), sweep_max_outer=k, sweep_outer_convergence_delta=1e-12)
    calmod = importlib.import_module("rigel.calibration.calibrate")
    orig = calmod.node_sweep
    cap = {}

    def wrapped(*a, **kw):
        b, d = orig(*a, **kw)
        cap["b"] = b
        return b, d

    calmod.node_sweep = wrapped
    try:
        calibrate(
            payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
            gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=cfg,
        )
    finally:
        calmod.node_sweep = orig
    b = cap["b"]
    return np.asarray(b.f_g), np.asarray(b.f_pos), np.asarray(b.f_neg)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cond = args[0] if args else DEFAULT
    n = int(sys.argv[sys.argv.index("--n") + 1]) if "--n" in sys.argv else 6
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    print(f"\n=== {cond} : outer-loop BULK convergence ===")
    print(f"  {'outer':>5} {'max|Δfg|':>9} {'mean|Δfg|':>9} {'#fg>.05':>8} {'#fg>.20':>8} {'Σf_g':>10}")
    prev = None
    for k in range(1, n + 1):
        fg, fp, fn = belief_at(ra, blob, k)
        if prev is not None:
            d = np.abs(fg - prev)
            print(f"  {k:>5} {d.max():>9.3f} {d.mean():>9.4f} {int((d>0.05).sum()):>8} "
                  f"{int((d>0.20).sum()):>8} {fg.sum():>10,.1f}")
        prev = fg


if __name__ == "__main__":
    main()
