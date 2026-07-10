"""Characterize the Gauss-Seidel sweep convergence (or lack of it).

node_sweep's max-Δf per pass does NOT shrink (oscillates ~0.3-0.6 out to 25 passes). This probe digs in:
runs calibrate() with max_passes = 1..N (deterministic init ⇒ run-k == state after pass k of one long run),
captures the full belief (f_g, f_pos, f_neg), and reports the convergence DISTRIBUTION per step:
  max-Δf_g, mean-Δf_g, # nodes with |Δf_g|>0.05, and whether the SAME nodes keep moving (flip-floppers).

`--freeze-vm`: freeze the RNA var~mean after its first fit (return the cached curve) to test whether the
per-pass var~mean REFIT (a moving precision target) drives the oscillation, vs the GS update structure itself.

    OMP_NUM_THREADS=1 python scripts/debug/gs_convergence_probe.py [cond] [--freeze-vm] [--n N]
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

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

DEFAULT = "gdna_gdna300_ss_0.50_nrna_none_capture_on"


def _belief_at(cond, k, ra, blob, freeze_vm):
    """Run calibrate with max_passes=k; return (f_g, f_pos, f_neg) of the swept belief."""
    cfg = dataclasses.replace(CalibrationConfig(), sweep_max_passes=k, sweep_convergence_delta=1e-12)
    cap = {}
    calmod = importlib.import_module("rigel.calibration.calibrate")  # the MODULE (pkg shadows the name)
    orig_sweep = calmod.node_sweep

    def wrapped(*a, **kw):
        belief, deltas = orig_sweep(*a, **kw)
        cap["b"] = belief
        return belief, deltas

    calmod.node_sweep = wrapped
    orig_fit = bp.fit_rna_varmean
    if freeze_vm:
        store = {}

        def frozen(*a, **kw):
            if "vm" not in store:
                store["vm"] = orig_fit(*a, **kw)
            return store["vm"]

        bp.fit_rna_varmean = frozen
    try:
        calibrate(
            payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
            gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=cfg,
        )
    finally:
        calmod.node_sweep = orig_sweep
        bp.fit_rna_varmean = orig_fit
    b = cap["b"]
    return np.asarray(b.f_g), np.asarray(b.f_pos), np.asarray(b.f_neg)


def run(cond, n, freeze_vm):
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    print(f"\n=== {cond} {'[FROZEN var~mean]' if freeze_vm else '[per-pass refit]'} ===")
    print(f"  {'pass':>4} {'max|Δfg|':>9} {'mean|Δfg|':>9} {'#fg>.05':>8} {'#fg>.20':>8} "
          f"{'max|Δfp|':>9} {'max|Δfn|':>9} {'persist':>8}")
    prev = None
    prev_movers = None
    for k in range(1, n + 1):
        fg, fp, fn = _belief_at(cond, k, ra, blob, freeze_vm)
        if prev is not None:
            dfg = np.abs(fg - prev[0])
            dfp = np.abs(fp - prev[1])
            dfn = np.abs(fn - prev[2])
            movers = set(np.where(dfg > 0.05)[0].tolist())
            persist = len(movers & prev_movers) if prev_movers is not None else 0
            print(f"  {k:>4} {dfg.max():>9.3f} {dfg.mean():>9.4f} {int((dfg>0.05).sum()):>8} "
                  f"{int((dfg>0.20).sum()):>8} {dfp.max():>9.3f} {dfn.max():>9.3f} {persist:>8}")
            prev_movers = movers
        prev = (fg, fp, fn)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    freeze = "--freeze-vm" in sys.argv
    n = 15
    if "--n" in sys.argv:
        n = int(sys.argv[sys.argv.index("--n") + 1])
    cond = args[0] if args else DEFAULT
    run(cond, n, freeze)


if __name__ == "__main__":
    main()
