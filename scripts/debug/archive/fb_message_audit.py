"""Audit the FB message system for OUTRIGHT DEFECTS (before any new weights): out-of-range modes, precision
pathologies, NaN/inf, and the gDNA-emission gate. A message is a Gaussian (mode, prec) on the dst's fraction —
``mode`` should be a valid fraction in [0,1] and ``prec`` finite ≥ 0. Anything else is a bug in the message,
not a modelling choice.

For each component (gDNA / RNA-pos / RNA-neg) and direction (forward α / backward β) scans every ACTIVE message
(prec > 0): counts mode<0, mode>1, NaN/inf; reports the ranges; dumps examples. Also tallies the gDNA-emission
gate (how much crossing mass is suppressed at strand-discontinuous boundaries).

    OMP_NUM_THREADS=1 python scripts/debug/fb_message_audit.py [cond ...]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

DEFAULT = ["gdna_gdna300_ss_0.50_nrna_none_capture_on", "gdna_none_ss_0.50_nrna_rnd_capture_off"]
_EPS = 1e-9


def capture(cond):
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    calmod = importlib.import_module("rigel.calibration.calibrate")
    orig = calmod.node_sweep
    cap = {}
    calmod.node_sweep = lambda *a, **k: orig(*a, _capture=cap, **k)
    try:
        calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
                  gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    finally:
        calmod.node_sweep = orig
    return cap


def scan_component(name, mode, prec):
    active = prec > _EPS
    n = int(active.sum())
    if n == 0:
        print(f"    {name:>10}: (no active messages)")
        return
    m = mode[active]
    bad_lo = int((m < -1e-6).sum())
    bad_hi = int((m > 1 + 1e-6).sum())
    nan = int((~np.isfinite(mode) | ~np.isfinite(prec))[active].sum())
    print(f"    {name:>10}: n={n:>5}  mode[min={m.min():+.3f} max={m.max():+.3f}]  "
          f"mode<0: {bad_lo:>4}  mode>1: {bad_hi:>4}  NaN/inf: {nan:>3}  "
          f"prec[max={prec[active].max():.1f}]")


def main():
    for cond in (sys.argv[1:] or DEFAULT):
        cap = capture(cond)
        a, b = cap["a_fwd"], cap["b_bwd"]
        print(f"\n=== {cond} : message defect scan ===")
        for dname, msgs in (("forward α", a), ("backward β", b)):
            print(f"  {dname}:")
            scan_component("gDNA", msgs[0], msgs[1])
            scan_component("RNA-pos", msgs[2], msgs[3])
            scan_component("RNA-neg", msgs[4], msgs[5])

        # gDNA-emission gate tally: crossing mass that CANNOT emit gDNA because the boundary is
        # strand-discontinuous (solvable=False), vs. the mass that can.
        solv = np.asarray(cap["solvable"], bool)
        ml, mr = cap["mass_l"], cap["mass_r"]
        cross = ml + mr
        emit_ok = solv & (cross > _EPS)
        suppressed = (~solv) & (cross > _EPS)
        print("  gDNA-emission gate (strand-criterion `solvable`):")
        print(f"    crossing mass that CAN emit gDNA  (solvable)   = {cross[emit_ok].sum():>12,.0f}")
        print(f"    crossing mass SUPPRESSED (not solvable)        = {cross[suppressed].sum():>12,.0f}"
              f"   ({int(suppressed.sum())} nodes)")


if __name__ == "__main__":
    main()
