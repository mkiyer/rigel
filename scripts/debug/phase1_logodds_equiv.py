"""Phase 1 verification: the log-odds λ-grid LOCAL solve vs the 2-simplex lattice LOCAL solve, on
SINGLE-STRAND solvable nodes. Identical ψ inputs (same global on σ(λ), same strand/spliced/Jeffreys);
only the grid differs (linear lattice ↔ log-odds). Bit-comparable ⇒ the 1-D reduction + the log-odds
grid lose no accuracy. Reports per-node |Δf_g| + the aggregate gDNA-mass agreement.

    OMP_NUM_THREADS=1 python scripts/debug/phase1_logodds_equiv.py [cond1 ...]
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

CONDS = [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna300_ss_0.50_nrna_rnd_capture_on",
    "gdna_none_ss_0.50_nrna_none_capture_off",
]


def run(cond: str):
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    calmod = importlib.import_module("rigel.calibration.calibrate")
    orig = calmod.node_sweep
    cap = {}
    calmod.node_sweep = lambda *a, **k: orig(*a, _capture=cap, **k)
    ng = int(os.environ.get("RIGEL_NGRID", "60"))
    cfg = CalibrationConfig(sweep_n_grid=ng)
    try:
        calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
                  gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=cfg)
    finally:
        calmod.node_sweep = orig

    fg_lat = np.asarray(cap["fg_loc"], float)
    fg_lo = np.asarray(cap["fg_loc_logodds"], float)
    fg_lo_amb = np.asarray(cap["fg_loc_logodds_ambig"], float)
    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    solv = np.asarray(cap["solvable"], bool)
    mass = np.asarray(cap["mass_global"], float)
    print(f"\n{cond}")

    def _report(label, mask, fg_lo_arr):
        d = np.abs(fg_lo_arr[mask] - fg_lat[mask])
        al = float(np.sum(fg_lat[mask] * mass[mask]))
        ao = float(np.sum(fg_lo_arr[mask] * mass[mask]))
        print(f"  [{label}] solvable n = {int(mask.sum()):,}")
        print(f"    |Δf_g|: median={np.median(d):.4f} mean={d.mean():.4f} "
              f"p95={np.percentile(d,95):.4f} p99={np.percentile(d,99):.4f} max={d.max():.4f}")
        print(f"    aggregate gDNA mass: lattice={al:,.0f} log-odds={ao:,.0f} "
              f"Δ={100*(ao-al)/max(abs(al),1.0):+.3f}%")
        return d, al, ao

    ss = (fp ^ fn) & solv          # single-strand solvable (Phase 1)
    amb = (fp & fn) & solv         # AMBIG solvable (Phase 2)
    d_ss, al_ss, ao_ss = _report("single-strand", ss, fg_lo)
    d_amb = al_amb = ao_amb = None
    if bool(amb.any()):
        d_amb, al_amb, ao_amb = _report("AMBIG (2-D)", amb, fg_lo_amb)
        ai = np.where(amb)[0]
        worst = ai[np.argsort(-np.abs(fg_lo_amb[amb] - fg_lat[amb]))[:5]]
        print("    worst AMBIG |Δf_g| (lat → lo, mass):")
        for n in worst:
            print(f"      node {n}: {fg_lat[n]:.4f} → {fg_lo_amb[n]:.4f}  "
                  f"Δ={abs(fg_lo_amb[n]-fg_lat[n]):.4f}  mass={mass[n]:,.0f}")
    return d_ss, al_ss, ao_ss, d_amb, al_amb, ao_amb


def main():
    conds = sys.argv[1:] or CONDS
    mx_ss = mx_amb = 0.0
    for c in conds:
        try:
            d_ss, _, _, d_amb, _, _ = run(c)
            mx_ss = max(mx_ss, float(d_ss.max()))
            if d_amb is not None and d_amb.size:
                mx_amb = max(mx_amb, float(d_amb.max()))
        except Exception as e:  # noqa: BLE001
            print(f"\n!! {c}: {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
    print(f"\n{'='*70}\nOVERALL max |Δf_g|:  single-strand={mx_ss:.4f}   AMBIG={mx_amb:.4f}")


if __name__ == "__main__":
    main()
