"""Verify pass-1 prior-free (FLOOR_OFF) robustness: FP (gdna_none), unstranded (ss0.50), low-count stability.

For each condition: run production calibrate() baseline vs FLOOR_OFF, report contained gDNA error vs the
by-origin oracle, the gdna_none FALSE-POSITIVE gDNA (must stay ~0), and low-count node stability (finite,
no blow-up). The floor's stated job was "keep low/zero-gDNA nodes finite" — this checks removing it is safe.

    OMP_NUM_THREADS=1 python scripts/debug/pass1_priorfree_robustness.py
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402

CONDS = [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",   # baseline (stranded, gDNA, capture)
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",   # UNSTRANDED + gDNA + capture (the hard case)
    "gdna_none_ss_0.99_nrna_none_capture_on",       # FP test, stranded
    "gdna_none_ss_0.50_nrna_none_capture_on",       # FP test, unstranded (hardest FP)
]

def run(blob, ra, env):
    for k in ["RIGEL_DBG_FLOOR_OFF"]:
        os.environ.pop(k, None)
    os.environ.update(env)
    from rigel.calibration.calibrate import calibrate
    from rigel.config import CalibrationConfig
    return calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
                     gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())

for cond in CONDS:
    print("=" * 78)
    print(cond)
    try:
        index, blob = build_or_load_cache(cond, False)
    except Exception as e:
        print(f"  cache build FAILED: {e}"); continue
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    sg = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    sr = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
    go = np.asarray(sg.contained.mass_unspliced, float)
    ro = np.asarray(sr.contained.mass_unspliced, float)
    sf = CalibrationSubstrate.from_payload(blob["payload_full"], ra)
    mass = np.asarray(sf.contained.mass_unspliced, float)
    lowc = mass < 20.0  # low-count region nodes (stability target)
    for name, env in [("baseline", {}), ("FLOOR_OFF", {"RIGEL_DBG_FLOOR_OFF": "1"})]:
        c = run(blob, ra, env)
        g = np.asarray(c.mass_gdna_contained, float)
        r = np.asarray(c.mass_rna_contained, float)
        nfin = int((~np.isfinite(g)).sum() + (~np.isfinite(r)).sum())
        cont_err = (g - go).sum()
        fp = g.sum()  # for gdna_none this is pure false-positive gDNA
        # low-count node stability: f_g spread + any extreme
        fg_low = g[lowc] / (g[lowc] + r[lowc] + 1e-9)
        line = (f"  {name:10s} prod_gDNA={g.sum():>12,.0f}  oracle_gDNA={go.sum():>12,.0f}  "
                f"cont_err={cont_err:>+12,.0f}  nonfinite={nfin}")
        if "none" in cond:
            line += f"  [FP gDNA={fp:,.0f}]"
        line += f"  | lowcount(n={int(lowc.sum())}) mean_fg={np.nanmean(fg_low) if lowc.any() else float('nan'):.3f}"
        print(line)
