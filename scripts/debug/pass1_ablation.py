"""Ablate pass-1 crushers (floor/global prior + imputation messages) on the AMBIG-exon nodes.

Runs production calibrate() under each toggle combo and reports the mass-weighted gDNA fraction over
the 92 AMBIG-exon region nodes (oracle ~0.757) plus a few representative nodes. Isolates which pass-1
component crushes the strand signal.

    OMP_NUM_THREADS=1 python scripts/debug/pass1_ablation.py
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, subprocess
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature  # noqa: E402

COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
index, blob = build_or_load_cache(COND, False)
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
sig = index.region_df["signature"].to_numpy()
scls = np.array([coarse_strand_from_signature(int(s)) for s in sig])
tcls = np.array([coarse_type_from_signature(int(s)) for s in sig])
ambig_exon = np.array([r for r in range(len(sig)) if scls[r] == 3 and tcls[r] == 2])
sub_g = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
sub_r = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
g_or = np.asarray(sub_g.contained.mass_unspliced, float)
r_or = np.asarray(sub_r.contained.mass_unspliced, float)

def run():
    from rigel.calibration.calibrate import calibrate
    from rigel.config import CalibrationConfig
    return calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
                     gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())

combos = [
    ("baseline",         {}),
    ("GLOBAL_OFF",       {"RIGEL_DBG_GLOBAL_OFF": "1"}),
    ("MSG_OFF",          {"RIGEL_DBG_MSG_OFF": "1"}),
    ("GLOBAL+MSG_OFF",   {"RIGEL_DBG_GLOBAL_OFF": "1", "RIGEL_DBG_MSG_OFF": "1"}),
]
ae = ambig_exon
M = np.asarray(sub_g.contained.mass_unspliced + sub_r.contained.mass_unspliced, float)[ae]
orc = (g_or[ae]) / (g_or[ae] + r_or[ae] + 1e-9)
orc_mw = (M * orc).sum() / M.sum()
print(f"AMBIG-exon nodes: {len(ae)}   oracle mass-weighted f_g = {orc_mw:.3f}   (mass {M.sum():,.0f})")
print(f"{'combo':>16}  {'mw_fg':>7}  {'vs_orc':>7}  {'gDNA_mass':>11}  | nodes 242/224/672/1083")
for name, env in combos:
    for k in ("RIGEL_DBG_GLOBAL_OFF", "RIGEL_DBG_MSG_OFF"):
        os.environ.pop(k, None)
    os.environ.update(env)
    cal = run()
    gpr = np.asarray(cal.mass_gdna_contained, float)
    rpr = np.asarray(cal.mass_rna_contained, float)
    fg = gpr[ae] / (gpr[ae] + rpr[ae] + 1e-9)
    mw = (M * fg).sum() / M.sum()
    d = dict(zip(ae.tolist(), fg.tolist()))
    reps = "  ".join(f"{d.get(r, float('nan')):.3f}" for r in [242, 224, 672, 1083])
    print(f"{name:>16}  {mw:>7.3f}  {mw-orc_mw:>+7.3f}  {(M*fg).sum():>11,.0f}  | {reps}")
