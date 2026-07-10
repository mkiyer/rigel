"""Verify the IPR theory: G/eff_len = mass-weighted MEAN gDNA density, and under capture the exon
(competition) density exceeds that mean because introns/crossings retain baseline gDNA that dilutes it.

Scans+calibrates one sweep condition, then dumps the gDNA node densities by class (exon-contained /
intron-contained / intergenic / crossing-seam) with mass + support + density, and compares:
  - the mass-weighted mean density  Σρ_n m_n / Σm_n   (== G/eff_len, the density the EM uses)
  - the exon-contained density       (the density where gDNA competes with mature RNA)
The ratio exon/mean is the predicted ~2x under-weighting.

    OMP_NUM_THREADS=1 python scripts/debug/ipr_theory_probe.py [condition_dir]
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import numpy as np

from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.calibrate import calibrate
from rigel.calibration.signature import coarse_type_from_signature
from rigel.splice import SpliceType

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/efflen_binding_sweep")
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna100_ss_0.99_nrna_none_capture_b200"

idx = TranscriptIndex.load(SUITE / "rigel_index")
bam = SUITE / COND / "sim_oracle.bam"
st, sm, flm, buf, pl = scan_and_buffer(str(bam), idx, BamScanConfig())
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
cal = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, CalibrationConfig())

sig = idx.region_df["signature"].to_numpy()
tcls = np.array([coarse_type_from_signature(int(s)) for s in sig])  # 0=intergenic 1=intron 2=exon
gc = np.asarray(cal.mass_gdna_contained, float)
Sc = np.maximum(np.asarray(cal.gdna_region_eff_len, float), 1e-9)   # contained support
# crossing-seam nodes (boundary sides): mass + support
gl = np.asarray(cal.mass_gdna_left, float); grr = np.asarray(cal.mass_gdna_right, float)
Sb = np.maximum(np.asarray(cal.gdna_boundary_len, float), 1e-9)

print(f"=== {COND} ===")
print(f"{'class':>16} {'n':>5} {'gDNA_mass':>12} {'support':>12} {'density(m/S)':>12} {'mass_share':>10}")
Gtot = gc.sum() + gl.sum() + grr.sum()
node_rho = []  # (rho_n, m_n) for the mass-weighted mean
def row(name, m, S):
    m = np.asarray(m, float); S = np.asarray(S, float)
    tot = m.sum()
    dens = tot / max(S[m > 0].sum(), 1e-9)  # class-aggregate density
    print(f"{name:>16} {int((m>0).sum()):>5} {tot:>12,.0f} {S[m>0].sum():>12,.0f} {dens:>12.3f} {100*tot/max(Gtot,1):>9.1f}%")
    # per-node densities for the mean
    ok = m > 0
    for mi, Si in zip(m[ok], S[ok]):
        node_rho.append((mi / Si, mi))
for name, tc in [("exon-contained", 2), ("intron-contained", 1), ("intergenic", 0)]:
    msk = tcls == tc
    row(name, np.where(msk, gc, 0.0), Sc)
# seams: attribute to their own class (crossing); support Sb applies per side
row("crossing-seam", np.concatenate([gl, grr]), np.concatenate([Sb, Sb]))

node_rho = np.array(node_rho)
mean_dens = float((node_rho[:, 0] * node_rho[:, 1]).sum() / node_rho[:, 1].sum())  # Σρm/Σm
exon_dens = gc[tcls == 2].sum() / max(Sc[(tcls == 2) & (gc > 0)].sum(), 1e-9)
print(f"\n  mass-weighted MEAN density Σρ_n·m_n/Σm_n  = {mean_dens:.3f}   (== G/eff_len, what the EM uses)")
print(f"  EXON-contained density (competition zone)  = {exon_dens:.3f}")
print(f"  under-weighting ratio exon/mean            = {exon_dens/max(mean_dens,1e-9):.2f}   (predicted ~2x)")
