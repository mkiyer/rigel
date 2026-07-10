"""Diagnose the calibration's per-node gDNA-vs-RNA DECONVOLUTION accuracy against the sim oracle.

The deconvolution splits each node's unspliced mass into gDNA vs RNA. This dumps, per region, the OBSERVED
(calibrated) vs TRUE (oracle) gDNA CONTAINED mass + gDNA FRACTION, classified by signature (exon/intron/
intergenic) and by true enrichment, so we can localise WHERE the deconvolution under-calls captured-exon
gDNA (the ~2× under-call blocking the enriched-mode reference). True masses come from the oracle read-name
origins (mrna/nrna/gdna + true span); observed from calibrate().

    OMP_NUM_THREADS=1 python scripts/debug/calib_gdna_accuracy.py <condition_dir> [out.tsv]
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import numpy as np, pandas as pd, pysam

from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.calibrate import calibrate
from rigel.calibration.signature import coarse_type_from_signature
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

cond = Path(sys.argv[1]); suite = cond.parent
out_tsv = Path(sys.argv[2]) if len(sys.argv) > 2 else None
idx = TranscriptIndex.load(suite / "rigel_index")
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
starts = np.asarray(ra.start, np.int64); ends = np.asarray(ra.end, np.int64)
ref_off = np.asarray(ra.ref_offsets, np.int64); name_to_id = idx.ref_name_to_id
n = ra.n_regions

# TRUE per-region contained mass by origin (fragment span from the read name; contained = within one region)
true_g = np.zeros(n); true_r = np.zeros(n)
with pysam.AlignmentFile(str(cond / "sim_oracle.bam"), "rb") as f:
    default_ref = f.references[0] if f.references else None
    seen: set[str] = set()
    for read in f:
        q = read.query_name
        if q in seen:
            continue
        seen.add(q)
        o = parse_origin(q)
        if o.start is None:
            continue
        rid = name_to_id.get(str(o.ref if o.ref is not None else default_ref))
        if rid is None:
            continue
        lo0, hi0 = int(ref_off[rid]), int(ref_off[rid + 1])
        a, b = int(o.start), int(o.end)
        r = lo0 + int(np.searchsorted(ends[lo0:hi0], a, side="right"))
        if r < hi0 and starts[r] <= a and b <= ends[r]:  # contained
            (true_g if o.kind == "gdna" else true_r)[r] += 1.0

# OBSERVED (deconvolved) per-region contained mass
st, sm, flm, buf, pl = scan_and_buffer(str(cond / "sim_oracle.bam"), idx, BamScanConfig())
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
cal = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, CalibrationConfig())
obs_g = np.asarray(cal.mass_gdna_contained, float)
obs_r = np.asarray(cal.mass_rna_contained, float)
eff_len = np.maximum(np.asarray(cal.gdna_region_eff_len, float), 1e-9)

sig = idx.region_df["signature"].to_numpy()
cls = np.array(["ig", "intron", "exon"])[[coarse_type_from_signature(int(s)) for s in sig]]
df = pd.DataFrame(dict(
    cls=cls, S=eff_len,
    true_g=true_g, true_r=true_r, obs_g=obs_g, obs_r=obs_r,
    true_gfrac=true_g / np.maximum(true_g + true_r, 1e-9),
    obs_gfrac=obs_g / np.maximum(obs_g + obs_r, 1e-9),
    true_dens=true_g / eff_len, obs_dens=obs_g / eff_len,
))
df["total"] = df.true_g + df.true_r
active = df[df.total >= 20]  # nodes with enough fragments to judge

print(f"=== {cond.name} ===")
print(f"totals: obs_g={obs_g.sum():,.0f} true_g={true_g.sum():,.0f}  obs_r={obs_r.sum():,.0f} true_r={true_r.sum():,.0f}")
print(f"\ngDNA-FRACTION error (obs_gfrac - true_gfrac; NEGATIVE = deconv UNDER-calls gDNA), by class"
      f" [mass-weighted over active nodes]:")
for c in ["exon", "intron", "ig"]:
    s = active[active.cls == c]
    if len(s):
        wgfrac_err = np.average(s.obs_gfrac - s.true_gfrac, weights=s.total)
        print(f"  {c:7s} n={len(s):4d}  Δgfrac={wgfrac_err:+.3f}  "
              f"true_gfrac={np.average(s.true_gfrac,weights=s.total):.3f} obs_gfrac={np.average(s.obs_gfrac,weights=s.total):.3f}  "
              f"true_gmass={s.true_g.sum():,.0f} obs_gmass={s.obs_g.sum():,.0f} ({s.obs_g.sum()/max(s.true_g.sum(),1):.2f}×)")
# EXON nodes split by whether EXPRESSED (has RNA) vs not — the key hypothesis (RNA masks gDNA)
ex = active[active.cls == "exon"]
for lab, s in [("exon EXPRESSED (true_r>0)", ex[ex.true_r > 0]), ("exon SILENT (true_r≈0)", ex[ex.true_r == 0])]:
    if len(s):
        print(f"  {lab:26s} n={len(s):4d}  Δgfrac={np.average(s.obs_gfrac-s.true_gfrac,weights=s.total):+.3f}  "
              f"obs/true gmass={s.obs_g.sum()/max(s.true_g.sum(),1):.2f}×")
# by enrichment quartile (true gDNA density) among expressed exons
exx = ex[ex.true_r > 0].copy()
if len(exx) >= 8:
    exx["q"] = pd.qcut(exx.true_dens, 4, labels=["Q1(low)", "Q2", "Q3", "Q4(hi)"], duplicates="drop")
    print("  expressed-exon under-call by TRUE-gDNA-density quartile:")
    for q, s in exx.groupby("q", observed=True):
        print(f"    {q:8s} true_dens≈{s.true_dens.median():6.2f}  Δgfrac={np.average(s.obs_gfrac-s.true_gfrac,weights=s.total):+.3f}  obs/true={s.obs_g.sum()/max(s.true_g.sum(),1):.2f}×")
if out_tsv:
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nper-region table -> {out_tsv}")
