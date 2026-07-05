"""Diagnose whether the calibration deconvolution SMOOTHS the per-node gDNA density (mass-conserving but
peak-lowering) vs the sim oracle. gDNA read-name coords ARE genomic, so true per-region gDNA is exact.

Compares OBSERVED (calibrated) vs TRUE per-region gDNA contained mass + density, and quantifies smoothing:
per-node density correlation, the peak ratio (obs/true on the hottest true nodes), and whether mass flows
DOWN the density gradient (obs<true on enriched nodes, obs>true on depleted neighbours). This is the
mechanism behind the enriched-mode reference reading low (obs peak ~55 vs true ~103).

    OMP_NUM_THREADS=1 python scripts/debug/calib_gdna_smoothing.py <condition_dir>
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
idx = TranscriptIndex.load(suite / "rigel_index")
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
starts = np.asarray(ra.start, np.int64); ends = np.asarray(ra.end, np.int64)
ref_off = np.asarray(ra.ref_offsets, np.int64); name_to_id = idx.ref_name_to_id
n = ra.n_regions

true_g = np.zeros(n)  # TRUE gDNA contained (genomic read-name span, within one region)
with pysam.AlignmentFile(str(cond / "sim_oracle.bam"), "rb") as f:
    default_ref = f.references[0] if f.references else None
    seen: set[str] = set()
    for read in f:
        q = read.query_name
        if q in seen:
            continue
        seen.add(q)
        o = parse_origin(q)
        if o.kind != "gdna" or o.start is None:
            continue
        rid = name_to_id.get(str(o.ref if o.ref is not None else default_ref))
        if rid is None:
            continue
        lo0, hi0 = int(ref_off[rid]), int(ref_off[rid + 1])
        a, b = int(o.start), int(o.end)
        r = lo0 + int(np.searchsorted(ends[lo0:hi0], a, side="right"))
        if r < hi0 and starts[r] <= a and b <= ends[r]:
            true_g[r] += 1.0

st, sm, flm, buf, pl = scan_and_buffer(str(cond / "sim_oracle.bam"), idx, BamScanConfig())
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
cal = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, CalibrationConfig())
obs_g = np.asarray(cal.mass_gdna_contained, float)
S = np.maximum(np.asarray(cal.gdna_region_eff_len, float), 1e-9)
sig = idx.region_df["signature"].to_numpy()
cls = np.array(["ig", "intron", "exon"])[[coarse_type_from_signature(int(s)) for s in sig]]

td, od = true_g / S, obs_g / S  # densities
m = (true_g + obs_g) > 5  # nodes with some gDNA
print(f"=== {cond.name} ===")
print(f"total gDNA contained: obs={obs_g.sum():,.0f} true={true_g.sum():,.0f} ({obs_g.sum()/max(true_g.sum(),1):.3f}×)")
print(f"per-node density: corr(obs,true)={np.corrcoef(od[m],td[m])[0,1]:.3f}  max: obs={od.max():.1f} true={td.max():.1f} ({od.max()/max(td.max(),1e-9):.2f}×)")
# peak smoothing: on the hottest TRUE-density nodes, obs/true
order = np.argsort(td)[::-1]
for k in (5, 20, 50):
    top = order[:k]
    print(f"  top-{k:2d} true-density nodes: mean true_dens={td[top].mean():6.1f} obs_dens={od[top].mean():6.1f}  obs/true={od[top].mean()/max(td[top].mean(),1e-9):.2f}  mass obs/true={obs_g[top].sum()/max(true_g[top].sum(),1):.2f}")
# mass flow down the gradient: split nodes into enriched (top) vs depleted; where did obs mass go vs true?
enr = td > np.median(td[m])
print(f"  ENRICHED nodes (td>median): obs_gmass/true_gmass={obs_g[enr].sum()/max(true_g[enr].sum(),1):.2f}× "
      f"(<1 ⇒ peak mass LOST)")
print(f"  DEPLETED nodes (td<=median): obs_gmass/true_gmass={obs_g[~enr & m].sum()/max(true_g[~enr & m].sum(),1):.2f}× "
      f"(>1 ⇒ received smoothed mass)")
# by class: the peak-density under-read where it matters (exons)
for c in ["exon", "intron", "ig"]:
    cm = m & (cls == c)
    if cm.sum() >= 5:
        o5 = np.argsort(np.where(cm, td, -1))[::-1][:max(5, cm.sum() // 10)]
        print(f"  {c:7s}: peak obs/true dens (top-decile)={od[o5].mean()/max(td[o5].mean(),1e-9):.2f}  "
              f"corr={np.corrcoef(od[cm],td[cm])[0,1]:.3f}")
