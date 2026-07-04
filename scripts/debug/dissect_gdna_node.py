"""Rank calibration nodes by ABSOLUTE gDNA deconvolution error (|obs_g - true_g|) and dump the top culprits
with the deconvolution's actual inputs — for single-node dissection of the gDNA↔RNA leak.

Per region (contained/unspliced): raw = total unspliced contained mass (what the deconv splits); true_g =
oracle gDNA (genomic read-name); true_r = raw - true_g (unspliced RNA); obs_g = calibrated gDNA. The strand
tilt (n_pos/n_neg) is the only intrinsic gDNA/RNA signal. Ranks by |obs_g - true_g| and prints the worst
nodes + the #1's neighbourhood (the cross-node imputation context).

    OMP_NUM_THREADS=1 python scripts/debug/dissect_gdna_node.py <condition_dir> [topN]
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import numpy as np, pandas as pd, pysam

from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.calibrate import calibrate
from rigel.calibration.signature import coarse_type_from_signature
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

cond = Path(sys.argv[1]); suite = cond.parent
topN = int(sys.argv[2]) if len(sys.argv) > 2 else 15
idx = TranscriptIndex.load(suite / "rigel_index")
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
starts = np.asarray(ra.start, np.int64); ends = np.asarray(ra.end, np.int64)
roff = np.asarray(ra.ref_offsets, np.int64); refid = np.asarray(ra.ref_id)
n = ra.n_regions

true_g = np.zeros(n)
with pysam.AlignmentFile(str(cond / "sim_oracle.bam"), "rb") as f:
    dref = f.references[0]; seen: set[str] = set()
    for r in f:
        q = r.query_name
        if q in seen:
            continue
        seen.add(q); o = parse_origin(q)
        if o.kind != "gdna" or o.start is None:
            continue
        rid = idx.ref_name_to_id.get(str(o.ref if o.ref is not None else dref))
        if rid is None:
            continue
        lo0, hi0 = int(roff[rid]), int(roff[rid + 1]); a, b = int(o.start), int(o.end)
        rr = lo0 + int(np.searchsorted(ends[lo0:hi0], a, side="right"))
        if rr < hi0 and starts[rr] <= a and b <= ends[rr]:
            true_g[rr] += 1.0

st, sm, flm, buf, pl = scan_and_buffer(str(cond / "sim_oracle.bam"), idx, BamScanConfig())
sub = CalibrationSubstrate.from_payload(pl, ra)
raw = np.asarray(sub.contained.mass_unspliced, float)         # total unspliced contained (gDNA+RNA)
npos = np.asarray(sub.contained.n_unspliced_pos, float)
nneg = np.asarray(sub.contained.n_unspliced_neg, float)
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
cal = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, CalibrationConfig())
obs_g = np.asarray(cal.mass_gdna_contained, float)
S = np.maximum(np.asarray(cal.gdna_region_eff_len, float), 1e-9)
sig = idx.region_df["signature"].to_numpy()
cls = np.array(["ig", "intron", "exon"])[[coarse_type_from_signature(int(s)) for s in sig]]

true_g = np.minimum(true_g, raw)  # true gDNA can't exceed the contained total
err = obs_g - true_g              # negative = deconv UNDER-calls gDNA (the leak of gDNA→RNA)
df = pd.DataFrame(dict(
    region=np.arange(n), ref=refid, start=starts, end=ends, cls=cls, S=S,
    raw=raw, true_g=true_g, obs_g=obs_g, err=err,
    true_gfrac=true_g / np.maximum(raw, 1e-9), obs_gfrac=obs_g / np.maximum(raw, 1e-9),
    sense_frac=npos / np.maximum(npos + nneg, 1e-9), npos=npos, nneg=nneg,
    true_dens=true_g / S, obs_dens=obs_g / S,
))
print(f"=== {cond.name} ===  total |err|={np.abs(err).sum():,.0f}  (Σ obs_g={obs_g.sum():,.0f} true_g={true_g.sum():,.0f})")
print(f"net err by class:  " + "  ".join(f"{c}:{df[df.cls==c].err.sum():+,.0f}" for c in ["exon","intron","ig"]))
top = df.reindex(df.err.abs().sort_values(ascending=False).index).head(topN)
pd.set_option("display.width", 240)
print(f"\nTOP {topN} nodes by |gDNA error|:")
print(top[["region","cls","start","end","S","raw","true_g","obs_g","err","true_gfrac","obs_gfrac","sense_frac"]]
      .to_string(index=False, float_format=lambda x: f"{x:,.1f}"))
# dissect #1 + its genomic neighbours (the imputation context)
w = int(top.iloc[0].region)
print(f"\n=== #1 CULPRIT: region {w} ({cls[w]}) [{starts[w]}-{ends[w]}) ref={refid[w]} ===")
print(f"  raw={raw[w]:,.0f}  true_g={true_g[w]:,.0f} (gfrac {true_g[w]/max(raw[w],1):.3f})  "
      f"obs_g={obs_g[w]:,.0f} (gfrac {obs_g[w]/max(raw[w],1):.3f})  ERR={err[w]:+,.0f}")
print(f"  strand: n_pos={npos[w]:,.0f} n_neg={nneg[w]:,.0f} sense_frac={npos[w]/max(npos[w]+nneg[w],1):.3f}  "
      f"(gDNA→0.5 sense; mature→rna_sense_frac={cal.rna_sense_frac:.2f})")
print(f"  eff_len S={S[w]:.0f}  true_dens={true_g[w]/S[w]:.2f} obs_dens={obs_g[w]/S[w]:.2f}  "
      f"gdna_density_global={cal.gdna_density_global:.4g}")
print(f"  NEIGHBOURS (imputation context):")
for d in (-2, -1, 1, 2):
    j = w + d
    if 0 <= j < n and refid[j] == refid[w]:
        print(f"    r{j:+d} ({cls[j]:6s}) raw={raw[j]:8,.0f} true_g={true_g[j]:8,.0f} obs_dens={obs_g[j]/S[j]:7.2f} true_dens={true_g[j]/S[j]:7.2f} sense={npos[j]/max(npos[j]+nneg[j],1):.2f}")
