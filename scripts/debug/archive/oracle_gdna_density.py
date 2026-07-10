"""Precompute the TRUE per-region gDNA contained density from a sim oracle BAM (for the ρ* oracle study).

Each gDNA fragment's true genomic span is encoded in its read name (``gdna:[ref:]start-end:...``), so we
parse origins directly (no alignment reconstruction). A fragment CONTAINED in a single region increments
that region's true-gDNA count; the true density is count / gdna_region_eff_len (the same effective support
the observed ρ* uses). Saves ``oracle_gdna_contained.npy`` (per-region true gDNA contained mass) next to the
condition, consumed by priors' RIGEL_RHOSTAR=oracle_* modes.

    OMP_NUM_THREADS=1 python scripts/debug/oracle_gdna_density.py <condition_dir>
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import numpy as np, pysam

from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.calibrate import calibrate
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

cond = Path(sys.argv[1])
suite = cond.parent
idx = TranscriptIndex.load(suite / "rigel_index")
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
starts = np.asarray(ra.start, np.int64); ends = np.asarray(ra.end, np.int64)
ref_off = np.asarray(ra.ref_offsets, np.int64)
name_to_id = idx.ref_name_to_id
n = ra.n_regions
true_contained = np.zeros(n, np.float64)

bam = cond / "sim_oracle.bam"
seen: set[str] = set()
with pysam.AlignmentFile(str(bam), "rb") as f:
    default_ref = f.references[0] if f.references else None
    for read in f:
        q = read.query_name
        if q in seen:
            continue
        seen.add(q)
        o = parse_origin(q)
        if o.kind != "gdna" or o.start is None:
            continue
        ref = o.ref if o.ref is not None else default_ref
        rid = name_to_id.get(str(ref))
        if rid is None:
            continue
        lo0, hi0 = int(ref_off[rid]), int(ref_off[rid + 1])
        a, b = int(o.start), int(o.end)
        # region containing a: first region with end > a; contained iff b <= that region's end
        r = lo0 + int(np.searchsorted(ends[lo0:hi0], a, side="right"))
        if r < hi0 and starts[r] <= a and b <= ends[r]:
            true_contained[r] += 1.0

# true density needs gdna_region_eff_len — get it from a calibration run (same support the observed ρ* uses)
st, sm, flm, buf, pl = scan_and_buffer(str(bam), idx, BamScanConfig())
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
cal = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, CalibrationConfig())
eff_len = np.maximum(np.asarray(cal.gdna_region_eff_len, np.float64), 1e-9)
true_density = true_contained / eff_len

out = cond / "rigel_out" / "oracle_gdna_contained.npy"
out.parent.mkdir(exist_ok=True)
np.save(out, true_contained)
obs = np.asarray(cal.mass_gdna_contained, np.float64)
m = (true_contained > 0) | (obs > 0)
print(f"regions with gDNA: true={int((true_contained>0).sum())} obs={int((obs>0).sum())}")
print(f"true contained total={true_contained.sum():.0f}  obs contained total={obs.sum():.0f}")
print(f"true density  (nonzero): median={np.median(true_density[true_contained>0]):.3f} max={true_density.max():.3f}")
print(f"obs  density  (nonzero): median={np.median((obs/eff_len)[obs>0]):.3f} max={(obs/eff_len).max():.3f}")
print(f"saved {out}")
