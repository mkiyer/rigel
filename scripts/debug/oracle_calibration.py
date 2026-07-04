"""Build an ORACLE calibration (perfect per-node gDNA/RNA deconvolution) from a sim oracle BAM.

Each fragment's true origin (mrna/nrna/gdna) + genomic span is in its read name, so we deposit the TRUE
per-node masses directly — contained fragments to their region, crossing fragments split across the
boundary sides they span — replacing calibration's deconvolved split with ground truth. Saves the mass
arrays to ``rigel_out/oracle_calib.npz``; calibrate() loads + overrides them when RIGEL_ORACLE_CALIB is set,
so the whole pipeline runs on perfect calibration and any residual leak is PURELY the eff-len formula.

    OMP_NUM_THREADS=1 python scripts/debug/oracle_calibration.py <condition_dir>
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import numpy as np, pysam

from rigel.index import TranscriptIndex
from rigel.calibration.region_arrays import RegionArrays
from rigel.sim.read_name import parse_origin

cond = Path(sys.argv[1]); suite = cond.parent
idx = TranscriptIndex.load(suite / "rigel_index")
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
starts = np.asarray(ra.start, np.int64); ends = np.asarray(ra.end, np.int64)
ref_off = np.asarray(ra.ref_offsets, np.int64); ref_id = np.asarray(ra.ref_id)
name_to_id = idx.ref_name_to_id; n = ra.n_regions

g_c = np.zeros(n); r_c = np.zeros(n)  # contained gDNA / RNA
g_l = np.zeros(n); r_l = np.zeros(n)  # left-side (right side of left boundary)
g_r = np.zeros(n); r_rt = np.zeros(n)  # right-side (left side of right boundary)
r_spl = np.zeros(n)  # spliced RNA (contained region of the spliced fragment's dominant block)

with pysam.AlignmentFile(str(cond / "sim_oracle.bam"), "rb") as f:
    default_ref = f.references[0] if f.references else None
    seen: set[str] = set()
    is_spliced_by_q: dict[str, bool] = {}
    for read in f:
        q = read.query_name
        spliced = read.cigartuples is not None and any(op == 3 for op, _ in read.cigartuples)  # N = splice
        is_spliced_by_q[q] = is_spliced_by_q.get(q, False) or spliced
        if q in seen:
            continue
        seen.add(q)
    # second pass would re-read; instead reuse: re-open to deposit with the spliced flag known
with pysam.AlignmentFile(str(cond / "sim_oracle.bam"), "rb") as f:
    default_ref = f.references[0] if f.references else None
    seen = set()
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
        is_g = o.kind == "gdna"
        # regions overlapping [a,b): first with end>a through last with start<b
        lo = lo0 + int(np.searchsorted(ends[lo0:hi0], a, side="right"))
        hi = lo0 + int(np.searchsorted(starts[lo0:hi0], b, side="left"))
        if hi <= lo:
            continue
        if hi - lo == 1:  # CONTAINED in one region
            (g_c if is_g else r_c)[lo] += 1.0
            if not is_g and is_spliced_by_q.get(q, False):
                r_spl[lo] += 1.0
        else:  # CROSSING: split 1.0 equally across the (hi-1-lo) internal boundaries it spans, 0.5/side
            nb = hi - 1 - lo
            w = 0.5 / nb
            for rr in range(lo, hi - 1):
                if ref_id[rr] != ref_id[rr + 1]:
                    continue
                (g_r if is_g else r_rt)[rr] += w
                (g_l if is_g else r_l)[rr + 1] += w
            if not is_g and is_spliced_by_q.get(q, False):
                r_spl[lo] += 1.0  # attribute spliced mass to its first region

out = cond / "rigel_out" / "oracle_calib.npz"
out.parent.mkdir(exist_ok=True)
np.savez(out, mass_gdna_contained=g_c, mass_rna_contained=r_c, mass_gdna_left=g_l, mass_rna_left=r_l,
         mass_gdna_right=g_r, mass_rna_right=r_rt, mass_rna_spliced=r_spl)
print(f"gDNA contained total={g_c.sum():.0f} crossing={g_l.sum()+g_r.sum():.0f}  "
      f"RNA contained={r_c.sum():.0f} crossing={r_l.sum()+r_rt.sum():.0f} spliced={r_spl.sum():.0f}")
print(f"saved {out}")
