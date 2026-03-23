#!/usr/bin/env python3
"""Investigate the 3.6% 'dominated' cases where MANE loses to a competitor.

Gap = exactly log(0.01) = -4.6052, which corresponds to EITHER:
  (a) 1 base of overhang (oh_log_pen * 1 = -4.6052), OR
  (b) 2 NM mismatches (mm_log_pen * 2 = 2 * -2.303 = -4.6052)

We need to figure out which scenario is happening and whether it's a bug.

Also investigate the bimodal tx_start distribution in dominated cases:
  - Group at tx_start ~ 0-200  (near the 5' end of MANE)
  - Group at tx_start ~ 600    (somewhere in the middle of MANE)
"""
import sys
sys.path.insert(0, "src")

import math
import gc
import numpy as np
import pysam
from collections import defaultdict, Counter

from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer, _setup_geometry_and_estimator, _score_fragments
from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig
from rigel.stats import PipelineStats

BASE = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
IDX  = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
MM2_BAM = f"{BASE}/sim_minimap2.bam"

print("Loading index...")
tx_index = TranscriptIndex.load(IDX)
tx_df = tx_index.t_df

t_info = {}
for _, r in tx_df.iterrows():
    t_info[int(r["t_index"])] = {
        "t_id": r["t_id"], "g_name": r.get("g_name", "?"),
        "length": int(r.get("length", 0)),
        "is_basic": bool(r.get("is_basic", False)),
        "is_mane": bool(r.get("is_mane", False)),
        "is_nrna": bool(r.get("is_synthetic_nrna", False)),
        "start": int(r.get("start", 0)), "end": int(r.get("end", 0)),
    }

mane_tidx   = next(ti for ti, inf in t_info.items() if inf["is_mane"] and inf["g_name"] == "GAPDH")
gapdh_mrna_set = {ti for ti, inf in t_info.items() if inf["g_name"] == "GAPDH" and not inf["is_nrna"]}

# Map frag_id → BAM qname (BAM is name-sorted, buffer processes in order)
print("Mapping frag_id → qname...")
frag_id_to_qname = {}
seen = set()
fid = 0
with pysam.AlignmentFile(MM2_BAM, "rb") as bam:
    for r in bam:
        if r.query_name not in seen:
            seen.add(r.query_name)
            frag_id_to_qname[fid] = r.query_name
            fid += 1

print(f"  Mapped {len(frag_id_to_qname)} qnames")

print("\nScanning + scoring MM2 BAM...")
scan_config = BamScanConfig(include_multimap=True)
pipeline_stats, strand_models, fl_models, buffer, _, _ = scan_and_buffer(
    MM2_BAM, tx_index, scan=scan_config
)
em_config = EMConfig()
scoring_config = FragmentScoringConfig()
geometry, estimator = _setup_geometry_and_estimator(tx_index, fl_models, em_config)
stats_obj = PipelineStats()
sf = _score_fragments(
    buffer, tx_index, strand_models, fl_models, stats_obj, estimator,
    scoring_config, log_every=200000, annotations=None
)
print(f"  EM units: {sf.n_units:,}")

offsets   = sf.offsets
t_indices = sf.t_indices
log_liks  = sf.log_liks
tx_starts = sf.tx_starts
tx_ends   = sf.tx_ends
frag_ids  = sf.frag_ids

# ---- Collect dominated + tied cases with full detail ----
TOL = 1e-4
OH_PEN  = math.log(0.01)   # -4.6052
NM_PEN  = math.log(0.1)    # -2.3026

dominated_units = []
tied_units   = []
best_units   = []

for unit_i in range(sf.n_units):
    s = int(offsets[unit_i])
    e = int(offsets[unit_i + 1])
    if e <= s:
        continue
    unit_t  = [int(t)  for t in t_indices[s:e]]
    unit_ll = [float(l) for l in log_liks[s:e]]
    unit_ts = [int(ts) for ts in tx_starts[s:e]]
    unit_te = [int(te) for te in tx_ends[s:e]]

    gapdh_idx = [j for j, t in enumerate(unit_t) if t in gapdh_mrna_set]
    if not gapdh_idx:
        continue
    if mane_tidx not in [unit_t[j] for j in gapdh_idx]:
        continue
    n_gapdh = len(gapdh_idx)
    if n_gapdh <= 1:
        continue

    g_t  = [unit_t[j]  for j in gapdh_idx]
    g_ll = [unit_ll[j] for j in gapdh_idx]
    g_ts = [unit_ts[j] for j in gapdh_idx]
    g_te = [unit_te[j] for j in gapdh_idx]

    mane_pos = g_t.index(mane_tidx)
    mane_ll  = g_ll[mane_pos]
    mane_ts_ = g_ts[mane_pos]
    mane_te_ = g_te[mane_pos]
    other_best = max(ll for j, ll in enumerate(g_ll) if j != mane_pos)
    gap = other_best - mane_ll
    fid = int(frag_ids[unit_i])

    rec = {
        "unit_i": unit_i, "fid": fid,
        "mane_ll": mane_ll, "other_best": other_best, "gap": gap,
        "mane_ts": mane_ts_, "mane_te": mane_te_,
        "cands": [(g_t[j], g_ll[j], g_ts[j], g_te[j]) for j in range(len(g_t))],
    }
    if abs(gap) <= TOL:
        tied_units.append(rec)
    elif gap > TOL:
        dominated_units.append(rec)
    else:
        best_units.append(rec)

print(f"\nCases: dominated={len(dominated_units)}, tied={len(tied_units)}, best={len(best_units)}")

# ---- Analyze dominated cases ----
dom_gaps = [r["gap"] for r in dominated_units]
print(f"\nDominated gap distribution:")
for g in sorted(set(round(x,4) for x in dom_gaps)):
    cnt = sum(1 for x in dom_gaps if abs(x-g) < 0.01)
    print(f"  gap={g:.4f}  count={cnt}")
    # Is it 1 overhang? 2 NM? or combination?
    oh_match = abs(g - abs(OH_PEN)) < 0.01  # 1 base overhang
    nm2_match = abs(g - 2 * abs(NM_PEN)) < 0.01  # 2 NM penalties
    print(f"    = 1 overhang base: {oh_match}")
    print(f"    = 2 NM penalties:  {nm2_match}")

# ---- Show detailed dominated cases ----
print(f"\n{'='*70}")
print("DOMINATED CASES: full candidate list + BAM qname")
print(f"{'='*70}")

for rec in dominated_units[:15]:
    qname = frag_id_to_qname.get(rec["fid"], "?")
    print(f"\n  unit={rec['unit_i']:6d} fid={rec['fid']:6d} "
          f"qname={qname}")
    print(f"  gap={rec['gap']:.4f}  mane_ts={rec['mane_ts']}  "
          f"mane_te={rec['mane_te']}")
    for ti, ll, ts, te in sorted(rec["cands"], key=lambda x: -x[1]):
        inf = t_info.get(ti, {})
        marker = " <<< MANE" if ti == mane_tidx else ""
        flen = te - ts
        print(f"    [{ti:6d}] {inf.get('t_id','?'):25s} ll={ll:9.4f} "
              f"tx=[{ts:4d},{te:4d}] flen={flen:3d} "
              f"gstart={inf.get('start',0):>10d}{marker}")

# ---- Key question: do dominated cases cluster at specific positions? ----
print(f"\n{'='*70}")
print("TX_START bimodal analysis for dominated cases")
print("Looking for two groups: near-5'-end and mid-transcript")
print(f"{'='*70}")

dom_ts = [r["mane_ts"] for r in dominated_units]
dom_ts_arr = np.array(dom_ts)
print(f"tx_start: mean={dom_ts_arr.mean():.0f} "
      f"p25={np.percentile(dom_ts_arr,25):.0f} "
      f"p50={np.median(dom_ts_arr):.0f} "
      f"p75={np.percentile(dom_ts_arr,75):.0f}")

bins = [0, 50, 100, 200, 400, 500, 550, 600, 650, 700, 800, 900, 1000, 1285]
hist, _ = np.histogram(dom_ts_arr, bins=bins)
for i, cnt in enumerate(hist):
    print(f"  tx_start [{bins[i]:5d},{bins[i+1]:5d}): {cnt:6,d}")

# ---- Look at the BAM records for a dominanted fragment ----
print(f"\n{'='*70}")
print("BAM INSPECTION: trace a dominated fragment through the BAM")
print(f"{'='*70}")

# Pick first few dominated cases and show ALL their BAM alignments
dom_qnames = set()
for rec in dominated_units[:20]:
    qn = frag_id_to_qname.get(rec["fid"])
    if qn:
        dom_qnames.add(qn)

# Collect all BAM records for these qnames
bam_records = defaultdict(list)
with pysam.AlignmentFile(MM2_BAM, "rb") as bam:
    for r in bam:
        if r.query_name in dom_qnames:
            nm = r.get_tag("NM") if r.has_tag("NM") else -1
            bam_records[r.query_name].append({
                "chrom": r.reference_name,
                "pos": r.reference_start,
                "cigar": r.cigarstring,
                "flag": r.flag,
                "nm": nm,
                "is_r1": bool(r.flag & 0x40),
                "is_r2": bool(r.flag & 0x80),
                "is_sec": bool(r.flag & 0x100),
                "is_rev": bool(r.flag & 0x10),
                "mate_pos": r.next_reference_start,
                "mate_chrom": r.next_reference_name,
            })

shown = 0
for rec in dominated_units[:8]:
    qn = frag_id_to_qname.get(rec["fid"])
    if not qn or qn not in bam_records:
        continue
    print(f"\n  qname={qn}")
    print(f"  EM gap={rec['gap']:.4f}, MANE ll={rec['mane_ll']:.4f}, "
          f"other_best={rec['other_best']:.4f}")
    print(f"  MANE tx=[{rec['mane_ts']},{rec['mane_te']}]")
    for al in bam_records[qn]:
        tag = "R1" if al["is_r1"] else "R2"
        sec = "SEC" if al["is_sec"] else "PRI"
        rev = "-" if al["is_rev"] else "+"
        print(f"    {tag} {sec} {al['chrom']}:{al['pos']:>10d} {rev}  "
              f"NM={al['nm']:2d}  cigar={al['cigar']}")
    shown += 1
    if shown >= 5:
        break

# ---- Cross-check: are dominated cases from cross-paired hits? ----
print(f"\n{'='*70}")
print("Fragment source analysis (from qname)")
print(f"{'='*70}")
src_cnts = Counter()
for rec in dominated_units:
    qn = frag_id_to_qname.get(rec["fid"], "?")
    src = qn.split(":")[0] if ":" in qn else qn
    src_cnts[src] += 1
print("Source transcripts for dominated units:")
for src, cnt in src_cnts.most_common(10):
    print(f"  {src:30s}: {cnt}")

# ---- PAIRED END OVERLAP DOUBLE-COUNTING ANALYSIS ----
# The TODO mentions: PE reads overlapping when flen < 2*readlen
# This could cause double-counted NM
print(f"\n{'='*70}")
print("DOUBLE-COUNTING CHECK: PE overlap NM double-counting?")
print(f"{'='*70}")

# For dominated fragments, check: do R1 and R2 overlap on chr12?
# If flen < 2*readlen, they overlap. Mismatches in the overlap
# would be counted in BOTH R1_NM and R2_NM → NM inflated for MANE

print("Checking PE overlap in dominated fragments:")
n_with_overlap = 0
n_no_overlap = 0
for rec in dominated_units[:200]:
    qn = frag_id_to_qname.get(rec["fid"])
    if not qn or qn not in bam_records:
        continue
    reads = bam_records[qn]
    pri_r1 = [r for r in reads if r["is_r1"] and not r["is_sec"] and r["chrom"] == "chr12"]
    pri_r2 = [r for r in reads if r["is_r2"] and not r["is_sec"] and r["chrom"] == "chr12"]
    if not pri_r1 or not pri_r2:
        continue
    # Check if cigar indicates they overlap
    # Read length ~ 150; flen = rec mane_te - mane_ts (in tx coords)
    tx_flen = rec["mane_te"] - rec["mane_ts"]
    # PE reads of ~150bp overlap if flen < 300bp
    # But tx_flen is in transcript space ~ similar to genomic flen for spliced
    r1 = pri_r1[0]
    r2 = pri_r2[0]
    # Estimate genomic flen from positions
    if r1["pos"] < r2["pos"]:
        gflen_approx = r2["pos"] - r1["pos"] + 150
    else:
        gflen_approx = r1["pos"] - r2["pos"] + 150
    r1_nm = r1["nm"]
    r2_nm = r2["nm"]
    if gflen_approx < 300:
        n_with_overlap += 1

frac_overlap = n_with_overlap / min(200, len(dominated_units)) if dominated_units else 0
print(f"  Dominated cases with PE overlap (flen<300): {n_with_overlap} / {min(200,len(dominated_units))} = {frac_overlap:.1%}")
