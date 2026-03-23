#!/usr/bin/env python3
"""Pinpoint the splice junction mismatch: WHY is MANE absent?

We know:
- mane_absent reads ARE spliced (100% have N in CIGAR)
- They have oh≥0 for whatever transcript they score best with
- MANE must be getting HIGH overhang (and pruned) due to an incorrect splice junction

Hypothesis: minimap2 uses a slightly different splice acceptor/donor position than
MANE's annotated junction. A few-base offset at the acceptor site makes MANE's
exon_bp lower (higher oh) than the competitor, and MANE gets pruned.

This script:
1. Gets the MANE transcript exon coordinates from the index
2. For a sample of mane_absent units, shows exactly which genomic junction
   was used and whether it matches MANE's annotation
3. Computes the oh values that MANE and the best competitor would get
"""
import sys
sys.path.insert(0, "src")

import math
import pysam
import numpy as np
from collections import defaultdict

from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer, _setup_geometry_and_estimator, _score_fragments
from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig
from rigel.stats import PipelineStats

BASE = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
IDX  = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
MM2_BAM    = f"{BASE}/sim_minimap2.bam"
ORACLE_BAM = f"{BASE}/sim_oracle_nsorted.bam"
TOL = 1e-4

print("Loading index...")
tx_index = TranscriptIndex.load(IDX)
tx_df = tx_index.t_df

t_info = {}
for _, r in tx_df.iterrows():
    t_info[int(r["t_index"])] = {
        "t_id": r["t_id"], "g_name": r.get("g_name", "?"),
        "length": int(r.get("length", 0)),
        "is_mane": bool(r.get("is_mane", False)),
        "is_nrna": bool(r.get("is_synthetic_nrna", False)),
        "chrom": str(r.get("chrom", "?")),
        "start": int(r.get("start", 0)),
        "end": int(r.get("end", 0)),
    }

mane_tidx = next(ti for ti, inf in t_info.items() if inf["is_mane"] and inf["g_name"] == "GAPDH")
gapdh_mrna_set = {ti for ti, inf in t_info.items() if inf["g_name"] == "GAPDH" and not inf["is_nrna"]}
mane_t_id = t_info[mane_tidx]["t_id"]
print(f"  MANE={mane_t_id} t_index={mane_tidx}")
print(f"  MANE genomic: chrom={t_info[mane_tidx]['chrom']} "
      f"start={t_info[mane_tidx]['start']} end={t_info[mane_tidx]['end']}")

# Get exon intervals from index
# tx_df should have exon start/end info, or we get it from the exon table
print("\nLooking for exon table in index...")
if hasattr(tx_index, 'exon_df'):
    exon_df = tx_index.exon_df
    mane_exons = exon_df[exon_df["t_index"] == mane_tidx].sort_values("start")
    print(f"  MANE exons from exon_df:")
    for _, e in mane_exons.iterrows():
        print(f"    exon: {e.get('chrom','?')}:{e['start']}-{e['end']}")
elif hasattr(tx_index, 'e_df'):
    exon_df = tx_index.e_df
    mane_exons = exon_df[exon_df["t_index"] == mane_tidx].sort_values("start")
    print(f"  MANE exons from e_df ({len(mane_exons)} exons):")
    for _, e in mane_exons.iterrows():
        print(f"    exon: {e.get('chrom','?')}:{int(e['start'])}-{int(e['end'])}")
else:
    print("  No exon table found, checking index attributes:")
    for attr in dir(tx_index):
        if not attr.startswith('_') and 'exon' in attr.lower():
            print(f"    {attr}")
    print("  Available attrs:", [a for a in dir(tx_index) if not a.startswith('_') and 'df' in a.lower()])
    mane_exons = None


def extract_source_tid(qn): return qn.split(":")[0]


def extract_junctions_from_cigar(chrom: str, pos: int, cigar: str) -> list:
    """Extract (donor_pos, acceptor_pos) for each N in cigar."""
    if not cigar:
        return []
    junctions = []
    cur_pos = pos
    import re
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    for length, op in ops:
        length = int(length)
        if op in ('M', 'D', '=', 'X'):
            cur_pos += length
        elif op == 'N':  # intron
            donor = cur_pos
            acceptor = cur_pos + length
            junctions.append((chrom, donor, acceptor, length))
            cur_pos += length
        elif op == 'S':  # soft clip doesn't advance
            pass
        elif op == 'I':  # insertion doesn't advance ref
            pass
        elif op == 'H':  # hard clip
            pass
    return junctions


# Score mm2 and find mane_absent units
print("\nScoring minimap2 BAM...")
# Build fid→qname
fid_to_qname = {}
seen = set()
fid = 0
with pysam.AlignmentFile(MM2_BAM, "rb") as bam:
    for r in bam:
        if r.query_name not in seen:
            seen.add(r.query_name)
            fid_to_qname[fid] = r.query_name
            fid += 1

scan_config = BamScanConfig(include_multimap=True)
_, strand_models, fl_models, buffer, _, _ = scan_and_buffer(MM2_BAM, tx_index, scan=scan_config)

# Capture raw buffer data before scoring
fid_raw = {}
for chunk in buffer.iter_chunks():
    for row_idx in range(chunk.size):
        fid = int(chunk.frag_id[row_idx])
        rl  = int(chunk.read_length[row_idx])
        nm  = int(chunk.nm[row_idx])
        t_start = int(chunk.t_offsets[row_idx])
        t_end   = int(chunk.t_offsets[row_idx + 1])
        ebp_map = {int(chunk.t_indices[k]): int(chunk.exon_bp[k]) for k in range(t_start, t_end)}
        fid_raw[fid] = (rl, nm, ebp_map)

em_config = EMConfig()
scoring_config = FragmentScoringConfig()
geometry, estimator = _setup_geometry_and_estimator(tx_index, fl_models, em_config)
stats_obj = PipelineStats()
sf = _score_fragments(buffer, tx_index, strand_models, fl_models, stats_obj, estimator,
                      scoring_config, log_every=500000, annotations=None)

# Find mane_absent units
mane_absent = []
for unit_i in range(sf.n_units):
    s = int(sf.offsets[unit_i]); e = int(sf.offsets[unit_i + 1])
    if e <= s: continue
    unit_t = [int(t) for t in sf.t_indices[s:e]]
    unit_ll = [float(l) for l in sf.log_liks[s:e]]
    gapdh_idx = [j for j, t in enumerate(unit_t) if t in gapdh_mrna_set]
    if not gapdh_idx: continue
    fid = int(sf.frag_ids[unit_i])
    g_t = [unit_t[j] for j in gapdh_idx]
    g_ll = [unit_ll[j] for j in gapdh_idx]
    if mane_tidx in g_t: continue
    best_j = max(range(len(g_t)), key=lambda j: g_ll[j])
    mane_absent.append({
        "unit_i": unit_i, "fid": fid,
        "qn": fid_to_qname.get(fid, "?"),
        "best_tidx": g_t[best_j], "best_ll": g_ll[best_j],
    })
print(f"  mane_absent: {len(mane_absent):,}")

# Build BAM record index for mane_absent qnames
key_qnames = {r["qn"] for r in mane_absent[:500]}
bam_recs = defaultdict(list)
with pysam.AlignmentFile(MM2_BAM, "rb") as bam:
    for r in bam:
        if r.query_name in key_qnames:
            bam_recs[r.query_name].append({
                "chrom": r.reference_name, "pos": r.reference_start,
                "pos_end": r.reference_end, "cigar": r.cigarstring,
                "is_r1": bool(r.flag & 0x40), "is_sec": bool(r.flag & 0x100),
                "nm": r.get_tag("NM") if r.has_tag("NM") else -1,
            })

# Also get oracle BAM records for same qnames to compare
oracle_recs = defaultdict(list)
with pysam.AlignmentFile(ORACLE_BAM, "rb") as bam:
    for r in bam:
        if r.query_name in key_qnames:
            oracle_recs[r.query_name].append({
                "chrom": r.reference_name, "pos": r.reference_start,
                "pos_end": r.reference_end, "cigar": r.cigarstring,
                "is_r1": bool(r.flag & 0x40), "is_sec": bool(r.flag & 0x100),
                "nm": r.get_tag("NM") if r.has_tag("NM") else -1,
            })

# Collect junction statistics for mane_absent units
print(f"\n{'='*80}")
print("JUNCTION ANALYSIS: what junctions does minimap2 use for mane_absent units?")
print(f"{'='*80}")

all_mm2_junctions = defaultdict(int)   # (chrom, donor, acceptor) → count
junction_vs_mane_absent = defaultdict(int)

for rec in mane_absent[:500]:
    qn = rec["qn"]
    mm2_aligns = [a for a in bam_recs.get(qn, []) if not a["is_sec"]]
    for a in mm2_aligns:
        junctions = extract_junctions_from_cigar(a["chrom"], a["pos"], a["cigar"])
        for j in junctions:
            all_mm2_junctions[j] += 1

print(f"\n  Top splice junctions used by mane_absent fragments in mm2 (first 500 mane_absent):")
print(f"  {'chrom':8s}  {'donor':>12s}  {'acceptor':>12s}  {'intron_len':>10s}  {'count':>8s}")
for (chrom, donor, acceptor, intron_len), cnt in sorted(all_mm2_junctions.items(), key=lambda x: -x[1])[:20]:
    print(f"  {chrom:8s}  {donor:>12,d}  {acceptor:>12,d}  {intron_len:>10,d}  {cnt:>8,d}")

# Compare to oracle junctions for same qnames
print(f"\n{'='*80}")
print("JUNCTION COMPARISON: oracle vs minimap2 for SAME queries")
print(f"{'='*80}")
print("\n  Fragment-level comparison (mane_absent qnames, primary alignments only):")
junc_matches = 0
junc_mismatches = 0
examples_shown = 0

for rec in mane_absent[:200]:
    qn = rec["qn"]
    src = extract_source_tid(qn)
    mm2_aligns = [a for a in bam_recs.get(qn, []) if not a["is_sec"]]
    oracle_aligns = [a for a in oracle_recs.get(qn, []) if not a["is_sec"]]

    mm2_junctions = []
    for a in mm2_aligns:
        mm2_junctions += extract_junctions_from_cigar(a["chrom"], a["pos"], a["cigar"])

    oracle_junctions = []
    for a in oracle_aligns:
        oracle_junctions += extract_junctions_from_cigar(a["chrom"], a["pos"], a["cigar"])

    # Only look at chr12 junctions
    mm2_j = [(d, ac, il) for (ch, d, ac, il) in mm2_junctions if ch == "chr12"]
    ora_j = [(d, ac, il) for (ch, d, ac, il) in oracle_junctions if ch == "chr12"]

    if mm2_j == ora_j:
        junc_matches += 1
    else:
        junc_mismatches += 1
        if examples_shown < 20:
            best_t = t_info.get(rec["best_tidx"], {})
            raw = fid_raw.get(rec["fid"], (0, 0, {}))
            ebp_map = raw[2]
            print(f"\n  qn={qn}  src={src}")
            print(f"  best_cand={best_t.get('t_id','?')}  mm2_ll={rec['best_ll']:.4f}")
            print(f"  buffer t_indices present: {sorted(ti for ti in ebp_map.keys() if ti in gapdh_mrna_set)}")
            print(f"  buffer ebp values: { {ti: ebp_map[ti] for ti in sorted(ebp_map.keys()) if ti in gapdh_mrna_set} }")
            print(f"  mm2 chr12 junctions:    {mm2_j}")
            print(f"  oracle chr12 junctions: {ora_j}")
            examples_shown += 1

pct_mismatch = 100 * junc_mismatches / max(junc_matches + junc_mismatches, 1)
print(f"\n  Junction comparison for first 200 mane_absent units:")
print(f"    same junctions (mm2 == oracle):      {junc_matches:5,d}  "
      f"({100*junc_matches/(junc_matches+junc_mismatches):.1f}%)")
print(f"    DIFFERENT junctions (mm2 != oracle): {junc_mismatches:5,d}  "
      f"({pct_mismatch:.1f}%)")

# Buffer: does MANE appear with high overhang?
print(f"\n{'='*80}")
print("RAW BUFFER ANALYSIS: for mane_absent, what overhang does MANE have in the buffer?")
print(f"{'='*80}")
print("(MANE in buffer = pruned due to high oh; MANE absent from buffer = not even found in cgranges)")

mane_in_buf = 0
mane_absent_from_buf = 0
oh_dist = defaultdict(int)

for rec in mane_absent:
    raw = fid_raw.get(rec["fid"], (0, 0, {}))
    rl, nm, ebp_map = raw
    if mane_tidx in ebp_map:
        mane_in_buf += 1
        mane_ebp = ebp_map[mane_tidx]
        oh = rl - mane_ebp
        oh_dist[oh] += 1
    else:
        mane_absent_from_buf += 1

total = mane_in_buf + mane_absent_from_buf
print(f"\n  MANE present in buffer (found by cgranges, then pruned): "
      f"{mane_in_buf:,d} ({100*mane_in_buf/max(total,1):.1f}%)")
print(f"  MANE absent from buffer (cgranges didn't find it):       "
      f"{mane_absent_from_buf:,d} ({100*mane_absent_from_buf/max(total,1):.1f}%)")

if oh_dist:
    print(f"\n  OH distribution for MANE (when present in buffer, before pruning):")
    for oh, cnt in sorted(oh_dist.items())[:20]:
        penalty = oh * math.log(0.01)
        print(f"    oh={oh:3d} (ll_pen={penalty:8.2f}):  {cnt:5,d}")

print("\nDONE")
