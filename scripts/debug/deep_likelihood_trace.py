#!/usr/bin/env python3
"""Deep per-candidate likelihood analysis for GAPDH multimappers.

Key questions:
1. Are MANE and shorter isoforms truly identical in log-likelihood for
   fragments in the shared 3' region?
2. What fraction of fragments CAN be distinguished (MANE has unique best LL)?
3. What scoring component (strand, FL, overhang, NM) drives any differences?
4. Why does ENST00000396859.5 capture 6x more fragments than it should?

(The original version below; this entire file is replaced by the rewrite)
-----------------------------------------------------------------------
ORIGINAL DOCSTRING:
Deep fragment likelihood tracer for GAPDH multimapper analysis.

Questions to answer:
1. For fragments from MANE's unique 5' region: does overhang correctly
   penalize shorter isoforms?
2. For fragments from the shared 3' region: are all isoforms truly tied
   in log-likelihood (identical scores)?
3. What is the NM landscape for primary vs secondary hits?
4. What does the EM see: can it distinguish MANE from ENST00000396859.5?

Strategy:
- Run rigel scan pipeline on mm2 BAM to get ScoredFragments
- Inspect per-unit log-likelihoods for GAPDH transcripts
- Use qname ordering (BAM is name-sorted) to map frag_ids back to source
"""
import sys
sys.path.insert(0, "src")

import math
import numpy as np
import pandas as pd
import pysam
from collections import defaultdict, Counter

# Rigel imports
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer, quant_from_buffer
from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig, PipelineConfig
from rigel.stats import PipelineStats

BASE = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
IDX = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
MM2_BAM = f"{BASE}/sim_minimap2.bam"
ORACLE_BAM = f"{BASE}/sim_oracle_nsorted.bam"

# ---- 1. Load index ----
print("Loading index...")
tx_index = TranscriptIndex.load(IDX)
tx_df = tx_index.t_df

# t_index → metadata
t_info = {}
for _, r in tx_df.iterrows():
    t_info[r["t_index"]] = {
        "t_id": r["t_id"],
        "g_name": r.get("g_name", "?"),
        "length": r.get("length", 0),
        "is_basic": bool(r.get("is_basic", False)),
        "is_mane": bool(r.get("is_mane", False)),
        "is_nrna": bool(r.get("is_synthetic_nrna", False)),
        "ref": r.get("ref", "?"),
        "start": int(r.get("start", 0)),
        "end": int(r.get("end", 0)),
    }

# GAPDH mRNA transcript indices
gapdh_tidx = [ti for ti, inf in t_info.items()
              if inf["g_name"] == "GAPDH" and not inf["is_nrna"]]
gapdh_tidx.sort()
mane_tidx = [ti for ti, inf in t_info.items() if inf["is_mane"]][0]
target_tidx = {47812, 47813, 47815, 47818, 47820, 47821}  # basic isoforms

print(f"\nGAPDH mRNA transcripts ({len(gapdh_tidx)} total):")
for ti in gapdh_tidx:
    inf = t_info[ti]
    star = "* MANE" if inf["is_mane"] else ("  basic" if inf["is_basic"] else "       ")
    print(f"  [{ti}] {star} {inf['t_id']:25s} len={inf['length']:5d} "
          f"start={inf['start']:>10d} end={inf['end']:>10d}")

# ---- 2. Build qname → truth source from mm2 BAM ----
# BAM is name-sorted; qnames encode the source transcript
# Format: {transcript_id}:{tx_start}-{tx_end}:{strand}:{frag_num}
print("\nParsing mm2 BAM qnames...")
qnames_ordered = []  # in order of first encounter (= buffer order for primary reads)
seen = set()
with pysam.AlignmentFile(MM2_BAM, "rb") as bam:
    for r in bam:
        qn = r.query_name
        if qn not in seen:
            seen.add(qn)
            qnames_ordered.append(qn)
print(f"  {len(qnames_ordered)} unique qnames")

def parse_source_tid(qname):
    """Extract source transcript-id from sim qname."""
    return qname.split(":")[0] if ":" in qname else qname

# Map qname → source transcript_id
qname_to_src = {qn: parse_source_tid(qn) for qn in qnames_ordered}
# frag_id in buffer = sequential index (0-based) in qnames_ordered
frag_id_to_qname = {i: qn for i, qn in enumerate(qnames_ordered)}
frag_id_to_src = {i: qname_to_src[qn] for i, qn in enumerate(qnames_ordered)}

# Reverse map: t_id → t_index
tid_to_tidx = {inf["t_id"]: ti for ti, inf in t_info.items()}

src_count = Counter(frag_id_to_src.values())
print("Top source transcripts (from qnames):")
for src, cnt in src_count.most_common(10):
    print(f"  {src}: {cnt}")

# ---- 3. Run full pipeline scan to get ScoredFragments ----
print("\nRunning scan_and_buffer on mm2 BAM...")
scan_config = BamScanConfig(include_multimap=True)
pipeline_stats, strand_models, fl_models, buffer, _, _ = scan_and_buffer(
    MM2_BAM, tx_index, scan=scan_config
)
print(f"  Buffer size: {buffer.total_fragments:,} fragments")

# Run scoring to get ScoredFragments
from rigel.scoring import FragmentScorer
from rigel.estimator import AbundanceEstimator
from rigel.scored_fragments import ScoredFragments
from rigel.scan import FragmentRouter

print("\nScoring fragments...")
em_config = EMConfig()
scoring_config = FragmentScoringConfig()
estimator = AbundanceEstimator(tx_index, em_config)
scorer = FragmentScorer.from_models(
    strand_models, fl_models, tx_index, scoring=scoring_config
)
router = FragmentRouter(scorer, estimator, PipelineStats(), tx_index, strand_models)
sf = router.scan(buffer, log_every=100000)
print(f"  EM units: {sf.n_units:,}")
print(f"  Candidates: {sf.n_candidates:,}")

# ---- 4. Analyze per-unit log-likelihoods ----
print("\nAnalyzing per-unit log-likelihoods...")

# For each EM unit, extract GAPDH candidates
# ScoredFragments layout: offsets[i]..offsets[i+1] = candidates for unit i
offsets = sf.offsets
t_indices = sf.t_indices
log_liks = sf.log_liks
frag_ids = sf.frag_ids

# Accumulate statistics
n_units = sf.n_units

# Per-unit stats for GAPDH
mane_has_best = 0
mane_tied_best = 0
mane_not_best = 0
mane_not_candidate = 0
mane_only_candidate = 0

# Detailed per-source analysis
src_stats = defaultdict(lambda: {
    "n_units": 0,
    "mane_best": 0,
    "mane_tied": 0,  # MANE tied with other GAPDH transcripts
    "mane_dominated": 0,
    "mane_not_candidate": 0,
    "mane_only": 0,
    "ll_gaps": [],  # best - mane_ll for cases where mane IS candidate
})

# Sample detailed units for inspection
sample_units = {"mane_dominated": [], "mane_tied": [], "mane_not_candidate": []}
MAX_SAMPLE = 5

total_gapdh_units = 0

for unit_i in range(n_units):
    start = offsets[unit_i]
    end = offsets[unit_i + 1]
    fid = int(frag_ids[unit_i])
    src_tid = frag_id_to_src.get(fid, "unknown")

    unit_t = t_indices[start:end]
    unit_ll = log_liks[start:end]

    # Find GAPDH mRNA candidates (exclude nRNA)
    gapdh_mask = np.array([
        (int(ti) in target_tidx) and not t_info.get(int(ti), {}).get("is_nrna", True)
        for ti in unit_t
    ])
    if not gapdh_mask.any():
        continue

    total_gapdh_units += 1
    gapdh_t = unit_t[gapdh_mask]
    gapdh_ll = unit_ll[gapdh_mask]
    best_ll_all = unit_ll.max()
    best_ll_gapdh = gapdh_ll.max()

    mane_in_candidates = int(mane_tidx) in [int(t) for t in gapdh_t]

    sst = src_stats[src_tid]
    sst["n_units"] += 1

    if not mane_in_candidates:
        mane_not_candidate += 1
        sst["mane_not_candidate"] += 1
        if len(sample_units["mane_not_candidate"]) < MAX_SAMPLE:
            sample_units["mane_not_candidate"].append({
                "unit_i": unit_i, "src": src_tid, "fid": fid,
                "candidates": [(int(t), float(ll)) for t, ll in zip(gapdh_t, gapdh_ll)],
            })
        continue

    mane_idx_in_arr = [j for j, t in enumerate(gapdh_t) if int(t) == int(mane_tidx)][0]
    mane_ll = float(gapdh_ll[mane_idx_in_arr])

    if len(gapdh_t) == 1:
        mane_only_candidate += 1
        sst["mane_only"] += 1
        continue

    # Compute gap between best and MANE
    other_lls = np.array([float(ll) for j, (t, ll) in enumerate(zip(gapdh_t, gapdh_ll))
                          if int(t) != int(mane_tidx)])
    best_other_ll = other_lls.max()
    gap = best_other_ll - mane_ll  # positive => another isoform beats MANE
    sst["ll_gaps"].append(gap)

    TOL = 1e-6
    if gap > TOL:
        mane_not_best += 1
        sst["mane_dominated"] += 1
        if len(sample_units["mane_dominated"]) < MAX_SAMPLE:
            sample_units["mane_dominated"].append({
                "unit_i": unit_i, "src": src_tid, "fid": fid,
                "mane_ll": mane_ll, "best_other_ll": best_other_ll, "gap": gap,
                "candidates": [(int(t), float(ll), t_info.get(int(t), {}).get("t_id","?"))
                               for t, ll in zip(gapdh_t, gapdh_ll)],
            })
    elif gap < -TOL:
        mane_has_best += 1
        sst["mane_best"] += 1
    else:
        mane_tied_best += 1
        sst["mane_tied"] += 1
        if len(sample_units["mane_tied"]) < MAX_SAMPLE:
            sample_units["mane_tied"].append({
                "unit_i": unit_i, "src": src_tid, "fid": fid,
                "mane_ll": mane_ll, "best_other_ll": best_other_ll,
                "candidates": [(int(t), float(ll), t_info.get(int(t), {}).get("t_id","?"))
                               for t, ll in zip(gapdh_t, gapdh_ll)],
            })

print(f"\n=== EM UNIT ANALYSIS ({total_gapdh_units:,} units with GAPDH candidates) ===")
print(f"  MANE only candidate (deterministic): {mane_only_candidate:,}")
print(f"  MANE has best LL (strictly):         {mane_has_best:,}")
print(f"  MANE tied for best LL:               {mane_tied_best:,}")
print(f"  MANE dominated by another isoform:   {mane_not_best:,}")
print(f"  MANE not even a candidate:           {mane_not_candidate:,}")
print()

# Distribution of ll_gaps for source=MANE fragments
mane_src_stats = src_stats.get("ENST00000229239.10", {})
if mane_src_stats and mane_src_stats.get("ll_gaps"):
    gaps = np.array(mane_src_stats["ll_gaps"])
    print(f"=== For src=MANE fragments ({mane_src_stats['n_units']:,} units with GAPDH cands) ===")
    print(f"  MANE only:           {mane_src_stats['mane_only']:,}")
    print(f"  MANE best:           {mane_src_stats['mane_best']:,}")
    print(f"  MANE tied:           {mane_src_stats['mane_tied']:,}")
    print(f"  MANE dominated:      {mane_src_stats['mane_dominated']:,}")
    print(f"  MANE not candidate:  {mane_src_stats['mane_not_candidate']:,}")
    print(f"  LL gap stats (best_other - MANE_LL):")
    print(f"    mean={gaps.mean():.4f} std={gaps.std():.4f} "
          f"min={gaps.min():.4f} max={gaps.max():.4f}")
    if (gaps > 0).any():
        print(f"    cases where another isoform has HIGHER LL: {(gaps>1e-6).sum()}")
    if (np.abs(gaps) < 1e-6).any():
        print(f"    cases where all GAPDH isoforms TIED: {(np.abs(gaps)<1e-6).sum()}")

print()
print("=== Per-source breakdown ===")
for src, st in sorted(src_stats.items(), key=lambda x: -x[1]["n_units"]):
    n = st["n_units"]
    if n < 10:
        continue
    gaps = np.array(st.get("ll_gaps", []))
    frac_tied = st["mane_tied"] / n if n else 0
    frac_dom = st["mane_dominated"] / n if n else 0
    frac_only = st["mane_only"] / n if n else 0
    frac_best = st["mane_best"] / n if n else 0
    frac_not = st["mane_not_candidate"] / n if n else 0
    print(f"  src={src:30s} n={n:6d} "
          f"only={frac_only:.2f} best={frac_best:.2f} "
          f"tied={frac_tied:.2f} dom={frac_dom:.2f} nc={frac_not:.2f}")

# ---- 5. Show sample units ----
print()
print("=== SAMPLE: MANE DOMINATED units (another isoform scores higher) ===")
for s in sample_units["mane_dominated"][:3]:
    print(f"\n  fid={s['fid']} src={s['src']} gap={s['gap']:.4f}")
    print(f"  MANE ll={s['mane_ll']:.4f} best_other={s['best_other_ll']:.4f}")
    for ti, ll, tid in sorted(s["candidates"], key=lambda x: -x[1]):
        marker = "<<MANE" if ti == int(mane_tidx) else ""
        print(f"    [{ti}] {tid:25s} ll={ll:8.4f} {marker}")

print()
print("=== SAMPLE: MANE TIED units ===")
for s in sample_units["mane_tied"][:2]:
    print(f"\n  fid={s['fid']} src={s['src']}")
    for ti, ll, tid in sorted(s["candidates"], key=lambda x: -x[1]):
        marker = "<<MANE" if ti == int(mane_tidx) else ""
        inf = t_info.get(ti, {})
        print(f"    [{ti}] {tid:25s} ll={ll:8.4f} len={inf.get('length','?'):5} {marker}")

# ---- 6. Analyze log-likelihood score components ----
# To understand if FL or strand or overhang is the discriminator,
# let me look at the MANE transcript vs closest competitor
print()
print("=== LL COMPONENT BREAKDOWN (for tied fragments) ===")
print("Need to check: are tied fragments in the shared 3' region?")
print("MANE-unique region: genomic 6534516-6534752 (before ENST00000396859.5 start)")
print("MANE+396859 shared: genomic 6534753-6538339")
print("MANE+619601 shared: genomic 6536604-6538374 (619601 starts much later!)")

# ---- 7. Analyze per-transcript fragment length distributions ----
# For MANE-source fragments that end up tied, where are they genomically?
# Use the frag_id → qname → parse tx coords to understand position
print()
print("=== Checking frag_id → position mapping for tied MANE fragments ===")
n_show = 10
shown = 0
for unit_i in range(n_units):
    if shown >= n_show:
        break
    start = offsets[unit_i]
    end = offsets[unit_i + 1]
    fid = int(frag_ids[unit_i])
    src_tid = frag_id_to_src.get(fid, "unknown")
    if src_tid != "ENST00000229239.10":
        continue

    unit_t = t_indices[start:end]
    unit_ll = log_liks[start:end]

    gapdh_mRNA_t = [(int(t), float(ll), t_info.get(int(t), {}).get("t_id","?"))
                    for t, ll in zip(unit_t, unit_ll)
                    if int(t) in target_tidx and not t_info.get(int(t), {}).get("is_nrna", True)]

    if len(gapdh_mRNA_t) < 2:
        continue

    # Check if tied
    lls = [ll for _, ll, _ in gapdh_mRNA_t]
    if max(lls) - min(lls) > 0.01:
        continue  # not tied

    qname = frag_id_to_qname.get(fid, "?")
    print(f"\n  fid={fid} qname={qname}")
    for ti, ll, tid in sorted(gapdh_mRNA_t, key=lambda x: -x[1]):
        marker = "<<MANE" if ti == int(mane_tidx) else ""
        inf = t_info.get(ti, {})
        print(f"    [{ti}] {tid:25s} ll={ll:8.4f} start={inf.get('start','?')} {marker}")
    shown += 1
