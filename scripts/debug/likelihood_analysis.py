#!/usr/bin/env python3
"""Deep per-candidate likelihood analysis for GAPDH multimappers.

Key questions:
1. Are MANE and shorter isoforms truly identical in log-likelihood for
   fragments in the shared 3' region?
2. What fraction of fragments CAN be distinguished (MANE has unique best LL)?
3. What scoring component drives any differences?
4. Why does ENST00000396859.5 absorb 6x more fragments than it should?
"""
import sys
sys.path.insert(0, "src")

import math
import gc
import numpy as np
from collections import defaultdict, Counter

from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer, _setup_geometry_and_estimator, _score_fragments
from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig
from rigel.stats import PipelineStats

BASE = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
IDX  = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
MM2_BAM = f"{BASE}/sim_minimap2.bam"

# ---- Load index ----
print("Loading index...")
tx_index = TranscriptIndex.load(IDX)
tx_df = tx_index.t_df

t_info = {}
for _, r in tx_df.iterrows():
    t_info[int(r["t_index"])] = {
        "t_id": r["t_id"],
        "g_name": r.get("g_name", "?"),
        "length": int(r.get("length", 0)),
        "is_basic": bool(r.get("is_basic", False)),
        "is_mane": bool(r.get("is_mane", False)),
        "is_nrna": bool(r.get("is_synthetic_nrna", False)),
        "start": int(r.get("start", 0)),
        "end": int(r.get("end", 0)),
    }

mane_tidx = next(ti for ti, inf in t_info.items()
                 if inf["is_mane"] and inf["g_name"] == "GAPDH")
gapdh_mrna_set = {ti for ti, inf in t_info.items()
                  if inf["g_name"] == "GAPDH" and not inf["is_nrna"]}
gapdh_basic_set = {ti for ti in gapdh_mrna_set if t_info[ti]["is_basic"]}

print(f"\nMANE t_index: {mane_tidx}")
print(f"\nGAPDH isoform genomic start positions:")
for ti in sorted(gapdh_mrna_set):
    inf = t_info[ti]
    print(f"  [{ti:6d}] {inf['t_id']:25s} start={inf['start']:>10d} "
          f"len={inf['length']:5d}  "
          f"{'MANE' if inf['is_mane'] else ('basic' if inf['is_basic'] else '')}")

# ---- Run pipeline scan ----
print("\nScanning MM2 BAM...")
scan_config = BamScanConfig(include_multimap=True)
pipeline_stats, strand_models, fl_models, buffer, _, _ = scan_and_buffer(
    MM2_BAM, tx_index, scan=scan_config
)
print(f"  Buffer size: {buffer.total_fragments:,} (includes secondary hits)")

em_config = EMConfig()
scoring_config = FragmentScoringConfig()
geometry, estimator = _setup_geometry_and_estimator(tx_index, fl_models, em_config)

print("Scoring fragments...")
stats_obj = PipelineStats()
sf = _score_fragments(
    buffer, tx_index, strand_models, fl_models, stats_obj, estimator,
    scoring_config, log_every=200000, annotations=None
)
print(f"  EM units: {sf.n_units:,}, Candidates: {sf.n_candidates:,}")

offsets   = sf.offsets
t_indices = sf.t_indices
log_liks  = sf.log_liks
tx_starts = sf.tx_starts
tx_ends   = sf.tx_ends
frag_ids  = sf.frag_ids

# ---- Per-unit LL analysis ----
print("\nAnalyzing per-unit log-likelihoods...")

TOL = 1e-4

case_counts = Counter()
n_gapdh_units_multi = 0  # MANE present + multiple GAPDH candidates

# Store for distributions
mane_ts_by_case = defaultdict(list)
gap_by_nCands = defaultdict(list)

samples = defaultdict(list)
MAX_SAMPLE = 6

mane_ll_vs_best_list = []

for unit_i in range(sf.n_units):
    s = int(offsets[unit_i])
    e = int(offsets[unit_i + 1])
    if e <= s:
        continue

    unit_t  = [int(t) for t in t_indices[s:e]]
    unit_ll = [float(ll) for ll in log_liks[s:e]]
    unit_ts = [int(ts) for ts in tx_starts[s:e]]
    unit_te = [int(te) for te in tx_ends[s:e]]

    # GAPDH mRNA candidates
    gapdh_idx = [j for j, t in enumerate(unit_t) if t in gapdh_mrna_set]
    if not gapdh_idx:
        continue
    if mane_tidx not in [unit_t[j] for j in gapdh_idx]:
        continue

    n_gapdh = len(gapdh_idx)
    if n_gapdh == 1:
        case_counts["mane_only"] += 1
        continue

    n_gapdh_units_multi += 1
    g_t   = [unit_t[j]  for j in gapdh_idx]
    g_ll  = [unit_ll[j] for j in gapdh_idx]
    g_ts  = [unit_ts[j] for j in gapdh_idx]
    g_te  = [unit_te[j] for j in gapdh_idx]

    mane_pos = g_t.index(mane_tidx)
    mane_ll  = g_ll[mane_pos]
    mane_ts_ = g_ts[mane_pos]
    mane_te_ = g_te[mane_pos]

    best_other_ll = max(ll for j, ll in enumerate(g_ll) if j != mane_pos)
    gap = best_other_ll - mane_ll

    mane_ll_vs_best_list.append((mane_ll, best_other_ll, gap, n_gapdh, mane_ts_))

    if abs(gap) <= TOL:
        case = "tied"
    elif gap < 0:
        case = "mane_best"
    else:
        case = "mane_dominated"

    case_counts[case] += 1
    mane_ts_by_case[case].append(mane_ts_)
    gap_by_nCands[n_gapdh].append(gap)

    if len(samples[case]) < MAX_SAMPLE:
        cand_info = []
        for j, (ti, ll, ts, te) in enumerate(zip(g_t, g_ll, g_ts, g_te)):
            flen = te - ts
            cand_info.append({
                "ti": ti, "ll": ll, "ts": ts, "te": te, "flen": flen,
                "t_id": t_info.get(ti, {}).get("t_id","?"),
                "g_start": t_info.get(ti, {}).get("start", 0),
                "t_len": t_info.get(ti, {}).get("length", 0),
                "is_mane": (ti == mane_tidx),
            })
        samples[case].append({
            "unit_i": unit_i,
            "cands": sorted(cand_info, key=lambda x: -x["ll"]),
            "gap": gap,
        })

total_with_mane = n_gapdh_units_multi + case_counts["mane_only"]
print(f"\n{'='*65}")
print(f"Units with MANE as candidate: {total_with_mane:,}")
print(f"  mane_only (deterministic):  {case_counts['mane_only']:,} "
      f"({100*case_counts['mane_only']/total_with_mane:.1f}%)")
print(f"  mane_best (MANE uniquely best): {case_counts['mane_best']:,} "
      f"({100*case_counts['mane_best']/total_with_mane:.1f}%)")
print(f"  tied (MANE = best other):   {case_counts['tied']:,} "
      f"({100*case_counts['tied']/total_with_mane:.1f}%)")
print(f"  mane_dominated:             {case_counts['mane_dominated']:,} "
      f"({100*case_counts['mane_dominated']/total_with_mane:.1f}%)")

arr = np.array(mane_ll_vs_best_list)
gap_arr = arr[:, 2]
print(f"\nLL gap stats (best_other - MANE_LL):")
print(f"  Mean={gap_arr.mean():.4f} Std={gap_arr.std():.4f} "
      f"Min={gap_arr.min():.4f} Max={gap_arr.max():.4f}")

# Histogram
edges = [-20, -10, -5, -2, -1, -0.5, -0.1, -TOL, TOL, 0.1, 0.5, 1, 2, 5, 10, 20]
hist, _ = np.histogram(gap_arr, bins=edges)
print("\nGap histogram (best_other - MANE_LL):")
for i, cnt in enumerate(hist):
    lo, hi = edges[i], edges[i+1]
    bar = "#" * min(int(cnt/100), 50)
    print(f"  [{lo:6.2f},{hi:5.2f}): {cnt:7,d} {bar}")

# ---- TX_START distributions ----
mane_len = t_info[mane_tidx]["length"]
print(f"\n{'='*65}")
print(f"MANE tx_start distributions by case  (MANE length={mane_len}nt)")
print("(tx_start=0 is 5' end of MANE; unique region = ~first 100-200nt)")
print(f"{'='*65}")
for case in ["mane_best", "tied", "mane_dominated"]:
    tsvs = mane_ts_by_case.get(case, [])
    if not tsvs:
        continue
    ts_arr = np.array(tsvs)
    print(f"\n  {case} (n={len(ts_arr):,})")
    print(f"    mean={ts_arr.mean():.0f} p25={np.percentile(ts_arr,25):.0f} "
          f"p50={np.median(ts_arr):.0f} p75={np.percentile(ts_arr,75):.0f} "
          f"max={ts_arr.max():.0f}")
    for thresh in [50, 100, 200, 400, 700]:
        frac = (ts_arr < thresh).mean()
        print(f"    tx_start < {thresh:4d}: {frac:.1%} ({(ts_arr<thresh).sum():,})")

# ---- Sample details ----
print(f"\n{'='*65}")
print("SAMPLE: TIED units — IDENTICAL LL for all GAPDH isoforms")
print("These are the isoform identifiability failures")
print(f"{'='*65}")
for s in samples["tied"][:3]:
    print(f"\n  unit={s['unit_i']} gap={s['gap']:.4f}")
    for c in s["cands"]:
        marker = " <<< MANE" if c["is_mane"] else ""
        print(f"    [{c['ti']:6d}] {c['t_id']:25s} ll={c['ll']:9.4f} "
              f"tx=[{c['ts']},{c['te']}] flen={c['flen']:4d} "
              f"g_start={c['g_start']:>10d} tlen={c['t_len']}{marker}")

print(f"\n{'='*65}")
print("SAMPLE: MANE_DOMINATED — another isoform has strictly higher LL")
print(f"{'='*65}")
for s in samples["mane_dominated"][:3]:
    print(f"\n  unit={s['unit_i']} gap={s['gap']:.4f}")
    for c in s["cands"]:
        marker = " <<< MANE" if c["is_mane"] else ""
        print(f"    [{c['ti']:6d}] {c['t_id']:25s} ll={c['ll']:9.4f} "
              f"tx=[{c['ts']},{c['te']}] flen={c['flen']:4d} "
              f"tlen={c['t_len']}{marker}")

print(f"\n{'='*65}")
print("SAMPLE: MANE_BEST — MANE uniquely has the highest LL")
print(f"{'='*65}")
for s in samples["mane_best"][:3]:
    print(f"\n  unit={s['unit_i']} gap={s['gap']:.4f}")
    for c in s["cands"]:
        marker = " <<< MANE" if c["is_mane"] else ""
        print(f"    [{c['ti']:6d}] {c['t_id']:25s} ll={c['ll']:9.4f} "
              f"tx=[{c['ts']},{c['te']}] flen={c['flen']:4d} "
              f"tlen={c['t_len']}{marker}")

# ---- LL gap vs n_candidates ----
print(f"\n{'='*65}")
print("LL GAP BY NUMBER OF GAPDH CANDIDATES")
print(f"{'='*65}")
for nc in sorted(gap_by_nCands.keys()):
    gaps = np.array(gap_by_nCands[nc])
    pct_tied = (np.abs(gaps) <= TOL).mean()
    pct_dom = (gaps > TOL).mean()
    print(f"  n_cands={nc:3d}: {len(gaps):6,d} units  "
          f"tied={pct_tied:.1%}  dominated={pct_dom:.1%}  "
          f"mean_gap={gaps.mean():.3f}")

# ---- What fraction of MANE fragments have ALL isoforms tied? ----
print(f"\n{'='*65}")
print("FUNDAMENTAL QUESTION: What fraction of GAPDH fragments")
print("cannot be assigned to the correct isoform?")
print(f"{'='*65}")
total_non_trivial = n_gapdh_units_multi
n_tied = case_counts["tied"]
n_mane_best = case_counts["mane_best"]
n_dom = case_counts["mane_dominated"]
print(f"\n  Total units with multiple GAPDH candidates: {total_non_trivial:,}")
print(f"  Identifiable (MANE uniquely best):           {n_mane_best:,} ({100*n_mane_best/total_non_trivial:.1f}%)")
print(f"  Not identifiable (tied or dominated):        {n_tied+n_dom:,} ({100*(n_tied+n_dom)/total_non_trivial:.1f}%)")
print(f"    Of which tied (can't distinguish):         {n_tied:,} ({100*n_tied/total_non_trivial:.1f}%)")
print(f"    Of which dominated (MANE loses LL):        {n_dom:,} ({100*n_dom/total_non_trivial:.1f}%)")
print()
print("CONCLUSION:")
print("  If most units are TIED → the EM is effectively guessing, distributing")
print("  fragments across all co-eligible isoforms equally (proportional to prior).")
print("  This explains why shorter isoforms absorb MANE fragments.")
print()
print("  If most units are MANE_BEST → there is a scoring BUG causing MANE to lose.")
