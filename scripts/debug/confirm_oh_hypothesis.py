#!/usr/bin/env python3
"""Confirm that the 4.6052 gap in dominated cases comes from oh=1, not NM.

Since nm is per-FRAGMENT (same for all candidates) and exon_bp is per-CANDIDATE,
the only source of LL variation between MANE and a competitor for the same
fragment is different exon_bp values, which drives oh = read_length - exon_bp.

This script directly inspects the raw buffer arrays for the first N dominated
units to show:
  - read_length[fid]: total aligned read length (fragment-level)
  - nm[fid]: NM count (fragment-level, same for all candidates)
  - exon_bp for MANE: leads to oh_mane = read_length - exon_bp_mane
  - exon_bp for best competitor: leads to oh_win = read_length - exon_bp_win
  - oh_mane should be oh_win + 1
"""
import sys
sys.path.insert(0, "src")

import math
import numpy as np
from collections import Counter

from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer, _setup_geometry_and_estimator, _score_fragments
from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig
from rigel.stats import PipelineStats

BASE = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
IDX  = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
MM2_BAM = f"{BASE}/sim_minimap2.bam"

OH_PEN = math.log(0.01)   # -4.6052
TOL = 1e-4

print("Loading index...")
tx_index = TranscriptIndex.load(IDX)
tx_df = tx_index.t_df

t_info = {}
for _, r in tx_df.iterrows():
    t_info[int(r["t_index"])] = {
        "t_id": r["t_id"],
        "g_name": r.get("g_name", "?"),
        "length": int(r.get("length", 0)),
        "is_mane": bool(r.get("is_mane", False)),
        "is_nrna": bool(r.get("is_synthetic_nrna", False)),
        "start": int(r.get("start", 0)),
        "end": int(r.get("end", 0)),
    }

mane_tidx = next(ti for ti, inf in t_info.items() if inf["is_mane"] and inf["g_name"] == "GAPDH")
gapdh_mrna_set = {ti for ti, inf in t_info.items() if inf["g_name"] == "GAPDH" and not inf["is_nrna"]}
print(f"  MANE t_index: {mane_tidx}")
print(f"  GAPDH mRNA transcripts: {sorted(gapdh_mrna_set)}")

print("\nScanning + buffering MM2 BAM...")
scan_config = BamScanConfig(include_multimap=True)
pipeline_stats, strand_models, fl_models, buffer, _, _ = scan_and_buffer(
    MM2_BAM, tx_index, scan=scan_config
)

# Build frag_id → raw fragment data cache BEFORE scoring releases the buffer
# _score_fragments calls buffer.release() at end, so we must cache here
print("Building frag_id → raw data cache from buffer (before scoring)...")
fid_raw = {}  # fid → (rl, nm, {t_index: exon_bp})
for chunk in buffer.iter_chunks():
    for row_idx in range(chunk.size):
        fid = int(chunk.frag_id[row_idx])
        rl  = int(chunk.read_length[row_idx])
        nm  = int(chunk.nm[row_idx])
        t_start = int(chunk.t_offsets[row_idx])
        t_end   = int(chunk.t_offsets[row_idx + 1])
        ebp_map = {
            int(chunk.t_indices[k]): int(chunk.exon_bp[k])
            for k in range(t_start, t_end)
        }
        fid_raw[fid] = (rl, nm, ebp_map)
print(f"  Cached {len(fid_raw):,} fragments")

print("\nScoring (buffer will be released after)...")
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

# ---- Find dominated units ----
print("\nFinding dominated units...")
dominated = []
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
    if not gapdh_idx or mane_tidx not in [unit_t[j] for j in gapdh_idx] or len(gapdh_idx) <= 1:
        continue

    g_t  = [unit_t[j]  for j in gapdh_idx]
    g_ll = [unit_ll[j] for j in gapdh_idx]
    g_ts = [unit_ts[j] for j in gapdh_idx]
    g_te = [unit_te[j] for j in gapdh_idx]

    mane_pos = g_t.index(mane_tidx)
    mane_ll  = g_ll[mane_pos]
    other_best_ll = max(ll for j, ll in enumerate(g_ll) if j != mane_pos)
    best_j = next(j for j, ll in enumerate(g_ll) if j != mane_pos and abs(ll - other_best_ll) < TOL)
    gap = other_best_ll - mane_ll

    if gap > TOL:  # MANE is dominated
        dominated.append({
            "unit_i": unit_i,
            "fid": int(frag_ids[unit_i]),
            "gap": gap,
            "mane_ll": mane_ll,
            "best_ll": other_best_ll,
            "best_tidx": g_t[best_j],
            "mane_ts": g_ts[mane_pos],
            "mane_te": g_te[mane_pos],
            "best_ts": g_ts[best_j],
            "best_te": g_te[best_j],
            "all_cands": list(zip(g_t, g_ll, g_ts, g_te)),
        })

print(f"  Found {len(dominated):,} dominated units")

# ---- For each dominated unit: extract raw buffer arrays ----
print(f"\n{'='*80}")
print("RAW BUFFER DATA for first 30 dominated units")
print("  rl  = read_length (fragment-level, same for all candidates)")
print("  nm  = NM count (fragment-level, same for all candidates)")
print("  ebp = exon_bp (candidate-level, varies by transcript)")
print("  oh  = rl - ebp")
print(f"{'='*80}")

oh_diffs = []
for rec in dominated[:200]:
    fid = rec["fid"]
    if fid not in fid_raw:
        continue
    rl, nm, ebp_map = fid_raw[fid]

    mane_ebp = ebp_map.get(mane_tidx)
    best_ebp = ebp_map.get(rec["best_tidx"])

    if mane_ebp is None or best_ebp is None:
        continue

    oh_mane = rl - mane_ebp
    oh_best = rl - best_ebp
    oh_diff = oh_mane - oh_best
    oh_diffs.append((oh_mane, oh_best, oh_diff, nm))

    if len(oh_diffs) <= 30:
        mane_t = t_info.get(mane_tidx, {})
        best_t = t_info.get(rec["best_tidx"], {})
        print(f"\n  unit={rec['unit_i']:5d} fid={fid:5d}  gap={rec['gap']:.4f}")
        print(f"    rl={rl}  nm={nm}  nm_pen={nm * math.log(0.1):.4f}")
        print(f"    MANE  ({mane_tidx:5d}) {mane_t.get('t_id','?'):25s}: "
              f"ebp={mane_ebp:3d}  oh={oh_mane}")
        print(f"    BEST  ({rec['best_tidx']:5d}) {best_t.get('t_id','?'):25s}: "
              f"ebp={best_ebp:3d}  oh={oh_best}")
        print(f"    oh_diff(MANE-BEST)={oh_diff}  → expected ll gap = {oh_diff * OH_PEN:.4f}")

# ---- Summary statistics ----
print(f"\n{'='*80}")
print("SUMMARY: oh_diff distribution across ALL dominated units")
print(f"{'='*80}")

oh_diff_counts = Counter(oh_diff for _, _, oh_diff, _ in oh_diffs)
for diff, cnt in sorted(oh_diff_counts.items()):
    print(f"  oh_diff = {diff:+d}  count={cnt:5,d}  "
          f"(pred gap={diff * abs(OH_PEN):.4f})")

nm_vals = [nm for _, _, _, nm in oh_diffs]
print(f"\n  NM distribution in dominated units:")
for nm_v, cnt in sorted(Counter(nm_vals).items()):
    print(f"    nm={nm_v}: {cnt:5,d}")

print(f"\n  Total dominated units analyzed: {len(oh_diffs):,}")
frac_oh1 = sum(1 for _, _, d, _ in oh_diffs if d == 1) / max(len(oh_diffs), 1)
print(f"  Fraction with oh_diff=+1 (MANE has 1 extra overhang): {frac_oh1:.1%}")
frac_oh2 = sum(1 for _, _, d, _ in oh_diffs if d == 2) / max(len(oh_diffs), 1)
print(f"  Fraction with oh_diff=+2 (MANE has 2 extra overhang): {frac_oh2:.1%}")

print("\nDONE")
