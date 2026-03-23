#!/usr/bin/env python3
"""Debug: why are only 24/1656 dominated units found via buffer t_indices?"""
import sys
sys.path.insert(0, "src")

import math
import numpy as np

from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer, _setup_geometry_and_estimator, _score_fragments
from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig
from rigel.stats import PipelineStats

BASE = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
IDX  = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
MM2_BAM = f"{BASE}/sim_minimap2.bam"

TOL = 1e-4

print("Loading index...")
tx_index = TranscriptIndex.load(IDX)
tx_df = tx_index.t_df

t_info = {}
for _, r in tx_df.iterrows():
    t_info[int(r["t_index"])] = {
        "t_id": r["t_id"],
        "g_name": r.get("g_name", "?"),
        "is_mane": bool(r.get("is_mane", False)),
        "is_nrna": bool(r.get("is_synthetic_nrna", False)),
    }

mane_tidx = next(ti for ti, inf in t_info.items() if inf["is_mane"] and inf["g_name"] == "GAPDH")
gapdh_mrna_set = {ti for ti, inf in t_info.items() if inf["g_name"] == "GAPDH" and not inf["is_nrna"]}
print(f"  MANE t_index: {mane_tidx}")
print(f"  GAPDH mRNA indices: {sorted(gapdh_mrna_set)}")

print("\nScanning + buffering...")
scan_config = BamScanConfig(include_multimap=True)
pipeline_stats, strand_models, fl_models, buffer, _, _ = scan_and_buffer(
    MM2_BAM, tx_index, scan=scan_config
)

# Before scoring: capture buffer data for a FEW specific fids
# We'll identify dominated fids AFTER scoring, so first capture EVERYTHING
# But we can't store everything — let's capture t_indices per fid
print("Capturing buffer fid→t_indices+ebp map before scoring...")
fid_tindices = {}  # fid → set of t_indices in buffer
fid_raw = {}       # fid → (rl, nm, {ti: ebp})
for chunk in buffer.iter_chunks():
    for row_idx in range(chunk.size):
        fid = int(chunk.frag_id[row_idx])
        rl  = int(chunk.read_length[row_idx])
        nm  = int(chunk.nm[row_idx])
        t_start = int(chunk.t_offsets[row_idx])
        t_end   = int(chunk.t_offsets[row_idx + 1])
        tis = {int(chunk.t_indices[k]): int(chunk.exon_bp[k]) for k in range(t_start, t_end)}
        fid_tindices[fid] = set(tis.keys())
        fid_raw[fid] = (rl, nm, tis)
print(f"  {len(fid_raw):,} fragments in buffer")

# Check presence of typical GAPDH indices
print("\nSample check: how many buffer fragments have MANE in t_indices?")
has_mane = sum(1 for s in fid_tindices.values() if mane_tidx in s)
print(f"  Fragments with MANE ({mane_tidx}) in buffer t_indices: {has_mane:,}")
has_any_gapdh = sum(1 for s in fid_tindices.values() if s & gapdh_mrna_set)
print(f"  Fragments with any GAPDH mRNA:                          {has_any_gapdh:,}")

# Show range of t_indices used in buffer
all_tis = set()
for s in fid_tindices.values():
    all_tis.update(s)
print(f"  Unique t_indices in buffer: {len(all_tis):,}")
gapdh_in_buf = all_tis & gapdh_mrna_set
print(f"  GAPDH mRNA t_indices present in buffer: {sorted(gapdh_in_buf)}")

print("\nScoring...")
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

# Find dominated units
dominated = []
for unit_i in range(sf.n_units):
    s = int(offsets[unit_i])
    e = int(offsets[unit_i + 1])
    if e <= s:
        continue
    unit_t  = [int(t) for t in t_indices[s:e]]
    unit_ll = [float(l) for l in log_liks[s:e]]
    gapdh_idx = [j for j, t in enumerate(unit_t) if t in gapdh_mrna_set]
    if not gapdh_idx or mane_tidx not in [unit_t[j] for j in gapdh_idx] or len(gapdh_idx) <= 1:
        continue
    g_ll = [unit_ll[j] for j in gapdh_idx]
    mane_pos_in_g = [unit_t[j] for j in gapdh_idx].index(mane_tidx)
    mane_ll = g_ll[mane_pos_in_g]
    other_best_ll = max(ll for j, ll in enumerate(g_ll) if j != mane_pos_in_g)
    gap = other_best_ll - mane_ll
    if gap > TOL:
        fid = int(frag_ids[unit_i])
        best_j = next(j for j, ll in enumerate(g_ll) if j != mane_pos_in_g and abs(ll - other_best_ll) < TOL)
        best_tidx = [unit_t[j] for j in gapdh_idx][best_j]
        dominated.append({"unit_i": unit_i, "fid": fid, "gap": gap, "best_tidx": best_tidx})

print(f"  Dominated units: {len(dominated):,}")

# Key diagnostic: for dominated units, check what's in fid_raw
print("\nDIAGNOSTIC: First 20 dominated units — buffer vs ScoredFragments t_indices")
for rec in dominated[:20]:
    fid = rec["fid"]
    buf_tis = fid_raw.get(fid, (None, None, {}))[2].keys() if fid in fid_raw else "NOT FOUND"
    mane_in_buf = mane_tidx in (fid_raw.get(fid, (None, None, {}))[2] if fid in fid_raw else {})
    best_in_buf = rec["best_tidx"] in (fid_raw.get(fid, (None, None, {}))[2] if fid in fid_raw else {})
    if fid in fid_raw:
        buf_gapdh = sorted(gapdh_mrna_set & set(fid_raw[fid][2].keys()))
    else:
        buf_gapdh = "NOT FOUND"
    print(f"  unit={rec['unit_i']:5d} fid={fid:6d}  gap={rec['gap']:.4f}")
    print(f"    MANE({mane_tidx}) in buf: {mane_in_buf}  BEST({rec['best_tidx']}) in buf: {best_in_buf}")
    print(f"    GAPDH mRNA in buf: {buf_gapdh}")

# Check if there is a t_index shift
# Maybe ScoredFragments t_indices are offset by n_nrna or some geometry factor
print("\nGeometry info about MANE t_index:")
print(f"  mane_tidx in ScoredFragments: {mane_tidx}")
print(f"  t_info[mane_tidx] = {t_info.get(mane_tidx)}")
# Find all indices near mane_tidx that are in buffer for dominated units
for rec in dominated[:5]:
    fid = rec["fid"]
    if fid not in fid_raw:
        continue
    buf_keys = sorted(fid_raw[fid][2].keys())
    close_to_mane = [k for k in buf_keys if abs(k - mane_tidx) <= 100]
    print(f"\n  fid={fid}: buffer t_indices near {mane_tidx}: {close_to_mane}")
    for k in close_to_mane:
        print(f"    t_index={k}: {t_info.get(k, '?')}")
