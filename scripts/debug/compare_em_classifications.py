#!/usr/bin/env python3
"""Compare oracle vs minimap2 EM unit classifications.

Key question: oracle achieves ~70% TPM for MANE while minimap2 gets only ~40%.
Both have ~35k MANE-source frags in the GAPDH EM, so the difference MUST come
from the likelihood structure — how many units are tied vs MANE-unique vs dominated.

This script runs the mane_only/mane_best/tied/mane_dominated analysis on
BOTH oracle and minimap2 and reports side-by-side, then samples the 'tied'
units in oracle to understand why MANE wins anyway.
"""
import sys
sys.path.insert(0, "src")

import math
import gc
import numpy as np
from collections import Counter

from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer, _setup_geometry_and_estimator, _score_fragments
from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig
from rigel.stats import PipelineStats

BASE = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
IDX  = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
MM2_BAM    = f"{BASE}/sim_minimap2.bam"
ORACLE_BAM = f"{BASE}/sim_oracle_nsorted.bam"

TOL = 1e-4
OH_PEN = math.log(0.01)

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
    }

mane_tidx = next(ti for ti, inf in t_info.items() if inf["is_mane"] and inf["g_name"] == "GAPDH")
gapdh_mrna_set = {ti for ti, inf in t_info.items() if inf["g_name"] == "GAPDH" and not inf["is_nrna"]}
print(f"  MANE t_index={mane_tidx}, GAPDH mRNA count={len(gapdh_mrna_set)}")


def classify_units(sf, label: str):
    """Classify all GAPDH EM units by mane_only/best/tied/dominated."""
    offsets   = sf.offsets
    t_indices = sf.t_indices
    log_liks  = sf.log_liks
    tx_starts = sf.tx_starts
    tx_ends   = sf.tx_ends
    frag_ids  = sf.frag_ids

    mane_only = []
    mane_best = []
    tied = []
    mane_dominated = []
    mane_absent = []  # mane not a candidate at all in a GAPDH unit

    for unit_i in range(sf.n_units):
        s = int(offsets[unit_i])
        e = int(offsets[unit_i + 1])
        if e <= s:
            continue
        unit_t  = [int(t) for t in t_indices[s:e]]
        unit_ll = [float(ll) for ll in log_liks[s:e]]
        unit_ts = [int(ts) for ts in tx_starts[s:e]]
        unit_te = [int(te) for te in tx_ends[s:e]]

        gapdh_idx = [j for j, t in enumerate(unit_t) if t in gapdh_mrna_set]
        if not gapdh_idx:
            continue

        g_t  = [unit_t[j] for j in gapdh_idx]
        g_ll = [unit_ll[j] for j in gapdh_idx]
        g_ts = [unit_ts[j] for j in gapdh_idx]
        g_te = [unit_te[j] for j in gapdh_idx]
        fid  = int(frag_ids[unit_i])

        if mane_tidx not in g_t:
            mane_absent.append({"unit_i": unit_i, "fid": fid, "n": len(g_t),
                                 "best_tidx": g_t[max(range(len(g_t)), key=lambda j: g_ll[j])]})
            continue

        n = len(g_t)
        if n == 1:
            mane_only.append({"unit_i": unit_i, "fid": fid, "ts": g_ts[0], "te": g_te[0]})
            continue

        mane_pos = g_t.index(mane_tidx)
        mane_ll  = g_ll[mane_pos]
        best_ll  = max(g_ll)
        gap = best_ll - mane_ll

        rec = {
            "unit_i": unit_i, "fid": fid,
            "mane_ll": mane_ll, "best_ll": best_ll, "gap": gap,
            "mane_ts": g_ts[mane_pos], "mane_te": g_te[mane_pos],
            "n_gapdh_cands": n,
        }
        if gap < -TOL:       mane_best.append(rec)
        elif gap <= TOL:     tied.append(rec)
        else:                mane_dominated.append(rec)

    total = len(mane_only) + len(mane_best) + len(tied) + len(mane_dominated) + len(mane_absent)
    print(f"\n=== {label} GAPDH EM unit classification ===")
    print(f"  Total units with any GAPDH candidate: {total:,}")
    print(f"  mane_only (MANE is sole candidate):   {len(mane_only):6,d}  ({100*len(mane_only)/max(total,1):.1f}%)")
    print(f"  mane_best (MANE wins uniquely):        {len(mane_best):6,d}  ({100*len(mane_best)/max(total,1):.1f}%)")
    print(f"  tied      (MANE = tied best):          {len(tied):6,d}  ({100*len(tied)/max(total,1):.1f}%)")
    print(f"  mane_dominated (MANE loses):           {len(mane_dominated):6,d}  ({100*len(mane_dominated)/max(total,1):.1f}%)")
    print(f"  mane_absent (MANE not a candidate):    {len(mane_absent):6,d}  ({100*len(mane_absent)/max(total,1):.1f}%)")

    # tx_start distribution for mane_best (what positions are informative?)
    if mane_best:
        ts_arr = np.array([r["mane_ts"] for r in mane_best])
        print(f"\n  mane_best tx_start: mean={ts_arr.mean():.0f} "
              f"min={ts_arr.min()} max={ts_arr.max()} "
              f"p25={np.percentile(ts_arr,25):.0f} p75={np.percentile(ts_arr,75):.0f}")
        bins = list(range(0, 1300, 50))
        hist, _ = np.histogram(ts_arr, bins=bins)
        print(f"  mane_best tx_start histogram (bin=50):")
        for i, cnt in enumerate(hist):
            if cnt > 0:
                print(f"    [{bins[i]:4d},{bins[i+1]:4d}): {cnt:5,d}")

    # gap distribution for mane_dominated
    if mane_dominated:
        gaps = [r["gap"] for r in mane_dominated]
        gap_counts = Counter(round(g, 4) for g in gaps)
        print(f"\n  mane_dominated gap distribution:")
        for g, cnt in sorted(gap_counts.items(), key=lambda x: x[0]):
            print(f"    gap={g:.4f}  count={cnt:,}")
        dom_ts = np.array([r["mane_ts"] for r in mane_dominated])
        print(f"  mane_dominated tx_start: mean={dom_ts.mean():.0f} "
              f"p25={np.percentile(dom_ts,25):.0f} p75={np.percentile(dom_ts,75):.0f}")

    # mane_absent: which best candidate?
    if mane_absent:
        absent_best = Counter(t_info.get(r["best_tidx"], {}).get("t_id", "?")
                              for r in mane_absent)
        print(f"\n  mane_absent: top best candidates (note: MANE not in these units):")
        for tid, cnt in absent_best.most_common(10):
            print(f"    {tid:35s}  {cnt:5,d}")

    return {
        "mane_only": len(mane_only),
        "mane_best": len(mane_best),
        "tied": len(tied),
        "mane_dominated": len(mane_dominated),
        "mane_absent": len(mane_absent),
        "total": total,
        "mane_best_data": mane_best,
        "tied_data": tied,
    }


def run_and_classify(bam_path, label):
    scan_config = BamScanConfig(include_multimap=True)
    _, strand_models, fl_models, buffer, _, _ = scan_and_buffer(bam_path, tx_index, scan=scan_config)
    em_config = EMConfig()
    scoring_config = FragmentScoringConfig()
    geometry, estimator = _setup_geometry_and_estimator(tx_index, fl_models, em_config)
    stats_obj = PipelineStats()
    sf = _score_fragments(buffer, tx_index, strand_models, fl_models, stats_obj, estimator,
                          scoring_config, log_every=500000, annotations=None)
    print(f"  {label}: {sf.n_units:,} EM units")
    result = classify_units(sf, label)
    del sf; gc.collect()
    return result


print("\n--- Running oracle ---")
oracle_result = run_and_classify(ORACLE_BAM, "ORACLE")

print("\n--- Running minimap2 ---")
mm2_result = run_and_classify(MM2_BAM, "MINIMAP2")

# ---- Side-by-side comparison ----
print(f"\n{'='*70}")
print("SIDE-BY-SIDE COMPARISON (GAPDH EM units)")
print(f"{'='*70}")
for key in ["mane_only", "mane_best", "tied", "mane_dominated", "mane_absent", "total"]:
    ov = oracle_result[key]
    mv = mm2_result[key]
    diff = mv - ov
    pct_o = 100.0 * ov / max(oracle_result["total"], 1)
    pct_m = 100.0 * mv / max(mm2_result["total"], 1)
    print(f"  {key:20s}  oracle={ov:7,d} ({pct_o:5.1f}%)   mm2={mv:7,d} ({pct_m:5.1f}%)   diff={diff:+7,d}")

print(f"\nKey insight: oracle 'mane_best' = {oracle_result['mane_best']:,} vs "
      f"mm2 'mane_best' = {mm2_result['mane_best']:,}")
print(f"  Difference: {mm2_result['mane_best'] - oracle_result['mane_best']:+,d} "
      f"units where MANE uniquely wins")
print(f"  mm2 converts {oracle_result['mane_best'] - mm2_result['mane_best']:,} oracle 'mane_best' "
      f"→ 'tied' or 'dominated'")

print("\nDONE")
