#!/usr/bin/env python3
"""Root cause analysis: WHY does minimap2 fail where oracle succeeds?

Hypothesis: The 96.3% tied EM units are NOT the issue. The problem is:
  (A) MANE-source fragments that ARE present in oracle EM but ABSENT from
      minimap2 GAPDH EM (aligned to pseudogenes instead)
  (B) Non-GAPDH-source (pseudogene) fragments that ENTER the minimap2
      GAPDH EM and inflate shorter isoforms

Since simulated fragment qnames encode source transcript (e.g.
ENST00000229239.10:start-end:strand:id), we can trace exactly:
  - For each GAPDH EM unit: what is the source transcript?
  - How many MANE-source fragments HAVE MANE as a candidate?
  - How many MANE-source fragments are ABSENT from GAPDH EM entirely?
  - What is in the GAPDH EM in minimap2 that shouldn't be there?
"""
import sys
sys.path.insert(0, "src")

import math
import gc
import numpy as np
import pysam
from collections import defaultdict, Counter
from pathlib import Path

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
        "t_id": r["t_id"],
        "g_name": r.get("g_name", "?"),
        "is_mane": bool(r.get("is_mane", False)),
        "is_nrna": bool(r.get("is_synthetic_nrna", False)),
        "chrom": r.get("chrom", "?"),
        "start": int(r.get("start", 0)),
        "end": int(r.get("end", 0)),
    }

mane_tidx = next(ti for ti, inf in t_info.items() if inf["is_mane"] and inf["g_name"] == "GAPDH")
gapdh_mrna_set = {ti for ti, inf in t_info.items() if inf["g_name"] == "GAPDH" and not inf["is_nrna"]}
tid_to_tidx = {inf["t_id"]: ti for ti, inf in t_info.items()}
print(f"  MANE t_index={mane_tidx} ({t_info[mane_tidx]['t_id']})")


def extract_source_tid(qname: str) -> str:
    """Extract source transcript ID from simulated qname like ENST00000229239.10:start-end:strand:id"""
    return qname.split(":")[0]


def get_gapdh_source_fids(bam_path: str) -> tuple[dict, set, set]:
    """Returns:
    - fid_to_qname: all fragment ids mapped to qname
    - gapdh_fids: fragment ids where source is any GAPDH mRNA transcript
    - mane_fids: fragment ids where source is MANE specifically
    """
    gapdh_t_ids = {t_info[ti]["t_id"] for ti in gapdh_mrna_set}
    fid_to_qname = {}
    gapdh_fids = set()
    mane_fids = set()
    seen = set()
    fid = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            qn = r.query_name
            if qn not in seen:
                seen.add(qn)
                fid_to_qname[fid] = qn
                src = extract_source_tid(qn)
                if src in gapdh_t_ids:
                    gapdh_fids.add(fid)
                if src == t_info[mane_tidx]["t_id"]:
                    mane_fids.add(fid)
                fid += 1
    return fid_to_qname, gapdh_fids, mane_fids


def score_bam(bam_path: str, label: str):
    """Score a BAM and return (sf, fid_to_qname, gapdh_fids, mane_fids)."""
    print(f"\n{'='*70}")
    print(f"Scanning + scoring {label}")
    print(f"{'='*70}")

    fid_to_qname, gapdh_fids, mane_fids = get_gapdh_source_fids(bam_path)
    print(f"  Total fragments: {len(fid_to_qname):,}")
    print(f"  MANE-source fragments: {len(mane_fids):,}")
    print(f"  Any-GAPDH-source fragments: {len(gapdh_fids):,}")

    scan_config = BamScanConfig(include_multimap=True)
    _, strand_models, fl_models, buffer, _, _ = scan_and_buffer(bam_path, tx_index, scan=scan_config)
    em_config = EMConfig()
    scoring_config = FragmentScoringConfig()
    geometry, estimator = _setup_geometry_and_estimator(tx_index, fl_models, em_config)
    stats_obj = PipelineStats()
    sf = _score_fragments(buffer, tx_index, strand_models, fl_models, stats_obj, estimator,
                          scoring_config, log_every=500000, annotations=None)
    print(f"  EM units: {sf.n_units:,}")
    return sf, fid_to_qname, gapdh_fids, mane_fids


def analyze_gapdh_em(sf, fid_to_qname, gapdh_fids, mane_fids, label: str):
    """Analyze EM unit distribution for GAPDH fragments."""
    offsets   = sf.offsets
    t_indices = sf.t_indices
    log_liks  = sf.log_liks
    frag_ids  = sf.frag_ids

    # Categorize each EM unit
    # For each EM unit containing any GAPDH mRNA t_index:
    #   - What is the source transcript of the fragment?
    #   - Is MANE a candidate?
    #   - Which GAPDH isoforms are candidates?

    unit_has_gapdh = []  # (fid, src_tid, has_mane, n_gapdh_cands, best_tidx)
    fids_in_gapdh_em = set()

    for unit_i in range(sf.n_units):
        s = int(offsets[unit_i])
        e = int(offsets[unit_i + 1])
        if e <= s:
            continue
        unit_t  = [int(t) for t in t_indices[s:e]]
        unit_ll = [float(ll) for ll in log_liks[s:e]]
        gapdh_idx = [j for j, t in enumerate(unit_t) if t in gapdh_mrna_set]
        if not gapdh_idx:
            continue
        fid = int(frag_ids[unit_i])
        fids_in_gapdh_em.add(fid)
        g_t  = [unit_t[j] for j in gapdh_idx]
        g_ll = [unit_ll[j] for j in gapdh_idx]
        has_mane = mane_tidx in g_t
        best_j = max(range(len(g_t)), key=lambda j: g_ll[j])
        best_tidx = g_t[best_j]
        src_qname = fid_to_qname.get(fid, "?")
        src_tid = extract_source_tid(src_qname)
        unit_has_gapdh.append({
            "fid": fid, "src_tid": src_tid, "has_mane": has_mane,
            "n_gapdh_cands": len(g_t), "best_tidx": best_tidx,
        })

    total_gapdh_units = len(unit_has_gapdh)
    mane_source_in_gapdh_em = sum(1 for u in unit_has_gapdh
                                   if extract_source_tid(fid_to_qname.get(u["fid"],"?")) == t_info[mane_tidx]["t_id"])
    mane_cand_in_gapdh_em = sum(1 for u in unit_has_gapdh if u["has_mane"])

    print(f"\n--- {label}: GAPDH EM analysis ---")
    print(f"  EM units with any GAPDH mRNA candidate: {total_gapdh_units:,}")
    print(f"  Unique fragment IDs in GAPDH EM:         {len(fids_in_gapdh_em):,}")
    print(f"  MANE-source frags in GAPDH EM:           {mane_source_in_gapdh_em:,}  "
          f"(of {len(mane_fids):,} total)")
    print(f"  MISSING MANE-source frags from GAPDH EM: {len(mane_fids - fids_in_gapdh_em):,}")
    print(f"  EM units where MANE IS a candidate:      {mane_cand_in_gapdh_em:,}")

    # Source breakdown for ALL units in GAPDH EM
    src_counts = Counter(u["src_tid"] for u in unit_has_gapdh)
    print(f"\n  Source transcript breakdown (top 15):")
    for src, cnt in src_counts.most_common(15):
        pct = 100.0 * cnt / total_gapdh_units
        gapdh_src = src in {t_info[ti]["t_id"] for ti in gapdh_mrna_set}
        marker = "  ← GAPDH" if gapdh_src else "  ← PSEUDOGENE/OTHER"
        print(f"    {src:35s} {cnt:6,d}  ({pct:5.1f}%){marker}")

    # Where do MANE-source fragments end up in this BAM?
    mane_src_units = [u for u in unit_has_gapdh
                      if extract_source_tid(fid_to_qname.get(u["fid"],"?")) == t_info[mane_tidx]["t_id"]]
    mane_missing_from_em = mane_fids - fids_in_gapdh_em
    print(f"\n  MANE-source fragments in GAPDH EM: {len(mane_src_units):,}")
    print(f"  Of those with MANE as candidate:   {sum(1 for u in mane_src_units if u['has_mane']):,}")
    print(f"  Of those WITHOUT MANE as candidate:{sum(1 for u in mane_src_units if not u['has_mane']):,}")
    if mane_src_units:
        best_counts = Counter(t_info.get(u["best_tidx"], {}).get("t_id", "?")
                              for u in mane_src_units if not u["has_mane"])
        print(f"  Where MANE-source goes when MANE absent:")
        for tid, cnt in best_counts.most_common(10):
            print(f"    best={tid:35s}  {cnt:5,d}")

    # Non-GAPDH source fragments that show up in GAPDH EM (spurious)
    spurious = [u for u in unit_has_gapdh
                if extract_source_tid(fid_to_qname.get(u["fid"],"?")) not in
                   {t_info[ti]["t_id"] for ti in gapdh_mrna_set}]
    print(f"\n  Non-GAPDH-source fragments in GAPDH EM (spurious): {len(spurious):,}")
    if spurious:
        spur_src = Counter(u["src_tid"] for u in spurious)
        print(f"  Top sources of spurious fragments:")
        for src, cnt in spur_src.most_common(10):
            print(f"    {src:35s}  {cnt:5,d}")
        spur_best = Counter(t_info.get(u["best_tidx"], {}).get("t_id", "?") for u in spurious)
        print(f"  Which GAPDH isoforms get the spurious fragments:")
        for tid, cnt in spur_best.most_common(10):
            print(f"    best_cand={tid:35s}  {cnt:5,d}")

    return fids_in_gapdh_em, mane_source_in_gapdh_em, len(spurious)


# ---- Run both BAMs ----
sf_oracle, fid2q_oracle, gapdh_fids_oracle, mane_fids_oracle = score_bam(ORACLE_BAM, "ORACLE")
em_fids_oracle, mane_in_oracle, spur_oracle = analyze_gapdh_em(
    sf_oracle, fid2q_oracle, gapdh_fids_oracle, mane_fids_oracle, "ORACLE")

del sf_oracle; gc.collect()

sf_mm2, fid2q_mm2, gapdh_fids_mm2, mane_fids_mm2 = score_bam(MM2_BAM, "MINIMAP2")
em_fids_mm2, mane_in_mm2, spur_mm2 = analyze_gapdh_em(
    sf_mm2, fid2q_mm2, gapdh_fids_mm2, mane_fids_mm2, "MINIMAP2")

# ---- Summary comparison ----
print(f"\n{'='*70}")
print("COMPARISON SUMMARY")
print(f"{'='*70}")
print(f"{'Metric':45s}  {'Oracle':>10}  {'Minimap2':>10}  {'Diff':>8}")
print("-" * 80)
print(f"{'MANE-source frags total':45s}  {len(mane_fids_oracle):>10,}  {len(mane_fids_mm2):>10,}")
print(f"{'MANE-source frags in GAPDH EM':45s}  {mane_in_oracle:>10,}  {mane_in_mm2:>10,}  "
      f"  {mane_in_mm2 - mane_in_oracle:>+6,}")
print(f"{'MANE frags missing from GAPDH EM':45s}  "
      f"{len(mane_fids_oracle - em_fids_oracle):>10,}  "
      f"{len(mane_fids_mm2 - em_fids_mm2):>10,}  "
      f"  {(len(mane_fids_mm2 - em_fids_mm2)) - (len(mane_fids_oracle - em_fids_oracle)):>+6,}")
print(f"{'Spurious (non-GAPDH) in GAPDH EM':45s}  {spur_oracle:>10,}  {spur_mm2:>10,}  "
      f"  {spur_mm2 - spur_oracle:>+6,}")

print("\nDONE")
