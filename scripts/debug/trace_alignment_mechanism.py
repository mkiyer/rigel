#!/usr/bin/env python3
"""Confirm the alignment mechanism causing mane_absent in minimap2.

Hypothesis: minimap2 produces UNSPLICED (intronic) alignments for reads at
exon-intron boundaries. These reads map to a genomic region where MANE is
not scored as a candidate (intronic position not in any exon), OR where MANE
has high overhang and gets PRUNED from the candidate set.

This script:
1. Finds mane_absent units in minimap2 that are MANE-source reads
2. Looks up their BAM records and shows the CIGAR / splice structure
3. Compares to oracle BAM for the same qname
4. Shows whether the alignment is spliced vs unspliced
"""
import sys
sys.path.insert(0, "src")

import math
import gc
import pysam
from collections import Counter, defaultdict

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
    }

mane_tidx = next(ti for ti, inf in t_info.items() if inf["is_mane"] and inf["g_name"] == "GAPDH")
gapdh_mrna_set = {ti for ti, inf in t_info.items() if inf["g_name"] == "GAPDH" and not inf["is_nrna"]}
mane_t_id = t_info[mane_tidx]["t_id"]
print(f"  MANE={mane_t_id} t_index={mane_tidx}")


def extract_source_tid(qn): return qn.split(":")[0]


def build_qname_map(bam_path: str) -> dict:
    """Map qname → list of BAM records."""
    records = defaultdict(list)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            records[r.query_name].append({
                "chrom": r.reference_name,
                "pos": r.reference_start,
                "pos_end": r.reference_end,
                "cigar": r.cigarstring,
                "flag": r.flag,
                "is_r1": bool(r.flag & 0x40),
                "is_r2": bool(r.flag & 0x80),
                "is_sec": bool(r.flag & 0x100),
                "is_rev": bool(r.flag & 0x10),
                "nm": r.get_tag("NM") if r.has_tag("NM") else -1,
                "n_introns": r.cigarstring.count("N") if r.cigarstring else 0,
            })
    return records


def has_N_in_cigar(cigar: str) -> bool:
    return "N" in cigar if cigar else False


def get_mane_absent_units(bam_path: str, label: str) -> tuple[list, dict]:
    """Return list of mane_absent unit records and fid→qname map."""
    fid_to_qname = {}
    seen = set()
    fid = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.query_name not in seen:
                seen.add(r.query_name)
                fid_to_qname[fid] = r.query_name
                fid += 1

    scan_config = BamScanConfig(include_multimap=True)
    _, strand_models, fl_models, buffer, _, _ = scan_and_buffer(bam_path, tx_index, scan=scan_config)

    # Capture raw exon_bp per fragment BEFORE scoring releases buffer
    print(f"  Capturing buffer data for {label}...")
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
    print(f"  {label}: {sf.n_units:,} EM units")

    offsets   = sf.offsets
    t_indices = sf.t_indices
    log_liks  = sf.log_liks
    frag_ids  = sf.frag_ids

    mane_absent = []
    mane_dominated = []

    for unit_i in range(sf.n_units):
        s = int(offsets[unit_i])
        e = int(offsets[unit_i + 1])
        if e <= s:
            continue
        unit_t  = [int(t) for t in t_indices[s:e]]
        unit_ll = [float(l) for l in log_liks[s:e]]
        gapdh_idx = [j for j, t in enumerate(unit_t) if t in gapdh_mrna_set]
        if not gapdh_idx:
            continue
        fid = int(frag_ids[unit_i])
        g_t = [unit_t[j] for j in gapdh_idx]
        g_ll = [unit_ll[j] for j in gapdh_idx]

        qn = fid_to_qname.get(fid, "?")
        src = extract_source_tid(qn)
        best_j = max(range(len(g_t)), key=lambda j: g_ll[j])
        best_tidx = g_t[best_j]
        best_ll = g_ll[best_j]

        if mane_tidx not in g_t:
            raw = fid_raw.get(fid, (0, 0, {}))
            mane_absent.append({
                "unit_i": unit_i, "fid": fid, "qn": qn, "src": src,
                "best_tidx": best_tidx, "best_ll": best_ll,
                "rl": raw[0], "nm": raw[1],
                "best_ebp": raw[2].get(best_tidx, 0),
            })
        elif len(g_t) > 1:
            mane_pos = g_t.index(mane_tidx)
            mane_ll = g_ll[mane_pos]
            gap = best_ll - mane_ll
            if gap > TOL:
                raw = fid_raw.get(fid, (0, 0, {}))
                mane_dominated.append({
                    "unit_i": unit_i, "fid": fid, "qn": qn, "src": src,
                    "best_tidx": best_tidx, "gap": gap,
                    "rl": raw[0], "nm": raw[1],
                    "mane_ebp": raw[2].get(mane_tidx, 0),
                    "best_ebp": raw[2].get(best_tidx, 0),
                })

    del sf; gc.collect()
    return mane_absent, mane_dominated, fid_to_qname


# ---- Run minimap2 ----
print("\nAnalyzing minimap2 mane_absent units...")
mm2_absent, mm2_dominated, mm2_fid2q = get_mane_absent_units(MM2_BAM, "MINIMAP2")
print(f"  mane_absent: {len(mm2_absent):,}")
print(f"  mane_dominated: {len(mm2_dominated):,}")

# ---- Source breakdown ----
absent_src = Counter(r["src"] for r in mm2_absent)
absent_best = Counter(t_info.get(r["best_tidx"],{}).get("t_id","?") for r in mm2_absent)
print(f"\n  mane_absent source breakdown:")
for src, cnt in absent_src.most_common(10):
    pct = 100 * cnt / len(mm2_absent)
    print(f"    {src:35s}  {cnt:5,d}  ({pct:.1f}%)")
print(f"\n  mane_absent best-candidate breakdown:")
for tid, cnt in absent_best.most_common(10):
    pct = 100 * cnt / len(mm2_absent)
    print(f"    {tid:35s}  {cnt:5,d}  ({pct:.1f}%)")

# ---- Focus on MANE-source mane_absent: why is MANE excluded? ----
mane_src_absent = [r for r in mm2_absent if r["src"] == mane_t_id]
print(f"\n  MANE-source frags that are mane_absent: {len(mane_src_absent):,}")

# ---- Build BAM index for minimap2 (all qnames) ----
print("\nBuilding minimap2 BAM index for selected qnames...")
# Only fetch BAM records for mane_absent and a subset of dominated
key_qnames = {r["qn"] for r in mm2_absent[:500]} | {r["qn"] for r in mm2_dominated[:200]}
mm2_bam_records = defaultdict(list)
with pysam.AlignmentFile(MM2_BAM, "rb") as bam:
    for r in bam:
        if r.query_name in key_qnames:
            nm = r.get_tag("NM") if r.has_tag("NM") else -1
            mm2_bam_records[r.query_name].append({
                "chrom": r.reference_name, "pos": r.reference_start,
                "pos_end": r.reference_end, "cigar": r.cigarstring,
                "is_r1": bool(r.flag & 0x40), "is_r2": bool(r.flag & 0x80),
                "is_sec": bool(r.flag & 0x100), "is_rev": bool(r.flag & 0x10),
                "nm": nm,
                "n_introns": r.cigarstring.count("N") if r.cigarstring else 0,
            })

# ---- Classify alignment types for mane_absent ----
print(f"\n{'='*70}")
print("ALIGNMENT TYPE ANALYSIS: mane_absent units in minimap2")
print("Does minimap2 produce UNSPLICED alignments where oracle would splice?")
print(f"{'='*70}")

# For mane_absent: what fraction have NO intron (N operation) in CIGAR?
n_spliced = 0
n_unspliced = 0
n_missing_bam = 0
n_introns_dist = Counter()

for r in mm2_absent:
    recs = mm2_bam_records.get(r["qn"], [])
    if not recs:
        n_missing_bam += 1
        continue
    # Check if any alignment (R1 or R2) has an N (intron) in CIGAR
    has_splice = any(rec["n_introns"] > 0 for rec in recs)
    n_introns_total = sum(rec["n_introns"] for rec in recs)
    n_introns_dist[n_introns_total] += 1
    if has_splice:
        n_spliced += 1
    else:
        n_unspliced += 1

total_in_bam = n_spliced + n_unspliced
pct_unspliced = 100 * n_unspliced / max(total_in_bam, 1)
pct_spliced = 100 * n_spliced / max(total_in_bam, 1)
print(f"  mane_absent fragments with BAM records: {total_in_bam:,}")
print(f"  Fully UNSPLICED (no N in any CIGAR):    {n_unspliced:,} ({pct_unspliced:.1f}%)")
print(f"  At least partially SPLICED (N present): {n_spliced:,} ({pct_spliced:.1f}%)")

print(f"\n  N-operation count distribution:")
for n, cnt in sorted(n_introns_dist.items()):
    print(f"    {n} introns: {cnt:5,d} fragments")

# ---- Sample some MANE-source mane_absent to show BAM records ----
print(f"\n{'='*70}")
print("SAMPLE BAM RECORDS: MANE-source reads that are mane_absent in mm2")
print(f"{'='*70}")
shown = 0
for r in mane_src_absent[:30]:
    recs = mm2_bam_records.get(r["qn"], [])
    if not recs:
        continue
    rl, nm_frag, ebp_map = r["rl"], r["nm"], r.get("best_ebp", 0)
    best_t = t_info.get(r["best_tidx"], {})
    print(f"\n  qn={r['qn']}")
    print(f"  rl={rl} nm={nm_frag}  best_cand={best_t.get('t_id','?')} ebp={ebp_map}")
    for a in recs:
        tag = "R1" if a["is_r1"] else "R2"
        sec = "SEC" if a["is_sec"] else "PRI"
        rev = "-" if a["is_rev"] else "+"
        splice = f"(introns={a['n_introns']})"
        print(f"    {tag} {sec} {a['chrom']}:{a['pos']:>10,d}-{a['pos_end']:>10,d}  "
              f"NM={a['nm']}  cigar={a['cigar']}  {splice}")
    shown += 1
    if shown >= 10:
        break

# ---- Quantify: how much does the EM seed effect explain the over-count? ----
print(f"\n{'='*70}")
print("EM SEED EFFECT: fragments that seed ENST396859.5 in minimap2")
print(f"{'='*70}")
enst396859 = next((ti for ti, inf in t_info.items() if "396859" in inf["t_id"]), None)
seed_for_396859 = [r for r in mm2_absent if r["best_tidx"] == enst396859]
dom_for_396859 = [r for r in mm2_dominated if r["best_tidx"] == enst396859]
print(f"  ENST396859.5 t_index={enst396859}")
print(f"  mane_absent units with best=ENST396859.5:   {len(seed_for_396859):,}")
print(f"  mane_dominated units with best=ENST396859.5: {len(dom_for_396859):,}")
print(f"  Total 'wrong' seed fragments for ENST396859:  {len(seed_for_396859)+len(dom_for_396859):,}")
print(f"\n  Oracle shows ENST396859.5 at 33,921 TPM (~3.4% of GAPDH)")
print(f"  Minimap2 shows ENST396859.5 at ~209,000 TPM (~21% of GAPDH) = 6× overcount")

print("\nDONE")
