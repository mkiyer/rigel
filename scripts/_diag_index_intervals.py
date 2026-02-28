#!/usr/bin/env python3
"""Diagnostic: Analyze the HulkIndex interval structure.

Reports:
1. Total intervals by type (EXON, INTRON, INTERGENIC)
2. Boundary-based redundancy: how many intervals share exact (ref, start, end)
3. Per-transcript interval counts
4. Simulated query hit distribution (sample genomic positions)
5. Equivalence class potential: intervals with identical (ref, start, end)
   but different t_index values
"""

import sys
import os
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from hulkrna.index import HulkIndex, INTERVALS_FEATHER, TRANSCRIPTS_FEATHER
from hulkrna.types import IntervalType


INDEX_DIR = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/hulkrna_index")


def main():
    print(f"Loading index from {INDEX_DIR}")
    index = HulkIndex.load(INDEX_DIR)

    # Read raw interval table
    iv_df = pd.read_feather(os.path.join(str(INDEX_DIR), INTERVALS_FEATHER))
    t_df = pd.read_feather(os.path.join(str(INDEX_DIR), TRANSCRIPTS_FEATHER))

    print(f"\n{'='*70}")
    print(f"INDEX SUMMARY")
    print(f"{'='*70}")
    print(f"Transcripts:  {len(t_df)}")
    print(f"Genes:        {index.num_genes}")
    print(f"Total intervals: {len(iv_df)}")

    # --- 1. Interval counts by type ---
    print(f"\n--- Interval counts by type ---")
    type_counts = Counter(iv_df["interval_type"].values)
    for itype_val, count in sorted(type_counts.items()):
        try:
            name = IntervalType(itype_val).name
        except ValueError:
            name = f"UNKNOWN({itype_val})"
        print(f"  {name:15s}: {count:>8,d}")

    # --- 2. Interval sizes by type ---
    print(f"\n--- Interval sizes (bp) by type ---")
    iv_df["size"] = iv_df["end"] - iv_df["start"]
    for itype_val in sorted(type_counts.keys()):
        try:
            name = IntervalType(itype_val).name
        except ValueError:
            name = f"UNKNOWN({itype_val})"
        mask = iv_df["interval_type"] == itype_val
        sizes = iv_df.loc[mask, "size"]
        print(f"  {name:15s}: median={sizes.median():>8,.0f}  "
              f"mean={sizes.mean():>8,.0f}  "
              f"min={sizes.min():>6,d}  max={sizes.max():>10,d}  "
              f"total={sizes.sum():>12,d}")

    # --- 3. Per-transcript interval counts ---
    print(f"\n--- Per-transcript interval counts ---")
    genic_mask = iv_df["t_index"] >= 0
    genic = iv_df[genic_mask]
    t_iv_counts = genic.groupby("t_index").size()
    print(f"  Transcripts with intervals: {len(t_iv_counts)}")
    print(f"  Intervals per transcript: "
          f"median={t_iv_counts.median():.0f}  "
          f"mean={t_iv_counts.mean():.1f}  "
          f"min={t_iv_counts.min()}  max={t_iv_counts.max()}")

    # Breakdown by type per transcript
    for itype_val in sorted(type_counts.keys()):
        if itype_val == int(IntervalType.INTERGENIC):
            continue
        try:
            name = IntervalType(itype_val).name
        except ValueError:
            name = f"UNKNOWN({itype_val})"
        mask2 = genic["interval_type"] == itype_val
        per_t = genic[mask2].groupby("t_index").size()
        if len(per_t) > 0:
            print(f"  {name:15s}: per-transcript median={per_t.median():.0f}  "
                  f"mean={per_t.mean():.1f}  max={per_t.max()}")

    # --- 4. Boundary redundancy: exact (ref, start, end) duplicates ---
    print(f"\n--- Boundary redundancy (exact ref, start, end matches) ---")
    genic_keys = genic[["ref", "start", "end"]].apply(tuple, axis=1)
    key_counts = genic_keys.value_counts()
    total_genic = len(genic)
    unique_boundaries = len(key_counts)
    max_dup = key_counts.max()
    mean_dup = key_counts.mean()
    print(f"  Total genic intervals:    {total_genic:>8,d}")
    print(f"  Unique (ref,start,end):   {unique_boundaries:>8,d}")
    print(f"  Redundancy ratio:         {total_genic/unique_boundaries:.2f}x")
    print(f"  Max duplicates at one boundary: {max_dup}")
    print(f"  Mean duplicates per boundary:   {mean_dup:.2f}")

    # Distribution of duplicate counts
    dup_dist = key_counts.value_counts().sort_index()
    print(f"  Duplicate-count distribution:")
    for count, n_boundaries in dup_dist.head(20).items():
        print(f"    {count} transcript(s) share boundary: {n_boundaries} unique intervals")

    # --- 4b. With type: exact (ref, start, end, type) duplicates ---
    print(f"\n--- Boundary redundancy with type (ref, start, end, type) ---")
    genic_keys_typed = genic[["ref", "start", "end", "interval_type"]].apply(tuple, axis=1)
    key_typed_counts = genic_keys_typed.value_counts()
    unique_typed = len(key_typed_counts)
    print(f"  Unique (ref,start,end,type): {unique_typed:>8,d}")
    print(f"  Redundancy ratio:            {total_genic/unique_typed:.2f}x")

    # --- 5. Equivalence classes: same (ref, start, end) → set of t_indices ---
    print(f"\n--- Equivalence classes (boundary → transcript set) ---")
    # Group genic exon intervals by (ref, start, end)
    exon_mask = genic["interval_type"] == int(IntervalType.EXON)
    exon_genic = genic[exon_mask]
    exon_groups = exon_genic.groupby(["ref", "start", "end"])["t_index"].apply(frozenset)

    # How many unique transcript sets?
    unique_t_sets = set(exon_groups.values)
    print(f"  Unique exon boundaries:     {len(exon_groups):>8,d}")
    print(f"  Unique transcript sets:     {len(unique_t_sets):>8,d}")

    # Distribution of transcript set sizes
    set_sizes = [len(s) for s in exon_groups.values]
    set_size_counter = Counter(set_sizes)
    print(f"  Transcript-set size distribution (exon boundaries):")
    for size, count in sorted(set_size_counter.items())[:15]:
        print(f"    size {size}: {count} boundaries")

    # Same for intron intervals
    intron_mask = genic["interval_type"] == int(IntervalType.INTRON)
    intron_genic = genic[intron_mask]
    intron_groups = intron_genic.groupby(["ref", "start", "end"])["t_index"].apply(frozenset)
    unique_intron_t_sets = set(intron_groups.values)
    print(f"\n  Unique intron boundaries:   {len(intron_groups):>8,d}")
    print(f"  Unique intron transcript sets: {len(unique_intron_t_sets):>8,d}")

    intron_set_sizes = [len(s) for s in intron_groups.values]
    intron_size_counter = Counter(intron_set_sizes)
    print(f"  Transcript-set size distribution (intron boundaries):")
    for size, count in sorted(intron_size_counter.items())[:15]:
        print(f"    size {size}: {count} boundaries")

    # --- 6. Simulated query hit distribution ---
    print(f"\n--- Simulated query hit distribution (100bp windows) ---")
    # Sample 10k random positions from the first reference
    refs = iv_df["ref"].unique()
    ref = refs[0]
    ref_mask = iv_df["ref"] == ref
    ref_start = iv_df.loc[ref_mask, "start"].min()
    ref_end = iv_df.loc[ref_mask, "end"].max()
    print(f"  Reference: {ref}, span: {ref_start:,d} - {ref_end:,d}")

    np.random.seed(42)
    n_queries = 10000
    query_starts = np.random.randint(ref_start, max(ref_end - 100, ref_start + 1), size=n_queries)
    hit_counts = []
    exon_hit_counts = []
    for qs in query_starts:
        hits = list(index.cr.overlap(ref, int(qs), int(qs + 100)))
        hit_counts.append(len(hits))
        n_exon = sum(1 for _, _, label in hits if int(index._iv_type[label]) == int(IntervalType.EXON))
        exon_hit_counts.append(n_exon)

    hit_arr = np.array(hit_counts)
    exon_arr = np.array(exon_hit_counts)
    print(f"  Total hits per query (100bp): "
          f"median={np.median(hit_arr):.0f}  "
          f"mean={hit_arr.mean():.1f}  "
          f"min={hit_arr.min()}  max={hit_arr.max()}")
    print(f"  Exon hits per query:          "
          f"median={np.median(exon_arr):.0f}  "
          f"mean={exon_arr.mean():.1f}  "
          f"min={exon_arr.min()}  max={exon_arr.max()}")

    # Distribution
    print(f"  Hit count distribution:")
    for bucket in [0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500]:
        n = (hit_arr <= bucket).sum()
        print(f"    ≤{bucket:>3d}: {n:>6d} ({100*n/len(hit_arr):>5.1f}%)")

    # --- 7. Transcript span analysis ---
    print(f"\n--- Transcript span analysis ---")
    # For each transcript, compute the span (end - start)
    # vs number of intervals (exons + introns)
    t_spans = t_df["end"] - t_df["start"]
    t_n_exons = exon_genic.groupby("t_index").size()
    t_n_introns = intron_genic.groupby("t_index").size()

    # Single-exon transcripts (have 0 introns)
    n_single_exon = (t_n_exons == 1).sum()
    n_multi_exon = (t_n_exons > 1).sum()
    print(f"  Single-exon transcripts: {n_single_exon}")
    print(f"  Multi-exon transcripts:  {n_multi_exon}")

    # What fraction of intervals come from multi-exon transcripts?
    multi_exon_tids = set(t_n_exons[t_n_exons > 1].index)
    n_intervals_multi = len(genic[genic["t_index"].isin(multi_exon_tids)])
    print(f"  Intervals from multi-exon: {n_intervals_multi} ({100*n_intervals_multi/total_genic:.1f}%)")

    # Intron intervals as fraction of total genic
    n_intron_total = len(intron_genic)
    n_exon_total = len(exon_genic)
    print(f"  Exon intervals:   {n_exon_total:>6,d} ({100*n_exon_total/total_genic:.1f}%)")
    print(f"  Intron intervals: {n_intron_total:>6,d} ({100*n_intron_total/total_genic:.1f}%)")

    # --- 8. Potential savings from collapsing ---
    print(f"\n--- Potential savings from interval collapsing ---")
    # If we stored each unique (ref, start, end, type) once
    # with a list of t_indices, how many intervals would we have?
    print(f"  Current interval count (total): {len(iv_df):>8,d}")
    print(f"  If collapsed to unique (ref,start,end,type): {unique_typed:>8,d}")
    print(f"  Savings: {len(iv_df) - unique_typed:>8,d} intervals "
          f"({100*(len(iv_df) - unique_typed)/len(iv_df):.1f}%)")

    # If we stored transcripts as spans instead of per-intron
    # Each transcript would contribute: 1 span interval + N exon intervals
    # instead of N exon intervals + (N-1) intron intervals
    n_current_genic = total_genic
    n_span_model = n_exon_total + len(t_n_exons)  # exons + 1 span per transcript
    print(f"\n  Span model (exons + 1 span per transcript):")
    print(f"    Current genic intervals: {n_current_genic:>8,d}")
    print(f"    Span model intervals:    {n_span_model:>8,d}")
    print(f"    Savings: {n_current_genic - n_span_model:>8,d} "
          f"({100*(n_current_genic - n_span_model)/n_current_genic:.1f}%)")


if __name__ == "__main__":
    main()
