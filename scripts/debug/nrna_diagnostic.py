#!/usr/bin/env python3
"""nRNA Redesign Diagnostic: Evaluate tolerance-based merging on a real index.

Loads a rigel index and computes statistics about nRNA span merging at
multiple tolerance levels. This informs the nRNA redesign plan by
quantifying how many synthetic nRNA transcripts are needed vs how many
are already covered by annotated single-exon equivalents.

Usage:
    conda activate rigel
    python scripts/debug/nrna_diagnostic.py /path/to/rigel_index
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


def load_index_data(index_dir: str) -> pd.DataFrame:
    """Load transcript DataFrame from feather file."""
    t_df = pd.read_feather(Path(index_dir) / "transcripts.feather")
    # Derive is_single_exon: exonic length == genomic span
    t_df["span"] = t_df["end"] - t_df["start"]
    t_df["is_single_exon"] = t_df["length"] == t_df["span"]
    return t_df


def cluster_coordinates(coords: np.ndarray, tolerance: int) -> np.ndarray:
    """Cluster sorted coordinates within tolerance, returning representative for each.

    For starts: representative = min of cluster (outer left envelope).
    For ends: representative = max of cluster (outer right envelope).

    Returns array of same length as coords with cluster representative index.
    """
    n = len(coords)
    if n == 0:
        return np.empty(0, dtype=np.int64)
    cluster_ids = np.empty(n, dtype=np.int64)
    cluster_id = 0
    cluster_start_idx = 0
    for i in range(n):
        if coords[i] - coords[cluster_start_idx] > tolerance:
            cluster_id += 1
            cluster_start_idx = i
        cluster_ids[i] = cluster_id
    return cluster_ids


def compute_merged_spans(
    multi_exon_df: pd.DataFrame,
    tolerance: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Cluster TSS/TES within tolerance and compute merged nRNA spans.

    Returns:
        merged_df: DataFrame with columns [ref, strand, start, end, n_transcripts]
        t_to_merged: array mapping each row in multi_exon_df to its merged span index
    """
    n = len(multi_exon_df)
    if n == 0:
        return pd.DataFrame(columns=["ref", "strand", "start", "end", "n_transcripts"]), np.array([], dtype=np.int32)

    # Work group by group: (ref, strand)
    t_merged_start = np.empty(n, dtype=np.int64)
    t_merged_end = np.empty(n, dtype=np.int64)

    for (_ref, _strand), grp in multi_exon_df.groupby(["ref", "strand"]):
        idx = grp.index  # positions in multi_exon_df

        # Cluster starts
        starts = grp["start"].values
        order_s = np.argsort(starts, kind="mergesort")
        sorted_starts = starts[order_s]
        cluster_ids_s = cluster_coordinates(sorted_starts, tolerance)
        # Representative start = min per cluster (already the first in each cluster)
        # Build cluster_id → min_start
        unique_cids_s = np.unique(cluster_ids_s)
        cid_to_min_start = {}
        for cid in unique_cids_s:
            mask = cluster_ids_s == cid
            cid_to_min_start[cid] = sorted_starts[mask].min()
        # Map back to original order
        rep_starts = np.empty(len(idx), dtype=np.int64)
        for i, oi in enumerate(order_s):
            rep_starts[oi] = cid_to_min_start[cluster_ids_s[i]]

        # Cluster ends
        ends = grp["end"].values
        order_e = np.argsort(ends, kind="mergesort")
        sorted_ends = ends[order_e]
        cluster_ids_e = cluster_coordinates(sorted_ends, tolerance)
        unique_cids_e = np.unique(cluster_ids_e)
        cid_to_max_end = {}
        for cid in unique_cids_e:
            mask = cluster_ids_e == cid
            cid_to_max_end[cid] = sorted_ends[mask].max()
        rep_ends = np.empty(len(idx), dtype=np.int64)
        for i, oi in enumerate(order_e):
            rep_ends[oi] = cid_to_max_end[cluster_ids_e[i]]

        for j, pos in enumerate(idx):
            t_merged_start[pos] = rep_starts[j]
            t_merged_end[pos] = rep_ends[j]

    # Build merged span key for each transcript
    refs = multi_exon_df["ref"].values
    strands = multi_exon_df["strand"].values
    keys = list(zip(refs, strands, t_merged_start, t_merged_end))

    # Dedup merged spans
    span_to_idx: dict[tuple, int] = {}
    merged_rows: list[dict] = []
    t_to_merged = np.empty(n, dtype=np.int32)

    for i, key in enumerate(keys):
        if key not in span_to_idx:
            span_to_idx[key] = len(merged_rows)
            merged_rows.append({
                "ref": key[0],
                "strand": key[1],
                "start": key[2],
                "end": key[3],
                "n_transcripts": 0,
            })
        midx = span_to_idx[key]
        merged_rows[midx]["n_transcripts"] += 1
        t_to_merged[i] = midx

    merged_df = pd.DataFrame(merged_rows)
    return merged_df, t_to_merged


def find_annotated_equivalents(
    merged_df: pd.DataFrame,
    single_exon_df: pd.DataFrame,
) -> tuple[np.ndarray, list[dict]]:
    """Find merged spans covered by annotated single-exon transcripts.

    Returns:
        covered: boolean array[len(merged_df)] — True if an annotated equiv exists
        examples: list of dicts with example matches
    """
    n_merged = len(merged_df)
    covered = np.zeros(n_merged, dtype=bool)
    examples = []

    if len(single_exon_df) == 0 or n_merged == 0:
        return covered, examples

    # Build lookup: (ref, strand) → sorted list of (start, end, g_name, t_id)
    se_by_loc: dict[tuple, list] = defaultdict(list)
    for _, row in single_exon_df.iterrows():
        se_by_loc[(row["ref"], row["strand"])].append(
            (row["start"], row["end"], row.get("g_name", ""), row.get("t_id", ""))
        )
    # Sort each group by start
    for key in se_by_loc:
        se_by_loc[key].sort()

    for midx, mrow in merged_df.iterrows():
        m_ref, m_strand = mrow["ref"], mrow["strand"]
        m_start, m_end = mrow["start"], mrow["end"]
        candidates = se_by_loc.get((m_ref, m_strand), [])
        if not candidates:
            continue

        # Phase 1: exact match
        for s_start, s_end, s_gname, s_tid in candidates:
            if s_start == m_start and s_end == m_end:
                covered[midx] = True
                if len(examples) < 20:
                    examples.append({
                        "type": "exact",
                        "merged_span": f"{m_ref}:{'+' if m_strand == 1 else '-'}:{m_start}-{m_end}",
                        "single_exon_tid": s_tid,
                        "gene_name": s_gname,
                    })
                break

        if covered[midx]:
            continue

        # Phase 2: full containment (S.start <= M.start AND M.end <= S.end)
        for s_start, s_end, s_gname, s_tid in candidates:
            if s_start > m_start:
                break  # sorted by start, no more candidates
            if s_start <= m_start and m_end <= s_end:
                covered[midx] = True
                if len(examples) < 20:
                    examples.append({
                        "type": "containment",
                        "merged_span": f"{m_ref}:{'+' if m_strand == 1 else '-'}:{m_start}-{m_end}",
                        "single_exon_tid": s_tid,
                        "gene_name": s_gname,
                        "s_span": f"{s_start}-{s_end}",
                    })
                break

    return covered, examples


def analyze_tolerance(
    t_df: pd.DataFrame,
    tolerance: int,
) -> dict:
    """Run full analysis at a given tolerance level."""
    multi_exon = t_df[~t_df["is_single_exon"]].reset_index(drop=True)
    single_exon = t_df[t_df["is_single_exon"]].reset_index(drop=True)

    n_multi = len(multi_exon)
    n_single = len(single_exon)

    # Current exact dedup (tolerance=0 baseline for comparison)
    exact_spans = multi_exon.groupby(["ref", "strand", "start", "end"]).size()
    n_exact_spans = len(exact_spans)

    # Tolerance-based merging
    merged_df, t_to_merged = compute_merged_spans(multi_exon, tolerance)
    n_merged_spans = len(merged_df)

    # Find annotated equivalents
    covered, examples = find_annotated_equivalents(merged_df, single_exon)
    n_covered = covered.sum()
    n_synthetic = n_merged_spans - n_covered

    # Cluster size distribution
    if len(merged_df) > 0:
        nt = merged_df["n_transcripts"].values
        hist_1 = (nt == 1).sum()
        hist_2_5 = ((nt >= 2) & (nt <= 5)).sum()
        hist_6_20 = ((nt >= 6) & (nt <= 20)).sum()
        hist_20p = (nt > 20).sum()
    else:
        hist_1 = hist_2_5 = hist_6_20 = hist_20p = 0

    return {
        "tolerance": tolerance,
        "n_multi_exon_tx": n_multi,
        "n_single_exon_tx": n_single,
        "n_exact_spans": n_exact_spans,
        "n_merged_spans": n_merged_spans,
        "reduction_pct": (1 - n_merged_spans / max(n_exact_spans, 1)) * 100,
        "n_annotated_equiv": n_covered,
        "n_synthetic_needed": n_synthetic,
        "cluster_1": hist_1,
        "cluster_2_5": hist_2_5,
        "cluster_6_20": hist_6_20,
        "cluster_20p": hist_20p,
        "examples": examples,
    }


def main():
    parser = argparse.ArgumentParser(
        description="nRNA redesign diagnostic: evaluate tolerance-based merging"
    )
    parser.add_argument("index_dir", help="Path to rigel index directory")
    parser.add_argument(
        "--tolerances",
        type=int,
        nargs="+",
        default=[0, 5, 10, 20, 50, 100],
        help="Tolerance levels in bp (default: 0 5 10 20 50 100)",
    )
    args = parser.parse_args()

    print(f"Loading index from {args.index_dir}...")
    t_df = load_index_data(args.index_dir)

    n_total = len(t_df)
    n_single = t_df["is_single_exon"].sum()
    n_multi = n_total - n_single
    print(f"  Total transcripts: {n_total:,}")
    print(f"  Single-exon:       {n_single:,} ({n_single/n_total*100:.1f}%)")
    print(f"  Multi-exon:        {n_multi:,} ({n_multi/n_total*100:.1f}%)")
    print()

    # Current nRNA count (from nrna.feather)
    nrna_path = Path(args.index_dir) / "nrna.feather"
    if nrna_path.exists():
        nrna_df = pd.read_feather(nrna_path)
        print(f"  Current nRNA entities: {len(nrna_df):,}")
    print()

    # Run analysis at each tolerance
    results = []
    all_examples = []
    for tol in args.tolerances:
        print(f"Analyzing tolerance = {tol} bp...")
        r = analyze_tolerance(t_df, tol)
        results.append(r)
        if r["examples"]:
            all_examples.extend(r["examples"])

    # Summary table
    print("\n" + "=" * 100)
    print("SUMMARY: nRNA Span Statistics by Tolerance Level")
    print("=" * 100)
    header = (
        f"{'Tol':>4s}  {'Exact':>8s}  {'Merged':>8s}  {'Reduct%':>8s}  "
        f"{'Annot.Eq':>8s}  {'Synth':>8s}  "
        f"{'Cl=1':>6s}  {'Cl2-5':>6s}  {'Cl6-20':>7s}  {'Cl>20':>6s}"
    )
    print(header)
    print("-" * 100)
    for r in results:
        print(
            f"{r['tolerance']:>4d}  {r['n_exact_spans']:>8,d}  {r['n_merged_spans']:>8,d}  "
            f"{r['reduction_pct']:>7.1f}%  "
            f"{r['n_annotated_equiv']:>8,d}  {r['n_synthetic_needed']:>8,d}  "
            f"{r['cluster_1']:>6,d}  {r['cluster_2_5']:>6,d}  {r['cluster_6_20']:>7,d}  "
            f"{r['cluster_20p']:>6,d}"
        )
    print()

    # Current vs proposed component counts
    print("=" * 80)
    print("COMPONENT COUNT COMPARISON")
    print("=" * 80)
    print(f"  Current architecture:")
    print(f"    mRNA components:  {n_total:>8,d} (one per transcript)")
    if nrna_path.exists():
        print(f"    nRNA components:  {len(nrna_df):>8,d} (one per unique span)")
        print(f"    Total (+ gDNA):   {n_total + len(nrna_df):>8,d} (+ 1 gDNA per locus)")
    print()
    for r in results:
        tol = r["tolerance"]
        n_synth = r["n_synthetic_needed"]
        print(f"  Proposed (tol={tol}bp):")
        print(f"    mRNA components:  {n_total:>8,d} (original transcripts)")
        print(f"    Synth nRNA:       {n_synth:>8,d} (new synthetic transcripts)")
        print(f"    Total (+ gDNA):   {n_total + n_synth:>8,d} (+ 1 gDNA per locus)")
        if nrna_path.exists():
            delta = (n_total + n_synth) - (n_total + len(nrna_df))
            sign = "+" if delta >= 0 else ""
            print(f"    Delta vs current: {sign}{delta:>8,d} components")
        print()

    # Annotated equivalent examples
    if all_examples:
        # Deduplicate
        seen = set()
        unique_examples = []
        for ex in all_examples:
            key = ex.get("single_exon_tid", "") + ex.get("merged_span", "")
            if key not in seen:
                seen.add(key)
                unique_examples.append(ex)

        print("=" * 80)
        print("ANNOTATED EQUIVALENT EXAMPLES (first 20)")
        print("=" * 80)
        for ex in unique_examples[:20]:
            match_type = ex["type"]
            gene = ex.get("gene_name", "?")
            span = ex["merged_span"]
            tid = ex.get("single_exon_tid", "?")
            extra = ""
            if match_type == "containment":
                extra = f" (SE span: {ex.get('s_span', '?')})"
            print(f"  [{match_type:>12s}] {gene:>15s}  {span}  ← {tid}{extra}")
        print()

    # Near-duplicate merging examples
    print("=" * 80)
    print("NEAR-DUPLICATE MERGING EXAMPLES (tolerance=20bp)")
    print("=" * 80)
    # Show specific cases at tol=20
    if 20 in args.tolerances:
        multi_exon = t_df[~t_df["is_single_exon"]].reset_index(drop=True)

        # Tolerance=0 (exact)
        exact_keys = multi_exon.groupby(["ref", "strand", "start", "end"]).first().reset_index()
        exact_keys["exact_span"] = exact_keys.apply(
            lambda r: (r["ref"], r["strand"], r["start"], r["end"]), axis=1
        )

        # Tolerance=20
        merged_20, t_to_merged_20 = compute_merged_spans(multi_exon, 20)

        # Find merged spans that consolidate multiple exact spans
        multi_exon_copy = multi_exon.copy()
        multi_exon_copy["merged_idx"] = t_to_merged_20
        multi_exon_copy["exact_key"] = list(
            zip(multi_exon["ref"], multi_exon["strand"], multi_exon["start"], multi_exon["end"])
        )

        # For each merged span, count distinct exact spans
        merged_groups = multi_exon_copy.groupby("merged_idx").agg(
            n_exact_spans=("exact_key", "nunique"),
            n_transcripts=("t_id", "count"),
            gene_names=("g_name", lambda x: ", ".join(sorted(set(str(v) for v in x if pd.notna(v)))[:3])),
            example_starts=("start", lambda x: sorted(set(x))[:5]),
            example_ends=("end", lambda x: sorted(set(x))[:5]),
        )
        multi_merged = merged_groups[merged_groups["n_exact_spans"] > 1].sort_values(
            "n_exact_spans", ascending=False
        )

        if len(multi_merged) > 0:
            print(f"  {len(multi_merged):,d} merged spans consolidate multiple exact spans")
            print()
            for midx, row in multi_merged.head(15).iterrows():
                mspan = merged_20.loc[midx]
                ref, strand, ms, me = mspan["ref"], mspan["strand"], mspan["start"], mspan["end"]
                strand_s = "+" if strand == 1 else "-"
                print(
                    f"  Merged: {ref}:{strand_s}:{ms}-{me}  "
                    f"({row['n_exact_spans']} exact spans, {row['n_transcripts']} tx)"
                )
                print(f"    Genes: {row['gene_names']}")
                starts = row["example_starts"]
                ends = row["example_ends"]
                print(f"    Starts: {starts}")
                print(f"    Ends:   {ends}")
                print()
        else:
            print("  No spans were consolidated at tolerance=20bp")
            print()


if __name__ == "__main__":
    main()
