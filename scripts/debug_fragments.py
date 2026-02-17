#!/usr/bin/env python3
"""
Fragment-level debugging script for hulkrna pipeline.

Traces each simulated fragment through the entire pipeline:
  1. Simulation truth (FASTQ read names → true transcript)
  2. Alignment (BAM records, mapping quality, splice junctions)
  3. Fragment resolution (candidate transcript set, merge criteria)
  4. Equivalence classes (groups of fragments with same candidate sets)

Reports discrepancies between truth and resolution at every stage.

Usage:
  PYTHONPATH=src python scripts/debug_fragments.py \
      --region FGFR2 \
      --bench-dir /Users/mkiyer/Downloads/hulkrna_runs/bench_zero_gdna_v2
"""
from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hulkrna.bam import parse_bam_file
from hulkrna.categories import CountCategory
from hulkrna.fragment import Fragment
from hulkrna.index import HulkIndex
from hulkrna.resolution import resolve_fragment
from hulkrna.types import MergeCriteria, Strand

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class FragmentTrace:
    """Complete trace of one simulated fragment through the pipeline."""
    read_name: str
    truth_tid: str  # True transcript from simulation
    truth_pos: str  # Simulated fragment coordinates

    # Alignment stage
    is_mapped: bool = True
    is_primary: bool = True
    is_secondary: bool = False
    is_supplementary: bool = False
    mapq: int = 0
    n_cigar_ops: int = 0
    has_splice: bool = False  # BAM record has N (splice) cigar op
    alignment_ref: str = ""
    alignment_start: int = 0
    alignment_end: int = 0

    # Resolution stage
    resolved: bool = False
    candidate_tids: frozenset[str] = field(default_factory=frozenset)
    n_candidates: int = 0
    truth_in_candidates: bool = False
    count_cat: str = ""
    merge_criteria: str = ""
    is_chimeric: bool = False
    insert_size: int = -1
    n_genes: int = 0
    exon_strand: str = ""

    # Classification
    classification: str = ""  # correct_unique, correct_ambig, wrong, lost


def parse_truth_from_read_name(qname: str) -> tuple[str, str]:
    """Extract truth transcript ID and position from FASTQ read name.

    Format: TRANSCRIPT_ID:START-END:STRAND:INDEX
    """
    # Remove /1 or /2 suffix
    if qname.endswith("/1") or qname.endswith("/2"):
        qname = qname[:-2]

    parts = qname.split(":")
    tid = parts[0]
    pos_info = ":".join(parts[1:]) if len(parts) > 1 else ""
    return tid, pos_info


def load_truth_counts(per_tx_csv: Path) -> dict[str, int]:
    """Load truth fragment counts from per_transcript_counts.csv."""
    truth = {}
    with open(per_tx_csv) as f:
        for row in csv.DictReader(f):
            truth[row["transcript_id"]] = int(row["truth"])
    return truth


def load_index_transcript_ids(index: HulkIndex) -> dict[int, str]:
    """Map transcript index → transcript ID."""
    return dict(enumerate(index.t_df["t_id"].tolist()))


def load_index_tid_to_idx(index: HulkIndex) -> dict[str, int]:
    """Map transcript ID → index."""
    return {tid: i for i, tid in enumerate(index.t_df["t_id"].tolist())}


def candidate_tids_from_indices(
    t_inds: frozenset[int], t_idx_to_tid: dict[int, str]
) -> frozenset[str]:
    """Convert transcript indices to transcript IDs."""
    return frozenset(t_idx_to_tid.get(i, f"?{i}") for i in t_inds)


# ---------------------------------------------------------------------------
# Main tracing
# ---------------------------------------------------------------------------

def trace_fragments(
    bam_path: Path,
    index: HulkIndex,
    sj_strand_tag: str = "ts",
) -> list[FragmentTrace]:
    """Trace all fragments through alignment → resolution."""

    t_idx_to_tid = load_index_transcript_ids(index)
    traces: list[FragmentTrace] = []

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        bam_stats: dict = {
            "n_reads": 0,
            "n_paired": 0,
            "n_unmapped": 0,
            "n_duplicate": 0,
            "n_supplementary": 0,
            "n_secondary": 0,
            "n_read_names": 0,
            "n_singletons": 0,
            "n_proper_pairs": 0,
        }

        for nh, hits in parse_bam_file(
            bam, bam_stats,
            skip_duplicates=True,
            include_multimap=True,
        ):
            num_hits = max(nh, len(hits))

            for r1_reads, r2_reads in hits:
                # Get truth from read name
                qname = ""
                if r1_reads:
                    qname = r1_reads[0].query_name
                elif r2_reads:
                    qname = r2_reads[0].query_name

                truth_tid, truth_pos = parse_truth_from_read_name(qname)

                trace = FragmentTrace(
                    read_name=qname,
                    truth_tid=truth_tid,
                    truth_pos=truth_pos,
                )

                # Alignment info from primary read
                primary = r1_reads[0] if r1_reads else (r2_reads[0] if r2_reads else None)
                if primary:
                    trace.is_mapped = not primary.is_unmapped
                    trace.mapq = primary.mapping_quality
                    trace.alignment_ref = primary.reference_name or ""
                    trace.alignment_start = primary.reference_start
                    trace.alignment_end = primary.reference_end or 0
                    trace.has_splice = any(
                        op == 3 for op, _ in (primary.cigartuples or [])
                    )

                # Build fragment and resolve
                frag = Fragment.from_reads(
                    r1_reads, r2_reads, sj_strand_tag=sj_strand_tag
                )

                if not frag.exons:
                    trace.resolved = False
                    trace.classification = "no_exons"
                    traces.append(trace)
                    continue

                result = resolve_fragment(frag, index)

                if result is None:
                    trace.resolved = False
                    trace.classification = "intergenic"
                    traces.append(trace)
                    continue

                trace.resolved = True

                if result.is_chimeric:
                    trace.is_chimeric = True
                    trace.classification = "chimeric"
                    traces.append(trace)
                    continue

                # Resolution details
                trace.candidate_tids = candidate_tids_from_indices(
                    result.t_inds, t_idx_to_tid
                )
                trace.n_candidates = len(result.t_inds)
                trace.truth_in_candidates = (truth_tid in trace.candidate_tids)
                trace.count_cat = CountCategory(result.count_cat).name
                trace.merge_criteria = MergeCriteria(result.merge_criteria).name
                trace.insert_size = result.insert_size
                trace.n_genes = result.n_genes
                trace.exon_strand = Strand(result.exon_strand).name if result.exon_strand else "NONE"

                # Classify
                if trace.truth_in_candidates:
                    if trace.n_candidates == 1:
                        trace.classification = "correct_unique"
                    else:
                        trace.classification = "correct_ambig"
                else:
                    trace.classification = "wrong_candidates"

                traces.append(trace)

    return traces


# ---------------------------------------------------------------------------
# Analysis & reporting
# ---------------------------------------------------------------------------

def analyze_traces(
    traces: list[FragmentTrace],
    truth_counts: dict[str, int],
    index: HulkIndex,
) -> None:
    """Comprehensive analysis of fragment traces."""

    t_idx_to_tid = load_index_transcript_ids(index)
    tid_to_idx = load_index_tid_to_idx(index)
    all_index_tids = set(t_idx_to_tid.values())
    truth_tids = set(truth_counts.keys())

    # =====================================================================
    # Section 1: Overview
    # =====================================================================
    print("\n" + "=" * 100)
    print("FRAGMENT-LEVEL DEBUG REPORT")
    print("=" * 100)

    n_total = len(traces)
    classifications = Counter(t.classification for t in traces)

    print(f"\nTotal fragments: {n_total}")
    print(f"\nClassification breakdown:")
    for cls, n in sorted(classifications.items(), key=lambda x: -x[1]):
        pct = 100.0 * n / n_total if n_total else 0
        print(f"  {cls:<25} {n:>8,} ({pct:>6.2f}%)")

    # =====================================================================
    # Section 2: Transcript set comparison (Index vs Truth)
    # =====================================================================
    print(f"\n{'=' * 100}")
    print("INDEX vs TRUTH TRANSCRIPT SETS")
    print(f"{'=' * 100}")

    in_both = truth_tids & all_index_tids
    in_truth_only = truth_tids - all_index_tids
    in_index_only = all_index_tids - truth_tids

    print(f"  Index transcripts:     {len(all_index_tids):>5}")
    print(f"  Truth transcripts:     {len(truth_tids):>5}")
    print(f"  In both:               {len(in_both):>5}")
    print(f"  In truth only (MISSING from index): {len(in_truth_only):>5}")
    print(f"  In index only (PHANTOM - not expressed): {len(in_index_only):>5}")

    if in_truth_only:
        print(f"\n  MISSING transcripts (in truth but not in index):")
        for tid in sorted(in_truth_only):
            print(f"    {tid} (truth={truth_counts.get(tid, 0)})")

    if in_index_only:
        print(f"\n  PHANTOM transcripts (in index but not expressed):")
        for tid in sorted(in_index_only):
            print(f"    {tid}")

    # =====================================================================
    # Section 3: Resolution accuracy per truth transcript
    # =====================================================================
    print(f"\n{'=' * 100}")
    print("PER-TRANSCRIPT RESOLUTION ACCURACY")
    print(f"{'=' * 100}")

    # Group traces by truth transcript
    by_truth: dict[str, list[FragmentTrace]] = defaultdict(list)
    for t in traces:
        by_truth[t.truth_tid].append(t)

    print(f"\n{'transcript_id':<28} {'truth':>6} {'resolved':>8} {'correct':>8} "
          f"{'unique':>7} {'ambig':>7} {'wrong':>7} {'lost':>6} "
          f"{'avg_cand':>9} {'cat_dist':<30}")
    print("-" * 140)

    rows = []
    for tid in sorted(truth_counts.keys(), key=lambda t: -truth_counts[t]):
        frags = by_truth.get(tid, [])
        n = truth_counts[tid]
        n_resolved = sum(1 for f in frags if f.resolved)
        n_correct = sum(1 for f in frags if f.truth_in_candidates)
        n_unique = sum(1 for f in frags if f.classification == "correct_unique")
        n_ambig = sum(1 for f in frags if f.classification == "correct_ambig")
        n_wrong = sum(1 for f in frags if f.classification == "wrong_candidates")
        n_lost = n - len(frags)  # Not in BAM at all

        avg_cand = np.mean([f.n_candidates for f in frags if f.resolved]) if any(f.resolved for f in frags) else 0

        # Count category distribution
        cat_counts = Counter(f.count_cat for f in frags if f.resolved)
        cat_str = ", ".join(f"{c}={v}" for c, v in sorted(cat_counts.items()))

        print(f"{tid:<28} {n:>6} {n_resolved:>8} {n_correct:>8} "
              f"{n_unique:>7} {n_ambig:>7} {n_wrong:>7} {n_lost:>6} "
              f"{avg_cand:>9.1f} {cat_str:<30}")

        rows.append({
            "tid": tid, "truth": n, "resolved": n_resolved,
            "correct": n_correct, "unique": n_unique, "ambig": n_ambig,
            "wrong": n_wrong, "lost": n_lost, "avg_cand": avg_cand,
        })

    # =====================================================================
    # Section 4: Equivalence classes
    # =====================================================================
    print(f"\n{'=' * 100}")
    print("EQUIVALENCE CLASSES (top 30 by fragment count)")
    print(f"{'=' * 100}")

    # Group fragments by candidate set
    eq_classes: dict[frozenset[str], list[FragmentTrace]] = defaultdict(list)
    for t in traces:
        if t.resolved and not t.is_chimeric and t.candidate_tids:
            eq_classes[t.candidate_tids].append(t)

    sorted_eqs = sorted(eq_classes.items(), key=lambda x: -len(x[1]))

    print(f"\nTotal equivalence classes: {len(sorted_eqs)}")
    print(f"\n{'#':>4} {'n_frags':>8} {'n_txns':>7} {'truth_dist':<50} {'transcripts':<60}")
    print("-" * 140)

    for i, (cand_set, frags) in enumerate(sorted_eqs[:30]):
        truth_dist = Counter(f.truth_tid for f in frags)
        truth_str = ", ".join(f"{t}={c}" for t, c in truth_dist.most_common(5))
        if len(truth_dist) > 5:
            truth_str += f" +{len(truth_dist)-5} more"

        cand_list = sorted(cand_set)
        cand_str = ", ".join(cand_list[:4])
        if len(cand_list) > 4:
            cand_str += f" +{len(cand_list)-4} more"

        print(f"{i+1:>4} {len(frags):>8} {len(cand_set):>7} {truth_str:<50} {cand_str:<60}")

    # =====================================================================
    # Section 5: Wrong candidate analysis
    # =====================================================================
    print(f"\n{'=' * 100}")
    print("WRONG CANDIDATE ANALYSIS (truth transcript NOT in candidate set)")
    print(f"{'=' * 100}")

    wrong_frags = [t for t in traces if t.classification == "wrong_candidates"]
    if wrong_frags:
        # Group by truth transcript
        wrong_by_truth = defaultdict(list)
        for f in wrong_frags:
            wrong_by_truth[f.truth_tid].append(f)

        print(f"\nTotal wrong: {len(wrong_frags)} fragments from {len(wrong_by_truth)} transcripts")
        print(f"\n{'truth_tid':<28} {'n_wrong':>8} {'truth':>6} {'candidate_sets (top 3)':<80}")
        print("-" * 130)

        for tid in sorted(wrong_by_truth.keys(), key=lambda t: -len(wrong_by_truth[t])):
            frags = wrong_by_truth[tid]
            # Most common wrong candidate sets
            cand_counter = Counter(f.candidate_tids for f in frags)
            top_sets = cand_counter.most_common(3)
            sets_str = " | ".join(
                f"[{','.join(sorted(s)[:3])}{'...' if len(s) > 3 else ''}]x{c}"
                for s, c in top_sets
            )
            print(f"{tid:<28} {len(frags):>8} {truth_counts.get(tid, 0):>6} {sets_str:<80}")

    # =====================================================================
    # Section 6: Phantom transcript absorption
    # =====================================================================
    print(f"\n{'=' * 100}")
    print("PHANTOM TRANSCRIPT ABSORPTION")
    print(f"  (How many fragments have phantom transcripts in their candidate sets)")
    print(f"{'=' * 100}")

    phantom_tids = in_index_only
    if phantom_tids:
        n_with_phantom = 0
        n_only_phantom = 0  # All candidates are phantoms
        phantom_appearance = Counter()  # How often each phantom appears

        for t in traces:
            if not t.resolved or t.is_chimeric:
                continue
            phantoms_in_cands = t.candidate_tids & phantom_tids
            if phantoms_in_cands:
                n_with_phantom += 1
                phantom_appearance.update(phantoms_in_cands)
                real_cands = t.candidate_tids - phantom_tids
                if not real_cands:
                    n_only_phantom += 1

        n_resolved = sum(1 for t in traces if t.resolved and not t.is_chimeric)
        print(f"\n  Fragments with ≥1 phantom candidate: {n_with_phantom:>8} / {n_resolved} "
              f"({100*n_with_phantom/n_resolved:.1f}%)")
        print(f"  Fragments with ONLY phantom candidates: {n_only_phantom:>8}")

        print(f"\n  Phantom transcript frequency in candidate sets:")
        print(f"  {'phantom_tid':<28} {'appearances':>12} {'pct_of_resolved':>16}")
        print(f"  {'-'*56}")
        for ptid, count in phantom_appearance.most_common(20):
            print(f"  {ptid:<28} {count:>12} {100*count/n_resolved:>15.1f}%")

    # =====================================================================
    # Section 7: Splice junction analysis
    # =====================================================================
    print(f"\n{'=' * 100}")
    print("COUNT CATEGORY DISTRIBUTION")
    print(f"{'=' * 100}")

    resolved_traces = [t for t in traces if t.resolved and not t.is_chimeric]
    cat_dist = Counter(t.count_cat for t in resolved_traces)
    for cat, n in sorted(cat_dist.items(), key=lambda x: -x[1]):
        pct = 100.0 * n / len(resolved_traces) if resolved_traces else 0
        print(f"  {cat:<25} {n:>8} ({pct:>6.2f}%)")

    # By merge criteria
    merge_dist = Counter(t.merge_criteria for t in resolved_traces)
    print(f"\nMerge criteria distribution:")
    for mc, n in sorted(merge_dist.items(), key=lambda x: -x[1]):
        pct = 100.0 * n / len(resolved_traces) if resolved_traces else 0
        print(f"  {mc:<30} {n:>8} ({pct:>6.2f}%)")

    # =====================================================================
    # Section 8: Detailed wrong fragment examples
    # =====================================================================
    if wrong_frags:
        print(f"\n{'=' * 100}")
        print("DETAILED WRONG FRAGMENT EXAMPLES (first 20)")
        print(f"{'=' * 100}")

        for i, f in enumerate(wrong_frags[:20]):
            print(f"\n  Fragment #{i+1}: {f.read_name}")
            print(f"    Truth: {f.truth_tid} (pos={f.truth_pos})")
            print(f"    Alignment: ref={f.alignment_ref}, start={f.alignment_start}, "
                  f"end={f.alignment_end}, mapq={f.mapq}, splice={f.has_splice}")
            print(f"    Resolution: cat={f.count_cat}, merge={f.merge_criteria}, "
                  f"n_genes={f.n_genes}, insert_size={f.insert_size}")
            print(f"    Candidates ({f.n_candidates}): {sorted(f.candidate_tids)}")
            print(f"    Strand: {f.exon_strand}")


def main():
    parser = argparse.ArgumentParser(description="Fragment-level pipeline debugger")
    parser.add_argument("--region", required=True, help="Region name")
    parser.add_argument("--bench-dir", required=True, type=Path,
                        help="Benchmark output directory")
    parser.add_argument("--sj-strand-tag", default="ts",
                        help="Splice junction strand tag (default: ts for minimap2)")
    parser.add_argument("--max-frags", type=int, default=0,
                        help="Max fragments to trace (0=all)")
    args = parser.parse_args()

    reg_dir = args.bench_dir / args.region
    if not reg_dir.exists():
        logger.error(f"Region directory not found: {reg_dir}")
        sys.exit(1)

    bam_path = reg_dir / "align" / "reads.bam"
    index_dir = reg_dir / "hulkrna_index"
    per_tx_csv = reg_dir / "per_transcript_counts.csv"

    logger.info(f"Loading index from {index_dir}")
    index = HulkIndex.load(index_dir)
    truth_counts = load_truth_counts(per_tx_csv)

    logger.info(f"Tracing fragments from {bam_path}")
    traces = trace_fragments(bam_path, index, sj_strand_tag=args.sj_strand_tag)

    logger.info(f"Traced {len(traces)} fragments, analyzing...")
    analyze_traces(traces, truth_counts, index)


if __name__ == "__main__":
    main()
