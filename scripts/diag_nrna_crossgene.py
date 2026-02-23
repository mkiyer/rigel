#!/usr/bin/env python3
"""Diagnose cross-gene intronic evidence feeding phantom nRNA init.

For each region in the pristine benchmark, runs the hulkrna pipeline and
examines per-transcript intronic sense/antisense counts, nrna_init values,
and nRNA EM counts. Checks whether affected transcripts overlap with other
genes (whose exonic reads would appear as intronic evidence for this one).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from hulkrna.index import HulkIndex
from hulkrna.pipeline import run_pipeline

BENCH_DIR = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_pristine_10_regions")


def analyse_region(reg_name: str) -> dict:
    reg_dir = BENCH_DIR / reg_name
    index_dir = reg_dir / "hulkrna_index"
    cond_dir = reg_dir / "gdna_none_nrna_none_ss_1.00"
    bam_ns = cond_dir / "align" / "reads_namesort.bam"

    if not bam_ns.exists():
        print(f"  SKIP (no BAM)")
        return {}

    index = HulkIndex.load(index_dir)
    n_t = index.num_transcripts

    result = run_pipeline(
        str(bam_ns), index,
        include_multimap=True,
        em_pseudocount=0.5,
        seed=101,
    )
    counter = result.counter

    # Key arrays
    intronic_sense = counter.transcript_intronic_sense
    intronic_anti = counter.transcript_intronic_antisense
    nrna_init = counter.nrna_init
    nrna_em = counter.nrna_em_counts
    mrna_total = counter.unique_counts.sum(axis=1) + counter.em_counts.sum(axis=1)
    exonic_lengths = index.t_df["length"].values.astype(np.float64)
    transcript_spans = (
        index.t_df["end"].values - index.t_df["start"].values
    ).astype(np.float64)

    print(f"\n  Transcripts: {n_t}, Genes: {index.num_genes}")
    print(f"  Total intronic_sense: {intronic_sense.sum():.1f}")
    print(f"  Total intronic_antisense: {intronic_anti.sum():.1f}")
    print(f"  Total nrna_init: {nrna_init.sum():.1f}")
    print(f"  Total nRNA EM: {nrna_em.sum():.1f}")

    # Find transcripts with non-trivial nRNA evidence
    nrna_thresh = 1.0
    affected = np.where(nrna_em > nrna_thresh)[0]
    if len(affected) == 0:
        print(f"  No transcripts with nRNA EM > {nrna_thresh}")
        return {}

    print(f"\n  Transcripts with nRNA EM > {nrna_thresh}: {len(affected)}")

    # Build gene coordinate lookup for overlap detection
    g_starts = index.g_df["start"].values
    g_ends = index.g_df["end"].values
    g_ids = index.g_df["gene_id"].values

    total_nrna_from_overlap = 0.0
    total_nrna_no_overlap = 0.0

    for t_idx in affected:
        t_id = index.t_df.iloc[t_idx]["transcript_id"]
        g_idx = int(index.t_to_g_arr[t_idx])
        g_id = g_ids[g_idx]
        t_start = int(index.t_df.iloc[t_idx]["start"])
        t_end = int(index.t_df.iloc[t_idx]["end"])

        # Find overlapping genes (different gene)
        overlapping = []
        for other_g in range(index.num_genes):
            if other_g == g_idx:
                continue
            if g_starts[other_g] < t_end and g_ends[other_g] > t_start:
                overlap_start = max(t_start, g_starts[other_g])
                overlap_end = min(t_end, g_ends[other_g])
                overlap_bp = overlap_end - overlap_start
                overlapping.append((g_ids[other_g], overlap_bp))

        nrna_val = float(nrna_em[t_idx])
        has_overlap = len(overlapping) > 0

        if has_overlap:
            total_nrna_from_overlap += nrna_val
        else:
            total_nrna_no_overlap += nrna_val

        # Print top cases
        overlap_str = ", ".join(f"{gid}({bp}bp)" for gid, bp in overlapping) if overlapping else "NONE"
        print(f"    {t_id} ({g_id}): "
              f"nRNA_EM={nrna_val:.1f}, "
              f"init={nrna_init[t_idx]:.1f}, "
              f"int_sense={intronic_sense[t_idx]:.1f}, "
              f"int_anti={intronic_anti[t_idx]:.1f}, "
              f"span={transcript_spans[t_idx]:.0f}, exonic={exonic_lengths[t_idx]:.0f}, "
              f"overlaps=[{overlap_str}]")

    total_nrna = nrna_em.sum()
    pct_overlap = 100.0 * total_nrna_from_overlap / total_nrna if total_nrna > 0 else 0
    print(f"\n  nRNA from gene-overlap transcripts: {total_nrna_from_overlap:.1f} "
          f"({pct_overlap:.1f}% of total {total_nrna:.1f})")
    print(f"  nRNA from non-overlap transcripts: {total_nrna_no_overlap:.1f}")

    return {
        "total_nrna": float(total_nrna),
        "nrna_overlap": float(total_nrna_from_overlap),
        "nrna_no_overlap": float(total_nrna_no_overlap),
        "pct_overlap": pct_overlap,
    }


def main():
    regions = sorted(
        d.name for d in BENCH_DIR.iterdir()
        if d.is_dir() and not d.name.startswith("aggregate")
    )

    grand_total = 0.0
    grand_overlap = 0.0

    for reg in regions:
        print(f"\n{'='*60}")
        print(f"Region: {reg}")
        print(f"{'='*60}")
        try:
            info = analyse_region(reg)
            if info:
                grand_total += info["total_nrna"]
                grand_overlap += info["nrna_overlap"]
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback; traceback.print_exc()

    if grand_total > 0:
        print(f"\n{'='*60}")
        print(f"GRAND TOTAL")
        print(f"{'='*60}")
        print(f"  Total nRNA EM: {grand_total:.1f}")
        print(f"  From gene-overlap: {grand_overlap:.1f} ({100*grand_overlap/grand_total:.1f}%)")
        print(f"  From non-overlap: {grand_total - grand_overlap:.1f}")


if __name__ == "__main__":
    main()
