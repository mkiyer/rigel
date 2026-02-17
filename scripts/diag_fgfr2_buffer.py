#!/usr/bin/env python3
from collections import Counter
from pathlib import Path

import pysam

from hulkrna.categories import CountCategory
from hulkrna.index import HulkIndex
from hulkrna.pipeline import scan_and_buffer
from hulkrna.types import MergeCriteria


def main() -> int:
    base = Path(
        "/Users/mkiyer/Downloads/hulkrna_runs/fgfr2_depth_runs/seed_303/chr10_121476640_121601584"
    )
    index = HulkIndex.load(base / "hulkrna_index")
    bam = base / "align" / "reads.bam"

    with pysam.AlignmentFile(str(bam), "rb") as bam_iter:
        stats, _, _, buf = scan_and_buffer(
            bam_iter,
            index,
            skip_duplicates=True,
            include_multimap=False,
            sj_strand_tag="auto",
            log_every=10_000_000,
        )

    cat_counts = Counter()
    merge_counts = Counter()
    spliced_annot_merge = Counter()
    size_by_cat = Counter()

    for bf in buf:
        cat = int(bf.count_cat)
        mc = int(bf.merge_criteria)
        n = len(bf.t_inds)
        cat_counts[cat] += 1
        merge_counts[mc] += 1
        spliced_annot_merge[mc] += int(cat == int(CountCategory.SPLICED_ANNOT))
        size_by_cat[(cat, min(n, 10))] += 1

    print(f"stats.n_fragments={stats.n_fragments}")
    print(f"stats.n_with_annotated_sj={stats.n_with_annotated_sj}")
    print(f"stats.n_with_unannotated_sj={stats.n_with_unannotated_sj}")
    print(f"stats.n_with_exon={stats.n_with_exon}")
    print(f"stats.n_with_intron_fallback={stats.n_with_intron_fallback}")
    print(f"stats.n_unique_gene={stats.n_unique_gene}")
    print(f"stats.n_multi_gene={stats.n_multi_gene}")

    print("\nCountCategory counts:")
    for c in sorted(cat_counts):
        print(f"{CountCategory(c).name}: {cat_counts[c]}")

    print("\nMergeCriteria counts:")
    for m in sorted(merge_counts):
        print(f"{MergeCriteria(m).name}: {merge_counts[m]}")

    print("\nSPLICED_ANNOT merge criteria:")
    for m in sorted(spliced_annot_merge):
        if spliced_annot_merge[m] > 0:
            print(f"{MergeCriteria(m).name}: {spliced_annot_merge[m]}")

    print("\nCandidate-size buckets by category (bucket=min(len(t_inds),10)):")
    for (c, b), v in sorted(size_by_cat.items()):
        if b in (1, 2, 3, 4, 5, 10):
            print(f"{CountCategory(c).name} n={b}: {v}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
