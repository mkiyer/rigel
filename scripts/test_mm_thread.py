#!/usr/bin/env python3
"""Quick test: compare 1-thread vs 4-thread scan with include_multimap=True."""
import sys
import time

from rigel._bam_impl import BamScanner
from rigel.index import TranscriptIndex


def main():
    idx_path = sys.argv[1] if len(sys.argv) > 1 else "/Users/mkiyer/Downloads/rigel_runs/rigel_index"
    bam_path = sys.argv[2] if len(sys.argv) > 2 else "/Users/mkiyer/Downloads/rigel_runs/sim_ccle_hela_cervix/gdna_high_ss_0.90_nrna_low/sim_oracle.bam"
    include_mm = True

    idx = TranscriptIndex.load(idx_path)

    rc = idx._resolver
    rc.set_gene_strands(idx.g_to_strand_arr.tolist())
    rc.set_transcript_strands(idx.t_to_strand_arr.tolist())

    # Single-threaded
    s1 = BamScanner(rc, "none", True, include_mm)
    t0 = time.time()
    r1 = s1.scan(bam_path)
    dt1 = time.time() - t0
    st1 = r1["stats"]
    print(f"1-thread: acc={r1['accumulator_size']}, n_frags={st1['n_fragments']}, "
          f"unique={st1['unique']}, multi={st1['multimapping']}, "
          f"mm_groups={st1['n_multimapper_groups']}, mm_alns={st1['n_multimapper_alignments']}, "
          f"time={dt1:.1f}s")

    # Multi-threaded
    s2 = BamScanner(rc, "none", True, include_mm)
    t0 = time.time()
    r2 = s2.scan(bam_path, n_workers=4)
    dt2 = time.time() - t0
    st2 = r2["stats"]
    print(f"4-thread: acc={r2['accumulator_size']}, n_frags={st2['n_fragments']}, "
          f"unique={st2['unique']}, multi={st2['multimapping']}, "
          f"mm_groups={st2['n_multimapper_groups']}, mm_alns={st2['n_multimapper_alignments']}, "
          f"time={dt2:.1f}s")

    # Print comparison
    keys = [
        "n_read_names", "unique", "multimapping", "n_fragments",
        "n_chimeric", "n_intergenic_unspliced", "n_intergenic_spliced",
        "n_with_exon", "n_with_annotated_sj", "n_with_unannotated_sj",
        "n_same_strand", "n_ambig_strand",
        "n_multimapper_groups", "n_multimapper_alignments",
    ]
    print("\nStat comparison:")
    any_diff = False
    for k in keys:
        a, b = st1.get(k, -1), st2.get(k, -1)
        tag = "" if a == b else " *** MISMATCH ***"
        if a != b:
            any_diff = True
        print(f"  {k:40s}: {a:>12} vs {b:>12}{tag}")

    acc_match = r1["accumulator_size"] == r2["accumulator_size"]
    print(f"\nACC MATCH: {acc_match}")
    print(f"ANY STAT DIFF: {any_diff}")

    return 0 if acc_match and not any_diff else 1


if __name__ == "__main__":
    sys.exit(main())
