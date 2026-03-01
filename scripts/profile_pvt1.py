#!/usr/bin/env python3
"""Profile hulkrna pipeline on PVT1 oracle BAM with cProfile.

Usage:
    PYTHONPATH=src python scripts/profile_pvt1.py [--gdna-nrna]
"""
import argparse
import cProfile
import pstats
import time
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gdna-nrna", action="store_true",
                        help="Use gDNA 10%% + nRNA 10%% condition if available")
    args = parser.parse_args()

    # Paths
    base = "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC"
    index_dir = f"{base}/hulkrna_index"
    cond = "gdna_none_nrna_none_ss_1.00"
    bam_path = f"{base}/{cond}/align_oracle/reads_namesort.bam"

    from hulkrna.index import HulkIndex
    from hulkrna.pipeline import run_pipeline

    print(f"BAM: {bam_path}")
    print(f"Index: {index_dir}")
    print()

    # Load index (not profiled — one-time cost)
    index = HulkIndex.load(index_dir)
    print(f"Index loaded: {index.num_transcripts} transcripts, {index.num_genes} genes")
    print()

    # --- Wall-clock timing ---
    t0 = time.perf_counter()
    profiler = cProfile.Profile()
    profiler.enable()
    result = run_pipeline(bam_path, index)
    profiler.disable()
    wall = time.perf_counter() - t0

    print(f"\n{'='*70}")
    print(f"Pipeline wall time: {wall:.2f}s")
    print(f"{'='*70}\n")

    # --- Stage timings (from stats) ---
    stats = result.stats
    print(f"Fragments processed: {stats.total}")
    print()

    # --- cProfile output ---
    ps = pstats.Stats(profiler)
    ps.sort_stats("cumulative")
    print("=== Top 40 by cumulative time ===")
    ps.print_stats(40)

    ps.sort_stats("tottime")
    print("\n=== Top 40 by self time ===")
    ps.print_stats(40)

    # --- Function call count ---
    total_calls = sum(v[0] for v in ps.stats.values())
    print(f"\nTotal function calls: {total_calls:,}")

    # Save to file for later analysis
    profiler.dump_stats("/tmp/hulkrna_pvt1_profile.prof")
    print("\nProfile saved to /tmp/hulkrna_pvt1_profile.prof")


if __name__ == "__main__":
    main()
