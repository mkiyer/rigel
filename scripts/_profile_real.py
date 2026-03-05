#!/usr/bin/env python3
"""Profile hulkrna on the real dataset.

Usage:
    cd /Users/mkiyer/proj/hulkrna
    PYTHONPATH=src python scripts/_profile_real.py
"""
import cProfile
import pstats
import time
import sys
import io

BAM = "/Users/mkiyer/Downloads/hulkrna_runs/star.collate.markdup.bam"
INDEX = "/Users/mkiyer/Downloads/hulkrna_runs/hulkrna_index"


def main():
    from hulkrna.index import TranscriptIndex
    from hulkrna.pipeline import run_pipeline

    # --- Load index ---
    t0 = time.perf_counter()
    index = TranscriptIndex.load(INDEX)
    t_index = time.perf_counter() - t0
    print(f"Index loaded in {t_index:.2f}s: "
          f"{index.num_transcripts} transcripts, {index.num_genes} genes")

    # --- Profile pipeline ---
    print(f"\n{'='*70}")
    print("Profiling hulkrna pipeline on real BAM...")
    print(f"{'='*70}\n")

    t0 = time.perf_counter()
    profiler = cProfile.Profile()
    profiler.enable()
    result = run_pipeline(BAM, index)
    profiler.disable()
    wall = time.perf_counter() - t0

    print(f"\n{'='*70}")
    print(f"Pipeline wall time: {wall:.2f}s")
    print(f"{'='*70}\n")

    # --- Stats ---
    stats = result.stats
    print(f"Fragments processed: {stats.total:,}")
    frags_per_sec = stats.total / wall if wall > 0 else 0
    print(f"Throughput: {frags_per_sec:,.0f} fragments/sec")
    print()

    # --- cProfile output ---
    ps = pstats.Stats(profiler, stream=sys.stdout)

    ps.sort_stats("cumulative")
    print("=== Top 60 by cumulative time ===")
    ps.print_stats(60)

    ps.sort_stats("tottime")
    print("\n=== Top 60 by self time ===")
    ps.print_stats(60)

    # --- Function call count ---
    total_calls = sum(v[0] for v in ps.stats.values())
    print(f"\nTotal function calls: {total_calls:,}")

    # Save profile
    prof_path = "/tmp/hulkrna_real_profile.prof"
    profiler.dump_stats(prof_path)
    print(f"\nProfile saved to {prof_path}")


if __name__ == "__main__":
    main()
