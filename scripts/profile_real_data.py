#!/usr/bin/env python3
"""Profile hulkrna on a real STAR-aligned BAM file.

Usage:
    cd /Users/mkiyer/proj/hulkrna
    conda run -n hulkrna python scripts/profile_real_data.py
"""
import cProfile
import pstats
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BAM_PATH = Path(
    "/Users/mkiyer/Downloads/hulkrna_runs/"
    "mctp_LBX0069_SI_42153_HFFFMDRX7/bam/star.collate.markdup.bam"
)
INDEX_DIR = Path(
    "/Users/mkiyer/Downloads/hulkrna_runs/refs/human/hulkrna_index"
)
PROF_OUT = Path("/tmp/hulkrna_real_profile.prof")


def main():
    from hulkrna.index import TranscriptIndex
    from hulkrna.pipeline import run_pipeline

    # --- Validate paths ---
    if not BAM_PATH.exists():
        sys.exit(f"BAM not found: {BAM_PATH}")
    if not INDEX_DIR.is_dir():
        sys.exit(f"Index not found: {INDEX_DIR}")

    # --- Load index ---
    print(f"Loading index from {INDEX_DIR} ...")
    t0 = time.perf_counter()
    index = TranscriptIndex.load(str(INDEX_DIR))
    t_index = time.perf_counter() - t0
    print(f"  {index.num_transcripts:,} transcripts, "
          f"{index.num_genes:,} genes  ({t_index:.1f}s)")

    # --- Profile pipeline ---
    print(f"\n{'='*70}")
    print(f"Profiling hulkrna pipeline on real data")
    print(f"  BAM: {BAM_PATH.name}")
    print(f"{'='*70}\n")

    t0 = time.perf_counter()
    profiler = cProfile.Profile()
    profiler.enable()
    result = run_pipeline(str(BAM_PATH), index)
    profiler.disable()
    wall = time.perf_counter() - t0

    # --- Summary ---
    print(f"\n{'='*70}")
    print(f"Pipeline wall time: {wall:.2f}s")
    print(f"{'='*70}\n")

    stats = result.stats
    sd = stats.to_dict()
    for key in sorted(sd):
        val = sd[key]
        if isinstance(val, (int, float)):
            print(f"  {key}: {val:,.2f}" if isinstance(val, float)
                  else f"  {key}: {val:,}")

    total_frags = stats.total
    frags_per_sec = total_frags / wall if wall > 0 else 0
    print(f"\nThroughput: {frags_per_sec:,.0f} fragments/sec")

    # --- cProfile output ---
    ps = pstats.Stats(profiler)

    ps.sort_stats("cumulative")
    print(f"\n{'='*70}")
    print("Top 60 by cumulative time")
    print(f"{'='*70}")
    ps.print_stats(60)

    ps.sort_stats("tottime")
    print(f"\n{'='*70}")
    print("Top 60 by self time")
    print(f"{'='*70}")
    ps.print_stats(60)

    # --- Function call count ---
    total_calls = sum(v[0] for v in ps.stats.values())
    print(f"\nTotal function calls: {total_calls:,}")

    # Save profile
    profiler.dump_stats(str(PROF_OUT))
    print(f"Profile saved to {PROF_OUT}")


if __name__ == "__main__":
    main()
