#!/usr/bin/env python3
"""Diagnostic: Profile actual cgranges query patterns during pipeline execution.

Instruments query_exon_with_coords to collect per-query hit statistics,
then reports the distribution and Python overhead breakdown.
"""

import sys
import time
from pathlib import Path
from collections import Counter

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from hulkrna.index import HulkIndex
from hulkrna.types import IntervalType

INDEX_DIR = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/hulkrna_index")
BAM_PATH = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/gdna_none_nrna_none_ss_1.00/align_oracle/reads_namesort.bam")


def main():
    print("Loading index...")
    index = HulkIndex.load(INDEX_DIR)

    # Instrument: collect stats on cgranges overlap calls
    cr = index.cr
    iv_t_index = index._iv_t_index
    iv_g_index = index._iv_g_index
    iv_type = index._iv_type

    # Parse BAM to get actual fragment exon blocks
    print("Parsing BAM to collect actual query patterns...")
    import pysam
    from hulkrna.bam import parse_bam_file
    from hulkrna.fragment import Fragment

    bam_fh = pysam.AlignmentFile(str(BAM_PATH), "rb")
    stats = {}

    hit_counts = []
    exon_hit_counts = []
    intron_hit_counts = []
    intergenic_hit_counts = []
    unique_t_per_query = []
    tuple_construction_count = 0
    total_queries = 0

    # Time the raw cgranges query vs Python tuple construction
    t_cgranges = 0.0
    t_tuple = 0.0
    t_total = 0.0

    max_frags = 20000  # Sample for speed
    n_frags = 0

    for nh, hits, sec_r1, sec_r2 in parse_bam_file(bam_fh, stats, include_multimap=True):
        for r1_reads, r2_reads in hits:
            frag = Fragment.from_reads(r1_reads, r2_reads)
            if frag is None or not frag.exons:
                continue
            n_frags += 1
            if n_frags > max_frags:
                break

        for exon_block in frag.exons:
            total_queries += 1

            # Time raw cgranges
            t0 = time.perf_counter()
            raw_hits = list(cr.overlap(exon_block.ref, exon_block.start, exon_block.end))
            t1 = time.perf_counter()
            t_cgranges += t1 - t0

            # Time tuple construction
            t0 = time.perf_counter()
            hits = []
            for h_start, h_end, label in raw_hits:
                hits.append((
                    int(iv_t_index[label]),
                    int(iv_g_index[label]),
                    int(iv_type[label]),
                    h_start,
                    h_end,
                ))
            t1 = time.perf_counter()
            t_tuple += t1 - t0

            tuple_construction_count += len(hits)
            hit_counts.append(len(raw_hits))

            # Count by type
            n_exon = n_intron = n_intergenic = 0
            t_set = set()
            for _, _, label in raw_hits:
                itype = int(iv_type[label])
                if itype == int(IntervalType.EXON):
                    n_exon += 1
                elif itype == int(IntervalType.INTRON):
                    n_intron += 1
                else:
                    n_intergenic += 1
                t_idx = int(iv_t_index[label])
                if t_idx >= 0:
                    t_set.add(t_idx)
            exon_hit_counts.append(n_exon)
            intron_hit_counts.append(n_intron)
            intergenic_hit_counts.append(n_intergenic)
            unique_t_per_query.append(len(t_set))

    print(f"\nProcessed {n_frags} fragments, {total_queries} exon block queries")

    hit_arr = np.array(hit_counts)
    exon_arr = np.array(exon_hit_counts)
    intron_arr = np.array(intron_hit_counts)
    ig_arr = np.array(intergenic_hit_counts)
    ut_arr = np.array(unique_t_per_query)

    print(f"\n{'='*70}")
    print(f"QUERY HIT STATISTICS (from actual BAM fragments)")
    print(f"{'='*70}")

    print(f"\n--- Total hits per query ---")
    print(f"  median={np.median(hit_arr):.0f}  mean={hit_arr.mean():.1f}  "
          f"min={hit_arr.min()}  max={hit_arr.max()}  total={hit_arr.sum():,d}")
    print(f"  p90={np.percentile(hit_arr, 90):.0f}  "
          f"p95={np.percentile(hit_arr, 95):.0f}  "
          f"p99={np.percentile(hit_arr, 99):.0f}")

    print(f"\n--- Exon hits per query ---")
    print(f"  median={np.median(exon_arr):.0f}  mean={exon_arr.mean():.1f}  "
          f"max={exon_arr.max()}")

    print(f"\n--- Intron hits per query ---")
    print(f"  median={np.median(intron_arr):.0f}  mean={intron_arr.mean():.1f}  "
          f"max={intron_arr.max()}")

    print(f"\n--- Intergenic hits per query ---")
    print(f"  median={np.median(ig_arr):.0f}  mean={ig_arr.mean():.1f}  "
          f"max={ig_arr.max()}")

    print(f"\n--- Unique transcripts per query ---")
    print(f"  median={np.median(ut_arr):.0f}  mean={ut_arr.mean():.1f}  "
          f"max={ut_arr.max()}")

    print(f"\n--- Timing breakdown (extrapolated to 100k frags) ---")
    scale = 100000 / n_frags
    print(f"  Raw cgranges queries:   {t_cgranges*scale:.2f}s")
    print(f"  Tuple construction:     {t_tuple*scale:.2f}s")
    print(f"  Total tuples built:     {int(tuple_construction_count*scale):,d}")

    # --- Collapsed-index simulation ---
    print(f"\n--- Collapsed-index simulation ---")
    # If we collapsed intervals so identical (ref,start,end) return one label
    # with a transcript set, how many fewer hits per query?
    collapsed_hit_counts = []
    for exon_block_idx in range(min(total_queries, len(hit_counts))):
        # We can't replay the exact query, but we can estimate
        # from the hit stats: unique boundaries vs total hits
        pass

    # Instead, re-run with dedup
    collapsed_hits_total = 0
    bam_fh2 = pysam.AlignmentFile(str(BAM_PATH), "rb")
    stats2 = {}
    n2 = 0
    for nh, hits, sec_r1, sec_r2 in parse_bam_file(bam_fh2, stats2, include_multimap=True):
        for r1_reads, r2_reads in hits:
            frag = Fragment.from_reads(r1_reads, r2_reads)
            if frag is None or not frag.exons:
                continue
            n2 += 1
            if n2 > max_frags:
                break
        for exon_block in frag.exons:
            raw_hits = list(cr.overlap(exon_block.ref, exon_block.start, exon_block.end))
            # How many unique (start, end) boundaries?
            boundaries = set()
            for h_start, h_end, label in raw_hits:
                boundaries.add((h_start, h_end))
            collapsed_hits_total += len(boundaries)
            collapsed_hit_counts.append(len(boundaries))

    c_arr = np.array(collapsed_hit_counts)
    print(f"  Current avg hits/query:    {hit_arr.mean():.1f}")
    print(f"  Collapsed avg hits/query:  {c_arr.mean():.1f}")
    print(f"  Hit reduction:             {100*(1 - c_arr.mean()/hit_arr.mean()):.1f}%")
    print(f"  Collapsed total hits:      {c_arr.sum():,d} (was {hit_arr.sum():,d})")


if __name__ == "__main__":
    main()
