#!/usr/bin/env python3
"""Compare gap correction effectiveness before and after the overlap fix.

Uses the same minimap2 BAM but loads the current (post-fix) code to trace
how many gaps are now being corrected.
"""
import sys
sys.path.insert(0, "src")

import pysam
import numpy as np
from collections import Counter, defaultdict
from rigel.index import TranscriptIndex

BASE = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v2"
COND = "gdna_none_ss_0.95_nrna_none"
BAM = f"{BASE}/{COND}/align_minimap2/reads_namesort.bam"
INDEX = f"{BASE}/rigel_index"

# Load index 
idx = TranscriptIndex.load(INDEX)

# Get the SJ data for manual inspection
import pyarrow.feather as pf
sj_df = pf.read_table(f"{INDEX}/sj.feather").to_pandas()
print(f"SJ index: {len(sj_df)} introns")

# Build a lookup: (chrom, start, end) -> list of t_indices
sj_lookup = {}
for _, row in sj_df.iterrows():
    key = (row['ref'], row['start'], row['end'])
    if key not in sj_lookup:
        sj_lookup[key] = []
    sj_lookup[key].append(row['t_index'])

print(f"Unique SJ coordinates: {len(sj_lookup)}")

# Now process BAM fragments and compute FL both ways (overlap vs containment)
# Read pairs and compute FL
n_pairs = 0
n_sampled = 0
n_max = 200000

# Counters for FL comparison
fl_old = []  # containment-based FL
fl_new = []  # overlap-based FL (current code)
fl_both_valid = 0
fl_only_new_valid = 0

# Sample paired reads
pairs = {}
with pysam.AlignmentFile(BAM, "rb") as bam:
    count = 0
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in pairs:
            pairs[qname] = read
        else:
            r1 = pairs.pop(qname)
            r2 = read
            n_pairs += 1
            # Only process multi-exon fragments (those with gaps between reads)
            # Each read contributes exon blocks from CIGAR
            count += 1
            if count >= n_max:
                break

print(f"\nProcessed {n_pairs:,} pairs, {count:,} sampled")

# Instead of manual FL computation, let me just compare the FL model
# outputs directly by looking at the output stats
print("\n=== POST-FIX FL Model ===")
print("SPLICED_ANNOT: mode=294, mean=298.0, n_obs=2,892,014")
print("UNSPLICED:     mode=295, mean=288.6, n_obs=3,230,946")  
print("Global:        mode=294, mean=293.0, n_obs=6,123,020")
print("\n(PRE-FIX had SPLICED_ANNOT mode=1000)")

# Calculate what percentage of all mapped fragments are affected
total_frags = 9_999_650
print(f"\nTotal fragments: {total_frags:,}")
print(f"Unambiguous FL observations: 6,123,020 ({6123020/total_frags*100:.1f}%)")
print(f"Ambiguous FL: 1,060,687 ({1060687/total_frags*100:.1f}%)")
print(f"FL used for training: {6123020 + 1060687:,} ({(6123020+1060687)/total_frags*100:.1f}%)")

# The key question: how many unique mappers with CIGAR introns now have correct FL?
# Let's count from the BAM
print("\n=== Checking unique mapper FL validity ===")
pairs2 = {}
stats = Counter()
fl_vals = []

with pysam.AlignmentFile(BAM, "rb") as bam:
    count = 0
    for read in bam:
        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in pairs2:
            pairs2[qname] = [read]
        else:
            pairs2[qname].append(read)
            if len(pairs2[qname]) == 2:
                r1, r2 = pairs2[qname]
                del pairs2[qname]
                
                # Check NH tag
                try:
                    nh1 = r1.get_tag("NH")
                except KeyError:
                    nh1 = 1
                
                if nh1 > 1:
                    stats["multi"] += 1
                    continue
                
                # Unique mapper - check for CIGAR introns
                has_intron_r1 = any(op == 3 for op, _ in r1.cigartuples) if r1.cigartuples else False
                has_intron_r2 = any(op == 3 for op, _ in r2.cigartuples) if r2.cigartuples else False
                
                if has_intron_r1 or has_intron_r2:
                    stats["unique_spliced"] += 1
                else:
                    stats["unique_unspliced"] += 1
                    # For unspliced, compute genomic span
                    if r1.reference_name == r2.reference_name:
                        span = abs(max(r1.reference_end, r2.reference_end) - 
                                  min(r1.reference_start, r2.reference_start))
                        if span > 500:
                            stats["unspliced_large_span"] += 1
                        else:
                            stats["unspliced_normal_span"] += 1
                            fl_vals.append(span)
                
                count += 1
                if count >= 500000:
                    break

print(f"\nUnique mapper stats ({count:,} pairs):")
for k, v in sorted(stats.items()):
    print(f"  {k}: {v:,}")

if fl_vals:
    arr = np.array(fl_vals)
    print(f"\nUnspliced normal-span FL (n={len(arr):,}):")
    print(f"  mean: {arr.mean():.1f}, median: {np.median(arr):.1f}")
    print(f"  mode: {Counter(arr).most_common(1)[0]}")
