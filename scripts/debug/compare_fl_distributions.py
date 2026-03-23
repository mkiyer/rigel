#!/usr/bin/env python3
"""Compare FL distributions between pre-fix and post-fix annotated BAMs.

Focus on reads that resolved to transcripts (ZP != intergenic).
"""
import pysam
import numpy as np
from collections import Counter

BASE_V1 = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine"
BASE_V2 = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v2"
COND = "gdna_none_ss_0.95_nrna_none"

def extract_fl_stats(bam_path, n_max=2000000):
    """Extract ZL (fragment length) stats from annotated BAM."""
    fl_values = []
    fl_counter = Counter()
    zp_counter = Counter()
    n_reads = 0
    n_r1 = 0
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            n_reads += 1
            if n_reads > n_max:
                break
            
            if not read.is_read2:  # R1 or unpaired
                n_r1 += 1
                try:
                    zp = read.get_tag("ZP")
                    zl = read.get_tag("ZL")
                    zp_counter[zp] += 1
                    if zl > 0:
                        fl_values.append(zl)
                        if zl <= 500:
                            fl_counter["<=500"] += 1
                        elif zl <= 1000:
                            fl_counter["501-1000"] += 1
                        elif zl <= 2000:
                            fl_counter["1001-2000"] += 1
                        else:
                            fl_counter[">2000"] += 1
                    else:
                        fl_counter["<=0_or_missing"] += 1
                except KeyError:
                    fl_counter["no_ZL_tag"] += 1
    
    return fl_values, fl_counter, zp_counter, n_r1

for label, base in [("PRE-FIX", BASE_V1), ("POST-FIX", BASE_V2)]:
    bam_path = f"{base}/{COND}/align_minimap2/annotated_default.bam"
    print(f"\n{'='*60}")
    print(f"=== {label} ===")
    print(f"{'='*60}")
    
    fl_values, fl_counter, zp_counter, n_r1 = extract_fl_stats(bam_path)
    
    print(f"\nTotal R1 reads processed: {n_r1:,}")
    print(f"\nPartition (ZP) distribution:")
    for zp, count in sorted(zp_counter.items(), key=lambda x: -x[1]):
        print(f"  {zp}: {count:,} ({100*count/n_r1:.1f}%)")
    
    print(f"\nFL (ZL) distribution:")
    for bucket, count in sorted(fl_counter.items()):
        print(f"  {bucket}: {count:,}")
    
    if fl_values:
        arr = np.array(fl_values)
        print(f"\nFL statistics (ZL > 0):")
        print(f"  n: {len(arr):,}")
        print(f"  mean: {arr.mean():.1f}")
        print(f"  median: {np.median(arr):.1f}")
        print(f"  mode (binned): {Counter(arr).most_common(1)[0]}")
        print(f"  std: {arr.std():.1f}")
        print(f"  p5: {np.percentile(arr, 5):.0f}")
        print(f"  p25: {np.percentile(arr, 25):.0f}")
        print(f"  p50: {np.percentile(arr, 50):.0f}")
        print(f"  p75: {np.percentile(arr, 75):.0f}")
        print(f"  p95: {np.percentile(arr, 95):.0f}")
        print(f"  p99: {np.percentile(arr, 99):.0f}")
        print(f"  max: {arr.max():.0f}")
        print(f"  fraction > 500: {(arr > 500).mean():.4f}")
        print(f"  fraction > 1000: {(arr > 1000).mean():.4f}")
