#!/usr/bin/env python3
"""Deeper analysis of gap types and correction failures."""
import numpy as np
import pandas as pd
import pysam
from collections import defaultdict

INDEX_DIR = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/rigel_index"
BAM_MM2 = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/gdna_none_ss_0.95_nrna_none/align_minimap2/reads_namesort.bam"
sj_df = pd.read_feather(f"{INDEX_DIR}/sj.feather")

# Pre-build SJ lookup
sj_by_ref = {}
for ref_name, group in sj_df.groupby("ref"):
    sj_by_ref[ref_name] = group[["start", "end", "t_index"]].values

def get_exon_blocks(read):
    blocks = []
    pos = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            blocks.append((pos, pos + length))
            pos += length
        elif op in (2, 3):
            pos += length
    return blocks

def get_introns(read):
    introns = []
    pos = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8, 2):
            pos += length
        elif op == 3:
            introns.append((pos, pos + length))
            pos += length
    return introns

mm2 = pysam.AlignmentFile(BAM_MM2, "rb")
pair_data = defaultdict(lambda: {"r1": [], "r2": []})
count = 0
for read in mm2:
    if read.is_unmapped or read.is_supplementary:
        continue
    count += 1
    if count > 1000000:
        break
    key = "r1" if read.is_read1 else "r2"
    pair_data[read.query_name][key].append(read)
mm2.close()

# Classify gaps
gap_types = {
    "small_same_exon": 0,        # gap < 50bp, no SJs — same exon gap
    "large_no_sj": 0,            # gap >= 50bp, no SJs at all — should be rare
    "sj_contained_exact": 0,     # SJ fully contained with tol=0
    "sj_near_miss_1bp": 0,       # SJ fails by 1bp
    "sj_near_miss_2_5bp": 0,     # SJ fails by 2-5bp
    "sj_near_miss_6_50bp": 0,    # SJ fails by 6-50bp
    "sj_far_miss": 0,            # SJ overlaps but extends far beyond gap
}

# Collect detailed stats for each category
gap_sizes_by_type = defaultdict(list)
upper_before_correction = []
upper_after_correction = []
no_correction_uppers = []

examples = {"sj_near_miss_1bp": [], "sj_far_miss": [], "large_no_sj": []}

for qname, pair in pair_data.items():
    r1_list = [r for r in pair["r1"] if not r.is_secondary]
    r2_list = [r for r in pair["r2"] if not r.is_secondary]
    if not r1_list or not r2_list:
        continue
    r1, r2 = r1_list[0], r2_list[0]
    if r1.reference_id != r2.reference_id:
        continue
    
    r1_blocks = get_exon_blocks(r1)
    r2_blocks = get_exon_blocks(r2)
    all_blocks = sorted(r1_blocks + r2_blocks)
    if len(all_blocks) < 2:
        continue
    
    all_introns = get_introns(r1) + get_introns(r2)
    cigar_intron_set = set((s, e) for s, e in all_introns)
    
    footprint = all_blocks[-1][1] - all_blocks[0][0]
    total_cigar_intron = sum(e - s for s, e in all_introns)
    upper = footprint - total_cigar_intron
    
    for i in range(len(all_blocks) - 1):
        gs = all_blocks[i][1]
        ge = all_blocks[i+1][0]
        if gs >= ge or (gs, ge) in cigar_intron_set:
            continue
        
        gap_size = ge - gs
        ref_name = r1.reference_name
        sj_data = sj_by_ref.get(ref_name)
        
        if sj_data is None or len(sj_data) == 0:
            if gap_size < 50:
                gap_types["small_same_exon"] += 1
            else:
                gap_types["large_no_sj"] += 1
                no_correction_uppers.append(upper)
            continue
        
        # Find overlapping SJs
        overlapping = sj_data[(sj_data[:, 0] < ge) & (sj_data[:, 1] > gs)]
        
        if len(overlapping) == 0:
            if gap_size < 50:
                gap_types["small_same_exon"] += 1
            else:
                gap_types["large_no_sj"] += 1
                no_correction_uppers.append(upper)
                if len(examples["large_no_sj"]) < 5:
                    examples["large_no_sj"].append({
                        "qname": qname, "ref": ref_name,
                        "gap": (gs, ge), "gap_size": gap_size,
                        "upper": upper, "footprint": footprint,
                    })
            continue
        
        # Check containment
        contained_exact = overlapping[
            (overlapping[:, 0] >= gs) & (overlapping[:, 1] <= ge)
        ]
        
        if len(contained_exact) > 0:
            gap_types["sj_contained_exact"] += 1
            # Compute correction
            max_correction = max(r[1] - r[0] for r in contained_exact)
            upper_before_correction.append(upper)
            upper_after_correction.append(upper - max_correction)
            continue
        
        # Find the near-miss: minimum overshoot
        start_overshoot = np.maximum(0, gs - overlapping[:, 0])
        end_overshoot = np.maximum(0, overlapping[:, 1] - ge)
        total_overshoot = start_overshoot + end_overshoot
        min_overshoot = total_overshoot.min()
        
        if min_overshoot <= 1:
            gap_types["sj_near_miss_1bp"] += 1
            no_correction_uppers.append(upper)
            if len(examples["sj_near_miss_1bp"]) < 5:
                best_idx = total_overshoot.argmin()
                examples["sj_near_miss_1bp"].append({
                    "qname": qname, "ref": ref_name,
                    "gap": (gs, ge), "gap_size": gap_size,
                    "sj": (int(overlapping[best_idx, 0]), int(overlapping[best_idx, 1])),
                    "start_over": int(start_overshoot[best_idx]),
                    "end_over": int(end_overshoot[best_idx]),
                    "upper": upper,
                })
        elif min_overshoot <= 5:
            gap_types["sj_near_miss_2_5bp"] += 1
            no_correction_uppers.append(upper)
        elif min_overshoot <= 50:
            gap_types["sj_near_miss_6_50bp"] += 1
            no_correction_uppers.append(upper)
        else:
            gap_types["sj_far_miss"] += 1
            no_correction_uppers.append(upper)
            if len(examples["sj_far_miss"]) < 5:
                best_idx = total_overshoot.argmin()
                examples["sj_far_miss"].append({
                    "qname": qname, "ref": ref_name,
                    "gap": (gs, ge), "gap_size": gap_size,
                    "sj": (int(overlapping[best_idx, 0]), int(overlapping[best_idx, 1])),
                    "start_over": int(start_overshoot[best_idx]),
                    "end_over": int(end_overshoot[best_idx]),
                    "upper": upper,
                })

print("=" * 70)
print("GAP CLASSIFICATION RESULTS")
print("=" * 70)
total = sum(gap_types.values())
for cat, cnt in gap_types.items():
    print(f"  {cat:30s}: {cnt:8,} ({100*cnt/max(total,1):5.1f}%)")
print(f"  {'TOTAL':30s}: {total:8,}")

print(f"\n## Uncorrected fragment length impact:")
if no_correction_uppers:
    nc = np.array(no_correction_uppers)
    print(f"  Uncorrected fragments: {len(nc):,}")
    print(f"  Uncorrected upper (FL proxy) stats:")
    print(f"    Mean: {nc.mean():.0f}")
    print(f"    Median: {np.median(nc):.0f}")
    print(f"    > 500: {(nc > 500).sum():,}")
    print(f"    > 1000: {(nc > 1000).sum():,}")

if upper_before_correction:
    bc = np.array(upper_before_correction)
    ac = np.array(upper_after_correction)
    print(f"\n  Corrected fragments: {len(bc):,}")
    print(f"  Before correction mean: {bc.mean():.0f}")
    print(f"  After correction mean: {ac.mean():.0f}")

print(f"\n## Examples of 1bp near-misses:")
for ex in examples["sj_near_miss_1bp"]:
    print(f"  {ex['qname'][:50]}...")
    print(f"    gap=[{ex['gap'][0]}, {ex['gap'][1]}) size={ex['gap_size']}")
    print(f"    sj=[{ex['sj'][0]}, {ex['sj'][1]}) size={ex['sj'][1]-ex['sj'][0]}")
    print(f"    start_over={ex['start_over']} end_over={ex['end_over']}")
    print(f"    upper={ex['upper']}")

print(f"\n## Examples of far misses:")
for ex in examples["sj_far_miss"]:
    print(f"  {ex['qname'][:50]}...")
    print(f"    gap=[{ex['gap'][0]}, {ex['gap'][1]}) size={ex['gap_size']}")
    print(f"    sj=[{ex['sj'][0]}, {ex['sj'][1]}) size={ex['sj'][1]-ex['sj'][0]}")
    print(f"    start_over={ex['start_over']} end_over={ex['end_over']}")
    print(f"    upper={ex['upper']}")

print(f"\n## Examples of large gaps with no SJ:")
for ex in examples["large_no_sj"]:
    print(f"  {ex['qname'][:50]}...")
    print(f"    gap=[{ex['gap'][0]}, {ex['gap'][1]}) size={ex['gap_size']}")
    print(f"    upper={ex['upper']} footprint={ex['footprint']}")
