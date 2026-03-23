#!/usr/bin/env python3
"""Verify the 1bp containment failure pattern for SJ gap correction."""
import numpy as np
import pandas as pd

INDEX_DIR = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/rigel_index"
sj_df = pd.read_feather(f"{INDEX_DIR}/sj.feather")

# The bug: SJ [8946405, 8949487) with gap [8946405, 8949486)
# The SJ end (8949487) > gap end (8949486) — fails containment by 1bp
# 
# Why? The gap is computed as:
#   gap_start = exon_block[i].end     (R1 exon block end)
#   gap_end   = exon_block[i+1].start (R2 exon block start)
#
# For a spliced alignment, the intron (SJ) is:
#   sj_start = exon_block.end         (same as gap_start — match!)
#   sj_end   = next_exon_block.start  (same as gap_end if read spans junction)
#
# BUT when R2 starts inside the next exon (not at the junction boundary),
# the R2 block start may be INSIDE the exon, not at the exon boundary.
# 
# In this case:
#   Exon in annotation ends at 8946405 (intron starts at 8946405)
#   Intron ends at 8949487 (next exon starts at 8949487)
#   R1 ends at 8946405 (matches intron start perfectly)
#   R2 starts at 8949486 (1bp BEFORE the intron end!)
#
# Wait — R2 starts at 8949486 but the intron ends at 8949487?
# That means R2 starts INSIDE the intron! That shouldn't be possible
# unless minimap2 is aligning R2 slightly wrong, or the read overlaps
# the junction boundary by 1bp but minimap2 doesn't report an N-skip.
#
# Actually, let me reconsider. The SJ coordinates in the index are
# [start, end) meaning the intron spans from position 'start' to 'end' (exclusive).
# The exon AFTER the intron starts at position 'end'.
# If R2 aligns starting at position 8949486, and the next exon starts at 8949487,
# then R2 starts 1bp BEFORE the exon boundary.
#
# This means either:
# a) Minimap2 is soft-clipping the first base or misaligning by 1bp
# b) The R2 mate actually does overlap into the intron by 1bp
# c) There's an off-by-one in coordinate systems
#
# Let's check several examples to see if this is a consistent pattern.

print("=" * 70)
print("CONTAINMENT FAILURE ANALYSIS")
print("=" * 70)

# Let's check: for the SJs of ENST00000000412.8, what are the exon boundaries?
t_df = pd.read_feather(f"{INDEX_DIR}/transcripts.feather")
t_row = t_df[t_df["t_id"] == "ENST00000000412.8"].iloc[0]
t_index = t_row["t_index"]

# Get exon coordinates from the intervals file
iv_df = pd.read_feather(f"{INDEX_DIR}/intervals.feather")
print(f"\nInterval columns: {list(iv_df.columns)}")
print(iv_df.head())

# The SJ for the target is [8946405, 8949487).
# If the exon before ends at 8946405 and the exon after starts at 8949487,
# then a read aligned from 8949486 overlaps into the intron by 1bp.
# This means the gap [gap_start=8946405, gap_end=8949486) is 1bp SHORT
# of containing the intron [8946405, 8949487).

# The CONTAINMENT CHECK in resolve_context.h is:
#   if (hs >= gs && he <= ge)
# where hs = sj start, he = sj end, gs = gap start, ge = gap end
#
# For our case: hs=8946405 >= gs=8946405 (TRUE) AND he=8949487 <= ge=8949486 (FALSE!)
# Fails by exactly 1bp.

# This is NOT a coordinate system bug. It's a fundamental issue:
# The gap between mates doesn't fully contain the intron because the R2 mate
# starts at a position that's within the last 1-2 bases of the intron (soft-clipped
# or not detected by minimap2).

# The real question: should the containment check use a TOLERANCE?
# i.e., hs >= gs - tol && he <= ge + tol
# where tol could be a small number like 5-10bp

# Let's quantify: how many SJ gap corrections would succeed with
# various tolerances?

import pysam
from collections import defaultdict

BAM_MM2 = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/gdna_none_ss_0.95_nrna_none/align_minimap2/reads_namesort.bam"

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
        if op in (0, 7, 8):
            pos += length
        elif op == 2:
            pos += length
        elif op == 3:
            introns.append((pos, pos + length))
            pos += length
    return introns

# Sample fragments and check gap containment
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

# For each fragment with a gap between mates and no CIGAR introns in the gap:
# Check the SJ index for matching introns and measure the containment failure
tolerances = [0, 1, 2, 3, 5, 10]
corrections_by_tol = {t: 0 for t in tolerances}
total_gapped = 0
total_large_gapped = 0  # genomic span > 500
sj_chr12 = sj_df[sj_df["ref"] == "chr12"]

# Pre-build SJ lookup for all chromosomes
sj_by_ref = {}
for ref_name, group in sj_df.groupby("ref"):
    sj_by_ref[ref_name] = group[["start", "end", "t_index"]].values

gap_miss_reasons = {"no_sj_overlap": 0, "sj_not_contained": 0, 
                     "sj_contained_no_t": 0, "corrected": 0}
near_miss_sizes = []  # how many bp off for "almost contained" SJs

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
    all_blocks = r1_blocks + r2_blocks
    if not all_blocks:
        continue
    
    r1_introns = get_introns(r1)
    r2_introns = get_introns(r2)
    all_introns = r1_introns + r2_introns
    
    # Sort blocks
    sorted_blocks = sorted(all_blocks)
    
    # Find gaps between consecutive blocks that aren't CIGAR introns
    cigar_intron_set = set((s, e) for s, e in all_introns)
    
    for i in range(len(sorted_blocks) - 1):
        gs = sorted_blocks[i][1]   # end of block i
        ge = sorted_blocks[i+1][0]  # start of block i+1
        if gs >= ge:
            continue
        if (gs, ge) in cigar_intron_set:
            continue
        
        gap_size = ge - gs
        total_gapped += 1
        
        genomic_span = sorted_blocks[-1][1] - sorted_blocks[0][0]
        total_intron_size = sum(e - s for s, e in all_introns)
        upper = genomic_span - total_intron_size
        if upper > 500:
            total_large_gapped += 1
        
        # Check SJ index for this gap
        ref_name = r1.reference_name
        sj_data = sj_by_ref.get(ref_name)
        if sj_data is None:
            gap_miss_reasons["no_sj_overlap"] += 1
            continue
        
        # Find SJs that overlap this gap
        overlapping = sj_data[
            (sj_data[:, 0] < ge) & (sj_data[:, 1] > gs)
        ]
        
        if len(overlapping) == 0:
            gap_miss_reasons["no_sj_overlap"] += 1
            continue
        
        # Check containment at various tolerances
        found_any = False
        for tol in tolerances:
            contained = overlapping[
                (overlapping[:, 0] >= gs - tol) & (overlapping[:, 1] <= ge + tol)
            ]
            if len(contained) > 0:
                corrections_by_tol[tol] += 1
                if tol == 0:
                    gap_miss_reasons["corrected"] += 1
                found_any = True
        
        if not found_any:
            gap_miss_reasons["sj_not_contained"] += 1
        elif corrections_by_tol[0] < corrections_by_tol[1]:
            # This gap had SJs that were contained with tol=1 but not tol=0
            # Record the near-miss
            for sj_row in overlapping:
                start_miss = max(0, gs - sj_row[0])  # how much sj extends before gap
                end_miss = max(0, sj_row[1] - ge)    # how much sj extends past gap
                if start_miss <= 10 or end_miss <= 10:
                    near_miss_sizes.append((start_miss, end_miss))

print(f"\n## Gap Analysis Results (from {len(pair_data):,} pairs)")
print(f"  Total gaps between mate blocks (no CIGAR intron): {total_gapped:,}")
print(f"  Large gaps (upper > 500): {total_large_gapped:,}")
print(f"\n  Gap miss reasons:")
for reason, cnt in gap_miss_reasons.items():
    print(f"    {reason}: {cnt:,}")
print(f"\n  Corrections by tolerance:")
for tol, cnt in sorted(corrections_by_tol.items()):
    print(f"    tol={tol:2d}bp: {cnt:,} corrections ({100*cnt/max(total_gapped,1):.1f}%)")

if near_miss_sizes:
    nm = np.array(near_miss_sizes)
    print(f"\n  Near-miss analysis ({len(nm)} cases):")
    print(f"    Start overshoot (SJ extends before gap): "
          f"mean={nm[:, 0].mean():.1f}, max={nm[:, 0].max()}")
    print(f"    End overshoot (SJ extends past gap): "
          f"mean={nm[:, 1].mean():.1f}, max={nm[:, 1].max()}")

print("\n" + "=" * 70)
print("END OF CONTAINMENT ANALYSIS")
print("=" * 70)
