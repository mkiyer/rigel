#!/usr/bin/env python3
"""Trace minimap2-aligned reads with inflated fragment lengths.

Looks at the annotated BAM from the benchmark to find fragments
where the fragment length is vastly inflated compared to the oracle.
"""
import pysam
import sys
from collections import Counter, defaultdict

BAM_MM2 = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/gdna_none_ss_0.95_nrna_none/align_minimap2/reads_namesort.bam"
BAM_ORACLE = "/Users/mkiyer/Downloads/rigel_runs/sim_pristine/gdna_none_ss_0.95_nrna_none/sim_oracle.bam"

def get_exon_blocks(read):
    """Extract aligned exon blocks from CIGAR."""
    blocks = []
    pos = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):  # M, =, X — consumes both
            blocks.append((pos, pos + length))
            pos += length
        elif op == 2:  # D — consumes ref
            pos += length
        elif op == 3:  # N — ref skip (intron)
            pos += length
        elif op in (1, 4, 5):  # I, S, H — don't consume ref
            pass
    return blocks

def get_introns(read):
    """Extract N-skip introns from CIGAR."""
    introns = []
    pos = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            pos += length
        elif op == 2:  # D
            pos += length
        elif op == 3:  # N — intron
            introns.append((pos, pos + length))
            pos += length
        elif op in (1, 4, 5):
            pass
    return introns

def compute_read_length(read):
    """Sum of M/D/=/X ops (aligned exon span)."""
    total = 0
    for op, length in read.cigartuples:
        if op in (0, 2, 7, 8):  # M, D, =, X
            total += length
    return total

print("=" * 80)
print("TRACING MINIMAP2 FRAGMENT LENGTH INFLATION")
print("=" * 80)

# ── 1. Sample some reads from minimap2 BAM ──
print("\n## 1. Sampling minimap2 reads to understand CIGAR patterns...")

mm2 = pysam.AlignmentFile(BAM_MM2, "rb")

# Collect paired reads grouped by qname
pair_data = defaultdict(lambda: {"r1": [], "r2": []})
count = 0
max_reads = 500000

for read in mm2:
    if read.is_unmapped or read.is_supplementary:
        continue
    count += 1
    if count > max_reads:
        break
    key = "r1" if read.is_read1 else "r2"
    pair_data[read.query_name][key].append(read)

mm2.close()

# Analyze fragment geometry
print(f"  Collected {len(pair_data)} read pairs")

# For primary alignments only, compute fragment geometry
n_split_cigar = 0
n_no_split = 0
n_both_split = 0
n_one_split = 0

genomic_spans = []
intron_counts = []
gap_sizes_no_intron = []
gap_sizes_with_intron = []

for qname, pair in list(pair_data.items())[:200000]:
    r1_list = pair["r1"]
    r2_list = pair["r2"]
    
    # Use primary alignments only
    r1_primary = [r for r in r1_list if not r.is_secondary]
    r2_primary = [r for r in r2_list if not r.is_secondary]
    
    if not r1_primary or not r2_primary:
        continue
    
    r1 = r1_primary[0]
    r2 = r2_primary[0]
    
    if r1.reference_id != r2.reference_id:
        continue
    
    r1_blocks = get_exon_blocks(r1)
    r2_blocks = get_exon_blocks(r2)
    r1_introns = get_introns(r1)
    r2_introns = get_introns(r2)
    
    all_introns = r1_introns + r2_introns
    
    r1_has_n = len(r1_introns) > 0
    r2_has_n = len(r2_introns) > 0
    
    if r1_has_n and r2_has_n:
        n_both_split += 1
    elif r1_has_n or r2_has_n:
        n_one_split += 1
    else:
        n_no_split += 1
    
    # Compute genomic span
    all_blocks = r1_blocks + r2_blocks
    if not all_blocks:
        continue
    min_start = min(b[0] for b in all_blocks)
    max_end = max(b[1] for b in all_blocks)
    genomic_span = max_end - min_start
    
    total_intron = sum(e - s for s, e in all_introns)
    upper = genomic_span - total_intron
    
    genomic_spans.append(genomic_span)
    intron_counts.append(len(all_introns))
    
    # Check for gap between R1 and R2 ends
    r1_extent = (min(b[0] for b in r1_blocks), max(b[1] for b in r1_blocks))
    r2_extent = (min(b[0] for b in r2_blocks), max(b[1] for b in r2_blocks))
    
    # Determine gap
    gap_start = min(r1_extent[1], r2_extent[1])
    gap_end = max(r1_extent[0], r2_extent[0])
    
    if gap_end > gap_start:
        gap_size = gap_end - gap_start
        if all_introns:
            gap_sizes_with_intron.append(gap_size)
        else:
            gap_sizes_no_intron.append(gap_size)

print(f"  Both reads have N-skip: {n_both_split:,}")
print(f"  One read has N-skip: {n_one_split:,}")
print(f"  Neither read has N-skip: {n_no_split:,}")
print(f"  Total pairs analyzed: {n_both_split + n_one_split + n_no_split:,}")

import numpy as np

gs = np.array(genomic_spans)
print(f"\n  Genomic span stats:")
print(f"    Mean: {gs.mean():.0f}")
print(f"    Median: {np.median(gs):.0f}")
print(f"    Mode (approx): {np.bincount(gs.clip(0, 2000)).argmax()}")
print(f"    {(gs > 500).sum():,} / {len(gs):,} have genomic span > 500")
print(f"    {(gs > 1000).sum():,} / {len(gs):,} have genomic span > 1000")

ic = np.array(intron_counts)
print(f"\n  Intron count per fragment:")
for n in range(4):
    cnt = (ic == n).sum()
    print(f"    {n} introns: {cnt:,} ({100*cnt/len(ic):.1f}%)")

# ── 2. Focus on fragments with large genomic span but no introns in CIGAR ──
print("\n## 2. Fragments with large genomic span but NO CIGAR introns...")
print("  (These are the fragments where minimap2 doesn't report the intron)")

large_no_intron = 0
example_frags = []
count = 0

mm2 = pysam.AlignmentFile(BAM_MM2, "rb")
pair_data2 = defaultdict(lambda: {"r1": [], "r2": []})

for read in mm2:
    if read.is_unmapped or read.is_supplementary:
        continue
    count += 1
    if count > 500000:
        break
    key = "r1" if read.is_read1 else "r2"
    pair_data2[read.query_name][key].append(read)

mm2.close()

for qname, pair in pair_data2.items():
    r1_list = pair["r1"]
    r2_list = pair["r2"]
    r1_primary = [r for r in r1_list if not r.is_secondary]
    r2_primary = [r for r in r2_list if not r.is_secondary]
    if not r1_primary or not r2_primary:
        continue
    r1 = r1_primary[0]
    r2 = r2_primary[0]
    if r1.reference_id != r2.reference_id:
        continue
    
    r1_blocks = get_exon_blocks(r1)
    r2_blocks = get_exon_blocks(r2)
    all_blocks = r1_blocks + r2_blocks
    if not all_blocks:
        continue
    
    all_introns = get_introns(r1) + get_introns(r2)
    
    min_start = min(b[0] for b in all_blocks)
    max_end = max(b[1] for b in all_blocks)
    genomic_span = max_end - min_start
    total_intron = sum(e - s for s, e in all_introns)
    upper = genomic_span - total_intron
    
    if upper > 500 and len(all_introns) == 0:
        large_no_intron += 1
        if len(example_frags) < 20:
            # Get the original transcript from qname  
            t_id = qname.split(":")[0]
            example_frags.append({
                "qname": qname,
                "t_id": t_id,
                "ref": r1.reference_name,
                "r1_blocks": r1_blocks,
                "r2_blocks": r2_blocks,
                "r1_cigar": r1.cigarstring,
                "r2_cigar": r2.cigarstring,
                "genomic_span": genomic_span,
                "upper": upper,
                "r1_pos": r1.reference_start,
                "r2_pos": r2.reference_start,
            })

print(f"  Large span (>500) with no CIGAR introns: {large_no_intron:,}")
print(f"\n  Examples:")
for f in example_frags[:10]:
    gap_start = max(b[1] for b in f["r1_blocks"])
    gap_end = min(b[0] for b in f["r2_blocks"])
    if gap_start > gap_end:
        gap_start = max(b[1] for b in f["r2_blocks"])
        gap_end = min(b[0] for b in f["r1_blocks"])
    gap_size = max(0, gap_end - gap_start)
    print(f"\n    {f['qname'][:60]}...")
    print(f"      Source transcript: {f['t_id']}")
    print(f"      Ref: {f['ref']}")
    print(f"      R1 CIGAR: {f['r1_cigar']}, blocks: {f['r1_blocks']}")
    print(f"      R2 CIGAR: {f['r2_cigar']}, blocks: {f['r2_blocks']}")
    print(f"      Genomic span: {f['genomic_span']}")
    print(f"      Upper (no intron subtracted): {f['upper']}")
    print(f"      Gap between mates: [{gap_start}, {gap_end}) size={gap_size}")

print("\n" + "=" * 80)
print("END OF TRACE")
print("=" * 80)
