#!/usr/bin/env python3
"""Directly prove the nRNA hard overhang gate bug by simulating the scoring step.

For a fragment where minimap2 doesn't detect the splice junction:
1. The read extends a few bp past an mRNA exon boundary into intronic space
2. For the mRNA transcript: e_bp = (read_within_exon), oh = rl - e_bp > 0
3. For the nRNA transcript: e_bp = full read (it's all "exon"), oh = 0
4. hard overhang gate keeps oh=0 only → nRNA wins → mRNA discarded

This script takes actual minimap2 reads and resolves them against both
mRNA and nRNA intervals to prove the mechanism.
"""
import sys
sys.path.insert(0, "src")

import pysam
import numpy as np
import pyarrow.feather as pf
from collections import defaultdict

BASE = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v2"
COND = "gdna_none_ss_0.95_nrna_none"
RAW_BAM = f"{BASE}/{COND}/align_minimap2/reads_namesort.bam"
INDEX = f"{BASE}/rigel_index"

# Load transcript and interval data
t_df = pf.read_table(f"{INDEX}/transcripts.feather").to_pandas()
iv_df = pf.read_table(f"{INDEX}/intervals.feather").to_pandas()

# Identify nRNA transcripts
nrna_set = set(t_df.loc[t_df["is_synthetic_nrna"] == True, "t_index"].values)
print(f"nRNA transcripts: {len(nrna_set)}")

# Build a simple exon interval index for a specific transcript
# Let's pick ENST00000000233.10 which we saw in examples
target_t_id = "ENST00000000233.10"
target_row = t_df[t_df["t_id"] == target_t_id]
if len(target_row) == 0:
    print(f"Transcript {target_t_id} not found!")
    sys.exit(1)
target_t_idx = target_row.iloc[0]["t_index"]
target_ref = target_row.iloc[0]["ref"]
print(f"\nTarget: {target_t_id} (t_index={target_t_idx}, ref={target_ref})")

# Get exon intervals for this transcript
target_exons = iv_df[(iv_df["t_index"] == target_t_idx) & (iv_df["interval_type"] == 0)]
target_exons = target_exons.sort_values("start")
print(f"Exons ({len(target_exons)}):")
for _, row in target_exons.iterrows():
    print(f"  [{row['start']}, {row['end']}) = {row['end'] - row['start']}bp")

# Get the nRNA transcript that covers this region
# Find nRNA transcripts on the same ref/strand
target_strand = target_row.iloc[0]["strand"]
nrna_candidates = t_df[
    (t_df["is_synthetic_nrna"] == True) &
    (t_df["ref"] == target_ref) &
    (t_df["strand"] == target_strand)
]
print(f"\nnRNA transcripts on {target_ref} strand {target_strand}: {len(nrna_candidates)}")

# Find the nRNA that overlaps this transcript's genomic region
target_start = target_exons.iloc[0]["start"]
target_end = target_exons.iloc[-1]["end"]
overlapping_nrna = nrna_candidates[
    (nrna_candidates["start"] <= target_end) &
    (nrna_candidates["end"] >= target_start)
]
print(f"Overlapping nRNA transcripts:")
for _, row in overlapping_nrna.iterrows():
    print(f"  {row['t_id']} (t_index={row['t_index']}): [{row['start']}, {row['end']})")

# Get nRNA exon intervals
if len(overlapping_nrna) > 0:
    nrna_t_idx = overlapping_nrna.iloc[0]["t_index"]
    nrna_exons = iv_df[(iv_df["t_index"] == nrna_t_idx) & (iv_df["interval_type"] == 0)]
    print(f"\nnRNA exons:")
    for _, row in nrna_exons.iterrows():
        print(f"  [{row['start']}, {row['end']}) = {row['end'] - row['start']}bp")

# Now trace actual reads that should map to this transcript
print("\n" + "="*60)
print("TRACING READS FROM MINIMAP2 BAM")
print("="*60)

# Read from raw BAM
examples_no_nsplice = []
examples_with_nsplice = []

with pysam.AlignmentFile(RAW_BAM, "rb") as bam:
    # Fetch reads aligned to the target region
    # Since name-sorted, we can't fetch by region directly
    # Instead, iterate and filter by read name prefix
    pairs = {}
    count = 0
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        # Filter for reads from our target transcript
        if not qname.startswith(target_t_id + ":"):
            if qname not in pairs:
                pairs[qname] = None  # skip
            continue
        
        if qname not in pairs:
            pairs[qname] = read
        elif pairs[qname] is not None:
            r1 = pairs.pop(qname)
            r2 = read
            
            # Parse CIGAR for N-skips
            def has_n_skip(r):
                return any(op == 3 for op, _ in r.cigartuples) if r.cigartuples else False
            
            r1_nsplice = has_n_skip(r1)
            r2_nsplice = has_n_skip(r2)
            
            # Compute exon blocks from alignment
            def get_blocks(r):
                """Get aligned blocks (M/=/X operations)."""
                blocks = []
                pos = r.reference_start
                for op, length in r.cigartuples:
                    if op in (0, 7, 8):  # M, =, X
                        blocks.append((pos, pos + length))
                        pos += length
                    elif op in (2, 3):  # D, N
                        pos += length
                    elif op in (1, 4, 5):  # I, S, H
                        pass
                return blocks
            
            def compute_exon_overlap(blocks, exon_start, exon_end):
                """Compute how many bases of the blocks overlap [exon_start, exon_end)."""
                total = 0
                for bs, be in blocks:
                    overlap = min(be, exon_end) - max(bs, exon_start)
                    if overlap > 0:
                        total += overlap
                return total
            
            r1_blocks = get_blocks(r1)
            r2_blocks = get_blocks(r2)
            all_blocks = sorted(r1_blocks + r2_blocks)
            
            # Compute e_bp for mRNA vs nRNA
            # mRNA: sum overlap with each exon
            mrna_ebp = 0
            for _, row in target_exons.iterrows():
                for bs, be in all_blocks:
                    overlap = min(be, row["end"]) - max(bs, row["start"])
                    if overlap > 0:
                        mrna_ebp += overlap
            
            # nRNA: sum overlap with nRNA "exon" (entire span)
            nrna_ebp = 0
            if len(overlapping_nrna) > 0:
                nrna_start = overlapping_nrna.iloc[0]["start"]
                nrna_end = overlapping_nrna.iloc[0]["end"]
                for bs, be in all_blocks:
                    overlap = min(be, nrna_end) - max(bs, nrna_start)
                    if overlap > 0:
                        nrna_ebp += overlap
            
            # Read length (sum of aligned bases from CIGAR)
            def read_len_from_cigar(r):
                return sum(l for op, l in r.cigartuples if op in (0, 7, 8))  # M/=/X
            
            rl = read_len_from_cigar(r1)  # per-read (R1 used for scoring)
            
            info = {
                "qname": qname,
                "r1_ref": r1.reference_name,
                "r1_pos": r1.reference_start,
                "r1_cigar": r1.cigarstring,
                "r2_cigar": r2.cigarstring,
                "r1_nsplice": r1_nsplice,
                "r2_nsplice": r2_nsplice,
                "rl": rl,
                "mrna_ebp": mrna_ebp,
                "nrna_ebp": nrna_ebp,
                "mrna_oh": max(rl - mrna_ebp, 0),
                "nrna_oh": max(rl - nrna_ebp, 0),
                "all_blocks": all_blocks,
            }
            
            if not r1_nsplice and not r2_nsplice:
                examples_no_nsplice.append(info)
            else:
                examples_with_nsplice.append(info)
            
            count += 1
        
        if count >= 1000:
            break

print(f"\nFragments from {target_t_id}: {count}")
print(f"  With N-skip in either read: {len(examples_with_nsplice)}")
print(f"  No N-skip in either read: {len(examples_no_nsplice)}")

# Show examples of fragments WITHOUT N-skip (the bug cases)
print(f"\n--- Fragments WITHOUT N-skip (nRNA should win hard overhang gate) ---")
if examples_no_nsplice:
    for ex in examples_no_nsplice[:10]:
        print(f"\n  {ex['qname']}:")
        print(f"    R1: {ex['r1_ref']}:{ex['r1_pos']} CIGAR={ex['r1_cigar']}")
        print(f"    R2 CIGAR: {ex['r2_cigar']}")
        print(f"    rl={ex['rl']}, mrna_ebp={ex['mrna_ebp']}, nrna_ebp={ex['nrna_ebp']}")
        print(f"    mrna_oh={ex['mrna_oh']}, nrna_oh={ex['nrna_oh']}")
        if ex['mrna_oh'] > 0 and ex['nrna_oh'] == 0:
            print(f"    >>> nRNA WINS hard overhang gate (mrna_oh={ex['mrna_oh']} > nrna_oh=0)")
        elif ex['mrna_oh'] == 0 and ex['nrna_oh'] == 0:
            print(f"    >>> TIE (both oh=0, both survive)")
        else:
            print(f"    >>> mrna wins or both equal")

# Show examples WITH N-skip for comparison
print(f"\n--- Fragments WITH N-skip (normal case) ---")
if examples_with_nsplice:
    for ex in examples_with_nsplice[:5]:
        print(f"\n  {ex['qname']}:")
        print(f"    R1: {ex['r1_ref']}:{ex['r1_pos']} CIGAR={ex['r1_cigar']}")
        print(f"    R2 CIGAR: {ex['r2_cigar']}")
        print(f"    rl={ex['rl']}, mrna_ebp={ex['mrna_ebp']}, nrna_ebp={ex['nrna_ebp']}")
        print(f"    mrna_oh={ex['mrna_oh']}, nrna_oh={ex['nrna_oh']}")

# Summary statistics
print(f"\n--- hard overhang gate OUTCOME SUMMARY (no N-skip fragments) ---")
nrna_wins = sum(1 for ex in examples_no_nsplice if ex['mrna_oh'] > 0 and ex['nrna_oh'] == 0)
tie = sum(1 for ex in examples_no_nsplice if ex['mrna_oh'] == 0 and ex['nrna_oh'] == 0)
mrna_wins = sum(1 for ex in examples_no_nsplice if ex['mrna_oh'] == 0 and ex['nrna_oh'] > 0)
other = len(examples_no_nsplice) - nrna_wins - tie - mrna_wins
print(f"  nRNA wins hard overhang gate: {nrna_wins} ({100*nrna_wins/max(len(examples_no_nsplice),1):.1f}%)")
print(f"  Tie (both oh=0): {tie} ({100*tie/max(len(examples_no_nsplice),1):.1f}%)")
print(f"  mRNA wins: {mrna_wins}")
print(f"  Other: {other}")
