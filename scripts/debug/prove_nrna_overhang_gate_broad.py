#!/usr/bin/env python3
"""Find fragments where minimap2 missed a splice junction, causing intronic overhang.

For EACH read (not fragment), check:
1. Read has no N-skip in CIGAR
2. Read's aligned block extends past an exon boundary of a matching transcript
3. This creates overhang for mRNA but not for nRNA

The key case: R1 is at [exon_end - 148, exon_end + 2] = 150M
but should have been [exon_end - 148, exon_end] + N-skip + [next_exon_start, next_exon_start + 2]
"""
import sys
sys.path.insert(0, "src")

import pysam
import numpy as np
import pyarrow.feather as pf
from collections import defaultdict, Counter

BASE = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v2"
COND = "gdna_none_ss_0.95_nrna_none"
RAW_BAM = f"{BASE}/{COND}/align_minimap2/reads_namesort.bam"
INDEX = f"{BASE}/rigel_index"

# Load transcript data
t_df = pf.read_table(f"{INDEX}/transcripts.feather").to_pandas()
iv_df = pf.read_table(f"{INDEX}/intervals.feather").to_pandas()

nrna_set = set(t_df.loc[t_df["is_synthetic_nrna"] == True, "t_index"].values)
mrna_t_df = t_df[t_df["is_synthetic_nrna"] != True].copy()

# Build exon lookup: for each (ref, strand), sorted list of exon boundaries
# with transcript index
print("Building exon boundary index...")
exon_ivs = iv_df[(iv_df["interval_type"] == 0) & (iv_df["t_index"] >= 0)].copy()
mrna_exon_ivs = exon_ivs[~exon_ivs["t_index"].isin(nrna_set)]

# Group by transcript, sort exons
t_exons = {}  # t_index -> sorted list of (start, end)
for _, row in mrna_exon_ivs.iterrows():
    ti = row["t_index"]
    if ti not in t_exons:
        t_exons[ti] = []
    t_exons[ti].append((row["start"], row["end"]))
for ti in t_exons:
    t_exons[ti].sort()

# Build intron boundaries for each transcript
t_introns = {}  # t_index -> list of (start, end)
for ti, exons in t_exons.items():
    introns = []
    for i in range(len(exons) - 1):
        introns.append((exons[i][1], exons[i+1][0]))
    if introns:
        t_introns[ti] = introns

# Build ref -> transcript index mapping
t_ref = dict(zip(t_df["t_index"], t_df["ref"]))

print(f"mRNA transcripts with introns: {len(t_introns)}")

# For efficiency, build an interval tree of exon-intron boundaries per chromosome
# Simple approach: for each chromosome, collect all exon boundaries
print("Building exon boundary set...")
from bisect import bisect_left, bisect_right

# Actually, let's just sample reads and check them against the transcript they came from
# We know the truth from the read name: ENST...:start-end:strand:fragid
print("\nScanning BAM for reads with intronic overhang...")

n_pairs = 0
n_has_overhang = 0
n_no_nsplice = 0
oh_examples = []
oh_stats = Counter()
mrna_oh_vals = []
nrna_oh_vals = []

pairs = {}
with pysam.AlignmentFile(RAW_BAM, "rb") as bam:
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        
        qname = read.query_name
        if qname not in pairs:
            pairs[qname] = read
            continue
        
        r1 = pairs.pop(qname)
        r2 = read
        n_pairs += 1
        
        # Parse truth transcript from read name
        parts = qname.split(":")
        if len(parts) < 4:
            continue
        truth_t_id = parts[0]
        
        # Find truth transcript index
        truth_rows = t_df[t_df["t_id"] == truth_t_id]
        if len(truth_rows) == 0:
            continue
        truth_t_idx = truth_rows.iloc[0]["t_index"]
        
        if truth_t_idx not in t_exons:
            continue
        if truth_t_idx not in t_introns:
            continue  # single-exon transcript
        
        exons = t_exons[truth_t_idx]
        truth_ref = t_ref[truth_t_idx]
        
        # Check each read individually
        for r in [r1, r2]:
            if r.reference_name != truth_ref:
                continue
            if not r.cigartuples:
                continue
            
            has_nsplice = any(op == 3 for op, _ in r.cigartuples)
            if has_nsplice:
                continue  # Skip reads with N-skip
            
            # Get aligned blocks
            blocks = []
            pos = r.reference_start
            for op, length in r.cigartuples:
                if op in (0, 7, 8):  # M/=/X
                    blocks.append((pos, pos + length))
                    pos += length
                elif op in (2, 3):  # D/N
                    pos += length
            
            if not blocks:
                continue
            
            # Check if any block extends past an exon into an intron
            read_start = blocks[0][0]
            read_end = blocks[-1][1]
            
            # Find overlapping exons
            exon_bp = 0
            for es, ee in exons:
                for bs, be in blocks:
                    overlap = min(be, ee) - max(bs, es)
                    if overlap > 0:
                        exon_bp += overlap
            
            rl = sum(be - bs for bs, be in blocks)
            oh = rl - exon_bp
            
            if oh > 0:
                n_has_overhang += 1
                mrna_oh_vals.append(oh)
                
                # Compute nRNA overhang (for the nRNA spanning this locus)
                # nRNA exon covers entire transcript span
                tx_start = exons[0][0]
                tx_end = exons[-1][1]
                # Find nRNA spanning this region
                nrna_ebp = 0
                for nrna_row in t_df[
                    (t_df["is_synthetic_nrna"] == True) & 
                    (t_df["ref"] == truth_ref)
                ].itertuples():
                    if nrna_row.start <= read_start and nrna_row.end >= read_end:
                        # This nRNA covers the read
                        nrna_ebp = rl  # entirely within nRNA "exon"
                        break
                
                nrna_oh = max(rl - nrna_ebp, 0)
                nrna_oh_vals.append(nrna_oh)
                
                if len(oh_examples) < 20:
                    oh_examples.append({
                        "qname": qname,
                        "truth_t": truth_t_id,
                        "ref": r.reference_name,
                        "pos": r.reference_start,
                        "cigar": r.cigarstring,
                        "rl": rl,
                        "exon_bp": exon_bp,
                        "mrna_oh": oh,
                        "nrna_oh": nrna_oh,
                        "read_span": f"[{read_start}, {read_end})",
                        "nearby_exon": str([(s,e) for s,e in exons if abs(s-read_start)<1000 or abs(e-read_end)<1000]),
                    })
                
                oh_stats[oh] += 1
        
        if n_pairs >= 500_000:
            break

print(f"\nTotal pairs processed: {n_pairs:,}")
print(f"Reads with intronic overhang (no N-skip): {n_has_overhang:,}")

if mrna_oh_vals:
    arr = np.array(mrna_oh_vals)
    print(f"\nmRNA overhang distribution:")
    print(f"  mean: {arr.mean():.1f}")
    print(f"  median: {np.median(arr):.0f}")
    print(f"  p90: {np.percentile(arr, 90):.0f}")
    print(f"  max: {arr.max()}")
    print(f"\n  Most common overhang values:")
    for val, cnt in sorted(oh_stats.items(), key=lambda x: -x[1])[:15]:
        print(f"    oh={val}: {cnt:,}")

if nrna_oh_vals:
    arr = np.array(nrna_oh_vals)
    print(f"\nnRNA overhang for same reads:")
    print(f"  all zero: {(arr == 0).all()}")
    print(f"  mean: {arr.mean():.1f}")

if oh_examples:
    print(f"\n=== EXAMPLES of intronic overhang ===")
    for ex in oh_examples[:10]:
        print(f"\n  {ex['qname']} (truth: {ex['truth_t']}):")
        print(f"    pos={ex['ref']}:{ex['pos']}, cigar={ex['cigar']}")
        print(f"    span={ex['read_span']}, rl={ex['rl']}")
        print(f"    exon_bp={ex['exon_bp']}, mrna_oh={ex['mrna_oh']}, nrna_oh={ex['nrna_oh']}")
        print(f"    nearby exons: {ex['nearby_exon']}")
        if ex['mrna_oh'] > 0 and ex['nrna_oh'] == 0:
            print(f"    >>> BUG: nRNA wins hard overhang gate (mrna_oh={ex['mrna_oh']}, nrna_oh=0)")
