#!/usr/bin/env python3
"""Optimized proof of nRNA hard overhang gate bug.

For reads where minimap2 missed a splice junction, the read extends into
intronic space. mRNA has overhang > 0, nRNA has overhang = 0 → nRNA wins hard overhang gate.
"""
import sys
sys.path.insert(0, "src")

import pysam
import numpy as np
import pyarrow.feather as pf
from collections import Counter

BASE = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v2"
COND = "gdna_none_ss_0.95_nrna_none"
RAW_BAM = f"{BASE}/{COND}/align_minimap2/reads_namesort.bam"
INDEX = f"{BASE}/rigel_index"

# Load transcript data
t_df = pf.read_table(f"{INDEX}/transcripts.feather").to_pandas()
iv_df = pf.read_table(f"{INDEX}/intervals.feather").to_pandas()

nrna_set = set(t_df.loc[t_df["is_synthetic_nrna"] == True, "t_index"].values)

# Build fast lookup: t_id -> t_index
t_id_to_idx = dict(zip(t_df["t_id"], t_df["t_index"]))

# Build exon lookup: t_index -> sorted list of (start, end)
print("Building exon index...")
exon_ivs = iv_df[(iv_df["interval_type"] == 0) & (iv_df["t_index"] >= 0)]
mrna_exon_ivs = exon_ivs[~exon_ivs["t_index"].isin(nrna_set)]

t_exons = {}
for ti, start, end in zip(mrna_exon_ivs["t_index"].values,
                           mrna_exon_ivs["start"].values,
                           mrna_exon_ivs["end"].values):
    if ti not in t_exons:
        t_exons[ti] = []
    t_exons[ti].append((int(start), int(end)))
for ti in t_exons:
    t_exons[ti].sort()

# Keep only multi-exon transcripts
t_exons = {ti: exons for ti, exons in t_exons.items() if len(exons) > 1}

# Build nRNA lookup: (ref, strand) -> sorted list of (start, end, t_index)
nrna_spans = {}
nrna_df = t_df[t_df["is_synthetic_nrna"] == True]
for ref, strand, start, end, ti in zip(
    nrna_df["ref"].values, nrna_df["strand"].values,
    nrna_df["start"].values, nrna_df["end"].values,
    nrna_df["t_index"].values
):
    key = (ref, int(strand))
    if key not in nrna_spans:
        nrna_spans[key] = []
    nrna_spans[key].append((int(start), int(end), int(ti)))
for key in nrna_spans:
    nrna_spans[key].sort()

# Transcript ref/strand lookup
t_ref = dict(zip(t_df["t_index"], t_df["ref"]))
t_strand = dict(zip(t_df["t_index"], t_df["strand"]))

print(f"Multi-exon mRNA transcripts: {len(t_exons)}")
print(f"nRNA span groups: {len(nrna_spans)}")

# Scan BAM
print("\nScanning BAM...")
n_pairs = 0
n_overhang_reads = 0
oh_vals = []
oh_examples = []
oh_counter = Counter()
gate_outcomes = Counter()

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
        colon = qname.find(":")
        if colon < 0:
            continue
        truth_t_id = qname[:colon]
        truth_t_idx = t_id_to_idx.get(truth_t_id, -1)
        if truth_t_idx < 0 or truth_t_idx not in t_exons:
            continue
        
        exons = t_exons[truth_t_idx]
        ref = t_ref[truth_t_idx]
        strand = t_strand[truth_t_idx]
        
        # Check each read
        for r in [r1, r2]:
            if not r.cigartuples or r.reference_name != ref:
                continue
            
            has_nsplice = any(op == 3 for op, _ in r.cigartuples)
            if has_nsplice:
                continue
            
            # Get aligned blocks
            blocks = []
            pos = r.reference_start
            for op, length in r.cigartuples:
                if op in (0, 7, 8):
                    blocks.append((pos, pos + length))
                    pos += length
                elif op in (2, 3):
                    pos += length
            
            if not blocks:
                continue
            
            # Compute exon overlap
            exon_bp = 0
            for es, ee in exons:
                for bs, be in blocks:
                    ov = min(be, ee) - max(bs, es)
                    if ov > 0:
                        exon_bp += ov
            
            rl = sum(be - bs for bs, be in blocks)
            mrna_oh = max(rl - exon_bp, 0)
            
            if mrna_oh > 0:
                n_overhang_reads += 1
                oh_vals.append(mrna_oh)
                oh_counter[mrna_oh] += 1
                
                # Check nRNA overhang
                read_start = blocks[0][0]
                read_end = blocks[-1][1]
                nrna_oh = rl  # worst case: no nRNA covers it
                key = (ref, strand)
                if key in nrna_spans:
                    for ns, ne, nti in nrna_spans[key]:
                        if ns <= read_start and ne >= read_end:
                            nrna_oh = 0
                            break
                
                if mrna_oh > 0 and nrna_oh == 0:
                    gate_outcomes["nrna_wins"] += 1
                elif mrna_oh == 0 and nrna_oh == 0:
                    gate_outcomes["tie"] += 1
                elif mrna_oh > 0 and nrna_oh > 0:
                    gate_outcomes["both_oh"] += 1
                else:
                    gate_outcomes["mrna_wins"] += 1
                
                if len(oh_examples) < 20 and mrna_oh > 0 and nrna_oh == 0:
                    oh_examples.append({
                        "qname": qname,
                        "truth_t": truth_t_id,
                        "pos": f"{r.reference_name}:{r.reference_start}",
                        "cigar": r.cigarstring,
                        "rl": rl,
                        "exon_bp": exon_bp,
                        "mrna_oh": mrna_oh,
                        "nrna_oh": nrna_oh,
                        "blocks": blocks,
                        "nearby_exons": [(s,e) for s,e in exons 
                                        if abs(s-read_start)<500 or abs(e-read_end)<500 
                                        or (s<=read_end and e>=read_start)],
                    })
        
        if n_pairs >= 200_000:
            break
        if n_pairs % 100_000 == 0:
            print(f"  {n_pairs:,} pairs, {n_overhang_reads:,} overhang reads")

print(f"\nTotal pairs: {n_pairs:,}")
print(f"Reads with intronic overhang (no N-skip): {n_overhang_reads:,}")

if oh_vals:
    arr = np.array(oh_vals)
    print(f"\nmRNA overhang distribution:")
    print(f"  mean: {arr.mean():.1f}")
    print(f"  median: {np.median(arr):.0f}")
    print(f"  p10: {np.percentile(arr, 10):.0f}")
    print(f"  p50: {np.percentile(arr, 50):.0f}")
    print(f"  p90: {np.percentile(arr, 90):.0f}")
    print(f"  max: {arr.max()}")
    
    print(f"\n  Top overhang values:")
    for val, cnt in sorted(oh_counter.items(), key=lambda x: -x[1])[:15]:
        print(f"    oh={val}: {cnt:,}")

print(f"\npruning outcomes for overhang reads:")
total_oh = sum(gate_outcomes.values())
for outcome, cnt in sorted(gate_outcomes.items(), key=lambda x: -x[1]):
    print(f"  {outcome}: {cnt:,} ({100*cnt/max(total_oh,1):.1f}%)")

if oh_examples:
    print(f"\n=== EXAMPLES: nRNA wins hard overhang gate due to intronic overhang ===")
    for ex in oh_examples[:10]:
        print(f"\n  {ex['qname']} (truth: {ex['truth_t']}):")
        print(f"    pos={ex['pos']}, cigar={ex['cigar']}")
        print(f"    rl={ex['rl']}, exon_bp={ex['exon_bp']}, mrna_oh={ex['mrna_oh']}, nrna_oh={ex['nrna_oh']}")
        print(f"    read blocks: {ex['blocks']}")
        print(f"    nearby exons: {ex['nearby_exons']}")
