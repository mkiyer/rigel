#!/usr/bin/env python3
"""Prove nRNA hard overhang gate bug — optimized with sorted nRNA spans."""
import sys
sys.path.insert(0, "src")

import pysam
import numpy as np
import pyarrow.feather as pf
from collections import Counter
from bisect import bisect_left

BASE = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v2"
COND = "gdna_none_ss_0.95_nrna_none"
RAW_BAM = f"{BASE}/{COND}/align_minimap2/reads_namesort.bam"
INDEX = f"{BASE}/rigel_index"

t_df = pf.read_table(f"{INDEX}/transcripts.feather").to_pandas()
iv_df = pf.read_table(f"{INDEX}/intervals.feather").to_pandas()
nrna_set = set(t_df.loc[t_df["is_synthetic_nrna"] == True, "t_index"].values)
t_id_to_idx = dict(zip(t_df["t_id"], t_df["t_index"]))

# Build exon lookup
print("Building exon index...")
exon_ivs = iv_df[(iv_df["interval_type"] == 0) & (iv_df["t_index"] >= 0)]
mrna_exon_ivs = exon_ivs[~exon_ivs["t_index"].isin(nrna_set)]
t_exons = {}
for ti, start, end in zip(mrna_exon_ivs["t_index"].values,
                           mrna_exon_ivs["start"].values,
                           mrna_exon_ivs["end"].values):
    ti = int(ti)
    if ti not in t_exons:
        t_exons[ti] = []
    t_exons[ti].append((int(start), int(end)))
for ti in t_exons:
    t_exons[ti].sort()
t_exons = {ti: exons for ti, exons in t_exons.items() if len(exons) > 1}

# Pre-compute: for each transcript, is there a covering nRNA?
# Build nRNA lookup per transcript
t_ref = dict(zip(t_df["t_index"], t_df["ref"]))
t_strand = dict(zip(t_df["t_index"], t_df["strand"]))

# Build: for each mRNA t_index -> covering nRNA span(s)
# nRNA covers mRNA if nrna.start <= mRNA.start and nrna.end >= mRNA.end
print("Building nRNA coverage lookup...")
nrna_df = t_df[t_df["is_synthetic_nrna"] == True]
# For each nRNA, the span is [start, end)
nrna_by_ref_strand = {}
for ref, strand, start, end in zip(
    nrna_df["ref"].values, nrna_df["strand"].values,
    nrna_df["start"].values, nrna_df["end"].values
):
    key = (ref, int(strand))
    if key not in nrna_by_ref_strand:
        nrna_by_ref_strand[key] = []
    nrna_by_ref_strand[key].append((int(start), int(end)))
# Sort by start
for key in nrna_by_ref_strand:
    nrna_by_ref_strand[key].sort()

# For each mRNA transcript, find if there's a covering nRNA
t_has_nrna = set()
for ti in t_exons:
    ref = t_ref.get(ti)
    strand = t_strand.get(ti)
    if ref is None or strand is None:
        continue
    key = (ref, int(strand))
    nrnas = nrna_by_ref_strand.get(key, [])
    exons = t_exons[ti]
    t_start = exons[0][0]
    t_end = exons[-1][1]
    for ns, ne in nrnas:
        if ns <= t_start and ne >= t_end:
            t_has_nrna.add(ti)
            break

print(f"Multi-exon mRNA with covering nRNA: {len(t_has_nrna)} / {len(t_exons)}")

# Scan BAM — only check transcripts with covering nRNA
print("\nScanning BAM...")
n_pairs = 0
n_overhang = 0
oh_vals = []
oh_examples = []
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
        
        colon = qname.find(":")
        if colon < 0:
            continue
        truth_t_id = qname[:colon]
        truth_t_idx = t_id_to_idx.get(truth_t_id, -1)
        if truth_t_idx < 0 or truth_t_idx not in t_exons:
            continue
        if truth_t_idx not in t_has_nrna:
            continue  # Skip transcripts without covering nRNA
        
        exons = t_exons[truth_t_idx]
        ref = t_ref[truth_t_idx]
        
        for r in [r1, r2]:
            if not r.cigartuples or r.reference_name != ref:
                continue
            
            has_nsplice = any(op == 3 for op, _ in r.cigartuples)
            if has_nsplice:
                continue
            
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
            
            exon_bp = 0
            for es, ee in exons:
                for bs, be in blocks:
                    ov = min(be, ee) - max(bs, es)
                    if ov > 0:
                        exon_bp += ov
            
            rl = sum(be - bs for bs, be in blocks)
            mrna_oh = max(rl - exon_bp, 0)
            
            if mrna_oh > 0:
                n_overhang += 1
                oh_vals.append(mrna_oh)
                # nRNA always has oh=0 since its single exon covers everything
                gate_outcomes["nrna_wins"] += 1
                
                if len(oh_examples) < 20:
                    read_start = blocks[0][0]
                    read_end = blocks[-1][1]
                    nearby = [(s,e) for s,e in exons if s<=read_end and e>=read_start]
                    oh_examples.append({
                        "qname": qname,
                        "truth_t": truth_t_id,
                        "pos": f"{r.reference_name}:{r.reference_start}",
                        "cigar": r.cigarstring,
                        "rl": rl,
                        "exon_bp": exon_bp,
                        "mrna_oh": mrna_oh,
                        "blocks": blocks,
                        "nearby_exons": nearby,
                    })
            else:
                gate_outcomes["tie"] += 1
        
        if n_pairs >= 500_000:
            break
        if n_pairs % 100_000 == 0:
            print(f"  {n_pairs:,} pairs, {n_overhang:,} overhang reads")

print(f"\nTotal pairs: {n_pairs:,}")
print(f"Reads with intronic overhang (no N-skip): {n_overhang:,}")
print(f"Rate: {n_overhang/max(n_pairs,1)*100:.2f}% of pairs have at least one overhang read")

if oh_vals:
    arr = np.array(oh_vals)
    print(f"\nmRNA overhang distribution:")
    print(f"  mean: {arr.mean():.1f}")
    print(f"  median: {np.median(arr):.0f}")
    print(f"  p10/p50/p90: {np.percentile(arr,10):.0f}/{np.percentile(arr,50):.0f}/{np.percentile(arr,90):.0f}")
    print(f"  max: {arr.max()}")
    
    oh_counter = Counter(oh_vals)
    print(f"\n  Top overhang values:")
    for val, cnt in sorted(oh_counter.items(), key=lambda x: -x[1])[:15]:
        print(f"    oh={val}: {cnt:,}")

print(f"\npruning outcomes for reads with/without overhang:")
for outcome, cnt in sorted(gate_outcomes.items(), key=lambda x: -x[1]):
    total = sum(gate_outcomes.values())
    print(f"  {outcome}: {cnt:,} ({100*cnt/max(total,1):.1f}%)")

if oh_examples:
    print(f"\n=== SMOKING GUN: nRNA wins hard overhang gate due to intronic overhang ===")
    for ex in oh_examples[:10]:
        print(f"\n  {ex['qname']} (truth: {ex['truth_t']}):")
        print(f"    pos={ex['pos']}, cigar={ex['cigar']}")
        print(f"    rl={ex['rl']}, exon_bp={ex['exon_bp']}")
        print(f"    mrna_oh={ex['mrna_oh']}, nrna_oh=0 → nRNA WINS hard overhang gate")
        print(f"    read blocks: {ex['blocks']}")
        print(f"    overlapping exons: {ex['nearby_exons']}")
