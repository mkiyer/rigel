#!/usr/bin/env python3
"""Prove that nRNA synthetic transcripts win hard overhang gate over mRNA due to overhang.

For fragments where minimap2 doesn't detect the splice junction, the read
extends a few bp into intronic space. For the mRNA transcript, this creates
overhang (oh > 0). For the synthetic nRNA transcript (single exon spanning
the whole locus including introns), oh = 0. The hard overhang gate keeps min-overhang
winners only, so the nRNA wins and the mRNA is discarded.

This script traces actual fragments from the annotated BAM to prove this.
"""
import sys
sys.path.insert(0, "src")

import pysam
import numpy as np
from collections import Counter, defaultdict
import pyarrow.feather as pf

BASE = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v2"
COND = "gdna_none_ss_0.95_nrna_none"
BAM = f"{BASE}/{COND}/align_minimap2/annotated_default.bam"
INDEX = f"{BASE}/rigel_index"

# Load transcript info
t_df = pf.read_table(f"{INDEX}/transcripts.feather").to_pandas()
iv_df = pf.read_table(f"{INDEX}/intervals.feather").to_pandas()

# Identify nRNA synthetic transcripts
nrna_mask = t_df["is_synthetic_nrna"] == True
nrna_t_indices = set(t_df.loc[nrna_mask, "t_index"].values)
print(f"Total transcripts: {len(t_df)}")
print(f"Synthetic nRNA transcripts: {len(nrna_t_indices)}")

# Count interval types for nRNA vs mRNA
nrna_iv = iv_df[iv_df["t_index"].isin(nrna_t_indices)]
mrna_iv = iv_df[(iv_df["t_index"] >= 0) & (~iv_df["t_index"].isin(nrna_t_indices))]
print(f"\nnRNA intervals: {len(nrna_iv)}")
for itype, count in nrna_iv["interval_type"].value_counts().items():
    print(f"  type {itype}: {count}")
print(f"\nmRNA intervals: {len(mrna_iv)}")
for itype, count in mrna_iv["interval_type"].value_counts().items():
    print(f"  type {itype}: {count}")

# nRNA exon interval sizes
nrna_exon = nrna_iv[nrna_iv["interval_type"] == 0]  # EXON
if len(nrna_exon) > 0:
    nrna_exon_sizes = (nrna_exon["end"] - nrna_exon["start"]).values
    print(f"\nnRNA EXON interval sizes:")
    print(f"  mean: {nrna_exon_sizes.mean():.0f}")
    print(f"  median: {np.median(nrna_exon_sizes):.0f}")
    print(f"  min: {nrna_exon_sizes.min()}")
    print(f"  max: {nrna_exon_sizes.max()}")

# For each nRNA transcript, check if it's single-exon
print(f"\nnRNA transcript exon counts:")
nrna_exon_counts = nrna_exon.groupby("t_index").size()
print(f"  all single-exon: {(nrna_exon_counts == 1).all()}")
print(f"  distribution: {nrna_exon_counts.value_counts().to_dict()}")

# Now trace actual fragments from annotated BAM
# The ZP tag tells us the partition (mRNA, intergenic, chimeric)
# The ZL tag gives fragment length
# ZT = transcript match info
print("\n" + "="*60)
print("TRACING ANNOTATED BAM FRAGMENTS")
print("="*60)

# Read annotated BAM and look for fragments classified as mRNA vs intergenic
partition_counts = Counter()
zl_by_partition = defaultdict(list)
oh_examples = []

n_reads = 0
n_r1 = 0
with pysam.AlignmentFile(BAM, "rb") as bam:
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        n_reads += 1
        if not read.is_read2:
            n_r1 += 1
            try:
                zp = read.get_tag("ZP")
                zl = read.get_tag("ZL")
                partition_counts[zp] += 1
                if zl > 0:
                    zl_by_partition[zp].append(zl)
                
                # Collect some examples of mRNA reads with large FL
                if zp == "mRNA" and zl > 500 and len(oh_examples) < 20:
                    zt = read.get_tag("ZT") if read.has_tag("ZT") else "?"
                    oh_examples.append({
                        "qname": read.query_name,
                        "ref": read.reference_name,
                        "pos": read.reference_start,
                        "zl": zl,
                        "zt": zt,
                        "cigar": read.cigarstring,
                    })
            except KeyError:
                partition_counts["no_tags"] += 1
        
        if n_reads >= 10_000_000:
            break

print(f"\nTotal R1 reads: {n_r1:,}")
print(f"\nPartition distribution:")
for p, count in sorted(partition_counts.items(), key=lambda x: -x[1]):
    print(f"  {p}: {count:,} ({100*count/n_r1:.2f}%)")

for p, fls in zl_by_partition.items():
    arr = np.array(fls)
    if len(arr) > 0:
        print(f"\n  {p} FL (n={len(arr):,}):")
        print(f"    mean={arr.mean():.1f}, median={np.median(arr):.1f}")
        print(f"    >500: {(arr>500).sum():,} ({(arr>500).mean()*100:.1f}%)")
        print(f"    >1000: {(arr>1000).sum():,} ({(arr>1000).mean()*100:.1f}%)")

if oh_examples:
    print(f"\nExamples of mRNA reads with FL > 500:")
    for ex in oh_examples[:10]:
        print(f"  {ex['qname']}: ref={ex['ref']}:{ex['pos']}, FL={ex['zl']}, ZT={ex['zt']}, CIGAR={ex['cigar']}")
