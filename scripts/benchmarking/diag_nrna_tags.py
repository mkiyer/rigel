#!/usr/bin/env python3
"""Quick check: what ZT and ZG tags are on nRNA-assigned reads."""
import pysam

BAM = (
    "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/"
    "gdna_none_nrna_none_ss_1.00/align_oracle/annotated_default.sorted.bam"
)

nrna_count = 0
sample = []
with pysam.AlignmentFile(BAM, "rb") as bam:
    for read in bam.fetch(until_eof=True):
        if read.is_read2:
            continue
        pool = read.get_tag("ZP")
        if pool == "nRNA":
            nrna_count += 1
            if len(sample) < 20:
                zt = read.get_tag("ZT")
                zg = read.get_tag("ZG")
                sample.append((read.query_name[:50], zt, zg, read.get_tag("ZW")))

print(f"Total nRNA R1 reads: {nrna_count}")
print()
print("Sample nRNA reads:")
for qn, zt, zg, zw in sample:
    print(f"  qname={qn:50s}  ZT={zt:30s}  ZG={zg:30s}  ZW={zw:.4f}")

# Count by gene
from collections import Counter
gene_counts = Counter()
with pysam.AlignmentFile(BAM, "rb") as bam:
    for read in bam.fetch(until_eof=True):
        if read.is_read2:
            continue
        if read.get_tag("ZP") == "nRNA":
            gene_counts[read.get_tag("ZG")] += 1

print()
print("nRNA by gene:")
for g, c in gene_counts.most_common():
    print(f"  {g:30s}  {c}")

# Also check: what are the ZT values for nRNA?
tx_counts = Counter()
with pysam.AlignmentFile(BAM, "rb") as bam:
    for read in bam.fetch(until_eof=True):
        if read.is_read2:
            continue
        if read.get_tag("ZP") == "nRNA":
            tx_counts[read.get_tag("ZT")] += 1

print()
print("nRNA by transcript (top 10):")
for t, c in tx_counts.most_common(10):
    print(f"  {t:30s}  {c}")
