#!/usr/bin/env python3
"""Diagnose nRNA leakage in annotated BAM from chr8 PVT1/MYC benchmark."""
import pysam
from collections import Counter, defaultdict

BAM = (
    "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/"
    "gdna_none_nrna_none_ss_1.00/align_oracle/annotated_default.sorted.bam"
)

pool_counts = Counter()
nrna_by_gene = Counter()
nrna_by_tx = Counter()
nrna_by_class = Counter()
nrna_by_splice = Counter()
nrna_examples = []  # first 20 nRNA fragments

seen_names = set()
n_records = 0
n_frags = 0

with pysam.AlignmentFile(BAM, "rb") as bam:
    for read in bam.fetch(until_eof=True):
        n_records += 1
        qname = read.query_name

        # Only process R1 to count fragments once
        if read.is_read2:
            continue

        n_frags += 1
        pool = read.get_tag("ZP")
        pool_counts[pool] += 1

        if pool == "nRNA":
            gene = read.get_tag("ZG")
            tx = read.get_tag("ZT")
            fclass = read.get_tag("ZC")
            splice = read.get_tag("ZS")
            posterior = read.get_tag("ZW")
            n_cand = read.get_tag("ZN")

            nrna_by_gene[gene] += 1
            nrna_by_tx[tx] += 1
            nrna_by_class[fclass] += 1
            nrna_by_splice[splice] += 1

            if len(nrna_examples) < 20:
                nrna_examples.append({
                    "qname": qname,
                    "gene": gene,
                    "tx": tx,
                    "class": fclass,
                    "splice": splice,
                    "posterior": posterior,
                    "n_cand": n_cand,
                    "pos": read.reference_start,
                    "cigar": read.cigarstring,
                    "is_primary": read.get_tag("ZH"),
                })

print(f"Total records: {n_records:,}")
print(f"Total fragments (R1 only): {n_frags:,}")
print()

print("=== POOL DISTRIBUTION ===")
for pool, cnt in pool_counts.most_common():
    print(f"  {pool:15s} {cnt:7,d}  ({100*cnt/n_frags:.2f}%)")

print()
print(f"=== nRNA FRAGMENTS BY GENE (top 15) ===")
for gene, cnt in nrna_by_gene.most_common(15):
    print(f"  {gene:30s} {cnt:6,d}")

print()
print(f"=== nRNA FRAGMENTS BY TRANSCRIPT (top 15) ===")
for tx, cnt in nrna_by_tx.most_common(15):
    print(f"  {tx:30s} {cnt:6,d}")

print()
print(f"=== nRNA FRAGMENTS BY FRAGMENT CLASS ===")
for fclass, cnt in nrna_by_class.most_common():
    print(f"  {fclass:20s} {cnt:6,d}")

print()
print(f"=== nRNA FRAGMENTS BY SPLICE TYPE ===")
for splice, cnt in nrna_by_splice.most_common():
    print(f"  {splice:20s} {cnt:6,d}")

print()
print("=== FIRST 20 nRNA FRAGMENT EXAMPLES ===")
for ex in nrna_examples:
    print(
        f"  {ex['qname']:30s}  gene={ex['gene']:25s}  tx={ex['tx']:25s}  "
        f"class={ex['class']:15s}  splice={ex['splice']:15s}  "
        f"post={ex['posterior']:.4f}  n_cand={ex['n_cand']}  "
        f"pos={ex['pos']}  cigar={ex['cigar']}  primary={ex['is_primary']}"
    )
