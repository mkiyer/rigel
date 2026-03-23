#!/usr/bin/env python3
"""Analyze PE secondary pairing to understand chimera classification.

For each fragment (qname) in the minimap2 BAM, group all records,
identify secondaries, and check if cross-chromosome R1/R2 secondaries exist.
"""
import pysam
from collections import defaultdict, Counter

bam = '/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file/sim_minimap2.bam'

# Group all records by qname
qname_records = defaultdict(list)
with pysam.AlignmentFile(bam, 'rb') as f:
    for r in f:
        qname_records[r.query_name].append(r)

n_fragments = len(qname_records)
print(f"Total fragments (qnames): {n_fragments}")

# Analyze each fragment
cross_chrom_count = 0
singleton_sec_count = 0
same_chrom_sec_count = 0
fragments_with_secs = 0
fragments_with_cross_chrom = 0

# For detailed analysis, track per-fragment
per_frag_sec_chroms = Counter()

for qname, records in qname_records.items():
    # Separate primary and secondary
    pri_r1 = [r for r in records if r.flag & 0x40 and not r.flag & 0x100]
    pri_r2 = [r for r in records if r.flag & 0x80 and not r.flag & 0x100]
    sec_r1 = [r for r in records if r.flag & 0x40 and r.flag & 0x100]
    sec_r2 = [r for r in records if r.flag & 0x80 and r.flag & 0x100]

    if not sec_r1 and not sec_r2:
        continue

    fragments_with_secs += 1

    sec_r1_chroms = set(r.reference_name for r in sec_r1)
    sec_r2_chroms = set(r.reference_name for r in sec_r2)

    # Check for cross-chromosome secondary pairs
    # The CROSS-PAIR step in pair_multimapper_reads would pair unmatched R1/R2
    # across chromosomes
    common_chroms = sec_r1_chroms & sec_r2_chroms
    only_r1_chroms = sec_r1_chroms - sec_r2_chroms
    only_r2_chroms = sec_r2_chroms - sec_r1_chroms

    has_cross = False
    if only_r1_chroms and only_r2_chroms:
        # R1 on some chrom, R2 on different chrom → potential cross-pair
        has_cross = True
        fragments_with_cross_chrom += 1
    elif sec_r1 and not sec_r2:
        singleton_sec_count += 1
    elif sec_r2 and not sec_r1:
        singleton_sec_count += 1
    elif common_chroms:
        same_chrom_sec_count += 1

    if not common_chroms and (only_r1_chroms or only_r2_chroms):
        cross_chrom_count += 1

    per_frag_sec_chroms[f"r1={'|'.join(sorted(sec_r1_chroms))},r2={'|'.join(sorted(sec_r2_chroms))}"] += 1

print(f"\nFragments with secondaries: {fragments_with_secs}")
print(f"Fragments with potential cross-chrom pairing: {fragments_with_cross_chrom}")
print(f"Fragments with only R1 or only R2 secs: {singleton_sec_count}")
print(f"Fragments with same-chrom secs: {same_chrom_sec_count}")
print(f"\nTop secondary chrom patterns:")
for pattern, count in per_frag_sec_chroms.most_common(20):
    print(f"  {count:>6d}  {pattern}")

# Detailed: for fragments with cross-chrom, show structure
print(f"\n--- Sample cross-chrom fragments ---")
shown = 0
for qname, records in qname_records.items():
    sec_r1 = [r for r in records if r.flag & 0x40 and r.flag & 0x100]
    sec_r2 = [r for r in records if r.flag & 0x80 and r.flag & 0x100]
    if not sec_r1 or not sec_r2:
        continue
    sec_r1_chroms = set(r.reference_name for r in sec_r1)
    sec_r2_chroms = set(r.reference_name for r in sec_r2)
    only_r1 = sec_r1_chroms - sec_r2_chroms
    only_r2 = sec_r2_chroms - sec_r1_chroms
    if only_r1 and only_r2:
        print(f"\n  {qname}")
        for r in sorted(records, key=lambda x: (x.flag & 0x100, x.flag & 0x80)):
            tag = 'R1' if r.flag & 0x40 else 'R2'
            sec = 'sec' if r.flag & 0x100 else 'pri'
            nm = r.get_tag('NM') if r.has_tag('NM') else -1
            print(f"    {tag} {sec} {r.reference_name}:{r.reference_start} NM={nm} cigar={r.cigarstring}")
        shown += 1
        if shown >= 3:
            break
