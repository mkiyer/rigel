"""Trace individual GAPDH fragments through oracle vs minimap2 resolution.

For a sample of fragments, show:
- Oracle: which transcripts/how many candidates
- Minimap2: which transcripts/how many candidates, which are on other chroms
- NM tag values across hits
"""
import sys
sys.path.insert(0, "src")

import pysam
import pandas as pd
from collections import defaultdict

base = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
oracle_bam = f"{base}/sim_oracle_nsorted.bam"
mm2_bam = f"{base}/sim_minimap2.bam"

# Parse all reads from both BAMs, grouped by qname
def parse_bam_by_qname(bam_path, max_qnames=None):
    """Parse BAM, return dict[qname] -> list of (chrom, pos, cigar, flag, mapq, NM, is_sec)"""
    groups = defaultdict(list)
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam:
            qname = read.query_name
            if max_qnames and len(groups) >= max_qnames and qname not in groups:
                continue
            nm = read.get_tag("NM") if read.has_tag("NM") else -1
            groups[qname].append({
                "chrom": read.reference_name,
                "pos": read.reference_start,
                "cigar": read.cigarstring,
                "flag": read.flag,
                "mapq": read.mapping_quality,
                "nm": nm,
                "is_r1": bool(read.flag & 0x40),
                "is_r2": bool(read.flag & 0x80),
                "is_sec": bool(read.flag & 0x100),
                "is_sup": bool(read.flag & 0x800),
                "is_rev": bool(read.flag & 0x10),
            })
    return groups

print("Parsing oracle BAM...")
oracle_groups = parse_bam_by_qname(oracle_bam)
print(f"  {len(oracle_groups)} qnames")

print("Parsing minimap2 BAM...")
mm2_groups = parse_bam_by_qname(mm2_bam)
print(f"  {len(mm2_groups)} qnames")

# Classify minimap2 fragments
unique_frags = []
multi_frags = []
chimeric_frags = []

for qname, reads in mm2_groups.items():
    primary = [r for r in reads if not r["is_sec"] and not r["is_sup"]]
    secondary = [r for r in reads if r["is_sec"]]
    
    # Check if primary pair maps to different chroms
    pri_chroms = set(r["chrom"] for r in primary)
    if len(pri_chroms) > 1:
        chimeric_frags.append(qname)
    elif len(secondary) == 0:
        unique_frags.append(qname)
    else:
        multi_frags.append(qname)

print(f"\nMinimap2 fragment classification:")
print(f"  Unique (no secondaries): {len(unique_frags)}")
print(f"  Multimapping (has secondaries): {len(multi_frags)}")
print(f"  Chimeric (primary on diff chroms): {len(chimeric_frags)}")

# Analyze multimappers: where do the secondary alignments go?
print(f"\n{'='*80}")
print(f"Multimapper Analysis")
print(f"{'='*80}")

# Count how many secondaries per fragment
n_sec_hist = defaultdict(int)
sec_chrom_counts = defaultdict(int)

for qname in multi_frags:
    reads = mm2_groups[qname]
    secondary = [r for r in reads if r["is_sec"]]
    n_sec = len(secondary)
    n_sec_hist[n_sec] += 1
    for r in secondary:
        sec_chrom_counts[r["chrom"]] += 1

print("\nSecondary alignment count distribution:")
for n, cnt in sorted(n_sec_hist.items()):
    print(f"  {n} secondaries: {cnt} fragments")

print(f"\nTop chromosomes in secondary alignments:")
for chrom, cnt in sorted(sec_chrom_counts.items(), key=lambda x: -x[1])[:20]:
    print(f"  {chrom}: {cnt}")

# Detailed trace of a few representative fragments
print(f"\n{'='*80}")
print(f"Detailed Fragment Traces (first 10 multimappers)")
print(f"{'='*80}")

for i, qname in enumerate(multi_frags[:10]):
    # Extract truth transcript from qname
    parts = qname.split(":")
    truth_tx = parts[0]
    
    oracle_reads = oracle_groups.get(qname, [])
    mm2_reads = mm2_groups[qname]
    
    primary_mm2 = [r for r in mm2_reads if not r["is_sec"] and not r["is_sup"]]
    secondary_mm2 = [r for r in mm2_reads if r["is_sec"]]
    
    print(f"\n--- Fragment {i+1}: {qname} ---")
    print(f"  Truth transcript: {truth_tx}")
    print(f"  Oracle: {len(oracle_reads)} reads")
    for r in oracle_reads:
        tag = "R1" if r["is_r1"] else "R2"
        strand = "-" if r["is_rev"] else "+"
        print(f"    {tag} {r['chrom']}:{r['pos']} {strand} CIGAR={r['cigar']} MQ={r['mapq']} NM={r['nm']}")
    
    print(f"  MM2 primary: {len(primary_mm2)} reads")
    for r in primary_mm2:
        tag = "R1" if r["is_r1"] else "R2"
        strand = "-" if r["is_rev"] else "+"
        print(f"    {tag} {r['chrom']}:{r['pos']} {strand} CIGAR={r['cigar']} MQ={r['mapq']} NM={r['nm']}")
    
    print(f"  MM2 secondary: {len(secondary_mm2)} reads")
    for r in secondary_mm2:
        tag = "R1" if r["is_r1"] else "R2"
        strand = "-" if r["is_rev"] else "+"
        print(f"    {tag} {r['chrom']}:{r['pos']} {strand} CIGAR={r['cigar']} MQ={r['mapq']} NM={r['nm']}")

# NM distribution for primary vs secondary
print(f"\n{'='*80}")
print(f"NM Tag Analysis (Edit Distance)")
print(f"{'='*80}")

primary_nm = []
secondary_nm = []
for qname in multi_frags:
    reads = mm2_groups[qname]
    for r in reads:
        if r["is_sec"]:
            secondary_nm.append(r["nm"])
        elif not r["is_sup"]:
            primary_nm.append(r["nm"])

print(f"\nPrimary alignments NM distribution:")
nm_hist = defaultdict(int)
for nm in primary_nm:
    nm_hist[nm] += 1
for nm, cnt in sorted(nm_hist.items()):
    print(f"  NM={nm}: {cnt} ({100*cnt/len(primary_nm):.1f}%)")

print(f"\nSecondary alignments NM distribution:")
nm_hist2 = defaultdict(int)
for nm in secondary_nm:
    nm_hist2[nm] += 1
for nm, cnt in sorted(nm_hist2.items()):
    print(f"  NM={nm}: {cnt} ({100*cnt/len(secondary_nm):.1f}%)")

# Key question: how many multimappers have NM=0 primary AND NM>0 secondary?
# (These are the ones where edit distance COULD help)
print(f"\n{'='*80}")
print(f"Edit Distance Discriminability")
print(f"{'='*80}")

n_pri_nm0_sec_nm_gt0 = 0
n_pri_nm0_sec_nm0 = 0
n_pri_nm_gt0 = 0

for qname in multi_frags:
    reads = mm2_groups[qname]
    primary = [r for r in reads if not r["is_sec"] and not r["is_sup"]]
    secondary = [r for r in reads if r["is_sec"]]
    
    pri_nm = [r["nm"] for r in primary]
    sec_nm = [r["nm"] for r in secondary]
    
    if all(nm == 0 for nm in pri_nm):
        if all(nm == 0 for nm in sec_nm):
            n_pri_nm0_sec_nm0 += 1
        elif any(nm > 0 for nm in sec_nm):
            n_pri_nm0_sec_nm_gt0 += 1
    else:
        n_pri_nm_gt0 += 1

print(f"\nMultimappers where primary NM=0, ALL secondary NM=0: {n_pri_nm0_sec_nm0} (indistinguishable)")
print(f"Multimappers where primary NM=0, SOME secondary NM>0: {n_pri_nm0_sec_nm_gt0} (edit distance helps)")
print(f"Multimappers where primary NM>0: {n_pri_nm_gt0}")
print(f"Total multimappers: {len(multi_frags)}")
