#!/usr/bin/env python3
"""Check intergenic reads in annotated BAM to understand classification."""
import pysam
import numpy as np
from collections import Counter

BAM = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v2/gdna_none_ss_0.95_nrna_none/align_minimap2/annotated_default.bam"

stats = []
zn_vals = []
n = 0
with pysam.AlignmentFile(BAM, "rb") as bam:
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.is_read2:
            continue
        n += 1
        try:
            zp = read.get_tag("ZP")
        except KeyError:
            continue

        if zp == "intergenic":
            try:
                zn = read.get_tag("ZN")
                zn_vals.append(zn)
            except KeyError:
                pass
            if len(stats) < 10:
                tags = {t: v for t, v in read.get_tags()}
                stats.append((read.query_name, {
                    "ZP": tags.get("ZP"),
                    "ZW": tags.get("ZW"),
                    "ZN": tags.get("ZN"),
                    "ZH": tags.get("ZH"),
                    "ZS": tags.get("ZS"),
                    "ZL": tags.get("ZL"),
                    "ZT": tags.get("ZT"),
                    "ZG": tags.get("ZG"),
                    "ZC": tags.get("ZC"),
                    "ref": read.reference_name,
                    "pos": read.reference_start,
                    "cigar": read.cigarstring,
                }))

        if n >= 2_000_000:
            break

print("Intergenic read examples:")
for qname, tags in stats[:5]:
    print(f"  {qname}:")
    for k, v in tags.items():
        print(f"    {k} = {v}")

if zn_vals:
    arr = np.array(zn_vals)
    print(f"\nZN (n_candidates) for intergenic reads:")
    print(f"  n={len(arr):,}")
    cn = Counter(arr)
    for val, cnt in sorted(cn.items())[:10]:
        print(f"  ZN={val}: {cnt:,}")
