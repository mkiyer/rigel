#!/usr/bin/env python3
"""Compare FL distributions between pre-fix and post-fix annotated BAMs.

The annotated BAM stores per-transcript FL in tags. We extract these
to see if the overlap fix improved the FL distribution.
"""
import pysam
import numpy as np
from collections import Counter

BASE_V1 = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine"
BASE_V2 = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v2"
COND = "gdna_none_ss_0.95_nrna_none"

def extract_fl_from_bam(bam_path, n_max=500000):
    """Extract fragment lengths from annotated BAM.
    
    Reads the 'fl' tag (integer fragment length) and the
    'Zf' tag (comma-sep per-transcript FL).
    """
    fl_list = []  # per-fragment FLs
    read1_tags = {}  # qname -> tags
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        count = 0
        for read in bam:
            if read.is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            
            qname = read.query_name
            tags = dict(read.get_tags())
            
            if qname not in read1_tags:
                read1_tags[qname] = tags
            else:
                # Have both reads, combine
                # fl tag = fragment length used by model
                for tdict in [read1_tags.pop(qname), tags]:
                    if 'fl' in tdict:
                        fl_list.append(int(tdict['fl']))
                        break
                    elif 'Zf' in tdict:
                        # Comma-sep per-transcript FLs
                        fls = [int(x) for x in str(tdict['Zf']).split(',') if x]
                        fl_list.extend(fls)
                        break
                
                count += 1
                if count >= n_max:
                    break
    
    return np.array(fl_list)


# First check which tags are available in the annotated BAM
for label, base in [("PRE-FIX", BASE_V1), ("POST-FIX", BASE_V2)]:
    bam_path = f"{base}/{COND}/align_minimap2/annotated_default.bam"
    print(f"\n=== {label}: {bam_path} ===")
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for i, read in enumerate(bam):
            if read.is_unmapped or read.is_secondary:
                continue
            tags = read.get_tags()
            print(f"  Read {read.query_name}: tags = {[t[0] for t in tags]}")
            if len(tags) > 0:
                for t, v in tags:
                    if isinstance(v, str) and len(v) > 50:
                        print(f"    {t} = {v[:100]}...")
                    else:
                        print(f"    {t} = {v}")
            if i >= 5:
                break

# Extract unique mapper FL from raw minimap2 BAM 
print("\n\n--- Extracting unique mapper FL from raw minimap2 BAM ---")
raw_bam = f"{BASE_V1}/{COND}/align_minimap2/reads_namesort.bam"
print(f"BAM: {raw_bam}")

# Count unique vs multi reads
nh_counts = Counter()
with pysam.AlignmentFile(raw_bam, "rb") as bam:
    n = 0
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        try:
            nh = read.get_tag("NH")
        except KeyError:
            nh = 1
        nh_counts[nh] += 1
        n += 1
        if n >= 1000000:
            break
print(f"  NH distribution (first 1M primary reads): {dict(sorted(nh_counts.items())[:10])}")
