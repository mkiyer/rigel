#!/usr/bin/env python3
"""Trace FL model and gap correction quality after the overlap fix."""
import sys
sys.path.insert(0, "src")

import pysam
import numpy as np
from collections import defaultdict, Counter
from rigel.index import TranscriptIndex

BASE_V1 = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine"
BASE_V2 = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v2"
COND = "gdna_none_ss_0.95_nrna_none"

# Load index (same for both)
idx = TranscriptIndex.load(f"{BASE_V2}/rigel_index")
ctx = idx.resolve_context

# Load annotated BAMs from both runs
bam_v1 = f"{BASE_V1}/{COND}/align_minimap2/annotated.bam"
bam_v2 = f"{BASE_V2}/{COND}/align_minimap2/annotated.bam"

def extract_fl_from_annotated_bam(bam_path, n_max=200000):
    """Extract fragment lengths from the annotated BAM's fl tag."""
    fl_values = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        seen = set()
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            if qname in seen:
                continue
            seen.add(qname)
            # Try to get FL from tags
            try:
                fl = read.get_tag("fl")
                if fl > 0:
                    fl_values.append(fl)
            except KeyError:
                pass
            if len(fl_values) >= n_max:
                break
    return np.array(fl_values)

def compute_fl_directly(bam_path, idx, n_max=200000):
    """Compute fragment lengths using the resolve context (current code)."""
    from rigel.native import ResolverScratch, ExonBlock, IntronBlock
    ctx = idx.resolve_context
    
    fl_all = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        seen = {}
        count = 0
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            if qname not in seen:
                seen[qname] = read
                continue
            
            r1 = seen.pop(qname)
            r2 = read
            
            # Get NH tag to check if unique
            try:
                nh = r1.get_tag("NH")
            except KeyError:
                nh = 1
            
            if nh != 1:
                continue
                
            count += 1
            if count >= n_max:
                break
    
    return count

# First, let's look at FL from the BAM tags
print("Checking annotated BAM for FL tags...")
for label, bam_path in [("PRE-FIX", bam_v1), ("POST-FIX", bam_v2)]:
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            tags_found = set()
            for i, read in enumerate(bam):
                if read.is_unmapped:
                    continue
                for tag, val in read.get_tags():
                    tags_found.add(tag)
                if i > 100:
                    break
            print(f"  {label}: tags = {sorted(tags_found)}")
    except Exception as e:
        print(f"  {label}: error - {e}")

# Instead, let's use the minimap2 BAM directly and compute FL
print("\n--- Computing FL from minimap2 BAM with current code ---")
mm2_bam = f"{BASE_V1}/{COND}/align_minimap2/sorted_by_name.bam"
print(f"BAM: {mm2_bam}")

# Actually let's just compute FL stats from the scan phase output 
# Look for scan stats or model training data
import os, json, glob

for label, base in [("PRE-FIX", BASE_V1), ("POST-FIX", BASE_V2)]:
    print(f"\n=== {label} ===")
    rigel_dir = glob.glob(f"{base}/{COND}/*rigel*") + glob.glob(f"{base}/{COND}/align_*/")
    for d in rigel_dir:
        for f in glob.glob(f"{d}/**", recursive=True):
            if f.endswith('.json') or 'stats' in os.path.basename(f).lower() or 'model' in os.path.basename(f).lower():
                print(f"  Found: {f}")

# Actually, the quant output should have the FL model. Let me search more broadly
print("\n--- Searching for rigel output files ---")
for label, base in [("PRE-FIX", BASE_V1), ("POST-FIX", BASE_V2)]:
    print(f"\n{label}:")
    for f in glob.glob(f"{base}/{COND}/**/*", recursive=True):
        bn = os.path.basename(f)
        if any(x in bn.lower() for x in ['frag', 'model', 'stats', 'quant', '.json']):
            sz = os.path.getsize(f) if os.path.isfile(f) else 0
            print(f"  {os.path.relpath(f, base)} ({sz} bytes)")
