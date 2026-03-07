#!/usr/bin/env python
"""Check whether locus structure (transcript-to-locus assignment) differs between runs."""
import pyarrow.feather as feather
import numpy as np

base = "/Users/mkiyer/Downloads/rigel_runs/baseline_output"
orig = "/Users/mkiyer/Downloads/rigel_runs/orig_run1"

df1 = feather.read_table(f"{base}/quant.feather").to_pandas()
df2 = feather.read_table(f"{orig}/quant.feather").to_pandas()

# Check transcript ordering is the same
assert (df1["transcript_id"].values == df2["transcript_id"].values).all(), "transcript order differs!"
print(f"Total transcripts: {len(df1)}")

# Check locus_id assignment
lid1 = df1["locus_id"].values
lid2 = df2["locus_id"].values

# Direct comparison (locus IDs may be arbitrary, but should be consistent)
same = (lid1 == lid2).sum()
diff = (lid1 != lid2).sum()
print(f"\nDirect locus_id comparison: {same} same, {diff} differ")

# But IDs could just be renumbered. Check structural equivalence:
# Two transcripts are in the same locus iff they share a locus_id
# Build equivalence: for each pair of adjacent transcripts, check if co-locus assignment is the same

# Simpler: group transcripts by locus_id and compare the groupings
from collections import defaultdict

def get_locus_groups(lids, tids):
    groups = defaultdict(set)
    for tid, lid in zip(tids, lids):
        groups[lid].add(tid)
    return set(frozenset(g) for g in groups.values())

groups1 = get_locus_groups(lid1, df1["transcript_id"].values)
groups2 = get_locus_groups(lid2, df2["transcript_id"].values)

print(f"\nNumber of loci: {len(groups1)} vs {len(groups2)}")

# Compare structurally
only_in_1 = groups1 - groups2
only_in_2 = groups2 - groups1
shared = groups1 & groups2

print(f"Structurally identical loci: {len(shared)}")
print(f"Loci only in run1: {len(only_in_1)}")
print(f"Loci only in run2: {len(only_in_2)}")

if only_in_1:
    # Show some examples
    examples = sorted(only_in_1, key=len, reverse=True)[:5]
    print(f"\nLargest loci that differ (run1 version):")
    for g in examples:
        print(f"  {len(g)} transcripts: {sorted(g)[:5]}...")
        # Find where these transcripts ended up in run2
        sample_tid = next(iter(g))
        idx = np.where(df2["transcript_id"].values == sample_tid)[0][0]
        lid2_val = lid2[idx]
        group2_tids = set(df2["transcript_id"].values[lid2 == lid2_val])
        print(f"    In run2, {sample_tid} is in locus {lid2_val} with {len(group2_tids)} transcripts")
        extra = group2_tids - g
        missing = g - group2_tids
        if extra:
            print(f"    Run2 adds: {sorted(extra)[:5]}")
        if missing:
            print(f"    Run2 loses: {sorted(missing)[:5]}")

# Also check loci.feather structural columns
ldf1 = feather.read_table(f"{base}/loci.feather").to_pandas()
ldf2 = feather.read_table(f"{orig}/loci.feather").to_pandas()

print(f"\n=== loci.feather ===")
print(f"n_transcripts: differ = {(ldf1['n_transcripts'].values != ldf2['n_transcripts'].values).sum()}/{len(ldf1)}")
print(f"n_genes: differ = {(ldf1['n_genes'].values != ldf2['n_genes'].values).sum()}/{len(ldf1)}")
print(f"n_em_fragments: differ = {(ldf1['n_em_fragments'].values != ldf2['n_em_fragments'].values).sum()}/{len(ldf1)}")

# Check if mrna_unambig differs (pre-EM, should be deterministic)
a = df1["mrna_unambig"].values
b = df2["mrna_unambig"].values
print(f"\nmrna_unambig: bit-exact = {np.array_equal(a, b, equal_nan=True)}")
ndiff = np.sum(a != b)
if ndiff > 0:
    diff = np.abs(a - b)
    print(f"  n_diff = {ndiff}/{len(a)}")
    maxidx = np.argmax(diff)
    print(f"  max_abs_diff = {diff[maxidx]:.6e}")
    print(f"  at transcript {df1['transcript_id'].iloc[maxidx]}, values: {a[maxidx]} vs {b[maxidx]}")
