#!/usr/bin/env python3
"""Investigate WHY gap SJ correction fails for minimap2 fragments.

For specific examples of inflated fragment lengths, check:
1. Does the SJ gap index contain the intron between the mates?
2. Does the fragment resolve to the correct transcript?
3. Is the containment check (hs >= gs && he <= ge) passing?
"""
import numpy as np
import pandas as pd
import pysam

# Load the rigel index to check SJ data
from rigel.index import TranscriptIndex

INDEX_DIR = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/rigel_index"
idx = TranscriptIndex.load(INDEX_DIR)

# Get SJ data from the index
sj_df = pd.read_feather(f"{INDEX_DIR}/sj.feather")
print(f"Total splice junctions in index: {len(sj_df):,}")
print(f"SJ columns: {list(sj_df.columns)}")
print(f"Sample:\n{sj_df.head()}")

# Focus on ENST00000000412.8 on chr12
# From trace: R1 on exon ending ~8946405, R2 on exon starting ~8949486
# Gap: [8946405, 8949486) — this should contain an intron

t_id = "ENST00000000412.8"
print(f"\n{'='*70}")
print(f"Investigating {t_id} on chr12")
print(f"{'='*70}")

# Find this transcript in the index
t_df = pd.read_feather(f"{INDEX_DIR}/transcripts.feather")
t_row = t_df[t_df["t_id"] == t_id]
if len(t_row) > 0:
    t_row = t_row.iloc[0]
    print(f"  t_index: {t_row['t_index']}")
    print(f"  ref: {t_row['ref']}")
    print(f"  strand: {t_row['strand']}")
    print(f"  start: {t_row['start']}")
    print(f"  end: {t_row['end']}")
    t_index = t_row["t_index"]
else:
    print(f"  NOT FOUND in index!")
    t_index = -1

# Find SJs for this transcript
t_sjs = sj_df[sj_df["t_index"] == t_index]
print(f"\n  Splice junctions for {t_id} (t_index={t_index}):")
for _, sj in t_sjs.iterrows():
    print(f"    ref={sj['ref']} [{sj['start']}, {sj['end']}) "
          f"strand={sj['strand']} size={sj['end']-sj['start']}")

# Check: is there an SJ in the gap [8946405, 8949486)?
gap_start = 8946405
gap_end = 8949486
print(f"\n  Gap region: [{gap_start}, {gap_end})")
matching_sjs = sj_df[
    (sj_df["ref"] == "chr12") &
    (sj_df["start"] >= gap_start) &
    (sj_df["end"] <= gap_end)
]
print(f"  SJs CONTAINED within gap (start >= {gap_start} AND end <= {gap_end}):")
for _, sj in matching_sjs.iterrows():
    is_target = sj["t_index"] == t_index
    print(f"    ref={sj['ref']} [{sj['start']}, {sj['end']}) "
          f"t_index={sj['t_index']} {'*** THIS TRANSCRIPT ***' if is_target else ''}")

# Also check for SJs that OVERLAP the gap but aren't fully contained
overlapping_sjs = sj_df[
    (sj_df["ref"] == "chr12") &
    (sj_df["start"] < gap_end) &
    (sj_df["end"] > gap_start)
]
print(f"\n  SJs OVERLAPPING gap (not necessarily contained):")
for _, sj in overlapping_sjs.iterrows():
    is_target = sj["t_index"] == t_index
    contained = sj["start"] >= gap_start and sj["end"] <= gap_end
    print(f"    ref={sj['ref']} [{sj['start']}, {sj['end']}) "
          f"t_index={sj['t_index']} contained={contained} "
          f"{'*** THIS TRANSCRIPT ***' if is_target else ''}")

# Check the genomic intervals for this region
iv_df = pd.read_feather(f"{INDEX_DIR}/intervals.feather")
print(f"\n  Genomic intervals covering gap region chr12:{gap_start}-{gap_end}:")
iv_overlap = iv_df[
    (iv_df["ref"] == "chr12") &
    (iv_df["start"] < gap_end) &
    (iv_df["end"] > gap_start)
]
for _, iv in iv_overlap.head(10).iterrows():
    print(f"    [{iv['start']}, {iv['end']}) type={iv.get('type', '?')}")

# Now check: what transcripts does the interval resolution return for this fragment?
# The fragment has exon blocks at [8946235, 8946405] and [8949486, 8949638]
# Let's see what the resolver gives us
print(f"\n  Checking resolver for exon blocks [8946235, 8946405] and [8949486, 8949638]:")

# Use the resolver to find matching transcripts
results1 = idx.query_interval("chr12", 8946235, 8946405)
results2 = idx.query_interval("chr12", 8949486, 8949638)

t_set_1 = set()
for _, _, _, tset in results1:
    t_set_1.update(tset)
t_set_2 = set()
for _, _, _, tset in results2:
    t_set_2.update(tset)

# Intersection
t_both = t_set_1 & t_set_2
print(f"  Transcripts matching block 1: {len(t_set_1)}")
print(f"  Transcripts matching block 2: {len(t_set_2)}")
print(f"  Transcripts matching BOTH blocks: {len(t_both)}")
print(f"  Target t_index {t_index} in block 1: {t_index in t_set_1}")
print(f"  Target t_index {t_index} in block 2: {t_index in t_set_2}")
print(f"  Target t_index {t_index} in intersection: {t_index in t_both}")

# Check if ANY of the overlapping SJs belong to transcripts in t_both
if len(overlapping_sjs) > 0:
    print(f"\n  SJ gap correction check:")
    for _, sj in overlapping_sjs.iterrows():
        in_t_set = sj["t_index"] in t_both
        contained = sj["start"] >= gap_start and sj["end"] <= gap_end
        would_correct = in_t_set and contained
        correction = sj["end"] - sj["start"] if would_correct else 0
        print(f"    SJ [{sj['start']}, {sj['end']}) t_index={sj['t_index']}: "
              f"in_t_set={in_t_set} contained={contained} "
              f"would_correct={would_correct} correction={correction}")

# Check another case — the mode=1000 suggests many fragments at exactly 1000
# which is the frag_max. Let's check fragments clamped to 1000
print(f"\n{'='*70}")
print(f"Checking for fragments clamped at 1000 (frag_max)")
print(f"{'='*70}")

# In the benchmark log, SPLICED_ANNOT mode=1000. This means the FL training
# is seeing many fragment lengths exactly at 1000. Let me check if there's
# clamping in the code.
print("  The FL model mode=1000 suggests that many spliced-annotated unique-mapper ")
print("  fragments have fragment_length computed as 1000 (the frag_max). This is ")
print("  likely because the FL is clamped to frag_max, and many uncorrected pairs ")
print("  have genomic spans > 1000 that get clamped.")

print("\n" + "=" * 70)
print("END OF INVESTIGATION")
print("=" * 70)
