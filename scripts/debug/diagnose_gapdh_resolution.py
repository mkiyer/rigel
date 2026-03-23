"""Diagnose why GAPDH reads get ZN=0 (intergenic) despite overlapping exons."""
import sys
sys.path.insert(0, "src")

import pandas as pd
from rigel.index import TranscriptIndex

idx_dir = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"

# 1) Load the full index (builds resolver internally)
print("=== Loading rigel index ===")
tx_index = TranscriptIndex.load(idx_dir)
resolver = tx_index.resolver
print(f"  Transcripts: {len(tx_index.t_to_g_arr)}")

# Get ref-to-id mapping
ref_to_id = resolver.get_ref_to_id()
chr12_id = ref_to_id.get("chr12", -1)
print(f"  chr12 internal ref_id = {chr12_id}")

# 2) Manually resolve a GAPDH read
# From the fragment TSV: chr12:6534806, CIGAR 55M1632N95M 
# Exon blocks: [6534806, 6534861) and [6536493, 6536588)
# Intron: [6534861, 6536493)

STRAND_NONE = 0
STRAND_POS = 1

print("\n=== Spliced GAPDH Read ===")
print("  chr12:6534806, CIGAR=55M1632N95M")
print("  Exon1: [6534806, 6534861)")
print("  Exon2: [6536493, 6536588)")
print("  Intron: [6534861, 6536493)")

result = resolver.resolve(
    exon_ref_ids=[chr12_id, chr12_id],
    exon_starts=[6534806, 6536493],
    exon_ends=[6534861, 6536588],
    exon_strands=[STRAND_NONE, STRAND_NONE],
    intron_ref_ids=[chr12_id],
    intron_starts=[6534861],
    intron_ends=[6536493],
    intron_strands=[STRAND_POS],
    genomic_footprint=6536588 - 6534806,
)

if result is None:
    print("  RESULT: None (intergenic) — BUG CONFIRMED!")
else:
    print(f"  RESULT: resolved with {len(result)} elements")
    for i, v in enumerate(result):
        print(f"    [{i}] = {v}")

# 3) Unspliced read in GAPDH exon
print("\n=== Unspliced GAPDH Read ===")
print("  chr12:6534809-6534959 (150bp in exon 2)")
result2 = resolver.resolve(
    exon_ref_ids=[chr12_id],
    exon_starts=[6534809],
    exon_ends=[6534959],
    exon_strands=[STRAND_NONE],
    intron_ref_ids=[],
    intron_starts=[],
    intron_ends=[],
    intron_strands=[],
    genomic_footprint=150,
)

if result2 is None:
    print("  RESULT: None (intergenic) — ALSO FAILS!")
else:
    print(f"  RESULT: resolved!")
    for i, v in enumerate(result2):
        print(f"    [{i}] = {v}")

# 4) Check GAPDH via Python cgranges 
print("\n=== Python cgranges query for GAPDH region ===")
cr = tx_index.cr
overlaps = list(cr.overlap("chr12", 6534806, 6534861))
print(f"  Overlaps for [6534806, 6534861): {len(overlaps)}")
for ref, start, end, label in overlaps[:20]:
    itype = tx_index._iv_type[label]
    tset = tx_index._iv_t_set[label]
    type_name = {0: "exon", 1: "transcript", 2: "intergenic"}
    print(f"    [{start}, {end}) type={type_name.get(itype, itype)} t_set={tset}")

overlaps2 = list(cr.overlap("chr12", 6536493, 6536588))
print(f"\n  Overlaps for [6536493, 6536588): {len(overlaps2)}")
for ref, start, end, label in overlaps2[:20]:
    itype = tx_index._iv_type[label]
    tset = tx_index._iv_t_set[label]
    type_name = {0: "exon", 1: "transcript", 2: "intergenic"}
    print(f"    [{start}, {end}) type={type_name.get(itype, itype)} t_set={tset}")
