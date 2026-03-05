#!/usr/bin/env python
"""Quick script to inspect data formats for full-genome benchmark."""
import pandas as pd

# Check quant.feather
counts = pd.read_feather(
    "/Users/mkiyer/Downloads/hulkrna_runs/mctp_LBX0069_SI_42153_HFFFMDRX7/"
    "hulkrna_output/quant.feather"
)
print("=== quant.feather ===")
print("Shape:", counts.shape)
print("Columns:", list(counts.columns))
print(counts.head(3).to_string())
print()
print("Sum of mrna_em for expressed:", counts[counts["mrna_em"] > 0]["mrna_em"].sum())
print("Expressed transcripts:", (counts["mrna_em"] > 0).sum())
print()

# Check hulkrna index transcripts
txdf = pd.read_feather(
    "/Users/mkiyer/Downloads/hulkrna_runs/refs/human/hulkrna_index/transcripts.feather"
)
print("=== transcripts.feather (index) ===")
print("Shape:", txdf.shape)
print("Columns:", list(txdf.columns))
print(txdf.head(3).to_string())
print()

# Check ref_lengths
rl = pd.read_feather(
    "/Users/mkiyer/Downloads/hulkrna_runs/refs/human/hulkrna_index/ref_lengths.feather"
)
print("=== ref_lengths.feather ===")
print("Shape:", rl.shape)
print("Columns:", list(rl.columns))
print(rl.head())
print("Num refs:", len(rl))
print()

# Check genome FASTA index
print("=== Genome FAI ===")
import pysam
fa = pysam.FastaFile(
    "/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genome_controls.fasta.bgz"
)
refs = fa.references
print(f"Num references: {len(refs)}")
print(f"First 10: {refs[:10]}")
print(f"Last 5: {refs[-5:]}")
fa.close()

# Check if salmon/kallisto indices exist 
import shutil
print()
print("=== Tool availability ===")
print(f"salmon: {shutil.which('salmon')}")
print(f"kallisto: {shutil.which('kallisto')}")
