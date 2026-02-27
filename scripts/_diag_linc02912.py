#!/usr/bin/env python3
"""Diagnose why ENST00000624314.1 nRNA component is alive in the locus EM."""
import sys
import numpy as np

sys.path.insert(0, "src")

from hulkrna.index import HulkIndex
from hulkrna.pipeline import run_pipeline
from hulkrna.config import PipelineConfig
from pathlib import Path

INDEX_DIR = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/hulkrna_index")
BAM_PATH = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/gdna_none_nrna_none_ss_1.00/align_oracle/reads_namesort.bam")

print("Loading index...")
index = HulkIndex.load(INDEX_DIR)

# Find ENST00000624314.1 in the index
t_df = index.t_df
target_tid = "ENST00000624314.1"
mask = t_df["t_id"] == target_tid
if mask.sum() == 0:
    print(f"ERROR: {target_tid} not found in index!")
    sys.exit(1)

t_idx = int(t_df.index[mask][0])
row = t_df.loc[t_idx]
print(f"\n=== {target_tid} in index ===")
print(f"  t_index = {t_idx}")
print(f"  gene = {row.g_id}")
print(f"  start = {row.start}, end = {row.end}")
print(f"  length (exonic) = {row.length}")
print(f"  strand = {row.strand}")
span = row.end - row.start
intronic = span - row.length
print(f"  genomic span = {span}")
print(f"  intronic span = {intronic}")
print(f"  span <= exonic_length? = {span <= row.length}")

print("\n--- Running pipeline to inspect internals ---")
cfg = PipelineConfig()
pipe = run_pipeline(BAM_PATH, index, config=cfg)

# Check nrna_init for this transcript
nrna_init = pipe.estimator.nrna_init
print(f"\n=== nRNA init ===")
print(f"  nrna_init[{t_idx}] = {nrna_init[t_idx]}")
print(f"  nrna_init range [{nrna_init.min():.4f}, {nrna_init.max():.4f}]")

# Check transcript_spans & exonic_lengths
est = pipe.estimator
print(f"\n=== Estimator internals ===")
if hasattr(est, '_transcript_spans') and est._transcript_spans is not None:
    print(f"  _transcript_spans[{t_idx}] = {est._transcript_spans[t_idx]}")
else:
    print(f"  _transcript_spans = None")
if hasattr(est, '_exonic_lengths') and est._exonic_lengths is not None:
    print(f"  _exonic_lengths[{t_idx}] = {est._exonic_lengths[t_idx]}")
else:
    print(f"  _exonic_lengths = None")

# Check unique counts
unique = est.unique_counts[t_idx].sum()
print(f"  unique_counts[{t_idx}] (sum) = {unique}")

# Check total counts
t_counts = est.t_counts
observed = float(t_counts[t_idx].sum())
print(f"  t_counts[{t_idx}] (sum) = {observed}")

# Check nrna_em_counts
nrna_em = est.nrna_em_counts
print(f"\n=== nRNA EM counts ===")
print(f"  nrna_em_counts[{t_idx}] = {nrna_em[t_idx]}")
print(f"  nrna_em total = {nrna_em.sum():.1f}")

# Show top nRNA EM counts
top_nrna = np.argsort(nrna_em)[::-1][:10]
print(f"\n  Top 10 nRNA EM counts:")
for i in top_nrna:
    if nrna_em[i] > 0:
        tid = t_df.loc[i, "t_id"]
        gid = t_df.loc[i, "g_id"]
        sp = t_df.loc[i, "end"] - t_df.loc[i, "start"]
        el = t_df.loc[i, "length"]
        print(f"    [{i}] {tid} ({gid}) span={sp} exonic_len={el} nrna_em={nrna_em[i]:.1f}")

# Check all nrna_init for this gene
gene_mask = t_df["g_id"] == row.g_id
gene_indices = t_df.index[gene_mask].tolist()
print(f"\n=== All transcripts in {row.g_id} ===")
for gi in gene_indices:
    gr = t_df.loc[gi]
    print(f"  [{gi}] {gr.t_id}: span={gr.end-gr.start}, exonic_len={gr.length}, nrna_init={nrna_init[gi]:.4f}")

print("\n=== Pipeline stats ===")
stats = pipe.stats
print(f"  n_fragments = {stats.n_fragments}")
print(f"  n_intergenic = {stats.n_intergenic}")
print(f"  n_gdna = {stats.n_gdna_total}")
print(f"  nrna_em total = {nrna_em.sum():.1f}")
print(f"  nrna_init total = {nrna_init.sum():.1f}")
