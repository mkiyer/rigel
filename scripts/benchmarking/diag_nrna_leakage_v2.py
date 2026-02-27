#!/usr/bin/env python
"""Diagnostic: run the pipeline on the existing BAM and check nRNA counts.

Specifically, check how many fragments are assigned to nRNA for
ENST00000624314.1 (LINC02912), PVT1 transcripts, and overall.
Also produce an annotated BAM for detailed inspection.
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

import numpy as np
from hulkrna.index import HulkIndex
from hulkrna.pipeline import run_pipeline
from hulkrna.config import PipelineConfig

# Paths
INDEX_DIR = "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/hulkrna_index"
BAM_PATH = "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/gdna_none_nrna_none_ss_1.00/align_oracle/reads_namesort.bam"
OUT_BAM = "/tmp/diag_annotated_v2.bam"

print("Loading index...")
index = HulkIndex.load(INDEX_DIR)
print(f"  {len(index.t_df)} transcripts, {len(index.g_df)} genes")

# Find ENST00000624314.1
t_ids = index.t_df["t_id"].values
t624_mask = np.array([str(t) == "ENST00000624314.1" for t in t_ids])
if t624_mask.any():
    t624_idx = int(np.where(t624_mask)[0][0])
    print(f"  ENST00000624314.1 at index {t624_idx}")
    g_idx_624 = int(index.t_to_g_arr[t624_idx])
    print(f"  Gene: {index.g_df['g_id'].values[g_idx_624]} (g_idx={g_idx_624})")
    print(f"  Gene name: {index.g_df['g_name'].values[g_idx_624]}")
    print(f"  t_strand: {index.t_df['strand'].values[t624_idx]}")
    t_start = index.t_df["start"].values[t624_idx]
    t_end = index.t_df["end"].values[t624_idx]
    t_len = index.t_df["length"].values[t624_idx]
    print(f"  start={t_start}, end={t_end}, length={t_len}")
    print(f"  span={t_end - t_start}, is_single_exon={t_end - t_start == t_len}")
else:
    t624_idx = -1
    print("  ENST00000624314.1 NOT FOUND")

# Find PVT1 gene
g_names = index.g_df["g_name"].values
pvt1_mask = np.array([str(g) == "PVT1" for g in g_names])
if pvt1_mask.any():
    pvt1_gidx = int(np.where(pvt1_mask)[0][0])
    pvt1_g_id = index.g_df["g_id"].values[pvt1_gidx]
    print(f"\n  PVT1 gene: {pvt1_g_id} (g_idx={pvt1_gidx})")
    pvt1_t_mask = index.t_to_g_arr == pvt1_gidx
    pvt1_t_indices = np.where(pvt1_t_mask)[0]
    print(f"  PVT1 has {len(pvt1_t_indices)} transcripts")
else:
    pvt1_gidx = -1
    pvt1_t_indices = np.array([])
    print("  PVT1 NOT FOUND")

# Run pipeline with annotation
print(f"\nRunning pipeline on {BAM_PATH}...")
config = PipelineConfig(annotated_bam_path=OUT_BAM)
result = run_pipeline(BAM_PATH, index, config=config)

counter = result.estimator
stats = result.stats
strand_models = result.strand_models

print(f"\n--- Pipeline Stats ---")
print(f"  Total fragments: {stats.total_fragments:,}")
print(f"  Deterministic unique: {stats.deterministic_unique_units:,}")
print(f"  EM-routed unique: {stats.em_routed_unique_units:,}")
print(f"  EM-routed isoform-ambig: {stats.em_routed_isoform_ambig_units:,}")
print(f"  EM-routed gene-ambig: {stats.em_routed_gene_ambig_units:,}")
print(f"  EM-routed multimapper: {stats.em_routed_multimapper_units:,}")
print(f"  nRNA EM count total: {counter.nrna_em_count:.2f}")
print(f"  gDNA EM count total: {counter.gdna_em_count:.2f}")

# Strand model details
sm = strand_models.exonic_spliced
print(f"\n--- Strand Model ---")
print(f"  SS = {sm.strand_specificity:.4f}")
print(f"  p_r1_sense = {sm.p_r1_sense:.6f}")
print(f"  anti_flag = {sm.anti_flag}")

# nRNA counts for ENST00000624314.1
if t624_idx >= 0:
    nrna_624 = counter.nrna_em_counts[t624_idx]
    mrna_624 = counter.t_counts[t624_idx].sum()
    unique_624 = counter.unique_counts[t624_idx].sum()
    em_624 = counter.em_counts[t624_idx].sum()
    print(f"\n--- ENST00000624314.1 (LINC02912) ---")
    print(f"  unique counts: {unique_624:.2f}")
    print(f"  EM counts: {em_624:.2f}")
    print(f"  total mRNA: {mrna_624:.2f}")
    print(f"  nRNA EM: {nrna_624:.2f}")

# nRNA counts for PVT1
if len(pvt1_t_indices) > 0:
    nrna_pvt1 = counter.nrna_em_counts[pvt1_t_indices].sum()
    mrna_pvt1 = counter.t_counts[pvt1_t_indices].sum()
    unique_pvt1 = counter.unique_counts[pvt1_t_indices].sum()
    em_pvt1 = counter.em_counts[pvt1_t_indices].sum()
    print(f"\n--- PVT1 (all isoforms) ---")
    print(f"  unique counts: {unique_pvt1:.2f}")
    print(f"  EM counts: {em_pvt1:.2f}")
    print(f"  total mRNA: {mrna_pvt1:.2f}")
    print(f"  nRNA EM: {nrna_pvt1:.2f}")
    # Show top nRNA isoforms
    pvt1_nrna = counter.nrna_em_counts[pvt1_t_indices]
    top_nrna = np.argsort(pvt1_nrna)[::-1][:5]
    if pvt1_nrna.sum() > 0:
        print(f"  Top nRNA PVT1 isoforms:")
        for k in top_nrna:
            if pvt1_nrna[k] > 0:
                tidx = pvt1_t_indices[k]
                print(f"    {t_ids[tidx]}: {pvt1_nrna[k]:.2f}")

# Overall nRNA
top_nrna_all = np.argsort(counter.nrna_em_counts)[::-1][:10]
print(f"\n--- Top 10 nRNA transcripts overall ---")
for tidx in top_nrna_all:
    ncount = counter.nrna_em_counts[tidx]
    if ncount > 0:
        gidx = int(index.t_to_g_arr[tidx])
        gname = str(index.g_df["g_name"].values[gidx])
        print(f"  {t_ids[tidx]} ({gname}): {ncount:.2f}")

# Check annotated BAM
print(f"\n--- Annotated BAM ---")
print(f"  Written to: {OUT_BAM}")

import pysam
nrna_reads = 0
nrna_by_zt = {}
nrna_by_zg = {}
bam = pysam.AlignmentFile(OUT_BAM, "rb")
for read in bam.fetch(until_eof=True):
    if not read.is_read1:
        continue
    zp = read.get_tag("ZP") if read.has_tag("ZP") else None
    if zp == "nRNA":
        nrna_reads += 1
        zt = read.get_tag("ZT") if read.has_tag("ZT") else "."
        zg = read.get_tag("ZG") if read.has_tag("ZG") else "."
        nrna_by_zt[zt] = nrna_by_zt.get(zt, 0) + 1
        nrna_by_zg[zg] = nrna_by_zg.get(zg, 0) + 1
bam.close()

print(f"  Total nRNA R1 reads: {nrna_reads}")
if nrna_by_zt:
    print(f"  nRNA by ZT (top 10):")
    for zt, cnt in sorted(nrna_by_zt.items(), key=lambda x: -x[1])[:10]:
        print(f"    {zt}: {cnt}")
if nrna_by_zg:
    print(f"  nRNA by ZG (top 10):")
    for zg, cnt in sorted(nrna_by_zg.items(), key=lambda x: -x[1])[:10]:
        print(f"    {zg}: {cnt}")

print("\nDone.")
