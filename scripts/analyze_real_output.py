#!/usr/bin/env python3
"""Analyze hulkrna output from real data run."""
import pandas as pd
import json
import numpy as np

outdir = '/Users/mkiyer/Downloads/hulkrna_runs/mctp_LBX0069_SI_42153_HFFFMDRX7/hulkrna_output'

# 1. Transcript-level counts
counts = pd.read_feather(f'{outdir}/quant.feather')
print('=== TRANSCRIPT COUNTS ===')
print(f'Shape: {counts.shape}')
print(f'Columns: {list(counts.columns)}')
print()

# Filter to expressed transcripts
expressed = counts[counts['mrna'] > 0]
print(f'Expressed transcripts: {len(expressed):,} / {len(counts):,}')
print(f'Total mRNA count: {expressed["mrna"].sum():,.1f}')

# Top 20 transcripts by mrna
top20 = expressed.nlargest(20, 'mrna')
print()
print('Top 20 transcripts by mrna:')
for _, row in top20.iterrows():
    t_id = row.get('transcript_id', '')
    gname = row.get('gene_name', '')
    print(f'  {t_id:30s}  {gname:15s}  mrna={row["mrna"]:,.1f}')

# nRNA columns
nrna_cols = [c for c in counts.columns if 'nrna' in c.lower() or 'nascent' in c.lower()]
print(f'\nnRNA columns: {nrna_cols}')
if nrna_cols:
    for col in nrna_cols:
        total = counts[col].sum()
        nonzero = (counts[col] > 0).sum()
        print(f'  {col} total: {total:,.1f}, nonzero: {nonzero:,}')

print()
print('=== GENE COUNTS ===')
gene_counts = pd.read_feather(f'{outdir}/gene_quant.feather')
print(f'Shape: {gene_counts.shape}')
print(f'Columns: {list(gene_counts.columns)}')
expressed_genes = gene_counts[gene_counts['mrna'] > 0]
print(f'Expressed genes: {len(expressed_genes):,} / {len(gene_counts):,}')
print(f'Total gene-level mRNA count: {expressed_genes["mrna"].sum():,.1f}')

# Gene-level nRNA
for col in gene_counts.columns:
    if 'nrna' in col.lower():
        total = gene_counts[col].sum()
        nonzero = (gene_counts[col] > 0).sum()
        print(f'  {col}: total={total:,.1f}, nonzero={nonzero:,}')

# Top 20 genes by mrna
top20g = expressed_genes.nlargest(20, 'mrna')
print()
print('Top 20 genes:')
for _, row in top20g.iterrows():
    gname = row.get('gene_name', '')
    print(f'  {gname:20s}  mrna={row["mrna"]:>8.1f}  nrna={row["nrna"]:>8.1f}  rna_total={row["rna_total"]:>8.1f}')

print()
print('=== SUMMARY ===')
with open(f'{outdir}/summary.json') as f:
    summary = json.load(f)

# Quantification rates
quant = summary.get('quantification', {})
frags = summary.get('fragments', {})
print(f'  n_fragments: {frags.get("n_fragments", 0):,}')
print(f'  mrna_total: {quant.get("mrna_total", 0):,.1f}')
print(f'  nrna_total: {quant.get("nrna_total", 0):,.1f}')
print(f'  gdna_total: {quant.get("gdna_total", 0):,.1f}')
print(f'  mrna_fraction: {quant.get("mrna_fraction", 0):.4f}')
print(f'  gdna_fraction: {quant.get("gdna_fraction", 0):.4f}')

# Strand model
print()
print('=== STRAND MODEL ===')
strand = summary.get('strand_models', {})
print(json.dumps(strand, indent=2)[:2000])

# Fragment length model
print()
print('=== FRAGMENT LENGTH MODEL ===')
fl = summary.get('frag_length_models', {})
for key in fl:
    if isinstance(fl[key], dict):
        info = fl[key]
        n_obs = info.get('n_observations', info.get('n_obs', '?'))
        mode = info.get('mode', '?')
        mean = info.get('mean', '?')
        print(f'  {key}: n_obs={n_obs}, mode={mode}, mean={mean}')

# Pipeline stats
print()
print('=== PIPELINE STATS ===')
pstats = summary.get('pipeline_stats', {})
for k in ['n_fragments', 'n_unique_gene', 'n_multi_gene', 'n_intergenic']:
    if k in pstats:
        print(f'  {k}: {pstats[k]:,}')

print()
print('=== LOCI ===')
import os
loci_path = f'{outdir}/loci.feather'
if os.path.exists(loci_path):
    loci = pd.read_feather(loci_path)
    print(f'Shape: {loci.shape}')
    print(f'Columns: {list(loci.columns)}')
    print(f'Total gdna: {loci["gdna"].sum():,.1f}')
    print(f'Loci with gdna > 0: {(loci["gdna"] > 0).sum():,}')
    top_gdna = loci.nlargest(10, 'gdna')
    print('Top 10 loci by gDNA:')
    for _, row in top_gdna.iterrows():
        print(f'  locus {int(row["locus_id"]):>5d}  chrom={row["chrom"]}  gdna={row["gdna"]:>8.1f}  gdna_rate={row["gdna_rate"]:.3f}')

print()
print('=== DETAIL COUNTS ===')
detail = pd.read_feather(f'{outdir}/quant_detail.feather')
print(f'Shape: {detail.shape}')
print(f'Columns: {list(detail.columns)}')

# Summary of numeric columns
for col in detail.columns:
    if detail[col].dtype in [np.float64, np.float32, np.int64, np.int32]:
        total = detail[col].sum()
        if total > 0:
            print(f'  {col}: total={total:,.1f}')

# nRNA analysis
print()
print('=== nRNA ANALYSIS ===')
if 'nrna' in counts.columns:
    nrna_total = counts['nrna'].sum()
    mrna_total = counts['mrna'].sum()
    print(f'Total nRNA count: {nrna_total:,.1f}')
    print(f'Total mRNA count: {mrna_total:,.1f}')
    print(f'nRNA / (mRNA + nRNA): {nrna_total / (mrna_total + nrna_total) * 100:.2f}%')

    # Top nRNA genes
    if 'nrna' in gene_counts.columns:
        top_nrna = gene_counts.nlargest(15, 'nrna')
        print()
        print('Top 15 genes by nRNA:')
        for _, row in top_nrna.iterrows():
            gname = row.get('gene_name', '')
            print(f'  {gname:20s}  nrna={row["nrna"]:>8.1f}  mrna={row["mrna"]:>8.1f}  '
                  f'nrna_frac={row["nrna"]/(row["mrna"]+row["nrna"]+1e-9):.3f}')
