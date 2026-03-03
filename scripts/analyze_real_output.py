#!/usr/bin/env python3
"""Analyze hulkrna output from real data run."""
import pandas as pd
import json
import numpy as np

outdir = '/Users/mkiyer/Downloads/hulkrna_runs/mctp_LBX0069_SI_42153_HFFFMDRX7/hulkrna_output'

# 1. Transcript-level counts
counts = pd.read_feather(f'{outdir}/counts.feather')
print('=== TRANSCRIPT COUNTS ===')
print(f'Shape: {counts.shape}')
print(f'Columns: {list(counts.columns)}')
print()

# Filter to expressed transcripts
expressed = counts[counts['count'] > 0]
print(f'Expressed transcripts: {len(expressed):,} / {len(counts):,}')
print(f'Total mRNA count: {expressed["count"].sum():,.1f}')

# Top 20 transcripts by count
top20 = expressed.nlargest(20, 'count')
print()
print('Top 20 transcripts by count:')
for _, row in top20.iterrows():
    t_id = row.get('transcript_id', '')
    gname = row.get('gene_name', '')
    print(f'  {t_id:30s}  {gname:15s}  count={row["count"]:,.1f}')

# nRNA columns
nrna_cols = [c for c in counts.columns if 'nrna' in c.lower() or 'nascent' in c.lower()]
print(f'\nnRNA columns: {nrna_cols}')
if nrna_cols:
    for col in nrna_cols:
        total = counts[col].sum()
        nonzero = (counts[col] > 0).sum()
        print(f'  {col} total: {total:,.1f}, nonzero: {nonzero:,}')

# gDNA columns
gdna_cols = [c for c in counts.columns if 'gdna' in c.lower()]
print(f'\ngDNA columns: {gdna_cols}')
if gdna_cols:
    for col in gdna_cols:
        total = counts[col].sum()
        print(f'  {col} total: {total:,.1f}')

print()
print('=== GENE COUNTS ===')
gene_counts = pd.read_feather(f'{outdir}/gene_counts.feather')
print(f'Shape: {gene_counts.shape}')
print(f'Columns: {list(gene_counts.columns)}')
expressed_genes = gene_counts[gene_counts['count'] > 0]
print(f'Expressed genes: {len(expressed_genes):,} / {len(gene_counts):,}')
print(f'Total gene-level mRNA count: {expressed_genes["count"].sum():,.1f}')

# Gene-level nRNA and gDNA
for col in gene_counts.columns:
    if 'nrna' in col.lower() or 'gdna' in col.lower():
        total = gene_counts[col].sum()
        nonzero = (gene_counts[col] > 0).sum()
        print(f'  {col}: total={total:,.1f}, nonzero={nonzero:,}')

# Top 20 genes by count
top20g = expressed_genes.nlargest(20, 'count')
print()
print('Top 20 genes:')
for _, row in top20g.iterrows():
    g_id = row.get('gene_id', '')
    gname = row.get('gene_name', '')
    cnt = row['count']
    nrna = row.get('nrna_count', 0)
    gdna = row.get('gdna_count', 0)
    print(f'  {gname:20s}  count={cnt:>8.1f}  nRNA={nrna:>8.1f}  gDNA={gdna:>8.1f}')

print()
print('=== STATS ===')
with open(f'{outdir}/stats.json') as f:
    stats = json.load(f)
# Key stats
key_stats = ['n_fragments', 'n_unique_gene', 'n_multi_gene',
             'n_intergenic', 'n_gdna_total', 'n_gdna_unique', 'n_gdna_em',
             'duplicate', 'total', 'unique', 'multimapping']
for k in key_stats:
    if k in stats:
        print(f'  {k}: {stats[k]:,}')

# Compute rates
n_frags = stats.get('n_fragments', 1)
n_gdna = stats.get('n_gdna_total', 0)
n_intergenic = stats.get('n_intergenic', 0)
print(f'\n  gDNA rate: {n_gdna / n_frags:.4f} ({n_gdna / n_frags * 100:.2f}%)')
print(f'  Intergenic rate: {n_intergenic / n_frags:.4f} ({n_intergenic / n_frags * 100:.2f}%)')

print()
print('=== STRAND MODEL ===')
with open(f'{outdir}/strand_model.json') as f:
    strand = json.load(f)
print(json.dumps(strand, indent=2)[:2000])

print()
print('=== FRAGMENT LENGTH MODEL ===')
with open(f'{outdir}/frag_length_models.json') as f:
    fl = json.load(f)
# Just summary
for key in fl:
    if isinstance(fl[key], dict):
        info = fl[key]
        n_obs = info.get('n_observations', info.get('n_obs', '?'))
        mode = info.get('mode', '?')
        mean = info.get('mean', '?')
        print(f'  {key}: n_obs={n_obs}, mode={mode}, mean={mean}')

print()
print('=== DETAIL COUNTS ===')
detail = pd.read_feather(f'{outdir}/counts_detail.feather')
print(f'Shape: {detail.shape}')
print(f'Columns: {list(detail.columns)}')

# Summary of count categories
for col in detail.columns:
    if 'count' in col.lower() or col in ['mrna', 'nrna', 'gdna']:
        if detail[col].dtype in [np.float64, np.float32, np.int64, np.int32]:
            total = detail[col].sum()
            if total > 0:
                print(f'  {col}: total={total:,.1f}')

# nRNA analysis
print()
print('=== nRNA ANALYSIS ===')
if 'nrna_count' in counts.columns:
    nrna_total = counts['nrna_count'].sum()
    mrna_total = counts['count'].sum()
    print(f'Total nRNA count: {nrna_total:,.1f}')
    print(f'Total mRNA count: {mrna_total:,.1f}')
    print(f'nRNA / (mRNA + nRNA): {nrna_total / (mrna_total + nrna_total) * 100:.2f}%')
    print(f'nRNA rate (nRNA / total counted): {nrna_total / (mrna_total + nrna_total + n_gdna):.4f}')

    # Top nRNA genes
    if 'nrna_count' in gene_counts.columns:
        top_nrna = gene_counts.nlargest(15, 'nrna_count')
        print()
        print('Top 15 genes by nRNA count:')
        for _, row in top_nrna.iterrows():
            gname = row.get('gene_name', '')
            print(f'  {gname:20s}  nRNA={row["nrna_count"]:>8.1f}  mRNA={row["count"]:>8.1f}  '
                  f'nRNA_frac={row["nrna_count"]/(row["count"]+row["nrna_count"]+1e-9):.3f}')
