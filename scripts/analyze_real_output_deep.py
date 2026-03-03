#!/usr/bin/env python3
"""Deep dive on nRNA, gDNA, and strand specificity from hulkrna output."""
import pandas as pd
import json
import numpy as np

outdir = '/Users/mkiyer/Downloads/hulkrna_runs/mctp_LBX0069_SI_42153_HFFFMDRX7/hulkrna_output'

gene_counts = pd.read_feather(f'{outdir}/gene_quant.feather')
counts = pd.read_feather(f'{outdir}/quant.feather')

# Key summary metrics
total_mrna = gene_counts['count'].sum()
total_nrna = gene_counts['n_nrna'].sum()
total_gdna = gene_counts['n_gdna'].sum()
total_counted = total_mrna + total_nrna + total_gdna

print('=== GLOBAL COMPOSITION ===')
print(f'Total mRNA:  {total_mrna:>12,.1f}  ({total_mrna/total_counted*100:>5.1f}%)')
print(f'Total nRNA:  {total_nrna:>12,.1f}  ({total_nrna/total_counted*100:>5.1f}%)')
print(f'Total gDNA:  {total_gdna:>12,.1f}  ({total_gdna/total_counted*100:>5.1f}%)')
print(f'Grand total: {total_counted:>12,.1f}')
print()

# nRNA rate as fraction of transcriptomic signal
nrna_rate = total_nrna / (total_mrna + total_nrna)
print(f'nRNA rate (nRNA / RNA): {nrna_rate:.4f} ({nrna_rate*100:.2f}%)')
gdna_rate = total_gdna / total_counted
print(f'gDNA rate (gDNA / total): {gdna_rate:.4f} ({gdna_rate*100:.2f}%)')
print()

# Strand specificity from strand model
print('=== STRAND SPECIFICITY ===')
with open(f'{outdir}/strand_model.json') as f:
    strand = json.load(f)

for model_name, model in strand['strand_models'].items():
    ss = model['probabilities']['strand_specificity']
    n_obs = model['observations']['total']
    print(f'  {model_name:25s}  strand_spec={ss:.6f}  n_obs={n_obs:>12,}')

# The relevant strand specificity is the spliced one (most reliable)
spliced_ss = strand['strand_models']['exonic_spliced']['probabilities']['strand_specificity']
print(f'\nUsing exonic_spliced strand_specificity: {spliced_ss:.6f}')

# Top nRNA genes
print('\n=== TOP 30 GENES BY nRNA COUNT ===')
top_nrna = gene_counts.nlargest(30, 'n_nrna')
print(f'{"gene_name":20s} {"count":>10s} {"n_nrna":>10s} {"n_gdna":>10s} {"nrna_frac":>10s}')
for _, row in top_nrna.iterrows():
    gname = row.get('gene_name', '')
    total = row['count'] + row['n_nrna']
    frac = row['n_nrna'] / total if total > 0 else 0
    print(f'{gname:20s} {row["count"]:>10.1f} {row["n_nrna"]:>10.1f} {row["n_gdna"]:>10.1f} {frac:>10.3f}')

# nRNA rate distribution
print('\n=== nRNA RATE DISTRIBUTION (per gene, expressed genes only) ===')
expressed = gene_counts[(gene_counts['count'] > 0) | (gene_counts['n_nrna'] > 0)].copy()
expressed['total_rna'] = expressed['count'] + expressed['n_nrna']
expressed['nrna_frac'] = expressed['n_nrna'] / expressed['total_rna']
expressed = expressed[expressed['total_rna'] > 1]  # filter noise

print(f'Genes with nRNA > 0: {(expressed["n_nrna"] > 0).sum():,} / {len(expressed):,}')
print(f'Median nRNA fraction: {expressed["nrna_frac"].median():.4f}')
print(f'Mean nRNA fraction: {expressed["nrna_frac"].mean():.4f}')
print(f'75th percentile: {expressed["nrna_frac"].quantile(0.75):.4f}')
print(f'90th percentile: {expressed["nrna_frac"].quantile(0.90):.4f}')
print(f'95th percentile: {expressed["nrna_frac"].quantile(0.95):.4f}')

# gDNA rate distribution
print('\n=== gDNA RATE PER GENE ===')
print(f'Genes with gDNA > 0: {(gene_counts["n_gdna"] > 0).sum():,}')
avg_gdna_rate = gene_counts.loc[gene_counts['n_gdna'] > 0, 'gdna_rate'].mean()
print(f'Mean per-gene gDNA rate: {avg_gdna_rate:.4f}')

# Fragment length model
print('\n=== FRAGMENT LENGTH MODEL ===')
with open(f'{outdir}/frag_length_models.json') as f:
    fl = json.load(f)
for key, val in fl.items():
    if isinstance(val, dict) and 'observations' in val:
        obs = val['observations']
        n = sum(obs.values()) if isinstance(obs, dict) else len(obs)
        print(f'  {key}: {n} items in observations')
    elif isinstance(val, dict):
        for k2, v2 in val.items():
            if isinstance(v2, (int, float)):
                print(f'  {key}.{k2}: {v2}')
            elif isinstance(v2, dict):
                print(f'  {key}.{k2}: {list(v2.keys())[:5]}...')
            elif isinstance(v2, list) and len(v2) < 5:
                print(f'  {key}.{k2}: {v2}')

# Intronic vs exonic counts for sense/antisense
print('\n=== SENSE/ANTISENSE BREAKDOWN (gene level) ===')
cols = ['n_antisense', 'n_sense_all', 'n_antisense_all', 'n_intronic_sense', 'n_intronic_antisense']
for col in cols:
    if col in gene_counts.columns:
        total = gene_counts[col].sum()
        nonzero = (gene_counts[col] > 0).sum()
        print(f'  {col:25s}  total={total:>12,.1f}  nonzero={nonzero:>8,}')

print('\n=== SUMMARY FOR SIMULATION CONFIG ===')
print(f'gDNA rate:           {gdna_rate:.4f}')
print(f'nRNA rate:           {nrna_rate:.4f}')
print(f'Strand specificity:  {spliced_ss:.6f}')
