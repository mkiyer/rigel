#!/usr/bin/env python3
"""Check POLR2A, OAZ1, and RPL41 locus-level transcript counts."""
import sys
import pandas as pd
import pyarrow.feather as pf

BASE = '/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine'
SUBDIR = f'{BASE}/gdna_none_ss_0.95_nrna_none'

# per_transcript_counts_oracle.csv has: transcript_id, gene_id, gene_name, mrna_abundance, nrna_abundance, mrna_truth, nrna_truth, kallisto, rigel_oracle, salmon
oracle_df = pd.read_csv(f'{SUBDIR}/per_transcript_counts_oracle.csv')
# per_transcript_counts_minimap2.csv has: transcript_id, gene_id, gene_name, mrna_abundance, nrna_abundance, mrna_truth, nrna_truth, rigel_minimap2
mm2_df = pd.read_csv(f'{SUBDIR}/per_transcript_counts_minimap2.csv')
txp = pf.read_feather(f'{BASE}/rigel_index/transcripts.feather')

def check_gene(gene_name=None, gene_id=None):
    if gene_name:
        ids = txp[txp['g_name'] == gene_name]['t_id'].tolist()
        gtypes = txp[txp['g_name'] == gene_name][['t_id', 'g_type']].set_index('t_id')['g_type'].to_dict()
    else:
        ids = txp[txp['g_id'].str.startswith(gene_id)]['t_id'].tolist()
        gtypes = txp[txp['g_id'].str.startswith(gene_id)][['t_id', 'g_type']].set_index('t_id')['g_type'].to_dict()

    orr = oracle_df[oracle_df['transcript_id'].isin(ids)][
        ['transcript_id', 'gene_name', 'mrna_truth', 'rigel_oracle', 'kallisto', 'salmon']].copy()
    mr = mm2_df[mm2_df['transcript_id'].isin(ids)][['transcript_id', 'rigel_minimap2']]
    df = orr.merge(mr, on='transcript_id', how='left')
    df['g_type'] = df['transcript_id'].map(gtypes)
    df = df.sort_values('mrna_truth', ascending=False)
    print(df.rename(columns={'mrna_truth': 'truth', 'rigel_oracle': 'oracle'}).to_string(index=False))
    print()


print("=== POLR2A - check all transcripts at this locus ===")
check_gene(gene_name='POLR2A')

print("=== OAZ1 - check all transcripts at this locus ===")
check_gene(gene_name='OAZ1')

print("=== RPL41 - check all transcripts at this locus ===")
check_gene(gene_name='RPL41')

# Check pseudogenes absorption for POLR2A - look for anything with POLR2A-like name
print("=== Genes with 'POLR2' in name ===")
polr = txp[txp['g_name'].str.contains('POLR2A', na=False)][['t_id','g_id','g_name','g_type']].drop_duplicates('g_name')
print(polr.to_string(index=False))
print()

# OAZ paralogs
print("=== Genes with 'OAZ' in name ===")
oaz = txp[txp['g_name'].str.contains('^OAZ', na=False)][['t_id','g_id','g_name','g_type']].drop_duplicates('g_name')
print(oaz.to_string(index=False))

# Check OAZ2 and OAZ3 counts
print()
print("=== OAZ2 locus ===")
check_gene(gene_name='OAZ2')

print("=== OAZ3 locus ===")
check_gene(gene_name='OAZ3')

# Also check ENSG00000280800 locus and ENSG00000278996 (overlapping lncRNA)
print("=== ENSG00000280800 locus ===")
check_gene(gene_id='ENSG00000280800')

print("=== ENSG00000278996 (nearby lncRNA on chr21) ===")
check_gene(gene_id='ENSG00000278996')
