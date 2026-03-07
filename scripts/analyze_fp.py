#!/usr/bin/env python3
"""Investigate rigel false positives vs salmon."""
import pandas as pd
import numpy as np

base = "/Users/mkiyer/Downloads/rigel_runs/benchmark_output/gdna_none_ss_1.00_nrna_none"
df = pd.read_csv(f"{base}/per_transcript_counts_oracle.csv")

truth = df['mrna_truth'].astype(float)
rigel = df['rigel_oracle'].astype(float)
salmon = df['salmon'].astype(float)

# False positives: truth == 0 but pred > 0
fp_rigel = df[(truth == 0) & (rigel > 0)].copy()
fp_salmon = df[(truth == 0) & (salmon > 0)].copy()

print(f"Rigel false positives: {len(fp_rigel)}")
print(f"Salmon false positives: {len(fp_salmon)}")
print()

# How big are the false positive counts?
print("=== Rigel FP count distribution ===")
fp_rigel_counts = fp_rigel['rigel_oracle']
for pct in [50, 75, 90, 95, 99, 99.9, 100]:
    print(f"  p{pct}: {np.percentile(fp_rigel_counts, pct):.4f}")
print(f"  mean: {fp_rigel_counts.mean():.4f}")
print(f"  sum:  {fp_rigel_counts.sum():.0f}")
print(f"  count > 1: {(fp_rigel_counts > 1).sum()}")
print(f"  count > 10: {(fp_rigel_counts > 10).sum()}")
print(f"  count > 100: {(fp_rigel_counts > 100).sum()}")
print()

print("=== Salmon FP count distribution ===")
fp_salmon_counts = fp_salmon['salmon']
for pct in [50, 75, 90, 95, 99, 99.9, 100]:
    print(f"  p{pct}: {np.percentile(fp_salmon_counts, pct):.4f}")
print(f"  mean: {fp_salmon_counts.mean():.4f}")
print(f"  sum:  {fp_salmon_counts.sum():.0f}")
print(f"  count > 1: {(fp_salmon_counts > 1).sum()}")
print(f"  count > 10: {(fp_salmon_counts > 10).sum()}")
print(f"  count > 100: {(fp_salmon_counts > 100).sum()}")
print()

# Are rigel FPs mostly in multi-transcript genes?
# Check top rigel false positives
top_fp = fp_rigel.nlargest(20, 'rigel_oracle')[['transcript_id', 'gene_id', 'gene_name', 'rigel_oracle', 'salmon']]
print("=== Top 20 rigel false positives ===")
print(top_fp.to_string(index=False))
print()

# Check if these FPs are in genes that DO have expressed transcripts
expressed_genes = set(df[truth > 0]['gene_id'].unique())
fp_in_expressed_gene = fp_rigel['gene_id'].isin(expressed_genes).sum()
fp_not_in_expressed_gene = len(fp_rigel) - fp_in_expressed_gene
print(f"FPs in genes with expressed transcripts: {fp_in_expressed_gene}")
print(f"FPs in fully unexpressed genes: {fp_not_in_expressed_gene}")
print()

# Spearman only on expressed transcripts (truth > 0 AND pred > 0)
def spearman_r(x, y):
    rx = np.empty_like(x)
    ry = np.empty_like(y)
    rx[np.argsort(x)] = np.arange(len(x), dtype=float)
    ry[np.argsort(y)] = np.arange(len(y), dtype=float)
    d = rx - ry
    n = len(x)
    return 1.0 - 6.0 * np.sum(d**2) / (n * (n**2 - 1))

# On expressed only
expressed = truth > 0
t_exp = truth[expressed].values
r_exp = rigel[expressed].values
s_exp = salmon[expressed].values

print("=== Metrics on EXPRESSED transcripts only ===")
print(f"N expressed: {expressed.sum()}")
print(f"Rigel Pearson:  {np.corrcoef(t_exp, r_exp)[0,1]:.6f}")
print(f"Salmon Pearson: {np.corrcoef(t_exp, s_exp)[0,1]:.6f}")
print(f"Rigel Spearman:  {spearman_r(t_exp, r_exp):.6f}")
print(f"Salmon Spearman: {spearman_r(t_exp, s_exp):.6f}")
print(f"Rigel MAE:  {np.mean(np.abs(t_exp - r_exp)):.4f}")
print(f"Salmon MAE: {np.mean(np.abs(t_exp - s_exp)):.4f}")
print(f"Rigel RMSE:  {np.sqrt(np.mean((t_exp - r_exp)**2)):.4f}")
print(f"Salmon RMSE: {np.sqrt(np.mean((t_exp - s_exp)**2)):.4f}")
