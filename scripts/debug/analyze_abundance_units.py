#!/usr/bin/env python3
"""Investigate transcript abundance units and MARD computation."""
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr

base = "/Users/mkiyer/Downloads/rigel_runs/sim_synthetic"
truth = pd.read_csv(f"{base}/truth_abundances.tsv", sep="\t")
quant = pd.read_feather(f"{base}/gdna_none_ss_0.99_nrna_none/rigel_out/quant.feather")

merged = truth.merge(
    quant[["transcript_id", "count", "count_em", "tpm"]],
    on="transcript_id", how="left"
).fillna(0)
expressed = merged[merged["mrna_abundance"] > 0].copy()

print("TRUTH COLUMN STATS:")
print(f"  mrna_abundance range: [{expressed['mrna_abundance'].min():.2f}, {expressed['mrna_abundance'].max():.2f}]")
print(f"  mrna_abundance sum (expressed): {expressed['mrna_abundance'].sum():.2f}")
total_all = truth["mrna_abundance"].sum()
print(f"  mrna_abundance sum (all): {total_all:.2f}")
print(f"  n_expressed: {len(expressed)}")
print()

print("RIGEL OUTPUT STATS:")
print(f"  tpm range: [{expressed['tpm'].min():.2f}, {expressed['tpm'].max():.2f}]")
print(f"  tpm sum (expressed): {expressed['tpm'].sum():.2f}")
print(f"  count range: [{expressed['count'].min():.2f}, {expressed['count'].max():.2f}]")
print(f"  count sum (expressed): {expressed['count'].sum():.2f}")
print()

# The truth "mrna_abundance" is a RELATIVE abundance (like TPM without length normalization)
# Rigel TPM is length-normalized. They're in different units!
# To compare properly:
# Option A: Compare TPM-to-TPM (normalize truth by length)
# Option B: Compare count-to-expected-count

# Option A: truth abundance → TPM (divide by spliced_length, normalize to 1M)
expressed["truth_rpk"] = expressed["mrna_abundance"] / expressed["spliced_length"] * 1000
total_rpk = expressed["truth_rpk"].sum()
expressed["truth_tpm"] = expressed["truth_rpk"] / total_rpk * 1e6

re_tpm = np.abs(expressed["tpm"] - expressed["truth_tpm"]) / (expressed["truth_tpm"] + 1e-6)
sp_r, _ = spearmanr(expressed["truth_tpm"], expressed["tpm"])
pe_r, _ = pearsonr(np.log2(expressed["truth_tpm"] + 1), np.log2(expressed["tpm"] + 1))

print("=" * 80)
print("COMPARISON A: Rigel TPM vs Truth TPM (length-normalized)")
print("=" * 80)
print(f"  Spearman: {sp_r:.4f}")
print(f"  Pearson (log2): {pe_r:.4f}")
print(f"  MARD: {re_tpm.mean():.3f}")
print(f"  Median RE: {re_tpm.median():.3f}")
print()

# Option B: Compare counts to expected counts
# Expected count = abundance * length / sum(abundance * length) * N_total
expressed["eff_weight"] = expressed["mrna_abundance"] * expressed["spliced_length"]
total_weight = expressed["eff_weight"].sum()
expressed["expected_count"] = expressed["eff_weight"] / total_weight * 1_000_000

re_count = np.abs(expressed["count"] - expressed["expected_count"]) / (expressed["expected_count"] + 1e-6)
sp_c, _ = spearmanr(expressed["expected_count"], expressed["count"])
pe_c, _ = pearsonr(np.log2(expressed["expected_count"] + 1), np.log2(expressed["count"] + 1))

print("=" * 80)
print("COMPARISON B: Rigel count vs Expected count (abundance×length)")
print("=" * 80)
print(f"  Spearman: {sp_c:.4f}")
print(f"  Pearson (log2): {pe_c:.4f}")
print(f"  MARD: {re_count.mean():.3f}")
print(f"  Median RE: {re_count.median():.3f}")
print()

# Show worst cases
print("TOP 10 WORST TPM ERRORS:")
expressed["re_tpm"] = re_tpm
worst = expressed.nlargest(10, "re_tpm")
print(f"  {'transcript_id':<18} {'truth_tpm':>10} {'rigel_tpm':>10} {'RE':>8} {'n_exons':>7} {'spl_len':>7}")
for _, row in worst.iterrows():
    print(f"  {row['transcript_id']:<18} {row['truth_tpm']:>10.1f} {row['tpm']:>10.1f} {row['re_tpm']:>8.2f} {int(row['n_exons']):>7} {int(row['spliced_length']):>7}")

print()
print("TOP 10 WORST COUNT ERRORS:")
expressed["re_count"] = re_count
worst_c = expressed.nlargest(10, "re_count")
print(f"  {'transcript_id':<18} {'exp_count':>10} {'rigel_count':>11} {'RE':>8} {'n_exons':>7} {'spl_len':>7}")
for _, row in worst_c.iterrows():
    print(f"  {row['transcript_id']:<18} {row['expected_count']:>10.1f} {row['count']:>11.1f} {row['re_count']:>8.2f} {int(row['n_exons']):>7} {int(row['spliced_length']):>7}")

# Option C: raw comparison (what the original MARD was measuring)
re_raw = np.abs(expressed["tpm"] - expressed["mrna_abundance"]) / (expressed["mrna_abundance"] + 1e-6)
print()
print("=" * 80)
print("COMPARISON C: Rigel TPM vs RAW truth abundance (WRONG - different units!)")
print("=" * 80)
print(f"  MARD: {re_raw.mean():.3f}")
print(f"  Median RE: {re_raw.median():.3f}")
print("  → This is what the original analysis computed — MEANINGLESS because")
print("    truth 'mrna_abundance' is NOT TPM (not length-normalized)!")
