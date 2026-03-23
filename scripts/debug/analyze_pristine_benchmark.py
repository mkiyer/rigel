#!/usr/bin/env python3
"""Deep analysis of pristine benchmark results.

Compares rigel (oracle + minimap2), salmon, and kallisto at
the per-transcript level to identify systematic error patterns.
"""
import numpy as np
import pandas as pd

BASE = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/gdna_none_ss_0.95_nrna_none"

# ── Load oracle comparison (rigel vs salmon vs kallisto) ──────
oracle_df = pd.read_csv(f"{BASE}/per_transcript_counts_oracle.csv")
mm2_df = pd.read_csv(f"{BASE}/per_transcript_counts_minimap2.csv")

print("=" * 80)
print("PRISTINE BENCHMARK ANALYSIS (zero gDNA, zero nRNA, SS=0.95)")
print("=" * 80)

# ── 1. High-level summary ────────────────────────────────────
print("\n## 1. Aggregate Metrics (Oracle Alignment)")
print(f"  Transcripts in annotation: {len(oracle_df):,}")
expressed = oracle_df[oracle_df["mrna_truth"] > 0]
print(f"  Expressed transcripts (truth > 0): {len(expressed):,}")
print(f"  Unexpressed transcripts: {len(oracle_df) - len(expressed):,}")

total_truth = oracle_df["mrna_truth"].sum()
print(f"  Total mRNA truth fragments: {total_truth:,.0f}")

for tool in ["rigel_oracle", "salmon", "kallisto"]:
    pred = oracle_df[tool].sum()
    mae = (oracle_df[tool] - oracle_df["mrna_truth"]).abs().mean()
    rmse = np.sqrt(((oracle_df[tool] - oracle_df["mrna_truth"]) ** 2).mean())
    # Pearson, Spearman on expressed only
    from scipy.stats import pearsonr, spearmanr
    mask = oracle_df["mrna_truth"] > 0
    r_p, _ = pearsonr(oracle_df.loc[mask, "mrna_truth"], oracle_df.loc[mask, tool])
    r_s, _ = spearmanr(oracle_df.loc[mask, "mrna_truth"], oracle_df.loc[mask, tool])
    print(f"\n  {tool}:")
    print(f"    Total predicted: {pred:,.1f}")
    print(f"    Total − truth:   {pred - total_truth:+,.1f} ({100*(pred - total_truth)/total_truth:+.3f}%)")
    print(f"    MAE:  {mae:.3f}")
    print(f"    RMSE: {rmse:.3f}")
    print(f"    Pearson (expressed):  {r_p:.6f}")
    print(f"    Spearman (expressed): {r_s:.6f}")

# ── 2. False positive analysis (unexpressed transcripts) ──────
print("\n\n## 2. False Positives (unexpressed transcripts predicted > 0)")
unexpressed = oracle_df[oracle_df["mrna_truth"] == 0]
for tool in ["rigel_oracle", "salmon", "kallisto"]:
    fp = unexpressed[unexpressed[tool] > 0]
    fp_total = unexpressed[tool].sum()
    print(f"\n  {tool}:")
    print(f"    Transcripts with FP > 0: {len(fp):,} / {len(unexpressed):,}")
    print(f"    Total FP fragments: {fp_total:,.1f}")
    print(f"    Max single FP: {unexpressed[tool].max():.1f}")
    if len(fp) > 0:
        print(f"    Top 5 FP transcripts:")
        top_fp = fp.nlargest(5, tool)
        for _, row in top_fp.iterrows():
            print(f"      {row['transcript_id']} ({row['gene_name']}): {row[tool]:.1f}")

# ── 3. Worst absolute errors (expressed transcripts) ─────────
print("\n\n## 3. Worst Absolute Errors (expressed transcripts)")
for tool in ["rigel_oracle", "salmon", "kallisto"]:
    expressed_copy = expressed.copy()
    expressed_copy[f"{tool}_err"] = expressed_copy[tool] - expressed_copy["mrna_truth"]
    expressed_copy[f"{tool}_abs_err"] = expressed_copy[f"{tool}_err"].abs()
    expressed_copy[f"{tool}_rel_err"] = (
        expressed_copy[f"{tool}_abs_err"] / expressed_copy["mrna_truth"].clip(lower=1)
    )

    print(f"\n  {tool} — Top 10 by absolute error:")
    worst = expressed_copy.nlargest(10, f"{tool}_abs_err")
    for _, row in worst.iterrows():
        print(f"    {row['transcript_id']:25s}  truth={row['mrna_truth']:8.0f}  "
              f"pred={row[tool]:8.1f}  err={row[f'{tool}_err']:+8.1f}  "
              f"rel={row[f'{tool}_rel_err']:+.3f}")

# ── 4. Systematic bias: over vs under estimation ─────────────
print("\n\n## 4. Systematic Bias (expressed transcripts)")
for tool in ["rigel_oracle", "salmon", "kallisto"]:
    err = expressed[tool] - expressed["mrna_truth"]
    n_over = (err > 0.5).sum()
    n_under = (err < -0.5).sum()
    n_accurate = ((err >= -0.5) & (err <= 0.5)).sum()
    sum_over = err[err > 0].sum()
    sum_under = err[err < 0].sum()
    print(f"\n  {tool}:")
    print(f"    Overestimated: {n_over:,} ({100*n_over/len(expressed):.1f}%)")
    print(f"    Underestimated: {n_under:,} ({100*n_under/len(expressed):.1f}%)")
    print(f"    Accurate (±0.5): {n_accurate:,} ({100*n_accurate/len(expressed):.1f}%)")
    print(f"    Sum overestimation: {sum_over:+,.1f}")
    print(f"    Sum underestimation: {sum_under:+,.1f}")

# ── 5. Error by expression level bins ────────────────────────
print("\n\n## 5. Error by Expression Level (expressed transcripts)")
bins = [0, 1, 10, 100, 1000, 10000, np.inf]
labels = ["(0,1]", "(1,10]", "(10,100]", "(100,1K]", "(1K,10K]", ">10K"]
expressed_copy = expressed.copy()
expressed_copy["bin"] = pd.cut(expressed_copy["mrna_truth"], bins=bins, labels=labels)

for tool in ["rigel_oracle", "salmon", "kallisto"]:
    print(f"\n  {tool}:")
    print(f"    {'Bin':12s} {'Count':>7s} {'MAE':>10s} {'MAPE%':>9s} {'Bias':>10s}")
    for lbl in labels:
        sub = expressed_copy[expressed_copy["bin"] == lbl]
        if len(sub) == 0:
            continue
        err = sub[tool] - sub["mrna_truth"]
        mae = err.abs().mean()
        mape = (err.abs() / sub["mrna_truth"].clip(lower=0.01)).mean() * 100
        bias = err.mean()
        print(f"    {lbl:12s} {len(sub):7d} {mae:10.3f} {mape:8.1f}% {bias:+10.3f}")

# ── 6. Rigel-specific: nRNA and gDNA attribution ─────────────
print("\n\n## 6. Rigel nRNA/gDNA Attribution (should be ~0 in pristine data)")
print("\n  Oracle alignment:")
rigel_oracle_total = oracle_df["rigel_oracle"].sum()
truth_total = oracle_df["mrna_truth"].sum()
nrna_leak = truth_total - rigel_oracle_total  # fragments attributed to nRNA/gDNA
print(f"    mRNA (truth): {truth_total:,.0f}")
print(f"    mRNA (rigel): {rigel_oracle_total:,.1f}")
print(f"    Fragments leaked to nRNA/gDNA: {nrna_leak:,.1f} ({100*nrna_leak/truth_total:.3f}%)")

# ── 7. Minimap2 alignment analysis ───────────────────────────
print("\n\n## 7. Minimap2 Alignment Effects")
for tool in ["rigel_minimap2"]:
    pred = mm2_df[tool].sum()
    mae = (mm2_df[tool] - mm2_df["mrna_truth"]).abs().mean()
    rmse = np.sqrt(((mm2_df[tool] - mm2_df["mrna_truth"]) ** 2).mean())
    mask = mm2_df["mrna_truth"] > 0
    r_p, _ = pearsonr(mm2_df.loc[mask, "mrna_truth"], mm2_df.loc[mask, tool])
    r_s, _ = spearmanr(mm2_df.loc[mask, "mrna_truth"], mm2_df.loc[mask, tool])
    print(f"\n  {tool}:")
    print(f"    Total predicted: {pred:,.1f}")
    print(f"    Total − truth:   {pred - total_truth:+,.1f} ({100*(pred - total_truth)/total_truth:+.3f}%)")
    print(f"    MAE:  {mae:.3f}")
    print(f"    RMSE: {rmse:.3f}")
    print(f"    Pearson (expressed):  {r_p:.6f}")
    print(f"    Spearman (expressed): {r_s:.6f}")

# ── 8. Head-to-head: rigel_oracle vs salmon vs kallisto ──────
print("\n\n## 8. Head-to-Head Comparison (Oracle)")
for tool_a, tool_b in [("rigel_oracle", "salmon"), ("rigel_oracle", "kallisto"), ("salmon", "kallisto")]:
    err_a = (expressed[tool_a] - expressed["mrna_truth"]).abs()
    err_b = (expressed[tool_b] - expressed["mrna_truth"]).abs()
    a_better = (err_a < err_b - 0.5).sum()
    b_better = (err_b < err_a - 0.5).sum()
    tied = len(expressed) - a_better - b_better
    print(f"\n  {tool_a} vs {tool_b} (expressed transcripts):")
    print(f"    {tool_a} wins: {a_better:,} ({100*a_better/len(expressed):.1f}%)")
    print(f"    {tool_b} wins: {b_better:,} ({100*b_better/len(expressed):.1f}%)")
    print(f"    Tied (within ±0.5): {tied:,} ({100*tied/len(expressed):.1f}%)")

# ── 9. Quantile error analysis ───────────────────────────────
print("\n\n## 9. Error Distribution (expressed transcripts)")
for tool in ["rigel_oracle", "salmon", "kallisto"]:
    err = (expressed[tool] - expressed["mrna_truth"]).abs()
    print(f"\n  {tool} absolute error quantiles:")
    for q in [0.5, 0.75, 0.90, 0.95, 0.99, 1.0]:
        print(f"    {q*100:5.1f}%ile: {err.quantile(q):10.2f}")

# ── 10. Worst rigel-specific errors (where rigel is much worse than salmon) ──
print("\n\n## 10. Rigel-Specific Failures (oracle, where |rigel_err| > 2*|salmon_err| and |rigel_err| > 5)")
rigel_err = (expressed["rigel_oracle"] - expressed["mrna_truth"]).abs()
salmon_err = (expressed["salmon"] - expressed["mrna_truth"]).abs()
rigel_worse = expressed[(rigel_err > 2 * salmon_err) & (rigel_err > 5)].copy()
rigel_worse["rigel_abs_err"] = rigel_err.loc[rigel_worse.index]
rigel_worse["salmon_abs_err"] = salmon_err.loc[rigel_worse.index]
rigel_worse = rigel_worse.sort_values("rigel_abs_err", ascending=False)
print(f"  Count: {len(rigel_worse):,} transcripts")
if len(rigel_worse) > 0:
    print(f"\n  Top 20:")
    for i, (_, row) in enumerate(rigel_worse.head(20).iterrows()):
        print(f"    {row['transcript_id']:25s} gene={row['gene_name']:20s}  "
              f"truth={row['mrna_truth']:8.0f}  rigel={row['rigel_oracle']:8.1f}  "
              f"salmon={row['salmon']:8.1f}  r_err={row['rigel_abs_err']:+8.1f}  "
              f"s_err={row['salmon_abs_err']:+8.1f}")

print("\n\n## 11. Salmon-Specific Failures (oracle, where |salmon_err| > 2*|rigel_err| and |salmon_err| > 5)")
salmon_worse = expressed[(salmon_err > 2 * rigel_err) & (salmon_err > 5)].copy()
salmon_worse["rigel_abs_err"] = rigel_err.loc[salmon_worse.index]
salmon_worse["salmon_abs_err"] = salmon_err.loc[salmon_worse.index]
salmon_worse = salmon_worse.sort_values("salmon_abs_err", ascending=False)
print(f"  Count: {len(salmon_worse):,} transcripts")
if len(salmon_worse) > 0:
    print(f"\n  Top 20:")
    for i, (_, row) in enumerate(salmon_worse.head(20).iterrows()):
        print(f"    {row['transcript_id']:25s} gene={row['gene_name']:20s}  "
              f"truth={row['mrna_truth']:8.0f}  rigel={row['rigel_oracle']:8.1f}  "
              f"salmon={row['salmon']:8.1f}  r_err={row['rigel_abs_err']:+8.1f}  "
              f"s_err={row['salmon_abs_err']:+8.1f}")

print("\n" + "=" * 80)
print("END OF ANALYSIS")
print("=" * 80)
