#!/usr/bin/env python3
"""Deep analysis of VCaP benchmark results - VBEM vs MAP-EM comparison.

Investigates why VBEM performs significantly worse than MAP-EM:
- Pearson R: 0.8145 (VBEM) vs 0.9862 (MAP)
- WARE: 0.2516 (VBEM) vs 0.0796 (MAP)
- Very-high expression Pearson R: 0.34 (VBEM) vs 0.91 (MAP)
"""

import pandas as pd
import numpy as np
import json
import sys
from pathlib import Path


def load_data():
    """Load per-transcript detail data."""
    detail = pd.read_parquet("results/vcap/per_transcript_detail.parquet")
    print(f"Per-transcript detail: {len(detail)} rows")
    print(f"Columns: {list(detail.columns)}")
    print(f"Tools: {detail['tool'].unique()}")
    return detail


def analyze_vbem_vs_map(detail):
    """Compare VBEM and MAP transcript-by-transcript."""
    vbem = detail[detail["tool"] == "rigel/vbem"].copy()
    map_ = detail[detail["tool"] == "rigel/map"].copy()

    # Use actual column names: 'predicted', 'truth_tpm', 'residual', 'abs_error', 'rel_error'
    # Rename for clarity
    vbem = vbem.rename(columns={"predicted": "pred_tpm"})
    map_ = map_.rename(columns={"predicted": "pred_tpm"})

    print(f"\n{'='*80}")
    print("VBEM vs MAP-EM COMPARISON")
    print(f"{'='*80}")
    print(f"VBEM transcripts: {len(vbem)}")
    print(f"MAP transcripts: {len(map_)}")

    # Compute log2 values (with pseudocount)
    pc = 0.1
    for df in [vbem, map_]:
        df["log2_pred"] = np.log2(df["pred_tpm"] + pc)
        df["log2_truth"] = np.log2(df["truth_tpm"] + pc)
        df["log2_error"] = df["log2_pred"] - df["log2_truth"]

    # Merge on transcript_id
    merged = vbem.merge(
        map_[["transcript_id", "pred_tpm", "log2_pred", "log2_truth", "log2_error"]],
        on="transcript_id",
        suffixes=("_vbem", "_map"),
    )
    print(f"Merged rows: {len(merged)}")

    # Basic stats
    print(f"\nVBEM pred_tpm stats:")
    print(merged["pred_tpm_vbem"].describe())
    print(f"\nMAP pred_tpm stats:")
    print(merged["pred_tpm_map"].describe())

    # Compute abs error difference
    merged["vbem_abs_error"] = np.abs(merged["log2_error_vbem"])
    merged["map_abs_error"] = np.abs(merged["log2_error_map"])
    merged["error_diff"] = merged["vbem_abs_error"] - merged["map_abs_error"]

    # How many transcripts does VBEM do worse on?
    vbem_worse = (merged["error_diff"] > 0.1).sum()
    vbem_better = (merged["error_diff"] < -0.1).sum()
    similar = len(merged) - vbem_worse - vbem_better
    print(f"\n--- Error comparison (log2 scale, |diff| > 0.1) ---")
    print(f"VBEM worse: {vbem_worse} ({100*vbem_worse/len(merged):.1f}%)")
    print(f"VBEM better: {vbem_better} ({100*vbem_better/len(merged):.1f}%)")
    print(f"Similar: {similar} ({100*similar/len(merged):.1f}%)")

    # Worst VBEM transcripts (biggest error_diff)
    print(f"\n--- Top 30 transcripts where VBEM is worst vs MAP ---")
    worst = merged.nlargest(30, "error_diff")
    for _, row in worst.iterrows():
        print(
            f"  {row['transcript_id']}: "
            f"truth={row['truth_tpm']:.2f} TPM, "
            f"vbem={row['pred_tpm_vbem']:.2f}, "
            f"map={row['pred_tpm_map']:.2f}, "
            f"log2err vbem={row['log2_error_vbem']:.2f}, "
            f"map={row['log2_error_map']:.2f}"
        )

    # Analyze by expression level
    print(f"\n--- Error by expression bin ---")
    bins = [0, 1, 10, 100, 1000, np.inf]
    labels = ["zero/low (0-1)", "low (1-10)", "mid (10-100)", "high (100-1000)", "very_high (1000+)"]
    merged["expr_bin"] = pd.cut(merged["truth_tpm"], bins=bins, labels=labels, right=False)
    for bin_label in labels:
        subset = merged[merged["expr_bin"] == bin_label]
        if len(subset) == 0:
            continue
        vbem_mae = subset["vbem_abs_error"].mean()
        map_mae = subset["map_abs_error"].mean()
        vbem_worse_n = (subset["error_diff"] > 0.1).sum()
        print(
            f"  {bin_label} (n={len(subset)}): "
            f"VBEM mean|log2err|={vbem_mae:.3f}, "
            f"MAP mean|log2err|={map_mae:.3f}, "
            f"VBEM worse in {vbem_worse_n}/{len(subset)}"
        )

    return merged


def analyze_vbem_sparsification(detail):
    """Check if VBEM is over-sparsifying (zero-forcing too aggressively)."""
    vbem = detail[detail["tool"] == "rigel/vbem"].copy()
    map_ = detail[detail["tool"] == "rigel/map"].copy()

    print(f"\n{'='*80}")
    print("VBEM SPARSIFICATION ANALYSIS")
    print(f"{'='*80}")

    # Check zero predictions
    vbem_zeros = (vbem["predicted"] == 0).sum()
    map_zeros = (map_["predicted"] == 0).sum()
    vbem_near_zero = (vbem["predicted"] < 0.01).sum()
    map_near_zero = (map_["predicted"] < 0.01).sum()

    print(f"VBEM zero predictions: {vbem_zeros}")
    print(f"MAP zero predictions: {map_zeros}")
    print(f"VBEM near-zero (<0.01 TPM): {vbem_near_zero}")
    print(f"MAP near-zero (<0.01 TPM): {map_near_zero}")

    # True positives that got zeroed
    vbem_expressed = vbem[vbem["truth_tpm"] > 1.0]
    map_expressed = map_[map_["truth_tpm"] > 1.0]
    vbem_false_zero = (vbem_expressed["predicted"] < 0.01).sum()
    map_false_zero = (map_expressed["predicted"] < 0.01).sum()
    print(f"\nFalse zeros (truth > 1 TPM, pred < 0.01):")
    print(f"  VBEM: {vbem_false_zero} / {len(vbem_expressed)}")
    print(f"  MAP: {map_false_zero} / {len(map_expressed)}")

    # Check if VBEM is concentrating mass on fewer transcripts
    vbem_active = (vbem["predicted"] > 0.1).sum()
    map_active = (map_["predicted"] > 0.1).sum()
    print(f"\nActive transcripts (> 0.1 TPM):")
    print(f"  VBEM: {vbem_active}")
    print(f"  MAP: {map_active}")

    # Gini coefficient of predictions
    def gini(x):
        x = np.sort(x[x > 0])
        if len(x) == 0:
            return 0
        n = len(x)
        index = np.arange(1, n + 1)
        return (2 * np.sum(index * x) / (n * np.sum(x))) - (n + 1) / n

    print(f"\nGini coefficient (higher = more concentrated):")
    print(f"  VBEM: {gini(vbem['predicted'].values):.4f}")
    print(f"  MAP: {gini(map_['predicted'].values):.4f}")
    print(f"  Truth: {gini(vbem['truth_tpm'].values):.4f}")


def analyze_extreme_errors(detail):
    """Analyze the most extreme VBEM errors."""
    vbem = detail[detail["tool"] == "rigel/vbem"].copy()

    print(f"\n{'='*80}")
    print("EXTREME VBEM ERRORS")
    print(f"{'='*80}")

    # Compute log2 values
    pc = 0.1
    vbem["log2_pred"] = np.log2(vbem["predicted"] + pc)
    vbem["log2_truth"] = np.log2(vbem["truth_tpm"] + pc)
    vbem["log2_error"] = vbem["log2_pred"] - vbem["log2_truth"]

    # Large over-estimates
    vbem["log2_abs_error"] = np.abs(vbem["log2_error"])
    overestimates = vbem[vbem["log2_error"] > 2.0].nlargest(20, "log2_error")
    print(f"\n--- Top 20 VBEM over-estimates (log2_error > 2) ---")
    print(f"Total count: {(vbem['log2_error'] > 2.0).sum()}")
    for _, row in overestimates.iterrows():
        print(
            f"  {row['transcript_id']} (gene={row.get('gene_id', '?')}): "
            f"truth={row['truth_tpm']:.2f}, pred={row['predicted']:.2f}, "
            f"log2err={row['log2_error']:.2f}"
        )

    # Large under-estimates
    underestimates = vbem[vbem["log2_error"] < -2.0].nsmallest(20, "log2_error")
    print(f"\n--- Top 20 VBEM under-estimates (log2_error < -2) ---")
    print(f"Total count: {(vbem['log2_error'] < -2.0).sum()}")
    for _, row in underestimates.iterrows():
        print(
            f"  {row['transcript_id']} (gene={row.get('gene_id', '?')}): "
            f"truth={row['truth_tpm']:.2f}, pred={row['predicted']:.2f}, "
            f"log2err={row['log2_error']:.2f}"
        )


def analyze_mass_redistribution(detail):
    """Analyze how VBEM redistributes mass compared to MAP."""
    vbem = detail[detail["tool"] == "rigel/vbem"].copy()
    map_ = detail[detail["tool"] == "rigel/map"].copy()

    print(f"\n{'='*80}")
    print("MASS REDISTRIBUTION ANALYSIS")
    print(f"{'='*80}")

    # Check total TPM
    print(f"Total TPM - truth: {vbem['truth_tpm'].sum():.0f}")
    print(f"Total TPM - VBEM: {vbem['predicted'].sum():.0f}")
    print(f"Total TPM - MAP: {map_['predicted'].sum():.0f}")

    # Merge
    merged = vbem[["transcript_id", "truth_tpm", "predicted"]].merge(
        map_[["transcript_id", "predicted"]],
        on="transcript_id",
        suffixes=("_vbem", "_map"),
    )

    # Where does VBEM put more mass than MAP?
    merged["vbem_excess"] = merged["predicted_vbem"] - merged["predicted_map"]
    excess_positive = merged[merged["vbem_excess"] > 1.0]
    excess_negative = merged[merged["vbem_excess"] < -1.0]
    print(f"\nTranscripts where VBEM > MAP by >1 TPM: {len(excess_positive)}")
    print(f"  Total excess TPM: {excess_positive['vbem_excess'].sum():.0f}")
    print(f"Transcripts where VBEM < MAP by >1 TPM: {len(excess_negative)}")
    print(f"  Total deficit TPM: {excess_negative['vbem_excess'].sum():.0f}")

    # Top transcripts receiving excess mass from VBEM
    top_excess = excess_positive.nlargest(20, "vbem_excess")
    print(f"\n--- Top 20 transcripts where VBEM assigns MORE mass than MAP ---")
    for _, row in top_excess.iterrows():
        print(
            f"  {row['transcript_id']}: "
            f"truth={row['truth_tpm']:.2f}, "
            f"vbem={row['predicted_vbem']:.2f}, "
            f"map={row['predicted_map']:.2f}, "
            f"excess={row['vbem_excess']:.2f}"
        )

    # Top transcripts losing mass in VBEM
    top_deficit = excess_negative.nsmallest(20, "vbem_excess")
    print(f"\n--- Top 20 transcripts where VBEM assigns LESS mass than MAP ---")
    for _, row in top_deficit.iterrows():
        print(
            f"  {row['transcript_id']}: "
            f"truth={row['truth_tpm']:.2f}, "
            f"vbem={row['predicted_vbem']:.2f}, "
            f"map={row['predicted_map']:.2f}, "
            f"deficit={row['vbem_excess']:.2f}"
        )


def analyze_nrna_siphon_effect(detail):
    """Analyze whether nRNA siphon is worse for VBEM."""
    print(f"\n{'='*80}")
    print("nRNA SIPHON COMPARISON")
    print(f"{'='*80}")

    pool = pd.read_csv("results/vcap/pool_summary.csv")
    print(pool.to_string())

    # nRNA truth is 0 for nrna_none condition
    # Any nRNA prediction is pure siphon
    for _, row in pool.iterrows():
        tool = row["tool"]
        mrna_err = row.get("mrna_rel_error", "N/A")
        nrna_pred = row.get("nrna_pred", 0)
        gdna_err = row.get("gdna_rel_error", "N/A")
        print(f"\n{tool}:")
        print(f"  mRNA rel error: {mrna_err}")
        print(f"  nRNA siphon (should be 0): {nrna_pred:.0f}")
        print(f"  gDNA rel error: {gdna_err}")


def analyze_correlation_breakdown(detail):
    """Analyze at what expression levels VBEM correlation breaks down."""
    vbem = detail[detail["tool"] == "rigel/vbem"].copy()
    map_ = detail[detail["tool"] == "rigel/map"].copy()

    print(f"\n{'='*80}")
    print("CORRELATION BREAKDOWN BY TPM QUANTILE")
    print(f"{'='*80}")

    # Use truth TPM quantiles
    expressed_vbem = vbem[vbem["truth_tpm"] > 0].copy()
    expressed_map = map_[map_["truth_tpm"] > 0].copy()

    quantiles = [0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0]
    thresholds = expressed_vbem["truth_tpm"].quantile(quantiles)

    for i in range(len(quantiles) - 1):
        lo = thresholds.iloc[i]
        hi = thresholds.iloc[i + 1]
        label = f"Q{quantiles[i]:.0%}-{quantiles[i+1]:.0%} (TPM {lo:.1f}-{hi:.1f})"

        v_sub = expressed_vbem[
            (expressed_vbem["truth_tpm"] >= lo) & (expressed_vbem["truth_tpm"] < hi)
        ]
        m_sub = expressed_map[
            (expressed_map["truth_tpm"] >= lo) & (expressed_map["truth_tpm"] < hi)
        ]

        if len(v_sub) < 5 or len(m_sub) < 5:
            continue

        from scipy.stats import pearsonr, spearmanr

        v_pearson = pearsonr(v_sub["truth_tpm"], v_sub["predicted"])[0]
        m_pearson = pearsonr(m_sub["truth_tpm"], m_sub["predicted"])[0]
        v_spearman = spearmanr(v_sub["truth_tpm"], v_sub["predicted"])[0]
        m_spearman = spearmanr(m_sub["truth_tpm"], m_sub["predicted"])[0]

        print(
            f"  {label} (n={len(v_sub)}): "
            f"VBEM r={v_pearson:.4f} ρ={v_spearman:.4f} | "
            f"MAP r={m_pearson:.4f} ρ={m_spearman:.4f}"
        )

    # Also check: top 100 most expressed
    print(f"\n--- Top 100 most expressed transcripts ---")
    top100 = expressed_vbem.nlargest(100, "truth_tpm")
    top100_map = expressed_map[expressed_map["transcript_id"].isin(top100["transcript_id"])]

    from scipy.stats import pearsonr
    v_r = pearsonr(top100["truth_tpm"], top100["predicted"])[0]
    m_r = pearsonr(top100_map["truth_tpm"], top100_map["predicted"])[0]
    print(f"  VBEM Pearson R: {v_r:.4f}")
    print(f"  MAP Pearson R: {m_r:.4f}")

    # Top 10 most expressed
    print(f"\n--- Top 10 most expressed transcripts ---")
    top10 = expressed_vbem.nlargest(10, "truth_tpm")
    top10_map = expressed_map[expressed_map["transcript_id"].isin(top10["transcript_id"])]
    for _, row in top10.iterrows():
        map_row = top10_map[top10_map["transcript_id"] == row["transcript_id"]]
        map_pred = map_row["predicted"].values[0] if len(map_row) > 0 else "N/A"
        print(
            f"  {row['transcript_id']}: "
            f"truth={row['truth_tpm']:.1f}, "
            f"vbem={row['predicted']:.1f}, "
            f"map={map_pred:.1f}"
        )


def main():
    detail = load_data()
    analyze_vbem_vs_map(detail)
    analyze_vbem_sparsification(detail)
    analyze_extreme_errors(detail)
    analyze_mass_redistribution(detail)
    analyze_nrna_siphon_effect(detail)
    analyze_correlation_breakdown(detail)


if __name__ == "__main__":
    main()
