#!/usr/bin/env python3
"""Deep-dive analysis of benchmark v5: rigel vs salmon vs kallisto.

Reads the summary.csv and per-transcript CSVs to produce:
1. Summary comparison tables (transcript-level and gene-level)
2. Per-condition breakdown by gDNA level and strand specificity
3. Expression-stratified accuracy (low/mid/high expressors)
4. Per-transcript error distributions and outlier identification
5. Root cause analysis of where each tool excels/fails
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
import sys

OUTDIR = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v5_all")


def load_summary():
    df = pd.read_csv(OUTDIR / "summary.csv")
    return df


def load_per_transcript(condition: str, aligner: str) -> pd.DataFrame:
    p = OUTDIR / condition / f"per_transcript_counts_{aligner}.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def section(title: str):
    bar = "=" * 80
    print(f"\n{bar}\n{title}\n{bar}\n")


def subsection(title: str):
    print(f"\n--- {title} ---\n")


# ═══════════════════════════════════════════════════════════
# 1. SUMMARY OVERVIEW
# ═══════════════════════════════════════════════════════════
def analyze_summary(df: pd.DataFrame):
    section("1. TRANSCRIPT-LEVEL SUMMARY (rigel_oracle vs salmon vs kallisto)")

    tx = df[df["level"] == "transcript"].copy()
    # Focus on tools that appear consistently
    tools_of_interest = ["rigel_oracle", "salmon", "kallisto", "rigel_minimap2"]
    tx = tx[tx["tool"].isin(tools_of_interest)]

    # Average across all conditions
    avg = (
        tx.groupby("tool")[["mean_abs_error", "rmse", "pearson", "spearman"]]
        .agg(["mean", "std"])
        .round(4)
    )
    print("Average across all 9 conditions:")
    print(avg.to_string())
    print()

    # By gDNA level
    subsection("By gDNA contamination level")
    for gdna in ["none", "low", "high"]:
        subset = tx[tx["gdna_label"] == gdna]
        avg_gdna = (
            subset.groupby("tool")[["mean_abs_error", "rmse", "pearson", "spearman"]]
            .mean()
            .round(4)
        )
        print(f"\ngDNA = {gdna}:")
        print(avg_gdna.to_string())

    # By strand specificity
    subsection("By strand specificity")
    for ss in [0.5, 0.9, 1.0]:
        subset = tx[tx["strand_specificity"] == ss]
        avg_ss = (
            subset.groupby("tool")[["mean_abs_error", "rmse", "pearson", "spearman"]]
            .mean()
            .round(4)
        )
        print(f"\nSS = {ss}:")
        print(avg_ss.to_string())

    section("2. GENE-LEVEL SUMMARY")
    gene = df[df["level"] == "gene"].copy()
    gene = gene[gene["tool"].isin(tools_of_interest)]
    avg_gene = (
        gene.groupby("tool")[["mean_abs_error", "rmse", "pearson", "spearman"]]
        .agg(["mean", "std"])
        .round(4)
    )
    print("Average across all 9 conditions:")
    print(avg_gene.to_string())


# ═══════════════════════════════════════════════════════════
# 2. POOL-LEVEL (mRNA/nRNA/gDNA totals)
# ═══════════════════════════════════════════════════════════
def analyze_pool_level(df: pd.DataFrame):
    section("3. POOL-LEVEL mRNA ACCURACY (rigel only)")

    tx = df[df["level"] == "transcript"].copy()
    rigel_tools = tx[tx["tool"].isin(["rigel_oracle", "rigel_minimap2"])]

    print(f"{'Condition':<45} {'Tool':<18} {'mRNA_truth':>12} {'mRNA_pred':>12} "
          f"{'Err%':>8} {'nRNA_truth':>12} {'nRNA_pred':>12} {'gDNA_truth':>10} {'gDNA_pred':>10}")
    print("-" * 160)

    for _, row in rigel_tools.iterrows():
        err_pct = (row["total_observed"] - row["total_truth"]) / row["total_truth"] * 100
        print(
            f"{row['dataset']:<45} {row['tool']:<18} "
            f"{row['total_truth']:>12.0f} {row['total_observed']:>12.0f} "
            f"{err_pct:>+7.2f}%"
        )


# ═══════════════════════════════════════════════════════════
# 3. EXPRESSION-STRATIFIED ANALYSIS
# ═══════════════════════════════════════════════════════════
def analyze_expression_strata(conditions=None):
    section("4. EXPRESSION-STRATIFIED ACCURACY")

    if conditions is None:
        conditions = [
            ("gdna_none_ss_0.90_nrna_default", "oracle"),
            ("gdna_low_ss_0.90_nrna_default", "oracle"),
            ("gdna_high_ss_0.90_nrna_default", "oracle"),
        ]

    for cond, aligner in conditions:
        ptx = load_per_transcript(cond, aligner)
        if ptx is None:
            print(f"  [SKIP] {cond}/{aligner}: no per-transcript CSV")
            continue

        subsection(f"{cond} / {aligner}")

        # Tools present
        tool_cols = [c for c in ptx.columns if c not in [
            "transcript_id", "gene_id", "gene_name",
            "mrna_abundance", "nrna_abundance", "mrna_truth", "nrna_truth"
        ]]
        print(f"Tools: {tool_cols}")

        truth_col = "mrna_truth"
        truth = ptx[truth_col].values.astype(float)

        # Define strata by expression quantiles (among expressed)
        expressed_mask = truth > 0
        n_expressed = expressed_mask.sum()
        n_total = len(truth)
        print(f"Transcripts: {n_total} total, {n_expressed} expressed")

        if n_expressed < 10:
            continue

        truth_expr = truth[expressed_mask]
        q25 = np.percentile(truth_expr, 25)
        q75 = np.percentile(truth_expr, 75)

        strata = {
            "low (0, Q25]": (truth > 0) & (truth <= q25),
            f"mid (Q25, Q75]": (truth > q25) & (truth <= q75),
            f"high (>Q75)": truth > q75,
            "zero (truth=0)": truth == 0,
        }
        print(f"Quantiles: Q25={q25:.1f}, Q75={q75:.1f}")

        for tool in tool_cols:
            obs = ptx[tool].values.astype(float)
            print(f"\n  Tool: {tool}")
            print(f"  {'Stratum':<25} {'N':>8} {'MAE':>10} {'MedAE':>10} {'RMSE':>10} {'Pearson':>10}")

            for name, mask in strata.items():
                n = mask.sum()
                if n == 0:
                    continue
                t = truth[mask]
                o = obs[mask]
                ae = np.abs(t - o)
                mae = ae.mean()
                medae = np.median(ae)
                rmse = np.sqrt(np.mean((t - o) ** 2))
                if t.std() > 0 and o.std() > 0:
                    pearson = np.corrcoef(t, o)[0, 1]
                else:
                    pearson = float("nan")
                print(f"  {name:<25} {n:>8} {mae:>10.2f} {medae:>10.2f} {rmse:>10.2f} {pearson:>10.4f}")


# ═══════════════════════════════════════════════════════════
# 4. PER-TRANSCRIPT ERROR DISTRIBUTION
# ═══════════════════════════════════════════════════════════
def analyze_error_distribution():
    section("5. ERROR DISTRIBUTION ANALYSIS (representative condition)")

    cond = "gdna_low_ss_0.90_nrna_default"
    aligner = "oracle"
    ptx = load_per_transcript(cond, aligner)
    if ptx is None:
        print("[SKIP] No per-transcript CSV")
        return

    truth = ptx["mrna_truth"].values.astype(float)
    tool_cols = [c for c in ptx.columns if c not in [
        "transcript_id", "gene_id", "gene_name",
        "mrna_abundance", "nrna_abundance", "mrna_truth", "nrna_truth"
    ]]

    expressed = truth > 0
    n_expr = expressed.sum()
    print(f"Condition: {cond} / {aligner}")
    print(f"Expressed transcripts: {n_expr}")
    print()

    for tool in tool_cols:
        obs = ptx[tool].values.astype(float)
        errors = obs[expressed] - truth[expressed]
        rel_errors = errors / truth[expressed]
        abs_errors = np.abs(errors)

        print(f"Tool: {tool}")
        print(f"  Absolute error percentiles: "
              f"P50={np.percentile(abs_errors, 50):.2f}, "
              f"P90={np.percentile(abs_errors, 90):.2f}, "
              f"P95={np.percentile(abs_errors, 95):.2f}, "
              f"P99={np.percentile(abs_errors, 99):.2f}")
        print(f"  Relative error percentiles: "
              f"P50={np.percentile(rel_errors, 50):.4f}, "
              f"P90={np.percentile(rel_errors, 90):.4f}, "
              f"P95={np.percentile(rel_errors, 95):.4f}, "
              f"P99={np.percentile(rel_errors, 99):.4f}")
        print(f"  Signed error: mean={errors.mean():.2f}, median={np.median(errors):.2f}")

        # False positives: predicted > 0 when truth = 0
        zero_truth = truth == 0
        fp = (obs[zero_truth] > 0.5).sum()
        fp_total = zero_truth.sum()
        print(f"  False positives: {fp}/{fp_total} ({100*fp/fp_total:.2f}% of zero-truth transcripts)")

        # False negatives: predicted ~0 when truth > threshold
        high_truth = truth > 10
        fn = (obs[high_truth] < 0.5).sum()
        fn_total = high_truth.sum()
        print(f"  False negatives (truth>10): {fn}/{fn_total} ({100*fn/fn_total:.2f}%)")
        print()


# ═══════════════════════════════════════════════════════════
# 5. HEAD-TO-HEAD: rigel_oracle vs salmon detailed
# ═══════════════════════════════════════════════════════════
def head_to_head_analysis():
    section("6. HEAD-TO-HEAD: rigel_oracle vs salmon vs kallisto")

    cond = "gdna_low_ss_0.90_nrna_default"
    aligner = "oracle"
    ptx = load_per_transcript(cond, aligner)
    if ptx is None:
        print("[SKIP]")
        return

    truth = ptx["mrna_truth"].values.astype(float)

    tools = ["rigel_oracle", "salmon", "kallisto"]
    available = [t for t in tools if t in ptx.columns]

    if len(available) < 2:
        print("Not enough tools for comparison")
        return

    expr_mask = truth > 0

    subsection(f"Condition: {cond} — expressed transcripts only (n={expr_mask.sum()})")

    # Compute per-tool absolute errors
    tool_errors = {}
    for tool in available:
        obs = ptx[tool].values.astype(float)
        tool_errors[tool] = np.abs(obs[expr_mask] - truth[expr_mask])

    # Pairwise: rigel_oracle beats salmon/kallisto on how many transcripts?
    if "rigel_oracle" in available:
        rigel_ae = tool_errors["rigel_oracle"]
        for other in ["salmon", "kallisto"]:
            if other not in available:
                continue
            other_ae = tool_errors[other]
            rigel_wins = (rigel_ae < other_ae).sum()
            other_wins = (other_ae < rigel_ae).sum()
            ties = (rigel_ae == other_ae).sum()
            total = len(rigel_ae)
            print(f"rigel_oracle vs {other}: "
                  f"rigel wins {rigel_wins}/{total} ({100*rigel_wins/total:.1f}%), "
                  f"{other} wins {other_wins}/{total} ({100*other_wins/total:.1f}%), "
                  f"ties {ties}/{total}")

    # Where does salmon beat rigel? (largest margin)
    if "rigel_oracle" in available and "salmon" in available:
        subsection("Transcripts where salmon beats rigel_oracle (top 20 by margin)")
        rigel_obs = ptx["rigel_oracle"].values.astype(float)[expr_mask]
        salmon_obs = ptx["salmon"].values.astype(float)[expr_mask]
        truth_e = truth[expr_mask]
        tx_ids = ptx["transcript_id"].values[expr_mask]
        gene_names = ptx["gene_name"].values[expr_mask]

        rigel_ae = np.abs(rigel_obs - truth_e)
        salmon_ae = np.abs(salmon_obs - truth_e)
        margin = rigel_ae - salmon_ae  # positive = salmon is better

        top_idx = np.argsort(margin)[::-1][:20]
        print(f"{'Transcript':<25} {'Gene':<15} {'Truth':>8} {'Rigel':>8} {'Salmon':>8} "
              f"{'Rigel_AE':>10} {'Salmon_AE':>10} {'Margin':>10}")
        for i in top_idx:
            print(f"{tx_ids[i]:<25} {gene_names[i]:<15} {truth_e[i]:>8.1f} {rigel_obs[i]:>8.1f} "
                  f"{salmon_obs[i]:>8.1f} {rigel_ae[i]:>10.1f} {salmon_ae[i]:>10.1f} {margin[i]:>10.1f}")

    # Where does rigel beat salmon? (largest margin)
    if "rigel_oracle" in available and "salmon" in available:
        subsection("Transcripts where rigel_oracle beats salmon (top 20 by margin)")
        bottom_idx = np.argsort(margin)[:20]
        print(f"{'Transcript':<25} {'Gene':<15} {'Truth':>8} {'Rigel':>8} {'Salmon':>8} "
              f"{'Rigel_AE':>10} {'Salmon_AE':>10} {'Margin':>10}")
        for i in bottom_idx:
            print(f"{tx_ids[i]:<25} {gene_names[i]:<15} {truth_e[i]:>8.1f} {rigel_obs[i]:>8.1f} "
                  f"{salmon_obs[i]:>8.1f} {rigel_ae[i]:>10.1f} {salmon_ae[i]:>10.1f} {margin[i]:>10.1f}")


# ═══════════════════════════════════════════════════════════
# 6. nRNA IMPACT ANALYSIS
# ═══════════════════════════════════════════════════════════
def nrna_impact_analysis():
    section("7. nRNA IMPACT: How does nRNA abundance affect tool accuracy?")

    cond = "gdna_none_ss_0.90_nrna_default"
    aligner = "oracle"
    ptx = load_per_transcript(cond, aligner)
    if ptx is None:
        print("[SKIP]")
        return

    truth = ptx["mrna_truth"].values.astype(float)
    nrna_truth = ptx["nrna_truth"].values.astype(float)
    expr_mask = truth > 0

    # Compute nrna:mrna ratio for expressed transcripts
    nrna_ratio = np.zeros_like(truth)
    nrna_ratio[expr_mask] = nrna_truth[expr_mask] / truth[expr_mask]

    tools = ["rigel_oracle", "salmon", "kallisto"]
    available = [t for t in tools if t in ptx.columns]

    # Stratify by nRNA ratio
    ratio_bins = [
        ("no nRNA (ratio=0)", nrna_ratio == 0),
        ("low nRNA (0 < ratio ≤ 1)", (nrna_ratio > 0) & (nrna_ratio <= 1)),
        ("mid nRNA (1 < ratio ≤ 5)", (nrna_ratio > 1) & (nrna_ratio <= 5)),
        ("high nRNA (ratio > 5)", nrna_ratio > 5),
    ]

    subsection(f"Condition: {cond}")
    print(f"{'Stratum':<30} {'N':>6}", end="")
    for tool in available:
        print(f"  {tool:>16}_MAE", end="")
    print()

    for name, mask in ratio_bins:
        combined = expr_mask & mask
        n = combined.sum()
        if n == 0:
            continue
        t = truth[combined]
        print(f"{name:<30} {n:>6}", end="")
        for tool in available:
            obs = ptx[tool].values.astype(float)[combined]
            mae = np.abs(obs - t).mean()
            print(f"  {mae:>20.2f}", end="")
        print()

    # Also show relative error
    subsection("Relative error (% of truth) by nRNA ratio stratum")
    print(f"{'Stratum':<30} {'N':>6}", end="")
    for tool in available:
        print(f"  {tool:>16}_RE%", end="")
    print()

    for name, mask in ratio_bins:
        combined = expr_mask & mask
        n = combined.sum()
        if n == 0:
            continue
        t = truth[combined]
        print(f"{name:<30} {n:>6}", end="")
        for tool in available:
            obs = ptx[tool].values.astype(float)[combined]
            re = np.abs(obs - t) / t
            median_re = np.median(re) * 100
            print(f"  {median_re:>20.1f}", end="")
        print()


# ═══════════════════════════════════════════════════════════
# 7. ISOFORM RESOLUTION ANALYSIS
# ═══════════════════════════════════════════════════════════
def isoform_analysis():
    section("8. ISOFORM RESOLUTION: Multi-isoform gene accuracy")

    cond = "gdna_none_ss_0.90_nrna_default"
    aligner = "oracle"
    ptx = load_per_transcript(cond, aligner)
    if ptx is None:
        print("[SKIP]")
        return

    truth = ptx["mrna_truth"].values.astype(float)
    tools = ["rigel_oracle", "salmon", "kallisto"]
    available = [t for t in tools if t in ptx.columns]

    # Count isoforms per gene
    gene_isoform_counts = ptx.groupby("gene_id").size()
    ptx_copy = ptx.copy()
    ptx_copy["n_isoforms"] = ptx_copy["gene_id"].map(gene_isoform_counts)

    # Strata: single-isoform vs multi-isoform
    iso_bins = [
        ("1 isoform", ptx_copy["n_isoforms"] == 1),
        ("2 isoforms", ptx_copy["n_isoforms"] == 2),
        ("3-5 isoforms", (ptx_copy["n_isoforms"] >= 3) & (ptx_copy["n_isoforms"] <= 5)),
        ("6-10 isoforms", (ptx_copy["n_isoforms"] >= 6) & (ptx_copy["n_isoforms"] <= 10)),
        (">10 isoforms", ptx_copy["n_isoforms"] > 10),
    ]

    expr_mask = truth > 0

    print(f"{'Isoform stratum':<20} {'N_expr':>8}", end="")
    for tool in available:
        print(f"  {tool:>16}_MAE", end="")
    print()

    for name, mask in iso_bins:
        combined = expr_mask & mask.values
        n = combined.sum()
        if n == 0:
            continue
        t = truth[combined]
        print(f"{name:<20} {n:>8}", end="")
        for tool in available:
            obs = ptx[tool].values.astype(float)[combined]
            mae = np.abs(obs - t).mean()
            print(f"  {mae:>20.2f}", end="")
        print()


# ═══════════════════════════════════════════════════════════
# 8. CONDITION SENSITIVITY HEATMAP
# ═══════════════════════════════════════════════════════════
def condition_heatmap(df: pd.DataFrame):
    section("9. CONDITION SENSITIVITY: MAE by (gDNA, SS) for each tool")

    tx = df[(df["level"] == "transcript")].copy()
    tools = ["rigel_oracle", "salmon", "kallisto", "rigel_minimap2"]
    tx = tx[tx["tool"].isin(tools)]

    for tool in tools:
        subset = tx[tx["tool"] == tool]
        pivot = subset.pivot_table(
            values="mean_abs_error",
            index="gdna_label",
            columns="strand_specificity",
            aggfunc="mean"
        )
        # Reorder rows
        pivot = pivot.reindex(["none", "low", "high"])
        print(f"\n{tool} — MAE by (gDNA, SS):")
        print(pivot.round(2).to_string())

    subsection("Pearson correlation by (gDNA, SS)")
    for tool in tools:
        subset = tx[tx["tool"] == tool]
        pivot = subset.pivot_table(
            values="pearson",
            index="gdna_label",
            columns="strand_specificity",
            aggfunc="mean"
        )
        pivot = pivot.reindex(["none", "low", "high"])
        print(f"\n{tool} — Pearson by (gDNA, SS):")
        print(pivot.round(4).to_string())


# ═══════════════════════════════════════════════════════════
# 9. COMPETITIVE SUMMARY
# ═══════════════════════════════════════════════════════════
def competitive_summary(df: pd.DataFrame):
    section("10. COMPETITIVE SUMMARY")

    tx = df[df["level"] == "transcript"].copy()

    # For each condition, rank tools by MAE
    conditions = tx["dataset"].unique()
    tool_ranks = []
    for cond in conditions:
        cond_data = tx[tx["dataset"] == cond]
        ranked = cond_data.sort_values("mean_abs_error")
        for rank, (_, row) in enumerate(ranked.iterrows(), 1):
            tool_ranks.append({"tool": row["tool"], "rank": rank, "condition": cond})

    ranks_df = pd.DataFrame(tool_ranks)
    avg_rank = ranks_df.groupby("tool")["rank"].mean().sort_values()
    print("Average rank across all conditions (lower is better):")
    print(avg_rank.round(2).to_string())
    print()

    # Win counts
    wins = ranks_df[ranks_df["rank"] == 1].groupby("tool").size()
    print("Number of conditions where each tool ranks #1:")
    print(wins.to_string())
    print()

    # Speed comparison
    subsection("Speed comparison (average quant time across conditions)")
    tools = ["rigel_oracle", "salmon", "kallisto", "rigel_minimap2"]
    for tool in tools:
        subset = tx[tx["tool"] == tool]
        avg_time = subset["elapsed_sec"].mean()
        avg_rss = subset["peak_rss_mb"].mean()
        print(f"  {tool:<20}: {avg_time:.1f}s avg, {avg_rss:.0f} MB RSS")


# ═══════════════════════════════════════════════════════════
# 10. ROOT CAUSE: Where salmon beats rigel
# ═══════════════════════════════════════════════════════════
def root_cause_salmon_advantage():
    section("11. ROOT CAUSE: Conditions where salmon/kallisto approach or beat rigel")

    df = load_summary()
    tx = df[df["level"] == "transcript"].copy()

    # For each condition, compare rigel_oracle MAE vs salmon MAE
    conditions = tx["dataset"].unique()

    print(f"{'Condition':<50} {'rigel_MAE':>10} {'salmon_MAE':>10} {'kall_MAE':>10} {'rigel/salmon':>12}")
    print("-" * 100)
    for cond in sorted(conditions):
        cond_data = tx[tx["dataset"] == cond]
        rigel_row = cond_data[cond_data["tool"] == "rigel_oracle"]
        salmon_row = cond_data[cond_data["tool"] == "salmon"]
        kallisto_row = cond_data[cond_data["tool"] == "kallisto"]

        if rigel_row.empty or salmon_row.empty:
            continue

        r_mae = rigel_row["mean_abs_error"].values[0]
        s_mae = salmon_row["mean_abs_error"].values[0]
        k_mae = kallisto_row["mean_abs_error"].values[0] if not kallisto_row.empty else float("nan")
        ratio = r_mae / s_mae

        print(f"{cond:<50} {r_mae:>10.2f} {s_mae:>10.2f} {k_mae:>10.2f} {ratio:>12.4f}")


# ═══════════════════════════════════════════════════════════
# 11. SALMON GENE-LEVEL ANALYSIS
# ═══════════════════════════════════════════════════════════
def gene_level_salmon_advantage():
    section("12. GENE-LEVEL: rigel vs salmon vs kallisto")

    df = load_summary()
    gene = df[df["level"] == "gene"].copy()
    tools = ["rigel_oracle", "salmon", "kallisto"]
    gene = gene[gene["tool"].isin(tools)]

    avg = gene.groupby("tool")[["mean_abs_error", "rmse", "pearson", "spearman"]].mean().round(4)
    print("Average gene-level metrics across all conditions:")
    print(avg.to_string())
    print()

    conditions = gene["dataset"].unique()
    print(f"\n{'Condition':<50} {'rigel_MAE':>10} {'salmon_MAE':>10} {'kall_MAE':>10}")
    print("-" * 90)
    for cond in sorted(conditions):
        cond_data = gene[gene["dataset"] == cond]
        r = cond_data[cond_data["tool"] == "rigel_oracle"]["mean_abs_error"].values
        s = cond_data[cond_data["tool"] == "salmon"]["mean_abs_error"].values
        k = cond_data[cond_data["tool"] == "kallisto"]["mean_abs_error"].values
        r_mae = r[0] if len(r) > 0 else float("nan")
        s_mae = s[0] if len(s) > 0 else float("nan")
        k_mae = k[0] if len(k) > 0 else float("nan")
        print(f"{cond:<50} {r_mae:>10.2f} {s_mae:>10.2f} {k_mae:>10.2f}")


# ═══════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════
def main():
    print("=" * 80)
    print("BENCHMARK V5 DEEP-DIVE ANALYSIS")
    print(f"Data: {OUTDIR}")
    print("=" * 80)

    df = load_summary()
    analyze_summary(df)
    analyze_pool_level(df)
    analyze_expression_strata()
    analyze_error_distribution()
    head_to_head_analysis()
    nrna_impact_analysis()
    isoform_analysis()
    condition_heatmap(df)
    competitive_summary(df)
    root_cause_salmon_advantage()
    gene_level_salmon_advantage()

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
