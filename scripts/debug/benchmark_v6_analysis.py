#!/usr/bin/env python3
"""Deep-dive analysis of benchmark v6: pruning ON/OFF × oracle/minimap2 vs salmon/kallisto.

Reads the summary.csv and per-transcript CSVs to produce:
 1. Summary comparison tables (transcript-level and gene-level)
 2. Pruning impact analysis (default vs no_prune, FN/FP rates)
 3. Realistic comparison: minimap2 rigel vs salmon vs kallisto
 4. Oracle comparison: upper bound accuracy
 5. Per-condition breakdown by gDNA level and strand specificity
 6. Expression-stratified accuracy (low/mid/high expressors)
 7. Error distribution analysis
 8. nRNA impact analysis
 9. Isoform resolution analysis
10. Head-to-head per-transcript comparisons
11. Pool-level mRNA/nRNA/gDNA accuracy (rigel only)
12. Performance (speed/memory)
13. Competitive ranking
"""

import pandas as pd
import numpy as np
from pathlib import Path

OUTDIR = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v6_prune")

# Expected tool names (rigel_{config}_{aligner} when multi_rigel=True)
RIGEL_ORACLE_DEFAULT = "rigel_default_oracle"
RIGEL_ORACLE_NOPRUNE = "rigel_no_prune_oracle"
RIGEL_MM2_DEFAULT = "rigel_default_minimap2"
RIGEL_MM2_NOPRUNE = "rigel_no_prune_minimap2"
SALMON = "salmon"
KALLISTO = "kallisto"

ALL_TOOLS = [
    RIGEL_ORACLE_DEFAULT, RIGEL_ORACLE_NOPRUNE,
    RIGEL_MM2_DEFAULT, RIGEL_MM2_NOPRUNE,
    SALMON, KALLISTO,
]

# Groupings for analysis
ORACLE_TOOLS = [RIGEL_ORACLE_DEFAULT, RIGEL_ORACLE_NOPRUNE]
MM2_TOOLS = [RIGEL_MM2_DEFAULT, RIGEL_MM2_NOPRUNE]
RIGEL_TOOLS = ORACLE_TOOLS + MM2_TOOLS
EXTERNAL_TOOLS = [SALMON, KALLISTO]
# "Realistic" comparison: minimap2 rigel variants + external tools
REALISTIC_TOOLS = MM2_TOOLS + EXTERNAL_TOOLS

# Short display names
SHORT = {
    RIGEL_ORACLE_DEFAULT: "rigel_orc_def",
    RIGEL_ORACLE_NOPRUNE: "rigel_orc_nop",
    RIGEL_MM2_DEFAULT: "rigel_mm2_def",
    RIGEL_MM2_NOPRUNE: "rigel_mm2_nop",
    SALMON: "salmon",
    KALLISTO: "kallisto",
}


def load_summary():
    return pd.read_csv(OUTDIR / "summary.csv")


def load_per_transcript(condition: str, aligner: str) -> pd.DataFrame | None:
    p = OUTDIR / condition / f"per_transcript_counts_{aligner}.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def section(title: str):
    bar = "=" * 90
    print(f"\n{bar}\n{title}\n{bar}\n")


def subsection(title: str):
    print(f"\n--- {title} ---\n")


def tool_cols(ptx: pd.DataFrame) -> list[str]:
    """Get tool column names from per-transcript CSV."""
    meta = {
        "transcript_id", "gene_id", "gene_name",
        "mrna_abundance", "nrna_abundance", "mrna_truth", "nrna_truth",
    }
    return [c for c in ptx.columns if c not in meta]


def compute_fn_fp(truth: np.ndarray, obs: np.ndarray, truth_threshold: float = 10.0):
    """Compute false negative and false positive rates."""
    # FP: predicted > 0.5 when truth = 0
    zero_mask = truth == 0
    n_zero = zero_mask.sum()
    fp = (obs[zero_mask] > 0.5).sum() if n_zero > 0 else 0
    fp_rate = fp / n_zero if n_zero > 0 else 0

    # FN: predicted < 0.5 when truth > threshold
    high_mask = truth > truth_threshold
    n_high = high_mask.sum()
    fn = (obs[high_mask] < 0.5).sum() if n_high > 0 else 0
    fn_rate = fn / n_high if n_high > 0 else 0

    return {
        "fp": fp, "fp_total": n_zero, "fp_rate": fp_rate,
        "fn": fn, "fn_total": n_high, "fn_rate": fn_rate,
    }


# ═══════════════════════════════════════════════════════════
# 1. SUMMARY OVERVIEW
# ═══════════════════════════════════════════════════════════
def analyze_summary(df: pd.DataFrame):
    section("1. TRANSCRIPT-LEVEL SUMMARY (all 6 tools)")

    tx = df[df["level"] == "transcript"].copy()
    present = [t for t in ALL_TOOLS if t in tx["tool"].unique()]
    tx = tx[tx["tool"].isin(present)]

    avg = (
        tx.groupby("tool")[["mean_abs_error", "rmse", "pearson", "spearman"]]
        .agg(["mean", "std"])
        .round(4)
    )
    print("Average across all 9 conditions:")
    print(avg.to_string())

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
    gene = gene[gene["tool"].isin(present)]
    avg_gene = (
        gene.groupby("tool")[["mean_abs_error", "rmse", "pearson", "spearman"]]
        .agg(["mean", "std"])
        .round(4)
    )
    print("Average across all 9 conditions:")
    print(avg_gene.to_string())


# ═══════════════════════════════════════════════════════════
# 2. PRUNING IMPACT (KEY NEW ANALYSIS)
# ═══════════════════════════════════════════════════════════
def analyze_pruning_impact(df: pd.DataFrame):
    section("3. PRUNING IMPACT: default (prune_threshold=0.1) vs no_prune")

    tx = df[df["level"] == "transcript"].copy()

    subsection("3a. Oracle: pruning ON vs OFF (transcript-level MAE)")
    pairs = [
        (RIGEL_ORACLE_DEFAULT, RIGEL_ORACLE_NOPRUNE, "oracle"),
        (RIGEL_MM2_DEFAULT, RIGEL_MM2_NOPRUNE, "minimap2"),
    ]
    for default_tool, noprune_tool, label in pairs:
        t_def = tx[tx["tool"] == default_tool]
        t_nop = tx[tx["tool"] == noprune_tool]
        if t_def.empty or t_nop.empty:
            print(f"  [{label}] SKIP — tool not found")
            continue

        print(f"\n  {label.upper()}: pruning ON vs OFF")
        print(f"  {'Condition':<50} {'prune_ON':>10} {'prune_OFF':>10} {'Δ':>10} {'Δ%':>8}")
        print("  " + "-" * 90)

        for _, row_def in t_def.iterrows():
            cond = row_def["dataset"]
            row_nop = t_nop[t_nop["dataset"] == cond]
            if row_nop.empty:
                continue
            mae_def = row_def["mean_abs_error"]
            mae_nop = row_nop["mean_abs_error"].values[0]
            delta = mae_nop - mae_def
            pct = 100 * delta / mae_def if mae_def > 0 else 0
            print(f"  {cond:<50} {mae_def:>10.4f} {mae_nop:>10.4f} {delta:>+10.4f} {pct:>+7.2f}%")

        # Averages
        avg_def = t_def["mean_abs_error"].mean()
        avg_nop = t_nop["mean_abs_error"].mean()
        delta_avg = avg_nop - avg_def
        pct_avg = 100 * delta_avg / avg_def if avg_def > 0 else 0
        print(f"  {'AVERAGE':<50} {avg_def:>10.4f} {avg_nop:>10.4f} {delta_avg:>+10.4f} {pct_avg:>+7.2f}%")

    # Gene-level pruning impact
    subsection("3b. Gene-level pruning impact")
    gene = df[df["level"] == "gene"].copy()
    for default_tool, noprune_tool, label in pairs:
        g_def = gene[gene["tool"] == default_tool]
        g_nop = gene[gene["tool"] == noprune_tool]
        if g_def.empty or g_nop.empty:
            continue
        avg_def = g_def["mean_abs_error"].mean()
        avg_nop = g_nop["mean_abs_error"].mean()
        pct = 100 * (avg_nop - avg_def) / avg_def if avg_def > 0 else 0
        print(f"  {label}: gene MAE  prune_ON={avg_def:.4f}  prune_OFF={avg_nop:.4f}  Δ={pct:+.2f}%")

    # Pearson comparison
    subsection("3c. Pearson correlation: pruning ON vs OFF")
    for default_tool, noprune_tool, label in pairs:
        t_def = tx[tx["tool"] == default_tool]
        t_nop = tx[tx["tool"] == noprune_tool]
        if t_def.empty or t_nop.empty:
            continue
        avg_p_def = t_def["pearson"].mean()
        avg_p_nop = t_nop["pearson"].mean()
        print(f"  {label}: Pearson  prune_ON={avg_p_def:.6f}  prune_OFF={avg_p_nop:.6f}")


# ═══════════════════════════════════════════════════════════
# 3. FALSE NEGATIVE / FALSE POSITIVE ANALYSIS (pruning focus)
# ═══════════════════════════════════════════════════════════
def analyze_fn_fp_pruning():
    section("4. FALSE NEGATIVE / FALSE POSITIVE RATES: pruning ON vs OFF")

    conditions = sorted([
        d.name for d in OUTDIR.iterdir()
        if d.is_dir() and d.name.startswith("gdna_")
    ])

    for aligner in ["oracle", "minimap2"]:
        subsection(f"Aligner: {aligner}")
        print(f"  {'Condition':<45} {'Tool':<22} {'FP':>6} {'FP%':>8} "
              f"{'FN(>10)':>8} {'FN%':>8}")
        print("  " + "-" * 100)

        all_fn = {}  # tool -> list of fn counts
        all_fp = {}
        for cond in conditions:
            ptx = load_per_transcript(cond, aligner)
            if ptx is None:
                continue
            truth = ptx["mrna_truth"].values.astype(float)
            tools = tool_cols(ptx)

            for tool in sorted(tools):
                obs = ptx[tool].values.astype(float)
                r = compute_fn_fp(truth, obs)
                all_fn.setdefault(tool, []).append(r["fn"])
                all_fp.setdefault(tool, []).append(r["fp"])
                print(
                    f"  {cond:<45} {tool:<22} "
                    f"{r['fp']:>6}/{r['fp_total']:<6} {100*r['fp_rate']:>6.2f}% "
                    f"{r['fn']:>6}/{r['fn_total']:<6} {100*r['fn_rate']:>6.2f}%"
                )

        # Averages
        subsection(f"Average FN/FP across conditions ({aligner})")
        for tool in sorted(all_fn.keys()):
            avg_fn = np.mean(all_fn[tool])
            avg_fp = np.mean(all_fp[tool])
            print(f"  {tool:<22}: avg FP={avg_fp:.0f}  avg FN(>10)={avg_fn:.0f}")


# ═══════════════════════════════════════════════════════════
# 4. PRUNING IMPACT ON FALSE NEGATIVES (detailed)
# ═══════════════════════════════════════════════════════════
def analyze_pruning_fn_detail():
    section("5. PRUNING FN DETAIL: Which transcripts are rescued by disabling pruning?")

    # Use representative condition
    cond = "gdna_low_ss_0.90_nrna_default"
    for aligner in ["oracle", "minimap2"]:
        ptx = load_per_transcript(cond, aligner)
        if ptx is None:
            continue

        tools = tool_cols(ptx)
        # Find default and no_prune tools for this aligner
        default_tool = [t for t in tools if "default" in t]
        noprune_tool = [t for t in tools if "no_prune" in t]
        if not default_tool or not noprune_tool:
            continue
        default_tool = default_tool[0]
        noprune_tool = noprune_tool[0]

        subsection(f"Condition: {cond} / {aligner}")

        truth = ptx["mrna_truth"].values.astype(float)
        obs_def = ptx[default_tool].values.astype(float)
        obs_nop = ptx[noprune_tool].values.astype(float)

        # FN in default but NOT in no_prune (rescued by disabling pruning)
        high_truth = truth > 10
        fn_def = high_truth & (obs_def < 0.5)
        fn_nop = high_truth & (obs_nop < 0.5)
        rescued = fn_def & ~fn_nop
        still_fn = fn_def & fn_nop
        new_fn = ~fn_def & fn_nop  # became FN when pruning disabled (shouldn't happen)

        n_fn_def = fn_def.sum()
        n_fn_nop = fn_nop.sum()
        n_rescued = rescued.sum()
        n_still = still_fn.sum()
        n_new = new_fn.sum()

        print(f"  FN with pruning ON:  {n_fn_def}")
        print(f"  FN with pruning OFF: {n_fn_nop}")
        print(f"  Rescued (FN→detected): {n_rescued}")
        print(f"  Still FN in both: {n_still}")
        print(f"  New FN (only in no_prune): {n_new}")

        if n_rescued > 0:
            subsection(f"Top rescued transcripts (pruning OFF recovers them) [{aligner}]")
            ptx_rescued = ptx[rescued].copy()
            ptx_rescued["truth"] = truth[rescued]
            ptx_rescued["pred_default"] = obs_def[rescued]
            ptx_rescued["pred_noprune"] = obs_nop[rescued]
            ptx_rescued = ptx_rescued.sort_values("truth", ascending=False)
            print(f"  {'Transcript':<25} {'Gene':<15} {'Truth':>8} {'Default':>10} {'NoPrune':>10}")
            for _, row in ptx_rescued.head(30).iterrows():
                print(f"  {row['transcript_id']:<25} {row['gene_name']:<15} "
                      f"{row['truth']:>8.0f} {row['pred_default']:>10.1f} {row['pred_noprune']:>10.1f}")

        # Characterize rescued transcripts: isoform count
        if n_rescued > 0:
            gene_counts = ptx.groupby("gene_id").size()
            rescued_genes = ptx.loc[rescued, "gene_id"]
            rescued_iso_counts = rescued_genes.map(gene_counts)
            print(f"\n  Rescued transcript isoform distribution:")
            print(f"    Median isoforms per gene: {rescued_iso_counts.median():.0f}")
            print(f"    Mean isoforms per gene: {rescued_iso_counts.mean():.1f}")
            print(f"    Single-isoform: {(rescued_iso_counts == 1).sum()}")
            print(f"    Multi-isoform: {(rescued_iso_counts > 1).sum()}")


# ═══════════════════════════════════════════════════════════
# 5. REALISTIC COMPARISON: minimap2 rigel vs salmon vs kallisto
# ═══════════════════════════════════════════════════════════
def analyze_realistic_comparison(df: pd.DataFrame):
    section("6. REALISTIC COMPARISON: rigel_minimap2 vs salmon vs kallisto")

    tx = df[df["level"] == "transcript"].copy()
    present = [t for t in REALISTIC_TOOLS if t in tx["tool"].unique()]
    tx = tx[tx["tool"].isin(present)]

    subsection("6a. Overall averages (minimap2 rigel variants + salmon + kallisto)")
    avg = (
        tx.groupby("tool")[["mean_abs_error", "rmse", "pearson", "spearman"]]
        .mean().round(4)
    )
    print(avg.to_string())

    subsection("6b. Per-condition MAE comparison")
    conditions = sorted(tx["dataset"].unique())
    header = f"{'Condition':<50}"
    for t in present:
        header += f" {SHORT.get(t, t):>14}"
    print(header)
    print("-" * (50 + 15 * len(present)))

    for cond in conditions:
        row_str = f"{cond:<50}"
        for t in present:
            val = tx[(tx["dataset"] == cond) & (tx["tool"] == t)]["mean_abs_error"]
            row_str += f" {val.values[0]:>14.4f}" if len(val) > 0 else f" {'—':>14}"
        print(row_str)

    # Averages
    row_str = f"{'AVERAGE':<50}"
    for t in present:
        val = tx[tx["tool"] == t]["mean_abs_error"].mean()
        row_str += f" {val:>14.4f}"
    print(row_str)

    subsection("6c. MAE by gDNA level (realistic tools)")
    for gdna in ["none", "low", "high"]:
        subset = tx[tx["gdna_label"] == gdna]
        row_str = f"  gDNA={gdna:<6}"
        for t in present:
            val = subset[subset["tool"] == t]["mean_abs_error"].mean()
            row_str += f" {SHORT.get(t, t):>14}={val:.4f}"
        print(row_str)

    subsection("6d. MAE by strand specificity (realistic tools)")
    for ss in [0.5, 0.9, 1.0]:
        subset = tx[tx["strand_specificity"] == ss]
        row_str = f"  SS={ss:<6}"
        for t in present:
            val = subset[subset["tool"] == t]["mean_abs_error"].mean()
            row_str += f" {SHORT.get(t, t):>14}={val:.4f}"
        print(row_str)

    # Gene-level
    subsection("6e. Gene-level realistic comparison")
    gene = df[(df["level"] == "gene")].copy()
    gene = gene[gene["tool"].isin(present)]
    avg_gene = gene.groupby("tool")[["mean_abs_error", "rmse", "pearson"]].mean().round(4)
    print(avg_gene.to_string())


# ═══════════════════════════════════════════════════════════
# 6. EXPRESSION-STRATIFIED ACCURACY
# ═══════════════════════════════════════════════════════════
def analyze_expression_strata():
    section("7. EXPRESSION-STRATIFIED ACCURACY")

    representative = [
        ("gdna_none_ss_0.90_nrna_default", "oracle"),
        ("gdna_low_ss_0.90_nrna_default", "oracle"),
        ("gdna_low_ss_0.90_nrna_default", "minimap2"),
    ]

    for cond, aligner in representative:
        ptx = load_per_transcript(cond, aligner)
        if ptx is None:
            print(f"  [SKIP] {cond}/{aligner}")
            continue

        subsection(f"{cond} / {aligner}")
        truth = ptx["mrna_truth"].values.astype(float)
        tools = tool_cols(ptx)
        expressed = truth > 0
        n_expr = expressed.sum()

        truth_expr = truth[expressed]
        q25 = np.percentile(truth_expr, 25)
        q75 = np.percentile(truth_expr, 75)

        strata = {
            f"low (0, {q25:.0f}]": (truth > 0) & (truth <= q25),
            f"mid ({q25:.0f}, {q75:.0f}]": (truth > q25) & (truth <= q75),
            f"high (>{q75:.0f})": truth > q75,
            "zero (truth=0)": truth == 0,
        }
        print(f"Expressed: {n_expr}, Q25={q25:.0f}, Q75={q75:.0f}")

        for tool in sorted(tools):
            obs = ptx[tool].values.astype(float)
            print(f"\n  Tool: {tool}")
            print(f"  {'Stratum':<25} {'N':>8} {'MAE':>10} {'MedAE':>10} {'RMSE':>10}")
            for name, mask in strata.items():
                n = mask.sum()
                if n == 0:
                    continue
                t = truth[mask]
                o = obs[mask]
                ae = np.abs(t - o)
                print(f"  {name:<25} {n:>8} {ae.mean():>10.2f} {np.median(ae):>10.2f} "
                      f"{np.sqrt(np.mean((t - o)**2)):>10.2f}")


# ═══════════════════════════════════════════════════════════
# 7. ERROR DISTRIBUTION
# ═══════════════════════════════════════════════════════════
def analyze_error_distribution():
    section("8. ERROR DISTRIBUTION ANALYSIS")

    for aligner in ["oracle", "minimap2"]:
        cond = "gdna_low_ss_0.90_nrna_default"
        ptx = load_per_transcript(cond, aligner)
        if ptx is None:
            continue

        subsection(f"{cond} / {aligner}")
        truth = ptx["mrna_truth"].values.astype(float)
        tools = tool_cols(ptx)
        expressed = truth > 0

        for tool in sorted(tools):
            obs = ptx[tool].values.astype(float)
            errors = obs[expressed] - truth[expressed]
            abs_errors = np.abs(errors)
            r = compute_fn_fp(truth, obs)

            print(f"Tool: {tool}")
            print(f"  Abs error P50={np.percentile(abs_errors, 50):.2f}, "
                  f"P90={np.percentile(abs_errors, 90):.2f}, "
                  f"P95={np.percentile(abs_errors, 95):.2f}, "
                  f"P99={np.percentile(abs_errors, 99):.2f}")
            print(f"  Signed error: mean={errors.mean():.2f}, median={np.median(errors):.2f}")
            print(f"  FP: {r['fp']}/{r['fp_total']} ({100*r['fp_rate']:.2f}%), "
                  f"FN(>10): {r['fn']}/{r['fn_total']} ({100*r['fn_rate']:.2f}%)")
            print()


# ═══════════════════════════════════════════════════════════
# 8. nRNA IMPACT
# ═══════════════════════════════════════════════════════════
def analyze_nrna_impact():
    section("9. nRNA IMPACT ON ACCURACY")

    cond = "gdna_none_ss_0.90_nrna_default"
    for aligner in ["oracle", "minimap2"]:
        ptx = load_per_transcript(cond, aligner)
        if ptx is None:
            continue

        subsection(f"{cond} / {aligner}")
        truth = ptx["mrna_truth"].values.astype(float)
        nrna_truth = ptx["nrna_truth"].values.astype(float)
        expr_mask = truth > 0
        tools = tool_cols(ptx)

        nrna_ratio = np.zeros_like(truth)
        nrna_ratio[expr_mask] = nrna_truth[expr_mask] / truth[expr_mask]

        ratio_bins = [
            ("no nRNA (=0)", nrna_ratio == 0),
            ("low (0,1]", (nrna_ratio > 0) & (nrna_ratio <= 1)),
            ("mid (1,5]", (nrna_ratio > 1) & (nrna_ratio <= 5)),
            ("high (>5)", nrna_ratio > 5),
        ]

        header = f"{'Stratum':<20} {'N':>6}"
        for t in sorted(tools):
            header += f" {SHORT.get(t, t):>16}"
        print(header)

        for name, mask in ratio_bins:
            combined = expr_mask & mask
            n = combined.sum()
            if n == 0:
                continue
            t = truth[combined]
            row_str = f"{name:<20} {n:>6}"
            for tool in sorted(tools):
                obs = ptx[tool].values.astype(float)[combined]
                mae = np.abs(obs - t).mean()
                row_str += f" {mae:>16.2f}"
            print(row_str)


# ═══════════════════════════════════════════════════════════
# 9. ISOFORM RESOLUTION
# ═══════════════════════════════════════════════════════════
def analyze_isoforms():
    section("10. ISOFORM RESOLUTION: Multi-isoform gene accuracy")

    cond = "gdna_none_ss_0.90_nrna_default"
    for aligner in ["oracle", "minimap2"]:
        ptx = load_per_transcript(cond, aligner)
        if ptx is None:
            continue

        subsection(f"{cond} / {aligner}")
        truth = ptx["mrna_truth"].values.astype(float)
        tools = tool_cols(ptx)
        gene_counts = ptx.groupby("gene_id").size()
        ptx_copy = ptx.copy()
        ptx_copy["n_isoforms"] = ptx_copy["gene_id"].map(gene_counts)
        expr_mask = truth > 0

        iso_bins = [
            ("1 isoform", ptx_copy["n_isoforms"] == 1),
            ("2 isoforms", ptx_copy["n_isoforms"] == 2),
            ("3-5 isoforms", (ptx_copy["n_isoforms"] >= 3) & (ptx_copy["n_isoforms"] <= 5)),
            ("6-10 isoforms", (ptx_copy["n_isoforms"] >= 6) & (ptx_copy["n_isoforms"] <= 10)),
            (">10 isoforms", ptx_copy["n_isoforms"] > 10),
        ]

        header = f"{'Iso stratum':<20} {'N':>8}"
        for t in sorted(tools):
            header += f" {SHORT.get(t, t):>16}"
        print(header)

        for name, mask in iso_bins:
            combined = expr_mask & mask.values
            n = combined.sum()
            if n == 0:
                continue
            t = truth[combined]
            row_str = f"{name:<20} {n:>8}"
            for tool in sorted(tools):
                obs = ptx[tool].values.astype(float)[combined]
                mae = np.abs(obs - t).mean()
                row_str += f" {mae:>16.2f}"
            print(row_str)


# ═══════════════════════════════════════════════════════════
# 10. HEAD-TO-HEAD
# ═══════════════════════════════════════════════════════════
def analyze_head_to_head():
    section("11. HEAD-TO-HEAD PER-TRANSCRIPT COMPARISONS")

    cond = "gdna_low_ss_0.90_nrna_default"

    # Oracle: pruning ON vs OFF
    ptx = load_per_transcript(cond, "oracle")
    if ptx is not None:
        subsection("Oracle: pruning ON vs OFF win rate")
        truth = ptx["mrna_truth"].values.astype(float)
        tools = tool_cols(ptx)
        expr = truth > 0
        def_tool = [t for t in tools if "default" in t]
        nop_tool = [t for t in tools if "no_prune" in t]
        if def_tool and nop_tool:
            obs_def = ptx[def_tool[0]].values.astype(float)
            obs_nop = ptx[nop_tool[0]].values.astype(float)
            ae_def = np.abs(obs_def[expr] - truth[expr])
            ae_nop = np.abs(obs_nop[expr] - truth[expr])
            def_wins = (ae_def < ae_nop).sum()
            nop_wins = (ae_nop < ae_def).sum()
            ties = (ae_def == ae_nop).sum()
            total = expr.sum()
            print(f"  {def_tool[0]} wins: {def_wins}/{total} ({100*def_wins/total:.1f}%)")
            print(f"  {nop_tool[0]} wins: {nop_wins}/{total} ({100*nop_wins/total:.1f}%)")
            print(f"  Ties: {ties}/{total} ({100*ties/total:.1f}%)")

    # Realistic: minimap2 rigel vs salmon vs kallisto
    ptx_orc = load_per_transcript(cond, "oracle")
    ptx_mm2 = load_per_transcript(cond, "minimap2")

    if ptx_orc is not None:
        subsection("Oracle: rigel_default vs salmon vs kallisto win rate")
        truth = ptx_orc["mrna_truth"].values.astype(float)
        expr = truth > 0
        tools_orc = tool_cols(ptx_orc)
        rigel_def = [t for t in tools_orc if "default" in t]
        for base_tool in ["salmon", "kallisto"]:
            if base_tool not in tools_orc or not rigel_def:
                continue
            obs_r = ptx_orc[rigel_def[0]].values.astype(float)
            obs_o = ptx_orc[base_tool].values.astype(float)
            ae_r = np.abs(obs_r[expr] - truth[expr])
            ae_o = np.abs(obs_o[expr] - truth[expr])
            r_wins = (ae_r < ae_o).sum()
            o_wins = (ae_o < ae_r).sum()
            total = expr.sum()
            print(f"  {rigel_def[0]} vs {base_tool}: "
                  f"rigel wins {r_wins}/{total} ({100*r_wins/total:.1f}%), "
                  f"{base_tool} wins {o_wins}/{total} ({100*o_wins/total:.1f}%)")

    if ptx_mm2 is not None and ptx_orc is not None:
        subsection("Realistic: minimap2 rigel_default vs salmon vs kallisto win rate (using oracle truth)")
        truth = ptx_orc["mrna_truth"].values.astype(float)
        expr = truth > 0
        # minimap2 ptx has different number of rows potentially? No, same transcripts
        tools_mm2 = tool_cols(ptx_mm2)
        rigel_mm2_def = [t for t in tools_mm2 if "default" in t]

        # For salmon/kallisto, they are in the oracle ptx
        for base_tool in ["salmon", "kallisto"]:
            if base_tool not in ptx_orc.columns or not rigel_mm2_def:
                continue
            # Need to align by transcript_id
            merged = ptx_orc[["transcript_id", "mrna_truth", base_tool]].merge(
                ptx_mm2[["transcript_id", rigel_mm2_def[0]]],
                on="transcript_id", how="inner",
            )
            t = merged["mrna_truth"].values.astype(float)
            e = t > 0
            obs_r = merged[rigel_mm2_def[0]].values.astype(float)
            obs_o = merged[base_tool].values.astype(float)
            ae_r = np.abs(obs_r[e] - t[e])
            ae_o = np.abs(obs_o[e] - t[e])
            r_wins = (ae_r < ae_o).sum()
            o_wins = (ae_o < ae_r).sum()
            total = e.sum()
            print(f"  {rigel_mm2_def[0]} vs {base_tool}: "
                  f"rigel wins {r_wins}/{total} ({100*r_wins/total:.1f}%), "
                  f"{base_tool} wins {o_wins}/{total} ({100*o_wins/total:.1f}%)")


# ═══════════════════════════════════════════════════════════
# 11. POOL-LEVEL mRNA ACCURACY
# ═══════════════════════════════════════════════════════════
def analyze_pool_level(df: pd.DataFrame):
    section("12. POOL-LEVEL mRNA TOTALS (rigel only)")

    tx = df[df["level"] == "transcript"].copy()
    rigel = tx[tx["tool"].isin(RIGEL_TOOLS)]

    print(f"{'Condition':<45} {'Tool':<22} {'mRNA_truth':>12} {'mRNA_pred':>12} {'Err%':>8}")
    print("-" * 105)

    for _, row in rigel.sort_values(["dataset", "tool"]).iterrows():
        if row["total_truth"] == 0:
            continue
        err_pct = (row["total_observed"] - row["total_truth"]) / row["total_truth"] * 100
        print(f"{row['dataset']:<45} {row['tool']:<22} "
              f"{row['total_truth']:>12.0f} {row['total_observed']:>12.0f} {err_pct:>+7.2f}%")


# ═══════════════════════════════════════════════════════════
# 12. CONDITION SENSITIVITY HEATMAP
# ═══════════════════════════════════════════════════════════
def analyze_condition_heatmap(df: pd.DataFrame):
    section("13. CONDITION SENSITIVITY: MAE heatmap by (gDNA, SS)")

    tx = df[df["level"] == "transcript"].copy()
    present = [t for t in ALL_TOOLS if t in tx["tool"].unique()]

    for tool in present:
        subset = tx[tx["tool"] == tool]
        pivot = subset.pivot_table(
            values="mean_abs_error", index="gdna_label", columns="strand_specificity",
        )
        pivot = pivot.reindex(["none", "low", "high"])
        print(f"\n{tool}:")
        print(pivot.round(4).to_string())


# ═══════════════════════════════════════════════════════════
# 13. PERFORMANCE
# ═══════════════════════════════════════════════════════════
def analyze_performance(df: pd.DataFrame):
    section("14. PERFORMANCE (speed, memory)")

    tx = df[df["level"] == "transcript"].copy()
    present = [t for t in ALL_TOOLS if t in tx["tool"].unique()]

    print(f"{'Tool':<25} {'Avg Time(s)':>12} {'Avg RSS(MB)':>12} {'Avg Throughput':>16}")
    print("-" * 70)
    for tool in present:
        subset = tx[tx["tool"] == tool]
        avg_t = subset["elapsed_sec"].mean()
        avg_rss = subset["peak_rss_mb"].mean()
        avg_tp = subset["throughput_frags_per_sec"].mean()
        print(f"{tool:<25} {avg_t:>12.1f} {avg_rss:>12.0f} {avg_tp:>16.0f}")


# ═══════════════════════════════════════════════════════════
# 14. COMPETITIVE RANKING
# ═══════════════════════════════════════════════════════════
def analyze_ranking(df: pd.DataFrame):
    section("15. COMPETITIVE RANKING")

    tx = df[df["level"] == "transcript"].copy()
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

    wins = ranks_df[ranks_df["rank"] == 1].groupby("tool").size()
    print("Conditions where each tool ranks #1:")
    print(wins.to_string())

    # Realistic ranking (exclude oracle tools)
    subsection("Realistic ranking (minimap2 rigel + salmon + kallisto only)")
    realistic = tx[tx["tool"].isin(REALISTIC_TOOLS)]
    tool_ranks_r = []
    for cond in conditions:
        cond_data = realistic[realistic["dataset"] == cond]
        if cond_data.empty:
            continue
        ranked = cond_data.sort_values("mean_abs_error")
        for rank, (_, row) in enumerate(ranked.iterrows(), 1):
            tool_ranks_r.append({"tool": row["tool"], "rank": rank})
    if tool_ranks_r:
        ranks_r = pd.DataFrame(tool_ranks_r)
        avg_r = ranks_r.groupby("tool")["rank"].mean().sort_values()
        print(avg_r.round(2).to_string())


# ═══════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════
def main():
    print("=" * 90)
    print("BENCHMARK V6 DEEP-DIVE ANALYSIS: Pruning ON/OFF × Oracle/Minimap2")
    print(f"Data: {OUTDIR}")
    print("=" * 90)

    df = load_summary()
    print(f"\nSummary: {len(df)} rows, tools: {sorted(df['tool'].unique())}")

    analyze_summary(df)
    analyze_pruning_impact(df)
    analyze_fn_fp_pruning()
    analyze_pruning_fn_detail()
    analyze_realistic_comparison(df)
    analyze_expression_strata()
    analyze_error_distribution()
    analyze_nrna_impact()
    analyze_isoforms()
    analyze_head_to_head()
    analyze_pool_level(df)
    analyze_condition_heatmap(df)
    analyze_performance(df)
    analyze_ranking(df)

    print("\n" + "=" * 90)
    print("ANALYSIS COMPLETE")
    print("=" * 90)


if __name__ == "__main__":
    main()
