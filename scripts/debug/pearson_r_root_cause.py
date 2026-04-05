#!/usr/bin/env python3
"""
Root cause analysis: Why Rigel's very-high-expression Pearson R (0.91)
is lower than Kallisto's (0.98) in pristine simulation data.

Approach:
1. Load per_transcript_detail.parquet (all tools)
2. Pivot to get per-transcript side-by-side comparison
3. Identify transcripts where Rigel is uniquely bad
4. If oracle results exist, classify errors: aligner-caused vs inherent
5. Drill into worst-case transcripts: gene structure, locus membership
6. Look for systematic patterns (isoform confusion, mega-loci, etc.)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def load_detail(parquet_path: Path) -> pd.DataFrame:
    """Load per-transcript detail and filter to the single pristine condition."""
    df = pd.read_parquet(parquet_path)
    # Filter to pristine condition only
    conds = df["condition"].unique()
    pristine = [c for c in conds if "gdna_none" in c and "nrna_none" in c]
    if not pristine:
        print(f"Available conditions: {conds}")
        raise ValueError("No pristine condition found")
    cond = pristine[0]
    print(f"Using condition: {cond}")
    return df[df["condition"] == cond].copy()


def pivot_tools(df: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot so each transcript has one row with columns per tool.
    Returns: transcript_id, gene_id, gene_name, truth_tpm,
             pred_<tool>, err_<tool>, rel_err_<tool> for each tool.
    """
    tools = sorted(df["tool"].unique())
    print(f"Tools available: {tools}")

    # Start with truth info (same across tools)
    truth_cols = ["transcript_id", "gene_id", "gene_name", "mrna_abundance",
                  "nrna_abundance", "truth_tpm"]
    # Take first occurrence of truth per transcript
    base = df.drop_duplicates("transcript_id")[truth_cols].copy()

    for tool in tools:
        tdf = df[df["tool"] == tool][["transcript_id", "predicted", "residual",
                                       "abs_error", "rel_error"]].copy()
        safe = tool.replace("/", "_")
        tdf = tdf.rename(columns={
            "predicted": f"pred_{safe}",
            "residual": f"resid_{safe}",
            "abs_error": f"abs_err_{safe}",
            "rel_error": f"rel_err_{safe}",
        })
        base = base.merge(tdf, on="transcript_id", how="outer")

    return base


def expression_bin(tpm: float) -> str:
    if tpm <= 0:
        return "zero"
    elif tpm < 10:
        return "low (0-10)"
    elif tpm < 100:
        return "med (10-100)"
    elif tpm < 1000:
        return "high (100-1K)"
    else:
        return "very high (1K+)"


def compute_pearson_by_bin(pv: pd.DataFrame, tool_col: str, truth_col: str = "truth_tpm"):
    """Compute Pearson R for each expression bin."""
    pv = pv.copy()
    pv["expr_bin"] = pv[truth_col].apply(expression_bin)
    results = []
    for bin_name, grp in pv.groupby("expr_bin"):
        mask = grp[truth_col].notna() & grp[tool_col].notna()
        g = grp[mask]
        if len(g) < 3:
            continue
        r, p = stats.pearsonr(g[truth_col], g[tool_col])
        results.append({
            "bin": bin_name,
            "n": len(g),
            "pearson_r": r,
            "mean_truth": g[truth_col].mean(),
            "mean_pred": g[tool_col].mean(),
        })
    return pd.DataFrame(results)


def find_rigel_unique_errors(pv: pd.DataFrame, rigel_tool: str = "rigel_vbem",
                              kallisto_tool: str = "kallisto",
                              min_truth_tpm: float = 10.0) -> pd.DataFrame:
    """
    Find transcripts where Rigel is much worse than Kallisto.
    Returns sorted by (abs_err_rigel - abs_err_kallisto) descending.
    """
    rcol = f"abs_err_{rigel_tool}"
    kcol = f"abs_err_{kallisto_tool}"

    mask = (pv["truth_tpm"] >= min_truth_tpm) & pv[rcol].notna() & pv[kcol].notna()
    sub = pv[mask].copy()
    sub["rigel_minus_kallisto"] = sub[rcol] - sub[kcol]
    sub["rigel_rel_err"] = sub[f"rel_err_{rigel_tool}"]
    sub["kallisto_rel_err"] = sub[f"rel_err_{kallisto_tool}"]

    return sub.sort_values("rigel_minus_kallisto", ascending=False)


def classify_errors_with_oracle(pv: pd.DataFrame) -> pd.DataFrame:
    """
    If oracle results exist, classify each transcript error into:
    - aligner_caused: Rigel(minimap2) bad, Rigel(oracle) good
    - inherent: Rigel(oracle) also bad
    - kallisto_shared: Both Kallisto and Rigel bad
    - rigel_better: Rigel actually does better than Kallisto
    """
    has_oracle = any("oracle" in c for c in pv.columns)
    if not has_oracle:
        print("No oracle results found — skipping aligner vs inherent classification")
        return pv

    # Find the oracle column
    oracle_col = [c for c in pv.columns if c.startswith("abs_err_") and "oracle" in c]
    if not oracle_col:
        return pv
    oracle_col = oracle_col[0]
    oracle_safe = oracle_col.replace("abs_err_", "")

    rigel_col = "abs_err_rigel_vbem"
    kallisto_col = "abs_err_kallisto"

    mask = pv["truth_tpm"] > 10
    sub = pv[mask].copy()

    def classify(row):
        r_err = row.get(rigel_col, np.nan)
        k_err = row.get(kallisto_col, np.nan)
        o_err = row.get(oracle_col, np.nan)
        truth = row["truth_tpm"]

        if pd.isna(r_err) or pd.isna(k_err):
            return "missing"

        r_rel = r_err / truth if truth > 0 else np.nan
        k_rel = k_err / truth if truth > 0 else np.nan
        o_rel = o_err / truth if truth > 0 and not pd.isna(o_err) else np.nan

        # Threshold: >20% relative error is "bad"
        r_bad = r_rel > 0.2
        k_bad = k_rel > 0.2

        if not r_bad:
            return "rigel_good"
        if r_bad and k_bad:
            return "shared_error"
        if r_bad and not k_bad:
            if not pd.isna(o_rel) and o_rel <= 0.2:
                return "aligner_caused"
            elif not pd.isna(o_rel) and o_rel > 0.2:
                return "inherent_rigel"
            else:
                return "rigel_worse_no_oracle"
        return "other"

    sub["error_class"] = sub.apply(classify, axis=1)
    return sub


def isoform_confusion_analysis(pv: pd.DataFrame, rigel_tool: str = "rigel_vbem"):
    """
    Look for isoform confusion: within a gene, one isoform gets too much signal
    while a sibling gets too little (anti-correlated errors).
    """
    pred_col = f"pred_{rigel_tool}"
    mask = pv["truth_tpm"] > 0
    sub = pv[mask].copy()
    sub["resid"] = sub[pred_col] - sub["truth_tpm"]

    # For each gene, check if there are offsetting errors
    results = []
    for gene_id, grp in sub.groupby("gene_id"):
        if len(grp) < 2:
            continue
        resids = grp["resid"].values
        truth = grp["truth_tpm"].values
        preds = grp[pred_col].values

        # Check for anti-correlation: some positive, some negative residuals
        pos_resid = resids[resids > 0].sum()
        neg_resid = resids[resids < 0].sum()
        total_truth = truth.sum()

        if total_truth < 10:
            continue

        # Gene-level prediction vs truth
        gene_pred = preds.sum()
        gene_truth = truth.sum()
        gene_rel_err = abs(gene_pred - gene_truth) / gene_truth

        # Isoform-level total absolute error
        isoform_total_err = np.abs(resids).sum()

        # Confusion index: how much of the error is internal redistribution vs total error
        # If gene-level is correct but isoform-level is wrong, it's pure confusion
        if isoform_total_err > 0:
            confusion_ratio = 1.0 - (abs(gene_pred - gene_truth) / isoform_total_err)
        else:
            confusion_ratio = 0.0

        results.append({
            "gene_id": gene_id,
            "gene_name": grp["gene_name"].iloc[0],
            "n_isoforms": len(grp),
            "gene_truth_tpm": gene_truth,
            "gene_pred_tpm": gene_pred,
            "gene_rel_err": gene_rel_err,
            "pos_resid_sum": pos_resid,
            "neg_resid_sum": neg_resid,
            "isoform_abs_err_sum": isoform_total_err,
            "confusion_ratio": confusion_ratio,
        })

    return pd.DataFrame(results).sort_values("isoform_abs_err_sum", ascending=False)


def compare_gene_vs_transcript_accuracy(pv: pd.DataFrame, tools: list[str]):
    """
    For each tool, compare gene-level vs transcript-level Pearson R
    in the very-high expression bin.
    """
    mask = pv["truth_tpm"] >= 1000
    sub = pv[mask].copy()

    results = []
    for tool in tools:
        pred_col = f"pred_{tool}"
        if pred_col not in sub.columns:
            continue

        # Transcript-level
        tmask = sub["truth_tpm"].notna() & sub[pred_col].notna()
        tg = sub[tmask]
        if len(tg) < 3:
            continue
        tx_r, _ = stats.pearsonr(tg["truth_tpm"], tg[pred_col])

        # Gene-level: aggregate
        gene = tg.groupby("gene_id").agg(
            gene_truth=("truth_tpm", "sum"),
            gene_pred=(pred_col, "sum"),
        ).reset_index()
        gm = gene[(gene["gene_truth"] > 0) & (gene["gene_pred"].notna())]
        if len(gm) < 3:
            continue
        gene_r, _ = stats.pearsonr(gm["gene_truth"], gm["gene_pred"])

        results.append({
            "tool": tool,
            "n_tx_vh": len(tg),
            "tx_pearson_r": tx_r,
            "n_genes_vh": len(gm),
            "gene_pearson_r": gene_r,
            "r_diff": gene_r - tx_r,
        })

    return pd.DataFrame(results)


def print_section(title: str):
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}\n")


def main():
    parser = argparse.ArgumentParser(description="Root cause analysis of Pearson R gap")
    parser.add_argument("--parquet", type=Path, required=True,
                        help="Path to per_transcript_detail.parquet")
    parser.add_argument("--min-truth-tpm", type=float, default=10.0,
                        help="Minimum truth TPM for error analysis")
    parser.add_argument("--top-n", type=int, default=30,
                        help="Number of top discrepant transcripts to show")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Directory to write output CSVs")
    args = parser.parse_args()

    # Load and pivot data
    print_section("Loading Data")
    df = load_detail(args.parquet)
    pv = pivot_tools(df)
    print(f"Pivoted table: {len(pv)} transcripts")

    # Determine available tools
    tool_cols = [c.replace("pred_", "") for c in pv.columns if c.startswith("pred_")]
    print(f"Tool columns: {tool_cols}")

    # ── Section 1: Pearson R by expression bin ──
    print_section("1. Pearson R by Expression Bin")
    for tool in tool_cols:
        print(f"\n--- {tool} ---")
        bindf = compute_pearson_by_bin(pv, f"pred_{tool}")
        print(bindf.to_string(index=False))

    # ── Section 2: Gene vs Transcript accuracy ──
    print_section("2. Gene-Level vs Transcript-Level Accuracy (VH Expression)")
    gt_df = compare_gene_vs_transcript_accuracy(pv, tool_cols)
    if len(gt_df):
        print(gt_df.to_string(index=False))
    else:
        print("(insufficient data)")

    # ── Section 3: Rigel-unique errors ──
    print_section("3. Transcripts Where Rigel is Uniquely Worse Than Kallisto")
    rigel_tool = "rigel_vbem" if "rigel_vbem" in tool_cols else tool_cols[0]
    kallisto_tool = "kallisto" if "kallisto" in tool_cols else None

    if kallisto_tool:
        errors = find_rigel_unique_errors(pv, rigel_tool, kallisto_tool,
                                           min_truth_tpm=args.min_truth_tpm)
        display_cols = [
            "transcript_id", "gene_name", "truth_tpm",
            f"pred_{rigel_tool}", f"pred_{kallisto_tool}",
            "rigel_rel_err", "kallisto_rel_err", "rigel_minus_kallisto"
        ]
        # Add oracle columns if available
        for col in pv.columns:
            if "oracle" in col and col.startswith("pred_"):
                oracle_safe = col.replace("pred_", "")
                display_cols.insert(5, col)

        print(f"Top {args.top_n} transcripts where Rigel VBEM >> Kallisto (truth >= {args.min_truth_tpm} TPM):")
        print(errors[display_cols].head(args.top_n).to_string(index=False))

        # Summary stats
        print(f"\nTotal transcripts with truth >= {args.min_truth_tpm}: {len(errors)}")
        print(f"Rigel worse: {(errors['rigel_minus_kallisto'] > 0).sum()}")
        print(f"Kallisto worse: {(errors['rigel_minus_kallisto'] < 0).sum()}")
        print(f"Mean diff (rigel - kallisto abs_err): {errors['rigel_minus_kallisto'].mean():.2f}")
        print(f"Median diff: {errors['rigel_minus_kallisto'].median():.2f}")

        # Large discrepancies
        threshold = 50  # TPM of absolute error difference
        big = errors[errors["rigel_minus_kallisto"] > threshold]
        print(f"\nTranscripts with Rigel error > Kallisto error by > {threshold} TPM: {len(big)}")
        if len(big):
            print(big[display_cols].to_string(index=False))
    else:
        print("Kallisto results not available")

    # ── Section 4: Isoform confusion ──
    print_section("4. Isoform Confusion Analysis")
    confusion = isoform_confusion_analysis(pv, rigel_tool)
    # Filter to genes with high confusion (> 50% of error is internal redistribution)
    high_confusion = confusion[
        (confusion["confusion_ratio"] > 0.7) & (confusion["isoform_abs_err_sum"] > 50)
    ].head(args.top_n)
    print(f"Genes with high isoform confusion (>70% internal, >50 TPM error):")
    print(high_confusion.to_string(index=False))

    # Same analysis for Kallisto
    if kallisto_tool:
        print(f"\n--- Kallisto isoform confusion for comparison ---")
        k_confusion = isoform_confusion_analysis(pv, kallisto_tool)
        k_high = k_confusion[
            (k_confusion["confusion_ratio"] > 0.7) & (k_confusion["isoform_abs_err_sum"] > 50)
        ].head(args.top_n)
        print(k_high.to_string(index=False))

    # ── Section 5: Error classification with oracle ──
    print_section("5. Error Classification (Aligner vs Inherent)")
    classified = classify_errors_with_oracle(pv)
    if "error_class" in classified.columns:
        print("\nError class distribution:")
        print(classified["error_class"].value_counts().to_string())

        # For each class, show top examples
        for cls in ["aligner_caused", "inherent_rigel", "shared_error"]:
            subset = classified[classified["error_class"] == cls]
            if len(subset) == 0:
                continue
            print(f"\n--- Top {cls} errors ---")
            display = subset.sort_values(f"abs_err_{rigel_tool}", ascending=False).head(10)
            cols = ["transcript_id", "gene_name", "truth_tpm",
                    f"pred_{rigel_tool}", f"pred_{kallisto_tool}"]
            # Add oracle pred if available
            for c in display.columns:
                if "oracle" in c and c.startswith("pred_"):
                    cols.append(c)
            print(display[cols].to_string(index=False))

    # ── Section 6: Contribution to Pearson R ──
    print_section("6. Contribution to Pearson R Degradation")
    if kallisto_tool and rigel_tool:
        vh = pv[pv["truth_tpm"] >= 1000].copy()
        rp = f"pred_{rigel_tool}"
        kp = f"pred_{kallisto_tool}"

        if len(vh) > 3:
            r_all, _ = stats.pearsonr(vh["truth_tpm"], vh[rp])
            k_all, _ = stats.pearsonr(vh["truth_tpm"], vh[kp])
            print(f"VH Pearson R — Rigel: {r_all:.4f}, Kallisto: {k_all:.4f}, gap: {k_all - r_all:.4f}")

            # Leave-one-out: which transcript's removal most improves Rigel's R?
            improvements = []
            for idx, row in vh.iterrows():
                vh_minus = vh.drop(idx)
                if len(vh_minus) < 3:
                    continue
                r_minus, _ = stats.pearsonr(vh_minus["truth_tpm"], vh_minus[rp])
                improvements.append({
                    "transcript_id": row["transcript_id"],
                    "gene_name": row["gene_name"],
                    "truth_tpm": row["truth_tpm"],
                    f"pred_{rigel_tool}": row[rp],
                    f"pred_{kallisto_tool}": row[kp],
                    "r_without": r_minus,
                    "r_improvement": r_minus - r_all,
                })

            imp_df = pd.DataFrame(improvements).sort_values("r_improvement", ascending=False)
            print(f"\nTop {args.top_n} transcripts whose removal most improves Rigel VH Pearson R:")
            print(imp_df.head(args.top_n).to_string(index=False))

            # Cumulative removal: remove top-N worst offenders
            print("\nCumulative removal effect:")
            top_ids = imp_df.head(20)["transcript_id"].tolist()
            for n in [1, 2, 3, 5, 10, 15, 20]:
                if n > len(top_ids):
                    break
                remove = set(top_ids[:n])
                vh_sub = vh[~vh["transcript_id"].isin(remove)]
                if len(vh_sub) < 3:
                    break
                r_n, _ = stats.pearsonr(vh_sub["truth_tpm"], vh_sub[rp])
                k_n, _ = stats.pearsonr(vh_sub["truth_tpm"], vh_sub[kp])
                print(f"  Remove top {n:2d}: Rigel R={r_n:.4f}, Kallisto R={k_n:.4f}, gap={k_n-r_n:.4f}")

    # ── Section 7: Systematic patterns ──
    print_section("7. Systematic Patterns in Rigel Errors")
    if rigel_tool:
        expressed = pv[pv["truth_tpm"] > 0].copy()
        rp = f"pred_{rigel_tool}"
        expressed["log2_truth"] = np.log2(expressed["truth_tpm"] + 1)
        expressed["log2fc"] = np.log2(expressed[rp] + 1) - np.log2(expressed["truth_tpm"] + 1)

        # Is there a systematic bias (over/under estimation)?
        vh = expressed[expressed["truth_tpm"] >= 1000]
        print(f"Very-high expression (>= 1K TPM): {len(vh)} transcripts")
        print(f"  Mean log2FC (pred/truth): {vh['log2fc'].mean():.4f}")
        print(f"  Median log2FC: {vh['log2fc'].median():.4f}")
        print(f"  Fraction overestimated: {(vh['log2fc'] > 0).mean():.3f}")
        print(f"  Fraction >2x overestimated: {(vh['log2fc'] > 1).mean():.3f}")
        print(f"  Fraction <0.5x underestimated: {(vh['log2fc'] < -1).mean():.3f}")

        # Distribution of errors
        for pct in [50, 75, 90, 95, 99]:
            print(f"  {pct}th percentile |log2FC|: {np.percentile(np.abs(vh['log2fc']), pct):.4f}")

    # ── Section 8: Write outputs ──
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        pv.to_parquet(args.output_dir / "pivoted_comparison.parquet", index=False)
        if kallisto_tool:
            errors.to_parquet(args.output_dir / "rigel_vs_kallisto_errors.parquet", index=False)
        confusion.to_parquet(args.output_dir / "isoform_confusion.parquet", index=False)
        if "error_class" in classified.columns:
            classified.to_parquet(args.output_dir / "error_classification.parquet", index=False)
        print(f"\nOutputs written to {args.output_dir}")


if __name__ == "__main__":
    main()
