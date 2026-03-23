#!/usr/bin/env python3
"""Deep root-cause analysis: oracle vs minimap2 benchmark_pristine.

Usage:
    python3 scripts/debug/analyze_minimap2_errors.py
"""
import pandas as pd
import numpy as np
from pathlib import Path

COND = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/gdna_none_ss_0.95_nrna_none"


def load_data():
    oracle_df = pd.read_csv(f"{COND}/per_transcript_counts_oracle.csv")
    mm2_df    = pd.read_csv(f"{COND}/per_transcript_counts_minimap2.csv")
    df = oracle_df.merge(mm2_df[["transcript_id", "rigel_minimap2"]], on="transcript_id", how="left")
    df["err_oracle"] = np.abs(df["mrna_truth"] - df["rigel_oracle"])
    df["err_mm2"]    = np.abs(df["mrna_truth"] - df["rigel_minimap2"])
    df["err_diff"]   = df["err_mm2"] - df["err_oracle"]
    return df


def overall_summary(df):
    print("=" * 80)
    print("OVERALL TRANSCRIPT-LEVEL METRICS")
    print("=" * 80)
    metrics = {
        "oracle": ("rigel_oracle",),
        "mm2   ": ("rigel_minimap2",),
    }
    for label, (col,) in metrics.items():
        truth = df["mrna_truth"].values.astype(float)
        obs   = df[col].values.astype(float)
        err   = np.abs(truth - obs)
        # Pearson
        if truth.std() > 0 and obs.std() > 0:
            pearson = float(np.corrcoef(truth, obs)[0, 1])
        else:
            pearson = 0.0
        mae  = float(err.mean())
        rmse = float(np.sqrt(np.mean((truth-obs)**2)))
        print(f"  {label}  MAE={mae:8.3f}  RMSE={rmse:8.3f}  Pearson={pearson:.4f}"
              f"  total_err={err.sum():.0f}")
    # Ratio
    mae_or  = df["err_oracle"].mean()
    mae_mm2 = df["err_mm2"].mean()
    print(f"\n  mm2/oracle MAE ratio: {mae_mm2/mae_or:.3f}×")


def error_distribution(df):
    print("\n" + "=" * 80)
    print("ERROR DISTRIBUTION (mm2 additional error per transcript)")
    print("=" * 80)
    diff = df["err_diff"]
    for pct in [50, 75, 90, 95, 99, 99.9]:
        print(f"  P{pct:.0f}: {np.percentile(diff, pct):.1f}")
    n_worse = (diff > 1).sum()
    n_total = len(diff)
    print(f"\n  Transcripts with mm2 err > oracle by >1 count: {n_worse:,} / {n_total:,} ({100*n_worse/n_total:.1f}%)")

    # Distribution by truth abundance bucket
    print("\n  Error by truth abundance bucket:")
    print(f"  {'Bucket':20s}  {'N':>8}  {'oracle_MAE':>10}  {'mm2_MAE':>10}  {'ratio':>8}")
    bins = [(0, 0, "zero"), (1, 10, "1-10"), (11, 100, "11-100"),
            (101, 1000, "101-1000"), (1001, 10000, "1001-10k"),
            (10001, 999999999, ">10k")]
    for lo, hi, label in bins:
        mask = (df["mrna_truth"] >= lo) & (df["mrna_truth"] <= hi)
        sub = df[mask]
        if len(sub) == 0:
            continue
        mae_or = sub["err_oracle"].mean()
        mae_m2 = sub["err_mm2"].mean()
        ratio  = mae_m2 / mae_or if mae_or > 0 else float("nan")
        print(f"  {label:20s}  {len(sub):>8,}  {mae_or:>10.3f}  {mae_m2:>10.3f}  {ratio:>8.3f}")


def gene_error_table(df, top_n=20):
    print("\n" + "=" * 80)
    print(f"TOP {top_n} GENES BY ADDITIONAL ERROR (mm2 vs oracle)")
    print("=" * 80)
    df["gene_id_str"] = df["gene_id"].astype(str)
    gene_err = df.groupby("gene_id_str").agg(
        gene_name    =("gene_name", "first"),
        truth_total  =("mrna_truth", "sum"),
        oracle_total =("rigel_oracle", "sum"),
        mm2_total    =("rigel_minimap2", "sum"),
        err_oracle   =("err_oracle", "sum"),
        err_mm2      =("err_mm2", "sum"),
        n_tx         =("transcript_id", "count"),
    ).reset_index()
    gene_err["err_diff"] = gene_err["err_mm2"] - gene_err["err_oracle"]
    gene_err["mm2_ratio"] = gene_err["mm2_total"] / gene_err["truth_total"].replace(0, np.nan)
    gene_err = gene_err.sort_values("err_diff", ascending=False)

    hdr = f"  {'gene_id':25s}  {'name':18s}  {'truth':>8}  {'oracle':>8}  {'mm2':>8}  {'err_diff':>9}  {'mm2/truth':>9}  ntx"
    print(hdr)
    print("  " + "-" * (len(hdr)-2))
    for _, row in gene_err.head(top_n).iterrows():
        ratio_str = f"{row.mm2_ratio:.3f}" if not np.isnan(row.mm2_ratio) else "   nan"
        print(
            f"  {row.gene_id_str:25s}  {str(row.gene_name):18s}  {row.truth_total:8.0f}  "
            f"{row.oracle_total:8.0f}  {row.mm2_total:8.0f}  {row.err_diff:9.0f}  "
            f"{ratio_str:>9}  {row.n_tx}"
        )

    top_5_err = gene_err.head(5)["err_mm2"].sum()
    total_err  = df["err_mm2"].sum()
    print(f"\n  Top 5 genes: {top_5_err:.0f} / {total_err:.0f} = {100*top_5_err/total_err:.1f}% of total mm2 error")
    return gene_err


def top_gene_deep_dive(df, gene_err, n_genes=5):
    print("\n" + "=" * 80)
    print("DEEP DIVE: TOP 5 ERROR GENES")
    print("=" * 80)
    top_genes = gene_err.head(n_genes)[["gene_id_str", "gene_name"]].values
    for gid, gname in top_genes:
        gene_df = df[df["gene_id"].astype(str) == gid].copy()
        gene_df = gene_df.sort_values("err_diff", ascending=False)
        print(f"\n  {gid}  ({gname}):")
        print(f"  {'transcript_id':42s}  {'truth':>8}  {'oracle':>8}  {'mm2':>8}  {'err_diff':>9}")
        for _, row in gene_df.iterrows():
            print(
                f"    {row.transcript_id:40s}  {row.mrna_truth:8.0f}  "
                f"{row.rigel_oracle:8.0f}  {row.rigel_minimap2:8.0f}  {row.err_diff:9.0f}"
            )


def pseudogene_pattern(df):
    """Identify systematic pseudogene absorption: truth=0, mm2>100."""
    print("\n" + "=" * 80)
    print("PSEUDOGENE / SPURIOUS ABSORPTION (truth=0 but mm2>100)")
    print("=" * 80)
    spurious = df[(df["mrna_truth"] == 0) & (df["rigel_minimap2"] > 100)].copy()
    spurious = spurious.sort_values("rigel_minimap2", ascending=False)
    print(f"  {len(spurious)} transcripts with truth=0 but mm2>100 (total spurious counts: {spurious.rigel_minimap2.sum():.0f})")
    print(f"  {'transcript_id':42s}  {'gene_id':25s}  {'gene_name':18s}  {'mm2':>8}")
    for _, row in spurious.head(20).iterrows():
        print(f"    {row.transcript_id:40s}  {str(row.gene_id):25s}  {str(row.gene_name):18s}  {row.rigel_minimap2:8.0f}")


def lost_fragments_analysis(df):
    """Identify where true fragments are going wrong."""
    print("\n" + "=" * 80)
    print("FRAGMENT LOSS: high-truth transcripts severely under-counted by mm2")
    print("=" * 80)
    lost = df[(df["mrna_truth"] > 500) & (df["rigel_minimap2"] / df["mrna_truth"].replace(0, np.nan) < 0.5)].copy()
    lost["mm2_ratio"] = lost["rigel_minimap2"] / lost["mrna_truth"]
    lost = lost.sort_values("err_diff", ascending=False)
    print(f"  {len(lost)} transcripts where truth>500 and mm2/truth < 0.5")
    print(f"  {'transcript_id':42s}  {'gene_name':18s}  {'truth':>8}  {'oracle':>8}  {'mm2':>8}  {'mm2/truth':>9}  {'err_diff':>9}")
    for _, row in lost.head(20).iterrows():
        print(
            f"    {row.transcript_id:40s}  {str(row.gene_name):18s}  {row.mrna_truth:8.0f}  "
            f"{row.rigel_oracle:8.0f}  {row.rigel_minimap2:8.0f}  {row.mm2_ratio:9.3f}  {row.err_diff:9.0f}"
        )


if __name__ == "__main__":
    df = load_data()
    overall_summary(df)
    error_distribution(df)
    gene_err = gene_error_table(df, top_n=20)
    top_gene_deep_dive(df, gene_err, n_genes=5)
    pseudogene_pattern(df)
    lost_fragments_analysis(df)
