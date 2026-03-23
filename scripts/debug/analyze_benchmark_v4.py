#!/usr/bin/env python3
"""Deep analysis of benchmark v4 vs v3 results."""

import pandas as pd
import numpy as np
from pathlib import Path

V3 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v3/gdna_none_ss_0.95_nrna_none")
V4 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/gdna_none_ss_0.95_nrna_none")


def load_oracle(path):
    df = pd.read_csv(path / "per_transcript_counts_oracle.csv")
    return df


def load_minimap2(path):
    df = pd.read_csv(path / "per_transcript_counts_minimap2.csv")
    return df


def mae(pred, truth):
    return np.mean(np.abs(pred - truth))


def rmse(pred, truth):
    return np.sqrt(np.mean((pred - truth) ** 2))


def spearman_active(pred, truth):
    from scipy.stats import spearmanr
    mask = truth > 0
    if mask.sum() < 3:
        return np.nan
    return spearmanr(pred[mask], truth[mask]).correlation


def report_tool_metrics(df, tool_col, truth_col="mrna_truth", label=""):
    pred = df[tool_col].values
    truth = df[truth_col].values
    active = truth > 0
    n_active = active.sum()
    n_expressed = (pred > 0).sum()

    abs_err = np.abs(pred - truth)
    signed_err = pred - truth

    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    print(f"  Transcripts (total / active / predicted>0): {len(truth)} / {n_active} / {n_expressed}")
    print(f"  MAE:     {mae(pred, truth):.2f}")
    print(f"  RMSE:    {rmse(pred, truth):.2f}")
    print(f"  Pearson: {np.corrcoef(pred, truth)[0,1]:.6f}")
    print(f"  Spearman (active only): {spearman_active(pred, truth):.4f}")

    # Among active transcripts
    print(f"\n  --- Active transcripts only (truth > 0, n={n_active}) ---")
    print(f"  MAE:     {mae(pred[active], truth[active]):.2f}")
    print(f"  RMSE:    {rmse(pred[active], truth[active]):.2f}")
    print(f"  Median abs error: {np.median(abs_err[active]):.2f}")
    print(f"  Total truth: {truth[active].sum():.0f}")
    print(f"  Total pred:  {pred[active].sum():.0f}")
    print(f"  Leakage (pred > 0 where truth == 0): {((pred > 0) & (~active)).sum()}")

    # Relative error for active transcripts with truth >= 10
    big = truth >= 10
    if big.sum() > 0:
        rel_err = np.abs(pred[big] - truth[big]) / truth[big]
        print(f"\n  --- Transcripts with truth >= 10 (n={big.sum()}) ---")
        print(f"  Median relative error: {np.median(rel_err):.4f}")
        print(f"  Mean relative error:   {np.mean(rel_err):.4f}")
        print(f"  95th pctile rel error: {np.percentile(rel_err, 95):.4f}")
        print(f"  Max relative error:    {np.max(rel_err):.4f}")

    # Error distribution by expression level
    print(f"\n  --- Error by expression bucket ---")
    bins = [(0, 0, "zero"), (1, 10, "1-10"), (11, 100, "11-100"),
            (101, 1000, "101-1k"), (1001, 10000, "1k-10k"), (10001, np.inf, "10k+")]
    for lo, hi, name in bins:
        mask = (truth >= lo) & (truth <= hi)
        n = mask.sum()
        if n == 0:
            continue
        me = mae(pred[mask], truth[mask])
        bias = np.mean(signed_err[mask])
        print(f"    [{name:>7s}] n={n:>6d}  MAE={me:>8.2f}  bias={bias:>+8.2f}")

    return abs_err, signed_err


def top_errors(df, tool_col, truth_col="mrna_truth", n=20):
    """Print the top N largest absolute errors."""
    pred = df[tool_col].values
    truth = df[truth_col].values
    abs_err = np.abs(pred - truth)
    idx = np.argsort(abs_err)[::-1][:n]

    print(f"\n  Top {n} errors:")
    print(f"  {'transcript_id':>30s}  {'gene_name':>15s}  {'truth':>8s}  {'pred':>8s}  {'error':>8s}  {'rel_err':>8s}")
    for i in idx:
        t = truth[i]
        p = pred[i]
        e = abs_err[i]
        r = e / t if t > 0 else np.inf
        print(f"  {df.iloc[i]['transcript_id']:>30s}  {str(df.iloc[i]['gene_name']):>15s}  {t:>8.0f}  {p:>8.1f}  {e:>8.1f}  {r:>8.3f}")


def regression_analysis(df_v3, df_v4, tool_col, truth_col="mrna_truth"):
    """Compare v3 vs v4 at the transcript level."""
    # Merge on transcript_id
    m = df_v3[["transcript_id", truth_col, tool_col]].merge(
        df_v4[["transcript_id", tool_col]],
        on="transcript_id", suffixes=("_v3", "_v4"),
    )
    truth = m[truth_col].values
    pred_v3 = m[f"{tool_col}_v3"].values
    pred_v4 = m[f"{tool_col}_v4"].values

    err_v3 = np.abs(pred_v3 - truth)
    err_v4 = np.abs(pred_v4 - truth)
    delta_err = err_v4 - err_v3  # positive = v4 is worse

    print(f"\n{'='*60}")
    print(f"  Regression: v4 vs v3 for {tool_col}")
    print(f"{'='*60}")
    print(f"  Transcripts compared: {len(m)}")
    print(f"  v3 MAE: {np.mean(err_v3):.2f}")
    print(f"  v4 MAE: {np.mean(err_v4):.2f}")
    print(f"  Delta MAE (v4 - v3): {np.mean(delta_err):+.2f}")
    print(f"  v3 RMSE: {np.sqrt(np.mean(err_v3**2)):.2f}")
    print(f"  v4 RMSE: {np.sqrt(np.mean(err_v4**2)):.2f}")

    improved = delta_err < -0.5
    regressed = delta_err > 0.5
    unchanged = ~improved & ~regressed
    print(f"\n  Transcripts improved (|err| decreased by >0.5): {improved.sum()}")
    print(f"  Transcripts regressed (|err| increased by >0.5): {regressed.sum()}")
    print(f"  Transcripts unchanged: {unchanged.sum()}")

    # Top regressions
    if regressed.sum() > 0:
        idx = np.argsort(delta_err)[::-1][:20]
        print(f"\n  Top 20 regressions (v4 worse than v3):")
        print(f"  {'transcript_id':>30s}  {'truth':>8s}  {'v3_pred':>8s}  {'v4_pred':>8s}  {'v3_err':>8s}  {'v4_err':>8s}  {'delta':>8s}")
        for i in idx:
            if delta_err[i] <= 0.5:
                break
            print(f"  {m.iloc[i]['transcript_id']:>30s}  {truth[i]:>8.0f}  {pred_v3[i]:>8.1f}  {pred_v4[i]:>8.1f}  {err_v3[i]:>8.1f}  {err_v4[i]:>8.1f}  {delta_err[i]:>+8.1f}")

    # Top improvements
    if improved.sum() > 0:
        idx = np.argsort(delta_err)[:20]
        print(f"\n  Top 20 improvements (v4 better than v3):")
        print(f"  {'transcript_id':>30s}  {'truth':>8s}  {'v3_pred':>8s}  {'v4_pred':>8s}  {'v3_err':>8s}  {'v4_err':>8s}  {'delta':>8s}")
        for i in idx:
            if delta_err[i] >= -0.5:
                break
            print(f"  {m.iloc[i]['transcript_id']:>30s}  {truth[i]:>8.0f}  {pred_v3[i]:>8.1f}  {pred_v4[i]:>8.1f}  {err_v3[i]:>8.1f}  {err_v4[i]:>8.1f}  {delta_err[i]:>+8.1f}")

    return m, delta_err


def nrna_leakage_analysis(df, label=""):
    """Analyze how much nRNA leakage exists (since ground truth has zero nRNA)."""
    nrna = df["nrna_abundance"].values if "nrna_abundance" in df.columns else None
    if nrna is None:
        return
    total_nrna = nrna.sum()
    total_mrna = df["mrna_abundance"].values.sum() if "mrna_abundance" in df.columns else 0
    print(f"\n  --- nRNA Leakage ({label}) ---")
    print(f"  Total nRNA predicted: {total_nrna:.0f}")
    print(f"  Total mRNA predicted: {total_mrna:.0f}")
    print(f"  nRNA fraction:        {total_nrna / (total_mrna + total_nrna + 1e-10):.6f}")
    if total_nrna > 0:
        nrna_per_tx = nrna[nrna > 0]
        print(f"  Transcripts with nRNA > 0: {len(nrna_per_tx)}")
        print(f"  Median nRNA (where > 0):   {np.median(nrna_per_tx):.2f}")
        print(f"  Max nRNA:                  {np.max(nrna_per_tx):.2f}")
        top_idx = np.argsort(nrna)[::-1][:10]
        print(f"  Top 10 nRNA transcripts:")
        for i in top_idx:
            if nrna[i] <= 0:
                break
            print(f"    {df.iloc[i]['transcript_id']:>30s}  nRNA={nrna[i]:>8.1f}  mRNA={df.iloc[i]['mrna_abundance']:>8.1f}  truth={df.iloc[i]['mrna_truth']:>8.0f}")


def main():
    print("Loading v3 data...")
    v3_oracle = load_oracle(V3)
    v3_minimap2 = load_minimap2(V3)

    print("Loading v4 data...")
    v4_oracle = load_oracle(V4)
    v4_minimap2 = load_minimap2(V4)

    # ===== V4 ORACLE ANALYSIS =====
    report_tool_metrics(v4_oracle, "rigel_oracle", label="v4 Rigel Oracle")
    top_errors(v4_oracle, "rigel_oracle")
    nrna_leakage_analysis(v4_oracle, label="v4 Rigel Oracle")

    report_tool_metrics(v4_oracle, "salmon", label="v4 Salmon")
    report_tool_metrics(v4_oracle, "kallisto", label="v4 Kallisto")

    # ===== V4 MINIMAP2 ANALYSIS =====
    report_tool_metrics(v4_minimap2, "rigel_minimap2", label="v4 Rigel Minimap2")
    top_errors(v4_minimap2, "rigel_minimap2")
    nrna_leakage_analysis(v4_minimap2, label="v4 Rigel Minimap2")

    # ===== REGRESSION: V4 vs V3 =====
    print("\n\n" + "#" * 60)
    print("#  REGRESSION ANALYSIS: V4 vs V3")
    print("#" * 60)

    regression_analysis(v3_oracle, v4_oracle, "rigel_oracle")
    regression_analysis(v3_minimap2, v4_minimap2, "rigel_minimap2")

    # ===== MINIMAP2 vs ORACLE GAP =====
    print("\n\n" + "#" * 60)
    print("#  MINIMAP2 vs ORACLE GAP ANALYSIS")
    print("#" * 60)

    m = v4_oracle[["transcript_id", "mrna_truth", "rigel_oracle"]].merge(
        v4_minimap2[["transcript_id", "rigel_minimap2"]],
        on="transcript_id",
    )
    truth = m["mrna_truth"].values
    ora_err = np.abs(m["rigel_oracle"].values - truth)
    mm2_err = np.abs(m["rigel_minimap2"].values - truth)
    gap = mm2_err - ora_err

    active = truth > 0
    print(f"\n  Oracle MAE (active):  {np.mean(ora_err[active]):.2f}")
    print(f"  Minimap2 MAE (active): {np.mean(mm2_err[active]):.2f}")
    print(f"  Gap (mm2 - ora):       {np.mean(gap[active]):.2f}")

    # Where does the gap come from?
    big_gap = gap > 5
    print(f"\n  Transcripts with gap > 5: {big_gap.sum()}")
    print(f"  Error mass from gap > 5:  {np.sum(gap[big_gap]):.0f} (out of total gap {np.sum(gap[active]):.0f})")

    idx = np.argsort(gap)[::-1][:30]
    print(f"\n  Top 30 minimap2-oracle gaps:")
    print(f"  {'transcript_id':>30s}  {'truth':>8s}  {'oracle':>8s}  {'mm2':>8s}  {'gap':>8s}")
    for i in idx:
        print(f"  {m.iloc[i]['transcript_id']:>30s}  {truth[i]:>8.0f}  {m.iloc[i]['rigel_oracle']:>8.1f}  {m.iloc[i]['rigel_minimap2']:>8.1f}  {gap[i]:>+8.1f}")


if __name__ == "__main__":
    main()
