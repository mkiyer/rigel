#!/usr/bin/env python3
"""Deep analysis of VCaP benchmark after VBEM_CLAMP_FLOOR fix.

Compares VBEM vs MAP-EM accuracy, analyzes locus-level convergence from
locus_stats.feather, checks zero-forcing regression, and produces a
detailed diagnostic report.

Usage:
    conda activate rigel
    python scripts/debug/vcap_clamp_floor_analysis.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

# Paths
BENCHMARK_DIR = Path(
    "/scratch/mkiyer_root/mkiyer0/shared_data/rigel_benchmarks/ccle_vcap_prostate"
)
CONDITION = "gdna_high_ss_0.90_nrna_none"
VBEM_DIR = BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "vbem"
MAP_DIR = BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "map"
TRUTH_PATH = BENCHMARK_DIR / "truth_abundances_nrna_none.tsv"
ANALYSIS_DIR = Path("results/vcap")


def load_quant(tool_dir: Path) -> pd.DataFrame:
    """Load transcript-level quantification."""
    return pd.read_feather(tool_dir / "quant.feather")


def load_gene_quant(tool_dir: Path) -> pd.DataFrame:
    """Load gene-level quantification."""
    return pd.read_feather(tool_dir / "gene_quant.feather")


def load_truth() -> pd.DataFrame:
    """Load ground truth abundances."""
    return pd.read_csv(TRUTH_PATH, sep="\t")


def load_locus_stats(tool_dir: Path) -> pd.DataFrame | None:
    """Load locus stats if available."""
    locus_stats_path = tool_dir / "locus_stats.feather"
    if locus_stats_path.exists():
        return pd.read_feather(locus_stats_path)
    return None


def load_summary(tool_dir: Path) -> dict:
    """Load summary.json."""
    with open(tool_dir / "summary.json") as f:
        return json.load(f)


def compute_metrics(truth_col: np.ndarray, pred_col: np.ndarray) -> dict:
    """Compute correlation and error metrics."""
    mask = truth_col > 0
    t = truth_col[mask]
    p = pred_col[mask]

    # Log2 transform with pseudocount
    log2_t = np.log2(t + 1)
    log2_p = np.log2(p + 1)

    pearson_r = np.corrcoef(log2_t, log2_p)[0, 1]
    spearman_r = pd.Series(log2_t).corr(pd.Series(log2_p), method="spearman")

    # MAE and RMSE in log2 space
    residuals = log2_p - log2_t
    mae = np.mean(np.abs(residuals))
    rmse = np.sqrt(np.mean(residuals**2))

    # WARE (weighted absolute relative error)
    total_truth = np.sum(truth_col)
    if total_truth > 0:
        weights = truth_col / total_truth
        ware = np.sum(weights * np.abs(pred_col - truth_col) / np.maximum(truth_col, 1))
    else:
        ware = np.nan

    # False zeros: truth > 1 TPM but prediction == 0
    tpm_threshold = 1.0
    false_zeros = np.sum((truth_col > tpm_threshold) & (pred_col == 0))

    return {
        "pearson_r": pearson_r,
        "spearman_r": spearman_r,
        "mae": mae,
        "rmse": rmse,
        "ware": ware,
        "false_zeros": int(false_zeros),
    }


def analyze_convergence(locus_stats: pd.DataFrame) -> dict:
    """Analyze EM convergence from locus stats."""
    n = len(locus_stats)
    max_sq_iters = 333  # 1000 / SQUAREM_BUDGET_DIVISOR

    locus_stats = locus_stats.copy()
    locus_stats["hit_max"] = locus_stats["squarem_iterations"] >= max_sq_iters
    locus_stats["converged"] = ~locus_stats["hit_max"]
    locus_stats["squarem_time_s"] = locus_stats["squarem_us"] / 1e6
    locus_stats["total_time_s"] = locus_stats["total_us"] / 1e6

    n_converged = int(locus_stats["converged"].sum())
    n_hit_max = int(locus_stats["hit_max"].sum())

    # Mega-locus stats
    mega = locus_stats[locus_stats["is_mega_locus"]]

    # Iteration distribution
    iter_percentiles = {
        "p50": float(locus_stats["squarem_iterations"].quantile(0.5)),
        "p90": float(locus_stats["squarem_iterations"].quantile(0.9)),
        "p95": float(locus_stats["squarem_iterations"].quantile(0.95)),
        "p99": float(locus_stats["squarem_iterations"].quantile(0.99)),
        "max": int(locus_stats["squarem_iterations"].max()),
    }

    # Time distribution
    time_percentiles = {
        "total_s": float(locus_stats["squarem_time_s"].sum()),
        "mean_s": float(locus_stats["squarem_time_s"].mean()),
        "p90_s": float(locus_stats["squarem_time_s"].quantile(0.9)),
        "p99_s": float(locus_stats["squarem_time_s"].quantile(0.99)),
        "max_s": float(locus_stats["squarem_time_s"].max()),
    }

    # Top 10 most expensive loci
    top_expensive = locus_stats.nlargest(10, "squarem_time_s")[
        ["locus_idx", "n_transcripts", "n_components", "n_equiv_classes",
         "squarem_iterations", "squarem_time_s", "is_mega_locus"]
    ].to_dict("records")

    result = {
        "n_loci": n,
        "n_converged": n_converged,
        "n_hit_max_iters": n_hit_max,
        "convergence_rate": n_converged / n if n > 0 else 0,
        "n_mega_loci": int(mega["is_mega_locus"].sum()),
        "iteration_percentiles": iter_percentiles,
        "time_distribution": time_percentiles,
        "top_10_expensive": top_expensive,
    }

    if len(mega) > 0:
        result["mega_locus_details"] = mega[
            ["locus_idx", "n_transcripts", "n_components", "n_equiv_classes",
             "squarem_iterations", "squarem_time_s", "total_time_s",
             "ec_total_elements"]
        ].to_dict("records")

    return result


def analyze_false_zeros(
    truth: pd.DataFrame,
    vbem_quant: pd.DataFrame,
    map_quant: pd.DataFrame,
) -> dict:
    """Analyze false zeros: transcripts with truth > threshold but estimate == 0."""
    merged_v = truth.merge(
        vbem_quant[["transcript_id", "count"]],
        on="transcript_id", how="left", suffixes=("", "_vbem"),
    )
    merged_v["count"] = merged_v["count"].fillna(0)

    merged_m = truth.merge(
        map_quant[["transcript_id", "count"]],
        on="transcript_id", how="left", suffixes=("", "_map"),
    )
    merged_m["count"] = merged_m["count"].fillna(0)

    thresholds = [0.1, 1.0, 10.0, 100.0]
    result = {}
    for thresh in thresholds:
        expressed = truth["mrna_abundance"] > thresh
        vbem_zeros = (merged_v["count"] == 0) & expressed
        map_zeros = (merged_m["count"] == 0) & expressed

        # Transcripts zero-forced in VBEM but not MAP
        vbem_only = vbem_zeros & ~map_zeros

        result[f"truth_gt_{thresh}"] = {
            "n_expressed": int(expressed.sum()),
            "vbem_false_zeros": int(vbem_zeros.sum()),
            "map_false_zeros": int(map_zeros.sum()),
            "vbem_only_zeros": int(vbem_only.sum()),
        }

        # Top zero-forced in VBEM (by truth abundance)
        if vbem_only.sum() > 0:
            top_zeros = merged_v[vbem_only].nlargest(
                min(20, int(vbem_only.sum())), "mrna_abundance"
            )[["transcript_id", "gene_id", "gene_name", "mrna_abundance"]].to_dict("records")
            result[f"truth_gt_{thresh}"]["top_vbem_false_zeros"] = top_zeros

    return result


def compare_tools(
    truth: pd.DataFrame,
    vbem_quant: pd.DataFrame,
    map_quant: pd.DataFrame,
) -> dict:
    """Side-by-side comparison of VBEM vs MAP-EM."""
    # Merge
    merged = truth[["transcript_id", "mrna_abundance"]].merge(
        vbem_quant[["transcript_id", "count"]].rename(columns={"count": "vbem"}),
        on="transcript_id", how="left",
    ).merge(
        map_quant[["transcript_id", "count"]].rename(columns={"count": "map"}),
        on="transcript_id", how="left",
    )
    merged["vbem"] = merged["vbem"].fillna(0)
    merged["map"] = merged["map"].fillna(0)

    truth_arr = merged["mrna_abundance"].values
    vbem_arr = merged["vbem"].values
    map_arr = merged["map"].values

    vbem_metrics = compute_metrics(truth_arr, vbem_arr)
    map_metrics = compute_metrics(truth_arr, map_arr)

    # Expression-stratified comparison
    bins = [0, 1, 10, 100, 1000, np.inf]
    labels = ["zero_to_1", "1_to_10", "10_to_100", "100_to_1000", "1000+"]
    merged["expr_bin"] = pd.cut(merged["mrna_abundance"], bins=bins, labels=labels, right=False)

    stratified = {}
    for label in labels:
        sub = merged[merged["expr_bin"] == label]
        if len(sub) == 0:
            continue
        t = sub["mrna_abundance"].values
        v = sub["vbem"].values
        m = sub["map"].values
        stratified[label] = {
            "n": len(sub),
            "vbem": compute_metrics(t, v),
            "map": compute_metrics(t, m),
        }

    return {
        "overall": {
            "vbem": vbem_metrics,
            "map": map_metrics,
        },
        "stratified": stratified,
    }


def main():
    print("=" * 70)
    print("VCaP VBEM_CLAMP_FLOOR Fix — Deep Analysis Report")
    print("=" * 70)
    print()

    # Load data
    truth = load_truth()
    vbem_quant = load_quant(VBEM_DIR)
    map_quant = load_quant(MAP_DIR)
    vbem_summary = load_summary(VBEM_DIR)
    map_summary = load_summary(MAP_DIR)
    vbem_locus_stats = load_locus_stats(VBEM_DIR)

    print("== 1. OVERALL ACCURACY COMPARISON ==")
    print()
    comparison = compare_tools(truth, vbem_quant, map_quant)

    print(f"{'Metric':<25} {'VBEM (new)':>15} {'MAP-EM':>15} {'Old VBEM':>15}")
    print("-" * 70)
    old_vbem = {"pearson_r": 0.8145, "ware": 0.4037, "false_zeros": 4865}
    for m in ["pearson_r", "spearman_r", "rmse", "mae", "ware", "false_zeros"]:
        v_val = comparison["overall"]["vbem"][m]
        m_val = comparison["overall"]["map"][m]
        old = old_vbem.get(m, "—")
        if isinstance(v_val, float):
            print(f"{m:<25} {v_val:>15.4f} {m_val:>15.4f} {old:>15}")
        else:
            print(f"{m:<25} {v_val:>15} {m_val:>15} {old:>15}")

    print()
    print("== 2. EXPRESSION-STRATIFIED METRICS ==")
    print()
    print(f"{'Bin':<15} {'N':>6} {'VBEM R':>10} {'MAP R':>10} {'VBEM WARE':>12} {'MAP WARE':>12} {'VBEM FZ':>8} {'MAP FZ':>8}")
    print("-" * 85)
    for label, data in comparison["stratified"].items():
        print(
            f"{label:<15} {data['n']:>6} "
            f"{data['vbem']['pearson_r']:>10.4f} "
            f"{data['map']['pearson_r']:>10.4f} "
            f"{data['vbem']['ware']:>12.4f} "
            f"{data['map']['ware']:>12.4f} "
            f"{data['vbem']['false_zeros']:>8} "
            f"{data['map']['false_zeros']:>8}"
        )

    print()
    print("== 3. FALSE ZEROS ANALYSIS ==")
    print()
    false_zeros = analyze_false_zeros(truth, vbem_quant, map_quant)
    for thresh_key, data in false_zeros.items():
        print(f"  {thresh_key}:")
        print(f"    Expressed: {data['n_expressed']}")
        print(f"    VBEM false zeros: {data['vbem_false_zeros']}")
        print(f"    MAP false zeros:  {data['map_false_zeros']}")
        print(f"    VBEM-only zeros:  {data['vbem_only_zeros']}")
        if "top_vbem_false_zeros" in data and data["top_vbem_false_zeros"]:
            print(f"    Top VBEM-only zeros (by truth):")
            for tz in data["top_vbem_false_zeros"][:10]:
                print(f"      {tz['transcript_id']:25s} {tz.get('gene_name',''):15s} truth={tz['mrna_abundance']:.1f}")
        print()

    print()
    print("== 4. POOL-LEVEL SUMMARY ==")
    print()
    for label, summary in [("VBEM", vbem_summary), ("MAP", map_summary)]:
        quant = summary.get("quantification", {})
        print(f"  {label}:")
        print(f"    Total mRNA:  {quant.get('total_mrna', 'N/A'):>15}")
        print(f"    Total nRNA:  {quant.get('total_nrna', 'N/A'):>15}")
        print(f"    Total gDNA:  {quant.get('total_gdna', 'N/A'):>15}")
        print(f"    N loci:      {quant.get('n_loci', 'N/A'):>15}")
        print()

    if vbem_locus_stats is not None:
        print()
        print("== 5. CONVERGENCE ANALYSIS (VBEM) ==")
        print()
        conv = analyze_convergence(vbem_locus_stats)

        print(f"  Total loci: {conv['n_loci']}")
        print(f"  Converged:  {conv['n_converged']} ({conv['convergence_rate']*100:.1f}%)")
        print(f"  Hit max:    {conv['n_hit_max_iters']}")
        print(f"  Mega-loci:  {conv['n_mega_loci']}")
        print()

        print("  Iteration percentiles:")
        for k, v in conv["iteration_percentiles"].items():
            print(f"    {k}: {v}")
        print()

        print("  Timing:")
        for k, v in conv["time_distribution"].items():
            print(f"    {k}: {v:.3f}")
        print()

        print("  Top 10 most expensive loci:")
        print(f"  {'Idx':>6} {'TX':>8} {'Comp':>8} {'ECs':>8} {'Iters':>6} {'Time(s)':>10} {'Mega':>5}")
        print("  " + "-" * 60)
        for loc in conv["top_10_expensive"]:
            print(
                f"  {loc['locus_idx']:>6} {loc['n_transcripts']:>8} "
                f"{loc['n_components']:>8} {loc['n_equiv_classes']:>8} "
                f"{loc['squarem_iterations']:>6} {loc['squarem_time_s']:>10.2f} "
                f"{'Y' if loc['is_mega_locus'] else 'N':>5}"
            )

        if "mega_locus_details" in conv:
            print()
            print("  Mega-locus details:")
            for ml in conv["mega_locus_details"]:
                print(f"    Locus {ml['locus_idx']}: "
                      f"{ml['n_transcripts']} transcripts, "
                      f"{ml['n_components']} components, "
                      f"{ml['n_equiv_classes']} ECs, "
                      f"{ml['squarem_iterations']} iters, "
                      f"{ml['squarem_time_s']:.1f}s SQUAREM, "
                      f"{ml['total_time_s']:.1f}s total")
    else:
        print()
        print("WARNING: No locus_stats.feather found for VBEM")

    # Delta analysis: VBEM vs MAP per transcript
    print()
    print("== 6. VBEM vs MAP DELTA ANALYSIS ==")
    print()
    merged = truth[["transcript_id", "gene_name", "mrna_abundance"]].merge(
        vbem_quant[["transcript_id", "count"]].rename(columns={"count": "vbem"}),
        on="transcript_id", how="left",
    ).merge(
        map_quant[["transcript_id", "count"]].rename(columns={"count": "map"}),
        on="transcript_id", how="left",
    )
    merged["vbem"] = merged["vbem"].fillna(0)
    merged["map"] = merged["map"].fillna(0)
    merged["vbem_err"] = merged["vbem"] - merged["mrna_abundance"]
    merged["map_err"] = merged["map"] - merged["mrna_abundance"]
    merged["delta"] = np.abs(merged["vbem_err"]) - np.abs(merged["map_err"])

    # Where VBEM is worse than MAP (positive delta = VBEM has larger error)
    vbem_worse = merged[merged["delta"] > 1].sort_values("delta", ascending=False)
    # Where VBEM is better than MAP
    vbem_better = merged[merged["delta"] < -1].sort_values("delta")

    print(f"  Transcripts where |VBEM error| > |MAP error| + 1: {len(vbem_worse)}")
    print(f"  Transcripts where |MAP error| > |VBEM error| + 1: {len(vbem_better)}")
    print()

    if len(vbem_worse) > 0:
        print("  Top 15 where VBEM is worse:")
        print(f"  {'transcript_id':30s} {'gene':15s} {'truth':>10} {'vbem':>10} {'map':>10} {'delta':>10}")
        print("  " + "-" * 85)
        for _, r in vbem_worse.head(15).iterrows():
            print(
                f"  {r['transcript_id']:30s} {str(r.get('gene_name','')):15s} "
                f"{r['mrna_abundance']:>10.1f} {r['vbem']:>10.1f} "
                f"{r['map']:>10.1f} {r['delta']:>10.1f}"
            )
        print()

    if len(vbem_better) > 0:
        print("  Top 15 where VBEM is better:")
        print(f"  {'transcript_id':30s} {'gene':15s} {'truth':>10} {'vbem':>10} {'map':>10} {'delta':>10}")
        print("  " + "-" * 85)
        for _, r in vbem_better.head(15).iterrows():
            print(
                f"  {r['transcript_id']:30s} {str(r.get('gene_name','')):15s} "
                f"{r['mrna_abundance']:>10.1f} {r['vbem']:>10.1f} "
                f"{r['map']:>10.1f} {r['delta']:>10.1f}"
            )

    print()
    print("== 7. SUMMARY ==")
    print()
    v_r = comparison["overall"]["vbem"]["pearson_r"]
    m_r = comparison["overall"]["map"]["pearson_r"]
    v_fz = comparison["overall"]["vbem"]["false_zeros"]
    m_fz = comparison["overall"]["map"]["false_zeros"]
    v_ware = comparison["overall"]["vbem"]["ware"]
    m_ware = comparison["overall"]["map"]["ware"]

    print(f"  VBEM Pearson R:     {v_r:.4f} (was 0.8145, now {v_r:.4f})")
    print(f"  MAP  Pearson R:     {m_r:.4f}")
    print(f"  VBEM WARE:          {v_ware:.4f} (was 0.4037, now {v_ware:.4f})")
    print(f"  MAP  WARE:          {m_ware:.4f}")
    print(f"  VBEM false zeros:   {v_fz} (was 4865, now {v_fz})")
    print(f"  MAP  false zeros:   {m_fz}")
    print()
    if v_r > 0.98:
        print("  RESULT: VBEM accuracy RESTORED. Clamp floor fix is effective.")
    else:
        print("  WARNING: VBEM accuracy not fully restored.")


if __name__ == "__main__":
    main()
