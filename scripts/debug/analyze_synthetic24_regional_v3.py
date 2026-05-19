#!/usr/bin/env python3
"""Compare synthetic-24 regional-v3 outputs against preserved baseline outputs."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy.stats import pearsonr, spearmanr

from rigel.sim.truth import parse_origin


DEFAULT_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24")
POOL_ORDER = ("mrna", "nrna", "gdna", "unresolved")
RUNS = {
    "baseline": ("rigel_out", "annotated.bam"),
    "regional_v3": ("rigel_regional_v3_out", "annotated_regional_v3.bam"),
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def load_frame(path: Path) -> pd.DataFrame:
    if path.with_suffix(".feather").exists():
        return pd.read_feather(path.with_suffix(".feather"))
    if path.with_suffix(".tsv").exists():
        return pd.read_csv(path.with_suffix(".tsv"), sep="\t")
    return pd.DataFrame()


def load_manifest(base: Path) -> tuple[list[dict], dict[str, dict]]:
    manifest = load_json(base / "manifest.json")
    conditions = list(manifest.get("conditions", []))
    return conditions, {str(row["name"]): row for row in conditions}


def true_pool_from_origin(query_name: str) -> str:
    origin = parse_origin(query_name)
    if origin.kind == "gdna":
        return "gdna"
    if origin.kind == "nrna":
        return "nrna"
    return "mrna"


def predicted_pool(read: pysam.AlignedSegment) -> str:
    try:
        zf = int(read.get_tag("ZF"))
    except KeyError:
        zf = 0
    try:
        zc = str(read.get_tag("ZC"))
    except KeyError:
        zc = ""

    if zf & 0x04 or zc == "intergenic":
        return "gdna"
    if zf & 0x08:
        return "nrna"
    if zf & 0x02:
        return "mrna"
    return "unresolved"


def collect_confusion(annotated_bam: Path) -> Counter[tuple[str, str]]:
    counts: Counter[tuple[str, str]] = Counter()
    with pysam.AlignmentFile(str(annotated_bam), "rb") as bam:
        for read in bam:
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue
            counts[(true_pool_from_origin(read.query_name), predicted_pool(read))] += 1
    return counts


def abundance_metrics(base: Path, condition: str, out_dir: Path, condition_meta: dict) -> dict[str, float]:
    quant = load_frame(out_dir / "quant")
    truth_name = str(condition_meta.get("truth_abundances", "truth_abundances.tsv"))
    truth = pd.read_csv(base / truth_name, sep="\t")
    if quant.empty:
        return {"spearman": np.nan, "pearson_log": np.nan, "mard": np.nan, "median_re": np.nan}

    merged = truth.merge(quant[["transcript_id", "tpm", "count", "count_em"]], on="transcript_id", how="left")
    merged = merged.fillna(0.0)
    expressed = merged[merged["mrna_abundance"] > 0].copy()
    if len(expressed) <= 2:
        return {"spearman": np.nan, "pearson_log": np.nan, "mard": np.nan, "median_re": np.nan}

    spearman = float(spearmanr(expressed["mrna_abundance"], expressed["tpm"]).statistic)
    pearson = float(
        pearsonr(
            np.log2(expressed["mrna_abundance"] + 1.0),
            np.log2(expressed["tpm"] + 1.0),
        ).statistic
    )
    rel_err = np.abs(expressed["tpm"] - expressed["mrna_abundance"]) / (
        expressed["mrna_abundance"] + 1.0e-6
    )
    return {
        "spearman": spearman,
        "pearson_log": pearson,
        "mard": float(rel_err.mean()),
        "median_re": float(rel_err.median()),
    }


def run_condition_metrics(base: Path, condition: str, condition_meta: dict, run: str) -> dict[str, object]:
    out_name, annotated_name = RUNS[run]
    cond_dir = base / condition
    out_dir = cond_dir / out_name
    summary = load_json(out_dir / "summary.json")
    loci = load_frame(out_dir / "loci")
    calib = summary.get("calibration", {})
    exposure = calib.get("regional_exposure", {})
    application = calib.get("regional_weighting_stats", {})

    mrna = float(loci["mrna"].sum()) if "mrna" in loci else np.nan
    nrna = float(loci["nrna"].sum()) if "nrna" in loci else 0.0
    gdna = float(loci["gdna"].sum()) if "gdna" in loci else np.nan
    total_locus = mrna + nrna + gdna if np.isfinite(mrna) and np.isfinite(gdna) else np.nan
    gdna_rate = gdna / total_locus if total_locus and total_locus > 0 else np.nan
    nrna_rate = nrna / total_locus if total_locus and total_locus > 0 else np.nan

    ratio = loci["gdna_eff_len_weight_ratio"] if "gdna_eff_len_weight_ratio" in loci else pd.Series(dtype=float)
    signals = [float(v.get("signal", 0.0)) for v in exposure.get("per_class", {}).values()]

    metrics = {
        "condition": condition,
        "run": run,
        "n_mrna_true": int(condition_meta.get("n_mrna", 0)),
        "n_nrna_true": int(condition_meta.get("n_nrna", 0)),
        "n_gdna_true": int(condition_meta.get("n_gdna", 0)),
        "strand_specificity": float(condition_meta.get("strand_specificity", np.nan)),
        "gdna_label": str(condition_meta.get("gdna_label", "")),
        "nrna_label": str(condition_meta.get("nrna_label", "")),
        "status": "ok" if (out_dir / "quant.feather").exists() else "missing",
        "mrna_locus": mrna,
        "nrna_locus": nrna,
        "gdna_locus": gdna,
        "total_locus": total_locus,
        "gdna_locus_rate": gdna_rate,
        "nrna_locus_rate": nrna_rate,
        "exposure_mode": str(exposure.get("mode", "none")),
        "exposure_rho_ref": float(exposure.get("rho_ref", np.nan)) if exposure else np.nan,
        "exposure_max_signal": max(signals) if signals else 0.0,
        "exposure_n_at_floor": int(exposure.get("n_at_floor", 0)) if exposure else 0,
        "weighted_units": int(application.get("n_units_weighted", 0)) if application else 0,
        "seen_units": int(application.get("n_units_seen", 0)) if application else 0,
        "weight_ratio_min": float(ratio.min()) if len(ratio) else np.nan,
        "weight_ratio_p05": float(ratio.quantile(0.05)) if len(ratio) else np.nan,
        "weight_ratio_median": float(ratio.median()) if len(ratio) else np.nan,
        "weight_ratio_p95": float(ratio.quantile(0.95)) if len(ratio) else np.nan,
    }
    metrics.update(abundance_metrics(base, condition, out_dir, condition_meta))
    return metrics


def condition_confusion_rows(base: Path, condition: str, run: str) -> list[dict[str, object]]:
    out_name, annotated_name = RUNS[run]
    annotated_bam = base / condition / annotated_name
    if not annotated_bam.exists():
        return []
    counts = collect_confusion(annotated_bam)
    rows = []
    for true_pool in POOL_ORDER[:-1]:
        for pred_pool in POOL_ORDER:
            rows.append({
                "condition": condition,
                "run": run,
                "true_pool": true_pool,
                "pred_pool": pred_pool,
                "count": counts.get((true_pool, pred_pool), 0),
            })
    return rows


def summarize_confusion(confusion: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (condition, run), grp in confusion.groupby(["condition", "run"]):
        pivot = grp.pivot_table(index="true_pool", columns="pred_pool", values="count", aggfunc="sum", fill_value=0)
        for pool in POOL_ORDER:
            if pool not in pivot.columns:
                pivot[pool] = 0
        total_mrna = float(pivot.loc["mrna"].sum()) if "mrna" in pivot.index else 0.0
        total_nrna = float(pivot.loc["nrna"].sum()) if "nrna" in pivot.index else 0.0
        total_gdna = float(pivot.loc["gdna"].sum()) if "gdna" in pivot.index else 0.0
        gdna_to_rna = float(pivot.loc["gdna", ["mrna", "nrna"]].sum()) if total_gdna else 0.0
        mrna_to_gdna = float(pivot.loc["mrna", "gdna"]) if total_mrna else 0.0
        nrna_to_gdna = float(pivot.loc["nrna", "gdna"]) if total_nrna else 0.0
        rows.append({
            "condition": condition,
            "run": run,
            "total_mrna_fragments": total_mrna,
            "total_nrna_fragments": total_nrna,
            "total_gdna_fragments": total_gdna,
            "mrna_to_gdna": mrna_to_gdna,
            "mrna_to_gdna_rate": mrna_to_gdna / total_mrna if total_mrna else 0.0,
            "nrna_to_gdna": nrna_to_gdna,
            "nrna_to_gdna_rate": nrna_to_gdna / total_nrna if total_nrna else 0.0,
            "gdna_to_rna": gdna_to_rna,
            "gdna_to_rna_rate": gdna_to_rna / total_gdna if total_gdna else 0.0,
            "gdna_correct_rate": float(pivot.loc["gdna", "gdna"]) / total_gdna if total_gdna else np.nan,
        })
    return pd.DataFrame(rows)


def comparison_table(metrics: pd.DataFrame, confusion_summary: pd.DataFrame) -> pd.DataFrame:
    merged = metrics.merge(confusion_summary, on=["condition", "run"], how="left")
    base = merged[merged["run"] == "baseline"].set_index("condition")
    regional = merged[merged["run"] == "regional_v3"].set_index("condition")
    rows = []
    for condition in regional.index:
        if condition not in base.index:
            continue
        left = base.loc[condition]
        right = regional.loc[condition]
        row = {"condition": condition}
        for col in [
            "gdna_locus_rate",
            "nrna_locus_rate",
            "spearman",
            "pearson_log",
            "mard",
            "median_re",
            "gdna_to_rna_rate",
            "mrna_to_gdna_rate",
            "nrna_to_gdna_rate",
        ]:
            row[f"baseline_{col}"] = float(left.get(col, np.nan))
            row[f"regional_{col}"] = float(right.get(col, np.nan))
            row[f"delta_{col}"] = row[f"regional_{col}"] - row[f"baseline_{col}"]
        for col in [
            "exposure_mode",
            "exposure_max_signal",
            "weight_ratio_min",
            "weight_ratio_p05",
            "weight_ratio_median",
            "weighted_units",
            "seen_units",
        ]:
            row[col] = right.get(col, np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def aggregate_confusion_matrix(confusion: pd.DataFrame, run: str) -> pd.DataFrame:
    sub = confusion[confusion["run"] == run]
    matrix = sub.pivot_table(index="true_pool", columns="pred_pool", values="count", aggfunc="sum", fill_value=0)
    for pool in POOL_ORDER:
        if pool not in matrix.columns:
            matrix[pool] = 0
    return matrix[list(POOL_ORDER)].reindex(POOL_ORDER[:-1], fill_value=0)


def markdown_table(frame: pd.DataFrame, *, index: bool = True) -> str:
    """Render a small DataFrame as GitHub-flavored Markdown without tabulate."""
    shown = frame.reset_index() if index else frame.reset_index(drop=True)
    columns = [str(col) for col in shown.columns]

    def cell(value: object) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, float):
            return f"{value:.6g}"
        return str(value)

    lines = ["| " + " | ".join(columns) + " |"]
    lines.append("| " + " | ".join("---" for _ in columns) + " |")
    for _, row in shown.iterrows():
        lines.append("| " + " | ".join(cell(row[col]) for col in shown.columns) + " |")
    return "\n".join(lines)


def write_report(out_dir: Path, metrics: pd.DataFrame, confusion: pd.DataFrame, comparison: pd.DataFrame) -> Path:
    baseline_matrix = aggregate_confusion_matrix(confusion, "baseline")
    regional_matrix = aggregate_confusion_matrix(confusion, "regional_v3")
    regional_metrics = metrics[metrics["run"] == "regional_v3"]
    regional_count = int((regional_metrics["exposure_mode"] == "regional").sum())
    uniform_count = int((regional_metrics["exposure_mode"] == "uniform").sum())

    def fmt(value: float, digits: int = 4) -> str:
        if pd.isna(value):
            return "n/a"
        return f"{value:.{digits}f}"

    max_abs_gdna_rate = float(comparison["delta_gdna_locus_rate"].abs().max())
    max_abs_leak = float(comparison["delta_gdna_to_rna_rate"].abs().max())
    max_abs_spearman = float(comparison["delta_spearman"].abs().max())
    max_abs_mard = float(comparison["delta_mard"].abs().max())
    worst_leak = comparison.reindex(comparison["delta_gdna_to_rna_rate"].abs().sort_values(ascending=False).index).head(8)
    worst_mard = comparison.reindex(comparison["delta_mard"].abs().sort_values(ascending=False).index).head(8)

    lines = [
        "# Synthetic-24 Regional Exposure Assessment",
        "",
        "Compared fresh `rigel_regional_v3_out` outputs against preserved `rigel_out` baseline outputs.",
        "",
        "## Run Summary",
        "",
        f"- Conditions analyzed: {len(regional_metrics)}",
        f"- Regional exposure modes: {regional_count} regional, {uniform_count} uniform",
        f"- Max absolute locus gDNA-rate delta vs baseline: {fmt(max_abs_gdna_rate)}",
        f"- Max absolute gDNA-to-RNA leak-rate delta vs baseline: {fmt(max_abs_leak)}",
        f"- Max absolute Spearman delta vs baseline: {fmt(max_abs_spearman)}",
        f"- Max absolute transcript MARD delta vs baseline: {fmt(max_abs_mard, 3)}",
        "",
        "## Aggregate Pool Confusion - Baseline",
        "",
        markdown_table(baseline_matrix),
        "",
        "## Aggregate Pool Confusion - Regional v3",
        "",
        markdown_table(regional_matrix),
        "",
        "## Largest gDNA-to-RNA Leak Deltas",
        "",
        markdown_table(worst_leak[[
            "condition",
            "baseline_gdna_to_rna_rate",
            "regional_gdna_to_rna_rate",
            "delta_gdna_to_rna_rate",
            "weight_ratio_median",
            "weight_ratio_p05",
            "exposure_max_signal",
        ]], index=False),
        "",
        "## Largest Transcript MARD Deltas",
        "",
        markdown_table(worst_mard[[
            "condition",
            "baseline_mard",
            "regional_mard",
            "delta_mard",
            "baseline_spearman",
            "regional_spearman",
            "delta_spearman",
        ]], index=False),
        "",
        "## Exposure Weight Summary",
        "",
        markdown_table(regional_metrics[[
            "condition",
            "exposure_mode",
            "exposure_max_signal",
            "weight_ratio_min",
            "weight_ratio_p05",
            "weight_ratio_median",
            "weight_ratio_p95",
            "weighted_units",
            "seen_units",
        ]], index=False),
        "",
    ]
    report_path = out_dir / "synthetic24_regional_v3_assessment.md"
    report_path.write_text("\n".join(lines) + "\n")
    return report_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--skip-confusion", action="store_true")
    parser.add_argument("--reuse-confusion", action="store_true")
    args = parser.parse_args()

    base = args.base.resolve()
    out_dir = args.out_dir or (base / "regional_v3_assessment")
    out_dir.mkdir(parents=True, exist_ok=True)
    conditions, condition_map = load_manifest(base)
    names = [str(row["name"]) for row in conditions]

    metric_rows = []
    for condition in names:
        for run in RUNS:
            metric_rows.append(run_condition_metrics(base, condition, condition_map[condition], run))
    metrics = pd.DataFrame(metric_rows)

    confusion_path = out_dir / "confusion_long.tsv"
    confusion_rows = []
    if args.reuse_confusion and confusion_path.exists():
        confusion = pd.read_csv(confusion_path, sep="\t")
    elif not args.skip_confusion:
        for condition in names:
            for run in RUNS:
                print(f"confusion {run} {condition}", flush=True)
                confusion_rows.extend(condition_confusion_rows(base, condition, run))
        confusion = pd.DataFrame(confusion_rows)
    else:
        confusion = pd.DataFrame(confusion_rows)
    if confusion.empty:
        confusion_summary = pd.DataFrame()
    else:
        confusion_summary = summarize_confusion(confusion)
    comparison = comparison_table(metrics, confusion_summary) if not confusion_summary.empty else pd.DataFrame()

    metrics.to_csv(out_dir / "condition_metrics.tsv", sep="\t", index=False)
    if not confusion.empty:
        confusion.to_csv(out_dir / "confusion_long.tsv", sep="\t", index=False)
        confusion_summary.to_csv(out_dir / "confusion_summary.tsv", sep="\t", index=False)
    if not comparison.empty:
        comparison.to_csv(out_dir / "baseline_vs_regional_v3.tsv", sep="\t", index=False)
        report_path = write_report(out_dir, metrics, confusion, comparison)
        print(f"Report: {report_path}")
        print(
            "Summary: "
            f"max_abs_delta_gdna_rate={comparison['delta_gdna_locus_rate'].abs().max():.6f} "
            f"max_abs_delta_gdna_to_rna={comparison['delta_gdna_to_rna_rate'].abs().max():.6f} "
            f"max_abs_delta_spearman={comparison['delta_spearman'].abs().max():.6f} "
            f"max_abs_delta_mard={comparison['delta_mard'].abs().max():.6f}"
        )
    else:
        print(f"Metrics: {out_dir / 'condition_metrics.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())