#!/usr/bin/env python3
"""Pool-confusion analysis for the 24-condition synthetic Rigel suite."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pysam

from rigel.sim.manifest import condition_manifest_map, load_manifest
from rigel.sim.truth import parse_origin


AF_MRNA_BIT = 0x02
AF_GDNA_BIT = 0x04
AF_NRNA_BIT = 0x08
AF_INTERGENIC_BIT = 0x20

DEFAULT_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24")
DEFAULT_OUT = Path("results/synthetic_24_pool_confusion_2026-05-15")
DEFAULT_REPORT = Path("docs/benchmarks/synthetic_24_pool_confusion_2026-05-15.md")

POOLS = ["mrna", "nrna", "gdna", "unresolved"]


def safe_div(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if denominator else float("nan")


def pct(value: float, digits: int = 2) -> str:
    if value is None or math.isnan(float(value)):
        return "n/a"
    return f"{100.0 * float(value):.{digits}f}%"


def fmt(value: float, digits: int = 3) -> str:
    if value is None or math.isnan(float(value)):
        return "n/a"
    return f"{float(value):.{digits}f}"


def sci(value: float, digits: int = 3) -> str:
    if value is None or math.isnan(float(value)):
        return "n/a"
    return f"{float(value):.{digits}g}"


def fmt_count(value: float) -> str:
    if value is None or math.isnan(float(value)):
        return "n/a"
    return f"{float(value):,.0f}"


def read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".feather":
        return pd.read_feather(path)
    return pd.read_csv(path, sep="\t")


def table_path(out_dir: Path, stem: str) -> Path:
    feather_path = out_dir / f"{stem}.feather"
    if feather_path.exists():
        return feather_path
    tsv_path = out_dir / f"{stem}.tsv"
    if tsv_path.exists():
        return tsv_path
    raise FileNotFoundError(f"Missing {stem}.feather/.tsv in {out_dir}")


def get_tag(read: pysam.AlignedSegment, tag: str, default: Any) -> Any:
    try:
        return read.get_tag(tag)
    except KeyError:
        return default


def predicted_pool_from_zf(zf_value: int) -> str:
    if zf_value & AF_GDNA_BIT:
        return "gdna"
    if zf_value & AF_NRNA_BIT:
        return "nrna"
    if zf_value & AF_MRNA_BIT:
        return "mrna"
    return "unresolved"


def markdown_table(records: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    if not records:
        return "_No rows._"
    headers = [header for header, _ in columns]
    rows = [[str(record.get(key, "")) for _, key in columns] for record in records]
    widths = [len(header) for header in headers]
    for row in rows:
        for column_index, cell in enumerate(row):
            widths[column_index] = max(widths[column_index], len(cell))
    header_line = "| " + " | ".join(
        header.ljust(widths[column_index]) for column_index, header in enumerate(headers)
    ) + " |"
    separator = "| " + " | ".join("-" * width for width in widths) + " |"
    body = [
        "| " + " | ".join(
            cell.ljust(widths[column_index]) for column_index, cell in enumerate(row)
        ) + " |"
        for row in rows
    ]
    return "\n".join([header_line, separator, *body])


def load_truth_table(base: Path, meta: dict[str, Any]) -> pd.DataFrame:
    truth_name = str(meta.get("truth_abundances", "truth_abundances.tsv"))
    return pd.read_csv(base / truth_name, sep="\t")


def scan_annotated_bam(
    bam_path: Path,
    truth_gene_by_tx: dict[str, str],
    condition: str,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pool_counts: Counter[tuple[str, str]] = Counter()
    zc_counts: Counter[tuple[str, str, str, str]] = Counter()
    locus_counts: Counter[tuple[int, str, str]] = Counter()
    zf_counts: Counter[tuple[str, str, int]] = Counter()
    posterior_sums: Counter[tuple[str, str]] = Counter()
    metrics: Counter[str] = Counter()

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue

            origin = parse_origin(read.query_name)
            true_pool = origin.kind
            zf_value = int(get_tag(read, "ZF", 0))
            predicted_pool = predicted_pool_from_zf(zf_value)
            zc_value = str(get_tag(read, "ZC", "."))
            zs_value = str(get_tag(read, "ZS", "."))
            zw_value = float(get_tag(read, "ZW", float("nan")))
            locus_id = int(get_tag(read, "ZL", -1))
            assigned_tx = str(get_tag(read, "ZT", ""))
            assigned_gene = str(get_tag(read, "ZG", ""))

            metrics["total_fragments"] += 1
            metrics[f"truth_{true_pool}"] += 1
            metrics[f"pred_{predicted_pool}"] += 1
            pool_counts[(true_pool, predicted_pool)] += 1
            zc_counts[(true_pool, predicted_pool, zc_value, zs_value)] += 1
            locus_counts[(locus_id, true_pool, predicted_pool)] += 1
            zf_counts[(true_pool, predicted_pool, zf_value)] += 1
            if not math.isnan(zw_value):
                posterior_sums[(true_pool, predicted_pool)] += zw_value

            if true_pool == "mrna":
                true_tx = origin.transcript_id or ""
                true_gene = truth_gene_by_tx.get(true_tx, true_tx.rsplit(".", 1)[0])
                if predicted_pool == "mrna" and assigned_tx == true_tx:
                    metrics["mrna_exact_tx"] += 1
                elif predicted_pool == "mrna" and assigned_gene == true_gene:
                    metrics["mrna_same_gene_wrong_tx"] += 1
                elif predicted_pool == "mrna":
                    metrics["mrna_wrong_gene"] += 1

    total_fragments = int(metrics["total_fragments"])
    confusion_rows = []
    for true_pool in POOLS:
        truth_total = sum(pool_counts[(true_pool, predicted_pool)] for predicted_pool in POOLS)
        for predicted_pool in POOLS:
            count = int(pool_counts[(true_pool, predicted_pool)])
            posterior_key = (true_pool, predicted_pool)
            confusion_rows.append(
                {
                    "condition": condition,
                    "true_pool": true_pool,
                    "pred_pool": predicted_pool,
                    "count": count,
                    "row_fraction": safe_div(count, truth_total),
                    "overall_fraction": safe_div(count, total_fragments),
                    "mean_winner_posterior": safe_div(posterior_sums[posterior_key], count),
                }
            )
    confusion_df = pd.DataFrame(confusion_rows)

    zc_df = pd.DataFrame(
        [
            {
                "condition": condition,
                "true_pool": true_pool,
                "pred_pool": predicted_pool,
                "zc": zc_value,
                "zs": zs_value,
                "count": int(count),
            }
            for (true_pool, predicted_pool, zc_value, zs_value), count in zc_counts.items()
        ]
    )
    locus_df = pd.DataFrame(
        [
            {
                "condition": condition,
                "locus_id": int(locus_id),
                "true_pool": true_pool,
                "pred_pool": predicted_pool,
                "count": int(count),
            }
            for (locus_id, true_pool, predicted_pool), count in locus_counts.items()
        ]
    )
    zf_df = pd.DataFrame(
        [
            {
                "condition": condition,
                "true_pool": true_pool,
                "pred_pool": predicted_pool,
                "zf": int(zf_value),
                "count": int(count),
            }
            for (true_pool, predicted_pool, zf_value), count in zf_counts.items()
        ]
    )
    return dict(metrics), confusion_df, zc_df, locus_df, zf_df


def add_rates(row: dict[str, Any]) -> dict[str, Any]:
    truth_mrna = row.get("truth_mrna", 0)
    truth_nrna = row.get("truth_nrna", 0)
    truth_gdna = row.get("truth_gdna", 0)
    total = row.get("total_fragments", 0)

    for true_pool in ["mrna", "nrna", "gdna"]:
        truth_total = row.get(f"truth_{true_pool}", 0)
        for predicted_pool in ["mrna", "nrna", "gdna", "unresolved"]:
            key = f"{true_pool}_to_{predicted_pool}"
            row[f"{key}_rate"] = safe_div(row.get(key, 0), truth_total)

    row["assigned_mrna_fraction"] = safe_div(row.get("pred_mrna", 0), total)
    row["assigned_nrna_fraction"] = safe_div(row.get("pred_nrna", 0), total)
    row["assigned_gdna_fraction"] = safe_div(row.get("pred_gdna", 0), total)
    row["assigned_unresolved_fraction"] = safe_div(row.get("pred_unresolved", 0), total)
    row["truth_mrna_fraction"] = safe_div(truth_mrna, total)
    row["truth_nrna_fraction"] = safe_div(truth_nrna, total)
    row["truth_gdna_fraction"] = safe_div(truth_gdna, total)
    row["mrna_exact_tx_rate"] = safe_div(row.get("mrna_exact_tx", 0), truth_mrna)
    row["mrna_gene_rate"] = safe_div(
        row.get("mrna_exact_tx", 0) + row.get("mrna_same_gene_wrong_tx", 0),
        truth_mrna,
    )
    row["mrna_wrong_gene_rate"] = safe_div(row.get("mrna_wrong_gene", 0), truth_mrna)
    row["nrna_rna_compatible_rate"] = safe_div(
        row.get("nrna_to_nrna", 0) + row.get("nrna_to_mrna", 0),
        truth_nrna,
    )
    row["gdna_to_rna_rate"] = safe_div(
        row.get("gdna_to_mrna", 0) + row.get("gdna_to_nrna", 0),
        truth_gdna,
    )
    row["rna_to_gdna_rate"] = safe_div(
        row.get("mrna_to_gdna", 0) + row.get("nrna_to_gdna", 0),
        truth_mrna + truth_nrna,
    )
    row["overall_pool_accuracy"] = safe_div(
        row.get("mrna_to_mrna", 0) + row.get("nrna_to_nrna", 0) + row.get("gdna_to_gdna", 0),
        total,
    )
    return row


def summarize_condition(
    *,
    base: Path,
    condition: str,
    meta: dict[str, Any],
    out_name: str,
    annotated_name: str,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out_dir = base / condition / out_name
    bam_path = base / condition / annotated_name
    if not (out_dir / "summary.json").exists():
        raise FileNotFoundError(f"Missing Rigel output for {condition}: {out_dir}")
    if not bam_path.exists():
        raise FileNotFoundError(f"Missing annotated BAM for {condition}: {bam_path}")

    truth = load_truth_table(base, meta)
    truth_gene_by_tx = dict(zip(truth["transcript_id"], truth["gene_id"]))
    scan_metrics, confusion_df, zc_df, locus_df, zf_df = scan_annotated_bam(
        bam_path,
        truth_gene_by_tx,
        condition,
    )

    summary = read_json(out_dir / "summary.json")
    loci = read_table(table_path(out_dir, "loci"))
    quant = read_table(table_path(out_dir, "quant"))
    nrna_quant = read_table(table_path(out_dir, "nrna_quant"))
    gene_quant = read_table(table_path(out_dir, "gene_quant"))

    calibration = summary.get("calibration", {})
    densities = calibration.get("global_densities", {})
    intergenic = densities.get("INTERGENIC", {})
    intron = densities.get("INTRON", {})
    exon_intron = densities.get("EXON-INTRON", {})
    fl_models = calibration.get("fl_models", {})
    strand_model = summary.get("strand_model", {})

    truth_mrna = int(meta.get("n_mrna", scan_metrics.get("truth_mrna", 0)))
    truth_nrna = int(meta.get("n_nrna", scan_metrics.get("truth_nrna", 0)))
    truth_gdna = int(meta.get("n_gdna", scan_metrics.get("truth_gdna", 0)))
    total_truth = truth_mrna + truth_nrna + truth_gdna

    row: dict[str, Any] = {
        "condition": condition,
        "gdna_label": meta.get("gdna_label"),
        "gdna_rate": float(meta.get("gdna_rate", np.nan)),
        "nrna_label": meta.get("nrna_label"),
        "strand_specificity": float(meta.get("strand_specificity", np.nan)),
        "truth_mrna": truth_mrna,
        "truth_nrna": truth_nrna,
        "truth_gdna": truth_gdna,
        "truth_total_manifest": int(meta.get("n_total", total_truth)),
        "loci_mrna": float(loci["mrna"].sum()) if "mrna" in loci.columns else np.nan,
        "loci_nrna": float(loci["nrna"].sum()) if "nrna" in loci.columns else np.nan,
        "loci_gdna": float(loci["gdna"].sum()) if "gdna" in loci.columns else np.nan,
        "quant_mrna_total": float(quant["count"].sum()) if "count" in quant.columns else np.nan,
        "nrna_quant_total": float(nrna_quant["count"].sum()) if "count" in nrna_quant.columns else 0.0,
        "gene_quant_total": float(gene_quant["count"].sum()) if "count" in gene_quant.columns else np.nan,
        "strand_model_estimate": float(strand_model.get("strand_specificity", np.nan)),
        "rna_fl_mean": float(fl_models.get("rna_fl_mean", np.nan)),
        "gdna_fl_mean": float(fl_models.get("gdna_fl_mean", np.nan)),
        "gdna_fl_rel_error": safe_div(float(fl_models.get("gdna_fl_mean", np.nan)) - 350.0, 350.0),
        "rho_intergenic": float(intergenic.get("rho", np.nan)),
        "rho_intron": float(intron.get("rho", np.nan)),
        "rho_exon_intron": float(exon_intron.get("rho", np.nan)),
        "rho_intron_uncorrected": float(intron.get("rho_uncorrected", np.nan)),
        "rho_exon_intron_uncorrected": float(exon_intron.get("rho_uncorrected", np.nan)),
        "rho_intron_estimated_fragments": float(intron.get("n_fragments_estimated", np.nan)),
        "rho_exon_intron_estimated_fragments": float(exon_intron.get("n_fragments_estimated", np.nan)),
        "rho_intron_strand_active": bool(intron.get("strand_active", False)),
        "rho_exon_intron_strand_active": bool(exon_intron.get("strand_active", False)),
        "rho_intron_over_intergenic": safe_div(
            float(intron.get("rho", np.nan)),
            float(intergenic.get("rho", np.nan)),
        ),
        "rho_exon_intron_over_intergenic": safe_div(
            float(exon_intron.get("rho", np.nan)),
            float(intergenic.get("rho", np.nan)),
        ),
    }
    row.update(scan_metrics)
    for true_pool in ["mrna", "nrna", "gdna"]:
        for predicted_pool in POOLS:
            row[f"{true_pool}_to_{predicted_pool}"] = int(
                confusion_df.loc[
                    (confusion_df["true_pool"] == true_pool)
                    & (confusion_df["pred_pool"] == predicted_pool),
                    "count",
                ].sum()
            )
    row = add_rates(row)
    row["loci_mrna_recovery"] = safe_div(row["loci_mrna"], truth_mrna)
    row["loci_nrna_recovery"] = safe_div(row["loci_nrna"], truth_nrna)
    row["loci_gdna_recovery"] = safe_div(row["loci_gdna"], truth_gdna)
    row["winner_nrna_recall"] = row.get("nrna_to_nrna_rate")
    row["winner_gdna_recall"] = row.get("gdna_to_gdna_rate")
    row["winner_mrna_recall"] = row.get("mrna_to_mrna_rate")
    return row, confusion_df, zc_df, locus_df, zf_df


def build_worst_locus_table(locus_df: pd.DataFrame, base: Path, out_name: str) -> pd.DataFrame:
    if locus_df.empty:
        return pd.DataFrame()
    rows = []
    for condition, condition_locus in locus_df.groupby("condition"):
        loci_path = table_path(base / condition / out_name, "loci")
        loci = read_table(loci_path)
        useful_cols = [
            column
            for column in [
                "locus_id",
                "locus_span_bp",
                "n_transcripts",
                "n_nrna_entities",
                "n_genes",
                "mrna",
                "nrna",
                "gdna",
                "gdna_rate",
                "gdna_prior",
                "gdna_eff_len_per_bp",
            ]
            if column in loci.columns
        ]
        errors = condition_locus[
            (condition_locus["locus_id"] >= 0)
            & (condition_locus["true_pool"] != condition_locus["pred_pool"])
        ]
        if errors.empty:
            continue
        error_summary = errors.groupby(["condition", "locus_id"], as_index=False)["count"].sum()
        error_summary = error_summary.rename(columns={"count": "pool_error_count"})
        locus_truth = condition_locus.pivot_table(
            index=["condition", "locus_id"],
            columns="true_pool",
            values="count",
            aggfunc="sum",
            fill_value=0,
        ).reset_index()
        for pool in ["mrna", "nrna", "gdna"]:
            if pool not in locus_truth.columns:
                locus_truth[pool] = 0
        merged = error_summary.merge(locus_truth, on=["condition", "locus_id"], how="left")
        merged = merged.merge(loci[useful_cols], on="locus_id", how="left", suffixes=("_truth", ""))
        merged["pool_error_fraction"] = merged["pool_error_count"] / merged[
            ["mrna_truth", "nrna_truth", "gdna_truth"]
        ].sum(axis=1).replace(0, np.nan)
        rows.append(merged)
    if not rows:
        return pd.DataFrame()
    result = pd.concat(rows, ignore_index=True)
    return result.sort_values("pool_error_count", ascending=False)


def build_no_locus_error_table(locus_df: pd.DataFrame) -> pd.DataFrame:
    """Summarize false assignments stamped with ZL=-1 in the annotated BAM.

    In current annotated BAMs, gDNA winners are stamped without a locus id because
    their winning transcript id is the gDNA sentinel. This table makes that bucket
    explicit rather than presenting it as a biological locus.
    """
    if locus_df.empty:
        return pd.DataFrame()
    subset = locus_df[(locus_df["locus_id"] < 0) & (locus_df["true_pool"] != locus_df["pred_pool"])]
    if subset.empty:
        return pd.DataFrame()
    grouped = subset.groupby(["condition", "true_pool", "pred_pool"], as_index=False)[
        "count"
    ].sum()
    return grouped.sort_values("count", ascending=False)


def build_strand_pair_delta(condition_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metric_columns = [
        "overall_pool_accuracy",
        "mrna_to_mrna_rate",
        "mrna_to_gdna_rate",
        "mrna_to_nrna_rate",
        "nrna_to_nrna_rate",
        "nrna_to_mrna_rate",
        "nrna_to_gdna_rate",
        "nrna_rna_compatible_rate",
        "gdna_to_gdna_rate",
        "gdna_to_rna_rate",
        "assigned_gdna_fraction",
        "assigned_nrna_fraction",
        "loci_gdna_recovery",
        "loci_nrna_recovery",
        "gdna_fl_mean",
        "rho_intron_over_intergenic",
        "rho_exon_intron_over_intergenic",
    ]
    for (gdna_label, nrna_label), group in condition_df.groupby(["gdna_label", "nrna_label"]):
        unstranded = group[group["strand_specificity"] < 0.75]
        stranded = group[group["strand_specificity"] > 0.75]
        if unstranded.empty or stranded.empty:
            continue
        left = unstranded.iloc[0]
        right = stranded.iloc[0]
        row: dict[str, Any] = {
            "gdna_label": gdna_label,
            "nrna_label": nrna_label,
            "unstranded_condition": left["condition"],
            "stranded_condition": right["condition"],
            "truth_mrna": right["truth_mrna"],
            "truth_nrna": right["truth_nrna"],
            "truth_gdna": right["truth_gdna"],
        }
        for metric in metric_columns:
            row[f"unstranded_{metric}"] = left.get(metric, np.nan)
            row[f"stranded_{metric}"] = right.get(metric, np.nan)
            row[f"delta_{metric}"] = right.get(metric, np.nan) - left.get(metric, np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def condition_report_rows(condition_df: pd.DataFrame) -> list[dict[str, str]]:
    sorted_df = condition_df.sort_values(["nrna_label", "gdna_rate", "strand_specificity"])
    rows = []
    for row in sorted_df.itertuples(index=False):
        rows.append(
            {
                "condition": row.condition,
                "truth_gdna": pct(row.truth_gdna_fraction),
                "pred_gdna": pct(row.assigned_gdna_fraction),
                "truth_nrna": pct(row.truth_nrna_fraction),
                "nrna_to_nrna": pct(row.nrna_to_nrna_rate),
                "nrna_to_gdna": pct(row.nrna_to_gdna_rate),
                "gdna_to_gdna": pct(row.gdna_to_gdna_rate),
                "gdna_to_rna": pct(row.gdna_to_rna_rate),
                "mrna_to_gdna": pct(row.mrna_to_gdna_rate),
                "pool_acc": pct(row.overall_pool_accuracy),
            }
        )
    return rows


def strand_delta_report_rows(delta_df: pd.DataFrame) -> list[dict[str, str]]:
    rows = []
    if delta_df.empty:
        return rows
    sorted_df = delta_df.sort_values(["nrna_label", "gdna_label"])
    for row in sorted_df.itertuples(index=False):
        rows.append(
            {
                "gdna": row.gdna_label,
                "nrna": row.nrna_label,
                "d_pool_acc": pct(row.delta_overall_pool_accuracy),
                "d_nrna_recall": pct(row.delta_nrna_to_nrna_rate),
                "d_nrna_to_gdna": pct(row.delta_nrna_to_gdna_rate),
                "d_gdna_recall": pct(row.delta_gdna_to_gdna_rate),
                "d_gdna_to_rna": pct(row.delta_gdna_to_rna_rate),
                "d_assigned_gdna": pct(row.delta_assigned_gdna_fraction),
            }
        )
    return rows


def zc_error_report_rows(zc_df: pd.DataFrame, condition: str, max_rows: int = 8) -> list[dict[str, str]]:
    errors = zc_df[(zc_df["condition"] == condition) & (zc_df["true_pool"] != zc_df["pred_pool"])]
    if errors.empty:
        return []
    grouped = errors.groupby(["true_pool", "pred_pool", "zc", "zs"], as_index=False)["count"].sum()
    grouped = grouped.sort_values("count", ascending=False).head(max_rows)
    return [
        {
            "true": row.true_pool,
            "pred": row.pred_pool,
            "ZC": row.zc,
            "ZS": row.zs,
            "count": fmt_count(row.count),
        }
        for row in grouped.itertuples(index=False)
    ]


def calibration_focus_rows(condition_df: pd.DataFrame) -> list[dict[str, str]]:
    focus = condition_df[condition_df["nrna_label"] == "high"].sort_values(
        ["gdna_rate", "strand_specificity"]
    )
    rows = []
    for row in focus.itertuples(index=False):
        rows.append(
            {
                "condition": row.condition,
                "ss": fmt(row.strand_specificity, 2),
                "gdna_fl": fmt(row.gdna_fl_mean, 1),
                "rho_ig": sci(row.rho_intergenic),
                "rho_in": sci(row.rho_intron),
                "rho_ex_in": sci(row.rho_exon_intron),
                "pred_gdna": pct(row.assigned_gdna_fraction),
                "nrna_to_gdna": pct(row.nrna_to_gdna_rate),
            }
        )
    return rows


def write_report(
    *,
    report_path: Path,
    base: Path,
    out_name: str,
    annotated_name: str,
    out_dir: Path,
    condition_df: pd.DataFrame,
    confusion_df: pd.DataFrame,
    zc_df: pd.DataFrame,
    strand_delta_df: pd.DataFrame,
    worst_locus_df: pd.DataFrame,
    no_locus_error_df: pd.DataFrame,
) -> None:
    stranded_df = condition_df[condition_df["strand_specificity"] > 0.75]
    unstranded_df = condition_df[condition_df["strand_specificity"] < 0.75]

    def pick_condition(condition_name: str) -> pd.Series:
        matches = condition_df[condition_df["condition"] == condition_name]
        if not matches.empty:
            return matches.iloc[0]
        return condition_df.iloc[0]

    zero_gdna_high_stranded = pick_condition("gdna_zero_ss_0.99_nrna_high")
    zero_gdna_high_unstranded = pick_condition("gdna_zero_ss_0.50_nrna_high")

    stranded_gdna = stranded_df[stranded_df["truth_gdna"] > 0]
    stranded_nrna = stranded_df[stranded_df["truth_nrna"] > 0]
    zero_nrna_stranded = stranded_df[stranded_df["nrna_label"] == "zero"]

    mean_stranded_gdna_recall = stranded_gdna["gdna_to_gdna_rate"].mean()
    mean_stranded_gdna_to_rna = stranded_gdna["gdna_to_rna_rate"].mean()
    mean_stranded_nrna_recall = stranded_nrna["nrna_to_nrna_rate"].mean()
    mean_stranded_nrna_to_gdna = stranded_nrna["nrna_to_gdna_rate"].mean()
    mean_unstranded_nrna_recall = unstranded_df[unstranded_df["truth_nrna"] > 0][
        "nrna_to_nrna_rate"
    ].mean()
    mean_unstranded_nrna_to_gdna = unstranded_df[unstranded_df["truth_nrna"] > 0][
        "nrna_to_gdna_rate"
    ].mean()
    mean_zero_nrna_stranded_mrna_to_gdna = zero_nrna_stranded["mrna_to_gdna_rate"].mean()

    worst_locus_rows = []
    if not worst_locus_df.empty:
        focus = worst_locus_df.head(12)
        for row in focus.itertuples(index=False):
            worst_locus_rows.append(
                {
                    "condition": row.condition,
                    "locus": str(int(row.locus_id)),
                    "errors": fmt_count(row.pool_error_count),
                    "err_frac": pct(row.pool_error_fraction),
                    "truth_mrna": fmt_count(getattr(row, "mrna_truth", np.nan)),
                    "truth_nrna": fmt_count(getattr(row, "nrna_truth", np.nan)),
                    "truth_gdna": fmt_count(getattr(row, "gdna_truth", np.nan)),
                    "pred_gdna_rate": pct(getattr(row, "gdna_rate", np.nan)),
                    "n_tx": fmt_count(getattr(row, "n_transcripts", np.nan)),
                }
            )

    no_locus_rows = []
    if not no_locus_error_df.empty:
        for row in no_locus_error_df.head(12).itertuples(index=False):
            no_locus_rows.append(
                {
                    "condition": row.condition,
                    "true": row.true_pool,
                    "pred": row.pred_pool,
                    "count": fmt_count(row.count),
                }
            )

    condition_rows = condition_report_rows(condition_df)
    delta_rows = strand_delta_report_rows(strand_delta_df)
    calibration_rows = calibration_focus_rows(condition_df)
    zero_high_error_rows = zc_error_report_rows(zc_df, "gdna_zero_ss_0.99_nrna_high")
    high_gdna_error_rows = zc_error_report_rows(zc_df, "gdna_high_ss_0.99_nrna_high")

    report = f"""# Synthetic 24 Pool Confusion Report - 2026-05-15

Base directory: `{base}`

Rigel output namespace: `<condition>/{out_name}`

Annotated BAM namespace: `<condition>/{annotated_name}`

Analysis artifacts: `{out_dir}`

## Executive Summary

This rerun completed all 24 synthetic scenarios with the current effective-length
normalization code. The pool-level confusion matrices use oracle read names as truth
and Rigel's annotated-BAM `ZF` winner flags as the predicted pool. Counts are therefore
hard winner assignments; the fractional EM pool masses are also written to the metrics
table from `loci.feather`.

In strand-specific cases (`ss=0.99`), Rigel recovers true gDNA well when gDNA is present:
mean gDNA winner recall across nonzero-gDNA stranded scenarios is
{pct(mean_stranded_gdna_recall)}, with mean gDNA-to-RNA leakage of
{pct(mean_stranded_gdna_to_rna)}. Mature RNA is also rarely mistaken for gDNA in clean
stranded, no-nRNA settings: mean mRNA-to-gDNA is
{pct(mean_zero_nrna_stranded_mrna_to_gdna)}.

The remaining hard problem is not ordinary mature mRNA versus gDNA. It is genic,
unspliced nRNA versus genic gDNA. Across stranded scenarios with true nRNA, mean nRNA
winner recall is {pct(mean_stranded_nrna_recall)} and mean nRNA-to-gDNA leakage is
{pct(mean_stranded_nrna_to_gdna)}. This is a major improvement over unstranded nRNA
scenarios, where mean nRNA recall is {pct(mean_unstranded_nrna_recall)} and mean
nRNA-to-gDNA leakage is {pct(mean_unstranded_nrna_to_gdna)}, but it is still not the
dramatic separation we want.

The zero-gDNA, high-nRNA contrast is the cleanest diagnostic. In the unstranded run,
Rigel sends {pct(zero_gdna_high_unstranded.nrna_to_gdna_rate)} of true nRNA fragments
to gDNA. In the stranded run, this falls to
{pct(zero_gdna_high_stranded.nrna_to_gdna_rate)}, and nRNA recall rises from
{pct(zero_gdna_high_unstranded.nrna_to_nrna_rate)} to
{pct(zero_gdna_high_stranded.nrna_to_nrna_rate)}. Strand information is being used, but
genic unspliced evidence remains partially gDNA-like.

## Per-Scenario Pool Confusion

The table below reports row-normalized winner assignments for each true pool. `nRNA to
gDNA` means true nascent fragments assigned to the gDNA pool. `gDNA to RNA` combines
gDNA assigned to mature or nascent RNA.

{markdown_table(condition_rows, [
        ("Condition", "condition"),
        ("Truth gDNA", "truth_gdna"),
        ("Pred gDNA", "pred_gdna"),
        ("Truth nRNA", "truth_nrna"),
        ("nRNA to nRNA", "nrna_to_nrna"),
        ("nRNA to gDNA", "nrna_to_gdna"),
        ("gDNA to gDNA", "gdna_to_gdna"),
        ("gDNA to RNA", "gdna_to_rna"),
        ("mRNA to gDNA", "mrna_to_gdna"),
        ("Pool acc", "pool_acc"),
    ])}

## Strand-Specific Delta

Each row compares `ss=0.99` minus `ss=0.50` for the same gDNA/nRNA setting. Positive
`d nRNA recall` and negative `d nRNA to gDNA` are improvements.

{markdown_table(delta_rows, [
        ("gDNA", "gdna"),
        ("nRNA", "nrna"),
        ("d pool acc", "d_pool_acc"),
        ("d nRNA recall", "d_nrna_recall"),
        ("d nRNA to gDNA", "d_nrna_to_gdna"),
        ("d gDNA recall", "d_gdna_recall"),
        ("d gDNA to RNA", "d_gdna_to_rna"),
        ("d assigned gDNA", "d_assigned_gdna"),
    ])}

## What The Strand Signal Fixes

Strand-specific data sharply reduces intronic nRNA masquerading as gDNA. nRNA is
unspliced and spans introns, while mature mRNA is confined to exons and exon-exon splice
junctions. In stranded data, opposite-strand genic-unspliced evidence can be discounted
as gDNA-like, so the calibration correction suppresses much of the intron-only false
gDNA signal. This is visible in the zero-gDNA high-nRNA pair above and in the large
negative `d nRNA to gDNA` values for low/high nRNA settings.

## What Remains Confused

The strand signal does not fully resolve genic unspliced fragments that lack splice
junction evidence. Some are exon-intron boundary fragments that can feed the
EXON-INTRON calibration channel; others are exon-confined mature-RNA fragments, which
are compatible with both a mature exon and the corresponding nascent transcript span.
In zero-gDNA high-nRNA stranded data, the top error categories are:

{markdown_table(zero_high_error_rows, [
        ("True", "true"),
        ("Pred", "pred"),
        ("ZC", "ZC"),
        ("ZS", "ZS"),
        ("Count", "count"),
    ])}

When true gDNA and high nRNA coexist in stranded data, the largest residual errors are:

{markdown_table(high_gdna_error_rows, [
        ("True", "true"),
        ("Pred", "pred"),
        ("ZC", "ZC"),
        ("ZS", "ZS"),
        ("Count", "count"),
    ])}

## Calibration Diagnostics In High-nRNA Scenarios

When true gDNA is absent or low, the gDNA FL model and genic gDNA densities can be
pulled toward nRNA-like evidence. `rho_ig` is the intergenic anchor; `rho_in` and
`rho_ex_in` are genic unspliced channels.

{markdown_table(calibration_rows, [
        ("Condition", "condition"),
        ("SS", "ss"),
        ("gDNA FL", "gdna_fl"),
        ("rho ig", "rho_ig"),
        ("rho intron", "rho_in"),
        ("rho exon-intron", "rho_ex_in"),
        ("Pred gDNA", "pred_gdna"),
        ("nRNA to gDNA", "nrna_to_gdna"),
    ])}

## Predicted-gDNA / No-transcript Bucket

Current annotated BAMs stamp gDNA winners with `ZT=.` and `ZL=-1`, so false gDNA
assignments cannot yet be localized back to the EM locus from the BAM alone. The table
below separates that diagnostic bucket from true biological loci. Preserving the source
locus id for gDNA winners would make future root-cause analysis much sharper.

{markdown_table(no_locus_rows, [
        ("Condition", "condition"),
        ("True", "true"),
        ("Pred", "pred"),
        ("Count", "count"),
    ])}

## Worst EM Loci By Non-gDNA Pool Error

{markdown_table(worst_locus_rows, [
        ("Condition", "condition"),
        ("Locus", "locus"),
        ("Errors", "errors"),
        ("Err frac", "err_frac"),
        ("Truth mRNA", "truth_mrna"),
        ("Truth nRNA", "truth_nrna"),
        ("Truth gDNA", "truth_gdna"),
        ("Pred gDNA rate", "pred_gdna_rate"),
        ("n tx", "n_tx"),
    ])}

## Root-Cause Interpretation

1. Strand-specific data works as intended for mature RNA versus gDNA. Mature fragments
   carry exonic or splice-junction evidence, and true gDNA has intergenic support plus
   continuous genomic compatibility. In stranded no-nRNA cases, mRNA-to-gDNA leakage is
   very small.
2. The main residual ambiguity is genic unspliced signal. nRNA and genic gDNA are both
    continuous over the transcript span and both cover introns. Mature mRNA also produces
    unspliced exon-confined fragments that are compatible with nRNA unless a splice
    junction or intronic evidence separates them. Strand correction removes much of the
    intron-only nRNA-as-gDNA signal, but it does not create a positive nRNA state during
    calibration; remaining genic-unspliced signal can still seed gDNA priors.
3. The gDNA fragment-length model is vulnerable when nRNA is abundant. If genic
   unspliced nRNA is admitted into gDNA training sources, the fitted gDNA FL can become
   RNA-like, which then reinforces nRNA/gDNA confusion in EM.
4. Pool-level mRNA behavior is better than transcript-level isoform behavior. Many mRNA
   mistakes are same-gene or same-locus isoform allocation issues rather than wrong-pool
   errors. That is a separate identifiability problem from the nRNA/gDNA confusion.

## Further Exploration

1. Build per-locus nRNA/gDNA truth summaries for the worst loci above and inspect the
   local transcript structures, especially loci with many nRNA entities or shared
   transcript spans.
2. Stratify nRNA errors by read geometry: intron-only, exon-intron boundary,
   exon-confined unspliced, and splice-junction-supported. The BAM `ZC` and `ZS` tags are
   a start, but a dedicated geometry label would make this much sharper.
3. Compare gDNA FL estimates trained from intergenic-only fragments versus the current
   pooled gDNA source in high-nRNA stranded cases.
4. Add synthetic fixtures for zero-gDNA plus high-nRNA at `ss=0.99`; this should be a
   regression guard for false gDNA from nascent RNA.
5. Preserve the source locus id for gDNA EM winners in annotated BAM output so false
    gDNA calls can be traced directly to their competing transcript/nRNA components.

## Improvements To Consider

1. Add an explicit nRNA-aware calibration state before projecting genic-unspliced signal
   into gDNA priors. Intergenic evidence should anchor true gDNA; genic intron/exon
   signal without intergenic support should not automatically become gDNA.
2. Train or regularize the gDNA FL model primarily from intergenic-only fragments when
   enough evidence exists, then use genic channels only after strand/nRNA correction.
3. Gate genic gDNA priors by consistency with intergenic density. If INTERGENIC density
   is near zero but INTRON or EXON-INTRON density is high, prefer nRNA or downweight the
   gDNA prior.
4. Introduce a locus-level competition feature between nRNA span coverage and gDNA: nRNA
   should be favored when unspliced coverage is transcript-span-local and lacks flanking
   intergenic support.
5. Report pool confusion diagnostics from annotated BAMs as a first-class benchmark
   artifact, because aggregate quant tables can hide where fragments moved.
6. Preserve locus ids for gDNA winners in annotated BAMs; this is a diagnostics
    improvement rather than a model change, but it would substantially shorten future
    root-cause analysis loops.

## Artifacts

- Condition metrics: `{out_dir / 'condition_metrics.tsv'}`
- Pool confusion counts: `{out_dir / 'pool_confusion_counts.tsv'}`
- Pool confusion matrices: `{out_dir / 'pool_confusion_matrices.tsv'}`
- Strand-pair deltas: `{out_dir / 'strand_pair_delta.tsv'}`
- ZC/ZS breakdown: `{out_dir / 'zc_zs_breakdown.tsv'}`
- ZF breakdown: `{out_dir / 'zf_breakdown.tsv'}`
- Locus pool breakdown: `{out_dir / 'locus_pool_breakdown.tsv'}`
- Predicted-gDNA/no-transcript bucket: `{out_dir / 'no_locus_pool_errors.tsv'}`
- Worst EM locus errors: `{out_dir / 'worst_locus_pool_errors.tsv'}`
"""
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--out-name", default="rigel_strand_v2_out")
    parser.add_argument("--annotated-name", default="annotated_strand_v2.bam")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--conditions", nargs="*", default=None)
    args = parser.parse_args()

    base = args.base.resolve()
    out_dir = args.out_dir.resolve()
    report_path = args.report.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = load_manifest(base)
    condition_meta = condition_manifest_map(manifest)
    conditions = args.conditions or list(condition_meta)

    condition_rows: list[dict[str, Any]] = []
    confusion_tables: list[pd.DataFrame] = []
    zc_tables: list[pd.DataFrame] = []
    locus_tables: list[pd.DataFrame] = []
    zf_tables: list[pd.DataFrame] = []

    for condition in conditions:
        meta = condition_meta[condition]
        row, confusion_df, zc_df, locus_df, zf_df = summarize_condition(
            base=base,
            condition=condition,
            meta=meta,
            out_name=args.out_name,
            annotated_name=args.annotated_name,
        )
        condition_rows.append(row)
        confusion_tables.append(confusion_df)
        zc_tables.append(zc_df)
        locus_tables.append(locus_df)
        zf_tables.append(zf_df)
        print(f"analyzed {condition}", flush=True)

    condition_df = pd.DataFrame(condition_rows)
    confusion_df = pd.concat(confusion_tables, ignore_index=True)
    zc_df = pd.concat(zc_tables, ignore_index=True) if zc_tables else pd.DataFrame()
    locus_df = pd.concat(locus_tables, ignore_index=True) if locus_tables else pd.DataFrame()
    zf_df = pd.concat(zf_tables, ignore_index=True) if zf_tables else pd.DataFrame()
    strand_delta_df = build_strand_pair_delta(condition_df)
    worst_locus_df = build_worst_locus_table(locus_df, base, args.out_name)
    no_locus_error_df = build_no_locus_error_table(locus_df)

    matrix_df = confusion_df.pivot_table(
        index=["condition", "true_pool"],
        columns="pred_pool",
        values="count",
        aggfunc="sum",
        fill_value=0,
    ).reset_index()

    condition_df.to_csv(out_dir / "condition_metrics.tsv", sep="\t", index=False)
    confusion_df.to_csv(out_dir / "pool_confusion_counts.tsv", sep="\t", index=False)
    matrix_df.to_csv(out_dir / "pool_confusion_matrices.tsv", sep="\t", index=False)
    strand_delta_df.to_csv(out_dir / "strand_pair_delta.tsv", sep="\t", index=False)
    zc_df.to_csv(out_dir / "zc_zs_breakdown.tsv", sep="\t", index=False)
    zf_df.to_csv(out_dir / "zf_breakdown.tsv", sep="\t", index=False)
    locus_df.to_csv(out_dir / "locus_pool_breakdown.tsv", sep="\t", index=False)
    worst_locus_df.to_csv(out_dir / "worst_locus_pool_errors.tsv", sep="\t", index=False)
    no_locus_error_df.to_csv(out_dir / "no_locus_pool_errors.tsv", sep="\t", index=False)

    write_report(
        report_path=report_path,
        base=base,
        out_name=args.out_name,
        annotated_name=args.annotated_name,
        out_dir=out_dir,
        condition_df=condition_df,
        confusion_df=confusion_df,
        zc_df=zc_df,
        strand_delta_df=strand_delta_df,
        worst_locus_df=worst_locus_df,
        no_locus_error_df=no_locus_error_df,
    )
    print(f"wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())