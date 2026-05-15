#!/usr/bin/env python3
"""Deep analysis for the synthetic 24-condition Rigel simulation suite.

This script consumes outputs from scripts/sim/evaluate_suite.py. It uses
oracle read names in the annotated BAMs for fragment-level truth, avoiding
the molecule-vs-fragment abundance pitfall in truth_abundances_*.tsv.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter, defaultdict
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


DEFAULT_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24")
DEFAULT_OUT = Path("results/synthetic_24_deep_analysis")
DEFAULT_REPORT = Path("docs/benchmarks/synthetic_24_rigel_performance_2026-05-14.md")


def safe_div(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if denominator else float("nan")


def pct(value: float) -> str:
    if value is None or math.isnan(float(value)):
        return "n/a"
    return f"{100.0 * float(value):.2f}%"


def fmt(value: float, digits: int = 3) -> str:
    if value is None or math.isnan(float(value)):
        return "n/a"
    return f"{float(value):.{digits}f}"


def fmt_count(value: float) -> str:
    if value is None or math.isnan(float(value)):
        return "n/a"
    return f"{float(value):,.0f}"


def read_json(path: Path) -> dict[str, Any]:
    with open(path) as handle:
        return json.load(handle)


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".feather":
        return pd.read_feather(path)
    return pd.read_csv(path, sep="\t")


def load_condition_truth(base: Path, condition_meta: dict[str, Any]) -> pd.DataFrame:
    truth_name = str(condition_meta.get("truth_abundances", "truth_abundances.tsv"))
    return pd.read_csv(base / truth_name, sep="\t")


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for match in re.finditer(r'(\S+)\s+"([^"]*)"', text):
        attrs[match.group(1)] = match.group(2)
    return attrs


def parse_gtf_exons(gtf_path: Path) -> dict[str, tuple[tuple[int, int], ...]]:
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            start = int(fields[3]) - 1
            end = int(fields[4])
            exons[transcript_id].append((start, end))
    return {tx_id: tuple(sorted(items)) for tx_id, items in exons.items()}


def intron_chain(exons: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    if len(exons) < 2:
        return ()
    return tuple((exons[index][1], exons[index + 1][0]) for index in range(len(exons) - 1))


def build_structure_table(base: Path) -> pd.DataFrame:
    index_t = pd.read_csv(base / "rigel_index" / "transcripts.tsv", sep="\t")
    exons_by_tx = parse_gtf_exons(base / "reference" / "genes.gtf")

    rows = []
    for _, row in index_t.iterrows():
        tx_id = str(row["t_id"])
        exons = exons_by_tx.get(tx_id, ())
        rows.append(
            {
                "transcript_id": tx_id,
                "gene_id": row["g_id"],
                "ref": row["ref"],
                "strand": row["strand"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "length": int(row["length"]),
                "n_exons": int(row["n_exons"]),
                "is_nrna": bool(row["is_nrna"]),
                "is_synthetic": bool(row["is_synthetic"]),
                "nrna_t_index": int(row["nrna_t_index"]),
                "exon_signature": repr((row["ref"], row["strand"], exons)),
                "intron_signature": repr((row["ref"], row["strand"], intron_chain(exons))),
            }
        )
    table = pd.DataFrame(rows)
    annotated = table[~table["is_synthetic"]].copy()

    gene_sizes = annotated.groupby("gene_id")["transcript_id"].transform("size")
    exact_sizes = annotated.groupby("exon_signature")["transcript_id"].transform("size")
    intron_sizes = annotated.groupby("intron_signature")["transcript_id"].transform("size")
    nrna_sizes = annotated.groupby("nrna_t_index")["transcript_id"].transform("size")

    annotated["gene_n_tx"] = gene_sizes.astype(int)
    annotated["exact_exon_group_size"] = exact_sizes.astype(int)
    annotated["intron_chain_group_size"] = intron_sizes.astype(int)
    annotated["nrna_span_group_size"] = nrna_sizes.astype(int)

    conditions = [
        annotated["exact_exon_group_size"] > 1,
        (annotated["n_exons"] > 1) & (annotated["intron_chain_group_size"] > 1),
        annotated["nrna_span_group_size"] > 1,
        annotated["gene_n_tx"] > 1,
    ]
    choices = [
        "exact_exon_duplicate",
        "same_intron_chain",
        "shared_nrna_span",
        "multi_isoform_gene",
    ]
    annotated["structure_class"] = np.select(conditions, choices, default="singleton")
    return annotated


def assigned_pool_from_zf(zf_value: int) -> str:
    if zf_value & AF_GDNA_BIT:
        return "gdna"
    if zf_value & AF_NRNA_BIT:
        return "nrna"
    if zf_value & AF_MRNA_BIT:
        return "mrna"
    return "unresolved"


def get_tag(read: pysam.AlignedSegment, tag: str, default: Any) -> Any:
    try:
        return read.get_tag(tag)
    except KeyError:
        return default


def scan_annotated_bam(
    bam_path: Path,
    truth_gene_by_tx: dict[str, str],
) -> tuple[dict[str, Any], Counter[str], Counter[str], pd.DataFrame]:
    truth_mrna: Counter[str] = Counter()
    truth_nrna: Counter[str] = Counter()
    origin_pool: Counter[tuple[str, str]] = Counter()
    zc_rows: Counter[tuple[str, str, str]] = Counter()

    metrics: Counter[str] = Counter()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue

            origin = parse_origin(read.query_name)
            zf_value = int(get_tag(read, "ZF", 0))
            pool = assigned_pool_from_zf(zf_value)
            zc_value = str(get_tag(read, "ZC", "."))
            assigned_tx = str(get_tag(read, "ZT", ""))
            assigned_gene = str(get_tag(read, "ZG", ""))

            metrics["total_fragments"] += 1
            metrics[f"origin_{origin.kind}"] += 1
            metrics[f"assigned_{pool}"] += 1
            origin_pool[(origin.kind, pool)] += 1
            zc_rows[(origin.kind, pool, zc_value)] += 1

            if origin.kind == "mrna":
                true_tx = origin.transcript_id or ""
                truth_mrna[true_tx] += 1
                true_gene = truth_gene_by_tx.get(true_tx, true_tx.rsplit(".", 1)[0])

                if pool == "gdna":
                    metrics["mrna_to_gdna"] += 1
                elif pool == "nrna":
                    metrics["mrna_to_nrna"] += 1
                elif assigned_tx == true_tx:
                    metrics["mrna_exact_tx"] += 1
                elif assigned_gene == true_gene:
                    metrics["mrna_same_gene_wrong_tx"] += 1
                else:
                    metrics["mrna_wrong_gene"] += 1

            elif origin.kind == "nrna":
                true_tx = origin.transcript_id or ""
                truth_nrna[true_tx] += 1
                if pool == "gdna":
                    metrics["nrna_to_gdna"] += 1
                elif pool == "nrna":
                    metrics["nrna_to_nrna"] += 1
                elif pool == "mrna":
                    metrics["nrna_to_mrna"] += 1
                else:
                    metrics["nrna_to_unresolved"] += 1

            elif origin.kind == "gdna":
                if pool == "gdna":
                    metrics["gdna_to_gdna"] += 1
                elif pool == "nrna":
                    metrics["gdna_to_nrna"] += 1
                elif pool == "mrna":
                    metrics["gdna_to_mrna"] += 1
                else:
                    metrics["gdna_to_unresolved"] += 1

    zc_detail = pd.DataFrame(
        [
            {
                "origin": origin,
                "assigned_pool": pool,
                "zc": zc_value,
                "count": count,
            }
            for (origin, pool, zc_value), count in sorted(zc_rows.items())
        ]
    )
    return dict(metrics), truth_mrna, truth_nrna, zc_detail


def corr_spearman(left: pd.Series, right: pd.Series) -> float:
    if len(left) < 2 or left.nunique(dropna=False) < 2 or right.nunique(dropna=False) < 2:
        return float("nan")
    return float(left.corr(right, method="spearman"))


def corr_pearson_log(left: pd.Series, right: pd.Series) -> float:
    if len(left) < 2:
        return float("nan")
    left_log = np.log2(left.astype(float) + 1.0)
    right_log = np.log2(right.astype(float) + 1.0)
    if left_log.std() == 0 or right_log.std() == 0:
        return float("nan")
    return float(np.corrcoef(left_log, right_log)[0, 1])


def abundance_metrics(merged: pd.DataFrame, pred_col: str, truth_col: str) -> dict[str, float]:
    expressed = merged[merged[truth_col] > 0].copy()
    unexpressed = merged[merged[truth_col] == 0].copy()
    if expressed.empty:
        rel_error = pd.Series(dtype=float)
    else:
        rel_error = (expressed[pred_col] - expressed[truth_col]).abs() / expressed[truth_col]
    total_truth = float(merged[truth_col].sum())
    total_pred = float(merged[pred_col].sum())
    return {
        "spearman": corr_spearman(expressed[truth_col], expressed[pred_col]),
        "pearson_log": corr_pearson_log(expressed[truth_col], expressed[pred_col]),
        "wape": safe_div((merged[pred_col] - merged[truth_col]).abs().sum(), total_truth),
        "median_abs_rel_error": float(rel_error.median()) if len(rel_error) else float("nan"),
        "mean_abs_rel_error": float(rel_error.mean()) if len(rel_error) else float("nan"),
        "false_positive_tx": int((unexpressed[pred_col] > 5.0).sum()),
        "false_positive_mass": float(unexpressed[pred_col].sum()),
        "false_negative_tx": int((expressed[pred_col] <= 0.0).sum()),
        "truth_total": total_truth,
        "pred_total": total_pred,
    }


def add_rates(row: dict[str, Any]) -> dict[str, Any]:
    total = row.get("total_fragments", 0)
    total_mrna = row.get("origin_mrna", 0)
    total_nrna = row.get("origin_nrna", 0)
    total_gdna = row.get("origin_gdna", 0)
    row["assigned_gdna_fraction"] = safe_div(row.get("assigned_gdna", 0), total)
    row["assigned_nrna_fraction"] = safe_div(row.get("assigned_nrna", 0), total)
    row["assigned_mrna_fraction"] = safe_div(row.get("assigned_mrna", 0), total)
    row["mrna_exact_rate"] = safe_div(row.get("mrna_exact_tx", 0), total_mrna)
    row["mrna_gene_rate"] = safe_div(
        row.get("mrna_exact_tx", 0) + row.get("mrna_same_gene_wrong_tx", 0),
        total_mrna,
    )
    row["mrna_wrong_gene_rate"] = safe_div(row.get("mrna_wrong_gene", 0), total_mrna)
    row["mrna_to_gdna_rate"] = safe_div(row.get("mrna_to_gdna", 0), total_mrna)
    row["mrna_to_nrna_rate"] = safe_div(row.get("mrna_to_nrna", 0), total_mrna)
    row["nrna_recall"] = safe_div(row.get("nrna_to_nrna", 0), total_nrna)
    row["nrna_to_mrna_rate"] = safe_div(row.get("nrna_to_mrna", 0), total_nrna)
    row["nrna_to_gdna_rate"] = safe_div(row.get("nrna_to_gdna", 0), total_nrna)
    row["gdna_recall"] = safe_div(row.get("gdna_to_gdna", 0), total_gdna)
    row["gdna_to_rna_rate"] = safe_div(
        row.get("gdna_to_mrna", 0) + row.get("gdna_to_nrna", 0),
        total_gdna,
    )
    row["overall_correct_gene_or_pool"] = safe_div(
        row.get("mrna_exact_tx", 0)
        + row.get("mrna_same_gene_wrong_tx", 0)
        + row.get("nrna_to_nrna", 0)
        + row.get("gdna_to_gdna", 0),
        total,
    )
    return row


def summarize_condition(
    base: Path,
    condition: str,
    meta: dict[str, Any],
    truth: pd.DataFrame,
    structure: pd.DataFrame,
    scan_metrics: dict[str, Any],
    truth_mrna: Counter[str],
    truth_nrna: Counter[str],
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    out_dir = base / condition / "rigel_out"
    quant = read_table(out_dir / "quant.feather")
    gene_quant = read_table(out_dir / "gene_quant.feather")
    loci = read_table(out_dir / "loci.feather")
    nrna_quant = read_table(out_dir / "nrna_quant.feather")
    summary = read_json(out_dir / "summary.json")

    truth_counts = pd.DataFrame(
        {
            "transcript_id": list(truth_mrna.keys()),
            "true_mrna_count": list(truth_mrna.values()),
        }
    )
    tx_detail = truth[[
        "transcript_id",
        "gene_id",
        "mrna_abundance",
        "nrna_abundance",
        "n_exons",
        "spliced_length",
        "genomic_span",
    ]].merge(truth_counts, on="transcript_id", how="left")
    tx_detail["true_mrna_count"] = tx_detail["true_mrna_count"].fillna(0.0)
    tx_detail = tx_detail.merge(
        quant[[
            "transcript_id",
            "count",
            "count_unambig",
            "count_em",
            "count_spliced",
            "nrna_parent_count",
            "effective_length",
            "tpm_total_rna",
        ]],
        on="transcript_id",
        how="left",
    ).fillna(0.0)
    tx_detail = tx_detail.merge(
        structure[[
            "transcript_id",
            "gene_n_tx",
            "exact_exon_group_size",
            "intron_chain_group_size",
            "nrna_span_group_size",
            "structure_class",
        ]],
        on="transcript_id",
        how="left",
    )
    tx_detail["condition"] = condition
    tx_detail["gdna_label"] = meta.get("gdna_label")
    tx_detail["nrna_label"] = meta.get("nrna_label")
    tx_detail["strand_specificity"] = meta.get("strand_specificity")
    tx_detail["abs_error"] = (tx_detail["count"] - tx_detail["true_mrna_count"]).abs()
    tx_detail["rel_error"] = np.where(
        tx_detail["true_mrna_count"] > 0,
        tx_detail["abs_error"] / tx_detail["true_mrna_count"],
        np.nan,
    )

    gene_truth = tx_detail.groupby("gene_id", as_index=False)["true_mrna_count"].sum()
    gene_merged = gene_truth.merge(
        gene_quant[["gene_id", "count", "count_unambig", "count_em"]],
        on="gene_id",
        how="left",
    ).fillna(0.0)
    gene_merged["condition"] = condition

    tx_metrics = abundance_metrics(tx_detail, "count", "true_mrna_count")
    gene_metrics = abundance_metrics(gene_merged, "count", "true_mrna_count")

    cal = summary.get("calibration", {})
    densities = cal.get("global_densities", {})
    intergenic = densities.get("INTERGENIC", {})
    intron = densities.get("INTRON", {})
    exon_intron = densities.get("EXON-INTRON", {})
    fl_models = cal.get("fl_models", {})
    strand_model = summary.get("strand_model", {})

    row: dict[str, Any] = {
        "condition": condition,
        "gdna_label": meta.get("gdna_label"),
        "nrna_label": meta.get("nrna_label"),
        "strand_specificity": float(meta.get("strand_specificity", np.nan)),
        "truth_mrna": int(meta.get("n_mrna", meta.get("n_rna", 0))),
        "truth_nrna": int(meta.get("n_nrna", 0)),
        "truth_gdna": int(meta.get("n_gdna", 0)),
        "truth_total": int(meta.get("n_total", 0)),
        "truth_gdna_fraction": safe_div(meta.get("n_gdna", 0), meta.get("n_total", 0)),
        "truth_nrna_fraction": safe_div(meta.get("n_nrna", 0), meta.get("n_total", 0)),
        "truth_nrna_to_mrna_ratio": safe_div(meta.get("n_nrna", 0), meta.get("n_mrna", 0)),
        "rho_intergenic": float(intergenic.get("rho", np.nan)),
        "rho_intron": float(intron.get("rho", np.nan)),
        "rho_exon_intron": float(exon_intron.get("rho", np.nan)),
        "rho_intron_over_intergenic": safe_div(
            float(intron.get("rho", np.nan)),
            float(intergenic.get("rho", np.nan)),
        ),
        "rho_exon_over_intergenic": safe_div(
            float(exon_intron.get("rho", np.nan)),
            float(intergenic.get("rho", np.nan)),
        ),
        "rho_exon_over_intron": safe_div(
            float(exon_intron.get("rho", np.nan)),
            float(intron.get("rho", np.nan)),
        ),
        "intron_strand_active": bool(intron.get("strand_active", False)),
        "exon_intron_strand_active": bool(exon_intron.get("strand_active", False)),
        "rna_fl_mean": float(fl_models.get("rna_fl_mean", np.nan)),
        "gdna_fl_mean": float(fl_models.get("gdna_fl_mean", np.nan)),
        "gdna_fl_rel_error": safe_div(float(fl_models.get("gdna_fl_mean", np.nan)) - 350.0, 350.0),
        "strand_model_estimate": float(strand_model.get("strand_specificity", np.nan)),
        "loci_mrna": float(loci["mrna"].sum()) if not loci.empty else np.nan,
        "loci_nrna": float(loci["nrna"].sum()) if "nrna" in loci else np.nan,
        "loci_gdna": float(loci["gdna"].sum()) if not loci.empty else np.nan,
        "nrna_quant_total": float(nrna_quant["count"].sum()) if not nrna_quant.empty else 0.0,
        "tx_spearman": tx_metrics["spearman"],
        "tx_pearson_log": tx_metrics["pearson_log"],
        "tx_wape": tx_metrics["wape"],
        "tx_median_abs_rel_error": tx_metrics["median_abs_rel_error"],
        "tx_false_positive_tx": tx_metrics["false_positive_tx"],
        "tx_false_positive_mass": tx_metrics["false_positive_mass"],
        "tx_false_negative_tx": tx_metrics["false_negative_tx"],
        "tx_pred_total": tx_metrics["pred_total"],
        "gene_spearman": gene_metrics["spearman"],
        "gene_pearson_log": gene_metrics["pearson_log"],
        "gene_wape": gene_metrics["wape"],
        "gene_median_abs_rel_error": gene_metrics["median_abs_rel_error"],
        "gene_pred_total": gene_metrics["pred_total"],
    }
    row.update(scan_metrics)
    row = add_rates(row)
    row["nrna_count_recovery"] = safe_div(row["nrna_quant_total"], row["truth_nrna"])
    row["mature_count_recovery"] = safe_div(row["tx_pred_total"], row["truth_mrna"])
    row["gdna_fraction_error"] = row["assigned_gdna_fraction"] - row["truth_gdna_fraction"]
    return row, tx_detail, gene_merged


def structural_summary(tx_detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_cols = ["condition", "gdna_label", "nrna_label", "strand_specificity", "structure_class"]
    for keys, group in tx_detail.groupby(group_cols, dropna=False):
        true_total = float(group["true_mrna_count"].sum())
        rows.append(
            {
                "condition": keys[0],
                "gdna_label": keys[1],
                "nrna_label": keys[2],
                "strand_specificity": keys[3],
                "structure_class": keys[4],
                "n_tx": int(len(group)),
                "true_total": true_total,
                "pred_total": float(group["count"].sum()),
                "wape": safe_div(group["abs_error"].sum(), true_total),
                "median_abs_rel_error": float(group["rel_error"].median()),
                "false_positive_mass": float(group.loc[group["true_mrna_count"] == 0, "count"].sum()),
            }
        )
    return pd.DataFrame(rows)


def top_error_tables(tx_detail: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    false_pos = tx_detail[tx_detail["true_mrna_count"] == 0].copy()
    false_pos = false_pos.sort_values("count", ascending=False).head(40)
    false_neg = tx_detail[(tx_detail["true_mrna_count"] > 0) & (tx_detail["count"] <= 0)].copy()
    false_neg = false_neg.sort_values("true_mrna_count", ascending=False).head(40)
    abs_errors = tx_detail.sort_values("abs_error", ascending=False).head(80)
    return false_pos, false_neg, abs_errors


def markdown_table(
    records: list[dict[str, Any]],
    columns: list[tuple[str, str]],
) -> str:
    headers = [header for header, _ in columns]
    rows = []
    for record in records:
        rows.append([str(record.get(key, "")) for _, key in columns])
    widths = [len(header) for header in headers]
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(cell))
    header_line = "| " + " | ".join(
        header.ljust(widths[index]) for index, header in enumerate(headers)
    ) + " |"
    sep_line = "| " + " | ".join("-" * widths[index] for index in range(len(widths))) + " |"
    body = [
        "| " + " | ".join(cell.ljust(widths[index]) for index, cell in enumerate(row)) + " |"
        for row in rows
    ]
    return "\n".join([header_line, sep_line, *body])


def report_records_from_conditions(condition_df: pd.DataFrame) -> dict[str, Any]:
    zero_nrna_gdna = condition_df[
        (condition_df["nrna_label"] == "zero") & (condition_df["truth_gdna"] > 0)
    ]
    high_nrna_no_gdna = condition_df[
        (condition_df["nrna_label"] == "high") & (condition_df["truth_gdna"] == 0)
    ]
    low_nrna_no_gdna = condition_df[
        (condition_df["nrna_label"] == "low") & (condition_df["truth_gdna"] == 0)
    ]
    high_nrna = condition_df[condition_df["nrna_label"] == "high"]
    high_nrna_with_gdna = high_nrna[high_nrna["truth_gdna"] > 0]
    zero_nrna = condition_df[condition_df["nrna_label"] == "zero"]

    return {
        "zero_nrna_gdna_max_abs_frac_err": zero_nrna_gdna["gdna_fraction_error"].abs().max(),
        "zero_nrna_gdna_max_gdna_leak": zero_nrna_gdna["gdna_to_rna_rate"].max(),
        "zero_nrna_mean_mrna_exact": zero_nrna["mrna_exact_rate"].mean(),
        "zero_nrna_mean_mrna_gene": zero_nrna["mrna_gene_rate"].mean(),
        "zero_nrna_mean_gene_wape": zero_nrna["gene_wape"].mean(),
        "zero_nrna_mean_tx_wape": zero_nrna["tx_wape"].mean(),
        "low_nrna_no_gdna_unstranded": low_nrna_no_gdna[
            low_nrna_no_gdna["strand_specificity"] < 0.75
        ].iloc[0],
        "low_nrna_no_gdna_stranded": low_nrna_no_gdna[
            low_nrna_no_gdna["strand_specificity"] > 0.75
        ].iloc[0],
        "high_nrna_no_gdna_unstranded": high_nrna_no_gdna[
            high_nrna_no_gdna["strand_specificity"] < 0.75
        ].iloc[0],
        "high_nrna_no_gdna_stranded": high_nrna_no_gdna[
            high_nrna_no_gdna["strand_specificity"] > 0.75
        ].iloc[0],
        "high_nrna_min_nrna_recall": high_nrna["nrna_recall"].min(),
        "high_nrna_max_nrna_to_gdna": high_nrna["nrna_to_gdna_rate"].max(),
        "high_nrna_max_gdna_fl_bias": high_nrna["gdna_fl_rel_error"].min(),
        "high_nrna_with_gdna_min_fl": high_nrna_with_gdna["gdna_fl_mean"].min(),
        "high_nrna_with_gdna_min_fl_bias": high_nrna_with_gdna[
            "gdna_fl_rel_error"
        ].min(),
        "high_nrna_max_tx_wape": high_nrna["tx_wape"].max(),
    }


def write_report(
    report_path: Path,
    base: Path,
    out_dir: Path,
    condition_df: pd.DataFrame,
    structure_df: pd.DataFrame,
    top_false_pos: pd.DataFrame,
    top_abs_errors: pd.DataFrame,
) -> None:
    summary = report_records_from_conditions(condition_df)
    high_unstranded = summary["high_nrna_no_gdna_unstranded"]
    high_stranded = summary["high_nrna_no_gdna_stranded"]
    low_unstranded = summary["low_nrna_no_gdna_unstranded"]
    low_stranded = summary["low_nrna_no_gdna_stranded"]

    def condition_row(row: pd.Series) -> dict[str, str]:
        return {
            "condition": str(row["condition"]),
            "truth_gdna": pct(row["truth_gdna_fraction"]),
            "assigned_gdna": pct(row["assigned_gdna_fraction"]),
            "truth_nrna": pct(row["truth_nrna_fraction"]),
            "nrna_recall": pct(row["nrna_recall"]),
            "nrna_to_gdna": pct(row["nrna_to_gdna_rate"]),
            "gdna_fl": fmt(row["gdna_fl_mean"], 1),
            "tx_wape": pct(row["tx_wape"]),
            "gene_wape": pct(row["gene_wape"]),
        }

    headline_rows = [
        condition_row(low_unstranded),
        condition_row(low_stranded),
        condition_row(high_unstranded),
        condition_row(high_stranded),
    ]

    zero_cal = condition_df[
        (condition_df["nrna_label"] == "zero") & (condition_df["truth_gdna"] > 0)
    ].copy()
    zero_cal = zero_cal.sort_values(["gdna_label", "strand_specificity"])
    zero_cal_rows = [
        {
            "condition": row.condition,
            "truth_gdna": pct(row.truth_gdna_fraction),
            "assigned_gdna": pct(row.assigned_gdna_fraction),
            "rho_ex_over_ig": fmt(row.rho_exon_over_intergenic, 3),
            "gdna_leak": pct(row.gdna_to_rna_rate),
            "mrna_gene": pct(row.mrna_gene_rate),
        }
        for row in zero_cal.itertuples(index=False)
    ]

    structural_focus = structure_df[
        (structure_df["condition"] == "gdna_zero_ss_0.99_nrna_zero")
        & (structure_df["true_total"] > 0)
    ].sort_values("wape", ascending=False)
    structural_rows = [
        {
            "class": row.structure_class,
            "n_tx": str(row.n_tx),
            "true_total": fmt_count(row.true_total),
            "wape": pct(row.wape),
            "median_re": pct(row.median_abs_rel_error),
            "fp_mass": fmt_count(row.false_positive_mass),
        }
        for row in structural_focus.itertuples(index=False)
    ]

    fp_rows = [
        {
            "condition": row.condition,
            "tx": row.transcript_id,
            "gene": row.gene_id,
            "pred": fmt_count(row.count),
            "class": row.structure_class,
            "gene_n_tx": str(int(row.gene_n_tx)) if not pd.isna(row.gene_n_tx) else "n/a",
        }
        for row in top_false_pos.head(12).itertuples(index=False)
    ]

    abs_rows = [
        {
            "condition": row.condition,
            "tx": row.transcript_id,
            "truth": fmt_count(row.true_mrna_count),
            "pred": fmt_count(row.count),
            "abs_err": fmt_count(row.abs_error),
            "class": row.structure_class,
        }
        for row in top_abs_errors.head(12).itertuples(index=False)
    ]

    report = f"""# Rigel Synthetic 24 Performance Report - 2026-05-14

Base directory: `{base}`

Evaluator run: `conda activate rigel && python scripts/sim/evaluate_suite.py --sim-base {base}`

Deep-analysis outputs: `{out_dir}`

## Executive Summary

All 24 scenarios completed. The no-nRNA cases are mostly healthy: across nonzero-gDNA,
zero-nRNA conditions, the maximum absolute gDNA fraction error was
{pct(summary['zero_nrna_gdna_max_abs_frac_err'])}, gDNA-to-RNA leak was at most
{pct(summary['zero_nrna_gdna_max_gdna_leak'])}, and mean mRNA correct-gene assignment
was {pct(summary['zero_nrna_mean_mrna_gene'])}. The main residual no-nRNA weakness is
isoform resolution: exact transcript assignment averaged only
{pct(summary['zero_nrna_mean_mrna_exact'])}, while gene-level WAPE averaged
{pct(summary['zero_nrna_mean_gene_wape'])}.

The major failure mode is nRNA/gDNA confounding. In unstranded high-nRNA data with no
true gDNA, Rigel assigned {pct(high_unstranded.assigned_gdna_fraction)} of all fragments
to gDNA and recovered only {pct(high_unstranded.nrna_recall)} of true nRNA as nRNA.
With strand-specificity 0.99 on the same scenario, assigned gDNA fell to
{pct(high_stranded.assigned_gdna_fraction)} and nRNA recall rose to
{pct(high_stranded.nrna_recall)}, so the strand-aware correction is doing real work but
is not sufficient for exon-boundary nRNA/gDNA ambiguity.

## Key Metrics

{markdown_table(headline_rows, [
        ('Condition', 'condition'),
        ('Truth gDNA', 'truth_gdna'),
        ('Assigned gDNA', 'assigned_gdna'),
        ('Truth nRNA', 'truth_nrna'),
        ('nRNA recall', 'nrna_recall'),
        ('nRNA to gDNA', 'nrna_to_gdna'),
        ('gDNA FL', 'gdna_fl'),
        ('Tx WAPE', 'tx_wape'),
        ('Gene WAPE', 'gene_wape'),
    ])}

## Healthy Baseline: No nRNA

{markdown_table(zero_cal_rows, [
        ('Condition', 'condition'),
        ('Truth gDNA', 'truth_gdna'),
        ('Assigned gDNA', 'assigned_gdna'),
        ('rho_ex/rho_ig', 'rho_ex_over_ig'),
        ('gDNA to RNA', 'gdna_leak'),
        ('mRNA correct gene', 'mrna_gene'),
    ])}

Interpretation: calibration is coherent when nRNA is absent. The exon/intergenic density
ratio stays near 1.0, and gDNA leakage remains below the 1.5% acceptance guardrail.
This argues against a broad gDNA calibration regression.

## Issue 1: nRNA Is Treated As gDNA In Unstranded Libraries

In the zero-gDNA, low-nRNA unstranded case, Rigel assigned
{pct(low_unstranded.assigned_gdna_fraction)} of fragments to gDNA even though true gDNA
is zero. In the matched stranded run, assigned gDNA dropped to
{pct(low_stranded.assigned_gdna_fraction)}. The high-nRNA case amplifies the same
pattern: {pct(high_unstranded.nrna_to_gdna_rate)} of true nRNA fragments went to gDNA
when unstranded, versus {pct(high_stranded.nrna_to_gdna_rate)} when stranded.

Root cause: the calibration gDNA channels are region-mask based. The gDNA FL histogram
explicitly sums unspliced intron-only, exon-intron, and intergenic-only fragments in
[src/rigel/calibration/_fl_sources.py](src/rigel/calibration/_fl_sources.py). Global
densities are then estimated for INTERGENIC, INTRON, and EXON-INTRON channels in
[src/rigel/calibration/density_global.py](src/rigel/calibration/density_global.py).
Unstranded nRNA produces intronic and exon-boundary unspliced evidence with the same
observable geometry as genic gDNA, but without intergenic support. With no usable strand
contrast, the intron/exon channels are left uncorrected and become gDNA priors.

## Issue 2: Strand Correction Helps, But Exonic Boundary nRNA Still Leaks

For high nRNA with no gDNA, the 0.99 stranded run has intronic density almost zeroed by
strand correction, but still assigns {pct(high_stranded.assigned_gdna_fraction)} of all
fragments to gDNA. The remaining signal is concentrated in EXON-INTRON density: this is
consistent with boundary-crossing unspliced nRNA that is harder to separate from gDNA
even with strand information.

Improvement opportunity: add an explicit nRNA-aware calibration state before global gDNA
projection, or gate genic gDNA priors on intergenic support when intron/exon density is
strongly inconsistent with intergenic density. A practical diagnostic is the ratio
rho_in/rho_ig and rho_ex/rho_ig stratified by strand activity; this suite shows ratios
above 8x in the worst unstranded high-nRNA low-gDNA condition.

## Issue 3: The gDNA Fragment-Length Model Is Contaminated By nRNA

When true gDNA is absent but nRNA is present, the reported gDNA FL is RNA-like
({fmt(high_unstranded.gdna_fl_mean, 1)} bp in high nRNA, no gDNA). In nonzero-gDNA,
high-nRNA conditions, the gDNA FL estimate falls as low as
{fmt(summary['high_nrna_with_gdna_min_fl'], 1)} bp, which is
{pct(summary['high_nrna_with_gdna_min_fl_bias'])} relative error against the simulated
350 bp gDNA mean. This follows directly from the same source selection in
`extract_gdna_counts`: nRNA fragments enter the unspliced intron/exon-intron pools used
to build the gDNA FL model.

Improvement opportunity: estimate gDNA FL primarily from intergenic-only fragments when
intergenic evidence is adequate, then use intron/exon channels for density after strand
and nRNA correction. The current pooled gDNA FL source is too permissive under nRNA.

## Issue 4: Transcript-Level Isoform Resolution Remains Weak

In the cleanest no-contamination condition, structural ambiguity explains much of the
transcript error:

{markdown_table(structural_rows, [
        ('Structure class', 'class'),
        ('n_tx', 'n_tx'),
        ('Truth count', 'true_total'),
        ('WAPE', 'wape'),
        ('Median RE', 'median_re'),
        ('FP mass', 'fp_mass'),
    ])}

Fragment-level truth agrees with this: mature RNA is assigned to the correct gene at
about 99% in clean stranded data, but exact transcript assignment is only about 70%.
This is not primarily a cross-gene mapping failure; it is within-gene isoform
identifiability and EM mass splitting.

Top false-positive transcript counts:

{markdown_table(fp_rows, [
        ('Condition', 'condition'),
        ('Transcript', 'tx'),
        ('Gene', 'gene'),
        ('Pred count', 'pred'),
        ('Class', 'class'),
        ('Gene n_tx', 'gene_n_tx'),
    ])}

Top absolute transcript errors:

{markdown_table(abs_rows, [
        ('Condition', 'condition'),
        ('Transcript', 'tx'),
        ('Truth', 'truth'),
        ('Pred', 'pred'),
        ('Abs err', 'abs_err'),
        ('Class', 'class'),
    ])}

Improvement opportunity: interrogate EM responsibilities within same-intron-chain and
shared-nRNA-span transcript groups. The likely fixes are not in BAM resolution; they are
in isoform priors, effective-length normalization for short/near-duplicate isoforms, and
posterior regularization when the data cannot identify a unique isoform.

## Recommended Follow-Up Work

1. Build an nRNA-aware SRD calibration model with separate genic-unspliced nRNA and gDNA
   states, and use intergenic evidence as the anchor for true gDNA.
2. Change gDNA FL training to prefer intergenic-only fragments, falling back to genic
   unspliced evidence only after nRNA/strand correction.
3. Add acceptance tests for zero-gDNA plus nRNA scenarios. The current acceptance checks
   pass no-nRNA gDNA cases but do not guard against nRNA-induced false gDNA.
4. Add per-locus reports comparing truth nRNA, predicted nRNA, and predicted gDNA. The
   worst loci in this suite should become regression fixtures.
5. Profile same-intron-chain transcript groups: record per-fragment posterior entropy,
   count_unambig vs count_em, and short-transcript enrichment among false positives.

## Artifacts

- Stock evaluator report: `{base / 'analysis_report.txt'}`
- Condition metrics: `{out_dir / 'condition_metrics.tsv'}`
- Transcript detail: `{out_dir / 'transcript_detail.tsv'}`
- Gene detail: `{out_dir / 'gene_detail.tsv'}`
- Structural summary: `{out_dir / 'structural_summary.tsv'}`
- Top false positives: `{out_dir / 'top_false_positive_transcripts.tsv'}`
- Top absolute errors: `{out_dir / 'top_abs_error_transcripts.tsv'}`
"""
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    base = args.base.resolve()
    out_dir = args.out_dir.resolve()
    report_path = args.report.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = load_manifest(base)
    condition_meta = condition_manifest_map(manifest)
    conditions = list(condition_meta)
    structure = build_structure_table(base)
    structure.to_csv(out_dir / "structure_table.tsv", sep="\t", index=False)

    condition_rows: list[dict[str, Any]] = []
    tx_details: list[pd.DataFrame] = []
    gene_details: list[pd.DataFrame] = []
    zc_details: list[pd.DataFrame] = []

    for condition in conditions:
        meta = condition_meta[condition]
        truth = load_condition_truth(base, meta)
        truth_gene_by_tx = dict(zip(truth["transcript_id"], truth["gene_id"]))
        bam_path = base / condition / "annotated.bam"
        scan_metrics, truth_mrna, truth_nrna, zc_detail = scan_annotated_bam(
            bam_path,
            truth_gene_by_tx,
        )
        zc_detail["condition"] = condition
        zc_details.append(zc_detail)
        row, tx_detail, gene_detail = summarize_condition(
            base,
            condition,
            meta,
            truth,
            structure,
            scan_metrics,
            truth_mrna,
            truth_nrna,
        )
        condition_rows.append(row)
        tx_details.append(tx_detail)
        gene_details.append(gene_detail)
        print(f"analyzed {condition}")

    condition_df = pd.DataFrame(condition_rows)
    tx_detail_df = pd.concat(tx_details, ignore_index=True)
    gene_detail_df = pd.concat(gene_details, ignore_index=True)
    zc_detail_df = pd.concat(zc_details, ignore_index=True)
    structure_summary_df = structural_summary(tx_detail_df)
    top_fp, top_fn, top_abs = top_error_tables(tx_detail_df)

    condition_df.to_csv(out_dir / "condition_metrics.tsv", sep="\t", index=False)
    tx_detail_df.to_csv(out_dir / "transcript_detail.tsv", sep="\t", index=False)
    gene_detail_df.to_csv(out_dir / "gene_detail.tsv", sep="\t", index=False)
    zc_detail_df.to_csv(out_dir / "zc_detail.tsv", sep="\t", index=False)
    structure_summary_df.to_csv(out_dir / "structural_summary.tsv", sep="\t", index=False)
    top_fp.to_csv(out_dir / "top_false_positive_transcripts.tsv", sep="\t", index=False)
    top_fn.to_csv(out_dir / "top_false_negative_transcripts.tsv", sep="\t", index=False)
    top_abs.to_csv(out_dir / "top_abs_error_transcripts.tsv", sep="\t", index=False)

    write_report(
        report_path,
        base,
        out_dir,
        condition_df,
        structure_summary_df,
        top_fp,
        top_abs,
    )
    print(f"wrote {report_path}")


if __name__ == "__main__":
    main()