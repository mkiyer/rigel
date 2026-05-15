#!/usr/bin/env python3
"""Compare synthetic 24 baseline Rigel outputs against a new output directory."""

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

DEFAULT_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24")
DEFAULT_OUT = Path("results/synthetic_24_efflen_2026-05-15")


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
    with path.open() as handle:
        return json.load(handle)


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".feather":
        return pd.read_feather(path)
    return pd.read_csv(path, sep="\t")


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


def scan_assignments(bam_path: Path) -> dict[str, float]:
    metrics: Counter[str] = Counter()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue
            origin = parse_origin(read.query_name)
            origin_kind = "mrna" if origin.kind == "rna" else origin.kind
            zf = int(get_tag(read, "ZF", 0))
            assigned_pool = assigned_pool_from_zf(zf)

            metrics["total_fragments"] += 1
            metrics[f"origin_{origin_kind}"] += 1
            metrics[f"assigned_{assigned_pool}"] += 1
            metrics[f"{origin_kind}_to_{assigned_pool}"] += 1

    row = dict(metrics)
    total = row.get("total_fragments", 0)
    origin_mrna = row.get("origin_mrna", 0)
    origin_nrna = row.get("origin_nrna", 0)
    origin_gdna = row.get("origin_gdna", 0)
    row["assigned_gdna_fraction"] = safe_div(row.get("assigned_gdna", 0), total)
    row["assigned_nrna_fraction"] = safe_div(row.get("assigned_nrna", 0), total)
    row["assigned_mrna_fraction"] = safe_div(row.get("assigned_mrna", 0), total)
    row["mrna_to_gdna_rate"] = safe_div(row.get("mrna_to_gdna", 0), origin_mrna)
    row["mrna_to_nrna_rate"] = safe_div(row.get("mrna_to_nrna", 0), origin_mrna)
    row["nrna_recall"] = safe_div(row.get("nrna_to_nrna", 0), origin_nrna)
    row["nrna_to_mrna_rate"] = safe_div(row.get("nrna_to_mrna", 0), origin_nrna)
    row["nrna_to_gdna_rate"] = safe_div(row.get("nrna_to_gdna", 0), origin_nrna)
    row["nrna_rna_compatible_rate"] = safe_div(
        row.get("nrna_to_nrna", 0) + row.get("nrna_to_mrna", 0),
        origin_nrna,
    )
    row["gdna_recall"] = safe_div(row.get("gdna_to_gdna", 0), origin_gdna)
    row["gdna_to_rna_rate"] = safe_div(
        row.get("gdna_to_mrna", 0) + row.get("gdna_to_nrna", 0),
        origin_gdna,
    )
    return row


def abundance_wape(truth: pd.DataFrame, quant: pd.DataFrame) -> float:
    merged = truth[["transcript_id", "mrna_abundance"]].merge(
        quant[["transcript_id", "count"]],
        on="transcript_id",
        how="left",
    )
    merged["count"] = merged["count"].fillna(0.0)
    return safe_div(
        (merged["count"] - merged["mrna_abundance"]).abs().sum(), merged["mrna_abundance"].sum()
    )


def summarize_output(
    *,
    base: Path,
    condition: str,
    meta: dict[str, Any],
    out_name: str,
    annotated_name: str,
    label: str,
) -> dict[str, Any]:
    out_dir = base / condition / out_name
    if not (out_dir / "quant.feather").exists():
        raise FileNotFoundError(f"Missing quant output: {out_dir}")

    truth = pd.read_csv(base / str(meta["truth_abundances"]), sep="\t")
    quant = read_table(out_dir / "quant.feather")
    loci = read_table(out_dir / "loci.feather")
    nrna = read_table(out_dir / "nrna_quant.feather")
    summary = read_json(out_dir / "summary.json")
    bam_metrics = scan_assignments(base / condition / annotated_name)

    pred_mrna = float(loci["mrna"].sum()) if "mrna" in loci else float("nan")
    pred_nrna = float(loci["nrna"].sum()) if "nrna" in loci else 0.0
    pred_gdna = float(loci["gdna"].sum()) if "gdna" in loci else float("nan")
    pred_total = pred_mrna + pred_nrna + pred_gdna

    gdna_eff = summary.get("gdna_eff_len", {})
    gdna_eff_value = gdna_eff.get("value", {}) if isinstance(gdna_eff, dict) else {}
    gdna_eff_per_bp = gdna_eff.get("per_bp", {}) if isinstance(gdna_eff, dict) else {}

    row: dict[str, Any] = {
        "label": label,
        "condition": condition,
        "gdna_label": meta.get("gdna_label"),
        "nrna_label": meta.get("nrna_label"),
        "strand_specificity": float(meta.get("strand_specificity", np.nan)),
        "truth_mrna": int(meta.get("n_mrna", 0)),
        "truth_nrna": int(meta.get("n_nrna", 0)),
        "truth_gdna": int(meta.get("n_gdna", 0)),
        "truth_total": int(meta.get("n_total", 0)),
        "pred_mrna": pred_mrna,
        "pred_nrna": pred_nrna,
        "pred_gdna": pred_gdna,
        "pred_total": pred_total,
        "pred_nrna_fraction": safe_div(pred_nrna, pred_total),
        "pred_gdna_fraction": safe_div(pred_gdna, pred_total),
        "nrna_quant_total": float(nrna["count"].sum()) if not nrna.empty else 0.0,
        "nrna_recovery_loci": safe_div(pred_nrna, meta.get("n_nrna", 0)),
        "gdna_recovery_loci": safe_div(pred_gdna, meta.get("n_gdna", 0)),
        "mrna_recovery_loci": safe_div(pred_mrna, meta.get("n_mrna", 0)),
        "tx_wape": abundance_wape(truth, quant),
        "gdna_eff_len_median": float(gdna_eff_value.get("median", np.nan)),
        "gdna_eff_len_per_bp_median": float(gdna_eff_per_bp.get("median", np.nan)),
    }
    row.update(bam_metrics)
    return row


def delta_table(metrics: pd.DataFrame, baseline_label: str, new_label: str) -> pd.DataFrame:
    baseline = metrics[metrics["label"] == baseline_label].set_index("condition")
    new = metrics[metrics["label"] == new_label].set_index("condition")
    rows = []
    delta_cols = [
        "pred_mrna",
        "pred_nrna",
        "pred_gdna",
        "pred_nrna_fraction",
        "pred_gdna_fraction",
        "nrna_recovery_loci",
        "gdna_recovery_loci",
        "mrna_recovery_loci",
        "tx_wape",
        "mrna_to_gdna_rate",
        "nrna_recall",
        "nrna_to_mrna_rate",
        "nrna_to_gdna_rate",
        "nrna_rna_compatible_rate",
        "gdna_recall",
        "gdna_to_rna_rate",
        "assigned_gdna_fraction",
        "assigned_nrna_fraction",
    ]
    for condition in new.index:
        left = baseline.loc[condition]
        right = new.loc[condition]
        row = {
            "condition": condition,
            "gdna_label": right["gdna_label"],
            "nrna_label": right["nrna_label"],
            "strand_specificity": right["strand_specificity"],
            "truth_mrna": right["truth_mrna"],
            "truth_nrna": right["truth_nrna"],
            "truth_gdna": right["truth_gdna"],
        }
        for col in delta_cols:
            row[f"baseline_{col}"] = left[col]
            row[f"new_{col}"] = right[col]
            row[f"delta_{col}"] = right[col] - left[col]
        rows.append(row)
    return pd.DataFrame(rows)


def markdown_table(records: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    headers = [header for header, _ in columns]
    rows = [[str(record.get(key, "")) for _, key in columns] for record in records]
    widths = [len(header) for header in headers]
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(cell))
    header = "| " + " | ".join(name.ljust(widths[i]) for i, name in enumerate(headers)) + " |"
    sep = "| " + " | ".join("-" * width for width in widths) + " |"
    body = [
        "| " + " | ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)) + " |"
        for row in rows
    ]
    return "\n".join([header, sep, *body])


def grouped_delta(delta: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    cols = [
        "delta_pred_nrna",
        "delta_pred_gdna",
        "delta_nrna_to_gdna_rate",
        "delta_nrna_rna_compatible_rate",
        "delta_gdna_to_rna_rate",
        "delta_mrna_to_gdna_rate",
        "delta_tx_wape",
    ]
    return delta.groupby(group_cols, as_index=False)[cols].mean(numeric_only=True)


def write_report(
    *,
    path: Path,
    base: Path,
    baseline_name: str,
    new_name: str,
    delta: pd.DataFrame,
) -> None:
    key_conditions = [
        "gdna_zero_ss_0.50_nrna_high",
        "gdna_zero_ss_0.99_nrna_high",
        "gdna_low_ss_0.50_nrna_high",
        "gdna_high_ss_0.50_nrna_high",
    ]
    key_rows = []
    for condition in key_conditions:
        if condition not in set(delta["condition"]):
            continue
        row = delta[delta["condition"] == condition].iloc[0]
        key_rows.append(
            {
                "Condition": condition,
                "nRNA pred old": fmt_count(row["baseline_pred_nrna"]),
                "nRNA pred new": fmt_count(row["new_pred_nrna"]),
                "gDNA pred old": fmt_count(row["baseline_pred_gdna"]),
                "gDNA pred new": fmt_count(row["new_pred_gdna"]),
                "nRNA->gDNA old": pct(row["baseline_nrna_to_gdna_rate"]),
                "nRNA->gDNA new": pct(row["new_nrna_to_gdna_rate"]),
                "gDNA->RNA old": pct(row["baseline_gdna_to_rna_rate"]),
                "gDNA->RNA new": pct(row["new_gdna_to_rna_rate"]),
            }
        )

    group = grouped_delta(delta, ["nrna_label", "strand_specificity"])
    group_rows = []
    for _, row in group.iterrows():
        group_rows.append(
            {
                "nRNA": row["nrna_label"],
                "SS": fmt(row["strand_specificity"], 2),
                "d nRNA pred": fmt_count(row["delta_pred_nrna"]),
                "d gDNA pred": fmt_count(row["delta_pred_gdna"]),
                "d nRNA->gDNA": pct(row["delta_nrna_to_gdna_rate"]),
                "d nRNA RNA-compat": pct(row["delta_nrna_rna_compatible_rate"]),
                "d gDNA->RNA": pct(row["delta_gdna_to_rna_rate"]),
                "d tx WAPE": fmt(row["delta_tx_wape"], 4),
            }
        )

    high_unstranded = delta[(delta["nrna_label"] == "high") & (delta["strand_specificity"] == 0.5)]
    zero_gdna_high_unstranded = delta[delta["condition"] == "gdna_zero_ss_0.50_nrna_high"]
    if not zero_gdna_high_unstranded.empty:
        z = zero_gdna_high_unstranded.iloc[0]
        headline = (
            f"For gdna_zero_ss_0.50_nrna_high, predicted nRNA changed from "
            f"{fmt_count(z['baseline_pred_nrna'])} to {fmt_count(z['new_pred_nrna'])}, "
            f"while predicted gDNA changed from {fmt_count(z['baseline_pred_gdna'])} "
            f"to {fmt_count(z['new_pred_gdna'])}. True nRNA is {fmt_count(z['truth_nrna'])} "
            "and true gDNA is 0."
        )
    else:
        headline = "Key zero-gDNA high-nRNA unstranded condition was not available."

    if not high_unstranded.empty:
        high_summary = (
            "Across high-nRNA unstranded conditions, mean delta nRNA->gDNA rate was "
            f"{pct(high_unstranded['delta_nrna_to_gdna_rate'].mean())}, and mean "
            f"delta gDNA->RNA rate was {pct(high_unstranded['delta_gdna_to_rna_rate'].mean())}."
        )
    else:
        high_summary = "High-nRNA unstranded aggregate was not available."

    report = f"""# Synthetic 24 gDNA Effective-Length Rerun

Baseline output: `{baseline_name}`
New output: `{new_name}`
Base directory: `{base}`

## Headline

{headline}

{high_summary}

## Key Conditions

{
        markdown_table(
            key_rows,
            [
                ("Condition", "Condition"),
                ("nRNA pred old", "nRNA pred old"),
                ("nRNA pred new", "nRNA pred new"),
                ("gDNA pred old", "gDNA pred old"),
                ("gDNA pred new", "gDNA pred new"),
                ("nRNA->gDNA old", "nRNA->gDNA old"),
                ("nRNA->gDNA new", "nRNA->gDNA new"),
                ("gDNA->RNA old", "gDNA->RNA old"),
                ("gDNA->RNA new", "gDNA->RNA new"),
            ],
        )
    }

## Mean Deltas By nRNA Level And Strand Specificity

Positive deltas mean the new effective-length implementation increased that metric.

{
        markdown_table(
            group_rows,
            [
                ("nRNA", "nRNA"),
                ("SS", "SS"),
                ("d nRNA pred", "d nRNA pred"),
                ("d gDNA pred", "d gDNA pred"),
                ("d nRNA->gDNA", "d nRNA->gDNA"),
                ("d nRNA RNA-compat", "d nRNA RNA-compat"),
                ("d gDNA->RNA", "d gDNA->RNA"),
                ("d tx WAPE", "d tx WAPE"),
            ],
        )
    }

## Artifacts

- `condition_metrics.tsv`: stacked baseline/new metrics.
- `condition_delta.tsv`: one row per condition with old/new/delta columns.
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(report)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--baseline-name", default="rigel_out")
    parser.add_argument("--new-name", default="rigel_efflen_out")
    parser.add_argument("--baseline-annotated", default="annotated.bam")
    parser.add_argument("--new-annotated", default="annotated_efflen.bam")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    base = args.base.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(base)
    condition_meta = condition_manifest_map(manifest)

    rows = []
    for condition, meta in condition_meta.items():
        print(f"baseline {condition}")
        rows.append(
            summarize_output(
                base=base,
                condition=condition,
                meta=meta,
                out_name=args.baseline_name,
                annotated_name=args.baseline_annotated,
                label="baseline",
            )
        )
        print(f"new      {condition}")
        rows.append(
            summarize_output(
                base=base,
                condition=condition,
                meta=meta,
                out_name=args.new_name,
                annotated_name=args.new_annotated,
                label="new",
            )
        )

    metrics = pd.DataFrame(rows)
    delta = delta_table(metrics, "baseline", "new")
    metrics.to_csv(out_dir / "condition_metrics.tsv", sep="\t", index=False)
    delta.to_csv(out_dir / "condition_delta.tsv", sep="\t", index=False)
    write_report(
        path=out_dir / "summary.md",
        base=base,
        baseline_name=args.baseline_name,
        new_name=args.new_name,
        delta=delta,
    )
    print(f"wrote {out_dir / 'condition_metrics.tsv'}")
    print(f"wrote {out_dir / 'condition_delta.tsv'}")
    print(f"wrote {out_dir / 'summary.md'}")


if __name__ == "__main__":
    main()
