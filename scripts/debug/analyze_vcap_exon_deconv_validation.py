#!/usr/bin/env python3
"""Compare VCaP RNA/gDNA confusion across recent Rigel calibration runs."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m")
DEFAULT_RESULTS = Path("results/vcap_exon_strand_deconv_validation_2026-05-20")
DEFAULT_REPORT = Path("docs/benchmarks/vcap_exon_strand_deconv_validation_2026-05-20.md")


@dataclass(frozen=True)
class RunSpec:
    label: str
    quant_dir: Path
    confusion_dir: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def fmt_int(value: float | int) -> str:
    return f"{int(round(float(value))):,}"


def fmt_float(value: float, digits: int = 4) -> str:
    return f"{float(value):.{digits}g}"


def fmt_pct(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    return f"{100.0 * float(value):.{digits}f}%"


def markdown_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    if not rows:
        return "_No rows._"
    header = "| " + " | ".join(label for label, _key in columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = ["| " + " | ".join(str(row.get(key, "")) for _label, key in columns) + " |" for row in rows]
    return "\n".join([header, sep, *body])


def read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def row_value(df: pd.DataFrame, true_source: str, predicted_pool: str) -> int:
    hit = df[(df["true_source"] == true_source) & (df["predicted_pool"] == predicted_pool)]
    if hit.empty:
        return 0
    return int(hit.iloc[0]["count"])


def row_rate(df: pd.DataFrame, true_source: str, predicted_pool: str) -> float:
    hit = df[(df["true_source"] == true_source) & (df["predicted_pool"] == predicted_pool)]
    if hit.empty:
        return 0.0
    return float(hit.iloc[0]["row_rate"])


def read_metrics(spec: RunSpec) -> dict[str, Any]:
    summary = read_json(spec.quant_dir / "summary.json")
    detail = pd.read_csv(spec.confusion_dir / "detailed_predicted_pool.tsv", sep="\t")
    matrix = pd.read_csv(spec.confusion_dir / "confusion_matrix.tsv", sep="\t")
    loci = pd.read_feather(spec.quant_dir / "loci.feather")

    gdna_total = sum(row_value(detail, "gdna", pool) for pool in ["mrna", "nrna", "gdna", "unresolved"])
    rna_total = sum(row_value(detail, "rna", pool) for pool in ["mrna", "nrna", "gdna", "unresolved"])

    quant = summary.get("quantification", {})
    cal = summary.get("calibration", {}) or {}
    regional = cal.get("regional_exposure", {}) or {}
    gdens = cal.get("global_densities", {}) or {}
    exon_composite = gdens.get("EXON-COMPOSITE", {}) or {}
    exon_contained = gdens.get("EXON-CONTAINED", {}) or {}

    return {
        "run": spec.label,
        "quant_dir": str(spec.quant_dir),
        "confusion_dir": str(spec.confusion_dir),
        "gdna_total_truth": gdna_total,
        "rna_total_truth": rna_total,
        "gdna_to_mrna": row_value(detail, "gdna", "mrna"),
        "gdna_to_nrna": row_value(detail, "gdna", "nrna"),
        "gdna_to_rna": row_value(matrix, "gdna", "rna"),
        "gdna_to_gdna": row_value(detail, "gdna", "gdna"),
        "gdna_unresolved": row_value(detail, "gdna", "unresolved"),
        "gdna_to_mrna_rate": row_rate(detail, "gdna", "mrna"),
        "gdna_to_nrna_rate": row_rate(detail, "gdna", "nrna"),
        "gdna_to_rna_rate": row_rate(matrix, "gdna", "rna"),
        "gdna_recall": row_rate(matrix, "gdna", "gdna"),
        "rna_to_gdna": row_value(detail, "rna", "gdna"),
        "rna_to_gdna_rate": row_rate(detail, "rna", "gdna"),
        "rna_recall": row_rate(matrix, "rna", "rna"),
        "em_mrna_total": float(quant.get("mrna_total", 0.0)),
        "em_nrna_total": float(quant.get("nrna_total", 0.0)),
        "em_gdna_total": float(quant.get("gdna_total", 0.0)),
        "mean_pi_gdna": float(cal.get("mean_pi_gdna", 0.0)),
        "gdna_eff_len_median": float(summary.get("gdna_eff_len", {}).get("value", {}).get("median", 0.0)),
        "gdna_weight_median": float(loci.get("gdna_em_exposure_weight", pd.Series([1.0])).median()),
        "gdna_weight_p05": float(loci.get("gdna_em_exposure_weight", pd.Series([1.0])).quantile(0.05)),
        "gdna_weight_p95": float(loci.get("gdna_em_exposure_weight", pd.Series([1.0])).quantile(0.95)),
        "regional_mode": regional.get("mode", "missing"),
        "regional_rho_global": float(regional.get("rho_global", 0.0)),
        "regional_rho_ref": float(regional.get("rho_ref", 0.0)),
        "regional_observed_log_spread": float(regional.get("observed_log_spread", 0.0)),
        "regional_null_log_spread": float(regional.get("null_log_spread", 0.0)),
        "regional_n_at_floor": int(regional.get("n_at_floor", 0)),
        "regional_negative_floors": int(regional.get("n_negative_rho_floored", 0)),
        "exon_intron_rho": float((gdens.get("EXON-INTRON", {}) or {}).get("rho", 0.0)),
        "exon_contained_rho": float(exon_contained.get("rho", 0.0)),
        "exon_contained_estimated": float(exon_contained.get("n_fragments_estimated", 0.0)),
        "exon_composite_rho": float(exon_composite.get("rho", 0.0)),
        "exon_composite_boundary_precision": float(exon_composite.get("boundary_precision", 0.0)),
        "exon_composite_contained_precision": float(exon_composite.get("contained_precision", 0.0)),
        "exon_composite_strand_power": float(exon_composite.get("strand_power", 0.0)),
    }


def build_specs(base_dir: Path, out_dir: Path) -> list[RunSpec]:
    return [
        RunSpec(
            "v4_3_with_mm",
            base_dir / "v4_3_with_mm",
            Path("results/vcap_rna20m_gdna20m_v4_3_confusion_2026-05-19"),
        ),
        RunSpec(
            "kappa_units_fix",
            base_dir / "kappa_units_fix",
            Path("results/vcap_rna20m_gdna20m_kappa_units_fix_confusion_2026-05-19"),
        ),
        RunSpec(
            "exon_strand_deconv_v1",
            base_dir / "exon_strand_deconv_v1",
            out_dir / "confusion",
        ),
    ]


def build_locus_delta(base_dir: Path, out_dir: Path) -> dict[str, pd.DataFrame]:
    old = pd.read_feather(base_dir / "kappa_units_fix" / "loci.feather")
    new = pd.read_feather(base_dir / "exon_strand_deconv_v1" / "loci.feather")
    keep = [
        "locus_id",
        "n_em_fragments",
        "mrna",
        "nrna",
        "gdna",
        "gdna_rate",
        "gdna_prior_count",
        "gdna_prior_count_em",
        "gdna_eff_len",
        "gdna_eff_len_unweighted",
        "gdna_eff_len_weight_ratio",
        "gdna_em_exposure_weight",
    ]
    merged = old[keep].merge(new[keep], on="locus_id", suffixes=("_kappa", "_new"))
    for col in ["mrna", "nrna", "gdna", "gdna_prior_count", "gdna_prior_count_em", "gdna_eff_len"]:
        merged[f"delta_{col}"] = merged[f"{col}_new"] - merged[f"{col}_kappa"]
    merged["delta_rna"] = merged["delta_mrna"] + merged["delta_nrna"]
    merged["gdna_weight_ratio_change"] = (
        merged["gdna_em_exposure_weight_new"] / merged["gdna_em_exposure_weight_kappa"].clip(lower=1e-12)
    )
    out = {
        "top_gdna_gain": merged.sort_values("delta_gdna", ascending=False).head(30),
        "top_gdna_loss": merged.sort_values("delta_gdna", ascending=True).head(30),
        "top_prior_gain": merged.sort_values("delta_gdna_prior_count", ascending=False).head(30),
        "all_locus_delta": merged,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in out.items():
        df.to_csv(out_dir / f"{name}.tsv", sep="\t", index=False)
    return out


def report_metric_rows(metrics: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows = []
    for m in metrics:
        rows.append(
            {
                "run": m["run"],
                "gdna_recall": fmt_pct(m["gdna_recall"]),
                "gdna_to_rna": f"{fmt_int(m['gdna_to_rna'])} ({fmt_pct(m['gdna_to_rna_rate'])})",
                "gdna_to_mrna": f"{fmt_int(m['gdna_to_mrna'])} ({fmt_pct(m['gdna_to_mrna_rate'])})",
                "gdna_to_nrna": f"{fmt_int(m['gdna_to_nrna'])} ({fmt_pct(m['gdna_to_nrna_rate'])})",
                "rna_to_gdna": f"{fmt_int(m['rna_to_gdna'])} ({fmt_pct(m['rna_to_gdna_rate'])})",
                "em_gdna": fmt_int(m["em_gdna_total"]),
                "em_mrna": fmt_int(m["em_mrna_total"]),
                "em_nrna": fmt_int(m["em_nrna_total"]),
            }
        )
    return rows


def report_delta_rows(metrics: list[dict[str, Any]], baseline: str) -> list[dict[str, str]]:
    by_run = {m["run"]: m for m in metrics}
    base = by_run[baseline]
    rows = []
    for m in metrics:
        if m["run"] == baseline:
            continue
        rows.append(
            {
                "comparison": f"{m['run']} - {baseline}",
                "delta_gdna_recall_pp": f"{100.0 * (m['gdna_recall'] - base['gdna_recall']):+.2f}",
                "delta_gdna_to_rna": f"{m['gdna_to_rna'] - base['gdna_to_rna']:+,}",
                "delta_gdna_to_mrna": f"{m['gdna_to_mrna'] - base['gdna_to_mrna']:+,}",
                "delta_gdna_to_nrna": f"{m['gdna_to_nrna'] - base['gdna_to_nrna']:+,}",
                "delta_rna_to_gdna": f"{m['rna_to_gdna'] - base['rna_to_gdna']:+,}",
                "delta_em_gdna": f"{m['em_gdna_total'] - base['em_gdna_total']:+,.0f}",
            }
        )
    return rows


def report_regional_rows(metrics: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows = []
    for m in metrics:
        rows.append(
            {
                "run": m["run"],
                "mode": m["regional_mode"],
                "rho_ref": fmt_float(m["regional_rho_ref"]),
                "spread": fmt_float(m["regional_observed_log_spread"]),
                "null": fmt_float(m["regional_null_log_spread"]),
                "floors": fmt_int(m["regional_n_at_floor"]),
                "neg_rate_floors": fmt_int(m["regional_negative_floors"]),
                "w_median": fmt_float(m["gdna_weight_median"]),
                "rho_ei": fmt_float(m["exon_intron_rho"]),
                "rho_ec": fmt_float(m["exon_contained_rho"]),
                "rho_comp": fmt_float(m["exon_composite_rho"]),
            }
        )
    return rows


def report_locus_rows(df: pd.DataFrame) -> list[dict[str, str]]:
    rows = []
    for row in df.head(12).to_dict("records"):
        rows.append(
            {
                "locus_id": int(row["locus_id"]),
                "n_em": fmt_int(row["n_em_fragments_new"]),
                "delta_gdna": f"{row['delta_gdna']:+,.0f}",
                "delta_mrna": f"{row['delta_mrna']:+,.0f}",
                "delta_nrna": f"{row['delta_nrna']:+,.0f}",
                "prior_kappa": fmt_float(row["gdna_prior_count_kappa"]),
                "prior_new": fmt_float(row["gdna_prior_count_new"]),
                "w_kappa": fmt_float(row["gdna_em_exposure_weight_kappa"]),
                "w_new": fmt_float(row["gdna_em_exposure_weight_new"]),
            }
        )
    return rows


def write_report(
    report: Path,
    metrics: list[dict[str, Any]],
    locus_tables: dict[str, pd.DataFrame],
    out_dir: Path,
) -> None:
    new = next(m for m in metrics if m["run"] == "exon_strand_deconv_v1")
    kappa = next(m for m in metrics if m["run"] == "kappa_units_fix")
    text = f"""# VCaP EXON Strand-Deconvolution Validation - 2026-05-20

Input BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam`

Fresh run: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1`

Truth source is flowcell-derived: RNA=`C6EL5ANXX`, gDNA=`H7MFFDSXY`. Counts are one primary read1 per fragment from Rigel annotated BAMs.

## Headline

Relative to `kappa_units_fix`, EXON strand deconvolution changed true-gDNA recall from {fmt_pct(kappa['gdna_recall'])} to {fmt_pct(new['gdna_recall'])}. True gDNA called RNA changed by {new['gdna_to_rna'] - kappa['gdna_to_rna']:+,} fragments; true RNA called gDNA changed by {new['rna_to_gdna'] - kappa['rna_to_gdna']:+,} fragments.

## Confusion Summary

{markdown_table(report_metric_rows(metrics), [
        ('Run', 'run'),
        ('gDNA recall', 'gdna_recall'),
        ('gDNA -> RNA', 'gdna_to_rna'),
        ('gDNA -> mRNA', 'gdna_to_mrna'),
        ('gDNA -> nRNA', 'gdna_to_nrna'),
        ('RNA -> gDNA', 'rna_to_gdna'),
        ('EM gDNA', 'em_gdna'),
        ('EM mRNA', 'em_mrna'),
        ('EM nRNA', 'em_nrna'),
    ])}

## Deltas

{markdown_table(report_delta_rows(metrics, 'kappa_units_fix'), [
        ('Comparison', 'comparison'),
        ('Δ gDNA recall pp', 'delta_gdna_recall_pp'),
        ('Δ gDNA -> RNA', 'delta_gdna_to_rna'),
        ('Δ gDNA -> mRNA', 'delta_gdna_to_mrna'),
        ('Δ gDNA -> nRNA', 'delta_gdna_to_nrna'),
        ('Δ RNA -> gDNA', 'delta_rna_to_gdna'),
        ('Δ EM gDNA', 'delta_em_gdna'),
    ])}

## Calibration And Exposure

{markdown_table(report_regional_rows(metrics), [
        ('Run', 'run'),
        ('Mode', 'mode'),
        ('rho_ref', 'rho_ref'),
        ('spread', 'spread'),
        ('null', 'null'),
        ('floor rows', 'floors'),
        ('neg floors', 'neg_rate_floors'),
        ('median weight', 'w_median'),
        ('rho EI', 'rho_ei'),
        ('rho EC', 'rho_ec'),
        ('rho composite', 'rho_comp'),
    ])}

## Largest Locus-Level gDNA Gains Versus kappa_units_fix

{markdown_table(report_locus_rows(locus_tables['top_gdna_gain']), [
        ('locus_id', 'locus_id'),
        ('n EM', 'n_em'),
        ('Δ gDNA', 'delta_gdna'),
        ('Δ mRNA', 'delta_mrna'),
        ('Δ nRNA', 'delta_nrna'),
        ('prior old', 'prior_kappa'),
        ('prior new', 'prior_new'),
        ('w old', 'w_kappa'),
        ('w new', 'w_new'),
    ])}

## Largest Locus-Level gDNA Losses Versus kappa_units_fix

{markdown_table(report_locus_rows(locus_tables['top_gdna_loss']), [
        ('locus_id', 'locus_id'),
        ('n EM', 'n_em'),
        ('Δ gDNA', 'delta_gdna'),
        ('Δ mRNA', 'delta_mrna'),
        ('Δ nRNA', 'delta_nrna'),
        ('prior old', 'prior_kappa'),
        ('prior new', 'prior_new'),
        ('w old', 'w_kappa'),
        ('w new', 'w_new'),
    ])}

## Artifacts

- Run metrics: `{out_dir / 'run_metrics.tsv'}`
- Detailed confusion comparison: `{out_dir / 'confusion_detail_compare.tsv'}`
- Locus deltas: `{out_dir / 'all_locus_delta.tsv'}`
- Top gDNA gains/losses: `{out_dir / 'top_gdna_gain.tsv'}`, `{out_dir / 'top_gdna_loss.tsv'}`
"""
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(text)


def main() -> int:
    args = parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    specs = build_specs(args.base_dir, out_dir)
    metrics = [read_metrics(spec) for spec in specs]
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(out_dir / "run_metrics.tsv", sep="\t", index=False)

    detail_frames = []
    for spec in specs:
        detail = pd.read_csv(spec.confusion_dir / "detailed_predicted_pool.tsv", sep="\t")
        detail.insert(0, "run", spec.label)
        detail_frames.append(detail)
    pd.concat(detail_frames, ignore_index=True).to_csv(
        out_dir / "confusion_detail_compare.tsv",
        sep="\t",
        index=False,
    )

    locus_tables = build_locus_delta(args.base_dir, out_dir)
    write_report(args.report, metrics, locus_tables, out_dir)
    print(f"wrote {args.report}")
    print(f"wrote {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())