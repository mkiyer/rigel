#!/usr/bin/env python3
"""Compare VCaP confusion for global high-tail exposure reference quantiles."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m")
DEFAULT_OUT = Path("results/vcap_global_high_tail_exposure_2026-05-21")
DEFAULT_REPORT = Path("docs/benchmarks/vcap_global_high_tail_exposure_validation_2026-05-21.md")
BASELINE_CONFUSION = Path("results/vcap_exon_strand_deconv_validation_2026-05-20/confusion")
Q999_CONFUSION = Path("results/vcap_global_q999_confusion_2026-05-21")
Q9995_CONFUSION = Path("results/vcap_global_q9995_confusion_2026-05-21")
DETAIL_POOLS = ("mrna", "nrna", "gdna", "unresolved")
BINARY_POOLS = ("rna", "gdna", "unresolved")


@dataclass(frozen=True)
class RunSpec:
    label: str
    quantile_label: str
    quant_dir: Path
    confusion_dir: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def fmt_int(value: float | int) -> str:
    return f"{int(round(float(value))):,}"


def fmt_float(value: float | int, digits: int = 4) -> str:
    return f"{float(value):.{digits}g}"


def fmt_pct(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    return f"{100.0 * float(value):.{digits}f}%"


def fmt_delta_int(value: float | int) -> str:
    return f"{int(round(float(value))):+,}"


def markdown_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    if not rows:
        return "_No rows._"
    header = "| " + " | ".join(label for label, _key in columns) + " |"
    sep = "| " + " | ".join("---" for _label, _key in columns) + " |"
    body = [
        "| " + " | ".join(str(row.get(key, "")) for _label, key in columns) + " |"
        for row in rows
    ]
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


def build_specs(base_dir: Path) -> list[RunSpec]:
    return [
        RunSpec("q95_baseline", "Q95", base_dir / "exon_strand_deconv_v1", BASELINE_CONFUSION),
        RunSpec("global_q999", "Q99.9", base_dir / "global_q999_v1", Q999_CONFUSION),
        RunSpec("global_q9995", "Q99.95", base_dir / "global_q9995_v1", Q9995_CONFUSION),
    ]


def read_metrics(spec: RunSpec) -> dict[str, Any]:
    summary = read_json(spec.quant_dir / "summary.json")
    detail = pd.read_csv(spec.confusion_dir / "detailed_predicted_pool.tsv", sep="\t")
    binary = pd.read_csv(spec.confusion_dir / "confusion_matrix.tsv", sep="\t")
    loci = pd.read_feather(spec.quant_dir / "loci.feather")
    locus3 = loci.loc[loci["locus_id"] == 3].iloc[0]

    gdna_total = sum(row_value(detail, "gdna", pool) for pool in DETAIL_POOLS)
    rna_total = sum(row_value(detail, "rna", pool) for pool in DETAIL_POOLS)
    cal = summary.get("calibration", {}) or {}
    quant = summary.get("quantification", {}) or {}
    regional = cal.get("regional_exposure", {}) or {}

    return {
        "run": spec.label,
        "quantile": spec.quantile_label,
        "quant_dir": str(spec.quant_dir),
        "confusion_dir": str(spec.confusion_dir),
        "reference_quantile": float(regional.get("reference_quantile", 0.95)),
        "rho_ref": float(regional.get("rho_ref", 0.0)),
        "rho_global": float(regional.get("rho_global", 0.0)),
        "n_at_floor": int(regional.get("n_at_floor", 0)),
        "gdna_truth_total": gdna_total,
        "rna_truth_total": rna_total,
        "gdna_to_mrna": row_value(detail, "gdna", "mrna"),
        "gdna_to_nrna": row_value(detail, "gdna", "nrna"),
        "gdna_to_gdna": row_value(detail, "gdna", "gdna"),
        "gdna_unresolved": row_value(detail, "gdna", "unresolved"),
        "gdna_to_rna": row_value(binary, "gdna", "rna"),
        "gdna_to_mrna_rate": row_rate(detail, "gdna", "mrna"),
        "gdna_to_nrna_rate": row_rate(detail, "gdna", "nrna"),
        "gdna_to_rna_rate": row_rate(binary, "gdna", "rna"),
        "gdna_recall": row_rate(binary, "gdna", "gdna"),
        "rna_to_gdna": row_value(detail, "rna", "gdna"),
        "rna_to_gdna_rate": row_rate(detail, "rna", "gdna"),
        "rna_recall": row_rate(binary, "rna", "rna"),
        "em_mrna_total": float(quant.get("mrna_total", 0.0)),
        "em_nrna_total": float(quant.get("nrna_total", 0.0)),
        "em_gdna_total": float(quant.get("gdna_total", 0.0)),
        "locus3_gdna": float(locus3["gdna"]),
        "locus3_mrna": float(locus3["mrna"]),
        "locus3_nrna": float(locus3["nrna"]),
        "locus3_gdna_eff_len": float(locus3["gdna_eff_len"]),
        "locus3_weight": float(locus3["gdna_em_exposure_weight"]),
        "locus3_prior_em": float(locus3["gdna_prior_count_em"]),
    }


def metric_rows(metrics: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows = []
    for item in metrics:
        rows.append(
            {
                "run": item["quantile"],
                "rho_ref": fmt_float(item["rho_ref"]),
                "gdna_recall": fmt_pct(item["gdna_recall"]),
                "gdna_to_rna": f"{fmt_int(item['gdna_to_rna'])} ({fmt_pct(item['gdna_to_rna_rate'])})",
                "gdna_to_mrna": f"{fmt_int(item['gdna_to_mrna'])} ({fmt_pct(item['gdna_to_mrna_rate'])})",
                "gdna_to_nrna": f"{fmt_int(item['gdna_to_nrna'])} ({fmt_pct(item['gdna_to_nrna_rate'])})",
                "rna_to_gdna": f"{fmt_int(item['rna_to_gdna'])} ({fmt_pct(item['rna_to_gdna_rate'])})",
                "em_gdna": fmt_int(item["em_gdna_total"]),
                "em_mrna": fmt_int(item["em_mrna_total"]),
                "em_nrna": fmt_int(item["em_nrna_total"]),
            }
        )
    return rows


def delta_rows(metrics: list[dict[str, Any]]) -> list[dict[str, str]]:
    base = next(item for item in metrics if item["run"] == "q95_baseline")
    rows = []
    for item in metrics:
        if item["run"] == "q95_baseline":
            continue
        leakage_delta = item["gdna_to_rna"] - base["gdna_to_rna"]
        leakage_reduction = -leakage_delta / base["gdna_to_rna"] if base["gdna_to_rna"] else 0.0
        rows.append(
            {
                "comparison": f"{item['quantile']} - Q95",
                "delta_gdna_recall_pp": f"{100.0 * (item['gdna_recall'] - base['gdna_recall']):+.2f}",
                "delta_gdna_to_rna": fmt_delta_int(leakage_delta),
                "leakage_reduction": fmt_pct(leakage_reduction),
                "delta_gdna_to_mrna": fmt_delta_int(item["gdna_to_mrna"] - base["gdna_to_mrna"]),
                "delta_gdna_to_nrna": fmt_delta_int(item["gdna_to_nrna"] - base["gdna_to_nrna"]),
                "delta_rna_to_gdna": fmt_delta_int(item["rna_to_gdna"] - base["rna_to_gdna"]),
                "delta_em_gdna": fmt_delta_int(item["em_gdna_total"] - base["em_gdna_total"]),
            }
        )
    return rows


def locus3_rows(metrics: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows = []
    for item in metrics:
        rows.append(
            {
                "run": item["quantile"],
                "weight": fmt_float(item["locus3_weight"], 6),
                "gdna_eff_len": fmt_int(item["locus3_gdna_eff_len"]),
                "prior_em": fmt_int(item["locus3_prior_em"]),
                "gdna": fmt_int(item["locus3_gdna"]),
                "mrna": fmt_int(item["locus3_mrna"]),
                "nrna": fmt_int(item["locus3_nrna"]),
            }
        )
    return rows


def interpretation(metrics: list[dict[str, Any]]) -> str:
    base = next(item for item in metrics if item["run"] == "q95_baseline")
    high_tail = [item for item in metrics if item["run"] != "q95_baseline"]
    lines = []
    for item in high_tail:
        total_delta = item["gdna_to_rna"] - base["gdna_to_rna"]
        mrna_delta = item["gdna_to_mrna"] - base["gdna_to_mrna"]
        nrna_delta = item["gdna_to_nrna"] - base["gdna_to_nrna"]
        rna_to_gdna_delta = item["rna_to_gdna"] - base["rna_to_gdna"]
        if total_delta > 0:
            verdict = "increased"
        elif total_delta < 0:
            verdict = "decreased"
        else:
            verdict = "did not change"
        lines.append(
            f"- {item['quantile']} {verdict} total true-gDNA -> RNA calls by "
            f"{fmt_delta_int(total_delta)} fragments versus Q95. The split was "
            f"{fmt_delta_int(mrna_delta)} mRNA and {fmt_delta_int(nrna_delta)} nRNA, "
            f"with RNA -> gDNA changing by {fmt_delta_int(rna_to_gdna_delta)}."
        )
    return "\n".join(lines)


def write_report(report: Path, metrics: list[dict[str, Any]], out_dir: Path) -> None:
    base = next(item for item in metrics if item["run"] == "q95_baseline")
    best = min((item for item in metrics if item["run"] != "q95_baseline"), key=lambda x: x["gdna_to_rna"])
    best_delta = best["gdna_to_rna"] - base["gdna_to_rna"]
    best_reduction = -best_delta / base["gdna_to_rna"] if base["gdna_to_rna"] else 0.0
    if best_delta > 0:
        headline = "increases gDNA-to-RNA leakage"
        headline_change = f"a {fmt_pct(best_delta / base['gdna_to_rna'])} increase relative to Q95"
    elif best_reduction >= 0.5:
        headline = "eliminates the majority of gDNA-to-RNA leakage"
        headline_change = f"a {fmt_pct(best_reduction)} reduction relative to Q95"
    elif best_delta < 0:
        headline = "reduces gDNA-to-RNA leakage, but not by a majority"
        headline_change = f"a {fmt_pct(best_reduction)} reduction relative to Q95"
    else:
        headline = "leaves gDNA-to-RNA leakage unchanged"
        headline_change = "no change relative to Q95"
    text = f"""# VCaP Global High-Tail Exposure Validation - 2026-05-21

Input run family: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m`

Compared runs:

- Q95 baseline: `{DEFAULT_BASE / 'exon_strand_deconv_v1'}`
- Q99.9: `{DEFAULT_BASE / 'global_q999_v1'}`
- Q99.95: `{DEFAULT_BASE / 'global_q9995_v1'}`

Truth source is flowcell-derived: RNA=`C6EL5ANXX`, gDNA=`H7MFFDSXY`. Counts are one primary read1 per fragment from Rigel annotated BAMs.

## Headline

The best high-tail run, `{best['quantile']}`, {headline}: true-gDNA fragments called RNA changed from {fmt_int(base['gdna_to_rna'])} to {fmt_int(best['gdna_to_rna'])}, {headline_change}.

## Confusion Summary

{markdown_table(metric_rows(metrics), [
        ('Run', 'run'),
        ('rho_ref', 'rho_ref'),
        ('gDNA recall', 'gdna_recall'),
        ('gDNA -> RNA', 'gdna_to_rna'),
        ('gDNA -> mRNA', 'gdna_to_mrna'),
        ('gDNA -> nRNA', 'gdna_to_nrna'),
        ('RNA -> gDNA', 'rna_to_gdna'),
        ('EM gDNA', 'em_gdna'),
        ('EM mRNA', 'em_mrna'),
        ('EM nRNA', 'em_nrna'),
    ])}

## Deltas Versus Q95 Baseline

{markdown_table(delta_rows(metrics), [
        ('Comparison', 'comparison'),
        ('Delta gDNA recall pp', 'delta_gdna_recall_pp'),
        ('Delta gDNA -> RNA', 'delta_gdna_to_rna'),
        ('Leakage reduction', 'leakage_reduction'),
        ('Delta gDNA -> mRNA', 'delta_gdna_to_mrna'),
        ('Delta gDNA -> nRNA', 'delta_gdna_to_nrna'),
        ('Delta RNA -> gDNA', 'delta_rna_to_gdna'),
        ('Delta EM gDNA', 'delta_em_gdna'),
    ])}

## Locus 3 / FLG2 Mega-Locus

{markdown_table(locus3_rows(metrics), [
        ('Run', 'run'),
        ('Exposure weight', 'weight'),
        ('L_gDNA', 'gdna_eff_len'),
        ('Prior EM', 'prior_em'),
        ('EM gDNA', 'gdna'),
        ('EM mRNA', 'mrna'),
        ('EM nRNA', 'nrna'),
    ])}

## Interpretation

{interpretation(metrics)}

The high-tail global reference quantile is therefore not sufficient as a standalone fix for the dominant VCaP leakage mode. The key diagnostic is the mRNA/nRNA split: if mRNA leakage falls while nRNA leakage rises, the higher reference scale is mainly changing the RNA sub-state balance rather than removing the false-RNA assignment pressure.

## Artifacts

- Metrics table: `{out_dir / 'run_metrics.tsv'}`
- Confusion detail comparison: `{out_dir / 'confusion_detail_compare.tsv'}`
"""
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(text)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    specs = build_specs(args.base_dir)
    metrics = [read_metrics(spec) for spec in specs]
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(args.out_dir / "run_metrics.tsv", sep="\t", index=False)

    details = []
    for spec in specs:
        detail = pd.read_csv(spec.confusion_dir / "detailed_predicted_pool.tsv", sep="\t")
        detail.insert(0, "run", spec.label)
        detail.insert(1, "quantile", spec.quantile_label)
        details.append(detail)
    pd.concat(details, ignore_index=True).to_csv(
        args.out_dir / "confusion_detail_compare.tsv",
        sep="\t",
        index=False,
    )
    write_report(args.report, metrics, args.out_dir)
    print(f"wrote {args.report}")
    print(f"wrote {args.out_dir / 'run_metrics.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
