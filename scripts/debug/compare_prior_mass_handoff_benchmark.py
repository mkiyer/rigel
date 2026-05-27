"""Compare synthetic capture benchmark metrics before/after prior-mass handoff."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


KEY_COLUMNS = [
    "condition",
    "true_gdna_fraction",
    "est_gdna_fraction",
    "gdna_delta",
    "gdna_to_rna_rate",
    "mrna_to_gdna_rate",
    "est_gdna_total",
    "est_rna",
    "prior_global_gdna",
    "prior_global_rna",
    "mRNA_count_mard",
    "mRNA_count_spearman",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sim-base",
        type=Path,
        default=Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb"),
    )
    return parser.parse_args()


def fmt_float(value: float, digits: int = 4) -> str:
    if not np.isfinite(value):
        return "n/a"
    return f"{value:.{digits}f}"


def fmt_pct(value: float, digits: int = 2) -> str:
    if not np.isfinite(value):
        return "n/a"
    return f"{100.0 * value:.{digits}f}%"


def fmt_count(value: float) -> str:
    if not np.isfinite(value):
        return "n/a"
    return f"{value:,.0f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    widths = [len(header) for header in headers]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))

    def line(values: list[str]) -> str:
        return "| " + " | ".join(value.ljust(widths[idx]) for idx, value in enumerate(values)) + " |"

    divider = "| " + " | ".join("-" * width for width in widths) + " |"
    return "\n".join([line(headers), divider, *[line(row) for row in rows]])


def load_metrics(path: Path, label: str) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    missing = [column for column in KEY_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"{label} metrics missing columns: {missing}")
    return table[KEY_COLUMNS].copy()


def build_comparison(before: pd.DataFrame, after: pd.DataFrame) -> pd.DataFrame:
    merged = before.merge(after, on="condition", suffixes=("_before", "_after"))
    merged["abs_gdna_delta_before"] = merged["gdna_delta_before"].abs()
    merged["abs_gdna_delta_after"] = merged["gdna_delta_after"].abs()
    merged["abs_gdna_delta_change"] = (
        merged["abs_gdna_delta_after"] - merged["abs_gdna_delta_before"]
    )
    merged["gdna_fraction_abs_error_before"] = (
        merged["est_gdna_fraction_before"] - merged["true_gdna_fraction_before"]
    ).abs()
    merged["gdna_fraction_abs_error_after"] = (
        merged["est_gdna_fraction_after"] - merged["true_gdna_fraction_after"]
    ).abs()
    merged["gdna_fraction_abs_error_change"] = (
        merged["gdna_fraction_abs_error_after"] - merged["gdna_fraction_abs_error_before"]
    )
    merged["gdna_to_rna_rate_change"] = (
        merged["gdna_to_rna_rate_after"] - merged["gdna_to_rna_rate_before"]
    )
    merged["mrna_to_gdna_rate_change"] = (
        merged["mrna_to_gdna_rate_after"] - merged["mrna_to_gdna_rate_before"]
    )
    merged["prior_global_gdna_change"] = (
        merged["prior_global_gdna_after"] - merged["prior_global_gdna_before"]
    )
    return merged.sort_values("condition")


def write_report(comparison: pd.DataFrame, out_dir: Path) -> None:
    overview_rows: list[list[str]] = []
    for _, row in comparison.iterrows():
        overview_rows.append(
            [
                str(row["condition"]),
                fmt_pct(row["true_gdna_fraction_after"]),
                fmt_pct(row["est_gdna_fraction_before"]),
                fmt_pct(row["est_gdna_fraction_after"]),
                fmt_count(row["gdna_delta_before"]),
                fmt_count(row["gdna_delta_after"]),
                fmt_count(row["abs_gdna_delta_change"]),
                fmt_pct(row["gdna_to_rna_rate_before"]),
                fmt_pct(row["gdna_to_rna_rate_after"]),
                fmt_pct(row["mrna_to_gdna_rate_before"]),
                fmt_pct(row["mrna_to_gdna_rate_after"]),
            ]
        )

    prior_rows: list[list[str]] = []
    for _, row in comparison.iterrows():
        prior_rows.append(
            [
                str(row["condition"]),
                fmt_count(row["prior_global_gdna_before"]),
                fmt_count(row["prior_global_gdna_after"]),
                fmt_count(row["prior_global_gdna_change"]),
                fmt_count(row["prior_global_rna_before"]),
                fmt_count(row["prior_global_rna_after"]),
            ]
        )

    abundance_rows: list[list[str]] = []
    for _, row in comparison.iterrows():
        abundance_rows.append(
            [
                str(row["condition"]),
                fmt_pct(row["mRNA_count_mard_before"]),
                fmt_pct(row["mRNA_count_mard_after"]),
                fmt_float(row["mRNA_count_spearman_before"], 4),
                fmt_float(row["mRNA_count_spearman_after"], 4),
            ]
        )

    lines = ["# Prior Mass Handoff Benchmark", ""]
    lines.append("## Pool And Assignment Metrics")
    lines.append("")
    lines.append(
        markdown_table(
            [
                "Condition",
                "True gDNA",
                "Before Rigel gDNA",
                "After Rigel gDNA",
                "Before gDNA delta",
                "After gDNA delta",
                "Abs delta change",
                "Before gDNA->RNA",
                "After gDNA->RNA",
                "Before mRNA->gDNA",
                "After mRNA->gDNA",
            ],
            overview_rows,
        )
    )
    lines.append("")
    lines.append("Negative `Abs delta change` means improved absolute pool-level gDNA error.")
    lines.append("")
    lines.append("## Adaptive Prior Global Counts")
    lines.append("")
    lines.append(
        markdown_table(
            [
                "Condition",
                "Before prior gDNA",
                "After prior gDNA",
                "Prior gDNA change",
                "Before prior RNA",
                "After prior RNA",
            ],
            prior_rows,
        )
    )
    lines.append("")
    lines.append("## Transcript Count Metrics")
    lines.append("")
    lines.append(
        markdown_table(
            [
                "Condition",
                "Before mRNA MARD",
                "After mRNA MARD",
                "Before Spearman",
                "After Spearman",
            ],
            abundance_rows,
        )
    )
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    improved = int((comparison["abs_gdna_delta_change"] < 0.0).sum())
    worsened = int((comparison["abs_gdna_delta_change"] > 0.0).sum())
    same = int((comparison["abs_gdna_delta_change"] == 0.0).sum())
    lines.append(f"Pool-level absolute gDNA error improved in {improved}, worsened in {worsened}, unchanged in {same} scenarios.")
    lines.append(
        "The main intended rescue is the high-gDNA SS 0.99 capture-on scenario; "
        "unstranded capture-on remains a calibration-stage failure rather than a prior-handoff failure."
    )
    (out_dir / "prior_mass_handoff_comparison.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    sim_base = args.sim_base.resolve()
    out_dir = sim_base / "benchmarks" / "prior_mass_handoff_after"
    out_dir.mkdir(parents=True, exist_ok=True)
    before_path = sim_base / "benchmarks" / "prior_mass_handoff_before" / "condition_metrics.tsv"
    after_path = sim_base / "condition_metrics.tsv"
    before = load_metrics(before_path, "before")
    after = load_metrics(after_path, "after")
    comparison = build_comparison(before, after)
    comparison.to_csv(out_dir / "prior_mass_handoff_comparison.tsv", sep="\t", index=False)
    write_report(comparison, out_dir)
    print(f"[compare] Wrote {out_dir}")


if __name__ == "__main__":
    main()