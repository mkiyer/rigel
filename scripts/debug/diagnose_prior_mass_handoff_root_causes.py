"""Diagnose prior-mass handoff behavior across synthetic capture benchmark snapshots."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


CONDITIONS = (
    "gdna_none_ss_0.99_nrna_none_capture_off",
    "gdna_none_ss_0.99_nrna_none_capture_on",
    "gdna_none_ss_0.50_nrna_none_capture_off",
    "gdna_none_ss_0.50_nrna_none_capture_on",
    "gdna_high_ss_0.99_nrna_none_capture_off",
    "gdna_high_ss_0.99_nrna_none_capture_on",
    "gdna_high_ss_0.50_nrna_none_capture_off",
    "gdna_high_ss_0.50_nrna_none_capture_on",
)

SNAPSHOTS = {
    "baseline_state_split": ("benchmarks", "prior_mass_handoff_before"),
    "all_region_prior_mass": ("benchmarks", "prior_mass_handoff_all_region_after"),
    "capture_gated_prior_mass": (),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sim-base",
        type=Path,
        default=Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb"),
    )
    parser.add_argument("--out-dir", type=Path, default=None)
    return parser.parse_args()


def snapshot_dir(sim_base: Path, snapshot: str) -> Path:
    parts = SNAPSHOTS[snapshot]
    if not parts:
        return sim_base
    return sim_base.joinpath(*parts)


def condition_output_dir(sim_base: Path, snapshot: str, condition: str) -> Path:
    root = snapshot_dir(sim_base, snapshot)
    archived = root / "rigel_out" / condition
    if archived.exists():
        return archived
    return root / condition / "rigel_out"


def load_metrics(sim_base: Path) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for snapshot in SNAPSHOTS:
        path = snapshot_dir(sim_base, snapshot) / "condition_metrics.tsv"
        table = pd.read_csv(path, sep="\t").copy()
        table.insert(0, "snapshot", snapshot)
        rows.append(table)
    return pd.concat(rows, ignore_index=True)


def load_truth_counts(sim_base: Path) -> pd.DataFrame:
    truth = pd.read_csv(sim_base / "truth_abundances_nrna_none.tsv", sep="\t")
    total = float(truth["mrna_abundance"].sum())
    truth = truth[["transcript_id", "gene_id", "gene_name", "mrna_abundance"]].copy()
    truth["truth_mrna_count"] = np.where(
        total > 0.0,
        truth["mrna_abundance"] / total * 100000.0,
        0.0,
    )
    return truth


def load_region_diagnostics(sim_base: Path) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for condition in CONDITIONS:
        path = sim_base / "diagnostics" / condition / "per_region_calibration_truth.tsv"
        table = pd.read_csv(path, sep="\t")
        state_unexpressed = table["prior_total"] * (
            table["p_state_unexpressed_offtarget"] + table["p_state_unexpressed_capture"]
        )
        state_expressed = table["prior_total"] * (
            table["p_state_expressed_capture"] + table["p_state_expressed_offtarget"]
        )
        table = table.assign(
            condition=condition,
            capture_label="on" if condition.endswith("capture_on") else "off",
            strand_label="ss0.99" if "ss_0.99" in condition else "ss0.50",
            gdna_label="high" if condition.startswith("gdna_high") else "none",
            state_implied_unexpressed_mass=state_unexpressed,
            state_implied_expressed_mass=state_expressed,
            probe_exon=table["has_probe_overlap"].astype(bool) & (table["region_type"] == "exon"),
            expressed_state=table["state_name"].isin(
                ["expressed_capture", "expressed_offtarget"]
            ),
        )
        rows.append(table)
    return pd.concat(rows, ignore_index=True)


def load_loci(sim_base: Path) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for snapshot in SNAPSHOTS:
        for condition in CONDITIONS:
            path = condition_output_dir(sim_base, snapshot, condition) / "loci.tsv"
            table = pd.read_csv(path, sep="\t")
            table.insert(0, "condition", condition)
            table.insert(0, "snapshot", snapshot)
            rows.append(table)
    return pd.concat(rows, ignore_index=True)


def load_quant(sim_base: Path, truth: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for snapshot in SNAPSHOTS:
        for condition in CONDITIONS:
            path = condition_output_dir(sim_base, snapshot, condition) / "quant.tsv"
            table = pd.read_csv(path, sep="\t")
            table = table.merge(truth, on=["transcript_id", "gene_id", "gene_name"], how="left")
            table["truth_mrna_count"] = table["truth_mrna_count"].fillna(0.0)
            table["error"] = table["count"] - table["truth_mrna_count"]
            table["abs_error"] = table["error"].abs()
            table.insert(0, "condition", condition)
            table.insert(0, "snapshot", snapshot)
            rows.append(table)
    return pd.concat(rows, ignore_index=True)


def safe_div(numer, denom) -> np.ndarray:
    n = np.asarray(numer, dtype=np.float64)
    d = np.asarray(denom, dtype=np.float64)
    return np.divide(n, d, out=np.full_like(n, np.nan, dtype=np.float64), where=d != 0.0)


def summarize_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "snapshot",
        "condition",
        "true_gdna_fraction",
        "est_gdna_fraction",
        "gdna_delta",
        "gdna_to_rna_rate",
        "mrna_to_gdna_rate",
        "prior_global_gdna",
        "prior_global_rna",
        "mRNA_count_mard",
        "mRNA_count_spearman",
        "mRNA_count_abs_error_total",
    ]
    return metrics[cols].sort_values(["condition", "snapshot"])


def summarize_regional_truth(regions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for condition, group in regions.groupby("condition", sort=False):
        probe_exon = group["probe_exon"]
        expressed = group["expressed_state"]
        true_gdna = float(group["true_gdna_mass"].sum())
        true_rna = float(group["true_rna_mass"].sum())
        prior_gdna = float(group["prior_gdna"].sum())
        state_unexpressed = float(group["state_implied_unexpressed_mass"].sum())
        probe_true_gdna = float(group.loc[probe_exon, "true_gdna_mass"].sum())
        probe_prior_gdna = float(group.loc[probe_exon, "prior_gdna"].sum())
        probe_state_unexpressed = float(
            group.loc[probe_exon, "state_implied_unexpressed_mass"].sum()
        )
        rows.append(
            {
                "condition": condition,
                "true_gdna": true_gdna,
                "true_rna": true_rna,
                "prior_mass_gdna": prior_gdna,
                "state_implied_unexpressed_mass": state_unexpressed,
                "prior_mass_gdna_error": prior_gdna - true_gdna,
                "state_unexpressed_minus_true_gdna": state_unexpressed - true_gdna,
                "gdna_in_expressed_states": float(group.loc[expressed, "true_gdna_mass"].sum()),
                "probe_exon_true_gdna": probe_true_gdna,
                "probe_exon_prior_gdna": probe_prior_gdna,
                "probe_exon_state_unexpressed_mass": probe_state_unexpressed,
            }
        )
    return pd.DataFrame(rows)


def summarize_loci(loci: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "alpha_gdna_add",
        "alpha_rna_add",
        "prior_n_local_gdna",
        "prior_n_local_rna",
        "prior_n_other_gdna",
        "prior_n_other_rna",
        "mrna",
        "nrna",
        "gdna",
    ]
    out = loci.groupby(["snapshot", "condition"], dropna=False)[cols].sum().reset_index()
    out["alpha_gdna_share"] = safe_div(
        out["alpha_gdna_add"], out["alpha_gdna_add"] + out["alpha_rna_add"]
    )
    out["local_prior_gdna_share"] = safe_div(
        out["prior_n_local_gdna"], out["prior_n_local_gdna"] + out["prior_n_local_rna"]
    )
    out["em_gdna_share"] = safe_div(out["gdna"], out["mrna"] + out["nrna"] + out["gdna"])
    return out.sort_values(["condition", "snapshot"])


def transcript_delta_focus(quant: pd.DataFrame, condition: str) -> pd.DataFrame:
    subset = quant[quant["condition"] == condition]
    pivot = subset.pivot_table(
        index=["transcript_id", "gene_id", "gene_name", "truth_mrna_count"],
        columns="snapshot",
        values=["count", "abs_error", "error"],
        aggfunc="first",
    )
    pivot.columns = [f"{field}_{snapshot}" for field, snapshot in pivot.columns]
    pivot = pivot.reset_index()
    pivot["abs_error_change_all_region_vs_baseline"] = (
        pivot["abs_error_all_region_prior_mass"] - pivot["abs_error_baseline_state_split"]
    )
    pivot["abs_error_change_gated_vs_baseline"] = (
        pivot["abs_error_capture_gated_prior_mass"] - pivot["abs_error_baseline_state_split"]
    )
    return pivot.sort_values("abs_error_change_all_region_vs_baseline", ascending=False)


def markdown_table(table: pd.DataFrame, max_rows: int | None = None) -> str:
    view = table if max_rows is None else table.head(max_rows)
    rows = [[format_value(value) for value in row] for row in view.to_numpy()]
    headers = [str(column) for column in view.columns]
    widths = [len(header) for header in headers]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))

    def line(values: list[str]) -> str:
        return "| " + " | ".join(value.ljust(widths[idx]) for idx, value in enumerate(values)) + " |"

    divider = "| " + " | ".join("-" * width for width in widths) + " |"
    return "\n".join([line(headers), divider, *[line(row) for row in rows]])


def format_value(value: object) -> str:
    if isinstance(value, float):
        if not np.isfinite(value):
            return "n/a"
        return f"{value:.5g}"
    return str(value)


def write_markdown(
    metrics: pd.DataFrame,
    regional: pd.DataFrame,
    locus_summary: pd.DataFrame,
    focus_tx: pd.DataFrame,
    out_dir: Path,
) -> None:
    lines = ["# Prior Mass Handoff Root-Cause Diagnostics", ""]
    lines.append("## Snapshot Metrics")
    lines.append("")
    lines.append(markdown_table(metrics))
    lines.append("")
    lines.append("## Regional Truth Versus Calibration Signals")
    lines.append("")
    lines.append(markdown_table(regional))
    lines.append("")
    lines.append("## Locus Prior Summary")
    lines.append("")
    lines.append(markdown_table(locus_summary))
    lines.append("")
    lines.append("## Capture-Off Transcript Regression Focus")
    lines.append("")
    lines.append(markdown_table(focus_tx, max_rows=20))
    lines.append("")
    (out_dir / "prior_mass_handoff_root_causes.md").write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    sim_base = args.sim_base.resolve()
    out_dir = args.out_dir or (sim_base / "diagnostics" / "prior_mass_handoff_root_causes")
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics = summarize_metrics(load_metrics(sim_base))
    regional = summarize_regional_truth(load_region_diagnostics(sim_base))
    locus_summary = summarize_loci(load_loci(sim_base))
    truth = load_truth_counts(sim_base)
    quant = load_quant(sim_base, truth)
    focus_tx = transcript_delta_focus(quant, "gdna_high_ss_0.99_nrna_none_capture_off")

    metrics.to_csv(out_dir / "snapshot_metrics.tsv", sep="\t", index=False)
    regional.to_csv(out_dir / "regional_truth_vs_signals.tsv", sep="\t", index=False)
    locus_summary.to_csv(out_dir / "locus_prior_summary.tsv", sep="\t", index=False)
    focus_tx.to_csv(out_dir / "capture_off_transcript_regression_focus.tsv", sep="\t", index=False)
    write_markdown(metrics, regional, locus_summary, focus_tx, out_dir)
    print(f"[diagnose] Wrote {out_dir}")


if __name__ == "__main__":
    main()