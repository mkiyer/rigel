"""Audit latent state labels against oracle region truth in the capture suite."""

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sim-base",
        type=Path,
        default=Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb"),
    )
    return parser.parse_args()


def safe_div(numer, denom):
    n = np.asarray(numer, dtype=np.float64)
    d = np.asarray(denom, dtype=np.float64)
    return np.divide(n, d, out=np.full_like(n, np.nan, dtype=np.float64), where=d != 0.0)


def fmt_pct(value: float) -> str:
    if not np.isfinite(value):
        return "n/a"
    return f"{100.0 * value:.2f}%"


def fmt_count(value: float) -> str:
    if not np.isfinite(value):
        return "n/a"
    return f"{value:,.0f}"


def fmt_float(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "n/a"
    return f"{value:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    widths = [len(header) for header in headers]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))

    def line(values: list[str]) -> str:
        return "| " + " | ".join(value.ljust(widths[idx]) for idx, value in enumerate(values)) + " |"

    divider = "| " + " | ".join("-" * width for width in widths) + " |"
    return "\n".join([line(headers), divider, *[line(row) for row in rows]])


def load_region_table(sim_base: Path, condition: str) -> pd.DataFrame:
    path = sim_base / "diagnostics" / condition / "per_region_calibration_truth.tsv"
    if not path.exists():
        raise FileNotFoundError(path)
    table = pd.read_csv(path, sep="\t")
    diagnostics = pd.DataFrame(
        {
            "condition": condition,
            "capture_label": "on" if condition.endswith("capture_on") else "off",
            "strand_label": "ss0.99" if "ss_0.99" in condition else "ss0.50",
            "gdna_label": "high" if condition.startswith("gdna_high") else "none",
            "state_implied_gdna": table["prior_total"]
            * (table["p_state_background"] + table["p_state_gdna_only_capture"]),
            "state_implied_rna": table["prior_total"]
            * (table["p_state_expressed_capture"] + table["p_state_expressed_offtarget"]),
            "mixed_truth": (table["true_gdna_mass"] > 0.0) & (table["true_rna_mass"] > 0.0),
            "expressed_state": table["state_name"].isin(
                ["expressed_capture", "expressed_offtarget"]
            ),
            "gdna_state": table["state_name"].isin(["background", "gdna_only_capture"]),
            "probe_exon": table["has_probe_overlap"].astype(bool)
            & (table["region_type"] == "exon"),
        }
    )
    return pd.concat([table, diagnostics], axis=1).copy()


def summarize_conditions(table: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for condition, group in table.groupby("condition", sort=False):
        true_gdna_total = float(group["true_gdna_mass"].sum())
        true_rna_total = float(group["true_rna_mass"].sum())
        prior_gdna_total = float(group["prior_gdna"].sum())
        state_implied_gdna_total = float(group["state_implied_gdna"].sum())
        expressed = group["expressed_state"]
        mixed = group["mixed_truth"]
        probe_exon = group["probe_exon"]
        rows.append(
            {
                "condition": condition,
                "true_gdna_total": true_gdna_total,
                "true_rna_total": true_rna_total,
                "prior_gdna_total": prior_gdna_total,
                "state_implied_gdna_total": state_implied_gdna_total,
                "prior_gdna_error": prior_gdna_total - true_gdna_total,
                "state_implied_gdna_error": state_implied_gdna_total - true_gdna_total,
                "true_gdna_in_expressed_states": float(group.loc[expressed, "true_gdna_mass"].sum()),
                "fraction_gdna_in_expressed_states": float(
                    safe_div(group.loc[expressed, "true_gdna_mass"].sum(), true_gdna_total)
                ),
                "mixed_region_true_gdna": float(group.loc[mixed, "true_gdna_mass"].sum()),
                "mixed_region_true_rna": float(group.loc[mixed, "true_rna_mass"].sum()),
                "mixed_gdna_in_expressed_states": float(
                    group.loc[mixed & expressed, "true_gdna_mass"].sum()
                ),
                "probe_exon_true_gdna": float(group.loc[probe_exon, "true_gdna_mass"].sum()),
                "probe_exon_true_rna": float(group.loc[probe_exon, "true_rna_mass"].sum()),
                "probe_exon_prior_gdna": float(group.loc[probe_exon, "prior_gdna"].sum()),
                "probe_exon_state_implied_gdna": float(
                    group.loc[probe_exon, "state_implied_gdna"].sum()
                ),
                "probe_exon_p_expressed_weighted": weighted_mean(
                    group.loc[probe_exon, "p_expressed"],
                    group.loc[probe_exon, "true_total_mass"],
                ),
                "probe_exon_p_captured_weighted": weighted_mean(
                    group.loc[probe_exon, "p_captured"],
                    group.loc[probe_exon, "true_total_mass"],
                ),
                "probe_exon_spliced_count": float(group.loc[probe_exon, "spliced_count"].sum()),
            }
        )
    return pd.DataFrame(rows)


def summarize_state_groups(table: pd.DataFrame) -> pd.DataFrame:
    keys = ["condition", "state_name", "region_type", "has_probe_overlap"]
    cols = [
        "true_gdna_mass",
        "true_rna_mass",
        "prior_gdna",
        "prior_rna",
        "state_implied_gdna",
        "state_implied_rna",
        "observed_compatible_count",
        "spliced_count",
    ]
    out = table.groupby(keys, dropna=False)[cols].sum().reset_index()
    out["true_gdna_fraction"] = safe_div(
        out["true_gdna_mass"],
        out["true_gdna_mass"] + out["true_rna_mass"],
    )
    out["prior_gdna_fraction"] = safe_div(out["prior_gdna"], out["prior_gdna"] + out["prior_rna"])
    out["state_implied_gdna_fraction"] = safe_div(
        out["state_implied_gdna"],
        out["state_implied_gdna"] + out["state_implied_rna"],
    )
    return out.sort_values(["condition", "true_gdna_mass"], ascending=[True, False])


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    v = np.asarray(values, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    total = float(np.sum(w))
    if total <= 0.0:
        return float("nan")
    return float(np.sum(v * w) / total)


def write_report(summary: pd.DataFrame, state_groups: pd.DataFrame, out_dir: Path) -> None:
    rows: list[list[str]] = []
    for _, row in summary.iterrows():
        rows.append(
            [
                str(row["condition"]),
                fmt_count(row["true_gdna_total"]),
                fmt_count(row["prior_gdna_total"]),
                fmt_count(row["state_implied_gdna_total"]),
                fmt_pct(row["fraction_gdna_in_expressed_states"]),
                fmt_count(row["probe_exon_true_gdna"]),
                fmt_count(row["probe_exon_prior_gdna"]),
                fmt_count(row["probe_exon_state_implied_gdna"]),
                fmt_float(row["probe_exon_p_expressed_weighted"], 3),
                fmt_float(row["probe_exon_p_captured_weighted"], 3),
                fmt_count(row["probe_exon_spliced_count"]),
            ]
        )

    focus = state_groups[
        (state_groups["condition"].isin(
            [
                "gdna_high_ss_0.99_nrna_none_capture_on",
                "gdna_high_ss_0.50_nrna_none_capture_on",
            ]
        ))
        & (state_groups["region_type"] == "exon")
        & (state_groups["has_probe_overlap"])
    ].copy()
    focus_rows: list[list[str]] = []
    for _, row in focus.iterrows():
        focus_rows.append(
            [
                str(row["condition"]),
                str(row["state_name"]),
                fmt_count(row["true_gdna_mass"]),
                fmt_count(row["true_rna_mass"]),
                fmt_pct(row["true_gdna_fraction"]),
                fmt_count(row["prior_gdna"]),
                fmt_pct(row["prior_gdna_fraction"]),
                fmt_count(row["state_implied_gdna"]),
                fmt_pct(row["state_implied_gdna_fraction"]),
                fmt_count(row["spliced_count"]),
            ]
        )

    lines = ["# Latent State Label Audit", ""]
    lines.append("## Condition Summary")
    lines.append("")
    lines.append(
        markdown_table(
            [
                "Condition",
                "True gDNA",
                "Prior gDNA",
                "State-implied gDNA",
                "gDNA in expressed states",
                "Probe-exon true gDNA",
                "Probe-exon prior gDNA",
                "Probe-exon state gDNA",
                "Probe-exon p_expr",
                "Probe-exon p_cap",
                "Probe-exon spliced",
            ],
            rows,
        )
    )
    lines.append("")
    lines.append("## Probe-Overlap Exon Focus")
    lines.append("")
    lines.append(
        markdown_table(
            [
                "Condition",
                "State",
                "True gDNA",
                "True RNA",
                "True gDNA frac",
                "Prior gDNA",
                "Prior gDNA frac",
                "State-implied gDNA",
                "State gDNA frac",
                "Spliced count",
            ],
            focus_rows,
        )
    )
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append(
        "The latent labels are inaccurate as pool labels because the four-state model "
        "does not include a mixed expressed-plus-gDNA state. Spliced/RNA-lower evidence "
        "pushes mixed exonic regions into expressed states, while the capture/density "
        "factors decide capture/offtarget status. The mass split must therefore come "
        "from `prior_mass`, not from `state_name` or `p_state_background + p_state_gdna_only_capture`."
    )
    (out_dir / "latent_state_label_audit.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    sim_base = args.sim_base.resolve()
    out_dir = sim_base / "diagnostics" / "latent_state_label_audit"
    out_dir.mkdir(parents=True, exist_ok=True)
    table = pd.concat([load_region_table(sim_base, condition) for condition in CONDITIONS], ignore_index=True)
    summary = summarize_conditions(table)
    state_groups = summarize_state_groups(table)
    table.to_csv(out_dir / "latent_state_region_table.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "latent_state_condition_summary.tsv", sep="\t", index=False)
    state_groups.to_csv(out_dir / "latent_state_group_summary.tsv", sep="\t", index=False)
    write_report(summary, state_groups, out_dir)
    print(f"[latent] Wrote {out_dir}")


if __name__ == "__main__":
    main()