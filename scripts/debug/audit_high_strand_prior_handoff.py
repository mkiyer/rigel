"""Audit calibration-to-EM prior handoff for high-gDNA strand-specific capture runs."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


CONDITIONS = (
    "gdna_high_ss_0.99_nrna_none_capture_off",
    "gdna_high_ss_0.99_nrna_none_capture_on",
)

STATE_COLUMNS = (
    "p_state_background",
    "p_state_gdna_only_capture",
    "p_state_expressed_capture",
    "p_state_expressed_offtarget",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sim-base",
        type=Path,
        default=Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb"),
    )
    parser.add_argument("--out-dir", type=Path, default=None)
    return parser.parse_args()


def entropy_weight(states: np.ndarray) -> np.ndarray:
    p = np.asarray(states, dtype=np.float64)
    p = np.clip(p, 0.0, 1.0)
    row_sum = p.sum(axis=1)
    p = np.divide(p, row_sum[:, None], out=np.zeros_like(p), where=row_sum[:, None] > 0.0)
    xlogx = np.zeros_like(p)
    positive = p > 0.0
    xlogx[positive] = p[positive] * np.log(p[positive])
    entropy = -np.sum(xlogx, axis=1)
    return np.clip(1.0 - entropy / math.log(4.0), 0.0, 1.0)


def load_region_table(sim_base: Path, condition: str) -> pd.DataFrame:
    path = sim_base / "diagnostics" / condition / "per_region_calibration_truth.tsv"
    if not path.exists():
        raise FileNotFoundError(path)
    table = pd.read_csv(path, sep="\t")
    missing = [col for col in STATE_COLUMNS if col not in table.columns]
    if missing:
        raise ValueError(f"{path} missing state columns: {missing}")
    states = table.loc[:, STATE_COLUMNS].to_numpy(dtype=np.float64)
    weight = entropy_weight(states)
    q_gdna_state = states[:, 0] + states[:, 1]
    q_rna_state = states[:, 2] + states[:, 3]
    q_total = q_gdna_state + q_rna_state
    q_gdna_state = np.divide(q_gdna_state, q_total, out=np.zeros_like(q_gdna_state), where=q_total > 0.0)
    q_rna_state = np.divide(q_rna_state, q_total, out=np.zeros_like(q_rna_state), where=q_total > 0.0)
    unspliced = table["prior_total"].to_numpy(dtype=np.float64)
    diagnostics = pd.DataFrame(
        {
            "entropy_weight": weight,
            "state_split_gdna_mass": unspliced * weight * q_gdna_state,
            "state_split_rna_mass": unspliced * weight * q_rna_state,
            "strand_split_gdna_mass": table["prior_gdna"].to_numpy(dtype=np.float64) * weight,
            "strand_split_rna_mass": table["prior_rna"].to_numpy(dtype=np.float64) * weight,
            "probe_bin": np.where(table["has_probe_overlap"], "probe_overlap", "no_probe"),
            "captured_bin": np.where(
                table["p_captured"] >= 0.5,
                "p_captured>=0.5",
                "p_captured<0.5",
            ),
            "condition": condition,
        }
    )
    diagnostics["discarded_strand_gdna_mass"] = (
        diagnostics["strand_split_gdna_mass"] - diagnostics["state_split_gdna_mass"]
    )
    return pd.concat([table, diagnostics], axis=1).copy()


def summarize_regions(table: pd.DataFrame) -> pd.DataFrame:
    keys = ["condition", "probe_bin", "captured_bin", "state_name", "region_type"]
    cols = [
        "true_gdna_mass",
        "true_rna_mass",
        "prior_gdna",
        "prior_rna",
        "strand_split_gdna_mass",
        "strand_split_rna_mass",
        "state_split_gdna_mass",
        "state_split_rna_mass",
        "discarded_strand_gdna_mass",
        "observed_compatible_count",
        "true_gdna_assigned_rna_major",
        "true_gdna_assigned_gdna_major",
        "true_rna_assigned_gdna_major",
        "true_rna_assigned_rna_major",
    ]
    grouped = table.groupby(keys, dropna=False)[cols].sum().reset_index()
    grouped.insert(0, "n_regions", table.groupby(keys, dropna=False).size().to_numpy())
    grouped["state_vs_strand_gdna_ratio"] = safe_div(
        grouped["state_split_gdna_mass"], grouped["strand_split_gdna_mass"]
    )
    grouped["gdna_to_rna_major_rate"] = safe_div(
        grouped["true_gdna_assigned_rna_major"],
        grouped["true_gdna_assigned_rna_major"] + grouped["true_gdna_assigned_gdna_major"],
    )
    return grouped.sort_values(
        ["discarded_strand_gdna_mass", "true_gdna_mass"], ascending=False
    )


def load_locus_table(sim_base: Path, condition: str) -> pd.DataFrame:
    path = sim_base / condition / "rigel_out" / "loci.tsv"
    if not path.exists():
        raise FileNotFoundError(path)
    table = pd.read_csv(path, sep="\t")
    table["condition"] = condition
    return table


def summarize_loci(table: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "n_em_fragments",
        "count_unambig",
        "mrna",
        "nrna",
        "gdna",
        "alpha_gdna_add",
        "alpha_rna_add",
        "prior_unspliced_total",
        "prior_n_local_gdna",
        "prior_n_local_rna",
        "prior_n_other_gdna",
        "prior_n_other_rna",
        "prior_ess_final",
    ]
    grouped = table.groupby("condition", dropna=False)[cols].sum().reset_index()
    grouped["alpha_gdna_share"] = safe_div(
        grouped["alpha_gdna_add"], grouped["alpha_gdna_add"] + grouped["alpha_rna_add"]
    )
    grouped["local_gdna_share"] = safe_div(
        grouped["prior_n_local_gdna"], grouped["prior_n_local_gdna"] + grouped["prior_n_local_rna"]
    )
    grouped["em_gdna_share"] = safe_div(grouped["gdna"], grouped["gdna"] + grouped["mrna"] + grouped["nrna"])
    return grouped


def safe_div(numer, denom):
    n = np.asarray(numer, dtype=np.float64)
    d = np.asarray(denom, dtype=np.float64)
    return np.divide(n, d, out=np.zeros_like(n, dtype=np.float64), where=d != 0.0)


def markdown_table(df: pd.DataFrame, max_rows: int) -> str:
    view = df.head(max_rows).copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda value: f"{value:.5g}")
    columns = [str(col) for col in view.columns]
    rows = [[str(value) for value in row] for row in view.to_numpy()]
    widths = [len(col) for col in columns]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))

    def fmt(row: list[str]) -> str:
        return "| " + " | ".join(value.ljust(widths[idx]) for idx, value in enumerate(row)) + " |"

    return "\n".join(
        [fmt(columns), "| " + " | ".join("-" * width for width in widths) + " |"]
        + [fmt(row) for row in rows]
    )


def main() -> None:
    args = parse_args()
    sim_base = args.sim_base.resolve()
    out_dir = args.out_dir or (sim_base / "diagnostics" / "high_strand_prior_handoff")
    out_dir.mkdir(parents=True, exist_ok=True)

    region_tables = [load_region_table(sim_base, condition) for condition in CONDITIONS]
    region_table = pd.concat(region_tables, ignore_index=True)
    region_summary = summarize_regions(region_table)
    region_summary.to_csv(out_dir / "region_state_vs_strand_handoff.tsv", sep="\t", index=False)

    locus_tables = [load_locus_table(sim_base, condition) for condition in CONDITIONS]
    locus_table = pd.concat(locus_tables, ignore_index=True)
    locus_summary = summarize_loci(locus_table)
    locus_summary.to_csv(out_dir / "locus_prior_summary.tsv", sep="\t", index=False)
    locus_table.to_csv(out_dir / "locus_prior_details.tsv", sep="\t", index=False)

    lines = ["# High-gDNA SS 0.99 Prior Handoff Audit", ""]
    lines.append("## Locus Prior Summary")
    lines.append("")
    lines.append(markdown_table(locus_summary, 10))
    lines.append("")
    lines.append("## Region State Split vs Strand Split")
    lines.append("")
    lines.append(markdown_table(region_summary, 16))
    lines.append("")
    lines.append("## Interpretation Aid")
    lines.append("")
    lines.append(
        "`strand_split_gdna_mass` is what the strand-aware calibration prior mass says after "
        "the same entropy weight used by adaptive priors. `state_split_gdna_mass` is what "
        "the current adaptive prior code reconstructs from latent state probabilities."
    )
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n")
    print(f"[audit] Wrote {out_dir}")


if __name__ == "__main__":
    main()