"""Review the prior-noise v2 proposal against synthetic capture artifacts."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_SIM_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb")

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

SNAPSHOT_DIRS = {
    "baseline_state_split": "prior_mass_handoff_before",
    "capture_gated_prior_mass": ".",
    "all_region_prior_mass": "prior_mass_handoff_all_region_after",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sim-base", type=Path, default=DEFAULT_SIM_BASE)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--null-sigma", type=float, default=3.0)
    parser.add_argument("--prior-data-fraction", type=float, default=0.25)
    return parser.parse_args()


def safe_div(numer, denom):
    n = np.asarray(numer, dtype=np.float64)
    d = np.asarray(denom, dtype=np.float64)
    return np.divide(n, d, out=np.full_like(n, np.nan, dtype=np.float64), where=d > 0.0)


def finite_quantile(values: pd.Series | np.ndarray, q: float) -> float:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.quantile(arr, q))


def load_condition_metrics(sim_base: Path) -> pd.DataFrame:
    return pd.read_csv(sim_base / "condition_metrics.tsv", sep="\t").set_index("condition")


def load_region_table(sim_base: Path, condition: str) -> pd.DataFrame:
    path = sim_base / "diagnostics" / condition / "per_region_calibration_truth.tsv"
    table = pd.read_csv(path, sep="\t")
    if "p_state_unexpressed_offtarget" not in table.columns and "p_state_background" in table:
        table = table.rename(
            columns={
                "p_state_background": "p_state_unexpressed_offtarget",
                "p_state_gdna_only_capture": "p_state_unexpressed_capture",
            }
        )
    return table


def load_loci(sim_base: Path, snapshot: str, condition: str) -> pd.DataFrame:
    run_dir = SNAPSHOT_DIRS[snapshot]
    if run_dir == ".":
        path = sim_base / condition / "rigel_out" / "loci.tsv"
    else:
        path = sim_base / "benchmarks" / run_dir / "rigel_out" / condition / "loci.tsv"
    return pd.read_csv(path, sep="\t")


def load_quant(sim_base: Path, snapshot: str, condition: str) -> pd.DataFrame:
    run_dir = SNAPSHOT_DIRS[snapshot]
    if run_dir == ".":
        path = sim_base / condition / "rigel_out" / "quant.tsv"
    else:
        path = sim_base / "benchmarks" / run_dir / "rigel_out" / condition / "quant.tsv"
    return pd.read_csv(path, sep="\t")


def exposure_review(sim_base: Path, metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for condition in CONDITIONS:
        table = load_region_table(sim_base, condition)
        rho = float(metrics.loc[condition, "rho_off"])
        off_expected = rho * np.maximum(table["contained_leff"].to_numpy(np.float64), 0.0)
        ratios = pd.DataFrame(
            {
                "off_expected": off_expected,
                "current_a_r": table["A_r"].to_numpy(np.float64),
                "observed_ratio": safe_div(table["observed_compatible_count"], off_expected),
                "prior_gdna_ratio": safe_div(table["prior_gdna"], off_expected),
                "true_gdna_ratio": safe_div(table["true_gdna_mass"], off_expected),
            }
        )
        table = pd.concat([table.copy(), ratios], axis=1)
        masks = {
            "all_regions": np.ones(len(table), dtype=bool),
            "probe_exons": table["has_probe_overlap"].astype(bool)
            & (table["region_type"] == "exon"),
            "probe_exons_true_gdna": table["has_probe_overlap"].astype(bool)
            & (table["region_type"] == "exon")
            & (table["true_gdna_mass"] > 0.0),
            "probe_exons_rna_only": table["has_probe_overlap"].astype(bool)
            & (table["region_type"] == "exon")
            & (table["true_gdna_mass"] <= 0.0)
            & (table["true_rna_mass"] > 0.0),
            "expressed_no_gdna": (table["true_rna_mass"] > 0.0)
            & (table["true_gdna_mass"] <= 0.0),
        }
        for group_name, mask in masks.items():
            group = table.loc[mask].copy()
            if group.empty:
                continue
            rows.append(
                {
                    "condition": condition,
                    "group": group_name,
                    "n_regions": int(len(group)),
                    "true_gdna": float(group["true_gdna_mass"].sum()),
                    "true_rna": float(group["true_rna_mass"].sum()),
                    "off_expected": float(group["off_expected"].sum()),
                    "observed_ratio_sum": float(
                        group["observed_compatible_count"].sum()
                        / max(group["off_expected"].sum(), 1.0e-12)
                    ),
                    "prior_gdna_ratio_sum": float(
                        group["prior_gdna"].sum() / max(group["off_expected"].sum(), 1.0e-12)
                    ),
                    "true_gdna_ratio_sum": float(
                        group["true_gdna_mass"].sum()
                        / max(group["off_expected"].sum(), 1.0e-12)
                    ),
                    "current_a_r_p50": finite_quantile(group["current_a_r"], 0.50),
                    "current_a_r_p95": finite_quantile(group["current_a_r"], 0.95),
                    "observed_ratio_p50": finite_quantile(group["observed_ratio"], 0.50),
                    "observed_ratio_p95": finite_quantile(group["observed_ratio"], 0.95),
                    "prior_gdna_ratio_p50": finite_quantile(group["prior_gdna_ratio"], 0.50),
                    "prior_gdna_ratio_p95": finite_quantile(group["prior_gdna_ratio"], 0.95),
                    "true_gdna_ratio_p50": finite_quantile(group["true_gdna_ratio"], 0.50),
                    "true_gdna_ratio_p95": finite_quantile(group["true_gdna_ratio"], 0.95),
                }
            )
    return pd.DataFrame(rows)


def prior_null_cap_review(
    sim_base: Path,
    metrics: pd.DataFrame,
    *,
    null_sigma: float,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for condition in CONDITIONS:
        table = load_region_table(sim_base, condition)
        strand_specificity = float(metrics.loc[condition, "strand_model_est"])
        q_error = min(strand_specificity, 1.0 - strand_specificity)
        kappa = metrics.loc[condition, "strand_channel_kappa_d"]
        kappa = float(kappa) if pd.notna(kappa) else 2.0
        sigma_inflation = math.sqrt(1.0 + 1.0 / max(kappa, 1.0))
        cap = (
            2.0
            * table["prior_total"].to_numpy(np.float64)
            * max(q_error, 1.0e-12)
            * (1.0 + null_sigma * sigma_inflation)
        )
        prior_gdna = table["prior_gdna"].to_numpy(np.float64)
        ineligible = (table["strand_flags"].fillna(0).astype(np.uint16).to_numpy() & 0x1) != 0
        expressed = table["p_expressed"].to_numpy(np.float64) >= 0.8
        capped_everywhere = np.minimum(prior_gdna, cap)
        capped_eligible_only = np.where(ineligible, prior_gdna, np.minimum(prior_gdna, cap))
        capped_expressed_only = np.where(expressed, np.minimum(prior_gdna, cap), prior_gdna)
        cap_binds = prior_gdna > cap
        rows.append(
            {
                "condition": condition,
                "q_error": q_error,
                "cap_factor_per_count": 2.0 * q_error * (1.0 + null_sigma * sigma_inflation),
                "true_gdna": float(table["true_gdna_mass"].sum()),
                "prior_gdna": float(np.sum(prior_gdna)),
                "capped_everywhere_gdna": float(np.sum(capped_everywhere)),
                "capped_eligible_only_gdna": float(np.sum(capped_eligible_only)),
                "capped_expressed_only_gdna": float(np.sum(capped_expressed_only)),
                "n_regions_cap_binds": int(np.count_nonzero(cap_binds)),
                "true_gdna_in_cap_bound_regions": float(table.loc[cap_binds, "true_gdna_mass"].sum()),
                "rna_in_cap_bound_regions": float(table.loc[cap_binds, "true_rna_mass"].sum()),
            }
        )
    return pd.DataFrame(rows)


def ess_and_floor_review(sim_base: Path, *, prior_data_fraction: float) -> pd.DataFrame:
    targets = [
        ("gdna_high_ss_0.99_nrna_none_capture_on", 1),
        ("gdna_high_ss_0.99_nrna_none_capture_off", 7),
    ]
    rows: list[dict[str, object]] = []
    for condition, locus_id in targets:
        for snapshot in SNAPSHOT_DIRS:
            loci = load_loci(sim_base, snapshot, condition)
            row = loci.loc[loci["locus_id"] == locus_id]
            if row.empty:
                continue
            item = row.iloc[0]
            data_cap = prior_data_fraction * float(item["n_em_fragments"])
            unspliced_cap = float(item["prior_unspliced_total"])
            proposed_cap_no_precision = min(data_cap, unspliced_cap)
            k_rna = int(item.get("n_transcripts", 0))
            rna_floor = float(item["alpha_rna_add"]) / max(k_rna, 1)
            rows.append(
                {
                    "condition": condition,
                    "locus_id": locus_id,
                    "snapshot": snapshot,
                    "n_em_fragments": float(item["n_em_fragments"]),
                    "prior_n_local_gdna": float(item["prior_n_local_gdna"]),
                    "prior_n_local_rna": float(item["prior_n_local_rna"]),
                    "alpha_gdna_add": float(item["alpha_gdna_add"]),
                    "alpha_rna_add": float(item["alpha_rna_add"]),
                    "prior_ess_final": float(item["prior_ess_final"]),
                    "data_fraction_cap": data_cap,
                    "proposed_cap_no_precision": proposed_cap_no_precision,
                    "n_rna_components": k_rna,
                    "dirichlet_floor_per_rna_component": rna_floor,
                    "gdna_eff_len_weight_ratio": float(item["gdna_eff_len_weight_ratio"]),
                    "gdna": float(item["gdna"]),
                    "mrna": float(item["mrna"]),
                    "nrna": float(item["nrna"]),
                }
            )
    return pd.DataFrame(rows)


def quant_focus(sim_base: Path) -> pd.DataFrame:
    condition = "gdna_high_ss_0.99_nrna_none_capture_off"
    transcripts = {"GENE0008.2", "GENE0008.3", "GENE0008.4"}
    rows: list[pd.DataFrame] = []
    for snapshot in SNAPSHOT_DIRS:
        quant = load_quant(sim_base, snapshot, condition)
        focus = quant.loc[quant["transcript_id"].isin(transcripts)].copy()
        focus.insert(0, "snapshot", snapshot)
        rows.append(focus)
    return pd.concat(rows, ignore_index=True)


def fmt_float(value: object, digits: int = 4) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(number):
        return "nan"
    return f"{number:.{digits}g}"


def markdown_table(df: pd.DataFrame, max_rows: int = 20) -> str:
    view = df.head(max_rows).copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(fmt_float)
    columns = [str(col) for col in view.columns]
    rows = [[str(value) for value in row] for row in view.to_numpy()]
    widths = [len(col) for col in columns]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))

    def line(values: list[str]) -> str:
        return "| " + " | ".join(value.ljust(widths[idx]) for idx, value in enumerate(values)) + " |"

    divider = "| " + " | ".join("-" * width for width in widths) + " |"
    return "\n".join([line(columns), divider, *[line(row) for row in rows]])


def write_report(
    out_dir: Path,
    exposure: pd.DataFrame,
    caps: pd.DataFrame,
    ess: pd.DataFrame,
    quant: pd.DataFrame,
) -> None:
    lines = ["# Prior Noise v2 Review Diagnostics", ""]
    lines.append("## Exposure Counterfactuals")
    lines.append("")
    focus = exposure[
        exposure["group"].isin(["probe_exons", "probe_exons_rna_only", "probe_exons_true_gdna"])
    ].copy()
    lines.append(markdown_table(focus, 40))
    lines.append("")
    lines.append("## Proposed Null Cap Counterfactual")
    lines.append("")
    lines.append(markdown_table(caps, 20))
    lines.append("")
    lines.append("## ESS and Dirichlet Floor Magnitudes")
    lines.append("")
    lines.append(markdown_table(ess, 20))
    lines.append("")
    lines.append("## GENE0008 Transcript Focus")
    lines.append("")
    cols = [col for col in ["snapshot", "transcript_id", "gene_id", "count", "locus_id"] if col in quant]
    lines.append(markdown_table(quant[cols], 20))
    lines.append("")
    lines.append("## Main Readouts")
    lines.append("")
    lines.append(
        "- Raw observed-mass exposure is high in RNA-only probe exons, so it is not "
        "source-specific enough to use directly as gDNA exposure."
    )
    lines.append(
        "- Per-region prior-gDNA exposure is source-specific but still too spiky in "
        "no-gDNA SS=0.99 probe exons; exposure needs aggregate/Gamma shrinkage and "
        "pool evidence, not raw region ratios."
    )
    lines.append(
        "- The proposed hard null cap would bind regions containing real gDNA in high-gDNA "
        "SS=0.99 conditions; it should become an excess-over-null shrinkage rule, not a "
        "global min(gdna, cap)."
    )
    lines.append(
        "- Locus-scaled ESS is directionally right; the data-fraction cap raises Locus 1 "
        "from 3k toward 25k but only modestly changes Locus 7."
    )
    lines.append(
        "- A symmetric RNA floor prevents exact zeros, but its guaranteed mass is only "
        "alpha_rna_add / K_rna; it is a stability guard, not a full isoform-resolution fix."
    )
    (out_dir / "prior_noise_v2_review.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    sim_base = args.sim_base.resolve()
    out_dir = args.out_dir or (sim_base / "diagnostics" / "prior_noise_v2_review")
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics = load_condition_metrics(sim_base)
    exposure = exposure_review(sim_base, metrics)
    caps = prior_null_cap_review(sim_base, metrics, null_sigma=float(args.null_sigma))
    ess = ess_and_floor_review(sim_base, prior_data_fraction=float(args.prior_data_fraction))
    quant = quant_focus(sim_base)

    exposure.to_csv(out_dir / "exposure_counterfactuals.tsv", sep="\t", index=False)
    caps.to_csv(out_dir / "prior_null_cap_counterfactual.tsv", sep="\t", index=False)
    ess.to_csv(out_dir / "ess_and_floor_counterfactual.tsv", sep="\t", index=False)
    quant.to_csv(out_dir / "gene0008_quant_focus.tsv", sep="\t", index=False)
    write_report(out_dir, exposure, caps, ess, quant)
    print(f"[review] wrote {out_dir}")


if __name__ == "__main__":
    main()