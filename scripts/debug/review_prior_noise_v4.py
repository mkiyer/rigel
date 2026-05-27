"""Diagnostic facts for the prior-noise v4 redesign."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sim-base", type=Path, default=DEFAULT_SIM_BASE)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--tail-prob", type=float, default=0.999)
    return parser.parse_args()


def load_condition_metrics(sim_base: Path) -> pd.DataFrame:
    return pd.read_csv(sim_base / "condition_metrics.tsv", sep="\t").set_index("condition")


def load_region_table(sim_base: Path, condition: str) -> pd.DataFrame:
    return pd.read_csv(sim_base / "diagnostics" / condition / "per_region_calibration_truth.tsv", sep="\t")


def finite_quantile(values: np.ndarray | pd.Series, q: float) -> float:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.quantile(arr, q))


def safe_ratio(numer: float, denom: float) -> float:
    return float(numer / denom) if denom > 0.0 else float("nan")


def kappa_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "strand_specificity",
        "strand_model_est",
        "strand_channel_kappa_d",
        "background_n_seed_regions",
        "boundary_sweep_regions_with_evidence",
        "prior_global_gdna",
        "prior_global_rna",
        "est_gdna_fraction",
        "mRNA_count_mard",
    ]
    return metrics.loc[list(CONDITIONS), cols].reset_index()


def exposure_geometry(sim_base: Path, metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for condition in CONDITIONS:
        if "capture_on" not in condition:
            continue
        table = load_region_table(sim_base, condition)
        rho = float(metrics.loc[condition, "rho_off"])
        off_expected = rho * np.maximum(table["contained_leff"].to_numpy(np.float64), 0.0)
        mask = table["has_probe_overlap"].astype(bool) & (table["region_type"] == "exon")
        group = table.loc[mask]
        off_group = off_expected[mask.to_numpy()]
        lengths = (group["end"].to_numpy(np.float64) - group["start"].to_numpy(np.float64)).clip(
            min=0.0
        )
        rows.append(
            {
                "condition": condition,
                "n_probe_exon_regions": int(mask.sum()),
                "probe_bp_sum": float(lengths.sum()),
                "probe_bp_p50": finite_quantile(lengths, 0.50),
                "probe_bp_p95": finite_quantile(lengths, 0.95),
                "off_expected_sum": float(off_group.sum()),
                "prior_gdna_sum": float(group["prior_gdna"].sum()),
                "true_gdna_sum": float(group["true_gdna_mass"].sum()),
                "prior_gdna_to_off_expected": safe_ratio(
                    float(group["prior_gdna"].sum()), float(off_group.sum())
                ),
                "true_gdna_to_off_expected": safe_ratio(
                    float(group["true_gdna_mass"].sum()), float(off_group.sum())
                ),
                "current_a_r_p50": finite_quantile(group["A_r"], 0.50),
                "current_a_r_p95": finite_quantile(group["A_r"], 0.95),
            }
        )
    return pd.DataFrame(rows)


def low_n_tail_table(tail_prob: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for q_error in (0.01, 0.05, 0.10):
        for n in (1, 5, 10, 20, 50, 100, 200):
            mean = n * q_error
            sd = math.sqrt(n * q_error * (1.0 - q_error))
            normal_3sigma = mean + 3.0 * sd
            rows.append(
                {
                    "n": n,
                    "q_error": q_error,
                    "mean": mean,
                    "normal_3sigma": normal_3sigma,
                    "binom_ppf_tail": float(binom.ppf(tail_prob, n, q_error)),
                }
            )
    return pd.DataFrame(rows)


def ess_scale_table() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for n_em in (100, 1_000, 10_000, 100_000, 1_000_000):
        rows.append(
            {
                "n_em_fragments": n_em,
                "current_cap": 3000.0,
                "linear_0p25N": 0.25 * n_em,
                "sqrt_N": math.sqrt(n_em),
                "sqrt_3000_at_100k": 3000.0 * math.sqrt(n_em / 100_000.0),
            }
        )
    return pd.DataFrame(rows)


def fmt_float(value: object, digits: int = 4) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(number):
        return "nan"
    return f"{number:.{digits}g}"


def markdown_table(df: pd.DataFrame, max_rows: int = 30) -> str:
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
    kappa: pd.DataFrame,
    exposure: pd.DataFrame,
    tails: pd.DataFrame,
    ess: pd.DataFrame,
) -> None:
    lines = ["# Prior Noise v4 Diagnostics", ""]
    lines.append("## Kappa and Source-Mass Summary")
    lines.append("")
    lines.append(markdown_table(kappa, 20))
    lines.append("")
    lines.append("## Capture Probe Exposure Geometry")
    lines.append("")
    lines.append(markdown_table(exposure, 20))
    lines.append("")
    lines.append("## Exact Binomial Tail vs Normal Approximation")
    lines.append("")
    lines.append(markdown_table(tails, 30))
    lines.append("")
    lines.append("## ESS Scale Comparison")
    lines.append("")
    lines.append(markdown_table(ess, 20))
    lines.append("")
    lines.append("## Main Readouts")
    lines.append("")
    lines.append("- The high-gDNA conditions report kappa_d at the 1e6 maximum, consistent with binomial or sub-binomial simulated gDNA strand balance.")
    lines.append("- Probe-exon off-target expectations can be tiny, so per-region or per-panel ratios are unstable unless source mass is first made trustworthy.")
    lines.append("- Exact binomial tails differ materially from the continuous 3-sigma approximation at low N; the v4 RNA-noise guard should use exact predictive tails.")
    lines.append("- A 0.25N ESS cap is intentionally large in mega-loci and must be benchmarked as a behavioral change, not treated as a theorem.")
    (out_dir / "prior_noise_v4_diagnostics.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    sim_base = args.sim_base.resolve()
    out_dir = args.out_dir or (sim_base / "diagnostics" / "prior_noise_v4_review")
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics = load_condition_metrics(sim_base)
    kappa = kappa_summary(metrics)
    exposure = exposure_geometry(sim_base, metrics)
    tails = low_n_tail_table(float(args.tail_prob))
    ess = ess_scale_table()
    kappa.to_csv(out_dir / "kappa_summary.tsv", sep="\t", index=False)
    exposure.to_csv(out_dir / "capture_exposure_geometry.tsv", sep="\t", index=False)
    tails.to_csv(out_dir / "exact_binomial_tails.tsv", sep="\t", index=False)
    ess.to_csv(out_dir / "ess_scale_comparison.tsv", sep="\t", index=False)
    write_report(out_dir, kappa, exposure, tails, ess)
    print(f"[review] wrote {out_dir}")


if __name__ == "__main__":
    main()