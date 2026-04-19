"""Summarise the calibration v4 sweep results.

Reads ``results.tsv`` produced by ``calibration_sweep.py`` and writes:

* ``summary_by_ss_gdna.tsv`` — mean / std / median λ̂ relative-error,
  grouped by (ss, gdna_fraction).  Reports both the pipeline-realistic
  and the truth-SS-injected calibrations.
* ``per_pathway.tsv`` — density / strand pathway accuracy in isolation.
* ``calibration_v4_summary.md`` — markdown table snippets ready to
  paste into ``docs/benchmarks/calibration_v4_baseline.md``.

Usage::

    python scripts/benchmark/analyze_calibration_sweep.py \\
        -i scripts/benchmark/results/calibration_v4
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _agg(g: pd.DataFrame, col: str) -> dict:
    s = g[col].dropna()
    if len(s) == 0:
        return {f"{col}_mean": np.nan, f"{col}_std": np.nan,
                f"{col}_median": np.nan, f"{col}_n": 0}
    return {
        f"{col}_mean": float(s.mean()),
        f"{col}_std": float(s.std(ddof=0)),
        f"{col}_median": float(s.median()),
        f"{col}_n": int(len(s)),
    }


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """Group by (ss, gdna_fraction) and compute summary stats."""
    rows = []
    for (ss, gf), g in df.groupby(["ss", "gdna_fraction"]):
        row = {"ss": ss, "gdna_fraction": gf, "n_replicates": len(g)}
        for col in [
            "lambda_relative_error",
            "lambda_log2_ratio",
            "gdna_em_relative_error",
            "truthss_lambda_relative_error",
            "truthss_strand_relative_error",
            "consistency_chi2",
            "truthss_consistency_chi2",
        ]:
            if col in g.columns:
                row.update(_agg(g, col))
        # Truth side
        row["lambda_truth_mean"] = float(g["lambda_truth"].mean())
        row["n_gdna_truth_mean"] = float(g["n_gdna_truth"].mean())
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["ss", "gdna_fraction"])


def md_table_pipeline(summary: pd.DataFrame) -> str:
    """Pipeline-realistic λ̂ recovery table (uses estimated SS)."""
    lines = [
        "| SS | gDNA fraction | n | truth λ_G | λ̂_pool / λ_G | rel. err. (median, [P25, P75]) |",
        "|----|---------------|---|-----------|---------------|--------------------------------|",
    ]
    for _, r in summary.iterrows():
        rel = r["lambda_relative_error_median"]
        ratio = (1.0 + rel) if not np.isnan(rel) else np.nan
        lines.append(
            f"| {r['ss']:.2f} | {r['gdna_fraction']:.2f} | "
            f"{int(r['n_replicates'])} | {r['lambda_truth_mean']:.3e} | "
            f"{ratio:.3f}× | "
            f"{rel:+.1%} ± {r['lambda_relative_error_std']:.1%} |"
        )
    return "\n".join(lines)


def md_table_truth_ss(summary: pd.DataFrame) -> str:
    """Strand-pathway-isolated table (truth SS injected)."""
    lines = [
        "| SS | gDNA fraction | n | truthSS λ̂_pool err | strand-only err |",
        "|----|---------------|---|---------------------|-----------------|",
    ]
    for _, r in summary.iterrows():
        if r["ss"] <= 0.5 + 1e-6:
            continue  # strand pathway disabled at SS=0.5
        rel_p = r.get("truthss_lambda_relative_error_median", np.nan)
        rel_s = r.get("truthss_strand_relative_error_median", np.nan)
        std_p = r.get("truthss_lambda_relative_error_std", np.nan)
        std_s = r.get("truthss_strand_relative_error_std", np.nan)
        lines.append(
            f"| {r['ss']:.2f} | {r['gdna_fraction']:.2f} | "
            f"{int(r['n_replicates'])} | "
            f"{rel_p:+.2%} ± {std_p:.2%} | "
            f"{rel_s:+.2%} ± {std_s:.2%} |"
        )
    return "\n".join(lines)


def md_table_zero_gdna(df: pd.DataFrame) -> str:
    """Sanity check at gDNA=0: λ̂ should be ≪ λ̂ at the next-up gDNA level."""
    z = df[df["gdna_fraction"] == 0.0]
    if z.empty:
        return "_(no zero-gDNA points)_"
    lines = [
        "| SS | seed | λ̂_pool (pipeline) | λ̂_density | λ̂_strand (truth-SS) |",
        "|----|------|---------------------|-------------|----------------------|",
    ]
    for _, r in z.iterrows():
        lines.append(
            f"| {r['ss']:.2f} | {int(r['seed'])} | "
            f"{r.get('lambda_pool', np.nan):.3e} | "
            f"{r.get('lambda_density', np.nan):.3e} | "
            f"{r.get('truthss_lambda_strand', np.nan)!s:>20} |"
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--in", dest="indir", type=Path, required=True,
                        help="Sweep output directory (contains results.tsv).")
    parser.add_argument("-o", "--out", dest="outdir", type=Path, default=None,
                        help="Output directory (default: same as --in).")
    args = parser.parse_args()

    indir = args.indir
    outdir = args.outdir or indir
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(indir / "results.tsv", sep="\t")
    print(f"Loaded {len(df)} rows from {indir / 'results.tsv'}")
    if "ok" in df.columns:
        n_failed = int((~df["ok"]).sum())
        if n_failed > 0:
            print(f"  WARNING: {n_failed} pipelines failed; excluding")
            df = df[df["ok"]]

    summary = summarise(df)
    summary_path = outdir / "summary_by_ss_gdna.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"  wrote {summary_path}  ({len(summary)} rows)")

    md_path = outdir / "calibration_v4_summary.md"
    parts = [
        "## Calibration v4 — synthetic sweep summary\n",
        "### Pipeline-realistic λ̂ recovery\n",
        "*Pipeline runs with strand-model-estimated SS (which on these "
        "oracle BAMs collapses to 0.5, so this is effectively the "
        "**density-only** pathway).*\n",
        md_table_pipeline(summary),
        "\n\n### Truth-SS-injected calibration\n",
        "*Calibration re-run with the simulated SS injected, exercising "
        "the strand pathway in isolation.*\n",
        md_table_truth_ss(summary),
        "\n\n### Zero-gDNA sanity check\n",
        "*With no contamination, λ̂ should collapse toward zero "
        "(bounded only by noise / RNA-leakage).*\n",
        md_table_zero_gdna(df),
        "\n",
    ]
    md_path.write_text("\n".join(parts))
    print(f"  wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
