"""Analyse the mini-genome calibration stress sweep.

Reads ``<out>/sweep_results.tsv`` plus the per-cell ``per_region.feather``
and optional ``per_read.feather`` files, and prints:

  * Sweep summary table (λ_ratio, π, π_soft, μ_R, σ_R by ss × gdna)
  * γ-distribution breakdown per block class
  * Per-read confusion matrix (truth vs block_kind)
  * Cell with the worst λ_ratio — drill-down by region

Usage
-----
    python -m scripts.calibration.stress_mini.analyze /tmp/mini_stress/sweep
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def summary(sweep: pd.DataFrame) -> None:
    cols = ["label", "ss_true", "gdna_fraction",
            "lambda_truth", "lambda_est", "lambda_ratio",
            "ss_est", "mu_R", "sigma_R",
            "mixing_pi", "mixing_pi_soft",
            "strand_used", "em_n_iter"]
    with pd.option_context("display.float_format", "{:.4g}".format,
                            "display.width", 220):
        print("=" * 80)
        print("SWEEP SUMMARY")
        print("=" * 80)
        print(sweep[cols].to_string(index=False))

    # Pivot table: λ_ratio grid
    pivot = sweep.pivot(index="ss_true", columns="gdna_fraction",
                        values="lambda_ratio")
    print("\nλ_est / λ_truth grid (NaN when λ_truth=0):")
    print(pivot.round(3).to_string())


def per_region_breakdown(cells_dir: Path) -> None:
    """Aggregate per-region statistics across all cells."""
    rows = []
    for pr_path in sorted(cells_dir.glob("*/per_region.feather")):
        label = pr_path.parent.name
        pr = pd.read_feather(pr_path)
        agg = (pr.groupby("block_kind", dropna=False)
                 .agg(n=("gamma", "size"),
                      sum_len=("length", "sum"),
                      sum_E=("mappable_effective_length", "sum"),
                      sum_ntot=("n_total", "sum"),
                      sum_ngdna_true=("n_gdna_true", "sum"),
                      sum_Egdna_est=("E_gdna", "sum"),
                      mean_gamma=("gamma", "mean"),
                      frac_gamma_gt_0p5=("gamma",
                          lambda g: float(np.mean(np.asarray(g, dtype=float) > 0.5))))
                 .reset_index())
        agg.insert(0, "label", label)
        rows.append(agg)
    if not rows:
        return
    df = pd.concat(rows, ignore_index=True)
    df["gdna_recall"] = np.where(
        df["sum_ngdna_true"] > 0,
        df["sum_Egdna_est"] / df["sum_ngdna_true"],
        np.nan,
    )
    with pd.option_context("display.float_format", "{:.4g}".format,
                            "display.width", 220,
                            "display.max_rows", None):
        print("\n" + "=" * 80)
        print("PER-REGION by block_kind (aggregated within each cell)")
        print("=" * 80)
        print(df.to_string(index=False))


def worst_cell_drilldown(sweep: pd.DataFrame, cells_dir: Path) -> None:
    """Show per-region details for the cell with largest |λ_ratio - 1|."""
    s = sweep.dropna(subset=["lambda_ratio"]).copy()
    if s.empty:
        return
    s["dev"] = np.abs(s["lambda_ratio"] - 1.0)
    worst = s.sort_values("dev", ascending=False).iloc[0]
    label = worst["label"]
    print("\n" + "=" * 80)
    print(f"WORST CELL: {label} — λ_ratio = {worst.lambda_ratio:.3f}")
    print("=" * 80)
    pr_path = cells_dir / label / "per_region.feather"
    if not pr_path.exists():
        print(f"  per_region.feather missing at {pr_path}")
        return
    pr = pd.read_feather(pr_path)
    # Top regions by contribution to total E_gdna
    top = (pr.assign(contrib=pr["E_gdna"])
             .sort_values("contrib", ascending=False)
             .head(15))
    cols = ["region_id", "start", "end", "length", "mappable_effective_length",
            "block_kind", "n_total", "n_mrna_true", "n_nrna_true",
            "n_gdna_true", "gamma", "E_gdna"]
    avail = [c for c in cols if c in top.columns]
    with pd.option_context("display.float_format", "{:.4g}".format,
                            "display.width", 220):
        print("Top 15 regions by E_gdna:")
        print(top[avail].to_string(index=False))

    # γ histogram summary
    g = pr["gamma"].to_numpy()
    g_fin = g[np.isfinite(g)]
    bins = [0, 0.01, 0.1, 0.5, 0.9, 0.99, 1.0 + 1e-9]
    labels = ["<0.01", "0.01-0.1", "0.1-0.5", "0.5-0.9", "0.9-0.99", "≥0.99"]
    counts, _ = np.histogram(g_fin, bins=bins)
    print("\nγ distribution:")
    for lbl, c in zip(labels, counts):
        print(f"  {lbl:>10}: {c:5d}")


def per_read_confusion(cells_dir: Path) -> None:
    """Confusion matrix: truth × block_kind across all cells (if present)."""
    all_tabs = []
    for p in sorted(cells_dir.glob("*/per_read.feather")):
        df = pd.read_feather(p)
        label = p.parent.name
        tab = pd.crosstab(df["truth"], df["block_kind"], margins=True)
        tab.insert(0, "cell", label)
        all_tabs.append(tab)
    if not all_tabs:
        return
    print("\n" + "=" * 80)
    print("PER-READ CONFUSION (truth × block_kind) per cell")
    print("=" * 80)
    for t in all_tabs:
        lbl = t.iloc[0, 0]
        print(f"\n--- {lbl} ---")
        print(t.drop(columns=["cell"]).to_string())


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("sweep_dir", type=Path,
                   help="Directory passed to run_sweep.py --out")
    p.add_argument("--no-per-read", action="store_true",
                   help="Skip the per-read confusion matrices")
    args = p.parse_args(argv)

    sweep = pd.read_csv(args.sweep_dir / "sweep_results.tsv", sep="\t")
    cells_dir = args.sweep_dir / "cells"

    summary(sweep)
    per_region_breakdown(cells_dir)
    worst_cell_drilldown(sweep, cells_dir)
    if not args.no_per_read:
        per_read_confusion(cells_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
