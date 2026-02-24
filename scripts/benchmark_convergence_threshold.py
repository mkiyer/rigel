#!/usr/bin/env python3
"""Benchmark EM convergence threshold impact on accuracy and speed.

Runs the hulkrna pipeline on all 10 pristine benchmark regions at
multiple convergence thresholds (1e-6, 1e-5, 1e-4), comparing
per-transcript counts against ground truth and against the tightest
threshold (1e-6) as baseline.

Usage
-----
    PYTHONPATH=src python scripts/benchmark_convergence_threshold.py \
        --bench-dir /path/to/bench_pristine_10_regions
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# -- Make hulkrna importable when run from repo root --
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import hulkrna.estimator as estimator_mod
from hulkrna.index import HulkIndex
from hulkrna.pipeline import run_pipeline


def _discover_regions(bench_dir: Path) -> list[Path]:
    """Find all region directories with an oracle BAM + hulkrna_index."""
    regions = []
    for d in sorted(bench_dir.iterdir()):
        if not d.is_dir():
            continue
        idx_dir = d / "hulkrna_index"
        cond_dir = d / "gdna_none_nrna_none_ss_1.00"
        bam = cond_dir / "align_oracle" / "reads_namesort.bam"
        truth = cond_dir / "per_transcript_counts.csv"
        if idx_dir.is_dir() and bam.is_file() and truth.is_file():
            regions.append(d)
    return regions


def _load_truth(region_dir: Path) -> pd.DataFrame:
    """Load per-transcript truth counts."""
    csv_path = (
        region_dir / "gdna_none_nrna_none_ss_1.00" / "per_transcript_counts.csv"
    )
    df = pd.read_csv(csv_path)
    return df[["transcript_id", "truth"]].copy()


def _run_at_threshold(
    index: HulkIndex,
    bam_path: Path,
    threshold: float,
) -> tuple[dict[str, float], float, int]:
    """Run pipeline with a given EM convergence threshold.

    Returns (transcript_id -> count dict, elapsed_sec, n_sq_iters).
    """
    # Monkey-patch the module-level convergence threshold.
    orig = estimator_mod._EM_CONVERGENCE_DELTA
    estimator_mod._EM_CONVERGENCE_DELTA = threshold
    try:
        t0 = time.monotonic()
        pipe = run_pipeline(bam_path, index, include_multimap=True)
        elapsed = time.monotonic() - t0
    finally:
        estimator_mod._EM_CONVERGENCE_DELTA = orig

    counts_df = pipe.estimator.get_counts_df(index)
    counts = {
        row.transcript_id: float(row.count)
        for row in counts_df.itertuples(index=False)
    }
    return counts, elapsed


def _accuracy_metrics(
    pred: dict[str, float],
    truth: dict[str, float],
) -> dict[str, float]:
    """Compute accuracy metrics between prediction and truth."""
    all_keys = sorted(set(pred) | set(truth))
    y_pred = np.array([pred.get(k, 0.0) for k in all_keys])
    y_true = np.array([truth.get(k, 0.0) for k in all_keys])

    diff = y_pred - y_true
    mae = float(np.abs(diff).mean())
    rmse = float(np.sqrt((diff ** 2).mean()))
    max_abs = float(np.abs(diff).max())

    # Pearson correlation
    if np.std(y_pred) > 0 and np.std(y_true) > 0:
        corr = float(np.corrcoef(y_pred, y_true)[0, 1])
    else:
        corr = 1.0

    # Mean absolute percentage error (on transcripts with truth > 0)
    mask = y_true > 0
    if mask.any():
        mape = float(np.abs(diff[mask] / y_true[mask]).mean()) * 100
    else:
        mape = 0.0

    # Sum accuracy (total count preservation)
    sum_pred = float(y_pred.sum())
    sum_true = float(y_true.sum())
    sum_diff = sum_pred - sum_true

    return {
        "mae": mae,
        "rmse": rmse,
        "max_abs_err": max_abs,
        "pearson_r": corr,
        "mape_pct": mape,
        "sum_pred": sum_pred,
        "sum_true": sum_true,
        "sum_diff": sum_diff,
    }


def _pairwise_metrics(
    pred_a: dict[str, float],
    pred_b: dict[str, float],
) -> dict[str, float]:
    """Compute difference metrics between two predictions."""
    all_keys = sorted(set(pred_a) | set(pred_b))
    a = np.array([pred_a.get(k, 0.0) for k in all_keys])
    b = np.array([pred_b.get(k, 0.0) for k in all_keys])

    diff = a - b
    mae = float(np.abs(diff).mean())
    max_abs = float(np.abs(diff).max())
    rmse = float(np.sqrt((diff ** 2).mean()))

    # Relative difference on transcripts with nonzero baseline
    mask = b > 0
    if mask.any():
        rel_diff = float(np.abs(diff[mask] / b[mask]).mean()) * 100
    else:
        rel_diff = 0.0

    return {
        "mae_vs_baseline": mae,
        "max_abs_vs_baseline": max_abs,
        "rmse_vs_baseline": rmse,
        "mean_rel_diff_pct": rel_diff,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bench-dir",
        type=Path,
        required=True,
        help="Path to bench_pristine_10_regions directory.",
    )
    parser.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=[1e-6, 1e-5, 1e-4],
        help="Convergence thresholds to test (default: 1e-6 1e-5 1e-4).",
    )
    args = parser.parse_args()

    regions = _discover_regions(args.bench_dir)
    if not regions:
        print(f"ERROR: No valid regions found in {args.bench_dir}", file=sys.stderr)
        sys.exit(1)

    thresholds = sorted(args.thresholds)
    baseline_thr = thresholds[0]  # tightest threshold as baseline

    print(f"Regions: {len(regions)}")
    print(f"Thresholds: {thresholds}")
    print(f"Baseline: {baseline_thr}")
    print()

    # Collect all results
    all_rows = []

    for region_dir in regions:
        region_name = region_dir.name
        idx_dir = region_dir / "hulkrna_index"
        bam_path = (
            region_dir / "gdna_none_nrna_none_ss_1.00"
            / "align_oracle" / "reads_namesort.bam"
        )

        print(f"{'='*70}")
        print(f"Region: {region_name}")
        print(f"{'='*70}")

        index = HulkIndex.load(idx_dir)
        truth_df = _load_truth(region_dir)
        truth_dict = dict(zip(truth_df.transcript_id, truth_df.truth))

        results_by_thr: dict[float, dict[str, float]] = {}
        timings: dict[float, float] = {}

        for thr in thresholds:
            print(f"  threshold={thr:.0e} ... ", end="", flush=True)
            counts, elapsed = _run_at_threshold(index, bam_path, thr)
            results_by_thr[thr] = counts
            timings[thr] = elapsed
            print(f"{elapsed:.2f}s")

        # Compute metrics for each threshold
        for thr in thresholds:
            counts = results_by_thr[thr]
            vs_truth = _accuracy_metrics(counts, truth_dict)

            row = {
                "region": region_name,
                "threshold": thr,
                "elapsed_s": timings[thr],
                **{f"truth_{k}": v for k, v in vs_truth.items()},
            }

            # Pairwise vs baseline (skip for baseline itself)
            if thr != baseline_thr:
                vs_base = _pairwise_metrics(counts, results_by_thr[baseline_thr])
                row.update(vs_base)
            else:
                row.update({
                    "mae_vs_baseline": 0.0,
                    "max_abs_vs_baseline": 0.0,
                    "rmse_vs_baseline": 0.0,
                    "mean_rel_diff_pct": 0.0,
                })

            all_rows.append(row)

    # Build results DataFrame
    df = pd.DataFrame(all_rows)

    # ── Summary tables ──
    print()
    print("=" * 80)
    print("RESULTS: Accuracy vs Ground Truth (per threshold)")
    print("=" * 80)

    for thr in thresholds:
        sub = df[df.threshold == thr]
        print(f"\n  Threshold: {thr:.0e}")
        print(f"  {'Region':<35} {'MAE':>8} {'RMSE':>8} {'MaxErr':>8} "
              f"{'Pearson':>8} {'MAPE%':>8} {'Time(s)':>8}")
        print(f"  {'-'*35} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
        for _, r in sub.iterrows():
            print(f"  {r.region:<35} {r.truth_mae:>8.2f} {r.truth_rmse:>8.2f} "
                  f"{r.truth_max_abs_err:>8.1f} {r.truth_pearson_r:>8.6f} "
                  f"{r.truth_mape_pct:>8.2f} {r.elapsed_s:>8.2f}")
        print(f"  {'MEAN':<35} {sub.truth_mae.mean():>8.2f} "
              f"{sub.truth_rmse.mean():>8.2f} {sub.truth_max_abs_err.mean():>8.1f} "
              f"{sub.truth_pearson_r.mean():>8.6f} "
              f"{sub.truth_mape_pct.mean():>8.2f} {sub.elapsed_s.mean():>8.2f}")

    print()
    print("=" * 80)
    print("RESULTS: Deviation from Baseline (threshold={:.0e})".format(baseline_thr))
    print("=" * 80)

    for thr in thresholds:
        if thr == baseline_thr:
            continue
        sub = df[df.threshold == thr]
        print(f"\n  Threshold: {thr:.0e} vs baseline {baseline_thr:.0e}")
        print(f"  {'Region':<35} {'MAE':>10} {'MaxAbs':>10} {'RMSE':>10} "
              f"{'RelDiff%':>10} {'dTime(s)':>10}")
        print(f"  {'-'*35} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
        base_sub = df[df.threshold == baseline_thr].set_index("region")
        for _, r in sub.iterrows():
            dt = r.elapsed_s - base_sub.loc[r.region, "elapsed_s"]
            print(f"  {r.region:<35} {r.mae_vs_baseline:>10.4f} "
                  f"{r.max_abs_vs_baseline:>10.4f} {r.rmse_vs_baseline:>10.4f} "
                  f"{r.mean_rel_diff_pct:>10.4f} {dt:>+10.2f}")
        print(f"  {'MEAN':<35} {sub.mae_vs_baseline.mean():>10.4f} "
              f"{sub.max_abs_vs_baseline.mean():>10.4f} "
              f"{sub.rmse_vs_baseline.mean():>10.4f} "
              f"{sub.mean_rel_diff_pct.mean():>10.4f} "
              f"{(sub.elapsed_s.mean() - df[df.threshold==baseline_thr].elapsed_s.mean()):>+10.2f}")

    # ── Timing summary ──
    print()
    print("=" * 80)
    print("TIMING SUMMARY")
    print("=" * 80)
    for thr in thresholds:
        sub = df[df.threshold == thr]
        print(f"  threshold={thr:.0e}: mean={sub.elapsed_s.mean():.2f}s, "
              f"max={sub.elapsed_s.max():.2f}s, "
              f"total={sub.elapsed_s.sum():.1f}s")

    # Save CSV
    csv_path = args.bench_dir / "convergence_threshold_benchmark.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nDetailed results saved to: {csv_path}")


if __name__ == "__main__":
    main()
