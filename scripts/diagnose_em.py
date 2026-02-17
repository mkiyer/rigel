#!/usr/bin/env python3
"""
Diagnostic script:  re-run hulkrna on a single benchmark region
with varying settings to isolate the root cause of EM misallocation.

Tests:
  1. Default settings (gdna_threshold=0.5, vbem_prior=0.01)
  2. No gDNA shadows (gdna_threshold=0.0)
  3. No sparsity (vbem_prior=1.0, making VBEM ≈ standard EM)
  4. No gDNA + no sparsity (gdna_threshold=0.0, vbem_prior=1.0)
  5. Uniform effective lengths (disable eff_len normalization)

Usage:
  python scripts/diagnose_em.py \
      --region HBB_cluster \
      --bench-dir /Users/mkiyer/Downloads/hulkrna_runs/bench_zero_gdna_v2
"""

import argparse
import csv
import logging
import sys
from pathlib import Path

import numpy as np

# Ensure hulkrna is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hulkrna.index import HulkIndex
from hulkrna.pipeline import run_pipeline

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
)
logger = logging.getLogger(__name__)


def load_truth(per_tx_csv: Path) -> dict[str, int]:
    """Load truth counts from per_transcript_counts.csv."""
    truth = {}
    with open(per_tx_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            truth[row["transcript_id"]] = int(row["truth"])
    return truth


def run_hulkrna_variant(
    bam_path: Path,
    index: HulkIndex,
    *,
    label: str,
    gdna_threshold: float = 0.5,
    vbem_prior: float = 0.01,
    uniform_eff_len: bool = False,
    seed: int = 42,
) -> dict[str, float]:
    """Run hulkrna with custom settings and return transcript counts."""
    logger.info(f"\n{'='*60}")
    logger.info(f"Running: {label}")
    logger.info(f"  gdna_threshold={gdna_threshold}, vbem_prior={vbem_prior}, "
                f"uniform_eff_len={uniform_eff_len}")
    logger.info(f"{'='*60}")

    pipe = run_pipeline(
        bam_path,
        index,
        seed=seed,
        sj_strand_tag="auto",
        gdna_threshold=gdna_threshold,
        em_pseudocount=vbem_prior,
    )

    counts_df = pipe.counter.get_counts_df(index)
    counts = {
        row.transcript_id: float(row.count)
        for row in counts_df.itertuples(index=False)
    }

    # Report stats
    stats = pipe.stats
    gdna_rate = pipe.counter.gdna_contamination_rate
    gdna_total = pipe.counter.gdna_total
    logger.info(f"  gDNA: total={gdna_total:.1f}, rate={gdna_rate:.2%}")
    logger.info(f"  Stats: frags={stats.n_fragments}, "
                f"chimeric={stats.n_chimeric}, intergenic={stats.n_intergenic}, "
                f"unique_gene={stats.n_unique_gene}, multi_gene={stats.n_multi_gene}, "
                f"gdna_em={pipe.counter.gdna_em_count:.0f}")

    return counts


def compute_mae(truth: dict, predicted: dict) -> float:
    """Compute MAE across all truth transcripts."""
    errors = []
    for tid, t in truth.items():
        p = predicted.get(tid, 0.0)
        errors.append(abs(t - p))
    return np.mean(errors) if errors else 0.0


def compute_pearson(truth: dict, predicted: dict) -> float:
    """Compute Pearson correlation."""
    tids = sorted(truth.keys())
    t = np.array([truth.get(tid, 0) for tid in tids], dtype=np.float64)
    p = np.array([predicted.get(tid, 0.0) for tid in tids], dtype=np.float64)
    if t.std() == 0 or p.std() == 0:
        return 0.0
    return float(np.corrcoef(t, p)[0, 1])


def main():
    parser = argparse.ArgumentParser(description="Diagnose EM misallocation")
    parser.add_argument("--region", required=True, help="Region name (e.g., HBB_cluster)")
    parser.add_argument("--bench-dir", required=True, type=Path,
                        help="Benchmark output directory")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    reg_dir = args.bench_dir / args.region
    if not reg_dir.exists():
        logger.error(f"Region directory not found: {reg_dir}")
        sys.exit(1)

    bam_path = reg_dir / "align" / "reads.bam"
    index_dir = reg_dir / "hulkrna_index"
    per_tx_csv = reg_dir / "per_transcript_counts.csv"

    if not bam_path.exists():
        logger.error(f"BAM not found: {bam_path}")
        sys.exit(1)

    truth = load_truth(per_tx_csv)
    index = HulkIndex.load(index_dir)

    # Load salmon/kallisto from CSV for comparison
    salmon_counts = {}
    kallisto_counts = {}
    with open(per_tx_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            salmon_counts[row["transcript_id"]] = float(row["salmon"])
            kallisto_counts[row["transcript_id"]] = float(row["kallisto"])

    # Define test configurations
    configs = [
        ("1_default", dict(gdna_threshold=0.5, vbem_prior=0.01)),
        ("2_no_gdna", dict(gdna_threshold=0.0, vbem_prior=0.01)),
        ("3_no_sparsity", dict(gdna_threshold=0.5, vbem_prior=1.0)),
        ("4_no_gdna_no_sparsity", dict(gdna_threshold=0.0, vbem_prior=1.0)),
        ("5_large_prior", dict(gdna_threshold=0.0, vbem_prior=0.1)),
    ]

    results = {}
    for label, kwargs in configs:
        counts = run_hulkrna_variant(
            bam_path, index, label=label, seed=args.seed, **kwargs
        )
        results[label] = counts

    # Print comparison table
    print("\n" + "="*100)
    print(f"DIAGNOSTIC RESULTS: {args.region}")
    print("="*100)

    # Summary metrics
    print(f"\n{'Config':<30} {'MAE':>10} {'Pearson':>10} {'Total':>12} {'Total_err':>12}")
    print("-"*74)

    for label, counts in results.items():
        mae = compute_mae(truth, counts)
        pearson = compute_pearson(truth, counts)
        total = sum(counts.get(tid, 0) for tid in truth)
        total_err = abs(sum(truth.values()) - total)
        print(f"{label:<30} {mae:>10.2f} {pearson:>10.5f} {total:>12.1f} {total_err:>12.1f}")

    # Reference: salmon and kallisto
    print("-"*74)
    mae_s = compute_mae(truth, salmon_counts)
    pearson_s = compute_pearson(truth, salmon_counts)
    total_s = sum(salmon_counts.get(tid, 0) for tid in truth)
    print(f"{'salmon':<30} {mae_s:>10.2f} {pearson_s:>10.5f} {total_s:>12.1f}")

    mae_k = compute_mae(truth, kallisto_counts)
    pearson_k = compute_pearson(truth, kallisto_counts)
    total_k = sum(kallisto_counts.get(tid, 0) for tid in truth)
    print(f"{'kallisto':<30} {mae_k:>10.2f} {pearson_k:>10.5f} {total_k:>12.1f}")

    # Per-transcript detail for worst misallocations
    print(f"\n{'='*130}")
    print("Per-transcript detail (sorted by default error)")
    print(f"{'='*130}")
    print(f"{'transcript_id':<25} {'truth':>8} {'default':>10} {'no_gdna':>10} "
          f"{'no_sparse':>10} {'no_both':>10} {'salmon':>10} {'kallisto':>10}")
    print("-"*130)

    # Sort by absolute error in default config
    tids_by_err = sorted(truth.keys(),
                         key=lambda t: abs(truth[t] - results["1_default"].get(t, 0)),
                         reverse=True)

    for tid in tids_by_err:
        t = truth[tid]
        d = results["1_default"].get(tid, 0)
        ng = results["2_no_gdna"].get(tid, 0)
        ns = results["3_no_sparsity"].get(tid, 0)
        nb = results["4_no_gdna_no_sparsity"].get(tid, 0)
        s = salmon_counts.get(tid, 0)
        k = kallisto_counts.get(tid, 0)
        print(f"{tid:<25} {t:>8} {d:>10.1f} {ng:>10.1f} {ns:>10.1f} "
              f"{nb:>10.1f} {s:>10.1f} {k:>10.1f}")


if __name__ == "__main__":
    main()
