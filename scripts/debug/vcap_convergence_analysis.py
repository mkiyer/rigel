#!/usr/bin/env python3
"""VCaP VBEM convergence deep-analysis.

Runs rigel quant via the Python API with ``emit_locus_stats=True`` to
capture per-locus EM convergence information, then produces a detailed
convergence report.

Usage:
    conda activate rigel
    python scripts/debug/vcap_convergence_analysis.py \
        --bam /path/to/aligned.bam \
        --index /path/to/rigel_index \
        --em-mode vbem \
        --threads 8 \
        -o results/vcap_convergence/
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def run_with_locus_stats(
    bam_path: str,
    index_dir: str,
    em_mode: str,
    threads: int,
    output_dir: str,
    seed: int = 42,
) -> tuple:
    """Run rigel pipeline with emit_locus_stats=True."""
    from rigel.config import EMConfig, PipelineConfig, BamScanConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import run_pipeline

    logging.info(f"Loading index from {index_dir}")
    index = TranscriptIndex.load(index_dir)
    logging.info(f"Loaded: {index.num_transcripts} transcripts")

    em_config = EMConfig(
        mode=em_mode,
        seed=seed,
    )
    scan_config = BamScanConfig(
        n_scan_threads=threads,
        sj_strand_tag="auto",
    )
    config = PipelineConfig(
        em=em_config,
        scan=scan_config,
        emit_locus_stats=True,
    )

    logging.info(f"Running pipeline with em_mode={em_mode}, threads={threads}")
    t0 = time.monotonic()
    result = run_pipeline(bam_path, index, config=config)
    elapsed = time.monotonic() - t0
    logging.info(f"Pipeline completed in {elapsed:.1f}s")

    return result, elapsed


def analyze_locus_stats(
    locus_stats: list[dict],
    em_mode: str,
    output_dir: Path,
) -> pd.DataFrame:
    """Analyze per-locus convergence statistics."""
    if not locus_stats:
        logging.warning("No locus stats available")
        return pd.DataFrame()

    df = pd.DataFrame(locus_stats)
    logging.info(f"Analyzing {len(df)} loci")

    # Compute derived metrics
    df["squarem_time_s"] = df["squarem_us"] / 1e6
    df["total_time_s"] = df["total_us"] / 1e6
    df["em_fraction"] = df["squarem_us"] / df["total_us"].clip(lower=1)

    # Determine max iterations (budget / 3 for SQUAREM)
    # Default: max_iterations=1000, so max_sq_iters = 333
    max_sq_iters = 333
    df["hit_max_iters"] = df["squarem_iterations"] >= max_sq_iters
    df["converged"] = ~df["hit_max_iters"]

    # Save full locus stats
    df.to_csv(output_dir / f"locus_stats_{em_mode}.csv", index=False)

    # Summary statistics
    n_loci = len(df)
    n_converged = df["converged"].sum()
    n_hit_max = df["hit_max_iters"].sum()

    summary = {
        "em_mode": em_mode,
        "n_loci": n_loci,
        "n_converged": int(n_converged),
        "n_hit_max_iters": int(n_hit_max),
        "convergence_rate": float(n_converged / n_loci) if n_loci > 0 else 0,
        "median_squarem_iters": float(df["squarem_iterations"].median()),
        "mean_squarem_iters": float(df["squarem_iterations"].mean()),
        "max_squarem_iters": int(df["squarem_iterations"].max()),
        "p90_squarem_iters": float(df["squarem_iterations"].quantile(0.9)),
        "p95_squarem_iters": float(df["squarem_iterations"].quantile(0.95)),
        "p99_squarem_iters": float(df["squarem_iterations"].quantile(0.99)),
        "total_squarem_time_s": float(df["squarem_time_s"].sum()),
        "mean_squarem_time_s": float(df["squarem_time_s"].mean()),
        "max_squarem_time_s": float(df["squarem_time_s"].max()),
        "total_n_components": int(df["n_components"].sum()),
        "max_n_components": int(df["n_components"].max()),
        "n_mega_loci": int(df["is_mega_locus"].sum()),
    }

    # Mega-locus specific stats
    mega = df[df["is_mega_locus"]]
    if len(mega) > 0:
        summary["mega_locus_details"] = []
        for _, row in mega.iterrows():
            summary["mega_locus_details"].append({
                "locus_idx": int(row["locus_idx"]),
                "n_transcripts": int(row["n_transcripts"]),
                "n_units": int(row["n_units"]),
                "n_components": int(row["n_components"]),
                "n_equiv_classes": int(row["n_equiv_classes"]),
                "ec_total_elements": int(row["ec_total_elements"]),
                "squarem_iterations": int(row["squarem_iterations"]),
                "squarem_time_s": float(row["squarem_time_s"]),
                "total_time_s": float(row["total_time_s"]),
                "hit_max_iters": bool(row["hit_max_iters"]),
                "estep_threads_used": int(row["estep_threads_used"]),
            })

    # Distribution of iterations
    iter_bins = [0, 1, 5, 10, 20, 50, 100, 200, 333, 334]
    iter_labels = ["0", "1-4", "5-9", "10-19", "20-49", "50-99",
                   "100-199", "200-333", ">333"]
    df["iter_bin"] = pd.cut(df["squarem_iterations"], bins=iter_bins,
                            labels=iter_labels, right=False)
    iter_dist = df["iter_bin"].value_counts().sort_index().to_dict()
    summary["iteration_distribution"] = {str(k): int(v) for k, v in iter_dist.items()}

    # Top 20 most expensive loci
    top_expensive = df.nlargest(20, "squarem_time_s")
    summary["top_20_expensive_loci"] = []
    for _, row in top_expensive.iterrows():
        summary["top_20_expensive_loci"].append({
            "locus_idx": int(row["locus_idx"]),
            "n_transcripts": int(row["n_transcripts"]),
            "n_components": int(row["n_components"]),
            "n_equiv_classes": int(row["n_equiv_classes"]),
            "squarem_iterations": int(row["squarem_iterations"]),
            "squarem_time_s": float(row["squarem_time_s"]),
            "hit_max_iters": bool(row["hit_max_iters"]),
        })

    # Correlation between locus size and iterations
    if len(df) > 10:
        from numpy import corrcoef
        valid = df[df["n_components"] > 1]
        if len(valid) > 10:
            corr = corrcoef(
                np.log1p(valid["n_components"].values),
                valid["squarem_iterations"].values,
            )[0, 1]
            summary["log_components_vs_iters_corr"] = float(corr)

    with open(output_dir / f"convergence_summary_{em_mode}.json", "w") as f:
        json.dump(summary, f, indent=2)

    logging.info(f"Convergence summary for {em_mode}:")
    logging.info(f"  Loci: {n_loci}, converged: {n_converged} "
                 f"({100*n_converged/n_loci:.1f}%), hit max: {n_hit_max}")
    logging.info(f"  Median iters: {summary['median_squarem_iters']:.0f}, "
                 f"P99: {summary['p99_squarem_iters']:.0f}")
    logging.info(f"  Total SQUAREM time: {summary['total_squarem_time_s']:.1f}s")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="VCaP VBEM convergence deep-analysis"
    )
    parser.add_argument("--bam", required=True, help="Input BAM path")
    parser.add_argument("--index", required=True, help="Rigel index directory")
    parser.add_argument("--em-mode", default="vbem", choices=["vbem", "map"],
                        help="EM mode (default: vbem)")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for convergence analysis")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    result, elapsed = run_with_locus_stats(
        bam_path=args.bam,
        index_dir=args.index,
        em_mode=args.em_mode,
        threads=args.threads,
        output_dir=str(output_dir),
        seed=args.seed,
    )

    # Extract locus stats
    locus_stats = result.estimator.locus_stats
    if locus_stats is None:
        logging.error("No locus stats captured - emit_locus_stats may be False")
        sys.exit(1)

    logging.info(f"Captured {len(locus_stats)} locus stat entries")

    # Analyze convergence
    df = analyze_locus_stats(locus_stats, args.em_mode, output_dir)

    # Save pipeline timing
    timing = {
        "em_mode": args.em_mode,
        "total_pipeline_time_s": elapsed,
        "bam": args.bam,
    }
    with open(output_dir / f"timing_{args.em_mode}.json", "w") as f:
        json.dump(timing, f, indent=2)

    logging.info(f"Results written to {output_dir}")


if __name__ == "__main__":
    main()
