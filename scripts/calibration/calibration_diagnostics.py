#!/usr/bin/env python3
"""Calibration Diagnostics — rendered markdown + plots for gDNA deconvolution.

Runs the calibration EM on a BAM file and produces a diagnostic report with:

1. **Input Data Summary** — region counts, coverage, strand, and density distributions
2. **Initialization Quality** — seed partition, expressed vs gDNA seeds
3. **EM Convergence** — per-iteration parameter traces and histogram evolution
4. **Final Results** — posterior distributions, density partition, strand separation
5. **Ground Truth** (optional) — classification accuracy when simulated truth is available

Usage::

    python scripts/calibration_diagnostics.py \\
        --bam sample.bam \\
        --index index_dir/ \\
        -o diagnostics/ \\
        [--truth]              # parse ground truth from simulated read names
        [--ss 0.95]            # override strand specificity (auto-estimated if omitted)
        [--max-iter 50]
        [--convergence-tol 1e-4]

Output: ``diagnostics/calibration_report.md`` + ``diagnostics/plots/*.png``
"""
from __future__ import annotations

import argparse
import base64
import logging
import re
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Ensure src/ is importable when run from scripts/
# ---------------------------------------------------------------------------
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

from rigel.calibration import (
    CalibrationResult,
    _DENSITY_BINS,
    _STRAND_BINS,
    build_density_histogram,
    build_strand_histogram,
    calibrate_gdna,
    compute_log_density,
    compute_region_stats,
    compute_sense_fraction,
    estimate_kappa_sym,
)
from rigel.index import TranscriptIndex
from rigel.region_evidence import count_region_evidence

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Matplotlib setup (Agg backend for headless rendering)
# ---------------------------------------------------------------------------
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    "figure.dpi": 150,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
})


# ===================================================================
# Ground truth parsing (simulated BAMs)
# ===================================================================


def parse_ground_truth_labels(
    bam_path: str,
    region_counts: pd.DataFrame,
) -> np.ndarray | None:
    """Parse ground truth gDNA labels from simulated BAM read names.

    Returns boolean array (n_regions,) where True = region is purely gDNA,
    or None if truth can't be parsed.

    Convention: read names starting with 'gdna:' are gDNA fragments.
    """
    try:
        import pysam
    except ImportError:
        logger.warning("pysam not available — cannot parse ground truth")
        return None

    n_regions = len(region_counts)

    # Count gDNA vs RNA fragments per region by scanning BAM read names.
    # This is approximate: we count unique read names that map to each region.
    # For a more precise approach, we'd need to resolve fragments to regions.
    # Instead, we classify regions by whether they have ANY spliced reads (RNA)
    # or only unspliced reads from gDNA read names.

    # Simpler approach: a region is "truly gDNA" if it has zero spliced reads
    # and all its unspliced reads come from gDNA fragments.
    # We approximate this from the count data + a BAM scan of read names.

    gdna_frag_names: set[str] = set()
    rna_frag_names: set[str] = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            if qname.startswith("gdna:"):
                gdna_frag_names.add(qname)
            else:
                rna_frag_names.add(qname)

    total_gdna = len(gdna_frag_names)
    total_rna = len(rna_frag_names)

    if total_gdna == 0:
        logger.warning("No gDNA fragments found in BAM — truth labels unavailable")
        return None

    logger.info(
        f"Ground truth: {total_gdna:,} gDNA fragments, {total_rna:,} RNA fragments"
    )

    # For region-level truth: use spliced reads as definitive RNA marker,
    # then classify unspliced-only regions by their dominant read type.
    # This is coarse but sufficient for diagnostic purposes.
    has_spliced = (
        region_counts["n_spliced_pos"].values
        + region_counts["n_spliced_neg"].values
    ) > 0
    has_unspliced = (
        region_counts["n_unspliced_pos"].values
        + region_counts["n_unspliced_neg"].values
    ) > 0

    # Regions with spliced reads are definitively RNA (truth = not gDNA)
    # Regions with only unspliced reads: need per-region read name scan.
    # For now: approximate as "true gDNA" if no spliced and no coverage
    # or very low coverage (ambiguous). Mark as unknown.
    truth_gdna = np.full(n_regions, np.nan)
    truth_gdna[has_spliced] = 0.0  # Definitely has RNA
    truth_gdna[~has_unspliced & ~has_spliced] = np.nan  # No data

    # For unspliced-only regions, we can't resolve per-region without
    # the full fragment-to-region mapping. Return what we have.
    # The ground_truth_fraction gives the global gDNA fraction.
    gdna_frac = total_gdna / max(total_gdna + total_rna, 1)

    # Store as metadata for the report
    truth_meta = {
        "n_gdna_frags": total_gdna,
        "n_rna_frags": total_rna,
        "gdna_fraction": gdna_frac,
        "per_region": truth_gdna,
        "has_spliced": has_spliced,
    }
    return truth_meta


# ===================================================================
# Strand specificity estimation from spliced-region strand ratios
# ===================================================================


def estimate_strand_specificity(stats: dict[str, np.ndarray]) -> float:
    """Estimate strand specificity from regions with spliced reads.

    Spliced reads are definitively RNA. For + genes with spliced reads,
    RNA follows R1-antisense convention: p(+) = 1 - SS.

    Returns SS ∈ [0.5, 1.0].
    """
    n_spliced = stats["n_spliced"]
    n_pos = stats["n_pos"]
    n_unspliced = stats["n_unspliced"]
    gene_strand = stats["gene_strand"]

    # Use regions with substantial spliced reads and clear gene strand
    has_splice = n_spliced > 5
    has_unspliced = n_unspliced >= 10
    valid = has_splice & has_unspliced

    # For + genes: p(+) = 1 - SS → SS = 1 - p(+) = 1 - (n_pos / n_unspliced)
    plus_mask = valid & (gene_strand == 1)
    # For - genes: p(+) = SS → SS = n_pos / n_unspliced
    minus_mask = valid & (gene_strand == -1)

    ss_estimates = []
    if plus_mask.any():
        strand_ratio = n_pos[plus_mask] / n_unspliced[plus_mask]
        ss_estimates.extend((1.0 - strand_ratio).tolist())
    if minus_mask.any():
        strand_ratio = n_pos[minus_mask] / n_unspliced[minus_mask]
        ss_estimates.extend(strand_ratio.tolist())

    if not ss_estimates:
        logger.warning("Cannot estimate SS — no spliced regions with known gene strand")
        return 0.5

    ss = float(np.median(ss_estimates))
    ss = np.clip(ss, 0.5, 1.0)
    logger.info(f"Estimated strand specificity: {ss:.4f} (from {len(ss_estimates)} regions)")
    return ss


# ===================================================================
# Plot helpers
# ===================================================================


def _hist_bar(ax, bin_edges, density, color, label, alpha=0.7):
    """Plot a histogram as a bar chart from bin_edges and density arrays."""
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    width = bin_edges[1] - bin_edges[0]
    ax.bar(centers, density, width=width, color=color, alpha=alpha, label=label)


def _save_fig(fig, plot_dir: Path, name: str) -> str:
    """Save figure and return relative path for markdown embedding."""
    path = plot_dir / f"{name}.png"
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return f"plots/{name}.png"


# ===================================================================
# Section 1: Input Data Summary
# ===================================================================


def plot_input_summary(
    stats: dict[str, np.ndarray],
    eligible: np.ndarray,
    plot_dir: Path,
) -> list[str]:
    """Generate input data summary plots."""
    paths = []
    n_total = stats["n_total"]
    region_length = stats["region_length"]
    density = stats["density"]

    # --- Coverage distribution ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # Total fragment count distribution
    ax = axes[0]
    ax.hist(n_total[eligible], bins=50, color="steelblue", alpha=0.8, edgecolor="white")
    ax.set_xlabel("Total fragments per region")
    ax.set_ylabel("Count")
    ax.set_title("Fragment Count Distribution")
    ax.axvline(np.median(n_total[eligible]), color="red", ls="--", label=f"median={np.median(n_total[eligible]):.0f}")
    ax.legend(fontsize=7)

    # Region length distribution
    ax = axes[1]
    ax.hist(region_length[eligible], bins=50, color="seagreen", alpha=0.8, edgecolor="white")
    ax.set_xlabel("Region length (bp)")
    ax.set_ylabel("Count")
    ax.set_title("Region Length Distribution")

    # Density distribution (log scale)
    ax = axes[2]
    d_pos = density[eligible & (density > 0)]
    if len(d_pos) > 0:
        ax.hist(np.log10(d_pos), bins=50, color="darkorange", alpha=0.8, edgecolor="white")
    ax.set_xlabel("log₁₀(density)")
    ax.set_ylabel("Count")
    ax.set_title("Density Distribution (eligible)")

    fig.suptitle("Input Data Summary", fontsize=12, fontweight="bold")
    fig.tight_layout()
    paths.append(_save_fig(fig, plot_dir, "01_input_summary"))

    # --- Strand ratio vs splice rate scatter ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    ax = axes[0]
    sr = stats["strand_ratio"][eligible]
    spl = stats["splice_rate"][eligible]
    valid = np.isfinite(sr)
    sc = ax.scatter(sr[valid], spl[valid], c=np.log10(n_total[eligible][valid] + 1),
                    s=5, alpha=0.5, cmap="viridis")
    plt.colorbar(sc, ax=ax, label="log₁₀(n_total + 1)")
    ax.set_xlabel("Strand ratio (n_pos / n_unspliced)")
    ax.set_ylabel("Splice rate")
    ax.set_title("Strand Ratio vs Splice Rate")

    ax = axes[1]
    ax.hist(spl[spl > 0], bins=50, color="mediumpurple", alpha=0.8, edgecolor="white")
    ax.set_xlabel("Splice rate (among spliced regions)")
    ax.set_ylabel("Count")
    ax.set_title("Splice Rate Distribution")

    fig.tight_layout()
    paths.append(_save_fig(fig, plot_dir, "02_strand_vs_splice"))

    return paths


# ===================================================================
# Section 2: Empirical Histograms (pre-EM)
# ===================================================================


def plot_empirical_histograms(
    sense_frac: np.ndarray,
    log_d: np.ndarray,
    eligible: np.ndarray,
    stats: dict[str, np.ndarray],
    plot_dir: Path,
) -> list[str]:
    """Plot the raw empirical distributions that feed into the EM."""
    paths = []

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # --- Sense fraction histogram ---
    ax = axes[0]
    sf_valid = sense_frac[eligible & np.isfinite(sense_frac)]
    if len(sf_valid) > 0:
        ax.hist(sf_valid, bins=_STRAND_BINS, color="steelblue", alpha=0.8, edgecolor="white")
        ax.axvline(0.5, color="red", ls="--", lw=1, label="gDNA expected (0.5)")
        ax.axvline(np.median(sf_valid), color="green", ls="--", lw=1, label=f"median={np.median(sf_valid):.3f}")
        ax.legend(fontsize=7)
    ax.set_xlabel("Sense Fraction (gene-strand normalised)")
    ax.set_ylabel("Count")
    ax.set_title("Sense Fraction Distribution (eligible)")

    # --- Log density histogram ---
    ax = axes[1]
    ld = log_d[eligible]
    if len(ld) > 0:
        ax.hist(ld, bins=_DENSITY_BINS, color="darkorange", alpha=0.8, edgecolor="white")
        ax.axvline(np.median(ld), color="green", ls="--", lw=1, label=f"median={np.median(ld):.2f}")
        ax.legend(fontsize=7)
    ax.set_xlabel("log(n/L + ε)")
    ax.set_ylabel("Count")
    ax.set_title("Log Density Distribution (eligible)")

    fig.suptitle("Raw Empirical Distributions", fontsize=12, fontweight="bold")
    fig.tight_layout()
    paths.append(_save_fig(fig, plot_dir, "03_empirical_histograms"))

    # --- Sense fraction split by spliced/unspliced ---
    fig, ax = plt.subplots(figsize=(8, 4))
    has_splice = stats["n_spliced"][eligible] > 0
    sf_elig = sense_frac[eligible]
    sf_spliced = sf_elig[has_splice & np.isfinite(sf_elig)]
    sf_unspliced = sf_elig[~has_splice & np.isfinite(sf_elig)]
    if len(sf_spliced) > 0:
        ax.hist(sf_spliced, bins=_STRAND_BINS, alpha=0.6, color="green",
                label=f"Spliced (n={len(sf_spliced)})", density=True, edgecolor="white")
    if len(sf_unspliced) > 0:
        ax.hist(sf_unspliced, bins=_STRAND_BINS, alpha=0.6, color="gray",
                label=f"Unspliced-only (n={len(sf_unspliced)})", density=True, edgecolor="white")
    ax.set_xlabel("Sense Fraction")
    ax.set_ylabel("Density")
    ax.set_title("Sense Fraction: Spliced vs Unspliced Regions")
    ax.legend()
    fig.tight_layout()
    paths.append(_save_fig(fig, plot_dir, "04_strand_spliced_vs_unspliced"))

    return paths


# ===================================================================
# Section 3: Initialization Quality
# ===================================================================


def plot_initialization(
    result: CalibrationResult,
    log_d: np.ndarray,
    sense_frac: np.ndarray,
    eligible: np.ndarray,
    plot_dir: Path,
) -> list[str]:
    """Visualize the two-phase seed partition."""
    paths = []
    init_diag = result.init_diagnostics
    if init_diag is None:
        return paths

    # Get the first iteration's gamma as the post-init state
    if result.iteration_history and len(result.iteration_history) > 0:
        init_gamma = result.iteration_history[0].get("gamma")
    else:
        init_gamma = result.region_posteriors

    if init_gamma is None:
        init_gamma = result.region_posteriors

    stats = result.region_stats
    n_spliced = stats["n_spliced"]

    fig = plt.figure(figsize=(14, 8))
    gs = GridSpec(2, 3, figure=fig)

    # 1) Log density colored by seed label
    ax = fig.add_subplot(gs[0, 0])
    expressed_seed = eligible & (n_spliced > 0)
    gdna_seed = eligible & (init_gamma == 1.0) & ~expressed_seed
    ambiguous = eligible & ~expressed_seed & ~gdna_seed

    ld_expressed = log_d[expressed_seed]
    ld_gdna = log_d[gdna_seed]
    ld_ambig = log_d[ambiguous]

    for vals, label, color in [
        (ld_expressed, f"Expressed seed (n={len(ld_expressed)})", "green"),
        (ld_gdna, f"gDNA seed (n={len(ld_gdna)})", "red"),
        (ld_ambig, f"Ambiguous (n={len(ld_ambig)})", "gray"),
    ]:
        if len(vals) > 0:
            ax.hist(vals, bins=30, alpha=0.5, color=color, label=label, density=True, edgecolor="white")
    ax.set_xlabel("log(n/L + ε)")
    ax.set_ylabel("Density")
    ax.set_title("Seed Partition by Log Density")
    ax.legend(fontsize=6)

    # 2) Sense fraction colored by seed label
    ax = fig.add_subplot(gs[0, 1])
    for mask, label, color in [
        (expressed_seed, "Expressed", "green"),
        (gdna_seed, "gDNA seed", "red"),
        (ambiguous, "Ambiguous", "gray"),
    ]:
        sf_vals = sense_frac[mask]
        sf_vals = sf_vals[np.isfinite(sf_vals)]
        if len(sf_vals) > 0:
            ax.hist(sf_vals, bins=30, alpha=0.5, color=color, label=label, density=True, edgecolor="white")
    ax.set_xlabel("Sense Fraction")
    ax.set_ylabel("Density")
    ax.set_title("Seed Partition by Sense Fraction")
    ax.legend(fontsize=6)

    # 3) Log density vs sense fraction scatter colored by seed
    ax = fig.add_subplot(gs[0, 2])
    for mask, label, color, marker in [
        (expressed_seed, "Expressed", "green", "o"),
        (gdna_seed, "gDNA seed", "red", "s"),
        (ambiguous, "Ambiguous", "gray", "."),
    ]:
        sf_m = sense_frac[mask]
        ld_m = log_d[mask]
        valid = np.isfinite(sf_m)
        if valid.any():
            ax.scatter(sf_m[valid], ld_m[valid], c=color, s=8, alpha=0.4,
                      label=label, marker=marker)
    ax.set_xlabel("Sense Fraction")
    ax.set_ylabel("log(n/L + ε)")
    ax.set_title("Seed Partition (2D)")
    ax.legend(fontsize=6, markerscale=2)

    # 4) Initialization summary text
    ax = fig.add_subplot(gs[1, 0])
    ax.axis("off")
    text_lines = [
        f"n_expressed_seed: {init_diag.get('n_expressed_seed', 'N/A')}",
        f"n_gdna_seed: {init_diag.get('n_gdna_seed', 'N/A')}",
        f"pi_init: {init_diag.get('pi_init', 'N/A'):.4f}" if isinstance(init_diag.get('pi_init'), float) else f"pi_init: {init_diag.get('pi_init', 'N/A')}",
        f"pristine_sample: {init_diag.get('pristine_sample', 'N/A')}",
        f"expressed_fallback: {init_diag.get('expressed_fallback', False)}",
        f"n_eligible: {int(eligible.sum())}",
        f"n_total_regions: {len(eligible)}",
    ]
    ax.text(0.05, 0.95, "\n".join(text_lines), transform=ax.transAxes,
            verticalalignment="top", fontfamily="monospace", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
    ax.set_title("Initialization Summary")

    # 5) Initial density histograms (from first M-step)
    if result.iteration_history and len(result.iteration_history) > 0:
        h = result.iteration_history[0]
        ax = fig.add_subplot(gs[1, 1])
        _hist_bar(ax, *h["gdna_density_hist"], "red", "gDNA", alpha=0.5)
        _hist_bar(ax, *h["rna_density_hist"], "green", "RNA", alpha=0.5)
        ax.set_xlabel("log(n/L + ε)")
        ax.set_ylabel("Probability")
        ax.set_title("Density Histograms (iter 1)")
        ax.legend(fontsize=7)

    # 6) Initial strand histograms (from first M-step)
    if result.iteration_history and len(result.iteration_history) > 0:
        h = result.iteration_history[0]
        ax = fig.add_subplot(gs[1, 2])
        _hist_bar(ax, *h["gdna_strand_hist"], "red", "gDNA", alpha=0.5)
        _hist_bar(ax, *h["rna_strand_hist"], "green", "RNA", alpha=0.5)
        ax.set_xlabel("Sense Fraction")
        ax.set_ylabel("Probability")
        ax.set_title("Strand Histograms (iter 1)")
        ax.legend(fontsize=7)

    fig.suptitle("Initialization Quality", fontsize=12, fontweight="bold")
    fig.tight_layout()
    paths.append(_save_fig(fig, plot_dir, "05_initialization"))

    return paths


# ===================================================================
# Section 4: EM Convergence
# ===================================================================


def plot_em_convergence(
    result: CalibrationResult,
    eligible: np.ndarray,
    plot_dir: Path,
) -> list[str]:
    """Plot EM iteration traces and histogram evolution."""
    paths = []
    history = result.iteration_history
    if not history or len(history) == 0:
        return paths

    iters = list(range(1, len(history) + 1))

    # --- Convergence traces ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(iters, [h["pi"] for h in history], "o-", color="steelblue", markersize=3)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("π (mixing proportion)")
    ax.set_title("π Convergence")
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(iters, [h["delta_pi"] for h in history], "o-", color="crimson", markersize=3)
    ax.set_yscale("log")
    ax.axhline(result._convergence_tol if hasattr(result, "_convergence_tol") else 1e-4,
               color="gray", ls="--", lw=1, label="tolerance")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("|Δπ|")
    ax.set_title("Convergence Rate")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(iters, [h["lambda_G"] for h in history], "o-", color="red", markersize=3, label="λ_G")
    ax.plot(iters, [h["lambda_E"] for h in history], "o-", color="green", markersize=3, label="λ_E")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Density (frags/bp)")
    ax.set_title("Component Densities")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(iters, [h["mean_gamma"] for h in history], "o-", color="darkorange", markersize=3)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Mean γ (eligible)")
    ax.set_title("Mean Posterior")
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"EM Convergence ({result.n_iterations} iters, "
        f"{'converged' if result.converged else 'NOT converged'})",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout()
    paths.append(_save_fig(fig, plot_dir, "06_convergence_traces"))

    # --- Histogram evolution (first, middle, last iteration) ---
    if len(history) >= 2:
        pick = [0]
        if len(history) > 2:
            pick.append(len(history) // 2)
        pick.append(len(history) - 1)

        fig, axes = plt.subplots(len(pick), 2, figsize=(12, 4 * len(pick)))
        if len(pick) == 1:
            axes = axes[np.newaxis, :]

        for row, idx in enumerate(pick):
            h = history[idx]
            it = idx + 1

            # Density histograms
            ax = axes[row, 0]
            _hist_bar(ax, *h["gdna_density_hist"], "red", "gDNA", alpha=0.5)
            _hist_bar(ax, *h["rna_density_hist"], "green", "RNA", alpha=0.5)
            ax.set_xlabel("log(n/L + ε)")
            ax.set_ylabel("Probability")
            ax.set_title(f"Density Histograms (iter {it})")
            ax.legend(fontsize=7)

            # Strand histograms
            ax = axes[row, 1]
            _hist_bar(ax, *h["gdna_strand_hist"], "red", "gDNA", alpha=0.5)
            _hist_bar(ax, *h["rna_strand_hist"], "green", "RNA", alpha=0.5)
            ax.set_xlabel("Sense Fraction")
            ax.set_ylabel("Probability")
            ax.set_title(f"Strand Histograms (iter {it})")
            ax.legend(fontsize=7)

        fig.suptitle("Histogram Evolution Across EM", fontsize=12, fontweight="bold")
        fig.tight_layout()
        paths.append(_save_fig(fig, plot_dir, "07_histogram_evolution"))

    # --- Posterior distribution evolution ---
    if len(history) >= 2:
        pick = [0]
        if len(history) > 2:
            pick.append(len(history) // 2)
        pick.append(len(history) - 1)

        fig, axes = plt.subplots(1, len(pick), figsize=(5 * len(pick), 4))
        if len(pick) == 1:
            axes = [axes]
        for i, idx in enumerate(pick):
            h = history[idx]
            gamma = h["gamma"][eligible]
            axes[i].hist(gamma, bins=50, color="steelblue", alpha=0.8, edgecolor="white")
            axes[i].set_xlabel("γ (posterior)")
            axes[i].set_ylabel("Count")
            axes[i].set_title(f"Iteration {idx + 1} (π={h['pi']:.3f})")
        fig.suptitle("Posterior Distribution Evolution", fontsize=12, fontweight="bold")
        fig.tight_layout()
        paths.append(_save_fig(fig, plot_dir, "08_posterior_evolution"))

    return paths


# ===================================================================
# Section 5: Final Results
# ===================================================================


def plot_final_results(
    result: CalibrationResult,
    log_d: np.ndarray,
    sense_frac: np.ndarray,
    eligible: np.ndarray,
    stats: dict[str, np.ndarray],
    plot_dir: Path,
) -> list[str]:
    """Plot final EM results."""
    paths = []
    gamma = result.region_posteriors

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig)

    # 1) Final posterior distribution
    ax = fig.add_subplot(gs[0, 0])
    g_elig = gamma[eligible]
    ax.hist(g_elig, bins=50, color="steelblue", alpha=0.8, edgecolor="white")
    ax.axvline(0.5, color="red", ls="--", lw=1)
    n_gdna = (g_elig > 0.5).sum()
    n_rna = (g_elig <= 0.5).sum()
    ax.set_xlabel("γ (posterior P(not expressed))")
    ax.set_ylabel("Count")
    ax.set_title(f"Final Posteriors (gDNA={n_gdna}, RNA={n_rna})")

    # 2) Log density colored by posterior
    ax = fig.add_subplot(gs[0, 1])
    ld_elig = log_d[eligible]
    sc = ax.scatter(ld_elig, g_elig, c=g_elig, cmap="RdYlGn_r", s=6, alpha=0.6,
                    vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="γ")
    ax.set_xlabel("log(n/L + ε)")
    ax.set_ylabel("γ")
    ax.set_title("Posterior vs Log Density")

    # 3) Sense fraction colored by posterior
    ax = fig.add_subplot(gs[0, 2])
    sf_elig = sense_frac[eligible]
    valid = np.isfinite(sf_elig)
    if valid.any():
        sc = ax.scatter(sf_elig[valid], g_elig[valid], c=g_elig[valid],
                        cmap="RdYlGn_r", s=6, alpha=0.6, vmin=0, vmax=1)
        plt.colorbar(sc, ax=ax, label="γ")
    ax.set_xlabel("Sense Fraction")
    ax.set_ylabel("γ")
    ax.set_title("Posterior vs Sense Fraction")

    # 4) Density distributions for classified regions
    ax = fig.add_subplot(gs[1, 0])
    is_gdna = eligible & (gamma > 0.5)
    is_rna = eligible & (gamma <= 0.5) & (stats["n_spliced"] == 0)
    is_spliced = eligible & (stats["n_spliced"] > 0)
    for mask, label, color in [
        (is_gdna, f"gDNA (γ>0.5, n={is_gdna.sum()})", "red"),
        (is_rna, f"RNA unspliced (γ≤0.5, n={is_rna.sum()})", "blue"),
        (is_spliced, f"RNA spliced (n={is_spliced.sum()})", "green"),
    ]:
        vals = log_d[mask]
        if len(vals) > 0:
            ax.hist(vals, bins=30, alpha=0.4, color=color, label=label, density=True, edgecolor="white")
    ax.set_xlabel("log(n/L + ε)")
    ax.set_ylabel("Density")
    ax.set_title("Density by Classification")
    ax.legend(fontsize=6)

    # 5) Sense fraction distributions for classified regions
    ax = fig.add_subplot(gs[1, 1])
    for mask, label, color in [
        (is_gdna, "gDNA (γ>0.5)", "red"),
        (is_rna, "RNA unspliced (γ≤0.5)", "blue"),
        (is_spliced, "RNA spliced", "green"),
    ]:
        vals = sense_frac[mask]
        vals = vals[np.isfinite(vals)]
        if len(vals) > 0:
            ax.hist(vals, bins=30, alpha=0.4, color=color, label=label, density=True, edgecolor="white")
    ax.set_xlabel("Sense Fraction")
    ax.set_ylabel("Density")
    ax.set_title("Sense Fraction by Classification")
    ax.legend(fontsize=6)

    # 6) Per-chromosome gDNA density
    ax = fig.add_subplot(gs[1, 2])
    per_ref = result.gdna_density_per_ref
    if per_ref:
        refs = sorted(per_ref.keys())
        densities = [per_ref[r] for r in refs]
        ax.barh(refs, densities, color="coral", alpha=0.8)
        ax.axvline(result.gdna_density_global, color="black", ls="--", lw=1,
                   label=f"global={result.gdna_density_global:.4f}")
        ax.set_xlabel("gDNA density (frags/bp)")
        ax.set_title("Per-Chromosome gDNA Density")
        ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "No per-ref density data", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title("Per-Chromosome gDNA Density")

    fig.suptitle("Final Calibration Results", fontsize=12, fontweight="bold")
    fig.tight_layout()
    paths.append(_save_fig(fig, plot_dir, "09_final_results"))

    # --- Final density and strand histograms ---
    if result.iteration_history and len(result.iteration_history) > 0:
        h = result.iteration_history[-1]
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

        ax = axes[0]
        _hist_bar(ax, *h["gdna_density_hist"], "red", "gDNA", alpha=0.5)
        _hist_bar(ax, *h["rna_density_hist"], "green", "RNA", alpha=0.5)
        ax.set_xlabel("log(n/L + ε)")
        ax.set_ylabel("Probability")
        ax.set_title("Final Density Histograms")
        ax.legend()

        ax = axes[1]
        _hist_bar(ax, *h["gdna_strand_hist"], "red", "gDNA", alpha=0.5)
        _hist_bar(ax, *h["rna_strand_hist"], "green", "RNA", alpha=0.5)
        ax.set_xlabel("Sense Fraction")
        ax.set_ylabel("Probability")
        ax.set_title("Final Strand Histograms")
        ax.legend()

        fig.suptitle("Final Component Histograms", fontsize=12, fontweight="bold")
        fig.tight_layout()
        paths.append(_save_fig(fig, plot_dir, "10_final_histograms"))

    return paths


# ===================================================================
# Section 6: Ground Truth Comparison
# ===================================================================


def plot_ground_truth(
    result: CalibrationResult,
    truth_meta: dict,
    eligible: np.ndarray,
    log_d: np.ndarray,
    sense_frac: np.ndarray,
    plot_dir: Path,
) -> list[str]:
    """Plot ground truth comparison plots (simulation only)."""
    paths = []
    gamma = result.region_posteriors

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # 1) True gDNA fraction vs estimated π
    ax = axes[0]
    true_frac = truth_meta["gdna_fraction"]
    est_pi = result.mixing_proportion
    ax.bar(["True gDNA\nfraction", "Estimated π"], [true_frac, est_pi],
           color=["steelblue", "coral"], alpha=0.8, edgecolor="black")
    ax.set_ylabel("Fraction")
    ax.set_title(f"gDNA Fraction: True={true_frac:.3f}, Est={est_pi:.3f}")
    ax.set_ylim(0, max(true_frac, est_pi) * 1.3 + 0.01)
    for i, v in enumerate([true_frac, est_pi]):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=9)

    # 2) Posterior distribution colored by truth (spliced = definitely RNA)
    ax = axes[1]
    has_spliced = truth_meta["has_spliced"]
    no_spliced = ~has_spliced & eligible
    g_rna = gamma[has_spliced & eligible]
    g_unknown = gamma[no_spliced]
    if len(g_rna) > 0:
        ax.hist(g_rna, bins=30, alpha=0.5, color="green",
                label=f"Spliced (RNA, n={len(g_rna)})", density=True, edgecolor="white")
    if len(g_unknown) > 0:
        ax.hist(g_unknown, bins=30, alpha=0.5, color="gray",
                label=f"Unspliced-only (n={len(g_unknown)})", density=True, edgecolor="white")
    ax.set_xlabel("γ (posterior)")
    ax.set_ylabel("Density")
    ax.set_title("Posteriors by Truth Category")
    ax.legend(fontsize=7)

    # 3) Summary statistics
    ax = axes[2]
    ax.axis("off")
    stats_text = [
        f"True gDNA frags:  {truth_meta['n_gdna_frags']:,}",
        f"True RNA frags:   {truth_meta['n_rna_frags']:,}",
        f"True gDNA frac:   {truth_meta['gdna_fraction']:.4f}",
        "",
        f"Estimated π:      {result.mixing_proportion:.4f}",
        f"Estimated λ_G:    {result.gdna_density_global:.4e}",
        f"Estimated λ_E:    {result.expressed_density:.4e}",
        f"κ_sym:            {result.kappa_sym:.1f}",
        "",
        f"Converged:        {result.converged}",
        f"Iterations:       {result.n_iterations}",
        f"n_spliced_regions:{has_spliced.sum()}",
        f"n_eligible:       {eligible.sum()}",
    ]
    ax.text(0.05, 0.95, "\n".join(stats_text), transform=ax.transAxes,
            verticalalignment="top", fontfamily="monospace", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
    ax.set_title("Ground Truth Summary")

    fig.suptitle("Ground Truth Comparison", fontsize=12, fontweight="bold")
    fig.tight_layout()
    paths.append(_save_fig(fig, plot_dir, "11_ground_truth"))

    return paths


# ===================================================================
# Markdown Report Generator
# ===================================================================


def generate_report(
    result: CalibrationResult,
    bam_path: str,
    ss: float,
    plot_paths: dict[str, list[str]],
    truth_meta: dict | None = None,
    run_params: dict | None = None,
) -> str:
    """Generate the full markdown report."""
    stats = result.region_stats
    eligible = result.eligible
    n_total_regions = len(result.region_posteriors)
    n_eligible = int(eligible.sum()) if eligible is not None else n_total_regions
    gamma = result.region_posteriors

    lines = []
    lines.append("# Calibration Diagnostics Report\n")
    lines.append(f"**BAM**: `{bam_path}`\n")
    lines.append(f"**Strand Specificity**: {ss:.4f}\n")
    lines.append(f"**Date**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")
    lines.append("")

    # --- Run Parameters ---
    if run_params:
        lines.append("## Run Parameters\n")
        lines.append("| Parameter | Value |")
        lines.append("|---|---|")
        for k, v in sorted(run_params.items()):
            lines.append(f"| {k} | `{v}` |")
        lines.append("")

    # --- Summary Table ---
    lines.append("## Summary\n")
    lines.append("| Parameter | Value |")
    lines.append("|---|---|")
    lines.append(f"| Total regions | {n_total_regions:,} |")
    lines.append(f"| Eligible regions | {n_eligible:,} |")
    lines.append(f"| Spliced regions | {int((stats['n_spliced'] > 0).sum()):,} |" if stats else "")
    lines.append(f"| Mixing proportion (π) | {result.mixing_proportion:.4f} |")
    lines.append(f"| gDNA density (λ_G) | {result.gdna_density_global:.4e} |")
    lines.append(f"| Expressed density (λ_E) | {result.expressed_density:.4e} |")
    lines.append(f"| κ_sym | {result.kappa_sym:.1f} |")
    lines.append(f"| Converged | {result.converged} |")
    lines.append(f"| Iterations | {result.n_iterations} |")

    if eligible is not None:
        g_elig = gamma[eligible]
        lines.append(f"| Classified gDNA (γ > 0.5) | {int((g_elig > 0.5).sum()):,} |")
        lines.append(f"| Classified RNA (γ ≤ 0.5) | {int((g_elig <= 0.5).sum()):,} |")
        lines.append(f"| Mean γ (eligible) | {float(g_elig.mean()):.4f} |")
    lines.append("")

    if truth_meta:
        lines.append("### Ground Truth\n")
        lines.append("| Parameter | Value |")
        lines.append("|---|---|")
        lines.append(f"| True gDNA fragments | {truth_meta['n_gdna_frags']:,} |")
        lines.append(f"| True RNA fragments | {truth_meta['n_rna_frags']:,} |")
        lines.append(f"| True gDNA fraction | {truth_meta['gdna_fraction']:.4f} |")
        lines.append(f"| Estimated π | {result.mixing_proportion:.4f} |")
        pi_err = abs(result.mixing_proportion - truth_meta['gdna_fraction'])
        lines.append(f"| π absolute error | {pi_err:.4f} |")
        lines.append("")

    # --- Section 1: Input Data ---
    lines.append("## 1. Input Data Summary\n")
    for p in plot_paths.get("input", []):
        lines.append(f"![Input Summary]({p})\n")

    # --- Section 2: Empirical Histograms ---
    lines.append("## 2. Empirical Distributions (Pre-EM)\n")
    for p in plot_paths.get("empirical", []):
        lines.append(f"![Empirical]({p})\n")

    # --- Section 3: Initialization ---
    lines.append("## 3. Initialization Quality\n")
    if result.init_diagnostics:
        init = result.init_diagnostics
        lines.append("| Parameter | Value |")
        lines.append("|---|---|")
        lines.append(f"| Expressed seed regions | {init.get('n_expressed_seed', 'N/A')} |")
        lines.append(f"| gDNA seed regions | {init.get('n_gdna_seed', 'N/A')} |")
        pi_init = init.get('pi_init')
        lines.append(f"| Initial π | {pi_init:.4f} |" if isinstance(pi_init, float) else f"| Initial π | {pi_init} |")
        lines.append(f"| Pristine sample | {init.get('pristine_sample', 'N/A')} |")
        lines.append(f"| Expressed fallback | {init.get('expressed_fallback', False)} |")
        lines.append("")
    for p in plot_paths.get("init", []):
        lines.append(f"![Initialization]({p})\n")

    # --- Section 4: EM Convergence ---
    lines.append("## 4. EM Convergence\n")
    if result.iteration_history:
        lines.append("### Iteration Table\n")
        lines.append("| Iter | π | λ_G | λ_E | Mean γ | |Δπ| |")
        lines.append("|---:|---:|---:|---:|---:|---:|")
        for i, h in enumerate(result.iteration_history):
            lines.append(
                f"| {i + 1} | {h['pi']:.4f} | {h['lambda_G']:.4e} | "
                f"{h['lambda_E']:.4e} | {h['mean_gamma']:.4f} | "
                f"{h['delta_pi']:.2e} |"
            )
        lines.append("")
    for p in plot_paths.get("convergence", []):
        lines.append(f"![Convergence]({p})\n")

    # --- Section 5: Final Results ---
    lines.append("## 5. Final Results\n")

    # Per-chromosome density table
    if result.gdna_density_per_ref:
        lines.append("### Per-Chromosome gDNA Density\n")
        lines.append("| Chromosome | Density |")
        lines.append("|---|---:|")
        for ref, d in sorted(result.gdna_density_per_ref.items()):
            lines.append(f"| {ref} | {d:.4e} |")
        lines.append("")

    for p in plot_paths.get("final", []):
        lines.append(f"![Final Results]({p})\n")

    # --- Section 6: Ground Truth (optional) ---
    if truth_meta:
        lines.append("## 6. Ground Truth Comparison\n")
        for p in plot_paths.get("truth", []):
            lines.append(f"![Ground Truth]({p})\n")

    return "\n".join(lines)


# ===================================================================
# Self-contained HTML generation
# ===================================================================

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Calibration Diagnostics Report</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI",
         Helvetica, Arial, sans-serif; line-height: 1.6;
         max-width: 1100px; margin: 0 auto; padding: 20px 40px;
         color: #24292f; background: #fff; }}
  h1 {{ border-bottom: 2px solid #d0d7de; padding-bottom: 8px; }}
  h2 {{ border-bottom: 1px solid #d0d7de; padding-bottom: 4px;
       margin-top: 32px; }}
  h3 {{ margin-top: 24px; }}
  table {{ border-collapse: collapse; margin: 12px 0; }}
  th, td {{ border: 1px solid #d0d7de; padding: 6px 12px;
           text-align: left; }}
  th {{ background: #f6f8fa; }}
  tr:nth-child(even) {{ background: #f6f8fa; }}
  img {{ max-width: 100%; height: auto; margin: 8px 0;
        border: 1px solid #d0d7de; border-radius: 4px; }}
  code {{ background: #f6f8fa; padding: 2px 6px; border-radius: 3px;
         font-size: 0.9em; }}
  strong {{ font-weight: 600; }}
  p {{ margin: 8px 0; }}
</style>
</head>
<body>
{body}
</body>
</html>
"""


def _md_to_html_simple(md: str, output_dir: Path) -> str:
    """Convert markdown report to HTML with base64-embedded images.

    Handles the subset of markdown used by generate_report():
    headings, tables, bold, code, image links, and paragraphs.
    """
    import html as html_mod

    out_lines: list[str] = []
    in_table = False
    is_header_row = False

    for line in md.split("\n"):
        stripped = line.strip()

        # Close table if we leave table context
        if in_table and not stripped.startswith("|"):
            out_lines.append("</table>")
            in_table = False

        # Headings
        if stripped.startswith("#"):
            m = re.match(r"^(#{1,6})\s+(.*)", stripped)
            if m:
                level = len(m.group(1))
                text = m.group(2)
                text = _inline_fmt(text)
                out_lines.append(f"<h{level}>{text}</h{level}>")
                continue

        # Image: ![alt](path)
        img_m = re.match(r"^!\[([^\]]*)\]\(([^)]+)\)\s*$", stripped)
        if img_m:
            alt = html_mod.escape(img_m.group(1))
            img_rel = img_m.group(2)
            img_path = output_dir / img_rel
            if img_path.exists():
                b64 = base64.b64encode(img_path.read_bytes()).decode("ascii")
                out_lines.append(
                    f'<img src="data:image/png;base64,{b64}" alt="{alt}">'
                )
            else:
                out_lines.append(f'<p><em>[Missing image: {alt}]</em></p>')
            continue

        # Table rows
        if stripped.startswith("|"):
            cells = [c.strip() for c in stripped.split("|")[1:-1]]
            # Skip separator rows (|---|---|)
            if all(re.match(r"^-+:?$|^:?-+:?$", c) for c in cells):
                is_header_row = False
                continue
            if not in_table:
                out_lines.append("<table>")
                in_table = True
                is_header_row = True
            tag = "th" if is_header_row else "td"
            row = "".join(f"<{tag}>{_inline_fmt(c)}</{tag}>" for c in cells)
            out_lines.append(f"<tr>{row}</tr>")
            is_header_row = False
            continue

        # Empty line
        if not stripped:
            continue

        # Paragraph
        out_lines.append(f"<p>{_inline_fmt(stripped)}</p>")

    if in_table:
        out_lines.append("</table>")

    return "\n".join(out_lines)


def _inline_fmt(text: str) -> str:
    """Apply inline markdown formatting: bold, code, italic."""
    import html as html_mod
    # Code spans first (so bold inside code isn't processed)
    text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
    # Bold
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    # Italic
    text = re.sub(r"\*([^*]+)\*", r"<em>\1</em>", text)
    return text


def generate_html(md_report: str, output_dir: Path) -> str:
    """Convert a markdown report to a self-contained HTML file.

    All ``![alt](plots/foo.png)`` references are replaced with
    base64-encoded data URIs so the HTML file is fully portable.
    """
    body = _md_to_html_simple(md_report, output_dir)
    return _HTML_TEMPLATE.format(body=body)


# ===================================================================
# Main entry point
# ===================================================================


def run_diagnostics(
    bam_path: str,
    index_dir: str,
    output_dir: str,
    *,
    strand_specificity: float | None = None,
    max_iterations: int = 50,
    convergence_tol: float = 1e-4,
    parse_truth: bool = False,
    emit_html: bool = True,
    run_params: dict | None = None,
) -> Path:
    """Run calibration diagnostics and write report.

    Parameters
    ----------
    bam_path : str
        Path to name-sorted BAM file.
    index_dir : str
        Path to rigel index directory.
    output_dir : str
        Output directory for report and plots.
    strand_specificity : float or None
        Override SS. If None, estimated from spliced reads.
    max_iterations : int
        Maximum EM iterations.
    convergence_tol : float
        |Δπ| convergence tolerance.
    parse_truth : bool
        Parse ground truth labels from simulated read names.
    emit_html : bool
        Also write a self-contained HTML report (default True).
    run_params : dict or None
        Optional dict of run parameters to include in the report.

    Returns
    -------
    Path to the generated markdown report.
    """
    out = Path(output_dir)
    plot_dir = out / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # --- Load index ---
    logger.info(f"Loading index from {index_dir}")
    index = TranscriptIndex.load(index_dir)
    if index.region_df is None or index.region_cr is None:
        raise RuntimeError(
            "Index does not contain calibration regions. "
            "Rebuild with a version that generates region_df."
        )

    # --- Count region evidence ---
    logger.info(f"Counting region evidence from {bam_path}")
    region_counts, fl_table = count_region_evidence(bam_path, index)
    logger.info(
        f"Region evidence: {len(region_counts)} regions, "
        f"{len(fl_table)} FL observations"
    )

    # --- Compute stats and estimate SS ---
    stats = compute_region_stats(region_counts, index.region_df)
    eligible = (stats["n_total"] > 0) & (stats["region_length"] > 0)

    if strand_specificity is None:
        strand_specificity = estimate_strand_specificity(stats)
    ss = strand_specificity

    # --- Run calibration with full diagnostics ---
    logger.info(f"Running calibration (SS={ss:.4f}, max_iter={max_iterations})")
    result = calibrate_gdna(
        region_counts, fl_table, index.region_df,
        strand_specificity=ss,
        max_iterations=max_iterations,
        convergence_tol=convergence_tol,
        diagnostics=True,
    )
    logger.info(
        f"Calibration complete: π={result.mixing_proportion:.4f}, "
        f"converged={result.converged}, iters={result.n_iterations}"
    )

    # Use diagnostic fields
    log_d = result.log_density
    sense_frac = result.sense_frac
    eligible = result.eligible
    assert log_d is not None and sense_frac is not None and eligible is not None

    # --- Ground truth ---
    truth_meta = None
    if parse_truth:
        logger.info("Parsing ground truth from BAM read names")
        truth_meta = parse_ground_truth_labels(bam_path, region_counts)

    # --- Generate plots ---
    plot_paths: dict[str, list[str]] = {}
    logger.info("Generating plots...")

    plot_paths["input"] = plot_input_summary(stats, eligible, plot_dir)
    plot_paths["empirical"] = plot_empirical_histograms(
        sense_frac, log_d, eligible, stats, plot_dir,
    )
    plot_paths["init"] = plot_initialization(
        result, log_d, sense_frac, eligible, plot_dir,
    )
    plot_paths["convergence"] = plot_em_convergence(result, eligible, plot_dir)
    plot_paths["final"] = plot_final_results(
        result, log_d, sense_frac, eligible, stats, plot_dir,
    )

    if truth_meta:
        plot_paths["truth"] = plot_ground_truth(
            result, truth_meta, eligible, log_d, sense_frac, plot_dir,
        )

    # --- Build run_params if not provided ---
    if run_params is None:
        run_params = {}
    run_params.setdefault("bam_path", str(bam_path))
    run_params.setdefault("index_dir", str(index_dir))
    run_params.setdefault("strand_specificity", f"{ss:.4f}")
    run_params.setdefault("max_iterations", str(max_iterations))
    run_params.setdefault("convergence_tol", str(convergence_tol))
    run_params.setdefault("parse_truth", str(parse_truth))

    # --- Generate markdown report ---
    logger.info("Generating markdown report...")
    report = generate_report(
        result, bam_path, ss, plot_paths, truth_meta, run_params,
    )
    report_path = out / "calibration_report.md"
    report_path.write_text(report)
    logger.info(f"Report written to {report_path}")

    # --- Generate self-contained HTML ---
    if emit_html:
        logger.info("Generating self-contained HTML report...")
        html_content = generate_html(report, out)
        html_path = out / "calibration_report.html"
        html_path.write_text(html_content)
        logger.info(f"HTML report written to {html_path}")

    return report_path


# ===================================================================
# CLI
# ===================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Calibration diagnostics for gDNA deconvolution",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              python scripts/calibration_diagnostics.py --bam sample.bam --index index/ -o diag/
              python scripts/calibration_diagnostics.py --bam sim.bam --index index/ -o diag/ --truth
              python scripts/calibration_diagnostics.py --bam sample.bam --index index/ -o diag/ --ss 0.95
        """),
    )
    parser.add_argument("--bam", required=True, help="Path to name-sorted BAM file")
    parser.add_argument("--index", required=True, help="Path to rigel index directory")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--ss", type=float, default=None,
                        help="Override strand specificity (auto-estimated if omitted)")
    parser.add_argument("--truth", action="store_true",
                        help="Parse ground truth from simulated read names")
    parser.add_argument("--max-iter", type=int, default=50,
                        help="Maximum EM iterations (default: 50)")
    parser.add_argument("--convergence-tol", type=float, default=1e-4,
                        help="Convergence tolerance for |Δπ| (default: 1e-4)")
    parser.add_argument("--no-html", action="store_true",
                        help="Skip HTML report generation")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose logging")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )

    run_diagnostics(
        bam_path=args.bam,
        index_dir=args.index,
        output_dir=args.output_dir,
        strand_specificity=args.ss,
        max_iterations=args.max_iter,
        convergence_tol=args.convergence_tol,
        parse_truth=args.truth,
        emit_html=not args.no_html,
    )


if __name__ == "__main__":
    main()
