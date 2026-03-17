#!/usr/bin/env python
"""Calibration stress-test framework.

Generates large synthetic region ensembles with known ground truth and
runs ``calibrate_gdna`` to characterize classification accuracy across
a grid of conditions:

    - strand specificity (0.5 → 1.0)
    - gDNA contamination level (low → dominant)
    - fragment-length overlap (identical → well-separated)
    - nascent-RNA fraction (0 → high)

Outputs a summary table and optional per-scenario diagnostic plots.

Usage
-----
    conda activate rigel
    python scripts/debug/calibration_stress_test.py          # default grid
    python scripts/debug/calibration_stress_test.py --quick   # 3×3 quick scan
    python scripts/debug/calibration_stress_test.py --outdir results/stress
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

# ── Rigel imports ──────────────────────────────────────────────────
from rigel.calibration import (
    GDNACalibration,
    calibrate_gdna,
    compute_region_stats,
    compute_sense_fraction,
)

# ===================================================================
# Scenario configuration
# ===================================================================


@dataclass
class ScenarioConfig:
    """Configuration for a single synthetic calibration scenario."""

    # -- region counts --
    n_rna_regions: int = 500
    n_gdna_regions: int = 500
    n_nrna_regions: int = 0

    # -- strand --
    strand_specificity: float = 0.95
    # Fraction of genes on + strand (rest on -)
    plus_strand_frac: float = 0.5

    # -- fragment-length --
    rna_fl_mean: float = 250.0
    rna_fl_std: float = 50.0
    gdna_fl_mean: float = 250.0
    gdna_fl_std: float = 80.0

    # -- read counts --
    rna_count_mean: float = 200.0  # mean total fragments per RNA region
    rna_count_std: float = 60.0
    gdna_count_mean: float = 50.0  # mean total fragments per gDNA region
    gdna_count_std: float = 20.0

    # -- nRNA --
    nrna_splice_rate: float = 0.0  # fraction of nRNA fragments that are spliced
    nrna_count_mean: float = 80.0
    nrna_count_std: float = 30.0
    # nRNA strand: same as RNA (reads on antisense strand)

    # -- region geometry --
    region_length_mean: float = 2000.0
    region_length_std: float = 500.0

    # -- splice signal --
    rna_splice_rate: float = 0.3  # fraction of RNA frags that are spliced
    gdna_splice_rate: float = 0.0  # gDNA should never be spliced

    # -- misc --
    seed: int = 42
    label: str = ""

    @property
    def true_pi(self) -> float:
        """True fraction of not-expressed (gDNA-only) regions."""
        total = self.n_rna_regions + self.n_gdna_regions + self.n_nrna_regions
        return self.n_gdna_regions / total if total > 0 else 0.0


# ===================================================================
# Synthetic data generator
# ===================================================================


def generate_scenario(cfg: ScenarioConfig) -> dict:
    """Generate synthetic region_counts, region_df, fl_table, and truth labels.

    Returns dict with keys:
        region_counts, region_df, fl_table, truth_labels, config
    where truth_labels is a 1D array: 0=RNA, 1=gDNA, 2=nRNA
    """
    rng = np.random.default_rng(cfg.seed)
    n_total = cfg.n_rna_regions + cfg.n_gdna_regions + cfg.n_nrna_regions

    # -- truth labels --
    truth = np.concatenate([
        np.zeros(cfg.n_rna_regions, dtype=np.int8),       # 0 = RNA
        np.ones(cfg.n_gdna_regions, dtype=np.int8),        # 1 = gDNA
        np.full(cfg.n_nrna_regions, 2, dtype=np.int8),    # 2 = nRNA
    ])

    # -- gene strand --
    is_plus = rng.random(n_total) < cfg.plus_strand_frac
    tx_pos = is_plus.copy()
    tx_neg = ~is_plus

    # -- region lengths --
    lengths = rng.normal(cfg.region_length_mean, cfg.region_length_std, n_total)
    lengths = np.clip(lengths, 200, 20000).astype(np.int32)

    # -- generate counts per component --
    n_unspliced_pos = np.zeros(n_total, dtype=np.float32)
    n_unspliced_neg = np.zeros(n_total, dtype=np.float32)
    n_spliced_pos = np.zeros(n_total, dtype=np.float32)
    n_spliced_neg = np.zeros(n_total, dtype=np.float32)
    fl_region_ids = []
    fl_frag_lens = []

    for i in range(n_total):
        component = truth[i]
        is_plus_gene = is_plus[i]

        if component == 0:
            # RNA region
            n_frags = max(1, int(rng.normal(cfg.rna_count_mean, cfg.rna_count_std)))
            n_spliced = int(rng.binomial(n_frags, cfg.rna_splice_rate))
            n_unspliced = n_frags - n_spliced
            ss = cfg.strand_specificity
        elif component == 1:
            # gDNA region
            n_frags = max(1, int(rng.normal(cfg.gdna_count_mean, cfg.gdna_count_std)))
            n_spliced = int(rng.binomial(n_frags, cfg.gdna_splice_rate))
            n_unspliced = n_frags - n_spliced
            ss = 0.5  # gDNA is unstranded
        else:
            # nRNA region
            n_frags = max(1, int(rng.normal(cfg.nrna_count_mean, cfg.nrna_count_std)))
            n_spliced = int(rng.binomial(n_frags, cfg.nrna_splice_rate))
            n_unspliced = n_frags - n_spliced
            ss = cfg.strand_specificity  # nRNA is stranded like RNA

        # Strand assignment for unspliced reads
        # R1-antisense convention: RNA on + gene → reads on - strand
        # So sense direction reads = n_unspliced * (1 - ss)
        # For + gene: n_pos = antisense portion (1 - ss), n_neg = sense portion (ss)
        # Wait, let's be precise:
        #   + gene: sense strand is -, so sense reads are neg
        #   + gene with SS: fraction (ss) of reads on antisense = neg strand
        #   Actually for R1-antisense: RNA from + gene → R1 maps to - strand
        #   So for + gene RNA: n_neg = ss * n_unspliced, n_pos = (1-ss) * n_unspliced
        #   For - gene RNA: n_pos = ss * n_unspliced, n_neg = (1-ss) * n_unspliced
        #   For gDNA (ss=0.5): symmetric
        if is_plus_gene:
            # + gene: RNA reads land on neg strand (R1-antisense)
            p_pos = 1.0 - ss
        else:
            # - gene: RNA reads land on pos strand
            p_pos = ss

        n_pos_us = int(rng.binomial(n_unspliced, p_pos))
        n_neg_us = n_unspliced - n_pos_us

        # Spliced reads have same strand bias
        n_pos_sp = int(rng.binomial(n_spliced, p_pos))
        n_neg_sp = n_spliced - n_pos_sp

        n_unspliced_pos[i] = n_pos_us
        n_unspliced_neg[i] = n_neg_us
        n_spliced_pos[i] = n_pos_sp
        n_spliced_neg[i] = n_neg_sp

        # Fragment lengths for unspliced reads
        if n_unspliced > 0:
            if component == 1:
                fls = rng.normal(cfg.gdna_fl_mean, cfg.gdna_fl_std, n_unspliced)
            else:
                fls = rng.normal(cfg.rna_fl_mean, cfg.rna_fl_std, n_unspliced)
            fls = np.clip(fls, 50, 1000).astype(np.int32)
            fl_region_ids.extend([i] * n_unspliced)
            fl_frag_lens.extend(fls.tolist())

    # -- assemble DataFrames --
    region_counts = pd.DataFrame({
        "region_id": np.arange(n_total, dtype=np.int32),
        "n_unspliced_pos": n_unspliced_pos,
        "n_unspliced_neg": n_unspliced_neg,
        "n_spliced_pos": n_spliced_pos,
        "n_spliced_neg": n_spliced_neg,
    })

    region_df = pd.DataFrame({
        "region_id": np.arange(n_total, dtype=np.int32),
        "tx_pos": tx_pos,
        "tx_neg": tx_neg,
        "length": lengths,
    })

    fl_table = pd.DataFrame({
        "region_id": np.array(fl_region_ids, dtype=np.int32),
        "frag_len": np.array(fl_frag_lens, dtype=np.int32),
    })

    return {
        "region_counts": region_counts,
        "region_df": region_df,
        "fl_table": fl_table,
        "truth_labels": truth,
        "config": cfg,
    }


# ===================================================================
# Evaluation metrics
# ===================================================================


@dataclass
class ScenarioResult:
    """Metrics from a single scenario run."""

    label: str
    # Config
    strand_specificity: float
    true_pi: float
    n_rna: int
    n_gdna: int
    n_nrna: int
    rna_fl_mean: float
    gdna_fl_mean: float

    # EM convergence
    converged: bool
    n_iterations: int
    estimated_pi: float
    pi_error: float  # estimated - true

    # Classification accuracy (at threshold 0.5)
    gdna_precision: float
    gdna_recall: float
    gdna_f1: float
    rna_precision: float
    rna_recall: float
    rna_f1: float

    # Posterior quality
    mean_gdna_posterior: float   # mean γ for true-gDNA regions
    mean_rna_posterior: float    # mean γ for true-RNA regions
    mean_nrna_posterior: float   # mean γ for true-nRNA regions (if any)
    auroc: float                 # ROC-AUC for gDNA vs all-expressed
    separation: float            # mean_gdna_posterior - mean_rna_posterior

    # Kappa
    kappa_sym: float
    gdna_density: float
    expressed_density: float

    # Timing
    elapsed_s: float


def evaluate_scenario(
    result: GDNACalibration,
    truth: np.ndarray,
    cfg: ScenarioConfig,
    elapsed: float,
) -> ScenarioResult:
    """Compute classification metrics from calibration result."""
    gamma = result.region_posteriors

    # Binary classification: gDNA (truth==1) vs expressed (truth==0 or 2)
    is_gdna_true = truth == 1
    is_expressed_true = ~is_gdna_true

    pred_gdna = gamma >= 0.5
    pred_expressed = ~pred_gdna

    # Precision / recall / F1 for gDNA class
    tp_g = (pred_gdna & is_gdna_true).sum()
    fp_g = (pred_gdna & is_expressed_true).sum()
    fn_g = (pred_expressed & is_gdna_true).sum()
    gdna_prec = tp_g / max(tp_g + fp_g, 1)
    gdna_rec = tp_g / max(tp_g + fn_g, 1)
    gdna_f1 = 2 * gdna_prec * gdna_rec / max(gdna_prec + gdna_rec, 1e-12)

    # Precision / recall / F1 for RNA/expressed class
    tp_r = (pred_expressed & is_expressed_true).sum()
    fp_r = (pred_expressed & is_gdna_true).sum()
    fn_r = (pred_gdna & is_expressed_true).sum()
    rna_prec = tp_r / max(tp_r + fp_r, 1)
    rna_rec = tp_r / max(tp_r + fn_r, 1)
    rna_f1 = 2 * rna_prec * rna_rec / max(rna_prec + rna_rec, 1e-12)

    # Posterior means per group
    mean_gdna_post = float(gamma[truth == 1].mean()) if (truth == 1).any() else np.nan
    mean_rna_post = float(gamma[truth == 0].mean()) if (truth == 0).any() else np.nan
    mean_nrna_post = float(gamma[truth == 2].mean()) if (truth == 2).any() else np.nan

    # AUC-ROC (manual, no sklearn dependency)
    auroc = _compute_auroc(gamma, is_gdna_true)
    separation = mean_gdna_post - mean_rna_post if not np.isnan(mean_rna_post) else np.nan

    return ScenarioResult(
        label=cfg.label,
        strand_specificity=cfg.strand_specificity,
        true_pi=cfg.true_pi,
        n_rna=cfg.n_rna_regions,
        n_gdna=cfg.n_gdna_regions,
        n_nrna=cfg.n_nrna_regions,
        rna_fl_mean=cfg.rna_fl_mean,
        gdna_fl_mean=cfg.gdna_fl_mean,
        converged=result.converged,
        n_iterations=result.n_iterations,
        estimated_pi=result.mixing_proportion,
        pi_error=result.mixing_proportion - cfg.true_pi,
        gdna_precision=float(gdna_prec),
        gdna_recall=float(gdna_rec),
        gdna_f1=float(gdna_f1),
        rna_precision=float(rna_prec),
        rna_recall=float(rna_rec),
        rna_f1=float(rna_f1),
        mean_gdna_posterior=mean_gdna_post,
        mean_rna_posterior=mean_rna_post,
        mean_nrna_posterior=mean_nrna_post,
        auroc=auroc,
        separation=separation,
        kappa_sym=result.kappa_sym,
        gdna_density=result.gdna_density_global,
        expressed_density=result.expressed_density,
        elapsed_s=elapsed,
    )


def _compute_auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Wilcoxon-Mann-Whitney AUC estimate. O(n log n)."""
    pos = scores[labels]
    neg = scores[~labels]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    # U-statistic approach
    all_scores = np.concatenate([pos, neg])
    all_labels = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    order = np.argsort(all_scores)
    all_labels = all_labels[order]
    # Sum of ranks for positives (1-indexed)
    ranks = np.arange(1, len(all_scores) + 1, dtype=np.float64)
    pos_rank_sum = ranks[all_labels == 1].sum()
    n_pos = len(pos)
    n_neg = len(neg)
    u = pos_rank_sum - n_pos * (n_pos + 1) / 2
    return float(u / (n_pos * n_neg))


# ===================================================================
# Scenario grid builders
# ===================================================================


def build_default_grid() -> list[ScenarioConfig]:
    """Full factorial stress test grid.

    Axes:
        strand_specificity : [0.50, 0.60, 0.75, 0.90, 0.95, 0.99]
        gDNA level         : [low, medium, equal, dominant]
        FL overlap          : [identical, moderate, well_separated]
        nRNA                : [none, moderate, high]

    Total: 6 × 4 × 3 × 3 = 216 scenarios
    """
    ss_values = [0.50, 0.60, 0.75, 0.90, 0.95, 0.99]

    # gDNA levels: (n_rna, n_gdna, gdna_count_mean, label_suffix)
    gdna_levels = [
        (800, 200, 30.0, "low_gdna"),        # 20% gDNA, sparse
        (600, 400, 60.0, "med_gdna"),         # 40% gDNA, moderate counts
        (500, 500, 200.0, "equal_gdna"),      # 50% gDNA, equal to RNA median
        (300, 700, 300.0, "dominant_gdna"),   # 70% gDNA, higher than RNA
    ]

    # FL overlap: (gdna_fl_mean, gdna_fl_std, label_suffix)
    fl_configs = [
        (250.0, 50.0, "fl_identical"),         # same as RNA
        (300.0, 80.0, "fl_moderate"),          # partial overlap
        (400.0, 100.0, "fl_separated"),        # well separated
    ]

    # nRNA: (n_nrna, nrna_count_mean, nrna_splice_rate, label_suffix)
    nrna_configs = [
        (0, 0.0, 0.0, "no_nrna"),
        (150, 80.0, 0.05, "mod_nrna"),
        (300, 120.0, 0.02, "high_nrna"),
    ]

    grid = []
    for ss in ss_values:
        for (n_rna, n_gdna, gdna_cm, gdna_lbl) in gdna_levels:
            for (gdna_flm, gdna_fls, fl_lbl) in fl_configs:
                for (n_nrna, nrna_cm, nrna_sr, nrna_lbl) in nrna_configs:
                    label = f"ss{ss:.2f}_{gdna_lbl}_{fl_lbl}_{nrna_lbl}"
                    grid.append(ScenarioConfig(
                        n_rna_regions=n_rna,
                        n_gdna_regions=n_gdna,
                        n_nrna_regions=n_nrna,
                        strand_specificity=ss,
                        rna_fl_mean=250.0,
                        rna_fl_std=50.0,
                        gdna_fl_mean=gdna_flm,
                        gdna_fl_std=gdna_fls,
                        gdna_count_mean=gdna_cm,
                        nrna_count_mean=nrna_cm,
                        nrna_splice_rate=nrna_sr,
                        seed=42,
                        label=label,
                    ))
    return grid


def build_quick_grid() -> list[ScenarioConfig]:
    """Minimal 3×3×2×2 = 36 scenario grid for fast iteration."""
    ss_values = [0.50, 0.75, 0.95]
    gdna_levels = [
        (800, 200, 30.0, "low_gdna"),
        (500, 500, 200.0, "equal_gdna"),
        (300, 700, 300.0, "dominant_gdna"),
    ]
    fl_configs = [
        (250.0, 50.0, "fl_identical"),
        (400.0, 100.0, "fl_separated"),
    ]
    nrna_configs = [
        (0, 0.0, 0.0, "no_nrna"),
        (200, 100.0, 0.03, "mod_nrna"),
    ]
    grid = []
    for ss in ss_values:
        for (n_rna, n_gdna, gdna_cm, gdna_lbl) in gdna_levels:
            for (gdna_flm, gdna_fls, fl_lbl) in fl_configs:
                for (n_nrna, nrna_cm, nrna_sr, nrna_lbl) in nrna_configs:
                    label = f"ss{ss:.2f}_{gdna_lbl}_{fl_lbl}_{nrna_lbl}"
                    grid.append(ScenarioConfig(
                        n_rna_regions=n_rna,
                        n_gdna_regions=n_gdna,
                        n_nrna_regions=n_nrna,
                        strand_specificity=ss,
                        gdna_fl_mean=gdna_flm,
                        gdna_fl_std=gdna_fls,
                        gdna_count_mean=gdna_cm,
                        nrna_count_mean=nrna_cm,
                        nrna_splice_rate=nrna_sr,
                        seed=42,
                        label=label,
                    ))
    return grid


def build_edge_case_grid() -> list[ScenarioConfig]:
    """Targeted edge-case scenarios that stress boundary conditions."""
    cases = []

    # 1. All gDNA, no RNA — can the EM handle π → 1?
    cases.append(ScenarioConfig(
        n_rna_regions=0, n_gdna_regions=1000, strand_specificity=0.95,
        label="edge_all_gdna",
    ))

    # 2. All RNA, no gDNA — π → 0
    cases.append(ScenarioConfig(
        n_rna_regions=1000, n_gdna_regions=0, strand_specificity=0.95,
        label="edge_all_rna",
    ))

    # 3. Very few regions total
    cases.append(ScenarioConfig(
        n_rna_regions=10, n_gdna_regions=10, strand_specificity=0.95,
        label="edge_tiny_n",
    ))

    # 4. Extreme gDNA dominance (95%)
    cases.append(ScenarioConfig(
        n_rna_regions=50, n_gdna_regions=950, strand_specificity=0.95,
        gdna_count_mean=300.0, label="edge_95pct_gdna",
    ))

    # 5. Very low gDNA (2%)
    cases.append(ScenarioConfig(
        n_rna_regions=980, n_gdna_regions=20, strand_specificity=0.95,
        gdna_count_mean=20.0, label="edge_2pct_gdna",
    ))

    # 6. Completely unstranded library
    cases.append(ScenarioConfig(
        n_rna_regions=500, n_gdna_regions=500, strand_specificity=0.50,
        label="edge_unstranded",
    ))

    # 7. FL distributions perfectly identical, unstranded
    cases.append(ScenarioConfig(
        n_rna_regions=500, n_gdna_regions=500,
        strand_specificity=0.50,
        gdna_fl_mean=250.0, gdna_fl_std=50.0,
        gdna_count_mean=200.0,
        label="edge_no_signal",
    ))

    # 8. Only nRNA, no mRNA
    cases.append(ScenarioConfig(
        n_rna_regions=0, n_gdna_regions=500, n_nrna_regions=500,
        strand_specificity=0.95, nrna_splice_rate=0.02,
        nrna_count_mean=100.0,
        label="edge_nrna_only",
    ))

    # 9. Very high count RNA with very low count gDNA
    cases.append(ScenarioConfig(
        n_rna_regions=500, n_gdna_regions=500,
        strand_specificity=0.95,
        rna_count_mean=1000.0, rna_count_std=200.0,
        gdna_count_mean=5.0, gdna_count_std=2.0,
        label="edge_extreme_count_diff",
    ))

    # 10. Very low count RNA with very high count gDNA
    cases.append(ScenarioConfig(
        n_rna_regions=500, n_gdna_regions=500,
        strand_specificity=0.95,
        rna_count_mean=10.0, rna_count_std=5.0,
        gdna_count_mean=500.0, gdna_count_std=100.0,
        label="edge_gdna_high_count_rna_low",
    ))

    # 11. All ambiguous strand regions (tx_pos AND tx_neg both True)
    cases.append(ScenarioConfig(
        n_rna_regions=500, n_gdna_regions=500,
        strand_specificity=0.95,
        plus_strand_frac=1.0,  # all "+" — gene_strand = +1
        label="edge_all_plus_strand",
    ))

    # 12. Nearly perfect strand specificity
    cases.append(ScenarioConfig(
        n_rna_regions=500, n_gdna_regions=500,
        strand_specificity=0.999,
        label="edge_perfect_strand",
    ))

    # 13. Heavy nRNA contamination with identical FL
    cases.append(ScenarioConfig(
        n_rna_regions=200, n_gdna_regions=300, n_nrna_regions=500,
        strand_specificity=0.90,
        gdna_fl_mean=250.0, gdna_fl_std=50.0,
        nrna_count_mean=150.0, nrna_splice_rate=0.0,
        label="edge_heavy_nrna_identical_fl",
    ))

    # 14. Bimodal RNA count distribution
    cases.append(ScenarioConfig(
        n_rna_regions=500, n_gdna_regions=500,
        strand_specificity=0.95,
        rna_count_mean=50.0, rna_count_std=150.0,  # very high variance
        label="edge_bimodal_rna_counts",
    ))

    return cases


# ===================================================================
# Runner
# ===================================================================


def run_scenario(cfg: ScenarioConfig) -> ScenarioResult:
    """Generate data, run calibrate_gdna, evaluate."""
    data = generate_scenario(cfg)

    t0 = time.perf_counter()
    result = calibrate_gdna(
        data["region_counts"],
        data["fl_table"],
        data["region_df"],
        strand_specificity=cfg.strand_specificity,
        max_iterations=100,
        convergence_tol=1e-5,
        diagnostics=False,
    )
    elapsed = time.perf_counter() - t0

    return evaluate_scenario(result, data["truth_labels"], cfg, elapsed)


def run_grid(
    grid: list[ScenarioConfig],
    *,
    verbose: bool = True,
) -> pd.DataFrame:
    """Run all scenarios and return results DataFrame."""
    results = []
    n = len(grid)
    for i, cfg in enumerate(grid):
        if verbose:
            print(f"  [{i + 1}/{n}] {cfg.label} ...", end="", flush=True)
        try:
            r = run_scenario(cfg)
            results.append(asdict(r))
            if verbose:
                sym = "✓" if r.converged else "✗"
                print(
                    f" {sym}  π_err={r.pi_error:+.3f}"
                    f"  gDNA_F1={r.gdna_f1:.3f}"
                    f"  RNA_F1={r.rna_f1:.3f}"
                    f"  AUC={r.auroc:.3f}"
                    f"  sep={r.separation:.3f}"
                    f"  ({r.elapsed_s:.2f}s)"
                )
        except Exception as e:
            if verbose:
                print(f" FAILED: {e}")
            # Record failure
            results.append({
                "label": cfg.label,
                "strand_specificity": cfg.strand_specificity,
                "true_pi": cfg.true_pi,
                "n_rna": cfg.n_rna_regions,
                "n_gdna": cfg.n_gdna_regions,
                "n_nrna": cfg.n_nrna_regions,
                "rna_fl_mean": cfg.rna_fl_mean,
                "gdna_fl_mean": cfg.gdna_fl_mean,
                "converged": False,
                "error": str(e),
            })

    return pd.DataFrame(results)


# ===================================================================
# Analysis & reporting
# ===================================================================


def summarize_by_axis(df: pd.DataFrame) -> str:
    """Produce a text summary table grouped by each experimental axis."""
    lines = []
    lines.append("=" * 100)
    lines.append("CALIBRATION STRESS TEST SUMMARY")
    lines.append("=" * 100)
    lines.append("")

    # Overall stats
    if "converged" in df.columns:
        n_conv = df["converged"].sum()
        lines.append(f"Total scenarios: {len(df)}")
        lines.append(f"Converged: {n_conv}/{len(df)}")
        if "error" in df.columns:
            n_err = df["error"].notna().sum()
            lines.append(f"Errors: {n_err}")
        lines.append("")

    numeric = df.select_dtypes(include=[np.number])
    if "auroc" in numeric.columns:
        lines.append(f"Overall AUC:    mean={df['auroc'].mean():.3f}  "
                      f"min={df['auroc'].min():.3f}  max={df['auroc'].max():.3f}")
    if "gdna_f1" in numeric.columns:
        lines.append(f"Overall gDNA F1: mean={df['gdna_f1'].mean():.3f}  "
                      f"min={df['gdna_f1'].min():.3f}  max={df['gdna_f1'].max():.3f}")
    if "rna_f1" in numeric.columns:
        lines.append(f"Overall RNA F1:  mean={df['rna_f1'].mean():.3f}  "
                      f"min={df['rna_f1'].min():.3f}  max={df['rna_f1'].max():.3f}")
    if "pi_error" in numeric.columns:
        lines.append(f"Overall π error: mean={df['pi_error'].mean():+.4f}  "
                      f"std={df['pi_error'].std():.4f}")
    lines.append("")

    # By strand specificity
    if "strand_specificity" in df.columns:
        lines.append("-" * 80)
        lines.append("BY STRAND SPECIFICITY")
        lines.append("-" * 80)
        for ss in sorted(df["strand_specificity"].unique()):
            sub = df[df["strand_specificity"] == ss]
            cols = {}
            for c in ["auroc", "gdna_f1", "rna_f1", "pi_error", "separation"]:
                if c in sub.columns:
                    cols[c] = f"{sub[c].mean():.3f}"
            lines.append(f"  SS={ss:.2f} (n={len(sub)}): " +
                          "  ".join(f"{k}={v}" for k, v in cols.items()))
        lines.append("")

    # By gDNA level
    if "true_pi" in df.columns:
        lines.append("-" * 80)
        lines.append("BY gDNA CONTAMINATION LEVEL")
        lines.append("-" * 80)
        bins = [(0, 0.15, "pristine(<15%)"),
                (0.15, 0.35, "low(15-35%)"),
                (0.35, 0.55, "medium(35-55%)"),
                (0.55, 0.80, "high(55-80%)"),
                (0.80, 1.01, "dominant(>80%)")]
        for lo, hi, lbl in bins:
            sub = df[(df["true_pi"] >= lo) & (df["true_pi"] < hi)]
            if len(sub) == 0:
                continue
            cols = {}
            for c in ["auroc", "gdna_f1", "rna_f1", "pi_error"]:
                if c in sub.columns:
                    cols[c] = f"{sub[c].mean():.3f}"
            lines.append(f"  {lbl} (n={len(sub)}): " +
                          "  ".join(f"{k}={v}" for k, v in cols.items()))
        lines.append("")

    # By FL overlap
    if "gdna_fl_mean" in df.columns:
        lines.append("-" * 80)
        lines.append("BY FRAGMENT-LENGTH OVERLAP")
        lines.append("-" * 80)
        for flm in sorted(df["gdna_fl_mean"].unique()):
            sub = df[df["gdna_fl_mean"] == flm]
            cols = {}
            for c in ["auroc", "gdna_f1", "rna_f1"]:
                if c in sub.columns:
                    cols[c] = f"{sub[c].mean():.3f}"
            lines.append(f"  gDNA_FL_mean={flm:.0f} (n={len(sub)}): " +
                          "  ".join(f"{k}={v}" for k, v in cols.items()))
        lines.append("")

    # By nRNA presence
    if "n_nrna" in df.columns:
        lines.append("-" * 80)
        lines.append("BY NASCENT-RNA PRESENCE")
        lines.append("-" * 80)
        for nn in sorted(df["n_nrna"].unique()):
            sub = df[df["n_nrna"] == nn]
            cols = {}
            for c in ["auroc", "gdna_f1", "rna_f1", "mean_nrna_posterior"]:
                if c in sub.columns and sub[c].notna().any():
                    cols[c] = f"{sub[c].mean():.3f}"
            lines.append(f"  n_nrna={nn} (n={len(sub)}): " +
                          "  ".join(f"{k}={v}" for k, v in cols.items()))
        lines.append("")

    # Worst scenarios (bottom 10 by AUC)
    if "auroc" in df.columns:
        lines.append("-" * 80)
        lines.append("WORST 15 SCENARIOS (by AUC)")
        lines.append("-" * 80)
        worst = df.nsmallest(15, "auroc")
        for _, row in worst.iterrows():
            flag = ""
            if row.get("auroc", 1.0) < 0.6:
                flag = " *** CRITICAL"
            elif row.get("auroc", 1.0) < 0.75:
                flag = " ** WARNING"
            lines.append(
                f"  {row['label']}: "
                f"AUC={row.get('auroc', 'N/A'):.3f}  "
                f"gF1={row.get('gdna_f1', 'N/A'):.3f}  "
                f"rF1={row.get('rna_f1', 'N/A'):.3f}  "
                f"π_err={row.get('pi_error', 'N/A'):+.3f}  "
                f"sep={row.get('separation', 'N/A'):.3f}"
                f"{flag}"
            )
        lines.append("")

    # Best scenarios (top 10 by AUC)
    if "auroc" in df.columns:
        lines.append("-" * 80)
        lines.append("BEST 10 SCENARIOS (by AUC)")
        lines.append("-" * 80)
        best = df.nlargest(10, "auroc")
        for _, row in best.iterrows():
            lines.append(
                f"  {row['label']}: "
                f"AUC={row.get('auroc', 'N/A'):.3f}  "
                f"gF1={row.get('gdna_f1', 'N/A'):.3f}  "
                f"rF1={row.get('rna_f1', 'N/A'):.3f}  "
                f"π_err={row.get('pi_error', 'N/A'):+.3f}"
            )
        lines.append("")

    return "\n".join(lines)


def failure_analysis(df: pd.DataFrame, threshold: float = 0.75) -> str:
    """Detailed analysis of scenarios where AUC < threshold."""
    lines = []
    bad = df[df.get("auroc", pd.Series(dtype=float)) < threshold] if "auroc" in df.columns else df.head(0)
    if len(bad) == 0:
        lines.append(f"No scenarios below AUC threshold {threshold}.")
        return "\n".join(lines)

    lines.append(f"\n{'=' * 80}")
    lines.append(f"FAILURE ANALYSIS: {len(bad)} scenarios below AUC={threshold}")
    lines.append(f"{'=' * 80}")

    # Group by dominant failure mode
    for _, row in bad.iterrows():
        lines.append(f"\n  >>> {row['label']}")
        lines.append(f"      AUC={row.get('auroc', 'N/A'):.3f}  "
                      f"gF1={row.get('gdna_f1', 'N/A'):.3f}  "
                      f"rF1={row.get('rna_f1', 'N/A'):.3f}")
        lines.append(f"      π_true={row.get('true_pi', 'N/A'):.3f}  "
                      f"π_est={row.get('estimated_pi', 'N/A'):.3f}  "
                      f"π_err={row.get('pi_error', 'N/A'):+.3f}")
        lines.append(f"      mean_γ_gdna={row.get('mean_gdna_posterior', 'N/A'):.3f}  "
                      f"mean_γ_rna={row.get('mean_rna_posterior', 'N/A'):.3f}  "
                      f"sep={row.get('separation', 'N/A'):.3f}")
        lines.append(f"      SS={row.get('strand_specificity', 'N/A'):.2f}  "
                      f"gDNA_FL={row.get('gdna_fl_mean', 'N/A'):.0f}  "
                      f"n_nrna={row.get('n_nrna', 0)}")

        # Diagnose failure mode
        sep = row.get("separation", 0)
        pi_err = row.get("pi_error", 0)
        ss = row.get("strand_specificity", 0.5)
        gdna_fl = row.get("gdna_fl_mean", 250)

        modes = []
        if abs(sep) < 0.2:
            modes.append("POOR_SEPARATION: posteriors overlap heavily")
        if abs(pi_err) > 0.15:
            modes.append(f"PI_BIAS: π off by {pi_err:+.3f}")
        if ss <= 0.55 and gdna_fl == 250.0:
            modes.append("NO_SIGNAL: unstranded + identical FL → no discriminating features")
        if ss <= 0.55:
            modes.append("WEAK_STRAND: strand channel contributes ~zero")
        if gdna_fl == 250.0:
            modes.append("IDENTICAL_FL: no FL discrimination")
        if row.get("n_nrna", 0) > 0:
            modes.append(f"NRNA_CONFOUND: {row.get('n_nrna', 0)} nRNA regions")
        if not modes:
            modes.append("UNKNOWN")

        lines.append(f"      Failure modes: {'; '.join(modes)}")

    return "\n".join(lines)


# ===================================================================
# Main
# ===================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Calibration stress test for region partition + EM",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--quick", action="store_true",
        help="Run minimal 36-scenario grid",
    )
    parser.add_argument(
        "--edge", action="store_true",
        help="Run edge-case scenarios only",
    )
    parser.add_argument(
        "--full", action="store_true",
        help="Run full 216-scenario grid (default)",
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Run full grid + edge cases",
    )
    parser.add_argument(
        "--outdir", type=str, default=None,
        help="Output directory for results (default: print to stdout)",
    )
    parser.add_argument(
        "--auc-threshold", type=float, default=0.75,
        help="AUC threshold for failure analysis (default: 0.75)",
    )
    args = parser.parse_args()

    # Select grid
    grids: dict[str, list[ScenarioConfig]] = {}
    if args.all:
        grids["full"] = build_default_grid()
        grids["edge"] = build_edge_case_grid()
    elif args.edge:
        grids["edge"] = build_edge_case_grid()
    elif args.quick:
        grids["quick"] = build_quick_grid()
    else:
        grids["full"] = build_default_grid()

    all_results = []
    for name, grid in grids.items():
        print(f"\n{'=' * 60}")
        print(f"Running {name} grid: {len(grid)} scenarios")
        print(f"{'=' * 60}")
        df = run_grid(grid, verbose=True)
        df["grid"] = name
        all_results.append(df)

    combined = pd.concat(all_results, ignore_index=True)

    # Summary
    summary = summarize_by_axis(combined)
    print("\n" + summary)

    # Failure analysis
    analysis = failure_analysis(combined, threshold=args.auc_threshold)
    print(analysis)

    # Save results
    if args.outdir:
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        combined.to_csv(outdir / "stress_test_results.csv", index=False)
        (outdir / "summary.txt").write_text(summary + "\n" + analysis)

        # Save config as JSON
        configs = []
        for name, grid in grids.items():
            for cfg in grid:
                d = asdict(cfg)
                d["grid"] = name
                configs.append(d)
        (outdir / "configs.json").write_text(json.dumps(configs, indent=2))

        print(f"\nResults saved to {outdir}/")


if __name__ == "__main__":
    main()
