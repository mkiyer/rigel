"""rigel.calibration — Iterative seed-extend gDNA calibration.

Estimates global gDNA nuisance parameters (κ_sym strand symmetry
concentration, fragment-length distribution) from region-level evidence
produced by :mod:`rigel.region_evidence`.

Algorithm overview:

1. Compute per-region summary statistics from counts + annotation.
2. Select *seed* regions (high-confidence gDNA-only: zero spliced,
   strand-symmetric for stranded data, low-density for unstranded).
3. Estimate initial κ_sym and gDNA FL from seed.
4. Score all regions → per-region gDNA weights w_r ∈ [0, 1].
5. Re-estimate parameters using weights; iterate until convergence.

See ``docs/rigel_revision_impl/gdna_calibration_plan.md`` for full design.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .frag_length_model import FragmentLengthModel

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

#: Minimum unspliced fragment count for a region to be eligible.
MIN_COVERAGE: int = 10

#: z-score threshold for strand symmetry in seed selection.
#: 1.5 corresponds roughly to a two-sided binomial p ≈ 0.13.
SEED_SYMMETRY_ZSCORE: float = 1.5

#: Percentile threshold for coverage density in seed selection.
#: Applied among zero-spliced, symmetric, sufficient-coverage regions.
SEED_COVERAGE_PERCENTILE: float = 50.0

#: Iteration limits.
MAX_ITERATIONS: int = 20
CONVERGENCE_TOL: float = 0.01

#: κ bounds.
KAPPA_MIN: float = 1.0
KAPPA_MAX: float = 10_000.0
KAPPA_FALLBACK: float = 50.0

#: Minimum number of weighted regions for κ estimation.
KAPPA_MIN_OBS: int = 10

#: Small constant to avoid division by zero.
_EPS: float = 1e-12


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GDNACalibration:
    """Results from iterative gDNA calibration."""

    kappa_sym: float
    gdna_density: float
    gdna_fl_model: FragmentLengthModel
    region_weights: np.ndarray  # (n_regions,) float64 ∈ [0, 1]
    n_iterations: int
    converged: bool
    # Diagnostics (populated when diagnostics=True in calibrate_gdna)
    region_stats: dict[str, np.ndarray] | None = None
    seed_mask: np.ndarray | None = None
    kappa_history: list[float] | None = None
    density_history: list[float] | None = None


# ---------------------------------------------------------------------------
# Per-region summary statistics
# ---------------------------------------------------------------------------


def compute_region_stats(
    region_counts: pd.DataFrame,
    region_df: pd.DataFrame,
) -> dict[str, np.ndarray]:
    """Compute per-region summary statistics.

    Parameters
    ----------
    region_counts : pd.DataFrame
        Columns: ``n_unspliced_pos``, ``n_unspliced_neg``,
        ``n_spliced_pos``, ``n_spliced_neg`` (float32).
    region_df : pd.DataFrame
        Region metadata with ``tx_pos``, ``tx_neg``, ``length`` columns.

    Returns
    -------
    dict with keys: n_pos, n_neg, n_unspliced, n_spliced, n_total,
    strand_ratio, splice_rate, density, gene_strand, region_length.
    All arrays have shape ``(n_regions,)`` and dtype float64.
    """
    n_pos = region_counts["n_unspliced_pos"].values.astype(np.float64)
    n_neg = region_counts["n_unspliced_neg"].values.astype(np.float64)
    n_spliced = (
        region_counts["n_spliced_pos"].values
        + region_counts["n_spliced_neg"].values
    ).astype(np.float64)

    n_unspliced = n_pos + n_neg
    n_total = n_unspliced + n_spliced

    strand_ratio = np.full(len(n_pos), np.nan)
    has_unspliced = n_unspliced > 0
    strand_ratio[has_unspliced] = n_pos[has_unspliced] / n_unspliced[has_unspliced]

    splice_rate = np.zeros(len(n_pos), dtype=np.float64)
    has_total = n_total > 0
    splice_rate[has_total] = n_spliced[has_total] / n_total[has_total]

    region_length = region_df["length"].values.astype(np.float64)
    density = np.zeros(len(n_pos), dtype=np.float64)
    has_length = region_length > 0
    density[has_length] = n_total[has_length] / region_length[has_length]

    # Gene strand context: +1 only‑sense, -1 only‑antisense, 0 ambiguous.
    tx_pos = region_df["tx_pos"].values.astype(bool)
    tx_neg = region_df["tx_neg"].values.astype(bool)
    gene_strand = np.zeros(len(region_df), dtype=np.int8)
    gene_strand[tx_pos & ~tx_neg] = 1
    gene_strand[~tx_pos & tx_neg] = -1

    return {
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_unspliced": n_unspliced,
        "n_spliced": n_spliced,
        "n_total": n_total,
        "strand_ratio": strand_ratio,
        "splice_rate": splice_rate,
        "density": density,
        "gene_strand": gene_strand,
        "region_length": region_length,
    }


# ---------------------------------------------------------------------------
# Seed selection
# ---------------------------------------------------------------------------


def select_seed(
    stats: dict[str, np.ndarray],
    strand_specificity: float,
    *,
    min_coverage: int = MIN_COVERAGE,
    symmetry_zscore: float = SEED_SYMMETRY_ZSCORE,
    coverage_percentile: float = SEED_COVERAGE_PERCENTILE,
) -> np.ndarray:
    """Select high-confidence gDNA seed regions.

    Always applies ALL three criteria simultaneously:

    1. Zero spliced fragments + sufficient unspliced coverage.
    2. Strand symmetry (binomial z-test for p = 0.5).  When the
       library is unstranded (SS ≈ 0.5) this test is ~vacuous
       (nearly all regions pass), which is the correct behaviour.
    3. Low coverage density (bottom percentile among regions that
       pass criteria 1–2).  This is the primary discriminator
       for unstranded data and provides an additional safety
       margin for stranded data.

    No hard threshold on strand specificity is needed — the
    informativeness of each criterion emerges naturally from the
    data.

    Parameters
    ----------
    stats : dict
        Output of :func:`compute_region_stats`.
    strand_specificity : float
        Library strand specificity ∈ [0.5, 1.0].  Documented for
        API consistency; not used directly (symmetry test is
        model-free).

    Returns
    -------
    np.ndarray
        Boolean mask of shape ``(n_regions,)`` — True for seed regions.
    """
    n_unspliced = stats["n_unspliced"]
    n_spliced = stats["n_spliced"]
    n_pos = stats["n_pos"]
    density = stats["density"]

    # 1. Zero spliced + minimum coverage
    mask = (n_spliced == 0) & (n_unspliced >= min_coverage)

    # 2. Strand symmetry: z = |2k - n| / sqrt(n), require z < threshold
    z = np.abs(2.0 * n_pos - n_unspliced) / np.sqrt(
        np.maximum(n_unspliced, 1.0)
    )
    mask = mask & (z < symmetry_zscore)

    # 3. Low coverage density: bottom percentile among qualifying
    qualifying = mask & (density > 0)
    if qualifying.any():
        threshold = np.percentile(
            density[qualifying], coverage_percentile
        )
        mask = mask & (density <= threshold)

    return mask


# ---------------------------------------------------------------------------
# κ_sym estimation (MoM, forced mean = 0.5)
# ---------------------------------------------------------------------------


def estimate_kappa_sym(
    stats: dict[str, np.ndarray],
    weights: np.ndarray,
    *,
    kappa_min: float = KAPPA_MIN,
    kappa_max: float = KAPPA_MAX,
    kappa_fallback: float = KAPPA_FALLBACK,
    kappa_min_obs: int = KAPPA_MIN_OBS,
) -> float:
    """Weighted Method-of-Moments estimation of κ_sym.

    Models strand ratios as Beta(κ/2, κ/2): mean forced to 0.5,
    variance = 0.25 / (κ + 1).

    Corrects for within-region sampling variance:
    ``observed_var = true_var + 0.25 * E[1/n]``.

    Parameters
    ----------
    stats : dict
        Output of :func:`compute_region_stats`.
    weights : np.ndarray
        Per-region weights ∈ [0, 1].

    Returns
    -------
    float
        Estimated κ_sym, clamped to [kappa_min, kappa_max].
    """
    n_unspliced = stats["n_unspliced"]
    strand_ratio = stats["strand_ratio"]

    # Only use regions with sufficient data and positive weight
    valid = (n_unspliced >= 2) & (weights > 0) & np.isfinite(strand_ratio)
    n_valid = valid.sum()

    if n_valid < kappa_min_obs:
        return kappa_fallback

    w = weights[valid]
    p = strand_ratio[valid]
    n = n_unspliced[valid]
    w_sum = w.sum()

    if w_sum < _EPS:
        return kappa_fallback

    # Weighted variance of strand ratios around 0.5
    observed_var = np.sum(w * (p - 0.5) ** 2) / w_sum

    # Sampling variance correction: E_w[0.25 / n]
    sampling_var = 0.25 * np.sum(w / n) / w_sum

    true_var = max(observed_var - sampling_var, _EPS)

    # κ from Var[Beta(κ/2, κ/2)] = 0.25 / (κ + 1)
    kappa = 0.25 / true_var - 1.0

    return float(np.clip(kappa, kappa_min, kappa_max))


# ---------------------------------------------------------------------------
# gDNA density estimation
# ---------------------------------------------------------------------------


def estimate_gdna_density(
    stats: dict[str, np.ndarray],
    weights: np.ndarray,
) -> float:
    """Weighted average gDNA density (fragments / bp) from region weights.

    Returns 0.0 if no valid regions.
    """
    n_unspliced = stats["n_unspliced"]
    region_length = stats["region_length"]

    valid = (weights > 0) & (region_length > 0) & (n_unspliced > 0)
    if not valid.any():
        return 0.0

    w = weights[valid]
    # Weight by both gDNA weight and region length (longer regions
    # contribute more to the density estimate).
    weighted_frags = np.sum(w * n_unspliced[valid])
    weighted_bp = np.sum(w * region_length[valid])

    if weighted_bp < _EPS:
        return 0.0

    return float(weighted_frags / weighted_bp)


# ---------------------------------------------------------------------------
# Region scoring
# ---------------------------------------------------------------------------


def score_regions(
    stats: dict[str, np.ndarray],
    kappa_sym: float,
    strand_specificity: float,
    gdna_density: float = 0.0,
) -> np.ndarray:
    """Score each region for gDNA purity.

    Returns per-region weights w_r ∈ [0, 1] indicating estimated
    proportion of unspliced fragments attributable to gDNA.

    Combines **all** available signals additively in log-likelihood
    ratio space:

    .. math::

        \\text{llr} = \\text{llr}_{\\text{strand}} + \\text{llr}_{\\text{density}}

    * **Strand**: normal-approximation LLR comparing symmetric gDNA
      model vs SS-informed RNA model.  Contributes ≈ 0 when
      SS ≈ 0.5 (unstranded), dominant when SS → 1.
    * **Density**: Poisson z-score comparing observed coverage to
      expected gDNA background.  Contributes 0 when *gdna_density*
      is unknown (≤ 0).
    * **Splice**: applied multiplicatively as ``(1 − splice_rate)``
      since spliced reads are definitively RNA.

    No hard threshold on strand specificity — the informativeness of
    each signal emerges naturally from the data and the math.

    Parameters
    ----------
    stats : dict
        Output of :func:`compute_region_stats`.
    kappa_sym : float
        Current κ_sym estimate.
    strand_specificity : float
        Library strand specificity ∈ [0.5, 1.0].
    gdna_density : float
        Expected gDNA density (frags/bp).  Pass 0.0 to disable the
        density component.
    """
    n_unspliced = stats["n_unspliced"]
    n_pos = stats["n_pos"]
    splice_rate = stats["splice_rate"]
    gene_strand = stats["gene_strand"]
    region_length = stats["region_length"]

    n_regions = len(n_unspliced)
    has_data = n_unspliced >= 2

    # --- Strand LLR (always computed; ≈0 when SS ≈ 0.5) ---
    llr_strand = _compute_strand_llr(
        n_pos, n_unspliced, has_data, kappa_sym,
        strand_specificity, gene_strand, n_regions,
    )

    # --- Density LLR (always computed; 0 when density unknown) ---
    llr_density = _compute_density_llr(
        n_unspliced, region_length, has_data,
        gdna_density, n_regions,
    )

    # --- Combined LLR → gDNA weight for unspliced fragments ---
    llr = np.clip(llr_strand + llr_density, -500, 500)

    w_unspliced = np.full(n_regions, 0.5, dtype=np.float64)
    w_unspliced[has_data] = 1.0 / (1.0 + np.exp(-llr[has_data]))

    # Splice adjustment: spliced reads are definitively RNA.
    weights = (1.0 - splice_rate) * w_unspliced

    return np.clip(weights, 0.0, 1.0)


def _compute_strand_llr(
    n_pos: np.ndarray,
    n_unspliced: np.ndarray,
    has_data: np.ndarray,
    kappa_sym: float,
    strand_specificity: float,
    gene_strand: np.ndarray,
    n_regions: int,
) -> np.ndarray:
    """Strand log-likelihood ratio (gDNA vs RNA).

    gDNA model: strand_ratio ~ N(0.5, var_gdna)
        where var_gdna = 0.25/(κ+1) + 0.25/n
    RNA model:  strand_ratio ~ N(p_rna, var_rna)
        where p_rna = SS for + gene, 1-SS for - gene, 0.5 for ambiguous
        and var_rna = p_rna*(1-p_rna)/n

    When SS ≈ 0.5, both models predict the same strand ratio → LLR ≈ 0.
    """
    llr = np.zeros(n_regions, dtype=np.float64)

    p = np.full(n_regions, 0.5)
    p[has_data] = n_pos[has_data] / n_unspliced[has_data]

    # gDNA variance: overdispersion + sampling
    var_gdna = np.full(n_regions, 1.0)
    var_gdna[has_data] = (
        0.25 / (kappa_sym + 1.0) + 0.25 / n_unspliced[has_data]
    )

    # RNA expected strand ratio and variance
    p_rna = np.full(n_regions, 0.5)
    p_rna[gene_strand == 1] = strand_specificity
    p_rna[gene_strand == -1] = 1.0 - strand_specificity

    var_rna = np.full(n_regions, 1.0)
    var_rna[has_data] = np.maximum(
        p_rna[has_data] * (1.0 - p_rna[has_data]) / n_unspliced[has_data],
        _EPS,
    )

    # Normal log-likelihood: -0.5 * ((x - mu)^2 / var + ln(var))
    ll_gdna = -0.5 * (
        (p[has_data] - 0.5) ** 2 / var_gdna[has_data]
        + np.log(var_gdna[has_data])
    )
    ll_rna = -0.5 * (
        (p[has_data] - p_rna[has_data]) ** 2 / var_rna[has_data]
        + np.log(var_rna[has_data])
    )

    llr[has_data] = ll_gdna - ll_rna
    return llr


def _compute_density_llr(
    n_unspliced: np.ndarray,
    region_length: np.ndarray,
    has_data: np.ndarray,
    gdna_density: float,
    n_regions: int,
) -> np.ndarray:
    """Density log-likelihood ratio via Poisson z-score.

    Under gDNA the expected unspliced count is
    ``μ = gdna_density × region_length``.  The z-score
    ``z = (k − μ) / √max(μ, 1)`` measures excess/deficit:

    * z > 0 (excess coverage) → LLR < 0 → evidence for RNA.
    * z < 0 (deficit) → LLR > 0 → evidence for gDNA.
    * z ≈ 0 (matches expected) → LLR ≈ 0 → no information.

    The z-score magnitude grows as O(√n) while the strand LLR
    grows as O(n), so strand naturally dominates for large samples
    where it is most reliable.

    Returns all zeros when *gdna_density* ≤ 0 (density unknown).
    """
    llr = np.zeros(n_regions, dtype=np.float64)

    if gdna_density <= 0:
        return llr

    mu = gdna_density * region_length
    valid = has_data & (mu > 0)

    z = (
        (n_unspliced[valid] - mu[valid])
        / np.sqrt(np.maximum(mu[valid], 1.0))
    )
    llr[valid] = -z  # excess → negative LLR (RNA)

    return llr


# ---------------------------------------------------------------------------
# gDNA fragment-length model
# ---------------------------------------------------------------------------


def build_gdna_fl_model(
    fl_region_ids: np.ndarray,
    fl_frag_lens: np.ndarray,
    region_weights: np.ndarray,
    max_fl: int = 1000,
) -> FragmentLengthModel:
    """Build gDNA fragment-length model weighted by region gDNA weights.

    Parameters
    ----------
    fl_region_ids : np.ndarray, shape (n_obs,)
        Region ID for each FL observation.
    fl_frag_lens : np.ndarray, shape (n_obs,)
        Fragment length for each observation.
    region_weights : np.ndarray, shape (n_regions,)
        Per-region gDNA weights ∈ [0, 1].
    max_fl : int
        Maximum fragment length for the model.

    Returns
    -------
    FragmentLengthModel
        Finalized gDNA FL model.
    """
    model = FragmentLengthModel(max_size=max_fl)

    if len(fl_region_ids) == 0:
        model.finalize()
        return model

    # Look up per-observation weight from region weights
    fl_weights = region_weights[fl_region_ids]

    # Filter valid observations
    valid = (fl_weights > 0) & (fl_frag_lens > 0) & (fl_frag_lens <= max_fl)
    if not valid.any():
        model.finalize()
        return model

    # Vectorized histogram via bincount
    valid_fl = fl_frag_lens[valid].astype(np.intp)
    valid_w = fl_weights[valid].astype(np.float64)

    weighted_hist = np.bincount(valid_fl, weights=valid_w, minlength=max_fl + 1)

    # Feed into model
    for fl_val in range(1, max_fl + 1):
        if weighted_hist[fl_val] > 0:
            model.observe(fl_val, weight=float(weighted_hist[fl_val]))

    model.finalize()
    return model


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def calibrate_gdna(
    region_counts: pd.DataFrame,
    fl_table: pd.DataFrame,
    region_df: pd.DataFrame,
    strand_specificity: float,
    *,
    min_coverage: int = MIN_COVERAGE,
    max_iterations: int = MAX_ITERATIONS,
    convergence_tol: float = CONVERGENCE_TOL,
    seed_symmetry_zscore: float = SEED_SYMMETRY_ZSCORE,
    seed_coverage_percentile: float = SEED_COVERAGE_PERCENTILE,
    kappa_min: float = KAPPA_MIN,
    kappa_max: float = KAPPA_MAX,
    kappa_fallback: float = KAPPA_FALLBACK,
    diagnostics: bool = False,
) -> GDNACalibration:
    """Iteratively estimate gDNA parameters via seed-extend.

    All three signals (strand symmetry, coverage density, spliced-read
    evidence) are combined into a single additive log-likelihood ratio
    at every iteration — no hard threshold on strand specificity.

    Parameters
    ----------
    region_counts : pd.DataFrame
        Per-region fragment counts (from :func:`count_region_evidence`).
    fl_table : pd.DataFrame
        Fragment-length observations with ``region_id`` and ``frag_len``
        columns (from :func:`count_region_evidence`).
    region_df : pd.DataFrame
        Region metadata (from ``index.region_df``).
    strand_specificity : float
        Library strand specificity ∈ [0.5, 1.0].
    diagnostics : bool
        When True, populate ``region_stats``, ``seed_mask``,
        ``kappa_history``, and ``density_history`` on the result.

    Returns
    -------
    GDNACalibration
    """
    # ------------------------------------------------------------------
    # Phase 1: summary statistics
    # ------------------------------------------------------------------
    stats = compute_region_stats(region_counts, region_df)

    # ------------------------------------------------------------------
    # Phase 2: seed selection (all criteria applied simultaneously)
    # ------------------------------------------------------------------
    seed_mask = select_seed(
        stats,
        strand_specificity,
        min_coverage=min_coverage,
        symmetry_zscore=seed_symmetry_zscore,
        coverage_percentile=seed_coverage_percentile,
    )
    seed_weights = seed_mask.astype(np.float64)

    logger.info(
        "gDNA calibration: %d seed regions (of %d total)",
        int(seed_mask.sum()),
        len(seed_mask),
    )

    # ------------------------------------------------------------------
    # Phase 3: initial parameter estimates from seed
    # ------------------------------------------------------------------
    kappa_sym = estimate_kappa_sym(
        stats,
        seed_weights,
        kappa_min=kappa_min,
        kappa_max=kappa_max,
        kappa_fallback=kappa_fallback,
    )
    gdna_density = estimate_gdna_density(stats, seed_weights)

    logger.info(
        "gDNA calibration: initial κ_sym=%.1f, density=%.2e",
        kappa_sym,
        gdna_density,
    )

    # Diagnostics: per-iteration history
    kappa_hist: list[float] = [kappa_sym] if diagnostics else []
    density_hist: list[float] = [gdna_density] if diagnostics else []

    # ------------------------------------------------------------------
    # Phase 4: iterate score → re-estimate
    # ------------------------------------------------------------------
    fl_region_ids = fl_table["region_id"].values.astype(np.intp)
    fl_frag_lens = fl_table["frag_len"].values.astype(np.intp)

    weights = seed_weights.copy()
    converged = False
    n_iter = 0

    for iteration in range(max_iterations):
        n_iter = iteration + 1

        # Score all regions (unified: strand + density + splice)
        weights = score_regions(
            stats,
            kappa_sym,
            strand_specificity,
            gdna_density,
        )

        # Re-estimate κ
        new_kappa = estimate_kappa_sym(
            stats,
            weights,
            kappa_min=kappa_min,
            kappa_max=kappa_max,
            kappa_fallback=kappa_fallback,
        )

        # Re-estimate density
        gdna_density = estimate_gdna_density(stats, weights)

        if diagnostics:
            kappa_hist.append(new_kappa)
            density_hist.append(gdna_density)

        # Convergence check
        rel_change = abs(new_kappa - kappa_sym) / max(kappa_sym, 1.0)

        logger.debug(
            "gDNA calibration iter %d: κ=%.1f → %.1f (Δ=%.4f), "
            "density=%.2e, mean_weight=%.3f",
            n_iter,
            kappa_sym,
            new_kappa,
            rel_change,
            gdna_density,
            weights.mean(),
        )

        kappa_sym = new_kappa

        if rel_change < convergence_tol:
            converged = True
            break

    logger.info(
        "gDNA calibration: %s after %d iterations, κ_sym=%.1f",
        "converged" if converged else "did not converge",
        n_iter,
        kappa_sym,
    )

    # ------------------------------------------------------------------
    # Phase 5: final FL model
    # ------------------------------------------------------------------
    gdna_fl = build_gdna_fl_model(fl_region_ids, fl_frag_lens, weights)

    return GDNACalibration(
        kappa_sym=kappa_sym,
        gdna_density=gdna_density,
        gdna_fl_model=gdna_fl,
        region_weights=weights,
        n_iterations=n_iter,
        converged=converged,
        region_stats=stats if diagnostics else None,
        seed_mask=seed_mask if diagnostics else None,
        kappa_history=kappa_hist if diagnostics else None,
        density_history=density_hist if diagnostics else None,
    )
