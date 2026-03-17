"""rigel.calibration — Aggregate-First EM for gDNA deconvolution.

Classifies genomic regions as *not expressed* (gDNA only) or *expressed*
(gDNA + RNA) using a two-component mixture model.

Three convergent signals drive the classification:

1. **Density** — Gaussian model on log-density (log(n/L + ε)).
   gDNA and RNA components each have (μ, σ²) estimated via EM.
2. **Strand balance** — Binomial model.  gDNA: Binom(k | n, 0.5);
   RNA: Binom(k | n, SS).  LLR vanishes when SS = 0.5 (no magic
   thresholds) and scales naturally with sample size.
3. **Fragment length** — FL_E frozen from spliced reads (gold
   standard); FL_G built iteratively from γ-weighted fragments.

Hard constraint: ``n_spliced > 0 → γ = 0`` (definitively expressed).
Zero-count regions are excluded from estimation; their posterior
defaults to the prior π.

Outputs:

* Per-region posteriors γ_r ∈ [0, 1] — P(not expressed | data).
* Global and per-chromosome gDNA density.
* Strand symmetry concentration κ (post-hoc diagnostic).
* gDNA fragment-length model.
* Mixing proportion π.

See ``docs/rigel_revision_plan/calibration_v3_aggregate_em.md``.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .frag_length_model import FragmentLengthModel

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Small constant to avoid division by zero.
# ---------------------------------------------------------------------------
_EPS: float = 1e-12

# ---------------------------------------------------------------------------
# Initialization thresholds.
GDNA_INIT_DENSITY_PERCENTILE: float = 0.10
GDNA_INIT_MIN_REGIONS: int = 100

# ---------------------------------------------------------------------------
# Vectorized log-gamma (for Beta-Binomial MLE)
# ---------------------------------------------------------------------------
_lgamma_vec = np.frompyfunc(math.lgamma, 1, 1)


def _vec_lgamma(x: np.ndarray) -> np.ndarray:
    """Vectorized log-gamma via ``math.lgamma``."""
    return _lgamma_vec(x).astype(np.float64)


# ---------------------------------------------------------------------------
# 1-D golden-section maximiser (no scipy)
# ---------------------------------------------------------------------------


def _golden_section_max(
    f,
    a: float,
    b: float,
    tol: float = 1e-4,
    max_iter: int = 100,
) -> float:
    """Argmax of *f* on [*a*, *b*] via golden-section search."""
    gr = (math.sqrt(5.0) + 1.0) / 2.0
    c = b - (b - a) / gr
    d = a + (b - a) / gr
    fc, fd = f(c), f(d)
    for _ in range(max_iter):
        if abs(b - a) < tol:
            break
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - (b - a) / gr
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + (b - a) / gr
            fd = f(d)
    return (a + b) / 2.0


def _beta_binom_loglik(
    kappa: float,
    k: np.ndarray,
    n: np.ndarray,
    w: np.ndarray,
) -> float:
    """γ-weighted log-likelihood of the symmetric Beta-Binomial.

    BetaBinom(k | n, κ/2, κ/2), summed over observations with weights *w*.
    Terms independent of κ are dropped.
    """
    alpha = kappa / 2.0
    per_obs = (
        _vec_lgamma(k + alpha)
        + _vec_lgamma(n - k + alpha)
        - _vec_lgamma(n + kappa)
        + math.lgamma(kappa)
        - 2.0 * math.lgamma(alpha)
    )
    return float(np.sum(w * per_obs))


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GDNACalibration:
    """Results from EM two-component gDNA deconvolution."""

    # Per-region posteriors: γ_r = P(not expressed | data)
    region_posteriors: np.ndarray  # (n_regions,) float64 ∈ [0, 1]

    # Global gDNA density (frags / bp)
    gdna_density_global: float

    # Per-reference-sequence gDNA density
    gdna_density_per_ref: dict[str, float]

    # Strand symmetry concentration: Beta(κ/2, κ/2)
    kappa_sym: float

    # gDNA fragment-length model
    gdna_fl_model: FragmentLengthModel

    # Mixing proportion: π = P(not expressed)
    mixing_proportion: float

    # Mean density of expressed regions (frags / bp)
    expressed_density: float

    # Convergence metadata
    n_iterations: int
    converged: bool

    # Diagnostics (populated when diagnostics=True)
    region_stats: dict[str, np.ndarray] | None = None
    iteration_history: list[dict] | None = None
    init_diagnostics: dict | None = None
    log_density: np.ndarray | None = None
    sense_frac: np.ndarray | None = None
    eligible: np.ndarray | None = None
    density_bin_edges: np.ndarray | None = None

    # --- Backward-compatible aliases ---
    @property
    def gdna_density(self) -> float:
        """Alias for ``gdna_density_global``."""
        return self.gdna_density_global

    @property
    def region_weights(self) -> np.ndarray:
        """Alias for ``region_posteriors``."""
        return self.region_posteriors


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
    strand_ratio, splice_rate, density, gene_strand, region_length,
    tx_pos, tx_neg, exon_pos, exon_neg, ref.
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

    # Gene strand context: +1 only-sense, -1 only-antisense, 0 ambiguous.
    tx_pos = region_df["tx_pos"].values.astype(bool)
    tx_neg = region_df["tx_neg"].values.astype(bool)
    gene_strand = np.zeros(len(region_df), dtype=np.int8)
    gene_strand[tx_pos & ~tx_neg] = 1
    gene_strand[~tx_pos & tx_neg] = -1

    # Annotation flags.
    exon_pos = (
        region_df["exon_pos"].values.astype(bool)
        if "exon_pos" in region_df.columns
        else np.zeros(len(region_df), dtype=bool)
    )
    exon_neg = (
        region_df["exon_neg"].values.astype(bool)
        if "exon_neg" in region_df.columns
        else np.zeros(len(region_df), dtype=bool)
    )
    ref = (
        region_df["ref"].values
        if "ref" in region_df.columns
        else np.full(len(region_df), "unknown", dtype=object)
    )

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
        "tx_pos": tx_pos,
        "tx_neg": tx_neg,
        "exon_pos": exon_pos,
        "exon_neg": exon_neg,
        "ref": ref,
    }


# ---------------------------------------------------------------------------
# Empirical strand histogram
# ---------------------------------------------------------------------------


def compute_sense_fraction(
    stats: dict[str, np.ndarray],
) -> np.ndarray:
    """Convert per-region strand_ratio to gene-strand-normalised sense fraction.

    R1-antisense convention (dUTP / TruSeq Stranded):

    * ``gene_strand == +1``: RNA reads land on − strand, so
      ``sense_frac = 1 − strand_ratio``.
    * ``gene_strand == -1``: RNA reads land on + strand, so
      ``sense_frac = strand_ratio``.
    * ``gene_strand == 0``: ambiguous — set to NaN (no strand signal).

    For gDNA, sense_frac is ~0.5 regardless of gene_strand.
    For RNA, sense_frac ≈ SS (close to 1.0 for stranded libraries).

    Returns an array of shape ``(n_regions,)`` with NaN where strand
    signal is unavailable (gene_strand == 0 or n_unspliced < 2).
    """
    strand_ratio = stats["strand_ratio"]
    gene_strand = stats["gene_strand"]
    n_unspliced = stats["n_unspliced"]

    sf = np.full(len(strand_ratio), np.nan, dtype=np.float64)

    valid = np.isfinite(strand_ratio) & (n_unspliced >= 2)
    plus = valid & (gene_strand == 1)
    minus = valid & (gene_strand == -1)

    sf[plus] = 1.0 - strand_ratio[plus]
    sf[minus] = strand_ratio[minus]

    return sf


def _compute_strand_llr_binomial(
    stats: dict[str, np.ndarray],
    strand_specificity: float,
    n_regions: int,
) -> np.ndarray:
    """Strand LLR using simple Binomial model.

    gDNA: Binom(k_sense | n, 0.5)  — gDNA is strand-symmetric.
    RNA:  Binom(k_sense | n, SS)   — RNA follows strand specificity.

    LLR = k·log(0.5/SS) + (n−k)·log(0.5/(1−SS)).

    When SS = 0.5 every term is zero — no magic threshold needed.
    Signal scales naturally with sample size: low-count regions
    contribute minimal LLR.

    Regions with ``gene_strand == 0`` or ``n_unspliced < 2`` get LLR = 0.
    """
    llr = np.zeros(n_regions, dtype=np.float64)

    n_unspliced = stats["n_unspliced"]
    n_pos = stats["n_pos"]
    gene_strand = stats["gene_strand"]

    valid = (n_unspliced >= 2) & (gene_strand != 0)
    if not valid.any():
        return llr

    # Clip SS away from 0 and 1 for numerical safety
    ss = float(np.clip(strand_specificity, _EPS, 1.0 - _EPS))

    # k_sense: number of gene-sense reads per region
    k_sense = np.empty(n_regions, dtype=np.float64)
    k_sense[gene_strand == 1] = (
        n_unspliced[gene_strand == 1] - n_pos[gene_strand == 1]
    )
    k_sense[gene_strand == -1] = n_pos[gene_strand == -1]

    k = k_sense[valid]
    n = n_unspliced[valid]

    log_ratio_sense = np.log(0.5 / ss)
    log_ratio_anti = np.log(0.5 / (1.0 - ss))

    llr[valid] = k * log_ratio_sense + (n - k) * log_ratio_anti
    return llr


# ---------------------------------------------------------------------------
# κ estimation (Beta-Binomial MLE, forced mean = 0.5)
# ---------------------------------------------------------------------------


def estimate_kappa_sym(
    stats: dict[str, np.ndarray],
    weights: np.ndarray,
) -> float | None:
    """Beta-Binomial MLE estimation of κ_sym.

    Maximises the γ-weighted log-likelihood of the symmetric
    Beta-Binomial model BetaBinom(k | n, κ/2, κ/2) via golden-section
    search.  Unlike Method-of-Moments, this remains accurate at low
    per-region counts because the Beta-Binomial generative model
    inherently accounts for binomial sampling variance.

    Returns ``None`` when fewer than 3 valid (weighted, n ≥ 2) regions
    are available.  Callers should disable the strand signal when κ is
    ``None``.
    """
    n_unspliced = stats["n_unspliced"]
    n_pos = stats["n_pos"]
    strand_ratio = stats["strand_ratio"]

    valid = (n_unspliced >= 2) & (weights > 0) & np.isfinite(strand_ratio)
    if valid.sum() < 3:
        return None

    k = n_pos[valid]
    n = n_unspliced[valid]
    w = weights[valid]

    kappa = _golden_section_max(
        lambda kap: _beta_binom_loglik(kap, k, n, w),
        0.01, 500.0,
    )
    return max(kappa, 0.0)


# ---------------------------------------------------------------------------
# Per-chromosome gDNA density
# ---------------------------------------------------------------------------


def _compute_per_ref_density(
    stats: dict[str, np.ndarray],
    posteriors: np.ndarray,
) -> dict[str, float]:
    """Compute γ-weighted gDNA density per reference sequence."""
    ref = stats["ref"]
    n_total = stats["n_total"]
    region_length = stats["region_length"]

    result: dict[str, float] = {}
    for ref_name in np.unique(ref):
        mask = (ref == ref_name) & (posteriors > 0) & (region_length > 0)
        if not mask.any():
            result[str(ref_name)] = 0.0
            continue
        w = posteriors[mask]
        frags = np.sum(w * n_total[mask])
        bp = np.sum(w * region_length[mask])
        result[str(ref_name)] = float(frags / bp) if bp > _EPS else 0.0

    return result


# ---------------------------------------------------------------------------
# gDNA fragment-length model
# ---------------------------------------------------------------------------


def build_gdna_fl_model(
    fl_region_ids: np.ndarray,
    fl_frag_lens: np.ndarray,
    region_weights: np.ndarray,
    max_fl: int = 1000,
) -> FragmentLengthModel:
    """Build gDNA fragment-length model weighted by region posteriors.

    Parameters
    ----------
    fl_region_ids : np.ndarray, shape (n_obs,)
        Region ID for each FL observation.
    fl_frag_lens : np.ndarray, shape (n_obs,)
        Fragment length for each observation.
    region_weights : np.ndarray, shape (n_regions,)
        Per-region posteriors γ_r ∈ [0, 1].
    max_fl : int
        Maximum fragment length for the model.
    """
    model = FragmentLengthModel(max_size=max_fl)

    if len(fl_region_ids) == 0:
        model.finalize()
        return model

    fl_weights = region_weights[fl_region_ids]
    valid = (fl_weights > 0) & (fl_frag_lens > 0) & (fl_frag_lens <= max_fl)
    if not valid.any():
        model.finalize()
        return model

    valid_fl = fl_frag_lens[valid].astype(np.intp)
    valid_w = fl_weights[valid].astype(np.float64)

    weighted_hist = np.bincount(valid_fl, weights=valid_w, minlength=max_fl + 1)

    for fl_val in range(1, max_fl + 1):
        if weighted_hist[fl_val] > 0:
            model.observe(fl_val, weight=float(weighted_hist[fl_val]))

    model.finalize()
    return model


# ---------------------------------------------------------------------------
# Density histogram (log-density channel)
# ---------------------------------------------------------------------------


def compute_log_density(
    stats: dict[str, np.ndarray],
    eligible: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Compute log-density = log(n/L + ε) for all regions.

    Parameters
    ----------
    stats : dict
        Region summary statistics from :func:`compute_region_stats`.
    eligible : np.ndarray, bool
        Mask for regions participating in the EM.

    Returns
    -------
    (log_d, epsilon) where log_d has shape (n_regions,) and epsilon is
    the global pseudocount.  Ineligible regions get log_d = 0.
    """
    n_total = stats["n_total"]
    region_length = stats["region_length"]
    n_regions = len(n_total)

    # Global ε = 1 / median(L) over eligible regions
    med_L = float(np.median(region_length[eligible]))
    epsilon = 1.0 / max(med_L, 1.0)

    log_d = np.zeros(n_regions, dtype=np.float64)
    valid = eligible & (region_length > 0)
    log_d[valid] = np.log(
        n_total[valid] / region_length[valid] + epsilon
    )
    return log_d, epsilon


def _compute_density_llr_gaussian(
    log_d: np.ndarray,
    mu_g: float,
    var_g: float,
    mu_r: float,
    var_r: float,
    eligible: np.ndarray,
    n_regions: int,
) -> np.ndarray:
    """Density LLR from Gaussian model on log-density.

    LLR = log N(x | μ_G, σ_G) − log N(x | μ_R, σ_R) for each region.
    Ineligible regions get LLR = 0.
    """
    llr = np.zeros(n_regions, dtype=np.float64)
    if not eligible.any():
        return llr

    x = log_d[eligible]
    llr[eligible] = -0.5 * (
        (x - mu_g) ** 2 / var_g
        - (x - mu_r) ** 2 / var_r
        + np.log(var_g)
        - np.log(var_r)
    )
    return llr


# ---------------------------------------------------------------------------
# Initialization: two-phase seed partition
# ---------------------------------------------------------------------------


def _seed_initial_partition(
    stats: dict[str, np.ndarray],
    log_d: np.ndarray,
    sense_frac: np.ndarray,
    strand_specificity: float,
    eligible: np.ndarray,
    *,
    density_percentile: float = GDNA_INIT_DENSITY_PERCENTILE,
    min_gdna_regions: int = GDNA_INIT_MIN_REGIONS,
) -> tuple[np.ndarray, float, dict]:
    """Two-phase initialization for the EM.

    Phase 1: expressed seed from ``n_spliced > 0``.
    Phase 2: gDNA seed from density eCDF + optional strand filter.

    Returns (gamma, pi_init, diagnostics_dict).
    """
    n_total = stats["n_total"]
    n_spliced = stats["n_spliced"]
    gene_strand = stats["gene_strand"]
    n_regions = len(n_total)

    gamma = np.full(n_regions, 0.5, dtype=np.float64)
    diag: dict = {}

    # Phase 1: expressed seed (definitively expressed)
    expressed_seed = (n_spliced > 0) & eligible
    gamma[expressed_seed] = 0.0

    # Build expressed log-density eCDF
    expressed_log_d = log_d[expressed_seed]
    n_expressed = expressed_seed.sum()

    if n_expressed < 2:
        # Fallback: not enough spliced regions
        # Use top 50th percentile of eligible regions
        log_d_eligible = log_d[eligible]
        if len(log_d_eligible) > 0:
            threshold = float(np.median(log_d_eligible))
            expressed_seed = eligible & (log_d >= threshold)
            gamma[expressed_seed] = 0.0
            expressed_log_d = log_d[expressed_seed]
            n_expressed = expressed_seed.sum()
        diag["expressed_fallback"] = True

    # Phase 2: gDNA seed from density eCDF
    # For each unspliced-only eligible region, compute its quantile in
    # the expressed log-density distribution.
    unspliced_only = eligible & (n_spliced == 0)
    gdna_seed = np.zeros(n_regions, dtype=bool)

    if n_expressed > 0 and unspliced_only.any():
        expressed_sorted = np.sort(expressed_log_d)
        # eCDF quantile for each unspliced-only region
        p_expressed = np.searchsorted(
            expressed_sorted, log_d[unspliced_only], side="right",
        ) / len(expressed_sorted)

        # Regions below the density percentile threshold
        below_threshold = p_expressed < density_percentile
        candidate_indices = np.where(unspliced_only)[0]
        gdna_candidates = candidate_indices[below_threshold]

        # Minimum seed size guarantee
        if len(gdna_candidates) < min_gdna_regions:
            # Sort all unspliced-only by p_expressed, take bottom N
            sorted_order = np.argsort(p_expressed)
            n_take = min(min_gdna_regions, len(sorted_order))
            gdna_candidates = candidate_indices[sorted_order[:n_take]]

        gdna_seed[gdna_candidates] = True

    # Supplementary strand filter (when stranded library)
    if strand_specificity > 0.55:
        strand_symmetric = (
            unspliced_only
            & np.isfinite(sense_frac)
            & (np.abs(sense_frac - 0.5) < 0.1)
            & (gene_strand != 0)
        )
        gdna_seed |= strand_symmetric

    gamma[gdna_seed] = 1.0

    # Ensure expressed and gdna seeds don't overlap
    gamma[expressed_seed] = 0.0

    # π from seed counts
    n_eligible = eligible.sum()
    n_gdna = gdna_seed.sum()
    if n_eligible > 0:
        pi_init = float(max(n_gdna, 1)) / float(n_eligible)
    else:
        pi_init = 0.5

    # Zero-count regions get the prior
    gamma[~eligible] = pi_init

    diag["n_expressed_seed"] = int(n_expressed)
    diag["n_gdna_seed"] = int(n_gdna)
    diag["pi_init"] = pi_init
    diag["pristine_sample"] = (
        pi_init < 0.02 or n_gdna < min_gdna_regions // 2
    )

    return gamma, pi_init, diag


# ---------------------------------------------------------------------------
# E-step
# ---------------------------------------------------------------------------


def _e_step(
    stats: dict[str, np.ndarray],
    pi: float,
    log_d: np.ndarray,
    eligible: np.ndarray,
    strand_specificity: float,
    mu_g: float,
    var_g: float,
    mu_r: float,
    var_r: float,
    fl_region_ids: np.ndarray | None = None,
    fl_frag_lens: np.ndarray | None = None,
    gdna_fl_model: FragmentLengthModel | None = None,
    rna_fl_model: FragmentLengthModel | None = None,
) -> np.ndarray:
    """Compute per-region posterior γ_r = P(not expressed | data).

    Hard constraints:
    - n_spliced > 0 → γ = 0 (definitively expressed).
    - not eligible (n_total = 0 or L = 0) → γ = π (prior).

    Density channel uses Gaussian model (μ, σ²) per component.
    Strand channel uses Binomial(k | n, p) — no histograms.
    """
    n_spliced = stats["n_spliced"]
    n_regions = len(stats["n_total"])

    gamma = np.full(n_regions, pi, dtype=np.float64)

    # Hard constraint: spliced → expressed
    has_splice = n_spliced > 0
    gamma[has_splice] = 0.0

    # Soft classification for eligible unspliced regions
    soft = eligible & ~has_splice
    if not soft.any():
        return gamma

    # --- Log prior odds ---
    pi_safe = np.clip(pi, _EPS, 1.0 - _EPS)
    log_prior_odds = np.log(pi_safe / (1.0 - pi_safe))

    # --- Density LLR (Gaussian) ---
    llr_density = _compute_density_llr_gaussian(
        log_d, mu_g, var_g, mu_r, var_r, soft, n_regions,
    )

    # --- Strand LLR (Binomial) ---
    llr_strand = _compute_strand_llr_binomial(
        stats, strand_specificity, n_regions,
    )

    # --- Fragment length LLR (optional) ---
    llr_fl = np.zeros(n_regions, dtype=np.float64)
    if (
        gdna_fl_model is not None
        and rna_fl_model is not None
        and fl_region_ids is not None
        and fl_frag_lens is not None
        and len(fl_region_ids) > 0
        and gdna_fl_model._total_weight > 0
        and rna_fl_model._total_weight > 0
    ):
        llr_fl = _compute_fl_llr(
            fl_region_ids, fl_frag_lens,
            gdna_fl_model, rna_fl_model, n_regions,
        )

    # --- Combined → posterior ---
    log_odds = log_prior_odds + llr_density + llr_strand + llr_fl
    log_odds_soft = np.clip(log_odds[soft], -500, 500)
    gamma[soft] = 1.0 / (1.0 + np.exp(-log_odds_soft))

    return gamma


def _compute_fl_llr(
    fl_region_ids: np.ndarray,
    fl_frag_lens: np.ndarray,
    gdna_fl_model: FragmentLengthModel,
    rna_fl_model: FragmentLengthModel,
    n_regions: int,
) -> np.ndarray:
    """Sum of per-fragment FL log-likelihood ratios per region.

    Uses **shape-normalized** densities: both histograms are converted
    to proportions and smoothed with a shared Dirichlet prior
    (α = 1/M per bin, i.e. one total pseudo-observation).  This
    ensures only distributional *shape* differences drive the LLR —
    not the total amount of data behind each model, which would
    otherwise introduce a systematic per-fragment bias when the two
    components have different total weights.
    """
    llr = np.zeros(n_regions, dtype=np.float64)
    if len(fl_region_ids) == 0:
        return llr

    max_size = gdna_fl_model.max_size
    M = max_size + 1

    # Normalize histograms to proportions → shape only.
    T_g = max(gdna_fl_model._total_weight, 1.0)
    T_r = max(rna_fl_model._total_weight, 1.0)

    # Shared symmetric Dirichlet smoothing: α = 1/M per bin
    # (one total pseudo-observation spread uniformly).
    alpha = 1.0 / M
    q_g = gdna_fl_model.counts / T_g + alpha
    q_r = rna_fl_model.counts / T_r + alpha

    # Pre-compute log-ratio lookup (both q vectors sum to 1 + Mα = 2,
    # so the shared denominator cancels in the ratio).
    log_ratio = np.log(q_g) - np.log(q_r)

    for i in range(len(fl_region_ids)):
        rid = int(fl_region_ids[i])
        fl = int(fl_frag_lens[i])
        if rid < n_regions and fl > 0:
            idx = min(fl, max_size)
            llr[rid] += log_ratio[idx]

    return llr


# ---------------------------------------------------------------------------
# M-step
# ---------------------------------------------------------------------------


def _m_step(
    stats: dict[str, np.ndarray],
    gamma: np.ndarray,
    log_d: np.ndarray,
    eligible: np.ndarray,
) -> tuple[float, float, float, float, float, float, float]:
    """Update parameters from posteriors.

    Returns ``(pi, lambda_G, lambda_E, mu_g, var_g, mu_r, var_r)``.

    Density Gaussian parameters use **region-level** weights (γ alone,
    not γ × n_total) because each region is one observation of
    log-density.
    """
    n_total = stats["n_total"]
    region_length = stats["region_length"]

    n_eligible = eligible.sum()
    if n_eligible == 0:
        return 0.5, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0

    # --- Mixing proportion (eligible regions only) ---
    pi = float(gamma[eligible].sum() / n_eligible)
    pi = np.clip(pi, _EPS, 1.0 - _EPS)

    # --- Scalar density summaries (for output / diagnostics) ---
    w_g = gamma[eligible]
    w_e = 1.0 - gamma[eligible]
    nt = n_total[eligible]
    rl = region_length[eligible]

    g_frags = np.sum(w_g * nt)
    g_bp = np.sum(w_g * rl)
    lambda_G = float(g_frags / g_bp) if g_bp > _EPS else 0.0

    e_frags = np.sum(w_e * nt)
    e_bp = np.sum(w_e * rl)
    lambda_E = float(e_frags / e_bp) if e_bp > _EPS else 0.0

    # --- Density Gaussian parameters (region-level weights) ---
    ld = log_d[eligible]
    mu_overall = float(np.mean(ld))
    var_overall = max(float(np.var(ld)), _EPS)

    w_g_sum = float(w_g.sum())
    if w_g_sum >= 1.0:
        mu_g = float(np.sum(w_g * ld) / w_g_sum)
        var_g = max(float(np.sum(w_g * (ld - mu_g) ** 2) / w_g_sum), _EPS)
    else:
        mu_g, var_g = mu_overall, var_overall

    w_e_sum = float(w_e.sum())
    if w_e_sum >= 1.0:
        mu_r = float(np.sum(w_e * ld) / w_e_sum)
        var_r = max(float(np.sum(w_e * (ld - mu_r) ** 2) / w_e_sum), _EPS)
    else:
        mu_r, var_r = mu_overall, var_overall

    return float(pi), lambda_G, lambda_E, mu_g, var_g, mu_r, var_r


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def calibrate_gdna(
    region_counts: pd.DataFrame,
    fl_table: pd.DataFrame,
    region_df: pd.DataFrame,
    strand_specificity: float,
    *,
    max_iterations: int = 50,
    convergence_tol: float = 1e-4,
    diagnostics: bool = False,
    density_percentile: float = GDNA_INIT_DENSITY_PERCENTILE,
    min_gdna_regions: int = GDNA_INIT_MIN_REGIONS,
) -> GDNACalibration:
    """Estimate gDNA parameters via Aggregate-First EM.

    Classifies regions as *not expressed* (gDNA only) or *expressed*
    (gDNA + RNA) using Gaussian density and Binomial strand models.

    Parameters
    ----------
    region_counts : pd.DataFrame
        Per-region fragment counts (from :func:`count_region_evidence`).
    fl_table : pd.DataFrame
        Fragment-length observations with ``region_id`` and ``frag_len``
        columns.
    region_df : pd.DataFrame
        Region metadata with ``tx_pos``, ``tx_neg``, ``length`` columns.
    strand_specificity : float
        Library strand specificity ∈ [0.5, 1.0].
    max_iterations : int
        Maximum EM iterations.
    convergence_tol : float
        Convergence tolerance for |Δπ|.
    diagnostics : bool
        Populate ``region_stats`` and ``iteration_history``.
    density_percentile : float
        Quantile threshold for gDNA seed initialization.
    min_gdna_regions : int
        Minimum number of gDNA seed regions.
    """
    # ------------------------------------------------------------------
    # Phase 1: summary statistics
    # ------------------------------------------------------------------
    stats = compute_region_stats(region_counts, region_df)

    n_regions = len(stats["n_total"])
    fl_region_ids = fl_table["region_id"].values.astype(np.intp)
    fl_frag_lens = fl_table["frag_len"].values.astype(np.intp)

    # Eligible regions: n_total > 0 and L > 0
    eligible = (stats["n_total"] > 0) & (stats["region_length"] > 0)

    # ------------------------------------------------------------------
    # Bailout: no data at all
    # ------------------------------------------------------------------
    if not eligible.any():
        logger.warning("gDNA calibration: no eligible regions")
        empty_fl = build_gdna_fl_model(fl_region_ids, fl_frag_lens,
                                        np.zeros(n_regions))
        return GDNACalibration(
            region_posteriors=np.ones(n_regions, dtype=np.float64),
            gdna_density_global=0.0,
            gdna_density_per_ref={},
            kappa_sym=0.0,
            gdna_fl_model=empty_fl,
            mixing_proportion=1.0,
            expressed_density=0.0,
            n_iterations=0,
            converged=True,
            region_stats=stats if diagnostics else None,
            iteration_history=[] if diagnostics else None,
            init_diagnostics=None,
            log_density=None,
            sense_frac=None,
            eligible=eligible if diagnostics else None,
            density_bin_edges=None,
        )

    # ------------------------------------------------------------------
    # Phase 2: preprocessing
    # ------------------------------------------------------------------
    sense_frac = compute_sense_fraction(stats)
    log_d, epsilon = compute_log_density(stats, eligible)

    # ------------------------------------------------------------------
    # Phase 3: initialization (two-phase seed partition)
    # ------------------------------------------------------------------
    gamma, pi, init_diag = _seed_initial_partition(
        stats, log_d, sense_frac, strand_specificity, eligible,
        density_percentile=density_percentile,
        min_gdna_regions=min_gdna_regions,
    )

    # Compute initial Gaussian density parameters from seed weights.
    ld_elig = log_d[eligible]
    w_g_init = gamma[eligible]
    w_e_init = 1.0 - gamma[eligible]

    w_g_sum = max(float(w_g_init.sum()), _EPS)
    mu_g = float(np.sum(w_g_init * ld_elig) / w_g_sum)
    var_g = max(float(np.sum(w_g_init * (ld_elig - mu_g) ** 2) / w_g_sum), _EPS)

    w_e_sum = max(float(w_e_init.sum()), _EPS)
    mu_r = float(np.sum(w_e_init * ld_elig) / w_e_sum)
    var_r = max(float(np.sum(w_e_init * (ld_elig - mu_r) ** 2) / w_e_sum), _EPS)

    # Build initial FL models from seeds
    gdna_fl = build_gdna_fl_model(fl_region_ids, fl_frag_lens, gamma)
    rna_fl = build_gdna_fl_model(fl_region_ids, fl_frag_lens, 1.0 - gamma)

    # Compute pi_soft: mixing proportion over unspliced ("soft") regions.
    # The E-step prior should reflect P(gDNA | unspliced), not P(gDNA)
    # globally, because hard-constraint spliced regions (gamma=0) are
    # already fixed and should not dilute the prior for soft regions.
    has_splice = stats["n_spliced"] > 0
    soft = eligible & ~has_splice
    n_soft = int(soft.sum())
    n_eligible = int(eligible.sum())
    if n_soft > 0:
        pi_soft = float(np.clip(
            pi * n_eligible / n_soft, _EPS, 1.0 - _EPS,
        ))
    else:
        pi_soft = pi

    logger.info(
        "gDNA calibration init: π=%.3f, π_soft=%.3f, n_expressed=%d, "
        "n_gdna=%d, n_soft=%d, pristine=%s, ε=%.2e, "
        "μ_G=%.2f, σ²_G=%.2f, μ_R=%.2f, σ²_R=%.2f",
        pi, pi_soft, init_diag["n_expressed_seed"], init_diag["n_gdna_seed"],
        n_soft, init_diag["pristine_sample"], epsilon,
        mu_g, var_g, mu_r, var_r,
    )

    # ------------------------------------------------------------------
    # Phase 4: EM iteration
    # ------------------------------------------------------------------
    history: list[dict] = []
    converged = False
    n_iter = 0
    lambda_G = 0.0
    lambda_E = 0.0

    for iteration in range(max_iterations):
        n_iter = iteration + 1

        # E-step: use pi_soft as prior (conditioned on unspliced)
        gamma = _e_step(
            stats, pi_soft, log_d, eligible, strand_specificity,
            mu_g, var_g, mu_r, var_r,
            fl_region_ids=fl_region_ids,
            fl_frag_lens=fl_frag_lens,
            gdna_fl_model=gdna_fl,
            rna_fl_model=rna_fl,
        )

        # M-step: update pi, lambdas, and Gaussian density params
        new_pi, new_lG, new_lE, mu_g, var_g, mu_r, var_r = _m_step(
            stats, gamma, log_d, eligible,
        )

        # Update FL models
        gdna_fl = build_gdna_fl_model(
            fl_region_ids, fl_frag_lens, gamma,
        )
        rna_fl = build_gdna_fl_model(
            fl_region_ids, fl_frag_lens, 1.0 - gamma,
        )

        # Update pi_soft for next E-step
        if n_soft > 0:
            new_pi_soft = float(np.clip(
                gamma[soft].sum() / n_soft, _EPS, 1.0 - _EPS,
            ))
        else:
            new_pi_soft = new_pi

        # Convergence check on pi_soft (the E-step prior)
        delta_pi = abs(new_pi_soft - pi_soft)

        if diagnostics:
            history.append({
                "pi": new_pi,
                "pi_soft": new_pi_soft,
                "lambda_G": new_lG,
                "lambda_E": new_lE,
                "mean_gamma": float(gamma[eligible].mean()),
                "delta_pi": delta_pi,
                "gamma": gamma.copy(),
                "mu_g": mu_g,
                "var_g": var_g,
                "mu_r": mu_r,
                "var_r": var_r,
            })

        logger.debug(
            "gDNA calibration iter %d: π=%.4f, π_soft=%.4f, λ_G=%.2e, "
            "λ_E=%.2e, mean_γ=%.3f, Δπ_soft=%.6f, "
            "μ_G=%.2f, σ²_G=%.2f, μ_R=%.2f, σ²_R=%.2f",
            n_iter, new_pi, new_pi_soft, new_lG, new_lE,
            gamma[eligible].mean(), delta_pi,
            mu_g, var_g, mu_r, var_r,
        )

        pi = new_pi
        pi_soft = new_pi_soft
        lambda_G = new_lG
        lambda_E = new_lE

        if delta_pi < convergence_tol:
            converged = True
            break

    logger.info(
        "gDNA calibration: %s after %d iterations. "
        "π=%.3f, λ_G=%.2e, λ_E=%.2e",
        "converged" if converged else "did not converge",
        n_iter, pi, lambda_G, lambda_E,
    )

    # ------------------------------------------------------------------
    # Phase 5: final posteriors and outputs
    # ------------------------------------------------------------------
    final_gamma = _e_step(
        stats, pi_soft, log_d, eligible, strand_specificity,
        mu_g, var_g, mu_r, var_r,
        fl_region_ids=fl_region_ids,
        fl_frag_lens=fl_frag_lens,
        gdna_fl_model=gdna_fl,
        rna_fl_model=rna_fl,
    )

    per_ref = _compute_per_ref_density(stats, final_gamma)
    final_fl = build_gdna_fl_model(fl_region_ids, fl_frag_lens, final_gamma)
    kappa_posthoc = estimate_kappa_sym(stats, final_gamma)

    return GDNACalibration(
        region_posteriors=final_gamma,
        gdna_density_global=lambda_G,
        gdna_density_per_ref=per_ref,
        kappa_sym=kappa_posthoc if kappa_posthoc is not None else 0.0,
        gdna_fl_model=final_fl,
        mixing_proportion=pi,
        expressed_density=lambda_E,
        n_iterations=n_iter,
        converged=converged,
        region_stats=stats if diagnostics else None,
        iteration_history=history if diagnostics else None,
        init_diagnostics=init_diag if diagnostics else None,
        log_density=log_d if diagnostics else None,
        sense_frac=sense_frac if diagnostics else None,
        eligible=eligible if diagnostics else None,
        density_bin_edges=None,
    )
