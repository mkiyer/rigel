"""rigel.priors — Empirical Bayes prior computation for the EM solver.

All prior functions that run between the scan pass and the locus EM
loop live here:

- ``compute_global_gdna_density`` — global gDNA density via strand
  correction.
- ``estimate_kappa`` — Method-of-Moments Beta concentration κ.
- ``compute_nrna_frac_priors`` — 3-tier density + strand EB
  hierarchical shrinkage for per-nRNA nascent fraction priors.

These are pure functions operating on numpy arrays and the
``AbundanceEstimator`` accumulator state.  They mutate
``estimator.nrna_frac_alpha`` / ``nrna_frac_beta`` in-place.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from .config import EMConfig

if TYPE_CHECKING:
    from .estimator import AbundanceEstimator

logger = logging.getLogger(__name__)


# ======================================================================
# Constants
# ======================================================================

#: Minimum denominator for the strand weight ``2s − 1`` below which
#: the strand component of the hybrid estimator is zeroed out.
#: Also imported by locus.py.
STRAND_DENOM_EPS: float = 0.01

#: Default mean fragment length when no observations are available.
DEFAULT_MEAN_FRAG: float = 200.0

# Default kappa bounds from EMConfig (single source of truth).
_EM_DEFAULTS = EMConfig()


# ======================================================================
# Hybrid density + strand nrna_frac estimator
# ======================================================================

def _compute_hybrid_nrna_frac_vec(
    exon_sense: np.ndarray,
    exon_anti: np.ndarray,
    intron_sense: np.ndarray,
    intron_anti: np.ndarray,
    L_exonic: np.ndarray,
    L_intronic: np.ndarray,
    strand_spec: float,
    gdna_density: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute nrna_frac estimates and evidence from per-unit count arrays.

    Uses the unified density + strand hybrid model:

    * **Density subtraction** isolates RNA rates from spatial zone
      densities (exonic vs intronic vs intergenic).
    * **Strand subtraction** isolates RNA from gDNA via the strand
      specificity.
    * The two are combined through inverse-variance weighting
      ``W = (2s − 1)²``.

    Returns
    -------
    nrna_frac : ndarray
        Nascent fraction estimates, same shape as inputs.
    evidence : ndarray
        Total exonic + intronic fragment count per unit.
    """
    exon_total = exon_sense + exon_anti
    intron_total = intron_sense + intron_anti
    evidence = exon_total + intron_total

    # --- Density estimator ---
    with np.errstate(divide="ignore", invalid="ignore"):
        D_exon = np.where(L_exonic > 0, exon_total / L_exonic, 0.0)
        D_intron = np.where(
            L_intronic > 0, intron_total / L_intronic, gdna_density,
        )
    den_nrna = np.maximum(0.0, D_intron - gdna_density)
    den_mrna = np.maximum(0.0, D_exon - D_intron)

    # --- Strand estimator ---
    denom = 2.0 * strand_spec - 1.0
    inv_var_weight = denom ** 2 if denom > STRAND_DENOM_EPS else 0.0

    if inv_var_weight > 0:
        exon_rna = np.maximum(0.0, (exon_sense - exon_anti) / denom)
        str_exon_density = exon_rna / np.maximum(1.0, L_exonic)
        intron_rna = np.maximum(0.0, (intron_sense - intron_anti) / denom)
        with np.errstate(divide="ignore", invalid="ignore"):
            str_intron_density = np.where(
                L_intronic > 0, intron_rna / L_intronic, 0.0,
            )
        str_nrna = str_intron_density
        str_mrna = np.maximum(0.0, str_exon_density - str_intron_density)
    else:
        str_nrna = den_nrna
        str_mrna = den_mrna

    # --- Weighted combination ---
    final_nrna = inv_var_weight * str_nrna + (1.0 - inv_var_weight) * den_nrna
    final_mrna = inv_var_weight * str_mrna + (1.0 - inv_var_weight) * den_mrna

    # --- nrna_frac = nRNA / (mRNA + nRNA) ---
    total_rate = final_mrna + final_nrna
    with np.errstate(divide="ignore", invalid="ignore"):
        nrna_frac = np.where(total_rate > 0, final_nrna / total_rate, 0.0)

    return nrna_frac, evidence


# ======================================================================
# Global gDNA density
# ======================================================================

def compute_global_gdna_density(
    estimator: AbundanceEstimator,
    strand_specificity: float,
) -> float:
    """Estimate global gDNA density (fragments / bp) via strand correction.

    Uses the same strand-based gDNA isolation as the EB prior system:
    from global unspliced sense/antisense totals, compute the gDNA
    fraction, then convert to a per-bp density using total exonic span.

    Returns 0.0 when strand specificity is too low for reliable
    correction.
    """
    from .locus import compute_gdna_density_from_strand

    total_sense = float(estimator.transcript_unspliced_sense.sum())
    total_anti = float(estimator.transcript_unspliced_antisense.sum())

    if estimator._exonic_lengths is not None:
        total_exonic_bp = float(estimator._exonic_lengths.sum())
    else:
        total_exonic_bp = 0.0

    if total_exonic_bp > 0:
        return compute_gdna_density_from_strand(
            total_sense, total_anti, total_exonic_bp, strand_specificity,
        )
    return 0.0


# ======================================================================
# Method-of-Moments κ estimation
# ======================================================================

def estimate_kappa(
    nrna_frac_array: np.ndarray,
    evidence_array: np.ndarray,
    min_evidence: float = 30.0,
    *,
    kappa_min: float = _EM_DEFAULTS.nrna_frac_kappa_min,
    kappa_max: float = _EM_DEFAULTS.nrna_frac_kappa_max,
    kappa_fallback: float = _EM_DEFAULTS.nrna_frac_kappa_fallback,
    kappa_min_obs: int = _EM_DEFAULTS.nrna_frac_kappa_min_obs,
) -> float:
    """Estimate Beta concentration κ via Method of Moments.

    Selects features whose evidence ≥ *min_evidence*, computes the
    mean (μ) and variance (σ²) of their raw nrna_frac values, and
    solves ``κ = μ(1 − μ)/σ² − 1``, clamped to [kappa_min, kappa_max].

    Returns *kappa_fallback* if fewer than *kappa_min_obs* features
    pass the evidence filter.
    """
    valid = evidence_array >= min_evidence
    if valid.sum() < kappa_min_obs:
        return kappa_fallback

    confident_nrna_fracs = nrna_frac_array[valid]
    mu = float(np.mean(confident_nrna_fracs))
    sigma2 = float(np.var(confident_nrna_fracs))

    if sigma2 <= 0.0:
        return kappa_max  # zero variance → maximum precision

    kappa = mu * (1.0 - mu) / sigma2 - 1.0
    return float(np.clip(kappa, kappa_min, kappa_max))


# ======================================================================
# Group aggregation helper
# ======================================================================

def _aggregate_nrna_frac_by_group(
    group_ids: np.ndarray,
    n_groups: int,
    exon_sense: np.ndarray,
    exon_anti: np.ndarray,
    intron_sense: np.ndarray,
    intron_anti: np.ndarray,
    L_exonic: np.ndarray,
    L_intronic: np.ndarray,
    strand_spec: float,
    gdna_density: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Aggregate counts by group ID, compute nrna_frac, map back per-transcript.

    Returns
    -------
    (group_nrna_frac, group_evidence, per_t_nrna_frac, per_t_evidence)
    """
    g_exon_s = np.zeros(n_groups, dtype=np.float64)
    g_exon_a = np.zeros(n_groups, dtype=np.float64)
    g_int_s = np.zeros(n_groups, dtype=np.float64)
    g_int_a = np.zeros(n_groups, dtype=np.float64)
    g_L_ex = np.zeros(n_groups, dtype=np.float64)
    g_L_in = np.zeros(n_groups, dtype=np.float64)

    valid = group_ids >= 0
    if valid.any():
        ids = group_ids[valid]
        np.add.at(g_exon_s, ids, exon_sense[valid])
        np.add.at(g_exon_a, ids, exon_anti[valid])
        np.add.at(g_int_s, ids, intron_sense[valid])
        np.add.at(g_int_a, ids, intron_anti[valid])
        np.add.at(g_L_ex, ids, L_exonic[valid])
        np.add.at(g_L_in, ids, L_intronic[valid])

    g_nrna_frac, g_evidence = _compute_hybrid_nrna_frac_vec(
        g_exon_s, g_exon_a, g_int_s, g_int_a,
        g_L_ex, g_L_in, strand_spec, gdna_density,
    )

    n_t = len(group_ids)
    per_t_nrna_frac = np.zeros(n_t, dtype=np.float64)
    per_t_evidence = np.zeros(n_t, dtype=np.float64)
    if valid.any():
        per_t_nrna_frac[valid] = g_nrna_frac[ids]
        per_t_evidence[valid] = g_evidence[ids]

    return g_nrna_frac, g_evidence, per_t_nrna_frac, per_t_evidence


# ======================================================================
# 3-tier nrna_frac prior computation
# ======================================================================

def compute_nrna_frac_priors(
    estimator: AbundanceEstimator,
    nrna_strands: np.ndarray,
    nrna_spans: np.ndarray,
    strand_specificity: float,
    gdna_density: float,
    kappa_global: float | None = None,
    kappa_locus: float | None = None,
    kappa_nrna: float = _EM_DEFAULTS.nrna_frac_kappa_nrna,
    mom_min_evidence_global: float = 50.0,
    mom_min_evidence_locus: float = 20.0,
    kappa_min: float = _EM_DEFAULTS.nrna_frac_kappa_min,
    kappa_max: float = _EM_DEFAULTS.nrna_frac_kappa_max,
    kappa_fallback: float = _EM_DEFAULTS.nrna_frac_kappa_fallback,
    kappa_min_obs: int = _EM_DEFAULTS.nrna_frac_kappa_min_obs,
) -> None:
    """Compute per-nRNA nrna_frac (nascent fraction) Beta priors.

    Uses a 3-tier **density + strand hybrid** model with smooth
    Empirical-Bayes hierarchical shrinkage:

    1. **Locus-strand nrna_frac** ← shrink toward empirical global
       nrna_frac with strength *kappa_global*.
    2. **Per-nRNA nrna_frac** ← shrink toward shrunk locus-strand
       nrna_frac with strength *kappa_locus*.

    The final Beta prior uses a constant pseudo-count *kappa_nrna*
    (default 5.0) so the EM solver is guided but not overwhelmed.

    Mutates ``estimator.nrna_frac_alpha`` and ``estimator.nrna_frac_beta``
    in-place.
    """
    n_nrna = estimator.num_nrna
    if n_nrna == 0:
        return

    t_to_nrna = estimator._t_to_nrna
    n_t = estimator.num_transcripts

    # === Phase 1: nRNA-level data aggregation ===

    nrna_exon_sense = np.zeros(n_nrna, dtype=np.float64)
    nrna_exon_anti = np.zeros(n_nrna, dtype=np.float64)
    if t_to_nrna is not None:
        np.add.at(nrna_exon_sense, t_to_nrna, estimator.transcript_exonic_sense)
        np.add.at(nrna_exon_anti, t_to_nrna, estimator.transcript_exonic_antisense)
    else:
        nrna_exon_sense[:] = estimator.transcript_exonic_sense
        nrna_exon_anti[:] = estimator.transcript_exonic_antisense

    nrna_intron_sense = estimator.transcript_intronic_sense.copy()
    nrna_intron_anti = estimator.transcript_intronic_antisense.copy()

    if estimator._exonic_lengths is not None:
        exonic_lengths = estimator._exonic_lengths
    else:
        exonic_lengths = np.ones(n_t, dtype=np.float64)

    t_exon_total = (
        estimator.transcript_exonic_sense
        + estimator.transcript_exonic_antisense
    )

    nrna_wt_len = np.zeros(n_nrna, dtype=np.float64)
    nrna_wt_sum = np.zeros(n_nrna, dtype=np.float64)
    if t_to_nrna is not None:
        np.add.at(nrna_wt_len, t_to_nrna, t_exon_total * exonic_lengths)
        np.add.at(nrna_wt_sum, t_to_nrna, t_exon_total)
    else:
        nrna_wt_len[:] = t_exon_total * exonic_lengths
        nrna_wt_sum[:] = t_exon_total

    nrna_len_fallback = np.zeros(n_nrna, dtype=np.float64)
    nrna_t_count = np.zeros(n_nrna, dtype=np.float64)
    if t_to_nrna is not None:
        np.add.at(nrna_len_fallback, t_to_nrna, exonic_lengths)
        np.add.at(nrna_t_count, t_to_nrna, 1.0)
    else:
        nrna_len_fallback[:] = exonic_lengths
        nrna_t_count[:] = 1.0
    nrna_len_fallback = np.where(
        nrna_t_count > 0, nrna_len_fallback / nrna_t_count, 1.0,
    )

    nrna_exonic_length = np.where(
        nrna_wt_sum > 0, nrna_wt_len / nrna_wt_sum, nrna_len_fallback,
    )
    nrna_intronic_length = np.maximum(1.0, nrna_spans - nrna_exonic_length)

    nrna_locus_ids = np.full(n_nrna, -1, dtype=np.int32)
    if t_to_nrna is not None:
        np.maximum.at(
            nrna_locus_ids, t_to_nrna, estimator.locus_id_per_transcript,
        )
    else:
        nrna_locus_ids[:] = estimator.locus_id_per_transcript

    strand_spec = strand_specificity

    # === Phase 2: 3-tier shrinkage pipeline ===

    nrna_eta, nrna_evidence = _compute_hybrid_nrna_frac_vec(
        nrna_exon_sense, nrna_exon_anti,
        nrna_intron_sense, nrna_intron_anti,
        nrna_exonic_length, nrna_intronic_length,
        strand_spec, gdna_density,
    )

    # Locus-strand η (aggregate nRNA-level by locus × strand)
    has_locus = nrna_locus_ids >= 0
    max_locus = (
        int(nrna_locus_ids.max()) + 1 if has_locus.any() else 0
    )
    ls_key = np.full(n_nrna, -1, dtype=np.int64)
    if max_locus > 0:
        ls_key[has_locus] = (
            nrna_locus_ids[has_locus].astype(np.int64) * 4
            + nrna_strands[has_locus].astype(np.int64)
        )
    n_ls = max(max_locus * 4, 1)
    ls_eta_flat, ls_ev_flat, ls_eta_per_nrna, ls_ev_per_nrna = (
        _aggregate_nrna_frac_by_group(
            ls_key, n_ls,
            nrna_exon_sense, nrna_exon_anti,
            nrna_intron_sense, nrna_intron_anti,
            nrna_exonic_length, nrna_intronic_length,
            strand_spec, gdna_density,
        )
    )

    # === Estimate κ hyperparameters via Method of Moments ===
    _mom_kw = dict(
        kappa_min=kappa_min, kappa_max=kappa_max,
        kappa_fallback=kappa_fallback, kappa_min_obs=kappa_min_obs,
    )
    if kappa_global is None:
        kappa_global = estimate_kappa(
            ls_eta_flat, ls_ev_flat,
            min_evidence=mom_min_evidence_global, **_mom_kw,
        )
    if kappa_locus is None:
        kappa_locus = estimate_kappa(
            nrna_eta, nrna_evidence,
            min_evidence=mom_min_evidence_locus, **_mom_kw,
        )

    # === Empirical global nrna_frac ===
    total_ls_evidence = float(ls_ev_flat.sum())
    if total_ls_evidence > 0:
        nrna_frac_global = float(
            np.dot(ls_ev_flat, ls_eta_flat) / total_ls_evidence,
        )
    else:
        nrna_frac_global = 0.0

    logger.debug(
        f"nrna_frac EB hyperparameters: nrna_frac_global={nrna_frac_global:.4f}, "
        f"κ_global={kappa_global:.1f}, κ_locus={kappa_locus:.1f}, "
        f"κ_nrna={kappa_nrna:.1f}"
    )

    # === Smooth Empirical-Bayes hierarchical shrinkage ===

    nrna_frac_ls = (
        (ls_ev_per_nrna * ls_eta_per_nrna + kappa_global * nrna_frac_global)
        / (ls_ev_per_nrna + kappa_global)
    )

    nrna_frac_est = (
        (nrna_evidence * nrna_eta + kappa_locus * nrna_frac_ls)
        / (nrna_evidence + kappa_locus)
    )

    nrna_frac_est = np.clip(nrna_frac_est, 1e-4, 1.0 - 1e-4)

    estimator.nrna_frac_alpha = nrna_frac_est * kappa_nrna
    estimator.nrna_frac_beta = (1.0 - nrna_frac_est) * kappa_nrna

    # --- Diagnostic logging ---
    median_eta = float(np.median(nrna_frac_est))
    mean_nrna_weight = float(
        np.mean(nrna_evidence / (nrna_evidence + kappa_locus)),
    )
    n_zero = int(np.sum(nrna_frac_est <= 1e-4))
    n_high = int(np.sum(nrna_frac_est >= 0.1))
    logger.info(
        f"nrna_frac priors: global={nrna_frac_global:.4f}, "
        f"median={median_eta:.4f}, "
        f"κ=[{kappa_global:.1f}, {kappa_locus:.1f}, {kappa_nrna:.1f}], "
        f"nRNA_wt={mean_nrna_weight:.2f}, "
        f"n_zero={n_zero}/{n_nrna}, n_high(≥0.1)={n_high}/{n_nrna}"
    )
