"""rigel.priors — Empirical Bayes prior computation for the EM solver.

All prior functions that run between the scan pass and the locus EM
loop live here:

- ``compute_global_gdna_density`` — global gDNA density via strand
  correction.
- ``estimate_kappa`` — Method-of-Moments Beta concentration κ.

These are pure functions operating on numpy arrays and the
``AbundanceEstimator`` accumulator state.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .estimator import AbundanceEstimator

logger = logging.getLogger(__name__)


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
            total_sense,
            total_anti,
            total_exonic_bp,
            strand_specificity,
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
    kappa_min: float = 2.0,
    kappa_max: float = 200.0,
    kappa_fallback: float = 5.0,
    kappa_min_obs: int = 20,
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
