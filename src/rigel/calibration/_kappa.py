"""rigel.calibration._kappa — Per-region NB Method-of-Moments overdispersion.

A single global negative-binomial $\\kappa$ estimate per density type.
Used by the M6 locoregional EB shrinkage.  See
``docs/calibration/calibration_v6_plan.md`` §2.6.3 and
``docs/calibration/m4_implementation_plan.md`` §4.

The estimator is the per-region MoM solution under the model
$N_R \\sim \\mathrm{NB}(\\mu_R, \\kappa)$ with $\\mu_R = \\hat\\rho \\cdot
L_{\\mathrm{eff}}(R)$:

$$\\hat\\kappa = \\frac{\\sum_R \\mu_R^2}{\\sum_R (N_R - \\mu_R)^2 - \\sum_R \\mu_R}$$

Per-region (rather than pooled) MoM is required because exposures
$L_{\\mathrm{eff}}(R)$ span 3-4 orders of magnitude (short exons vs whole
introns); a single pooled mean conflates true overdispersion with
exposure heterogeneity.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


KAPPA_DEFAULT = 100.0
KAPPA_MIN = 1.0
KAPPA_MAX = 1.0e6
MIN_REGIONS_FOR_MOM = 5


@dataclass(frozen=True, slots=True)
class KappaEstimate:
    """One NB overdispersion estimate plus a fallback diagnostic."""

    value: float          # always finite, always in [KAPPA_MIN, KAPPA_MAX]
    n_regions: int        # rows that entered the MoM sum
    fallback_used: bool   # True ⇒ value == KAPPA_DEFAULT
    fallback_reason: str  # "" iff fallback_used is False


def estimate_kappa(
    counts: np.ndarray,
    eff_lengths: np.ndarray,
    rho_hat: float,
) -> KappaEstimate:
    """Per-region NB Method-of-Moments overdispersion estimate.

    Parameters
    ----------
    counts
        ``(R,)`` int64 per-region observed counts (the same regions
        whose densities went into ``rho_hat``).
    eff_lengths
        ``(R,)`` float64 per-region $L_{\\mathrm{eff}}$.
    rho_hat
        Global density (fragments / bp) for this region type.

    Falls back to ``KAPPA_DEFAULT`` (recorded with a reason) when:

    * fewer than ``MIN_REGIONS_FOR_MOM`` regions, or
    * ``rho_hat <= 0``, or all ``eff_lengths == 0``, or
    * the MoM excess variance is non-positive (Poisson or
      under-dispersed input).

    The returned value is always clipped to ``[KAPPA_MIN, KAPPA_MAX]``.
    """
    if counts.shape != eff_lengths.shape:
        raise ValueError(
            f"estimate_kappa: counts.shape ({counts.shape}) != "
            f"eff_lengths.shape ({eff_lengths.shape})."
        )
    n = int(counts.size)
    if n < MIN_REGIONS_FOR_MOM:
        return _fallback(n, "n_regions < MIN_REGIONS_FOR_MOM")
    if rho_hat <= 0.0:
        return _fallback(n, "rho_hat <= 0")
    if not bool(np.any(eff_lengths > 0)):
        return _fallback(n, "all eff_lengths == 0")

    mu = rho_hat * eff_lengths.astype(np.float64, copy=False)
    diff = counts.astype(np.float64, copy=False) - mu
    excess = float(np.sum(diff * diff) - np.sum(mu))
    if excess <= 0.0:
        return _fallback(n, "excess variance <= 0 (Poisson or under-dispersed)")

    kappa = float(np.sum(mu * mu) / excess)
    kappa = float(np.clip(kappa, KAPPA_MIN, KAPPA_MAX))
    return KappaEstimate(
        value=kappa, n_regions=n, fallback_used=False, fallback_reason=""
    )


def _fallback(n_regions: int, reason: str) -> KappaEstimate:
    return KappaEstimate(
        value=KAPPA_DEFAULT,
        n_regions=n_regions,
        fallback_used=True,
        fallback_reason=reason,
    )
