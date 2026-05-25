"""Symmetric beta-binomial gDNA strand-balance estimator.

Moved from ``calibration.density_global`` in v4 Phase 0
(``docs/fineregions/density_model_impl_plan_v4.md`` §4).
The strand-balance estimator is a strand-model concern that has nothing
to do with density modelling; keeping it in ``density_global`` was a
historical accident from the v3 monolith.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


__all__ = [
    "STRAND_KAPPA_DEFAULT",
    "STRAND_KAPPA_MIN",
    "STRAND_KAPPA_MAX",
    "MIN_REGIONS_FOR_STRAND_MOM",
    "StrandBalanceEstimate",
    "estimate_strand_balance",
]


STRAND_KAPPA_DEFAULT: float = 1.0e6
STRAND_KAPPA_MIN: float = 1.0e-3
STRAND_KAPPA_MAX: float = 1.0e6
MIN_REGIONS_FOR_STRAND_MOM: int = 2


@dataclass(frozen=True, slots=True)
class StrandBalanceEstimate:
    """Symmetric beta-binomial strand overdispersion estimate for gDNA."""

    kappa: float
    n_regions: int
    n_fragments: float
    n_pos: float
    n_neg: float
    observed_pos_fraction: float
    residual_sum: float
    binomial_variance_sum: float
    max_overdispersed_variance_sum: float
    overdispersion_factor: float
    fallback_used: bool = False
    fallback_reason: str = ""

    @property
    def alpha(self) -> float:
        return float(self.kappa / 2.0)

    @property
    def beta(self) -> float:
        return float(self.kappa / 2.0)

    def to_summary_dict(self) -> dict[str, object]:
        return {
            "model": "symmetric_beta_binomial",
            "mean_fixed": 0.5,
            "kappa": float(self.kappa),
            "alpha": self.alpha,
            "beta": self.beta,
            "n_regions": int(self.n_regions),
            "n_fragments": float(self.n_fragments),
            "n_pos": float(self.n_pos),
            "n_neg": float(self.n_neg),
            "observed_pos_fraction": float(self.observed_pos_fraction),
            "residual_sum": float(self.residual_sum),
            "binomial_variance_sum": float(self.binomial_variance_sum),
            "max_overdispersed_variance_sum": float(self.max_overdispersed_variance_sum),
            "overdispersion_factor": float(self.overdispersion_factor),
            "fallback_used": bool(self.fallback_used),
            "fallback_reason": str(self.fallback_reason),
        }


def estimate_strand_balance(
    pos_counts: np.ndarray,
    neg_counts: np.ndarray,
    region_mask: np.ndarray,
) -> StrandBalanceEstimate:
    """Estimate beta-binomial strand overdispersion with fixed mean 0.5."""
    pos = np.asarray(pos_counts, dtype=np.float64)
    neg = np.asarray(neg_counts, dtype=np.float64)
    if pos.shape != neg.shape:
        raise ValueError(
            f"estimate_strand_balance: pos.shape ({pos.shape}) != neg.shape ({neg.shape})."
        )

    mask = np.asarray(region_mask, dtype=bool)
    if mask.shape != pos.shape:
        raise ValueError(
            f"estimate_strand_balance: region_mask.shape ({mask.shape}) != pos.shape ({pos.shape})."
        )

    total = pos + neg
    eligible = mask & (total > 0.0)
    n_regions = int(eligible.sum())
    if n_regions < MIN_REGIONS_FOR_STRAND_MOM:
        return _strand_fallback(
            pos[eligible],
            neg[eligible],
            n_regions,
            "n_regions < MIN_REGIONS_FOR_STRAND_MOM",
        )

    pos_e = pos[eligible]
    neg_e = neg[eligible]
    total_e = total[eligible]
    residual = pos_e - 0.5 * total_e
    residual_sum = float(np.sum(residual * residual))
    binomial_variance_sum = float(0.25 * np.sum(total_e))
    max_variance_sum = float(0.25 * np.sum(total_e * total_e))
    n_fragments = float(total_e.sum())
    n_pos = float(pos_e.sum())
    n_neg = float(neg_e.sum())
    observed_pos_fraction = float(n_pos / n_fragments) if n_fragments > 0.0 else 0.5
    overdispersion_factor = (
        float(residual_sum / binomial_variance_sum) if binomial_variance_sum > 0.0 else 0.0
    )

    if binomial_variance_sum <= 0.0:
        return _strand_fallback(pos_e, neg_e, n_regions, "no positive strand exposure")
    if residual_sum <= binomial_variance_sum:
        return StrandBalanceEstimate(
            kappa=STRAND_KAPPA_MAX,
            n_regions=n_regions,
            n_fragments=n_fragments,
            n_pos=n_pos,
            n_neg=n_neg,
            observed_pos_fraction=observed_pos_fraction,
            residual_sum=residual_sum,
            binomial_variance_sum=binomial_variance_sum,
            max_overdispersed_variance_sum=max_variance_sum,
            overdispersion_factor=overdispersion_factor,
            fallback_used=True,
            fallback_reason="residual variance <= binomial expectation",
        )

    numerator = max_variance_sum - residual_sum
    denominator = residual_sum - binomial_variance_sum
    if denominator <= 0.0 or not np.isfinite(denominator):
        kappa = STRAND_KAPPA_MAX
        fallback_used = True
        fallback_reason = "invalid MoM denominator"
    else:
        kappa = numerator / denominator
        fallback_used = False
        fallback_reason = ""
    if not np.isfinite(kappa):
        kappa = STRAND_KAPPA_MAX
        fallback_used = True
        fallback_reason = "non-finite MoM estimate"
    kappa = float(np.clip(kappa, STRAND_KAPPA_MIN, STRAND_KAPPA_MAX))

    return StrandBalanceEstimate(
        kappa=kappa,
        n_regions=n_regions,
        n_fragments=n_fragments,
        n_pos=n_pos,
        n_neg=n_neg,
        observed_pos_fraction=observed_pos_fraction,
        residual_sum=residual_sum,
        binomial_variance_sum=binomial_variance_sum,
        max_overdispersed_variance_sum=max_variance_sum,
        overdispersion_factor=overdispersion_factor,
        fallback_used=fallback_used,
        fallback_reason=fallback_reason,
    )


def _strand_fallback(
    pos: np.ndarray,
    neg: np.ndarray,
    n_regions: int,
    reason: str,
) -> StrandBalanceEstimate:
    n_pos = float(np.asarray(pos, dtype=np.float64).sum())
    n_neg = float(np.asarray(neg, dtype=np.float64).sum())
    n_fragments = n_pos + n_neg
    return StrandBalanceEstimate(
        kappa=STRAND_KAPPA_DEFAULT,
        n_regions=int(n_regions),
        n_fragments=n_fragments,
        n_pos=n_pos,
        n_neg=n_neg,
        observed_pos_fraction=float(n_pos / n_fragments) if n_fragments > 0.0 else 0.5,
        residual_sum=0.0,
        binomial_variance_sum=0.0,
        max_overdispersed_variance_sum=0.0,
        overdispersion_factor=0.0,
        fallback_used=True,
        fallback_reason=reason,
    )
