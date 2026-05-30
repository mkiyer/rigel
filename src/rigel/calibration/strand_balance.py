"""RNA strand-balance model (D2): κ_rna (mean) + ρ_r_bb (overdispersion).

The E-step's strand log-Bayes-factor (PR 4) compares an observed sense count
against an RNA Beta-Binomial ``BB(n, κ_rna, ρ_r_bb)`` and a gDNA
``BB(n, 0.5, ρ_d_bb)``. This module supplies the RNA parameters:

* **``κ_rna``** (the RNA sense mean) is **not fit here** — it is the live
  ``StrandModel``'s ``p_r1_sense`` (the MLE over annotated spliced unique
  mappers, an uncontaminated measure of library strand specificity; doc 01
  §4.3, PR03 §III.1).
* **``ρ_r_bb``** (the strand overdispersion) **is** fit here, by
  method-of-moments over integer ``(k_sense, n)`` spliced junction
  observations pooled from the substrate's three views (contained + per-side
  boundary). The substrate's spliced channels are motif-oriented at deposit
  (PR 2.5), so observations are valid even in strand-ambiguous regions — no
  ``ts_class`` filtering is needed.

Overdispersion is the robustness mechanism for the strand channel (it absorbs
imbalanced regions without an outlier cliff). Whether a library truly needs it
(vs. a Binomial, ``ρ_r_bb ≈ 0``) is to be confirmed on real data (PR 7).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..strand_model import StrandModels
    from .substrate import CalibrationSubstrate

# Numerical floor keeping the Beta-Binomial shape params α = κ(1−ρ)/ρ and
# β = (1−κ)(1−ρ)/ρ finite and positive — both κ_rna and ρ_r_bb are clamped to
# (_BB_FLOOR, 1 − _BB_FLOOR). Spec: doc 03 §8 `_BB_FLOOR`.
_BB_FLOOR = 1.0e-6

# No-data fallback overdispersion: a mild, non-zero default so the strand
# channel is not treated as a sharp (Binomial) discriminator when there is too
# little spliced evidence to estimate dispersion. Spec: doc 03 §7 init value
# (matches the placeholder's _RHO_R_BB_INIT).
_RHO_R_BB_FALLBACK = 0.01

# Minimum distinct spliced observations needed to estimate an overdispersion:
# the cross-observation variance is undefined for fewer than two. Not a
# tunable — a mathematical floor.
_MIN_OBS_FOR_OVERDISPERSION = 2


@dataclass(frozen=True, slots=True)
class StrandBalance:
    """Fitted RNA strand model: mean κ_rna + Beta-Binomial overdispersion ρ_r_bb."""

    kappa_rna: float  # RNA sense mean, in (0, 1)
    rho_r_bb: float  # RNA strand BB overdispersion, in (0, 1)
    n_observations: int  # number of pooled spliced junction observations (n > 0)
    n_total_reads: int  # total spliced reads pooled
    fallback_used: bool
    fallback_reason: str


def _mom_overdispersion(k_sense: np.ndarray, n: np.ndarray, kappa: float) -> float:
    """Method-of-moments overdispersion ρ for a Beta-Binomial of known mean κ.

    With ``Var(k_i) = n_i κ(1−κ)[1 + (n_i − 1)ρ]``, equate the pooled observed
    squared residuals to the pooled expected variance and solve for ρ::

        ρ = [ Σ(k_i − κ n_i)² / (κ(1−κ)) − Σ n_i ] / Σ n_i(n_i − 1)

    Returns 0.0 for the degenerate cases (κ at a bound, or every n_i == 1 so
    overdispersion is unseparable from the binomial term); the caller clamps.
    """
    var_unit = kappa * (1.0 - kappa)
    sum_n_nm1 = float(np.sum(n * (n - 1.0)))
    if var_unit <= 0.0 or sum_n_nm1 <= 0.0:
        return 0.0
    resid_ss = float(np.sum((k_sense - kappa * n) ** 2))
    return (resid_ss / var_unit - float(np.sum(n))) / sum_n_nm1


def fit_strand_balance(
    substrate: "CalibrationSubstrate",
    strand_model: "StrandModels",
) -> StrandBalance:
    """Fit the RNA strand-balance model from the substrate's spliced channels."""
    kappa = min(max(float(strand_model.p_r1_sense), _BB_FLOOR), 1.0 - _BB_FLOOR)

    # Pool integer (k_sense, n) spliced observations from all three views.
    k_parts: list[np.ndarray] = []
    n_parts: list[np.ndarray] = []
    for view in (substrate.contained, substrate.left, substrate.right):
        n = view.n_spliced
        mask = n > 0
        if np.any(mask):
            k_parts.append(view.n_spliced_sense[mask].astype(np.float64))
            n_parts.append(n[mask].astype(np.float64))

    if not n_parts:
        return StrandBalance(
            kappa_rna=kappa,
            rho_r_bb=_RHO_R_BB_FALLBACK,
            n_observations=0,
            n_total_reads=0,
            fallback_used=True,
            fallback_reason="no spliced observations",
        )

    k_sense = np.concatenate(k_parts)
    n = np.concatenate(n_parts)
    n_obs = int(k_sense.shape[0])
    n_total = int(n.sum())

    if n_obs < _MIN_OBS_FOR_OVERDISPERSION:
        return StrandBalance(
            kappa_rna=kappa,
            rho_r_bb=_RHO_R_BB_FALLBACK,
            n_observations=n_obs,
            n_total_reads=n_total,
            fallback_used=True,
            fallback_reason=(
                f"only {n_obs} spliced observation(s) (< {_MIN_OBS_FOR_OVERDISPERSION})"
            ),
        )

    rho = _mom_overdispersion(k_sense, n, kappa)
    rho = min(max(rho, _BB_FLOOR), 1.0 - _BB_FLOOR)
    return StrandBalance(
        kappa_rna=kappa,
        rho_r_bb=rho,
        n_observations=n_obs,
        n_total_reads=n_total,
        fallback_used=False,
        fallback_reason="",
    )


__all__ = ["StrandBalance", "fit_strand_balance"]
