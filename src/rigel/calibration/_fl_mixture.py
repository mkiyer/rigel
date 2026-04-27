"""SRD v1 Pass 2 — 1-D mixture EM for the unified gDNA pool.

Fits the well-posed mixture

    pool_FL(l) = π · gDNA_FL(l) + (1 − π) · RNA_FL(l)

where ``RNA_FL(l)`` is fixed (provided by the caller) and ``gDNA_FL(l)``
is a free Dirichlet-smoothed histogram. ``π`` is a scalar in [0, 1].

EM is convex given fixed ``RNA_FL`` and converges in <30 iterations.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FLMixtureResult:
    """Output of the 1-D pool mixture fit.

    Attributes
    ----------
    gdna_counts : np.ndarray[float64]
        Recovered ``gDNA_FL`` histogram, scaled to ``π · N_pool`` so a
        downstream Dirichlet shrinkage sees the correct effective
        sample size.
    pi : float
        Pool-level gDNA fraction.
    converged : bool
        True iff EM stopped because the change in ``π`` fell below
        ``tol``.
    n_iter : int
        Number of EM iterations actually run.
    n_pool : int
        Total fragment count in the input pool histogram.
    """

    gdna_counts: np.ndarray
    pi: float
    converged: bool
    n_iter: int
    n_pool: int


def fit_fl_mixture(
    pool_hist: np.ndarray,
    rna_fl_probs: np.ndarray,
    *,
    max_iter: int = 50,
    tol: float = 1e-4,
    pi_init: float = 0.1,
    smoothing: float = 1.0,
) -> FLMixtureResult:
    """Fit ``pool = π·gDNA + (1−π)·RNA`` with fixed RNA component.

    Parameters
    ----------
    pool_hist : np.ndarray[float64]
        Observed FL counts in the unified pool (length ``L``).
    rna_fl_probs : np.ndarray[float64]
        Fixed RNA FL probability vector (length ``L``); must sum to ~1.
    max_iter, tol
        EM stopping criteria (on ``|Δπ|``).
    pi_init
        Initial mixing fraction.
    smoothing
        Per-bin additive smoothing for the gDNA histogram each step,
        preventing zero-probability bins.
    """
    pool = np.asarray(pool_hist, dtype=np.float64)
    rna = np.asarray(rna_fl_probs, dtype=np.float64)
    n_bins = pool.size
    n_pool = int(pool.sum())

    if rna.size != n_bins:
        raise ValueError(
            f"rna_fl_probs length {rna.size} != pool_hist length {n_bins}"
        )

    if n_pool == 0:
        return FLMixtureResult(
            gdna_counts=np.zeros(n_bins, dtype=np.float64),
            pi=0.0,
            converged=False,
            n_iter=0,
            n_pool=0,
        )

    # Defensive normalize of RNA component.
    rna_sum = float(rna.sum())
    rna = rna / rna_sum if rna_sum > 0 else np.full(n_bins, 1.0 / n_bins)

    # Initialize gDNA component as the smoothed pool marginal.
    gdna = pool + smoothing
    gdna /= gdna.sum()

    pi = float(pi_init)
    converged = False
    n_iter = 0

    for it in range(1, max_iter + 1):
        n_iter = it
        num = pi * gdna
        den = num + (1.0 - pi) * rna
        with np.errstate(divide="ignore", invalid="ignore"):
            r_g = np.where(den > 0, num / den, 0.0)

        e_gdna = pool * r_g
        new_pi = float(e_gdna.sum() / max(pool.sum(), 1e-300))
        new_gdna = e_gdna + smoothing
        s = new_gdna.sum()
        new_gdna = new_gdna / s if s > 0 else np.full(n_bins, 1.0 / n_bins)

        if abs(new_pi - pi) < tol:
            pi = new_pi
            gdna = new_gdna
            converged = True
            break
        pi = new_pi
        gdna = new_gdna

    gdna_counts = gdna * (pi * float(pool.sum()))
    return FLMixtureResult(
        gdna_counts=gdna_counts,
        pi=pi,
        converged=converged,
        n_iter=n_iter,
        n_pool=n_pool,
    )
