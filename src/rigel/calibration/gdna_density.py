"""The gDNA background rate from per-object counts and their opportunity; layer 2.

:func:`one_sided_rate` estimates the rate on a pool that transcription may have added to, the one the
fragment-length model uses, and reports the naive pooled rate ``sum n / sum E`` beside it. No structural class of
objects is pure gDNA, but transcription can only add fragments to an object, so every observed density ``n_i / E_i``
overestimates the gDNA rate and the clean objects sit on the low side. The estimator assumes the direction of
contamination, which the biology gives, instead of a pure pool, which it does not.

With ``lam_i = rho * E_i`` and pure gDNA ``n_i ~ Poisson(lam_i)``, the negative part of the residual has a closed
form, De Moivre's mean-absolute-deviation identity::

    E[(lam - N)+]  =  1/2 E|N - lam|  =  lam^(floor(lam)+1) * exp(-lam) / floor(lam)!   =:  D(lam)

Contamination only raises ``n_i`` and so only lowers ``(lam_i - n_i)+``. Matching the observed one-sided mass to its
null expectation is one scalar equation, with no trim depth or threshold::

    F(rho)  =  sum_i (rho*E_i - n_i)+  -  sum_i D(rho*E_i)  =  0

Its one bias is one-signed: a heavily contaminated object adds 0 to the first sum but its full ``D(lam_i)`` to the
second, so the root is pushed up, never down.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

from ..config import CONSTANTS

__all__ = [
    "GdnaDensityFit",
    "one_sided_rate",
    "poisson_lower_mean",
]


@dataclass(frozen=True, slots=True)
class GdnaDensityFit:
    """A gDNA density estimate with its diagnostics.

    ``rate`` is the one-sided estimate and ``pooled_rate`` the naive one. ``rate_over_pooled`` says how much
    contamination the fit found; near 1 means the pool was already clean.
    """

    rate: float
    pooled_rate: float
    n_objects: int
    total_exposure: float
    bracket_ok: bool

    @property
    def rate_over_pooled(self) -> float:
        """``rate / pooled_rate`` — 1.0 when the pool was clean, below 1 when contamination was found."""
        return self.rate / self.pooled_rate if self.pooled_rate > 0.0 else float("nan")

    @property
    def informative(self) -> bool:
        """False when there was no support to fit on, or the bracket failed — the caller must fall back."""
        return bool(self.bracket_ok and self.n_objects > 0 and self.total_exposure > 0.0)


def poisson_lower_mean(lam: np.ndarray) -> np.ndarray:
    """``E[(lam - N)+]`` for ``N ~ Poisson(lam)``, exactly: ``lam^(k+1) e^-lam / k!`` with ``k = floor(lam)``.

    De Moivre's mean-absolute-deviation identity for the Poisson, halved, since ``E[N - lam] = 0`` makes the two
    one-sided means equal.
    """
    lam = np.asarray(lam, dtype=np.float64)
    out = np.zeros(lam.shape, dtype=np.float64)
    pos = lam > 0.0
    if not np.any(pos):
        return out
    lp = lam[pos]
    k = np.floor(lp)
    # in log space: a large `lam` overflows `lam**(k+1)` and `k!` separately while their ratio is O(sqrt(lam))
    out[pos] = np.exp((k + 1.0) * np.log(lp) - lp - _log_factorial(k))
    return out


def _log_factorial(k: np.ndarray) -> np.ndarray:
    """``gammaln(k + 1)`` for integer-valued ``k``, by table whenever the table is smaller than the data.

    The arguments are the integers up to the largest, usually far fewer than the objects, and the table holds the
    same function at the same points, so its values are identical, not approximated.
    """
    if k.size == 0:
        return gammaln(k + 1.0)
    kmax = float(k.max())
    if not np.isfinite(kmax) or kmax + 2.0 > k.size:
        return gammaln(k + 1.0)
    table = gammaln(np.arange(int(kmax) + 2, dtype=np.float64))
    return table[k.astype(np.int64) + 1]


def one_sided_rate(counts, exposure) -> GdnaDensityFit:
    """The gDNA density from the low side of the per-object density distribution.

    Solves ``F(rho) = sum (rho*E_i - n_i)+ - sum D(rho*E_i) = 0`` by bisection. ``F(0) = 0`` exactly, and just above
    zero ``F`` is negative: ``D(lam) -> lam`` as ``lam -> 0``, so ``F(0+) ~ -rho * sum_{n_i>0} E_i``. ``F`` then rises
    without bound, since the first sum is eventually linear in ``rho`` and the second grows as ``sqrt(rho)``, so there
    is exactly one root above zero. ``F`` is not monotone, but bisection needs only ``F(lo) <= 0 <= F(hi)``, and
    because ``F(0+) < 0`` whenever any object has a count (the early return ensures one does), the lower end walks
    off zero on its own.

    The loop ends when the bracket closes to adjacent floats: their midpoint rounds to one of them, so ``½(lo + hi)``
    cannot move again, and no iteration budget is needed.
    """
    n = np.asarray(counts, dtype=np.float64).ravel()
    e = np.asarray(exposure, dtype=np.float64).ravel()
    if n.shape != e.shape:
        raise ValueError(f"counts {n.shape} and exposure {e.shape} must be aligned per object")
    live = e > 0.0
    n, e = n[live], e[live]
    total_e = float(e.sum())
    total_n = float(n.sum())
    pooled = total_n / total_e if total_e > 0.0 else 0.0
    if not (total_e > 0.0 and total_n > 0.0):
        # no support or no fragments: there is no density to estimate, and a zero-gDNA library must not be
        # handed a fabricated rate
        return GdnaDensityFit(0.0, pooled, int(n.size), total_e, bracket_ok=False)

    def f(rho: float) -> float:
        lam = rho * e
        return float(np.clip(lam - n, 0.0, None).sum() - poisson_lower_mean(lam).sum())

    lo, hi = 0.0, CONSTANTS.calibration.bracket_headroom * pooled
    if not (f(lo) <= 0.0 <= f(hi)):
        return GdnaDensityFit(pooled, pooled, int(n.size), total_e, bracket_ok=False)
    while True:
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:  # the bracket is closed to the last float
            break
        if f(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return GdnaDensityFit(0.5 * (lo + hi), pooled, int(n.size), total_e, bracket_ok=True)
