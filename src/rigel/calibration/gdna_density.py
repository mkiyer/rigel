"""The gDNA background RATE, from counts and their opportunity.

:func:`one_sided_rate` estimates it for a pool that transcription may have added to, and it is the one the
fragment-length model uses. Beside its estimate it reports the naive pooled rate ``sum n / sum E``, which
assumes the pool it is handed is pure gDNA.

The one-sided argument, which is why the estimator exists: no structural class of objects is
pure gDNA. "Intergenic" is whatever the annotation leaves over, RNA sits inside introns wherever it
has not spliced, and pervasive transcription is real. But transcription can only ADD fragments to an
object, never remove them, so every object's observed density ``n_i / E_i`` is an overestimate and
the clean objects are the ones on the LOW side. Using the side a contaminant cannot reach is an
assumption about DIRECTION, which the biology gives, in place of an assumption about PURITY, which
it does not.

The estimator. With ``lam_i = rho * E_i``, under pure gDNA ``n_i ~ Poisson(lam_i)`` and ``E[n_i - lam_i]
= 0``, so the negative part of the residual has a closed form — De Moivre's mean-absolute-deviation
identity for the Poisson::

    E[(lam - N)+]  =  1/2 E|N - lam|  =  lam^(floor(lam)+1) * exp(-lam) / floor(lam)!   =:  D(lam)

Contamination only RAISES ``n_i``, hence only LOWERS ``(lam_i - n_i)+``. Matching the observed one-sided
mass to its exact null expectation is then one scalar equation::

    F(rho)  =  sum_i (rho*E_i - n_i)+  -  sum_i D(rho*E_i)  =  0

No trim depth, no quantile, no threshold, no chosen level. Keeping the lowest ``q`` of objects by
density and pooling those also works, but a safe ``q`` is a function of the library's contaminated
fraction — a property of the sample, not of the tool. The identity above is that estimator's limit
taken exactly.

The one bias is one-signed: a heavily contaminated object contributes ``0`` to the first sum but
its full ``D(lam_i)`` to the second, so the root is pushed UP. Against an origin-split oracle the
estimate runs a few percent high and is never low. How that bias scales with the contaminated
fraction is not established, so a library dirtier than any simulated panel may sit outside the
range the panels show.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

from ..config import CONSTANTS

__all__ = [
    "GdnaDensityFit",
    "contained_opportunity",
    "one_sided_rate",
    "poisson_lower_mean",
]


@dataclass(frozen=True, slots=True)
class GdnaDensityFit:
    """A gDNA density estimate with the diagnostics that make it auditable.

    ``rate`` is the one-sided estimate; ``pooled_rate`` is the naive one the pool would have given.
    ``rate_over_pooled`` is the number to read: it is how much contamination the fit found, and a
    value near 1 means the pool was already clean rather than that the estimator did nothing.
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

    De Moivre's mean-absolute-deviation identity for the Poisson, halved, since ``E[N - lam] = 0``
    makes the two one-sided means equal. Closed form at every ``lam``, so the estimator that uses it
    needs no simulation and no tabulation; gated against the exact truncated sum and Monte Carlo.
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

    The log-gamma above is evaluated at ``floor(lam) + 1``, which is an INTEGER, and `one_sided_rate`
    evaluates it over every object at every step of its bisection. The distinct arguments are therefore the
    integers up to the largest, which on a real library is orders of magnitude fewer than the objects. The
    table holds that same function at those same points, so its values are IDENTICAL rather than
    approximated, and the direct form stays for the case where the table would be the larger of the two.
    """
    if k.size == 0:
        return gammaln(k + 1.0)
    kmax = float(k.max())
    if not np.isfinite(kmax) or kmax + 2.0 > k.size:
        return gammaln(k + 1.0)
    table = gammaln(np.arange(int(kmax) + 2, dtype=np.float64))
    return table[k.astype(np.int64) + 1]


def one_sided_rate(counts, exposure) -> GdnaDensityFit:
    """The gDNA density from the LOW side of the per-object density distribution.

    Solves ``F(rho) = sum (rho*E_i - n_i)+ - sum D(rho*E_i) = 0`` by bisection.

    The bracket is structural, not guessed. ``F(0) = 0`` exactly, since both sums vanish, and just
    above zero ``F`` is strictly negative: ``D(lam) -> lam`` as ``lam -> 0``, so the second sum
    approaches ``rho * sum E`` while the first collects only the objects with ``n_i = 0``, leaving
    ``F(0+) ~ -rho * sum_{n_i>0} E_i``. ``F`` then rises without bound, because the first sum is
    eventually linear in ``rho`` while the second grows only as ``sqrt(rho)``. So there is exactly
    one root above zero.

    ``F`` is NOT monotone — it dips before it rises — but bisection does not need it to be, only
    ``F(lo) <= 0 <= F(hi)``. The root at ``rho = 0`` is not sticky either: because ``F(0+) < 0``
    whenever any object carries a count, the first midpoints test negative and the lower end walks
    off zero on its own. The early return above is what guarantees that premise, so no separate
    guard against the boundary root is needed.

    THE LOOP ENDS WHEN THE BRACKET DOES, and that is a property of float64 rather than a step count: the
    midpoint of two ADJACENT floats rounds to one of them, so once ``mid`` is no longer strictly inside
    ``(lo, hi)`` every further halving reassigns the value it already holds and the returned
    ``½(lo + hi)`` cannot move again. Testing for that is the whole termination rule — no iteration
    budget to choose and none to justify.
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
        # No support, or no fragments at all: there is no density to estimate and saying so is the
        # answer. A zero-gDNA library reaches this and must NOT be handed a fabricated rate.
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


def contained_opportunity(pmf, lengths) -> np.ndarray:
    """``E_f[(ell - L + 1)+]`` per object — the CONTAINED opportunity a uniformly placed fragment has.

    Evaluated by cumulative sums rather than an ``(objects x lengths)`` outer product, which at
    genome scale would be tens of billions of entries. Writing ``F`` and ``S`` for the pmf's
    cumulative mass and cumulative first moment, ``E[(ell-L+1)+] = (ell+1)*F(ell) - S(ell)``, and
    for ``ell`` beyond the pmf's support that is exactly ``ell + 1 - E[L]``. Cost is
    ``O(objects + max_length)``.
    """
    p = np.asarray(pmf, dtype=np.float64)
    ell = np.asarray(lengths, dtype=np.float64)
    if p.size == 0:
        return np.zeros(ell.shape, dtype=np.float64)
    L = np.arange(p.size, dtype=np.float64)
    cum_mass = np.cumsum(p)
    cum_first = np.cumsum(p * L)
    idx = np.clip(ell, 0, p.size - 1).astype(np.int64)
    out = (ell + 1.0) * cum_mass[idx] - cum_first[idx]
    return np.maximum(out, 0.0)
