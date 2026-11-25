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

Over-dispersed gDNA, each object's rate scattering as a Gamma of shape ``s`` about ``rho``, is
``n_i ~ NegBinom(mean lam_i, size s)``, and its lower mean deviation has the same kind of closed form::

    E[(lam - X)+]  =  (lam/s) * (floor(lam) + s) * P(X = floor(lam))

which telescopes as De Moivre's does and is his at ``s = inf``. :func:`one_sided_rate` takes ``s`` (Poisson by default,
as the fragment-length model reads its pool), and :func:`background_law` fits a pool's whole law, rate and spread, for
the intron abundance model.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import betaln, gammaln

from ..config import CONSTANTS
from .numeric import bisect_to_adjacent

__all__ = [
    "BackgroundLaw",
    "GdnaDensityFit",
    "background_law",
    "log_negbinom",
    "nb_lower_mean",
    "one_sided_rate",
    "poisson_lower_mean",
]


@dataclass(frozen=True, slots=True)
class GdnaDensityFit:
    """A gDNA density estimate with its diagnostics.

    ``rate`` is the one-sided estimate and ``pooled_rate`` the naive one; their ratio says how much contamination the
    fit found, near 1 when the pool was already clean.
    """

    rate: float
    pooled_rate: float
    n_objects: int
    total_exposure: float
    bracket_ok: bool

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
    out[pos] = np.exp((k + 1.0) * np.log(lp) - lp - log_factorial(k))
    return out


def log_negbinom(n, lam, size: float) -> np.ndarray:
    """``log P(X = n)`` for ``X ~ NegBinom(mean lam, size)`` at real ``n >= 0``; ``size = inf`` is the Poisson.

    Stable at any size: ``log Γ(n + size) − log Γ(size)`` is taken as ``log Γ(n) − log B(n, size)`` (0 at ``n = 0``),
    since the direct difference of two log-gammas loses every digit once ``size`` dwarfs ``n``, and ``log(size/(size +
    lam))`` as ``−log1p(lam/size)``."""
    n = np.asarray(n, dtype=np.float64)
    lam = np.maximum(np.asarray(lam, dtype=np.float64), CONSTANTS.calibration.rate_log_floor)
    if not np.isfinite(size):
        return n * np.log(lam) - lam - gammaln(n + 1.0)
    safe = np.where(n > 0.0, n, 1.0)
    rising = np.where(n > 0.0, gammaln(safe) - betaln(safe, size), 0.0)
    shrink = np.log1p(lam / size)
    return rising - gammaln(n + 1.0) - size * shrink + n * (np.log(lam) - np.log(size) - shrink)


def nb_lower_mean(lam: np.ndarray, size: float) -> np.ndarray:
    """``E[(lam - X)+]`` for ``X ~ NegBinom(mean lam, size)``, exactly: ``(lam/size)(m + size) P(X = m)``, ``m =
    floor(lam)``.

    The truncated sum ``sum_{k<=m} (lam - k) p_k`` telescopes: ``(lam - k) p_k = c_k - c_(k-1)`` with ``c_k =
    (lam/size)(k + size) p_k``, from the negative binomial's ratio ``k p_k = (k - 1 + size) (lam/(lam + size))
    p_(k-1)``. At ``size = inf`` it is De Moivre's identity, :func:`poisson_lower_mean`.
    """
    if not np.isfinite(size):
        return poisson_lower_mean(lam)
    lam = np.asarray(lam, dtype=np.float64)
    out = np.zeros(lam.shape, dtype=np.float64)
    pos = lam > 0.0
    if not np.any(pos):
        return out
    lp = lam[pos]
    m = np.floor(lp)
    out[pos] = lp * (1.0 + m / size) * np.exp(log_negbinom(m, lp, size))
    return out


def log_factorial(k: np.ndarray) -> np.ndarray:
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


def one_sided_rate(counts, exposure, size: float = np.inf) -> GdnaDensityFit:
    """The gDNA density from the low side of the per-object density distribution.

    Solves ``F(rho) = sum (rho*E_i - n_i)+ - sum D(rho*E_i) = 0`` by bisection, ``D`` the lower mean deviation of a
    pool whose objects scatter with negative binomial ``size`` (:func:`nb_lower_mean`; ``inf``, Poisson, by default). ``F(0) = 0`` exactly, and just above
    zero ``F`` is negative: ``D(lam) -> lam`` as ``lam -> 0``, so ``F(0+) ~ -rho * sum_{n_i>0} E_i``. ``F`` then rises
    without bound, since the first sum is eventually linear in ``rho`` and the second grows as ``sqrt(rho)``, so there
    is exactly one root above zero. ``F`` is not monotone, but bisection needs only ``F(lo) <= 0 <= F(hi)``, and
    because ``F(0+) < 0`` whenever any object has a count (the early return ensures one does), the lower end walks
    off zero on its own.

    The bisection closes the bracket to adjacent floats (:func:`~.numeric.bisect_to_adjacent`), so no iteration
    budget is needed.
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
        return float(np.clip(lam - n, 0.0, None).sum() - nb_lower_mean(lam, size).sum())

    lo, hi = 0.0, CONSTANTS.calibration.bracket_headroom * pooled
    if not (f(lo) <= 0.0 <= f(hi)):
        return GdnaDensityFit(pooled, pooled, int(n.size), total_e, bracket_ok=False)
    lo, hi = bisect_to_adjacent(lambda rho: not f(rho) > 0.0, lo, hi)
    return GdnaDensityFit(0.5 * (lo + hi), pooled, int(n.size), total_e, bracket_ok=True)


@dataclass(frozen=True, slots=True)
class BackgroundLaw:
    """The intergenic background's gDNA count law, the reference the intron abundance model reads an intron's count
    against: an object of gDNA opportunity ``E`` draws ``g ~ NegBinom(mean exp(log_rate)·E, size)``, except with
    prior probability ``departure``, when its gDNA is not at the background and any ``g`` in ``[0, C]`` is as likely.

    Fitted by :func:`background_law`, so that a region the law cannot explain does not set it. ``size`` joins the
    over-dispersion the pool supports and the rate's own posterior width, ``1/size = 1/s + 1/(sum w·n + ½)``, so an
    intron is one more draw from the law rather than a draw at a known rate. ``pooled_rate``, the naive ``sum n /
    sum E``, and ``kept``, the regions' summed weight under the law, report how much of the pool departed from it.
    The pool is fitted at the run's ``calibration.intron_departure_prior`` as ``departure``: it can measure its own
    spread but not how often a region leaves it. The introns read the law at their library's own ``departure`` where
    their strands measure one."""

    log_rate: float
    size: float
    n_objects: int
    total_exposure: float
    departure: float
    pooled_rate: float
    kept: float

    @property
    def informative(self) -> bool:
        """False when the pool has no opportunity to read a rate from; the model is then off."""
        return bool(self.total_exposure > CONSTANTS.calibration.opportunity_floor)


def _golden_max(f, lo: float, hi: float, width: float) -> float:
    """The maximiser of a unimodal ``f`` on ``[lo, hi]`` by golden-section search, to a bracket at most ``width``."""
    g = (5.0**0.5 - 1.0) / 2.0  # the golden section
    x1, x2 = hi - g * (hi - lo), lo + g * (hi - lo)
    f1, f2 = f(x1), f(x2)
    while hi - lo > width:
        if f1 >= f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - g * (hi - lo)
            f1 = f(x1)
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + g * (hi - lo)
            f2 = f(x2)
    return x1 if f1 >= f2 else x2


def background_law(counts, exposure, departure: float) -> BackgroundLaw:
    """The gDNA law of a pool of objects read as gDNA (the intergenic regions), from their counts and opportunity,
    fitted so that an object the law cannot explain (unannotated transcription, a copy-number change, a probe's reach)
    does not set it.

    Each object is read as the intron abundance model reads an intron: with probability ``1 - departure`` it is at
    the background and its count is a draw from ``NegBinom(rate·E, s)``, and otherwise any count is as likely,
    ``1/(n + 1)``. At each over-dispersion ``phi = 1/s`` the rate is the one-sided root (:func:`one_sided_rate`), which
    added fragments cannot move. ``phi`` is the value that maximises the pool's likelihood under that two-part law,
    searched by golden section between Poisson (``phi = 0``) and the pool's moment estimate, which a departing object
    only inflates. Each object's weight is its posterior probability of being at the background.
    """
    n = np.asarray(counts, dtype=np.float64).ravel()
    e = np.asarray(exposure, dtype=np.float64).ravel()
    if n.shape != e.shape:
        raise ValueError(f"counts {n.shape} and exposure {e.shape} must be aligned per object")
    c = CONSTANTS.calibration
    total_n, total_e = float(n.sum()), float(e.sum())
    pooled = total_n / total_e if total_e > 0.0 else 0.0
    if not total_e > c.opportunity_floor:
        return BackgroundLaw(-np.inf, np.inf, int(n.size), total_e, departure, pooled, 0.0)
    shape = float(c.background_rate_shape)
    if not total_n > c.mass_floor:  # no fragments: the Jeffreys posterior alone, Poisson
        return BackgroundLaw(
            float(np.log(shape) - np.log(total_e)),
            shape,
            int(n.size),
            total_e,
            departure,
            pooled,
            float(n.size),
        )
    keep, floor = float(np.log1p(-departure)), np.log(departure) - np.log(n + 1.0)

    def law_at(phi: float):
        size = 1.0 / phi if phi > 0.0 else np.inf
        rate = one_sided_rate(n, e, size).rate
        rate = rate if rate > 0.0 else pooled
        at_bg = keep + log_negbinom(n, rate * e, size)
        return float(np.logaddexp(at_bg, floor).sum()), rate, at_bg

    _, rate, _ = law_at(0.0)
    lam = rate * e
    sum_lam2 = float((lam * lam).sum())
    phi_hi = (
        (float(((n - lam) ** 2).sum()) - float(lam.sum())) / sum_lam2
        if sum_lam2 > c.mass_floor
        else 0.0
    )
    phi = 0.0
    if phi_hi > 0.0:
        found = _golden_max(
            lambda p: law_at(p)[0], 0.0, phi_hi, c.background_fit_tolerance * phi_hi
        )
        phi = found if law_at(found)[0] > law_at(0.0)[0] else 0.0
    _, rate, at_bg = law_at(phi)
    w = np.exp(at_bg - np.logaddexp(at_bg, floor))
    inv = phi + 1.0 / (float((w * n).sum()) + shape)
    return BackgroundLaw(
        float(np.log(rate)), 1.0 / inv, int(n.size), total_e, departure, pooled, float(w.sum())
    )
