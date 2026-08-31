"""rigel.calibration.gdna_density — the gDNA background RATE, from counts and their opportunity.

Two estimators of one quantity, over one signature ``(counts, exposure) -> rate``, so there is exactly
one implementation of each and no second home to drift:

* :func:`pooled_log_rate` — ``ln((sum g + shape) / sum E)``, the Gamma-posterior location. It assumes the
  pool it is handed is pure gDNA, which is why it is the naive one.
* :func:`one_sided_rate` — the same rate WITHOUT that assumption, for a pool that transcription may have
  added to. This is the one the fragment-length model uses.

⭐⭐ **THE ONE-SIDED ARGUMENT, which is why the second exists.** No structural class of objects is pure
gDNA — "intergenic" is whatever the annotation leaves over, nascent RNA sits inside introns by
definition, and pervasive transcription is real. But **transcription can only ADD fragments to an object,
never remove them**, so every object's observed density ``n_i / E_i`` is an OVERESTIMATE and the clean
objects are the ones on the LOW side. Using the side a contaminant cannot reach is an assumption about
DIRECTION, which the biology gives, in place of an assumption about PURITY, which it does not.

**The estimator.** With ``lam_i = rho * E_i``, under pure gDNA ``n_i ~ Poisson(lam_i)`` and ``E[n_i - lam_i]
= 0``, so the negative part of the residual has a closed form — De Moivre's mean-absolute-deviation
identity for the Poisson::

    E[(lam - N)+]  =  1/2 E|N - lam|  =  lam^(floor(lam)+1) * exp(-lam) / floor(lam)!   =:  D(lam)

Contamination only RAISES ``n_i``, hence only LOWERS ``(lam_i - n_i)+``. Matching the observed one-sided
mass to its exact null expectation is then one scalar equation::

    F(rho)  =  sum_i (rho*E_i - n_i)+  -  sum_i D(rho*E_i)  =  0

⭐ **No trim depth, no quantile, no threshold, no chosen level.** An earlier form of this estimator kept
the lowest ``q`` of objects by density and pooled those; it worked, but the safe ``q`` is a function of the
library's contaminated fraction — a property of the sample, not of the tool — and it broke at ``q = 0.80``.
The identity above is that estimator's limit taken exactly.

⚠ **The one bias, and it is one-signed.** A heavily contaminated object contributes ``0`` to the first sum
but its full ``D(lam_i)`` to the second, so the root is pushed UP. Measured against an origin-split oracle
it runs **+0.2 % to +6.7 %** high; it is never low. ⛔ Its scaling with the contaminated fraction is NOT
established, so a library dirtier than any simulated panel may sit outside that range.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

__all__ = [
    "GdnaDensityFit",
    "contained_opportunity",
    "one_sided_rate",
    "poisson_lower_mean",
    "pooled_log_rate",
    "region_lengths_from_partition",
]

#: Bisection steps. 200 halvings take a float64 bracket below its own resolution, so this is a
#: termination guarantee rather than a tuned iteration count.
_BISECT_STEPS = 200

#: The bracket's upper end as a multiple of the pooled rate. The one-sided root is always BELOW the
#: pooled rate (contamination only inflates it), so any multiple above 1 brackets it; the value is a
#: headroom guard, and ``bracket_ok`` reports if it ever failed to.
_BRACKET_HEADROOM = 10.0


@dataclass(frozen=True, slots=True)
class GdnaDensityFit:
    """A gDNA density estimate with the diagnostics that make it auditable.

    ``rate`` is the one-sided estimate; ``pooled_rate`` is the naive one the pool would have given. ⭐
    ``rate_over_pooled`` is the number to read: it is how much contamination the fit found, and a value
    near 1 means the pool was already clean rather than that the estimator did nothing.
    """

    rate: float
    pooled_rate: float
    n_objects: int
    total_counts: float
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

    De Moivre's mean-absolute-deviation identity for the Poisson, halved (``E[N - lam] = 0`` makes the two
    one-sided means equal). ⭐ Closed form at every ``lam``, so the estimator that uses it needs no
    simulation and no tabulation. Verified against the exact truncated sum and against Monte Carlo.
    """
    lam = np.asarray(lam, dtype=np.float64)
    out = np.zeros(lam.shape, dtype=np.float64)
    pos = lam > 0.0
    if not np.any(pos):
        return out
    lp = lam[pos]
    k = np.floor(lp)
    # in log space: a large `lam` overflows `lam**(k+1)` and `k!` separately while their ratio is O(sqrt(lam))
    out[pos] = np.exp((k + 1.0) * np.log(lp) - lp - gammaln(k + 1.0))
    return out


def pooled_log_rate(counts, exposure, *, shape: float = 0.0) -> float:
    """``ln((sum counts + shape) / sum exposure)`` — the naive pooled rate, in logs.

    ⛔ **It assumes the pool is pure**, which no structural class of objects is; that is what
    :func:`one_sided_rate` exists to avoid. It stays here because it is the other estimator of the same
    quantity and both belong in one module.

    ``shape`` is the caller's prior pseudo-count (a Gamma-posterior location uses the Jeffreys ``½``);
    the default of 0 is the plain ratio. Returns ``-inf`` when the pool has no support, which is the
    honest value rather than a raised error: a chain with no intergenic region is a real state.
    """
    se = float(np.asarray(exposure, dtype=np.float64).sum())
    if not se > 0.0:
        return -np.inf
    return float(np.log(float(np.asarray(counts, dtype=np.float64).sum()) + shape) - np.log(se))


def one_sided_rate(counts, exposure) -> GdnaDensityFit:
    """The gDNA density from the LOW side of the per-object density distribution.

    Solves ``F(rho) = sum (rho*E_i - n_i)+ - sum D(rho*E_i) = 0`` by bisection.

    ⭐ **The bracket is structural, not guessed.** ``F(0) = 0`` exactly (both sums vanish), and just above
    zero ``F`` is strictly negative — ``D(lam) -> lam`` as ``lam -> 0``, so the second sum approaches
    ``rho * sum E`` while the first collects only the objects with ``n_i = 0``, leaving
    ``F(0+) ~ -rho * sum_{n_i>0} E_i``. ``F`` then rises without bound, because the first sum is eventually
    linear in ``rho`` while the second grows only as ``sqrt(rho)``. So there is exactly one root above
    zero.

    ⚠ ``F`` is NOT monotone — it dips before it rises — but bisection does not need it to be, only
    ``F(lo) <= 0 <= F(hi)``. ⭐ **And the root at ``rho = 0`` is not sticky**: because ``F(0+) < 0``
    whenever any object carries a count, the first midpoints test negative and the lower end walks off
    zero on its own. The early return above is what guarantees that premise, so no separate guard against
    the boundary root is needed — one was written, and removed when no perturbation could make it fire.
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
        # answer. ⭐ A zero-gDNA library reaches this and must NOT be handed a fabricated rate.
        return GdnaDensityFit(0.0, pooled, int(n.size), total_n, total_e, bracket_ok=False)

    def f(rho: float) -> float:
        lam = rho * e
        return float(np.clip(lam - n, 0.0, None).sum() - poisson_lower_mean(lam).sum())

    lo, hi = 0.0, _BRACKET_HEADROOM * pooled
    if not (f(lo) <= 0.0 <= f(hi)):
        return GdnaDensityFit(pooled, pooled, int(n.size), total_n, total_e, bracket_ok=False)
    for _ in range(_BISECT_STEPS):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return GdnaDensityFit(0.5 * (lo + hi), pooled, int(n.size), total_n, total_e, bracket_ok=True)


def contained_opportunity(pmf, lengths) -> np.ndarray:
    """``E_f[(ell - L + 1)+]`` per object — the CONTAINED opportunity a uniformly placed fragment has.

    ⭐ Evaluated by cumulative sums rather than an ``(objects x lengths)`` outer product, which at genome
    scale would be tens of billions of entries. Writing ``F`` and ``S`` for the pmf's cumulative mass and
    cumulative first moment, ``E[(ell-L+1)+] = (ell+1)*F(ell) - S(ell)``, and for ``ell`` beyond the pmf's
    support that is exactly ``ell + 1 - E[L]``. Cost is ``O(objects + max_length)``.
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


def region_lengths_from_partition(region_bounds, ref_pos_offsets, n_regions: int) -> np.ndarray:
    """Region lengths from the scanner's own partition, differenced **PER REFERENCE**.

    ⛔ ``np.diff(region_bounds)`` straight through is WRONG. The bound positions are concatenated per
    reference (``ref_pos_offsets`` delimits them), and each reference's ``k`` regions contribute ``k + 1``
    positions — so a plain diff manufactures one phantom region spanning the junction between every pair
    of adjacent references. It is silent, plausible-looking, and lands on whichever reference happens to
    follow; caught here once by a gate rather than by a wrong answer.
    """
    b = np.asarray(region_bounds, dtype=np.float64).ravel()
    off = np.asarray(ref_pos_offsets, dtype=np.int64).ravel()
    out = np.zeros(int(n_regions), dtype=np.float64)
    at = 0
    for r in range(off.size - 1):
        lo, hi = int(off[r]), int(off[r + 1])
        k = hi - lo - 1  # k regions from k+1 bound positions
        if k <= 0:
            continue
        out[at : at + k] = np.diff(b[lo:hi])
        at += k
    return out
