"""The locks' departure: what a locked object's two strands can share, estimated per library (layer 5).

An object whose annotation admits no RNA strand is locked: an intergenic region, or a boundary beside one. Every
unspliced fragment it holds is gDNA by structure. gDNA is double-stranded, so its fragments land on either genome strand
evenly, while a transcribed molecule lands on the plus column at κ or ``1 − κ`` by its strand. A locked object's two
columns therefore say whether its count can be gDNA: what one column holds beyond the other cannot, and at a perfect
protocol its gDNA is at most twice the smaller column, the part the two strands can share.

A stranded library measures it, in the introns' departure prior's two groups (:mod:`.intron_departure`), with the lock
in the background's place. Given the object's count ``n``, its plus-column count ``u`` reads its gDNA share ``f``
through ``u ~ Binom(n, p(f))``, the strand likelihood the solve itself uses::

    null       f = 1, every fragment gDNA:                       a = Binom(u; n, ½)
    departed   a share 1 − f is transcribed from one strand, f uniform on [0, 1], the strand either one:
               p_s(f) = f/2 + (1 − f)·q_s,   q₊ = κ,  q₋ = 1 − κ
               b_s = ∫₀¹ Binom(u; n, p_s(f)) df,   b = (b₊ + b₋)/2     differences of regularised incomplete betas

    L_j(ε) = (1 − ε)·a_j + ε·b_j

ε is the library's rate of departed locked objects, and ε̂ its maximum-likelihood value on ``[0, 1]``: the log-likelihood
is concave in ε (each term the log of an affine function), so its slope has one root, found by bisection to adjacent
floats, or the maximum is at an end. It has no prior and no setting, because ε = 0 is the lock as it stands.

A locked object is gDNA until its own strands say otherwise. The two groups' posterior, ``(1 − ε̂)·a / L(ε̂)``, would
move an object whose strands say nothing (``b = a``: a single fragment) to ``1 − ε̂`` on the library's rate alone. So
the weight the lock keeps is that posterior over its prior::

    w = min(1, a / L(ε̂)) = 1 / (1 + ε̂·(b/a − 1)₊)

the lock holds an object unless the library's law of locked objects explains its strands better than the lock does.
Where ``b ≤ a`` the object is all gDNA at any ε̂. Where its strands favour the departed reading, ε̂ tempers how far
they are believed: one one-sided object among thousands of balanced ones departs less than the same object where many
do, and at ε̂ = 1 it departs by ``a/b``. Its departed weight, ``1 − w``, is the part of the posterior its own strands
added to the prior, ``(P(departed | u) − ε̂)/(1 − ε̂)``, and never exceeds the posterior.

Its gDNA share is ``φ = w + (1 − w)·f̄``, ``f̄`` the departed posterior mean of ``f``; the rest, ``1 − φ``, is
transcribed, split between the strands by their posterior weights. The two parts of φ are different things. ``w`` is
the part of the object that is still a lock: structurally pure gDNA, a count. ``(1 − w)·f̄`` is the departed reading's
gDNA, an estimate with a width, as any single-strand slot's is. ``Var(log φ)`` is the spread of ``log f`` over the two groups, the
departed one read at quantiles of its posterior: a lock whose strands share about one fragment, give or take one, is
wider than the landscape's location floor and does not train it.

Four properties hold by construction. A lock whose own strands do not favour the departed reading is the lock, whatever
the library's rate. At κ = ½ every ``p_s`` is ½, ``b = a``, the likelihood is flat, and every share
is 1: an unstranded library is the lock. Under the null ``E[b/a] = 1``, so an all-gDNA library's slope at ε = 0 has mean
zero and ε̂ is 0 or of sampling size: binomial noise is not read as departure, as the bare shared part
``2·min(u, n − u)`` would read it. And transcription on both strands in equal amounts stays in the shared part: the
strands cannot see it, and it remains gDNA.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import betainc, betaincc, betainccinv, betaincinv
from scipy.stats import binom

from ..config import CONSTANTS
from .intron_departure import _score
from .numeric import beta_span, bisect_to_adjacent

__all__ = ["LockFit", "fit_lock_departure"]


@dataclass(frozen=True, slots=True)
class LockFit:
    """The locks' departure for one library, per object, aligned with the fit's inputs.

    ``departure`` is ε̂ and ``n_objects`` the locked objects holding a fragment. ``evidence`` is what the strands say
    for ε̂ over the lock, in nats: 0 where the lock stands. ``share`` is each object's gDNA share φ, ``kept`` the part
    of it that is still a lock, ``w``, ``rna_pos`` and ``rna_neg`` the transcribed remainder by strand (they sum to 1
    with the share), and ``variance`` is ``Var(log φ)``: 0 where the object is all gDNA."""

    departure: float
    n_objects: int
    evidence: float
    share: np.ndarray
    kept: np.ndarray
    rna_pos: np.ndarray
    rna_neg: np.ndarray
    variance: np.ndarray


def _quantiles(a, b, x1: float, x2: float, level):
    """``(objects, nodes)``: the quantiles at ``level`` of a ``Beta(a, b)`` truncated to ``[x1, x2]``; through the
    complements where the range sits in the upper tail, as :func:`~.numeric.beta_span` does."""
    lo, hi = betainc(a, b, x1)[:, None], betainc(a, b, x2)[:, None]
    c_lo, c_hi = betaincc(a, b, x1)[:, None], betaincc(a, b, x2)[:, None]
    lower = betaincinv(a[:, None], b[:, None], lo + level[None, :] * (hi - lo))
    upper = betainccinv(a[:, None], b[:, None], c_lo - level[None, :] * (c_lo - c_hi))
    return np.clip(np.where(lo > 0.5, upper, lower), x1, x2)


def _departed(plus, count, q: float):
    """One strand's departed reading, its RNA on the plus column at rate ``q``: the marginal likelihood ``b_s``, the
    posterior mean of the gDNA share ``f``, and the first two moments of ``log f``. With ``p = f/2 + (1 − f)·q`` uniform
    on its range, ``p``'s posterior is a ``Beta(u + 1, n − u + 1)`` truncated to that range: the mass and the mean are
    ratios of incomplete betas.

    ``log f`` is read at quantiles of that posterior. It is unbounded where ``f → 0``, the end of the range at
    ``p = q``, so the quantile level runs as ``s²`` from that end over midpoints ``s``, with weight ``2s``: the
    integrand ``2s·logᵐ f`` is then bounded, and the midpoint rule is exact to a part in a thousand at 64 nodes where
    even levels are five parts in a hundred low."""
    x1, x2 = min(q, 0.5), max(q, 0.5)
    a, b = plus + 1.0, count - plus + 1.0
    m0 = beta_span(a, b, x1, x2)
    live = m0 > 0.0
    mass = m0 / ((count + 1.0) * (x2 - x1))
    p1 = a / (count + 2.0) * beta_span(a + 1.0, b, x1, x2) / np.where(live, m0, 1.0)
    mean = np.where(live, np.clip((p1 - q) / (0.5 - q), 0.0, 1.0), 0.0)
    nodes = int(CONSTANTS.calibration.departure_quadrature_nodes)
    s = (np.arange(nodes) + 0.5) / nodes
    level = s * s if q < 0.5 else 1.0 - s * s  # from the end where f = 0
    f = (_quantiles(a, b, x1, x2, level) - q) / (0.5 - q)
    with np.errstate(divide="ignore"):
        log_f = np.log(np.clip(f, 0.0, 1.0))
    weight = 2.0 * s / nodes
    first, second = (log_f * weight).sum(axis=1), (log_f**2 * weight).sum(axis=1)
    return mass, mean, np.where(live, first, 0.0), np.where(live, second, 0.0)


def _slope(d, eps: float) -> float:
    """``d/dε Σ log L_j`` from each object's ``d = log b − log a``, the introns' score summed; ``±∞`` at an end where a
    term's likelihood vanishes."""
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        return float(_score(d, eps).sum())


def _departure(d) -> float:
    """ε̂: the root of the log-likelihood's slope on ``(0, 1)``, or the end it rises or falls to."""
    if not _slope(d, float(np.finfo(np.float64).tiny)) > 0.0:
        return 0.0
    if _slope(d, 1.0) >= 0.0:
        return 1.0
    lo, hi = bisect_to_adjacent(lambda eps: _slope(d, eps) > 0.0, 0.0, 1.0)
    return 0.5 * (lo + hi)


def _gain(d, eps: float) -> float:
    """``Σ log(L_j(ε)/a_j)``, the strands' log-likelihood at ε over the lock: ``log((1 − ε) + ε·e^d)`` per object,
    written so that neither sign of ``d`` overflows."""
    e = np.exp(-np.abs(d))
    with np.errstate(divide="ignore"):
        term = np.where(d > 0.0, d + np.log(eps + (1.0 - eps) * e), np.log1p(eps * (e - 1.0)))
    return float(term.sum())


def fit_lock_departure(plus, minus, kappa: float) -> LockFit:
    """The locks' departure for a library's locked objects: their unspliced counts on the two genome-strand columns and
    the library's RNA sense fraction κ, ½ where the strand channel is dead. Returns :class:`LockFit`, every array
    aligned with the inputs; an object holding no fragment keeps the lock."""
    u = np.asarray(plus, np.float64).ravel()
    v = np.asarray(minus, np.float64).ravel()
    if u.shape != v.shape:
        raise ValueError(f"plus {u.shape} and minus {v.shape} must be aligned per locked object")
    if np.any(u != np.floor(u)) or np.any(v != np.floor(v)):
        # a fraction has no binomial probability
        raise ValueError("the strand counts must be whole fragments, as the scan tallies them")
    size = u.shape[0]
    share, kept, variance = np.ones(size), np.ones(size), np.zeros(size)
    rna_pos, rna_neg = np.zeros(size), np.zeros(size)
    held = (u + v) > 0.0
    n_objects = int(held.sum())
    if n_objects == 0 or kappa == 0.5:  # no fragment or no strand: the lock stands
        return LockFit(0.0, n_objects, 0.0, share, kept, rna_pos, rna_neg, variance)
    plus_held, count = u[held], (u + v)[held]
    log_a = binom.logpmf(plus_held, count, 0.5)
    b_pos, f_pos, l1_pos, l2_pos = _departed(plus_held, count, float(kappa))
    b_neg, f_neg, l1_neg, l2_neg = _departed(plus_held, count, 1.0 - float(kappa))
    total = b_pos + b_neg
    with np.errstate(divide="ignore"):
        # −∞ where the departed reading is impossible: the lock stands there
        d = np.log(0.5 * total) - log_a
    eps = _departure(d)
    if eps == 0.0:
        return LockFit(0.0, n_objects, 0.0, share, kept, rna_pos, rna_neg, variance)
    with np.errstate(over="ignore"):
        # gDNA until proven otherwise: the lock's posterior over its prior, min(1, a/L), so the library's rate tempers
        # a lock its own strands move and never moves one by itself
        w = 1.0 / (1.0 + eps * np.expm1(np.maximum(d, 0.0)))
    safe = np.where(total > 0.0, total, 1.0)
    s_pos, s_neg = b_pos / safe, b_neg / safe  # the departed reading's strand, given the counts
    mean = s_pos * f_pos + s_neg * f_neg

    def over_strands(pos, neg):
        # a strand's reading the counts rule out carries no weight, whatever its own moments are
        return np.where(s_pos > 0.0, s_pos * pos, 0.0) + np.where(s_neg > 0.0, s_neg * neg, 0.0)

    with np.errstate(invalid="ignore"):
        # log f over the two groups: 0 where the object is all gDNA, the departed posterior's own log f otherwise
        first = (1.0 - w) * over_strands(l1_pos, l1_neg)
        second = (1.0 - w) * over_strands(l2_pos, l2_neg)
        var = np.where(w < 1.0, np.maximum(second - first * first, 0.0), 0.0)
        var = np.where(np.isfinite(var), var, np.inf)
    share[held] = w + (1.0 - w) * mean
    kept[held] = w
    rna_pos[held] = (1.0 - w) * s_pos * (1.0 - f_pos)
    rna_neg[held] = (1.0 - w) * s_neg * (1.0 - f_neg)
    variance[held] = var
    return LockFit(
        float(eps), n_objects, max(_gain(d, eps), 0.0), share, kept, rna_pos, rna_neg, variance
    )
