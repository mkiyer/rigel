"""The introns' departure prior ε, estimated per library from their antisense counts (layer 5).

An intron's gDNA is at the intergenic background (``gdna_density.BackgroundLaw``) with probability ``1 − ε``; with
probability ε it departs (a probe's reach, a copy-number change, a mappability gap), and then any split of its count is
as likely. ε caps how far the background can argue against an intron's own strand evidence, and it is a base rate,
the library's own: in a capture library many introns beside a probed exon depart, in a plain one few do.

A stranded library measures it. On a single-strand intron gDNA lands on either strand evenly while RNA lands on the
antisense strand at ``1 − κ``, so the antisense count ``A``, given the intron's count ``C``, reads the intron's gDNA
share ``f`` through ``A ~ Binom(C, p(f))``, ``p(f) = f/2 + (1 − f)(1 − κ)``: the strand likelihood the solve itself
uses. At the background, ``f = g/C`` with ``g`` the background's count truncated to ``0…C``; departed, ``f`` is uniform
on ``[0, 1]``. Each intron's evidence is::

    L_i(ε) = (1 − ε)·a_i + ε·b_i
    a_i = E[Binom(A_i; C_i, p(g/C_i)) | g ≤ C_i]      the background, at its own midpoint quantiles
    b_i = ∫₀¹ Binom(A_i; C_i, p(f)) df                 departed: a difference of regularised incomplete betas

Both terms are conditioned on ``C``, so the estimate reads the strand alone: an intron's count, however far its RNA
puts it above the background, is no evidence here, and at κ = ½ every ``p(f) = ½``, ``a_i = b_i``, and the likelihood
is flat.

``a_i`` is the background's reading of the intron as the abundance row makes it, so ε̂ is the weight of the row's flat
part that the strands support: how often the background misreads an intron, for any reason. A departure of the
intron's gDNA is one reason. A sparse intron is the other: where the background expects less than one fragment, every
fragment the intron holds is an excess, and the row, which charges nothing for RNA, reads the excess as RNA. If those
fragments are gDNA, the strands say so, so a library whose sparse introns hold mostly gDNA reads near 1, and its rows
step aside for the strands.

The prior is ``Beta(2, 1/m0)``, ``m0`` the run's ``calibration.intron_departure_prior``. Its mode is ``m0`` and it holds
one departing intron's worth of evidence at that rate, so it says "about m0, within a few-fold". A library with no
informative intron returns ``m0``; one with thousands is decided by its data. ε̂ is the posterior mode: the log
posterior is concave in ε (each term the log of an affine function), so its derivative has one root, found by
bisection to adjacent floats.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from scipy.stats import binom, nbinom, poisson

from ..config import CONSTANTS
from .gdna_density import BackgroundLaw
from .numeric import beta_span, bisect_to_adjacent

__all__ = ["DepartureFit", "fit_departure_prior"]


@dataclass(frozen=True, slots=True)
class DepartureFit:
    """The introns' departure prior for one library, with the evidence it rests on.

    ``departure`` is ε̂, the posterior mode of a prior centred on the run setting. ``n_introns`` counts the
    single-strand introns holding a fragment. ``evidence`` is what the data say for ε̂ over the setting, in nats: the
    log-likelihood of the introns' strands at ε̂ less that at the setting. Near 0 the setting explains the strands as
    well, and the estimate is the prior's; at tens of nats the library decides."""

    departure: float
    n_introns: int
    evidence: float


def _background_quantiles(mu, size, count, nodes: int):
    """``(introns, nodes)``: midpoint quantiles of each intron's background count, truncated to ``0…count``."""
    u = (np.arange(nodes) + 0.5) / nodes
    dist = nbinom(size, size / (size + mu[:, None])) if np.isfinite(size) else poisson(mu[:, None])
    top = dist.cdf(count[:, None])
    return np.minimum(dist.ppf(u[None, :] * top), count[:, None])


def departure_evidence(sense, antisense, exposure, law: BackgroundLaw, kappa: float):
    """``log a_i, log b_i`` per intron: its antisense count's likelihood given its count, at the background and
    departed. Every intron must carry a fragment."""
    a = np.asarray(antisense, np.float64)
    c = np.asarray(sense, np.float64) + a
    mu = np.exp(law.log_rate) * np.asarray(exposure, np.float64)
    p_rna, p_gdna = 1.0 - kappa, 0.5  # p(f) runs from 1 − κ (all RNA) to ½ (all gDNA)
    g = _background_quantiles(
        mu, law.size, c, int(CONSTANTS.calibration.departure_quadrature_nodes)
    )
    f = g / c[:, None]
    log_a = binom.logpmf(a[:, None], c[:, None], p_gdna * f + p_rna * (1.0 - f))
    log_a = np.logaddexp.reduce(log_a, axis=1) - np.log(g.shape[1])
    # ∫ Binom(A; C, p) dp over [x1, x2] is (I_x2 − I_x1)(A + 1, C − A + 1)/(C + 1); over the span x2 − x1 it is the
    # uniform f's average. The upper tail goes through the complements, which keep their digits where both I near 1.
    x1, x2 = min(p_rna, p_gdna), max(p_rna, p_gdna)
    mass = beta_span(a + 1.0, c - a + 1.0, x1, x2)
    with np.errstate(divide="ignore"):
        log_b = np.log(mass) - np.log(c + 1.0) - np.log(x2 - x1)
    return log_a, log_b


def _log_mix(d, eps):
    """``log L_i(ε)`` less ``log max(a_i, b_i)``, from ``d = log b_i − log a_i``. It is finite at every ``d``, ``±∞``
    included, and the reference cancels in a difference between two ε: the intron's log-likelihood ratio."""
    e = np.exp(-np.abs(d))
    return np.where(d > 0.0, np.log(eps + (1.0 - eps) * e), np.log1p(-eps * (1.0 - e)))


def _score(d, eps):
    """``d/dε log L_i``, from ``d = log b_i − log a_i``, stable at either extreme of ``d``."""
    ed = np.exp(-np.abs(d))
    return np.where(
        d > 0.0, (1.0 - ed) / ((1.0 - eps) * ed + eps), (ed - 1.0) / ((1.0 - eps) + eps * ed)
    )


def fit_departure_prior(
    sense, antisense, exposure, law: BackgroundLaw, kappa: float, prior_mode: float
):
    """ε̂ for a library's single-strand introns: their sense and antisense unspliced counts (antisense relative to each
    intron's annotated strand) and gDNA opportunity, the background law, the library's RNA sense fraction κ, and the
    prior's mode, the run setting. Returns :class:`DepartureFit`."""
    m0 = float(prior_mode)
    n0 = 1.0 / m0
    s = np.asarray(sense, np.float64).ravel()
    a = np.asarray(antisense, np.float64).ravel()
    e = np.asarray(exposure, np.float64).ravel()
    if not (s.shape == a.shape == e.shape):
        raise ValueError(
            f"sense {s.shape}, antisense {a.shape} and exposure {e.shape} must be aligned per intron"
        )
    if np.any(s != np.floor(s)) or np.any(a != np.floor(a)):
        # a fraction has no binomial probability, so its intron would read as impossible at the background
        raise ValueError(
            "the sense and antisense counts must be whole fragments, as the scan tallies them"
        )
    live = (s + a) > 0.0
    s, a, e = s[live], a[live], e[live]
    n = int(s.size)
    if (
        n == 0 or kappa == 0.5 or not law.informative
    ):  # no fragment, no strand, or no background: no evidence
        return DepartureFit(m0, n, 0.0)
    log_a, log_b = departure_evidence(s, a, e, law, kappa)
    fa, fb = np.isfinite(log_a), np.isfinite(log_b)
    # d = log b − log a; an intron impossible under both carries no evidence about ε
    d = np.where(fa & fb, log_b - log_a, np.where(fb, np.inf, -np.inf))[fa | fb]
    if not np.any(d != 0.0):
        return DepartureFit(m0, n, 0.0)

    # the log posterior's derivative, Σ score_i + (α − 1)/ε − (β − 1)/(1 − ε) with α − 1 = 1 and β − 1 = n0 − 1,
    # falls from +∞ to −∞ across (0, 1)
    def slope(eps):
        return float(_score(d, eps).sum()) + 1.0 / eps - (n0 - 1.0) / (1.0 - eps)

    lo, hi = bisect_to_adjacent(lambda eps: slope(eps) > 0.0, 0.0, 1.0)
    eps = 0.5 * (lo + hi)
    # The mode gains at least the likelihood the prior gives up by leaving m0, so the data's gain is not negative; a
    # rounding below 0 is 0.
    evidence = float((_log_mix(d, eps) - _log_mix(d, m0)).sum())
    return DepartureFit(eps, n, max(evidence, 0.0))
