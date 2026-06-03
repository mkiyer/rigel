"""Phase 2 — the strand clue: per-node gDNA-fraction posterior from the sense/antisense split.

A region's unspliced fragments are a mix of gDNA (sense rate κ_d = 0.5, unstranded) and RNA
(sense rate κ_rna). Given the observed oriented split (S sense, A antisense; N = S+A), infer
the gDNA fraction ``π_g`` as a **posterior**, not a point estimate, so that **thin counts
stay diffuse** (the review's Risk A): with a flat prior, a region with 1–2 fragments yields
a wide posterior centred near 0.5 rather than a confident 0/1 call.

Model (oriented so the RNA sense rate is κ_rna for both strands; ``_sense_flux`` convention):

    P(oriented-sense | π_g) = ½·π_g + κ_rna·(1 − π_g)
    S ~ Binomial(N, P)              flat Beta prior on π_g

The Beta-Binomial **supporting distribution** enters as an optional strand overdispersion
``ρ`` that caps the effective sample size ``N_eff = N / (1 + ρ(N−1))`` (so the strand call
cannot become arbitrarily certain beyond ~1/ρ fragments). ``ρ = 0`` is the Binomial limit.

Output is the posterior mean ``π_g`` and its variance (the joint decode, Phase 3, weights
the two clues by precision = 1/variance). Only **POS/NEG** regions are strand-decodable;
AMBIG has no valid sense split, and intergenic (NONE) has no transcript — both are left to
the count clue / sweep.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .signature import TS_NEG, TS_POS

# Machine-scale guard keeping the grid off the {0,1} singularities of log P.
_PI_EPS = 1.0e-4
# Uniform-distribution variance — the posterior variance of a fully uninformative node.
_UNIFORM_VAR = 1.0 / 12.0


@dataclass(frozen=True, slots=True)
class StrandDecode:
    """Per-region gDNA-fraction posterior from the strand clue."""

    pi_g: np.ndarray  # float64[R] — posterior-mean gDNA fraction
    pi_g_var: np.ndarray  # float64[R] — posterior variance (precision = 1/var for the joint)
    decodable: np.ndarray  # bool[R] — strand-decodable (POS/NEG, identifiable κ_rna)
    n_decodable: int


def strand_decode(
    substrate,
    kappa_rna: float,
    *,
    strand_overdispersion: float = 0.0,
    n_grid: int = 200,
) -> StrandDecode:
    """Posterior-mean gDNA fraction per region from the contained sense/antisense split."""
    ts = np.asarray(substrate.ts_class)
    c = substrate.contained
    pos = c.n_unspliced_pos.astype(np.float64)
    neg = c.n_unspliced_neg.astype(np.float64)
    # Orient to transcript sense (cf. estep._sense_flux): NEG → genome '−', else '+'.
    sense = np.where(ts == TS_NEG, neg, pos)
    n = pos + neg
    antisense = n - sense
    r = ts.shape[0]

    pi_g = np.full(r, 0.5, dtype=np.float64)
    pi_g_var = np.full(r, _UNIFORM_VAR, dtype=np.float64)
    decodable = (ts == TS_POS) | (ts == TS_NEG)

    # Unstranded library → the strand clue carries no information (Binomial degenerate).
    if abs(2.0 * kappa_rna - 1.0) < _PI_EPS:
        return StrandDecode(pi_g, pi_g_var, np.zeros(r, dtype=bool), 0)

    grid = np.linspace(_PI_EPS, 1.0 - _PI_EPS, n_grid)
    p_plus = 0.5 * grid + kappa_rna * (1.0 - grid)  # P(oriented-sense | π_g)
    log_p = np.log(p_plus)
    log_1mp = np.log1p(-p_plus)

    for i in np.flatnonzero(decodable):
        ni = n[i]
        if ni < 0.5:
            continue  # no fragments → stays at the uninformative prior (0.5, wide)
        # Beta-Binomial cap on the effective sample size (ρ = 0 ⇒ Binomial).
        scale = 1.0 / (1.0 + strand_overdispersion * (ni - 1.0))
        loglik = scale * (sense[i] * log_p + antisense[i] * log_1mp)  # flat prior
        w = np.exp(loglik - loglik.max())
        w /= w.sum()
        mean = float(np.dot(grid, w))
        pi_g[i] = mean
        pi_g_var[i] = float(np.dot((grid - mean) ** 2, w))

    return StrandDecode(
        pi_g=pi_g, pi_g_var=pi_g_var, decodable=decodable, n_decodable=int(decodable.sum())
    )


def strand_loglik(
    pi_g: np.ndarray,
    sense: float,
    antisense: float,
    kappa_rna: float,
    *,
    strand_overdispersion: float = 0.0,
) -> np.ndarray:
    """Beta-Binomial-mixture strand log-likelihood of one node over a ``π_g`` grid (§6).

    Of ``N = sense + antisense`` discrete fragments, a fraction ``π_g`` are gDNA
    (oriented-sense rate ½) and ``1−π_g`` are RNA (rate ``κ_rna``). The per-fragment sense
    rate is Beta-distributed, so the count is over-dispersed: the variance is inflated by
    ``1 + ρ(N−1)`` (``ρ = 0`` ⇒ Binomial). Normal moment approximation (the joint's fast
    path; the exact mixture is reserved for very small ``N``)::

        p(π_g) = ½·π_g + κ_rna·(1−π_g);  μ = N·p;  σ² = N·p(1−p)·[1 + ρ(N−1)]
        loglik(π_g) = −½ (sense − μ)² / σ²  −  ½ log σ²
    """
    n = float(sense) + float(antisense)
    p = 0.5 * pi_g + kappa_rna * (1.0 - pi_g)
    mu = n * p
    var = n * p * (1.0 - p) * (1.0 + strand_overdispersion * max(n - 1.0, 0.0))
    var = np.maximum(var, 1.0e-9)
    return -0.5 * (float(sense) - mu) ** 2 / var - 0.5 * np.log(var)


__all__ = ["StrandDecode", "strand_decode", "strand_loglik"]
