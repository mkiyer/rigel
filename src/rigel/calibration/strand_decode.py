"""Strand likelihood — the Beta-Binomial-mixture strand log-likelihood of one node.

Used by the joint deconvolution (:mod:`joint_decode`). A region's unspliced fragments are a
mix of gDNA (oriented-sense rate ½, unstranded) and RNA (oriented-sense rate ``rna_sense_frac``).
Given the observed (sense, antisense) split, this returns the strand log-likelihood over a grid
of candidate gDNA fractions ``π_g``, which the joint decode multiplies into the count prior.

Only POS/NEG nodes carry a valid sense split; AMBIG (both strands) and intergenic (no
transcript) nodes are strand-uninformative and left to the count clue.
"""

from __future__ import annotations

import numpy as np


def strand_loglik(
    pi_g: np.ndarray,
    sense: float,
    antisense: float,
    rna_sense_frac: float,
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
    p = 0.5 * pi_g + rna_sense_frac * (1.0 - pi_g)
    mu = n * p
    var = n * p * (1.0 - p) * (1.0 + strand_overdispersion * max(n - 1.0, 0.0))
    var = np.maximum(var, 1.0e-9)
    return -0.5 * (float(sense) - mu) ** 2 / var - 0.5 * np.log(var)


__all__ = ["strand_loglik"]
