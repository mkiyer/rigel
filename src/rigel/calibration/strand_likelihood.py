"""Strand likelihood — the Beta-Binomial-mixture strand log-likelihood of one node.

Used by the joint deconvolution (:mod:`joint_deconv`). A region's unspliced fragments are a
mix of gDNA (oriented-sense rate ½, unstranded) and RNA (oriented-sense rate ``rna_sense_frac``).
Given the observed (sense, antisense) split, this returns the strand log-likelihood over a grid
of candidate gDNA fractions ``gdna_frac``, which the joint deconvolution multiplies into the count prior.

Only POS/NEG nodes carry a valid sense split; AMBIG (both strands) and intergenic (no
transcript) nodes are strand-uninformative and left to the count clue.
"""

from __future__ import annotations

import numpy as np


def strand_loglik(
    gdna_frac: np.ndarray,
    sense: float,
    antisense: float,
    rna_sense_frac: float,
    *,
    strand_overdispersion: float = 0.0,
) -> np.ndarray:
    """Beta-Binomial-mixture strand log-likelihood of one node over a ``gdna_frac`` grid (§6).

    Of ``N = sense + antisense`` discrete fragments, a fraction ``gdna_frac`` are gDNA
    (oriented-sense rate 1/2) and ``1−gdna_frac`` are RNA (rate ``rna_sense_frac``). The
    per-fragment sense rate is Beta-distributed, so the count is over-dispersed: the variance is
    inflated by ``1 + strand_overdispersion·(N−1)`` (``0`` ⇒ Binomial). Normal moment
    approximation (the joint's fast path; the exact mixture is reserved for very small ``N``)::

        p   = (1/2)·gdna_frac + rna_sense_frac·(1−gdna_frac);  mean = N·p
        var = N·p·(1−p)·[1 + strand_overdispersion·(N−1)]
        loglik(gdna_frac) = −(1/2)·(sense − mean)² / var  −  (1/2)·log(var)
    """
    n = float(sense) + float(antisense)
    p = 0.5 * gdna_frac + rna_sense_frac * (1.0 - gdna_frac)
    mean = n * p
    var = n * p * (1.0 - p) * (1.0 + strand_overdispersion * max(n - 1.0, 0.0))
    var = np.maximum(var, 1.0e-9)
    return -0.5 * (float(sense) - mean) ** 2 / var - 0.5 * np.log(var)


__all__ = ["strand_loglik"]
