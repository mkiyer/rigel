"""The simplex per-node solver — partition a node's unspliced mass into ``(f_rna₊, f_rna₋, f_g)``.

Spec: ``docs/calibration/propagation_simplex_plan.md`` (rev-2). Calibration models **only RNA vs gDNA**
(the per-locus EM separates nascent from mature downstream), so every node's unspliced mass is one pie
cut into three slices ``(f_rna₊, f_rna₋, f_g)`` summing to 1 — **no nascent slice, no mature slice**.
Evidence *pushes* the point on the 2-simplex; the answer is always a normalized partition, so
over-subtraction is structurally impossible and an under-constrained slice stays wide ("unknown").

This module is the **per-node solve only** (plan §3/§7). Its ``_mixture_strand_loglik`` is the reusable
three-component strand-likelihood primitive the ``simplex_logodds`` log-density solve builds on.

The strand term is the **three-component generalization** of :func:`strand_likelihood.strand_loglik`:
of ``N`` unspliced fragments a fraction ``f_g`` are gDNA (genomic-plus rate ½, overdispersion ``od_g``),
``f₊`` are ``+``-strand RNA (plus rate ``κ``, od ``od_r``) and ``f₋`` are ``−``-strand RNA (plus rate
``1−κ``, od ``od_r``). The mixture plus-strand probability is ``p = ½·f_g + κ·f₊ + (1−κ)·f₋ = ½ +
(κ−½)·(f₊−f₋)`` — it depends **only on the tilt** ``t = f₊−f₋`` (an identity, not an approximation), so
the strand constrains ``t`` and the count/prior constrain the RNA-vs-gDNA magnitude orthogonally. For a
single-strand node (``f₋≡0``) every term collapses back to ``strand_loglik`` exactly — the no-regression
guard.

Evidence terms (each a precision-weighted log term; plan §3):
  * **strand** — the 3-component mixture loglik above (precision = BB curvature ``N·(2κ−1)²``, intrinsic).
  * **count** — Gaussian pull of ``f_g`` toward ``count_gdna_frac`` with ``count_precision`` (the
    propagated count evidence; Poisson floor to start, loess/NB ``var~mean`` to follow — plan §4).
  * **spliced RNA lower bound (sided)** — a soft one-sided push ``f_s·U ≥ spliced_s`` (RNA can't be below
    the mature its junctions directly observe); precision = the spliced flux (Poisson). Single-exon ⇒ no
    junction ⇒ no push (correct: that RNA is honestly "unknown").
  * **gDNA prior** — a weak Dirichlet pseudo-count on ``f_g`` (RNA prior flat / "earned"); keeps ``f_g``
    off zero with no evidence and degrades gracefully at κ≈½. ``gdna_prior_count`` is the one sanctioned
    starting constant (placeholder pending derivation; plan §5).

The one-sided spliced shape (clipped-Gaussian) and the Dirichlet form are the toy-sweep-calibratable
*shape* choices flagged in plan §10 — no magic widths (every precision is a derived flux count).
"""

from __future__ import annotations

import numpy as np

__all__ = ["_mixture_strand_loglik"]

_EPS = 1.0e-9


def _mixture_strand_loglik(u_pos, n, f_g, f_pos, f_neg, kappa, od_g, od_r):
    """Three-component gDNA/RNA₊/RNA₋ strand loglik — :func:`strand_loglik` generalized to two RNA strands.

    Broadcasts ``(u_pos, n)`` of shape ``(nodes, 1)`` against the lattice ``(f_*)`` of shape
    ``(1, P)`` → ``(nodes, P)``. Mean ``N·p`` with ``p = ½·f_g + κ·f₊ + (1−κ)·f₋``; variance is the
    Binomial-mixture term plus each component's overdispersion excess scaled by its own ``μ(1−μ)``
    (gDNA ¼, RNA κ(1−κ) for both strands), exactly as in the two-component model.
    """
    p = 0.5 * f_g + kappa * f_pos + (1.0 - kappa) * f_neg
    mean = n * p
    rscale = kappa * (1.0 - kappa)  # κ(1−κ): each RNA strand's μ(1−μ)
    var = (
        n * p * (1.0 - p)
        + (n * f_g) ** 2 * 0.25 * od_g
        + (n * f_pos) ** 2 * rscale * od_r
        + (n * f_neg) ** 2 * rscale * od_r
    )
    var = np.maximum(var, _EPS)
    return -0.5 * (u_pos - mean) ** 2 / var - 0.5 * np.log(var)
