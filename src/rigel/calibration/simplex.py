"""The three-component strand-likelihood primitive shared by the log-density per-node solve.

Calibration models **only RNA vs gDNA** (the per-locus EM separates nascent from mature downstream), so every
node's unspliced mass is one pie ``(f_rna₊, f_rna₋, f_g)`` summing to 1 — no nascent/mature slice.

This module holds only :func:`_mixture_strand_loglik` — the reusable three-component strand-likelihood the
``simplex_logodds`` log-density solve builds on. It is the **three-component generalization** of
:func:`strand_likelihood.strand_loglik`: of ``N`` unspliced fragments a fraction ``f_g`` are gDNA (plus rate ½,
overdispersion ``od_g``), ``f₊`` are ``+``-strand RNA (plus rate ``κ``, od ``od_r``), ``f₋`` are ``−``-strand
RNA (plus rate ``1−κ``, od ``od_r``). The mixture plus-strand probability ``p = ½·f_g + κ·f₊ + (1−κ)·f₋ = ½ +
(κ−½)·(f₊−f₋)`` depends **only on the tilt** ``t = f₊−f₋`` (an identity), so the strand constrains ``t`` and the
count/prior constrain the RNA-vs-gDNA magnitude orthogonally. A single-strand node (``f₋≡0``) collapses to
``strand_loglik`` exactly (the no-regression guard).
"""

from __future__ import annotations

import numpy as np

__all__ = ["_mixture_strand_loglik"]

_EPS = 1.0e-9


def _mixture_strand_loglik(
    u_pos, n, f_g, f_pos, f_neg, kappa, od_g, od_r, f_g_ref, f_pos_ref, f_neg_ref
):
    """Three-component gDNA/RNA₊/RNA₋ strand loglik — :func:`strand_loglik` generalized to two RNA strands.

    Broadcasts ``(u_pos, n)`` of shape ``(nodes, 1)`` against the lattice ``(f_*)`` of shape
    ``(1, P)`` → ``(nodes, P)``. Mean ``N·p`` with ``p = ½·f_g + κ·f₊ + (1−κ)·f₋``.

    **Count-zero-information freeze** (`docs/calibration/archive/node_prior_design.md` §2): the mean stays LIVE in
    the solved composition ``(f_g, f_pos, f_neg)`` — the legitimate strand channel — but the variance is
    evaluated at the fixed REFERENCE composition ``(f_g_ref, f_pos_ref, f_neg_ref)`` (per-node scalars,
    broadcast). This keeps the heteroscedastic precision (the count still sets a composition-aware variance
    via the reference) while removing the ``f_g``-tilt of the normalizer, so the raw count can no longer
    manufacture a composition preference toward the variance-minimum when the mean degenerates (κ→½). The
    reference is a NEUTRAL structural default at init and the incoming belief in the sweep.
    """
    p = (
        0.5 * f_g + kappa * f_pos + (1.0 - kappa) * f_neg
    )  # LIVE mean channel (the composition solved)
    mean = n * p
    rscale = kappa * (1.0 - kappa)  # κ(1−κ): each RNA strand's μ(1−μ)
    # Variance at the REFERENCE composition (NOT the solved f_g) — the freeze.
    p_ref = 0.5 * f_g_ref + kappa * f_pos_ref + (1.0 - kappa) * f_neg_ref
    var = (
        n * p_ref * (1.0 - p_ref)
        + (n * f_g_ref) ** 2 * 0.25 * od_g
        + (n * f_pos_ref) ** 2 * rscale * od_r
        + (n * f_neg_ref) ** 2 * rscale * od_r
    )
    var = np.maximum(var, _EPS)
    return -0.5 * (u_pos - mean) ** 2 / var - 0.5 * np.log(var)
