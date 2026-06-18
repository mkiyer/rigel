"""The per-node ψ solve over the 2-simplex pie ``(f₊, f₋, f_g)`` — the shared core of the BP sweep.

A node's latent is the pie ``(f₊, f₋, f_g)`` on the triangular lattice (P points). Its evidence ``ψ_i``
(:func:`_local_loglik`) is the sum of: the 3-component strand mixture (`simplex._mixture_strand_loglik`, the
only count→precision path); the sided spliced lower bound; the node-class prior (Jeffreys at single-strand
nodes, the global gDNA prior at AMBIG / intergenic); and the optional gDNA / per-strand RNA imputation priors
(the neighbour messages). ``f_g`` is the posterior **median** of ``ψ_i`` over the lattice, ``f_g_var`` its
posterior **variance** (the per-node confidence), ``f±`` the posterior **means**; a forbidden strand axis is
masked off, so an intergenic node keeps only the ``f_g=1`` vertex.

This module is the **per-node primitives** (``_local_loglik`` + the marginal helpers + ``_solve_nodes``);
:mod:`rigel.calibration.bp_solver` drives them over the region↔boundary chain (init + the directional sweep).
"""

from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

from .simplex import _mixture_strand_loglik, _simplex_lattice
from .strand_deconv import NodeDeconv

__all__ = ["_solve_nodes", "_local_loglik", "_fg_median", "_fg_var", "_axis_mean", "_simplex_lattice"]

_EPS = 1.0e-9
_PRIOR_EPS = 1.0e-3  # Jeffreys edge floor (the exact lattice vertex with 1e-9 over-rewards the prior)
_STRAND_PRIOR = 0.5  # Beta(½,½) strand reference prior (matches strand_deconv._STRAND_PRIOR / the fusion)


def _local_loglik(u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, kappa, od_g, od_r,
                  lattice, strand_obs=None, global_mu=None, global_tau=0.0,
                  gdna_imp_frac=None, gdna_imp_precision=None,
                  rna_imp_frac=None, rna_imp_precision=None):
    """ψ_i over the lattice — strand mixture + sided spliced lower bound + node-class prior + imputation.

    The Bayesian hierarchy (CALIBRATION_ARCHITECTURE.md §1), with a **node-class-dependent** prior. The
    count enters ONLY through the strand mixture's overdispersed Fisher information (§0 invariant); there is
    no count prior and no RNA-magnitude prior — a fragment count carries no intrinsic gDNA/RNA composition
    information, only the strand mixture's precision. A node is solved by its strand likelihood + the
    node-class prior + the neighbour imputation priors:

    * **single-strand (strand-observable) nodes** carry the **Beta(½,½) strand reference prior** (the
      Jeffreys term `(½−1)·(log f_g + log f_pos + log f_neg)`). It concentrates ``f_g`` at the
      likelihood-favoured value against the mixture's overdispersion-induced spread (without it, balanced
      pure-gDNA nodes under-call gDNA). It is the correct reference for a binomial proportion and is NOT the
      harmful vertex-push (the likelihood already favours that vertex at a single-strand node).
    * **AMBIG / intergenic (non-strand-observable) nodes** carry the **global gDNA prior** (the foundation)
      instead — a Gaussian pull of ``f_g`` toward the per-node baseline
      ``global_mu = clip(ρ_global·eff_len/mass, 0, 1)`` with **per-node precision ``global_tau``**
      (``1/σ²_global``, σ²_global = the robust MAD spread of the per-node densities in ``calibrate``; a
      scalar/array, ≥ the 1-pseudo-observation foundation). Here the U-shaped Jeffreys *would* push a
      flat-likelihood node to the ``f_g=1`` vertex (phantom gDNA); the informative global baseline replaces
      it, so in a pure-RNA library (``ρ_global ≈ 0`` ⇒ ``global_mu ≈ 0``) an unanchored AMBIG node settles
      at ``f_g ≈ 0``. ``(n, P)``.
    """
    f_pos, f_neg, f_g = lattice
    n = u_pos + u_neg
    psi = _mixture_strand_loglik(
        u_pos[:, None], n[:, None], f_g[None, :], f_pos[None, :], f_neg[None, :], kappa, od_g, od_r
    )
    inv_n = np.where(n > 0.0, 1.0 / np.maximum(n, _EPS), 0.0)[:, None]
    short_p = np.maximum(spliced_pos[:, None] * inv_n - f_pos[None, :], 0.0)
    short_n = np.maximum(spliced_neg[:, None] * inv_n - f_neg[None, :], 0.0)
    psi = psi - 0.5 * spliced_pos[:, None] * short_p**2 - 0.5 * spliced_neg[:, None] * short_n**2
    # node-class prior: Jeffreys reference at strand-observable nodes, global gDNA prior elsewhere.
    if strand_obs is not None:
        jeff = (_STRAND_PRIOR - 1.0) * (
            np.log(np.clip(f_g, _PRIOR_EPS, 1.0))
            + np.log(np.clip(f_pos, _PRIOR_EPS, 1.0))
            + np.log(np.clip(f_neg, _PRIOR_EPS, 1.0))
        )
        psi = psi + strand_obs[:, None].astype(np.float64) * jeff[None, :]
        if global_mu is not None:
            gt = np.asarray(global_tau, dtype=np.float64)
            gt = gt[:, None] if gt.ndim else gt  # per-node array → column; scalar → broadcast
            gpen = -0.5 * gt * (f_g[None, :] - global_mu[:, None]) ** 2
            psi = psi + (~strand_obs)[:, None].astype(np.float64) * gpen
    # gDNA IMPUTATION prior: a Gaussian pull of f_g toward the neighbour-imputed fraction (the identity
    # mean) at the precision 1/(σ²_bio + predictor Poisson noise) — CALIBRATION_ARCHITECTURE §1.2 +
    # imputation_variance_model.md. τ_imp=0 (no observable flank) ⇒ a no-op; a sharp strand likelihood
    # dominates this diffuse prior by their honest precisions (emergent deference, no tuned weight).
    if gdna_imp_frac is not None and gdna_imp_precision is not None:
        ti = np.asarray(gdna_imp_precision, dtype=np.float64)[:, None]
        psi = psi - 0.5 * ti * (f_g[None, :] - np.asarray(gdna_imp_frac, dtype=np.float64)[:, None]) ** 2
    # RNA IMPUTATION prior: per-strand Gaussian pull of f± toward the neighbour-imputed strand-s fraction
    # (each strand INDEPENDENT — the unspliced mass is partitioned by the fraction state, never shared). Since
    # f₊+f₋+f_g=1, pulling f± toward their imputed RNA sharpens f_g. τ=0 (no strand-s RNA neighbour) ⇒ no-op.
    # CALIBRATION_ARCHITECTURE §1.2 + the R↔B↔R chain (the nascent crosses exon↔intron via the boundary node).
    if rna_imp_frac is not None and rna_imp_precision is not None:
        for f_axis, mu_s, tau_s in (
            (f_pos, rna_imp_frac[0], rna_imp_precision[0]),
            (f_neg, rna_imp_frac[1], rna_imp_precision[1]),
        ):
            tr = np.asarray(tau_s, dtype=np.float64)[:, None]
            psi = psi - 0.5 * tr * (f_axis[None, :] - np.asarray(mu_s, dtype=np.float64)[:, None]) ** 2
    forbid = ((~allow_pos)[:, None] & (f_pos[None, :] > _EPS)) | (
        (~allow_neg)[:, None] & (f_neg[None, :] > _EPS)
    )
    return np.where(forbid, -np.inf, psi)


def _fg_median(belief, f_g_g):
    """Per-node posterior **median** of ``f_g`` over the lattice (the belief marginalised to the ``f_g``
    axis). The median — not the mean — avoids the strand-posterior overdispersion skew that biased the
    earlier grid-MAP low (the +8.7 pt regression)."""
    post = np.exp(belief - logsumexp(belief, axis=1, keepdims=True))  # (m,P)
    order = np.argsort(f_g_g)
    fgs = f_g_g[order]
    cw = np.cumsum(post[:, order], axis=1)
    idx = np.clip((cw < 0.5).sum(axis=1), 0, fgs.size - 1)  # first lattice f_g with CDF ≥ ½
    return fgs[idx]


def _fg_var(belief, f_g_g):
    """Per-node posterior VARIANCE of ``f_g`` over the lattice belief — the per-node confidence (the
    combined strand+count+global precision; small ⇒ confident). Feeds ``NodeDeconv.gdna_frac_var``."""
    post = np.exp(belief - logsumexp(belief, axis=1, keepdims=True))  # (m,P)
    mean = post @ f_g_g
    return np.maximum(post @ (f_g_g**2) - mean**2, 0.0)


def _axis_mean(belief, axis_g):
    """Per-node posterior MEAN of a lattice axis (f_pos or f_neg) over the belief. Feeds the per-strand RNA
    imputation (the current-state partition of the unspliced mass — the strands never share it)."""
    post = np.exp(belief - logsumexp(belief, axis=1, keepdims=True))  # (m,P)
    return post @ axis_g


def _solve_nodes(
    u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, strand_obs,
    mass_unspl, mass_spliced, *, kappa, od_g, od_r, n_grid,
    mu_global=None, global_tau=0.0, gdna_imp_frac=None, gdna_imp_precision=None,
    rna_imp_frac=None, rna_imp_precision=None,
) -> NodeDeconv:
    """The node-agnostic per-node simplex solve — solves an array of nodes (regions OR boundaries) from their
    per-node sufficient statistics, with no knowledge of node type. Each node's belief IS its local evidence
    ``ψ_i`` (``_local_loglik``: strand mixture + sided spliced floor + node-class prior + imputation pulls);
    the cross-node imputation enters ψ_i as a prior, computed upstream. ``f_g`` = posterior median over the
    lattice, ``f_g_var`` = posterior variance, ``f_pos``/``f_neg`` = posterior MEANS (the current-state
    partition for the per-strand RNA imputation); a node with no fragments (``u_pos+u_neg == 0``) reports
    ``f_g = f_pos = f_neg = 0``. The region/boundary wrappers build the arrays + the global baseline and call
    this core.
    """
    f_pos_g, f_neg_g, f_g_g = _simplex_lattice(int(n_grid))
    lattice = (f_pos_g, f_neg_g, f_g_g)
    psi = _local_loglik(u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, kappa,
                        od_g, od_r, lattice,
                        strand_obs=strand_obs, global_mu=mu_global, global_tau=global_tau,
                        gdna_imp_frac=gdna_imp_frac, gdna_imp_precision=gdna_imp_precision,
                        rna_imp_frac=rna_imp_frac, rna_imp_precision=rna_imp_precision)
    f_g = _fg_median(psi, f_g_g)
    f_g_var = _fg_var(psi, f_g_g)
    f_pos = _axis_mean(psi, f_pos_g)
    f_neg = _axis_mean(psi, f_neg_g)
    active = (u_pos + u_neg) > 0.0
    f_g = np.where(active, np.clip(f_g, 0.0, 1.0), 0.0)
    f_g_var = np.where(active, f_g_var, 0.0)
    f_pos = np.where(active, np.clip(f_pos, 0.0, 1.0), 0.0)
    f_neg = np.where(active, np.clip(f_neg, 0.0, 1.0), 0.0)
    return NodeDeconv(
        gdna_mass=f_g * mass_unspl, rna_mass=(1.0 - f_g) * mass_unspl + mass_spliced,
        gdna_frac=f_g, gdna_frac_var=f_g_var, rna_pos_frac=f_pos, rna_neg_frac=f_neg,
    )
