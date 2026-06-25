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

__all__ = ["_solve_nodes", "_local_loglik", "_binom_pseudo", "_fg_median", "_fg_var", "_axis_mean",
           "_simplex_lattice"]

_EPS = 1.0e-9
_PRIOR_EPS = 1.0e-3  # Jeffreys edge floor (the exact lattice vertex with 1e-9 over-rewards the prior)
_STRAND_PRIOR = 0.5  # Beta(½,½) strand reference prior (matches strand_deconv._STRAND_PRIOR / the fusion)


def _binom_pseudo(f, mu, count):
    """Count-space (Beta/binomial) pseudo-count log-prior on ONE axis: ``α·log f + (count−α)·log(1−f)``,
    ``α = count·μ`` — i.e. ``count`` pseudo-fragments split (this-component vs rest) at fraction ``μ``.

    This is the count-space replacement for the Gaussian ``−½·τ·(f−μ)²`` prior, with NO ``(M/E)²`` Jacobian:
    the curvature at the mode is ``count/(μ(1−μ))`` — the binomial Fisher information of ``count``
    pseudo-fragments — so the prior's weight is the COUNT (competing one-for-one with the strand likelihood's
    real fragments), never the destination mass². Two-sided (mode at ``μ``); the complement term
    ``(count−α)·log(1−f)`` pins the wall when ``μ→0`` (the part the rejected one-sided ``α·log f`` lacked).
    ``f`` is a lattice row ``(1,P)``; ``mu``/``count`` are node columns ``(m,1)`` → ``(m,P)``; ``f`` is clipped
    to ``[_PRIOR_EPS, 1−_PRIOR_EPS]`` for the logs (same edge floor as the Jeffreys term)."""
    alpha = count * mu
    fc = np.clip(f, _PRIOR_EPS, 1.0 - _PRIOR_EPS)
    return alpha * np.log(fc) + (count - alpha) * np.log1p(-fc)


def _local_loglik(u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, kappa, od_g, od_r,
                  lattice, strand_obs=None, global_logprior=None,
                  gdna_imp_mode=None, gdna_imp_prec=None,
                  rna_imp_mode=None, rna_imp_prec=None):
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
    * **ALL solvable nodes** (single-strand AND AMBIG — the ``strand_obs`` fork is dissolved) carry the
      **population gDNA prior** as a precomputed **count-space Negative-Binomial** log-prior on ``f_g``
      (``global_logprior``, an ``(n,P)`` array built in ``bp_solver.node_sweep``). It is a count-currency
      term with NO ``(M/E)²`` Jacobian: its pin/defer is governed by counts + length + the seed-fit gDNA
      between-node spread. In a zero-gDNA library the seed-based ``ρ_global → 0`` so it pins ``f_g → 0``;
      under capture the spread is large so it defers to the messages. Applying it to single-strand nodes too
      counteracts the Jeffreys vertex-push at κ≈½ (the Bug-C fix). ``(n, P)``.
    * **imputation messages** enter as per-component **density-Gaussians** ``−½·prec·(f_c − mode)²``
      (`precision_state_design.md`): ``mode`` = the neighbour's imputed fraction, ``prec`` = the honest
      precision (the source's own `Var_own` blended with the communication `σ²_bio` + sampling). An unsure
      source ⇒ low ``prec`` ⇒ weak pull; no edge log-wall. The count still enters the LIKELIHOOD only through
      the strand mixture's overdispersed Fisher information (§0 invariant).
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
    # node-class prior #1: Jeffreys Beta(½,½) strand reference at single-strand (strand-observable) nodes.
    if strand_obs is not None:
        jeff = (_STRAND_PRIOR - 1.0) * (
            np.log(np.clip(f_g, _PRIOR_EPS, 1.0))
            + np.log(np.clip(f_pos, _PRIOR_EPS, 1.0))
            + np.log(np.clip(f_neg, _PRIOR_EPS, 1.0))
        )
        psi = psi + strand_obs[:, None].astype(np.float64) * jeff[None, :]
    # node-class prior #2: the population gDNA prior as a COUNT-SPACE Negative-Binomial log-prior on f_g,
    # precomputed per node (no (M/E)² Jacobian) and applied to ALL solvable nodes (fork dissolved — §3.3).
    if global_logprior is not None:
        psi = psi + np.asarray(global_logprior, dtype=np.float64)
    # gDNA IMPUTATION message — a DENSITY-GAUSSIAN on f_g: ``−½·prec·(f_g − mode)²`` (precision_state_design.md
    # §5/§7, Phase-0 decision). ``mode`` = the imputed fraction (ρ_src·E_dst−spliced)/M_dst; ``prec`` =
    # (M_dst/E_dst)²/(Var_own,src + σ²_bio + pois) — the source's own uncertainty is blended into the variance,
    # so an unsure source ⇒ low prec ⇒ weak pull (NO log-wall, NO cap). prec=0 (no flank) ⇒ a no-op.
    if gdna_imp_mode is not None and gdna_imp_prec is not None:
        m = np.asarray(gdna_imp_mode, dtype=np.float64)[:, None]
        p = np.asarray(gdna_imp_prec, dtype=np.float64)[:, None]
        psi = psi - 0.5 * p * (f_g[None, :] - m) ** 2
    # per-strand RNA IMPUTATION messages — the same density-Gaussian on each strand axis (each INDEPENDENT;
    # since f₊+f₋+f_g=1, pulling f± toward their imputed RNA sharpens f_g). prec=0 (no strand-s flank) ⇒ no-op.
    if rna_imp_mode is not None and rna_imp_prec is not None:
        for f_axis, m_s, p_s in (
            (f_pos, rna_imp_mode[0], rna_imp_prec[0]),
            (f_neg, rna_imp_mode[1], rna_imp_prec[1]),
        ):
            psi = psi - 0.5 * np.asarray(p_s, dtype=np.float64)[:, None] * (
                f_axis[None, :] - np.asarray(m_s, dtype=np.float64)[:, None]
            ) ** 2
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
    """Per-node posterior VARIANCE of a lattice axis (``f_g``, ``f_pos`` or ``f_neg`` — pass the matching grid)
    over the belief: ``E[f²] − E[f]²``. This is the **precision state** (`precision_state_design.md`): small ⇒
    confident, large ⇒ unresolved (e.g. a balanced AMBIG node's flat ``f_g``). Moment-matched (mean/variance),
    distinct from the robust ``_fg_median`` readout (the two operators, plan P2). Feeds ``NodeBelief.var_*`` and
    the honest message send ``Var_own = (M/E)²·Var(f_c)``."""
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
    global_logprior=None, gdna_imp_mode=None, gdna_imp_prec=None,
    rna_imp_mode=None, rna_imp_prec=None,
) -> NodeDeconv:
    """The node-agnostic per-node simplex solve — solves an array of nodes (regions OR boundaries) from their
    per-node sufficient statistics, with no knowledge of node type. Each node's belief IS its local evidence
    ``ψ_i`` (``_local_loglik``: strand mixture + sided spliced floor + node-class prior + imputation pulls);
    the cross-node imputation enters ψ_i as a prior, computed upstream. ``f_g`` = posterior median over the
    lattice, ``f_pos``/``f_neg`` = posterior MEANS (the current-state partition for the per-strand RNA
    imputation), and ``*_frac_var`` = the per-component posterior VARIANCES (the precision state, moment-matched
    — `precision_state_design.md`); a node with no fragments (``u_pos+u_neg == 0``) reports all of them 0. The
    region/boundary wrappers build the arrays + the global baseline and call this core.
    """
    f_pos_g, f_neg_g, f_g_g = _simplex_lattice(int(n_grid))
    lattice = (f_pos_g, f_neg_g, f_g_g)
    psi = _local_loglik(u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, kappa,
                        od_g, od_r, lattice,
                        strand_obs=strand_obs, global_logprior=global_logprior,
                        gdna_imp_mode=gdna_imp_mode, gdna_imp_prec=gdna_imp_prec,
                        rna_imp_mode=rna_imp_mode, rna_imp_prec=rna_imp_prec)
    f_g = _fg_median(psi, f_g_g)
    f_pos = _axis_mean(psi, f_pos_g)
    f_neg = _axis_mean(psi, f_neg_g)
    var_g = _fg_var(psi, f_g_g)
    var_pos = _fg_var(psi, f_pos_g)
    var_neg = _fg_var(psi, f_neg_g)
    active = (u_pos + u_neg) > 0.0
    f_g = np.where(active, np.clip(f_g, 0.0, 1.0), 0.0)
    f_pos = np.where(active, np.clip(f_pos, 0.0, 1.0), 0.0)
    f_neg = np.where(active, np.clip(f_neg, 0.0, 1.0), 0.0)
    var_g = np.where(active, var_g, 0.0)
    var_pos = np.where(active, var_pos, 0.0)
    var_neg = np.where(active, var_neg, 0.0)
    return NodeDeconv(
        gdna_mass=f_g * mass_unspl, rna_mass=(1.0 - f_g) * mass_unspl + mass_spliced,
        gdna_frac=f_g, rna_pos_frac=f_pos, rna_neg_frac=f_neg,
        gdna_frac_var=var_g, rna_pos_frac_var=var_pos, rna_neg_frac_var=var_neg,
    )
