"""The per-node grid solve over the 2-simplex pie ``(f₊, f₋, f_g)`` (CALIBRATION_PLAN §2/§4).

Each region node's latent is the pie ``(f₊, f₋, f_g)`` on the triangular lattice (P points). The per-node
evidence ``ψ_i`` (`_local_loglik`) = the 3-component strand mixture (`simplex._mixture_strand_loglik`) + the
sided spliced lower bound + the node-class prior (Jeffreys at single-strand nodes, the global gDNA prior at
AMBIG / intergenic). The node's ``f_g`` is the posterior **median** of ``ψ_i`` over the lattice; its posterior
**variance** is the per-node confidence. Intergenic → only the ``f_g=1`` vertex survives (the forbid mask).

Cross-node **imputation** — the odds-propagation that previously resolved AMBIG nodes from their same-strand
exon neighbours (`log(f_c/f_g)` coupling on a per-reference chain) — was removed in the **Step-1 precision
rebuild** (CALIBRATION_ARCHITECTURE §6.5: the `q_rna` magic edge-coupling). Its principled,
reliability-weighted replacement lands in **Step 2**; until then every node is solved by its own local
evidence + the global foundation. Boundary sides keep `strand_deconv.deconv_sides` (the flux transport,
post-solve).
"""

from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

from .signature import TS_AMBIG, TS_NEG, TS_POS
from .simplex import _mixture_strand_loglik, _simplex_lattice
from .strand_deconv import NodeDeconv

__all__ = ["deconv_regions_sweep"]

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
    no count prior and no RNA-magnitude prior (both were count→fraction-precision violations, removed in the
    Step-1 precision rebuild — §6.1). A node is solved by its strand likelihood + the node-class prior:

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


def deconv_regions_sweep(
    substrate, region_arrays, *, rna_sense_frac, gdna_strand_overdispersion=0.0,
    rna_strand_overdispersion=0.0, n_grid=20, rho_global=0.0, region_eff_len=None, global_tau=None,
    gdna_imp_frac=None, gdna_imp_precision=None,
    rna_imp_frac_pos=None, rna_imp_frac_neg=None, rna_imp_precision_pos=None, rna_imp_precision_neg=None,
):
    """Per-region gDNA fraction by the per-node grid solve (see module docstring).

    The node combines its **strand likelihood** (the only count→precision path — the overdispersed BB Fisher
    info) with the **node-class prior**. ``rho_global`` + ``region_eff_len`` + ``global_tau``: the **global
    gDNA prior** (foundation) — baseline fraction ``μ_global = clip(ρ_global·eff_len/mass, 0, 1)`` with
    per-node precision ``global_tau`` (``1/σ²_global`` from ``calibrate``; default 1-pseudo-observation; toy
    tests omit ``region_eff_len`` ⇒ no global prior). ``n_grid`` = lattice K. Returns the posterior-**median**
    ``f_g`` and its posterior **variance** (the per-node confidence). Cross-node imputation is deferred to
    Step 2 (see module docstring). See CALIBRATION_ARCHITECTURE.md §1/§4.
    """
    ts = np.asarray(region_arrays.strand_class)
    c = substrate.contained
    u_pos = c.n_unspliced_pos.astype(np.float64)
    u_neg = c.n_unspliced_neg.astype(np.float64)
    mass_unspl = np.asarray(c.mass_unspliced, dtype=np.float64)
    mass_spliced = np.asarray(c.mass_spliced, dtype=np.float64)
    spl_sense = c.n_spliced_sense.astype(np.float64)
    # sided spliced (oriented): single-strand regions only (AMBIG relies on the propagated odds — deferred).
    spliced_pos = np.where(ts == TS_POS, spl_sense, 0.0)
    spliced_neg = np.where(ts == TS_NEG, spl_sense, 0.0)
    allow_pos = (ts == TS_POS) | (ts == TS_AMBIG)
    allow_neg = (ts == TS_NEG) | (ts == TS_AMBIG)
    strand_obs = (ts == TS_POS) | (ts == TS_NEG)

    # GLOBAL gDNA prior (foundation): per-node baseline fraction μ_global = clip(ρ_global·eff_len/mass, 0, 1).
    # Weak (global_tau ≈ 1 pseudo-observation) so it only governs nodes the strand leaves undetermined
    # (κ=½ / AMBIG / thin); ρ_global ≈ 0 in a pure-RNA library ⇒ μ_global ≈ 0 ⇒ unanchored nodes settle at
    # f_g ≈ 0 (no phantom gDNA). Active only when region_eff_len is supplied (toy tests omit it ⇒ no prior).
    if region_eff_len is not None:
        eff = np.asarray(region_eff_len, dtype=np.float64)
        mu_global = np.clip(
            np.where(mass_unspl > 0.0, float(rho_global) * eff / np.maximum(mass_unspl, _EPS), 0.0),
            0.0, 1.0,
        )
        # per-node global precision from calibrate (1/σ²_global, σ²_global = the robust MAD spread of the
        # per-node densities); default to the 1-pseudo-observation foundation when not supplied (toy tests).
        gtau = 1.0 if global_tau is None else np.asarray(global_tau, dtype=np.float64)
    else:
        mu_global = None
        gtau = 0.0

    # Per-strand RNA imputation pulls (CALIBRATION_ARCHITECTURE §1.2; the R↔B↔R chain): (μ_pos, μ_neg) +
    # (τ_pos, τ_neg). None ⇒ no RNA prior (a no-op in _local_loglik).
    rna_imp_frac = None if rna_imp_frac_pos is None else (rna_imp_frac_pos, rna_imp_frac_neg)
    rna_imp_precision = (
        None if rna_imp_precision_pos is None else (rna_imp_precision_pos, rna_imp_precision_neg)
    )
    return _solve_nodes(
        u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, strand_obs,
        mass_unspl, mass_spliced,
        kappa=float(rna_sense_frac), od_g=gdna_strand_overdispersion, od_r=rna_strand_overdispersion,
        n_grid=n_grid, mu_global=mu_global, global_tau=gtau,
        gdna_imp_frac=gdna_imp_frac, gdna_imp_precision=gdna_imp_precision,
        rna_imp_frac=rna_imp_frac, rna_imp_precision=rna_imp_precision,
    )
