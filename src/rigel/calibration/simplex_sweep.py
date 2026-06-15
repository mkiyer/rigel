"""Phase-2b — the odds-propagation grid sum-product sweep (CALIBRATION_PLAN §2/§4).

Belief propagation on a per-reference chain of region nodes. Each node's latent is the pie
``(f₊, f₋, f_g)`` on the 2-simplex triangular lattice (P points). Messages are ``(P,)`` log-likelihood
vectors; the edge transition is ``(P,P)``. The propagated quantity is the **per-strand RNA:gDNA log-odds**
``log(f_c/f_g)`` — the only enrichment-AND-content-invariant signal (CALIBRATION_PLAN §0b): the edge
penalises the *difference* in log-odds, gated to same-strand exon stretches (continuous RNA) and decoupled
(`Q=∞`) at exon↔intron / silent-strand transitions. gDNA is the residual of the simplex (never propagated
*directly* — but it couples implicitly, since the propagated log-odds `log(f_c/f_g)` share `f_g` in the
denominator). Exact in two sweeps on a chain (forward + backward); per-reference chunked.

`ψ_i` (local evidence) = the 3-component strand mixture (`simplex._mixture_strand_loglik`) + the sided
spliced lower bound + a weak gDNA prior — i.e. `solve_node` minus the count term (the count is not local
evidence; the cross-node signal is the odds). Seeds emerge from `ψ_i` (intergenic → only the `f_g=1`
vertex survives the strand mixture + prior).

This is the grid sum-product that superseded the (now-removed) scalar-RTS scaffold. Boundary sides keep
`strand_deconv.deconv_sides` (the flux transport is unchanged, post-sweep). `Q_RNA` (the log-odds coupling
variance = the RNA relative variance) is currently a scalar from the Phase-2a fit median — a flagged
placeholder; the per-edge `var~mean` lookup is the refinement.
"""

from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

from .signature import BIT_EXON_NEG, BIT_EXON_POS, TS_AMBIG, TS_NEG, TS_NONE, TS_POS
from .simplex import _mixture_strand_loglik, _simplex_lattice
from .strand_deconv import NodeDeconv

__all__ = ["deconv_regions_sweep"]

_EPS = 1.0e-9
_PRIOR_EPS = 1.0e-3  # Jeffreys edge floor (the exact lattice vertex with 1e-9 over-rewards the prior)
_STRAND_PRIOR = 0.5  # Beta(½,½) strand reference prior (matches strand_deconv._STRAND_PRIOR / the fusion)


def _log_odds(f_pos, f_neg, f_g):
    """Per-lattice-point log-odds ``log(f_c/f_g)`` (clipped to bound the simplex edges)."""
    fp = np.clip(f_pos, _EPS, 1.0)
    fn = np.clip(f_neg, _EPS, 1.0)
    g = np.clip(f_g, _EPS, 1.0)
    return np.log(fp) - np.log(g), np.log(fn) - np.log(g)


def _local_loglik(u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, kappa, od_g, od_r,
                  lattice, count_frac=None, count_trust_beta=0.0, count_precision=None, strand_obs=None,
                  global_mu=None, global_tau=0.0):
    """ψ_i over the lattice — strand mixture + sided spliced lower bound + count prior + node-class prior.

    The Bayesian hierarchy (per_node_deconv_hierarchy_design.md), with a **node-class-dependent** prior:

    * **single-strand (strand-observable) nodes** carry the **Beta(½,½) strand reference prior** (the
      Jeffreys term `(½−1)·(log f_g + log f_pos + log f_neg)`) — exactly what the fusion's
      ``strand_posterior_gdna_frac`` uses, ungated. It concentrates ``f_g`` at the likelihood-favoured
      value against the mixture's overdispersion-induced spread (without it, balanced pure-gDNA nodes
      under-call gDNA). It is the correct reference for a binomial proportion and is NOT the harmful
      vertex-push (the likelihood already favours that vertex at a single-strand node).
    * **AMBIG / intergenic (non-strand-observable) nodes** carry the **global gDNA prior** instead — a
      Gaussian pull of ``f_g`` toward the per-node baseline ``global_mu = clip(ρ_global·eff_len/mass, 0, 1)``
      with **per-node precision ``global_tau``** (``1/σ²_global``, where ``σ²_global`` is the robust MAD
      spread of the per-node densities in ``calibrate``; a scalar/array, ≥ the 1-pseudo-observation
      foundation). Here the U-shaped Jeffreys *would* push a
      flat-likelihood node to the ``f_g=1`` vertex (phantom gDNA); the informative global baseline replaces
      it, so in a pure-RNA library (``ρ_global ≈ 0`` ⇒ tight ``global_tau``) an unanchored AMBIG node
      settles at ``f_g ≈ 0``.

    The count term is a Gaussian prior on ``f_g`` toward ``count_frac`` with per-node precision
    ``count_precision`` (= ``1/σ²_frac`` from the var~mean, capped at the node's own-count Poisson ceiling;
    falls back to the scalar β for the toy tests). On the simplex lattice each quadratic is a truncated
    Gaussian (renormalized over the in-simplex points). ``(n, P)``.
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
    if count_frac is not None:
        # count prior: Gaussian pull of f_g toward count_frac. Per-node precision τ_count = count_precision
        # (= 1/v_rel, the degeneracy-free count reliability) when supplied; else the scalar β (toy tests).
        if count_precision is not None:
            tau_c = np.asarray(count_precision, dtype=np.float64)[:, None]
            psi = psi - 0.5 * tau_c * (f_g[None, :] - count_frac[:, None]) ** 2
        elif count_trust_beta > 0.0:
            psi = psi - 0.5 * count_trust_beta * (f_g[None, :] - count_frac[:, None]) ** 2
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


def _edge_logphi(lo_pos, lo_neg, q_pos, q_neg):
    """``(P,P)`` log-coupling: ``−½[(Δlo₊)²/Q₊ + (Δlo₋)²/Q₋]``; ``Q=∞`` ⇒ that strand decoupled (0)."""
    p = lo_pos.size
    out = np.zeros((p, p), dtype=np.float64)
    if np.isfinite(q_pos) and q_pos > 0.0:
        d = lo_pos[None, :] - lo_pos[:, None]
        out = out - 0.5 * d * d / q_pos
    if np.isfinite(q_neg) and q_neg > 0.0:
        d = lo_neg[None, :] - lo_neg[:, None]
        out = out - 0.5 * d * d / q_neg
    return out


def _sweep_chain(psi, q_pos_edges, q_neg_edges, lo_pos, lo_neg):
    """Forward + backward grid sum-product over one chain. ``psi`` ``(m,P)``; returns belief ``(m,P)``."""
    m = psi.shape[0]
    if m == 1:
        return psi
    mats = [_edge_logphi(lo_pos, lo_neg, q_pos_edges[k], q_neg_edges[k]) for k in range(m - 1)]
    alpha = np.zeros_like(psi)
    for i in range(1, m):
        prev = psi[i - 1] + alpha[i - 1]  # (P,)
        alpha[i] = logsumexp(prev[:, None] + mats[i - 1], axis=0)  # over prev state → (P,)
    beta = np.zeros_like(psi)
    for i in range(m - 2, -1, -1):
        nxt = psi[i + 1] + beta[i + 1]
        beta[i] = logsumexp(mats[i] + nxt[None, :], axis=1)  # over next state → (P,)
    return psi + alpha + beta


def deconv_regions_sweep(
    substrate, region_arrays, *, rna_sense_frac, gdna_strand_overdispersion=0.0,
    rna_strand_overdispersion=0.0, count_gdna_frac=None, count_trust_beta=0.0, count_precision=None,
    q_rna=0.25, n_grid=20, rho_global=0.0, region_eff_len=None, info_scale=10.0, global_tau=None,
):
    """Per-region gDNA fraction by the odds-propagation grid sum-product (see module docstring).

    ``count_gdna_frac`` + ``count_precision`` (per-node ``τ_count``, from the var~mean; ``count_trust_beta``
    is the scalar fallback for the toy tests): the count prior toward the count fraction in ``ψ_i``, gated by
    ``(1−w)`` so it yields to the strand. ``rho_global`` + ``region_eff_len`` + ``global_tau``: the
    **global gDNA prior** (foundation) — baseline fraction ``μ_global = clip(ρ_global·eff_len/mass, 0, 1)``
    with per-node precision ``global_tau`` (``1/σ²_global`` from ``calibrate``; default 1-pseudo-observation;
    toy tests omit ``region_eff_len`` ⇒ no global prior). ``q_rna`` = log-odds coupling variance; ``n_grid``
    = lattice K. Returns the posterior-**median** ``f_g`` and its posterior **variance** (the per-node
    confidence). See per_node_deconv_hierarchy_design.md / CALIBRATION_PLAN_v2.md.
    """
    ts = np.asarray(region_arrays.strand_class)
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    ref_id = np.asarray(region_arrays.ref_id)
    kappa = float(rna_sense_frac)
    c = substrate.contained
    u_pos = c.n_unspliced_pos.astype(np.float64)
    u_neg = c.n_unspliced_neg.astype(np.float64)
    U = u_pos + u_neg
    mass_unspl = np.asarray(c.mass_unspliced, dtype=np.float64)
    mass_spliced = np.asarray(c.mass_spliced, dtype=np.float64)
    spl_sense = c.n_spliced_sense.astype(np.float64)
    # sided spliced (oriented): single-strand regions only (AMBIG relies on the propagated odds — deferred).
    spliced_pos = np.where(ts == TS_POS, spl_sense, 0.0)
    spliced_neg = np.where(ts == TS_NEG, spl_sense, 0.0)
    allow_pos = (ts == TS_POS) | (ts == TS_AMBIG)
    allow_neg = (ts == TS_NEG) | (ts == TS_AMBIG)

    # GLOBAL gDNA prior (foundation): per-node baseline fraction μ_global = clip(ρ_global·eff_len/mass, 0, 1)
    # — the count module's own count_gdna_frac form with the global pooled density substituted for the local
    # one. Weak (global_tau ≈ 1 pseudo-observation) so it only governs nodes the strand AND count leave
    # undetermined; ρ_global ≈ 0 in a pure-RNA library ⇒ μ_global ≈ 0 ⇒ unanchored nodes settle at f_g ≈ 0
    # (no phantom gDNA). Active only when region_eff_len is supplied (toy tests omit it ⇒ no global prior).
    strand_obs = (ts == TS_POS) | (ts == TS_NEG)
    # Count is a FALLBACK: down-weight its precision by (1−w), w = I/(I+I₀), I = N·(2κ−1)², so it YIELDS to
    # the strand+propagation where they are informative (single-strand, anchored AMBIG via the tilt) and
    # GOVERNS where they are silent (κ=½ / unanchored). w=0 at TS_NONE (intergenic — no transcript strand to
    # defer to ⇒ count governs). Without this the per-node count precision (~N) over-rides the propagation
    # at AMBIG and the complex-locus win is lost.
    info = U * (2.0 * kappa - 1.0) ** 2
    w_strand = np.where(ts == TS_NONE, 0.0, info / (info + float(info_scale)))
    eff_count_precision = (
        None if count_precision is None
        else (1.0 - w_strand) * np.asarray(count_precision, dtype=np.float64)
    )
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

    f_pos_g, f_neg_g, f_g_g = _simplex_lattice(int(n_grid))
    lattice = (f_pos_g, f_neg_g, f_g_g)
    lo_pos, lo_neg = _log_odds(f_pos_g, f_neg_g, f_g_g)
    cf = None if count_gdna_frac is None else np.clip(np.asarray(count_gdna_frac, np.float64), 0.0, 1.0)
    psi = _local_loglik(u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, kappa,
                        gdna_strand_overdispersion, rna_strand_overdispersion, lattice,
                        count_frac=cf, count_trust_beta=count_trust_beta,
                        count_precision=eff_count_precision,
                        strand_obs=strand_obs, global_mu=mu_global, global_tau=gtau)

    exon_pos = (sig & BIT_EXON_POS) != 0
    exon_neg = (sig & BIT_EXON_NEG) != 0
    f_g = np.zeros(U.shape, dtype=np.float64)
    f_g_var = np.zeros(U.shape, dtype=np.float64)
    # per-reference chains (genomic order); Q finite only on same-strand-exon edges, else ∞ (decouple).
    for ref in np.unique(ref_id):
        idx = np.flatnonzero(ref_id == ref)  # ascending = genomic order
        if idx.size == 0:
            continue
        qp = np.full(idx.size - 1, np.inf)
        qn = np.full(idx.size - 1, np.inf)
        for k in range(idx.size - 1):
            a, b = idx[k], idx[k + 1]
            if exon_pos[a] and exon_pos[b]:
                qp[k] = q_rna
            if exon_neg[a] and exon_neg[b]:
                qn[k] = q_rna
        belief = _sweep_chain(psi[idx], qp, qn, lo_pos, lo_neg)
        f_g[idx] = _fg_median(belief, f_g_g)
        f_g_var[idx] = _fg_var(belief, f_g_g)

    active = U > 0.0
    f_g = np.where(active, np.clip(f_g, 0.0, 1.0), 0.0)
    f_g_var = np.where(active, f_g_var, 0.0)
    return NodeDeconv(
        gdna_mass=f_g * mass_unspl, rna_mass=(1.0 - f_g) * mass_unspl + mass_spliced,
        gdna_frac=f_g, gdna_frac_var=f_g_var,
    )
