"""Phase-2b — the odds-propagation grid sum-product sweep (CALIBRATION_PLAN §2/§4).

Belief propagation on a per-reference chain of region nodes. Each node's latent is the pie
``(f₊, f₋, f_g)`` on the 2-simplex triangular lattice (P points). Messages are ``(P,)`` log-likelihood
vectors; the edge transition is ``(P,P)``. The propagated quantity is the **per-strand RNA:gDNA log-odds**
``log(f_c/f_g)`` — the only enrichment-AND-content-invariant signal (CALIBRATION_PLAN §0b): the edge
penalises the *difference* in log-odds, gated to same-strand exon stretches (continuous RNA) and decoupled
(`Q=∞`) at exon↔intron / silent-strand transitions. gDNA is the residual (it falls out of the simplex,
never propagated). Exact in two sweeps on a chain (forward + backward); per-reference chunked.

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

from .signature import BIT_EXON_NEG, BIT_EXON_POS, TS_AMBIG, TS_NEG, TS_POS
from .simplex import _mixture_strand_loglik, _simplex_lattice
from .strand_deconv import NodeDeconv

__all__ = ["deconv_regions_sweep"]

_EPS = 1.0e-9


def _log_odds(f_pos, f_neg, f_g):
    """Per-lattice-point log-odds ``log(f_c/f_g)`` (clipped to bound the simplex edges)."""
    fp = np.clip(f_pos, _EPS, 1.0)
    fn = np.clip(f_neg, _EPS, 1.0)
    g = np.clip(f_g, _EPS, 1.0)
    return np.log(fp) - np.log(g), np.log(fn) - np.log(g)


def _local_loglik(u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, kappa, od_g, od_r,
                  gdna_prior_count, lattice, count_frac=None, count_trust_beta=0.0):
    """ψ_i over the lattice — strand mixture + sided spliced lower bound + (β-trusted) count + gDNA prior.

    The count is **local** evidence (rev-3: the *odds* is the propagated message, not the count), at the
    fixed trust ``β`` (the gDNA-magnitude fallback where strand/odds are weak). ``(n, P)``.
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
    if count_frac is not None and count_trust_beta > 0.0:
        psi = psi - 0.5 * count_trust_beta * (f_g[None, :] - count_frac[:, None]) ** 2
    psi = psi + gdna_prior_count * np.log(f_g[None, :] + _EPS)
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
    rna_strand_overdispersion=0.0, count_gdna_frac=None, count_trust_beta=0.0, gdna_prior_count=0.0,
    q_rna=0.25, n_grid=20,
):
    """Per-region gDNA fraction by the odds-propagation grid sum-product (see module docstring).

    ``count_gdna_frac`` + ``count_trust_beta``: the local (β-trusted) count fallback in ``ψ_i`` (None ⇒
    off, e.g. the toy tests). ``q_rna`` is the log-odds coupling variance (RNA relative variance) — a scalar
    placeholder from the Phase-2a fit; finite only on same-strand exon edges. ``n_grid`` = the lattice
    resolution K (P points). ``f_g`` is read as the posterior **median** (skew-safe).
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

    f_pos_g, f_neg_g, f_g_g = _simplex_lattice(int(n_grid))
    lattice = (f_pos_g, f_neg_g, f_g_g)
    lo_pos, lo_neg = _log_odds(f_pos_g, f_neg_g, f_g_g)
    cf = None if count_gdna_frac is None else np.clip(np.asarray(count_gdna_frac, np.float64), 0.0, 1.0)
    psi = _local_loglik(u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, kappa,
                        gdna_strand_overdispersion, rna_strand_overdispersion, gdna_prior_count, lattice,
                        count_frac=cf, count_trust_beta=count_trust_beta)

    exon_pos = (sig & BIT_EXON_POS) != 0
    exon_neg = (sig & BIT_EXON_NEG) != 0
    f_g = np.zeros(U.shape, dtype=np.float64)
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

    f_g = np.where(U > 0.0, np.clip(f_g, 0.0, 1.0), 0.0)
    return NodeDeconv(
        gdna_mass=f_g * mass_unspl, rna_mass=(1.0 - f_g) * mass_unspl + mass_spliced,
        gdna_frac=f_g, gdna_frac_var=np.zeros_like(f_g),
    )
