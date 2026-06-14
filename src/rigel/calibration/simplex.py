"""The simplex per-node solver — partition a node's unspliced mass into ``(f_rna₊, f_rna₋, f_g)``.

Spec: ``docs/calibration/propagation_simplex_plan.md`` (rev-2). Calibration models **only RNA vs gDNA**
(the per-locus EM separates nascent from mature downstream), so every node's unspliced mass is one pie
cut into three slices ``(f_rna₊, f_rna₋, f_g)`` summing to 1 — **no nascent slice, no mature slice**.
Evidence *pushes* the point on the 2-simplex; the answer is always a normalized partition, so
over-subtraction is structurally impossible and an under-constrained slice stays wide ("unknown").

This module is the **per-node solve only** (plan §3/§7). Its ``_simplex_lattice`` and
``_mixture_strand_loglik`` are the reusable primitives the ``simplex_sweep`` grid sum-product builds on;
``solve_node`` (the grid-MAP combine) is retained as a primitive but is not the production combine.

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

from dataclasses import dataclass

import numpy as np

from .run_fill import same_ref_left_right
from .signature import TS_AMBIG, TS_NEG, TS_NONE, TS_POS

__all__ = ["SimplexSolution", "solve_node", "SignatureInit", "init_from_signature", "region_adjacency"]

_EPS = 1.0e-9


@dataclass(frozen=True, slots=True)
class SimplexSolution:
    """Per-node simplex partition (all float64[n]); ``f_rna_pos + f_rna_neg + f_g = 1`` where active."""

    f_rna_pos: np.ndarray
    f_rna_neg: np.ndarray
    f_g: np.ndarray
    f_g_var: np.ndarray  # posterior variance of f_g (the node's gDNA-fraction uncertainty / precision⁻¹)


@dataclass(frozen=True, slots=True)
class SignatureInit:
    """Per-node initialization from the signature (plan §2). All length-``R``.

    The signature is a *fact* about which slices can be nonzero: an absent transcript strand means that
    RNA slice is exactly 0 (certain), an active strand means that slice is unknown (to be earned by
    evidence), and gDNA is always admissible. Intergenic (``TS_NONE``) nodes carry no transcript, so
    they are fully specified seeds (``f_g = 1``).
    """

    allow_pos: np.ndarray  # bool — the +RNA slice may be nonzero (signature has a + transcript)
    allow_neg: np.ndarray  # bool — the −RNA slice may be nonzero
    f_g_fixed: np.ndarray  # float — pinned f_g (1.0 for intergenic seeds; NaN where not pinned)
    intergenic_seed: np.ndarray  # bool — fully specified at init (TS_NONE: f_g=1, no work)
    strand_seedable: np.ndarray  # bool — a single defined strand (POS/NEG): own strand can seed it


def init_from_signature(strand_class: np.ndarray) -> SignatureInit:
    """Map each node's transcript-strand class to its simplex init (plan §2).

    ``TS_NONE`` → ``(0,0,1)`` fully specified (intergenic seed); ``TS_POS`` → ``(nan,0,prior)``;
    ``TS_NEG`` → ``(0,nan,prior)``; ``TS_AMBIG`` → ``(nan,nan,prior)``. POS/NEG are *strand-seedable*
    (their own strand deconvolution can resolve them ab initio if the data is informative); AMBIG and
    NONE are not strand-seedable (NONE needs none; AMBIG is filled by propagation/count).
    """
    ts = np.asarray(strand_class)
    allow_pos = (ts == TS_POS) | (ts == TS_AMBIG)
    allow_neg = (ts == TS_NEG) | (ts == TS_AMBIG)
    intergenic_seed = ts == TS_NONE
    f_g_fixed = np.where(intergenic_seed, 1.0, np.nan)
    strand_seedable = (ts == TS_POS) | (ts == TS_NEG)
    return SignatureInit(
        allow_pos=allow_pos, allow_neg=allow_neg, f_g_fixed=f_g_fixed,
        intergenic_seed=intergenic_seed, strand_seedable=strand_seedable,
    )


def region_adjacency(ref_id: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Genomic-order ``(prev, next)`` region indices within each reference (−1 at a reference edge).

    The calibration region partition is sorted in genomic order per reference, so region ``i``'s chain
    neighbours are ``i−1`` / ``i+1`` when they share ``i``'s reference. This is the adjacency of the
    **chain** the propagation runs over (a locus is a linear run of regions/boundaries → a tree → exact
    two-sweep message passing; see ``docs/calibration/propagation_message_passing.md``).
    """
    ref = np.asarray(ref_id)
    r = ref.shape[0]
    left_same, right_same = same_ref_left_right(ref)
    prev = np.full(r, -1, dtype=np.int64)
    nxt = np.full(r, -1, dtype=np.int64)
    idx = np.arange(r)
    prev[left_same] = idx[left_same] - 1
    nxt[right_same] = idx[right_same] + 1
    return prev, nxt


def _simplex_lattice(n_grid: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The ``(f₊, f₋, f_g)`` lattice on the 2-simplex: all ``(i, j)`` with ``i+j ≤ n_grid``, ``/n_grid``.

    Returns three length-``P`` arrays (``P = (n_grid+1)(n_grid+2)/2``). The barycentric grid is uniform
    over the triangle, so the softmax posterior mean over it is an unbiased simplex integral.
    """
    g = int(n_grid)
    ii, jj = np.meshgrid(np.arange(g + 1), np.arange(g + 1), indexing="ij")
    keep = (ii + jj) <= g
    i = ii[keep].astype(np.float64)
    j = jj[keep].astype(np.float64)
    f_pos = i / g
    f_neg = j / g
    return f_pos, f_neg, 1.0 - f_pos - f_neg


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


def solve_node(
    u_pos,
    u_neg,
    *,
    kappa,
    spliced_pos=0.0,
    spliced_neg=0.0,
    count_gdna_frac=None,
    count_precision=0.0,
    allow_pos=True,
    allow_neg=True,
    strand_od_gdna=0.0,
    strand_od_rna=0.0,
    gdna_prior_count=0.5,
    rna_dirichlet_alpha=1.0,
    n_grid=60,
) -> SimplexSolution:
    """Solve each node's ``(f_rna₊, f_rna₋, f_g)`` pie by a precision-weighted MAP on the 2-simplex.

    All inputs broadcast to length ``n``. ``u_pos/u_neg`` are the node's plus/minus unspliced counts;
    ``spliced_pos/spliced_neg`` the sided spliced flux (0 for boundaries / single-exon); ``count_gdna_frac``
    + ``count_precision`` the propagated count evidence (``count_precision=0`` ⇒ no count pull);
    ``allow_pos/allow_neg`` the signature's active-strand mask (a disallowed slice is pinned to 0).
    Returns the posterior-mean partition + ``f_g`` posterior variance.

    Single-strand collapse (``allow_neg=False`` ⇒ ``f₋≡0``) reduces to today's strand posterior — the
    no-regression guard. Mass is conserved by construction (the lattice lies on the simplex).
    """
    u_pos = np.atleast_1d(np.asarray(u_pos, dtype=np.float64))
    u_neg = np.atleast_1d(np.asarray(u_neg, dtype=np.float64))
    n = u_pos + u_neg
    shape = np.broadcast(u_pos, u_neg).shape
    sp = np.broadcast_to(np.asarray(spliced_pos, dtype=np.float64), shape).astype(np.float64)
    sn = np.broadcast_to(np.asarray(spliced_neg, dtype=np.float64), shape).astype(np.float64)
    allow_p = np.broadcast_to(np.asarray(allow_pos, dtype=bool), shape)
    allow_n = np.broadcast_to(np.asarray(allow_neg, dtype=bool), shape)
    cg = (
        None if count_gdna_frac is None
        else np.broadcast_to(np.asarray(count_gdna_frac, dtype=np.float64), shape).astype(np.float64)
    )
    cprec = np.broadcast_to(np.asarray(count_precision, dtype=np.float64), shape).astype(np.float64)
    kappa = float(kappa)

    f_pos_g, f_neg_g, f_g_g = _simplex_lattice(n_grid)  # (P,)
    # Total mass the spliced lower bound divides into a fraction (skip empty nodes).
    n_col = n[:, None]
    inv_n = np.where(n_col > 0.0, 1.0 / np.maximum(n_col, _EPS), 0.0)

    logpost = _mixture_strand_loglik(
        u_pos[:, None], n_col, f_g_g[None, :], f_pos_g[None, :], f_neg_g[None, :],
        kappa, strand_od_gdna, strand_od_rna,
    )  # (n, P)

    # Count term: Gaussian pull of f_g toward count_gdna_frac, precision count_precision.
    if cg is not None:
        logpost = logpost - 0.5 * cprec[:, None] * (f_g_g[None, :] - cg[:, None]) ** 2

    # Spliced RNA lower bound (sided, soft one-sided): f_s ≥ spliced_s/N, precision = the spliced flux.
    lb_pos = sp * inv_n[:, 0]
    lb_neg = sn * inv_n[:, 0]
    short_pos = np.maximum(lb_pos[:, None] - f_pos_g[None, :], 0.0)
    short_neg = np.maximum(lb_neg[:, None] - f_neg_g[None, :], 0.0)
    logpost = logpost - 0.5 * sp[:, None] * short_pos**2 - 0.5 * sn[:, None] * short_neg**2

    # gDNA Dirichlet pseudo-count (RNA slices flat / "earned"); keeps f_g off zero with no evidence.
    logprior = (
        gdna_prior_count * np.log(f_g_g + _EPS)
        + (rna_dirichlet_alpha - 1.0) * (np.log(f_pos_g + _EPS) + np.log(f_neg_g + _EPS))
    )
    logpost = logpost + logprior[None, :]

    # Signature mask: pin a disallowed RNA strand to 0 (forbid lattice points using it).
    forbid = ((~allow_p)[:, None] & (f_pos_g[None, :] > _EPS)) | (
        (~allow_n)[:, None] & (f_neg_g[None, :] > _EPS)
    )
    logpost = np.where(forbid, -np.inf, logpost)

    post = np.exp(logpost - np.nanmax(np.where(np.isfinite(logpost), logpost, -np.inf), axis=1,
                                      keepdims=True))
    post = np.where(np.isfinite(logpost), post, 0.0)
    norm = post.sum(axis=1, keepdims=True)
    post = post / np.where(norm > 0.0, norm, 1.0)

    f_pos = post @ f_pos_g
    f_neg = post @ f_neg_g
    f_g = post @ f_g_g
    f_g_var = post @ (f_g_g**2) - f_g**2

    active = n > 0.0
    f_g = np.where(active, np.clip(f_g, 0.0, 1.0), 0.0)
    f_pos = np.where(active, np.clip(f_pos, 0.0, 1.0), 0.0)
    f_neg = np.where(active, np.clip(f_neg, 0.0, 1.0), 0.0)
    f_g_var = np.where(active, np.maximum(f_g_var, 0.0), 0.0)
    return SimplexSolution(f_rna_pos=f_pos, f_rna_neg=f_neg, f_g=f_g, f_g_var=f_g_var)
