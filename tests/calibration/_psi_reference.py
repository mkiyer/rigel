"""The readable statements of ψ's pieces that the gates recompute with — the oracles, not a second solver.

ψ is native (`native/psi_kernel.h`); the four readers at the end of this file (`psi_cube`, `gdna_arm`,
`posterior_median_fg`, `compose`) call its gate bindings. The gates that hold it to its derivations need
the derivations written down once, in numpy, small enough to read against `EQUATIONS.md`: the
three-component strand term with the variance frozen at a reference composition and its two-component
special case (the collapse `test_strand_likelihood_reference.py` gates), the two Jeffreys arms,
the delivered row's map onto the cube, a log-sum-exp, and the delivery table stated by hand. Nothing here
asserts and nothing in `src/` reads it.
"""

from __future__ import annotations

import numpy as np
from scipy.special import log_expit, logsumexp

from rigel.calibration.simplex_logodds import (
    _TILT_NODES,
    CubeRows,
    _arm,
    _cube_args,
    _prior,
    _reference_composition,
)
from rigel.native import psi_compose, psi_cube_native, psi_gdna_arm, psi_posterior_median


def strand_loglik_mixture(
    u_pos, n, f_g, f_pos, f_neg, kappa, od_g, od_r, f_g_ref, f_pos_ref, f_neg_ref
):
    """The three-component gDNA / RNA₊ / RNA₋ strand log-likelihood of the + column count ``u_pos`` at the
    mean ``n·p``, ``p = ½·f_g + κ·f₊ + (1−κ)·f₋``, with the variance frozen at the REFERENCE composition
    (the count-zero-information freeze: the count sets precision, never composition). Broadcasts, and keeps
    the dtype of its inputs (the replay's budget gate evaluates it in float32)."""
    p = 0.5 * f_g + kappa * f_pos + (1.0 - kappa) * f_neg
    mean = n * p
    rscale = kappa * (1.0 - kappa)
    p_ref = 0.5 * f_g_ref + kappa * f_pos_ref + (1.0 - kappa) * f_neg_ref
    var = (
        n * p_ref * (1.0 - p_ref)
        + (n * f_g_ref) ** 2 * 0.25 * od_g
        + (n * f_pos_ref) ** 2 * rscale * od_r
        + (n * f_neg_ref) ** 2 * rscale * od_r
    )
    var = np.maximum(var, 1.0e-9)
    return -0.5 * (u_pos - mean) ** 2 / var - 0.5 * np.log(var)


def strand_loglik(
    gdna_frac: np.ndarray,
    sense: float,
    antisense: float,
    rna_sense_frac: float,
    *,
    gdna_strand_overdispersion: float = 0.0,
    rna_strand_overdispersion: float = 0.0,
) -> np.ndarray:
    """The TWO-component gDNA/RNA strand log-likelihood of one region over a ``gdna_frac`` grid — the readable
    special case ψ's three-component term (`strand_loglik_mixture`, the kernel's `strand_term`) collapses onto
    when one RNA strand is dead; the gate is `test_strand_likelihood_reference.py`.

    Of ``N = sense + antisense`` discrete unspliced fragments, a fraction ``gdna_frac`` are gDNA
    (oriented-sense rate ½, intra-class correlation ``gdna_strand_overdispersion``) and
    ``1 − gdna_frac`` are RNA (oriented-sense rate ``rna_sense_frac``, intra-class correlation
    ``rna_strand_overdispersion``). The mixture sense count has mean ``N·p`` and a variance in
    three parts — the Binomial mixture variance plus the excess variance each component's shared
    per-region sense rate contributes, scaled by that component's own ``μ_c(1−μ_c)``: the
    ``N·gdna_frac`` gDNA fragments (mean ½) add ``(N·gdna_frac)²·¼·gdna_strand_overdispersion`` and
    the ``N·(1−gdna_frac)`` RNA fragments (mean κ) add
    ``(N·(1−gdna_frac))²·κ(1−κ)·rna_strand_overdispersion``. Normal moment approximation::

        p   = ½·gdna_frac + rna_sense_frac·(1 − gdna_frac);  mean = N·p
        var = N·p·(1 − p)
            + (N·gdna_frac)²·¼·gdna_strand_overdispersion
            + (N·(1 − gdna_frac))²·κ(1 − κ)·rna_strand_overdispersion
        loglik(gdna_frac) = −½·(sense − mean)² / var  −  ½·log(var)

    Limits: ``gdna_frac → 1`` ⇒ Beta-Binomial(½, od_gdna); ``gdna_frac → 0`` ⇒ Beta-Binomial(κ,
    od_rna); both ``od → 0`` ⇒ the Binomial mixture exactly. Symmetry: at ``κ = ½`` the RNA scale
    ``κ(1−κ)`` equals the gDNA ¼, so with
    ``od_gdna = od_rna`` the variance is flat in ``gdna_frac`` (the means coincide and only the
    ``g² + (1−g)²`` scaling depends on ``gdna_frac``) — an unstranded region is uninformative. Each
    component's excess variance uses its own mean ``μ_c(1−μ_c)``; the normal-moment vs exact-mixture
    discrepancy is small and ~constant in ``N``.
    """
    # Fully elementwise in (sense, antisense, gdna_frac), so it broadcasts: scalar (sense, antisense)
    # with a 1-D grid returns one region's curve; column (sense, antisense) of shape (K, 1) with a row
    # grid of shape (1, n_grid) returns the whole (K, n_grid) batch at once.
    n = sense + antisense
    p = 0.5 * gdna_frac + rna_sense_frac * (1.0 - gdna_frac)
    mean = n * p
    rna_var_scale = rna_sense_frac * (1.0 - rna_sense_frac)  # κ(1−κ); the RNA component's μ(1−μ)
    var = (
        n * p * (1.0 - p)
        + (n * gdna_frac) ** 2 * 0.25 * gdna_strand_overdispersion
        + (n * (1.0 - gdna_frac)) ** 2 * rna_var_scale * rna_strand_overdispersion
    )
    var = np.maximum(var, 1.0e-9)
    return -0.5 * (sense - mean) ** 2 / var - 0.5 * np.log(var)


def jeffreys_arms(lam, c_g: float = 0.5, c_r: float = 0.5):
    """The two reference arms over the λ grid, ``c_g·log f_g + c_r·log(1 − f_g)`` — ``½`` and ``½`` is
    ψ's `_JEFFREYS_REF`; other exponents are the gates' ablations, delivered to ψ as the DIFFERENCE from
    the reference on the λ-row channel (the arms are additive, so an exponent change is a row)."""
    lam = np.asarray(lam, np.float64)
    return c_g * log_expit(lam) + c_r * log_expit(-lam)


def row_at(row, fg, tau):
    """One delivered row (a record with ``profile_pos``, ``profile_neg``, ``u``, ``total``, ``opportunity``,
    ``rho_ref``) evaluated over ψ's cells — at each ``(λ, θ)`` the strand's share
    ``f_s = (1 − f_g)(1 ± τ)/2`` implies the density ``f_s·n/a_r``, and the held profile is read at
    ``log(ρ_s/ρ_ref)``; the strands' rows add and the result is max-normalised over the whole row (the
    kernel's map, `EQUATIONS.md` §9e). ``fg`` is ``(K,)``; ``tau`` ``(K_t,)`` or ``(K, K_t)``."""
    fg = np.asarray(fg, np.float64)
    tau = np.asarray(tau, np.float64)
    if tau.ndim == 1:
        tau = tau[None, :]
    f_act = (1.0 - fg)[:, None]
    out = np.zeros(np.broadcast_shapes(f_act.shape, tau.shape))
    scale = float(row.total) / float(row.opportunity)
    for prof, sign in ((row.profile_pos, 1.0), (row.profile_neg, -1.0)):
        if prof is None:
            continue
        prof = np.asarray(prof, np.float64)
        with np.errstate(divide="ignore"):
            u_s = np.log(f_act * (1.0 + sign * tau) / 2.0 * scale) - np.log(float(row.rho_ref))
        out += np.interp(u_s, np.asarray(row.u, np.float64), prof, left=prof[0], right=prof[-1])
    return out - out.max()


def lse(a, axis, keepdims=False):
    """log Σ exp along ``axis``; an all-``−∞`` slice gives ``−∞``."""
    return logsumexp(np.asarray(a, np.float64), axis=axis, keepdims=keepdims)


def cube_rows_of(rows: dict, u) -> CubeRows:
    """A `CubeRows` table from ``{slot: (profile_pos | None, profile_neg | None, total, opportunity,
    rho_ref)}`` on the grid ``u`` — how a gate states a delivery by hand."""
    u = np.asarray(u, np.float64)
    slots = sorted(rows)
    t = CubeRows.blank(len(slots), u)
    for r, k in enumerate(slots):
        pos, neg, total, opportunity, rho_ref = rows[k]
        t.slot[r] = k
        if pos is not None:
            t.profile_pos[r] = pos
            t.has_pos[r] = True
        if neg is not None:
            t.profile_neg[r] = neg
            t.has_neg[r] = True
        t.total[r], t.opportunity[r], t.rho_ref[r] = total, opportunity, rho_ref
    return t


class Row:
    """One delivered row as a record, for `row_at` — what a gate reads out of a `CubeRows` table."""

    def __init__(self, table: CubeRows, r: int):
        self.profile_pos = table.profile_pos[r] if table.has_pos[r] else None
        self.profile_neg = table.profile_neg[r] if table.has_neg[r] else None
        self.u = table.u
        self.total, self.opportunity, self.rho_ref = (
            table.total[r],
            table.opportunity[r],
            table.rho_ref[r],
        )


def row_record(profile_pos, profile_neg, u, total, opportunity, rho_ref) -> Row:
    """One delivered row stated by hand, for `row_at`."""
    return Row(cube_rows_of({0: (profile_pos, profile_neg, total, opportunity, rho_ref)}, u), 0)


# ── the kernel's own ψ pieces, read through its gate bindings (`native/solve_kernel.cpp`) ───────────────


def gdna_arm(log_rho, logP, lam, mass, eff):
    """The fitted gDNA arm as an ``(m, K)`` matrix — the kernel's own construction, for the gates: the curve
    ``(log_rho, logP)`` read at ``log σ(λ) + log M − log E`` for every slot and cell (`native.psi_gdna_arm`)."""
    lam = np.ascontiguousarray(lam, np.float64)
    mass = np.ascontiguousarray(mass, np.float64)
    out = np.zeros((mass.shape[0], lam.shape[0]))
    psi_gdna_arm(
        np.ascontiguousarray(log_rho, np.float64),
        np.ascontiguousarray(logP, np.float64),
        lam,
        mass,
        np.ascontiguousarray(eff, np.float64),
        out,
    )
    return out


def psi_cube(
    u_pos,
    u_neg,
    allow_pos,
    allow_neg,
    fg_ref,
    fpos_ref,
    fneg_ref,
    *,
    kappa,
    od_g,
    od_r,
    lam,
    ambig: bool,
    gdna_prior=None,
    gdna_support=None,
    row_logprior=None,
    cube_rows=None,
    n_tilt: int | None = None,
):
    """ψ ITSELF over the ``(λ, θ)`` cube for ``m`` slots of one class — the strand term + the two Jeffreys
    arms (``_JEFFREYS_REF``) + the fitted gDNA arm (the curve ``gdna_prior`` read on ``gdna_support`` at each
    cell's density) + the delivered composition rows, added per cell (+ the delivered cube rows) + the
    θ quadrature's log-weights — as ``(m, K, C)`` in float64, with the two strand-fraction grids it was
    evaluated on, ``(f_pos, f_neg)``, and the tilt ``tau`` they were built from. A single-strand call
    (``ambig=False``) has one column, the tilt of each slot's live strand (``τ = ±1``) and no weight. An
    AMBIG call places the θ nodes across each slot's strand term (``n_tilt`` of them, ``_TILT_NODES`` by
    derivation) for the MIXED hypothesis and appends the two PURE hypotheses as the columns ``τ = +1, −1``
    (the tilt atom: the three at equal reference weight, a delivered level on a strand ruling the other
    strand's atom out). The same code the solve runs (`native.psi_cube_native`); the gates read it here."""
    lam = np.ascontiguousarray(lam, np.float64)
    K = lam.shape[0]
    u_pos = np.ascontiguousarray(u_pos, np.float64)
    u_neg = np.ascontiguousarray(u_neg, np.float64)
    ap, an = np.ascontiguousarray(allow_pos, bool), np.ascontiguousarray(allow_neg, bool)
    fg_ref, fpos_ref, fneg_ref = _reference_composition(ap, an, fg_ref, fpos_ref, fneg_ref)
    n_tilt = int(_TILT_NODES if n_tilt is None else n_tilt)
    m, C = u_pos.shape[0], (n_tilt + 2 if ambig else 1)
    psi, f_pos, f_neg, tau = (np.zeros((m, K, C)) for _ in range(4))
    psi_cube_native(
        u_pos=u_pos,
        u_neg=u_neg,
        allow_pos=ap,
        allow_neg=an,
        fg_ref=fg_ref,
        fpos_ref=fpos_ref,
        fneg_ref=fneg_ref,
        kappa=float(kappa),
        od_g=float(od_g),
        od_r=float(od_r),
        lam=lam,
        gdna=_arm(gdna_prior, gdna_support, m),
        row_logprior=_prior(row_logprior, K),
        **_cube_args(cube_rows, lam),
        n_tilt=n_tilt,
        ambig=bool(ambig),
        out_psi=psi,
        out_fpos=f_pos,
        out_fneg=f_neg,
        out_tau=tau,
    )
    return psi, f_pos, f_neg, tau


def posterior_median_fg(post, lam):
    """Per-slot point estimate of ``f_g``: the posterior's ½-QUANTILE, read off the CDF on the uniform λ
    lattice — a continuous quantile (the grid mass as a histogram with edges at the midpoints, the crossing
    bin interpolated ON λ, then mapped through σ: median equivariance, which is why ``f_g`` is a median
    and not a mean). ``post``: ``(m, K)``; ``lam``: ``(K,)``. Returns ``(m,)``."""
    post = np.ascontiguousarray(post, np.float64)
    out = np.zeros(post.shape[0])
    psi_posterior_median(post, np.ascontiguousarray(lam, np.float64), out)
    return out


def compose(f_g, w_pos, allow_pos, allow_neg):
    """ψ's composition as the MAP from its two parameters — ``f_g`` (the level, a median) and ``w_pos``
    (the + strand's share of the RNA) — onto the simplex: ``f_pos = (1 − f_g)·w``, ``f_neg = (1 − f_g)·(1 −
    w)``, so closure is structural; the share is clamped to ``[0, 1]`` and restricted to the admissible
    strands (a single-strand slot's whole RNA sits on its live strand whatever the share says; a slot with
    neither strand has no RNA to place). Returns ``(f_pos, f_neg)``."""
    f_g = np.ascontiguousarray(f_g, np.float64)
    ap, an = np.ascontiguousarray(allow_pos, bool), np.ascontiguousarray(allow_neg, bool)
    f_pos, f_neg = np.zeros(f_g.shape[0]), np.zeros(f_g.shape[0])
    psi_compose(f_g, np.ascontiguousarray(w_pos, np.float64), ap, an, f_pos, f_neg)
    return f_pos, f_neg
