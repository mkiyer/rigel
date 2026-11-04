"""Readable statements of ψ's pieces that the tests recompute with: the oracles, not a second solver.

ψ is native (`native/psi_kernel.h`); the readers at the end of this module call it: `psi_solve` through the solve's own
binding with every option the kernel takes, and `psi_cube`, `gdna_arm`, `posterior_median_fg` and `compose` through its
test bindings. The tests that hold it to its derivations need those
derivations written once in numpy: the three-component strand term with the variance frozen at a reference composition
and its two-component special case (the collapse `test_strand_likelihood_reference.py` checks), the two Jeffreys arms,
the delivered row's map onto the cube, a log-sum-exp, and the delivery table stated by hand. Nothing here asserts and
nothing in `src/` reads it.
"""

from __future__ import annotations

import numpy as np
from scipy.special import log_expit, logsumexp

from rigel.calibration.region_chain import RegionDeconv
from rigel.calibration.simplex_logodds import (
    TILT_NODES,
    CubeRows,
    logodds_grid,
    reference_composition,
)
from rigel.native import psi_solve as psi_solve_native
from rigel._solve_impl.testing import compose as psi_compose
from rigel._solve_impl.testing import gdna_arm as psi_gdna_arm
from rigel._solve_impl.testing import posterior_median as psi_posterior_median
from rigel._solve_impl.testing import psi_cube as psi_cube_native


def strand_loglik_mixture(u_pos, n, f_g, f_pos, f_neg, kappa, f_g_ref, f_pos_ref, f_neg_ref):
    """The three-component gDNA / RNA₊ / RNA₋ strand log-likelihood of the + column count ``u_pos`` at the mean ``n·p``,
    ``p = ½·f_g + κ·f₊ + (1−κ)·f₋``, with the binomial variance frozen at the reference composition (the count sets
    precision, never composition). Broadcasts, and keeps the dtype of its inputs."""
    p = 0.5 * f_g + kappa * f_pos + (1.0 - kappa) * f_neg
    mean = n * p
    p_ref = 0.5 * f_g_ref + kappa * f_pos_ref + (1.0 - kappa) * f_neg_ref
    var = np.maximum(n * p_ref * (1.0 - p_ref), 1.0e-9)
    return -0.5 * (u_pos - mean) ** 2 / var - 0.5 * np.log(var)


def strand_loglik(
    gdna_frac: np.ndarray, sense: float, antisense: float, rna_sense_frac: float
) -> np.ndarray:
    """The two-component gDNA/RNA strand log-likelihood of one region over a ``gdna_frac`` grid: the special case ψ's
    three-component term (`strand_loglik_mixture`, the kernel's `strand_term`) collapses onto when one RNA strand is
    dead; `test_strand_likelihood_reference.py` checks the collapse.

    Of ``N = sense + antisense`` unspliced fragments, a fraction ``gdna_frac`` are gDNA (oriented-sense rate ½) and ``1
    − gdna_frac`` RNA (oriented-sense rate ``rna_sense_frac``). Normal moment approximation of the Binomial mixture::

        p   = ½·gdna_frac + rna_sense_frac·(1 − gdna_frac);  mean = N·p;  var = N·p·(1 − p)
        loglik(gdna_frac) = −½·(sense − mean)² / var  −  ½·log(var)
    """
    # Fully elementwise in (sense, antisense, gdna_frac), so it broadcasts: scalar (sense, antisense)
    # with a 1-D grid returns one region's curve; column (sense, antisense) of shape (K, 1) with a row
    # grid of shape (1, n_grid) returns the whole (K, n_grid) batch at once.
    n = sense + antisense
    p = 0.5 * gdna_frac + rna_sense_frac * (1.0 - gdna_frac)
    var = np.maximum(n * p * (1.0 - p), 1.0e-9)
    return -0.5 * (sense - n * p) ** 2 / var - 0.5 * np.log(var)


def jeffreys_arms(lam, c_g: float = 0.5, c_r: float = 0.5):
    """The two reference arms over the λ grid, ``c_g·log f_g + c_r·log(1 − f_g)``. ``½`` and ``½`` is ψ's
    `JEFFREYS_REF`; other exponents are the tests' ablations, delivered to ψ as the difference from the reference on the
    λ-row channel (the arms are additive, so an exponent change is a row)."""
    lam = np.asarray(lam, np.float64)
    return c_g * log_expit(lam) + c_r * log_expit(-lam)


def row_at(row, fg, tau):
    """One delivered row (a record with ``profile_pos``, ``profile_neg``, ``u``, ``total``, ``opportunity``,
    ``rho_ref``) evaluated over ψ's cells: at each ``(λ, θ)`` the strand's share ``f_s = (1 − f_g)(1 ± τ)/2`` implies
    the density ``f_s·n/a_r``, and the held profile is read at ``log(ρ_s/ρ_ref)``; the strands' rows add and the result
    is max-normalised over the whole row (the kernel's map). ``fg`` is ``(K,)``; ``tau`` ``(K_t,)`` or ``(K, K_t)``."""
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
    """A `CubeRows` table from ``{slot: (profile_pos | None, profile_neg | None, total, opportunity, rho_ref)}`` on the
    grid ``u``: how a test states a delivery by hand."""
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
    """One delivered row as a record, for `row_at`: what a test reads out of a `CubeRows` table."""

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


# ── the kernel's own ψ, with every option it takes ──────────────────────────────────────────────


def _cube_args(cube_rows, u) -> dict:
    """The delivered rows as ψ's kernel takes them — the table's arrays, or an empty table's on this grid."""
    if cube_rows is None:
        return CubeRows.blank(0, u).kernel_args()
    if cube_rows.u.shape != u.shape or not np.array_equal(cube_rows.u, u):
        raise ValueError("the delivered rows are not on the solve grid")
    return cube_rows.kernel_args()


def _reference(allow_pos, allow_neg, fg_ref, fpos_ref, fneg_ref):
    """The composition the strand variance is frozen at: the supplied one (the sweep's incoming belief), else the
    solve's own default (`simplex_logodds.reference_composition`)."""
    if fg_ref is None or fpos_ref is None or fneg_ref is None:
        return reference_composition(allow_pos, allow_neg)
    return tuple(np.ascontiguousarray(a, np.float64) for a in (fg_ref, fpos_ref, fneg_ref))


def _prior(rows, K: int):
    """A per-slot ``(m, K)`` log-prior row array on the solve grid as the kernel reads it, or ``None``."""
    if rows is None:
        return None
    rows = np.ascontiguousarray(rows, np.float64)
    if rows.ndim != 2 or rows.shape[1] != K:
        raise ValueError(
            f"a ψ prior must be (m, K) on the solve grid; got {rows.shape} for K = {K}"
        )
    return rows


def _arm(gdna_prior, gdna_support, m: int):
    """The fitted gDNA arm as the kernel takes it, ``(log_rho, logP, mass, eff)``, or ``None`` for the reference alone.

    ``(log_rho, logP)`` is the landscape's curve, a density in log-rate, and ``(mass, eff)`` the per-slot support it
    is read on; ψ evaluates it at each cell's ``log ρ_c = log f_g + log M − log E``, ends held constant.
    """
    if gdna_prior is None:
        if gdna_support is not None:
            raise ValueError(
                "a gDNA support without a curve: the arm is the curve read on the support"
            )
        return None
    if gdna_support is None:
        raise ValueError("the fitted gDNA arm needs its per-slot support (mass, eff) to be read on")
    log_rho, logP = (np.ascontiguousarray(a, np.float64) for a in gdna_prior)
    mass, eff = (np.ascontiguousarray(a, np.float64) for a in gdna_support)
    if log_rho.ndim != 1 or log_rho.shape != logP.shape or log_rho.shape[0] < 1:
        raise ValueError("the gDNA arm's curve must be two (G,) arrays")
    if mass.shape != (m,) or eff.shape != (m,):
        raise ValueError(
            f"the gDNA arm's support must be (mass, eff) of shape ({m},); got {mass.shape} and {eff.shape}"
        )
    return log_rho, logP, mass, eff


def psi_solve(
    u_pos,
    u_neg,
    allow_pos,
    allow_neg,
    *,
    kappa,
    n_grid,
    L: float,
    gdna_prior=None,
    gdna_support=None,
    row_logprior=None,
    fg_ref=None,
    fpos_ref=None,
    fneg_ref=None,
    cube_rows=None,
    n_tilt: int | None = None,
    n_threads: int = 1,
) -> RegionDeconv:
    """`simplex_logodds.solve_regions_logodds` with every option ψ's kernel takes (the binding `native.psi_solve`);
    with none of them it is that solve.

    ``gdna_prior`` is the landscape's curve ``(log_rho, logP)`` and ``gdna_support`` the per-slot ``(mass, eff)`` it is
    read on (`_arm`). ``row_logprior`` is the policy's delivered composition rows, ``(m, K)`` on the λ grid;
    ``cube_rows`` is the RNA level lanes' delivery (`CubeRows`). ``fg_ref``/``fpos_ref``/``fneg_ref`` is the reference
    composition the strand variance is frozen at (`_reference`). ``n_tilt`` defaults to ``TILT_NODES``.
    """
    u_pos = np.ascontiguousarray(u_pos, np.float64)
    u_neg = np.ascontiguousarray(u_neg, np.float64)
    m = u_pos.shape[0]
    ap, an = np.ascontiguousarray(allow_pos, bool), np.ascontiguousarray(allow_neg, bool)
    fg_ref, fpos_ref, fneg_ref = _reference(ap, an, fg_ref, fpos_ref, fneg_ref)
    slots = np.flatnonzero((ap | an) & ((u_pos + u_neg) > 0.0)).astype(np.int64)
    lam, _fg = logodds_grid(int(n_grid), L)
    out = {k: np.zeros(m) for k in ("fg", "fp", "fn", "vg")}
    if slots.size:
        psi_solve_native(
            slots=slots,
            u_pos=u_pos,
            u_neg=u_neg,
            allow_pos=ap,
            allow_neg=an,
            fg_ref=fg_ref,
            fpos_ref=fpos_ref,
            fneg_ref=fneg_ref,
            kappa=float(kappa),
            lam=lam,
            gdna=_arm(gdna_prior, gdna_support, m),
            row_logprior=_prior(row_logprior, lam.shape[0]),
            **_cube_args(cube_rows, lam),
            n_tilt=int(TILT_NODES if n_tilt is None else n_tilt),
            out_fg=out["fg"],
            out_fpos=out["fp"],
            out_fneg=out["fn"],
            out_var=out["vg"],
            n_threads=int(n_threads),
        )
    return RegionDeconv(
        gdna_frac=out["fg"],
        rna_pos_frac=out["fp"],
        rna_neg_frac=out["fn"],
        gdna_frac_var=out["vg"],
    )


# ── the kernel's own ψ pieces, read through its test bindings (`native/testing/solve.cpp`) ────────


def gdna_arm(log_rho, logP, lam, mass, eff):
    """The fitted gDNA arm as an ``(m, K)`` matrix, the kernel's own construction (`gdna_arm`): the curve ``(log_rho,
    logP)`` read at ``log σ(λ) + log M − log E`` for every slot and cell."""
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
    lam,
    ambig: bool,
    gdna_prior=None,
    gdna_support=None,
    row_logprior=None,
    cube_rows=None,
    n_tilt: int | None = None,
):
    """ψ itself over the ``(λ, θ)`` cube for ``m`` slots of one class: the strand term, the two Jeffreys arms
    (``JEFFREYS_REF``), the fitted gDNA arm (the curve ``gdna_prior`` read on ``gdna_support`` at each cell's density),
    the delivered composition rows added per cell, the delivered cube rows, and the θ quadrature's log-weights, as ``(m,
    K, C)`` in float64, with the two strand-fraction grids it was evaluated on, ``(f_pos, f_neg)``, and the tilt ``tau``
    they were built from.

    A single-strand call (``ambig=False``) has one column, the tilt of each slot's live strand (``τ = ±1``) and no
    weight. An AMBIG call places ``n_tilt`` θ nodes (``TILT_NODES`` by derivation) across each slot's strand term for
    the mixed hypothesis and appends the two pure hypotheses as the columns ``τ = +1, −1`` (the three at equal reference
    weight, a delivered level on a strand ruling the other strand's atom out). It is the code the solve runs (the
    kernel's `psi_cube`).
    """
    lam = np.ascontiguousarray(lam, np.float64)
    K = lam.shape[0]
    u_pos = np.ascontiguousarray(u_pos, np.float64)
    u_neg = np.ascontiguousarray(u_neg, np.float64)
    ap, an = np.ascontiguousarray(allow_pos, bool), np.ascontiguousarray(allow_neg, bool)
    fg_ref, fpos_ref, fneg_ref = _reference(ap, an, fg_ref, fpos_ref, fneg_ref)
    n_tilt = int(TILT_NODES if n_tilt is None else n_tilt)
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
    """Per-slot point estimate of ``f_g``: the posterior's ½-quantile, read off the CDF on the uniform λ lattice as a
    continuous quantile (the grid mass as a histogram with edges at the midpoints, the crossing bin interpolated on λ,
    then mapped through σ; the median is equivariant, which is why ``f_g`` is a median and not a mean). ``post``: ``(m,
    K)``; ``lam``: ``(K,)``. Returns ``(m,)``."""
    post = np.ascontiguousarray(post, np.float64)
    out = np.zeros(post.shape[0])
    psi_posterior_median(post, np.ascontiguousarray(lam, np.float64), out)
    return out


def compose(f_g, w_pos, allow_pos, allow_neg):
    """ψ's composition as the map from its two parameters, ``f_g`` (the level, a median) and ``w_pos`` (the + strand's
    share of the RNA), onto the simplex: ``f_pos = (1 − f_g)·w``, ``f_neg = (1 − f_g)·(1 − w)``, so closure is
    structural. The share is clamped to ``[0, 1]`` and restricted to the admissible strands (a single-strand slot's
    whole RNA sits on its live strand whatever the share says; a slot with neither strand has no RNA to place). Returns
    ``(f_pos, f_neg)``."""
    f_g = np.ascontiguousarray(f_g, np.float64)
    ap, an = np.ascontiguousarray(allow_pos, bool), np.ascontiguousarray(allow_neg, bool)
    f_pos, f_neg = np.zeros(f_g.shape[0]), np.zeros(f_g.shape[0])
    psi_compose(f_g, np.ascontiguousarray(w_pos, np.float64), ap, an, f_pos, f_neg)
    return f_pos, f_neg
