"""The readable statements of ψ's pieces that the gates recompute with — the oracles, not a second solver.

ψ is native (`native/psi_kernel.cpp`, read through `simplex_logodds.psi_cube`); the gates that hold it to its
derivations need the derivations written down once, in numpy, small enough to read against `EQUATIONS.md`:
the three-component strand term with the variance frozen at a reference composition, the two Jeffreys arms,
the delivered row's map onto the cube, a log-sum-exp, and the delivery table stated by hand. Nothing here
asserts and nothing in `src/` reads it.
"""

from __future__ import annotations

import numpy as np
from scipy.special import log_expit, logsumexp

from rigel.calibration.simplex_logodds import CubeRows


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
