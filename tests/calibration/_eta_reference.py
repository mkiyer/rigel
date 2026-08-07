"""THE ``eta`` TRANSFER KERNEL — the executable form of the transfer derivation.

⭐ ``eta = lambda - log(E_g/E_r) = log(rho_g/rho_R)`` is the FRAME-FREE composition coordinate: it crosses
any EDGE whose two endpoints draw on the SAME populations as the **identity**, so the destination's own
opportunities do the conversion exactly and belief-free.

This module is the small pile of algebra that claim rests on, kept in one place so the gates in
``test_eta_transfer.py`` and the ``scripts/design`` prototype cannot drift apart (``TRAPS.md`` A11 — one
home, every consumer imports it). It follows the precedent of ``_oracle.py``: a test-side module that
``scripts/design`` imports, which is what keeps ``scripts/`` un-importable-from-production while still
letting a prototype and its gates share one definition.

⛔ **This is a PROTOTYPE's kernel, not a fixture.** When the rebuild moves into ``src/`` this file is
deleted and the gates re-point at the production module.

⚠ **THREE POPULATIONS AND THERE IS NO FOURTH** — ``gDNA``, ``RNA+``, ``RNA-``. :func:`population_set`
is that axiom made executable; ``|T| <= 3`` is asserted, not assumed.
"""

from __future__ import annotations

import numpy as np
from scipy.special import digamma, polygamma

__all__ = [
    "GDNA",
    "RNA_NEG",
    "RNA_POS",
    "cross_composition",
    "eta_from_lambda",
    "lambda_from_eta",
    "log_rate_posterior",
    "population_set",
    "rebuild_densities",
    "residual_eta_scalar",
    "share_from_density",
    "shares_from_lambda",
    "tilt_angle",
]

#: Column order of the population indicator returned by :func:`population_set`.
GDNA, RNA_POS, RNA_NEG = 0, 1, 2

_EPS = 1.0e-300  # guards a log of an exact zero only; never a modelling threshold


def population_set(free_pos, free_neg, eff_g, eff_r) -> np.ndarray:
    """``T(slot)`` as an ``(n, 3)`` boolean indicator over ``[gDNA, RNA+, RNA-]``.

    gDNA is genomically continuous, so it is admissible wherever it has OPPORTUNITY — there is no gDNA
    analogue of a forbidden strand. Each RNA strand is admissible iff the annotation frees it there AND
    there is RNA opportunity. A slot shorter than a fragment has ``E = 0`` and holds no population at all.

    ⭐ These are exactly the solver's own ``live`` predicates, so this is a READ of two bits and two
    opportunities, never a new model — which is what makes ``|T| <= 3`` structural.
    """
    eff_g = np.asarray(eff_g, dtype=np.float64)
    eff_r = np.asarray(eff_r, dtype=np.float64)
    out = np.empty(eff_g.shape + (3,), dtype=bool)
    out[..., GDNA] = eff_g > 0.0
    out[..., RNA_POS] = np.asarray(free_pos, dtype=bool) & (eff_r > 0.0)
    out[..., RNA_NEG] = np.asarray(free_neg, dtype=bool) & (eff_r > 0.0)
    return out


def log_opportunity_ratio(eff_g, eff_r) -> np.ndarray:
    """``log(E_g/E_r)`` — the whole of the conversion between ``lambda`` and ``eta``.

    A known geometric constant read off the index: no belief, no count, no ratio of totals.
    """
    return np.log(np.maximum(np.asarray(eff_g, dtype=np.float64), _EPS)) - np.log(
        np.maximum(np.asarray(eff_r, dtype=np.float64), _EPS)
    )


def eta_from_lambda(lam, eff_g, eff_r) -> np.ndarray:
    """``eta = lambda - log(E_g/E_r) = log(rho_g/rho_R)`` — into the frame-free coordinate."""
    return np.asarray(lam, dtype=np.float64) - log_opportunity_ratio(eff_g, eff_r)


def lambda_from_eta(eta, eff_g, eff_r) -> np.ndarray:
    """``lambda = eta + log(E_g/E_r)`` — back into a slot's own share frame.

    ⚠ ``eta`` and ``lambda`` differ by a CONSTANT, so their variances are equal: a precision crosses an
    edge as the identity too, and there is no scale noise to charge for the conversion.
    """
    return np.asarray(eta, dtype=np.float64) + log_opportunity_ratio(eff_g, eff_r)


def shares_from_lambda(lam) -> tuple[np.ndarray, np.ndarray]:
    """``(phi_g, phi_R)`` from ``lambda`` — the logistic, written stably."""
    lam = np.asarray(lam, dtype=np.float64)
    phi_g = np.where(lam >= 0.0, 1.0 / (1.0 + np.exp(-lam)), np.exp(lam) / (1.0 + np.exp(lam)))
    return phi_g, 1.0 - phi_g


def rebuild_densities(eta, mass, eff_g, eff_r):
    """``(rho_g, rho_R)`` at a slot, from ``eta`` and the slot's OWN ``(M, E_g, E_r)``.

    ⭐⭐ ``sum_c rho_c E_c = M`` holds BY CONSTRUCTION for any ``eta`` — which is the whole reason the mass
    pin has nothing left to restore.
    """
    mass = np.asarray(mass, dtype=np.float64)
    phi_g, phi_r = shares_from_lambda(lambda_from_eta(eta, eff_g, eff_r))
    eff_g = np.maximum(np.asarray(eff_g, dtype=np.float64), _EPS)
    eff_r = np.maximum(np.asarray(eff_r, dtype=np.float64), _EPS)
    return phi_g * mass / eff_g, phi_r * mass / eff_r


def log_rate_posterior(count, eff):
    """The HONEST density claim: ``a`` fragments over opportunity ``E`` under the Jeffreys rate prior.

    ``rho ~ Gamma(a + 1/2, E)``, so ``E[log rho] = psi0(a + 1/2) - log E`` and
    ``Var(log rho) = psi1(a + 1/2)`` — finite at ``a = 0``, where the location is
    ``psi0(1/2) = -1.9635`` and the claim is ``rho = 0.1404 / E``.

    ⭐ Length lives in the LOCATION (``-log E``), never in the variance: a zero count is always "a factor
    of ~9 either way", and the length sets what it is a factor of nine around. No exact zero is ever
    formed, so a multiplicative transport has nothing unrepresentable to carry.
    """
    a = np.asarray(count, dtype=np.float64)
    eff = np.maximum(np.asarray(eff, dtype=np.float64), _EPS)
    return digamma(a + 0.5) - np.log(eff), polygamma(1, a + 0.5)


#: The saturation point of :func:`residual_eta_scalar` when the shared claim already exceeds the
#: destination's own count. 40 nats is ``f_g = 1 - 4e-18``, i.e. all-gDNA to well below float resolution —
#: a representable stand-in for the vertex, not a tuned threshold. ⛔ It is deliberately NOT a shift.
_SATURATED = 40.0


def residual_eta_scalar(rho_g_src, mass_dst, eff_g_dst, eff_r_dst) -> float:
    """The ``|N| >= 1`` branch for ONE hop, in plain floats.

    ⚠ The scalar twin of :func:`cross_composition`'s mismatch branch. It exists because the sweep's scan is
    a sequential Python loop where building a one-element array per hop costs more than the arithmetic —
    the same trade `bp_solver` records for its own scalar relay. ⛔ ``test_eta_transfer`` gates the two
    forms against each other over a randomised sweep, so ONE HOME holds the definition and the fast path
    cannot drift from it (`TRAPS.md` A11).
    """
    residual = float(mass_dst) - float(rho_g_src) * float(eff_g_dst)
    if residual <= 0.0:
        return _SATURATED
    return float(
        np.log(max(float(rho_g_src), _EPS)) - np.log(max(residual / float(eff_r_dst), _EPS))
    )


def tilt_angle(rho_pos, rho_neg, prec_pos, prec_neg):
    """``(theta, precision)`` — the strand tilt in **psi's own coordinate**, from two RNA densities.

    ⛔⛔ **THE TILT IS AN ANGLE AND NOT A LOG-ODDS, AND psi'S GRID SAYS SO.**
    ``simplex_logodds._tilt_grid`` spans exactly ``[-pi/2, +pi/2]``, because the coordinate is
    ``theta = arcsin(tau)`` with ``tau = (f_pos - f_neg)/f_R``. A mode outside that interval is not a
    claim psi can represent at all: the message term ``-1/2 * p * (theta - m)^2`` is then MONOTONE over
    the whole grid, so it pins the tilt at the boundary — ``tau = +-1``, "all RNA on one strand" — and
    the strand likelihood, no longer free to integrate the nuisance out, explains whatever it then cannot
    fit by calling the mass gDNA. Measured at ``g00 ss0.99 capture_on``: a raw log-odds delivered modes
    of ``+-4.6`` against a grid ending at 1.571, and muting that ONE channel at the psi boundary took
    ``Sum|err|`` from 1,751,145 to 452,326 — **74 %** of the arm's error there.

    ⛔⛔ **AND IT MUST BE BUILT FROM DECONVOLVED RNA DENSITIES, NOT FROM RAW COUNTS.** psi's ``tau`` is
    the tilt of the **RNA**; a count tilt is the tilt of the **MASS**, and F2 relates them by
    ``tau_obs = (2*kappa - 1) * (1 - f_g) * tau`` — which on an antisense protocol INVERTS THE SIGN
    (``kappa = 0.0101`` at ``ss = 0.99`` gives ``2*kappa - 1 = -0.98``). Measured: correcting the
    coordinate alone, still on raw counts, left ``g00`` at 1,568,412 — no better — because it then
    asserted the right angle with the wrong sign at a precision of ``n``. ⭐ ``rho_pos`` and ``rho_neg``
    are ``f_pos * M / E_r`` and ``f_neg * M / E_r``, i.e. psi's OWN coordinates up to a common factor
    that cancels in the ratio, so ``tau`` needs neither ``kappa`` nor a belief about ``f_g``. This is
    ``bp_solver``'s own ``tth`` formula, reached here by derivation rather than by copying it.

    ⭐ **The precision is the delta method and carries no new constant.** With
    ``Var(log rho_c) = 1/p_c`` independent across the two strands,
    ``d tau / d log rho_+- = +-(1 - tau^2)/2`` and ``d theta / d tau = (1 - tau^2)^(-1/2)``, so

        ``Var(theta) = ((1 - tau^2)/4) * (1/p_+ + 1/p_-)``.

    ⚠ The arcsin map is singular at ``tau = +-1`` (all RNA on one strand): the linearisation the delta
    method rests on is invalid there, and a variance going to zero would read as INFINITE confidence in
    a boundary tilt. ``live`` is the whole guard — both densities strictly positive, so ``|tau| < 1``,
    and both precisions strictly positive, so there is something to propagate. ⛔ ONE home: a second,
    redundant test on the variance is what made this gate VACUOUS under perturbation the first time it
    was written (`TRAPS.md` A2 found it, A11 names the shape).

    ⭐ And ``1 - tau^2`` is written as ``4*rho_+*rho_-/(rho_+ + rho_-)^2`` — algebraically the same
    quantity, but strictly positive whenever both densities are, where ``1 - tau*tau`` rounds to exactly
    zero as soon as one strand dominates the other by more than the float resolution. That rounding is
    the same infinite precision arriving by a different door.
    """
    rp = np.asarray(rho_pos, dtype=np.float64)
    rn = np.asarray(rho_neg, dtype=np.float64)
    pp = np.asarray(prec_pos, dtype=np.float64)
    pn = np.asarray(prec_neg, dtype=np.float64)
    live = (rp > 0.0) & (rn > 0.0) & (pp > 0.0) & (pn > 0.0)
    total = np.where(live, rp + rn, 1.0)
    tau = np.where(live, (rp - rn) / total, 0.0)
    theta = np.where(live, np.arcsin(np.clip(tau, -1.0, 1.0)), 0.0)
    var = (rp * rn) / (total * total) * (
        1.0 / np.where(live, pp, 1.0) + 1.0 / np.where(live, pn, 1.0)
    )
    return theta, np.where(live, 1.0 / np.where(live, var, 1.0), 0.0)


def share_from_density(rho, eff, mass, prec, *, name):
    """``log`` of a density measurement read as a claim on the destination's OWN share.

    A neighbour states a density; the destination converts it with its own ``(M, E)`` and nothing else —
    no belief crosses, so this is D4-clean.

    ⛔⛔ **AN ASSERTION, NOT A CLIP** (owner, 2026-08-06). A share is in ``[0, 1]`` by definition, so a
    value above 1 is an absurdity and not a number to be repaired. Clipping it to ``1.0`` states
    **maximal gDNA** — which at a zero-gDNA library is the worst answer available — and converts a loud
    contradiction into a silent, confident error (``TRAPS.md`` A6: a ``min()`` clip hid an exact factor
    of 2 for months). A prototype must be maximally brittle: what an over-unit share means is that the
    null is inconsistent with the destination's own count, and §5 already says how to read that.

    ⚠ Only checked where the claim is CARRIED (``prec > 0``). A muted channel's mode is never read, and
    a slot with no mass would otherwise raise on a ``0/eps`` it does not deliver.
    """
    rho = np.asarray(rho, dtype=np.float64)
    eff = np.asarray(eff, dtype=np.float64)
    mass = np.asarray(mass, dtype=np.float64)
    share = rho * eff / np.maximum(mass, _EPS)
    bad = (np.asarray(prec, dtype=np.float64) > 0.0) & (share > 1.0 + 1e-9)
    if bad.any():
        i = int(np.flatnonzero(bad)[np.argmax(share[np.flatnonzero(bad)])])
        raise AssertionError(
            f"{name}: a density measurement claims share {share[i]:.6g} > 1 at slot {i} "
            f"(rho={rho[i]:.6g}, eff={eff[i]:.6g}, mass={mass[i]:.6g}) — the null is inconsistent "
            f"with the destination's own count and must be read as §5 says, never clipped to 1.0"
        )
    return np.log(np.maximum(share, _EPS))


def cross_composition(eta_src, t_src, t_dst, *, rho_g_src, mass_dst, eff_g_dst, eff_r_dst):
    """Carry a composition claim from slot ``s`` to slot ``d`` — the three regimes of the derivation.

    Returns ``(eta_dst, determined)``. ``determined`` is False only where the destination holds a
    population the message cannot speak about at all.

    * ``T_s == T_d`` — the populations match, so ``eta`` crosses as the IDENTITY. Nothing else is needed:
      the destination's own opportunities convert it and its own ``M`` supplies the level.
    * ``T_s > T_d`` — drop the strands the destination cannot hold, then as above.
    * ``T_s < T_d`` — the shared densities cross UNCHANGED (a density is already frame-free), so the mass
      they account for at the destination is ``A = sum_{c in C} rho_c(s) E_c(d)`` and the newly-active
      population takes the residual ``M_d - A``. With one newly-active population that is DETERMINED —
      one unknown, one equation. With two (both strands newly live) the RNA TOTAL is still determined and
      only the TILT is not.

    ⚠ ``A > M_d`` is possible from sampling noise. The null is then inconsistent with the destination's
    own count and the honest reading is that the new population is absent; the excess is absorbed by the
    shared claim's own uncertainty. ⛔ No shift is introduced here — every rule that adds doubt at pass-0
    was priced on the full panel and refused by the zero control.
    """
    t_src = np.asarray(t_src, dtype=bool)
    t_dst = np.asarray(t_dst, dtype=bool)
    shared = t_src & t_dst
    newly = t_dst & ~t_src
    same = np.all(shared == t_dst, axis=-1)  # nothing newly active => eta crosses as the identity

    mass_dst = np.asarray(mass_dst, dtype=np.float64)
    eff_g_dst = np.maximum(np.asarray(eff_g_dst, dtype=np.float64), _EPS)
    eff_r_dst = np.maximum(np.asarray(eff_r_dst, dtype=np.float64), _EPS)

    # the shared populations' claim on the destination's own mass (only gDNA can be shared-and-carried
    # as a bare density here; a shared RNA strand travels inside eta, which is the `same` branch).
    acc = np.where(shared[..., GDNA], np.asarray(rho_g_src, dtype=np.float64) * eff_g_dst, 0.0)
    residual = mass_dst - acc
    # the null is inconsistent with the destination's own count => the new population is absent.
    residual = np.maximum(residual, 0.0)
    rho_new = residual / eff_r_dst
    with np.errstate(divide="ignore", invalid="ignore"):
        eta_resid = np.log(np.maximum(np.asarray(rho_g_src, dtype=np.float64), _EPS)) - np.log(
            np.maximum(rho_new, _EPS)
        )
    eta_out = np.where(same, np.asarray(eta_src, dtype=np.float64), eta_resid)
    # undetermined only in the TILT, and only when both strands are newly active; the gDNA-vs-RNA split
    # the tool exists to estimate is determined in every case.
    determined = same | (newly.sum(axis=-1) <= 1)
    return eta_out, determined
