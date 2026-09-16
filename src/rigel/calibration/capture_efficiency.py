"""The per-piece capture efficiency — the posterior mean of each piece's clipped gDNA density against
the fully captured level, under the fitted gDNA landscape, from every unspliced gDNA object whose
fragments overlap the piece.

gDNA is one template at a uniform rate before capture, so its density after capture on a piece of the
genome is that piece's capture efficiency up to the one unit ``ρ_ref`` — the located enriched mode of
the landscape (`abundance_landscape.located_enriched_mode`), the fully captured level. A piece's
efficiency is

    c̃_p = E[ min(ρ_p / ρ_ref, 1) | evidence ],

the expectation under the landscape ``P(log ρ)`` as the prior (the population's own statement of where
gDNA densities sit) and the Poisson counting rule as the likelihood. At high depth it is the plug-in
``min(k/S/ρ_ref, 1)``; at low depth it is the population's mixture weighted by the piece's own
likelihood; a piece with no evidence at all reads the population's clipped mean. No constant enters and
no floor: the multimapper floor ``C/(C+1)`` this replaces was a 30–40× bias on every unprobed
transcript (`ISSUES: ruler-multimapper-floor-caps-the-correction`), protecting against a plug-in's
exact zero that a posterior mean never produces.

THE EVIDENCE is two kinds of object. The piece's own CONTAINED count ``k_p`` on its support ``S_p``,
``E[k_p] = ρ_p S_p``. And the CROSSING counts at every boundary within a fragment's reach: the
accumulator counts a fragment at every boundary it crosses, and a crossing fragment's expected weight is
the mean efficiency over its bases, so ``E[k_e] = Σ_q A_eq ρ_q`` over the pieces ``q`` its fragments
cover, with ``A_eq`` the base-starts each piece takes (`effective_length.crossing_base_shares`). A piece
shorter than a fragment holds no contained fragment — its support and its own count are nil — and its
evidence is its edge crossings read against its neighbours' levels: an exon of 40 bp beside a 1 kb intron
whose own long count pins the intron takes the count the intron cannot explain. That is what makes a
transcript of tiny exons readable at all (the plan's worked locus; the test chromosome's tiny-exon block).

THE APPORTIONMENT. A crossing count is a Poisson sum over the pieces within reach, and the E-step of a
Poisson sum attributes it: ``z_eq = k_e · A_eq ρ̄_q / Σ_q' A_eq' ρ̄_q'`` with ``ρ̄`` the pieces'
own-count posterior means, each piece then reading its share as a count on its own exposure ``A_eq``.
One pass: the neighbours are held at what their own counts say. Iterating the attribution with the
updated means is EM on the joint and converges (13 passes on the test chromosome, 59 on the ladder, at
the grid step), and it changes nothing the truth instrument can see — identical on every row of the
test chromosome, +0.44 against +0.46 nat on the unprobed class of the ladder's `g05 ss.99 ON` row — so
the one pass ships; a joint update of neighbours instead (each conditioning on the other's mean with the
whole count) never settles, 64–66 pieces of the test chromosome and ~2,000 of the ladder flipping by 8 nat
every pass, and is refused.

A BOUNDARY has an efficiency of its own too — the posterior of its crossing count on its crossing support,
the object the calibration's crossing mass sits on. The transcript ruler never reads it (no boundary object
enters a length over bases), but the locus prior must: its COUNT is the calibration's masses on regions and
boundaries, and a length consistent with that count is those same objects at their own efficiencies.

Layer 5: a posterior under the density prior. The consumers are `calibrate`, which publishes the
efficiencies on the result, and through it `capture_eff_length` and `priors`. Gate:
``tests/calibration/test_capture_efficiency.py``.
"""

from __future__ import annotations

import numpy as np

from .landscape import DensityLandscape
from .simplex_logodds import _block_rows

__all__ = ["capture_efficiencies"]

_EPS = 1e-300


def _posterior(
    landscape: DensityLandscape,
    rho_ref: float,
    k_own: np.ndarray,
    S_own: np.ndarray,
    z_cross: np.ndarray,
    A_cross: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Per piece ``(E[min(ρ/ρ_ref, 1)], E[ρ])`` on the landscape's grid: the prior times the Poisson
    likelihood of the own count on ``S_own`` and of the attributed crossing count ``z_cross`` on the
    pooled exposure ``A_cross`` (the terms of one piece pool exactly into one Poisson term). Row tiles,
    as ψ tiles (`simplex_logodds._block_rows`): a million pieces never exist against the grid at once."""
    rho = np.exp(np.asarray(landscape.log_rho, dtype=np.float64))
    lp = np.asarray(landscape.logP, dtype=np.float64)
    clipped = np.minimum(rho / rho_ref, 1.0)
    n = k_own.shape[0]
    c_tilde = np.empty(n)
    e_rho = np.empty(n)
    rows = _block_rows(rho.size, 8)
    for r0 in range(0, n, rows):
        sl = slice(r0, min(r0 + rows, n))
        ll = np.zeros((sl.stop - sl.start, rho.size))
        for k, S in ((k_own[sl], S_own[sl]), (z_cross[sl], A_cross[sl])):
            lam = rho[None, :] * S[:, None]
            live = S[:, None] > 0.0
            ll += np.where(live, k[:, None] * np.log(np.maximum(lam, _EPS)) - lam, 0.0)
        lw = ll + lp[None, :]
        lw -= lw.max(1, keepdims=True)
        w = np.exp(lw)
        w /= w.sum(1, keepdims=True)
        c_tilde[sl] = w @ clipped
        e_rho[sl] = w @ rho
    return c_tilde, e_rho


def capture_efficiencies(
    landscape: DensityLandscape,
    rho_ref: float,
    count_region: np.ndarray,
    eff_region: np.ndarray,
    count_boundary: np.ndarray,
    eff_boundary: np.ndarray,
    shares: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """``(efficiency_region, efficiency_boundary)``: every piece's ``c̃_p`` from its own contained count
    and the crossings within reach apportioned to it, and every boundary's ``c̃_e`` from its own crossing
    count on its own support — the efficiency of the object the calibration's crossing mass sits on,
    which is what a locus prior's length needs beside its count (`priors.assemble_priors`).

    ``count_region`` / ``count_boundary`` are the deconvolved gDNA masses (a mass, so the Poisson enters
    through the gamma function as a continuation), ``eff_region`` / ``eff_boundary`` the contained and
    crossing supports, ``shares`` the ``(boundary, region, share)`` triples of
    `effective_length.crossing_base_shares`. The crossings are attributed with the own-count posterior
    means, once.
    """
    k_reg = np.maximum(np.asarray(count_region, dtype=np.float64), 0.0)
    S_reg = np.maximum(np.asarray(eff_region, dtype=np.float64), 0.0)
    k_bnd = np.maximum(np.asarray(count_boundary, dtype=np.float64), 0.0)
    S_bnd = np.maximum(np.asarray(eff_boundary, dtype=np.float64), 0.0)
    E, Q, A = shares
    n = k_reg.shape[0]

    zero = np.zeros(n)
    _c_own, rho_own = _posterior(landscape, rho_ref, k_reg, S_reg, zero, zero)
    total = np.zeros(k_bnd.shape[0])
    np.add.at(total, E, A * rho_own[Q])
    z = k_bnd[E] * A * rho_own[Q] / np.maximum(total[E], _EPS)
    z_pool = np.zeros(n)
    A_pool = np.zeros(n)
    np.add.at(z_pool, Q, z)
    np.add.at(A_pool, Q, A)
    c_region, _rho = _posterior(landscape, rho_ref, k_reg, S_reg, z_pool, A_pool)
    zb = np.zeros(k_bnd.shape[0])
    c_boundary, _rho_b = _posterior(landscape, rho_ref, k_bnd, S_bnd, zb, zb)
    return c_region, c_boundary
