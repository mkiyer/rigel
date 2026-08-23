"""The PASS-0 FAN-OUT policy — stage 3 of the first-pass redesign: solved single-stranded introns
speak to their flanking ``intron|exon`` boundaries, and nothing else says anything.

       Gate: ``tests/calibration/test_fanout_policy.py``

The message is the intron's OWN local solve delivered on the λ channel — mode ``logit(f_g)`` clipped
into the grid domain (`TRAPS: off-grid-message-mode`), precision ``tau_lam``, the strand +
intron-factory Fisher. The reference's location term does not enter the precision (the §9c.2 ruling:
a prior's curvature is not the data's information), and the transfer carries **no drift allowance —
measured, not assumed**: `hop_currency.py`'s excess-over-floor for ``B exon|intron[sj] ← R intron``
is 0.2–0.4 % on both panels, capture ON and OFF (2026-08-23), so the source posterior's own width is
the honest price of this hop.

⭐ **Why this hop is legitimate where `intron_phi` was refuted.** The premise — the intron's contained
fragments and the boundary's unspliced crossings are ONE population — is certified by the stage-0
confusion matrix (32/32 conditions, zero violating fragments), and the currency map says COMPOSITION
is what survives this hop under capture (the LEVEL dies by ~30×). `intron_phi` transported the same
composition but as an UNCONDITIONAL overwrite with no fuse; here the claim is one precision-weighted
factor among the destination's own channels, and the intron beliefs themselves are post-
`measured_intron_reference`.

⛔ The licence is STRUCTURAL and lives on the SOURCE alone: ``ss_intron_region``
(`structural_claims`, the stage-0 substrate). ⭐ The destination end is IMPLIED, not masked — a
perturbation sweep proved a destination mask does no work, and the derivation is why: a boundary
flanking a single-stranded intron is ALWAYS ``ss_intron_boundary`` (the intron flank carries no exon
bits, so no contiguous exon crosses; and the intron's transcript continues across the boundary, so
exactly its one strand is continuous there). Likewise an empty intron transports nothing by
arithmetic — its ``tau_lam`` IS the claim's precision, and a zero-precision claim is no claim
(`TRAPS: zero-the-precision-with-the-value`, stated once as the state pair). ⚠ When stage 4 widens
the fan-out to exon destinations, the destination question returns FOR REAL and gets its own mask.
Fan-out depth is ONE by design; the boundary→exon leg is stage 4, which carries its own measured
drift (5–10 % excess under capture-ON) and lands separately.
"""

from __future__ import annotations

import numpy as np

from . import NeighbourState, PsiMessage, StepContext

__all__ = ["FanOutPolicy"]

_EPS = 1.0e-12


class _FanOutRelay:
    def __init__(self, ctx: StepContext):
        claims = ctx.claims
        L = float(ctx.logodds_window)
        f = np.clip(np.asarray(ctx.own.f_g, np.float64), _EPS, 1.0 - _EPS)
        lam = np.clip(np.log(f) - np.log1p(-f), -L, L)
        tau = np.asarray(ctx.own.tau_lam, np.float64)
        source = np.asarray(claims.ss_intron_region, bool)
        # the per-slot claim each SOURCE publishes; a non-source publishes (0, 0) — no claim, no
        # precision, one statement (`TRAPS: zero-the-precision-with-the-value`).
        self._lam = np.where(source, lam, 0.0)
        self._tau = np.where(source, tau, 0.0)

    def scan(self, *, backward: bool):
        """Depth-ONE fan-out: the claims do not chain, so the running state IS the static per-slot
        claim and ``step`` has nothing to accumulate — the backbone's source-indexed gather does the
        entire delivery."""
        state = (self._lam, self._tau)

        def step(s: int, i: int) -> None:
            return None

        def publish():
            return state

        return step, publish

    def deliver(self, left: NeighbourState, right: NeighbourState) -> PsiMessage:
        ll, lt = (np.asarray(a, np.float64) for a in left.state)
        rl, rt = (np.asarray(a, np.float64) for a in right.state)
        lt = np.where(left.valid, lt, 0.0)
        rt = np.where(right.valid, rt, 0.0)
        prec = lt + rt
        # two licensed flanks (an intron|intron boundary) are two independent measurements of one
        # shared composition: precision-weighted mean, precisions ADD.
        mode = np.where(prec > 0.0, (lt * ll + rt * rl) / np.maximum(prec, _EPS), 0.0)
        return PsiMessage(lam_mode=mode, lam_prec=prec)


class FanOutPolicy:
    """Stage 3's message policy: intron → flanking ``intron|exon`` boundary, λ channel only."""

    name = "fanout"

    def prepare(self, ctx: StepContext) -> _FanOutRelay:
        return _FanOutRelay(ctx)


def flank_to_exon_lambda(
    lam_b,
    tau_b,
    U,
    S_face,
    eff_sj_face,
    eff_gx,
    eff_rx,
    eff_gc,
    eff_rc,
    logodds_window,
):
    """Stage 4's transfer — route-merge ∘ reframe, derived JOINTLY: the flank's composed unspliced
    claim ``(lam_b, tau_b)`` plus the spliced flux entering the exon, delivered as ``(lam, tau)`` in
    the EXON's frame.

    The routes merge as DENSITIES (each count over ITS OWN opportunity — the two footings differ by a
    measured 7–10×, so raw-count pooling is inadmissible), and the ratio is re-formed through the
    exon's opportunities, so the absolute level CANCELS — a uniform capture pull moves nothing, which
    is what the currency map certified about this hop::

        f     = sigma(lam_b)                       rho_g  = f·U/E_gx        rho_nu = (1−f)·U/E_rx
        rho_mu = S/E_sj                            rho_r  = rho_nu + rho_mu
        lam   = log(rho_g·E_gc) − log(rho_r·E_rc)          clipped into [−L, +L]

    ⛔ **ONE function, not route-merge THEN reframe through a (density, variance) interface — and that
    is `TRAPS: two-gaussians-one-latent` speaking**: ``rho_g`` and ``rho_nu`` share ``lam_b`` and
    ``U``, so a split interface loses their covariance and double-counts both. Carried jointly (the
    delta method on the three RAW statistics ``lam_b``, ``log U``, ``log S`` — each entering once,
    D.1's requirement)::

        Var(lam) = [(1−f) + f·rho_nu/rho_r]²/tau_b + (rho_mu/rho_r)²·(1/U + 1/S)

    whose limits are the falsification: at ``S = 0`` the transfer is a pure reparameterization —
    ``lam = lam_b + log(E_rx·E_gc/(E_gx·E_rc))`` at EXACTLY ``1/tau_b``, the count noise cancelling in
    the ratio — and at ``S ≫`` the gDNA side's Poisson noise carries. Monte-Carlo-gated.

    A claim requires ``tau_b > 0``, live opportunities, and live densities (``U = 0`` withdraws the claim through ``rho_g``) — anything else returns ``tau = 0``,
    value and precision withdrawn in one statement. ``S_face``/``eff_sj_face`` are the FACE of the sj
    flux on the exon's side (`RegionGeometry.sj_count_lo/hi`, ``eff_sj_lo/hi`` — the population of a
    message is direction-dependent).
    """
    lam_b = np.asarray(lam_b, np.float64)
    tau_b = np.asarray(tau_b, np.float64)
    U = np.asarray(U, np.float64)
    S = np.asarray(S_face, np.float64)
    E_sj = np.asarray(eff_sj_face, np.float64)
    gx = np.asarray(eff_gx, np.float64)
    rx = np.asarray(eff_rx, np.float64)
    gc = np.asarray(eff_gc, np.float64)
    rc = np.asarray(eff_rc, np.float64)
    L = float(logodds_window)

    # ⭐ no explicit U guard: with U = 0 the gDNA density is 0 and the ``rho_g > 0`` refinement below
    #   withdraws the claim — the perturbation sweep proved a separate conjunct did no work.
    live = (tau_b > 0.0) & (gx > 0.0) & (rx > 0.0) & (gc > 0.0) & (rc > 0.0)
    f = 1.0 / (1.0 + np.exp(-np.clip(lam_b, -L, L)))
    with np.errstate(divide="ignore", invalid="ignore"):
        rho_g = np.where(live, f * U / np.maximum(gx, _EPS), 0.0)
        rho_nu = np.where(live, (1.0 - f) * U / np.maximum(rx, _EPS), 0.0)
        rho_mu = np.where((S > 0.0) & (E_sj > 0.0), S / np.maximum(E_sj, _EPS), 0.0)
        rho_r = rho_nu + rho_mu
        live = live & (rho_g > 0.0) & (rho_r > 0.0)
        lam = np.where(
            live,
            np.clip(
                np.log(np.maximum(rho_g * gc, _EPS)) - np.log(np.maximum(rho_r * rc, _EPS)),
                -L,
                L,
            ),
            0.0,
        )
        w_mu = np.where(live, rho_mu / np.maximum(rho_r, _EPS), 0.0)
        a = np.where(live, (1.0 - f) + f * rho_nu / np.maximum(rho_r, _EPS), 0.0)
        var = a * a / np.maximum(tau_b, _EPS) + w_mu * w_mu * (
            1.0 / np.maximum(U, _EPS) + np.where(S > 0.0, 1.0 / np.maximum(S, _EPS), 0.0)
        )
        tau = np.where(live & (var > 0.0), 1.0 / np.maximum(var, _EPS), 0.0)
    return lam, tau
