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
