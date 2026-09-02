"""TransferPolicy — the intron -> intron|exon boundary COMPOSITION TRANSFER (rung 1 of the
ground-up message rebuild; owner ruling 2026-09-01).

       Gate: ``tests/calibration/test_transfer_policy.py``

**The one claim this policy makes, and its whole scope.** An intron REGION and its intron|exon
BOUNDARY share their unspliced population set — mature RNA cannot cross that boundary
contiguously, so both objects hold {gDNA, unspliced RNA} of the same transcripts — and a
composition is scale-free, so it crosses the hop with no reframe, no flank pair and no mass
budget. The policy therefore delivers, at each such boundary, its intron flank's OWN
density-deconvolution factor (the intron factory's NegBinom row — the shipped likelihood, whose
``alpha_eff`` already prices the fitted overdispersion and the background posterior's width)
as ``PsiMessage.lam_rows``, VERBATIM. Every other channel and every other slot is silent, the
exon-side message does not exist here, and ``scan`` relays nothing — ONE hop is structural, not
a discipline.

**The hop cost is MEASURED, priced, and zero beyond what the row already carries** — an
imputation must cost something every hop, and this one's cost was derived rather than waived:
the certified intron-vs-boundary composition dispersion is 0.000 off capture (32,534 ladder
pairs, counting noise subtracted) and 0.200 under capture, and an A/B of the row blurred by
that dispersion against the unblurred row moved every ladder condition by under half a percent
(the unblurred row was slightly BETTER at the largest rows), because ``alpha_eff``'s width
already dominates it. So the delivered row is the source's factor unmodified, and no dispersion
constant ships; if the pair dispersion is ever re-priced as material, it re-enters as a
runtime-fitted widening, never as a constant.

**Why a row factor and not a Gaussian.** Near the pure-gDNA vertex the intron's likelihood is
one-sided — flat toward +lambda with a cliff below — and a symmetric ``(mode, precision)`` pair
cannot carry a cliff: its honest second moment reads the flat top and collapses. The row factor
is psi's general evidence currency and keeps the cliff.

**The laws it honours by construction**: the sender publishes its claim unchanged (the row is
the source's own factor, never a tailoring to the destination); the recipient decides (psi
fuses the rows against the slot's own evidence in the FINAL solve only — the certified-flux
citizenship); a no-claim stays a no-claim (a flat source row, an absent factory, an
evidence-free provider all deliver silence, never a zero-filled channel).
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from . import NeighbourState, PsiMessage, StepContext

__all__ = ["TransferPolicy"]

_EPS = 1.0e-9


class TransferPolicy:
    """The composition-transfer policy. ``rows_at(n_grid, logodds_window)`` supplies the intron
    factory's per-slot rows on the sweep's own grid (``calibrate`` passes its memoized
    ``_intron_prior_at``; ``None`` means no factory evidence and the policy is silent — the
    rung-0 identity)."""

    name = "transfer"

    def __init__(self, rows_at: Callable):
        self._rows_at = rows_at

    def prepare(self, ctx: StepContext) -> "_PreparedTransfer":
        src = self._rows_at(int(ctx.n_grid), float(ctx.logodds_window))
        n = ctx.n_slots
        if src is None:
            return _PreparedTransfer(None)
        src = np.asarray(src, np.float64)
        # ── the shared-population pairs, derived from the annotation bits every sweep ────────
        # A BOUNDARY with an exon REGION on one flank and an RNA-admitting non-exon REGION
        # (an intron) on the other. Where a region is exonic for ANY transcript its signature
        # is exon, so an overlapping-isoform locus refuses the pair structurally.
        is_bnd = np.asarray(ctx.is_boundary, bool)
        is_exon = np.asarray(ctx.is_exon_region, bool)
        fp = np.asarray(ctx.free_pos, bool)
        fn = np.asarray(ctx.free_neg, bool)
        is_intron = ~is_bnd & ~is_exon & (fp | fn)
        left = np.asarray(ctx.left, np.int64)
        right = np.asarray(ctx.right, np.int64)
        li = np.clip(left, 0, n - 1)
        ri = np.clip(right, 0, n - 1)
        ok = is_bnd & (left >= 0) & (right >= 0)
        take_right = ok & is_exon[li] & is_intron[ri]
        take_left = ok & is_exon[ri] & is_intron[li] & ~take_right
        dst = np.flatnonzero(take_right | take_left)
        if dst.size == 0:
            return _PreparedTransfer(None)
        srcslot = np.where(take_right, ri, li)[dst]

        rows = np.zeros((n, src.shape[1]))
        live = False
        for b, j in zip(dst.tolist(), srcslot.tolist()):
            row = src[j]
            if np.ptp(row) <= _EPS:
                continue  # a flat factory row carries no claim
            rows[b] = row - row.max()
            live = True
        return _PreparedTransfer(rows if live else None)


class _PreparedTransfer:
    """One sweep's delivery: the blurred rows, or silence. ``scan`` relays nothing — the
    transfer is strictly one-hop by construction, never by discipline."""

    def __init__(self, rows: np.ndarray | None):
        self._rows = rows

    def scan(self, *, backward: bool):
        return None

    def deliver(self, left: NeighbourState, right: NeighbourState) -> PsiMessage:
        if self._rows is None:
            return PsiMessage.silent()
        return PsiMessage(lam_rows=self._rows)
