"""TransferPolicy — the COMPOSITION TRANSFER policy (rungs 1 and 2 of the ground-up message
rebuild; owner rulings 2026-09-01).

       Gate: ``tests/calibration/test_transfer_policy.py``

**RUNG 1 — the intron|exon BOUNDARY.** An intron REGION and its intron|exon BOUNDARY share
their unspliced population set — mature RNA cannot cross that boundary contiguously, so both
objects hold {gDNA, unspliced RNA} of the same transcripts — and a composition is scale-free,
so it crosses the hop with no reframe, no flank pair and no mass budget. The policy delivers,
at each such boundary, its intron flank's OWN density-deconvolution factor (the intron
factory's NegBinom row — the shipped likelihood, whose ``alpha_eff`` already prices the fitted
overdispersion and the background posterior's width) as ``PsiMessage.lam_rows``, VERBATIM.

**RUNG 2 — the EXON, by the FACE-COMPOSED transfer.** At a LICENSED face (no transcript
terminus, the same strand set on both sides — an UNMEASURED population change refuses the hop;
the measured one, the spliced route, joins the message instead), the exon's composition claim
is the intron row transported through the MONOTONE face map

    lam_e(lam_u) = log(n_u*sigma(lam_u)*E_g^e/A_g) − log((n_u*(1−sigma(lam_u))/A_r + s)*E_r^e)

(``n_u`` the face's unspliced crossing count, ``A``/``E`` the face's and exon's opportunities,
``s`` the face's route-summed certified spliced density) and widened by the face's OWN counting
variance ``trigamma(n_u+1/2) + trigamma(n_s+1/2)`` — the delta-method price of the map's
measured position, per face, no constants. Every ratio is formed WITHIN one locale, so no level
ever crosses a capture cliff (the anchor's refuted assumption), and the measured ``s``
structurally CAPS the exon's claimable gDNA share (the map saturates). A depleted face
(``n_u = 0``), a flat source row, or an unlicensed face delivers nothing. Everything else is
silent, and ``scan`` relays nothing — ONE hop is structural, not a discipline. ⚠ Two factors in
tension can otherwise ANNIHILATE into a flat posterior whose estimator wanders (measured: a
13-fragment face's unwidened cliff sat on the prior's 8-nat truth peak and the point estimate
landed at the centroid of a vacuous plateau) — the width is what prevents it.

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

from scipy.special import polygamma

from ..splice_graph import (
    FLAG_TES_NEG as _TES_NEG,
    FLAG_TES_POS as _TES_POS,
    FLAG_TSS_NEG as _TSS_NEG,
    FLAG_TSS_POS as _TSS_POS,
)
from . import NeighbourState, PsiMessage, StepContext

__all__ = ["TransferPolicy", "face_is_licensed", "face_map_lambda", "transport_row"]

_EPS = 1.0e-9
_TERM = _TSS_POS | _TSS_NEG | _TES_POS | _TES_NEG


def face_is_licensed(flags_b, fp_e, fn_e, fp_i, fn_i) -> bool:
    """The rung-2 face licence: NO transcript terminus at the boundary and the SAME strand set
    on both sides. An UNMEASURED population change (a terminus admits uncounted molecules; a
    strand flip changes membership) refuses the hop; the measured spliced route does not — it
    joins the message instead."""
    return (not (int(flags_b) & _TERM)) and bool(fp_e) == bool(fp_i) and bool(fn_e) == bool(fn_i)


def face_map_lambda(lam, n_u, a_g_b, a_r_b, e_g_e, e_r_e, s):
    """The face map ``lam_e(lam_u)`` on the grid — monotone nondecreasing, saturating at the
    certified-flux ceiling ``log(n_u*E_g^e/A_g) − log(s*E_r^e)`` when ``s > 0`` (a measured
    spliced density CAPS the claimable gDNA share), and a pure opportunity shift at ``s = 0``."""
    lam = np.asarray(lam, np.float64)
    sig = 1.0 / (1.0 + np.exp(-lam))
    g_arm = float(n_u) * sig / float(a_g_b) * float(e_g_e)
    r_arm = (float(n_u) * (1.0 - sig) / float(a_r_b) + float(s)) * float(e_r_e)
    return np.log(np.maximum(g_arm, 1e-300)) - np.log(np.maximum(r_arm, 1e-300))


def transport_row(row, lam, lam_e_of_u, n_u, n_s):
    """One face's transported row: the source row read at the map's preimage (a likelihood
    evaluated at the corresponding parameter point — NO Jacobian), then widened by the face's
    counting variance ``trigamma(n_u+1/2) + trigamma(n_s+1/2)`` and max-renormalized. Beyond
    the map's ceiling the row takes the flat limit (a SOFT cap — conservative)."""
    lam = np.asarray(lam, np.float64)
    r = np.asarray(row, np.float64)
    lam_u_of_x = np.interp(lam, np.asarray(lam_e_of_u, np.float64), lam, left=lam[0], right=lam[-1])
    out = np.interp(lam_u_of_x, lam, r - r.max())
    v = float(polygamma(1, float(n_u) + 0.5) + polygamma(1, float(n_s) + 0.5))
    if v > 0.0 and lam.shape[0] > 1:
        dlam = float(lam[1] - lam[0])
        half = max(int(np.ceil(4.0 * float(np.sqrt(v)) / dlam)), 1)
        x = np.arange(-half, half + 1) * dlam
        kern = np.exp(-0.5 * x * x / v)
        kern /= kern.sum()
        pr = np.convolve(np.pad(np.exp(out - out.max()), half, mode="edge"), kern, mode="valid")
        out = np.log(np.maximum(pr, 1.0e-300))
    return out - out.max()


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

        # ── RUNG 2: the face-composed transfer into exons ────────────────────────────────────
        # For each exon, each LICENSED intron|exon face contributes the intron row transported
        # through the face map and widened by the face's own counting variance; two licensed
        # faces sum as independent witnesses (left face reads the left intron, right the right).
        lam = np.linspace(-float(ctx.logodds_window), float(ctx.logodds_window), src.shape[1])
        flags = np.asarray(ctx.boundary_flags, np.uint16)
        n_u = np.asarray(ctx.n_slot, np.float64)
        a_g = np.asarray(ctx.eff_gdna_global, np.float64)
        a_r = np.asarray(ctx.eff_rna, np.float64)
        rr_lo = np.asarray(ctx.route_rate_lo, np.float64).sum(axis=1)
        rr_hi = np.asarray(ctx.route_rate_hi, np.float64).sum(axis=1)
        sc_lo = np.asarray(ctx.sj_count_lo, np.float64).sum(axis=1)
        sc_hi = np.asarray(ctx.sj_count_hi, np.float64).sum(axis=1)
        for e in np.flatnonzero(is_exon):
            add = None
            for b, entering_hi in ((left[e], True), (right[e], False)):
                if b < 0 or not is_bnd[b]:
                    continue
                i = left[b] if right[b] == e else right[b]
                if i < 0 or not is_intron[i]:
                    continue  # an exon|exon or terminus-adjacent face: not this rung's hop
                if not face_is_licensed(flags[b], fp[e], fn[e], fp[i], fn[i]):
                    continue  # an UNMEASURED population change refuses the hop
                row_i = src[i]
                if np.ptp(row_i) <= _EPS or not (
                    n_u[b] > 0 and a_g[b] > 0 and a_r[b] > 0 and a_g[e] > 0 and a_r[e] > 0
                ):
                    continue  # no split evidence, or a depleted/empty face: silence
                s = float((rr_hi if entering_hi else rr_lo)[b])
                n_s = float((sc_hi if entering_hi else sc_lo)[b])
                le = face_map_lambda(lam, n_u[b], a_g[b], a_r[b], a_g[e], a_r[e], s)
                r = transport_row(row_i, lam, le, n_u[b], n_s)
                add = r if add is None else add + r
            if add is not None:
                rows[e] = add - add.max()
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
