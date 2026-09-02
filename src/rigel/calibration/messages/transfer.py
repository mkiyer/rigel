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

**RUNG 1 COMPLETED — the EXON -> intron|exon BOUNDARY message (owner design 2026-09-02).** The
exon publishes its OWN evidence only — its strand row (`simplex_logodds.strand_row_logodds`, the
solver's own frozen-variance term, so an unstranded library's exon says nothing), and only when
the solver's derived strand DEADBAND declares the channel live (`own.tau_lam > 0`; no constant).
The exon holds one population more than the crossing — the mature RNA that arrived by the sj — so
the boundary reads the exon's row AT the splice-in map: the owner's rescale / subtract / rescale,
in which the enrichment ratio CANCELS and only the face's own spliced-to-unspliced ratio survives,
``f_b = f_E (U + S) / U``. ⛔ Both components must see ONE opportunity ratio between the two
locales for that cancellation to hold, so the map runs on the capture-blind GEOMETRIC opportunity
(``eff_gdna_global``) for gDNA and RNA alike, with the spliced density ``S / A_g^b``; a
capture-aware opportunity on one component alone re-introduced a level across locales and was
measured to fail on sparse probes. The width is the MARGINAL over the measured ratio,
``log rho ~ N(log S/U, 1/S + 1/U)`` on equal-probability nodes — wide exactly where the map is
sensitive (a thin face near the pure-gDNA vertex), no wider elsewhere; a uniform blur in ``lam``
under-states it there. The premise the message carries — spliced and unspliced fragments at the
same face share capture affinity — was MEASURED as a bias (``a ~ 1.3`` under benign capture, ~2.2
on junction probes, 1 off capture) with negligible spread; it is recorded, not corrected. Ladder:
unstranded rows byte-identical, stranded capture-ON rows 0.987–0.995x, capture-OFF within 33
fragments, the reversal fails everywhere the message acts.

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
from scipy.stats import norm

from ..splice_graph import (
    FLAG_TES_NEG as _TES_NEG,
    FLAG_TES_POS as _TES_POS,
    FLAG_TSS_NEG as _TSS_NEG,
    FLAG_TSS_POS as _TSS_POS,
)
from ..simplex_logodds import strand_row_logodds
from . import NeighbourState, PsiMessage, StepContext

__all__ = [
    "TransferPolicy",
    "edge_bound_row",
    "face_is_licensed",
    "face_map_lambda",
    "splice_out_row",
    "transport_row",
]

_EPS = 1.0e-9
_TERM = _TSS_POS | _TSS_NEG | _TES_POS | _TES_NEG
#: the marginal over ``log rho`` is taken on equal-probability nodes of the standard normal —
#: quadrature resolution, like ``n_grid``, not a model constant
_MARGINAL_NODES = norm.ppf((np.arange(9) + 0.5) / 9.0)


def edge_bound_row(lam, n_b, n_e, a_g_b, a_g_e):
    """RUNG 3 — the intergenic|exon EDGE's claim, LOWER BOUND ONLY (owner ruling 2026-09-02:
    the enrichment-ceiling upper side is refused as over-engineering, and the zero-edge
    residual at zero-gDNA libraries is an ACCEPTED error until something clearly better exists).

    The edge's crossing is structurally pure gDNA (certified 0.0 RNA on every panel), and by
    the documented tool-scope assumption — probe panels neither target intergenic sequence at
    gene boundaries nor trail past annotated ends — an exon is at-least-as-enriched as its own
    edge, so the bias direction is STRUCTURAL: the edge under-reads, never over-reads. The
    honest claim is therefore the PROFILE likelihood over the nuisance enrichment ``s >= 1``:
    ``sup_s Pois(n_b; c(lam)/s)`` with ``c = sigma(lam)·n_e·A_g^edge/E_g^exon`` — exactly 0
    wherever ``c >= n_b`` (some enrichment explains any excess) and the edge count's own
    one-sided Poisson tail ``n_b·log(c/n_b) − (c − n_b)`` below it. At ``n_b = 0`` the row is
    identically zero: a zero edge is VACUOUS, never a claim (measured: near-zero rows perturb
    flat posteriors through refit amplification, so nothing below the profile's own content may
    be delivered)."""
    lam = np.asarray(lam, np.float64)
    n_b = float(n_b)
    if not n_b > 0.0:
        return np.zeros_like(lam)
    sig = 1.0 / (1.0 + np.exp(-lam))
    c = sig * float(n_e) * float(a_g_b) / float(a_g_e)
    return np.where(c >= n_b, 0.0, n_b * np.log(np.maximum(c, 1e-300) / n_b) - (c - n_b))


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


def splice_out_row(row_e, lam, n_u, n_s, a_g_b, a_g_e):
    """The exon -> boundary claim: the exon's own row ``row_e`` read at the GEOMETRIC splice-in map
    (both components on the gDNA opportunity ``a_g``; spliced density ``n_s / a_g_b``), marginalised
    over the measured ratio ``log rho ~ N(log n_s/n_u, trigamma(n_s+1/2) + trigamma(n_u+1/2))`` — the
    same counting price rung 2 charges its ingredients; max-normalised. Vacuous
    (all zeros) at a depleted face or a flat row — silence, never a near-zero row."""
    lam = np.asarray(lam, np.float64)
    r = np.asarray(row_e, np.float64)
    if not (n_u > 0.0 and a_g_b > 0.0 and a_g_e > 0.0) or np.ptp(r) <= _EPS:
        return np.zeros_like(lam)
    r = r - r.max()
    sd = float(np.sqrt(polygamma(1, float(n_s) + 0.5) + polygamma(1, float(n_u) + 0.5)))
    acc = np.zeros_like(lam)
    for z in _MARGINAL_NODES:
        s = float(n_s) / float(a_g_b) * float(np.exp(z * sd))
        m = face_map_lambda(lam, n_u, a_g_b, a_g_b, a_g_e, a_g_e, s)
        acc += np.exp(np.interp(m, lam, r, left=r[0], right=r[-1]))
    out = np.log(np.maximum(acc / _MARGINAL_NODES.size, 1.0e-300))
    out -= out.max()
    return out if np.ptp(out) > _EPS else np.zeros_like(lam)


class TransferPolicy:
    """The composition-transfer policy. ``rows_at(n_grid, logodds_window)`` supplies the intron
    factory's per-slot rows on the sweep's own grid (``calibrate`` passes its memoized
    ``_intron_prior_at``; ``None`` means no factory evidence and the policy is silent — the
    rung-0 identity)."""

    name = "transfer"

    def __init__(self, rows_at: Callable, strand: tuple[float, float, float] | None = None):
        """``strand`` = ``(kappa, od_g, od_r)``, the library's fitted strand model, which the exon ->
        boundary message needs to state an exon's own row; ``None`` leaves that message off
        (rungs 1–3 exactly)."""
        self._rows_at = rows_at
        self._strand = None if strand is None else tuple(float(x) for x in strand)

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

        # ── the LICENSED intron|exon faces of an exon: (boundary, intron, the face's flux column) ──
        # A face is licensed when no terminus sits on it and both flanks admit the same strands (an
        # UNMEASURED population change refuses the hop). The spliced flux that enters the exon across
        # its LEFT face is the sj whose HIGH end is that boundary, and vice versa — the same column
        # serves the splice-in (rung 2) and the splice-out (item 1) directions of one face.
        lam = np.linspace(-float(ctx.logodds_window), float(ctx.logodds_window), src.shape[1])
        flags = np.asarray(ctx.boundary_flags, np.uint16)
        n_u = np.asarray(ctx.n_slot, np.float64)
        a_g = np.asarray(ctx.eff_gdna_global, np.float64)
        a_r = np.asarray(ctx.eff_rna, np.float64)
        rr = (
            np.asarray(ctx.route_rate_lo, np.float64).sum(axis=1),
            np.asarray(ctx.route_rate_hi, np.float64).sum(axis=1),
        )
        sc = (
            np.asarray(ctx.sj_count_lo, np.float64).sum(axis=1),
            np.asarray(ctx.sj_count_hi, np.float64).sum(axis=1),
        )

        def licensed_faces(e):
            for b, hi in ((left[e], 1), (right[e], 0)):
                if b < 0 or not is_bnd[b]:
                    continue
                i = left[b] if right[b] == e else right[b]
                if i < 0 or not is_intron[i]:
                    continue  # an exon|exon or terminus-adjacent face: not these rungs' hop
                if face_is_licensed(flags[b], fp[e], fn[e], fp[i], fn[i]):
                    yield b, i, hi

        # ── RUNG 2: the face-composed transfer into exons ────────────────────────────────────
        # Each licensed face contributes the intron row transported through the face map and
        # widened by the face's own counting variance; two faces sum as independent witnesses.
        for e in np.flatnonzero(is_exon):
            add = None
            for b, i, hi in licensed_faces(e):
                row_i = src[i]
                if np.ptp(row_i) <= _EPS or not (
                    n_u[b] > 0 and a_g[b] > 0 and a_r[b] > 0 and a_g[e] > 0 and a_r[e] > 0
                ):
                    continue  # no split evidence, or a depleted/empty face: silence
                le = face_map_lambda(lam, n_u[b], a_g[b], a_r[b], a_g[e], a_r[e], float(rr[hi][b]))
                r = transport_row(row_i, lam, le, n_u[b], float(sc[hi][b]))
                add = r if add is None else add + r
            # ── RUNG 3: the intergenic|exon EDGE faces, lower bound only ─────────────────
            is_intergenic = ~is_bnd & ~is_exon & ~fp & ~fn
            for b in (left[e], right[e]):
                if b < 0 or not is_bnd[b]:
                    continue
                o = left[b] if right[b] == e else right[b]
                if o < 0 or not is_intergenic[o]:
                    continue  # not a gene edge
                if not (n_u[b] > 0 and a_g[b] > 0 and a_g[e] > 0):
                    continue  # a zero edge is vacuous by the profile — silence, never ~0 rows
                r = edge_bound_row(lam, n_u[b], n_u[e], a_g[b], a_g[e])
                if np.ptp(r) <= _EPS:
                    continue
                r = r - r.max()
                add = r if add is None else add + r
            if add is not None:
                rows[e] = add - add.max()
                live = True

        # ── RUNG 1 COMPLETED: the EXON -> intron|exon BOUNDARY message (item 1) ─────────────────
        # A strand-live exon publishes its own strand row; each licensed face's boundary reads it
        # at the geometric splice-in map, marginalised over the face's measured ratio, and sums it
        # with the intron row already there — two independent witnesses of one composition.
        if self._strand is not None:
            kappa, od_g, od_r = self._strand
            tau = np.asarray(ctx.own.tau_lam, np.float64)
            cnt = np.asarray(ctx.unspliced_count, np.float64)
            belief = np.asarray(ctx.belief_fg, np.float64)
            for e in np.flatnonzero(is_exon & (fp != fn) & (tau > 0.0)):
                f_ref = float(belief[e]) if np.isfinite(belief[e]) else 0.5
                row_e = strand_row_logodds(
                    lam, cnt[e, 0], cnt[e, 1], bool(fp[e]), kappa, od_g, od_r, f_ref
                )
                if np.ptp(row_e) <= _EPS:
                    continue  # inside the deadband in all but name: silence
                for b, _i, hi in licensed_faces(e):
                    r = splice_out_row(row_e, lam, n_u[b], float(sc[hi][b]), a_g[b], a_g[e])
                    if np.ptp(r) <= _EPS:
                        continue
                    rows[b] = rows[b] + r
                    rows[b] -= rows[b].max()
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
