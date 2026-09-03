"""TransferPolicy — the COMPOSITION TRANSFER policy: one hop, nothing relayed, silence-not-zeros.

       Gate: ``tests/calibration/test_transfer_policy.py``

Every message is a composition — a scale-free statement of gDNA versus RNA — carried across ONE
face by the row constructors of `transfer_rows`, so no level ever crosses a capture cliff and no
constant is anywhere. Ten messages ship, each delivered as ``PsiMessage.lam_rows``:

* **intron → intron|exon boundary** (rung 1): the intron's own factory row, VERBATIM — the two
  objects share their unspliced population (mature RNA cannot cross that boundary), and the hop
  cost was measured at zero beyond the row's own width.
* **boundary → exon** (rung 2): the intron row through the splice-in face map at every LICENSED
  face (no terminus, the same strand set), widened by the face's counting variance; two faces sum
  as independent witnesses.
* **intergenic|exon edge → exon** (rung 3): the sign-certified LOWER bound; a zero edge is vacuous.
* **exon → intron|exon boundary** (rung 1 completed): the exon's OWN strand row, only where the
  solver's derived strand deadband declares the channel live, read at the splice-in map backwards
  on capture-blind opportunities and marginalised over the face's measured ratio.
* **intron|exon boundary → intron** (rung 2 completed): the boundary's OWN strand row, VERBATIM —
  the two objects share one unspliced population exactly, so there is no splice-out term and no
  premise, and the s = 0 opportunity map is the identity under the one-opportunity rule; gated by
  the same derived deadband at the boundary, where the solver's own λ-precision IS the strand
  Fisher information. Under capture a probed exon's faces carry many times the unprobed intron's
  own information; off capture they add a few percent.
* **the exon|exon TERMINUS boundary and its OUTSIDE exon** (rung 4, item 5): three messages. The
  terminating transcripts cover one flank; the OUTSIDE flank (read off the flag alone) holds exactly
  the population that crosses the boundary — counting the SPLICED crossing too, since a mature
  molecule that crosses contiguously and splices elsewhere lands in the spliced bank while a fragment
  contained in the outside piece has no junction by geometry. So the licence is item 1's map with the
  spliced crossing as S, ``f_b = f_O (U_b + S_b) / U_b`` (certified on ladder truth): the outside
  exon's own strand row and THE COMPOSED TRANSPORT — what rungs 2–3 delivered to the outside exon,
  carried one hop further, from a snapshot so the boundary never hears its own row back — reach the
  boundary through the splice-out map, and the boundary's own strand row reaches the outside exon
  through the face map with the spliced density. An empty outside piece sends and receives nothing;
  the inside flank is never a composition destination.
* **the INSIDE flank of a terminus** (rung 4, item 6 — the owner's ABUNDANCE-DISCREPANCY rule): where
  composition cannot cross, the two objects' TOTAL abundances differ by a measured ratio `r`, and two
  hypotheses explain it — enrichment (the gDNA abundance scales by `r`: composition transfers) or new
  RNA (the gDNA abundance is unchanged: the composition shifts by `r`) — so the inside's gDNA share is
  `f_c · s / r` with the step `s` between them and never above `f_c`. The boundary's OWN strand row
  travels that map; the step's spread is the hop's premise, FITTED from the served pairs' two
  witnesses beyond counting, zero where they agree — the new-RNA point exactly when the library shows
  no enrichment structure. The precision is set by the discrepancy: first do no harm.
* **the ALTERNATIVE SPLICE SITE** (rung 4, item 7): an exon|exon boundary carrying a junction and no
  terminus. The flank on the junction's intron side shares the boundary's full unspliced crossing
  (item 5's law with the spliced crossing); the exon-of-both flank holds the crossing plus the
  isoform that splices out here, MEASURED at the face as the route flux (item 1's law with the
  spliced crossing plus the flux) — both certified on the ladder at every flank length, the flanks
  swapped opening gaps that scale with the leaving isoform's share. Each flank's own strand row
  travels to the boundary and the boundary's own row to each flank, every message carrying THE
  OWNER'S DISCREPANCY RULE PER PAIR: where the pair's two witnesses — the boundary's own strand mode
  and the flank's mapped to it — disagree beyond counting, this pair's messages are widened by that
  excess; nothing is pooled across pairs and no mode is shifted (under capture the gDNA landscape
  tapers at a probed exon's edge by a locus-dependent amount, priced where it is seen).
  Deadband-gated both ways; unstranded libraries send nothing.

The laws the policy keeps: the sender publishes its claim unchanged; the recipient (psi, in the
final solve) fuses the rows against the slot's own evidence; a no-claim stays a no-claim — a flat
row, an absent factory, an evidence-free provider all deliver SILENCE, never a zero-filled channel;
``scan`` relays nothing, so one hop is structural rather than a discipline.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from ..simplex_logodds import strand_row_logodds
from . import NeighbourState, PsiMessage, StepContext
from .transfer_rows import (
    EPS,
    abundance_row,
    blur_row,
    boundary_shares_strand,
    edge_bound_row,
    face_is_licensed,
    face_map_lambda,
    junction_flanks,
    outside_flank,
    splice_out_row,
    transport_row,
)

__all__ = ["TransferPolicy"]


class TransferPolicy:
    """``rows_at(n_grid, logodds_window)`` supplies the intron factory's per-slot rows on the sweep's
    own grid (``calibrate`` passes its memoized ``_intron_prior_at``; ``None`` means no factory
    evidence and the policy is silent — the rung-0 identity). ``strand = (kappa, od_g, od_r)`` is the
    library's fitted strand model, which the exon → boundary and boundary → intron messages need to
    state a slot's own row; ``None`` leaves both off."""

    name = "transfer"

    def __init__(self, rows_at: Callable, strand: tuple[float, float, float] | None = None):
        self._rows_at = rows_at
        self._strand = None if strand is None else tuple(float(x) for x in strand)

    def prepare(self, ctx: StepContext) -> "_PreparedTransfer":
        src = self._rows_at(int(ctx.n_grid), float(ctx.logodds_window))
        if src is None:
            return _PreparedTransfer(None)
        src = np.asarray(src, np.float64)
        n = ctx.n_slots
        # ── the chain's structure, from the annotation bits, every sweep ─────────────────────
        is_bnd = np.asarray(ctx.is_boundary, bool)
        is_exon = np.asarray(ctx.is_exon_region, bool)
        fp = np.asarray(ctx.free_pos, bool)
        fn = np.asarray(ctx.free_neg, bool)
        is_intron = ~is_bnd & ~is_exon & (fp | fn)
        is_intergenic = ~is_bnd & ~is_exon & ~fp & ~fn
        left = np.asarray(ctx.left, np.int64)
        right = np.asarray(ctx.right, np.int64)
        lam = np.linspace(-float(ctx.logodds_window), float(ctx.logodds_window), src.shape[1])
        flags = np.asarray(ctx.boundary_flags, np.uint16)
        n_u = np.asarray(ctx.n_slot, np.float64)
        n_s = np.asarray(ctx.spliced_slot, np.float64)
        flux = np.asarray(ctx.sj_count, np.float64).sum(axis=1)
        a_g = np.asarray(ctx.eff_gdna_global, np.float64)
        a_r = np.asarray(ctx.eff_rna, np.float64)
        tau = np.asarray(ctx.own.tau_lam, np.float64)
        # the flux that crosses an exon's LEFT face belongs to the sj whose HIGH end is that boundary
        # and vice versa; one column serves both directions of a face
        rr = (
            np.asarray(ctx.route_rate_lo, np.float64).sum(axis=1),
            np.asarray(ctx.route_rate_hi, np.float64).sum(axis=1),
        )
        sc = (
            np.asarray(ctx.sj_count_lo, np.float64).sum(axis=1),
            np.asarray(ctx.sj_count_hi, np.float64).sum(axis=1),
        )
        rows = np.zeros((n, src.shape[1]))
        live = False

        def deliver(slot, row):
            nonlocal live
            if np.ptp(row) > EPS:
                rows[slot] += row - row.max()
                rows[slot] -= rows[slot].max()
                live = True

        def other_flank(b, e):
            return left[b] if right[b] == e else right[b]

        def licensed_faces(e):
            """(boundary, intron, flux column) for each licensed intron|exon face of exon ``e``."""
            for b, hi in ((left[e], 1), (right[e], 0)):
                if b < 0 or not is_bnd[b]:
                    continue
                i = other_flank(b, e)
                if (
                    i >= 0
                    and is_intron[i]
                    and face_is_licensed(flags[b], fp[e], fn[e], fp[i], fn[i])
                ):
                    yield b, i, hi

        # ── rung 1: the intron's factory row, verbatim, at every intron|exon boundary ─────────
        for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
            lo, hi_ = left[b], right[b]
            if is_exon[lo] and is_intron[hi_]:
                deliver(b, src[hi_])
            elif is_exon[hi_] and is_intron[lo]:
                deliver(b, src[lo])

        strand_live = np.zeros(n, bool)
        if self._strand is not None:
            strand_live = is_exon & (fp != fn) & (tau > 0.0)
            cnt = np.asarray(ctx.unspliced_count, np.float64)
            belief = np.asarray(ctx.belief_fg, np.float64)
            kappa, od_g, od_r = self._strand

            def own_row(x):
                """A slot's OWN strand row: its unspliced count split by genome strand, the variance
                frozen at its incoming belief (a source-side read); the mode is data only."""
                f_ref = float(belief[x]) if np.isfinite(belief[x]) else 0.5
                return strand_row_logodds(
                    lam, cnt[x, 0], cnt[x, 1], bool(fp[x]), kappa, od_g, od_r, f_ref
                )

        for e in np.flatnonzero(is_exon):
            # ── rung 2: each licensed face's intron row, through the face map, into the exon ──
            for b, i, hi in licensed_faces(e):
                if np.ptp(src[i]) <= EPS or not (
                    n_u[b] > 0 and a_g[b] > 0 and a_r[b] > 0 and a_g[e] > 0 and a_r[e] > 0
                ):
                    continue  # no split evidence, or a depleted face: silence
                le = face_map_lambda(lam, n_u[b], a_g[b], a_r[b], a_g[e], a_r[e], float(rr[hi][b]))
                deliver(e, transport_row(src[i], lam, le, n_u[b], float(sc[hi][b])))
            # ── rung 3: each intergenic|exon edge's lower bound, into the exon ──────────────
            for b in (left[e], right[e]):
                if b < 0 or not is_bnd[b]:
                    continue
                o = other_flank(b, e)
                if o >= 0 and is_intergenic[o] and n_u[b] > 0 and a_g[b] > 0 and a_g[e] > 0:
                    deliver(e, edge_bound_row(lam, n_u[b], n_u[e], a_g[b], a_g[e]))
            # ── rung 1 completed: the exon's own strand row, out to each licensed face ─────
            if strand_live[e]:
                row_e = own_row(e)
                for b, _i, hi in licensed_faces(e):
                    deliver(b, splice_out_row(row_e, lam, n_u[b], float(sc[hi][b]), a_g[b], a_g[e]))
        # ── rung 2 completed: each intron|exon boundary's own strand row, verbatim, into its intron ──
        # The intron and its boundary share ONE unspliced population exactly (mature RNA crosses
        # neither boundary), so the composition transfers with no splice-out term and no premise, and
        # the s = 0 opportunity map is the identity under the one-opportunity rule (the fitted fl gap
        # on the panel is ±0.02 nats — a refinement recorded for BOTH directions of this hop). The
        # row is the boundary's OWN evidence — its crossing count split by genome strand, the
        # variance frozen at the boundary's incoming belief (a source-side read; the mode is data
        # only) — never its belief, which already holds rung 1's intron row. With no factory at a
        # boundary, ``own.tau_lam`` there IS the strand Fisher information, so the deadband gate is
        # the derived one and an unstranded library sends nothing, structurally.
        if self._strand is not None:
            for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
                lo, hi_ = left[b], right[b]
                if is_exon[lo] and is_intron[hi_]:
                    i = hi_
                elif is_exon[hi_] and is_intron[lo]:
                    i = lo
                else:
                    continue
                if not (boundary_shares_strand(fp[b], fn[b], fp[i], fn[i]) and tau[b] > 0.0):
                    continue
                deliver(i, own_row(b))
        # ── rung 4, item 5: the exon|exon TERMINUS boundary and its OUTSIDE exon ───────────────
        # The terminating transcripts cover exactly one flank; the OUTSIDE flank holds exactly the
        # population that crosses the boundary — counting BOTH crossings, the unspliced U_b and the
        # SPLICED S_b (a mature molecule that crosses contiguously and splices elsewhere lands in the
        # spliced bank, while a fragment contained in the outside piece has no junction by geometry).
        # So the licence is item 1's map with S_b for the splice-out flux, f_b = f_O (U_b + S_b) / U_b,
        # certified on ladder truth (the residual +0.004 after the map, +0.055 before). Every message
        # here is a shipped constructor: the outside exon's own strand row to the boundary through
        # `splice_out_row`; the boundary's own strand row to the outside exon through the face map
        # with the spliced density; and THE COMPOSED TRANSPORT — what rungs 2–3 delivered to the
        # outside exon, carried one hop further through the same map. The no-echo law is
        # STRUCTURAL: this block reads ``rows`` as the earlier messages left them and accumulates
        # its own deliveries apart, so a boundary can never hear its own row back through the exon.
        # An empty outside piece (zero contained opportunity) sends nothing and is sent nothing; the
        # inside flank is never a destination.
        rows5 = np.zeros_like(rows)
        served5 = np.zeros(n, bool)

        def deliver5(slot, row):
            if np.ptp(row) > EPS:
                rows5[slot] += row - row.max()
                served5[slot] = True

        for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
            lo, hi_ = left[b], right[b]
            if not (is_exon[lo] and is_exon[hi_]):
                continue
            o, _i = outside_flank(flags[b], lo, hi_)
            if o is None or not boundary_shares_strand(fp[b], fn[b], fp[o], fn[o]):
                continue
            if not (n_u[b] > 0 and a_g[b] > 0 and a_g[o] > 0):
                continue
            if np.ptp(rows[o]) > EPS:
                deliver5(b, splice_out_row(rows[o], lam, n_u[b], n_s[b], a_g[b], a_g[o]))
            if self._strand is None:
                continue
            if strand_live[o]:
                deliver5(b, splice_out_row(own_row(o), lam, n_u[b], n_s[b], a_g[b], a_g[o]))
            if tau[b] > 0.0:
                le = face_map_lambda(lam, n_u[b], a_g[b], a_g[b], a_g[o], a_g[o], n_s[b] / a_g[b])
                deliver5(o, transport_row(own_row(b), lam, le, n_u[b], n_s[b]))
        for slot in np.flatnonzero(served5):
            deliver(slot, rows5[slot])
        if self._strand is None:
            return _PreparedTransfer(rows if live else None)
        # ── rung 4, item 6: THE ABUNDANCE-DISCREPANCY message into the INSIDE flank of a terminus ──
        # (owner design, 2026-09-02). Where composition cannot cross, the two objects' TOTAL
        # abundances differ by a measured ratio r, and two hypotheses explain it — enrichment (the gDNA
        # abundance scales by r: composition transfers) or new RNA (the gDNA abundance is unchanged:
        # the composition shifts by r) — so the inside's gDNA share is f_c · s / r with the step s
        # between them and never above f_c (new RNA can only add). The boundary's OWN strand row
        # travels the map; the step's spread is the hop's premise, FITTED from the served pairs' two
        # witnesses (the boundary's and the inside exon's strand modes, data only) by the method of
        # moments beyond counting — zero where they agree within counting, so the message is the
        # new-RNA point exactly when the library shows no enrichment structure. Delivered only where
        # the boundary's strand channel is live; the reverse direction and the forwarded multi-hop
        # arrivals are owed to the scan phase with their own premises.
        kappa = self._strand[0]
        ks_of = lambda x: kappa if fp[x] else 1.0 - kappa  # noqa: E731
        served = []
        for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
            lo, hi_ = left[b], right[b]
            if not (is_exon[lo] and is_exon[hi_]):
                continue
            _o, i = outside_flank(flags[b], lo, hi_)
            if i is None or not boundary_shares_strand(fp[b], fn[b], fp[i], fn[i]):
                continue
            if not (tau[b] > 0.0 and n_u[b] > 0 and n_u[i] > 0 and a_g[b] > 0 and a_g[i] > 0):
                continue
            served.append((b, i, (n_u[i] / a_g[i]) / ((n_u[b] + n_s[b]) / a_g[b])))
        v_step = 0.0
        if served:
            log_s, var_s = [], []
            for b, i, r in served:
                if not tau[i] > 0.0:
                    continue  # one witness only: nothing to fit from this pair
                ok = True
                modes = []
                for x in (b, i):
                    n = cnt[x].sum()
                    p = cnt[x, 0] / n
                    f = (p - ks_of(x)) / (0.5 - ks_of(x))
                    if not 0.0 < f < 1.0:
                        ok = False  # a vertex mode has no logarithm
                        break
                    modes.append((f, p * (1.0 - p) / n / (p - ks_of(x)) ** 2))
                if not ok:
                    continue
                (f_b, v_b), (f_i, v_i) = modes
                log_s.append(np.log(f_i / (f_b * n_u[b] / (n_u[b] + n_s[b]) / r)))
                var_s.append(v_b + v_i + 1.0 / n_u[i] + 1.0 / (n_u[b] + n_s[b]))
            if len(log_s) >= 2:
                ls, vs = np.asarray(log_s), np.asarray(var_s)
                w = 1.0 / vs
                w /= w.sum()
                mu = float(w @ ls)
                v_step = max(0.0, float(w @ (ls - mu) ** 2) - float(w @ vs))
        for b, i, r in served:
            deliver(i, abundance_row(own_row(b), lam, n_u[b], n_s[b], r, n_u[i], v_step))
        # ── rung 4, item 7: THE ALTERNATIVE SPLICE SITE — an exon|exon boundary with a junction ──
        # The flank on the junction's intron side, C, shares the boundary's full unspliced crossing
        # (item 5's law with the spliced crossing S_b); the exon-of-both flank, E, holds the crossing
        # plus the isoform that splices out here, MEASURED at the face as the route flux F (item 1's
        # law with S_b + F). Certified on the ladder at every flank length, the flanks swapped opening
        # 0.10–0.26 gaps in f. Each flank's own strand row travels to the boundary through
        # `splice_out_row`, the boundary's own row to each flank through the face map with the
        # matching spliced density — and every message carries THE OWNER'S DISCREPANCY RULE, PER PAIR:
        # where the pair's two witnesses (the boundary's own strand mode and the flank's mapped to it)
        # disagree beyond counting, this pair's messages are widened by that excess — nothing pooled
        # across pairs, no shift of any mode (under capture the gDNA landscape tapers at a probed
        # exon's edge by a locus-dependent amount; the disagreement it makes is priced where it is
        # seen). Both directions, deadband-gated.
        served7 = []
        for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
            lo, hi_ = left[b], right[b]
            if not (is_exon[lo] and is_exon[hi_]):
                continue
            c_side, e_side = junction_flanks(flags[b], lo, hi_)
            if c_side is None or not (n_u[b] > 0 and a_g[b] > 0):
                continue
            for x, s_out in ((e_side, n_s[b] + flux[b]), (c_side, n_s[b])):
                if boundary_shares_strand(fp[b], fn[b], fp[x], fn[x]) and a_g[x] > 0:
                    served7.append((b, x, s_out, None))
        pair_width = {}
        for b, x, s_out, _kind in served7:
            if not (tau[b] > 0.0 and tau[x] > 0.0):
                continue
            ok, lo_v = True, []
            for y in (b, x):
                n = cnt[y].sum()
                p = cnt[y, 0] / n
                f = (p - ks_of(y)) / (0.5 - ks_of(y))
                if not 0.0 < f < 1.0:
                    ok = False  # a vertex mode has no log-odds
                    break
                v_log = p * (1.0 - p) / n / (p - ks_of(y)) ** 2
                lo_v.append((np.log(f / (1.0 - f)), v_log / (1.0 - f) ** 2))
            if not ok:
                continue
            (lo_b, v_b), (lo_x, v_x) = lo_v
            v_ratio = s_out / (n_u[b] * (n_u[b] + s_out)) if s_out > 0 else 0.0
            d = lo_b - lo_x - np.log((n_u[b] + s_out) / n_u[b])
            pair_width[(b, x)] = max(0.0, d * d - (v_b + v_x + v_ratio))
        for b, x, s_out, _kind in served7:
            width = pair_width.get((b, x), 0.0)
            if strand_live[x]:
                row = splice_out_row(own_row(x), lam, n_u[b], s_out, a_g[b], a_g[x])
                if np.ptp(row) > EPS:
                    deliver(b, blur_row(row, lam, width))
            if tau[b] > 0.0:
                le = face_map_lambda(lam, n_u[b], a_g[b], a_g[b], a_g[x], a_g[x], s_out / a_g[b])
                row = transport_row(own_row(b), lam, le, n_u[b], s_out)
                if np.ptp(row) > EPS:
                    deliver(x, blur_row(row, lam, width))
        return _PreparedTransfer(rows if live else None)


class _PreparedTransfer:
    """One sweep's delivery: the rows, or silence. ``scan`` relays nothing."""

    def __init__(self, rows: np.ndarray | None):
        self._rows = rows

    def scan(self, *, backward: bool):
        return None

    def deliver(self, left: NeighbourState, right: NeighbourState) -> PsiMessage:
        if self._rows is None:
            return PsiMessage.silent()
        return PsiMessage(lam_rows=self._rows)
