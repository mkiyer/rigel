"""TransferPolicy — the COMPOSITION TRANSFER policy: one hop, nothing relayed, silence-not-zeros.

       Gate: ``tests/calibration/test_transfer_policy.py``

Every message is a composition — a scale-free statement of gDNA versus RNA — carried across ONE
face by the row constructors of `transfer_rows`, so no level ever crosses a capture cliff and no
constant is anywhere. Four messages ship, each delivered as ``PsiMessage.lam_rows``:

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
    edge_bound_row,
    face_is_licensed,
    face_map_lambda,
    splice_out_row,
    transport_row,
)

__all__ = ["TransferPolicy"]


class TransferPolicy:
    """``rows_at(n_grid, logodds_window)`` supplies the intron factory's per-slot rows on the sweep's
    own grid (``calibrate`` passes its memoized ``_intron_prior_at``; ``None`` means no factory
    evidence and the policy is silent — the rung-0 identity). ``strand = (kappa, od_g, od_r)`` is the
    library's fitted strand model, which the exon → boundary message needs to state an exon's own
    row; ``None`` leaves that message off."""

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
        a_g = np.asarray(ctx.eff_gdna_global, np.float64)
        a_r = np.asarray(ctx.eff_rna, np.float64)
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
            strand_live = is_exon & (fp != fn) & (np.asarray(ctx.own.tau_lam, np.float64) > 0.0)
            cnt = np.asarray(ctx.unspliced_count, np.float64)
            belief = np.asarray(ctx.belief_fg, np.float64)

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
                kappa, od_g, od_r = self._strand
                f_ref = float(belief[e]) if np.isfinite(belief[e]) else 0.5
                row_e = strand_row_logodds(
                    lam, cnt[e, 0], cnt[e, 1], bool(fp[e]), kappa, od_g, od_r, f_ref
                )
                for b, _i, hi in licensed_faces(e):
                    deliver(b, splice_out_row(row_e, lam, n_u[b], float(sc[hi][b]), a_g[b], a_g[e]))
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
