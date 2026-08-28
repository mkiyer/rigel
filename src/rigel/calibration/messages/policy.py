"""`MessagePolicy` — the foundation spec running on the calibration backbone.

WHAT IT IS. `MessagePolicy` plugs a `(PropagationModel, SolveModel)` pair from
`messages/foundation.py` — the ratified message architecture: one `Message` with provenance
lanes, the propagate/solve timepoints, and the laws the skeleton enforces — into the
backbone's existing `Policy` protocol. It supplies the parts that are the same whatever the
models do: each node's OWN message, built from banks the context already carries (the
self-solve densities in the unspliced lanes, the certified sj flux in the spliced lanes), the
two scans, and the delivery of two travelled messages plus the node's own to the solve.

WHAT IT IS FOR. Message propagation exists for the slots whose own solve has no composition
channel — unstranded data and AMBIG slots, where the strand likelihood is flat and the local
answer is a default rather than a measurement. ⭐ The bar is NOT to beat `SilentPolicy`: on
strand-specific data a sighted exon's own solve is excellent and messages can only disturb it.
The goal is to solve the unstranded case while doing as little harm as possible relative to
silence.

WHERE IT STARTS. With `PassThroughPropagation` and `SilentSolve` this policy is byte-identical
to `SilentPolicy` through the real backbone (gated in
`tests/calibration/test_message_policy.py`), because silence IS the solve ignoring its
messages — a node's belief only updates at the solve.

`SilentPolicy` is the measured floor; `RelayPolicy` is the shipped tool. Both are frozen and
kept for comparison.

⚠ The scan builds per-hop `Message` objects in a Python loop: correct and readable, not fast.
"""

from __future__ import annotations

import numpy as np

from . import NeighbourState, PsiMessage, StepContext
from .foundation import Claim, Hop, Message, PropagationModel, SolveModel
from .variance import count_logvar

__all__ = [
    "MessagePolicy",
    "RowsSolve",
    "SilentSolve",
]


class SilentSolve(SolveModel):
    """The owner's trivial silence: the solve IGNORES the messages — every delivered precision
    is zero — and nothing else needs to exist, because a node's belief only updates at the
    solve."""

    def solve_unspliced(self, own, forward, backward) -> PsiMessage:
        return PsiMessage.silent()

    def solve_spliced(self, own, forward, backward):
        return None


class RowsSolve(SilentSolve):
    """A spliced-lane solve that delivers a precomputed row-factor — the shipped certified-flux
    arithmetic riding the spec's seam. The identity anchor against the channels-off relay."""

    def __init__(self, rows):
        self._rows = np.asarray(rows, np.float64)

    def solve_spliced(self, own, forward, backward):
        return self._rows


class MessagePolicy:
    """The spec runner: `(PropagationModel, SolveModel)` on the backbone's Policy protocol."""

    name = "message"

    def __init__(self, propagation: PropagationModel, solve_model: SolveModel):
        self.propagation = propagation
        self.solve_model = solve_model

    def prepare(self, ctx: StepContext) -> "_PreparedMessage":
        self.propagation.prepare(ctx)
        self.solve_model.prepare(ctx)
        return _PreparedMessage(ctx, self.propagation, self.solve_model)


def _own_messages(ctx: StepContext) -> dict:
    """Each node's OWN message, from banks the context already carries: the self-solve densities
    and precisions in the unspliced lanes; the certified sj flux, as a counting-honest rate
    claim (precision = 1/trigamma(count + ½), the house counting law), in the spliced lanes."""
    own = ctx.own
    spl = np.asarray(ctx.sj_count, np.float64)
    esp = np.asarray(ctx.eff_sj, np.float64)
    # the ROUTE-SUMMED rate (both faces — a boundary's certified flux, correctly pooled at the
    # geometry source); the pooled ratio under-reads k-route faces ~k×
    rr = np.asarray(ctx.route_rate_lo, np.float64) + np.asarray(ctx.route_rate_hi, np.float64)

    def spliced(col: int):
        live = esp[:, col] > 0.0
        ab = np.where(live, rr[:, col], 0.0)
        pr = np.where(live, 1.0 / count_logvar(spl[:, col]), 0.0)
        return ab, pr, pr  # a certified count IS a measurement: both streams

    # MEASUREMENT SEEDS (the shipped relay's citizenship, ported): the gDNA lane's measured
    # stream is the struct-locked anchor's own counting precision — the one unspliced witness
    # that COUNTED something (an intergenic pure-gDNA region's fragments are gDNA by
    # structure). Deconvolved RNA beliefs seed no measurement: a solve's output is not a count.
    lock = np.asarray(own.struct_lock, bool)
    pg = np.asarray(own.prec_g, np.float64)
    zero = np.zeros_like(pg)
    sp_pos, sp_neg = spliced(0), spliced(1)
    return {
        "unspliced_gdna": (np.asarray(own.rho_g, np.float64), pg, np.where(lock, pg, 0.0)),
        "unspliced_rna_pos": (
            np.asarray(own.rho_pos, np.float64),
            np.asarray(own.prec_pos, np.float64),
            zero,
        ),
        "unspliced_rna_neg": (
            np.asarray(own.rho_neg, np.float64),
            np.asarray(own.prec_neg, np.float64),
            zero,
        ),
        "spliced_rna_pos": sp_pos,
        "spliced_rna_neg": sp_neg,
    }


def _own_spliced_faces(ctx: StepContext) -> dict:
    """Per-direction spliced lanes: forward (low→high) presents each boundary's HIGH (acceptor)
    face; backward presents the LOW (donor) face — so the delivered flank rate is exactly the
    exon-facing one, route-summed. Precision is counting-honest on THAT FACE's own flux count
    (`sj_count_lo`/`sj_count_hi`), never the slot's lo+hi total.

    ⚠ The two differ only where a position is a donor for one intron AND an acceptor for
    another, which is RARE — measured on the ladder's annotation: 1 of ~3,600 publishing faces
    per strand (0.03 %). But where it happens the total is not a small over-statement: the
    worst face there reads 176x its own count's precision, because trigamma is decreasing and
    the two faces' fluxes can differ by orders of magnitude. Pricing a face on its own count
    costs nothing and removes the tail.

    ⭐ THE PUBLICATION LICENCE (sender-side): a face is published only toward an exon whose
    flank is STRUCTURALLY COMPLETE — every route into it certified, no terminus admitting
    molecules the flux cannot see. An incomplete face's rate is a KNOWN underestimate; the
    honest statement is silence, not a truth claim at a value known to be low. Measured: the
    lane-liveness gate alone admitted 12,448 exons the structural claim refuses (3.6× the
    licensed population, 1.2M fragments of C at median C/mu 0.003) and those slots were the
    whole of the ③ zero-control regression."""
    n = int(np.asarray(ctx.n_slot).shape[0])
    left = np.clip(np.asarray(ctx.left, np.int64), 0, n - 1)
    right = np.clip(np.asarray(ctx.right, np.int64), 0, n - 1)
    has_l = np.asarray(ctx.left, np.int64) >= 0
    has_r = np.asarray(ctx.right, np.int64) >= 0
    comp_l = np.asarray(ctx.left_interface_certified, bool)
    comp_r = np.asarray(ctx.right_interface_certified, bool)
    # hi face feeds the exon on the boundary's RIGHT (its left flank); lo the mirror
    ok_hi = has_r & comp_l[right]
    ok_lo = has_l & comp_r[left]

    def face(rate, count, ok):
        out = {}
        for col, lane in ((0, "spliced_rna_pos"), (1, "spliced_rna_neg")):
            r = np.asarray(rate, np.float64)[:, col]
            c = np.asarray(count, np.float64)[:, col]
            live = ok & (r > 0.0)
            pr = np.where(live, 1.0 / count_logvar(c), 0.0)
            out[lane] = (np.where(live, r, 0.0), pr, pr)
        return out

    return {
        False: face(ctx.route_rate_hi, ctx.sj_count_hi, ok_hi),
        True: face(ctx.route_rate_lo, ctx.sj_count_lo, ok_lo),
    }


class _PreparedMessage:
    """One sweep's working object: own messages, the scan kernel, and the solve delivery."""

    def __init__(self, ctx: StepContext, propagation: PropagationModel, solve_model: SolveModel):
        self.ctx = ctx
        self._prop = propagation
        self._solve = solve_model
        self.own = _own_messages(ctx)
        self._own_faces = _own_spliced_faces(ctx)

    def _message_at(self, state: dict, i: int) -> Message:
        return Message(
            **{
                lane: Claim(
                    float(state[lane][0][i]),
                    float(state[lane][1][i]),
                    float(state[lane][2][i]),
                )
                for lane in Message.LANES
            }
        )

    def scan(self, *, backward: bool):
        # state starts as each node's OWN message; a hop overwrites the destination with the
        # propagated blend, so a reference terminal (skipped by the backbone) keeps its own.
        state = {lane: tuple(arr.copy() for arr in self.own[lane]) for lane in Message.LANES}
        # ⭐ the SPLICED lanes are direction-aware: travelling low→high a boundary presents its
        # HIGH (acceptor) face — the routes entering its rightward exon — and the mirror going
        # back, so the adjacent solve receives exactly the exon-facing flank rate.
        own_scan = dict(self.own)
        if self._own_faces is not None:
            face = self._own_faces[bool(backward)]
            for lane in ("spliced_rna_pos", "spliced_rna_neg"):
                state[lane] = tuple(arr.copy() for arr in face[lane])
                # ⛔ the scan's own-message must present the SAME masked face, or propagate's
                # fuse re-publishes the unmasked slot-total lane over the licence at every
                # visited slot (measured: 3,663 unlicensed exons received faces this way)
                own_scan[lane] = face[lane]

        def step(s: int, i: int) -> None:
            incoming = self._message_at(state, s)
            mine = self._message_at(own_scan, i)
            out = self._prop.propagate(mine, incoming, Hop(src=s, dst=i, backward=backward))
            for lane in Message.LANES:
                claim = out.lane(lane)
                state[lane][0][i] = claim.abundance
                state[lane][1][i] = claim.precision
                state[lane][2][i] = claim.measured

        def publish():
            return tuple(arr for lane in Message.LANES for arr in state[lane])

        return step, publish

    def _travelled(self, nb: NeighbourState) -> Message:
        valid = np.asarray(nb.valid, bool)
        lanes = {}
        for k, lane in enumerate(Message.LANES):
            ab = np.asarray(nb.state[3 * k], np.float64)
            pr = np.asarray(nb.state[3 * k + 1], np.float64)
            ms = np.asarray(nb.state[3 * k + 2], np.float64)
            lanes[lane] = Claim(
                np.where(valid, ab, 0.0), np.where(valid, pr, 0.0), np.where(valid, ms, 0.0)
            )
        return Message(**lanes)

    def deliver(self, left: NeighbourState, right: NeighbourState) -> PsiMessage:
        own = Message(**{lane: Claim(*self.own[lane]) for lane in Message.LANES})
        return self._solve.solve(own, self._travelled(left), self._travelled(right))
