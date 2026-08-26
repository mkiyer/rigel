"""THE BRIDGE — the foundation spec running on the real backbone (owner directive, 2026-08-26).

⭐⭐⭐ **This module is the CONVERGENCE TARGET, not another sandbox policy.** The Policy
architecture was a development vehicle — a way to build new message code without breaking the
shipped answer — and it is scheduled for retirement: the plan is ONE policy, with the existing
ones treated as DONORS whose measured-good parts are extracted here and whose files are then
deleted (converge-and-delete). `UnifiedPolicy` is the runner that will remain.

What it does: plugs a `(PropagationModel, SolveModel)` pair from the foundation spec into the
backbone's existing protocol. `prepare` builds each node's OWN message from the banks the
context already carries — the self-solve densities into the unspliced lanes, the certified sj
flux into the spliced lanes (which is what retires the separate certified-flux side channel:
the flux is simply part of what a boundary knows). The scans initialise every node's state to
its own message and let `propagate` integrate-and-relay hop by hop (a reference terminal keeps
its own message — the backbone skips terminal hops). `deliver` hands the solve model the two
travelled messages plus the node's own, and the solve is the only place a belief is formed.

The two byte-identity anchors (gated, `tests/calibration/test_unified_bridge.py`):

* pass-through propagation + `SilentSolve` ≡ `SilentPolicy`, byte-for-byte through the real
  backbone — silence IS the solve ignoring the messages, nothing else changes;
* pass-through propagation + `RowsSolve` (the shipped anchor arithmetic as a spliced-lane solve)
  ≡ the relay running only its certified-flux stream, byte-for-byte.

⚠ The scan builds per-hop `Message` objects in a Python loop — correct and honest, not yet fast.
Vectorising the hop kernel is deliberately deferred until a REAL propagation model exists, so
speed work optimises arithmetic that will survive.
"""

from __future__ import annotations

import numpy as np

from . import NeighbourState, PsiMessage, StepContext
from .foundation import Claim, Hop, Message, PropagationModel, SolveModel
from .variance import count_logvar

__all__ = ["FrameAwarePropagation", "RowsSolve", "SilentSolve", "UnifiedPolicy"]


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


class UnifiedPolicy:
    """The spec runner: `(PropagationModel, SolveModel)` on the backbone's Policy protocol."""

    name = "unified"

    def __init__(self, propagation: PropagationModel, solve_model: SolveModel):
        self.propagation = propagation
        self.solve_model = solve_model

    def prepare(self, ctx: StepContext) -> "_PreparedUnified":
        self.propagation.prepare(ctx)
        return _PreparedUnified(ctx, self.propagation, self.solve_model)


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
        return ab, pr

    sp_pos, sp_neg = spliced(0), spliced(1)
    return {
        "unspliced_gdna": (np.asarray(own.rho_g, np.float64), np.asarray(own.prec_g, np.float64)),
        "unspliced_rna_pos": (
            np.asarray(own.rho_pos, np.float64),
            np.asarray(own.prec_pos, np.float64),
        ),
        "unspliced_rna_neg": (
            np.asarray(own.rho_neg, np.float64),
            np.asarray(own.prec_neg, np.float64),
        ),
        "spliced_rna_pos": sp_pos,
        "spliced_rna_neg": sp_neg,
    }


class _PreparedUnified:
    """One sweep's working object: own messages, the scan kernel, and the solve delivery."""

    def __init__(self, ctx: StepContext, propagation: PropagationModel, solve_model: SolveModel):
        self.ctx = ctx
        self._prop = propagation
        self._solve = solve_model
        self.own = _own_messages(ctx)

    def _message_at(self, state: dict, i: int) -> Message:
        return Message(
            **{
                lane: Claim(float(state[lane][0][i]), float(state[lane][1][i]))
                for lane in Message.LANES
            }
        )

    def scan(self, *, backward: bool):
        # state starts as each node's OWN message; a hop overwrites the destination with the
        # propagated blend, so a reference terminal (skipped by the backbone) keeps its own.
        state = {
            lane: (self.own[lane][0].copy(), self.own[lane][1].copy()) for lane in Message.LANES
        }

        def step(s: int, i: int) -> None:
            incoming = self._message_at(state, s)
            mine = self._message_at(self.own, i)
            out = self._prop.propagate(mine, incoming, Hop(src=s, dst=i, backward=backward))
            for lane in Message.LANES:
                claim = out.lane(lane)
                state[lane][0][i] = claim.abundance
                state[lane][1][i] = claim.precision

        def publish():
            return tuple(arr for lane in Message.LANES for arr in state[lane])

        return step, publish

    def _travelled(self, nb: NeighbourState) -> Message:
        valid = np.asarray(nb.valid, bool)
        lanes = {}
        for k, lane in enumerate(Message.LANES):
            ab = np.asarray(nb.state[2 * k], np.float64)
            pr = np.asarray(nb.state[2 * k + 1], np.float64)
            lanes[lane] = Claim(np.where(valid, ab, 0.0), np.where(valid, pr, 0.0))
        return Message(**lanes)

    def deliver(self, left: NeighbourState, right: NeighbourState) -> PsiMessage:
        own = Message(
            **{lane: Claim(self.own[lane][0], self.own[lane][1]) for lane in Message.LANES}
        )
        return self._solve.solve(own, self._travelled(left), self._travelled(right))


class FrameAwarePropagation(PropagationModel):
    """THE FIRST REAL PROPAGATION MODEL — the donors' measured-good parts on the spec.

    Extracted from CurrencyPolicy (imported SIDEWAYS, physically migrated when the donor is
    deleted): the per-hop enrichment reading from belief-free face totals, the KNOB
    ``w = (log r)² / ((log r)² + v)`` interpolating between transporting the abundance (w → 0)
    and believing the frame change in full (w → 1) with no constant, its two derived costs
    (the knob's own estimate variance ``w·v`` and the residual disagreement ``((1−w)·log r)²``),
    and the library premise by precision-weighted moments. From the relay: the PER-STRAND
    population licence — a strand admissible on one side of a hop and not the other is a
    different RNA population, and that lane's claim is refused (value and precision in one
    statement). The gDNA lane always crosses: gDNA is genomically continuous.

    ⚠ Recorded, deliberately not in this first extraction: the per-hop-class premise (the
    population-equality bit is the natural class key — the donor's pooled fit is kept verbatim,
    dead mask and all, until the class-keyed fit is derived on the anchored tree); the sampling
    CAP (needs per-component source counts the lanes do not carry); and the SPLICED transit law
    — certified counts cross UNREFRAMED (capture-invariance, the anchor's measured property) and
    for now cost nothing, because their honest transit price is derived WITH the spliced solve.
    """

    def __init__(self):
        self._tables = None
        self._premise = 0.0

    def prepare(self, ctx) -> None:
        from .currency import enrichment_ratio, premise_logvar
        from .variance import count_logvar

        n = int(np.asarray(ctx.n_slot).shape[0])
        fp = np.asarray(ctx.free_pos, bool)
        fn = np.asarray(ctx.free_neg, bool)
        n_obs = np.asarray(ctx.n_slot, np.float64)
        self._tables = {}
        for backward in (False, True):
            nbr = np.asarray(ctx.right if backward else ctx.left, np.int64)
            src = np.clip(nbr, 0, n - 1)
            r = enrichment_ratio(ctx, backward=backward)
            self._tables[backward] = {
                "log_r": np.log(np.maximum(r, 1e-300)),
                "v_r": count_logvar(n_obs) + count_logvar(n_obs[src]),
                "fp": fp,
                "fn": fn,
                "src": src,
            }
        lr = np.concatenate([t["log_r"] for t in self._tables.values()])
        vr = np.concatenate([t["v_r"] for t in self._tables.values()])
        self._premise = float(premise_logvar(lr, vr))

    def attenuate(self, claim: Claim, lane: str, hop: "Hop | None") -> Claim:
        t = self._tables[hop.backward]
        i, s = hop.dst, int(t["src"][hop.dst])
        if lane.startswith("spliced"):
            return claim  # a measurement: never reframed; its transit price is step 3's derivation
        if lane == "unspliced_rna_pos" and bool(t["fp"][i]) != bool(t["fp"][s]):
            return Claim.silent()
        if lane == "unspliced_rna_neg" and bool(t["fn"][i]) != bool(t["fn"][s]):
            return Claim.silent()
        lr = float(t["log_r"][i])
        v = float(t["v_r"][i])
        lr2 = lr * lr
        w = 0.0 if lr2 <= 0.0 else (lr2 / (lr2 + v) if v > 0.0 else 1.0)
        value = float(claim.abundance) * float(np.exp(w * lr))
        resid = (1.0 - w) * lr
        v_hop = w * v + resid * resid + self._premise
        p = float(claim.precision)
        return Claim(abundance=value, precision=p / (1.0 + p * v_hop) if v_hop > 0.0 else p)
