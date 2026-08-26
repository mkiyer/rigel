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
from scipy.special import polygamma

from . import NeighbourState, PsiMessage, StepContext
from .foundation import Claim, Message, PropagationModel, SolveModel

__all__ = ["RowsSolve", "SilentSolve", "UnifiedPolicy"]


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
        return _PreparedUnified(ctx, self.propagation, self.solve_model)


def _own_messages(ctx: StepContext) -> dict:
    """Each node's OWN message, from banks the context already carries: the self-solve densities
    and precisions in the unspliced lanes; the certified sj flux, as a counting-honest rate
    claim (precision = 1/trigamma(count + ½), the house counting law), in the spliced lanes."""
    own = ctx.own
    spl = np.asarray(ctx.sj_count, np.float64)
    esp = np.asarray(ctx.eff_sj, np.float64)

    def spliced(col: int):
        live = esp[:, col] > 0.0
        ab = np.where(live, spl[:, col] / np.maximum(esp[:, col], 1e-300), 0.0)
        pr = np.where(live, 1.0 / polygamma(1, spl[:, col] + 0.5), 0.0)
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
        return Message(**{lane: Claim(float(state[lane][0][i]), float(state[lane][1][i]))
                          for lane in Message.LANES})

    def scan(self, *, backward: bool):
        # state starts as each node's OWN message; a hop overwrites the destination with the
        # propagated blend, so a reference terminal (skipped by the backbone) keeps its own.
        state = {lane: (self.own[lane][0].copy(), self.own[lane][1].copy())
                 for lane in Message.LANES}

        def step(s: int, i: int) -> None:
            incoming = self._message_at(state, s)
            mine = self._message_at(self.own, i)
            out = self._prop.propagate(mine, incoming)
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
        own = Message(**{lane: Claim(self.own[lane][0], self.own[lane][1])
                         for lane in Message.LANES})
        return self._solve.solve(own, self._travelled(left), self._travelled(right))
