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

__all__ = [
    "AllocationSolve",
    "FrameAwarePropagation",
    "RowsSolve",
    "SilentSolve",
    "UnifiedPolicy",
    "allocate",
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


class UnifiedPolicy:
    """The spec runner: `(PropagationModel, SolveModel)` on the backbone's Policy protocol."""

    name = "unified"

    def __init__(self, propagation: PropagationModel, solve_model: SolveModel):
        self.propagation = propagation
        self.solve_model = solve_model

    def prepare(self, ctx: StepContext) -> "_PreparedUnified":
        self.propagation.prepare(ctx)
        self.solve_model.prepare(ctx)
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


def allocate(mu, p, *, total, total_precision, absorber):
    """THE OWNER'S CONSERVATION ALLOCATION (ratified 2026-08-26): distribute the residual between
    the arriving component claims and the node's measured total IN PROPORTION TO VARIANCE — the
    certified/strong components hold their values, the weak absorb. The derived repair of the
    equal-weight rescale (a weak claim can no longer be amplified to satisfy a constraint).

    Soft form (the total is itself a measurement): minimise Σ p_i (x_i − μ_i)² +
    p_T (Σx − T)², giving x_i = μ_i + k/p_i with k = p_T·(T − Σμ) / (1 + p_T·Σ(1/p_i)).
    A licensed ABSORBER (terminus-admitted NEW RNA, precision → 0) takes any DEFICIT whole — the
    known components keep their values and the solve is in abundance space, never share space —
    but cannot absorb an EXCESS (new transcription cannot be negative). Clamped at zero with one
    variance-weighted redistribution pass."""
    mu = np.asarray(mu, np.float64)
    p = np.asarray(p, np.float64)
    live = p > 0.0
    if not live.any():
        return mu.copy()
    deficit = float(total) - float(mu[live].sum())
    if absorber and deficit > 0.0:
        return mu.copy()  # the absorber takes the whole deficit; every claim holds
    x = mu.copy()
    inv = np.where(live, 1.0 / np.where(live, p, 1.0), 0.0)
    for _ in range(2):  # one clamp-and-redistribute pass is enough for K = 5
        free = live & (x >= 0.0)
        s_inv = float(inv[free].sum())
        if s_inv <= 0.0:
            break
        d = float(total) - float(x[live].sum())
        k = float(total_precision) * d / (1.0 + float(total_precision) * s_inv)
        x[free] = x[free] + k * inv[free]
        if np.all(x >= 0.0):
            break
        x = np.maximum(x, 0.0)
    return np.maximum(x, 0.0)


class AllocationSolve(SolveModel):
    """THE UNSPLICED SOLVE (commit two of the build): reception → allocation → conversion.

    Reception: the two directions fuse as independent witnesses, then each fused claim YIELDS
    where the node's own evidence disagrees (the DerSimonian–Laird deflation, gap measured
    against the node's OWN lane — never a fused belief; silent own ⇒ no deflation, which is the
    split census's requirement: messages yield exactly where the slot is sighted). Allocation:
    the owner's conservation law across the five arriving claims against the node's model-free
    total, with terminus-licensed absorbers. Conversion: currency's delivery — log-share of the
    node's own measured mass, clipped into the grid domain. v1 emits the three component
    channels; the spliced rows arrive with the next commit (the spliced claims already
    participate in the allocation constraint)."""

    _UNSPLICED = ("unspliced_gdna", "unspliced_rna_pos", "unspliced_rna_neg")
    _ALL = (
        "spliced_rna_pos",
        "spliced_rna_neg",
        "unspliced_gdna",
        "unspliced_rna_pos",
        "unspliced_rna_neg",
    )

    def prepare(self, ctx) -> None:
        from ..region_geometry import terminus_flank_gain
        from ..simplex_logodds import _log_fg, _logodds_grid

        n = int(np.asarray(ctx.n_slot).shape[0])
        self._E = {
            "unspliced_gdna": np.asarray(ctx.eff_gdna_global, np.float64),
            "unspliced_rna_pos": np.asarray(ctx.eff_rna, np.float64),
            "unspliced_rna_neg": np.asarray(ctx.eff_rna, np.float64),
        }
        self._M = np.asarray(ctx.mass, np.float64)
        # the node's MODEL-FREE total density and its counting precision — the conservation total
        inv = np.asarray(ctx.inv_abundance, np.float64)
        self._T = (
            inv
            + np.asarray(ctx.inv_sj_lo, np.float64).sum(axis=1)
            + np.asarray(ctx.inv_sj_hi, np.float64).sum(axis=1)
        )
        n_tot = np.asarray(ctx.n_slot, np.float64) + np.asarray(ctx.spliced_slot, np.float64)
        self._pT = np.where(n_tot > 0, 1.0 / count_logvar(n_tot), 0.0)
        # terminus-licensed absorbers: NEW RNA may start/stop beside a flank whose population
        # grows — coarse v1 licence (any gain at an adjacent boundary), refinement recorded
        rgain, lgain = terminus_flank_gain(ctx.boundary_flags)
        left = np.clip(np.asarray(ctx.left, np.int64), 0, n - 1)
        right = np.clip(np.asarray(ctx.right, np.int64), 0, n - 1)
        gain = np.asarray(rgain, bool) | np.asarray(lgain, bool)
        self._absorber = gain[left] | gain[right] | gain
        lam, _ = _logodds_grid(int(ctx.n_grid), float(ctx.logodds_window))
        dom = _log_fg(lam)
        self._dom = (float(dom[0]), float(dom[-1]))

    @staticmethod
    def _fuse_dir(f, b):
        pf = np.asarray(f.precision, np.float64)
        pb = np.asarray(b.precision, np.float64)
        p = pf + pb
        v = np.where(
            p > 0,
            (pf * np.asarray(f.abundance, np.float64) + pb * np.asarray(b.abundance, np.float64))
            / np.where(p > 0, p, 1.0),
            0.0,
        )
        return v, p

    def solve_unspliced(self, own, forward, backward) -> PsiMessage:
        n = np.asarray(own.unspliced_gdna.precision).shape[0]
        val, prec = {}, {}
        for lane in self._ALL:
            v, p = self._fuse_dir(forward.lane(lane), backward.lane(lane))
            # RECEPTION — the yielding rule: where the node's OWN lane disagrees, the arriving
            # claim loses precision (DerSimonian–Laird; gap vs the own lane, never a fused
            # belief). Own silent ⇒ no deflation: messages yield exactly where the slot is
            # sighted (the split census's requirement).
            po = np.asarray(own.lane(lane).precision, np.float64)
            vo = np.asarray(own.lane(lane).abundance, np.float64)
            both = (p > 0) & (po > 0) & (v > 0) & (vo > 0)
            g2 = np.zeros(n)
            g2[both] = np.log(v[both] / vo[both]) ** 2
            b2 = np.maximum(
                g2
                - np.where(p > 0, 1.0 / np.where(p > 0, p, 1.0), 0.0)
                - np.where(po > 0, 1.0 / np.where(po > 0, po, 1.0), 0.0),
                0.0,
            )
            with np.errstate(divide="ignore"):
                deflated = 1.0 / (np.where(p > 0, 1.0 / np.where(p > 0, p, 1.0), np.inf) + b2)
            p = np.where(both & (b2 > 0), deflated, p)
            val[lane], prec[lane] = v, p
        # THE ALLOCATION — per slot, the five claims against the node's model-free total
        mu = np.stack([val[k] for k in self._ALL], axis=1)
        pp = np.stack([prec[k] for k in self._ALL], axis=1)
        live_rows = (pp > 0).any(axis=1) & (self._pT > 0)
        for i in np.flatnonzero(live_rows):
            mu[i] = allocate(
                mu[i],
                pp[i],
                total=self._T[i],
                total_precision=self._pT[i],
                absorber=bool(self._absorber[i]),
            )

        # CONVERSION — log-share of the node's own measured mass, clipped into the grid domain
        def channel(lane):
            x = mu[:, self._ALL.index(lane)]
            p = pp[:, self._ALL.index(lane)]
            live = (p > 0) & (self._M > 0) & (x > 0)
            share = np.where(live, x * self._E[lane] / np.where(self._M > 0, self._M, 1.0), 1.0)
            mode = np.clip(np.log(np.maximum(share, 1e-300)), self._dom[0], self._dom[1])
            return np.where(live, mode, 0.0), np.where(live, p, 0.0)

        mg, pg = channel("unspliced_gdna")
        mp_, pp_ = channel("unspliced_rna_pos")
        mn_, pn_ = channel("unspliced_rna_neg")
        return PsiMessage(gdna_mode=mg, gdna_prec=pg, rna_mode=(mp_, mn_), rna_prec=(pp_, pn_))

    def solve_spliced(self, own, forward, backward):
        return None
