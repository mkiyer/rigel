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

import dataclasses

import numpy as np

from . import NeighbourState, PsiMessage, StepContext
from .foundation import Claim, Hop, Message, PropagationModel, SolveModel
from .variance import count_logvar

#: numerical budget for the level solve's re-linearisation — a convergence tolerance and a
#: per-step log bound (stability, not model constants).
_LEVEL_STEPS = 8
_LEVEL_TOL = 1e-9
_LEVEL_MAX_STEP = 30.0

__all__ = [
    "AllocationSolve",
    "FrameAwarePropagation",
    "RowsSolve",
    "SilentSolve",
    "UnifiedPolicy",
    "allocate",
    "allocate_level",
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
    exon-facing one, route-summed. Precision is counting-honest on the slot's flux count
    (conservative: never sharper than the total flux supports).

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
        False: face(ctx.route_rate_hi, ctx.sj_count, ok_hi),
        True: face(ctx.route_rate_lo, ctx.sj_count, ok_lo),
    }


def count_from_precision(p) -> np.ndarray:
    """Invert the counting law: precision = 1/trigamma(count + ½) ⇒ count. Newton on the
    strictly-decreasing trigamma — float-exact in a handful of steps."""
    from scipy.special import polygamma

    p = np.asarray(p, np.float64)
    target = np.where(p > 0, 1.0 / np.where(p > 0, p, 1.0), np.inf)
    c = np.maximum(1.0 / np.maximum(target, 1e-12) - 0.5, 0.0)
    for _ in range(30):
        f = polygamma(1, c + 0.5) - target
        df = polygamma(2, c + 0.5)
        step = np.where(np.abs(df) > 0, f / np.where(np.abs(df) > 0, df, 1.0), 0.0)
        c = np.maximum(c - step, 0.0)
    return np.where(np.isfinite(target), c, 0.0)


def flank_complete(*, spliced_live, terminus_gain) -> np.ndarray:
    """The ruled completeness branch: a flank is complete iff it carries certified routes AND no
    terminus gain admits molecules the flux cannot see. v1 withholds the flux claim at incomplete
    flanks (the one-sided form is the recorded refinement)."""
    return np.asarray(spliced_live, bool) & ~np.asarray(terminus_gain, bool)


class _PreparedUnified:
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
        # the SHARED LEVEL variance travels beside the lanes: own beliefs are locally anchored
        # (0) and each hop adds its reframe's scale cost, diluted by the destination's own
        n_slots = int(np.asarray(self.ctx.n_slot).shape[0])
        level = np.zeros(n_slots, np.float64)
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
            incoming = dataclasses.replace(self._message_at(state, s), level_logvar=float(level[s]))
            mine = self._message_at(own_scan, i)
            out = self._prop.propagate(mine, incoming, Hop(src=s, dst=i, backward=backward))
            for lane in Message.LANES:
                claim = out.lane(lane)
                state[lane][0][i] = claim.abundance
                state[lane][1][i] = claim.precision
                state[lane][2][i] = claim.measured
            level[i] = out.level_logvar

        def publish():
            return tuple(arr for lane in Message.LANES for arr in state[lane]) + (level,)

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
        lv = np.asarray(nb.state[3 * len(Message.LANES)], np.float64)
        return Message(**lanes, level_logvar=np.where(valid, lv, 0.0))

    def deliver(self, left: NeighbourState, right: NeighbourState) -> PsiMessage:
        own = Message(**{lane: Claim(*self.own[lane]) for lane in Message.LANES})
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

    ⚠ Recorded, deliberately not in this extraction: the sampling
    CAP (needs per-component source counts the lanes do not carry); and the SPLICED transit law
    — certified counts cross UNREFRAMED (capture-invariance, the anchor's measured property) and
    for now cost nothing, because their honest transit price is derived WITH the spliced solve.
    """

    #: A/B flag (the factorial contract: each fix priced in isolation and together):
    #: relay's scan-time conservation projection at licensed hops.
    _mass_rescale = False

    def __init__(self):
        self._tables = None
        self._premise = 0.0

    def propagate(self, own, incoming, hop=None):
        out = super().propagate(own, incoming, hop)
        if not self._mass_rescale or hop is None:
            return out
        # THE MASS RESCALE (relay's scan-time conservation projection, ported): at a licensed
        # hop the running state must satisfy the local identity sum_c rho_c*E_c = M — an
        # OVERWRITE by k = M/S, never a fuse (the slot's observed mass is a fact), so the
        # message carries the COMPOSITION forward while the LEVEL is re-measured at every
        # licensed slot. This is what makes any level-transport rule safe: capture steps,
        # enrichment misjudgments and stale anchors are erased one hop later. Licensed = the
        # incoming state supplied both components (a composition), or the destination is
        # structurally pure-gDNA (its whole observed mass IS gDNA — AXIOM 0). Components the
        # message did not supply enter the budget at the destination's own density.
        t = self._tables[hop.backward]
        i = hop.dst
        g = out.lane("unspliced_gdna")
        p_ = out.lane("unspliced_rna_pos")
        n_ = out.lane("unspliced_rna_neg")
        # THE SUPPLY TEST (currency's form, the sharper of the two donors): the state must
        # carry precision on the gDNA lane AND on EVERY strand the destination admits — else
        # the identity has no right-hand side and the rescale would scale a PARTIAL claim up
        # to account for mass its missing components hold (measured: 73,728 invented gDNA
        # fragments at a zero-gDNA control under the any-RNA-precision form).
        supplied = (
            float(incoming.unspliced_gdna.precision) > 0.0
            and (float(incoming.unspliced_rna_pos.precision) > 0.0 or not bool(t["fp"][i]))
            and (float(incoming.unspliced_rna_neg.precision) > 0.0 or not bool(t["fn"][i]))
        )
        if not (supplied or bool(t["g1"][i])):
            return out
        eg, er, mass = float(t["E_g"][i]), float(t["E_r"][i]), float(t["M"][i])
        vg = float(g.abundance) if float(g.precision) > 0.0 else float(t["own_g"][i])
        vp = float(p_.abundance) if float(p_.precision) > 0.0 else float(t["own_p"][i])
        vn = float(n_.abundance) if float(n_.precision) > 0.0 else float(t["own_n"][i])
        budget = vg * eg + (vp + vn) * er
        if budget <= 1e-9 or mass <= 1e-9:
            return out
        k = mass / budget
        for lane, claim in (
            ("unspliced_gdna", g),
            ("unspliced_rna_pos", p_),
            ("unspliced_rna_neg", n_),
        ):
            if float(claim.precision) > 0.0:
                out = out.with_lane(
                    lane, Claim(float(claim.abundance) * k, claim.precision, claim.measured)
                )
        return out

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
                # the mass-rescale banks: the local conservation identity's ingredients
                "E_g": np.asarray(ctx.eff_gdna_global, np.float64),
                "E_r": np.asarray(ctx.eff_rna, np.float64),
                "M": np.asarray(ctx.mass, np.float64),
                "own_g": np.asarray(ctx.own.rho_g, np.float64),
                "own_p": np.asarray(ctx.own.rho_pos, np.float64),
                "own_n": np.asarray(ctx.own.rho_neg, np.float64),
                "g1": ~fp & ~fn,
            }
        for backward in (False, True):
            t = self._tables[backward]
            src = t["src"]
            # THE SPLICE-SHED RULE (relay's SPLICE OUT, extracted): RNA arriving at a
            # boundary from an EXON sheds the spliced share — the RNA that spliced out at this
            # junction skipped the intron, so it cannot also be inside it. One RNA, two routes:
            # the flux names the part that took the junction; the remainder continues unspliced.
            # The shed is the flux SPANNING the intron the message is entering:
            # forward (heading into the intron on the boundary's genomic-HIGH side) sheds the
            # routes whose LOW end is this boundary; backward sheds the HIGH-ended routes.
            face = np.asarray(ctx.route_rate_hi if backward else ctx.route_rate_lo, np.float64)
            src_is_exon = np.asarray(ctx.is_exon_region, bool)[src]
            gate = np.asarray(ctx.is_boundary, bool) & src_is_exon
            t["shed_pos"] = np.where(gate, face[:, 0], 0.0)
            t["shed_neg"] = np.where(gate, face[:, 1], 0.0)
        lr = np.concatenate([t["log_r"] for t in self._tables.values()])
        vr = np.concatenate([t["v_r"] for t in self._tables.values()])
        # ⛔ ONE pooled scalar, and that is now a MEASURED choice, not a deferral: the class-keyed
        # fit (region-end exon vs not, gated and perturbation-verified) was built and REFUTED on
        # the panel — lowering the sparse family's charge let composition misinformation through
        # (g00-OFF 91.5k → 394.4k, g05-unstranded 567.8k → 848.4k; only g00-ON improved, 27.2k
        # → 21.5k). The premise prices LEVEL transport; the residual damage is COMPOSITION
        # transport, which log-r moments cannot see — pricing that is the open lever, and it is
        # not a premise refinement.
        self._premise = float(premise_logvar(lr, vr))

    def level_cost(self, hop: "Hop | None") -> float:
        """The hop's SHARED SCALE cost: the knob's own estimate variance ``w·v`` plus the
        residual disagreement ``((1−w)·log r)²``. Both describe the reframe factor ``r^w``,
        which multiplies every lane identically — a LEVEL statement, never a composition one."""
        if hop is None or self._tables is None:
            return 0.0
        t = self._tables[hop.backward]
        lr = float(t["log_r"][hop.dst])
        v = float(t["v_r"][hop.dst])
        lr2 = lr * lr
        w = 0.0 if lr2 <= 0.0 else (lr2 / (lr2 + v) if v > 0.0 else 1.0)
        resid = (1.0 - w) * lr
        return w * v + resid * resid

    def attenuate(self, claim: Claim, lane: str, hop: "Hop | None") -> Claim:
        t = self._tables[hop.backward]
        i, s = hop.dst, int(t["src"][hop.dst])
        if lane == "unspliced_rna_pos" and bool(t["fp"][i]) != bool(t["fp"][s]):
            return Claim.silent()
        if lane == "unspliced_rna_neg" and bool(t["fn"][i]) != bool(t["fn"][s]):
            return Claim.silent()
        if lane == "unspliced_rna_pos" and t["shed_pos"][i] > 0.0:
            kept = max(float(claim.abundance) - float(t["shed_pos"][i]), 0.0)
            live = kept > 0.0
            claim = Claim(kept, claim.precision if live else 0.0, claim.measured if live else 0.0)
            if claim.precision <= 0.0:
                return Claim.silent()
        if lane == "unspliced_rna_neg" and t["shed_neg"][i] > 0.0:
            kept = max(float(claim.abundance) - float(t["shed_neg"][i]), 0.0)
            live = kept > 0.0
            claim = Claim(kept, claim.precision if live else 0.0, claim.measured if live else 0.0)
            if claim.precision <= 0.0:
                return Claim.silent()
        lr = float(t["log_r"][i])
        v = float(t["v_r"][i])
        lr2 = lr * lr
        w = 0.0 if lr2 <= 0.0 else (lr2 / (lr2 + v) if v > 0.0 else 1.0)
        # ⭐ THE SPLIT (2026-08-27): the knob's two costs — its own estimate variance w·v and
        # the residual disagreement ((1−w)·log r)² — are properties of the REFRAME, a single
        # scale applied to every lane alike, so they are charged to the message's SHARED LEVEL
        # (`level_cost`) and not to this claim. What remains here is the PREMISE: "a
        # neighbour's abundances apply at all", which is a per-lane composition price.
        # ⛔ THE gDNA-CONTINUITY RULE WAS BUILT, KEYED THREE WAYS, AND REFUTED ON THE PANEL
        # WITHOUT ITS SECOND HALF. The enrichment ratio mixes TWO effects: hybrid capture,
        # which binds nucleic acid regardless of origin and so enriches gDNA and RNA ALIKE at
        # probe positions (owner, 2026-08-26 — gDNA's density genuinely steps there), and
        # EXPRESSION, which moves RNA only. A gDNA level should ride the first and not the
        # second, and the ratio does not separate them; relay's answer is to transport an
        # unsupplied source's level UNSCALED (r_g = 1, may_share_composition) beside the
        # scan-time MASS RESCALE that overwrites the level to each licensed slot's own measured
        # total (region_geometry's capture-landscape ruling: the landscape — capture included —
        # is carried locally by the pure-gDNA population's own measurements).
        # Ported without the rescale: a static per-slot licence broke the knob's telescoping
        # cancellation (a value RATCHET, arriving gdna densities to 3.9e+32); the running-state
        # licence killed the ratchet but still lost capture-ON (g98-unstranded-ON 3.8M -> 5.3M)
        # while a fuse-based pure-gDNA re-anchor lattice measured too weak to repair the drift
        # (a fuse negotiates; the rescale overwrites). The upside is real — continuity halved
        # g05-unstranded-OFF — so porting the rescale, THEN continuity, is the recorded path.
        value = float(claim.abundance) * float(np.exp(w * lr))
        v_hop = self._premise
        p = float(claim.precision)
        ms = float(claim.measured)
        if v_hop > 0.0:
            p = p / (1.0 + p * v_hop)
            ms = ms / (1.0 + ms * v_hop)
        return Claim(abundance=value, precision=p, measured=ms)


def allocate_level(mu, v_own, *, total, total_precision, level_logvar, absorber):
    """⭐⭐⭐ THE LEVEL SOLVE (derived 2026-08-27) — ONE conservation with TWO variance kinds.

    A transported claim's uncertainty is not all of one kind. The reframe knob's costs
    multiply EVERY lane of a message identically, so they are a SHARED SCALE uncertainty
    ``V`` (a level statement); each component additionally carries its own ``v_c`` (a
    composition statement). Spending the shared part as per-component variance is what turns a
    common-mode level error into a composition error — the measured pathology where a
    near-zero gDNA claim eats an unstranded slot's unexplained RNA mass.

    Minimise ``a²/2V + Σ u_c²/2v_c`` subject to ``Σ μ_c e^(a+u_c) = M`` (a the shared log-level
    correction, u_c the per-component ones; the total carries its own counting precision
    ``p_M``, so the constraint is soft). Stationarity around the current point gives::

        D = M − S,  S = Σ y_c            W = V·S² + Σ v_c·y_c²
        a   = V  ·p_M·D·S   / (1 + p_M·W)
        u_c = v_c·p_M·D·y_c / (1 + p_M·W)

    applied as ``y_c ← y_c·exp(a + u_c)`` and re-linearised until the constraint is met — the
    LOG form, so positivity is structural and no clamp-and-redistribute pass exists.

    ⭐ **The three operators built separately are its limits**, which is why they fought when
    stacked: ``V ≫ v_c`` ⇒ ``a → log(M/S)``, ``u → 0`` — the multiplicative mass rescale, the
    composition held; ``V = 0`` ⇒ the additive precision-weighted allocation; one lane
    uninformative ⇒ the residual lands there, which is ``residual_level``'s law.

    A licensed ABSORBER (terminus-admitted or annotation-admitted-but-unwitnessed RNA) takes a
    DEFICIT whole — the known components hold — but cannot absorb an EXCESS.
    """
    mu = np.asarray(mu, np.float64)
    v = np.asarray(v_own, np.float64)
    live = (mu > 0.0) & np.isfinite(v)
    if not live.any():
        return mu.copy()
    y = mu.copy()
    S = float(y[live].sum())
    if absorber and float(total) - S > 0.0:
        return y
    pM = float(total_precision)
    V = float(level_logvar)
    for _ in range(_LEVEL_STEPS):
        S = float(y[live].sum())
        D = float(total) - S
        if abs(D) <= _LEVEL_TOL * max(float(total), 1.0):
            break
        W = V * S * S + float((v[live] * y[live] ** 2).sum())
        denom = 1.0 + pM * W
        if denom <= 0.0 or not np.isfinite(denom):
            break
        a = V * pM * D * S / denom
        u = np.where(live, v * pM * D * y / denom, 0.0)
        step = np.clip(a + u, -_LEVEL_MAX_STEP, _LEVEL_MAX_STEP)
        y = np.where(live, y * np.exp(step), y)
    return y


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
        # the count identity's own precision: M is the UNSPLICED count, so its counting basis
        # is n_slot alone (the spliced fragments are not in it)
        self._pM = 1.0 / count_logvar(np.asarray(ctx.n_slot, np.float64))
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
        # the node's OWN composition variance per arm — what the reception law is keyed on
        # (relay's mismatch_deflate regime): inf = composition-blind (AMBIG / unstranded),
        # messages pass untouched; small = the strand solve has this slot, disagreement kills
        from ..region_init import own_composition_logvar

        self._v_own_g, self._v_own_r = own_composition_logvar(
            ctx.own.f_g, ctx.own.tau_lam, ctx.own.struct_lock
        )
        self._free_pos = np.asarray(ctx.free_pos, bool)
        self._free_neg = np.asarray(ctx.free_neg, bool)
        # ── the spliced law's tables ─────────────────────────────────────────────────────────
        _, self._fg_grid = _logodds_grid(int(ctx.n_grid), float(ctx.logodds_window))
        self._n_slot = np.asarray(ctx.n_slot, np.float64)
        self._eff_rna = np.asarray(ctx.eff_rna, np.float64)
        self._is_exon = np.asarray(ctx.is_exon_region, bool)
        self._is_boundary = np.asarray(ctx.is_boundary, bool)
        self._ss_intron = np.asarray(ctx.ss_intron_boundary, bool)
        self._gain_l = gain[left]
        self._gain_r = gain[right]
        rc_hi = np.asarray(ctx.route_count_hi, np.int64).sum(axis=1)
        rc_lo = np.asarray(ctx.route_count_lo, np.int64).sum(axis=1)
        self._k_left = rc_hi[left]  # the left boundary's acceptor face feeds this exon
        self._k_right = rc_lo[right]
        self._nb_left_is_exon = self._is_exon[left]
        self._nb_right_is_exon = self._is_exon[right]

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
        ms = np.asarray(f.measured, np.float64) + np.asarray(b.measured, np.float64)
        return v, p, ms

    def solve_unspliced(self, own, forward, backward) -> PsiMessage:
        n = np.asarray(own.unspliced_gdna.precision).shape[0]
        val, prec, meas = {}, {}, {}
        for lane in self._ALL:
            # THE COVERAGE RULE: the allocation's components are exactly the claims whose
            # fragments the total counted. The spliced components are the node's OWN lanes —
            # a boundary's T contains its own certified flux (live own lanes cover it); a
            # region's T contains no flux (its own spliced lanes are structurally silent, and
            # the arriving faces stay out — their information enters through solve_spliced's
            # rows instead, which also removes a double-count).
            if lane not in self._UNSPLICED:
                o = own.lane(lane)
                v = np.asarray(o.abundance, np.float64)
                p = np.asarray(o.precision, np.float64)
                ms = np.asarray(o.measured, np.float64)
                val[lane], prec[lane], meas[lane] = v, p, ms
                continue
            v, p, ms = self._fuse_dir(forward.lane(lane), backward.lane(lane))
            if lane in self._UNSPLICED:
                # THE RECEPTION LAW (relay's mismatch_deflate, ported): per stream,
                # p_eff = 1/max(v_stream, G^2 - v_own) with G the log gap between the arriving
                # and the OWN lane value and v_own the node's own COMPOSITION variance —
                # composition-blind (v_own = inf: AMBIG, unstranded) passes untouched, which is
                # exactly what lets messages solve unstranded data while a sighted slot caps a
                # disagreeing claim at ~1/G^2. Exactly one side at zero with own composition
                # evidence is a CONTRADICTION and kills the claim outright, both streams.
                v_own = self._v_own_g if lane == "unspliced_gdna" else self._v_own_r
                vo = np.asarray(own.lane(lane).abundance, np.float64)
                seen = np.isfinite(v_own)
                both = (p > 0) & (v > 0) & (vo > 0) & seen
                g2 = np.zeros(n)
                g2[both] = np.log(v[both] / vo[both]) ** 2
                cap = np.maximum(g2 - np.where(seen, v_own, 0.0), 0.0)
                with np.errstate(divide="ignore"):
                    p = np.where(
                        both & (cap > 0),
                        1.0
                        / np.maximum(np.where(p > 0, 1.0 / np.where(p > 0, p, 1.0), np.inf), cap),
                        p,
                    )
                    ms = np.where(
                        both & (cap > 0),
                        1.0
                        / np.maximum(
                            np.where(ms > 0, 1.0 / np.where(ms > 0, ms, 1.0), np.inf), cap
                        ),
                        ms,
                    )
                ms = np.where(np.isfinite(ms), ms, 0.0)
                contra = seen & (((v > 0) & (p > 0) & (vo <= 0)) | ((v <= 0) & (p > 0) & (vo > 0)))
                p = np.where(contra, 0.0, p)
                ms = np.where(contra, 0.0, ms)
            val[lane], prec[lane], meas[lane] = v, p, ms
        if self._residual_level:
            # RESIDUAL_LEVEL (relay's boundary law, ported): everything the imputed gDNA level
            # cannot explain about a boundary's OWN observed mass is continuing RNA — the
            # generic density deconvolution with the gDNA prior supplied by the message. At a
            # low-gDNA boundary rho_nu -> M/E_r at near-counting precision: the positive half
            # of the pincer (the anchored near-zero gDNA level is the negative half). The
            # split across admissible strands follows the arriving RNA values, half each when
            # both are admissible and unclaimed; delivered as a MEASUREMENT (the donor feeds
            # its mode and measurement streams alike — the level is anchored to the slot's own
            # counted mass).
            from .variance import residual_level

            b = self._is_boundary
            pg_arr = prec["unspliced_gdna"]
            with np.errstate(divide="ignore"):
                v_g = np.where(pg_arr > 0, 1.0 / np.where(pg_arr > 0, pg_arr, 1.0), np.inf)
            rho_nu, v_log, _v_lin = residual_level(
                self._M,
                self._n_slot,
                val["unspliced_gdna"],
                self._E["unspliced_gdna"],
                self._E["unspliced_rna_pos"],
                v_g,
            )
            with np.errstate(divide="ignore"):
                p_nu = np.where(
                    np.isfinite(v_log) & (v_log > 0), 1.0 / np.maximum(v_log, 1e-12), 0.0
                )
            # deliver only the information ABOVE the estimator's own ignorance point: the
            # sigma_f >> 1 branch returns f_R ~ Uniform(0,1) at k = 3 — ZERO knowledge, by the
            # estimator's own declaration — and its precision 1/trigamma(3) must not arrive at
            # psi as counted evidence (measured: half-the-mass-is-RNA claims at every deep
            # gDNA-rich boundary whose imputed gdna claim arrived weak). The subtraction point
            # is the estimator's structure, not a tunable.
            from scipy.special import polygamma

            p_nu = np.maximum(p_nu - 1.0 / float(polygamma(1, 3.0)), 0.0)
            # THE RESIDUAL CLAIM IS AN IMPUTATION and passes through the same reception law
            # as any arriving claim: at a slot with own composition evidence a disagreeing
            # residual is capped at ~1/G^2 (the uninformative f_R = 1/2 branch cannot arrive
            # as a confident measurement against sighted strata); a composition-blind slot
            # passes it untouched — exactly where the claim wins.
            vo_r = np.asarray(own.unspliced_rna_pos.abundance, np.float64) + np.asarray(
                own.unspliced_rna_neg.abundance, np.float64
            )
            seen_r = np.isfinite(self._v_own_r)
            both_r = (p_nu > 0) & (rho_nu > 0) & (vo_r > 0) & seen_r
            g2r = np.zeros(rho_nu.shape[0])
            g2r[both_r] = np.log(rho_nu[both_r] / vo_r[both_r]) ** 2
            cap_r = np.maximum(g2r - np.where(seen_r, self._v_own_r, 0.0), 0.0)
            with np.errstate(divide="ignore"):
                p_nu = np.where(
                    both_r & (cap_r > 0),
                    1.0
                    / np.maximum(
                        np.where(p_nu > 0, 1.0 / np.where(p_nu > 0, p_nu, 1.0), np.inf), cap_r
                    ),
                    p_nu,
                )
            live_nu = b & (p_nu > 0)
            vp0, vn0 = val["unspliced_rna_pos"], val["unspliced_rna_neg"]
            tot0 = vp0 + vn0
            share_p = np.where(
                self._free_pos & ~self._free_neg,
                1.0,
                np.where(
                    ~self._free_pos & self._free_neg,
                    0.0,
                    np.where(tot0 > 0, vp0 / np.where(tot0 > 0, tot0, 1.0), 0.5),
                ),
            )
            for lane, sh, free in (
                ("unspliced_rna_pos", share_p, self._free_pos),
                ("unspliced_rna_neg", 1.0 - share_p, self._free_neg),
            ):
                m_ = live_nu & free & (sh > 0)
                add_v = rho_nu * sh
                add_p = p_nu * sh
                v0, p0, ms0 = val[lane], prec[lane], meas[lane]
                pf = p0 + np.where(m_, add_p, 0.0)
                vf = np.where(
                    m_ & (pf > 0),
                    (p0 * v0 + add_p * add_v) / np.where(pf > 0, pf, 1.0),
                    v0,
                )
                val[lane] = vf
                prec[lane] = np.where(m_, pf, p0)
                meas[lane] = np.where(m_, ms0 + add_p, ms0)
        # THE ALLOCATION — per slot, the five claims against the node's model-free total.
        # THE UNSEEN-COMPONENT ABSORBER: a deficit at a node whose annotation admits an RNA
        # strand NOBODY has evidence about (silent lane, admissible bit set) belongs to that
        # unseen component, never to the weakest live witness — without it, an unstranded
        # library's RNA mass (invisible to every witness: the strand solve is blind and the
        # neighbours' RNA lanes are silent) is force-fed to the near-zero gDNA anchor claim as
        # phantom gDNA (measured: the whole g05-unstranded family).
        # ⭐⭐⭐ THE CONSERVATION IS THE COUNT IDENTITY: sum_c rho_c * E_c = M — each density
        # weighted by ITS OWN opportunity, summing to the slot's OBSERVED unspliced count.
        # That is `region_gdna_geometry`'s stated contract and it is what makes a delivered
        # f_c a COUNT share, which is exactly what psi consumes (`channel` below divides by M).
        # ⛔ The earlier form summed raw DENSITIES against the reciprocal-opportunity total:
        # a different identity, exact at boundaries but WRONG at regions, where
        # `inv_abundance` reads rho * P(w <= ell) (region_geometry's contract; measured on the
        # ladder — median exon P = 0.452, p5 = 0.044, so the total was up to 23x too low at
        # half of all exon slots) — and inconsistent with both the delivery and the mass
        # rescale, which already used the count identity.
        # ⭐ The SPLICED lanes leave the allocation with it: their fragments are not in M
        # (they are counted separately), so they carry no information into the unspliced
        # conservation. Their evidence reaches psi through `solve_spliced`'s rows alone — the
        # owner's rule that spliced fragments impute the adjacent region and do nothing else.
        # The allocation runs in COUNT space (y_c = rho_c * E_c, a plain sum against M) and
        # converts back: `allocate`'s closed form is unchanged, the constraint is now exact.
        E = np.stack([self._E[k] for k in self._UNSPLICED], axis=1)
        mu = np.stack([val[k] for k in self._UNSPLICED], axis=1) * E
        pp = np.stack([prec[k] for k in self._UNSPLICED], axis=1)
        unseen_rna = (self._free_pos & (prec["unspliced_rna_pos"] <= 0)) | (
            self._free_neg & (prec["unspliced_rna_neg"] <= 0)
        )
        absorber = self._absorber | unseen_rna
        live_rows = (pp > 0).any(axis=1) & (self._pM > 0) & (self._M > 0)
        if self._level_solve:
            # ⭐⭐⭐ THE LEVEL SOLVE: one conservation, two variance kinds. V is the SHARED
            # scale uncertainty the arriving messages accumulated (blended by how much each
            # direction contributed, with the node's OWN belief entering at zero — it is
            # locally anchored); v_c is each component's own. The three operators built
            # separately (mass rescale, additive allocation, residual_level) are its limits.
            lf = np.asarray(forward.level_logvar, np.float64)
            lb = np.asarray(backward.level_logvar, np.float64)
            pf = sum(np.asarray(forward.lane(k).precision, np.float64) for k in self._UNSPLICED)
            pb = sum(np.asarray(backward.lane(k).precision, np.float64) for k in self._UNSPLICED)
            po = sum(np.asarray(own.lane(k).precision, np.float64) for k in self._UNSPLICED)
            wsum = pf + pb + po
            V = np.where(wsum > 0, (pf * lf + pb * lb) / np.where(wsum > 0, wsum, 1.0), 0.0)
            with np.errstate(divide="ignore"):
                v_own = np.where(pp > 0, 1.0 / np.where(pp > 0, pp, 1.0), np.inf)
            for i in np.flatnonzero(live_rows):
                mu[i] = allocate_level(
                    mu[i],
                    v_own[i],
                    total=self._M[i],
                    total_precision=self._pM[i],
                    level_logvar=float(V[i]),
                    absorber=bool(absorber[i]),
                )
        elif self._conserve_multiplicative:
            # THE COMMON-MODE FORM (relay's k = M/S, A/B flag): when the live claims are
            # imputations of comparable quality, a shortfall against the slot's mass is a LEVEL
            # error they share — scale them together and leave the COMPOSITION untouched.
            # Shifting the residual additively onto the weakest lane converts a level error
            # into a composition error, which is how a near-zero gDNA claim eats an unstranded
            # slot's unexplained RNA mass. A licensed absorber still takes a deficit whole.
            live = pp > 0
            S = np.where(live, mu, 0.0).sum(axis=1)
            hold = absorber & (S < self._M)
            k = np.where((S > 0) & live_rows & ~hold, self._M / np.where(S > 0, S, 1.0), 1.0)
            mu = np.where(live, mu * k[:, None], mu)
        else:
            for i in np.flatnonzero(live_rows):
                mu[i] = allocate(
                    mu[i],
                    pp[i],
                    total=self._M[i],
                    total_precision=self._pM[i],
                    absorber=bool(absorber[i]),
                )
        with np.errstate(divide="ignore", invalid="ignore"):
            mu = np.where(E > 0, mu / np.where(E > 0, E, 1.0), 0.0)

        # CONVERSION — log-share of the node's own measured mass, clipped into the grid domain
        def channel(lane):
            # MEASUREMENT CITIZENSHIP (relay's law): the delivered precision is the measured
            # stream — the belief stream weighted the value blend and the allocation, but only
            # counted witnesses may arrive at psi as precision
            x = mu[:, self._UNSPLICED.index(lane)]
            p = meas[lane]
            live = (p > 0) & (self._M > 0) & (x > 0)
            share = np.where(live, x * self._E[lane] / np.where(self._M > 0, self._M, 1.0), 1.0)
            mode = np.clip(np.log(np.maximum(share, 1e-300)), self._dom[0], self._dom[1])
            return np.where(live, mode, 0.0), np.where(live, p, 0.0)

        mg, pg = channel("unspliced_gdna")
        mp_, pp_ = channel("unspliced_rna_pos")
        mn_, pn_ = channel("unspliced_rna_neg")
        return PsiMessage(gdna_mode=mg, gdna_prec=pg, rna_mode=(mp_, mn_), rna_prec=(pp_, pn_))

    _spliced_enabled = True
    #: A/B flag (the factorial contract): relay's residual_level at boundary solves.
    _residual_level = False
    #: ⭐ A/B flag: THE LEVEL SOLVE — one conservation with two variance kinds, of which the
    #: additive allocation and the multiplicative rescale are limits (`allocate_level`).
    _level_solve = False
    #: A/B flag: the conservation OPERATOR — multiplicative (relay's common-mode k = M/S,
    #: composition-preserving) instead of the additive precision-weighted allocation.
    _conserve_multiplicative = False

    def solve_spliced(self, own, forward, backward):
        """THE SPLICED LAW, lane-native (commit three): flank rates from the two arriving
        spliced lanes (one-hop reach ⇒ each is exactly the adjacent flank's exon-facing
        route-summed rate); NB size recovered from the counting-honest precisions; the derived
        width law V_route(class) + (V_tail − V_route)+; unspliced-RNA nodes from the ARRIVING RNA
        claims (already deconvolved upstream — no background subtraction, the ruled blend); the
        guarded-Gaussian family at boundaries. Completeness withholds the claim where a terminus
        admits unseen molecules (the ruled v1 branch)."""
        if not self._spliced_enabled:
            return None
        from .. import rna_anchor as RA

        n = self._n_slot.shape[0]
        fg = self._fg_grid

        def lane_pair(msg):
            rp = np.asarray(msg.spliced_rna_pos.abundance)
            pp_ = np.asarray(msg.spliced_rna_pos.precision)
            rn = np.asarray(msg.spliced_rna_neg.abundance)
            pn_ = np.asarray(msg.spliced_rna_neg.precision)
            rate = np.where(pp_ > 0, rp, 0.0) + np.where(pn_ > 0, rn, 0.0)
            cnt = count_from_precision(pp_) + count_from_precision(pn_)
            live = (pp_ > 0) | (pn_ > 0)
            return rate, cnt, live

        rate_l, cnt_l, live_l = lane_pair(forward)
        rate_r, cnt_r, live_r = lane_pair(backward)
        comp_l = flank_complete(spliced_live=live_l, terminus_gain=self._gain_l)
        comp_r = flank_complete(spliced_live=live_r, terminus_gain=self._gain_r)

        def rna_blend(msg):
            a = np.asarray(msg.unspliced_rna_pos.abundance) + np.asarray(
                msg.unspliced_rna_neg.abundance
            )
            p = np.asarray(msg.unspliced_rna_pos.precision) + np.asarray(
                msg.unspliced_rna_neg.precision
            )
            return a, p

        af, pf = rna_blend(forward)
        ab_, pb = rna_blend(backward)
        # a BOUNDARY's unspliced-RNA prediction may read only its NON-exon side: the exon
        # side's RNA largely splices out at the junction, so it is not part of what crosses
        # this boundary unspliced. An exon's input blends both sides (its neighbours' states
        # are intron-side already, post the splice-shed rule)
        fwd_exonic = self._nb_left_is_exon
        bwd_exonic = self._nb_right_is_exon
        pf_b = np.where(self._is_boundary & fwd_exonic, 0.0, pf)
        pb_b = np.where(self._is_boundary & bwd_exonic, 0.0, pb)
        p_nas = pf_b + pb_b
        rho_nas = np.where(
            p_nas > 0, (pf_b * af + pb_b * ab_) / np.where(p_nas > 0, p_nas, 1.0), 0.0
        )

        rows = np.zeros((n, fg.shape[0]))
        C = self._n_slot

        # ── the exon law ─────────────────────────────────────────────────────────────────────
        # the licensed population and the measurability gates — matching the graduated anchor's
        # selection exactly (the sender's publication licence supplies completeness; these add
        # the recipient's own requirements: observed count, RNA opportunity)
        e = self._is_exon & (comp_l | comp_r) & (self._n_slot > 0) & (self._eff_rna > 0)
        if e.any():
            n_fl = comp_l[e].astype(float) + comp_r[e].astype(float)
            rate = (
                np.where(comp_l[e], rate_l[e], 0.0) + np.where(comp_r[e], rate_r[e], 0.0)
            ) / np.maximum(n_fl, 1.0)
            mu = rate * self._eff_rna[e]
            size = np.where(comp_l[e], cnt_l[e], 0.0) + np.where(comp_r[e], cnt_r[e], 0.0) + 0.5
            Ce = C[e]
            fit = RA.left_fit_center_spread(Ce, mu)
            if fit is not None:
                mu = mu * float(np.exp(fit[0]))
            both = comp_l[e] & comp_r[e]
            single = (self._k_left[e] <= 1) & (self._k_right[e] <= 1)
            V_route = {}
            for cls, m in (("single", single), ("multi", ~single)):
                mm = both & m & (rate_l[e] > 0) & (rate_r[e] > 0)
                V_route[cls] = RA.route_pair_log_variance(rate_l[e][mm], rate_r[e][mm])
            V_tail = RA.left_tail_log_variance(Ce, mu)
            v_arr = np.full(int(e.sum()), np.nan)
            for cls, m in (("single", single), ("multi", ~single)):
                vr = V_route[cls]
                if vr is None and V_tail is None:
                    v = None
                elif vr is None:
                    v = V_tail
                elif V_tail is None:
                    v = vr
                else:
                    v = vr + max(0.0, V_tail - vr)
                if v is not None:
                    v_arr[m] = v
            live_v = np.isfinite(v_arr)
            nodes = None
            has_nas = np.zeros(int(e.sum()), bool)  # ⛔ waypoint: the blend-driven
            # unspliced-RNA nodes are disabled pending their strand-aware derivation (measured
            # ~inert either way on the panel); the graduated anchor's intron-read win at the
            # unspliced-RNA stress rows is knowingly given back until the derivation lands
            if has_nas.any():
                sd = np.where(has_nas, 1.0 / np.sqrt(np.maximum(p_nas[e], 1e-300)), 0.0)
                z = np.array([-0.9674, 0.0, 0.9674])  # the 1/6, 1/2, 5/6 normal quantiles
                nodes = (
                    np.maximum(rho_nas[e], 0.0)[:, None]
                    * np.exp(sd[:, None] * z[None, :])
                    * self._eff_rna[e][:, None]
                )
                nodes[~has_nas] = 0.0
            if live_v.any():
                sub = np.flatnonzero(e)[live_v]
                rows[sub] += RA._quadrature_rows(
                    Ce[live_v],
                    mu[live_v],
                    size[live_v],
                    scatter_log_variance=float(np.nanmean(v_arr[live_v])),
                    nascent_count_nodes=None if nodes is None else nodes[live_v],
                    fg_grid=fg,
                )

        # ── the boundary family (guarded Gaussian; cliff-flat behaviour preserved) ──────────
        b = self._ss_intron & (p_nas > 0) & (C > 0) & (self._eff_rna > 0)
        if b.any():
            mu_b = np.maximum(rho_nas[b], 0.0) * self._eff_rna[b]
            fit_b = RA.left_fit_center_spread(C[b], mu_b)
            if fit_b is not None:
                mu_b = mu_b * float(np.exp(fit_b[0]))
            mad_b = RA.left_tail_log_variance(C[b], mu_b)
            # the certified route disagreement joins the width family (the graduated law's
            # guard): at capture-ON it carries the common mode the fit and the MAD are blind
            # to, and without it the boundary rows arrive overconfident at the capture cliff
            both_all = comp_l & comp_r & (rate_l > 0) & (rate_r > 0)
            V_pair_all = RA.route_pair_log_variance(rate_l[both_all], rate_r[both_all])
            ests = [
                v
                for v in (V_pair_all, (fit_b[1] if fit_b is not None else None), mad_b)
                if v is not None
            ]
            if ests:
                V_b = max(ests)
                claim_var = np.where(p_nas[b] > 0, 1.0 / p_nas[b], 0.0)
                var = mu_b + float(np.expm1(V_b)) * mu_b**2 + claim_var * mu_b**2
                rows[np.flatnonzero(b)] += RA._gaussian_rows(C[b], mu_b, var, fg)

        return rows if float(np.ptp(rows)) > 0.0 else None
