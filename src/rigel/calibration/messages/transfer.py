"""TransferPolicy — the COMPOSITION TRANSFER policy on the two-phase backbone: every node states its
own claim, every directed face carries a recipient's rule, and the two passes do the rest.

       Gate: ``tests/calibration/test_transfer_policy.py``

Every message is a COMPOSITION PROFILE — a max-normalised log-likelihood over the solve grid of the
destination's gDNA share — carried across ONE face by a derived map from `transfer_rows`, so no level
ever crosses a capture cliff and no constant is anywhere. The policy has three parts:

* **Every node's OWN CLAIM** (``prepare``): an intron's factory profile (its density against the
  intergenic background, the rows ``rows_at`` supplies); an exon's or a boundary's own strand profile
  where the solver's derived strand deadband declares the channel live (``strand``); an
  intergenic|exon edge's gDNA COUNT (the level lane: the edge's crossing is structurally pure gDNA).
  ⛔ A claim is data only — never a belief, which already holds the prior and the neighbours.
* **The RECIPIENT's rule per directed face** (``prepare``): absent = STOP (composition cannot cross:
  the recipient holds SILENCE); the identity = FORWARD (the two objects share one unspliced
  population exactly); a map = MODIFY (the face's arithmetic with its counting width and, where two
  witnesses exist, the pair's own discrepancy). The rules ARE the ten shipped messages:

  - intron ⇄ intron|exon boundary: FORWARD both ways (rung 1; item 2 needs a shared single strand).
  - boundary → exon at a LICENSED face (no terminus, the same strand set): the splice-in face map,
    the certified flux capping the claimable gDNA share, widened by the face's counting (rung 2).
  - exon → boundary at that face: the splice-in map read backwards, marginalised over the face's
    spliced-to-unspliced ratio (item 1).
  - intergenic|exon edge → exon: the sign-certified LOWER bound on the exon's gDNA from the edge's
    count — the one level-lane rule; a zero edge is vacuous (rung 3).
  - the exon|exon TERMINUS boundary ⇄ its OUTSIDE exon: the licence counts the spliced crossing,
    ``f_b = f_O (U_b + S_b) / U_b`` (item 5) — the outside exon's message travels the splice-out map,
    the boundary's the face map with the spliced density; the inside flank takes the owner's
    ABUNDANCE-DISCREPANCY map, the step's spread fitted from the served pairs' two witnesses (item 6).
  - the ALTERNATIVE SPLICE SITE ⇄ both flanks: the intron-side flank shares the full unspliced
    crossing, the exon-of-both flank the crossing plus the leaving isoform measured as the route flux;
    each pair widened by its own disagreement beyond counting, nothing pooled (item 7).

* **The two passes and the solve** (``propagate`` / ``solve``): what a node SENDS is its own claim
  composed with what it holds from its far side (a witness product: profiles add), and the recipient's
  rule maps or stops it — so a claim travels as far as the faces admit it, each hop charging its own
  counting width, and no node ever hears its own claim back (the forward pass composes only what came
  from the left, the backward pass only what came from the right). At the solve the two held profiles
  add, and ψ fuses them with the slot's own evidence and the prior.

The laws the policy keeps: the sender publishes its claim unchanged; the recipient decides; a no-claim
stays a no-claim — a flat profile, an absent factory, an evidence-free provider all deliver SILENCE,
never a zero-filled channel; a message is built from the source's claim and the recipient's constants
and observations, never the recipient's belief.

⚠ **What the passes changed against the hand-built one-hop policy they replace (2026-09-04).** A
boundary's own strand profile now also reaches the exon through the face map; an exon forwards what it
holds through EVERY face it has, and a boundary forwards the exon's splice-out profile into its intron;
the hand-built two-hop messages (rung 2, item 5's composed transport) are ordinary hops. Measured before
landing (`policy_prototype.py`, the record in the sandbox): the form wins both halves of the ladder
against the one-hop policy, 7/8 and 7/8, at pass zero and through the pipeline.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from ..simplex_logodds import strand_row_logodds
from . import SILENCE, Message, PsiMessage, StepContext
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


def _norm(row):
    row = np.asarray(row, np.float64)
    return row - row.max()


def _fuse(parts):
    """Independent witnesses about one slot: log-profiles add, then re-normalise."""
    out = None
    for p in parts:
        out = p if out is None else out + p
    return None if out is None else _norm(out)


def _forward(row):
    """FORWARD: the identity — a hop whose two objects share one population exactly."""
    return row


class TransferPolicy:
    """``rows_at(n_grid, logodds_window)`` supplies the intron factory's per-slot profiles on the sweep's
    own grid (``calibrate`` passes its memoized ``_intron_prior_at``; ``None`` means no factory
    evidence and the policy is silent — the rung-0 identity). ``strand = (kappa, od_g, od_r)`` is the
    library's fitted strand model, which the exon's and the boundary's own claims need; ``None``
    leaves every strand claim off."""

    name = "transfer"

    def __init__(self, rows_at: Callable, strand: tuple[float, float, float] | None = None):
        self._rows_at = rows_at
        self._strand = None if strand is None else tuple(float(x) for x in strand)

    # ── phase 0: every node's own claim, and the recipient's rule per directed face ──────────────
    def prepare(self, ctx: StepContext) -> "_PreparedTransfer":
        src = self._rows_at(int(ctx.n_grid), float(ctx.logodds_window))
        if src is None:
            return _PreparedTransfer(None, {}, 0)
        src = np.asarray(src, np.float64)
        n, K = ctx.n_slots, src.shape[1]
        is_bnd = np.asarray(ctx.is_boundary, bool)
        is_exon = np.asarray(ctx.is_exon_region, bool)
        fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
        is_intron = ~is_bnd & ~is_exon & (fp | fn)
        is_intergenic = ~is_bnd & ~is_exon & ~fp & ~fn
        left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
        lam = np.linspace(-float(ctx.logodds_window), float(ctx.logodds_window), K)
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

        # ── every node's OWN CLAIM: data only, never a belief ─────────────────────────────────────
        own: list = [None] * n
        strand_live = np.zeros(n, bool)
        own_row = None
        cnt = None
        if self._strand is not None:
            cnt = np.asarray(ctx.unspliced_count, np.float64)
            belief = np.asarray(ctx.belief_fg, np.float64)
            kappa, od_g, od_r = self._strand
            strand_live = is_exon & (fp != fn) & (tau > 0.0)

            def own_row(x):
                """A slot's OWN strand profile: its unspliced count split by genome strand, the variance
                frozen at its incoming belief (a source-side read); the mode is data only."""
                f_ref = float(belief[x]) if np.isfinite(belief[x]) else 0.5
                return _norm(
                    strand_row_logodds(
                        lam, cnt[x, 0], cnt[x, 1], bool(fp[x]), kappa, od_g, od_r, f_ref
                    )
                )

        for i in np.flatnonzero(is_intron):
            if np.ptp(src[i]) > EPS:
                own[i] = _norm(src[i])
        if own_row is not None:
            for e in np.flatnonzero(strand_live):
                own[e] = own_row(e)
            # a boundary's own profile is a statement about ONE live strand: an AMBIG boundary's split
            # constrains only the tilt, never the gDNA level (the Schur complement the local solve applies)
            for b in np.flatnonzero(is_bnd & (tau > 0.0) & (fp != fn)):
                own[b] = own_row(b)

        # ── the recipient's rule per DIRECTED face ────────────────────────────────────────────────
        rule: dict = {}

        def other_flank(b, e):
            return left[b] if right[b] == e else right[b]

        def intron_exon_pair(b):
            lo, hi_ = left[b], right[b]
            if lo < 0 or hi_ < 0:
                return None, None
            if is_exon[lo] and is_intron[hi_]:
                return lo, hi_
            if is_exon[hi_] and is_intron[lo]:
                return hi_, lo
            return None, None

        # rung 1 / item 2 / rung 2 / item 1: the intron|exon face
        for b in np.flatnonzero(is_bnd):
            e, i = intron_exon_pair(b)
            if e is None:
                continue
            rule[(int(i), int(b))] = _forward  # rung 1: one shared unspliced population
            if boundary_shares_strand(fp[b], fn[b], fp[i], fn[i]):
                rule[(int(b), int(i))] = (
                    _forward  # item 2: the identity under the one-opportunity rule
                )
            hi = 1 if left[e] == b else 0
            if not face_is_licensed(flags[b], fp[e], fn[e], fp[i], fn[i]):
                continue
            if not (n_u[b] > 0 and a_g[b] > 0 and a_r[b] > 0 and a_g[e] > 0 and a_r[e] > 0):
                continue  # a depleted face: silence
            le = face_map_lambda(lam, n_u[b], a_g[b], a_r[b], a_g[e], a_r[e], float(rr[hi][b]))

            def face_rule(row, le=le, b=b, s=float(sc[hi][b])):
                return transport_row(row, lam, le, n_u[b], s)

            def out_rule(row, b=b, e=e, s=float(sc[hi][b])):
                return splice_out_row(row, lam, n_u[b], s, a_g[b], a_g[e])

            rule[(int(b), int(e))] = face_rule  # rung 2
            rule[(int(e), int(b))] = out_rule  # item 1

        # rung 3: the intergenic|exon EDGE — the level lane's one rule. The edge's own claim is its
        # gDNA count (a marker profile: its value is read by the rule, not the row); the recipient
        # exon converts it with its own count into the one-sided profile-likelihood bound.
        for e in np.flatnonzero(is_exon):
            for b in (left[e], right[e]):
                if b < 0 or not is_bnd[b]:
                    continue
                o = other_flank(b, e)
                if o >= 0 and is_intergenic[o] and n_u[b] > 0 and a_g[b] > 0 and a_g[e] > 0:
                    own[b] = np.zeros(K)

                    def edge_rule(row, b=b, e=e):
                        return edge_bound_row(lam, n_u[b], n_u[e], a_g[b], a_g[e])

                    rule[(int(b), int(e))] = edge_rule

        # items 5 and 6: the exon|exon TERMINUS boundary — its outside exon by the spliced-crossing
        # licence, its inside exon by the abundance-discrepancy map
        served6 = []
        for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
            lo, hi_ = left[b], right[b]
            if not (is_exon[lo] and is_exon[hi_]):
                continue
            o, i = outside_flank(flags[b], lo, hi_)
            if o is None:
                continue
            if (
                boundary_shares_strand(fp[b], fn[b], fp[o], fn[o])
                and n_u[b] > 0
                and a_g[b] > 0
                and a_g[o] > 0
            ):

                def out5(row, b=b, o=o):
                    return splice_out_row(row, lam, n_u[b], n_s[b], a_g[b], a_g[o])

                le5 = face_map_lambda(lam, n_u[b], a_g[b], a_g[b], a_g[o], a_g[o], n_s[b] / a_g[b])

                def back5(row, le=le5, b=b):
                    return transport_row(row, lam, le, n_u[b], n_s[b])

                rule[(int(o), int(b))] = out5
                rule[(int(b), int(o))] = back5
            if (
                self._strand is not None
                and boundary_shares_strand(fp[b], fn[b], fp[i], fn[i])
                and n_u[b] > 0
                and n_u[i] > 0
                and a_g[b] > 0
                and a_g[i] > 0
            ):
                served6.append((b, i, (n_u[i] / a_g[i]) / ((n_u[b] + n_s[b]) / a_g[b])))
        if served6:
            # the step's spread: the hop's premise, FITTED from the served pairs' two witnesses (the
            # boundary's and the inside exon's strand modes, data only) by the method of moments
            # beyond counting — zero where they agree within counting
            kappa = self._strand[0]
            ks_of = lambda x: kappa if fp[x] else 1.0 - kappa  # noqa: E731
            log_s, var_s = [], []
            for b, i, r in served6:
                if not (tau[b] > 0.0 and tau[i] > 0.0):
                    continue  # one witness only: nothing to fit from this pair
                ok, modes = True, []
                for x in (b, i):
                    nn = cnt[x].sum()
                    p = cnt[x, 0] / nn
                    f = (p - ks_of(x)) / (0.5 - ks_of(x))
                    if not 0.0 < f < 1.0:
                        ok = False  # a vertex mode has no logarithm
                        break
                    modes.append((f, p * (1.0 - p) / nn / (p - ks_of(x)) ** 2))
                if not ok:
                    continue
                (f_b, v_b), (f_i, v_i) = modes
                log_s.append(np.log(f_i / (f_b * n_u[b] / (n_u[b] + n_s[b]) / r)))
                var_s.append(v_b + v_i + 1.0 / n_u[i] + 1.0 / (n_u[b] + n_s[b]))
            v_step = 0.0
            if len(log_s) >= 2:
                ls, vs = np.asarray(log_s), np.asarray(var_s)
                w = 1.0 / vs
                w /= w.sum()
                mu = float(w @ ls)
                v_step = max(0.0, float(w @ (ls - mu) ** 2) - float(w @ vs))
            for b, i, r in served6:

                def step6(row, b=b, i=i, r=r, v=v_step):
                    return abundance_row(row, lam, n_u[b], n_s[b], r, n_u[i], v)

                rule[(int(b), int(i))] = step6

        # item 7: the ALTERNATIVE SPLICE SITE, both flanks both ways, the discrepancy rule PER PAIR
        if self._strand is not None:
            kappa = self._strand[0]
            ks_of = lambda x: kappa if fp[x] else 1.0 - kappa  # noqa: E731
            for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
                lo, hi_ = left[b], right[b]
                if not (is_exon[lo] and is_exon[hi_]):
                    continue
                c_side, e_side = junction_flanks(flags[b], lo, hi_)
                if c_side is None or not (n_u[b] > 0 and a_g[b] > 0):
                    continue
                for x, s_out in ((e_side, n_s[b] + flux[b]), (c_side, n_s[b])):
                    if not (boundary_shares_strand(fp[b], fn[b], fp[x], fn[x]) and a_g[x] > 0):
                        continue
                    width = 0.0
                    if tau[b] > 0.0 and tau[x] > 0.0:
                        ok, lo_v = True, []
                        for y in (b, x):
                            nn = cnt[y].sum()
                            p = cnt[y, 0] / nn
                            f = (p - ks_of(y)) / (0.5 - ks_of(y))
                            if not 0.0 < f < 1.0:
                                ok = False  # a vertex mode has no log-odds
                                break
                            v_log = p * (1.0 - p) / nn / (p - ks_of(y)) ** 2
                            lo_v.append((np.log(f / (1.0 - f)), v_log / (1.0 - f) ** 2))
                        if ok:
                            (lo_b, v_b), (lo_x, v_x) = lo_v
                            v_ratio = s_out / (n_u[b] * (n_u[b] + s_out)) if s_out > 0 else 0.0
                            d = lo_b - lo_x - np.log((n_u[b] + s_out) / n_u[b])
                            width = max(0.0, d * d - (v_b + v_x + v_ratio))
                    le7 = face_map_lambda(
                        lam, n_u[b], a_g[b], a_g[b], a_g[x], a_g[x], s_out / a_g[b]
                    )

                    def out7(row, b=b, x=x, s=s_out, w=width):
                        return blur_row(splice_out_row(row, lam, n_u[b], s, a_g[b], a_g[x]), lam, w)

                    def back7(row, le=le7, b=b, s=s_out, w=width):
                        return blur_row(transport_row(row, lam, le, n_u[b], s), lam, w)

                    rule[(int(x), int(b))] = out7
                    rule[(int(b), int(x))] = back7

        return _PreparedTransfer(own, rule, K)


class _PreparedTransfer:
    """One sweep's working object: every node's own claim, the rules, the two passes' state."""

    def __init__(self, own, rule: dict, n_grid: int):
        self.own = own
        self.rule = rule
        self._K = int(n_grid)
        self.held: dict = {False: None, True: None}

    # ── phase 1: propagate — the recipient's kernel, run by the backbone in chain order ──────────
    def propagate(self, *, backward: bool):
        if self.own is None or not self.rule:
            return None  # nothing to say anywhere: every node holds SILENCE
        held: list = [None] * len(self.own)
        self.held[bool(backward)] = held

        def receive(s: int, i: int):
            """What ``s`` SENDS is its own claim composed with what it holds from its far side (written
            by this pass one step earlier); the rule for the face maps it (MODIFY), passes it (FORWARD)
            or is absent (STOP)."""
            fn = self.rule.get((int(s), int(i)))
            if fn is None:
                return SILENCE
            parts = []
            if self.own[s] is not None:
                parts.append(self.own[s])
            far = held[s]
            if far is not None and far.composition is not None:
                parts.append(far.composition)
            sending = _fuse(parts)
            if sending is None:
                return SILENCE
            out = fn(sending)
            if out is None or np.ptp(out) <= EPS:
                return SILENCE
            msg = Message(composition=_norm(out))
            held[i] = msg
            return msg

        return receive

    # ── phase 2: solve — the policy's half ───────────────────────────────────────────────────────
    def solve(self, from_left: list, from_right: list) -> PsiMessage:
        if self.own is None:
            return PsiMessage.silent()
        n = len(self.own)
        rows = np.zeros((n, self._K))
        live = False
        for i in range(n):
            parts = [
                m.composition
                for m in (from_left[i], from_right[i])
                if m is not None and m.composition is not None
            ]
            if parts:
                rows[i] = _fuse(parts)
                live = True
        return PsiMessage(lam_rows=rows) if live else PsiMessage.silent()
