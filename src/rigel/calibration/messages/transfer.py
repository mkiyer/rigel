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
  - intergenic|exon edge → exon: THE EDGE'S LEVEL, one-sided (rule 5 as a level, the owner's design
    2026-09-04): the exon has at least the edge's gDNA density, at the count's Poisson width; nothing
    above (no local witness prices capture's enrichment of the interior), a zero count vacuous
    (darkness is not absence).
  - the exon|exon TERMINUS boundary ⇄ its OUTSIDE exon: the licence counts the spliced crossing,
    ``f_b = f_O (U_b + S_b) / U_b`` (item 5) — the outside exon's message travels the splice-out map,
    the boundary's the face map with the spliced density.
  - a TERMINUS boundary → the region INSIDE it (exon|exon and exon|intron alike): THE LEVEL RULE
    (the owner's design, 2026-09-04). Composition cannot cross a terminus (new transcription starts
    or ends there), the gDNA LEVEL can: the boundary's OWN strand profile — a measurement, never what
    it holds — carried through the level-kept map (the inside's gDNA share is the boundary's share
    times its crossing density, over the inside's own total), its shape preserved, blurred by the two
    totals' counting and by the pair's own discrepancies: the excess of the totals' disagreement over
    counting, and, where the inside exon's strand channel is live, the excess of the two strand
    modes' disagreement over counting. Value kept, nothing pooled, no hypothesis chosen for the
    discrepancy. With no own claim the boundary still sends what every library measures: its
    crossing's total bounds the inside's gDNA density above.
  - the ALTERNATIVE SPLICE SITE ⇄ both flanks: the intron-side flank shares the full unspliced
    crossing, the exon-of-both flank the crossing plus the leaving isoform measured as the route flux;
    each pair widened by its own disagreement beyond counting, nothing pooled (item 7).
  - EVERY OTHER DIRECTED FACE — a strand change, termini pointing both ways, the AMBIG complex, and
    every face into or out of an EMPTY node (a piece with no total) — carries THE LEVEL LANE
    (2026-09-05): the gDNA level as an ABSOLUTE profile over the log density, which needs no map and
    no recipient. An empty node forwards it unchanged; a full node emits the intersection of its own
    level's lower side and the level it holds; the recipient prices the hop (both totals' counting
    plus the abundance discrepancy beyond it) and takes it as a LOWER bound — a level that crosses a
    face says "at least this much gDNA" and nothing more. The only faces with no rule at all lead into
    intergenic regions (structurally pure gDNA) or off the chain.
  - THE RNA LEVEL LANES (2026-09-08, the both-stranded locus): one lane per strand, its FACES from
    the flag bits (strand ``s`` crosses a face iff the boundary carries none of ``s``'s bits and both
    nodes admit ``s``; across ``s``'s own junction it enters ``s``'s INTRON — the intron test is per
    strand, from ``StepContext.exon_pos`` / ``exon_neg`` — and not ``s``'s exon; a terminus of ``s``
    stops ``s``), TWO-SIDED only between an intron of ``s`` and its own boundary (one shared unspliced
    population: the whole profile, counting-only price), lower-only and priced by the strand's own
    counts everywhere else. SOURCES: a single-strand node's own claim read as its live strand's RNA
    level; the certified flux at each of an exon's junctions as that strand's level at the exon (one
    hop, boundary → exon; two junctions pay their pair's disagreement beyond counting). DELIVERED at
    AMBIG nodes only, as one row over ψ's (λ, θ) cube (`PsiMessage.cube_rows`): a lower bound on RNA+
    is an upper bound on the gDNA share through the node's own strand counts — the side the gDNA lane
    cannot give (THE BRACKET THEOREM, gated). The tilt needs no lane of its own.

* **The two passes and the solve** (``propagate`` / ``solve``): a node SENDS two things apart — its
  own claim (a measurement) and what it holds from its far side (an imputation) — and the recipient's
  rule decides what to do with each: a composition rule composes them (a witness product: profiles
  add) and maps the product; a level rule reads the measurement only, because an imputation is never
  re-issued as a level; the lane intersects them as bounds; either may stop — so a claim travels as
  far as the faces admit it, each hop charging its own counting width, and no node ever hears its own
  claim back (the forward pass composes only what came from the left, the backward pass only what came
  from the right). At the solve the two held profiles add, the two held levels intersect and join them
  as a constraint, and ψ fuses the row with the slot's own evidence and the prior.

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

from ..simplex_logodds import _tilt_grid, strand_row_logodds
from . import SILENCE, Level, Message, PsiMessage, StepContext
from .transfer_rows import (
    EPS,
    SJ_FLAGS,
    blur_row,
    boundary_shares_strand,
    count_logvar,
    count_price,
    cube_row,
    edge_level_row,
    face_is_licensed,
    face_map_lambda,
    flux_level,
    hop_price,
    intersect,
    junction_exon_side,
    junction_flanks,
    level_bound_row,
    level_map_lambda,
    level_of_profile,
    level_row,
    lower_side,
    outside_flank,
    poisson_level,
    priced_level,
    profile_of_level,
    read_column,
    rna_level_of_profile,
    rna_row_of_level,
    splice_out_row,
    strand_bits,
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


def _forward(own, held):
    """FORWARD: the identity — a hop whose two objects share one population exactly."""
    return _fuse([r for r in (own, held) if r is not None])


def _composed(mapping):
    """A composition rule: the sender's own claim and what it holds compose (profiles add) and the
    face's map carries the product."""

    def rule(own, held):
        sending = _fuse([r for r in (own, held) if r is not None])
        return None if sending is None else mapping(sending)

    return rule


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

            rule[(int(b), int(e))] = _composed(face_rule)  # rung 2
            rule[(int(e), int(b))] = _composed(out_rule)  # item 1

        # rule 5: the intergenic|exon EDGE — a LEVEL, one-sided (the owner's design, 2026-09-04; the
        # form the ladder kept). The edge's own claim is its gDNA COUNT (a marker profile; the rule reads
        # the count), and the exon converts it through its own total: the exon has at least the edge's
        # gDNA density, at the count's own Poisson width; nothing above (no local witness prices a
        # probed interior's enrichment over its edge), and a zero count is vacuous (darkness under
        # capture is not absence).
        for e in np.flatnonzero(is_exon):
            for b in (left[e], right[e]):
                if b < 0 or not is_bnd[b]:
                    continue
                o = other_flank(b, e)
                if not (o >= 0 and is_intergenic[o] and a_g[b] > 0 and a_g[e] > 0 and n_u[e] > 0):
                    continue
                own[b] = np.zeros(K)  # the level claim (the rule reads the count, not the row)
                row_edge = edge_level_row(lam, n_u[b], n_u[e], a_g[b], a_g[e])

                def edge_rule(own_claim, held, row=row_edge):
                    """The edge's level: its count, converted at the exon; nothing is held at an edge."""
                    return row if np.ptp(row) > EPS else None

                rule[(int(b), int(e))] = edge_rule

        # item 5: the exon|exon TERMINUS boundary and its OUTSIDE exon — composition, by the
        # spliced-crossing licence; and THE LEVEL RULE into the region INSIDE every terminus
        for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
            lo, hi_ = left[b], right[b]
            o, i = outside_flank(flags[b], lo, hi_)
            if o is None:
                continue
            # the sj+terminus boundary (2026-09-08): a junction sharing the terminus does not change
            # which flank is inside; its measured flux belongs to the flank on the junction's exon side
            ex_side = junction_exon_side(flags[b], lo, hi_) if int(flags[b]) & SJ_FLAGS else None
            flux_b = float(flux[b]) if ex_side is not None else 0.0
            if (
                is_exon[lo]
                and is_exon[hi_]
                and boundary_shares_strand(fp[b], fn[b], fp[o], fn[o])
                and n_u[b] > 0
                and a_g[b] > 0
                and a_g[o] > 0
            ):
                # the outside flank holds the crossing plus, when the junction's exon lies outside, the
                # isoform that splices out there — item 7's arithmetic
                s_out = n_s[b] + (flux_b if ex_side == o else 0.0)

                def out5(row, b=b, o=o, s=s_out):
                    return splice_out_row(row, lam, n_u[b], s, a_g[b], a_g[o])

                le5 = face_map_lambda(lam, n_u[b], a_g[b], a_g[b], a_g[o], a_g[o], s_out / a_g[b])

                def back5(row, le=le5, b=b, s=s_out):
                    return transport_row(row, lam, le, n_u[b], s)

                rule[(int(o), int(b))] = _composed(out5)
                rule[(int(b), int(o))] = _composed(back5)
            # ── THE LEVEL RULE: the boundary's gDNA level into the inside region ─────────────────
            if not (is_exon[i] and (is_exon[o] or is_intron[o])):
                continue
            if not (n_u[b] > 0 and n_u[i] > 0 and a_g[b] > 0 and a_g[i] > 0):
                continue
            if not boundary_shares_strand(fp[b], fn[b], fp[i], fn[i]):
                continue  # a strand change is the plan's step B
            density_b = n_u[b] / a_g[b]
            total_b = n_u[b] + n_s[b]
            # the pair's discrepancies, each the excess of a disagreement over its counting:
            # (i) the totals per opportunity, the dampening the owner ruled — the boundary's total gains
            # the junction's flux when the junction's exon is the inside (the RNA joining there is
            # measured; only the terminus's own transcription is not)
            total_i = total_b + (flux_b if ex_side == i else 0.0)
            r = (n_u[i] / a_g[i]) / (total_i / a_g[b])
            v_pair = max(0.0, float(np.log(r)) ** 2 - (1.0 / n_u[i] + 1.0 / total_i))
            # (ii) predicts the inside's share from the crossing alone: the flux is not a crossing
            r_mode = (n_u[i] / a_g[i]) / (total_b / a_g[b])
            # (ii) the two strand modes, where both channels are live (data only, as item 7)
            if self._strand is not None and tau[b] > 0.0 and tau[i] > 0.0:
                kappa = self._strand[0]
                ok, modes = True, []
                for y in (b, i):
                    nn = cnt[y].sum()
                    p = cnt[y, 0] / nn
                    ks = kappa if fp[y] else 1.0 - kappa
                    f = (p - ks) / (0.5 - ks)
                    if not 0.0 < f < 1.0:
                        ok = False  # a vertex mode has no log-odds
                        break
                    v_log = p * (1.0 - p) / nn / (p - ks) ** 2
                    modes.append((f, v_log / (1.0 - f) ** 2))
                if ok:
                    (f_b, v_b), (f_i, v_i) = modes
                    f_pred = min(
                        f_b / r_mode, 1.0 - 1e-9
                    )  # the level kept: the boundary's share over r
                    d = np.log(f_i / (1.0 - f_i)) - np.log(f_pred / (1.0 - f_pred))
                    v_pair += max(0.0, float(d * d) - (v_b + v_i + 1.0 / n_u[i] + 1.0 / total_b))
            v_level = float(count_logvar(n_u[b]) + count_logvar(n_u[i])) + v_pair
            m_level = level_map_lambda(lam, density_b, a_g[i], n_u[i])
            bound = level_bound_row(lam, density_b, a_g[i], n_u[i], v_level)

            def level_rule(own, held, m=m_level, v=v_level, ub=bound):
                """A level is made from the sender's MEASUREMENT only: its own profile through the
                level-kept map; with none, the crossing total's upper bound. What it holds — an
                imputation — never crosses a level face."""
                return level_row(own, lam, m, v) if own is not None else ub

            rule[(int(b), int(i))] = level_rule

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

                    rule[(int(x), int(b))] = _composed(out7)
                    rule[(int(b), int(x))] = _composed(back7)

        # ── THE LEVEL LANE: the default rule of every directed face that has no composition rule ──
        # A gDNA level is absolute (a profile over u = log(rho / rho_ref)), so it needs no map and no
        # recipient: it crosses the faces composition cannot (strand changes, termini both ways, the
        # AMBIG complex) and the EMPTY node — a piece with no total, 52 % of the ladder's exon pieces —
        # which forwards it unchanged. A full node emits the product of its own level and the priced
        # level it holds; the recipient prices the hop and takes the level as a LOWER bound.
        empty = ~(n_u > 0.0) | ~(a_g > 0.0)
        pure = is_intergenic & (a_g > 0.0)
        # the lane's coordinate: the library's structurally pure gDNA density. It is a REFERENCE for a
        # log axis, so any positive density serves; when the intergenic count is zero (a zero-gDNA
        # library with no intergenic reads) the library-wide density stands in — measured 2026-09-09:
        # with the coordinate at zero the gDNA lane AND the RNA lanes that hang off it were never built,
        # and a both-stranded overlap at zero gDNA read 0.95 gDNA against a truth of 0.
        rho_ref = float(n_u[pure].sum() / a_g[pure].sum()) if a_g[pure].sum() > 0.0 else 0.0
        if not rho_ref > 0.0:
            rho_ref = float(n_u[a_g > 0.0].sum() / max(a_g[a_g > 0.0].sum(), EPS))
        lane = None
        if rho_ref > 0.0:
            u = lam
            gene_edge = np.zeros(n, bool)
            for b in np.flatnonzero(is_bnd):
                lo, hi_ = left[b], right[b]
                gene_edge[b] = (lo >= 0 and is_intergenic[lo]) or (hi_ >= 0 and is_intergenic[hi_])
            own_level: list = [None] * n
            for x in np.flatnonzero(~empty & ~is_intergenic):
                if gene_edge[x]:
                    own_level[x] = poisson_level(u, n_u[x], a_g[x], rho_ref)
                elif own[x] is not None and np.ptp(own[x]) > EPS:
                    own_level[x] = level_of_profile(own[x], lam, u, n_u[x], a_g[x], rho_ref)
            faces = {
                (int(x), int(y))
                for x in range(n)
                for y in (left[x], right[x])
                if y >= 0 and not is_intergenic[y] and (int(x), int(y)) not in rule
            }
            lane = _Lane(u, lam, rho_ref, n_u, a_g, empty, own_level, faces)
        rna = None if lane is None else self._rna_lanes(ctx, lane, own, is_intron, is_intergenic)
        cube = _CubeSite(
            fp & fn,
            n_u,
            a_r,
            int(ctx.n_tilt) if ctx.n_tilt else K,
            {"pos": fp, "neg": fn},
            left,
            right,
        )
        return _PreparedTransfer(own, rule, K, lane, rna, cube)

    # ── THE RNA LEVEL LANES (the both-stranded locus, phase 1; owner rulings 2026-09-08) ─────────
    def _rna_lanes(self, ctx, lane, own, is_intron, is_intergenic):
        """One lane per strand: FACES from the flag bits (strand ``s``'s level crosses a face iff the
        boundary carries none of ``s``'s four bits and both nodes admit ``s``; across ``s``'s OWN
        junction it enters ``s``'s intron — the crossing IS the intron's unspliced population — and
        not ``s``'s exon; a terminus of ``s`` stops ``s`` both ways), TWO-SIDED only between an intron
        of ``s`` and its own boundary (rung 1's one shared unspliced population; ⭐ the intron test is
        PER STRAND — a region that admits ``s`` and carries no exon of ``s`` is ``s``'s intron whatever
        the other strand does there), and SOURCES: a single-strand node's own claim read as its live
        strand's RNA level, and the certified flux at each of an exon's junctions as that strand's level
        at the exon (two junctions charged the pair's own disagreement beyond counting; one, counting
        alone). The coordinate ``rho_ref_s`` is the library's strand-``s`` unspliced density over its
        single-strand exons. Nothing pooled, no constant."""
        n = ctx.n_slots
        is_bnd = np.asarray(ctx.is_boundary, bool)
        is_exon = np.asarray(ctx.is_exon_region, bool)
        fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
        exon_of = {"pos": np.asarray(ctx.exon_pos, bool), "neg": np.asarray(ctx.exon_neg, bool)}
        left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
        flags = np.asarray(ctx.boundary_flags, np.uint16)
        n_u = np.asarray(ctx.n_slot, np.float64)
        a_r = np.asarray(ctx.eff_rna, np.float64)
        cnt = np.asarray(ctx.unspliced_count, np.float64)
        rr = (np.asarray(ctx.route_rate_lo, np.float64), np.asarray(ctx.route_rate_hi, np.float64))
        sc = (np.asarray(ctx.sj_count_lo, np.float64), np.asarray(ctx.sj_count_hi, np.float64))
        empty = ~(n_u > 0.0) | ~(a_r > 0.0)
        single = ~(fp & fn)
        kappa = None if self._strand is None else float(self._strand[0])
        # the split is a witness of a strand's RNA only where the library's strand channel is live —
        # the derived deadband's verdict, already on the context: a κ within its noise of ½ zeroes the
        # strand precision of every single-strand exon (`tau_lam` there is the strand term alone; at an
        # intron the factory's density precision joins it, so introns cannot stand for the channel)
        tau = np.asarray(ctx.own.tau_lam, np.float64)
        split_live = kappa is not None and bool(np.any(tau[is_exon & (fp != fn)] > 0.0))
        lanes = {}
        for name, free, col in (("pos", fp, 0), ("neg", fn, 1)):
            all_bits, _sj_bits, term_bits = strand_bits[name]
            # the lane's WITNESS is the count of the reads strand-``s`` RNA produces: its own genome-strand
            # column when the library reads sense, the other under an antisense protocol (`read_column`;
            # measured 2026-09-09 at κ ≈ 0.01: the own column held a handful of reads at every
            # single-strand exon, every hop was priced as counting on nothing, and the floors arrived
            # blurred to nothing at a zero-gDNA overlap that then read 0.35 gDNA)
            col_read = read_column(col, kappa)
            intron_s = ~is_bnd & free & ~exon_of[name]
            faces, two_sided = set(), set()
            for x in range(n):
                for y in (left[x], right[x]):
                    if y < 0 or is_intergenic[y] or not (free[x] and free[y]):
                        continue
                    b = x if is_bnd[x] else y
                    i = y if is_bnd[x] else x
                    f = int(flags[b])
                    if not (f & all_bits):
                        faces.add((int(x), int(y)))
                        if intron_s[i]:
                            two_sided.add((int(x), int(y)))
                    elif not (f & term_bits) and intron_s[i]:
                        faces.add((int(x), int(y)))
                        two_sided.add((int(x), int(y)))
            sel = is_exon & free & single & (a_r > 0.0)
            rho_ref = (
                float(cnt[sel, col_read].sum() / a_r[sel].sum()) if a_r[sel].sum() > 0.0 else 0.0
            )
            own_level: list = [None] * n
            flux_of: list = [None] * n
            if rho_ref > 0.0:
                for x in np.flatnonzero(~empty & free):
                    parts = []
                    if single[x] and own[x] is not None and np.ptp(own[x]) > EPS:
                        parts.append(
                            rna_level_of_profile(own[x], lane.lam, lane.u, n_u[x], a_r[x], rho_ref)
                        )
                    if is_exon[x]:
                        # the junction's estimate of the exon's RNA abundance, priced by THE NODE PAIR
                        # (owner ruling 2026-09-08): the junction's spliced count at its route rate
                        # against the exon's own count of that strand per RNA opportunity — read on the
                        # genome-strand column strand-s RNA READS on (its own at kappa >= 1/2, the other
                        # under an antisense protocol; the rate is in transcript-strand terms). Both
                        # counts' counting plus their disagreement beyond it; a probed junction beside
                        # an unprobed exon disagrees and goes weak. Kept LOWER-SIDED: the two-sided
                        # estimate over-claimed at the probe cliff (measured 2026-09-08).
                        for b in (left[x], right[x]):
                            if b < 0 or not is_bnd[b]:
                                continue
                            hi = 1 if left[x] == b else 0
                            c_j, r_j = float(sc[hi][b, col]), float(rr[hi][b, col])
                            if not (c_j > 0.0 and r_j > 0.0):
                                continue
                            v = count_price(c_j, c_j / r_j, cnt[x, col_read], a_r[x])
                            fl = flux_level(lane.u, c_j, r_j, rho_ref, v)
                            parts.append(fl)
                            if flux_of[x] is None:
                                flux_of[x] = {}
                            flux_of[x][int(b)] = fl
                    if parts:
                        own_level[x] = intersect(parts)
            lanes[name] = _RnaLane(
                lane.u,
                rho_ref,
                cnt[:, col_read],
                a_r,
                empty,
                own_level,
                frozenset(faces),
                frozenset(two_sided),
                flux_of,
                cnt[:, 1 - col_read] if split_live else None,
            )
        return lanes


class _Lane:
    """The level lane's static data for one sweep: the grid, the reference density, every node's total,
    opportunity, emptiness and own level, and the directed faces the lane serves."""

    __slots__ = ("u", "lam", "rho_ref", "n_u", "a_g", "empty", "own_level", "faces")

    def __init__(self, u, lam, rho_ref, n_u, a_g, empty, own_level, faces):
        self.u, self.lam, self.rho_ref = u, lam, float(rho_ref)
        self.n_u, self.a_g, self.empty = n_u, a_g, empty
        self.own_level, self.faces = own_level, faces

    def emit(self, s: int, held: Message | None) -> Level | None:
        """What ``s`` sends on the lane: an empty node forwards what it holds unchanged; a full node the
        INTERSECTION of its own level's lower side and the (already priced) level it holds. Bounds
        intersect, they do not multiply: a product of one-sided claims sharpens along a chain into a
        hard bound at the noisiest node's mode (measured: nine terminus boundaries with one true gDNA
        fragment each moved 2 → 29 on the ladder's `g05 ss.99 OFF`); the tighter bound winning at each
        density sharpens nothing."""
        far = None if held is None else held.level_gdna
        if self.empty[s]:
            return far
        own = self.own_level[s]
        parts = [
            p
            for p in (
                None if own is None else lower_side(own),
                None if far is None else far.profile,
            )
            if p is not None
        ]
        if not parts:
            return None
        return Level(intersect(parts), float(self.n_u[s]), float(self.a_g[s]))

    def receive(self, level: Level, x: int) -> Level:
        """What a FULL recipient holds: the level priced for this hop and taken as a lower bound, now
        standing at ``x``."""
        v = hop_price(level.n, level.a, self.n_u[x], self.a_g[x])
        return Level(priced_level(level.profile, self.u, v), float(self.n_u[x]), float(self.a_g[x]))

    def row(self, level: Level, x: int):
        """A held level as ``x``'s composition profile (a coordinate change: it was priced on arrival)."""
        return profile_of_level(
            level.profile, self.u, self.lam, self.n_u[x], self.a_g[x], self.rho_ref
        )


class _RnaLane:
    """ONE strand's RNA level lane: the coordinate ``u = log(rho_s / rho_ref_s)`` on the solve grid,
    each node's strand count (the lane's own witness) and RNA opportunity, its emptiness, its own
    level, the directed faces the lane serves and the faces it crosses TWO-SIDED."""

    __slots__ = (
        "u",
        "rho_ref",
        "count",
        "a",
        "empty",
        "own_level",
        "faces",
        "two_sided",
        "flux",
        "other",
    )

    def __init__(
        self, u, rho_ref, count, a, empty, own_level, faces, two_sided, flux=None, other=None
    ):
        self.u, self.rho_ref = u, float(rho_ref)
        self.count, self.a, self.empty = count, a, empty
        self.own_level, self.faces, self.two_sided = own_level, faces, two_sided
        #: per node, ``{boundary: level}`` — the junction's priced estimate of the exon's RNA, kept per
        #: FACE so the solve can tell which face's composition already carries it
        self.flux = [None] * len(own_level) if flux is None else flux
        #: the OTHER genome-strand column's count per node, where the library's strand channel is live
        #: (the derived deadband): the split's asymmetry ``count − other`` is the node's own estimate of
        #: this strand's RNA count, the witness the hop's price compares (`None`: the column count is)
        self.other = other

    def witness(self, y: int):
        """The strand's RNA count at ``y`` and its Poisson variance, read from the column split: the
        asymmetry ``count − other`` (gDNA splits evenly, so it cancels; the other strand's RNA reads on
        the other column) over the protocol's strand contrast ``|1 − 2κ|`` — a factor common to every
        node, so it cancels from every ratio the price takes and is left out. A non-positive asymmetry
        is a DARK node: no RNA of this strand is measurable there."""
        c_r, c_o = float(self.count[y]), float(self.other[y])
        return c_r - c_o, c_r + c_o

    def emit(self, s: int, x: int, far: Level | None) -> Level | None:
        """An empty node forwards; a full node the INTERSECTION of its own level and what it holds —
        the own level WHOLE across a two-sided face (its boundary shares its population exactly), its
        lower side everywhere else. Bounds intersect, they do not multiply."""
        if self.empty[s]:
            return far
        own = self.own_level[s]
        if own is not None and (int(s), int(x)) not in self.two_sided:
            own = lower_side(own)
        parts = [p for p in (own, None if far is None else far.profile) if p is not None]
        if not parts:
            return None
        rna_count = rna_var = None
        if self.other is not None:
            rna_count, rna_var = self.witness(s)
        return Level(intersect(parts), float(self.count[s]), float(self.a[s]), rna_count, rna_var)

    def receive(self, level: Level, s: int, x: int) -> Level:
        """What a FULL recipient holds: across a TWO-SIDED face the whole profile (an intron and its
        own boundary share one unspliced population), everywhere else its lower side — and on EVERY
        face the hop's price: both column counts' counting plus the disagreement, beyond its own
        counting, between the two nodes' estimates of THIS STRAND's abundance. ⭐ The witness is the
        column split's asymmetry (`witness`), not the column count: the column holds gDNA's half, which
        jumps with every probe edge whether or not this strand's RNA is there, while the asymmetry
        is this strand's RNA alone. So a dark recipient (no measurable RNA of the strand) agrees with
        a dark claim and the claim arrives whole — the perfectly dark host intron's "no RNA of mine
        here" that resolves the tilt at an antisense exon's boundaries under capture (the test
        chromosome's probed span loci, ~2,000 fragments a row) — while a lit recipient disagrees with a
        claim from a dimmer node by the cliff between them and blurs it away (the ladder's `g05 ss.99
        ON`, 2026-09-09: a + intron at 0.005 fragments per base whose nascent RNA the − gene's probe
        captured 170-fold, carried whole under a counting-only exemption, read a 93 % RNA junction as
        86 % gDNA). Where the library's strand channel is dead the column count is the witness
        (`count_price`). Measured against the exemption and against `count_price` on the column
        counts, three panels and the ladder, halves apart, both frames."""
        v = float(count_logvar(level.n) + count_logvar(self.count[x]))
        if self.other is None or level.rna_count is None:
            v = count_price(level.n, level.a, self.count[x], self.a[x])
        else:
            n_s, v_s = float(level.rna_count), float(level.rna_count_var)
            n_x, v_x = self.witness(x)
            if n_s > 0.0 and n_x > 0.0:
                r = (n_x / float(self.a[x])) / (n_s / float(level.a))
                v += max(0.0, float(np.log(r)) ** 2 - (v_s / (n_s * n_s) + v_x / (n_x * n_x)))
        p = level.profile if (int(s), int(x)) in self.two_sided else lower_side(level.profile)
        p = blur_row(p, self.u, v) if v > 0.0 else p
        rna_count = rna_var = None
        if self.other is not None:
            rna_count, rna_var = self.witness(x)
        return Level(p, float(self.count[x]), float(self.a[x]), rna_count, rna_var)


_RNA_FIELD = {"pos": "level_rna_pos", "neg": "level_rna_neg"}


class _CubeSite:
    """What the cube delivery reads at a node: the AMBIG mask, every node's total and RNA
    opportunity, and the tilt grid's size."""

    __slots__ = ("ambig", "n_u", "a_r", "n_tilt", "free", "left", "right")

    def __init__(self, ambig, n_u, a_r, n_tilt, free=None, left=None, right=None):
        self.ambig, self.n_u, self.a_r, self.n_tilt = ambig, n_u, a_r, int(n_tilt)
        self.free, self.left, self.right = free, left, right


class _PreparedTransfer:
    """One sweep's working object: every node's own claim, the rules, the gDNA lane, the two RNA
    lanes, the two passes' state."""

    def __init__(
        self,
        own,
        rule: dict,
        n_grid: int,
        lane: _Lane | None = None,
        rna: dict | None = None,
        cube: "_CubeSite | None" = None,
    ):
        self.own = own
        self.rule = rule
        self.lane = lane
        self.rna = rna
        self.cube = cube
        self._K = int(n_grid)
        self.held: dict = {False: None, True: None}

    # ── phase 1: propagate — the recipient's kernel, run by the backbone in chain order ──────────
    def propagate(self, *, backward: bool):
        if self.own is None or not (
            self.rule
            or (self.lane and self.lane.faces)
            or (self.rna and any(r.faces for r in self.rna.values()))
        ):
            return None  # nothing to say anywhere: every node holds SILENCE
        held: list = [None] * len(self.own)
        self.held[bool(backward)] = held
        lane = self.lane
        rna = self.rna

        def receive(s: int, i: int):
            """``s`` sends two things apart — its own claim and what it holds from its far side (written
            by this pass one step earlier) — and the face decides: a composition rule composes and maps
            (MODIFY), passes (FORWARD) or reads the measurement only (a level face); with no rule the
            LANE carries the gDNA level; absent both (a structural pure-gDNA neighbour) STOP."""
            far = held[s]
            comp = None
            fn = self.rule.get((int(s), int(i)))
            if fn is not None:
                out = fn(self.own[s], None if far is None else far.composition)
                if out is not None and np.ptp(out) > EPS:
                    comp = _norm(out)
            level = None
            if lane is not None and (int(s), int(i)) in lane.faces:
                level = lane.emit(s, far)
                if level is not None and not lane.empty[i]:
                    level = lane.receive(level, i)
            fields = {}
            if rna:
                for name, rl in rna.items():
                    if (int(s), int(i)) not in rl.faces:
                        continue
                    lv = rl.emit(s, i, None if far is None else getattr(far, _RNA_FIELD[name]))
                    if lv is not None and not rl.empty[i]:
                        lv = rl.receive(lv, s, i)
                    if lv is not None:
                        fields[_RNA_FIELD[name]] = lv
            if comp is None and level is None and not fields:
                return SILENCE
            msg = Message(composition=comp, level_gdna=level, **fields)
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
        lane = self.lane
        for i in range(n):
            parts, bounds = [], []
            for m in (from_left[i], from_right[i]):
                if m is None:
                    continue
                if m.composition is not None:
                    parts.append(m.composition)
                if m.level_gdna is not None and lane is not None and not lane.empty[i]:
                    bounds.append(lane.row(m.level_gdna, i))
            if bounds:
                parts.append(intersect(bounds))  # two bounds on one density: the tighter wins
            if parts:
                rows[i] = _fuse(parts)
                live = True
        live = self._ceilings(from_left, from_right, rows) or live
        cube = self._cube_rows(from_left, from_right)
        if not live and not cube:
            return PsiMessage.silent()
        return PsiMessage(lam_rows=rows if live else None, cube_rows=cube or None)

    def _ceilings(self, from_left, from_right, rows) -> bool:
        """THE UPPER SIDE AT SINGLE-STRAND NODES (phase 2, 2026-09-08): an RNA level of the node's live
        strand says "at most this much gDNA". Read ONLY from a face that sent no composition — a held
        level on that side, and the node's own junction flux at that face — because a licensed face's
        splice-in map already carries the flux as its cap and a composition already carries its
        sender's witnesses (reading them again counted them twice: the weak-κ zero control 42 → 1,540).
        Bounds intersect; the row joins the node's other witnesses. Returns whether anything was added."""
        rna, site, lane = self.rna, self.cube, self.lane
        if not rna or site is None or site.free is None:
            return False
        added = False
        for name, rl in rna.items():
            field = _RNA_FIELD[name]
            for i in np.flatnonzero(site.free[name] & ~site.ambig & ~rl.empty):
                bounds = []
                for side, m in ((site.left[i], from_left[i]), (site.right[i], from_right[i])):
                    if m is None or m.composition is not None:
                        continue
                    lv = getattr(m, field)
                    if lv is not None:
                        bounds.append(lv.profile)
                    fx = rl.flux[i]
                    if fx is not None and int(side) in fx:
                        bounds.append(fx[int(side)])
                if not bounds:
                    continue
                row = rna_row_of_level(
                    intersect(bounds), rl.u, lane.lam, site.n_u[i], site.a_r[i], rl.rho_ref
                )
                if np.ptp(row) <= EPS:
                    continue
                rows[i] = _fuse([rows[i], row]) if np.ptp(rows[i]) > EPS else row
                added = True
        return added

    def _cube_rows(self, from_left, from_right) -> dict:
        """THE DELIVERY AT AMBIG NODES: the held RNA levels — both sides intersected, plus the node's
        OWN flux level (the spliced claim's one hop, boundary → exon, read at the exon; an AMBIG node
        has no own strand claim, so its own level is the flux alone) — as one row over ψ's cube. The
        tilt needs no lane of its own: both strands' bounds constrain it through the shares."""
        rna, site = self.rna, self.cube
        if not rna or site is None:
            return {}
        n_u, a_r = site.n_u, site.a_r
        pos, neg = rna["pos"], rna["neg"]
        lam = self.lane.lam
        theta = _tilt_grid(site.n_tilt)
        rho = {"pos": pos.rho_ref, "neg": neg.rho_ref}
        out = {}
        for i in np.flatnonzero(site.ambig & ~pos.empty):
            profiles = {}
            for name, rl in (("pos", pos), ("neg", neg)):
                field = _RNA_FIELD[name]
                bounds = [
                    getattr(m, field).profile
                    for m in (from_left[i], from_right[i])
                    if m is not None and getattr(m, field) is not None
                ]
                if rl.own_level[i] is not None:
                    bounds.append(lower_side(rl.own_level[i]))
                if bounds:
                    profiles[name] = intersect(bounds)
            if profiles:
                out[int(i)] = cube_row(profiles, pos.u, lam, theta, n_u[i], a_r[i], rho)
        return out
