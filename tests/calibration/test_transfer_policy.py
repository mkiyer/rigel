"""Gates for the transfer policy's shape and for the two mechanisms that belong to the policy
rather than to any one face.

The shape: the backbone protocol; the installer and its refusal of an unknown name; the native pass
wired and read without a copy; the two passes against an independent recursive reference (single hops
of the native pass composed recursively); the no-echo law; the own claims' liveness under the
strand protocol decision; and the completion contract, that every directed face
carries a rule or a lane face or leads into structural pure gDNA. The mechanisms: the solve's
ceilings, read only from a face that sent no composition, and the gDNA level lane, where a received
level is a lower bound the hop widens, two bounds intersect rather than multiply, and an empty node
forwards what it holds. The per-face builders are gated by `test_transfer_faces.py`.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

import rigel.calibration.sweep as SW
from rigel.calibration.blocks import SweepCapture
from rigel.calibration.messages import Policy
from rigel.calibration.messages.silent import SilentPolicy
from rigel.native import transfer_rows as R
from _transfer_harness import (
    EDGE,
    FORWARD,
    LEVEL,
    NONE,
    RULE_NAMES,
    SPLICE_OUT,
    TRANSPORT,
    Faces,
    LevelLane,
    Prepared,
    Received,
    RowTable,
    _drive,
    _bits,
    _ctx_of,
    _dead_boundaries,
    _drive_the_backbone,
    _full_policy,
    _hop,
    _lane_prepared,
    _leaves,
    _pairs,
    _prepared,
    _rna,
    _rna_lanes_of,
    _rule,
    _with_alt_splice_sites,
    fuse,
    intersect,
)


def test_the_transfer_policy_satisfies_the_backbone_protocol():
    from rigel.calibration.messages.transfer import TransferPolicy

    assert isinstance(TransferPolicy(), Policy)


def test_the_policy_name_installs_the_transfer_policy(sweep_inputs):
    """`message_policy = "transfer"` must install `TransferPolicy` (an unreadable knob is worse
    than no knob), and the unknown-name refusal must survive the new branch."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.config import CalibrationConfig

    calibrate_mod = sys.modules["rigel.calibration.calibrate"]
    seen: list = []
    orig = SW.solve_chain

    def spy(*a, **kw):
        seen.append(kw.get("policy"))
        return orig(*a, **kw)

    calibrate_mod.solve_chain = spy
    try:
        calibrate_mod.calibrate(
            payload=sweep_inputs["payload"],
            config=_dc.replace(CalibrationConfig(), message_policy="transfer"),
            **sweep_inputs["calibrate_kw"],
        )
    finally:
        calibrate_mod.solve_chain = orig
    assert seen and all(isinstance(p, TransferPolicy) for p in seen)
    assert all(p.strand is not None for p in seen), (
        "the exon -> boundary message must be ON in production"
    )
    with pytest.raises(ValueError, match="unknown message_policy"):
        calibrate_mod.calibrate(
            payload=sweep_inputs["payload"],
            config=_dc.replace(CalibrationConfig(), message_policy="no-such-policy"),
            **sweep_inputs["calibrate_kw"],
        )


def _reference_rows(prepared, ctx):
    """An independent implementation of the two passes, recursive rather than sequential: the message
    into ``i`` from its neighbour ``s`` is the face's rule applied to ``s``'s own claim composed with
    the message into ``s`` from ``s``'s OTHER neighbour (one hop of the native pass, `_rule`, on a table
    holding only that composed message), and the rows are the two messages composed. Nothing here reads
    the policy's pass state; only its claims, rules and single hops."""
    left, right = list(ctx.left), list(ctx.right)
    own, faces = prepared.own, prepared.faces
    memo = {}

    def norm(r):
        return r - r.max()

    def into(s, i):
        key = (s, i)
        if key in memo:
            return memo[key]
        out = None
        if faces.has(s, i):
            far = left[s] if right[s] == i else right[s]
            m = into(far, s) if far >= 0 else None
            r = _rule(prepared, s, i, held=m)
            if r is not None and np.ptp(r) > R.EPS:
                out = norm(r)
        memo[key] = out
        return out

    lane = prepared.lanes.get("gdna")
    lmemo = {}

    def level_into(s, i):
        """The lane, recursively: nothing unless the face is a lane face; an EMPTY sender forwards the
        level that reaches it from its far side unchanged; a full sender sends the product of its own
        level and the level it holds; a full recipient prices the hop (both totals' counting plus the
        discrepancy beyond it) and takes the level as a lower bound."""
        key = (s, i)
        if key in lmemo:
            return lmemo[key]
        out = None
        if lane is not None and lane.serves(s, i):
            far = left[s] if right[s] == i else right[s]
            held = level_into(far, s) if far >= 0 else None
            if lane.empty[s]:
                out = held
            else:
                own_s = lane.own_level[s]
                if own_s is not None:
                    own_s = np.maximum.accumulate(own_s)  # its lower side
                parts = [q for q in (own_s, None if held is None else held[0]) if q is not None]
                if parts:
                    tight = parts[0] if len(parts) == 1 else np.minimum(parts[0], parts[1])
                    out = (norm(tight), float(lane.count[s]), float(lane.a[s]))
            if out is not None and not lane.empty[i]:
                v = R.hop_price(out[1], out[2], lane.count[i], lane.a[i])
                priced = R.lower_side(out[0])
                priced = R.blur_row(priced, lane.u, v) if v > 0.0 else priced
                out = (priced, float(lane.count[i]), float(lane.a[i]))
        lmemo[key] = out
        return out

    n = len(own)
    rows = np.zeros((n, int(ctx.n_grid)))
    for i in range(n):
        parts = [into(s, i) for s in (left[i], right[i]) if s >= 0]
        if lane is not None and not lane.empty[i]:
            bounds = []
            for s in (left[i], right[i]):
                lv = level_into(s, i) if s >= 0 else None
                if lv is not None:
                    bounds.append(
                        R.profile_of_level(
                            lv[0], lane.u, lane.lam, lane.total[i], lane.a[i], lane.rho_ref
                        )
                    )
            if bounds:  # two bounds on one density intersect: the pointwise tighter
                parts.append(
                    norm(bounds[0] if len(bounds) == 1 else np.minimum(bounds[0], bounds[1]))
                )
        parts = [p for p in parts if p is not None]
        if parts:
            rows[i] = norm(sum(parts))
    return rows


def test_the_backbone_rows_equal_an_independent_recursive_reference_of_the_passes(sweep_inputs):
    """The passes, gated against a second implementation: what the backbone's two sequential passes
    deliver equals the recursive definition — the message into a node is the face's rule applied to
    the sender's own claim composed with what reached the sender from ITS far side — on the live toy
    and on the toy with populated inside pieces and alternative splice sites (every rule family live)."""
    pol, _p, _g, _w = _full_policy(sweep_inputs)
    for ctx in (_ctx_of(sweep_inputs), _with_alt_splice_sites(_ctx_of(sweep_inputs))):
        prepared = _prepared(pol, ctx)
        rows = _drive_the_backbone(prepared, ctx)
        assert rows.any(), "the toy delivered nothing — this gate would prove nothing"
        np.testing.assert_allclose(rows, _reference_rows(prepared, ctx), rtol=0, atol=1e-10)


def test_PERTURBATION_no_node_ever_hears_its_own_claim_back(sweep_inputs):
    """The no-echo law, watched: replace ONE node's own claim by a distinctive profile and re-run the
    passes — what that node HOLDS from either side must not move (its claim never returns to it),
    while some other node's rows must (the claim did travel). Checked at an intron with two live faces
    and at a strand-live exon."""
    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    prepared = _prepared(pol, ctx)
    lam = np.linspace(-window, window, n_grid)
    probes = [
        i
        for i in range(ctx.n_slots)
        if prepared.own[i] is not None and np.ptp(prepared.own[i]) > R.EPS
    ]
    assert probes, "no node carries a claim — this gate would prove nothing"
    checked = 0
    for i in probes[:12]:
        base_rows, fl, br = _drive(prepared, ctx)
        held_before = [t.composition[i] if t.has_composition[i] else None for t in (fl, br)]
        saved = prepared.own[i].copy()  # a row of the table: copy before the row is rewritten
        prepared.own[i] = -0.5 * ((lam - 3.3) / 0.2) ** 2  # a spike nowhere near any real claim
        rows, fl, br = _drive(prepared, ctx)
        held_after = [t.composition[i] if t.has_composition[i] else None for t in (fl, br)]
        prepared.own[i] = saved
        for a, b in zip(held_before, held_after):
            if a is None or b is None:
                assert a is None and b is None
            else:
                np.testing.assert_array_equal(a, b, err_msg=f"slot {i} heard its own claim back")
        moved = [
            j for j in range(ctx.n_slots) if j != i and not np.array_equal(rows[j], base_rows[j])
        ]
        if moved:
            checked += 1
    assert checked >= 1, "no perturbed claim travelled anywhere — the gate could not have fired"


def test_a_dead_strand_channel_carries_no_own_claim(sweep_inputs):
    """The vacuity law, for every claim carried by the strand channel: where a node's strand channel
    is dead (``own.tau_lam == 0``), that node's OWN CLAIM is
    absent — an exon's and a boundary's alike — so nothing of its own can travel; and a policy built
    without strand parameters carries no strand claim anywhere."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer import TransferPolicy

    pol, _rows, _g, _w = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    n = int(ctx.n_slots)
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    live = _prepared(pol, ctx)
    assert any(live.own[e] is not None for e in np.flatnonzero(is_exon)), "no live exon claim"
    dead = _prepared(pol, _dc.replace(ctx, has_own_composition=np.zeros(n, bool)))
    for x in range(n):
        if is_exon[x]:
            assert dead.own[x] is None, f"a dead exon {x} carries a claim"
        elif is_bnd[x] and dead.own[x] is not None:
            assert not dead.own[x].any(), f"a dead boundary {x} carries a strand claim"
    no_strand = _prepared(TransferPolicy(), ctx)
    for x in range(n):
        if is_exon[x]:
            assert no_strand.own[x] is None
    # boundaries dead, exons live: the boundary's own claim is absent at every boundary
    dead_b = _prepared(pol, _dead_boundaries(ctx))
    for b in np.flatnonzero(is_bnd):
        assert dead_b.own[b] is None or not dead_b.own[b].any()


def test_every_directed_face_is_served_or_faces_structural_pure_gdna(sweep_inputs):
    """The completion contract, structurally: after `prepare`, every directed face (x → y) of the chain
    carries a composition rule or is a lane face, unless an end of it is an intergenic region — a
    TERMINAL: structurally pure gDNA, nothing to impute there and no own level to send, so it neither
    receives nor sends and the chain breaks at it. STOP by omission is impossible by construction."""
    ctx = _ctx_of(sweep_inputs)
    prepared = _prepared(_full_policy(sweep_inputs)[0], ctx)
    lane = prepared.lanes.get("gdna")
    assert lane is not None and lane.face.any()
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    intergenic = ~is_bnd & ~is_exon & ~fp & ~fn
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    unserved = [
        (int(x), int(y))
        for x in range(int(ctx.n_slots))
        for y in (left[x], right[x])
        if y >= 0
        and not intergenic[y]
        and not intergenic[x]
        and not prepared.faces.has(x, y)
        and not lane.serves(x, y)
    ]
    assert not unserved, f"faces with no rule and no lane: {unserved[:8]}"
    # and a lane face never doubles a composition rule (no witness counted twice)
    ruled, laned = set(prepared.faces.pairs()), _pairs(lane.face, ctx)
    assert not (ruled & laned)
    # a terminal is at neither end of any served face: it is where the chain breaks
    assert intergenic.any(), "the toy has no intergenic region — this clause would prove nothing"
    served = ruled | laned
    assert not [f for f in served if intergenic[f[0]] or intergenic[f[1]]]


# ── the solve's ceilings: an RNA level as an upper side on a single-strand node's gDNA share ─────


def test_an_rna_level_reads_as_a_ceiling_on_the_gdna_share_and_round_trips():
    """`rna_row_of_level`: a lower-only RNA level (non-decreasing in u) is NON-INCREASING in λ — "at least
    this much RNA" is "at most this much gDNA"; the level of a profile read back is the profile where the
    coordinate resolves. PERTURBATION: the gDNA map in its place does not round-trip."""

    K = 60
    lam = np.linspace(-10.0, 10.0, K)
    u = lam
    n, a_r, rho = 500.0, 100.0, 0.02
    floor = np.maximum.accumulate(-0.5 * ((u - 1.0) / 0.3) ** 2)
    row = R.rna_row_of_level(floor - floor.max(), u, lam, n, a_r, rho)
    assert np.all(np.diff(row) <= 1e-9) and row[0] > row[-1] + 5.0
    prof = -0.5 * ((lam - 1.5) / 0.7) ** 2
    back = R.rna_row_of_level(
        R.rna_level_of_profile(prof, lam, u, n, a_r, rho), u, lam, n, a_r, rho
    )
    f_r = 1.0 / (1.0 + np.exp(lam))
    u_of = np.log(f_r * n / (a_r * rho))
    du = np.abs(np.gradient(u_of))
    inside = (u_of > u[0]) & (u_of < u[-1]) & (du >= 0.5 * (u[1] - u[0]))
    tol = 0.15 + 0.02 * np.abs(prof[inside])
    assert np.all(np.abs(back[inside] - prof[inside]) <= tol)
    wrong = R.rna_row_of_level(R.level_of_profile(prof, lam, u, n, a_r, rho), u, lam, n, a_r, rho)
    assert not np.all(np.abs(wrong[inside] - prof[inside]) <= tol)


def test_the_ceiling_is_read_only_from_a_face_that_sent_no_composition():
    """The policy's solve on a hand-built single-strand exon — THE CEILING: the left face sent a
    COMPOSITION with an RNA level and a junction flux — nothing of it is read (the face map already
    carries them); the right face sent an RNA level and no composition, and the exon has a flux at that
    junction — both are read, intersected, and delivered as a non-increasing row. An AMBIG node, an
    empty node and a node whose live strand admits nothing get no ceiling. PERTURBATION: with the left
    message's composition removed, its level and flux join the intersection and the row changes."""
    K = 41
    lam = np.linspace(-6.0, 6.0, K)
    u = lam
    #             0: bnd   1: exon(+)   2: bnd   3: ambig   4: empty exon(+)
    n_u = np.array([50.0, 400.0, 50.0, 300.0, 0.0])
    a_r = np.full(5, 100.0)
    empty = ~(n_u > 0)
    fp = np.array([True, True, True, True, True])
    fn = np.array([False, False, False, True, False])
    left = np.array([-1, 0, 1, 2, 3])
    right = np.array([1, 2, 3, 4, -1])

    def floor(u_b):
        return -0.5 * np.maximum(0.0, (u_b - u) / 0.2) ** 2

    flux = RowTable((5, 2), K)  # exon 1's junctions: slot 0 to its left, 2 right
    flux[1, 0], flux[1, 1] = floor(0.4), floor(0.1)
    none = _bits(5, [])
    no_levels = RowTable(5, K)
    pos = LevelLane("pos", u, lam, 0.5, n_u / 2, a_r, empty, no_levels, none, total=n_u, flux=flux)
    neg = LevelLane("neg", u, lam, 0.4, n_u / 2, a_r, empty, no_levels, none, total=n_u)
    gd = LevelLane("gdna", u, lam, 0.5, n_u, a_r, empty, no_levels, none)
    prep = Prepared.hand_built(
        RowTable(5, K), Faces(lam, left, right), {"gdna": gd, "pos": pos, "neg": neg}, fp, fn
    )
    comp = -0.5 * ((lam - 1.0) / 0.5) ** 2
    held_l, held_r = floor(0.8), floor(-0.3)
    from_left, from_right = Received.empty(5, K), Received.empty(5, K)
    from_left.has_neighbour[1] = from_right.has_neighbour[1] = True
    from_left.composition[1], from_left.has_composition[1] = comp, True
    from_left.level_rna_pos.write(1, held_l, 25.0, 100.0)
    from_right.level_rna_pos.write(1, held_r, 25.0, 100.0)
    rows = prep.solve(from_left, from_right).lam_rows
    assert rows is not None
    ceiling = R.rna_row_of_level(intersect([held_r, flux[(1, 1)]]), u, lam, 400.0, 100.0, 0.5)
    assert np.all(np.diff(ceiling) <= 1e-9)
    # the delivered row is the held composition fused with the ceiling — the left face's level and
    # flux are NOT in it
    np.testing.assert_allclose(rows[1], fuse([comp - comp.max(), ceiling]), atol=1e-12)
    assert not rows[[0, 2, 3, 4]].any()
    # the perturbation: the left composition removed → its level and flux join
    from_left2 = Received.empty(5, K)
    from_left2.has_neighbour[1] = True
    from_left2.level_rna_pos.write(1, held_l, 25.0, 100.0)
    rows2 = prep.solve(from_left2, from_right).lam_rows
    assert rows2 is not None
    want2 = R.rna_row_of_level(
        intersect([held_l, flux[(1, 0)], held_r, flux[(1, 1)]]),
        u,
        lam,
        400.0,
        100.0,
        0.5,
    )
    np.testing.assert_allclose(rows2[1], want2, atol=1e-12)
    assert np.all(np.diff(rows2[1]) <= 1e-9)  # a ceiling alone: non-increasing in the gDNA share
    assert not np.allclose(rows2[1], rows[1])


def test_the_flux_is_kept_per_face_and_a_licensed_face_keeps_the_ceiling_out(sweep_inputs):
    """On the live toy: every exon with a junction flux carries it PER FACE in the lane; and at every
    single-strand exon whose faces all sent a composition the delivered row equals the row without the
    ceiling step (nothing added), so the flux inside a face map is never counted twice."""
    ctx, prepared = _rna_lanes_of(sweep_inputs)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    sc = np.asarray(ctx.sj_count_lo) + np.asarray(ctx.sj_count_hi)
    per_face = 0
    for name, lane in _rna(prepared).items():
        for x in np.flatnonzero(is_exon):
            fx = [
                (b, lane.flux_at(x, side))
                for side, b in ((0, int(ctx.left[x])), (1, int(ctx.right[x])))
                if lane.flux_at(x, side) is not None
            ]
            if not fx:
                continue
            for b, prof in fx:
                assert b >= 0 and sc[b].sum() > 0
                assert np.all(np.diff(prof) >= -1e-12)  # a lower bound
                per_face += 1
    assert per_face > 0
    rows, fl, br = _drive(prepared, ctx)
    without = np.zeros_like(rows)
    n = len(prepared.own)
    gd = prepared.lanes["gdna"]
    for i in range(n):
        parts, bounds = [], []
        for t in (fl, br):
            if t.has_composition[i]:
                parts.append(t.composition[i])
            if t.level_gdna.present[i] and not gd.empty[i]:
                bounds.append(
                    R.profile_of_level(
                        t.level_gdna.profile[i], gd.u, gd.lam, gd.total[i], gd.a[i], gd.rho_ref
                    )
                )
        if bounds:
            parts.append(intersect(bounds))
        if parts:
            without[i] = fuse(parts)
    fp, fn = prepared.free_pos, prepared.free_neg
    # every side that exists sent a composition
    all_comp = (fl.has_composition | ~fl.has_neighbour) & (br.has_composition | ~br.has_neighbour)
    single = (fp ^ fn) & ~prepared.ambig
    quiet = single & all_comp
    assert quiet.sum() > 0
    np.testing.assert_allclose(rows[quiet], without[quiet], atol=1e-12)


# ── the gDNA level lane: a level crosses a face as a lower bound the hop widens ──────────────────


def test_a_received_gdna_level_is_a_lower_bound_and_the_hop_widens_it():
    """The gDNA lane's hop into a full recipient (no two-sided face): what the recipient holds is
    non-decreasing in u (a level that crosses a face says "at least this much gDNA" and nothing more), a
    two-sided input loses only its upper side, and a hop across a density cliff — a larger price — widens
    it. The sender is EMPTY and holds the level with its witness, so it forwards and the recipient prices."""
    u = np.linspace(-10, 10, 60)
    two_sided = -0.5 * ((u - 1.0) / 0.4) ** 2
    two_sided -= two_sided.max()

    def received_after(a):
        lane = LevelLane(
            "gdna",
            u,
            u,
            0.05,
            np.array([1.0e6, 1.0e6]),
            a,
            np.array([True, False]),
            RowTable(2, u.shape[0]),
            _bits(2, [(0, 1)]),
        )
        held = Received.empty(2, u.shape[0])
        held.level_gdna.write(0, two_sided, 1.0e6, 1.0e4)  # what the sender holds, with its witness
        _hop(_lane_prepared(lane, [-1, 0], [1, -1]), 0, 1, held)
        assert held.level_gdna.present[1]
        return held.level_gdna.profile[1].copy()

    # one density on both sides: the price is the two totals' counting alone
    lower = received_after(np.array([1.0e4, 1.0e4]))
    assert np.all(np.diff(lower) >= -1e-12), "not a lower bound"
    below = u < 0.0
    np.testing.assert_allclose(lower[below], two_sided[below], atol=1e-3)
    assert np.all(lower[u > 2.0] > -1e-9), "the upper side was kept"
    # a cliff between the two nodes: the discrepancy beyond counting widens the bound
    wide = received_after(np.array([1.0e4, 4.0e4]))
    assert np.all(np.diff(wide) >= -1e-12)
    assert wide[np.argmin(np.abs(u + 1.0))] > lower[np.argmin(np.abs(u + 1.0))], "no widening"


def test_the_level_coordinates_round_trip_through_a_node_total():
    """`level_of_profile` and `profile_of_level` are one map read both ways: a composition profile
    taken to a level at (n, a) and back at the same (n, a) is itself wherever the level lies below
    the node's total; and the same level read at a node with a LARGER total is a smaller gDNA share
    (the level is kept, the share is not)."""

    lam = np.linspace(-10, 10, 601)  # a fine grid: the round trip interpolates twice
    u = lam
    rho_ref, n, a = 0.05, 400.0, 1000.0
    row = -0.5 * ((lam - np.log(0.3 / 0.7)) / 1.0) ** 2  # a share near 0.3 → 120 gDNA of 400
    level = R.level_of_profile(row, lam, u, n, a, rho_ref)
    back = R.profile_of_level(level, u, lam, n, a, rho_ref)
    core = np.abs(lam - np.log(0.3 / 0.7)) < 2.0
    np.testing.assert_allclose(back[core], (row - row.max())[core], atol=0.05)
    bigger = R.profile_of_level(level, u, lam, 4.0 * n, a, rho_ref)
    assert lam[np.argmax(bigger)] < lam[np.argmax(back)], "a larger total must read a smaller share"


def test_bounds_intersect_and_do_not_multiply():
    """`intersect`: the pointwise tighter of two lower bounds, never their product — two identical soft
    bounds intersect to themselves (a product would double the penalty), and a tight bound beats a
    loose one at every density."""

    u = np.linspace(-10, 10, 60)
    soft = R.lower_side(-0.5 * ((u - 1.0) / 2.0) ** 2)
    np.testing.assert_allclose(intersect([soft, soft]), soft, atol=1e-12)
    tight = R.lower_side(-0.5 * ((u - 1.0) / 0.5) ** 2)
    np.testing.assert_allclose(intersect([soft, tight]), tight, atol=1e-12)
    both = intersect([R.lower_side(-0.5 * ((u - 3.0) / 1.0) ** 2), soft])
    assert np.all(both <= soft + 1e-12) and np.all(np.diff(both) >= -1e-12)


def test_a_full_node_emits_the_intersection_of_its_own_lower_side_and_what_it_holds():
    """The ratchet gate: on a hand-built lane (three full nodes in a row, which the toy does not have),
    the level a full node emits is the pointwise tighter of its own lower side and the priced level it
    holds — never their sum, which sharpened a chain of nine one-fragment boundaries into a hard bound
    on the ladder. The recipient is EMPTY, so the hop stops at the emission and it is read as sent.
    PERTURBATION: a pass whose full node multiplies fails here."""
    u = np.linspace(-10, 10, 60)
    own = RowTable(3, u.shape[0])
    own[1] = -0.5 * ((u + 1.0) / 1.5) ** 2
    lane = LevelLane(
        "gdna",
        u,
        u,
        0.05,
        np.array([100.0, 120.0, 90.0]),
        np.array([50.0, 60.0, 45.0]),
        np.array([False, False, True]),
        own,
        _bits(3, [(0, 1), (1, 2)]),
    )
    held = Received.empty(3, u.shape[0])
    held.level_gdna.write(
        1, R.lower_side(-0.5 * ((u - 0.5) / 0.8) ** 2), 100.0, 50.0
    )  # node 1 holds
    _hop(_lane_prepared(lane, [-1, 0, 1], [1, 2, -1]), 1, 2, held)
    assert held.level_gdna.present[2]
    sent = held.level_gdna.profile[2]
    want = intersect([R.lower_side(own[1]), held.level_gdna.profile[1]])
    np.testing.assert_allclose(sent, want, atol=1e-12)
    product = R.lower_side(own[1]) + held.level_gdna.profile[1]
    assert np.max(np.abs(sent - (product - product.max()))) > 0.5, "the sum: the ratchet"
    assert (held.level_gdna.count[2], held.level_gdna.opportunity[2]) == (120.0, 60.0)


def test_the_hop_price_is_both_countings_and_the_discrepancy_beyond_them():
    from scipy.special import polygamma

    counting = float(polygamma(1, 50.5) + polygamma(1, 200.5))
    assert abs(R.hop_price(50, 100.0, 200, 400.0) - counting) < 1e-12, (
        "equal densities: counting only"
    )
    d = np.log(8.0) ** 2 - (1 / 50 + 1 / 200)
    assert abs(R.hop_price(50, 100.0, 200, 50.0) - (counting + d)) < 1e-12, (
        "the excess over counting"
    )
    assert R.hop_price(50, 100.0, 200, 50.0) == R.hop_price(50, 100.0, 200, 50.0)


def test_a_composition_arrives_only_through_a_face_with_a_composition_rule(sweep_inputs):
    """The table's ``has_composition`` against the face table, independently of the kernel that wrote
    it: on the live toy a node holds a composition from a side only where the directed face into it
    carries a composition rule, and a node that heard something through a face with NO rule heard a
    level. PERTURBATION: a kernel that marks a level as a composition fails here."""
    ctx = _ctx_of(sweep_inputs)
    prepared = _prepared(_full_policy(sweep_inputs)[0], ctx)
    _rows, fl, br = _drive(prepared, ctx)
    kind = prepared.faces.kind
    for side, t in ((0, fl), (1, br)):
        ruled = kind[:, side] != NONE
        assert not (t.has_composition & ~ruled).any(), (
            "a composition arrived through a face with no rule"
        )
        assert (t.heard & ~ruled & t.has_level).sum() == (t.heard & ~ruled).sum()
        assert (t.has_composition & ruled).any(), (
            "no composition crossed a ruled face: the gate is vacuous"
        )
        assert (t.heard & ~ruled).any(), "no level crossed an unruled face: the gate is vacuous"


def test_an_empty_node_forwards_the_level_it_holds_unchanged(sweep_inputs):
    """The empty node is transparent: the toy's inside pieces have no total, so what the boundary beyond
    such a piece holds from it must be exactly what the piece received — same profile, same (n, a) of
    the last full node — priced only at the full recipient. PERTURBATION: a policy whose empty nodes
    price the hop fails this gate."""

    ctx = _ctx_of(sweep_inputs)
    prepared = _prepared(_full_policy(sweep_inputs)[0], ctx)
    _rows, fl, br = _drive(prepared, ctx)
    lane = prepared.lanes["gdna"]
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    checked = 0
    for backward, nbr, held in ((False, left, fl.level_gdna), (True, right, br.level_gdna)):
        for e in np.flatnonzero(lane.empty):
            s = int(nbr[e])
            if s < 0 or not held.present[e]:
                continue
            x = int(right[e] if not backward else left[e])
            if x < 0 or not held.present[x] or lane.empty[x]:
                continue
            v = R.hop_price(held.count[e], held.opportunity[e], lane.count[x], lane.a[x])
            want = R.lower_side(held.profile[e])
            want = R.blur_row(want, lane.u, v) if v > 0.0 else want
            np.testing.assert_allclose(held.profile[x], want, atol=1e-9)
            assert (held.count[x], held.opportunity[x]) == (float(lane.count[x]), float(lane.a[x]))
            checked += 1
    assert checked > 0, "no level crossed an empty node on the toy: the gate proved nothing"


def test_the_face_table_holds_one_of_five_kinds_with_finite_parameters_at_real_faces(sweep_inputs):
    """`Faces`, the rules as typed tables: on the live toy every rule the builders wrote is one of the
    five kinds, its scalar parameters and its rows are finite, `at` reads back what the table holds,
    and every ruled face is a real face (the source is the destination's neighbour on that side).
    A rule at a face that does not exist, or a second rule at a face, is refused by the kernel's own table
    writer (`transfer_kernel.h`, ``FacesOut::set``), so a builder cannot address the wrong neighbour and
    the builders' faces stay disjoint."""
    ctx = _ctx_of(sweep_inputs)
    prepared = _prepared(_full_policy(sweep_inputs)[0], ctx)
    faces = prepared.faces
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    kinds = {FORWARD, TRANSPORT, SPLICE_OUT, EDGE, LEVEL}
    seen = set()
    for s, i in faces.pairs():
        assert s == (left[i] if s < i else right[i]), f"({s}, {i}) is not a face"
        f = faces.at(s, i)
        assert f.kind in kinds and f.kind != NONE
        seen.add(f.kind)
        for name in ("n_u", "n_s", "a_b", "a_x", "width", "var"):
            assert np.isfinite(getattr(f, name))
        for row in (f.row, f.row2):
            assert row is None or (row.shape == (int(ctx.n_grid),) and np.isfinite(row).all())
        if f.kind == TRANSPORT:
            assert f.row is not None and f.n_u > 0.0
    assert {FORWARD, TRANSPORT, SPLICE_OUT} <= seen, [RULE_NAMES[k] for k in sorted(seen)]
    assert len(RULE_NAMES) == 6


def test_the_layer_reads_a_row_only_under_its_mask(sweep_inputs):
    """THE TABLES ARE ALLOCATED UNFILLED, so this is the reader audit made executable on the kernel's own
    tables: every matrix the builders returned is poisoned with NaN where its mask says no row — the
    claims, each lane's own levels (the face rows and the flux store hold written rows only) — and both
    received tables' compositions and level profiles are poisoned whole before the passes; the pass's
    tables (every bit, count and witness, every PRESENT row) and the delivered channels are bit-identical
    to the unpoisoned run. One reader of an ABSENT row — in the pass or the solve — and NaN reaches a
    number. The block pipeline's arena holds the same unfilled matrices and reads them through the same
    code (`native/solve_kernel.cpp`), which the replay holds bit-identical. PERTURBATION: a kernel reading
    a row without its bit (`lane_emit` taking the own level unmasked; the solve reading a composition
    without its bit) fails here."""
    pol, _rows, _n_grid, _window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    n, K = int(ctx.n_slots), int(ctx.n_grid)

    def run(poison):
        prepared = _prepared(pol, ctx)
        if poison:
            prepared.own.rows[~prepared.own.mask] = np.nan
            for ln in prepared.lanes.values():
                ln.own_level.rows[~ln.own_level.mask] = np.nan
        order = np.arange(n, dtype=np.int64)
        tables = []
        for nbr, seq, backward in ((ctx.left, order, False), (ctx.right, order[::-1], True)):
            nbr = np.asarray(nbr, np.int64)
            received = Received.empty(n, K)
            if poison:
                received.composition.fill(np.nan)
                for lane in Received.LANES:
                    getattr(received, lane).profile.fill(np.nan)
            received.has_neighbour[seq] = nbr[seq] >= 0
            prepared.run_pass(received, seq, nbr, np.zeros(n, bool), backward=backward)
            tables.append(received)
        return tables, prepared.solve(*tables)

    plain, poisoned = run(False), run(True)
    assert plain[0][0].heard.any() and plain[1].lam_rows is not None, "the toy must carry messages"
    for a, b in zip(plain[0], poisoned[0]):
        for name, x, y in _leaves(a, b):
            if x.dtype != bool and x.ndim == 2:
                assert np.isfinite(x).all(), f"{name}: the plain run holds a non-finite row"
            assert np.array_equal(x, y, equal_nan=True), f"{name} moved under poison"
    rows, prows = plain[1].lam_rows, poisoned[1].lam_rows
    assert np.isfinite(rows).all() and np.array_equal(rows, prows), "the delivered λ rows moved"
    # the toy delivers no AMBIG cube rows, so the audit covers the λ rows alone
    assert plain[1].cube_rows is None and poisoned[1].cube_rows is None


def test_the_backbone_runs_the_sweep_in_one_native_call(sweep_inputs, monkeypatch):
    """The wiring: `sweep.solve_chain` reaches `native.solve_blocks` ONCE per sweep, whatever the policy —
    a spy counts the calls and reads the policy the kernel was told — and the transfer sweep's capture
    holds received tables that heard something, while the silent sweep runs no layer and holds none."""
    calls = []
    real = SW.solve_blocks

    def spy(**kw):
        calls.append(int(kw["policy"]))
        return real(**kw)

    monkeypatch.setattr(SW, "solve_blocks", spy)
    pol = _full_policy(sweep_inputs)[0]
    cap = SweepCapture()
    SW.solve_chain(*sweep_inputs["args"], **sweep_inputs["kw"], policy=pol, _capture=cap)
    assert calls == [1] and cap.from_left is not None
    assert Received.from_kernel(cap.from_left).heard.any() and cap.policy_name == "transfer"
    quiet = SweepCapture()
    SW.solve_chain(
        *sweep_inputs["args"], **sweep_inputs["kw"], policy=SilentPolicy(), _capture=quiet
    )
    assert calls == [1, 0] and quiet.from_left is None and quiet.policy_name == "silent"
