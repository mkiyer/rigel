"""Gates for the transfer policy's shape and for the two mechanisms that belong to the policy
rather than to any one face.

The shape: the backbone protocol; the installer and its refusal of an unknown name; the two passes
against an independent recursive reference; the no-echo law; the own claims' liveness under the
strand protocol decision; and the completion contract, that every directed face
carries a rule or a lane face or leads into structural pure gDNA. The mechanisms: phase 2's
ceilings, read only from a face that sent no composition, and the gDNA level lane, where a received
level is a lower bound the hop widens, two bounds intersect rather than multiply, and an empty node
forwards what it holds. The per-face builders are gated by `test_transfer_faces.py`.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

import rigel.calibration.sweep as SW
from rigel.calibration.messages import Policy
from _transfer_harness import (
    _drive,
    _bits,
    _ctx_of,
    _dead_boundaries,
    _drive_the_backbone,
    _full_policy,
    _pairs,
    _prepared,
    _rna,
    _rna_lanes_of,
    _with_alt_splice_sites,
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
    assert all(p._strand is not None for p in seen), (
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
    the message into ``s`` from ``s``'s OTHER neighbour, and the rows are the two messages composed.
    Nothing here reads the policy's pass state; only its claims and rules."""
    from rigel.calibration.messages.transfer_rows import EPS

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
            r = faces.apply(s, i, own[s], m)
            if r is not None and np.ptp(r) > EPS:
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
        from rigel.calibration.messages.transfer_rows import blur_row, hop_price, lower_side

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
                v = hop_price(out[1], out[2], lane.count[i], lane.a[i])
                priced = lower_side(out[0])
                priced = blur_row(priced, lane.u, v) if v > 0.0 else priced
                out = (priced, float(lane.count[i]), float(lane.a[i]))
        lmemo[key] = out
        return out

    n = len(own)
    rows = np.zeros((n, int(ctx.n_grid)))
    for i in range(n):
        parts = [into(s, i) for s in (left[i], right[i]) if s >= 0]
        if lane is not None and not lane.empty[i]:
            from rigel.calibration.messages.transfer_rows import profile_of_level

            bounds = []
            for s in (left[i], right[i]):
                lv = level_into(s, i) if s >= 0 else None
                if lv is not None:
                    bounds.append(
                        profile_of_level(
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
    from rigel.calibration.messages.transfer_rows import EPS

    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    prepared = _prepared(pol, ctx)
    lam = np.linspace(-window, window, n_grid)
    probes = [
        i
        for i in range(ctx.n_slots)
        if prepared.own[i] is not None and np.ptp(prepared.own[i]) > EPS
    ]
    assert probes, "no node carries a claim — this gate would prove nothing"
    checked = 0
    for i in probes[:12]:
        base_rows, fl, br = _drive(prepared, ctx)
        held_before = [t.composition[i] if t.has_composition[i] else None for t in (fl, br)]
        saved = prepared.own[i]
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


# ── phase 2's ceilings: an RNA level as an upper side on a single-strand node's gDNA share ───────


def test_an_rna_level_reads_as_a_ceiling_on_the_gdna_share_and_round_trips():
    """`rna_row_of_level`: a lower-only RNA level (non-decreasing in u) is NON-INCREASING in λ — "at least
    this much RNA" is "at most this much gDNA"; the level of a profile read back is the profile where the
    coordinate resolves. PERTURBATION: the gDNA map in its place does not round-trip."""
    from rigel.calibration.messages.transfer_rows import (
        level_of_profile,
        rna_level_of_profile,
        rna_row_of_level,
    )

    K = 60
    lam = np.linspace(-10.0, 10.0, K)
    u = lam
    n, a_r, rho = 500.0, 100.0, 0.02
    floor = np.maximum.accumulate(-0.5 * ((u - 1.0) / 0.3) ** 2)
    row = rna_row_of_level(floor - floor.max(), u, lam, n, a_r, rho)
    assert np.all(np.diff(row) <= 1e-9) and row[0] > row[-1] + 5.0
    prof = -0.5 * ((lam - 1.5) / 0.7) ** 2
    back = rna_row_of_level(rna_level_of_profile(prof, lam, u, n, a_r, rho), u, lam, n, a_r, rho)
    f_r = 1.0 / (1.0 + np.exp(lam))
    u_of = np.log(f_r * n / (a_r * rho))
    du = np.abs(np.gradient(u_of))
    inside = (u_of > u[0]) & (u_of < u[-1]) & (du >= 0.5 * (u[1] - u[0]))
    tol = 0.15 + 0.02 * np.abs(prof[inside])
    assert np.all(np.abs(back[inside] - prof[inside]) <= tol)
    wrong = rna_row_of_level(level_of_profile(prof, lam, u, n, a_r, rho), u, lam, n, a_r, rho)
    assert not np.all(np.abs(wrong[inside] - prof[inside]) <= tol)


def test_the_ceiling_is_read_only_from_a_face_that_sent_no_composition():
    """`_PreparedTransfer._ceilings` on a hand-built single-strand exon: the left face sent a
    COMPOSITION with an RNA level and a junction flux — nothing of it is read (the face map already
    carries them); the right face sent an RNA level and no composition, and the exon has a flux at that
    junction — both are read, intersected, and delivered as a non-increasing row. An AMBIG node, an
    empty node and a node whose live strand admits nothing get no ceiling. PERTURBATION: with the left
    message's composition removed, its level and flux join the intersection and the row changes."""
    from rigel.calibration.messages import Received
    from rigel.calibration.messages.faces import Faces
    from rigel.calibration.messages.lanes import LevelLane
    from rigel.calibration.messages.transfer import _PreparedTransfer, _SolveSite
    from rigel.calibration.messages.transfer_rows import intersect, rna_row_of_level

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

    flux = {
        (1, 0): floor(0.4),
        (1, 1): floor(0.1),
    }  # exon 1's junctions: slot 0 to its left, 2 right
    none = _bits(5, [])
    pos = LevelLane("pos", u, lam, 0.5, n_u / 2, a_r, empty, [None] * 5, none, total=n_u, flux=flux)
    neg = LevelLane("neg", u, lam, 0.4, n_u / 2, a_r, empty, [None] * 5, none, total=n_u)
    gd = LevelLane("gdna", u, lam, 0.5, n_u, a_r, empty, [None] * 5, none)
    site = _SolveSite(fp & fn, {"pos": fp, "neg": fn})
    prep = _PreparedTransfer(
        [None] * 5, Faces(lam, left, right), {"gdna": gd, "pos": pos, "neg": neg}, site
    )
    comp = -0.5 * ((lam - 1.0) / 0.5) ** 2
    held_l, held_r = floor(0.8), floor(-0.3)
    from_left, from_right = Received.empty(5, K), Received.empty(5, K)
    from_left.has_neighbour[1] = from_right.has_neighbour[1] = True
    from_left.composition[1], from_left.has_composition[1] = comp, True
    from_left.level_rna_pos.write(1, held_l, 25.0, 100.0)
    from_right.level_rna_pos.write(1, held_r, 25.0, 100.0)
    rows = np.zeros((5, K))
    assert prep._ceilings(from_left, from_right, rows)
    want = rna_row_of_level(intersect([held_r, flux[(1, 1)]]), u, lam, 400.0, 100.0, 0.5)
    np.testing.assert_allclose(rows[1], want, atol=1e-12)
    assert np.all(np.diff(rows[1]) <= 1e-9)
    assert not rows[[0, 2, 3, 4]].any()
    # the perturbation: the left composition removed → its level and flux join
    from_left2 = Received.empty(5, K)
    from_left2.has_neighbour[1] = True
    from_left2.level_rna_pos.write(1, held_l, 25.0, 100.0)
    rows2 = np.zeros((5, K))
    assert prep._ceilings(from_left2, from_right, rows2)
    want2 = rna_row_of_level(
        intersect([held_l, flux[(1, 0)], held_r, flux[(1, 1)]]),
        u,
        lam,
        400.0,
        100.0,
        0.5,
    )
    np.testing.assert_allclose(rows2[1], want2, atol=1e-12)
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
    from rigel.calibration.messages.faces import fuse
    from rigel.calibration.messages.transfer_rows import intersect

    for i in range(n):
        parts, bounds = [], []
        for t in (fl, br):
            if t.has_composition[i]:
                parts.append(t.composition[i])
            if t.level_gdna.present[i] and not prepared.lanes["gdna"].empty[i]:
                bounds.append(prepared.lanes["gdna"].row(t.level_gdna.profile[i], i))
        if bounds:
            parts.append(intersect(bounds))
        if parts:
            without[i] = fuse(parts)
    site = prepared.site
    fp, fn = site.free["pos"], site.free["neg"]
    # every side that exists sent a composition
    all_comp = (fl.has_composition | ~fl.has_neighbour) & (br.has_composition | ~br.has_neighbour)
    single = (fp ^ fn) & ~site.ambig
    quiet = single & all_comp
    assert quiet.sum() > 0
    np.testing.assert_allclose(rows[quiet], without[quiet], atol=1e-12)


# ── the gDNA level lane: a level crosses a face as a lower bound the hop widens ──────────────────


def test_a_received_gdna_level_is_a_lower_bound_and_the_hop_widens_it():
    """`LevelLane.receive` on the gDNA lane (no two-sided face): what a full recipient holds is
    non-decreasing in u (a level that crosses a face says "at least this much gDNA" and nothing more),
    a two-sided input loses only its upper side, and a hop across a density cliff — a larger price —
    widens it."""
    from rigel.calibration.messages import Levels
    from rigel.calibration.messages.lanes import LevelLane

    u = np.linspace(-10, 10, 60)
    two_sided = -0.5 * ((u - 1.0) / 0.4) ** 2
    two_sided -= two_sided.max()
    # one density on both sides: the price is the two totals' counting alone
    flat = LevelLane(
        "gdna",
        u,
        u,
        0.05,
        np.array([1.0e6, 1.0e6]),
        np.array([1.0e4, 1.0e4]),
        np.zeros(2, bool),
        [None, None],
        _bits(2, [(0, 1)]),
    )
    held = Levels.empty(2, u.shape[0])
    held.write(1, two_sided, 1.0e6, 1.0e4)  # what the hop wrote at 1 before the recipient priced it
    flat.receive(held, 0, 1)
    lower = held.profile[1].copy()  # a row of the table: copy before the row is rewritten
    assert np.all(np.diff(lower) >= -1e-12), "not a lower bound"
    below = u < 0.0
    np.testing.assert_allclose(lower[below], two_sided[below], atol=1e-3)
    assert np.all(lower[u > 2.0] > -1e-9), "the upper side was kept"
    # a cliff between the two nodes: the discrepancy beyond counting widens the bound
    cliff = LevelLane(
        "gdna",
        u,
        u,
        0.05,
        np.array([1.0e6, 1.0e6]),
        np.array([1.0e4, 4.0e4]),
        np.zeros(2, bool),
        [None, None],
        _bits(2, [(0, 1)]),
    )
    held.write(1, two_sided, 1.0e6, 1.0e4)
    cliff.receive(held, 0, 1)
    wide = held.profile[1]
    assert np.all(np.diff(wide) >= -1e-12)
    assert wide[np.argmin(np.abs(u + 1.0))] > lower[np.argmin(np.abs(u + 1.0))], "no widening"


def test_the_level_coordinates_round_trip_through_a_node_total():
    """`level_of_profile` and `profile_of_level` are one map read both ways: a composition profile
    taken to a level at (n, a) and back at the same (n, a) is itself wherever the level lies below
    the node's total; and the same level read at a node with a LARGER total is a smaller gDNA share
    (the level is kept, the share is not)."""
    from rigel.calibration.messages.transfer_rows import level_of_profile, profile_of_level

    lam = np.linspace(-10, 10, 601)  # a fine grid: the round trip interpolates twice
    u = lam
    rho_ref, n, a = 0.05, 400.0, 1000.0
    row = -0.5 * ((lam - np.log(0.3 / 0.7)) / 1.0) ** 2  # a share near 0.3 → 120 gDNA of 400
    level = level_of_profile(row, lam, u, n, a, rho_ref)
    back = profile_of_level(level, u, lam, n, a, rho_ref)
    core = np.abs(lam - np.log(0.3 / 0.7)) < 2.0
    np.testing.assert_allclose(back[core], (row - row.max())[core], atol=0.05)
    bigger = profile_of_level(level, u, lam, 4.0 * n, a, rho_ref)
    assert lam[np.argmax(bigger)] < lam[np.argmax(back)], "a larger total must read a smaller share"


def test_bounds_intersect_and_do_not_multiply():
    """`intersect`: the pointwise tighter of two lower bounds, never their product — two identical soft
    bounds intersect to themselves (a product would double the penalty), and a tight bound beats a
    loose one at every density."""
    from rigel.calibration.messages.transfer_rows import intersect, lower_side

    u = np.linspace(-10, 10, 60)
    soft = lower_side(-0.5 * ((u - 1.0) / 2.0) ** 2)
    np.testing.assert_allclose(intersect([soft, soft]), soft, atol=1e-12)
    tight = lower_side(-0.5 * ((u - 1.0) / 0.5) ** 2)
    np.testing.assert_allclose(intersect([soft, tight]), tight, atol=1e-12)
    both = intersect([lower_side(-0.5 * ((u - 3.0) / 1.0) ** 2), soft])
    assert np.all(both <= soft + 1e-12) and np.all(np.diff(both) >= -1e-12)


def test_a_full_node_emits_the_intersection_of_its_own_lower_side_and_what_it_holds():
    """The ratchet gate: on a hand-built lane (three full nodes in a row, which the toy does not have),
    the level a full node emits is the pointwise tighter of its own lower side and the priced level it
    holds — never their sum, which sharpened a chain of nine one-fragment boundaries into a hard bound
    on the ladder. PERTURBATION: an `emit` that multiplies fails here."""
    from rigel.calibration.messages import Levels
    from rigel.calibration.messages.lanes import LevelLane
    from rigel.calibration.messages.transfer_rows import intersect, lower_side

    u = np.linspace(-10, 10, 60)
    own = [None, -0.5 * ((u + 1.0) / 1.5) ** 2, None]
    lane = LevelLane(
        "gdna",
        u,
        u,
        0.05,
        np.array([100.0, 120.0, 90.0]),
        np.array([50.0, 60.0, 45.0]),
        np.zeros(3, bool),
        own,
        _bits(3, [(0, 1), (1, 2)]),
    )
    held = Levels.empty(3, u.shape[0])
    held.write(1, lower_side(-0.5 * ((u - 0.5) / 0.8) ** 2), 100.0, 50.0)  # what node 1 holds
    assert lane.emit(1, 2, held)
    sent = held.profile[2]
    want = intersect([lower_side(own[1]), held.profile[1]])
    np.testing.assert_allclose(sent, want, atol=1e-12)
    product = lower_side(own[1]) + held.profile[1]
    assert np.max(np.abs(sent - (product - product.max()))) > 0.5, "the sum: the ratchet"
    assert (held.count[2], held.opportunity[2]) == (120.0, 60.0)


def test_the_hop_price_is_both_countings_and_the_discrepancy_beyond_them():
    from scipy.special import polygamma

    from rigel.calibration.messages.transfer_rows import hop_price

    counting = float(polygamma(1, 50.5) + polygamma(1, 200.5))
    assert abs(hop_price(50, 100.0, 200, 400.0) - counting) < 1e-12, (
        "equal densities: counting only"
    )
    d = np.log(8.0) ** 2 - (1 / 50 + 1 / 200)
    assert abs(hop_price(50, 100.0, 200, 50.0) - (counting + d)) < 1e-12, "the excess over counting"
    assert hop_price(50, 100.0, 200, 50.0) == hop_price(50, 100.0, 200, 50.0)


def test_a_composition_arrives_only_through_a_face_with_a_composition_rule(sweep_inputs):
    """The table's ``has_composition`` against the face table, independently of the kernel that wrote
    it: on the live toy a node holds a composition from a side only where the directed face into it
    carries a composition rule, and a node that heard something through a face with NO rule heard a
    level. PERTURBATION: a kernel that marks a level as a composition fails here."""
    from rigel.calibration.messages.faces import NONE

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
    from rigel.calibration.messages.transfer_rows import blur_row, hop_price, lower_side

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
            v = hop_price(held.count[e], held.opportunity[e], lane.count[x], lane.a[x])
            want = lower_side(held.profile[e])
            want = blur_row(want, lane.u, v) if v > 0.0 else want
            np.testing.assert_allclose(held.profile[x], want, atol=1e-9)
            assert (held.count[x], held.opportunity[x]) == (float(lane.count[x]), float(lane.a[x]))
            checked += 1
    assert checked > 0, "no level crossed an empty node on the toy: the gate proved nothing"


def test_the_face_table_holds_one_of_five_kinds_with_finite_parameters_at_real_faces(sweep_inputs):
    """`Faces`, the rules as typed tables: on the live toy every rule the builders wrote is one of the
    five kinds, its scalar parameters and its rows are finite, `at` reads back what the table holds,
    and every ruled face is a real face (the source is the destination's neighbour on that side).
    PERTURBATION: a rule at a face that does not exist is refused, so a builder cannot address the
    wrong neighbour; a second rule at a face is refused, so the builders' faces stay disjoint."""
    import pytest

    from rigel.calibration.messages.faces import (
        EDGE,
        FORWARD,
        LEVEL,
        NONE,
        RULE_NAMES,
        SPLICE_OUT,
        TRANSPORT,
        Faces,
    )

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
        if f.kind == LEVEL:
            assert f.row is not None and f.row2 is not None and f.var > 0.0
    assert {FORWARD, TRANSPORT, SPLICE_OUT} <= seen, [RULE_NAMES[k] for k in sorted(seen)]
    assert len(RULE_NAMES) == 6
    # a rule can only be written at a face that exists
    n = int(ctx.n_slots)
    table = Faces(np.linspace(-1.0, 1.0, int(ctx.n_grid)), left, right)
    i = int(np.flatnonzero(left >= 0)[0])
    table.set(int(left[i]), i, FORWARD)
    assert table.kind_at(int(left[i]), i) == FORWARD and table.has(int(left[i]), i)
    stranger = int(left[i]) - 1 if int(left[i]) > 0 else (i + 2 if i + 2 < n else i - 2)
    with pytest.raises(ValueError, match="no face into"):
        table.set(stranger, i, FORWARD)
    # and a face carries one rule: the builders' faces are disjoint, so a second write is a builder
    # addressing another's face, refused rather than silently winning or losing by order
    with pytest.raises(ValueError, match="already carries"):
        table.set(int(left[i]), i, TRANSPORT)
