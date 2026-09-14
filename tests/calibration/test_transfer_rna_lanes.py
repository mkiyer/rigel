"""The RNA level lanes (`lanes.rna_lanes`, one `lanes.LevelLane` per strand) and the cube
delivery — the machinery the both-stranded locus needs.

Each gate was written against the prototype first and watched firing on its perturbation. They
cover the faces derived from the flag bits per strand, the level coordinate and its round trip, the
hop priced on the strand's own counts, the flux level, which nodes are sources, the no-echo law,
the two-sided hop at the pair's price, the split witness that reads a node's estimate of one
strand's RNA from its column asymmetry, a one-sided profile surviving the map onto the cube, the
bracket theorem (three lower bounds plus the strand equation give a two-sided gDNA share), and the
cube row as the intersected held levels plus the node's own flux.
"""

from __future__ import annotations

import numpy as np
from scipy.special import polygamma

from _transfer_harness import (
    _passes,
    _bits,
    _held_rna,
    _pairs,
    _prepared,
    _rna,
    _rna_lanes_of,
    _strand_intron,
    _two_sided,
)


def test_the_rna_faces_come_from_the_flag_bits_per_strand(sweep_inputs):
    """A strand's level crosses a face iff the boundary carries none of that strand's four bits and both
    nodes admit the strand; across the strand's OWN junction it enters the strand's intron (two-sided:
    the crossing IS the intron's unspliced population) and not its exon; a terminus of the strand stops
    it both ways; every two-sided face is an intron ↔ own-boundary face. PERTURBATION: with the junction
    bits masked out of the flags the derivation opens the exon face at every junction."""
    from rigel.calibration.messages.transfer_rows import strand_bits

    ctx, prepared = _rna_lanes_of(sweep_inputs)
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    checked = 0
    for name, lane in _rna(prepared).items():
        all_bits, sj_bits, term_bits = strand_bits[name]
        free = np.asarray(ctx.free_pos if name == "pos" else ctx.free_neg, bool)
        intron_s = _strand_intron(ctx, name)
        for x, y in _pairs(lane.face, ctx):
            assert free[x] and free[y], (name, x, y)
        for x, y in _pairs(lane.two_sided, ctx):
            i = y if is_bnd[x] else x
            assert intron_s[i] and lane.serves(x, y), (name, x, y)
        for b in np.flatnonzero(is_bnd):
            lo, hi = left[b], right[b]
            if lo < 0 or hi < 0:
                continue
            f = int(flags[b])
            for e, i in ((lo, hi), (hi, lo)):
                if not (is_exon[e] and intron_s[i]):
                    continue
                if (f & sj_bits) and not (f & term_bits):
                    checked += 1
                    assert not lane.serves(b, e) and not lane.serves(e, b)
                    assert _two_sided(lane, i, b) and _two_sided(lane, b, i)
                if f & term_bits:
                    assert not lane.serves(b, e) and not lane.serves(b, i)
    assert checked > 0, "no junction face on the toy — this gate would prove nothing"
    # the perturbation: junction bits masked → the exon faces at junctions open
    import dataclasses

    masked = np.asarray(
        flags & ~np.uint16(strand_bits["pos"][1] | strand_bits["neg"][1]), np.uint16
    )
    _c2, p2 = _rna_lanes_of(sweep_inputs, dataclasses.replace(ctx, boundary_flags=masked))
    opened = sum(
        1
        for name, lane in _rna(p2).items()
        for b in np.flatnonzero(is_bnd)
        for e in (left[b], right[b])
        if e >= 0
        and is_exon[e]
        and (int(flags[b]) & strand_bits[name][1])
        and not (int(flags[b]) & strand_bits[name][2])
        and lane.serves(b, e)
    )
    assert opened > 0


def test_the_rna_coordinate_round_trips_where_it_resolves():
    """`rna_level_of_profile` then its inverse reproduces the profile wherever one λ cell moves the level
    by at least half a grid cell (near ``f_r → 1`` the level saturates at the total and several λ cells
    share one u cell — the gDNA level's limit at ``f_g → 1``). PERTURBATION: the gDNA map in its place
    does not round-trip."""
    from rigel.calibration.messages.transfer_rows import level_of_profile, rna_level_of_profile

    K = 60
    lam = np.linspace(-10.0, 10.0, K)
    u = lam
    row = -0.5 * ((lam - 1.5) / 0.7) ** 2
    n, a_r, rho = 500.0, 100.0, 0.02
    f_r = 1.0 / (1.0 + np.exp(lam))
    u_of_lam = np.log(f_r * n / (a_r * rho))
    du = np.abs(np.gradient(u_of_lam))
    inside = (u_of_lam > u[0]) & (u_of_lam < u[-1]) & (du >= 0.5 * (u[1] - u[0]))
    assert inside.sum() >= K // 3
    tol = 0.15 + 0.02 * np.abs(row[inside])  # two linear interpolations, deeper in the tail

    def back(level):
        b = np.interp(u_of_lam, u, level)
        return b - b.max()

    good = back(rna_level_of_profile(row, lam, u, n, a_r, rho))
    assert np.all(np.abs(good[inside] - row[inside]) <= tol)
    wrong = back(level_of_profile(row, lam, u, n, a_r, rho))
    assert not np.all(np.abs(wrong[inside] - row[inside]) <= tol)


def test_the_rna_hop_is_priced_by_the_strands_own_counts():
    """`hop_price` on an RNA lane's own witness: two nodes with equal totals but strand counts 30 and 3
    pay the strand counts' price, not the totals'; a zero count on either side pays counting alone and
    stays finite."""
    from rigel.calibration.messages.transfer_rows import hop_price

    v_strand = hop_price(30.0, 100.0, 3.0, 100.0)
    v_total = hop_price(3000.0, 100.0, 3000.0, 100.0)
    assert v_strand > 10.0 * v_total
    v0 = hop_price(30.0, 100.0, 0.0, 100.0)
    assert np.isfinite(v0) and abs(v0 - float(polygamma(1, 30.5) + polygamma(1, 0.5))) < 1e-12


def test_the_flux_level_is_a_lower_bound_priced_by_the_node_pair():
    """`flux_level`: non-decreasing in u; unpriced, below the rate it IS the spliced count's Poisson
    likelihood on the rate's own opportunity; priced by the junction-to-exon disagreement it stays a
    lower bound and is WIDER (the truth pays less when the pair disagrees); a zero count claims nothing.
    `read_column`: the exon count a junction is priced against is the strand's own column when the
    library reads sense and the other column under an antisense protocol (PERTURBATION: the wrong column
    reads a 90 % transcript as 10 % and prices the floor away)."""
    from rigel.calibration.messages.transfer_rows import (
        flux_level,
        hop_price,
        poisson_level,
        read_column,
    )

    u = np.linspace(-10.0, 10.0, 60)
    rho = 0.02
    fl = flux_level(u, 40.0, 40.0 / 2000.0, rho)
    assert np.all(np.diff(fl) >= -1e-12)
    pl = poisson_level(u, 40.0, 2000.0, rho)
    below = u <= u[np.argmax(pl)]
    np.testing.assert_allclose(fl[below], pl[below] - pl.max())
    assert flux_level(u, 0.0, 0.0, rho) is None
    # a probed junction beside an unprobed exon: the exon's strand count per opportunity sits far below
    # the junction's rate, the pair disagrees, the floor widens and stays a floor
    v_match = hop_price(300.0, 200.0, 1500.0, 1000.0)  # 1.5/bp against 1.5/bp
    v_cliff = hop_price(300.0, 200.0, 150.0, 1000.0)  # 1.5/bp against 0.15/bp
    assert v_cliff > 10.0 * v_match
    wide = flux_level(u, 300.0, 1.5, rho, v_cliff)
    tight = flux_level(u, 300.0, 1.5, rho, v_match)
    assert np.all(np.diff(wide) >= -1e-12)
    at_half = np.argmin(np.abs(u - (np.log(1.5 / rho) - np.log(2.0))))  # half the rate
    assert wide[at_half] > tight[at_half] + 1.0
    # the column follows the protocol
    assert read_column(0, 0.99) == 0 and read_column(1, 0.99) == 1 and read_column(0, None) == 0
    assert read_column(0, 0.01) == 1 and read_column(1, 0.01) == 0
    # the perturbation: at kappa = 0.01 a 90 % − transcript's reads sit on the + column; pricing its
    # junction against the − column (10 %) disagrees 9x and prices the floor away
    v_right = hop_price(180.0, 215.0, 1345.0, 1785.0)
    v_wrong = hop_price(180.0, 215.0, 173.0, 1785.0)
    assert v_wrong > 4.0 and v_right < 0.5


def test_the_rna_sources_are_single_strand_claims_and_the_flux_at_the_exon_only(sweep_inputs):
    """On the live toy: an RNA level exists only at a node that admits the strand and is not empty; a
    single-strand node with a live own claim has one; a junction boundary with flux and no own claim has
    NONE (the spliced claim is one hop, boundary → exon); and some exon with flux carries the flux
    level (its lower side is a floor at the route rate)."""
    ctx, prepared = _rna_lanes_of(sweep_inputs)
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    sc = np.asarray(ctx.sj_count_lo) + np.asarray(ctx.sj_count_hi)
    n_src = 0
    for name, lane in _rna(prepared).items():
        free = np.asarray(ctx.free_pos if name == "pos" else ctx.free_neg, bool)
        for x in range(int(ctx.n_slots)):
            lv = lane.own_level[x]
            if lv is not None:
                n_src += 1
                assert free[x]
                # an EMPTY source is an exon piece beside a lit junction, its level the flux's
                # alone, travelling with the flux's witness
                if lane.empty[x]:
                    assert is_exon[x] and (lane.flux_row[x] >= 0).any()
                    assert lane.flux_witness[x] is not None
                    assert lane.flux_witness[x][0] > 0.0
        for b in np.flatnonzero(is_bnd):
            if sc[b].sum() > 0 and prepared.own[b] is None:
                assert lane.own_level[b] is None, b
    assert n_src > 0
    fluxed = [e for e in np.flatnonzero(is_exon) if sc[[ctx.left[e], ctx.right[e]]].sum() > 0]
    assert fluxed and any(
        prepared.lanes[k].own_level[e] is not None for e in fluxed for k in ("pos", "neg")
    )


def test_PERTURBATION_an_rna_level_never_returns_to_its_source(sweep_inputs):
    """The no-echo law on the RNA lanes: sharpen one node's own RNA level to a hard false floor; what
    that node HOLDS from either side must not move, while some neighbouring slot's held level must."""
    ctx, prepared = _rna_lanes_of(sweep_inputs)

    def passes(prep):
        return _passes(prep, ctx)

    lane = prepared.lanes["pos"]
    before = passes(prepared)
    cands = [
        x
        for x in range(int(ctx.n_slots))
        if lane.own_level[x] is not None
        and any(_held_rna(h, x, "level_rna_pos") is not None for h in before)
    ]
    assert cands, "no + node both sourcing and holding an RNA level on the toy"
    x0 = cands[len(cands) // 2]
    was = [_held_rna(h, x0, "level_rna_pos") for h in before]
    _c2, p2 = _rna_lanes_of(sweep_inputs)
    p2.lanes["pos"].own_level[x0] = np.where(lane.u > 3.0, 0.0, -50.0)
    after = passes(p2)
    now = [_held_rna(h, x0, "level_rna_pos") for h in after]
    for a, b in zip(was, now):
        assert (a is None and b is None) or (
            a is not None and b is not None and np.array_equal(a, b)
        )
    moved = 0
    for y in range(max(0, x0 - 4), min(int(ctx.n_slots), x0 + 5)):
        if y == x0:
            continue
        for h1, h2 in zip(before, after):
            a, b = _held_rna(h1, y, "level_rna_pos"), _held_rna(h2, y, "level_rna_pos")
            if (a is None) != (b is None) or (a is not None and not np.array_equal(a, b)):
                moved += 1
    assert moved > 0, (
        "the sharpened claim reached nobody — the lane is dead and this gate proves nothing"
    )


def test_the_two_sided_hop_keeps_the_whole_profile_and_pays_the_pairs_price():
    """`LevelLane.receive` on an RNA lane: across a two-sided face the profile keeps its upper side;
    across any other face it is lower-sided. BOTH pay the pair's price (`hop_price`: both counts' counting plus
    the disagreement of the strand's two count densities beyond it) — no face is exempt. Where the
    pair agrees the price is counting alone (the identity hop's natural price); where it disagrees
    by a cliff the upper side is blurred away — under a counting-only exemption a lit intron's sharp
    upper side crosses a probe cliff unpriced and reads a mostly-RNA junction as mostly gDNA."""
    from rigel.calibration.messages import Levels
    from rigel.calibration.messages.lanes import LevelLane
    from rigel.calibration.messages.transfer_rows import blur_row, hop_price, lower_side

    u = np.linspace(-10.0, 10.0, 60)
    prof = -0.5 * ((u - 0.0) / 0.4) ** 2

    def received(lane, count, opportunity):
        """Row 1 as the hop wrote it — ``prof`` with the sender's witness — then priced by ``lane``."""
        held = Levels.empty(2, u.shape[0])
        held.write(1, prof, count, opportunity)
        lane.receive(held, 0, 1)
        return held.profile[1]

    # the pair AGREES (one density on both sides): the price is counting alone and the upper side stands
    count = np.array([30.0, 300.0])
    a = np.array([300.0, 3000.0])
    empty, none = np.zeros(2, bool), [None, None]
    one_face = _bits(2, [(0, 1)])
    lane = LevelLane("pos", u, u, 0.5, count, a, empty, none, one_face, two_sided=one_face)
    two = received(lane, 30.0, 300.0)
    v_agree = float(polygamma(1, 30.5) + polygamma(1, 300.5))
    assert abs(hop_price(30.0, 300.0, 300.0, 3000.0) - v_agree) < 1e-12
    np.testing.assert_allclose(two, blur_row(prof, u, v_agree), atol=1e-9)
    assert two.max() - two[-1] > 1.0  # the upper side survives
    # the pair DISAGREES by a cliff (the strand's density 800x higher at the recipient): the whole
    # profile still crosses, but at the pair's price, not counting alone
    count = np.array([3.0, 229.0])
    a = np.array([2200.0, 200.0])
    cliff = LevelLane("pos", u, u, 0.5, count, a, empty, none, one_face, two_sided=one_face)
    two_cliff = received(cliff, 3.0, 2200.0)
    v_pair = hop_price(3.0, 2200.0, 229.0, 200.0)
    v_counting = float(polygamma(1, 3.5) + polygamma(1, 229.5))
    assert v_pair > v_counting + 10.0  # log(800)^2 ~ 45 nats^2 beyond counting
    np.testing.assert_allclose(two_cliff, blur_row(prof, u, v_pair), atol=1e-9)
    assert not np.allclose(two_cliff, blur_row(prof, u, v_counting), atol=1e-3), (
        "the exemption is back: the two-sided face charged counting alone across a cliff"
    )
    # any other face: the lower side at the same price
    open_lane = LevelLane("pos", u, u, 0.5, count, a, empty, none, one_face)
    one = received(open_lane, 3.0, 2200.0)
    np.testing.assert_allclose(one, blur_row(lower_side(prof), u, v_pair), atol=1e-9)
    assert np.all(np.diff(one) >= -1e-9)


def test_the_rna_hop_witness_is_the_split_where_the_strand_channel_is_live():
    """`LevelLane.receive` with the other column given (the library's strand channel live): the price's
    disagreement term reads each node's estimate of THIS strand's RNA from its column split's asymmetry
    (`witness`), not from the column count. (i) Two DARK nodes (no asymmetry) agree whatever their
    column densities do — the whole profile crosses a two-sided face at counting alone, where the
    column witness would have blurred it by a 100-fold cliff; (ii) two LIT nodes at a cliff disagree
    by the ratio of their asymmetries per opportunity, beyond the asymmetries' own Poisson counting;
    (iii) the witness travels with the level across an EMPTY node, so the next full recipient prices
    against the last full node and not the empty; (iv) with no other column the column count is the
    witness (`hop_price` on the column counts). Falsified by making the column count the witness again (i fires) and by
    the counting-only exemption (ii fires)."""
    from rigel.calibration.messages import Levels
    from rigel.calibration.messages.lanes import LevelLane
    from rigel.calibration.messages.transfer_rows import blur_row, hop_price, lower_side

    u = np.linspace(-10.0, 10.0, 60)
    prof = -0.5 * ((u - 0.0) / 0.4) ** 2
    # (i) dark → dark across a 100-fold column-density cliff: counting alone, the whole profile
    count = np.array([5.0, 367.0])  # this strand's read column
    other = np.array([8.0, 400.0])  # the other column: no asymmetry on either node
    a = np.array([277.0, 223.0])
    empty, none = np.zeros(2, bool), [None, None]
    one_face = _bits(2, [(0, 1)])
    lane = LevelLane(
        "pos", u, u, 0.5, count, a, empty, none, one_face, two_sided=one_face, other=other
    )
    held = Levels.empty(2, u.shape[0])
    held.write(1, prof, 5.0, 277.0, 5.0 - 8.0, 13.0)  # the sent level, with the sender's witness
    lane.receive(held, 0, 1)
    got = held.profile[1]
    v_counting = float(polygamma(1, 5.5) + polygamma(1, 367.5))
    np.testing.assert_allclose(got, blur_row(prof, u, v_counting), atol=1e-9)
    assert got.max() - got[-1] > 1.0, "the dark claim's upper side must survive"
    v_column = hop_price(5.0, 277.0, 367.0, 223.0)
    assert v_column > v_counting + 5.0
    assert not np.allclose(got, blur_row(prof, u, v_column), atol=1e-3), (
        "the column count is the witness again: two dark nodes were charged a cliff"
    )
    assert (
        held.has_witness[1]
        and held.rna_count[1] == 367.0 - 400.0
        and held.rna_count_var[1] == 767.0
    )
    # (ii) lit → lit at a cliff: the asymmetries' ratio per opportunity, beyond their counting
    count = np.array([69.0, 198.0])
    other = np.array([0.0, 10.0])
    a = np.array([14054.0, 227.0])
    lit = LevelLane(
        "pos", u, u, 0.5, count, a, empty, none, one_face, two_sided=one_face, other=other
    )
    held.write(1, prof, 69.0, 14054.0, 69.0, 69.0)
    lit.receive(held, 0, 1)
    got = held.profile[1]
    n_s, v_s, n_x, v_x = 69.0, 69.0, 188.0, 208.0
    r = (n_x / 227.0) / (n_s / 14054.0)
    v_lit = float(polygamma(1, 69.5) + polygamma(1, 198.5)) + max(
        0.0, np.log(r) ** 2 - (v_s / n_s**2 + v_x / n_x**2)
    )
    assert v_lit > 20.0  # log(170)^2 ~ 26 nats^2
    np.testing.assert_allclose(got, blur_row(prof, u, v_lit), atol=1e-9)
    assert not np.allclose(
        got, blur_row(prof, u, float(polygamma(1, 69.5) + polygamma(1, 198.5))), atol=1e-3
    ), "the exemption is back: a lit claim crossed a cliff at counting alone"
    # (iii) the witness rides across an EMPTY node
    count = np.array([69.0, 0.0, 198.0])
    other = np.array([0.0, 0.0, 10.0])
    a = np.array([14054.0, 50.0, 227.0])
    empty = np.array([False, True, False])
    chain = LevelLane(
        "pos",
        u,
        u,
        0.5,
        count,
        a,
        empty,
        [prof, None, None],
        _bits(3, [(0, 1), (1, 2)]),
        other=other,
    )
    held = Levels.empty(3, u.shape[0])
    assert chain.emit(0, 1, held)
    assert held.rna_count[1] == 69.0 and held.rna_count_var[1] == 69.0
    assert chain.emit(1, 2, held)
    assert np.array_equal(held.profile[2], held.profile[1]) and (
        held.count[2],
        held.opportunity[2],
        held.rna_count[2],
        held.rna_count_var[2],
    ) == (held.count[1], held.opportunity[1], 69.0, 69.0), (
        "an empty node forwards the level and its witness unchanged"
    )
    chain.receive(held, 1, 2)
    np.testing.assert_allclose(held.profile[2], blur_row(lower_side(prof), u, v_lit), atol=1e-9)
    # (iv) no other column: the column count is the witness
    dead = LevelLane(
        "pos",
        u,
        u,
        0.5,
        np.array([69.0, 198.0]),
        np.array([14054.0, 227.0]),
        np.zeros(2, bool),
        [None, None],
        one_face,
        two_sided=one_face,
    )
    held = Levels.empty(2, u.shape[0])
    held.write(1, prof, 69.0, 14054.0)
    dead.receive(held, 0, 1)
    np.testing.assert_allclose(
        held.profile[1], blur_row(prof, u, hop_price(69.0, 14054.0, 198.0, 227.0)), atol=1e-9
    )
    assert not held.has_witness[1]


def test_a_lower_only_profile_stays_one_sided_on_the_cube():
    """`CubeRow.at`: a lower-only RNA+ profile evaluated on the ``(λ, θ)`` cube is non-decreasing in θ
    (the + share rises with τ) and non-increasing in λ (it falls with the gDNA share) — one-sided through
    the map, no parametric summary. PERTURBATION: a two-sided profile is not monotone."""
    from rigel.calibration.simplex_logodds import CubeRow

    lam = np.linspace(-10.0, 10.0, 60)
    tau = np.sin(np.linspace(-0.5 * np.pi, 0.5 * np.pi, 60))
    fg = 1.0 / (1.0 + np.exp(-lam))
    u = lam
    floor = -0.5 * np.maximum(0.0, (0.0 - u) / 0.3) ** 2
    row = CubeRow(floor, None, u, 400.0, 100.0, 0.5).at(fg, tau)
    assert np.all(np.diff(row, axis=1) >= -1e-9) and np.all(np.diff(row, axis=0) <= 1e-9)
    two = CubeRow(-0.5 * (u / 0.3) ** 2, None, u, 400.0, 100.0, 0.5).at(fg, tau)
    assert not (np.all(np.diff(two, axis=1) >= -1e-9) and np.all(np.diff(two, axis=0) <= 1e-9))


def test_THE_BRACKET_THEOREM_three_lower_bounds_and_the_strand_equation_bracket_the_gdna_share():
    """On a hand-built AMBIG node (truth ``f_g`` 0.5, ``f_+`` 0.3, ``f_−`` 0.2, n = 400): three lower
    bounds — the gDNA level and both RNA levels, each a floor at its truth — plus the node's own
    strand counts give a two-sided gDNA share (a 90 % interval narrower than 0.15 that contains the
    truth), on stranded (κ = 0.99) and unstranded (κ = 0.5) data alike. Removing the gDNA bound
    opens the lower side and removing either RNA bound opens the upper side."""
    import rigel.calibration.simplex_logodds as sl
    from rigel.calibration.messages.transfer_rows import profile_of_level

    K = 60
    lam = np.linspace(-10.0, 10.0, K)
    u = lam
    theta = np.linspace(-0.5 * np.pi, 0.5 * np.pi, K)
    fg = 1.0 / (1.0 + np.exp(-lam))
    n, a_g, a_r = 400.0, 100.0, 100.0
    rho = {"g": 0.5, "pos": 0.5, "neg": 0.5}
    truth = dict(g=0.5, pos=0.3, neg=0.2)

    def floor(u_b, w=0.1):
        return -0.5 * np.maximum(0.0, (u_b - u) / w) ** 2

    def interval(kappa, drop=()):
        p_plus = 0.5 * truth["g"] + kappa * truth["pos"] + (1 - kappa) * truth["neg"]
        u_pos = round(p_plus * n)
        F = np.float32
        tau = np.sin(theta)
        fpk = ((1 - fg)[:, None] * (1 + tau)[None, :] / 2).astype(F)
        fnk = ((1 - fg)[:, None] * (1 - tau)[None, :] / 2).astype(F)
        psi = sl._mixture_strand_loglik(
            np.asarray([u_pos], F)[:, None, None],
            np.asarray([n], F)[:, None, None],
            fg.astype(F)[None, :, None],
            fpk[None],
            fnk[None],
            kappa,
            0.0,
            0.0,
            np.asarray([0.5], F)[:, None, None],
            np.asarray([0.25], F)[:, None, None],
            np.asarray([0.25], F)[:, None, None],
        )[0].astype(np.float64)
        psi += np.asarray(sl._gdna_arm(lam, None) + sl._rna_arm(lam), np.float64).reshape(-1)[
            :, None
        ]
        profiles = {
            s: floor(np.log(truth[s] * n / a_r / rho[s])) for s in ("pos", "neg") if s not in drop
        }
        if profiles:
            psi += sl.CubeRow(profiles.get("pos"), profiles.get("neg"), u, n, a_r, rho["pos"]).at(
                fg, tau
            )
        if "g" not in drop:
            lvl = floor(np.log(truth["g"] * n / a_g / rho["g"]))
            psi += profile_of_level(lvl, u, lam, n, a_g, rho["g"])[:, None]
        post = np.exp(psi - psi.max()).sum(axis=1)
        cdf = np.cumsum(post / post.sum())
        return float(fg[np.searchsorted(cdf, 0.05)]), float(
            fg[min(np.searchsorted(cdf, 0.95), K - 1)]
        )

    for kappa in (0.99, 0.5):
        lo, hi = interval(kappa)
        assert hi - lo < 0.15 and lo <= 0.5 <= hi + 0.05, (kappa, lo, hi)
        lo_g, _hi_g = interval(kappa, drop=("g",))
        assert lo - lo_g >= 0.05, (kappa, lo, lo_g)
        for s in ("pos", "neg"):
            _lo_s, hi_s = interval(kappa, drop=(s,))
            assert hi_s - hi >= 0.05, (kappa, s, hi, hi_s)


def test_the_cube_delivery_is_the_intersected_held_levels_and_the_own_flux(sweep_inputs):
    """`_PreparedTransfer._cube_rows` on a hand-built AMBIG node: the delivered `CubeRow` carries, per
    strand, the intersection of the two held levels and the node's own (flux) level's lower side, with the
    node's total and RNA opportunity and the lanes' reference densities; a non-AMBIG node, an empty node
    and a node holding nothing deliver no row."""
    from rigel.calibration.messages import Received
    from rigel.calibration.messages.faces import Faces
    from rigel.calibration.messages.lanes import LevelLane
    from rigel.calibration.messages.transfer import _PreparedTransfer, _SolveSite
    from rigel.calibration.messages.transfer_rows import intersect, lower_side

    K = 41
    lam = np.linspace(-6.0, 6.0, K)
    u = lam
    n_u = np.array([50.0, 400.0, 0.0, 300.0])
    a_r = np.array([100.0, 100.0, 100.0, 100.0])
    empty = ~(n_u > 0)
    own_pos = [
        None,
        -0.5 * np.maximum(0.0, (1.0 - u) / 0.3) ** 2 - 0.5 * (u / 5.0) ** 2,
        None,
        None,
    ]
    none = _bits(4, [])
    pos = LevelLane("pos", u, lam, 0.5, n_u / 2, a_r, empty, own_pos, none, total=n_u)
    neg = LevelLane("neg", u, lam, 0.4, n_u / 2, a_r, empty, [None] * 4, none, total=n_u)
    gd = LevelLane("gdna", u, lam, 0.5, n_u, a_r, empty, [None] * 4, none)
    ambig = np.array([False, True, True, True])
    free = {"pos": np.ones(4, bool), "neg": ambig}
    left, right = np.array([-1, 0, 1, 2]), np.array([1, 2, 3, -1])
    site = _SolveSite(ambig, free)
    prep = _PreparedTransfer(
        [None] * 4, Faces(lam, left, right), {"gdna": gd, "pos": pos, "neg": neg}, site
    )
    lv_l = -0.5 * np.maximum(0.0, (0.5 - u) / 0.2) ** 2
    lv_r = -0.5 * np.maximum(0.0, (0.0 - u) / 0.2) ** 2
    lv_n = -0.5 * np.maximum(0.0, (-1.0 - u) / 0.2) ** 2
    from_left, from_right = Received.empty(4, K), Received.empty(4, K)
    from_left.has_neighbour[1] = from_right.has_neighbour[1] = True
    from_left.level_rna_pos.write(1, lv_l, 25.0, 100.0)
    from_right.level_rna_pos.write(1, lv_r, 25.0, 100.0)
    from_right.level_rna_neg.write(1, lv_n, 25.0, 100.0)
    rows = prep._cube_rows(from_left, from_right)
    assert set(rows) == {1}
    got = rows[1]
    np.testing.assert_allclose(
        got.profile_pos, intersect([lv_l, lv_r, lower_side(own_pos[1])]), atol=1e-12
    )
    np.testing.assert_allclose(got.profile_neg, intersect([lv_n]), atol=1e-12)
    assert np.array_equal(got.u, u)
    assert (got.total, got.opportunity, got.rho_ref) == (400.0, 100.0, 0.5)


def _empty_piece_ctx(flux: float = 40.0, rate: float = 0.02):
    """A five-slot chain, hand-built: an intron of ``+`` | its ACCEPTOR carrying ``flux`` spliced
    crossings at route rate ``rate`` | an EMPTY exon piece of ``+`` (no unspliced fragment, no RNA
    opportunity — a piece shorter than a fragment) | a plain contiguity boundary (no bits: a face for
    every lane; empty, so it forwards) | a full single-strand ``+`` exon. The strand channel is live (``tau_lam > 0`` at the
    full exon, κ = 0.99)."""

    from rigel.calibration.messages import BlockContext
    from rigel.calibration.splice_graph import FLAG_ACCEPTOR_POS

    n = 5
    fp = np.array([True, True, True, True, True])
    fn = np.zeros(n, bool)
    flags = np.zeros(n, np.uint16)
    flags[1] = FLAG_ACCEPTOR_POS
    n_slot = np.array([300.0, 25.0, 0.0, 0.0, 400.0])
    cnt = np.stack([n_slot * 0.9, n_slot * 0.1], axis=1)
    a_g = np.array([3000.0, 200.0, 40.0, 200.0, 800.0])
    a_r = np.array([3000.0, 200.0, 0.0, 200.0, 800.0])
    sj_hi = np.zeros((n, 2))
    rr_hi = np.zeros((n, 2))
    sj_hi[1, 0], rr_hi[1, 0] = flux, rate
    return BlockContext(
        eff_gdna=a_g,
        eff_rna=a_r,
        sj_count=sj_hi.copy(),
        sj_count_lo=np.zeros((n, 2)),
        sj_count_hi=sj_hi,
        route_rate_lo=np.zeros((n, 2)),
        route_rate_hi=rr_hi,
        unspliced_count=cnt,
        spliced_count=np.zeros((n, 2)),
        left=np.array([-1, 0, 1, 2, 3]),
        right=np.array([1, 2, 3, 4, -1]),
        is_boundary=np.array([False, True, False, True, False]),
        is_exon_region=np.array([False, False, True, False, True]),
        free_pos=fp,
        free_neg=fn,
        exon_pos=np.array([False, False, True, False, True]),
        exon_neg=np.zeros(n, bool),
        boundary_flags=flags,
        has_own_composition=np.array([False, False, False, False, True]),
        belief_fg=np.full(n, 0.5),
        n_grid=41,
        logodds_window=10.0,
        factory_rows=np.zeros((n, 41)),  # a factory with nothing to say: the lanes alone
        strand_live=True,  # the deadband is open: the full exon's split is a witness
    )


def test_an_empty_exon_piece_beside_a_lit_junction_is_a_flux_source():
    """ISSUES: flux-source-skipped-at-an-empty-exon-piece. The junction's flux is a measurement of
    the exon's RNA whether or not the exon piece holds a fragment of its own: at an empty piece the
    flux level is built (lower-sided), priced by `hop_price` on the piece's zero count — counting
    alone, the same rule every hop pays — and emitted from the piece with the flux's own witness
    (the spliced count on the route rate's opportunity), so the next full node prices the hop as
    `flux_level` is priced at a full exon: the junction's rate against its own strand column. A
    silent junction builds nothing and the piece forwards as before."""
    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.calibration.messages.transfer_rows import blur_row, flux_level, hop_price, lower_side

    ctx = _empty_piece_ctx()
    pol = TransferPolicy(strand=(0.99, 0.02, 0.02))
    prepared = _prepared(pol, ctx)
    lane = prepared.lanes["pos"]
    assert lane.empty[2], "the piece must be EMPTY for this gate to say anything"
    assert lane.serves(2, 3) and lane.serves(3, 4)
    own = lane.own_level[2]
    assert own is not None, "no flux level at the empty piece: the source was skipped"
    assert np.all(np.diff(own) >= -1e-12), "the flux level is lower-sided"
    assert lane.flux_at(2, 0) is not None  # its junction, slot 1, is on its left
    # the price is `hop_price` on the piece's own (zero) count: counting alone, both counts
    v = hop_price(40.0, 40.0 / 0.02, 0.0, 0.0)
    np.testing.assert_allclose(own, flux_level(lane.u, 40.0, 0.02, lane.rho_ref, v), atol=1e-12)
    # the piece EMITS its level with the flux's witness, and forwards what it holds beside it
    from rigel.calibration.messages import Levels, Received

    sent = Levels.empty(5, lane.u.shape[0])
    assert lane.emit(2, 3, sent)
    assert (sent.count[3], sent.opportunity[3]) == (40.0, 40.0 / 0.02) and not sent.has_witness[3]
    np.testing.assert_allclose(sent.profile[3], lower_side(own), atol=1e-12)
    # ... through the empty boundary unchanged, and priced at the full exon by the junction's rate
    # against the exon's own column
    table = Received.empty(5, int(ctx.n_grid))
    receive = prepared.propagate(table, backward=False)
    receive(2, 3)
    assert table.level_rna_pos.present[3] and table.level_rna_pos.count[3] == 40.0
    receive(3, 4)
    assert table.level_rna_pos.present[4]
    got = table.level_rna_pos.profile[4]
    v4 = hop_price(40.0, 40.0 / 0.02, float(lane.count[4]), float(lane.a[4]))
    np.testing.assert_allclose(got, blur_row(lower_side(own), lane.u, v4), atol=1e-9)
    # PERTURBATION: a silent junction builds no source, and the empty piece forwards nothing
    quiet = _prepared(
        TransferPolicy(strand=(0.99, 0.02, 0.02)), _empty_piece_ctx(flux=0.0, rate=0.0)
    )
    assert quiet.lanes["pos"].own_level[2] is None
    assert not quiet.lanes["pos"].emit(2, 3, Levels.empty(5, lane.u.shape[0]))


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE RNA LANES EXIST WHENEVER THEIR OWN COORDINATE EXISTS — not when something else does
# ══════════════════════════════════════════════════════════════════════════════════════════════════════
#
# Found by the encompassing-transcript audit (2026-09-13): a single-exon TB− over a two-exon TA+ delivered
# no RNA level anywhere, because three unrelated things each silenced the lanes — the intron factory had
# no rows (no coarse intron in the chain), the gDNA lane had no coordinate (a gDNA-free library), and the
# + lane did not exist (TA+'s exons are both AMBIG, so no single-strand + exon anywhere). A whole
# chromosome never shows any of the three, which is why they shipped.


def _neg_only_ctx():
    """The empty-piece chain mirrored onto the − strand at slots 0..4 with a both-stranded exon
    appended: TB− has a single-strand exon (slot 4) whose own claim is a − level, and the + strand has
    NO single-strand exon anywhere, so the library's + coordinate is zero."""

    from rigel.calibration.messages import BlockContext
    from rigel.calibration.splice_graph import FLAG_ACCEPTOR_POS

    n = 7
    fn = np.array([True, True, True, True, True, True, True])
    flags = np.zeros(n, np.uint16)
    flags[1] = FLAG_ACCEPTOR_POS
    n_slot = np.array([300.0, 25.0, 0.0, 0.0, 400.0, 30.0, 500.0])
    cnt = np.stack([n_slot * 0.1, n_slot * 0.9], axis=1)  # the RNA reads on the − column
    a_g = np.array([3000.0, 200.0, 40.0, 200.0, 800.0, 200.0, 900.0])
    a_r = np.array([3000.0, 200.0, 0.0, 200.0, 800.0, 200.0, 900.0])
    return BlockContext(
        eff_gdna=a_g,
        eff_rna=a_r,
        sj_count=np.zeros((n, 2)),
        sj_count_lo=np.zeros((n, 2)),
        sj_count_hi=np.zeros((n, 2)),
        route_rate_lo=np.zeros((n, 2)),
        route_rate_hi=np.zeros((n, 2)),
        unspliced_count=cnt,
        spliced_count=np.zeros((n, 2)),
        left=np.array([-1, 0, 1, 2, 3, 4, 5]),
        right=np.array([1, 2, 3, 4, 5, 6, -1]),
        is_boundary=np.array([False, True, False, True, False, True, False]),
        is_exon_region=np.array([False, False, True, False, True, False, True]),
        free_pos=np.array([False, False, False, False, False, True, True]),
        free_neg=fn,
        exon_pos=np.array([False, False, False, False, False, False, True]),
        exon_neg=np.array([False, False, True, False, True, False, True]),
        boundary_flags=flags,
        has_own_composition=np.array([False, False, False, False, True, False, False]),
        belief_fg=np.full(n, 0.5),
        n_grid=41,
        logodds_window=10.0,
        factory_rows=np.zeros((n, 41)),
        strand_live=True,
    )


def test_the_lanes_are_built_when_the_intron_factory_has_no_rows():
    """A chain with no coarse intron has no factory rows; the message layer is still the message
    layer. PERTURBATION: the same context with rows builds the same lanes."""
    import dataclasses

    from rigel.calibration.messages.transfer import TransferPolicy

    pol = TransferPolicy(strand=(0.99, 0.02, 0.02))
    with_rows = _prepared(pol, _neg_only_ctx())
    no_rows = _prepared(pol, dataclasses.replace(_neg_only_ctx(), factory_rows=None))
    assert "neg" in with_rows.lanes, "the gate's premise: the − lane exists with rows"
    assert set(no_rows.lanes) == set(with_rows.lanes), (
        f"the lanes vanished with the factory rows: {sorted(no_rows.lanes)} against {sorted(with_rows.lanes)}"
    )


def test_the_rna_lanes_are_built_without_a_gdna_lane():
    """A gDNA-free library has no gDNA level coordinate and hence no gDNA lane; its RNA lanes have
    their own coordinates and must exist. PERTURBATION: with a positive gDNA density the gDNA lane
    joins them."""
    from rigel.calibration.messages.transfer import TransferPolicy, _Library

    pol = TransferPolicy(strand=(0.99, 0.02, 0.02))
    ctx = _neg_only_ctx()
    lib = pol.library(ctx)
    assert lib.rho_rna > 0.0, "the gate's premise: the RNA coordinate exists"
    gdna_free = _Library(0.0, lib.rho_rna, lib.split_live)
    prepared = pol.prepare(ctx, gdna_free)
    assert "gdna" not in prepared.lanes
    assert "neg" in prepared.lanes, "the − lane died with the gDNA lane"
    assert "gdna" in pol.prepare(ctx, lib).lanes


def test_a_strands_level_is_delivered_to_the_cube_when_the_other_strand_has_no_coordinate():
    """TB−'s level must reach the both-stranded exon (slot 6) although the + strand has no
    single-strand exon anywhere, hence no coordinate and no source: the delivery reads whichever RNA
    lane holds something."""
    from rigel.calibration.messages import Received
    from rigel.calibration.messages.transfer import TransferPolicy

    pol = TransferPolicy(strand=(0.99, 0.02, 0.02))
    ctx = _neg_only_ctx()
    prepared = _prepared(pol, ctx)
    assert prepared.lanes["pos"].rho_ref > 0.0, "one coordinate serves both strands"
    K = int(ctx.n_grid)
    from_left, from_right = Received.empty(7, K), Received.empty(7, K)
    from_left.has_neighbour[6] = True
    u = prepared.lanes["neg"].u
    from_left.level_rna_neg.write(6, -0.5 * np.maximum(0.0, (0.0 - u) / 0.3) ** 2, 400.0, 800.0)
    rows = prepared._cube_rows(from_left, from_right)
    assert 6 in rows, "the − level held at the both-stranded exon was not delivered"
    assert rows[6].profile_neg is not None and rows[6].profile_pos is None


def test_a_junctions_flux_is_a_source_when_the_strand_has_no_single_strand_exon():
    """TA+'s exons are both both-stranded, so + has no single-strand exon and, on the old per-strand
    coordinate, no coordinate and no source: its junction's certified flux was silently unused. A level
    is absolute and the coordinate only an origin, so one RNA coordinate serves both strands and the +
    flux level is built at the exon and delivered into its cube. PERTURBATION: the same context with no
    spliced fragments at the junction builds no + level."""
    import dataclasses

    from rigel.calibration.messages import Received
    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.calibration.splice_graph import FLAG_ACCEPTOR_POS

    pol = TransferPolicy(strand=(0.99, 0.02, 0.02))
    base = _neg_only_ctx()
    flags = base.boundary_flags.copy()
    flags[5] = FLAG_ACCEPTOR_POS  # the + junction's acceptor, its exon (slot 6) to the right
    sj_hi, rr_hi = np.zeros((7, 2)), np.zeros((7, 2))
    sj_hi[5, 0], rr_hi[5, 0] = 40.0, 0.02
    lit = dataclasses.replace(
        base, boundary_flags=flags, sj_count=sj_hi.copy(), sj_count_hi=sj_hi, route_rate_hi=rr_hi
    )
    prepared = _prepared(pol, lit)
    lane = prepared.lanes["pos"]
    assert lane.own_level[6] is not None, "the + junction's flux built no level at its exon"
    K = int(lit.n_grid)
    rows = prepared._cube_rows(Received.empty(7, K), Received.empty(7, K))
    assert 6 in rows and rows[6].profile_pos is not None, "the + flux level was not delivered"
    dark = _prepared(pol, base)
    assert dark.lanes["pos"].own_level[6] is None
