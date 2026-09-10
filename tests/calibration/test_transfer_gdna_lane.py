"""THE LEVEL LANE on the gDNA population (`transfer._gdna_lane`, `transfer._LevelLane`): a
received level is a lower bound the hop widens; the level coordinates round-trip through a node's
total; bounds INTERSECT and do not multiply (the ratchet); a full node emits the intersection of its
own lower side and what it holds; the hop price is both countings and the discrepancy beyond them;
and an empty node forwards the level it holds unchanged."""

from __future__ import annotations

import numpy as np
import pytest

from _transfer_harness import _ctx_of, _drive_the_backbone, _full_policy, capture_sweep_inputs


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    return capture_sweep_inputs(tmp_path_factory)


def test_a_received_gdna_level_is_a_lower_bound_and_the_hop_widens_it():
    """`_LevelLane.receive` on the gDNA lane (no two-sided face): what a full recipient holds is
    non-decreasing in u (a level that crosses a face says "at least this much gDNA" and nothing more),
    a two-sided input loses only its upper side, and a hop across a density cliff — a larger price —
    widens it."""
    from rigel.calibration.messages import Level
    from rigel.calibration.messages.transfer import _LevelLane

    u = np.linspace(-10, 10, 60)
    two_sided = -0.5 * ((u - 1.0) / 0.4) ** 2
    two_sided -= two_sided.max()
    # one density on both sides: the price is the two totals' counting alone
    flat = _LevelLane(
        "gdna",
        u,
        u,
        0.05,
        np.array([1.0e6, 1.0e6]),
        np.array([1.0e4, 1.0e4]),
        np.zeros(2, bool),
        [None, None],
        {(0, 1)},
    )
    lower = flat.receive(Level(two_sided, 1.0e6, 1.0e4), 0, 1).profile
    assert np.all(np.diff(lower) >= -1e-12), "not a lower bound"
    below = u < 0.0
    np.testing.assert_allclose(lower[below], two_sided[below], atol=1e-3)
    assert np.all(lower[u > 2.0] > -1e-9), "the upper side was kept"
    # a cliff between the two nodes: the discrepancy beyond counting widens the bound
    cliff = _LevelLane(
        "gdna",
        u,
        u,
        0.05,
        np.array([1.0e6, 1.0e6]),
        np.array([1.0e4, 4.0e4]),
        np.zeros(2, bool),
        [None, None],
        {(0, 1)},
    )
    wide = cliff.receive(Level(two_sided, 1.0e6, 1.0e4), 0, 1).profile
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
    """THE RATCHET GATE: on a hand-built lane (three full nodes in a row, which the toy does not have),
    the level a full node emits is the pointwise tighter of its own lower side and the priced level it
    holds — never their sum, which sharpened a chain of nine one-fragment boundaries into a hard bound
    on the ladder. PERTURBATION: an `emit` that multiplies fails here."""
    from rigel.calibration.messages import Level
    from rigel.calibration.messages.transfer import _LevelLane
    from rigel.calibration.messages.transfer_rows import intersect, lower_side

    u = np.linspace(-10, 10, 60)
    own = [None, -0.5 * ((u + 1.0) / 1.5) ** 2, None]
    lane = _LevelLane(
        "gdna",
        u,
        u,
        0.05,
        np.array([100.0, 120.0, 90.0]),
        np.array([50.0, 60.0, 45.0]),
        np.zeros(3, bool),
        own,
        {(0, 1), (1, 2)},
    )
    held = Level(lower_side(-0.5 * ((u - 0.5) / 0.8) ** 2), 100.0, 50.0)
    sent = lane.emit(1, 2, held)
    want = intersect([lower_side(own[1]), held.profile])
    np.testing.assert_allclose(sent.profile, want, atol=1e-12)
    product = lower_side(own[1]) + held.profile
    assert np.max(np.abs(sent.profile - (product - product.max()))) > 0.5, "the sum: the ratchet"
    assert (sent.n, sent.a) == (120.0, 60.0)


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


def test_an_empty_node_forwards_the_level_it_holds_unchanged(sweep_inputs):
    """THE EMPTY NODE IS TRANSPARENT: the toy's inside pieces have no total, so what the boundary beyond
    such a piece holds from it must be exactly what the piece received — same profile, same (n, a) of
    the last full node — priced only at the full recipient. PERTURBATION: a policy whose empty nodes
    price the hop fails this gate."""
    from rigel.calibration.messages.transfer_rows import blur_row, hop_price, lower_side

    ctx = _ctx_of(sweep_inputs)
    prepared = _full_policy(sweep_inputs)[0].prepare(ctx)
    _drive_the_backbone(prepared, ctx)
    lane = prepared.lanes["gdna"]
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    checked = 0
    for backward, nbr in ((False, left), (True, right)):
        held = prepared.held[backward]
        for e in np.flatnonzero(lane.empty):
            s = int(nbr[e])
            if s < 0 or held[e] is None or held[e].level_gdna is None:
                continue
            x = int(right[e] if not backward else left[e])
            if x < 0 or held[x] is None or held[x].level_gdna is None or lane.empty[x]:
                continue
            got, sent = held[x].level_gdna, held[e].level_gdna
            v = hop_price(sent.n, sent.a, lane.count[x], lane.a[x])
            want = lower_side(sent.profile)
            want = blur_row(want, lane.u, v) if v > 0.0 else want
            np.testing.assert_allclose(got.profile, want, atol=1e-9)
            assert (got.n, got.a) == (float(lane.count[x]), float(lane.a[x]))
            checked += 1
    assert checked > 0, "no level crossed an empty node on the toy: the gate proved nothing"
