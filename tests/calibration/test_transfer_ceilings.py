"""PHASE 2's CEILINGS at single-strand nodes (`transfer._PreparedTransfer._ceilings`, 2026-09-08):
an RNA level of the node's live strand reads as an upper side on its gDNA share and round-trips;
it is read ONLY from a face that sent no composition; and on the live toy the flux is kept per
face and a licensed face keeps the ceiling out, so nothing is counted twice."""

from __future__ import annotations

import numpy as np
import pytest

from _transfer_harness import _rna, _rna_lanes_of, capture_sweep_inputs


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    return capture_sweep_inputs(tmp_path_factory)


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
    from rigel.calibration.messages import Level, Message
    from rigel.calibration.messages.transfer import _LevelLane, _PreparedTransfer, _SolveSite
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

    flux = [None, {0: floor(0.4), 2: floor(0.1)}, None, None, None]
    pos = _LevelLane(
        "pos", u, lam, 0.5, n_u / 2, a_r, empty, [None] * 5, set(), total=n_u, flux=flux
    )
    neg = _LevelLane("neg", u, lam, 0.4, n_u / 2, a_r, empty, [None] * 5, set(), total=n_u)
    gd = _LevelLane("gdna", u, lam, 0.5, n_u, a_r, empty, [None] * 5, set())
    site = _SolveSite(fp & fn, {"pos": fp, "neg": fn}, left, right, 20)
    prep = _PreparedTransfer([None] * 5, {}, K, {"gdna": gd, "pos": pos, "neg": neg}, site)
    comp = -0.5 * ((lam - 1.0) / 0.5) ** 2
    held_l = Level(floor(0.8), 25.0, 100.0)
    held_r = Level(floor(-0.3), 25.0, 100.0)
    from_left = [None, Message(composition=comp, level_rna_pos=held_l), None, None, None]
    from_right = [None, Message(level_rna_pos=held_r), None, None, None]
    rows = np.zeros((5, K))
    assert prep._ceilings(from_left, from_right, rows)
    want = rna_row_of_level(intersect([held_r.profile, flux[1][2]]), u, lam, 400.0, 100.0, 0.5)
    np.testing.assert_allclose(rows[1], want, atol=1e-12)
    assert np.all(np.diff(rows[1]) <= 1e-9)
    assert not rows[[0, 2, 3, 4]].any()
    # the perturbation: the left composition removed → its level and flux join
    from_left2 = [None, Message(level_rna_pos=held_l), None, None, None]
    rows2 = np.zeros((5, K))
    assert prep._ceilings(from_left2, from_right, rows2)
    want2 = rna_row_of_level(
        intersect([held_l.profile, flux[1][0], held_r.profile, flux[1][2]]),
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
            fx = lane.flux[x]
            if fx is None:
                continue
            for b, prof in fx.items():
                assert b in (int(ctx.left[x]), int(ctx.right[x])) and sc[b].sum() > 0
                assert np.all(np.diff(prof) >= -1e-12)  # a lower bound
                per_face += 1
    assert per_face > 0
    from rigel.calibration.messages import SILENCE

    order = list(range(ctx.n_slots))
    held = []
    for nbr, seq, backward in (
        (np.asarray(ctx.left, np.int64), order, False),
        (np.asarray(ctx.right, np.int64), order[::-1], True),
    ):
        receive = prepared.propagate(backward=backward)
        got = [None] * len(order)
        for i in seq:
            s = int(nbr[i])
            if s >= 0:
                got[i] = SILENCE if receive is None else receive(s, i)
        held.append(got)
    msg = prepared.solve(*held)
    rows = np.asarray(msg.lam_rows)
    without = np.zeros_like(rows)
    n = len(prepared.own)
    from rigel.calibration.messages.transfer import _fuse
    from rigel.calibration.messages.transfer_rows import intersect

    for i in range(n):
        parts, bounds = [], []
        for m in (held[0][i], held[1][i]):
            if m is None:
                continue
            if m.composition is not None:
                parts.append(m.composition)
            if m.level_gdna is not None and not prepared.lanes["gdna"].empty[i]:
                bounds.append(prepared.lanes["gdna"].row(m.level_gdna.profile, i))
        if bounds:
            parts.append(intersect(bounds))
        if parts:
            without[i] = _fuse(parts)
    site = prepared.site
    fp, fn = site.free["pos"], site.free["neg"]
    all_comp = np.array(
        [
            all(m is None or m.composition is not None for m in (held[0][i], held[1][i]))
            for i in range(n)
        ]
    )
    single = (fp ^ fn) & ~site.ambig
    quiet = single & all_comp
    assert quiet.sum() > 0
    np.testing.assert_allclose(rows[quiet], without[quiet], atol=1e-12)
