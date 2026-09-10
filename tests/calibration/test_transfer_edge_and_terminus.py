"""THE EDGE'S LEVEL and THE TERMINUS (`transfer._edge_level`, `transfer._terminus_rules`): the
one-sided edge level with a vacuous zero count; the terminus orientation read from the flag alone;
the sj+terminus boundary placing the junction's flux on the junction's exon side; the level-kept
map and the crossing bound; and on the live toy the composition with the OUTSIDE exon at exactly
the terminus pairs (nowhere once the flags are cleared) and THE LEVEL RULE into every inside,
made from the sender's measurement alone."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import polygamma

from _transfer_harness import (
    _ctx_of,
    _full_policy,
    _strand_of,
    _strand_row_of,
    _with_populated_inside,
    capture_sweep_inputs,
)


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    return capture_sweep_inputs(tmp_path_factory)


def test_the_edge_level_is_one_sided_and_a_zero_count_is_vacuous():
    """THE EDGE'S LEVEL (the owner's design 2026-09-04, the form the ladder kept): below the edge's level the
    count's exact Poisson — the exon has at least the edge's gDNA density; NOTHING above it (no local
    witness prices capture's enrichment of the interior); a ZERO count is vacuous — darkness under
    capture is not absence."""
    from rigel.calibration.messages.transfer_rows import edge_level_row

    lam = np.linspace(-8, 8, 801)
    sig = 1 / (1 + np.exp(-lam))
    n_e, a_b, a_e = 400.0, 200.0, 800.0
    row = edge_level_row(
        lam, 25.0, n_e, a_b, a_e
    )  # the implied edge count is 100 f: the level at f = 0.25
    c = sig * n_e * a_b / a_e
    below = c < 25.0
    np.testing.assert_allclose(row[below], 25 * np.log(c[below] / 25) - (c[below] - 25), atol=1e-9)
    assert np.all(row[~below] == 0.0), "nothing above the edge's level"
    assert float(np.interp(np.log(0.1 / 0.9), lam, row)) < -3.0, (
        "below it the count's own Poisson charges"
    )
    assert not edge_level_row(lam, 0.0, n_e, a_b, a_e).any(), "a zero count is vacuous"


def _expected_terminus_rows(si, ctx, strand, lam, exon_rows):
    """Item 5 recomputed INDEPENDENTLY of the policy. At every exon|exon boundary carrying exactly one
    terminus direction and no splice junction, the OUTSIDE flank is read off the flag alone (TSS+ and
    TES− bodies extend right, outside = LEFT; TES+ and TSS− extend left, outside = RIGHT). The boundary
    receives (i) the splice-in and edge rows of its outside exon (``exon_rows``, the independent recompute) and
    (ii) the exon's own frozen-variance strand row when its channel is live, each through the
    splice-out map with S = the boundary's SPLICED crossing; the outside exon receives the boundary's
    own strand row through the face map with that spliced density, widened by the counting variance.
    Returns ``(rows_at_boundaries, rows_at_outside_exons)`` keyed by slot."""
    from rigel.calibration.messages.transfer_rows import (
        face_map_lambda,
        splice_out_row,
        transport_row,
    )
    from rigel.calibration.splice_graph import (
        FLAG_ACCEPTOR_NEG,
        FLAG_ACCEPTOR_POS,
        FLAG_DONOR_NEG,
        FLAG_DONOR_POS,
        FLAG_TES_NEG,
        FLAG_TES_POS,
        FLAG_TSS_NEG,
        FLAG_TSS_POS,
    )

    kappa, od_g, od_r = strand
    sj = FLAG_DONOR_POS | FLAG_DONOR_NEG | FLAG_ACCEPTOR_POS | FLAG_ACCEPTOR_NEG
    body_right, body_left = FLAG_TSS_POS | FLAG_TES_NEG, FLAG_TES_POS | FLAG_TSS_NEG
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    n_u = np.asarray(ctx.n_slot, np.float64)
    n_s = np.asarray(ctx.spliced_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna, np.float64)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    tau = np.asarray(ctx.own.tau_lam, np.float64)
    belief = np.asarray(ctx.belief_fg, np.float64)
    fg = 1.0 / (1.0 + np.exp(-lam))

    def strand_row(x):
        n = cnt[x].sum()
        ks = kappa if fp[x] else 1.0 - kappa
        f_ref = float(np.clip(belief[x], 1e-9, 1 - 1e-9))
        p = 0.5 * fg + ks * (1 - fg)
        p_ref = 0.5 * f_ref + ks * (1 - f_ref)
        var = max(
            n * p_ref * (1 - p_ref)
            + (n * f_ref) ** 2 * 0.25 * od_g
            + (n * (1 - f_ref)) ** 2 * ks * (1 - ks) * od_r,
            1e-9,
        )
        row = -0.5 * (cnt[x, 0] - n * p) ** 2 / var
        return row - row.max()

    at_b, at_o = {}, {}
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0 or not (is_exon[lo] and is_exon[hi]):
            continue
        f = int(flags[b])
        to_right, to_left = bool(f & body_right), bool(f & body_left)
        if (f & sj) or to_right == to_left:
            continue
        o = lo if to_right else hi
        if fp[b] != fp[o] or fn[b] != fn[o] or fp[b] == fn[b]:
            continue
        if not (n_u[b] > 0 and A_g[b] > 0 and A_g[o] > 0):
            continue
        acc = np.zeros(lam.shape[0])
        if int(o) in exon_rows:
            acc += splice_out_row(exon_rows[int(o)], lam, n_u[b], n_s[b], A_g[b], A_g[o])
        if tau[o] > 0.0:
            acc += splice_out_row(strand_row(o), lam, n_u[b], n_s[b], A_g[b], A_g[o])
        if np.ptp(acc) > 1e-9:
            at_b[int(b)] = acc
        if tau[b] > 0.0:
            le = face_map_lambda(lam, n_u[b], A_g[b], A_g[b], A_g[o], A_g[o], n_s[b] / A_g[b])
            at_o.setdefault(int(o), np.zeros(lam.shape[0]))
            at_o[int(o)] += transport_row(strand_row(b), lam, le, n_u[b], n_s[b])
    return at_b, at_o


def _item5_slots(ctx):
    """Item 5's delivery sites — every exon|exon boundary carrying one terminus direction and no sj
    whose outside flank shares its single strand, and that outside exon — derived from the flags
    alone, so the earlier gates can hand these slots to the terminus's own gate."""
    from rigel.calibration.messages.transfer_rows import boundary_shares_strand, outside_flank

    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    out = set()
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0 or not (is_exon[lo] and is_exon[hi]):
            continue
        o, _i = outside_flank(flags[b], lo, hi)
        if o is not None and boundary_shares_strand(fp[b], fn[b], fp[o], fn[o]):
            out |= {int(b), int(o)}
    return out


def _terminus_flags_cleared(ctx):
    """The context with every exon|exon boundary's terminus bits cleared — the terminus rules see none,
    every other message is unchanged (they read flags at intron|exon faces only)."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer_rows import TERMINUS

    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16).copy()
    for b in np.flatnonzero(is_bnd):
        if left[b] >= 0 and right[b] >= 0 and is_exon[left[b]] and is_exon[right[b]]:
            flags[b] &= ~TERMINUS
    return _dc.replace(ctx, boundary_flags=flags)


def test_the_terminus_orientation_reads_the_flag_alone():
    """The orientation table, gated DIRECTLY: a + start or a − end extends genomic-right, so the
    OUTSIDE flank is the left one; a + end or a − start extends left, so it is the right one; termini
    pointing both ways or no terminus at all give no side. ⭐ A splice junction sharing the boundary does
    NOT change the side (the sj+terminus case, 2026-09-08): the four ladder families resolve as their
    terminus does, and a junction with no terminus still gives no side (the perturbation)."""
    from rigel.calibration.messages.transfer_rows import junction_exon_side, outside_flank
    from rigel.calibration.splice_graph import (
        FLAG_ACCEPTOR_NEG,
        FLAG_ACCEPTOR_POS,
        FLAG_DONOR_NEG,
        FLAG_DONOR_POS,
        FLAG_TES_NEG,
        FLAG_TES_POS,
        FLAG_TSS_NEG,
        FLAG_TSS_POS,
    )

    assert outside_flank(FLAG_TSS_POS, 7, 9) == (7, 9)
    assert outside_flank(FLAG_TES_NEG, 7, 9) == (7, 9)
    assert outside_flank(FLAG_TES_POS, 7, 9) == (9, 7)
    assert outside_flank(FLAG_TSS_NEG, 7, 9) == (9, 7)
    assert outside_flank(FLAG_TSS_POS | FLAG_TES_POS, 7, 9) == (None, None), "both ways: no side"
    assert outside_flank(0, 7, 9) == (None, None), "no terminus: no side"
    # the sj+terminus families: the terminus decides, the junction says where its flux belongs
    assert outside_flank(FLAG_TSS_POS | FLAG_ACCEPTOR_POS, 7, 9) == (
        7,
        9,
    )  # a start at an exon's low edge
    assert outside_flank(FLAG_TES_POS | FLAG_DONOR_POS, 7, 9) == (
        9,
        7,
    )  # an end at an exon's high edge
    assert outside_flank(FLAG_TSS_NEG | FLAG_DONOR_NEG, 7, 9) == (
        9,
        7,
    )  # a − start at an exon's high edge
    assert outside_flank(FLAG_TES_NEG | FLAG_ACCEPTOR_NEG, 7, 9) == (
        7,
        9,
    )  # a − end at an exon's low edge
    assert (
        junction_exon_side(FLAG_TSS_POS | FLAG_ACCEPTOR_POS, 7, 9) == 9
    )  # ACC: the intron is left
    assert junction_exon_side(FLAG_TES_POS | FLAG_DONOR_POS, 7, 9) == 7  # DON: the intron is right
    assert junction_exon_side(FLAG_TSS_POS, 7, 9) is None
    assert junction_exon_side(FLAG_DONOR_POS | FLAG_ACCEPTOR_POS, 7, 9) is None, (
        "junctions both ways"
    )
    # the perturbations: a junction alone gives no side; termini both ways with a junction give none
    assert outside_flank(FLAG_DONOR_POS, 7, 9) == (None, None)
    assert outside_flank(FLAG_ACCEPTOR_NEG | FLAG_DONOR_NEG, 7, 9) == (None, None)
    assert outside_flank(FLAG_TSS_POS | FLAG_TES_POS | FLAG_ACCEPTOR_POS, 7, 9) == (None, None)


def test_the_sj_terminus_boundary_places_the_flux_where_the_junctions_exon_is(sweep_inputs):
    """On the toy with one boundary made an sj+terminus boundary (a TSS flag added to a licensed
    acceptor): the terminus rule now serves the inside exon and its totals' disagreement carries the
    junction's flux (a lower price than without it, since the RNA joining at the junction is measured),
    while the strand-mode prediction keeps the crossing alone; the junction rules (the splice-in map, the
    alternative splice site) leave
    that face. PERTURBATION: with the flux zeroed the price rises back to the plain form's."""
    import dataclasses

    from rigel.calibration.messages.transfer_rows import SJ_FLAGS, TERMINUS, junction_exon_side
    from rigel.calibration.splice_graph import FLAG_ACCEPTOR_POS, FLAG_TSS_POS

    pol, _p, _g, _w = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp = np.asarray(ctx.free_pos, bool)
    flags = np.asarray(ctx.boundary_flags, np.uint16).copy()
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    n_u = np.asarray(ctx.n_slot, np.float64)
    flux = np.asarray(ctx.sj_count, np.float64).sum(axis=1)
    # a + acceptor whose right flank is an exon with a total and a junction flux, no terminus yet
    cands = [
        b
        for b in np.flatnonzero(is_bnd)
        if int(flags[b]) == int(FLAG_ACCEPTOR_POS)
        and right[b] >= 0
        and is_exon[right[b]]
        and fp[b]
        and n_u[b] > 0
        and n_u[right[b]] > 0
        and flux[b] > 0
    ]
    assert cands, "no plain + acceptor with flux on the toy"
    b = int(cands[0])
    i = int(right[b])
    flags[b] = np.uint16(int(flags[b]) | int(FLAG_TSS_POS))
    ctx2 = dataclasses.replace(ctx, boundary_flags=flags)
    assert (int(flags[b]) & TERMINUS) and (int(flags[b]) & SJ_FLAGS)
    assert junction_exon_side(flags[b], left[b], right[b]) == i  # the junction's exon is the inside
    prep = pol.prepare(ctx2)
    rule = prep.rule.get((b, i))
    assert rule is not None and getattr(rule, "__name__", "") == "level_rule"
    v_with = rule.__defaults__[1]
    # the plain form: the same boundary with its flux zeroed
    sjc = np.asarray(ctx.sj_count, np.float64).copy()
    sjc[b] = 0.0
    prep0 = pol.prepare(dataclasses.replace(ctx2, sj_count=sjc))
    rule0 = prep0.rule.get((b, i))
    assert rule0 is not None and getattr(rule0, "__name__", "") == "level_rule"
    v_plain = rule0.__defaults__[1]
    assert v_with < v_plain, (v_with, v_plain)
    # the junction rules leave the face: the splice-in map into the inside and the splice-out map out of it are gone
    prep_before = pol.prepare(ctx)
    assert (b, i) in prep_before.rule and getattr(
        prep_before.rule[(b, i)], "__name__", ""
    ) != "level_rule"


def _expected_level_rows(si, ctx, strand, lam):
    """THE LEVEL RULE recomputed INDEPENDENTLY of the policy: at every terminus boundary with an inside
    EXON (the outside an exon or an intron), the boundary's own strand row through the level-kept map
    (the boundary's share times its crossing density, times the inside's opportunity, over the
    inside's total), blurred by both totals' counting plus the pair's discrepancies — the totals'
    disagreement beyond counting and, where both strand channels are live, the two strand modes'
    disagreement beyond counting. With no own row the crossing total's one-sided upper bound. Keyed by
    the inside slot; ``(rows, served pairs)``."""
    from rigel.calibration.messages.transfer_rows import (
        boundary_shares_strand,
        level_bound_row,
        level_map_lambda,
        level_row,
        outside_flank,
    )

    kappa, od_g, od_r = strand
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    is_intron = ~is_bnd & ~is_exon & (fp | fn)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    n_u = np.asarray(ctx.n_slot, np.float64)
    n_s = np.asarray(ctx.spliced_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna, np.float64)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    tau = np.asarray(ctx.own.tau_lam, np.float64)
    out, served = {}, []
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0:
            continue
        o, i = outside_flank(flags[b], lo, hi)
        if i is None or not is_exon[i] or not (is_exon[o] or is_intron[o]):
            continue
        if not boundary_shares_strand(fp[b], fn[b], fp[i], fn[i]):
            continue
        if not (n_u[b] > 0 and n_u[i] > 0 and A_g[b] > 0 and A_g[i] > 0):
            continue
        served.append((int(b), int(i), "exon|exon" if is_exon[o] else "exon|intron"))
        d_b, T_b = n_u[b] / A_g[b], n_u[b] + n_s[b]
        r = (n_u[i] / A_g[i]) / (T_b / A_g[b])
        v = max(0.0, np.log(r) ** 2 - (1.0 / n_u[i] + 1.0 / T_b))
        if tau[b] > 0.0 and tau[i] > 0.0:
            ok, modes = True, []
            for y in (b, i):
                n = cnt[y].sum()
                p = cnt[y, 0] / n
                ks = kappa if fp[y] else 1.0 - kappa
                f = (p - ks) / (0.5 - ks)
                if not 0.0 < f < 1.0:
                    ok = False
                    break
                modes.append((f, p * (1 - p) / n / (p - ks) ** 2 / (1 - f) ** 2))
            if ok:
                (f_b, v_b), (f_i, v_i) = modes
                f_pred = min(f_b / r, 1 - 1e-9)
                dd = np.log(f_i / (1 - f_i)) - np.log(f_pred / (1 - f_pred))
                v += max(0.0, dd * dd - (v_b + v_i + 1.0 / n_u[i] + 1.0 / T_b))
        v += float(polygamma(1, n_u[b] + 0.5) + polygamma(1, n_u[i] + 0.5))
        m = level_map_lambda(lam, d_b, A_g[i], n_u[i])
        if tau[b] > 0.0 and fp[b] != fn[b]:
            row = level_row(_strand_row_of(ctx, strand, lam, b), lam, m, v)
        else:
            row = level_bound_row(lam, d_b, A_g[i], n_u[i], v)
        if np.ptp(row) > 1e-9:
            out.setdefault(int(i), np.zeros(lam.shape[0]))
            out[int(i)] += row
    return out, served


def test_the_level_map_keeps_the_level_and_the_bound_is_vacuous_below_the_total():
    """THE LEVEL RULE's arithmetic, gated directly: through the level-kept map a boundary whose share is
    f_b and whose crossing density is d lands the inside at f_b * d * E_i / T_i (the gDNA level kept,
    the inside's own total supplying the rest); the map is monotone; a profile read through it keeps
    its peak there and its width grows with the dampening; the total's upper bound charges nothing
    below the crossing's density and charges above it."""
    from rigel.calibration.messages.transfer_rows import (
        level_bound_row,
        level_map_lambda,
        level_row,
    )

    lam = np.linspace(-8, 8, 801)
    sig = 1 / (1 + np.exp(-lam))
    d_b, E_i, T_i = (
        0.2,
        500.0,
        200.0,
    )  # the boundary crosses 0.2 fragments/base; the inside holds 200
    m = level_map_lambda(lam, d_b, E_i, T_i)
    assert np.all(np.diff(m) >= 0.0)
    f_b = 0.4
    row_b = -0.5 * ((lam - np.log(f_b / (1 - f_b))) / 0.05) ** 2
    want = f_b * d_b * E_i / T_i  # 0.2: the level kept
    tight = level_row(row_b, lam, m, 0.001)
    assert float(sig[np.argmax(tight)]) == pytest.approx(want, abs=0.01)

    def _var(r):
        w = np.exp(r - r.max())
        w /= w.sum()
        mu = w @ lam
        return float(w @ (lam * lam) - mu * mu)

    wide = level_row(row_b, lam, m, 0.5)
    assert _var(wide) > 3 * _var(tight), "the dampening must widen the delivered profile"
    assert not level_row(np.zeros_like(lam), lam, m, 0.1).any(), "a flat profile is vacuous"
    ub = level_bound_row(lam, d_b, E_i, T_i, 0.01)
    below = sig * T_i / E_i < d_b  # inside gDNA density below the crossing's total density
    assert np.all(ub[below] == 0.0) and np.all(ub[~below] <= 0.0) and np.any(ub[~below] < 0.0)


def test_the_terminus_rules_land_at_the_outside_pair_and_nowhere_when_the_flags_clear(sweep_inputs):
    """ITEM 5 as claims and rules: at every exon|exon terminus boundary the rule outside exon → boundary
    applied to the exon's claim is the independently recomputed splice-out row with the SPLICED
    crossing, the rule boundary → outside exon applied to the boundary's claim the recomputed face-map
    row; the rules exist exactly at the outside pairs and vanish when the terminus bits are cleared."""
    from rigel.calibration.messages.transfer_rows import (
        face_map_lambda,
        splice_out_row,
        transport_row,
    )
    from rigel.calibration.simplex_logodds import _logodds_grid

    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    prepared = pol.prepare(ctx)
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    n_u = np.asarray(ctx.n_slot, np.float64)
    n_s = np.asarray(ctx.spliced_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna, np.float64)
    tau = np.asarray(ctx.own.tau_lam, np.float64)
    sites = _item5_slots(ctx)
    served = 0
    for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
        if not (is_exon[left[b]] and is_exon[right[b]]):
            continue
        outs = [
            (s, i) for (s, i) in prepared.rule if (s == b and is_exon[i]) or (i == b and is_exon[s])
        ]
        if int(b) not in sites:
            assert not outs, f"boundary {b} carries exon|exon rules without a served terminus"
            continue
        served += 1
        (o,) = (
            {s if s != b else i for (s, i) in outs} & {int(left[b]), int(right[b])}
            if outs
            else (None,)
        )
        assert o is not None
        want = splice_out_row(
            _strand_row_of(ctx, strand, lam, o), lam, n_u[b], n_s[b], A_g[b], A_g[o]
        )
        if tau[o] > 0.0 and prepared.own[o] is not None:
            got = prepared.rule[(int(o), int(b))](prepared.own[o], None)
            np.testing.assert_allclose(got - got.max(), want - want.max(), rtol=0, atol=1e-10)
        if tau[b] > 0.0:
            le = face_map_lambda(lam, n_u[b], A_g[b], A_g[b], A_g[o], A_g[o], n_s[b] / A_g[b])
            want_o = transport_row(_strand_row_of(ctx, strand, lam, b), lam, le, n_u[b], n_s[b])
            got_o = prepared.rule[(int(b), int(o))](prepared.own[b], None)
            np.testing.assert_allclose(
                got_o - got_o.max(), want_o - want_o.max(), rtol=0, atol=1e-10
            )
    assert served >= 2, (
        "the toy must carry two served terminus boundaries or this gate proves nothing"
    )
    cleared = pol.prepare(_terminus_flags_cleared(ctx))
    for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
        if is_exon[left[b]] and is_exon[right[b]]:
            assert not any(s == b or i == b for (s, i) in cleared.rule), (
                f"rules survive at {b} with no terminus"
            )


def test_the_level_rule_serves_every_terminus_inside_from_the_measurement_alone(sweep_inputs):
    """THE LEVEL RULE (the owner's design 2026-09-04) as a rule: at every terminus boundary with an inside exon —
    exon|exon and exon|intron alike — the rule applied to the boundary's OWN claim equals the
    independent recompute; what the boundary HOLDS changes nothing (an imputation never becomes a
    level); with no own claim the crossing total's upper bound is what crosses; the rule vanishes when
    the terminus bits are cleared; and no library-wide parameter exists — changing another pair's
    counts leaves this pair's message unchanged."""
    import dataclasses as _dc

    from rigel.calibration.simplex_logodds import _logodds_grid

    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _with_populated_inside(_ctx_of(sweep_inputs))
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    prepared = pol.prepare(ctx)
    expected, served = _expected_level_rows(sweep_inputs, ctx, strand, lam)
    assert len(served) >= 2, (
        "the toy must carry two served terminus pairs or this gate proves nothing"
    )
    spike = -0.5 * ((lam - 2.0) / 0.1) ** 2  # a held imputation that must not cross a level face
    for b, i, _kind in served:
        fn = prepared.rule.get((b, i))
        assert fn is not None, f"no level rule at terminus pair ({b}, {i})"
        got = fn(prepared.own[b], None)
        np.testing.assert_allclose(
            got - got.max(), expected[i] - expected[i].max(), rtol=0, atol=1e-10
        )
        np.testing.assert_array_equal(
            fn(prepared.own[b], spike), got, err_msg="what is held crossed a level face"
        )
        bound = fn(None, spike)
        assert bound is not None and np.ptp(bound) > 0.0 and np.all(bound <= 0.0), (
            "no upper bound without a claim"
        )
    cleared = pol.prepare(_terminus_flags_cleared(ctx))
    for b, i, _kind in served:
        assert (b, i) not in cleared.rule, f"a level rule survives at ({b}, {i}) with no terminus"
    b0, i0, _k = served[0]
    b1, i1, _k = served[-1]
    cnt = np.asarray(ctx.unspliced_count, np.float64).copy()
    cnt[i1] = cnt[i1] * 3.0 + 7.0
    other = pol.prepare(_dc.replace(ctx, unspliced_count=cnt, n_slot=cnt.sum(axis=1)))
    np.testing.assert_array_equal(
        other.rule[(b0, i0)](other.own[b0], None), prepared.rule[(b0, i0)](prepared.own[b0], None)
    )
