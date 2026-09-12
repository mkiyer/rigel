"""Gates for the transfer policy's per-face message builders, every row recomputed on the live
toy by a second implementation so a policy bug cannot hide.

The intron|exon face: the licence refuses unmeasured population changes; the splice-in face map is
monotone and capped by a measured flux, with its width exactly the face's counting variance; the
splice-out row is the count form marginalised over the ratio's counting; the boundary-to-intron
licence needs one shared strand; and every row lands at exactly the licensed slots. The edge and
the terminus: the one-sided edge level with a vacuous zero count; the orientation read from the
flag alone; the sj+terminus boundary placing the junction's flux on the junction's exon side; the
level-kept map and the crossing bound; and the level rule into every inside. The alternative
splice site: the flanks from the flag kind, and both flanks' rules, priced per pair.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import polygamma

from _transfer_harness import (
    _ctx_of,
    _expected_pairs,
    _full_policy,
    _intron_mask,
    _prepared,
    _strand_of,
    _strand_row_of,
    _with_alt_splice_sites,
    _with_populated_inside,
)


# ── the intron|exon face: the licence, the splice-in map and the splice-out row ──────────────────


def _expected_exon_rows(si, ctx, src, lam, strand_of_si=None):
    """The splice-in exon rows recomputed INDEPENDENTLY of the policy: for every licensed
    intron|exon face of every exon, the intron row pushed through the face map and widened by
    the face's counting variance — a second implementation, so a policy bug cannot hide."""
    from scipy.special import polygamma

    from rigel.calibration.splice_graph import (
        FLAG_TES_NEG,
        FLAG_TES_POS,
        FLAG_TSS_NEG,
        FLAG_TSS_POS,
    )

    term = FLAG_TSS_POS | FLAG_TSS_NEG | FLAG_TES_POS | FLAG_TES_NEG
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    is_intron = ~is_bnd & ~is_exon & (fp | fn)
    left = np.asarray(ctx.left, np.int64)
    right = np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    n_u = np.asarray(ctx.n_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna, np.float64)
    A_r = np.asarray(ctx.eff_rna, np.float64)
    rr_lo = np.asarray(ctx.route_rate_lo, np.float64).sum(axis=1)
    rr_hi = np.asarray(ctx.route_rate_hi, np.float64).sum(axis=1)
    sc_lo = np.asarray(ctx.sj_count_lo, np.float64).sum(axis=1)
    sc_hi = np.asarray(ctx.sj_count_hi, np.float64).sum(axis=1)
    is_intergenic = ~is_bnd & ~is_exon & ~fp & ~fn
    sig = 1.0 / (1.0 + np.exp(-lam))
    dlam = float(lam[1] - lam[0])
    out = {}
    for e in np.flatnonzero(is_exon):
        add = np.zeros(lam.shape[0])
        live = False
        for b in (left[e], right[e]):
            # the intergenic|exon EDGE's level, lower bound only
            if b < 0 or not is_bnd[b]:
                continue
            o = left[b] if right[b] == e else right[b]
            if o < 0 or not is_intergenic[o]:
                continue
            if not (A_g[b] > 0 and A_g[e] > 0 and n_u[e] > 0):
                continue
            from rigel.calibration.messages.transfer_rows import edge_level_row

            r = edge_level_row(lam, n_u[b], n_u[e], A_g[b], A_g[e])
            if np.ptp(r) <= 1e-9:
                continue
            add += r - r.max()
            live = True
        for b, hi in ((left[e], True), (right[e], False)):
            if b < 0 or not is_bnd[b]:
                continue
            i = left[b] if right[b] == e else right[b]
            if i < 0 or not is_intron[i] or (flags[b] & term):
                continue
            if fp[e] != fp[i] or fn[e] != fn[i]:
                continue
            row_i = src[i]
            if np.ptp(row_i) <= 1e-9 or not (
                n_u[b] > 0 and A_g[b] > 0 and A_r[b] > 0 and A_g[e] > 0 and A_r[e] > 0
            ):
                continue
            s = float((rr_hi if hi else rr_lo)[b])
            n_s = float((sc_hi if hi else sc_lo)[b])
            le = np.log(n_u[b] * sig / A_g[b] * A_g[e]) - np.log(
                (n_u[b] * (1.0 - sig) / A_r[b] + s) * A_r[e]
            )
            r = np.interp(
                np.interp(lam, le, lam, left=lam[0], right=lam[-1]), lam, row_i - row_i.max()
            )
            v = float(polygamma(1, n_u[b] + 0.5) + polygamma(1, n_s + 0.5))
            half = max(int(np.ceil(4.0 * np.sqrt(v) / dlam)), 1)
            x = np.arange(-half, half + 1) * dlam
            kern = np.exp(-0.5 * x * x / v)
            kern /= kern.sum()
            pr = np.convolve(np.pad(np.exp(r - r.max()), half, mode="edge"), kern, mode="valid")
            r = np.log(np.maximum(pr, 1e-300))
            add += r - r.max()
            live = True
        if live:
            out[int(e)] = add - add.max()
    return out


def test_the_face_licence_refuses_unmeasured_population_changes():
    """The licence predicate, gated directly: the integration gate cannot falsify the terminus
    branch on this toy, because no terminus-flagged intron|exon face exists there, so the pure
    function carries the falsification instead."""
    from rigel.calibration.messages.transfer_rows import face_is_licensed
    from rigel.calibration.splice_graph import (
        FLAG_TES_NEG,
        FLAG_TES_POS,
        FLAG_TSS_NEG,
        FLAG_TSS_POS,
    )

    assert face_is_licensed(0, True, False, True, False)
    for flag in (FLAG_TSS_POS, FLAG_TSS_NEG, FLAG_TES_POS, FLAG_TES_NEG):
        assert not face_is_licensed(flag, True, False, True, False), f"terminus {flag} must refuse"
    assert not face_is_licensed(0, True, False, True, True), "a strand-set change must refuse"
    assert not face_is_licensed(0, True, True, False, True), "a strand-set change must refuse"
    assert face_is_licensed(0, True, True, True, True), "matching AMBIG sides are licensed"


def test_the_face_map_is_monotone_and_flux_capped():
    """The map lam_e(lam_u) is monotone nondecreasing; a measured spliced density CAPS the
    exon's claimable f_g at the closed-form ceiling; s = 0 degenerates to a pure shift by the
    opportunity ratio."""
    from rigel.calibration.messages.transfer_rows import face_map_lambda

    lam = np.linspace(-10, 10, 401)
    n_u, A_g_b, A_r_b, A_g_e, A_r_e, s = 40.0, 200.0, 210.0, 800.0, 790.0, 3.0
    le = face_map_lambda(lam, n_u, A_g_b, A_r_b, A_g_e, A_r_e, s)
    assert np.all(np.diff(le) >= -1e-12), "the map must be monotone nondecreasing"
    ceiling = np.log(n_u / A_g_b * A_g_e) - np.log(s * A_r_e)
    # the ceiling is a SUPREMUM: the grid's last point sits sigma(+L)-short of it (~5e-5 nats
    # at L = 10), so approach is asserted at that scale and the CAP is asserted hard.
    assert le[-1] == pytest.approx(ceiling, abs=1e-3)
    assert le.max() <= ceiling + 1e-9, "certified flux must CAP the claimable f_g"
    le0 = face_map_lambda(lam, n_u, A_g_b, A_r_b, A_g_e, A_r_e, 0.0)
    shift = np.log((A_g_e / A_g_b) / (A_r_e / A_r_b))
    np.testing.assert_allclose(le0, lam + shift, rtol=0, atol=1e-9)


def test_the_ingredient_width_adds_the_counting_variance():
    """A Gaussian factor transported through an identity-shaped map must widen by exactly the
    face's counting variance trigamma(n_u+1/2) + trigamma(n_s+1/2), the delta-method width — the
    analytic property an inverted or mis-scaled kernel cannot fake."""
    from scipy.special import polygamma

    from rigel.calibration.messages.transfer_rows import transport_row

    lam = np.linspace(-10, 10, 401)
    row = -0.5 * (lam - 0.8) ** 2  # unit-variance Gaussian factor
    n_u, n_s = 9.0, 3.0
    v = float(polygamma(1, n_u + 0.5) + polygamma(1, n_s + 0.5))
    out = transport_row(row, lam, lam.copy(), n_u, n_s)  # identity map isolates the width

    def _var(r):
        w = np.exp(r - r.max())
        w /= w.sum()
        m = w @ lam
        return float(w @ (lam * lam) - m * m)

    assert _var(row) == pytest.approx(1.0, rel=0.02)
    assert _var(out) == pytest.approx(1.0 + v, rel=0.02), "the width must add exactly the variance"
    assert abs(float(out.max())) < 1e-12


def _expected_splice_out_rows(si, ctx, strand, lam):
    """The splice-out boundary rows recomputed INDEPENDENTLY of the policy: for every licensed
    intron|exon face of every exon whose strand channel is live, the exon's own frozen-variance
    strand row read at the GEOMETRIC splice-in map (both components on ``eff_gdna``, the
    spliced density ``S / A_g^b``), marginalised over ``log rho ~ N(log S/U, trigamma(S+1/2) + trigamma(U+1/2))`` on nine
    equal-probability nodes — a second implementation, so a policy bug cannot hide."""
    from scipy.special import polygamma
    from scipy.stats import norm

    from rigel.calibration.splice_graph import (
        FLAG_TES_NEG,
        FLAG_TES_POS,
        FLAG_TSS_NEG,
        FLAG_TSS_POS,
    )

    kappa, od_g, od_r = strand
    term = FLAG_TSS_POS | FLAG_TSS_NEG | FLAG_TES_POS | FLAG_TES_NEG
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    is_intron = ~is_bnd & ~is_exon & (fp | fn)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    n_u = np.asarray(ctx.n_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna, np.float64)
    sc_lo = np.asarray(ctx.sj_count_lo, np.float64).sum(axis=1)
    sc_hi = np.asarray(ctx.sj_count_hi, np.float64).sum(axis=1)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    live = np.asarray(ctx.own_live, bool)
    belief = np.asarray(ctx.belief_fg, np.float64)
    fg = 1.0 / (1.0 + np.exp(-lam))
    nodes = norm.ppf((np.arange(9) + 0.5) / 9.0)
    out = {}
    for e in np.flatnonzero(is_exon):
        if fp[e] == fn[e] or not live[e]:
            continue
        n = cnt[e].sum()
        ks = kappa if fp[e] else 1.0 - kappa
        f_ref = float(np.clip(belief[e], 1e-3, 1 - 1e-3))
        p = 0.5 * fg + ks * (1 - fg)
        p_ref = 0.5 * f_ref + ks * (1 - f_ref)
        var = max(
            n * p_ref * (1 - p_ref)
            + (n * f_ref) ** 2 * 0.25 * od_g
            + (n * (1 - f_ref)) ** 2 * ks * (1 - ks) * od_r,
            1e-9,
        )
        row_e = -0.5 * (cnt[e, 0] - n * p) ** 2 / var
        row_e -= row_e.max()
        if np.ptp(row_e) <= 1e-9:
            continue
        for b, hi in ((left[e], True), (right[e], False)):
            if b < 0 or not is_bnd[b]:
                continue
            i = left[b] if right[b] == e else right[b]
            if i < 0 or not is_intron[i] or (flags[b] & term) or fp[e] != fp[i] or fn[e] != fn[i]:
                continue
            n_s = float((sc_hi if hi else sc_lo)[b])
            if not (n_u[b] > 0 and A_g[b] > 0 and A_g[e] > 0):
                continue
            v = float(polygamma(1, n_s + 0.5) + polygamma(1, n_u[b] + 0.5))
            acc = np.zeros(lam.shape[0])
            for z in nodes:
                s = n_s / A_g[b] * np.exp(z * np.sqrt(v))
                m = np.log(n_u[b] * fg / A_g[b] * A_g[e]) - np.log(
                    (n_u[b] * (1 - fg) / A_g[b] + s) * A_g[e]
                )
                acc += np.exp(np.interp(m, lam, row_e, left=row_e[0], right=row_e[-1]))
            r = np.log(np.maximum(acc / nodes.size, 1e-300))
            r -= r.max()
            if np.ptp(r) <= 1e-9:
                continue
            out.setdefault(int(b), np.zeros(lam.shape[0]))
            out[int(b)] += r
    return out


def test_the_splice_out_row_is_the_count_form_widened_by_the_marginal():
    """The splice-out row's three analytic properties: (i) with equal face and region geometry its
    mode sits at f_b = f_E (U+S)/U, the enrichment ratio having cancelled;
    (ii) more gDNA at the boundary than in the exon whenever S > 0 (the reversed subtraction cannot
    pass); (iii) the marginal over log rho is WIDER at a thin face than at a deep one."""
    from rigel.calibration.messages.transfer_rows import splice_out_row

    lam = np.linspace(-8, 8, 801)
    sig = 1 / (1 + np.exp(-lam))
    f_e = 0.20
    row_e = -0.5 * ((lam - np.log(f_e / (1 - f_e))) / 0.05) ** 2  # a sharp exon belief at f_E

    def _mode_f(r):
        return float(sig[np.argmax(r)])

    def _var(r):
        w = np.exp(r - r.max())
        w /= w.sum()
        m = w @ lam
        return float(w @ (lam * lam) - m * m)

    deep = splice_out_row(row_e, lam, n_u=2000.0, n_s=3000.0, a_g_b=210.0, a_g_e=210.0)
    assert _mode_f(deep) == pytest.approx(f_e * (2000 + 3000) / 2000, abs=0.01)
    assert _mode_f(deep) > f_e
    thin = splice_out_row(row_e, lam, n_u=6.0, n_s=9.0, a_g_b=210.0, a_g_e=210.0)
    assert _var(thin) > 3 * _var(deep), "a thin face must deliver a much wider claim"
    assert abs(float(deep.max())) < 1e-12 and abs(float(thin.max())) < 1e-12
    assert not splice_out_row(row_e, lam, n_u=0.0, n_s=9.0, a_g_b=210.0, a_g_e=210.0).any()


def _expected_boundary_rows(si, ctx, strand, lam):
    """The boundary→intron rows recomputed INDEPENDENTLY of the policy: for every intron|exon pair whose
    boundary admits the intron's SINGLE strand set and whose strand channel the solver declares live
    (``own.tau_lam > 0`` at the boundary — the strand Fisher information itself, there being no
    factory at a boundary), the boundary's own frozen-variance strand row VERBATIM, summed at the
    intron — a second implementation, so a policy bug cannot hide."""
    kappa, od_g, od_r = strand
    pairs, _n = _expected_pairs(si)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    live = np.asarray(ctx.own_live, bool)
    belief = np.asarray(ctx.belief_fg, np.float64)
    fg = 1.0 / (1.0 + np.exp(-lam))
    out = {}
    for b, i in pairs:
        if fp[b] != fp[i] or fn[b] != fn[i] or fp[b] == fn[b] or not live[b]:
            continue
        n = cnt[b].sum()
        ks = kappa if fp[b] else 1.0 - kappa
        f_ref = float(np.clip(belief[b], 1e-9, 1 - 1e-9))
        p = 0.5 * fg + ks * (1 - fg)
        p_ref = 0.5 * f_ref + ks * (1 - f_ref)
        var = max(
            n * p_ref * (1 - p_ref)
            + (n * f_ref) ** 2 * 0.25 * od_g
            + (n * (1 - f_ref)) ** 2 * ks * (1 - ks) * od_r,
            1e-9,
        )
        row = -0.5 * (cnt[b, 0] - n * p) ** 2 / var
        row -= row.max()
        if np.ptp(row) <= 1e-9:
            continue
        out.setdefault(int(i), np.zeros(lam.shape[0]))
        out[int(i)] += row
    return out


def test_the_boundary_strand_licence_requires_one_shared_strand():
    """The pair predicate, gated DIRECTLY — the integration toy cannot falsify it (every pair there
    shares one strand), so the pure function carries the falsification: one strand, shared by both
    flanks, is licensed; an AMBIG pair (whose strand split constrains the tilt, never f_g), a
    strand-set change and an empty boundary refuse."""
    from rigel.calibration.messages.transfer_rows import boundary_shares_strand

    assert boundary_shares_strand(True, False, True, False)
    assert boundary_shares_strand(False, True, False, True)
    assert not boundary_shares_strand(True, True, True, True), "an AMBIG pair must refuse"
    assert not boundary_shares_strand(True, False, True, True), "a strand-set change must refuse"
    assert not boundary_shares_strand(False, False, True, False), "an empty boundary must refuse"


def test_the_intron_face_carries_the_pair_identity_and_the_face_map(sweep_inputs):
    """The intron's face, as claims and rules: its own claim is its factory profile and its rule into
    each boundary is the identity (FORWARD); the rule from a licensed face into the exon, applied to
    the intron's claim and summed with the edge's level rule, equals the independently recomputed
    splice-in + edge rows of that exon; an unlicensed face has no rule into the exon."""
    from rigel.calibration.messages.transfer import FORWARD, TransferPolicy
    from rigel.calibration.simplex_logodds import _logodds_grid

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    pairs, _n = _expected_pairs(sweep_inputs)
    assert pairs, "the toy must carry at least one intron|exon pair or this gate proves nothing"
    ctx = _ctx_of(sweep_inputs)
    prepared = _prepared(TransferPolicy(), ctx)
    src = ctx.factory_rows
    lam, _ = _logodds_grid(n_grid, window)
    for b, j in pairs:
        assert prepared.faces.kind_at(j, b) == FORWARD, f"pair ({j}, {b}) is not FORWARD"
        np.testing.assert_array_equal(prepared.own[j], src[j] - src[j].max())
    exon_rows = _expected_exon_rows(sweep_inputs, ctx, src, lam)
    assert exon_rows, "the toy must license at least one exon face or this gate proves nothing"
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    intron = _intron_mask(ctx)
    for e in np.flatnonzero(is_exon):
        acc = None
        for b in (left[e], right[e]):
            if b < 0 or not prepared.faces.has(b, e):
                continue
            i = left[b] if right[b] == e else right[b]
            claim = prepared.own[i] if (i >= 0 and intron[i]) else prepared.own[b]
            if claim is None:
                continue
            r = prepared.faces.apply(b, e, claim, None)
            if r is not None and np.ptp(r) > 1e-9:
                acc = (r - r.max()) if acc is None else acc + (r - r.max())
        if int(e) in exon_rows:
            assert acc is not None, f"exon {e} has no rule where the recompute expects rows"
            np.testing.assert_allclose(acc - acc.max(), exon_rows[int(e)], rtol=0, atol=1e-12)
        else:
            assert acc is None, f"exon {e} received a rule it is not licensed for"


def test_the_exon_and_boundary_own_claims_are_the_strand_rows_and_their_rules_the_maps(
    sweep_inputs,
):
    """The two strand-borne rules, as claims and rules: at every licensed face of a strand-live exon the rule
    exon → boundary applied to the exon's claim, summed per boundary, equals the independently
    recomputed splice-out rows; at every intron|exon pair sharing one strand with a live boundary,
    the rule boundary → intron is the identity and the boundary's claim, summed per intron, equals
    the independently recomputed strand rows."""
    from rigel.calibration.messages.transfer import FORWARD
    from rigel.calibration.simplex_logodds import _logodds_grid

    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    prepared = _prepared(pol, ctx)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    expected = _expected_splice_out_rows(sweep_inputs, ctx, strand, lam)
    assert expected, (
        "the toy must carry a strand-live exon with a licensed face or this gate proves nothing"
    )
    intron = _intron_mask(ctx)
    got = {}
    for e in np.flatnonzero(is_exon):
        if prepared.own[e] is None:
            continue
        for b in (left[e], right[e]):
            if b < 0 or not prepared.faces.has(e, b):
                continue
            other = left[b] if right[b] == e else right[b]
            if other < 0 or not intron[other]:
                continue  # the terminus and alternative-splice rules leave through exon|exon faces: their own gates
            r = prepared.faces.apply(e, b, prepared.own[e], None)
            if r is not None and np.ptp(r) > 1e-9:
                got[int(b)] = got.get(int(b), 0.0) + (r - r.max())
    assert set(got) == set(expected), set(got) ^ set(expected)
    for b in expected:
        np.testing.assert_allclose(
            got[b] - got[b].max(), expected[b] - expected[b].max(), rtol=0, atol=1e-10
        )
    expected_i = _expected_boundary_rows(sweep_inputs, ctx, strand, lam)
    assert expected_i, (
        "the toy must carry a strand-live intron|exon pair or this gate proves nothing"
    )
    pairs, _n = _expected_pairs(sweep_inputs)
    got_i = {}
    for b, i in pairs:
        if not prepared.faces.has(b, i) or prepared.own[b] is None:
            continue
        assert prepared.faces.kind_at(b, i) == FORWARD, (
            f"the boundary → intron rule at ({b}, {i}) is not the identity"
        )
        got_i[int(i)] = got_i.get(int(i), 0.0) + prepared.own[b]
    assert set(got_i) == set(expected_i), set(got_i) ^ set(expected_i)
    for i in expected_i:
        np.testing.assert_allclose(
            got_i[i] - got_i[i].max(), expected_i[i] - expected_i[i].max(), rtol=0, atol=1e-10
        )


# ── the edge's one-sided level and the terminus rules ────────────────────────────────────────────


def test_the_edge_level_is_one_sided_and_a_zero_count_is_vacuous():
    """The edge's level, in the form the ladder kept: below the edge's level the count's exact
    Poisson — the exon has at least the edge's gDNA density; nothing above it, because no local
    witness prices capture's enrichment of the interior; and a zero count is vacuous, because
    darkness under capture is not absence."""
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
    """The terminus rows recomputed INDEPENDENTLY of the policy. At every exon|exon boundary carrying exactly one
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
    live = np.asarray(ctx.own_live, bool)
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
        if live[o]:
            acc += splice_out_row(strand_row(o), lam, n_u[b], n_s[b], A_g[b], A_g[o])
        if np.ptp(acc) > 1e-9:
            at_b[int(b)] = acc
        if live[b]:
            le = face_map_lambda(lam, n_u[b], A_g[b], A_g[b], A_g[o], A_g[o], n_s[b] / A_g[b])
            at_o.setdefault(int(o), np.zeros(lam.shape[0]))
            at_o[int(o)] += transport_row(strand_row(b), lam, le, n_u[b], n_s[b])
    return at_b, at_o


def _item5_slots(ctx):
    """The terminus rules' delivery sites — every exon|exon boundary carrying one terminus direction and no sj
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
    pointing both ways or no terminus at all give no side. A splice junction sharing the boundary does
    NOT change the side, which is the sj+terminus case: the four ladder families resolve as their
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
    from rigel.calibration.messages.transfer import LEVEL

    prep = _prepared(pol, ctx2)
    face = prep.faces.at(b, i)
    assert face.kind == LEVEL
    v_with = face.var
    # the plain form: the same boundary with its flux zeroed
    sjc = np.asarray(ctx.sj_count, np.float64).copy()
    sjc[b] = 0.0
    prep0 = _prepared(pol, dataclasses.replace(ctx2, sj_count=sjc))
    face0 = prep0.faces.at(b, i)
    assert face0.kind == LEVEL
    v_plain = face0.var
    assert v_with < v_plain, (v_with, v_plain)
    # the junction rules leave the face: the splice-in map into the inside and the splice-out map out of it are gone
    prep_before = _prepared(pol, ctx)
    assert prep_before.faces.has(b, i) and prep_before.faces.kind_at(b, i) != LEVEL


def _expected_level_rows(si, ctx, strand, lam):
    """The level rule recomputed INDEPENDENTLY of the policy: at every terminus boundary with an inside
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
    live = np.asarray(ctx.own_live, bool)
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
        if live[b] and live[i]:
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
        if live[b] and fp[b] != fn[b]:
            row = level_row(_strand_row_of(ctx, strand, lam, b), lam, m, v)
        else:
            row = level_bound_row(lam, d_b, A_g[i], n_u[i], v)
        if np.ptp(row) > 1e-9:
            out.setdefault(int(i), np.zeros(lam.shape[0]))
            out[int(i)] += row
    return out, served


def test_the_level_map_keeps_the_level_and_the_bound_is_vacuous_below_the_total():
    """The level rule's arithmetic, gated directly: through the level-kept map a boundary whose share is
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
    """The terminus rules, as claims and rules: at every exon|exon terminus boundary the rule outside exon → boundary
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
    prepared = _prepared(pol, ctx)
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    n_u = np.asarray(ctx.n_slot, np.float64)
    n_s = np.asarray(ctx.spliced_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna, np.float64)
    live = np.asarray(ctx.own_live, bool)
    sites = _item5_slots(ctx)
    served = 0
    for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
        if not (is_exon[left[b]] and is_exon[right[b]]):
            continue
        outs = [
            (s, i)
            for (s, i) in prepared.faces.pairs()
            if (s == b and is_exon[i]) or (i == b and is_exon[s])
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
        if live[o] and prepared.own[o] is not None:
            got = prepared.faces.apply(o, b, prepared.own[o], None)
            np.testing.assert_allclose(got - got.max(), want - want.max(), rtol=0, atol=1e-10)
        if live[b]:
            le = face_map_lambda(lam, n_u[b], A_g[b], A_g[b], A_g[o], A_g[o], n_s[b] / A_g[b])
            want_o = transport_row(_strand_row_of(ctx, strand, lam, b), lam, le, n_u[b], n_s[b])
            got_o = prepared.faces.apply(b, o, prepared.own[b], None)
            np.testing.assert_allclose(
                got_o - got_o.max(), want_o - want_o.max(), rtol=0, atol=1e-10
            )
    assert served >= 2, (
        "the toy must carry two served terminus boundaries or this gate proves nothing"
    )
    cleared = _prepared(pol, _terminus_flags_cleared(ctx))
    for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
        if is_exon[left[b]] and is_exon[right[b]]:
            assert not any(s == b or i == b for (s, i) in cleared.faces.pairs()), (
                f"rules survive at {b} with no terminus"
            )


def test_the_level_rule_serves_every_terminus_inside_from_the_measurement_alone(sweep_inputs):
    """The level rule, as a rule: at every terminus boundary with an inside exon —
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
    prepared = _prepared(pol, ctx)
    expected, served = _expected_level_rows(sweep_inputs, ctx, strand, lam)
    assert len(served) >= 2, (
        "the toy must carry two served terminus pairs or this gate proves nothing"
    )
    spike = -0.5 * ((lam - 2.0) / 0.1) ** 2  # a held imputation that must not cross a level face
    for b, i, _kind in served:
        assert prepared.faces.has(b, i), f"no level rule at terminus pair ({b}, {i})"
        got = prepared.faces.apply(b, i, prepared.own[b], None)
        np.testing.assert_allclose(
            got - got.max(), expected[i] - expected[i].max(), rtol=0, atol=1e-10
        )
        np.testing.assert_array_equal(
            prepared.faces.apply(b, i, prepared.own[b], spike),
            got,
            err_msg="what is held crossed a level face",
        )
        bound = prepared.faces.apply(b, i, None, spike)
        assert bound is not None and np.ptp(bound) > 0.0 and np.all(bound <= 0.0), (
            "no upper bound without a claim"
        )
    cleared = _prepared(pol, _terminus_flags_cleared(ctx))
    for b, i, _kind in served:
        assert not cleared.faces.has(b, i), f"a level rule survives at ({b}, {i}) with no terminus"
    b0, i0, _k = served[0]
    b1, i1, _k = served[-1]
    cnt = np.asarray(ctx.unspliced_count, np.float64).copy()
    cnt[i1] = cnt[i1] * 3.0 + 7.0
    other = _prepared(pol, _dc.replace(ctx, unspliced_count=cnt, n_slot=cnt.sum(axis=1)))
    np.testing.assert_array_equal(
        other.faces.apply(b0, i0, other.own[b0], None),
        prepared.faces.apply(b0, i0, prepared.own[b0], None),
    )


# ── the alternative splice site: both flanks, priced per pair ────────────────────────────────────


def test_the_junction_flanks_read_the_flag_kind_alone():
    """The flank predicate, gated DIRECTLY: a DONOR bit puts the intron to the right on either strand
    (the flags are genomic-order), an ACCEPTOR bit to the left; a terminus on the boundary, no junction,
    or junctions both ways give no flanks."""
    from rigel.calibration.messages.transfer_rows import junction_flanks
    from rigel.calibration.splice_graph import (
        FLAG_ACCEPTOR_NEG,
        FLAG_ACCEPTOR_POS,
        FLAG_DONOR_NEG,
        FLAG_DONOR_POS,
        FLAG_TSS_POS,
    )

    assert junction_flanks(FLAG_DONOR_POS, 7, 9) == (9, 7)
    assert junction_flanks(FLAG_DONOR_NEG, 7, 9) == (9, 7)
    assert junction_flanks(FLAG_ACCEPTOR_POS, 7, 9) == (7, 9)
    assert junction_flanks(FLAG_ACCEPTOR_NEG, 7, 9) == (7, 9)
    assert junction_flanks(FLAG_DONOR_POS | FLAG_TSS_POS, 7, 9) == (None, None), (
        "sj+terminus: its own item"
    )
    assert junction_flanks(0, 7, 9) == (None, None)
    assert junction_flanks(FLAG_DONOR_POS | FLAG_ACCEPTOR_NEG, 7, 9) == (None, None), (
        "junctions both ways"
    )


def _expected_alt_splice_rows(si, ctx, strand, lam):
    """The alternative-splice rows recomputed INDEPENDENTLY of the policy on the patched context: at
    each junction boundary the two flanks' own strand rows through the splice-out map (E with
    S_b + F, C with S_b) into the boundary, and the boundary's own row through the face map into
    each flank, each blurred by the pair's OWN disagreement beyond counting — per pair, nothing
    pooled. Keyed by slot; the pair widths returned beside."""
    from rigel.calibration.messages.transfer_rows import (
        blur_row,
        boundary_shares_strand,
        face_map_lambda,
        junction_flanks,
        splice_out_row,
        transport_row,
    )

    kappa, od_g, od_r = strand
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    n_u = np.asarray(ctx.n_slot, np.float64)
    n_s = np.asarray(ctx.spliced_slot, np.float64)
    flux = np.asarray(ctx.sj_count, np.float64).sum(axis=1)
    A_g = np.asarray(ctx.eff_gdna, np.float64)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    live = np.asarray(ctx.own_live, bool)
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

    served = []
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0 or not (is_exon[lo] and is_exon[hi]):
            continue
        c_side, e_side = junction_flanks(flags[b], lo, hi)
        if c_side is None or not (n_u[b] > 0 and A_g[b] > 0):
            continue
        for x, s_out, kind in ((e_side, n_s[b] + flux[b], "E"), (c_side, n_s[b], "C")):
            if boundary_shares_strand(fp[b], fn[b], fp[x], fn[x]) and A_g[x] > 0:
                served.append((int(b), int(x), float(s_out), kind))
    width = {}
    for b, x, s_out, _kind in served:
        if not (live[b] and live[x]):
            continue
        lo_v = []
        for y in (b, x):
            n = cnt[y].sum()
            ks = kappa if fp[y] else 1.0 - kappa
            p = cnt[y, 0] / n
            f = (p - ks) / (0.5 - ks)
            if not 0.0 < f < 1.0:
                break
            lo_v.append((np.log(f / (1 - f)), p * (1 - p) / n / (p - ks) ** 2 / (1 - f) ** 2))
        if len(lo_v) < 2:
            continue
        v_ratio = s_out / (n_u[b] * (n_u[b] + s_out)) if s_out > 0 else 0.0
        d = lo_v[0][0] - lo_v[1][0] - np.log((n_u[b] + s_out) / n_u[b])
        width[(b, x)] = max(0.0, d * d - (lo_v[0][1] + lo_v[1][1] + v_ratio))
    out = {}
    for b, x, s_out, _kind in served:
        w = width.get((b, x), 0.0)
        if live[x] and fp[x] != fn[x]:
            row = splice_out_row(strand_row(x), lam, n_u[b], s_out, A_g[b], A_g[x])
            if np.ptp(row) > 1e-9:
                out.setdefault(b, np.zeros(lam.shape[0]))
                out[b] += blur_row(row, lam, w)
        if live[b]:
            le = face_map_lambda(lam, n_u[b], A_g[b], A_g[b], A_g[x], A_g[x], s_out / A_g[b])
            row = transport_row(strand_row(b), lam, le, n_u[b], s_out)
            if np.ptp(row) > 1e-9:
                out.setdefault(x, np.zeros(lam.shape[0]))
                out[x] += blur_row(row, lam, w)
    return out, width


def test_the_alt_splice_rules_carry_both_flanks_with_the_pair_width(sweep_inputs):
    """The alternative splice site, as claims and rules on the patched toy: the rule flank → junction boundary applied to the
    flank's claim, summed over the two flanks, equals the independently recomputed rows at the
    boundary; the rule boundary → flank applied to the boundary's claim equals the recomputed row at
    each flank — every message blurred by its own pair's disagreement beyond counting."""
    from rigel.calibration.simplex_logodds import _logodds_grid

    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _with_alt_splice_sites(_ctx_of(sweep_inputs))
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    prepared = _prepared(pol, ctx)
    expected, widths = _expected_alt_splice_rows(sweep_inputs, ctx, strand, lam)
    assert expected and any(w > 0.0 for w in widths.values()), (
        "the patched toy must carry served junctions with a live pair width or this gate proves nothing"
    )
    is_bnd = np.asarray(ctx.is_boundary, bool)
    got = {}
    for s, d in prepared.faces.pairs():
        if not (
            (is_bnd[s] and s in widths_keys(widths)) or (is_bnd[d] and d in widths_keys(widths))
        ):
            continue
        if prepared.own[s] is None:
            continue
        r = prepared.faces.apply(s, d, prepared.own[s], None)
        if r is not None and np.ptp(r) > 1e-9:
            got[int(d)] = got.get(int(d), 0.0) + (r - r.max())
    for slot, want in expected.items():
        assert slot in got, f"slot {slot} has no live item-7 rule"
        np.testing.assert_allclose(
            got[slot] - got[slot].max(), want - want.max(), rtol=0, atol=1e-10
        )


def widths_keys(widths):
    """The junction boundaries the recompute served (its width keys are ``(boundary, flank)`` pairs)."""
    return {int(k[0]) if isinstance(k, tuple) else int(k) for k in widths}
