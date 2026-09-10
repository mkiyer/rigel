"""THE INTRON|EXON FACE (`transfer._splice_faces`): the licence refuses unmeasured population
changes; the splice-in face map is monotone and flux-capped and its width is exactly the face's
counting variance; the splice-out row is the count form marginalised over the ratio's counting;
the boundary → intron licence needs one shared strand; and on the live toy the intron's pair
identity, the face map into the exon and the strand rows at exactly the licensed slots, each
recomputed by a second implementation so a policy bug cannot hide."""

from __future__ import annotations

import numpy as np
import pytest

from _transfer_harness import (
    _ctx_of,
    _expected_pairs,
    _full_policy,
    _intron_mask,
    _live_provider,
    _strand_of,
    capture_sweep_inputs,
)


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    return capture_sweep_inputs(tmp_path_factory)


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
            # the intergenic|exon EDGE's level, lower bound only (owner ruling 2026-09-02)
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
    """The licence predicate, gated DIRECTLY — the integration gate cannot falsify the terminus
    branch on this toy (no terminus-flagged intron|exon face exists there; watched 2026-09-01),
    so the pure function carries the falsification instead."""
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
    """A Gaussian factor transported through an identity-shaped map must widen by EXACTLY the
    face's counting variance trigamma(n_u+1/2) + trigamma(n_s+1/2) (the delta-method width) —
    the analytic property an inverted or mis-scaled kernel cannot fake (the forward hop's lesson,
    watched 2026-09-01)."""
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
    """The item-1 boundary rows recomputed INDEPENDENTLY of the policy: for every licensed
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
    tau = np.asarray(ctx.own.tau_lam, np.float64)
    belief = np.asarray(ctx.belief_fg, np.float64)
    fg = 1.0 / (1.0 + np.exp(-lam))
    nodes = norm.ppf((np.arange(9) + 0.5) / 9.0)
    out = {}
    for e in np.flatnonzero(is_exon):
        if fp[e] == fn[e] or not tau[e] > 0.0:
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
    mode sits at f_b = f_E (U+S)/U — the owner's arithmetic with the enrichment ratio cancelled;
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
    """The item-2 intron rows recomputed INDEPENDENTLY of the policy: for every intron|exon pair whose
    boundary admits the intron's SINGLE strand set and whose strand channel the solver declares live
    (``own.tau_lam > 0`` at the boundary — the strand Fisher information itself, there being no
    factory at a boundary), the boundary's own frozen-variance strand row VERBATIM, summed at the
    intron — a second implementation, so a policy bug cannot hide."""
    kappa, od_g, od_r = strand
    pairs, _n = _expected_pairs(si)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    tau = np.asarray(ctx.own.tau_lam, np.float64)
    belief = np.asarray(ctx.belief_fg, np.float64)
    fg = 1.0 / (1.0 + np.exp(-lam))
    out = {}
    for b, i in pairs:
        if fp[b] != fp[i] or fn[b] != fn[i] or fp[b] == fn[b] or not tau[b] > 0.0:
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
    """RUNGS 1–3 as claims and rules: an intron's own claim is its factory profile and its rule into
    each boundary is the identity (FORWARD); the rule from a licensed face into the exon, applied to
    the intron's claim and summed with the edge's level rule, equals the independently recomputed
    splice-in + edge rows of that exon; an unlicensed face has no rule into the exon."""
    from rigel.calibration.messages.transfer import TransferPolicy, _forward
    from rigel.calibration.simplex_logodds import _logodds_grid

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    pairs, _n = _expected_pairs(sweep_inputs)
    assert pairs, "the toy must carry at least one intron|exon pair or this gate proves nothing"
    ctx = _ctx_of(sweep_inputs)
    prepared = TransferPolicy(provider).prepare(ctx)
    src = provider(n_grid, window)
    lam, _ = _logodds_grid(n_grid, window)
    for b, j in pairs:
        assert prepared.rule.get((int(j), int(b))) is _forward, f"pair ({j}, {b}) is not FORWARD"
        np.testing.assert_array_equal(prepared.own[j], src[j] - src[j].max())
    exon_rows = _expected_exon_rows(sweep_inputs, ctx, src, lam)
    assert exon_rows, "the toy must license at least one exon face or this gate proves nothing"
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    intron = _intron_mask(ctx)
    for e in np.flatnonzero(is_exon):
        acc = None
        for b in (left[e], right[e]):
            fn = prepared.rule.get((int(b), int(e)))
            if fn is None:
                continue
            i = left[b] if right[b] == e else right[b]
            claim = prepared.own[i] if (i >= 0 and intron[i]) else prepared.own[b]
            if claim is None:
                continue
            r = fn(claim, None)
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
    """ITEMS 1 and 2 as claims and rules: at every licensed face of a strand-live exon the rule
    exon → boundary applied to the exon's claim, summed per boundary, equals the independently
    recomputed splice-out rows; at every intron|exon pair sharing one strand with a live boundary,
    the rule boundary → intron is the identity and the boundary's claim, summed per intron, equals
    the independently recomputed strand rows."""
    from rigel.calibration.messages.transfer import _forward
    from rigel.calibration.simplex_logodds import _logodds_grid

    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    prepared = pol.prepare(ctx)
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
            fn = prepared.rule.get((int(e), int(b)))
            if fn is None or b < 0:
                continue
            other = left[b] if right[b] == e else right[b]
            if other < 0 or not intron[other]:
                continue  # the terminus and alternative-splice rules leave through exon|exon faces: their own gates
            r = fn(prepared.own[e], None)
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
        fn = prepared.rule.get((int(b), int(i)))
        if fn is None or prepared.own[b] is None:
            continue
        assert fn is _forward, f"the boundary → intron rule at ({b}, {i}) is not the identity"
        got_i[int(i)] = got_i.get(int(i), 0.0) + prepared.own[b]
    assert set(got_i) == set(expected_i), set(got_i) ^ set(expected_i)
    for i in expected_i:
        np.testing.assert_allclose(
            got_i[i] - got_i[i].max(), expected_i[i] - expected_i[i].max(), rtol=0, atol=1e-10
        )
