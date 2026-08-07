"""Gates for the ``eta`` transfer derivation (`_eta_reference`).

One test per property the rebuild actually rests on, so a regression names itself.

⛔ **A2 — each of these was run against a deliberately-broken kernel and watched to fire**, and the
perturbation that matters is named in each docstring. The single most important one:
:func:`test_eta_crosses_as_the_identity_and_the_shares_do_not` fails under "composition transfer is the
identity on ``lambda``", which is the error the first draft of the derivation actually made.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import digamma, polygamma

from _eta_reference import (
    GDNA,
    RNA_NEG,
    RNA_POS,
    cross_composition,
    eta_from_lambda,
    lambda_from_eta,
    log_rate_posterior,
    population_set,
    rebuild_densities,
    residual_eta_scalar,
    share_from_density,
    shares_from_lambda,
    tilt_angle,
)

#: (E_g, E_r) at source -> destination. Every pair has a DIFFERENT opportunity ratio at the two ends,
#: which is the only situation in which eta and lambda are distinguishable at all.
OPPORTUNITY_PAIRS = [
    ((100.0, 200.0), (250.0, 100.0)),
    ((3000.0, 1500.0), (254.0, 254.0)),
    ((50.0, 400.0), (900.0, 300.0)),
]


# ── the frame-free coordinate ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("src,dst", OPPORTUNITY_PAIRS)
def test_eta_crosses_as_the_identity_and_the_shares_do_not(src, dst):
    """⭐⭐⭐ THE CENTREPIECE. What transfers is the DENSITY RATIO, not the share.

    "gDNA is uniform" and "the same transcripts run through both slots" are statements about DENSITIES.
    Shares are derived, and they move whenever the two slots' opportunity ratio differs — which between a
    NODE and an EDGE it always does.

    ⛔ The second assertion is not decoration. A gate that only checked ``eta_s == eta_d`` would also pass
    under "lambda crosses unchanged" whenever the two ratios happen to agree, and would ACTIVELY ENSHRINE
    the first draft's error on a fixture where they do not. Asserting the shares are NOT equal is what
    makes this a gate on the right invariant.
    """
    (eg_s, er_s), (eg_d, er_d) = src, dst
    rho_g, rho_r = 3.0, 7.0  # one physical pair of densities, seen from two frames
    lam_s = np.log(rho_g * eg_s / (rho_r * er_s))
    lam_d = np.log(rho_g * eg_d / (rho_r * er_d))

    eta_s = eta_from_lambda(lam_s, eg_s, er_s)
    eta_d = eta_from_lambda(lam_d, eg_d, er_d)
    assert eta_s == pytest.approx(np.log(rho_g / rho_r), abs=1e-12)
    assert eta_d == pytest.approx(eta_s, abs=1e-12)

    # the delivered lambda differs from the source's by EXACTLY the geometric constant
    geo = np.log(eg_d / er_d) - np.log(eg_s / er_s)
    assert lambda_from_eta(eta_s, eg_d, er_d) == pytest.approx(lam_s + geo, abs=1e-12)

    # ...and the shares genuinely move, so the identity above is a real claim about a real difference.
    phi_s, _ = shares_from_lambda(lam_s)
    phi_d, _ = shares_from_lambda(lam_d)
    assert abs(phi_s - phi_d) > 1e-3, "fixture is degenerate: the two frames agree on the share"


def test_eta_lambda_roundtrip_is_exact():
    """The conversion is a known constant, so it must be invertible to machine precision in both
    directions and at every opportunity ratio — including ``E_g == E_r``, where it is the identity."""
    rng = np.random.default_rng(0)
    lam = rng.normal(0.0, 3.0, 512)
    eg = np.exp(rng.normal(6.0, 2.0, 512))
    er = np.exp(rng.normal(6.0, 2.0, 512))
    back = lambda_from_eta(eta_from_lambda(lam, eg, er), eg, er)
    assert np.max(np.abs(back - lam)) < 1e-12


# ── the mass identity ───────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("src,dst", OPPORTUNITY_PAIRS)
def test_mass_identity_holds_by_construction(src, dst):
    """⭐⭐ ``sum_c rho_c E_c = M`` EXACTLY, at the destination, for ANY eta.

    This is the whole reason the mass pin has nothing left to restore: rebuild from eta with the
    destination's own ``(M, E_g, E_r)`` and the identity is arithmetic, not something to be enforced.

    ⛔ Perturbation: dividing ``rho_R`` by ``E_g`` instead of ``E_r`` (the classic frame slip) breaks this
    at 1e-2 on every pair, while leaving the eta-invariance gate above completely green.
    """
    (_eg_s, _er_s), (eg_d, er_d) = src, dst
    mass = 4321.0
    for eta in (-8.0, -1.0, 0.0, 0.37, 5.0, 20.0):
        rho_g, rho_r = rebuild_densities(eta, mass, eg_d, er_d)
        assert rho_g * eg_d + rho_r * er_d == pytest.approx(mass, rel=1e-10)


def test_mass_identity_holds_vectorised_over_a_whole_chain():
    """The same identity over many slots at once — the form the sweep actually calls it in."""
    rng = np.random.default_rng(7)
    n = 4096
    eta = rng.normal(0.0, 6.0, n)
    mass = np.exp(rng.normal(4.0, 3.0, n))
    eg = np.exp(rng.normal(6.0, 2.0, n))
    er = np.exp(rng.normal(6.0, 2.0, n))
    rho_g, rho_r = rebuild_densities(eta, mass, eg, er)
    assert np.max(np.abs(rho_g * eg + rho_r * er - mass) / mass) < 1e-10


# ── the reference constants ψ is built on ───────────────────────────────────────────────────────────


def test_jeffreys_reference_constants():
    """``psi0(1/2) = -gamma - 2 log 2`` and ``psi1(1/2) = pi^2/2``, to 15 digits.

    These are not tunables — they are the digamma and trigamma of the Jeffreys 1/2 the solver's reference
    is already built on, and every number in the rebuild is one of them or a geometric quantity.
    """
    assert digamma(0.5) == pytest.approx(-np.euler_gamma - 2.0 * np.log(2.0), abs=1e-15)
    assert digamma(0.5) == pytest.approx(-1.9635100260214235, abs=1e-15)
    assert polygamma(1, 0.5) == pytest.approx(np.pi**2 / 2.0, abs=1e-15)
    assert polygamma(1, 0.5) == pytest.approx(4.934802200544679, abs=1e-15)


def test_var_lambda_under_the_beta_half_half_reference_is_pi_squared():
    """``Var(lambda)`` under the reference prior ``f_g ~ Beta(1/2, 1/2)`` is ``2 psi1(1/2) = pi^2``.

    ``lambda = logit(f_g)``, and for ``f ~ Beta(a, b)`` the logit has variance ``psi1(a) + psi1(b)``.
    Checked against a direct high-resolution quadrature so the closed form is not merely restated.
    """
    closed = 2.0 * polygamma(1, 0.5)
    assert closed == pytest.approx(np.pi**2, abs=1e-13)
    assert closed == pytest.approx(9.86960440108936, abs=1e-12)

    # quadrature on the logit axis: f = sigmoid(l), Beta(1/2,1/2) density times the Jacobian.
    lam = np.linspace(-60.0, 60.0, 2_000_001)
    logw = -0.5 * np.log1p(np.exp(-lam)) - 0.5 * np.log1p(np.exp(lam))  # log[f^-.5 (1-f)^-.5 * f(1-f)]
    w = np.exp(logw - logw.max())
    w /= np.trapezoid(w, lam)
    assert np.trapezoid(w * lam, lam) == pytest.approx(0.0, abs=1e-8)
    assert np.trapezoid(w * lam * lam, lam) == pytest.approx(closed, rel=1e-6)


def test_zero_count_claim_is_finite_and_length_sets_only_the_location():
    """⭐ A zero count is a MEASUREMENT of a density, and its relative confidence does not depend on the
    length: ``psi1(1/2)`` is the log-variance at ``a = 0`` whatever ``E`` is, while the location moves by
    ``-log E``. 38.4 kb and 82 b are 470x apart in the claim and identical in confidence.

    ⛔ Perturbation: the shipped ``log(a/E)`` location is ``-inf`` here, which is what pairs a finite
    variance with an impossible mode — this gate fires on it immediately.
    """
    loc_long, var_long = log_rate_posterior(0.0, 38_400.0)
    loc_short, var_short = log_rate_posterior(0.0, 82.0)
    assert np.isfinite(loc_long) and np.isfinite(loc_short)
    assert var_long == pytest.approx(polygamma(1, 0.5), abs=1e-15)
    assert var_short == pytest.approx(var_long, abs=1e-15)
    assert np.exp(loc_long) == pytest.approx(np.exp(digamma(0.5)) / 38_400.0, rel=1e-12)
    assert loc_short - loc_long == pytest.approx(np.log(38_400.0 / 82.0), abs=1e-12)
    # ...and no exact zero is ever formed, so a multiplicative transport has nothing it cannot carry.
    assert np.exp(loc_long) > 0.0


def test_log_rate_variance_falls_monotonically_with_the_count():
    """More fragments, tighter claim — and never the ``1/n = inf`` that made a zero-mass slot a relay
    BARRIER. The finite value at ``a = 0`` is the whole point."""
    a = np.arange(0.0, 64.0)
    _loc, var = log_rate_posterior(a, 1000.0)
    assert np.all(np.diff(var) < 0.0)
    assert np.isfinite(var[0]) and var[0] == pytest.approx(polygamma(1, 0.5), abs=1e-15)


# ── AXIOM 0, made executable ────────────────────────────────────────────────────────────────────────


def test_population_set_never_exceeds_three_members():
    """⛔ THREE POPULATIONS AND THERE IS NO FOURTH. Over every combination of the two annotation bits and
    the two opportunities, ``|T| <= 3`` — and the indicator has exactly three columns, so a fourth
    population is not expressible rather than merely unlikely."""
    fp, fn, eg, er = (
        np.array(v)
        for v in zip(
            *[
                (p, n, g, r)
                for p in (False, True)
                for n in (False, True)
                for g in (0.0, 900.0)
                for r in (0.0, 700.0)
            ],
            strict=True,
        )
    )
    t = population_set(fp, fn, eg, er)
    assert t.shape[-1] == 3
    assert np.all(t.sum(axis=-1) <= 3)
    assert t.sum(axis=-1).max() == 3, "fixture never reaches |T| = 3; the bound is untested"


def test_population_set_is_a_read_of_two_bits_and_two_opportunities():
    """gDNA is admissible wherever it has opportunity — there is no gDNA analogue of a forbidden strand —
    and a slot shorter than a fragment (``E = 0``) holds no population at all."""
    t = population_set([True, True, False], [False, True, False], [10.0, 10.0, 0.0], [10.0, 10.0, 0.0])
    assert list(t[0]) == [True, True, False]
    assert list(t[1]) == [True, True, True]
    assert list(t[2]) == [False, False, False], "no opportunity => no population"

    # gDNA needs no free strand, and an RNA strand needs BOTH its bit and RNA opportunity.
    t = population_set([True], [True], [500.0], [0.0])
    assert list(t[0]) == [True, False, False]


# ── the transfer rule ───────────────────────────────────────────────────────────────────────────────


def test_matched_populations_transfer_eta_untouched():
    """``T_s == T_d`` is the ``|N| = 0`` case: the message is complete and eta is delivered verbatim,
    whatever the two slots' geometry does."""
    t = population_set([True], [False], [100.0], [200.0])
    eta = np.array([-1.234])
    out, determined = cross_composition(
        eta, t, t, rho_g_src=np.array([0.05]), mass_dst=np.array([900.0]),
        eff_g_dst=np.array([250.0]), eff_r_dst=np.array([100.0]),
    )
    assert out[0] == pytest.approx(-1.234, abs=1e-12)
    assert bool(determined[0])


def test_a_newly_active_strand_takes_the_residual_and_is_determined():
    """⭐ ``T_s < T_d`` with ONE newly-active population: the shared gDNA density crosses UNCHANGED, the
    destination's own count fixes the residual, and the split is DETERMINED — one unknown, one equation.
    No interval, no free parameter, no scale to integrate out.

    ⛔ Perturbation: letting the shared density be re-scaled by any ratio of totals (the retired reframe
    ``r``) moves the delivered eta by that ratio and this gate fires.
    """
    t_src = population_set([False], [False], [400.0], [0.0])  # {g}
    t_dst = population_set([True], [False], [250.0], [100.0])  # {g, R+}
    rho_g, mass, eg_d, er_d = 0.02, 1000.0, 250.0, 100.0
    out, determined = cross_composition(
        np.array([0.0]), t_src, t_dst, rho_g_src=np.array([rho_g]),
        mass_dst=np.array([mass]), eff_g_dst=np.array([eg_d]), eff_r_dst=np.array([er_d]),
    )
    assert bool(determined[0])
    # the destination's own observations do the conversion, exactly as the derivation states
    rho_r = (mass - rho_g * eg_d) / er_d
    assert out[0] == pytest.approx(np.log(rho_g / rho_r), abs=1e-12)
    # ...and the rebuilt densities still account for the destination's whole mass
    g, r = rebuild_densities(out[0], mass, eg_d, er_d)
    assert g * eg_d + r * er_d == pytest.approx(mass, rel=1e-10)
    assert g == pytest.approx(rho_g, rel=1e-10), "the shared density must cross UNCHANGED"


def test_the_shared_density_crosses_unchanged_over_a_WHOLE_CHAIN():
    """⛔ A14 — this gate exists because its absence was caught by a perturbation that fired NOTHING.

    Every other transfer fixture here is a SINGLE slot, so any defect expressed relative to the batch (a
    ratio of totals, a mean, a normalisation — i.e. exactly the retired reframe ``r``) is identically 1
    on them and hides completely. The sweep calls this vectorised over the whole chain, so the gate must
    be too, with slots whose masses genuinely differ.

    The invariant: the shared gDNA density arrives UNCHANGED at every destination, whatever that
    destination's own mass and opportunities are. Nothing in the transfer may consult a ratio of totals.
    """
    rng = np.random.default_rng(23)
    n = 512
    t_src = population_set(np.zeros(n, bool), np.zeros(n, bool), np.full(n, 400.0), np.zeros(n))
    t_dst = population_set(np.ones(n, bool), np.zeros(n, bool), np.full(n, 250.0), np.full(n, 100.0))
    rho_g = np.exp(rng.normal(-4.0, 1.0, n))
    eg_d = np.exp(rng.normal(5.5, 1.0, n))
    er_d = np.exp(rng.normal(5.0, 1.0, n))
    # masses spanning four decades, so a batch-relative rescale cannot be the identity
    mass = rho_g * eg_d * np.exp(rng.uniform(0.5, 9.0, n))

    out, determined = cross_composition(
        np.zeros(n), t_src, t_dst, rho_g_src=rho_g, mass_dst=mass,
        eff_g_dst=eg_d, eff_r_dst=er_d,
    )
    assert np.all(determined)
    g, r = rebuild_densities(out, mass, eg_d, er_d)
    assert np.max(np.abs(g - rho_g) / rho_g) < 1e-9, "a ratio of totals leaked into the transfer"
    assert np.max(np.abs(g * eg_d + r * er_d - mass) / mass) < 1e-10


def test_two_newly_active_strands_determine_the_total_and_not_the_tilt():
    """⭐⭐⭐ The ONLY thing a population mismatch ever leaves undetermined is the TILT, and only when both
    strands are newly active. The gDNA-vs-RNA split — the quantity the tool exists to estimate — is
    determined in every case."""
    t_src = population_set([False], [False], [400.0], [0.0])  # {g}
    t_dst = population_set([True], [True], [250.0], [100.0])  # {g, R+, R-}
    out, determined = cross_composition(
        np.array([0.0]), t_src, t_dst, rho_g_src=np.array([0.02]),
        mass_dst=np.array([1000.0]), eff_g_dst=np.array([250.0]), eff_r_dst=np.array([100.0]),
    )
    assert not bool(determined[0]), "both strands newly live => the tilt is not determined"
    assert np.isfinite(out[0]), "but the gDNA-vs-RNA split still is"


def test_an_inconsistent_shared_claim_resolves_to_an_absent_new_population_without_a_shift():
    """``A > M_d`` — the shared populations account for more than the destination observed. The null is
    inconsistent with the destination's own count, and the honest reading is that the new population is
    absent. ⛔ No shift is introduced: every rule that added doubt at pass-0 was priced on the full panel
    and refused by the zero control."""
    t_src = population_set([False], [False], [400.0], [0.0])
    t_dst = population_set([True], [False], [250.0], [100.0])
    out, _ = cross_composition(
        np.array([0.0]), t_src, t_dst, rho_g_src=np.array([99.0]),  # far more gDNA than M_d admits
        mass_dst=np.array([10.0]), eff_g_dst=np.array([250.0]), eff_r_dst=np.array([100.0]),
    )
    assert out[0] > 0.0 and np.isfinite(out[0]), "must saturate toward all-gDNA, not blow up"


def test_the_scalar_fast_path_agrees_with_the_vectorised_definition():
    """⛔ A11 — ONE HOME for the definition, and a gate on the fast path so the duplication cannot drift.

    The sweep's scan is a sequential loop and calls the scalar twin per hop. Over a randomised sweep that
    spans the saturating branch as well as the ordinary one, the two must agree exactly.
    """
    rng = np.random.default_rng(31)
    n = 4000
    rho_g = np.exp(rng.normal(-4.0, 2.5, n))
    eg = np.exp(rng.normal(5.5, 1.0, n))
    er = np.exp(rng.normal(5.0, 1.0, n))
    # a deliberately wide spread of masses, so ~a third of the draws trip the inconsistent-claim branch
    mass = rho_g * eg * np.exp(rng.normal(1.0, 3.0, n))

    t_src = population_set(np.zeros(n, bool), np.zeros(n, bool), np.full(n, 400.0), np.zeros(n))
    t_dst = population_set(np.ones(n, bool), np.zeros(n, bool), eg, er)
    vec, _det = cross_composition(
        np.zeros(n), t_src, t_dst, rho_g_src=rho_g, mass_dst=mass, eff_g_dst=eg, eff_r_dst=er
    )
    sca = np.array([residual_eta_scalar(rho_g[i], mass[i], eg[i], er[i]) for i in range(n)])

    saturated = mass - rho_g * eg <= 0.0
    assert saturated.sum() > n // 20, "fixture never reaches the saturating branch; it is untested"
    # off the saturating branch the two forms are the same arithmetic and must agree to machine precision
    assert np.max(np.abs(vec[~saturated] - sca[~saturated])) < 1e-9
    # on it, both must report all-gDNA rather than blowing up or changing sign
    assert np.all(sca[saturated] > 0.0) and np.all(np.isfinite(sca[saturated]))
    assert np.all(vec[saturated] > 0.0) and np.all(np.isfinite(vec[saturated]))


def test_the_mirror_gate_a_palindromic_fixture_delivers_mirror_image_answers():
    """⭐ A13's shape — an invariance the code cannot fake, rather than a reconstruction of a value the
    code may decide not to use.

    Swap the two opportunities and negate eta and every downstream quantity must mirror exactly: the
    shares exchange, the densities exchange, and the delivered eta negates. A gate that recomputed the
    transfer and compared would pass on a path where the composition never enters; this one cannot.
    """
    rng = np.random.default_rng(11)
    eta = rng.normal(0.0, 4.0, 256)
    mass = np.exp(rng.normal(5.0, 2.0, 256))
    eg = np.exp(rng.normal(6.0, 1.5, 256))
    er = np.exp(rng.normal(6.0, 1.5, 256))

    g, r = rebuild_densities(eta, mass, eg, er)
    g_m, r_m = rebuild_densities(-eta, mass, er, eg)  # the palindrome: swap roles, negate the coordinate
    assert np.max(np.abs(g_m - r)) < 1e-9 * np.max(np.abs(r))
    assert np.max(np.abs(r_m - g)) < 1e-9 * np.max(np.abs(g))

    phi_g, phi_r = shares_from_lambda(lambda_from_eta(eta, eg, er))
    phi_g_m, phi_r_m = shares_from_lambda(lambda_from_eta(-eta, er, eg))
    assert np.max(np.abs(phi_g_m - phi_r)) < 1e-12
    assert np.max(np.abs(phi_r_m - phi_g)) < 1e-12


def test_the_mirror_gate_holds_through_a_population_mismatch():
    """The mirror must survive the ``|N| = 1`` branch too, otherwise the residual rule has a handedness
    the derivation does not — the failure mode A13 exists to catch."""
    t_src = population_set([False], [False], [400.0], [0.0])
    t_dst_pos = population_set([True], [False], [250.0], [100.0])
    t_dst_neg = population_set([False], [True], [250.0], [100.0])
    kw = dict(
        rho_g_src=np.array([0.02]), mass_dst=np.array([1000.0]),
        eff_g_dst=np.array([250.0]), eff_r_dst=np.array([100.0]),
    )
    a, da = cross_composition(np.array([0.0]), t_src, t_dst_pos, **kw)
    b, db = cross_composition(np.array([0.0]), t_src, t_dst_neg, **kw)
    assert a[0] == pytest.approx(b[0], abs=1e-12), "the residual rule must not prefer a strand"
    assert bool(da[0]) and bool(db[0])


def test_eta_precision_crosses_unchanged_because_the_offset_is_a_constant():
    """eta and lambda differ by a CONSTANT, so a variance is invariant under the conversion. That is why
    a belief-free exact conversion has no scale noise to charge and ``transfer_logvar``'s non-graft branch
    has nothing left to price."""
    rng = np.random.default_rng(3)
    lam = rng.normal(0.0, 2.0, 1000)
    eg, er = 250.0, 100.0
    eta = eta_from_lambda(lam, eg, er)
    assert np.var(eta) == pytest.approx(np.var(lam), rel=1e-12)
    assert np.var(lambda_from_eta(eta, 900.0, 300.0)) == pytest.approx(np.var(lam), rel=1e-12)


def test_population_column_order_is_the_documented_one():
    """A13 again, in miniature: the column constants must not silently permute under a refactor."""
    assert (GDNA, RNA_POS, RNA_NEG) == (0, 1, 2)
    t = population_set([True], [False], [1.0], [1.0])
    assert t[0, GDNA] and t[0, RNA_POS] and not t[0, RNA_NEG]


# ── ψ's TILT COORDINATE, and the assertion that replaced a clip ─────────────────────────────────────
#
# ⛔⛔⛔ These gates exist because ONE unit error carried 74 % of the ``eta`` arm's error at the
# zero-gDNA control. The channel delivered a raw log-odds into ψ's ``theta`` slot, whose grid is the
# ANGLE ``arcsin(tau)`` and spans exactly ``[-pi/2, +pi/2]``. Nothing in the suite could see it: no gate
# read the tilt at all, and the arm's own debug capture did not publish it.


def test_the_tilt_is_an_angle_and_never_leaves_psis_grid():
    """⭐⭐⭐ THE GATE THE ``g00`` BLOW-UP NEEDED. ``theta`` must lie inside ψ's own tilt grid.

    ``simplex_logodds._tilt_grid`` spans ``[-pi/2, +pi/2]`` because ``theta = arcsin(tau)`` and
    ``tau`` is a tilt in ``[-1, 1]``. A mode outside that interval is not a claim ψ can represent: its
    Gaussian ``-1/2 p (theta - m)^2`` becomes MONOTONE across the whole grid, pinning the tilt at the
    boundary and destroying the nuisance integration that protects ``f_g`` at an AMBIG slot.

    ⛔ PERTURBATION: return ``log(rho_pos) - log(rho_neg)`` — the refuted form — and this fires on the
    extreme fixtures, where the log-odds reaches 6.9 against a bound of 1.5708.
    """
    rng = np.random.default_rng(20260806)
    rho_p = 10.0 ** rng.uniform(-3.0, 3.0, size=20_000)
    rho_n = 10.0 ** rng.uniform(-3.0, 3.0, size=20_000)
    prec = rng.uniform(0.5, 500.0, size=20_000)
    theta, _ = tilt_angle(rho_p, rho_n, prec, prec)
    assert np.all(np.isfinite(theta))
    assert np.abs(theta).max() <= np.pi / 2 + 1e-12, (
        f"a tilt of {np.abs(theta).max()} is outside ψ's grid; the mode is not representable"
    )
    # and the fixture must actually REACH the region a log-odds would have escaped from — A14: a gate
    # on a bounded fixture would pass with the refuted form in place.
    log_odds = np.log(rho_p) - np.log(rho_n)
    assert np.abs(log_odds).max() > np.pi / 2, "the fixture cannot distinguish the two forms"


def test_the_tilt_is_the_rna_tilt_and_carries_the_sign_of_the_density_difference():
    """⛔ The SOURCE, not only the coordinate. ``theta > 0`` iff there is more ``+``-strand RNA.

    Correcting the coordinate while still reading raw COUNTS was measured and did NOT help ``g00``,
    because F2 gives ``tau_obs = (2 kappa - 1)(1 - f_g) tau`` and an antisense protocol
    (``kappa = 0.0101``) flips the sign. Densities that are already ψ's own ``f_pos``/``f_neg`` up to a
    common factor carry no such factor.

    ⛔ PERTURBATION: swap the two density arguments and this fires on every row.
    """
    theta, _ = tilt_angle(
        np.array([9.0, 1.0, 5.0]), np.array([1.0, 9.0, 5.0]),
        np.array([100.0, 100.0, 100.0]), np.array([100.0, 100.0, 100.0]),
    )
    assert theta[0] > 0.0 and theta[1] < 0.0
    assert theta[0] == pytest.approx(-theta[1])
    assert theta[2] == pytest.approx(0.0)


def test_the_tilt_precision_is_the_delta_method_and_needs_no_constant():
    """The variance is ``((1 - tau^2)/4)(1/p_+ + 1/p_-)`` — checked against a numeric Jacobian.

    ⛔ PERTURBATION: drop the ``(1 - tau^2)`` factor, or the ``1/4``, and this fires.
    """
    rng = np.random.default_rng(11)
    rp = rng.uniform(0.5, 40.0, size=500)
    rn = rng.uniform(0.5, 40.0, size=500)
    pp = rng.uniform(2.0, 300.0, size=500)
    pn = rng.uniform(2.0, 300.0, size=500)
    _, prec = tilt_angle(rp, rn, pp, pn)

    h = 1e-6  # d(theta)/d(log rho_c) by central difference, one column at a time
    def theta_of(a, b):
        return tilt_angle(a, b, pp, pn)[0]

    d_p = (theta_of(rp * np.exp(h), rn) - theta_of(rp * np.exp(-h), rn)) / (2 * h)
    d_n = (theta_of(rp, rn * np.exp(h)) - theta_of(rp, rn * np.exp(-h))) / (2 * h)
    var_numeric = d_p**2 / pp + d_n**2 / pn
    np.testing.assert_allclose(1.0 / prec, var_numeric, rtol=1e-5)


def test_a_slot_with_rna_on_one_strand_only_states_no_tilt_rather_than_an_infinite_one():
    """⚠ The arcsin map is singular at ``tau = +-1``. The honest output there is precision 0.

    ⛔ PERTURBATION: drop the ``rho > 0`` requirement and the precision goes to ``inf``, which is a slot
    asserting a boundary tilt with unbounded confidence — the exact failure this whole family is about.
    """
    theta, prec = tilt_angle(
        np.array([5.0, 0.0, 5.0]), np.array([0.0, 5.0, 5.0]),
        np.array([50.0, 50.0, 50.0]), np.array([50.0, 50.0, 50.0]),
    )
    assert prec[0] == 0.0 and prec[1] == 0.0
    assert prec[2] > 0.0 and np.isfinite(prec[2])
    assert np.all(np.isfinite(theta))


# ── MONOTONICITY: more gDNA must never deliver less gDNA (owner, 2026-08-06) ─────────────────────────


def test_more_gdna_at_the_source_never_delivers_less_gdna_at_the_destination():
    """⭐⭐ THE OWNER'S MONOTONICITY GATE, at the one place a density becomes a composition.

    ``residual_eta_scalar`` is where a neighbour's gDNA density is turned into the destination's
    ``eta``. Raising that density must not lower the delivered ``eta``, hence must not lower ``f_g``:
    the mapping is monotone by derivation (a larger shared claim leaves a smaller RNA residual), and an
    inversion anywhere in the chain would show up here first. Cheap, and it catches a whole class of
    future sign and reciprocal errors.

    ⛔ PERTURBATION: swap the two ``log`` terms in ``residual_eta_scalar`` and this fires on every pair.
    """
    rng = np.random.default_rng(7)
    for _ in range(200):
        mass = float(rng.uniform(1.0, 1e5))
        eff_g = float(rng.uniform(1.0, 5e3))
        eff_r = float(rng.uniform(1.0, 5e3))
        rho = np.sort(10.0 ** rng.uniform(-6.0, np.log10(mass / eff_g), size=25))
        eta = np.array([residual_eta_scalar(r, mass, eff_g, eff_r) for r in rho])
        assert np.all(np.diff(eta) >= -1e-9), "a larger source gDNA density delivered a smaller eta"
        f_g = 1.0 / (1.0 + np.exp(-(eta + np.log(eff_g / eff_r))))
        assert np.all(np.diff(f_g) >= -1e-12), "a larger source gDNA density delivered a smaller f_g"


def test_monotonicity_holds_through_the_vectorised_transfer_too():
    """B14 — the scalar fast path and the vectorised definition are TWO places, so both are gated."""
    t_src = np.tile(np.array([True, False, False]), (30, 1))
    t_dst = np.tile(np.array([True, True, False]), (30, 1))
    rho = np.sort(10.0 ** np.linspace(-6.0, -1.0, 30))
    eta, determined = cross_composition(
        np.zeros(30), t_src, t_dst,
        rho_g_src=rho, mass_dst=np.full(30, 5_000.0),
        eff_g_dst=np.full(30, 800.0), eff_r_dst=np.full(30, 400.0),
    )
    assert determined.all()
    assert np.all(np.diff(eta) >= -1e-9)


# ── the ASSERTION that replaced the [0,1] clip ──────────────────────────────────────────────────────


def test_an_over_unit_share_raises_rather_than_being_clipped():
    """⭐⭐⭐ THE OWNER'S PRESCRIPTION, AND THE GATE THAT IT FIRES (`TRAPS.md` A2, A6).

    A share is in ``[0, 1]`` by definition. Clipping ``5.4e9`` to ``1.0`` states MAXIMAL gDNA — the worst
    available answer at a zero-gDNA library — turning a loud absurdity into a silent confident error. A
    prototype must be brittle.

    ⛔ This is the INJECTED out-of-range value the owner asked for: the assertion must FIRE, and the
    message must name the slot so the failure is actionable rather than a bare traceback.
    """
    rho = np.array([0.1, 0.1, 50.0])          # the third claims 50 frag/base
    eff = np.array([100.0, 100.0, 100.0])     # over 100 b => 5,000 fragments
    mass = np.array([1_000.0, 1_000.0, 10.0])  # at a slot holding 10
    with pytest.raises(AssertionError, match="share 500 > 1 at slot 2"):
        share_from_density(rho, eff, mass, np.array([1.0, 1.0, 1.0]), name="injected")


def test_a_muted_channel_is_not_asserted_on():
    """⚠ Only a claim that is CARRIED can be absurd. A precision of 0 is never read by ψ, and a slot
    with no mass would otherwise raise on a ratio nothing delivers.

    ⛔ PERTURBATION: drop the ``prec > 0`` conjunct and this fires — which would make the assertion
    unusable on any real chain and is how a brittle check gets weakened back into a clip.
    """
    out = share_from_density(
        np.array([50.0]), np.array([100.0]), np.array([10.0]), np.array([0.0]), name="muted"
    )
    assert np.isfinite(out).all()


def test_a_legal_share_is_the_log_of_exactly_the_share():
    """The happy path is not decoration: it is what says the assertion has not eaten the value."""
    out = share_from_density(
        np.array([0.25, 1.0]), np.array([100.0, 100.0]), np.array([1_000.0, 100.0]),
        np.array([3.0, 3.0]), name="legal",
    )
    np.testing.assert_allclose(np.exp(out), [0.025, 1.0], rtol=1e-12)
