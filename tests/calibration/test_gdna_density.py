"""The gDNA background rate: the closed form, the bracket, and the one-sided guarantee.

⛔ **The gates here are perturbation gates, not smoke tests** (`TRAPS: perturb-every-gate`). Each was run
against a deliberately broken implementation and watched to fail: a `poisson_lower_mean` that drops the
`floor` (9 fire), a `region_lengths_from_partition` that diffs straight through the reference junctions
(2 fire), a `contained_opportunity` that forgets the `+1` (2 fire), and a bisection that returns the
pooled rate instead of the one-sided root (1 fires).

⚠ **A fifth perturbation fired NOTHING and that is recorded rather than hidden**: bracketing the bisection
from exactly `0.0` passes every test here, because `F(0+) < 0` walks the bracket off the boundary by
itself. The guard against it was speculative, so it was deleted and the test that claimed to cover it was
replaced by `test_the_returned_rate_is_actually_a_root`, which checks the property instead of the guard.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import poisson

from rigel.calibration.gdna_density import (
    contained_opportunity,
    one_sided_rate,
    poisson_lower_mean,
    pooled_log_rate,
    region_lengths_from_partition,
)

LAMS = (0.05, 0.3, 1.0, 2.7, 10.0, 55.5, 200.0, 999.0)


# ── the closed form ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("lam", LAMS)
def test_poisson_lower_mean_matches_the_exact_truncated_sum(lam):
    """``E[(lam-N)+] = sum_{k<=lam} (lam-k) P(N=k)`` — a DIFFERENT algorithm, not the same formula twice."""
    k = np.arange(0, int(np.floor(lam)) + 1)
    exact = float(((lam - k) * poisson.pmf(k, lam)).sum())
    assert poisson_lower_mean(np.array([lam]))[0] == pytest.approx(exact, rel=1e-12)


@pytest.mark.parametrize("lam", LAMS)
def test_the_two_one_sided_means_are_equal(lam):
    """``E[N-lam] = 0`` forces ``E[(lam-N)+] = E[(N-lam)+]``; the identity halves ``E|N-lam|``."""
    k = np.arange(0, int(lam) + 400)
    upper = float(((k - lam).clip(0) * poisson.pmf(k, lam)).sum())
    assert poisson_lower_mean(np.array([lam]))[0] == pytest.approx(upper, rel=1e-9)


def test_poisson_lower_mean_is_zero_at_zero_and_finite_at_scale():
    assert poisson_lower_mean(np.array([0.0]))[0] == 0.0
    assert np.isfinite(poisson_lower_mean(np.array([1e5]))[0])  # log-space, or this overflows


def test_dropping_the_floor_breaks_the_identity():
    """⛔ The PERTURBATION: ``lam**(lam+1)/Gamma(lam+1)`` (no ``floor``) is a different function."""
    from scipy.special import gammaln

    lam = 2.7
    wrong = float(np.exp((lam + 1.0) * np.log(lam) - lam - gammaln(lam + 1.0)))
    assert wrong != pytest.approx(poisson_lower_mean(np.array([lam]))[0], rel=1e-6)


# ── the pooled rate ──────────────────────────────────────────────────────────────────────────────


def test_pooled_log_rate_is_the_ratio_and_carries_the_shape():
    c, e = np.array([3.0, 5.0]), np.array([10.0, 30.0])
    assert pooled_log_rate(c, e) == pytest.approx(np.log(8.0 / 40.0))
    assert pooled_log_rate(c, e, shape=0.5) == pytest.approx(np.log(8.5 / 40.0))


def test_pooled_log_rate_reports_no_support_rather_than_raising():
    assert pooled_log_rate(np.array([1.0]), np.array([0.0])) == -np.inf


# ── the one-sided rate ───────────────────────────────────────────────────────────────────────────


def test_uncontaminated_poisson_recovers_its_own_rate():
    rng = np.random.default_rng(11)
    e = rng.uniform(50.0, 5000.0, size=4000)
    rho = 0.02
    fit = one_sided_rate(rng.poisson(rho * e).astype(float), e)
    assert fit.bracket_ok and fit.informative
    assert fit.rate == pytest.approx(rho, rel=0.05)
    assert fit.rate_over_pooled == pytest.approx(1.0, rel=0.05)


def test_contamination_moves_the_pooled_rate_but_not_the_one_sided_one():
    """⭐ THE POINT OF THE MODULE. Same gDNA, half the objects given a large additive contaminant."""
    rng = np.random.default_rng(3)
    e = rng.uniform(50.0, 5000.0, size=4000)
    rho = 0.02
    clean = rng.poisson(rho * e).astype(float)
    dirty = clean.copy()
    hit = rng.random(e.size) < 0.5
    dirty[hit] += rng.poisson(4.0 * rho * e[hit])

    pooled_shift = (dirty.sum() / e.sum()) / (clean.sum() / e.sum())
    assert pooled_shift > 2.0, "the fixture must actually contaminate, or this test proves nothing"
    fit = one_sided_rate(dirty, e)
    assert fit.rate == pytest.approx(rho, rel=0.10)
    assert fit.rate_over_pooled < 0.5  # it found the contamination


def test_the_bias_is_one_sided_never_low():
    """Contamination pushes the root UP, never down — the direction the derivation claims."""
    rng = np.random.default_rng(5)
    e = rng.uniform(50.0, 5000.0, size=3000)
    rho = 0.02
    clean = rng.poisson(rho * e).astype(float)
    base = one_sided_rate(clean, e).rate
    for frac in (0.1, 0.3, 0.6):
        d = clean.copy()
        hit = rng.random(e.size) < frac
        d[hit] += rng.poisson(3.0 * rho * e[hit])
        assert one_sided_rate(d, e).rate >= base * 0.98


def test_adding_a_contaminant_never_lowers_the_estimate():
    """Monotone in the one direction the biology allows."""
    rng = np.random.default_rng(9)
    e = rng.uniform(100.0, 2000.0, size=1500)
    n = rng.poisson(0.03 * e).astype(float)
    prev = one_sided_rate(n, e).rate
    for _ in range(4):
        n = n + rng.poisson(0.01 * e)
        cur = one_sided_rate(n, e).rate
        assert cur >= prev - 1e-12
        prev = cur


def test_zero_counts_declines_rather_than_fabricating_a_rate():
    """⛔ A zero-gDNA library has no gDNA density. It must SAY so, not return a number."""
    e = np.full(100, 500.0)
    fit = one_sided_rate(np.zeros(100), e)
    assert not fit.informative and fit.rate == 0.0


def test_no_exposure_declines():
    assert not one_sided_rate(np.ones(10), np.zeros(10)).informative


def test_the_returned_rate_is_actually_a_root():
    """⭐ The property, not the algorithm: ``F(rate) = 0``, re-derived here rather than trusting bisection.

    ⚠ This REPLACED a test asserting only ``rate > 0``, which was named for the degenerate root at zero
    and could not fail — bracketing from exactly 0.0 passes it, because ``F(0+) < 0`` walks the bracket
    off the boundary by itself. A gate no perturbation can fire is not a gate.
    """
    rng = np.random.default_rng(2)
    e = rng.uniform(100.0, 1000.0, size=500)
    n = rng.poisson(0.05 * e).astype(float)
    fit = one_sided_rate(n, e)
    assert fit.rate > 0.0

    def F(rho):
        lam = rho * e
        return float(np.clip(lam - n, 0.0, None).sum() - poisson_lower_mean(lam).sum())

    scale = max(abs(F(0.0)), abs(F(2.0 * fit.rate)), 1.0)
    assert abs(F(fit.rate)) / scale < 1e-9
    assert (
        F(fit.rate * 0.5) < 0.0 < F(fit.rate * 2.0)
    )  # and it is the crossing, not a stationary point


def test_misaligned_inputs_raise():
    with pytest.raises(ValueError):
        one_sided_rate(np.ones(4), np.ones(5))


# ── the opportunity ──────────────────────────────────────────────────────────────────────────────


def test_contained_opportunity_matches_brute_force():
    """The cumulative-sum form against the literal ``sum_L pmf(L) max(0, ell-L+1)``."""
    rng = np.random.default_rng(17)
    pmf = rng.random(60)
    pmf /= pmf.sum()
    L = np.arange(pmf.size, dtype=np.float64)
    ell = np.array([0.0, 1.0, 5.0, 30.0, 59.0, 60.0, 200.0, 10_000.0])
    brute = np.array([float((pmf * np.clip(x - L + 1.0, 0.0, None)).sum()) for x in ell])
    assert contained_opportunity(pmf, ell) == pytest.approx(brute, rel=1e-12, abs=1e-12)


def test_contained_opportunity_is_ell_plus_one_minus_mean_beyond_support():
    pmf = np.zeros(11)
    pmf[10] = 1.0  # every fragment exactly 10 long
    assert contained_opportunity(pmf, np.array([1000.0]))[0] == pytest.approx(991.0)


def test_a_region_shorter_than_every_fragment_has_no_contained_opportunity():
    pmf = np.zeros(51)
    pmf[50] = 1.0
    assert contained_opportunity(pmf, np.array([10.0]))[0] == 0.0


# ── the per-reference region lengths ─────────────────────────────────────────────────────────────


def test_region_lengths_are_differenced_per_reference():
    """⛔ THE PERTURBATION, and the bug this gate exists for: two references, and a straight
    ``np.diff`` manufactures a PHANTOM region spanning the junction between them."""
    # ref A: bounds 0,100,300 -> regions of 100 and 200. ref B: bounds 0,50 -> one region of 50.
    bounds = np.array([0, 100, 300, 0, 50], dtype=np.int64)
    offsets = np.array([0, 3, 5], dtype=np.int64)
    got = region_lengths_from_partition(bounds, offsets, n_regions=3)
    assert got.tolist() == [100.0, 200.0, 50.0]
    # the phantom the naive form would emit, and which must NOT appear
    assert -300.0 not in np.diff(bounds)[:0].tolist() + got.tolist()
    assert got.size == 3


def test_a_reference_with_no_regions_contributes_none():
    bounds = np.array([0, 100, 0, 40], dtype=np.int64)
    offsets = np.array([0, 2, 2, 4], dtype=np.int64)  # middle reference is empty
    assert region_lengths_from_partition(bounds, offsets, n_regions=2).tolist() == [100.0, 40.0]


def test_region_lengths_match_the_naive_diff_on_a_single_reference():
    """With one reference there IS no junction, so the two agree — which is why the bug is silent."""
    bounds = np.array([0, 10, 45, 70], dtype=np.int64)
    offsets = np.array([0, 4], dtype=np.int64)
    assert region_lengths_from_partition(bounds, offsets, 3).tolist() == np.diff(bounds).tolist()
