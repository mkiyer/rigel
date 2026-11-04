"""The gDNA background rate: the closed form `poisson_lower_mean` and the contamination-robust one-sided root."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import poisson

from rigel.calibration.gdna_density import (
    one_sided_rate,
    poisson_lower_mean,
)


LAMS = (0.05, 0.3, 1.0, 2.7, 10.0, 55.5, 200.0, 999.0)


# ── the closed form ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("lam", LAMS)
def test_poisson_lower_mean_matches_the_exact_truncated_sum(lam):
    """``E[(lam-N)+] = sum_{k<=lam} (lam-k) P(N=k)``: a different algorithm, not the same formula twice."""
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
    """``lam**(lam+1)/Gamma(lam+1)``, with no ``floor``, is a different function."""
    from scipy.special import gammaln

    lam = 2.7
    wrong = float(np.exp((lam + 1.0) * np.log(lam) - lam - gammaln(lam + 1.0)))
    assert wrong != pytest.approx(poisson_lower_mean(np.array([lam]))[0], rel=1e-6)


def test_the_log_gamma_table_is_the_direct_form_TO_THE_BIT():
    """The log-gamma inside the Poisson identity is evaluated at integers, so it is tabulated, and a table of the same
    function at the same points is exact: over a dense grid the two forms must be equal, not close. Both branches are
    covered (the table is used only when it is smaller than the data, so short arrays fall back to the direct call and
    long ones take the table); a table built at ``k`` instead of ``k + 1`` would fail on every row.
    """
    from scipy.special import gammaln

    from rigel.calibration.gdna_density import _log_factorial

    for size in (1, 5, 1000, 200000):
        rng = np.random.default_rng(size)
        for hi in (0.5, 3.0, 40.0, 5000.0):
            k = np.floor(rng.uniform(0.0, hi, size=size))
            assert np.array_equal(_log_factorial(k), gammaln(k + 1.0)), f"size={size} hi={hi}"
    # the degenerate ends: an empty array, every k zero, and one huge k that must take the fallback
    assert np.array_equal(_log_factorial(np.zeros(0)), gammaln(np.zeros(0) + 1.0))
    assert np.array_equal(_log_factorial(np.zeros(7)), gammaln(np.zeros(7) + 1.0))
    big = np.array([1e9, 2.0, 3.0])
    assert np.array_equal(_log_factorial(big), gammaln(big + 1.0))


# ── the one-sided rate ───────────────────────────────────────────────────────────────────────────


def test_uncontaminated_poisson_recovers_its_own_rate():
    rng = np.random.default_rng(11)
    e = rng.uniform(50.0, 5000.0, size=4000)
    rho = 0.02
    fit = one_sided_rate(rng.poisson(rho * e).astype(float), e)
    assert fit.bracket_ok and fit.informative
    assert fit.rate == pytest.approx(rho, rel=0.05)
    assert fit.rate_over_pooled == pytest.approx(1.0, rel=0.05)


def test_the_bisection_on_its_closed_bracket_is_the_fixed_200_step_loop_TO_THE_BIT():
    """The bisection stops when its float64 bracket closes, about sixty halvings, and must equal a fixed 200-step loop
    to the bit, since a closed bracket can no longer change the value. Stopping one halving early (``mid <= lo`` relaxed
    to a relative width) would move the rate.
    """
    from rigel.calibration.gdna_density import poisson_lower_mean as _plm

    def fixed_200(n, e):
        total_e, total_n = float(e.sum()), float(n.sum())
        pooled = total_n / total_e

        def f(rho):
            lam = rho * e
            return float(np.clip(lam - n, 0.0, None).sum() - _plm(lam).sum())

        lo, hi = 0.0, 10.0 * pooled
        assert f(lo) <= 0.0 <= f(hi)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0.0:
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    rng = np.random.default_rng(5)
    for scale, size in ((0.002, 500), (0.02, 4000), (2.0, 300), (1e-5, 1000)):
        e = rng.uniform(50.0, 5000.0, size=size)
        n = rng.poisson(scale * e).astype(float)
        if n.sum() == 0.0:
            continue
        got = one_sided_rate(n, e)
        assert got.bracket_ok
        assert got.rate == fixed_200(n, e), f"the closed bracket moved the rate at rho={scale}"


def test_contamination_moves_the_pooled_rate_but_not_the_one_sided_one():
    """Same gDNA, half the objects given a large additive contaminant: the pooled rate moves, the one-sided one does
    not."""
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
    """Contamination pushes the root up, never down, as the derivation claims."""
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
    """A zero-gDNA library has no gDNA density and must say so rather than return a number."""
    e = np.full(100, 500.0)
    fit = one_sided_rate(np.zeros(100), e)
    assert not fit.informative and fit.rate == 0.0


def test_no_exposure_declines():
    assert not one_sided_rate(np.ones(10), np.zeros(10)).informative


def test_the_returned_rate_is_actually_a_root():
    """``F(rate) = 0``, re-derived here rather than trusting the bisection. A check of ``rate > 0`` alone could not
    fail, because a bracket starting at exactly 0.0 already walks off the boundary (``F(0+) < 0``), so the root is
    checked directly.
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
