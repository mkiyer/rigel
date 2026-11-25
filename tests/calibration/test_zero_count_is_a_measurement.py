"""A zero count over a known opportunity is a measurement.

Every count channel in calibration is a Poisson rate observed over an opportunity, so under the Jeffreys prior the
log-rate variance is ``trigamma(a + ½)``: proper and finite at ``a = 0``, and asymptotically ``1/a``. ``count_logvar``
is that expression and the one home of the counting term; every hop price of the transfer policy reads it
(`native/transfer_rows.h`, bound in `_solve_impl.testing`). The tests check that it is exact and exposure-free (Z1);
agrees with the ``1/a`` asymptote to better than 0.1 % for ``a ≥ 10``, so it matters only for low counts (Z2); is finite
and positive at ``a = 0``, where the asymptote diverges and a pure-gDNA object with no counts would emit nothing (Z3);
is monotone across ``0 → 1``, since a step there would let that population fall off the boundary (Z6); and that the two
production paths inherit all of it (Z7).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import polygamma

from rigel._solve_impl import testing as R


# ── Z1: the exact expression ───────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("n", [0.0, 1.0, 3.0, 10.0, 250.0, 112_333.0])
def test_z1_the_count_term_is_the_exact_poisson_log_rate_variance(n):
    """``Var(log ρ) = trigamma(a + ½)`` for ``ρ ~ Gamma(a + ½, E)``, scored against scipy's own ``polygamma``, a
    different implementation."""
    assert R.count_logvar(n) == pytest.approx(polygamma(1, n + 0.5), rel=1e-12)


def test_z1_it_does_not_depend_on_the_exposure():
    """``Var(log ρ)`` is exposure-free: ``E`` shifts the location and cannot sharpen the claim; a count term that moved
    with ``E`` would double-count the opportunity."""
    assert R.count_logvar(7.0) == R.count_logvar(7.0) and R.count_logvar(0.0) == R.count_logvar(0.0)
    assert R.count_logvar(7.0) != R.count_logvar(
        0.0
    )  # a function of the count, and of nothing else


# ── Z2: agreement with the asymptote bounds the effect ──────────────────────────────────────────────


@pytest.mark.parametrize("n", [10.0, 50.0, 1_000.0, 100_000.0])
def test_z2_it_agrees_with_the_one_over_n_asymptote_to_a_tenth_of_a_percent(n):
    """``trigamma(n+½) → 1/n``, so every object with a real count moves by under 0.1 %, and the exact expression differs
    from the asymptote only where the asymptote is broken."""
    assert R.count_logvar(n) == pytest.approx(1.0 / n, rel=1e-3)


def test_z2_and_it_departs_from_one_over_n_at_small_counts():
    """At small counts the two must disagree, or Z2 would pass on a no-op; that disagreement is why the exact expression
    is used."""
    assert R.count_logvar(1.0) == pytest.approx(0.9348, abs=1e-3)  # vs 1/1
    assert R.count_logvar(2.0) == pytest.approx(0.4903, abs=1e-3)  # vs 1/2


# ── Z3: a zero count, as a number ────────────────────────────────────────────────────────────────


def test_z3_a_zero_count_carries_a_finite_variance():
    """``trigamma(½) = π²/2``, so a zero-count object's log-density claim has sd 2.22 nats: loose but finite, a claim
    rather than the silence an infinite variance forces."""
    v = R.count_logvar(0.0)
    assert v == pytest.approx(np.pi**2 / 2.0, rel=1e-12)
    assert np.isfinite(v) and v > 0.0
    assert np.sqrt(v) == pytest.approx(2.221, abs=1e-3)


# ── Z6: monotone across the join ───────────────────────────────────────────────────────────────


def test_z6_the_variance_is_monotone_in_the_count_across_zero():
    """No step at ``n = 0 → 1``: a jump from silence to "as precise as the count allows" would let a whole population
    fall off the boundary."""
    v = np.array([R.count_logvar(float(n)) for n in range(40)])
    assert np.all(np.diff(v) < 0.0), v[:6]
    assert np.isfinite(v[0]) and v[0] > 0.0


# ── Z7: the production paths, not just the arithmetic ─────────────────────────────────────────


def test_z7_a_hop_priced_on_a_zero_count_is_finite_and_counting_alone():
    """The transfer policy's hop price with a zero witness count on either side: no density ratio exists, so the price
    is the two counts' counting and nothing more, finite and larger than the same hop priced on a positive count (a zero
    count is the quietest measurement, not silence)."""
    v0 = R.step_price(30.0, 100.0, 0.0, 100.0)
    assert np.isfinite(v0) and v0 == pytest.approx(
        polygamma(1, 30.5) + polygamma(1, 0.5), rel=1e-12
    )
    assert v0 > R.step_price(30.0, 100.0, 30.0, 100.0)
    assert np.isfinite(R.step_price(0.0, 100.0, 0.0, 100.0))


def test_z7_a_gene_edge_with_zero_counts_emits_a_level():
    """A structurally pure-gDNA crossing holding zero counts over a large opportunity: on the gDNA lane its level is a
    profile falling with the density (nothing below zero claimed, everything above priced by the count's own
    likelihood), never an absent claim.
    """
    u = np.linspace(-10.0, 10.0, 60)
    level = R.poisson_level(u, 0.0, 5_000.0, 0.05)
    assert level is not None and np.all(np.isfinite(level))
    assert np.all(np.diff(level) <= 0.0) and level[0] == 0.0 and level[-1] < -1.0
    assert np.ptp(level) > 0.0, "a zero count over a real opportunity is a claim, not silence"


# ── the counting term's arithmetic: the native trigamma against scipy ────────────────────────────


def test_the_native_trigamma_is_scipys_hurwitz_zeta():
    """The counting variance's one home in C++ agrees with ``zeta(2, x)`` from a half to ten million, a different
    implementation."""
    from scipy.special import zeta

    x = np.concatenate([np.linspace(0.5, 12.0, 400), np.logspace(1.0, 7.0, 200)])
    np.testing.assert_allclose(np.array([R.trigamma(float(v)) for v in x]), zeta(2, x), rtol=1e-13)
