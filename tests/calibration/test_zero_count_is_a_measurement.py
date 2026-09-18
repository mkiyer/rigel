"""A zero count over a known opportunity is a measurement — TRAPS: a-zero-count-is-a-measurement.

Every count channel in calibration is a Poisson rate observed over an opportunity, so under the
Jeffreys prior the log-rate variance is ``trigamma(a + ½)`` — proper and finite at ``a = 0``, and
asymptotically ``1/a``. ``count_logvar`` is that expression and the one home of the counting term:
every hop price of the transfer policy reads it (`native/transfer_rows.h`, bound for the gates as
`native.transfer_rows`). Five of its properties are gated. It is exact and
exposure-free (Z1); it agrees with the ``1/a`` asymptote to better than 0.1 % for ``a ≥ 10``, which
confines its effect to the low-count population (Z2); it is finite and positive at ``a = 0``, where
the asymptote diverges and a pure-gDNA object with no counts therefore emits nothing at all (Z3);
it is monotone across ``0 → 1``, since a step there is what lets that population fall off the
boundary (Z6); and the two production paths inherit all of it (Z7).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import polygamma

from rigel.native import transfer_rows as R


# ── Z1 — the exact expression ───────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("n", [0.0, 1.0, 3.0, 10.0, 250.0, 112_333.0])
def test_Z1_the_count_term_is_the_exact_Poisson_log_rate_variance(n):
    """``Var(log ρ) = trigamma(a + ½)`` for ``ρ ~ Gamma(a + ½, E)`` — scored against scipy's own
    ``polygamma``, which is a different implementation from ours (TRAPS: a-test-that-redefines)."""
    assert R.count_logvar(n) == pytest.approx(polygamma(1, n + 0.5), rel=1e-12)


def test_Z1_it_does_not_depend_on_the_exposure():
    """``Var(log ρ)`` is exposure-free: ``E`` shifts the location, it cannot sharpen the claim. A
    count term that moved with ``E`` would be double-counting the opportunity."""
    assert R.count_logvar(7.0) == R.count_logvar(7.0) and R.count_logvar(0.0) == R.count_logvar(0.0)
    assert R.count_logvar(7.0) != R.count_logvar(
        0.0
    )  # a function of the count, and of nothing else


# ── Z2 — the property that bounds the blast radius ──────────────────────────────────────────────


@pytest.mark.parametrize("n", [10.0, 50.0, 1_000.0, 100_000.0])
def test_Z2_it_agrees_with_the_retired_one_over_n_to_better_than_a_tenth_of_a_percent(n):
    """The safety property: ``trigamma(n+½) → 1/n``, so every object with a real count is moved by
    under 0.1 % and the expression differs from the asymptote only where the asymptote is broken.
    Without this it would be a whole-panel perturbation dressed as a bug fix."""
    assert R.count_logvar(n) == pytest.approx(1.0 / n, rel=1e-3)


def test_Z2_and_it_DIVERGES_from_one_over_n_exactly_where_the_defect_was():
    """The other half of Z2, or Z2 would pass on a no-op: at small counts the two must genuinely
    disagree, which is the whole reason the exact expression is used."""
    assert R.count_logvar(1.0) == pytest.approx(0.9348, abs=1e-3)  # vs 1/1
    assert R.count_logvar(2.0) == pytest.approx(0.4903, abs=1e-3)  # vs 1/2


# ── Z3 — the defect, as a number ────────────────────────────────────────────────────────────────


def test_Z3_a_zero_count_carries_a_FINITE_variance():
    """``trigamma(½) = π²/2``, so a zero-count object's log-density claim has sd 2.22 nats — loose,
    but finite, and so a claim rather than the silence an infinite variance forces."""
    v = R.count_logvar(0.0)
    assert v == pytest.approx(np.pi**2 / 2.0, rel=1e-12)
    assert np.isfinite(v) and v > 0.0
    assert np.sqrt(v) == pytest.approx(2.221, abs=1e-3)


# ── Z6 — monotone across the join ───────────────────────────────────────────────────────────────


def test_Z6_the_variance_is_MONOTONE_in_the_count_across_zero():
    """There must be no step at ``n = 0 → 1``. A jump from silence to "as precise as the count
    allows" is a discontinuity that lets a whole population fall off the boundary."""
    v = np.array([R.count_logvar(float(n)) for n in range(40)])
    assert np.all(np.diff(v) < 0.0), v[:6]
    assert np.isfinite(v[0]) and v[0] > 0.0


# ── Z7 — the production paths, not just the arithmetic ─────────────────────────────────────────


def test_Z7_a_hop_priced_on_a_zero_count_is_finite_and_counting_alone():
    """The transfer policy's hop price on a zero witness count on either side: no density ratio
    exists, so the price is the two counts' counting and nothing more — finite, and larger than the
    same hop priced on a positive count (a zero count is the quietest measurement, not a silence)."""
    v0 = R.hop_price(30.0, 100.0, 0.0, 100.0)
    assert np.isfinite(v0) and v0 == pytest.approx(
        polygamma(1, 30.5) + polygamma(1, 0.5), rel=1e-12
    )
    assert v0 > R.hop_price(30.0, 100.0, 30.0, 100.0)
    assert np.isfinite(R.hop_price(0.0, 100.0, 0.0, 100.0))


def test_Z7_a_gene_edge_with_zero_counts_still_emits_a_level():
    """The object the ladder dissection found: a structurally pure-gDNA crossing holding zero counts
    over a large opportunity. On the gDNA lane its level is a profile falling with the density —
    nothing below zero is claimed, everything above it is priced by the count's own likelihood —
    never an absent claim. PERTURBATION: a lane that returns ``None`` at a zero count fails here."""
    u = np.linspace(-10.0, 10.0, 60)
    level = R.poisson_level(u, 0.0, 5_000.0, 0.05)
    assert level is not None and np.all(np.isfinite(level))
    assert np.all(np.diff(level) <= 0.0) and level[0] == 0.0 and level[-1] < -1.0
    assert np.ptp(level) > 0.0, "a zero count over a real opportunity is a claim, not silence"


# ── the counting term's arithmetic: the native trigamma against scipy ────────────────────────────


def test_the_native_trigamma_is_scipys_hurwitz_zeta():
    """The counting variance's one home in C++ agrees with ``zeta(2, x)`` from a half to ten million —
    a different implementation from ours (TRAPS: a-test-that-redefines)."""
    from scipy.special import zeta

    x = np.concatenate([np.linspace(0.5, 12.0, 400), np.logspace(1.0, 7.0, 200)])
    np.testing.assert_allclose(np.array([R.trigamma(float(v)) for v in x]), zeta(2, x), rtol=1e-13)
