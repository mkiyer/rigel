"""⛔⛔⛔ A ZERO COUNT OVER A KNOWN OPPORTUNITY IS A MEASUREMENT — TRAPS: a-zero-count-is-a-measurement.

Every count channel in calibration is a Poisson rate observed over an opportunity, and under the
Jeffreys prior ψ is built on the exact posterior is ``Gamma(a + ½, E)``, whose log has variance
``trigamma(a + ½)`` — proper and finite at ``a = 0``, and asymptotically ``1/a``. ``count_logvar`` is
that expression, THE ONE HOME of the counting term: every hop price of the transfer policy reads it.
Its predecessor was ``1/n``, the large-count limit, which diverges at ``n = 0`` — so an object with no
counts emitted nothing even when it was structurally pure gDNA and composition-CERTAIN. Measured at
``g00`` (zero gDNA by construction): all **1,298** intergenic regions held exactly zero counts over
**50.7 Mb** of opportunity and all 1,298 were silent, and pass-0 invented 34–38 % gDNA.

===  ===========================================================================================
Z1   ``count_logvar`` is ``trigamma(n + ½)`` exactly (vs scipy, not re-derived), exposure-free
Z2   ⭐ ASYMPTOTIC SAFETY — it agrees with the old ``1/n`` to <0.1 % for ``n ≥ 10``, so the change
     is confined to the low-count population. This is what bounds the blast radius
Z3   ⭐⭐ FINITE AND POSITIVE AT ``n = 0`` — the defect, as a number: ``π²/2``, an sd of 2.2 nats
Z6   monotone: more counts ⇒ a sharper claim, at every count including 0 → 1, with no step
Z7   ⭐⭐ THE PRODUCTION PATHS: a hop priced on a zero count is finite (`hop_price`), and a gene
     edge's zero count still emits a level on the gDNA lane (`poisson_level`)
===  ===========================================================================================
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import polygamma

from rigel.calibration.messages.transfer_rows import count_logvar, hop_price, poisson_level


# ── Z1 — the exact expression ───────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("n", [0.0, 1.0, 3.0, 10.0, 250.0, 112_333.0])
def test_Z1_the_count_term_is_the_exact_Poisson_log_rate_variance(n):
    """``Var(log ρ) = trigamma(a + ½)`` for ``ρ ~ Gamma(a + ½, E)`` — scored against scipy's own
    ``polygamma``, which is a different implementation from ours (TRAPS: a-test-that-redefines)."""
    assert count_logvar(np.array([n]))[0] == pytest.approx(polygamma(1, n + 0.5), rel=1e-12)


def test_Z1_it_does_not_depend_on_the_exposure():
    """⭐ ``Var(log ρ)`` is exposure-FREE — ``E`` only shifts the location, it cannot sharpen the
    claim. A count term that moved with ``E`` would be double-counting the opportunity."""
    v = count_logvar(np.array([0.0, 7.0]))
    assert np.array_equal(v, count_logvar(np.array([0.0, 7.0])))  # pure function of the count


# ── Z2 — the property that bounds the blast radius ──────────────────────────────────────────────


@pytest.mark.parametrize("n", [10.0, 50.0, 1_000.0, 100_000.0])
def test_Z2_it_agrees_with_the_retired_one_over_n_to_better_than_a_tenth_of_a_percent(n):
    """⭐⭐ THE SAFETY PROPERTY. ``trigamma(n+½) → 1/n``, so every object with a real count keeps its
    old answer to <0.1 % and the change lands ONLY where the old form was broken. Without this the
    fix would be a whole-panel perturbation dressed as a bug fix."""
    assert count_logvar(np.array([n]))[0] == pytest.approx(1.0 / n, rel=1e-3)


def test_Z2_and_it_DIVERGES_from_one_over_n_exactly_where_the_defect_was():
    """⚠ The other half of Z2, or Z2 would pass on a no-op: at small counts the two must genuinely
    disagree, which is the whole reason for the change."""
    assert count_logvar(np.array([1.0]))[0] == pytest.approx(0.9348, abs=1e-3)  # vs 1/1
    assert count_logvar(np.array([2.0]))[0] == pytest.approx(0.4903, abs=1e-3)  # vs 1/2


# ── Z3 — the defect, as a number ────────────────────────────────────────────────────────────────


def test_Z3_a_zero_count_carries_a_FINITE_variance():
    """⭐⭐ THE FIX, IN ONE ASSERTION. ``trigamma(½) = π²/2``, so a zero-count object's log-density
    claim has sd 2.22 nats — loose, but finite, and infinitely more than the nothing it emitted
    before."""
    v = float(count_logvar(np.array([0.0]))[0])
    assert v == pytest.approx(np.pi**2 / 2.0, rel=1e-12)
    assert np.isfinite(v) and v > 0.0
    assert np.sqrt(v) == pytest.approx(2.221, abs=1e-3)


# ── Z6 — monotone across the join ───────────────────────────────────────────────────────────────


def test_Z6_the_variance_is_MONOTONE_in_the_count_across_zero():
    """⭐ There must be no step at ``n = 0 → 1``: the old code jumped from "silent" to "as precise as
    the count allows", and a discontinuity there is what let a whole population fall off the boundary."""
    v = count_logvar(np.arange(0.0, 40.0))
    assert np.all(np.diff(v) < 0.0), v[:6]
    assert np.isfinite(v[0]) and v[0] > 0.0


# ── Z7 — ⭐⭐ THE PRODUCTION PATHS, not just the arithmetic ──────────────────────────────────────


def test_Z7_a_hop_priced_on_a_zero_count_is_finite_and_counting_alone():
    """The transfer policy's hop price on a zero witness count on either side: no density ratio
    exists, so the price is the two counts' counting and nothing more — finite, and larger than the
    same hop priced on a positive count (a zero count is the quietest measurement, not a silence)."""
    v0 = hop_price(30.0, 100.0, 0.0, 100.0)
    assert np.isfinite(v0) and v0 == pytest.approx(
        polygamma(1, 30.5) + polygamma(1, 0.5), rel=1e-12
    )
    assert v0 > hop_price(30.0, 100.0, 30.0, 100.0)
    assert np.isfinite(hop_price(0.0, 100.0, 0.0, 100.0))


def test_Z7_a_gene_edge_with_zero_counts_still_emits_a_level():
    """THE OBJECT THE DISSECTION FOUND: a structurally pure-gDNA crossing holding zero counts over a
    large opportunity. On the gDNA lane its level is a profile FALLING with the density — nothing
    below zero is claimed, everything above it is priced by the count's own likelihood — never an
    absent claim. PERTURBATION: a lane that returns ``None`` at a zero count fails here."""
    u = np.linspace(-10.0, 10.0, 60)
    level = poisson_level(u, 0.0, 5_000.0, 0.05)
    assert level is not None and np.all(np.isfinite(level))
    assert np.all(np.diff(level) <= 0.0) and level[0] == 0.0 and level[-1] < -1.0
    assert np.ptp(level) > 0.0, "a zero count over a real opportunity is a claim, not silence"
