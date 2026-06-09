"""Count overdispersion fit — the count-side twin of the strand Beta-Binomial overdispersion.

Checks the NB common-dispersion method of moments recovers a known α, the effective-count limiter
is monotone and saturates at 1/α, and the per-type fit/shrinkage behaves: small α on uniform
seeds, large α on heterogeneous (capture-like) seeds, shrink toward the pooled trend, and the
geometric α₀ fallback when no seeds exist.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.count_dispersion import (
    CountDispersionModel,
    effective_count,
    fit_gdna_count_overdispersion,
)


def _nb_sample(rng, mu, alpha, size):
    """Draw NB counts with mean ``mu`` and ``Var = mu + alpha·mu²`` (gamma-Poisson)."""
    if alpha <= 0.0:
        return rng.poisson(mu, size=size).astype(np.float64)
    r = 1.0 / alpha  # shape; Var = mu + mu²/r = mu + alpha·mu²
    p = r / (r + mu)
    return rng.negative_binomial(r, p, size=size).astype(np.float64)


def test_effective_count_monotone_and_saturates():
    n = np.array([0.0, 1.0, 10.0, 100.0, 1e6])
    eff = effective_count(n, alpha=0.5)
    assert np.all(np.diff(eff) >= 0.0)  # monotone non-decreasing in N
    assert eff[0] == 0.0
    assert eff[-1] == pytest.approx(2.0, rel=1e-3)  # → 1/α = 2 for large N


def test_effective_count_poisson_limit():
    # α = 0 ⇒ N_eff = N (the raw Poisson precision, no limiting).
    n = np.array([3.0, 50.0, 1234.0])
    assert np.allclose(effective_count(n, alpha=0.0), n)


def test_geometric_alpha_one_is_one_effective_obs():
    # α₀ = 1 (geometric) ⇒ N_eff → 1: a count is worth ~1 pseudo-observation absent evidence.
    assert effective_count(np.array([1e9]), alpha=1.0)[0] == pytest.approx(1.0, rel=1e-6)


def test_mom_recovers_known_alpha():
    # Uniform-mean NB seeds: the pooled MoM should recover the injected dispersion. Use equal eff
    # lengths so μ = ρ̄·L is constant across seeds and the trend ρ̄ matches the sampling mean.
    rng = np.random.default_rng(7)
    true_alpha = 0.3
    mu = 200.0
    counts = _nb_sample(rng, mu, true_alpha, size=4000)
    eff_len = np.full(counts.shape, 10.0)  # ρ̄ = mean(counts)/10 ⇒ μ = ρ̄·10 ≈ mu
    model = fit_gdna_count_overdispersion(
        counts, eff_len, np.array([]), np.array([]), prior_alpha=1.0, prior_weight=0.0
    )
    assert model.alpha_contained == pytest.approx(true_alpha, rel=0.15)


def test_uniform_seeds_small_alpha_heterogeneous_large():
    rng = np.random.default_rng(11)
    eff_len = np.full(2000, 10.0)
    uniform = _nb_sample(rng, 200.0, 0.001, size=2000)  # tight (≈Poisson)
    hetero = _nb_sample(rng, 200.0, 1.5, size=2000)  # capture-like overdispersion
    a_uniform = fit_gdna_count_overdispersion(
        uniform, eff_len, np.array([]), np.array([]), prior_alpha=1.0, prior_weight=0.0
    ).alpha_contained
    a_hetero = fit_gdna_count_overdispersion(
        hetero, eff_len, np.array([]), np.array([]), prior_alpha=1.0, prior_weight=0.0
    ).alpha_contained
    assert a_uniform < 0.1
    assert a_hetero > 0.5
    assert a_hetero > a_uniform


def test_two_types_distinct_and_shrink_toward_pooled():
    # Contained uniform (α≈0), crossing heterogeneous (α large). With prior_weight 0 the two types
    # stay distinct; a large prior_weight blends each toward the global pooled trend.
    rng = np.random.default_rng(13)
    cl = np.full(1500, 10.0)
    xl = np.full(1500, 5.0)
    contained = _nb_sample(rng, 200.0, 0.002, size=1500)
    crossing = _nb_sample(rng, 100.0, 1.2, size=1500)
    sharp = fit_gdna_count_overdispersion(
        contained, cl, crossing, xl, prior_alpha=1.0, prior_weight=0.0
    )
    assert sharp.alpha_contained < 0.1
    assert sharp.alpha_crossing > 0.5
    # Heavy shrink pulls the two per-type estimates toward each other (the pooled trend).
    blended = fit_gdna_count_overdispersion(
        contained, cl, crossing, xl, prior_alpha=1.0, prior_weight=1e6
    )
    assert blended.alpha_contained == pytest.approx(blended.alpha_crossing, rel=2e-2)
    assert blended.alpha_contained == pytest.approx(blended.alpha_pooled, rel=2e-2)


def test_no_seeds_falls_back_to_prior_alpha():
    model = fit_gdna_count_overdispersion(
        np.array([]), np.array([]), np.array([]), np.array([]), prior_alpha=1.0, prior_weight=30.0
    )
    assert isinstance(model, CountDispersionModel)
    assert model.fallback_used
    assert model.alpha_contained == pytest.approx(1.0)
    assert model.alpha_crossing == pytest.approx(1.0)
    assert model.n_contained_seeds == 0 and model.n_crossing_seeds == 0


def test_alpha_is_non_negative():
    # Under-dispersed (Var < μ) seeds would give a negative MoM; the floor clamps α to 0.
    # A constant count is maximally under-dispersed (zero variance).
    counts = np.full(500, 200.0)
    eff_len = np.full(500, 10.0)
    model = fit_gdna_count_overdispersion(
        counts, eff_len, np.array([]), np.array([]), prior_alpha=1.0, prior_weight=0.0
    )
    assert model.alpha_contained == 0.0
