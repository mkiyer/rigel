"""RNA strand-balance fit: posterior-predictive Beta-Binomial from the spliced 2x2 (PR 9)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from rigel.calibration.strand_balance import StrandBalance, fit_strand_balance


def _strand_model(p_r1_sense: float, n_observations: int):
    """Minimal StrandModels stand-in: the posterior-predictive needs only these two."""
    return SimpleNamespace(p_r1_sense=p_r1_sense, n_observations=n_observations)


def test_posterior_mean_and_overdispersion():
    # n_obs=10, p_r1_sense=0.8 → n_same=8 → Beta(9, 3): κ=9/12=0.75; ρ=1/(10+3)=1/13.
    sb = fit_strand_balance(_strand_model(0.8, 10))
    np.testing.assert_allclose(sb.rna_sense_frac, 9.0 / 12.0)
    np.testing.assert_allclose(sb.rna_strand_overdispersion, 1.0 / 13.0)
    assert sb.n_observations == 10
    assert not sb.fallback_used


def test_dense_collapses_toward_binomial():
    # Abundant spliced reads → overdispersion → 0 (the Binomial limit) and κ → p.
    sb = fit_strand_balance(_strand_model(0.8, 200_000))
    assert sb.rna_strand_overdispersion < 1e-4
    np.testing.assert_allclose(sb.rna_sense_frac, 0.8, atol=1e-4)


def test_sparse_is_overdispersed():
    # Few spliced reads → large overdispersion (robust — not a razor-sharp Binomial).
    sb1 = fit_strand_balance(_strand_model(1.0, 1))  # n_same=1 → Beta(2, 1)
    sb5 = fit_strand_balance(_strand_model(1.0, 5))
    np.testing.assert_allclose(sb1.rna_strand_overdispersion, 1.0 / 4.0)
    np.testing.assert_allclose(sb5.rna_strand_overdispersion, 1.0 / 8.0)
    assert sb1.rna_strand_overdispersion > sb5.rna_strand_overdispersion  # sparser → wider


def test_overdispersion_monotone_decreasing_in_power():
    rhos = [
        fit_strand_balance(_strand_model(0.7, n)).rna_strand_overdispersion
        for n in (2, 10, 100, 1000)
    ]
    assert all(a > b for a, b in zip(rhos, rhos[1:]))


def test_kappa_never_at_the_bound():
    # Even a "perfectly stranded" point estimate is pulled off 0/1 by the Beta(1,1) prior.
    for p in (0.0, 1.0):
        sb = fit_strand_balance(_strand_model(p, 10))
        assert 0.0 < sb.rna_sense_frac < 1.0
    np.testing.assert_allclose(
        fit_strand_balance(_strand_model(1.0, 10)).rna_sense_frac, 11.0 / 12.0
    )


def test_zero_spliced_is_symmetric_fallback():
    # No spliced reads → Beta(1,1): κ=0.5 (channel neutral), ρ=1/3, fallback flagged.
    sb = fit_strand_balance(_strand_model(0.5, 0))
    np.testing.assert_allclose(sb.rna_sense_frac, 0.5)
    np.testing.assert_allclose(sb.rna_strand_overdispersion, 1.0 / 3.0)
    assert sb.fallback_used and sb.n_observations == 0


def test_returns_strand_balance_type():
    sb = fit_strand_balance(_strand_model(0.8, 20))
    assert isinstance(sb, StrandBalance)
    assert 0.0 < sb.rna_sense_frac < 1.0
    assert 0.0 < sb.rna_strand_overdispersion < 1.0
