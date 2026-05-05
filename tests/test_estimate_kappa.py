"""Tests for ``rigel.calibration._kappa.estimate_kappa``.

Pins the per-region NB MoM estimator and its fallback policy.

See ``docs/calibration/m4_implementation_plan.md`` §5.2.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration._kappa import (
    KAPPA_DEFAULT,
    KAPPA_MAX,
    KAPPA_MIN,
    MIN_REGIONS_FOR_MOM,
    KappaEstimate,
    estimate_kappa,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _draw_nb(mu: np.ndarray, kappa: float, rng: np.random.Generator) -> np.ndarray:
    """Draw NB(mu, kappa) counts using numpy's (n=kappa, p=kappa/(kappa+mu))."""
    p = kappa / (kappa + mu)
    return rng.negative_binomial(n=kappa, p=p).astype(np.int64)


# ---------------------------------------------------------------------------
# Degenerate inputs → fallback
# ---------------------------------------------------------------------------


class TestFallbacks:
    def test_too_few_regions(self):
        n = MIN_REGIONS_FOR_MOM - 1
        ke = estimate_kappa(
            np.zeros(n, dtype=np.int64),
            np.full(n, 1000.0),
            rho_hat=1e-4,
        )
        assert ke.fallback_used
        assert ke.value == KAPPA_DEFAULT
        assert "MIN_REGIONS_FOR_MOM" in ke.fallback_reason
        assert ke.n_regions == n

    def test_rho_zero(self):
        ke = estimate_kappa(
            np.zeros(50, dtype=np.int64),
            np.full(50, 1000.0),
            rho_hat=0.0,
        )
        assert ke.fallback_used
        assert ke.value == KAPPA_DEFAULT
        assert "rho_hat" in ke.fallback_reason

    def test_rho_negative(self):
        ke = estimate_kappa(
            np.zeros(50, dtype=np.int64),
            np.full(50, 1000.0),
            rho_hat=-1e-3,
        )
        assert ke.fallback_used
        assert "rho_hat" in ke.fallback_reason

    def test_all_zero_eff_lengths(self):
        ke = estimate_kappa(
            np.zeros(50, dtype=np.int64),
            np.zeros(50, dtype=np.float64),
            rho_hat=1e-4,
        )
        assert ke.fallback_used
        assert "eff_lengths" in ke.fallback_reason

    def test_underdispersed_constant_counts(self):
        # All counts equal to mu (Poisson-like with σ²=0 < μ) → excess <= 0.
        eff = np.full(50, 1000.0)
        rho = 0.01
        mu = rho * eff[0]
        counts = np.full(50, int(round(mu)), dtype=np.int64)
        ke = estimate_kappa(counts, eff, rho)
        assert ke.fallback_used
        assert "excess variance" in ke.fallback_reason


# ---------------------------------------------------------------------------
# Recovery from heterogeneous-exposure NB draws
# ---------------------------------------------------------------------------


class TestRecovery:
    @pytest.mark.parametrize("kappa_true", [2.0, 20.0, 200.0])
    def test_recovery_within_30pct(self, kappa_true):
        # Per-region mean count μ_R must be large enough that the MoM
        # signal Σμ²/κ stands out above Poisson noise on Σ(N-μ)².  With
        # 500 regions and μ ~ 50 the signal/noise ratio is ample for all
        # three κ values; with μ ~ 0.1 (the original spec) κ=200 is
        # unidentifiable in finite samples.
        rng = np.random.default_rng(seed=int(kappa_true * 1000))
        n = 500
        eff = rng.lognormal(mean=np.log(50_000.0), sigma=0.8, size=n)
        rho = 1e-3
        mu = rho * eff
        counts = _draw_nb(mu, kappa_true, rng)
        ke = estimate_kappa(counts, eff, rho)
        assert not ke.fallback_used, ke.fallback_reason
        rel = abs(ke.value - kappa_true) / kappa_true
        assert rel < 0.30, (
            f"kappa_true={kappa_true} estimate={ke.value:.3f} rel_err={rel:.3f}"
        )

    def test_clipping_at_max(self):
        # Tight Poisson-like data should give very large κ — clip to KAPPA_MAX.
        rng = np.random.default_rng(seed=7)
        n = 500
        eff = np.full(n, 50_000.0)
        rho = 1e-3
        mu = rho * eff
        counts = rng.poisson(mu).astype(np.int64)
        ke = estimate_kappa(counts, eff, rho)
        if not ke.fallback_used:
            assert ke.value <= KAPPA_MAX
            assert ke.value >= KAPPA_MIN

    def test_clipping_at_min(self):
        # Heavily over-dispersed data (κ → 0); MoM may go below 1 → clipped.
        rng = np.random.default_rng(seed=9)
        n = 500
        eff = np.full(n, 50_000.0)
        rho = 1e-3
        mu = rho * eff
        counts = _draw_nb(mu, kappa=0.5, rng=rng)
        ke = estimate_kappa(counts, eff, rho)
        if not ke.fallback_used:
            assert ke.value >= KAPPA_MIN
            assert ke.value <= KAPPA_MAX


# ---------------------------------------------------------------------------
# Schema invariants
# ---------------------------------------------------------------------------


class TestSchema:
    def test_determinism(self):
        rng = np.random.default_rng(seed=42)
        eff = rng.lognormal(mean=np.log(1000.0), sigma=0.5, size=200)
        rho = 1e-4
        counts = _draw_nb(rho * eff, kappa=20.0, rng=rng)
        a = estimate_kappa(counts, eff, rho)
        b = estimate_kappa(counts, eff, rho)
        assert a == b

    def test_no_fallback_means_empty_reason(self):
        rng = np.random.default_rng(seed=1)
        eff = rng.lognormal(mean=np.log(1000.0), sigma=0.5, size=200)
        rho = 1e-4
        counts = _draw_nb(rho * eff, kappa=20.0, rng=rng)
        ke = estimate_kappa(counts, eff, rho)
        if not ke.fallback_used:
            assert ke.fallback_reason == ""

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="shape"):
            estimate_kappa(
                np.zeros(10, dtype=np.int64),
                np.zeros(11, dtype=np.float64),
                rho_hat=1e-4,
            )

    def test_returns_kappaestimate(self):
        ke = estimate_kappa(
            np.zeros(2, dtype=np.int64),
            np.zeros(2, dtype=np.float64),
            rho_hat=1e-4,
        )
        assert isinstance(ke, KappaEstimate)
        assert isinstance(ke.value, float)
        assert isinstance(ke.n_regions, int)
        assert isinstance(ke.fallback_used, bool)
        assert isinstance(ke.fallback_reason, str)
