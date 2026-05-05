"""Tests for ``rigel.calibration.density_loco.shrink_to_loco``."""

from __future__ import annotations

import math

import pytest

from rigel.calibration.density_loco import shrink_to_loco


def test_kappa_zero_returns_pure_local():
    # κ → 0 ⇒ N / L_eff
    assert shrink_to_loco(10.0, 100.0, 0.05, 0.0) == pytest.approx(0.1)


def test_kappa_large_returns_global():
    # κ → ∞ ⇒ ρ_global
    out = shrink_to_loco(10.0, 100.0, 0.05, 1.0e9)
    assert out == pytest.approx(0.05, rel=1e-6)


def test_zero_leff_returns_global():
    assert shrink_to_loco(0.0, 0.0, 0.07, 0.0) == 0.07
    assert shrink_to_loco(5.0, 0.0, 0.07, 0.0) == 0.07


def test_zero_local_count_shrinks_toward_global():
    # N=0, L_eff>0 ⇒ ρ_global · κ / (L_eff + κ)  (closed form)
    out = shrink_to_loco(0.0, 100.0, 0.05, 50.0)
    assert out == pytest.approx(0.05 * 50.0 / (100.0 + 50.0))


def test_monotone_in_n():
    a = shrink_to_loco(10.0, 100.0, 0.05, 50.0)
    b = shrink_to_loco(20.0, 100.0, 0.05, 50.0)
    assert b > a


def test_negative_input_raises():
    with pytest.raises(ValueError):
        shrink_to_loco(-1.0, 100.0, 0.05, 50.0)
    with pytest.raises(ValueError):
        shrink_to_loco(0.0, -1.0, 0.05, 50.0)
    with pytest.raises(ValueError):
        shrink_to_loco(0.0, 100.0, -0.05, 50.0)
    with pytest.raises(ValueError):
        shrink_to_loco(0.0, 100.0, 0.05, -1.0)


def test_nan_propagates_via_validation():
    with pytest.raises(ValueError):
        shrink_to_loco(math.nan, 100.0, 0.05, 50.0)
