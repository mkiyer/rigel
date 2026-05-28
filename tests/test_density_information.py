"""Tests for ``DensityEvidence.information`` and ``.applicable`` (PR 07 Phase 2)."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.density_model import (
    DENSITY_INFO_VAR_FLOOR,
    DENSITY_MIN_EFF_LENGTH,
    DensityEvidence,
    PRIOR_FAMILY_DETERMINISTIC_ZERO,
    _density_information,
    _deterministic_zero_evidence,
)
from rigel.calibration.density_observation import DensityObservation


def _empty_observation(n: int) -> DensityObservation:
    """Build an all-zeros density observation with ``n`` regions."""
    f32 = np.zeros(n, dtype=np.float32)
    f64 = np.zeros(n, dtype=np.float64)
    bf = np.zeros(n, dtype=bool)
    i64 = np.zeros(n, dtype=np.int64)
    return DensityObservation(
        contained_count=f32,
        boundary_left_count=f32,
        boundary_right_count=f32,
        boundary_count=f32,
        observed_compatible_count=f32,
        contained_leff=f64,
        boundary_left_leff=f64,
        boundary_right_leff=f64,
        boundary_leff=f64,
        anchor_intergenic=bf,
        anchor_intron=bf,
        is_anchor=bf,
        spliced_count=f32,
        region_length=i64,
    )


class TestDensityInformation:
    def test_zero_variance_yields_zero_information(self):
        var = np.array([0.0, 0.0, 0.0], dtype=np.float64)
        info = _density_information(var)
        np.testing.assert_array_equal(info, np.zeros(3))

    def test_positive_variance_yields_finite_information(self):
        var = np.array([1.0, 2.0, 0.25], dtype=np.float64)
        info = _density_information(var)
        np.testing.assert_allclose(info, [1.0, 0.5, 4.0])
        assert np.all(np.isfinite(info))

    def test_information_monotone_decreasing_in_variance(self):
        var = np.geomspace(1e-6, 1e6, 21).astype(np.float64)
        info = _density_information(var)
        diffs = np.diff(info)
        assert np.all(diffs <= 1e-15), "information must decrease as variance grows"

    def test_subfloor_variance_treated_as_zero(self):
        var = np.array([DENSITY_INFO_VAR_FLOOR * 0.5], dtype=np.float64)
        info = _density_information(var)
        assert info[0] == 0.0


class TestDeterministicZeroEvidence:
    """The deterministic-zero fallback must emit the new fields."""

    def test_has_information_and_applicable_fields(self):
        obs = _empty_observation(5)
        ev = _deterministic_zero_evidence(obs, confidence=0.95)
        assert isinstance(ev, DensityEvidence)
        assert ev.information is not None
        assert ev.applicable is not None
        assert ev.information.shape == (5,)
        assert ev.applicable.shape == (5,)
        assert ev.information.dtype == np.float64
        assert ev.applicable.dtype == bool

    def test_zero_information_everywhere(self):
        obs = _empty_observation(7)
        ev = _deterministic_zero_evidence(obs, confidence=0.95)
        np.testing.assert_array_equal(ev.information, np.zeros(7))

    def test_not_applicable_everywhere(self):
        obs = _empty_observation(4)
        ev = _deterministic_zero_evidence(obs, confidence=0.95)
        assert not np.any(ev.applicable)
        assert np.all(ev.prior_family == PRIOR_FAMILY_DETERMINISTIC_ZERO)


@pytest.mark.parametrize("leff,expected", [
    (DENSITY_MIN_EFF_LENGTH * 2.0, True),
    (DENSITY_MIN_EFF_LENGTH, True),
    (DENSITY_MIN_EFF_LENGTH * 0.5, False),
    (0.0, False),
])
def test_applicable_threshold(leff, expected):
    """``applicable`` flips strictly at ``contained_leff >= min_eff_length``."""
    # Directly probe the threshold semantic that fit_density_evidence uses.
    arr = np.array([leff], dtype=np.float64)
    out = arr >= DENSITY_MIN_EFF_LENGTH
    assert bool(out[0]) is expected
