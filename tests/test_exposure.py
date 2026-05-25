"""Phase 5 tests for :class:`rigel.calibration.exposure.RegionExposure`.

Covers the uniform constructor, the density constructor, ineligibility
wiring, ``max_exposure`` clipping, and the summary dict (including the new
``A_p99`` field added in Phase 5).

The detailed eligibility / flag handling of ``from_density`` is also
exercised in
:mod:`tests.test_region_exposure_from_density`; this file keeps the
contractual tests close to the class for editor navigation.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.density_model import (
    FLAG_LOW_BOUNDARY_OPPORTUNITY,
    DensityEvidence,
)
from rigel.calibration.exposure import RegionExposure


def _make_evidence(
    *,
    relative_exposure: np.ndarray,
    rho_post: np.ndarray,
    flags: np.ndarray,
    rho_ref: float = 0.04,
) -> DensityEvidence:
    R = int(relative_exposure.size)
    return DensityEvidence(
        rho_post=np.asarray(rho_post, dtype=np.float64),
        relative_exposure=np.asarray(relative_exposure, dtype=np.float64),
        mean_unbounded=np.zeros(R, dtype=np.float64),
        upper_unbounded=np.zeros(R, dtype=np.float64),
        prior_family=np.zeros(R, dtype=np.uint8),
        fallback_depth=np.zeros(R, dtype=np.uint8),
        flags=np.asarray(flags, dtype=np.uint8),
        confidence=0.95,
        priors={},
        rho_ref=float(rho_ref),
        rho_ref_source="INTERGENIC",
    )


class TestUniformConstructor:
    def test_shapes_and_dtypes(self) -> None:
        exposure = RegionExposure.uniform(5)
        assert exposure.mode == "uniform"
        assert exposure.A_r.dtype == np.float32
        assert exposure.rho_r.dtype == np.float32
        assert exposure.flags.dtype == np.uint8
        np.testing.assert_array_equal(exposure.A_r, np.ones(5, dtype=np.float32))
        np.testing.assert_array_equal(exposure.rho_r, np.zeros(5, dtype=np.float32))
        assert exposure.eligible.dtype == bool
        assert exposure.eligible.all()
        assert exposure.rho_ref == 0.0
        assert exposure.reference_quantile == 0.0

    def test_zero_regions(self) -> None:
        exposure = RegionExposure.uniform(0)
        assert exposure.A_r.size == 0
        assert exposure.eligible.size == 0

    def test_negative_R_raises(self) -> None:
        with pytest.raises(ValueError):
            RegionExposure.uniform(-1)


class TestFromDensity:
    def test_alignment(self) -> None:
        rel = np.array([0.4, 1.2, 2.0], dtype=np.float64)
        rho = np.array([0.016, 0.048, 0.080], dtype=np.float64)
        ev = _make_evidence(
            relative_exposure=rel,
            rho_post=rho,
            flags=np.zeros(3, dtype=np.uint8),
        )
        exposure = RegionExposure.from_density(ev)
        assert exposure.mode == "density"
        np.testing.assert_array_equal(exposure.A_r, rel.astype(np.float32))
        np.testing.assert_array_equal(exposure.rho_r, rho.astype(np.float32))
        assert exposure.rho_ref == pytest.approx(0.04)

    def test_ineligibility_preserves_A_r(self) -> None:
        rel = np.array([0.5, 2.0], dtype=np.float64)
        flags = np.array(
            [0, FLAG_LOW_BOUNDARY_OPPORTUNITY], dtype=np.uint8
        )
        ev = _make_evidence(
            relative_exposure=rel,
            rho_post=np.array([0.02, 0.08], dtype=np.float64),
            flags=flags,
        )
        exposure = RegionExposure.from_density(ev)
        assert exposure.eligible.tolist() == [True, False]
        # Ineligible row keeps posterior A_r — no overwrite to 1.0.
        np.testing.assert_array_equal(exposure.A_r, rel.astype(np.float32))

    def test_max_exposure_clips(self) -> None:
        rel = np.array([0.5, 1.5, 3.0], dtype=np.float64)
        ev = _make_evidence(
            relative_exposure=rel,
            rho_post=np.array([0.02, 0.06, 0.12], dtype=np.float64),
            flags=np.zeros(3, dtype=np.uint8),
        )
        exposure = RegionExposure.from_density(ev, max_exposure=1.5)
        np.testing.assert_array_equal(
            exposure.A_r,
            np.array([0.5, 1.5, 1.5], dtype=np.float32),
        )

    def test_max_exposure_invalid(self) -> None:
        ev = _make_evidence(
            relative_exposure=np.ones(1, dtype=np.float64),
            rho_post=np.zeros(1, dtype=np.float64),
            flags=np.zeros(1, dtype=np.uint8),
        )
        with pytest.raises(ValueError):
            RegionExposure.from_density(ev, max_exposure=0.0)
        with pytest.raises(ValueError):
            RegionExposure.from_density(ev, max_exposure=-1.0)
        with pytest.raises(ValueError):
            RegionExposure.from_density(ev, max_exposure=float("nan"))


class TestSummaryDict:
    def test_uniform_summary(self) -> None:
        summary = RegionExposure.uniform(4).to_summary_dict()
        assert summary["mode"] == "uniform"
        assert summary["n_regions"] == 4
        assert summary["n_regions_eligible"] == 4
        assert summary["A_min"] == 1.0
        assert summary["A_mean"] == 1.0
        assert summary["A_p99"] == 1.0
        assert summary["A_max"] == 1.0

    def test_density_summary_includes_p99(self) -> None:
        rel = np.linspace(0.1, 3.0, 100, dtype=np.float64)
        ev = _make_evidence(
            relative_exposure=rel,
            rho_post=rel * 0.04,
            flags=np.zeros(100, dtype=np.uint8),
        )
        summary = RegionExposure.from_density(ev).to_summary_dict()
        assert summary["mode"] == "density"
        assert summary["A_min"] == pytest.approx(0.1, abs=1e-6)
        assert summary["A_max"] == pytest.approx(3.0, abs=1e-6)
        # p99 lies between mean and max for an evenly spaced grid.
        assert summary["A_mean"] < summary["A_p99"] <= summary["A_max"]
        # p99 of linspace(0.1, 3.0, 100) is the 99th-percentile sample.
        assert summary["A_p99"] == pytest.approx(
            float(np.quantile(rel, 0.99)), abs=1e-6
        )

    def test_zero_region_summary(self) -> None:
        summary = RegionExposure.uniform(0).to_summary_dict()
        assert summary["n_regions"] == 0
        assert summary["A_p99"] == 1.0
