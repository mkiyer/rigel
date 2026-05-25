"""Phase 4 tests for ``RegionExposure.from_density``.

The classmethod lifts a transient :class:`DensityEvidence` into the
consumer-side region-exposure surface that calibration hands to locus EM.

Plan reference: ``docs/fineregions/density_model_impl_plan_v4.md`` §8.3 / §8.6.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.density_model import (
    FLAG_LOW_BOUNDARY_OPPORTUNITY,
    FLAG_PRIOR_DOMINATED,
    DensityEvidence,
)
from rigel.calibration.exposure import RegionExposure


def _make_evidence(
    *,
    relative_exposure: np.ndarray,
    rho_post: np.ndarray,
    flags: np.ndarray,
    rho_ref: float = 0.05,
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


class TestRegionExposureFromDensity:
    def test_alignment_with_relative_exposure(self) -> None:
        rel = np.array([0.5, 1.0, 1.7, 0.0], dtype=np.float64)
        rho = np.array([0.025, 0.05, 0.085, 0.0], dtype=np.float64)
        ev = _make_evidence(
            relative_exposure=rel,
            rho_post=rho,
            flags=np.zeros(4, dtype=np.uint8),
        )

        exposure = RegionExposure.from_density(ev)

        assert exposure.mode == "density"
        assert exposure.A_r.dtype == np.float32
        assert exposure.rho_r.dtype == np.float32
        np.testing.assert_array_equal(
            exposure.A_r, rel.astype(np.float32)
        )
        np.testing.assert_array_equal(
            exposure.rho_r, rho.astype(np.float32)
        )
        assert exposure.rho_ref == pytest.approx(0.05)
        assert exposure.reference_quantile == 0.0

    def test_max_exposure_clips_high_regions(self) -> None:
        rel = np.array([0.5, 1.7, 3.0], dtype=np.float64)
        ev = _make_evidence(
            relative_exposure=rel,
            rho_post=np.array([0.025, 0.085, 0.15], dtype=np.float64),
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
            RegionExposure.from_density(ev, max_exposure=float("nan"))

    def test_eligibility_uses_low_boundary_opportunity(self) -> None:
        rel = np.array([0.4, 1.2, 0.9, 2.0], dtype=np.float64)
        flags = np.array(
            [
                0,
                FLAG_PRIOR_DOMINATED,
                FLAG_LOW_BOUNDARY_OPPORTUNITY,
                FLAG_LOW_BOUNDARY_OPPORTUNITY | FLAG_PRIOR_DOMINATED,
            ],
            dtype=np.uint8,
        )
        ev = _make_evidence(
            relative_exposure=rel,
            rho_post=np.array([0.02, 0.06, 0.045, 0.1], dtype=np.float64),
            flags=flags,
        )

        exposure = RegionExposure.from_density(ev)

        assert exposure.eligible.tolist() == [True, True, False, False]
        # Ineligible regions retain their posterior A_r (no overwrite).
        np.testing.assert_array_equal(
            exposure.A_r, rel.astype(np.float32)
        )
        np.testing.assert_array_equal(exposure.flags, flags)

    def test_mode_is_density(self) -> None:
        ev = _make_evidence(
            relative_exposure=np.ones(2, dtype=np.float64),
            rho_post=np.full(2, 0.05, dtype=np.float64),
            flags=np.zeros(2, dtype=np.uint8),
        )
        exposure = RegionExposure.from_density(ev)
        assert exposure.mode == "density"
