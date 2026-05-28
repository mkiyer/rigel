"""Tests for PR04 region exposure estimation."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._arrays import RegionArrays
from rigel.calibration._exposure import (
    MIN_EXPOSURE_FACTOR,
    REGION_EXPOSURE_FLOORED,
    component_bp_weighted_exposure,
)
from rigel.calibration.background_model import BackgroundModel
from rigel.calibration.calibration_iteration import (
    METHOD_BACKGROUND_FALLBACK,
    METHOD_BOUNDARY,
    METHOD_STRAND,
    BackgroundDensity,
    RegionUnsplicedMass,
)
from rigel.calibration.exposure import (
    FLAG_EXPOSURE_BOOTSTRAP_NEUTRAL,
    FLAG_EXPOSURE_IMPUTED_TIER3,
    FLAG_EXPOSURE_NOT_TAU2_POOL,
    FLAG_EXPOSURE_NO_SUPPORT,
    RegionExposure,
    estimate_region_exposure,
)


def _background_model() -> BackgroundModel:
    return BackgroundModel(
        rho_off_alpha=1.0,
        rho_off_beta=100.0,
        rho_off_mean=0.01,
        seed_mask=np.ones(3, dtype=bool),
        top_t_exclusion_mask=np.zeros(3, dtype=bool),
        n_seed_regions=3,
        n_fragments=1.0,
        eff_length=100.0,
        fit_status="ok",
        flags=np.zeros(3, dtype=np.uint16),
    )


def _density(*, fit_status: str = "converged") -> BackgroundDensity:
    return BackgroundDensity(
        rho0_mean=0.01,
        alpha0=1.0,
        beta0=100.0,
        log_dispersion=float(np.log(10.0)),
        n_effective_regions=3.0,
        n_regions_in_pool=2,
        method_histogram=(2, 0, 1),
        fit_status=fit_status,
    )


def _mass(
    *,
    gdna: list[float] | None = None,
    counts: list[int] | None = None,
    method: list[int] | None = None,
    precision: list[float] | None = None,
) -> RegionUnsplicedMass:
    gdna_arr = np.asarray(gdna if gdna is not None else [1.0, 9.0, 100.0], dtype=np.float64)
    region_count = int(gdna_arr.size)
    return RegionUnsplicedMass(
        total_mass=np.maximum(gdna_arr, 100.0),
        gdna_mass=gdna_arr,
        rna_mass=np.maximum(gdna_arr, 100.0) - gdna_arr,
        region_size_bp=np.full(region_count, 100.0, dtype=np.float64),
        unspliced_counts=np.asarray(
            counts if counts is not None else [100, 100, 100], dtype=np.uint64
        ),
        method=np.asarray(
            method
            if method is not None
            else [METHOD_STRAND, METHOD_BOUNDARY, METHOD_BACKGROUND_FALLBACK],
            dtype=np.uint8,
        ),
        precision=np.asarray(
            precision if precision is not None else [1.0, 1.0, 1.0], dtype=np.float64
        ),
        flags=np.zeros(region_count, dtype=np.uint16),
    )


def _previous_tau2(tau2: float) -> RegionExposure:
    omega = np.ones(3, dtype=np.float64)
    return RegionExposure(
        omega=omega,
        log_omega=np.log(omega),
        raw_ratio=omega,
        log_raw_ratio=np.log(omega),
        shrink_weight=np.zeros(3, dtype=np.float64),
        v_obs=np.ones(3, dtype=np.float64),
        lambda_global=np.ones(3, dtype=np.float64),
        rho0=0.01,
        tau2=tau2,
        tau2_hat=tau2,
        support_count=np.ones(3, dtype=np.uint64),
        tau2_pool_size=3,
        tau2_method="moment",
        flags=np.zeros(3, dtype=np.uint16),
    )


def test_bootstrap_density_is_exposure_neutral() -> None:
    exposure = estimate_region_exposure(
        _mass(counts=[0, 100, 100]),
        BackgroundDensity.from_bootstrap(_background_model()),
        np.ones(3, dtype=np.float64),
    )

    assert exposure.tau2_method == "bootstrap_neutral"
    assert exposure.tau2 == pytest.approx(0.0)
    np.testing.assert_allclose(exposure.omega, np.ones(3))
    np.testing.assert_allclose(exposure.shrink_weight, np.zeros(3))
    assert np.all((exposure.flags & FLAG_EXPOSURE_BOOTSTRAP_NEUTRAL) != 0)
    assert (exposure.flags[0] & FLAG_EXPOSURE_NO_SUPPORT) != 0


def test_no_pool_density_is_exposure_neutral_but_keeps_raw_diagnostics() -> None:
    density = BackgroundDensity(
        rho0_mean=0.01,
        alpha0=1.0,
        beta0=100.0,
        log_dispersion=float(np.log(10.0)),
        n_effective_regions=0.0,
        n_regions_in_pool=0,
        method_histogram=(0, 0, 3),
        fit_status="fallback_bootstrap",
    )

    exposure = estimate_region_exposure(
        _mass(
            method=[
                METHOD_BACKGROUND_FALLBACK,
                METHOD_BACKGROUND_FALLBACK,
                METHOD_BACKGROUND_FALLBACK,
            ]
        ),
        density,
        np.ones(3, dtype=np.float64),
    )

    assert exposure.tau2_method == "no_pool_neutral"
    np.testing.assert_allclose(exposure.omega, np.ones(3))
    assert exposure.raw_ratio[2] > 1.0
    assert np.all((exposure.flags & FLAG_EXPOSURE_IMPUTED_TIER3) != 0)
    assert np.all((exposure.flags & FLAG_EXPOSURE_NOT_TAU2_POOL) != 0)


def test_method_of_moments_shrinks_pool_rows_and_neutralizes_tier3() -> None:
    exposure = estimate_region_exposure(
        _mass(),
        _density(),
        np.ones(3, dtype=np.float64),
    )

    assert exposure.tau2_method == "moment"
    assert exposure.tau2_pool_size == 2
    assert exposure.tau2 > 0.0
    assert exposure.omega[0] == pytest.approx(1.0)
    assert 1.0 < exposure.omega[1] < exposure.raw_ratio[1]
    assert exposure.shrink_weight[1] > 0.0
    assert exposure.omega[2] == pytest.approx(1.0)
    assert exposure.shrink_weight[2] == pytest.approx(0.0)
    assert (exposure.flags[2] & FLAG_EXPOSURE_IMPUTED_TIER3) != 0
    assert (exposure.flags[2] & FLAG_EXPOSURE_NOT_TAU2_POOL) != 0


def test_tau2_damping_uses_previous_real_fit() -> None:
    undamped = estimate_region_exposure(
        _mass(),
        _density(),
        np.ones(3, dtype=np.float64),
    )
    damped = estimate_region_exposure(
        _mass(),
        _density(),
        np.ones(3, dtype=np.float64),
        previous=_previous_tau2(10.0),
        tau2_damping=0.25,
    )

    assert damped.tau2_method == "moment_damped"
    assert damped.tau2_hat == pytest.approx(undamped.tau2_hat)
    assert damped.tau2 == pytest.approx(0.75 * 10.0 + 0.25 * undamped.tau2_hat)


def test_invalid_scalar_parameters_raise() -> None:
    with pytest.raises(ValueError, match="alpha_floor"):
        estimate_region_exposure(
            _mass(),
            _density(),
            np.ones(3, dtype=np.float64),
            alpha_floor=0.0,
        )
    with pytest.raises(ValueError, match="omega_floor"):
        estimate_region_exposure(
            _mass(),
            _density(),
            np.ones(3, dtype=np.float64),
            omega_floor=2.0,
        )


def _region_arrays(rows: list[tuple[int, int]]) -> RegionArrays:
    return RegionArrays.from_region_df(
        pd.DataFrame(
            {
                "ref_name": pd.array(["chr1"] * len(rows), dtype="string"),
                "start": np.array([start for start, _end in rows], dtype=np.int64),
                "end": np.array([end for _start, end in rows], dtype=np.int64),
                "signature": np.zeros(len(rows), dtype=np.uint8),
            }
        ),
        {"chr1": 0},
    )


def test_component_bp_weighted_exposure_averages_over_block_overlap() -> None:
    exposure = component_bp_weighted_exposure(
        block_ref_ids=np.array([0], dtype=np.int32),
        block_starts=np.array([25], dtype=np.int64),
        block_ends=np.array([75], dtype=np.int64),
        component_ids=np.array([0], dtype=np.int64),
        n_components=1,
        region_arrays=_region_arrays([(0, 50), (50, 100)]),
        omega=np.array([2.0, 4.0], dtype=np.float64),
    )

    np.testing.assert_allclose(exposure.covered_bp, [50.0])
    np.testing.assert_allclose(exposure.weighted_bp, [150.0])
    np.testing.assert_allclose(exposure.exposure_factor, [3.0])


def test_component_bp_weighted_exposure_requires_strict_region_coverage() -> None:
    with pytest.raises(ValueError, match="does not cover"):
        component_bp_weighted_exposure(
            block_ref_ids=np.array([0], dtype=np.int32),
            block_starts=np.array([0], dtype=np.int64),
            block_ends=np.array([100], dtype=np.int64),
            component_ids=np.array([0], dtype=np.int64),
            n_components=1,
            region_arrays=_region_arrays([(0, 50)]),
            omega=np.array([1.0], dtype=np.float64),
        )


def test_component_bp_weighted_exposure_floors_positive_small_values() -> None:
    exposure = component_bp_weighted_exposure(
        block_ref_ids=np.array([0], dtype=np.int32),
        block_starts=np.array([0], dtype=np.int64),
        block_ends=np.array([100], dtype=np.int64),
        component_ids=np.array([0], dtype=np.int64),
        n_components=1,
        region_arrays=_region_arrays([(0, 100)]),
        omega=np.array([MIN_EXPOSURE_FACTOR / 10.0], dtype=np.float64),
    )

    np.testing.assert_allclose(exposure.exposure_factor, [MIN_EXPOSURE_FACTOR])
    assert (exposure.flags[0] & REGION_EXPOSURE_FLOORED) != 0
