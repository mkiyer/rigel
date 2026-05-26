"""Tests for v6 boundary evidence sweeps."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.background_model import BackgroundModel
from rigel.calibration.boundaries import BoundaryTable
from rigel.calibration.boundary_model import BoundaryLocalPosterior, build_boundary_local_posterior
from rigel.calibration.boundary_sweep import (
    FLAG_SWEEP_FROM_LEFT,
    compute_boundary_transfer_weight,
    run_boundary_sweep,
)
from rigel.calibration.density_observation import DensityObservation


def _background(alpha: float = 1.0, beta: float = 100.0) -> BackgroundModel:
    return BackgroundModel(
        rho_off_alpha=alpha,
        rho_off_beta=beta,
        rho_off_mean=alpha / beta,
        seed_mask=np.ones(3, dtype=bool),
        top_t_exclusion_mask=np.zeros(3, dtype=bool),
        n_seed_regions=3,
        n_fragments=alpha,
        eff_length=beta,
        fit_status="ok",
        flags=np.zeros(3, dtype=np.uint16),
    )


def _observation() -> DensityObservation:
    contained_leff = np.full(3, 1000.0, dtype=np.float64)
    count = np.zeros(3, dtype=np.float32)
    zero_leff = np.zeros(3, dtype=np.float64)
    mask = np.ones(3, dtype=bool)
    return DensityObservation(
        contained_count=count.copy(),
        boundary_left_count=count.copy(),
        boundary_right_count=count.copy(),
        boundary_count=count.copy(),
        observed_compatible_count=count.copy(),
        contained_leff=contained_leff,
        boundary_left_leff=zero_leff.copy(),
        boundary_right_leff=zero_leff.copy(),
        boundary_leff=zero_leff.copy(),
        anchor_intergenic=mask.copy(),
        anchor_intron=np.zeros_like(mask, dtype=bool),
        is_anchor=mask,
        spliced_count=count.copy(),
        region_length=contained_leff.copy(),
    )


def _boundaries(
    left_unspliced: list[float] | None = None,
    right_spliced: list[float] | None = None,
) -> BoundaryTable:
    boundary_count = 4
    zeros32 = np.zeros(boundary_count, dtype=np.float32)
    left_mass = zeros32.copy() if left_unspliced is None else np.asarray(left_unspliced, dtype=np.float32)
    spliced = zeros32.copy() if right_spliced is None else np.asarray(right_spliced, dtype=np.float32)
    leff = np.array([0.0, 100.0, 100.0, 0.0], dtype=np.float64)
    return BoundaryTable(
        boundary_pos=np.array([0, 100, 200, 300], dtype=np.int64),
        ref_id=np.zeros(boundary_count, dtype=np.int32),
        ref_region_offsets=np.array([0, 3], dtype=np.int64),
        ref_boundary_offsets=np.array([0, 4], dtype=np.int64),
        is_terminal=np.array([True, False, False, True]),
        left_region_unspliced_pos=left_mass,
        left_region_unspliced_neg=zeros32.copy(),
        left_region_spliced_pos=zeros32.copy(),
        left_region_spliced_neg=zeros32.copy(),
        right_region_unspliced_pos=zeros32.copy(),
        right_region_unspliced_neg=zeros32.copy(),
        right_region_spliced_pos=spliced,
        right_region_spliced_neg=zeros32.copy(),
        left_region_boundary_leff=leff,
        right_region_boundary_leff=leff.copy(),
    )


def _local(alpha: list[float], beta: list[float]) -> BoundaryLocalPosterior:
    zeros = np.zeros(3, dtype=np.float32)
    return BoundaryLocalPosterior(
        alpha_excess=np.asarray(alpha, dtype=np.float64),
        beta_excess=np.asarray(beta, dtype=np.float64),
        mu_local=zeros.copy(),
        upper_local=zeros.copy(),
        flags=np.zeros(3, dtype=np.uint16),
    )


def test_zero_transfer_weights_make_sweep_equal_local_excess() -> None:
    local = _local([10.0, 20.0, 30.0], [1.0, 2.0, 3.0])
    result = run_boundary_sweep(
        local,
        _boundaries(),
        _observation(),
        _background(),
        transfer_weight=np.zeros(4, dtype=np.float64),
    )

    np.testing.assert_allclose(result.alpha_excess, local.alpha_excess)
    np.testing.assert_allclose(result.beta_excess, local.beta_excess)
    np.testing.assert_allclose(result.forward_alpha_excess, 0.0)
    np.testing.assert_allclose(result.reverse_alpha_excess, 0.0)


def test_boundary_sweep_propagates_attenuated_evidence_to_neighbor() -> None:
    local = _local([100.0, 0.0, 0.0], [10.0, 0.0, 0.0])
    result = run_boundary_sweep(
        local,
        _boundaries(),
        _observation(),
        _background(),
        transfer_weight=np.array([0.0, 0.5, 0.0, 0.0], dtype=np.float64),
    )

    assert result.forward_alpha_excess.tolist() == [0.0, 50.0, 0.0]
    assert result.forward_beta_excess.tolist() == [0.0, 5.0, 0.0]
    assert result.alpha_excess.tolist() == [100.0, 50.0, 0.0]
    assert result.mu_sweep[1] > 10.0
    assert (result.flags[1] & FLAG_SWEEP_FROM_LEFT) != 0


def test_compute_transfer_weight_is_zero_at_terminals_and_tracks_boundary_information() -> None:
    boundaries = _boundaries(
        left_unspliced=[0.0, 100.0, 1.0, 0.0],
        right_spliced=[0.0, 0.0, 20.0, 0.0],
    )

    weights = compute_boundary_transfer_weight(boundaries, None, _background())

    assert weights[0] == pytest.approx(0.0)
    assert weights[3] == pytest.approx(0.0)
    assert weights[1] > weights[2]
    assert 0.0 <= weights[1] <= 1.0


def test_sweep_from_built_local_posterior_moves_internal_exon_like_boundary_signal() -> None:
    boundaries = _boundaries(left_unspliced=[0.0, 80.0, 0.0, 0.0])
    local = build_boundary_local_posterior(_observation(), boundaries, _background())

    result = run_boundary_sweep(
        local,
        boundaries,
        _observation(),
        _background(),
        transfer_weight=np.array([0.0, 0.25, 0.25, 0.0], dtype=np.float64),
    )

    assert local.alpha_excess.tolist() == [80.0, 0.0, 0.0]
    assert result.alpha_excess[1] == pytest.approx(20.0)
    assert result.mu_sweep[1] > local.mu_local[1]