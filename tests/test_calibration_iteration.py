"""Tests for Phase III four-state calibration iteration."""

from __future__ import annotations

import numpy as np

from rigel.calibration._arrays import RegionArrays
from rigel.calibration.background_model import BackgroundModel
from rigel.calibration.boundaries import BoundaryTable
from rigel.calibration.boundary_model import BoundaryLocalPosterior
from rigel.calibration.calibration_iteration import (
    FLAG_BOUNDARY_SPARSE,
    RegionCalibration,
    calibration_e_step,
    run_calibration_iteration,
)
from rigel.calibration.density_observation import DensityObservation
from rigel.calibration.fractional_evidence import transcript_strand_class
from rigel.calibration.latent_states import (
    STATE_BACKGROUND,
    STATE_EXPRESSED_CAPTURE,
    STATE_EXPRESSED_OFFTARGET,
    STATE_GDNA_ONLY_CAPTURE,
)
from rigel.calibration.signature import pack_signature


def _region_arrays() -> RegionArrays:
    signatures = np.asarray(
        [
            pack_signature(),
            pack_signature(),
            pack_signature(exon_pos=True),
            pack_signature(exon_pos=True),
        ],
        dtype=np.uint8,
    )
    starts = np.arange(4, dtype=np.int64) * 100
    return RegionArrays(
        ref_id=np.zeros(4, dtype=np.int32),
        start=starts,
        end=starts + 100,
        signature=signatures,
        ts_class=transcript_strand_class(signatures),
        ref_offsets=np.array([0, 4], dtype=np.int32),
        order=np.arange(4, dtype=np.int64),
        n_refs=1,
    )


def _observation() -> DensityObservation:
    contained = np.asarray([0.0, 30.0, 0.0, 30.0], dtype=np.float32)
    leff = np.full(4, 1000.0, dtype=np.float64)
    spliced = np.asarray([0.0, 0.0, 20.0, 20.0], dtype=np.float32)
    zeros_count = np.zeros(4, dtype=np.float32)
    zeros_leff = np.zeros(4, dtype=np.float64)
    return DensityObservation(
        contained_count=contained,
        boundary_left_count=zeros_count.copy(),
        boundary_right_count=zeros_count.copy(),
        boundary_count=zeros_count.copy(),
        observed_compatible_count=contained.copy(),
        contained_leff=leff,
        boundary_left_leff=zeros_leff.copy(),
        boundary_right_leff=zeros_leff.copy(),
        boundary_leff=zeros_leff.copy(),
        anchor_intergenic=np.array([True, True, False, False]),
        anchor_intron=np.zeros(4, dtype=bool),
        is_anchor=np.array([True, True, False, False]),
        spliced_count=spliced,
        region_length=leff.copy(),
    )


def _background() -> BackgroundModel:
    return BackgroundModel(
        rho_off_alpha=1.0,
        rho_off_beta=100.0,
        rho_off_mean=0.01,
        seed_mask=np.array([True, False, False, False]),
        top_t_exclusion_mask=np.zeros(4, dtype=bool),
        n_seed_regions=1,
        n_fragments=1.0,
        eff_length=100.0,
        fit_status="sparse",
        flags=np.zeros(4, dtype=np.uint16),
    )


def _boundaries() -> BoundaryTable:
    boundary_count = 5
    zeros32 = np.zeros(boundary_count, dtype=np.float32)
    zeros64 = np.zeros(boundary_count, dtype=np.float64)
    return BoundaryTable(
        boundary_pos=np.array([0, 100, 200, 300, 400], dtype=np.int64),
        ref_id=np.zeros(boundary_count, dtype=np.int32),
        ref_region_offsets=np.array([0, 4], dtype=np.int64),
        ref_boundary_offsets=np.array([0, 5], dtype=np.int64),
        is_terminal=np.array([True, False, False, False, True]),
        left_region_unspliced_pos=zeros32.copy(),
        left_region_unspliced_neg=zeros32.copy(),
        left_region_spliced_pos=zeros32.copy(),
        left_region_spliced_neg=zeros32.copy(),
        right_region_unspliced_pos=zeros32.copy(),
        right_region_unspliced_neg=zeros32.copy(),
        right_region_spliced_pos=zeros32.copy(),
        right_region_spliced_neg=zeros32.copy(),
        left_region_boundary_leff=zeros64.copy(),
        right_region_boundary_leff=zeros64.copy(),
    )


def _local_posterior() -> BoundaryLocalPosterior:
    return BoundaryLocalPosterior(
        alpha_excess=np.asarray([0.0, 5.0, 0.0, 5.0], dtype=np.float64),
        beta_excess=np.asarray([0.0, 100.0, 0.0, 100.0], dtype=np.float64),
        mu_local=np.asarray([10.0, 30.0, 10.0, 30.0], dtype=np.float32),
        upper_local=np.asarray([20.0, 45.0, 20.0, 45.0], dtype=np.float32),
        flags=np.zeros(4, dtype=np.uint16),
    )


def test_calibration_e_step_identifies_four_archetypes() -> None:
    result = calibration_e_step(
        _region_arrays(),
        _observation(),
        _boundaries(),
        _background(),
        local_posterior=_local_posterior(),
        transfer_weight=np.zeros(5, dtype=np.float64),
    )

    dominant = np.argmax(result.p_states, axis=1).tolist()
    assert dominant == [
        STATE_BACKGROUND,
        STATE_GDNA_ONLY_CAPTURE,
        STATE_EXPRESSED_OFFTARGET,
        STATE_EXPRESSED_CAPTURE,
    ]
    np.testing.assert_allclose(result.p_states.sum(axis=1), 1.0, rtol=1.0e-6, atol=1.0e-6)
    assert np.all(result.mu_gdna >= 0.0)
    assert np.all(result.upper_gdna >= result.mu_gdna)
    assert np.all(np.isfinite(result.A_r))
    assert np.all(result.A_r >= 0.0)
    assert np.all(np.isfinite(result.gamma_r))
    assert (result.flags[0] & FLAG_BOUNDARY_SPARSE) != 0


def test_calibration_iteration_converges_with_fixed_scalar_parameters() -> None:
    calibration = run_calibration_iteration(
        _region_arrays(),
        _observation(),
        _boundaries(),
        _background(),
        local_posterior=_local_posterior(),
        transfer_weight=np.zeros(5, dtype=np.float64),
        max_calibration_passes=5,
        damping=0.0,
    )

    assert isinstance(calibration, RegionCalibration)
    assert calibration.converged
    assert calibration.n_passes == 2
    assert len(calibration.pass_diagnostics) == calibration.n_passes
    assert calibration.pass_diagnostics[-1]["converged"] is True
    np.testing.assert_allclose(calibration.p_states.sum(axis=1), 1.0, rtol=1.0e-6, atol=1.0e-6)
    dominant = np.argmax(calibration.p_states, axis=1).tolist()
    assert dominant == [
        STATE_BACKGROUND,
        STATE_GDNA_ONLY_CAPTURE,
        STATE_EXPRESSED_OFFTARGET,
        STATE_EXPRESSED_CAPTURE,
    ]
    assert np.all(calibration.upper_gdna >= calibration.mu_gdna)
    assert np.all(np.isfinite(calibration.A_r))
    assert np.all(calibration.A_r >= 0.0)
