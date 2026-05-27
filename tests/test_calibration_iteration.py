"""Tests for Phase III four-state calibration iteration."""

from __future__ import annotations

import numpy as np

from rigel.calibration._arrays import RegionArrays
from rigel.calibration.background_model import BackgroundModel
from rigel.calibration.boundaries import BoundaryTable
from rigel.calibration.boundary_model import BoundaryLocalPosterior
from rigel.calibration.calibration_iteration import (
    FLAG_BOUNDARY_SPARSE,
    PRIOR_MASS_METHOD_DENSITY,
    PRIOR_MASS_METHOD_STRAND,
    RegionCalibration,
    build_prior_mass_deconvolution,
    calibration_e_step,
    run_calibration_iteration,
)
from rigel.calibration.density_observation import DensityObservation
from rigel.calibration.fractional_evidence import transcript_strand_class
from rigel.calibration.latent_states import (
    STATE_EXPRESSED_CAPTURE,
    STATE_EXPRESSED_OFFTARGET,
    STATE_UNEXPRESSED_CAPTURE,
    STATE_UNEXPRESSED_OFFTARGET,
)
from rigel.calibration.signature import pack_signature
from rigel.calibration.strand_deconv import RegionGdnaChannelEstimate


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
        STATE_UNEXPRESSED_OFFTARGET,
        STATE_UNEXPRESSED_CAPTURE,
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
        STATE_UNEXPRESSED_OFFTARGET,
        STATE_UNEXPRESSED_CAPTURE,
        STATE_EXPRESSED_OFFTARGET,
        STATE_EXPRESSED_CAPTURE,
    ]
    assert np.all(calibration.upper_gdna >= calibration.mu_gdna)
    assert np.all(np.isfinite(calibration.A_r))
    assert np.all(calibration.A_r >= 0.0)


def test_prior_mass_deconvolution_excludes_spliced_mass_in_density_fallback() -> None:
    observation = DensityObservation(
        contained_count=np.array([10.0], dtype=np.float32),
        boundary_left_count=np.array([20.0], dtype=np.float32),
        boundary_right_count=np.array([30.0], dtype=np.float32),
        boundary_count=np.array([50.0], dtype=np.float32),
        observed_compatible_count=np.array([60.0], dtype=np.float32),
        contained_leff=np.array([100.0], dtype=np.float64),
        boundary_left_leff=np.array([10.0], dtype=np.float64),
        boundary_right_leff=np.array([10.0], dtype=np.float64),
        boundary_leff=np.array([20.0], dtype=np.float64),
        anchor_intergenic=np.array([False]),
        anchor_intron=np.array([False]),
        is_anchor=np.array([False]),
        spliced_count=np.array([999.0], dtype=np.float32),
        region_length=np.array([100.0], dtype=np.float64),
    )

    prior_mass = build_prior_mass_deconvolution(
        observation,
        mu_gdna=np.array([12.0], dtype=np.float32),
        strand_channels=None,
    )

    np.testing.assert_allclose(prior_mass.unspliced_total, [60.0])
    np.testing.assert_allclose(prior_mass.gdna_unspliced_mean, [12.0])
    np.testing.assert_allclose(prior_mass.rna_unspliced_mean, [48.0])
    assert prior_mass.method.tolist() == [PRIOR_MASS_METHOD_DENSITY]
    assert float(prior_mass.mass_conservation_error()[0]) == 0.0


def test_prior_mass_deconvolution_uses_compartment_strand_means_when_available() -> None:
    observation = DensityObservation(
        contained_count=np.array([10.0], dtype=np.float32),
        boundary_left_count=np.array([20.0], dtype=np.float32),
        boundary_right_count=np.array([30.0], dtype=np.float32),
        boundary_count=np.array([50.0], dtype=np.float32),
        observed_compatible_count=np.array([60.0], dtype=np.float32),
        contained_leff=np.array([100.0], dtype=np.float64),
        boundary_left_leff=np.array([10.0], dtype=np.float64),
        boundary_right_leff=np.array([10.0], dtype=np.float64),
        boundary_leff=np.array([20.0], dtype=np.float64),
        anchor_intergenic=np.array([False]),
        anchor_intron=np.array([False]),
        is_anchor=np.array([False]),
        spliced_count=np.array([999.0], dtype=np.float32),
        region_length=np.array([100.0], dtype=np.float64),
    )
    strand_channels = RegionGdnaChannelEstimate(
        contained_mean=np.array([1.0], dtype=np.float32),
        contained_upper=np.array([2.0], dtype=np.float32),
        contained_rna_lower=np.array([8.0], dtype=np.float32),
        contained_precision=np.array([1.0], dtype=np.float32),
        boundary_left_mean=np.array([2.0], dtype=np.float32),
        boundary_left_upper=np.array([3.0], dtype=np.float32),
        boundary_left_rna_lower=np.array([17.0], dtype=np.float32),
        boundary_left_precision=np.array([1.0], dtype=np.float32),
        boundary_right_mean=np.array([3.0], dtype=np.float32),
        boundary_right_upper=np.array([4.0], dtype=np.float32),
        boundary_right_rna_lower=np.array([26.0], dtype=np.float32),
        boundary_right_precision=np.array([1.0], dtype=np.float32),
        flags=np.zeros(1, dtype=np.uint16),
        kappa_d=10.0,
        p_r1_sense=0.95,
        internal_rna_lower_ci=0.95,
        contained_reliability=np.ones(1, dtype=np.float32),
        contained_log_bayes_factor=np.zeros(1, dtype=np.float32),
        boundary_left_reliability=np.ones(1, dtype=np.float32),
        boundary_left_log_bayes_factor=np.zeros(1, dtype=np.float32),
        boundary_right_reliability=np.ones(1, dtype=np.float32),
        boundary_right_log_bayes_factor=np.zeros(1, dtype=np.float32),
    )

    prior_mass = build_prior_mass_deconvolution(
        observation,
        mu_gdna=np.array([99.0], dtype=np.float32),
        strand_channels=strand_channels,
    )

    np.testing.assert_allclose(prior_mass.unspliced_total, [60.0])
    np.testing.assert_allclose(prior_mass.gdna_unspliced_mean, [6.0])
    np.testing.assert_allclose(prior_mass.rna_unspliced_mean, [54.0])
    assert prior_mass.method.tolist() == [PRIOR_MASS_METHOD_STRAND]
    assert float(prior_mass.mass_conservation_error()[0]) == 0.0


def test_prior_mass_deconvolution_weights_strand_means_by_reliability() -> None:
    observation = DensityObservation(
        contained_count=np.array([100.0], dtype=np.float32),
        boundary_left_count=np.array([0.0], dtype=np.float32),
        boundary_right_count=np.array([0.0], dtype=np.float32),
        boundary_count=np.array([0.0], dtype=np.float32),
        observed_compatible_count=np.array([100.0], dtype=np.float32),
        contained_leff=np.array([100.0], dtype=np.float64),
        boundary_left_leff=np.array([0.0], dtype=np.float64),
        boundary_right_leff=np.array([0.0], dtype=np.float64),
        boundary_leff=np.array([0.0], dtype=np.float64),
        anchor_intergenic=np.array([False]),
        anchor_intron=np.array([False]),
        is_anchor=np.array([False]),
        spliced_count=np.array([0.0], dtype=np.float32),
        region_length=np.array([100.0], dtype=np.float64),
    )
    strand_channels = RegionGdnaChannelEstimate(
        contained_mean=np.array([80.0], dtype=np.float32),
        contained_upper=np.array([90.0], dtype=np.float32),
        contained_rna_lower=np.array([10.0], dtype=np.float32),
        contained_precision=np.array([0.8], dtype=np.float32),
        boundary_left_mean=np.array([0.0], dtype=np.float32),
        boundary_left_upper=np.array([0.0], dtype=np.float32),
        boundary_left_rna_lower=np.array([0.0], dtype=np.float32),
        boundary_left_precision=np.array([0.0], dtype=np.float32),
        boundary_right_mean=np.array([0.0], dtype=np.float32),
        boundary_right_upper=np.array([0.0], dtype=np.float32),
        boundary_right_rna_lower=np.array([0.0], dtype=np.float32),
        boundary_right_precision=np.array([0.0], dtype=np.float32),
        flags=np.zeros(1, dtype=np.uint16),
        kappa_d=10.0,
        p_r1_sense=0.95,
        internal_rna_lower_ci=0.95,
        contained_reliability=np.array([0.25], dtype=np.float32),
        contained_log_bayes_factor=np.zeros(1, dtype=np.float32),
        boundary_left_reliability=np.zeros(1, dtype=np.float32),
        boundary_left_log_bayes_factor=np.zeros(1, dtype=np.float32),
        boundary_right_reliability=np.zeros(1, dtype=np.float32),
        boundary_right_log_bayes_factor=np.zeros(1, dtype=np.float32),
    )

    prior_mass = build_prior_mass_deconvolution(
        observation,
        mu_gdna=np.array([100.0], dtype=np.float32),
        strand_channels=strand_channels,
    )

    np.testing.assert_allclose(prior_mass.unspliced_total, [100.0])
    np.testing.assert_allclose(prior_mass.gdna_unspliced_mean, [20.0])
    np.testing.assert_allclose(prior_mass.rna_unspliced_mean, [80.0])
    np.testing.assert_allclose(prior_mass.precision, [0.2])
    assert float(prior_mass.mass_conservation_error()[0]) == 0.0


def test_prior_mass_deconvolution_zero_reliability_blocks_strand_gdna() -> None:
    observation = DensityObservation(
        contained_count=np.array([50.0], dtype=np.float32),
        boundary_left_count=np.array([0.0], dtype=np.float32),
        boundary_right_count=np.array([0.0], dtype=np.float32),
        boundary_count=np.array([0.0], dtype=np.float32),
        observed_compatible_count=np.array([50.0], dtype=np.float32),
        contained_leff=np.array([100.0], dtype=np.float64),
        boundary_left_leff=np.array([0.0], dtype=np.float64),
        boundary_right_leff=np.array([0.0], dtype=np.float64),
        boundary_leff=np.array([0.0], dtype=np.float64),
        anchor_intergenic=np.array([False]),
        anchor_intron=np.array([False]),
        is_anchor=np.array([False]),
        spliced_count=np.array([0.0], dtype=np.float32),
        region_length=np.array([100.0], dtype=np.float64),
    )
    strand_channels = RegionGdnaChannelEstimate(
        contained_mean=np.array([50.0], dtype=np.float32),
        contained_upper=np.array([50.0], dtype=np.float32),
        contained_rna_lower=np.array([0.0], dtype=np.float32),
        contained_precision=np.array([1.0], dtype=np.float32),
        boundary_left_mean=np.array([0.0], dtype=np.float32),
        boundary_left_upper=np.array([0.0], dtype=np.float32),
        boundary_left_rna_lower=np.array([0.0], dtype=np.float32),
        boundary_left_precision=np.array([0.0], dtype=np.float32),
        boundary_right_mean=np.array([0.0], dtype=np.float32),
        boundary_right_upper=np.array([0.0], dtype=np.float32),
        boundary_right_rna_lower=np.array([0.0], dtype=np.float32),
        boundary_right_precision=np.array([0.0], dtype=np.float32),
        flags=np.zeros(1, dtype=np.uint16),
        kappa_d=10.0,
        p_r1_sense=0.95,
        internal_rna_lower_ci=0.95,
        contained_reliability=np.array([0.0], dtype=np.float32),
        contained_log_bayes_factor=np.zeros(1, dtype=np.float32),
        boundary_left_reliability=np.zeros(1, dtype=np.float32),
        boundary_left_log_bayes_factor=np.zeros(1, dtype=np.float32),
        boundary_right_reliability=np.zeros(1, dtype=np.float32),
        boundary_right_log_bayes_factor=np.zeros(1, dtype=np.float32),
    )

    prior_mass = build_prior_mass_deconvolution(
        observation,
        mu_gdna=np.array([50.0], dtype=np.float32),
        strand_channels=strand_channels,
    )

    np.testing.assert_allclose(prior_mass.gdna_unspliced_mean, [0.0])
    np.testing.assert_allclose(prior_mass.rna_unspliced_mean, [50.0])
    np.testing.assert_allclose(prior_mass.precision, [0.0])
