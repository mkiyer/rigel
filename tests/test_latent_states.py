"""Tests for Phase III four-state latent tensor helpers."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration._arrays import RegionArrays
from rigel.calibration.background_model import BackgroundModel
from rigel.calibration.boundary_sweep import BoundarySweepResult
from rigel.calibration.density_observation import DensityObservation
from rigel.calibration.fractional_evidence import transcript_strand_class
from rigel.calibration.latent_states import (
    N_STATES,
    STATE_BACKGROUND,
    STATE_EXPRESSED_CAPTURE,
    STATE_EXPRESSED_OFFTARGET,
    STATE_GDNA_ONLY_CAPTURE,
    build_logbf_capture,
    build_logbf_expression,
    build_logbf_gdna_density,
    build_logbf_strand,
    build_state_log_prior,
    build_state_tensor,
)
from rigel.calibration.signature import pack_signature
from rigel.calibration.strand_deconv import RegionGdnaChannelEstimate


def _region_arrays(signatures: list[int]) -> RegionArrays:
    signature_array = np.asarray(signatures, dtype=np.uint8)
    region_count = int(signature_array.shape[0])
    starts = np.arange(region_count, dtype=np.int64) * 100
    return RegionArrays(
        ref_id=np.zeros(region_count, dtype=np.int32),
        start=starts,
        end=starts + 100,
        signature=signature_array,
        ts_class=transcript_strand_class(signature_array),
        ref_offsets=np.array([0, region_count], dtype=np.int32),
        order=np.arange(region_count, dtype=np.int64),
        n_refs=1,
    )


def _background(region_count: int, seed_mask: list[bool] | None = None) -> BackgroundModel:
    if seed_mask is None:
        seed_mask = [False] * region_count
    return BackgroundModel(
        rho_off_alpha=1.0,
        rho_off_beta=100.0,
        rho_off_mean=0.01,
        seed_mask=np.asarray(seed_mask, dtype=bool),
        top_t_exclusion_mask=np.zeros(region_count, dtype=bool),
        n_seed_regions=int(np.count_nonzero(seed_mask)),
        n_fragments=1.0,
        eff_length=100.0,
        fit_status="ok",
        flags=np.zeros(region_count, dtype=np.uint16),
    )


def _observation(
    contained: list[float],
    leff: list[float],
    spliced: list[float],
) -> DensityObservation:
    contained_count = np.asarray(contained, dtype=np.float32)
    contained_leff = np.asarray(leff, dtype=np.float64)
    spliced_count = np.asarray(spliced, dtype=np.float32)
    zeros_count = np.zeros_like(contained_count, dtype=np.float32)
    zeros_leff = np.zeros_like(contained_leff, dtype=np.float64)
    anchor = np.ones_like(contained_count, dtype=bool)
    return DensityObservation(
        contained_count=contained_count,
        boundary_left_count=zeros_count.copy(),
        boundary_right_count=zeros_count.copy(),
        boundary_count=zeros_count.copy(),
        observed_compatible_count=contained_count.copy(),
        contained_leff=contained_leff,
        boundary_left_leff=zeros_leff.copy(),
        boundary_right_leff=zeros_leff.copy(),
        boundary_leff=zeros_leff.copy(),
        anchor_intergenic=anchor.copy(),
        anchor_intron=np.zeros_like(anchor, dtype=bool),
        is_anchor=anchor,
        spliced_count=spliced_count,
        region_length=contained_leff.copy(),
    )


def _sweep(mu: list[float], alpha: list[float], beta: list[float]) -> BoundarySweepResult:
    mu_array = np.asarray(mu, dtype=np.float32)
    region_count = int(mu_array.shape[0])
    alpha_array = np.asarray(alpha, dtype=np.float64)
    beta_array = np.asarray(beta, dtype=np.float64)
    return BoundarySweepResult(
        alpha_excess=alpha_array,
        beta_excess=beta_array,
        forward_alpha_excess=np.zeros(region_count, dtype=np.float64),
        forward_beta_excess=np.zeros(region_count, dtype=np.float64),
        reverse_alpha_excess=np.zeros(region_count, dtype=np.float64),
        reverse_beta_excess=np.zeros(region_count, dtype=np.float64),
        mu_sweep=mu_array,
        upper_sweep=np.maximum(mu_array, 0.0).astype(np.float32),
        transfer_weight=np.zeros(region_count + 1, dtype=np.float64),
        flags=np.zeros(region_count, dtype=np.uint16),
    )


def _strand_channels(rna_lower: list[float]) -> RegionGdnaChannelEstimate:
    rna = np.asarray(rna_lower, dtype=np.float32)
    region_count = int(rna.shape[0])
    zeros = np.zeros(region_count, dtype=np.float32)
    return RegionGdnaChannelEstimate(
        contained_mean=zeros.copy(),
        contained_upper=zeros.copy(),
        contained_rna_lower=rna,
        contained_precision=np.ones(region_count, dtype=np.float32),
        boundary_left_mean=zeros.copy(),
        boundary_left_upper=zeros.copy(),
        boundary_left_rna_lower=zeros.copy(),
        boundary_left_precision=zeros.copy(),
        boundary_right_mean=zeros.copy(),
        boundary_right_upper=zeros.copy(),
        boundary_right_rna_lower=zeros.copy(),
        boundary_right_precision=zeros.copy(),
        flags=np.zeros(region_count, dtype=np.uint16),
        kappa_d=10.0,
        p_r1_sense=0.95,
        internal_rna_lower_ci=0.95,
    )


def test_state_log_prior_favors_background_for_intergenic_and_intronic_regions() -> None:
    region_arrays = _region_arrays(
        [pack_signature(), pack_signature(intron_pos=True), pack_signature(exon_pos=True)]
    )
    background = _background(3, seed_mask=[True, False, False])

    prior = build_state_log_prior(region_arrays, background, pass_index=0)

    assert prior.shape == (3, N_STATES)
    assert prior[0, STATE_BACKGROUND] == pytest.approx(2.0)
    assert prior[0, STATE_GDNA_ONLY_CAPTURE] == pytest.approx(0.0)
    assert prior[0, STATE_EXPRESSED_CAPTURE] == pytest.approx(-2.0)
    assert prior[0, STATE_EXPRESSED_OFFTARGET] == pytest.approx(-2.0)
    assert prior[1, STATE_BACKGROUND] == pytest.approx(0.75)
    assert prior[1, STATE_EXPRESSED_OFFTARGET] == pytest.approx(-1.0)
    np.testing.assert_allclose(prior[2], 0.0)

    late_prior = build_state_log_prior(region_arrays, background, pass_index=2)
    assert late_prior[0, STATE_BACKGROUND] == pytest.approx(1.0)


def test_likelihood_terms_have_expected_direction() -> None:
    observation = _observation([0.0, 30.0], [1000.0, 1000.0], [0.0, 5.0])
    background = _background(2)
    sweep = _sweep([10.0, 30.0], [0.0, 5.0], [0.0, 100.0])

    expression = build_logbf_expression(observation)
    capture = build_logbf_capture(sweep, observation, background)
    density = build_logbf_gdna_density(observation, background, sweep)

    assert expression[1, 0] < expression[1, 1]
    assert capture[0, 1] < capture[0, 0]
    assert capture[1, 1] > capture[1, 0]
    assert density[1, 1] > density[1, 0]


def test_strand_summary_penalizes_not_expressed_states_when_rna_lower_bound_exists() -> None:
    logbf = build_logbf_strand(_strand_channels([0.0, 5.0]), 2)

    assert logbf[0, STATE_BACKGROUND] == pytest.approx(0.0)
    assert logbf[1, STATE_BACKGROUND] < 0.0
    assert logbf[1, STATE_GDNA_ONLY_CAPTURE] < 0.0
    assert logbf[1, STATE_EXPRESSED_CAPTURE] == pytest.approx(0.0)
    assert logbf[1, STATE_EXPRESSED_OFFTARGET] == pytest.approx(0.0)


def test_state_tensor_rows_sum_to_one_and_identify_combined_evidence() -> None:
    state_log_prior = np.zeros((2, N_STATES), dtype=np.float64)
    expression = np.array([[0.0, 0.0], [-10.0, 0.0]], dtype=np.float64)
    capture = np.array([[0.0, -10.0], [-10.0, 0.0]], dtype=np.float64)
    density = np.zeros((2, 2), dtype=np.float64)
    strand = np.zeros((2, N_STATES), dtype=np.float64)

    p_states = build_state_tensor(state_log_prior, expression, capture, density, strand)

    np.testing.assert_allclose(p_states.sum(axis=1), 1.0, rtol=1.0e-6, atol=1.0e-6)
    assert int(np.argmax(p_states[1])) == STATE_EXPRESSED_CAPTURE
