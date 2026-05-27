"""Tests for v6 boundary-to-contained gDNA imputation."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.background_model import BackgroundModel
from rigel.calibration.boundaries import BoundaryTable
from rigel.calibration.boundary_model import (
    FLAG_BOUNDARY_NO_EVIDENCE,
    FLAG_BOUNDARY_STRAND_RNA,
    build_boundary_local_posterior,
    predict_contained_gdna_from_excess,
)
from rigel.calibration.density_observation import DensityObservation
from rigel.calibration.strand_deconv import RegionGdnaChannelEstimate


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


def _observation(leff: list[float]) -> DensityObservation:
    contained_leff = np.asarray(leff, dtype=np.float64)
    count = np.zeros(contained_leff.shape, dtype=np.float32)
    zero_leff = np.zeros_like(contained_leff, dtype=np.float64)
    mask = np.ones(contained_leff.shape, dtype=bool)
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
    *,
    left_unspliced: list[float] | None = None,
    right_unspliced: list[float] | None = None,
    left_leff: list[float] | None = None,
    right_leff: list[float] | None = None,
) -> BoundaryTable:
    boundary_count = 4
    zeros32 = np.zeros(boundary_count, dtype=np.float32)
    zeros64 = np.zeros(boundary_count, dtype=np.float64)
    left_mass = zeros32.copy() if left_unspliced is None else np.asarray(left_unspliced, dtype=np.float32)
    right_mass = (
        zeros32.copy() if right_unspliced is None else np.asarray(right_unspliced, dtype=np.float32)
    )
    left_opp = zeros64.copy() if left_leff is None else np.asarray(left_leff, dtype=np.float64)
    right_opp = zeros64.copy() if right_leff is None else np.asarray(right_leff, dtype=np.float64)
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
        right_region_unspliced_pos=right_mass,
        right_region_unspliced_neg=zeros32.copy(),
        right_region_spliced_pos=zeros32.copy(),
        right_region_spliced_neg=zeros32.copy(),
        left_region_boundary_leff=left_opp,
        right_region_boundary_leff=right_opp,
    )


def _strand_channels(boundary_right_mean: list[float]) -> RegionGdnaChannelEstimate:
    mean = np.zeros(3, dtype=np.float32)
    boundary_right = np.asarray(boundary_right_mean, dtype=np.float32)
    zeros = np.zeros(3, dtype=np.float32)
    return RegionGdnaChannelEstimate(
        contained_mean=mean.copy(),
        contained_upper=mean.copy(),
        contained_rna_lower=zeros.copy(),
        contained_precision=zeros.copy(),
        boundary_left_mean=zeros.copy(),
        boundary_left_upper=zeros.copy(),
        boundary_left_rna_lower=zeros.copy(),
        boundary_left_precision=zeros.copy(),
        boundary_right_mean=boundary_right,
        boundary_right_upper=boundary_right.copy(),
        boundary_right_rna_lower=np.array([0.0, 2.0, 0.0], dtype=np.float32),
        boundary_right_precision=zeros.copy(),
        flags=np.zeros(3, dtype=np.uint16),
        kappa_d=10.0,
        p_r1_sense=0.95,
        internal_rna_lower_ci=0.95,
    )


def test_predict_contained_gdna_uses_background_when_no_excess() -> None:
    mean, upper = predict_contained_gdna_from_excess(
        _background(),
        np.array([1000.0], dtype=np.float64),
        np.array([0.0], dtype=np.float64),
        np.array([0.0], dtype=np.float64),
        confidence=0.95,
    )

    assert mean[0] == pytest.approx(10.0)
    assert upper[0] >= mean[0]


def test_local_boundary_prediction_increases_with_strong_boundary_evidence() -> None:
    boundaries = _boundaries(
        left_unspliced=[0.0, 50.0, 0.0, 0.0],
        left_leff=[0.0, 100.0, 0.0, 0.0],
    )
    posterior = build_boundary_local_posterior(_observation([1000.0, 1000.0, 1000.0]), boundaries, _background())

    assert posterior.alpha_excess.tolist() == [50.0, 0.0, 0.0]
    assert posterior.beta_excess.tolist() == [100.0, 0.0, 0.0]
    assert posterior.mu_local[0] == pytest.approx(255.0)
    assert posterior.mu_local[1] == pytest.approx(10.0)
    assert (posterior.flags[1] & FLAG_BOUNDARY_NO_EVIDENCE) != 0


def test_local_boundary_prediction_uses_strand_deconvolved_boundary_counts() -> None:
    boundaries = _boundaries(
        left_unspliced=[0.0, 100.0, 0.0, 0.0],
        left_leff=[0.0, 100.0, 0.0, 0.0],
    )
    channels = _strand_channels([5.0, 0.0, 0.0])

    posterior = build_boundary_local_posterior(
        _observation([1000.0, 1000.0, 1000.0]),
        boundaries,
        _background(),
        strand_channels=channels,
    )

    assert posterior.alpha_excess[0] == pytest.approx(5.0)
    assert posterior.mu_local[0] == pytest.approx(30.0)
    assert (posterior.flags[1] & FLAG_BOUNDARY_STRAND_RNA) != 0