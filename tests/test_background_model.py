"""Tests for the v6 off-target background model."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.background_model import (
    FLAG_BACKGROUND_SEED,
    FLAG_STRAND_RNA_EXCLUDED,
    FLAG_TOP_T_EXCLUDED,
    fit_background_model,
)
from rigel.calibration.density_observation import DensityObservation
from rigel.calibration.strand_deconv import RegionGdnaChannelEstimate


def _observation(
    contained: list[float],
    leff: list[float],
    *,
    anchor: list[bool] | None = None,
    spliced: list[float] | None = None,
) -> DensityObservation:
    counts = np.asarray(contained, dtype=np.float32)
    contained_leff = np.asarray(leff, dtype=np.float64)
    if anchor is None:
        anchor = [True] * counts.size
    if spliced is None:
        spliced = [0.0] * counts.size
    anchor_mask = np.asarray(anchor, dtype=bool)
    zeros_count = np.zeros_like(counts, dtype=np.float32)
    zeros_leff = np.zeros_like(contained_leff, dtype=np.float64)
    return DensityObservation(
        contained_count=counts,
        boundary_left_count=zeros_count.copy(),
        boundary_right_count=zeros_count.copy(),
        boundary_count=zeros_count.copy(),
        observed_compatible_count=counts.copy(),
        contained_leff=contained_leff,
        boundary_left_leff=zeros_leff.copy(),
        boundary_right_leff=zeros_leff.copy(),
        boundary_leff=zeros_leff.copy(),
        anchor_intergenic=anchor_mask.copy(),
        anchor_intron=np.zeros_like(anchor_mask, dtype=bool),
        is_anchor=anchor_mask,
        spliced_count=np.asarray(spliced, dtype=np.float32),
        region_length=contained_leff.copy(),
    )


def _strand_channels(
    contained_mean: list[float],
    contained_rna_lower: list[float],
) -> RegionGdnaChannelEstimate:
    mean = np.asarray(contained_mean, dtype=np.float32)
    rna_lower = np.asarray(contained_rna_lower, dtype=np.float32)
    zeros = np.zeros_like(mean, dtype=np.float32)
    return RegionGdnaChannelEstimate(
        contained_mean=mean,
        contained_upper=mean.copy(),
        contained_rna_lower=rna_lower,
        contained_precision=zeros.copy(),
        boundary_left_mean=zeros.copy(),
        boundary_left_upper=zeros.copy(),
        boundary_left_rna_lower=zeros.copy(),
        boundary_left_precision=zeros.copy(),
        boundary_right_mean=zeros.copy(),
        boundary_right_upper=zeros.copy(),
        boundary_right_rna_lower=zeros.copy(),
        boundary_right_precision=zeros.copy(),
        flags=np.zeros(mean.shape, dtype=np.uint16),
        kappa_d=10.0,
        p_r1_sense=0.95,
        rna_lower_confidence=0.95,
    )


def test_background_model_fits_gamma_posterior_from_anchor_regions() -> None:
    obs = _observation([10.0, 20.0, 40.0, 100.0], [1000.0, 2000.0, 4000.0, 1000.0])

    model = fit_background_model(
        obs,
        top_t_fraction=0.0,
        min_seed_regions=3,
        alpha_floor=1.0,
        beta_floor=1.0,
    )

    assert model.fit_status == "ok"
    assert model.seed_mask.tolist() == [True, True, True, True]
    assert model.n_seed_regions == 4
    assert model.rho_off_alpha == pytest.approx(171.0)
    assert model.rho_off_beta == pytest.approx(8001.0)
    assert model.rho_off_mean == pytest.approx(171.0 / 8001.0)
    assert np.all((model.flags[model.seed_mask] & FLAG_BACKGROUND_SEED) != 0)


def test_background_model_excludes_spliced_nonanchors_low_opportunity_and_top_tail() -> None:
    obs = _observation(
        [10.0, 12.0, 500.0, 8.0, 9.0],
        [1000.0, 1000.0, 1000.0, 1000.0, 0.25],
        anchor=[True, True, True, False, True],
        spliced=[0.0, 3.0, 0.0, 0.0, 0.0],
    )

    model = fit_background_model(obs, top_t_fraction=0.5, min_seed_regions=1)

    assert model.seed_mask.tolist() == [True, False, False, False, False]
    assert model.top_t_exclusion_mask.tolist() == [False, False, True, False, False]
    assert (model.flags[2] & FLAG_TOP_T_EXCLUDED) != 0


def test_background_model_uses_strand_channel_gdna_and_excludes_rna_lower_bound() -> None:
    obs = _observation([100.0, 100.0, 100.0], [1000.0, 1000.0, 1000.0])
    channels = _strand_channels([20.0, 30.0, 40.0], [0.0, 5.0, 0.0])

    model = fit_background_model(
        obs,
        strand_channels=channels,
        top_t_fraction=0.0,
        min_seed_regions=2,
        alpha_floor=1.0,
        beta_floor=1.0,
    )

    assert model.seed_mask.tolist() == [True, False, True]
    assert model.fit_status == "ok"
    assert model.n_fragments == pytest.approx(60.0)
    assert model.rho_off_mean == pytest.approx(61.0 / 2001.0)
    assert (model.flags[1] & FLAG_STRAND_RNA_EXCLUDED) != 0