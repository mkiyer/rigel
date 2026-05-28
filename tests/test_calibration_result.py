"""Tests for the v6 ``CalibrationResult`` handoff contract."""

from __future__ import annotations

import dataclasses
import json

import numpy as np
import pytest

from rigel.calibration._diagnostics import Diagnostics
from rigel.calibration._result import CalibrationResult, build_calibration_result
from rigel.calibration.background_model import BackgroundModel
from rigel.calibration.boundary_model import BoundaryLocalPosterior
from rigel.calibration.boundary_sweep import BoundarySweepResult
from rigel.calibration.calibration_iteration import RegionCalibration
from rigel.calibration.scan_payload import FL_HIST_N_BINS, CalibrationScanPayload
from rigel.calibration.signature import N_CHANNELS, N_FL_POOLS, N_SIGNATURES
from rigel.frag_length_model import FragmentLengthModels
from rigel.splice import SpliceType


def _payload() -> CalibrationScanPayload:
    region_counts = np.zeros((1, N_CHANNELS), dtype=np.float32)
    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    return CalibrationScanPayload(
        region_counts=region_counts,
        channel_mass=region_counts.sum(axis=0, dtype=np.float64),
        signature_mass=signature_mass,
        fl_pool_mass=fl_pool_mass,
        fl_pool_total=fl_pool_mass.sum(axis=1),
        region_unspliced_support=np.zeros(1, dtype=np.uint64),
        region_spliced_support=np.zeros(1, dtype=np.uint64),
        n_observed=0,
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_excluded_strand_ambig=0,
        n_unobserved=0,
        n_unannotated_ref=0,
        n_fl_unavailable=0,
        resolver_splicing_anchor_tolerance=0,
        n_regions=1,
    )


def _scan_trained() -> FragmentLengthModels:
    models = FragmentLengthModels(max_size=200)
    models.global_model.counts[100] = 100.0
    models.global_model._total_weight = 100.0
    spliced = models.category_models[SpliceType.SPLICED_ANNOT]
    spliced.counts[100] = 50.0
    spliced._total_weight = 50.0
    return models


def _background() -> BackgroundModel:
    return BackgroundModel(
        rho_off_alpha=1.0,
        rho_off_beta=100.0,
        rho_off_mean=0.01,
        seed_mask=np.array([True], dtype=bool),
        top_t_exclusion_mask=np.array([False], dtype=bool),
        n_seed_regions=1,
        n_fragments=0.0,
        eff_length=100.0,
        fit_status="sparse",
        flags=np.zeros(1, dtype=np.uint16),
    )


def _local() -> BoundaryLocalPosterior:
    return BoundaryLocalPosterior(
        alpha_excess=np.zeros(1, dtype=np.float64),
        beta_excess=np.zeros(1, dtype=np.float64),
        mu_local=np.zeros(1, dtype=np.float32),
        upper_local=np.zeros(1, dtype=np.float32),
        flags=np.zeros(1, dtype=np.uint16),
    )


def _sweep() -> BoundarySweepResult:
    return BoundarySweepResult(
        alpha_excess=np.zeros(1, dtype=np.float64),
        beta_excess=np.zeros(1, dtype=np.float64),
        forward_alpha_excess=np.zeros(1, dtype=np.float64),
        forward_beta_excess=np.zeros(1, dtype=np.float64),
        reverse_alpha_excess=np.zeros(1, dtype=np.float64),
        reverse_beta_excess=np.zeros(1, dtype=np.float64),
        mu_sweep=np.zeros(1, dtype=np.float32),
        upper_sweep=np.zeros(1, dtype=np.float32),
        transfer_weight=np.zeros(2, dtype=np.float64),
        flags=np.zeros(1, dtype=np.uint16),
    )


def _region_calibration() -> RegionCalibration:
    from rigel.calibration.calibration_iteration import (
        METHOD_STRAND,
        BackgroundDensity,
        RegionUnsplicedMass,
    )

    return RegionCalibration(
        p_states=np.array([[1.0, 0.0]], dtype=np.float32),
        mu_gdna=np.zeros(1, dtype=np.float32),
        upper_gdna=np.zeros(1, dtype=np.float32),
        rna_lower=np.zeros(1, dtype=np.float32),
        region_unspliced_mass=RegionUnsplicedMass(
            total_mass=np.zeros(1, dtype=np.float64),
            gdna_mass=np.zeros(1, dtype=np.float64),
            rna_mass=np.zeros(1, dtype=np.float64),
            region_size_bp=np.full(1, 100.0, dtype=np.float64),
            unspliced_counts=np.zeros(1, dtype=np.uint64),
            method=np.full(1, METHOD_STRAND, dtype=np.uint8),
            precision=np.zeros(1, dtype=np.float32),
            flags=np.zeros(1, dtype=np.uint16),
        ),
        background_density=BackgroundDensity.from_bootstrap(_background()),
        A_r=np.ones(1, dtype=np.float32),
        kappa_d=None,
        n_passes=1,
        converged=True,
        flags=np.zeros(1, dtype=np.uint16),
        pass_diagnostics=({"pass_index": 0, "converged": True},),
        background_model=_background(),
        boundary_local=_local(),
        boundary_sweep=_sweep(),
    )


def test_build_calibration_result_basic() -> None:
    result = build_calibration_result(
        payload=_payload(),
        scan_trained=_scan_trained(),
        region_calibration=_region_calibration(),
        region_signature=np.zeros(1, dtype=np.uint8),
    )

    assert isinstance(result, CalibrationResult)
    assert isinstance(result.diagnostics, Diagnostics)
    assert result.region_calibration.p_states.shape == (1, 2)
    assert result.background_model.rho_off_mean == pytest.approx(0.01)
    assert not hasattr(result, "density_evidence")
    assert not hasattr(result, "region_gdna")
    assert not hasattr(result, "region_exposure")
    assert not hasattr(result, "fused_region_gdna")


def test_build_calibration_result_requires_region_calibration() -> None:
    with pytest.raises(ValueError, match="region_calibration"):
        build_calibration_result(payload=_payload(), scan_trained=_scan_trained())


def test_calibration_result_is_frozen() -> None:
    result = build_calibration_result(
        payload=_payload(),
        scan_trained=_scan_trained(),
        region_calibration=_region_calibration(),
        region_signature=np.zeros(1, dtype=np.uint8),
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.n_multi_loci = 99  # type: ignore[misc]


def test_to_summary_dict_is_json_serialisable() -> None:
    result = build_calibration_result(
        payload=_payload(),
        scan_trained=_scan_trained(),
        region_calibration=_region_calibration(),
        region_signature=np.zeros(1, dtype=np.uint8),
    )

    summary = result.to_summary_dict()
    assert summary["region_calibration"]["n_regions"] == 1
    assert summary["background_model"]["n_regions"] == 1
    assert summary["boundary_local"]["n_regions"] == 1
    assert summary["boundary_sweep"]["n_regions"] == 1
    assert summary["prior_table"] is None
    json.dumps(summary)
