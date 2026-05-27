"""Acceptance tests for the v3 grouped RNA/gDNA prior contract."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.calibration_iteration import (
    PRIOR_MASS_METHOD_DENSITY,
    PriorMassDeconvolution,
    RegionCalibration,
)
from rigel.calibration.prior import PriorTable, assemble_priors, enable_gdna_for_multilocus
from rigel.config import EMConfig
from rigel.frag_length_model import FragmentLengthModel
from rigel.locus import Locus, MultiLocus


def _delta_fl(length: int, *, max_size: int = 512) -> FragmentLengthModel:
    counts = np.zeros(max_size + 1, dtype=np.float64)
    counts[length] = 10_000.0
    return FragmentLengthModel.from_counts(counts, max_size=max_size)


def _index() -> SimpleNamespace:
    return SimpleNamespace(
        region_df=pd.DataFrame(
            {
                "region_id": np.array([0], dtype=np.int64),
                "ref_name": pd.array(["chr1"], dtype="string"),
                "start": np.array([0], dtype=np.int64),
                "end": np.array([100], dtype=np.int64),
                "signature": np.zeros(1, dtype=np.uint8),
            }
        ),
        ref_name_to_id={"chr1": 0},
        ref_lengths={"chr1": 10_000},
    )


def _region_calibration(*, gdna: float, rna: float) -> RegionCalibration:
    total = np.array([gdna + rna], dtype=np.float32)
    gdna_arr = np.array([gdna], dtype=np.float32)
    rna_arr = np.array([rna], dtype=np.float32)
    return RegionCalibration(
        p_states=np.array([[1.0, 0.0, 0.0, 0.0]], dtype=np.float32),
        mu_gdna=gdna_arr.copy(),
        upper_gdna=gdna_arr.copy(),
        rna_lower=np.zeros(1, dtype=np.float32),
        prior_mass=PriorMassDeconvolution(
            unspliced_total=total,
            gdna_unspliced_mean=gdna_arr,
            rna_unspliced_mean=rna_arr,
            method=np.array([PRIOR_MASS_METHOD_DENSITY], dtype=np.uint8),
            precision=np.zeros(1, dtype=np.float32),
            flags=np.zeros(1, dtype=np.uint16),
        ),
        A_r=np.ones(1, dtype=np.float32),
        gamma_r=np.ones(1, dtype=np.float32),
        rho_off=0.01,
        kappa_d=None,
        capture_enrichment_target=1.0,
        n_passes=1,
        converged=True,
        flags=np.zeros(1, dtype=np.uint16),
        pass_diagnostics=(),
    )


def _calibration(region_calibration: RegionCalibration) -> SimpleNamespace:
    return SimpleNamespace(
        region_calibration=region_calibration,
        fl_models=SimpleNamespace(gdna=_delta_fl(50)),
    )


def _em_data(*, is_spliced: list[bool], gdna_log_liks: list[float]) -> SimpleNamespace:
    return SimpleNamespace(
        is_spliced=np.asarray(is_spliced, dtype=bool),
        gdna_log_liks=np.asarray(gdna_log_liks, dtype=np.float64),
        genomic_midpoint=np.array([50] * len(is_spliced), dtype=np.int64),
    )


def _multi_locus(unit_indices: list[int]) -> MultiLocus:
    return MultiLocus(
        multi_locus_id=0,
        transcript_indices=np.array([0], dtype=np.int32),
        unit_indices=np.asarray(unit_indices, dtype=np.int32),
        gdna_span=100,
        loci=(Locus(ref="chr1", ref_id=0, start=0, end=100),),
    )


def test_grouped_prior_table_exposes_paired_mass_and_bounded_alpha() -> None:
    priors = assemble_priors(
        multi_loci=[_multi_locus([0])],
        em_data=_em_data(is_spliced=[False], gdna_log_liks=[-1.0]),
        index=_index(),
        calibration=_calibration(_region_calibration(gdna=5.0, rna=1000.0)),
        em_config=EMConfig(),
    )

    assert isinstance(priors, PriorTable)
    assert priors.gdna_expected_count[0] == pytest.approx(5.0)
    assert priors.rna_expected_count[0] == pytest.approx(1000.0)
    assert priors.prior_unspliced_total[0] == pytest.approx(1005.0)
    assert priors.alpha_gdna_add[0] == pytest.approx(priors.gdna_prior_count_em[0])
    assert priors.alpha_rna_add[0] > priors.alpha_gdna_add[0]
    assert priors.prior_budget[0] <= EMConfig().aggregate_prior_max_count
    assert priors.prior_mass_conservation_error[0] == pytest.approx(0.0)


def test_strength_zero_preserves_projected_mass_but_disables_em_alpha() -> None:
    priors = assemble_priors(
        multi_loci=[_multi_locus([0])],
        em_data=_em_data(is_spliced=[False], gdna_log_liks=[-1.0]),
        index=_index(),
        calibration=_calibration(_region_calibration(gdna=10.0, rna=20.0)),
        em_config=EMConfig(aggregate_prior_strength=0.0),
    )

    assert priors.gdna_expected_count[0] == pytest.approx(10.0)
    assert priors.rna_expected_count[0] == pytest.approx(20.0)
    assert priors.alpha_gdna_add[0] == pytest.approx(0.0)
    assert priors.alpha_rna_add[0] == pytest.approx(0.0)
    assert priors.gdna_prior_count_em[0] == pytest.approx(0.0)


def test_enable_gdna_helper_is_structural_diagnostic_only() -> None:
    locus = _multi_locus([0, 1])
    em_data = _em_data(is_spliced=[True, False], gdna_log_liks=[-np.inf, -2.0])
    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=_index(),
        calibration=_calibration(_region_calibration(gdna=0.0, rna=30.0)),
        em_config=EMConfig(),
    )

    assert enable_gdna_for_multilocus(locus, em_data) is True
    assert priors.enable_gdna[0] == 1
    assert priors.alpha_gdna_add[0] == pytest.approx(0.0)
    assert priors.alpha_rna_add[0] > 0.0


def test_all_spliced_locus_reports_no_structural_gdna_candidate() -> None:
    locus = _multi_locus([0, 1])
    em_data = _em_data(is_spliced=[True, True], gdna_log_liks=[0.0, 0.0])
    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=_index(),
        calibration=_calibration(_region_calibration(gdna=20.0, rna=0.0)),
        em_config=EMConfig(),
    )

    assert enable_gdna_for_multilocus(locus, em_data) is False
    assert priors.enable_gdna[0] == 0
    assert priors.gdna_expected_count[0] > 0.0
    assert priors.alpha_gdna_add[0] == pytest.approx(0.0)
    assert priors.alpha_rna_add[0] == pytest.approx(0.0)
