"""Tests for fused regional prior allocation onto MultiLocus EM inputs."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.exposure import RegionExposure
from rigel.calibration.integration import FusedRegionGdnaEvidence
from rigel.calibration.prior import (
    COUNT_ALLOC_GEOMETRY,
    assemble_priors,
    enable_gdna_for_multilocus,
)
from rigel.frag_length_model import FragmentLengthModel
from rigel.locus import Locus, MultiLocus


def _delta_fl(length: int, *, max_size: int = 512) -> FragmentLengthModel:
    counts = np.zeros(max_size + 1, dtype=np.float64)
    counts[length] = 10_000.0
    return FragmentLengthModel.from_counts(counts, max_size=max_size)


def _index(region_rows: list[tuple[int, int]]) -> SimpleNamespace:
    return SimpleNamespace(
        region_df=pd.DataFrame(
            {
                "region_id": np.arange(len(region_rows), dtype=np.int64),
                "ref_name": pd.array(["chr1"] * len(region_rows), dtype="string"),
                "start": np.array([start for start, _end in region_rows], dtype=np.int64),
                "end": np.array([end for _start, end in region_rows], dtype=np.int64),
                "signature": np.zeros(len(region_rows), dtype=np.uint8),
            }
        ),
        ref_name_to_id={"chr1": 0},
        ref_lengths={"chr1": 10_000},
    )


def _fused(mean: list[float], upper: list[float] | None = None) -> FusedRegionGdnaEvidence:
    mean_arr = np.asarray(mean, dtype=np.float32)
    upper_arr = np.asarray(upper if upper is not None else mean, dtype=np.float32)
    zeros = np.zeros(mean_arr.size, dtype=np.float32)
    return FusedRegionGdnaEvidence(
        mean_count=mean_arr,
        upper_count=upper_arr,
        variance_count=zeros.copy(),
        rna_lower_count=zeros.copy(),
        observed_compatible_count=mean_arr.copy(),
        density_weight=np.ones(mean_arr.size, dtype=np.float32),
        strand_weight=zeros.copy(),
        density_applicable=np.ones(mean_arr.size, dtype=bool),
        strand_applicable=np.zeros(mean_arr.size, dtype=bool),
        tail_probability=zeros.copy(),
        expected_tail_count=zeros.copy(),
        flags=np.zeros(mean_arr.size, dtype=np.uint8),
    )


def _calibration(fused: FusedRegionGdnaEvidence) -> SimpleNamespace:
    return SimpleNamespace(
        fused_region_gdna=fused,
        region_exposure=RegionExposure.uniform(int(fused.mean_count.size)),
        fl_models=SimpleNamespace(gdna=_delta_fl(50)),
    )


def _em_data(
    *,
    is_spliced: list[bool],
    gdna_log_liks: list[float],
    genomic_midpoint: list[int] | None = None,
) -> SimpleNamespace:
    if genomic_midpoint is None:
        genomic_midpoint = [10] * len(is_spliced)
    return SimpleNamespace(
        is_spliced=np.asarray(is_spliced, dtype=bool),
        gdna_log_liks=np.asarray(gdna_log_liks, dtype=np.float64),
        genomic_midpoint=np.asarray(genomic_midpoint, dtype=np.int64),
    )


def _ml(lid: int, start: int, end: int, units: list[int]) -> MultiLocus:
    locus = Locus(ref="chr1", ref_id=0, start=start, end=end)
    return MultiLocus(
        multi_locus_id=lid,
        transcript_indices=np.array([lid], dtype=np.int32),
        unit_indices=np.asarray(units, dtype=np.int32),
        gdna_span=max(end - start, 1),
        loci=(locus,),
    )


def test_geometry_allocation_conserves_mass_for_disjoint_loci() -> None:
    index = _index([(0, 100), (100, 200)])
    loci = [_ml(0, 0, 100, [0]), _ml(1, 100, 200, [1])]
    em_data = _em_data(is_spliced=[False, False], gdna_log_liks=[-1.0, -2.0])

    priors = assemble_priors(
        multi_loci=loci,
        em_data=em_data,
        index=index,
        calibration=_calibration(_fused([10.0, 20.0], [12.0, 24.0])),
    )

    np.testing.assert_allclose(priors.gdna_expected_count, np.array([10.0, 20.0]))
    np.testing.assert_allclose(priors.gdna_upper_count, np.array([12.0, 24.0]))
    assert priors.unallocated_expected_count == pytest.approx(0.0)
    assert np.sum(priors.gdna_prior_count_em) == pytest.approx(30.0)
    np.testing.assert_array_equal(priors.count_allocation_mode, np.full(2, COUNT_ALLOC_GEOMETRY))


def test_geometry_allocation_reports_regions_touching_multiple_loci() -> None:
    index = _index([(0, 100)])
    loci = [_ml(0, 0, 60, [0]), _ml(1, 40, 100, [1])]
    em_data = _em_data(is_spliced=[False, False], gdna_log_liks=[-1.0, -2.0])

    priors = assemble_priors(
        multi_loci=loci,
        em_data=em_data,
        index=index,
        calibration=_calibration(_fused([10.0])),
    )

    np.testing.assert_allclose(priors.gdna_expected_count, np.array([5.0, 5.0]))
    np.testing.assert_allclose(priors.multi_locus_region_mass, np.array([5.0, 5.0]))
    assert np.sum(priors.gdna_expected_count) == pytest.approx(10.0)


def test_geometry_allocation_reports_partial_coverage_mass() -> None:
    index = _index([(0, 100)])
    loci = [_ml(0, 0, 50, [0])]
    em_data = _em_data(is_spliced=[False], gdna_log_liks=[-1.0])

    priors = assemble_priors(
        multi_loci=loci,
        em_data=em_data,
        index=index,
        calibration=_calibration(_fused([10.0])),
    )

    assert priors.gdna_expected_count[0] == pytest.approx(10.0)
    assert priors.partial_coverage_region_mass[0] == pytest.approx(10.0)


def test_zero_expected_prior_with_positive_upper_still_enables_gdna() -> None:
    index = _index([(0, 100)])
    locus = _ml(0, 0, 100, [0, 1])
    em_data = _em_data(is_spliced=[True, False], gdna_log_liks=[-np.inf, -1.0])

    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=index,
        calibration=_calibration(_fused([0.0], [2.0])),
    )

    assert priors.gdna_prior_count_em[0] == pytest.approx(0.0)
    assert priors.enable_gdna[0] == 1
    assert enable_gdna_for_multilocus(locus, em_data) is True


def test_zero_upper_bound_disables_gdna_even_with_candidate() -> None:
    index = _index([(0, 100)])
    locus = _ml(0, 0, 100, [0, 1])
    em_data = _em_data(is_spliced=[True, False], gdna_log_liks=[-np.inf, -1.0])

    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=index,
        calibration=_calibration(_fused([0.0], [0.0])),
    )

    assert priors.gdna_prior_count_em[0] == pytest.approx(0.0)
    assert priors.gdna_upper_count[0] == pytest.approx(0.0)
    assert priors.enable_gdna[0] == 0
    assert enable_gdna_for_multilocus(locus, em_data) is True


def test_all_spliced_locus_is_ineligible_even_with_positive_prior() -> None:
    index = _index([(0, 100)])
    locus = _ml(0, 0, 100, [0, 1])
    em_data = _em_data(is_spliced=[True, True], gdna_log_liks=[-1.0, -2.0])

    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=index,
        calibration=_calibration(_fused([10.0])),
    )

    assert priors.gdna_prior_count_em[0] > 0.0
    assert priors.enable_gdna[0] == 0


def test_uniform_exposure_keeps_weighted_denominator_equal_to_unweighted() -> None:
    index = _index([(0, 100)])
    locus = _ml(0, 0, 100, [0])
    em_data = _em_data(is_spliced=[False], gdna_log_liks=[-1.0])

    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=index,
        calibration=_calibration(_fused([10.0])),
    )

    assert priors.gdna_eff_len[0] == pytest.approx(priors.gdna_eff_len_unweighted[0])
    assert priors.gdna_em_exposure_weight[0] == pytest.approx(1.0)
