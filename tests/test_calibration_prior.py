"""Tests for adaptive prior projection onto MultiLocus EM inputs."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.adaptive_prior import PRIOR_STRUCTURAL_GATED
from rigel.calibration.calibration_iteration import (
    PRIOR_MASS_METHOD_DENSITY,
    PriorMassDeconvolution,
    RegionCalibration,
)
from rigel.calibration.prior import assemble_priors, enable_gdna_for_multilocus
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


def _region_calibration(
    unspliced_total: list[float],
    *,
    p_states: np.ndarray | None = None,
    exposure: list[float] | None = None,
) -> RegionCalibration:
    unspliced = np.asarray(unspliced_total, dtype=np.float32)
    region_count = int(unspliced.size)
    if p_states is None:
        p_states = np.zeros((region_count, 4), dtype=np.float32)
        p_states[:, 0] = 1.0
    exposure_arr = np.asarray(
        exposure if exposure is not None else np.ones(region_count, dtype=np.float32),
        dtype=np.float32,
    )
    return RegionCalibration(
        p_states=np.asarray(p_states, dtype=np.float32),
        mu_gdna=unspliced.copy(),
        upper_gdna=unspliced.copy(),
        rna_lower=np.zeros(region_count, dtype=np.float32),
        prior_mass=PriorMassDeconvolution(
            unspliced_total=unspliced,
            gdna_unspliced_mean=unspliced.copy(),
            rna_unspliced_mean=np.zeros(region_count, dtype=np.float32),
            method=np.full(region_count, PRIOR_MASS_METHOD_DENSITY, dtype=np.uint8),
            precision=np.zeros(region_count, dtype=np.float32),
            flags=np.zeros(region_count, dtype=np.uint16),
        ),
        A_r=exposure_arr,
        gamma_r=exposure_arr.copy(),
        rho_off=0.01,
        kappa_d=None,
        capture_enrichment_target=1.0,
        n_passes=1,
        converged=True,
        flags=np.zeros(region_count, dtype=np.uint16),
        pass_diagnostics=(),
    )


def _calibration(region_calibration: RegionCalibration) -> SimpleNamespace:
    return SimpleNamespace(
        region_calibration=region_calibration,
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
        calibration=_calibration(_region_calibration([10.0, 20.0])),
    )

    np.testing.assert_allclose(priors.prior_n_local_gdna, np.array([10.0, 20.0]))
    np.testing.assert_allclose(priors.alpha_gdna_add, np.array([10.0, 20.0]))
    np.testing.assert_allclose(priors.alpha_rna_add, np.zeros(2))
    assert priors.unallocated_unspliced_count == pytest.approx(0.0)
    assert priors.global_n_gdna == pytest.approx(30.0)
    np.testing.assert_array_equal(priors.n_regions_touched, np.array([1, 1], dtype=np.int32))


def test_geometry_allocation_reports_regions_touching_multiple_loci() -> None:
    index = _index([(0, 100)])
    loci = [_ml(0, 0, 60, [0]), _ml(1, 40, 100, [1])]
    em_data = _em_data(is_spliced=[False, False], gdna_log_liks=[-1.0, -2.0])

    priors = assemble_priors(
        multi_loci=loci,
        em_data=em_data,
        index=index,
        calibration=_calibration(_region_calibration([10.0])),
    )

    np.testing.assert_allclose(priors.prior_n_local_gdna, np.array([5.0, 5.0]))
    np.testing.assert_allclose(priors.multi_locus_region_mass, np.array([5.0, 5.0]))
    assert np.sum(priors.alpha_gdna_add) == pytest.approx(10.0)


def test_geometry_allocation_reports_partial_coverage_mass() -> None:
    index = _index([(0, 100)])
    locus = _ml(0, 0, 50, [0])
    em_data = _em_data(is_spliced=[False], gdna_log_liks=[-1.0])

    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=index,
        calibration=_calibration(_region_calibration([10.0])),
    )

    assert priors.prior_n_local_gdna[0] == pytest.approx(10.0)
    assert priors.partial_coverage_region_mass[0] == pytest.approx(10.0)


def test_zero_unspliced_mass_keeps_native_gdna_eligible_with_candidate() -> None:
    index = _index([(0, 100)])
    locus = _ml(0, 0, 100, [0, 1])
    em_data = _em_data(is_spliced=[True, False], gdna_log_liks=[-np.inf, -1.0])

    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=index,
        calibration=_calibration(_region_calibration([0.0])),
    )

    assert priors.alpha_gdna_add[0] == pytest.approx(0.0)
    assert priors.alpha_rna_add[0] == pytest.approx(0.0)
    assert priors.enable_gdna[0] == 1
    assert enable_gdna_for_multilocus(locus, em_data) is True


def test_all_spliced_locus_is_structurally_gated_even_with_positive_prior() -> None:
    index = _index([(0, 100)])
    locus = _ml(0, 0, 100, [0, 1])
    em_data = _em_data(is_spliced=[True, True], gdna_log_liks=[-1.0, -2.0])

    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=index,
        calibration=_calibration(_region_calibration([10.0])),
    )

    assert priors.prior_n_local_gdna[0] > 0.0
    assert priors.alpha_gdna_add[0] == pytest.approx(0.0)
    assert priors.alpha_rna_add[0] == pytest.approx(0.0)
    assert (priors.prior_flags[0] & PRIOR_STRUCTURAL_GATED) != 0
    assert priors.enable_gdna[0] == 0


def test_uniform_exposure_keeps_weighted_denominator_equal_to_unweighted() -> None:
    index = _index([(0, 100)])
    locus = _ml(0, 0, 100, [0])
    em_data = _em_data(is_spliced=[False], gdna_log_liks=[-1.0])

    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=index,
        calibration=_calibration(_region_calibration([10.0])),
    )

    assert priors.gdna_eff_len[0] == pytest.approx(priors.gdna_eff_len_unweighted[0])
    assert priors.gdna_em_exposure_weight[0] == pytest.approx(1.0)


def test_region_calibration_exposure_weights_denominator() -> None:
    index = _index([(0, 100)])
    locus = _ml(0, 0, 100, [0])
    em_data = _em_data(is_spliced=[False], gdna_log_liks=[-1.0])

    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=index,
        calibration=_calibration(_region_calibration([10.0], exposure=[2.5])),
    )

    assert priors.gdna_eff_len[0] == pytest.approx(
        priors.gdna_eff_len_unweighted[0] * 2.5
    )
    assert priors.gdna_em_exposure_weight[0] == pytest.approx(2.5)
