"""Acceptance tests for the adaptive v5/v6 grouped prior contract."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.adaptive_prior import PRIOR_STRUCTURAL_GATED
from rigel.calibration.background_model import BackgroundModel
from rigel.calibration.calibration_iteration import (
    METHOD_STRAND,
    BackgroundDensity,
    RegionCalibration,
    RegionUnsplicedMass,
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


def _region_calibration(
    *,
    unspliced: float,
    p_state: list[float],
    gdna_unspliced: float | None = None,
    rna_unspliced: float | None = None,
) -> RegionCalibration:
    total = np.array([unspliced], dtype=np.float32)
    gdna = np.array(
        [unspliced if gdna_unspliced is None else gdna_unspliced],
        dtype=np.float32,
    )
    rna = np.array(
        [0.0 if rna_unspliced is None else rna_unspliced],
        dtype=np.float32,
    )
    return RegionCalibration(
        p_states=np.array([p_state], dtype=np.float32),
        mu_gdna=total.copy(),
        upper_gdna=total.copy(),
        rna_lower=np.zeros(1, dtype=np.float32),
        region_unspliced_mass=RegionUnsplicedMass(
            total_mass=total.astype(np.float64),
            gdna_mass=gdna.astype(np.float64),
            rna_mass=rna.astype(np.float64),
            region_size_bp=np.array([100.0], dtype=np.float64),
            unspliced_counts=np.array([max(int(round(unspliced)), 0)], dtype=np.uint64),
            method=np.array([METHOD_STRAND], dtype=np.uint8),
            precision=np.zeros(1, dtype=np.float32),
            flags=np.zeros(1, dtype=np.uint16),
        ),
        background_density=BackgroundDensity.from_bootstrap(
            BackgroundModel(
                rho_off_alpha=1.0,
                rho_off_beta=99.0,
                rho_off_mean=0.01,
                seed_mask=np.ones(1, dtype=bool),
                top_t_exclusion_mask=np.zeros(1, dtype=bool),
                n_seed_regions=1,
                n_fragments=1.0,
                eff_length=100.0,
                fit_status="ok",
                flags=np.zeros(1, dtype=np.uint16),
            )
        ),
        A_r=np.ones(1, dtype=np.float32),
        kappa_d=None,
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


def test_adaptive_prior_table_exposes_paired_mass_and_summary() -> None:
    priors = assemble_priors(
        multi_loci=[_multi_locus([0])],
        em_data=_em_data(is_spliced=[False], gdna_log_liks=[-1.0]),
        index=_index(),
        calibration=_calibration(
            _region_calibration(
                unspliced=100.0,
                p_state=[0.25, 0.75],
                gdna_unspliced=25.0,
                rna_unspliced=75.0,
            )
        ),
        em_config=EMConfig(),
    )

    assert isinstance(priors, PriorTable)
    assert priors.prior_unspliced_total[0] == pytest.approx(100.0)
    assert priors.alpha_gdna_add[0] > 0.0
    assert priors.alpha_rna_add[0] > priors.alpha_gdna_add[0]
    assert priors.prior_ess_final[0] == pytest.approx(
        priors.alpha_gdna_add[0] + priors.alpha_rna_add[0]
    )
    summary = priors.to_summary_dict()
    assert summary["name"] == "p_unexpressed_soft_gate_interim_no_exposure_gating"
    assert summary["rna_call_bias"] == pytest.approx(0.5)
    assert summary["n_loci_with_prior_mass"] == 1


def test_rna_call_bias_shifts_split_without_changing_ess() -> None:
    kwargs = dict(
        multi_loci=[_multi_locus([0])],
        em_data=_em_data(is_spliced=[False], gdna_log_liks=[-1.0]),
        index=_index(),
        calibration=_calibration(
            _region_calibration(
                unspliced=100.0,
                p_state=[0.5, 0.5],
                gdna_unspliced=50.0,
                rna_unspliced=50.0,
            )
        ),
    )

    conservative = assemble_priors(**kwargs, em_config=EMConfig(rna_call_bias=0.25))
    liberal = assemble_priors(**kwargs, em_config=EMConfig(rna_call_bias=0.75))

    assert conservative.prior_ess_final[0] == pytest.approx(liberal.prior_ess_final[0])
    assert conservative.prior_rna_share_final[0] < conservative.prior_rna_share_v5[0]
    assert liberal.prior_rna_share_final[0] > liberal.prior_rna_share_v5[0]
    assert conservative.alpha_rna_add[0] < liberal.alpha_rna_add[0]


def test_enable_gdna_helper_is_structural_diagnostic_only() -> None:
    locus = _multi_locus([0, 1])
    em_data = _em_data(is_spliced=[True, False], gdna_log_liks=[-np.inf, -2.0])
    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=_index(),
        calibration=_calibration(
            _region_calibration(
                unspliced=30.0,
                p_state=[0.0, 1.0],
                gdna_unspliced=0.0,
                rna_unspliced=30.0,
            )
        ),
        em_config=EMConfig(),
    )

    assert enable_gdna_for_multilocus(locus, em_data) is True
    assert priors.enable_gdna[0] == 1
    assert priors.alpha_gdna_add[0] == pytest.approx(0.0)
    assert priors.alpha_rna_add[0] == pytest.approx(0.0)
    assert priors.prior_region_weight[0] == pytest.approx(0.0)


def test_all_spliced_locus_reports_no_structural_gdna_candidate() -> None:
    locus = _multi_locus([0, 1])
    em_data = _em_data(is_spliced=[True, True], gdna_log_liks=[0.0, 0.0])
    priors = assemble_priors(
        multi_loci=[locus],
        em_data=em_data,
        index=_index(),
        calibration=_calibration(
            _region_calibration(unspliced=20.0, p_state=[1.0, 0.0])
        ),
        em_config=EMConfig(),
    )

    assert enable_gdna_for_multilocus(locus, em_data) is False
    assert priors.enable_gdna[0] == 0
    assert priors.prior_n_local_gdna[0] > 0.0
    assert priors.alpha_gdna_add[0] == pytest.approx(0.0)
    assert priors.alpha_rna_add[0] == pytest.approx(0.0)
    assert (priors.prior_flags[0] & PRIOR_STRUCTURAL_GATED) != 0
