"""Current per-MultiLocus adaptive gDNA prior projection tests.

The retired per-locus gDNA mass estimator was replaced by grouped adaptive
priors assembled from region calibration state probabilities and projected
onto ``MultiLocus`` geometry.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.calibration_iteration import PriorMassDeconvolution, RegionCalibration
from rigel.calibration.latent_states import N_STATES, STATE_IS_EXPRESSED
from rigel.calibration.prior import assemble_priors, enable_gdna_for_multilocus
from rigel.calibration.signature import pack_signature
from rigel.config import EMConfig
from rigel.frag_length_model import FragmentLengthModel
from rigel.locus import Locus, MultiLocus
from rigel.scored_fragments import ScoredFragments


GDNA_STATE = [1.0, 0.0, 0.0, 0.0]
RNA_STATE = [0.0, 0.0, 1.0, 0.0]
UNIFORM_STATE = [0.25, 0.25, 0.25, 0.25]


def _region_df() -> pd.DataFrame:
    rows = [
        ("chr1", 0, 100, pack_signature()),
        ("chr1", 100, 200, pack_signature(intron_pos=True)),
        ("chr1", 200, 300, pack_signature(exon_pos=True)),
    ]
    df = pd.DataFrame(rows, columns=["ref_name", "start", "end", "signature"])
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    df.index = df["region_id"].to_numpy()
    return df


def _delta_fl(length: int = 50) -> FragmentLengthModel:
    counts = np.zeros(128, dtype=np.float64)
    counts[length] = 1000.0
    return FragmentLengthModel.from_counts(counts, max_size=127)


def _region_calibration(
    p_states: list[list[float]],
    unspliced_total: list[float],
    *,
    exposure: list[float] | None = None,
) -> RegionCalibration:
    p = np.asarray(p_states, dtype=np.float32)
    unspliced = np.asarray(unspliced_total, dtype=np.float32)
    if exposure is None:
        exposure = [1.0] * len(unspliced)
    if p.ndim == 2 and p.shape[1] == N_STATES and p.shape[0] == len(unspliced):
        gdna_share = p[:, ~STATE_IS_EXPRESSED].sum(axis=1)
        rna_share = p[:, STATE_IS_EXPRESSED].sum(axis=1)
        total = gdna_share + rna_share
        gdna_mean = unspliced * np.divide(gdna_share, total, out=np.zeros_like(total), where=total > 0.0)
        rna_mean = unspliced - gdna_mean
    else:
        gdna_mean = np.zeros(len(unspliced), dtype=np.float32)
        rna_mean = unspliced.copy()
    prior_mass = PriorMassDeconvolution(
        unspliced_total=unspliced,
        gdna_unspliced_mean=gdna_mean,
        rna_unspliced_mean=rna_mean,
        method=np.ones(len(unspliced), dtype=np.uint8),
        precision=np.ones(len(unspliced), dtype=np.float32),
        flags=np.zeros(len(unspliced), dtype=np.uint16),
    )
    return RegionCalibration(
        p_states=p,
        mu_gdna=prior_mass.gdna_unspliced_mean,
        upper_gdna=prior_mass.gdna_unspliced_mean,
        rna_lower=prior_mass.rna_unspliced_mean,
        prior_mass=prior_mass,
        A_r=np.asarray(exposure, dtype=np.float32),
        gamma_r=np.ones(len(unspliced), dtype=np.float32),
        rho_off=0.0,
        kappa_d=None,
        capture_enrichment_target=1.0,
        n_passes=1,
        converged=True,
        flags=np.zeros(len(unspliced), dtype=np.uint16),
        pass_diagnostics=(),
    )


def _multi_locus(locus_id: int, start: int, end: int, *, unit_indices=None) -> MultiLocus:
    if unit_indices is None:
        unit_indices = [locus_id]
    return MultiLocus(
        multi_locus_id=locus_id,
        transcript_indices=np.array([locus_id], dtype=np.int32),
        unit_indices=np.asarray(unit_indices, dtype=np.int32),
        gdna_span=max(end - start, 1),
        loci=(Locus(ref="chr1", ref_id=0, start=start, end=end),),
    )


def _scored_fragments(*, gdna_log_liks: list[float], is_spliced: list[bool]) -> ScoredFragments:
    n_units = len(gdna_log_liks)
    return ScoredFragments(
        offsets=np.arange(n_units + 1, dtype=np.int64),
        t_indices=np.arange(n_units, dtype=np.int32),
        log_liks=np.zeros(n_units, dtype=np.float64),
        count_cols=np.zeros(n_units, dtype=np.uint8),
        coverage_weights=np.ones(n_units, dtype=np.float64),
        locus_t_indices=np.arange(n_units, dtype=np.int32),
        locus_count_cols=np.zeros(n_units, dtype=np.uint8),
        is_spliced=np.asarray(is_spliced, dtype=bool),
        gdna_log_liks=np.asarray(gdna_log_liks, dtype=np.float64),
        frag_ids=np.arange(n_units, dtype=np.int64),
        frag_class=np.zeros(n_units, dtype=np.int8),
        splice_type=np.zeros(n_units, dtype=np.uint8),
        genomic_midpoint=np.arange(n_units, dtype=np.int64),
        n_units=n_units,
        n_candidates=n_units,
    )


def _index_and_calibration(region_calibration: RegionCalibration):
    index = SimpleNamespace(
        region_df=_region_df(),
        ref_name_to_id={"chr1": 0},
        ref_lengths={"chr1": 1000, 0: 1000},
    )
    calibration = SimpleNamespace(
        region_calibration=region_calibration,
        fl_models=SimpleNamespace(gdna=_delta_fl()),
    )
    return index, calibration


def test_enable_gdna_requires_unspliced_unit_with_finite_gdna_likelihood() -> None:
    locus = _multi_locus(0, 0, 100, unit_indices=[0, 1, 2])
    em_data = _scored_fragments(
        gdna_log_liks=[-1.0, np.inf, -3.0],
        is_spliced=[True, False, False],
    )

    assert enable_gdna_for_multilocus(locus, em_data)

    no_candidate = _scored_fragments(
        gdna_log_liks=[-1.0, np.inf, np.nan],
        is_spliced=[True, False, False],
    )
    assert not enable_gdna_for_multilocus(locus, no_candidate)


def test_assemble_priors_projects_region_state_mass_to_matching_loci() -> None:
    region_calibration = _region_calibration(
        [GDNA_STATE, RNA_STATE, UNIFORM_STATE],
        [10.0, 20.0, 30.0],
    )
    index, calibration = _index_and_calibration(region_calibration)
    multi_loci = [_multi_locus(0, 0, 100), _multi_locus(1, 100, 200), _multi_locus(2, 200, 300)]
    em_data = _scored_fragments(
        gdna_log_liks=[-1.0, -1.0, -1.0],
        is_spliced=[False, False, False],
    )

    prior_table = assemble_priors(
        multi_loci=multi_loci,
        em_data=em_data,
        index=index,
        calibration=calibration,
        em_config=EMConfig(rna_call_bias=0.5),
    )

    np.testing.assert_allclose(prior_table.prior_n_local_gdna, [10.0, 0.0, 0.0])
    np.testing.assert_allclose(prior_table.prior_n_local_rna, [0.0, 20.0, 0.0])
    np.testing.assert_allclose(prior_table.prior_unspliced_total, [10.0, 20.0, 30.0])
    assert prior_table.alpha_gdna_add[0] == pytest.approx(10.0)
    assert prior_table.alpha_rna_add[1] == pytest.approx(20.0)
    assert prior_table.enable_gdna.tolist() == [1, 1, 1]
    assert prior_table.n_units_used_for_diagnostics.tolist() == [1, 1, 1]


def test_assemble_priors_structural_gate_zeroes_locus_without_gdna_candidate() -> None:
    region_calibration = _region_calibration(
        [GDNA_STATE, RNA_STATE, GDNA_STATE],
        [10.0, 20.0, 30.0],
    )
    index, calibration = _index_and_calibration(region_calibration)
    multi_loci = [_multi_locus(0, 0, 100), _multi_locus(1, 100, 200), _multi_locus(2, 200, 300)]
    em_data = _scored_fragments(
        gdna_log_liks=[-1.0, -1.0, np.nan],
        is_spliced=[False, False, False],
    )

    prior_table = assemble_priors(
        multi_loci=multi_loci,
        em_data=em_data,
        index=index,
        calibration=calibration,
    )

    assert prior_table.enable_gdna.tolist() == [1, 1, 0]
    assert prior_table.alpha_gdna_add[2] == pytest.approx(0.0)
    assert prior_table.alpha_rna_add[2] == pytest.approx(0.0)


def test_assemble_priors_applies_region_exposure_to_gdna_eff_len() -> None:
    region_calibration = _region_calibration(
        [GDNA_STATE, RNA_STATE, UNIFORM_STATE],
        [10.0, 20.0, 30.0],
        exposure=[0.25, 0.5, 2.0],
    )
    index, calibration = _index_and_calibration(region_calibration)
    multi_loci = [_multi_locus(0, 0, 100), _multi_locus(1, 100, 200), _multi_locus(2, 200, 300)]
    em_data = _scored_fragments(
        gdna_log_liks=[-1.0, -1.0, -1.0],
        is_spliced=[False, False, False],
    )

    prior_table = assemble_priors(
        multi_loci=multi_loci,
        em_data=em_data,
        index=index,
        calibration=calibration,
    )

    np.testing.assert_allclose(prior_table.gdna_em_exposure_weight, [0.25, 0.5, 2.0])
    np.testing.assert_allclose(
        prior_table.gdna_eff_len,
        np.maximum(prior_table.gdna_eff_len_unweighted * prior_table.gdna_em_exposure_weight, 1.0),
    )


def test_region_calibration_rejects_prior_mass_that_does_not_conserve_unspliced() -> None:
    with pytest.raises(ValueError, match="conserve"):
        PriorMassDeconvolution(
            unspliced_total=np.array([10.0], dtype=np.float32),
            gdna_unspliced_mean=np.array([4.0], dtype=np.float32),
            rna_unspliced_mean=np.array([7.0], dtype=np.float32),
            method=np.ones(1, dtype=np.uint8),
            precision=np.ones(1, dtype=np.float32),
            flags=np.zeros(1, dtype=np.uint16),
        )


def test_region_calibration_state_rows_must_match_latent_state_count() -> None:
    assert N_STATES == 4
    with pytest.raises(ValueError, match="p_states"):
        _region_calibration([[1.0, 0.0, 0.0]], [1.0])