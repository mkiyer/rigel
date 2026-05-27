"""Production-path tests for the native batch locus EM wrapper."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator


def _partition(
    *,
    n_units: int,
    log_liks: tuple[float, ...],
    gdna_log_lik: float = -np.inf,
    is_spliced: bool = False,
):
    n_t = len(log_liks)
    offsets = np.arange(n_units + 1, dtype=np.int64) * n_t
    return (
        offsets,
        np.tile(np.arange(n_t, dtype=np.int32), n_units),
        np.tile(np.asarray(log_liks, dtype=np.float64), n_units),
        np.ones(n_units * n_t, dtype=np.float64),
        np.zeros(n_units * n_t, dtype=np.uint8),
        np.full(n_units, 1 if is_spliced else 0, dtype=np.uint8),
        np.full(n_units, gdna_log_lik, dtype=np.float64),
        np.zeros(n_units, dtype=np.int32),
        np.zeros(n_units, dtype=np.uint8),
    )


def _estimator(n_t: int, *, mode: str = "vbem") -> AbundanceEstimator:
    est = AbundanceEstimator(
        num_transcripts=n_t,
        em_config=EMConfig(
            mode=mode,
            iterations=200,
            convergence_delta=1e-9,
            assignment_mode="fractional",
            seed=0,
        ),
    )
    est._t_eff_len = np.ones(n_t, dtype=np.float64)
    return est


def test_equal_rna_likelihoods_split_evenly():
    n_units = 100
    est = _estimator(2, mode="map")
    total_gdna, locus_rna, locus_gdna = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(0.0, 0.0))],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        enable_gdna=np.array([0], dtype=np.uint8),
    )

    assert total_gdna == pytest.approx(0.0)
    assert locus_gdna[0] == pytest.approx(0.0)
    assert locus_rna[0] == pytest.approx(n_units)
    np.testing.assert_allclose(est.em_counts.sum(axis=1), [50.0, 50.0], atol=1e-6)


def test_enabled_gdna_component_absorbs_likelihood_mass_without_prior_count():
    n_units = 80
    est = _estimator(2, mode="map")
    total_gdna, locus_rna, locus_gdna = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(-1.0, -1.0), gdna_log_lik=0.0)],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    assert total_gdna > 0.8 * n_units
    assert locus_gdna[0] == pytest.approx(total_gdna)
    assert locus_rna[0] + locus_gdna[0] == pytest.approx(n_units)


def test_gdna_effective_length_downweights_gdna_component():
    n_units = 100
    partition = _partition(n_units=n_units, log_liks=(0.0,), gdna_log_lik=0.0)

    est_short = _estimator(1, mode="map")
    _total_short, _rna_short, gdna_short = est_short.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        gdna_eff_len=np.array([1.0], dtype=np.float64),
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    est_long = _estimator(1, mode="map")
    _total_long, _rna_long, gdna_long = est_long.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        gdna_eff_len=np.array([100.0], dtype=np.float64),
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    assert gdna_long[0] < gdna_short[0]


def test_aggregate_rna_prior_reduces_gdna_share_without_isoform_floor():
    n_units = 100
    partition = _partition(n_units=n_units, log_liks=(0.0,), gdna_log_lik=0.0)

    est_unprior = _estimator(1, mode="map")
    _total_unprior, _rna_unprior, gdna_unprior = est_unprior.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        rna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        gdna_eff_len=np.array([1.0], dtype=np.float64),
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    est_rna = _estimator(1, mode="map")
    _total_rna, _rna_rna, gdna_rna = est_rna.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        rna_prior_count=np.array([10.0], dtype=np.float64),
        index=None,
        gdna_eff_len=np.array([1.0], dtype=np.float64),
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    assert gdna_rna[0] < gdna_unprior[0]


def test_grouped_priors_inactive_without_structural_gdna_candidate():
    n_units = 40
    est = _estimator(1, mode="map")
    total_gdna, locus_rna, locus_gdna = est.run_batch_locus_em_partitioned(
        partition_tuples=[
            _partition(n_units=n_units, log_liks=(0.0,), gdna_log_lik=0.0, is_spliced=True)
        ],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([100.0], dtype=np.float64),
        rna_prior_count=np.array([100.0], dtype=np.float64),
        index=None,
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    assert total_gdna == pytest.approx(0.0)
    assert locus_gdna[0] == pytest.approx(0.0)
    assert locus_rna[0] == pytest.approx(n_units)


def test_assignment_outputs_follow_partition_units():
    n_units = 5
    est = _estimator(2, mode="map")
    result = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(0.0, -4.0))],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        enable_gdna=np.array([0], dtype=np.uint8),
        emit_assignments=True,
    )

    total_gdna, _locus_rna, _locus_gdna, winner_tid, winner_post, n_candidates = result
    assert total_gdna == pytest.approx(0.0)
    np.testing.assert_array_equal(winner_tid, np.zeros(n_units, dtype=np.int32))
    assert np.all(winner_post > 0.98)
    np.testing.assert_array_equal(n_candidates, np.full(n_units, 2, dtype=np.int16))
