"""Phase 0.5 (Bayesian-prior redesign): ``prior_weight_rna`` policy lock.

``prior_weight_rna`` is an **RNA-prior-allocation weight**. It distributes
a positive ``alpha_rna`` budget across RNA components proportional to
coverage × weight. With ``alpha_rna == 0`` (the asymmetric-prior default),
this allocation path collapses and ``prior_weight_rna`` must be
objective-neutral.

This test pins that contract:

    Varying ``prior_weight_rna`` while ``alpha_rna == 0`` produces
    bit-identical EM output.

If a future change accidentally couples ``prior_weight_rna`` to the
likelihood (a tempting but theoretically wrong way to keep nRNA
suppression alive under the asymmetric prior), this test will fail.

See ``docs/bayesian_prior/bayesian_prior_plan_v3.md`` §5 Phase 0.5.
"""

from __future__ import annotations

import numpy as np

from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator


def _make_one_locus_partition(n_units: int, n_transcripts: int):
    n_cands_per_unit = n_transcripts
    total = n_units * n_cands_per_unit

    offsets = np.arange(n_units + 1, dtype=np.int64) * n_cands_per_unit
    t_indices = np.tile(np.arange(n_transcripts, dtype=np.int32), n_units)
    # Asymmetric likelihoods so the EM has work to do.
    log_liks = np.tile(
        np.array([-1.0, -2.0, -1.5][:n_transcripts], dtype=np.float64),
        n_units,
    )
    coverage_weights = np.ones(total, dtype=np.float64)
    count_cols = np.zeros(total, dtype=np.uint8)
    is_spliced = np.zeros(n_units, dtype=np.uint8)
    gdna_log_liks = np.full(n_units, -1.0, dtype=np.float64)
    locus_t_indices = np.zeros(n_units, dtype=np.int32)
    locus_count_cols = np.zeros(n_units, dtype=np.uint8)

    return (
        offsets, t_indices, log_liks, coverage_weights, count_cols,
        is_spliced, gdna_log_liks, locus_t_indices, locus_count_cols,
    )


def _make_estimator(n_transcripts: int) -> AbundanceEstimator:
    rc = AbundanceEstimator(
        num_transcripts=n_transcripts,
        em_config=EMConfig(iterations=200, convergence_delta=1e-12, seed=0),
    )
    rc._t_eff_len = np.full(n_transcripts, 1000.0, dtype=np.float64)
    return rc


def _run(rc: AbundanceEstimator, partition, weights: np.ndarray | None,
         alpha_rna: float):
    locus_t_lists = [np.arange(rc.num_transcripts, dtype=np.int32)]
    pwr = [weights] if weights is not None else None
    total_gdna, locus_rna, locus_gdna = rc.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=locus_t_lists,
        alpha_gdna=np.array([0.0], dtype=np.float64),
        alpha_rna=np.array([alpha_rna], dtype=np.float64),
        index=None,
        locus_prior_weight_rna=pwr,
    )
    return total_gdna, np.asarray(locus_rna), np.asarray(locus_gdna), \
        rc.em_counts.copy()


class TestPriorWeightRnaPolicyLock:
    """``prior_weight_rna`` is objective-neutral when ``alpha_rna == 0``."""

    def test_varying_weights_bitidentical_under_zero_alpha_rna(self):
        n_t = 3
        partition = _make_one_locus_partition(n_units=80, n_transcripts=n_t)

        # Three different weight vectors: uniform, all-zeros, skewed.
        # Under ``alpha_rna == 0`` they must all produce the *same* EM
        # output. ``len = n_t + 1`` (RNA components + gDNA slot).
        weights_a = np.ones(n_t + 1, dtype=np.float32)
        weights_b = np.zeros(n_t + 1, dtype=np.float32)
        weights_c = np.array([0.1, 5.0, 1.0, 0.0], dtype=np.float32)
        weights_none = None  # nullptr path → all-ones inside C++

        results = []
        for w in (weights_a, weights_b, weights_c, weights_none):
            rc = _make_estimator(n_t)
            results.append(_run(rc, partition, w, alpha_rna=0.0))

        gdna_ref, rna_ref, gdna_arr_ref, em_ref = results[0]
        for gdna, rna, gdna_arr, em in results[1:]:
            assert gdna == gdna_ref, (
                "prior_weight_rna must not affect total gDNA when alpha_rna=0"
            )
            np.testing.assert_array_equal(rna, rna_ref)
            np.testing.assert_array_equal(gdna_arr, gdna_arr_ref)
            np.testing.assert_array_equal(em, em_ref)
