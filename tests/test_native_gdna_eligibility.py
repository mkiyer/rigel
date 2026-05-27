"""Grouped prior redesign: native gDNA availability is structural.

These tests pin the new native semantics:

* ``gdna_prior_count == 0`` with finite unspliced gDNA candidates keeps the
    gDNA component active.
* A compatibility ``enable_gdna`` array no longer acts as a modeling gate;
    native derives candidate availability from the partition itself.

The tests construct a minimal one-locus partition by hand to avoid the
calibration / scoring machinery; they exercise the public Python wrapper
``AbundanceEstimator.run_batch_locus_em_partitioned``.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator


# ---------------------------------------------------------------------------
# Minimal one-locus partition factory
# ---------------------------------------------------------------------------

def _make_one_locus_partition(
    *,
    n_units: int,
    n_transcripts: int,
    log_lik_t0: float = -1.0,
    log_lik_t1: float = -2.0,
    gdna_log_lik: float = -1.5,
    is_spliced: bool = False,
):
    """Build a one-locus partition tuple with all units identical.

    Each unit carries a candidate for transcript 0, transcript 1, and
    (when ``is_spliced`` is False and ``gdna_log_lik`` is finite) gDNA.
    """
    assert n_transcripts >= 2, "need at least two transcripts to exercise EM"

    n_cands_per_unit = n_transcripts
    total_cands = n_units * n_cands_per_unit

    offsets = np.arange(n_units + 1, dtype=np.int64) * n_cands_per_unit
    t_indices = np.tile(
        np.arange(n_transcripts, dtype=np.int32), n_units
    )
    log_liks = np.tile(
        np.array(
            [log_lik_t0] + [log_lik_t1] * (n_transcripts - 1),
            dtype=np.float64,
        ),
        n_units,
    )
    coverage_weights = np.ones(total_cands, dtype=np.float64)
    count_cols = np.zeros(total_cands, dtype=np.uint8)

    is_spliced_arr = np.full(
        n_units, 1 if is_spliced else 0, dtype=np.uint8
    )
    gdna_log_liks_arr = np.full(n_units, gdna_log_lik, dtype=np.float64)
    locus_t_indices = np.zeros(n_units, dtype=np.int32)  # all anchored on t0
    locus_count_cols = np.zeros(n_units, dtype=np.uint8)

    return (
        offsets,
        t_indices,
        log_liks,
        coverage_weights,
        count_cols,
        is_spliced_arr,
        gdna_log_liks_arr,
        locus_t_indices,
        locus_count_cols,
    )


def _make_estimator(n_transcripts: int) -> AbundanceEstimator:
    rc = AbundanceEstimator(
        num_transcripts=n_transcripts,
        em_config=EMConfig(iterations=200, convergence_delta=1e-9, seed=0),
    )
    rc._t_eff_len = np.full(n_transcripts, 1000.0, dtype=np.float64)
    return rc


# ---------------------------------------------------------------------------
# Phase 0 invariants
# ---------------------------------------------------------------------------

class TestEligibilityDecoupling:
    """``gdna_prior_count`` no longer disables the gDNA component."""

    def test_zero_gdna_prior_enabled_keeps_gdna_active(self):
        """``gdna_prior_count == 0`` and ``enable_gdna == True`` ⇒ gDNA absorbs mass.

        With finite gDNA likelihoods, the gDNA component is eligible and
        absorbs a non-trivial share of mass even when the calibration prior
        count is zero.
        """
        n_units = 100
        n_t = 2
        rc = _make_estimator(n_t)

        partition = _make_one_locus_partition(
            n_units=n_units, n_transcripts=n_t,
            log_lik_t0=-1.0, log_lik_t1=-1.0,   # identical RNA likelihoods
            gdna_log_lik=-0.5,                   # gDNA wins the likelihood
        )
        locus_t_lists = [np.arange(n_t, dtype=np.int32)]

        total_gdna, _rna, _g = rc.run_batch_locus_em_partitioned(
            partition_tuples=[partition],
            locus_transcript_indices=locus_t_lists,
            gdna_prior_count=np.array([0.0], dtype=np.float64),
            index=None,
            enable_gdna=np.array([1], dtype=np.uint8),
        )

        # gDNA likelihood dominates ⇒ most mass should land on gDNA.
        assert total_gdna > 0.5 * n_units, (
            f"expected gDNA to absorb >50% with finite likelihoods and "
            f"enable_gdna=True; got total_gdna={total_gdna} of {n_units} units"
        )

    def test_compat_enable_false_does_not_disable_structural_candidate(self):
        """A compatibility ``enable_gdna=False`` input is ignored by native v3."""
        n_units = 50
        n_t = 2
        rc = _make_estimator(n_t)

        partition = _make_one_locus_partition(
            n_units=n_units, n_transcripts=n_t,
            log_lik_t0=-1.0, log_lik_t1=-1.0,
            gdna_log_lik=-0.5,
        )
        locus_t_lists = [np.arange(n_t, dtype=np.int32)]

        total_gdna, _rna, _g = rc.run_batch_locus_em_partitioned(
            partition_tuples=[partition],
            locus_transcript_indices=locus_t_lists,
            gdna_prior_count=np.array([100.0], dtype=np.float64),
            index=None,
            enable_gdna=np.array([0], dtype=np.uint8),
        )

        assert total_gdna > 0.5 * n_units, (
            "native should derive gDNA availability from finite unspliced candidates, "
            f"not the compatibility enable_gdna array; got total_gdna={total_gdna}"
        )

    def test_default_enable_gdna_inferred_from_partition(self):
        """When ``enable_gdna`` is None, the wrapper computes it from the
        partition (any unspliced unit with a finite gDNA log-lik ⇒ enabled).
        """
        n_t = 2
        rc = _make_estimator(n_t)

        # All-spliced partition ⇒ no per-unit gDNA candidate ⇒ enable=0.
        spliced_part = _make_one_locus_partition(
            n_units=20, n_transcripts=n_t,
            is_spliced=True,
        )
        # Unspliced partition with finite gDNA log-liks ⇒ enable=1.
        unspliced_part = _make_one_locus_partition(
            n_units=20, n_transcripts=n_t,
            is_spliced=False, gdna_log_lik=-0.5,
        )
        # Unspliced but non-finite gDNA log-lik ⇒ enable=0.
        nogdna_part = _make_one_locus_partition(
            n_units=20, n_transcripts=n_t,
            is_spliced=False, gdna_log_lik=-np.inf,
        )

        rc2 = _make_estimator(n_t)
        rc3 = _make_estimator(n_t)
        locus_t_lists = [np.arange(n_t, dtype=np.int32)]

        # All three with gdna_prior_count=0. Only the unspliced+finite
        # case produces gDNA assignments.
        g_spl, _, _ = rc.run_batch_locus_em_partitioned(
            [spliced_part], locus_t_lists,
            np.zeros(1), index=None,
        )
        g_uns, _, _ = rc2.run_batch_locus_em_partitioned(
            [unspliced_part], locus_t_lists,
            np.zeros(1), index=None,
        )
        g_no, _, _ = rc3.run_batch_locus_em_partitioned(
            [nogdna_part], locus_t_lists,
            np.zeros(1), index=None,
        )

        assert g_spl == 0.0, "spliced partition has no gDNA candidates"
        assert g_no == 0.0, "non-finite gDNA log-liks ⇒ no gDNA candidates"
        assert g_uns > 0.0, (
            "unspliced+finite gDNA log-lik must enable component "
            "even when gdna_prior_count == 0"
        )


class TestGdnaPriorCount:
    """gDNA prior counts are accepted without any RNA-prior companion."""

    def test_positive_gdna_prior_produces_finite_outputs(self):
        """A positive gDNA prior count produces finite, conserved outputs."""
        n_t = 2
        rc = _make_estimator(n_t)
        partition = _make_one_locus_partition(
            n_units=30, n_transcripts=n_t, gdna_log_lik=-1.0,
        )
        locus_t_lists = [np.arange(n_t, dtype=np.int32)]

        total_gdna, locus_rna, locus_gdna = rc.run_batch_locus_em_partitioned(
            partition_tuples=[partition],
            locus_transcript_indices=locus_t_lists,
            gdna_prior_count=np.array([5.0], dtype=np.float64),
            index=None,
        )

        # Outputs must be finite and sum to total fragments.
        assert np.isfinite(total_gdna)
        assert np.isfinite(locus_rna).all() and np.isfinite(locus_gdna).all()
        assert total_gdna >= 0.0
        # Conservation up to assignment fractional remainder.
        assigned = float(locus_rna[0]) + float(locus_gdna[0])
        assert assigned == pytest.approx(30.0, rel=1e-9)
