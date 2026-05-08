"""Phase 0 (Bayesian-prior redesign): gDNA component eligibility is now an
explicit per-locus boolean (``locus_enable_gdna``), decoupled from
``alpha_gdna > 0``.

These tests pin the new native semantics:

* ``alpha_gdna == 0`` with ``enable_gdna == True`` keeps the gDNA component
  active (it can absorb posterior mass under the EM likelihood).
* ``enable_gdna == False`` disables the gDNA component regardless of the
  ``alpha_gdna`` value.
* When ``alpha_rna == 0`` (the asymmetric-prior default), the warm-start
  ratio override does not fire and the coverage-derived warm start is
  preserved.

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
    """``alpha_gdna`` no longer disables the gDNA component."""

    def test_zero_alpha_gdna_enabled_keeps_gdna_active(self):
        """``alpha_gdna == 0`` and ``enable_gdna == True`` ⇒ gDNA absorbs mass.

        Pre-Phase-0 behavior: ``alpha_gdna == 0`` zeroed the gDNA prior in
        the extractor, leaving the gDNA component permanently off and forcing
        all posterior mass onto transcripts.
        Post-Phase-0 behavior: with finite gDNA likelihoods, the gDNA
        component is eligible and absorbs a non-trivial share of mass.
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
            alpha_gdna=np.array([0.0], dtype=np.float64),
            alpha_rna=np.array([0.0], dtype=np.float64),
            index=None,
            enable_gdna=np.array([1], dtype=np.uint8),
        )

        # gDNA likelihood dominates ⇒ most mass should land on gDNA.
        assert total_gdna > 0.5 * n_units, (
            f"expected gDNA to absorb >50% with finite likelihoods and "
            f"enable_gdna=True; got total_gdna={total_gdna} of {n_units} units"
        )

    def test_disabled_gdna_yields_zero_gdna_mass(self):
        """``enable_gdna == False`` zeros gDNA mass regardless of alpha_gdna."""
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
            alpha_gdna=np.array([100.0], dtype=np.float64),  # large prior!
            alpha_rna=np.array([1.0], dtype=np.float64),
            index=None,
            enable_gdna=np.array([0], dtype=np.uint8),
        )

        assert total_gdna == 0.0, (
            f"expected zero gDNA mass with enable_gdna=False; got {total_gdna}"
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

        # All three with alpha_gdna=0 — pre-Phase-0, all three would have
        # gDNA disabled. Post-Phase-0, only the unspliced+finite case
        # produces gDNA assignments.
        g_spl, _, _ = rc.run_batch_locus_em_partitioned(
            [spliced_part], locus_t_lists,
            np.zeros(1), np.zeros(1), index=None,
        )
        g_uns, _, _ = rc2.run_batch_locus_em_partitioned(
            [unspliced_part], locus_t_lists,
            np.zeros(1), np.zeros(1), index=None,
        )
        g_no, _, _ = rc3.run_batch_locus_em_partitioned(
            [nogdna_part], locus_t_lists,
            np.zeros(1), np.zeros(1), index=None,
        )

        assert g_spl == 0.0, "spliced partition has no gDNA candidates"
        assert g_no == 0.0, "non-finite gDNA log-liks ⇒ no gDNA candidates"
        assert g_uns > 0.0, (
            "unspliced+finite gDNA log-lik must enable component "
            "even when alpha_gdna == 0 (Phase 0 semantics)"
        )


class TestWarmStartWithZeroAlphaRna:
    """``alpha_rna == 0`` no longer triggers the gDNA warm-start override."""

    def test_zero_alpha_rna_keeps_coverage_warm_start(self):
        """With ``alpha_rna == 0``, the gDNA warm-start ratio override is
        skipped. The result is sensitive to the starting point only via
        EM convergence, but at minimum the call must succeed and produce
        finite outputs (no NaN from a divide-by-zero ratio).
        """
        n_t = 2
        rc = _make_estimator(n_t)
        partition = _make_one_locus_partition(
            n_units=30, n_transcripts=n_t, gdna_log_lik=-1.0,
        )
        locus_t_lists = [np.arange(n_t, dtype=np.int32)]

        total_gdna, locus_rna, locus_gdna = rc.run_batch_locus_em_partitioned(
            partition_tuples=[partition],
            locus_transcript_indices=locus_t_lists,
            alpha_gdna=np.array([5.0], dtype=np.float64),  # nonzero prior
            alpha_rna=np.array([0.0], dtype=np.float64),    # asymmetric default
            index=None,
        )

        # Outputs must be finite and sum to total fragments.
        assert np.isfinite(total_gdna)
        assert np.isfinite(locus_rna).all() and np.isfinite(locus_gdna).all()
        assert total_gdna >= 0.0
        # Conservation up to assignment fractional remainder.
        assigned = float(locus_rna[0]) + float(locus_gdna[0])
        assert assigned == pytest.approx(30.0, rel=1e-9)
