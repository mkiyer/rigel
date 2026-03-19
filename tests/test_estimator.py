"""Tests for rigel.estimator — AbundanceEstimator with locus-level EM assignment."""

import numpy as np
import pandas as pd
import pytest

from rigel.types import Strand
from rigel.splice import (
    SpliceType,
    SpliceStrandCol,
)
from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator
from rigel.scored_fragments import ScoredFragments

from conftest import _UNSPLICED_SENSE, _make_locus_em_data, _run_and_assign


# =====================================================================
# Helpers — lightweight mock index + fixtures
# =====================================================================


class MockIndex:
    """Minimal mock of TranscriptIndex with the arrays AbundanceEstimator needs."""

    def __init__(
        self,
        num_transcripts,
        num_genes,
        t_to_g,
        t_to_strand,
        g_to_strand,
        t_ids=None,
        g_ids=None,
        g_names=None,
        t_gnames=None,
    ):
        self.num_transcripts = num_transcripts
        self.num_genes = num_genes
        self.t_to_g_arr = np.array(t_to_g, dtype=np.int64)
        self.t_to_strand_arr = np.array(t_to_strand, dtype=np.int64)
        self.g_to_strand_arr = np.array(g_to_strand, dtype=np.int64)

        # DataFrames for output methods
        if t_ids is None:
            t_ids = [f"t{i}" for i in range(num_transcripts)]
        if g_ids is None:
            g_ids = [f"g{i}" for i in range(num_genes)]
        if g_names is None:
            g_names = [f"Gene{i}" for i in range(num_genes)]
        if t_gnames is None:
            t_gnames = [g_names[g] for g in t_to_g]

        self.t_df = pd.DataFrame(
            {
                "t_id": t_ids,
                "g_id": [g_ids[g] for g in t_to_g],
                "g_name": t_gnames,
            }
        )
        self.g_df = pd.DataFrame(
            {
                "g_id": g_ids,
                "g_name": g_names,
            }
        )


def _make_index():
    """3 transcripts, 2 genes.

    t0, t1 → g0 (+strand)
    t2     → g1 (-strand)
    """
    return MockIndex(
        num_transcripts=3,
        num_genes=2,
        t_to_g=[0, 0, 1],
        t_to_strand=[int(Strand.POS), int(Strand.POS), int(Strand.NEG)],
        g_to_strand=[int(Strand.POS), int(Strand.NEG)],
    )


def _make_frag_length_models():
    from rigel.frag_length_model import FragmentLengthModels

    im = FragmentLengthModels()
    im.observe(250, SpliceType.UNSPLICED)
    return im


def _make_em_data(
    t_indices_per_unit,
    log_liks_per_unit=None,
    count_cols_per_unit=None,
    num_transcripts=None,
):
    """Build ScoredFragments from a list of per-unit candidate lists.

    The global ScoredFragments contains mRNA + nRNA candidates only (no gDNA).
    """
    offsets = [0]
    flat_t = []
    flat_lk = []
    flat_cc = []

    for u, t_list in enumerate(t_indices_per_unit):
        for j, t_idx in enumerate(t_list):
            flat_t.append(t_idx)
            if log_liks_per_unit is not None:
                flat_lk.append(log_liks_per_unit[u][j])
            else:
                flat_lk.append(0.0)
            if count_cols_per_unit is not None:
                flat_cc.append(count_cols_per_unit[u][j])
            else:
                flat_cc.append(_UNSPLICED_SENSE)
        offsets.append(len(flat_t))

    n_units = len(t_indices_per_unit)
    n_candidates = len(flat_t)

    if num_transcripts is None:
        num_transcripts = (max(flat_t) + 1) if flat_t else 0
    nrna_base = num_transcripts

    # Build locus tracking arrays
    locus_t = np.full(n_units, -1, dtype=np.int32)
    locus_cc = np.zeros(n_units, dtype=np.uint8)
    for u, t_list in enumerate(t_indices_per_unit):
        cc_list = count_cols_per_unit[u] if count_cols_per_unit else None
        for j, t_idx in enumerate(t_list):
            if t_idx < nrna_base:
                locus_t[u] = t_idx
                locus_cc[u] = cc_list[j] if cc_list else _UNSPLICED_SENSE
                break

    return ScoredFragments(
        offsets=np.array(offsets, dtype=np.int64),
        t_indices=np.array(flat_t, dtype=np.int32),
        log_liks=np.array(flat_lk, dtype=np.float64),
        count_cols=np.array(flat_cc, dtype=np.uint8),
        coverage_weights=np.ones(n_candidates, dtype=np.float64),
        tx_starts=np.zeros(n_candidates, dtype=np.int32),
        tx_ends=np.ones(n_candidates, dtype=np.int32),
        locus_t_indices=locus_t,
        locus_count_cols=locus_cc,
        is_spliced=np.zeros(n_units, dtype=bool),
        gdna_log_liks=np.full(n_units, -np.inf, dtype=np.float64),
        frag_ids=np.arange(n_units, dtype=np.int64),
        frag_class=np.zeros(n_units, dtype=np.int8),
        splice_type=np.zeros(n_units, dtype=np.uint8),
        n_units=n_units,
        n_candidates=n_candidates,
        nrna_base_index=nrna_base,
        genomic_footprints=np.full(n_units, 200, dtype=np.int32),
    )


# =====================================================================
# ScoredFragments construction
# =====================================================================


class TestScoredFragments:
    def test_empty_em_data(self):
        em = _make_em_data([])
        assert em.n_units == 0
        assert em.n_candidates == 0
        assert len(em.offsets) == 1
        assert em.offsets[0] == 0

    def test_single_unit_two_candidates(self):
        em = _make_em_data([[0, 1]])
        assert em.n_units == 1
        assert em.n_candidates == 2
        assert list(em.offsets) == [0, 2]
        assert list(em.t_indices) == [0, 1]

    def test_multiple_units(self):
        em = _make_em_data([[0, 1], [2], [0, 1, 2]])
        assert em.n_units == 3
        assert em.n_candidates == 6
        assert list(em.offsets) == [0, 2, 3, 6]

    def test_new_fields_present(self):
        """New locus-EM fields are populated."""
        em = _make_em_data([[0, 1]])
        assert em.is_spliced.shape == (1,)
        assert em.gdna_log_liks.shape == (1,)
        assert em.gdna_log_liks[0] == -np.inf  # default unspliced


# =====================================================================
# Locus EM — convergence tests
# =====================================================================


class TestLocusEM:
    def test_empty_locus_em(self):
        """Empty ScoredFragments produces zero em_counts."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data([], num_transcripts=3)
        _run_and_assign(rc, bundle, em_iterations=5)
        assert rc.em_counts.sum() == 0.0

    def test_uniform_priors_stay_uniform(self):
        """With no unambig counts and equal likelihoods, mRNA em_counts stay uniform."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data([[0, 1]] * 1000, num_transcripts=3)
        _run_and_assign(rc, bundle, em_iterations=10)

        # mRNA transcripts 0 and 1 should receive approximately equal counts
        t0_em = rc.em_counts[0].sum()
        t1_em = rc.em_counts[1].sum()
        assert t0_em == pytest.approx(t1_em, rel=0.05)
        # t2 not in any unit → no counts
        assert rc.em_counts[2].sum() < t0_em * 0.01

    def test_unambig_counts_bias_em(self):
        """Unique counts on t0 should bias EM toward t0."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 100.0

        bundle = _make_locus_em_data(
            [[0, 1]] * 500,
            num_transcripts=3,
            rc=rc,
        )
        _run_and_assign(rc, bundle, em_iterations=10)

        assert rc.em_counts[0].sum() > rc.em_counts[1].sum()

    def test_zero_iterations_uses_unambig_priors(self):
        """With em_iterations=0, assignment comes from unambig counts + prior only."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 10.0
        rc.unambig_counts[1, 0] = 5.0

        bundle = _make_locus_em_data(
            [[0, 1]] * 100,
            num_transcripts=3,
            rc=rc,
        )
        _run_and_assign(rc, bundle, em_iterations=0)

        # t0 gets more prior weight than t1
        assert rc.em_counts[0].sum() > rc.em_counts[1].sum()
        assert rc.em_counts.sum() == pytest.approx(100.0, abs=1.0)

    def test_likelihood_influences_em(self):
        """Candidates with higher likelihoods should attract more mass."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data(
            [[0, 1]] * 1000,
            log_liks_per_unit=[[0.0, -10.0]] * 1000,
            num_transcripts=3,
        )
        _run_and_assign(rc, bundle, em_iterations=10)

        assert rc.em_counts[0].sum() > rc.em_counts[1].sum()

    def test_convergence_detected(self):
        """EM should converge early for simple problems."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 100.0
        rc.unambig_counts[1, 0] = 100.0

        bundle = _make_locus_em_data(
            [[0, 1]] * 100,
            num_transcripts=2,
            rc=rc,
        )
        _run_and_assign(rc, bundle, em_iterations=100)
        # Just verify it runs without error and produces sensible output
        assert rc.em_counts.sum() == pytest.approx(100.0, abs=1.0)


# =====================================================================
# Locus assignment — posterior expected-count assignment
# =====================================================================


class TestLocusAssignment:
    def test_single_candidate_gets_full_count(self):
        """Unit with one candidate → count goes entirely to that transcript."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data([[0]], num_transcripts=3)
        _run_and_assign(rc, bundle, em_iterations=1)

        assert rc.em_counts[0].sum() == pytest.approx(1.0)
        assert rc.em_counts.sum() == pytest.approx(1.0)

    def test_two_candidates_assigns_one(self):
        """Each unit contributes exactly 1.0 total count."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data([[0, 1]], num_transcripts=3)
        _run_and_assign(rc, bundle, em_iterations=1)

        assert rc.em_counts.sum() == pytest.approx(1.0)

    def test_n_units_equals_n_counts(self):
        """N ambiguous units → N total em_counts."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data([[0, 1]] * 100, num_transcripts=3)
        _run_and_assign(rc, bundle, em_iterations=5)

        assert rc.em_counts.sum() == pytest.approx(100.0, abs=1.0)

    def test_total_counts_equals_unambig_plus_em(self):
        """t_counts = unambig_counts + em_counts."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))

        # 10 unambig (set directly — deterministic assignment now handled in C++)
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 10.0

        # 5 ambiguous
        bundle = _make_locus_em_data([[0, 1]] * 5, num_transcripts=3, rc=rc)
        _run_and_assign(rc, bundle, em_iterations=5)

        assert rc.unambig_counts.sum() == 10.0
        assert rc.em_counts.sum() == pytest.approx(5.0, abs=0.5)
        np.testing.assert_array_almost_equal(
            rc.t_counts, rc.unambig_counts + rc.em_counts, decimal=1
        )

    def test_distribution_follows_priors(self):
        """Over many fragments, expected-count assignments follow unambig priors."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 90.0
        rc.unambig_counts[1, _UNSPLICED_SENSE] = 10.0

        bundle = _make_locus_em_data(
            [[0, 1]] * 10000,
            num_transcripts=3,
            rc=rc,
        )
        _run_and_assign(rc, bundle, em_iterations=10)

        t0_em = rc.em_counts[0].sum()
        t1_em = rc.em_counts[1].sum()
        assert t0_em > t1_em
        # Both should receive counts (not winner-take-all)
        assert t1_em > 500  # ~10% of 10000

    def test_writes_to_em_not_unambig(self):
        """Batch EM only writes to em_counts."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data([[0, 1]] * 50, num_transcripts=3)

        unique_before = rc.unambig_counts.copy()
        _run_and_assign(rc, bundle, em_iterations=5)

        np.testing.assert_array_equal(rc.unambig_counts, unique_before)
        assert rc.em_counts.sum() == pytest.approx(50.0, abs=1.0)

    def test_empty_em_data_does_nothing(self):
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data([], num_transcripts=3)
        _run_and_assign(rc, bundle, em_iterations=5)

        assert rc.em_counts.sum() == 0.0

    def test_splice_strand_col_respected(self):
        """Count goes to the correct column (category×strand)."""
        cc = int(SpliceStrandCol.SPLICED_ANNOT_ANTISENSE)
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data(
            [[0]],
            count_cols_per_unit=[[cc]],
            num_transcripts=3,
        )
        _run_and_assign(rc, bundle, em_iterations=1)

        assert rc.em_counts[0, cc] == pytest.approx(1.0, abs=0.1)
        assert rc.em_counts.sum() == pytest.approx(1.0, abs=0.1)

    def test_deterministic_with_same_seed(self):
        """Expected-count assignment is deterministic."""
        results = []
        for _ in range(3):
            rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
            rc.unambig_counts[0, _UNSPLICED_SENSE] = 50.0
            rc.unambig_counts[1, _UNSPLICED_SENSE] = 50.0
            bundle = _make_locus_em_data(
                [[0, 1]] * 100,
                num_transcripts=3,
                rc=rc,
            )
            _run_and_assign(rc, bundle, em_iterations=10)
            results.append(rc.em_counts.copy())

        np.testing.assert_array_equal(results[0], results[1])
        np.testing.assert_array_equal(results[1], results[2])

    def test_different_seeds_same_expected_counts(self):
        """Different seeds do not affect expected-count assignment."""
        counts = []
        for seed in [1, 2]:
            rc = AbundanceEstimator(3, em_config=EMConfig(seed=seed))
            rc.unambig_counts[0, _UNSPLICED_SENSE] = 50.0
            rc.unambig_counts[1, _UNSPLICED_SENSE] = 50.0
            bundle = _make_locus_em_data(
                [[0, 1]] * 100,
                num_transcripts=3,
                rc=rc,
            )
            _run_and_assign(rc, bundle, em_iterations=10)
            counts.append(rc.em_counts.copy())

        np.testing.assert_array_equal(counts[0], counts[1])

    def test_fractional_counts(self):
        """EM counts may be fractional under expected-count assignment."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 50.0
        rc.unambig_counts[1, _UNSPLICED_SENSE] = 30.0

        bundle = _make_locus_em_data(
            [[0, 1]] * 200,
            num_transcripts=3,
            rc=rc,
        )
        _run_and_assign(rc, bundle, em_iterations=10)

        frac = rc.em_counts - np.floor(rc.em_counts)
        assert np.any(frac > 0.0)

    def test_equal_candidates_share_counts(self):
        """Equal-probability candidates all receive counts.

        With 4 candidates of equal probability (~25%), all should
        receive a substantial share over many fragments.
        """
        rc = AbundanceEstimator(4, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data(
            [[0, 1, 2, 3]] * 10000,
            num_transcripts=4,
        )
        _run_and_assign(rc, bundle, em_iterations=10)

        per_t = rc.em_counts.sum(axis=1)
        # Each should get ~2500 counts; verify all get at least 2000
        for t in range(4):
            assert per_t[t] > 2000, f"t{t} got {per_t[t]} counts — expected ~2500"
        assert rc.em_counts.sum() == pytest.approx(10000.0, abs=10.0)


# =====================================================================
# Posterior mean
# =====================================================================


class TestPosteriorMean:
    def test_no_em_returns_nan(self):
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        pm = rc.posterior_mean()
        assert np.all(np.isnan(pm))

    def test_single_candidate_posterior_one(self):
        """One candidate per unit → posterior is 1.0, mean is 1.0."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data([[0]] * 10, num_transcripts=3)
        _run_and_assign(rc, bundle, em_iterations=1)

        pm = rc.posterior_mean()
        assert pm[0] == pytest.approx(1.0)
        assert np.isnan(pm[1])  # no assignments
        assert np.isnan(pm[2])

    def test_strong_prior_high_posterior(self):
        """Strong prior → assigned units have high mean posterior."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 100.0
        bundle = _make_locus_em_data(
            [[0, 1]] * 100,
            num_transcripts=2,
            rc=rc,
        )
        _run_and_assign(rc, bundle, em_iterations=10)

        pm = rc.posterior_mean()
        # t0 gets most assignments with high posterior
        assert pm[0] > 0.8


# =====================================================================
# Simultaneous resolution — no phase ordering bias
# =====================================================================


class TestSimultaneousResolution:
    """Verify that isoform-ambig and gene-ambig fragments
    are resolved simultaneously, not sequentially."""

    def test_same_counts_regardless_of_order(self):
        """Shuffling ambiguous units produces identical EM counts.

        EM convergence is order-independent.  Expected counts agree.
        """
        rc1 = AbundanceEstimator(4, em_config=EMConfig(seed=42))
        rc1.unambig_counts[0, 0] = 50.0
        rc1.unambig_counts[2, 0] = 30.0

        units = [[0, 1]] * 200 + [[0, 2]] * 200
        bundle1 = _make_locus_em_data(units, num_transcripts=4, rc=rc1)
        _run_and_assign(rc1, bundle1, em_iterations=10)

        rc2 = AbundanceEstimator(4, em_config=EMConfig(seed=42))
        rc2.unambig_counts[0, 0] = 50.0
        rc2.unambig_counts[2, 0] = 30.0

        units_rev = [[0, 2]] * 200 + [[0, 1]] * 200
        bundle2 = _make_locus_em_data(units_rev, num_transcripts=4, rc=rc2)
        _run_and_assign(rc2, bundle2, em_iterations=10)

        np.testing.assert_allclose(rc1.em_counts, rc2.em_counts, atol=1e-10)
        # Total counts should match (400 units each)
        assert rc1.em_counts.sum() == pytest.approx(400.0, abs=1.0)
        assert rc2.em_counts.sum() == pytest.approx(400.0, abs=1.0)


# =====================================================================
# Multimapper-like units (multi-alignment molecules)
# =====================================================================


class TestMultimapperEM:
    def test_multimapper_molecule_one_count(self):
        """Multimapper molecule (many candidates) → exactly 1.0 total count."""
        rc = AbundanceEstimator(4, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data([[0, 1, 2, 3]], num_transcripts=4)
        _run_and_assign(rc, bundle, em_iterations=5)

        assert rc.em_counts.sum() == pytest.approx(1.0, abs=0.1)

    def test_multimapper_and_ambig_together(self):
        """Mixed ambiguous + multimapper units in same EM."""
        rc = AbundanceEstimator(4, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 100.0

        units = [[0, 1]] * 50 + [[0, 2]] * 50 + [[0, 1, 2, 3]] * 25
        bundle = _make_locus_em_data(units, num_transcripts=4, rc=rc)
        _run_and_assign(rc, bundle, em_iterations=10)

        assert rc.em_counts.sum() == pytest.approx(125.0, abs=2.0)


# =====================================================================
# Primary counts DataFrame output
# =====================================================================


class TestCountsOutput:
    def test_get_counts_df_columns(self):
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        df = rc.get_counts_df(index)
        assert df.shape[0] == 3
        expected_cols = [
            "transcript_id",
            "gene_id",
            "gene_name",
            "locus_id",
            "effective_length",
            "mrna",
            "mrna_unambig",
            "mrna_em",
            "mrna_high_conf",
            "mrna_spliced",
            "nrna",
            "rna_total",
            "tpm",
            "posterior_mean",
        ]
        assert list(df.columns) == expected_cols

    def test_get_gene_counts_df_columns(self):
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        df = rc.get_gene_counts_df(index)
        assert df.shape[0] == 2
        expected_cols = [
            "gene_id",
            "gene_name",
            "locus_id",
            "effective_length",
            "mrna",
            "mrna_unambig",
            "mrna_em",
            "mrna_high_conf",
            "mrna_spliced",
            "nrna",
            "rna_total",
            "tpm",
        ]
        assert list(df.columns) == expected_cols

    def test_counts_include_both_sources(self):
        """Total count sums unambig + em counts."""
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 5.0
        rc.em_counts[0, _UNSPLICED_SENSE] = 3.0

        df = rc.get_counts_df(index)
        assert df.loc[0, "mrna"] == 8.0
        assert df.loc[0, "mrna_unambig"] == 5.0
        assert df.loc[0, "mrna_em"] == 3.0

    def test_gene_counts_aggregate(self):
        """Gene counts aggregate across transcripts."""
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 5.0
        rc.unambig_counts[1, _UNSPLICED_SENSE] = 3.0  # same gene
        rc.unambig_counts[2, _UNSPLICED_SENSE] = 7.0  # different gene

        df = rc.get_gene_counts_df(index)
        assert df.loc[0, "mrna"] == 8.0  # g0 = t0 + t1
        assert df.loc[1, "mrna"] == 7.0  # g1 = t2

    def test_gene_effective_length_abundance_weighted(self):
        """Gene effective length is abundance-weighted mean of transcript eff lens."""
        index = _make_index()  # t0,t1 → g0; t2 → g1
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        # Set different effective lengths per transcript
        rc._t_eff_len = np.array([500.0, 1000.0, 300.0])
        # t0 gets 90% of counts, t1 gets 10% → weighted toward t0
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 90.0
        rc.unambig_counts[1, _UNSPLICED_SENSE] = 10.0
        # t2 has counts too
        rc.unambig_counts[2, _UNSPLICED_SENSE] = 50.0

        df = rc.get_gene_counts_df(index)
        # g0: (90*500 + 10*1000) / (90+10) = 55000/100 = 550
        assert df.loc[0, "effective_length"] == pytest.approx(550.0)
        # g1: single transcript → 300
        assert df.loc[1, "effective_length"] == pytest.approx(300.0)

    def test_gene_effective_length_zero_counts_uses_mean(self):
        """Zero-count genes use unweighted mean of transcript effective lengths."""
        index = _make_index()  # t0,t1 → g0; t2 → g1
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc._t_eff_len = np.array([500.0, 1000.0, 300.0])
        # No counts at all → both genes are zero-count

        df = rc.get_gene_counts_df(index)
        # g0: mean(500, 1000) = 750
        assert df.loc[0, "effective_length"] == pytest.approx(750.0)
        # g1: mean(300) = 300
        assert df.loc[1, "effective_length"] == pytest.approx(300.0)

    def test_spliced_counts(self):
        """mrna_spliced captures both SPLICED_ANNOT and SPLICED_UNANNOT."""
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        # SPLICED_ANNOT sense
        rc.unambig_counts[0, SpliceStrandCol.SPLICED_ANNOT_SENSE] = 3.0
        # SPLICED_UNANNOT antisense
        rc.unambig_counts[0, SpliceStrandCol.SPLICED_UNANNOT_ANTISENSE] = 2.0
        # UNSPLICED (should not be in mrna_spliced)
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 10.0

        df = rc.get_counts_df(index)
        assert df.loc[0, "mrna_spliced"] == 5.0
        assert df.loc[0, "mrna"] == 15.0

    def test_identifiers_present(self):
        """Output includes transcript/gene identifiers."""
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        df = rc.get_counts_df(index)
        assert list(df["transcript_id"]) == ["t0", "t1", "t2"]
        assert list(df["gene_id"]) == ["g0", "g0", "g1"]


# =====================================================================
# Detail DataFrame output (long format)
# =====================================================================


class TestDetailOutput:
    def test_detail_empty(self):
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        df = rc.get_detail_df(index)
        assert len(df) == 0
        assert "transcript_id" in df.columns
        assert "category" in df.columns
        assert "source" in df.columns

    def test_detail_unambig_only(self):
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 5.0
        rc.unambig_counts[2, SpliceStrandCol.SPLICED_ANNOT_SENSE] = 3.0

        df = rc.get_detail_df(index)
        assert len(df) == 2
        assert set(df["source"]) == {"unambig"}

    def test_detail_both_sources(self):
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 5.0
        rc.em_counts[0, _UNSPLICED_SENSE] = 3.0

        df = rc.get_detail_df(index)
        assert len(df) == 2
        assert set(df["source"]) == {"unambig", "em"}
        assert df["count"].sum() == 8.0

    def test_detail_has_category(self):
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 1.0
        rc.unambig_counts[0, SpliceStrandCol.SPLICED_ANNOT_SENSE] = 1.0

        df = rc.get_detail_df(index)
        assert set(df["category"]) == {"unspliced", "spliced_annot"}


# =====================================================================
# gDNA in locus EM — unspliced compete with gDNA shadow
# =====================================================================


class TestGDNAInLocusEM:
    """Verify gDNA shadow competes with mRNA in locus EM."""

    def test_gdna_absorbs_when_init_high(self):
        """With large gdna_init and equal likelihoods, gDNA takes share."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data(
            [[0]] * 100,
            num_transcripts=2,
            include_gdna=True,
            locus_gamma=0.5,
            gdna_log_lik=0.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        # gDNA should absorb some fragments
        assert gdna_count > 0

    def test_strong_rna_beats_gdna(self):
        """When transcript likelihood >> gDNA, most go to transcript."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 500.0
        bundle = _make_locus_em_data(
            [[0]] * 100,
            log_liks_per_unit=[[0.0]] * 100,
            num_transcripts=2,
            rc=rc,
            include_gdna=True,
            locus_gamma=0.1,
            gdna_log_lik=-20.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert rc.em_counts[0].sum() > 90
        assert gdna_count < 10

    def test_total_counts_preserved_with_gdna(self):
        """em_counts + nrna_em + gdna == n_units."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data(
            [[0, 1]] * 200,
            num_transcripts=2,
            include_nrna=True,
            include_gdna=True,
            locus_gamma=0.3,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        total = rc.em_counts.sum() + rc.nrna_em_counts.sum() + gdna_count
        assert total == pytest.approx(200.0, abs=1.0)

    def test_no_gdna_candidate_means_no_gdna_assignment(self):
        """Without gDNA candidate, all counts go to RNA."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data(
            [[0, 1]] * 100,
            num_transcripts=2,
            include_gdna=False,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=5)
        gdna_count = pool_counts["gdna"]

        assert gdna_count == 0.0
        assert rc.em_counts.sum() == pytest.approx(100.0, abs=1.0)
