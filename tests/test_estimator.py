"""The locus EM — `AbundanceEstimator`, the native batch wrapper it drives, and the warm start that
seeds it. The first block runs the estimator over `ScoredFragments` and a mock index: locus assignment,
posterior means, simultaneous resolution, multimappers, the counts and detail tables, partitioned
effective lengths, the gDNA component and discrete assignment. The second calls
`run_batch_locus_em_partitioned` directly on the production path, gating structural gDNA eligibility,
the gDNA effective length, the aggregate RNA prior, and the rule that every RNA component receives its
share of that prior in proportion to its evidence. The third gates `EMConfig.warm_start`: what the initial theta is
derived from, and that an unknown name is refused.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.types import Strand
from rigel.splice import (
    SpliceStrandCol,
)
from rigel.config import EMConfig, TranscriptGeometry
from rigel.estimator import AbundanceEstimator
from rigel.scored_fragments import ScoredFragments

from _em_harness import _UNSPLICED_SENSE, _make_locus_em_data, _run_and_assign


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
                "g_type": ["protein_coding"] * num_transcripts,
                "ref": ["chr1"] * num_transcripts,
                "strand": t_to_strand,
                "start": list(range(0, num_transcripts * 1000, 1000)),
                "end": list(range(500, num_transcripts * 1000 + 500, 1000)),
                "length": [500] * num_transcripts,
                "is_basic": [True] * num_transcripts,
                "is_mane": [False] * num_transcripts,
                "is_nrna": [False] * num_transcripts,
                "is_synthetic": [False] * num_transcripts,
                "nrna_t_index": [-1] * num_transcripts,
                "nrna_n_contributors": [0] * num_transcripts,
            }
        )
        self.g_df = pd.DataFrame(
            {
                "g_id": g_ids,
                "g_name": g_names,
                "g_type": ["protein_coding"] * num_genes,
                "ref": ["chr1"] * num_genes,
                "strand": g_to_strand,
                "start": list(range(0, num_genes * 2000, 2000)),
                "end": list(range(1500, num_genes * 2000 + 1500, 2000)),
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


def _lengths(values):
    """One effective length per transcript, the same for the output and the EM (capture-off)."""
    a = np.asarray(values, dtype=np.float64)
    return TranscriptGeometry(effective_lengths=a, effective_lengths_em=a)


def _make_em_data(t_indices_per_unit):
    """Build ScoredFragments from a list of per-unit candidate lists, every candidate an unspliced
    sense hit at log-likelihood 0.

    The global ScoredFragments contains mRNA + nRNA candidates only (no gDNA).
    """
    offsets = [0]
    flat_t = []
    for t_list in t_indices_per_unit:
        flat_t.extend(t_list)
        offsets.append(len(flat_t))

    n_units = len(t_indices_per_unit)
    n_candidates = len(flat_t)

    # Build locus tracking arrays: each unit's first candidate
    locus_t = np.full(n_units, -1, dtype=np.int32)
    locus_cc = np.zeros(n_units, dtype=np.uint8)
    for u, t_list in enumerate(t_indices_per_unit):
        if t_list:
            locus_t[u] = t_list[0]
            locus_cc[u] = _UNSPLICED_SENSE

    return ScoredFragments(
        offsets=np.array(offsets, dtype=np.int64),
        t_indices=np.array(flat_t, dtype=np.int32),
        log_liks=np.zeros(n_candidates, dtype=np.float64),
        count_cols=np.full(n_candidates, _UNSPLICED_SENSE, dtype=np.uint8),
        coverage_weights=np.ones(n_candidates, dtype=np.float64),
        locus_t_indices=locus_t,
        locus_count_cols=locus_cc,
        is_spliced=np.zeros(n_units, dtype=bool),
        gdna_log_liks=np.full(n_units, -np.inf, dtype=np.float64),
        frag_ids=np.arange(n_units, dtype=np.int64),
        frag_class=np.zeros(n_units, dtype=np.int8),
        splice_type=np.zeros(n_units, dtype=np.uint8),
        n_units=n_units,
        n_candidates=n_candidates,
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
# Locus assignment — the default sampled draw, and fractional where a test sets it
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

        # 10 unambig, set directly
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
        """Over many fragments, sampled assignments follow the unambig priors."""
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
        # Both should receive counts (not hard-gated to one)
        assert t1_em > 500  # ~10% of 10000

    def test_writes_to_em_not_unambig(self):
        """Batch EM only writes to em_counts."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data([[0, 1]] * 50, num_transcripts=3)

        unique_before = rc.unambig_counts.copy()
        _run_and_assign(rc, bundle, em_iterations=5)

        np.testing.assert_array_equal(rc.unambig_counts, unique_before)
        assert rc.em_counts.sum() == pytest.approx(50.0, abs=1.0)

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
        """The sampled assignment is deterministic for a fixed seed."""
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
        """Different seeds do not affect fractional (expected-count) assignment."""
        counts = []
        for seed in [1, 2]:
            rc = AbundanceEstimator(3, em_config=EMConfig(seed=seed, assignment_mode="fractional"))
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
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42, assignment_mode="fractional"))
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
        rc1 = AbundanceEstimator(4, em_config=EMConfig(seed=42, assignment_mode="fractional"))
        rc1.unambig_counts[0, 0] = 50.0
        rc1.unambig_counts[2, 0] = 30.0

        units = [[0, 1]] * 200 + [[0, 2]] * 200
        bundle1 = _make_locus_em_data(units, num_transcripts=4, rc=rc1)
        _run_and_assign(rc1, bundle1, em_iterations=10)

        rc2 = AbundanceEstimator(4, em_config=EMConfig(seed=42, assignment_mode="fractional"))
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
            "gene_type",
            "ref",
            "strand",
            "start",
            "end",
            "length",
            "effective_length",
            "em_effective_length",
            "locus_id",
            "nrna_id",
            "is_basic",
            "is_mane",
            "is_nrna",
            "count",
            "count_unambig",
            "count_em",
            "count_spliced",
            "nrna_parent_count",
            "tpm",
            "tpm_total_rna",
            "posterior_mean",
        ]
        assert list(df.columns) == expected_cols

    def test_counts_df_keeps_raw_and_em_effective_lengths_separate(self):
        index = _make_index()
        geometry = TranscriptGeometry(
            effective_lengths=np.array([100.0, 200.0, 300.0]),
            effective_lengths_em=np.array([50.0, 100.0, 300.0]),
        )
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42), geometry=geometry)

        df = rc.get_counts_df(index)
        assert df.loc[0, "effective_length"] == pytest.approx(100.0)
        assert df.loc[0, "em_effective_length"] == pytest.approx(50.0)

    def test_get_gene_counts_df_columns(self):
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        df = rc.get_gene_counts_df(index)
        assert df.shape[0] == 2
        expected_cols = [
            "gene_id",
            "gene_name",
            "gene_type",
            "ref",
            "strand",
            "start",
            "end",
            "n_transcripts",
            "locus_id",
            "effective_length",
            "count",
            "mature_count",
            "nascent_count",
            "count_unambig",
            "count_em",
            "count_spliced",
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
        assert df.loc[0, "count"] == 8.0
        assert df.loc[0, "count_unambig"] == 5.0
        assert df.loc[0, "count_em"] == 3.0

    def test_gene_counts_aggregate(self):
        """Gene counts aggregate across transcripts."""
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 5.0
        rc.unambig_counts[1, _UNSPLICED_SENSE] = 3.0  # same gene
        rc.unambig_counts[2, _UNSPLICED_SENSE] = 7.0  # different gene

        df = rc.get_gene_counts_df(index)
        assert df.loc[0, "count"] == 8.0  # g0 = t0 + t1
        assert df.loc[1, "count"] == 7.0  # g1 = t2

    def test_gene_effective_length_abundance_weighted(self):
        """Gene effective length is abundance-weighted mean of transcript eff lens."""
        index = _make_index()  # t0,t1 → g0; t2 → g1
        rc = AbundanceEstimator(
            3, em_config=EMConfig(seed=42), geometry=_lengths([500.0, 1000.0, 300.0])
        )
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
        rc = AbundanceEstimator(
            3, em_config=EMConfig(seed=42), geometry=_lengths([500.0, 1000.0, 300.0])
        )
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
        assert df.loc[0, "count_spliced"] == 5.0
        assert df.loc[0, "count"] == 15.0

    def test_identifiers_present(self):
        """Output includes transcript/gene identifiers."""
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        df = rc.get_counts_df(index)
        assert list(df["transcript_id"]) == ["t0", "t1", "t2"]
        assert list(df["gene_id"]) == ["g0", "g0", "g1"]


class TestPartitionedEffectiveLength:
    """Production partitioned EM must consume per-transcript L̃ values."""

    def test_partitioned_em_uses_nontrivial_effective_lengths(self):
        units = [[0, 1]] * 100

        rc_equal = AbundanceEstimator(
            2,
            em_config=EMConfig(seed=42, assignment_mode="fractional"),
        )
        _run_and_assign(
            rc_equal,
            _make_locus_em_data(units, num_transcripts=2),
            em_iterations=100,
        )
        equal_counts = rc_equal.em_counts.sum(axis=1)
        np.testing.assert_allclose(equal_counts[0], equal_counts[1], rtol=1e-3, atol=1e-3)

        rc_len = AbundanceEstimator(
            2,
            em_config=EMConfig(seed=42, assignment_mode="fractional"),
            geometry=_lengths([10.0, 1.0]),
        )
        _run_and_assign(
            rc_len,
            _make_locus_em_data(units, num_transcripts=2),
            em_iterations=100,
        )
        length_counts = rc_len.em_counts.sum(axis=1)

        assert length_counts[1] > length_counts[0] * 5.0


# =====================================================================
# gDNA in locus EM — unspliced compete with gDNA shadow
# =====================================================================


class TestGDNAInLocusEM:
    """Verify gDNA shadow competes with mRNA in locus EM."""

    def test_gdna_absorbs_when_init_high(self):
        """With a gDNA prior count of 5 and equal likelihoods, gDNA takes a share."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42, assignment_mode="fractional"))
        bundle = _make_locus_em_data(
            [[0]] * 100,
            num_transcripts=2,
            include_gdna=True,
            gdna_prior_count=5.0,
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
            gdna_prior_count=1.0,
            gdna_log_lik=-20.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert rc.em_counts[0].sum() > 90
        assert gdna_count < 10

    def test_total_counts_preserved_with_gdna(self):
        """em_counts + gdna == n_units."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data(
            [[0, 1]] * 200,
            num_transcripts=2,
            include_nrna=True,
            include_gdna=True,
            gdna_prior_count=3.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        total = rc.em_counts.sum() + gdna_count
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


# =====================================================================
# Discrete assignment modes (map / sample)
# =====================================================================


class TestDiscreteAssignment:
    """Test post-EM discrete fragment assignment modes."""

    def test_map_mode_produces_integer_counts(self):
        """MAP assignment assigns each fragment to exactly one component."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42, assignment_mode="map"))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 50.0
        rc.unambig_counts[1, _UNSPLICED_SENSE] = 30.0
        bundle = _make_locus_em_data(
            [[0, 1]] * 200,
            num_transcripts=3,
            rc=rc,
        )
        _run_and_assign(rc, bundle, em_iterations=10)
        # All EM counts should be integers
        frac = rc.em_counts - np.floor(rc.em_counts)
        assert not np.any(frac > 0.0), "MAP mode should produce only integer counts"
        assert rc.em_counts.sum() == pytest.approx(200.0, abs=1.0)

    def test_sample_mode_produces_integer_counts(self):
        """Sample assignment assigns each fragment to exactly one component."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42, assignment_mode="sample"))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 50.0
        rc.unambig_counts[1, _UNSPLICED_SENSE] = 30.0
        bundle = _make_locus_em_data(
            [[0, 1]] * 200,
            num_transcripts=3,
            rc=rc,
        )
        _run_and_assign(rc, bundle, em_iterations=10)
        frac = rc.em_counts - np.floor(rc.em_counts)
        assert not np.any(frac > 0.0), "Sample mode should produce only integer counts"
        assert rc.em_counts.sum() == pytest.approx(200.0, abs=1.0)

    def test_sample_mode_deterministic_with_same_seed(self):
        """Same seed → same sample assignment results."""
        results = []
        for _ in range(2):
            rc = AbundanceEstimator(3, em_config=EMConfig(seed=42, assignment_mode="sample"))
            rc.unambig_counts[0, _UNSPLICED_SENSE] = 50.0
            rc.unambig_counts[1, _UNSPLICED_SENSE] = 30.0
            bundle = _make_locus_em_data(
                [[0, 1]] * 200,
                num_transcripts=3,
                rc=rc,
            )
            _run_and_assign(rc, bundle, em_iterations=10)
            results.append(rc.em_counts.copy())
        np.testing.assert_array_equal(results[0], results[1])

    def test_sample_mode_different_seeds_differ(self):
        """Different seeds → different sample assignment results."""
        results = []
        for seed in [42, 99]:
            rc = AbundanceEstimator(3, em_config=EMConfig(seed=seed, assignment_mode="sample"))
            rc.unambig_counts[0, _UNSPLICED_SENSE] = 50.0
            rc.unambig_counts[1, _UNSPLICED_SENSE] = 30.0
            bundle = _make_locus_em_data(
                [[0, 1]] * 200,
                num_transcripts=3,
                rc=rc,
            )
            _run_and_assign(rc, bundle, em_iterations=10)
            results.append(rc.em_counts.copy())
        # With different seeds, counts should differ (probabilistic but overwhelming)
        assert not np.array_equal(results[0], results[1])

    def test_map_mode_winner_matches_highest_posterior(self):
        """MAP assigns all ambiguous fragments to the component with higher prior."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42, assignment_mode="map"))
        # Transcript 0 has much more unambiguous support
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 500.0
        rc.unambig_counts[1, _UNSPLICED_SENSE] = 10.0
        bundle = _make_locus_em_data(
            [[0, 1]] * 100,
            num_transcripts=3,
            rc=rc,
        )
        _run_and_assign(rc, bundle, em_iterations=10)
        # In MAP mode, the stronger transcript should win all or most
        assert rc.em_counts[0].sum() >= 90, "MAP should favor the transcript with higher posterior"

    def test_fractional_mode_preserves_posteriors(self):
        """Fractional mode produces non-integer counts (original behavior)."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42, assignment_mode="fractional"))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 50.0
        rc.unambig_counts[1, _UNSPLICED_SENSE] = 30.0
        bundle = _make_locus_em_data(
            [[0, 1]] * 200,
            num_transcripts=3,
            rc=rc,
        )
        _run_and_assign(rc, bundle, em_iterations=10)
        frac = rc.em_counts - np.floor(rc.em_counts)
        assert np.any(frac > 0.0), "Fractional mode should produce non-integer counts"

    def test_total_counts_preserved_all_modes(self):
        """All three modes preserve total fragment count."""
        for mode in ["fractional", "map", "sample"]:
            rc = AbundanceEstimator(3, em_config=EMConfig(seed=42, assignment_mode=mode))
            rc.unambig_counts[0, _UNSPLICED_SENSE] = 50.0
            bundle = _make_locus_em_data(
                [[0, 1]] * 300,
                num_transcripts=3,
                rc=rc,
            )
            _run_and_assign(rc, bundle, em_iterations=10)
            total = rc.em_counts.sum()
            assert total == pytest.approx(300.0, abs=1.0), f"mode={mode} lost fragments"

    def test_min_posterior_threshold(self):
        """Components below min_posterior threshold are excluded from assignment."""
        # With a high min_posterior, only the dominant component should win
        rc = AbundanceEstimator(
            4, em_config=EMConfig(seed=42, assignment_mode="map", assignment_min_posterior=0.3)
        )
        # Strongly favor transcript 0
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 500.0
        bundle = _make_locus_em_data(
            [[0, 1, 2, 3]] * 200,
            num_transcripts=4,
            rc=rc,
        )
        _run_and_assign(rc, bundle, em_iterations=10)
        # Transcript 0 should get all 200 (others below threshold)
        assert rc.em_counts[0].sum() >= 190


# ── The native batch locus EM wrapper, on the production path ────────────────────────────────────


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
    return est


def test_equal_rna_likelihoods_split_evenly():
    n_units = 100
    est = _estimator(2, mode="map")
    total_gdna, locus_rna, locus_gdna = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(0.0, 0.0))],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
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
        gdna_eff_len=np.array([1.0], dtype=np.float64),
    )

    est_long = _estimator(1, mode="map")
    _total_long, _rna_long, gdna_long = est_long.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        gdna_eff_len=np.array([100.0], dtype=np.float64),
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
        gdna_eff_len=np.array([1.0], dtype=np.float64),
    )

    est_rna = _estimator(1, mode="map")
    _total_rna, _rna_rna, gdna_rna = est_rna.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        rna_prior_count=np.array([10.0], dtype=np.float64),
        gdna_eff_len=np.array([1.0], dtype=np.float64),
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
        emit_assignments=True,
    )

    total_gdna, _locus_rna, _locus_gdna, winner_tid, winner_post, n_candidates = result
    assert total_gdna == pytest.approx(0.0)
    np.testing.assert_array_equal(winner_tid, np.zeros(n_units, dtype=np.int32))
    assert np.all(winner_post > 0.98)
    np.testing.assert_array_equal(n_candidates, np.full(n_units, 2, dtype=np.int16))


# ---------------------------------------------------------------------------
# Structural gDNA eligibility: native derives candidate availability from the
# partition itself, and ``gdna_prior_count`` does not gate the gDNA component.
# ---------------------------------------------------------------------------


def test_gdna_candidates_are_derived_from_the_partition():
    """A locus admits gDNA iff some unspliced unit carries a finite gDNA log-lik, whatever the
    prior count."""
    n_units = 20
    locus_t_lists = [np.array([0, 1], dtype=np.int32)]

    # All-spliced partition ⇒ no per-unit gDNA candidate.
    spliced_part = _partition(
        n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-0.5, is_spliced=True
    )
    # Unspliced partition with finite gDNA log-liks ⇒ candidates.
    unspliced_part = _partition(
        n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-0.5, is_spliced=False
    )
    # Unspliced but non-finite gDNA log-lik ⇒ no candidate.
    nogdna_part = _partition(
        n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-np.inf, is_spliced=False
    )

    # All three with gdna_prior_count=0: only the unspliced+finite case produces gDNA assignments.
    g_spl, _, _ = _estimator(2, mode="map").run_batch_locus_em_partitioned(
        [spliced_part],
        locus_t_lists,
        np.zeros(1),
    )
    g_uns, _, _ = _estimator(2, mode="map").run_batch_locus_em_partitioned(
        [unspliced_part],
        locus_t_lists,
        np.zeros(1),
    )
    g_no, _, _ = _estimator(2, mode="map").run_batch_locus_em_partitioned(
        [nogdna_part],
        locus_t_lists,
        np.zeros(1),
    )

    assert g_spl == 0.0, "spliced partition has no gDNA candidates"
    assert g_no == 0.0, "non-finite gDNA log-liks ⇒ no gDNA candidates"
    assert g_uns > 0.0, (
        "unspliced+finite gDNA log-lik must admit the component even when gdna_prior_count == 0"
    )


def test_positive_gdna_prior_produces_finite_outputs():
    """A positive gDNA prior count produces finite, conserved outputs."""
    n_units = 30
    est = _estimator(2, mode="map")
    total_gdna, locus_rna, locus_gdna = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-1.0)],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([5.0], dtype=np.float64),
    )

    # Outputs must be finite and sum to total fragments.
    assert np.isfinite(total_gdna)
    assert np.isfinite(locus_rna).all() and np.isfinite(locus_gdna).all()
    assert total_gdna >= 0.0
    # Conservation up to assignment fractional remainder.
    assigned = float(locus_rna[0]) + float(locus_gdna[0])
    assert assigned == pytest.approx(30.0, rel=1e-9)


# ── EVERY RNA COMPONENT GETS ITS SHARE OF THE RNA PRIOR ──────────────────────────────────────────
#
# The RNA pseudocount is distributed over the locus's RNA components in proportion to the evidence
# each already carries, and NO component is singled out for zero (`EQUATIONS.md` §9b). RNA is RNA:
# whether the annotation happens to assert a given RNA component is not a fact about this locus's
# composition, so the allocation does not read it.
#
# Withholding the share from any component — SYNTHETIC nascent entities, say, the spans the index
# manufactured — would make the prior's factor un-common over the pool, so the prior ALONE would
# redistribute RNA between entities the data cannot tell apart. These gates pin the rule and the one
# guard it keeps: a component with no evidence cannot be revived by prior mass, because the weights
# are the evidence.


def _run(est, partition, t_idx, *, rna_prior=0.0, gdna_prior=0.0):
    return est.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.asarray(t_idx, dtype=np.int32)],
        gdna_prior_count=np.array([gdna_prior], dtype=np.float64),
        rna_prior_count=np.array([rna_prior], dtype=np.float64),
    )


def _mixed_partition(unique: tuple[int, ...], shared: int, gdna_log_lik: float = -np.inf):
    """One locus: ``unique[t]`` fragments compatible with transcript ``t`` alone, then ``shared``
    fragments compatible with every transcript, all at log-likelihood 0."""
    cands = [[t] for t, n in enumerate(unique) for _ in range(n)]
    cands += [list(range(len(unique)))] * shared
    offsets = np.cumsum([0] + [len(c) for c in cands]).astype(np.int64)
    flat = np.concatenate([np.asarray(c, dtype=np.int32) for c in cands])
    n_units, n_cand = len(cands), flat.shape[0]
    return (
        offsets,
        flat,
        np.zeros(n_cand, dtype=np.float64),
        np.ones(n_cand, dtype=np.float64),
        np.zeros(n_cand, dtype=np.uint8),
        np.zeros(n_units, dtype=np.uint8),
        np.full(n_units, gdna_log_lik, dtype=np.float64),
        np.zeros(n_units, dtype=np.int32),
        np.zeros(n_units, dtype=np.uint8),
    )


def _yield_estimator(yields, *, mode: str) -> AbundanceEstimator:
    return AbundanceEstimator(
        num_transcripts=len(yields),
        em_config=EMConfig(
            mode=mode, iterations=200, convergence_delta=1e-9, assignment_mode="fractional", seed=0
        ),
        geometry=_lengths(yields),
    )


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_the_split_is_invariant_to_a_common_thinning_of_every_yield(mode):
    """Capture thins the yield of every component in a locus; thinned by ONE factor, no count may move,
    because the E-step reads each component's count against its yield and only the ratios within the
    locus decide. A floor on the yield breaks it: clamping a 300 bp yield thinned a thousandfold to
    1 bp moves the split of the shared fragments."""

    def counts(scale):
        est = _yield_estimator([300.0 * scale, 3000.0 * scale], mode=mode)
        _run(est, _mixed_partition((100, 100), 200), [0, 1])
        return est.em_counts.sum(axis=1)

    full, thinned = counts(1.0), counts(1e-3)
    assert full[0] > 100.0 and full[1] > 100.0, full  # the shared fragments are split, not dumped
    np.testing.assert_allclose(thinned, full, rtol=1e-9)


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_transcript_with_no_start_position_cannot_emit(mode):
    """A transcript shorter than every fragment has a yield of exactly 0 and produced nothing: its share
    of every fragment is 0 — not the share of a floored 1 bp yield, which would make it the densest
    component in the locus and hand it every shared fragment. The fragments go to the components that can emit; a
    fragment no component can emit is left unassigned, and nothing is NaN."""
    est = _yield_estimator([0.0, 1000.0], mode=mode)
    est.run_batch_locus_em_partitioned(
        partition_tuples=[_mixed_partition((100, 0), 200)],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        rna_prior_count=np.array([0.0], dtype=np.float64),
        emit_locus_stats=True,
    )
    np.testing.assert_allclose(est.em_counts.sum(axis=1), [0.0, 200.0], atol=1e-9)
    assert np.isfinite(
        est.locus_stats[0]["final_data_loglik"]
    )  # the unemittable rows are outside it

    est = _yield_estimator([0.0], mode=mode)
    _t, rna, gdna = est.run_batch_locus_em_partitioned(
        partition_tuples=[_mixed_partition((0,), 100, gdna_log_lik=0.0)],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        gdna_eff_len=np.array([100.0], dtype=np.float64),
    )
    assert gdna[0] == pytest.approx(100.0) and rna[0] == pytest.approx(0.0)

    est = _yield_estimator([0.0], mode=mode)
    _t, rna, gdna = _run(est, _mixed_partition((0,), 100), [0])
    assert rna[0] == 0.0 and gdna[0] == 0.0
    assert np.all(np.isfinite(est.em_counts)) and est.em_counts.sum() == 0.0


_SHAPES = {
    # every component holds fragments only IT can explain, so none is pruned and the within-RNA
    # split is genuinely at stake
    "tied2": ((100, 100), 0),
    "shared2": ((60, 60), 80),
    "skewed3": ((60, 40, 20), 80),
    "wide4": ((50, 40, 30, 20), 60),
}


def _counts_against_prior(shape, *, mode, priors=(0.0, 500.0)):
    """The same locus solved at two RNA-prior magnitudes. No gDNA candidate anywhere
    (`_mixed_partition`'s default `gdna_log_lik = -inf`), so `theta` normalises over the RNA
    components alone and the prior's only possible effect is the one under test."""
    unique, shared = shape
    part = _mixed_partition(unique, shared)
    out = []
    for rna_prior in priors:
        est = _estimator(len(unique), mode=mode)
        _run(est, part, list(range(len(unique))), rna_prior=rna_prior)
        out.append(est.em_counts.sum(axis=1).copy())
    return out


@pytest.mark.parametrize("shape", list(_SHAPES), ids=list(_SHAPES))
def test_the_RNA_prior_moves_NO_component_s_SHARE_of_the_RNA_pool_under_MAP(shape):
    """The rule's promise, stated where it is EXACTLY true.

    The RNA pseudocount reaches every RNA component as the same factor `(1 + rna_prior/rna_count)`,
    and under MAP `theta` is proportional to those counts — so the factor cancels and no prior
    magnitude whatsoever can move the answer. That is what "distributed over the RNA components, none
    singled out for zero" MEANS, and it is exactly what a per-component eligibility test breaks: hold
    one component out and the factor stops being common, so the prior ALONE redistributes RNA between
    entities the data cannot tell apart.

    An eligibility test does not merely tilt the split: the withheld factor compounds once per M-step,
    so on a tied two-component locus the held-out component loses essentially all of its mass.

    ⛔ The equality is EXACT, not approximate. One common multiply leaves `theta` bit for bit where it
    was; anything that reads a component would not.
    """
    none, large = _counts_against_prior(_SHAPES[shape], mode="map")
    np.testing.assert_array_equal(
        large, none, err_msg=f"the RNA prior redistributed the RNA pool: {none} -> {large}"
    )


@pytest.mark.parametrize("shape", list(_SHAPES), ids=list(_SHAPES))
def test_under_VBEM_the_prior_moves_the_split_only_by_the_DIGAMMA_CORRECTION(shape):
    """The same claim under the shipped mode, where it is exact only in the limit — and the bound is
    DERIVED rather than chosen.

    VBEM's M-step is `theta_i ∝ exp(psi(alpha_i))`, which is NOT scale-equivariant, so a common factor
    `c = 1 + rna_prior/rna_count` on every RNA alpha does not cancel the way it does under MAP. With
    `psi(x) = log x − 1/(2x) + O(x^-2)`,

        psi(c·alpha_i) − psi(alpha_i) = log c + (1 − 1/c)/(2·alpha_i) + O(alpha_i^-2)

    and `log c` IS common, so it normalises away. What is left is a per-component residual
    `(1 − 1/c)/(2·alpha_i)`, largest at the smallest alpha and bounded by `1/(2·alpha_min)` since
    `c ≥ 1`. So the uniformity of the allocation survives; only the M-step's own nonlinearity moves
    the shares, by an amount that vanishes as the locus deepens.

    ⭐ It separates this rule from an eligibility test by orders of magnitude, which is what a gate is
    for: on `skewed3` the residual moves a share by 1.4e-3 against a bound of 1.5e-2, while an
    eligibility test moves a component by 100 % of its mass.
    """
    none, large = _counts_against_prior(_SHAPES[shape], mode="vbem")
    share_none, share_large = none / none.sum(), large / large.sum()
    bound = 1.0 / (2.0 * none.min())
    moved = np.abs(share_large / share_none - 1.0)
    assert moved.max() <= bound, (
        f"the prior moved a share by {moved.max():.3e}, beyond the digamma residual's "
        f"bound 1/(2·alpha_min) = {bound:.3e}: {none} -> {large}"
    )


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_SHADOW_SPAN_STILL_LOSES_to_the_transcript_it_shadows(mode):
    """What guards against a zombie.

    The RNA prior reaches every RNA component as the same factor, so it does nothing to make a shadow
    decay; the decay rests on the LIKELIHOOD alone, and the likelihood is sufficient. A shadow span is longer than the
    transcript it shadows, so at equal per-fragment likelihood its yield ratio `kappa = w_N/w_T` is
    strictly below 1 and a component with no evidence of its own decays geometrically at that rate
    (`EQUATIONS.md` §9b).

    Component 1 here explains nothing component 0 does not, and its yield is 10× longer. It must lose
    the shared mass even with a large RNA prior in play.
    """
    est = _yield_estimator([1000.0, 10000.0], mode=mode)
    _run(
        est,
        _partition(n_units=200, log_liks=(0.0, 0.0), gdna_log_lik=-2.0),
        [0, 1],
        rna_prior=100.0,
        gdna_prior=10.0,
    )
    counts = est.em_counts.sum(axis=1)
    assert counts[1] < 0.01 * counts.sum(), (
        f"a shadow span with no evidence of its own kept mass: {counts}"
    )


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_component_the_DATA_SUPPORTS_keeps_its_mass(mode):
    """The gate that stops the one above becoming "kill every long component". A rule that zeroes
    anything with a long yield passes it; only this separates the two. Component 1 is the strictly
    better explanation here — a 10× longer yield against a likelihood advantage of `e^8` — and it
    must keep the mass."""
    est = _yield_estimator([1000.0, 10000.0], mode=mode)
    _run(
        est,
        _partition(n_units=200, log_liks=(-8.0, 0.0), gdna_log_lik=-8.0),
        [0, 1],
        rna_prior=100.0,
        gdna_prior=10.0,
    )
    counts = est.em_counts.sum(axis=1)
    assert counts[1] > 0.9 * counts.sum(), (
        f"a component with decisive likelihood support was suppressed anyway: {counts}"
    )


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_ZERO_EVIDENCE_component_is_NOT_revived_by_the_prior(mode):
    """The guard the rule keeps, and the reason the weights are the evidence rather than a
    flat share (`EQUATIONS.md` §9b.1). `out[i]` is proportional to `raw[i]`, so `out[i] = 0` is an
    ABSORBING STATE that no prior magnitude escapes — and the RNA pseudocount is a FRAGMENT COUNT,
    tens to thousands on an expressed locus, so a flat share would hand every component far more than
    the ~0.16–0.47 alpha units at which one activates.

    Component 1 is unreachable: it is a candidate for no unit at all, so its warm start is 0.
    """
    part = _mixed_partition((200, 0), 0, gdna_log_lik=-2.0)
    est = _yield_estimator([1000.0, 1000.0], mode=mode)
    _run(est, part, [0, 1], rna_prior=5000.0, gdna_prior=10.0)
    counts = est.em_counts.sum(axis=1)
    assert counts[1] == 0.0, f"a component with no evidence was revived by prior mass: {counts}"


# The per-M-step identity this design rests on is gated in C++, not here.
#
# `apply_grouped_prior_update` guarantees, FOR A GIVEN raw_counts VECTOR, that the RNA components sum
# to `rna_count + rna_prior`, so the RNA prior is redistributed strictly WITHIN the RNA pool and the
# library gDNA fraction cannot move by this rule. The function is `static` in `em_solver.cpp`; the
# test-only binding `_apply_grouped_prior_update_test` reaches it and
# `tests/native/test_grouped_prior_update.py` holds the identity.
#
# ⚠ It is a PER-M-STEP algebraic identity, NOT an end-to-end one, and conflating the two produces a
# wrong test. End to end the gDNA total legitimately DEPENDS on `rna_prior`: setting the
# gDNA:RNA split is the prior's whole purpose (`gdna_total = gdna_count + gdna_prior`,
# `rna_total = rna_count + rna_prior`), and a larger RNA prior shifts theta, hence the E-step, hence
# the next iteration's `gdna_count`. A test asserting "the library gDNA fraction is invariant to
# rna_prior" is asserting something FALSE BY DESIGN.


# ── ``EMConfig.warm_start`` — what the EM's initial ``theta`` is derived FROM ─────────────────────
#
# The shipped warm start seeds ``theta`` with each component's unambiguous total plus a
# coverage-weighted share of the ambiguous fragments, and then projects that seed through the
# calibration prior. Seed and prior are two different methods, and the projection multiplies them — so
# a coverage-weighted share scaled by a per-transcript allocation derived some other way is neither
# method's answer. ``warm_start="prior"`` zeroes the seed, so a prior can be tested on its own terms.
# It is only meaningful with a per-transcript weight, and that is a property rather than a caveat:
# under the shipped evidence-proportional rule ``out[i]`` is proportional to ``raw[i]``, so an all-zero
# seed yields an all-zero RNA pool and the whole locus goes to gDNA — ``theta = 0`` is an absorbing
# state. The tests below pin both halves, because the failure is silent: a starved locus still
# converges. The helpers carry a ``_warm_start_`` prefix because the block above defines its own
# ``_partition`` / ``_estimator`` / ``_run`` with different defaults.


def _warm_start_partition(
    *, n_units: int, log_liks: tuple[float, ...], gdna_log_lik: float = -50.0
):
    """``gdna_log_lik`` defaults to a FINITE (but hopeless) value on purpose. With ``-inf`` no unit
    carries a gDNA candidate, and ``apply_grouped_prior_update`` then discards the RNA prior along with
    the gDNA one — so a fixture built that way measures a prior that was never applied."""
    n_t = len(log_liks)
    offsets = np.arange(n_units + 1, dtype=np.int64) * n_t
    return (
        offsets,
        np.tile(np.arange(n_t, dtype=np.int32), n_units),
        np.tile(np.asarray(log_liks, dtype=np.float64), n_units),
        np.ones(n_units * n_t, dtype=np.float64),
        np.zeros(n_units * n_t, dtype=np.uint8),
        np.zeros(n_units, dtype=np.uint8),
        np.full(n_units, gdna_log_lik, dtype=np.float64),
        np.zeros(n_units, dtype=np.int32),
        np.zeros(n_units, dtype=np.uint8),
    )


def _warm_start_estimator(n_t: int, *, warm_start: str = "coverage", mode: str = "vbem"):
    est = AbundanceEstimator(
        num_transcripts=n_t,
        em_config=EMConfig(
            mode=mode,
            iterations=500,
            convergence_delta=1e-10,
            assignment_mode="fractional",
            seed=0,
            warm_start=warm_start,
        ),
    )
    return est


def _warm_start_run(
    est,
    partition,
    t_idx,
    *,
    weight=None,
    rna_prior=0.0,
    gdna_prior=0.0,
    iterations=500,
):
    # `em_iterations` is an EXPLICIT argument of the estimator entry point, NOT read from
    # `em_config` — a helper that sets it on the config alone silently runs the default 1000.
    est.run_batch_locus_em_partitioned(
        em_iterations=iterations,
        partition_tuples=[partition],
        locus_transcript_indices=[np.asarray(t_idx, dtype=np.int32)],
        gdna_prior_count=np.array([gdna_prior], dtype=np.float64),
        rna_prior_count=np.array([rna_prior], dtype=np.float64),
        rna_prior_weight=None if weight is None else np.asarray(weight, dtype=np.float64),
    )
    return est.em_counts.sum(axis=1)


# ──────────────────────────────────────────────────────────────────────────────
# The default is untouched
# ──────────────────────────────────────────────────────────────────────────────


def test_the_SHIPPED_default_is_coverage_and_is_BYTE_IDENTICAL_to_no_switch():
    """The control. A config field that changes the default is not a switch, it is a release."""
    assert EMConfig().warm_start == "coverage"
    part = _warm_start_partition(n_units=100, log_liks=(0.0, -0.4), gdna_log_lik=-1.0)
    a = _warm_start_run(_warm_start_estimator(2), part, [0, 1], rna_prior=30.0, gdna_prior=10.0)
    b = _warm_start_run(
        _warm_start_estimator(2, warm_start="coverage"),
        part,
        [0, 1],
        rna_prior=30.0,
        gdna_prior=10.0,
    )
    np.testing.assert_array_equal(a, b)


# ──────────────────────────────────────────────────────────────────────────────
# theta starts from the prior alone
# ──────────────────────────────────────────────────────────────────────────────


def test_the_WEIGHT_moves_the_converged_split_into_the_weights_ratio():
    """Two components whose likelihoods are IDENTICAL, so the data cannot separate them and the only
    thing saying where mass belongs is the prior. At the fixed point ``theta_i`` is proportional to
    ``raw_i + a_i``; with flat likelihoods ``raw_i`` is itself proportional to ``theta_i``, so the
    algebra collapses to ``theta_i = a_i / Σa`` — the WEIGHTS' ratio."""
    part = _warm_start_partition(n_units=100, log_liks=(0.0, 0.0))
    flat = _warm_start_run(_warm_start_estimator(2), part, [0, 1], rna_prior=40.0)
    tilted = _warm_start_run(
        _warm_start_estimator(2), part, [0, 1], weight=[3.0, 1.0], rna_prior=40.0
    )

    assert flat[0] == pytest.approx(flat[1], rel=1e-6), (
        "flat likelihoods must split evenly with no weight"
    )
    assert tilted[0] / tilted[1] == pytest.approx(3.0, rel=0.05), (
        "the weight did not move the converged split into its own ratio"
    )


def test_the_split_is_MONOTONE_in_the_weight_ratio():
    """A PROPERTY rather than a number. The exact converged ratio is perturbed by the gDNA component
    sharing the locus, so pinning it to `w0/w1` across a sweep would be fitting the fixture. What the
    allocation must guarantee is direction and order: equal weights split evenly, and a heavier weight
    always takes strictly more."""
    part = _warm_start_partition(n_units=200, log_liks=(0.0, 0.0))
    ratios = [
        _warm_start_run(_warm_start_estimator(2), part, [0, 1], weight=[w, 1.0], rna_prior=25.0)
        for w in (0.25, 1.0, 3.0, 9.0)
    ]
    got = [float(o[0] / o[1]) for o in ratios]
    assert got[1] == pytest.approx(1.0, rel=1e-6), "equal weights did not split evenly"
    assert got == sorted(got), f"the split is not monotone in the weight ratio: {got}"
    assert got[0] < 1.0 < got[2] < got[3]


def test_the_prior_only_warm_start_CHANGES_THE_SEED_but_not_a_UNIQUE_fixed_point():
    """What zeroing the warm start actually does: it changes where the EM STARTS, not where it ends.
    When the fixed point is unique the two arms converge to the same answer, which is correct
    behaviour rather than evidence the switch is inert.

    So the switch is demonstrated on the ITERATION, not on the limit: truncated to a single pass the
    two seeds give visibly different answers, and at convergence they agree. A test that only compared
    converged values could not tell a live switch from a dead one
    (``TRAPS: could-the-arm-have-fired``).
    """
    part = _warm_start_partition(n_units=100, log_liks=(0.0, -0.7))
    kw = dict(weight=[3.0, 1.0], rna_prior=40.0, gdna_prior=5.0)

    def _one(ws, iters):
        return _warm_start_run(
            _warm_start_estimator(2, warm_start=ws), part, [0, 1], iterations=iters, **kw
        )

    early_cov, early_pri = _one("coverage", 1), _one("prior", 1)
    assert not np.allclose(early_cov, early_pri), "the warm-start switch never reached the solver"

    late_cov, late_pri = _one("coverage", 500), _one("prior", 500)
    np.testing.assert_allclose(late_cov, late_pri, rtol=1e-6)


def test_the_RNA_PRIOR_REACHES_a_locus_with_NO_gDNA_CANDIDATE():
    """Reading BOTH pseudocounts through ``has_gdna`` discards the whole RNA prior::

        gdna_prior = has_gdna ? ... : 0.0
        rna_prior  = has_gdna ? ... : 0.0     // <- discards the RNA prior too

    A locus with no gDNA candidate is one whose fragments are ALL SPLICED — a gDNA candidate is
    appended to every unspliced unit — so that withholds the RNA prior from exactly the loci whose RNA
    is most certain.

    It is invisible under the shipped weights, which make it a no-op: the evidence-proportional rule
    enters as a COMMON factor over the eligible components, and a common factor cancels under
    normalisation. An informative per-component weight cancels nothing, which is what surfaces it, so
    the gate below runs both.
    """
    flat = (0.0, 0.0)
    part = _warm_start_partition(n_units=100, log_liks=flat, gdna_log_lik=-np.inf)
    weighted = _warm_start_run(
        _warm_start_estimator(2), part, [0, 1], weight=[3.0, 1.0], rna_prior=40.0
    )
    assert weighted[0] / weighted[1] == pytest.approx(3.0, rel=0.05), (
        "the RNA prior was discarded at a locus with no gDNA candidate"
    )
    # and the SHIPPED weights cancel there, which is what makes the defect invisible by default
    shipped = _warm_start_run(_warm_start_estimator(2), part, [0, 1], rna_prior=40.0)
    assert shipped[0] == pytest.approx(shipped[1], rel=1e-9)


def test_an_unknown_warm_start_is_REFUSED():
    """``Literal`` is a type-checker annotation, not a runtime constraint — so a typo must be refused
    explicitly, or it silently falls through to the shipped path."""
    with pytest.raises(ValueError, match="Unknown warm start"):
        EMConfig(warm_start="zero")


# ── The shadow-vs-gDNA contest: theta_n = 0 is stable only while L_n >= L_g ──────────────────────


def _shadow_locus(n_sense: int, n_anti: int, ss: float = 0.99):
    """One locus of unspliced fragments that ONLY gDNA and the gene's shadow span can explain — an
    intronic position, where no mature isoform reaches. Every fragment is gDNA in truth. gDNA is
    strandless (½ either way); the shadow carries the library's strandedness."""
    ti, ll, gll = [], [], []
    for n, p_rna in ((n_sense, ss), (n_anti, 1.0 - ss)):
        for _ in range(int(n)):
            ti += [0, 1]  # [0] the mature isoform, out of reach here; [1] the shadow
            ll += [-np.inf, float(np.log(p_rna))]
            gll.append(float(np.log(0.5)))
    n_units = int(n_sense + n_anti)
    return (
        np.arange(n_units + 1, dtype=np.int64) * 2,
        np.asarray(ti, dtype=np.int32),
        np.asarray(ll, dtype=np.float64),
        np.ones(n_units * 2, dtype=np.float64),
        np.zeros(n_units * 2, dtype=np.uint8),
        np.zeros(n_units, dtype=np.uint8),
        np.asarray(gll, dtype=np.float64),
        np.zeros(n_units, dtype=np.int32),
        np.zeros(n_units, dtype=np.uint8),
    )


def _shadow_share(ratio: float, *, n: int = 20_000, gdna_prior: float = 0.0) -> float:
    """The shadow's share of ``n`` fragments that are ALL gDNA in truth, at ``L_g / L_n = ratio``."""
    L_n = 10_000.0
    est = AbundanceEstimator(
        2,
        em_config=EMConfig(
            mode="vbem",
            iterations=400,
            convergence_delta=1e-12,
            assignment_mode="fractional",
            seed=0,
        ),
        # the mature isoform's own length is irrelevant here — it reaches no fragment in this locus
        geometry=_lengths([L_n * 10.0, L_n]),
    )
    est.run_batch_locus_em_partitioned(
        partition_tuples=[_shadow_locus(n // 2, n // 2)],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([gdna_prior], dtype=np.float64),
        gdna_eff_len=np.array([L_n * ratio], dtype=np.float64),
    )
    return float(est.em_counts[1].sum()) / n


@pytest.mark.parametrize("ratio", [0.5, 1.0])
def test_a_shadow_entity_holding_NOTHING_decays_to_ZERO_while_it_is_the_LONGER_component(ratio):
    """⛔ THE THRESHOLD IS `L_g / L_n = 1`, AND IT IS DERIVED (`EQUATIONS.md` §9b), NOT CHOSEN.

    Every fragment here is gDNA, so the only correct answer is that the shadow holds none of them.
    While the shadow's effective length is at least the gDNA component's, `theta_n = 0` is the fixed
    point the EM reaches and stays at: the shadow's density is multiplied by `L_g/L_n <= 1` each
    iteration and collapses.
    """
    assert _shadow_share(ratio) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize(
    "ratio,expected",
    [(2.0, 0.3519), (6.223, 0.5269), (20.0, 0.8204)],
)
def test_a_shadow_entity_holding_NOTHING_GROWS_once_it_is_the_SHORTER_component(ratio, expected):
    """⛔⛔ THE SIPHON'S MECHANISM (`ISSUES: nascent-siphons-gdna-under-capture`). Past the threshold
    `theta_n = 0` is UNSTABLE — the shadow's density is multiplied by `L_g/L_n > 1` every iteration and
    climbs off zero — and it settles where the strand channel alone stops it. The shares here are the
    closed-form fixed point of `EQUATIONS.md` §9b, so this gate pins the SOLVER against the derivation
    rather than against a recorded run.

    ⭐ `L_g` is the whole MultiLocus's opportunity and `L_n` one gene's span, so `L_g > L_n` is
    structural and grows with the component's gene count: on the ladder at `g50 ss.99 ON` the ratio's
    mass-weighted mean over the leaked fragments is 9.7.

    Perturbation: dividing the shadow by the locus's `L_g` instead of its own length reads 0.0 on every
    row. Dropping the strand term leaves the contest fully degenerate and it goes to the CORNER — the
    shorter component takes everything: 1.0 at every ratio above 1, exactly 0.5 at 1, 0.0 below. So the
    strand channel is not what opens this channel, only what bounds it.
    """
    assert _shadow_share(ratio) == pytest.approx(expected, abs=2e-3)


def test_the_gDNA_PSEUDOCOUNT_bounds_the_shadow_but_does_NOT_close_the_channel():
    """The calibration's own gDNA count anchors `theta_g`, so it damps the climb — which is why the
    ladder leaks ~10-20 % of the contested pool where the bare contest gives ~50 %. It is a bound, not
    a fix: the shadow still takes a large share, and at `g98`, where the anchor is strongest, the
    capture-ON rows still over-call nascent by 100x.
    """
    free = _shadow_share(6.223, gdna_prior=0.0)
    anchored = _shadow_share(6.223, gdna_prior=10_000.0)  # half the locus's fragments
    assert anchored < free, "the gDNA pseudocount did not damp the shadow at all"
    assert anchored > 0.3, "the pseudocount closed a channel it can only bound"
