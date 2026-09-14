"""The locus EM — `AbundanceEstimator`, the native batch wrapper it drives, and the warm start that
seeds it. The first block runs the estimator over `ScoredFragments` and a mock index: locus assignment,
posterior means, simultaneous resolution, multimappers, the counts and detail tables, partitioned
effective lengths, the gDNA component and discrete assignment. The second calls
`run_batch_locus_em_partitioned` directly on the production path, gating structural gDNA eligibility,
the gDNA effective length and log-odds bias, the aggregate RNA prior, and the rule that withholds that
prior from synthetic nascent entities. The third gates `EMConfig.warm_start`: what the initial theta is
derived from, and that an unknown name is refused.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.types import Strand
from rigel.splice import (
    SpliceType,
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

    # Build locus tracking arrays
    locus_t = np.full(n_units, -1, dtype=np.int32)
    locus_cc = np.zeros(n_units, dtype=np.uint8)
    for u, t_list in enumerate(t_indices_per_unit):
        cc_list = count_cols_per_unit[u] if count_cols_per_unit else None
        for j, t_idx in enumerate(t_list):
            locus_t[u] = t_idx
            locus_cc[u] = cc_list[j] if cc_list else _UNSPLICED_SENSE
            break

    return ScoredFragments(
        offsets=np.array(offsets, dtype=np.int64),
        t_indices=np.array(flat_t, dtype=np.int32),
        log_liks=np.array(flat_lk, dtype=np.float64),
        count_cols=np.array(flat_cc, dtype=np.uint8),
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
        """With large gdna_init and equal likelihoods, gDNA takes share."""
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
        """em_counts + nrna_em + gdna == n_units."""
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


def _estimator(
    n_t: int, *, mode: str = "vbem", gdna_em_llr_bias: float = 0.0
) -> AbundanceEstimator:
    est = AbundanceEstimator(
        num_transcripts=n_t,
        em_config=EMConfig(
            mode=mode,
            iterations=200,
            convergence_delta=1e-9,
            assignment_mode="fractional",
            seed=0,
            gdna_em_llr_bias=gdna_em_llr_bias,
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


def test_gdna_em_llr_bias_favors_gdna_assignment():
    # An ambiguous fragment (transcript and gDNA equally likely, equal eff-lens):
    # neutral (bias=0) splits ~evenly; a positive gDNA LLR bias shifts mass toward
    # gDNA (the FP-aversion knob — fewer gDNA->RNA leaks, more RNA->gDNA siphons);
    # a negative bias shifts toward RNA. Monotone in the bias.
    n_units = 100
    partition = _partition(n_units=n_units, log_liks=(0.0,), gdna_log_lik=0.0)
    kw = dict(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        gdna_eff_len=np.array([1.0], dtype=np.float64),
        enable_gdna=np.array([1], dtype=np.uint8),
    )
    _t, _r, gdna_neg = _estimator(1, gdna_em_llr_bias=-2.0).run_batch_locus_em_partitioned(**kw)
    _t, _r, gdna_neutral = _estimator(1, gdna_em_llr_bias=0.0).run_batch_locus_em_partitioned(**kw)
    _t, _r, gdna_pos = _estimator(1, gdna_em_llr_bias=2.0).run_batch_locus_em_partitioned(**kw)

    assert gdna_neg[0] < gdna_neutral[0] < gdna_pos[0]
    # neutral splits an exactly-tied 1-transcript-vs-gDNA fragment ~50/50
    assert gdna_neutral[0] == pytest.approx(n_units / 2, rel=0.05)
    # a +2-nat (~7.4:1 odds) bias pushes the gDNA share well above half
    assert gdna_pos[0] > 0.8 * n_units


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


# ---------------------------------------------------------------------------
# Structural gDNA eligibility: ``gdna_prior_count`` no longer gates the gDNA
# component, and the compatibility ``enable_gdna`` array is not a modeling gate
# — native derives candidate availability from the partition itself.
# ---------------------------------------------------------------------------


def test_compat_enable_false_does_not_disable_structural_candidate():
    """A compatibility ``enable_gdna=False`` input is ignored by native v3."""
    n_units = 50
    est = _estimator(2, mode="map")
    total_gdna, _rna, _g = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(-1.0, -1.0), gdna_log_lik=-0.5)],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([100.0], dtype=np.float64),
        index=None,
        enable_gdna=np.array([0], dtype=np.uint8),
    )

    assert total_gdna > 0.5 * n_units, (
        "native should derive gDNA availability from finite unspliced candidates, "
        f"not the compatibility enable_gdna array; got total_gdna={total_gdna}"
    )


def test_default_enable_gdna_inferred_from_partition():
    """When ``enable_gdna`` is None, the wrapper computes it from the partition
    (any unspliced unit with a finite gDNA log-lik ⇒ enabled).
    """
    n_units = 20
    locus_t_lists = [np.array([0, 1], dtype=np.int32)]

    # All-spliced partition ⇒ no per-unit gDNA candidate ⇒ enable=0.
    spliced_part = _partition(
        n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-0.5, is_spliced=True
    )
    # Unspliced partition with finite gDNA log-liks ⇒ enable=1.
    unspliced_part = _partition(
        n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-0.5, is_spliced=False
    )
    # Unspliced but non-finite gDNA log-lik ⇒ enable=0.
    nogdna_part = _partition(
        n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-np.inf, is_spliced=False
    )

    # All three with gdna_prior_count=0 and enable_gdna omitted (inferred).
    # Only the unspliced+finite case produces gDNA assignments.
    g_spl, _, _ = _estimator(2, mode="map").run_batch_locus_em_partitioned(
        [spliced_part],
        locus_t_lists,
        np.zeros(1),
        index=None,
    )
    g_uns, _, _ = _estimator(2, mode="map").run_batch_locus_em_partitioned(
        [unspliced_part],
        locus_t_lists,
        np.zeros(1),
        index=None,
    )
    g_no, _, _ = _estimator(2, mode="map").run_batch_locus_em_partitioned(
        [nogdna_part],
        locus_t_lists,
        np.zeros(1),
        index=None,
    )

    assert g_spl == 0.0, "spliced partition has no gDNA candidates"
    assert g_no == 0.0, "non-finite gDNA log-liks ⇒ no gDNA candidates"
    assert g_uns > 0.0, (
        "unspliced+finite gDNA log-lik must enable component even when gdna_prior_count == 0"
    )


def test_positive_gdna_prior_produces_finite_outputs():
    """A positive gDNA prior count produces finite, conserved outputs."""
    n_units = 30
    est = _estimator(2, mode="map")
    total_gdna, locus_rna, locus_gdna = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-1.0)],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
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


# ── SYNTHETIC NASCENT ENTITIES GET NO RNA PRIOR ──────────────────────────────────────────────────
#
# A synthetic nascent entity is a shadow span the INDEX manufactured; no annotation asserts it exists.
# The null hypothesis is therefore that it is ABSENT, and it earns mass only from fragments the data
# cannot explain any other way. These gates pin that behaviour and, just as importantly, pin the two
# things it must NOT do: kill an entity the data supports, and disturb a locus that has none.


class _StubIndex:
    """The minimum `run_batch_locus_em_partitioned` reads: `t_df["is_synthetic"]`."""

    def __init__(self, flags):
        import pandas as pd

        self.t_df = pd.DataFrame({"is_synthetic": np.asarray(flags, dtype=bool)})


def _run(est, partition, t_idx, *, index=None, rna_prior=0.0, gdna_prior=0.0, enable_gdna=1):
    return est.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.asarray(t_idx, dtype=np.int32)],
        gdna_prior_count=np.array([gdna_prior], dtype=np.float64),
        rna_prior_count=np.array([rna_prior], dtype=np.float64),
        index=index,
        enable_gdna=np.array([enable_gdna], dtype=np.uint8),
    )


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_locus_with_no_synthetic_component_is_BIT_IDENTICAL(mode):
    """The control. The rule must be invisible where it does not apply — and BIT-identical, not
    close: the C++ skips the annotated/synthetic split entirely when the mask is empty, and that is
    what makes every other number here attributable (TRAPS: byte-identity-gate)."""
    outs = []
    for index in (None, _StubIndex([False, False])):
        est = _estimator(2, mode=mode)
        _run(
            est,
            _partition(n_units=100, log_liks=(0.0, -0.5), gdna_log_lik=-1.0),
            [0, 1],
            index=index,
            rna_prior=50.0,
            gdna_prior=20.0,
        )
        outs.append(est.em_counts.sum(axis=1).copy())
    np.testing.assert_array_equal(outs[0], outs[1])


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_SHARED_ONLY_synthetic_entity_LOSES_mass_to_the_annotated_one(mode):
    """The zombie. Every fragment it holds is equally well explained by the annotated transcript,
    so it has no evidence of its own and the prior no longer props it up."""
    part = _partition(n_units=200, log_liks=(0.0, 0.0), gdna_log_lik=-2.0)
    counts = {}
    for tag, index in (("shipped", None), ("gated", _StubIndex([False, True]))):
        est = _estimator(2, mode=mode)
        _run(est, part, [0, 1], index=index, rna_prior=100.0, gdna_prior=10.0)
        counts[tag] = est.em_counts.sum(axis=1).copy()
    assert counts["gated"][1] < counts["shipped"][1], (
        f"the synthetic component did not lose mass: {counts}"
    )
    assert counts["gated"][0] > counts["shipped"][0], (
        f"the annotated component did not gain what the synthetic lost: {counts}"
    )


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_synthetic_entity_the_DATA_SUPPORTS_still_survives(mode):
    """The gate that stops this becoming "kill all nascent RNA". A rule that zeroes everything
    passes the zombie test above; only this one separates the two. The synthetic component is the
    strictly better explanation here, and it must keep the mass."""
    est = _estimator(2, mode=mode)
    _run(
        est,
        _partition(n_units=200, log_liks=(-8.0, 0.0), gdna_log_lik=-8.0),
        [0, 1],
        index=_StubIndex([False, True]),
        rna_prior=100.0,
        gdna_prior=10.0,
    )
    counts = est.em_counts.sum(axis=1)
    assert counts[1] > 0.9 * counts.sum(), (
        f"a synthetic entity with decisive likelihood support was suppressed anyway: {counts}"
    )


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_an_ALL_SYNTHETIC_locus_takes_NO_rna_prior(mode):
    """The degenerate locus of the derivation: no annotated component means no eligible recipient,
    so the RNA prior is 0 and the pool must outcompete gDNA unaided. It must not silently fall back
    to handing the prior to the synthetic components after all."""
    part = _partition(n_units=100, log_liks=(0.0, 0.0), gdna_log_lik=0.0)
    res = {}
    for tag, rna_prior in (("no_prior", 0.0), ("big_prior", 500.0)):
        est = _estimator(2, mode=mode)
        _total, _rna, gdna = _run(
            est, part, [0, 1], index=_StubIndex([True, True]), rna_prior=rna_prior, gdna_prior=10.0
        )
        res[tag] = (est.em_counts.sum(axis=1).copy(), float(gdna[0]))
    np.testing.assert_allclose(res["no_prior"][0], res["big_prior"][0], rtol=1e-9, atol=1e-9)
    assert res["no_prior"][1] == pytest.approx(res["big_prior"][1], rel=1e-9), (
        "the RNA prior reached an all-synthetic locus and moved the gDNA split"
    )


def test_is_synthetic_is_read_NOT_is_nrna():
    """A single-exon annotated transcript carries `is_nrna = True` and is simultaneously the
    nascent and the mature form of a REAL gene. It must keep its prior. The estimator must key on
    `is_synthetic` alone, so an index carrying only `is_nrna` changes nothing."""
    import pandas as pd

    class _NrnaOnlyIndex:
        def __init__(self):
            self.t_df = pd.DataFrame({"is_nrna": [False, True]})

    est = AbundanceEstimator(num_transcripts=2, em_config=EMConfig(mode="map"))
    np.testing.assert_array_equal(est._t_is_synthetic(_NrnaOnlyIndex(), 2), np.zeros(0, np.uint8))
    np.testing.assert_array_equal(
        est._t_is_synthetic(_StubIndex([False, False]), 2), np.zeros(0, np.uint8)
    )
    np.testing.assert_array_equal(
        est._t_is_synthetic(_StubIndex([False, True]), 2), np.array([0, 1], np.uint8)
    )


# The one invariant this design rests on has no test, and cannot have one from Python.
#
# `apply_grouped_prior_update` guarantees, FOR A GIVEN raw_counts VECTOR, that the RNA components sum
# to `rna_count + rna_prior` — so withholding the prior from synthetic components redistributes mass
# strictly WITHIN the RNA pool. It is asserted only in a C++ comment.
#
# It is a PER-M-STEP algebraic identity, NOT an end-to-end one, and conflating the two has produced a
# wrong test twice. End to end the gDNA total legitimately DEPENDS on `rna_prior`: setting the
# gDNA:RNA split is the prior's whole purpose (`gdna_total = gdna_count + gdna_prior`,
# `rna_total = rna_count + rna_prior`), and a larger RNA prior shifts theta, hence the E-step, hence
# the next iteration's `gdna_count`. A test asserting "the gDNA total is invariant to rna_prior" is
# therefore asserting something false by design, and it fails on the no-synthetic configuration too.
#
# The function is `static` in em_solver.cpp, so nothing can call it directly. Exposing it behind a
# test-only binding — the executable-specification pattern this repo already uses for the accumulator
# — is the way to gate the identity, and it is a prerequisite for per-transcript prior work, since a
# per-transcript vector makes the identity harder rather than easier to hold.


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
    enable_gdna=1,
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
        index=None,
        enable_gdna=np.array([enable_gdna], dtype=np.uint8),
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
        _warm_start_estimator(2), part, [0, 1], weight=[3.0, 1.0], rna_prior=40.0, enable_gdna=0
    )
    assert weighted[0] / weighted[1] == pytest.approx(3.0, rel=0.05), (
        "the RNA prior was discarded at a locus with no gDNA candidate"
    )
    # and the SHIPPED weights cancel there, which is what makes the defect invisible by default
    shipped = _warm_start_run(_warm_start_estimator(2), part, [0, 1], rna_prior=40.0, enable_gdna=0)
    assert shipped[0] == pytest.approx(shipped[1], rel=1e-9)


def test_an_unknown_warm_start_is_REFUSED():
    """``Literal`` is a type-checker annotation, not a runtime constraint — so a typo must be refused
    explicitly, or it silently falls through to the shipped path."""
    with pytest.raises(ValueError, match="Unknown warm start"):
        EMConfig(warm_start="zero")
