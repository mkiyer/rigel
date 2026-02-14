"""Tests for hulkrna.counter — ReadCounter with unified EM assignment."""

import numpy as np
import pytest

from hulkrna.types import Strand, MergeCriteria
from hulkrna.categories import (
    AssignmentSource,
    CountCategory,
    CountStrand,
    CountType,
    NUM_COUNT_TYPES,
)
from hulkrna.resolution import ResolvedFragment
from hulkrna.counter import ReadCounter, EMData
from hulkrna.strand_model import StrandModel, StrandModels


# =====================================================================
# Helpers — lightweight mock index + fixtures
# =====================================================================


class MockIndex:
    """Minimal mock of HulkIndex with the arrays ReadCounter needs."""

    def __init__(self, num_transcripts, num_genes, t_to_g, t_to_strand, g_to_strand):
        self.num_transcripts = num_transcripts
        self.num_genes = num_genes
        self.t_to_g_arr = np.array(t_to_g, dtype=np.int64)
        self.t_to_strand_arr = np.array(t_to_strand, dtype=np.int64)
        self.g_to_strand_arr = np.array(g_to_strand, dtype=np.int64)


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


def _make_strand_model_fr():
    """FR-strand model: p_r1_sense ≈ 0.98."""
    sm = StrandModel()
    for _ in range(100):
        sm.observe(Strand.POS, Strand.POS)
    for _ in range(2):
        sm.observe(Strand.POS, Strand.NEG)
    return sm


def _make_strand_models_fr():
    """StrandModels container with FR-strand model in all slots."""
    sm = _make_strand_model_fr()
    return StrandModels(
        exonic_spliced=sm,
        exonic=sm,
        intronic=sm,
        intergenic=StrandModel(),
    )


def _make_resolved(**kwargs):
    defaults = dict(
        t_inds=frozenset({0}),
        n_genes=1,
        count_cat=CountCategory.UNSPLICED,
        exon_strand=Strand.POS,
        sj_strand=Strand.NONE,
        insert_size=250,
        merge_criteria=MergeCriteria.INTERSECTION,
        num_hits=1,
    )
    defaults.update(kwargs)
    return ResolvedFragment(**defaults)


def _make_insert_models():
    from hulkrna.insert_model import InsertSizeModels
    im = InsertSizeModels()
    im.observe(250, CountCategory.UNSPLICED)
    return im


def _make_em_data(
    t_indices_per_unit,
    log_liks_per_unit=None,
    count_types_per_unit=None,
    gdna_base_index=None,
):
    """Build EMData from a list of per-unit candidate lists.

    Parameters
    ----------
    t_indices_per_unit : list[list[int]]
        Candidate transcript indices for each unit.
    log_liks_per_unit : list[list[float]] or None
        Log-likelihoods per candidate (default: all zeros = equal).
    count_types_per_unit : list[list[int]] or None
        Count type column index per candidate (default: all 4 = UNSPLICED_SENSE).
    gdna_base_index : int or None
        Base gDNA shadow index. Default: ``max(all_t_indices) + 1`` or 0
        if there are no candidates.
    """
    offsets = [0]
    flat_t = []
    flat_lk = []
    flat_ct = []

    for u, t_list in enumerate(t_indices_per_unit):
        for j, t_idx in enumerate(t_list):
            flat_t.append(t_idx)
            if log_liks_per_unit is not None:
                flat_lk.append(log_liks_per_unit[u][j])
            else:
                flat_lk.append(0.0)
            if count_types_per_unit is not None:
                flat_ct.append(count_types_per_unit[u][j])
            else:
                flat_ct.append(int(CountType.UNSPLICED_SENSE))
        offsets.append(len(flat_t))

    n_units = len(t_indices_per_unit)
    n_candidates = len(flat_t)

    # Derive gdna_base_index sentinel (one past max transcript index)
    if gdna_base_index is None:
        gdna_base_index = (max(flat_t) + 1) if flat_t else 0

    # Build locus tracking arrays: best non-gDNA transcript per unit.
    # For each unit pick the first non-gDNA candidate (or -1).
    locus_t = np.full(n_units, -1, dtype=np.int32)
    locus_ct = np.zeros(n_units, dtype=np.uint8)
    for u, t_list in enumerate(t_indices_per_unit):
        ct_list = count_types_per_unit[u] if count_types_per_unit else None
        for j, t_idx in enumerate(t_list):
            if t_idx < gdna_base_index:
                locus_t[u] = t_idx
                locus_ct[u] = ct_list[j] if ct_list else int(CountType.UNSPLICED_SENSE)
                break

    return EMData(
        offsets=np.array(offsets, dtype=np.int64),
        t_indices=np.array(flat_t, dtype=np.int32),
        log_liks=np.array(flat_lk, dtype=np.float64),
        count_types=np.array(flat_ct, dtype=np.uint8),
        locus_t_indices=locus_t,
        locus_count_types=locus_ct,
        n_units=n_units,
        n_candidates=n_candidates,
        gdna_base_index=gdna_base_index,
    )


# =====================================================================
# classify_strand
# =====================================================================


class TestClassifyStrand:
    def test_high_likelihood_returns_sense(self):
        sm = _make_strand_model_fr()
        cs = ReadCounter.classify_strand(Strand.POS, int(Strand.POS), sm)
        assert cs == CountStrand.SENSE

    def test_low_likelihood_returns_antisense(self):
        sm = _make_strand_model_fr()
        cs = ReadCounter.classify_strand(Strand.POS, int(Strand.NEG), sm)
        assert cs == CountStrand.ANTISENSE

    def test_ambiguous_for_unstranded(self):
        sm = StrandModel()  # no observations → p_r1_sense = 0.5
        cs = ReadCounter.classify_strand(Strand.POS, int(Strand.POS), sm)
        assert cs == CountStrand.AMBIGUOUS


# =====================================================================
# assign_unique
# =====================================================================


class TestAssignUnique:
    def test_single_transcript_single_gene(self):
        index = _make_index()
        sms = _make_strand_models_fr()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)

        resolved = _make_resolved(
            t_inds=frozenset({0}),
            n_genes=1,
            exon_strand=Strand.POS,
            count_cat=CountCategory.UNSPLICED,
        )
        rc.assign_unique(resolved, index, sms)

        col = int(CountType.UNSPLICED_SENSE)
        assert rc.unique_counts[0, col] == 1.0
        assert rc.t_counts[0, col] == 1.0

    def test_empty_t_inds_does_nothing(self):
        index = _make_index()
        sms = _make_strand_models_fr()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)

        resolved = _make_resolved(t_inds=frozenset(), n_genes=0)
        rc.assign_unique(resolved, index, sms)

        assert rc.unique_counts.sum() == 0.0

    def test_always_assigns_one_count(self):
        index = _make_index()
        sms = _make_strand_models_fr()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)

        resolved = _make_resolved(
            t_inds=frozenset({0}),
            n_genes=1,
            exon_strand=Strand.POS,
            num_hits=4,
        )
        rc.assign_unique(resolved, index, sms)

        col = int(CountType.UNSPLICED_SENSE)
        assert rc.unique_counts[0, col] == 1.0

    def test_writes_to_unique_not_em(self):
        index = _make_index()
        sms = _make_strand_models_fr()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)

        resolved = _make_resolved(t_inds=frozenset({0}), n_genes=1)
        rc.assign_unique(resolved, index, sms)

        assert rc.unique_counts.sum() == 1.0
        assert rc.em_counts.sum() == 0.0


# =====================================================================
# EMData construction
# =====================================================================


class TestEMData:
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


# =====================================================================
# run_em — EM convergence tests
# =====================================================================


class TestRunEM:
    def test_empty_em_data_sets_theta(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data([])
        rc.run_em(em, em_iterations=5)
        assert rc._converged_theta is not None
        assert rc._converged_theta.shape == (5,)  # N_t + N_g (3 + 2)
        assert rc._converged_theta.sum() == pytest.approx(1.0)

    def test_uniform_priors_stay_uniform(self):
        """With no unique counts and equal likelihoods, theta stays uniform."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        # 1000 fragments, each with candidates [0, 1], equal log-liks
        em = _make_em_data([[0, 1]] * 1000)
        rc.run_em(em, em_iterations=10)

        theta = rc._converged_theta
        # t0 and t1 should be roughly equal, t2 only has alpha
        assert theta[0] == pytest.approx(theta[1], rel=0.01)
        assert theta[2] < theta[0]  # t2 has no evidence

    def test_unique_counts_bias_em(self):
        """Unique counts on t0 should bias EM toward t0."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        # Give t0 many unique counts
        rc.unique_counts[0, int(CountType.UNSPLICED_SENSE)] = 100.0

        # 500 ambiguous fragments between t0 and t1
        em = _make_em_data([[0, 1]] * 500)
        rc.run_em(em, em_iterations=10)

        theta = rc._converged_theta
        assert theta[0] > theta[1]

    def test_zero_iterations_uses_unique_priors(self):
        """With em_iterations=0, theta is just unique counts + alpha."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, 0] = 10.0
        rc.unique_counts[1, 0] = 5.0

        em = _make_em_data([[0, 1]] * 100)
        rc.run_em(em, em_iterations=0)

        theta = rc._converged_theta
        # theta has N_t + N_g elements:
        # t0=10+1=11, t1=5+1=6, t2=0+1=1, shadow_g0=0+1=1, shadow_g1=0+1=1
        # sum=20
        assert theta[0] == pytest.approx(11.0 / 20.0)
        assert theta[1] == pytest.approx(6.0 / 20.0)
        assert theta[2] == pytest.approx(1.0 / 20.0)
        assert theta[3] == pytest.approx(1.0 / 20.0)  # shadow g0
        assert theta[4] == pytest.approx(1.0 / 20.0)  # shadow g1

    def test_likelihood_influences_em(self):
        """Candidates with higher likelihoods should attract more mass."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        # 1000 fragments: candidate t0 has high log-lik, t1 has low
        em = _make_em_data(
            [[0, 1]] * 1000,
            log_liks_per_unit=[[0.0, -10.0]] * 1000,
        )
        rc.run_em(em, em_iterations=10)

        theta = rc._converged_theta
        assert theta[0] > theta[1]

    def test_convergence_detected(self):
        """EM should converge early for simple problems."""
        rc = ReadCounter(num_transcripts=2, num_genes=1, seed=42)
        rc.unique_counts[0, 0] = 100.0
        rc.unique_counts[1, 0] = 100.0

        em = _make_em_data([[0, 1]] * 100)
        rc.run_em(em, em_iterations=100)  # allow many iterations
        # Should converge way before 100 iterations
        theta = rc._converged_theta
        assert theta is not None


# =====================================================================
# assign_ambiguous — stochastic assignment using EM posteriors
# =====================================================================


class TestAssignAmbiguous:
    def test_single_candidate_gets_full_count(self):
        """Unit with one candidate → count goes to that transcript."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data([[0]])
        rc.run_em(em, em_iterations=1)
        rc.assign_ambiguous(em)

        assert rc.em_counts[0].sum() == 1.0
        assert rc.em_counts.sum() == 1.0

    def test_two_candidates_picks_one(self):
        """Each unit contributes exactly 1 count."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data([[0, 1]])
        rc.run_em(em, em_iterations=1)
        rc.assign_ambiguous(em)

        assert rc.em_counts.sum() == 1.0

    def test_n_units_equals_n_counts(self):
        """N ambiguous units → N total em_counts."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data([[0, 1]] * 100)
        rc.run_em(em, em_iterations=5)
        rc.assign_ambiguous(em)

        assert rc.em_counts.sum() == 100.0

    def test_total_counts_equals_unique_plus_em(self):
        """t_counts = unique_counts + em_counts."""
        index = _make_index()
        sms = _make_strand_models_fr()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)

        # 10 unique
        for _ in range(10):
            rc.assign_unique(
                _make_resolved(t_inds=frozenset({0})),
                index, sms,
            )

        # 5 ambiguous
        em = _make_em_data([[0, 1]] * 5)
        rc.run_em(em, em_iterations=5)
        rc.assign_ambiguous(em)

        assert rc.unique_counts.sum() == 10.0
        assert rc.em_counts.sum() == 5.0
        assert rc.t_counts.sum() == 15.0
        np.testing.assert_array_equal(
            rc.t_counts, rc.unique_counts + rc.em_counts
        )

    def test_distribution_follows_priors(self):
        """Over many fragments, EM + assignment follows unique priors."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, int(CountType.UNSPLICED_SENSE)] = 90.0
        rc.unique_counts[1, int(CountType.UNSPLICED_SENSE)] = 10.0

        em = _make_em_data([[0, 1]] * 10000)
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        # t0 should get roughly 90% of em_counts
        t0_em = rc.em_counts[0].sum()
        t1_em = rc.em_counts[1].sum()
        assert t0_em > t1_em

    def test_writes_to_em_not_unique(self):
        """assign_ambiguous only writes to em_counts."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data([[0, 1]] * 50)
        rc.run_em(em, em_iterations=5)

        unique_before = rc.unique_counts.copy()
        rc.assign_ambiguous(em)

        np.testing.assert_array_equal(rc.unique_counts, unique_before)
        assert rc.em_counts.sum() == 50.0

    def test_empty_em_data_does_nothing(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data([])
        rc.run_em(em, em_iterations=5)
        rc.assign_ambiguous(em)

        assert rc.em_counts.sum() == 0.0

    def test_count_type_respected(self):
        """Count goes to the correct CountType column."""
        ct_col = int(CountType.SPLICED_ANNOT_ANTISENSE)
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data(
            [[0]],
            count_types_per_unit=[[ct_col]],
        )
        rc.run_em(em, em_iterations=1)
        rc.assign_ambiguous(em)

        assert rc.em_counts[0, ct_col] == 1.0
        assert rc.em_counts.sum() == 1.0

    def test_requires_run_em_first(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data([[0, 1]])
        with pytest.raises(RuntimeError, match="run_em"):
            rc.assign_ambiguous(em)


# =====================================================================
# Simultaneous resolution — no phase ordering bias
# =====================================================================


class TestSimultaneousResolution:
    """Verify that isoform-ambig and gene-ambig fragments
    are resolved simultaneously, not sequentially."""

    def test_same_result_regardless_of_order(self):
        """Shuffling ambiguous units produces statistically similar results.

        This tests the key property: no order dependence.
        """
        rc1 = ReadCounter(num_transcripts=4, num_genes=2, seed=42)
        rc1.unique_counts[0, 0] = 50.0
        rc1.unique_counts[2, 0] = 30.0

        # Mix of isoform-ambig (within gene 0: t0, t1) and
        # gene-ambig (across genes: t0, t2)
        units = [[0, 1]] * 200 + [[0, 2]] * 200
        em1 = _make_em_data(units)
        rc1.run_em(em1, em_iterations=10)
        rc1.assign_ambiguous(em1)

        # Same but reversed order of units
        rc2 = ReadCounter(num_transcripts=4, num_genes=2, seed=42)
        rc2.unique_counts[0, 0] = 50.0
        rc2.unique_counts[2, 0] = 30.0

        units_rev = [[0, 2]] * 200 + [[0, 1]] * 200
        em2 = _make_em_data(units_rev)
        rc2.run_em(em2, em_iterations=10)
        rc2.assign_ambiguous(em2)

        # Converged theta should be identical (EM is order-independent)
        np.testing.assert_allclose(
            rc1._converged_theta, rc2._converged_theta, atol=1e-10
        )


# =====================================================================
# Multimapper-like units (multi-alignment molecules)
# =====================================================================


class TestMultimapperEM:
    def test_multimapper_molecule_one_count(self):
        """Multimapper molecule (many candidates) → exactly 1 count."""
        rc = ReadCounter(num_transcripts=4, num_genes=2, seed=42)
        # Simulate molecule with 3 alignments, each to different transcripts
        em = _make_em_data([[0, 1, 2, 3]])
        rc.run_em(em, em_iterations=5)
        rc.assign_ambiguous(em)

        assert rc.em_counts.sum() == 1.0

    def test_multimapper_and_ambig_together(self):
        """Mixed ambiguous + multimapper units in same EM."""
        rc = ReadCounter(num_transcripts=4, num_genes=2, seed=42)
        rc.unique_counts[0, 0] = 100.0

        # 50 isoform-ambig, 50 gene-ambig, 25 multimapper
        units = [[0, 1]] * 50 + [[0, 2]] * 50 + [[0, 1, 2, 3]] * 25
        em = _make_em_data(units)
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        assert rc.em_counts.sum() == 125.0


# =====================================================================
# Wide DataFrame output
# =====================================================================


class TestCounterOutput:
    def test_get_t_counts_df(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        df = rc.get_t_counts_df()
        assert df.shape == (3, NUM_COUNT_TYPES)
        assert list(df.columns) == list(CountType.columns())

    def test_get_g_counts_df(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        t_to_g_arr = np.array([0, 0, 1], dtype=np.int64)
        df = rc.get_g_counts_df(t_to_g_arr)
        assert df.shape == (2, NUM_COUNT_TYPES)
        assert list(df.columns) == list(CountType.columns())

    def test_wide_counts_include_both_sources(self):
        """Wide output sums unique + em counts."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, 0] = 5.0
        rc.em_counts[0, 0] = 3.0

        df = rc.get_t_counts_df()
        assert df.iloc[0, 0] == 8.0


# =====================================================================
# Sparse DataFrame output
# =====================================================================


class TestSparseOutput:
    def test_sparse_empty(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        df = rc.get_sparse_t_counts_df()
        assert len(df) == 0
        assert list(df.columns) == [
            "transcript_idx", "category", "strand", "source", "count"
        ]

    def test_sparse_t_counts_unique_only(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, int(CountType.UNSPLICED_SENSE)] = 5.0
        rc.unique_counts[2, int(CountType.SPLICED_ANNOT_SENSE)] = 3.0

        df = rc.get_sparse_t_counts_df()
        assert len(df) == 2
        assert set(df["source"]) == {"unique"}
        assert set(df["transcript_idx"]) == {0, 2}
        assert df["count"].sum() == 8

    def test_sparse_t_counts_both_sources(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, int(CountType.UNSPLICED_SENSE)] = 5.0
        rc.em_counts[0, int(CountType.UNSPLICED_SENSE)] = 3.0

        df = rc.get_sparse_t_counts_df()
        assert len(df) == 2
        assert set(df["source"]) == {"unique", "em"}
        assert df["count"].sum() == 8

    def test_sparse_g_counts(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, int(CountType.UNSPLICED_SENSE)] = 5.0
        rc.unique_counts[1, int(CountType.UNSPLICED_SENSE)] = 3.0
        t_to_g_arr = np.array([0, 0, 1], dtype=np.int64)

        df = rc.get_sparse_g_counts_df(t_to_g_arr)
        # g0 gets 5+3=8 from unique
        assert len(df) == 1
        row = df.iloc[0]
        assert row["gene_idx"] == 0
        assert row["count"] == 8
        assert row["source"] == "unique"

    def test_sparse_columns_are_categorical(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, 0] = 1.0

        df = rc.get_sparse_t_counts_df()
        assert hasattr(df["category"].dtype, "categories")
        assert hasattr(df["strand"].dtype, "categories")
