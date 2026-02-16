"""Tests for hulkrna.counter — ReadCounter with unified EM assignment."""

import numpy as np
import pandas as pd
import pytest

from hulkrna.types import Strand, MergeCriteria
from hulkrna.categories import (
    CountCategory,
    CountCol,
    NUM_COUNT_COLS,
)
from hulkrna.resolution import ResolvedFragment
from hulkrna.counter import ReadCounter, EMData
from hulkrna.strand_model import StrandModel, StrandModels


# =====================================================================
# Helpers — lightweight mock index + fixtures
# =====================================================================


class MockIndex:
    """Minimal mock of HulkIndex with the arrays ReadCounter needs."""

    def __init__(self, num_transcripts, num_genes, t_to_g, t_to_strand, g_to_strand,
                 t_ids=None, g_ids=None, g_names=None, t_gnames=None):
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

        self.t_df = pd.DataFrame({
            "t_id": t_ids,
            "g_id": [g_ids[g] for g in t_to_g],
            "g_name": t_gnames,
        })
        self.g_df = pd.DataFrame({
            "g_id": g_ids,
            "g_name": g_names,
        })


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


# Default column for UNSPLICED_SENSE
_UNSPLICED_SENSE = int(CountCol.UNSPLICED_SENSE)


def _make_em_data(
    t_indices_per_unit,
    log_liks_per_unit=None,
    count_cols_per_unit=None,
    gdna_base_index=None,
):
    """Build EMData from a list of per-unit candidate lists.

    Parameters
    ----------
    t_indices_per_unit : list[list[int]]
        Candidate transcript indices for each unit.
    log_liks_per_unit : list[list[float]] or None
        Log-likelihoods per candidate (default: all zeros = equal).
    count_cols_per_unit : list[list[int]] or None
        Count column index per candidate (default: all UNSPLICED_SENSE = 2).
    gdna_base_index : int or None
        Base gDNA shadow index. Default: ``max(all_t_indices) + 1`` or 0
        if there are no candidates.
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

    # Derive gdna_base_index sentinel (one past max transcript index)
    if gdna_base_index is None:
        gdna_base_index = (max(flat_t) + 1) if flat_t else 0

    # Build locus tracking arrays: best non-gDNA transcript per unit.
    locus_t = np.full(n_units, -1, dtype=np.int32)
    locus_cc = np.zeros(n_units, dtype=np.uint8)
    for u, t_list in enumerate(t_indices_per_unit):
        cc_list = count_cols_per_unit[u] if count_cols_per_unit else None
        for j, t_idx in enumerate(t_list):
            if t_idx < gdna_base_index:
                locus_t[u] = t_idx
                locus_cc[u] = cc_list[j] if cc_list else _UNSPLICED_SENSE
                break

    return EMData(
        offsets=np.array(offsets, dtype=np.int64),
        t_indices=np.array(flat_t, dtype=np.int32),
        log_liks=np.array(flat_lk, dtype=np.float64),
        count_cols=np.array(flat_cc, dtype=np.uint8),
        locus_t_indices=locus_t,
        locus_count_cols=locus_cc,
        n_units=n_units,
        n_candidates=n_candidates,
        gdna_base_index=gdna_base_index,
    )


# =====================================================================
# is_antisense
# =====================================================================


class TestIsAntisense:
    def test_high_likelihood_returns_false(self):
        """Good strand match → sense → not antisense."""
        sm = _make_strand_model_fr()
        result = ReadCounter.is_antisense(Strand.POS, int(Strand.POS), sm)
        assert result is False

    def test_low_likelihood_returns_true(self):
        """Poor strand match → antisense."""
        sm = _make_strand_model_fr()
        result = ReadCounter.is_antisense(Strand.POS, int(Strand.NEG), sm)
        assert result is True

    def test_unstranded_returns_false(self):
        """Unstranded model (p=0.5) is treated as non-antisense."""
        sm = StrandModel()  # no observations → p = 0.5
        result = ReadCounter.is_antisense(Strand.POS, int(Strand.POS), sm)
        assert result is False


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

        col = _UNSPLICED_SENSE
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

        col = _UNSPLICED_SENSE
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
        em = _make_em_data([[0, 1]] * 1000)
        rc.run_em(em, em_iterations=10)

        theta = rc._converged_theta
        assert theta[0] == pytest.approx(theta[1], rel=0.01)
        assert theta[2] < theta[0]

    def test_unique_counts_bias_em(self):
        """Unique counts on t0 should bias EM toward t0."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, _UNSPLICED_SENSE] = 100.0

        em = _make_em_data([[0, 1]] * 500)
        rc.run_em(em, em_iterations=10)

        theta = rc._converged_theta
        assert theta[0] > theta[1]

    def test_zero_iterations_uses_unique_priors(self):
        """With em_iterations=0, theta comes from unique counts + prior only."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, 0] = 10.0
        rc.unique_counts[1, 0] = 5.0

        em = _make_em_data([[0, 1]] * 100)
        rc.run_em(em, em_iterations=0)

        theta = rc._converged_theta
        # No warm-up or VBEM iterations → alpha = unique_totals + prior
        # alpha = [10+0.01, 5+0.01, 0+0.01, 0+1.0, 0+1.0]
        # Shadows use _SHADOW_PRIOR=1.0, transcripts use vbem_prior=0.01
        assert theta[0] > theta[1] > theta[2]
        assert theta is not None
        assert theta.sum() == pytest.approx(1.0)

    def test_likelihood_influences_em(self):
        """Candidates with higher likelihoods should attract more mass."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
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
        rc.run_em(em, em_iterations=100)
        theta = rc._converged_theta
        assert theta is not None


# =====================================================================
# assign_ambiguous — MAP assignment using EM posteriors
# =====================================================================


class TestAssignAmbiguous:
    def test_single_candidate_gets_full_count(self):
        """Unit with one candidate → count goes entirely to that transcript."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data([[0]])
        rc.run_em(em, em_iterations=1)
        rc.assign_ambiguous(em)

        assert rc.em_counts[0].sum() == 1.0
        assert rc.em_counts.sum() == 1.0

    def test_two_candidates_assigns_one(self):
        """Each unit contributes exactly 1.0 total count to one transcript."""
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
        np.testing.assert_array_almost_equal(
            rc.t_counts, rc.unique_counts + rc.em_counts
        )

    def test_distribution_follows_priors(self):
        """Over many fragments, sampled assignments follow unique priors."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, _UNSPLICED_SENSE] = 90.0
        rc.unique_counts[1, _UNSPLICED_SENSE] = 10.0

        em = _make_em_data([[0, 1]] * 10000)
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        t0_em = rc.em_counts[0].sum()
        t1_em = rc.em_counts[1].sum()
        assert t0_em > t1_em
        # Both should receive counts (not winner-take-all)
        assert t1_em > 500  # ~10% of 10000

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

    def test_count_col_respected(self):
        """Count goes to the correct column (category×strand)."""
        cc = int(CountCol.SPLICED_ANNOT_ANTISENSE)
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data(
            [[0]],
            count_cols_per_unit=[[cc]],
        )
        rc.run_em(em, em_iterations=1)
        rc.assign_ambiguous(em)

        assert rc.em_counts[0, cc] == 1.0
        assert rc.em_counts.sum() == 1.0

    def test_requires_run_em_first(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data([[0, 1]])
        with pytest.raises(RuntimeError, match="run_em"):
            rc.assign_ambiguous(em)

    def test_deterministic_with_same_seed(self):
        """Same seed → identical sampled counts."""
        results = []
        for _ in range(3):
            rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
            rc.unique_counts[0, _UNSPLICED_SENSE] = 50.0
            rc.unique_counts[1, _UNSPLICED_SENSE] = 50.0
            em = _make_em_data([[0, 1]] * 100)
            rc.run_em(em, em_iterations=10)
            rc.assign_ambiguous(em)
            results.append(rc.em_counts.copy())

        np.testing.assert_array_equal(results[0], results[1])
        np.testing.assert_array_equal(results[1], results[2])

    def test_different_seeds_differ(self):
        """Different seeds produce different sampled counts."""
        counts = []
        for seed in [1, 2]:
            rc = ReadCounter(num_transcripts=3, num_genes=2, seed=seed)
            rc.unique_counts[0, _UNSPLICED_SENSE] = 50.0
            rc.unique_counts[1, _UNSPLICED_SENSE] = 50.0
            em = _make_em_data([[0, 1]] * 100)
            rc.run_em(em, em_iterations=10)
            rc.assign_ambiguous(em)
            counts.append(rc.em_counts.copy())

        # Very unlikely to be identical with different seeds
        assert not np.array_equal(counts[0], counts[1])

    def test_integer_counts(self):
        """All EM counts are integer-valued (each fragment → 1 transcript)."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, _UNSPLICED_SENSE] = 50.0
        rc.unique_counts[1, _UNSPLICED_SENSE] = 30.0

        em = _make_em_data([[0, 1]] * 200)
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        # Every element should be an exact integer
        np.testing.assert_array_equal(rc.em_counts, np.floor(rc.em_counts))

    def test_equal_candidates_share_counts(self):
        """Equal-probability candidates all receive counts (not winner-take-all).

        With 4 candidates of equal probability (each ~25%), all should
        receive a substantial share over many fragments.
        """
        rc = ReadCounter(num_transcripts=4, num_genes=2, seed=42)
        em = _make_em_data([[0, 1, 2, 3]] * 10000)
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        per_t = rc.em_counts.sum(axis=1)
        # Each should get ~2500 counts; verify all get at least 2000
        for t in range(4):
            assert per_t[t] > 2000, (
                f"t{t} got {per_t[t]} counts — expected ~2500"
            )
        assert rc.em_counts.sum() == 10000.0


# =====================================================================
# Posterior mean
# =====================================================================


class TestPosteriorMean:
    def test_no_em_returns_nan(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        pm = rc.posterior_mean()
        assert np.all(np.isnan(pm))

    def test_single_candidate_posterior_one(self):
        """One candidate per unit → posterior is 1.0, mean is 1.0."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        em = _make_em_data([[0]] * 10)
        rc.run_em(em, em_iterations=1)
        rc.assign_ambiguous(em)

        pm = rc.posterior_mean()
        assert pm[0] == pytest.approx(1.0)
        assert np.isnan(pm[1])  # no assignments
        assert np.isnan(pm[2])

    def test_strong_prior_high_posterior(self):
        """Strong prior → assigned units have high mean posterior."""
        rc = ReadCounter(num_transcripts=2, num_genes=1, seed=42)
        rc.unique_counts[0, _UNSPLICED_SENSE] = 100.0
        em = _make_em_data([[0, 1]] * 100)
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        pm = rc.posterior_mean()
        # t0 gets most assignments with high posterior
        assert pm[0] > 0.8


# =====================================================================
# Simultaneous resolution — no phase ordering bias
# =====================================================================


class TestSimultaneousResolution:
    """Verify that isoform-ambig and gene-ambig fragments
    are resolved simultaneously, not sequentially."""

    def test_same_theta_regardless_of_order(self):
        """Shuffling ambiguous units produces identical EM theta.

        EM convergence is order-independent.  Sampled counts will
        differ (different draw sequences) but abundances agree.
        """
        rc1 = ReadCounter(num_transcripts=4, num_genes=2, seed=42)
        rc1.unique_counts[0, 0] = 50.0
        rc1.unique_counts[2, 0] = 30.0

        units = [[0, 1]] * 200 + [[0, 2]] * 200
        em1 = _make_em_data(units)
        rc1.run_em(em1, em_iterations=10)
        rc1.assign_ambiguous(em1)

        rc2 = ReadCounter(num_transcripts=4, num_genes=2, seed=42)
        rc2.unique_counts[0, 0] = 50.0
        rc2.unique_counts[2, 0] = 30.0

        units_rev = [[0, 2]] * 200 + [[0, 1]] * 200
        em2 = _make_em_data(units_rev)
        rc2.run_em(em2, em_iterations=10)
        rc2.assign_ambiguous(em2)

        np.testing.assert_allclose(
            rc1._converged_theta, rc2._converged_theta, atol=1e-10
        )
        # Total counts should match (400 units each)
        assert rc1.em_counts.sum() == rc2.em_counts.sum() == 400.0


# =====================================================================
# Multimapper-like units (multi-alignment molecules)
# =====================================================================


class TestMultimapperEM:
    def test_multimapper_molecule_one_count(self):
        """Multimapper molecule (many candidates) → exactly 1.0 total count."""
        rc = ReadCounter(num_transcripts=4, num_genes=2, seed=42)
        em = _make_em_data([[0, 1, 2, 3]])
        rc.run_em(em, em_iterations=5)
        rc.assign_ambiguous(em)

        assert rc.em_counts.sum() == 1.0

    def test_multimapper_and_ambig_together(self):
        """Mixed ambiguous + multimapper units in same EM."""
        rc = ReadCounter(num_transcripts=4, num_genes=2, seed=42)
        rc.unique_counts[0, 0] = 100.0

        units = [[0, 1]] * 50 + [[0, 2]] * 50 + [[0, 1, 2, 3]] * 25
        em = _make_em_data(units)
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        assert rc.em_counts.sum() == 125.0


# =====================================================================
# Primary counts DataFrame output
# =====================================================================


class TestCountsOutput:
    def test_get_counts_df_columns(self):
        index = _make_index()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        df = rc.get_counts_df(index)
        assert df.shape[0] == 3
        expected_cols = [
            "transcript_id", "gene_id", "gene_name",
            "count", "count_unique", "count_spliced",
            "count_em", "count_high_conf", "n_gdna",
            "posterior_mean",
        ]
        assert list(df.columns) == expected_cols

    def test_get_gene_counts_df_columns(self):
        index = _make_index()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        df = rc.get_gene_counts_df(index)
        assert df.shape[0] == 2
        expected_cols = [
            "gene_id", "gene_name",
            "count", "count_unique", "count_spliced",
            "count_em", "count_high_conf", "n_gdna",
            "n_antisense", "gdna_rate",
        ]
        assert list(df.columns) == expected_cols

    def test_counts_include_both_sources(self):
        """Total count sums unique + em counts."""
        index = _make_index()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, _UNSPLICED_SENSE] = 5.0
        rc.em_counts[0, _UNSPLICED_SENSE] = 3.0

        df = rc.get_counts_df(index)
        assert df.loc[0, "count"] == 8.0
        assert df.loc[0, "count_unique"] == 5.0
        assert df.loc[0, "count_em"] == 3.0

    def test_gene_counts_aggregate(self):
        """Gene counts aggregate across transcripts."""
        index = _make_index()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, _UNSPLICED_SENSE] = 5.0
        rc.unique_counts[1, _UNSPLICED_SENSE] = 3.0  # same gene
        rc.unique_counts[2, _UNSPLICED_SENSE] = 7.0  # different gene

        df = rc.get_gene_counts_df(index)
        assert df.loc[0, "count"] == 8.0  # g0 = t0 + t1
        assert df.loc[1, "count"] == 7.0  # g1 = t2

    def test_spliced_counts(self):
        """count_spliced captures both SPLICED_ANNOT and SPLICED_UNANNOT."""
        index = _make_index()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        # SPLICED_ANNOT sense
        rc.unique_counts[0, CountCol.SPLICED_ANNOT_SENSE] = 3.0
        # SPLICED_UNANNOT antisense
        rc.unique_counts[0, CountCol.SPLICED_UNANNOT_ANTISENSE] = 2.0
        # UNSPLICED (should not be in count_spliced)
        rc.unique_counts[0, _UNSPLICED_SENSE] = 10.0

        df = rc.get_counts_df(index)
        assert df.loc[0, "count_spliced"] == 5.0
        assert df.loc[0, "count"] == 15.0

    def test_identifiers_present(self):
        """Output includes transcript/gene identifiers."""
        index = _make_index()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        df = rc.get_counts_df(index)
        assert list(df["transcript_id"]) == ["t0", "t1", "t2"]
        assert list(df["gene_id"]) == ["g0", "g0", "g1"]


# =====================================================================
# Detail DataFrame output (long format)
# =====================================================================


class TestDetailOutput:
    def test_detail_empty(self):
        index = _make_index()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        df = rc.get_detail_df(index)
        assert len(df) == 0
        assert "transcript_id" in df.columns
        assert "category" in df.columns
        assert "source" in df.columns

    def test_detail_unique_only(self):
        index = _make_index()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, _UNSPLICED_SENSE] = 5.0
        rc.unique_counts[2, CountCol.SPLICED_ANNOT_SENSE] = 3.0

        df = rc.get_detail_df(index)
        assert len(df) == 2
        assert set(df["source"]) == {"unique"}

    def test_detail_both_sources(self):
        index = _make_index()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, _UNSPLICED_SENSE] = 5.0
        rc.em_counts[0, _UNSPLICED_SENSE] = 3.0

        df = rc.get_detail_df(index)
        assert len(df) == 2
        assert set(df["source"]) == {"unique", "em"}
        assert df["count"].sum() == 8.0

    def test_detail_has_category(self):
        index = _make_index()
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.unique_counts[0, _UNSPLICED_SENSE] = 1.0
        rc.unique_counts[0, CountCol.SPLICED_ANNOT_SENSE] = 1.0

        df = rc.get_detail_df(index)
        assert set(df["category"]) == {"unspliced", "spliced_annot"}


# =====================================================================
# gDNA summary output
# =====================================================================


class TestGDNASummaryOutput:
    def test_gdna_summary_dict(self):
        shadow = np.array([10.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
            shadow_init=shadow,
        )
        rc.gdna_em_counts[0] = 5.0
        rc.unique_counts[0, 0] = 80.0
        rc.em_counts[0, 0] = 5.0

        summary = rc.gdna_summary()
        assert summary["gdna_unique_count"] == 10.0
        assert summary["gdna_em_count"] == 5.0
        assert summary["gdna_total"] == 15.0
        assert summary["rna_unique_total"] == 80.0
        assert summary["rna_em_total"] == 5.0
        assert summary["gdna_contamination_rate"] == pytest.approx(0.15, abs=1e-6)
