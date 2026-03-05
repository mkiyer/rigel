"""Tests for hulkrna.estimator — AbundanceEstimator with locus-level EM assignment."""

import numpy as np
import pandas as pd
import pytest

from hulkrna.types import Strand, MergeOutcome
from hulkrna.splice import (
    SpliceType,
    SpliceStrandCol,
)
from hulkrna.config import EMConfig
from hulkrna.estimator import AbundanceEstimator, ScoredFragments, Locus, LocusEMInput, compute_hybrid_nrna_frac_priors, estimate_kappa
from hulkrna.strand_model import StrandModel, StrandModels

from conftest import _UNSPLICED_SENSE, _make_locus_em_data, _run_and_assign


# =====================================================================
# Helpers — lightweight mock index + fixtures
# =====================================================================


class MockIndex:
    """Minimal mock of TranscriptIndex with the arrays AbundanceEstimator needs."""

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
    """StrandModels container with FR-strand model."""
    sm = _make_strand_model_fr()
    return StrandModels(
        exonic_spliced=sm,
    )


def _make_resolved(**kwargs):
    from types import SimpleNamespace
    defaults = dict(
        t_inds=frozenset({0}),
        ambig_strand=0,
        splice_type=SpliceType.UNSPLICED,
        exon_strand=Strand.POS,
        sj_strand=Strand.NONE,
        frag_lengths={0: 250},
        merge_criteria=MergeOutcome.INTERSECTION,
        num_hits=1,
        genomic_footprint=250,
        genomic_start=1000,
    )
    defaults.update(kwargs)
    ns = SimpleNamespace(**defaults)
    ns.is_chimeric = False
    ns.is_same_strand = not ns.ambig_strand
    return ns


def _make_frag_length_models():
    from hulkrna.frag_length_model import FragmentLengthModels
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
# is_antisense
# =====================================================================


class TestIsAntisense:
    def test_high_likelihood_returns_false(self):
        """Good strand match → sense → not antisense."""
        sm = _make_strand_model_fr()
        result = AbundanceEstimator.is_antisense(Strand.POS, int(Strand.POS), sm)
        assert result is False

    def test_low_likelihood_returns_true(self):
        """Poor strand match → antisense."""
        sm = _make_strand_model_fr()
        result = AbundanceEstimator.is_antisense(Strand.POS, int(Strand.NEG), sm)
        assert result is True

    def test_unstranded_returns_false(self):
        """Unstranded model (p=0.5) is treated as non-antisense."""
        sm = StrandModel()  # no observations → p = 0.5
        result = AbundanceEstimator.is_antisense(Strand.POS, int(Strand.POS), sm)
        assert result is False


# =====================================================================
# assign_unambig
# =====================================================================


class TestAssignUnique:
    def test_single_transcript_single_gene(self):
        index = _make_index()
        sms = _make_strand_models_fr()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))

        resolved = _make_resolved(
            t_inds=frozenset({0}),
            ambig_strand=0,
            exon_strand=Strand.POS,
            splice_type=SpliceType.UNSPLICED,
        )
        rc.assign_unambig(resolved, index, sms)

        col = _UNSPLICED_SENSE
        assert rc.unambig_counts[0, col] == 1.0
        assert rc.t_counts[0, col] == 1.0

    def test_empty_t_inds_does_nothing(self):
        index = _make_index()
        sms = _make_strand_models_fr()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))

        resolved = _make_resolved(t_inds=frozenset(), ambig_strand=0)
        rc.assign_unambig(resolved, index, sms)

        assert rc.unambig_counts.sum() == 0.0

    def test_always_assigns_one_count(self):
        index = _make_index()
        sms = _make_strand_models_fr()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))

        resolved = _make_resolved(
            t_inds=frozenset({0}),
            ambig_strand=0,
            exon_strand=Strand.POS,
            num_hits=4,
        )
        rc.assign_unambig(resolved, index, sms)

        col = _UNSPLICED_SENSE
        assert rc.unambig_counts[0, col] == 1.0

    def test_writes_to_unambig_not_em(self):
        index = _make_index()
        sms = _make_strand_models_fr()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))

        resolved = _make_resolved(t_inds=frozenset({0}), ambig_strand=0)
        rc.assign_unambig(resolved, index, sms)

        assert rc.unambig_counts.sum() == 1.0
        assert rc.em_counts.sum() == 0.0


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
        """Empty LocusEMInput produces valid theta."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data([], num_transcripts=3)
        theta, alpha = rc.run_locus_em(lem, em_iterations=5)
        assert theta is not None
        assert theta.shape == (lem.n_components,)
        assert theta.sum() == pytest.approx(1.0)

    def test_uniform_priors_stay_uniform(self):
        """With no unambig counts and equal likelihoods, mRNA theta stays uniform."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data([[0, 1]] * 1000, num_transcripts=3)
        theta, _alpha = rc.run_locus_em(lem, em_iterations=10)

        # mRNA components 0 and 1 should be approximately equal
        assert theta[0] == pytest.approx(theta[1], rel=0.01)
        # t2 not in any unit → smaller
        assert theta[2] < theta[0]

    def test_unambig_counts_bias_em(self):
        """Unique counts on t0 should bias EM toward t0."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 100.0

        lem = _make_locus_em_data(
            [[0, 1]] * 500, num_transcripts=3, rc=rc,
        )
        theta, _alpha = rc.run_locus_em(lem, em_iterations=10)

        assert theta[0] > theta[1]

    def test_zero_iterations_uses_unambig_priors(self):
        """With em_iterations=0, theta comes from unambig counts + prior only."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 10.0
        rc.unambig_counts[1, 0] = 5.0

        lem = _make_locus_em_data(
            [[0, 1]] * 100, num_transcripts=3, rc=rc,
        )
        theta, _alpha = rc.run_locus_em(lem, em_iterations=0)

        # t0 gets more prior weight than t1
        assert theta[0] > theta[1]
        assert theta is not None
        assert theta.sum() == pytest.approx(1.0)

    def test_likelihood_influences_em(self):
        """Candidates with higher likelihoods should attract more mass."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data(
            [[0, 1]] * 1000,
            log_liks_per_unit=[[0.0, -10.0]] * 1000,
            num_transcripts=3,
        )
        theta, _alpha = rc.run_locus_em(lem, em_iterations=10)

        assert theta[0] > theta[1]

    def test_convergence_detected(self):
        """EM should converge early for simple problems."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 100.0
        rc.unambig_counts[1, 0] = 100.0

        lem = _make_locus_em_data(
            [[0, 1]] * 100, num_transcripts=2, rc=rc,
        )
        theta, _alpha = rc.run_locus_em(lem, em_iterations=100)
        assert theta is not None


# =====================================================================
# Locus assignment — posterior expected-count assignment
# =====================================================================


class TestLocusAssignment:
    def test_single_candidate_gets_full_count(self):
        """Unit with one candidate → count goes entirely to that transcript."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data([[0]], num_transcripts=3)
        _run_and_assign(rc, lem, em_iterations=1)

        assert rc.em_counts[0].sum() == pytest.approx(1.0)
        assert rc.em_counts.sum() == pytest.approx(1.0)

    def test_two_candidates_assigns_one(self):
        """Each unit contributes exactly 1.0 total count."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data([[0, 1]], num_transcripts=3)
        _run_and_assign(rc, lem, em_iterations=1)

        assert rc.em_counts.sum() == pytest.approx(1.0)

    def test_n_units_equals_n_counts(self):
        """N ambiguous units → N total em_counts."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data([[0, 1]] * 100, num_transcripts=3)
        _run_and_assign(rc, lem, em_iterations=5)

        assert rc.em_counts.sum() == pytest.approx(100.0)

    def test_total_counts_equals_unambig_plus_em(self):
        """t_counts = unambig_counts + em_counts."""
        index = _make_index()
        sms = _make_strand_models_fr()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))

        # 10 unambig
        for _ in range(10):
            rc.assign_unambig(
                _make_resolved(t_inds=frozenset({0})),
                index, sms,
            )

        # 5 ambiguous
        lem = _make_locus_em_data([[0, 1]] * 5, num_transcripts=3, rc=rc)
        _run_and_assign(rc, lem, em_iterations=5)

        assert rc.unambig_counts.sum() == 10.0
        assert rc.em_counts.sum() == pytest.approx(5.0)
        assert rc.t_counts.sum() == pytest.approx(15.0)
        np.testing.assert_array_almost_equal(
            rc.t_counts, rc.unambig_counts + rc.em_counts
        )

    def test_distribution_follows_priors(self):
        """Over many fragments, expected-count assignments follow unambig priors."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 90.0
        rc.unambig_counts[1, _UNSPLICED_SENSE] = 10.0

        lem = _make_locus_em_data(
            [[0, 1]] * 10000, num_transcripts=3, rc=rc,
        )
        _run_and_assign(rc, lem, em_iterations=10)

        t0_em = rc.em_counts[0].sum()
        t1_em = rc.em_counts[1].sum()
        assert t0_em > t1_em
        # Both should receive counts (not winner-take-all)
        assert t1_em > 500  # ~10% of 10000

    def test_writes_to_em_not_unambig(self):
        """assign_locus_ambiguous only writes to em_counts."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data([[0, 1]] * 50, num_transcripts=3)

        unique_before = rc.unambig_counts.copy()
        _run_and_assign(rc, lem, em_iterations=5)

        np.testing.assert_array_equal(rc.unambig_counts, unique_before)
        assert rc.em_counts.sum() == pytest.approx(50.0)

    def test_empty_em_data_does_nothing(self):
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data([], num_transcripts=3)
        _run_and_assign(rc, lem, em_iterations=5)

        assert rc.em_counts.sum() == 0.0

    def test_splice_strand_col_respected(self):
        """Count goes to the correct column (category×strand)."""
        cc = int(SpliceStrandCol.SPLICED_ANNOT_ANTISENSE)
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data(
            [[0]],
            count_cols_per_unit=[[cc]],
            num_transcripts=3,
        )
        _run_and_assign(rc, lem, em_iterations=1)

        assert rc.em_counts[0, cc] == pytest.approx(1.0)
        assert rc.em_counts.sum() == pytest.approx(1.0)

    def test_deterministic_with_same_seed(self):
        """Expected-count assignment is deterministic."""
        results = []
        for _ in range(3):
            rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
            rc.unambig_counts[0, _UNSPLICED_SENSE] = 50.0
            rc.unambig_counts[1, _UNSPLICED_SENSE] = 50.0
            lem = _make_locus_em_data(
                [[0, 1]] * 100, num_transcripts=3, rc=rc,
            )
            _run_and_assign(rc, lem, em_iterations=10)
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
            lem = _make_locus_em_data(
                [[0, 1]] * 100, num_transcripts=3, rc=rc,
            )
            _run_and_assign(rc, lem, em_iterations=10)
            counts.append(rc.em_counts.copy())

        np.testing.assert_array_equal(counts[0], counts[1])

    def test_fractional_counts(self):
        """EM counts may be fractional under expected-count assignment."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 50.0
        rc.unambig_counts[1, _UNSPLICED_SENSE] = 30.0

        lem = _make_locus_em_data(
            [[0, 1]] * 200, num_transcripts=3, rc=rc,
        )
        _run_and_assign(rc, lem, em_iterations=10)

        frac = rc.em_counts - np.floor(rc.em_counts)
        assert np.any(frac > 0.0)

    def test_equal_candidates_share_counts(self):
        """Equal-probability candidates all receive counts.

        With 4 candidates of equal probability (~25%), all should
        receive a substantial share over many fragments.
        """
        rc = AbundanceEstimator(4, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data(
            [[0, 1, 2, 3]] * 10000, num_transcripts=4,
        )
        _run_and_assign(rc, lem, em_iterations=10)

        per_t = rc.em_counts.sum(axis=1)
        # Each should get ~2500 counts; verify all get at least 2000
        for t in range(4):
            assert per_t[t] > 2000, (
                f"t{t} got {per_t[t]} counts — expected ~2500"
            )
        assert rc.em_counts.sum() == pytest.approx(10000.0)


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
        lem = _make_locus_em_data([[0]] * 10, num_transcripts=3)
        _run_and_assign(rc, lem, em_iterations=1)

        pm = rc.posterior_mean()
        assert pm[0] == pytest.approx(1.0)
        assert np.isnan(pm[1])  # no assignments
        assert np.isnan(pm[2])

    def test_strong_prior_high_posterior(self):
        """Strong prior → assigned units have high mean posterior."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 100.0
        lem = _make_locus_em_data(
            [[0, 1]] * 100, num_transcripts=2, rc=rc,
        )
        _run_and_assign(rc, lem, em_iterations=10)

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

        EM convergence is order-independent.  Expected counts agree.
        """
        rc1 = AbundanceEstimator(4, em_config=EMConfig(seed=42))
        rc1.unambig_counts[0, 0] = 50.0
        rc1.unambig_counts[2, 0] = 30.0

        units = [[0, 1]] * 200 + [[0, 2]] * 200
        lem1 = _make_locus_em_data(units, num_transcripts=4, rc=rc1)
        theta1, _ = _run_and_assign(rc1, lem1, em_iterations=10)

        rc2 = AbundanceEstimator(4, em_config=EMConfig(seed=42))
        rc2.unambig_counts[0, 0] = 50.0
        rc2.unambig_counts[2, 0] = 30.0

        units_rev = [[0, 2]] * 200 + [[0, 1]] * 200
        lem2 = _make_locus_em_data(units_rev, num_transcripts=4, rc=rc2)
        theta2, _ = _run_and_assign(rc2, lem2, em_iterations=10)

        np.testing.assert_allclose(theta1, theta2, atol=1e-10)
        # Total counts should match (400 units each)
        assert rc1.em_counts.sum() == pytest.approx(400.0)
        assert rc2.em_counts.sum() == pytest.approx(400.0)


# =====================================================================
# Multimapper-like units (multi-alignment molecules)
# =====================================================================


class TestMultimapperEM:
    def test_multimapper_molecule_one_count(self):
        """Multimapper molecule (many candidates) → exactly 1.0 total count."""
        rc = AbundanceEstimator(4, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data([[0, 1, 2, 3]], num_transcripts=4)
        _run_and_assign(rc, lem, em_iterations=5)

        assert rc.em_counts.sum() == pytest.approx(1.0)

    def test_multimapper_and_ambig_together(self):
        """Mixed ambiguous + multimapper units in same EM."""
        rc = AbundanceEstimator(4, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 100.0

        units = [[0, 1]] * 50 + [[0, 2]] * 50 + [[0, 1, 2, 3]] * 25
        lem = _make_locus_em_data(units, num_transcripts=4, rc=rc)
        _run_and_assign(rc, lem, em_iterations=10)

        assert rc.em_counts.sum() == pytest.approx(125.0, abs=1e-6)


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
            "transcript_id", "gene_id", "gene_name",
            "locus_id",
            "effective_length",
            "mrna", "mrna_unambig", "mrna_em",
            "mrna_high_conf", "mrna_spliced",
            "nrna", "rna_total", "tpm",
            "posterior_mean",
        ]
        assert list(df.columns) == expected_cols

    def test_get_gene_counts_df_columns(self):
        index = _make_index()
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        df = rc.get_gene_counts_df(index)
        assert df.shape[0] == 2
        expected_cols = [
            "gene_id", "gene_name",
            "locus_id",
            "effective_length",
            "mrna", "mrna_unambig", "mrna_em",
            "mrna_high_conf", "mrna_spliced",
            "nrna", "rna_total", "tpm",
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
        lem = _make_locus_em_data(
            [[0]] * 100,
            num_transcripts=2,
            include_gdna=True,
            gdna_init=500.0,
            gdna_log_lik=0.0,
        )
        theta, pool_counts = _run_and_assign(rc, lem, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        # gDNA should absorb some fragments
        assert gdna_count > 0

    def test_strong_rna_beats_gdna(self):
        """When transcript likelihood >> gDNA, most go to transcript."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, _UNSPLICED_SENSE] = 500.0
        lem = _make_locus_em_data(
            [[0]] * 100,
            log_liks_per_unit=[[0.0]] * 100,
            num_transcripts=2,
            rc=rc,
            include_gdna=True,
            gdna_init=1.0,
            gdna_log_lik=-20.0,
        )
        theta, pool_counts = _run_and_assign(rc, lem, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert rc.em_counts[0].sum() > 95
        assert gdna_count < 5

    def test_total_counts_preserved_with_gdna(self):
        """em_counts + nrna_em + gdna == n_units."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data(
            [[0, 1]] * 200,
            num_transcripts=2,
            include_nrna=True,
            include_gdna=True,
            gdna_init=50.0,
        )
        theta, pool_counts = _run_and_assign(rc, lem, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        total = rc.em_counts.sum() + rc.nrna_em_counts.sum() + gdna_count
        assert total == pytest.approx(200.0)

    def test_no_gdna_candidate_means_no_gdna_assignment(self):
        """Without gDNA candidate, all counts go to RNA."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data(
            [[0, 1]] * 100,
            num_transcripts=2,
            include_gdna=False,
        )
        theta, pool_counts = _run_and_assign(rc, lem, em_iterations=5)
        gdna_count = pool_counts["gdna"]

        assert gdna_count == 0.0
        assert rc.em_counts.sum() == pytest.approx(100.0)


# =====================================================================
# nrna_frac (nascent fraction) prior cascade tests
# =====================================================================


class TestEstimateKappa:
    """Tests for estimate_kappa() Method-of-Moments helper."""

    def test_zero_variance_returns_max(self):
        """All nrna_frac identical → σ²=0 → returns maximum κ."""
        nrna_frac = np.full(30, 0.3)
        ev = np.full(30, 50.0)
        k = estimate_kappa(nrna_frac, ev, min_evidence=20.0)
        assert k == 200.0  # _KAPPA_MAX

    def test_too_few_observations_returns_fallback(self):
        """Fewer than 20 valid features → returns fallback."""
        rng = np.random.default_rng(0)
        nrna_frac = rng.uniform(0.1, 0.9, size=5)
        ev = np.full(5, 100.0)
        k = estimate_kappa(nrna_frac, ev, min_evidence=20.0)
        assert k == 5.0  # _KAPPA_FALLBACK

    def test_known_beta_distribution(self):
        """Samples from Beta(2, 8) → μ=0.2, σ²≈0.0145 → κ≈10."""
        rng = np.random.default_rng(42)
        nrna_frac = rng.beta(2.0, 8.0, size=1000)
        ev = np.full(1000, 100.0)
        k = estimate_kappa(nrna_frac, ev, min_evidence=20.0)
        # True κ = α + β = 10
        assert k == pytest.approx(10.0, abs=1.5)

    def test_high_variance_gives_low_kappa(self):
        """Uniform(0,1) → σ²≈0.083 → κ close to 2 (our floor)."""
        rng = np.random.default_rng(7)
        nrna_frac = rng.uniform(0.0, 1.0, size=500)
        ev = np.full(500, 100.0)
        k = estimate_kappa(nrna_frac, ev, min_evidence=20.0)
        assert k == pytest.approx(2.0, abs=0.5)

    def test_evidence_filter_excludes_noisy(self):
        """Only features above min_evidence are used."""
        rng = np.random.default_rng(99)
        # 30 confident features with tight distribution → high κ
        nrna_frac_conf = np.full(30, 0.15)
        ev_conf = np.full(30, 100.0)
        # 200 noisy features with random nrna_frac but low evidence
        nrna_frac_noisy = rng.uniform(0.0, 1.0, size=200)
        ev_noisy = np.full(200, 5.0)
        nrna_frac = np.concatenate([nrna_frac_conf, nrna_frac_noisy])
        ev = np.concatenate([ev_conf, ev_noisy])
        k = estimate_kappa(nrna_frac, ev, min_evidence=50.0)
        # Should only see the 30 tight features → high κ
        assert k >= 100.0

    def test_clamp_lower_bound(self):
        """κ is clamped to at least 2.0."""
        nrna_frac = np.array([0.0, 1.0] * 20, dtype=np.float64)
        ev = np.full(40, 100.0)
        k = estimate_kappa(nrna_frac, ev, min_evidence=20.0)
        assert k >= 2.0


class TestEtaPriorShrinkage:
    """Tests for compute_hybrid_nrna_frac_priors() smooth EB shrinkage.

    The global prior nrna_frac_global is now *empirical* — the evidence-weighted
    mean of all locus-strand nrna_frac values.  When there is no evidence,
    nrna_frac_global = 0 (conservative: assume no nRNA).  Tests that check
    exact values pass explicit κ values.
    """

    KG, KL, KT = 2.0, 10.0, 5.0  # standard explicit kappas

    def _make_estimator(self, n_t, L_exonic=None, L_intronic=None):
        """Create a minimal AbundanceEstimator with geometry arrays."""
        rc = AbundanceEstimator(n_t)
        if L_exonic is not None:
            rc._exonic_lengths = np.array(L_exonic, dtype=np.float64)
        if L_intronic is not None:
            spans = rc._exonic_lengths + np.array(L_intronic, dtype=np.float64)
            rc._transcript_spans = spans
        return rc

    def _kw(self):
        return dict(kappa_global=self.KG, kappa_locus=self.KL, kappa_tss=self.KT)

    # -- Zero evidence → empirical global prior (nrna_frac=0) --

    def test_no_evidence_shrinks_to_global(self):
        """No data → nrna_frac_global=0 (empirical), κ=5 → α≈0, β≈5."""
        rc = self._make_estimator(3, [1000]*3, [4000]*3)
        locus_ids = np.full(3, -1, dtype=np.int32)
        strands = np.array([1, 1, 2], dtype=np.int32)
        compute_hybrid_nrna_frac_priors(
            rc, None, strands, locus_ids,
            strand_specificity=0.95, gdna_density=0.0, **self._kw(),
        )
        # nrna_frac_global=0 (no evidence), all nrna_frac→0, κ=5 → α≈0, β≈5
        np.testing.assert_allclose(rc.nrna_frac_alpha, 0.0, atol=0.01)
        np.testing.assert_allclose(rc.nrna_frac_beta, 5.0, atol=0.01)

    # -- Transcript level with shrinkage --

    def test_density_mode_transcript_level(self):
        """Unstranded (s=0.5): density nrna_frac shrunk toward empirical global."""
        rc = self._make_estimator(1, [1000], [4000])
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_SENSE] = 10.0
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_ANTISENSE] = 10.0
        rc.transcript_exonic_sense[0] = 10.0
        rc.transcript_exonic_antisense[0] = 10.0
        rc.transcript_intronic_sense[0] = 8.0
        rc.transcript_intronic_antisense[0] = 8.0

        locus_ids = np.full(1, -1, dtype=np.int32)
        strands = np.array([1], dtype=np.int32)
        compute_hybrid_nrna_frac_priors(
            rc, None, strands, locus_ids,
            strand_specificity=0.5, gdna_density=0.0, **self._kw(),
        )
        # Single transcript: kappa = kappa_tss = 5.0 (constant, weak prior)
        # nrna_frac ≈ 0.176 → α = 5 * 0.176 ≈ 0.878
        assert rc.nrna_frac_alpha[0] == pytest.approx(0.878, abs=0.05)
        kappa = rc.nrna_frac_alpha[0] + rc.nrna_frac_beta[0]
        assert kappa == pytest.approx(5.0, abs=0.01)

    def test_strand_mode_transcript_level(self):
        """Perfectly stranded (s=1.0): strand nrna_frac shrunk toward empirical global."""
        rc = self._make_estimator(1, [1000], [4000])
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_SENSE] = 20.0
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_ANTISENSE] = 0.0
        rc.transcript_exonic_sense[0] = 20.0
        rc.transcript_exonic_antisense[0] = 0.0
        rc.transcript_intronic_sense[0] = 8.0
        rc.transcript_intronic_antisense[0] = 0.0

        locus_ids = np.full(1, -1, dtype=np.int32)
        strands = np.array([1], dtype=np.int32)
        compute_hybrid_nrna_frac_priors(
            rc, None, strands, locus_ids,
            strand_specificity=1.0, gdna_density=0.0, **self._kw(),
        )
        # Single transcript: kappa = kappa_tss = 5.0 (constant, weak prior)
        # nrna_frac ≈ 0.085 → α = 5 * 0.085 ≈ 0.424
        assert rc.nrna_frac_alpha[0] == pytest.approx(0.424, abs=0.05)
        kappa = rc.nrna_frac_alpha[0] + rc.nrna_frac_beta[0]
        assert kappa == pytest.approx(5.0, abs=0.01)

    def test_gdna_subtraction_reduces_nrna(self):
        """Positive gDNA density reduces the nRNA estimate (directional)."""
        rc = self._make_estimator(1, [1000], [4000])
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_SENSE] = 10.0
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_ANTISENSE] = 10.0
        rc.transcript_exonic_sense[0] = 10.0
        rc.transcript_exonic_antisense[0] = 10.0
        rc.transcript_intronic_sense[0] = 8.0
        rc.transcript_intronic_antisense[0] = 8.0

        locus_ids = np.full(1, -1, dtype=np.int32)
        strands = np.array([1], dtype=np.int32)
        compute_hybrid_nrna_frac_priors(
            rc, None, strands, locus_ids,
            strand_specificity=0.5, gdna_density=0.0, **self._kw(),
        )
        alpha_no_gdna = rc.nrna_frac_alpha[0]

        # Fresh estimator for gDNA run (nrna_frac_global is re-derived each call)
        rc2 = self._make_estimator(1, [1000], [4000])
        rc2.unambig_counts[0, SpliceStrandCol.UNSPLICED_SENSE] = 10.0
        rc2.unambig_counts[0, SpliceStrandCol.UNSPLICED_ANTISENSE] = 10.0
        rc2.transcript_exonic_sense[0] = 10.0
        rc2.transcript_exonic_antisense[0] = 10.0
        rc2.transcript_intronic_sense[0] = 8.0
        rc2.transcript_intronic_antisense[0] = 8.0
        compute_hybrid_nrna_frac_priors(
            rc2, None, strands, locus_ids,
            strand_specificity=0.5, gdna_density=0.002, **self._kw(),
        )
        alpha_with_gdna = rc2.nrna_frac_alpha[0]

        # gDNA subtraction should lower nrna_frac and thus α
        assert alpha_with_gdna < alpha_no_gdna
        assert alpha_with_gdna == pytest.approx(0.488, abs=0.05)

    def test_single_exon_nrna_suppressed(self):
        """Single-exon transcript (L_intronic=0): nrna_frac pulled toward 0."""
        rc = self._make_estimator(1, [1000], [0])
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_SENSE] = 10.0
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_ANTISENSE] = 10.0
        rc.transcript_exonic_sense[0] = 10.0
        rc.transcript_exonic_antisense[0] = 10.0

        locus_ids = np.full(1, -1, dtype=np.int32)
        strands = np.array([1], dtype=np.int32)
        compute_hybrid_nrna_frac_priors(
            rc, None, strands, locus_ids,
            strand_specificity=0.5, gdna_density=0.001, **self._kw(),
        )
        nrna_frac = rc.nrna_frac_alpha[0] / (rc.nrna_frac_alpha[0] + rc.nrna_frac_beta[0])
        # Raw nrna_frac=0, nrna_frac_global=0; single-exon → nRNA strongly suppressed
        assert nrna_frac < 0.01

    # -- Hierarchical shrinkage --

    def test_tss_group_shrinkage(self):
        """Low-evidence transcript borrows strength from TSS group."""
        rc = self._make_estimator(3, [1000]*3, [4000]*3)
        # t0: low evidence (6 frags)
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_SENSE] = 2.0
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_ANTISENSE] = 2.0
        rc.transcript_exonic_sense[0] = 2.0
        rc.transcript_exonic_antisense[0] = 2.0
        rc.transcript_intronic_sense[0] = 1.0
        rc.transcript_intronic_antisense[0] = 1.0
        # t1: strong evidence (36 frags)
        rc.unambig_counts[1, SpliceStrandCol.UNSPLICED_SENSE] = 10.0
        rc.unambig_counts[1, SpliceStrandCol.UNSPLICED_ANTISENSE] = 10.0
        rc.transcript_exonic_sense[1] = 10.0
        rc.transcript_exonic_antisense[1] = 10.0
        rc.transcript_intronic_sense[1] = 8.0
        rc.transcript_intronic_antisense[1] = 8.0
        # t2: no evidence

        tss_groups = np.array([0, 0, 1], dtype=np.int32)
        locus_ids = np.full(3, -1, dtype=np.int32)
        strands = np.array([1, 1, 1], dtype=np.int32)
        compute_hybrid_nrna_frac_priors(
            rc, tss_groups, strands, locus_ids,
            strand_specificity=0.5, gdna_density=0.0, **self._kw(),
        )
        # t0: low evidence, borrows from TSS group; α≈0.685, κ=5
        assert rc.nrna_frac_alpha[0] == pytest.approx(0.685, abs=0.05)
        kappa0 = rc.nrna_frac_alpha[0] + rc.nrna_frac_beta[0]
        assert kappa0 == pytest.approx(5.0, abs=0.01)

        # t1: strong evidence, own data dominates; α≈0.970
        assert rc.nrna_frac_alpha[1] == pytest.approx(0.970, abs=0.05)

        # t2: zero evidence, TSS group 1 also zero → falls to global (nrna_frac≈0)
        assert rc.nrna_frac_alpha[2] < 0.01

    def test_locus_strand_shrinkage(self):
        """No TSS groups: locus-strand provides parent information."""
        rc = self._make_estimator(3, [1000]*3, [4000]*3)
        # t0: low evidence (6 frags)
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_SENSE] = 2.0
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_ANTISENSE] = 2.0
        rc.transcript_exonic_sense[0] = 2.0
        rc.transcript_exonic_antisense[0] = 2.0
        rc.transcript_intronic_sense[0] = 1.0
        rc.transcript_intronic_antisense[0] = 1.0
        # t1: no evidence
        # t2: strong evidence
        rc.unambig_counts[2, SpliceStrandCol.UNSPLICED_SENSE] = 10.0
        rc.unambig_counts[2, SpliceStrandCol.UNSPLICED_ANTISENSE] = 10.0
        rc.transcript_exonic_sense[2] = 10.0
        rc.transcript_exonic_antisense[2] = 10.0
        rc.transcript_intronic_sense[2] = 8.0
        rc.transcript_intronic_antisense[2] = 8.0

        locus_ids = np.array([0, 0, 0], dtype=np.int32)
        strands = np.array([1, 1, 1], dtype=np.int32)
        compute_hybrid_nrna_frac_priors(
            rc, None, strands, locus_ids,
            strand_specificity=0.5, gdna_density=0.0, **self._kw(),
        )
        # t0: low evidence, inherits from locus-strand; α≈0.767
        assert rc.nrna_frac_alpha[0] == pytest.approx(0.767, abs=0.05)

        # t1: zero own data → inherits parent (locus-strand nrna_frac); α≈0.938
        assert rc.nrna_frac_alpha[1] == pytest.approx(0.938, abs=0.05)

        # t2: strong evidence, own data dominates; α≈0.992
        assert rc.nrna_frac_alpha[2] == pytest.approx(0.992, abs=0.05)

    # -- Strand separation --

    def test_different_strands_separate_locus_pools(self):
        """Locus-strand groups transcripts by strand."""
        rc = self._make_estimator(2, [1000]*2, [4000]*2)
        # t0: +strand
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_SENSE] = 10.0
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_ANTISENSE] = 10.0
        rc.transcript_exonic_sense[0] = 10.0
        rc.transcript_exonic_antisense[0] = 10.0
        rc.transcript_intronic_sense[0] = 6.0
        rc.transcript_intronic_antisense[0] = 6.0
        # t1: -strand
        rc.unambig_counts[1, SpliceStrandCol.UNSPLICED_SENSE] = 15.0
        rc.unambig_counts[1, SpliceStrandCol.UNSPLICED_ANTISENSE] = 15.0
        rc.transcript_exonic_sense[1] = 15.0
        rc.transcript_exonic_antisense[1] = 15.0
        rc.transcript_intronic_sense[1] = 10.0
        rc.transcript_intronic_antisense[1] = 10.0

        locus_ids = np.array([0, 0], dtype=np.int32)
        strands = np.array([1, 2], dtype=np.int32)
        compute_hybrid_nrna_frac_priors(
            rc, None, strands, locus_ids,
            strand_specificity=0.5, gdna_density=0.0, **self._kw(),
        )
        nrna_frac0 = rc.nrna_frac_alpha[0] / (rc.nrna_frac_alpha[0] + rc.nrna_frac_beta[0])
        nrna_frac1 = rc.nrna_frac_alpha[1] / (rc.nrna_frac_alpha[1] + rc.nrna_frac_beta[1])
        # Raw nrna_frac0≈0.150, nrna_frac1≈0.167; both shrink slightly
        assert nrna_frac0 == pytest.approx(0.150, abs=0.01)
        assert nrna_frac1 == pytest.approx(0.167, abs=0.01)
        # κ = kappa_tss (constant, weak prior)
        kappa0 = rc.nrna_frac_alpha[0] + rc.nrna_frac_beta[0]
        kappa1 = rc.nrna_frac_alpha[1] + rc.nrna_frac_beta[1]
        assert kappa0 == pytest.approx(5.0, abs=0.01)
        assert kappa1 == pytest.approx(5.0, abs=0.01)

    # -- Shrinkage properties --

    def test_constant_kappa_regardless_of_evidence(self):
        """Kappa is constant (kappa_tss) regardless of evidence level."""
        # With the weakly informative prior, kappa = kappa_tss always.
        # The EM itself sees all data, so kappa need not scale with evidence.
        for ne in [5, 30, 100]:
            rc = self._make_estimator(1, [1000], [4000])
            rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_SENSE] = float(ne)
            rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_ANTISENSE] = float(ne)
            rc.transcript_exonic_sense[0] = float(ne)
            rc.transcript_exonic_antisense[0] = float(ne)
            rc.transcript_intronic_sense[0] = float(ne // 2)
            rc.transcript_intronic_antisense[0] = float(ne // 2)
            locus_ids = np.array([0], dtype=np.int32)
            strands = np.array([1], dtype=np.int32)
            compute_hybrid_nrna_frac_priors(
                rc, None, strands, locus_ids,
                strand_specificity=0.5, gdna_density=0.0, **self._kw(),
            )
            kappa = rc.nrna_frac_alpha[0] + rc.nrna_frac_beta[0]
            # κ is always kappa_tss, independent of fragment count
            assert kappa == pytest.approx(5.0, abs=0.01)
            # nrna_frac mean is valid
            nf = rc.nrna_frac_alpha[0] / kappa
            assert 0 < nf < 1

    # -- Auto MoM estimation (None defaults) --

    def test_auto_kappa_sparse_uses_fallback(self):
        """With < 20 features, MoM falls back to κ=5; no evidence → nrna_frac≈0."""
        rc = self._make_estimator(3, [1000]*3, [4000]*3)
        locus_ids = np.full(3, -1, dtype=np.int32)
        strands = np.array([1, 1, 2], dtype=np.int32)
        # No kappas → auto-estimate → sparse → all fallback to 5.0
        compute_hybrid_nrna_frac_priors(
            rc, None, strands, locus_ids,
            strand_specificity=0.5, gdna_density=0.0,
        )
        # With no evidence: nrna_frac_global=0, all κ=5 → α≈0, β≈5
        np.testing.assert_allclose(rc.nrna_frac_alpha, 0.0, atol=0.01)
        np.testing.assert_allclose(rc.nrna_frac_beta, 5.0, atol=0.01)

    def test_auto_kappa_produces_finite_output(self):
        """Auto-estimation produces valid finite priors regardless."""
        rc = self._make_estimator(1, [1000], [4000])
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_SENSE] = 10.0
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_ANTISENSE] = 10.0
        rc.transcript_exonic_sense[0] = 10.0
        rc.transcript_exonic_antisense[0] = 10.0
        rc.transcript_intronic_sense[0] = 8.0
        rc.transcript_intronic_antisense[0] = 8.0
        locus_ids = np.full(1, -1, dtype=np.int32)
        strands = np.array([1], dtype=np.int32)
        # Leave κ as None → MoM auto
        compute_hybrid_nrna_frac_priors(
            rc, None, strands, locus_ids,
            strand_specificity=0.5, gdna_density=0.0,
        )
        assert np.all(np.isfinite(rc.nrna_frac_alpha))
        assert np.all(np.isfinite(rc.nrna_frac_beta))
        assert np.all(rc.nrna_frac_alpha > 0)
        assert np.all(rc.nrna_frac_beta > 0)

    # -- Output types --

    def test_nrna_frac_alpha_beta_are_float64(self):
        """Output arrays should be float64."""
        rc = self._make_estimator(2, [1000]*2, [4000]*2)
        locus_ids = np.full(2, -1, dtype=np.int32)
        strands = np.array([1, 2], dtype=np.int32)
        compute_hybrid_nrna_frac_priors(
            rc, None, strands, locus_ids,
            strand_specificity=0.95, gdna_density=0.0,
        )
        assert rc.nrna_frac_alpha.dtype == np.float64
        assert rc.nrna_frac_beta.dtype == np.float64

    def test_no_geometry_falls_back_gracefully(self):
        """Estimator without geometry arrays still produces valid priors."""
        rc = AbundanceEstimator(2)
        rc.unambig_counts[0, SpliceStrandCol.UNSPLICED_SENSE] = 20.0
        rc.transcript_exonic_sense[0] = 20.0
        locus_ids = np.full(2, -1, dtype=np.int32)
        strands = np.array([1, 1], dtype=np.int32)
        compute_hybrid_nrna_frac_priors(
            rc, None, strands, locus_ids,
            strand_specificity=0.5, gdna_density=0.0,
        )
        assert np.all(np.isfinite(rc.nrna_frac_alpha))
        assert np.all(np.isfinite(rc.nrna_frac_beta))
