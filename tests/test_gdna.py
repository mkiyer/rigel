"""Tests for gDNA contamination modeling in EM.

Verifies that the gDNA pseudo-component correctly competes with
transcripts in the EM, and that gDNA assignments are tracked
in the separate gdna_em_count / gdna_locus_counts arrays.
"""

import numpy as np
import pytest

from hulkrna.categories import CountCategory, CountType, NUM_COUNT_TYPES
from hulkrna.counter import ReadCounter, EMData
from hulkrna.pipeline import _score_gdna_candidate, _GDNA_SPLICE_PENALTIES
from hulkrna.insert_model import InsertSizeModels
from hulkrna.strand_model import StrandModel, StrandModels
from hulkrna.types import Strand


# =====================================================================
# Helpers
# =====================================================================


def _make_insert_models(n_obs=200, size=250):
    """InsertSizeModels with observations for intergenic + unspliced."""
    im = InsertSizeModels()
    for _ in range(n_obs):
        im.observe(size, CountCategory.UNSPLICED)
        im.observe(size, None)  # intergenic (count_cat=None)
    return im


def _make_strand_models_default():
    """StrandModels with untrained intergenic model (P_strand=0.5)."""
    return StrandModels()


def _make_em_data_with_gdna(
    t_indices_per_unit,
    log_liks_per_unit,
    count_types_per_unit=None,
    gdna_base_index=None,
):
    """Build EMData where some candidates are gDNA (using gdna_index).

    Like the test_counter helper, but explicitly supports gDNA index
    in the candidate lists.
    """
    offsets = [0]
    flat_t = []
    flat_lk = []
    flat_ct = []

    for u, t_list in enumerate(t_indices_per_unit):
        for j, t_idx in enumerate(t_list):
            flat_t.append(t_idx)
            flat_lk.append(log_liks_per_unit[u][j])
            if count_types_per_unit is not None:
                flat_ct.append(count_types_per_unit[u][j])
            else:
                flat_ct.append(int(CountType.UNSPLICED_SENSE))
        offsets.append(len(flat_t))

    n_units = len(t_indices_per_unit)
    n_candidates = len(flat_t)

    if gdna_base_index is None:
        gdna_base_index = (max(flat_t) + 1) if flat_t else 0

    # Build locus tracking: best non-gDNA candidate per unit
    locus_t = np.full(n_units, -1, dtype=np.int32)
    locus_ct = np.zeros(n_units, dtype=np.uint8)
    for u, t_list in enumerate(t_indices_per_unit):
        ct_list = count_types_per_unit[u] if count_types_per_unit else None
        best_ll = -np.inf
        for j, t_idx in enumerate(t_list):
            if t_idx < gdna_base_index and log_liks_per_unit[u][j] > best_ll:
                best_ll = log_liks_per_unit[u][j]
                locus_t[u] = t_idx
                locus_ct[u] = (
                    ct_list[j] if ct_list else int(CountType.UNSPLICED_SENSE)
                )

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
# _score_gdna_candidate
# =====================================================================


class TestScoreGDNA:
    """Tests for pipeline._score_gdna_candidate."""

    def test_unspliced_no_penalty(self):
        """UNSPLICED gets splice_penalty=1.0 → log(1)=0."""
        im = _make_insert_models()
        sms = _make_strand_models_default()
        score = _score_gdna_candidate(
            int(Strand.POS), int(CountCategory.UNSPLICED), 250,
            sms, im,
        )
        # log(0.5) + log(P_insert) + log(1.0)
        assert score < 0  # log(0.5) is negative
        assert score == pytest.approx(
            np.log(0.5) + im.intergenic.log_likelihood(250)
        )

    def test_intron_no_penalty(self):
        """INTRON gets splice_penalty=1.0."""
        im = _make_insert_models()
        sms = _make_strand_models_default()
        score = _score_gdna_candidate(
            int(Strand.POS), int(CountCategory.INTRON), 250,
            sms, im,
        )
        expected = np.log(0.5) + im.intergenic.log_likelihood(250)
        assert score == pytest.approx(expected)

    def test_spliced_unannot_penalty(self):
        """SPLICED_UNANNOT gets heavy penalty (default 0.01)."""
        im = _make_insert_models()
        sms = _make_strand_models_default()
        score_unannot = _score_gdna_candidate(
            int(Strand.POS), int(CountCategory.SPLICED_UNANNOT), 250,
            sms, im,
        )
        score_unspliced = _score_gdna_candidate(
            int(Strand.POS), int(CountCategory.UNSPLICED), 250,
            sms, im,
        )
        # SPLICED_UNANNOT should have a much lower score
        assert score_unannot < score_unspliced
        penalty_diff = score_unspliced - score_unannot
        assert penalty_diff == pytest.approx(-np.log(0.01), rel=0.01)

    def test_custom_penalties(self):
        """Custom penalty dict overrides defaults."""
        im = _make_insert_models()
        sms = _make_strand_models_default()
        custom = {CountCategory.UNSPLICED: 0.5}
        score = _score_gdna_candidate(
            int(Strand.POS), int(CountCategory.UNSPLICED), 250,
            sms, im,
            gdna_splice_penalties=custom,
        )
        expected = np.log(0.5) + im.intergenic.log_likelihood(250) + np.log(0.5)
        assert score == pytest.approx(expected)

    def test_zero_insert_size(self):
        """insert_size=0 → P_insert=1.0 (log=0)."""
        im = _make_insert_models()
        sms = _make_strand_models_default()
        score = _score_gdna_candidate(
            int(Strand.POS), int(CountCategory.UNSPLICED), 0,
            sms, im,
        )
        # log(0.5) + 0 + log(1.0) = log(0.5)
        assert score == pytest.approx(np.log(0.5))


# =====================================================================
# gDNA theta vector and EM
# =====================================================================


class TestGDNAThetaVector:
    """Verify per-gene gDNA shadows are properly represented in EM theta."""

    def test_theta_has_shadow_elements(self):
        """Theta has N_t + N_g elements (per-gene shadows)."""
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        assert rc.gdna_base_index == 3

        # Shadow for gene 0 at index 3, gene 1 at index 4
        em = _make_em_data_with_gdna(
            [[0, 3]], [[0.0, 0.0]], gdna_base_index=3,
        )
        rc.run_em(em, em_iterations=1)
        assert rc._converged_theta.shape == (5,)  # 3 transcripts + 2 genes
        assert rc._converged_theta.sum() == pytest.approx(1.0)

    def test_shadow_init_feeds_prior(self):
        """shadow_init biases the per-gene gDNA prior."""
        shadow = np.array([100.0])  # 1 gene
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
            shadow_init=shadow,
        )
        # Shadow index for gene 0 = gdna_base_index + 0 = 2
        em = _make_em_data_with_gdna(
            [[0, 2]] * 100, [[0.0, 0.0]] * 100, gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)

        theta = rc._converged_theta
        # gDNA shadow has 100 shadow_init; t0 has alpha=1 only
        assert theta[2] > theta[0], (
            f"gDNA theta ({theta[2]:.4f}) should dominate over t0 ({theta[0]:.4f})"
        )

    def test_no_shadow_init_minimal_absorption(self):
        """With zero shadow_init, gDNA absorbs only alpha."""
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
        )
        rc.unique_counts[0, 0] = 1000.0

        em = _make_em_data_with_gdna(
            [[0, 2]] * 100, [[0.0, -5.0]] * 100, gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)

        theta = rc._converged_theta
        # t0 has 1000 unique counts; shadow has only alpha=1 and worse likelihood
        assert theta[2] < 0.01, f"gDNA theta too high: {theta[2]:.4f}"


# =====================================================================
# gDNA assignment routing
# =====================================================================


class TestGDNAAssignment:
    """Verify that gDNA assignments go to gdna_em_counts, not em_counts."""

    def test_gdna_assignment_not_in_em_counts(self):
        """Fragments assigned to gDNA are NOT in em_counts."""
        # 2 transcripts, 1 gene → shadow index = 2
        shadow_idx = 2
        shadow = np.array([1000.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
            shadow_init=shadow,
        )
        # 100 units, each with t0 (bad lik) and shadow (good lik)
        em = _make_em_data_with_gdna(
            [[0, shadow_idx]] * 100,
            [[-10.0, 0.0]] * 100,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        # Most should go to gDNA
        assert rc.gdna_em_count > 0
        total_assigned = rc.em_counts.sum() + rc.gdna_em_count
        assert total_assigned == 100.0

    def test_total_counts_preserved(self):
        """em_counts + gdna_em_count == n_units for all assigned fragments."""
        # 3 transcripts, 2 genes → shadow indices 3, 4
        shadow = np.array([25.0, 25.0])
        rc = ReadCounter(
            num_transcripts=3, num_genes=2, seed=42,
            shadow_init=shadow,
        )
        # Use shadow for gene 0 (index 3)
        em = _make_em_data_with_gdna(
            [[0, 1, 3]] * 200,
            [[0.0, 0.0, 0.0]] * 200,
            gdna_base_index=3,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        assert rc.em_counts.sum() + rc.gdna_em_count == 200.0

    def test_strong_transcript_signal_beats_gdna(self):
        """When transcript likelihood >> gDNA, most go to transcript."""
        # 2 transcripts, 1 gene → shadow index = 2
        shadow_idx = 2
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
        )
        rc.unique_counts[0, 0] = 500.0

        # t0 has log_lik=0.0, shadow has log_lik=-20.0 (terrible)
        em = _make_em_data_with_gdna(
            [[0, shadow_idx]] * 500,
            [[0.0, -20.0]] * 500,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        # Almost all should go to t0
        assert rc.em_counts[0].sum() > 490
        assert rc.gdna_em_count < 10

    def test_high_shadow_prior_absorbs_fragments(self):
        """With large shadow_init, gDNA takes larger share."""
        # 2 transcripts, 1 gene → shadow index = 2
        shadow_idx = 2
        shadow = np.array([500.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
            shadow_init=shadow,
        )
        # Equal likelihoods — prior drives assignment
        em = _make_em_data_with_gdna(
            [[0, shadow_idx]] * 500,
            [[0.0, 0.0]] * 500,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        # gDNA should get a meaningful share
        assert rc.gdna_em_count > 100, (
            f"Expected gDNA to absorb many fragments, got {rc.gdna_em_count}"
        )


# =====================================================================
# gDNA locus attribution
# =====================================================================


class TestGDNALocusAttribution:
    """When a fragment goes to gDNA shadow, its locus transcript gets credited."""

    def test_locus_counts_populated(self):
        """gDNA assignments populate gdna_locus_counts."""
        # 2 transcripts, 1 gene → shadow index = 2
        shadow_idx = 2
        shadow = np.array([1000.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
            shadow_init=shadow,
        )
        ct_col = int(CountType.UNSPLICED_SENSE)
        # All units: t0 as transcript candidate, shadow with better score
        em = _make_em_data_with_gdna(
            [[0, shadow_idx]] * 100,
            [[-10.0, 0.0]] * 100,
            count_types_per_unit=[[ct_col, ct_col]] * 100,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        # locus_t_indices point to t0, so gdna_locus_counts[0] should be nonzero
        assert rc.gdna_locus_counts[0].sum() > 0

    def test_locus_counts_zero_when_no_gdna(self):
        """No gDNA assignments → locus counts stay zero."""
        # 2 transcripts, 1 gene → shadow index = 2
        shadow_idx = 2
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
        )
        rc.unique_counts[0, 0] = 1000.0
        # shadow has terrible likelihood → none assigned
        em = _make_em_data_with_gdna(
            [[0, shadow_idx]] * 50,
            [[0.0, -50.0]] * 50,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        assert rc.gdna_locus_counts.sum() == 0.0


# =====================================================================
# gDNA contamination rate and properties
# =====================================================================


class TestGDNAProperties:
    def test_gdna_total_includes_unique_and_em(self):
        shadow = np.array([50.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
            shadow_init=shadow,
        )
        rc.gdna_em_counts[0] = 30.0
        assert rc.gdna_total == 80.0

    def test_gdna_contamination_rate(self):
        shadow = np.array([10.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
            shadow_init=shadow,
        )
        rc.unique_counts[0, 0] = 80.0
        rc.em_counts[0, 0] = 10.0
        # total RNA = 80 + 10 = 90; gDNA_total = 10 + 0 = 10; total = 100
        assert rc.gdna_contamination_rate == pytest.approx(0.1)

    def test_contamination_rate_zero_when_empty(self):
        rc = ReadCounter(num_transcripts=2, num_genes=1, seed=42)
        assert rc.gdna_contamination_rate == 0.0


# =====================================================================
# gDNA output methods
# =====================================================================


class TestGDNAOutput:
    def test_get_t_gdna_df_shape(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        df = rc.get_t_gdna_df()
        assert df.shape == (3, NUM_COUNT_TYPES)
        assert list(df.columns) == list(CountType.columns())

    def test_get_g_gdna_df_aggregates(self):
        rc = ReadCounter(num_transcripts=3, num_genes=2, seed=42)
        rc.gdna_locus_counts[0, 0] = 5.0
        rc.gdna_locus_counts[1, 0] = 3.0  # same gene as t0
        rc.gdna_locus_counts[2, 0] = 7.0  # different gene

        t_to_g = np.array([0, 0, 1], dtype=np.int64)
        df = rc.get_g_gdna_df(t_to_g)
        assert df.shape == (2, NUM_COUNT_TYPES)
        assert df.iloc[0, 0] == 8.0  # g0 = t0 + t1
        assert df.iloc[1, 0] == 7.0  # g1 = t2

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


# =====================================================================
# Stats integration
# =====================================================================


class TestStatsGDNA:
    """Verify PipelineStats gDNA field changes."""

    def test_intergenic_split(self):
        from hulkrna.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 100
        stats.n_intergenic_spliced = 10
        assert stats.n_intergenic == 110

    def test_gdna_unique_property(self):
        from hulkrna.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 80
        stats.n_intergenic_spliced = 20
        assert stats.n_gdna_unique == 100

    def test_gdna_total_property(self):
        from hulkrna.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 80
        stats.n_intergenic_spliced = 20
        stats.n_gdna_em = 50
        assert stats.n_gdna_total == 150

    def test_to_dict_includes_gdna(self):
        from hulkrna.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 10
        stats.n_intergenic_spliced = 5
        stats.n_gdna_em = 3

        d = stats.to_dict()
        assert d["n_intergenic"] == 15
        assert d["n_gdna_unique"] == 15
        assert d["n_gdna_total"] == 18


# =====================================================================
# compute_shadow_init (empirical Bayes shrinkage)
# =====================================================================


class TestComputeShadowInit:
    """Tests for pipeline.compute_shadow_init."""

    def test_zero_counts_returns_zeros(self):
        """No counts at all → all shadow_init = 0."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((3, 12), dtype=np.float64)
        t_to_g = np.array([0, 0, 1], dtype=np.int64)
        result = compute_shadow_init(unique_counts, t_to_g, 2, 0)
        np.testing.assert_array_equal(result, [0.0, 0.0])

    def test_no_genic_reads_gives_zeros(self):
        """With only intergenic fragments (no genic reads), shadow = 0."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((2, 12), dtype=np.float64)
        t_to_g = np.array([0, 1], dtype=np.int64)
        # 100 intergenic fragments but zero genic reads → depth=0 → shadow=0
        result = compute_shadow_init(unique_counts, t_to_g, 2, 100)
        np.testing.assert_array_equal(result, [0.0, 0.0])

    def test_intergenic_with_gene_depth_sets_floor(self):
        """Genes with reads get a depth-scaled floor even without antisense."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((2, 12), dtype=np.float64)
        t_to_g = np.array([0, 1], dtype=np.int64)
        # Both genes have 100 sense reads, plus 100 intergenic
        unique_counts[0, 1] = 100.0  # INTRON_SENSE
        unique_counts[1, 1] = 100.0  # INTRON_SENSE
        result = compute_shadow_init(unique_counts, t_to_g, 2, 100, kappa=1.0)
        # n_total = 200 + 100 = 300
        # estimated_gdna = 100 + 0 = 100
        # θ_global = 100 / 300 = 1/3
        # depth = [100, 100]
        # shadow = [0 + 1.0*100*(1/3), 0 + 1.0*100*(1/3)] = [33.33, 33.33]
        assert result[0] == pytest.approx(100.0 / 3)
        assert result[1] == pytest.approx(100.0 / 3)

    def test_antisense_boosts_gene_shadow(self):
        """Genes with antisense counts get proportionally higher shadow."""
        from hulkrna.pipeline import compute_shadow_init
        from hulkrna.categories import CountStrand

        unique_counts = np.zeros((3, 12), dtype=np.float64)
        t_to_g = np.array([0, 0, 1], dtype=np.int64)
        # Gene 0: 50 antisense via t0
        as_col = int(CountStrand.ANTISENSE)  # col 2
        unique_counts[0, as_col] = 50.0
        # Gene 1: 50 sense reads via t2 (to give it depth)
        unique_counts[2, 1] = 50.0  # INTRON_SENSE
        result = compute_shadow_init(unique_counts, t_to_g, 2, 0, kappa=1.0)
        # n_total = 100, estimated_gdna = 0 + 2*50 = 100
        # θ_global = 100 / 100 = 1.0
        # depth = [50, 50], g_antisense = [50, 0]
        # gene 0: 2*50 + 1.0*50*1.0 = 150.0
        # gene 1: 0 + 1.0*50*1.0 = 50.0
        assert result[0] == pytest.approx(150.0)
        assert result[1] == pytest.approx(50.0)
        assert result[0] > result[1]  # antisense evidence boosts gene 0

    def test_kappa_scales_global_shrinkage(self):
        """Higher kappa → more shrinkage toward global rate."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((2, 12), dtype=np.float64)
        t_to_g = np.array([0, 1], dtype=np.int64)
        # Give genes depth so the global term is non-zero
        unique_counts[0, 1] = 100.0  # INTRON_SENSE
        unique_counts[1, 1] = 100.0  # INTRON_SENSE
        r_lo = compute_shadow_init(unique_counts, t_to_g, 2, 100, kappa=0.5)
        r_hi = compute_shadow_init(unique_counts, t_to_g, 2, 100, kappa=5.0)
        assert r_hi[0] > r_lo[0]

    def test_depth_scaling(self):
        """Deeper gene gets larger absolute shadow floor."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((2, 12), dtype=np.float64)
        t_to_g = np.array([0, 1], dtype=np.int64)
        # Gene 0: 1000 sense reads (deeply sequenced)
        # Gene 1: 10 sense reads (lowly sequenced)
        unique_counts[0, 1] = 1000.0
        unique_counts[1, 1] = 10.0
        result = compute_shadow_init(unique_counts, t_to_g, 2, 100, kappa=1.0)
        # Both have 0 antisense, so shadow comes from depth*θ_global
        # θ_global = 100 / 1110 ≈ 0.09
        # gene 0: 0 + 1.0*1000*0.09 ≈ 90
        # gene 1: 0 + 1.0*10*0.09 ≈ 0.9
        assert result[0] > result[1] * 50  # ~100x ratio
        assert result[0] > 0
        assert result[1] > 0

    def test_output_shape(self):
        """Returns array of shape (num_genes,)."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((5, 12), dtype=np.float64)
        t_to_g = np.array([0, 0, 1, 2, 2], dtype=np.int64)
        result = compute_shadow_init(unique_counts, t_to_g, 3, 10)
        assert result.shape == (3,)
        assert result.dtype == np.float64


# =====================================================================
# Per-gene shadow EM integration
# =====================================================================


class TestPerGeneShadowEM:
    """Integration tests for per-gene shadow behavior in EM."""

    def test_two_genes_independent_shadows(self):
        """Two genes get independent shadow components in theta."""
        # 3 transcripts (t0, t1 → gene 0; t2 → gene 1), 2 genes
        # Shadow indices: 3 (gene 0), 4 (gene 1)
        shadow = np.array([100.0, 0.0])  # gene 0 has gDNA evidence
        rc = ReadCounter(
            num_transcripts=3, num_genes=2, seed=42,
            shadow_init=shadow,
        )
        # Units: each competes t0 vs shadow_gene0 (index 3)
        em = _make_em_data_with_gdna(
            [[0, 3]] * 200,
            [[0.0, 0.0]] * 200,
            gdna_base_index=3,
        )
        rc.run_em(em, em_iterations=10)

        theta = rc._converged_theta
        assert theta.shape == (5,)
        # Shadow for gene 0 (index 3) should have significant mass
        assert theta[3] > 0.1
        # Shadow for gene 1 (index 4) should have minimal mass (only alpha)
        assert theta[4] < theta[3]

    def test_per_gene_em_counts_tracking(self):
        """Verify gdna_em_counts tracks per-gene gDNA assignments."""
        # 2 transcripts, 2 genes → shadow indices 2, 3
        shadow = np.array([500.0, 500.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=2, seed=42,
            shadow_init=shadow,
        )
        # 200 units: t0 vs shadow_gene0, 200 units: t1 vs shadow_gene1
        em = _make_em_data_with_gdna(
            [[0, 2]] * 200 + [[1, 3]] * 200,
            [[0.0, 0.0]] * 200 + [[0.0, 0.0]] * 200,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        # Both genes should have some gDNA assignments
        assert rc.gdna_em_counts[0] > 0
        assert rc.gdna_em_counts[1] > 0
        # Total should be preserved
        total = rc.em_counts.sum() + rc.gdna_em_counts.sum()
        assert total == 400.0

    def test_gene_with_high_antisense_absorbs_more(self):
        """Gene with high shadow_init absorbs more ambiguous fragments."""
        # 2 transcripts (1 per gene), 2 genes → shadow indices 2, 3
        shadow = np.array([200.0, 1.0])  # gene 0 has high gDNA evidence
        rc = ReadCounter(
            num_transcripts=2, num_genes=2, seed=42,
            shadow_init=shadow,
        )
        rc.unique_counts[0, 0] = 10.0
        rc.unique_counts[1, 0] = 10.0

        # 500 units each: t0 vs shadow_g0, t1 vs shadow_g1
        em = _make_em_data_with_gdna(
            [[0, 2]] * 500 + [[1, 3]] * 500,
            [[0.0, 0.0]] * 500 + [[0.0, 0.0]] * 500,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        # Gene 0's shadow should absorb more than gene 1's
        assert rc.gdna_em_counts[0] > rc.gdna_em_counts[1]


# =====================================================================
# StrandModels container
# =====================================================================


class TestStrandModelsContainer:
    """Tests for the StrandModels multi-region container."""

    def test_default_construction(self):
        from hulkrna.strand_model import StrandModels

        sm = StrandModels()
        assert sm.exonic_spliced.n_observations == 0
        assert sm.exonic.n_observations == 0
        assert sm.intronic.n_observations == 0
        assert sm.intergenic.n_observations == 0

    def test_delegation_to_exonic_spliced(self):
        from hulkrna.strand_model import StrandModels
        from hulkrna.types import Strand

        sm = StrandModels()
        for _ in range(100):
            sm.exonic_spliced.observe(Strand.POS, Strand.POS)
        for _ in range(2):
            sm.exonic_spliced.observe(Strand.POS, Strand.NEG)

        assert sm.p_r1_sense == sm.exonic_spliced.p_r1_sense
        assert sm.strand_specificity == sm.exonic_spliced.strand_specificity
        assert sm.read1_sense == sm.exonic_spliced.read1_sense
        assert sm.n_observations == sm.exonic_spliced.n_observations
        assert sm.strand_likelihood(Strand.POS, Strand.POS) == sm.exonic_spliced.strand_likelihood(Strand.POS, Strand.POS)

    def test_to_dict_has_all_regions(self):
        from hulkrna.strand_model import StrandModels

        sm = StrandModels()
        d = sm.to_dict()
        assert "exonic_spliced" in d
        assert "exonic" in d
        assert "intronic" in d
        assert "intergenic" in d
        assert "observations" in d["exonic_spliced"]

    def test_write_json(self, tmp_path):
        from hulkrna.strand_model import StrandModels
        from hulkrna.types import Strand
        import json

        sm = StrandModels()
        for _ in range(50):
            sm.exonic_spliced.observe(Strand.POS, Strand.POS)

        out = tmp_path / "strand_models.json"
        sm.write_json(out)
        with open(out) as f:
            data = json.load(f)
        assert "strand_models" in data
        assert "exonic_spliced" in data["strand_models"]
        assert "exonic" in data["strand_models"]
        assert "intronic" in data["strand_models"]
        assert "intergenic" in data["strand_models"]

    def test_independent_models(self):
        """Each region's model is independent."""
        from hulkrna.strand_model import StrandModels
        from hulkrna.types import Strand

        sm = StrandModels()
        sm.exonic_spliced.observe(Strand.POS, Strand.POS)
        sm.exonic.observe(Strand.POS, Strand.NEG)
        sm.intronic.observe(Strand.NEG, Strand.POS)
        sm.intergenic.observe(Strand.NEG, Strand.NEG)

        assert sm.exonic_spliced.n_observations == 1
        assert sm.exonic.n_observations == 1
        assert sm.intronic.n_observations == 1
        assert sm.intergenic.n_observations == 1

    def test_model_for_category(self):
        """model_for_category selects the right sub-model."""
        from hulkrna.strand_model import StrandModels
        from hulkrna.categories import CountCategory

        sm = StrandModels()
        assert sm.model_for_category(int(CountCategory.SPLICED_ANNOT)) is sm.exonic_spliced
        assert sm.model_for_category(int(CountCategory.SPLICED_UNANNOT)) is sm.exonic
        assert sm.model_for_category(int(CountCategory.UNSPLICED)) is sm.exonic
        assert sm.model_for_category(int(CountCategory.INTRON)) is sm.intronic
