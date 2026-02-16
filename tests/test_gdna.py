"""Tests for gDNA contamination modeling in EM.

Verifies that the gDNA pseudo-component correctly competes with
transcripts in the EM, and that gDNA assignments are tracked
in the separate gdna_em_count / gdna_locus_counts arrays.
"""

import numpy as np
import pytest

from hulkrna.categories import CountCategory, CountCol, NUM_COUNT_COLS
from hulkrna.counter import ReadCounter, EMData
from hulkrna.pipeline import _score_gdna_candidate, _GDNA_SPLICE_PENALTIES
from hulkrna.insert_model import InsertSizeModels
from hulkrna.strand_model import StrandModel, StrandModels
from hulkrna.types import Strand


# Default column for UNSPLICED_SENSE
_UNSPLICED_SENSE = int(CountCol.UNSPLICED_SENSE)


# =====================================================================
# Helpers
# =====================================================================


def _make_insert_models(n_obs=200, size=250):
    """InsertSizeModels with observations for all categories + intergenic."""
    im = InsertSizeModels()
    for _ in range(n_obs):
        im.observe(size, CountCategory.SPLICED_ANNOT)
        im.observe(size, CountCategory.UNSPLICED)
        im.observe(size, CountCategory.SPLICED_UNANNOT)
        im.observe(size, CountCategory.INTRON)
        im.observe(size, None)  # intergenic
    return im


def _make_strand_models_default():
    """StrandModels with untrained intergenic model (P_strand=0.5)."""
    return StrandModels()


def _make_em_data_with_gdna(
    t_indices_per_unit,
    log_liks_per_unit,
    count_cols_per_unit=None,
    gdna_base_index=None,
):
    """Build EMData where some candidates are gDNA (using gdna_index)."""
    offsets = [0]
    flat_t = []
    flat_lk = []
    flat_cc = []

    for u, t_list in enumerate(t_indices_per_unit):
        for j, t_idx in enumerate(t_list):
            flat_t.append(t_idx)
            flat_lk.append(log_liks_per_unit[u][j])
            if count_cols_per_unit is not None:
                flat_cc.append(count_cols_per_unit[u][j])
            else:
                flat_cc.append(_UNSPLICED_SENSE)
        offsets.append(len(flat_t))

    n_units = len(t_indices_per_unit)
    n_candidates = len(flat_t)

    if gdna_base_index is None:
        gdna_base_index = (max(flat_t) + 1) if flat_t else 0

    # Build locus tracking: best non-gDNA candidate per unit
    locus_t = np.full(n_units, -1, dtype=np.int32)
    locus_cc = np.zeros(n_units, dtype=np.uint8)
    for u, t_list in enumerate(t_indices_per_unit):
        cc_list = count_cols_per_unit[u] if count_cols_per_unit else None
        best_ll = -np.inf
        for j, t_idx in enumerate(t_list):
            if t_idx < gdna_base_index and log_liks_per_unit[u][j] > best_ll:
                best_ll = log_liks_per_unit[u][j]
                locus_t[u] = t_idx
                locus_cc[u] = (
                    cc_list[j] if cc_list else _UNSPLICED_SENSE
                )

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
# _score_gdna_candidate
# =====================================================================


class TestScoreGDNA:
    """Tests for pipeline._score_gdna_candidate."""

    def test_unspliced_no_penalty(self):
        """UNSPLICED gets splice_penalty=1.0 → log(1)=0."""
        im = _make_insert_models()
        sms = _make_strand_models_default()
        cat = int(CountCategory.UNSPLICED)
        score = _score_gdna_candidate(
            int(Strand.POS), cat, 250,
            sms, im,
        )
        isize_model = im.category_models.get(cat, im.global_model)
        assert score < 0
        assert score == pytest.approx(
            np.log(0.5) + isize_model.log_likelihood(250)
        )

    def test_intron_no_penalty(self):
        """INTRON gets splice_penalty=1.0."""
        im = _make_insert_models()
        sms = _make_strand_models_default()
        cat = int(CountCategory.INTRON)
        score = _score_gdna_candidate(
            int(Strand.POS), cat, 250,
            sms, im,
        )
        isize_model = im.category_models.get(cat, im.global_model)
        expected = np.log(0.5) + isize_model.log_likelihood(250)
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
        assert score_unannot < score_unspliced
        penalty_diff = score_unspliced - score_unannot
        assert penalty_diff == pytest.approx(-np.log(0.01), rel=0.01)

    def test_custom_penalties(self):
        """Custom penalty dict overrides defaults."""
        im = _make_insert_models()
        sms = _make_strand_models_default()
        cat = int(CountCategory.UNSPLICED)
        custom = {CountCategory.UNSPLICED: 0.5}
        score = _score_gdna_candidate(
            int(Strand.POS), cat, 250,
            sms, im,
            gdna_splice_penalties=custom,
        )
        isize_model = im.category_models.get(cat, im.global_model)
        expected = np.log(0.5) + isize_model.log_likelihood(250) + np.log(0.5)
        assert score == pytest.approx(expected)

    def test_zero_insert_size(self):
        """insert_size=0 → P_insert=1.0 (log=0)."""
        im = _make_insert_models()
        sms = _make_strand_models_default()
        score = _score_gdna_candidate(
            int(Strand.POS), int(CountCategory.UNSPLICED), 0,
            sms, im,
        )
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

        em = _make_em_data_with_gdna(
            [[0, 3]], [[0.0, 0.0]], gdna_base_index=3,
        )
        rc.run_em(em, em_iterations=1)
        assert rc._converged_theta.shape == (5,)
        assert rc._converged_theta.sum() == pytest.approx(1.0)

    def test_shadow_init_feeds_prior(self):
        """shadow_init biases the per-gene gDNA prior."""
        shadow = np.array([100.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
            shadow_init=shadow,
        )
        em = _make_em_data_with_gdna(
            [[0, 2]] * 100, [[0.0, 0.0]] * 100, gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)

        theta = rc._converged_theta
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
        assert theta[2] < 0.01, f"gDNA theta too high: {theta[2]:.4f}"


# =====================================================================
# gDNA assignment routing
# =====================================================================


class TestGDNAAssignment:
    """Verify that gDNA assignments go to gdna_em_counts, not em_counts."""

    def test_gdna_assignment_not_in_em_counts(self):
        """Fragments assigned to gDNA are NOT in em_counts."""
        shadow_idx = 2
        shadow = np.array([1000.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
            shadow_init=shadow,
        )
        em = _make_em_data_with_gdna(
            [[0, shadow_idx]] * 100,
            [[-10.0, 0.0]] * 100,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        assert rc.gdna_em_count > 0
        total_assigned = rc.em_counts.sum() + rc.gdna_em_count
        assert total_assigned == 100.0

    def test_total_counts_preserved(self):
        """em_counts + gdna_em_count == n_units for all assigned fragments."""
        shadow = np.array([25.0, 25.0])
        rc = ReadCounter(
            num_transcripts=3, num_genes=2, seed=42,
            shadow_init=shadow,
        )
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
        shadow_idx = 2
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
        )
        rc.unique_counts[0, 0] = 500.0

        em = _make_em_data_with_gdna(
            [[0, shadow_idx]] * 500,
            [[0.0, -20.0]] * 500,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        assert rc.em_counts[0].sum() > 490
        assert rc.gdna_em_count < 10

    def test_high_shadow_prior_absorbs_fragments(self):
        """With large shadow_init, gDNA takes larger share."""
        shadow_idx = 2
        shadow = np.array([500.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
            shadow_init=shadow,
        )
        em = _make_em_data_with_gdna(
            [[0, shadow_idx]] * 500,
            [[0.0, 0.0]] * 500,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

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
        shadow_idx = 2
        shadow = np.array([1000.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
            shadow_init=shadow,
        )
        cc = _UNSPLICED_SENSE
        em = _make_em_data_with_gdna(
            [[0, shadow_idx]] * 100,
            [[-10.0, 0.0]] * 100,
            count_cols_per_unit=[[cc, cc]] * 100,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        assert rc.gdna_locus_counts[0].sum() > 0

    def test_locus_counts_zero_when_no_gdna(self):
        """No gDNA assignments → locus counts stay zero."""
        shadow_idx = 2
        rc = ReadCounter(
            num_transcripts=2, num_genes=1, seed=42,
        )
        rc.unique_counts[0, 0] = 1000.0
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
        assert rc.gdna_contamination_rate == pytest.approx(0.1)

    def test_contamination_rate_zero_when_empty(self):
        rc = ReadCounter(num_transcripts=2, num_genes=1, seed=42)
        assert rc.gdna_contamination_rate == 0.0


# =====================================================================
# gDNA output methods
# =====================================================================


class TestGDNAOutput:
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

        unique_counts = np.zeros((3, NUM_COUNT_COLS), dtype=np.float64)
        t_to_g = np.array([0, 0, 1], dtype=np.int64)
        result = compute_shadow_init(unique_counts, t_to_g, 2, 0)
        np.testing.assert_array_equal(result, [0.0, 0.0])

    def test_no_genic_reads_gives_zeros(self):
        """With only intergenic fragments (no genic reads), shadow = 0."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((2, NUM_COUNT_COLS), dtype=np.float64)
        t_to_g = np.array([0, 1], dtype=np.int64)
        result = compute_shadow_init(unique_counts, t_to_g, 2, 100)
        np.testing.assert_array_equal(result, [0.0, 0.0])

    def test_intergenic_with_gene_depth_sets_floor(self):
        """Genes with reads get a depth-scaled floor even without antisense."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((2, NUM_COUNT_COLS), dtype=np.float64)
        t_to_g = np.array([0, 1], dtype=np.int64)
        # INTRON_SENSE column
        intron_sense_col = int(CountCol.INTRON_SENSE)
        unique_counts[0, intron_sense_col] = 100.0
        unique_counts[1, intron_sense_col] = 100.0
        result = compute_shadow_init(unique_counts, t_to_g, 2, 100, kappa=1.0)
        # n_total = 200 + 100 = 300
        # estimated_gdna = 100 + 0 = 100
        # θ_global = 100 / 300 = 1/3
        assert result[0] == pytest.approx(100.0 / 3)
        assert result[1] == pytest.approx(100.0 / 3)

    def test_antisense_boosts_gene_shadow(self):
        """Genes with antisense counts get proportionally higher shadow."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((3, NUM_COUNT_COLS), dtype=np.float64)
        t_to_g = np.array([0, 0, 1], dtype=np.int64)
        # Gene 0: 50 antisense via t0 (e.g. INTRON_ANTISENSE = col 1)
        intron_anti_col = int(CountCol.INTRON_ANTISENSE)
        unique_counts[0, intron_anti_col] = 50.0
        # Gene 1: 50 sense reads via t2
        intron_sense_col = int(CountCol.INTRON_SENSE)
        unique_counts[2, intron_sense_col] = 50.0
        result = compute_shadow_init(unique_counts, t_to_g, 2, 0, kappa=1.0)
        # n_total = 100, estimated_gdna = 0 + 2*50 = 100
        # θ_global = 100 / 100 = 1.0
        # gene 0: 2*50 + 1.0*50*1.0 = 150.0
        # gene 1: 0 + 1.0*50*1.0 = 50.0
        assert result[0] == pytest.approx(150.0)
        assert result[1] == pytest.approx(50.0)
        assert result[0] > result[1]

    def test_kappa_scales_global_shrinkage(self):
        """Higher kappa → more shrinkage toward global rate."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((2, NUM_COUNT_COLS), dtype=np.float64)
        t_to_g = np.array([0, 1], dtype=np.int64)
        intron_sense = int(CountCol.INTRON_SENSE)
        unique_counts[0, intron_sense] = 100.0
        unique_counts[1, intron_sense] = 100.0
        r_lo = compute_shadow_init(unique_counts, t_to_g, 2, 100, kappa=0.5)
        r_hi = compute_shadow_init(unique_counts, t_to_g, 2, 100, kappa=5.0)
        assert r_hi[0] > r_lo[0]

    def test_depth_scaling(self):
        """Deeper gene gets larger absolute shadow floor."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((2, NUM_COUNT_COLS), dtype=np.float64)
        t_to_g = np.array([0, 1], dtype=np.int64)
        intron_sense = int(CountCol.INTRON_SENSE)
        unique_counts[0, intron_sense] = 1000.0
        unique_counts[1, intron_sense] = 10.0
        result = compute_shadow_init(unique_counts, t_to_g, 2, 100, kappa=1.0)
        assert result[0] > result[1] * 50
        assert result[0] > 0
        assert result[1] > 0

    def test_output_shape(self):
        """Returns array of shape (num_genes,)."""
        from hulkrna.pipeline import compute_shadow_init

        unique_counts = np.zeros((5, NUM_COUNT_COLS), dtype=np.float64)
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
        shadow = np.array([100.0, 0.0])
        rc = ReadCounter(
            num_transcripts=3, num_genes=2, seed=42,
            shadow_init=shadow,
        )
        em = _make_em_data_with_gdna(
            [[0, 3]] * 200,
            [[0.0, 0.0]] * 200,
            gdna_base_index=3,
        )
        rc.run_em(em, em_iterations=10)

        theta = rc._converged_theta
        assert theta.shape == (5,)
        assert theta[3] > 0.1
        assert theta[4] < theta[3]

    def test_per_gene_em_counts_tracking(self):
        """Verify gdna_em_counts tracks per-gene gDNA assignments."""
        shadow = np.array([500.0, 500.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=2, seed=42,
            shadow_init=shadow,
        )
        em = _make_em_data_with_gdna(
            [[0, 2]] * 200 + [[1, 3]] * 200,
            [[0.0, 0.0]] * 200 + [[0.0, 0.0]] * 200,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

        assert rc.gdna_em_counts[0] > 0
        assert rc.gdna_em_counts[1] > 0
        total = rc.em_counts.sum() + rc.gdna_em_counts.sum()
        assert total == 400.0

    def test_gene_with_high_antisense_absorbs_more(self):
        """Gene with high shadow_init absorbs more ambiguous fragments."""
        shadow = np.array([200.0, 1.0])
        rc = ReadCounter(
            num_transcripts=2, num_genes=2, seed=42,
            shadow_init=shadow,
        )
        rc.unique_counts[0, 0] = 10.0
        rc.unique_counts[1, 0] = 10.0

        em = _make_em_data_with_gdna(
            [[0, 2]] * 500 + [[1, 3]] * 500,
            [[0.0, 0.0]] * 500 + [[0.0, 0.0]] * 500,
            gdna_base_index=2,
        )
        rc.run_em(em, em_iterations=10)
        rc.assign_ambiguous(em)

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
        """model_for_category returns the best pure-RNA model."""
        from hulkrna.strand_model import StrandModels, StrandModel
        from hulkrna.categories import CountCategory
        from hulkrna.types import Strand

        # --- Case 1: exonic_spliced has enough observations ---
        sm1 = StrandModels()
        for _ in range(20):
            sm1.exonic_spliced.observe(Strand.POS, Strand.POS)
        m = sm1.model_for_category(int(CountCategory.SPLICED_ANNOT))
        assert m is sm1.exonic_spliced
        assert sm1.model_for_category(int(CountCategory.UNSPLICED)) is m
        assert sm1.model_for_category(int(CountCategory.INTRON)) is m

        # --- Case 2: exonic_spliced insufficient, exonic sufficient ---
        sm2 = StrandModels()
        for _ in range(20):
            sm2.exonic.observe(Strand.POS, Strand.POS)
            sm2.intergenic.observe(Strand.POS, Strand.POS)
        m2 = sm2.model_for_category(int(CountCategory.UNSPLICED))
        assert m2 is not sm2.exonic_spliced
        assert sm2.model_for_category(int(CountCategory.INTRON)) is m2

        # --- Case 3: default (no observations) → falls back to exonic ---
        sm3 = StrandModels()
        m3 = sm3.model_for_category(int(CountCategory.UNSPLICED))
        assert m3 is sm3.exonic
