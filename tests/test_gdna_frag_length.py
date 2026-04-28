"""Tests for gDNA/RNA fragment length distribution framework.

Covers:
- FragmentLengthModel.from_counts() factory for pre-built distributions
- build_scoring_models() correctly sets rna_model = SPLICED_ANNOT
- Integration with scoring (gdna_model vs rna_model likelihood separation)
- Serialization of model fields
"""

import math

import numpy as np
import pytest

from rigel.frag_length_model import (
    FragmentLengthModel,
    FragmentLengthModels,
)
from rigel.splice import SpliceType


# =====================================================================
# FragmentLengthModel.from_counts() factory
# =====================================================================


class TestFromCounts:
    """Tests for the from_counts() factory method."""

    def test_basic_creation(self):
        """Create a model from a simple histogram."""
        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        counts[250] = 50.0
        model = FragmentLengthModel.from_counts(counts)

        assert model._finalized
        assert model.max_size == 500
        assert model.total_weight == pytest.approx(150.0)
        assert model.n_observations == 150
        assert model._log_prob is not None

    def test_max_size_inferred(self):
        """max_size is inferred as len(counts) - 1."""
        counts = np.ones(101)
        model = FragmentLengthModel.from_counts(counts)
        assert model.max_size == 100

    def test_explicit_max_size(self):
        """Explicit max_size with shorter counts array zero-fills."""
        counts = np.array([0, 0, 10, 20, 10])
        model = FragmentLengthModel.from_counts(counts, max_size=1000)
        assert model.max_size == 1000
        assert model.total_weight == pytest.approx(40.0)
        assert model.counts[2] == pytest.approx(10.0)
        assert model.counts[999] == pytest.approx(0.0)

    def test_log_likelihood_matches_trained(self):
        """from_counts() produces same log_likelihood as observe+finalize."""
        trained = FragmentLengthModel(max_size=500)
        for _ in range(100):
            trained.observe(200)
        for _ in range(50):
            trained.observe(300)
        trained.finalize()

        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        counts[300] = 50.0
        factory = FragmentLengthModel.from_counts(counts)

        for length in [0, 100, 200, 250, 300, 500]:
            assert factory.log_likelihood(length) == pytest.approx(trained.log_likelihood(length))

    def test_tail_decay_works(self):
        """Queries beyond max_size use exponential tail decay."""
        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        model = FragmentLengthModel.from_counts(counts)

        ll_500 = model.log_likelihood(500)
        ll_600 = model.log_likelihood(600)
        assert ll_600 < ll_500
        assert ll_600 == pytest.approx(model._tail_base + (600 - 500) * math.log(0.99))

    def test_empty_counts(self):
        """Zero-count histogram yields uniform distribution."""
        counts = np.zeros(101, dtype=np.float64)
        model = FragmentLengthModel.from_counts(counts)
        assert model._finalized
        ll_0 = model.log_likelihood(0)
        ll_50 = model.log_likelihood(50)
        assert ll_0 == pytest.approx(ll_50)
        assert ll_0 == pytest.approx(-np.log(101))

    def test_single_peak(self):
        """Single-peak histogram: mode has highest likelihood."""
        counts = np.zeros(501, dtype=np.float64)
        counts[250] = 1000.0
        model = FragmentLengthModel.from_counts(counts)
        assert model.log_likelihood(250) > model.log_likelihood(100)
        assert model.log_likelihood(250) > model.log_likelihood(400)
        assert model.mode == 250

    def test_statistics(self):
        """Mean, std, median, mode are computed correctly."""
        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        counts[300] = 100.0
        model = FragmentLengthModel.from_counts(counts)
        assert model.mean == pytest.approx(250.0)
        assert model.median == pytest.approx(200.0)
        assert model.mode in (200, 300)

    def test_from_counts_matches_normal_training(self):
        """Manually constructed histogram via from_counts matches
        the same histogram built through observe() calls."""
        trained = FragmentLengthModel(max_size=500)
        for _ in range(100):
            trained.observe(180)
        for _ in range(200):
            trained.observe(220)
        for _ in range(50):
            trained.observe(260)
        trained.finalize()

        counts = np.zeros(501, dtype=np.float64)
        counts[180] = 100
        counts[220] = 200
        counts[260] = 50
        factory = FragmentLengthModel.from_counts(counts)

        for length in range(501):
            assert factory.log_likelihood(length) == pytest.approx(
                trained.log_likelihood(length)
            ), f"Mismatch at length={length}"


# =====================================================================
# build_scoring_models()
# =====================================================================


class TestBuildScoringModels:
    """Verify build_scoring_models sets rna_model = SPLICED_ANNOT."""

    def test_rna_model_equals_spliced_annot(self):
        """rna_model histogram matches SPLICED_ANNOT after build."""
        m = FragmentLengthModels(max_size=500)
        for fl in [180, 200, 220, 200, 210]:
            m.observe(fl, splice_type=SpliceType.SPLICED_ANNOT)
        for fl in [300, 310, 320]:
            m.observe(fl, splice_type=SpliceType.UNSPLICED)
        for fl in [400, 410]:
            m.observe(fl, splice_type=None)  # intergenic

        m.build_scoring_models()

        spliced = m.category_models[SpliceType.SPLICED_ANNOT]
        np.testing.assert_array_equal(m.rna_model.counts, spliced.counts)
        assert m.rna_model.total_weight == pytest.approx(spliced.total_weight)

    def test_rna_model_independent_of_intergenic(self):
        """Intergenic observations do NOT affect rna_model."""
        m = FragmentLengthModels(max_size=500)
        for fl in [200] * 100:
            m.observe(fl, splice_type=SpliceType.SPLICED_ANNOT)
        for fl in [400] * 500:
            m.observe(fl, splice_type=None)

        m.build_scoring_models()
        assert m.rna_model.total_weight == pytest.approx(100.0)
        assert m.rna_model.mean == pytest.approx(200.0)

    def test_rna_model_independent_of_unspliced(self):
        """Unspliced genic observations do NOT affect rna_model."""
        m = FragmentLengthModels(max_size=500)
        for fl in [200] * 100:
            m.observe(fl, splice_type=SpliceType.SPLICED_ANNOT)
        for fl in [350] * 200:
            m.observe(fl, splice_type=SpliceType.UNSPLICED)

        m.build_scoring_models()
        assert m.rna_model.total_weight == pytest.approx(100.0)
        assert m.rna_model.mean == pytest.approx(200.0)

    def test_empty_spliced_uses_global_prior(self):
        """No spliced observations -> rna_model has 0 weight but gets global prior."""
        m = FragmentLengthModels(max_size=100)
        for fl in [50] * 10:
            m.observe(fl, splice_type=SpliceType.UNSPLICED)
        m.build_scoring_models()
        # rna_model stays at 0 category-specific weight;
        # global prior (applied during finalize) provides the distribution.
        assert m.rna_model.total_weight == pytest.approx(0.0)
        m.finalize()
        assert m.rna_model._finalized
        assert np.isfinite(m.rna_model.log_likelihood(50))

    def test_finalize_after_build(self):
        """build -> finalize produces valid lookup tables."""
        m = FragmentLengthModels(max_size=500)
        for fl in [200] * 200:
            m.observe(fl, splice_type=SpliceType.SPLICED_ANNOT)
        for fl in [350] * 100:
            m.observe(fl, splice_type=None)

        m.build_scoring_models()
        m.finalize()

        assert m.rna_model._finalized
        assert m.rna_model._log_prob is not None
        assert np.isfinite(m.rna_model.log_likelihood(200))


# =====================================================================
# Scoring integration (gdna injected externally)
# =====================================================================


class TestScoringIntegration:
    """Integration: build -> inject gdna -> finalize -> scoring."""

    def test_rna_gdna_discriminate(self):
        """RNA and gDNA models assign higher probability to their own peak."""
        m = FragmentLengthModels(max_size=500)
        for fl in [150] * 300:
            m.observe(fl, splice_type=SpliceType.SPLICED_ANNOT)
        for fl in [400] * 200:
            m.observe(fl, splice_type=None)

        m.build_scoring_models()

        gdna_counts = np.zeros(501, dtype=np.float64)
        gdna_counts[400] = 500.0
        m.gdna_model = FragmentLengthModel.from_counts(gdna_counts)

        m.finalize()

        assert m.rna_model.log_likelihood(150) > m.gdna_model.log_likelihood(150)
        assert m.gdna_model.log_likelihood(400) > m.rna_model.log_likelihood(400)

    def test_custom_gdna_model_override(self):
        """Install a from_counts() model as gdna_model on a container."""
        m = FragmentLengthModels(max_size=500)
        for _ in range(200):
            m.observe(200, SpliceType.SPLICED_ANNOT)
            m.observe(200, None)

        m.build_scoring_models()

        gdna_counts = np.zeros(501, dtype=np.float64)
        gdna_counts[350] = 500.0
        gdna_counts[400] = 300.0
        m.gdna_model = FragmentLengthModel.from_counts(gdna_counts)
        m.finalize()

        assert m.gdna_model.log_likelihood(350) > m.gdna_model.log_likelihood(200)
        assert m.rna_model.log_likelihood(200) > m.rna_model.log_likelihood(350)

    def test_log_likelihood_is_proper_distribution(self):
        """Exponentiated log-likelihoods approximately sum to 1."""
        m = FragmentLengthModels(max_size=500)
        for fl in [200] * 500:
            m.observe(fl, splice_type=SpliceType.SPLICED_ANNOT)
        for fl in [350] * 300:
            m.observe(fl, splice_type=None)

        m.build_scoring_models()

        gdna_counts = np.zeros(501, dtype=np.float64)
        gdna_counts[350] = 300.0
        m.gdna_model = FragmentLengthModel.from_counts(gdna_counts)
        m.finalize()

        for model in [m.rna_model, m.gdna_model]:
            probs = np.exp([model.log_likelihood(k) for k in range(501)])
            assert abs(probs.sum() - 1.0) < 0.01


# =====================================================================
# Serialization
# =====================================================================


class TestToDict:
    """Verify serialization includes expected keys."""

    def test_to_dict_has_all_keys(self):
        m = FragmentLengthModels(max_size=500)
        for fl in [200] * 50:
            m.observe(fl, splice_type=SpliceType.SPLICED_ANNOT)
        for fl in [300] * 20:
            m.observe(fl, splice_type=None)

        m.build_scoring_models()
        m.finalize()

        d = m.to_dict()
        assert "rna" in d
        assert "gdna" in d
        assert "global" in d
        assert "unspliced_same_strand" not in d
        assert "unspliced_opp_strand" not in d

    def test_to_dict_summary_accurate(self):
        """Serialized summary stats match live model."""
        m = FragmentLengthModels(max_size=500)
        for fl in [200] * 100:
            m.observe(fl, splice_type=SpliceType.SPLICED_ANNOT)
        for fl in [350] * 50:
            m.observe(fl, splice_type=None)

        m.build_scoring_models()
        m.finalize()

        d = m.to_dict()
        rna_summary = d["rna"]["summary"]
        assert rna_summary["total_weight"] == pytest.approx(m.rna_model.total_weight, rel=0.01)
