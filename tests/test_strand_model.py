"""Tests for hulkrna.strand_model — Bayesian strand model."""

import json
import math

import pytest

from hulkrna.types import Strand
from hulkrna.strand_model import StrandModel, StrandModels


class TestStrandModelObserve:
    """Test the 2×2 contingency table and derived counts."""

    def test_default_counts_are_zero(self):
        sm = StrandModel()
        assert sm.pos_pos == 0
        assert sm.pos_neg == 0
        assert sm.neg_pos == 0
        assert sm.neg_neg == 0
        assert sm.n_observations == 0

    def test_observe_pos_pos(self):
        sm = StrandModel()
        sm.observe(Strand.POS, Strand.POS)
        assert sm.pos_pos == 1
        assert sm.n_same == 1
        assert sm.n_opposite == 0

    def test_observe_pos_neg(self):
        sm = StrandModel()
        sm.observe(Strand.POS, Strand.NEG)
        assert sm.pos_neg == 1
        assert sm.n_same == 0
        assert sm.n_opposite == 1

    def test_observe_neg_pos(self):
        sm = StrandModel()
        sm.observe(Strand.NEG, Strand.POS)
        assert sm.neg_pos == 1
        assert sm.n_opposite == 1

    def test_observe_neg_neg(self):
        sm = StrandModel()
        sm.observe(Strand.NEG, Strand.NEG)
        assert sm.neg_neg == 1
        assert sm.n_same == 1

    def test_n_observations(self):
        sm = StrandModel()
        sm.observe(Strand.POS, Strand.POS)
        sm.observe(Strand.POS, Strand.NEG)
        sm.observe(Strand.NEG, Strand.POS)
        assert sm.n_observations == 3


class TestStrandModelPosterior:
    """Test Beta posterior computation."""

    def test_prior_only(self):
        sm = StrandModel(prior_alpha=1.0, prior_beta=1.0)
        assert sm.alpha == 1.0
        assert sm.beta == 1.0
        assert sm.p_r1_sense == 0.5

    def test_strong_fr_library(self):
        """Simulate a strongly FR-stranded library (most same-direction)."""
        sm = StrandModel()
        for _ in range(95):
            sm.observe(Strand.POS, Strand.POS)
        for _ in range(5):
            sm.observe(Strand.POS, Strand.NEG)
        # p_r1_sense ≈ (1 + 95) / (1 + 95 + 1 + 5) = 96/102 ≈ 0.941
        assert sm.p_r1_sense > 0.9
        assert sm.strand_specificity > 0.9
        assert sm.read1_sense is True

    def test_strong_rf_library(self):
        """Simulate a strongly RF-stranded library (most opposite-direction)."""
        sm = StrandModel()
        for _ in range(5):
            sm.observe(Strand.POS, Strand.POS)
        for _ in range(95):
            sm.observe(Strand.POS, Strand.NEG)
        assert sm.p_r1_sense < 0.1
        assert sm.p_r1_antisense > 0.9
        assert sm.strand_specificity > 0.9
        assert sm.read1_sense is False


class TestStrandModelLikelihood:
    """Test the strand_likelihood() method."""

    def test_same_direction_returns_p_r1_sense(self):
        sm = StrandModel()
        for _ in range(50):
            sm.observe(Strand.POS, Strand.POS)
        p = sm.strand_likelihood(Strand.POS, Strand.POS)
        assert p == pytest.approx(sm.p_r1_sense)

    def test_opposite_direction_returns_p_r1_antisense(self):
        sm = StrandModel()
        for _ in range(50):
            sm.observe(Strand.POS, Strand.POS)
        p = sm.strand_likelihood(Strand.POS, Strand.NEG)
        assert p == pytest.approx(sm.p_r1_antisense)

    def test_ambiguous_exon_returns_half(self):
        sm = StrandModel()
        sm.observe(Strand.POS, Strand.POS)
        p = sm.strand_likelihood(Strand.AMBIGUOUS, Strand.POS)
        assert p == 0.5

    def test_none_gene_strand_returns_half(self):
        sm = StrandModel()
        sm.observe(Strand.POS, Strand.POS)
        p = sm.strand_likelihood(Strand.POS, Strand.NONE)
        assert p == 0.5


class TestStrandModelSerialization:
    """Test to_dict and write_json."""

    def test_to_dict_structure(self):
        sm = StrandModel()
        sm.observe(Strand.POS, Strand.POS)
        sm.observe(Strand.NEG, Strand.POS)
        d = sm.to_dict()
        assert "observations" in d
        assert "posterior" in d
        assert "probabilities" in d
        assert "protocol" in d
        assert d["observations"]["total"] == 2
        assert d["observations"]["n_same"] == 1
        assert d["observations"]["n_opposite"] == 1

    def test_write_json_creates_file(self, tmp_path):
        sm = StrandModel()
        for _ in range(20):
            sm.observe(Strand.POS, Strand.POS)
        path = tmp_path / "strand.json"
        sm.write_json(path)
        assert path.exists()
        data = json.loads(path.read_text())
        assert "strand_model" in data
        assert "posterior" in data["strand_model"]

    def test_write_json_includes_ci_with_enough_observations(self, tmp_path):
        sm = StrandModel()
        for _ in range(50):
            sm.observe(Strand.POS, Strand.POS)
        for _ in range(5):
            sm.observe(Strand.POS, Strand.NEG)
        path = tmp_path / "strand.json"
        sm.write_json(path)
        data = json.loads(path.read_text())
        ci = data["strand_model"]["posterior"].get("ci_95")
        assert ci is not None
        assert len(ci) == 2
        assert ci[0] < ci[1]


class TestStrandModelProperties:
    def test_posterior_variance(self):
        sm = StrandModel(prior_alpha=2.0, prior_beta=3.0)
        # No observations: alpha=2, beta=3
        # variance = (2*3) / ((5)^2 * 6) = 6/150 = 0.04
        assert sm.posterior_variance() == pytest.approx(0.04)

    def test_posterior_95ci(self):
        sm = StrandModel()
        for _ in range(100):
            sm.observe(Strand.POS, Strand.POS)
        lo, hi = sm.posterior_95ci()
        assert lo < hi
        assert lo > 0.9  # should be heavily skewed towards 1.0
        assert hi <= 1.0


class TestStrandModelsFallback:
    def test_min_spliced_observations_threshold_controls_primary_model(self):
        models = StrandModels()
        models.min_spliced_observations = 5

        for _ in range(5):
            models.exonic_spliced.observe(Strand.POS, Strand.POS)
        for _ in range(40):
            models.exonic.observe(Strand.POS, Strand.NEG)

        best = models.model_for_category(0)
        assert best is models.exonic_spliced

    def test_pooled_fallback_includes_spliced_and_warns(self, caplog):
        models = StrandModels()
        models.min_spliced_observations = 10

        for _ in range(3):
            models.exonic_spliced.observe(Strand.POS, Strand.POS)
        for _ in range(20):
            models.exonic.observe(Strand.POS, Strand.POS)

        with caplog.at_level("WARNING"):
            best = models.model_for_category(0)

        assert best is not models.exonic
        assert best is not models.exonic_spliced
        assert best.n_observations == (
            models.exonic.n_observations + models.exonic_spliced.n_observations
        )
        assert "using pooled fallback model" in caplog.text

    def test_raises_when_pooled_observations_below_threshold(self):
        models = StrandModels()
        models.min_spliced_observations = 10

        for _ in range(2):
            models.exonic_spliced.observe(Strand.POS, Strand.POS)
        for _ in range(3):
            models.exonic.observe(Strand.POS, Strand.POS)

        with pytest.raises(RuntimeError, match="Insufficient strand"):
            models.model_for_category(0)
