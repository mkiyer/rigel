"""Tests for rigel.strand_model — strand model."""

import json

import pytest

from rigel.types import Strand
from rigel.strand_model import StrandModel, StrandModels


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
    """Test MLE probability computation."""

    def test_no_observations(self):
        sm = StrandModel()
        assert sm.p_r1_sense == 0.5

    def test_strong_fr_library(self):
        """Simulate a strongly R1-sense library (most same-direction)."""
        sm = StrandModel()
        for _ in range(95):
            sm.observe(Strand.POS, Strand.POS)
        for _ in range(5):
            sm.observe(Strand.POS, Strand.NEG)
        # p_r1_sense = 95 / 100 = 0.95
        assert sm.p_r1_sense == pytest.approx(0.95)
        assert sm.strand_specificity > 0.9
        assert sm.read1_sense is True

    def test_strong_rf_library(self):
        """Simulate a strongly R1-antisense library (most opposite-direction)."""
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
        sm.finalize()
        p = sm.strand_likelihood(Strand.POS, Strand.POS)
        assert p == pytest.approx(sm.p_r1_sense)

    def test_opposite_direction_returns_p_r1_antisense(self):
        sm = StrandModel()
        for _ in range(50):
            sm.observe(Strand.POS, Strand.POS)
        sm.finalize()
        p = sm.strand_likelihood(Strand.POS, Strand.NEG)
        assert p == pytest.approx(sm.p_r1_antisense)

    def test_ambiguous_exon_returns_half(self):
        sm = StrandModel()
        sm.observe(Strand.POS, Strand.POS)
        sm.finalize()
        p = sm.strand_likelihood(Strand.AMBIGUOUS, Strand.POS)
        assert p == 0.5

    def test_none_gene_strand_returns_half(self):
        sm = StrandModel()
        sm.observe(Strand.POS, Strand.POS)
        sm.finalize()
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
        assert "estimate" in d
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
        assert "estimate" in data["strand_model"]

    def test_write_json_includes_ci_with_enough_observations(self, tmp_path):
        sm = StrandModel()
        for _ in range(50):
            sm.observe(Strand.POS, Strand.POS)
        for _ in range(5):
            sm.observe(Strand.POS, Strand.NEG)
        path = tmp_path / "strand.json"
        sm.write_json(path)
        data = json.loads(path.read_text())
        ci = data["strand_model"]["estimate"].get("ci_95")
        assert ci is not None
        assert len(ci) == 2
        assert ci[0] < ci[1]


class TestStrandModelProperties:
    def test_posterior_variance_with_data(self):
        sm = StrandModel()
        for _ in range(80):
            sm.observe(Strand.POS, Strand.POS)
        for _ in range(20):
            sm.observe(Strand.POS, Strand.NEG)
        # p = 80/100 = 0.8, variance = 0.8 * 0.2 / 100 = 0.0016
        assert sm.posterior_variance() == pytest.approx(0.0016)

    def test_posterior_variance_no_observations(self):
        sm = StrandModel()
        assert sm.posterior_variance() == 0.25

    def test_posterior_95ci(self):
        sm = StrandModel()
        for _ in range(95):
            sm.observe(Strand.POS, Strand.POS)
        for _ in range(5):
            sm.observe(Strand.POS, Strand.NEG)
        lo, hi = sm.posterior_95ci()
        assert lo < hi
        assert lo > 0.88  # should be heavily skewed towards 1.0
        assert hi <= 1.0


class TestStrandModelsContainer:
    """Test the simplified StrandModels container."""

    def test_delegation_to_exonic_spliced(self):
        models = StrandModels()
        for _ in range(50):
            models.exonic_spliced.observe(Strand.POS, Strand.POS)
        models.finalize()
        assert models.strand_specificity == models.exonic_spliced.strand_specificity
        assert models.p_r1_sense == models.exonic_spliced.p_r1_sense
        assert models.read1_sense == models.exonic_spliced.read1_sense
        assert models.n_observations == models.exonic_spliced.n_observations

    def test_mle_on_finalize(self):
        models = StrandModels()
        for _ in range(10):
            models.exonic_spliced.observe(Strand.POS, Strand.POS)
        models.finalize()
        # MLE: p_r1_sense = 10/10 = 1.0
        assert models.p_r1_sense == pytest.approx(1.0)

    def test_zero_observations_warns(self, caplog):
        models = StrandModels()
        with caplog.at_level("WARNING"):
            models.finalize()
        assert "No spliced strand observations" in caplog.text
        assert models.p_r1_sense == 0.5

    def test_low_observations_warns(self, caplog):
        models = StrandModels()
        for _ in range(5):
            models.exonic_spliced.observe(Strand.POS, Strand.POS)
        with caplog.at_level("WARNING"):
            models.finalize()
        assert "Only 5 spliced strand observations" in caplog.text

    def test_diagnostic_models_finalized(self):
        models = StrandModels()
        for _ in range(10):
            models.exonic.observe(Strand.POS, Strand.POS)
        for _ in range(10):
            models.intergenic.observe(Strand.POS, Strand.POS)
        models.finalize()
        assert models.exonic._finalized
        assert models.intergenic._finalized

    def test_to_dict_structure(self):
        models = StrandModels()
        for _ in range(20):
            models.exonic_spliced.observe(Strand.POS, Strand.POS)
        models.finalize()
        d = models.to_dict()
        assert "exonic_spliced" in d
        assert "diagnostics" in d
        assert "exonic" in d["diagnostics"]
        assert "intergenic" in d["diagnostics"]
