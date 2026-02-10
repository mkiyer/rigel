"""Tests for hulkrna.insert_model — Insert size distribution model."""

import numpy as np
import pytest

from hulkrna.insert_model import InsertSizeModel


class TestInsertSizeModelBasic:
    """Test basic InsertSizeModel functionality."""

    def test_default_construction(self):
        model = InsertSizeModel()
        assert model.max_size == 1000
        assert model.n_observations == 0
        assert model.counts.shape == (1001,)
        assert model.total_weight == 0.0

    def test_custom_max_size(self):
        model = InsertSizeModel(max_size=500)
        assert model.max_size == 500
        assert model.counts.shape == (501,)

    def test_observe_single(self):
        model = InsertSizeModel()
        model.observe(200)
        assert model.n_observations == 1
        assert model.counts[200] == 1.0
        assert model.total_weight == 1.0

    def test_observe_multiple(self):
        model = InsertSizeModel()
        model.observe(200)
        model.observe(200)
        model.observe(300)
        assert model.n_observations == 3
        assert model.counts[200] == 2.0
        assert model.counts[300] == 1.0

    def test_observe_weighted(self):
        model = InsertSizeModel()
        model.observe(200, weight=2.5)
        assert model.n_observations == 1
        assert model.counts[200] == 2.5
        assert model.total_weight == 2.5

    def test_observe_clamped_to_max(self):
        model = InsertSizeModel(max_size=500)
        model.observe(600)
        assert model.counts[500] == 1.0  # clamped to overflow bin
        assert model.counts[499] == 0.0

    def test_observe_clamped_negative(self):
        model = InsertSizeModel()
        model.observe(-5)
        assert model.counts[0] == 1.0  # clamped to 0

    def test_observe_at_max_size(self):
        model = InsertSizeModel(max_size=100)
        model.observe(100)
        assert model.counts[100] == 1.0


class TestInsertSizeModelStatistics:
    """Test summary statistics."""

    def test_mean_simple(self):
        model = InsertSizeModel()
        model.observe(200)
        model.observe(300)
        assert model.mean == 250.0

    def test_mean_empty(self):
        model = InsertSizeModel()
        assert model.mean == 0.0

    def test_std_simple(self):
        model = InsertSizeModel()
        model.observe(200)
        model.observe(300)
        assert model.std == pytest.approx(50.0, abs=0.01)

    def test_std_empty(self):
        model = InsertSizeModel()
        assert model.std == 0.0

    def test_median_simple(self):
        model = InsertSizeModel()
        for _ in range(10):
            model.observe(200)
        for _ in range(5):
            model.observe(300)
        # Median should be 200 (more weight there)
        assert model.median == 200.0

    def test_median_empty(self):
        model = InsertSizeModel()
        assert model.median == 0.0

    def test_mode(self):
        model = InsertSizeModel()
        model.observe(200)
        model.observe(200)
        model.observe(300)
        assert model.mode == 200

    def test_mode_empty(self):
        model = InsertSizeModel()
        assert model.mode == 0


class TestInsertSizeModelLikelihood:
    """Test log-likelihood computation."""

    def test_log_likelihood_uniform_when_empty(self):
        model = InsertSizeModel(max_size=100)
        ll = model.log_likelihood(50)
        assert ll == pytest.approx(-np.log(101), abs=1e-10)

    def test_log_likelihood_observed_size(self):
        model = InsertSizeModel(max_size=100)
        for _ in range(100):
            model.observe(200)  # clamped to 100
        ll_at_100 = model.log_likelihood(100)
        ll_at_50 = model.log_likelihood(50)
        # Size 100 (overflow) observed 100 times, size 50 never
        assert ll_at_100 > ll_at_50

    def test_log_likelihood_never_negative_inf(self):
        """Laplace smoothing ensures no -inf."""
        model = InsertSizeModel(max_size=100)
        model.observe(50)
        ll = model.log_likelihood(99)  # never observed
        assert np.isfinite(ll)
        assert ll < 0


class TestInsertSizeModelSerialization:
    """Test to_dict and write_json."""

    def test_to_dict_empty(self):
        model = InsertSizeModel()
        d = model.to_dict()
        assert d["summary"]["n_observations"] == 0
        assert d["summary"]["mean"] == 0.0
        assert d["histogram"]["values"] == []

    def test_to_dict_with_data(self):
        model = InsertSizeModel()
        model.observe(200)
        model.observe(200)
        model.observe(300)
        d = model.to_dict()
        assert d["summary"]["n_observations"] == 3
        assert d["summary"]["mode"] == 200
        assert d["histogram"]["range"][0] == 200
        assert d["histogram"]["range"][1] == 300
        # Values should cover bins 200..300
        assert len(d["histogram"]["values"]) == 101
        assert d["histogram"]["values"][0] == 2.0   # bin 200
        assert d["histogram"]["values"][-1] == 1.0  # bin 300

    def test_write_json(self, tmp_path):
        import json
        model = InsertSizeModel()
        model.observe(250)
        model.observe(250)
        model.observe(300)
        json_path = tmp_path / "insert_model.json"
        model.write_json(json_path)
        assert json_path.exists()
        with open(json_path) as f:
            data = json.load(f)
        assert "insert_size_model" in data
        assert data["insert_size_model"]["summary"]["n_observations"] == 3

    def test_to_dict_native_types(self):
        """All values should be native Python types (not numpy)."""
        model = InsertSizeModel()
        model.observe(200)
        d = model.to_dict()
        # Check that everything is serializable
        import json
        json.dumps(d)  # should not raise


# ---------------------------------------------------------------------------
# Tests: InsertSizeModels (per-category container)
# ---------------------------------------------------------------------------

class TestInsertSizeModels:
    """Test the InsertSizeModels container class."""

    def test_default_construction(self):
        from hulkrna.insert_model import InsertSizeModels
        models = InsertSizeModels()
        assert models.max_size == 1000
        assert models.n_observations == 0
        assert models.global_model.n_observations == 0
        assert models.intergenic.n_observations == 0

    def test_observe_with_category(self):
        from hulkrna.insert_model import InsertSizeModels
        from hulkrna.core import CountCategory

        models = InsertSizeModels()
        models.observe(200, count_cat=CountCategory.SPLICED_ANNOT)
        models.observe(300, count_cat=CountCategory.UNSPLICED)

        assert models.n_observations == 2
        assert models.global_model.n_observations == 2
        assert models.category_models[CountCategory.SPLICED_ANNOT].n_observations == 1
        assert models.category_models[CountCategory.UNSPLICED].n_observations == 1
        assert models.category_models[CountCategory.INTRON].n_observations == 0
        assert models.intergenic.n_observations == 0

    def test_observe_intergenic(self):
        from hulkrna.insert_model import InsertSizeModels

        models = InsertSizeModels()
        models.observe(150, count_cat=None)

        assert models.n_observations == 1
        assert models.global_model.n_observations == 1
        assert models.intergenic.n_observations == 1
        assert models.intergenic.counts[150] == 1.0

    def test_observe_routing(self):
        """Each observation goes to global + exactly one sub-model."""
        from hulkrna.insert_model import InsertSizeModels
        from hulkrna.core import CountCategory

        models = InsertSizeModels()
        models.observe(200, count_cat=CountCategory.INTRON)
        models.observe(300, count_cat=CountCategory.SPLICED_UNANNOT)
        models.observe(100, count_cat=None)  # intergenic

        assert models.n_observations == 3
        # Each category model has exactly 1 or 0
        assert models.category_models[CountCategory.INTRON].n_observations == 1
        assert models.category_models[CountCategory.SPLICED_UNANNOT].n_observations == 1
        assert models.category_models[CountCategory.SPLICED_ANNOT].n_observations == 0
        assert models.category_models[CountCategory.UNSPLICED].n_observations == 0
        assert models.intergenic.n_observations == 1

    def test_to_dict(self):
        from hulkrna.insert_model import InsertSizeModels
        from hulkrna.core import CountCategory

        models = InsertSizeModels()
        models.observe(200, count_cat=CountCategory.SPLICED_ANNOT)
        d = models.to_dict()

        assert "global" in d
        assert "intergenic" in d
        assert "spliced_annot" in d
        assert "intron" in d
        assert "unspliced" in d
        assert "spliced_unannot" in d
        assert d["global"]["summary"]["n_observations"] == 1
        assert d["spliced_annot"]["summary"]["n_observations"] == 1

    def test_write_json(self, tmp_path):
        import json
        from hulkrna.insert_model import InsertSizeModels
        from hulkrna.core import CountCategory

        models = InsertSizeModels()
        models.observe(200, count_cat=CountCategory.UNSPLICED)
        json_path = tmp_path / "insert_models.json"
        models.write_json(json_path)

        assert json_path.exists()
        with open(json_path) as f:
            data = json.load(f)
        assert "insert_size_models" in data
        assert "global" in data["insert_size_models"]
        assert "intergenic" in data["insert_size_models"]
        assert "unspliced" in data["insert_size_models"]
