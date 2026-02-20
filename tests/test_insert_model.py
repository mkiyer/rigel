"""Tests for hulkrna.insert_model — insert size distribution model."""

import json
import math

import numpy as np
import pytest

from hulkrna.categories import SpliceType
from hulkrna.insert_model import InsertSizeModel, InsertSizeModels


# =====================================================================
# InsertSizeModel
# =====================================================================


class TestInsertSizeModelBasic:
    def test_default_construction(self):
        m = InsertSizeModel()
        assert m.max_size == 1000
        assert m.n_observations == 0
        assert m.counts.shape == (1001,)
        assert m.total_weight == 0.0

    def test_custom_max_size(self):
        m = InsertSizeModel(max_size=500)
        assert m.counts.shape == (501,)

    def test_observe_accumulates(self):
        m = InsertSizeModel(max_size=100)
        m.observe(50)
        m.observe(50)
        m.observe(75)
        assert m.n_observations == 3
        assert m.counts[50] == 2.0
        assert m.counts[75] == 1.0

    def test_observe_weight(self):
        m = InsertSizeModel(max_size=100)
        m.observe(50, weight=2.5)
        assert m.n_observations == 1
        assert m.counts[50] == 2.5
        assert m.total_weight == 2.5

    def test_observe_clamps_negative(self):
        m = InsertSizeModel(max_size=100)
        m.observe(-10)
        assert m.counts[0] == 1.0

    def test_observe_clamps_overflow(self):
        m = InsertSizeModel(max_size=100)
        m.observe(999)
        assert m.counts[100] == 1.0  # overflow bin


class TestInsertSizeModelStatistics:
    @pytest.fixture
    def model_with_data(self):
        m = InsertSizeModel(max_size=500)
        # Add observations at specific sizes
        for size in [200, 200, 200, 250, 250, 300]:
            m.observe(size)
        return m

    def test_mean(self, model_with_data):
        # (200*3 + 250*2 + 300*1) / 6 = (600+500+300)/6 = 1400/6 ≈ 233.33
        assert model_with_data.mean == pytest.approx(1400 / 6, rel=1e-6)

    def test_std(self, model_with_data):
        mu = 1400 / 6
        var = (3 * (200 - mu) ** 2 + 2 * (250 - mu) ** 2 + 1 * (300 - mu) ** 2) / 6
        assert model_with_data.std == pytest.approx(math.sqrt(var), rel=1e-6)

    def test_median(self, model_with_data):
        # Sorted: 200,200,200,250,250,300 → median at position 3 → 200 or 250
        # cumsum: [0..199]=0, [200]=3, ..., [250]=5, ... total=6, half=3
        # searchsorted finds first index where cumsum >= 3, which is 200
        assert model_with_data.median == 200.0

    def test_mode(self, model_with_data):
        assert model_with_data.mode == 200

    def test_empty_statistics(self):
        m = InsertSizeModel(max_size=100)
        assert m.mean == 0.0
        assert m.std == 0.0
        assert m.median == 0.0
        assert m.mode == 0


class TestInsertSizeModelLikelihood:
    def test_log_likelihood_uniform_when_empty(self):
        m = InsertSizeModel(max_size=100)
        ll = m.log_likelihood(50)
        assert ll == pytest.approx(-math.log(101))

    def test_log_likelihood_with_data(self):
        m = InsertSizeModel(max_size=100)
        m.observe(50, weight=9.0)
        # count[50] = 9, total = 9, laplace: log((9+1)/(9+101)) = log(10/110)
        expected = math.log(10.0 / 110.0)
        assert m.log_likelihood(50) == pytest.approx(expected)

    def test_log_likelihood_unseen_size(self):
        m = InsertSizeModel(max_size=100)
        m.observe(50, weight=9.0)
        # count[30] = 0, total = 9, laplace: log((0+1)/(9+101)) = log(1/110)
        expected = math.log(1.0 / 110.0)
        assert m.log_likelihood(30) == pytest.approx(expected)

    def test_log_likelihood_is_negative(self):
        m = InsertSizeModel(max_size=100)
        m.observe(50)
        assert m.log_likelihood(50) < 0


class TestInsertSizeModelSerialization:
    def test_to_dict_structure(self):
        m = InsertSizeModel(max_size=100)
        m.observe(50)
        m.observe(60)
        d = m.to_dict()
        assert "summary" in d
        assert "histogram" in d
        assert d["summary"]["n_observations"] == 2
        assert d["summary"]["max_size"] == 100

    def test_to_dict_histogram_trimmed(self):
        m = InsertSizeModel(max_size=100)
        m.observe(50)
        m.observe(60)
        d = m.to_dict()
        assert d["histogram"]["range"] == [50, 60]
        assert len(d["histogram"]["values"]) == 11  # 50 through 60 inclusive

    def test_to_dict_empty_histogram(self):
        m = InsertSizeModel(max_size=100)
        d = m.to_dict()
        assert d["histogram"]["values"] == []

    def test_write_json(self, tmp_path):
        m = InsertSizeModel(max_size=100)
        m.observe(50)
        path = tmp_path / "insert.json"
        m.write_json(path)
        assert path.exists()
        data = json.loads(path.read_text())
        assert "insert_size_model" in data


# =====================================================================
# InsertSizeModels (per-category container)
# =====================================================================


class TestInsertSizeModels:
    def test_construction(self):
        models = InsertSizeModels(max_size=200)
        assert models.global_model.max_size == 200
        assert models.intergenic.max_size == 200
        assert len(models.category_models) == len(SpliceType)

    def test_observe_routes_to_global_and_category(self):
        models = InsertSizeModels()
        models.observe(250, splice_type=SpliceType.SPLICED_ANNOT)
        assert models.global_model.n_observations == 1
        assert models.category_models[SpliceType.SPLICED_ANNOT].n_observations == 1
        assert models.category_models[SpliceType.UNSPLICED].n_observations == 0
        assert models.intergenic.n_observations == 0

    def test_observe_intergenic(self):
        models = InsertSizeModels()
        models.observe(300, splice_type=None)
        assert models.global_model.n_observations == 1
        assert models.intergenic.n_observations == 1
        for cat in SpliceType:
            assert models.category_models[cat].n_observations == 0

    def test_n_observations_delegates_to_global(self):
        models = InsertSizeModels()
        models.observe(100, splice_type=SpliceType.UNSPLICED)
        models.observe(200, splice_type=None)
        assert models.n_observations == 2

    def test_to_dict(self):
        models = InsertSizeModels()
        models.observe(250, splice_type=SpliceType.SPLICED_ANNOT)
        d = models.to_dict()
        assert "global" in d
        assert "intergenic" in d
        assert "spliced_annot" in d
        assert "unspliced" in d

    def test_write_json(self, tmp_path):
        models = InsertSizeModels()
        models.observe(250, splice_type=SpliceType.SPLICED_ANNOT)
        path = tmp_path / "insert_models.json"
        models.write_json(path)
        assert path.exists()
        data = json.loads(path.read_text())
        assert "insert_size_models" in data
