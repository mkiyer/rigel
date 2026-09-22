"""`rigel.frag_length_model.FragmentLengthModel` in the two shapes production builds: `from_pmf`, the
scoring and effective-length model built from `calibration.fl`'s pmfs, and the raw-histogram model the
QC report summarises (`FragmentLengthModel(counts=...)`). Gated: construction, the statistics reported
(the censored overflow bin excluded where it must be), the effective length implied, the scoring table,
and the serialization.
"""

import math

import numpy as np
import pytest

from rigel.frag_length_model import FragmentLengthModel


# =====================================================================
# FragmentLengthModel
# =====================================================================


class TestFragmentLengthModelBasic:
    def test_default_construction(self):
        m = FragmentLengthModel()
        assert m.max_size == 1000
        assert m.n_observations == 0
        assert m.counts.shape == (1001,)
        assert m.total_weight == 0.0

    def test_custom_max_size(self):
        m = FragmentLengthModel(max_size=500)
        assert m.counts.shape == (501,)


class TestFragmentLengthModelStatistics:
    @pytest.fixture
    def model_with_data(self):
        counts = np.zeros(501, dtype=np.float64)
        counts[200], counts[250], counts[300] = 3.0, 2.0, 1.0  # 200,200,200,250,250,300
        return FragmentLengthModel(max_size=500, counts=counts)

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
        m = FragmentLengthModel(max_size=100)
        assert m.mean == 0.0
        assert m.std == 0.0
        assert m.median == 0.0
        assert m.mode == 0

    def test_overflow_bin_excluded_from_mean_mode_std(self):
        """mean/std/mode ignore the >= max_size overflow bin; median does not."""
        counts = np.zeros(101, dtype=np.float64)
        counts[50] = 90.0  # in-range bulk
        counts[100] = 10.0  # overflow bin (>= max_size=100), a censored long tail
        m = FragmentLengthModel(max_size=100, counts=counts)  # raw QC view (unfinalized)

        # mean/mode/std describe the in-range bulk (all at 50), NOT dragged to 100.
        assert m.mode == 50
        assert m.mean == pytest.approx(50.0)
        assert m.std == pytest.approx(0.0)
        # median is a rank statistic over the full sample: with only 10% censored,
        # the 50th percentile still lands at the in-range bulk.
        assert m.median == 50.0

    def test_overflow_reported_in_to_dict(self):
        """The censored overflow bin is surfaced as a QC diagnostic."""
        counts = np.zeros(101, dtype=np.float64)
        counts[50] = 90.0
        counts[100] = 10.0
        m = FragmentLengthModel(max_size=100, counts=counts)  # unfinalized (raw QC view)
        summ = m.to_dict()["summary"]
        assert summ["mode"] == 50
        assert summ["mean"] == pytest.approx(50.0, abs=0.01)
        assert summ["overflow"]["count"] == pytest.approx(10.0)
        assert summ["overflow"]["fraction"] == pytest.approx(0.1)


class TestFragmentLengthModelEffectiveLength:
    @staticmethod
    def _oracle_eff_len(probs, length):
        sizes = np.arange(len(probs), dtype=np.float64)
        usable = (sizes > 0) & (sizes <= length)
        return float(np.sum(probs[usable] * (length - sizes[usable] + 1.0)))

    def test_terminal_bin_is_not_double_counted_for_long_transcripts(self):
        pmf = np.zeros(11, dtype=np.float64)
        pmf[10] = 1.0
        m = FragmentLengthModel.from_pmf(pmf, max_size=10)

        probs = np.exp(m._log_prob)
        expected = self._oracle_eff_len(probs, length=20)
        actual = m.compute_all_transcript_eff_lens(np.array([20]), min_value=0.0)[0]

        assert actual == pytest.approx(expected)

    def test_from_pmf_preserves_scoring_distribution(self):
        # mass at length 0 and at the overflow bin, so the zero-length term and the tail base are both
        # distinguishable from their neighbours
        pmf = np.zeros(11, dtype=np.float64)
        pmf[0], pmf[2], pmf[8], pmf[10] = 0.1, 0.2, 0.5, 0.2

        m = FragmentLengthModel.from_pmf(pmf, max_size=10)

        np.testing.assert_allclose(m.pmf, pmf, rtol=0.0, atol=1e-15)
        # the table the scorer reads, and the base its tail decays from
        assert m._log_prob[2] == pytest.approx(math.log(0.2))
        assert m._log_prob[8] == pytest.approx(math.log(0.5))
        assert m._tail_base == pytest.approx(math.log(0.2))
        expected = self._oracle_eff_len(pmf, length=8)
        actual = m.compute_all_transcript_eff_lens(np.array([8]), min_value=0.0)[0]
        assert actual == pytest.approx(expected)


class TestFragmentLengthModelSerialization:
    @staticmethod
    def _two_fragments():
        counts = np.zeros(101, dtype=np.float64)
        counts[50] = counts[60] = 1.0
        return FragmentLengthModel(max_size=100, counts=counts)

    def test_to_dict_structure(self):
        d = self._two_fragments().to_dict()
        assert "summary" in d
        assert "histogram" in d
        assert d["summary"]["n_observations"] == 2
        assert d["summary"]["max_size"] == 100

    def test_to_dict_histogram_trimmed(self):
        d = self._two_fragments().to_dict()
        assert d["histogram"]["range"] == [50, 60]
        assert len(d["histogram"]["values"]) == 11  # 50 through 60 inclusive

    def test_to_dict_empty_histogram(self):
        m = FragmentLengthModel(max_size=100)
        d = m.to_dict()
        assert d["histogram"]["values"] == []


# The report's fragment-length categories come from `FLModels` rather than from a per-splice-category
# histogram of the scanner's, and are gated in tests/test_summary_report.py.
