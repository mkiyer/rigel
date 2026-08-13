"""Tests for rigel.strand_model — the strand model and its per-sj SJ strand table."""

import numpy as np
import pytest

from rigel.types import Strand
from rigel.strand_model import SJStrandTable, StrandModel, StrandModels


def _labels(pairs):
    """Expand ``[(align, sj, count), ...]`` into the C++ scanner's two parallel label arrays."""
    align, sj = [], []
    for a, s, n in pairs:
        align += [int(a)] * n
        sj += [int(s)] * n
    return np.asarray(align, dtype=np.int8), np.asarray(sj, dtype=np.int8)


def _model(pairs) -> StrandModel:
    return StrandModel.from_labels(*_labels(pairs))


def _table(rows) -> SJStrandTable:
    """Build a table from ``[(motif_strand, n_sense, n_antisense), ...]``."""
    motif = [int(m) for m, _, _ in rows]
    return SJStrandTable(
        ref_id=np.zeros(len(rows), dtype=np.int32),
        start=np.arange(len(rows), dtype=np.int64) * 1000,
        end=np.arange(len(rows), dtype=np.int64) * 1000 + 100,
        motif_strand=np.asarray(motif, dtype=np.int8),
        n_sense=np.asarray([s for _, s, _ in rows], dtype=np.int64),
        n_antisense=np.asarray([a for _, _, a in rows], dtype=np.int64),
    )


class TestStrandModelCounts:
    """The 2×2 contingency table and derived counts."""

    def test_default_counts_are_zero(self):
        sm = StrandModel()
        assert sm.pos_pos == 0
        assert sm.pos_neg == 0
        assert sm.neg_pos == 0
        assert sm.neg_neg == 0
        assert sm.n_observations == 0
        assert sm.sj_table is None

    def test_from_labels_fills_each_cell(self):
        sm = _model(
            [
                (Strand.POS, Strand.POS, 3),
                (Strand.POS, Strand.NEG, 5),
                (Strand.NEG, Strand.POS, 7),
                (Strand.NEG, Strand.NEG, 11),
            ]
        )
        assert (sm.pos_pos, sm.pos_neg, sm.neg_pos, sm.neg_neg) == (3, 5, 7, 11)
        assert sm.n_same == 3 + 11
        assert sm.n_opposite == 5 + 7
        assert sm.n_observations == 26

    def test_is_immutable(self):
        sm = _model([(Strand.POS, Strand.POS, 1)])
        with pytest.raises(Exception):
            sm.pos_pos = 99


class TestSJStrandTable:
    """The per-sj refinement and the marginal identity that licenses it."""

    def test_empty(self):
        t = SJStrandTable.empty()
        assert t.n_sj == 0
        assert t.n_observations == 0
        assert t.contingency() == (0, 0, 0, 0)

    def test_marginal_is_the_2x2(self):
        """⭐ THE correctness argument: the 2×2 is exactly the table's marginal."""
        t = _table([(Strand.POS, 3, 7), (Strand.POS, 10, 2), (Strand.NEG, 5, 4)])
        pos_pos, pos_neg, neg_pos, neg_neg = t.contingency()
        assert pos_pos == 3 + 10  # sense on motif-POS sj
        assert neg_pos == 7 + 2  # antisense on motif-POS sj
        assert neg_neg == 5  # sense on motif-NEG sj
        assert pos_neg == 4  # antisense on motif-NEG sj

    def test_from_sj_table_agrees_with_from_labels(self):
        """Both constructors, one population: the same fragments give the same 2×2 and κ."""
        rows = [(Strand.POS, 3, 7), (Strand.POS, 10, 2), (Strand.NEG, 5, 4)]
        from_table = StrandModel.from_sj_table(_table(rows))
        # The same fragments as raw labels: motif POS + sense ⇒ (align POS, sj POS), etc.
        from_labels = _model(
            [
                (Strand.POS, Strand.POS, 13),  # sense on motif POS
                (Strand.NEG, Strand.POS, 9),  # antisense on motif POS
                (Strand.NEG, Strand.NEG, 5),  # sense on motif NEG
                (Strand.POS, Strand.NEG, 4),  # antisense on motif NEG
            ]
        )
        assert (from_table.pos_pos, from_table.pos_neg, from_table.neg_pos, from_table.neg_neg) == (
            from_labels.pos_pos,
            from_labels.pos_neg,
            from_labels.neg_pos,
            from_labels.neg_neg,
        )
        assert from_table.p_r1_sense == pytest.approx(from_labels.p_r1_sense)

    def test_depth_and_totals(self):
        t = _table([(Strand.POS, 3, 7), (Strand.NEG, 100, 0)])
        np.testing.assert_array_equal(t.depth, [10, 100])
        assert t.n_observations == 110
        assert t.n_sj == 2
        assert StrandModel.from_sj_table(t).p_r1_sense == pytest.approx(103 / 110)

    def test_to_dict_reports_deep_sj(self):
        t = _table([(Strand.POS, 1, 1), (Strand.POS, 60, 60), (Strand.NEG, 900, 200)])
        d = t.to_dict()
        assert d["n_sj"] == 3
        assert d["n_observations"] == 2 + 120 + 1100
        assert d["n_sj_depth_ge_100"] == 2
        assert d["n_sj_depth_ge_1000"] == 1
        assert d["depth_max"] == 1100

    def test_model_carries_table_through(self):
        t = _table([(Strand.POS, 4, 1)])
        sm = StrandModel.from_sj_table(t)
        assert sm.sj_table is t
        assert sm.n_observations == 5


class TestStrandModelPosterior:
    """MLE probability computation."""

    def test_no_observations(self):
        assert StrandModel().p_r1_sense == 0.5

    def test_strong_fr_library(self):
        """A strongly R1-sense library (most same-direction)."""
        sm = _model([(Strand.POS, Strand.POS, 95), (Strand.POS, Strand.NEG, 5)])
        assert sm.p_r1_sense == pytest.approx(0.95)
        assert sm.strand_specificity > 0.9
        assert sm.read1_sense is True

    def test_strong_rf_library(self):
        """A strongly R1-antisense library (most opposite-direction)."""
        sm = _model([(Strand.POS, Strand.POS, 5), (Strand.POS, Strand.NEG, 95)])
        assert sm.p_r1_sense < 0.1
        assert sm.p_r1_antisense > 0.9
        assert sm.strand_specificity > 0.9
        assert sm.read1_sense is False


class TestStrandModelProperties:
    def test_posterior_variance_with_data(self):
        sm = _model([(Strand.POS, Strand.POS, 80), (Strand.POS, Strand.NEG, 20)])
        # p = 80/100 = 0.8, variance = 0.8 * 0.2 / 100 = 0.0016
        assert sm.posterior_variance() == pytest.approx(0.0016)

    def test_posterior_variance_no_observations(self):
        assert StrandModel().posterior_variance() == 0.25

    def test_posterior_95ci(self):
        sm = _model([(Strand.POS, Strand.POS, 95), (Strand.POS, Strand.NEG, 5)])
        lo, hi = sm.posterior_95ci()
        assert lo < hi
        assert lo > 0.88  # should be heavily skewed towards 1.0
        assert hi <= 1.0


class TestStrandModelsContainer:
    """The StrandModels container and its construction from a scanner dict."""

    def _scan_dict(self, rows, exonic=()):
        t = _table(rows)
        align, sj = _labels(exonic) if exonic else (np.empty(0, np.int8), np.empty(0, np.int8))
        return {
            "sj_ref_id": t.ref_id,
            "sj_start": t.start,
            "sj_end": t.end,
            "sj_motif_strand": t.motif_strand,
            "sj_n_sense": t.n_sense,
            "sj_n_antisense": t.n_antisense,
            "exonic_obs": align,
            "exonic_truth": sj,
        }

    def test_delegation_to_exonic_spliced(self):
        models = StrandModels.from_scan(self._scan_dict([(Strand.POS, 50, 0)]))
        assert models.strand_specificity == models.exonic_spliced.strand_specificity
        assert models.p_r1_sense == models.exonic_spliced.p_r1_sense
        assert models.read1_sense == models.exonic_spliced.read1_sense
        assert models.n_observations == models.exonic_spliced.n_observations

    def test_mle_from_scan(self):
        models = StrandModels.from_scan(self._scan_dict([(Strand.POS, 10, 0)]))
        assert models.p_r1_sense == pytest.approx(1.0)

    def test_spliced_2x2_is_the_table_marginal(self):
        """The container's spliced model has ONE source of truth — the sj table."""
        rows = [(Strand.POS, 30, 3), (Strand.NEG, 12, 5)]
        models = StrandModels.from_scan(self._scan_dict(rows))
        assert models.exonic_spliced.contingency_matches_table()
        assert models.sj_table.n_sj == 2

    def test_zero_observations_warns(self, caplog):
        with caplog.at_level("WARNING"):
            models = StrandModels.from_scan(self._scan_dict([]))
        assert "No spliced strand observations" in caplog.text
        assert models.p_r1_sense == 0.5

    def test_low_observations_warns(self, caplog):
        with caplog.at_level("WARNING"):
            StrandModels.from_scan(self._scan_dict([(Strand.POS, 5, 0)]))
        assert "Only 5 spliced strand observations" in caplog.text

    def test_diagnostic_model_built_from_labels_and_has_no_table(self):
        models = StrandModels.from_scan(
            self._scan_dict([(Strand.POS, 1, 0)], exonic=[(Strand.POS, Strand.POS, 10)])
        )
        assert models.exonic.n_observations == 10
        assert models.exonic.sj_table is None

    def test_default_container_has_an_empty_table(self):
        assert StrandModels().sj_table.n_sj == 0
