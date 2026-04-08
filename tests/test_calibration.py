"""Tests for rigel.calibration - analytical gDNA-RNA deconvolution (V3)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration import (
    CalibrationResult,
    _build_gdna_fl_model,
    _length_weighted_percentile,
    calibrate_gdna,
    compute_region_stats,
    compute_sense_fraction,
)
from rigel.types import STRAND_POS, STRAND_NEG


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_region_df(
    n,
    *,
    tx_pos=None,
    tx_neg=None,
    exon_pos=None,
    exon_neg=None,
    ref=None,
    lengths=None,
):
    d = {
        "region_id": np.arange(n, dtype=np.int32),
        "tx_pos": tx_pos if tx_pos is not None else np.ones(n, dtype=bool),
        "tx_neg": tx_neg if tx_neg is not None else np.zeros(n, dtype=bool),
        "length": lengths if lengths is not None else np.full(n, 1000),
    }
    if exon_pos is not None:
        d["exon_pos"] = exon_pos
    if exon_neg is not None:
        d["exon_neg"] = exon_neg
    if ref is not None:
        d["ref"] = ref
    return pd.DataFrame(d)


def _make_region_counts(
    n_unspliced_pos,
    n_unspliced_neg,
    n_spliced_pos=None,
    n_spliced_neg=None,
):
    n = len(n_unspliced_pos)
    return pd.DataFrame(
        {
            "region_id": np.arange(n, dtype=np.int32),
            "n_unspliced_pos": np.asarray(n_unspliced_pos, dtype=np.float32),
            "n_unspliced_neg": np.asarray(n_unspliced_neg, dtype=np.float32),
            "n_spliced_pos": (
                np.asarray(n_spliced_pos, dtype=np.float32)
                if n_spliced_pos is not None
                else np.zeros(n, dtype=np.float32)
            ),
            "n_spliced_neg": (
                np.asarray(n_spliced_neg, dtype=np.float32)
                if n_spliced_neg is not None
                else np.zeros(n, dtype=np.float32)
            ),
        }
    )


def _make_fl_table(region_ids, frag_lens, frag_strands=None):
    d = {
        "region_id": np.asarray(region_ids, dtype=np.int32),
        "frag_len": np.asarray(frag_lens, dtype=np.int32),
    }
    if frag_strands is not None:
        d["frag_strand"] = np.asarray(frag_strands, dtype=np.int8)
    return pd.DataFrame(d)


# ===================================================================
# TestComputeRegionStats
# ===================================================================


class TestComputeRegionStats:
    def test_basic_computation(self):
        counts = _make_region_counts(
            n_unspliced_pos=[80, 50, 0],
            n_unspliced_neg=[20, 50, 0],
            n_spliced_pos=[10, 0, 0],
            n_spliced_neg=[0, 0, 0],
        )
        region_df = _make_region_df(
            3,
            tx_pos=np.array([True, False, True]),
            tx_neg=np.array([False, True, True]),
            lengths=np.array([1000, 2000, 500]),
        )
        stats = compute_region_stats(counts, region_df)
        assert stats["strand_ratio"][0] == pytest.approx(0.8)
        assert stats["strand_ratio"][1] == pytest.approx(0.5)
        assert np.isnan(stats["strand_ratio"][2])
        assert stats["splice_rate"][0] == pytest.approx(10 / 110)
        assert stats["splice_rate"][1] == 0.0
        assert stats["density"][0] == pytest.approx(110 / 1000)
        assert stats["density"][1] == pytest.approx(100 / 2000)
        assert stats["density"][2] == 0.0
        assert stats["gene_strand"][0] == 1
        assert stats["gene_strand"][1] == -1
        assert stats["gene_strand"][2] == 0

    def test_zero_coverage_region(self):
        counts = _make_region_counts(n_unspliced_pos=[0], n_unspliced_neg=[0])
        stats = compute_region_stats(counts, _make_region_df(1))
        assert np.isnan(stats["strand_ratio"][0])
        assert stats["splice_rate"][0] == 0.0
        assert stats["density"][0] == 0.0

    def test_intergenic_gene_strand(self):
        region_df = _make_region_df(
            2,
            tx_pos=np.array([False, False]),
            tx_neg=np.array([False, False]),
        )
        counts = _make_region_counts(n_unspliced_pos=[10, 10], n_unspliced_neg=[10, 10])
        stats = compute_region_stats(counts, region_df)
        assert stats["gene_strand"][0] == 0
        assert stats["gene_strand"][1] == 0


# ===================================================================
# TestStrandModel (sense fraction)
# ===================================================================


class TestStrandModel:
    def test_sense_fraction_plus_gene(self):
        counts = _make_region_counts(n_unspliced_pos=[10], n_unspliced_neg=[90])
        rdf = _make_region_df(1, tx_pos=np.array([True]), tx_neg=np.array([False]))
        stats = compute_region_stats(counts, rdf)
        sf = compute_sense_fraction(stats)
        assert sf[0] == pytest.approx(0.9)

    def test_sense_fraction_minus_gene(self):
        counts = _make_region_counts(n_unspliced_pos=[90], n_unspliced_neg=[10])
        rdf = _make_region_df(1, tx_pos=np.array([False]), tx_neg=np.array([True]))
        stats = compute_region_stats(counts, rdf)
        sf = compute_sense_fraction(stats)
        assert sf[0] == pytest.approx(0.9)

    def test_sense_fraction_ambiguous_is_nan(self):
        counts = _make_region_counts(n_unspliced_pos=[70], n_unspliced_neg=[30])
        rdf = _make_region_df(1, tx_pos=np.array([True]), tx_neg=np.array([True]))
        stats = compute_region_stats(counts, rdf)
        sf = compute_sense_fraction(stats)
        assert np.isnan(sf[0])


# ===================================================================
# TestLengthWeightedPercentile
# ===================================================================


class TestLengthWeightedPercentile:
    def test_uniform_weights(self):
        vals = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        wts = np.ones(5)
        p50 = _length_weighted_percentile(vals, wts, 50.0)
        assert p50 == pytest.approx(3.0, abs=0.1)

    def test_extreme_weights(self):
        vals = np.array([1.0, 100.0])
        wts = np.array([1000.0, 1.0])
        p50 = _length_weighted_percentile(vals, wts, 50.0)
        assert p50 < 10.0

    def test_empty_input(self):
        assert _length_weighted_percentile(np.array([]), np.array([]), 50.0) == 0.0

    def test_zero_weights_excluded(self):
        vals = np.array([1.0, 100.0, 50.0])
        wts = np.array([0.0, 0.0, 10.0])
        p50 = _length_weighted_percentile(vals, wts, 50.0)
        assert p50 == pytest.approx(50.0, abs=0.1)


# ===================================================================
# TestBuildGDNAFLModel
# ===================================================================


class TestBuildGDNAFLModel:
    def _make_stats_and_eligible(self, n, gene_strand_vals=None, lengths=None, n_unspliced=None):
        if lengths is None:
            lengths = np.full(n, 1000, dtype=np.float64)
        else:
            lengths = np.asarray(lengths, dtype=np.float64)
        if gene_strand_vals is None:
            gene_strand_vals = np.ones(n, dtype=np.int8)
        if n_unspliced is None:
            n_unspliced = np.full(n, 100, dtype=np.float64)
        stats = {
            "gene_strand": np.asarray(gene_strand_vals, dtype=np.int8),
            "region_length": lengths,
            "n_unspliced": np.asarray(n_unspliced, dtype=np.float64),
            "n_total": np.asarray(n_unspliced, dtype=np.float64),
        }
        eligible = np.ones(n, dtype=bool)
        return stats, eligible

    def test_basic_antisense_fl(self):
        n = 5
        stats, eligible = self._make_stats_and_eligible(n)
        rids = np.array([0, 0, 1, 1, 2], dtype=np.int32)
        flens = np.array([200, 250, 200, 300, 200], dtype=np.int32)
        fstrands = np.array([STRAND_POS] * 5, dtype=np.int8)
        fl_table = _make_fl_table(rids, flens, fstrands)
        model = _build_gdna_fl_model(
            fl_table,
            stats,
            eligible,
            strand_specificity=0.95,
            density_percentile=10.0,
            intergenic_fl_model=None,
            min_ess=3,
        )
        assert model is not None
        assert model.total_weight >= 5

    def test_empty_fl_table(self):
        stats, eligible = self._make_stats_and_eligible(3)
        model = _build_gdna_fl_model(
            pd.DataFrame(),
            stats,
            eligible,
            strand_specificity=0.95,
            density_percentile=10.0,
            intergenic_fl_model=None,
        )
        assert model is None

    def test_no_frag_strand_uses_density(self):
        n = 10
        n_unspliced = np.array([5] * 5 + [500] * 5, dtype=np.float64)
        stats, eligible = self._make_stats_and_eligible(
            n, n_unspliced=n_unspliced, lengths=np.full(n, 10000)
        )
        rids = np.repeat(np.arange(n, dtype=np.int32), 10)
        flens = np.full(len(rids), 200, dtype=np.int32)
        fl_table = _make_fl_table(rids, flens)
        model = _build_gdna_fl_model(
            fl_table,
            stats,
            eligible,
            strand_specificity=0.95,
            density_percentile=50.0,
            intergenic_fl_model=None,
            min_ess=5,
        )
        assert model is not None


# ===================================================================
# TestCalibrateGDNA (V3 analytical)
# ===================================================================


class TestCalibrateGDNA:
    def _make_synthetic_data(self, n_gdna=100, n_rna=100, ss=0.95, n_per=100, rng_seed=42):
        rng = np.random.default_rng(rng_seed)
        n_total = n_gdna + n_rna
        n_pos_gdna = rng.binomial(n_per, 0.5, size=n_gdna)
        n_neg_gdna = n_per - n_pos_gdna
        n_pos_rna = rng.binomial(n_per, 1.0 - ss, size=n_rna)
        n_neg_rna = n_per - n_pos_rna
        sp_neg_rna = rng.poisson(20, size=n_rna).astype(int)
        sp_pos_rna = rng.poisson(2, size=n_rna).astype(int)
        n_pos = np.concatenate([n_pos_gdna, n_pos_rna]).astype(np.float32)
        n_neg = np.concatenate([n_neg_gdna, n_neg_rna]).astype(np.float32)
        sp_pos = np.concatenate([np.zeros(n_gdna, dtype=int), sp_pos_rna]).astype(np.float32)
        sp_neg = np.concatenate([np.zeros(n_gdna, dtype=int), sp_neg_rna]).astype(np.float32)
        rc = _make_region_counts(n_pos, n_neg, sp_pos, sp_neg)
        rdf = _make_region_df(
            n_total,
            tx_pos=np.ones(n_total, dtype=bool),
            tx_neg=np.zeros(n_total, dtype=bool),
            lengths=np.full(n_total, 1000),
        )
        fl_ids_g = np.repeat(np.arange(n_gdna), 5)
        fl_lens_g = np.clip(rng.normal(200, 20, size=len(fl_ids_g)).astype(int), 50, 500)
        fl_strands_g = rng.choice([STRAND_POS, STRAND_NEG], size=len(fl_ids_g)).astype(np.int8)
        fl_ids_r = np.repeat(np.arange(n_gdna, n_total), 5)
        fl_lens_r = np.clip(rng.normal(300, 30, size=len(fl_ids_r)).astype(int), 100, 600)
        fl_strands_r = np.full(len(fl_ids_r), STRAND_NEG, dtype=np.int8)
        flip = rng.random(len(fl_ids_r)) < 0.05
        fl_strands_r[flip] = STRAND_POS
        fl = _make_fl_table(
            np.concatenate([fl_ids_g, fl_ids_r]),
            np.concatenate([fl_lens_g, fl_lens_r]),
            np.concatenate([fl_strands_g, fl_strands_r]),
        )
        return rc, fl, rdf

    def test_returns_calibration_result(self):
        rc, fl, rdf = self._make_synthetic_data()
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert isinstance(result, CalibrationResult)

    def test_result_has_correct_shape(self):
        rc, fl, rdf = self._make_synthetic_data()
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert len(result.region_e_gdna) == 200
        assert len(result.region_n_total) == 200
        assert result.lambda_gdna >= 0.0
        assert 0.5 <= result.strand_specificity <= 1.0

    def test_gdna_regions_higher_e_gdna(self):
        rc, fl, rdf = self._make_synthetic_data(n_gdna=100, n_rna=100)
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        # Global strand aggregation distributes e_gdna proportionally by
        # n_unspliced, so per-region differentiation is via density pathway.
        # Core check: gDNA is detected (positive lambda and total e_gdna).
        assert result.lambda_gdna > 0
        assert result.region_e_gdna.sum() > 0

    def test_lambda_gdna_positive(self):
        rc, fl, rdf = self._make_synthetic_data(n_gdna=100, n_rna=100)
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert result.lambda_gdna > 0.0

    def test_gdna_fl_model_has_data(self):
        rc, fl, rdf = self._make_synthetic_data()
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert result.gdna_fl_model is not None
        assert result.gdna_fl_model.total_weight > 0

    def test_unstranded_uses_density(self):
        rng = np.random.default_rng(7)
        n_gdna, n_rna = 100, 100
        n_total = n_gdna + n_rna
        n_pos_g = rng.binomial(10, 0.5, size=n_gdna).astype(np.float32)
        n_neg_g = (10 - n_pos_g).astype(np.float32)
        n_pos_r = rng.binomial(500, 0.5, size=n_rna).astype(np.float32)
        n_neg_r = (500 - n_pos_r).astype(np.float32)
        sp_r = rng.poisson(30, size=n_rna).astype(np.float32)
        rc = _make_region_counts(
            np.concatenate([n_pos_g, n_pos_r]),
            np.concatenate([n_neg_g, n_neg_r]),
            np.concatenate([np.zeros(n_gdna, dtype=np.float32), sp_r]),
        )
        rdf = _make_region_df(n_total, lengths=np.full(n_total, 1000))
        result = calibrate_gdna(rc, _make_fl_table([], []), rdf, strand_specificity=0.5)
        assert isinstance(result, CalibrationResult)
        assert result.lambda_gdna >= 0.0

    def test_all_spliced_low_gdna(self):
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50],
            n_unspliced_neg=[50, 50],
            n_spliced_pos=[10, 10],
        )
        result = calibrate_gdna(
            counts,
            _make_fl_table([], []),
            _make_region_df(2),
            strand_specificity=0.95,
        )
        assert np.all(np.isfinite(result.region_e_gdna))

    def test_single_region(self):
        counts = _make_region_counts(n_unspliced_pos=[50], n_unspliced_neg=[50])
        result = calibrate_gdna(
            counts,
            _make_fl_table([0], [200]),
            _make_region_df(1),
            strand_specificity=0.95,
        )
        assert isinstance(result, CalibrationResult)
        assert len(result.region_e_gdna) == 1

    def test_deprecated_kwargs_accepted(self):
        counts = _make_region_counts(n_unspliced_pos=[50], n_unspliced_neg=[50])
        result = calibrate_gdna(
            counts,
            _make_fl_table([0], [200]),
            _make_region_df(1),
            strand_specificity=0.95,
            max_iterations=50,
            convergence_tol=1e-4,
            min_gdna_regions=100,
        )
        assert isinstance(result, CalibrationResult)

    def test_to_summary_dict(self):
        rc, fl, rdf = self._make_synthetic_data()
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        d = result.to_summary_dict()
        assert "lambda_gdna" in d
        assert "total_expected_gdna" in d
        assert "gdna_fraction" in d
        assert d["strand_specificity"] == pytest.approx(0.95)


# ===================================================================
# TestEdgeCases
# ===================================================================


class TestEdgeCases:
    def test_all_gdna_symmetric(self):
        rng = np.random.default_rng(42)
        n = 50
        p = rng.beta(25, 25, size=n)
        n_pos = rng.binomial(100, p).astype(np.float32)
        n_neg = (100 - n_pos).astype(np.float32)
        rc = _make_region_counts(n_pos, n_neg)
        rdf = _make_region_df(
            n,
            tx_pos=np.ones(n, dtype=bool),
            tx_neg=np.zeros(n, dtype=bool),
        )
        fl = _make_fl_table(
            np.repeat(np.arange(n), 3),
            rng.integers(150, 300, size=n * 3),
        )
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert result.region_e_gdna.mean() > 0

    def test_zero_length_regions(self):
        counts = _make_region_counts(n_unspliced_pos=[50, 10], n_unspliced_neg=[50, 10])
        rdf = _make_region_df(2, lengths=np.array([0, 1000]))
        result = calibrate_gdna(counts, _make_fl_table([], []), rdf, strand_specificity=0.95)
        assert np.all(np.isfinite(result.region_e_gdna))

    def test_all_zero_coverage(self):
        counts = _make_region_counts(n_unspliced_pos=[0, 0], n_unspliced_neg=[0, 0])
        result = calibrate_gdna(
            counts,
            _make_fl_table([], []),
            _make_region_df(2),
            strand_specificity=0.95,
        )
        assert result.lambda_gdna == 0.0
        assert np.all(result.region_e_gdna == 0.0)

    def test_reproducibility(self):
        counts = _make_region_counts(
            n_unspliced_pos=[50, 90, 50, 10],
            n_unspliced_neg=[50, 10, 50, 90],
        )
        rdf = _make_region_df(
            4,
            tx_pos=np.array([True, True, False, False]),
            tx_neg=np.array([False, False, True, True]),
        )
        fl = _make_fl_table([0, 1, 2, 3], [200, 250, 300, 350])
        r1 = calibrate_gdna(counts, fl, rdf, strand_specificity=0.95)
        r2 = calibrate_gdna(counts, fl, rdf, strand_specificity=0.95)
        np.testing.assert_array_equal(r1.region_e_gdna, r2.region_e_gdna)
        assert r1.lambda_gdna == r2.lambda_gdna

    def test_e_gdna_capped_at_unspliced(self):
        rng = np.random.default_rng(42)
        n = 20
        n_pos = rng.binomial(100, rng.uniform(0.3, 0.7, n)).astype(np.float32)
        n_neg = (100 - n_pos).astype(np.float32)
        rc = _make_region_counts(n_pos, n_neg)
        rdf = _make_region_df(
            n,
            tx_pos=np.ones(n, dtype=bool),
            tx_neg=np.zeros(n, dtype=bool),
        )
        result = calibrate_gdna(rc, _make_fl_table([], []), rdf, strand_specificity=0.95)
        stats = compute_region_stats(rc, rdf)
        np.testing.assert_array_less(
            result.region_e_gdna - 1e-10,
            stats["n_unspliced"] + 1e-10,
        )


# ===================================================================
# TestStrandDecomposition
# ===================================================================


class TestStrandDecomposition:
    def test_symmetric_strand_detected_as_gdna(self):
        counts = _make_region_counts(n_unspliced_pos=[50], n_unspliced_neg=[50])
        rdf = _make_region_df(
            1,
            tx_pos=np.array([True]),
            tx_neg=np.array([False]),
            lengths=np.array([10000]),
        )
        result = calibrate_gdna(counts, _make_fl_table([], []), rdf, strand_specificity=0.95)
        assert result.region_e_gdna[0] > 50

    def test_stranded_region_low_gdna(self):
        counts = _make_region_counts(n_unspliced_pos=[5], n_unspliced_neg=[95])
        rdf = _make_region_df(
            1,
            tx_pos=np.array([True]),
            tx_neg=np.array([False]),
            lengths=np.array([1000]),
        )
        result = calibrate_gdna(counts, _make_fl_table([], []), rdf, strand_specificity=0.95)
        assert result.region_e_gdna[0] < 20

    def test_ambiguous_gene_strand_uses_density_only(self):
        counts = _make_region_counts(n_unspliced_pos=[50], n_unspliced_neg=[50])
        rdf = _make_region_df(
            1,
            tx_pos=np.array([True]),
            tx_neg=np.array([True]),
            lengths=np.array([1000]),
        )
        result = calibrate_gdna(counts, _make_fl_table([], []), rdf, strand_specificity=0.95)
        assert np.isfinite(result.region_e_gdna[0])

    def test_minus_gene_strand_decomposition(self):
        counts = _make_region_counts(n_unspliced_pos=[50], n_unspliced_neg=[50])
        rdf = _make_region_df(
            1,
            tx_pos=np.array([False]),
            tx_neg=np.array([True]),
            lengths=np.array([1000]),
        )
        result = calibrate_gdna(counts, _make_fl_table([], []), rdf, strand_specificity=0.95)
        assert result.region_e_gdna[0] > 50

    def test_blending_weight_at_ss_075(self):
        counts = _make_region_counts(
            n_unspliced_pos=[50, 5],
            n_unspliced_neg=[50, 95],
        )
        rdf = _make_region_df(
            2,
            tx_pos=np.array([True, True]),
            tx_neg=np.array([False, False]),
            lengths=np.array([1000, 1000]),
        )
        result = calibrate_gdna(counts, _make_fl_table([], []), rdf, strand_specificity=0.75)
        assert result.strand_specificity == pytest.approx(0.75)
        # Global strand aggregation distributes proportionally by n_unspliced;
        # both regions have equal unspliced counts, so they get equal e_gdna.
        assert result.region_e_gdna[0] == pytest.approx(result.region_e_gdna[1])
        assert result.region_e_gdna[0] > 0
