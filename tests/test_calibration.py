"""Tests for rigel.calibration — IRLS + Tukey biweight (v4)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration import (
    CalibrationResult,
    build_gdna_fl_model,
    calibrate_gdna,
    compute_region_stats,
    compute_sense_fraction,
)
from rigel.types import STRAND_NEG, STRAND_POS


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
    mappable_effective_length=None,
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
    if mappable_effective_length is not None:
        d["mappable_effective_length"] = np.asarray(mappable_effective_length, dtype=np.float32)
    return pd.DataFrame(d)


def _make_region_counts(
    n_unspliced_pos, n_unspliced_neg, n_spliced_pos=None, n_spliced_neg=None
):
    n = len(n_unspliced_pos)
    return pd.DataFrame(
        {
            "region_id": np.arange(n, dtype=np.int32),
            "n_unspliced_pos": np.asarray(n_unspliced_pos, dtype=np.float32),
            "n_unspliced_neg": np.asarray(n_unspliced_neg, dtype=np.float32),
            "n_spliced_pos": np.asarray(
                n_spliced_pos if n_spliced_pos is not None else np.zeros(n),
                dtype=np.float32,
            ),
            "n_spliced_neg": np.asarray(
                n_spliced_neg if n_spliced_neg is not None else np.zeros(n),
                dtype=np.float32,
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


def _calibrate(
    rc, fl, rdf, *, strand_specificity=0.95, mean_frag_len=200.0,
    intergenic_fl_model=None, fl_prior_ess=1000.0, diagnostics=False,
):
    return calibrate_gdna(
        rc, fl, rdf,
        strand_specificity,
        mean_frag_len=mean_frag_len,
        intergenic_fl_model=intergenic_fl_model,
        fl_prior_ess=fl_prior_ess,
        diagnostics=diagnostics,
    )


# ===========================================================================
# Region statistics
# ===========================================================================


class TestComputeRegionStats:
    def test_basic_computation(self):
        counts = _make_region_counts(
            n_unspliced_pos=[80, 50, 0],
            n_unspliced_neg=[20, 50, 0],
            n_spliced_pos=[10, 0, 0],
            n_spliced_neg=[0, 0, 0],
        )
        rdf = _make_region_df(
            3,
            tx_pos=np.array([True, False, True]),
            tx_neg=np.array([False, True, True]),
            lengths=np.array([1000, 2000, 500]),
        )
        s = compute_region_stats(counts, rdf)
        assert s["strand_ratio"][0] == pytest.approx(0.8)
        assert s["strand_ratio"][1] == pytest.approx(0.5)
        assert np.isnan(s["strand_ratio"][2])
        assert s["splice_rate"][0] == pytest.approx(10 / 110)
        assert s["density"][0] == pytest.approx(110 / 1000)
        assert s["tx_strand"][0] == 1
        assert s["tx_strand"][1] == -1
        assert s["tx_strand"][2] == 0

    def test_exposure_uses_mappable_when_present(self):
        counts = _make_region_counts(n_unspliced_pos=[10], n_unspliced_neg=[10])
        rdf = _make_region_df(
            1, lengths=np.array([1000]), mappable_effective_length=np.array([400.0])
        )
        s = compute_region_stats(counts, rdf)
        assert s["mappable_bp"][0] == pytest.approx(400.0)
        assert s["region_length"][0] == pytest.approx(1000.0)

    def test_exposure_falls_back_to_length(self):
        counts = _make_region_counts(n_unspliced_pos=[10], n_unspliced_neg=[10])
        rdf = _make_region_df(1, lengths=np.array([800]))
        s = compute_region_stats(counts, rdf)
        assert s["mappable_bp"][0] == pytest.approx(800.0)


class TestSenseFraction:
    def test_plus_gene(self):
        counts = _make_region_counts(n_unspliced_pos=[10], n_unspliced_neg=[90])
        rdf = _make_region_df(1, tx_pos=np.array([True]), tx_neg=np.array([False]))
        sf = compute_sense_fraction(compute_region_stats(counts, rdf))
        assert sf[0] == pytest.approx(0.9)

    def test_minus_gene(self):
        counts = _make_region_counts(n_unspliced_pos=[90], n_unspliced_neg=[10])
        rdf = _make_region_df(1, tx_pos=np.array([False]), tx_neg=np.array([True]))
        sf = compute_sense_fraction(compute_region_stats(counts, rdf))
        assert sf[0] == pytest.approx(0.9)

    def test_ambiguous_is_nan(self):
        counts = _make_region_counts(n_unspliced_pos=[70], n_unspliced_neg=[30])
        rdf = _make_region_df(1, tx_pos=np.array([True]), tx_neg=np.array([True]))
        sf = compute_sense_fraction(compute_region_stats(counts, rdf))
        assert np.isnan(sf[0])


# ===========================================================================
# gDNA FL model
# ===========================================================================


class TestBuildGDNAFLModel:
    def _stats(self, n, *, tx_strand=None, length=1000.0, n_unspliced=100.0):
        return {
            "tx_strand": np.full(n, 1, dtype=np.int8) if tx_strand is None
                            else np.asarray(tx_strand, dtype=np.int8),
            "region_length": np.full(n, length, dtype=np.float64),
            "n_unspliced": np.full(n, n_unspliced, dtype=np.float64),
            "n_total": np.full(n, n_unspliced, dtype=np.float64),
        }

    def test_returns_model_with_antisense_data(self):
        n = 3
        stats = self._stats(n)
        rids = np.array([0, 0, 1, 1, 2], dtype=np.int32)
        flens = np.array([200, 250, 200, 300, 200], dtype=np.int32)
        fstrands = np.full(5, STRAND_POS, dtype=np.int8)  # antisense for +1 genes
        fl = _make_fl_table(rids, flens, fstrands)
        model = build_gdna_fl_model(
            fl, stats, region_weight=np.ones(n),
            strand_specificity=0.95, intergenic_fl_model=None, fl_prior_ess=0.0,
        )
        assert model is not None
        assert model.total_weight >= 5

    def test_weighted_by_region_weight(self):
        # Region 0 has weight 0, so its fragments contribute zero.
        stats = self._stats(2)
        rids = np.array([0, 0, 1], dtype=np.int32)
        flens = np.array([200, 200, 200], dtype=np.int32)
        fstrands = np.full(3, STRAND_POS, dtype=np.int8)
        fl = _make_fl_table(rids, flens, fstrands)
        model = build_gdna_fl_model(
            fl, stats, region_weight=np.array([0.0, 1.0]),
            strand_specificity=0.95, intergenic_fl_model=None, fl_prior_ess=0.0,
        )
        # Only region 1's single fragment retained.
        assert model.total_weight == pytest.approx(1.0)

    def test_returns_uniform_when_no_data_no_prior(self):
        stats = self._stats(2)
        model = build_gdna_fl_model(
            _make_fl_table([], []), stats, region_weight=np.ones(2),
            strand_specificity=0.95, intergenic_fl_model=None, fl_prior_ess=0.0,
        )
        # No-cliff philosophy: a uniform fallback is always returned.
        assert model is not None
        assert model.total_weight == 0.0


# ===========================================================================
# calibrate_gdna — orchestration
# ===========================================================================


class TestCalibrateGDNA:
    def _synth(self, n_gdna=100, n_rna=100, ss=0.95, n_per=100, rng_seed=42):
        rng = np.random.default_rng(rng_seed)
        n_total = n_gdna + n_rna
        n_pos_g = rng.binomial(n_per, 0.5, size=n_gdna)
        n_neg_g = n_per - n_pos_g
        n_pos_r = rng.binomial(n_per, 1.0 - ss, size=n_rna)
        n_neg_r = n_per - n_pos_r
        sp_pos_r = rng.poisson(2, size=n_rna).astype(int)
        sp_neg_r = rng.poisson(20, size=n_rna).astype(int)
        rc = _make_region_counts(
            np.concatenate([n_pos_g, n_pos_r]),
            np.concatenate([n_neg_g, n_neg_r]),
            np.concatenate([np.zeros(n_gdna, int), sp_pos_r]),
            np.concatenate([np.zeros(n_gdna, int), sp_neg_r]),
        )
        rdf = _make_region_df(n_total, lengths=np.full(n_total, 1000))
        fl_ids_g = np.repeat(np.arange(n_gdna), 5)
        fl_lens_g = np.clip(rng.normal(200, 20, size=len(fl_ids_g)), 50, 500).astype(int)
        fl_strands_g = rng.choice([STRAND_POS, STRAND_NEG], size=len(fl_ids_g)).astype(np.int8)
        fl = _make_fl_table(fl_ids_g, fl_lens_g, fl_strands_g)
        return rc, fl, rdf

    def test_returns_calibration_result(self):
        rc, fl, rdf = self._synth()
        r = _calibrate(rc, fl, rdf)
        assert isinstance(r, CalibrationResult)

    def test_shapes_and_bounds(self):
        rc, fl, rdf = self._synth()
        r = _calibrate(rc, fl, rdf)
        assert len(r.region_e_gdna) == 200
        assert len(r.region_n_total) == 200
        assert r.lambda_gdna >= 0
        assert 0.5 <= r.strand_specificity <= 1.0

    def test_gdna_detected_in_mixture(self):
        rc, fl, rdf = self._synth(n_gdna=100, n_rna=100)
        r = _calibrate(rc, fl, rdf)
        assert r.lambda_gdna > 0
        assert r.region_e_gdna.sum() > 0

    def test_diagnostics_populated(self):
        rc, fl, rdf = self._synth()
        r = _calibrate(rc, fl, rdf, strand_specificity=0.95, diagnostics=True)
        # v4 EM diagnostics.
        assert r.region_gamma is not None
        assert r.region_gamma.shape == (200,)
        assert r.mu_R is not None
        assert r.sigma_R is not None and r.sigma_R > 0
        assert r.mixing_pi is not None and 0.0 <= r.mixing_pi <= 1.0
        assert r.mixing_pi_soft is not None and 0.0 <= r.mixing_pi_soft <= 1.0
        assert r.em_n_iter >= 1
        assert r.n_eligible > 0

    def test_unstranded_still_converges(self):
        rc, fl, rdf = self._synth()
        r = _calibrate(rc, fl, rdf, strand_specificity=0.5, diagnostics=True)
        # EM still runs; strand channel contributes 0 but count+FL carry signal.
        assert r.em_n_iter >= 1
        assert r.lambda_gdna >= 0

    def test_e_gdna_capped_at_n_total(self):
        rng = np.random.default_rng(0)
        n = 30
        n_pos = rng.binomial(100, rng.uniform(0.3, 0.7, n)).astype(np.float32)
        n_neg = (100 - n_pos).astype(np.float32)
        rc = _make_region_counts(n_pos, n_neg)
        rdf = _make_region_df(n, lengths=np.full(n, 1000))
        r = _calibrate(rc, _make_fl_table([], []), rdf)
        s = compute_region_stats(rc, rdf)
        np.testing.assert_array_less(r.region_e_gdna - 1e-10, s["n_total"] + 1e-10)

    def test_to_summary_dict(self):
        rc, fl, rdf = self._synth()
        r = _calibrate(rc, fl, rdf)
        d = r.to_summary_dict()
        for key in (
            "lambda_gdna", "mu_R", "sigma_R", "mixing_pi", "mixing_pi_soft",
            "strand_used", "strand_z",
            "em_n_iter", "em_converged", "n_eligible", "n_soft",
            "n_spliced_hard", "total_expected_gdna",
            "gdna_fraction", "strand_specificity",
        ):
            assert key in d

    def test_reproducibility(self):
        rc, fl, rdf = self._synth(rng_seed=123)
        r1 = _calibrate(rc, fl, rdf)
        r2 = _calibrate(rc, fl, rdf)
        np.testing.assert_array_equal(r1.region_e_gdna, r2.region_e_gdna)
        assert r1.lambda_gdna == r2.lambda_gdna

    def test_all_zero_coverage(self):
        rc = _make_region_counts(n_unspliced_pos=[0, 0], n_unspliced_neg=[0, 0])
        r = _calibrate(rc, _make_fl_table([], []), _make_region_df(2))
        assert r.lambda_gdna == 0.0
        assert np.all(r.region_e_gdna == 0.0)

    def test_short_regions_excluded_from_density(self):
        # Region 0 length < mean_frag_len is ineligible; region 1 is eligible.
        rc = _make_region_counts(n_unspliced_pos=[5, 5], n_unspliced_neg=[5, 5])
        rdf = _make_region_df(2, lengths=np.array([50, 5000]))
        r = _calibrate(rc, _make_fl_table([], []), rdf, mean_frag_len=200.0)
        assert np.all(np.isfinite(r.region_e_gdna))


# ===========================================================================
# Edge cases
# ===========================================================================


class TestEdgeCases:
    def test_zero_length_regions(self):
        rc = _make_region_counts(n_unspliced_pos=[50, 10], n_unspliced_neg=[50, 10])
        rdf = _make_region_df(2, lengths=np.array([0, 1000]))
        r = _calibrate(rc, _make_fl_table([], []), rdf)
        assert np.all(np.isfinite(r.region_e_gdna))

    def test_single_region(self):
        rc = _make_region_counts(n_unspliced_pos=[50], n_unspliced_neg=[50])
        rdf = _make_region_df(1, lengths=np.array([5000]))
        r = _calibrate(rc, _make_fl_table([0], [200]), rdf)
        assert isinstance(r, CalibrationResult)
        assert len(r.region_e_gdna) == 1


# ===========================================================================
# Direction tests on the new estimator
# ===========================================================================


class TestDirectional:
    def test_higher_density_yields_higher_lambda(self):
        # Two scenarios with identical structure but different density.
        n = 50
        rdf = _make_region_df(n, lengths=np.full(n, 5000))
        # Symmetric strand → both pathways agree gDNA-like.
        rc_low = _make_region_counts(
            n_unspliced_pos=np.full(n, 5, dtype=np.float32),
            n_unspliced_neg=np.full(n, 5, dtype=np.float32),
        )
        rc_hi = _make_region_counts(
            n_unspliced_pos=np.full(n, 50, dtype=np.float32),
            n_unspliced_neg=np.full(n, 50, dtype=np.float32),
        )
        r_low = _calibrate(rc_low, _make_fl_table([], []), rdf)
        r_hi = _calibrate(rc_hi, _make_fl_table([], []), rdf)
        assert r_hi.lambda_gdna > r_low.lambda_gdna

    def test_stranded_clean_rna_yields_lower_lambda(self):
        # Highly stranded RNA-like data should give λ̂ near zero.
        n = 50
        rdf = _make_region_df(n, lengths=np.full(n, 5000))
        # Strong sense bias = RNA. With +1 gene, sense = NEG strand.
        rc = _make_region_counts(
            n_unspliced_pos=np.full(n, 1, dtype=np.float32),
            n_unspliced_neg=np.full(n, 99, dtype=np.float32),
        )
        r = _calibrate(rc, _make_fl_table([], []), rdf, strand_specificity=0.99)
        # With strand LLR active, expressed regions pushed away from γ=1;
        # λ̂ absorbs only the antisense leak ≈ (1-SS) * rate, which is small.
        assert r.lambda_gdna < 0.005
