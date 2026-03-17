"""Tests for rigel.calibration - EM two-component gDNA deconvolution."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration import (
    GDNACalibration,
    _compute_density_llr_gaussian,
    _compute_strand_llr_betabinom,
    _compute_strand_llr_binomial,
    _e_step,
    _m_step,
    _seed_initial_partition,
    build_gdna_fl_model,
    calibrate_gdna,
    compute_log_density,
    compute_region_stats,
    compute_sense_fraction,
    estimate_kappa_marginal,
)


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
    return pd.DataFrame({
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
    })


def _make_fl_table(region_ids, frag_lens):
    return pd.DataFrame({
        "region_id": np.asarray(region_ids, dtype=np.int32),
        "frag_len": np.asarray(frag_lens, dtype=np.int32),
    })


def _density_setup(stats):
    """Helper: compute log_d and eligible mask."""
    eligible = (stats["n_total"] > 0) & (stats["region_length"] > 0)
    log_d, epsilon = compute_log_density(stats, eligible)
    return log_d, eligible


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
            2, tx_pos=np.array([False, False]), tx_neg=np.array([False, False]),
        )
        counts = _make_region_counts(n_unspliced_pos=[10, 10], n_unspliced_neg=[10, 10])
        stats = compute_region_stats(counts, region_df)
        assert stats["gene_strand"][0] == 0
        assert stats["gene_strand"][1] == 0


# ===================================================================
# TestStrandModel (empirical histograms + binomial bootstrap)
# ===================================================================


class TestStrandModel:

    def test_sense_fraction_plus_gene(self):
        """+ gene: RNA on − strand → sense_frac = 1 − strand_ratio."""
        counts = _make_region_counts(
            n_unspliced_pos=[10], n_unspliced_neg=[90],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, rdf)
        sf = compute_sense_fraction(stats)
        # strand_ratio = 10/100 = 0.1; gene_strand = +1 → sf = 1 − 0.1 = 0.9
        assert sf[0] == pytest.approx(0.9)

    def test_sense_fraction_minus_gene(self):
        """− gene: RNA on + strand → sense_frac = strand_ratio."""
        counts = _make_region_counts(
            n_unspliced_pos=[90], n_unspliced_neg=[10],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([False]), tx_neg=np.array([True]),
        )
        stats = compute_region_stats(counts, rdf)
        sf = compute_sense_fraction(stats)
        # strand_ratio = 90/100 = 0.9; gene_strand = -1 → sf = 0.9
        assert sf[0] == pytest.approx(0.9)

    def test_sense_fraction_ambiguous_is_nan(self):
        counts = _make_region_counts(
            n_unspliced_pos=[70], n_unspliced_neg=[30],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([True]),
        )
        stats = compute_region_stats(counts, rdf)
        sf = compute_sense_fraction(stats)
        # gene_strand=0 (ambiguous) → sense_frac is NaN (no strand signal)
        assert np.isnan(sf[0])

    def test_binomial_llr_zero_when_unstranded(self):
        """When SS=0.5, Binomial LLR is exactly zero for all regions."""
        counts = _make_region_counts(
            n_unspliced_pos=[50], n_unspliced_neg=[50],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, rdf)
        llr = _compute_strand_llr_binomial(stats, 0.5, 1)
        assert llr[0] == pytest.approx(0.0, abs=1e-10)

    def test_binomial_symmetric_favors_gdna(self):
        """Symmetric counts (50/50) at SS=0.95 → positive LLR (gDNA)."""
        counts = _make_region_counts(
            n_unspliced_pos=[50], n_unspliced_neg=[50],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, rdf)
        llr = _compute_strand_llr_binomial(stats, 0.95, 1)
        assert llr[0] > 0

    def test_binomial_biased_favors_rna(self):
        """Strong antisense bias at SS=0.95 → negative LLR (RNA)."""
        counts = _make_region_counts(
            n_unspliced_pos=[5], n_unspliced_neg=[95],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, rdf)
        llr = _compute_strand_llr_binomial(stats, 0.95, 1)
        assert llr[0] < 0

    def test_binomial_scales_with_n(self):
        """LLR magnitude should increase with sample size."""
        # Symmetric strand (sense_frac ≈ 0.5), different n
        counts_small = _make_region_counts(
            n_unspliced_pos=[50], n_unspliced_neg=[50],
        )
        counts_large = _make_region_counts(
            n_unspliced_pos=[500], n_unspliced_neg=[500],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([False]),
        )
        stats_s = compute_region_stats(counts_small, rdf)
        stats_l = compute_region_stats(counts_large, rdf)
        llr_s = _compute_strand_llr_binomial(stats_s, 0.95, 1)
        llr_l = _compute_strand_llr_binomial(stats_l, 0.95, 1)
        # Both positive (gDNA-like: sf near 0.5), larger n → stronger signal
        assert llr_s[0] > 0
        assert llr_l[0] > 0
        assert abs(llr_l[0]) > abs(llr_s[0])


# ===================================================================
# TestDensityLLR — empirical histogram density channel
# ===================================================================


class TestLogDensityAndGaussianLLR:
    """Tests for compute_log_density and _compute_density_llr_gaussian."""

    def test_compute_log_density_basic(self):
        # Two regions: n_total = [100, 0], lengths = [1000, 1000].
        counts = _make_region_counts(
            n_unspliced_pos=[50, 0], n_unspliced_neg=[50, 0],
        )
        stats = compute_region_stats(counts, _make_region_df(2, lengths=np.array([1000, 1000])))
        eligible = stats["n_total"] > 0
        log_d, eps = compute_log_density(stats, eligible)
        # ε = 1/1000 = 0.001
        assert eps == pytest.approx(0.001, rel=0.01)
        # Region 0: log(100/1000 + 0.001) = log(0.101)
        assert log_d[0] == pytest.approx(np.log(0.101), abs=0.01)
        # Region 1: n_total=0 → not eligible → log_d = 0
        assert log_d[1] == 0.0

    def test_compute_log_density_epsilon_from_median(self):
        # 3 regions with lengths [100, 500, 1000] → median = 500 → ε = 0.002
        counts = _make_region_counts(
            n_unspliced_pos=[10, 10, 10], n_unspliced_neg=[10, 10, 10],
        )
        stats = compute_region_stats(
            counts, _make_region_df(3, lengths=np.array([100, 500, 1000])),
        )
        eligible = stats["n_total"] > 0
        _, eps = compute_log_density(stats, eligible)
        assert eps == pytest.approx(1.0 / 500.0)

    def test_zero_count_regions_share_same_baseline(self):
        # Two zero-count regions with different lengths.
        # Under ε = 1/median(L), both get log(ε).
        counts = _make_region_counts(
            n_unspliced_pos=[0, 0, 50], n_unspliced_neg=[0, 0, 50],
        )
        stats = compute_region_stats(
            counts, _make_region_df(3, lengths=np.array([100, 5000, 1000])),
        )
        eligible = stats["n_total"] > 0
        log_d, eps = compute_log_density(stats, eligible)
        # Zero-count regions are not eligible → log_d stays 0.
        # This is correct: they don't participate in the EM.
        assert log_d[0] == 0.0
        assert log_d[1] == 0.0

    def test_density_llr_low_density_favors_gdna(self):
        # Gaussian: gDNA at low density, RNA at high density.
        # Query a low-density region → should favor gDNA (LLR > 0)
        log_d = np.array([-4.9])
        eligible = np.array([True])
        llr = _compute_density_llr_gaussian(
            log_d, mu_g=-5.0, var_g=0.1, mu_r=0.0, var_r=0.5,
            eligible=eligible, n_regions=1,
        )
        assert llr[0] > 0

    def test_density_llr_high_density_favors_rna(self):
        # Query at the RNA peak → should favor RNA (LLR < 0)
        log_d = np.array([0.5])
        eligible = np.array([True])
        llr = _compute_density_llr_gaussian(
            log_d, mu_g=-5.0, var_g=0.5, mu_r=0.0, var_r=0.5,
            eligible=eligible, n_regions=1,
        )
        assert llr[0] < 0

    def test_density_llr_ineligible_returns_zero(self):
        log_d = np.array([-3.0, 0.0])
        eligible = np.array([False, True])
        llr = _compute_density_llr_gaussian(
            log_d, mu_g=-5.0, var_g=0.5, mu_r=0.0, var_r=0.5,
            eligible=eligible, n_regions=2,
        )
        assert llr[0] == 0.0
        assert llr[1] != 0.0

    def test_extreme_density_finite(self):
        # Very high expression → should still produce finite LLR
        log_d = np.array([3.0])
        eligible = np.array([True])
        llr = _compute_density_llr_gaussian(
            log_d, mu_g=-5.0, var_g=1.0, mu_r=0.0, var_r=1.0,
            eligible=eligible, n_regions=1,
        )
        assert np.isfinite(llr[0])


# ===================================================================
# TestEstimateKappaMarginal
# ===================================================================


class TestEstimateKappaMarginal:

    def test_high_symmetry_high_kappa(self):
        rng = np.random.default_rng(42)
        n_regions = 200
        p_true = rng.beta(50, 50, size=n_regions)
        n_pos = rng.binomial(100, p_true).astype(np.float64)
        n_neg = (100 - n_pos).astype(np.float64)
        stats = compute_region_stats(
            _make_region_counts(n_pos, n_neg), _make_region_df(n_regions),
        )
        # All-gDNA scenario: γ = 1 everywhere
        kappa = estimate_kappa_marginal(stats, np.ones(n_regions), 0.95)
        assert kappa is not None
        assert 30 < kappa < 300

    def test_high_overdispersion_low_kappa(self):
        rng = np.random.default_rng(123)
        n_regions = 200
        p_true = rng.beta(2, 2, size=n_regions)
        n_pos = rng.binomial(100, p_true).astype(np.float64)
        n_neg = (100 - n_pos).astype(np.float64)
        stats = compute_region_stats(
            _make_region_counts(n_pos, n_neg), _make_region_df(n_regions),
        )
        kappa = estimate_kappa_marginal(stats, np.ones(n_regions), 0.95)
        assert kappa is not None
        assert 1 < kappa < 15

    def test_all_rna_scenario(self):
        """When γ = 0 everywhere, κ should still converge (not degenerate)."""
        rng = np.random.default_rng(77)
        n_regions = 200
        kappa_true = 50.0
        ss = 0.95
        # RNA-like: k_sense ~ BetaBin(n, κ·SS, κ·(1-SS)), centered at SS.
        # With gene_strand=+1: k_sense = n_unspliced - n_pos
        p_sense = rng.beta(kappa_true * ss, kappa_true * (1 - ss), size=n_regions)
        k_sense = rng.binomial(100, p_sense)
        n_pos = (100 - k_sense).astype(np.float64)
        n_neg = k_sense.astype(np.float64)
        stats = compute_region_stats(
            _make_region_counts(n_pos, n_neg), _make_region_df(n_regions),
        )
        kappa = estimate_kappa_marginal(stats, np.zeros(n_regions), ss)
        assert kappa is not None
        assert kappa > 10  # should fit the RNA distribution, not collapse

    def test_mixed_scenario(self):
        """Mixed γ values: κ lies between pure-component optima."""
        rng = np.random.default_rng(99)
        ng = 100
        # gDNA-like: symmetric, high dispersion
        p_g = rng.beta(5, 5, size=ng)
        n_pos_g = rng.binomial(100, p_g)
        # RNA-like: stranded, high concentration
        p_e = rng.beta(50, 50, size=ng)
        n_pos_e = rng.binomial(100, p_e)
        n_pos = np.concatenate([n_pos_g, n_pos_e]).astype(np.float64)
        n_neg = np.concatenate([100 - n_pos_g, 100 - n_pos_e]).astype(np.float64)
        stats = compute_region_stats(
            _make_region_counts(n_pos, n_neg), _make_region_df(2 * ng),
        )
        gamma = np.concatenate([np.ones(ng), np.zeros(ng)])
        kappa = estimate_kappa_marginal(stats, gamma, 0.95)
        assert kappa is not None
        assert kappa > 1

    def test_insufficient_data_returns_none(self):
        counts = _make_region_counts(n_unspliced_pos=[50, 50], n_unspliced_neg=[50, 50])
        stats = compute_region_stats(counts, _make_region_df(2))
        kappa = estimate_kappa_marginal(stats, np.ones(2), 0.95)
        assert kappa is None

    def test_unstranded_ss_reduces_to_symmetric(self):
        """With SS=0.5, marginal reduces to symmetric BetaBin (stable)."""
        rng = np.random.default_rng(42)
        n_regions = 200
        p_true = rng.beta(50, 50, size=n_regions)
        n_pos = rng.binomial(100, p_true).astype(np.float64)
        n_neg = (100 - n_pos).astype(np.float64)
        stats = compute_region_stats(
            _make_region_counts(n_pos, n_neg), _make_region_df(n_regions),
        )
        # SS=0.5: both components are identical, γ should be irrelevant
        kappa_g1 = estimate_kappa_marginal(stats, np.ones(n_regions), 0.5)
        kappa_g0 = estimate_kappa_marginal(stats, np.zeros(n_regions), 0.5)
        assert kappa_g1 is not None and kappa_g0 is not None
        np.testing.assert_allclose(kappa_g1, kappa_g0, rtol=0.01)

    def test_sampling_variance_correction(self):
        rng = np.random.default_rng(7)
        n_regions = 500
        p_true = rng.beta(50, 50, size=n_regions)
        n_pos = rng.binomial(10, p_true).astype(np.float64)
        n_neg = (10 - n_pos).astype(np.float64)
        stats = compute_region_stats(
            _make_region_counts(n_pos, n_neg), _make_region_df(n_regions),
        )
        kappa = estimate_kappa_marginal(stats, np.ones(n_regions), 0.95)
        assert kappa is not None
        assert kappa > 20


# ===================================================================
# TestBetaBinomialStrandLLR
# ===================================================================


class TestBetaBinomialStrandLLR:
    """Tests for _compute_strand_llr_betabinom."""

    def test_symmetric_ss_gives_zero_llr(self):
        """With SS=0.5, LLR must be exactly zero."""
        counts = _make_region_counts(
            n_unspliced_pos=[70, 30, 50],
            n_unspliced_neg=[30, 70, 50],
        )
        rdf = _make_region_df(
            3, tx_pos=np.array([True, True, True]),
            tx_neg=np.array([False, False, False]),
        )
        stats = compute_region_stats(counts, rdf)
        kappa = 20.0
        llr = _compute_strand_llr_betabinom(stats, 0.5, kappa, 3)
        np.testing.assert_allclose(llr, 0.0, atol=1e-10)

    def test_positive_llr_for_symmetric_strand(self):
        """A region with 50/50 strand ratio should favor gDNA (LLR > 0) at high SS."""
        counts = _make_region_counts(
            n_unspliced_pos=[50], n_unspliced_neg=[50],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, rdf)
        llr = _compute_strand_llr_betabinom(stats, 0.95, 20.0, 1)
        assert llr[0] > 0  # favors gDNA

    def test_negative_llr_for_stranded_region(self):
        """A strongly stranded region should favor RNA (LLR < 0) at high SS."""
        # gene_strand=+1, TruSeq: sense = neg strand → k_sense = n_unspliced - n_pos
        # n_pos=5, n_unspliced=100 → k_sense=95 (high sense fraction → RNA)
        counts = _make_region_counts(
            n_unspliced_pos=[5], n_unspliced_neg=[95],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, rdf)
        llr = _compute_strand_llr_betabinom(stats, 0.95, 20.0, 1)
        assert llr[0] < 0  # favors RNA

    def test_ambiguous_gene_strand_gives_zero_llr(self):
        """Regions with gene_strand=0 should produce LLR=0."""
        counts = _make_region_counts(
            n_unspliced_pos=[80], n_unspliced_neg=[20],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([True]),  # both → strand=0
        )
        stats = compute_region_stats(counts, rdf)
        llr = _compute_strand_llr_betabinom(stats, 0.95, 20.0, 1)
        assert llr[0] == 0.0

    def test_low_count_regions_produce_weak_llr(self):
        """Regions with very few reads should produce near-zero LLR."""
        counts = _make_region_counts(
            n_unspliced_pos=[1], n_unspliced_neg=[1],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, rdf)
        llr = _compute_strand_llr_betabinom(stats, 0.95, 20.0, 1)
        assert abs(llr[0]) < 2.0  # should be modest, not extreme

    def test_higher_kappa_makes_llr_more_extreme(self):
        """Higher κ (less overdispersion) should give stronger strand signal."""
        counts = _make_region_counts(
            n_unspliced_pos=[5], n_unspliced_neg=[95],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, rdf)
        llr_low = _compute_strand_llr_betabinom(stats, 0.95, 10.0, 1)
        llr_high = _compute_strand_llr_betabinom(stats, 0.95, 100.0, 1)
        # Both negative (RNA), higher kappa should be more negative
        assert llr_high[0] < llr_low[0]

    def test_betabinom_vs_binomial_agreement_at_high_kappa(self):
        """At very high κ, Beta-Binomial should approximate Binomial."""
        counts = _make_region_counts(
            n_unspliced_pos=[20, 80, 50],
            n_unspliced_neg=[80, 20, 50],
        )
        rdf = _make_region_df(
            3, tx_pos=np.array([True, True, True]),
            tx_neg=np.array([False, False, False]),
        )
        stats = compute_region_stats(counts, rdf)
        big_kappa = 1e6
        llr_bb = _compute_strand_llr_betabinom(stats, 0.9, big_kappa, 3)
        llr_bin = _compute_strand_llr_binomial(stats, 0.9, 3)
        np.testing.assert_allclose(llr_bb, llr_bin, rtol=0.05)

    def test_e_step_uses_betabinom_when_kappa_provided(self):
        """E-step should use Beta-Binomial when kappa is given."""
        counts = _make_region_counts(
            n_unspliced_pos=[50, 5], n_unspliced_neg=[50, 95],
        )
        rdf = _make_region_df(
            2, tx_pos=np.array([True, True]),
            tx_neg=np.array([False, False]),
        )
        stats = compute_region_stats(counts, rdf)
        log_d, eligible = _density_setup(stats)
        ld_elig = log_d[eligible]
        mu = float(np.mean(ld_elig)) if eligible.any() else 0.0
        var = max(float(np.var(ld_elig)), 1e-12) if eligible.any() else 1.0

        # Without kappa (Binomial fallback)
        gamma_bin = _e_step(
            stats, 0.5, log_d, eligible, 0.95,
            mu_g=mu, var_g=var, mu_r=mu, var_r=var,
        )
        # With kappa (Beta-Binomial)
        gamma_bb = _e_step(
            stats, 0.5, log_d, eligible, 0.95,
            mu_g=mu, var_g=var, mu_r=mu, var_r=var,
            kappa=20.0,
        )
        # Both should classify similarly (gDNA-like vs RNA-like patterns)
        # but Beta-Binomial should be less extreme (more regularized)
        assert gamma_bin[0] > gamma_bin[1]  # region 0 more gDNA-like
        assert gamma_bb[0] > gamma_bb[1]


# ===================================================================
# TestInitialization
# ===================================================================


class TestInitialization:
    """Tests for _seed_initial_partition."""

    def test_spliced_regions_seed_expressed(self):
        # Regions with spliced reads → gamma = 0 (expressed seed).
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50, 50, 50],
            n_unspliced_neg=[50, 50, 50, 50],
            n_spliced_pos=[10, 0, 5, 0],
        )
        stats = compute_region_stats(counts, _make_region_df(4, lengths=np.full(4, 1000)))
        log_d, eligible = _density_setup(stats)
        sf = compute_sense_fraction(stats)
        gamma, pi_init, diag = _seed_initial_partition(
            stats, log_d, sf, 0.95, eligible,
        )
        # Spliced regions get gamma = 0
        assert gamma[0] == 0.0
        assert gamma[2] == 0.0
        assert diag["n_expressed_seed"] >= 2

    def test_gdna_seed_from_low_density(self):
        # Low-density unspliced regions should be seeded as gDNA.
        rng = np.random.default_rng(42)
        n = 200
        # 100 low-density + 100 high-density
        n_pos = np.concatenate([
            rng.binomial(5, 0.5, 100),      # low count
            rng.binomial(100, 0.5, 100),     # high count
        ]).astype(np.float32)
        n_neg = np.concatenate([
            rng.binomial(5, 0.5, 100),
            rng.binomial(100, 0.5, 100),
        ]).astype(np.float32)
        # Give high-density regions spliced reads (expressed seed)
        sp = np.concatenate([
            np.zeros(100, dtype=np.float32),
            rng.poisson(5, 100).astype(np.float32),
        ])
        counts = _make_region_counts(n_pos, n_neg, sp)
        stats = compute_region_stats(counts, _make_region_df(n, lengths=np.full(n, 1000)))
        log_d, eligible = _density_setup(stats)
        sf = compute_sense_fraction(stats)
        gamma, pi_init, diag = _seed_initial_partition(
            stats, log_d, sf, 0.95, eligible,
        )
        assert diag["n_gdna_seed"] > 0
        assert 0.0 < pi_init < 1.0

    def test_min_gdna_regions_guarantee(self):
        # Even when density percentile selects fewer, min_gdna_regions applies.
        counts = _make_region_counts(
            n_unspliced_pos=[50] * 200, n_unspliced_neg=[50] * 200,
            n_spliced_pos=[10] * 100 + [0] * 100,
        )
        stats = compute_region_stats(
            counts, _make_region_df(200, lengths=np.full(200, 1000)),
        )
        log_d, eligible = _density_setup(stats)
        sf = compute_sense_fraction(stats)
        _, _, diag = _seed_initial_partition(
            stats, log_d, sf, 0.95, eligible, min_gdna_regions=50,
        )
        assert diag["n_gdna_seed"] >= 50

    def test_pristine_detection(self):
        # All regions have spliced reads → almost no gDNA seed.
        counts = _make_region_counts(
            n_unspliced_pos=[50] * 20, n_unspliced_neg=[50] * 20,
            n_spliced_pos=[10] * 20,
        )
        stats = compute_region_stats(
            counts, _make_region_df(20, lengths=np.full(20, 1000)),
        )
        log_d, eligible = _density_setup(stats)
        sf = compute_sense_fraction(stats)
        _, pi_init, diag = _seed_initial_partition(
            stats, log_d, sf, 0.95, eligible,
        )
        # All spliced → pi very small or flagged pristine
        assert pi_init < 0.1 or diag["pristine_sample"]

    def test_strand_filter_supplements_density(self):
        # Strand-symmetric unspliced regions get added to gDNA seed.
        # + genes: sense_frac = 1 - strand_ratio; for sf ≈ 0.5, need strand_ratio ≈ 0.5
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50, 5, 5],
            n_unspliced_neg=[50, 50, 95, 95],
            n_spliced_pos=[0, 0, 10, 10],
        )
        rdf = _make_region_df(
            4, tx_pos=np.array([True, True, True, True]),
            tx_neg=np.array([False, False, False, False]),
            lengths=np.full(4, 1000),
        )
        stats = compute_region_stats(counts, rdf)
        log_d, eligible = _density_setup(stats)
        sf = compute_sense_fraction(stats)
        gamma, _, diag = _seed_initial_partition(
            stats, log_d, sf, 0.95, eligible,
        )
        # Regions 0,1 are strand-symmetric → seeded as gDNA
        assert gamma[0] == 1.0 or gamma[1] == 1.0


# ===================================================================
# TestEStep
# ===================================================================


class TestEStep:

    def _call_e_step(self, stats, pi=0.5, ss=0.95):
        """Helper that computes density Gaussian + calls _e_step."""
        log_d, eligible = _density_setup(stats)
        n = len(stats["n_total"])
        # Use broad uninformative Gaussian (equal for both components)
        ld_elig = log_d[eligible]
        if eligible.any():
            mu = float(np.mean(ld_elig))
            var = max(float(np.var(ld_elig)), 1e-12)
        else:
            mu, var = 0.0, 1.0
        return _e_step(
            stats, pi, log_d, eligible, ss,
            mu_g=mu, var_g=var, mu_r=mu, var_r=var,
        )

    def test_spliced_regions_get_gamma_zero(self):
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50], n_unspliced_neg=[50, 50],
            n_spliced_pos=[10, 0],
        )
        stats = compute_region_stats(counts, _make_region_df(2))
        gamma = self._call_e_step(stats, pi=0.5)
        assert gamma[0] == 0.0
        assert gamma[1] > 0.0

    def test_zero_coverage_gets_gamma_pi(self):
        # Zero-count → ineligible → γ = π (prior).
        counts = _make_region_counts(n_unspliced_pos=[0, 50], n_unspliced_neg=[0, 50])
        stats = compute_region_stats(counts, _make_region_df(2))
        gamma = self._call_e_step(stats, pi=0.3)
        assert gamma[0] == pytest.approx(0.3)

    def test_gdna_like_region_high_gamma(self):
        # + gene, symmetric strand (sf ≈ 0.5), SS=0.95 → gDNA-like
        counts = _make_region_counts(n_unspliced_pos=[50], n_unspliced_neg=[50])
        rdf = _make_region_df(1, tx_pos=np.array([True]), tx_neg=np.array([False]),
                              lengths=np.array([10000]))
        stats = compute_region_stats(counts, rdf)
        gamma = self._call_e_step(stats, pi=0.5)
        assert gamma[0] > 0.5

    def test_rna_like_region_low_gamma(self):
        # + gene, SS=0.95: RNA reads on − strand → n_pos LOW (R1-antisense)
        counts = _make_region_counts(n_unspliced_pos=[5], n_unspliced_neg=[95])
        rdf = _make_region_df(1, tx_pos=np.array([True]), tx_neg=np.array([False]),
                              lengths=np.array([1000]))
        stats = compute_region_stats(counts, rdf)
        gamma = self._call_e_step(stats, pi=0.5, ss=0.95)
        assert gamma[0] < 0.3

    def test_strand_signal_classifies_rna_correctly(self):
        """Binomial strand LLR correctly identifies RNA-like regions."""
        # + gene at SS=0.95: RNA reads on − strand → n_pos LOW
        counts = _make_region_counts(n_unspliced_pos=[5], n_unspliced_neg=[95])
        rdf = _make_region_df(1, tx_pos=np.array([True]), tx_neg=np.array([False]),
                              lengths=np.array([1000]))
        stats = compute_region_stats(counts, rdf)
        gamma = self._call_e_step(stats, pi=0.5, ss=0.95)
        # This RNA-like region should have low gamma
        assert gamma[0] < 0.5

    def test_prior_odds_affect_gamma(self):
        # Ambiguous gene_strand → no strand signal; density neutral; prior matters
        counts = _make_region_counts(n_unspliced_pos=[50], n_unspliced_neg=[50])
        rdf = _make_region_df(1, tx_pos=np.ones(1, dtype=bool),
                              tx_neg=np.ones(1, dtype=bool), lengths=np.array([1000]))
        stats = compute_region_stats(counts, rdf)
        g_low = self._call_e_step(stats, pi=0.2)
        g_high = self._call_e_step(stats, pi=0.8)
        assert g_high[0] > g_low[0]


# ===================================================================
# TestMStep
# ===================================================================


class TestMStep:

    def test_basic_m_step(self):
        # gDNA regions: low count.  RNA regions: high count.
        counts = _make_region_counts(
            n_unspliced_pos=[10, 10, 200, 200],
            n_unspliced_neg=[10, 10, 200, 200],
        )
        stats = compute_region_stats(counts, _make_region_df(4, lengths=np.full(4, 1000)))
        log_d, eligible = _density_setup(stats)
        gamma = np.array([1.0, 1.0, 0.0, 0.0])
        pi, lG, lE, mu_g, var_g, mu_r, var_r = _m_step(
            stats, gamma, log_d, eligible,
        )
        assert pi == pytest.approx(0.5)
        assert lG == pytest.approx(0.02)
        assert lE > lG
        # Gaussian params: gDNA mean should be lower than RNA mean
        assert mu_g < mu_r
        assert var_g > 0
        assert var_r > 0

    def test_density_output_types(self):
        counts = _make_region_counts(n_unspliced_pos=[50, 50], n_unspliced_neg=[50, 50])
        stats = compute_region_stats(counts, _make_region_df(2, lengths=np.full(2, 1000)))
        log_d, eligible = _density_setup(stats)
        gamma = np.array([1.0, 1.0])
        pi, lG, lE, mu_g, var_g, mu_r, var_r = _m_step(
            stats, gamma, log_d, eligible,
        )
        assert isinstance(lG, float)
        assert isinstance(lE, float)
        assert isinstance(mu_g, float)
        assert isinstance(var_g, float)

    def test_all_expressed_m_step(self):
        counts = _make_region_counts(n_unspliced_pos=[50, 50], n_unspliced_neg=[50, 50])
        stats = compute_region_stats(counts, _make_region_df(2, lengths=np.full(2, 1000)))
        log_d, eligible = _density_setup(stats)
        gamma = np.array([0.0, 0.0])
        pi, lG, lE, _, _, _, _ = _m_step(
            stats, gamma, log_d, eligible,
        )
        assert pi < 0.1
        assert lE >= lG


# ===================================================================
# TestBuildGDNAFLModel
# ===================================================================


class TestBuildGDNAFLModel:

    def test_basic_fl_model(self):
        model = build_gdna_fl_model(
            np.array([0, 0, 1, 1, 2], dtype=np.int32),
            np.array([200, 250, 200, 300, 200], dtype=np.int32),
            np.array([1.0, 0.5, 0.0]),
        )
        assert model.counts[200] == pytest.approx(1.5)
        assert model.counts[250] == pytest.approx(1.0)
        assert model.counts[300] == pytest.approx(0.5)

    def test_empty_fl_table(self):
        model = build_gdna_fl_model(
            np.array([], dtype=np.int32), np.array([], dtype=np.int32), np.ones(5),
        )
        assert model._finalized
        assert model.total_weight == 0.0

    def test_zero_weights_excluded(self):
        model = build_gdna_fl_model(
            np.array([0, 1], dtype=np.int32),
            np.array([200, 300], dtype=np.int32),
            np.array([0.0, 0.0]),
        )
        assert model.total_weight == 0.0

    def test_fl_out_of_range_excluded(self):
        model = build_gdna_fl_model(
            np.array([0, 0, 0], dtype=np.int32),
            np.array([0, 500, 1500], dtype=np.int32),
            np.array([1.0]),
            max_fl=1000,
        )
        assert model.counts[500] == pytest.approx(1.0)
        assert model.counts[0] == 0.0


# ===================================================================
# TestCalibrateGDNA (integration)
# ===================================================================


class TestCalibrateGDNA:

    def _make_synthetic_data(
        self, n_gdna=100, n_rna=100, kappa_true=50.0, ss=0.95,
        n_per=100, rng_seed=42,
    ):
        rng = np.random.default_rng(rng_seed)
        n_total = n_gdna + n_rna
        alpha = kappa_true / 2.0
        p_gdna = rng.beta(alpha, alpha, size=n_gdna)
        n_pos_gdna = rng.binomial(n_per, p_gdna)
        n_neg_gdna = n_per - n_pos_gdna

        # R1-antisense convention: + gene → RNA reads on − strand → p(+) = 1−SS
        n_pos_rna = rng.binomial(n_per, 1.0 - ss, size=n_rna)
        n_neg_rna = n_per - n_pos_rna
        # Spliced reads follow same antisense convention
        sp_neg_rna = rng.poisson(20, size=n_rna).astype(int)
        sp_pos_rna = rng.poisson(2, size=n_rna).astype(int)

        n_pos = np.concatenate([n_pos_gdna, n_pos_rna]).astype(np.float32)
        n_neg = np.concatenate([n_neg_gdna, n_neg_rna]).astype(np.float32)
        sp_pos = np.concatenate([np.zeros(n_gdna, dtype=int), sp_pos_rna]).astype(np.float32)
        sp_neg = np.concatenate([np.zeros(n_gdna, dtype=int), sp_neg_rna]).astype(np.float32)

        rc = _make_region_counts(n_pos, n_neg, sp_pos, sp_neg)
        rdf = _make_region_df(
            n_total, tx_pos=np.ones(n_total, dtype=bool),
            tx_neg=np.zeros(n_total, dtype=bool), lengths=np.full(n_total, 1000),
        )

        fl_ids_g = np.repeat(np.arange(n_gdna), 5)
        fl_lens_g = np.clip(rng.normal(200, 20, size=len(fl_ids_g)).astype(int), 50, 500)
        fl_ids_r = np.repeat(np.arange(n_gdna, n_total), 5)
        fl_lens_r = np.clip(rng.normal(300, 30, size=len(fl_ids_r)).astype(int), 100, 600)

        fl = _make_fl_table(
            np.concatenate([fl_ids_g, fl_ids_r]),
            np.concatenate([fl_lens_g, fl_lens_r]),
        )
        return rc, fl, rdf, kappa_true

    def test_convergence(self):
        rc, fl, rdf, _ = self._make_synthetic_data()
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert isinstance(result, GDNACalibration)
        assert result.converged
        assert result.n_iterations <= 50

    def test_kappa_recovery(self):
        rc, fl, rdf, kappa_true = self._make_synthetic_data(kappa_true=50.0)
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert kappa_true / 3 < result.kappa < kappa_true * 3
        assert result.kappa > 0

    def test_posteriors_gdna_high_rna_low(self):
        rc, fl, rdf, _ = self._make_synthetic_data(n_gdna=100, n_rna=100)
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        gdna_gamma = result.region_posteriors[:100]
        rna_gamma = result.region_posteriors[100:]
        assert gdna_gamma.mean() > rna_gamma.mean()
        assert (gdna_gamma > 0.5).sum() > 50
        assert (rna_gamma < 0.5).sum() > 50

    def test_mixing_proportion_recovered(self):
        rc, fl, rdf, _ = self._make_synthetic_data(n_gdna=100, n_rna=100)
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert 0.2 < result.mixing_proportion < 0.8

    def test_per_ref_density(self):
        rc, fl, rdf, _ = self._make_synthetic_data(n_gdna=50, n_rna=50)
        rdf = rdf.copy()
        rdf["ref"] = ["chr1"] * 50 + ["chr2"] * 50
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert isinstance(result.gdna_density_per_ref, dict)
        assert "chr1" in result.gdna_density_per_ref
        assert "chr2" in result.gdna_density_per_ref

    def test_gdna_fl_model_has_data(self):
        rc, fl, rdf, _ = self._make_synthetic_data()
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert result.gdna_fl_model.total_weight > 0
        assert result.gdna_fl_model._finalized

    def test_unstranded_convergence(self):
        rng = np.random.default_rng(7)
        n_pos_g = rng.binomial(100, 0.5, size=100).astype(np.float32)
        n_neg_g = (100 - n_pos_g).astype(np.float32)
        n_pos_r = rng.binomial(500, 0.5, size=100).astype(np.float32)
        n_neg_r = (500 - n_pos_r).astype(np.float32)
        sp_r = rng.poisson(30, size=100).astype(np.float32)
        rc = _make_region_counts(
            np.concatenate([n_pos_g, n_pos_r]),
            np.concatenate([n_neg_g, n_neg_r]),
            np.concatenate([np.zeros(100, dtype=np.float32), sp_r]),
        )
        rdf = _make_region_df(200, lengths=np.full(200, 1000))
        result = calibrate_gdna(rc, _make_fl_table([], []), rdf, strand_specificity=0.5)
        assert isinstance(result, GDNACalibration)
        assert result.n_iterations >= 1

    def test_all_spliced_no_gdna(self):
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50], n_unspliced_neg=[50, 50],
            n_spliced_pos=[10, 10],
        )
        result = calibrate_gdna(counts, _make_fl_table([], []),
                                _make_region_df(2), strand_specificity=0.95)
        assert result.region_posteriors.max() == 0.0

    def test_single_region(self):
        counts = _make_region_counts(n_unspliced_pos=[50], n_unspliced_neg=[50])
        result = calibrate_gdna(counts, _make_fl_table([0], [200]),
                                _make_region_df(1), strand_specificity=0.95)
        assert isinstance(result, GDNACalibration)

    def test_diagnostics_populated(self):
        rc, fl, rdf, _ = self._make_synthetic_data()
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95, diagnostics=True)
        assert result.region_stats is not None
        assert isinstance(result.iteration_history, list)
        assert len(result.iteration_history) == result.n_iterations


# ===================================================================
# TestEdgeCases
# ===================================================================


class TestEdgeCases:

    def test_all_gdna(self):
        rng = np.random.default_rng(42)
        n = 50
        p = rng.beta(25, 25, size=n)
        n_pos = rng.binomial(100, p).astype(np.float32)
        n_neg = (100 - n_pos).astype(np.float32)
        rc = _make_region_counts(n_pos, n_neg)
        rdf = _make_region_df(n, tx_pos=np.ones(n, dtype=bool),
                              tx_neg=np.zeros(n, dtype=bool))
        fl = _make_fl_table(np.repeat(np.arange(n), 3), rng.integers(150, 300, size=n * 3))
        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert result.region_posteriors.mean() > 0.5

    def test_zero_length_regions(self):
        counts = _make_region_counts(n_unspliced_pos=[50, 10], n_unspliced_neg=[50, 10])
        rdf = _make_region_df(2, lengths=np.array([0, 1000]))
        result = calibrate_gdna(counts, _make_fl_table([], []), rdf, strand_specificity=0.95)
        assert np.all(np.isfinite(result.region_posteriors))

    def test_all_zero_coverage(self):
        counts = _make_region_counts(n_unspliced_pos=[0, 0], n_unspliced_neg=[0, 0])
        result = calibrate_gdna(counts, _make_fl_table([], []),
                                _make_region_df(2), strand_specificity=0.95)
        assert result.gdna_density_global == 0.0
        assert result.converged
        assert result.n_iterations == 0

    def test_reproducibility(self):
        counts = _make_region_counts(
            n_unspliced_pos=[50, 90, 50, 10],
            n_unspliced_neg=[50, 10, 50, 90],
        )
        rdf = _make_region_df(
            4, tx_pos=np.array([True, True, False, False]),
            tx_neg=np.array([False, False, True, True]),
        )
        fl = _make_fl_table([0, 1, 2, 3], [200, 250, 300, 350])
        r1 = calibrate_gdna(counts, fl, rdf, strand_specificity=0.95)
        r2 = calibrate_gdna(counts, fl, rdf, strand_specificity=0.95)
        assert r1.kappa == r2.kappa
        assert r1.n_iterations == r2.n_iterations
        np.testing.assert_array_equal(r1.region_posteriors, r2.region_posteriors)


# ===================================================================
# TestSignalCombination
# ===================================================================


class TestSignalCombination:

    def _call_e_step(self, stats, pi=0.5, ss=0.95):
        """Helper that computes density Gaussian + calls _e_step."""
        log_d, eligible = _density_setup(stats)
        n = len(stats["n_total"])
        ld_elig = log_d[eligible]
        if eligible.any():
            mu = float(np.mean(ld_elig))
            var = max(float(np.var(ld_elig)), 1e-12)
        else:
            mu, var = 0.0, 1.0
        return _e_step(
            stats, pi, log_d, eligible, ss,
            mu_g=mu, var_g=var, mu_r=mu, var_r=var,
        )

    def test_strand_contributes_zero_when_unstranded(self):
        """With SS=0.5, Binomial strand LLR should be exactly zero."""
        counts = _make_region_counts(
            n_unspliced_pos=[50], n_unspliced_neg=[50],
        )
        rdf = _make_region_df(
            1, tx_pos=np.array([True]), tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, rdf)
        llr = _compute_strand_llr_binomial(stats, 0.5, 1)
        assert abs(llr[0]) < 1e-10

    def test_strand_and_density_reinforce(self):
        # Region 0: symmetric strand + low density on + gene → gDNA-like (high γ)
        # Region 1: biased strand (R1-antisense: n_pos LOW for + gene RNA) + high count → RNA-like
        # Need many regions so histograms carry information
        rng = np.random.default_rng(42)
        n_gdna, n_rna = 50, 50
        # gDNA: symmetric strand, low count
        n_pos_g = rng.binomial(20, 0.5, n_gdna).astype(np.float32)
        n_neg_g = (20 - n_pos_g).astype(np.float32)
        # RNA (+ gene, SS=0.95): n_pos LOW (antisense), high count
        n_pos_r = rng.binomial(200, 0.05, n_rna).astype(np.float32)
        n_neg_r = (200 - n_pos_r).astype(np.float32)
        sp_r = rng.poisson(10, n_rna).astype(np.float32)
        counts = _make_region_counts(
            np.concatenate([n_pos_g, n_pos_r]),
            np.concatenate([n_neg_g, n_neg_r]),
            np.concatenate([np.zeros(n_gdna, dtype=np.float32), sp_r]),
        )
        n_total = n_gdna + n_rna
        rdf = _make_region_df(
            n_total, tx_pos=np.ones(n_total, dtype=bool),
            tx_neg=np.zeros(n_total, dtype=bool),
            lengths=np.full(n_total, 1000),
        )
        result = calibrate_gdna(counts, _make_fl_table([], []), rdf, strand_specificity=0.95)
        # gDNA regions should have higher γ than RNA regions on average
        assert result.region_posteriors[:n_gdna].mean() > result.region_posteriors[n_gdna:].mean()

    def test_full_integration_three_groups(self):
        rng = np.random.default_rng(42)
        p_a = rng.beta(25, 25, size=100)
        npa = rng.binomial(100, p_a).astype(np.float32)
        nna = (100 - npa).astype(np.float32)
        spa = np.zeros(100, dtype=np.float32)

        npb = rng.binomial(300, 0.95, size=100).astype(np.float32)
        nnb = (300 - npb).astype(np.float32)
        spb = rng.poisson(30, size=100).astype(np.float32)

        p_c = rng.beta(25, 25, size=100)
        npc = rng.binomial(100, p_c).astype(np.float32)
        nnc = (100 - npc).astype(np.float32)
        spc = rng.poisson(5, size=100).astype(np.float32)

        rc = _make_region_counts(
            np.concatenate([npa, npb, npc]),
            np.concatenate([nna, nnb, nnc]),
            np.concatenate([spa, spb, spc]),
        )
        rdf = _make_region_df(
            300, tx_pos=np.ones(300, dtype=bool),
            tx_neg=np.zeros(300, dtype=bool), lengths=np.full(300, 1000),
        )
        result = calibrate_gdna(rc, _make_fl_table([], []), rdf, strand_specificity=0.95)
        wa = result.region_posteriors[:100].mean()
        wb = result.region_posteriors[100:200].mean()
        wc = result.region_posteriors[200:].mean()
        assert wa > wc > wb
