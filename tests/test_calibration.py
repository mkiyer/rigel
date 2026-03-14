"""Tests for rigel.calibration — iterative seed-extend gDNA calibration."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration import (
    GDNACalibration,
    build_gdna_fl_model,
    calibrate_gdna,
    compute_region_stats,
    estimate_gdna_density,
    estimate_kappa_sym,
    score_regions,
    select_seed,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_region_df(
    n: int,
    *,
    tx_pos: np.ndarray | None = None,
    tx_neg: np.ndarray | None = None,
    lengths: np.ndarray | None = None,
) -> pd.DataFrame:
    """Build a synthetic region_df with n regions."""
    return pd.DataFrame({
        "region_id": np.arange(n, dtype=np.int32),
        "tx_pos": tx_pos if tx_pos is not None else np.ones(n, dtype=bool),
        "tx_neg": tx_neg if tx_neg is not None else np.zeros(n, dtype=bool),
        "length": lengths if lengths is not None else np.full(n, 1000),
    })


def _make_region_counts(
    n_unspliced_pos: np.ndarray,
    n_unspliced_neg: np.ndarray,
    n_spliced_pos: np.ndarray | None = None,
    n_spliced_neg: np.ndarray | None = None,
) -> pd.DataFrame:
    """Build a synthetic region_counts DataFrame."""
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
    """Build a synthetic FL table."""
    return pd.DataFrame({
        "region_id": np.asarray(region_ids, dtype=np.int32),
        "frag_len": np.asarray(frag_lens, dtype=np.int32),
    })


# ===================================================================
# TestComputeRegionStats
# ===================================================================


class TestComputeRegionStats:
    """Tests for compute_region_stats."""

    def test_basic_computation(self):
        """Correct strand_ratio, splice_rate, density, gene_strand."""
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

        # Strand ratios
        assert stats["strand_ratio"][0] == pytest.approx(0.8)
        assert stats["strand_ratio"][1] == pytest.approx(0.5)
        assert np.isnan(stats["strand_ratio"][2])

        # Splice rates
        assert stats["splice_rate"][0] == pytest.approx(10 / 110)
        assert stats["splice_rate"][1] == 0.0
        assert stats["splice_rate"][2] == 0.0

        # Density
        assert stats["density"][0] == pytest.approx(110 / 1000)
        assert stats["density"][1] == pytest.approx(100 / 2000)
        assert stats["density"][2] == 0.0

        # Gene strand
        assert stats["gene_strand"][0] == 1   # tx_pos only
        assert stats["gene_strand"][1] == -1  # tx_neg only
        assert stats["gene_strand"][2] == 0   # both

    def test_zero_coverage_region(self):
        """Regions with zero counts produce safe defaults."""
        counts = _make_region_counts(
            n_unspliced_pos=[0],
            n_unspliced_neg=[0],
        )
        region_df = _make_region_df(1)
        stats = compute_region_stats(counts, region_df)

        assert np.isnan(stats["strand_ratio"][0])
        assert stats["splice_rate"][0] == 0.0
        assert stats["density"][0] == 0.0

    def test_intergenic_gene_strand(self):
        """Intergenic regions (no tx flags) get gene_strand = 0."""
        region_df = _make_region_df(
            2,
            tx_pos=np.array([False, False]),
            tx_neg=np.array([False, False]),
        )
        counts = _make_region_counts(
            n_unspliced_pos=[10, 10],
            n_unspliced_neg=[10, 10],
        )
        stats = compute_region_stats(counts, region_df)
        assert stats["gene_strand"][0] == 0
        assert stats["gene_strand"][1] == 0


# ===================================================================
# TestSelectSeed
# ===================================================================


class TestSelectSeed:
    """Tests for select_seed (stranded and unstranded)."""

    def test_stranded_symmetric_selected(self):
        """Symmetric, zero-spliced, sufficient coverage → selected."""
        counts = _make_region_counts(
            n_unspliced_pos=[50],
            n_unspliced_neg=[50],
        )
        region_df = _make_region_df(1)
        stats = compute_region_stats(counts, region_df)

        mask = select_seed(stats, strand_specificity=0.95)
        assert mask[0]

    def test_stranded_asymmetric_excluded(self):
        """Highly asymmetric strand ratio → excluded."""
        counts = _make_region_counts(
            n_unspliced_pos=[95],
            n_unspliced_neg=[5],
        )
        region_df = _make_region_df(1)
        stats = compute_region_stats(counts, region_df)

        mask = select_seed(stats, strand_specificity=0.95)
        assert not mask[0]

    def test_stranded_spliced_excluded(self):
        """Any spliced reads → excluded from seed."""
        counts = _make_region_counts(
            n_unspliced_pos=[50],
            n_unspliced_neg=[50],
            n_spliced_pos=[1],
        )
        region_df = _make_region_df(1)
        stats = compute_region_stats(counts, region_df)

        mask = select_seed(stats, strand_specificity=0.95)
        assert not mask[0]

    def test_stranded_low_coverage_excluded(self):
        """Below min_coverage → excluded."""
        counts = _make_region_counts(
            n_unspliced_pos=[3],
            n_unspliced_neg=[3],
        )
        region_df = _make_region_df(1)
        stats = compute_region_stats(counts, region_df)

        mask = select_seed(stats, strand_specificity=0.95, min_coverage=10)
        assert not mask[0]

    def test_stranded_mix(self):
        """Mix of symmetric/asymmetric/spliced → correct selection."""
        counts = _make_region_counts(
            # Region 0: symmetric, no splice, high coverage → seed
            # Region 1: asymmetric, no splice → not seed
            # Region 2: symmetric, has splice → not seed
            # Region 3: symmetric, low coverage → not seed
            n_unspliced_pos=[50, 95, 48, 3],
            n_unspliced_neg=[50, 5, 52, 3],
            n_spliced_pos=[0, 0, 5, 0],
        )
        region_df = _make_region_df(4)
        stats = compute_region_stats(counts, region_df)

        mask = select_seed(stats, strand_specificity=0.95, min_coverage=10)
        assert mask[0] is np.True_
        assert mask[1] is np.False_
        assert mask[2] is np.False_
        assert mask[3] is np.False_

    def test_unstranded_low_density_selected(self):
        """Unstranded mode: low-density, zero-spliced → selected."""
        counts = _make_region_counts(
            # Region 0: low density (10 frags / 10000 bp = 0.001)
            # Region 1: high density (100 frags / 1000 bp = 0.1)
            n_unspliced_pos=[5, 50],
            n_unspliced_neg=[5, 50],
        )
        region_df = _make_region_df(
            2, lengths=np.array([10000, 1000]),
        )
        stats = compute_region_stats(counts, region_df)

        # At 50th percentile, only region 0 (lower density) passes
        mask = select_seed(
            stats, strand_specificity=0.5, coverage_percentile=50.0,
        )
        assert mask[0]
        assert not mask[1]

    def test_unstranded_spliced_excluded(self):
        """Unstranded: zero coverage but has spliced → excluded."""
        counts = _make_region_counts(
            n_unspliced_pos=[10],
            n_unspliced_neg=[10],
            n_spliced_pos=[5],
        )
        region_df = _make_region_df(1, lengths=np.array([10000]))
        stats = compute_region_stats(counts, region_df)

        mask = select_seed(stats, strand_specificity=0.5)
        assert not mask[0]

    def test_empty_seed(self):
        """All regions have spliced reads → empty seed."""
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50],
            n_unspliced_neg=[50, 50],
            n_spliced_pos=[1, 1],
        )
        region_df = _make_region_df(2)
        stats = compute_region_stats(counts, region_df)

        mask = select_seed(stats, strand_specificity=0.95)
        assert mask.sum() == 0


# ===================================================================
# TestEstimateKappaSym
# ===================================================================


class TestEstimateKappaSym:
    """Tests for estimate_kappa_sym (MoM)."""

    def test_high_symmetry_high_kappa(self):
        """Regions all near p=0.5 → high κ (low overdispersion)."""
        rng = np.random.default_rng(42)
        n_regions = 200
        n_per = 100

        # BetaBinomial: p ~ Beta(50, 50) → κ = 100
        # Generate strand ratios close to 0.5
        p_true = rng.beta(50, 50, size=n_regions)
        n_pos = rng.binomial(n_per, p_true).astype(np.float64)
        n_neg = (n_per - n_pos).astype(np.float64)

        counts = _make_region_counts(n_pos, n_neg)
        region_df = _make_region_df(n_regions)
        stats = compute_region_stats(counts, region_df)
        weights = np.ones(n_regions)

        kappa = estimate_kappa_sym(stats, weights)
        # Should recover κ ~ 100 reasonably (within factor of 2)
        assert 30 < kappa < 300

    def test_high_overdispersion_low_kappa(self):
        """Widely varying strand ratios → low κ."""
        rng = np.random.default_rng(123)
        n_regions = 200
        n_per = 100

        # BetaBinomial: p ~ Beta(2, 2) → κ = 4
        p_true = rng.beta(2, 2, size=n_regions)
        n_pos = rng.binomial(n_per, p_true).astype(np.float64)
        n_neg = (n_per - n_pos).astype(np.float64)

        counts = _make_region_counts(n_pos, n_neg)
        region_df = _make_region_df(n_regions)
        stats = compute_region_stats(counts, region_df)
        weights = np.ones(n_regions)

        kappa = estimate_kappa_sym(stats, weights)
        # Should recover κ ~ 4 reasonably
        assert 1 < kappa < 15

    def test_weighted_estimation(self):
        """High-weight regions dominate the κ estimate."""
        rng = np.random.default_rng(99)
        n_per_group = 100
        n_per = 100

        # Group A: κ ~ 100 (tight), weight = 1.0
        p_a = rng.beta(50, 50, size=n_per_group)
        n_pos_a = rng.binomial(n_per, p_a)
        n_neg_a = n_per - n_pos_a

        # Group B: κ ~ 4 (wide), weight = 0.01
        p_b = rng.beta(2, 2, size=n_per_group)
        n_pos_b = rng.binomial(n_per, p_b)
        n_neg_b = n_per - n_pos_b

        n_pos = np.concatenate([n_pos_a, n_pos_b]).astype(np.float64)
        n_neg = np.concatenate([n_neg_a, n_neg_b]).astype(np.float64)

        counts = _make_region_counts(n_pos, n_neg)
        region_df = _make_region_df(2 * n_per_group)
        stats = compute_region_stats(counts, region_df)

        weights = np.concatenate([
            np.ones(n_per_group),
            np.full(n_per_group, 0.01),
        ])

        kappa = estimate_kappa_sym(stats, weights)
        # Should be close to group A's κ ~ 100, not group B's ~ 4
        assert kappa > 30

    def test_fallback_on_insufficient_data(self):
        """Too few weighted regions → return fallback."""
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50],
            n_unspliced_neg=[50, 50],
        )
        region_df = _make_region_df(2)
        stats = compute_region_stats(counts, region_df)
        weights = np.ones(2)

        kappa = estimate_kappa_sym(
            stats, weights, kappa_fallback=42.0, kappa_min_obs=10,
        )
        assert kappa == 42.0

    def test_sampling_variance_correction(self):
        """Low-coverage regions: correction prevents κ underestimate."""
        rng = np.random.default_rng(7)
        n_regions = 500
        n_per = 10  # low coverage → high sampling variance

        # True κ = 100 → tight Beta(50, 50)
        p_true = rng.beta(50, 50, size=n_regions)
        n_pos = rng.binomial(n_per, p_true).astype(np.float64)
        n_neg = (n_per - n_pos).astype(np.float64)

        counts = _make_region_counts(n_pos, n_neg)
        region_df = _make_region_df(n_regions)
        stats = compute_region_stats(counts, region_df)
        weights = np.ones(n_regions)

        kappa = estimate_kappa_sym(stats, weights)
        # Without correction, κ would be much lower.
        # With correction, should be > 20 at least.
        assert kappa > 20


# ===================================================================
# TestEstimateGDNADensity
# ===================================================================


class TestEstimateGDNADensity:
    """Tests for estimate_gdna_density."""

    def test_basic(self):
        """Correct density from uniform regions."""
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50],
            n_unspliced_neg=[50, 50],
        )
        region_df = _make_region_df(2, lengths=np.array([1000, 1000]))
        stats = compute_region_stats(counts, region_df)
        weights = np.ones(2)

        density = estimate_gdna_density(stats, weights)
        # Total unspliced = 200, total bp = 2000
        assert density == pytest.approx(0.1)

    def test_weighted_density(self):
        """Weights affect density estimate."""
        counts = _make_region_counts(
            # Region 0: 100 frags, 1000 bp → density 0.1
            # Region 1: 10 frags, 1000 bp → density 0.01
            n_unspliced_pos=[50, 5],
            n_unspliced_neg=[50, 5],
        )
        region_df = _make_region_df(2, lengths=np.array([1000, 1000]))
        stats = compute_region_stats(counts, region_df)

        # Weight only region 1
        weights = np.array([0.0, 1.0])
        density = estimate_gdna_density(stats, weights)
        assert density == pytest.approx(0.01)

    def test_zero_weights(self):
        """All zero weights → density 0."""
        counts = _make_region_counts(
            n_unspliced_pos=[50],
            n_unspliced_neg=[50],
        )
        region_df = _make_region_df(1)
        stats = compute_region_stats(counts, region_df)
        weights = np.zeros(1)

        density = estimate_gdna_density(stats, weights)
        assert density == 0.0


# ===================================================================
# TestScoreRegions
# ===================================================================


class TestScoreRegions:
    """Tests for score_regions."""

    def test_stranded_gdna_region(self):
        """Symmetric strand ratio + no splicing → high gDNA weight."""
        counts = _make_region_counts(
            n_unspliced_pos=[50],
            n_unspliced_neg=[50],
        )
        # Region is on + strand gene, SS=0.95
        region_df = _make_region_df(
            1,
            tx_pos=np.array([True]),
            tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, region_df)

        # Under gDNA: strand ratio should be ~0.5 ✓
        # Under RNA (+gene, SS=0.95): strand ratio should be ~0.95 ✗
        weights = score_regions(stats, kappa_sym=50.0, strand_specificity=0.95)
        assert weights[0] > 0.8

    def test_stranded_rna_region(self):
        """Biased strand ratio + spliced reads → low gDNA weight."""
        counts = _make_region_counts(
            n_unspliced_pos=[90],
            n_unspliced_neg=[10],
            n_spliced_pos=[20],
        )
        # + strand gene, SS=0.95 → RNA expects ~0.95 strand ratio
        region_df = _make_region_df(
            1,
            tx_pos=np.array([True]),
            tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, region_df)

        weights = score_regions(stats, kappa_sym=50.0, strand_specificity=0.95)
        assert weights[0] < 0.2

    def test_stranded_mixed_regions(self):
        """Verify ordering: gDNA-like > ambiguous > RNA-like."""
        counts = _make_region_counts(
            # Region 0: symmetric, no splice → gDNA-like
            # Region 1: biased, no splice → RNA-like
            # Region 2: symmetric, some splice → mixed
            n_unspliced_pos=[50, 90, 50],
            n_unspliced_neg=[50, 10, 50],
            n_spliced_pos=[0, 0, 20],
        )
        region_df = _make_region_df(
            3,
            tx_pos=np.array([True, True, True]),
            tx_neg=np.array([False, False, False]),
        )
        stats = compute_region_stats(counts, region_df)

        weights = score_regions(stats, kappa_sym=50.0, strand_specificity=0.95)
        assert weights[0] > weights[2] > weights[1]

    def test_unstranded_density_scoring(self):
        """Unstranded mode: low density → high weight, high density → low."""
        counts = _make_region_counts(
            # Region 0: low density (10 frags / 10000 bp = 0.001)
            # Region 1: high density (500 frags / 1000 bp = 0.5)
            n_unspliced_pos=[5, 250],
            n_unspliced_neg=[5, 250],
        )
        region_df = _make_region_df(
            2,
            lengths=np.array([10000, 1000]),
        )
        stats = compute_region_stats(counts, region_df)

        # gdna_density = 0.001 (like region 0)
        weights = score_regions(
            stats, kappa_sym=50.0, strand_specificity=0.5,
            gdna_density=0.001,
        )
        # Region 0 (matches expected density) should have higher weight
        # than region 1 (huge excess)
        assert weights[0] > weights[1]

    def test_splice_rate_adjustment(self):
        """Splice rate reduces gDNA weight proportionally."""
        # Two identical regions, one with spliced reads
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50],
            n_unspliced_neg=[50, 50],
            n_spliced_pos=[0, 50],
        )
        region_df = _make_region_df(
            2,
            tx_pos=np.array([True, True]),
            tx_neg=np.array([False, False]),
        )
        stats = compute_region_stats(counts, region_df)

        weights = score_regions(stats, kappa_sym=50.0, strand_specificity=0.95)
        # Region 1 has splice_rate = 50/150 ≈ 0.33
        # So w1 ≈ w0 * (1 - 0.33)
        assert weights[1] < weights[0]
        ratio = weights[1] / max(weights[0], 1e-10)
        assert ratio == pytest.approx(1.0 - 50 / 150, abs=0.05)

    def test_no_data_regions(self):
        """Regions with no data get weight 0 (from splice adjustment)."""
        counts = _make_region_counts(
            n_unspliced_pos=[0],
            n_unspliced_neg=[0],
        )
        region_df = _make_region_df(1)
        stats = compute_region_stats(counts, region_df)

        weights = score_regions(stats, kappa_sym=50.0, strand_specificity=0.95)
        # n_unspliced < 2 → has_data=False → strand_weight=0.5
        # splice_rate=0 → weight = 0.5
        assert weights[0] == pytest.approx(0.5)

    def test_antisense_gene_scoring(self):
        """Antisense gene: RNA expects strand_ratio ~ (1-SS)."""
        # For - strand gene with SS=0.95, RNA expects strand_ratio ~ 0.05
        # Region with 5% on + strand → consistent with RNA
        counts = _make_region_counts(
            n_unspliced_pos=[5],
            n_unspliced_neg=[95],
        )
        region_df = _make_region_df(
            1,
            tx_pos=np.array([False]),
            tx_neg=np.array([True]),
        )
        stats = compute_region_stats(counts, region_df)

        weights = score_regions(stats, kappa_sym=50.0, strand_specificity=0.95)
        # Strand ratio 0.05 matches RNA expectation of 1-SS=0.05 → low gDNA
        assert weights[0] < 0.3


# ===================================================================
# TestBuildGDNAFLModel
# ===================================================================


class TestBuildGDNAFLModel:
    """Tests for build_gdna_fl_model."""

    def test_basic_fl_model(self):
        """FL observations weighted by region weights."""
        region_ids = np.array([0, 0, 1, 1, 2], dtype=np.int32)
        frag_lens = np.array([200, 250, 200, 300, 200], dtype=np.int32)
        weights = np.array([1.0, 0.5, 0.0])

        model = build_gdna_fl_model(region_ids, frag_lens, weights)

        # Region 0: weight 1.0 → FL 200 (w=1), 250 (w=1)
        # Region 1: weight 0.5 → FL 200 (w=0.5), 300 (w=0.5)
        # Region 2: weight 0.0 → excluded
        assert model.counts[200] == pytest.approx(1.5)
        assert model.counts[250] == pytest.approx(1.0)
        assert model.counts[300] == pytest.approx(0.5)

    def test_empty_fl_table(self):
        """Empty FL table → finalized but empty model."""
        model = build_gdna_fl_model(
            np.array([], dtype=np.int32),
            np.array([], dtype=np.int32),
            np.ones(5),
        )
        assert model._finalized
        assert model.total_weight == 0.0

    def test_zero_weights_excluded(self):
        """FL observations from zero-weight regions are excluded."""
        region_ids = np.array([0, 1], dtype=np.int32)
        frag_lens = np.array([200, 300], dtype=np.int32)
        weights = np.array([0.0, 0.0])

        model = build_gdna_fl_model(region_ids, frag_lens, weights)
        assert model.total_weight == 0.0

    def test_fl_out_of_range_excluded(self):
        """FL > max_fl or <= 0 excluded."""
        region_ids = np.array([0, 0, 0], dtype=np.int32)
        frag_lens = np.array([0, 500, 1500], dtype=np.int32)
        weights = np.array([1.0])

        model = build_gdna_fl_model(region_ids, frag_lens, weights, max_fl=1000)
        # Only fl=500 should be counted
        assert model.counts[500] == pytest.approx(1.0)
        assert model.counts[0] == 0.0


# ===================================================================
# TestCalibrateGDNA (integration)
# ===================================================================


class TestCalibrateGDNA:
    """Integration tests for calibrate_gdna."""

    def _make_synthetic_data(
        self,
        n_gdna: int = 100,
        n_rna: int = 100,
        kappa_true: float = 50.0,
        ss: float = 0.95,
        n_per: int = 100,
        rng_seed: int = 42,
    ):
        """Generate synthetic region data with known gDNA/RNA split.

        Returns (region_counts, fl_table, region_df, expected_kappa).
        """
        rng = np.random.default_rng(rng_seed)
        n_total = n_gdna + n_rna

        # gDNA regions: symmetric strand, no splicing, FL ~ 200±50
        alpha = kappa_true / 2.0
        p_gdna = rng.beta(alpha, alpha, size=n_gdna)
        n_pos_gdna = rng.binomial(n_per, p_gdna)
        n_neg_gdna = n_per - n_pos_gdna
        sp_pos_gdna = np.zeros(n_gdna, dtype=int)
        sp_neg_gdna = np.zeros(n_gdna, dtype=int)

        # RNA regions: biased strand (+ gene, SS=0.95), has splicing
        n_pos_rna = rng.binomial(n_per, ss, size=n_rna)
        n_neg_rna = n_per - n_pos_rna
        sp_pos_rna = rng.poisson(20, size=n_rna).astype(int)
        sp_neg_rna = rng.poisson(2, size=n_rna).astype(int)

        n_pos = np.concatenate([n_pos_gdna, n_pos_rna]).astype(np.float32)
        n_neg = np.concatenate([n_neg_gdna, n_neg_rna]).astype(np.float32)
        sp_pos = np.concatenate([sp_pos_gdna, sp_pos_rna]).astype(np.float32)
        sp_neg = np.concatenate([sp_neg_gdna, sp_neg_rna]).astype(np.float32)

        region_counts = _make_region_counts(n_pos, n_neg, sp_pos, sp_neg)

        # All regions are + strand gene exons
        region_df = _make_region_df(
            n_total,
            tx_pos=np.ones(n_total, dtype=bool),
            tx_neg=np.zeros(n_total, dtype=bool),
            lengths=np.full(n_total, 1000),
        )

        # FL observations: gDNA regions get FL ~ 200, RNA ~ 300
        fl_ids_g = np.repeat(np.arange(n_gdna), 5)
        fl_lens_g = rng.normal(200, 20, size=len(fl_ids_g)).astype(int)
        fl_lens_g = np.clip(fl_lens_g, 50, 500)

        fl_ids_r = np.repeat(np.arange(n_gdna, n_total), 5)
        fl_lens_r = rng.normal(300, 30, size=len(fl_ids_r)).astype(int)
        fl_lens_r = np.clip(fl_lens_r, 100, 600)

        fl_table = _make_fl_table(
            np.concatenate([fl_ids_g, fl_ids_r]),
            np.concatenate([fl_lens_g, fl_lens_r]),
        )

        return region_counts, fl_table, region_df, kappa_true

    def test_convergence(self):
        """Algorithm converges within max_iterations."""
        rc, fl, rdf, _ = self._make_synthetic_data()

        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert isinstance(result, GDNACalibration)
        assert result.converged
        assert result.n_iterations <= 20

    def test_kappa_recovery(self):
        """Recovered κ is in the right ballpark of true κ."""
        rc, fl, rdf, kappa_true = self._make_synthetic_data(kappa_true=50.0)

        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        # Within factor of 3
        assert kappa_true / 3 < result.kappa_sym < kappa_true * 3

    def test_weights_gdna_high_rna_low(self):
        """gDNA regions get higher weights than RNA regions."""
        rc, fl, rdf, _ = self._make_synthetic_data(n_gdna=100, n_rna=100)

        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)

        gdna_weights = result.region_weights[:100]
        rna_weights = result.region_weights[100:]

        # gDNA mean weight should be substantially higher
        assert gdna_weights.mean() > rna_weights.mean()
        # gDNA weights should mostly be > 0.5
        assert (gdna_weights > 0.5).sum() > 50
        # RNA weights should mostly be < 0.5
        assert (rna_weights < 0.5).sum() > 50

    def test_gdna_fl_model_has_data(self):
        """The gDNA FL model has observations."""
        rc, fl, rdf, _ = self._make_synthetic_data()

        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        assert result.gdna_fl_model.total_weight > 0
        assert result.gdna_fl_model._finalized

    def test_unstranded_convergence(self):
        """Algorithm works with unstranded data (SS=0.5)."""
        rng = np.random.default_rng(7)
        n_regions = 200
        n_per = 100

        # gDNA: symmetric, no splice, low density
        n_pos_g = rng.binomial(n_per, 0.5, size=100).astype(np.float32)
        n_neg_g = (n_per - n_pos_g).astype(np.float32)
        sp_g = np.zeros(100, dtype=np.float32)

        # RNA: symmetric (unstranded!), but has splicing and higher expression
        n_pos_r = rng.binomial(n_per * 5, 0.5, size=100).astype(np.float32)
        n_neg_r = (n_per * 5 - n_pos_r).astype(np.float32)
        sp_r = rng.poisson(30, size=100).astype(np.float32)

        rc = _make_region_counts(
            np.concatenate([n_pos_g, n_pos_r]),
            np.concatenate([n_neg_g, n_neg_r]),
            np.concatenate([sp_g, sp_r]),
        )
        rdf = _make_region_df(
            n_regions,
            lengths=np.full(n_regions, 1000),
        )
        fl = _make_fl_table([], [])

        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.5)
        assert isinstance(result, GDNACalibration)
        assert result.n_iterations >= 1

    def test_all_spliced_fallback(self):
        """If every region has spliced reads → fallback κ, empty seed."""
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50],
            n_unspliced_neg=[50, 50],
            n_spliced_pos=[10, 10],
        )
        rdf = _make_region_df(2)
        fl = _make_fl_table([], [])

        result = calibrate_gdna(counts, fl, rdf, strand_specificity=0.95)
        assert result.kappa_sym > 0  # should be fallback
        assert result.n_iterations >= 1

    def test_single_region(self):
        """Single region doesn't crash."""
        counts = _make_region_counts(
            n_unspliced_pos=[50],
            n_unspliced_neg=[50],
        )
        rdf = _make_region_df(1)
        fl = _make_fl_table([0], [200])

        result = calibrate_gdna(counts, fl, rdf, strand_specificity=0.95)
        assert isinstance(result, GDNACalibration)

    def test_custom_parameters(self):
        """Custom tuning parameters are respected."""
        rc, fl, rdf, _ = self._make_synthetic_data()

        result = calibrate_gdna(
            rc, fl, rdf,
            strand_specificity=0.95,
            min_coverage=5,
            max_iterations=3,
            convergence_tol=1e-6,
            kappa_min=10.0,
            kappa_max=500.0,
        )
        assert result.n_iterations <= 3
        assert 10.0 <= result.kappa_sym <= 500.0


# ===================================================================
# Edge cases
# ===================================================================


class TestEdgeCases:
    """Edge case tests for robustness."""

    def test_all_gdna(self):
        """Dataset with only gDNA (no spliced reads anywhere)."""
        rng = np.random.default_rng(42)
        n = 50
        n_per = 100
        p = rng.beta(25, 25, size=n)
        n_pos = rng.binomial(n_per, p).astype(np.float32)
        n_neg = (n_per - n_pos).astype(np.float32)

        rc = _make_region_counts(n_pos, n_neg)
        rdf = _make_region_df(n, tx_pos=np.ones(n, dtype=bool),
                              tx_neg=np.zeros(n, dtype=bool))
        fl = _make_fl_table(
            np.repeat(np.arange(n), 3),
            rng.integers(150, 300, size=n * 3),
        )

        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)
        # All regions should get high weight
        assert result.region_weights.mean() > 0.5

    def test_zero_length_regions(self):
        """Zero-length regions don't cause division by zero."""
        counts = _make_region_counts(
            n_unspliced_pos=[50, 10],
            n_unspliced_neg=[50, 10],
        )
        rdf = _make_region_df(2, lengths=np.array([0, 1000]))
        fl = _make_fl_table([], [])

        result = calibrate_gdna(counts, fl, rdf, strand_specificity=0.95)
        assert np.all(np.isfinite(result.region_weights))

    def test_reproducibility(self):
        """Same input → same output (deterministic)."""
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

        assert r1.kappa_sym == r2.kappa_sym
        assert r1.n_iterations == r2.n_iterations
        np.testing.assert_array_equal(r1.region_weights, r2.region_weights)


# ===================================================================
# TestUnifiedScoring — tests for combined strand + density + splice
# ===================================================================


class TestUnifiedScoring:
    """Tests that strand, density, and splice signals combine correctly."""

    def test_strand_contributes_zero_when_unstranded(self):
        """When SS=0.5, strand LLR is small — density dominates."""
        # Two regions with identical counts but different region lengths
        # (so different density)
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50],
            n_unspliced_neg=[50, 50],
        )
        # Intergenic (gene_strand=0)
        region_df = _make_region_df(
            2,
            tx_pos=np.array([False, False]),
            tx_neg=np.array([False, False]),
            lengths=np.array([10000, 100]),
        )
        stats = compute_region_stats(counts, region_df)

        # Without density: both regions have identical strand → equal weights
        w_no_density = score_regions(
            stats, kappa_sym=50.0, strand_specificity=0.5,
        )
        assert w_no_density[0] == pytest.approx(w_no_density[1], abs=0.01)

        # With density: R0 (100/10000=0.01) matches gDNA density,
        # R1 (100/100=1.0) has extreme excess
        w_with_density = score_regions(
            stats, kappa_sym=50.0, strand_specificity=0.5,
            gdna_density=0.01,
        )
        # Density differentiates: R0 matches background, R1 excess → RNA
        assert w_with_density[0] > w_with_density[1]

    def test_density_supplements_strand_for_stranded_data(self):
        """Density provides additional discrimination beyond strand."""
        # Two regions: same strand behavior (symmetric, gDNA-like)
        # but different density via different region lengths.
        counts = _make_region_counts(
            n_unspliced_pos=[50, 50],
            n_unspliced_neg=[50, 50],
        )
        region_df = _make_region_df(
            2,
            tx_pos=np.array([True, True]),
            tx_neg=np.array([False, False]),
            lengths=np.array([10000, 100]),  # very different lengths
        )
        stats = compute_region_stats(counts, region_df)

        # Without density: both have same strand → same weight
        w_no_density = score_regions(
            stats, kappa_sym=50.0, strand_specificity=0.70,
        )
        assert w_no_density[0] == pytest.approx(w_no_density[1], abs=0.01)

        # With density: R0 (100/10000=0.01) matches gDNA,
        # R1 (100/100=1.0) has extreme excess
        w_with_density = score_regions(
            stats, kappa_sym=50.0, strand_specificity=0.70,
            gdna_density=0.01,
        )
        # Density differentiates: R0 matches, R1 huge excess → RNA
        assert w_with_density[0] > w_with_density[1]

    def test_strand_and_density_reinforce(self):
        """Strand + density in agreement → stronger classification."""
        # Region 0: symmetric strand + low density → gDNA (both agree)
        # Region 1: biased strand + high density → RNA (both agree)
        counts = _make_region_counts(
            n_unspliced_pos=[50, 90],
            n_unspliced_neg=[50, 10],
        )
        region_df = _make_region_df(
            2,
            tx_pos=np.array([True, True]),
            tx_neg=np.array([False, False]),
            lengths=np.array([10000, 1000]),
        )
        stats = compute_region_stats(counts, region_df)

        w = score_regions(
            stats, kappa_sym=50.0, strand_specificity=0.95,
            gdna_density=0.01,  # Region 0: expected=100, actual=100 → match!
                                # Region 1: expected=10, actual=100 → excess!
        )
        assert w[0] > 0.8   # strong gDNA: strand + density agree
        assert w[1] < 0.1   # strong RNA: strand + density agree

    def test_strand_and_density_conflict(self):
        """When strand says gDNA but density says RNA, weight is reduced."""
        # Symmetric strand (gDNA-like) but high density (RNA-like).
        # Use moderate SS=0.70 so strand signal isn't overwhelming.
        counts = _make_region_counts(
            n_unspliced_pos=[50],
            n_unspliced_neg=[50],
        )
        region_df = _make_region_df(
            1,
            tx_pos=np.array([True]),
            tx_neg=np.array([False]),
            lengths=np.array([1000]),
        )
        stats = compute_region_stats(counts, region_df)

        # Strand only: symmetric → gDNA
        w_strand_only = score_regions(
            stats, kappa_sym=50.0, strand_specificity=0.70,
        )
        assert w_strand_only[0] > 0.9

        # With density: expected 10, observed 100 → strong excess
        w_combined = score_regions(
            stats, kappa_sym=50.0, strand_specificity=0.70,
            gdna_density=0.01,
        )
        # Density evidence reduces weight vs strand-only
        assert w_combined[0] < w_strand_only[0]

    def test_moderate_ss_uses_both_signals(self):
        """At intermediate SS (e.g., 0.75), both strand and density contribute."""
        # Region 0: gDNA-like on both signals
        # Region 1: RNA-like on both signals
        counts = _make_region_counts(
            n_unspliced_pos=[50, 75],
            n_unspliced_neg=[50, 25],
        )
        region_df = _make_region_df(
            2,
            tx_pos=np.array([True, True]),
            tx_neg=np.array([False, False]),
            lengths=np.array([10000, 1000]),
        )
        stats = compute_region_stats(counts, region_df)

        w = score_regions(
            stats, kappa_sym=50.0, strand_specificity=0.75,
            gdna_density=0.01,
        )
        # Region 0: symmetric strand (gDNA) + matches density → gDNA
        # Region 1: biased strand (RNA at SS=0.75) + excess density → RNA
        assert w[0] > w[1]

    def test_seed_applies_all_criteria(self):
        """Seed selection applies symmetry AND density for ALL SS values."""
        # Construct 4 regions:
        # R0: symmetric, zero splice, low density → seed ✓
        # R1: symmetric, zero splice, high density → fails density
        # R2: asymmetric, zero splice, low density → fails symmetry
        # R3: symmetric, has splice, low density → fails splice
        counts = _make_region_counts(
            n_unspliced_pos=[50, 250, 90, 50],
            n_unspliced_neg=[50, 250, 10, 50],
            n_spliced_pos=[0, 0, 0, 5],
        )
        region_df = _make_region_df(
            4,
            lengths=np.array([10000, 1000, 10000, 10000]),
        )
        stats = compute_region_stats(counts, region_df)

        # Test with stranded data
        mask_stranded = select_seed(
            stats, strand_specificity=0.95, coverage_percentile=50.0,
        )
        assert mask_stranded[0]    # passes all
        assert not mask_stranded[1]  # fails density
        assert not mask_stranded[2]  # fails symmetry
        assert not mask_stranded[3]  # fails splice

        # Same criteria for unstranded data — same results since
        # symmetry test still applies (symmetric regions pass regardless of SS)
        mask_unstranded = select_seed(
            stats, strand_specificity=0.5, coverage_percentile=50.0,
        )
        assert mask_unstranded[0]
        assert not mask_unstranded[1]
        assert not mask_unstranded[2]
        assert not mask_unstranded[3]

    def test_no_ss_threshold_boundary(self):
        """No discontinuity at any SS value — results vary smoothly."""
        counts = _make_region_counts(
            n_unspliced_pos=[80],
            n_unspliced_neg=[20],
        )
        region_df = _make_region_df(
            1,
            tx_pos=np.array([True]),
            tx_neg=np.array([False]),
        )
        stats = compute_region_stats(counts, region_df)

        # Score at densely-spaced SS values around old threshold (0.65)
        ss_values = np.arange(0.50, 0.96, 0.01)
        weights = np.array([
            score_regions(
                stats, kappa_sym=50.0, strand_specificity=ss,
                gdna_density=0.1,
            )[0]
            for ss in ss_values
        ])

        # Key property: no sudden jump between adjacent SS values.
        # A hard threshold would produce a near-1.0 jump; the smooth
        # unified model has bounded gradient.
        diffs = np.abs(np.diff(weights))
        assert np.all(diffs < 0.5), (
            f"Max jump = {diffs.max():.3f} at SS={ss_values[diffs.argmax()]:.2f}"
        )

    def test_all_three_signals_combined_integration(self):
        """Full integration: strand + density + splice all contribute."""
        rng = np.random.default_rng(42)
        n_regions = 300

        # Group A (100): pure gDNA — symmetric, no splice, background density
        # Group B (100): pure RNA — biased, spliced, high density
        # Group C (100): ambiguous — symmetric, some splice, moderate density
        n_per = 100

        # Group A
        p_a = rng.beta(25, 25, size=100)
        npa = rng.binomial(n_per, p_a).astype(np.float32)
        nna = (n_per - npa).astype(np.float32)
        spa = np.zeros(100, dtype=np.float32)

        # Group B
        npb = rng.binomial(n_per * 3, 0.95, size=100).astype(np.float32)
        nnb = (n_per * 3 - npb).astype(np.float32)
        spb = rng.poisson(30, size=100).astype(np.float32)

        # Group C
        p_c = rng.beta(25, 25, size=100)
        npc = rng.binomial(n_per, p_c).astype(np.float32)
        nnc = (n_per - npc).astype(np.float32)
        spc = rng.poisson(5, size=100).astype(np.float32)

        rc = _make_region_counts(
            np.concatenate([npa, npb, npc]),
            np.concatenate([nna, nnb, nnc]),
            np.concatenate([spa, spb, spc]),
        )
        rdf = _make_region_df(
            n_regions,
            tx_pos=np.ones(n_regions, dtype=bool),
            tx_neg=np.zeros(n_regions, dtype=bool),
            lengths=np.full(n_regions, 1000),
        )
        fl = _make_fl_table([], [])

        result = calibrate_gdna(rc, fl, rdf, strand_specificity=0.95)

        # Group A (gDNA) should have highest mean weight
        # Group B (RNA) should have lowest
        # Group C (ambiguous) should be in between
        wa = result.region_weights[:100].mean()
        wb = result.region_weights[100:200].mean()
        wc = result.region_weights[200:].mean()

        assert wa > wc > wb
