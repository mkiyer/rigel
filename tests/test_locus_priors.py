"""Tests for compute_locus_priors() — per-locus Dirichlet prior computation.

Validates mode-aware Dirichlet prior splitting:
  - γ = Σ(E[gDNA]_r) / Σ(N_total_r) per locus
  - α_gDNA = γ × c_base
  - α_RNA  = (1−γ) × c_base
  - α_gDNA + α_RNA = c_base  always
  - Fallback when no regions overlap a locus
  - c_base scales linearly
  - Priors independent of N_locus (no κ×N)
  - Pure RNA → α_gDNA = 0
  - Pure gDNA → α_gDNA ≈ c_base
  - No-cgranges fallback pathway
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import pytest

from rigel.calibration import CalibrationResult
from rigel.locus import compute_locus_priors
from rigel.scored_fragments import Locus


# ---------------------------------------------------------------------------
# Minimal mock index with region_cr
# ---------------------------------------------------------------------------


@dataclass
class _MockRegionIndex:
    """Minimal TranscriptIndex substitute with region_cr for testing."""

    region_df: pd.DataFrame | None
    _region_cr: object = None

    @property
    def region_cr(self):
        if self._region_cr is not None:
            return self._region_cr
        if self.region_df is None or len(self.region_df) == 0:
            return None
        from rigel._cgranges_impl import cgranges

        cr = cgranges()
        for _, row in self.region_df.iterrows():
            cr.add(row["ref"], int(row["start"]), int(row["end"]), int(row["region_id"]))
        cr.index()
        self._region_cr = cr
        return cr


def _make_calibration(
    region_e_gdna: np.ndarray,
    region_n_total: np.ndarray,
    ss: float = 0.95,
) -> CalibrationResult:
    """Build a minimal CalibrationResult."""
    return CalibrationResult(
        region_e_gdna=np.asarray(region_e_gdna, dtype=np.float64),
        region_n_total=np.asarray(region_n_total, dtype=np.float64),
        gdna_fl_model=None,
        lambda_gdna=0.001,
        strand_specificity=ss,
    )


def _make_locus(
    locus_id: int,
    ref: str,
    start: int,
    end: int,
    n_units: int = 1,
) -> Locus:
    """Build a Locus with a single merged interval."""
    return Locus(
        locus_id=locus_id,
        transcript_indices=np.array([0], dtype=np.int32),
        unit_indices=np.arange(n_units, dtype=np.int32),
        gdna_span=end - start,
        merged_intervals=[(ref, start, end)],
    )


def _make_region_df(regions: list[tuple[str, int, int]]) -> pd.DataFrame:
    """Build a region DataFrame from (ref, start, end) tuples."""
    return pd.DataFrame(
        {
            "ref": [r[0] for r in regions],
            "start": [r[1] for r in regions],
            "end": [r[2] for r in regions],
            "region_id": list(range(len(regions))),
        }
    )


# ======================================================================
# Basic c_base prior computation
# ======================================================================


class TestCBaseComputation:
    """Verify α_gDNA = γ × c_base, α_RNA = (1−γ) × c_base."""

    def test_basic_split(self):
        """α_gDNA + α_RNA = c_base."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([100.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal, c_base=5.0)

        # γ = 100/1000 = 0.1
        assert alpha_g[0] == pytest.approx(0.5)   # 0.1 × 5.0
        assert alpha_r[0] == pytest.approx(4.5)   # 0.9 × 5.0
        assert alpha_g[0] + alpha_r[0] == pytest.approx(5.0)

    def test_pure_rna(self):
        """γ=0 → α_gDNA=0, α_RNA=c_base."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([0.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal, c_base=5.0)

        assert alpha_g[0] == pytest.approx(0.0)
        assert alpha_r[0] == pytest.approx(5.0)

    def test_pure_gdna(self):
        """γ=1 → α_gDNA=c_base, α_RNA=0."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([1000.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal, c_base=5.0)

        assert alpha_g[0] == pytest.approx(5.0)
        assert alpha_r[0] == pytest.approx(0.0)

    def test_c_base_scales_linearly(self):
        """c_base=10 gives 2× priors vs c_base=5."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([100.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        ag5, ar5 = compute_locus_priors(loci, index, cal, c_base=5.0)
        ag10, ar10 = compute_locus_priors(loci, index, cal, c_base=10.0)

        assert ag10[0] == pytest.approx(2.0 * ag5[0])
        assert ar10[0] == pytest.approx(2.0 * ar5[0])
        # Ratio unchanged
        assert ag5[0] / ar5[0] == pytest.approx(ag10[0] / ar10[0])

    def test_no_dependence_on_n_locus(self):
        """Priors independent of locus size (no κ×N)."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([100.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci_small = [_make_locus(0, "chr1", 0, 1000, n_units=10)]
        loci_large = [_make_locus(0, "chr1", 0, 1000, n_units=10000)]

        ag_s, ar_s = compute_locus_priors(loci_small, index, cal, c_base=5.0)
        ag_l, ar_l = compute_locus_priors(loci_large, index, cal, c_base=5.0)

        # Identical priors regardless of locus size
        assert ag_s[0] == pytest.approx(ag_l[0])
        assert ar_s[0] == pytest.approx(ar_l[0])

    def test_no_dependence_on_ss(self):
        """Priors independent of strand specificity (no κ interpolation)."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        results = []
        for ss in [0.50, 0.65, 0.90, 1.0]:
            cal = _make_calibration(
                region_e_gdna=np.array([100.0]),
                region_n_total=np.array([1000.0]),
                ss=ss,
            )
            ag, ar = compute_locus_priors(loci, index, cal, c_base=5.0)
            results.append((ag[0], ar[0]))

        # All SS values give identical priors
        for ag, ar in results:
            assert ag == pytest.approx(results[0][0])
            assert ar == pytest.approx(results[0][1])


# ======================================================================
# Region overlap
# ======================================================================


class TestRegionOverlap:
    """Verify per-locus γ from cgranges overlap."""

    def test_single_region(self):
        """One locus overlapping one region → exact fraction formula."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([100.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal, c_base=5.0)

        # γ = 100/1000 = 0.1
        assert alpha_g[0] == pytest.approx(0.5)
        assert alpha_r[0] == pytest.approx(4.5)

    def test_multiple_regions(self):
        """Locus spanning two regions → sums E[gDNA] and N from both."""
        regions = [("chr1", 0, 500), ("chr1", 500, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([30.0, 70.0]),
            region_n_total=np.array([500.0, 500.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal, c_base=5.0)

        # γ = 100/1000 = 0.1
        assert alpha_g[0] == pytest.approx(0.5)
        assert alpha_r[0] == pytest.approx(4.5)

    def test_no_overlap_fallback(self):
        """Locus on chr2 with regions only on chr1 → global γ fallback."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([200.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr2", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal, c_base=5.0)

        # Fallback γ = 200/1000 = 0.2
        assert alpha_g[0] == pytest.approx(1.0)
        assert alpha_r[0] == pytest.approx(4.0)

    def test_no_cgranges_fallback(self):
        """Index without region_cr → global γ fallback."""
        cal = _make_calibration(
            region_e_gdna=np.array([100.0, 200.0]),
            region_n_total=np.array([500.0, 500.0]),
        )
        index = _MockRegionIndex(region_df=None)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal, c_base=5.0)

        # Fallback γ = 300/1000 = 0.3
        assert alpha_g[0] == pytest.approx(1.5)
        assert alpha_r[0] == pytest.approx(3.5)

    def test_multiple_loci_different_gamma(self):
        """Two loci with different overlaps get different γ."""
        regions = [("chr1", 0, 1000), ("chr1", 2000, 3000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([10.0, 90.0]),
            region_n_total=np.array([500.0, 500.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [
            _make_locus(0, "chr1", 0, 1000),
            _make_locus(1, "chr1", 2000, 3000),
        ]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal, c_base=5.0)

        # Locus 0: γ = 10/500 = 0.02, α_gDNA = 0.1
        assert alpha_g[0] == pytest.approx(0.1)
        # Locus 1: γ = 90/500 = 0.18, α_gDNA = 0.9
        assert alpha_g[1] == pytest.approx(0.9)
        assert alpha_g[1] > alpha_g[0]
        # Both sum to c_base
        assert alpha_g[0] + alpha_r[0] == pytest.approx(5.0)
        assert alpha_g[1] + alpha_r[1] == pytest.approx(5.0)

    def test_ratio_independent_of_depth(self):
        """γ depends on fraction, not absolute count; ratio stable across depths."""
        regions = [("chr1", 0, 1000), ("chr1", 2000, 3000)]
        region_df = _make_region_df(regions)
        # Same γ (10%) but very different depths
        cal = _make_calibration(
            region_e_gdna=np.array([1.0, 1000.0]),
            region_n_total=np.array([10.0, 10000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [
            _make_locus(0, "chr1", 0, 1000),
            _make_locus(1, "chr1", 2000, 3000),
        ]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal, c_base=5.0)

        # Both loci have γ = 0.1 despite 1000× depth difference
        assert alpha_g[0] / (alpha_g[0] + alpha_r[0]) == pytest.approx(0.1)
        assert alpha_g[1] / (alpha_g[1] + alpha_r[1]) == pytest.approx(0.1)


# ======================================================================
# Edge cases
# ======================================================================


class TestEdgeCases:
    def test_empty_loci(self):
        """Zero loci → empty arrays."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([50.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)

        alpha_g, alpha_r = compute_locus_priors([], index, cal)

        assert len(alpha_g) == 0
        assert len(alpha_r) == 0

    def test_zero_total_fragments(self):
        """All regions with n_total=0 → fallback γ=0 (division guard)."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([0.0]),
            region_n_total=np.array([0.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal, c_base=5.0)

        # γ = 0/max(0,1) = 0
        assert alpha_g[0] == pytest.approx(0.0)
        assert alpha_r[0] == pytest.approx(5.0)

    def test_high_gdna_produces_large_alpha_ratio(self):
        """If nearly all fragments are gDNA, α_gDNA >> α_RNA."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([9500.0]),
            region_n_total=np.array([10000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal, c_base=5.0)

        # γ = 0.95
        assert alpha_g[0] == pytest.approx(4.75)
        assert alpha_r[0] == pytest.approx(0.25)
        assert alpha_g[0] > 10 * alpha_r[0]

    def test_default_c_base(self):
        """Default c_base=5.0 produces expected values."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([500.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal)

        # γ = 0.5, default c_base = 5.0
        assert alpha_g[0] == pytest.approx(2.5)
        assert alpha_r[0] == pytest.approx(2.5)
