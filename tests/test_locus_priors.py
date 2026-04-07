"""Tests for compute_locus_priors() — per-locus Dirichlet prior computation.

Validates fraction-based prior splitting:
  - γ = Σ(E[gDNA]_r) / Σ(N_total_r) per locus
  - α_gDNA = γ × C,  α_RNA = (1−γ) × C
  - α_gDNA + α_RNA = C  always
  - Fallback when no regions overlap a locus
  - Budget C sensitivity
  - Pure RNA → α_gDNA = 0, α_RNA = C
  - Pure gDNA → α_gDNA ≈ C, α_RNA ≈ 0
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
) -> Locus:
    """Build a Locus with a single merged interval."""
    return Locus(
        locus_id=locus_id,
        transcript_indices=np.array([0], dtype=np.int32),
        unit_indices=np.array([0], dtype=np.int32),
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
# Basic fraction-based prior computation
# ======================================================================


class TestBasicPriors:
    """Verify γ = E[gDNA]/N and budget splitting α_gDNA + α_rna = C."""

    def test_single_locus_single_region(self):
        """One locus overlapping one region → exact fraction formula."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([100.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal)

        # γ = 100/1000 = 0.1, C = 1.0
        assert alpha_g[0] == pytest.approx(0.1)
        assert alpha_r[0] == pytest.approx(0.9)
        assert alpha_g[0] + alpha_r[0] == pytest.approx(1.0)

    def test_multiple_overlapping_regions(self):
        """Locus spanning two regions → sums E[gDNA] and N from both."""
        regions = [("chr1", 0, 500), ("chr1", 500, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([30.0, 70.0]),
            region_n_total=np.array([500.0, 500.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal)

        # γ = 100/1000 = 0.1
        assert alpha_g[0] == pytest.approx(0.1)
        assert alpha_r[0] == pytest.approx(0.9)

    def test_no_overlap_uses_fallback(self):
        """Locus on chr2 with regions only on chr1 → global γ fallback."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([200.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr2", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal)

        # Fallback γ = 200/1000 = 0.2
        assert alpha_g[0] == pytest.approx(0.2)
        assert alpha_r[0] == pytest.approx(0.8)

    def test_multiple_loci(self):
        """Two loci with different overlaps get different γ but same total."""
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

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal)

        # Locus 0: γ = 10/500 = 0.02
        assert alpha_g[0] == pytest.approx(0.02)
        # Locus 1: γ = 90/500 = 0.18
        assert alpha_g[1] == pytest.approx(0.18)
        assert alpha_g[1] > alpha_g[0]
        # α_gDNA + α_rna = C always
        assert alpha_g[0] + alpha_r[0] == pytest.approx(1.0)
        assert alpha_g[1] + alpha_r[1] == pytest.approx(1.0)


# ======================================================================
# Budget conservation: α_gDNA + α_RNA = C
# ======================================================================


class TestBudgetConservation:
    """α_gDNA + α_RNA must always equal the total pseudocount C."""

    def test_budget_sums_to_C_default(self):
        """With default C=1.0, budgets sum to 1 for all loci."""
        regions = [("chr1", 0, 1000), ("chr1", 2000, 3000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([10.0, 90.0]),
            region_n_total=np.array([500.0, 2000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [
            _make_locus(0, "chr1", 0, 1000),
            _make_locus(1, "chr1", 2000, 3000),
        ]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal)

        for i in range(len(loci)):
            assert alpha_g[i] + alpha_r[i] == pytest.approx(1.0)

    def test_budget_sums_to_C_custom(self):
        """With custom C, budgets sum to C."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([50.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(
            loci, index, cal, total_pseudocount=5.0,
        )

        assert alpha_g[0] + alpha_r[0] == pytest.approx(5.0)
        # γ = 0.05, so α_gDNA = 0.25, α_rna = 4.75
        assert alpha_g[0] == pytest.approx(0.25)
        assert alpha_r[0] == pytest.approx(4.75)

    def test_ratio_independent_of_locus_depth(self):
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

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal)

        # Both loci have γ = 0.1 despite 1000× depth difference
        assert alpha_g[0] / (alpha_g[0] + alpha_r[0]) == pytest.approx(0.1)
        assert alpha_g[1] / (alpha_g[1] + alpha_r[1]) == pytest.approx(0.1)


# ======================================================================
# Budget C sensitivity
# ======================================================================


class TestBudgetSensitivity:
    """C controls total magnitude; γ ratio is unaffected."""

    def test_doubling_C_doubles_both_alphas(self):
        """Doubling C doubles both α_gDNA and α_RNA."""
        regions = [("chr1", 0, 1000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([100.0]),
            region_n_total=np.array([1000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        ag1, ar1 = compute_locus_priors(loci, index, cal, total_pseudocount=1.0)
        ag2, ar2 = compute_locus_priors(loci, index, cal, total_pseudocount=2.0)

        assert ag2[0] == pytest.approx(2.0 * ag1[0])
        assert ar2[0] == pytest.approx(2.0 * ar1[0])
        # Ratio unchanged
        assert ag1[0] / ar1[0] == pytest.approx(ag2[0] / ar2[0])


# ======================================================================
# Pure RNA → α_gDNA = 0
# ======================================================================


class TestPureRNA:
    """In a pure RNA sample, E[gDNA] = 0 → γ = 0 → α_gDNA = 0."""

    def test_zero_e_gdna_gives_zero_alpha(self):
        """If calibration reports E[gDNA] = 0 everywhere, α_gDNA = 0."""
        regions = [("chr1", 0, 1000), ("chr1", 2000, 3000)]
        region_df = _make_region_df(regions)
        cal = _make_calibration(
            region_e_gdna=np.array([0.0, 0.0]),
            region_n_total=np.array([5000.0, 3000.0]),
        )
        index = _MockRegionIndex(region_df)
        loci = [
            _make_locus(0, "chr1", 0, 1000),
            _make_locus(1, "chr1", 2000, 3000),
        ]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal)

        assert alpha_g[0] == pytest.approx(0.0)
        assert alpha_g[1] == pytest.approx(0.0)
        # All budget goes to RNA
        assert alpha_r[0] == pytest.approx(1.0)
        assert alpha_r[1] == pytest.approx(1.0)


# ======================================================================
# Pure gDNA → α_gDNA dominates
# ======================================================================


class TestPureGDNA:
    """In a pure gDNA sample, E[gDNA] ≈ N_total → α_gDNA ≈ C."""

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

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal)

        # γ = 0.95, so α_gDNA = 0.95, α_RNA = 0.05
        assert alpha_g[0] == pytest.approx(0.95)
        assert alpha_r[0] == pytest.approx(0.05)
        assert alpha_g[0] > 10 * alpha_r[0]


# ======================================================================
# Fallback pathways
# ======================================================================


class TestFallback:
    """Tests for the no-cgranges and no-overlap fallback paths."""

    def test_no_region_cr_uses_global_fallback(self):
        """Index without region_cr → global γ fallback."""
        cal = _make_calibration(
            region_e_gdna=np.array([100.0, 200.0]),
            region_n_total=np.array([500.0, 500.0]),
        )
        index = _MockRegionIndex(region_df=None)
        loci = [_make_locus(0, "chr1", 0, 1000)]

        alpha_g, alpha_r = compute_locus_priors(loci, index, cal)

        # Fallback γ = 300/1000 = 0.3
        assert alpha_g[0] == pytest.approx(0.3)
        assert alpha_r[0] == pytest.approx(0.7)

    def test_empty_loci_returns_empty(self):
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
