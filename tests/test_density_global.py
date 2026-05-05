"""Tests for ``rigel.calibration.density_global.compute_global_densities``.

Pins the three locked formulas (INTERGENIC, INTRON, EXON-INTRON), the
overlap geometry $L_{\\mathrm{eff}} = L + \\bar L_{gDNA} - 1$, and the
schema/validation contract.

See ``docs/calibration/m4_implementation_plan.md`` §5.1.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._kappa import KAPPA_DEFAULT
from rigel.calibration.density_global import (
    GlobalDensityTable,
    GlobalGdnaDensity,
    compute_global_densities,
    l_eff_overlap,
)
from rigel.calibration.regions import RegionType
from rigel.calibration.scan_payload import (
    FL_HIST_N_BINS,
    MASK_N_STATES,
    CalibrationScanPayload,
)


_MASK_EXON = 0b001
_MASK_INTRON = 0b010
_MASK_INTERGENIC = 0b100


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------


def _make_region_df(rows: list[dict]) -> pd.DataFrame:
    """Build a region DataFrame from row dicts (with sensible defaults)."""
    df = pd.DataFrame(rows)
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    for col in ("start", "end"):
        df[col] = df[col].astype(np.int64)
    df["type"] = df["type"].astype(np.uint8)
    df.index = df["region_id"].to_numpy()
    return df


def _empty_payload(n_regions: int) -> dict:
    return {
        "global_counts": np.zeros(MASK_N_STATES, dtype=np.int64),
        "per_region_counts": np.zeros((n_regions, MASK_N_STATES), dtype=np.int64),
        "fl_hist": np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64),
        "u_left": np.zeros(n_regions, dtype=np.int64),
        "u_right": np.zeros(n_regions, dtype=np.int64),
        "n_observed": 0,
        "n_excluded_multimap": 0,
        "n_excluded_chimera": 0,
        "n_excluded_artifact": 0,
        "n_unobserved": 0,
        "n_unannotated_ref": 0,
    }


def _wrap_payload(d: dict) -> CalibrationScanPayload:
    n_obs = int(d["global_counts"].sum())
    d["n_observed"] = n_obs
    # Mirror n_obs into fl_hist[mask=1] bin 0 so the validator's
    # `sum(fl_hist) == n_observed` check passes — unused by M4.
    d["fl_hist"][_MASK_EXON, 0] = n_obs
    return CalibrationScanPayload.from_scan_dict(d, n_total=n_obs)


# ---------------------------------------------------------------------------
# Hand-counted scenarios
# ---------------------------------------------------------------------------


class TestHandCounted:
    def test_three_region_scenario(self):
        """One INTERGENIC, one INTRON, one EXON region with hand counts."""
        df = _make_region_df([
            {"ref_name": "chr1", "start": 0,    "end": 1000,
             "type": int(RegionType.INTERGENIC),
             "boundary_flux_left": False, "boundary_flux_right": False},
            {"ref_name": "chr1", "start": 1000, "end": 1500,
             "type": int(RegionType.EXON),
             "boundary_flux_left": False, "boundary_flux_right": True},
            {"ref_name": "chr1", "start": 1500, "end": 2500,
             "type": int(RegionType.INTRON),
             "boundary_flux_left": False, "boundary_flux_right": False},
        ])
        d = _empty_payload(3)
        d["per_region_counts"][0, _MASK_INTERGENIC] = 50    # INTERGENIC
        d["per_region_counts"][2, _MASK_INTRON] = 20         # INTRON
        d["u_right"][1] = 7                                  # right edge of EXON
        d["global_counts"][_MASK_INTERGENIC] = 50
        d["global_counts"][_MASK_INTRON] = 20
        d["global_counts"][_MASK_EXON | _MASK_INTRON] = 7    # boundary fragments

        payload = _wrap_payload(d)
        fl_mean = 200.0
        out = compute_global_densities(df, payload, gdna_fl_mean=fl_mean)

        # INTERGENIC: 50 frags / (1000 + 200 - 1) bp
        leff_ig = 1000.0 + fl_mean - 1.0
        assert out.intergenic.n_fragments == 50
        assert out.intergenic.eff_length_bp == pytest.approx(leff_ig)
        assert out.intergenic.rho == pytest.approx(50.0 / leff_ig)
        # INTRON: 20 / (1000 + 200 - 1)
        leff_in = 1000.0 + fl_mean - 1.0
        assert out.intron.n_fragments == 20
        assert out.intron.eff_length_bp == pytest.approx(leff_in)
        assert out.intron.rho == pytest.approx(20.0 / leff_in)
        # EXON-INTRON: u_R contributes only because bf_right is True.
        # Numerator = 7 * 1; denominator = (0 + 1) * 200 = 200.
        assert out.exon_intron.n_fragments == 7
        assert out.exon_intron.eff_length_bp == pytest.approx(200.0)
        assert out.exon_intron.rho == pytest.approx(7.0 / 200.0)
        assert out.exon_intron.n_regions_used == 1

    def test_pure_mrna_library(self):
        """All counts in EXON_ONLY mask → all three gDNA densities zero."""
        df = _make_region_df([
            {"ref_name": "chr1", "start": 0, "end": 500,
             "type": int(RegionType.INTERGENIC),
             "boundary_flux_left": False, "boundary_flux_right": False},
            {"ref_name": "chr1", "start": 500, "end": 1000,
             "type": int(RegionType.EXON),
             "boundary_flux_left": True, "boundary_flux_right": True},
            {"ref_name": "chr1", "start": 1000, "end": 2000,
             "type": int(RegionType.INTRON),
             "boundary_flux_left": False, "boundary_flux_right": False},
        ])
        d = _empty_payload(3)
        d["per_region_counts"][1, _MASK_EXON] = 100
        d["global_counts"][_MASK_EXON] = 100

        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl_mean=350.0)
        assert out.intergenic.rho == 0.0
        assert out.intron.rho == 0.0
        assert out.exon_intron.rho == 0.0
        # Empty numerators → kappa fallback.
        for d_ in (out.intergenic, out.intron, out.exon_intron):
            assert d_.kappa.value == KAPPA_DEFAULT
            assert d_.kappa.fallback_used

    def test_intron_only_library(self):
        df = _make_region_df([
            {"ref_name": "chr1", "start": 0, "end": 1000,
             "type": int(RegionType.INTRON),
             "boundary_flux_left": False, "boundary_flux_right": False},
            {"ref_name": "chr1", "start": 1000, "end": 2000,
             "type": int(RegionType.INTRON),
             "boundary_flux_left": False, "boundary_flux_right": False},
        ])
        d = _empty_payload(2)
        d["per_region_counts"][0, _MASK_INTRON] = 30
        d["per_region_counts"][1, _MASK_INTRON] = 70
        d["global_counts"][_MASK_INTRON] = 100
        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl_mean=300.0)
        assert out.intergenic.rho == 0.0
        assert out.intron.n_fragments == 100
        leff = (1000.0 + 300.0 - 1.0) * 2
        assert out.intron.rho == pytest.approx(100.0 / leff)
        assert out.exon_intron.rho == 0.0


class TestExonIntronEligibility:
    def test_all_flags_false_zero_density(self):
        df = _make_region_df([
            {"ref_name": "chr1", "start": 0, "end": 500,
             "type": int(RegionType.EXON),
             "boundary_flux_left": False, "boundary_flux_right": False},
        ])
        d = _empty_payload(1)
        d["u_left"][0] = 99
        d["u_right"][0] = 99
        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl_mean=200.0)
        assert out.exon_intron.rho == 0.0
        assert out.exon_intron.n_regions_used == 0
        assert out.exon_intron.eff_length_bp == 0.0

    def test_single_side_eligible(self):
        """bf_left=True, bf_right=False → u_R must NOT contribute."""
        df = _make_region_df([
            {"ref_name": "chr1", "start": 0, "end": 500,
             "type": int(RegionType.EXON),
             "boundary_flux_left": True, "boundary_flux_right": False},
        ])
        d = _empty_payload(1)
        d["u_left"][0] = 5
        d["u_right"][0] = 999     # ineligible side; must be ignored
        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl_mean=200.0)
        assert out.exon_intron.n_fragments == 5
        assert out.exon_intron.eff_length_bp == pytest.approx(200.0)
        assert out.exon_intron.rho == pytest.approx(5.0 / 200.0)


class TestEmptyPerType:
    def test_no_intron_regions(self):
        """Single-exon-genome: no INTRON rows at all → intron.rho == 0."""
        df = _make_region_df([
            {"ref_name": "chr1", "start": 0, "end": 1000,
             "type": int(RegionType.INTERGENIC),
             "boundary_flux_left": False, "boundary_flux_right": False},
            {"ref_name": "chr1", "start": 1000, "end": 1500,
             "type": int(RegionType.EXON),
             "boundary_flux_left": False, "boundary_flux_right": False},
        ])
        d = _empty_payload(2)
        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl_mean=200.0)
        assert out.intron.rho == 0.0
        assert out.intron.n_regions_used == 0
        assert out.intron.eff_length_bp == 0.0


class TestLeffGeometryLock:
    """The most important test in this file: pins the overlap formula."""

    def test_l_eff_overlap_helper(self):
        # 100 bp region, FL=350 → 100 + 350 - 1 = 449.
        assert l_eff_overlap(100.0, 350.0) == pytest.approx(449.0)
        # vector form
        spans = np.array([100.0, 200.0, 50.0])
        out = l_eff_overlap(spans, 350.0)
        np.testing.assert_allclose(out, np.array([449.0, 549.0, 399.0]))

    def test_pinned_through_orchestrator(self):
        df = _make_region_df([
            {"ref_name": "chr1", "start": 0, "end": 100,
             "type": int(RegionType.INTERGENIC),
             "boundary_flux_left": False, "boundary_flux_right": False},
        ])
        d = _empty_payload(1)
        d["per_region_counts"][0, _MASK_INTERGENIC] = 10
        d["global_counts"][_MASK_INTERGENIC] = 10
        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl_mean=350.0)
        assert out.intergenic.eff_length_bp == pytest.approx(449.0)
        assert out.intergenic.rho == pytest.approx(10.0 / 449.0)


# ---------------------------------------------------------------------------
# Schema / validation
# ---------------------------------------------------------------------------


class TestValidation:
    def test_region_count_mismatch(self):
        df = _make_region_df([
            {"ref_name": "chr1", "start": 0, "end": 100,
             "type": int(RegionType.INTERGENIC),
             "boundary_flux_left": False, "boundary_flux_right": False},
        ])
        d = _empty_payload(2)  # mismatched
        payload = _wrap_payload(d)
        with pytest.raises(ValueError, match="Rebuild the index"):
            compute_global_densities(df, payload, gdna_fl_mean=200.0)

    def test_gdna_fl_mean_must_be_positive(self):
        df = _make_region_df([
            {"ref_name": "chr1", "start": 0, "end": 100,
             "type": int(RegionType.INTERGENIC),
             "boundary_flux_left": False, "boundary_flux_right": False},
        ])
        d = _empty_payload(1)
        payload = _wrap_payload(d)
        with pytest.raises(ValueError, match="gdna_fl_mean"):
            compute_global_densities(df, payload, gdna_fl_mean=0.0)
        with pytest.raises(ValueError, match="gdna_fl_mean"):
            compute_global_densities(df, payload, gdna_fl_mean=-1.0)

    def test_to_summary_dict_round_trip(self):
        df = _make_region_df([
            {"ref_name": "chr1", "start": 0, "end": 100,
             "type": int(RegionType.INTERGENIC),
             "boundary_flux_left": False, "boundary_flux_right": False},
        ])
        d = _empty_payload(1)
        d["per_region_counts"][0, _MASK_INTERGENIC] = 5
        d["global_counts"][_MASK_INTERGENIC] = 5
        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl_mean=200.0)
        assert isinstance(out, GlobalDensityTable)
        sd = out.to_summary_dict()
        assert sd["gdna_fl_mean"] == 200.0
        for key in ("INTERGENIC", "INTRON", "EXON-INTRON"):
            block = sd[key]
            for col in (
                "rho", "n_fragments", "eff_length_bp", "n_regions_used",
                "kappa", "kappa_n_regions", "kappa_fallback", "kappa_reason",
            ):
                assert col in block, f"{key} missing {col}"

    def test_returns_globalgdnadensity_per_type(self):
        df = _make_region_df([
            {"ref_name": "chr1", "start": 0, "end": 100,
             "type": int(RegionType.INTERGENIC),
             "boundary_flux_left": False, "boundary_flux_right": False},
        ])
        d = _empty_payload(1)
        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl_mean=200.0)
        for d_ in (out.intergenic, out.intron, out.exon_intron):
            assert isinstance(d_, GlobalGdnaDensity)
        assert out.intergenic.type == "INTERGENIC"
        assert out.intron.type == "INTRON"
        assert out.exon_intron.type == "EXON-INTRON"
