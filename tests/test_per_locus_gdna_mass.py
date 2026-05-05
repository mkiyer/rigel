"""Tests for ``estimate_locus_gdna`` (per-Locus gDNA mass)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration._kappa import KappaEstimate
from rigel.calibration._region_index_py import RegionIndexPy
from rigel.calibration.density_global import GlobalDensityTable, GlobalGdnaDensity
from rigel.calibration.locus_prior import (
    FLAG_EXON_INTRON_NO_ELIGIBLE,
    FLAG_PI_CLIPPED,
    estimate_locus_gdna,
)
from rigel.calibration.regions import RegionType
from rigel.calibration.scan_payload import (
    FL_HIST_N_BINS,
    MASK_N_STATES,
    CalibrationScanPayload,
)
from rigel.locus import Locus


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

def _kappa(value: float = 100.0) -> KappaEstimate:
    return KappaEstimate(value=value, n_regions=10, fallback_used=False, fallback_reason="")


def _gdt(
    rho_ig: float = 0.0,
    rho_in: float = 0.0,
    rho_ei: float = 0.0,
    fl_mean: float = 200.0,
) -> GlobalDensityTable:
    return GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC", rho=rho_ig, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=0, kappa=_kappa(),
        ),
        intron=GlobalGdnaDensity(
            type="INTRON", rho=rho_in, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=0, kappa=_kappa(),
        ),
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON", rho=rho_ei, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=0, kappa=_kappa(),
        ),
        gdna_fl_mean=fl_mean,
    )


def _build_arrays(
    regions: list[tuple[str, int, int, int, bool, bool]],
    counts_intergenic: list[int] | None = None,
    counts_intron: list[int] | None = None,
    u_left: list[int] | None = None,
    u_right: list[int] | None = None,
):
    """regions: list of (ref, start, end, type, bf_left, bf_right)."""
    n = len(regions)
    df = pd.DataFrame(
        {
            "region_id": np.arange(n, dtype=np.int64),
            "ref_name": [r[0] for r in regions],
            "start": np.array([r[1] for r in regions], dtype=np.int64),
            "end": np.array([r[2] for r in regions], dtype=np.int64),
            "type": np.array([r[3] for r in regions], dtype=np.uint8),
            "strand": np.zeros(n, dtype=np.uint8),
            "tx_pos_bp": np.zeros(n, dtype=np.int64),
            "tx_neg_bp": np.zeros(n, dtype=np.int64),
            "exon_pos_bp": np.zeros(n, dtype=np.int64),
            "exon_neg_bp": np.zeros(n, dtype=np.int64),
            "boundary_flux_left": np.array([r[4] for r in regions], dtype=bool),
            "boundary_flux_right": np.array([r[5] for r in regions], dtype=bool),
        }
    )
    refs = sorted({r[0] for r in regions})
    ref_name_to_id = {name: i for i, name in enumerate(refs)}
    region_arrays = RegionArrays.from_region_df(df, ref_name_to_id)

    per_region = np.zeros((n, MASK_N_STATES), dtype=np.int64)
    if counts_intergenic is not None:
        per_region[:, 0b100] = np.array(counts_intergenic, dtype=np.int64)
    if counts_intron is not None:
        per_region[:, 0b010] = np.array(counts_intron, dtype=np.int64)
    payload = CalibrationScanPayload(
        global_counts=np.zeros(MASK_N_STATES, dtype=np.int64),
        per_region_counts=per_region,
        fl_hist=np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64),
        u_left=np.array(u_left if u_left is not None else [0] * n, dtype=np.int64),
        u_right=np.array(u_right if u_right is not None else [0] * n, dtype=np.int64),
        n_observed=0,
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_unobserved=0,
        n_unannotated_ref=0,
    )
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    region_index = RegionIndexPy(arrays=region_arrays)
    return region_arrays, payload_arrays, region_index, ref_name_to_id


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_single_intergenic_locus_exact_count():
    # One intergenic region [0, 1000); 50 fragments; FL mean 200.
    # κ = 0 ⇒ ρ_loco = 50 / (1000 + 199) = exact local density.
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        counts_intergenic=[50],
    )
    gdt = _gdt(rho_ig=0.0, fl_mean=200.0)
    # Force κ → 0 via a bespoke kappa
    gdt = GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC", rho=0.0, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=0, kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        intron=gdt.intron, exon_intron=gdt.exon_intron, gdna_fl_mean=200.0,
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    est = estimate_locus_gdna(
        locus=locus, n_obs=50, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl_mean=200.0,
    )
    leff_full = 1000 + 199.0
    assert est.n_gdna_intergenic == pytest.approx(50.0)
    assert est.n_gdna_intron == 0.0
    assert est.n_gdna_exon_intron == 0.0
    assert est.n_gdna == pytest.approx(50.0)
    assert est.leff_loco[0] == pytest.approx(leff_full)
    assert est.pi_gdna == pytest.approx(1.0)


def test_intron_zero_count_shrinks_to_global():
    # One intron region [0, 1000) with 0 counts; ρ_global = 0.05; κ = 100.
    # ρ_loco = (0 + 100·0.05) / (1199 + 100) = 5 / 1299
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 1000, int(RegionType.INTRON), False, False)],
        counts_intron=[0],
    )
    gdt = GlobalDensityTable(
        intergenic=_gdt().intergenic,
        intron=GlobalGdnaDensity(
            type="INTRON", rho=0.05, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=10, kappa=_kappa(100.0),
        ),
        exon_intron=_gdt().exon_intron, gdna_fl_mean=200.0,
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    est = estimate_locus_gdna(
        locus=locus, n_obs=10, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl_mean=200.0,
    )
    expected_rho = (0 + 100 * 0.05) / (1199.0 + 100.0)
    assert est.rho_loco[1] == pytest.approx(expected_rho)
    assert est.n_gdna_intron == pytest.approx(expected_rho * 1199.0)


def test_exon_only_no_eligible_uses_global():
    # One EXON region [0, 500) with NO eligible boundaries; ρ_global = 0.03.
    # n_eligible == 0 ⇒ ρ_loco = ρ_global, n_gdna = ρ · L_eff_full.
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 500, int(RegionType.EXON), False, False)],
        u_left=[0], u_right=[0],
    )
    gdt = GlobalDensityTable(
        intergenic=_gdt().intergenic,
        intron=_gdt().intron,
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON", rho=0.03, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=10, kappa=_kappa(50.0),
        ),
        gdna_fl_mean=200.0,
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=500)
    est = estimate_locus_gdna(
        locus=locus, n_obs=100, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl_mean=200.0,
    )
    leff_full = 500 + 199.0
    assert est.n_eligible_boundaries == 0
    assert est.fallback_flags & FLAG_EXON_INTRON_NO_ELIGIBLE
    assert est.rho_loco[2] == pytest.approx(0.03)
    assert est.n_gdna_exon_intron == pytest.approx(0.03 * leff_full)


def test_exon_with_both_eligible_boundaries_exact():
    # EXON region [0, 500) with both boundaries eligible; u_L=10, u_R=15.
    # κ = 0 ⇒ ρ_loco = (10 + 15) / (2 · 200) = 25 / 400 = 0.0625
    # n_gdna = 0.0625 · L_eff_full
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 500, int(RegionType.EXON), True, True)],
        u_left=[10], u_right=[15],
    )
    gdt = GlobalDensityTable(
        intergenic=_gdt().intergenic,
        intron=_gdt().intron,
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON", rho=0.0, n_fragments=25, eff_length_bp=400.0,
            n_regions_used=1,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        gdna_fl_mean=200.0,
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=500)
    est = estimate_locus_gdna(
        locus=locus, n_obs=100, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl_mean=200.0,
    )
    leff_full = 500 + 199.0
    assert est.n_eligible_boundaries == 2
    assert est.rho_loco[2] == pytest.approx(25.0 / 400.0)
    assert est.n_gdna_exon_intron == pytest.approx(25.0 / 400.0 * leff_full)


def test_short_locus_overlap_leff_regression():
    """Locus shorter than gdna_fl_mean must produce L_eff > 0.

    Master regression for the L_eff overlap formula
    (`/memories/repo/calibration-overlap-leff-2026-04.md`).
    """
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 50, int(RegionType.INTERGENIC), False, False)],
        counts_intergenic=[5],
    )
    gdt = _gdt(rho_ig=0.05, fl_mean=200.0)
    locus = Locus(ref="chr1", ref_id=0, start=0, end=50)
    est = estimate_locus_gdna(
        locus=locus, n_obs=10, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl_mean=200.0,
    )
    # Containment formula would give L_eff = max(0, 50 - 200 + 1) = 0 ⇒ n_gdna = 0.
    # Overlap formula gives L_eff = 50 + 199 = 249 ⇒ n_gdna > 0.
    expected_leff = 50 + 199.0
    assert est.leff_loco[0] == pytest.approx(expected_leff)
    assert est.n_gdna > 0.0


def test_all_three_region_types_summed():
    ra, pa, idx, _ = _build_arrays(
        regions=[
            ("chr1", 0, 100, int(RegionType.INTERGENIC), False, False),
            ("chr1", 100, 600, int(RegionType.EXON), True, False),
            ("chr1", 600, 1000, int(RegionType.INTRON), False, False),
        ],
        counts_intergenic=[5, 0, 0],
        counts_intron=[0, 0, 10],
        u_left=[0, 7, 0], u_right=[0, 0, 0],
    )
    gdt = GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC", rho=0.0, n_fragments=5, eff_length_bp=299.0,
            n_regions_used=1,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        intron=GlobalGdnaDensity(
            type="INTRON", rho=0.0, n_fragments=10, eff_length_bp=599.0,
            n_regions_used=1,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON", rho=0.0, n_fragments=7, eff_length_bp=200.0,
            n_regions_used=1,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        gdna_fl_mean=200.0,
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    est = estimate_locus_gdna(
        locus=locus, n_obs=22, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl_mean=200.0,
    )
    # κ = 0 ⇒ pure-local densities reproduce input fragments exactly:
    #   intergenic: n_gdna = 5
    #   intron:     n_gdna = 10
    #   exon-intron: ρ = 7/(1·200), L_eff_full = 500 + 199 = 699 ⇒ 7/200 · 699
    assert est.n_gdna_intergenic == pytest.approx(5.0)
    assert est.n_gdna_intron == pytest.approx(10.0)
    assert est.n_gdna_exon_intron == pytest.approx(7.0 / 200.0 * (500 + 199.0))
    assert est.n_gdna == pytest.approx(
        5.0 + 10.0 + 7.0 / 200.0 * 699.0
    )


def test_pi_clipped_when_overestimate():
    # Force massive intergenic over-prediction → pi > 1.
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        counts_intergenic=[1000],
    )
    gdt = GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC", rho=0.0, n_fragments=1000, eff_length_bp=1199.0,
            n_regions_used=1,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        intron=_gdt().intron, exon_intron=_gdt().exon_intron, gdna_fl_mean=200.0,
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    est = estimate_locus_gdna(
        locus=locus, n_obs=50,  # only 50 observed but predicted = 1000
        region_index=idx, region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl_mean=200.0,
    )
    assert est.pi_gdna == pytest.approx(1.0)
    assert est.fallback_flags & FLAG_PI_CLIPPED


def test_no_region_overlap_raises():
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        counts_intergenic=[5],
    )
    gdt = _gdt(rho_ig=0.05, fl_mean=200.0)
    # Locus on a different reference id (out of range) ⇒ empty overlap ⇒ raise.
    locus = Locus(ref="chrX", ref_id=99, start=0, end=1000)
    with pytest.raises(RuntimeError, match="overlaps no regions"):
        estimate_locus_gdna(
            locus=locus, n_obs=10, region_index=idx,
            region_arrays=ra, payload_arrays=pa,
            global_densities=gdt, gdna_fl_mean=200.0,
        )
