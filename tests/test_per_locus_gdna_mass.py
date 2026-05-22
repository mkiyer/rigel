"""Tests for ``estimate_locus_gdna`` (per-Locus gDNA mass)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration._exposure import boundary_crossing_exposure
from rigel.calibration._kappa import KappaEstimate
from rigel.calibration._region_index_py import RegionIndexPy
from rigel.calibration._orient import ORIENT_OPP, ORIENT_SAME, StrandSummary
from rigel.calibration.density_global import (
    ExonCompositeDensity,
    GlobalDensityTable,
    GlobalGdnaDensity,
)
from rigel.calibration.locus_prior import (
    FLAG_EXON_INTRON_NO_ELIGIBLE,
    FLAG_PI_CLIPPED,
    estimate_locus_gdna,
    expected_gdna_count_global,
)
from rigel.calibration.regions import RegionStrand, RegionType
from rigel.calibration.scan_payload import (
    FL_HIST_N_BINS,
    MASK_N_STATES,
    CalibrationScanPayload,
)
from rigel.frag_length_model import FragmentLengthModel
from rigel.locus import Locus


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------


def _delta_fl(length: int, *, max_size: int = 1024) -> FragmentLengthModel:
    """FL model peaked at ``length`` (Laplace smoothing remains, so use
    :func:`_oracle_leff` for the exact expected denominator)."""
    counts = np.zeros(max_size + 1, dtype=np.float64)
    counts[length] = 10_000.0
    return FragmentLengthModel.from_counts(counts, max_size=max_size)


def _oracle_leff(spans, gdna_fl: FragmentLengthModel) -> np.ndarray:
    """Ground-truth FL-PMF-weighted containment L_eff (no salmon floor)."""
    return gdna_fl.compute_all_transcript_eff_lens(np.asarray(spans, dtype=np.int64), min_value=0.0)


def _kappa(value: float = 100.0) -> KappaEstimate:
    return KappaEstimate(value=value, n_regions=10, fallback_used=False, fallback_reason="")


def _gdt(
    rho_ig: float = 0.0,
    rho_in: float = 0.0,
    rho_ei: float = 0.0,
    fl_mean: int = 200,
) -> GlobalDensityTable:
    return GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC",
            rho=rho_ig,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=0,
            kappa=_kappa(),
        ),
        intron=GlobalGdnaDensity(
            type="INTRON",
            rho=rho_in,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=0,
            kappa=_kappa(),
        ),
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON",
            rho=rho_ei,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=0,
            kappa=_kappa(),
        ),
        gdna_fl=_delta_fl(fl_mean),
    )


def _build_arrays(
    regions: list[tuple[str, int, int, int, bool, bool]],
    counts_intergenic: list[int] | None = None,
    counts_intron: list[int] | None = None,
    u_left: list[int] | None = None,
    u_right: list[int] | None = None,
    strands: list[int] | None = None,
    exon_contained_same: list[int] | None = None,
    exon_contained_opp: list[int] | None = None,
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
            "strand": np.array(
                strands if strands is not None else [int(RegionStrand.NONE)] * n,
                dtype=np.uint8,
            ),
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
    intron_by_orient = np.zeros((n, 3), dtype=np.int64)
    intron_by_orient[:, 2] = per_region[:, 0b010]
    u_left_arr = np.array(u_left if u_left is not None else [0] * n, dtype=np.int64)
    u_right_arr = np.array(u_right if u_right is not None else [0] * n, dtype=np.int64)
    u_left_by_orient = np.zeros((n, 3), dtype=np.int64)
    u_right_by_orient = np.zeros((n, 3), dtype=np.int64)
    u_left_by_orient[:, 2] = u_left_arr
    u_right_by_orient[:, 2] = u_right_arr
    exon_contained_by_orient = np.zeros((n, 3), dtype=np.int64)
    if exon_contained_same is not None:
        exon_contained_by_orient[:, ORIENT_SAME] = np.array(
            exon_contained_same,
            dtype=np.int64,
        )
    if exon_contained_opp is not None:
        exon_contained_by_orient[:, ORIENT_OPP] = np.array(
            exon_contained_opp,
            dtype=np.int64,
        )
    payload = CalibrationScanPayload(
        global_counts=np.zeros(MASK_N_STATES, dtype=np.int64),
        per_region_counts=per_region,
        fl_hist=np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64),
        u_left=u_left_arr,
        u_right=u_right_arr,
        intron_counts_by_orient=intron_by_orient,
        exon_contained_counts_by_orient=exon_contained_by_orient,
        u_left_by_orient=u_left_by_orient,
        u_right_by_orient=u_right_by_orient,
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
    # κ = 0 ⇒ ρ_loco = 50 / L_eff_contained = exact local density.
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        counts_intergenic=[50],
    )
    gdna_fl = _delta_fl(200)
    gdt = _gdt(rho_ig=0.0, fl_mean=200)
    # Force κ → 0 via a bespoke kappa
    gdt = GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC",
            rho=0.0,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=0,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        intron=gdt.intron,
        exon_intron=gdt.exon_intron,
        gdna_fl=gdna_fl,
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    est = estimate_locus_gdna(
        locus=locus,
        n_obs=50,
        region_index=idx,
        region_arrays=ra,
        payload_arrays=pa,
        global_densities=gdt,
        gdna_fl=gdna_fl,
    )
    leff_full = float(_oracle_leff([1000], gdna_fl)[0])
    assert est.n_gdna_intergenic == pytest.approx(50.0)
    assert est.n_gdna_intron == 0.0
    assert est.n_gdna_boundary_observed == 0.0
    assert est.n_gdna_exon_only == 0.0
    assert est.n_gdna == pytest.approx(50.0)
    assert est.leff_loco[0] == pytest.approx(leff_full)
    assert est.pi_gdna == pytest.approx(1.0)


def test_intron_zero_count_shrinks_to_global():
    # One intron region [0, 1000) with 0 counts; ρ_global = 0.05; κ = 100.
    # ρ_loco = (0 + κ·ρ_global) / (L_eff_contained + κ)
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 1000, int(RegionType.INTRON), False, False)],
        counts_intron=[0],
    )
    gdna_fl = _delta_fl(200)
    gdt = GlobalDensityTable(
        intergenic=_gdt().intergenic,
        intron=GlobalGdnaDensity(
            type="INTRON",
            rho=0.05,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=10,
            kappa=_kappa(100.0),
        ),
        exon_intron=_gdt().exon_intron,
        gdna_fl=gdna_fl,
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    est = estimate_locus_gdna(
        locus=locus,
        n_obs=10,
        region_index=idx,
        region_arrays=ra,
        payload_arrays=pa,
        global_densities=gdt,
        gdna_fl=gdna_fl,
    )
    leff = float(_oracle_leff([1000], gdna_fl)[0])
    expected_rho = (0 + 100 * 0.05) / (leff + 100.0)
    assert est.rho_loco[1] == pytest.approx(expected_rho)
    assert est.n_gdna_intron == pytest.approx(expected_rho * leff)


def test_exon_only_no_eligible_uses_global():
    # One EXON region [0, 500) with NO eligible boundaries; ρ_global = 0.03.
    # n_eligible == 0 ⇒ ρ_loco = ρ_global, n_exon_only = ρ · L_eff_full.
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 500, int(RegionType.EXON), False, False)],
        u_left=[0],
        u_right=[0],
    )
    gdna_fl = _delta_fl(200)
    gdt = GlobalDensityTable(
        intergenic=_gdt().intergenic,
        intron=_gdt().intron,
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON",
            rho=0.03,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=10,
            kappa=_kappa(50.0),
        ),
        gdna_fl=gdna_fl,
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=500)
    est = estimate_locus_gdna(
        locus=locus,
        n_obs=100,
        region_index=idx,
        region_arrays=ra,
        payload_arrays=pa,
        global_densities=gdt,
        gdna_fl=gdna_fl,
    )
    leff_full = float(_oracle_leff([500], gdna_fl)[0])
    assert est.n_eligible_boundaries == 0
    assert est.fallback_flags & FLAG_EXON_INTRON_NO_ELIGIBLE
    assert est.rho_loco[2] == pytest.approx(0.03)
    assert est.n_gdna_boundary_observed == 0.0
    assert est.n_gdna_exon_only == pytest.approx(0.03 * leff_full)


def test_exon_with_both_eligible_boundaries_exact():
    # EXON region [0, 500) with both boundaries eligible inside the locus;
    # u_L=10, u_R=15. κ=0 ⇒ ρ_loco = 25 / (2 · B_cross).
    # Decomposition: n_gdna_boundary_observed = 25 (raw events),
    # n_gdna_exon_only = ρ_loco · L_eff_contained(500).
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 500, int(RegionType.EXON), True, True)],
        u_left=[10],
        u_right=[15],
    )
    gdna_fl = _delta_fl(200)
    b_cross = boundary_crossing_exposure(gdna_fl)
    gdt = GlobalDensityTable(
        intergenic=_gdt().intergenic,
        intron=_gdt().intron,
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON",
            rho=0.0,
            n_fragments=25,
            eff_length_bp=2.0 * b_cross,
            n_regions_used=1,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        gdna_fl=gdna_fl,
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=500)
    est = estimate_locus_gdna(
        locus=locus,
        n_obs=100,
        region_index=idx,
        region_arrays=ra,
        payload_arrays=pa,
        global_densities=gdt,
        gdna_fl=gdna_fl,
    )
    leff_full = float(_oracle_leff([500], gdna_fl)[0])
    expected_rho = 25.0 / (2.0 * b_cross)
    assert est.n_eligible_boundaries == 2
    assert est.rho_loco[2] == pytest.approx(expected_rho)
    assert est.n_gdna_boundary_observed == pytest.approx(25.0)
    assert est.n_gdna_exon_only == pytest.approx(expected_rho * leff_full)
    # Decomposition invariant.
    assert est.n_gdna == pytest.approx(
        est.n_gdna_intergenic
        + est.n_gdna_intron
        + est.n_gdna_boundary_observed
        + est.n_gdna_exon_only
    )


def test_exon_contained_only_strand_evidence_can_replace_boundary_fallback():
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 500, int(RegionType.EXON), False, False)],
        strands=[int(RegionStrand.POS)],
        exon_contained_same=[50],
        exon_contained_opp=[50],
    )
    gdna_fl = _delta_fl(200)
    gdt = GlobalDensityTable(
        intergenic=_gdt().intergenic,
        intron=_gdt().intron,
        exon_intron=_gdt().exon_intron,
        gdna_fl=gdna_fl,
        exon_contained=GlobalGdnaDensity(
            type="EXON-CONTAINED",
            rho=0.0,
            n_fragments=100,
            eff_length_bp=0.0,
            n_regions_used=1,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=500)

    est = estimate_locus_gdna(
        locus=locus,
        n_obs=100,
        region_index=idx,
        region_arrays=ra,
        payload_arrays=pa,
        global_densities=gdt,
        gdna_fl=gdna_fl,
        strand_summary=StrandSummary(p_r1_sense=1.0, n_observations=1000),
    )

    leff_full = float(_oracle_leff([500], gdna_fl)[0])
    expected_rho = 100.0 / leff_full
    assert est.n_eligible_boundaries == 0
    assert est.fallback_flags & FLAG_EXON_INTRON_NO_ELIGIBLE
    assert est.n_gdna_boundary_observed == 0.0
    assert est.rho_loco[2] == pytest.approx(expected_rho)
    assert est.n_gdna_exon_only == pytest.approx(100.0)
    assert est.exon_gdna.rho_exon_boundary == pytest.approx(0.0)
    assert est.exon_gdna.rho_exon_contained == pytest.approx(expected_rho)
    assert est.exon_gdna.n_exon_contained_observed == pytest.approx(100.0)


def test_global_prior_keeps_boundary_density_but_uses_composite_for_contained_exon():
    ra, _, idx, _ = _build_arrays(
        regions=[("chr1", 0, 500, int(RegionType.EXON), True, True)],
    )
    gdna_fl = _delta_fl(200)
    b_cross = boundary_crossing_exposure(gdna_fl)
    gdt = GlobalDensityTable(
        intergenic=_gdt().intergenic,
        intron=_gdt().intron,
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON",
            rho=0.01,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=1,
            kappa=_kappa(),
        ),
        gdna_fl=gdna_fl,
        exon_contained=GlobalGdnaDensity(
            type="EXON-CONTAINED",
            rho=0.1,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=1,
            kappa=_kappa(),
            strand_active=True,
        ),
        exon_composite=ExonCompositeDensity(
            rho=0.1,
            boundary_rho=0.01,
            contained_rho=0.1,
            boundary_precision=0.0,
            contained_precision=1.0,
            strand_power=1.0,
            contained_active=True,
        ),
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=500)

    parts = expected_gdna_count_global(
        locus=locus,
        region_index=idx,
        region_arrays=ra,
        global_densities=gdt,
        gdna_fl=gdna_fl,
        b_cross=b_cross,
    )

    leff_full = float(_oracle_leff([500], gdna_fl)[0])
    assert parts.boundary_crossing_expected == pytest.approx(0.01 * 2.0 * b_cross)
    assert parts.exon_contained_expected == pytest.approx(0.1 * leff_full)
    assert parts.total == pytest.approx(
        parts.boundary_crossing_expected + parts.exon_contained_expected
    )


def test_short_locus_contained_leff_semantics():
    """Short locus + delta-FL: containment correctly gives L_eff = 0
    and n_gdna = 0 (a 200-bp fragment cannot fit in a 50-bp region).

    Realistic FL distributions have mass at small lengths so L_eff > 0
    in practice; this test pins the mathematical contract that the
    *contained* numerator and denominator measure the same event.

    Supersedes the pre-R1 test that expected the (incorrect) overlap
    formula L_eff = span + fl_mean - 1; see
    ``docs/calibration/density_eff_length_fix_2026_05.md``.
    """
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 50, int(RegionType.INTERGENIC), False, False)],
        counts_intergenic=[5],
    )
    gdna_fl = _delta_fl(200)
    gdt = _gdt(rho_ig=0.05, fl_mean=200)
    locus = Locus(ref="chr1", ref_id=0, start=0, end=50)
    est = estimate_locus_gdna(
        locus=locus,
        n_obs=10,
        region_index=idx,
        region_arrays=ra,
        payload_arrays=pa,
        global_densities=gdt,
        gdna_fl=gdna_fl,
    )
    expected_leff = float(_oracle_leff([50], gdna_fl)[0])
    assert est.leff_loco[0] == pytest.approx(expected_leff)
    # With smoothing the leff is small but nonzero; n_gdna may be small
    # but the locoregional shrinkage falls back to the global rho when
    # leff is tiny relative to kappa.
    expected_rho = (5 + 100 * 0.05) / (expected_leff + 100)
    assert est.rho_loco[0] == pytest.approx(expected_rho)


def test_all_three_region_types_summed():
    ra, pa, idx, _ = _build_arrays(
        regions=[
            ("chr1", 0, 100, int(RegionType.INTERGENIC), False, False),
            ("chr1", 100, 600, int(RegionType.EXON), True, False),
            ("chr1", 600, 1000, int(RegionType.INTRON), False, False),
        ],
        counts_intergenic=[5, 0, 0],
        counts_intron=[0, 0, 10],
        u_left=[0, 7, 0],
        u_right=[0, 0, 0],
    )
    gdna_fl = _delta_fl(200)
    b_cross = boundary_crossing_exposure(gdna_fl)
    gdt = GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC",
            rho=0.0,
            n_fragments=5,
            eff_length_bp=299.0,
            n_regions_used=1,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        intron=GlobalGdnaDensity(
            type="INTRON",
            rho=0.0,
            n_fragments=10,
            eff_length_bp=599.0,
            n_regions_used=1,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON",
            rho=0.0,
            n_fragments=7,
            eff_length_bp=b_cross,
            n_regions_used=1,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        gdna_fl=gdna_fl,
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    est = estimate_locus_gdna(
        locus=locus,
        n_obs=22,
        region_index=idx,
        region_arrays=ra,
        payload_arrays=pa,
        global_densities=gdt,
        gdna_fl=gdna_fl,
    )
    # κ = 0 ⇒ pure-local densities. With proration ratio == 1 (clip == full):
    #   intergenic: n_mass = 5
    #   intron:     n_mass = 10
    #   boundary:   only left side eligible (bf_left=True, bf_right=False),
    #               but the exon region is [100,600); start=100 lies in
    #               [0,1000] window ⇒ 1 eligible side. ρ = 7 / (1 · B_cross)
    #               n_boundary_observed = 7
    #               n_exon_only = ρ · L_eff_contained(500)
    leff_exon = float(_oracle_leff([500], gdna_fl)[0])
    expected_rho_b = 7.0 / b_cross
    assert est.n_gdna_intergenic == pytest.approx(5.0)
    assert est.n_gdna_intron == pytest.approx(10.0)
    assert est.n_eligible_boundaries == 1
    assert est.n_gdna_boundary_observed == pytest.approx(7.0)
    assert est.n_gdna_exon_only == pytest.approx(expected_rho_b * leff_exon)
    assert est.n_gdna == pytest.approx(5.0 + 10.0 + 7.0 + expected_rho_b * leff_exon)


def test_pi_clipped_when_overestimate():
    # Force massive intergenic over-prediction → pi > 1.
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        counts_intergenic=[1000],
    )
    gdt = GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC",
            rho=0.0,
            n_fragments=1000,
            eff_length_bp=1199.0,
            n_regions_used=1,
            kappa=KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason=""),
        ),
        intron=_gdt().intron,
        exon_intron=_gdt().exon_intron,
        gdna_fl=_delta_fl(200),
    )
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    est = estimate_locus_gdna(
        locus=locus,
        n_obs=50,  # only 50 observed but predicted = 1000
        region_index=idx,
        region_arrays=ra,
        payload_arrays=pa,
        global_densities=gdt,
        gdna_fl=_delta_fl(200),
    )
    assert est.pi_gdna == pytest.approx(1.0)
    assert est.fallback_flags & FLAG_PI_CLIPPED


def test_no_region_overlap_raises():
    ra, pa, idx, _ = _build_arrays(
        regions=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        counts_intergenic=[5],
    )
    gdt = _gdt(rho_ig=0.05, fl_mean=200)
    # Locus on a different reference id (out of range) ⇒ empty overlap ⇒ raise.
    locus = Locus(ref="chrX", ref_id=99, start=0, end=1000)
    with pytest.raises(RuntimeError, match="overlaps no regions"):
        estimate_locus_gdna(
            locus=locus,
            n_obs=10,
            region_index=idx,
            region_arrays=ra,
            payload_arrays=pa,
            global_densities=gdt,
            gdna_fl=_delta_fl(200),
        )
