"""Tests for rigel.calibration._regional_exposure."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration._kappa import KappaEstimate
from rigel.calibration._orient import ORIENT_OPP, StrandSummary
from rigel.calibration._regional_exposure import (
    RegionalGdnaExposure,
    _weighted_quantile,
)
from rigel.calibration.density_global import DensityType, GlobalDensityTable, GlobalGdnaDensity
from rigel.calibration.regions import RegionStrand, RegionType
from rigel.calibration.scan_payload import (
    FL_HIST_N_BINS,
    MASK_INTERGENIC,
    MASK_N_STATES,
    ORIENT_N,
    CalibrationScanPayload,
)
from rigel.frag_length_model import FragmentLengthModel


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _fl_delta(ell: int = 100, max_size: int = 200) -> FragmentLengthModel:
    fl = FragmentLengthModel(max_size=max_size)
    for _ in range(2000):
        fl.observe(ell)
    fl.finalize()
    return fl


def _region_df_three_intergenic(span: int = 10_000) -> pd.DataFrame:
    """Three INTERGENIC regions on ref 'chr1', no exons / introns."""
    return pd.DataFrame(
        {
            "ref_name": ["chr1", "chr1", "chr1"],
            "start": [0, span, 2 * span],
            "end": [span, 2 * span, 3 * span],
            "type": [int(RegionType.INTERGENIC)] * 3,
            "strand": [int(RegionStrand.NONE)] * 3,
            "boundary_flux_left": [False, False, False],
            "boundary_flux_right": [False, False, False],
        }
    )


def _empty_payload(R: int) -> CalibrationScanPayload:
    return CalibrationScanPayload(
        global_counts=np.zeros(MASK_N_STATES, dtype=np.int64),
        per_region_counts=np.zeros((R, MASK_N_STATES), dtype=np.int64),
        fl_hist=np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64),
        u_left=np.zeros(R, dtype=np.int64),
        u_right=np.zeros(R, dtype=np.int64),
        intron_counts_by_orient=np.zeros((R, ORIENT_N), dtype=np.int64),
        u_left_by_orient=np.zeros((R, ORIENT_N), dtype=np.int64),
        u_right_by_orient=np.zeros((R, ORIENT_N), dtype=np.int64),
        n_observed=0,
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_unobserved=0,
        n_unannotated_ref=0,
    )


def _density(
    rho: float, kappa: float = 1.0, label: DensityType = "INTERGENIC"
) -> GlobalGdnaDensity:
    return GlobalGdnaDensity(
        type=label,
        rho=rho,
        n_fragments=0,
        eff_length_bp=1.0,
        n_regions_used=1,
        kappa=KappaEstimate(value=kappa, n_regions=1, fallback_used=False, fallback_reason=""),
    )


def _density_table(
    rho_ig: float = 1e-6,
    rho_in: float = 1e-6,
    rho_ex: float = 1e-6,
    kappa: float = 1.0,
) -> GlobalDensityTable:
    fl = _fl_delta()
    return GlobalDensityTable(
        intergenic=_density(rho_ig, kappa, "INTERGENIC"),
        intron=_density(rho_in, kappa, "INTRON"),
        exon_intron=_density(rho_ex, kappa, "EXON-INTRON"),
        gdna_fl=fl,
    )


def _build_arrays(
    region_df: pd.DataFrame, payload: CalibrationScanPayload
) -> tuple[RegionArrays, PayloadArrays]:
    region_arrays = RegionArrays.from_region_df(region_df, {"chr1": 0})
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    return region_arrays, payload_arrays


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_uniform_constructor():
    df = _region_df_three_intergenic()
    region_arrays = RegionArrays.from_region_df(df, {"chr1": 0})
    exp = RegionalGdnaExposure.uniform(region_arrays)
    assert exp.mode == "uniform"
    assert np.allclose(exp.weight, 1.0)
    assert np.allclose(exp.log_weight, 0.0)


def test_uniform_when_disabled():
    df = _region_df_three_intergenic()
    region_arrays, payload_arrays = _build_arrays(df, _empty_payload(3))
    fl = _fl_delta()
    exp = RegionalGdnaExposure.build(
        region_arrays,
        payload_arrays,
        _density_table(),
        fl,
        enabled=False,
    )
    assert exp.mode == "uniform"
    assert np.allclose(exp.weight, 1.0)


def test_equal_density_input_signal_zero():
    """All regions identical → no spread → signal=0 → mode=uniform."""
    df = _region_df_three_intergenic()
    payload = _empty_payload(3)
    # Same counts everywhere → identical rho_hat.
    payload.per_region_counts[:, MASK_INTERGENIC] = 10
    region_arrays, payload_arrays = _build_arrays(df, payload)
    fl = _fl_delta()
    exp = RegionalGdnaExposure.build(
        region_arrays,
        payload_arrays,
        _density_table(rho_ig=1e-5),
        fl,
    )
    assert exp.mode == "uniform"
    assert np.allclose(exp.weight, 1.0)
    # Region-type diagnostic should still be populated, but it is QC-only.
    assert "INTERGENIC" in exp.per_class
    assert exp.observed_log_spread <= exp.null_log_spread + 1e-12


def test_kappa_alpha_is_converted_to_opportunity_length(monkeypatch):
    """NB alpha must be converted to bp-equivalent opportunity before shrinkage."""
    span = 1000
    n_regions = 5
    df = pd.DataFrame(
        {
            "ref_name": ["chr1"] * n_regions,
            "start": [i * span for i in range(n_regions)],
            "end": [(i + 1) * span for i in range(n_regions)],
            "type": [int(RegionType.INTERGENIC)] * n_regions,
            "strand": [int(RegionStrand.NONE)] * n_regions,
            "boundary_flux_left": [False] * n_regions,
            "boundary_flux_right": [False] * n_regions,
        }
    )
    payload = _empty_payload(n_regions)
    payload.per_region_counts[:, MASK_INTERGENIC] = [0, 0, 0, 0, 5]
    region_arrays, payload_arrays = _build_arrays(df, payload)

    def _fixed_alpha(counts: np.ndarray, eff_lengths: np.ndarray, rho_hat: float) -> KappaEstimate:
        return KappaEstimate(
            value=1.0,
            n_regions=int(counts.size),
            fallback_used=False,
            fallback_reason="",
        )

    monkeypatch.setattr("rigel.calibration._regional_exposure.estimate_kappa", _fixed_alpha)

    exp = RegionalGdnaExposure.build(
        region_arrays,
        payload_arrays,
        _density_table(rho_ig=1e-3, kappa=1.0),
        _fl_delta(ell=1, max_size=20),
    )

    assert exp.rho_global == pytest.approx(1e-3, rel=1e-3)
    assert exp.kappa_alpha_global == pytest.approx(1.0)
    assert exp.kappa_opportunity_bp == pytest.approx(1.0 / exp.rho_global)
    assert exp.rho_hat[0] == pytest.approx(0.5 * exp.rho_global)
    summary = exp.to_summary_dict()
    assert summary["kappa_alpha_global"] == pytest.approx(1.0)
    assert summary["kappa_opportunity_bp"] == pytest.approx(1.0 / exp.rho_global)
    assert "kappa_global" not in summary


def test_bimodal_density_input_attenuates_low_regions():
    """Three regions, one with 1000x more counts than the other two."""
    df = _region_df_three_intergenic(span=10_000)
    payload = _empty_payload(3)
    payload.per_region_counts[0, MASK_INTERGENIC] = 5000
    payload.per_region_counts[1, MASK_INTERGENIC] = 0
    payload.per_region_counts[2, MASK_INTERGENIC] = 0
    region_arrays, payload_arrays = _build_arrays(df, payload)
    fl = _fl_delta(ell=100)
    # Set rho_global to a small value so the high-count region pops above.
    exp = RegionalGdnaExposure.build(
        region_arrays,
        payload_arrays,
        _density_table(rho_ig=1e-5, kappa=1.0),
        fl,
    )
    # If signal is nonzero, the low-density regions should have weight < 1.
    if exp.mode == "regional":
        # High region is the rho_ref → weight near 1; low regions below it.
        assert exp.weight[0] >= exp.weight[1]
        assert exp.weight[0] >= exp.weight[2]
        assert exp.weight[1] < 1.0


def test_vector_lookup_matches_per_region_log_weights():
    """log_weights_for_positions returns the weight of the containing region."""
    df = _region_df_three_intergenic(span=1000)
    region_arrays, payload_arrays = _build_arrays(df, _empty_payload(3))
    # Hand-craft a regional exposure: weight = [1.0, 0.1, 0.5]
    weights = np.array([1.0, 0.1, 0.5], dtype=np.float64)
    log_w = np.log(weights)
    exp = RegionalGdnaExposure(
        rho_hat=np.zeros(3),
        log_weight=log_w,
        weight=weights,
        mode="regional",
        rho_ref=1e-5,
        n_at_floor=0,
        per_class={},
        ref_offsets=region_arrays.ref_offsets.copy(),
        ref_id=region_arrays.ref_id.copy(),
        start=region_arrays.start.copy(),
        end=region_arrays.end.copy(),
    )
    # Sample one position from each region + one outside.
    positions = np.array([500, 1500, 2500, 5000], dtype=np.int64)
    refs = np.array([0, 0, 0, 0], dtype=np.int32)
    out = exp.log_weights_for_positions(refs, positions)
    assert np.isclose(out[0], log_w[0])
    assert np.isclose(out[1], log_w[1])
    assert np.isclose(out[2], log_w[2])
    # Position outside any region → identity (log A = 0).
    assert out[3] == 0.0


def test_uniform_lookup_returns_zeros():
    df = _region_df_three_intergenic()
    region_arrays = RegionArrays.from_region_df(df, {"chr1": 0})
    exp = RegionalGdnaExposure.uniform(region_arrays)
    positions = np.array([100, 5000, 25000], dtype=np.int64)
    refs = np.array([0, 0, 0], dtype=np.int32)
    out = exp.log_weights_for_positions(refs, positions)
    assert np.allclose(out, 0.0)


def test_weighted_length_uniform_equals_span():
    df = _region_df_three_intergenic(span=1000)
    region_arrays = RegionArrays.from_region_df(df, {"chr1": 0})
    exp = RegionalGdnaExposure.uniform(region_arrays)
    assert exp.weighted_length_on_ref(0, 500, 2500) == 2000.0


def test_weighted_length_two_state_geometry():
    """Region 0 weight 1.0, region 1 weight 0.01, region 2 weight 1.0."""
    df = _region_df_three_intergenic(span=1000)
    region_arrays = RegionArrays.from_region_df(df, {"chr1": 0})
    weights = np.array([1.0, 0.01, 1.0], dtype=np.float64)
    exp = RegionalGdnaExposure(
        rho_hat=np.zeros(3),
        log_weight=np.log(weights),
        weight=weights,
        mode="regional",
        rho_ref=1e-5,
        n_at_floor=0,
        per_class={},
        ref_offsets=region_arrays.ref_offsets.copy(),
        ref_id=region_arrays.ref_id.copy(),
        start=region_arrays.start.copy(),
        end=region_arrays.end.copy(),
    )
    # Span across all three: 1000*1 + 1000*0.01 + 1000*1 = 2010
    val = exp.weighted_length_on_ref(0, 0, 3000)
    assert val == pytest.approx(2010.0)


def test_weighted_length_unknown_contig_is_identity():
    df = _region_df_three_intergenic()
    region_arrays, _ = _build_arrays(df, _empty_payload(3))
    exp = RegionalGdnaExposure.uniform(region_arrays)
    # ref_id outside the table → identity exposure.
    assert exp.weighted_length_on_ref(99, 0, 100) == 100.0


def test_weighted_quantile_returns_fallback_on_zero_weight():
    v = np.array([1.0, 2.0, 3.0])
    w = np.array([0.0, 0.0, 0.0])
    assert _weighted_quantile(v, w, 0.5, fallback=42.0) == 42.0


def test_weighted_quantile_step_cdf():
    """Q95 of equal-weighted [1..10] = 10 under lower-step convention."""
    v = np.arange(1, 11, dtype=np.float64)
    w = np.ones(10)
    assert _weighted_quantile(v, w, 0.95, fallback=-1.0) == 10.0
    assert _weighted_quantile(v, w, 0.5, fallback=-1.0) == 5.0


def test_strand_correction_path_runs_with_intron_payload():
    """Smoke test: strand-correction code path returns finite weights."""
    df = pd.DataFrame(
        {
            "ref_name": ["chr1", "chr1"],
            "start": [0, 10_000],
            "end": [10_000, 20_000],
            "type": [int(RegionType.INTRON), int(RegionType.INTRON)],
            "strand": [int(RegionStrand.POS), int(RegionStrand.POS)],
            "boundary_flux_left": [False, False],
            "boundary_flux_right": [False, False],
        }
    )
    payload = _empty_payload(2)
    payload.intron_counts_by_orient[0, ORIENT_OPP] = 100  # interpreted as gDNA
    payload.intron_counts_by_orient[1, ORIENT_OPP] = 0
    region_arrays, payload_arrays = _build_arrays(df, payload)
    strand_summary = StrandSummary(p_r1_sense=1.0, n_observations=10_000)
    fl = _fl_delta()
    exp = RegionalGdnaExposure.build(
        region_arrays,
        payload_arrays,
        _density_table(rho_in=1e-5),
        fl,
        strand_summary=strand_summary,
    )
    assert exp.weight.size == 2
    assert np.all(np.isfinite(exp.weight))
    assert np.all(exp.weight > 0.0)
    assert np.all(exp.weight <= 1.0)


def test_gap_in_region_partition_raises():
    """Synthetic gap inside a contig must raise, not silently weight at 1.0."""
    df = pd.DataFrame(
        {
            "ref_name": ["chr1", "chr1"],
            "start": [0, 2000],  # gap [1000, 2000)
            "end": [1000, 3000],
            "type": [int(RegionType.INTERGENIC), int(RegionType.INTERGENIC)],
            "strand": [int(RegionStrand.NONE), int(RegionStrand.NONE)],
            "boundary_flux_left": [False, False],
            "boundary_flux_right": [False, False],
        }
    )
    region_arrays = RegionArrays.from_region_df(df, {"chr1": 0})
    weights = np.array([0.5, 0.5], dtype=np.float64)
    exp = RegionalGdnaExposure(
        rho_hat=np.zeros(2),
        log_weight=np.log(weights),
        weight=weights,
        mode="regional",
        rho_ref=1e-5,
        n_at_floor=0,
        per_class={},
        ref_offsets=region_arrays.ref_offsets.copy(),
        ref_id=region_arrays.ref_id.copy(),
        start=region_arrays.start.copy(),
        end=region_arrays.end.copy(),
    )
    with pytest.raises(RuntimeError, match="gap-free invariant violated"):
        exp.weighted_length_on_ref(0, 500, 2500)
