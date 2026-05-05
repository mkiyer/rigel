"""Tests for ``rigel.calibration._result_v6.CalibrationResult`` (M7)."""

from __future__ import annotations

import dataclasses
import json

import numpy as np
import pytest

from rigel.calibration._fl_pool import compute_pool_fl_models
from rigel.calibration._kappa import KappaEstimate
from rigel.calibration._result_v6 import (
    CalibrationResult,
    build_calibration_result,
    build_multi_locus_prior_df,
    build_per_locus_gdna_df,
)
from rigel.calibration.density_global import GlobalDensityTable, GlobalGdnaDensity
from rigel.calibration.locus_prior import (
    LocusGdnaEstimate,
    MultiLocusPrior,
    PriorTable,
)
from rigel.calibration.scan_payload import (
    FL_HIST_N_BINS,
    MASK_N_STATES,
    CalibrationScanPayload,
)
from rigel.frag_length_model import FragmentLengthModels
from rigel.locus import Locus


# ---------------------------------------------------------------------------
# Locked-column schemas (must match _result_v6 module constants)
# ---------------------------------------------------------------------------

EXPECTED_MULTI_LOCUS_COLUMNS = [
    "multi_locus_id", "n_obs", "n_gdna", "n_rna", "pi_gdna", "n_loci",
]
EXPECTED_PER_LOCUS_COLUMNS = [
    "multi_locus_id", "ref", "start", "end", "span",
    "n_obs", "n_gdna",
    "n_gdna_intergenic", "n_gdna_intron", "n_gdna_exon_intron",
    "pi_gdna", "n_eligible_boundaries", "fallback_flags",
]


# ---------------------------------------------------------------------------
# Fakes
# ---------------------------------------------------------------------------

def _kappa_zero() -> KappaEstimate:
    return KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason="")


def _gdt(fl_mean: float = 200.0) -> GlobalDensityTable:
    return GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC", rho=0.0, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=0, kappa=_kappa_zero(),
        ),
        intron=GlobalGdnaDensity(
            type="INTRON", rho=0.0, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=0, kappa=_kappa_zero(),
        ),
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON", rho=0.0, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=0, kappa=_kappa_zero(),
        ),
        gdna_fl_mean=fl_mean,
    )


def _payload() -> CalibrationScanPayload:
    return CalibrationScanPayload(
        global_counts=np.zeros(MASK_N_STATES, dtype=np.int64),
        per_region_counts=np.zeros((0, MASK_N_STATES), dtype=np.int64),
        fl_hist=np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64),
        u_left=np.zeros(0, dtype=np.int64),
        u_right=np.zeros(0, dtype=np.int64),
        n_observed=100, n_excluded_multimap=10, n_excluded_chimera=2,
        n_excluded_artifact=1, n_unobserved=5, n_unannotated_ref=3,
    )


def _make_estimate(
    ref: str, start: int, end: int, *, n_obs: int, n_gdna: float,
    intergenic: float = 0.0, intron: float = 0.0, exon_intron: float = 0.0,
    flags: int = 0,
) -> LocusGdnaEstimate:
    return LocusGdnaEstimate(
        locus=Locus(ref=ref, ref_id=0, start=start, end=end),
        n_obs=n_obs,
        n_gdna_intergenic=intergenic,
        n_gdna_intron=intron,
        n_gdna_exon_intron=exon_intron,
        n_gdna=n_gdna,
        pi_gdna=(n_gdna / n_obs) if n_obs > 0 else 0.0,
        rho_loco=(0.0, 0.0, 0.0),
        leff_loco=(0.0, 0.0, 0.0),
        n_eligible_boundaries=0,
        fallback_flags=flags,
    )


def _mlp(mid: int, ests: tuple[LocusGdnaEstimate, ...]) -> MultiLocusPrior:
    n_obs = sum(e.n_obs for e in ests)
    n_gdna = sum(e.n_gdna for e in ests)
    return MultiLocusPrior(
        multi_locus_id=mid, n_obs=n_obs, n_gdna=n_gdna,
        n_rna=max(0.0, n_obs - n_gdna),
        pi_gdna=(n_gdna / n_obs) if n_obs > 0 else 0.0,
        per_locus=ests,
    )


def _prior_table(mlps: tuple[MultiLocusPrior, ...]) -> PriorTable:
    n = len(mlps)
    alpha_g = np.array([m.pi_gdna * 10.0 for m in mlps], dtype=np.float64)
    alpha_r = np.array([(1.0 - m.pi_gdna) * 10.0 for m in mlps], dtype=np.float64)
    pwr = [np.ones(2, dtype=np.float32) for _ in range(n)]
    return PriorTable(
        multi_locus_priors=mlps,
        alpha_gdna=alpha_g, alpha_rna=alpha_r,
        prior_weight_rna=pwr, c_base_value=10.0,
    )


def _pool():
    """Build a real PoolFLModels via compute_pool_fl_models (no fakes)."""
    flm = FragmentLengthModels(max_size=1000)
    flm.global_model.counts[200:300] = 100.0
    flm.global_model._total_weight = float(flm.global_model.counts.sum())
    flm.global_model.finalize()
    flm.rna_model.counts[250:280] = 50.0
    flm.rna_model._total_weight = float(flm.rna_model.counts.sum())
    flm.rna_model.finalize()
    return compute_pool_fl_models(_payload(), flm)


# ---------------------------------------------------------------------------
# 1. Builder + roundtrip
# ---------------------------------------------------------------------------

def test_build_calibration_result_basic():
    est = _make_estimate("chr1", 0, 1000, n_obs=10, n_gdna=3.0, intergenic=3.0)
    mlp = _mlp(0, (est,))
    pt = _prior_table((mlp,))
    res = build_calibration_result(
        payload=_payload(), global_densities=_gdt(), pool=_pool(), prior_table=pt,
    )
    assert isinstance(res, CalibrationResult)
    assert res.n_multi_loci == 1
    assert res.payload_summary["n_observed"] == 100
    assert res.payload_summary["n_excluded_multimap"] == 10
    assert res.payload_summary["n_unannotated_ref"] == 3


def test_calibration_result_is_frozen():
    pt = _prior_table((_mlp(0, (_make_estimate("chr1", 0, 1, n_obs=1, n_gdna=0.0),)),))
    res = build_calibration_result(
        payload=_payload(), global_densities=_gdt(), pool=_pool(), prior_table=pt,
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        res.n_multi_loci = 99  # type: ignore[misc]


# ---------------------------------------------------------------------------
# 2. Convenience accessors (zero-copy / identity)
# ---------------------------------------------------------------------------

def test_convenience_accessors_alias_prior_table():
    est = _make_estimate("chr1", 0, 1000, n_obs=10, n_gdna=4.0, intergenic=4.0)
    pt = _prior_table((_mlp(0, (est,)),))
    res = build_calibration_result(
        payload=_payload(), global_densities=_gdt(fl_mean=185.0), pool=_pool(), prior_table=pt,
    )
    assert res.alpha_gdna is pt.alpha_gdna
    assert res.alpha_rna is pt.alpha_rna
    assert res.prior_weight_rna is pt.prior_weight_rna
    assert res.gdna_fl_mean == res.pool.gdna_fl_mean


# ---------------------------------------------------------------------------
# 3. with_priors → new instance, dataframes rebuilt
# ---------------------------------------------------------------------------

def test_with_priors_returns_new_instance_and_rebuilds_dfs():
    est0 = _make_estimate("chr1", 0, 1000, n_obs=10, n_gdna=2.0, intergenic=2.0)
    pt0 = _prior_table((_mlp(0, (est0,)),))
    res0 = build_calibration_result(
        payload=_payload(), global_densities=_gdt(), pool=_pool(), prior_table=pt0,
    )

    est1a = _make_estimate("chr2", 0, 500, n_obs=20, n_gdna=5.0, intergenic=5.0)
    est1b = _make_estimate("chr2", 1000, 1500, n_obs=5, n_gdna=1.0, intron=1.0)
    pt1 = _prior_table((_mlp(0, (est1a, est1b)),))
    res1 = res0.with_priors(pt1)

    assert res1 is not res0
    assert res0.prior_table is pt0      # original untouched
    assert res1.prior_table is pt1
    assert res1.n_multi_loci == 1
    assert len(res1.per_locus_gdna_df) == 2
    assert res1.multi_locus_prior_df.iloc[0]["n_loci"] == 2


# ---------------------------------------------------------------------------
# 4. Diagnostic dataframes — locked schemas
# ---------------------------------------------------------------------------

def test_multi_locus_prior_df_schema_and_content():
    est = _make_estimate("chr1", 0, 1000, n_obs=10, n_gdna=3.0, intergenic=3.0)
    df = build_multi_locus_prior_df((_mlp(0, (est,)),))
    assert list(df.columns) == EXPECTED_MULTI_LOCUS_COLUMNS
    row = df.iloc[0]
    assert row["multi_locus_id"] == 0
    assert row["n_obs"] == 10
    assert row["n_gdna"] == pytest.approx(3.0)
    assert row["n_rna"] == pytest.approx(7.0)
    assert row["pi_gdna"] == pytest.approx(0.3)
    assert row["n_loci"] == 1


def test_per_locus_gdna_df_schema_and_content():
    e0 = _make_estimate("chr1", 0, 100, n_obs=10, n_gdna=2.0, intergenic=2.0, flags=1)
    e1 = _make_estimate("chr2", 500, 700, n_obs=5, n_gdna=1.0, intron=1.0)
    df = build_per_locus_gdna_df((_mlp(0, (e0, e1)),))
    assert list(df.columns) == EXPECTED_PER_LOCUS_COLUMNS
    assert len(df) == 2
    assert df.iloc[0]["ref"] == "chr1"
    assert df.iloc[0]["span"] == 100
    assert df.iloc[0]["n_gdna_intergenic"] == pytest.approx(2.0)
    assert df.iloc[0]["fallback_flags"] == 1
    assert df.iloc[1]["ref"] == "chr2"
    assert df.iloc[1]["n_gdna_intron"] == pytest.approx(1.0)


def test_dataframes_handle_empty_input():
    df_m = build_multi_locus_prior_df(())
    df_p = build_per_locus_gdna_df(())
    assert list(df_m.columns) == EXPECTED_MULTI_LOCUS_COLUMNS
    assert list(df_p.columns) == EXPECTED_PER_LOCUS_COLUMNS
    assert len(df_m) == 0
    assert len(df_p) == 0


# ---------------------------------------------------------------------------
# 5. Summary serialisation
# ---------------------------------------------------------------------------

def test_to_summary_dict_is_json_serialisable():
    est = _make_estimate("chr1", 0, 1000, n_obs=10, n_gdna=3.0, intergenic=3.0)
    pt = _prior_table((_mlp(0, (est,)),))
    res = build_calibration_result(
        payload=_payload(), global_densities=_gdt(), pool=_pool(), prior_table=pt,
    )
    blob = json.dumps(res.to_summary_dict())
    assert "global_densities" in blob
    assert "pool" in blob
    assert "payload_summary" in blob
    assert "mean_pi_gdna" in blob
    assert "c_base" in blob


def test_to_summary_dict_mean_pi_gdna_handles_empty_priors():
    pt = PriorTable(
        multi_locus_priors=(),
        alpha_gdna=np.zeros(0, dtype=np.float64),
        alpha_rna=np.zeros(0, dtype=np.float64),
        prior_weight_rna=[],
        c_base_value=10.0,
    )
    res = build_calibration_result(
        payload=_payload(), global_densities=_gdt(), pool=_pool(), prior_table=pt,
    )
    s = res.to_summary_dict()
    assert s["mean_pi_gdna"] == 0.0
    assert s["n_multi_loci"] == 0
