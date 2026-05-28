"""Tests for ``rigel.calibration.fl`` (M7 v2)."""

from __future__ import annotations

import json

import numpy as np
import pytest

from rigel.calibration._diagnostics import Diagnostics
from rigel.calibration._fl_sources import (
    extract_gdna_counts,
    extract_global_counts,
    extract_rna_counts,
)
from rigel.calibration.fl import (
    POOL_SCORING_PRIOR_ESS,
    POOL_QUALITY_GOOD_THRESHOLD,
    POOL_QUALITY_WEAK_THRESHOLD,
    build_fl_models,
)
from rigel.calibration.scan_payload import FL_HIST_N_BINS, CalibrationScanPayload
from rigel.calibration.signature import (
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_BOUNDARY_LEFT,
    COMPARTMENT_CONTAINED,
    FL_POOL_EXONIC_CONTAINED,
    FL_POOL_INTERGENIC_BOUNDARY,
    FL_POOL_INTERGENIC_CONTAINED,
    FL_POOL_INTRONIC_BOUNDARY,
    FL_POOL_INTRONIC_CONTAINED,
    N_CHANNELS,
    N_FL_POOLS,
    N_SIGNATURES,
    SPLICE_SPLICED,
    SPLICE_UNSPLICED,
    channel_index,
    pack_signature,
)
from rigel.frag_length_model import FragmentLengthModel, FragmentLengthModels
from rigel.splice import SpliceType


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

MAX_SIZE = 1000
N_BINS = MAX_SIZE + 1
assert FL_HIST_N_BINS >= N_BINS


def _peaked(center: int, total: int) -> np.ndarray:
    """Sharp Gaussian-ish FL histogram centered at ``center`` summing to ``total``."""
    x = np.arange(N_BINS)
    raw = np.exp(-((x - center) ** 2) / (2 * 30.0**2))
    raw /= raw.sum()
    counts = np.round(raw * total).astype(np.int64)
    diff = total - int(counts.sum())
    counts[center] += diff
    return counts


def _pmf_mean(pmf: np.ndarray) -> float:
    return float(np.dot(np.arange(pmf.size, dtype=np.float64), pmf))


def _payload(
    fl_pool_rows: dict[int, np.ndarray] | None = None,
    channel_mass: np.ndarray | None = None,
    signature_mass: np.ndarray | None = None,
    n_observed: int | None = None,
) -> CalibrationScanPayload:
    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    if fl_pool_rows is not None:
        for pool_idx, vals in fl_pool_rows.items():
            values = np.asarray(vals, dtype=np.float64)
            fl_pool_mass[pool_idx, : values.size] = values

    channel = (
        np.asarray(channel_mass, dtype=np.float64)
        if channel_mass is not None
        else np.zeros(N_CHANNELS, dtype=np.float64)
    )
    signature = (
        np.asarray(signature_mass, dtype=np.float64)
        if signature_mass is not None
        else np.zeros(N_SIGNATURES, dtype=np.float64)
    )
    if n_observed is None:
        n_observed = int(max(channel.sum(), signature.sum(), fl_pool_mass.sum()))

    return CalibrationScanPayload(
        region_counts=np.zeros((0, N_CHANNELS), dtype=np.float32),
        channel_mass=channel,
        signature_mass=signature,
        fl_pool_mass=fl_pool_mass,
        fl_pool_total=fl_pool_mass.sum(axis=1),
        region_unspliced_support=np.zeros(0, dtype=np.uint64),
        region_spliced_support=np.zeros(0, dtype=np.uint64),
        n_observed=n_observed,
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_excluded_strand_ambig=0,
        n_unobserved=0,
        n_unannotated_ref=0,
        n_fl_unavailable=0,
        resolver_splicing_anchor_tolerance=0,
        n_regions=0,
    )


def _scan_trained(
    global_hist: np.ndarray,
    spliced_hist: np.ndarray | None = None,
) -> FragmentLengthModels:
    flm = FragmentLengthModels(max_size=MAX_SIZE)
    flm.global_model.counts[: global_hist.size] = global_hist.astype(np.float64)
    flm.global_model._total_weight = float(global_hist.sum())
    if spliced_hist is not None:
        spliced_model = flm.category_models[SpliceType.SPLICED_ANNOT]
        spliced_model.counts[: spliced_hist.size] = spliced_hist.astype(np.float64)
        spliced_model._total_weight = float(spliced_hist.sum())
    return flm


# ---------------------------------------------------------------------------
# Source extractors
# ---------------------------------------------------------------------------


def test_extract_global_counts_returns_float64_view_of_global_model():
    global_hist = _peaked(300, 12_345)
    scan = _scan_trained(global_hist)
    out = extract_global_counts(scan)
    assert out.dtype == np.float64
    assert out.sum() == 12_345


def test_extract_rna_counts_pulls_from_spliced_annot_bin():
    global_hist = _peaked(300, 50_000)
    spliced_hist = _peaked(250, 7_777)
    scan = _scan_trained(global_hist, spliced_hist=spliced_hist)
    out = extract_rna_counts(scan)
    assert out.dtype == np.float64
    assert out.sum() == 7_777


def test_extract_gdna_counts_sums_only_intergenic_and_intronic_fl_pools():
    payload = _payload(
        {
            FL_POOL_EXONIC_CONTAINED: np.array([1_000_000]),
            FL_POOL_INTRONIC_CONTAINED: np.array([1000]),
            FL_POOL_INTRONIC_BOUNDARY: np.array([2000]),
            FL_POOL_INTERGENIC_CONTAINED: np.array([3000]),
            FL_POOL_INTERGENIC_BOUNDARY: np.array([4000]),
        }
    )
    out = extract_gdna_counts(payload)
    assert out.dtype == np.float64
    assert int(out.sum()) == 1000 + 2000 + 3000 + 4000


# ---------------------------------------------------------------------------
# build_fl_models - quality classifier branches
# ---------------------------------------------------------------------------


def test_good_branch_pure_empirical():
    n_good = POOL_QUALITY_GOOD_THRESHOLD + 100
    fl = build_fl_models(
        global_counts=_peaked(300, 50_000),
        rna_counts=_peaked(280, n_good),
        gdna_counts=_peaked(180, n_good),
        max_size=MAX_SIZE,
    )
    assert fl.rna_quality == "good"
    assert fl.gdna_quality == "good"
    assert abs(fl.rna.mean - 280.0) < 5.0
    assert abs(fl.gdna.mean - 180.0) < 5.0


def test_good_branch_scoring_pmf_stays_empirical_for_short_fl():
    fl = build_fl_models(
        global_counts=_peaked(300, 50_000),
        rna_counts=_peaked(70, 6_000),
        gdna_counts=_peaked(70, 6_000),
        max_size=MAX_SIZE,
    )
    assert fl.rna_quality == "good"
    assert fl.gdna_quality == "good"
    assert abs(_pmf_mean(fl.rna.pmf) - 70.0) < 3.0
    assert abs(_pmf_mean(fl.gdna.pmf) - 70.0) < 3.0


def test_weak_branch_eb_borrows_from_global():
    n_weak = (POOL_QUALITY_GOOD_THRESHOLD + POOL_QUALITY_WEAK_THRESHOLD) // 2
    global_counts = _peaked(400, 50_000)
    fl_weak = build_fl_models(
        global_counts=global_counts,
        rna_counts=_peaked(180, n_weak),
        gdna_counts=_peaked(180, n_weak),
        max_size=MAX_SIZE,
    )
    fl_good = build_fl_models(
        global_counts=global_counts,
        rna_counts=_peaked(180, POOL_QUALITY_GOOD_THRESHOLD + 1000),
        gdna_counts=_peaked(180, POOL_QUALITY_GOOD_THRESHOLD + 1000),
        max_size=MAX_SIZE,
    )
    assert fl_weak.rna_quality == "weak"
    assert fl_weak.gdna_quality == "weak"
    assert fl_weak.rna.log_likelihood(400) > fl_good.rna.log_likelihood(400)
    assert fl_weak.gdna.log_likelihood(400) > fl_good.gdna.log_likelihood(400)


def test_low_nonzero_pool_builds_weak_model_without_erasing_channel_signal():
    fl = build_fl_models(
        global_counts=_peaked(300, 50_000),
        rna_counts=_peaked(70, 50),
        gdna_counts=_peaked(70, 50),
        max_size=MAX_SIZE,
    )
    assert fl.rna_quality == "weak"
    assert fl.gdna_quality == "weak"
    assert fl.rna is not fl.global_
    assert fl.gdna is not fl.global_
    assert abs(fl.rna.mean - 70.0) < 35.0
    assert abs(fl.gdna.mean - 70.0) < 35.0


def test_fallback_branch_identity_shares_global():
    n_fallback = POOL_QUALITY_WEAK_THRESHOLD - 1
    fl = build_fl_models(
        global_counts=_peaked(300, 50_000),
        rna_counts=_peaked(180, n_fallback),
        gdna_counts=_peaked(180, n_fallback),
        max_size=MAX_SIZE,
    )
    assert fl.rna_quality == "fallback"
    assert fl.gdna_quality == "fallback"
    assert fl.rna is fl.global_
    assert fl.gdna is fl.global_


def test_empty_pools_fall_back():
    fl = build_fl_models(
        global_counts=_peaked(300, 10_000),
        rna_counts=np.zeros(N_BINS, dtype=np.int64),
        gdna_counts=np.zeros(N_BINS, dtype=np.int64),
        max_size=MAX_SIZE,
    )
    assert fl.rna_quality == "fallback"
    assert fl.gdna_quality == "fallback"
    assert fl.n_rna == 0
    assert fl.n_gdna == 0


# ---------------------------------------------------------------------------
# Single-policy invariants
# ---------------------------------------------------------------------------


def test_eb_symmetry_rna_and_gdna():
    global_counts = _peaked(400, 50_000)
    pool = _peaked(180, 1500)
    fl = build_fl_models(
        global_counts=global_counts,
        rna_counts=pool,
        gdna_counts=pool,
        max_size=MAX_SIZE,
    )
    assert fl.rna_quality == fl.gdna_quality == "weak"
    grid = np.array([100, 180, 300, 400, 500, 800])
    rna_ll = np.array([fl.rna.log_likelihood(int(x)) for x in grid])
    gdna_ll = np.array([fl.gdna.log_likelihood(int(x)) for x in grid])
    np.testing.assert_array_equal(rna_ll, gdna_ll)


def test_scoring_surfaces_are_neutral_when_gdna_fallback():
    global_counts = _peaked(200, 50_000)
    rna_counts = _peaked(280, 39)
    gdna_counts = np.zeros(N_BINS, dtype=np.int64)

    fl = build_fl_models(
        global_counts=global_counts,
        rna_counts=rna_counts,
        gdna_counts=gdna_counts,
        max_size=MAX_SIZE,
    )

    assert fl.rna_quality == "weak"
    assert fl.gdna_quality == "fallback"
    assert fl.rna_fl_reliability == pytest.approx(39.0 / (39.0 + POOL_SCORING_PRIOR_ESS))
    assert fl.gdna_fl_reliability == 0.0
    assert fl.fl_contrast_weight == 0.0
    np.testing.assert_allclose(fl.rna_scoring.pmf, fl.global_.pmf, rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(fl.gdna_scoring.pmf, fl.global_.pmf, rtol=0.0, atol=1e-15)


def test_scoring_surfaces_use_joint_contrast_reliability():
    global_counts = _peaked(250, 50_000)
    rna_counts = _peaked(100, 50)
    gdna_counts = _peaked(400, 10_000)

    fl = build_fl_models(
        global_counts=global_counts,
        rna_counts=rna_counts,
        gdna_counts=gdna_counts,
        max_size=MAX_SIZE,
    )

    rna_rel = float(rna_counts.sum()) / (float(rna_counts.sum()) + POOL_SCORING_PRIOR_ESS)
    gdna_rel = float(gdna_counts.sum()) / (float(gdna_counts.sum()) + POOL_SCORING_PRIOR_ESS)
    expected_weight = min(rna_rel, gdna_rel)
    assert fl.fl_contrast_weight == pytest.approx(expected_weight)

    global_pmf = fl.global_.pmf
    rna_emp = rna_counts.astype(np.float64) / float(rna_counts.sum())
    gdna_emp = gdna_counts.astype(np.float64) / float(gdna_counts.sum())
    expected_rna = expected_weight * rna_emp + (1.0 - expected_weight) * global_pmf
    expected_gdna = expected_weight * gdna_emp + (1.0 - expected_weight) * global_pmf
    expected_rna /= expected_rna.sum()
    expected_gdna /= expected_gdna.sum()

    np.testing.assert_allclose(fl.rna_scoring.pmf, expected_rna, rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(fl.gdna_scoring.pmf, expected_gdna, rtol=1e-12, atol=1e-15)
    assert fl.rna_scoring.mean > 100.0
    assert fl.rna_scoring.mean < fl.global_.mean
    assert fl.gdna_scoring.mean > fl.global_.mean
    assert fl.gdna_scoring.mean < 400.0


def test_scoring_prior_ess_parameter_controls_contrast_weight():
    global_counts = _peaked(250, 50_000)
    rna_counts = _peaked(100, 50)
    gdna_counts = _peaked(400, 10_000)

    weakly_damped = build_fl_models(
        global_counts=global_counts,
        rna_counts=rna_counts,
        gdna_counts=gdna_counts,
        max_size=MAX_SIZE,
        scoring_prior_ess=50.0,
    )
    strongly_damped = build_fl_models(
        global_counts=global_counts,
        rna_counts=rna_counts,
        gdna_counts=gdna_counts,
        max_size=MAX_SIZE,
        scoring_prior_ess=500.0,
    )

    assert weakly_damped.fl_contrast_weight == pytest.approx(0.5)
    assert strongly_damped.fl_contrast_weight == pytest.approx(50.0 / 550.0)
    assert weakly_damped.fl_contrast_weight > strongly_damped.fl_contrast_weight
    assert abs(weakly_damped.rna_scoring.mean - weakly_damped.global_.mean) > abs(
        strongly_damped.rna_scoring.mean - strongly_damped.global_.mean
    )


def test_global_is_unconditional_anchor_no_prior():
    fl = build_fl_models(
        global_counts=_peaked(300, 100_000),
        rna_counts=np.zeros(N_BINS, dtype=np.int64),
        gdna_counts=np.zeros(N_BINS, dtype=np.int64),
        max_size=MAX_SIZE,
    )
    assert abs(fl.global_.mean - 300.0) < 1.0
    assert fl.rna is fl.global_
    assert fl.gdna is fl.global_


# ---------------------------------------------------------------------------
# EXONIC isolation: counted, not modeled as gDNA
# ---------------------------------------------------------------------------


def test_exonic_fl_pool_does_not_affect_gdna_model_via_extractor():
    intron_only = _peaked(180, 10_000)
    payload_clean = _payload({FL_POOL_INTRONIC_CONTAINED: intron_only})
    payload_with_exon = _payload(
        {
            FL_POOL_INTRONIC_CONTAINED: intron_only,
            FL_POOL_EXONIC_CONTAINED: _peaked(800, 1_000_000),
        }
    )
    scan = _scan_trained(_peaked(300, 50_000))

    fl_clean = build_fl_models(
        global_counts=extract_global_counts(scan),
        rna_counts=extract_rna_counts(scan),
        gdna_counts=extract_gdna_counts(payload_clean),
        max_size=MAX_SIZE,
    )
    fl_with = build_fl_models(
        global_counts=extract_global_counts(scan),
        rna_counts=extract_rna_counts(scan),
        gdna_counts=extract_gdna_counts(payload_with_exon),
        max_size=MAX_SIZE,
    )
    grid = np.array([100, 180, 300, 500, 800])
    for x in grid:
        assert fl_clean.gdna.log_likelihood(int(x)) == fl_with.gdna.log_likelihood(int(x))


# ---------------------------------------------------------------------------
# Diagnostics - accountability + no mask-int leakage
# ---------------------------------------------------------------------------


def test_diagnostics_total_equals_n_observed_when_no_exclusions():
    payload = _payload(n_observed=31)
    diag = Diagnostics.from_payload(payload)
    assert diag.total() == 31
    assert diag.total() == payload.n_observed


def test_diagnostics_named_fields_are_human_readable():
    channel = np.zeros(N_CHANNELS, dtype=np.float64)
    channel[channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)] = 2
    channel[channel_index(COMPARTMENT_BOUNDARY_LEFT, SPLICE_SPLICED, CHANNEL_STRAND_NEG)] = 3

    signature = np.zeros(N_SIGNATURES, dtype=np.float64)
    signature[0] = 5
    signature[pack_signature(intron_pos=True)] = 7
    signature[pack_signature(exon_neg=True)] = 11

    payload = _payload(
        fl_pool_rows={
            FL_POOL_INTERGENIC_CONTAINED: np.array([4]),
            FL_POOL_INTRONIC_BOUNDARY: np.array([6]),
        },
        channel_mass=channel,
        signature_mass=signature,
    )
    diag = Diagnostics.from_payload(payload)
    summary = diag.to_summary_dict()

    assert diag.mass_by_compartment["CONTAINED"] == 2
    assert diag.mass_by_compartment["BOUNDARY_LEFT"] == 3
    assert diag.mass_by_splice["UNSPLICED"] == 2
    assert diag.mass_by_splice["SPLICED"] == 3
    assert diag.mass_by_strand["POS"] == 2
    assert diag.mass_by_strand["NEG"] == 3
    assert diag.mass_by_coarse_class["INTERGENIC"] == 5
    assert diag.mass_by_coarse_class["INTRON"] == 7
    assert diag.mass_by_coarse_class["EXON"] == 11
    assert summary["fl_pool_total"]["INTERGENIC_CONTAINED"] == 4
    assert summary["fl_pool_total"]["INTRONIC_BOUNDARY"] == 6


def test_summary_dicts_contain_no_mask_int_keys():
    diag = Diagnostics.from_payload(
        _payload(
            channel_mass=np.ones(N_CHANNELS, dtype=np.float64),
            signature_mass=np.ones(N_SIGNATURES, dtype=np.float64),
        )
    )
    fl = build_fl_models(
        global_counts=_peaked(300, 10_000),
        rna_counts=_peaked(280, 5000),
        gdna_counts=_peaked(180, 5000),
        max_size=MAX_SIZE,
    )
    blob = json.dumps(diag.to_summary_dict()) + json.dumps(fl.to_summary_dict())
    for forbidden in (
        "mask_0",
        "mask_1",
        "mask_2",
        "mask_3",
        "mask_4",
        "mask_5",
        "mask_6",
        "mask_7",
    ):
        assert forbidden not in blob


# ---------------------------------------------------------------------------
# Frozen contract + finalized models
# ---------------------------------------------------------------------------


def test_fl_models_is_frozen():
    fl = build_fl_models(
        global_counts=_peaked(300, 10_000),
        rna_counts=_peaked(280, 5000),
        gdna_counts=_peaked(180, 5000),
        max_size=MAX_SIZE,
    )
    with pytest.raises(Exception):
        fl.rna_quality = "weak"  # type: ignore[misc]


def test_returned_models_are_finalized_and_callable():
    fl = build_fl_models(
        global_counts=_peaked(300, 10_000),
        rna_counts=_peaked(280, 5000),
        gdna_counts=_peaked(180, 5000),
        max_size=MAX_SIZE,
    )
    for model in (fl.global_, fl.rna, fl.gdna):
        assert isinstance(model, FragmentLengthModel)
        _ = model.log_likelihood(180)
        _ = model.log_likelihood(800)
