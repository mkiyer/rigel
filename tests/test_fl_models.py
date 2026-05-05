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
    POOL_QUALITY_GOOD_THRESHOLD,
    POOL_QUALITY_WEAK_THRESHOLD,
    build_fl_models,
)
from rigel.calibration.scan_payload import (
    FL_HIST_N_BINS,
    MASK_EXON,
    MASK_INTERGENIC,
    MASK_INTRON,
    MASK_N_STATES,
    CalibrationScanPayload,
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
    """Sharp Gaussian-ish FL histogram centred at ``center`` summing to ``total``."""
    x = np.arange(N_BINS)
    raw = np.exp(-((x - center) ** 2) / (2 * 30.0 ** 2))
    raw /= raw.sum()
    counts = np.round(raw * total).astype(np.int64)
    diff = total - int(counts.sum())
    counts[center] += diff
    return counts


def _payload(
    fl_hist_rows: dict[int, np.ndarray] | None = None,
    global_counts: np.ndarray | None = None,
) -> CalibrationScanPayload:
    h = np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64)
    if fl_hist_rows is not None:
        for mask, vals in fl_hist_rows.items():
            v = np.asarray(vals, dtype=np.int64)
            h[mask, : v.size] = v
    gc = (
        np.asarray(global_counts, dtype=np.int64)
        if global_counts is not None
        else np.zeros(MASK_N_STATES, dtype=np.int64)
    )
    n_observed = int(gc.sum())
    return CalibrationScanPayload(
        global_counts=gc,
        per_region_counts=np.zeros((0, MASK_N_STATES), dtype=np.int64),
        fl_hist=h,
        u_left=np.zeros(0, dtype=np.int64),
        u_right=np.zeros(0, dtype=np.int64),
        n_observed=n_observed, n_excluded_multimap=0, n_excluded_chimera=0,
        n_excluded_artifact=0, n_unobserved=0, n_unannotated_ref=0,
    )


def _scan_trained(global_hist: np.ndarray, spliced_hist: np.ndarray | None = None) -> FragmentLengthModels:
    flm = FragmentLengthModels(max_size=MAX_SIZE)
    flm.global_model.counts[: global_hist.size] = global_hist.astype(np.float64)
    flm.global_model._total_weight = float(global_hist.sum())
    if spliced_hist is not None:
        sm = flm.category_models[SpliceType.SPLICED_ANNOT]
        sm.counts[: spliced_hist.size] = spliced_hist.astype(np.float64)
        sm._total_weight = float(spliced_hist.sum())
    return flm


# ---------------------------------------------------------------------------
# Source extractors
# ---------------------------------------------------------------------------

def test_extract_global_counts_returns_int64_view_of_global_model():
    g = _peaked(300, 12_345)
    scan = _scan_trained(g)
    out = extract_global_counts(scan)
    assert out.dtype == np.int64
    assert out.sum() == 12_345


def test_extract_rna_counts_pulls_from_spliced_annot_bin():
    g = _peaked(300, 50_000)
    s = _peaked(250, 7_777)
    scan = _scan_trained(g, spliced_hist=s)
    out = extract_rna_counts(scan)
    assert out.dtype == np.int64
    assert out.sum() == 7_777


def test_extract_gdna_counts_sums_only_unspliced_pool_masks():
    """Mask 0b001 (EXON_ONLY) MUST NOT contribute — spliced-genomic-span trap."""
    payload = _payload({
        MASK_EXON:                       np.array([1_000_000]),  # ignored
        MASK_INTRON:                     np.array([1000]),
        MASK_EXON | MASK_INTRON:         np.array([2000]),
        MASK_INTERGENIC:                 np.array([3000]),
        MASK_EXON | MASK_INTERGENIC:     np.array([100]),       # ignored
        MASK_INTRON | MASK_INTERGENIC:   np.array([100]),       # ignored
        MASK_EXON | MASK_INTRON | MASK_INTERGENIC: np.array([100]),  # ignored
    })
    out = extract_gdna_counts(payload)
    assert out.dtype == np.int64
    assert int(out.sum()) == 1000 + 2000 + 3000


# ---------------------------------------------------------------------------
# build_fl_models — quality classifier branches
# ---------------------------------------------------------------------------

def test_good_branch_pure_empirical():
    n_g = POOL_QUALITY_GOOD_THRESHOLD + 100
    fl = build_fl_models(
        global_counts=_peaked(300, 50_000),
        rna_counts=_peaked(280, n_g),
        gdna_counts=_peaked(180, n_g),
        max_size=MAX_SIZE,
    )
    assert fl.rna_quality == "good"
    assert fl.gdna_quality == "good"
    assert abs(fl.rna.mean - 280.0) < 5.0
    assert abs(fl.gdna.mean - 180.0) < 5.0


def test_weak_branch_eb_borrows_from_global():
    n_g = (POOL_QUALITY_GOOD_THRESHOLD + POOL_QUALITY_WEAK_THRESHOLD) // 2
    global_counts = _peaked(400, 50_000)
    fl_weak = build_fl_models(
        global_counts=global_counts,
        rna_counts=_peaked(180, n_g),
        gdna_counts=_peaked(180, n_g),
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
    # Under EB, the weak pool borrows mass from global (peak at 400),
    # giving strictly higher likelihood there than the pure-empirical good fit.
    assert fl_weak.rna.log_likelihood(400) > fl_good.rna.log_likelihood(400)
    assert fl_weak.gdna.log_likelihood(400) > fl_good.gdna.log_likelihood(400)


def test_fallback_branch_identity_shares_global():
    n = POOL_QUALITY_WEAK_THRESHOLD - 1
    fl = build_fl_models(
        global_counts=_peaked(300, 50_000),
        rna_counts=_peaked(180, n),
        gdna_counts=_peaked(180, n),
        max_size=MAX_SIZE,
    )
    assert fl.rna_quality == "fallback"
    assert fl.gdna_quality == "fallback"
    # Identity-share — no recomputation.
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
    """The single EB policy must produce IDENTICAL rna and gdna models when
    fed identical raw counts.  This is the structural guard against
    accidental policy divergence."""
    global_counts = _peaked(400, 50_000)
    pool = _peaked(180, 1500)   # weak branch for both
    fl = build_fl_models(
        global_counts=global_counts,
        rna_counts=pool, gdna_counts=pool,
        max_size=MAX_SIZE,
    )
    assert fl.rna_quality == fl.gdna_quality == "weak"
    grid = np.array([100, 180, 300, 400, 500, 800])
    rna_ll = np.array([fl.rna.log_likelihood(int(x)) for x in grid])
    gdna_ll = np.array([fl.gdna.log_likelihood(int(x)) for x in grid])
    np.testing.assert_array_equal(rna_ll, gdna_ll)


def test_global_is_unconditional_anchor_no_prior():
    """``global_`` is finalized with NO prior (it IS the prior).
    Verify by feeding only global counts: rna/gdna fall back, but the
    fallback ``is global_`` and global mean equals the empirical mean.
    """
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
# EXON_ONLY isolation: counted, not modelled
# ---------------------------------------------------------------------------

def test_exon_only_does_not_affect_gdna_model_via_extractor():
    """One million EXON_ONLY observations in the gDNA-pool histogram
    column must NOT shift any FL-bin log-likelihood of ``gdna``."""
    intron_only = _peaked(180, 10_000)
    payload_clean = _payload({MASK_INTRON: intron_only})
    payload_with_exon = _payload({
        MASK_INTRON: intron_only,
        MASK_EXON:   _peaked(800, 1_000_000),
    })
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
# Diagnostics — accountability + no mask-int leakage
# ---------------------------------------------------------------------------

def test_diagnostics_total_equals_n_observed():
    gc = np.array([7, 11, 13, 17, 19, 23, 29, 31], dtype=np.int64)
    payload = _payload(global_counts=gc)
    diag = Diagnostics.from_payload(payload)
    assert diag.total() == int(gc.sum())
    assert diag.total() == payload.n_observed


def test_diagnostics_named_fields_match_mask_codes():
    gc = np.zeros(MASK_N_STATES, dtype=np.int64)
    gc[0] = 1
    gc[MASK_EXON] = 2
    gc[MASK_INTRON] = 3
    gc[MASK_EXON | MASK_INTRON] = 4
    gc[MASK_INTERGENIC] = 5
    gc[MASK_EXON | MASK_INTERGENIC] = 6
    gc[MASK_INTRON | MASK_INTERGENIC] = 7
    gc[MASK_EXON | MASK_INTRON | MASK_INTERGENIC] = 8
    diag = Diagnostics.from_payload(_payload(global_counts=gc))
    assert diag.n_unannotated == 1
    assert diag.n_exon_only == 2
    assert diag.n_intron_only == 3
    assert diag.n_exon_intron == 4
    assert diag.n_intergenic_only == 5
    assert diag.n_exon_intergenic == 6
    assert diag.n_intron_intergenic == 7
    assert diag.n_exon_intron_intergenic == 8


def test_summary_dicts_contain_no_mask_int_keys():
    """Regression guard: no public surface may name mask integers."""
    diag = Diagnostics.from_payload(_payload(global_counts=np.ones(MASK_N_STATES, dtype=np.int64)))
    fl = build_fl_models(
        global_counts=_peaked(300, 10_000),
        rna_counts=_peaked(280, 5000),
        gdna_counts=_peaked(180, 5000),
        max_size=MAX_SIZE,
    )
    blob = json.dumps(diag.to_summary_dict()) + json.dumps(fl.to_summary_dict())
    for forbidden in ("mask_0", "mask_1", "mask_2", "mask_3", "mask_4",
                      "mask_5", "mask_6", "mask_7"):
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
    for m in (fl.global_, fl.rna, fl.gdna):
        assert isinstance(m, FragmentLengthModel)
        # log_likelihood must answer without raising on any in-range FL.
        _ = m.log_likelihood(180)
        _ = m.log_likelihood(800)
