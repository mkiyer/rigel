"""Tests for ``rigel.calibration._fl_pool`` (M7)."""

from __future__ import annotations

import json

import numpy as np
import pytest

from rigel.calibration._fl_pool import (
    POOL_EB_PRIOR_ESS,
    POOL_QUALITY_GOOD_THRESHOLD,
    POOL_QUALITY_WEAK_THRESHOLD,
    PoolFLModels,
    compute_pool_fl_models,
)
from rigel.calibration.scan_payload import (
    FL_HIST_N_BINS,
    MASK_N_STATES,
    CalibrationScanPayload,
)
from rigel.frag_length_model import FragmentLengthModel, FragmentLengthModels


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

MAX_SIZE = 1000  # FragmentLengthModel default
N_BINS = MAX_SIZE + 1
assert FL_HIST_N_BINS >= N_BINS


def _fl_hist(rows: dict[int, np.ndarray]) -> np.ndarray:
    h = np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64)
    for mask, vals in rows.items():
        v = np.asarray(vals, dtype=np.int64)
        h[mask, : v.size] = v
    return h


def _payload(fl_hist: np.ndarray) -> CalibrationScanPayload:
    return CalibrationScanPayload(
        global_counts=np.zeros(MASK_N_STATES, dtype=np.int64),
        per_region_counts=np.zeros((0, MASK_N_STATES), dtype=np.int64),
        fl_hist=fl_hist,
        u_left=np.zeros(0, dtype=np.int64),
        u_right=np.zeros(0, dtype=np.int64),
        n_observed=0, n_excluded_multimap=0, n_excluded_chimera=0,
        n_excluded_artifact=0, n_unobserved=0, n_unannotated_ref=0,
    )


def _scan_trained_with_global(
    global_counts: np.ndarray,
    rna_counts: np.ndarray | None = None,
) -> FragmentLengthModels:
    flm = FragmentLengthModels(max_size=MAX_SIZE)
    flm.global_model.counts[: global_counts.size] = global_counts.astype(np.float64)
    flm.global_model._total_weight = float(global_counts.sum())
    if rna_counts is not None:
        flm.rna_model.counts[: rna_counts.size] = rna_counts.astype(np.float64)
        flm.rna_model._total_weight = float(rna_counts.sum())
    flm.global_model.finalize()
    flm.rna_model.finalize()
    return flm


def _hist_peaked_at(center: int, total: int) -> np.ndarray:
    """Sharp Gaussian-ish histogram centred at ``center`` summing to ``total``."""
    x = np.arange(N_BINS)
    raw = np.exp(-((x - center) ** 2) / (2 * 30.0 ** 2))
    raw /= raw.sum()
    counts = np.round(raw * total).astype(np.int64)
    diff = total - int(counts.sum())
    counts[center] += diff
    return counts


# ---------------------------------------------------------------------------
# 1. Quality classifier branches
# ---------------------------------------------------------------------------

def test_good_branch_pure_empirical():
    n_g = POOL_QUALITY_GOOD_THRESHOLD + 100
    intron = _hist_peaked_at(180, n_g)
    fl_hist = _fl_hist({2: intron})
    payload = _payload(fl_hist)
    scan = _scan_trained_with_global(_hist_peaked_at(250, 10_000))

    pool = compute_pool_fl_models(payload, scan)

    assert pool.gdna_quality == "good"
    assert pool.gdna_n_fragments == n_g
    assert pool.gdna_used_global_fallback is False
    assert pool.gdna_eb_ess == 0.0
    # Mean stays close to the empirical centre (180), not the global centre (250).
    assert abs(pool.gdna_fl_mean - 180.0) < 5.0


def test_weak_branch_eb_shrinks_toward_global():
    """In the weak branch, the global histogram is folded into the
    posterior predictive log-prob (via Dirichlet-Multinomial smoothing).
    The empirical histogram is unchanged (so ``mean`` is unchanged), but
    ``log_likelihood`` at a point where the global has mass and the
    empirical has none must be measurably higher under ``weak`` than
    under ``good``.
    """
    n_g_weak = (POOL_QUALITY_GOOD_THRESHOLD + POOL_QUALITY_WEAK_THRESHOLD) // 2
    n_g_good = POOL_QUALITY_GOOD_THRESHOLD + 1000
    intron_weak = _hist_peaked_at(180, n_g_weak)
    intron_good = _hist_peaked_at(180, n_g_good)
    global_counts = _hist_peaked_at(400, 50_000)

    pool_weak = compute_pool_fl_models(
        _payload(_fl_hist({2: intron_weak})),
        _scan_trained_with_global(global_counts),
    )
    pool_good = compute_pool_fl_models(
        _payload(_fl_hist({2: intron_good})),
        _scan_trained_with_global(global_counts),
    )

    assert pool_weak.gdna_quality == "weak"
    assert pool_weak.gdna_eb_ess == POOL_EB_PRIOR_ESS
    assert pool_weak.gdna_used_global_fallback is False
    # 400 is the global peak, far from the empirical peak (180).
    # Under EB shrinkage, the weak branch borrows mass from the global,
    # giving a strictly higher log_likelihood at 400.
    assert pool_weak.gdna.log_likelihood(400) > pool_good.gdna.log_likelihood(400)


def test_fallback_branch_uses_global():
    n_g = POOL_QUALITY_WEAK_THRESHOLD - 1
    intron = _hist_peaked_at(180, n_g)
    fl_hist = _fl_hist({2: intron})
    payload = _payload(fl_hist)
    global_counts = _hist_peaked_at(400, 50_000)
    scan = _scan_trained_with_global(global_counts)

    pool = compute_pool_fl_models(payload, scan)

    assert pool.gdna_quality == "fallback"
    assert pool.gdna_used_global_fallback is True
    assert pool.gdna_n_fragments == n_g
    # Identity-ish: pool gDNA mean ≈ global mean.
    assert abs(pool.gdna_fl_mean - 400.0) < 1.0


def test_empty_pool_falls_back():
    fl_hist = _fl_hist({})
    payload = _payload(fl_hist)
    scan = _scan_trained_with_global(_hist_peaked_at(300, 10_000))
    pool = compute_pool_fl_models(payload, scan)
    assert pool.gdna_quality == "fallback"
    assert pool.gdna_n_fragments == 0
    assert pool.gdna_used_global_fallback is True


# ---------------------------------------------------------------------------
# 2. Mask aggregation correctness
# ---------------------------------------------------------------------------

def test_gdna_pool_sums_masks_2_3_4_only():
    """Pool counts must include masks 2, 3, 4 and exclude every other mask."""
    fl_hist = _fl_hist({
        0: np.array([100]),
        1: np.array([1_000_000]),     # SPLICED-GENOMIC-SPAN TRAP — must be ignored.
        2: np.array([1000]),
        3: np.array([2000]),
        4: np.array([3000]),
        5: np.array([100]),
        6: np.array([100]),
        7: np.array([100]),
    })
    payload = _payload(fl_hist)
    scan = _scan_trained_with_global(_hist_peaked_at(300, 10_000))
    pool = compute_pool_fl_models(payload, scan)
    assert pool.gdna_n_fragments == 1000 + 2000 + 3000


def test_spliced_mask_does_not_leak_into_gdna():
    """Even a giant mask-1 column must not affect ``gdna``."""
    intron = _hist_peaked_at(180, 10_000)
    fl_hist = _fl_hist({1: _hist_peaked_at(800, 1_000_000), 2: intron})
    payload = _payload(fl_hist)
    scan = _scan_trained_with_global(_hist_peaked_at(300, 50_000))
    pool = compute_pool_fl_models(payload, scan)
    assert abs(pool.gdna_fl_mean - 180.0) < 5.0


# ---------------------------------------------------------------------------
# 3. Passthrough identity
# ---------------------------------------------------------------------------

def test_rna_and_global_passthrough_object_identity():
    payload = _payload(_fl_hist({2: _hist_peaked_at(180, 10_000)}))
    scan = _scan_trained_with_global(
        _hist_peaked_at(300, 10_000), rna_counts=_hist_peaked_at(250, 5_000)
    )
    pool = compute_pool_fl_models(payload, scan)
    assert pool.rna_spliced is scan.rna_model
    assert pool.global_ is scan.global_model


# ---------------------------------------------------------------------------
# 4. Annotation-gap diagnostic
# ---------------------------------------------------------------------------

def test_annotation_gap_counters_cover_non_pool_masks():
    fl_hist = _fl_hist({
        0: np.array([10]),
        1: np.array([20]),
        2: np.array([1000]),
        5: np.array([30]),
        6: np.array([40]),
        7: np.array([50]),
    })
    payload = _payload(fl_hist)
    scan = _scan_trained_with_global(_hist_peaked_at(300, 10_000))
    pool = compute_pool_fl_models(payload, scan)
    assert pool.n_pool_annotation_gap == {
        "mask_0": 10, "mask_1": 20, "mask_5": 30, "mask_6": 40, "mask_7": 50,
    }


# ---------------------------------------------------------------------------
# 5. Frozen / immutability + summary serialisation
# ---------------------------------------------------------------------------

def test_pool_fl_models_is_frozen():
    payload = _payload(_fl_hist({2: _hist_peaked_at(180, 10_000)}))
    scan = _scan_trained_with_global(_hist_peaked_at(300, 10_000))
    pool = compute_pool_fl_models(payload, scan)
    with pytest.raises(Exception):
        pool.gdna_quality = "weak"  # type: ignore[misc]


def test_summary_dict_is_json_serialisable():
    payload = _payload(_fl_hist({2: _hist_peaked_at(180, 10_000), 4: np.array([100])}))
    scan = _scan_trained_with_global(_hist_peaked_at(300, 10_000))
    pool = compute_pool_fl_models(payload, scan)
    s = pool.to_summary_dict()
    blob = json.dumps(s)
    assert "gdna_quality" in blob
    assert "n_pool_annotation_gap" in blob


def test_returned_models_are_finalized_FragmentLengthModel():
    payload = _payload(_fl_hist({2: _hist_peaked_at(180, 10_000)}))
    scan = _scan_trained_with_global(_hist_peaked_at(300, 10_000))
    pool = compute_pool_fl_models(payload, scan)
    assert isinstance(pool.gdna, FragmentLengthModel)
    # Finalized models can answer log_likelihood without raising.
    _ = pool.gdna.log_likelihood(180)
    _ = pool.rna_spliced.log_likelihood(250)
    _ = pool.global_.log_likelihood(300)
