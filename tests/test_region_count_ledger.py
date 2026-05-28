"""Tests for `rigel.calibration.region_count_ledger` (v4 Phase 1 §5)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration.region_count_ledger import (
    RegionCountLedger,
    build_region_count_ledger,
)
from rigel.calibration.scan_payload import FL_HIST_N_BINS, CalibrationScanPayload
from rigel.calibration.signature import (
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_BOUNDARY_LEFT,
    COMPARTMENT_BOUNDARY_RIGHT,
    COMPARTMENT_CONTAINED,
    N_CHANNELS,
    N_FL_POOLS,
    N_SIGNATURES,
    SPLICE_SPLICED,
    SPLICE_UNSPLICED,
    channel_index,
    pack_signature,
)


# ---- fixtures --------------------------------------------------------------


def _region_df() -> pd.DataFrame:
    rows = [
        ("chr1", 0, 500, pack_signature()),
        ("chr1", 500, 1000, pack_signature(exon_pos=True)),
        ("chr1", 1000, 2000, pack_signature(intron_pos=True)),
        ("chr1", 2000, 2600, pack_signature(intron_pos=True, intron_neg=True)),
        ("chr1", 2600, 3000, pack_signature(intron_pos=True, exon_pos=True)),
    ]
    df = pd.DataFrame(rows, columns=["ref_name", "start", "end", "signature"])
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    df.index = df["region_id"].to_numpy()
    return df


def _payload(region_counts: np.ndarray) -> CalibrationScanPayload:
    n_regions = region_counts.shape[0]
    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    for row, signature in enumerate(_region_df()["signature"].to_numpy()):
        signature_mass[int(signature)] += float(region_counts[row].sum(dtype=np.float64))
    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    return CalibrationScanPayload(
        region_counts=region_counts,
        channel_mass=region_counts.sum(axis=0, dtype=np.float64),
        signature_mass=signature_mass,
        fl_pool_mass=fl_pool_mass,
        fl_pool_total=fl_pool_mass.sum(axis=1),
        region_contained_unspliced_support=np.zeros(n_regions, dtype=np.uint32),
        region_boundary_left_unspliced_support=np.zeros(n_regions, dtype=np.uint32),
        region_boundary_right_unspliced_support=np.zeros(n_regions, dtype=np.uint32),
        region_contained_spliced_support=np.zeros(n_regions, dtype=np.uint32),
        region_boundary_left_spliced_support=np.zeros(n_regions, dtype=np.uint32),
        region_boundary_right_spliced_support=np.zeros(n_regions, dtype=np.uint32),
        region_unspliced_support=np.zeros(n_regions, dtype=np.uint64),
        region_spliced_support=np.zeros(n_regions, dtype=np.uint64),
        n_observed=int(region_counts.sum(dtype=np.float64)),
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_excluded_strand_ambig=0,
        n_unobserved=0,
        n_unannotated_ref=0,
        n_fl_unavailable=0,
        resolver_splicing_anchor_tolerance=0,
        n_regions=n_regions,
    )


def _arrays() -> tuple[RegionArrays, PayloadArrays, np.ndarray]:
    region_df = _region_df()
    region_arrays = RegionArrays.from_region_df(region_df, {"chr1": 0})
    rng = np.random.default_rng(20260524)
    rc = rng.integers(0, 50, size=(len(region_df), N_CHANNELS)).astype(np.float32)
    payload = _payload(rc)
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    return region_arrays, payload_arrays, payload_arrays.region_counts_sorted


# ---- tests -----------------------------------------------------------------


def test_ledger_views_match_payload_channel_slices() -> None:
    _, payload_arrays, rc = _arrays()
    ledger = build_region_count_ledger(payload_arrays)
    assert isinstance(ledger, RegionCountLedger)
    np.testing.assert_array_equal(
        ledger.contained_unspliced_pos,
        rc[:, channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)],
    )
    np.testing.assert_array_equal(
        ledger.contained_unspliced_neg,
        rc[:, channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)],
    )
    np.testing.assert_array_equal(
        ledger.boundary_left_spliced_pos,
        rc[:, channel_index(COMPARTMENT_BOUNDARY_LEFT, SPLICE_SPLICED, CHANNEL_STRAND_POS)],
    )
    np.testing.assert_array_equal(
        ledger.boundary_right_spliced_neg,
        rc[:, channel_index(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_SPLICED, CHANNEL_STRAND_NEG)],
    )


def test_ledger_holds_views_not_copies() -> None:
    _, payload_arrays, _ = _arrays()
    ledger = build_region_count_ledger(payload_arrays)
    rc = payload_arrays.region_counts_sorted
    for arr in (
        ledger.contained_unspliced_pos,
        ledger.contained_unspliced_neg,
        ledger.boundary_left_unspliced_pos,
        ledger.boundary_left_unspliced_neg,
        ledger.boundary_right_unspliced_pos,
        ledger.boundary_right_unspliced_neg,
        ledger.contained_spliced_pos,
        ledger.contained_spliced_neg,
        ledger.boundary_left_spliced_pos,
        ledger.boundary_left_spliced_neg,
        ledger.boundary_right_spliced_pos,
        ledger.boundary_right_spliced_neg,
    ):
        assert np.shares_memory(arr, rc)


def test_unspliced_totals_match_direct_sums() -> None:
    _, payload_arrays, rc = _arrays()
    ledger = build_region_count_ledger(payload_arrays)

    def chan(c: int, s: int, k: int) -> np.ndarray:
        return rc[:, channel_index(c, s, k)]

    np.testing.assert_array_equal(
        ledger.contained_unspliced_total(),
        chan(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
        + chan(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG),
    )
    np.testing.assert_array_equal(
        ledger.boundary_left_unspliced_total(),
        chan(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
        + chan(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG),
    )
    np.testing.assert_array_equal(
        ledger.boundary_right_unspliced_total(),
        chan(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
        + chan(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG),
    )
    np.testing.assert_array_equal(
        ledger.boundary_unspliced_total(),
        ledger.boundary_left_unspliced_total() + ledger.boundary_right_unspliced_total(),
    )


def test_spliced_totals_match_direct_sums() -> None:
    _, payload_arrays, rc = _arrays()
    ledger = build_region_count_ledger(payload_arrays)

    def chan(c: int, s: int, k: int) -> np.ndarray:
        return rc[:, channel_index(c, s, k)]

    contained_spliced = chan(
        COMPARTMENT_CONTAINED, SPLICE_SPLICED, CHANNEL_STRAND_POS
    ) + chan(COMPARTMENT_CONTAINED, SPLICE_SPLICED, CHANNEL_STRAND_NEG)
    boundary_spliced = (
        chan(COMPARTMENT_BOUNDARY_LEFT, SPLICE_SPLICED, CHANNEL_STRAND_POS)
        + chan(COMPARTMENT_BOUNDARY_LEFT, SPLICE_SPLICED, CHANNEL_STRAND_NEG)
        + chan(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_SPLICED, CHANNEL_STRAND_POS)
        + chan(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_SPLICED, CHANNEL_STRAND_NEG)
    )

    np.testing.assert_array_equal(ledger.contained_spliced_total(), contained_spliced)
    np.testing.assert_array_equal(ledger.boundary_spliced_total(), boundary_spliced)
    np.testing.assert_array_equal(
        ledger.spliced_total(), contained_spliced + boundary_spliced
    )


def test_payload_arrays_no_longer_carries_totals() -> None:
    _, payload_arrays, _ = _arrays()
    for name in (
        "contained_unspliced_total",
        "boundary_left_unspliced_total",
        "boundary_right_unspliced_total",
    ):
        assert not hasattr(payload_arrays, name), f"PayloadArrays still exposes {name!r}"


def test_build_rejects_non_2d_region_counts() -> None:
    bad = PayloadArrays.__new__(PayloadArrays)
    object.__setattr__(bad, "region_counts_sorted", np.zeros(5, dtype=np.float32))
    with pytest.raises(ValueError, match="region_counts_sorted must be 2D"):
        build_region_count_ledger(bad)
