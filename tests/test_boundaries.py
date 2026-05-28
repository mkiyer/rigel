"""Tests for v6 region-boundary table construction."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration.boundaries import build_boundary_table, validate_boundary_table
from rigel.calibration.region_count_ledger import build_region_count_ledger
from rigel.calibration.scan_payload import FL_HIST_N_BINS, CalibrationScanPayload
from rigel.calibration.signature import (
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_BOUNDARY_LEFT,
    COMPARTMENT_BOUNDARY_RIGHT,
    N_CHANNELS,
    N_FL_POOLS,
    N_SIGNATURES,
    SPLICE_SPLICED,
    SPLICE_UNSPLICED,
    channel_index,
    pack_signature,
)


def _payload(region_counts: np.ndarray, signatures: np.ndarray) -> CalibrationScanPayload:
    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    for row_idx, signature in enumerate(signatures):
        signature_mass[int(signature)] += float(region_counts[row_idx].sum(dtype=np.float64))
    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    return CalibrationScanPayload(
        region_counts=region_counts,
        channel_mass=region_counts.sum(axis=0, dtype=np.float64),
        signature_mass=signature_mass,
        fl_pool_mass=fl_pool_mass,
        fl_pool_total=fl_pool_mass.sum(axis=1),
        region_contained_unspliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_boundary_left_unspliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_boundary_right_unspliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_contained_spliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_boundary_left_spliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_boundary_right_spliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_unspliced_support=np.zeros(region_counts.shape[0], dtype=np.uint64),
        region_spliced_support=np.zeros(region_counts.shape[0], dtype=np.uint64),
        n_observed=int(region_counts.sum(dtype=np.float64)),
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_excluded_strand_ambig=0,
        n_unobserved=0,
        n_unannotated_ref=0,
        n_fl_unavailable=0,
        resolver_splicing_anchor_tolerance=0,
        n_regions=region_counts.shape[0],
    )


def _build_arrays() -> tuple[RegionArrays, PayloadArrays]:
    rows = [
        ("chr1", 0, 10, pack_signature()),
        ("chr1", 10, 20, pack_signature()),
        ("chr1", 20, 30, pack_signature()),
        ("chr2", 5, 15, pack_signature()),
    ]
    region_df = pd.DataFrame(rows, columns=["ref_name", "start", "end", "signature"])
    region_df["region_id"] = np.arange(len(region_df), dtype=np.int64)
    region_df.index = region_df["region_id"].to_numpy()

    rc = np.zeros((len(rows), N_CHANNELS), dtype=np.float32)
    brp = channel_index(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    brn = channel_index(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)
    blp = channel_index(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    bln = channel_index(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)
    brsp = channel_index(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_SPLICED, CHANNEL_STRAND_POS)
    blsn = channel_index(COMPARTMENT_BOUNDARY_LEFT, SPLICE_SPLICED, CHANNEL_STRAND_NEG)

    rc[0, brp] = 1.0
    rc[0, brn] = 2.0
    rc[0, brsp] = 3.0
    rc[1, blp] = 4.0
    rc[1, bln] = 5.0
    rc[1, brp] = 6.0
    rc[2, bln] = 7.0
    rc[2, blsn] = 8.0

    region_arrays = RegionArrays.from_region_df(region_df, {"chr1": 0, "chr2": 1})
    payload_arrays = PayloadArrays.from_payload(
        _payload(rc, region_df["signature"].to_numpy()), region_arrays
    )
    return region_arrays, payload_arrays


def test_boundary_table_maps_internal_left_and_right_sides() -> None:
    region_arrays, payload_arrays = _build_arrays()
    ledger = build_region_count_ledger(payload_arrays)
    side_leff = np.array([10.0, 20.0, 30.0, 40.0], dtype=np.float64)

    boundaries = build_boundary_table(region_arrays, ledger, side_leff)
    validate_boundary_table(boundaries, region_arrays, ledger)

    assert boundaries.ref_boundary_offsets.tolist() == [0, 4, 6]
    assert boundaries.boundary_pos.tolist() == [0, 10, 20, 30, 5, 15]
    assert boundaries.is_terminal.tolist() == [True, False, False, True, True, True]

    assert boundaries.left_region_unspliced_pos[1] == pytest.approx(1.0)
    assert boundaries.left_region_unspliced_neg[1] == pytest.approx(2.0)
    assert boundaries.left_region_spliced_pos[1] == pytest.approx(3.0)
    assert boundaries.right_region_unspliced_pos[1] == pytest.approx(4.0)
    assert boundaries.right_region_unspliced_neg[1] == pytest.approx(5.0)

    assert boundaries.left_region_unspliced_pos[2] == pytest.approx(6.0)
    assert boundaries.right_region_unspliced_neg[2] == pytest.approx(7.0)
    assert boundaries.right_region_spliced_neg[2] == pytest.approx(8.0)
    assert boundaries.left_region_boundary_leff.tolist() == [0.0, 10.0, 20.0, 0.0, 0.0, 0.0]
    assert boundaries.right_region_boundary_leff.tolist() == [0.0, 20.0, 30.0, 0.0, 0.0, 0.0]


def test_terminal_boundaries_are_zero_filled() -> None:
    region_arrays, payload_arrays = _build_arrays()
    ledger = build_region_count_ledger(payload_arrays)
    boundaries = build_boundary_table(region_arrays, ledger, np.ones(4, dtype=np.float64))

    terminals = np.where(boundaries.is_terminal)[0]
    for arr in (
        boundaries.left_region_unspliced_pos,
        boundaries.left_region_unspliced_neg,
        boundaries.left_region_spliced_pos,
        boundaries.left_region_spliced_neg,
        boundaries.right_region_unspliced_pos,
        boundaries.right_region_unspliced_neg,
        boundaries.right_region_spliced_pos,
        boundaries.right_region_spliced_neg,
        boundaries.left_region_boundary_leff,
        boundaries.right_region_boundary_leff,
    ):
        np.testing.assert_allclose(arr[terminals], 0.0)


def test_boundary_region_index_helpers_skip_terminals() -> None:
    region_arrays, payload_arrays = _build_arrays()
    ledger = build_region_count_ledger(payload_arrays)
    boundaries = build_boundary_table(region_arrays, ledger, np.ones(4, dtype=np.float64))

    assert boundaries.left_region_index().tolist() == [-1, 0, 1, -1, -1, -1]
    assert boundaries.right_region_index().tolist() == [-1, 1, 2, -1, -1, -1]


def test_boundary_table_rejects_noncontiguous_regions() -> None:
    rows = [
        ("chr1", 0, 10, pack_signature()),
        ("chr1", 15, 20, pack_signature()),
    ]
    region_df = pd.DataFrame(rows, columns=["ref_name", "start", "end", "signature"])
    region_df["region_id"] = np.arange(len(region_df), dtype=np.int64)
    region_df.index = region_df["region_id"].to_numpy()
    region_arrays = RegionArrays.from_region_df(region_df, {"chr1": 0})
    rc = np.zeros((2, N_CHANNELS), dtype=np.float32)
    payload_arrays = PayloadArrays.from_payload(
        _payload(rc, region_df["signature"].to_numpy()), region_arrays
    )
    ledger = build_region_count_ledger(payload_arrays)

    with pytest.raises(ValueError, match="contiguous"):
        build_boundary_table(region_arrays, ledger, np.ones(2, dtype=np.float64))