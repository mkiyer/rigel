"""Tests for pipeline wiring helpers."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd

from rigel.calibration.regions import BoundaryKind, RegionStrand, RegionType, SIGNATURE_SENTINEL
from rigel.pipeline import _wire_calibration_regions


class _ScannerSpy:
    def __init__(self) -> None:
        self.calls = []

    def set_regions(self, *args):
        self.calls.append(args)


def _region_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "ref_name": pd.array(["chrA", "chrB"], dtype="string"),
            "start": np.array([0, 0], dtype=np.int64),
            "end": np.array([100, 200], dtype=np.int64),
            "signature": np.array([0x0, 0x2], dtype=np.uint8),
            "left_signature": np.array([SIGNATURE_SENTINEL, SIGNATURE_SENTINEL], dtype=np.uint8),
            "right_signature": np.array([SIGNATURE_SENTINEL, SIGNATURE_SENTINEL], dtype=np.uint8),
            "boundary_kind_left": np.array(
                [int(BoundaryKind.NONE), int(BoundaryKind.NONE)], dtype=np.uint8
            ),
            "boundary_kind_right": np.array(
                [int(BoundaryKind.NONE), int(BoundaryKind.NONE)], dtype=np.uint8
            ),
            "type": np.array([int(RegionType.INTERGENIC), int(RegionType.EXON)], dtype=np.uint8),
            "strand": np.array([int(RegionStrand.NONE), int(RegionStrand.POS)], dtype=np.uint8),
        }
    )


def test_wire_calibration_regions_uses_index_ref_ids_not_resolver_map():
    scanner = _ScannerSpy()
    index = SimpleNamespace(ref_name_to_id={"chrA": 0, "chrB": 1})

    _wire_calibration_regions(scanner, index, _region_df())

    assert len(scanner.calls) == 1
    (
        ref_ids,
        starts,
        ends,
        signatures,
        left_signatures,
        right_signatures,
        boundary_kind_left,
        boundary_kind_right,
        type_masks,
        strands,
        n_refs,
        splicing_anchor_tolerance,
    ) = scanner.calls[0]
    np.testing.assert_array_equal(ref_ids, np.array([0, 1], dtype=np.int32))
    np.testing.assert_array_equal(starts, np.array([0, 0], dtype=np.int64))
    np.testing.assert_array_equal(ends, np.array([100, 200], dtype=np.int64))
    np.testing.assert_array_equal(signatures, np.array([0x0, 0x2], dtype=np.uint8))
    np.testing.assert_array_equal(
        left_signatures, np.array([SIGNATURE_SENTINEL, SIGNATURE_SENTINEL], dtype=np.uint8)
    )
    np.testing.assert_array_equal(
        right_signatures, np.array([SIGNATURE_SENTINEL, SIGNATURE_SENTINEL], dtype=np.uint8)
    )
    np.testing.assert_array_equal(boundary_kind_left, np.array([0, 0], dtype=np.uint8))
    np.testing.assert_array_equal(boundary_kind_right, np.array([0, 0], dtype=np.uint8))
    np.testing.assert_array_equal(type_masks, np.array([0b100, 0b001], dtype=np.uint8))
    np.testing.assert_array_equal(strands, np.array([0, 1], dtype=np.uint8))
    assert n_refs == 2
    assert splicing_anchor_tolerance == 0


def test_wire_calibration_regions_drops_unknown_refs(caplog):
    scanner = _ScannerSpy()
    index = SimpleNamespace(ref_name_to_id={"chrA": 0})

    _wire_calibration_regions(scanner, index, _region_df())

    assert len(scanner.calls) == 1
    (
        ref_ids,
        starts,
        ends,
        signatures,
        left_signatures,
        right_signatures,
        boundary_kind_left,
        boundary_kind_right,
        type_masks,
        strands,
        n_refs,
        _k,
    ) = scanner.calls[0]
    np.testing.assert_array_equal(ref_ids, np.array([0], dtype=np.int32))
    np.testing.assert_array_equal(starts, np.array([0], dtype=np.int64))
    np.testing.assert_array_equal(ends, np.array([100], dtype=np.int64))
    np.testing.assert_array_equal(signatures, np.array([0x0], dtype=np.uint8))
    np.testing.assert_array_equal(left_signatures, np.array([SIGNATURE_SENTINEL], dtype=np.uint8))
    np.testing.assert_array_equal(right_signatures, np.array([SIGNATURE_SENTINEL], dtype=np.uint8))
    np.testing.assert_array_equal(boundary_kind_left, np.array([0], dtype=np.uint8))
    np.testing.assert_array_equal(boundary_kind_right, np.array([0], dtype=np.uint8))
    np.testing.assert_array_equal(type_masks, np.array([0b100], dtype=np.uint8))
    np.testing.assert_array_equal(strands, np.array([0], dtype=np.uint8))
    assert n_refs == 1
    assert "Dropping 1 calibration regions" in caplog.text
