"""Native RegionIndex metadata carriage tests."""

from __future__ import annotations

import numpy as np
import pytest

from rigel import native
from rigel.calibration import signature as sig
from rigel.calibration.regions import SIGNATURE_SENTINEL


def _type_mask_from_signature(signature: int) -> int:
    return int(np.uint8(1) << np.uint8(2 - sig.coarse_type_from_signature(signature)))


def _region_index_arrays() -> tuple[np.ndarray, ...]:
    signatures = np.arange(sig.N_SIGNATURES, dtype=np.uint8)
    starts = np.arange(sig.N_SIGNATURES, dtype=np.int64) * 10
    ends = starts + 10
    left_signatures = np.empty(sig.N_SIGNATURES, dtype=np.uint8)
    right_signatures = np.empty(sig.N_SIGNATURES, dtype=np.uint8)
    left_signatures[0] = SIGNATURE_SENTINEL
    left_signatures[1:] = signatures[:-1]
    right_signatures[:-1] = signatures[1:]
    right_signatures[-1] = SIGNATURE_SENTINEL
    boundary_kind_left = (np.arange(sig.N_SIGNATURES, dtype=np.uint8) % 6).astype(np.uint8)
    boundary_kind_right = ((np.arange(sig.N_SIGNATURES, dtype=np.uint8) + 1) % 6).astype(np.uint8)
    type_masks = np.array([_type_mask_from_signature(int(s)) for s in signatures], dtype=np.uint8)
    strands = np.array(
        [sig.coarse_strand_from_signature(int(s)) for s in signatures], dtype=np.uint8
    )
    ref_ids = np.zeros(sig.N_SIGNATURES, dtype=np.int32)
    return (
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
    )


def test_native_region_index_carries_phase2_metadata():
    arrays = _region_index_arrays()
    index = native.RegionIndex()
    index.set(*arrays, n_refs=1)

    (
        _ref_ids,
        starts,
        ends,
        signatures,
        left_signatures,
        right_signatures,
        boundary_kind_left,
        boundary_kind_right,
        type_masks,
        strands,
    ) = arrays

    assert index.n_regions() == sig.N_SIGNATURES
    assert index.n_refs() == 1
    for rid in range(sig.N_SIGNATURES):
        assert index.start(rid) == int(starts[rid])
        assert index.end(rid) == int(ends[rid])
        assert index.signature(rid) == int(signatures[rid])
        assert index.left_signature(rid) == int(left_signatures[rid])
        assert index.right_signature(rid) == int(right_signatures[rid])
        assert index.boundary_kind_left(rid) == int(boundary_kind_left[rid])
        assert index.boundary_kind_right(rid) == int(boundary_kind_right[rid])
        assert index.type_mask(rid) == int(type_masks[rid])
        assert index.strand(rid) == int(strands[rid])

    assert index.overlap(0, 15, 35) == [1, 2, 3]


def test_native_region_index_rejects_invalid_metadata():
    arrays = list(_region_index_arrays())
    arrays[3] = arrays[3].copy()
    arrays[3][0] = 16

    index = native.RegionIndex()
    with pytest.raises(ValueError, match="signature"):
        index.set(*arrays, n_refs=1)

    arrays = list(_region_index_arrays())
    arrays[6] = arrays[6].copy()
    arrays[6][0] = 6

    index = native.RegionIndex()
    with pytest.raises(ValueError, match="boundary kind"):
        index.set(*arrays, n_refs=1)
