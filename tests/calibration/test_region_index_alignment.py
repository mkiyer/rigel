"""Region geometry ↔ index/payload alignment."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.errors import CalibrationSubstrateError
from rigel.calibration.region_arrays import RegionArrays
from rigel.pipeline import _check_region_payload_alignment


def _region_arrays(mini_index):
    return RegionArrays.from_region_df(mini_index.region_df, mini_index.ref_name_to_id)


def test_region_arrays_align_with_index(mini_index):
    ra = _region_arrays(mini_index)

    # R_total == len(region_df).
    assert ra.n_regions == len(mini_index.region_df)

    # Per-ref offsets equal the grouped region counts in ref_names order.
    counts = [int((mini_index.region_df["ref_name"] == ref).sum()) for ref in mini_index.ref_names]
    expected = np.concatenate([[0], np.cumsum(counts)]).astype(np.int32)
    np.testing.assert_array_equal(ra.ref_offsets, expected)


def test_neighbour_differs_per_ref(mini_index):
    ra = _region_arrays(mini_index)
    for f in range(ra.n_refs):
        lo, hi = int(ra.ref_offsets[f]), int(ra.ref_offsets[f + 1])
        sigs = ra.signature[lo:hi]
        if len(sigs) > 1:
            assert (sigs[:-1] != sigs[1:]).all(), f"ref {f} has adjacent equal signatures"


def test_alignment_guard_accepts_matching_payload(mini_index):
    ra = _region_arrays(mini_index)
    payload = SimpleNamespace(
        r_total=ra.n_regions,
        ref_region_offsets=ra.ref_offsets.astype(np.int64),
    )
    # Should not raise.
    _check_region_payload_alignment(ra, payload)


def test_alignment_guard_rejects_count_mismatch(mini_index):
    ra = _region_arrays(mini_index)
    payload = SimpleNamespace(
        r_total=ra.n_regions + 1,
        ref_region_offsets=ra.ref_offsets.astype(np.int64),
    )
    with pytest.raises(CalibrationSubstrateError):
        _check_region_payload_alignment(ra, payload)


def test_alignment_guard_rejects_offset_mismatch(mini_index):
    ra = _region_arrays(mini_index)
    bad = ra.ref_offsets.astype(np.int64).copy()
    bad[1] += 1  # perturb an offset so it no longer matches the geometry
    payload = SimpleNamespace(r_total=ra.n_regions, ref_region_offsets=bad)
    with pytest.raises(CalibrationSubstrateError):
        _check_region_payload_alignment(ra, payload)


def test_alignment_guard_rejects_none_payload(mini_index):
    ra = _region_arrays(mini_index)
    with pytest.raises(CalibrationSubstrateError):
        _check_region_payload_alignment(ra, None)
