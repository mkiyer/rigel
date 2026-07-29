"""Region geometry ↔ index/payload alignment."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.errors import CalibrationSubstrateError
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate


def _region_arrays(mini_index):
    """The LIVE geometry — the v8 node partition the scanner deposits into (plan W1b)."""
    return RegionArrays.from_index(mini_index)


def test_region_arrays_align_with_index(mini_index):
    ra = _region_arrays(mini_index)

    # R_total == len(nodes_df).
    assert ra.n_regions == len(mini_index.nodes_df)

    # Per-ref offsets equal the grouped node counts in ref_names order.
    counts = [int((mini_index.nodes_df["ref_name"] == ref).sum()) for ref in mini_index.ref_names]
    expected = np.concatenate([[0], np.cumsum(counts)]).astype(np.int32)
    np.testing.assert_array_equal(ra.ref_offsets, expected)


def test_alignment_guard_accepts_matching_payload(mini_index):
    ra = _region_arrays(mini_index)
    payload = SimpleNamespace(
        r_total=ra.n_regions,
        ref_region_offsets=ra.ref_offsets.astype(np.int64),
    )
    # Should not raise.
    CalibrationSubstrate._check_alignment(payload, ra)


def test_alignment_guard_rejects_count_mismatch(mini_index):
    ra = _region_arrays(mini_index)
    payload = SimpleNamespace(
        r_total=ra.n_regions + 1,
        ref_region_offsets=ra.ref_offsets.astype(np.int64),
    )
    with pytest.raises(CalibrationSubstrateError):
        CalibrationSubstrate._check_alignment(payload, ra)


def test_alignment_guard_rejects_offset_mismatch(mini_index):
    ra = _region_arrays(mini_index)
    bad = ra.ref_offsets.astype(np.int64).copy()
    bad[1] += 1  # perturb an offset so it no longer matches the geometry
    payload = SimpleNamespace(r_total=ra.n_regions, ref_region_offsets=bad)
    with pytest.raises(CalibrationSubstrateError):
        CalibrationSubstrate._check_alignment(payload, ra)


def test_alignment_guard_rejects_none_payload(mini_index):
    ra = _region_arrays(mini_index)
    with pytest.raises(CalibrationSubstrateError):
        CalibrationSubstrate._check_alignment(None, ra)
