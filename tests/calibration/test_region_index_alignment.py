"""Region geometry ↔ index/payload alignment."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.errors import CalibrationSubstrateError
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate


def _region_arrays(mini_index):
    """The LIVE geometry — the v8 region partition the scanner deposits into (plan W1b)."""
    return RegionArrays.from_index(mini_index)


def test_region_arrays_align_with_index(mini_index):
    ra = _region_arrays(mini_index)

    # R_total == len(regions_df).
    assert ra.n_regions == len(mini_index.regions_df)

    # Per-ref offsets equal the grouped region counts in ref_names order.
    counts = [int((mini_index.regions_df["ref_name"] == ref).sum()) for ref in mini_index.ref_names]
    expected = np.concatenate([[0], np.cumsum(counts)]).astype(np.int32)
    np.testing.assert_array_equal(ra.ref_offsets, expected)


def _payload(n_regions, ref_region_offsets):
    """The two fields the alignment guard reads. ⚠ Both are on the REGION axis: the payload's boundary and
    sj axes are sized from it (``E = N − n_refs``), so a region-axis mismatch is the one that has
    to be caught at the door."""
    return SimpleNamespace(n_regions=n_regions, ref_region_offsets=np.asarray(ref_region_offsets, np.int64))


def test_alignment_guard_accepts_matching_payload(mini_index):
    ra = _region_arrays(mini_index)
    # Should not raise.
    CalibrationSubstrate._check_alignment(_payload(ra.n_regions, ra.ref_offsets), ra)


def test_alignment_guard_rejects_count_mismatch(mini_index):
    ra = _region_arrays(mini_index)
    with pytest.raises(CalibrationSubstrateError):
        CalibrationSubstrate._check_alignment(_payload(ra.n_regions + 1, ra.ref_offsets), ra)


def test_alignment_guard_rejects_offset_mismatch(mini_index):
    ra = _region_arrays(mini_index)
    bad = ra.ref_offsets.astype(np.int64).copy()
    bad[1] += 1  # perturb an offset so it no longer matches the geometry
    with pytest.raises(CalibrationSubstrateError):
        CalibrationSubstrate._check_alignment(_payload(ra.n_regions, bad), ra)


def test_alignment_guard_rejects_none_payload(mini_index):
    ra = _region_arrays(mini_index)
    with pytest.raises(CalibrationSubstrateError):
        CalibrationSubstrate._check_alignment(None, ra)
