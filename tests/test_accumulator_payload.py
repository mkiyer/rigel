"""Schema invariants for :class:`rigel.scan_payload.AccumulatorPayload`."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.scan_payload import AccumulatorPayload, N_CHANNELS


def _toy_calibration_dict(
    boundaries_per_ref: list[list[int]],
) -> dict:
    """Build a flat calibration-dict that mirrors the C++ scanner ABI.

    ``boundaries_per_ref[f]`` is the sorted boundary positions for ref ``f``.
    A ref with zero positions contributes zero regions and zero boundaries.
    A ref with ``k+1`` positions contributes ``k`` regions and ``k+1``
    boundary objects.
    """
    n_refs = len(boundaries_per_ref)

    flat_positions: list[int] = []
    ref_pos_offsets = np.zeros(n_refs + 1, dtype=np.int64)
    ref_region_offsets = np.zeros(n_refs + 1, dtype=np.int64)
    ref_boundary_offsets = np.zeros(n_refs + 1, dtype=np.int64)

    for f, positions in enumerate(boundaries_per_ref):
        flat_positions.extend(int(p) for p in positions)
        ref_pos_offsets[f + 1] = ref_pos_offsets[f] + len(positions)
        if len(positions) == 0:
            ref_region_offsets[f + 1] = ref_region_offsets[f]
            ref_boundary_offsets[f + 1] = ref_boundary_offsets[f]
        else:
            n_regions = len(positions) - 1
            ref_region_offsets[f + 1] = ref_region_offsets[f] + n_regions
            ref_boundary_offsets[f + 1] = ref_boundary_offsets[f] + len(positions)

    r_total = int(ref_region_offsets[-1])
    b_total = int(ref_boundary_offsets[-1])
    return {
        "boundaries": np.asarray(flat_positions, dtype=np.int64),
        "ref_pos_offsets": ref_pos_offsets,
        "ref_region_offsets": ref_region_offsets,
        "ref_boundary_offsets": ref_boundary_offsets,
        "region_contained": np.zeros(r_total * N_CHANNELS, dtype=np.uint32),
        "boundary_mass_left": np.zeros(b_total * N_CHANNELS, dtype=np.float32),
        "boundary_mass_right": np.zeros(b_total * N_CHANNELS, dtype=np.float32),
        "boundary_flux_left": np.zeros(b_total * N_CHANNELS, dtype=np.uint32),
        "boundary_flux_right": np.zeros(b_total * N_CHANNELS, dtype=np.uint32),
        "n_channels": N_CHANNELS,
        "n_refs": n_refs,
    }


class TestAccumulatorPayloadFromScanResult:
    def test_single_ref_three_regions(self):
        cal = _toy_calibration_dict([[0, 100, 250, 500]])
        payload = AccumulatorPayload.from_scan_result({"calibration": cal})

        assert payload.n_refs == 1
        assert payload.n_channels == N_CHANNELS
        assert payload.r_total == 3
        assert payload.b_obj_total == 4
        assert payload.region_contained.shape == (3, N_CHANNELS)
        assert payload.boundary_mass_left.shape == (4, N_CHANNELS)
        assert payload.boundary_mass_right.shape == (4, N_CHANNELS)
        assert payload.boundary_flux_left.shape == (4, N_CHANNELS)
        assert payload.boundary_flux_right.shape == (4, N_CHANNELS)
        assert payload.region_contained.dtype == np.uint32
        assert payload.boundary_mass_left.dtype == np.float32
        assert payload.boundary_mass_right.dtype == np.float32
        assert payload.boundary_flux_left.dtype == np.uint32
        assert payload.boundary_flux_right.dtype == np.uint32

    def test_multi_ref_mixed_zero(self):
        # ref0: 2 regions, ref1: empty, ref2: 1 region.
        cal = _toy_calibration_dict([[0, 50, 200], [], [0, 1000]])
        payload = AccumulatorPayload.from_scan_result({"calibration": cal})

        assert payload.n_refs == 3
        assert payload.r_total == 3  # 2 + 0 + 1
        assert payload.b_obj_total == 5  # 3 + 0 + 2
        np.testing.assert_array_equal(payload.ref_pos_offsets, [0, 3, 3, 5])
        np.testing.assert_array_equal(payload.ref_region_offsets, [0, 2, 2, 3])
        np.testing.assert_array_equal(payload.ref_boundary_offsets, [0, 3, 3, 5])

    def test_all_refs_empty(self):
        cal = _toy_calibration_dict([[], [], []])
        payload = AccumulatorPayload.from_scan_result({"calibration": cal})
        assert payload.r_total == 0
        assert payload.b_obj_total == 0
        assert payload.boundaries.size == 0

    def test_contiguity(self):
        cal = _toy_calibration_dict([[0, 10, 20]])
        payload = AccumulatorPayload.from_scan_result({"calibration": cal})
        for arr in (
            payload.boundaries,
            payload.ref_pos_offsets,
            payload.ref_region_offsets,
            payload.ref_boundary_offsets,
            payload.region_contained,
            payload.boundary_mass_left,
            payload.boundary_mass_right,
            payload.boundary_flux_left,
            payload.boundary_flux_right,
        ):
            assert arr.flags["C_CONTIGUOUS"]

    def test_immutable_dataclass(self):
        cal = _toy_calibration_dict([[0, 10]])
        payload = AccumulatorPayload.from_scan_result({"calibration": cal})
        with pytest.raises((AttributeError, TypeError)):
            payload.n_refs = 99  # frozen

    def test_missing_calibration_raises(self):
        with pytest.raises(ValueError, match="calibration.*None"):
            AccumulatorPayload.from_scan_result({"calibration": None})

    def test_wrong_n_channels_raises(self):
        cal = _toy_calibration_dict([[0, 10]])
        cal["n_channels"] = 8
        with pytest.raises(ValueError, match="n_channels"):
            AccumulatorPayload.from_scan_result({"calibration": cal})

    def test_pos_offsets_length_mismatch_raises(self):
        cal = _toy_calibration_dict([[0, 10]])
        cal["ref_pos_offsets"] = np.array([0, 1, 2], dtype=np.int64)
        with pytest.raises(ValueError, match="ref_pos_offsets"):
            AccumulatorPayload.from_scan_result({"calibration": cal})
