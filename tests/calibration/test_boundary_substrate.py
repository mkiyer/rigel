"""BoundarySubstrate — the boundary-node view of the payload (PLAN v6 §2).

The boundary nodes carry the SAME information as ``CalibrationSubstrate.left/right`` (which project the
boundary sides onto regions), re-keyed by boundary. These tests pin that re-indexing identity — the Phase-A
gate that the bipartite boundary substrate is a faithful, lossless re-view (no payload reshape, no C++ change).
"""

from __future__ import annotations

import numpy as np

from _synthetic import make_synthetic_payload

from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate


def _views_equal(a: "object", i: int, b: "object", j: int) -> bool:
    """Row-wise equality of two SubstrateView rows across all six fields."""
    return all(
        np.array_equal(getattr(a, f)[i], getattr(b, f)[j])
        for f in (
            "n_unspliced_pos",
            "n_unspliced_neg",
            "n_spliced_sense",
            "n_spliced_antisense",
            "mass_unspliced",
            "mass_spliced",
        )
    )


def test_boundary_adjacency_indices():
    payload, _ = make_synthetic_payload()
    bsub = BoundarySubstrate.from_payload(payload)
    # 3 regions ⇒ 4 boundaries (2 terminal + 2 internal); terminals have one -1 flank.
    assert bsub.n_boundaries == 4
    np.testing.assert_array_equal(bsub.left_region, [-1, 0, 1, 2])
    np.testing.assert_array_equal(bsub.right_region, [0, 1, 2, -1])


def test_boundary_sides_read_payload_channels():
    payload, _ = make_synthetic_payload()
    bsub = BoundarySubstrate.from_payload(payload)
    # b1 LEFT side = boundary_flux_left[1]=[4,1,0,0] / mass_left[1]=[1.5,0.5,0,0]
    assert bsub.left.n_unspliced_pos[1] == 4 and bsub.left.n_unspliced_neg[1] == 1
    np.testing.assert_allclose(bsub.left.mass_unspliced[1], 2.0)
    # b1 RIGHT side = boundary_flux_right[1]=[6,2,0,0]
    assert bsub.right.n_unspliced_pos[1] == 6 and bsub.right.n_unspliced_neg[1] == 2
    # b2 spliced channels read correctly (sense on the left side, antisense on the right side).
    assert bsub.left.n_spliced_sense[2] == 1 and bsub.right.n_spliced_antisense[2] == 1


def test_reindexing_identity_vs_region_substrate():
    """The load-bearing Phase-A gate: a boundary's two sides equal the flanking regions' projected views.

    For an internal boundary ``b`` with ``(left_region=lr, right_region=rr)``:
        BoundarySubstrate.left[b]  == CalibrationSubstrate.right[lr]   (the side inside the left region)
        BoundarySubstrate.right[b] == CalibrationSubstrate.left[rr]    (the side inside the right region)
    """
    payload, ra = make_synthetic_payload()
    rsub = CalibrationSubstrate.from_payload(payload, ra)
    bsub = BoundarySubstrate.from_payload(payload)

    checked = 0
    for b in range(bsub.n_boundaries):
        lr, rr = int(bsub.left_region[b]), int(bsub.right_region[b])
        if lr >= 0:
            assert _views_equal(bsub.left, b, rsub.right, lr), f"boundary {b} left ≠ region {lr} right"
            checked += 1
        if rr >= 0:
            assert _views_equal(bsub.right, b, rsub.left, rr), f"boundary {b} right ≠ region {rr} left"
            checked += 1
    # 2 internal boundaries × 2 sides + 2 terminals × 1 inner side = 6 checked adjacencies.
    assert checked == 6


def test_terminal_offedge_side_is_empty():
    payload, _ = make_synthetic_payload()
    bsub = BoundarySubstrate.from_payload(payload)
    # b0 is a reference-start terminal (left_region == -1): its LEFT (off-edge) side carries no mass.
    assert bsub.left_region[0] == -1
    assert bsub.left.mass_unspliced[0] == 0.0 and bsub.left.n_unspliced[0] == 0
    # b3 is a reference-end terminal (right_region == -1): its RIGHT (off-edge) side carries no mass.
    assert bsub.right_region[3] == -1
    assert bsub.right.mass_unspliced[3] == 0.0 and bsub.right.n_unspliced[3] == 0
