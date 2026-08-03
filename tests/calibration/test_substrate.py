"""CalibrationSubstrate — ONE type, the payload's own axes, and the fixed point decoded exactly once.

        (S5.d)

⛔ **What this replaces.** ``CalibrationSubstrate`` held three per-REGION views (contained / left / right)
and ``BoundarySubstrate`` held the same numbers re-keyed by boundary. Two classes, one set of numbers,
two keyings — and they existed solely because a boundary had two sides. A contiguous edge is a 0-bp line
with ONE set of numbers, so both the second class and the whole left/right axis dissolve.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.errors import CalibrationSubstrateError
from rigel.calibration.substrate import CalibrationSubstrate, PopulationView

from _synthetic import INV_LENGTH_SCALE, make_synthetic_payload


@pytest.fixture
def substrate():
    payload, region_arrays = make_synthetic_payload()
    return CalibrationSubstrate.from_payload(payload, region_arrays), payload, region_arrays


# ---------------------------------------------------------------------------
# the five populations, on their own axes
# ---------------------------------------------------------------------------


def test_every_population_is_present_on_the_RIGHT_axis(substrate):
    """Nodes, contiguous edges and junctions are three axes off by one per reference.

    A population read against the wrong axis is the defect class that once dropped 476,719 of 476,732
    fragments while every golden test passed, so the shapes are asserted rather than assumed.
    """
    sub, payload, _ = substrate
    for view, n in (
        (sub.node_contained, payload.n_nodes),
        (sub.node_spanning, payload.n_nodes),
        (sub.edge_unspliced, payload.n_edges),
        (sub.edge_spliced, payload.n_edges),
        (sub.junction, payload.n_sj),
    ):
        assert view.count.shape == (n, 2)
        assert view.inv_length_sum.shape == (n, 2)
        assert view.length_sum.shape == (n, 2)
    assert payload.n_nodes != payload.n_edges, "the fixture must not let an axis mix-up pass"


def test_the_columns_are_GENOME_STRAND_and_nothing_is_re_oriented(substrate):
    """⭐ ONE convention. Sense/antisense is transcript-relative, derived by the consumer from a
    junction's own strand, and never stored — so no field here may be named for it."""
    sub, payload, _ = substrate
    np.testing.assert_array_equal(sub.node_contained.count, payload.node_contained_count)
    for name in dir(sub):
        assert "sense" not in name, (
            f"{name} names a transcript-relative concept the schema does not store"
        )


def test_no_population_is_a_VIEW_OF_ANOTHER(substrate):
    """The old substrate's `left`/`right` were the same numbers twice. Nothing here may alias."""
    sub, _, _ = substrate
    banks = [sub.node_contained, sub.node_spanning, sub.edge_unspliced, sub.edge_spliced]
    totals = [int(b.count.sum()) for b in banks]
    assert len(set(totals)) == len(totals), "the fixture gives every bank a distinct total"


# ---------------------------------------------------------------------------
# the fixed point
# ---------------------------------------------------------------------------


def test_the_fixed_point_is_DECODED_HERE_and_only_here(substrate):
    """⭐ ``inv_length_sum`` leaves the payload as ``round(2^32 / placements)`` and must arrive at a
    consumer as a real number. Decoding at the boundary is what stops every downstream module needing to
    know about ``2^32`` — one place to be wrong instead of many."""
    sub, payload, _ = substrate
    np.testing.assert_allclose(
        sub.node_contained.inv_length_sum,
        payload.node_contained_inv_length_sum.astype(np.float64) / INV_LENGTH_SCALE,
        rtol=0,
        atol=0,
    )
    assert sub.node_contained.inv_length_sum.dtype == np.float64


def test_a_decoded_sum_recovers_the_reciprocal_placements_it_was_built_from(substrate):
    """The fixture deposited ``n`` fragments at 50 placements into the contained bank, so the decoded
    sum must read ``n / 50`` — the round trip, not just a division."""
    sub, payload, _ = substrate
    counts = payload.node_contained_count.astype(np.float64)
    # ⚠ rtol above the fixed point's own quantisation, which the spec bounds at 6.9e-8 relative over
    # L in [40, 1000] — asserting tighter would be asserting the rounding, not the decode.
    np.testing.assert_allclose(sub.node_contained.inv_length_sum, counts / 50.0, rtol=1e-7)


def test_length_sum_is_NOT_scaled(substrate):
    """It is a plain sum of integers — no fixed point, no scale. Applying the density decode to it would
    divide by 4.3 billion and read as zero everywhere."""
    sub, payload, _ = substrate
    np.testing.assert_array_equal(
        sub.node_contained.length_sum, payload.node_contained_length_sum.astype(np.float64)
    )


# ---------------------------------------------------------------------------
# derived quantities
# ---------------------------------------------------------------------------


def test_mean_length_is_length_sum_over_count(substrate):
    """The quantity `length_sum` exists for: an object's own mean fragment length, model-free."""
    sub, _, _ = substrate
    view = sub.node_contained
    expected = view.length_sum.sum(axis=1) / view.count.sum(axis=1)
    np.testing.assert_allclose(view.mean_length, expected)
    assert np.allclose(view.mean_length, 50.0), "the fixture deposited only 50 bp fragments here"


def test_mean_length_at_ZERO_count_is_NAN_not_zero():
    """⛔ An object with no fragments has no mean length, and that is not 0.

    an object with no opportunity must emit NOTHING, never a floored
    value — a "no data" default of 100 % gDNA was actively seeding false gDNA into neighbouring exons.
    Zero here would read as "the fragments here are infinitely short".
    """
    view = PopulationView(
        count=np.zeros((2, 2), np.int64),
        inv_length_sum=np.zeros((2, 2), np.float64),
        length_sum=np.zeros((2, 2), np.float64),
    )
    assert np.all(np.isnan(view.mean_length))


def test_total_count_sums_the_two_strands(substrate):
    sub, payload, _ = substrate
    np.testing.assert_array_equal(
        sub.node_contained.total_count, payload.node_contained_count.sum(axis=1)
    )


# ---------------------------------------------------------------------------
# the alignment check — the one piece of the old substrate that survives
# ---------------------------------------------------------------------------


def test_a_geometry_with_the_WRONG_OBJECT_COUNT_is_refused():
    """Load-bearing, and kept verbatim in spirit from the old substrate: a geometry/payload mismatch
    otherwise surfaces as a shape error deep in the solver, pointing nowhere near its cause."""
    payload, region_arrays = make_synthetic_payload()
    import dataclasses

    trimmed = dataclasses.replace(
        region_arrays,
        start=region_arrays.start[:2],
        end=region_arrays.end[:2],
        ref_id=region_arrays.ref_id[:2],
        signature=region_arrays.signature[:2],
        strand_class=region_arrays.strand_class[:2],
        region_size_bp=region_arrays.region_size_bp[:2],
    )
    with pytest.raises(CalibrationSubstrateError, match="2 objects but the payload has 3"):
        CalibrationSubstrate.from_payload(payload, trimmed)


def test_a_geometry_with_DRIFTED_PER_REFERENCE_OFFSETS_is_refused():
    """⛔ The count can match while the per-reference slicing does not. That is the exact defect that
    once dropped 476,719 of 476,732 real fragments inside deposit() with every golden test green, so
    matching totals are NOT sufficient evidence."""
    import dataclasses

    payload, region_arrays = make_synthetic_payload()
    drifted = dataclasses.replace(region_arrays, ref_offsets=np.array([0, 2], dtype=np.int32))
    with pytest.raises(CalibrationSubstrateError, match="per-reference offsets"):
        CalibrationSubstrate.from_payload(payload, drifted)


def test_a_None_payload_is_refused_by_NAME():
    with pytest.raises(CalibrationSubstrateError, match="set_regions"):
        CalibrationSubstrate.from_payload(None, make_synthetic_payload()[1])


# ---------------------------------------------------------------------------
# what died
# ---------------------------------------------------------------------------


def test_BoundarySubstrate_and_the_left_right_axis_are_GONE():
    """Two classes holding one set of numbers in two keyings existed only because a boundary had two
    sides. An edge does not, so both go — and with them ``_make_view`` and the re-keying identity."""
    from rigel.calibration import substrate as mod

    assert not hasattr(mod, "BoundarySubstrate")
    assert not hasattr(mod, "SubstrateView")
    assert not hasattr(mod, "_make_view")
    sub = CalibrationSubstrate.from_payload(*make_synthetic_payload())
    for dead in ("left", "right", "contained"):
        assert not hasattr(sub, dead), f"{dead} is a per-face concept and must not survive"
