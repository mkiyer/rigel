"""Schema invariants for :class:`rigel.scan_payload.AccumulatorPayload`.

    Spec: ``tests/native/_accumulator_reference.py``   ·   Plan: ``docs/IMPLEMENTATION_PLAN.md`` §3.5

The payload is the boundary where the C++ tally becomes Python. Its field names **are** the
specification's ``Tally`` field names, character for character, so the tests below check the schema
against `Tally` itself rather than against a hand-written list — a list would be free to drift, and the
whole point of sharing one vocabulary is that it cannot.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rigel.scan_payload import N_STRAND_COLUMNS, AccumulatorPayload, ScanQC

from native._accumulator_reference import Tally


#: Two references: chr1 with 4 cuts (3 nodes, 2 lines) and chr2 with 3 (2 nodes, 1 line). A third
#: reference has NO cuts, which is legal and contributes nothing — the case where per-reference offset
#: arithmetic goes wrong if it is written as a plain subtraction.
CUTS_PER_REF = [[0, 100, 200, 600], [0, 500, 900], []]
MAX_LENGTH = 12


def _deposited_lengths(total: int) -> np.ndarray:
    """A uint32[MAX_LENGTH + 1] histogram summing to exactly ``total``."""
    out = np.zeros(MAX_LENGTH + 1, dtype=np.uint32)
    out[MAX_LENGTH // 2] = total
    return out


def _calibration_dict(**overrides) -> dict:
    """A flat calibration dict shaped exactly as ``BamScanner::build_result`` emits one."""
    ref_cut_offsets = np.zeros(len(CUTS_PER_REF) + 1, np.int64)
    ref_node_offsets = np.zeros(len(CUTS_PER_REF) + 1, np.int64)
    ref_edge_offsets = np.zeros(len(CUTS_PER_REF) + 1, np.int64)
    for f, cuts in enumerate(CUTS_PER_REF):
        ref_cut_offsets[f + 1] = ref_cut_offsets[f] + len(cuts)
        ref_node_offsets[f + 1] = ref_node_offsets[f] + max(len(cuts) - 1, 0)
        ref_edge_offsets[f + 1] = ref_edge_offsets[f] + max(len(cuts) - 2, 0)

    n_nodes, n_edges, n_sj = (
        int(ref_node_offsets[-1]),
        int(ref_edge_offsets[-1]),
        3,
    )
    ref_sj_offsets = np.asarray([0, 2, 3, 3], np.int64)

    cal = {
        "cut_positions": np.asarray([c for cuts in CUTS_PER_REF for c in cuts], np.int64),
        "ref_cut_offsets": ref_cut_offsets,
        "ref_node_offsets": ref_node_offsets,
        "ref_edge_offsets": ref_edge_offsets,
        "ref_sj_offsets": ref_sj_offsets,
        "node_contained_count": np.arange(n_nodes * 2, dtype=np.uint32),
        "node_contained_inv_length_sum": np.arange(n_nodes * 2, dtype=np.uint64) * 7,
        "node_contained_length_sum": np.arange(n_nodes * 2, dtype=np.uint64) * 11,
        "node_spanning_count": np.arange(n_nodes * 2, dtype=np.uint32) * 2,
        "node_spanning_inv_length_sum": np.arange(n_nodes * 2, dtype=np.uint64) * 3,
        "node_spanning_length_sum": np.arange(n_nodes * 2, dtype=np.uint64) * 13,
        "node_start_count": np.arange(n_nodes, dtype=np.uint32),
        "edge_unspliced_count": np.arange(n_edges * 2, dtype=np.uint32),
        "edge_unspliced_inv_length_sum": np.arange(n_edges * 2, dtype=np.uint64),
        "edge_unspliced_length_sum": np.arange(n_edges * 2, dtype=np.uint64) * 17,
        "edge_spliced_count": np.arange(n_edges * 2, dtype=np.uint32),
        "edge_spliced_inv_length_sum": np.arange(n_edges * 2, dtype=np.uint64),
        "edge_spliced_length_sum": np.arange(n_edges * 2, dtype=np.uint64) * 19,
        "sj_count": np.arange(n_sj * 2, dtype=np.uint32),
        "sj_inv_length_sum": np.arange(n_sj * 2, dtype=np.uint64),
        "sj_length_sum": np.arange(n_sj * 2, dtype=np.uint64) * 23,
        "pool_lengths": np.arange(5 * (MAX_LENGTH + 1), dtype=np.int64),
        # ⭐ C1: the unconditional histogram must bin EXACTLY the deposited fragments, so this fixture
        # can no longer carry an arbitrary array — 41 here, matching qc.deposited below. That coupling is
        # the invariant doing its job at the door.
        "deposited_lengths": _deposited_lengths(41),
        "qc": {
            "deposited": 41,
            "dropped_too_long": 1,
            "dropped_empty": 2,
            "dropped_strand_undefined": 3,
            "dropped_ambiguous_path": 4,
            "unannotated_introns": 5,
            "contradictory_sj_strand": 6,
            "sj_implicit_fragments": 7,
            "introns_absorbed": 8,
        },
        "n_strand_columns": N_STRAND_COLUMNS,
        "n_fragment_pools": 5,
        "max_length": MAX_LENGTH,
        "n_refs": len(CUTS_PER_REF),
    }
    cal.update(overrides)
    return {"calibration": cal}


def _payload(**overrides) -> AccumulatorPayload:
    return AccumulatorPayload.from_scan_result(_calibration_dict(**overrides))


# ---------------------------------------------------------------------------
# the schema IS the specification's Tally
# ---------------------------------------------------------------------------


def test_the_payload_carries_every_field_of_the_specifications_Tally():
    """⛔ Read off ``Tally``, never written out here.

    The payload, the reference and the parity gate share one vocabulary precisely so that no mapping
    table exists to drift. A hand-written list in this test would be that table.
    """
    payload = _payload()
    for field in dataclasses.fields(Tally):
        assert hasattr(payload, field.name), (
            f"the payload has no {field.name!r}; it is a field of the specification's Tally, so either "
            f"the payload dropped it or the two vocabularies have diverged"
        )


def test_the_two_column_banks_are_reshaped_and_the_one_column_ones_are_not():
    payload = _payload()
    n_nodes, n_edges, n_sj = payload.n_nodes, payload.n_edges, payload.n_sj
    assert (n_nodes, n_edges, n_sj) == (5, 3, 3)
    for name in (
        "node_contained_count",
        "node_contained_inv_length_sum",
        "node_contained_length_sum",
        "node_spanning_count",
        "node_spanning_inv_length_sum",
        "node_spanning_length_sum",
    ):
        assert getattr(payload, name).shape == (n_nodes, N_STRAND_COLUMNS), name
    for name in (
        "edge_unspliced_count",
        "edge_unspliced_inv_length_sum",
        "edge_unspliced_length_sum",
        "edge_spliced_count",
        "edge_spliced_inv_length_sum",
        "edge_spliced_length_sum",
    ):
        assert getattr(payload, name).shape == (n_edges, N_STRAND_COLUMNS), name
    for name in ("sj_count", "sj_inv_length_sum", "sj_length_sum"):
        assert getattr(payload, name).shape == (n_sj, N_STRAND_COLUMNS), name
    assert payload.node_start_count.shape == (n_nodes,)
    assert payload.pool_lengths.shape == (5, MAX_LENGTH + 1)


def test_the_dtypes_are_the_specifications_dtypes():
    """⚠ Counts are uint32 and densities uint64, and the payload must not silently widen either.

    A count that arrives as int64 compares equal to the specification's uint32 by value, so a value-only
    check would pass while the schema had changed underneath it.
    """
    payload = _payload()
    reference = Tally.zeros(n_nodes=5, n_edges=3, n_sj=3, max_length=MAX_LENGTH)
    for field in dataclasses.fields(Tally):
        if field.name == "qc":
            continue
        expected = getattr(reference, field.name).dtype
        assert getattr(payload, field.name).dtype == expected, field.name


def test_a_WRONG_dtype_is_REJECTED_rather_than_coerced():
    """⛔ Checking the output dtype is not enough, and a perturbation proved it.

    ``ascontiguousarray(x, dtype=uint32)`` will happily narrow an int64 array, so a payload that *coerces*
    still reports the right dtype and passes a check on its own output. What that hides is a C++ side that
    has stopped agreeing with the schema — the dtype is part of byte-identity, and a widened count compares
    equal by value all the way down. So the wrong dtype has to arrive and be refused.
    """
    with pytest.raises(ValueError, match="node_contained_count.*dtype"):
        _payload(node_contained_count=np.arange(10, dtype=np.int64))
    with pytest.raises(ValueError, match="sj_inv_length_sum.*dtype"):
        _payload(sj_inv_length_sum=np.arange(6, dtype=np.uint32))


def test_a_MISSING_qc_denominator_is_REJECTED():
    """⛔ Also found by perturbation: nothing was feeding an incomplete qc block.

    Design §10.3 requires every one of these to be emitted, because every conservation statement
    downstream has to be able to name what it excluded. A denominator that silently arrives absent is a
    statement that cannot.
    """
    qc = dict(_calibration_dict()["calibration"]["qc"])
    del qc["dropped_ambiguous_path"]
    with pytest.raises(ValueError, match="dropped_ambiguous_path"):
        _payload(qc=qc)


def test_a_reference_with_no_cuts_contributes_nothing_to_any_axis():
    """The third reference is empty. Per-reference offsets written as a plain subtraction go negative."""
    payload = _payload()
    assert payload.ref_node_offsets[-1] == payload.ref_node_offsets[-2]
    assert payload.ref_edge_offsets[-1] == payload.ref_edge_offsets[-2]
    assert payload.n_refs == 3


# ---------------------------------------------------------------------------
# the QC denominators
# ---------------------------------------------------------------------------


def test_qc_is_typed_so_a_misspelled_denominator_fails_loudly():
    """⚠ Design §10.3 requires these to be EMITTED, and every conservation statement to name its
    denominator. A dict would answer a typo with a KeyError at the call site; a dataclass answers at the
    boundary, and the field names are the specification's own."""
    payload = _payload()
    assert isinstance(payload.qc, ScanQC)
    assert payload.qc.deposited == 41
    assert payload.qc.dropped_ambiguous_path == 4
    assert payload.qc.sj_implicit_fragments == 7
    with pytest.raises(AttributeError):
        _ = payload.qc.dropped_ambigous_path  # noqa: B018 — the typo IS the test


def test_the_qc_fields_are_exactly_the_specifications_qc_keys():
    """Same argument as the Tally check: one vocabulary, and no list here to drift from it."""
    reference_keys = set(Tally.zeros(1, 0, 0, 1).qc)
    assert {f.name for f in dataclasses.fields(ScanQC)} == reference_keys


def test_the_start_count_invariant_is_checkable_from_the_payload_alone():
    """``Σ node_start_count == deposited`` — the accumulator's one non-tautological invariant.

    It replaced three "conservation identities" whose right-hand sides could only be evaluated by
    re-running the deposit, so a deliberately broken replay satisfied all three while 91 % of the
    crossings were junk. This one is checkable against a number the deposit counted independently, and
    the payload must carry both halves of it.
    """
    counts = np.zeros(5, np.uint32)
    counts[:3] = [10, 20, 11]
    payload = _payload(node_start_count=counts, qc=dict(_calibration_dict()["calibration"]["qc"]))
    assert int(payload.node_start_count.sum()) == payload.qc.deposited


# ---------------------------------------------------------------------------
# validation: every failure names both the observed and the expected value
# ---------------------------------------------------------------------------


def test_a_missing_calibration_block_is_an_ERROR_not_an_empty_payload():
    with pytest.raises(ValueError, match="set_regions"):
        AccumulatorPayload.from_scan_result({"calibration": None})


def test_a_wrong_column_count_is_rejected():
    with pytest.raises(ValueError, match="n_strand_columns"):
        _payload(n_strand_columns=4)


def test_an_offset_array_of_the_wrong_length_is_rejected():
    with pytest.raises(ValueError, match="ref_node_offsets"):
        _payload(ref_node_offsets=np.zeros(2, np.int64))


def test_an_array_that_does_not_divide_by_its_axis_is_rejected():
    """⛔ The failure mode this catches is a payload silently reshaped to the wrong number of rows."""
    with pytest.raises(ValueError, match="node_contained_count"):
        _payload(node_contained_count=np.arange(7, dtype=np.uint32))


def test_the_offsets_must_agree_with_the_cut_axis_they_describe():
    """A reference contributing ``c`` cuts owns ``c-1`` nodes and ``c-2`` lines. Re-derived here rather
    than trusted, because an offset array that merely has the right LENGTH can still be inconsistent."""
    bad = np.asarray(
        [0, 3, 5, 5], np.int64
    )  # chr1 claims 3 nodes; 4 cuts means 3 — chr2 claims 2, ok
    bad[1] = 99
    with pytest.raises(ValueError, match="ref_node_offsets"):
        _payload(ref_node_offsets=bad)


# ---------------------------------------------------------------------------
# ownership
# ---------------------------------------------------------------------------


def test_the_payload_holds_VIEWS_and_does_not_copy():
    """⚠ A live footgun, documented because someone will be tempted to 'add a cast for safety'.

    ``np.ascontiguousarray(x, dtype=D)`` is a **no-op** when the array already has dtype ``D``, so the
    payload holds views over the capsule-owned C++ heap and the payload object is the keep-alive. Adding
    a defensive cast would silently double peak memory on a 1.04 M-node partition.
    """
    cal = _calibration_dict()
    payload = AccumulatorPayload.from_scan_result(cal)
    source = cal["calibration"]["node_start_count"]
    assert payload.node_start_count.base is source or payload.node_start_count is source, (
        "the payload copied an array whose dtype already matched — that doubles peak memory"
    )


def test_a_deposited_lengths_HISTOGRAM_THAT_DOES_NOT_BIN_EVERY_FRAGMENT_IS_REJECTED():
    """⭐ **C1's invariant, refused at the door.** ``Σ deposited_lengths`` must equal ``qc.deposited``.

    This histogram is about to become the empirical-Bayes anchor for **every** fragment-length model in
    the tool (`docs/FRAGMENT_LENGTH_AUDIT.md`), so an off-by-N is not a cosmetic error — it silently
    re-weights the anchor against the pools it is supposed to anchor, which is a subtler version of the
    frame mismatch C1 exists to remove.

    ⚠ The check has to live at the payload boundary and not only in the accumulator's own tests, because
    the payload is what a **cached** scan is rebuilt from — and a cache can be truncated, partially
    written, or produced by a build whose schema digest happened to collide.
    """
    n = MAX_LENGTH + 1
    with pytest.raises(ValueError, match="deposited_lengths sums to"):
        _payload(deposited_lengths=_deposited_lengths(40))  # one short of qc.deposited = 41
    with pytest.raises(ValueError, match="deposited_lengths sums to"):
        _payload(deposited_lengths=_deposited_lengths(42))  # one too many
    with pytest.raises(ValueError, match="deposited_lengths has shape"):
        _payload(deposited_lengths=np.zeros(n - 1, dtype=np.uint32))
    # the self-consistent one is accepted
    assert int(_payload().deposited_lengths.sum()) == 41
