"""Schema invariants for :class:`rigel.scan_payload.AccumulatorPayload`.

    Spec: ``tests/native/_accumulator_reference.py``

The payload is the boundary where the C++ tally becomes Python. Its field names **are** the
specification's ``Tally`` field names, character for character, so the tests below check the schema
against `Tally` itself rather than against a hand-written list — a list would be free to drift, and the
whole point of sharing one vocabulary is that it cannot.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rigel.scan_payload import N_STRAND_COLUMNS, AccumulatorPayload, GapCensus, ScanQC
from rigel.types import Strand

from native._accumulator_reference import (
    Accumulator,
    GapHypothesis,
    Partition,
    Tally,
)


#: Two references: chr1 with 4 region_bounds (3 regions, 2 boundaries) and chr2 with 3 (2 regions, 1 boundary). A third
#: reference has NO region_bounds, which is legal and contributes nothing — the case where per-reference offset
#: arithmetic goes wrong if it is written as a plain subtraction.
REGION_BOUNDS_PER_REF = [[0, 100, 200, 600], [0, 500, 900], []]
MAX_LENGTH = 12

#: ⭐ The deferred bank in the fixture is produced by the SPECIFICATION rather than hand-written, so the
#: nested CSR the payload validates is one the reference actually emits. A hand-written bank is a second
#: encoding of the same layout, and the two would be free to drift — which is the whole reason the payload's
#: field names are read off ``Tally`` instead of listed.
#:
#: ⚠ Four fragments, because ``qc.deferred_undetermined_gap`` below is 4 and the payload refuses a bank
#: whose record count disagrees with the counter that describes it.
_DEFERRED_REGION_BOUNDS = [0, 100, 200, 300, 400, 500, 600]


def _deferred_bank(n_fragments: int = 4) -> dict[str, np.ndarray]:
    """``n_fragments`` deferred records, flattened exactly as ``Tally.deferred_arrays()`` specifies."""
    partition = Partition.from_region_bounds(
        [_DEFERRED_REGION_BOUNDS], region_types=[[0, 2, 1, 2, 1, 0]]
    )
    acc = Accumulator(partition, max_fragment_length=1000)
    for i in range(n_fragments):
        outcome = acc.deposit(
            0,
            100 + i,
            500,
            hypotheses=(
                GapHypothesis(((200, 300),), sj_strand=Strand.POS, supporting_t_inds=(i,)),
                GapHypothesis(),
            ),
        )
        assert outcome.value == "deferred_undetermined_gap", outcome
    return acc.tally.deferred_arrays()


def _gap_census(deferred: int = 4) -> dict[str, int]:
    """A census whose three ``gap_deferred_*`` sum to ``deferred`` — the partition the payload checks."""
    return {
        "gap_resolved_spliced": 7,
        "gap_deferred_rna_or_gdna": deferred,
        "gap_deferred_which_introns": 0,
        "gap_deferred_both": 0,
    }


def _deposited_lengths(total: int) -> np.ndarray:
    """A uint32[MAX_LENGTH + 1] histogram summing to exactly ``total``."""
    out = np.zeros(MAX_LENGTH + 1, dtype=np.uint32)
    out[MAX_LENGTH // 2] = total
    return out


def _calibration_dict(**overrides) -> dict:
    """A flat calibration dict shaped exactly as ``BamScanner::build_result`` emits one."""
    ref_region_bound_offsets = np.zeros(len(REGION_BOUNDS_PER_REF) + 1, np.int64)
    ref_region_offsets = np.zeros(len(REGION_BOUNDS_PER_REF) + 1, np.int64)
    ref_boundary_offsets = np.zeros(len(REGION_BOUNDS_PER_REF) + 1, np.int64)
    for f, region_bounds in enumerate(REGION_BOUNDS_PER_REF):
        ref_region_bound_offsets[f + 1] = ref_region_bound_offsets[f] + len(region_bounds)
        ref_region_offsets[f + 1] = ref_region_offsets[f] + max(len(region_bounds) - 1, 0)
        ref_boundary_offsets[f + 1] = ref_boundary_offsets[f] + max(len(region_bounds) - 2, 0)

    n_regions, n_boundaries, n_sj = (
        int(ref_region_offsets[-1]),
        int(ref_boundary_offsets[-1]),
        3,
    )
    ref_sj_offsets = np.asarray([0, 2, 3, 3], np.int64)

    cal = {
        "region_bounds": np.asarray(
            [c for region_bounds in REGION_BOUNDS_PER_REF for c in region_bounds], np.int64
        ),
        "ref_region_bound_offsets": ref_region_bound_offsets,
        "ref_region_offsets": ref_region_offsets,
        "ref_boundary_offsets": ref_boundary_offsets,
        "ref_sj_offsets": ref_sj_offsets,
        "region_contained_count": np.arange(n_regions * 2, dtype=np.uint32),
        "region_contained_inv_opportunity_sum": np.arange(n_regions, dtype=np.float64) * 7,
        "region_start_count": np.arange(n_regions * 2, dtype=np.uint32),
        "region_end_count": np.arange(n_regions * 2, dtype=np.uint32) * 3,
        "region_span_count": np.arange(n_regions * 2, dtype=np.uint32) * 5,
        "boundary_unspliced_count": np.arange(n_boundaries * 2, dtype=np.uint32),
        "boundary_unspliced_inv_length_sum": np.arange(n_boundaries, dtype=np.float64),
        "boundary_spliced_count": np.arange(n_boundaries * 2, dtype=np.uint32),
        # ⭐ The conserved masses: ONE value per boundary, so `n_boundaries` and not `n_boundaries * 2`. A fixture
        # that gave them a strand axis would pass every value check and disagree with the emitter.
        "boundary_unspliced_mass": np.arange(n_boundaries, dtype=np.float64) * 29,
        "boundary_spliced_mass": np.arange(n_boundaries, dtype=np.float64) * 31,
        "sj_count": np.arange(n_sj * 2, dtype=np.uint32),
        "sj_inv_length_sum": np.arange(n_sj, dtype=np.float64),
        "sj_mass": np.arange(n_sj * 2, dtype=np.float64),
        "pool_lengths": np.arange(5 * (MAX_LENGTH + 1), dtype=np.int64),
        # ⭐ TRAPS: a-purity-filter-is-a-length-filter: the unconditional histogram must bin EXACTLY the deposited fragments, so this fixture
        # can no longer carry an arbitrary array — 41 here, matching qc.deposited below. That coupling is
        # the invariant doing its job at the door.
        "deposited_lengths": _deposited_lengths(41),
        "qc": {
            "deposited": 41,
            "dropped_too_long": 1,
            "dropped_empty": 2,
            "dropped_strand_undefined": 3,
            "deferred_undetermined_gap": 4,
            "unannotated_introns": 5,
            "contradictory_sj_strand": 6,
            "introns_absorbed": 8,
        },
        "deferred": _deferred_bank(),
        "gap_resolution": _gap_census(),
        "n_strand_columns": N_STRAND_COLUMNS,
        "n_fragment_pools": 5,
        "max_length": MAX_LENGTH,
        "n_refs": len(REGION_BOUNDS_PER_REF),
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
    n_regions, n_boundaries, n_sj = payload.n_regions, payload.n_boundaries, payload.n_sj
    assert (n_regions, n_boundaries, n_sj) == (5, 3, 3)
    # ⭐ The COUNTS keep both genome-strand columns — the strand model is a Beta-Binomial over them.
    for name, rows in (
        ("region_contained_count", n_regions),
        ("boundary_unspliced_count", n_boundaries),
        ("boundary_spliced_count", n_boundaries),
        ("sj_count", n_sj),
        ("sj_mass", n_sj),
    ):
        assert getattr(payload, name).shape == (rows, N_STRAND_COLUMNS), name
    # ⛔ The length moments and the conserved masses carry ONE column: which strand a read aligned to
    # says nothing about whether the molecule was gDNA or RNA, and every consumer summed the two.
    for name, rows in (
        ("region_contained_inv_opportunity_sum", n_regions),
        ("boundary_unspliced_inv_length_sum", n_boundaries),
        ("boundary_unspliced_mass", n_boundaries),
        ("boundary_spliced_mass", n_boundaries),
        ("sj_inv_length_sum", n_sj),
    ):
        assert getattr(payload, name).shape == (rows,), name
    assert payload.region_start_count.shape == (n_regions, 2)
    assert payload.region_end_count.shape == (n_regions, 2)
    assert payload.region_span_count.shape == (n_regions, 2)
    assert payload.pool_lengths.shape == (5, MAX_LENGTH + 1)


def test_the_dtypes_are_the_specifications_dtypes():
    """⚠ Counts are uint32 and densities uint64, and the payload must not silently widen either.

    A count that arrives as int64 compares equal to the specification's uint32 by value, so a value-only
    check would pass while the schema had changed underneath it.

    ⚠ Three ``Tally`` fields are not arrays — ``qc`` and ``gap_resolution`` are dicts of counters and
    ``deferred`` is a list of records — and each is checked by its own test below. They are skipped by
    ASKING THE REFERENCE what type it holds, never by naming them here: a name would let a field that
    stopped being an array drop silently out of this gate.
    """
    payload = _payload()
    reference = Tally.zeros(n_regions=5, n_boundaries=3, n_sj=3, max_length=MAX_LENGTH)
    checked = 0
    for field in dataclasses.fields(Tally):
        expected = getattr(reference, field.name)
        if not isinstance(expected, np.ndarray):
            continue
        assert getattr(payload, field.name).dtype == expected.dtype, field.name
        checked += 1
    # ⚠ The floor moves only when the SCHEMA moves, and then deliberately. 18 arrays → 20 when the two
    # conserved masses landed → 14 when the six dead banks were removed (three ``region_spanning_*``, the
    # two spliced-boundary length moments, ``sj_length_sum``) → **12** on 2026-08-13 with
    # ``region_contained_length_sum`` and ``boundary_unspliced_length_sum``, whose stated justification did not
    # survive measurement (`scan_payload`'s docstring has the retraction). A floor that drifted down on
    # its own would be this gate quietly narrowing, which is the one thing it exists to catch.
    # ⚠ ``sj_mass`` going per-strand in that same change did NOT move this floor: it is one array either
    # way, and its SHAPE is gated by the test above rather than here.
    assert checked >= 12, f"only {checked} arrays compared; the gate has narrowed"


def test_the_DEFERRED_bank_is_int64_throughout_and_carries_the_specifications_arrays():
    """⭐ The side buffer's own schema check, since it is not one array and cannot join the loop above.

    ⚠ One dtype for the whole bank, and it is ``int64`` even for the two strand columns — which are
    ``int32`` everywhere else in the scanner. The parity gate compares dtypes, so a narrowing conversion at
    the ABI would compare equal by value and hide the change.
    """
    payload = _payload()
    reference = Tally.zeros(n_regions=5, n_boundaries=3, n_sj=3, max_length=MAX_LENGTH)
    expected = set(reference.deferred_arrays())
    actual = {f.name for f in dataclasses.fields(payload.deferred)}
    assert actual == expected, (
        f"the payload's deferred bank carries {sorted(actual)}; the specification's flattening emits "
        f"{sorted(expected)}"
    )
    for name in sorted(expected):
        array = getattr(payload.deferred, name)
        assert array.dtype == np.int64, f"deferred.{name} has dtype {array.dtype}, expected int64"
    assert payload.deferred.n_fragments == payload.qc.deferred_undetermined_gap == 4
    assert payload.deferred.n_hypotheses == 8, "two hypotheses per deferred fragment"


def test_a_WRONG_dtype_is_REJECTED_rather_than_coerced():
    """⛔ Checking the output dtype is not enough, and a perturbation proved it.

    ``ascontiguousarray(x, dtype=uint32)`` will happily narrow an int64 array, so a payload that *coerces*
    still reports the right dtype and passes a check on its own output. What that hides is a C++ side that
    has stopped agreeing with the schema — the dtype is part of byte-identity, and a widened count compares
    equal by value all the way down. So the wrong dtype has to arrive and be refused.
    """
    with pytest.raises(ValueError, match="region_contained_count.*dtype"):
        _payload(region_contained_count=np.arange(10, dtype=np.int64))
    with pytest.raises(ValueError, match="sj_inv_length_sum.*dtype"):
        _payload(sj_inv_length_sum=np.arange(6, dtype=np.uint32))


def test_a_MISSING_qc_denominator_is_REJECTED():
    """⛔ Also found by perturbation: nothing was feeding an incomplete qc block.

    Design §10.3 requires every one of these to be emitted, because every conservation statement
    downstream has to be able to name what it excluded. A denominator that silently arrives absent is a
    statement that cannot.
    """
    qc = dict(_calibration_dict()["calibration"]["qc"])
    del qc["deferred_undetermined_gap"]
    with pytest.raises(ValueError, match="deferred_undetermined_gap"):
        _payload(qc=qc)


def test_a_reference_with_no_region_bounds_contributes_nothing_to_any_axis():
    """The third reference is empty. Per-reference offsets written as a plain subtraction go negative."""
    payload = _payload()
    assert payload.ref_region_offsets[-1] == payload.ref_region_offsets[-2]
    assert payload.ref_boundary_offsets[-1] == payload.ref_boundary_offsets[-2]
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
    assert payload.qc.deferred_undetermined_gap == 4
    with pytest.raises(AttributeError):
        _ = payload.qc.dropped_ambigous_path  # noqa: B018 — the typo IS the test


def test_the_qc_fields_are_exactly_the_specifications_qc_keys():
    """Same argument as the Tally check: one vocabulary, and no list here to drift from it."""
    reference_keys = set(Tally.zeros(1, 0, 0, 1).qc)
    assert {f.name for f in dataclasses.fields(ScanQC)} == reference_keys


def test_the_gap_census_fields_are_exactly_the_specifications_keys():
    """⛔ Including the ABSENCE of ``gap_resolved_unspliced``.

    That class existed and no fragment could enter it: a spliced hypothesis region_bounds bases the unspliced one
    keeps, so the unspliced path is always the longest and can never be the sole survivor. Reading the key
    set off the specification is what stops it reappearing on one side only.
    """
    reference_keys = set(Tally.zeros(1, 0, 0, 1).gap_resolution)
    assert {f.name for f in dataclasses.fields(GapCensus)} == reference_keys
    assert "gap_resolved_unspliced" not in reference_keys


def test_a_MISSING_gap_census_subclass_is_REJECTED():
    """The subclasses are exhaustive by construction, so a missing one is a partition that does not close."""
    census = _gap_census()
    del census["gap_deferred_both"]
    with pytest.raises(ValueError, match="gap_deferred_both"):
        _payload(gap_resolution=census)


def test_a_DEFERRED_BANK_THAT_DISAGREES_WITH_ITS_OWN_COUNTER_IS_REJECTED():
    """⭐ **The conservation half, refused at the door.**

    ``deposited + deferred + dropped_* == offered`` is worth nothing if the deferred term is a number with
    no fragments behind it. ⚠ The check has to live at the payload boundary and not only in the
    accumulator's tests, because the payload is what a **cached** scan is rebuilt from — and a cache can be
    truncated, partially written, or produced by a build whose schema digest happened to collide.
    """
    with pytest.raises(ValueError, match="deferred bank holds 3 fragments"):
        _payload(deferred=_deferred_bank(3))
    with pytest.raises(ValueError, match="deferred bank holds 5 fragments"):
        _payload(deferred=_deferred_bank(5))
    with pytest.raises(ValueError, match="partition that does not close"):
        _payload(gap_resolution=_gap_census(deferred=3))
    assert _payload().deferred.n_fragments == 4  # the self-consistent one is accepted


def test_a_TRUNCATED_deferred_CSR_is_REJECTED_rather_than_indexed_off_the_end():
    """⛔ The second pass indexes every one of these arrays.

    A bank whose offsets outrun its values does not fail loudly when it is read — it scores one fragment
    against another fragment's hypotheses, or reads zeros, and returns a plausible answer. So the CSR is
    re-derived at the door, exactly as ``ref_region_offsets`` is.
    """
    bank = _deferred_bank()
    truncated = dict(bank) | {"hypothesis_introns": bank["hypothesis_introns"][:-2]}
    with pytest.raises(ValueError, match="hypothesis_intron_offsets.*ends at"):
        _payload(deferred=truncated)

    short_record = dict(bank) | {"start": bank["start"][:-1]}
    with pytest.raises(ValueError, match="deferred\\['start'\\] has 3 entries"):
        _payload(deferred=short_record)

    backwards = dict(bank) | {"hypothesis_offsets": bank["hypothesis_offsets"][::-1].copy()}
    with pytest.raises(ValueError, match="hypothesis_offsets"):
        _payload(deferred=backwards)

    widened = dict(bank) | {"ref": bank["ref"].astype(np.int32)}
    with pytest.raises(ValueError, match="deferred\\['ref'\\] has dtype int32"):
        _payload(deferred=widened)


def test_a_DEFERRED_RECORD_WITH_FEWER_THAN_TWO_HYPOTHESES_IS_REJECTED():
    """⭐ A fragment is deferred BECAUSE two or more hypotheses survived.

    A record carrying one is a bank that lost the others, and the second pass would then "choose" from a
    set of one and deposit an answer nothing supported — the exact outcome the deferral exists to prevent.
    """
    bank = _deferred_bank()
    # Drop the LAST record's second (unspliced) hypothesis, which carries neither introns nor supporting
    # transcripts — so every other array stays self-consistent and the ONLY thing wrong is the run length.
    # ⚠ Written this way on purpose: a bank that also broke the CSR would be caught by the test above and
    # this one would pass for the wrong reason.
    lone = dict(bank) | {
        "hypothesis_offsets": np.asarray([0, 2, 4, 6, 7], np.int64),
        "hypothesis_sj_strand": bank["hypothesis_sj_strand"][:-1],
        "hypothesis_intron_offsets": bank["hypothesis_intron_offsets"][:-1],
        "hypothesis_t_offsets": bank["hypothesis_t_offsets"][:-1],
    }
    with pytest.raises(ValueError, match="carries 1 hypotheses"):
        _payload(deferred=lone)


def test_the_start_count_invariant_is_checkable_from_the_payload_alone():
    """``Σ region_start_count == deposited`` — the accumulator's one non-tautological invariant.

    It replaced three "conservation identities" whose right-hand sides could only be evaluated by
    re-running the deposit, so a deliberately broken replay satisfied all three while 91 % of the
    crossings were junk. This one is checkable against a number the deposit counted independently, and
    the payload must carry both halves of it.
    """
    counts = np.zeros(10, np.uint32)  # 5 regions x 2 strand columns, flat as the C++ emits it
    counts[:3] = [10, 20, 11]
    payload = _payload(region_start_count=counts, qc=dict(_calibration_dict()["calibration"]["qc"]))
    assert int(payload.region_start_count.sum()) == payload.qc.deposited
    # ⭐ the ledger closes TWICE over since 2026-08-21: the END bank is its mirror
    ends = np.zeros(10, np.uint32)
    ends[2:6] = [5, 6, 7, 23]
    payload = _payload(region_end_count=ends, qc=dict(_calibration_dict()["calibration"]["qc"]))
    assert int(payload.region_end_count.sum()) == payload.qc.deposited


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
    with pytest.raises(ValueError, match="ref_region_offsets"):
        _payload(ref_region_offsets=np.zeros(2, np.int64))


def test_an_array_that_does_not_divide_by_its_axis_is_rejected():
    """⛔ The failure mode this catches is a payload silently reshaped to the wrong number of rows."""
    with pytest.raises(ValueError, match="region_contained_count"):
        _payload(region_contained_count=np.arange(7, dtype=np.uint32))


def test_the_offsets_must_agree_with_the_region_bound_axis_they_describe():
    """A reference contributing ``c`` region_bounds owns ``c-1`` regions and ``c-2`` boundaries. Re-derived here rather
    than trusted, because an offset array that merely has the right LENGTH can still be inconsistent."""
    bad = np.asarray(
        [0, 3, 5, 5], np.int64
    )  # chr1 claims 3 regions; 4 region_bounds means 3 — chr2 claims 2, ok
    bad[1] = 99
    with pytest.raises(ValueError, match="ref_region_offsets"):
        _payload(ref_region_offsets=bad)


# ---------------------------------------------------------------------------
# ownership
# ---------------------------------------------------------------------------


def test_the_payload_holds_VIEWS_and_does_not_copy():
    """⚠ A live footgun, documented because someone will be tempted to 'add a cast for safety'.

    ``np.ascontiguousarray(x, dtype=D)`` is a **no-op** when the array already has dtype ``D``, so the
    payload holds views over the capsule-owned C++ heap and the payload object is the keep-alive. Adding
    a defensive cast would silently double peak memory on a 1.04 M-region partition.
    """
    cal = _calibration_dict()
    payload = AccumulatorPayload.from_scan_result(cal)
    source = cal["calibration"]["region_start_count"]
    assert payload.region_start_count.base is source or payload.region_start_count is source, (
        "the payload copied an array whose dtype already matched — that doubles peak memory"
    )


def test_a_deposited_lengths_HISTOGRAM_THAT_DOES_NOT_BIN_EVERY_FRAGMENT_IS_REJECTED():
    """⭐ **TRAPS: a-purity-filter-is-a-length-filter's invariant, refused at the door.** ``Σ deposited_lengths`` must equal ``qc.deposited``.

    This histogram is about to become the empirical-Bayes anchor for **every** fragment-length model in
    the tool, so an off-by-N is not a cosmetic error — it silently
    re-weights the anchor against the pools it is supposed to anchor, which is a subtler version of the
    frame mismatch TRAPS: a-purity-filter-is-a-length-filter exists to remove.

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
