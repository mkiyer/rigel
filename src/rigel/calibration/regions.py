"""rigel.calibration.regions — Fine-grained calibration region schema.

This module defines:

- ``RegionType`` — the 3-way partition (INTERGENIC / INTRON / EXON).
- ``RegionStrand`` — bitwise OR-able strand flags.
- ``BoundaryKind`` — a small enum classifying adjacent fine-region edges.
- ``build_fine_region_table`` — vectorized-style event sweep that emits the
    v4 ``regions.feather`` table keyed by a 4-bit fine signature.
- ``load_regions(path)`` — load and dtype-coerce ``regions.feather``.
- ``validate_against_ref_lengths(region_df, ref_lengths)`` — enforce the
    load-time invariants from the fine-region migration plan.

The genome is tiled per-reference into maximal intervals with identical
``{intron_pos, intron_neg, exon_pos, exon_neg}`` annotation state. Coarse
``type`` and ``strand`` fields remain as derived bridge columns until the
fractional accumulator cutover removes the old consumers.
"""

from __future__ import annotations

from enum import IntEnum, IntFlag
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

from ..transcript import Transcript
from ..types import Strand


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


class RegionType(IntEnum):
    """The three region types in the genome partition."""

    INTERGENIC = 0
    INTRON = 1
    EXON = 2


class RegionStrand(IntFlag):
    """Strand flags for a region.

    ``AMBIG`` is the bitwise OR ``POS | NEG`` and indicates that the region
    contains transcripts (or exons, for ``EXON`` regions) on both strands.
    """

    NONE = 0
    POS = 1
    NEG = 2
    AMBIG = 3  # POS | NEG


class BoundaryKind(IntEnum):
    """Adjacent-region boundary classification for fine-region indexes."""

    NONE = 0
    EXON_INTRON = 1
    EXON_INTERGENIC = 2
    INTRON_INTERGENIC = 3
    EXON_EXON_STRAND_SWITCH = 4
    OTHER = 5


# ---------------------------------------------------------------------------
# v4 row schema
# ---------------------------------------------------------------------------

SIGNATURE_SENTINEL = np.uint8(0xFF)

REGION_COLUMNS = [
    "region_id",
    "ref_name",
    "start",
    "end",
    "length",
    "signature",
    "intron_pos",
    "intron_neg",
    "exon_pos",
    "exon_neg",
    "type",
    "strand",
    "boundary_flux_left",
    "boundary_flux_right",
    "left_signature",
    "right_signature",
    "boundary_kind_left",
    "boundary_kind_right",
]

REGION_COLUMN_DTYPES: dict[str, type | np.dtype] = {
    "region_id": np.int64,
    "start": np.int64,
    "end": np.int64,
    "length": np.int64,
    "signature": np.uint8,
    "intron_pos": np.bool_,
    "intron_neg": np.bool_,
    "exon_pos": np.bool_,
    "exon_neg": np.bool_,
    "type": np.uint8,
    "strand": np.uint8,
    "boundary_flux_left": np.bool_,
    "boundary_flux_right": np.bool_,
    "left_signature": np.uint8,
    "right_signature": np.uint8,
    "boundary_kind_left": np.uint8,
    "boundary_kind_right": np.uint8,
}


# ---------------------------------------------------------------------------
# Fine-region builder
# ---------------------------------------------------------------------------


def classify_boundary_kind(left_signature: int, right_signature: int) -> int:
    """Classify the boundary between two adjacent fine signatures.

    The classification is symmetric: ``classify(a, b) == classify(b, a)``.
    ``SIGNATURE_SENTINEL`` represents a reference end and yields
    :class:`BoundaryKind.NONE`.
    """
    from . import signature as sig

    left = int(left_signature)
    right = int(right_signature)
    if left == int(SIGNATURE_SENTINEL) or right == int(SIGNATURE_SENTINEL):
        return int(BoundaryKind.NONE)

    left_exon_pos = bool(left & sig.BIT_EXON_POS)
    left_exon_neg = bool(left & sig.BIT_EXON_NEG)
    left_intron_pos = bool(left & sig.BIT_INTRON_POS)
    left_intron_neg = bool(left & sig.BIT_INTRON_NEG)

    right_exon_pos = bool(right & sig.BIT_EXON_POS)
    right_exon_neg = bool(right & sig.BIT_EXON_NEG)
    right_intron_pos = bool(right & sig.BIT_INTRON_POS)
    right_intron_neg = bool(right & sig.BIT_INTRON_NEG)

    same_strand_exon_intron = (
        (left_exon_pos and right_intron_pos)
        or (left_intron_pos and right_exon_pos)
        or (left_exon_neg and right_intron_neg)
        or (left_intron_neg and right_exon_neg)
    )
    if same_strand_exon_intron:
        return int(BoundaryKind.EXON_INTRON)

    left_has_exon = left_exon_pos or left_exon_neg
    right_has_exon = right_exon_pos or right_exon_neg
    left_has_intron = left_intron_pos or left_intron_neg
    right_has_intron = right_intron_pos or right_intron_neg
    left_intergenic = left == 0
    right_intergenic = right == 0

    if (left_intergenic and right_has_exon) or (right_intergenic and left_has_exon):
        return int(BoundaryKind.EXON_INTERGENIC)
    if (left_intergenic and right_has_intron) or (right_intergenic and left_has_intron):
        return int(BoundaryKind.INTRON_INTERGENIC)

    same_strand_exon_exon = (left_exon_pos and right_exon_pos) or (left_exon_neg and right_exon_neg)
    opposite_strand_exon_exon = (left_exon_pos and right_exon_neg) or (
        left_exon_neg and right_exon_pos
    )
    if left_has_exon and right_has_exon and opposite_strand_exon_exon and not same_strand_exon_exon:
        return int(BoundaryKind.EXON_EXON_STRAND_SWITCH)

    return int(BoundaryKind.OTHER)


def _add_interval_events(
    events: list[tuple[int, int, int]],
    *,
    ref_name: str,
    ref_length: int,
    transcript_id: str | None,
    start: int,
    end: int,
    bit: int,
) -> None:
    """Append start/end delta events for one validated half-open interval."""
    start = int(start)
    end = int(end)
    if end <= start:
        return
    if start < 0 or end > ref_length:
        label = transcript_id if transcript_id is not None else "<unknown>"
        raise ValueError(
            f"Transcript {label} has interval [{start}, {end}) outside "
            f"reference {ref_name!r} length {ref_length}."
        )
    events.append((start, bit, +1))
    events.append((end, bit, -1))


def _pack_counts(
    intron_pos_count: int,
    intron_neg_count: int,
    exon_pos_count: int,
    exon_neg_count: int,
) -> int:
    from .signature import pack_signature

    return pack_signature(
        intron_pos=intron_pos_count > 0,
        intron_neg=intron_neg_count > 0,
        exon_pos=exon_pos_count > 0,
        exon_neg=exon_neg_count > 0,
    )


def _region_row(ref_name: str, start: int, end: int, signature: int) -> dict[str, object]:
    from . import signature as sig

    signature = sig.validate_signature(signature)
    return {
        "region_id": -1,
        "ref_name": ref_name,
        "start": int(start),
        "end": int(end),
        "length": int(end) - int(start),
        "signature": signature,
        "intron_pos": bool(signature & sig.BIT_INTRON_POS),
        "intron_neg": bool(signature & sig.BIT_INTRON_NEG),
        "exon_pos": bool(signature & sig.BIT_EXON_POS),
        "exon_neg": bool(signature & sig.BIT_EXON_NEG),
        "type": sig.coarse_type_from_signature(signature),
        "strand": sig.coarse_strand_from_signature(signature),
        "boundary_flux_left": False,
        "boundary_flux_right": False,
        "left_signature": int(SIGNATURE_SENTINEL),
        "right_signature": int(SIGNATURE_SENTINEL),
        "boundary_kind_left": int(BoundaryKind.NONE),
        "boundary_kind_right": int(BoundaryKind.NONE),
    }


def _build_reference_rows(
    ref_name: str,
    ref_length: int,
    events: list[tuple[int, int, int]],
) -> list[dict[str, object]]:
    from . import signature as sig

    if ref_length < 0:
        raise ValueError(f"Reference {ref_name!r} has negative length {ref_length}.")
    if ref_length == 0:
        return []

    events.sort(key=lambda event: event[0])
    rows: list[dict[str, object]] = []
    intron_pos_count = 0
    intron_neg_count = 0
    exon_pos_count = 0
    exon_neg_count = 0
    previous_position = 0
    event_index = 0
    n_events = len(events)

    def emit_segment(segment_start: int, segment_end: int) -> None:
        if segment_end <= segment_start:
            return
        fine_signature = _pack_counts(
            intron_pos_count,
            intron_neg_count,
            exon_pos_count,
            exon_neg_count,
        )
        if rows and rows[-1]["end"] == segment_start and rows[-1]["signature"] == fine_signature:
            rows[-1]["end"] = segment_end
            rows[-1]["length"] = segment_end - int(rows[-1]["start"])
            return
        rows.append(_region_row(ref_name, segment_start, segment_end, fine_signature))

    while event_index < n_events:
        position = events[event_index][0]
        if position > previous_position:
            emit_segment(previous_position, position)
            previous_position = position

        while event_index < n_events and events[event_index][0] == position:
            _, bit, delta = events[event_index]
            if bit == sig.BIT_INTRON_POS:
                intron_pos_count += delta
            elif bit == sig.BIT_INTRON_NEG:
                intron_neg_count += delta
            elif bit == sig.BIT_EXON_POS:
                exon_pos_count += delta
            elif bit == sig.BIT_EXON_NEG:
                exon_neg_count += delta
            else:  # pragma: no cover - internal guard
                raise RuntimeError(f"unknown fine-region bit {bit}")
            event_index += 1

    if ref_length > previous_position:
        emit_segment(previous_position, ref_length)
    return rows


def _finalize_region_rows(rows: list[dict[str, object]]) -> pd.DataFrame:
    if not rows:
        return _coerce_region_dtypes(pd.DataFrame(columns=REGION_COLUMNS))

    for region_id, row in enumerate(rows):
        row["region_id"] = region_id

    frame = pd.DataFrame(rows, columns=REGION_COLUMNS)

    for _, ref_index in frame.groupby("ref_name", sort=False).groups.items():
        positions = np.asarray(ref_index, dtype=np.int64)
        signatures = frame.loc[positions, "signature"].to_numpy(np.uint8, copy=False)
        types = frame.loc[positions, "type"].to_numpy(np.uint8, copy=False)

        left_signatures = np.full(len(positions), int(SIGNATURE_SENTINEL), dtype=np.uint8)
        right_signatures = np.full(len(positions), int(SIGNATURE_SENTINEL), dtype=np.uint8)
        if len(positions) > 1:
            left_signatures[1:] = signatures[:-1]
            right_signatures[:-1] = signatures[1:]

        boundary_kind_left = np.full(len(positions), int(BoundaryKind.NONE), dtype=np.uint8)
        boundary_kind_right = np.full(len(positions), int(BoundaryKind.NONE), dtype=np.uint8)
        for offset in range(len(positions) - 1):
            kind = classify_boundary_kind(int(signatures[offset]), int(signatures[offset + 1]))
            boundary_kind_right[offset] = kind
            boundary_kind_left[offset + 1] = kind

        frame.loc[positions, "left_signature"] = left_signatures
        frame.loc[positions, "right_signature"] = right_signatures
        frame.loc[positions, "boundary_kind_left"] = boundary_kind_left
        frame.loc[positions, "boundary_kind_right"] = boundary_kind_right
        frame.loc[positions, "boundary_flux_left"] = types == int(RegionType.EXON)
        frame.loc[positions, "boundary_flux_right"] = types == int(RegionType.EXON)
        if len(positions) > 0:
            frame.loc[positions[0], "boundary_flux_left"] = False
            frame.loc[positions[-1], "boundary_flux_right"] = False

    return _coerce_region_dtypes(frame)


def build_fine_region_table(
    transcripts: list[Transcript],
    ref_lengths: Mapping[str, int],
) -> pd.DataFrame:
    """Build the v4 fine-grained calibration region table.

    The builder sweeps exon and transcript-intron interval events per
    reference. Synthetic transcripts are ignored, and every reference in
    ``ref_lengths`` is represented by a tiling row set unless its length is
    zero.
    """
    from . import signature as sig

    normalized_ref_lengths = {
        str(ref_name): int(ref_length) for ref_name, ref_length in ref_lengths.items()
    }
    events_by_ref: dict[str, list[tuple[int, int, int]]] = {
        ref_name: [] for ref_name in normalized_ref_lengths
    }

    for transcript in transcripts:
        if transcript.is_synthetic or not transcript.exons:
            continue
        if str(transcript.ref) not in normalized_ref_lengths:
            raise ValueError(
                f"Transcript {transcript.t_id} has reference {transcript.ref!r} "
                "not found in ref_lengths."
            )
        ref_name = str(transcript.ref)
        ref_length = normalized_ref_lengths[ref_name]
        events = events_by_ref[ref_name]
        if transcript.strand == Strand.POS:
            exon_bit = sig.BIT_EXON_POS
            intron_bit = sig.BIT_INTRON_POS
        elif transcript.strand == Strand.NEG:
            exon_bit = sig.BIT_EXON_NEG
            intron_bit = sig.BIT_INTRON_NEG
        else:
            continue

        for exon in transcript.exons:
            _add_interval_events(
                events,
                ref_name=ref_name,
                ref_length=ref_length,
                transcript_id=transcript.t_id,
                start=int(exon.start),
                end=int(exon.end),
                bit=exon_bit,
            )

        for intron_start, intron_end in transcript.introns():
            _add_interval_events(
                events,
                ref_name=ref_name,
                ref_length=ref_length,
                transcript_id=transcript.t_id,
                start=int(intron_start),
                end=int(intron_end),
                bit=intron_bit,
            )

    rows: list[dict[str, object]] = []
    for ref_name, ref_length in normalized_ref_lengths.items():
        rows.extend(_build_reference_rows(ref_name, ref_length, events_by_ref[ref_name]))

    return _finalize_region_rows(rows)


# ---------------------------------------------------------------------------
# Load + validation
# ---------------------------------------------------------------------------


def _missing_schema_error(path: str | Path, missing: set[str]) -> ValueError:
    if "signature" in missing:
        return ValueError(
            f"regions.feather at {path} is missing 'signature' (pre-v4 index). "
            "Rebuild the index (rigel index --fasta ... --gtf ...)."
        )
    return ValueError(
        f"regions.feather at {path} is missing required columns: "
        f"{sorted(missing)}. Rebuild the index "
        "(rigel index --fasta ... --gtf ...)."
    )


def _coerce_region_dtypes(region_df: pd.DataFrame) -> pd.DataFrame:
    region_df = region_df.copy()
    for column in REGION_COLUMNS:
        if column not in region_df.columns:
            region_df[column] = pd.Series(dtype="string" if column == "ref_name" else object)
    region_df = region_df.loc[:, REGION_COLUMNS]
    for column, dtype in REGION_COLUMN_DTYPES.items():
        if region_df[column].dtype != dtype:
            region_df[column] = region_df[column].astype(dtype)
    region_df["ref_name"] = region_df["ref_name"].astype("string")
    return region_df


def load_regions(path: str | Path) -> pd.DataFrame:
    """Read ``regions.feather`` and coerce dtypes to the v4 schema.

    The returned DataFrame has its row index reset to ``region_id``
    (which is itself ``0..N-1``).  Raises :class:`ValueError` if the
    file is missing required columns or is a pre-v4 region table.
    """
    df = pd.read_feather(str(path))
    required = set(REGION_COLUMNS)
    missing = required - set(df.columns)
    if missing:
        raise _missing_schema_error(path, missing)
    return _coerce_region_dtypes(df).reset_index(drop=True)


def _require_region_columns(region_df: pd.DataFrame) -> None:
    missing = set(REGION_COLUMNS) - set(region_df.columns)
    if missing:
        raise _missing_schema_error("<in-memory>", missing)


def _expected_signature_flags(signatures: np.ndarray) -> tuple[np.ndarray, ...]:
    from . import signature as sig

    return (
        (signatures & sig.BIT_INTRON_POS) != 0,
        (signatures & sig.BIT_INTRON_NEG) != 0,
        (signatures & sig.BIT_EXON_POS) != 0,
        (signatures & sig.BIT_EXON_NEG) != 0,
    )


def _validate_derived_columns(region_df: pd.DataFrame) -> None:
    from . import signature as sig

    signatures = region_df["signature"].to_numpy(np.uint8, copy=False)
    if ((signatures > 0x0F) & (signatures != int(SIGNATURE_SENTINEL))).any():
        bad_row = int(
            np.flatnonzero((signatures > 0x0F) & (signatures != int(SIGNATURE_SENTINEL)))[0]
        )
        raise ValueError(
            f"regions.feather: region {bad_row} has invalid signature "
            f"{int(signatures[bad_row])}. Rebuild the index."
        )
    if (signatures == int(SIGNATURE_SENTINEL)).any():
        bad_row = int(np.flatnonzero(signatures == int(SIGNATURE_SENTINEL))[0])
        raise ValueError(
            f"regions.feather: region {bad_row} uses reserved signature 0xFF. Rebuild the index."
        )

    flag_columns = ("intron_pos", "intron_neg", "exon_pos", "exon_neg")
    for column, expected in zip(flag_columns, _expected_signature_flags(signatures), strict=True):
        actual = region_df[column].to_numpy(bool, copy=False)
        if not np.array_equal(actual, expected):
            bad_row = int(np.flatnonzero(actual != expected)[0])
            raise ValueError(
                f"regions.feather: region {bad_row} column {column!r} does "
                "not match 'signature'. Rebuild the index."
            )

    expected_type = np.fromiter(
        (sig.coarse_type_from_signature(int(value)) for value in signatures),
        dtype=np.uint8,
        count=len(signatures),
    )
    actual_type = region_df["type"].to_numpy(np.uint8, copy=False)
    if not np.array_equal(actual_type, expected_type):
        bad_row = int(np.flatnonzero(actual_type != expected_type)[0])
        raise ValueError(
            f"regions.feather: region {bad_row} has type {int(actual_type[bad_row])}, "
            f"expected {int(expected_type[bad_row])} from signature. Rebuild the index."
        )

    expected_strand = np.fromiter(
        (sig.coarse_strand_from_signature(int(value)) for value in signatures),
        dtype=np.uint8,
        count=len(signatures),
    )
    actual_strand = region_df["strand"].to_numpy(np.uint8, copy=False)
    if not np.array_equal(actual_strand, expected_strand):
        bad_row = int(np.flatnonzero(actual_strand != expected_strand)[0])
        raise ValueError(
            f"regions.feather: region {bad_row} has strand "
            f"{int(actual_strand[bad_row])}, expected {int(expected_strand[bad_row])} "
            "from signature. Rebuild the index."
        )


def validate_against_ref_lengths(
    region_df: pd.DataFrame,
    ref_lengths: Mapping[str, int],
) -> None:
    """Enforce v4 fine-region table invariants.

    Raises :class:`ValueError` with an actionable message on the first
    violation. Checks (in order):

    1. ``region_id`` column equals row index.
    2. All region lengths are positive and match ``end - start``.
    3. Fine signature-derived columns are internally consistent.
    4. Every ``ref_name`` exists in ``ref_lengths``.
    5. Per-reference: regions are sorted by start, non-overlapping, and
       cover ``[0, ref_lengths[ref])`` with no gaps.
    """
    _require_region_columns(region_df)
    if not (region_df.index.to_numpy() == region_df["region_id"].to_numpy()).all():
        raise ValueError(
            "regions.feather: row index does not match 'region_id' column. Rebuild the index."
        )

    lengths = region_df["end"].to_numpy() - region_df["start"].to_numpy()
    if (lengths <= 0).any():
        bad = int(np.flatnonzero(lengths <= 0)[0])
        raise ValueError(
            f"regions.feather: region {bad} has non-positive length "
            f"(end - start = {int(lengths[bad])}). Rebuild the index."
        )
    stored_lengths = region_df["length"].to_numpy(np.int64, copy=False)
    if not np.array_equal(stored_lengths, lengths):
        bad = int(np.flatnonzero(stored_lengths != lengths)[0])
        raise ValueError(
            f"regions.feather: region {bad} has length {int(stored_lengths[bad])}, "
            f"expected {int(lengths[bad])}. Rebuild the index."
        )

    _validate_derived_columns(region_df)

    refs_in_table = set(region_df["ref_name"].unique().tolist())
    unknown = refs_in_table - set(ref_lengths.keys())
    if unknown:
        raise ValueError(
            f"regions.feather references {sorted(unknown)} not found in "
            f"ref_lengths.feather. Rebuild the index."
        )

    for ref, ref_len in ref_lengths.items():
        sub = region_df[region_df["ref_name"] == ref]
        if sub.empty:
            if ref_len > 0:
                raise ValueError(
                    f"regions.feather: reference '{ref}' (length {ref_len}) "
                    f"has no regions. Rebuild the index."
                )
            continue
        starts = sub["start"].to_numpy()
        ends = sub["end"].to_numpy()
        order = np.argsort(starts, kind="stable")
        if not np.array_equal(order, np.arange(len(starts))):
            raise ValueError(
                f"regions.feather: reference '{ref}' regions are not in "
                f"sorted start order. Rebuild the index."
            )
        if starts[0] != 0:
            raise ValueError(
                f"regions.feather: reference '{ref}' first region starts "
                f"at {int(starts[0])}, expected 0. Rebuild the index."
            )
        if ends[-1] != ref_len:
            raise ValueError(
                f"regions.feather: reference '{ref}' last region ends at "
                f"{int(ends[-1])}, expected {ref_len}. Rebuild the index."
            )
        if len(starts) > 1 and not np.array_equal(ends[:-1], starts[1:]):
            mismatch_index = int(np.flatnonzero(ends[:-1] != starts[1:])[0])
            raise ValueError(
                f"regions.feather: reference '{ref}' has a gap or overlap "
                f"between region ending at {int(ends[mismatch_index])} and region "
                f"starting at {int(starts[mismatch_index + 1])}. Rebuild the index."
            )

        signatures = sub["signature"].to_numpy(np.uint8, copy=False)
        if len(signatures) > 1 and (signatures[:-1] == signatures[1:]).any():
            same_index = int(np.flatnonzero(signatures[:-1] == signatures[1:])[0])
            raise ValueError(
                f"regions.feather: reference '{ref}' has adjacent regions with "
                f"the same signature at rows {int(sub.index[same_index])} and "
                f"{int(sub.index[same_index + 1])}. Rebuild the index."
            )

        expected_left = np.full(len(signatures), int(SIGNATURE_SENTINEL), dtype=np.uint8)
        expected_right = np.full(len(signatures), int(SIGNATURE_SENTINEL), dtype=np.uint8)
        expected_kind_left = np.full(len(signatures), int(BoundaryKind.NONE), dtype=np.uint8)
        expected_kind_right = np.full(len(signatures), int(BoundaryKind.NONE), dtype=np.uint8)
        if len(signatures) > 1:
            expected_left[1:] = signatures[:-1]
            expected_right[:-1] = signatures[1:]
            for offset in range(len(signatures) - 1):
                kind = classify_boundary_kind(int(signatures[offset]), int(signatures[offset + 1]))
                expected_kind_right[offset] = kind
                expected_kind_left[offset + 1] = kind

        checks = (
            ("left_signature", expected_left),
            ("right_signature", expected_right),
            ("boundary_kind_left", expected_kind_left),
            ("boundary_kind_right", expected_kind_right),
        )
        for column, expected in checks:
            actual = sub[column].to_numpy(np.uint8, copy=False)
            if not np.array_equal(actual, expected):
                bad_offset = int(np.flatnonzero(actual != expected)[0])
                bad_row = int(sub.index[bad_offset])
                raise ValueError(
                    f"regions.feather: region {bad_row} column {column!r} "
                    "does not match neighboring signatures. Rebuild the index."
                )

        expected_bf_left = sub["type"].to_numpy(np.uint8, copy=False) == int(RegionType.EXON)
        expected_bf_right = expected_bf_left.copy()
        expected_bf_left[0] = False
        expected_bf_right[-1] = False
        for column, expected in (
            ("boundary_flux_left", expected_bf_left),
            ("boundary_flux_right", expected_bf_right),
        ):
            actual = sub[column].to_numpy(bool, copy=False)
            if not np.array_equal(actual, expected):
                bad_offset = int(np.flatnonzero(actual != expected)[0])
                bad_row = int(sub.index[bad_offset])
                raise ValueError(
                    f"regions.feather: region {bad_row} column {column!r} "
                    "does not match the Phase 1 coarse bridge rule. Rebuild the index."
                )


def load_ref_lengths(path: str | Path) -> dict[str, int]:
    """Read ``ref_lengths.feather`` into an insertion-ordered dict.

    Iteration order matches the on-disk row order, which is the canonical
    ``ref_id`` assignment per calibration plan §2.11.
    """
    df = pd.read_feather(str(path))
    return {str(r): int(L) for r, L in zip(df["ref"], df["length"], strict=True)}
