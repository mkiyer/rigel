"""rigel.calibration.regions — Per-reference region partition.

The genome is tiled per-reference into *regions*: maximal intervals over
which the 4-bit annotation :mod:`signature <rigel.calibration.signature>`
(exon/intron × strand) is constant. By construction **a region and each of
its neighbours have different signatures** — adjacent segments with identical
signatures are merged — so the boundaries between regions sit exactly at
signature transitions. The accumulator deposits per-fragment evidence into
these regions and the boundaries between them.

Public API
----------
- ``build_region_partition(transcripts, ref_lengths)`` — event-sweep builder
  returning a typed DataFrame with one row per region (carrying its
  ``signature``).
- ``build_region_partition_arrays(index)`` — flatten the partition to the
  ``(boundary_positions, ref_pos_offsets, region_types)`` ABI expected by
  :py:meth:`rigel.native.BamScanner.set_regions` (``region_types`` carries the
  per-region coarse type for the gDNA FL pools, PR 4c).
- ``load_regions(path)`` — read ``regions.feather`` and coerce dtypes.
- ``validate_against_ref_lengths(region_df, ref_lengths)`` — enforce the
  partition invariants (tiling, ordering, signature range, neighbour-differs).
- ``load_ref_lengths(path)`` — read ``ref_lengths.feather`` into an
  insertion-ordered dict keyed by reference name.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

from ..transcript import Transcript
from ..types import Strand
from .signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    N_SIGNATURES,
    coarse_type_array,
    pack_signature,
)


REGION_COLUMNS = ["region_id", "ref_name", "start", "end", "length", "signature"]

REGION_COLUMN_DTYPES: dict[str, type | np.dtype] = {
    "region_id": np.int64,
    "start": np.int64,
    "end": np.int64,
    "length": np.int64,
    "signature": np.uint8,
}


# ---------------------------------------------------------------------------
# Builder (event sweep with adjacent-equal-signature merge)
# ---------------------------------------------------------------------------


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
    """Append +1/-1 delta events for one validated half-open interval."""
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


def _build_reference_rows(
    ref_name: str,
    ref_length: int,
    events: list[tuple[int, int, int]],
) -> list[dict[str, object]]:
    """Sweep one reference's events into merged constant-signature segments."""
    if ref_length < 0:
        raise ValueError(f"Reference {ref_name!r} has negative length {ref_length}.")
    if ref_length == 0:
        return []

    events.sort(key=lambda event: event[0])
    rows: list[dict[str, object]] = []
    counts = {BIT_INTRON_POS: 0, BIT_INTRON_NEG: 0, BIT_EXON_POS: 0, BIT_EXON_NEG: 0}
    previous_position = 0
    event_index = 0
    n_events = len(events)

    def emit_segment(segment_start: int, segment_end: int) -> None:
        if segment_end <= segment_start:
            return
        signature = pack_signature(
            intron_pos=counts[BIT_INTRON_POS] > 0,
            intron_neg=counts[BIT_INTRON_NEG] > 0,
            exon_pos=counts[BIT_EXON_POS] > 0,
            exon_neg=counts[BIT_EXON_NEG] > 0,
        )
        # Merge: extend the previous region when it abuts and shares a signature
        # (this is what enforces the neighbour-differs invariant).
        if rows and rows[-1]["end"] == segment_start and rows[-1]["signature"] == signature:
            rows[-1]["end"] = segment_end
            rows[-1]["length"] = segment_end - int(rows[-1]["start"])
            return
        rows.append(
            {
                "region_id": -1,
                "ref_name": ref_name,
                "start": segment_start,
                "end": segment_end,
                "length": segment_end - segment_start,
                "signature": signature,
            }
        )

    while event_index < n_events:
        position = events[event_index][0]
        if position > previous_position:
            emit_segment(previous_position, position)
            previous_position = position
        while event_index < n_events and events[event_index][0] == position:
            _, bit, delta = events[event_index]
            counts[bit] += delta
            event_index += 1

    if ref_length > previous_position:
        emit_segment(previous_position, ref_length)
    return rows


def build_region_partition(
    transcripts: list[Transcript],
    ref_lengths: Mapping[str, int],
) -> pd.DataFrame:
    """Build the per-reference region partition.

    For each reference, exon and intron interval events from non-synthetic
    transcripts are swept into maximal constant-signature segments; adjacent
    segments with identical signatures are merged. ``region_id`` is assigned
    globally in genomic order (``ref_lengths`` iteration order, then start).
    References with zero length contribute no rows.
    """
    normalized_ref_lengths = {str(name): int(length) for name, length in ref_lengths.items()}

    events_by_ref: dict[str, list[tuple[int, int, int]]] = {
        ref_name: [] for ref_name in normalized_ref_lengths
    }

    for tx in transcripts:
        if tx.is_synthetic or not tx.exons:
            continue
        ref_name = str(tx.ref)
        if ref_name not in normalized_ref_lengths:
            raise ValueError(
                f"Transcript {tx.t_id} has reference {tx.ref!r} not found in ref_lengths."
            )
        if tx.strand == Strand.POS:
            exon_bit, intron_bit = BIT_EXON_POS, BIT_INTRON_POS
        elif tx.strand == Strand.NEG:
            exon_bit, intron_bit = BIT_EXON_NEG, BIT_INTRON_NEG
        else:
            continue

        ref_length = normalized_ref_lengths[ref_name]
        events = events_by_ref[ref_name]
        for exon in tx.exons:
            _add_interval_events(
                events,
                ref_name=ref_name,
                ref_length=ref_length,
                transcript_id=tx.t_id,
                start=int(exon.start),
                end=int(exon.end),
                bit=exon_bit,
            )
        for intron_start, intron_end in tx.introns():
            _add_interval_events(
                events,
                ref_name=ref_name,
                ref_length=ref_length,
                transcript_id=tx.t_id,
                start=int(intron_start),
                end=int(intron_end),
                bit=intron_bit,
            )

    rows: list[dict[str, object]] = []
    for ref_name, ref_length in normalized_ref_lengths.items():
        rows.extend(_build_reference_rows(ref_name, ref_length, events_by_ref[ref_name]))

    if not rows:
        return _coerce_region_dtypes(pd.DataFrame(columns=REGION_COLUMNS))

    for region_id, row in enumerate(rows):
        row["region_id"] = region_id
    return _coerce_region_dtypes(pd.DataFrame(rows, columns=REGION_COLUMNS))


def build_region_partition_arrays(index) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flatten the index's region partition into BamScanner.set_regions arrays.

    For each reference ``f`` in ``index.ref_names``, the boundary positions
    are the sorted set ``[r0.start, r0.end, r1.end, ...]`` of the per-ref
    region rows (which tile the reference contiguously). References with zero
    regions contribute zero positions.

    Returns
    -------
    boundaries : int64[B_pos_total]
        Concatenated boundary positions, ref-major.
    ref_pos_offsets : int64[n_refs + 1]
        Offsets into ``boundaries``; ``ref_pos_offsets[n_refs] == B_pos_total``.
    region_types : uint8[R_total]
        Per-region coarse type (0=intergenic, 1=intron, 2=exon), ref-major in
        the same region order as the partition — the gDNA FL-pool region axis
        (PR 4c). Aligns 1:1 with the accumulator's regions.
    """
    region_df = index.region_df
    ref_names = index.ref_names
    n_refs = len(ref_names)

    by_ref: dict[str, pd.DataFrame] = {
        ref: grp for ref, grp in region_df.groupby("ref_name", sort=False)
    }

    per_ref_positions: list[np.ndarray] = []
    per_ref_types: list[np.ndarray] = []
    offsets = np.zeros(n_refs + 1, dtype=np.int64)
    for i, ref in enumerate(ref_names):
        grp = by_ref.get(ref)
        if grp is None or len(grp) == 0:
            per_ref_positions.append(np.empty(0, dtype=np.int64))
            per_ref_types.append(np.empty(0, dtype=np.uint8))
            offsets[i + 1] = offsets[i]
            continue
        starts = grp["start"].to_numpy(np.int64, copy=False)
        ends = grp["end"].to_numpy(np.int64, copy=False)
        positions = np.empty(len(starts) + 1, dtype=np.int64)
        positions[:-1] = starts
        positions[-1] = ends[-1]
        per_ref_positions.append(positions)
        per_ref_types.append(coarse_type_array(grp["signature"].to_numpy()))
        offsets[i + 1] = offsets[i] + positions.shape[0]

    boundaries = (
        np.concatenate(per_ref_positions) if per_ref_positions else np.empty(0, dtype=np.int64)
    )
    region_types = np.concatenate(per_ref_types) if per_ref_types else np.empty(0, dtype=np.uint8)
    return boundaries, offsets, region_types


# ---------------------------------------------------------------------------
# Load + validation
# ---------------------------------------------------------------------------


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


def _missing_schema_error(path: str | Path, missing: set[str]) -> ValueError:
    if "signature" in missing:
        return ValueError(
            f"regions.feather at {path} is missing 'signature' (pre-signature index). "
            "Rebuild the index (rigel index --fasta ... --gtf ...)."
        )
    return ValueError(
        f"regions.feather at {path} is missing required columns: "
        f"{sorted(missing)}. Rebuild the index (rigel index --fasta ... --gtf ...)."
    )


def load_regions(path: str | Path) -> pd.DataFrame:
    """Read ``regions.feather`` and coerce dtypes to the partition schema.

    The returned DataFrame has its row index reset to ``region_id`` (which is
    itself ``0..N-1``). Raises :class:`ValueError` if required columns are
    missing.
    """
    df = pd.read_feather(str(path))
    missing = set(REGION_COLUMNS) - set(df.columns)
    if missing:
        raise _missing_schema_error(path, missing)
    df = df.loc[:, list(REGION_COLUMNS)]
    return _coerce_region_dtypes(df).reset_index(drop=True)


def validate_against_ref_lengths(
    region_df: pd.DataFrame,
    ref_lengths: Mapping[str, int],
) -> None:
    """Enforce the region partition invariants.

    Raises :class:`ValueError` with an actionable message on the first
    violation. Checks (in order):

    1. ``region_id`` column equals row index.
    2. All region lengths are positive and match ``end - start``.
    3. Every ``signature`` is a valid 4-bit value.
    4. Every ``ref_name`` exists in ``ref_lengths``.
    5. Per-reference: regions are sorted by start, non-overlapping, cover
       ``[0, ref_lengths[ref])`` with no gaps, and **no two adjacent regions
       share a signature** (the merge invariant).
    """
    missing = set(REGION_COLUMNS) - set(region_df.columns)
    if missing:
        raise _missing_schema_error("<in-memory>", missing)

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

    signatures = region_df["signature"].to_numpy(np.uint8, copy=False)
    if (signatures >= N_SIGNATURES).any():
        bad = int(np.flatnonzero(signatures >= N_SIGNATURES)[0])
        raise ValueError(
            f"regions.feather: region {bad} has invalid signature "
            f"{int(signatures[bad])}. Rebuild the index."
        )

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

        ref_signatures = sub["signature"].to_numpy(np.uint8, copy=False)
        if len(ref_signatures) > 1 and (ref_signatures[:-1] == ref_signatures[1:]).any():
            same = int(np.flatnonzero(ref_signatures[:-1] == ref_signatures[1:])[0])
            raise ValueError(
                f"regions.feather: reference '{ref}' has adjacent regions with the "
                f"same signature at rows {int(sub.index[same])} and "
                f"{int(sub.index[same + 1])} (should have been merged). Rebuild the index."
            )


def load_ref_lengths(path: str | Path) -> dict[str, int]:
    """Read ``ref_lengths.feather`` into an insertion-ordered dict.

    Iteration order matches the on-disk row order, which is the canonical
    ``ref_id`` assignment.
    """
    df = pd.read_feather(str(path))
    return {str(r): int(L) for r, L in zip(df["ref"], df["length"], strict=True)}


__all__ = [
    "REGION_COLUMNS",
    "REGION_COLUMN_DTYPES",
    "build_region_partition",
    "build_region_partition_arrays",
    "load_regions",
    "validate_against_ref_lengths",
    "load_ref_lengths",
]
