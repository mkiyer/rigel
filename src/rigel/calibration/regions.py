"""rigel.calibration.regions — Per-reference region partition.

This module defines the per-reference region partition consumed by the
fractional accumulator. The genome is tiled per-reference into maximal
intervals whose exon/intron annotation state is constant across the
interval; the accumulator then deposits per-fragment evidence into these
regions (and the boundaries between them) without distinguishing strand
or splicing classes at partition time.

Public API
----------
- ``build_region_partition(transcripts, ref_lengths)`` — event-sweep
  builder that returns a typed DataFrame with one row per region.
- ``build_region_partition_arrays(index)`` — flatten the partition to
  the ``(boundary_positions, ref_pos_offsets)`` ABI expected by
  :py:meth:`rigel.native.BamScanner.set_regions`.
- ``load_regions(path)`` — read ``regions.feather`` and coerce dtypes.
- ``validate_against_ref_lengths(region_df, ref_lengths)`` — enforce the
  partition invariants (tiling, ordering, non-negative lengths).
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


REGION_COLUMNS = ["region_id", "ref_name", "start", "end", "length"]

REGION_COLUMN_DTYPES: dict[str, type | np.dtype] = {
    "region_id": np.int64,
    "start": np.int64,
    "end": np.int64,
    "length": np.int64,
}


# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------


def _collect_breakpoints(
    transcripts: list[Transcript],
    ref_name: str,
    ref_length: int,
) -> np.ndarray:
    """Return sorted unique breakpoints for a reference, clipped to [0, L]."""
    breakpoints: list[int] = [0, int(ref_length)]
    for tx in transcripts:
        if tx.is_synthetic or not tx.exons:
            continue
        if str(tx.ref) != ref_name:
            continue
        if tx.strand not in (Strand.POS, Strand.NEG):
            continue
        for exon in tx.exons:
            breakpoints.append(int(exon.start))
            breakpoints.append(int(exon.end))
        for intron_start, intron_end in tx.introns():
            breakpoints.append(int(intron_start))
            breakpoints.append(int(intron_end))

    pts = np.asarray(breakpoints, dtype=np.int64)
    np.clip(pts, 0, int(ref_length), out=pts)
    pts = np.unique(pts)
    return pts


def build_region_partition(
    transcripts: list[Transcript],
    ref_lengths: Mapping[str, int],
) -> pd.DataFrame:
    """Build the per-reference region partition.

    For each reference, the partition tiles ``[0, ref_length)`` with
    non-overlapping intervals whose endpoints are the union of exon and
    intron boundaries from non-synthetic transcripts on that reference.
    References with zero length contribute no rows.
    """
    normalized_ref_lengths = {
        str(name): int(length) for name, length in ref_lengths.items()
    }

    for tx in transcripts:
        if tx.is_synthetic or not tx.exons:
            continue
        if str(tx.ref) not in normalized_ref_lengths:
            raise ValueError(
                f"Transcript {tx.t_id} has reference {tx.ref!r} "
                "not found in ref_lengths."
            )

    rows_ref: list[str] = []
    rows_start: list[int] = []
    rows_end: list[int] = []

    for ref_name, ref_length in normalized_ref_lengths.items():
        if ref_length <= 0:
            continue
        pts = _collect_breakpoints(transcripts, ref_name, ref_length)
        for i in range(len(pts) - 1):
            start = int(pts[i])
            end = int(pts[i + 1])
            if end <= start:
                continue
            rows_ref.append(ref_name)
            rows_start.append(start)
            rows_end.append(end)

    if not rows_ref:
        return _coerce_region_dtypes(pd.DataFrame(columns=REGION_COLUMNS))

    starts = np.asarray(rows_start, dtype=np.int64)
    ends = np.asarray(rows_end, dtype=np.int64)
    frame = pd.DataFrame(
        {
            "region_id": np.arange(len(rows_ref), dtype=np.int64),
            "ref_name": pd.array(rows_ref, dtype="string"),
            "start": starts,
            "end": ends,
            "length": ends - starts,
        }
    )
    return _coerce_region_dtypes(frame)


def build_region_partition_arrays(index) -> tuple[np.ndarray, np.ndarray]:
    """Flatten the index's region partition into BamScanner.set_regions arrays.

    For each reference ``f`` in ``index.ref_names``, the boundary positions
    are the sorted unique set ``[r0.start, r0.end, r1.end, ...]`` of the
    per-ref region rows (which tile the reference contiguously).
    References with zero regions contribute zero positions.

    Returns
    -------
    boundaries : int64[B_pos_total]
        Concatenated boundary positions, ref-major.
    ref_pos_offsets : int64[n_refs + 1]
        Offsets into ``boundaries``; ``ref_pos_offsets[n_refs] == B_pos_total``.
    """
    region_df = index.region_df
    ref_names = index.ref_names
    n_refs = len(ref_names)

    by_ref: dict[str, pd.DataFrame] = {
        ref: grp for ref, grp in region_df.groupby("ref_name", sort=False)
    }

    per_ref_positions: list[np.ndarray] = []
    offsets = np.zeros(n_refs + 1, dtype=np.int64)
    for i, ref in enumerate(ref_names):
        grp = by_ref.get(ref)
        if grp is None or len(grp) == 0:
            per_ref_positions.append(np.empty(0, dtype=np.int64))
            offsets[i + 1] = offsets[i]
            continue
        starts = grp["start"].to_numpy(np.int64, copy=False)
        ends = grp["end"].to_numpy(np.int64, copy=False)
        positions = np.empty(len(starts) + 1, dtype=np.int64)
        positions[:-1] = starts
        positions[-1] = ends[-1]
        per_ref_positions.append(positions)
        offsets[i + 1] = offsets[i] + positions.shape[0]

    boundaries = (
        np.concatenate(per_ref_positions)
        if per_ref_positions
        else np.empty(0, dtype=np.int64)
    )
    return boundaries, offsets


# ---------------------------------------------------------------------------
# Load + validation
# ---------------------------------------------------------------------------


def _coerce_region_dtypes(region_df: pd.DataFrame) -> pd.DataFrame:
    region_df = region_df.copy()
    for column in REGION_COLUMNS:
        if column not in region_df.columns:
            region_df[column] = pd.Series(
                dtype="string" if column == "ref_name" else object
            )
    region_df = region_df.loc[:, REGION_COLUMNS]
    for column, dtype in REGION_COLUMN_DTYPES.items():
        if region_df[column].dtype != dtype:
            region_df[column] = region_df[column].astype(dtype)
    region_df["ref_name"] = region_df["ref_name"].astype("string")
    return region_df


def _missing_schema_error(path: str | Path, missing: set[str]) -> ValueError:
    return ValueError(
        f"regions.feather at {path} is missing required columns: "
        f"{sorted(missing)}. Rebuild the index "
        "(rigel index --fasta ... --gtf ...)."
    )


def load_regions(path: str | Path) -> pd.DataFrame:
    """Read ``regions.feather`` and coerce dtypes to the partition schema.

    The returned DataFrame has its row index reset to ``region_id``
    (which is itself ``0..N-1``). Raises :class:`ValueError` if the file
    is missing required columns.
    """
    df = pd.read_feather(str(path))
    required = set(REGION_COLUMNS)
    missing = required - set(df.columns)
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
    3. Every ``ref_name`` exists in ``ref_lengths``.
    4. Per-reference: regions are sorted by start, non-overlapping, and
       cover ``[0, ref_lengths[ref])`` with no gaps.
    """
    missing = set(REGION_COLUMNS) - set(region_df.columns)
    if missing:
        raise _missing_schema_error("<in-memory>", missing)

    if not (region_df.index.to_numpy() == region_df["region_id"].to_numpy()).all():
        raise ValueError(
            "regions.feather: row index does not match 'region_id' column. "
            "Rebuild the index."
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
