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


REGION_CORE_COLUMNS = ["region_id", "ref_name", "start", "end", "length", "signature"]
# Per-region mature-eligibility (spliced-exon coverage) on each strand — the structural input to the
# calibration mature/nascent message split (docs/calibration/three_component_mature_nascent_design.md §4).
REGION_MATURE_COLUMNS = ["mature_eligible_pos", "mature_eligible_neg"]
REGION_COLUMNS = REGION_CORE_COLUMNS + REGION_MATURE_COLUMNS

REGION_COLUMN_DTYPES: dict[str, type | np.dtype] = {
    "region_id": np.int64,
    "start": np.int64,
    "end": np.int64,
    "length": np.int64,
    "signature": np.uint8,
    "mature_eligible_pos": np.bool_,
    "mature_eligible_neg": np.bool_,
}

# ---------------------------------------------------------------------------
# Boundary partition (one row per region interface; the calibration boundary nodes)
# ---------------------------------------------------------------------------
# A reference with k regions has k+1 boundaries B0 R0 B1 … R(k-1) Bk (node_chain.build_node_chain): B0 at the
# reference start, Bk at the end, and internal Bi at region (i-1)'s end == region i's start. The row order
# (ref-major, i = 0..k in genomic order) aligns by construction with the accumulator's boundary indexing.
BOUNDARY_COLUMNS = [
    "boundary_id", "ref_name", "position",
    "is_tss", "is_tes", "is_splice_junction", "genomic_sj_strand",
]
BOUNDARY_COLUMN_DTYPES: dict[str, type | np.dtype] = {
    "boundary_id": np.int64,
    "position": np.int64,
    "is_tss": np.bool_,
    "is_tes": np.bool_,
    "is_splice_junction": np.bool_,
    "genomic_sj_strand": np.int8,  # Strand: 0 none / 1 POS / 2 NEG / 3 both (coincident opposite junctions)
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
    region_df = pd.DataFrame(rows, columns=REGION_CORE_COLUMNS)
    pos, neg = _compute_mature_eligible(region_df, transcripts)
    region_df["mature_eligible_pos"] = pos
    region_df["mature_eligible_neg"] = neg
    return _coerce_region_dtypes(region_df)


def _compute_mature_eligible(
    region_df: pd.DataFrame, transcripts: list[Transcript]
) -> tuple[np.ndarray, np.ndarray]:
    """Per-region mature-eligibility on each strand: does the region **overlap** an exon of a
    **multi-exon** (spliced) transcript on that strand?

    This is the structural statement of "spliced ⇒ mature"
    (`docs/calibration/three_component_mature_nascent_design.md` §4): a single-exon transcript's exon and a
    retained-intron region are correctly NOT mature-eligible (their RNA is nascent). **Overlap, not full
    coverage** — a region has a constant signature but its multi-exon-vs-single-exon coverage can vary within
    it (a multi-exon terminal exon merges with an adjacent single-exon exon when no intron separates them);
    overlap over-approximates safely (it never *under*-calls mature-eligibility, so it can never re-introduce a
    mature→intron bleed). Returns ``(pos, neg)`` bool arrays aligned to ``region_df`` rows.
    """
    n = len(region_df)
    pos = np.zeros(n, dtype=bool)
    neg = np.zeros(n, dtype=bool)
    if n == 0:
        return pos, neg

    # Merged multi-exon exon intervals per (ref, strand).
    exons_by: dict[tuple[str, int], list[tuple[int, int]]] = {}
    for tx in transcripts:
        if tx.is_synthetic or len(tx.exons) < 2:
            continue
        if tx.strand not in (Strand.POS, Strand.NEG):
            continue
        key = (str(tx.ref), int(tx.strand))
        lst = exons_by.setdefault(key, [])
        for exon in tx.exons:
            lst.append((int(exon.start), int(exon.end)))
    merged: dict[tuple[str, int], tuple[np.ndarray, np.ndarray]] = {}
    for key, intervals in exons_by.items():
        intervals.sort()
        starts: list[int] = []
        ends: list[int] = []
        for s0, e0 in intervals:
            if starts and s0 <= ends[-1]:
                ends[-1] = max(ends[-1], e0)
            else:
                starts.append(s0)
                ends.append(e0)
        merged[key] = (np.asarray(starts, np.int64), np.asarray(ends, np.int64))

    refs = region_df["ref_name"].to_numpy()
    r_start = region_df["start"].to_numpy(np.int64)
    r_end = region_df["end"].to_numpy(np.int64)
    for (ref, strand_int), (m_start, m_end) in merged.items():
        out = pos if strand_int == int(Strand.POS) else neg
        sel = np.flatnonzero(refs == ref)
        if sel.size == 0:
            continue
        rs, re = r_start[sel], r_end[sel]
        # Merged intervals are sorted & disjoint. Region [rs,re) overlaps one iff the last interval starting
        # STRICTLY before re ends after rs (earlier intervals end even earlier; later ones start at/after re).
        # ``side="left"`` excludes an interval that merely touches at m_start == re (half-open, no overlap).
        idx = np.searchsorted(m_start, re, side="left") - 1
        ok = idx >= 0
        hit = np.zeros(sel.size, dtype=bool)
        hit[ok] = m_end[idx[ok]] > rs[ok]
        out[sel] = hit
    return pos, neg


def build_boundary_partition(
    region_df: pd.DataFrame,
    transcripts: list[Transcript],
    ref_lengths: Mapping[str, int],
) -> pd.DataFrame:
    """Build the per-reference boundary partition — one row per region interface, carrying the annotation
    structural flags the calibration mature/nascent split reads
    (`docs/calibration/three_component_mature_nascent_design.md` §4).

    For a reference with ``k`` regions there are ``k+1`` boundaries (``node_chain``: B0 R0 … R(k-1) Bk); the
    boundary positions are ``[r0.start, r0.end, r1.end, …, r(k-1).end]``, ref-major in genomic order — aligning
    by construction with the accumulator's boundary indexing. Per boundary position ``p``:

    * ``is_splice_junction`` / ``genomic_sj_strand`` — ``p`` is a transcript **intron endpoint** (donor or
      acceptor); the motif strand is the transcript's strand (``3`` = coincident opposite-strand junctions).
    * ``is_tss`` / ``is_tes`` — ``p`` is a transcript's 5′ / 3′ terminus (orientation by strand; a ``+``
      transcript's TSS is its first-exon start, its TES its last-exon end; ``−`` is mirrored).
    """
    normalized_ref_lengths = {str(name): int(length) for name, length in ref_lengths.items()}
    # Per-ref annotation event positions (sets), by role and strand.
    sj_pos: dict[str, set[int]] = {}
    sj_neg: dict[str, set[int]] = {}
    tss: dict[str, set[int]] = {}
    tes: dict[str, set[int]] = {}
    for tx in transcripts:
        if tx.is_synthetic or not tx.exons or tx.strand not in (Strand.POS, Strand.NEG):
            continue
        ref = str(tx.ref)
        sj = sj_pos if tx.strand == Strand.POS else sj_neg
        s = sj.setdefault(ref, set())
        for intron_start, intron_end in tx.introns():
            s.add(int(intron_start))
            s.add(int(intron_end))
        first, last = int(tx.exons[0].start), int(tx.exons[-1].end)
        tss_p = first if tx.strand == Strand.POS else last
        tes_p = last if tx.strand == Strand.POS else first
        tss.setdefault(ref, set()).add(tss_p)
        tes.setdefault(ref, set()).add(tes_p)

    by_ref: dict[str, pd.DataFrame] = {
        ref: grp for ref, grp in region_df.groupby("ref_name", sort=False)
    }
    rows: list[dict[str, object]] = []
    for ref in normalized_ref_lengths:
        grp = by_ref.get(ref)
        if grp is None or len(grp) == 0:
            continue
        starts = grp["start"].to_numpy(np.int64, copy=False)
        ends = grp["end"].to_numpy(np.int64, copy=False)
        positions = np.empty(len(starts) + 1, dtype=np.int64)
        positions[:-1] = starts
        positions[-1] = ends[-1]
        sp, sn = sj_pos.get(ref, set()), sj_neg.get(ref, set())
        tss_r, tes_r = tss.get(ref, set()), tes.get(ref, set())
        for p in positions.tolist():
            in_pos, in_neg = p in sp, p in sn
            sj_strand = int(Strand.POS) * in_pos + int(Strand.NEG) * in_neg
            if in_pos and in_neg:
                sj_strand = int(Strand.AMBIGUOUS)
            rows.append(
                {
                    "boundary_id": -1,
                    "ref_name": ref,
                    "position": int(p),
                    "is_tss": p in tss_r,
                    "is_tes": p in tes_r,
                    "is_splice_junction": in_pos or in_neg,
                    "genomic_sj_strand": sj_strand,
                }
            )
    if not rows:
        return _coerce_boundary_dtypes(pd.DataFrame(columns=BOUNDARY_COLUMNS))
    for boundary_id, row in enumerate(rows):
        row["boundary_id"] = boundary_id
    return _coerce_boundary_dtypes(pd.DataFrame(rows, columns=BOUNDARY_COLUMNS))


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
            dtype = "string" if column == "ref_name" else REGION_COLUMN_DTYPES.get(column, object)
            region_df[column] = pd.Series(dtype=dtype)
    region_df = region_df.loc[:, REGION_COLUMNS]
    for column, dtype in REGION_COLUMN_DTYPES.items():
        if region_df[column].dtype != dtype:
            region_df[column] = region_df[column].astype(dtype)
    region_df["ref_name"] = region_df["ref_name"].astype("string")
    return region_df


def _coerce_boundary_dtypes(boundary_df: pd.DataFrame) -> pd.DataFrame:
    boundary_df = boundary_df.copy()
    for column in BOUNDARY_COLUMNS:
        if column not in boundary_df.columns:
            dtype = "string" if column == "ref_name" else BOUNDARY_COLUMN_DTYPES.get(column, object)
            boundary_df[column] = pd.Series(dtype=dtype)
    boundary_df = boundary_df.loc[:, BOUNDARY_COLUMNS]
    for column, dtype in BOUNDARY_COLUMN_DTYPES.items():
        if boundary_df[column].dtype != dtype:
            boundary_df[column] = boundary_df[column].astype(dtype)
    boundary_df["ref_name"] = boundary_df["ref_name"].astype("string")
    return boundary_df


def load_boundaries(path: str | Path) -> pd.DataFrame:
    """Read ``boundaries.feather`` and coerce dtypes to the boundary-partition schema. Raises
    :class:`ValueError` (rebuild the index) if required columns are missing (an old, pre-boundary index)."""
    df = pd.read_feather(str(path))
    missing = set(BOUNDARY_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(
            f"boundaries.feather at {path} is missing required columns: {sorted(missing)}. "
            "Rebuild the index (rigel index --fasta ... --gtf ...)."
        )
    df = df.loc[:, list(BOUNDARY_COLUMNS)]
    return _coerce_boundary_dtypes(df).reset_index(drop=True)


def validate_boundaries_against_regions(
    boundary_df: pd.DataFrame, region_df: pd.DataFrame
) -> None:
    """Enforce the boundary-partition invariant: each reference has ``regions + 1`` boundaries, and each
    boundary position equals the corresponding region interface (the ``[start…, last end]`` sequence). Raises
    :class:`ValueError` (rebuild the index) on the first violation."""
    missing = set(BOUNDARY_COLUMNS) - set(boundary_df.columns)
    if missing:
        raise ValueError(
            f"boundaries.feather is missing required columns: {sorted(missing)}. Rebuild the index."
        )
    reg_by_ref = {ref: grp for ref, grp in region_df.groupby("ref_name", sort=False)}
    bnd_by_ref = {ref: grp for ref, grp in boundary_df.groupby("ref_name", sort=False)}
    for ref, rgrp in reg_by_ref.items():
        bgrp = bnd_by_ref.get(ref)
        n_expected = len(rgrp) + 1
        if bgrp is None or len(bgrp) != n_expected:
            raise ValueError(
                f"boundaries.feather: reference '{ref}' has {0 if bgrp is None else len(bgrp)} boundaries, "
                f"expected {n_expected} (regions + 1). Rebuild the index."
            )
        starts = rgrp["start"].to_numpy(np.int64, copy=False)
        ends = rgrp["end"].to_numpy(np.int64, copy=False)
        expected_pos = np.empty(len(starts) + 1, dtype=np.int64)
        expected_pos[:-1] = starts
        expected_pos[-1] = ends[-1]
        if not np.array_equal(bgrp["position"].to_numpy(np.int64, copy=False), expected_pos):
            raise ValueError(
                f"boundaries.feather: reference '{ref}' boundary positions do not match the region "
                f"interfaces. Rebuild the index."
            )


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
    "REGION_CORE_COLUMNS",
    "REGION_MATURE_COLUMNS",
    "REGION_COLUMN_DTYPES",
    "BOUNDARY_COLUMNS",
    "BOUNDARY_COLUMN_DTYPES",
    "build_region_partition",
    "build_region_partition_arrays",
    "build_boundary_partition",
    "load_regions",
    "load_boundaries",
    "validate_against_ref_lengths",
    "validate_boundaries_against_regions",
    "load_ref_lengths",
]
