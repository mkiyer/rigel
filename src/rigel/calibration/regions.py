"""rigel.calibration.regions — Region partition schema and emitter.

This module defines:

- ``RegionType`` — the 3-way partition (INTERGENIC / INTRON / EXON).
- ``RegionStrand`` — bitwise OR-able strand flags.
- ``RegionRecord`` — the on-disk row schema for ``regions.feather``.
- ``emit_regions(ref, layout)`` — generator that consumes the
  ``_iter_reference_layout`` stream from :mod:`rigel.index` and emits one
  ``RegionRecord`` per genomic region.
- ``load_regions(path)`` — load and dtype-coerce ``regions.feather``.
- ``validate_against_ref_lengths(region_df, ref_lengths)`` — enforce the
  load-time invariants from the calibration plan §2.2.

The genome is tiled per-reference into INTERGENIC spans (gaps between
genic spans) and per-genic-span EXON / INTRON regions, where EXON is the
union of all-strand exons within the genic span and INTRON is the
interstice. EXON wins at every position with any-strand exon coverage.

See ``docs/calibration/calibration_v6_plan.md`` §2.1 / §2.2.
"""

from __future__ import annotations

from enum import IntEnum, IntFlag
from pathlib import Path
from typing import TYPE_CHECKING, Iterator, Mapping, NamedTuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from ..index import _GenicSpan, _IntergenicSpan


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
    AMBIG = 3   # POS | NEG


# ---------------------------------------------------------------------------
# Row schema
# ---------------------------------------------------------------------------

class RegionRecord(NamedTuple):
    """One row of ``regions.feather``.

    Columns match the locked schema in the calibration plan §2.2.
    """

    region_id: int           # int64; assigned globally in genomic order
    ref_name: str
    start: int               # int64, 0-based inclusive
    end: int                 # int64, 0-based exclusive
    type: int                # uint8 (RegionType)
    strand: int              # uint8 (RegionStrand)
    tx_pos_bp: int           # int64, bp of (+)-strand transcript overlap
    tx_neg_bp: int           # int64, bp of (-)-strand transcript overlap
    exon_pos_bp: int         # int64, bp of (+)-strand exon overlap
    exon_neg_bp: int         # int64, bp of (-)-strand exon overlap
    boundary_flux_left: bool
    boundary_flux_right: bool


# ---------------------------------------------------------------------------
# Emitter
# ---------------------------------------------------------------------------

def emit_regions(
    ref: str,
    layout: Iterator["_IntergenicSpan | _GenicSpan"],
) -> Iterator[RegionRecord]:
    """Convert the per-reference layout stream into RegionRecord rows.

    Each ``_IntergenicSpan`` becomes exactly one INTERGENIC region.
    Each ``_GenicSpan`` is decomposed into a sequence of EXON and INTRON
    regions.

    ``region_id`` is left as ``-1`` here; the caller assigns globally
    monotonic ids after collecting all references.
    """
    # Local import to avoid a circular dependency with rigel.index.
    from ..index import _GenicSpan, _IntergenicSpan
    from ..types import Strand

    for span in layout:
        if isinstance(span, _IntergenicSpan):
            yield RegionRecord(
                region_id=-1,
                ref_name=ref,
                start=int(span.start),
                end=int(span.end),
                type=int(RegionType.INTERGENIC),
                strand=int(RegionStrand.NONE),
                tx_pos_bp=0,
                tx_neg_bp=0,
                exon_pos_bp=0,
                exon_neg_bp=0,
                boundary_flux_left=False,
                boundary_flux_right=False,
            )
            continue

        # Genic span: build per-strand exon mark arrays over the span.
        assert isinstance(span, _GenicSpan)
        s_start, s_end = int(span.start), int(span.end)
        L = s_end - s_start

        # Boolean coverage: (POS-exon, NEG-exon, POS-tx, NEG-tx) per bp.
        # int8 marks for any/either-strand exon, used to find region boundaries.
        # Note: we allocate L bytes per genic span; for typical loci this is
        # << 1MB. Mega-loci (HLA, etc.) may hit ~5MB which is still fine.
        import numpy as np
        exon_pos = np.zeros(L, dtype=np.bool_)
        exon_neg = np.zeros(L, dtype=np.bool_)
        tx_pos = np.zeros(L, dtype=np.bool_)
        tx_neg = np.zeros(L, dtype=np.bool_)

        for t in span.transcripts:
            if t.is_synthetic:
                continue
            t_lo = max(0, int(t.start) - s_start)
            t_hi = min(L, int(t.end) - s_start)
            if t_hi <= t_lo:
                continue
            if t.strand == Strand.POS:
                tx_pos[t_lo:t_hi] = True
            elif t.strand == Strand.NEG:
                tx_neg[t_lo:t_hi] = True
            for e in t.exons:
                e_lo = max(0, int(e.start) - s_start)
                e_hi = min(L, int(e.end) - s_start)
                if e_hi <= e_lo:
                    continue
                if t.strand == Strand.POS:
                    exon_pos[e_lo:e_hi] = True
                elif t.strand == Strand.NEG:
                    exon_neg[e_lo:e_hi] = True

        any_exon = exon_pos | exon_neg
        # Run-length decompose any_exon into alternating EXON/INTRON runs.
        # The first run is INTRON iff any_exon[0] is False.
        # Genic span by construction contains at least one exon, so any_exon
        # has at least one True position somewhere.
        # Find boundary indices via diff:
        diffs = np.diff(any_exon.astype(np.int8))
        boundary_idx = np.flatnonzero(diffs) + 1   # positions where state flips
        run_starts = np.concatenate(([0], boundary_idx))
        run_ends = np.concatenate((boundary_idx, [L]))

        # Pre-compute genic-span-relative cumulative sums for fast bp counts.
        cum = {
            "exon_pos": np.concatenate(([0], np.cumsum(exon_pos.astype(np.int64)))),
            "exon_neg": np.concatenate(([0], np.cumsum(exon_neg.astype(np.int64)))),
            "tx_pos": np.concatenate(([0], np.cumsum(tx_pos.astype(np.int64)))),
            "tx_neg": np.concatenate(([0], np.cumsum(tx_neg.astype(np.int64)))),
        }

        n_runs = run_starts.size
        for i in range(n_runs):
            r_lo = int(run_starts[i])
            r_hi = int(run_ends[i])
            is_exon = bool(any_exon[r_lo])

            ep = int(cum["exon_pos"][r_hi] - cum["exon_pos"][r_lo])
            en = int(cum["exon_neg"][r_hi] - cum["exon_neg"][r_lo])
            tp = int(cum["tx_pos"][r_hi] - cum["tx_pos"][r_lo])
            tn = int(cum["tx_neg"][r_hi] - cum["tx_neg"][r_lo])

            if is_exon:
                rtype = RegionType.EXON
                strand = RegionStrand.NONE
                if ep > 0:
                    strand |= RegionStrand.POS
                if en > 0:
                    strand |= RegionStrand.NEG
                # Boundary-flux: True iff the boundary touches an INTRON of
                # the same genic span (i.e., not the first/last run).
                bfl = i > 0
                bfr = i < n_runs - 1
            else:
                rtype = RegionType.INTRON
                strand = RegionStrand.NONE
                if tp > 0:
                    strand |= RegionStrand.POS
                if tn > 0:
                    strand |= RegionStrand.NEG
                bfl = False
                bfr = False

            yield RegionRecord(
                region_id=-1,
                ref_name=ref,
                start=s_start + r_lo,
                end=s_start + r_hi,
                type=int(rtype),
                strand=int(strand),
                tx_pos_bp=tp,
                tx_neg_bp=tn,
                exon_pos_bp=ep,
                exon_neg_bp=en,
                boundary_flux_left=bfl,
                boundary_flux_right=bfr,
            )


# ---------------------------------------------------------------------------
# Load + validation
# ---------------------------------------------------------------------------

#: Locked dtypes for the on-disk schema. Must match what
#: ``build_index_artifacts`` writes.
REGION_COLUMN_DTYPES: dict[str, type | np.dtype] = {
    "region_id": np.int64,
    "start": np.int64,
    "end": np.int64,
    "type": np.uint8,
    "strand": np.uint8,
    "tx_pos_bp": np.int64,
    "tx_neg_bp": np.int64,
    "exon_pos_bp": np.int64,
    "exon_neg_bp": np.int64,
    "boundary_flux_left": np.bool_,
    "boundary_flux_right": np.bool_,
}


def load_regions(path: str | Path) -> pd.DataFrame:
    """Read ``regions.feather`` and coerce dtypes to the locked schema.

    The returned DataFrame has its row index reset to ``region_id``
    (which is itself ``0..N-1``).  Raises :class:`ValueError` if the
    file is missing required columns.
    """
    df = pd.read_feather(str(path))
    required = set(RegionRecord._fields)
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"regions.feather at {path} is missing required columns: "
            f"{sorted(missing)}. Rebuild the index "
            f"(rigel index --fasta ... --gtf ...)."
        )
    for col, dt in REGION_COLUMN_DTYPES.items():
        if df[col].dtype != dt:
            df[col] = df[col].astype(dt)
    df["ref_name"] = df["ref_name"].astype("string")
    df = df.reset_index(drop=True)
    return df


def validate_against_ref_lengths(
    region_df: pd.DataFrame,
    ref_lengths: Mapping[str, int],
) -> None:
    """Enforce calibration plan §2.2 invariants.

    Raises :class:`ValueError` with an actionable message on the first
    violation. Checks (in order):

    1. ``region_id`` column equals row index.
    2. All region lengths are positive.
    3. Every ``ref_name`` exists in ``ref_lengths``.
    4. Per-reference: regions are sorted by start, non-overlapping, and
       cover ``[0, ref_lengths[ref])`` with no gaps.
    """
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
            i = int(np.flatnonzero(ends[:-1] != starts[1:])[0])
            raise ValueError(
                f"regions.feather: reference '{ref}' has a gap or overlap "
                f"between region ending at {int(ends[i])} and region "
                f"starting at {int(starts[i + 1])}. Rebuild the index."
            )


def load_ref_lengths(path: str | Path) -> dict[str, int]:
    """Read ``ref_lengths.feather`` into an insertion-ordered dict.

    Iteration order matches the on-disk row order, which is the canonical
    ``ref_id`` assignment per calibration plan §2.11.
    """
    df = pd.read_feather(str(path))
    return {str(r): int(L) for r, L in zip(df["ref"], df["length"], strict=True)}
