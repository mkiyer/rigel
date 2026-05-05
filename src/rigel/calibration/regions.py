"""rigel.calibration.regions — Region partition schema and emitter.

This module defines:

- ``RegionType`` — the 3-way partition (INTERGENIC / INTRON / EXON).
- ``RegionStrand`` — bitwise OR-able strand flags.
- ``RegionRecord`` — the on-disk row schema for ``regions.feather``.
- ``emit_regions(ref, layout)`` — generator that consumes the
  ``_iter_reference_layout`` stream from :mod:`rigel.index` and emits one
  ``RegionRecord`` per genomic region.

The genome is tiled per-reference into INTERGENIC spans (gaps between
genic spans) and per-genic-span EXON / INTRON regions, where EXON is the
union of all-strand exons within the genic span and INTRON is the
interstice. EXON wins at every position with any-strand exon coverage.

See ``docs/calibration/calibration_v6_plan.md`` §2.1 / §2.2.
"""

from __future__ import annotations

from enum import IntEnum, IntFlag
from typing import TYPE_CHECKING, Iterator, NamedTuple

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
