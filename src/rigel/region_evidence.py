"""
rigel.region_evidence — Fragment-to-region fractional counting (Stage 1 prototype).

Pure-Python implementation using pysam for BAM iteration and cgranges for
region overlap.  Produces two outputs:

1. **Region evidence** — a DataFrame with four float32 count columns
   (unspliced_pos, unspliced_neg, spliced_pos, spliced_neg) per region,
   accumulated via fractional overlap of aligned blocks.

2. **Fragment-length observation table** — a (region_id, frag_len) DataFrame
   for unspliced unique-mappers fully contained within a single region.

See ``docs/rigel_revision_impl/fragment_region_counting_plan.md`` for the
full design rationale.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from rigel.types import STRAND_NEG, STRAND_POS

logger = logging.getLogger(__name__)

# Count column indices
_COL_UNSPLICED_POS = 0
_COL_UNSPLICED_NEG = 1
_COL_SPLICED_POS = 2
_COL_SPLICED_NEG = 3

# CIGAR operation codes (from htslib/pysam)
_CIGAR_M = 0   # alignment match (sequence match or mismatch)
_CIGAR_I = 1   # insertion to the reference
_CIGAR_D = 2   # deletion from the reference
_CIGAR_N = 3   # skipped region from the reference (intron)
_CIGAR_S = 4   # soft clipping
_CIGAR_H = 5   # hard clipping
_CIGAR_EQ = 7  # sequence match
_CIGAR_X = 8   # sequence mismatch


# ---------------------------------------------------------------------------
# Fragment assembly
# ---------------------------------------------------------------------------

def assemble_fragment(
    records: list[pysam.AlignedSegment],
) -> tuple[list[tuple[int, int, int, int]], bool]:
    """Build merged aligned blocks from BAM records sharing a query name.

    Parameters
    ----------
    records : list[pysam.AlignedSegment]
        All BAM records for a single query name (name-sorted BAM).

    Returns
    -------
    blocks : list[tuple[int, int, int, int]]
        Merged aligned blocks as ``(ref_id, start, end, strand)``.
        Coordinates are 0-based half-open.  ``ref_id`` is the integer
        reference index from the BAM header.
    is_spliced : bool
        True if any record contains a CIGAR N (skip/intron) operation.
    """
    # Key: (ref_id, strand) → list of (start, end) intervals
    exon_dict: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    has_intron = False

    for rec in records:
        if rec.is_unmapped or rec.is_qcfail or rec.is_duplicate:
            continue
        if rec.cigartuples is None:
            continue

        ref = rec.reference_id

        # Strand: R2 gets flipped (Illumina PE converging-strand convention)
        if rec.is_read2:
            strand = STRAND_POS if rec.is_reverse else STRAND_NEG
        else:
            strand = STRAND_NEG if rec.is_reverse else STRAND_POS

        # Walk CIGAR to extract aligned blocks.
        # M/=/X and D ops are part of the same contiguous aligned block;
        # N (intron) splits blocks.
        pos = rec.reference_start
        block_start = None
        for op, length in rec.cigartuples:
            if op in (_CIGAR_M, _CIGAR_EQ, _CIGAR_X):
                if block_start is None:
                    block_start = pos
                pos += length
            elif op == _CIGAR_D:
                pos += length
            elif op == _CIGAR_N:
                if block_start is not None:
                    exon_dict[(ref, strand)].append((block_start, pos))
                    block_start = None
                has_intron = True
                pos += length
            # I, S, H: no reference movement
        if block_start is not None:
            exon_dict[(ref, strand)].append((block_start, pos))

    # Merge overlapping/adjacent intervals per (ref_id, strand) group
    blocks: list[tuple[int, int, int, int]] = []
    for (ref, strand), intervals in exon_dict.items():
        intervals.sort()
        merged_start, merged_end = intervals[0]
        for s, e in intervals[1:]:
            if s <= merged_end:
                merged_end = max(merged_end, e)
            else:
                blocks.append((ref, merged_start, merged_end, strand))
                merged_start, merged_end = s, e
        blocks.append((ref, merged_start, merged_end, strand))

    return blocks, has_intron


# ---------------------------------------------------------------------------
# Fractional counting
# ---------------------------------------------------------------------------

def count_fragment(
    blocks: list[tuple[int, int, int, int]],
    is_spliced: bool,
    region_cr,
    ref_names: tuple[str, ...] | list[str],
    counts: np.ndarray,
    fl_region_ids: list[int] | None = None,
    fl_frag_lens: list[int] | None = None,
) -> None:
    """Accumulate fractional region counts for one fragment (in-place).

    Parameters
    ----------
    blocks : list[tuple[int, int, int, int]]
        Merged aligned blocks ``(ref_id, start, end, strand)``.
    is_spliced : bool
        Whether the fragment has intron (CIGAR N) operations.
    region_cr
        cgranges index of the region partition.
    ref_names : tuple[str, ...] or list[str]
        Mapping from integer reference ID to reference name string
        (used for cgranges queries).
    counts : np.ndarray
        Shape ``(n_regions, 4)`` float64 accumulator, modified in-place.
    fl_region_ids, fl_frag_lens : list[int] | None
        If provided, fragment-length observations for unspliced fragments
        fully contained in a single region are appended.
    """
    if not blocks:
        return

    # Determine fragment strand — must be unanimous (non-chimeric)
    strand_set: set[int] = set()
    ref_set: set[int] = set()
    for ref, _s, _e, strand in blocks:
        strand_set.add(strand)
        ref_set.add(ref)

    # Chimeric: multiple strands or multiple references
    if len(strand_set) != 1 or len(ref_set) != 1:
        return
    frag_strand = strand_set.pop()
    if frag_strand not in (STRAND_POS, STRAND_NEG):
        return

    # Total aligned length (denominator)
    total_aligned = sum(end - start for _, start, end, _ in blocks)
    if total_aligned <= 0:
        return

    # Column index: unspliced=0,1  spliced=2,3  +POS=0,+NEG=1
    col = (_COL_SPLICED_POS if is_spliced else _COL_UNSPLICED_POS) + (
        0 if frag_strand == STRAND_POS else 1
    )

    # Accumulate fractional overlap per region
    for ref_id, start, end, _ in blocks:
        ref_name = ref_names[ref_id]
        for rgn_start, rgn_end, region_id in region_cr.overlap(ref_name, start, end):
            overlap_bp = min(end, rgn_end) - max(start, rgn_start)
            counts[region_id, col] += overlap_bp / total_aligned

    # Fragment-length tabulation: unspliced, single-region containment
    if fl_region_ids is not None and not is_spliced:
        frag_ref = ref_names[blocks[0][0]]
        genomic_start = min(b[1] for b in blocks)
        genomic_end = max(b[2] for b in blocks)
        hits = region_cr.overlap(frag_ref, genomic_start, genomic_end)
        if len(hits) == 1:
            fl_region_ids.append(hits[0][2])
            fl_frag_lens.append(genomic_end - genomic_start)


# ---------------------------------------------------------------------------
# Process a qname group
# ---------------------------------------------------------------------------

def _process_group(
    records: list[pysam.AlignedSegment],
    region_cr,
    ref_names: tuple[str, ...],
    counts: np.ndarray,
    fl_region_ids: list[int],
    fl_frag_lens: list[int],
) -> None:
    """Process all BAM records for a single query name."""
    # Determine NH (number of hits) — skip multi-mappers
    for rec in records:
        if rec.is_unmapped:
            continue
        nh = rec.get_tag("NH") if rec.has_tag("NH") else 1
        if nh > 1:
            return
        break  # only need one mapped record to check NH

    blocks, is_spliced = assemble_fragment(records)
    count_fragment(blocks, is_spliced, region_cr, ref_names, counts,
                   fl_region_ids, fl_frag_lens)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def count_region_evidence(
    bam_path: str | Path,
    index,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Count per-region fragment evidence from a name-sorted BAM.

    Parameters
    ----------
    bam_path : str or Path
        Path to a name-sorted BAM file.
    index : TranscriptIndex
        Loaded index with ``region_df`` and ``region_cr``.

    Returns
    -------
    region_counts : pd.DataFrame
        One row per region with columns ``region_id``, ``n_unspliced_pos``,
        ``n_unspliced_neg``, ``n_spliced_pos``, ``n_spliced_neg`` (float32).
    fl_table : pd.DataFrame
        Fragment-length observations for unspliced single-region fragments.
        Columns: ``region_id`` (int32), ``frag_len`` (int32).
    """
    n_regions = len(index.region_df)
    counts = np.zeros((n_regions, 4), dtype=np.float64)

    fl_region_ids: list[int] = []
    fl_frag_lens: list[int] = []

    n_groups = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        ref_names = bam.references
        current_qname: str | None = None
        current_group: list[pysam.AlignedSegment] = []

        for rec in bam:
            if rec.query_name != current_qname:
                if current_group:
                    _process_group(
                        current_group, index.region_cr, ref_names, counts,
                        fl_region_ids, fl_frag_lens,
                    )
                    n_groups += 1
                current_qname = rec.query_name
                current_group = [rec]
            else:
                current_group.append(rec)

        # Process last group
        if current_group:
            _process_group(
                current_group, index.region_cr, ref_names, counts,
                fl_region_ids, fl_frag_lens,
            )
            n_groups += 1

    logger.info(f"Processed {n_groups} query-name groups")

    # Build region evidence DataFrame (float32 for storage efficiency)
    region_counts = pd.DataFrame({
        "region_id": index.region_df["region_id"].values.astype(np.int32),
        "n_unspliced_pos": counts[:, _COL_UNSPLICED_POS].astype(np.float32),
        "n_unspliced_neg": counts[:, _COL_UNSPLICED_NEG].astype(np.float32),
        "n_spliced_pos": counts[:, _COL_SPLICED_POS].astype(np.float32),
        "n_spliced_neg": counts[:, _COL_SPLICED_NEG].astype(np.float32),
    })

    # Build fragment-length observation table
    fl_table = pd.DataFrame({
        "region_id": np.array(fl_region_ids, dtype=np.int32),
        "frag_len": np.array(fl_frag_lens, dtype=np.int32),
    })

    logger.info(
        f"Region evidence: {n_regions} regions, "
        f"{fl_table.shape[0]} FL observations"
    )
    return region_counts, fl_table
