"""Shared genomic-interval helpers for the sim package.

Two operations the capture sampler and the probe designer both need: merging overlapping genomic
intervals, and projecting one genomic block onto a transcript's spliced coordinates.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..transcript import Transcript
from ..types import Strand


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Sort and union overlapping/touching ``(start, end)`` intervals."""
    ordered = sorted(intervals)
    if not ordered:
        return []
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def project_genomic_block_to_transcript(
    transcript: Transcript,
    block_start: int,
    block_end: int,
) -> list[tuple[int, int]] | None:
    """Project one genomic block onto transcript coordinates (None if not fully exonic).

    Returns the transcript-space intervals the genomic ``[block_start, block_end)`` maps to,
    oriented to the transcript's strand. ``None`` if any base of the block falls outside the
    transcript's exons (i.e. the block is not fully contained in the spliced transcript).
    """
    tx_len = int(transcript.length or transcript.compute_length())
    consumed = 0
    projected: list[tuple[int, int]] = []
    mapped_bp = 0

    for exon in transcript.exons:
        exon_len = exon.end - exon.start
        overlap_start = max(block_start, exon.start)
        overlap_end = min(block_end, exon.end)
        if overlap_start < overlap_end:
            tx_start = consumed + (overlap_start - exon.start)
            tx_end = consumed + (overlap_end - exon.start)
            if transcript.strand == Strand.NEG:
                tx_start, tx_end = tx_len - tx_end, tx_len - tx_start
            projected.append((tx_start, tx_end))
            mapped_bp += overlap_end - overlap_start
        consumed += exon_len

    if mapped_bp != block_end - block_start:
        return None
    return sorted(projected)
