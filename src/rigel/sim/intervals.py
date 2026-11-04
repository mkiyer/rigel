"""Interval helpers shared by the capture sampler and the probe designer: merging intervals, and projecting genomic
blocks onto a transcript's spliced coordinates."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

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


def clip_blocks_to_transcript(
    transcript: Transcript, blocks: Sequence[tuple[int, int]]
) -> list[tuple[int, int]]:
    """The parts of genomic ``blocks`` inside ``transcript``'s exons, in strand-oriented transcript coordinates, block
    by block and exon by exon. A block need not be fully exonic: a probe hanging off an exon edge still binds the
    bases that are there."""
    tx_len = int(transcript.length or transcript.compute_length())
    out: list[tuple[int, int]] = []
    for block_start, block_end in blocks:
        consumed = 0
        for exon in transcript.exons:
            lo = max(block_start, exon.start)
            hi = min(block_end, exon.end)
            if lo < hi:
                tx_start = consumed + (lo - exon.start)
                tx_end = consumed + (hi - exon.start)
                if transcript.strand == Strand.NEG:
                    tx_start, tx_end = tx_len - tx_end, tx_len - tx_start
                out.append((tx_start, tx_end))
            consumed += exon.end - exon.start
    return out
