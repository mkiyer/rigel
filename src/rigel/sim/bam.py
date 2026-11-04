"""BAM record construction and coordinate projection, shared by every simulator writer.

:func:`transcript_to_genomic_blocks` maps a strand-oriented transcript interval to the genomic exon blocks it
covers; :func:`blocks_to_cigar` writes the gaps between blocks as ``N`` operations; :func:`make_aligned_segment`
fills in the simulator's defaults. The flags and the base quality live here so the writers cannot disagree.
"""

from __future__ import annotations

import pysam

from ..transcript import Transcript
from ..types import Strand

__all__ = [
    "FLAG_PAIRED",
    "FLAG_PROPER_PAIR",
    "FLAG_REVERSE",
    "FLAG_MATE_REVERSE",
    "FLAG_READ1",
    "FLAG_READ2",
    "BASE_R1_FLAG",
    "BASE_R2_FLAG",
    "BASE_QUALITY",
    "blocks_to_cigar",
    "make_aligned_segment",
    "transcript_to_genomic_blocks",
]


FLAG_PAIRED = 0x1
FLAG_PROPER_PAIR = 0x2
FLAG_REVERSE = 0x10
FLAG_MATE_REVERSE = 0x20
FLAG_READ1 = 0x40
FLAG_READ2 = 0x80

BASE_R1_FLAG = FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_READ1
BASE_R2_FLAG = FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_READ2

#: The quality of every simulated base, in the FASTQ and the oracle BAM: Phred+33 ``I``, Q40.
BASE_QUALITY = "I"


def transcript_to_genomic_blocks(
    frag_start: int,
    frag_end: int,
    transcript: Transcript,
) -> list[tuple[int, int]]:
    """Map a strand-oriented transcript interval to genomic exon blocks, ascending."""
    exons = transcript.exons

    if transcript.strand == Strand.NEG:
        tx_len = sum(exon.end - exon.start for exon in exons)
        frag_start, frag_end = tx_len - frag_end, tx_len - frag_start

    blocks: list[tuple[int, int]] = []
    consumed = 0
    for exon in exons:
        exon_len = exon.end - exon.start
        exon_tx_start = consumed
        exon_tx_end = consumed + exon_len

        overlap_start = max(frag_start, exon_tx_start)
        overlap_end = min(frag_end, exon_tx_end)
        if overlap_start < overlap_end:
            offset_start = overlap_start - exon_tx_start
            offset_end = overlap_end - exon_tx_start
            blocks.append((exon.start + offset_start, exon.start + offset_end))

        consumed += exon_len
        if consumed >= frag_end:
            break

    return blocks


def blocks_to_cigar(blocks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Convert genomic blocks to pysam CIGAR tuples."""
    cigar: list[tuple[int, int]] = []
    for index, (block_start, block_end) in enumerate(blocks):
        if index > 0:
            prev_end = blocks[index - 1][1]
            intron_len = block_start - prev_end
            if intron_len > 0:
                cigar.append((pysam.CREF_SKIP, intron_len))
        match_len = block_end - block_start
        if match_len > 0:
            cigar.append((pysam.CMATCH, match_len))
    return cigar


def make_aligned_segment(
    header: pysam.AlignmentHeader,
    query_name: str,
    query_sequence: str,
    flag: int,
    reference_id: int,
    reference_start: int,
    cigar: list[tuple[int, int]],
    mate_reference_id: int,
    mate_reference_start: int,
    template_length: int,
    mapping_quality: int = 255,
    tags: list | None = None,
) -> pysam.AlignedSegment:
    """Build a pysam ``AlignedSegment`` with simulator-standard defaults."""
    segment = pysam.AlignedSegment(header)
    segment.query_name = query_name
    segment.query_sequence = query_sequence
    segment.flag = flag
    segment.reference_id = reference_id
    segment.reference_start = reference_start
    segment.cigar = cigar
    segment.mapping_quality = mapping_quality
    segment.query_qualities = pysam.qualitystring_to_array(BASE_QUALITY * len(query_sequence))
    segment.next_reference_id = mate_reference_id
    segment.next_reference_start = mate_reference_start
    segment.template_length = template_length
    if tags:
        segment.set_tags(tags)
    return segment
