"""Test-only scoring helpers extracted from rigel.scoring.

These functions are used by unit tests to verify scoring logic
without constructing a full FragmentScorer.  They are NOT part of
the production pipeline.
"""

import math

import numpy as np

from rigel.scoring import (
    GDNA_SPLICE_PENALTIES,
    LOG_HALF,
    LOG_SAFE_FLOOR,
)
from rigel.types import STRAND_NEG


def genomic_to_transcript_pos(
    genomic_pos: int,
    exon_intervals: np.ndarray,
    strand: int,
    transcript_length: int,
) -> int:
    """Map a genomic position to a transcript-relative 5'→3' coordinate.

    Walks through sorted exon intervals ``(start, end)`` accumulating
    spliced offset until the exon containing *genomic_pos* is found.
    Positions falling in introns are mapped to the preceding exon
    boundary.  Negative-strand transcripts are flipped so that
    offset 0 is the 5' end.

    Parameters
    ----------
    genomic_pos : int
        Genomic coordinate to map.
    exon_intervals : np.ndarray
        ``(n_exons, 2)`` int32 array of ``[start, end)`` intervals,
        sorted by genomic start position.
    strand : int
        Transcript strand (``Strand.POS`` or ``Strand.NEG``).
    transcript_length : int
        Total spliced exonic length of the transcript.

    Returns
    -------
    int
        Transcript-relative position in ``[0, transcript_length]``.
    """
    offset = 0
    n = len(exon_intervals)
    for i in range(n):
        ex_start = int(exon_intervals[i, 0])
        ex_end = int(exon_intervals[i, 1])
        if genomic_pos < ex_start:
            # Position is upstream of (or in an intron before) this exon
            break
        if genomic_pos < ex_end:
            # Position is inside this exon
            offset += genomic_pos - ex_start
            break
        offset += ex_end - ex_start
    else:
        # Position is past the last exon
        offset = transcript_length

    offset = max(0, min(offset, transcript_length))
    if strand == STRAND_NEG:
        offset = transcript_length - offset
    return offset


def score_gdna_standalone(
    exon_strand: int,
    splice_type: int,
    frag_length: int,
    frag_length_models,
    gdna_splice_penalties: dict | None = None,
) -> float:
    """Compute gDNA log-likelihood for a fragment.

    gDNA has no strand bias: strand probability is always 0.5.

    Intended for unit tests and external callers that do not construct
    a ``FragmentScorer``.
    """
    penalties = gdna_splice_penalties or GDNA_SPLICE_PENALTIES
    splice_pen = penalties.get(splice_type, 1.0)

    log_p_insert = (
        frag_length_models.gdna_model.log_likelihood(frag_length)
        if frag_length > 0
        else 0.0
    )

    return (
        LOG_HALF
        + log_p_insert
        + math.log(max(splice_pen, LOG_SAFE_FLOOR))
    )
