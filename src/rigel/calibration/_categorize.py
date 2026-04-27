"""SRD Pass 0 — per-fragment categorization.

Walks a ``_FinalizedChunk`` once and assigns each fragment one of seven
categories (see :class:`FragmentCategory`). Categorization is fully
vectorized over the chunk's CSR arrays and the C++ scanner's already-computed
per-candidate ``exon_bp`` column.

The exon-fit test asks the only question that actually distinguishes RNA
from gDNA at the geometric level: *does at least one candidate transcript's
annotated exons cover (almost) the entire read?* Concretely an unspliced
fragment "fits" iff
``read_length - max(exon_bp[i] for i in candidates) <= exon_fit_tolerance_bp``.

Per-candidate ``intron_bp`` is intentionally **not** used — it counts
``transcript_span_bp - exon_bp``, which excludes any read bases that
overhang the transcript boundary into surrounding intergenic regions.
A small ``intron_bp`` does NOT prove an exonic fit; only a large
``exon_bp`` does.

Strand classification reuses the scanner's ``ambig_strand`` flag (mixed-strand
exon overlap) and the ``exon_strand`` field. tx_strand for sense/antisense
determination is taken from the candidate transcripts; when the consensus
``exon_strand`` is unknown or the overlap is mixed-strand, the fragment is
``UNSPLICED_EXONIC_AMBIG``.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np

from ..buffer import _FinalizedChunk
from ..splice import SpliceType
from ..types import Strand


class FragmentCategory(IntEnum):
    """SRD per-fragment category. Values are stable; persisted in arrays."""

    SPLICED = 0
    UNSPLICED_SENSE_EXONIC = 1
    UNSPLICED_ANTISENSE_EXONIC = 2
    UNSPLICED_EXONIC_AMBIG = 3
    EXON_INCOMPATIBLE = 4
    INTRONIC = 5
    INTERGENIC = 6


N_CATEGORIES = 7


@dataclass(frozen=True)
class CategorizedChunk:
    """Per-fragment SRD categorization output for one buffer chunk.

    All arrays have length ``chunk.size`` (parallel to the input chunk).

    Attributes
    ----------
    category : np.ndarray[uint8]
        :class:`FragmentCategory` value per fragment.
    frag_length : np.ndarray[int32]
        Representative fragment length per fragment (taken from the first
        transcript candidate's per-candidate length, or from
        ``read_length`` for fragments with no candidates).
    ref_id : np.ndarray[int32]
        Reference (chrom) integer id per fragment, derived from the first
        transcript candidate. ``-1`` for fragments with no candidates.
    g0 : np.ndarray[int32]
        Fragment genomic start (inclusive). Copied from
        ``chunk.genomic_start``. ``-1`` if unknown.
    g1 : np.ndarray[int32]
        Fragment genomic end (exclusive: ``g0 + genomic_footprint``).
        ``-1`` if unknown.
    """

    category: np.ndarray
    frag_length: np.ndarray
    ref_id: np.ndarray
    g0: np.ndarray
    g1: np.ndarray

    @property
    def n_per_category(self) -> np.ndarray:
        """Per-category count vector (length ``N_CATEGORIES``)."""
        return np.bincount(self.category, minlength=N_CATEGORIES).astype(np.int64)


def categorize_chunk(
    chunk: _FinalizedChunk,
    t_to_strand: np.ndarray,
    t_to_ref_id: np.ndarray,
    *,
    exon_fit_tolerance_bp: int = 5,
) -> CategorizedChunk:
    """Categorize every fragment in *chunk*.

    Parameters
    ----------
    chunk
        A finalized buffer chunk from ``FragmentBuffer.iter_chunks()``.
    t_to_strand
        Per-transcript strand array (``int8[n_transcripts]``) — values from
        :class:`rigel.types.Strand` (POS=1, NEG=2, NONE=0, ANY=3).
    t_to_ref_id
        Per-transcript reference id array (``int32[n_transcripts]``).
    exon_fit_tolerance_bp
        Tolerance (bp) for the exon-fit test. A fragment "fits" iff
        ``read_length - max(exon_bp) <= tol`` across its candidate
        transcripts.

    Returns
    -------
    CategorizedChunk
    """
    n = chunk.size
    cats = np.empty(n, dtype=np.uint8)
    ref_id = np.full(n, -1, dtype=np.int32)
    frag_length = np.zeros(n, dtype=np.int32)

    splice_type = chunk.splice_type
    exon_strand = chunk.exon_strand
    ambig_strand = chunk.ambig_strand
    t_offsets = chunk.t_offsets
    t_indices = chunk.t_indices
    intron_bp = chunk.intron_bp
    exon_bp = chunk.exon_bp
    cand_frag_lengths = chunk.frag_lengths

    n_cands = np.diff(t_offsets).astype(np.intp)
    has_cands = n_cands > 0

    # First-candidate index per fragment (-1 if none). Used to derive
    # ref_id and a representative fragment length.
    first_cand_idx = np.where(has_cands, t_offsets[:-1], -1)

    valid_first = first_cand_idx >= 0
    if valid_first.any():
        idx = first_cand_idx[valid_first]
        first_t = t_indices[idx]
        ref_id[valid_first] = t_to_ref_id[first_t]
        frag_length[valid_first] = cand_frag_lengths[idx]
    # For fragments with no transcript candidates, fall back to read_length
    # so Pool B has a sensible length for histogram construction.
    no_cands = ~valid_first
    if no_cands.any():
        frag_length[no_cands] = chunk.read_length[no_cands].astype(np.int32)

    # SPLICED — anything with N in the CIGAR.
    spliced_mask = (splice_type == SpliceType.SPLICED_ANNOT) | (
        splice_type == SpliceType.SPLICED_UNANNOT
    )
    cats[spliced_mask] = FragmentCategory.SPLICED

    # Unspliced fragments only from here on.
    unspliced_mask = ~spliced_mask

    # INTERGENIC = unspliced + no transcript candidates.
    intergenic_mask = unspliced_mask & ~has_cands
    cats[intergenic_mask] = FragmentCategory.INTERGENIC

    # Per-fragment max(exon_bp) across candidates. This is the only quantity
    # that actually answers "is this read covered by an exon?".
    #
    # NOTE on per-candidate semantics (see resolve_context.h:1055-1067):
    #   exon_bp[i]   = bp of the read covered by transcript i's annotated EXONS.
    #   intron_bp[i] = bp of the read inside transcript i's TRANSCRIPT SPAN
    #                  but NOT in any of its exons (i.e. introns of i).
    # Crucially: exon_bp + intron_bp <= read_length.  Bases of the read that
    # fall OUTSIDE transcript i's span entirely (overhanging the boundary
    # into a neighboring intron/intergenic region) are counted in NEITHER.
    # The earlier rule "min(intron_bp) <= tol" therefore mis-classified any
    # read with a small overhang as exonic — including unique gDNA reads
    # that splash a few bp past a transcript edge.
    #
    # The correct test for "this read fits in an exon of at least one
    # candidate" is:  read_length - max(exon_bp) <= tol.
    has_cands_unspliced = unspliced_mask & has_cands
    if has_cands_unspliced.any():
        starts = t_offsets[:-1]
        safe_starts = np.where(n_cands > 0, starts, 0).astype(np.intp)
        if t_indices.size > 0:
            max_exon = np.maximum.reduceat(exon_bp, safe_starts)
        else:
            max_exon = np.zeros(n, dtype=exon_bp.dtype)

        # Per-fragment read length (already populated for has_cands frags).
        rl = chunk.read_length.astype(np.int32, copy=False)
        tol = int(exon_fit_tolerance_bp)

        # exon_fit: ≥1 candidate transcript's annotated exons cover the
        # read (within `tol` bp of slop).
        exon_fit = has_cands_unspliced & ((rl - max_exon) <= tol)

        # INTRONIC = no exon overlap at all on any candidate. The read
        # overlaps a transcript span but lands entirely in introns (and
        # possibly overhanging into intergenic territory).
        intronic_mask = has_cands_unspliced & ~exon_fit & (max_exon == 0)
        cats[intronic_mask] = FragmentCategory.INTRONIC

        # EXON_INCOMPATIBLE = some exon overlap but the union of exonic
        # bp on every candidate falls short of the read by more than tol.
        # Includes reads that overhang transcript edges into surrounding
        # genomic territory.
        incompat_mask = has_cands_unspliced & ~exon_fit & (max_exon > 0)
        cats[incompat_mask] = FragmentCategory.EXON_INCOMPATIBLE

        # UNSPLICED_EXONIC{_SENSE,_ANTISENSE,_AMBIG}.
        # tx_strand for sense/antisense determination is taken via the
        # consensus: if all candidate transcripts agree on a strand and
        # ambig_strand == 0, use that; otherwise AMBIG.
        if exon_fit.any():
            # Consensus tx_strand per fragment: min == max across candidates.
            # Reuse reduceat pattern.
            t_strand_per_cand = t_to_strand[t_indices].astype(np.int8)
            if t_indices.size > 0:
                min_str = np.minimum.reduceat(t_strand_per_cand, safe_starts)
                max_str = np.maximum.reduceat(t_strand_per_cand, safe_starts)
            else:
                min_str = np.zeros(n, dtype=np.int8)
                max_str = np.zeros(n, dtype=np.int8)
            consensus_strand = np.where(min_str == max_str, min_str, np.int8(0))

            # An unambiguous tx_strand is POS or NEG and matches the C++
            # scanner's mixed-strand flag (ambig_strand == 0).
            unambig_tx = (
                exon_fit
                & (ambig_strand == 0)
                & ((consensus_strand == int(Strand.POS)) | (consensus_strand == int(Strand.NEG)))
            )

            # Aligned read strand: from chunk.exon_strand (POS/NEG/NONE).
            sense_mask = (
                unambig_tx
                & ((exon_strand == int(Strand.POS)) | (exon_strand == int(Strand.NEG)))
                & (exon_strand == consensus_strand)
            )
            antisense_mask = (
                unambig_tx
                & ((exon_strand == int(Strand.POS)) | (exon_strand == int(Strand.NEG)))
                & (exon_strand != consensus_strand)
            )
            ambig_mask = exon_fit & ~sense_mask & ~antisense_mask

            cats[sense_mask] = FragmentCategory.UNSPLICED_SENSE_EXONIC
            cats[antisense_mask] = FragmentCategory.UNSPLICED_ANTISENSE_EXONIC
            cats[ambig_mask] = FragmentCategory.UNSPLICED_EXONIC_AMBIG

    return CategorizedChunk(
        category=cats,
        frag_length=frag_length,
        ref_id=ref_id,
        g0=chunk.genomic_start.astype(np.int32, copy=False),
        g1=(chunk.genomic_start.astype(np.int32) + chunk.genomic_footprint.astype(np.int32)),
    )


def build_t_to_ref_id(t_df) -> np.ndarray:
    """Build per-transcript reference (chrom) integer id array.

    Parameters
    ----------
    t_df
        :attr:`rigel.index.TranscriptIndex.t_df` (the transcript DataFrame
        with a categorical ``ref`` column).

    Returns
    -------
    np.ndarray[int32]
        ``t_to_ref_id[t_idx]`` = integer id of the chromosome
        (categorical code) for transcript ``t_idx``.
    """
    refs = t_df["ref"]
    # Pandas categorical => integer codes are stable per-load.
    if hasattr(refs, "cat"):
        return refs.cat.codes.to_numpy(dtype=np.int32, copy=True)
    # Fallback: factorize.
    codes, _ = np.unique(np.asarray(refs), return_inverse=True)
    return codes.astype(np.int32)
