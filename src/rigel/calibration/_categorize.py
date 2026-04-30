"""SRD v3 Phase 1 — per-fragment geometric categorization with
**transcript-frame** strand labeling.

Walks a ``_FinalizedChunk`` once and labels each unspliced fragment with one
of four geometric categories (see :class:`FragmentCategory`) and a
transcript-frame strand label (see :class:`FragmentStrand`).  Fully
vectorised over the chunk's strand-aware overlap columns produced by the
C++ scanner:

* ``exon_bp_pos`` / ``exon_bp_neg``  — bp covered by ANY (+/−)-strand
  transcript's annotated exons.
* ``tx_bp_pos``  / ``tx_bp_neg``    — bp covered by ANY (+/−)-strand
  transcript's full span (transcript-tile interval).
* ``exon_strand``                   — the read's inferred genomic
  strand (POS/NEG); AMBIGUOUS or UNKNOWN (R1/R2 disagree or
  unresolved) fragments are dropped from labeling.

Per-strand intronic overlap is *derived* in Python as
``intron_bp_strand = max(tx_bp_strand − exon_bp_strand, 0)`` because the
cgranges index has no explicit ``ITYPE_INTRON`` interval type.

The categorization is geometric and library-agnostic.  The transcript-
frame strand label compares the read's *protocol-corrected* template
strand against the strand(s) of overlapping transcripts.  The C++
scanner stores ``chunk.exon_strand`` as the read's genomic strand
(after the dUTP-style R2 flip); the categorizer applies an additional
flip when the library is R1-antisense (``read1_sense=False``) so that
SENSE always means "fragment derives from this transcript's template":

* SENSE     — fragment template strand matches the (single)
              overlapping transcript strand.
* ANTISENSE — fragment template strand is opposite to the (single)
              overlapping transcript strand.
* AMBIG     — transcripts on **both** strands overlap (the fragment is
              sense to one and antisense to another); cannot classify.

For INTERGENIC fragments (no transcript overlap) the same enum values
are reused as a pure naming convention — SENSE = template on (+)
genomic strand, ANTISENSE = template on (−).  Downstream consumers
should treat INTERGENIC labels as a 50/50 sanity diagnostic on
genomic-strand assignment, not as biological sense/antisense.

Categories
----------
Per `docs/calibration/srd_v2_phase2plus_handoff.md` §2 / §4.

* ``INTERGENIC``      — ``tx_bp_pos == 0`` and ``tx_bp_neg == 0``
* ``INTRONIC``        — ``exon_bp_pos == 0`` and ``exon_bp_neg == 0``
                        and not INTERGENIC
* ``EXON_CONTAINED``  — there exists at least one candidate transcript
                        ``T`` whose exon set covers all read blocks of
                        the fragment, i.e.
                        ``read_length - max_T(exon_bp[T])
                        <= exon_fit_tolerance_bp``.  Per-transcript
                        ``exon_bp[T]`` (CSR ``exon_bp`` column) sums
                        per-read-block coverage by ``T``'s annotated
                        exons, so ``exon_bp[T] == read_length``
                        means every block lies inside an exon of ``T``
                        (and, for an UNSPLICED fragment, the inner
                        mate gap necessarily lies in the same exon).
                        Aggregating across transcripts via the scalar
                        ``exon_bp_pos / _neg`` columns is unsafe — those
                        sum per-cgranges-hit per-block bp and double-count
                        overlapping exons of distinct transcripts; the
                        check must be per-T to honour the user-stated
                        rule "contained in an exon of *any* single
                        transcript".
* ``EXON_INCOMPATIBLE`` — none of the above

Pre-filter
----------
Categorization is only meaningful for **unique-mapper, truly-unspliced
fragments with a known read strand** — ``num_hits == 1`` AND
``splice_type == UNSPLICED`` AND ``exon_strand ∈ {POS, NEG}``.  All
other fragments (multi-mapper, spliced, AMBIGUOUS read strand) are
held out of calibration entirely.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np

from ..buffer import _FinalizedChunk
from ..splice import SPLICE_UNSPLICED
from ..types import STRAND_NEG, STRAND_POS


class FragmentCategory(IntEnum):
    """SRD v2 per-fragment geometric category.  Values are stable;
    persisted as the row index into ``CalibrationResult.category_counts``.
    """

    INTERGENIC = 0
    INTRONIC = 1
    EXON_CONTAINED = 2
    EXON_INCOMPATIBLE = 3


N_CATEGORIES = len(FragmentCategory)


class FragmentStrand(IntEnum):
    """SRD v3 transcript-frame strand label of a categorized fragment.

    Values are stable; persisted as the column index into
    ``CalibrationResult.category_counts``.

    SENSE / ANTISENSE are defined relative to the strand of the
    transcript(s) overlapping the fragment.  For INTERGENIC fragments
    (no transcript overlap) the labels are reused as a pure naming
    convention — SENSE = read on (+) genomic strand, ANTISENSE = read
    on (−).  AMBIG marks fragments where transcripts on BOTH strands
    overlap (fragment is sense to one and antisense to another).
    """

    SENSE = 0
    ANTISENSE = 1
    AMBIG = 2


N_FRAGMENT_STRANDS = len(FragmentStrand)


@dataclass(frozen=True)
class CategorizedChunk:
    """Per-fragment SRD v2 categorization output for one buffer chunk.

    All arrays have length ``chunk.size`` (parallel to the input chunk).
    Fragments outside the unique-mapper-UNSPLICED-known-read-strand
    pre-filter receive ``category = 255`` (sentinel) and ``strand = 0``;
    callers should mask them out before tallying.

    Attributes
    ----------
    category : np.ndarray[uint8]
        :class:`FragmentCategory` value (or ``255`` sentinel) per fragment.
    strand : np.ndarray[uint8]
        :class:`FragmentStrand` value per fragment (``0`` for
        filtered-out — coincides with ``SENSE`` numerically; mask with
        ``keep`` before reading).
    keep : np.ndarray[bool]
        ``True`` iff the fragment passed the pre-filter; the
        strand/category fields are valid only where ``keep`` is set.
    frag_length : np.ndarray[int32]
        Fragment length for the calibration pool — sourced from
        ``genomic_footprint`` (always >= 0; SRD v2 Phase 1).
    """

    category: np.ndarray
    strand: np.ndarray
    keep: np.ndarray
    frag_length: np.ndarray

    @property
    def counts_by_category_strand(self) -> np.ndarray:
        """``(N_CATEGORIES, N_FRAGMENT_STRANDS)`` int64 count matrix.

        Counts only ``keep == True`` fragments.
        """
        out = np.zeros((N_CATEGORIES, N_FRAGMENT_STRANDS), dtype=np.int64)
        if not self.keep.any():
            return out
        cat = self.category[self.keep].astype(np.intp)
        st = self.strand[self.keep].astype(np.intp)
        flat = cat * N_FRAGMENT_STRANDS + st
        flat_counts = np.bincount(
            flat, minlength=N_CATEGORIES * N_FRAGMENT_STRANDS
        ).astype(np.int64)
        return flat_counts.reshape((N_CATEGORIES, N_FRAGMENT_STRANDS))


def _compute_fragment_strand(
    bp_pos: np.ndarray,
    bp_neg: np.ndarray,
    read_strand: np.ndarray,
    intergenic_mask: np.ndarray,
    *,
    read1_sense: bool = True,
) -> np.ndarray:
    """Return per-fragment :class:`FragmentStrand` value (uint8).

    Labels are in **biology / transcript frame**: SENSE means the
    fragment was transcribed from the overlapping transcript's
    template; ANTISENSE means it was transcribed from the opposite
    template (or, for gDNA, originated from the opposite genomic
    strand).

    The C++ scanner stores ``chunk.exon_strand`` as the read's
    *genomic* strand (after the dUTP-style R2 flip).  In an
    ``R1-sense`` library (``read1_sense=True``, e.g. KAPA Stranded),
    that genomic strand equals the transcript's template strand for
    sense fragments, so a direct strand comparison gives the
    biology-frame label.  In an ``R1-antisense`` library
    (``read1_sense=False``, e.g. Illumina TruSeq dUTP), the read's
    genomic strand is the **opposite** of the template strand for
    sense fragments — so we flip the read strand before comparing.

    Without this protocol-aware flip the labels are inverted on every
    R1-antisense library (the most common stranded protocol).  The
    discovery is documented in the SRD v3 oracle-BAM deep-dive
    (April 2026); see also ``docs/calibration/srd_v3_phase1_strand_relabel.md``
    §3 (Operation A) for the original framing.

    Parameters
    ----------
    bp_pos, bp_neg
        Discriminating per-strand bp counts. For exon classes use
        ``exon_bp_*``; for INTRONIC use the derived ``intron_bp_*``;
        for INTERGENIC the values are unused (forced to the
        SENSE/ANTISENSE-by-genomic-strand convention via
        ``intergenic_mask``).
    read_strand
        ``chunk.exon_strand`` — read's inferred genomic strand. The
        caller must have already filtered to ``{STRAND_POS, STRAND_NEG}``
        via the ``keep`` mask; behavior on AMBIG/UNKNOWN read-strand
        rows is unspecified.
    intergenic_mask
        Bool array marking INTERGENIC fragments.  For these the label
        is the SENSE/ANTISENSE convention based purely on
        the (protocol-corrected) read strand.
    read1_sense
        Library protocol from :class:`StrandModel`.  ``True`` for
        R1-sense protocols, ``False`` for R1-antisense.  Defaults to
        ``True`` for backwards compatibility with synthetic test
        fixtures that build R1-sense data.
    """
    out = np.zeros(bp_pos.shape, dtype=np.uint8)

    # Transcripts on both strands overlap → AMBIG (sense to one,
    # antisense to another). Wins over the read-strand comparison.
    both = (bp_pos > 0) & (bp_neg > 0)

    # Effective "template strand" of the fragment: in R1-sense
    # libraries this equals the read's genomic strand; in R1-antisense
    # libraries it is the opposite.
    if read1_sense:
        tpl_pos = read_strand == np.uint8(STRAND_POS)
        tpl_neg = read_strand == np.uint8(STRAND_NEG)
    else:
        tpl_pos = read_strand == np.uint8(STRAND_NEG)
        tpl_neg = read_strand == np.uint8(STRAND_POS)

    sense = (~both) & (
        (tpl_pos & (bp_pos > 0))
        | (tpl_neg & (bp_neg > 0))
    )
    antisense = (~both) & (
        (tpl_pos & (bp_neg > 0))
        | (tpl_neg & (bp_pos > 0))
    )

    out[sense] = int(FragmentStrand.SENSE)
    out[antisense] = int(FragmentStrand.ANTISENSE)
    out[both] = int(FragmentStrand.AMBIG)

    # INTERGENIC: pure naming convention — SENSE if the fragment's
    # template strand is (+), ANTISENSE if it is (−). Override
    # anything the above wrote (which is all-zero since
    # bp_pos == bp_neg == 0 by definition).
    out[intergenic_mask & tpl_pos] = int(FragmentStrand.SENSE)
    out[intergenic_mask & tpl_neg] = int(FragmentStrand.ANTISENSE)

    return out


def _per_fragment_max_exon_bp(
    chunk: _FinalizedChunk,
    t_strand_arr: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-fragment max ``exon_bp[T]`` over candidate transcripts ``T``,
    split by transcript strand.

    For each fragment ``i``, returns ``(max_pos[i], max_neg[i])`` where
    ``max_pos[i]`` = ``max{ chunk.exon_bp[k] : k ∈ [t_offsets[i], t_offsets[i+1]),
    t_strand_arr[chunk.t_indices[k]] == STRAND_POS }`` (and 0 when the set
    is empty).  Same for ``max_neg`` with ``STRAND_NEG``.

    The CSR ``exon_bp[k]`` column carries per-(fragment, candidate-transcript)
    annotated-exon coverage in basepairs (sum over read blocks).  For an
    UNSPLICED fragment whose mate pair lies inside a single exon of
    transcript ``T``, ``exon_bp[k_T] == read_length``.  This per-T
    maximum is therefore the precise check for "contained in an exon of
    *any* single transcript".

    Implementation uses ``numpy.maximum.at`` (unbuffered scatter-max)
    over a per-row group index built from ``t_offsets``.
    """
    n = chunk.size
    max_pos = np.zeros(n, dtype=np.int32)
    max_neg = np.zeros(n, dtype=np.int32)

    if chunk.exon_bp.size == 0 or n == 0:
        return max_pos, max_neg

    counts = np.diff(chunk.t_offsets).astype(np.int64, copy=False)
    if counts.sum() == 0:
        return max_pos, max_neg

    # Per-candidate transcript strand lookup.
    t_inds = chunk.t_indices
    valid = (t_inds >= 0) & (t_inds < t_strand_arr.size)
    cand_strand = np.zeros(t_inds.size, dtype=np.int32)
    cand_strand[valid] = t_strand_arr[t_inds[valid]].astype(np.int32, copy=False)

    bp = chunk.exon_bp.astype(np.int32, copy=False)
    bp_pos = np.where(cand_strand == np.int32(STRAND_POS), bp, np.int32(0))
    bp_neg = np.where(cand_strand == np.int32(STRAND_NEG), bp, np.int32(0))

    # Per-candidate → per-fragment row index.
    group_idx = np.repeat(np.arange(n, dtype=np.int64), counts)

    np.maximum.at(max_pos, group_idx, bp_pos)
    np.maximum.at(max_neg, group_idx, bp_neg)
    return max_pos, max_neg


def categorize_chunk(
    chunk: _FinalizedChunk,
    *,
    t_strand_arr: np.ndarray,
    exon_fit_tolerance_bp: int = 0,
    read1_sense: bool = True,
) -> CategorizedChunk:
    """Categorize every fragment in *chunk*.

    Pre-filter (``keep == True``):

    * ``num_hits == 1``                         — unique mapper
    * ``splice_type == UNSPLICED``              — truly unspliced
    * ``exon_strand ∈ {STRAND_POS, STRAND_NEG}`` — known read strand

    Fragments failing any of the above are flagged ``keep == False`` and
    should be ignored by downstream tallies.

    Parameters
    ----------
    chunk
        A finalised buffer chunk from ``FragmentBuffer.iter_chunks()``.
        Must carry the SRD v2 strand-bp columns (``exon_bp_pos``/``_neg``,
        ``tx_bp_pos``/``_neg``) and ``exon_strand``.
    t_strand_arr
        ``index.t_to_strand_arr`` — transcript strand lookup table
        (``STRAND_POS`` / ``STRAND_NEG`` / ``STRAND_NONE``) indexed by
        global ``t_index``.  Required for the per-transcript
        EXON_CONTAINED check (see ``_per_fragment_max_exon_bp``).
    exon_fit_tolerance_bp
        Tolerance (bp) for the EXON_CONTAINED test.  Defaults to 0:
        with the per-transcript ``max(exon_bp[T])`` formulation, oracle
        UNSPLICED fragments produced without alignment noise satisfy
        the check exactly, so any non-zero tolerance only obscures
        genuine alignment artefacts in real data.

    Returns
    -------
    CategorizedChunk
    """
    n = chunk.size

    # Pre-filter: unique-mapper, UNSPLICED, AND read strand known
    # (POS or NEG).  AMBIGUOUS / UNKNOWN read-strand fragments are
    # excluded entirely so that the transcript-frame strand label is
    # well-defined for every kept fragment.
    read_strand_known = (
        (chunk.exon_strand == np.uint8(STRAND_POS))
        | (chunk.exon_strand == np.uint8(STRAND_NEG))
    )
    keep = (
        (chunk.num_hits == 1)
        & (chunk.splice_type == np.uint8(SPLICE_UNSPLICED))
        & read_strand_known
    )

    cats = np.full(n, np.uint8(255), dtype=np.uint8)  # sentinel
    strand = np.zeros(n, dtype=np.uint8)
    frag_length = chunk.genomic_footprint.astype(np.int32, copy=False)

    if not keep.any():
        return CategorizedChunk(
            category=cats, strand=strand, keep=keep, frag_length=frag_length
        )

    rl = chunk.read_length.astype(np.int32, copy=False)
    ebp_pos = chunk.exon_bp_pos
    ebp_neg = chunk.exon_bp_neg
    tbp_pos = chunk.tx_bp_pos
    tbp_neg = chunk.tx_bp_neg
    read_strand = chunk.exon_strand

    # Derived per-strand intronic overlap (no ITYPE_INTRON in cgranges).
    ibp_pos = np.maximum(tbp_pos - ebp_pos, 0)
    ibp_neg = np.maximum(tbp_neg - ebp_neg, 0)

    # Per-transcript max exon coverage (the correct EXON_CONTAINED metric).
    # ``exon_bp[T]`` sums per-read-block coverage of fragment by ``T``'s
    # annotated exons; ``max_T exon_bp[T] >= read_length - tol`` iff
    # some single transcript ``T`` has its exon set fully cover every
    # block (and, since exons are contiguous, the inner mate gap too).
    # Aggregating ``exon_bp_pos/_neg`` is unsafe — those columns sum
    # per-cgranges-hit, double-counting overlapping exons of distinct
    # transcripts and (for paired-end fragments with non-zero gap)
    # collapsing to ``read_length`` regardless of containment.
    max_t_exon_bp_pos, max_t_exon_bp_neg = _per_fragment_max_exon_bp(
        chunk, t_strand_arr
    )
    max_t_exon_bp = np.maximum(max_t_exon_bp_pos, max_t_exon_bp_neg)

    tol = np.int32(int(exon_fit_tolerance_bp))

    # Mutually exclusive in priority order: INTERGENIC > INTRONIC >
    # EXON_CONTAINED > EXON_INCOMPATIBLE.
    intergenic_mask = (tbp_pos == 0) & (tbp_neg == 0)
    intronic_mask = (ebp_pos == 0) & (ebp_neg == 0) & ~intergenic_mask
    exon_contained_mask = (
        ~intergenic_mask
        & ~intronic_mask
        & ((rl - max_t_exon_bp) <= tol)
    )
    exon_incompat_mask = (
        ~intergenic_mask & ~intronic_mask & ~exon_contained_mask
    )

    cats[keep & intergenic_mask] = int(FragmentCategory.INTERGENIC)
    cats[keep & intronic_mask] = int(FragmentCategory.INTRONIC)
    cats[keep & exon_contained_mask] = int(FragmentCategory.EXON_CONTAINED)
    cats[keep & exon_incompat_mask] = int(FragmentCategory.EXON_INCOMPATIBLE)

    # Transcript-frame strand label per category — exon_bp_* for exon
    # classes, intron_bp_* for INTRONIC, SENSE/ANTISENSE-by-genomic-
    # strand convention for INTERGENIC (handled inside
    # _compute_fragment_strand via intergenic_mask).
    exon_strand_lbl = _compute_fragment_strand(
        ebp_pos, ebp_neg, read_strand,
        intergenic_mask=intergenic_mask, read1_sense=read1_sense,
    )
    intron_strand_lbl = _compute_fragment_strand(
        ibp_pos, ibp_neg, read_strand,
        intergenic_mask=intergenic_mask, read1_sense=read1_sense,
    )

    strand[keep & intergenic_mask] = exon_strand_lbl[keep & intergenic_mask]
    strand[keep & intronic_mask] = intron_strand_lbl[keep & intronic_mask]
    strand[keep & exon_contained_mask] = exon_strand_lbl[keep & exon_contained_mask]
    strand[keep & exon_incompat_mask] = exon_strand_lbl[keep & exon_incompat_mask]

    return CategorizedChunk(
        category=cats, strand=strand, keep=keep, frag_length=frag_length
    )
