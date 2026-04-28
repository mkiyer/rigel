"""SRD v2 Pass 0 — per-fragment geometric categorization.

Walks a ``_FinalizedChunk`` once and labels each unspliced fragment with one
of four geometric categories (see :class:`FragmentCategory`) and a
strand sub-label (see :class:`StrandLabel`).  Fully vectorised over the
chunk's strand-aware overlap columns produced by the C++ scanner:

* ``exon_bp_pos`` / ``exon_bp_neg``  — bp covered by ANY (+/−)-strand
  transcript's annotated exons.
* ``tx_bp_pos``  / ``tx_bp_neg``    — bp covered by ANY (+/−)-strand
  transcript's full span (transcript-tile interval).

Per-strand intronic overlap is *derived* in Python as
``intron_bp_strand = max(tx_bp_strand − exon_bp_strand, 0)`` because the
cgranges index has no explicit ``ITYPE_INTRON`` interval type.

The categorization is geometric and library-agnostic — strand
specificity affects which fragments populate which (category, strand)
cell, but the membership rule itself does not branch on SS.

Categories
----------
Per `docs/calibration/srd_v2_phase2plus_handoff.md` §2 / §4.

* ``INTERGENIC``      — ``tx_bp_pos == 0`` and ``tx_bp_neg == 0``
* ``INTRONIC``        — ``exon_bp_pos == 0`` and ``exon_bp_neg == 0``
                        and not INTERGENIC (i.e. inside a transcript span
                        but no exon overlap on either strand)
* ``EXON_CONTAINED``  — ``genomic_footprint − max(exon_bp_pos, exon_bp_neg)
                        <= exon_fit_tolerance_bp`` (the fragment fits
                        inside the exonic envelope of one strand within
                        a small tolerance)
* ``EXON_INCOMPATIBLE`` — none of the above (some exon overlap but
                        material genomic span lands outside the exonic
                        envelope of every strand)

Strand sub-label uses the strand-bp counts of the discriminating overlap
type (``exon_bp_*`` for exon classes, ``intron_bp_*`` derived for
INTRONIC, none for INTERGENIC).  Possible labels: ``NONE``, ``POS``,
``NEG``, ``AMBIG`` (both strands have overlap — strand-overlap zone).

Pre-filter
----------
Categorization is only meaningful for **unique-mapper, truly-unspliced**
fragments — ``num_hits == 1`` and ``splice_type == UNSPLICED``.  All
other splice types (annotated/unannotated/implicit/artifact) are
held out of calibration entirely; this filter is applied by the caller
in :mod:`rigel.calibration._simple`.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np

from ..buffer import _FinalizedChunk
from ..splice import SPLICE_UNSPLICED


class FragmentCategory(IntEnum):
    """SRD v2 per-fragment geometric category.  Values are stable;
    persisted as the row index into ``CalibrationResult.category_counts``.
    """

    INTERGENIC = 0
    INTRONIC = 1
    EXON_CONTAINED = 2
    EXON_INCOMPATIBLE = 3


N_CATEGORIES = len(FragmentCategory)


class StrandLabel(IntEnum):
    """Strand sub-label of a categorized fragment.  Values are stable;
    persisted as the column index into ``CalibrationResult.category_counts``.
    """

    NONE = 0
    POS = 1
    NEG = 2
    AMBIG = 3


N_STRAND_LABELS = len(StrandLabel)


@dataclass(frozen=True)
class CategorizedChunk:
    """Per-fragment SRD v2 categorization output for one buffer chunk.

    All arrays have length ``chunk.size`` (parallel to the input chunk).
    Fragments outside the ``UNSPLICED`` pre-filter receive
    ``category = 255`` (sentinel) and ``strand = 0``; callers should mask
    them out before tallying.

    Attributes
    ----------
    category : np.ndarray[uint8]
        :class:`FragmentCategory` value (or ``255`` sentinel) per fragment.
    strand : np.ndarray[uint8]
        :class:`StrandLabel` value per fragment (``0`` for filtered-out).
    keep : np.ndarray[bool]
        ``True`` iff the fragment passed the unspliced unique-mapper
        pre-filter; the strand/category fields are valid only where
        ``keep`` is set.
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
        """``(N_CATEGORIES, N_STRAND_LABELS)`` int64 count matrix.

        Counts only ``keep == True`` fragments.
        """
        out = np.zeros((N_CATEGORIES, N_STRAND_LABELS), dtype=np.int64)
        if not self.keep.any():
            return out
        cat = self.category[self.keep].astype(np.intp)
        st = self.strand[self.keep].astype(np.intp)
        # Flatten (cat, strand) into a single bincount index.
        flat = cat * N_STRAND_LABELS + st
        flat_counts = np.bincount(
            flat, minlength=N_CATEGORIES * N_STRAND_LABELS
        ).astype(np.int64)
        return flat_counts.reshape((N_CATEGORIES, N_STRAND_LABELS))


def _strand_label(bp_pos: np.ndarray, bp_neg: np.ndarray) -> np.ndarray:
    """Return per-fragment :class:`StrandLabel` value (uint8) from a pair
    of strand-bp arrays."""
    out = np.zeros(bp_pos.shape, dtype=np.uint8)
    p = bp_pos > 0
    n = bp_neg > 0
    out[p & ~n] = int(StrandLabel.POS)
    out[~p & n] = int(StrandLabel.NEG)
    out[p & n] = int(StrandLabel.AMBIG)
    return out


def categorize_chunk(
    chunk: _FinalizedChunk,
    *,
    exon_fit_tolerance_bp: int = 5,
) -> CategorizedChunk:
    """Categorize every fragment in *chunk*.

    Only fragments with ``num_hits == 1`` AND ``splice_type == UNSPLICED``
    receive a meaningful category; all others are flagged
    ``keep == False`` and should be ignored by downstream tallies.

    Parameters
    ----------
    chunk
        A finalised buffer chunk from ``FragmentBuffer.iter_chunks()``.
        Must carry the SRD v2 strand-bp columns (``exon_bp_pos``/``_neg``,
        ``tx_bp_pos``/``_neg``).
    exon_fit_tolerance_bp
        Tolerance (bp) for the EXON_CONTAINED test.  A fragment is
        contained iff
        ``genomic_footprint − max(exon_bp_pos, exon_bp_neg) <= tol``.

    Returns
    -------
    CategorizedChunk
    """
    n = chunk.size

    # Pre-filter: unique-mapper AND truly UNSPLICED.  All other splice
    # types (annotated, unannotated, implicit, artifact) are excluded
    # from the calibration pool entirely.
    keep = (chunk.num_hits == 1) & (chunk.splice_type == np.uint8(SPLICE_UNSPLICED))

    cats = np.full(n, np.uint8(255), dtype=np.uint8)  # sentinel
    strand = np.zeros(n, dtype=np.uint8)
    frag_length = chunk.genomic_footprint.astype(np.int32, copy=False)

    if not keep.any():
        return CategorizedChunk(
            category=cats, strand=strand, keep=keep, frag_length=frag_length
        )

    fp = chunk.genomic_footprint
    ebp_pos = chunk.exon_bp_pos
    ebp_neg = chunk.exon_bp_neg
    tbp_pos = chunk.tx_bp_pos
    tbp_neg = chunk.tx_bp_neg

    # Derived per-strand intronic overlap (no ITYPE_INTRON in cgranges).
    ibp_pos = np.maximum(tbp_pos - ebp_pos, 0)
    ibp_neg = np.maximum(tbp_neg - ebp_neg, 0)

    exon_either_lower = np.maximum(ebp_pos, ebp_neg).astype(np.int32, copy=False)
    tol = np.int32(int(exon_fit_tolerance_bp))

    # Mutually exclusive in priority order: INTERGENIC > INTRONIC >
    # EXON_CONTAINED > EXON_INCOMPATIBLE.  Strand-overlap zones (a
    # genomic position covered by both + and - transcripts) contribute
    # bp to BOTH ebp_pos and ebp_neg; those are surfaced via the AMBIG
    # strand sub-label.
    intergenic_mask = (tbp_pos == 0) & (tbp_neg == 0)
    intronic_mask = (ebp_pos == 0) & (ebp_neg == 0) & ~intergenic_mask
    exon_contained_mask = (
        ~intergenic_mask
        & ~intronic_mask
        & ((fp - exon_either_lower) <= tol)
    )
    exon_incompat_mask = (
        ~intergenic_mask & ~intronic_mask & ~exon_contained_mask
    )

    cats[keep & intergenic_mask] = int(FragmentCategory.INTERGENIC)
    cats[keep & intronic_mask] = int(FragmentCategory.INTRONIC)
    cats[keep & exon_contained_mask] = int(FragmentCategory.EXON_CONTAINED)
    cats[keep & exon_incompat_mask] = int(FragmentCategory.EXON_INCOMPATIBLE)

    # Strand sub-label per category — uses the discriminating overlap
    # type for each:  exon_bp_* for exon classes, intron_bp_* for
    # INTRONIC, NONE (zero by construction) for INTERGENIC.
    exon_strand_lbl = _strand_label(ebp_pos, ebp_neg)
    intron_strand_lbl = _strand_label(ibp_pos, ibp_neg)

    strand[keep & intronic_mask] = intron_strand_lbl[keep & intronic_mask]
    strand[keep & exon_contained_mask] = exon_strand_lbl[keep & exon_contained_mask]
    strand[keep & exon_incompat_mask] = exon_strand_lbl[keep & exon_incompat_mask]
    # INTERGENIC keeps strand = NONE (0) by construction.

    return CategorizedChunk(
        category=cats, strand=strand, keep=keep, frag_length=frag_length
    )
