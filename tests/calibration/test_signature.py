"""The region signature — packing, deriving and the strand classes read off it.

A signature is four bits (exon/intron × +/−) and every downstream strand question is a function of
them, so these gates hold the encoding itself: the bit layout, the range check, the coarse strand a
signature collapses to, the vectorised twin agreeing with the scalar form, the one shared numbering
with ``rigel.types.Strand``, and the nascent/mature active-strand classifiers a boundary's
taxonomy is then built from.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    N_SIGNATURES,
    TS_AMBIG,
    TS_NEG,
    TS_NONE,
    TS_POS,
    RegionStrand,
    mrna_active_strands,
    nrna_active_strands,
    transcript_strand_class,
)


def test_transcript_strand_class_array():
    sig = np.array(
        [
            0,  # intergenic → NONE
            BIT_EXON_POS,  # pos → POS
            BIT_INTRON_POS,  # pos → POS
            BIT_EXON_NEG,  # neg → NEG
            BIT_EXON_POS | BIT_EXON_NEG,  # overlap → AMBIG
            BIT_INTRON_POS | BIT_EXON_NEG,  # both strands → AMBIG
        ],
        dtype=np.uint8,
    )
    ts = transcript_strand_class(sig)
    assert ts.dtype == np.int8
    np.testing.assert_array_equal(ts, [TS_NONE, TS_POS, TS_POS, TS_NEG, TS_AMBIG, TS_AMBIG])


def test_strand_convention_unified():
    """``TS_*`` is one convention with ``RegionStrand`` (== ``rigel.types.Strand``): NONE=0, POS=1,
    NEG=2, AMBIG=3. Two numberings for one concept is how a sign error becomes invisible."""
    from rigel.types import Strand

    assert (TS_NONE, TS_POS, TS_NEG, TS_AMBIG) == (0, 1, 2, 3)
    assert (
        int(RegionStrand.NONE),
        int(RegionStrand.POS),
        int(RegionStrand.NEG),
        int(RegionStrand.AMBIG),
    ) == (
        TS_NONE,
        TS_POS,
        TS_NEG,
        TS_AMBIG,
    )
    assert (int(Strand.POS), int(Strand.NEG)) == (TS_POS, TS_NEG)


# ---------------------------------------------------------------------------
# nRNA/mRNA-active classifier
# ---------------------------------------------------------------------------


def test_nrna_active_strands():
    """Nascent-active = a transcript is present (exon OR intron) on the strand."""
    sig = np.array(
        [
            0,  # intergenic
            BIT_EXON_POS,  # exon+
            BIT_INTRON_POS,  # intron+
            BIT_EXON_NEG,  # exon−
            BIT_INTRON_NEG,  # intron−
            BIT_EXON_POS | BIT_EXON_NEG,  # ambig exon
            BIT_INTRON_POS | BIT_EXON_NEG,  # +intron, −exon
        ],
        dtype=np.int64,
    )
    pos, neg = nrna_active_strands(sig)
    np.testing.assert_array_equal(pos, [False, True, True, False, False, True, True])
    np.testing.assert_array_equal(neg, [False, False, False, True, True, True, True])


def test_mrna_active_strands():
    """Mature-active = an EXON is present on the strand (introns carry no mature RNA)."""
    sig = np.array(
        [
            0,  # intergenic
            BIT_EXON_POS,  # exon+
            BIT_INTRON_POS,  # intron+ → NO mature
            BIT_EXON_NEG,  # exon−
            BIT_INTRON_NEG,  # intron− → NO mature
            BIT_EXON_POS | BIT_EXON_NEG,  # ambig exon
            BIT_INTRON_POS | BIT_EXON_NEG,  # +intron (no mature), −exon (mature)
        ],
        dtype=np.int64,
    )
    pos, neg = mrna_active_strands(sig)
    np.testing.assert_array_equal(pos, [False, True, False, False, False, True, False])
    np.testing.assert_array_equal(neg, [False, False, False, True, False, True, True])


def test_mrna_implies_nrna_all_signatures():
    """Across all 16 signatures, mature-active ⇒ nascent-active on each strand (an exon carries both;
    an intron carries only nascent)."""
    sig = np.arange(N_SIGNATURES, dtype=np.int64)
    nrp, nrn = nrna_active_strands(sig)
    mrp, mrn = mrna_active_strands(sig)
    assert np.all(~mrp | nrp)  # mrna_pos ⇒ nrna_pos
    assert np.all(~mrn | nrn)


def test_boundary_taxonomy_from_flank_helpers():
    """The four boundary types as the AND of the two flanks' helper masks.
    A strand crosses (nascent) iff both flanks are nrna-active; it is mature-capable iff both are
    mrna-active."""
    # (left_sig, right_sig) per type, tested on the + strand unless noted.
    intergenic, exon_p, intron_p = 0, BIT_EXON_POS, BIT_INTRON_POS
    ambig = BIT_EXON_POS | BIT_EXON_NEG

    def cross_mature(left, right):
        L, R = np.array([left]), np.array([right])
        nlp, _ = nrna_active_strands(L)
        nrp, _ = nrna_active_strands(R)
        mlp, mln = mrna_active_strands(L)
        mrp, mrn = mrna_active_strands(R)
        return bool((nlp & nrp)[0]), bool((mlp & mrp)[0])

    # 1) intergenic ↔ exon: no + transcript on the left flank ⇒ no crossing ⇒ gDNA sink.
    assert cross_mature(intergenic, exon_p) == (False, False)
    # 2) intron ↔ exon: + transcript on both, but the intron has no exon ⇒ nascent-only.
    assert cross_mature(intron_p, exon_p) == (True, False)
    # 3) exon ↔ exon: contiguous + exon ⇒ mature-capable.
    assert cross_mature(exon_p, exon_p) == (True, True)
    # 4) ambig ↔ ambig: both strands cross; on + it is mature-capable.
    assert cross_mature(ambig, ambig) == (True, True)
    # ambig ↔ ambig also crosses on the − strand (the AMBIG 2-D region).
    _, nln = nrna_active_strands(np.array([ambig]))
    _, nrn = nrna_active_strands(np.array([ambig]))
    assert bool((nln & nrn)[0])
