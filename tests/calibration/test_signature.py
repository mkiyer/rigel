"""The region signature — packing, deriving and the strand classes read off it.

A signature is four bits (exon/intron × +/−) and every downstream strand question is a function of
them, so these gates hold the encoding itself: the bit layout, the range check, the coarse strand a
signature collapses to, the vectorised twin agreeing with the scalar form, the one shared numbering
with ``rigel.types.Strand``, and the active-strand classifier a boundary's taxonomy is built from.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    TS_AMBIG,
    TS_NEG,
    TS_NONE,
    TS_POS,
    RegionStrand,
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
# the RNA-active classifier
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


def test_boundary_taxonomy_from_flank_helpers():
    """The boundary types as the AND of the two flanks' helper masks: a strand crosses iff both flanks
    are nrna-active."""
    # (left_sig, right_sig) per type, tested on the + strand unless noted.
    intergenic, exon_p, intron_p = 0, BIT_EXON_POS, BIT_INTRON_POS
    ambig = BIT_EXON_POS | BIT_EXON_NEG

    def crosses(left, right):
        nlp, _ = nrna_active_strands(np.array([left]))
        nrp, _ = nrna_active_strands(np.array([right]))
        return bool((nlp & nrp)[0])

    # 1) intergenic ↔ exon: no + transcript on the left flank ⇒ no crossing ⇒ gDNA sink.
    assert crosses(intergenic, exon_p) is False
    # 2) intron ↔ exon, 3) exon ↔ exon, 4) ambig ↔ ambig: a + transcript on both flanks ⇒ it crosses.
    assert crosses(intron_p, exon_p) is True
    assert crosses(exon_p, exon_p) is True
    assert crosses(ambig, ambig) is True
    # ambig ↔ ambig also crosses on the − strand (the AMBIG 2-D region).
    _, nln = nrna_active_strands(np.array([ambig]))
    _, nrn = nrna_active_strands(np.array([ambig]))
    assert bool((nln & nrn)[0])
