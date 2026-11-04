"""The region signature: the per-strand RNA masks read off it.

A signature is four bits (exon/intron × +/−), and every downstream strand question is a function of them, so these tests
hold the active-strand classifier a region's and a boundary's ``free_*`` gate is built from.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    rna_active_strands,
)


# ---------------------------------------------------------------------------
# the RNA-active classifier
# ---------------------------------------------------------------------------


def test_rna_active_strands():
    """RNA-active means a transcript is present (exon or intron) on the strand."""
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
    pos, neg = rna_active_strands(sig)
    np.testing.assert_array_equal(pos, [False, True, True, False, False, True, True])
    np.testing.assert_array_equal(neg, [False, False, False, True, True, True, True])


def test_boundary_taxonomy_from_flank_helpers():
    """The boundary types as the AND of the two flanks' helper masks: a strand crosses iff both flanks
    are RNA-active."""
    # (left_sig, right_sig) per type, tested on the + strand unless noted.
    intergenic, exon_p, intron_p = 0, BIT_EXON_POS, BIT_INTRON_POS
    ambig = BIT_EXON_POS | BIT_EXON_NEG

    def crosses(left, right):
        left_pos, _ = rna_active_strands(np.array([left]))
        right_pos, _ = rna_active_strands(np.array([right]))
        return bool((left_pos & right_pos)[0])

    # 1) intergenic ↔ exon: no + transcript on the left flank, so no crossing: a gDNA sink.
    assert crosses(intergenic, exon_p) is False
    # 2) intron ↔ exon, 3) exon ↔ exon, 4) ambig ↔ ambig: a + transcript on both flanks, so it crosses.
    assert crosses(intron_p, exon_p) is True
    assert crosses(exon_p, exon_p) is True
    assert crosses(ambig, ambig) is True
    # ambig ↔ ambig also crosses on the − strand (the AMBIG 2-D region).
    _, left_neg = rna_active_strands(np.array([ambig]))
    _, right_neg = rna_active_strands(np.array([ambig]))
    assert bool((left_neg & right_neg)[0])
