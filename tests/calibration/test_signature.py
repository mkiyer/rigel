"""Region signature encoding — pack/derive/strand-class."""

from __future__ import annotations

import numpy as np
import pytest

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
    RegionType,
    coarse_strand_from_signature,
    coarse_type_from_signature,
    is_ambiguous_signature,
    pack_signature,
    transcript_strand_class,
    validate_signature,
)


def test_pack_signature_bits():
    assert pack_signature() == 0
    assert pack_signature(exon_pos=True) == BIT_EXON_POS
    assert pack_signature(intron_neg=True) == BIT_INTRON_NEG
    # Overlapping +/- exons: both exon bits set.
    assert pack_signature(exon_pos=True, exon_neg=True) == (BIT_EXON_POS | BIT_EXON_NEG)
    # Full 4-bit signature.
    assert pack_signature(intron_pos=True, intron_neg=True, exon_pos=True, exon_neg=True) == (
        BIT_INTRON_POS | BIT_INTRON_NEG | BIT_EXON_POS | BIT_EXON_NEG
    )


def test_validate_signature_range():
    assert validate_signature(0) == 0
    assert validate_signature(15) == 15
    with pytest.raises(ValueError):
        validate_signature(16)
    with pytest.raises(ValueError):
        validate_signature(-1)


def test_coarse_type_exon_wins_over_intron():
    assert coarse_type_from_signature(0) == int(RegionType.INTERGENIC)
    assert coarse_type_from_signature(BIT_INTRON_POS) == int(RegionType.INTRON)
    assert coarse_type_from_signature(BIT_EXON_NEG) == int(RegionType.EXON)
    # Mixed exon + intron → EXON (exon wins).
    assert coarse_type_from_signature(BIT_EXON_POS | BIT_INTRON_POS) == int(RegionType.EXON)


def test_coarse_strand_from_signature():
    assert coarse_strand_from_signature(0) == int(RegionStrand.NONE)
    assert coarse_strand_from_signature(BIT_EXON_POS) == int(RegionStrand.POS)
    assert coarse_strand_from_signature(BIT_INTRON_NEG) == int(RegionStrand.NEG)
    # Both strands present → AMBIG.
    assert coarse_strand_from_signature(BIT_EXON_POS | BIT_INTRON_NEG) == int(RegionStrand.AMBIG)
    assert is_ambiguous_signature(BIT_EXON_POS | BIT_EXON_NEG)
    assert not is_ambiguous_signature(BIT_EXON_POS)


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
    # transcript_strand_class is the vectorised twin of coarse_strand_from_signature.
    np.testing.assert_array_equal(ts, [coarse_strand_from_signature(int(s)) for s in sig])


def test_strand_convention_unified():
    """TS_* is ONE convention with RegionStrand (== rigel.types.Strand): NONE=0,
    POS=1, NEG=2, AMBIG=3. (Regression: TS_NEG was historically -1, TS_AMBIG 2.)"""
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
