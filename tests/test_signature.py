"""Tests for fine-region signature and 12-channel layout helpers."""

from __future__ import annotations

import pytest

from rigel.calibration import signature as sig
from rigel.calibration.regions import RegionStrand, RegionType


COARSE_TABLE = [
    (0x0, RegionType.INTERGENIC, RegionStrand.NONE, False),
    (0x1, RegionType.EXON, RegionStrand.NEG, False),
    (0x2, RegionType.EXON, RegionStrand.POS, False),
    (0x3, RegionType.EXON, RegionStrand.AMBIG, True),
    (0x4, RegionType.INTRON, RegionStrand.NEG, False),
    (0x5, RegionType.EXON, RegionStrand.NEG, False),
    (0x6, RegionType.EXON, RegionStrand.AMBIG, True),
    (0x7, RegionType.EXON, RegionStrand.AMBIG, True),
    (0x8, RegionType.INTRON, RegionStrand.POS, False),
    (0x9, RegionType.EXON, RegionStrand.AMBIG, True),
    (0xA, RegionType.EXON, RegionStrand.POS, False),
    (0xB, RegionType.EXON, RegionStrand.AMBIG, True),
    (0xC, RegionType.INTRON, RegionStrand.AMBIG, True),
    (0xD, RegionType.EXON, RegionStrand.AMBIG, True),
    (0xE, RegionType.EXON, RegionStrand.AMBIG, True),
    (0xF, RegionType.EXON, RegionStrand.AMBIG, True),
]


@pytest.mark.parametrize("signature, coarse_type, coarse_strand, ambiguous", COARSE_TABLE)
def test_coarse_derivation_table(signature, coarse_type, coarse_strand, ambiguous):
    assert sig.coarse_type_from_signature(signature) == int(coarse_type)
    assert sig.coarse_strand_from_signature(signature) == int(coarse_strand)
    assert sig.is_ambiguous_signature(signature) is ambiguous


@pytest.mark.parametrize("signature, _coarse_type, _coarse_strand, _ambiguous", COARSE_TABLE)
def test_pack_signature_round_trip(signature, _coarse_type, _coarse_strand, _ambiguous):
    packed = sig.pack_signature(
        intron_pos=bool(signature & sig.BIT_INTRON_POS),
        intron_neg=bool(signature & sig.BIT_INTRON_NEG),
        exon_pos=bool(signature & sig.BIT_EXON_POS),
        exon_neg=bool(signature & sig.BIT_EXON_NEG),
    )
    assert packed == signature


def test_channel_index_round_trip():
    seen = set()
    for compartment in range(sig.N_COMPARTMENTS):
        for splice_idx in range(sig.N_SPLICE_STATES):
            for strand_idx in range(sig.N_CHANNEL_STRANDS):
                channel = sig.channel_index(compartment, splice_idx, strand_idx)
                assert sig.channel_tuple(channel) == (compartment, splice_idx, strand_idx)
                seen.add(channel)
    assert seen == set(range(sig.N_CHANNELS))


def test_channel_layout_examples():
    assert sig.channel_index(sig.CHAN_CONTAINED, sig.SPLICE_UNSPLICED, sig.CHANNEL_STRAND_POS) == 0
    assert sig.channel_index(sig.CHAN_CONTAINED, sig.SPLICE_UNSPLICED, sig.CHANNEL_STRAND_NEG) == 1
    assert sig.channel_index(sig.CHAN_CONTAINED, sig.SPLICE_SPLICED, sig.CHANNEL_STRAND_POS) == 2
    assert (
        sig.channel_index(sig.CHAN_BOUNDARY_LEFT, sig.SPLICE_UNSPLICED, sig.CHANNEL_STRAND_POS) == 4
    )
    assert (
        sig.channel_index(sig.CHAN_BOUNDARY_RIGHT, sig.SPLICE_SPLICED, sig.CHANNEL_STRAND_NEG) == 11
    )


@pytest.mark.parametrize("bad_signature", [-1, 16, 99])
def test_invalid_signature_rejected(bad_signature):
    with pytest.raises(ValueError, match="signature"):
        sig.coarse_type_from_signature(bad_signature)


@pytest.mark.parametrize(
    "args",
    [
        (-1, 0, 0),
        (3, 0, 0),
        (0, -1, 0),
        (0, 2, 0),
        (0, 0, -1),
        (0, 0, 2),
    ],
)
def test_invalid_channel_args_rejected(args):
    with pytest.raises(ValueError):
        sig.channel_index(*args)
