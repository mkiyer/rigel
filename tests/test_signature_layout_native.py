"""Python-to-native layout parity for fine-region signature helpers."""

from __future__ import annotations

from rigel import native
from rigel.calibration import signature as sig


def test_signature_constants_match_native():
    assert native.REGION_SIG_INTRON_POS == sig.BIT_INTRON_POS
    assert native.REGION_SIG_INTRON_NEG == sig.BIT_INTRON_NEG
    assert native.REGION_SIG_EXON_POS == sig.BIT_EXON_POS
    assert native.REGION_SIG_EXON_NEG == sig.BIT_EXON_NEG
    assert native.REGION_SIG_N_STATES == sig.N_SIGNATURES


def test_channel_constants_match_native():
    assert native.REGION_CHAN_CONTAINED == sig.CHAN_CONTAINED
    assert native.REGION_CHAN_BOUNDARY_LEFT == sig.CHAN_BOUNDARY_LEFT
    assert native.REGION_CHAN_BOUNDARY_RIGHT == sig.CHAN_BOUNDARY_RIGHT
    assert native.REGION_N_CHANNELS == sig.N_CHANNELS


def test_pack_and_coarse_derivation_match_native():
    for signature in range(sig.N_SIGNATURES):
        assert (
            native.region_pack_signature(
                bool(signature & sig.BIT_INTRON_POS),
                bool(signature & sig.BIT_INTRON_NEG),
                bool(signature & sig.BIT_EXON_POS),
                bool(signature & sig.BIT_EXON_NEG),
            )
            == signature
        )
        assert native.region_coarse_type_from_signature(
            signature
        ) == sig.coarse_type_from_signature(signature)
        assert native.region_coarse_strand_from_signature(
            signature
        ) == sig.coarse_strand_from_signature(signature)


def test_channel_index_matches_native():
    for compartment in range(sig.N_COMPARTMENTS):
        for splice_idx in range(sig.N_SPLICE_STATES):
            for strand_idx in range(sig.N_CHANNEL_STRANDS):
                assert native.region_channel_index(
                    compartment, splice_idx, strand_idx
                ) == sig.channel_index(compartment, splice_idx, strand_idx)
