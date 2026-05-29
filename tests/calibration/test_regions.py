"""Region partition builder — merged constant-signature segments."""

from __future__ import annotations

from rigel.calibration.regions import build_region_partition, validate_against_ref_lengths
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_POS,
)
from rigel.transcript import Transcript
from rigel.types import Interval, Strand

REF_LEN = 1000
REF_LENGTHS = {"chr1": REF_LEN}


def _tx(strand, exons, t_id):
    return Transcript(
        ref="chr1",
        strand=strand,
        exons=[Interval(s, e) for s, e in exons],
        t_id=t_id,
    )


def test_single_transcript_partition():
    # One + transcript: exons [100,200),[400,500), intron [200,400).
    tx = _tx(Strand.POS, [(100, 200), (400, 500)], "t1")
    df = build_region_partition([tx], REF_LENGTHS)

    rows = list(zip(df["start"], df["end"], df["signature"]))
    assert rows == [
        (0, 100, 0),
        (100, 200, BIT_EXON_POS),
        (200, 400, BIT_INTRON_POS),
        (400, 500, BIT_EXON_POS),
        (500, REF_LEN, 0),
    ]
    # region_id is sequential and the partition validates.
    assert list(df["region_id"]) == list(range(len(df)))
    validate_against_ref_lengths(df.set_index("region_id", drop=False), REF_LENGTHS)


def test_overlapping_same_strand_exons_merge():
    # tx1 exons [100,200),[400,500); tx2 single exon [150,250) overlaps tx1.
    tx1 = _tx(Strand.POS, [(100, 200), (400, 500)], "t1")
    tx2 = _tx(Strand.POS, [(150, 250)], "t2")
    df = build_region_partition([tx1, tx2], REF_LENGTHS)

    rows = list(zip(df["start"], df["end"], df["signature"]))
    # The inner breakpoint at 150 does NOT split [100,200): both sides are
    # exon_pos, so the merge collapses them into one region.
    assert rows == [
        (0, 100, 0),
        (100, 200, BIT_EXON_POS),
        (200, 250, BIT_EXON_POS | BIT_INTRON_POS),
        (250, 400, BIT_INTRON_POS),
        (400, 500, BIT_EXON_POS),
        (500, REF_LEN, 0),
    ]


def test_neighbour_differs_invariant():
    # Opposite-strand exons that abut must NOT merge (signatures differ).
    tx_pos = _tx(Strand.POS, [(100, 300)], "tp")
    tx_neg = _tx(Strand.NEG, [(300, 500)], "tn")
    df = build_region_partition([tx_pos, tx_neg], REF_LENGTHS)

    sigs = df["signature"].to_numpy()
    # No two adjacent regions share a signature.
    assert (sigs[:-1] != sigs[1:]).all()
    assert BIT_EXON_POS in set(sigs)
    assert BIT_EXON_NEG in set(sigs)
    validate_against_ref_lengths(df.set_index("region_id", drop=False), REF_LENGTHS)


def test_no_transcripts_single_intergenic_region():
    df = build_region_partition([], REF_LENGTHS)
    rows = list(zip(df["start"], df["end"], df["signature"]))
    assert rows == [(0, REF_LEN, 0)]
