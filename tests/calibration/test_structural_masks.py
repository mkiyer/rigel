"""Structural eligibility masks for the three-channel (mature/nascent/gDNA) calibration design.

These assert the index-build structural properties (`docs/calibration/three_component_mature_nascent_design.md`
§4) on hand-built transcript topologies: per-region `mature_eligible_{pos,neg}` (multi-exon exon coverage) and
the per-boundary annotation flags (`is_tss`/`is_tes`/`is_splice_junction`/`genomic_sj_strand`). The masks are the
foundation of the mature/nascent message split, so every hard case — single-exon, intron-retention, opposite-
strand overlap (AMBIG), and the multi-transcript "mature can splice or continue" seam — is pinned here.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.regions import (
    build_boundary_partition,
    build_region_partition,
    validate_boundaries_against_regions,
)
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _tx(t_id: str, strand: Strand, exons: list[tuple[int, int]]) -> Transcript:
    return Transcript(ref="chr1", strand=strand, exons=[Interval(s, e) for s, e in exons], t_id=t_id)


def _region_map(region_df, col: str) -> dict[tuple[int, int], bool]:
    return {
        (int(s), int(e)): bool(v)
        for s, e, v in zip(region_df["start"], region_df["end"], region_df[col])
    }


def _boundary_map(boundary_df, col):
    return {int(p): v for p, v in zip(boundary_df["position"], boundary_df[col])}


def test_three_transcript_seam_mature_and_junctions():
    """TA/TB both multi-exon; the pure intron (9000,16000) is the sole nascent source; every other genic
    region is mature-eligible. All four junctions are + splice junctions; the outer termini are TSS/TES."""
    ta = _tx("TA", Strand.POS, [(1000, 5000), (20000, 24000)])
    tb = _tx("TB", Strand.POS, [(1000, 9000), (16000, 24000)])
    ref_lengths = {"chr1": 30000}
    rdf = build_region_partition([ta, tb], ref_lengths)
    bdf = build_boundary_partition(rdf, [ta, tb], ref_lengths)
    validate_boundaries_against_regions(bdf, rdf)

    me = _region_map(rdf, "mature_eligible_pos")
    assert me == {
        (0, 1000): False,       # intergenic
        (1000, 5000): True,     # TA/TB exon
        (5000, 9000): True,     # TA intron + TB exon → mature via TB (multi-exon)
        (9000, 16000): False,   # pure intron → nascent source
        (16000, 20000): True,   # TA intron + TB exon
        (20000, 24000): True,   # TA/TB exon
        (24000, 30000): False,  # intergenic
    }
    assert not rdf["mature_eligible_neg"].any()  # no − transcripts

    sj = _boundary_map(bdf, "is_splice_junction")
    assert [p for p, v in sj.items() if v] == [5000, 9000, 16000, 20000]
    strand = _boundary_map(bdf, "genomic_sj_strand")
    assert all(strand[p] == int(Strand.POS) for p in (5000, 9000, 16000, 20000))
    assert [p for p, v in _boundary_map(bdf, "is_tss").items() if v] == [1000]
    assert [p for p, v in _boundary_map(bdf, "is_tes").items() if v] == [24000]


def test_intron_retention_region_is_not_mature_eligible():
    """A retained-intron region (TA's intron + a single-exon TB spanning it) is NOT mature-eligible — its RNA
    is nascent (TB single-exon ≡ nascent + TA's nascent), even though it carries an EX+ signature bit."""
    ta = _tx("TA", Strand.POS, [(1000, 5000), (10000, 14000)])
    tb = _tx("TB", Strand.POS, [(1000, 14000)])  # single-exon, retains TA's intron
    rdf = build_region_partition([ta, tb], {"chr1": 20000})
    me = _region_map(rdf, "mature_eligible_pos")
    assert me[(5000, 10000)] is False  # TA intron + TB single-exon → nascent
    assert me[(1000, 5000)] is True    # TA/TB exon, TA is multi-exon
    assert me[(10000, 14000)] is True


def test_single_exon_transcript_never_mature_eligible():
    """A single-exon transcript has no splicing (mature ≡ nascent) → no region is mature-eligible and no
    boundary is a splice junction; only its two termini are TSS/TES (orientation by strand)."""
    tc = _tx("TC", Strand.NEG, [(2000, 6000)])
    ref_lengths = {"chr1": 10000}
    rdf = build_region_partition([tc], ref_lengths)
    bdf = build_boundary_partition(rdf, [tc], ref_lengths)
    assert not rdf["mature_eligible_pos"].any()
    assert not rdf["mature_eligible_neg"].any()
    assert not bdf["is_splice_junction"].any()
    # − strand: TSS is the genomically-higher terminus (last exon end), TES the lower (first exon start).
    assert [p for p, v in _boundary_map(bdf, "is_tss").items() if v] == [6000]
    assert [p for p, v in _boundary_map(bdf, "is_tes").items() if v] == [2000]


def test_opposite_strand_overlap_ambig_both_strands_mature_eligible():
    """Overlapping multi-exon transcripts on opposite strands: mature-eligibility is independent per strand."""
    td = _tx("TD", Strand.POS, [(1000, 3000), (5000, 7000)])
    te = _tx("TE", Strand.NEG, [(2000, 4000), (6000, 8000)])
    rdf = build_region_partition([td, te], {"chr1": 10000})
    assert rdf["mature_eligible_pos"].any()
    assert rdf["mature_eligible_neg"].any()
    # A region covered only by TD's exon is mature-eligible on + but not −, and vice-versa.
    mp = _region_map(rdf, "mature_eligible_pos")
    mn = _region_map(rdf, "mature_eligible_neg")
    assert mp[(1000, 2000)] is True and mn[(1000, 2000)] is False  # TD exon only


def test_boundary_count_matches_regions_plus_one():
    ta = _tx("TA", Strand.POS, [(1000, 5000), (20000, 24000)])
    tb = _tx("TB", Strand.POS, [(1000, 9000), (16000, 24000)])
    ref_lengths = {"chr1": 30000}
    rdf = build_region_partition([ta, tb], ref_lengths)
    bdf = build_boundary_partition(rdf, [ta, tb], ref_lengths)
    assert len(bdf) == len(rdf) + 1
    # positions are exactly the region interfaces [start…, last end]
    starts = rdf["start"].to_numpy(np.int64)
    ends = rdf["end"].to_numpy(np.int64)
    expected = np.append(starts, ends[-1])
    assert np.array_equal(bdf["position"].to_numpy(np.int64), expected)


def test_coincident_opposite_strand_junction_marked_both():
    """Two transcripts on opposite strands whose introns share an endpoint → that boundary is a splice
    junction on both strands (genomic_sj_strand == AMBIGUOUS/3), the documented rare edge case."""
    # + intron 5000-9000 (donor at 5000); − intron 3000-5000 (its endpoint 5000 is an acceptor on −).
    tp = _tx("TP", Strand.POS, [(1000, 5000), (9000, 12000)])
    tm = _tx("TM", Strand.NEG, [(2000, 3000), (5000, 8000)])
    rdf = build_region_partition([tp, tm], {"chr1": 15000})
    bdf = build_boundary_partition(rdf, [tp, tm], {"chr1": 15000})
    strand = _boundary_map(bdf, "genomic_sj_strand")
    assert strand[5000] == int(Strand.AMBIGUOUS)  # + donor and − endpoint coincide
