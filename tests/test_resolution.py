"""Fragment resolution: turning an aligned pair into a fragment with a candidate transcript set.

Splice-junction narrowing, the per-candidate exonic bases, candidate filtering and isoform discrimination, all
through the C++ ``FragmentResolver`` on indexes built from GTF strings.
"""

import textwrap

import pytest
from _index_builder import build_test_index

from rigel.types import Strand
from rigel.splice import SpliceType
from _resolution_reference import (
    GenomicInterval,
    make_fragment,
    resolve_fragment,
)


# =====================================================================
# resolve_fragment through the C++ FragmentResolver
# =====================================================================
#
# Each test class uses a session-scoped fixture that builds the
# TranscriptIndex it needs; the shared mini_index fixture (from
# conftest.py) serves when the MINI_GTF layout suffices.
# =====================================================================


def _t_map(index):
    """Return dict mapping transcript ID string to integer index."""
    return dict(zip(index.t_df["t_id"], index.t_df["t_index"]))


def _exon(ref, start, end, strand=Strand.POS):
    return GenomicInterval(ref, start, end, strand)


def _ambiguous_multiexon_gtf(n_transcripts=8):
    lines = []
    for t_idx in range(n_transcripts):
        exon_coords = [
            (100, 150),
            (200, 250),
            (300, 350),
            (400, 450),
            (600 + t_idx, 650 + t_idx),
        ]
        for start, end in exon_coords:
            lines.append(
                f"chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\t"
                f'gene_id "g1"; transcript_id "t{t_idx}"; '
                f'gene_name "G1"; gene_type "protein_coding"; tag "basic";'
            )
    return "\n".join(lines) + "\n"


# GTF for overlap tests: 3 transcripts in one gene
# 0-based half-open after parse:
#   t_two_exon: exons (100,200), (350,400)   transcript span (100,400)
#   t_one_left: exon  (100,200)              transcript span (100,200)
#   t_one_right: exon (350,500)              transcript span (350,500)
OVERLAP_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t_two_exon"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t351\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t_two_exon"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t_one_left"; gene_name "G1"; gene_type "protein_coding";
    chr1\ttest\texon\t351\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t_one_right"; gene_name "G1"; gene_type "protein_coding";
""")


@pytest.fixture(scope="session")
def overlap_index(tmp_path_factory):
    """3-transcript overlap-test index."""
    return build_test_index(tmp_path_factory, OVERLAP_GTF, name="overlap_idx")


# Two single-exon genes on opposite strands, overlapping on [300, 500) after parse.
ANTISENSE_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t101\t500\t.\t+\t.\tgene_id "gp"; transcript_id "tp"; gene_name "GP"; gene_type "protein_coding";
    chr1\ttest\texon\t301\t700\t.\t-\t.\tgene_id "gm"; transcript_id "tm"; gene_name "GM"; gene_type "protein_coding";
""")


@pytest.fixture(scope="session")
def antisense_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, ANTISENSE_GTF, name="antisense_idx")


class TestTheIndexInstallsTheTranscriptStrands:
    """The resolver reads each transcript's strand, and a fragment whose candidates disagree on it is
    ``ambig_strand``. The strands are an index fact, so the loaded index installs them and a fragment resolved
    without a scan sees them too."""

    def test_an_antisense_overlap_is_strand_ambiguous_without_a_scan(self, antisense_index):
        tm = _t_map(antisense_index)
        result = resolve_fragment(make_fragment(exons=(_exon("chr1", 350, 450),)), antisense_index)
        assert result is not None
        assert {tm["tp"], tm["tm"]} <= result.t_inds
        assert result.ambig_strand == 1

    def test_a_same_strand_overlap_is_not(self, overlap_index):
        tm = _t_map(overlap_index)
        result = resolve_fragment(make_fragment(exons=(_exon("chr1", 150, 190),)), overlap_index)
        assert result is not None
        assert {tm["t_two_exon"], tm["t_one_left"]} <= result.t_inds
        assert result.ambig_strand == 0


class TestResolveFragment:
    """Annotated SJ narrows transcript candidates via C++ FragmentResolver.

    Uses MINI_GTF:
      t1(idx=0): exons (99,200),(299,400),(499,600) → SJs (200,299), (400,499)
      t2(idx=1): exons (99,200),(499,600)           → SJ  (200,499)
    """

    def test_annotated_sj_prioritizes_sj_supported_transcripts(self, mini_index):
        """Fragment spanning SJ (200,299) resolves to t1 only."""
        tm = _t_map(mini_index)
        frag = make_fragment(
            exons=(
                _exon("chr1", 150, 200),
                _exon("chr1", 299, 350),
            ),
            introns=(GenomicInterval("chr1", 200, 299, Strand.POS),),
        )
        result = resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_ANNOT)
        assert result.t_inds == frozenset({tm["t1"]})

    def test_annotated_sj_matches_when_intron_strand_unknown(self, mini_index):
        """SJ match still works when fragment's intron strand is NONE."""
        tm = _t_map(mini_index)
        frag = make_fragment(
            exons=(
                _exon("chr1", 150, 200),
                _exon("chr1", 299, 350),
            ),
            introns=(GenomicInterval("chr1", 200, 299, Strand.NONE),),
        )
        result = resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_ANNOT)
        assert result.t_inds == frozenset({tm["t1"]})


class TestOverlapProfileViaResolve:
    """Test the per-candidate exonic bases produced by resolve_fragment.

    Uses overlap_index (OVERLAP_GTF):
      t_two_exon:  exons (100,200),(350,400)  — 2-exon, intron (200,350)
      t_one_left:  exon  (100,200)
      t_one_right: exon  (350,500)

    The 2-exon transcript's intron (200,350) is unambiguous: no other
    transcript's exon covers that region.
    """

    def test_full_exon_overlap_single_block(self, overlap_index):
        """Fragment fully inside shared exon → exon_bp=read_length."""
        tm = _t_map(overlap_index)
        frag = make_fragment(exons=(_exon("chr1", 150, 200),), introns=())
        result = resolve_fragment(frag, overlap_index)
        assert result is not None
        assert result.read_length == 50
        assert result.exon_bp[tm["t_two_exon"]] == 50
        assert result.exon_bp[tm["t_one_left"]] == 50

    def test_intronic_overlap(self, overlap_index):
        """Fragment in intronic region → exon_bp=0."""
        tm = _t_map(overlap_index)
        frag = make_fragment(exons=(_exon("chr1", 220, 320),), introns=())
        result = resolve_fragment(frag, overlap_index)
        assert result is not None
        assert result.read_length == 100
        # Only t_two_exon's transcript span covers 220-320
        t2e = tm["t_two_exon"]
        assert t2e in result.t_inds
        assert result.exon_bp[t2e] == 0

    def test_mixed_exon_intron_overlap(self, overlap_index):
        """Fragment spanning exon-intron boundary → partial exon bp."""
        tm = _t_map(overlap_index)
        frag = make_fragment(exons=(_exon("chr1", 180, 220),), introns=())
        result = resolve_fragment(frag, overlap_index)
        assert result is not None
        assert result.read_length == 40
        # t_two_exon: exon(100-200) covers 180-200=20bp
        assert result.exon_bp[tm["t_two_exon"]] == 20
        # t_one_left: exon(100-200) covers 180-200=20bp
        assert result.exon_bp[tm["t_one_left"]] == 20

    def test_multi_block_overlap_sums(self, overlap_index):
        """Two exon blocks accumulate exon overlap: 50+50=100bp for t_two_exon."""
        tm = _t_map(overlap_index)
        frag = make_fragment(
            exons=(
                _exon("chr1", 150, 200),
                _exon("chr1", 350, 400),
            ),
            introns=(GenomicInterval("chr1", 200, 350, Strand.POS),),
        )
        result = resolve_fragment(frag, overlap_index)
        assert result is not None
        # SJ (200,350) matches t_two_exon → SPLICED_ANNOT, narrowed to {t_two_exon}
        assert result.splice_type == int(SpliceType.SPLICED_ANNOT)
        t2e = tm["t_two_exon"]
        assert result.t_inds == frozenset({t2e})
        # exon(150-200)=50 + exon(350-400)=50 = 100bp exon
        assert result.exon_bp[t2e] == 100

    def test_intronic_cross_transcript_exons(self, tmp_path_factory):
        """Fragment in t0's intronic region covered by t1's exon → unambig=0.

        t0: exons (100,200),(400,500) — intron (200,400)
        t1: exon  (100,500)           — covers t0's entire intron
        Fragment (250,350): t0 exon_bp=0, t1 exon_bp=100.
        """
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t101\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "G1"; gene_type "protein_coding";
        """)
        idx = build_test_index(tmp_path_factory, gtf, name="cross_tx")
        tm = _t_map(idx)
        frag = make_fragment(exons=(_exon("chr1", 250, 350),), introns=())
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 100
        assert result.exon_bp[tm["t0"]] == 0
        assert result.exon_bp[tm["t1"]] == 100

    def test_intronic_partial_exon_overlap(self, tmp_path_factory):
        """Partial overlap: some intronic bp are exonic for another transcript.

        t0: exons (100,200),(400,500) — intron (200,400)
        t1: exon  (300,400) — covers only part of t0's intron
        Fragment (250,350): t0 has exon_bp=0, t1 has exon_bp=50
        """
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t301\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "G1"; gene_type "protein_coding";
        """)
        idx = build_test_index(tmp_path_factory, gtf, name="partial_ovl")
        tm = _t_map(idx)
        frag = make_fragment(exons=(_exon("chr1", 250, 350),), introns=())
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 100
        # t0: wholly intronic
        assert result.exon_bp[tm["t0"]] == 0
        # t1: exon(300-400) clipped to (300,350)=50bp
        assert result.exon_bp[tm["t1"]] == 50


def test_ambiguous_multiexon_candidate_snapshot(tmp_path_factory):
    """Many candidates across many exon blocks resolve by sorted-vector intersection."""
    idx = build_test_index(
        tmp_path_factory,
        _ambiguous_multiexon_gtf(),
        genome_size=1000,
        name="ambig_multiexon",
    )
    frag = make_fragment(
        exons=(
            _exon("chr1", 110, 140),
            _exon("chr1", 210, 240),
            _exon("chr1", 310, 340),
            _exon("chr1", 410, 440),
        ),
        introns=(),
    )

    result = resolve_fragment(frag, idx)
    assert result is not None
    by_index = dict(zip(idx.t_df["t_index"], idx.t_df["t_id"]))
    resolved_ids = {by_index[t] for t in result.t_inds}

    assert resolved_ids == {
        "t0",
        "t1",
        "t2",
        "t3",
        "t4",
        "t5",
        "t6",
        "t7",
        "RIGEL_NRNA_chr1_1_99_657",
    }
    assert result.splice_type == int(SpliceType.SPLICED_IMPLICIT)
    assert result.read_length == 120
    assert all(result.exon_bp[t] == 120 for t in result.t_inds)


# =====================================================================
# Overlap filtering integration in resolve_fragment
# =====================================================================


class TestOverlapFiltering:
    """Verify overlap filtering behavior for spliced vs unspliced fragments."""

    def test_unspliced_fragment_retains_all_candidates(self, overlap_index):
        """An unspliced fragment spanning exon + intron retains all candidates.

        Fragment (150,250) overlaps t_two_exon and t_one_left via their
        shared exon (100,200).  Both are kept even though t_one_left has
        less overlap — the EM handles weighting.
        """
        tm = _t_map(overlap_index)
        frag = make_fragment(exons=(_exon("chr1", 150, 250),), introns=())
        result = resolve_fragment(frag, overlap_index)
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)
        assert tm["t_two_exon"] in result.t_inds
        assert tm["t_one_left"] in result.t_inds

    def test_equal_overlap_keeps_all(self, overlap_index):
        """When all candidates have equal exon overlap, all are retained."""
        tm = _t_map(overlap_index)
        frag = make_fragment(exons=(_exon("chr1", 100, 200),), introns=())
        result = resolve_fragment(frag, overlap_index)
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)
        # Both t_two_exon and t_one_left share exon (100,200) fully
        assert tm["t_two_exon"] in result.t_inds
        assert tm["t_one_left"] in result.t_inds


# =====================================================================
# Fragment-length discrimination (3-exon gene, 2 isoforms)
# =====================================================================


class TestFragLengthDiscrimination:
    """Verify SJ-based isoform discrimination with the MINI_GTF 3-exon gene.

    MINI_GTF layout:
      t1(idx=0): exons (99,200),(299,400),(499,600) → SJs (200,299),(400,499)
      t2(idx=1): exons (99,200),(499,600)           → SJ  (200,499)

    A fragment's splice junction uniquely identifies which isoform
    produced it when the SJs differ.
    """

    def test_spliced_read_t2_specific_sj(self, mini_index):
        """Read with SJ (200,499) — specific to t2 (2-exon, skips middle)."""
        tm = _t_map(mini_index)
        frag = make_fragment(
            exons=(_exon("chr1", 150, 200), _exon("chr1", 499, 550)),
            introns=(GenomicInterval("chr1", 200, 499, Strand.POS),),
        )
        result = resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_ANNOT)
        assert result.t_inds == frozenset({tm["t2"]})

    def test_unspliced_read_in_shared_exon(self, mini_index):
        """Unspliced read fully in shared exon1 → t1, t2, and the nascent span."""
        tm = _t_map(mini_index)
        frag = make_fragment(exons=(_exon("chr1", 120, 180),), introns=())
        result = resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)
        assert {tm["t1"], tm["t2"]}.issubset(result.t_inds)

    def test_unspliced_read_in_t1_only_exon(self, mini_index):
        """Unspliced read in exon2 (299,400) — exonic for t1 only.

        t2's transcript span (99,600) also covers this region via
        intronic overlap, and the nascent span covers the full
        span, so all three are candidates; scoring differentiates.
        """
        tm = _t_map(mini_index)
        frag = make_fragment(exons=(_exon("chr1", 320, 380),), introns=())
        result = resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)
        assert {tm["t1"], tm["t2"]}.issubset(result.t_inds)

    def test_frag_length_differs_between_isoforms(self, mini_index):
        """Fragment spanning exon1→exon3 via SJ (200,499) → t2 only."""
        tm = _t_map(mini_index)
        frag = make_fragment(
            exons=(_exon("chr1", 150, 200), _exon("chr1", 499, 550)),
            introns=(GenomicInterval("chr1", 200, 499, Strand.POS),),
        )
        result = resolve_fragment(frag, mini_index)
        assert result is not None
        # SJ (200,499) → t2 specifically
        assert result.t_inds == frozenset({tm["t2"]})
