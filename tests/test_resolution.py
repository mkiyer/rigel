"""Tests for rigel.resolution — chimera detection and fragment resolution."""

import textwrap

import pytest
from conftest import build_test_index

from rigel.types import ChimeraType, MergeOutcome, Strand, GenomicInterval
from rigel.splice import SpliceType
from _resolution_reference import (
    _detect_intrachromosomal_chimera,
    make_fragment,
    resolve_fragment,
)


# =====================================================================
# _detect_intrachromosomal_chimera
# =====================================================================


class TestDetectIntrachromosomalChimera:
    """Transcript-set disjointness detection for intrachromosomal chimeras.

    ⚠ ``max_fragment_length`` is passed EXPLICITLY on every call, and the value is part of each case.
    Since 2026-08-19 compatibility is checked before a candidate is called a chimera
    (`TestGenomicCompatibilityIsCheckedBeforeChimera`), so a SAME-STRAND candidate is only a chimera when
    its implied fragment length exceeds the limit — those cases pass a limit BELOW their own span so they
    still exercise the disjointness mechanism. The DIFFERENT-STRAND and already-``None`` cases are decided
    before the compatibility check and the value cannot reach them."""

    def test_single_block_none(self):
        """One annotated block cannot be chimeric."""
        blocks = (GenomicInterval("chr1", 100, 200, Strand.POS),)
        t_sets = [frozenset({0, 1})]
        assert _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=1000) is None

    def test_empty_t_sets_none(self):
        """All empty transcript sets → not chimeric."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 300, 400, Strand.POS),
        )
        t_sets = [frozenset(), frozenset()]
        assert _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=1000) is None

    def test_one_empty_one_annotated(self):
        """One annotated, one intergenic → not chimeric."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 300, 400, Strand.POS),
        )
        t_sets = [frozenset({0}), frozenset()]
        assert _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=1000) is None

    def test_shared_transcript_none(self):
        """Two blocks sharing a transcript → connected → not chimeric."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 300, 400, Strand.POS),
        )
        t_sets = [frozenset({0, 1}), frozenset({1, 2})]
        assert _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=1000) is None

    def test_disjoint_same_strand_chimera(self):
        """Two blocks with disjoint transcript sets, same strand → CIS_STRAND_SAME."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 1000, 1100, Strand.POS),
        )
        t_sets = [frozenset({0, 1}), frozenset({2, 3})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=500)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_SAME
        assert chimera_gap == 800  # 1000 - 200

    def test_disjoint_diff_strand_chimera(self):
        """Two blocks with disjoint sets, different strands → CIS_STRAND_DIFF."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 500, 600, Strand.NEG),
        )
        t_sets = [frozenset({0}), frozenset({5})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=1000)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_DIFF
        assert chimera_gap == 300  # 500 - 200

    def test_overlapping_disjoint_blocks(self):
        """Disjoint transcript sets with overlapping genomic coordinates."""
        blocks = (
            GenomicInterval("chr1", 100, 500, Strand.POS),
            GenomicInterval("chr1", 300, 700, Strand.NEG),
        )
        t_sets = [frozenset({0}), frozenset({1})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=1000)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_DIFF
        assert chimera_gap == 0  # overlapping

    def test_three_blocks_two_components(self):
        """Three blocks: two connected (shared tx), one disjoint."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 300, 400, Strand.POS),
            GenomicInterval("chr1", 2000, 2100, Strand.NEG),
        )
        # Block 0 and 1 share transcript 1 → connected
        # Block 2 is disjoint from both → chimeric
        t_sets = [frozenset({0, 1}), frozenset({1, 2}), frozenset({5})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=1000)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_DIFF
        # Gap: min(2000-200, 2000-400) = 1600
        assert chimera_gap == 1600

    def test_three_disjoint_blocks(self):
        """Three blocks, all pairwise disjoint → three components → chimeric."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 500, 600, Strand.POS),
            GenomicInterval("chr1", 1000, 1100, Strand.POS),
        )
        t_sets = [frozenset({0}), frozenset({1}), frozenset({2})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=500)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_SAME
        # Gap: min(500-200, 1000-200, 1000-600) = 300
        assert chimera_gap == 300

    def test_connected_via_transitive_overlap(self):
        """A-B overlap, B-C overlap, but A-C are disjoint → all connected."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 300, 400, Strand.POS),
            GenomicInterval("chr1", 500, 600, Strand.POS),
        )
        # A={0,1}, B={1,2}, C={2,3} → A-B share 1, B-C share 2 → all connected
        t_sets = [frozenset({0, 1}), frozenset({1, 2}), frozenset({2, 3})]
        assert _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=1000) is None

    def test_adjacent_blocks_zero_gap(self):
        """Disjoint transcript sets with blocks that are adjacent (gap=0)."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 200, 300, Strand.POS),
        )
        t_sets = [frozenset({0}), frozenset({1})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=100)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_SAME
        assert chimera_gap == 0  # adjacent, no gap

    def test_mixed_annotated_intergenic_blocks(self):
        """Annotated blocks are disjoint, intergenic block between them."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 500, 600, Strand.NONE),  # intergenic
            GenomicInterval("chr1", 1000, 1100, Strand.NEG),
        )
        # Block 1 is intergenic (empty set), blocks 0 and 2 are disjoint
        t_sets = [frozenset({0}), frozenset(), frozenset({5})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets, max_fragment_length=1000)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_DIFF
        assert chimera_gap == 800  # 1000 - 200


# =====================================================================
# resolve_fragment — real C++ FragmentResolver tests
# =====================================================================
#
# All tests below use real TranscriptIndex instances built from GTF strings.
# The C++ FragmentResolver handles interval queries, SJ matching, and
# overlap computation.  Mock indexes have been removed.
#
# Each test class defines a session-scoped fixture that builds the
# specific TranscriptIndex it needs.  A shared mini_index fixture (from
# conftest.py) is reused when the MINI_GTF layout suffices.
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
    """Test per-candidate overlap profiles produced by resolve_fragment.

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
        assert result.overlap_bp[tm["t_two_exon"]] == (50, 0)
        assert result.overlap_bp[tm["t_one_left"]] == (50, 0)

    def test_intronic_overlap(self, overlap_index):
        """Fragment in intronic region → exon_bp=0, intron_bp=read_length."""
        tm = _t_map(overlap_index)
        frag = make_fragment(exons=(_exon("chr1", 220, 320),), introns=())
        result = resolve_fragment(frag, overlap_index)
        assert result is not None
        assert result.read_length == 100
        # Only t_two_exon's transcript span covers 220-320
        t2e = tm["t_two_exon"]
        assert t2e in result.t_inds
        assert result.overlap_bp[t2e] == (0, 100)

    def test_mixed_exon_intron_overlap(self, overlap_index):
        """Fragment spanning exon-intron boundary → partial exon + intron bp."""
        tm = _t_map(overlap_index)
        frag = make_fragment(exons=(_exon("chr1", 180, 220),), introns=())
        result = resolve_fragment(frag, overlap_index)
        assert result is not None
        assert result.read_length == 40
        # t_two_exon: exon(100-200) covers 180-200=20bp; tx span covers all 40bp
        #   intron_bp=40-20=20
        assert result.overlap_bp[tm["t_two_exon"]] == (20, 20)
        # t_one_left: exon(100-200) covers 180-200=20bp; tx ends at 200 → 20bp tx
        #   intron_bp=0
        assert result.overlap_bp[tm["t_one_left"]] == (20, 0)

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
        assert result.overlap_bp[t2e] == (100, 0)

    def test_intronic_cross_transcript_exons(self, tmp_path_factory):
        """Fragment in t0's intronic region covered by t1's exon → unambig=0.

        t0: exons (100,200),(400,500) — intron (200,400)
        t1: exon  (100,500)           — covers t0's entire intron
        Fragment (250,350): t0 intron_bp=100, t1 exon_bp=100.
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
        assert result.overlap_bp[tm["t0"]] == (0, 100)
        assert result.overlap_bp[tm["t1"]] == (100, 0)

    def test_intronic_partial_exon_overlap(self, tmp_path_factory):
        """Partial overlap: some intronic bp are exonic for another transcript.

        t0: exons (100,200),(400,500) — intron (200,400)
        t1: exon  (300,400) — covers only part of t0's intron
        Fragment (250,350): t0 has intron_bp=100, t1 has exon_bp=50
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
        # t0: intron_bp=100
        assert result.overlap_bp[tm["t0"]] == (0, 100)
        # t1: exon(300-400) clipped to (300,350)=50bp
        assert result.overlap_bp[tm["t1"]] == (50, 0)


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
    assert result.merge_criteria == int(MergeOutcome.UNION)
    assert result.splice_type == int(SpliceType.SPLICED_IMPLICIT)
    assert result.chimera_type == int(ChimeraType.NONE)
    assert result.read_length == 120
    assert all(result.overlap_bp[t] == (120, 0) for t in result.t_inds)


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

    def test_spliced_fragment_sj_narrows_candidates(self, mini_index):
        """Spliced fragments with annotated SJs: SJ match narrows first."""
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
        # SJ (200,299) → t1 only
        assert result.t_inds == frozenset({tm["t1"]})

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

    def test_spliced_read_t1_specific_sj(self, mini_index):
        """Read with SJ (200,299) — specific to t1 (3-exon)."""
        tm = _t_map(mini_index)
        frag = make_fragment(
            exons=(_exon("chr1", 150, 200), _exon("chr1", 299, 350)),
            introns=(GenomicInterval("chr1", 200, 299, Strand.POS),),
        )
        result = resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_ANNOT)
        assert result.t_inds == frozenset({tm["t1"]})

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
        """Unspliced read fully in shared exon1 → t1, t2, and synthetic nRNA."""
        tm = _t_map(mini_index)
        frag = make_fragment(exons=(_exon("chr1", 120, 180),), introns=())
        result = resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)
        assert {tm["t1"], tm["t2"]}.issubset(result.t_inds)

    def test_unspliced_read_in_t1_only_exon(self, mini_index):
        """Unspliced read in exon2 (299,400) — exonic for t1 only.

        t2's transcript span (99,600) also covers this region via
        intronic overlap, and the synthetic nRNA exon covers the full
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


# =====================================================================
# Intronic bp regression tests
# =====================================================================


class TestIntronicBpAccumulation:
    """Verify resolve_fragment correctly computes exon_bp and intron_bp."""

    def test_intron_bp_for_lone_intron(self, tmp_path_factory):
        """Fragment in a transcript's intron (single transcript).

        t0: exons (100,200),(400,500) → intron (200,400)
        Fragment (250,350) → intron_bp = 100.
        """
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
        """)
        idx = build_test_index(tmp_path_factory, gtf, name="lone_intron")
        frag = make_fragment(exons=(_exon("chr1", 250, 350),), introns=())
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 100
        assert result.overlap_bp[0] == (0, 100)

    def test_intron_bp_when_cross_exon_covers(self, tmp_path_factory):
        """Fragment in t0's intron entirely covered by t1's exon.

        t0: exons (100,200),(400,500) — intron (200,400)
        t1: exon  (200,400) — covers t0's entire intron
        Fragment (250,350) → t0 intron_bp=100, t1 exon_bp=100.
        """
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "G1"; gene_type "protein_coding";
        """)
        idx = build_test_index(tmp_path_factory, gtf, name="cross_exon")
        tm = _t_map(idx)
        frag = make_fragment(exons=(_exon("chr1", 250, 350),), introns=())
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 100
        assert result.overlap_bp[tm["t0"]] == (0, 100)
        assert result.overlap_bp[tm["t1"]] == (100, 0)

    def test_intron_partial(self, tmp_path_factory):
        """Fragment partially intronic for t0, partially exonic for t1.

        t0: exons (100,200),(400,500) — intron (200,400)
        t1: exon  (300,400) — covers only part of t0's intron
        Fragment (250,350) → t0 intron_bp=100, t1 exon_bp=50
        """
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t301\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "G1"; gene_type "protein_coding";
        """)
        idx = build_test_index(tmp_path_factory, gtf, name="partial_unambig")
        tm = _t_map(idx)
        frag = make_fragment(exons=(_exon("chr1", 250, 350),), introns=())
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 100
        assert result.overlap_bp[tm["t0"]] == (0, 100)
        assert result.overlap_bp[tm["t1"]] == (50, 0)

    def test_exon_only_fragment_has_zero_intron_bp(self, tmp_path_factory):
        """Fragment fully within an exon → intron_bp=0.

        t0: exons (100,300),(400,500) — intron (300,400)
        Fragment (150,250) → fully in first exon.
        """
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t0"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
        """)
        idx = build_test_index(tmp_path_factory, gtf, name="exon_only")
        frag = make_fragment(exons=(_exon("chr1", 150, 250),), introns=())
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 100
        assert result.overlap_bp[0] == (100, 0)


# =====================================================================
# ⭐⭐⭐ COMPATIBILITY BEFORE CHIMERA — owner ruling, 2026-08-19
# =====================================================================


class TestGenomicCompatibilityIsCheckedBeforeChimera:
    """⛔⛔ **A CONTIGUOUS GENOMIC MOLECULE IS NOT A CHIMERA, HOWEVER MANY TRANSCRIPTS IT SPANS.**

    Transcript-set disjointness alone called ``CHIMERA_CIS_STRAND_SAME`` on ordinary gDNA: a genomic
    fragment is contiguous and routinely covers two transcripts that share nothing, which is what the
    annotation looks like there rather than evidence of a rearrangement. Every such fragment was a proper
    inward-facing pair with no junction anywhere in it — a ``100M``/``100M`` ``NH:i:1`` pair — so there was
    no splice evidence for a chimera to be inferred FROM.

    ⭐ **The rule (owner):** a fragment is a chimera only if the mates are genomically INCOMPATIBLE —
    a different reference, an orientation that is not facing inward, or an implied fragment length beyond
    the library's ``max_frag_length``. Compatibility is checked FIRST.

    ⚠ **Measured cost of the old predicate**, on ``gdna_g98_ss_0.99_nrna_none_capture_off``: **4,087**
    gDNA fragments dropped per condition — 0.04 % of fragments carrying **2.4 % of every boundary
    crossing**, because the predicate fires exactly where transcripts are short and dense, which is
    exactly where a fragment crosses many boundaries (3.66 crossings lost each, against 1.65 for an
    average crosser). That is the whole boundary-crossing deficit that failed FIELD certification on the
    8 dense capture-OFF conditions.

    ⭐ The orientation half needs no new predicate: ``build_fragment`` keys blocks by
    ``(ref_id, ref_strand)`` with R2's orientation flipped, so both mates of an inward-facing pair carry
    the SAME ``strand`` and ``unique_strands.size() == 1`` already MEANS "facing inward".
    """

    #: two blocks, disjoint transcript sets, same strand — implied fragment span 100..1100 = 1000 bp
    BLOCKS = (
        GenomicInterval("chr1", 100, 200, Strand.POS),
        GenomicInterval("chr1", 1000, 1100, Strand.POS),
    )
    T_SETS = [frozenset({0, 1}), frozenset({2, 3})]

    def test_a_compatible_pair_is_not_a_chimera(self):
        """⭐ Same reference, facing inward, span within the maximum ⇒ an ordinary genomic molecule."""
        assert (
            _detect_intrachromosomal_chimera(self.BLOCKS, self.T_SETS, max_fragment_length=2000)
            is None
        ), "a contiguous, inward-facing, in-range genomic fragment was called a chimera"

    def test_an_impossible_span_is_still_a_chimera(self):
        """⛔ The length half of the rule: beyond ``max_frag_length`` no molecule explains the pair."""
        result = _detect_intrachromosomal_chimera(self.BLOCKS, self.T_SETS, max_fragment_length=500)
        assert result is not None, "a 1,000 bp span passed a 500 bp limit"
        assert result[0] == ChimeraType.CIS_STRAND_SAME

    def test_the_span_is_the_OUTER_extent_not_the_gap(self):
        """⚠ The quantity is the implied FRAGMENT LENGTH — outermost start to outermost end — not the
        ``chimera_gap`` between the blocks. The two differ by the blocks' own lengths (here 1000 vs 800),
        and using the gap would admit a fragment 200 bp longer than the library can contain."""
        assert (
            _detect_intrachromosomal_chimera(self.BLOCKS, self.T_SETS, max_fragment_length=900)
            is not None
        ), "the gap (800) was used where the implied fragment length (1000) is the rule"

    def test_an_outward_facing_pair_is_a_chimera_at_any_span(self):
        """⛔ The orientation half: opposite ``strand`` on the two components means the mates do not face
        inward, and no fragment length rescues that."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 1000, 1100, Strand.NEG),
        )
        result = _detect_intrachromosomal_chimera(blocks, self.T_SETS, max_fragment_length=100_000)
        assert result is not None, "an outward-facing pair was rescued by the length check"
        assert result[0] == ChimeraType.CIS_STRAND_DIFF

    def test_the_NATIVE_resolver_agrees(self, tmp_path_factory):
        """⭐⭐ The C++ kernel, not the reference — this is the path the scanner actually runs.

        Two transcripts 200 bp apart sharing nothing. A gDNA fragment covering one block in each is a
        perfectly ordinary molecule and must resolve with ``chimera_type == NONE``."""
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "gA"; transcript_id "tA"; gene_name "GA"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id "gB"; transcript_id "tB"; gene_name "GB"; gene_type "protein_coding"; tag "basic";
        """)
        idx = build_test_index(tmp_path_factory, gtf, name="compat_chimera")
        idx.resolver.set_max_fragment_length(1000)
        frag = make_fragment(
            exons=(
                GenomicInterval("chr1", 120, 200, Strand.POS),
                GenomicInterval("chr1", 400, 480, Strand.POS),
            )
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.chimera_type == int(ChimeraType.NONE), (
            "the native resolver called a 360 bp inward-facing genomic fragment a chimera"
        )

    def test_the_NATIVE_resolver_measures_the_SPAN_not_the_gap(self, tmp_path_factory):
        """⛔ The native twin of ``test_the_span_is_the_OUTER_extent_not_the_gap``, and it exists because
        breaking the C++ to use ``min_gap`` fired NO native gate (verified by perturbation).

        Blocks (100,200) and (1000,1100): implied fragment length **1000**, gap between them **800**. At a
        limit of 900 the fragment is impossible and the gap is not — so a gap-based rule rescues a molecule
        100 bp longer than the library can contain."""
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "gA"; transcript_id "tA"; gene_name "GA"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t1001\t1100\t.\t+\t.\tgene_id "gB"; transcript_id "tB"; gene_name "GB"; gene_type "protein_coding"; tag "basic";
        """)
        idx = build_test_index(tmp_path_factory, gtf, name="compat_span_vs_gap")
        idx.resolver.set_max_fragment_length(900)
        frag = make_fragment(
            exons=(
                GenomicInterval("chr1", 100, 200, Strand.POS),
                GenomicInterval("chr1", 1000, 1100, Strand.POS),
            )
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.chimera_type == int(ChimeraType.CIS_STRAND_SAME), (
            "a 1,000 bp implied fragment passed a 900 bp limit — the GAP (800) was measured, not the span"
        )

    def test_the_NATIVE_resolver_still_refuses_an_impossible_span(self, tmp_path_factory):
        """⛔ The control: the same geometry with a limit the span cannot meet stays a chimera, so the
        gate above is not passing because the rescue is unconditional."""
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "gA"; transcript_id "tA"; gene_name "GA"; gene_type "protein_coding"; tag "basic";
            chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id "gB"; transcript_id "tB"; gene_name "GB"; gene_type "protein_coding"; tag "basic";
        """)
        idx = build_test_index(tmp_path_factory, gtf, name="compat_chimera_ctl")
        idx.resolver.set_max_fragment_length(100)
        frag = make_fragment(
            exons=(
                GenomicInterval("chr1", 120, 200, Strand.POS),
                GenomicInterval("chr1", 400, 480, Strand.POS),
            )
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.chimera_type == int(ChimeraType.CIS_STRAND_SAME)
