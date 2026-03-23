"""Tests for transcript-space fragment length computation.

Validates that compute_frag_lengths() correctly projects fragment genomic
endpoints to transcript coordinate space and computes FL as the distance.

Tests cover:
  - Basic spliced and unspliced fragments
  - Transcript boundary overhang (before first exon, past last exon)
  - Endpoint overhang into internal introns
  - Internal intronic overhang (safe — not at endpoints)
  - Both-end overhang
  - Multi-intron gap correction (large unsequenced gap spanning multiple introns)
  - nRNA (single-exon) candidates
  - Negative-strand transcripts
  - Multiple candidate transcripts with different exon structure
"""

import textwrap

import pytest

from rigel.types import Strand, GenomicInterval
from rigel.resolution import make_fragment, resolve_fragment
from conftest import build_test_index


def _t_map(index):
    """Return dict mapping transcript ID string to integer index."""
    return dict(zip(index.t_df["t_id"], index.t_df["t_index"]))


def _exon(ref, start, end, strand=Strand.POS):
    return GenomicInterval(ref, start, end, strand)


# =====================================================================
# GTF definitions for test indexes
# =====================================================================

# Two-exon transcript for basic gap tests
# 0-based half-open after parse:
#   t_basic: exons (999, 2000), (4999, 6000)
#   transcript length = 1001 + 1001 = 2002
#   intron = (2000, 4999)
BASIC_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t1000\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t_basic"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t5000\t6000\t.\t+\t.\tgene_id "g1"; transcript_id "t_basic"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
""")

# Three-exon transcript for multi-intron gap tests
# 0-based half-open after parse:
#   t_three: exons (999, 2000), (4999, 6000), (8999, 10000)
#   transcript length = 1001 + 1001 + 1001 = 3003
#   introns = (2000, 4999), (6000, 8999)
THREE_EXON_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t1000\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t_three"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t5000\t6000\t.\t+\t.\tgene_id "g1"; transcript_id "t_three"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t9000\t10000\t.\t+\t.\tgene_id "g1"; transcript_id "t_three"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
""")

# Four-exon transcript for even larger gap tests
# 0-based half-open:
#   t_four: exons (999,2000), (4999,6000), (8999,10000), (12999,14000)
#   transcript length = 1001 * 4 = 4004
#   introns = (2000,4999), (6000,8999), (10000,12999)
FOUR_EXON_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t1000\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t_four"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t5000\t6000\t.\t+\t.\tgene_id "g1"; transcript_id "t_four"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t9000\t10000\t.\t+\t.\tgene_id "g1"; transcript_id "t_four"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t13000\t14000\t.\t+\t.\tgene_id "g1"; transcript_id "t_four"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
""")

# Two isoforms with different exon structures
# 0-based half-open after parse:
#   t_long:  exons (999, 2000), (2999, 4000), (4999, 6000) — 3 exons
#   t_short: exons (999, 2000), (4999, 6000)               — 2 exons (skips middle)
#   transcript lengths: t_long = 3003, t_short = 2002
ISOFORM_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t1000\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t_long"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t3000\t4000\t.\t+\t.\tgene_id "g1"; transcript_id "t_long"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t5000\t6000\t.\t+\t.\tgene_id "g1"; transcript_id "t_long"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t1000\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t_short"; gene_name "G1"; gene_type "protein_coding";
    chr1\ttest\texon\t5000\t6000\t.\t+\t.\tgene_id "g1"; transcript_id "t_short"; gene_name "G1"; gene_type "protein_coding";
""")

# Negative-strand transcript
# 0-based half-open after parse:
#   t_neg: exons (999, 2000), (4999, 6000)
NEG_STRAND_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t1000\t2000\t.\t-\t.\tgene_id "g1"; transcript_id "t_neg"; gene_name "G1"; gene_type "lncRNA";
    chr1\ttest\texon\t5000\t6000\t.\t-\t.\tgene_id "g1"; transcript_id "t_neg"; gene_name "G1"; gene_type "lncRNA";
""")


# =====================================================================
# Fixtures
# =====================================================================


@pytest.fixture(scope="session")
def basic_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, BASIC_GTF, genome_size=20000, name="fl_basic")


@pytest.fixture(scope="session")
def three_exon_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, THREE_EXON_GTF, genome_size=20000, name="fl_three")


@pytest.fixture(scope="session")
def four_exon_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, FOUR_EXON_GTF, genome_size=20000, name="fl_four")


@pytest.fixture(scope="session")
def isoform_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, ISOFORM_GTF, genome_size=20000, name="fl_isoform")


@pytest.fixture(scope="session")
def neg_strand_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, NEG_STRAND_GTF, genome_size=20000, name="fl_neg")


# =====================================================================
# Helper to extract FL for a specific transcript from resolve result
# =====================================================================


def _get_fl(result, index, t_id):
    """Get fragment length for transcript t_id from a ResolvedFragment."""
    tm = _t_map(index)
    t_idx = tm[t_id]
    fl_dict = result.frag_lengths
    return fl_dict.get(t_idx)


# =====================================================================
# Basic fragment length tests
# =====================================================================


class TestBasicFragLength:
    """Basic FL computation for simple fragment/transcript configurations."""

    def test_single_block_fragment(self, basic_index):
        """Single-block fragment: FL = block length, same for all candidates."""
        frag = make_fragment(exons=(_exon("chr1", 1200, 1400),), introns=())
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 200

    def test_spliced_fragment_observed_intron(self, basic_index):
        """Spliced fragment with observed intron: FL = exonic bases only.

        t_basic exons: (999,2000) and (4999,6000), intron (2000,4999)
        Fragment: blocks (1500,2000) + (4999,5500), intron (2000,4999)
        Genomic footprint = 5500 - 1500 = 4000
        Expected FL = 500 + 501 = 1001 (only exonic bases)
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1500, 2000),
                _exon("chr1", 4999, 5500),
            ),
            introns=(GenomicInterval("chr1", 2000, 4999, Strand.POS),),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 1001

    def test_unspliced_fragment_within_single_exon(self, basic_index):
        """Unspliced paired-end in single exon: FL = block span.

        Fragment: blocks (1200,1350) + (1600,1750), both within exon (999,2000)
        Genomic footprint = 1750 - 1200 = 550
        Expected FL = 550 (all exonic, no intron to subtract)
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1200, 1350),
                _exon("chr1", 1600, 1750),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 550  # 1750 - 1200


# =====================================================================
# Transcript boundary overhang tests
# =====================================================================


class TestTranscriptBoundaryOverhang:
    """FL must count bases that overhang beyond transcript exon boundaries."""

    def test_start_overhang_before_first_exon(self, basic_index):
        """Fragment starts before the first exon of the transcript.

        t_basic exons: (999, 2000), (4999, 6000)
        Fragment: blocks (994, 1100) + (5500, 5600)
        gstart = 994, gend = 5600
        tx_pos(994) = 994 - 999 = -5 (5 bp before first exon)
        tx_pos(5600) = 1001 + 601 = 1602 (1001 from exon1 + 601 into exon2)
        FL = |1602 - (-5)| = 1607
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 994, 1100),
                _exon("chr1", 5500, 5600),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 1607

    def test_end_overhang_past_last_exon(self, basic_index):
        """Fragment extends past the last exon of the transcript.

        t_basic exons: (999, 2000), (4999, 6000)
        Fragment: blocks (1500, 1600) + (5900, 6005)
        gstart = 1500, gend = 6005
        tx_pos(1500) = 501
        tx_pos(6005) = 2002 + 5 = 2007 (t_len=2002, 5bp past last exon)
        FL = 2007 - 501 = 1506
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1500, 1600),
                _exon("chr1", 5900, 6005),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 1506

    def test_both_end_overhang(self, basic_index):
        """Fragment overhangs both before first and past last exon.

        t_basic exons: (999, 2000), (4999, 6000)
        Fragment: blocks (996, 1050) + (5980, 6003)
        gstart = 996, gend = 6003
        tx_pos(996) = 996 - 999 = -3
        tx_pos(6003) = 2002 + 3 = 2005
        FL = |2005 - (-3)| = 2008
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 996, 1050),
                _exon("chr1", 5980, 6003),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 2008


# =====================================================================
# Endpoint vs internal intronic overhang tests
# =====================================================================


class TestIntronicOverhang:
    """Endpoint intronic overhang must be counted; internal is safe to ignore."""

    def test_endpoint_overhang_into_intron(self, basic_index):
        """Fragment gend extends past exon1 into the intron.

        t_basic exons: (999, 2000), (4999, 6000)
        Fragment: blocks (1600, 1750) + (1855, 2005)
        gstart = 1600, gend = 2005
        tx_pos(1600) = 601
        tx_pos(2005) = 1001 + 5 = 1006 (exon1 ends at tx offset 1001, +5 overhang)
        FL = 1006 - 601 = 405
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1600, 1750),
                _exon("chr1", 1855, 2005),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 405

    def test_internal_intronic_overhang_not_inflating(self, basic_index):
        """Internal block overhangs into intron (not at gstart/gend) don't inflate FL.

        t_basic exons: (999, 2000), (4999, 6000), intron (2000, 4999)
        Fragment: blocks (1855, 2005), (4995, 5149)
        The first block extends 5bp into the intron (at 2005) — but this is
        NOT gend (gend=5149). The second block starts 4bp before exon2
        (at 4995) — but this is NOT gstart (gstart=1855).
        gstart = 1855, gend = 5149
        tx_pos(1855) = 856 (inside exon1)
        tx_pos(5149) = 1001 + 150 = 1151 (inside exon2)
        FL = 1151 - 856 = 295
        Internal block overhangs at 2005 and 4995 are irrelevant.
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1855, 2005),
                _exon("chr1", 4995, 5149),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 295

    def test_gstart_in_intron(self, basic_index):
        """Fragment gstart is in the intron between exons.

        t_basic exons: (999, 2000), (4999, 6000), intron (2000, 4999)
        Fragment: blocks (2003, 2100) + (5200, 5400)
        gstart = 2003, gend = 5400
        tx_pos(2003) = 1001 + 3 = 1004 (3bp past exon1 end into intron)
        tx_pos(5400) = 1001 + 401 = 1402 (inside exon2)
        FL = 1402 - 1004 = 398
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 2003, 2100),
                _exon("chr1", 5200, 5400),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 398


# =====================================================================
# Multi-intron gap correction tests
# =====================================================================


class TestMultiIntronGap:
    """Large unsequenced gap spanning multiple introns."""

    def test_gap_spanning_one_intron_three_exon(self, three_exon_index):
        """Paired-end with gap spanning one intron of a 3-exon transcript.

        t_three exons: (999,2000), (4999,6000), (8999,10000)
        Fragment: blocks (1500, 1600) + (5500, 5600)
        gstart=1500, gend=5600
        tx_pos(1500) = 501
        tx_pos(5600) = 1001 + 601 = 1602
        FL = 1602 - 501 = 1101
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1500, 1600),
                _exon("chr1", 5500, 5600),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, three_exon_index)
        assert result is not None
        fl = _get_fl(result, three_exon_index, "t_three")
        assert fl == 1101

    def test_gap_spanning_two_introns(self, three_exon_index):
        """Paired-end with gap spanning TWO introns (exon1 → exon3).

        t_three exons: (999,2000), (4999,6000), (8999,10000)
        Fragment: blocks (1800, 1900) + (9200, 9300)
        gstart=1800, gend=9300
        Genomic footprint = 9300 - 1800 = 7500
        tx_pos(1800) = 801
        tx_pos(9300) = 1001 + 1001 + 301 = 2303
        FL = 2303 - 801 = 1502
        Without correction, genomic FL would be 7500 — wrong by ~6000.
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1800, 1900),
                _exon("chr1", 9200, 9300),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, three_exon_index)
        assert result is not None
        fl = _get_fl(result, three_exon_index, "t_three")
        assert fl == 1502

    def test_gap_spanning_three_introns(self, four_exon_index):
        """Paired-end with gap spanning THREE introns (exon1 → exon4).

        t_four exons: (999,2000), (4999,6000), (8999,10000), (12999,14000)
        Fragment: blocks (1900, 1950) + (13500, 13550)
        gstart=1900, gend=13550
        Genomic footprint = 13550 - 1900 = 11650
        tx_pos(1900) = 901
        tx_pos(13550) = 1001 + 1001 + 1001 + 551 = 3554
        FL = 3554 - 901 = 2653
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1900, 1950),
                _exon("chr1", 13500, 13550),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, four_exon_index)
        assert result is not None
        fl = _get_fl(result, four_exon_index, "t_four")
        assert fl == 2653

    def test_spliced_fragment_plus_gap_intron(self, three_exon_index):
        """Read1 spliced (exon1→exon2), Read2 in exon3 — gap spans intron 2.

        t_three exons: (999,2000), (4999,6000), (8999,10000)
        Fragment: blocks (1800, 2000) + (4999, 5100) + (9200, 9300)
        introns: (2000, 4999) observed
        gstart=1800, gend=9300
        tx_pos(1800) = 801
        tx_pos(9300) = 1001 + 1001 + 301 = 2303
        FL = 2303 - 801 = 1502
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1800, 2000),
                _exon("chr1", 4999, 5100),
                _exon("chr1", 9200, 9300),
            ),
            introns=(GenomicInterval("chr1", 2000, 4999, Strand.POS),),
        )
        result = resolve_fragment(frag, three_exon_index)
        assert result is not None
        fl = _get_fl(result, three_exon_index, "t_three")
        assert fl == 1502

    def test_spliced_fragment_across_all_three_introns(self, four_exon_index):
        """Read1 spliced (exon1→exon2→exon3), Read2 in exon4 via gap.

        t_four exons: (999,2000), (4999,6000), (8999,10000), (12999,14000)
        Fragment: blocks (1900,2000) + (4999,6000) + (8999,9100) + (13500,13550)
        introns: (2000,4999), (6000,8999) observed
        gstart=1900, gend=13550
        tx_pos(1900) = 901
        tx_pos(13550) = 1001 + 1001 + 1001 + 551 = 3554
        FL = 3554 - 901 = 2653
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1900, 2000),
                _exon("chr1", 4999, 6000),
                _exon("chr1", 8999, 9100),
                _exon("chr1", 13500, 13550),
            ),
            introns=(
                GenomicInterval("chr1", 2000, 4999, Strand.POS),
                GenomicInterval("chr1", 6000, 8999, Strand.POS),
            ),
        )
        result = resolve_fragment(frag, four_exon_index)
        assert result is not None
        fl = _get_fl(result, four_exon_index, "t_four")
        assert fl == 2653


# =====================================================================
# Different FL per isoform (multi-transcript)
# =====================================================================


class TestIsoformFL:
    """FL differs between candidate transcripts with different exon structures."""

    def test_gap_fragment_different_fl_per_isoform(self, isoform_index):
        """Fragment with gap: FL depends on whether middle exon is present.

        t_long exons:  (999,2000), (2999,4000), (4999,6000) — 3 exons
        t_short exons: (999,2000), (4999,6000) — 2 exons (skips middle)

        Fragment: blocks (1900, 1950) + (5050, 5100)
        gstart=1900, gend=5100

        For t_long:
          tx_pos(1900) = 901
          tx_pos(5100) = 1001 + 1001 + 101 = 2103
          FL = 2103 - 901 = 1202

        For t_short:
          tx_pos(1900) = 901
          tx_pos(5100) = 1001 + 101 = 1102
          FL = 1102 - 901 = 201
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1900, 1950),
                _exon("chr1", 5050, 5100),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, isoform_index)
        assert result is not None
        fl_long = _get_fl(result, isoform_index, "t_long")
        fl_short = _get_fl(result, isoform_index, "t_short")
        assert fl_long == 1202
        assert fl_short == 201

    def test_spliced_fragment_narrows_to_one_isoform(self, isoform_index):
        """Fragment with SJ specific to t_long resolves to one isoform.

        t_long SJs: (2000,2999), (4000,4999)
        t_short SJ: (2000,4999)

        Fragment: blocks (1900,2000) + (2999,3100)
        intron: (2000,2999) → matches t_long only
        gstart=1900, gend=3100
        tx_pos(1900) = 901
        tx_pos(3100) = 1001 + 101 = 1102
        FL = 1102 - 901 = 201
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1900, 2000),
                _exon("chr1", 2999, 3100),
            ),
            introns=(GenomicInterval("chr1", 2000, 2999, Strand.POS),),
        )
        result = resolve_fragment(frag, isoform_index)
        assert result is not None
        fl_long = _get_fl(result, isoform_index, "t_long")
        assert fl_long == 201


# =====================================================================
# Negative-strand transcript
# =====================================================================


class TestNegativeStrand:
    """FL is strand-agnostic (physical length of the molecule)."""

    def test_neg_strand_basic_fl(self, neg_strand_index):
        """Fragment on negative-strand transcript: FL is same as positive strand.

        t_neg exons: (999, 2000), (4999, 6000) — same coords, negative strand
        Fragment: blocks (1500, 1600) + (5500, 5600)
        gstart=1500, gend=5600
        tx_pos(1500) = 501
        tx_pos(5600) = 1001 + 601 = 1602
        FL = 1602 - 501 = 1101
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1500, 1600, Strand.NEG),
                _exon("chr1", 5500, 5600, Strand.NEG),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, neg_strand_index)
        assert result is not None
        fl = _get_fl(result, neg_strand_index, "t_neg")
        assert fl == 1101

    def test_neg_strand_with_overhang(self, neg_strand_index):
        """Negative-strand fragment with start overhang.

        t_neg exons: (999, 2000), (4999, 6000)
        Fragment: blocks (994, 1100) + (5500, 5600)
        gstart=994, gend=5600
        tx_pos(994) = 994 - 999 = -5
        tx_pos(5600) = 1001 + 601 = 1602
        FL = |1602 - (-5)| = 1607
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 994, 1100, Strand.NEG),
                _exon("chr1", 5500, 5600, Strand.NEG),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, neg_strand_index)
        assert result is not None
        fl = _get_fl(result, neg_strand_index, "t_neg")
        assert fl == 1607


# =====================================================================
# nRNA (nascent RNA) — single synthetic exon
# =====================================================================


class TestNRNAFragLength:
    """nRNA uses a single synthetic exon covering the full genomic interval.

    For nRNA candidates, FL should equal the genomic span between
    fragment endpoints (no introns to skip).
    """

    def test_nrna_fl_equals_genomic_span(self, basic_index):
        """Fragment in a region covered by nRNA synthetic transcript.

        The build_test_index creates synthetic nRNA transcripts.
        nRNA for g1 covers the full gene span as a single exon.
        For a fragment entirely within the nRNA span,
        FL = genomic distance between endpoints.

        Fragment: blocks (1500, 1600) + (5500, 5600)
        gstart=1500, gend=5600
        For nRNA (single exon), all positions map linearly → FL = 5600-1500 = 4100
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1500, 1600),
                _exon("chr1", 5500, 5600),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        # Find the nRNA candidate (synthetic nRNA has "is_synthetic_nrna" flag)
        t_df = basic_index.t_df
        nrna_mask = t_df["is_synthetic_nrna"].values if "is_synthetic_nrna" in t_df.columns else None
        if nrna_mask is not None:
            nrna_tidxs = set(t_df.loc[nrna_mask.astype(bool), "t_index"].tolist())
            fl_dict = result.frag_lengths
            for t_idx in nrna_tidxs:
                if t_idx in fl_dict:
                    # nRNA: FL = genomic footprint between endpoints
                    assert fl_dict[t_idx] == 4100


# =====================================================================
# Combined overhang + multi-intron gap tests
# =====================================================================


class TestCombinedOverhangGap:
    """Overhang + gap spanning introns simultaneously."""

    def test_start_overhang_plus_two_intron_gap(self, three_exon_index):
        """Fragment starts before first exon and spans two introns.

        t_three exons: (999,2000), (4999,6000), (8999,10000)
        Fragment: blocks (996, 1050) + (9200, 9300)
        gstart=996, gend=9300
        tx_pos(996) = 996 - 999 = -3
        tx_pos(9300) = 1001 + 1001 + 301 = 2303
        FL = |2303 - (-3)| = 2306
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 996, 1050),
                _exon("chr1", 9200, 9300),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, three_exon_index)
        assert result is not None
        fl = _get_fl(result, three_exon_index, "t_three")
        assert fl == 2306

    def test_end_overhang_plus_intron_endpoint(self, three_exon_index):
        """Fragment gend extends past last exon AND gstart in intron.

        t_three exons: (999,2000), (4999,6000), (8999,10000)
        Fragment: blocks (2003, 2100) + (9950, 10007)
        gstart = 2003, gend = 10007
        tx_pos(2003) = 1001 + 3 = 1004 (3bp into intron after exon1)
        tx_pos(10007) = 1001 + 1001 + 1001 + 7 = 3010 (7bp past last exon)
        FL = 3010 - 1004 = 2006
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 2003, 2100),
                _exon("chr1", 9950, 10007),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, three_exon_index)
        assert result is not None
        fl = _get_fl(result, three_exon_index, "t_three")
        assert fl == 2006


# =====================================================================
# Edge cases
# =====================================================================


class TestEdgeCases:
    """Degenerate and boundary scenarios."""

    def test_fragment_exactly_at_exon_boundaries(self, basic_index):
        """Fragment spans exactly from exon1 start to exon2 end.

        t_basic exons: (999, 2000), (4999, 6000)
        Fragment: blocks (999, 1050) + (5950, 6000)
        gstart=999, gend=6000
        tx_pos(999) = 0 (start of exon1)
        tx_pos(6000) = 1001 + 1001 = 2002 (end of exon2 = t_len)
        FL = 2002
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 999, 1050),
                _exon("chr1", 5950, 6000),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 2002

    def test_fragment_one_bp_overhang_start(self, basic_index):
        """Fragment starts exactly 1bp before first exon.

        t_basic exons: (999, 2000), (4999, 6000)
        Fragment: blocks (998, 1050) + (5950, 6000)
        gstart=998, gend=6000
        tx_pos(998) = 998 - 999 = -1
        tx_pos(6000) = 2002
        FL = |2002 - (-1)| = 2003
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 998, 1050),
                _exon("chr1", 5950, 6000),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 2003

    def test_fragment_one_bp_overhang_end(self, basic_index):
        """Fragment ends exactly 1bp past last exon.

        t_basic exons: (999, 2000), (4999, 6000)
        Fragment: blocks (999, 1050) + (5950, 6001)
        gstart=999, gend=6001
        tx_pos(999) = 0
        tx_pos(6001) = 2002 + 1 = 2003
        FL = 2003
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 999, 1050),
                _exon("chr1", 5950, 6001),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 2003

    def test_fragment_one_bp_overhang_into_intron(self, basic_index):
        """Fragment gend extends exactly 1bp past exon1 into intron.

        t_basic exons: (999, 2000), (4999, 6000)
        Fragment: blocks (1900, 1950) + (1990, 2001)
        gstart=1900, gend=2001
        tx_pos(1900) = 901
        tx_pos(2001) = 1001 + 1 = 1002
        FL = 1002 - 901 = 101
        Genomic footprint = 2001 - 1900 = 101 (matches, as expected for single-exon case)
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1900, 1950),
                _exon("chr1", 1990, 2001),
            ),
            introns=(),
        )
        result = resolve_fragment(frag, basic_index)
        assert result is not None
        fl = _get_fl(result, basic_index, "t_basic")
        assert fl == 101
