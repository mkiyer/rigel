"""A fragment's length on each candidate transcript is the length of the path the candidate implies: its
genomic span less the junctions the reads observed and less the candidate's own introns lying in an
unsequenced gap, within K of the gap's ends. An intron the reads cover is the molecule's own sequence and is
never cut. It is the length calibration deposits for the same path (``Accumulator.length_under``).

Spliced and unspliced fragments; overhang past either end of the transcript and into an intron; both ends
overhanging at once; a large unsequenced mate gap spanning one, two or three introns, so the cut has to
remove every intron in the gap rather than the first; nascent single-exon candidates, whose length is the
genomic span; negative-strand transcripts, where the length is strand-agnostic; and several candidates of
different exon structure, which imply DIFFERENT lengths for the same fragment. That last point is why the
length is per candidate and not per fragment. ``TestOneLengthRule`` reads one molecule three ways and gets
one length. ``TestLengthRuleByEnumeration`` holds every candidate of thousands of fragments, whose mates
stop short of introns, overrun them and run through them, to a base-by-base count and to the accumulator's
own length.
"""

import contextlib
import textwrap

import numpy as np
import pytest

from rigel.native import Accumulator as NativeAccumulator
from rigel.types import Strand, GenomicInterval
from _resolution_reference import make_fragment, resolve_fragment
from _index_builder import build_test_index
from native._accumulator_reference import GapHypothesis


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

# Introns shorter than a mate, so mates run into them and through them
# 0-based half-open after parse:
#   g1: t_a exons (100,200), (260,400), (430,600) — introns (200,260), (400,430)
#       t_b exons (100,200), (400,600)             — intron (200,400) holds t_a's middle exon
#       t_c exons (150,230), (260,400)             — its first exon runs into t_a's first intron
#   g2: T2  exons (1000,1200), (1260,1500)         — intron (1200,1260)
SHORT_INTRON_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t_a"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t261\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t_a"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t431\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t_a"; gene_name "G1"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t_b"; gene_name "G1"; gene_type "protein_coding";
    chr1\ttest\texon\t401\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t_b"; gene_name "G1"; gene_type "protein_coding";
    chr1\ttest\texon\t151\t230\t.\t+\t.\tgene_id "g1"; transcript_id "t_c"; gene_name "G1"; gene_type "protein_coding";
    chr1\ttest\texon\t261\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t_c"; gene_name "G1"; gene_type "protein_coding";
    chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "g2"; transcript_id "T2"; gene_name "G2"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t1261\t1500\t.\t+\t.\tgene_id "g2"; transcript_id "T2"; gene_name "G2"; gene_type "protein_coding"; tag "basic";
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


@pytest.fixture(scope="session")
def short_intron_index(tmp_path_factory):
    return build_test_index(
        tmp_path_factory, SHORT_INTRON_GTF, genome_size=2000, name="fl_short_intron"
    )


@contextlib.contextmanager
def _tuned(index, K, max_frag=0):
    """The index's resolver with intron-overhang tolerance ``K`` and maximum fragment length ``max_frag``,
    restored to the session fixtures' defaults (both 0) on exit."""
    resolver = index.resolver
    resolver.set_intron_overhang_tolerance(K)
    resolver.set_max_fragment_length(max_frag)
    try:
        yield
    finally:
        resolver.set_intron_overhang_tolerance(0)
        resolver.set_max_fragment_length(0)


# =====================================================================
# Helper to extract FL for a specific transcript from resolve result
# =====================================================================


def _get_fl(result, index, t_id):
    """Get fragment length for transcript t_id from a ResolvedFragment."""
    tm = _t_map(index)
    t_idx = tm[t_id]
    fl_dict = result.frag_lengths
    return fl_dict.get(t_idx)


def _nascent(index, start, end):
    """The t_index of the nascent span over ``[start, end)``: the synthetic single-exon transcript."""
    t = index.t_df
    hit = t[t["is_synthetic"] & (t["start"] == start) & (t["end"] == end)]
    assert len(hit) == 1, hit
    return int(hit["t_index"].iloc[0])


# =====================================================================
# One molecule, one length
# =====================================================================


class TestOneLengthRule:
    """One 170-bp molecule, 1150 to 1320 across T2's 60-bp intron (1200, 1260), read three ways.

    An intron the reads cover is the molecule's own sequence: whether the mates overlap (one block) or each
    runs into the intron (25 and 15 bases, more than K = 3), T2's length is the span, 170. Only an intron in
    the unsequenced gap is cut, and only by the candidate that has it: T2 110, the nascent span 170. The
    endpoint rule cut every intron of the candidate between the fragment's ends, and read the mates that run
    into the intron as 110 on T2.
    """

    @pytest.mark.parametrize(
        "blocks, on_t2, on_nascent",
        [
            pytest.param(((1150, 1320),), 170, 170, id="mates-overlap"),
            pytest.param(((1150, 1225), (1245, 1320)), 170, 170, id="mates-run-into-the-intron"),
            pytest.param(((1150, 1200), (1260, 1320)), 110, 170, id="the-intron-in-the-gap"),
        ],
    )
    def test_one_molecule_has_one_length(self, short_intron_index, blocks, on_t2, on_nascent):
        frag = make_fragment(exons=tuple(_exon("chr1", s, e) for s, e in blocks))
        with _tuned(short_intron_index, K=3, max_frag=1000):
            result = resolve_fragment(frag, short_intron_index)
        assert result is not None
        assert _get_fl(result, short_intron_index, "T2") == on_t2
        assert result.frag_lengths.get(_nascent(short_intron_index, 1000, 1500)) == on_nascent


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
    """An end overhanging into an intron counts; an intron the mates run into is cut only within K."""

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

    def test_an_intron_the_mates_run_into_is_kept(self, basic_index):
        """The mates run 5 and 4 bases into t_basic's intron (2000, 4999), more than K = 3.

        The reads show those bases, so on t_basic the intron is the molecule's own sequence and is not cut: the
        length is the span, 5149 - 1855 = 3294, and the EM prices the covered bases as overhang. The endpoint
        rule cut the intron and gave 295, a plausible length for a molecule t_basic cannot have made.
        """
        frag = make_fragment(
            exons=(
                _exon("chr1", 1855, 2005),
                _exon("chr1", 4995, 5149),
            ),
            introns=(),
        )
        with _tuned(basic_index, K=3):
            result = resolve_fragment(frag, basic_index)
        assert result is not None
        assert _get_fl(result, basic_index, "t_basic") == 3294

    def test_an_overrun_within_K_is_cut(self, basic_index):
        """The mates run 3 and 2 bases into t_basic's intron, within K = 3: the aligner ran past the donor and
        the acceptor, so the intron lies in the gap and is cut whole, the overrun bases with it.
        (5151 - 1853) - 2999 = 299."""
        frag = make_fragment(
            exons=(
                _exon("chr1", 1853, 2003),
                _exon("chr1", 4997, 5151),
            ),
            introns=(),
        )
        with _tuned(basic_index, K=3):
            result = resolve_fragment(frag, basic_index)
        assert result is not None
        assert _get_fl(result, basic_index, "t_basic") == 299

    def test_gstart_in_intron(self, basic_index):
        """Fragment gstart is in the intron between exons.

        t_basic exons: (999, 2000), (4999, 6000), intron (2000, 4999)
        Fragment: blocks (2003, 2100) + (5200, 5400)
        gstart = 2003 is 3 bp into the intron, so the first read lies wholly in it and the fragment
        jumps no intron: under t_basic the molecule is 4999 - 2003 = 2996 bp of overhang before exon2's
        start plus 401 exonic bases, FL = 3397 = 5400 - 2003, its genomic span. Measured forward from
        exon1's END instead, the start lands one intron length downstream and reads |3397 - 2999| = 398.
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
        assert fl == 3397


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
        # Find the nRNA candidate (synthetic nRNA has "is_synthetic" flag)
        t_df = basic_index.t_df
        nrna_mask = t_df["is_synthetic"].values if "is_synthetic" in t_df.columns else None
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

        t_three exons: (999,2000), (4999,6000), (8999,10000); introns (2000,4999), (6000,8999)
        Fragment: blocks (2003, 2100) + (9950, 10007)
        gstart = 2003 is in the first intron, which the fragment does not jump; the second lies wholly
        inside the span and is jumped.
        FL = (10007 - 2003) - 2999 = 5005: 2996 bp of start overhang + 1001 + 1001 exonic + 7 bp past
        the last exon.
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
        assert fl == 5005


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


# =====================================================================
# The rule, against a base-by-base count
# =====================================================================

#: The indexes the enumeration draws fragments on: short introns that mates run through, and long ones that a
#: mate gap spans.
_ENUMERATED = ("short_intron_index", "three_exon_index", "isoform_index")
#: Fragments drawn per index.
_N_FRAGMENTS = 2000
#: The longest mate drawn: longer than every intron of ``short_intron_index``, so mates run through introns.
_MAX_MATE = 250
#: How far from an exon boundary a mate's end is drawn: past K = 3, so an overrun lies within K or beyond it.
_REACH = 5
#: The bases a spliced mate aligns on each side of its junction.
_ANCHORS = (1, 2, 3, 4, 5, 20, 75)
#: The resolver's maximum fragment length while enumerating, so that a fragment joining the two genes of
#: ``short_intron_index`` is an ordinary molecule rather than a chimera.
_MAX_FRAG = 1000


def _all_exons(index):
    """Every transcript's sorted ``[start, end)`` exons, by t_index, nascent spans included."""
    return {int(t): np.asarray(index.get_exon_intervals(int(t))) for t in index.t_df["t_index"]}


def _introns(exons):
    """A transcript's introns as ``(start, end)`` pairs, from its sorted ``[start, end)`` exons."""
    return [(int(a), int(b)) for a, b in zip(exons[:-1, 1], exons[1:, 0])]


def _points(index):
    """Every position within ``_REACH`` of an exon boundary of any transcript, nascent spans included, and
    each exon's midpoint."""
    points = set()
    for t in index.t_df["t_index"]:
        for s, e in index.get_exon_intervals(int(t)):
            for b in (int(s), int(e)):
                points.update(range(b - _REACH, b + _REACH + 1))
            points.add((int(s) + int(e)) // 2)
    return np.asarray(sorted(points))


def _mate(rng, points):
    """One mate's aligned interval: two drawn positions at most ``_MAX_MATE`` apart."""
    while True:
        start = int(rng.choice(points))
        ends = points[(points > start) & (points <= start + _MAX_MATE)]
        if ends.size:
            return start, int(rng.choice(ends))


def _merge(intervals):
    """Sorted ``(start, end)`` blocks, overlapping and abutting ones merged, as the scanner assembles mates."""
    blocks = []
    for s, e in sorted(intervals):
        if blocks and s <= blocks[-1][1]:
            blocks[-1] = (blocks[-1][0], max(blocks[-1][1], e))
        else:
            blocks.append((s, e))
    return blocks


def _fragments(index, seed):
    """``_N_FRAGMENTS`` fragments on ``index`` as ``(blocks, observed junctions)``: two mates whose ends are
    drawn from `_points`, so they stop short of introns, overrun them by every distance up to and past K, and
    run through them whole. One in three has its first mate spliced across an annotated junction, each anchor
    drawn from ``_ANCHORS``, and its second mate clear of that junction."""
    rng = np.random.default_rng(seed)
    points = _points(index)
    junctions = sorted({j for exons in _all_exons(index).values() for j in _introns(exons)})
    for _ in range(_N_FRAGMENTS):
        if rng.random() < 1 / 3:
            a, b = junctions[int(rng.integers(len(junctions)))]
            left, right = (int(x) for x in rng.choice(_ANCHORS, 2))
            mate = _mate(rng, points)
            while mate[1] > a and mate[0] < b:
                mate = _mate(rng, points)
            yield _merge([(a - left, a), (b, b + right), mate]), ((a, b),)
        else:
            yield _merge([_mate(rng, points), _mate(rng, points)]), ()


def _resolve_blocks(index, blocks, observed):
    """The fragment with these blocks and observed junctions, resolved; ``None`` when it has no lengths (no
    result, or a chimera)."""
    frag = make_fragment(
        exons=tuple(_exon("chr1", s, e) for s, e in blocks),
        introns=tuple(GenomicInterval("chr1", a, b, Strand.POS) for a, b in observed),
    )
    result = resolve_fragment(frag, index)
    return None if result is None or result.chimera_type else result


def _gaps(blocks, observed):
    """The unsequenced gaps: the holes between consecutive blocks that are not an observed junction."""
    return [(x[1], y[0]) for x, y in zip(blocks, blocks[1:]) if (x[1], y[0]) not in observed]


def _implied(exons, gaps, K):
    """The candidate's introns that overlap an unsequenced gap and lie within ``K`` of its ends: what its path
    cuts beyond the observed junctions."""
    return [
        (a, b)
        for a, b in _introns(exons)
        if b > a and any(a < g1 and b > g0 and a >= g0 - K and b <= g1 + K for g0, g1 in gaps)
    ]


def _path_length(blocks, observed, implied):
    """The bases of the fragment's span, first block's start to last block's end, outside the observed
    junctions and the implied introns, counted one by one."""
    pos = np.arange(blocks[0][0], blocks[-1][1])
    kept = np.ones(pos.size, dtype=bool)
    for a, b in (*observed, *implied):
        kept &= (pos < a) | (pos >= b)
    return int(kept.sum())


def _overlap(x, y):
    """The bases two ``(start, end)`` intervals share; zero or less when they share none."""
    return min(x[1], y[1]) - max(x[0], y[0])


class TestLengthRuleByEnumeration:
    """Every candidate's length, for thousands of drawn fragments on three indexes, is the base-by-base count
    of the path it implies (`_path_length`), and is the accumulator's length for that path.

    The mates' ends sit at and beside exon boundaries, so the mates stop short of introns, overrun them within
    K and past it, and run through them whole; a third of the fragments carry an observed junction. The
    endpoint rule this replaced cut every intron of the candidate between the fragment's two ends, so it
    shortened a candidate whose intron the reads cover: it fails here wherever a mate runs into or through an
    intron by more than K.
    """

    @pytest.mark.parametrize("K", [0, 3])
    def test_every_candidate_length_is_its_path_length(self, K, request):
        n_checked = n_spliced = n_covered = n_within_k = 0
        wrong = []
        for name in _ENUMERATED:
            index = request.getfixturevalue(name)
            exons = _all_exons(index)
            with _tuned(index, K, _MAX_FRAG):
                for blocks, observed in _fragments(index, seed=0):
                    result = _resolve_blocks(index, blocks, observed)
                    if result is None:
                        continue
                    gaps = _gaps(blocks, observed)
                    for t in result.t_inds:
                        ex = exons[int(t)]
                        implied = _implied(ex, gaps, K)
                        want = _path_length(blocks, observed, implied)
                        got = result.frag_lengths.get(
                            int(t), 0
                        )  # a missing length is absent from the dict
                        n_checked += 1
                        n_spliced += bool(observed)
                        n_covered += any(
                            _overlap(i, x) > K
                            for i in _introns(ex)
                            if i not in observed
                            for x in blocks
                        )
                        n_within_k += any(_overlap(i, x) > 0 for i in implied for x in blocks)
                        if got != want:
                            wrong.append((name, blocks, observed, int(t), got, want))
        # ⚠ The enumeration must reach each case it exists for, or a pass means nothing: observed junctions,
        # introns the reads cover by more than K, and (K > 0) introns cut although a mate runs into them.
        assert n_spliced > 100, n_spliced
        assert n_covered > 100, n_covered
        assert n_within_k > 100 or K == 0, n_within_k
        assert not wrong, (
            f"{len(wrong)} of {n_checked} lengths wrong (index, blocks, observed, t, got, want): {wrong[:4]}"
        )

    def test_every_candidate_length_is_the_deposits_length(self, request):
        """The resolver's length on a candidate is ``Accumulator.length_under`` for the path the candidate
        implies: one formula for L in C++, which the EM and calibration share. The obvious second formula, the
        span less the summed introns, disagrees wherever a cut intron runs past the fragment's end."""
        K = 3
        n_checked = n_clipped = 0
        wrong = []
        for name in _ENUMERATED:
            index = request.getfixturevalue(name)
            exons = _all_exons(index)
            accumulator = NativeAccumulator(
                region_bounds=np.asarray([0, 1 << 30], dtype=np.int64),
                region_types=np.zeros(1, dtype=np.uint8),
                max_length=_MAX_FRAG,
                ref=0,
            )
            with _tuned(index, K, _MAX_FRAG):
                for blocks, observed in _fragments(index, seed=1):
                    result = _resolve_blocks(index, blocks, observed)
                    if result is None:
                        continue
                    start, end = blocks[0][0], blocks[-1][1]
                    gaps = _gaps(blocks, observed)
                    for t in result.t_inds:
                        implied = _implied(exons[int(t)], gaps, K)
                        want = accumulator.length_under(
                            start=start,
                            end=end,
                            observed_introns=observed,
                            hypotheses=(GapHypothesis(tuple(implied)),),
                        )
                        got = result.frag_lengths.get(int(t), 0)
                        n_checked += 1
                        n_clipped += any(a < start or b > end for a, b in implied)
                        if got != want:
                            wrong.append((name, blocks, observed, int(t), got, want))
        assert n_clipped > 10, n_clipped  # ⚠ the cases a second formula gets wrong are reached
        assert not wrong, (
            f"{len(wrong)} of {n_checked} lengths wrong (index, blocks, observed, t, got, want): {wrong[:4]}"
        )
