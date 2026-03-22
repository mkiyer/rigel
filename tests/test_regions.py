"""Tests for the calibration region partition (flattened index).

Uses the shared mini_index fixture from conftest.py:

  chr1 (2000 bp):
    g1 (+): t0 exons (99,200),(299,400),(499,600)   span [99,600)
            t1 exons (99,200),(499,600)              span [99,600)
    g2 (-): t2 exons (999,1100),(1199,1300)          span [999,1300)

Note: mini_index uses build_test_index which creates only chr1
(no chr2 in the FASTA).

Expected region boundaries on chr1 (sorted unique):
  [0, 99, 200, 299, 400, 499, 600, 999, 1100, 1199, 1300, 2000]
  → 11 atomic regions
"""

import numpy as np
import pandas as pd
import pytest

from rigel.index import (
    TranscriptIndex,
    build_region_table,
    read_transcripts,
    load_reference_lengths,
    REGIONS_FEATHER,
)
from rigel.types import Strand


# ═════════════════════════════════════════════════════════════════════════════
# Helper: expected regions for the mini GTF
# ═════════════════════════════════════════════════════════════════════════════

# fmt: off
EXPECTED_CHR1 = [
    # (start, end, exon_pos, exon_neg, tx_pos, tx_neg)
    # Region partition built from annotated transcripts only
    # (synthetic nRNA transcripts are excluded).
    (   0,   99, False, False, False, False),  # intergenic
    (  99,  200, True,  False, True,  False),  # g1 exon 1 (t0+t1)
    ( 200,  299, False, False, True,  False),  # g1 intron
    ( 299,  400, True,  False, True,  False),  # g1 exon 2 (t0 only)
    ( 400,  499, False, False, True,  False),  # g1 intron
    ( 499,  600, True,  False, True,  False),  # g1 exon 3 (t0+t1)
    ( 600,  999, False, False, False, False),  # intergenic
    ( 999, 1100, False, True,  False, True ),  # g2 exon 1 (t2)
    (1100, 1199, False, False, False, True ),  # g2 intron
    (1199, 1300, False, True,  False, True ),  # g2 exon 2 (t2)
    (1300, 2000, False, False, False, False),  # intergenic
]
# fmt: on


# ═════════════════════════════════════════════════════════════════════════════
# Core partition tests (use mini_index which runs build + load)
# ═════════════════════════════════════════════════════════════════════════════


class TestRegionPartitionShape:
    """Basic shape and schema checks on the region table."""

    def test_region_df_exists(self, mini_index):
        assert mini_index.region_df is not None

    def test_region_cr_exists(self, mini_index):
        assert mini_index.region_cr is not None

    def test_total_region_count(self, mini_index):
        assert len(mini_index.region_df) == 11

    def test_columns(self, mini_index):
        expected_cols = {
            "region_id", "ref", "start", "end", "length",
            "exon_pos", "exon_neg", "tx_pos", "tx_neg",
        }
        assert set(mini_index.region_df.columns) == expected_cols

    def test_region_id_sequential(self, mini_index):
        ids = mini_index.region_df["region_id"].values
        np.testing.assert_array_equal(ids, np.arange(len(ids)))


class TestReferenceComplete:
    """Every reference is fully covered with no overlaps or gaps."""

    def test_chr1_total_length(self, mini_index):
        df = mini_index.region_df
        chr1 = df[df["ref"] == "chr1"]
        assert chr1["length"].sum() == 2000

    def test_only_chr1(self, mini_index):
        """build_test_index creates only chr1."""
        refs = mini_index.region_df["ref"].unique().tolist()
        assert refs == ["chr1"]

    def test_no_overlaps(self, mini_index):
        """Adjacent regions must be contiguous with no overlap."""
        df = mini_index.region_df
        for ref, grp in df.groupby("ref"):
            grp = grp.sort_values("start")
            starts = grp["start"].values
            ends = grp["end"].values
            # Each region start must equal the previous region end
            if len(starts) > 1:
                np.testing.assert_array_equal(starts[1:], ends[:-1])

    def test_starts_at_zero(self, mini_index):
        df = mini_index.region_df
        for ref, grp in df.groupby("ref"):
            assert grp["start"].min() == 0

    def test_length_column_consistent(self, mini_index):
        df = mini_index.region_df
        computed = df["end"] - df["start"]
        pd.testing.assert_series_equal(
            df["length"].reset_index(drop=True),
            computed.reset_index(drop=True),
            check_names=False,
        )


class TestFlagCorrectness:
    """Verify the four boolean flags for every region in the mini GTF."""

    def test_chr1_regions(self, mini_index):
        df = mini_index.region_df
        chr1 = df[df["ref"] == "chr1"].sort_values("start").reset_index(drop=True)
        assert len(chr1) == len(EXPECTED_CHR1)
        for i, (start, end, e_pos, e_neg, t_pos, t_neg) in enumerate(EXPECTED_CHR1):
            row = chr1.iloc[i]
            assert row["start"] == start, f"region {i}: start"
            assert row["end"] == end, f"region {i}: end"
            assert row["exon_pos"] == e_pos, f"region {i}: exon_pos"
            assert row["exon_neg"] == e_neg, f"region {i}: exon_neg"
            assert row["tx_pos"] == t_pos, f"region {i}: tx_pos"
            assert row["tx_neg"] == t_neg, f"region {i}: tx_neg"

    def test_intergenic_flags(self, mini_index):
        """Leading and trailing intergenic regions have all flags False."""
        df = mini_index.region_df
        first = df.iloc[0]
        assert first["start"] == 0
        assert not first["exon_pos"]
        assert not first["exon_neg"]
        assert not first["tx_pos"]
        assert not first["tx_neg"]


class TestFlagConsistency:
    """exon_* implies tx_* on the same strand."""

    def test_exon_pos_implies_tx_pos(self, mini_index):
        df = mini_index.region_df
        mask = df["exon_pos"]
        assert df.loc[mask, "tx_pos"].all()

    def test_exon_neg_implies_tx_neg(self, mini_index):
        df = mini_index.region_df
        mask = df["exon_neg"]
        assert df.loc[mask, "tx_neg"].all()


class TestDerivedAnnotations:
    """Verify derived context categories."""

    def test_intergenic_regions(self, mini_index):
        df = mini_index.region_df
        intergenic = df[
            ~df["exon_pos"] & ~df["exon_neg"] & ~df["tx_pos"] & ~df["tx_neg"]
        ]
        # chr1: [0,99), [600,999), [1300,2000)
        assert len(intergenic) == 3

    def test_intronic_pos_regions(self, mini_index):
        df = mini_index.region_df
        intronic_pos = df[df["tx_pos"] & ~df["exon_pos"]]
        # g1 has 2 intronic regions: [200,299) and [400,499)
        assert len(intronic_pos) == 2

    def test_intronic_neg_regions(self, mini_index):
        df = mini_index.region_df
        intronic_neg = df[df["tx_neg"] & ~df["exon_neg"]]
        # g2 has 1 intronic region: [1100,1199)
        assert len(intronic_neg) == 1


# ═════════════════════════════════════════════════════════════════════════════
# cgranges lookup tests
# ═════════════════════════════════════════════════════════════════════════════


class TestRegionCgranges:
    """Verify the cgranges index returns correct region_ids."""

    def test_intergenic_query(self, mini_index):
        """Query a point in the leading intergenic region."""
        hits = list(mini_index.region_cr.overlap("chr1", 10, 11))
        assert len(hits) == 1
        region_id = hits[0][2]
        row = mini_index.region_df[mini_index.region_df["region_id"] == region_id]
        assert row.iloc[0]["start"] == 0
        assert row.iloc[0]["end"] == 99

    def test_exon_query(self, mini_index):
        """Query a point inside the first exon."""
        hits = list(mini_index.region_cr.overlap("chr1", 150, 151))
        assert len(hits) == 1
        region_id = hits[0][2]
        row = mini_index.region_df[mini_index.region_df["region_id"] == region_id]
        assert row.iloc[0]["exon_pos"] is True or row.iloc[0]["exon_pos"] == True

    def test_boundary_query(self, mini_index):
        """Query that spans a region boundary returns two regions."""
        # spans [99,600) and [600,999)
        hits = list(mini_index.region_cr.overlap("chr1", 599, 601))
        assert len(hits) == 2

    def test_trailing_intergenic_query(self, mini_index):
        """Query a point in the trailing intergenic region."""
        hits = list(mini_index.region_cr.overlap("chr1", 1500, 1501))
        assert len(hits) == 1
        region_id = hits[0][2]
        row = mini_index.region_df[
            mini_index.region_df["region_id"] == region_id
        ]
        assert row.iloc[0]["start"] == 1300
        assert row.iloc[0]["end"] == 2000


# ═════════════════════════════════════════════════════════════════════════════
# Strand-ambiguous region test (separate fixture)
# ═════════════════════════════════════════════════════════════════════════════


class TestStrandAmbiguous:
    """When a pos and neg transcript overlap, the shared region is ambiguous."""

    @pytest.fixture
    def ambig_region_df(self, mini_gtf_file, mini_fasta_file, tmp_path):
        """Build a region table with overlapping pos/neg transcripts."""
        import textwrap
        # Override GTF: two transcripts on opposite strands that overlap
        ambig_gtf = textwrap.dedent("""\
            chr1\ttest\texon\t100\t300\t.\t+\t.\tgene_id "gA"; transcript_id "tA"; gene_name "A"; gene_type "pc";
            chr1\ttest\texon\t200\t400\t.\t-\t.\tgene_id "gB"; transcript_id "tB"; gene_name "B"; gene_type "pc";
        """)
        gtf_path = tmp_path / "ambig.gtf"
        gtf_path.write_text(ambig_gtf)

        transcripts = read_transcripts(gtf_path)
        ref_lengths = {"chr1": 1000}
        return build_region_table(transcripts, ref_lengths)

    def test_overlap_region_is_ambiguous(self, ambig_region_df):
        """The overlap [200,300) should have both pos and neg flags."""
        df = ambig_region_df
        overlap = df[(df["start"] == 199) & (df["end"] == 300)]
        assert len(overlap) == 1
        row = overlap.iloc[0]
        assert row["exon_pos"] == True
        assert row["exon_neg"] == True
        assert row["tx_pos"] == True
        assert row["tx_neg"] == True

    def test_pos_only_region(self, ambig_region_df):
        """[99,199) should be positive-strand only."""
        df = ambig_region_df
        region = df[(df["start"] == 99) & (df["end"] == 199)]
        assert len(region) == 1
        row = region.iloc[0]
        assert row["exon_pos"] == True
        assert row["exon_neg"] == False

    def test_neg_only_region(self, ambig_region_df):
        """[300,400) should be negative-strand only."""
        df = ambig_region_df
        region = df[(df["start"] == 300) & (df["end"] == 400)]
        assert len(region) == 1
        row = region.iloc[0]
        assert row["exon_pos"] == False
        assert row["exon_neg"] == True


# ═════════════════════════════════════════════════════════════════════════════
# Empty-reference test (build_region_table directly)
# ═════════════════════════════════════════════════════════════════════════════


class TestEmptyReference:
    """References with no transcripts produce a single intergenic region."""

    def test_empty_ref_single_region(self):
        ref_lengths = {"chr1": 1000}
        df = build_region_table([], ref_lengths)
        assert len(df) == 1
        row = df.iloc[0]
        assert row["ref"] == "chr1"
        assert row["start"] == 0
        assert row["end"] == 1000
        assert not row["exon_pos"]
        assert not row["tx_pos"]

    def test_multi_ref_with_empty(self, mini_gtf_file):
        """When a second ref has no transcripts, it gets one intergenic region."""
        transcripts = read_transcripts(mini_gtf_file)
        ref_lengths = {"chr1": 2000, "chr2": 500}
        df = build_region_table(transcripts, ref_lengths)
        chr2 = df[df["ref"] == "chr2"]
        assert len(chr2) == 1
        assert chr2.iloc[0]["start"] == 0
        assert chr2.iloc[0]["end"] == 500
        assert not chr2.iloc[0]["exon_pos"]


# ═════════════════════════════════════════════════════════════════════════════
# Multi-chromosome end-to-end test (full build + load)
# ═════════════════════════════════════════════════════════════════════════════

# GTF with genes on two chromosomes, third chromosome empty.
#   chr1 (+): tA — single exon 100-200
#   chr2 (-): tB — two exons 100-200, 300-400
#   chr3: no genes (800 bp, purely intergenic)
#
# 0-based half-open after GTF parse:
#   tA: exon (99,200), span [99,200)
#   tB: exons (99,200),(299,400), span [99,400)

import textwrap

MULTI_CHR_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "gA"; transcript_id "tA"; gene_name "A"; gene_type "pc";
    chr2\ttest\texon\t100\t200\t.\t-\t.\tgene_id "gB"; transcript_id "tB"; gene_name "B"; gene_type "pc";
    chr2\ttest\texon\t300\t400\t.\t-\t.\tgene_id "gB"; transcript_id "tB"; gene_name "B"; gene_type "pc";
""")


@pytest.fixture(scope="module")
def multi_chr_index(tmp_path_factory):
    """Build + load an index with genes on chr1 and chr2, chr3 empty."""
    import pysam

    base = tmp_path_factory.mktemp("multi_chr_idx")
    gtf_path = base / "test.gtf"
    gtf_path.write_text(MULTI_CHR_GTF)

    fasta_path = base / "genome.fa"
    with open(fasta_path, "w") as f:
        for name, size in [("chr1", 500), ("chr2", 600), ("chr3", 800)]:
            f.write(f">{name}\n")
            seq = "N" * size
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")
    pysam.faidx(str(fasta_path))

    idx_dir = base / "index"
    TranscriptIndex.build(fasta_path, gtf_path, idx_dir, write_tsv=True)
    return TranscriptIndex.load(idx_dir)


class TestMultiChromosome:
    """Full build+load end-to-end with genes on two chromosomes + one empty."""

    def test_all_three_refs_present(self, multi_chr_index):
        refs = sorted(multi_chr_index.region_df["ref"].unique().tolist())
        assert refs == ["chr1", "chr2", "chr3"]

    def test_region_id_globally_sequential(self, multi_chr_index):
        ids = multi_chr_index.region_df["region_id"].values
        np.testing.assert_array_equal(ids, np.arange(len(ids)))

    def test_reference_complete_all_refs(self, multi_chr_index):
        df = multi_chr_index.region_df
        expected = {"chr1": 500, "chr2": 600, "chr3": 800}
        for ref, expected_len in expected.items():
            total = df.loc[df["ref"] == ref, "length"].sum()
            assert total == expected_len, f"{ref}: {total} != {expected_len}"

    def test_no_overlaps_all_refs(self, multi_chr_index):
        df = multi_chr_index.region_df
        for ref, grp in df.groupby("ref"):
            grp = grp.sort_values("start")
            starts = grp["start"].values
            ends = grp["end"].values
            if len(starts) > 1:
                np.testing.assert_array_equal(starts[1:], ends[:-1])

    def test_chr1_regions(self, multi_chr_index):
        """chr1: [0,99) intergenic, [99,200) exon+, [200,500) intergenic."""
        df = multi_chr_index.region_df
        c1 = df[df["ref"] == "chr1"].sort_values("start").reset_index(drop=True)
        assert len(c1) == 3
        # Leading intergenic
        assert c1.iloc[0]["start"] == 0
        assert c1.iloc[0]["end"] == 99
        assert not c1.iloc[0]["tx_pos"]
        # Exon +
        assert c1.iloc[1]["start"] == 99
        assert c1.iloc[1]["end"] == 200
        assert c1.iloc[1]["exon_pos"] == True
        assert c1.iloc[1]["tx_pos"] == True
        assert c1.iloc[1]["exon_neg"] == False
        # Trailing intergenic
        assert c1.iloc[2]["start"] == 200
        assert c1.iloc[2]["end"] == 500

    def test_chr2_regions(self, multi_chr_index):
        """chr2: intergenic, exon-, intron-, exon-, intergenic."""
        df = multi_chr_index.region_df
        c2 = df[df["ref"] == "chr2"].sort_values("start").reset_index(drop=True)
        assert len(c2) == 5
        # [0,99) intergenic
        assert c2.iloc[0]["start"] == 0
        assert not c2.iloc[0]["tx_neg"]
        # [99,200) exon-
        assert c2.iloc[1]["start"] == 99
        assert c2.iloc[1]["end"] == 200
        assert c2.iloc[1]["exon_neg"] == True
        assert c2.iloc[1]["tx_neg"] == True
        # [200,299) intron-
        assert c2.iloc[2]["start"] == 200
        assert c2.iloc[2]["end"] == 299
        assert c2.iloc[2]["exon_neg"] == False
        assert c2.iloc[2]["tx_neg"] == True
        # [299,400) exon-
        assert c2.iloc[3]["start"] == 299
        assert c2.iloc[3]["end"] == 400
        assert c2.iloc[3]["exon_neg"] == True
        assert c2.iloc[3]["tx_neg"] == True
        # [400,600) intergenic
        assert c2.iloc[4]["start"] == 400
        assert c2.iloc[4]["end"] == 600

    def test_chr3_single_intergenic(self, multi_chr_index):
        """chr3 has no genes — should be one intergenic region."""
        df = multi_chr_index.region_df
        c3 = df[df["ref"] == "chr3"]
        assert len(c3) == 1
        row = c3.iloc[0]
        assert row["start"] == 0
        assert row["end"] == 800
        assert not row["exon_pos"]
        assert not row["exon_neg"]
        assert not row["tx_pos"]
        assert not row["tx_neg"]

    def test_cgranges_cross_ref_queries(self, multi_chr_index):
        """cgranges lookup returns correct regions across chromosomes."""
        cr = multi_chr_index.region_cr
        df = multi_chr_index.region_df

        # Query exon on chr1
        hits = list(cr.overlap("chr1", 150, 151))
        assert len(hits) == 1
        rid = hits[0][2]
        assert df.loc[df["region_id"] == rid, "exon_pos"].iloc[0] == True

        # Query in gene intron on chr2 — intronic (not exonic) region
        hits = list(cr.overlap("chr2", 250, 251))
        assert len(hits) == 1
        rid = hits[0][2]
        row = df[df["region_id"] == rid].iloc[0]
        assert row["tx_neg"] == True
        assert row["exon_neg"] == False

        # Query intergenic on chr3
        hits = list(cr.overlap("chr3", 400, 401))
        assert len(hits) == 1
        rid = hits[0][2]
        row = df[df["region_id"] == rid].iloc[0]
        assert not row["tx_pos"]
        assert not row["tx_neg"]

    def test_exon_implies_tx_all_refs(self, multi_chr_index):
        df = multi_chr_index.region_df
        assert df.loc[df["exon_pos"], "tx_pos"].all()
        assert df.loc[df["exon_neg"], "tx_neg"].all()


# ═════════════════════════════════════════════════════════════════════════════
# Adjacent-region merging test
# ═════════════════════════════════════════════════════════════════════════════


class TestAdjacentMerging:
    """Adjacent bins with identical flags are merged into one region."""

    def test_adjacent_exons_merged(self):
        """Two adjacent exons on the same strand from different transcripts
        should produce one merged exon region."""
        import textwrap
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "A"; gene_type "pc";
            chr1\ttest\texon\t200\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t2"; gene_name "A"; gene_type "pc";
        """)
        # After GTF parse: t1 exon (99,200), t2 exon (199,300)
        # Boundaries: {0, 99, 199, 200, 300, 1000}
        # Bins: [0,99) [99,199) [199,200) [200,300) [300,1000)
        # Flags: [99,199) = exon+/tx+, [199,200) = exon+/tx+, [200,300) = exon+/tx+
        # These three should merge into one [99,300) exon+ region
        from pathlib import Path
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            gtf_path = Path(td) / "test.gtf"
            gtf_path.write_text(gtf)
            transcripts = read_transcripts(gtf_path)
        ref_lengths = {"chr1": 1000}
        df = build_region_table(transcripts, ref_lengths)
        # Should have 3 regions: [0,99) intergenic, [99,300) exon+, [300,1000) intergenic
        assert len(df) == 3
        mid = df.iloc[1]
        assert mid["start"] == 99
        assert mid["end"] == 300
        assert mid["exon_pos"] == True
        assert mid["tx_pos"] == True

    def test_non_adjacent_not_merged(self):
        """Non-adjacent regions with same flags stay separate."""
        import textwrap
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "A"; gene_type "pc";
            chr1\ttest\texon\t400\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "A"; gene_type "pc";
        """)
        # t1: exons (99,200) and (399,500), tx span [99,500)
        # The intronic region [200,399) has tx_pos but NOT exon_pos → different flags
        # So exon regions can't merge across it
        from pathlib import Path
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            gtf_path = Path(td) / "test.gtf"
            gtf_path.write_text(gtf)
            transcripts = read_transcripts(gtf_path)
        ref_lengths = {"chr1": 1000}
        df = build_region_table(transcripts, ref_lengths)
        exon_regions = df[df["exon_pos"]]
        assert len(exon_regions) == 2  # two separate exon regions

    def test_merge_reduces_count(self):
        """Merging should reduce region count when there are adjacent
        identically-flagged bins."""
        import textwrap
        # Three transcripts all with single exons that tile [100,400)
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "A"; gene_type "pc";
            chr1\ttest\texon\t200\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t2"; gene_name "A"; gene_type "pc";
            chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t3"; gene_name "A"; gene_type "pc";
        """)
        from pathlib import Path
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            gtf_path = Path(td) / "test.gtf"
            gtf_path.write_text(gtf)
            transcripts = read_transcripts(gtf_path)
        ref_lengths = {"chr1": 1000}
        df = build_region_table(transcripts, ref_lengths)
        # Should have: [0,99) intergenic, [99,400) exon+, [400,1000) intergenic
        assert len(df) == 3

    def test_reference_complete_after_merge(self):
        """Total length is preserved after merging."""
        import textwrap
        gtf = textwrap.dedent("""\
            chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "A"; gene_type "pc";
            chr1\ttest\texon\t200\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t2"; gene_name "A"; gene_type "pc";
        """)
        from pathlib import Path
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            gtf_path = Path(td) / "test.gtf"
            gtf_path.write_text(gtf)
            transcripts = read_transcripts(gtf_path)
        ref_lengths = {"chr1": 1000}
        df = build_region_table(transcripts, ref_lengths)
        assert df["length"].sum() == 1000
