"""Comprehensive tests for TranscriptIndex.load() correctness.

Verifies that every data structure produced by load() is correct after
the vectorized optimization (replacing pandas groupby/itertuples with
numpy sort + boundary detection).

Uses the shared mini_index fixture from conftest.py:
  g1 (+): t0 exons (99,200),(299,400),(499,600)  — 3 exons
          t1 exons (99,200),(499,600)              — 2 exons
  g2 (-): t2 exons (999,1100),(1199,1300)          — 2 exons
  genome: chr1 = 2000bp
"""

import numpy as np
import pytest

from rigel.index import TranscriptIndex
from rigel.types import GenomicInterval, IntervalType


# ═════════════════════════════════════════════════════════════════════
# Transcript & gene table integrity
# ═════════════════════════════════════════════════════════════════════


class TestTranscriptTable:
    """Verify t_df and derived arrays."""

    def test_num_transcripts(self, mini_index):
        assert mini_index.num_transcripts == 5  # 3 annotated + 2 synthetic nRNA

    def test_num_genes(self, mini_index):
        # Option A: synthetic nRNA transcripts get dedicated gene-neutral
        # rows in g_df. ``num_genes`` is the kernel-facing count (annotated
        # + synthetic), ``num_annotated_genes`` is user-facing.
        assert mini_index.num_annotated_genes == 2
        assert mini_index.num_genes == 4  # 2 annotated + 2 synthetic

    def test_t_index_matches_row_index(self, mini_index):
        t_df = mini_index.t_df
        assert (t_df.index == t_df["t_index"]).all()

    def test_transcript_ids(self, mini_index):
        t_ids = mini_index.t_df["t_id"].tolist()
        assert t_ids[:3] == ["t1", "t2", "t3"]
        # 2 synthetic nRNA transcripts appended
        assert len(t_ids) == 5
        assert all(tid.startswith("RIGEL_NRNA_") for tid in t_ids[3:])

    def test_gene_ids(self, mini_index):
        g_ids = mini_index.t_df["g_id"].tolist()
        assert g_ids[:3] == ["g1", "g1", "g2"]
        # Synthetic nRNA transcripts are gene-neutral: their g_id equals
        # their own t_id (``RIGEL_NRNA_...``).
        assert all(gid.startswith("RIGEL_NRNA_") for gid in g_ids[3:])
        assert g_ids[3:] == mini_index.t_df["t_id"].tolist()[3:]

    def test_t_to_g_arr(self, mini_index):
        # Synthetics now claim their own g_index (2 and 3); which maps to
        # which depends on merge order but each annotated gene still has
        # its two transcripts at indices 0, 1, 2.
        arr = mini_index.t_to_g_arr
        assert list(arr[:3]) == [0, 0, 1]
        assert sorted(arr[3:].tolist()) == [2, 3]

    def test_t_to_strand_arr(self, mini_index):
        arr = mini_index.t_to_strand_arr
        # t0, t1 on + (1); t2 on - (2)
        assert arr[0] == 1
        assert arr[1] == 1
        assert arr[2] == 2

    def test_g_to_strand_arr(self, mini_index):
        arr = mini_index.g_to_strand_arr
        assert arr[0] == 1  # g1 +
        assert arr[1] == 2  # g2 -


class TestGeneTable:
    """Verify derived gene table."""

    def test_gene_table_shape(self, mini_index):
        # 2 annotated + 2 synthetic gene-neutral rows
        assert len(mini_index.g_df) == 4
        assert "is_synthetic" in mini_index.g_df.columns
        assert int(mini_index.g_df["is_synthetic"].sum()) == 2

    def test_gene_g_index_matches(self, mini_index):
        g_df = mini_index.g_df
        assert (g_df.index == g_df["g_index"]).all()

    def test_gene_boundaries(self, mini_index):
        g_df = mini_index.g_df
        # g1: spans t0 and t1 — min start=99, max end=600
        g1 = g_df.loc[0]
        assert g1["start"] == 99
        assert g1["end"] == 600
        # g2: t2 only — start=999, end=1300
        g2 = g_df.loc[1]
        assert g2["start"] == 999
        assert g2["end"] == 1300

    def test_gene_num_transcripts(self, mini_index):
        g_df = mini_index.g_df
        # Annotated gene rows now contain only annotated transcripts;
        # synthetics get their own dedicated rows (num_transcripts == 1).
        assert g_df.loc[0, "num_transcripts"] == 2  # t1, t2
        assert g_df.loc[1, "num_transcripts"] == 1  # t3
        assert g_df.loc[2, "num_transcripts"] == 1  # synthetic
        assert g_df.loc[3, "num_transcripts"] == 1  # synthetic


# ═════════════════════════════════════════════════════════════════════
# Collapsed interval index (cgranges)
# ═════════════════════════════════════════════════════════════════════


class TestCollapsedIntervals:
    """Verify the cgranges overlap index + _iv_type + _iv_t_set."""

    def test_query_exon_t0_t1_shared(self, mini_index):
        """Query the first exon (99-200) shared by t0 and t1."""
        gi = GenomicInterval("chr1", 99, 200, 1)
        hits = mini_index.query(gi)
        exon_hits = [h for h in hits if h[2] == int(IntervalType.EXON)]
        # Should find at least one exon hit containing {0, 1}
        tsets = [h[3] for h in exon_hits]
        merged = frozenset().union(*tsets) if tsets else frozenset()
        assert 0 in merged
        assert 1 in merged

    def test_query_exon_t0_only(self, mini_index):
        """Query the middle exon (299-400) unique to t0."""
        gi = GenomicInterval("chr1", 299, 400, 1)
        hits = mini_index.query(gi)
        exon_hits = [h for h in hits if h[2] == int(IntervalType.EXON)]
        tsets = [h[3] for h in exon_hits]
        merged = frozenset().union(*tsets) if tsets else frozenset()
        assert 0 in merged
        assert 1 not in merged  # t1 doesn't have this exon

    def test_query_exon_t2_neg_strand(self, mini_index):
        """Query t2's first exon (999-1100) on neg strand."""
        gi = GenomicInterval("chr1", 999, 1100, 2)
        hits = mini_index.query(gi)
        exon_hits = [h for h in hits if h[2] == int(IntervalType.EXON)]
        tsets = [h[3] for h in exon_hits]
        merged = frozenset().union(*tsets) if tsets else frozenset()
        assert 2 in merged
        assert 0 not in merged
        assert 1 not in merged

    def test_query_intergenic(self, mini_index):
        """Query a region between gene clusters — should be intergenic."""
        gi = GenomicInterval("chr1", 700, 800, 0)
        hits = mini_index.query(gi)
        types = {h[2] for h in hits}
        # Should have an intergenic hit
        assert int(IntervalType.INTERGENIC) in types

    def test_transcript_span_hits(self, mini_index):
        """Query should include TRANSCRIPT-type intervals for full spans."""
        gi = GenomicInterval("chr1", 200, 350, 1)
        hits = mini_index.query(gi)
        tx_hits = [h for h in hits if h[2] == int(IntervalType.TRANSCRIPT)]
        # t0 spans 99-600, this query overlaps that
        tsets = [h[3] for h in tx_hits]
        merged = frozenset().union(*tsets) if tsets else frozenset()
        assert 0 in merged

    def test_iv_type_and_iv_t_set_lengths_match(self, mini_index):
        """_iv_type and _iv_t_set must have the same length."""
        assert len(mini_index._iv_type) == len(mini_index._iv_t_set)

    def test_iv_t_set_are_frozensets(self, mini_index):
        """All transcript sets should be frozensets."""
        for tset in mini_index._iv_t_set:
            assert isinstance(tset, frozenset)

    def test_all_iv_types_are_valid(self, mini_index):
        """All interval types should be valid IntervalType values."""
        valid = {int(e) for e in IntervalType}
        for itype in mini_index._iv_type:
            assert itype in valid, f"Unknown interval type: {itype}"


# ═════════════════════════════════════════════════════════════════════
# Per-transcript exon intervals
# ═════════════════════════════════════════════════════════════════════


class TestExonIntervals:
    """Verify _t_exon_intervals dict."""

    def test_all_transcripts_have_entries(self, mini_index):
        """Every transcript should have an exon interval entry."""
        for t_idx in range(mini_index.num_transcripts):
            assert t_idx in mini_index._t_exon_intervals, (
                f"Missing exon intervals for t_idx={t_idx}"
            )

    def test_t0_has_three_exons(self, mini_index):
        """t0 should have 3 exon intervals."""
        exons = mini_index._t_exon_intervals[0]
        assert exons.shape == (3, 2)
        assert exons.dtype == np.int32

    def test_t1_has_two_exons(self, mini_index):
        """t1 should have 2 exon intervals."""
        exons = mini_index._t_exon_intervals[1]
        assert exons.shape == (2, 2)

    def test_t2_has_two_exons(self, mini_index):
        """t2 should have 2 exon intervals."""
        exons = mini_index._t_exon_intervals[2]
        assert exons.shape == (2, 2)

    def test_t0_exon_coordinates(self, mini_index):
        """t0 exons should be (99,200), (299,400), (499,600)."""
        exons = mini_index._t_exon_intervals[0]
        expected = np.array([[99, 200], [299, 400], [499, 600]], dtype=np.int32)
        np.testing.assert_array_equal(exons, expected)

    def test_t1_exon_coordinates(self, mini_index):
        """t1 exons should be (99,200), (499,600)."""
        exons = mini_index._t_exon_intervals[1]
        expected = np.array([[99, 200], [499, 600]], dtype=np.int32)
        np.testing.assert_array_equal(exons, expected)

    def test_t2_exon_coordinates(self, mini_index):
        """t2 exons should be (999,1100), (1199,1300)."""
        exons = mini_index._t_exon_intervals[2]
        expected = np.array([[999, 1100], [1199, 1300]], dtype=np.int32)
        np.testing.assert_array_equal(exons, expected)

    def test_exons_sorted_by_start(self, mini_index):
        """Exon intervals for every transcript must be sorted by start coord."""
        for t_idx, exons in mini_index._t_exon_intervals.items():
            starts = exons[:, 0]
            assert np.all(starts[1:] >= starts[:-1]), (
                f"Exons not sorted for t_idx={t_idx}: {exons}"
            )

    def test_get_exon_intervals_method(self, mini_index):
        """The get_exon_intervals(t_idx) method should return same data."""
        for t_idx in range(mini_index.num_transcripts):
            result = mini_index.get_exon_intervals(t_idx)
            expected = mini_index._t_exon_intervals[t_idx]
            np.testing.assert_array_equal(result, expected)

    def test_get_exon_intervals_missing(self, mini_index):
        """get_exon_intervals for non-existent transcript returns None."""
        assert mini_index.get_exon_intervals(999) is None

    def test_no_negative_t_index(self, mini_index):
        """No negative t_index keys in the exon intervals dict."""
        for t_idx in mini_index._t_exon_intervals:
            assert t_idx >= 0


# ═════════════════════════════════════════════════════════════════════
# Splice junction data structures
# ═════════════════════════════════════════════════════════════════════


class TestSpliceJunctions:
    """Verify sj_map."""

    def test_sj_map_not_empty(self, mini_index):
        """mini_index has spliced transcripts so SJ map should be non-empty."""
        # t0: 2 introns (200-299, 400-499), t1: 1 intron (200-499), t2: 1 intron (1100-1199)
        assert len(mini_index.sj_map) > 0

    def test_sj_map_keys_are_tuples(self, mini_index):
        """SJ map keys should be (ref, start, end, strand) tuples."""
        for key in mini_index.sj_map:
            assert isinstance(key, tuple)
            assert len(key) == 4
            ref, start, end, strand = key
            assert isinstance(ref, str)
            assert isinstance(start, (int, np.integer))
            assert isinstance(end, (int, np.integer))
            assert isinstance(strand, (int, np.integer))

    def test_sj_map_values_are_frozensets(self, mini_index):
        """SJ map values should be frozensets of ints."""
        for val in mini_index.sj_map.values():
            assert isinstance(val, frozenset)
            for t_idx in val:
                assert isinstance(t_idx, (int, np.integer))
                assert t_idx >= 0

    def test_t0_introns_in_sj_map(self, mini_index):
        """t0's introns should be in the SJ map."""
        # t0 intron 1: 200-299 (between exons 99-200 and 299-400)
        # t0 intron 2: 400-499 (between exons 299-400 and 499-600)
        found_t0 = set()
        for (ref, start, end, strand), tset in mini_index.sj_map.items():
            if ref == "chr1" and 0 in tset:
                found_t0.add((start, end))
        assert len(found_t0) >= 2, f"Expected at least 2 SJs for t0, got {found_t0}"

    def test_t2_intron_in_sj_map(self, mini_index):
        """t2's intron (1100-1199) should be in the SJ map."""
        found = False
        for (ref, start, end, strand), tset in mini_index.sj_map.items():
            if ref == "chr1" and 2 in tset and start == 1100 and end == 1199:
                found = True
                assert strand == 2  # neg strand
                break
        assert found, "t2 intron 1100-1199 not found in sj_map"


# ═════════════════════════════════════════════════════════════════════
# C++ FragmentResolver
# ═════════════════════════════════════════════════════════════════════


class TestFragmentResolver:
    """Verify the native FragmentResolver is wired up."""

    def test_resolver_exists(self, mini_index):
        assert hasattr(mini_index, "resolver")
        assert mini_index.resolver is not None


# ═════════════════════════════════════════════════════════════════════
# Round-trip: build → load → query consistency
# ═════════════════════════════════════════════════════════════════════


class TestBuildLoadRoundTrip:
    """Build an index, load it, verify all structures agree."""

    @pytest.fixture
    def round_trip_index(self, tmp_path):
        """Build and load an index from GENCODE-style GTF."""
        import pysam
        from conftest import MINI_GTF

        gtf_path = tmp_path / "test.gtf"
        gtf_path.write_text(MINI_GTF)

        fasta_path = tmp_path / "genome.fa"
        with open(fasta_path, "w") as f:
            f.write(">chr1\n")
            seq = "N" * 2000
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")
        pysam.faidx(str(fasta_path))

        idx_dir = tmp_path / "idx"
        TranscriptIndex.build(fasta_path, gtf_path, idx_dir, write_tsv=False)
        return TranscriptIndex.load(idx_dir, retain_test_structures=True)

    def test_two_loads_produce_same_data(self, tmp_path):
        """Loading twice from the same dir should yield identical structures."""
        import pysam
        from conftest import MINI_GTF

        gtf_path = tmp_path / "test.gtf"
        gtf_path.write_text(MINI_GTF)
        fasta_path = tmp_path / "genome.fa"
        with open(fasta_path, "w") as f:
            f.write(">chr1\n")
            seq = "N" * 2000
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")
        pysam.faidx(str(fasta_path))

        idx_dir = tmp_path / "idx"
        TranscriptIndex.build(fasta_path, gtf_path, idx_dir, write_tsv=False)

        idx1 = TranscriptIndex.load(idx_dir, retain_test_structures=True)
        idx2 = TranscriptIndex.load(idx_dir, retain_test_structures=True)

        # Transcript tables
        assert idx1.t_df.equals(idx2.t_df)
        assert idx1.g_df.equals(idx2.g_df)
        np.testing.assert_array_equal(idx1.t_to_g_arr, idx2.t_to_g_arr)
        np.testing.assert_array_equal(idx1.t_to_strand_arr, idx2.t_to_strand_arr)
        np.testing.assert_array_equal(idx1.g_to_strand_arr, idx2.g_to_strand_arr)

        # Collapsed intervals
        assert len(idx1._iv_type) == len(idx2._iv_type)
        assert idx1._iv_type == idx2._iv_type
        assert len(idx1._iv_t_set) == len(idx2._iv_t_set)
        for s1, s2 in zip(idx1._iv_t_set, idx2._iv_t_set):
            assert s1 == s2

        # Exon intervals
        assert set(idx1._t_exon_intervals.keys()) == set(idx2._t_exon_intervals.keys())
        for t_idx in idx1._t_exon_intervals:
            np.testing.assert_array_equal(
                idx1._t_exon_intervals[t_idx],
                idx2._t_exon_intervals[t_idx],
            )

        # Splice junctions
        assert idx1.sj_map == idx2.sj_map

    def test_query_all_exons_round_trip(self, round_trip_index):
        """Query each exon from t_exon_intervals and confirm EXON hits."""
        idx = round_trip_index
        for t_idx, exons in idx._t_exon_intervals.items():
            for start, end in exons:
                gi = GenomicInterval("chr1", int(start), int(end), 0)
                hits = idx.query(gi)
                exon_hits = [h for h in hits if h[2] == int(IntervalType.EXON)]
                tsets = [h[3] for h in exon_hits]
                merged = frozenset().union(*tsets) if tsets else frozenset()
                assert t_idx in merged, (
                    f"t_idx={t_idx} not in exon query result for "
                    f"({start},{end}): {merged}"
                )



# ═════════════════════════════════════════════════════════════════════
# Stress test: larger synthetic index
# ═════════════════════════════════════════════════════════════════════


class TestLargerIndex:
    """Build a larger index with many overlapping transcripts and verify."""

    COMPLEX_GTF = ""

    @pytest.fixture(scope="class")
    def complex_index(self, tmp_path_factory):
        """Build an index with 10 genes, ~30 transcripts, overlapping exons."""
        from conftest import build_test_index

        # Build a more complex GTF with multiple genes and overlapping exons
        lines = []
        # Gene A (+ strand): 3 transcripts with varying exon structures
        for exon in [(100, 300), (500, 700), (900, 1100)]:
            lines.append(
                f'chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t+\t.\t'
                f'gene_id "gA"; transcript_id "tA1"; gene_name "GeneA"; '
                f'gene_type "protein_coding"; tag "basic";'
            )
        for exon in [(100, 300), (900, 1100)]:
            lines.append(
                f'chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t+\t.\t'
                f'gene_id "gA"; transcript_id "tA2"; gene_name "GeneA"; '
                f'gene_type "protein_coding";'
            )
        for exon in [(100, 300), (500, 700), (900, 1100), (1300, 1500)]:
            lines.append(
                f'chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t+\t.\t'
                f'gene_id "gA"; transcript_id "tA3"; gene_name "GeneA"; '
                f'gene_type "protein_coding";'
            )

        # Gene B (- strand): 2 transcripts
        for exon in [(2000, 2200), (2500, 2700)]:
            lines.append(
                f'chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t-\t.\t'
                f'gene_id "gB"; transcript_id "tB1"; gene_name "GeneB"; '
                f'gene_type "lncRNA";'
            )
        for exon in [(2000, 2200), (2500, 2700), (3000, 3200)]:
            lines.append(
                f'chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t-\t.\t'
                f'gene_id "gB"; transcript_id "tB2"; gene_name "GeneB"; '
                f'gene_type "lncRNA";'
            )

        # Gene C (+ strand): overlapping with Gene A on same strand
        for exon in [(800, 1100), (1300, 1600)]:
            lines.append(
                f'chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t+\t.\t'
                f'gene_id "gC"; transcript_id "tC1"; gene_name "GeneC"; '
                f'gene_type "protein_coding";'
            )

        # Gene D (+ strand): single-exon gene (unspliced)
        lines.append(
            'chr1\ttest\texon\t4000\t5000\t.\t+\t.\t'
            'gene_id "gD"; transcript_id "tD1"; gene_name "GeneD"; '
            'gene_type "protein_coding"; tag "basic";'
        )

        gtf_text = "\n".join(lines) + "\n"
        return build_test_index(
            tmp_path_factory, gtf_text, genome_size=6000, name="complex_idx"
        )

    def test_transcript_count(self, complex_index):
        # 7 annotated + 5 synthetic nRNA (gA:2, gB:2, gC:1)
        assert complex_index.num_transcripts == 12

    def test_gene_count(self, complex_index):
        # 4 annotated (gA, gB, gC, gD) + 5 synthetic gene-neutral rows
        assert complex_index.num_annotated_genes == 4
        assert complex_index.num_genes == 9

    def test_all_transcripts_have_exon_intervals(self, complex_index):
        for t_idx in range(complex_index.num_transcripts):
            assert t_idx in complex_index._t_exon_intervals

    def test_exon_intervals_shape_per_transcript(self, complex_index):
        t_df = complex_index.t_df
        expected_exon_counts = {
            "tA1": 3, "tA2": 2, "tA3": 4,
            "tB1": 2, "tB2": 3,
            "tC1": 2, "tD1": 1,
        }
        for _, row in t_df.iterrows():
            t_idx = row["t_index"]
            t_id = row["t_id"]
            exons = complex_index._t_exon_intervals[t_idx]
            if t_id.startswith("RIGEL_NRNA_"):
                # synthetic nRNA transcripts are single-exon
                assert exons.shape[0] == 1, (
                    f"{t_id} (t_idx={t_idx}): synthetic nRNA should have "
                    f"1 exon, got {exons.shape[0]}"
                )
            else:
                assert exons.shape[0] == expected_exon_counts[t_id], (
                    f"{t_id} (t_idx={t_idx}): expected {expected_exon_counts[t_id]} "
                    f"exons, got {exons.shape[0]}"
                )

    def test_unspliced_gene_has_no_sj(self, complex_index):
        """tD1 is single-exon → should have no splice junctions."""
        t_df = complex_index.t_df
        td1_idx = int(t_df.loc[t_df["t_id"] == "tD1", "t_index"].iloc[0])
        for (ref, start, end, strand), tset in complex_index.sj_map.items():
            assert td1_idx not in tset, (
                f"Single-exon transcript tD1 (t_idx={td1_idx}) should not "
                f"appear in SJ map, but found in ({ref},{start},{end},{strand})"
            )

    def test_overlapping_exons_both_found(self, complex_index):
        """Query exon region 900-1100 shared by gA and gC transcripts."""
        gi = GenomicInterval("chr1", 899, 1100, 0)
        hits = complex_index.query(gi)
        exon_hits = [h for h in hits if h[2] == int(IntervalType.EXON)]
        tsets = [h[3] for h in exon_hits]
        merged = frozenset().union(*tsets) if tsets else frozenset()
        # tA1, tA2, tA3 all have exons at 900-1100; tC1 has exon at 800-1100
        t_df = complex_index.t_df
        for _, row in t_df.iterrows():
            t_id = row["t_id"]
            t_idx = int(row["t_index"])
            if t_id in ("tA1", "tA2", "tA3", "tC1"):
                assert t_idx in merged, (
                    f"{t_id} (t_idx={t_idx}) should overlap query 899-1100"
                )

    def test_intergenic_between_gene_clusters(self, complex_index):
        """Region 3500-3900 should be intergenic."""
        gi = GenomicInterval("chr1", 3500, 3900, 0)
        hits = complex_index.query(gi)
        types = {h[2] for h in hits}
        assert int(IntervalType.INTERGENIC) in types

    def test_exon_interval_values_are_int32(self, complex_index):
        for t_idx, exons in complex_index._t_exon_intervals.items():
            assert exons.dtype == np.int32, (
                f"t_idx={t_idx}: expected int32, got {exons.dtype}"
            )


# ═════════════════════════════════════════════════════════════════════
# Edge cases
# ═════════════════════════════════════════════════════════════════════


class TestEdgeCases:
    """Edge cases: empty SJs, single-exon transcripts."""

    @pytest.fixture
    def single_exon_index(self, tmp_path_factory):
        """Index with only single-exon transcripts (no splice junctions)."""
        from conftest import build_test_index

        gtf = (
            'chr1\ttest\texon\t100\t500\t.\t+\t.\t'
            'gene_id "g1"; transcript_id "t1"; gene_name "G1"; '
            'gene_type "protein_coding"; tag "basic";\n'
            'chr1\ttest\texon\t1000\t1500\t.\t-\t.\t'
            'gene_id "g2"; transcript_id "t2"; gene_name "G2"; '
            'gene_type "lncRNA";\n'
        )
        return build_test_index(
            tmp_path_factory, gtf, genome_size=2000, name="single_exon_idx"
        )

    def test_single_exon_loads_ok(self, single_exon_index):
        assert single_exon_index.num_transcripts == 2

    def test_single_exon_empty_sj_map(self, single_exon_index):
        assert len(single_exon_index.sj_map) == 0

    def test_single_exon_exon_intervals(self, single_exon_index):
        for t_idx in range(2):
            exons = single_exon_index._t_exon_intervals[t_idx]
            assert exons.shape == (1, 2)

    def test_single_exon_query(self, single_exon_index):
        gi = GenomicInterval("chr1", 99, 500, 1)
        hits = single_exon_index.query(gi)
        exon_hits = [h for h in hits if h[2] == int(IntervalType.EXON)]
        assert len(exon_hits) > 0

    def test_single_exon_resolver_exists(self, single_exon_index):
        assert single_exon_index.resolver is not None
