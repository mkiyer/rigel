"""`rigel.index` end to end: what `TranscriptIndex.load()` produces, what it refuses, the duplicate-transcript
guard the reader applies, and the layout iterator underneath the interval and region tables.

The first block verifies every structure `load()` builds (the transcript and gene tables and their derived arrays,
exon intervals, splice junctions, the fragment resolver) and the build/load round trip. The second gates what
`build()` and `load()` refuse: a reordered transcript table or a gapped gene index, an index missing a file `load()`
reads, a setting `IndexConfig` refuses and a GTF with no transcripts (both refused before an existing index is
touched), and an index whose rebuild stopped part-way. The third gates the duplicate-exon-structure guard in
`read_transcripts` and its collapse flag; the fourth `iter_reference_layout`, which both `intervals.feather` and
`regions.feather` are derived from.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rigel.index import (
    EDGES_FEATHER,
    INTERVALS_FEATHER,
    REF_LENGTHS_FEATHER,
    REGIONS_FEATHER,
    SJ_FEATHER,
    TRANSCRIPT_ALIASES_FEATHER,
    TRANSCRIPTS_FEATHER,
    IndexLoadError,
    TranscriptIndex,
    GenicSpan,
    IntergenicSpan,
    _sj_csr,
    iter_reference_layout,
    read_transcripts,
)
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _mini_sources(base: Path) -> tuple[Path, Path]:
    """Write ``MINI_GTF`` and an indexed 2,000-bp chr1 of known sequence into *base*; return (FASTA, GTF)."""
    from _index_builder import MINI_GTF, write_fasta

    fasta = write_fasta(base / "genome.fa", {"chr1": 2000})
    gtf = base / "test.gtf"
    gtf.write_text(MINI_GTF)
    return fasta, gtf


# The first block drives `conftest.py`'s shared `mini_index` fixture, whose geometry the numbers
# below are read off:
#   g1 (+): t0 exons (99,200),(299,400),(499,600)  — 3 exons
#           t1 exons (99,200),(499,600)              — 2 exons
#   g2 (-): t2 exons (999,1100),(1199,1300)          — 2 exons
#   genome: chr1 = 2000bp


# ═════════════════════════════════════════════════════════════════════
# Transcript & gene table integrity
# ═════════════════════════════════════════════════════════════════════


class TestTranscriptTable:
    """Verify t_df and derived arrays."""

    def test_num_transcripts(self, mini_index):
        assert mini_index.num_transcripts == 5  # 3 annotated + 2 synthetic spans

    def test_num_genes(self, mini_index):
        # Synthetic spans get dedicated gene-neutral rows in g_df.
        # ``num_genes`` is the kernel-facing count (annotated + synthetic),
        # ``num_annotated_genes`` is user-facing.
        assert mini_index.num_annotated_genes == 2
        assert mini_index.num_genes == 4  # 2 annotated + 2 synthetic

    def test_t_index_matches_row_index(self, mini_index):
        t_df = mini_index.t_df
        assert (t_df.index == t_df["t_index"]).all()

    def test_transcript_ids(self, mini_index):
        t_ids = mini_index.t_df["t_id"].tolist()
        assert t_ids[:3] == ["t1", "t2", "t3"]
        # 2 synthetic spans appended
        assert len(t_ids) == 5
        assert all(tid.startswith("RIGEL_NASCENT_") for tid in t_ids[3:])

    def test_gene_ids(self, mini_index):
        g_ids = mini_index.t_df["g_id"].tolist()
        assert g_ids[:3] == ["g1", "g1", "g2"]
        # Synthetic spans are gene-neutral: their g_id equals
        # their own t_id (``RIGEL_NASCENT_...``).
        assert all(gid.startswith("RIGEL_NASCENT_") for gid in g_ids[3:])
        assert g_ids[3:] == mini_index.t_df["t_id"].tolist()[3:]

    def test_t_to_g_arr(self, mini_index):
        # Synthetics claim their own g_index (2 and 3); which maps to
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
        # Annotated gene rows contain only annotated transcripts;
        # synthetics get their own dedicated rows (num_transcripts == 1).
        assert g_df.loc[0, "num_transcripts"] == 2  # t1, t2
        assert g_df.loc[1, "num_transcripts"] == 1  # t3
        assert g_df.loc[2, "num_transcripts"] == 1  # synthetic
        assert g_df.loc[3, "num_transcripts"] == 1  # synthetic


# ═════════════════════════════════════════════════════════════════════
# Per-transcript exon intervals
# ═════════════════════════════════════════════════════════════════════


class TestExonIntervals:
    """Verify each transcript's exons, its slice of the exon CSR (`get_exon_intervals`)."""

    def test_all_transcripts_have_exons(self, mini_index):
        """Every transcript, a synthetic span too, has at least one exon."""
        for t_idx in range(mini_index.num_transcripts):
            assert len(mini_index.get_exon_intervals(t_idx)) > 0, f"No exons for t_idx={t_idx}"

    def test_t0_has_three_exons(self, mini_index):
        """t0 should have 3 exon intervals."""
        exons = mini_index.get_exon_intervals(0)
        assert exons.shape == (3, 2)
        assert exons.dtype == np.int32

    def test_t1_has_two_exons(self, mini_index):
        """t1 should have 2 exon intervals."""
        exons = mini_index.get_exon_intervals(1)
        assert exons.shape == (2, 2)

    def test_t2_has_two_exons(self, mini_index):
        """t2 should have 2 exon intervals."""
        exons = mini_index.get_exon_intervals(2)
        assert exons.shape == (2, 2)

    def test_t0_exon_coordinates(self, mini_index):
        """t0 exons should be (99,200), (299,400), (499,600)."""
        exons = mini_index.get_exon_intervals(0)
        expected = np.array([[99, 200], [299, 400], [499, 600]], dtype=np.int32)
        np.testing.assert_array_equal(exons, expected)

    def test_t1_exon_coordinates(self, mini_index):
        """t1 exons should be (99,200), (499,600)."""
        exons = mini_index.get_exon_intervals(1)
        expected = np.array([[99, 200], [499, 600]], dtype=np.int32)
        np.testing.assert_array_equal(exons, expected)

    def test_t2_exon_coordinates(self, mini_index):
        """t2 exons should be (999,1100), (1199,1300)."""
        exons = mini_index.get_exon_intervals(2)
        expected = np.array([[999, 1100], [1199, 1300]], dtype=np.int32)
        np.testing.assert_array_equal(exons, expected)

    def test_exons_sorted_by_start(self, mini_index):
        """Exon intervals for every transcript must be sorted by start coord."""
        for t_idx in range(mini_index.num_transcripts):
            exons = mini_index.get_exon_intervals(t_idx)
            starts = exons[:, 0]
            assert np.all(starts[1:] >= starts[:-1]), f"Exons not sorted for t_idx={t_idx}: {exons}"

    def test_get_exon_intervals_refuses_outside_the_table(self, mini_index):
        """A transcript index outside the table raises; -1 would otherwise read an empty slice as a transcript."""
        for t_idx in (mini_index.num_transcripts, 999, -1):
            with pytest.raises(IndexError):
                mini_index.get_exon_intervals(t_idx)


# ═════════════════════════════════════════════════════════════════════
# Splice junction data structures
# ═════════════════════════════════════════════════════════════════════


class TestSpliceJunctions:
    """Verify the annotated junctions the index installs in the resolver (`_sj_map`)."""

    def test_sj_map_not_empty(self, mini_index):
        """mini_index has spliced transcripts so SJ map should be non-empty."""
        # t0: 2 introns (200-299, 400-499), t1: 1 intron (200-499), t2: 1 intron (1100-1199)
        assert len(_sj_map(mini_index)) > 0

    def test_sj_map_keys_are_tuples(self, mini_index):
        """SJ map keys should be (ref, start, end, strand) tuples."""
        for key in _sj_map(mini_index):
            assert isinstance(key, tuple)
            assert len(key) == 4
            ref, start, end, strand = key
            assert isinstance(ref, str)
            assert isinstance(start, (int, np.integer))
            assert isinstance(end, (int, np.integer))
            assert isinstance(strand, (int, np.integer))

    def test_sj_map_values_are_frozensets(self, mini_index):
        """SJ map values should be frozensets of ints."""
        for val in _sj_map(mini_index).values():
            assert isinstance(val, frozenset)
            for t_idx in val:
                assert isinstance(t_idx, (int, np.integer))
                assert t_idx >= 0

    def test_t0_introns_in_sj_map(self, mini_index):
        """t0's introns should be in the SJ map."""
        # t0 intron 1: 200-299 (between exons 99-200 and 299-400)
        # t0 intron 2: 400-499 (between exons 299-400 and 499-600)
        found_t0 = set()
        for (ref, start, end, strand), tset in _sj_map(mini_index).items():
            if ref == "chr1" and 0 in tset:
                found_t0.add((start, end))
        assert len(found_t0) >= 2, f"Expected at least 2 SJs for t0, got {found_t0}"

    def test_t2_intron_in_sj_map(self, mini_index):
        """t2's intron (1100-1199) should be in the SJ map."""
        found = False
        for (ref, start, end, strand), tset in _sj_map(mini_index).items():
            if ref == "chr1" and 2 in tset and start == 1100 and end == 1199:
                found = True
                assert strand == 2  # neg strand
                break
        assert found, "t2 intron 1100-1199 not found in the junction map"

    def test_the_resolver_and_calibration_hold_the_same_junctions(self, mini_index):
        assert _resolver_junctions(mini_index) == _calibration_junctions(mini_index)


def _sj_map(index) -> dict[tuple, frozenset[int]]:
    """``(ref, start, end, strand) → frozenset(t_index)``: the junctions :meth:`TranscriptIndex.load` installs in the
    resolver (`_sj_csr` of the index's junction table), each with the transcripts holding it."""
    refs, starts, ends, strands, t_flat, t_offsets = _sj_csr(
        pd.read_feather(Path(index.index_dir) / SJ_FEATHER)
    )
    return {
        (ref, start, end, strand): frozenset(t_flat[lo:hi])
        for ref, start, end, strand, lo, hi in zip(
            refs, starts, ends, strands, t_offsets[:-1], t_offsets[1:]
        )
    }


def _resolver_junctions(index) -> set[tuple[str, int, int, int]]:
    """The junctions the resolver matches a read's against (`_sj_map`), as (reference, start, end, strand)."""
    return {
        (str(ref), int(start), int(end), int(strand)) for ref, start, end, strand in _sj_map(index)
    }


def _calibration_junctions(index) -> set[tuple[str, int, int, int]]:
    """The junctions calibration's deposit matches a read's against: the splice graph's sj edges on the flat
    region-bound axis (`build_sj_arrays`), turned back into coordinates. Built along another path than `_sj_map`; a
    junction one holds and the other lacks would be spliced RNA to the EM and unspliced to calibration."""
    from rigel.calibration.splice_graph import build_region_partition_arrays, build_sj_arrays

    bounds, ref_offsets, _ = build_region_partition_arrays(index)
    sj = build_sj_arrays(index)
    donor = np.repeat(np.arange(len(sj.offsets) - 1), np.diff(sj.offsets))
    ref = np.searchsorted(ref_offsets, donor, side="right") - 1
    return {
        (index.ref_names[r], int(bounds[d]), int(bounds[a]), int(s))
        for r, d, a, s in zip(ref, donor, sj.boundary_right, sj.strand)
    }


# ═════════════════════════════════════════════════════════════════════
# C++ FragmentResolver
# ═════════════════════════════════════════════════════════════════════


class TestFragmentResolver:
    """Verify the native FragmentResolver is wired up."""

    def test_resolver_exists(self, mini_index):
        assert hasattr(mini_index, "resolver")
        assert mini_index.resolver is not None

    def test_the_builders_name_only_the_references_the_index_installed(self):
        """The reference ids are the index's (`set_ref_names`), and a BAM reference outside them maps to -1, off the
        index. A builder that minted an id for a reference it named would give that reference a place no partition
        holds, so an interval or a junction on one is an index defect, and refused."""
        from rigel.native import FragmentResolver

        resolver = FragmentResolver()
        resolver.set_ref_names(["chr1"])
        with pytest.raises(ValueError, match="chrDecoy"):
            resolver.build_overlap_index(
                ["chr1", "chrDecoy"], [0, 0], [10, 10], [0, 0], [0, 0], [0, 1, 2]
            )
        with pytest.raises(ValueError, match="chrDecoy"):
            resolver.build_sj_map(["chrDecoy"], [10], [20], [1], [0], [0, 1])
        resolver.build_sj_map(["chr1"], [10], [20], [1], [0], [0, 1])


# ═════════════════════════════════════════════════════════════════════
# Round-trip: build → load consistency
# ═════════════════════════════════════════════════════════════════════


class TestBuildLoadRoundTrip:
    """Build an index, load it, verify all structures agree."""

    def test_two_loads_produce_same_data(self, tmp_path):
        """Loading twice from the same dir should yield identical structures."""
        fasta_path, gtf_path = _mini_sources(tmp_path)
        idx_dir = tmp_path / "idx"
        TranscriptIndex.build(fasta_path, gtf_path, idx_dir, write_tsv=False)

        idx1 = TranscriptIndex.load(idx_dir)
        idx2 = TranscriptIndex.load(idx_dir)

        # Transcript tables
        assert idx1.t_df.equals(idx2.t_df)
        assert idx1.g_df.equals(idx2.g_df)
        np.testing.assert_array_equal(idx1.t_to_g_arr, idx2.t_to_g_arr)
        np.testing.assert_array_equal(idx1.t_to_strand_arr, idx2.t_to_strand_arr)

        # Exon CSR
        for a, b in zip(idx1.build_exon_csr(), idx2.build_exon_csr(), strict=True):
            np.testing.assert_array_equal(a, b)


# ═════════════════════════════════════════════════════════════════════
# Stress test: larger synthetic index
# ═════════════════════════════════════════════════════════════════════


class TestLargerIndex:
    """Build a larger index with many overlapping transcripts and verify."""

    COMPLEX_GTF = ""

    @pytest.fixture(scope="class")
    def complex_index(self, tmp_path_factory):
        """Build an index with 10 genes, ~30 transcripts, overlapping exons."""
        from _index_builder import build_test_index

        # Build a more complex GTF with multiple genes and overlapping exons
        lines = []
        # Gene A (+ strand): 3 transcripts with varying exon structures
        for exon in [(100, 300), (500, 700), (900, 1100)]:
            lines.append(
                f"chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t+\t.\t"
                f'gene_id "gA"; transcript_id "tA1"; gene_name "GeneA"; '
                f'gene_type "protein_coding"; tag "basic";'
            )
        for exon in [(100, 300), (900, 1100)]:
            lines.append(
                f"chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t+\t.\t"
                f'gene_id "gA"; transcript_id "tA2"; gene_name "GeneA"; '
                f'gene_type "protein_coding";'
            )
        for exon in [(100, 300), (500, 700), (900, 1100), (1300, 1500)]:
            lines.append(
                f"chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t+\t.\t"
                f'gene_id "gA"; transcript_id "tA3"; gene_name "GeneA"; '
                f'gene_type "protein_coding";'
            )

        # Gene B (- strand): 2 transcripts
        for exon in [(2000, 2200), (2500, 2700)]:
            lines.append(
                f"chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t-\t.\t"
                f'gene_id "gB"; transcript_id "tB1"; gene_name "GeneB"; '
                f'gene_type "lncRNA";'
            )
        for exon in [(2000, 2200), (2500, 2700), (3000, 3200)]:
            lines.append(
                f"chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t-\t.\t"
                f'gene_id "gB"; transcript_id "tB2"; gene_name "GeneB"; '
                f'gene_type "lncRNA";'
            )

        # Gene C (+ strand): overlapping with Gene A on same strand
        for exon in [(800, 1100), (1300, 1600)]:
            lines.append(
                f"chr1\ttest\texon\t{exon[0]}\t{exon[1]}\t.\t+\t.\t"
                f'gene_id "gC"; transcript_id "tC1"; gene_name "GeneC"; '
                f'gene_type "protein_coding";'
            )

        # Gene D (+ strand): single-exon gene (unspliced)
        lines.append(
            "chr1\ttest\texon\t4000\t5000\t.\t+\t.\t"
            'gene_id "gD"; transcript_id "tD1"; gene_name "GeneD"; '
            'gene_type "protein_coding"; tag "basic";'
        )

        gtf_text = "\n".join(lines) + "\n"
        return build_test_index(tmp_path_factory, gtf_text, genome_size=6000, name="complex_idx")

    def test_transcript_count(self, complex_index):
        # 7 annotated + 5 synthetic spans (gA:2, gB:2, gC:1)
        assert complex_index.num_transcripts == 12

    def test_gene_count(self, complex_index):
        # 4 annotated (gA, gB, gC, gD) + 5 synthetic gene-neutral rows
        assert complex_index.num_annotated_genes == 4
        assert complex_index.num_genes == 9

    def test_all_transcripts_have_exon_intervals(self, complex_index):
        for t_idx in range(complex_index.num_transcripts):
            assert len(complex_index.get_exon_intervals(t_idx)) > 0

    def test_exon_intervals_shape_per_transcript(self, complex_index):
        t_df = complex_index.t_df
        expected_exon_counts = {
            "tA1": 3,
            "tA2": 2,
            "tA3": 4,
            "tB1": 2,
            "tB2": 3,
            "tC1": 2,
            "tD1": 1,
        }
        for _, row in t_df.iterrows():
            t_idx = row["t_index"]
            t_id = row["t_id"]
            exons = complex_index.get_exon_intervals(t_idx)
            if t_id.startswith("RIGEL_NASCENT_"):
                # synthetic spans are single-exon
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
        for (ref, start, end, strand), tset in _sj_map(complex_index).items():
            assert td1_idx not in tset, (
                f"Single-exon transcript tD1 (t_idx={td1_idx}) should not "
                f"appear in SJ map, but found in ({ref},{start},{end},{strand})"
            )

    def test_the_resolver_and_calibration_hold_the_same_junctions(self, complex_index):
        assert _resolver_junctions(complex_index) == _calibration_junctions(complex_index)

    def test_exon_interval_values_are_int32(self, complex_index):
        for t_idx in range(complex_index.num_transcripts):
            exons = complex_index.get_exon_intervals(t_idx)
            assert exons.dtype == np.int32, f"t_idx={t_idx}: expected int32, got {exons.dtype}"


# ═════════════════════════════════════════════════════════════════════
# Edge cases
# ═════════════════════════════════════════════════════════════════════


class TestEdgeCases:
    """Edge cases: empty SJs, single-exon transcripts."""

    @pytest.fixture
    def single_exon_index(self, tmp_path_factory):
        """Index with only single-exon transcripts (no splice junctions)."""
        from _index_builder import build_test_index

        gtf = (
            "chr1\ttest\texon\t100\t500\t.\t+\t.\t"
            'gene_id "g1"; transcript_id "t1"; gene_name "G1"; '
            'gene_type "protein_coding"; tag "basic";\n'
            "chr1\ttest\texon\t1000\t1500\t.\t-\t.\t"
            'gene_id "g2"; transcript_id "t2"; gene_name "G2"; '
            'gene_type "lncRNA";\n'
        )
        return build_test_index(tmp_path_factory, gtf, genome_size=2000, name="single_exon_idx")

    def test_single_exon_loads_ok(self, single_exon_index):
        assert single_exon_index.num_transcripts == 2

    def test_single_exon_empty_sj_map(self, single_exon_index):
        assert len(_sj_map(single_exon_index)) == 0

    def test_single_exon_exon_intervals(self, single_exon_index):
        for t_idx in range(2):
            exons = single_exon_index.get_exon_intervals(t_idx)
            assert exons.shape == (1, 2)

    def test_single_exon_resolver_exists(self, single_exon_index):
        assert single_exon_index.resolver is not None


# ── What ``TranscriptIndex`` refuses: a reordered table, a missing file, a bad setting, an empty GTF, a half rebuild ─


def test_load_refuses_a_reordered_transcript_table(tmp_path: Path):
    """Every per-transcript array is keyed by row position, so reordered rows would mis-map them."""
    fasta, gtf = _mini_sources(tmp_path)
    idx_dir = tmp_path / "idx"
    TranscriptIndex.build(fasta, gtf, idx_dir, write_tsv=False)
    t_df = pd.read_feather(idx_dir / TRANSCRIPTS_FEATHER)
    t_df.iloc[::-1].reset_index(drop=True).to_feather(idx_dir / TRANSCRIPTS_FEATHER)

    with pytest.raises(IndexLoadError, match="row index does not match 't_index'"):
        TranscriptIndex.load(idx_dir)


def test_load_refuses_a_gene_index_with_a_gap(tmp_path: Path):
    """The gene table is keyed by row position too, so a ``g_index`` that skips a value would mis-map it."""
    fasta, gtf = _mini_sources(tmp_path)
    idx_dir = tmp_path / "idx"
    TranscriptIndex.build(fasta, gtf, idx_dir, write_tsv=False)
    t_df = pd.read_feather(idx_dir / TRANSCRIPTS_FEATHER)
    t_df["g_index"] += 1
    t_df.to_feather(idx_dir / TRANSCRIPTS_FEATHER)

    with pytest.raises(IndexLoadError, match="does not match 'g_index'"):
        TranscriptIndex.load(idx_dir)


@pytest.mark.parametrize(
    "required",
    [
        REF_LENGTHS_FEATHER,
        TRANSCRIPTS_FEATHER,
        INTERVALS_FEATHER,
        REGIONS_FEATHER,
        EDGES_FEATHER,
        SJ_FEATHER,
    ],
)
def test_load_refuses_an_index_missing_a_file_it_reads_by_name(tmp_path: Path, required):
    """`rigel quant` reports a load refusal as an error, so a missing file is refused by name rather than
    escaping as pyarrow's ``FileNotFoundError`` or a bare ``RuntimeError``, a traceback with no remedy."""
    fasta, gtf = _mini_sources(tmp_path)
    idx_dir = tmp_path / "idx"
    TranscriptIndex.build(fasta, gtf, idx_dir, write_tsv=False)
    (idx_dir / required).unlink()

    with pytest.raises(IndexLoadError, match=rf"(?s)is missing {re.escape(required)}\b.*Rebuild"):
        TranscriptIndex.load(idx_dir)


def test_load_does_not_need_the_alias_table(tmp_path: Path):
    """``transcript_aliases.feather`` is written for the user and never read, so an index without it loads."""
    fasta, gtf = _mini_sources(tmp_path)
    idx_dir = tmp_path / "idx"
    TranscriptIndex.build(fasta, gtf, idx_dir, write_tsv=False)
    (idx_dir / TRANSCRIPT_ALIASES_FEATHER).unlink()

    assert TranscriptIndex.load(idx_dir).num_transcripts == 5


@pytest.mark.parametrize(
    "setting, refusal",
    [
        ({"end_tolerance": -5}, "end_tolerance must be >= 0; got -5"),
        ({"feather_compression": "gzip"}, "Unknown Feather compression: 'gzip'"),
        ({"gtf_parse_mode": "lenient"}, "Unknown GTF parse mode: 'lenient'"),
    ],
    ids=["negative-tolerance", "compression", "parse-mode"],
)
def test_build_refuses_a_setting_its_config_refuses_before_writing(
    tmp_path: Path, setting, refusal
):
    """`build` validates its settings as one `IndexConfig`, for every caller: a negative tolerance would cluster as 0
    does while the manifest records it as given."""
    fasta, gtf = _mini_sources(tmp_path)
    idx_dir = tmp_path / "idx"

    with pytest.raises(ValueError, match=re.escape(refusal)):
        TranscriptIndex.build(fasta, gtf, idx_dir, write_tsv=False, **setting)
    assert not idx_dir.exists()


def test_build_refuses_a_gtf_without_transcripts_before_writing(tmp_path: Path):
    fasta, gtf = _mini_sources(tmp_path)
    gtf.write_text("")
    idx_dir = tmp_path / "idx"

    with pytest.raises(ValueError, match="no transcripts"):
        TranscriptIndex.build(fasta, gtf, idx_dir, write_tsv=False)
    assert list(idx_dir.iterdir()) == []


def test_a_reference_longer_than_an_int32_coordinate_is_refused_by_name(tmp_path):
    """A coordinate is an int32 from the index to the scan, as it is in BAM, so a reference longer than 2^31 - 1 bp
    cannot be addressed. pysam reports a length as a C int, clipped at that largest value, so the build reads the
    FASTA's own index: the reference is refused by name at its true length, never recorded 2,147,483,647 bp long."""
    too_long = int(np.iinfo(np.int32).max) + 1
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n>chrHuge\nACGT\n")
    # the index of a FASTA whose chrHuge is one base too long; only the index is read for a length
    (tmp_path / "genome.fa.fai").write_text(
        f"chr1\t1000\t6\t1000\t1001\nchrHuge\t{too_long}\t1016\t4\t5\n"
    )
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text('chr1\ttest\texon\t100\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
    idx_dir = tmp_path / "idx"
    with pytest.raises(ValueError, match=f"chrHuge \\({too_long:,} bp\\)"):
        TranscriptIndex.build(fasta, gtf, idx_dir, write_tsv=False)
    assert list(idx_dir.iterdir()) == []


def test_an_empty_gtf_rebuild_leaves_the_existing_index_loadable(tmp_path: Path):
    """The refusal comes before any old file is deleted, the manifest among them, so a rebuild refused
    for an empty GTF leaves the index it would have replaced untouched and loadable."""
    fasta, gtf = _mini_sources(tmp_path)
    idx_dir = tmp_path / "idx"
    TranscriptIndex.build(fasta, gtf, idx_dir, write_tsv=False)
    files_before = {p.name: p.read_bytes() for p in idx_dir.iterdir()}
    n_transcripts = TranscriptIndex.load(idx_dir).num_transcripts

    empty = tmp_path / "empty.gtf"
    empty.write_text("")
    with pytest.raises(ValueError, match="no transcripts"):
        TranscriptIndex.build(fasta, empty, idx_dir, write_tsv=False)

    assert {p.name: p.read_bytes() for p in idx_dir.iterdir()} == files_before
    assert TranscriptIndex.load(idx_dir).num_transcripts == n_transcripts


def test_an_interrupted_rebuild_leaves_no_manifest(tmp_path: Path, monkeypatch):
    """The manifest is written last, so a rebuild that stops part-way must not leave the old one
    vouching for a mix of old and new files. The failure is injected at the first write, so an old
    manifest removed only after an index file is written fails here too."""
    fasta, gtf = _mini_sources(tmp_path)
    idx_dir = tmp_path / "idx"
    TranscriptIndex.build(fasta, gtf, idx_dir, write_tsv=False)

    def interrupted(*args, **kwargs):
        raise RuntimeError("interrupted")

    monkeypatch.setattr(pd.DataFrame, "to_feather", interrupted)
    with pytest.raises(RuntimeError, match="interrupted"):
        TranscriptIndex.build(fasta, gtf, idx_dir, write_tsv=False)

    assert (idx_dir / TRANSCRIPTS_FEATHER).exists()
    with pytest.raises(RuntimeError, match="no manifest.json"):
        TranscriptIndex.load(idx_dir)


# ── The duplicate-transcript guard in ``read_transcripts``, and its collapse flag ────────────────


def _write_gtf(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n")


def _gtf_line(
    feature: str, start: int, end: int, strand: str, gid: str, tid: str, ref: str = "chr1"
) -> str:
    return (
        f"{ref}\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t"
        f'gene_id "{gid}"; transcript_id "{tid}";'
    )


def test_duplicate_exon_structure_raises(tmp_path: Path):
    """Two transcripts with identical exon coordinates must raise ValueError."""
    gtf = tmp_path / "dup.gtf"
    _write_gtf(
        gtf,
        [
            _gtf_line("transcript", 100, 500, "+", "G1", "T1"),
            _gtf_line("exon", 100, 200, "+", "G1", "T1"),
            _gtf_line("exon", 400, 500, "+", "G1", "T1"),
            # Identical exon set, different transcript_id
            _gtf_line("transcript", 100, 500, "+", "G1", "T2"),
            _gtf_line("exon", 100, 200, "+", "G1", "T2"),
            _gtf_line("exon", 400, 500, "+", "G1", "T2"),
        ],
    )
    with pytest.raises(ValueError, match=r"identical exon coordinates"):
        read_transcripts(gtf)


def test_same_intron_chain_different_utrs_allowed(tmp_path: Path):
    """Transcripts that share an intron chain but differ in UTR boundaries
    are biologically realistic and must not trigger the guard."""
    gtf = tmp_path / "utrs.gtf"
    _write_gtf(
        gtf,
        [
            _gtf_line("transcript", 100, 500, "+", "G1", "T1"),
            _gtf_line("exon", 100, 200, "+", "G1", "T1"),
            _gtf_line("exon", 400, 500, "+", "G1", "T1"),
            # Same intron (200,400) but different 5' and 3' UTR
            _gtf_line("transcript", 120, 480, "+", "G1", "T2"),
            _gtf_line("exon", 120, 200, "+", "G1", "T2"),
            _gtf_line("exon", 400, 480, "+", "G1", "T2"),
        ],
    )
    txs = read_transcripts(gtf)
    assert len(txs) == 2
    assert {t.t_id for t in txs} == {"T1", "T2"}


def test_different_strand_not_a_duplicate(tmp_path: Path):
    """Identical coords on opposite strands are distinct transcripts."""
    gtf = tmp_path / "strand.gtf"
    _write_gtf(
        gtf,
        [
            _gtf_line("transcript", 100, 500, "+", "G1", "T1"),
            _gtf_line("exon", 100, 200, "+", "G1", "T1"),
            _gtf_line("exon", 400, 500, "+", "G1", "T1"),
            _gtf_line("transcript", 100, 500, "-", "G2", "T2"),
            _gtf_line("exon", 100, 200, "-", "G2", "T2"),
            _gtf_line("exon", 400, 500, "-", "G2", "T2"),
        ],
    )
    txs = read_transcripts(gtf)
    assert len(txs) == 2


def test_duplicate_message_lists_offending_ids(tmp_path: Path):
    gtf = tmp_path / "dup.gtf"
    _write_gtf(
        gtf,
        [
            _gtf_line("transcript", 100, 500, "+", "G1", "ALPHA"),
            _gtf_line("exon", 100, 500, "+", "G1", "ALPHA"),
            _gtf_line("transcript", 100, 500, "+", "G1", "BETA"),
            _gtf_line("exon", 100, 500, "+", "G1", "BETA"),
        ],
    )
    with pytest.raises(ValueError) as excinfo:
        read_transcripts(gtf)
    msg = str(excinfo.value)
    assert "ALPHA" in msg
    assert "BETA" in msg


def test_duplicate_error_mentions_collapse_flag(tmp_path: Path):
    """The error should point users at the --collapse-duplicate-transcripts escape hatch."""
    gtf = tmp_path / "dup.gtf"
    _write_gtf(
        gtf,
        [
            _gtf_line("transcript", 100, 500, "+", "G1", "T1"),
            _gtf_line("exon", 100, 500, "+", "G1", "T1"),
            _gtf_line("transcript", 100, 500, "+", "G1", "T2"),
            _gtf_line("exon", 100, 500, "+", "G1", "T2"),
        ],
    )
    with pytest.raises(ValueError, match=r"--collapse-duplicate-transcripts"):
        read_transcripts(gtf)


def test_collapse_keeps_lexicographically_smallest_id(tmp_path: Path):
    """With collapse enabled, a duplicate group keeps only the lexicographically
    smallest transcript ID; distinct transcripts are untouched and t_index is
    reassigned contiguously."""
    gtf = tmp_path / "collapse.gtf"
    _write_gtf(
        gtf,
        [
            # duplicate pair (identical exons), IDs deliberately out of lexical order
            _gtf_line("transcript", 100, 500, "+", "G1", "ENST00000002"),
            _gtf_line("exon", 100, 200, "+", "G1", "ENST00000002"),
            _gtf_line("exon", 400, 500, "+", "G1", "ENST00000002"),
            _gtf_line("transcript", 100, 500, "+", "G1", "ENST00000001"),
            _gtf_line("exon", 100, 200, "+", "G1", "ENST00000001"),
            _gtf_line("exon", 400, 500, "+", "G1", "ENST00000001"),
            # a distinct transcript (different exons) — must survive
            _gtf_line("transcript", 1000, 1500, "+", "G2", "ENST00000009"),
            _gtf_line("exon", 1000, 1500, "+", "G2", "ENST00000009"),
        ],
    )
    txs = read_transcripts(gtf, collapse_duplicate_transcripts=True)
    assert {t.t_id for t in txs} == {"ENST00000001", "ENST00000009"}
    # dropped transcript is gone; kept indices are contiguous 0..N-1
    assert sorted(t.t_index for t in txs) == [0, 1]


def test_collapse_three_way_group_keeps_one(tmp_path: Path):
    """A 3-way duplicate group collapses to the single lexicographically smallest ID."""
    gtf = tmp_path / "collapse3.gtf"
    lines: list[str] = []
    for tid in ("ENST0000C", "ENST0000A", "ENST0000B"):
        lines += [
            _gtf_line("transcript", 100, 500, "+", "G1", tid),
            _gtf_line("exon", 100, 200, "+", "G1", tid),
            _gtf_line("exon", 400, 500, "+", "G1", tid),
        ]
    _write_gtf(gtf, lines)
    txs = read_transcripts(gtf, collapse_duplicate_transcripts=True)
    assert {t.t_id for t in txs} == {"ENST0000A"}


def test_collapsed_duplicates_are_written_as_an_alias_map(tmp_path: Path):
    """A user whose annotation names a collapsed transcript finds the one that stands for it; with nothing
    collapsed the table is written empty."""
    from rigel.index import TRANSCRIPT_ALIASES_FEATHER

    fasta, _ = _mini_sources(tmp_path)
    gtf = tmp_path / "collapse.gtf"
    _write_gtf(
        gtf,
        [
            _gtf_line("transcript", 100, 500, "+", "G1", "ENST00000002"),
            _gtf_line("exon", 100, 200, "+", "G1", "ENST00000002"),
            _gtf_line("exon", 400, 500, "+", "G1", "ENST00000002"),
            _gtf_line("transcript", 100, 500, "+", "G1", "ENST00000001"),
            _gtf_line("exon", 100, 200, "+", "G1", "ENST00000001"),
            _gtf_line("exon", 400, 500, "+", "G1", "ENST00000001"),
            _gtf_line("transcript", 1000, 1500, "+", "G2", "ENST00000009"),
            _gtf_line("exon", 1000, 1500, "+", "G2", "ENST00000009"),
        ],
    )
    TranscriptIndex.build(fasta, gtf, tmp_path / "idx", collapse_duplicate_transcripts=True)
    aliases = pd.read_feather(tmp_path / "idx" / TRANSCRIPT_ALIASES_FEATHER)
    assert aliases.to_dict("records") == [
        {"dropped_t_id": "ENST00000002", "kept_t_id": "ENST00000001"}
    ]

    (tmp_path / "plain").mkdir()
    _, plain_gtf = _mini_sources(tmp_path / "plain")
    TranscriptIndex.build(fasta, plain_gtf, tmp_path / "plain_idx")
    plain = pd.read_feather(tmp_path / "plain_idx" / TRANSCRIPT_ALIASES_FEATHER)
    assert list(plain.columns) == ["dropped_t_id", "kept_t_id"] and plain.empty


def test_the_index_writes_no_column_nothing_reads(tmp_path: Path):
    """The simulator's abundances and the exon count in ``transcripts.feather``, and the constant interval type
    in ``sj.feather``, are not written: nothing reads them."""
    from rigel.index import SJ_FEATHER

    fasta, gtf = _mini_sources(tmp_path)
    TranscriptIndex.build(fasta, gtf, tmp_path / "idx")
    t_columns = set(pd.read_feather(tmp_path / "idx" / TRANSCRIPTS_FEATHER).columns)
    assert not t_columns & {"abundance", "nrna_abundance", "n_exons"}
    assert "interval_type" not in pd.read_feather(tmp_path / "idx" / SJ_FEATHER).columns
    assert TranscriptIndex.load(tmp_path / "idx").num_transcripts > 0


def test_a_single_exon_transcript_covering_two_spans_counts_the_contributors_of_both():
    """The nascent entity a single-exon transcript stands for holds every multi-exon transcript of every span it
    covers: one span of one transcript and one of two, three in all. An overwrite would keep only the last span's
    two."""
    from rigel.index import create_nascent_spans

    def tx(t_index, t_id, exons):
        return Transcript(
            ref="chr1",
            strand=Strand.POS,
            exons=[Interval(s, e) for s, e in exons],
            t_id=t_id,
            g_id="G1",
            t_index=t_index,
            g_index=0,
        )

    covering = tx(0, "SE", [(100, 2000)])
    transcripts = [
        covering,
        tx(1, "T1", [(200, 300), (400, 500)]),
        tx(2, "T2", [(1000, 1100), (1200, 1300)]),
        tx(3, "T3", [(1000, 1100), (1250, 1300)]),
    ]
    synthetics, _t_to_span_key, carriers = create_nascent_spans(transcripts, tolerance=20)
    assert synthetics == [] and len(carriers) == 2
    assert all(entity is covering for entity in carriers.values())
    assert covering.span_n_contributors == 3


def _span_entities(transcripts, tolerance=20):
    """The nascent entities' extents: the synthetics' and the carrying single-exon transcripts', and the carriers."""
    from rigel.index import create_nascent_spans

    synthetics, _t_to_span_key, carriers = create_nascent_spans(transcripts, tolerance=tolerance)
    extents = {(s.start, s.end) for s in synthetics} | {(t.start, t.end) for t in carriers.values()}
    return sorted(extents), carriers


def _multi_exon(t_index, start, end):
    """A two-exon transcript over [start, end) on chr1, +."""
    return Transcript(
        ref="chr1",
        strand=Strand.POS,
        exons=[Interval(start, start + 50), Interval(end - 50, end)],
        t_id=f"T{t_index}",
        g_id="G1",
        t_index=t_index,
        g_index=0,
    )


def _within(a, b, tolerance=20):
    return abs(a[0] - b[0]) <= tolerance and abs(a[1] - b[1]) <= tolerance


def test_no_two_spans_lie_within_the_tolerance_at_both_ends():
    """Ends 1,000, 1,020 and 1,021: a cluster anchored at its first coordinate holds 1,000 and 1,020 and starts
    another at 1,021, leaving two spans 1 bp apart. Within the tolerance at both ends is one entity."""
    spans, _ = _span_entities([_multi_exon(i, 0, end) for i, end in enumerate((1000, 1020, 1021))])
    assert spans == [(0, 1021)]


def test_a_chain_within_the_tolerance_is_one_span():
    """Ends 1,000, 1,015 and 1,030, each within 20 bp of the next: one entity over their envelope."""
    spans, _ = _span_entities([_multi_exon(i, 0, end) for i, end in enumerate((1000, 1015, 1030))])
    assert spans == [(0, 1030)]


def test_envelopes_within_the_tolerance_are_one_span():
    """A (100, 1,000), X (110, 1,015) and B (115, 1,030) chain into one entity over (100, 1,030); D (85, 1,045) is
    near none of them but within 20 bp of their envelope at both ends, so it joins: no two spans left near."""
    points = [(100, 1000), (110, 1015), (115, 1030), (85, 1045)]
    spans, _ = _span_entities([_multi_exon(i, s, e) for i, (s, e) in enumerate(points)])
    assert spans == [(85, 1045)]


def test_a_span_keeps_its_start_when_its_neighbour_ends_far_away():
    """Starts 0 and 15, ends 1,000 and 5,000: two entities, each its own extent. Clustering starts and ends
    separately would give the second the first's start."""
    spans, _ = _span_entities([_multi_exon(0, 0, 1000), _multi_exon(1, 15, 5000)])
    assert spans == [(0, 1000), (15, 5000)]


def test_no_span_is_made_equal_to_a_single_exon_transcript():
    """A single-exon transcript is its own nascent RNA: one 2 bp and 7 bp inside a span's ends (TMEM121B's case)
    is that span's entity, and no synthetic span is made beside it."""
    single = Transcript(
        ref="chr1",
        strand=Strand.POS,
        exons=[Interval(2, 993)],
        t_id="SE",
        g_id="G1",
        t_index=0,
        g_index=0,
    )
    spans, covered = _span_entities([single, _multi_exon(1, 0, 1000)])
    assert spans == [(2, 993)] and list(covered.values()) == [single]


def test_a_single_exon_transcript_near_a_chains_member_but_not_its_span_stands_for_nothing():
    """TDRG1's case: contributors chain from start 100 to start 160, so the span starts at 100; a single-exon
    transcript within 20 bp of the last contributor but 60 bp short of the span is not the span's equal, so the span is
    a synthetic of its own (it does not contain it either)."""
    single = Transcript(
        ref="chr1",
        strand=Strand.POS,
        exons=[Interval(165, 1000)],
        t_id="SE",
        g_id="G1",
        t_index=0,
        g_index=0,
    )
    chain = [_multi_exon(i + 1, s, 1000) for i, s in enumerate((100, 115, 130, 145, 160))]
    spans, covered = _span_entities([single, *chain])
    assert covered == {} and spans == [(100, 1000)]


def test_the_tolerance_is_the_setting():
    """Ends 1,000 and 1,030 are one entity at a tolerance of 30 and two at 20."""
    transcripts = [_multi_exon(0, 0, 1000), _multi_exon(1, 0, 1030)]
    assert _span_entities(transcripts, tolerance=30)[0] == [(0, 1030)]
    assert _span_entities(transcripts, tolerance=20)[0] == [(0, 1000), (0, 1030)]


def test_collapse_is_noop_without_duplicates(tmp_path: Path):
    """collapse=True must not drop transcripts that are merely similar (shared intron,
    different UTRs)."""
    gtf = tmp_path / "nodup.gtf"
    _write_gtf(
        gtf,
        [
            _gtf_line("transcript", 100, 500, "+", "G1", "T1"),
            _gtf_line("exon", 100, 200, "+", "G1", "T1"),
            _gtf_line("exon", 400, 500, "+", "G1", "T1"),
            _gtf_line("transcript", 120, 480, "+", "G1", "T2"),
            _gtf_line("exon", 120, 200, "+", "G1", "T2"),
            _gtf_line("exon", 400, 480, "+", "G1", "T2"),
        ],
    )
    txs = read_transcripts(gtf, collapse_duplicate_transcripts=True)
    assert {t.t_id for t in txs} == {"T1", "T2"}


# ── ``rigel.index.iter_reference_layout`` — the tiling both feathers are derived from ───────────
#
# The layout iterator is the single source of truth for ``intervals.feather`` and
# ``regions.feather`` alike, so what it yields must be a perfect tiling of every reference:
# alternating INTERGENIC / GENIC spans, no overlaps, no gaps, synthetic spans excluded.


def _mk_tx(
    t_idx: int, ref: str, strand: Strand, exons: list[tuple[int, int]], is_synthetic: bool = False
) -> Transcript:
    """Construct a minimal Transcript for layout tests."""
    return Transcript(
        ref=ref,
        strand=strand,
        exons=[Interval(s, e) for s, e in exons],
        t_id=f"t{t_idx}",
        g_id=f"g{t_idx}",
        t_index=t_idx,
        g_index=t_idx,
        is_synthetic=is_synthetic,
    )


def _assert_tiles(spans, ref_length: int):
    """Assert that ``spans`` exactly tile [0, ref_length) with alternating types."""
    assert spans, "expected at least one span"
    cursor = 0
    last_type = None
    for s in spans:
        assert s.start == cursor, f"gap or overlap before {s} (cursor={cursor})"
        assert s.end > s.start
        cur_type = type(s).__name__
        assert cur_type != last_type, f"two adjacent {cur_type} spans"
        last_type = cur_type
        cursor = s.end
    assert cursor == ref_length, f"layout ends at {cursor}, expected {ref_length}"


def test_empty_reference_yields_single_intergenic():
    spans = list(iter_reference_layout(1000, []))
    assert spans == [IntergenicSpan(0, 1000)]


def test_zero_length_reference_yields_nothing():
    assert list(iter_reference_layout(0, [])) == []


def test_single_transcript_in_middle_three_spans():
    t = _mk_tx(0, "chr1", Strand.POS, [(100, 200), (300, 400)])
    spans = list(iter_reference_layout(1000, [t]))
    _assert_tiles(spans, 1000)
    assert len(spans) == 3
    assert spans[0] == IntergenicSpan(0, 100)
    assert isinstance(spans[1], GenicSpan)
    assert (spans[1].start, spans[1].end) == (100, 400)
    assert spans[1].transcripts == (t,)
    assert spans[2] == IntergenicSpan(400, 1000)


def test_transcript_at_left_boundary_no_left_intergenic():
    t = _mk_tx(0, "chr1", Strand.POS, [(0, 200)])
    spans = list(iter_reference_layout(500, [t]))
    _assert_tiles(spans, 500)
    assert len(spans) == 2
    assert isinstance(spans[0], GenicSpan)
    assert spans[1] == IntergenicSpan(200, 500)


def test_transcript_at_right_boundary_no_right_intergenic():
    t = _mk_tx(0, "chr1", Strand.POS, [(800, 1000)])
    spans = list(iter_reference_layout(1000, [t]))
    _assert_tiles(spans, 1000)
    assert len(spans) == 2
    assert spans[0] == IntergenicSpan(0, 800)
    assert isinstance(spans[1], GenicSpan)


def test_overlapping_transcripts_strand_agnostic_coalesce():
    """Two transcripts on opposite strands that overlap → one genic span."""
    t1 = _mk_tx(0, "chr1", Strand.POS, [(100, 300)])
    t2 = _mk_tx(1, "chr1", Strand.NEG, [(250, 500)])
    spans = list(iter_reference_layout(1000, sorted([t1, t2], key=lambda x: (x.start, x.end))))
    _assert_tiles(spans, 1000)
    genic = [s for s in spans if isinstance(s, GenicSpan)]
    assert len(genic) == 1
    assert (genic[0].start, genic[0].end) == (100, 500)
    assert {t.t_index for t in genic[0].transcripts} == {0, 1}


def test_disjoint_transcripts_yield_separate_genic_spans():
    t1 = _mk_tx(0, "chr1", Strand.POS, [(100, 200)])
    t2 = _mk_tx(1, "chr1", Strand.POS, [(500, 700)])
    spans = list(iter_reference_layout(1000, [t1, t2]))
    _assert_tiles(spans, 1000)
    assert len(spans) == 5  # int, gen, int, gen, int
    assert sum(isinstance(s, GenicSpan) for s in spans) == 2


def test_synthetic_excluded_from_layout():
    """Synthetic spans must not coalesce or extend genic spans."""
    real = _mk_tx(0, "chr1", Strand.POS, [(100, 200)])
    syn = _mk_tx(1, "chr1", Strand.POS, [(150, 800)], is_synthetic=True)
    spans = list(iter_reference_layout(1000, sorted([real, syn], key=lambda x: (x.start, x.end))))
    _assert_tiles(spans, 1000)
    # The synthetic [150,800) must not extend the real [100,200) genic span.
    genic = [s for s in spans if isinstance(s, GenicSpan)]
    assert len(genic) == 1
    assert (genic[0].start, genic[0].end) == (100, 200)
    assert genic[0].transcripts == (real,)


def test_touching_transcripts_coalesce():
    """t1.end == t2.start → coalesced into one genic span (touching boundary)."""
    t1 = _mk_tx(0, "chr1", Strand.POS, [(100, 200)])
    t2 = _mk_tx(1, "chr1", Strand.POS, [(200, 300)])
    spans = list(iter_reference_layout(1000, [t1, t2]))
    _assert_tiles(spans, 1000)
    genic = [s for s in spans if isinstance(s, GenicSpan)]
    assert len(genic) == 1
    assert (genic[0].start, genic[0].end) == (100, 300)


def test_the_resolver_inputs_are_the_per_run_specification():
    """The index load builds the resolver's inputs vectorised: the collapsed overlap index (one interval per distinct
    ``(ref, start, end, type)``, its distinct non-negative transcripts sorted), the junction map (each distinct
    ``(ref, start, end, strand)`` with its distinct transcripts sorted) and the exon CSR (each transcript's exons by
    start, a synthetic's own span replacing anything the table holds for it). Each must equal a run-by-run construction
    on a table with duplicate rows, the intergenic sentinel, several references and synthetic transcripts.

    Perturbation: keeping the sentinel, losing a run's sort, or letting a synthetic keep the table's exons fires it."""
    from rigel.index import _collapsed_intervals, _exon_csr, _sj_csr
    from rigel.types import IntervalType

    rng = np.random.default_rng(17)
    n, n_t = 4_000, 300
    iv = pd.DataFrame(
        {
            "ref": rng.choice(["chr10", "chr2", "chrX"], n),
            "start": rng.integers(0, 500, n),
            "t_index": rng.integers(-1, n_t, n),
            "interval_type": rng.choice([int(t) for t in IntervalType], n),
        }
    )
    iv["end"] = iv["start"] + rng.integers(1, 50, n)
    iv = pd.concat([iv, iv.sample(500, random_state=1)], ignore_index=True)

    def runs(df, keys):
        groups: dict[tuple, set] = {}
        for row in df.itertuples(index=False):
            groups.setdefault(tuple(getattr(row, k) for k in keys), set()).add(int(row.t_index))
        return sorted(groups.items())

    collapsed = _collapsed_intervals(iv)
    want = runs(iv, ("ref", "start", "end", "interval_type"))
    assert [tuple(x) for x in zip(*collapsed[:4])] == [k for k, _ in want]
    members = [sorted(t for t in ts if t >= 0) for _, ts in want]
    assert collapsed[5] == list(np.cumsum([0] + [len(m) for m in members]))
    assert collapsed[4] == [t for m in members for t in m]
    assert any(-1 in ts for _, ts in want) and any(len(m) > 1 for m in members)

    sj = iv.rename(columns={"interval_type": "strand"})
    sj["t_index"] = sj["t_index"].abs()
    csr = _sj_csr(sj)
    want = runs(sj, ("ref", "start", "end", "strand"))
    assert [tuple(x) for x in zip(*csr[:4])] == [k for k, _ in want]
    assert csr[4] == [t for _, ts in want for t in sorted(ts)]

    synthetic = np.zeros(n_t, dtype=bool)
    synthetic[rng.choice(n_t, 40, replace=False)] = True
    t_df = pd.DataFrame(
        {
            "is_synthetic": synthetic,
            "start": rng.integers(0, 500, n_t),
            "end": rng.integers(600, 900, n_t),
        }
    )
    offsets, starts, ends, before = _exon_csr(iv, t_df)
    exon = iv[(iv["interval_type"] == int(IntervalType.EXON)) & (iv["t_index"] >= 0)]
    for t in range(n_t):
        if synthetic[t]:
            expected = [(int(t_df["start"][t]), int(t_df["end"][t]))]
        else:
            rows = exon[exon["t_index"] == t].sort_values("start", kind="stable")
            expected = list(zip(rows["start"].tolist(), rows["end"].tolist()))
        lo, hi = int(offsets[t]), int(offsets[t + 1])
        assert list(zip(starts[lo:hi].tolist(), ends[lo:hi].tolist())) == expected, (
            f"transcript {t}"
        )
        lengths = [e - s for s, e in expected]
        assert before[lo:hi].tolist() == np.cumsum([0] + lengths)[:-1].tolist(), f"transcript {t}"
    assert synthetic[exon["t_index"]].any(), (
        "some synthetic must hold table exons, which its span replaces"
    )
