"""Regression tests for unambiguous intronic interval computation.

These tests guard against the nRNA phantom-count bug where ambiguous
intronic regions (intronic for transcript A but exonic for transcript B)
were incorrectly counted as evidence for nascent RNA, causing massive
false nRNA assignment.

Tests cover three levels:
1. Index build: UNAMBIG_INTRON intervals generated correctly
2. Resolution: unambig_intron_bp accumulated correctly from cgranges
3. End-to-end: zero-nRNA scenario produces zero phantom nRNA
"""

import pytest

from hulkrna.index import (
    _gen_cluster_unambig_intron_intervals,
    _merge_exon_intervals,
    _subtract_from_interval,
    _gen_genomic_intervals,
    build_genomic_intervals,
    TranscriptIndex,
)
from hulkrna.sim import Scenario, SimConfig, run_benchmark
from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
from hulkrna.pipeline import run_pipeline
from hulkrna.transcript import Transcript
from hulkrna.types import IntervalType, Strand, Interval


# =====================================================================
# Helper interval math
# =====================================================================


class TestMergeExonIntervals:
    def test_empty(self):
        assert _merge_exon_intervals([]) == []

    def test_single(self):
        assert _merge_exon_intervals([(10, 20)]) == [(10, 20)]

    def test_non_overlapping(self):
        assert _merge_exon_intervals([(10, 20), (30, 40)]) == [
            (10, 20), (30, 40),
        ]

    def test_overlapping(self):
        assert _merge_exon_intervals([(10, 30), (20, 40)]) == [(10, 40)]

    def test_adjacent(self):
        assert _merge_exon_intervals([(10, 20), (20, 30)]) == [(10, 30)]

    def test_multiple_merge(self):
        assert _merge_exon_intervals([(10, 25), (20, 35), (50, 60)]) == [
            (10, 35), (50, 60),
        ]

    def test_unsorted_input(self):
        assert _merge_exon_intervals([(30, 40), (10, 20)]) == [
            (10, 20), (30, 40),
        ]


class TestSubtractFromInterval:
    def test_no_subtraction(self):
        result = _subtract_from_interval(10, 50, [])
        assert result == [(10, 50)]

    def test_full_subtraction(self):
        result = _subtract_from_interval(10, 50, [(0, 100)])
        assert result == []

    def test_partial_subtraction_left(self):
        result = _subtract_from_interval(10, 50, [(0, 30)])
        assert result == [(30, 50)]

    def test_partial_subtraction_right(self):
        result = _subtract_from_interval(10, 50, [(30, 100)])
        assert result == [(10, 30)]

    def test_middle_subtraction(self):
        result = _subtract_from_interval(10, 50, [(20, 30)])
        assert result == [(10, 20), (30, 50)]

    def test_multiple_subtractions(self):
        result = _subtract_from_interval(10, 100, [(20, 30), (50, 60)])
        assert result == [(10, 20), (30, 50), (60, 100)]

    def test_empty_interval(self):
        result = _subtract_from_interval(10, 10, [(5, 15)])
        assert result == []


# =====================================================================
# Cluster-level UNAMBIG_INTRON generation
# =====================================================================


def _make_transcript(ref, strand, exons, t_index):
    """Create a Transcript with specified exons."""
    t = Transcript()
    t.ref = ref
    t.strand = strand
    t.exons = [Interval(s, e) for s, e in exons]
    t.t_index = t_index
    t.g_index = 0
    t.t_id = f"t{t_index}"
    t.g_id = "g0"
    return t


class TestGenClusterUnambigIntronIntervals:
    def test_single_multi_exon_no_overlap(self):
        """One 2-exon transcript alone → full intron is unambiguous."""
        t0 = _make_transcript("chr1", Strand.POS,
                              [(100, 200), (400, 500)], t_index=0)
        intervals = list(_gen_cluster_unambig_intron_intervals([t0]))
        assert len(intervals) == 1
        iv = intervals[0]
        assert iv.interval_type == IntervalType.UNAMBIG_INTRON
        assert iv.t_index == 0
        assert iv.start == 200
        assert iv.end == 400

    def test_single_exon_transcript_no_introns(self):
        """Single-exon transcript → no UNAMBIG_INTRON intervals."""
        t0 = _make_transcript("chr1", Strand.POS,
                              [(100, 500)], t_index=0)
        intervals = list(_gen_cluster_unambig_intron_intervals([t0]))
        assert len(intervals) == 0

    def test_two_transcripts_full_exon_cover(self):
        """t1's exon fully covers t0's intron → no unambig introns for t0."""
        t0 = _make_transcript("chr1", Strand.POS,
                              [(100, 200), (400, 500)], t_index=0)
        t1 = _make_transcript("chr1", Strand.POS,
                              [(100, 500)], t_index=1)
        intervals = list(_gen_cluster_unambig_intron_intervals([t0, t1]))
        # t0's intron 200-400 is fully covered by t1's exon 100-500
        # t1 is single-exon, no introns
        assert len(intervals) == 0

    def test_two_transcripts_partial_exon_cover(self):
        """t1's exon partially covers t0's intron → reduced unambig intron."""
        t0 = _make_transcript("chr1", Strand.POS,
                              [(100, 200), (400, 500)], t_index=0)
        t1 = _make_transcript("chr1", Strand.POS,
                              [(300, 400)], t_index=1)
        intervals = list(_gen_cluster_unambig_intron_intervals([t0, t1]))
        # t0's intron 200-400, t1's exon 300-400
        # → unambig: 200-300
        assert len(intervals) == 1
        assert intervals[0].start == 200
        assert intervals[0].end == 300
        assert intervals[0].t_index == 0

    def test_three_exon_transcript_two_introns(self):
        """3-exon transcript → two introns, each evaluated independently."""
        t0 = _make_transcript("chr1", Strand.POS,
                              [(100, 200), (400, 500), (700, 800)], t_index=0)
        intervals = list(_gen_cluster_unambig_intron_intervals([t0]))
        assert len(intervals) == 2
        starts = [iv.start for iv in intervals]
        ends = [iv.end for iv in intervals]
        assert 200 in starts and 500 in starts
        assert 400 in ends and 700 in ends

    def test_multiple_transcripts_complex(self):
        """Two multi-exon transcripts with overlapping introns.

        t0: exon 100-200, exon 500-600 (intron 200-500)
        t1: exon 150-300, exon 450-550 (intron 300-450)
        Global exon union: [(100, 300), (450, 600)]
        t0 intron 200-500 minus global exons → 300-450
        t1 intron 300-450 minus global exons → 300-450
        """
        t0 = _make_transcript("chr1", Strand.POS,
                              [(100, 200), (500, 600)], t_index=0)
        t1 = _make_transcript("chr1", Strand.POS,
                              [(150, 300), (450, 550)], t_index=1)
        intervals = list(_gen_cluster_unambig_intron_intervals([t0, t1]))
        # t0: intron 200-500, minus global exons [(100,300),(450,600)] → 300-450
        # t1: intron 300-450, minus global exons [(100,300),(450,600)] → 300-450
        assert len(intervals) == 2
        for iv in intervals:
            assert iv.start == 300
            assert iv.end == 450

    def test_empty_cluster(self):
        """Empty cluster → no intervals."""
        intervals = list(_gen_cluster_unambig_intron_intervals([]))
        assert len(intervals) == 0


# =====================================================================
# Full index build with UNAMBIG_INTRON
# =====================================================================


class TestBuildGenomicIntervalsWithUnambigIntron:
    def test_unambig_intron_in_built_intervals(self):
        """build_genomic_intervals includes UNAMBIG_INTRON for multi-exon."""
        t0 = _make_transcript("chr1", Strand.POS,
                              [(100, 200), (400, 500)], t_index=0)
        ref_lengths = {"chr1": 1000}
        df = build_genomic_intervals([t0], ref_lengths)
        unambig = df[df["interval_type"] == int(IntervalType.UNAMBIG_INTRON)]
        assert len(unambig) == 1
        row = unambig.iloc[0]
        assert row["start"] == 200
        assert row["end"] == 400
        assert row["t_index"] == 0

    def test_no_unambig_intron_when_fully_covered(self):
        """No UNAMBIG_INTRON when intron is covered by another exon."""
        t0 = _make_transcript("chr1", Strand.POS,
                              [(100, 200), (400, 500)], t_index=0)
        t1 = _make_transcript("chr1", Strand.POS,
                              [(100, 500)], t_index=1)
        ref_lengths = {"chr1": 1000}
        df = build_genomic_intervals([t0, t1], ref_lengths)
        unambig = df[df["interval_type"] == int(IntervalType.UNAMBIG_INTRON)]
        assert len(unambig) == 0


# =====================================================================
# End-to-end regression: zero nRNA in pristine conditions
# =====================================================================


class TestNRNAPhantomRegression:
    """End-to-end tests that zero-nRNA scenarios produce near-zero nRNA.

    These are the critical regression tests for the phantom nRNA bug.
    A scenario with overlapping genes (the most problematic case for
    ambiguous intron quantification) should not produce phantom nRNA counts
    when no nascent RNA is simulated.
    """

    def test_overlapping_genes_zero_nrna(self, tmp_path):
        """Overlapping genes, nRNA=0 → pipeline nRNA should be near zero.

        This is the exact scenario that exposed the phantom nRNA bug:
        - g1 (+): 2-exon transcript (has intron overlapping g2's exon)
        - g2 (-): 1-exon transcript spanning g1's intron
        With nRNA=0, the pipeline should NOT assign fragments to nRNA pool.
        """
        sc = Scenario(
            "nrna_phantom_regression", genome_length=5000,
            seed=42, work_dir=tmp_path / "nrna_phantom",
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1500, 1800)],
             "abundance": 100},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(400, 1600)],
             "abundance": 50},
        ])
        try:
            sim_cfg = SimConfig(
                frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                read_length=100, strand_specificity=1.0, seed=42,
            )
            result = sc.build_oracle(
                n_fragments=500, sim_config=sim_cfg,
                nrna_abundance=0,
            )
            config = PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
            )
            pr = run_pipeline(result.bam_path, result.index, config=config)
            bench = run_benchmark(result, pr,
                                  scenario_name="nrna_phantom_regression")
            # With nRNA=0 truth, pipeline should produce very little nRNA
            # Allow up to 3% of total fragments as tolerance
            max_phantom = 0.03 * bench.n_fragments
            assert bench.n_nrna_pipeline <= max_phantom, (
                f"Phantom nRNA detected: {bench.n_nrna_pipeline:.0f} "
                f"(limit={max_phantom:.0f}, "
                f"n_fragments={bench.n_fragments})"
            )
        finally:
            sc.cleanup()

    def test_dense_overlapping_genes_zero_nrna(self, tmp_path):
        """Densely overlapping gene cluster, nRNA=0 → near-zero nRNA.

        Three genes with extensive exon/intron overlap — stress tests
        the UNAMBIG_INTRON computation.
        """
        sc = Scenario(
            "nrna_dense_regression", genome_length=8000,
            seed=42, work_dir=tmp_path / "nrna_dense",
        )
        # g1 (+): long 3-exon transcript with wide introns
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(200, 500), (1500, 1800), (3000, 3300)],
             "abundance": 100},
        ])
        # g2 (-): single-exon spanning g1's first intron
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(400, 1600)],
             "abundance": 50},
        ])
        # g3 (+): 2-exon transcript overlapping g1's second intron
        sc.add_gene("g3", "+", [
            {"t_id": "t3",
             "exons": [(1700, 2000), (2500, 3100)],
             "abundance": 30},
        ])
        try:
            sim_cfg = SimConfig(
                frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                read_length=100, strand_specificity=1.0, seed=42,
            )
            result = sc.build_oracle(
                n_fragments=800, sim_config=sim_cfg,
                nrna_abundance=0,
            )
            config = PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
            )
            pr = run_pipeline(result.bam_path, result.index, config=config)
            bench = run_benchmark(result, pr,
                                  scenario_name="nrna_dense_regression")
            max_phantom = 0.05 * bench.n_fragments
            assert bench.n_nrna_pipeline <= max_phantom, (
                f"Phantom nRNA detected: {bench.n_nrna_pipeline:.0f} "
                f"(limit={max_phantom:.0f}, "
                f"n_fragments={bench.n_fragments})"
            )
        finally:
            sc.cleanup()

    def test_non_overlapping_multi_exon_zero_nrna(self, tmp_path):
        """Non-overlapping multi-exon genes with nRNA=0.

        Even with unambiguous introns (all introns are unambiguous when
        no overlapping exons exist), pipeline should not call nRNA when
        the truth is zero nRNA.
        """
        sc = Scenario(
            "nrna_non_overlap_regression", genome_length=5000,
            seed=42, work_dir=tmp_path / "nrna_non_overlap",
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(200, 400), (800, 1000)],
             "abundance": 100},
        ])
        sc.add_gene("g2", "+", [
            {"t_id": "t2",
             "exons": [(2000, 2200), (2600, 2800)],
             "abundance": 80},
        ])
        try:
            sim_cfg = SimConfig(
                frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                read_length=100, strand_specificity=1.0, seed=42,
            )
            result = sc.build_oracle(
                n_fragments=500, sim_config=sim_cfg,
                nrna_abundance=0,
            )
            config = PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
            )
            pr = run_pipeline(result.bam_path, result.index, config=config)
            bench = run_benchmark(result, pr,
                                  scenario_name="nrna_non_overlap")
            max_phantom = 0.03 * bench.n_fragments
            assert bench.n_nrna_pipeline <= max_phantom, (
                f"Phantom nRNA detected: {bench.n_nrna_pipeline:.0f} "
                f"(limit={max_phantom:.0f}, "
                f"n_fragments={bench.n_fragments})"
            )
        finally:
            sc.cleanup()
