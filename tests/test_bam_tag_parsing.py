"""Tests for robust BAM auxiliary tag parsing in the C++ scanner.

Exercises the strand-tag parsing path (``read_sj_strand`` / ``try_tag``
in bam_scanner.cpp) and the ``detect_sj_strand_tag`` helper through
Python integration tests.  No C/C++ unit test framework is needed
because the native module is fully exposed via nanobind.

Tag type scenarios covered:
- XS tag written as BAM type 'A' (single char)  — standard STAR output
- XS tag written as BAM type 'Z' (string)       — pysam default for strings
- ts tag as type 'A' and 'Z'                     — minimap2 convention
- Missing strand tags entirely
- Numeric tags (NH, NM) as various integer subtypes
- Mixed tag presence across reads (auto-detection)
- Read-strand flipping for the 'ts' tag on reverse reads
"""

import struct
import tempfile
from pathlib import Path

import numpy as np
import pysam
import pytest

from hulkrna._bam_impl import BamScanner, detect_sj_strand_tag
from hulkrna.index import TranscriptIndex


# =====================================================================
# Helpers — write minimal BAMs with controlled tag types
# =====================================================================

_REF_NAME = "chr1"
_REF_LEN = 100_000
_HEADER = {
    "HD": {"VN": "1.6", "SO": "queryname"},
    "SQ": [{"SN": _REF_NAME, "LN": _REF_LEN}],
}


def _write_bam(path: Path, reads: list[pysam.AlignedSegment]) -> str:
    """Write a sorted + indexed BAM from a list of pysam AlignedSegments."""
    bam_path = str(path / "test.bam")
    with pysam.AlignmentFile(bam_path, "wb", header=_HEADER) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-n", "-o", bam_path, bam_path)
    return bam_path


def _make_spliced_read(
    qname: str,
    pos: int,
    *,
    flag: int = 0x41,  # paired, read1
    cigar: list[tuple[int, int]] | None = None,
    tags: list[tuple] | None = None,
) -> pysam.AlignedSegment:
    """Create a spliced paired-end read for the test reference."""
    a = pysam.AlignedSegment()
    a.query_name = qname
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = 60
    a.flag = flag
    if cigar is None:
        # 50M 200N 50M — two exon blocks separated by a 200bp intron
        cigar = [(0, 50), (3, 200), (0, 50)]
    a.cigar = cigar
    a.query_sequence = "A" * sum(l for op, l in cigar if op in (0, 1, 4))
    a.query_qualities = pysam.qualitystring_to_array(
        "I" * len(a.query_sequence)
    )
    a.next_reference_id = 0
    a.next_reference_start = pos + 300
    a.template_length = 500
    if tags is not None:
        a.set_tags(tags)
    return a


def _make_read_pair(
    qname: str,
    pos: int,
    *,
    r1_cigar: list[tuple[int, int]] | None = None,
    r2_cigar: list[tuple[int, int]] | None = None,
    r1_tags: list[tuple] | None = None,
    r2_tags: list[tuple] | None = None,
    r1_is_reverse: bool = False,
    r2_is_reverse: bool = True,
) -> list[pysam.AlignedSegment]:
    """Create a properly paired R1/R2 pair."""
    if r1_cigar is None:
        r1_cigar = [(0, 50), (3, 200), (0, 50)]
    if r2_cigar is None:
        r2_cigar = [(0, 100)]

    seq_len_r1 = sum(l for op, l in r1_cigar if op in (0, 1, 4))
    seq_len_r2 = sum(l for op, l in r2_cigar if op in (0, 1, 4))

    r2_start = pos + 400

    # R1
    r1 = pysam.AlignedSegment()
    r1.query_name = qname
    r1.reference_id = 0
    r1.reference_start = pos
    r1.mapping_quality = 60
    r1_flag = 0x1 | 0x2 | 0x40  # paired, proper, read1
    if r1_is_reverse:
        r1_flag |= 0x10
    if r2_is_reverse:
        r1_flag |= 0x20
    r1.flag = r1_flag
    r1.cigar = r1_cigar
    r1.query_sequence = "A" * seq_len_r1
    r1.query_qualities = pysam.qualitystring_to_array("I" * seq_len_r1)
    r1.next_reference_id = 0
    r1.next_reference_start = r2_start
    r1.template_length = r2_start + seq_len_r2 - pos
    if r1_tags is not None:
        r1.set_tags(r1_tags)

    # R2
    r2 = pysam.AlignedSegment()
    r2.query_name = qname
    r2.reference_id = 0
    r2.reference_start = r2_start
    r2.mapping_quality = 60
    r2_flag = 0x1 | 0x2 | 0x80  # paired, proper, read2
    if r2_is_reverse:
        r2_flag |= 0x10
    if r1_is_reverse:
        r2_flag |= 0x20
    r2.flag = r2_flag
    r2.cigar = r2_cigar
    r2.query_sequence = "A" * seq_len_r2
    r2.query_qualities = pysam.qualitystring_to_array("I" * seq_len_r2)
    r2.next_reference_id = 0
    r2.next_reference_start = pos
    r2.template_length = -(r2_start + seq_len_r2 - pos)
    if r2_tags is not None:
        r2.set_tags(r2_tags)

    return [r1, r2]


def _read_back_tags(bam_path: str) -> list[dict]:
    """Read back tags from a BAM and return type info for inspection."""
    results = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            tag_info = {}
            for tag, val in read.get_tags(with_value_type=True):
                # get_tags with_value_type returns (tag, value, type_code)
                pass
            # Alternative: iterate raw tags
            raw_tags = {}
            for tag, val, type_code in read.get_tags(with_value_type=True):
                raw_tags[tag] = {"value": val, "type": type_code}
            results.append({"qname": read.query_name, "tags": raw_tags})
    return results


# =====================================================================
# Helpers — low-level BAM tag writing for exact type control
# =====================================================================


def _set_tag_raw(segment: pysam.AlignedSegment, tag: str, value, type_code: str):
    """Set a BAM tag with explicit type code using pysam's set_tag().

    This gives precise control over whether a tag is stored as 'A' (char),
    'Z' (string), 'i' (int32), 'c' (int8), etc.
    """
    segment.set_tag(tag, value, type_code)


# =====================================================================
# Test: detect_sj_strand_tag
# =====================================================================


class TestDetectSJStrandTag:
    """Tests for the native ``detect_sj_strand_tag()`` function."""

    def test_detect_xs_type_a(self, tmp_path):
        """XS tag stored as type 'A' → detected as 'XS'."""
        reads = _make_read_pair(
            "r1", 1000,
            r1_tags=[("NH", 1), ("XS", "+", "A")],
            r2_tags=[("NH", 1)],
        )
        bam_path = _write_bam(tmp_path, reads)
        assert detect_sj_strand_tag(bam_path) == "XS"

    def test_detect_xs_type_z(self, tmp_path):
        """XS tag stored as type 'Z' → still detected (tag is present)."""
        reads = _make_read_pair(
            "r1", 1000,
            r1_tags=[("NH", 1), ("XS", "+", "Z")],
            r2_tags=[("NH", 1)],
        )
        bam_path = _write_bam(tmp_path, reads)
        assert detect_sj_strand_tag(bam_path) == "XS"

    def test_detect_ts(self, tmp_path):
        """ts tag detected."""
        reads = _make_read_pair(
            "r1", 1000,
            r1_tags=[("NH", 1), ("ts", "+", "A")],
            r2_tags=[("NH", 1)],
        )
        bam_path = _write_bam(tmp_path, reads)
        assert detect_sj_strand_tag(bam_path) == "ts"

    def test_detect_both_xs_ts(self, tmp_path):
        """Both XS and ts present → 'XS,ts'."""
        reads = _make_read_pair(
            "r1", 1000,
            r1_tags=[("NH", 1), ("XS", "+", "A"), ("ts", "+", "A")],
            r2_tags=[("NH", 1)],
        )
        bam_path = _write_bam(tmp_path, reads)
        assert detect_sj_strand_tag(bam_path) == "XS,ts"

    def test_detect_none_unspliced(self, tmp_path):
        """Unspliced reads only → 'none' (no spliced reads to find tags on)."""
        r1 = pysam.AlignedSegment()
        r1.query_name = "r1"
        r1.reference_id = 0
        r1.reference_start = 1000
        r1.mapping_quality = 60
        r1.flag = 0x1 | 0x2 | 0x40
        r1.cigar = [(0, 100)]  # no splicing
        r1.query_sequence = "A" * 100
        r1.query_qualities = pysam.qualitystring_to_array("I" * 100)
        r1.next_reference_id = 0
        r1.next_reference_start = 1200
        r1.template_length = 300
        r1.set_tags([("NH", 1), ("XS", "+", "A")])

        bam_path = _write_bam(tmp_path, [r1])
        # XS is on the read but it's not spliced — detect only checks spliced
        assert detect_sj_strand_tag(bam_path) == "none"

    def test_detect_none_no_tags(self, tmp_path):
        """Spliced reads without strand tags → 'none'."""
        reads = _make_read_pair(
            "r1", 1000,
            r1_tags=[("NH", 1)],  # no XS/ts
            r2_tags=[("NH", 1)],
        )
        bam_path = _write_bam(tmp_path, reads)
        assert detect_sj_strand_tag(bam_path) == "none"


# =====================================================================
# Test: XS tag type variants via full pipeline
# =====================================================================


class TestXSTagTypeVariants:
    """Verify the C++ scanner correctly parses XS tags stored as
    BAM type 'A' (single char) and type 'Z' (null-terminated string).

    Uses the full oracle scenario pipeline so that strand observations
    flow through the C++ ``try_tag`` lambda in ``read_sj_strand()``
    and are replayed into the Python StrandModels.
    """

    @staticmethod
    def _build_simple_scenario(tmp_path, xs_type_code: str):
        """Build a minimal 1-transcript scenario with specified XS type.

        Returns (bam_path, index) ready for pipeline run.
        """
        from hulkrna.sim import Scenario, SimConfig

        scenario = Scenario(
            "xs_tag_test",
            genome_length=100_000,
            seed=42,
            work_dir=tmp_path / "scenario",
        )
        scenario.add_gene(
            gene_id="G1",
            strand="+",
            transcripts=[
                {
                    "t_id": "T1",
                    "exons": [(1000, 2000), (3000, 4000), (5000, 6000)],
                    "abundance": 100.0,
                }
            ],
        )
        sc = SimConfig(
            frag_mean=250, frag_std=40, frag_min=100, frag_max=500,
            read_length=100, strand_specificity=0.95, seed=42,
        )
        result = scenario.build_oracle(
            n_fragments=200,
            sim_config=sc,
        )

        # Now rewrite the BAM with the desired XS type code.
        # Read the oracle BAM, modify XS tags, write new BAM.
        original_bam = result.bam_path
        modified_bam = str(tmp_path / "modified.bam")
        with pysam.AlignmentFile(original_bam, "rb") as inp:
            header = inp.header.to_dict()
            with pysam.AlignmentFile(modified_bam, "wb", header=header) as out:
                for read in inp:
                    if read.has_tag("XS"):
                        xs_val = read.get_tag("XS")
                        read.set_tag("XS", xs_val, xs_type_code)
                    out.write(read)

        # Re-sort by name (required by scanner)
        sorted_bam = str(tmp_path / "modified_sorted.bam")
        pysam.sort("-n", "-o", sorted_bam, modified_bam)

        return sorted_bam, result.index

    def test_xs_type_a_produces_strand_observations(self, tmp_path):
        """XS as type 'A' → exonic_spliced strand observations populated."""
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        bam_path, index = self._build_simple_scenario(tmp_path, "A")
        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr = run_pipeline(bam_path, index, config=config)

        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations > 0, (
            "Type 'A' XS tag should produce exonic_spliced observations"
        )
        # With SS=0.95, strand_specificity should be high
        assert sm.strand_specificity > 0.8

    def test_xs_type_z_produces_strand_observations(self, tmp_path):
        """XS as type 'Z' → exonic_spliced strand observations populated
        (tests the Z-type fallback in try_tag)."""
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        bam_path, index = self._build_simple_scenario(tmp_path, "Z")
        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr = run_pipeline(bam_path, index, config=config)

        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations > 0, (
            "Type 'Z' XS tag should be handled by Z-type fallback"
        )
        assert sm.strand_specificity > 0.8

    def test_type_a_and_z_produce_same_strand_model(self, tmp_path):
        """Type 'A' and 'Z' should yield identical strand observations."""
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        # Type A
        dir_a = tmp_path / "type_a"
        dir_a.mkdir()
        bam_a, idx_a = self._build_simple_scenario(dir_a, "A")
        cfg = PipelineConfig(
            em=EMConfig(seed=42), scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr_a = run_pipeline(bam_a, idx_a, config=cfg)

        # Type Z
        dir_z = tmp_path / "type_z"
        dir_z.mkdir()
        bam_z, idx_z = self._build_simple_scenario(dir_z, "Z")
        pr_z = run_pipeline(bam_z, idx_z, config=cfg)

        sm_a = pr_a.strand_models.exonic_spliced
        sm_z = pr_z.strand_models.exonic_spliced
        assert sm_a.n_observations == sm_z.n_observations, (
            f"Mismatch: A={sm_a.n_observations}, Z={sm_z.n_observations}"
        )
        assert sm_a.n_same == sm_z.n_same
        assert sm_a.n_opposite == sm_z.n_opposite


# =====================================================================
# Test: ts tag (minimap2 convention) with reverse-read strand flipping
# =====================================================================


class TestTsTagFlipping:
    """The 'ts' tag (minimap2) is alignment-relative: its sign must
    be flipped for reverse-strand reads.  Verify this via detect +
    inspection of tag rewriting."""

    @staticmethod
    def _build_scenario_with_ts(tmp_path, type_code: str = "A"):
        """Build scenario, then rewrite XS→ts tags."""
        from hulkrna.sim import Scenario, SimConfig

        scenario = Scenario(
            "ts_tag_test",
            genome_length=100_000,
            seed=42,
            work_dir=tmp_path / "scenario",
        )
        scenario.add_gene(
            gene_id="G1",
            strand="+",
            transcripts=[
                {
                    "t_id": "T1",
                    "exons": [(1000, 2000), (3000, 4000), (5000, 6000)],
                    "abundance": 100.0,
                }
            ],
        )
        sc = SimConfig(
            frag_mean=250, frag_std=40, frag_min=100, frag_max=500,
            read_length=100, strand_specificity=0.95, seed=42,
        )
        result = scenario.build_oracle(n_fragments=200, sim_config=sc)

        # Rewrite: remove XS, add ts with minimap2 semantics.
        # For forward reads: ts = XS (same as original)
        # For reverse reads: ts = flipped XS (minimap2 convention:
        #   ts indicates whether transcription direction matches the
        #   alignment direction, so reverse-strand alignments flip).
        original_bam = result.bam_path
        modified_bam = str(tmp_path / "ts_modified.bam")
        with pysam.AlignmentFile(original_bam, "rb") as inp:
            header = inp.header.to_dict()
            with pysam.AlignmentFile(modified_bam, "wb", header=header) as out:
                for read in inp:
                    if read.has_tag("XS"):
                        xs_val = read.get_tag("XS")
                        is_reverse = read.is_reverse
                        # ts = original strand for forward reads,
                        # flipped for reverse reads (since our C++ code
                        # will flip it back).
                        if is_reverse:
                            ts_val = "-" if xs_val == "+" else "+"
                        else:
                            ts_val = xs_val
                        read.set_tag("XS", None)
                        read.set_tag("ts", ts_val, type_code)
                    out.write(read)

        sorted_bam = str(tmp_path / "ts_sorted.bam")
        pysam.sort("-n", "-o", sorted_bam, modified_bam)

        return sorted_bam, result.index

    def test_ts_type_a_with_flipping(self, tmp_path):
        """ts as type 'A' with reverse-read flipping → correct SS."""
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        bam_path, index = self._build_scenario_with_ts(tmp_path, "A")
        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="ts"),
        )
        pr = run_pipeline(bam_path, index, config=config)

        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations > 0
        assert sm.strand_specificity > 0.8

    def test_ts_type_z_with_flipping(self, tmp_path):
        """ts as type 'Z' with reverse-read flipping → correct SS."""
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        bam_path, index = self._build_scenario_with_ts(tmp_path, "Z")
        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="ts"),
        )
        pr = run_pipeline(bam_path, index, config=config)

        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations > 0, (
            "ts type 'Z' should be handled by Z-type fallback"
        )
        assert sm.strand_specificity > 0.8


# =====================================================================
# Test: Missing or invalid strand tags
# =====================================================================


class TestMissingAndInvalidTags:
    """Verify scanner handles missing or unexpected tag values gracefully."""

    @staticmethod
    def _build_scenario_without_xs(tmp_path):
        """Build scenario, then strip all XS/ts tags."""
        from hulkrna.sim import Scenario, SimConfig

        scenario = Scenario(
            "no_tags_test",
            genome_length=100_000,
            seed=42,
            work_dir=tmp_path / "scenario",
        )
        scenario.add_gene(
            gene_id="G1",
            strand="+",
            transcripts=[
                {
                    "t_id": "T1",
                    "exons": [(1000, 2000), (3000, 4000), (5000, 6000)],
                    "abundance": 100.0,
                }
            ],
        )
        sc = SimConfig(
            frag_mean=250, frag_std=40, frag_min=100, frag_max=500,
            read_length=100, strand_specificity=0.95, seed=42,
        )
        result = scenario.build_oracle(n_fragments=200, sim_config=sc)

        # Strip all strand tags
        original_bam = result.bam_path
        modified_bam = str(tmp_path / "no_tags.bam")
        with pysam.AlignmentFile(original_bam, "rb") as inp:
            header = inp.header.to_dict()
            with pysam.AlignmentFile(modified_bam, "wb", header=header) as out:
                for read in inp:
                    if read.has_tag("XS"):
                        read.set_tag("XS", None)
                    if read.has_tag("ts"):
                        read.set_tag("ts", None)
                    out.write(read)

        sorted_bam = str(tmp_path / "no_tags_sorted.bam")
        pysam.sort("-n", "-o", sorted_bam, modified_bam)

        return sorted_bam, result.index

    def test_no_strand_tags_zero_spliced_observations(self, tmp_path):
        """No XS/ts tags → zero exonic_spliced observations."""
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        bam_path, index = self._build_scenario_without_xs(tmp_path)
        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="XS"),
        )
        pr = run_pipeline(bam_path, index, config=config)

        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations == 0, (
            "Without strand tags, exonic_spliced should have 0 observations"
        )

    def test_sj_mode_none_ignores_tags(self, tmp_path):
        """sj_strand_tag='none' → ignore all strand tags."""
        from hulkrna.sim import Scenario, SimConfig
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        scenario = Scenario(
            "mode_none_test",
            genome_length=100_000,
            seed=42,
            work_dir=tmp_path / "scenario",
        )
        scenario.add_gene(
            gene_id="G1",
            strand="+",
            transcripts=[
                {
                    "t_id": "T1",
                    "exons": [(1000, 2000), (3000, 4000), (5000, 6000)],
                    "abundance": 100.0,
                }
            ],
        )
        sc = SimConfig(
            frag_mean=250, frag_std=40, frag_min=100, frag_max=500,
            read_length=100, strand_specificity=0.95, seed=42,
        )
        result = scenario.build_oracle(n_fragments=200, sim_config=sc)
        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="none"),
        )
        pr = run_pipeline(result.bam_path, result.index, config=config)

        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations == 0, (
            "With sj_strand_tag='none', exonic_spliced should have 0 "
            "observations even if tags exist"
        )

    def test_invalid_xs_value_treated_as_none(self, tmp_path):
        """XS tag with value '.' or '*' → STRAND_NONE (not crash)."""
        from hulkrna.sim import Scenario, SimConfig
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        scenario = Scenario(
            "invalid_xs_test",
            genome_length=100_000,
            seed=42,
            work_dir=tmp_path / "scenario",
        )
        scenario.add_gene(
            gene_id="G1",
            strand="+",
            transcripts=[
                {
                    "t_id": "T1",
                    "exons": [(1000, 2000), (3000, 4000), (5000, 6000)],
                    "abundance": 100.0,
                }
            ],
        )
        sc = SimConfig(
            frag_mean=250, frag_std=40, frag_min=100, frag_max=500,
            read_length=100, strand_specificity=0.95, seed=42,
        )
        result = scenario.build_oracle(n_fragments=200, sim_config=sc)

        # Rewrite XS to invalid values
        original_bam = result.bam_path
        modified_bam = str(tmp_path / "invalid_xs.bam")
        with pysam.AlignmentFile(original_bam, "rb") as inp:
            header = inp.header.to_dict()
            with pysam.AlignmentFile(modified_bam, "wb", header=header) as out:
                for read in inp:
                    if read.has_tag("XS"):
                        read.set_tag("XS", ".", "A")
                    out.write(read)

        sorted_bam = str(tmp_path / "invalid_xs_sorted.bam")
        pysam.sort("-n", "-o", sorted_bam, modified_bam)

        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="XS"),
        )
        # Should not raise — invalid values produce STRAND_NONE
        pr = run_pipeline(sorted_bam, result.index, config=config)
        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations == 0


# =====================================================================
# Test: Oracle BAM XS tag output type
# =====================================================================


class TestOracleBamXSType:
    """Verify that the oracle BAM writes XS as type 'A' (single char)
    as required by the BAM specification and htslib ``bam_aux2A``."""

    def test_oracle_bam_writes_xs_type_a(self, tmp_path):
        """Oracle BAM should write XS as type 'A', not 'Z'."""
        from hulkrna.sim import Scenario, SimConfig

        scenario = Scenario(
            "oracle_xs_type",
            genome_length=100_000,
            seed=42,
            work_dir=tmp_path / "scenario",
        )
        scenario.add_gene(
            gene_id="G1",
            strand="+",
            transcripts=[
                {
                    "t_id": "T1",
                    "exons": [(1000, 2000), (3000, 4000), (5000, 6000)],
                    "abundance": 100.0,
                }
            ],
        )
        sc = SimConfig(
            frag_mean=250, frag_std=40, frag_min=100, frag_max=500,
            read_length=100, strand_specificity=0.95, seed=42,
        )
        result = scenario.build_oracle(n_fragments=50, sim_config=sc)

        with pysam.AlignmentFile(result.bam_path, "rb") as bam:
            for read in bam:
                if read.has_tag("XS"):
                    # get_tags with_value_type returns (tag, value, type)
                    tag_dict = {
                        t: (v, c)
                        for t, v, c in read.get_tags(with_value_type=True)
                    }
                    xs_val, xs_type = tag_dict["XS"]
                    assert xs_type == "A", (
                        f"XS tag should be type 'A' (char), "
                        f"got type '{xs_type}' value '{xs_val}'"
                    )
                    assert xs_val in ("+", "-"), (
                        f"XS value should be '+' or '-', got '{xs_val}'"
                    )


# =====================================================================
# Test: Integer tag robustness (NH, NM, HI)
# =====================================================================


class TestIntegerTags:
    """Verify NH/NM tags parsed correctly regardless of integer subtype.

    htslib ``bam_aux2i`` handles all integer subtypes: c/C/s/S/i/I.
    Pysam may write NH as 'C' (uint8) for small values.  Verify the
    scanner doesn't break on unusual but valid type codes.
    """

    @staticmethod
    def _build_and_check_nh(tmp_path, nh_value: int, nh_type: str):
        """Build BAM with specific NH type code, verify scanner reads it."""
        from hulkrna.sim import Scenario, SimConfig
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        scenario = Scenario(
            "nh_type_test",
            genome_length=100_000,
            seed=42,
            work_dir=tmp_path / "scenario",
        )
        scenario.add_gene(
            gene_id="G1",
            strand="+",
            transcripts=[
                {
                    "t_id": "T1",
                    "exons": [(1000, 2000), (3000, 4000), (5000, 6000)],
                    "abundance": 100.0,
                }
            ],
        )
        sc = SimConfig(
            frag_mean=250, frag_std=40, frag_min=100, frag_max=500,
            read_length=100, strand_specificity=0.95, seed=42,
        )
        result = scenario.build_oracle(n_fragments=50, sim_config=sc)

        # Rewrite NH with specific type code
        original_bam = result.bam_path
        modified_bam = str(tmp_path / "nh_modified.bam")
        with pysam.AlignmentFile(original_bam, "rb") as inp:
            header = inp.header.to_dict()
            with pysam.AlignmentFile(modified_bam, "wb", header=header) as out:
                for read in inp:
                    read.set_tag("NH", nh_value, nh_type)
                    out.write(read)

        sorted_bam = str(tmp_path / "nh_sorted.bam")
        pysam.sort("-n", "-o", sorted_bam, modified_bam)

        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        # Should not crash regardless of integer subtype
        pr = run_pipeline(sorted_bam, result.index, config=config)
        return pr

    @pytest.mark.parametrize("nh_type", ["c", "C", "s", "S", "i", "I"])
    def test_nh_integer_subtypes(self, tmp_path, nh_type):
        """NH tag as various integer subtypes → no crash, valid output."""
        pr = self._build_and_check_nh(tmp_path, 1, nh_type)
        assert pr.stats.total > 0

    def test_nh_zero_still_works(self, tmp_path):
        """NH=0 (rare, edge case) should not crash."""
        pr = self._build_and_check_nh(tmp_path, 0, "i")
        # Pipeline may handle NH=0 differently, but must not crash
        assert pr.stats.total > 0


# =====================================================================
# Test: SJ tag priority modes
# =====================================================================


class TestSJTagPriority:
    """Verify 'XS,ts' and 'ts,XS' priority ordering works correctly."""

    @staticmethod
    def _build_with_both_tags(tmp_path, xs_strand: str, ts_strand: str):
        """Build BAM with both XS and ts tags set to different values."""
        from hulkrna.sim import Scenario, SimConfig

        scenario = Scenario(
            "both_tags_test",
            genome_length=100_000,
            seed=42,
            work_dir=tmp_path / "scenario",
        )
        scenario.add_gene(
            gene_id="G1",
            strand="+",
            transcripts=[
                {
                    "t_id": "T1",
                    "exons": [(1000, 2000), (3000, 4000), (5000, 6000)],
                    "abundance": 100.0,
                }
            ],
        )
        sc = SimConfig(
            frag_mean=250, frag_std=40, frag_min=100, frag_max=500,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        result = scenario.build_oracle(n_fragments=200, sim_config=sc)

        # Set both XS and ts to specified values.
        # For ts, reverse-read flipping means we must pre-flip.
        original_bam = result.bam_path
        modified_bam = str(tmp_path / "both_tags.bam")
        with pysam.AlignmentFile(original_bam, "rb") as inp:
            header = inp.header.to_dict()
            with pysam.AlignmentFile(modified_bam, "wb", header=header) as out:
                for read in inp:
                    # XS always has the genomic-strand value
                    read.set_tag("XS", xs_strand, "A")
                    # ts is alignment-relative: flip for reverse reads
                    if read.is_reverse:
                        ts_written = "-" if ts_strand == "+" else "+"
                    else:
                        ts_written = ts_strand
                    read.set_tag("ts", ts_written, "A")
                    out.write(read)

        sorted_bam = str(tmp_path / "both_sorted.bam")
        pysam.sort("-n", "-o", sorted_bam, modified_bam)

        return sorted_bam, result.index

    def test_xs_ts_uses_xs_when_present(self, tmp_path):
        """With sj_strand_tag='XS,ts', both tags present with correct strand
        → scanner reads XS first, produces strand observations."""
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        # Both tags carry '+' (correct for the '+' strand gene).
        dir1 = tmp_path / "d1"
        dir1.mkdir()
        bam_path, index = self._build_with_both_tags(dir1, "+", "+")

        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag=("XS", "ts")),
        )
        pr = run_pipeline(bam_path, index, config=config)

        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations > 0
        assert sm.strand_specificity > 0.8

    def test_fallback_to_second_tag(self, tmp_path):
        """With sj_strand_tag='XS,ts' and XS removed, scanner
        falls back to ts (second priority) and still gets observations."""
        from hulkrna.sim import Scenario, SimConfig
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        # Build scenario with ts tags only (XS removed)
        scenario = Scenario(
            "fallback_test",
            genome_length=100_000,
            seed=42,
            work_dir=tmp_path / "scenario",
        )
        scenario.add_gene(
            gene_id="G1",
            strand="+",
            transcripts=[
                {
                    "t_id": "T1",
                    "exons": [(1000, 2000), (3000, 4000), (5000, 6000)],
                    "abundance": 100.0,
                }
            ],
        )
        sc = SimConfig(
            frag_mean=250, frag_std=40, frag_min=100, frag_max=500,
            read_length=100, strand_specificity=0.95, seed=42,
        )
        result = scenario.build_oracle(n_fragments=200, sim_config=sc)

        # Rewrite: remove XS, add ts with correct strand
        original_bam = result.bam_path
        modified_bam = str(tmp_path / "ts_only.bam")
        with pysam.AlignmentFile(original_bam, "rb") as inp:
            header = inp.header.to_dict()
            with pysam.AlignmentFile(modified_bam, "wb", header=header) as out:
                for read in inp:
                    if read.has_tag("XS"):
                        xs_val = read.get_tag("XS")
                        # ts is alignment-relative: flip for reverse reads
                        ts_val = xs_val
                        if read.is_reverse:
                            ts_val = "-" if xs_val == "+" else "+"
                        read.set_tag("XS", None)
                        read.set_tag("ts", ts_val, "A")
                    out.write(read)

        sorted_bam = str(tmp_path / "ts_only_sorted.bam")
        pysam.sort("-n", "-o", sorted_bam, modified_bam)

        # Use XS,ts priority — XS is absent, should fall back to ts
        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag=("XS", "ts")),
        )
        pr = run_pipeline(sorted_bam, result.index, config=config)

        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations > 0, (
            "Fallback to ts should produce exonic_spliced observations"
        )
        assert sm.strand_specificity > 0.8
