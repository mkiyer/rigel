"""The C++ scanner reads BAM auxiliary tags the way real aligners write them.

One aligner writes ``XS`` as a single character and another as a string; ``ts`` follows a different
convention again and must be flipped on a reverse read; ``NH`` and ``NM`` arrive as any of several
integer subtypes. The same library tagged two ways must produce the same strand model; ``auto`` reads
both tags on every spliced record, an agreement giving the strand and a disagreement an ambiguous
motif, counted; an ordered spec takes the first tag present; an unparseable value must read as absent
rather than as a strand, and a library with no strand tags at all must yield zero spliced
observations rather than invent one. Driven through Python, because the native module is fully
exposed and needs no C++ harness.
"""

import pysam
import pytest


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
        from rigel.sim import Scenario, ReadSimConfig

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
        sc = ReadSimConfig(
            frag_mean=250,
            frag_std=40,
            frag_min=100,
            frag_max=500,
            read_length=100,
            strand_specificity=0.95,
            seed=42,
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
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        bam_path, index = self._build_simple_scenario(tmp_path, "A")
        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr = run_pipeline(bam_path, index, config=config)

        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations > 0, "Type 'A' XS tag should produce exonic_spliced observations"
        # With SS=0.95, strand_specificity should be high
        assert sm.strand_specificity > 0.8

    def test_xs_type_z_produces_strand_observations(self, tmp_path):
        """XS as type 'Z' → exonic_spliced strand observations populated
        (tests the Z-type fallback in try_tag)."""
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        bam_path, index = self._build_simple_scenario(tmp_path, "Z")
        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr = run_pipeline(bam_path, index, config=config)

        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations > 0, "Type 'Z' XS tag should be handled by Z-type fallback"
        assert sm.strand_specificity > 0.8

    def test_type_a_and_z_produce_same_strand_model(self, tmp_path):
        """Type 'A' and 'Z' should yield identical strand observations."""
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        # Type A
        dir_a = tmp_path / "type_a"
        dir_a.mkdir()
        bam_a, idx_a = self._build_simple_scenario(dir_a, "A")
        cfg = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="auto"),
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
        from rigel.sim import Scenario, ReadSimConfig

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
        sc = ReadSimConfig(
            frag_mean=250,
            frag_std=40,
            frag_min=100,
            frag_max=500,
            read_length=100,
            strand_specificity=0.95,
            seed=42,
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
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

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
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        bam_path, index = self._build_scenario_with_ts(tmp_path, "Z")
        config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="ts"),
        )
        pr = run_pipeline(bam_path, index, config=config)

        sm = pr.strand_models.exonic_spliced
        assert sm.n_observations > 0, "ts type 'Z' should be handled by Z-type fallback"
        assert sm.strand_specificity > 0.8


# =====================================================================
# Test: Missing or invalid strand tags
# =====================================================================


class TestMissingAndInvalidTags:
    """Verify scanner handles missing or unexpected tag values gracefully."""

    @staticmethod
    def _build_scenario_without_xs(tmp_path):
        """Build scenario, then strip all XS/ts tags."""
        from rigel.sim import Scenario, ReadSimConfig

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
        sc = ReadSimConfig(
            frag_mean=250,
            frag_std=40,
            frag_min=100,
            frag_max=500,
            read_length=100,
            strand_specificity=0.95,
            seed=42,
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
        """No XS/ts tags → zero exonic_spliced observations (scanner level).

        Asserts the scanner behaviour directly via ``scan_and_buffer``; what the full pipeline
        makes of such a library is ``test_a_library_without_spliced_reads_calibrates_as_unstranded``.
        """
        from rigel.config import BamScanConfig
        from rigel.pipeline import scan_and_buffer

        bam_path, index = self._build_scenario_without_xs(tmp_path)
        _, strand_models, _, _ = scan_and_buffer(bam_path, index, BamScanConfig(sj_strand_tag="XS"))

        sm = strand_models.exonic_spliced
        assert sm.n_observations == 0, (
            "Without strand tags, exonic_spliced should have 0 observations"
        )

    def test_a_library_without_spliced_reads_calibrates_as_unstranded(self, tmp_path):
        """Zero spliced strand observations → κ = ½, od = 0, a dead strand channel, and the run completes.

        With no spliced read there is no strand anchor, so the strand channel carries nothing and the
        deconvolution reads density alone; the fallback must work with no RNA (owner, 2026-09-30), so it
        is a valid library, never an error.
        """
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        bam_path, index = self._build_scenario_without_xs(tmp_path)
        config = PipelineConfig(em=EMConfig(seed=42), scan=BamScanConfig(sj_strand_tag="XS"))
        result = run_pipeline(bam_path, index, config=config)
        cal = result.calibration
        assert cal.rna_sense_frac == 0.5
        assert cal.gdna_strand_overdispersion == 0.0
        assert cal.rna_strand_overdispersion == 0.0

    def test_sj_mode_none_ignores_tags(self, tmp_path):
        """sj_strand_tag='none' → ignore all strand tags (scanner level)."""
        from rigel.sim import Scenario, ReadSimConfig
        from rigel.config import BamScanConfig
        from rigel.pipeline import scan_and_buffer

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
        sc = ReadSimConfig(
            frag_mean=250,
            frag_std=40,
            frag_min=100,
            frag_max=500,
            read_length=100,
            strand_specificity=0.95,
            seed=42,
        )
        result = scenario.build_oracle(n_fragments=200, sim_config=sc)
        _, strand_models, _, _ = scan_and_buffer(
            str(result.bam_path), result.index, BamScanConfig(sj_strand_tag="none")
        )

        sm = strand_models.exonic_spliced
        assert sm.n_observations == 0, (
            "With sj_strand_tag='none', exonic_spliced should have 0 "
            "observations even if tags exist"
        )

    def test_invalid_xs_value_treated_as_none(self, tmp_path):
        """XS tag with value '.' or '*' → STRAND_NONE (scanner level, not crash)."""
        from rigel.sim import Scenario, ReadSimConfig
        from rigel.config import BamScanConfig
        from rigel.pipeline import scan_and_buffer

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
        sc = ReadSimConfig(
            frag_mean=250,
            frag_std=40,
            frag_min=100,
            frag_max=500,
            read_length=100,
            strand_specificity=0.95,
            seed=42,
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

        # Invalid values produce STRAND_NONE (no crash); scanner-level check.
        _, strand_models, _, _ = scan_and_buffer(
            sorted_bam, result.index, BamScanConfig(sj_strand_tag="XS")
        )
        sm = strand_models.exonic_spliced
        assert sm.n_observations == 0


# =====================================================================
# Test: Oracle BAM XS tag output type
# =====================================================================


class TestOracleBamXSType:
    """Verify that the oracle BAM writes XS as type 'A' (single char)
    as required by the BAM specification and htslib ``bam_aux2A``."""

    def test_oracle_bam_writes_xs_type_a(self, tmp_path):
        """Oracle BAM should write XS as type 'A', not 'Z'."""
        from rigel.sim import Scenario, ReadSimConfig

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
        sc = ReadSimConfig(
            frag_mean=250,
            frag_std=40,
            frag_min=100,
            frag_max=500,
            read_length=100,
            strand_specificity=0.95,
            seed=42,
        )
        result = scenario.build_oracle(n_fragments=50, sim_config=sc)

        with pysam.AlignmentFile(result.bam_path, "rb") as bam:
            for read in bam:
                if read.has_tag("XS"):
                    # get_tags with_value_type returns (tag, value, type)
                    tag_dict = {t: (v, c) for t, v, c in read.get_tags(with_value_type=True)}
                    xs_val, xs_type = tag_dict["XS"]
                    assert xs_type == "A", (
                        f"XS tag should be type 'A' (char), got type '{xs_type}' value '{xs_val}'"
                    )
                    assert xs_val in ("+", "-"), f"XS value should be '+' or '-', got '{xs_val}'"


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
        from rigel.sim import Scenario, ReadSimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

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
        sc = ReadSimConfig(
            frag_mean=250,
            frag_std=40,
            frag_min=100,
            frag_max=500,
            read_length=100,
            strand_specificity=0.95,
            seed=42,
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
        """NH=0 (rare, boundary case) should not crash."""
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
        from rigel.sim import Scenario, ReadSimConfig

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
        sc = ReadSimConfig(
            frag_mean=250,
            frag_std=40,
            frag_min=100,
            frag_max=500,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
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
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

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
        from rigel.sim import Scenario, ReadSimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

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
        sc = ReadSimConfig(
            frag_mean=250,
            frag_std=40,
            frag_min=100,
            frag_max=500,
            read_length=100,
            strand_specificity=0.95,
            seed=42,
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
        assert sm.n_observations > 0, "Fallback to ts should produce exonic_spliced observations"
        assert sm.strand_specificity > 0.8


# =====================================================================
# Test: the tag read per record
# =====================================================================


@pytest.fixture(scope="module")
def plus_strand_library(tmp_path_factory):
    """A stranded library of one '+' gene: every true junction motif is '+'."""
    from rigel.sim import ReadSimConfig, Scenario

    work = tmp_path_factory.mktemp("per_record_tags")
    scenario = Scenario("per_record_tags", genome_length=100_000, seed=42, work_dir=work)
    scenario.add_gene(
        gene_id="G1",
        strand="+",
        transcripts=[
            {"t_id": "T1", "exons": [(1000, 2000), (3000, 4000), (5000, 6000)], "abundance": 100.0}
        ],
    )
    sim = ReadSimConfig(
        frag_mean=250,
        frag_std=40,
        frag_min=100,
        frag_max=500,
        read_length=100,
        strand_specificity=0.95,
        seed=42,
    )
    result = scenario.build_oracle(n_fragments=600, sim_config=sim)
    yield work, result
    scenario.cleanup()


def _retag(src, dst, *, xs, ts):
    """``src`` rewritten in its own record order, every spliced record given XS ``xs`` (None: no XS) and a ts naming
    GENOMIC strand ``ts`` (None: no ts), written read-relative as minimap2 writes it. Returns the spliced records."""
    flip = {"+": "-", "-": "+"}
    n = 0
    with (
        pysam.AlignmentFile(str(src), "rb") as inp,
        pysam.AlignmentFile(str(dst), "wb", template=inp) as out,
    ):
        for r in inp:
            if any(op == 3 for op, _ in (r.cigartuples or ())):
                n += 1
                r.set_tag("XS", None)
                r.set_tag("ts", None)
                if xs is not None:
                    r.set_tag("XS", xs, "A")
                if ts is not None:
                    r.set_tag("ts", flip[ts] if r.is_reverse else ts, "A")
            out.write(r)
    return n


class TestTagReadPerRecord:
    """``auto`` reads XS and ts on every spliced record. Every spliced record is counted exactly once: given a
    strand, given none, or, when its two tags disagree, ambiguous."""

    @staticmethod
    def _scan(bam, index, spec):
        from rigel.config import BamScanConfig
        from rigel.pipeline import scan_and_buffer

        stats, strand_models, _, _ = scan_and_buffer(
            str(bam), index, BamScanConfig(sj_strand_tag=spec)
        )
        counts = (stats.n_sj_tag_tagged, stats.n_sj_tag_untagged, stats.n_sj_tag_conflicting)
        return counts, stats.n_sj_tag_invalid, strand_models.exonic_spliced

    def test_agreeing_tags_give_the_strand(self, plus_strand_library):
        work, result = plus_strand_library
        n = _retag(result.bam_path, work / "agree.bam", xs="+", ts="+")
        counts, invalid, sm = self._scan(work / "agree.bam", result.index, "auto")
        assert counts == (n, 0, 0) and invalid == 0
        assert sm.n_observations > 0 and sm.strand_specificity > 0.8

    def test_disagreeing_tags_are_ambiguous_counted_and_train_no_strand(self, plus_strand_library):
        work, result = plus_strand_library
        n = _retag(result.bam_path, work / "conflict.bam", xs="+", ts="-")
        counts, invalid, sm = self._scan(work / "conflict.bam", result.index, "auto")
        assert counts == (0, 0, n) and invalid == 0
        assert sm.n_observations == 0

    def test_an_ordered_spec_takes_the_first_tag_present(self, plus_strand_library):
        """A named order is the user's choice: XS, the truth here, is read, and no disagreement is looked for."""
        work, result = plus_strand_library
        n = _retag(result.bam_path, work / "conflict_ordered.bam", xs="+", ts="-")
        counts, _, sm = self._scan(work / "conflict_ordered.bam", result.index, "XS,ts")
        assert counts == (n, 0, 0)
        assert sm.strand_specificity > 0.8

    def test_an_unreadable_xs_does_not_hide_a_readable_ts(self, plus_strand_library):
        work, result = plus_strand_library
        n = _retag(result.bam_path, work / "invalid_xs.bam", xs=".", ts="+")
        counts, invalid, sm = self._scan(work / "invalid_xs.bam", result.index, "auto")
        assert counts == (n, 0, 0) and invalid == n
        assert sm.n_observations > 0 and sm.strand_specificity > 0.8
        counts, invalid, _ = self._scan(work / "invalid_xs.bam", result.index, "XS")
        assert counts == (0, n, 0) and invalid == n
