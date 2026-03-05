"""
Tests for the annotated BAM output feature (hulkrna.annotate).

Verifies that:
1. The AnnotationTable collects correct per-fragment metadata.
2. write_annotated_bam produces a valid BAM with the expected tags.
3. The annotation pipeline integrates end-to-end with run_pipeline.
"""

import shutil

import numpy as np
import pytest

from hulkrna.annotate import (
    AnnotationTable,
    POOL_CODE,
    POOL_LABEL,
    POOL_MRNA,
    POOL_NRNA,
    POOL_GDNA,
    POOL_INTERGENIC,
    POOL_CHIMERIC,
    _FRAG_CLASS_LABELS,
    _splice_type_label,
)


# =====================================================================
# Unit tests: AnnotationTable
# =====================================================================


class TestAnnotationTable:
    """Tests for the AnnotationTable data structure."""

    def test_create_empty(self):
        tbl = AnnotationTable.create(100)
        assert tbl.size == 0
        assert tbl.capacity == 100

    def test_add_and_get(self):
        tbl = AnnotationTable.create(10)
        tbl.add(
            frag_id=42,
            best_tid=5,
            best_gid=2,
            pool=POOL_CODE[POOL_MRNA],
            posterior=0.95,
            frag_class=0,
            n_candidates=3,
            splice_type=1,
        )
        assert tbl.size == 1
        ann = tbl.get(42)
        assert ann is not None
        assert ann["best_tid"] == 5
        assert ann["best_gid"] == 2
        assert ann["pool"] == POOL_CODE[POOL_MRNA]
        assert abs(ann["posterior"] - 0.95) < 0.01
        assert ann["frag_class"] == 0
        assert ann["n_candidates"] == 3
        assert ann["splice_type"] == 1

    def test_get_missing(self):
        tbl = AnnotationTable.create(10)
        assert tbl.get(999) is None

    def test_grow(self):
        tbl = AnnotationTable.create(2)
        for i in range(10):
            tbl.add(frag_id=i, best_tid=i, best_gid=0, pool=0,
                    posterior=1.0, frag_class=0, n_candidates=1)
        assert tbl.size == 10
        assert tbl.capacity >= 10
        for i in range(10):
            ann = tbl.get(i)
            assert ann is not None
            assert ann["best_tid"] == i

    def test_pool_codes_roundtrip(self):
        """All pool labels have codes and vice versa."""
        for label, code in POOL_CODE.items():
            assert POOL_LABEL[code] == label

    def test_frag_class_labels(self):
        """Known fragment class codes have labels."""
        assert _FRAG_CLASS_LABELS[0] == "unambig"
        assert _FRAG_CLASS_LABELS[1] == "ambig_same_strand"
        assert _FRAG_CLASS_LABELS[3] == "multimapper"
        assert _FRAG_CLASS_LABELS[-1] == "intergenic"

    def test_splice_type_label(self):
        """SpliceType codes convert to lowercase labels."""
        from hulkrna.splice import SpliceType
        for st in SpliceType:
            label = _splice_type_label(int(st))
            assert label == st.name.lower()


# =====================================================================
# Integration test: annotated BAM via run_pipeline
# =====================================================================

pytestmark_scenario = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and/or samtools not found in PATH",
)


@pytestmark_scenario
class TestAnnotatedBamIntegration:
    """End-to-end: scenario → run_pipeline with annotated_bam_path."""

    @pytest.fixture
    def scenario(self, tmp_path):
        from hulkrna.sim import Scenario
        sc = Scenario(
            "annot_bam_test", genome_length=8000, seed=42,
            work_dir=tmp_path / "annot_bam",
        )
        # Spliced gene for strand training + deterministic unambig
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 1000), (1500, 2000)],
             "abundance": 80},
        ])
        # Second gene for EM-routed isoform ambiguity
        sc.add_gene("g2", "+", [
            {"t_id": "t2",
             "exons": [(3000, 3500), (4000, 4500)],
             "abundance": 40},
            {"t_id": "t3",
             "exons": [(3000, 3500), (4000, 4500), (5000, 5500)],
             "abundance": 40},
        ])
        yield sc
        sc.cleanup()

    def test_annotated_bam_produced(self, scenario, tmp_path):
        """run_pipeline with annotated_bam_path writes a valid BAM."""
        import pysam
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.sim import SimConfig, run_benchmark
        from hulkrna.pipeline import run_pipeline

        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        result = scenario.build(
            n_fragments=200, sim_config=sim_config,
        )

        annotated_bam = tmp_path / "annotated.bam"
        pr = run_pipeline(
            result.bam_path,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
                annotated_bam_path=annotated_bam,
            ),
        )

        # Verify the annotated BAM exists and has records
        assert annotated_bam.exists(), "Annotated BAM not created"
        bam = pysam.AlignmentFile(str(annotated_bam), "rb")
        records = list(bam.fetch(until_eof=True))
        bam.close()
        assert len(records) > 0, "Annotated BAM is empty"

        # Every record should have the expected tags
        expected_tags = {"ZT", "ZG", "ZP", "ZW", "ZC", "ZH", "ZN", "ZS"}
        for rec in records:
            tags_present = {t[0] for t in rec.get_tags()}
            missing = expected_tags - tags_present
            assert not missing, (
                f"Record {rec.query_name} missing tags: {missing}"
            )

        # Check tag value types
        rec = records[0]
        assert isinstance(rec.get_tag("ZT"), str)
        assert isinstance(rec.get_tag("ZG"), str)
        assert isinstance(rec.get_tag("ZP"), str)
        assert isinstance(rec.get_tag("ZW"), float)
        assert isinstance(rec.get_tag("ZC"), str)
        assert isinstance(rec.get_tag("ZH"), int)
        assert isinstance(rec.get_tag("ZN"), int)
        assert isinstance(rec.get_tag("ZS"), str)

        # Pool values should be valid
        valid_pools = {POOL_MRNA, POOL_NRNA, POOL_GDNA,
                       POOL_INTERGENIC, POOL_CHIMERIC}
        for rec in records:
            assert rec.get_tag("ZP") in valid_pools, (
                f"Invalid pool: {rec.get_tag('ZP')}"
            )

        # ZH should be 0 or 1
        for rec in records:
            assert rec.get_tag("ZH") in (0, 1)

        # ZW should be in [0, 1]
        for rec in records:
            w = rec.get_tag("ZW")
            assert 0.0 <= w <= 1.0 + 1e-6, f"ZW out of range: {w}"

        # At least some records should have mRNA assignment
        mrna_records = [r for r in records if r.get_tag("ZP") == POOL_MRNA]
        assert len(mrna_records) > 0, "No mRNA assignments found"

    def test_annotated_bam_counts_match(self, scenario, tmp_path):
        """Pipeline counts are identical with and without annotation."""
        from hulkrna.sim import SimConfig
        from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
        from hulkrna.pipeline import run_pipeline

        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        result = scenario.build(
            n_fragments=200, sim_config=sim_config,
        )

        base_config = PipelineConfig(
            em=EMConfig(seed=42),
            scan=BamScanConfig(sj_strand_tag="ts"),
        )

        # Run without annotation
        pr1 = run_pipeline(result.bam_path, result.index, config=base_config)
        counts1 = pr1.estimator.get_counts_df(result.index)

        # Run with annotation
        annotated_bam = tmp_path / "annotated.bam"
        from dataclasses import replace as _replace
        annot_config = _replace(base_config, annotated_bam_path=annotated_bam)
        pr2 = run_pipeline(
            result.bam_path, result.index, config=annot_config,
        )
        counts2 = pr2.estimator.get_counts_df(result.index)

        # Counts should be identical (annotation doesn't change quantification)
        np.testing.assert_array_almost_equal(
            counts1["mrna"].values,
            counts2["mrna"].values,
            decimal=6,
            err_msg="Annotation mode changed mrna totals",
        )
