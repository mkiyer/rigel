"""
Tests for the annotated BAM output feature (rigel.annotate).

Verifies that:
1. The AnnotationTable collects correct per-fragment metadata.
2. write_annotated_bam produces a valid BAM with the expected tags.
3. The annotation pipeline integrates end-to-end with run_pipeline.
"""

import shutil

import numpy as np
import pytest

from rigel.annotate import (
    AnnotationTable,
    AF_GDNA,
    AF_GDNA_RESOLVED,
    AF_NRNA,
    AF_NRNA_RESOLVED,
    AF_RESOLVED,
    AF_SYNTH_RESOLVED,
    AF_SYNTHETIC,
    AF_TRANSCRIPT,
    AF_UNRESOLVED,
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
            tx_flags=AF_TRANSCRIPT,
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
        assert ann["tx_flags"] == AF_TRANSCRIPT
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
            tbl.add(frag_id=i, best_tid=i, best_gid=0, tx_flags=AF_TRANSCRIPT,
                    posterior=1.0, frag_class=0, n_candidates=1)
        assert tbl.size == 10
        assert tbl.capacity >= 10
        for i in range(10):
            ann = tbl.get(i)
            assert ann is not None
            assert ann["best_tid"] == i

    def test_zf_flags_valid_values(self):
        """Valid ZF values follow the implication chain."""
        assert AF_UNRESOLVED == 0
        assert AF_TRANSCRIPT == AF_RESOLVED  # 1
        assert AF_GDNA_RESOLVED == (AF_RESOLVED | AF_GDNA)  # 3
        assert AF_NRNA_RESOLVED == (AF_RESOLVED | AF_NRNA)  # 5
        assert AF_SYNTH_RESOLVED == (AF_RESOLVED | AF_NRNA | AF_SYNTHETIC)  # 13
        # All resolved values are odd (bit 0 set)
        for v in (AF_TRANSCRIPT, AF_GDNA_RESOLVED, AF_NRNA_RESOLVED, AF_SYNTH_RESOLVED):
            assert v & AF_RESOLVED, f"ZF={v} should have is_resolved set"

    def test_frag_class_labels(self):
        """Known fragment class codes have labels."""
        assert _FRAG_CLASS_LABELS[0] == "unambig"
        assert _FRAG_CLASS_LABELS[1] == "ambig_same_strand"
        assert _FRAG_CLASS_LABELS[3] == "multimapper"
        assert _FRAG_CLASS_LABELS[-1] == "intergenic"

    def test_splice_type_label(self):
        """SpliceType codes convert to lowercase labels."""
        from rigel.splice import SpliceType
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
        from rigel.sim import Scenario
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
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.sim import SimConfig
        from rigel.pipeline import run_pipeline

        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        result = scenario.build(
            n_fragments=200, sim_config=sim_config,
        )

        annotated_bam = tmp_path / "annotated.bam"
        run_pipeline(
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
        expected_tags = {"ZT", "ZG", "ZR", "ZF", "ZW", "ZC", "ZH", "ZN", "ZS"}
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
        assert isinstance(rec.get_tag("ZR"), str)
        assert isinstance(rec.get_tag("ZF"), int)
        assert isinstance(rec.get_tag("ZW"), float)
        assert isinstance(rec.get_tag("ZC"), str)
        assert isinstance(rec.get_tag("ZH"), int)
        assert isinstance(rec.get_tag("ZN"), int)
        assert isinstance(rec.get_tag("ZS"), str)

        # ZF values should be valid (0, 1, 3, 5, or 13)
        valid_zf = {AF_UNRESOLVED, AF_TRANSCRIPT, AF_GDNA_RESOLVED,
                     AF_NRNA_RESOLVED, AF_SYNTH_RESOLVED}
        for rec in records:
            zf = rec.get_tag("ZF")
            assert zf in valid_zf, f"Invalid ZF: {zf}"

        # ZH should be 0 or 1
        for rec in records:
            assert rec.get_tag("ZH") in (0, 1)

        # ZW should be in [0, 1]
        for rec in records:
            w = rec.get_tag("ZW")
            assert 0.0 <= w <= 1.0 + 1e-6, f"ZW out of range: {w}"

        # At least some records should have transcript assignment (ZF & 1 set)
        transcript_records = [
            r for r in records if (r.get_tag("ZF") & AF_RESOLVED)
        ]
        assert len(transcript_records) > 0, "No resolved assignments found"

    def test_annotated_bam_counts_match(self, scenario, tmp_path):
        """Pipeline counts are identical with and without annotation."""
        from rigel.sim import SimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        result = scenario.build(
            n_fragments=200, sim_config=sim_config,
        )

        base_config = PipelineConfig(
            em=EMConfig(seed=42, assignment_mode="fractional"),
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
            counts1["count"].values,
            counts2["count"].values,
            decimal=6,
            err_msg="Annotation mode changed count totals",
        )

    def test_annotated_bam_preserves_records_and_collation(
        self, scenario, tmp_path
    ):
        """Annotated BAM MUST contain exactly the same records as input
        and MUST remain collated (qnames grouped contiguously).

        This is a contract for downstream tools: rigel accepts collated
        BAM in and produces collated BAM out with no record drops or
        duplications.
        """
        import pysam
        from collections import Counter
        from rigel.sim import SimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        result = scenario.build(n_fragments=500, sim_config=sim_config)

        annotated_bam = tmp_path / "annotated.bam"
        run_pipeline(
            result.bam_path,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
                annotated_bam_path=annotated_bam,
            ),
        )

        def _record_key(r):
            # Identity of a BAM record = qname + flag + ref + pos.
            # Allows counting duplicates even if tags differ.
            return (
                r.query_name,
                int(r.flag),
                int(r.reference_id),
                int(r.reference_start),
            )

        # ---- 1. Record-count parity ----
        with pysam.AlignmentFile(str(result.bam_path), "rb") as fin:
            in_records = [_record_key(r) for r in fin.fetch(until_eof=True)]
        with pysam.AlignmentFile(str(annotated_bam), "rb") as fout:
            out_records = [_record_key(r) for r in fout.fetch(until_eof=True)]

        assert len(out_records) == len(in_records), (
            f"Record count changed: in={len(in_records)} out={len(out_records)}"
        )

        # Same multiset (no drops, no duplications)
        assert Counter(in_records) == Counter(out_records), (
            "Annotated BAM does not contain exactly the same records as input"
        )

        # ---- 2. Collation preservation ----
        # A collated BAM has each qname appearing in a single contiguous
        # run. Equivalent to: number of distinct qname runs equals
        # number of distinct qnames.
        def _qname_runs(records):
            runs = 0
            prev = object()
            for r in records:
                qn = r[0]
                if qn != prev:
                    runs += 1
                    prev = qn
            return runs

        n_unique_qnames = len({r[0] for r in out_records})
        n_runs = _qname_runs(out_records)
        assert n_runs == n_unique_qnames, (
            f"Annotated BAM is not collated: "
            f"{n_runs} qname runs vs {n_unique_qnames} unique qnames"
        )

    def test_annotated_bam_preserves_records_with_multimappers(
        self, scenario, tmp_path
    ):
        """Record-identity preservation under NH>1 multimappers.

        Simulates a standard scenario, then fabricates NH=2 cross-mapped
        secondary alignments for a subset of read pairs to force the
        writer through the no-HI multimapper pairing path
        (``pair_multimapper_reads``).  The output BAM must still contain
        exactly the same multiset of records as the augmented input.
        """
        import pysam
        from collections import Counter
        from rigel.sim import SimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        result = scenario.build(n_fragments=200, sim_config=sim_config)

        # Augment the BAM: add NH=2 secondary alignments for a subset
        # of primary records, WITHOUT HI tags, at a different ref
        # position.  This exercises the non-HI multimapper path and
        # the cross-pair combinatorial block in pair_multimapper_reads.
        augmented = tmp_path / "augmented.bam"
        with pysam.AlignmentFile(str(result.bam_path), "rb") as fin:
            header = fin.header
            records = list(fin.fetch(until_eof=True))

        # Group by qname (input is collated) and set NH=2 on every
        # record of a subset, then append a secondary copy per record.
        from collections import defaultdict
        by_qname = defaultdict(list)
        for r in records:
            by_qname[r.query_name].append(r)
        qnames = list(by_qname.keys())
        # Pick 10 fragments to turn into multimappers
        mm_qnames = set(qnames[:10])

        augmented_records = []
        for qn, recs in by_qname.items():
            if qn in mm_qnames and len(recs) >= 2:
                # Stamp NH=2 on primaries and add a secondary copy
                # at an offset genomic position on the same ref.
                new_group = []
                for r in recs:
                    r.set_tag("NH", 2)
                    new_group.append(r)
                for r in recs:
                    sec = r.__copy__()
                    sec.flag = r.flag | 0x100  # mark secondary
                    sec.reference_start = max(
                        0, r.reference_start + 500
                    )
                    sec.set_tag("NH", 2)
                    new_group.append(sec)
                augmented_records.extend(new_group)
            else:
                augmented_records.extend(recs)

        with pysam.AlignmentFile(
            str(augmented), "wb", header=header
        ) as fout:
            for r in augmented_records:
                fout.write(r)

        annotated_bam = tmp_path / "annotated_mm.bam"
        run_pipeline(
            augmented,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(
                    sj_strand_tag="ts", include_multimap=True,
                ),
                annotated_bam_path=annotated_bam,
            ),
        )

        def _record_key(r):
            return (
                r.query_name,
                int(r.flag),
                int(r.reference_id),
                int(r.reference_start),
            )

        with pysam.AlignmentFile(str(augmented), "rb") as fin:
            in_keys = [_record_key(r) for r in fin.fetch(until_eof=True)]
        with pysam.AlignmentFile(str(annotated_bam), "rb") as fout:
            out_keys = [_record_key(r) for r in fout.fetch(until_eof=True)]

        assert len(out_keys) == len(in_keys), (
            f"Multimapper record count changed: "
            f"in={len(in_keys)} out={len(out_keys)}"
        )
        assert Counter(in_keys) == Counter(out_keys), (
            "Multimapper annotated BAM does not preserve input records "
            "exactly (drops or duplicates detected)."
        )

        # Collation check
        def _qname_runs(keys):
            runs, prev = 0, object()
            for k in keys:
                if k[0] != prev:
                    runs += 1
                    prev = k[0]
            return runs

        assert _qname_runs(out_keys) == len({k[0] for k in out_keys}), (
            "Multimapper annotated BAM is not collated"
        )

    # =================================================================
    # ZB:i — per-record splice-artifact tag
    # =================================================================

    def test_zb_tag_schema_always_present(self, scenario, tmp_path):
        """Every annotated record carries ZB:i, and its value is >= 0.

        This scenario has no blacklist loaded, so every ZB should be 0.
        """
        import pysam
        from rigel.sim import SimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        result = scenario.build(n_fragments=200, sim_config=sim_config)

        annotated_bam = tmp_path / "annotated_zb.bam"
        run_pipeline(
            result.bam_path,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
                annotated_bam_path=annotated_bam,
            ),
        )

        with pysam.AlignmentFile(str(annotated_bam), "rb") as fout:
            records = list(fout.fetch(until_eof=True))
        assert len(records) > 0

        # Every record that receives any annotation tag must also have ZB,
        # and its value must be a non-negative int.
        for rec in records:
            tags = {t[0] for t in rec.get_tags()}
            if "ZT" in tags:  # record was annotated by rigel
                assert "ZB" in tags, (
                    f"Record {rec.query_name} missing ZB tag"
                )
                zb = rec.get_tag("ZB")
                assert isinstance(zb, int)
                assert zb >= 0
                # No blacklist loaded in this scenario → ZB must be 0
                assert zb == 0, (
                    f"ZB={zb} with no blacklist loaded (expected 0)"
                )

    def test_zb_tag_marks_blacklisted_junctions(self, scenario, tmp_path):
        """Injecting a targeted blacklist entry produces ZB>=1 on the
        matching reads, and the sum of ZB across the output BAM equals
        the Pass-1 scanner stat ``n_sj_blacklisted`` (Pass-1/Pass-2
        parity invariant).
        """
        import pandas as pd
        import pysam
        from rigel.sim import SimConfig
        from rigel.index import TranscriptIndex, SJ_BLACKLIST_FEATHER
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        result = scenario.build(n_fragments=500, sim_config=sim_config)

        # Inject a blacklist covering gene g1's junction (1000, 1500)
        # with huge anchor envelopes so every spliced read over it is
        # dropped.  The scenario's ref name equals scenario.name by default.
        ref = scenario.ref_name
        bl_df = pd.DataFrame({
            "ref": [ref],
            "start": pd.array([1000], dtype="int32"),
            "end": pd.array([1500], dtype="int32"),
            "max_anchor_left": pd.array([1000], dtype="int32"),
            "max_anchor_right": pd.array([1000], dtype="int32"),
        })
        bl_df.to_feather(result.index_dir / SJ_BLACKLIST_FEATHER)

        # Reload so the C++ resolver picks up the new blacklist
        reloaded_index = TranscriptIndex.load(result.index_dir)

        annotated_bam = tmp_path / "annotated_zb_bl.bam"
        pr = run_pipeline(
            result.bam_path,
            reloaded_index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
                annotated_bam_path=annotated_bam,
            ),
        )

        with pysam.AlignmentFile(str(annotated_bam), "rb") as fout:
            records = list(fout.fetch(until_eof=True))

        zb_values = [
            rec.get_tag("ZB") for rec in records if rec.has_tag("ZB")
        ]
        n_marked = sum(1 for v in zb_values if v > 0)
        zb_sum = sum(zb_values)

        assert n_marked > 0, (
            "Expected at least one record with ZB>0 after injecting a "
            "blacklist, got none."
        )

        # Pass-1 / Pass-2 parity: the sum of ZB over the annotated BAM
        # equals the scanner's n_sj_blacklisted stat.
        assert pr.stats.n_sj_blacklisted > 0, (
            "Pass-1 stat n_sj_blacklisted should be > 0"
        )
        assert zb_sum == pr.stats.n_sj_blacklisted, (
            f"ZB sum ({zb_sum}) != Pass-1 n_sj_blacklisted "
            f"({pr.stats.n_sj_blacklisted})"
        )

    def test_zb_tag_absent_on_filtered_passthrough(self, scenario, tmp_path):
        """Filtered pass-through records (QCFAIL etc.) receive no Z*
        tags, including no ZB.  Augment the input BAM with a QCFAIL
        record and verify the output preserves it verbatim without tags.
        """
        import pysam
        from rigel.sim import SimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        result = scenario.build(n_fragments=100, sim_config=sim_config)

        # Clone a primary record and flip the QCFAIL bit.
        augmented = tmp_path / "aug_qcfail.bam"
        with pysam.AlignmentFile(str(result.bam_path), "rb") as fin:
            header = fin.header
            records = list(fin.fetch(until_eof=True))
        target = records[0]
        qc = target.__copy__()
        qc.query_name = "QCFAIL_SENTINEL"
        qc.flag = qc.flag | 0x200  # BAM_FQCFAIL
        # Prepend so it is clearly outside any existing qname group
        records.insert(0, qc)
        with pysam.AlignmentFile(
            str(augmented), "wb", header=header
        ) as fout:
            for r in records:
                fout.write(r)

        annotated_bam = tmp_path / "annotated_qcfail.bam"
        run_pipeline(
            augmented,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
                annotated_bam_path=annotated_bam,
            ),
        )

        with pysam.AlignmentFile(str(annotated_bam), "rb") as fout:
            qc_out = [
                r for r in fout.fetch(until_eof=True)
                if r.query_name == "QCFAIL_SENTINEL"
            ]
        assert len(qc_out) == 1, (
            f"QCFAIL sentinel not preserved: found {len(qc_out)} copies"
        )
        tags = {t[0] for t in qc_out[0].get_tags()}
        # Filtered pass-through records get no rigel annotations
        for tag in ("ZT", "ZG", "ZF", "ZB"):
            assert tag not in tags, (
                f"Filtered pass-through record carries {tag}; should not"
            )

