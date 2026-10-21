"""`rigel.annotate` — writing the scanner's per-fragment decisions back out as BAM tags.

The annotation table itself: what it stores per fragment, how it grows, the canonical ZF values,
and the category a transcript's win is stamped with. Then the writer, through a real
`run_pipeline`: the annotated BAM must be valid, its tag values must lie in their domains, a
transcript's win must read annotated or nascent as the output tables do, its per-class counts must
agree with the run's own, and it must preserve every input record and its collation, multimappers
included, because an annotated BAM that silently drops records is not the same library the numbers
describe. The blacklist tag is gated for presence, for the blacklisted case and for its absence on
a filtered passthrough.
"""

import shutil
from collections import defaultdict

import numpy as np
import pytest

from rigel.annotate import (
    AnnotationTable,
    AF_ANNOTATED,
    AF_ANNOTATED_BIT,
    AF_GDNA_BIT,
    AF_GDNA_EM,
    AF_GDNA_INTERGENIC,
    AF_INTERGENIC_BIT,
    AF_MULTIMAPPER_DROP,
    AF_MULTIMAPPER_DROP_BIT,
    AF_NASCENT,
    AF_NASCENT_BIT,
    AF_RESOLVED,
    AF_UNRESOLVED,
    winner_flag,
    winner_flags,
)


# Set of the only legitimate ZF outputs produced by rigel.
VALID_ZF_VALUES = {
    AF_UNRESOLVED,
    AF_ANNOTATED,
    AF_NASCENT,
    AF_GDNA_EM,
    AF_GDNA_INTERGENIC,
    AF_MULTIMAPPER_DROP,
}

#: Bits 3 and 6, retired: bit 3 marked a single-exon winner (0x09) or a synthetic span's (0x19) under an
#: earlier layout, and bit 6 the resolver's chimera call (0x40), deleted; a value carrying either can only
#: be read with that layout's meaning.
RETIRED_BITS = 0x08 | 0x40

#: Every bit the table names; a ZF value carries no other.
NAMED_BITS = (
    AF_RESOLVED
    | AF_ANNOTATED_BIT
    | AF_GDNA_BIT
    | AF_NASCENT_BIT
    | AF_INTERGENIC_BIT
    | AF_MULTIMAPPER_DROP_BIT
)


def _assert_zf_invariants(zf: int) -> None:
    """Assert the four ZF invariants `rigel.annotate` states, and that no bit outside its table is set."""
    is_resolved = bool(zf & AF_RESOLVED)
    is_annotated = bool(zf & AF_ANNOTATED_BIT)
    is_gdna = bool(zf & AF_GDNA_BIT)
    is_nascent = bool(zf & AF_NASCENT_BIT)
    is_interg = bool(zf & AF_INTERGENIC_BIT)
    is_mm_drop = bool(zf & AF_MULTIMAPPER_DROP_BIT)
    n_kind = int(is_annotated) + int(is_gdna) + int(is_nascent)

    # 1: resolved XOR mm_dropped — every stamped record is exactly one or the other (the unresolved
    #    sentinel ZF=0 also satisfies this: neither side is set).
    assert not (is_resolved and is_mm_drop), f"ZF={zf:#x}: resolved AND dropped bits both set"
    # 2: a winner-kind bit only when resolved, and then exactly one, or none for RNA no annotated
    #    transcript holds (spliced intergenic)
    if is_resolved:
        assert n_kind == 1 or (n_kind == 0 and is_interg), (
            f"ZF={zf:#x}: resolved but {n_kind} winner-kind bits set"
        )
    else:
        assert n_kind == 0, f"ZF={zf:#x}: winner-kind bit set without is_resolved"
    # 3: intergenic -> resolved, gDNA when unspliced and RNA (no winner-kind bit) when spliced
    if is_interg:
        assert is_resolved and not (is_annotated or is_nascent), (
            f"ZF={zf:#x}: is_intergenic misplaced"
        )
    # 4: the retired bits are never set, nor any bit the table does not name
    assert not zf & RETIRED_BITS, f"ZF={zf:#x}: a retired bit (0x08, 0x40) is set"
    assert not zf & ~NAMED_BITS, f"ZF={zf:#x}: a bit outside the table is set"


# =====================================================================
# Unit tests: AnnotationTable
# =====================================================================


class TestAnnotationTable:
    """Tests for the AnnotationTable data structure."""

    def test_create_empty(self):
        tbl = AnnotationTable.create(100)
        assert tbl.size == 0
        assert tbl.capacity == 100

    def test_add(self):
        tbl = AnnotationTable.create(10)
        tbl.add(
            frag_id=42,
            best_tid=5,
            best_gid=2,
            tx_flags=AF_ANNOTATED,
            posterior=0.95,
            frag_class=0,
            n_candidates=3,
            splice_type=1,
        )
        assert tbl.size == 1
        assert tbl.frag_ids[0] == 42
        assert tbl.best_tid[0] == 5
        assert tbl.best_gid[0] == 2
        assert tbl.tx_flags[0] == AF_ANNOTATED
        assert abs(tbl.posterior[0] - 0.95) < 0.01
        assert tbl.frag_class[0] == 0
        assert tbl.n_candidates[0] == 3
        assert tbl.splice_type[0] == 1

    def test_grow(self):
        tbl = AnnotationTable.create(2)
        for i in range(10):
            tbl.add(
                frag_id=i,
                best_tid=i,
                best_gid=0,
                tx_flags=AF_ANNOTATED,
                posterior=1.0,
                frag_class=0,
                n_candidates=1,
            )
        assert tbl.size == 10
        assert tbl.capacity >= 10
        np.testing.assert_array_equal(tbl.frag_ids[: tbl.size], np.arange(10))
        np.testing.assert_array_equal(tbl.best_tid[: tbl.size], np.arange(10))

    def test_zf_canonical_values(self):
        """All canonical ZF values satisfy the invariants, and none reuses a retired bit (0x08, 0x40)."""
        # Documented composed values
        assert AF_UNRESOLVED == 0x00
        assert AF_ANNOTATED == 0x03
        assert AF_NASCENT == 0x11
        assert AF_GDNA_EM == 0x05
        assert AF_GDNA_INTERGENIC == 0x25
        assert AF_MULTIMAPPER_DROP == 0x80
        # All canonical values satisfy the invariants, which include the retired bits' absence
        for v in VALID_ZF_VALUES:
            _assert_zf_invariants(v)
        # Resolved-ness is coherent
        assert AF_ANNOTATED & AF_RESOLVED
        assert AF_NASCENT & AF_RESOLVED
        assert AF_GDNA_EM & AF_RESOLVED
        assert AF_GDNA_INTERGENIC & AF_RESOLVED
        assert not (AF_MULTIMAPPER_DROP & AF_RESOLVED)

    def test_a_transcript_win_is_annotated_or_nascent(self):
        """A transcript's win reads as the output tables report its RNA: a synthetic span's is nascent
        (0x11), every reference transcript's annotated (0x03), single-exon included, so the helpers take
        only ``is_synthetic``. The scalar and the vectorized form agree."""
        assert winner_flag(False) == AF_ANNOTATED == 0x03
        assert winner_flag(True) == AF_NASCENT == 0x11
        is_synthetic = np.array([False, True, False, True])
        flags = winner_flags(is_synthetic)
        assert flags.dtype == np.uint8
        np.testing.assert_array_equal(flags, [0x03, 0x11, 0x03, 0x11])
        assert [winner_flag(bool(s)) for s in is_synthetic] == flags.tolist()
        for v in flags.tolist():
            _assert_zf_invariants(v)
        assert winner_flags(np.zeros(0, dtype=bool)).dtype == np.uint8


# =====================================================================
# Integration test: annotated BAM via run_pipeline
# =====================================================================

pytestmark_scenario = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and/or samtools not found in PATH",
)


@pytest.fixture
def scenario(tmp_path):
    from rigel.sim import Scenario

    sc = Scenario(
        "annot_bam_test",
        genome_length=8000,
        seed=42,
        work_dir=tmp_path / "annot_bam",
    )
    # Spliced gene for strand training + deterministic unambig
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(500, 1000), (1500, 2000)], "abundance": 80},
        ],
    )
    # Second gene for EM-routed isoform ambiguity
    sc.add_gene(
        "g2",
        "+",
        [
            {"t_id": "t2", "exons": [(3000, 3500), (4000, 4500)], "abundance": 40},
            {
                "t_id": "t3",
                "exons": [(3000, 3500), (4000, 4500), (5000, 5500)],
                "abundance": 40,
            },
        ],
    )
    yield sc
    sc.cleanup()


@pytest.fixture
def winner_kinds_scenario(tmp_path):
    """One of each kind of transcript a fragment can be won by: a multi-exon transcript nothing covers, so
    the index adds a synthetic span over it; a single-exon transcript covering its multi-exon sibling's
    span, so there is no synthetic span there; and a single-exon transcript alone."""
    from rigel.sim import Scenario

    sc = Scenario(
        "annot_winner_kinds",
        genome_length=8000,
        seed=42,
        work_dir=tmp_path / "winner_kinds",
    )
    sc.add_gene(
        "g1", "+", [{"t_id": "t_multi", "exons": [(500, 1000), (1500, 2000)], "abundance": 60}]
    )
    sc.add_gene(
        "g2",
        "+",
        [
            {"t_id": "t_sibling", "exons": [(3000, 3400), (3800, 4200)], "abundance": 30},
            {"t_id": "t_cover", "exons": [(3000, 4200)], "abundance": 30},
        ],
    )
    sc.add_gene("g3", "+", [{"t_id": "t_single", "exons": [(5500, 6500)], "abundance": 60}])
    yield sc
    sc.cleanup()


@pytestmark_scenario
class TestAnnotatedBamIntegration:
    """End-to-end: scenario → run_pipeline with annotated_bam_path."""

    def test_annotated_bam_produced(self, scenario, tmp_path):
        """run_pipeline with annotated_bam_path writes a valid BAM."""
        import pysam
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.sim import ReadSimConfig
        from rigel.pipeline import run_pipeline

        sim_config = ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
        )
        result = scenario.build(
            n_fragments=200,
            sim_config=sim_config,
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
            assert not missing, f"Record {rec.query_name} missing tags: {missing}"

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

        # ZF values must be one of the canonical outputs, and every
        # value must satisfy the redesign invariants.
        for rec in records:
            zf = rec.get_tag("ZF")
            assert zf in VALID_ZF_VALUES, f"Invalid ZF: {zf:#x}"
            _assert_zf_invariants(zf)

        # ZC values: pre-EM ambiguity enum + "." for non-EM records.
        valid_zc = {"unambig", "ambig_same_strand", "ambig_opp_strand", "multimapper", "."}
        for rec in records:
            zc = rec.get_tag("ZC")
            assert zc in valid_zc, f"Unexpected ZC value: {zc!r}"

        # ZS values: every SpliceType, lower-cased, or "unknown".
        from rigel.splice import SpliceType

        valid_zs = {st.name.lower() for st in SpliceType} | {"unknown"}
        for rec in records:
            zs = rec.get_tag("ZS")
            assert zs in valid_zs, f"Unexpected ZS value: {zs!r}"

        # ZH should be 0 or 1
        for rec in records:
            assert rec.get_tag("ZH") in (0, 1)

        # ZW should be in [0, 1]
        for rec in records:
            w = rec.get_tag("ZW")
            assert 0.0 <= w <= 1.0 + 1e-6, f"ZW out of range: {w}"

        # At least some records should have transcript assignment (ZF & 1 set)
        transcript_records = [r for r in records if (r.get_tag("ZF") & AF_RESOLVED)]
        assert len(transcript_records) > 0, "No resolved assignments found"

    def test_annotated_bam_counts_match(self, scenario, tmp_path):
        """Pipeline counts are identical with and without annotation."""
        from rigel.sim import ReadSimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
        )
        result = scenario.build(
            n_fragments=200,
            sim_config=sim_config,
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
            result.bam_path,
            result.index,
            config=annot_config,
        )
        counts2 = pr2.estimator.get_counts_df(result.index)

        # Counts should be identical (annotation doesn't change quantification)
        np.testing.assert_array_almost_equal(
            counts1["count"].values,
            counts2["count"].values,
            decimal=6,
            err_msg="Annotation mode changed count totals",
        )

    def test_annotated_bam_preserves_records_and_collation(self, scenario, tmp_path):
        """Annotated BAM MUST contain exactly the same records as input
        and MUST remain collated (qnames grouped contiguously).

        This is a contract for downstream tools: rigel accepts collated
        BAM in and produces collated BAM out with no record drops or
        duplications.
        """
        import pysam
        from collections import Counter
        from rigel.sim import ReadSimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
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
            f"Annotated BAM is not collated: {n_runs} qname runs vs {n_unique_qnames} unique qnames"
        )

    def test_annotated_bam_preserves_records_with_multimappers(self, scenario, tmp_path):
        """Record-identity preservation under NH>1 multimappers.

        Simulates a standard scenario, then fabricates NH=2 cross-mapped
        secondary alignments for a subset of read pairs, so the writer pairs
        every R1 place with every R2 place and a record joins several hits.
        The output BAM must still contain exactly the same multiset of
        records as the augmented input.
        """
        import pysam
        from collections import Counter
        from rigel.sim import ReadSimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
        )
        result = scenario.build(n_fragments=200, sim_config=sim_config)

        # Augment the BAM: add NH=2 secondary alignments for a subset
        # of primary records at a different ref position, so each R1
        # place pairs with each R2 place.
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
                    sec.reference_start = max(0, r.reference_start + 500)
                    sec.set_tag("NH", 2)
                    new_group.append(sec)
                augmented_records.extend(new_group)
            else:
                augmented_records.extend(recs)

        with pysam.AlignmentFile(str(augmented), "wb", header=header) as fout:
            for r in augmented_records:
                fout.write(r)

        annotated_bam = tmp_path / "annotated_mm.bam"
        run_pipeline(
            augmented,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(
                    sj_strand_tag="ts",
                    include_multimap=True,
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
            f"Multimapper record count changed: in={len(in_keys)} out={len(out_keys)}"
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
        from rigel.sim import ReadSimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
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
                assert "ZB" in tags, f"Record {rec.query_name} missing ZB tag"
                zb = rec.get_tag("ZB")
                assert isinstance(zb, int)
                assert zb >= 0
                # No blacklist loaded in this scenario → ZB must be 0
                assert zb == 0, f"ZB={zb} with no blacklist loaded (expected 0)"

    def test_zb_tag_marks_blacklisted_sj(self, scenario, tmp_path):
        """Injecting a targeted blacklist entry produces ZB>=1 on the
        matching reads, and the sum of ZB across the output BAM equals
        the Pass-1 scanner stat ``n_sj_blacklisted`` (Pass-1/Pass-2
        parity invariant).
        """
        import pandas as pd
        import pysam
        from rigel.sim import ReadSimConfig
        from _index_builder import rebuild_with_splice_blacklist
        from rigel.index import TranscriptIndex
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
        )
        result = scenario.build(n_fragments=500, sim_config=sim_config)

        # Inject a blacklist covering gene g1's sj (1000, 1500)
        # with huge anchor envelopes so every spliced read over it is
        # dropped.  The scenario's ref name equals scenario.name by default.
        ref = scenario.ref_name
        bl_df = pd.DataFrame(
            {
                "ref": [ref],
                "start": pd.array([1000], dtype="int32"),
                "end": pd.array([1500], dtype="int32"),
                "max_anchor_left": pd.array([1000], dtype="int32"),
                "max_anchor_right": pd.array([1000], dtype="int32"),
            }
        )
        rebuild_with_splice_blacklist(result.index_dir, bl_df)

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

        zb_values = [rec.get_tag("ZB") for rec in records if rec.has_tag("ZB")]
        n_marked = sum(1 for v in zb_values if v > 0)
        zb_sum = sum(zb_values)

        assert n_marked > 0, (
            "Expected at least one record with ZB>0 after injecting a blacklist, got none."
        )

        # Pass-1 / Pass-2 parity: the sum of ZB over the annotated BAM
        # equals the scanner's n_sj_blacklisted stat.
        assert pr.stats.n_sj_blacklisted > 0, "Pass-1 stat n_sj_blacklisted should be > 0"
        assert zb_sum == pr.stats.n_sj_blacklisted, (
            f"ZB sum ({zb_sum}) != Pass-1 n_sj_blacklisted ({pr.stats.n_sj_blacklisted})"
        )

    def test_zb_tag_absent_on_filtered_passthrough(self, scenario, tmp_path):
        """Filtered pass-through records (QCFAIL etc.) receive no Z*
        tags, including no ZB.  Augment the input BAM with a QCFAIL
        record and verify the output preserves it verbatim without tags.
        """
        import pysam
        from rigel.sim import ReadSimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
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
        with pysam.AlignmentFile(str(augmented), "wb", header=header) as fout:
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
            qc_out = [r for r in fout.fetch(until_eof=True) if r.query_name == "QCFAIL_SENTINEL"]
        assert len(qc_out) == 1, f"QCFAIL sentinel not preserved: found {len(qc_out)} copies"
        tags = {t[0] for t in qc_out[0].get_tags()}
        # Filtered pass-through records get no rigel annotations
        for tag in ("ZT", "ZG", "ZF", "ZB"):
            assert tag not in tags, f"Filtered pass-through record carries {tag}; should not"

    # =================================================================
    # ZF redesign — outcome-specific stamping tests
    # =================================================================

    def test_zf_invariants_and_zc_domain_hold_over_fixture(self, scenario, tmp_path):
        """Every annotated record satisfies the ZF invariants, and ZC is the input-ambiguity class or,
        on a template withheld with no usable alignment (ZF 0), its reason; never ``intergenic``.
        """
        import pysam
        from rigel.sim import ReadSimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
        )
        result = scenario.build(n_fragments=300, sim_config=sim_config)

        annotated_bam = tmp_path / "annotated_inv.bam"
        run_pipeline(
            result.bam_path,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
                annotated_bam_path=annotated_bam,
            ),
        )

        allowed_zc = {"unambig", "ambig_same_strand", "ambig_opp_strand", "multimapper", "."}
        withheld_zc = {"unpaired", "unannotated_splice", "off_index", "chimeric"}

        with pysam.AlignmentFile(str(annotated_bam), "rb") as fout:
            for rec in fout.fetch(until_eof=True):
                if not rec.has_tag("ZF"):
                    continue  # filtered pass-through
                zf = rec.get_tag("ZF")
                assert zf in VALID_ZF_VALUES, f"Unexpected ZF={zf:#x} on {rec.query_name}"
                _assert_zf_invariants(zf)
                zc = rec.get_tag("ZC")
                assert zc in allowed_zc or (zf == 0 and zc in withheld_zc), (
                    f"ZC={zc!r} with ZF={zf:#x} on {rec.query_name}"
                )
                assert zc != "intergenic", "ZC carries a terminal outcome ZF already states"

    @staticmethod
    def _annotate_with_intergenic_sentinel(scenario, tmp_path, *, spliced: bool):
        """Run the pipeline over the scenario's library plus one pair moved to 6,000, past every transcript
        (genome 8,000, genes end by 5,500), and return the sentinel's annotated records. The donor pair is
        spliced (a CIGAR with an N) or not, as asked."""
        from collections import defaultdict

        import pysam
        from rigel.sim import ReadSimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
        )
        result = scenario.build(n_fragments=200, sim_config=sim_config)
        with pysam.AlignmentFile(str(result.bam_path), "rb") as fin:
            header = fin.header
            records = list(fin.fetch(until_eof=True))
        by_qname = defaultdict(list)
        for r in records:
            by_qname[r.query_name].append(r)
        donor = next(
            rs
            for rs in by_qname.values()
            if len(rs) >= 2
            and not any(r.is_secondary for r in rs)
            and any("N" in (r.cigarstring or "") for r in rs) == spliced
        )
        sentinel = []
        for r in donor:
            c = r.__copy__()
            c.query_name = "INTERGENIC_SENTINEL"
            c.reference_start = 6000
            for tag in ("XS", "ts", "NH", "HI"):
                try:
                    c.set_tag(tag, None)
                except Exception:
                    pass
            sentinel.append(c)
        augmented = tmp_path / "aug_intergenic.bam"
        with pysam.AlignmentFile(str(augmented), "wb", header=header) as fout:
            for r in sentinel + records:
                fout.write(r)
        annotated_bam = tmp_path / "annotated_intergenic.bam"
        pr = run_pipeline(
            augmented,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
                annotated_bam_path=annotated_bam,
            ),
        )
        with pysam.AlignmentFile(str(annotated_bam), "rb") as fout:
            out = [r for r in fout.fetch(until_eof=True) if r.query_name == "INTERGENIC_SENTINEL"]
        assert len(out) >= 1, "Intergenic sentinel lost in output"
        return out, pr.stats

    def test_zf_intergenic_marked_gdna(self, scenario, tmp_path):
        """An unspliced read pair no transcript holds is stamped with the deterministic gDNA fast-path bits on
        ZF (is_resolved & is_gdna & is_intergenic) and ZC=".", and counted as unspliced intergenic."""
        sentinel, stats = self._annotate_with_intergenic_sentinel(scenario, tmp_path, spliced=False)
        for rec in sentinel:
            zf = rec.get_tag("ZF")
            assert zf == AF_GDNA_INTERGENIC, (
                f"Intergenic record should have ZF=AF_GDNA_INTERGENIC "
                f"({AF_GDNA_INTERGENIC:#x}), got {zf:#x}"
            )
            _assert_zf_invariants(zf)
            assert zf & AF_GDNA_BIT, "Intergenic should set is_gdna bit"
            assert zf & AF_INTERGENIC_BIT, "Intergenic should set is_intergenic bit"
            assert rec.get_tag("ZC") == "."
        assert stats.n_intergenic >= 1

    def test_zf_spliced_intergenic_withheld(self, scenario, tmp_path):
        """A spliced read pair that no annotated transcript holds is withheld (the annotated-path rule): ZF 0, ZC
        ``unannotated_splice``, and nothing counted intergenic."""
        sentinel, stats = self._annotate_with_intergenic_sentinel(scenario, tmp_path, spliced=True)
        for rec in sentinel:
            assert rec.get_tag("ZF") == AF_UNRESOLVED
            assert rec.get_tag("ZC") == "unannotated_splice"
        assert stats.n_withheld_unannotated_splice >= 1

    def test_zf_multimapper_dropped_without_include_multimap(self, scenario, tmp_path):
        """When --no-include-multimap drops a multimapper, ZF is
        AF_MULTIMAPPER_DROP and ZC is '.'.
        """
        import pysam
        from collections import defaultdict
        from rigel.sim import ReadSimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
        )
        result = scenario.build(n_fragments=100, sim_config=sim_config)

        # Augment: turn one qname group into an NH=2 multimapper.
        augmented = tmp_path / "aug_mm_drop.bam"
        with pysam.AlignmentFile(str(result.bam_path), "rb") as fin:
            header = fin.header
            records = list(fin.fetch(until_eof=True))

        by_qname = defaultdict(list)
        for r in records:
            by_qname[r.query_name].append(r)
        victim_qname = next(qn for qn, rs in by_qname.items() if len(rs) >= 2)
        out_records = []
        for qn, recs in by_qname.items():
            if qn == victim_qname:
                group = []
                for r in recs:
                    r.set_tag("NH", 2)
                    group.append(r)
                for r in recs:
                    sec = r.__copy__()
                    sec.flag = r.flag | 0x100
                    sec.reference_start = max(0, r.reference_start + 500)
                    sec.set_tag("NH", 2)
                    group.append(sec)
                out_records.extend(group)
            else:
                out_records.extend(recs)

        with pysam.AlignmentFile(str(augmented), "wb", header=header) as fout:
            for r in out_records:
                fout.write(r)

        annotated_bam = tmp_path / "annotated_mm_drop.bam"
        run_pipeline(
            augmented,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(
                    sj_strand_tag="ts",
                    include_multimap=False,
                ),
                annotated_bam_path=annotated_bam,
            ),
        )

        with pysam.AlignmentFile(str(annotated_bam), "rb") as fout:
            victim_out = [r for r in fout.fetch(until_eof=True) if r.query_name == victim_qname]

        assert len(victim_out) >= 1, "Dropped-multimapper qname lost from output"
        for rec in victim_out:
            zf = rec.get_tag("ZF")
            assert zf == AF_MULTIMAPPER_DROP, (
                f"Dropped multimapper should have ZF=AF_MULTIMAPPER_DROP "
                f"({AF_MULTIMAPPER_DROP:#x}), got {zf:#x}"
            )
            _assert_zf_invariants(zf)
            assert not (zf & AF_RESOLVED)
            assert rec.get_tag("ZC") == "."

    def test_zf_names_the_category_of_the_winning_transcript(self, winner_kinds_scenario, tmp_path):
        """ZF names a transcript's win as the output tables report its RNA: every reference transcript's
        fragments read annotated (0x03), single-exon included — alone, or covering its sibling's span, both
        0x09 before — and a synthetic span's read nascent (0x11, 0x19 before); no record sets the retired
        bit. The single-exon and span fragments are unspliced, so the EM decides them; the multi-exon
        transcripts' spliced ones are also decided on the deterministic path, so both stamp sites are read.
        """
        import pysam
        from rigel.sim import ReadSimConfig
        from rigel.config import EMConfig, PipelineConfig, BamScanConfig
        from rigel.pipeline import run_pipeline

        sim_config = ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=42,
        )
        result = winner_kinds_scenario.build(
            n_fragments=600, sim_config=sim_config, nrna_abundance=60
        )

        annotated_bam = tmp_path / "annotated_kinds.bam"
        run_pipeline(
            result.bam_path,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
                annotated_bam_path=annotated_bam,
            ),
        )

        t_df = result.index.t_df
        spans = set(t_df.loc[t_df["is_synthetic"].astype(bool), "t_id"])
        annotated = {"t_multi", "t_sibling", "t_cover", "t_single"}
        stamped = defaultdict(set)  # winning transcript -> the ZF values its records carry
        # winners the deterministic path stamped: only it writes ZC unambig with ZS spliced_annot
        deterministic = set()
        with pysam.AlignmentFile(str(annotated_bam), "rb") as fout:
            for rec in fout.fetch(until_eof=True):
                if not rec.has_tag("ZF"):
                    continue  # filtered pass-through
                zf = rec.get_tag("ZF")
                _assert_zf_invariants(zf)
                if rec.get_tag("ZI") >= 0:
                    stamped[rec.get_tag("ZT")].add(zf)
                    if (rec.get_tag("ZC"), rec.get_tag("ZS")) == ("unambig", "spliced_annot"):
                        deterministic.add(rec.get_tag("ZT"))

        won = set(stamped)
        assert won <= annotated | spans, f"unexpected winners: {sorted(won - annotated - spans)}"
        assert {"t_cover", "t_single"} <= won, (
            f"a single-exon transcript won nothing: {sorted(won)}"
        )
        assert won & spans, f"no fragment was won by a synthetic span: {sorted(won)}"
        assert deterministic & {"t_multi", "t_sibling"}, "the deterministic path stamped no winner"
        for t_id, values in stamped.items():
            expected = AF_NASCENT if t_id in spans else AF_ANNOTATED
            assert values == {expected}, (
                f"{t_id} won records stamped {sorted(hex(v) for v in values)}, not {expected:#x}"
            )


@pytestmark_scenario
class TestAnnotatedBamWriterRobustness:
    """What an aligner's or an annotation's ordinary output must not break in the annotated-BAM pass."""

    _DOCUMENTED_ZS = {
        "spliced_annot",
        "spliced_implicit",
        "splice_artifact",
        "unspliced",
        "unknown",
    }

    @staticmethod
    def _sim(scenario, n_fragments=300):
        from rigel.sim import ReadSimConfig

        return scenario.build(
            n_fragments=n_fragments,
            sim_config=ReadSimConfig(
                frag_mean=200,
                frag_std=30,
                frag_min=80,
                frag_max=450,
                read_length=100,
                strand_specificity=1.0,
                seed=42,
            ),
        )

    @staticmethod
    def _run(bam_path, index, annotated_bam, **kwargs):
        from rigel.config import BamScanConfig, EMConfig, PipelineConfig
        from rigel.pipeline import run_pipeline

        return run_pipeline(
            bam_path,
            index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
                annotated_bam_path=annotated_bam,
            ),
            **kwargs,
        )

    @staticmethod
    def _rewrite(src, dst, edit):
        """Copy a BAM record by record, in order, through ``edit(record, seen_names)``."""
        import pysam

        seen = set()
        with pysam.AlignmentFile(str(src), "rb") as fin:
            with pysam.AlignmentFile(str(dst), "wb", header=fin.header) as fout:
                for r in fin.fetch(until_eof=True):
                    edit(r, seen)
                    seen.add(r.query_name)
                    fout.write(r)

    def test_a_gene_without_a_name_is_annotated_with_a_dot(self, scenario, tmp_path):
        """A GTF need not name its genes (Ensembl leaves an unnamed gene without ``gene_name``; NCBI never writes it),
        so the index holds NaN there. The writer takes strings: it raised a TypeError after scan, calibration and EM,
        and the run was lost."""
        import pysam

        result = self._sim(scenario)
        g_df = result.index.g_df
        g_df["g_name"] = g_df["g_name"].astype(object)
        g_df.loc[g_df["g_id"] == "g1", "g_name"] = np.nan
        out = tmp_path / "annotated.bam"
        self._run(result.bam_path, result.index, out)
        with pysam.AlignmentFile(str(out), "rb") as f:
            names = {
                r.get_tag("ZR")
                for r in f.fetch(until_eof=True)
                if r.has_tag("ZG") and r.get_tag("ZG") == "g1"
            }
        assert names == {"."}

    def test_a_read_name_whose_first_record_is_filtered_stays_one_run(self, scenario, tmp_path):
        """A filtered record (here unmapped) is written through as it is read. When it opens a new read name while
        the previous read name's records are still held, writing it first split its own read name around the held
        group, and the output was no longer collated."""
        import pysam

        result = self._sim(scenario)
        src = tmp_path / "unmapped_first.bam"

        def unmap_first_mate_of_every_other_name(r, seen):
            if r.query_name not in seen and hash(r.query_name) % 2 == 0:
                r.flag |= 0x4

        self._rewrite(result.bam_path, src, unmap_first_mate_of_every_other_name)
        out = tmp_path / "annotated.bam"
        self._run(src, result.index, out)
        with pysam.AlignmentFile(str(out), "rb") as f:
            names = [r.query_name for r in f.fetch(until_eof=True)]
        runs = sum(1 for i, name in enumerate(names) if i == 0 or name != names[i - 1])
        assert runs == len(set(names))

    def test_an_aligners_tag_of_another_type_is_replaced(self, scenario, tmp_path):
        """HISAT2 writes ``ZS:i`` (its second-best alignment score). Rigel's ``ZS`` is a string, htslib will not
        update a tag in place across types, and the return was ignored, so the record kept HISAT2's integer."""
        import pysam

        result = self._sim(scenario)
        src = tmp_path / "zs_int.bam"
        self._rewrite(result.bam_path, src, lambda r, seen: r.set_tag("ZS", -12, value_type="i"))
        out = tmp_path / "annotated.bam"
        self._run(src, result.index, out)
        with pysam.AlignmentFile(str(out), "rb") as f:
            stamped = [r.get_tag("ZS") for r in f.fetch(until_eof=True) if r.has_tag("ZT")]
        assert stamped and set(stamped) <= self._DOCUMENTED_ZS

    def test_the_result_is_handed_over_before_the_annotated_bam_pass(self, scenario, tmp_path):
        """The CLI writes its tables in this hook, so a failure in the annotated-BAM pass cannot lose a
        finished run."""
        result = self._sim(scenario)
        out = tmp_path / "annotated.bam"
        seen = []
        self._run(
            result.bam_path,
            result.index,
            out,
            before_annotated_bam=lambda res: seen.append((out.exists(), res.estimator is not None)),
        )
        assert seen == [(False, True)] and out.exists()

    def test_the_cli_writes_its_tables_before_the_annotated_bam_pass(
        self, scenario, tmp_path, monkeypatch
    ):
        """The annotated-BAM pass re-reads the whole BAM after the quantification is done; if it fails, the finished
        tables must already be on disk rather than lost with it."""
        import rigel.annotate as annotate
        from rigel.cli import build_parser

        result = self._sim(scenario)

        def fail(*args, **kwargs):
            raise RuntimeError("the annotated-BAM pass failed")

        monkeypatch.setattr(annotate, "write_annotated_bam", fail)
        out = tmp_path / "out"
        args = build_parser().parse_args(
            [
                "quant",
                "--bam",
                str(result.bam_path),
                "--index",
                str(result.index_dir),
                "-o",
                str(out),
                "--annotated-bam",
                str(tmp_path / "annotated.bam"),
                "--sj-strand-tag",
                "ts",
            ]
        )
        with pytest.raises(RuntimeError, match="annotated-BAM pass"):
            args.func(args)
        assert (out / "quant.feather").exists() and (out / "summary.json").exists()
