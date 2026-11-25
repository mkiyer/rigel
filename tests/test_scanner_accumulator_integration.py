"""End-to-end check that BamScanner populates an AccumulatorPayload.

Builds a tiny oracle scenario, calls ``scan_and_buffer`` directly, and
verifies that the returned ``AccumulatorPayload`` has the expected shape
derived from the index's region partition and that fragments deposited
non-trivial fractional mass into at least one region or boundary.
"""

from __future__ import annotations


import numpy as np
import pytest

from rigel.calibration.splice_graph import build_region_partition_arrays
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import scan_and_buffer
from rigel.scan_payload import AccumulatorPayload
from rigel.sim import ReadSimConfig, Scenario

SEED = 1234


@pytest.fixture
def oracle(tmp_path):
    sc = Scenario(
        "acc_int",
        genome_length=5000,
        seed=SEED,
        work_dir=tmp_path / "acc_int",
    )
    sc.add_gene(
        "g1",
        "+",
        [{"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 50}],
    )
    sc.add_gene(
        "g2",
        "-",
        [{"t_id": "t2", "exons": [(2500, 2700), (3000, 3200)], "abundance": 30}],
    )
    sim_config = ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=1.0,
        seed=SEED,
    )
    result = sc.build_oracle(n_fragments=200, sim_config=sim_config)
    yield result
    sc.cleanup()


def _scan_full(result, index=None):
    """The scan's ``(stats, payload)`` — the two QC objects the census identity spans."""
    config = PipelineConfig(
        em=EMConfig(seed=SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    stats, _, _, payload = scan_and_buffer(
        str(result.bam_path), result.index if index is None else index, config.scan
    )
    return stats, payload


def _scan(result) -> AccumulatorPayload:
    return _scan_full(result)[1]


def _accounted(qc) -> int:
    """Every fragment the accumulator was offered, by its outcome: deposited, dropped by name, or held for the second
    pass on its gap or on its placements."""
    return (
        int(qc.deposited)
        + int(qc.dropped_too_long)
        + int(qc.dropped_empty)
        + int(qc.dropped_strand_undefined)
        + int(qc.deferred_undetermined_gap)
        + int(qc.deferred_multiple_placements)
    )


class TestScannerAccumulatorIntegration:
    def test_the_scan_refuses_to_run_without_the_region_partition(self, oracle, monkeypatch):
        """The calibration tally is not optional: a scan whose regions were never installed is refused, not run with
        an empty tally."""
        import rigel.pipeline

        monkeypatch.setattr(rigel.pipeline, "_wire_calibration_regions", lambda *a, **k: None)
        with pytest.raises(RuntimeError, match="set_regions was not called"):
            _scan(oracle)

    def test_payload_is_populated(self, oracle):
        payload = _scan(oracle)
        assert payload is not None
        assert isinstance(payload, AccumulatorPayload)

    def test_payload_shape_matches_index_partition(self, oracle):
        """The payload's three axes must be exactly what the index's partition implies.

        ``region_bounds`` are the region-bound positions: a reference with ``k`` of them owns
        ``k − 1`` regions and ``k − 2`` interior boundaries. Counting ``k`` boundary objects per
        reference instead (the interiors plus two data-free terminals) gives an axis two entries
        too long, and every per-boundary array then addresses the wrong rows.
        """
        index = oracle.index
        region_bounds, ref_region_bound_offsets, region_types = build_region_partition_arrays(index)
        payload = _scan(oracle)

        np.testing.assert_array_equal(payload.region_bounds, region_bounds)
        np.testing.assert_array_equal(payload.ref_region_bound_offsets, ref_region_bound_offsets)
        assert payload.n_refs == len(index.ref_names)
        assert region_types.shape == (payload.n_regions,)  # one type per region

        diffs = np.diff(ref_region_bound_offsets)
        expected_regions = int(np.sum(np.maximum(diffs - 1, 0)))
        expected_boundaries = int(np.sum(np.maximum(diffs - 2, 0)))
        assert payload.n_regions == expected_regions
        assert payload.n_boundaries == expected_boundaries
        # boundaries = regions − (non-empty refs), stated a second way: the two derivations must agree.
        n_live_refs = int(np.sum(diffs > 1))
        assert payload.n_boundaries == payload.n_regions - n_live_refs

    def test_fl_pools_emitted(self, oracle):
        """The scan emits the five pure fragment-length pools, binned at the same L as every other
        bank. The pools are integer counts on a ``(N_FRAGMENT_POOLS, max_length + 1)`` grid: there
        is no ``fl_pool_mass`` and no separate ``fl_max_size``, because a pool is a histogram of the
        same molecule length the accumulator deposits by."""
        from rigel.calibration.fl import gdna_fl_mass
        from rigel.scan_payload import N_FRAGMENT_POOLS

        payload = _scan(oracle)
        assert payload.pool_lengths.shape == (N_FRAGMENT_POOLS, payload.max_length + 1)
        # gDNA pool aggregation is well-formed (non-negative, right length).
        gdna = gdna_fl_mass(payload)
        assert gdna.shape == (payload.max_length + 1,)
        assert float(gdna.sum()) >= 0.0

    def test_at_least_some_mass_deposited(self, oracle):
        payload = _scan(oracle)
        # One tally answers this: region_start_count is incremented once per accepted fragment, so its
        # total is the deposit count.
        assert int(np.asarray(payload.region_start_count).sum()) > 0, "scanner deposited nothing"
        assert int(payload.qc.deposited) == int(np.asarray(payload.region_start_count).sum())


class TestFragmentLengthAnchor:
    """The empirical-Bayes anchor is the accumulator's unconditional histogram.

    ``build_fl_models`` EB-shrinks the accumulator's pure pools toward an anchor. An anchor taken
    from the scanner would measure fragment length by another rule over another population, so
    accumulator-frame pools would be shrunk toward a scanner-frame anchor.

    The oracle classifies fragments the accumulator does not deposit (those whose mate gap may hold an
    annotated intron are held for the second pass), so the competing definitions disagree on it and a
    value gate tells them apart.
    """

    def test_the_anchor_is_the_accumulators_own_deposited_histogram(self, oracle):
        """The anchor is ``deposited_lengths``, over exactly the population the pools are drawn from.

        "Unconditional given deposit": an anchor over a wider population would re-create the frame
        mismatch somewhere new, so the anchor's total must equal ``qc.deposited`` and not the count of
        everything the scanner classified.
        """
        from rigel.pipeline import library_fl_models
        from rigel.splice import SpliceType, census_field

        stats, payload = _scan_full(oracle)
        classified = sum(int(getattr(stats, census_field(st))) for st in SpliceType)
        assert classified > payload.qc.deposited, "every classified fragment deposited: no teeth"

        fl = library_fl_models(payload, oracle.index, RegionArrays.from_index(oracle.index))
        np.testing.assert_array_equal(fl.global_counts, np.asarray(payload.deposited_lengths))
        assert float(fl.global_counts.sum()) == float(payload.qc.deposited)
        np.testing.assert_allclose(
            fl.global_pmf, np.asarray(payload.deposited_lengths) / payload.qc.deposited
        )

    # No gate compares the anchor against a scanner histogram, because the scanner builds none.

    def test_a_foreign_anchor_cannot_be_passed_at_all(self):
        """The gate here is structural, not a value.

        The two tests above pin what ``build_fl_models`` returns; neither can catch a call site that
        hands it the wrong array, which is the mixed-frame defect at issue. So the public entry point
        takes the payload and all three histograms are read off that one object in that one frame:
        there is no ``global_counts`` parameter to get wrong.
        """
        import inspect

        from rigel.calibration.fl import build_fl_models

        params = inspect.signature(build_fl_models).parameters
        assert "global_counts" not in params, (
            "build_fl_models still accepts a raw anchor, so a caller can still mix frames"
        )
        assert "payload" in params
        with pytest.raises(TypeError):
            build_fl_models(global_counts=np.zeros(11), rna_counts=np.zeros(11))


class TestSpliceCensus:
    """The per-fragment splice breakdown is scanner QC, and it closes the books.

    Reading `rigel report`'s splice-type counts off the fragment-length category models would count
    only the fragments that contributed a length observation, a population gated by the
    transcript-space unanimity test and by the single-block rule on the intergenic path. The counts
    belong where the classification is made, and nothing is routed through the accumulator to obtain
    them: its subject is fragment length, the scanner's is what it saw.
    """

    def test_every_splice_type_is_censused(self, oracle):
        """Every :class:`SpliceType` reaches Python, by a name derived in both languages.

        The failure mode this exists for is silent: ``pipeline`` copies the scanner's dict onto
        ``PipelineStats`` with ``stats_dict.get(key, 0)``, so a C++ key that does not match a Python
        field reads **zero** rather than raising. A category added to the enum and forgotten in
        either language would be reported as "none of those were seen".
        """
        from rigel.splice import SpliceType, census_field

        stats, _ = _scan_full(oracle)
        for stype in SpliceType:
            assert hasattr(stats, census_field(stype)), (
                f"{stype.name} has no census field; PipelineStats and the C++ "
                f"splice_type_label table have drifted"
            )

    def test_the_census_accounts_for_every_fragment_offered_to_the_accumulator(self, oracle):
        """The invariant: every censused fragment either deposited, was named as a rejection, or was
        held for the second pass. Nothing is lost between the two subsystems.

            Σ census == qc.deposited + Σ qc.dropped_* + qc.deferred_undetermined_gap
                    + qc.deferred_multiple_placements

        Both count fragments: the census counts a template once, whatever its placements, as the accumulator's
        offer does (:meth:`test_a_multimapper_is_censused_once` holds it where they differ).

        This is the same externally-checkable form as ``sum(deposited_lengths) == qc.deposited`` and
        a different statement: that one says every deposited fragment was binned by length, this one
        says every classified fragment was accounted for on its way to the deposit.
        """
        from rigel.splice import SpliceType, census_field

        stats, payload = _scan_full(oracle)
        census = {st: int(getattr(stats, census_field(st))) for st in SpliceType}
        qc = payload.qc

        offered = sum(census.values())
        assert offered == _accounted(qc), (
            f"census {census} does not reconcile with deposited={qc.deposited} dropped=("
            f"too_long={qc.dropped_too_long}, empty={qc.dropped_empty}, "
            f"strand_undefined={qc.dropped_strand_undefined}, "
            f"deferred_undetermined_gap={qc.deferred_undetermined_gap})"
        )

    def test_the_census_has_teeth_on_this_fixture(self, oracle):
        """The identity above is satisfiable by an all-zero census. This fixture must populate at
        least two categories, so that dropping any single category's increment is detectable."""
        from rigel.splice import SpliceType, census_field

        stats, _ = _scan_full(oracle)
        census = {st.name: int(getattr(stats, census_field(st))) for st in SpliceType}
        assert sum(census.values()) > 0, f"nothing was censused: {census}"
        assert sum(1 for n in census.values() if n > 0) >= 2, (
            f"only one splice category is populated, so this fixture cannot detect a dropped "
            f"increment in the others: {census}"
        )

    def test_a_multimapper_is_censused_once(self, tmp_path):
        """The census counts fragments, as the accumulator's ``qc`` block does: a template with two usable hits is one
        fragment, censused once, as its first hit with a candidate (the hit ``with_annotated_sj`` reads), else its
        first hit. So the census closes against the ``qc`` block with multimappers in the library, and
        ``spliced_annotated`` and ``with_annotated_sj`` count one unit.

        A two-exon transcript on each of two references; 20 templates spliced on both, 10 unspliced on both, 5 whose
        first hit is intergenic and second spliced, and 5 unique unspliced pairs: 40 fragments on 75 placements."""
        from _bam_builder import R1_FWD, R2_REV, SECONDARY, record, write_bam
        from _index_builder import build_index

        from rigel.splice import SpliceType, census_field

        refs = {"chr1": 4000, "chr2": 4000}
        gtf = "".join(
            f'{ref}\tt\texon\t{s}\t{e}\t.\t+\t.\tgene_id "g_{ref}"; transcript_id "t_{ref}";\n'
            for ref in refs
            for s, e in ((1001, 1600), (2001, 2200))
        )
        index = build_index(tmp_path, gtf, refs)

        def pair(name, ref, hi, f_pos, f_cigar, r_pos, nh=2):
            extra = SECONDARY if hi == 2 else 0
            return [
                record(
                    name, R1_FWD | extra, f_pos, f_cigar, r_pos, ref=ref, nh=nh, hi=hi, nm=0, xs="+"
                ),
                record(
                    name, R2_REV | extra, r_pos, "100M", f_pos, ref=ref, nh=nh, hi=hi, nm=0, xs="+"
                ),
            ]

        reads = []
        for k in range(20):  # read 1 crosses the intron [1600, 2000) on both references
            reads += pair(f"spliced{k}", 0, 1, 1550, "50M400N50M", 2050)
            reads += pair(f"spliced{k}", 1, 2, 1550, "50M400N50M", 2050)
        for k in range(10):
            reads += pair(f"unspliced{k}", 0, 1, 1100, "100M", 1300)
            reads += pair(f"unspliced{k}", 1, 2, 1100, "100M", 1300)
        for k in range(5):  # the first hit lies past the transcript, so it has no candidate
            reads += pair(f"mixed{k}", 0, 1, 3000, "100M", 3200)
            reads += pair(f"mixed{k}", 1, 2, 1550, "50M400N50M", 2050)
        for k in range(5):
            reads += pair(f"unique{k}", 0, 1, 1100, "100M", 1300, nh=1)
        bam = write_bam(tmp_path / "nh2.bam", reads, refs)
        stats, _, _, payload = scan_and_buffer(str(bam), index, BamScanConfig(total_threads=1))

        census = {st: int(getattr(stats, census_field(st))) for st in SpliceType}
        assert stats.n_placed_multimapping == 35 and stats.n_placed_unique == 5
        assert sum(census.values()) == 40, f"the census counted placements: {census}"
        assert census[SpliceType.SPLICED_ANNOT] == 25 == stats.n_with_annotated_sj
        assert census[SpliceType.UNSPLICED] == 15
        assert payload.qc.deferred_multiple_placements > 0, "no template was held on its placements"
        assert sum(census.values()) == _accounted(payload.qc)

    # The census is a superset of the fragment-length observation population, and no separate
    # observation counter exists to compare it against; the population statement is the identity
    # above: every censused fragment either deposits, is a named rejection, or is held for the second
    # pass.
