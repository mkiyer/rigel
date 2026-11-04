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


@pytest.fixture
def blacklisted_oracle(tmp_path):
    """The same scenario twice: its index without a splice blacklist, and with one that names
    **every sj the simulator wrote**.

    This exists so the ``SPLICE_ARTIFACT`` term of the census identity is exercised rather than
    asserted over a zero: a gate that only ever sees 0 on the term it was built for checks nothing.
    """
    import pandas as pd
    import pysam

    from _index_builder import rebuild_with_splice_blacklist
    from rigel.index import TranscriptIndex

    sc = Scenario("bl_int", genome_length=5000, seed=SEED, work_dir=tmp_path / "bl_int")
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

    # The sj are read back out of the BAM rather than derived from the exon list above, so
    # the fixture cannot silently blacklist nothing if the simulator's coordinate convention moves.
    introns: set[tuple[str, int, int]] = set()
    with pysam.AlignmentFile(str(result.bam_path), "rb") as bam:
        for rec in bam:
            if not rec.cigartuples:
                continue
            pos = rec.reference_start
            for op, ln in rec.cigartuples:
                if op == 3:  # CIGAR N
                    introns.add((rec.reference_name, pos, pos + ln))
                if op in (0, 2, 3, 7, 8):
                    pos += ln
    assert introns, "the fixture simulated no spliced reads, so it cannot make an artifact"

    rows = sorted(introns)
    # Anchors far larger than the 100 bp reads, so every crossing read is blacklisted and the
    # promotion is not a function of where the fragment happened to start.
    rebuild_with_splice_blacklist(
        result.index_dir,
        pd.DataFrame(
            {
                "ref": [r[0] for r in rows],
                "start": np.asarray([r[1] for r in rows], dtype=np.int32),
                "end": np.asarray([r[2] for r in rows], dtype=np.int32),
                "max_anchor_left": np.full(len(rows), 500, dtype=np.int32),
                "max_anchor_right": np.full(len(rows), 500, dtype=np.int32),
            }
        ),
    )

    yield result, TranscriptIndex.load(result.index_dir)
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


class TestScannerAccumulatorIntegration:
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

    These tests use the blacklist fixture on purpose: on the plain oracle a perfect BAM with no
    ambiguity makes the competing definitions agree, so a value gate there passes whichever anchor
    is wired in.
    """

    @staticmethod
    def _scan_raw(bam_path, index):
        config = PipelineConfig(em=EMConfig(seed=SEED), scan=BamScanConfig(sj_strand_tag="auto"))
        return scan_and_buffer(str(bam_path), index, config.scan)

    def test_the_anchor_is_the_accumulators_own_deposited_histogram(self, blacklisted_oracle):
        """The anchor is ``deposited_lengths``, over exactly the population the pools are drawn from.

        "Unconditional given deposit": an anchor over a wider population would re-create the frame
        mismatch somewhere new, so the anchor's total must equal ``qc.deposited`` and not the count of
        everything the scanner classified.
        """
        from rigel.calibration.fl import build_fl_models

        result, index = blacklisted_oracle
        _, _, _, payload = self._scan_raw(result.bam_path, index)

        fl = build_fl_models(payload)
        np.testing.assert_array_equal(fl.global_counts, np.asarray(payload.deposited_lengths))
        assert float(fl.global_counts.sum()) == float(payload.qc.deposited)
        np.testing.assert_allclose(
            fl.global_pmf, np.asarray(payload.deposited_lengths) / payload.qc.deposited
        )

    # No gate compares the anchor against a scanner histogram, because the scanner builds none.

    def test_a_FOREIGN_ANCHOR_CANNOT_BE_PASSED_AT_ALL(self):
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
    them: its subject is fragment length, the scanner's is what it saw and what it held out.
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
        """The invariant: every censused fragment either deposited, was named as a rejection, was
        held for the second pass, or was named as a hold-out. Nothing is lost between the two
        subsystems.

            Σ census − census[SPLICE_ARTIFACT]
                == qc.deposited + Σ qc.dropped_* + qc.deferred_undetermined_gap

        Every fragment of this fixture has one placement, so the census, which counts placements,
        and the ``qc`` counters, which count fragments, count the same fragments.

        ``SPLICE_ARTIFACT`` is subtracted because it is censused and then **held out by the
        scanner**: a blacklisted CIGAR-N sj may be a real-but-rejected sj or a wholly
        incorrect alignment, so the span that would be deposited is derived from an alignment the
        scanner has already refused to believe. Identifying and filtering those is the scanner's
        job, and the census is where that decision becomes visible.

        This is the same externally-checkable form as ``sum(deposited_lengths) == qc.deposited`` and
        a different statement: that one says every deposited fragment was binned by length, this one
        says every classified fragment was accounted for on its way to the deposit.
        """
        from rigel.splice import SpliceType, census_field

        stats, payload = _scan_full(oracle)
        census = {st: int(getattr(stats, census_field(st))) for st in SpliceType}
        qc = payload.qc

        offered = sum(census.values())
        accounted = (
            int(qc.deposited)
            + int(qc.dropped_too_long)
            + int(qc.dropped_empty)
            + int(qc.dropped_strand_undefined)
            + int(qc.deferred_undetermined_gap)
        )
        assert offered - census[SpliceType.SPLICE_ARTIFACT] == accounted, (
            f"census {census} (artifacts held out: {census[SpliceType.SPLICE_ARTIFACT]}) does not "
            f"reconcile with deposited={qc.deposited} dropped=("
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

    def test_the_identity_holds_when_fragments_are_held_out_as_artifacts(self, blacklisted_oracle):
        """The identity, on a library where the hold-out term is the dominant one.

        Without this the ``− census[SPLICE_ARTIFACT]`` term is asserted only over zero, and any
        mistake in it (counting artifacts as deposited, as not-offered, or not at all) is invisible.
        """
        from rigel.splice import SpliceType, census_field

        result, index = blacklisted_oracle
        stats, payload = _scan_full(result, index)
        census = {st: int(getattr(stats, census_field(st))) for st in SpliceType}
        qc = payload.qc

        assert census[SpliceType.SPLICE_ARTIFACT] > 0, (
            "the blacklist held nothing out, so this test proves nothing about artifacts"
        )
        offered = sum(census.values())
        accounted = (
            int(qc.deposited)
            + int(qc.dropped_too_long)
            + int(qc.dropped_empty)
            + int(qc.dropped_strand_undefined)
            + int(qc.deferred_undetermined_gap)
        )
        assert offered - census[SpliceType.SPLICE_ARTIFACT] == accounted

    def test_the_artifact_census_is_confirmed_by_the_deposit_TOTAL(
        self, oracle, blacklisted_oracle
    ):
        """The artifact count, derived a second way, from a subsystem that has never heard of an
        artifact.

        The same BAM is scanned twice against the same annotation, differing only in whether the
        index carries a splice blacklist. Blacklisting relabels fragments ``SPLICED_ANNOT →
        SPLICE_ARTIFACT``; the scanner then holds exactly those out. So the accumulator's own
        ``qc.deposited``, which knows nothing of splice types, blacklists or artifacts, must fall
        by **exactly** the artifact census:

            deposited(no blacklist) − deposited(blacklisted) == census[SPLICE_ARTIFACT]

        This crosses a subsystem boundary, so a census placed after the hold-out return, reading a
        confident zero for the very category the report names, cannot survive it.

        The two supporting equalities are asserted as well, because the subtraction is only
        attributable to the relabel if nothing else moved: the census total is unchanged (the same
        fragments were classified), and the artifacts came out of ``SPLICED_ANNOT``.
        """
        from rigel.splice import SpliceType, census_field

        clean_stats, clean_payload = _scan_full(oracle)
        result, index = blacklisted_oracle
        bl_stats, bl_payload = _scan_full(result, index)

        clean = {st: int(getattr(clean_stats, census_field(st))) for st in SpliceType}
        bl = {st: int(getattr(bl_stats, census_field(st))) for st in SpliceType}

        assert sum(bl.values()) == sum(clean.values()), (
            f"the blacklist changed how many fragments were CLASSIFIED ({sum(bl.values())} vs "
            f"{sum(clean.values())}); the deposit difference is then not attributable to the relabel"
        )
        n_artifact = bl[SpliceType.SPLICE_ARTIFACT]
        assert n_artifact > 0
        assert clean[SpliceType.SPLICED_ANNOT] - bl[SpliceType.SPLICED_ANNOT] == n_artifact, (
            f"artifacts did not come out of SPLICED_ANNOT: {clean} -> {bl}"
        )
        assert int(clean_payload.qc.deposited) - int(bl_payload.qc.deposited) == n_artifact, (
            f"the scanner censused {n_artifact} artifacts but the accumulator's deposit total fell "
            f"by {int(clean_payload.qc.deposited) - int(bl_payload.qc.deposited)}"
        )

    # The census is a superset of the fragment-length observation population, and no separate
    # observation counter exists to compare it against; the population statement is the identity
    # above: every censused fragment either deposits, is a named rejection, is held for the second
    # pass, or is a named hold-out.
