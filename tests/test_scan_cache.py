"""Scan once, calibrate many times — and refuse a cache whose index has moved underneath it.

    TODO item 2 (the cached substrate)   ·   `docs/testing/testing_plan.md`

⭐ **WHY THIS EXISTS.** Scanning is the expensive step and calibration is the one being iterated on: a
real cfRNA run is index load ~8 s, BAM scan ~2 s, **calibration ~66 s** (`CARRY_FORWARD.md` §1 fact 22),
and a 5 M-fragment simulated condition costs far more than that to scan. Caching the scan's output took a
24-condition sweep from ~13 min to ~9 s on the old path. This is that, rebuilt against the S4 payload.

⛔ **THE KEY IS THE WHOLE POINT, AND IT NEEDS THREE PARTS, NOT ONE.**

* ``graph_hash`` — nodes plus the junction CSR. The payload already carries it, so the tally self-keys.
* ⭐ **a REACH digest.** `TODO.md` logs that ``reach`` is consumed by calibration and is covered by
  **neither** ``partition_hash`` **nor** ``graph_hash`` — correctly, since neither the scan nor the
  accumulator reads it — and that the gap "becomes live the moment something caches a calibration
  output". A cache *loaded against an index* is that moment. The 2026-07-30 index rebuild moved ~38 % of
  human contiguous reaches while **both** hashes stayed byte-identical, so a reach-blind key would have
  loaded a stale cache against a moved index and verified clean.
* **the scan config**, because two scans of one BAM under different settings are different tallies.

⚠ Anything derivable from the index is **rebuilt on load, never stored** — `RegionArrays.from_index` is
0.11 s and `build_boundary_flags_array` is 0.04 s against an 8.45 s index load that happens anyway.
Storing them is how a cache goes stale against the thing it describes (`CARRY_FORWARD.md` §3 trap 25).
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import numpy as np
import pytest

from rigel import scan_payload
from rigel.scan_cache import (
    REACH_COLUMNS,
    ScanCacheKeyError,
    calibration_inputs,
    index_derived_inputs,
    read_scan_cache,
    reach_digest,
    write_scan_cache,
)
from rigel.scan_payload import AccumulatorPayload


SEED = 1234


@pytest.fixture(scope="module")
def scanned(tmp_path_factory):
    """A REAL scan of a tiny oracle BAM: the payload, strand model and FL histograms it produced.

    ⚠ A real scan rather than a hand-built payload on purpose — the cache's whole job is to reproduce
    what the scanner emitted, including the per-junction strand table and the five length pools, and a
    fixture that never ran the scanner could not catch a field the writer forgot.
    """
    from rigel.config import BamScanConfig
    from rigel.pipeline import scan_and_buffer
    from rigel.sim import ReadSimConfig, Scenario

    work = tmp_path_factory.mktemp("scan_cache")
    scenario = Scenario("cache", genome_length=6000, seed=SEED, work_dir=work / "sim")
    scenario.add_gene(
        "g1",
        "+",
        [{"t_id": "t1", "exons": [(200, 400), (600, 800), (1200, 1500)], "abundance": 60}],
    )
    scenario.add_gene(
        "g2", "-", [{"t_id": "t2", "exons": [(2500, 2700), (3000, 3200)], "abundance": 40}]
    )
    scenario.add_gene("g3", "+", [{"t_id": "t3", "exons": [(4000, 4600)], "abundance": 25}])
    result = scenario.build_oracle(
        n_fragments=600,
        sim_config=ReadSimConfig(
            frag_mean=200,
            frag_std=40,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=0.99,
            seed=SEED,
        ),
    )
    scan = BamScanConfig(sj_strand_tag="auto")
    _stats, strand_model, fl_models, _buffer, payload = scan_and_buffer(
        str(result.bam_path), result.index, scan
    )
    yield result.index, str(result.bam_path), scan, payload, strand_model, fl_models
    scenario.cleanup()


def round_trip(scanned, tmp_path: Path) -> tuple:
    index, bam, scan, payload, strand_model, fl_models = scanned
    cache_dir = tmp_path / "cache"
    write_scan_cache(
        cache_dir,
        payload=payload,
        strand_model=strand_model,
        frag_length_models=fl_models,
        index=index,
        bam=bam,
        scan_config=scan,
    )
    return cache_dir, read_scan_cache(cache_dir, index)


class TestTheCachedTallyIsTheSCANNEDTally:
    """Byte-identity, field by field. ⚠ This is the gate that can run TODAY — `calibrate()` cannot,
    because `substrate.py` still reads payload fields S4 deleted, so an end-to-end calibration
    comparison is blocked on S5. Input identity is the strongest available statement."""

    def test_every_payload_array_survives_the_round_trip_byte_identical(self, scanned, tmp_path):
        _cache_dir, cache = round_trip(scanned, tmp_path)
        original: AccumulatorPayload = scanned[3]
        for field in dataclasses.fields(AccumulatorPayload):
            before = getattr(original, field.name)
            after = getattr(cache.payload, field.name)
            if isinstance(before, np.ndarray):
                assert after.dtype == before.dtype, f"{field.name} dtype moved"
                assert after.shape == before.shape, f"{field.name} shape moved"
                assert np.array_equal(after, before), f"{field.name} content moved"
            elif dataclasses.is_dataclass(before):
                assert dataclasses.asdict(after) == dataclasses.asdict(before), field.name
            else:
                assert after == before, field.name

    def test_the_strand_model_survives_including_its_per_junction_table(self, scanned, tmp_path):
        """⚠ The 2x2 is the MARGINAL of the per-junction table; the dispersion fit needs the table
        itself, so a cache that kept only the 2x2 would silently disable the dispersion estimate."""
        _cache_dir, cache = round_trip(scanned, tmp_path)
        original = scanned[4]
        for name in ("exonic_spliced", "exonic"):
            before, after = getattr(original, name), getattr(cache.strand_model, name)
            assert (after.pos_pos, after.pos_neg, after.neg_pos, after.neg_neg) == (
                before.pos_pos,
                before.pos_neg,
                before.neg_pos,
                before.neg_neg,
            ), f"{name} 2x2 moved"
            if before.sj_table is None:
                assert after.sj_table is None, f"{name} invented an sj_table"
            else:
                for column in ("ref_id", "start", "end", "motif_strand", "n_sense", "n_antisense"):
                    assert np.array_equal(
                        getattr(after.sj_table, column), getattr(before.sj_table, column)
                    ), f"{name}.sj_table.{column} moved"

    def test_the_fl_histograms_are_cached_RAW_not_as_derived_pmfs(self, scanned, tmp_path):
        """⭐ `build_fl_models` must stay the single source of truth. Freezing its output into the
        cache would mean a change to the FL model silently does not apply to cached scans."""
        _cache_dir, cache = round_trip(scanned, tmp_path)
        fl_models = scanned[5]
        from rigel.splice import SpliceType

        assert np.array_equal(cache.fl_global_counts, fl_models.global_model.counts)
        assert np.array_equal(
            cache.fl_rna_counts, fl_models.category_models[SpliceType.SPLICED_ANNOT].counts
        )
        assert cache.fl_max_size == fl_models.max_size


class TestTheKeyRefusesAMovedIndex:
    def test_a_changed_graph_hash_is_refused(self, scanned, tmp_path):
        cache_dir, _cache = round_trip(scanned, tmp_path)
        manifest = json.loads((cache_dir / "manifest.json").read_text())
        manifest["graph_hash"] = "0" * 16
        (cache_dir / "manifest.json").write_text(json.dumps(manifest))
        with pytest.raises(ScanCacheKeyError, match="graph_hash"):
            read_scan_cache(cache_dir, scanned[0])

    def test_a_changed_REACH_is_refused(self, scanned, tmp_path):
        """⭐ The one neither existing hash covers. A rebuild moved ~38 % of human contiguous reaches
        with `partition_hash` AND `graph_hash` byte-identical — this is the check that notices."""
        cache_dir, _cache = round_trip(scanned, tmp_path)
        manifest = json.loads((cache_dir / "manifest.json").read_text())
        manifest["reach_digest"] = "0" * 16
        (cache_dir / "manifest.json").write_text(json.dumps(manifest))
        with pytest.raises(ScanCacheKeyError, match="reach"):
            read_scan_cache(cache_dir, scanned[0])

    def test_a_cache_from_a_DIFFERENT_ACCUMULATOR_SCHEMA_is_refused(self, scanned, tmp_path):
        """⛔ The gap none of the other three keys covers: the ACCUMULATOR's own field list.

        `graph_hash` describes the index, `reach_digest` the reaches, `scan_config_digest` the scan
        settings — none of them moves when the accumulator changes. S5.a made that concrete by adding
        `length_sum` to every population; without this key a cache written the day before would be
        accepted and then die deep in `_payload_from_parts` with a bare `KeyError`, which reads as a bug
        rather than as a stale cache.
        """
        cache_dir, _cache = round_trip(scanned, tmp_path)
        manifest = json.loads((cache_dir / "manifest.json").read_text())
        manifest["payload_schema_digest"] = "0" * 16
        (cache_dir / "manifest.json").write_text(json.dumps(manifest))
        with pytest.raises(ScanCacheKeyError, match="schema"):
            read_scan_cache(cache_dir, scanned[0])

    def test_a_cache_predating_the_key_ENTIRELY_is_refused(self, scanned, tmp_path):
        """A manifest written before the key existed has no such entry at all. Absent must not read as
        agreement — that is how a stale cache verifies clean (`CARRY_FORWARD.md` §3 trap 25)."""
        cache_dir, _cache = round_trip(scanned, tmp_path)
        manifest = json.loads((cache_dir / "manifest.json").read_text())
        del manifest["payload_schema_digest"]
        (cache_dir / "manifest.json").write_text(json.dumps(manifest))
        with pytest.raises(ScanCacheKeyError, match="schema"):
            read_scan_cache(cache_dir, scanned[0])

    def test_the_schema_digest_moves_when_a_PAYLOAD_FIELD_moves(self):
        """Teeth on the digest itself: it must be a function of the field list, not a constant."""
        from rigel.scan_cache import payload_schema_digest

        before = payload_schema_digest()
        original = scan_payload.AccumulatorPayload.__dataclass_fields__
        trimmed = {k: v for k, v in original.items() if k != "sj_length_sum"}
        try:
            scan_payload.AccumulatorPayload.__dataclass_fields__ = trimmed
            assert payload_schema_digest() != before, "dropping a field must move the digest"
        finally:
            scan_payload.AccumulatorPayload.__dataclass_fields__ = original
        assert payload_schema_digest() == before

    def test_a_changed_scan_config_is_refused(self, scanned, tmp_path):
        cache_dir, _cache = round_trip(scanned, tmp_path)
        manifest = json.loads((cache_dir / "manifest.json").read_text())
        manifest["scan_config_digest"] = "0" * 16
        (cache_dir / "manifest.json").write_text(json.dumps(manifest))
        with pytest.raises(ScanCacheKeyError, match="scan"):
            read_scan_cache(cache_dir, scanned[0])

    @pytest.mark.parametrize("column", REACH_COLUMNS)
    def test_the_reach_digest_depends_on_EVERY_reach_column(self, scanned, column):
        """⛔ Teeth on the digest itself, one column at a time.

        The first version perturbed only `reach_lo_pos`, so a digest that covered just the POS pair
        passed — half the reach unguarded. Per strand and per side, each column must move it.
        """
        index = scanned[0]
        before = reach_digest(index)
        original = index.edges_df[column].copy()
        try:
            index.edges_df.loc[index.edges_df.index[0], column] = int(original.iloc[0]) + 1
            assert reach_digest(index) != before, f"{column} does not reach the digest"
        finally:
            index.edges_df[column] = original
        assert reach_digest(index) == before

    def test_a_cache_whose_graph_hash_disagrees_with_the_INDEX_is_refused(self, scanned, tmp_path):
        """⛔ Distinct from the manifest-tampering test above, which trips the self-consistency check
        first and so never reaches the index comparison. Move the manifest AND the payload together —
        a cache that is internally consistent but describes a different index."""
        cache_dir, _cache = round_trip(scanned, tmp_path)
        manifest = json.loads((cache_dir / "manifest.json").read_text())
        manifest["graph_hash"] = "0" * 16
        manifest["payload_scalars"]["graph_hash"] = "0" * 16
        (cache_dir / "manifest.json").write_text(json.dumps(manifest))
        with pytest.raises(ScanCacheKeyError, match="index graph_hash"):
            read_scan_cache(cache_dir, scanned[0])


class TestNothingDerivableFromTheIndexIsStored:
    def test_the_cache_holds_no_region_arrays_or_boundary_flags(self, scanned, tmp_path):
        """They are 0.15 s to rebuild and a stale copy of them is a silent wrong answer."""
        cache_dir, _cache = round_trip(scanned, tmp_path)
        stored = {p.name for p in cache_dir.iterdir()}
        assert not any("region" in name or "boundary" in name for name in stored), stored

    def test_the_index_derived_inputs_are_rebuilt_not_loaded(self, scanned, tmp_path):
        index = scanned[0]
        from rigel.calibration.region_arrays import RegionArrays
        from rigel.calibration.splice_graph import build_boundary_flags_array

        inputs = index_derived_inputs(index)
        assert inputs["region_arrays"].n_regions == RegionArrays.from_index(index).n_regions
        assert np.array_equal(inputs["boundary_flags"], build_boundary_flags_array(index))

    def test_the_index_derived_names_are_ones_calibrate_accepts(self, scanned):
        """⭐ Read off `calibrate`'s signature, never written out here — a renamed parameter fails this
        test rather than silently escaping the cache."""
        import inspect

        from rigel.calibration.calibrate import calibrate

        accepted = set(inspect.signature(calibrate).parameters)
        offered = set(index_derived_inputs(scanned[0]))
        assert offered <= accepted, f"not accepted by calibrate(): {offered - accepted}"

    def test_calibration_inputs_offers_every_argument_calibrate_needs(self, scanned, tmp_path):
        """✅ Unblocked by S5.b — `fl.py` now reads the payload's five pure pools.

        It was a strict xfail while `gdna_fl_mass` still reached for `payload.fl_pool_mass`, the field
        S4 replaced. Strict, so that it would FAIL the moment it started working rather than sit green
        and forgotten.
        """
        import inspect

        from rigel.calibration.calibrate import calibrate

        _cache_dir, cache = round_trip(scanned, tmp_path)
        inputs = calibration_inputs(cache, scanned[0])
        accepted = set(inspect.signature(calibrate).parameters)
        assert set(inputs) <= accepted, f"not accepted by calibrate(): {set(inputs) - accepted}"
        required = {"payload", "region_arrays", "strand_model", "gdna_fl_pmf", "rna_fl_pmf"}
        assert required <= set(inputs), f"missing: {required - set(inputs)}"
        assert inputs["gdna_fl_pmf"].sum() == pytest.approx(1.0, abs=1e-9)


@pytest.mark.xfail(
    reason=(
        "STEP 4 (the population-prior seed) is blocked on S5. Extracting "
        "InjectedCalibrationPriors requires calibrate() to run, and substrate.py still reads "
        "payload.region_contained / boundary_flux_*, which S4 deleted. This test documents the "
        "dependency rather than pretending the seed path is verified."
    ),
    strict=True,
)
def test_population_priors_can_be_extracted_from_a_cached_scan(scanned, tmp_path):
    from rigel.calibration.calibrate import calibrate
    from rigel.config import CalibrationConfig

    _cache_dir, cache = round_trip(scanned, tmp_path)
    debug: dict = {}
    calibrate(**calibration_inputs(cache, scanned[0]), config=CalibrationConfig(), _debug=debug)
    assert debug["calibration_priors"] is not None
