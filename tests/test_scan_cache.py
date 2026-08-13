"""Scan once, calibrate many times — and refuse a cache whose index has moved underneath it.

    TODO item 2 (the cached substrate)

⭐ **WHY THIS EXISTS.** Scanning is the expensive step and calibration is the one being iterated on: a
real cfRNA run is index load ~8 s, BAM scan ~2 s, **calibration ~66 s**,
and a 5 M-fragment simulated condition costs far more than that to scan. Caching the scan's output took a
24-condition sweep from ~13 min to ~9 s on the old path. This is that, rebuilt against the S4 payload.

⛔ **THE KEY IS THE WHOLE POINT, AND IT NEEDS THREE PARTS, NOT ONE.**

* ``graph_hash`` — regions plus the junction CSR. The payload already carries it, so the tally self-keys.
* ⭐ **a REACH digest.** logs that ``reach`` is consumed by calibration and is covered by
  **neither** ``partition_hash`` **nor** ``graph_hash`` — correctly, since neither the scan nor the
  accumulator reads it — and that the gap "becomes live the moment something caches a calibration
  output". A cache *loaded against an index* is that moment. The 2026-07-30 index rebuild moved ~38 % of
  human contiguous reaches while **both** hashes stayed byte-identical, so a reach-blind key would have
  loaded a stale cache against a moved index and verified clean.
* **the scan config**, because two scans of one BAM under different settings are different tallies.

⚠ Anything derivable from the index is **rebuilt on load, never stored** — `RegionArrays.from_index` is
0.11 s and `build_edge_flags_array` is 0.04 s against an 8.45 s index load that happens anyway.
Storing them is how a cache goes stale against the thing it describes.
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
    _stats, strand_model, _buffer, payload = scan_and_buffer(
        str(result.bam_path), result.index, scan
    )
    yield result.index, str(result.bam_path), scan, payload, strand_model
    scenario.cleanup()


def round_trip(scanned, tmp_path: Path) -> tuple:
    index, bam, scan, payload, strand_model = scanned
    cache_dir = tmp_path / "cache"
    write_scan_cache(
        cache_dir,
        payload=payload,
        strand_model=strand_model,
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
        """⭐ Field by field, and ONE LEVEL DOWN — a nested bank's arrays are the point, not its repr.

        ⛔ This compared nested banks with ``dataclasses.asdict(after) == dataclasses.asdict(before)``,
        which is a plain ``==`` over whatever they hold. That was adequate while every nested bank held only
        counters; ``DeferredFragments`` holds thirteen arrays, and the comparison first raised, then — had
        the arrays been scalars — would have compared two truncated string reprs and passed. Arrays are
        compared as arrays, with their dtype, at whatever depth they sit.
        """
        _cache_dir, cache = round_trip(scanned, tmp_path)
        original: AccumulatorPayload = scanned[3]

        def compare(before, after, name):
            if isinstance(before, np.ndarray):
                assert isinstance(after, np.ndarray), f"{name} came back as {type(after).__name__}"
                assert after.dtype == before.dtype, f"{name} dtype moved"
                assert after.shape == before.shape, f"{name} shape moved"
                assert np.array_equal(after, before), f"{name} content moved"
            elif dataclasses.is_dataclass(before):
                assert type(after) is type(before), (
                    f"{name} came back as {type(after).__name__}, not {type(before).__name__} — a bank "
                    f"rebuilt as a plain dict has lost its validation as well as its type"
                )
                for sub in dataclasses.fields(before):
                    compare(
                        getattr(before, sub.name), getattr(after, sub.name), f"{name}.{sub.name}"
                    )
            else:
                assert after == before, name

        for field in dataclasses.fields(AccumulatorPayload):
            compare(getattr(original, field.name), getattr(cache.payload, field.name), field.name)

    def test_THE_SIDE_BUFFER_SURVIVES_THE_ROUND_TRIP_AND_IS_NOT_EMPTY(self, scanned, tmp_path):
        """⭐ **S2's gate.** Byte-identity over an empty bank is free, so the fixture must defer something.

        ⛔ And the bank must come back as a ``DeferredFragments``, not as whatever JSON happened to hold.
        The write path splits a nested dataclass by the type of each sub-field precisely because
        ``json.dumps(..., default=str)`` would stringify an ndarray to a TRUNCATED repr — silently, and the
        second pass would then drain coordinates parsed out of text.
        """
        from rigel.scan_payload import DeferredFragments

        _cache_dir, cache = round_trip(scanned, tmp_path)
        original: AccumulatorPayload = scanned[3]
        assert original.deferred.n_fragments > 0, (
            "nothing was deferred, so this round trip is asserted over an empty bank — which survives any "
            "serialisation at all. The fixture must produce an undetermined gap."
        )
        assert isinstance(cache.payload.deferred, DeferredFragments)
        for field in dataclasses.fields(DeferredFragments):
            before = getattr(original.deferred, field.name)
            after = getattr(cache.payload.deferred, field.name)
            assert after.dtype == np.int64 == before.dtype, f"deferred.{field.name} dtype moved"
            assert np.array_equal(after, before), f"deferred.{field.name} content moved"
        assert cache.payload.gap_resolution == original.gap_resolution

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

    def test_THE_FL_HISTOGRAMS_ARE_THE_PAYLOAD_AND_ARE_NOT_CACHED_SEPARATELY(
        self, scanned, tmp_path
    ):
        """⭐ TRAPS: pure-and-length-censored: there is nothing left to cache beside the payload.

        The cache used to carry a second FL block — the scanner's own histogram, stored so it could
        serve as the empirical-Bayes anchor, plus a ``fl_rna_counts`` field (**TRAPS: no-prior-means-haldane**) that was written,
        read back, and consumed by nothing. Both are gone. Every fragment-length histogram is a FIELD
        of the payload now, so caching the payload caches them, in one frame, by construction.

        ⚠ ``build_fl_models`` still stays the single source of truth for the derived pmfs, which are
        still not cached — freezing its output would mean a change to the FL model silently does not
        apply to cached scans.
        """
        import dataclasses

        from rigel.calibration.fl import build_fl_models

        _cache_dir, cache = round_trip(scanned, tmp_path)

        fields = {f.name for f in dataclasses.fields(cache)}
        assert not {"fl_global_counts", "fl_rna_counts", "fl_max_size"} & fields, (
            f"a separate fragment-length block survives on the cache: {sorted(fields)}"
        )

        # And the models rebuild from the cached payload alone, identically to the live scan's.
        live = build_fl_models(scanned[3])
        cached = build_fl_models(cache.payload)
        assert np.array_equal(cached.global_counts, live.global_counts)
        assert np.array_equal(cached.rna_counts, live.rna_counts)
        assert np.array_equal(cached.gdna_counts, live.gdna_counts)


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
        agreement — that is how a stale cache verifies clean."""
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
        trimmed = {k: v for k, v in original.items() if k != "sj_inv_length_sum"}
        try:
            scan_payload.AccumulatorPayload.__dataclass_fields__ = trimmed
            assert payload_schema_digest() != before, "dropping a field must move the digest"
        finally:
            scan_payload.AccumulatorPayload.__dataclass_fields__ = original
        assert payload_schema_digest() == before

    def test_the_deposit_digest_is_STABLE_across_calls_and_across_PROCESSES(self):
        """⭐ A key that wobbles refuses every cache, which is as useless as one that never moves.

        Across processes as well as calls: every channel is an integer and integer addition is
        associative, so this is deterministic by construction rather than by luck.
        """
        import subprocess
        import sys

        from rigel.scan_cache import deposit_digest

        first = deposit_digest()
        assert first == deposit_digest()
        other = subprocess.run(
            [sys.executable, "-c", "from rigel.scan_cache import deposit_digest;print(deposit_digest())"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        assert other == first, f"digest differs across processes: {other} != {first}"

    def test_the_deposit_digest_ENTERS_the_schema_digest(self):
        """⛔⛔ **THE WIRING, and without it the digest is a number nobody consults.**

        ``payload_schema_digest`` hashes field NAMES and column counts, and a deposit-RULE change moves
        neither — so a cache written under the old rule was accepted by the key and silently served OLD
        VALUES to NEW CODE. That is ``TRAPS: a-hash-that-misses-its-artifact`` in the key written to
        prevent it, for the FOURTH time. This asserts the deposit digest is actually folded in.
        """
        from rigel import scan_cache

        before = scan_cache.payload_schema_digest()
        original = scan_cache.deposit_digest
        try:
            scan_cache.deposit_digest = lambda: "not-the-real-digest"
            assert scan_cache.payload_schema_digest() != before, (
                "the deposit digest does not reach the cache key, so a deposit-rule change would "
                "leave every stale cache accepted"
            )
        finally:
            scan_cache.deposit_digest = original
        assert scan_cache.payload_schema_digest() == before

    def test_a_changed_DEPOSIT_RULE_moves_the_digest(self):
        """⭐⭐⭐ **THE CLAIM ITSELF: the digest is a function of BEHAVIOUR, not of names.**

        Perturbs the deposit rule in the executable specification — which the C++ is held byte-identical
        to — and asserts the digest computed over it moves. ⛔ Every field name, dtype and shape is
        untouched by the perturbation, so ``payload_schema_digest``'s name/column half cannot see it;
        only a behavioural hash can.
        """
        from tests.native._accumulator_reference import Accumulator, Partition

        from tests.native._digest_fixture import reference_deposit_digest

        before = reference_deposit_digest(Accumulator, Partition)
        original = Accumulator.deposit
        try:
            # A rule change that renames nothing: halve every conserved-mass deposit.
            def perturbed(self, *args, **kwargs):
                outcome = original(self, *args, **kwargs)
                self.tally.sj_mass //= 2
                return outcome

            Accumulator.deposit = perturbed
            assert reference_deposit_digest(Accumulator, Partition) != before, (
                "a changed deposit rule left the digest where it was — it is hashing names, not "
                "behaviour, and every stale cache would still be accepted"
            )
        finally:
            Accumulator.deposit = original
        assert reference_deposit_digest(Accumulator, Partition) == before

    def test_the_deposit_digest_AGREES_between_the_NATIVE_and_the_SPECIFICATION(self):
        """⛔ The key certifies what the PRODUCTION scanner deposited, so it is computed from the C++.

        If the two ever drifted, the key would be certifying an artifact nothing writes. This is the
        gate that makes reading the native side safe.
        """
        from tests.native._accumulator_reference import Accumulator, Partition

        from rigel.scan_cache import deposit_digest

        from tests.native._digest_fixture import reference_deposit_digest

        assert reference_deposit_digest(Accumulator, Partition) == deposit_digest(), (
            "the specification and the C++ disagree on the digest fixture, so the cache key certifies "
            "an artifact the scanner does not produce"
        )

    def test_the_schema_digest_moves_when_a_BANK_CHANGES_SHAPE(self):
        """⛔⛔ **The gap this key had until 2026-08-08, and it is the one it exists to close.**

        A bank moving from two genome-strand columns to one keeps its name, its dtype and its axis — so a
        name-only digest does not move. And nothing downstream catches it: ``_bank`` validates the C++
        dict in ``from_scan_result``, while ``_payload_from_parts`` puts the ``.npz`` arrays STRAIGHT
        into the payload with no shape check. A stale cache would have been accepted by the key and then
        failed with a shape error pointing nowhere near its cause — ``TRAPS: a-hash-that-misses-its-artifact``.

        ⚠ Measured: collapsing five length moments to ``[n]`` left the digest at ``a025995ea3ce7d4f``,
        byte for byte.
        """
        from rigel.scan_cache import payload_schema_digest

        before = payload_schema_digest()
        original = scan_payload.SINGLE_COLUMN_AXES
        moved = tuple(row for row in original if row[0] != "sj_inv_length_sum")
        try:
            scan_payload.SINGLE_COLUMN_AXES = moved
            assert payload_schema_digest() != before, (
                "a bank changing its column count must move the digest — the field NAMES are identical, "
                "so nothing else can notice"
            )
        finally:
            scan_payload.SINGLE_COLUMN_AXES = original
        assert payload_schema_digest() == before

    def test_the_schema_digest_moves_when_a_NESTED_BANK_FIELD_moves(self):
        """⭐ **The defect this key was written to prevent, and for a long time did not.**

        The digest hashed ``AccumulatorPayload``'s top-level field names ALONE, so a change *inside* a
        nested bank was invisible to it: renaming a ``ScanQC`` field let a stale cache be **accepted by the
        key** and then fail deep in ``_payload_from_parts`` with a bare ``TypeError`` — precisely the
        failure mode the docstring above says the key exists to prevent. Two statements about one contract,
        disagreeing (the same shape as ``check_scan_config``'s).

        ⛔ S1 made the stakes much higher rather than lower: ``DeferredFragments`` puts **thirteen** array
        names inside one payload field and every one of them is an ``.npz`` key, so an unrecursed digest
        would wave through a cache missing an entire bank.

        ⚠ Tested one level down on each of the three nested banks, because "it recurses" is not the claim —
        the claim is that each specific bank is covered, and a recursion that skipped one would still
        recurse.
        """
        from rigel.scan_cache import _payload_field_types, _schema_names, payload_schema_digest

        before = payload_schema_digest()
        names = set(_schema_names())
        banks = {
            name: cls
            for name, cls in _payload_field_types().items()
            if dataclasses.is_dataclass(cls)
        }
        assert set(banks) == {"qc", "deferred", "gap_resolution"}, (
            f"the payload's nested banks are {sorted(banks)}; a new one must join this gate deliberately "
            f"rather than by a recursion that happens to reach it"
        )
        for bank, nested in banks.items():
            for sub in dataclasses.fields(nested):
                assert f"{bank}__{sub.name}" in names, (
                    f"{bank}.{sub.name} is not in the digest's name list, so a cache written before it "
                    f"existed would be accepted by the key"
                )
            original = nested.__dataclass_fields__
            dropped = next(iter(original))
            try:
                nested.__dataclass_fields__ = {k: v for k, v in original.items() if k != dropped}
                assert payload_schema_digest() != before, (
                    f"dropping {bank}.{dropped} did not move the digest — a change inside a nested bank "
                    f"is invisible to the key"
                )
            finally:
                nested.__dataclass_fields__ = original
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
    def test_the_cache_holds_no_region_arrays_or_edge_flags(self, scanned, tmp_path):
        """They are 0.15 s to rebuild and a stale copy of them is a silent wrong answer."""
        cache_dir, _cache = round_trip(scanned, tmp_path)
        stored = {p.name for p in cache_dir.iterdir()}
        assert not any("region" in name or "boundary" in name for name in stored), stored

    def test_the_index_derived_inputs_are_rebuilt_not_loaded(self, scanned, tmp_path):
        index = scanned[0]
        from rigel.calibration.region_arrays import RegionArrays
        from rigel.calibration.splice_graph import build_edge_flags_array

        inputs = index_derived_inputs(index)
        assert inputs["region_arrays"].n_regions == RegionArrays.from_index(index).n_regions
        assert np.array_equal(inputs["edge_flags"], build_edge_flags_array(index))

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
        # ⚠ EVERY argument without a usable default, not a hand-kept list. `junctions` was added to
        # `calibrate` at S5.f and this helper was not updated — a hand-kept `required` set cannot catch
        # that, because the omission is invisible until something calls it.
        required = {
            name
            for name, param in inspect.signature(calibrate).parameters.items()
            if param.default is inspect.Parameter.empty and name != "config"
        } | {"junctions"}
        assert required <= set(inputs), f"missing: {required - set(inputs)}"
        assert inputs["gdna_fl_pmf"].sum() == pytest.approx(1.0, abs=1e-9)

    def test_calibration_inputs_ACTUALLY_DRIVE_calibrate(self, scanned, tmp_path):
        """⭐ The signature check above cannot see a MIS-SIZED argument, only a missing name.

        `junctions` must address the same graph the payload was scanned on; an axis of the wrong length
        places every splice on the wrong line. Calling `calibrate` for real is the only check that
        covers it, and it is cheap on this fixture.
        """
        from rigel.calibration.calibrate import calibrate
        from rigel.config import CalibrationConfig

        _cache_dir, cache = round_trip(scanned, tmp_path)
        result = calibrate(**calibration_inputs(cache, scanned[0]), config=CalibrationConfig())
        assert result.n_regions > 0
        assert result.n_junctions == cache.payload.n_sj


def test_population_priors_can_be_extracted_from_a_cached_scan(scanned, tmp_path):
    """✅ Unblocked by S5.f — `calibrate()` runs, so the population-prior seed path is live.

    ⚠ It was a **strict** xfail naming S5 as the blocker. Strict is what made it honest: the moment it
    started working it would have failed loudly rather than sitting green and forgotten. ⛔ It was in
    fact still xfailing at S5.f for a DIFFERENT reason than the one recorded — `calibration_inputs` had
    not been given `junctions` — which is exactly why a stale xfail reason is worth nothing.
    """
    from rigel.calibration.calibrate import calibrate
    from rigel.config import CalibrationConfig

    _cache_dir, cache = round_trip(scanned, tmp_path)
    debug: dict = {}
    calibrate(**calibration_inputs(cache, scanned[0]), config=CalibrationConfig(), _debug=debug)
    assert debug["calibration_priors"] is not None


# ══ P3 · the drain and the cache ════════════════════════════════════════════════════════════════════


def test_the_schema_digest_sees_fields_nested_TWO_levels_deep():
    """⛔ The digest recursed exactly ONE level until 2026-08-02, and `DrainQC` nests a `GapCensus` inside
    itself — so `drain__census_before__*` was invisible to the key. That is X8's defect one level down.

    ⚠ **And `drain` is `DrainQC | None`**, for which `dataclasses.is_dataclass` is False, so a plain check
    dropped the whole bank rather than just its nesting. Both halves are gated here.
    """
    from rigel.scan_cache import _schema_names

    names = set(_schema_names())
    assert "drain" in names
    assert "drain__offered" in names, (
        "the drain bank's own fields are missing from the schema key — `DrainQC | None` is a union, and "
        "`is_dataclass` is False for it."
    )
    assert "drain__census_before__gap_resolved_spliced" in names, (
        "a field nested TWO levels deep is invisible to the schema key; the digest is what refuses a "
        "stale cache instead of failing deep in the loader."
    )


def test_the_cache_REFUSES_a_drained_payload(scanned, tmp_path):
    """⛔ The cache holds a SCAN. Caching a drained payload would bake one draw into it and destroy the
    property the whole second-pass design rests on — that one scan can be drained repeatedly at different
    seeds without re-reading the BAM (and P5/P6 both need it).

    ⚠ It would also serialise `DrainQC.census_before` through `json.dumps(default=str)` as a stringified
    repr, silently — X9's defect one level down.
    """
    import dataclasses

    from rigel.scan_cache import write_scan_cache
    from rigel.scan_payload import DrainQC, GapCensus

    index, bam, scan_config, payload, strand_model = scanned
    drained = dataclasses.replace(
        payload,
        drain=DrainQC(
            offered=0,
            deposited=0,
            dropped_too_long=0,
            dropped_empty=0,
            dropped_strand_undefined=0,
            chose_genomic=0,
            chose_spliced=0,
            census_before=GapCensus.zeros(),
        ),
    )
    with pytest.raises(ValueError, match="DRAINED"):
        write_scan_cache(
            tmp_path / "refused",
            payload=drained,
            strand_model=strand_model,
            index=index,
            bam=bam,
            scan_config=scan_config,
        )


def test_an_UNDRAINED_payload_still_round_trips_with_the_new_field(scanned, tmp_path):
    """⚠ The other side of the refusal: adding the field must not break the ordinary path, and `drain` must
    come back as `None` rather than as the string ``"None"``."""
    from rigel.scan_cache import read_scan_cache, write_scan_cache

    index, bam, scan_config, payload, strand_model = scanned
    assert payload.drain is None
    write_scan_cache(
        tmp_path / "ok",
        payload=payload,
        strand_model=strand_model,
        index=index,
        bam=bam,
        scan_config=scan_config,
    )
    loaded = read_scan_cache(tmp_path / "ok", index)
    assert loaded.payload.drain is None, (
        f"drain came back as {loaded.payload.drain!r}; a JSON round trip must not turn None into a string"
    )
    assert loaded.payload.deferred.n_fragments == payload.deferred.n_fragments
