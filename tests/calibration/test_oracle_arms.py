"""Tests for ``scripts/design/_oracle_arms.py``: the oracle cache and the per-object scorer ``calibration_vs_oracle.py``
reads.

Every number the scorer prints is a difference between two per-object arrays, so every way of getting it wrong yields a
plausible answer; each test breaks the thing it guards and asserts the guard notices. The toy serves the cache and
sum-to-full tests; the scoring tests run on synthetic arrays. A single-reference toy cannot judge the deposit path,
since it hides ref-id mismatches; ``tests/native/`` covers that path.
"""

from __future__ import annotations

import dataclasses
from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer
import _oracle_arms
from _prior_toy import build_toy


# ── the toy ──────────────────────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """gDNA + mature + nascent, with spliced reads. The staggered isoform boundary matters: ``boundary_spliced`` (a
    molecule that crossed a contiguous boundary having spliced elsewhere) can only be deposited where a region bound
    falls inside another transcript's exon, so a single-isoform gene would leave the bank zero and gate 1's perturbation
    would remove nothing.
    """
    return build_toy(tmp_path_factory)


# ── Gate 0: the oracle cache hits, is refused when stale, and is still validated ──────────────────


def test_the_oracle_cache_hits_and_the_cached_path_is_validated(toy, tmp_path, monkeypatch):
    """The cache lets a panel be re-measured without re-splitting every BAM, provided:

    1. a warm run does not touch the BAM (``split_by_origin`` is made to raise);
    2. it reproduces the cold build's arrays exactly;
    3. it still runs sum-to-full: a corrupted cached partition on disk must abort the warm load.
    """
    import _oracle

    from rigel.config import PipelineConfig as PC

    bam, index, cfg = str(toy.bam_path), toy.index, PC()
    cache = tmp_path / "orc_cache"
    *_, full = scan_and_buffer(bam, index, cfg.scan)

    cold = _oracle_arms.load_or_build_oracle(bam, index, cfg, tmp_path / "w1", "t", full, cache)

    def _boom(*a, **k):
        raise AssertionError("the warm path re-split the BAM; the cache did not hit")

    monkeypatch.setattr(_oracle_arms, "split_by_origin", _boom)
    warm = _oracle_arms.load_or_build_oracle(bam, index, cfg, tmp_path / "w2", "t", full, cache)
    for origin in _oracle.ORIGINS:
        np.testing.assert_array_equal(
            np.asarray(warm.parts[origin].region_contained_count),
            np.asarray(cold.parts[origin].region_contained_count),
        )

    # a forced build re-splits even though the cache would hit: after a change to what the scan offers a stale partition
    # still loads, and only a rebuild makes it describe the scan again
    calls = {"n": 0}
    real = _oracle.split_by_origin

    def counting(*a, **k):
        calls["n"] += 1
        return real(*a, **k)

    monkeypatch.setattr(_oracle_arms, "split_by_origin", counting)
    forced = _oracle_arms.load_or_build_oracle(
        bam, index, cfg, tmp_path / "wf", "t", full, cache, force=True
    )
    assert calls["n"] == 1, "a forced build reused the cache"
    for origin in _oracle.ORIGINS:
        np.testing.assert_array_equal(
            np.asarray(forced.parts[origin].region_contained_count),
            np.asarray(cold.parts[origin].region_contained_count),
        )
    monkeypatch.setattr(_oracle_arms, "split_by_origin", _boom)

    # the cached path still runs sum-to-full
    npz = cache / "t" / "gdna" / "payload.npz"
    data = {k: v for k, v in np.load(npz).items()}
    data["region_contained_count"] = data["region_contained_count"].copy()
    data["region_contained_count"][0, 0] += 1
    np.savez_compressed(npz, **data)
    with pytest.raises(AssertionError, match="oracle INVALID"):
        _oracle_arms.load_or_build_oracle(bam, index, cfg, tmp_path / "w3", "t", full, cache)


def test_a_cache_that_does_not_describe_this_scan_is_rebuilt_not_reused(toy, tmp_path, monkeypatch):
    """A cache keyed to a different scan configuration is a different tally: it must be rebuilt, never reused, and a
    miss is normal rather than an error. The cache is populated under the default scan config, then requested under a
    changed one, which must re-split the BAM.

    ``full`` is re-scanned under the changed config too: two scan configs are two tallies, so a default-config ``full``
    would fail sum-to-full. The cache key refuses a stale partition, and the identity independently refuses a ``full``
    that does not match the partitions.
    """
    from rigel.config import PipelineConfig as PC

    bam, index, cfg = str(toy.bam_path), toy.index, PC()
    cache = tmp_path / "orc_cache"
    *_, full = scan_and_buffer(bam, index, cfg.scan)
    _oracle_arms.load_or_build_oracle(bam, index, cfg, tmp_path / "w1", "t", full, cache)

    calls = {"n": 0}
    real = _oracle_arms.split_by_origin

    def counting(*a, **k):
        calls["n"] += 1
        return real(*a, **k)

    monkeypatch.setattr(_oracle_arms, "split_by_origin", counting)
    changed = dataclasses.replace(
        cfg, scan=dataclasses.replace(cfg.scan, max_frag_length=cfg.scan.max_frag_length // 2)
    )
    *_, changed_full = scan_and_buffer(bam, index, changed.scan)
    _oracle_arms.load_or_build_oracle(
        bam, index, changed, tmp_path / "w2", "t", changed_full, cache
    )
    assert calls["n"] == 1, "a cache from a different scan config was reused instead of rebuilt"

    # and the rebuilt cache now serves the changed config without re-splitting.
    _oracle_arms.load_or_build_oracle(
        bam, index, changed, tmp_path / "w3", "t", changed_full, cache
    )
    assert calls["n"] == 1, "the rebuilt cache did not hit on the next call"


# ── Gate 1: T's per-axis totals are the full payload's totals ─────────────────────────────────────


def test_the_truth_totals_equal_the_full_payload_per_axis(toy, tmp_path):
    """``check_same_basis`` must hold between T and the payload it partitions, per axis, never pooled: ``n_regions`` and
    ``n_boundaries`` differ by only ``n_refs``, so equal and opposite errors on the two axes could cancel. Dropping the
    spliced term from ``count_rna_boundary`` must be refused: ``chain_boundary_deconv`` builds ``rna = (1−f_g)·unspliced
    + spliced``, so such a T would be on a different basis from every P.
    """
    from _oracle import OracleTruth

    oracle = OracleTruth.from_bam(str(toy.bam_path), toy.index, PipelineConfig(), tmp_path, "t")
    ra = RegionArrays.from_frame(toy.index.regions_df, toy.index.ref_name_to_id)
    full = CalibrationSubstrate.from_payload(oracle.full, ra)
    truth = oracle.override_masses(ra)
    _oracle_arms.check_same_basis("T", SimpleNamespace(**truth), full)  # holds

    spliced = np.asarray(full.boundary_spliced.count, np.float64).sum(axis=1)
    assert spliced.sum() > 0, "the toy must EXERCISE the spliced bank or this perturbation is inert"
    broken = {**truth, "count_rna_boundary": truth["count_rna_boundary"] - spliced}
    with pytest.raises(ValueError, match="boundary"):
        _oracle_arms.check_same_basis("T", SimpleNamespace(**broken), full)


# ── Gate 2: an arm and its truth are on the same basis, per object ────────────────────────────────


def test_an_arm_on_a_different_basis_is_refused():
    """Every arm's per-object total must equal T's, which makes ``f_g`` comparable: fractions over different
    denominators subtract to nothing meaningful. Scoring against another axis' truth (a different object count) and
    scaling the arm's masses (a different denominator) are both refused.
    """
    g, r = np.array([4.0, 1.0, 3.0]), np.array([6.0, 2.0, 0.0])
    tg, tr = np.array([5.0, 0.0, 3.0]), np.array([5.0, 3.0, 0.0])
    _oracle_arms.score_axis(g, r, tg, tr)  # the same per-object totals: holds
    with pytest.raises(ValueError, match="different axes"):
        _oracle_arms.score_axis(g, r, tg[:-1], tr[:-1])
    with pytest.raises(ValueError, match="DIFFERENT BASES"):
        _oracle_arms.score_axis(g * 2.0, r * 2.0, tg, tr)


# ── Gate 3: no data is absent, never f_g = 0 ──────────────────────────────────────────────────────


def test_an_object_with_no_mass_is_absent_not_a_confident_zero():
    """No data must be inert, neither "100 % gDNA" nor "0 % gDNA". Most regions in a real index carry no fragments, so a
    scorer turning 0/0 into a number would report a good answer for most of the genome. The mass-weighted mean gives a
    zero-mass object zero weight and so cannot see this; the test is on the count of scored objects, where a floored 0
    would show.
    """
    g = np.array([4.0, 0.0, 0.0, 1.0])
    r = np.array([6.0, 0.0, 0.0, 0.0])
    frac, total = _oracle_arms.object_fractions(g, r)
    assert np.isnan(frac[1]) and np.isnan(frac[2]), "0/0 must be NaN, not 0.0"
    np.testing.assert_allclose(frac[[0, 3]], [0.4, 1.0])

    score = _oracle_arms.score_axis(
        g, r, np.array([2.0, 0.0, 0.0, 1.0]), np.array([8.0, 0.0, 0.0, 0.0])
    )
    assert score.n_scored == 2, "empty objects were counted as perfectly scored objects"
    assert score.mass == 11.0
    np.testing.assert_allclose(score.abs_err, 2.0)
    np.testing.assert_allclose(score.mwae, 2.0 / 11.0)


# ── Gate 4: the directional split is reported and does not cancel ─────────────────────────────────


def test_the_directional_split_is_reported_and_the_net_is_their_difference():
    """A library-level net hides a large under-call beside a large over-call, so the two directions are separate fields
    and their relationship is an identity. An arm whose errors all point one way must have one direction exactly zero,
    or the split would be measuring noise.
    """
    s = _oracle_arms.score_axis(
        np.array([4.0, 1.0, 3.0]),
        np.array([6.0, 2.0, 0.0]),
        np.array([5.0, 0.0, 3.0]),
        np.array([5.0, 3.0, 0.0]),
    )
    assert s.over_call > 0.0 and s.under_call > 0.0, "the fixture must err both ways"
    np.testing.assert_allclose(s.over_call - s.under_call, s.net_err, atol=1e-12)
    np.testing.assert_allclose(s.over_call + s.under_call, s.abs_err, rtol=1e-12)

    one_way = _oracle_arms.score_axis(
        np.array([5.0, 6.0]),
        np.array([5.0, 4.0]),
        np.array([3.0, 2.0]),
        np.array([7.0, 8.0]),
    )
    assert one_way.under_call == 0.0 and one_way.over_call > 0.0


# ── the certified table is keyed on the scan it was derived from ──────────────────────────────────


def _slot_truth(oracle_dir, **extra):
    """A two-slot certified table: one region, one boundary."""
    from rigel.calibration.region_chain import BOUNDARY, REGION

    oracle_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        oracle_dir / "slot_truth.npz",
        kind=np.array([REGION, BOUNDARY]),
        obj=np.array([0, 0]),
        n_gdna=np.array([5.0, 2.0]),
        **extra,
    )


def _scan(**key):
    """A scan cache's manifest, the provenance `read_scan_cache` returns and `certified_gdna` is keyed on."""
    manifest = {
        "graph_hash": "a" * 16,
        "reach_digest": "b" * 16,
        "payload_schema_digest": "c" * 16,
        "transcript_model_digest": "d" * 16,
        "scan_config": {"sj_strand_tag": "auto"},
        "bam": {"blake2b": "e" * 16},
    }
    return {**manifest, **key}


def test_a_certified_table_reads_against_the_scan_it_was_derived_from(tmp_path):
    from _scan_cache import scan_cache_key

    scan = _scan()
    _slot_truth(tmp_path / "c", scan_key=scan_cache_key(scan))
    truth = _oracle_arms.certified_gdna(tmp_path / "c", scan)
    assert truth["region"].tolist() == [5.0] and truth["boundary"].tolist() == [2.0]


def test_a_certified_table_from_another_scan_is_refused(tmp_path):
    """A re-scan (a new BAM, an index rebuild, a deposit-rule change) leaves `slot_truth.npz` behind until
    `calibration_oracle.py` re-certifies it; read against the new scan it would score the new estimate against the
    old truth."""
    from _scan_cache import ScanCacheKeyError, scan_cache_key

    _slot_truth(tmp_path / "c", scan_key=scan_cache_key(_scan()))
    with pytest.raises(ScanCacheKeyError, match="slot_truth"):
        _oracle_arms.certified_gdna(tmp_path / "c", _scan(graph_hash="f" * 16))
    with pytest.raises(ScanCacheKeyError, match="slot_truth"):
        _oracle_arms.certified_gdna(tmp_path / "c", _scan(bam={"blake2b": "0" * 16}))


def test_a_certified_table_that_records_no_scan_is_refused(tmp_path):
    from _scan_cache import ScanCacheKeyError

    _slot_truth(tmp_path / "c")
    with pytest.raises(ScanCacheKeyError, match="slot_truth"):
        _oracle_arms.certified_gdna(tmp_path / "c", _scan())
