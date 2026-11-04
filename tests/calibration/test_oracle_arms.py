"""Tests for ``scripts/design/_oracle_arms.py``: the oracle cache and the per-object scorer ``calibration_vs_oracle.py``
reads.

Every number the scorer prints is a difference between two per-object arrays, so every way of getting it wrong yields a
plausible answer; each test breaks the thing it guards and asserts the guard notices. The toy serves the cache and
sum-to-full tests; the scoring tests run on synthetic arrays. A single-reference toy cannot judge the deposit path,
since it hides ref-id mismatches; ``tests/native/`` covers that path.
"""

from __future__ import annotations

import dataclasses
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from _prior_toy import build_toy


def _load_sibling(name: str):
    """Import a ``scripts/design/`` instrument by path.

    ``scripts/`` is not a package: it is a toolkit of instruments, not a library with a public API, and loading by path
    keeps the dependency one-directional. The module goes into ``sys.modules`` before it runs, because ``@dataclass``
    resolves its class's ``__module__`` through that table and an unregistered module fails at class definition. The
    scripts' directory goes on ``sys.path`` too, as running a script puts it there: the instruments import their shared
    helper ``_shared`` by name.
    """
    import sys

    key = name[:-3]
    if key not in _MODULES:
        path = Path(__file__).resolve().parents[2] / "scripts" / "design" / name
        if str(path.parent) not in sys.path:
            sys.path.insert(0, str(path.parent))
        spec = importlib.util.spec_from_file_location(key, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        _MODULES[key] = module
        spec.loader.exec_module(module)
    return _MODULES[key]


_MODULES: dict = {}
P0 = _load_sibling("_oracle_arms.py")


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


def test_the_oracle_cache_hits_and_the_cached_path_is_STILL_VALIDATED(toy, tmp_path, monkeypatch):
    """The cache lets a panel be re-measured without re-splitting every BAM, provided:

    1. a warm run does not touch the BAM (``_split_bam`` is made to raise);
    2. it reproduces the cold build's arrays exactly;
    3. it still runs sum-to-full: a corrupted cached partition on disk must abort the warm load.
    """
    import _oracle

    import _oracle_arms as mod

    from rigel.config import PipelineConfig as PC

    bam, index, cfg = str(toy.bam_path), toy.index, PC()
    cache = tmp_path / "orc_cache"
    full = _oracle._scan_payload(bam, index, cfg)

    cold = mod.load_or_build_oracle(bam, index, cfg, tmp_path / "w1", "t", full, cache)

    def _boom(*a, **k):
        raise AssertionError("the warm path re-split the BAM; the cache did not hit")

    monkeypatch.setattr(mod, "_split_bam", _boom)
    warm = mod.load_or_build_oracle(bam, index, cfg, tmp_path / "w2", "t", full, cache)
    for origin in _oracle.ORIGINS:
        np.testing.assert_array_equal(
            np.asarray(warm.parts[origin].region_contained_count),
            np.asarray(cold.parts[origin].region_contained_count),
        )

    # a forced build re-splits even though the cache would hit: after a change to what the scan offers a stale partition
    # still loads, and only a rebuild makes it describe the scan again
    calls = {"n": 0}
    real = _oracle._split_bam

    def counting(*a, **k):
        calls["n"] += 1
        return real(*a, **k)

    monkeypatch.setattr(mod, "_split_bam", counting)
    forced = mod.load_or_build_oracle(
        bam, index, cfg, tmp_path / "wf", "t", full, cache, force=True
    )
    assert calls["n"] == 1, "a forced build reused the cache"
    for origin in _oracle.ORIGINS:
        np.testing.assert_array_equal(
            np.asarray(forced.parts[origin].region_contained_count),
            np.asarray(cold.parts[origin].region_contained_count),
        )
    monkeypatch.setattr(mod, "_split_bam", _boom)

    # the cached path still runs sum-to-full
    npz = cache / "t" / "gdna" / "payload.npz"
    data = {k: v for k, v in np.load(npz).items()}
    data["region_contained_count"] = data["region_contained_count"].copy()
    data["region_contained_count"][0, 0] += 1
    np.savez_compressed(npz, **data)
    with pytest.raises(AssertionError, match="oracle INVALID"):
        mod.load_or_build_oracle(bam, index, cfg, tmp_path / "w3", "t", full, cache)


def test_a_cache_that_does_not_describe_this_SCAN_is_rebuilt_not_reused(toy, tmp_path, monkeypatch):
    """A cache keyed to a different scan configuration is a different tally: it must be rebuilt, never reused, and a
    miss is normal rather than an error. The cache is populated under the default scan config, then requested under a
    changed one, which must re-split the BAM.

    ``full`` is re-scanned under the changed config too: two scan configs are two tallies, so a default-config ``full``
    would fail sum-to-full. The cache key refuses a stale partition, and the identity independently refuses a ``full``
    that does not match the partitions.
    """
    import _oracle

    import _oracle_arms as mod

    from rigel.config import PipelineConfig as PC

    bam, index, cfg = str(toy.bam_path), toy.index, PC()
    cache = tmp_path / "orc_cache"
    full = _oracle._scan_payload(bam, index, cfg)
    mod.load_or_build_oracle(bam, index, cfg, tmp_path / "w1", "t", full, cache)

    calls = {"n": 0}
    real = mod._split_bam

    def counting(*a, **k):
        calls["n"] += 1
        return real(*a, **k)

    monkeypatch.setattr(mod, "_split_bam", counting)
    changed = dataclasses.replace(
        cfg, scan=dataclasses.replace(cfg.scan, max_frag_length=cfg.scan.max_frag_length // 2)
    )
    changed_full = _oracle._scan_payload(bam, index, changed)
    mod.load_or_build_oracle(bam, index, changed, tmp_path / "w2", "t", changed_full, cache)
    assert calls["n"] == 1, "a cache from a different scan config was reused instead of rebuilt"

    # and the rebuilt cache now serves the changed config without re-splitting.
    mod.load_or_build_oracle(bam, index, changed, tmp_path / "w3", "t", changed_full, cache)
    assert calls["n"] == 1, "the rebuilt cache did not hit on the next call"


# ── Gate 1: T's per-axis totals are the full payload's totals ─────────────────────────────────────


def test_T_totals_equal_the_full_payload_PER_AXIS(toy, tmp_path):
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
    P0.check_same_basis("T", SimpleNamespace(**truth), full)  # holds

    spliced = np.asarray(full.boundary_spliced.count, np.float64).sum(axis=1)
    assert spliced.sum() > 0, "the toy must EXERCISE the spliced bank or this perturbation is inert"
    broken = {**truth, "count_rna_boundary": truth["count_rna_boundary"] - spliced}
    with pytest.raises(ValueError, match="boundary"):
        P0.check_same_basis("T", SimpleNamespace(**broken), full)


# ── Gate 2: an arm and its truth are on the same basis, per object ────────────────────────────────


def test_an_arm_on_a_DIFFERENT_BASIS_is_refused():
    """Every arm's per-object total must equal T's, which makes ``f_g`` comparable: fractions over different
    denominators subtract to nothing meaningful. Scoring against another axis' truth (a different object count) and
    scaling the arm's masses (a different denominator) are both refused.
    """
    g, r = np.array([4.0, 1.0, 3.0]), np.array([6.0, 2.0, 0.0])
    tg, tr = np.array([5.0, 0.0, 3.0]), np.array([5.0, 3.0, 0.0])
    P0.score_axis(g, r, tg, tr)  # the same per-object totals: holds
    with pytest.raises(ValueError, match="different axes"):
        P0.score_axis(g, r, tg[:-1], tr[:-1])
    with pytest.raises(ValueError, match="DIFFERENT BASES"):
        P0.score_axis(g * 2.0, r * 2.0, tg, tr)


# ── Gate 3: no data is absent, never f_g = 0 ──────────────────────────────────────────────────────


def test_an_object_with_no_mass_is_ABSENT_not_a_confident_zero():
    """No data must be inert, neither "100 % gDNA" nor "0 % gDNA". Most regions in a real index carry no fragments, so a
    scorer turning 0/0 into a number would report a good answer for most of the genome. The mass-weighted mean gives a
    zero-mass object zero weight and so cannot see this; the test is on the count of scored objects, where a floored 0
    would show.
    """
    g = np.array([4.0, 0.0, 0.0, 1.0])
    r = np.array([6.0, 0.0, 0.0, 0.0])
    frac, total = P0.object_fractions(g, r)
    assert np.isnan(frac[1]) and np.isnan(frac[2]), "0/0 must be NaN, not 0.0"
    np.testing.assert_allclose(frac[[0, 3]], [0.4, 1.0])

    score = P0.score_axis(g, r, np.array([2.0, 0.0, 0.0, 1.0]), np.array([8.0, 0.0, 0.0, 0.0]))
    assert score.n_scored == 2, "empty objects were counted as perfectly scored objects"
    assert score.mass == 11.0
    np.testing.assert_allclose(score.abs_err, 2.0)
    np.testing.assert_allclose(score.mwae, 2.0 / 11.0)


# ── Gate 4: the directional split is reported and does not cancel ─────────────────────────────────


def test_the_directional_split_is_reported_and_the_net_is_their_DIFFERENCE():
    """A library-level net hides a large under-call beside a large over-call, so the two directions are separate fields
    and their relationship is an identity. An arm whose errors all point one way must have one direction exactly zero,
    or the split would be measuring noise.
    """
    s = P0.score_axis(
        np.array([4.0, 1.0, 3.0]),
        np.array([6.0, 2.0, 0.0]),
        np.array([5.0, 0.0, 3.0]),
        np.array([5.0, 3.0, 0.0]),
    )
    assert s.over_call > 0.0 and s.under_call > 0.0, "the fixture must err both ways"
    np.testing.assert_allclose(s.over_call - s.under_call, s.net_err, atol=1e-12)
    np.testing.assert_allclose(s.over_call + s.under_call, s.abs_err, rtol=1e-12)

    one_way = P0.score_axis(
        np.array([5.0, 6.0]),
        np.array([5.0, 4.0]),
        np.array([3.0, 2.0]),
        np.array([7.0, 8.0]),
    )
    assert one_way.under_call == 0.0 and one_way.over_call > 0.0
