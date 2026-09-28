"""Falsification gates for ``scripts/design/_oracle_arms.py`` — the oracle cache and the per-object scorer
``calibration_vs_oracle.py`` reads.

Every number the scorer prints is a difference between two per-object arrays, so every way of getting
that wrong is a way of getting a *plausible* answer. These gates are the ways.

EACH GATE HERE CARRIES ITS OWN PERTURBATION, in the same test, because a gate that has never been
watched to fire has not been written yet (TRAPS: perturb-every-gate). Reading a gate is not evidence;
each test below breaks the thing it guards and asserts the guard notices.

The toy serves the cache and sum-to-full gates; the scoring gates run on synthetic arrays. A
single-reference toy is deliberately NOT enough to judge the deposit path, since a single-reference
index hides ref-id-space mismatches. The deposit path has its own gates in ``tests/native/`` and its
truth-scored instruments run on the panel.
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
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario


def _load_sibling(name: str):
    """Import a ``scripts/design/`` instrument by path.

    ``scripts/`` is not a package and must not become one — it is a toolkit of instruments, not an
    importable library, and putting an ``__init__.py`` in it would make every script a public API
    surface with a compatibility obligation. Loading by path keeps the dependency one-directional:
    the test knows about the script, the script knows nothing about the test.

    The module goes into ``sys.modules`` BEFORE it is executed. ``@dataclass`` resolves its own
    class's ``__module__`` through that table, so a module that is not registered fails at class
    definition with an ``AttributeError`` on ``None`` — nothing to do with the script. The scripts'
    directory goes on ``sys.path`` too, as running a script puts it there: the instruments import
    their shared helper ``_shared`` by name.
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
    """gDNA + mature + nascent, with spliced reads.

    The staggered isoform boundary is load-bearing: ``boundary_spliced`` — a molecule that crossed a
    contiguous boundary having spliced *elsewhere* — can only be deposited where a region bound falls
    INSIDE another transcript's exon. A single-isoform gene has no such bound, so its spliced-boundary
    bank is identically zero and GATE 1's perturbation removes nothing.
    """
    wd = tmp_path_factory.mktemp("p0_orc")
    sc = Scenario("p0", genome_length=9000, seed=17, work_dir=wd / "sim")
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(600, 1100), (1800, 2300)], "abundance": 60},
            # the stagger must sit CLOSE to the sj: the fragment has to reach the boundary
            # contiguously AND reach the sj, so a bound 700 bp away is one no fragment spans.
            {"t_id": "t1b", "exons": [(1000, 1100), (1800, 1900)], "abundance": 30},
        ],
    )
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(4000, 4500), (5200, 5700)], "abundance": 40}])
    sc.add_gene(
        "g3",
        "+",
        [
            {"t_id": "t3", "exons": [(6600, 7100), (7800, 8300)], "abundance": 25},
            {"t_id": "t3b", "exons": [(6600, 7100), (7400, 7440), (7800, 8300)], "abundance": 12},
        ],
    )
    return sc.build_oracle(
        n_rna_fragments=4000,
        gdna_fraction=1.0,
        nrna_abundance=15.0,
        sim_config=ReadSimConfig(
            frag_mean=180,
            frag_std=30,
            frag_min=80,
            frag_max=400,
            read_length=90,
            strand_specificity=0.99,
            seed=17,
        ),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=240, frag_std=45),
    )


# ── GATE 0: the oracle cache HITS, is REFUSED when stale, and is still validated ──────────────────


def test_the_oracle_cache_hits_and_the_cached_path_is_STILL_VALIDATED(toy, tmp_path, monkeypatch):
    """The cache exists so a solver-debugging campaign can re-measure the panel without re-splitting
    every BAM. Three things must hold or it is a liability rather than a saving:

    1. a warm run must not touch the BAM at all;
    2. it must reproduce the cold build's arrays exactly;
    3. it must STILL run sum-to-full — a cached oracle that skipped validation would be a silently
       wrong truth source feeding every number downstream.

    PERTURBATION (1): make ``_split_bam`` raise, so a warm run that touched the BAM cannot pass.
    PERTURBATION (3): corrupt a cached partition on disk and require the warm load to abort.
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

    # the cached path is NOT exempt from sum-to-full
    npz = cache / "t" / "gdna" / "payload.npz"
    data = {k: v for k, v in np.load(npz).items()}
    data["region_contained_count"] = data["region_contained_count"].copy()
    data["region_contained_count"][0, 0] += 1
    np.savez_compressed(npz, **data)
    with pytest.raises(AssertionError, match="oracle INVALID"):
        mod.load_or_build_oracle(bam, index, cfg, tmp_path / "w3", "t", full, cache)


def test_a_cache_that_does_not_describe_this_SCAN_is_rebuilt_not_reused(toy, tmp_path, monkeypatch):
    """A cache keyed to a different scan configuration is a different tally. It must be REBUILT,
    never silently reused — and never propagated as an error either, since a miss is normal.

    PERTURBATION: populate the cache under the default scan config, then ask for it under a changed
    one and require the BAM to be re-split.

    ``full`` is re-scanned under the changed config too, and that is not tidying: reusing the
    default-config ``full`` makes sum-to-full reject the result outright, because two scan configs
    are two different tallies, and the test would then be asserting the wrong thing. It also shows
    the guarantee is belt-and-braces — the cache KEY refuses a stale partition, and the IDENTITY
    independently refuses a ``full`` that does not match the partitions, even when the key was
    bypassed.
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

    # ...and the rebuilt cache now serves the CHANGED config without re-splitting again.
    mod.load_or_build_oracle(bam, index, changed, tmp_path / "w3", "t", changed_full, cache)
    assert calls["n"] == 1, "the rebuilt cache did not hit on the next call"


# ── GATE 1: T's per-axis totals are the FULL payload's totals ─────────────────────────────────────


def test_T_totals_equal_the_full_payload_PER_AXIS(toy, tmp_path):
    """``check_same_basis`` must hold between T and the payload it claims to partition, **per axis**.

    Per axis, never pooled: ``n_regions`` and ``n_boundaries`` differ by only ``n_refs``, so an error on
    one axis cancelling an equal and opposite one on the other is not far-fetched.

    PERTURBATION: drop the spliced term from ``count_rna_boundary``. That is the exact schema mistake
    ``override_masses`` exists to avoid — ``chain_boundary_deconv`` builds ``rna = (1−f_g)·unspliced +
    spliced``, so a T without the spliced term is on a different basis from every P.
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


# ── GATE 2: an arm and its truth are on the same basis, per object ────────────────────────────────


def test_an_arm_on_a_DIFFERENT_BASIS_is_refused():
    """Every arm's per-object total must equal T's per-object total. This is what makes ``f_g``
    comparable at all: two arrays of fractions over different denominators subtract to a number that
    means nothing.

    PERTURBATION (a): score against another axis' truth (a different object count).
    PERTURBATION (b): scale the arm's masses, i.e. put it on a different denominator.
    """
    g, r = np.array([4.0, 1.0, 3.0]), np.array([6.0, 2.0, 0.0])
    tg, tr = np.array([5.0, 0.0, 3.0]), np.array([5.0, 3.0, 0.0])
    P0.score_axis(g, r, tg, tr)  # the same per-object totals: holds
    with pytest.raises(ValueError, match="different axes"):
        P0.score_axis(g, r, tg[:-1], tr[:-1])
    with pytest.raises(ValueError, match="DIFFERENT BASES"):
        P0.score_axis(g * 2.0, r * 2.0, tg, tr)


# ── GATE 3: no data is ABSENT, never f_g = 0 ──────────────────────────────────────────────────────


def test_an_object_with_no_mass_is_ABSENT_not_a_confident_zero():
    """No data must be inert: never "100 % gDNA", and never its mirror "0 % gDNA" either. Most regions
    in any real index carry no fragments at all, so a scorer that turns 0/0 into a number reports a
    beautiful answer for the majority of the genome.

    PERTURBATION: the mass-weighted mean is *blind* to this by construction (a zero-mass object gets
    zero weight), so the gate is on the COUNT of scored objects — which is exactly where a floored 0
    would hide.
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


# ── GATE 4: the directional split is reported and does not cancel ─────────────────────────────────


def test_the_directional_split_is_reported_and_the_net_is_their_DIFFERENCE():
    """The library-level number looks far better than the per-object answer
    because a large under-call sits next to a large over-call. Reporting only the net is what makes
    that invisible, so the two directions are separate fields and their relationship is an identity.

    PERTURBATION: an arm whose errors all point one way must have one direction exactly zero — if
    both are always populated the split is measuring noise, not direction.
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
