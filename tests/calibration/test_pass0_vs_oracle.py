"""Falsification gates for ``scripts/design/pass0_vs_oracle.py`` — Stage B step 1's instrument.

The instrument answers: does calibration's PRIOR-FREE first solve find what the payload actually
contains, and where it does not, was the information destroyed by the accumulator or missed by the
solver? Every number it prints is a difference between two per-object arrays, so every way of getting
that wrong is a way of getting a *plausible* answer. These gates are the ways.

⛔ **EACH GATE HERE CARRIES ITS OWN PERTURBATION**, in the same test, because a gate that has never
been watched to fire has not been written yet (``test_accumulator_worker_determinism`` learned this the
expensive way, and so did the four-pool gDNA model at 8 perturbations / 1 blind gate). Reading a gate
is not evidence; each test below breaks the thing it guards and asserts the guard notices.

⚠ The scenario is a single-reference toy. That is enough for the SCORING logic these gates cover — it
is arithmetic over two per-object arrays — and it is deliberately NOT enough to judge the deposit path
(a single-reference index hides ref-id-space mismatches). The deposit path has its own gates in
``tests/native/`` and its truth-scored instruments run on the panel.
"""

from __future__ import annotations

import dataclasses
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from rigel.calibration.region_init import has_own_composition_evidence
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import CalibrationConfig, PipelineConfig
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario


def _load_sibling(name: str):
    """Import a ``scripts/design/`` instrument by path.

    ⚠ ``scripts/`` is not a package and must not become one — it is a toolkit of instruments, not an
    importable library, and putting an ``__init__.py`` in it would make every script a public API
    surface with a compatibility obligation. Loading by path keeps the dependency one-directional:
    the test knows about the script, the script knows nothing about the test.

    ⚠ The module goes into ``sys.modules`` BEFORE it is executed. ``@dataclass`` resolves its own
    class's ``__module__`` through that table, so a module that is not registered fails at class
    definition with an ``AttributeError`` on ``None`` — nothing to do with the script.
    """
    import sys

    key = name[:-3]
    if key not in _MODULES:
        path = Path(__file__).resolve().parents[2] / "scripts" / "design" / name
        spec = importlib.util.spec_from_file_location(key, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        _MODULES[key] = module
        spec.loader.exec_module(module)
    return _MODULES[key]


_MODULES: dict = {}
P0 = _load_sibling("pass0_vs_oracle.py")


# ── the toy ──────────────────────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """gDNA + mature + nascent, with spliced reads (``calibrate`` refuses a library without them).

    ⭐ Two structures here are load-bearing, and without them two gates below are vacuous:

    * **staggered isoform boundaries.** ``boundary_spliced`` — a molecule that crossed a contiguous boundary
      having spliced *elsewhere* — can only be deposited where a region_bound falls INSIDE another
      transcript's exon. A single-isoform gene has no such region_bound, so its spliced-boundary bank is
      identically zero and GATE 2's perturbation removes nothing.
    * **a region shorter than the minimum fragment length.** ``region_contained`` requires the fragment
      to fit, so a 40 bp region can never hold one. That is what gives the toy genuinely EMPTY objects
      — everything else here is covered ~100x deep — and GATE 6 is about exactly those.
    """
    wd = tmp_path_factory.mktemp("p0_orc")
    sc = Scenario("p0", genome_length=9000, seed=17, work_dir=wd / "sim")
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(600, 1100), (1800, 2300)], "abundance": 60},
            # ⚠ the stagger must sit CLOSE to the sj: the fragment has to reach the boundary
            # contiguously AND reach the sj, so a region_bound 700 bp away is one no fragment spans.
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


@pytest.fixture(scope="module")
def measured(toy, tmp_path_factory):
    """One full run of the instrument on the toy — the object every gate below interrogates.

    ⭐ The two ``truth`` pmfs handed in are DELIBERATELY separated far beyond anything the toy
    realises. C_input's job here is to prove the lever is wired, and a lever that changes nothing
    proves nothing (a capture A/B once "passed" with the feature removed because its two arms shared
    no random input).
    """
    return P0.measure_condition(
        bam=str(toy.bam_path),
        index=toy.index,
        pipeline_config=PipelineConfig(),
        calibration_config=CalibrationConfig(),
        work_dir=tmp_path_factory.mktemp("p0_split"),
        tag="toy",
        truth_pmfs=lambda size: (_spike_pmf(90, size + 1), _spike_pmf(330, size + 1)),
    )


def _spike_pmf(mean: int, size: int) -> np.ndarray:
    """A narrow triangular pmf centred on ``mean`` — a length model with a known, wide separation."""
    w = np.arange(size, dtype=np.float64)
    p = np.maximum(30.0 - np.abs(w - mean), 0.0)
    return p / p.sum()


# ── GATE 0: the oracle cache HITS, is REFUSED when stale, and is still validated ──────────────────


def test_the_oracle_cache_hits_and_the_cached_path_is_STILL_VALIDATED(toy, tmp_path, monkeypatch):
    """⭐ The cache exists so a solver-debugging campaign can re-measure the panel without re-splitting
    every BAM. Three things must hold or it is a liability rather than a saving:

    1. a warm run must not touch the BAM at all;
    2. it must reproduce the cold build's arrays exactly;
    3. ⛔ it must STILL run sum-to-full — a cached oracle that skipped validation would be a silently
       wrong truth source feeding every number downstream.

    PERTURBATION (1): make ``_split_bam`` raise, so a warm run that touched the BAM cannot pass.
    PERTURBATION (3): corrupt a cached partition on disk and require the warm load to abort.
    """
    import _oracle

    import pass0_vs_oracle as mod

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

    # ⛔ the cached path is NOT exempt from sum-to-full
    npz = cache / "t" / "gdna" / "payload.npz"
    data = {k: v for k, v in np.load(npz).items()}
    data["region_contained_count"] = data["region_contained_count"].copy()
    data["region_contained_count"][0, 0] += 1
    np.savez_compressed(npz, **data)
    with pytest.raises(AssertionError, match="oracle INVALID"):
        mod.load_or_build_oracle(bam, index, cfg, tmp_path / "w3", "t", full, cache)


def test_a_cache_that_does_not_describe_this_SCAN_is_rebuilt_not_reused(toy, tmp_path, monkeypatch):
    """⛔ A cache keyed to a different scan configuration is a different tally. It must be REBUILT,
    never silently reused — and never propagated as an error either, since a miss is normal.

    PERTURBATION: populate the cache under the default scan config, then ask for it under a changed
    one and require the BAM to be re-split.

    ⚠ ``full`` is re-scanned under the changed config too, and that is not tidying — the first draft
    of this test reused the default-config ``full`` and sum-to-full rejected the result outright.
    That is the identity doing its job (two scan configs are two different tallies), but it meant the
    test was asserting the wrong thing. It also shows the guarantee is belt-and-braces: the cache KEY
    refuses a stale partition, and the IDENTITY independently refuses a ``full`` that does not match
    the partitions, even when the key was bypassed.
    """
    import _oracle

    import pass0_vs_oracle as mod

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


# ── GATE 1: the oracle's sum-to-full identity actually runs on the condition under test ───────────


def test_a_corrupted_partition_ABORTS_the_measurement(toy, tmp_path, monkeypatch):
    """The oracle's trustworthiness rests entirely on sum-to-full. This asserts the instrument is
    behind that gate rather than beside it — that a broken partition stops the measurement instead of
    producing a plausible table.

    PERTURBATION: add one deposit to one partition's payload, so the parts no longer sum to the full.
    """
    import _oracle

    real = _oracle._scan_payload
    calls = {"n": 0}

    def corrupting(bam, index, cfg):
        payload = real(bam, index, cfg)
        calls["n"] += 1
        if calls["n"] == 1:  # the first partition scanned
            payload.region_contained_count[0, 0] += 1
        return payload

    monkeypatch.setattr(_oracle, "_scan_payload", corrupting)
    with pytest.raises(AssertionError, match="oracle INVALID"):
        P0.measure_condition(
            bam=str(toy.bam_path),
            index=toy.index,
            pipeline_config=PipelineConfig(),
            calibration_config=CalibrationConfig(),
            work_dir=tmp_path,
            tag="corrupt",
        )


# ── GATE 2: T's per-axis totals are the FULL payload's totals ─────────────────────────────────────


def test_T_totals_equal_the_full_payload_PER_AXIS(measured, toy):
    """``check_same_basis`` must hold between T and the payload it claims to partition, **per axis**.

    ⚠ Per axis, never pooled: ``n_regions`` and ``n_boundaries`` differ by only ``n_refs``, so an error on
    one axis cancelling an equal and opposite one on the other is not far-fetched.

    PERTURBATION: drop the spliced term from ``mass_rna_boundary``. That is the exact schema mistake
    ``override_masses`` exists to avoid — ``chain_boundary_deconv`` builds ``rna = (1−f_g)·unspliced +
    spliced``, so a T without the spliced term is on a different basis from every P.
    """
    ra = RegionArrays.from_frame(toy.index.regions_df, toy.index.ref_name_to_id)
    full = CalibrationSubstrate.from_payload(measured.oracle.full, ra)
    P0.check_same_basis("T", measured.truth, full)  # holds

    spliced = np.asarray(full.boundary_spliced.count, np.float64).sum(axis=1)
    assert spliced.sum() > 0, "the toy must EXERCISE the spliced bank or this perturbation is inert"
    broken = dataclasses.replace(
        measured.truth, mass_rna_boundary=measured.truth.mass_rna_boundary - spliced
    )
    with pytest.raises(ValueError, match="boundary"):
        P0.check_same_basis("T", broken, full)


# ── GATE 3: P and T are on the same basis, per object ─────────────────────────────────────────────


def test_P_and_T_are_on_the_SAME_BASIS_per_object(measured, toy):
    """Every arm's per-object total must equal T's per-object total on each axis. This is what makes
    ``f_g`` comparable at all: two arrays of fractions over different denominators subtract to a
    number that means nothing.

    PERTURBATION (a): score the region axis against the boundary axis' truth.
    PERTURBATION (b): scale one arm's masses, i.e. put P on a different denominator.
    """
    ra = RegionArrays.from_frame(toy.index.regions_df, toy.index.ref_name_to_id)
    full = CalibrationSubstrate.from_payload(measured.oracle.full, ra)
    for name, arm in measured.arms.items():
        P0.check_same_basis(name, arm, full)  # holds for every arm

    p = measured.arms["pass0"]
    with pytest.raises(ValueError):
        P0.score_axis(
            p.mass_gdna_region,
            p.mass_rna_region,
            measured.truth.mass_gdna_boundary,
            measured.truth.mass_rna_boundary,
        )
    with pytest.raises(ValueError):
        P0.score_axis(
            p.mass_gdna_region * 2.0,
            p.mass_rna_region * 2.0,
            measured.truth.mass_gdna_region,
            measured.truth.mass_rna_region,
        )


# ── GATE 4: calib_refit_iters=0 IS pass-0 ─────────────────────────────────────────────────────────


def test_refit_iters_zero_reproduces_debug_belief_pass0(measured, toy):
    """⭐ The config lever the instrument relies on must be the same quantity ``_debug`` exposes —
    checked ONCE, here, so the instrument can use the lever and nothing has to spelunk. Two ways of
    obtaining one quantity is how two modules come to disagree about it.

    PERTURBATION: ask for ONE refit iteration. The refit is the whole point of Phase 2, so a single
    iteration must move the answer; if it does not, the lever is inert and this test is vacuous.
    """
    from rigel.calibration.sweep import chain_region_deconv

    ra = RegionArrays.from_frame(toy.index.regions_df, toy.index.ref_name_to_id)
    substrate = CalibrationSubstrate.from_payload(measured.payload, ra)
    debug = measured.debug_final
    from_debug = chain_region_deconv(debug["chain"], debug["belief_pass0"], substrate).gdna_mass

    np.testing.assert_array_equal(measured.arms["pass0"].mass_gdna_region, from_debug)

    one = P0.calibrate_arm(
        measured.payload,
        measured.calibrate_kwargs,
        dataclasses.replace(CalibrationConfig(), calib_refit_iters=1),
    )
    assert not np.array_equal(one.mass_gdna_region, from_debug), (
        "one refit iteration left the answer byte-identical: the lever this instrument depends on "
        "does nothing, so the pass-0/final distinction it reports is fictional."
    )


# ── GATE 5: the undetermined class exists, is reported, and responds to the length gap ────────────


def test_the_UNDETERMINED_class_is_reported_and_tracks_the_length_gap(measured, toy):
    """C_info classifies each object as identified-or-not by the 2×2, and the undetermined class must
    be a reported CLASS carrying its own mass share — never averaged into the rest.

    PERTURBATIONS, in both directions:
      * two IDENTICAL length pmfs ⇒ the 2×2 has no separation anywhere ⇒ nothing is identified;
      * two widely separated pmfs ⇒ something is.
    """
    ra = RegionArrays.from_frame(toy.index.regions_df, toy.index.ref_name_to_id)
    substrate = CalibrationSubstrate.from_payload(measured.payload, ra)
    chain = measured.debug_final["chain"]

    size = int(measured.payload.max_length) + 1
    same = _spike_pmf(200, size)
    flat = P0.info_class_masks(chain, ra, substrate, same, same.copy())
    apart = P0.info_class_masks(chain, ra, substrate, _spike_pmf(150, size), _spike_pmf(330, size))

    for axis in ("region", "boundary"):
        live = ~flat[axis]["absent"]
        assert live.sum() > 0, f"the toy has no live {axis} objects; the gate would be vacuous"
        # equal means ⇒ the length channel carries EXACTLY zero information, at any depth.
        assert flat[axis]["identified"].sum() == 0
        assert (flat[axis]["undet_no_separation"] & live).sum() == int(live.sum())
        # a wide separation must rescue some of them, or the classifier is not reading the pmfs
        assert apart[axis]["identified"].sum() > 0

    # and the share is REPORTED, mass-weighted, for every class — an unreported class is an averaged
    # -away class.
    shares = measured.info_shares["region"]
    assert set(shares) == set(P0.INFO_CLASSES)
    assert abs(sum(shares.values()) - 1.0) < 1e-9


# ── GATE 6: no data is ABSENT, never f_g = 0 ──────────────────────────────────────────────────────


def test_an_object_with_no_mass_is_ABSENT_not_a_confident_zero(measured):
    """⭐ "No data" must be inert, never "100 % gDNA" — and its mirror, never "0 % gDNA", is just as
    wrong. Most regions in any real index carry no fragments at all, so a scorer that turns 0/0 into a
    number reports a beautiful answer for the majority of the genome.

    PERTURBATION: the mass-weighted mean is *blind* to this by construction (a zero-mass object gets
    zero weight), so the gate is on the COUNT of scored objects and on the class shares — which is
    exactly where a floored 0 would hide.
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

    # ...and on the real toy: exactly the objects with mass are scored, and there is at least one
    # without, or the boundary above proves nothing about the path the instrument actually runs.
    region = measured.scores["pass0"]["region"]["ALL"]
    live = (
        np.asarray(measured.truth.mass_gdna_region) + np.asarray(measured.truth.mass_rna_region)
    ) > 0
    assert region.n_scored == int(live.sum())
    assert (~live).sum() > 0, "the toy has no empty region; the gate would be vacuous"


# ── GATE 7: the solver classes are the SOLVER's, not a second opinion ─────────────────────────────


def test_the_solver_classes_match_the_composition_evidence_census(measured, toy):
    """⭐ ``composition_evidence_census.py`` already partitions slots into own-evidence / no-evidence /
    structurally-locked, reproducing ``region_init``'s own definitions. This instrument must ask the
    SAME question, or the project acquires two definitions of one class and they drift — which has
    already happened here once, between a docstring and the code it sat beside, for months.

    ⭐⭐ **RE-POINTED 2026-08-05, AND THE RE-POINTING IS THE FIX THIS GATE WAS ASKING FOR.** It used to
    compare two hand-written copies of ``1e-9`` (``P0._EPS == census._EPS``) — a drift DETECTOR for a
    definition that had three homes. ``region_init.has_own_composition_evidence`` is now the single home
    and every instrument imports it, so there is no longer anything to drift; the honest gate is that
    the census's partition IS the solver's predicate, evaluated by calling it. TRAPS: a-test-that-redefines' own
    prescription: one home, every consumer importing, rather than a gate policing copies.

    PERTURBATION: shift the evidence threshold off ``region_init``'s and watch the partition move.
    """
    census = _load_sibling("composition_evidence_census.py")
    cap = measured.debug_pass0["capture"]
    chain = measured.debug_pass0["chain"]

    tau = np.asarray(cap["_tau0_lam"], np.float64)
    is_region = np.asarray(chain.kind) == P0.REGION
    # ⛔ TRAPS: no-magic-numbers, on BOTH axes — `_type_belief` locks the class without consulting the axis. The earlier
    # `(~solvable) & is_region` filed every structurally-locked BOUNDARY as `relay_only`, i.e. as an object
    # whose answer came from its neighbours, when nothing was ever asked of it.
    census_lock = ~np.asarray(cap["free_pos"], bool) & ~np.asarray(cap["free_neg"], bool)
    # ⭐ the census asks the SOLVER's question by calling the solver's own predicate — imported, not
    #   restated. This boundary IS the census's boundary (`composition_evidence_census.census_one`).
    census_no_ev = ~census.has_own_composition_evidence(tau) & (~census_lock)

    slot = P0.solver_slot_classes(cap, chain)
    np.testing.assert_array_equal(slot["struct_lock"], census_lock)
    np.testing.assert_array_equal(slot["relay_only"], census_no_ev)
    np.testing.assert_array_equal(slot["own_evidence"], ~(census_lock | census_no_ev))
    # ⭐ and the census's predicate is production's, by identity rather than by an equal-constants
    #   assertion — the object itself, not two numbers that happen to agree.
    assert census.has_own_composition_evidence is has_own_composition_evidence

    # PERTURBATION 1: a threshold above every finite tau collapses own-evidence into relay-only.
    moved = P0.solver_slot_classes(cap, chain, eps=float(np.max(tau)) + 1.0)
    assert moved["own_evidence"].sum() == 0
    assert slot["own_evidence"].sum() > 0, "no slot has own evidence; the gate would be vacuous"

    # PERTURBATION 2: the region-only lock must be a DIFFERENT partition on this fixture, or the
    # correction above is untested and could silently revert.
    region_only_lock = (~np.asarray(cap["solvable"], bool)) & is_region
    assert not np.array_equal(region_only_lock, census_lock), (
        "the region-only and both-axes locks agree on this fixture, so it cannot demonstrate the "
        "defect — the scenario needs a G1 BOUNDARY carrying mass (an intergenic<->exon boundary)"
    )

    # ⭐ AND this partition must stay the SOLVER's gate, not a judgement about scoreability. The
    # audit answers the second question with a CURVE over sd(λ) = 1/√τ and no region_bound at all (a floor was
    # derived, implemented and refuted — τ is continuous). If a threshold ever appears here, the two
    # questions have been collapsed into one.
    sa = _load_sibling("solvability_audit.py")
    assert not hasattr(sa, "own_evidence_tau_floor"), (
        "a resolving-power THRESHOLD is back in the audit; it was refuted because tau is continuous "
        "across the region on 4 of 5 ladder conditions, so any floor is a tuned constant"
    )
    assert sa.SD_LAMBDA_DECADES[-1] == np.inf, (
        "the sd(lambda) curve must have an unbounded top band"
    )


# ── GATE 8: the classes partition the axis, and the error decomposes over them ────────────────────


def test_the_classes_PARTITION_the_mass_and_the_error(measured):
    """Per-class scoring is worthless if the classes overlap or leak: a class that double-counts an
    object reports its error twice, and a gap between classes hides error entirely. Both are silent.

    PERTURBATION: drop the LARGEST class from the sum and watch the identity break — dropping an
    arbitrary one is not a perturbation at all, because a class that happens to be empty on this
    scenario (``struct_lock`` is, on a toy with no intergenic pure-gDNA region) can be removed with
    no effect whatsoever, and the gate would then "pass" while testing nothing.
    """
    for kind, table in (("solver", measured.scores), ("info", measured.info_scores)):
        names = P0.SOLVER_CLASSES if kind == "solver" else P0.INFO_CLASSES
        for arm in table:
            for axis in ("region", "boundary"):
                per_class = table[arm][axis]
                whole = per_class["ALL"]
                classes = [per_class[c] for c in names]
                np.testing.assert_allclose(sum(c.mass for c in classes), whole.mass, rtol=1e-9)
                np.testing.assert_allclose(
                    sum(c.abs_err for c in classes), whole.abs_err, rtol=1e-9
                )
                np.testing.assert_allclose(
                    sum(c.net_err for c in classes), whole.net_err, atol=1e-6
                )
                assert sum(c.n_scored for c in classes) == whole.n_scored
                biggest = max(c.mass for c in classes)
                assert biggest > 0.0
                assert whole.mass - biggest != pytest.approx(whole.mass, rel=1e-9)


# ── GATE 9: the directional split is reported and does not cancel ─────────────────────────────────


def test_the_directional_split_is_reported_and_the_net_is_their_DIFFERENCE(measured):
    """⭐ The library-level number looks an order of magnitude better than the per-object answer
    because a large under-call sits next to a large over-call. Reporting only the net is what makes
    that invisible, so the two directions are separate fields and their relationship is an identity.

    PERTURBATION: an arm whose errors all point one way must have one direction exactly zero — if
    both are always populated the split is measuring noise, not direction.
    """
    for arm in measured.scores:
        for axis in ("region", "boundary"):
            s = measured.scores[arm][axis]["ALL"]
            np.testing.assert_allclose(s.over_call - s.under_call, s.net_err, atol=1e-6)
            np.testing.assert_allclose(s.over_call + s.under_call, s.abs_err, rtol=1e-9)

    one_way = P0.score_axis(
        np.array([5.0, 6.0]),
        np.array([5.0, 4.0]),
        np.array([3.0, 2.0]),
        np.array([7.0, 8.0]),
    )
    assert one_way.under_call == 0.0 and one_way.over_call > 0.0
