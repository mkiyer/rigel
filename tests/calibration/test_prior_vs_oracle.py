"""Falsification gates for ``scripts/design/prior_vs_oracle.py`` — the instrument that scores
calibration's ENDPOINT (``LocusPriors``) against the origin-split oracle.

Everything the instrument prints is a difference between two per-locus arrays of FRAGMENT COUNTS, so
every way of getting it wrong is a way of getting a *plausible* number: a lever that never fired, a
reference that quietly became the arm, a weight that rewards the solver for declining to answer, a
projection that loses mass off the end of a locus. These gates are those ways.

⛔ **EACH GATE CARRIES ITS OWN PERTURBATION**, in the same test. A gate that has never been watched to
fire has not been written yet — this project has found holes in five already-green gate sets that way,
once 7 of 9. Reading a gate is not evidence.

⚠ The scenario is a single-reference toy, which is enough for the SCORING and PLUMBING these gates
cover (arithmetic over two per-locus arrays, plus the override lever) and deliberately not enough to
judge the deposit path — that has its own gates in ``tests/native/`` and its truth-scored instruments
run on the panel.
"""

from __future__ import annotations

import dataclasses
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario

_MODULES: dict = {}


def _load_sibling(name: str):
    """Import a ``scripts/design/`` instrument by path.

    ⚠ ``scripts/`` is not a package and must not become one. ⚠ The module is registered in
    ``sys.modules`` BEFORE execution: ``@dataclass`` resolves its own class's ``__module__`` through
    that table and fails at class-definition time otherwise.
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


PV = _load_sibling("prior_vs_oracle.py")


# ── the toys ─────────────────────────────────────────────────────────────────────────────────────


def _scenario(name: str, seed: int, work_dir) -> Scenario:
    """Three genes, staggered isoforms, and one region shorter than the shortest fragment.

    ⭐ The stagger is load-bearing for the ``boundary_spliced`` bank (a contiguous crossing by a molecule
    that spliced elsewhere can only land where a region_bound falls inside another transcript's exon), and the
    short region is what gives the toy genuinely EMPTY objects — the population the NaN-not-zero gate is
    about. Both are the same structures ``test_pass0_vs_oracle`` relies on, for the same reasons.
    """
    sc = Scenario(name, genome_length=9000, seed=seed, work_dir=work_dir)
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(600, 1100), (1800, 2300)], "abundance": 60},
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
    return sc


_SIM = ReadSimConfig(
    frag_mean=180,
    frag_std=30,
    frag_min=80,
    frag_max=400,
    read_length=90,
    strand_specificity=0.99,
    seed=17,
)


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """Contaminated: gDNA at 1.0x the RNA fragment count, plus nascent RNA inside the introns."""
    sc = _scenario("pv", 17, tmp_path_factory.mktemp("pv_sim"))
    return sc.build_oracle(
        n_rna_fragments=4000,
        gdna_fraction=1.0,
        nrna_abundance=15.0,
        sim_config=_SIM,
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=240, frag_std=45),
    )


@pytest.fixture(scope="module")
def toy_zero_gdna(tmp_path_factory):
    """⛔ **THE OWNER-REQUIRED ZERO-gDNA CONTROL, at the instrument level.** Truth is exactly 0 at
    every locus, so there is nothing for a false positive to cancel against."""
    sc = _scenario("pv0", 23, tmp_path_factory.mktemp("pv0_sim"))
    return sc.build_oracle(
        n_rna_fragments=4000,
        gdna_fraction=0.0,
        nrna_abundance=15.0,
        sim_config=_SIM,
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=240, frag_std=45),
    )


def _measure(result, tmp_path, tag):
    return PV.measure_condition(
        str(result.bam_path),
        result.index,
        PipelineConfig(),
        tmp_path,
        tag,
        oracle_cache=None,
        drained_arm=False,
    )


@pytest.fixture(scope="module")
def measured(toy, tmp_path_factory):
    return _measure(toy, tmp_path_factory.mktemp("pv_run"), "toy")


@pytest.fixture(scope="module")
def measured_zero(toy_zero_gdna, tmp_path_factory):
    return _measure(toy_zero_gdna, tmp_path_factory.mktemp("pv0_run"), "toy_zero")


# ── GATE 1: the noop arm is byte-identical, AND it notices when it is not ────────────────────────


def _nudged_prior(measured, field, index, delta):
    """``assemble_priors`` re-run with ``+delta`` added to one element of one mass array."""
    cal = _rebuild_calibration(measured)
    arr = np.array(getattr(cal, field), dtype=np.float64, copy=True)
    arr[index] += delta
    return PV.PRIORS.assemble_priors(
        dataclasses.replace(cal, **{field: arr}), measured.region_arrays, measured.multi_loci
    )


def _moved(a, b) -> bool:
    return any(not np.array_equal(getattr(a, f), getattr(b, f)) for f in PV.PRIOR_FIELDS)


def test_the_noop_arm_is_byte_identical_and_the_lever_resolves_a_PICOFRAGMENT(measured):
    """⛔ TRAPS: byte-identity-gate. The whole instrument rests on ``dataclasses.replace`` being an inert way to swap
    the six mass arrays. If the replace dropped a field, O would differ from P for a reason that is
    not deconvolution error, and that bug would BE the headline number.

    ⭐ **The perturbation site had to be found by perturbation, and finding it is the gate's real
    content.** The first version nudged ``argmax(mass_gdna_region)`` and read "no effect" — because on
    any genome the largest gDNA region is INTERGENIC, and ``_project_regions_to_loci`` drops every
    region overlapping no locus. Correct behaviour, and it would have retired the gate as broken
    (TRAPS: could-the-arm-have-fired). So this asserts BOTH directions: in-locus moves, intergenic does not.

    ⚠ **Measured resolution: 1e-12 fragments at an in-locus region moves the prior by 1.0e-12; one ULP
    (≈3e-14 on a mass of 207) does not.** The projection is a plain summation, so a perturbation below
    the summand's own rounding is absorbed. 1e-12 fragments is twelve orders of magnitude below the
    unit the prior is denominated in, so the lever is not the limiting factor anywhere.
    """
    assert all(measured.noop_identical.values()), measured.noop_identical

    from rigel.calibration.priors import _RNA_SIGNATURE_BITS

    sig = np.asarray(measured.region_arrays.signature).astype(np.int64)
    in_locus = (sig & _RNA_SIGNATURE_BITS) != 0
    mass = np.asarray(measured.calibration.mass_gdna_region, np.float64)
    inside = int(np.argmax(np.where(in_locus, mass, -1.0)))
    outside = int(np.argmax(np.where(~in_locus, mass, -1.0)))
    assert mass[inside] > 0.0 and mass[outside] > 0.0, (
        "the toy has no gDNA at an in-locus region or no gDNA at an intergenic one — one half of this "
        "gate would be vacuous"
    )

    base = measured.priors["P"]
    assert _moved(_nudged_prior(measured, "mass_gdna_region", inside, 1e-12), base), (
        "1e-12 fragments at an in-locus region changed no prior — the lever cannot resolve an override"
    )
    # ⛔ The intergenic direction is asserted on the COUNT fields only. The locus projection drops
    #   intergenic regions from the counts — that is what this gate protects — but ``gdna_eff_len``
    #   may legitimately move: the eff-length contraction's ``_global_reference_density`` KDE reads
    #   EVERY region's gDNA density by design (its docstring: "detected from the data with no
    #   assumption about probe locations") and SNAPS to a real region's density, which can be the
    #   nudged intergenic region itself. Asserting all of ``PRIOR_FIELDS`` here conflated the two
    #   paths and held only while the snap happened to land elsewhere (exposed by the 2026-08-24
    #   reference-location deletion moving the mass landscape).
    nudged_out = _nudged_prior(measured, "mass_gdna_region", outside, 1.0)
    count_moved = any(
        not np.array_equal(getattr(nudged_out, f), getattr(base, f))
        for f in ("gdna_prior_count", "rna_prior_count")
    )
    assert not count_moved, (
        "a whole fragment at an INTERGENIC region changed a locus COUNT prior — the locus "
        "projection is no longer dropping regions that overlap no locus"
    )


def test_the_prior_reads_five_of_the_six_override_fields_and_provably_NOT_the_sj(measured):
    """⛔ TRAPS: an-ablation-that-never-ran, applied to the override itself. Five of the six arrays ``override_masses``
    writes must reach the prior — an override landing on a field nothing reads is an override that
    never ran, and it would read as "calibration is already correct on that channel".

    ⭐ **And the sixth must provably NOT reach it.** ``count_rna_sj`` is certified RNA: it is
    exported for QC and the prior deliberately does not read it, because the prior arbitrates only the
    UNSPLICED fragments and a locus whose RNA is fully spliced should have a near-zero
    ``rna_prior_count`` (owner ruling, 2026-07-30). That ruling lives in a docstring; this is what
    keeps it true.
    """
    base = measured.priors["P"]
    reads = {}
    for field in PV.OVERRIDE_FIELDS:
        site = _biggest_in_locus_site(measured, field)
        reads[field] = _moved(_nudged_prior(measured, field, site, 1.0), base)
    assert reads["count_rna_sj"] is False, (
        "the sj flux now reaches the prior — that is the owner ruling reversed, not a test "
        "failure to widen"
    )
    silent = [f for f, m in reads.items() if not m and f != "count_rna_sj"]
    assert not silent, f"override fields the prior never reads: {silent}"


def test_an_override_that_stops_writing_a_field_ABORTS_rather_than_scoring(measured, monkeypatch):
    """⛔ TRAPS: an-ablation-that-never-ran. If ``override_masses`` were changed to stop writing one of the six
    arrays, O would silently keep the SHIPPED value there and would be a hybrid of truth and estimate
    — which reads as "calibration is better than we thought" and is the most flattering possible bug.
    """
    cal = _rebuild_calibration(measured)
    full = measured.oracle.override_masses(measured.region_arrays)
    for dropped in PV.OVERRIDE_FIELDS:
        partial = {k: v for k, v in full.items() if k != dropped}
        monkeypatch.setattr(
            type(measured.oracle), "override_masses", lambda self, ra, _p=partial: dict(_p)
        )
        with pytest.raises(RuntimeError, match="OVERRIDE_FIELDS"):
            PV.oracle_priors(measured.oracle, cal, measured.region_arrays, measured.multi_loci)


# ── GATE 2: the lever COULD have fired ───────────────────────────────────────────────────────────


def test_the_oracle_lever_actually_MOVES_the_prior(measured):
    """⛔ TRAPS: could-the-arm-have-fired. "P equals O" would be the headline result of the whole campaign, so the
    one thing that must not produce it is a lever that did nothing. On a toy with real gDNA the two
    priors must differ at a substantial number of loci and by a substantial total.

    ⚠ Stated as a floor on the COUNT of differing loci as well as on the total, because a single
    enormous locus difference and a broad small one are different findings and only one of them
    proves the lever reaches the whole array.
    """
    p, o = measured.priors["P"], measured.priors["O"]
    differ = ~np.isclose(p.gdna_prior_count, o.gdna_prior_count, rtol=1e-9, atol=1e-9)
    assert differ.sum() >= 2, f"the oracle lever moved {differ.sum()} loci — it barely fired"
    assert o.gdna_prior_count.sum() > 0.0, "the oracle prior claims no gDNA on a contaminated toy"


# ── GATE 3: the two arms cannot be scored against a different locus partition ────────────────────


def test_scoring_against_a_DIFFERENT_locus_partition_raises(measured):
    """⛔ The locus partition is a function of the SCORING stage, not of the index — ``build_multi_loci``
    unions transcripts linked by scored fragments — so two runs of the pipeline can produce different
    numbers of loci. Index-aligning two such arrays would silently compare locus 7 of one run with
    locus 7 of another, which is not a small error.
    """
    p = measured.priors["P"].gdna_prior_count
    o = measured.priors["O"].gdna_prior_count
    PV.score_arm(p, o)  # the honest call
    with pytest.raises(ValueError, match="different locus partitions"):
        PV.score_arm(p, o[:-1])


# ── GATE 4: an empty prior makes NO claim, and that is NaN and not zero ──────────────────────────


def test_a_locus_with_no_prior_is_ABSENT_from_the_composition_not_a_confident_zero():
    """⛔ A ``(0, 0)`` prior is "this locus says nothing", not "this locus is pure RNA". Flooring it to
    ``phi = 0`` inflates the scored denominator with loci that have no answer to get wrong, and its
    mirror (``phi = 1`` from a ``0/0`` guarded the other way) reads as a confident all-gDNA claim —
    the exact shape that once seeded false gDNA into neighbouring exons.

    ⭐ The perturbation replaces the NaN with 0 and shows the scored COUNT moves, because that is the
    place the damage is visible: the mass-weighted mean is blind to it (a zero-scale locus carries
    zero weight), so a gate written on ``mwae_phi`` alone would pass with the bug in.
    """
    gdna = np.array([10.0, 0.0, 4.0])
    rna = np.array([30.0, 0.0, 0.0])
    phi, scale = PV.composition(gdna, rna)
    assert np.isnan(phi[1]) and scale[1] == 0.0
    assert phi[0] == pytest.approx(0.25) and phi[2] == pytest.approx(1.0)

    arm = _fake_priors(gdna, rna)
    ref = _fake_priors(np.array([12.0, 0.0, 3.0]), np.array([28.0, 0.0, 1.0]))
    honest = PV.score_composition(arm, ref)
    floored = PV.score_composition(
        _fake_priors(gdna, np.where((gdna + rna) == 0, 1.0, rna)),
        _fake_priors(
            ref.gdna_prior_count,
            np.where((ref.gdna_prior_count + ref.rna_prior_count) == 0, 1.0, ref.rna_prior_count),
        ),
    )
    assert honest["n_scored"] == 2
    assert floored["n_scored"] == 3, "flooring an empty prior must be visible in the scored COUNT"


# ── GATE 5: phi and scale are SEPARATE axes, and the weight is the REFERENCE's ───────────────────


def test_rescaling_the_arm_moves_the_SCALE_and_leaves_phi_UNTOUCHED():
    """⭐ A prior can be right about the RATIO and wrong about the STRENGTH, or the reverse, and one
    number cannot say which. This gate proves the two reported axes are actually independent:
    multiplying both of the arm's counts by ``k`` leaves ``phi`` exactly where it was and must move
    ``scale_log10_ratio`` by exactly ``log10(k)``.

    ⛔ It also pins the WEIGHT to the reference. TRAPS: honesty-metrics-reward-ignorance: if ``mwae_phi`` were weighted by
    the ARM's own scale, a mechanism could improve it by shrinking its prior to nothing exactly at the
    loci it gets wrong — an accuracy metric that rewards saying less. Here ``k = 1e-6`` is that
    shrinkage taken to its limit, and ``mwae_phi`` must not move at all.
    """
    gdna = np.array([10.0, 1.0, 40.0])
    rna = np.array([30.0, 99.0, 10.0])
    ref = _fake_priors(np.array([20.0, 5.0, 10.0]), np.array([20.0, 95.0, 40.0]))
    base = PV.score_composition(_fake_priors(gdna, rna), ref)
    for k in (1e-6, 1e3):
        scaled = PV.score_composition(_fake_priors(gdna * k, rna * k), ref)
        assert scaled["mwae_phi"] == pytest.approx(base["mwae_phi"], rel=1e-12)
        assert scaled["weight"] == pytest.approx(base["weight"], rel=1e-12)
        assert scaled["scale_log10_ratio"] == pytest.approx(
            base["scale_log10_ratio"] + np.log10(k), abs=1e-12
        )


# ── GATE 6: gdna_eff_len is weighted like its consumer ───────────────────────────────────────────


def test_eff_len_error_ignores_loci_with_no_gDNA_and_notices_loci_with_some():
    """⚠ TRAPS: weight-it-like-the-consumer. ``gdna_eff_len`` divides the gDNA component's abundance and nothing else,
    so at a locus with no gDNA it is a number nothing reads. Weighting it by the locus TOTAL would
    report the error of an inert array — and on this panel most loci are that.

    The perturbation is two-sided, which is the point: corrupting the eff-len where the reference puts
    no gDNA must change nothing, and corrupting it where the reference puts a lot must change it.
    """
    ref = _fake_priors(
        np.array([0.0, 1000.0]), np.array([50.0, 50.0]), eff_len=np.array([500.0, 500.0])
    )
    base = PV.score_eff_len(
        _fake_priors(np.zeros(2), np.zeros(2), eff_len=np.array([500.0, 500.0])), ref
    )
    inert = PV.score_eff_len(
        _fake_priors(np.zeros(2), np.zeros(2), eff_len=np.array([5.0, 500.0])), ref
    )
    live = PV.score_eff_len(
        _fake_priors(np.zeros(2), np.zeros(2), eff_len=np.array([500.0, 5.0])), ref
    )
    assert base["w_rel_err"] == pytest.approx(0.0)
    assert inert["w_rel_err"] == pytest.approx(0.0), "a zero-gDNA locus must carry zero weight"
    assert live["w_rel_err"] > 0.9, "a locus carrying all the gDNA must dominate the weighted error"


# ── GATE 7: F conserves every fragment, and the residue is named ─────────────────────────────────


def test_the_fragment_truth_projection_LOSES_NOTHING_it_does_not_report(measured):
    """⛔ ``_project_regions_to_loci`` DROPS every region overlapping no locus — that is correct (an
    intergenic fragment belongs to no prior) and it is also the one place F could quietly lose mass
    and read as a smaller assembler error. So the identity ``Σ F + dropped == Σ region_start_count``
    must hold EXACTLY, per origin, and ``dropped`` must be reported rather than absorbed.

    ⭐ The perturbation removes one locus from the projection and watches the residue absorb exactly
    that locus's count — proving the identity is measuring the projection and not just restating a sum.
    """
    for origin, arm, drop in (
        ("gdna", measured.f_gdna, measured.f_dropped["gdna"]),
        ("rna", measured.f_rna_upper, measured.f_dropped["rna"]),
    ):
        if origin == "gdna":
            total = float(np.asarray(measured.oracle.parts["gdna"].region_start_count).sum())
        else:
            total = float(
                np.asarray(measured.oracle.parts["mrna"].region_start_count).sum()
                + np.asarray(measured.oracle.parts["nrna"].region_start_count).sum()
            )
        assert arm.sum() + drop == pytest.approx(total, rel=1e-9), origin
        assert drop >= 0.0, f"{origin}: a NEGATIVE residue means the projection invented fragments"

    short = measured.multi_loci[:-1]
    g2, r2, drop2 = PV.fragment_truth(measured.oracle, measured.region_arrays, short)
    assert g2.shape[0] == len(short)
    assert drop2["gdna"] >= measured.f_dropped["gdna"], (
        "removing a locus did not increase the dropped residue — the identity is not measuring the "
        "projection"
    )


def test_the_gdna_partition_carries_NO_spliced_deposit_so_F_gdna_needs_no_subtraction(measured):
    """⭐ This is why F is EXACT on the gDNA arm and only a bound on the RNA arm, and it is physics
    rather than a convention: gDNA does not splice, so there is no spliced sub-population inside
    ``region_start_count`` for the gdna partition to withhold.

    ⛔ The perturbation writes a single spliced deposit into the gdna partition and asserts
    ``OracleTruth`` refuses the whole oracle — because if it did not, F_gdna would silently become a
    bound too and the instrument's strongest claim would be false.
    """
    from _oracle import OracleTruth

    g = measured.oracle.parts["gdna"]
    assert int(np.asarray(g.boundary_spliced_count).sum()) == 0
    assert int(np.asarray(g.sj_count).sum()) == 0

    poisoned = np.array(g.boundary_spliced_count, copy=True)
    poisoned[0, 0] += 1
    fake = dataclasses.replace(g, boundary_spliced_count=poisoned)
    with pytest.raises(AssertionError, match="boundary_spliced_count"):
        OracleTruth.from_parts(measured.oracle.full, {**measured.oracle.parts, "gdna": fake})


# ── GATE 8: the ZERO-gDNA control ────────────────────────────────────────────────────────────────


def test_at_zero_gDNA_the_ORACLE_prior_is_identically_zero_and_the_shipped_one_is_scored_against_it(
    measured_zero,
):
    """⛔⛔ **THE OWNER-REQUIRED ZERO CONTROL.** With no gDNA in the library the oracle's gDNA mass is
    exactly 0 at every object, so ``O.gdna_prior_count`` must be exactly 0 at every locus — not small,
    not floored, zero. Anything the SHIPPED prior puts there is a false positive with nothing to
    cancel it, which is the only reading of that arm that is unambiguous.

    ⭐ The perturbation is the other direction and it is what makes the assertion non-vacuous: hand
    the same assembler a single fabricated gDNA fragment and the prior must come off zero. A gate that
    only ever sees zeros cannot tell "correct" from "the array is not wired".
    """
    o = measured_zero.priors["O"]
    assert float(np.asarray(o.gdna_prior_count).sum()) == 0.0, (
        "the oracle prior claims gDNA in a library that has none"
    )
    assert float(np.asarray(o.rna_prior_count).sum()) > 0.0, (
        "the oracle prior claims no RNA either — this arm is not wired, not correct"
    )

    cal = _rebuild_calibration(measured_zero)
    truth = measured_zero.oracle.override_masses(measured_zero.region_arrays)
    seeded = np.array(truth["mass_gdna_region"], copy=True)
    # the largest RNA region — one with real opportunity, so the mass is not dropped by
    # ``_mass_where_there_is_opportunity``
    i = int(np.argmax(np.asarray(truth["mass_rna_region"])))
    seeded[i] = 1.0
    with_one = PV.PRIORS.assemble_priors(
        dataclasses.replace(cal, **{**truth, "mass_gdna_region": seeded}),
        measured_zero.region_arrays,
        measured_zero.multi_loci,
    )
    assert float(np.asarray(with_one.gdna_prior_count).sum()) > 0.0, (
        "one fabricated gDNA fragment produced a prior of exactly zero — the zero above is the "
        "array being dead, not the library being clean"
    )


# ── GATE 9: a capture that never happened is an error, not a zero ────────────────────────────────


def test_a_run_that_never_reaches_assemble_priors_RAISES(measured, monkeypatch):
    """⛔ TRAPS: an-ablation-that-never-ran. ``quant_from_buffer`` returns early when there are no EM units, and a
    silently-absent capture would read as "a condition with no loci and therefore no error" — the most
    flattering possible failure of the harness.
    """
    monkeypatch.setattr(PV, "quant_from_buffer", lambda *a, **k: (None, None))
    with pytest.raises(RuntimeError, match="never called"):
        PV.capture_priors(None, None, None, None, None, None, None, None, PipelineConfig())


# ── GATE 10: the aggregate re-derives its rates, never averages them ─────────────────────────────


def test_the_stratum_aggregate_is_a_RATIO_OF_SUMS_not_a_mean_of_ratios():
    """⛔ TRAPS: never-pool-the-strata's third way. A panel condition at 10 M fragments and one at 10 k are not two
    equally-informative opinions about a rate, and averaging their ``rel`` values gives the shallow
    one equal say. The aggregate must recompute from the summed totals.

    ⭐ The two constructed rows differ by 1,000x in depth and have opposite-signed errors, so the mean
    of ratios and the ratio of sums are far apart and a gate that confused them could not pass by
    coincidence.
    """
    deep = PV.ArmScore(
        n_loci=1,
        n_claiming=1,
        total_arm=1_100_000.0,
        total_ref=1_000_000.0,
        net_err=100_000.0,
        abs_err=100_000.0,
        over_call=100_000.0,
        under_call=0.0,
    )
    shallow = PV.ArmScore(
        n_loci=1,
        n_claiming=1,
        total_arm=100.0,
        total_ref=1_000.0,
        net_err=-900.0,
        abs_err=900.0,
        over_call=0.0,
        under_call=900.0,
    )
    agg = PV._agg([deep, shallow])
    assert agg.rel == pytest.approx(100_900.0 / 1_001_000.0)
    assert agg.rel != pytest.approx((deep.rel + shallow.rel) / 2.0, rel=1e-3)
    assert agg.net_err == pytest.approx(99_100.0)
    assert agg.over_call == pytest.approx(100_000.0) and agg.under_call == pytest.approx(900.0)


# ── GATE 11: the directional split is reported and reconciles ────────────────────────────────────


def test_over_and_under_call_are_reported_separately_and_reconcile(measured):
    """⭐ The library-level figure is ``|Σ(a − a*)|`` and the per-locus answer is ``Σ|a − a*|``; when a
    large under-call sits next to a large over-call the first flatters the second by whatever
    ``cancellation`` reports. Both halves must therefore exist and must reconcile exactly::

        over − under == net        over + under == abs
    """
    for arm in ("gdna_prior_count", "rna_prior_count"):
        s = PV.score_arm(getattr(measured.priors["P"], arm), getattr(measured.priors["O"], arm))
        assert s.over_call - s.under_call == pytest.approx(s.net_err, rel=1e-9, abs=1e-6)
        assert s.over_call + s.under_call == pytest.approx(s.abs_err, rel=1e-9, abs=1e-6)
        assert s.over_call >= 0.0 and s.under_call >= 0.0


# ── GATE 12: the frag_id join ALIGNS, and a one-fragment slip is loud ────────────────────────────


def test_the_frag_id_join_is_gated_by_a_COUNT_IDENTITY_and_it_REFUSES_a_walk_that_slipped(toy):
    """⛔⛔ **THE Fo ARM'S ONE SILENT FAILURE MODE, AND WHY THE GATE IS ARITHMETIC.** ``frag_origin`` is
    indexed by the scanner's ``frag_id``; the walk re-derives that counter from the BAM. Slip by a
    single fragment and every unit still gets a *plausible* origin label, every total still looks like
    a count, and nothing is out of range to raise on.

    ⭐ **The gate is therefore an identity against the scanner's own counters**, not a smell test:
    ``stats.total`` is every record it read and ``stats.n_read_names`` is incremented once per qname
    group inside its worker, so it IS the number of ``frag_id``\\ s issued. Two monotone counters over
    one file that agree on both totals cannot have disagreed in the middle.

    ⛔ Perturbed in three directions — one record too many, one group too many, one group too few — and
    the un-perturbed identity is asserted too, so a guard that refused everything would not pass.
    """
    from _oracle import check_walk_alignment, frag_id_origins

    cfg = PipelineConfig()
    stats, _sm, _buf, _payload = scan_and_buffer(
        str(toy.bam_path),
        toy.index,
        dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(str(toy.bam_path))),
    )
    walk, diag = frag_id_origins(str(toy.bam_path), cfg.scan)
    check_walk_alignment(diag, stats)
    assert diag["n_groups"] == int(stats.n_read_names) > 0
    assert walk.shape[0] == diag["n_groups"]

    for field, delta in (("n_records", 1), ("n_groups", 1), ("n_groups", -1)):
        with pytest.raises(RuntimeError, match="does NOT reproduce"):
            check_walk_alignment({**diag, field: diag[field] + delta}, stats)


def test_the_SPLICED_gDNA_diagnostic_fires_on_a_BLOCK_SIZED_slip_and_is_blind_to_a_SMALL_one(
    measured,
):
    """⭐ The join's secondary diagnostic: gDNA does not splice, so a spliced unit labelled ``gdna`` is
    impossible physics and its count reads out a gross misalignment.

    ⛔⛔ **AND ITS SENSITIVITY IS MEASURED HERE RATHER THAN ASSUMED, because the first version of this
    gate asserted the opposite and failed.** The simulator writes each population as a CONTIGUOUS BLOCK,
    so BAM order has a handful of origin transitions (15 on a 10 M-fragment panel condition) and a
    one-fragment roll mislabels only the fragments sitting on those boundaries — a couple, none of them
    necessarily spliced. ⭐ So this test pins BOTH halves: a roll of one is invisible, a roll across a
    block is loud. That is why the hard gate is the count identity and not this.
    """
    d = measured.overlap.diag
    assert d["spliced_gdna_units"] == 0, (
        "a spliced unit is labelled gdna on an UNPERTURBED run — the join is already misaligned"
    )
    assert d["spliced_rna_units"] > 0, (
        "no spliced units at all: the detector has nothing to detect with and this gate is inert"
    )
    assert d["walk"]["n_transitions"] < d["walk"]["n_groups"] // 100, (
        "the origins are INTERLEAVED on this substrate, not blocked — then a one-fragment roll would "
        "be visible and the reasoning below no longer describes the panel"
    )

    def rolled(shift):
        return PV.overlap_truth(
            measured.multi_loci,
            PV.unit_origins(measured.units["frag_ids"], np.roll(measured.frag_origin, shift)),
            measured.units["is_spliced"],
            measured.units["n_units"],
            d["walk"],
        ).diag["spliced_gdna_units"]

    assert rolled(1) == 0, (
        "a one-fragment roll DID show up — good news for the diagnostic, but then the blocked-origin "
        "reasoning in this docstring is wrong and must be rewritten, not widened"
    )
    half = d["walk"]["n_groups"] // 2
    assert rolled(half) > 0, (
        "rolling every origin label across a population block did not put gDNA on a single spliced "
        "unit — the diagnostic cannot see even a gross slip and is worth nothing"
    )


# ── GATE 13: a filtered record does NOT advance frag_id, and the config decides which ─────────────


def _rewrite_bam(src: Path, dst: Path, *, insert_after: int, flag: int):
    """``src`` with ONE synthetic record inserted after group ``insert_after``, carrying ``flag``.

    ⭐ A fresh, PARSEABLE qname, so the only difference between counting it and skipping it is the
    off-by-one — not a crash in ``parse_origin`` that would pass the test for the wrong reason.
    """
    import pysam

    with pysam.AlignmentFile(str(src), "rb") as fin:
        recs = list(fin)
        header = fin.header
    groups, seen = [], None
    for r in recs:
        if r.query_name != seen:
            seen = r.query_name
            groups.append([])
        groups[-1].append(r)
    ghost = recs[0].__copy__()
    ghost.query_name = "gdna:ref0:1000-1100:+:987654"
    ghost.flag = recs[0].flag | flag
    out = []
    for i, g in enumerate(groups):
        out += g
        if i == insert_after:
            out.append(ghost)
    with pysam.AlignmentFile(str(dst), "wb", header=header) as fo:
        for r in out:
            fo.write(r)
    return len(groups)


def test_a_FILTERED_record_does_not_advance_frag_id_and_skip_duplicates_decides_which_are(
    toy, tmp_path
):
    """⛔ The scanner rejects QC-fail / unmapped / duplicate records in pass 1 **before** it stamps a
    ``frag_id``, so a walk that counted them would shift every later fragment's label. And *which*
    records are rejected is a CONFIG question — ``skip_duplicates`` — which is why
    ``frag_id_origins`` takes the scan config rather than assuming.

    ⭐ Three arms over the same poisoned BAM: a QC-fail ghost (always filtered, mapping unchanged), the
    same ghost as a duplicate under ``skip_duplicates=True`` (filtered, unchanged), and under
    ``skip_duplicates=False`` (counted, and every later label shifts). The third arm is the
    perturbation: it proves the config argument is load-bearing and not decoration.
    """
    from _oracle import frag_id_origins

    scan = PipelineConfig().scan
    base, _ = frag_id_origins(str(toy.bam_path), scan)

    qcfail = tmp_path / "qcfail.bam"
    _rewrite_bam(Path(toy.bam_path), qcfail, insert_after=3, flag=0x200)
    got, diag = frag_id_origins(str(qcfail), scan)
    assert np.array_equal(got, base), "a QC-fail record advanced frag_id"
    assert diag["n_filtered"] == 1 and diag["n_groups"] == base.shape[0]

    dup = tmp_path / "dup.bam"
    _rewrite_bam(Path(toy.bam_path), dup, insert_after=3, flag=0x400)
    kept, _ = frag_id_origins(str(dup), dataclasses.replace(scan, skip_duplicates=True))
    assert np.array_equal(kept, base), "a duplicate advanced frag_id under skip_duplicates=True"

    counted, diag_c = frag_id_origins(str(dup), dataclasses.replace(scan, skip_duplicates=False))
    assert diag_c["n_filtered"] == 0
    assert counted.shape[0] == base.shape[0] + 1, (
        "skip_duplicates=False did not count the duplicate — the config argument is inert, and a walk "
        "that ignores it can disagree with the scan it is joined to"
    )
    assert not np.array_equal(counted[:5], base[:5]) or not np.array_equal(
        counted[4:], base[3:-1]
    ), "the extra group did not shift any label, so this BAM cannot detect a miscount"


def test_an_UNPAIRED_record_makes_the_walk_REFUSE_rather_than_count_it(toy, tmp_path):
    """⛔ The production scanner throws on an unpaired read, so a walk that tolerated one would be
    counting groups no scan ever made. ⭐ The perturbation clears the PAIRED bit on one record."""
    from _oracle import frag_id_origins

    single = tmp_path / "single.bam"
    _rewrite_bam(Path(toy.bam_path), single, insert_after=3, flag=0)
    import pysam

    with pysam.AlignmentFile(str(single), "rb") as fin:
        recs, header = list(fin), fin.header
    recs[0].flag = recs[0].flag & ~0x1
    with pysam.AlignmentFile(str(single), "wb", header=header) as fo:
        for r in recs:
            fo.write(r)
    with pytest.raises(AssertionError, match="unpaired"):
        frag_id_origins(str(single), PipelineConfig().scan)


# ── GATE 14: every unit is counted ONCE and the residue is named ──────────────────────────────────


def test_Fo_counts_every_unit_ONCE_and_the_non_candidate_residue_RECONCILES(measured):
    """⛔ ``Fo`` is a per-locus fragment COUNT, so the two ways to get it wrong are to count a unit
    twice (a unit claimed by two loci) and to lose one silently (a unit claimed by none). Both are
    checkable against totals the arm does not compute:

        Σ Fo[gdna] + Σ Fo[rna] + orphan_units == n_units          nothing double-counted, nothing lost
        Σ Fo[origin] + nonunit_fragments[origin] == the library's own total for that origin

    ⭐ The perturbation drops one locus and watches BOTH residues absorb exactly its units — an
    identity that merely restated a sum could not do that.
    """
    d = measured.overlap.diag
    total = measured.overlap.gdna.sum() + measured.overlap.rna_all.sum()
    assert total + d["orphan_units"] == d["n_units"]
    for origin, arm in (("gdna", measured.overlap.gdna), ("rna", measured.overlap.rna_all)):
        lib = (
            d["walk"]["totals"]["gdna"]
            if origin == "gdna"
            else d["walk"]["totals"]["mrna"] + d["walk"]["totals"]["nrna"]
        )
        assert arm.sum() + d["nonunit_fragments"][origin] == pytest.approx(lib, rel=1e-12)
        assert d["nonunit_fragments"][origin] >= 0.0, (
            f"{origin}: more units than fragments — a unit is being counted twice"
        )

    dropped = measured.multi_loci[-1]
    short = PV.overlap_truth(
        measured.multi_loci[:-1],
        PV.unit_origins(measured.units["frag_ids"], measured.frag_origin),
        measured.units["is_spliced"],
        measured.units["n_units"],
        d["walk"],
    )
    lost = len(dropped.unit_indices)
    assert lost > 0, "the dropped locus had no units, so this perturbation tests nothing"
    assert short.diag["orphan_units"] == d["orphan_units"] + lost
    assert (short.gdna.sum() + short.rna_all.sum()) == total - lost


# ── GATE 15: the RNA arm's two populations are the SPLICE BIT and nothing else ────────────────────


def test_the_RNA_arm_splits_on_is_spliced_and_the_two_populations_RECONCILE(measured):
    """⭐ ``rna_prior_count`` withholds spliced mass, so ``Fo`` reports two RNA arrays: the assembler's
    target (unspliced units) and the EM's own RNA evidence (all units). ⛔ They must differ by exactly
    the spliced RNA units and by nothing else.

    ⭐ The perturbation replaces ``is_spliced`` with all-False and then all-True: the first must collapse
    the two arrays onto each other element-wise, the second must empty the unspliced one. A split driven
    by anything other than that bit survives one of the two.
    """
    o = measured.overlap
    assert np.all(o.rna_unspliced <= o.rna_all)
    assert o.rna_all.sum() - o.rna_unspliced.sum() == pytest.approx(
        o.diag["spliced_rna_units"] - 0.0, rel=1e-12
    ), "the two RNA populations do not differ by the spliced unit count"

    args = (
        measured.multi_loci,
        PV.unit_origins(measured.units["frag_ids"], measured.frag_origin),
    )
    n = measured.units["n_units"]
    none_spliced = PV.overlap_truth(*args, np.zeros(n, bool), n, o.diag["walk"])
    assert np.array_equal(none_spliced.rna_unspliced, none_spliced.rna_all)
    assert none_spliced.diag["spliced_rna_units"] == 0
    all_spliced = PV.overlap_truth(*args, np.ones(n, bool), n, o.diag["walk"])
    assert all_spliced.rna_unspliced.sum() == 0.0
    assert np.array_equal(all_spliced.rna_all, o.rna_all), (
        "the splice bit moved the ALL-RNA array — it must only split it"
    )


# ── GATE 16: the join ABORTS on a frag_id the walk never issued ──────────────────────────────────


def test_a_unit_frag_id_the_WALK_NEVER_ISSUED_aborts_instead_of_indexing(measured):
    """⛔ ``frag_origin`` is indexed BY ``frag_id``. A walk of the wrong BAM, or one that grouped
    differently, yields an array of the wrong length — and numpy would wrap a negative index silently
    and raise a bare ``IndexError`` for a large one, neither of which says "the join is broken".

    ⭐ Falsified in both directions, and the in-range case is asserted too: a guard that rejected
    everything would also pass the two raises.
    """
    origins = np.asarray([2, 0, 1], np.int8)
    assert PV.unit_origins(np.asarray([0, 2]), origins).tolist() == [2, 1]
    with pytest.raises(RuntimeError, match="frag_id"):
        PV.unit_origins(np.asarray([0, 3]), origins)
    with pytest.raises(RuntimeError, match="frag_id"):
        PV.unit_origins(np.asarray([-1, 0]), origins)
    # and the real arrays are in range, so the guard is not the reason the arm looks healthy
    PV.unit_origins(measured.units["frag_ids"], measured.frag_origin)


# ── GATE 17: Fo follows the SHIPPED unit→locus map, and the prior does NOT ───────────────────────


def test_Fo_is_keyed_by_the_SHIPPED_unit_indices_and_a_SWAP_moves_the_counts(measured):
    """⭐ ``MultiLocus.unit_indices`` is the array ``locus_partition`` scatters by, so it — and not any
    genomic-overlap rule invented here — decides which locus's prior a fragment's evidence lands in.
    ⛔ The perturbation swaps two loci's unit sets and demands the counts swap with them. A tally
    driven by geometry instead would not move.
    """
    ml = measured.multi_loci
    order = np.argsort([-len(m.unit_indices) for m in ml])
    a, b = int(order[0]), int(order[1])
    assert len(ml[a].unit_indices) and len(ml[b].unit_indices)
    swapped = list(ml)
    swapped[a] = dataclasses.replace(ml[a], unit_indices=ml[b].unit_indices)
    swapped[b] = dataclasses.replace(ml[b], unit_indices=ml[a].unit_indices)
    got = PV.overlap_truth(
        swapped,
        PV.unit_origins(measured.units["frag_ids"], measured.frag_origin),
        measured.units["is_spliced"],
        measured.units["n_units"],
        measured.overlap.diag["walk"],
    )
    assert got.gdna[a] == measured.overlap.gdna[b]
    assert got.gdna[b] == measured.overlap.gdna[a]
    assert got.gdna[a] != got.gdna[b], "the two loci carry equal counts, so a swap proves nothing"


def test_assemble_priors_is_BLIND_to_unit_indices_so_Fo_is_not_circular(measured):
    """⛔⛔ ``Fo`` is built from ``unit_indices`` and scored against a prior built by
    ``assemble_priors``. If that function read a unit count, "the assembler reproduces the EM's own
    count" would be a tautology rather than a result.

    ⭐ Behavioural, not a source grep: every locus's ``unit_indices`` is emptied and the three prior
    arrays must come back BYTE-identical. ⚠ And the same perturbation is shown to move ``Fo`` to
    nothing, so the invariance is the assembler's and not the perturbation's failure to bite.
    """
    from rigel.calibration.priors import assemble_priors

    blinded = [
        dataclasses.replace(m, unit_indices=np.zeros(0, dtype=m.unit_indices.dtype))
        for m in measured.multi_loci
    ]
    ref = assemble_priors(measured.calibration, measured.region_arrays, measured.multi_loci)
    got = assemble_priors(measured.calibration, measured.region_arrays, blinded)
    for field in PV.PRIOR_FIELDS:
        assert np.array_equal(getattr(got, field), getattr(ref, field)), (
            f"{field} moved when unit_indices was emptied — assemble_priors READS the unit count and "
            "the Fo comparison is circular"
        )
    empty = PV.overlap_truth(
        blinded,
        PV.unit_origins(measured.units["frag_ids"], measured.frag_origin),
        measured.units["is_spliced"],
        measured.units["n_units"],
        measured.overlap.diag["walk"],
    )
    assert empty.gdna.sum() == 0.0 and empty.diag["orphan_units"] == measured.units["n_units"], (
        "emptying unit_indices did not move Fo either — the perturbation does not bite"
    )


# ── helpers ──────────────────────────────────────────────────────────────────────────────────────


def _in_locus_regions(measured) -> np.ndarray:
    """``bool[N]`` — regions the locus projection keeps. An intergenic region carries no exon/intron bit
    and is DROPPED, so a perturbation there is inert by design (see the noop gate)."""
    from rigel.calibration.priors import _RNA_SIGNATURE_BITS

    sig = np.asarray(measured.region_arrays.signature).astype(np.int64)
    return (sig & _RNA_SIGNATURE_BITS) != 0


def _biggest_in_locus_site(measured, field):
    """The largest element of ``field`` that the locus projection actually reaches.

    ⛔ Boundary-indexed arrays are selected through the SHIPPED ``_boundary_locus_shares`` rather than a local
    rule: a locus's boundaries are the boundaries that TOUCH its regions, which is exactly the decision that
    function exists to make. Restating it here would drift from the code under test.
    """
    from rigel.calibration.priors import _boundary_locus_shares

    arr = np.asarray(getattr(measured.calibration, field), np.float64)
    keep = _in_locus_regions(measured)
    if "boundary" in field:
        e_idx, _lid, _w = _boundary_locus_shares(
            measured.region_arrays, measured.multi_loci, len(measured.multi_loci)
        )
        keep = np.zeros(int(measured.calibration.n_boundaries), dtype=bool)
        keep[e_idx] = True
    if arr.shape[0] == keep.shape[0]:
        ranked = np.where(keep.reshape((-1,) + (1,) * (arr.ndim - 1)), arr, -np.inf)
    else:  # the sj axis — never projected, and never read by the prior
        ranked = arr
    site = np.unravel_index(int(np.argmax(ranked)), arr.shape)
    return site


def _fake_priors(gdna, rna, eff_len=None):
    """A ``LocusPriors`` with hand-chosen arrays — the scoring functions take the real type."""
    from rigel.calibration.priors import LocusPriors

    g = np.asarray(gdna, np.float64)
    return LocusPriors(
        gdna_prior_count=g,
        rna_prior_count=np.asarray(rna, np.float64),
        gdna_eff_len=np.ones_like(g) if eff_len is None else np.asarray(eff_len, np.float64),
    )


def _rebuild_calibration(m):
    """The condition's own ``CalibrationResult``. ⚠ Reconstructed from the SHIPPED masses the
    measurement kept rather than re-running ``calibrate``: the perturbations below only need an object
    whose six mass arrays are P's, and re-solving would take seconds per gate and could drift."""
    return m.calibration
