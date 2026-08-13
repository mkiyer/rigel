"""Falsification gates for ``scripts/design/solvability_audit.py``.

⛔ **THE INSTRUMENT EXISTS BECAUSE THE OBVIOUS MEASUREMENT WAS WRONG.** Scoring every object that
carries mass counts an object correctly saying "I cannot be solved without a prior" as an error, and on
the hard condition that buried a 0.0456 answer inside a 0.3150 one. Every gate below guards one of the
ways that mistake, or a neighbouring one, can come back.

⛔ Each gate carries its own perturbation.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from rigel.calibration.region_chain import REGION
from rigel.config import CalibrationConfig, PipelineConfig
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        path = Path(__file__).resolve().parents[2] / "scripts" / "design" / name
        spec = importlib.util.spec_from_file_location(key, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


SA = _sibling("solvability_audit.py")
P0 = _sibling("pass0_vs_oracle.py")


@pytest.fixture(scope="module")
def audited(tmp_path_factory):
    wd = tmp_path_factory.mktemp("sa_sim")
    sc = Scenario("sa", genome_length=9000, seed=29, work_dir=wd / "sim")
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(600, 1100), (1800, 2300)], "abundance": 70},
            {"t_id": "t1b", "exons": [(1000, 1100), (1800, 1900)], "abundance": 25},
        ],
    )
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(4000, 4500), (5200, 5700)], "abundance": 35}])
    sc.add_gene("g3", "+", [{"t_id": "t3", "exons": [(6600, 7100), (7800, 8300)], "abundance": 20}])
    res = sc.build_oracle(
        n_rna_fragments=4000,
        gdna_fraction=1.0,
        nrna_abundance=15.0,
        sim_config=ReadSimConfig(
            frag_mean=180,
            frag_std=30,
            frag_min=80,
            frag_max=400,
            # ⭐ UNSTRANDED on purpose: at kappa = 1/2 the strand channel is exactly silent,
            # which is what CREATES the undetermined population this instrument is about.
            # A stranded toy has strand evidence everywhere and GATE 1 is vacuous on it.
            read_length=90,
            strand_specificity=0.50,
            seed=29,
        ),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=240, frag_std=45),
    )
    cfg = CalibrationConfig()
    m = P0.measure_condition(
        bam=str(res.bam_path),
        index=res.index,
        pipeline_config=PipelineConfig(),
        calibration_config=cfg,
        work_dir=tmp_path_factory.mktemp("sa_split"),
        tag="sa",
    )
    return m, SA.audit(m, axis="region", config=cfg), cfg


@pytest.fixture(scope="module")
def audited_stranded(tmp_path_factory):
    """⭐ The STRANDED twin. The two populations this instrument distinguishes need OPPOSITE
    libraries to exist at all: the undetermined class requires kappa = 1/2 (strand silent), and a
    slot carrying two channels at once requires strand to be live. One scenario cannot show both,
    and asserting both against one fixture is how a gate ends up vacuous."""
    wd = tmp_path_factory.mktemp("sa_sim_ss")
    sc = Scenario("sass", genome_length=9000, seed=31, work_dir=wd / "sim")
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(600, 1100), (1800, 2300)], "abundance": 70}])
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(4000, 4500), (5200, 5700)], "abundance": 35}])
    res = sc.build_oracle(
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
            seed=31,
        ),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=240, frag_std=45),
    )
    cfg = CalibrationConfig()
    m = P0.measure_condition(
        bam=str(res.bam_path),
        index=res.index,
        pipeline_config=PipelineConfig(),
        calibration_config=cfg,
        work_dir=tmp_path_factory.mktemp("sa_split_ss"),
        tag="sass",
    )
    return m, SA.audit(m, axis="region", config=cfg), cfg


# ── GATE 1: the undetermined population is EXCLUDED, and it is not empty ─────────────────────────


def test_the_UNDETERMINED_population_is_excluded_and_is_NOT_EMPTY(audited):
    """⛔ The whole point. Undetermined and determined must partition the live objects with no overlap
    and no gap, and the undetermined set must be non-empty or the instrument is measuring nothing that
    the naive all-objects score did not already measure.

    PERTURBATION: fold the undetermined objects back in and require the scored error to CHANGE — if it
    does not, this scenario cannot demonstrate the distinction and the gate is vacuous.
    """
    m, a, _ = audited
    det, und, live = a["determined"], a["undetermined"], a["live"]
    assert not (det & und).any(), "an object is both determined and undetermined"
    np.testing.assert_array_equal(det | und, live)
    assert und.sum() > 0, "no undetermined objects; the gate would be vacuous"
    assert det.sum() > 0, "no determined objects; nothing would be scored"

    err = np.abs(a["err"])
    assert err[live].sum() != pytest.approx(err[det].sum(), rel=1e-9), (
        "including the undetermined population changes nothing here, so this scenario cannot "
        "demonstrate the distinction the instrument exists to make"
    )


# ── GATE 2: the channels partition, and 'none' is exactly the undetermined set ────────────────────


def test_the_own_evidence_CHANNELS_partition_and_match_the_undetermined_set(
    audited, audited_stranded
):
    """Solvability is "some own-evidence channel can speak here", so ``none`` must be exactly the
    undetermined set and the union of the rest exactly the determined set.

    ⚠ The channels themselves OVERLAP and must not be asserted to partition — ``tau_lam`` is the SUM
    of the strand and factory arms, so a single-stranded intron region has both. A first draft of this
    gate asserted a partition and failed on exactly those regions; the instrument was corrected to
    describe them as overlapping capabilities and to warn that per-channel mass double-counts.

    PERTURBATION: require the overlap to be NON-EMPTY, so the gate fails if the two definitions ever
    silently become exclusive again.
    """
    m, a, cfg = audited
    ch = a["channels"]
    live = a["live"]
    assert not (ch["none"] & (ch["locked"] | ch["strand"] | ch["factory"])).any()
    np.testing.assert_array_equal(ch["none"] & live, a["undetermined"])
    np.testing.assert_array_equal(
        (ch["locked"] | ch["strand"] | ch["factory"]) & live, a["determined"]
    )
    # ⭐ the overlap itself needs the STRANDED twin: at kappa = 1/2 strand is silent everywhere, so
    # no slot on `audited` can carry two channels and the assertion would be vacuous there.
    sch = audited_stranded[1]["channels"]
    assert (sch["strand"] & sch["factory"]).sum() > 0, (
        "no slot carries both strand and factory even on a stranded library — the two definitions "
        "have silently become exclusive"
    )


# ── GATE 3: the standardised discrepancy is in LOG space ──────────────────────────────────────────


def test_the_discrepancy_is_LOG_space_because_var_gdna_IS(audited):
    """⛔ ``var_gdna`` is ``Var(log f_g)`` while its name says "fraction variance". A linear error over
    a log-space variance is the defect that once moved a suite total 0.046 → 1.007 and inverted a
    per-class ranking, so this gate pins the space explicitly.

    PERTURBATION: the LINEAR form must give a materially different answer — if the two agreed, the
    distinction would not matter and this gate would be decoration.
    """
    grid = np.array([0.01, 0.5, 0.99])
    pred = np.array([0.10, 0.20, 0.80])
    true = np.array([0.20, 0.20, 0.40])
    var = np.array([0.25, 0.25, 0.25])  # sd = 0.5 in LOG space

    z, gap, sd = SA.standardised_discrepancy(pred, true, var, grid)
    np.testing.assert_allclose(gap, np.log(pred) - np.log(true))
    np.testing.assert_allclose(z, (np.log(pred) - np.log(true)) / 0.5)

    linear_z = (pred - true) / 0.5
    assert not np.allclose(z, linear_z), "log and linear agree here; the gate proves nothing"


def test_truth_is_clipped_to_the_SOLVERS_OWN_GRID_not_an_epsilon(audited):
    """⛔ Truth is exactly 0 at a pure-RNA object and exactly 1 at a structurally-pure-gDNA one — both
    common — and ``log 0`` is not a number. The clip must be the solver's own λ-grid support, because
    that is the best answer the solver could ever have given.

    PERTURBATION: a different grid must move the answer, proving the grid is actually consulted rather
    than a hard-coded epsilon standing in for it.
    """
    narrow = np.array([0.05, 0.5, 0.95])
    wide = np.array([1e-4, 0.5, 1 - 1e-4])
    pred = np.array([0.5, 0.5])
    true = np.array([0.0, 1.0])  # the two degenerate cases
    var = np.array([1.0, 1.0])

    z_n, gap_n, _ = SA.standardised_discrepancy(pred, true, var, narrow)
    z_w, gap_w, _ = SA.standardised_discrepancy(pred, true, var, wide)
    assert np.all(np.isfinite(gap_n)) and np.all(np.isfinite(gap_w)), "log 0 leaked through"
    assert not np.allclose(gap_n, gap_w), "the grid is ignored; a fixed epsilon is being used"
    assert np.abs(gap_w).sum() > np.abs(gap_n).sum(), "a wider grid must admit a larger discrepancy"


def test_a_CERTAIN_slot_that_is_wrong_is_infinitely_confident_not_a_division_by_zero(audited):
    """⭐ ``sd == 0`` is the solver calling an object CERTAIN. If it is also wrong that is the worst
    possible case, and it must sort to the top rather than crash or silently become 0.

    PERTURBATION: certain AND right must be 0, not inf — otherwise every locked object would be
    flagged.
    """
    grid = np.array([0.01, 0.5, 0.99])
    z, _, _ = SA.standardised_discrepancy(
        np.array([0.9, 0.5]), np.array([0.1, 0.5]), np.array([0.0, 0.0]), grid
    )
    assert np.isinf(z[0]), "a certain-and-wrong slot must be infinitely confident"
    assert z[1] == 0.0, "a certain-and-right slot must be 0, not inf"


# ── GATE 4: the ablation ladder is the solver's own arithmetic ────────────────────────────────────


def test_the_ladders_FINAL_rung_IS_the_arm_it_claims_to_be(audited):
    """⭐ The ladder's conclusions are differences between rungs, so a projection error would invent a
    channel effect out of nothing. The FINAL rung must be byte-equal to the arm's own per-object
    fraction — that is the check that the whole projection is right.

    PERTURBATION: the rungs must not all be identical, or the ladder cannot attribute anything.
    """
    m, a, _ = audited
    g = np.asarray(m.arms["pass0"].mass_gdna_region, np.float64)
    r = np.asarray(m.arms["pass0"].mass_rna_region, np.float64)
    frac, _ = P0.object_fractions(g, r)
    live = a["live"]
    np.testing.assert_allclose(a["ladder"]["f_g"][live], frac[live], rtol=0, atol=1e-12)

    rungs = [a["ladder"][k][live] for k in ("fg_strand", "fg_loc", "f_g")]
    assert not np.allclose(rungs[0], rungs[2]), "strand and final coincide; the ladder is flat"


# ── GATE 5: the z bands and the report cover every scored object ──────────────────────────────────


# ── GATE 6: own evidence must RESOLVE something on the solver's own λ grid ────────────────────────


def test_own_evidence_STRENGTH_is_reported_as_a_CURVE_because_tau_is_CONTINUOUS(audited):
    """⛔⛔ **A precision above zero is not evidence, and no threshold can say where evidence starts.**

    ``tau_lam`` is a Fisher precision on ``λ = log(f_g/f_R)``, so an object's own statement has sd
    ``1/√τ`` nats against a solver that represents only ``λ ∈ [−L, +L]``. The strand arm carries
    ``I(f_g) ∝ (2κ−1)²`` (`EQUATIONS.md` §5.2), **exactly zero at κ = ½** — but κ is FITTED, so on a
    genuinely unstranded library it misses ½ by a few 1e-4 and τ lands at ~1e-7 rather than at 0. The
    ``τ > 1e-9`` cut then scores that object as SOLVABLE while its own statement has sd ~10³ nats.

    ⛔ **AND A BETTER THRESHOLD IS NOT THE FIX — that was tried and REFUTED.** A resolving-power floor
    at ``1/(2L)²`` was derived and implemented, and its own insensitivity gate killed it: τ is
    CONTINUOUS across that region on 4 of 5 ladder conditions, so no interval is empty and any floor
    is a tuned constant. Only the unstranded capture-OFF row is bimodal, and there the two clusters
    are the silent strand arm (~1e-7) and the live intron factory (~1e-1).

    ⭐ So the instrument reports the CURVE. This gate pins that it covers the population exactly and
    resolves the two extremes into different rows — which is all a curve has to do.

    PERTURBATION: a single-decade band must NOT reproduce the curve's discrimination, proving the
    decades are load-bearing rather than decoration.
    """
    m, a, cfg = audited
    # ⛔ BEHAVIOURAL, not by name: the channel test must BE the solver's own gate, so that "has a
    # channel" and "the evidence is strong enough to score" stay two separate questions. A threshold
    # smuggled in here would silently re-partition every stratum table in the project.
    cap, chain = m.debug_pass0["capture"], m.debug_pass0["chain"]
    tau = np.asarray(cap["_tau0_lam"], np.float64)
    ch = SA.channel_masks(cap, chain, cfg)
    np.testing.assert_array_equal(
        ch["strand"] | ch["factory"],
        (tau > SA._EPS) & ~ch["locked"],
        err_msg="channel_masks no longer asks the solver's own question; a strength threshold has "
        "been introduced where a curve belongs",
    )
    # ⚠ …and on a CONSTRUCTED τ that actually spans the region, because the fixture's own-evidence
    # slots all sit at one strength and a floor placed below them would be invisible here.
    probe_tau = np.array([0.0, 1e-8, 1e-7, 1e-5, 1e-3, 2e-3, 1e-2, 1e-1, 1.0])
    probe = {
        "_tau0_lam": probe_tau,
        "free_pos": np.ones(probe_tau.size, bool),
        "free_neg": np.zeros(probe_tau.size, bool),
        "intron_prior": None,
    }
    pchain = type("C", (), {"n_slots": probe_tau.size})()
    pch = SA.channel_masks(probe, pchain, cfg)
    np.testing.assert_array_equal(
        pch["strand"] | pch["factory"],
        probe_tau > SA._EPS,
        err_msg="a strength threshold has been reintroduced in channel_masks — it was REFUTED "
        "because tau is continuous across this region, so any floor is a tuned constant",
    )

    # ⚠ and the capture must be checked against the chain it claims to describe: a mismatch would
    # shift every mask by one slot, which is invisible in aggregate.
    with pytest.raises(ValueError, match="different partitions"):
        SA.channel_masks(probe, type("C", (), {"n_slots": probe_tau.size + 1})(), cfg)

    det = a["determined"] & ~a["channels"]["locked"]
    rows = SA.resolving_power_rows(a, det)
    assert sum(n for _, n, _, _, _, _ in rows) == int(det.sum()), "the curve loses objects"
    np.testing.assert_allclose(
        sum(mass for _, _, mass, _, _, _ in rows), float(a["total"][det].sum()), rtol=1e-9
    )
    np.testing.assert_allclose(
        sum(e for _, _, _, e, _, _ in rows), float(np.abs(a["err"])[det].sum()), rtol=1e-9
    )
    assert SA.SD_LAMBDA_DECADES[-1] == np.inf, (
        "the top band must be unbounded or the objects with NO usable own evidence are dropped — "
        "which is the entire population the curve exists to make visible"
    )

    # ⭐ THE DISCRIMINATION is a property of the banding, so it is tested on a CONSTRUCTED sd(λ) that
    # spans the two populations the ladder shows. A 9 kb toy has three own-evidence slots, all at the
    # same strength, so asserting separation on the fixture would only prove the toy is small.
    synth = {
        "sd_lam": np.array([0.5, 3.0, 30.0, 300.0, 3000.0, np.inf]),
        "err": np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0]),
        "total": np.ones(6) * 10.0,
        "f_pred": np.full(6, 0.5),
        "f_true": np.full(6, 1.0),
    }
    mask = np.ones(6, bool)
    srows = [r for r in SA.resolving_power_rows(synth, mask) if r[1]]
    assert len(srows) == 5, f"the decades do not separate five orders of magnitude: {srows}"
    assert [r[1] for r in srows] == [1, 1, 1, 1, 2], "the top band must absorb 3000 AND inf"
    assert sum(r[3] for r in srows) == pytest.approx(63.0), "the error does not decompose"

    # PERTURBATION: collapse to one band and the discrimination must vanish.
    saved = SA.SD_LAMBDA_DECADES
    try:
        SA.SD_LAMBDA_DECADES = (np.inf,)
        flat = [r for r in SA.resolving_power_rows(synth, mask) if r[1]]
        assert len(flat) == 1 and flat[0][1] == 6, "the collapse did not collapse"
    finally:
        SA.SD_LAMBDA_DECADES = saved


def test_the_sd_LAMBDA_is_the_solvers_own_tau_and_locked_slots_are_CERTAIN(audited):
    """``sd_lam`` must be ``1/√τ`` from the solver's own ``_tau0_lam`` — not a re-derivation — with
    ``0`` at a structurally-locked slot (certain, nothing to resolve) and ``inf`` where no channel
    spoke at all. Those two encodings are what make the curve's ends meaningful.

    PERTURBATION: the mapping must be strictly monotone-decreasing in τ, so a stronger channel cannot
    land in a weaker band.
    """
    m, a, cfg = audited
    cap, chain = m.debug_pass0["capture"], m.debug_pass0["chain"]
    tau = np.asarray(cap["_tau0_lam"], np.float64)
    ch = SA.channel_masks(cap, chain, cfg)
    kind = np.asarray(chain.kind)
    axis_sel = kind == REGION  # `audited` is the region axis
    obj = np.asarray(chain.obj_idx, np.int64)

    sd = a["sd_lam"]
    for slot in np.flatnonzero(axis_sel):
        i = obj[slot]
        if ch["locked"][slot]:
            assert sd[i] == 0.0, "a structurally-locked slot must read as certain, not as uncertain"
        elif tau[slot] <= 0.0:
            assert np.isinf(sd[i]), "a slot with no channel must read as infinitely uncertain"
        else:
            assert sd[i] == pytest.approx(1.0 / np.sqrt(tau[slot]), rel=1e-9)

    # PERTURBATION: monotonicity. Two live slots, the stronger must have the smaller sd.
    live = axis_sel & ~ch["locked"] & (tau > 0)
    if live.sum() >= 2:
        idx = np.flatnonzero(live)
        order = np.argsort(tau[idx])
        s = sd[obj[idx[order]]]
        assert np.all(np.diff(s) <= 1e-12), "sd(λ) is not monotone-decreasing in τ"


def test_a_STRUCTURALLY_LOCKED_EDGE_is_as_DETERMINED_as_a_locked_region(audited):
    """⛔ ``locked`` was ``~solvable & (kind == REGION)``, so a G1 **edge** — an intergenic↔exon seam,
    where RNA cannot cross a gene boundary and the solver pins ``{0,0,1}`` at ``Var(log f_g) = 0`` —
    fell into ``none`` and was EXCLUDED from the scored population as "honest ignorance". It is the
    opposite: it is structurally certain, and correct. Structural certainty is a property of the
    signature, not of which axis the object lives on.

    PERTURBATION: the fixture must actually contain locked edges, and folding them in must move the
    scored mass — otherwise the gate is vacuous.
    """
    m, a, cfg = audited
    cap, chain = m.debug_pass0["capture"], m.debug_pass0["chain"]
    kind = np.asarray(chain.kind)
    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    g1 = ~fp & ~fn

    ch = SA.channel_masks(cap, chain, cfg)
    locked_edges = ch["locked"] & (kind != REGION)
    assert (g1 & (kind != REGION)).sum() > 0, "no G1 edges in the fixture; the gate is vacuous"
    assert locked_edges.sum() > 0, "a G1 edge is not classed locked"
    np.testing.assert_array_equal(ch["locked"], g1)
    assert not (ch["locked"] & ch["none"]).any(), "a locked slot is also 'none'"

    # PERTURBATION: on the EDGE axis those objects must carry mass, or excluding them cost nothing.
    edge_audit = SA.audit(m, axis="edge", config=cfg)
    assert edge_audit["channels"]["locked"].sum() > 0, "no locked objects on the edge axis"
    assert edge_audit["total"][edge_audit["channels"]["locked"] & edge_audit["live"]].sum() > 0, (
        "the locked edges carry no mass, so their misclassification cost nothing and this gate "
        "cannot demonstrate the defect"
    )


def test_the_UNDETERMINED_class_is_checked_for_the_ONE_thing_it_can_get_wrong(audited):
    """⛔⛔ The undetermined population is excluded from every error total — correctly, because
    ``f_g ≈ ½`` at zero precision is a true statement about itself. But `SUCCESS.md` names the failure
    mode that exclusion leaves open ("claiming a precision it has not earned") and **nothing was
    checking it**, so a 395,251-fragment systematic error on ``gdna_g25_ss_0.50_nrna_none_capture_off``
    — 87 exon regions driven to ``f_g = 0.829`` against a truth of 0.009 — was invisible, at 0.0 %
    scored, on a condition whose reported mwae is 0.0170.

    ⭐ The check needs no threshold: the class's correct answer is ½, so bucket it by ``|f_pred − ½|``
    and report the error and the claimed precision per bucket.

    PERTURBATION: the buckets must cover the class exactly and must SEPARATE a moved object from an
    unmoved one — a single bucket would hide exactly what this exists to show.
    """
    m, a, cfg = audited
    rows = SA.undetermined_overreach_rows(a)
    und = a["undetermined"]
    assert sum(n for _, n, _, _, _, _ in rows) == int(und.sum()), "the buckets lose objects"
    np.testing.assert_allclose(
        sum(mass for _, _, mass, _, _, _ in rows), float(a["total"][und].sum()), rtol=1e-9
    )
    np.testing.assert_allclose(
        sum(e for _, _, _, e, _, _ in rows), float(np.abs(a["err"])[und].sum()), rtol=1e-9
    )
    assert SA.OVERREACH_BANDS[-1] == 0.5, (
        "the top band must reach ½ — |f_pred − ½| cannot exceed it, so a lower ceiling drops the "
        "objects moved furthest, which are the whole point"
    )

    # ⭐ separation, on a CONSTRUCTED class spanning all four bands: ½ at sd = ∞ is the CORRECT
    # answer; 0.83 with a FINITE sd is the defect. They must not share a row.
    synth = {
        "undetermined": np.ones(5, bool),
        "f_pred": np.array([0.500, 0.600, 0.700, 0.830, 0.980]),  # |f−½| = 0, .10, .20, .33, .48
        "total": np.full(5, 10.0),
        "err": np.array([0.0, 1.0, 4.0, 8.0, 9.0]),
        "sd": np.array([np.inf, np.inf, 0.4, 0.4, 0.4]),
    }
    srows = [r for r in SA.undetermined_overreach_rows(synth) if r[1]]
    assert len(srows) == 4, f"the four bands do not separate five spread values: {srows}"
    assert srows[0][1] == 1 and srows[0][5] == pytest.approx(0.0), (
        "the object correctly AT ½ with sd = inf must sit alone in the lowest band and must not "
        "count as claiming precision"
    )
    assert srows[-1][5] == pytest.approx(1.0), "a finite sd must count as claiming precision"
    assert srows[-1][3] > srows[0][3], "the error must concentrate in the moved bands"

    # PERTURBATION: one bucket hides the defect entirely — the ½ object and the 0.83 object merge.
    saved = SA.OVERREACH_BANDS
    try:
        SA.OVERREACH_BANDS = (0.5,)
        flat = [r for r in SA.undetermined_overreach_rows(synth) if r[1]]
        assert len(flat) == 1 and flat[0][1] == 5, "the collapse did not collapse"
        assert 0.0 < flat[0][5] < 1.0, "the collapsed row averages away the precision claim"
    finally:
        SA.OVERREACH_BANDS = saved


def test_the_z_BANDS_account_for_every_solvable_object(audited):
    """A confidence profile that silently drops objects — NaN z, infinite z — would understate the
    confidently-wrong class, which is the one number the instrument exists to produce.

    PERTURBATION: bands that stopped short of infinity would drop the certain-and-wrong objects; the
    band set must therefore reach ``inf`` and the counts must sum to the whole population.
    """
    m, a, _ = audited
    det = a["determined"]
    rows = SA._band_table(a, det)
    assert sum(n for _, n, _, _ in rows) == int(det.sum()), "the bands lose objects"
    np.testing.assert_allclose(
        sum(mass for _, _, mass, _ in rows), float(a["total"][det].sum()), rtol=1e-9
    )
    assert SA.Z_BANDS[-1][1] == np.inf, "the top band must be unbounded or certain-wrong is dropped"
