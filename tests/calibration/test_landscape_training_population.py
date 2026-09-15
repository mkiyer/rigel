"""The gDNA landscape prior's training population: which slots are allowed to train it.

A slot that holds no COMPOSITION — no own composition channel, no composition row received from a
neighbour, not structurally locked — is not in the training population. A level, a ceiling (a lane's,
or the node's own flux's) and a cube row are all BOUNDS: the value the solve settles on inside the
admitted half-line is the prior's own, so training the prior on it is training the prior on its echo
(`ISSUES: gdna-landscape-trains-on-false-positives`). `sweep.solve_chain` publishes the predicate as
`RegionBelief.has_composition`, read off the held messages, and `calibrate._fit_gdna_hyperprior` selects on
it AND on whether the solve LOCATES the slot: a posterior wider than one nat² in ``log f_g``
(`landscape._LOCATED_VAR`, the estimator's one-fragment wall through ``Var(log c) = 1/c``) has no
location whatever its evidence — a strand term at a pure-RNA vertex, an empty intron's factory row, a
one-sided delivered row — and its median is the reference's under its bound. The zero-count anchor
trains regardless, being a structural statement rather than a solve. "Any non-flat λ-row" is NOT the
predicate — `PsiMessage.lam_rows` fuses compositions and bounds together, so that reading keeps exactly
the bound-only slots this rule excludes.

PERTURBATION, each watched: with `has_composition` forced True everywhere, with the held compositions
dropped from the predicate, and with "any non-flat row" in place of the held composition, the
identity gates fail; with the selector ignoring the predicate, the population gates fail.
"""

from __future__ import annotations

import sys
from types import SimpleNamespace

import numpy as np
import pytest

import rigel.calibration.sweep as SW
from rigel.calibration.blocks import SweepCapture
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.region_chain import REGION
from rigel.calibration.region_geometry import g1_locked
from rigel.calibration.region_init import has_own_composition_evidence
from rigel.calibration.signature import RegionType
from _transfer_harness import _ctx_of, _full_policy, _passes, _prepared

CAL = sys.modules["rigel.calibration.calibrate"]


def _expected_has_composition(sweep_inputs, policy, capture):
    """The predicate re-derived INDEPENDENTLY of the solve: the two passes driven here on the same
    prepared policy, a slot has a composition iff either held message carries a COMPOSITION, or its own channel
    is live, or it is structurally certain. A level, a ceiling or a cube row does not count."""
    ctx = _ctx_of(sweep_inputs)
    prepared = _prepared(policy, ctx)
    from_left, from_right = _passes(prepared, ctx)
    comp = from_left.has_composition | from_right.has_composition
    tau = np.asarray(capture.tau_lam, np.float64)
    return has_own_composition_evidence(tau) | g1_locked(capture.free_pos, capture.free_neg) | comp


def test_the_solve_publishes_the_informed_predicate_as_the_solve_used_it(sweep_inputs):
    policy, *_ = _full_policy(sweep_inputs)
    capture = SweepCapture()
    out = SW.solve_chain(
        *sweep_inputs["args"], **sweep_inputs["kw"], policy=policy, _capture=capture
    )
    has_composition = np.asarray(out.has_composition, bool)
    assert has_composition.shape == (int(sweep_inputs["args"][0].n_slots),)
    assert np.array_equal(has_composition, _expected_has_composition(sweep_inputs, policy, capture))
    # not vacuous: the toy has both kinds
    assert has_composition.any() and (~has_composition).any()


def test_a_bound_with_a_row_does_not_inform_but_a_composition_does(sweep_inputs):
    """THE DISTINCTION THE RULING TURNS ON, on two blind slots (no own channel, not locked): a stub
    policy delivers a LEVEL to one and a COMPOSITION to the other, and its solve writes a non-flat row
    at BOTH — so "any non-flat row" would give both a composition. Only the composition's slot may."""
    from rigel.calibration.messages import PsiMessage

    cap = SweepCapture()
    SW.solve_chain(*sweep_inputs["args"], **sweep_inputs["kw"], policy=SilentPolicy(), _capture=cap)
    blind = np.flatnonzero(
        ~has_own_composition_evidence(cap.tau_lam)
        & ~g1_locked(cap.free_pos, cap.free_neg)
        & (np.asarray(cap.left, np.int64) >= 0)
    )
    assert blind.size >= 2, "the toy has fewer than two blind slots with a left neighbour"
    lvl_slot, comp_slot = int(blind[0]), int(blind[1])
    K = int(sweep_inputs["kw"]["n_grid"])
    lam = np.linspace(-1.0, 1.0, K)
    row = -0.5 * lam**2

    class _Prepared:
        def __init__(self, n):
            self.n = n

        def propagate(self, received, *, backward):
            if backward:
                return None

            def receive(s, i):
                if i == lvl_slot:
                    received.level_gdna.write(i, row, 3.0, 100.0)
                if i == comp_slot:
                    received.composition[i] = row
                    received.has_composition[i] = True

            return receive

        def solve(self, from_left, from_right):
            rows = np.zeros((self.n, K))
            rows[lvl_slot] = row
            rows[comp_slot] = row
            return PsiMessage(lam_rows=rows)

    class _Stub:
        name = "bound-vs-composition-stub"

        def library(self, view):
            return None

        def prepare(self, ctx, library):
            return _Prepared(int(ctx.n_slots))

    out = SW.solve_chain(*sweep_inputs["args"], **sweep_inputs["kw"], policy=_Stub())
    has_composition = np.asarray(out.has_composition, bool)
    assert has_composition[comp_slot], "a received composition must inform"
    assert not has_composition[lvl_slot], "a level with a row is a bound only and must not inform"


def test_silence_shrinks_the_informed_set_to_own_evidence_and_certainty(sweep_inputs):
    """PERTURBATION: with no messages, the delivered rows vanish and the ``has_composition`` set must shrink to
    the own channel plus structural certainty — and at least one slot must change, or the message
    layer's contribution to the predicate is not being read."""
    policy, *_ = _full_policy(sweep_inputs)
    live = np.asarray(
        SW.solve_chain(*sweep_inputs["args"], **sweep_inputs["kw"], policy=policy).has_composition,
        bool,
    )
    cap = SweepCapture()
    silent = SW.solve_chain(
        *sweep_inputs["args"], **sweep_inputs["kw"], policy=SilentPolicy(), _capture=cap
    )
    silent_has_composition = np.asarray(silent.has_composition, bool)
    own = has_own_composition_evidence(cap.tau_lam) | g1_locked(cap.free_pos, cap.free_neg)
    assert np.array_equal(silent_has_composition, own)
    assert (live & ~silent_has_composition).any(), "no slot was has_composition by a message alone"
    assert not (silent_has_composition & ~live).any()


def _synthetic_population():
    """Seven chain slots: R0 empty non-exon (the anchor), B1, R2 expressed single-strand exon, B3,
    R4 expressed locked (intergenic), B5, R6 expressed single-strand intron."""
    kind = np.array([REGION, 1, REGION, 1, REGION, 1, REGION], np.int8)
    obj = np.array([0, 0, 1, 1, 2, 2, 3], np.int64)
    chain = SimpleNamespace(kind=kind, obj_idx=obj, n_slots=7)
    fp = np.array([0, 0, 1, 1, 0, 0, 1], bool)
    fn = np.array([0, 0, 0, 0, 0, 0, 0], bool)
    statics = SimpleNamespace(free_pos=fp, free_neg=fn)
    from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS, coarse_type_array

    # coarse types (intergenic, exon, intergenic, intron) — the exon is the one the anchor excludes
    sig = np.array([0, BIT_EXON_POS, 0, BIT_INTRON_POS], np.int64)
    assert coarse_type_array(sig).tolist() == [
        RegionType.INTERGENIC,
        RegionType.EXON,
        RegionType.INTERGENIC,
        RegionType.INTRON,
    ]
    region_arrays = SimpleNamespace(signature=sig)
    mass = np.array([0.0, 5.0, 100.0, 5.0, 80.0, 5.0, 60.0])
    eff = np.full(7, 1000.0)
    belief = SimpleNamespace(
        f_g=np.array([0.0, 0.5, 0.3, 0.5, 1.0, 0.5, 0.9]),
        var_gdna=np.array([np.inf, 1.0, 0.8, 1.0, 0.0, 1.0, 0.5]),
        has_composition=np.array([False, True, True, True, True, True, True]),
    )
    return chain, belief, statics, region_arrays, mass, eff


def _training_counts(monkeypatch, belief, parts):
    chain, _, statics, region_arrays, mass, eff = parts
    seen = {}

    def spy(count, mass_, eff_, var, *, anchor, knn_scale=0.5, domain=None, prev=None):
        seen["count"] = np.asarray(count, np.float64)
        seen["anchor"] = np.asarray(anchor, bool)
        seen["domain_mass"] = None if domain is None else np.asarray(domain[0], np.float64)
        return None

    monkeypatch.setattr(CAL, "fit_landscape", spy)
    monkeypatch.setattr(CAL, "_MIN_TRAIN", 1)
    CAL._fit_gdna_hyperprior(chain, belief, statics, region_arrays, mass, eff)
    return seen


def test_a_flat_likelihood_slot_is_not_in_the_training_population(monkeypatch):
    parts = _synthetic_population()
    chain, belief, *_ = parts
    full = _training_counts(monkeypatch, belief, parts)
    # all three expressed regions plus the anchor train when every slot has a composition
    assert full["count"].shape == (4,) and full["anchor"].sum() == 1
    # PERTURBATION: the exon's likelihood was flat -> it leaves the population; nothing else moves
    blind = SimpleNamespace(
        f_g=belief.f_g,
        var_gdna=belief.var_gdna,
        has_composition=belief.has_composition & ~(np.arange(7) == 2),
    )
    part = _training_counts(monkeypatch, blind, parts)
    assert part["count"].shape == (3,)
    assert np.array_equal(part["count"], np.array([0.0, 80.0, 54.0]))
    assert part["anchor"].sum() == 1
    # THE GRID IS THE CONSUMERS' DOMAIN: every slot with opportunity, boundaries included, reaches the
    # estimator as ``domain`` — the exon that left the training set is still on the axis the prior is
    # read at, and so is every boundary
    assert np.array_equal(part["domain_mass"], np.array([0.0, 5.0, 100.0, 5.0, 80.0, 5.0, 60.0]))


def test_a_slot_wider_than_one_nat_does_not_train_whatever_its_evidence(monkeypatch):
    """THE LOCATION FLOOR: with a composition and a solve narrower than a nat² the exon trains; with the
    same composition and a solve wider than that it leaves the population (a pure-RNA vertex read at its
    resolution, the poison of the gDNA-free stranded rows: 3,771 false fragments trained at the first
    refit); at exactly one nat² it stays, the floor being the count rule's inclusive wall; and a slot
    with a narrow solve but no composition still does not train — the no-echo rule and the floor are a
    conjunction, since a bound-only slot sharpened by the prior alone is the prior's echo. The intron
    (var 0.5) and the locked region (var 0) are untouched throughout."""
    parts = _synthetic_population()
    chain, belief, *_ = parts
    full = _training_counts(monkeypatch, belief, parts)
    assert full["count"].shape == (4,), "the exon at var 0.8 trains"

    def with_var(v, has=None):
        var = belief.var_gdna.copy()
        var[2] = v
        return SimpleNamespace(
            f_g=belief.f_g,
            var_gdna=var,
            has_composition=belief.has_composition if has is None else has,
        )

    wide = _training_counts(monkeypatch, with_var(2.0), parts)
    assert wide["count"].shape == (3,) and np.array_equal(
        wide["count"], np.array([0.0, 80.0, 54.0])
    )
    edge = _training_counts(monkeypatch, with_var(1.0), parts)
    assert edge["count"].shape == (4,), "the floor is inclusive: one nat² is the wall itself"
    echo = _training_counts(
        monkeypatch, with_var(0.3, belief.has_composition & ~(np.arange(7) == 2)), parts
    )
    assert echo["count"].shape == (3,), "a bound-only slot does not train however narrow its solve"


def test_the_substrate_guard_measures_the_domain_not_the_cut(monkeypatch):
    """With the cut leaving only the anchor, the prior still fits as long as the DOMAIN has enough
    substrate — otherwise a gDNA-free toy loses its refit and invents gDNA in its place."""
    parts = _synthetic_population()
    _, belief, *_ = parts
    monkeypatch.setattr(
        CAL, "_MIN_TRAIN", 4
    )  # the domain (4) passes, the cut population (1) would not
    none = SimpleNamespace(
        f_g=belief.f_g, var_gdna=belief.var_gdna, has_composition=np.zeros(7, bool)
    )
    chain, _, statics, region_arrays, mass, eff = parts
    calls = []
    monkeypatch.setattr(CAL, "fit_landscape", lambda *a, **k: calls.append(1))
    CAL._fit_gdna_hyperprior(chain, none, statics, region_arrays, mass, eff)
    assert calls, "the refit was refused although the domain has enough substrate"


def test_the_anchor_trains_whatever_the_predicate_says(monkeypatch):
    parts = _synthetic_population()
    _, belief, *_ = parts
    none = SimpleNamespace(
        f_g=belief.f_g, var_gdna=belief.var_gdna, has_composition=np.zeros(7, bool)
    )
    seen = _training_counts(monkeypatch, none, parts)
    assert seen["count"].shape == (1,) and seen["anchor"].all() and seen["count"][0] == 0.0


def test_a_belief_without_the_predicate_trains_the_old_population(monkeypatch):
    """An initial belief (`init_beliefs`, no solve yet) carries no predicate; the selector then reads
    the annotation alone — the population before this rule — so nothing upstream of the first sweep
    changes shape."""
    parts = _synthetic_population()
    _, belief, *_ = parts
    old = SimpleNamespace(f_g=belief.f_g, var_gdna=belief.var_gdna, has_composition=None)
    seen = _training_counts(monkeypatch, old, parts)
    assert seen["count"].shape == (4,)


# ── THE E-STEP ON THE KERNELS THAT HAVE NO LOCATION ──────────────────────────────────────────────────
#
# A region trained at less than one fragment has no location of its own: the estimator already centres
# it at its resolution wall (``max(count, 1)``), and its Poisson kernel is flat below that wall. Summing
# such a kernel normalised to unit mass spreads that mass uniformly under the wall — so short empty
# regions would deposit prior mass at exon densities. The refit loop already holds the previous fit, and
# the deconvolution's E-step places a location-free kernel where the population is: kernel × previous
# landscape, renormalised. Counted kernels keep their own location, so an enriched minority cannot be
# competed away.


def _prev_at_floor(grid_log10):
    """A previous landscape concentrated at the grid's floor."""
    from rigel.calibration.landscape import DensityLandscape

    p = np.exp(-0.5 * ((grid_log10 - grid_log10[0]) / 0.1) ** 2) + 1e-12
    return DensityLandscape(
        log_rho=grid_log10 * np.log(10.0),
        logP=np.log(p / p.sum()),
        n_train=1,
        centre=np.array([grid_log10[0] * np.log(10.0)]),
        width=np.array([0.1 * np.log(10.0)]),
    )


def test_the_estep_moves_a_location_free_kernel_to_the_population_and_leaves_a_counted_one_alone():
    from rigel.calibration import landscape as LS

    count = np.array([0.0, 0.4, 30.0, 0.0])
    mass = np.array([0.0, 2.0, 100.0, 0.0])
    eff = np.array([50.0, 50.0, 1000.0, 20000.0])
    grid = LS._grid(mass, eff)
    base = LS._poisson_kernels(count, eff, grid)
    prev = _prev_at_floor(grid)
    out = LS._estep_kernels(base.copy(), count, grid, prev)

    def mean(k):
        return float((k * grid).sum() / k.sum())

    # the two location-free kernels (0 and 0.4 fragments) move DOWN toward the floor
    assert mean(out[0]) < mean(base[0]) - 0.5
    assert mean(out[1]) < mean(base[1]) - 0.5
    # the counted kernel is bit-identical
    assert np.array_equal(out[2], base[2])
    # every kernel still has unit mass
    assert np.allclose(out.sum(1), 1.0)
    # PERTURBATION: no previous landscape -> nothing moves
    assert np.array_equal(LS._estep_kernels(base.copy(), count, grid, None), base)


def test_fit_landscape_without_a_previous_fit_is_the_shipped_fit():
    from rigel.calibration import landscape as LS

    rng = np.random.default_rng(3)
    n = 40
    eff = np.exp(rng.uniform(np.log(50), np.log(20000), n))
    count = np.where(rng.random(n) < 0.5, 0.0, 1.0 + rng.poisson(3.0, n).astype(float))
    mass = np.where(count > 0, count + rng.poisson(5.0, n), 0.0)  # the zero-count rows are ANCHORS
    var = np.where(count > 0, 0.5, np.inf)
    anchor = mass <= 0.0
    assert anchor.sum() >= 10 and (~anchor).sum() >= 10
    a = LS.fit_landscape(count, mass, eff, var, anchor=anchor)
    b = LS.fit_landscape(count, mass, eff, var, anchor=anchor, prev=None)
    assert np.array_equal(a.logP, b.logP)
    # and WITH a previous fit at the floor the mass above the floor's neighbourhood falls
    c = LS.fit_landscape(
        count, mass, eff, var, anchor=anchor, prev=_prev_at_floor(a.log_rho / np.log(10.0))
    )
    g = a.log_rho / np.log(10.0)
    above = g > g[0] + 1.0
    assert np.exp(c.logP)[above].sum() < np.exp(a.logP)[above].sum()


def test_the_refit_loop_hands_each_fit_the_previous_landscape(sweep_inputs, monkeypatch):
    """The first fit sees no previous landscape; every later fit sees the one before it."""
    from rigel.config import CalibrationConfig

    seen = []
    orig = CAL.fit_landscape

    def spy(*a, **k):
        seen.append(k.get("prev"))
        return orig(*a, **k)

    monkeypatch.setattr(CAL, "fit_landscape", spy)
    CAL.calibrate(
        payload=sweep_inputs["payload"],
        config=CalibrationConfig(calib_refit_iters=3),
        **sweep_inputs["calibrate_kw"],
    )
    assert len(seen) == 3
    assert seen[0] is None
    assert seen[1] is not None and seen[2] is not None
    assert seen[2] is not seen[1]


def test_the_result_publishes_the_last_landscapes_located_enriched_mode(sweep_inputs, monkeypatch):
    """The reference the ruler and the prior assembler read is the located enriched mode of the LAST
    refit's landscape, published on the result; with no refit there is no landscape and no reference.
    PERTURBATION: with the located-mode reader forced to answer, the result carries exactly that answer."""
    from rigel.calibration.abundance_landscape import AbundanceMode
    from rigel.config import CalibrationConfig

    res0 = CAL.calibrate(
        payload=sweep_inputs["payload"],
        config=CalibrationConfig(calib_refit_iters=0),
        **sweep_inputs["calibrate_kw"],
    )
    assert res0.gdna_reference_density is None

    seen = []
    orig = CAL.located_enriched_mode

    def spy(ls):
        out = orig(ls)
        seen.append(out)
        return out

    monkeypatch.setattr(CAL, "located_enriched_mode", spy)
    res = CAL.calibrate(
        payload=sweep_inputs["payload"],
        config=CalibrationConfig(calib_refit_iters=2),
        **sweep_inputs["calibrate_kw"],
    )
    assert len(seen) == 1, "the reader runs once, on the last landscape"
    if seen[0] is None:
        assert res.gdna_reference_density is None
    else:
        assert res.gdna_reference_density == pytest.approx(float(np.exp(seen[0].log_rho)))

    forced = AbundanceMode(log_rho=-2.0, basin_mass=0.3, width=0.1, lo=-3.0, hi=-1.0)
    monkeypatch.setattr(CAL, "located_enriched_mode", lambda ls: forced)
    res2 = CAL.calibrate(
        payload=sweep_inputs["payload"],
        config=CalibrationConfig(calib_refit_iters=1),
        **sweep_inputs["calibrate_kw"],
    )
    assert res2.gdna_reference_density == pytest.approx(float(np.exp(-2.0)))
