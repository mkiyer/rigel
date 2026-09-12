"""The gDNA landscape prior's training population: which slots are allowed to train it.

A slot that holds no COMPOSITION — no own composition channel, no composition row received from a
neighbour, not structurally locked — is not in the training population. A level, a ceiling (a lane's,
or the node's own flux's) and a cube row are all BOUNDS: the value the solve settles on inside the
admitted half-line is the prior's own, so training the prior on it is training the prior on its echo
(`ISSUES: gdna-landscape-trains-on-false-positives`). `sweep.solve_chain` publishes the predicate as
`RegionBelief.informed`, read off the held messages, and `calibrate._fit_gdna_hyperprior` selects on
it; the zero-count anchor trains regardless, being a structural statement rather than a solve. "Any
non-flat λ-row" is NOT the predicate — `PsiMessage.lam_rows` fuses compositions and bounds together,
so that reading keeps exactly the bound-only slots this rule excludes.

PERTURBATION, each watched: with `informed` forced True everywhere, with the held compositions
dropped from the predicate, and with "any non-flat row" in place of the held composition, the
identity gates fail; with the selector ignoring the predicate, the population gates fail.
"""

from __future__ import annotations

import sys
from types import SimpleNamespace

import numpy as np

import rigel.calibration.sweep as SW
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.region_chain import REGION
from rigel.calibration.region_geometry import g1_locked
from rigel.calibration.region_init import has_own_composition_evidence
from rigel.calibration.signature import RegionType
from _transfer_harness import _ctx_of, _full_policy, _prepared

CAL = sys.modules["rigel.calibration.calibrate"]


def _expected_informed(sweep_inputs, policy, capture):
    """The predicate re-derived INDEPENDENTLY of the solve: the two passes driven here on the same
    prepared policy, a slot informed iff either held message carries a COMPOSITION, or its own channel
    is live, or it is structurally certain. A level, a ceiling or a cube row does not count."""
    from rigel.calibration.messages import SILENCE

    ctx = _ctx_of(sweep_inputs)
    prepared = _prepared(policy, ctx)
    n = int(ctx.n_slots)
    comp = np.zeros(n, bool)
    order = list(range(n))
    for nbr, seq, backward in (
        (np.asarray(ctx.left, np.int64), order, False),
        (np.asarray(ctx.right, np.int64), order[::-1], True),
    ):
        receive = prepared.propagate(backward=backward)
        for i in seq:
            src = int(nbr[i])
            if src < 0:
                continue
            m = SILENCE if receive is None else receive(src, i)
            if m is not None and m.composition is not None:
                comp[i] = True
    tau = np.asarray(capture["_tau0_lam"], np.float64)
    return (
        has_own_composition_evidence(tau)
        | g1_locked(capture["free_pos"], capture["free_neg"])
        | comp
    )


def test_the_solve_publishes_the_informed_predicate_as_the_solve_used_it(sweep_inputs):
    policy, *_ = _full_policy(sweep_inputs)
    capture: dict = {}
    out = SW.solve_chain(
        *sweep_inputs["args"], **sweep_inputs["kw"], policy=policy, _capture=capture
    )
    informed = np.asarray(out.informed, bool)
    assert informed.shape == (int(sweep_inputs["args"][0].n_slots),)
    assert np.array_equal(informed, _expected_informed(sweep_inputs, policy, capture))
    # not vacuous: the toy has both kinds
    assert informed.any() and (~informed).any()


def test_a_bound_with_a_row_does_not_inform_but_a_composition_does(sweep_inputs):
    """THE DISTINCTION THE RULING TURNS ON, on two blind slots (no own channel, not locked): a stub
    policy delivers a LEVEL to one and a COMPOSITION to the other, and its solve writes a non-flat row
    at BOTH — so "any non-flat row" would call both informed. Only the composition's slot may."""
    from rigel.calibration.messages import Level, Message, PsiMessage

    cap: dict = {}
    SW.solve_chain(*sweep_inputs["args"], **sweep_inputs["kw"], policy=SilentPolicy(), _capture=cap)
    blind = np.flatnonzero(
        ~has_own_composition_evidence(cap["_tau0_lam"])
        & ~g1_locked(cap["free_pos"], cap["free_neg"])
        & (np.asarray(cap["left"], np.int64) >= 0)
    )
    assert blind.size >= 2, "the toy has fewer than two blind slots with a left neighbour"
    lvl_slot, comp_slot = int(blind[0]), int(blind[1])
    K = int(sweep_inputs["kw"]["n_grid"])
    lam = np.linspace(-1.0, 1.0, K)
    row = -0.5 * lam**2

    class _Prepared:
        def __init__(self, n):
            self.n = n

        def propagate(self, *, backward):
            if backward:
                return None

            def receive(s, i):
                if i == lvl_slot:
                    return Message(level_gdna=Level(profile=row, n=3.0, a=100.0))
                if i == comp_slot:
                    return Message(composition=row)
                return Message()

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
    informed = np.asarray(out.informed, bool)
    assert informed[comp_slot], "a received composition must inform"
    assert not informed[lvl_slot], "a level with a row is a bound only and must not inform"


def test_silence_shrinks_the_informed_set_to_own_evidence_and_certainty(sweep_inputs):
    """PERTURBATION: with no messages, the delivered rows vanish and the informed set must shrink to
    the own channel plus structural certainty — and at least one slot must change, or the message
    layer's contribution to the predicate is not being read."""
    policy, *_ = _full_policy(sweep_inputs)
    live = np.asarray(
        SW.solve_chain(*sweep_inputs["args"], **sweep_inputs["kw"], policy=policy).informed, bool
    )
    cap: dict = {}
    silent = SW.solve_chain(
        *sweep_inputs["args"], **sweep_inputs["kw"], policy=SilentPolicy(), _capture=cap
    )
    silent_informed = np.asarray(silent.informed, bool)
    own = has_own_composition_evidence(cap["_tau0_lam"]) | g1_locked(
        cap["free_pos"], cap["free_neg"]
    )
    assert np.array_equal(silent_informed, own)
    assert (live & ~silent_informed).any(), "no slot was informed by a message alone"
    assert not (silent_informed & ~live).any()


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
        var_gdna=np.array([np.inf, 1.0, 2.0, 1.0, 0.0, 1.0, 0.5]),
        informed=np.array([False, True, True, True, True, True, True]),
    )
    return chain, belief, statics, region_arrays, mass, eff


def _training_counts(monkeypatch, belief, parts):
    chain, _, statics, region_arrays, mass, eff = parts
    seen = {}

    def spy(
        count, mass_, eff_, var, *, anchor, strength=1.0, knn_scale=0.5, domain=None, prev=None
    ):
        seen["count"] = np.asarray(count, np.float64)
        seen["anchor"] = np.asarray(anchor, bool)
        seen["domain_mass"] = None if domain is None else np.asarray(domain[0], np.float64)
        return None

    monkeypatch.setattr(CAL, "fit_landscape", spy)
    monkeypatch.setattr(CAL, "_MIN_TRAIN", 1)
    CAL._fit_gdna_hyperprior(chain, belief, statics, region_arrays, mass, eff, strength=1.0)
    return seen


def test_a_flat_likelihood_slot_is_not_in_the_training_population(monkeypatch):
    parts = _synthetic_population()
    chain, belief, *_ = parts
    full = _training_counts(monkeypatch, belief, parts)
    # all three expressed regions plus the anchor train when every slot is informed
    assert full["count"].shape == (4,) and full["anchor"].sum() == 1
    # PERTURBATION: the exon's likelihood was flat -> it leaves the population; nothing else moves
    blind = SimpleNamespace(
        f_g=belief.f_g, var_gdna=belief.var_gdna, informed=belief.informed & ~(np.arange(7) == 2)
    )
    part = _training_counts(monkeypatch, blind, parts)
    assert part["count"].shape == (3,)
    assert np.array_equal(part["count"], np.array([0.0, 80.0, 54.0]))
    assert part["anchor"].sum() == 1
    # THE GRID IS THE CONSUMERS' DOMAIN: every slot with opportunity, boundaries included, reaches the
    # estimator as ``domain`` — the exon that left the training set is still on the axis the prior is
    # read at, and so is every boundary
    assert np.array_equal(part["domain_mass"], np.array([0.0, 5.0, 100.0, 5.0, 80.0, 5.0, 60.0]))


def test_the_substrate_guard_measures_the_domain_not_the_cut(monkeypatch):
    """With the cut leaving only the anchor, the prior still fits as long as the DOMAIN has enough
    substrate — otherwise a gDNA-free toy loses its refit and invents gDNA in its place."""
    parts = _synthetic_population()
    _, belief, *_ = parts
    monkeypatch.setattr(
        CAL, "_MIN_TRAIN", 4
    )  # the domain (4) passes, the cut population (1) would not
    none = SimpleNamespace(f_g=belief.f_g, var_gdna=belief.var_gdna, informed=np.zeros(7, bool))
    chain, _, statics, region_arrays, mass, eff = parts
    calls = []
    monkeypatch.setattr(CAL, "fit_landscape", lambda *a, **k: calls.append(1))
    CAL._fit_gdna_hyperprior(chain, none, statics, region_arrays, mass, eff, strength=1.0)
    assert calls, "the refit was refused although the domain has enough substrate"


def test_the_anchor_trains_whatever_the_predicate_says(monkeypatch):
    parts = _synthetic_population()
    _, belief, *_ = parts
    none = SimpleNamespace(f_g=belief.f_g, var_gdna=belief.var_gdna, informed=np.zeros(7, bool))
    seen = _training_counts(monkeypatch, none, parts)
    assert seen["count"].shape == (1,) and seen["anchor"].all() and seen["count"][0] == 0.0


def test_a_belief_without_the_predicate_trains_the_old_population(monkeypatch):
    """An initial belief (`init_beliefs`, no solve yet) carries no predicate; the selector then reads
    the annotation alone — the population before this rule — so nothing upstream of the first sweep
    changes shape."""
    parts = _synthetic_population()
    _, belief, *_ = parts
    old = SimpleNamespace(f_g=belief.f_g, var_gdna=belief.var_gdna, informed=None)
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
    return DensityLandscape(log_rho=grid_log10 * np.log(10.0), logP=np.log(p / p.sum()), n_train=1)


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
