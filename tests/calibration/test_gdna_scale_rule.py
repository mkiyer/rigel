"""⭐⭐⭐ THE gDNA SCALE RULE — the gates on how a gDNA LEVEL crosses a composition-mismatched step.

The reframe ``r = rho_tot(dst)/rho_tot(src)`` delivers ``phi_c(src)·rho_tot(dst)``: the source's density
SHARE applied to the destination's observed total. It is a composition imputation, and it is exact only
where the two objects share a composition. At an ``intergenic|exon boundary -> EXON`` step the source is
structurally 100 % gDNA and the destination is not, so ``phi_g(src) = 1`` and the delivered gDNA level
collapses to the destination's **own** total — TRAPS: a-message-from-the-destinations-belief, the message carrying zero information and
confirming the destination.

The rule under test: a COMPOSITION crosses by ``r``, licensed by the λ-emission gate's own predicate
(the source must have SUPPLIED both components of the pair); a gDNA LEVEL crosses UNSCALED, because gDNA
is uniform along the genome before capture. Capture needs no branch and no scale factor: it is carried by
the pure-gDNA population's OWN measurements, which the relay's mass pin restores at every such object,
so the level an exon receives is its flanking BOUNDARY's enriched measurement rather than the off-target floor.

⭐ Six gates. The base fixture is the owner's ``TA_single_exon`` geometry in unit-test form:
``intergenic | exon+ | intergenic``, chain ``N E N E N``, so the exon (slot 2) can ONLY be solved through
its two flanking gene-boundary BOUNDARIES (slots 1, 3) — a wrong exon there is a message-passing failure and
nothing else. Gate 5 extends it to an interior exon, which is where TRAPS: a-purity-filter-is-a-length-filter's panel mass actually lives.

⚠ **Perturbation-verified, and here is what they cover** (`scratchpad/perturb_gates.py`): reverting the
rule fires 7 of 8, reverting it in EITHER twin alone fires at least gate 5, licensing every step fires 7,
licensing none fires gate 6, and testing only the gDNA half of the licence fires 7. ⛔ **One perturbation
fires NOTHING and that is recorded rather than papered over**: dropping the gDNA conjunct from the licence
(testing RNA alone) is inert, because ``pg[src] == 0`` with RNA precision live implies the source's own
gDNA density is 0 and the delivered level is then 0 at any scale. See the note beside ``lend`` in
`messages.head`.
"""

from __future__ import annotations

import functools

import numpy as np
import pytest

from rigel.calibration.messages.head import HeadPolicy
from rigel.calibration.sweep import solve_chain

from rigel.calibration.effective_length import (
    UNBOUNDED_REACH,
    contained_eff_length,
    crossing_eff_length,
)
from rigel.calibration.region_geometry import init_beliefs
from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS

from _synthetic import make_chain_parts


#: ⚠ These gates exercise HEADPOLICY's operators, so the policy is named EXPLICITLY. ``solve_chain``
#: defaults to ``SilentPolicy``, which sends nothing — every assertion below would then be vacuous, which
#: is TRAPS: could-the-arm-have-fired exactly ("check the arm COULD have changed something").
region_sweep = functools.partial(solve_chain, policy=HeadPolicy())

#: the chain is ``N E N E N``: intergenic(0) BOUNDARY(1) EXON(2) BOUNDARY(3) intergenic(4)
IG_LEFT, BOUNDARY_LEFT, EXON, BOUNDARY_RIGHT, IG_RIGHT = 0, 1, 2, 3, 4
PURE = (IG_LEFT, BOUNDARY_LEFT, BOUNDARY_RIGHT, IG_RIGHT)


def _delta_pmf(length):
    p = np.zeros(length + 1)
    p[length] = 1.0
    return p


def _single_exon_chain(*, rho_region, rho_boundary, rna_counts, region_bp=1000.0):
    """``intergenic | exon+ | intergenic`` with a gDNA field laid down as ``rho x own opportunity``, plus
    ``rna_counts`` extra UNSTRANDED RNA counts on the exon only.

    ``rho_region`` / ``rho_boundary`` are the gDNA densities of the off-probe REGION class and the exon-adjacent
    BOUNDARY class. Equal ⇒ a flat landscape (capture OFF). Unequal ⇒ a capture step, with the exon's own
    gDNA laid at the on-probe level ``rho_boundary`` — the geometry hybrid capture produces.

    ⚠ The RNA is split evenly across the two GENOME strands, which is what an unstranded library
    deposits, so the strand channel carries exactly zero composition information (`EQUATIONS.md` §5.2)
    and the exon has no own answer — the state TRAPS: a-purity-filter-is-a-length-filter was measured in.
    """
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    region_eff = contained_eff_length(np.full(3, region_bp), gdna_fl)
    unb = np.full(1, UNBOUNDED_REACH)
    boundary_eff = float(crossing_eff_length(gdna_fl, unb, unb)[0])
    # gDNA: the flanks at the off-probe rate, the exon at the on-probe rate (== off-probe when flat).
    g_region = np.array([rho_region, rho_boundary, rho_region]) * region_eff
    g_boundary = rho_boundary * boundary_eff
    region_count = g_region + np.array([0.0, float(rna_counts), 0.0])
    return make_chain_parts(
        [0, BIT_EXON_POS, 0],
        region_size_bp=region_bp,
        region_pos=region_count / 2,
        region_neg=region_count / 2,
        boundary_pos=g_boundary / 2,
        boundary_neg=g_boundary / 2,
        gdna_fl=gdna_fl,
        rna_fl=rna_fl,
    )


#: ⚠ Not the 40 the older fixtures in `test_sweep.py` use. The λ lattice's spacing near ``f_g = ½`` is
#: ``2·window/n_grid`` in log-odds, i.e. ~0.12 in ``f_g`` at 40 points — larger than the residual under
#: test, so a coarse grid would be measuring its own discretization. At 200 the posterior-mean error on
#: this fixture is < 0.006 on every resolvable rung. It is a computational budget, not a modelling choice.
_N_GRID = 200


def _solve(parts, *, kappa=0.5):
    cap = {}
    final = region_sweep(
        parts.chain,
        parts.statics,
        parts.geometry,
        init_beliefs(
            parts.chain, parts.geometry, parts.statics, rna_sense_frac=kappa, n_grid=_N_GRID
        ),
        parts.region_arrays,
        rna_sense_frac=kappa,
        n_grid=_N_GRID,
        _capture=cap,
    )
    count = np.asarray(parts.geometry.unspliced_count, float).sum(axis=1)
    eff_g = np.asarray(parts.geometry.eff_gdna, float)
    return final, cap, final.f_g * count / eff_g


def _sweep(rho_region, rho_boundary, rna_ladder):
    """The delivered gDNA level at the exon, and its own f_g, as the exon's RNA content moves."""
    out = []
    for rna in rna_ladder:
        parts = _single_exon_chain(rho_region=rho_region, rho_boundary=rho_boundary, rna_counts=rna)
        final, cap, rho_g = _solve(parts)
        count = np.asarray(parts.geometry.unspliced_count, float).sum(axis=1)
        eff_g = np.asarray(parts.geometry.eff_gdna, float)
        out.append(
            {
                "rna": rna,
                # the MESSAGE: the fused delivered gDNA density at the exon
                "cg": float(np.asarray(cap["_uni"][-1]["cg"], float)[EXON]),
                # the CONSEQUENCE: the exon's own answer, and the truth the fixture laid down
                "pred_fg": float(final.f_g[EXON]),
                "true_fg": float(rho_boundary * eff_g[EXON] / count[EXON]),
                # the destination's OWN total density — what TRAPS: a-message-from-the-destinations-belief delivers, kept so a failure can name it
                "m_over_eg": float(count[EXON] / eff_g[EXON]),
                "rho_g_out": float(rho_g[EXON]),
                "pure_fg": final.f_g[list(PURE)].copy(),
            }
        )
    return out


# ── the RNA ladder: the exon goes from ~pure gDNA to 1000x RNA, gDNA PINNED ─────────────────────────
#: ~700 gDNA counts on the exon at rho = 1.0, so this spans true f_g 1.00 -> 0.001.
_RNA_LADDER = (0.0, 700.0, 7_000.0, 70_000.0, 700_000.0)


def test_delivered_gdna_level_is_independent_of_the_destinations_own_total():
    """⛔⛔ **GATE 1 — TRAPS: a-message-from-the-destinations-belief, stated as an invariance AND as a value.** The gDNA field is PINNED
    and only the exon's RNA content moves, over 1000x. The level DELIVERED to the exon is a claim about
    gDNA, so it must not move at all — and since the field is uniform it must equal the field.

    Pre-fix this fails by construction: the delivered level is exactly ``rho_tot(dst)``, so it tracks the
    exon's RNA 1:1 — measured spread **982.7x** over this ladder. ⭐ A prediction that moves with data it
    is not about is the same tell as one that does not move with data it IS about."""
    rows = _sweep(1.0, 1.0, _RNA_LADDER)
    cg = np.array([r["cg"] for r in rows])
    assert np.all(cg > 0.0), f"the gDNA message died on the ladder: {cg}"
    spread = cg.max() / cg.min()
    assert spread < 1.001, (
        "the delivered gDNA level tracks the destination's OWN total (TRAPS: a-message-from-the-destinations-belief): "
        f"spread {spread:.1f}x over a 1000x RNA sweep, {cg}"
    )
    assert cg == pytest.approx(1.0, rel=1e-9), f"delivered {cg} against a laid-down field of 1.0"


def test_exon_fg_tracks_its_truth_across_the_density_sweep():
    """⭐⭐ **GATE 2 — the exon's answer must MOVE with the data.** With the gDNA level delivered
    correctly, ``f_g = rho_g·E_g/M`` is arithmetic, so the exon's ``f_g`` must follow its truth down three
    decades. Pre-fix it sat flat at **0.9095** on every rung while truth spanned 1.00 -> 0.001.

    ⚠ **The zero-RNA rung is excluded from the accuracy bound, and that is not a tolerance dodge.** There
    the exon is structurally pure gDNA with no evidence of its own, and ψ's uninformative Jeffreys
    reference deliberately holds such a region off the ``f_g = 1`` vertex until the data earn it — the same
    residual `test_sweep.test_gdna_sweep_factor1_ambig_recovery` measures from the other side, with the
    same instruction not to attack it with damping. The delivered LEVEL is exact there too (gate 1); what
    is short is ψ's readout at the vertex, which is a different mechanism and a different fix."""
    rows = _sweep(1.0, 1.0, _RNA_LADDER)
    pred = np.array([r["pred_fg"] for r in rows])
    true = np.array([r["true_fg"] for r in rows])
    assert np.all(np.diff(pred) < 0.0), f"f_g is not monotone decreasing in RNA density: {pred}"
    assert pred.max() / max(pred.min(), 1e-12) > 100.0, f"f_g barely moves over 1000x RNA: {pred}"
    rna = np.array([r["rna"] for r in rows])
    err = np.abs(pred - true)[rna > 0.0]
    assert err.max() < 0.02, f"|f_g - truth| = {err} (pred {pred}, true {true})"
    # the vertex rung: short by ψ's reference, not by the scale rule, and it must not be far short.
    assert 0.85 < pred[0] < 1.0, f"the zero-RNA exon read {pred[0]:.4f}"


def test_every_structurally_pure_gdna_object_stays_exact():
    """⭐ **GATE 3 — a fix that breaks the anchors has traded one error for another.** Both intergenic
    regions and both gene-boundary BOUNDARIES are structurally pure gDNA; all four must read exactly 1.0 on
    every rung of the ladder."""
    for row in _sweep(1.0, 1.0, _RNA_LADDER):
        assert np.all(row["pure_fg"] == 1.0), f"rna={row['rna']}: pure objects {row['pure_fg']}"


@pytest.mark.parametrize("step", (1.0, 20.0, 200.0))
def test_capture_step_is_carried_and_the_off_probe_floor_is_not(step):
    """⭐⭐⭐ **GATE 4 — THE SAME EXPRESSION UNDER BOTH CAPTURE ARMS, AND THIS IS WHERE 'NO BRANCH' IS
    PROVED.** The off-probe flanks sit ``step`` x below the gene-boundary BOUNDARIES and the exon's own gDNA,
    which is the geometry hybrid capture produces; ``step = 1`` is capture-OFF and is the SAME row of the
    same parametrization, not a separate case.

    The exon must receive the **ON-PROBE** level — its flanking BOUNDARIES' own measurement — on every rung.
    That is a three-way discrimination and each wrong answer has a name:

    * ``step`` (correct) — the BOUNDARY's own enriched measurement;
    * ``1.0`` — the off-probe floor leaking through an unscaled relay that never re-anchors;
    * ``rho_tot(exon)`` — TRAPS: a-message-from-the-destinations-belief again.

    ⭐ Nothing scales the level here: it arrives unscaled and is re-anchored at the BOUNDARY by the relay's
    mass identity, which is what makes the capture landscape a MEASUREMENT rather than a model. Run at
    high RNA so ``rho_tot(exon)`` is far from every candidate."""
    rows = _sweep(1.0, step, (70_000.0,))
    cg, rho_tot = rows[0]["cg"], rows[0]["m_over_eg"]
    assert cg == pytest.approx(step, rel=0.02), (
        f"delivered {cg:.4g}: expected the on-probe level {step:g} "
        f"(1.0 = the off-probe floor; {rho_tot:.4g} = rho_tot(exon), i.e. TRAPS: a-message-from-the-destinations-belief)"
    )
    assert rows[0]["pred_fg"] == pytest.approx(rows[0]["true_fg"], abs=0.05)


def test_level_survives_two_hops_through_an_rna_rich_exon():
    """⭐⭐ **GATE 5 — the RELAY's half of the rule, and the INTERIOR exon.** Gates 1–4 all read the
    COMBINE (`_pin`/`_uni`), and on a five-slot chain every message crosses one step, so the sequential
    relay is never exercised. That is the twin-drift hole the twin note on ``messages.head``'s ``scan`` exists for:
    ``_relay`` and ``_transport`` are two hand-maintained copies of one transform, and a change landed in
    only one of them must fail something.

    The fixture is TRAPS: a-purity-filter-is-a-length-filter's real panel population — an exon with no intergenic neighbour at all:
    ``intergenic | exon+ | exon+ | exon+ | intergenic``. The middle exon carries NO RNA (truth
    ``f_g = 1``) and both its neighbours are RNA-rich, so the correct gDNA level can only reach it by
    relaying THROUGH an RNA-rich object. If any step re-imputes the destination's total on the way, the
    middle exon inherits its neighbours' crowding instead of the gDNA field.

    ⭐⭐ **THE BOUND IS EXACT, AND IT ONLY BECAME EXACT WHEN THE MASS PIN WAS LICENSED TOO.** It used to
    be 1.5x, because a second D4-family mechanism sat on this fixture: the relay's mass pin filled every
    component a message did NOT supply from the destination's OWN density, which at an evidence-free
    object is ψ's uninformative ``fg_loc ~ 1/2``. It therefore reserved about half the budget for RNA no
    message claimed — the running level inflated 1.96x at the RNA-rich exon and deflated again at the next
    evidence-free boundary, a multiplicative random walk landing at 0.80 against a field of 1.0. The pin is
    now gated by the same licence (`test_relay_mass_pin`), so nothing rescales the level on this chain and
    it arrives as the field, exactly.

    ⚠ ``f_g`` at the interior exon is 0.914 and not 1.000, and that residual is NOT this rule or the pin:
    it is ψ's uninformative reference holding an evidence-free region off the vertex by ~0.08, pinned as a
    strict xfail against the TRUTH in `test_relay_mass_pin`. The bound here is short of it deliberately.

    ⭐ **Every variant that re-imputes a total fails this gate**, and did so even at the old loose bound:

    ====================================  ==========  ==========
    arm                                   delivered   f_g (truth 1.000)
    ====================================  ==========  ==========
    the rule, pin licensed                    1.000       0.914
    the rule, pin unlicensed everywhere       0.796       0.853
    ``r_g = r`` (rule reverted)               0.477       0.487
    relay reverted only (twin drift)          0.634       0.658
    combine reverted only (twin drift)        0.599       0.635
    ====================================  ==========  ==========

    ⚠ And under the reverted arms the relay's running level reaches **56.4x** the field at the RNA-rich
    exon before the pin drags it back — the TRAPS: a-message-from-the-destinations-belief re-imputation, visible in flight."""
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    rho, bp_ = 1.0, 1_000.0
    region_eff = contained_eff_length(np.full(5, bp_), gdna_fl)
    unb = np.full(1, UNBOUNDED_REACH)
    boundary_eff = float(crossing_eff_length(gdna_fl, unb, unb)[0])
    # a UNIFORM gDNA field on every object, plus heavy RNA on the two OUTER exons only
    region_count = rho * region_eff + np.array([0.0, 60_000.0, 0.0, 60_000.0, 0.0])
    parts = make_chain_parts(
        [0, BIT_EXON_POS, BIT_EXON_POS, BIT_EXON_POS, 0],
        region_size_bp=bp_,
        region_pos=region_count / 2,
        region_neg=region_count / 2,
        boundary_pos=rho * boundary_eff / 2,
        boundary_neg=rho * boundary_eff / 2,
        gdna_fl=gdna_fl,
        rna_fl=rna_fl,
    )
    final, cap, _ = _solve(parts)
    mid = 4  # chain N E N E N E N E N → the middle exon
    cg = float(np.asarray(cap["_uni"][-1]["cg"], float)[mid])
    assert cg == pytest.approx(rho, rel=1e-9), (
        f"the level reaching the interior exon is {cg:.9g} against a field of {rho:g} — a factor of "
        f"{max(cg / rho, rho / cg):.4f}, so a total was re-imputed somewhere along the relay"
    )
    assert final.f_g[mid] > 0.90, f"the RNA-free interior exon read f_g = {final.f_g[mid]:.4f}"


def test_composition_shared_hops_keep_the_total_density_reframe():
    """⭐⭐ **GATE 7 — the rule is a LICENCE, not a replacement.** Where the source really does supply
    both components of the pair, ``r`` is the shared enrichment ratio and the imputation is sound: the
    gDNA arm must use it, unchanged and bit-identical. Where the source supplies one, it must not.

    The fixture is a STRANDED chain (κ = 0.95, library sample sizes large enough to open the strand
    deadband), so the exon and intron regions earn real composition evidence and can lend a composition,
    while the two gene-boundary BOUNDARIES and the intergenic flanks still cannot.

    ⚠ **Both populations must be non-empty**, and that is what makes the gate two-sided: forcing the
    licence always-on or always-off each empties one of them. Without this gate the single-exon fixture
    cannot tell those two lies apart, because it has no composition-shared step at all."""
    parts = make_chain_parts(  # intergenic | exon+ | intron+ | exon+ | intergenic
        [0, BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS, 0],
        region_size_bp=1000.0,
        region_pos=[100.0, 900.0, 400.0, 900.0, 100.0],  # sense-tilted RNA on the genic regions
        region_neg=[100.0, 50.0, 30.0, 50.0, 100.0],
        boundary_pos=[20.0, 60.0, 60.0, 20.0],
        boundary_neg=[20.0, 10.0, 10.0, 20.0],
    )
    cap = {}
    region_sweep(
        parts.chain,
        parts.statics,
        parts.geometry,
        init_beliefs(
            parts.chain, parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=_N_GRID
        ),
        parts.region_arrays,
        rna_sense_frac=0.95,
        n_rna_obs=10_000.0,
        n_gdna_obs=10_000.0,
        n_grid=_N_GRID,
        _capture=cap,
    )
    lend = np.concatenate([np.asarray(d["lend"], bool) for d in cap["_pin"]])
    r = np.concatenate([np.asarray(d["r"], float) for d in cap["_pin"]])
    r_g = np.concatenate([np.asarray(d["r_g"], float) for d in cap["_pin"]])
    valid = np.concatenate([np.asarray(d["valid"], bool) for d in cap["_pin"]])

    assert (lend & valid).any(), "no step is licensed — the rule has become a blanket replacement"
    assert (~lend & valid).any(), "every step is licensed — the licence is not being applied at all"
    shared = lend & valid
    assert np.array_equal(r_g[shared], r[shared]), (
        f"a composition-shared step changed scale: r_g {r_g[shared]} vs r {r[shared]}"
    )
    un = ~lend & valid
    assert np.all(r_g[un] == 1.0), f"an unlicensed step scaled the gDNA level: {r_g[un]}"
    # …and it must be a real test: `r` has to be far from 1 on some of those steps, or the two branches
    # agree by accident and the fixture proves nothing.
    assert np.abs(np.log(np.maximum(r[un], 1e-300))).max() > 1.0, (
        f"the fixture's r is already ~1 on every unlicensed step: {r[un]}"
    )


# ── ⭐⭐⭐ THE RNA LEVEL LICENCE — the population licence is PER STRAND (owner, 2026-08-18) ─────────
#
# Found by the dissection loop on `g05 ss0.50 capture_off`, in two measured halves. ① At a TSS− boundary
# the licence refused the hop and the gDNA level crossed unscaled — but the RNA arms were still
# multiplied by the FULL total-density ratio, crushing a continuing + population 7x (its density barely
# changes across the OTHER strand's terminus; the total drops because the ENDING population exits): true
# f_g 0.0023 solved to 0.9009. ② The v1 repair — every arm unscaled at a denied hop — fixed those slots
# and FLOODED same-strand tandem genes instead (+114,812 fragments, 52 % in ten slots flanked by
# TES+|TSS+ pairs): where THIS strand's population changes, the source does not measure the
# destination's population at all, and neither `r` (a self-confirming echo,
# TRAPS: a-message-from-the-destinations-belief) nor `1` (the neighbour gene's level, unbounded) is valid.
#
# The rule, one decision per arm: both populations intact ⇒ reframe by r; THIS arm intact but the OTHER
# changed ⇒ r = 1 (a DENSITY transfer); THIS arm changed ⇒ NO CLAIM — value and precisions zeroed.


def _transport_fixture(pop_pair: bool, pop_p: bool, pop_n: bool):
    """One hop through the vectorised `_transport` with explicit pair and per-strand licences."""
    import dataclasses
    from types import SimpleNamespace

    from test_sweep_backbone import _ctx

    from rigel.calibration.messages import NeighbourState

    ctx = _ctx(
        free_neg=np.ones(8, bool)
    )  # both strands admissible: the zeroing gate is not under test
    n = ctx.n_slots
    z = np.zeros(n)
    own = SimpleNamespace(
        rho_g=z,
        rho_pos=z,
        rho_neg=z,
        prec_g=z,
        prec_pos=z,
        prec_neg=z,
        tau_lam=z,
        struct_lock=np.zeros(n, bool),
        f_g=np.full(n, 0.5),
    )
    two = np.zeros((n, 2))
    geom = SimpleNamespace(
        n_slots=n,
        unspliced_count=np.full((n, 2), 50.0),
        eff_gdna=np.full(n, 200.0),
        eff_rna=np.full(n, 200.0),
        spliced_count=two,
        sj_count=two,
        eff_sj=np.full((n, 2), 200.0),
        sj_count_lo=two,
        sj_count_hi=two,
        eff_sj_lo=np.full((n, 2), 200.0),
        eff_sj_hi=np.full((n, 2), 200.0),
    )
    ctx = dataclasses.replace(ctx, own=own, geometry=geom)
    relay = HeadPolicy().prepare(ctx)
    src = np.arange(n) - 1
    valid = src >= 0
    src = np.clip(src, 0, n - 1)
    rg, rp, rn = np.full(n, 0.05), np.full(n, 1.60), np.full(n, 0.80)
    pg, pp, pn = np.full(n, 2.0), np.full(n, 3.0), np.full(n, 3.0)
    zero = np.zeros(n)
    nb = NeighbourState(
        state=(rg, rp, rn, pg, pp, pn, pg.copy(), pp.copy(), pn.copy(), zero), valid=valid, src=src
    )
    out = relay._transport(
        nb,
        np.full(n, 2.0),
        np.full(n, 10.0),  # a 5x frame step: r = 0.2 wherever framed
        np.full(n, pop_pair),
        np.full(n, pop_p),
        np.full(n, pop_n),
    )
    return valid, (rg, rp, rn), out


def test_a_terminus_on_one_strand_silences_that_arm_and_frees_the_other():
    """⛔ A − terminus: the − arm delivers NO CLAIM (value AND precisions zero); the + arm crosses as a
    DENSITY (r = 1, because the total ratio is corrupted by the − population's exit); the gDNA level
    crosses unscaled (the pair licence is refused)."""
    valid, (rg, rp, rn), out = _transport_fixture(False, True, False)
    tg, tp, tn, tpp, tpn = out[0], out[1], out[2], out[4], out[5]
    tmn = out[8]
    assert np.allclose(tg[valid], rg[valid]), "gDNA level scaled on an unlicensed pair"
    assert np.allclose(tp[valid], rp[valid]), "the intact + arm must cross as a DENSITY (r = 1)"
    assert np.all(tn[valid] == 0.0), "the changed − arm delivered a value"
    assert np.all(tpn[valid] == 0.0) and np.all(tmn[valid] == 0.0), (
        "the changed − arm delivered PRECISION — a zero at live precision is the claim 'there is no "
        "− RNA', which the source cannot make about a population it does not measure"
    )
    assert np.all(tpp[valid] > 0.0), "the intact + arm lost its precision"


def test_a_terminus_on_both_strands_silences_both_arms():
    """⛔ TSS+ and TES− at one boundary: neither RNA population crosses; both arms deliver nothing."""
    valid, (rg, _rp, _rn), out = _transport_fixture(False, False, False)
    tg, tp, tn, tpp, tpn = out[0], out[1], out[2], out[4], out[5]
    assert np.allclose(tg[valid], rg[valid])
    assert np.all(tp[valid] == 0.0) and np.all(tn[valid] == 0.0)
    assert np.all(tpp[valid] == 0.0) and np.all(tpn[valid] == 0.0)


def test_a_clean_hop_still_reframes_every_arm():
    """⭐ No population change anywhere ⇒ the reframe is licensed and EVERY arm scales by r — the rule
    must not turn the composition sharing off wholesale."""
    valid, (rg, rp, rn), out = _transport_fixture(True, True, True)
    tg, tp, tn = out[0], out[1], out[2]
    r = 0.2
    assert np.allclose(tg[valid], rg[valid] * r)
    assert np.allclose(tp[valid], rp[valid] * r)
    assert np.allclose(tn[valid], rn[valid] * r)


# ── ⭐⭐⭐ THE GRAFT AT A NO-CLAIM HOP — found by the zero-control dissection (2026-08-18, session 2) ────
#
# `g00 ss0.50 nrna_mid capture_on` (62 k → 228 k under the per-strand licence) walked to its mechanism:
# a zero RNA claim WITH LIVE PRECISION, relayed one hop, where the mass pin's licence (`_lend`: both
# components supplied) then rescales the whole budget onto gDNA — ``k = M/(tg·E_g)`` — and delivers
# ``M/E_g``, "all your mass is gDNA", into a 26 k-fragment exon (slot 3568: true 0.000 → 0.9037). One of
# its two feeders is the graft: at a hop the per-strand rule refuses, the arm's VALUE (which had the sj
# density ``gp`` folded in) is zeroed, and the graft block then adds the sj COUNT'S precision back onto
# the zeroed arm — a confident "there is no RNA on this strand", the exact claim
# `test_a_terminus_on_one_strand_silences_that_arm_and_frees_the_other` forbids, reached through the
# graft path its fixture (no sj) never exercised. Measured at 56584 → 56583: ``n = 0 @ pn = 280``.
#
# The rule stands as landed — THIS strand changed ⇒ NO CLAIM — and the graft is a claim about this strand.
# ψ reads ``rna_mode = log(cp·E_r/M)`` as a LEVEL, so delivering the sj density alone would be a partial
# level (RNA biased low ⇒ false gDNA), and the exon's own certified spliced count is already in its own
# evidence. So a refused arm delivers NOTHING, graft included: value AND precision zero, both twins.


def _graft_transport_fixture(pop_p: bool):
    """One hop through `_transport` where every even slot is an EXON and every odd slot a BOUNDARY
    carrying + sj flux, so the hop odd → even is a GRAFT; the + arm's licence is the argument."""
    import dataclasses
    from types import SimpleNamespace

    from test_sweep_backbone import _ctx

    from rigel.calibration.messages import NeighbourState

    ctx = _ctx(free_neg=np.ones(8, bool))
    n = ctx.n_slots
    is_bnd = np.asarray(ctx.is_boundary, bool)
    sj = np.zeros((n, 2))
    sj[is_bnd, 0] = 40.0  # + sj flux at every boundary
    ctx = dataclasses.replace(ctx, is_exon_region=~is_bnd, sj_count=sj)
    z = np.zeros(n)
    own = SimpleNamespace(
        rho_g=z,
        rho_pos=z,
        rho_neg=z,
        prec_g=z,
        prec_pos=z,
        prec_neg=z,
        tau_lam=z,
        struct_lock=np.zeros(n, bool),
        f_g=np.full(n, 0.5),
    )
    two = np.zeros((n, 2))
    geom = SimpleNamespace(
        n_slots=n,
        unspliced_count=np.full((n, 2), 50.0),
        eff_gdna=np.full(n, 200.0),
        eff_rna=np.full(n, 200.0),
        spliced_count=two,
        sj_count=sj,
        eff_sj=np.full((n, 2), 200.0),
        sj_count_lo=two,
        sj_count_hi=two,
        eff_sj_lo=np.full((n, 2), 200.0),
        eff_sj_hi=np.full((n, 2), 200.0),
    )
    ctx = dataclasses.replace(ctx, own=own, geometry=geom)
    relay = HeadPolicy().prepare(ctx)
    src = np.arange(n) - 1
    valid = src >= 0
    src = np.clip(src, 0, n - 1)
    rg, rp, rn = np.full(n, 0.05), np.full(n, 1.60), np.full(n, 0.80)
    pg, pp, pn = np.full(n, 2.0), np.full(n, 3.0), np.full(n, 3.0)
    zero = np.zeros(n)
    nb = NeighbourState(
        state=(rg, rp, rn, pg, pp, pn, pg.copy(), pp.copy(), pn.copy(), zero), valid=valid, src=src
    )
    out = relay._transport(
        nb,
        np.full(n, 2.0),
        np.full(n, 10.0),
        np.full(n, pop_p),
        np.full(n, pop_p),
        np.ones(n, bool),
    )
    graft = valid & ~is_bnd & is_bnd[src]
    assert graft.any(), "the fixture built no graft hop"
    return graft, out


def test_the_graft_delivers_no_precision_on_a_refused_arm_VECTORISED():
    """⛔ At a no-claim hop the graft's sj precision must go with the arm's value — a zero at the sj
    count's precision is 'there is no + RNA here', which the source cannot claim."""
    graft, out = _graft_transport_fixture(pop_p=False)
    tp, tpp, tmp = out[1], out[4], out[7]
    assert np.all(tp[graft] == 0.0), "a refused + arm delivered a value"
    assert np.all(tpp[graft] == 0.0), (
        f"the graft added MODE-FUSION precision back onto a refused arm: {tpp[graft]}"
    )
    assert np.all(tmp[graft] == 0.0), (
        f"the graft added MEASUREMENT precision back onto a refused arm: {tmp[graft]}"
    )
    # …and the control: with the arm licensed the SAME fixture's graft is live, so the gate above is not
    # passing vacuously (TRAPS: could-the-arm-have-fired).
    graft, out = _graft_transport_fixture(pop_p=True)
    assert np.all(out[4][graft] > 0.0) and np.all(out[7][graft] > 0.0), (
        "the control's graft is dead"
    )


def _graft_scan_fixture(refuse: bool):
    """The SCALAR twin: the same chain through `scan(backward=False)`, the + licence set from the
    boundary flags (TSS+ at every boundary ⇒ the right flank gains + RNA ⇒ the + arm into every exon is
    refused; 0 ⇒ licensed)."""
    import dataclasses
    from types import SimpleNamespace

    from test_sweep_backbone import _ctx

    from rigel.calibration.splice_graph import FLAG_TSS_POS

    ctx = _ctx(free_neg=np.ones(8, bool))
    n = ctx.n_slots
    is_bnd = np.asarray(ctx.is_boundary, bool)
    sj = np.zeros((n, 2))
    sj[is_bnd, 0] = 40.0
    flags = np.where(is_bnd, int(FLAG_TSS_POS) if refuse else 0, 0).astype(np.uint16)
    z = np.zeros(n)
    own = SimpleNamespace(
        rho_g=z,
        rho_pos=z,
        rho_neg=z,
        prec_g=z,
        prec_pos=z,
        prec_neg=z,
        tau_lam=z,
        struct_lock=np.zeros(n, bool),
        f_g=np.full(n, 0.5),
    )
    two = np.zeros((n, 2))
    geom = SimpleNamespace(
        n_slots=n,
        unspliced_count=np.full((n, 2), 50.0),
        eff_gdna=np.full(n, 200.0),
        eff_rna=np.full(n, 200.0),
        spliced_count=two,
        sj_count=sj,
        eff_sj=np.full((n, 2), 200.0),
        sj_count_lo=two,
        sj_count_hi=two,
        eff_sj_lo=np.full((n, 2), 200.0),
        eff_sj_hi=np.full((n, 2), 200.0),
    )
    ctx = dataclasses.replace(
        ctx, is_exon_region=~is_bnd, sj_count=sj, boundary_flags=flags, own=own, geometry=geom
    )
    relay = HeadPolicy().prepare(ctx)
    step, publish = relay.scan(backward=False)
    for i in range(1, n):
        step(i - 1, i)
    rg, rp, rn, pg, pp, pn, mg, mp, mn, tau = publish()
    exons = ~is_bnd & (np.arange(n) >= 2)  # exons whose left neighbour is a sj-carrying boundary
    return exons, rp, pp, mp


def test_the_graft_delivers_no_precision_on_a_refused_arm_SCALAR():
    """The scan twin of the gate above; the two must agree (`TRAPS: two twins`)."""
    exons, rp, pp, mp = _graft_scan_fixture(refuse=True)
    assert np.all(pp[exons] == 0.0), f"scan: graft precision survived a refused + arm: {pp[exons]}"
    assert np.all(mp[exons] == 0.0), f"scan: graft MEASUREMENT precision survived: {mp[exons]}"
    exons, rp, pp, mp = _graft_scan_fixture(refuse=False)
    assert np.all(pp[exons] > 0.0) and np.all(mp[exons] > 0.0), "the scan control's graft is dead"
