"""⭐⭐⭐ THE RELAY'S MASS PIN IS GATED BY THE COMPOSITION-IMPUTATION LICENCE.

The relay used to restore ``Sum_c rho_c·E_c = M`` at **every** slot, filling components the message did
NOT supply from the **destination's own** density (``messages.head``'s scan kernel, the ``_k`` block; the scalar
twin of ``_pin_v``). For a gDNA-only message that rescale had an exact closed form in two published
quantities::

    k = 1 / (phi_msg + R_own)        phi_msg = (level in) x E_g / M        R_own = (op + on) x E_r / M

``R_own`` is the RNA share of the destination's OWN message-free self-solve, and at an object with no
composition evidence it is **exactly 0.5000** — psi's uninformative reference gives the free strand
exactly half. So the pin reserved half the mass budget for RNA the message never claimed, the delivered
gDNA level was doubled, and because ``k`` depended on ``phi_msg`` the level then random-walked along a
chain: in fraction space the map ``phi -> phi / (phi + R_own)`` has fixed points ``0`` (repelling) and
``1 - R_own`` = 1/2 (attracting), so **the pin pushed every message's composition toward the
destination's own uninformative half.** In the limit ``R_own -> 0`` it was TRAPS: a-message-from-the-destinations-belief at full
strength: ``k -> 1/phi_msg`` and the delivered level became ``M/E_g``, the destination's own total
density with the message's value cancelled out algebraically.

⭐⭐⭐ **THE FIX IS A DERIVATION RATHER THAN A PATCH, AND IT IS ONE RULE: NO BELIEF MAY REACH THE
BUDGET.** ``Sum_c rho_c·E_c = M`` is an identity **under the imputation premise** — a matched reframe
delivers ``rho_c^msg = phi_c·rho_tot(dst) = rho_c^dst,true``, so the components account for exactly the
fragments the destination observed. The pin restores it with ``k = M/S``, and ``S`` fills every component
the message did not supply from the destination's own density. TRAPS: a-message-from-the-destinations-belief permits a message to use the
destination's CONSTANTS and its OBSERVATIONS and never its BELIEFS, so the pin is licensed in exactly the
two states where nothing believed reaches ``S``:

    (i)  the message SUPPLIED the composition — the λ-emission gate's own predicate (`EQUATIONS.md` §3.5),
         so the premise is granted and the identity is implied. Nothing is filled in;
    (ii) the destination is a structurally pure-gDNA object — there is no unsupplied component to fill in.
         ``f_g = 1`` there is STRUCTURE, so ``S = rho_g·E_g`` and the pin hands the object its own MEASURED
         density ``M/E_g``.

Anywhere else, the identity is not implied and there is nothing for the pin to restore.

⭐ **AND IT IS NOT A DELETION — case (ii) is load-bearing, not a caveat.** A pure-gDNA BOUNDARY measures the
gDNA density at its own capture stratum and the exon behind it has no other way to hear it, so dropping
that branch delivers the off-probe floor to every exon under capture (measured: it fires
`test_gdna_scale_rule`'s capture gate at 20x and 200x). It is also why a per-capture-class landscape ratio
built to do the same job measured inert (`region_geometry`'s deleted-landscape note). ⛔ And deleting the pin
outright is a third, worse thing: unanchored, the relay's residual compounded to ``Sum_c rho_c E_c / M``
p99 = 31–288x and max 519x, and rescaling all three components blindly regressed capture-OFF 3.6x.
`test_the_pin_fires_ONLY_where_no_belief_can_reach_its_budget` is what keeps this a licence.

⭐ Measured on the simulated toy chromosome (`toy_harness.py --spec nested_exons` — three NESTED
single-exon transcripts, so the gene is five exon REGIONS and four ``exon|exon`` BOUNDARIES with no intron
anywhere; donor `g75 ss0.50 capture_off`; gDNA uniform at 0.0769 counts/bp). Before, the delivered level
rose monotonically and symmetrically with step distance from the gene ends — **0.071 → 0.102 → 0.097 →
0.162 → 0.192 → 0.215**, reaching **2.8x** the truth, and the gene's mass-weighted ``|Δf_g|`` was
**0.2618**. Licensing the pin takes it to **0.2264**, and the whole 36-condition gDNA ladder moves
0.0415 → 0.0414 with **no condition regressing**.

⚠⚠ **AND THE REST OF THE WAY IS BLOCKED BY TRAPS: pure-and-length-censored/TRAPS: a-threshold-on-a-fitted-residue, NOT BY THIS.** On the same toy the pin still fires
through the gene interior, because ``lend`` is a PRECISION test with no floor and κ is FITTED at 0.500689
rather than exactly ½ on a nominally unstranded library — so ``I(f_g) ∝ (2κ−1)²`` lands at 1e-6…1e-4
instead of 0 and every evidence-free exon "supplies" a composition whose own statement is 10²–10³ nats
wide against a ±10-nat grid. Switching the pin off entirely reaches 0.0760 on this toy; the difference is
that leak. ⛔ A precision floor here would be a tuned constant and is not this change (`ROADMAP.md` §1 **the-capture-level-residual**).

⚠ **Perturbation-verified**, and here is the coverage (`scratchpad/perturb.sh`). ⭐ Note that dropping
the composition branch is caught by gates in OTHER files — so both branches have an independently
measured consequence and neither is decoration:

    ``if g1_l[i]``                (drop the composition branch)  2: test_sweep's
                                                                   mature_measurement_disagreement_silenced,
                                                                   test_toy_harness's intron_composition_dependence
    ``if _lend``                  (drop the structural branch)   3: test_gdna_scale_rule's capture gate at
                                                                   20x and 200x, + the own-measurement gate here
    ``if True``                   (the pre-fix behaviour)        7: every effect gate in this file
    ``if False``                  (delete the pin)               5: the union of the two above
    ``_struct`` for ``g1_locked`` (region-only, so not the BOUNDARY)   3: as ``if _lend``
    ``and`` for ``or``                                           5: the union
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
from rigel.calibration.region_geometry import g1_locked, init_beliefs
from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS

from _synthetic import make_chain_parts


#: ⚠ These gates exercise HEADPOLICY's operators, so the policy is named EXPLICITLY. ``solve_chain``
#: defaults to ``SilentPolicy``, which sends nothing — every assertion below would then be vacuous, which
#: is TRAPS: could-the-arm-have-fired exactly ("check the arm COULD have changed something").
region_sweep = functools.partial(solve_chain, policy=HeadPolicy())

_N_GRID = 200


def _delta_pmf(length):
    p = np.zeros(length + 1)
    p[length] = 1.0
    return p


def _uniform_field_chain(*, rho=1.0, rna, bp=1000.0, rho_first=None):
    """``intergenic | exon+ x k | intergenic`` with gDNA laid down at ``rho`` x each object's own
    opportunity, plus ``rna[i]`` extra unstranded RNA counts on region ``i``.

    ⭐ The gDNA field is UNIFORM by construction, so the truth is ``rho`` at every object and any
    departure of the relayed level from ``rho`` is the relay's own doing. ``rho_first`` overrides the
    field on the FIRST region only — the handle for the TRAPS: a-message-from-the-destinations-belief invariance test.

    ⚠ **No step on this chain is licensed**, and that is the point of it: the library is unstranded and
    there is no sj anywhere, so no slot ever earns RNA precision of its own and no source can lend
    a composition. It is therefore the pure-unlicensed fixture. The licensed population lives on
    `_stranded_chain`.
    """
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    n = len(rna)
    region_eff = contained_eff_length(np.full(n, bp), gdna_fl)
    unb = np.full(1, UNBOUNDED_REACH)
    boundary_eff = float(crossing_eff_length(gdna_fl, unb, unb)[0])
    field = np.full(n, float(rho))
    if rho_first is not None:
        field[0] = float(rho_first)
    region_count = field * region_eff + np.asarray(rna, float)
    return make_chain_parts(
        [0] + [BIT_EXON_POS] * (n - 2) + [0],
        region_size_bp=bp,
        region_pos=region_count / 2,
        region_neg=region_count / 2,
        boundary_pos=rho * boundary_eff / 2,
        boundary_neg=rho * boundary_eff / 2,
        gdna_fl=gdna_fl,
        rna_fl=rna_fl,
    )


def _stranded_chain():
    """``intergenic | exon+ | intron+ | exon+ | intergenic`` on a STRANDED library.

    ⭐ **The two-sided fixture**: the genic regions earn real composition evidence from the strand tilt and
    can therefore lend one, while the intergenic flanks and the gene-boundary BOUNDARIES cannot. So the same
    chain carries licensed and unlicensed steps, which is what makes "the pin fires iff licensed"
    falsifiable in both directions. Same geometry as `test_gdna_scale_rule`'s licence gate.
    """
    return make_chain_parts(
        [0, BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS, 0],
        region_size_bp=1000.0,
        region_pos=[100.0, 900.0, 400.0, 900.0, 100.0],  # sense-tilted RNA on the genic regions
        region_neg=[100.0, 50.0, 30.0, 50.0, 100.0],
        boundary_pos=[20.0, 60.0, 60.0, 20.0],
        boundary_neg=[20.0, 10.0, 10.0, 20.0],
    )


def _relay_walk(parts, *, kappa=0.5, n_obs=0.0):
    """Solve, then read the FORWARD relay's running gDNA level and everything the closed form needs.

    Returns a list of per-slot dicts describing the step INTO that slot from its left neighbour.

    * ``lend`` — the composition-imputation licence for this step, recomputed here from the relay's own
      running precisions at the source (``fwd_pg``/``fwd_pp``/``fwd_pn``), which is the state ``_lend``
      reads: the forward pass visits the chain in genomic order, so the source's entry is final before
      the step into ``i`` is taken.
    * ``r`` — the reframe ratio for the step, ``rho_tot(dst)/rho_tot(src)``. On an unlicensed step the
      gDNA arm's own scale ``r_g`` is 1 and the level crosses unscaled (`EQUATIONS.md` §3.5), so
      ``level_out == level_in`` exactly wherever the destination has no own gDNA precision to fuse in.
    * ``g1`` — is the DESTINATION a structurally pure-gDNA object? Case (ii) of the pin's licence.
    * ``pinned`` — the pin's licence for this step, ``lend or g1``: the two states in which no belief
      reaches the budget.
    * ``S_over_M`` — the destination's mass identity as the relay left it. 1 means the pin fired.
    """
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
        n_rna_obs=n_obs,
        n_gdna_obs=n_obs,
        n_grid=_N_GRID,
        _capture=cap,
    )
    st = cap["_uni_static"]
    M, EG, ER, fwd = st["M"], st["E_g"], st["E_r"], st["fwd_g"]
    op, on, pg_own = st["op"], st["on"], st["pg_own"]
    # ⭐ the reframe frame is a PAIR, one total per FLANK, and a hop pairs them by ROLE: this walk follows
    # the FORWARD relay, whose source is the genomic-LOW neighbour, so the destination presents its
    # LOW-flank total and the source its HIGH-flank one (`region_total_density`).
    rho_lo, rho_hi = st["rho_lo"], st["rho_hi"]
    left = np.asarray(parts.chain.left, np.int64)
    g1 = g1_locked(parts.statics.free_pos, parts.statics.free_neg)
    out = []
    for s in range(M.shape[0]):
        p = int(left[s])
        has_src = p >= 0
        level_in = float(fwd[p]) if has_src else float("nan")
        supplied = st["fwd_g"][s] * EG[s] + (st["fwd_p"][s] + st["fwd_n"][s]) * ER[s]
        lend = (
            bool(st["fwd_pg"][p] > 0.0 and (st["fwd_pp"][p] + st["fwd_pn"][p]) > 0.0)
            if has_src
            else None
        )
        out.append(
            {
                "slot": s,
                "is_region": bool(np.asarray(parts.chain.kind)[s] == 0),
                "M_over_Eg": float(M[s] / EG[s]),
                "R_own": float((op[s] + on[s]) * ER[s] / M[s]),
                "own_prec": float(pg_own[s]),
                "lend": lend,
                "g1": bool(g1[s]),
                "pinned": (bool(lend or g1[s]) if has_src else None),
                "r": float(rho_lo[s] / rho_hi[p])
                if (has_src and rho_hi[p] > 0.0)
                else float("nan"),
                "level_in": level_in,
                "phi_msg": level_in * EG[s] / M[s],
                "level_out": float(fwd[s]),
                "S_over_M": float(supplied / M[s]) if M[s] > 0.0 else float("nan"),
                "f_g": float(final.f_g[s]),
            }
        )
    return out, final, cap


#: three exons, RNA on the OUTER two, the middle one RNA-FREE. Truth is f_g = 1 at every object.
_RNA = (0.0, 60_000.0, 0.0, 60_000.0, 0.0)
_MID = 4  # chain N E N E N E N E N → the middle exon


def _steps(walk, *, pinned):
    """The steps whose delivered claim is READABLE and whose pin licence is ``pinned``.

    ⚠ A destination with own gDNA precision fuses its own belief into the level, so the running value
    there is not the delivered claim and no exact statement can be made about it. Those slots are
    excluded — from both arms, by the same rule.
    """
    return [
        r
        for r in walk
        if r["pinned"] is pinned
        and r["own_prec"] == 0.0
        and np.isfinite(r["level_in"])
        and r["level_in"] > 0.0
    ]


# ── 1. THE LICENCE ─────────────────────────────────────────────────────────────────────────────────


def test_an_unpinned_step_does_not_rescale_the_level_at_all():
    """⭐⭐ **THE FIX, as an exact statement.** Where no licence applies there is no imputation premise,
    so ``Sum_c rho_c·E_c = M`` is not implied and the pin has nothing to restore: the delivered gDNA
    level must arrive **exactly** as the source measured it. Not approximately — an unpinned gDNA arm
    applies no scale of any kind, so this is float equality.

    ⛔ Pre-fix this failed at every step: the pin rescaled by ``1/(phi_msg + R_own)``, which on this
    fixture is ~2 at an evidence-free destination and ``1/phi_msg`` at a pure-gDNA one."""
    walk, _, _ = _relay_walk(_uniform_field_chain(rna=_RNA))
    steps = _steps(walk, pinned=False)
    assert len(steps) >= 5, f"only {len(steps)} readable unpinned steps in the fixture"
    for row in steps:
        assert row["level_out"] == row["level_in"], (
            f"slot {row['slot']}: an UNPINNED step scaled the level "
            f"{row['level_in']:.9g} → {row['level_out']:.9g} "
            f"(k = {row['level_out'] / row['level_in']:.9f}; the pin's old closed form gives "
            f"{1.0 / (row['phi_msg'] + row['R_own']):.9f})"
        )


def test_the_pin_fires_ONLY_where_no_belief_can_reach_its_budget():
    """⭐⭐⭐ **THIS IS A LICENCE AND NOT A DELETION, and this gate is the whole difference.** The pin's
    budget fills every component the message did not supply from the destination's own density, so
    TRAPS: a-message-from-the-destinations-belief permits it in exactly two states — and this gate asserts both of them fire and that
    nothing else does.

    * **(ii) the destination is structurally pure gDNA.** There is no unsupplied component to fill in:
      ``f_g = 1`` is STRUCTURE, so the budget is ``rho_g·E_g`` and the pin hands the object its own
      MEASURED density ``M/E_g``. Both ``intergenic|exon`` BOUNDARIES of this fixture are that object.
    * **not licensed anywhere else on this chain** — unstranded, no sj, so no slot earns RNA
      precision and nothing can lend a composition. Every other step must leave the level alone.

    Two-sided on ONE fixture, because either lie empties one of the two populations: pinning everywhere
    (the pre-fix behaviour) puts a rescale on the second population, and pinning nowhere leaves the
    first with a violated identity.

    ⚠ Case (i) — a licensed source at a NON-pure-gDNA destination — cannot be READ on any fixture: a
    source only earns RNA precision where it has composition evidence, and on a chain where it does, the
    destinations have their own precision too and fuse it into the level, so the delivered claim is not
    observable. `test_the_lend_branch_of_the_licence_is_live` gates that branch structurally instead."""
    walk, _, _ = _relay_walk(_uniform_field_chain(rna=_RNA))
    pin, off = _steps(walk, pinned=True), _steps(walk, pinned=False)
    assert pin, "nothing is pinned — the licence has become a blanket removal"
    assert off, "everything is pinned — the licence is not being applied at all"
    assert all(r["g1"] and not r["lend"] for r in pin), (
        "this fixture is meant to pin ONLY on structural certainty; a source lent a composition on it"
    )
    for row in pin:
        assert row["level_out"] == pytest.approx(row["M_over_Eg"], rel=1e-12), (
            f"slot {row['slot']} is structurally pure gDNA, so the pin must hand it its own measured "
            f"density {row['M_over_Eg']:.9g}; it delivered {row['level_out']:.9g}"
        )
        assert row["S_over_M"] == pytest.approx(1.0, rel=1e-12), (
            f"slot {row['slot']}: the mass identity reads {row['S_over_M']:.9f} after the pin"
        )
    for row in off:
        assert row["level_out"] == row["level_in"], (
            f"slot {row['slot']} has no licence, yet the level was scaled "
            f"{row['level_in']:.9g} → {row['level_out']:.9g}"
        )
        assert abs(row["S_over_M"] - 1.0) > 0.05, (
            f"slot {row['slot']}: the identity happens to hold ({row['S_over_M']:.6f}) with the pin "
            "off, so the unpinned arm of this gate proves nothing"
        )


def test_the_lend_branch_of_the_licence_is_live():
    """⭐ **CASE (i) IS REACHABLE AND IS NOT SHADOWED BY CASE (ii).** The pin's licence is
    ``lend or g1``, and a gate set that only ever exercises the ``g1`` branch would pass with ``lend``
    deleted. This asserts the composition branch fires on steps whose destination is NOT structurally
    pure gDNA — i.e. that it decides something the structural branch does not.

    The fixture is STRANDED, so the genic regions earn composition evidence of their own and can lend one,
    while the intergenic flanks and the gene-boundary BOUNDARIES cannot. ⚠ The effect on the delivered level
    is not readable here (see the note in the gate above); what is asserted is the predicate."""
    walk, _, _ = _relay_walk(_stranded_chain(), kappa=0.95, n_obs=10_000.0)
    steps = [r for r in walk if r["pinned"] is not None]
    lend_only = [r for r in steps if r["lend"] and not r["g1"]]
    neither = [r for r in steps if not r["lend"] and not r["g1"]]
    assert lend_only, "no step is licensed by composition alone — `lend` decides nothing here"
    assert neither, "every step carries some licence — the fixture proves no exclusion"
    assert [r for r in steps if r["g1"]], "no structurally pure-gDNA destination in the fixture"


def test_an_evidence_free_object_still_reserves_half_a_budget_that_no_longer_rescales_anything():
    """⭐⭐ **WHY THE OLD FACTOR WAS ~2, kept because the quantity is unchanged and still true.** At an
    object with no composition evidence psi's uninformative reference puts **exactly one half** on the
    free RNA strand, so ``R_own`` is 0.5000 to machine precision and the pin's old factor
    ``1/(phi_msg + R_own)`` went to 2 as the message's gDNA fraction went to 0. Half of every
    evidence-free object's mass budget was reserved for RNA that no message claimed.

    ⭐ ``R_own`` is psi's *reference*, not a defect, and this gate asserts it is **still exactly a
    half** — the fix changed the pin's licence, not the reference. What it also asserts is that the half
    no longer reaches the delivered level: the budget is still computable, and nothing consumes it."""
    walk, _, _ = _relay_walk(_uniform_field_chain(rna=_RNA))
    free = [r for r in walk if r["own_prec"] == 0.0 and r["R_own"] > 0.0]
    assert free, "no evidence-free object in the fixture"
    for row in free:
        assert row["R_own"] == pytest.approx(0.5, abs=1e-9), (
            f"slot {row['slot']}: R_own = {row['R_own']:.9f}, not the reference's exact half"
        )
    # the old factor is still ~2 on the smallest-phi step — and is now applied to nothing.
    small = min(free, key=lambda r: r["phi_msg"])
    assert 1.0 / (small["phi_msg"] + small["R_own"]) == pytest.approx(2.0, rel=0.05), (
        f"slot {small['slot']}: phi_msg = {small['phi_msg']:.4g} so the old k would be ~2"
    )
    assert small["level_out"] == small["level_in"], (
        f"slot {small['slot']}: the half-budget still rescaled the level "
        f"({small['level_in']:.6g} → {small['level_out']:.6g})"
    )


# ── 2. THE CONSEQUENCE ─────────────────────────────────────────────────────────────────────────────


def test_a_uniform_gdna_field_relays_unchanged():
    """⭐⭐ **THE DEFECT, NOW THE GATE.** The gDNA field is uniform, so the level the relay carries must
    be that field at every slot — there is nothing for it to be but ``rho``.

    Measured before the licence: **1.000 → 1.000 → 1.955 → 0.796 → 0.614 → 0.551 → 1.089 → 1.000**, i.e.
    up to **1.96x** and down to **0.55x** the field within one locus. ⚠ The bound is 1e-9 rather than a
    few per cent because nothing on this chain is licensed and no other operator touches the gDNA arm:
    the relayed level is the field exactly, or a scale has crept back in."""
    rho = 1.0
    walk, _, _ = _relay_walk(_uniform_field_chain(rna=_RNA, rho=rho))
    delivered = np.array([r["level_out"] for r in walk if r["own_prec"] == 0.0])
    worst = float(np.max(np.maximum(delivered / rho, rho / delivered)))
    assert worst < 1.0 + 1e-9, (
        f"the relayed level departs from the uniform field by up to {worst:.6f}x: "
        f"{np.round(delivered, 6)}"
    )


def test_the_running_product_of_the_rescales_is_exactly_one_at_every_step():
    """⭐ **THE CAUSAL CHAIN, closed arithmetically so no source ablation is needed** — plus the reason
    the defect survived so long.

    The level the relay carries is the field times the *running* product of the per-step rescales. With
    the pin licensed, every step on this chain has a rescale of exactly 1, so that product is 1 at every
    slot and the level is the field. Nothing else in the relay touches the gDNA arm.

    ⚠⚠ **AND THIS IS WHY NO AGGREGATE COULD SEE THE OLD DEFECT.** The old product TELESCOPED back to 1
    at the far end of the locus (``_relay``'s own comment records the pairwise version, "exon→bnd 0.431 x
    bnd→exon 2.298 ≈ 1"), and the last step into a structurally-pure-gDNA object rewrote the level to
    that object's own total anyway, erasing the history. The excursion was entirely in the MIDDLE of the
    chain, so it was visible only per object and only away from a pure-gDNA object. ⭐ The gate on the
    fix is therefore that there is **no excursion to telescope** — a peak of 1, not a return to 1."""
    walk, _, _ = _relay_walk(_uniform_field_chain(rna=_RNA, rho=1.0))
    prod, field, excursion = 1.0, 1.0, []
    for row in _steps(walk, pinned=False) + _steps(walk, pinned=True):
        prod *= row["level_out"] / row["level_in"]
        excursion.append(prod)
        assert row["level_out"] == pytest.approx(field * prod, rel=1e-12), (
            f"slot {row['slot']}: the running level is not field x prod(k) — something other than the "
            "pin is rescaling it"
        )
    excursion = np.array(excursion)
    peak = float(np.max(np.maximum(excursion, 1.0 / excursion)))
    assert peak < 1.0 + 1e-9, (
        f"the running product leaves 1 (peak {peak:.6f}) — a rescale is firing on a chain where no step "
        f"is licensed: {np.round(excursion, 6)}"
    )


def test_the_delivered_level_tracks_the_gdna_FIELD_and_not_any_destinations_own_total():
    """⛔⛔ **TRAPS: a-message-from-the-destinations-belief, AS AN INVARIANCE, AT THE DESTINATIONS WHERE TRAPS: a-message-from-the-destinations-belief ACTUALLY BIT.** The old pin's fixed point
    was half the destination's OWN total observed density, so the delivered level was a function of each
    destination's crowding rather than of the gDNA field. The check needs no tolerance: **scale the whole
    gDNA field by 10x with the RNA held fixed, and every delivered level must scale by 10x exactly.**

    ⭐ **What makes this a discrimination and not a scaling identity**: the RNA is unchanged, so the
    RNA-carrying exons' own totals ``M/E_g`` move by only ~1.1x over the same 10x field, while the field
    itself moves 10x. On those slots the two hypotheses are an order of magnitude apart. ⚠ The BOUNDARIES of
    this fixture carry no RNA, so their own total and the field are the SAME quantity and they
    discriminate nothing — they are asserted on for completeness, and the gate requires that at least one
    genuinely discriminating slot exists."""
    walks = {rho: _relay_walk(_uniform_field_chain(rna=_RNA, rho=rho))[0] for rho in (1.0, 10.0)}
    lo, hi = walks[1.0], walks[10.0]
    read = [r["slot"] for r in lo if r["own_prec"] == 0.0 and np.isfinite(r["level_in"])]
    assert len(read) >= 5, f"only {len(read)} readable slots in the fixture"
    # the slots where the rival hypothesis is far away — without one of these the gate is an identity
    sharp = [s for s in read if hi[s]["M_over_Eg"] / lo[s]["M_over_Eg"] < 2.0]
    assert sharp, (
        "every readable destination's own total scales with the field here, so tracking the field and "
        "tracking a destination's total are not distinguishable on this fixture"
    )
    for s in read:
        assert hi[s]["level_out"] / lo[s]["level_out"] == pytest.approx(10.0, rel=1e-9), (
            f"slot {s}: the delivered level went {lo[s]['level_out']:.6g} → {hi[s]['level_out']:.6g} "
            f"for a 10x gDNA field (a ratio of {hi[s]['level_out'] / lo[s]['level_out']:.4f}); its own "
            f"total moved {hi[s]['M_over_Eg'] / lo[s]['M_over_Eg']:.4f}x"
        )


def test_a_structurally_pure_gdna_destination_IS_told_its_own_measurement():
    """⭐⭐⭐ **AND HERE THE DELIVERED LEVEL *IS* THE DESTINATION'S OWN NUMBER — DELIBERATELY, AND IT IS
    NOT TRAPS: a-message-from-the-destinations-belief.** At a structurally pure-gDNA object the composition is certain by STRUCTURE, so ``M/E_g`` is
    a direct OBSERVATION of the very quantity the message is about, not a belief about it. TRAPS: a-message-from-the-destinations-belief forbids a
    message built from the destination's beliefs and explicitly permits its constants and its data; this
    is the latter. So the pin fires here by case (ii) of its licence and overwrites the relayed level —
    and it must, because an upstream level that disagrees is simply worse evidence about the gDNA density
    at this position than the position's own count.

    ⭐⭐ **THIS IS THE OPERATOR THE CAPTURE LANDSCAPE TRAVELS ON.** An ``intergenic|exon`` BOUNDARY measures
    the gDNA density at its own capture stratum, and the exon behind it has no other way to hear it — a
    G1 BOUNDARY carries ``prec_g = 0`` and so cannot ORIGINATE a level through the fuse (`ROADMAP.md` §1 **reframe-and-level-together**=). Drop
    case (ii) and the off-probe intergenic floor leaks straight through to the exon: measured, and it
    fires `test_gdna_scale_rule.test_capture_step_is_carried_and_the_off_probe_floor_is_not` at 20x and
    200x. A per-capture-class landscape ratio built to do this job explicitly measured inert
    (`region_geometry`'s deleted-landscape note) because this pin already did it.

    The fixture disagrees on purpose: the field on the FIRST region alone is scaled 10x, so the level
    arriving at the BOUNDARY is 10x the BOUNDARY's own measurement, and the BOUNDARY's answer must be its own."""
    BOUNDARY = 1  # chain N E N E N E N E N → the intergenic|exon boundary, structurally pure gDNA
    levels = []
    for rho_first in (1.0, 10.0):
        walk, _, _ = _relay_walk(_uniform_field_chain(rna=_RNA, rho=1.0, rho_first=rho_first))
        row = walk[BOUNDARY]
        assert row["g1"] and not row["lend"], (
            f"slot {BOUNDARY} is not the structural case (g1={row['g1']}, lend={row['lend']})"
        )
        assert row["level_in"] == pytest.approx(rho_first, rel=0.05), (
            f"the upstream field did not reach the BOUNDARY: level_in {row['level_in']:.4g}"
        )
        assert row["level_out"] == pytest.approx(row["M_over_Eg"], rel=1e-12), (
            f"level_out {row['level_out']:.9g} != its own measured density {row['M_over_Eg']:.9g}"
        )
        levels.append(row["level_out"])
    assert levels[0] == pytest.approx(levels[1], rel=1e-12), (
        f"delivered {levels[0]:.9g} then {levels[1]:.9g}: a pure-gDNA object's answer moved when a "
        "DISAGREEING upstream level did, so it is no longer reporting its own measurement"
    )


@pytest.mark.xfail(
    reason=(
        "NOT the nested-exons probe, and not the mass pin: the level reaching the RNA-free interior exon is now the field "
        "EXACTLY (see the gate above it), and the exon still reads 0.914 rather than 1.000. What is "
        "left is psi's uninformative reference holding an evidence-free region off the f_g = 1 vertex by "
        "~0.08 on its own — test_sweep.test_gdna_sweep_factor1_ambig_recovery measures the same "
        "residual from the other side. Strict: the bound is the TRUTH and must not be widened to the "
        "number the solver currently reaches."
    ),
    strict=True,
)
def test_an_rna_free_interior_exon_reads_its_truth():
    """⛔ **WHAT IS STILL SHORT, pinned as a truth rather than as a tolerance.** The middle exon carries
    no RNA at all, so its truth is ``f_g = 1.000``. Both its neighbours are RNA-rich, so the level can
    only reach it by relaying through them — and it now does, exactly (the gate below). The remaining
    0.086 is psi's readout at the vertex, a different mechanism with a different fix."""
    _, final, _ = _relay_walk(_uniform_field_chain(rna=_RNA, rho=1.0))
    assert final.f_g[_MID] > 0.95, f"the RNA-free interior exon reads f_g = {final.f_g[_MID]:.4f}"


def test_the_rna_free_interior_exon_receives_the_field_exactly():
    """⭐⭐ **THE HALF OF THE COST THAT THE LICENCE DOES FIX, separated from the half it does not.** The
    middle exon's two neighbours are RNA-rich and evidence-free, so under the old pin the level was
    inflated 1.96x at the first of them and dragged back down to 0.80x by the time it arrived — a
    multiplicative random walk through psi's half-budget. The delivered level must now be the field, and
    an exact bound is available because nothing here is licensed.

    ⚠ Its ``f_g`` is still 0.914 against a truth of 1.000, and that residual is the xfail above. Keeping
    the two apart is the point: a level that is exact and a readout that is not are different defects."""
    walk, _, _ = _relay_walk(_uniform_field_chain(rna=_RNA, rho=1.0))
    assert walk[_MID]["own_prec"] == 0.0, (
        "the middle exon has own evidence; the fixture has drifted"
    )
    assert walk[_MID]["level_out"] == pytest.approx(1.0, rel=1e-9), (
        f"the level reaching the RNA-free interior exon is {walk[_MID]['level_out']:.9g} against a "
        "field of 1.0"
    )
