"""The backbone's four assertions, and the contract that keeps the shipped policy the shipped policy.

The four: every node ends the two passes holding one message from each neighbour it has; every
delivered row is one row per slot on the solve grid and finite; a slot's population set has at most
three members; and the write-back touches only solvable slots. TRAPS: perturb-every-gate is the
shape of the whole file — each assertion has a matching perturbation test that constructs a policy
committing exactly that defect and asserts the backbone refuses it, because a gate with no firing
perturbation has not been written yet, it has been typed. Byte-identity against the shipped solver
per condition is not gated here: it needs a real chain and a BAM, and it is
``scripts/design/rename_identity.py`` with ``scripts/profiling/sweep_replay.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration import sweep as SW
from rigel.calibration.blocks import SweepCapture, block_slice
from rigel.calibration.message_cache import MessageCache
from rigel.calibration.messages import BlockContext, PsiMessage, Received
from rigel.calibration.messages.silent import SilentPolicy


N = 8


def _ctx(*, free_pos=None, free_neg=None, n_grid=60) -> BlockContext:
    """A minimal BlockContext. Only the fields the assertions read need to be real."""
    ones = np.ones(N)
    fp = np.ones(N, bool) if free_pos is None else np.asarray(free_pos, bool)
    fn = np.zeros(N, bool) if free_neg is None else np.asarray(free_neg, bool)
    return BlockContext(
        eff_gdna=ones * 200.0,
        eff_rna=ones * 200.0,
        sj_count=np.zeros((N, 2)),
        sj_count_lo=np.zeros((N, 2)),
        sj_count_hi=np.zeros((N, 2)),
        route_rate_lo=np.zeros((N, 2)),
        route_rate_hi=np.zeros((N, 2)),
        unspliced_count=np.ones((N, 2)) * 5.0,
        spliced_count=np.zeros((N, 2)),
        left=np.arange(-1, N - 1),
        right=np.append(np.arange(1, N), -1),
        is_boundary=np.arange(N) % 2 == 1,
        is_exon_region=np.zeros(N, bool),
        free_pos=fp,
        free_neg=fn,
        exon_pos=np.zeros(N, bool),
        exon_neg=np.zeros(N, bool),
        boundary_flags=np.zeros(N, np.int64),
        has_own_composition=np.zeros(N, bool),
        belief_fg=ones,
        n_grid=n_grid,
        logodds_window=10.0,
    )


def _counts(msg: PsiMessage, ctx: BlockContext | None = None):
    c = SW.AssertionCounts()
    SW._check_message(msg, ctx if ctx is not None else _ctx(), c)
    return c


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 1 — THE TWO PHASES: every node holds a message from each neighbour
# it has; a real hop must ARRIVE; a missing neighbour is not silence; the kernel sees indices only.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


class _Echo:
    """A policy whose kernel records every hop and writes a gDNA level naming its source into the
    destination's row — so the pass's ORDER and SIDES are observable, and so the solve can be shown the
    two tables."""

    name = "echo"

    def __init__(self):
        self.hops = {False: [], True: []}
        self.held = None

    def library(self, view):
        return None

    def prepare(self, ctx, library):
        return self

    def run_pass(self, received, seq, nbr, terminal, *, backward: bool):
        for i in np.asarray(seq).tolist():
            s = int(nbr[i])
            if s < 0 or terminal[i]:
                continue
            self.hops[backward].append((s, i))
            received.level_gdna.write(i, np.zeros(received.composition.shape[1]), float(s), 0.0)

    def solve(self, from_left, from_right):
        self.held = (from_left, from_right)
        return PsiMessage.silent()


def _passes(pol, ctx, terminal=None):
    K = int(ctx.n_grid)
    prepared = pol.prepare(ctx, None)
    order = list(range(int(ctx.n_slots)))
    fl = SW._pass(order, list(ctx.left), prepared, K, backward=False, terminal=terminal)
    br = SW._pass(order[::-1], list(ctx.right), prepared, K, backward=True, terminal=terminal)
    return fl, br


def test_an_empty_table_is_the_two_states_and_nothing_heard():
    """`Received.empty`: no neighbour anywhere, nothing present, nothing heard; the two states the table
    expresses — NO NEIGHBOUR and SILENCE — partition the nodes that heard nothing, and a written level
    is present with or without a witness but is never a composition."""
    t = Received.empty(N, 7)
    assert t.composition.shape == (N, 7) and t.level_gdna.profile.shape == (N, 7)
    assert not t.has_neighbour.any() and not t.has_composition.any() and not t.heard.any()
    assert np.array_equal(t.no_neighbour, ~t.has_neighbour) and not t.silence.any()
    t.has_neighbour[1:] = True
    assert np.array_equal(t.silence, t.has_neighbour) and not (t.silence & t.no_neighbour).any()
    t.level_rna_pos.write(3, np.zeros(7), 5.0, 100.0)
    t.level_gdna.write(4, np.zeros(7), 6.0, 200.0, 2.0, 8.0)
    assert t.level_rna_pos.present[3] and not t.level_rna_pos.has_witness[3]
    assert t.level_gdna.present[4] and t.level_gdna.has_witness[4]
    assert (t.level_gdna.count[4], t.level_gdna.opportunity[4]) == (6.0, 200.0)
    assert t.heard[3] and t.heard[4] and not t.silence[3] and not t.has_composition.any()
    assert np.array_equal(t.has_level, t.heard)
    for lane in Received.LANES:
        one = Received.empty(2, 3)
        getattr(one, lane).write(1, np.zeros(3), 1.0, 1.0)
        assert one.heard[1] and not one.heard[0], f"a table with only {lane} read as nothing heard"


def test_every_node_holds_what_each_neighbour_it_has_sent_and_an_open_side_is_no_neighbour():
    """After the two passes every interior node holds one level naming its low neighbour and one naming
    its high; ``has_neighbour`` is exactly the chain's links, so the chain's two end nodes read NO
    NEIGHBOUR on the open side — which is not silence — and a level that arrived is not a composition."""
    ctx = _ctx()
    left, right = np.asarray(ctx.left), np.asarray(ctx.right)
    fl, br = _passes(_Echo(), ctx)
    assert np.array_equal(fl.has_neighbour, left >= 0) and np.array_equal(
        br.has_neighbour, right >= 0
    )
    assert np.array_equal(fl.level_gdna.present, left >= 0)
    assert np.array_equal(br.level_gdna.present, right >= 0)
    assert np.array_equal(fl.level_gdna.count[left >= 0], left[left >= 0].astype(float))
    assert np.array_equal(br.level_gdna.count[right >= 0], right[right >= 0].astype(float))
    assert fl.no_neighbour.sum() == 1 and br.no_neighbour.sum() == 1, "one open side each"
    assert not fl.silence.any() and not br.silence.any(), "every real hop carried a level"
    assert not fl.has_composition.any() and not br.has_composition.any(), (
        "a level is not a composition"
    )


def test_the_passes_run_in_chain_order_and_read_one_side_each():
    """The forward pass visits low→high reading each node's LOW neighbour; the backward pass the
    mirror — so what a source holds from its far side is written before it is asked to send."""
    pol = _Echo()
    _passes(pol, _ctx())
    assert pol.hops[False] == [(i - 1, i) for i in range(1, N)]
    assert pol.hops[True] == [(i + 1, i) for i in range(N - 2, -1, -1)]


def test_a_policy_that_sends_nothing_leaves_silence_at_every_node_with_a_neighbour():
    """A pass that writes nothing means every node holds SILENCE from that side — a neighbour and
    nothing present — distinguishable from the open side of the chain."""

    class _Quiet(_Echo):
        def run_pass(self, received, seq, nbr, terminal, *, backward: bool):
            return

    fl, br = _passes(_Quiet(), _ctx())
    for t in (fl, br):
        assert np.array_equal(t.silence, t.has_neighbour) and t.silence.sum() == N - 1
        assert np.array_equal(t.no_neighbour, ~t.has_neighbour) and t.no_neighbour.sum() == 1
        assert not t.heard.any()


def test_every_lane_written_at_a_hop_reaches_the_solve_in_the_same_table():
    """The lanes: a kernel that fills every lane — the composition and the three levels — leaves them in
    the destination's row, and the solve receives the very tables the passes filled. The backbone
    carries; it never reads a lane."""
    ctx = _ctx()
    K = int(ctx.n_grid)
    row = -0.5 * np.linspace(-1.0, 1.0, K) ** 2

    class _Full(_Echo):
        def run_pass(self, received, seq, nbr, terminal, *, backward: bool):
            for i in np.asarray(seq).tolist():
                if nbr[i] < 0 or terminal[i]:
                    continue
                received.composition[i] = row
                received.has_composition[i] = True
                received.level_gdna.write(i, row - 2.0, 0.1, 1.0)
                received.level_rna_pos.write(i, row - 3.0, 0.2, 1.0, 4.0, 4.0)
                received.level_rna_neg.write(i, row - 4.0, 0.3, 1.0)

    pol = _Full()
    fl, br = _passes(pol, ctx)
    pol.solve(fl, br)
    assert pol.held[0] is fl and pol.held[1] is br, (
        "the solve must receive the tables the passes filled"
    )
    for t, nbr in ((fl, np.asarray(ctx.left)), (br, np.asarray(ctx.right))):
        has = nbr >= 0
        assert np.array_equal(t.has_composition, has) and np.array_equal(t.heard, has)
        assert np.array_equal(t.composition[has], np.tile(row, (int(has.sum()), 1)))
        for lane, shift in zip(Received.LANES, (2.0, 3.0, 4.0)):
            lv = getattr(t, lane)
            assert np.array_equal(lv.present, has)
            assert np.array_equal(lv.profile[has], np.tile(row - shift, (int(has.sum()), 1)))
        assert (
            np.array_equal(t.level_rna_pos.has_witness, has) and not t.level_gdna.has_witness.any()
        )


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 2 — every delivered ROW is one row per slot on the solve grid, and finite.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_lambda_rows_are_checked_for_shape_and_finiteness():
    """The λ-row channel: one row per slot on the solve grid is accepted and counted finite; a row array
    of another shape is REFUSED outright; a non-finite row is COUNTED and, with no waiver, raises."""
    ctx = _ctx()
    K = int(ctx.n_grid)
    ok = _counts(PsiMessage(lam_rows=np.zeros((N, K))), ctx)
    assert ok["lam_rows_finite"] == {"violations": 0, "eligible": N}
    with pytest.raises(ValueError, match="lam_rows has shape"):
        _counts(PsiMessage(lam_rows=np.zeros((N + 1, K))), ctx)
    with pytest.raises(ValueError, match="lam_rows has shape"):
        _counts(PsiMessage(lam_rows=np.zeros(N)), ctx)
    bad = np.zeros((N, K))
    bad[2, 0] = np.nan
    with pytest.raises(AssertionError, match="lam_rows_finite"):
        _counts(PsiMessage(lam_rows=bad), ctx)


def test_a_waiver_is_never_silent():
    """Every waived assertion carries a written reason, so a reader learns the defect rather than the
    exemption. An empty reason would be a widened predicate wearing a waiver's clothes."""
    for name, why in SW._KNOWN_VIOLATIONS.items():
        assert len(why) > 80, f"{name}'s waiver does not say what the defect is"


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 3 — |T| <= 3.  AXIOM 0, made executable.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_population_set_is_at_most_three_because_it_is_a_function_of_two_bits():
    """``T(slot) = {gDNA} u {RNA+ if free_pos} u {RNA- if free_neg}``, so ``|T| = 1 + free_pos +
    free_neg`` and it is in ``{1, 2, 3}`` for every possible input. There are three populations and
    there is no fourth: "mature" and "nascent" are not species, and RNA inside an intron is RNA that
    has not spliced at that position. Structural rather than something to remember, which is the
    point."""
    for fp in (True, False):
        for fn in (True, False):
            ctx = _ctx(free_pos=np.full(N, fp), free_neg=np.full(N, fn))
            pop = ctx.population_size()
            assert set(np.unique(pop)) <= {1, 2, 3}
            assert np.all(pop == 1 + int(fp) + int(fn))


def test_PERTURBATION_a_fourth_population_is_REFUSED():
    """Axiom 0's tell, executable: a population set with more than three members. A derivation that
    opens with ``{gDNA, nascent+, nascent-, mature+, mature-}`` produces a wrong table every
    time."""
    counts = SW.AssertionCounts()
    with pytest.raises(AssertionError, match="population_at_most_three"):
        counts.note("population_at_most_three", np.array([4, 4, 5]) > 3, np.ones(3, bool))


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 4 — the write-back touches only `solvable` slots.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_PERTURBATION_a_writeback_outside_solvable_is_REFUSED():
    """Getting this wrong reads as a byte-identity failure of ``max|delta| = 1.0``
    (TRAPS: byte-identity-gate): a replay compares the solve's raw output against the shipped
    belief, and the two differ by exactly this mask. Reproducing a pipeline stage means reproducing
    its write-back.

    A locked slot — one with no admissible RNA strand — is never solved and keeps its
    signature-binary init, because RNA cannot cross a gene boundary so its unspliced mass is purely
    gDNA."""
    untouched = np.array([False, True, True, False])
    changed = np.array([False, False, True, False])
    counts = SW.AssertionCounts()
    with pytest.raises(AssertionError, match="writeback_only_solvable"):
        counts.note("writeback_only_solvable", untouched & changed, untouched)


def test_a_writeback_confined_to_solvable_is_accepted():
    counts = SW.AssertionCounts()
    untouched = np.array([False, True, True, False])
    counts.note("writeback_only_solvable", untouched & np.zeros(4, bool), untouched)
    assert counts["writeback_only_solvable"] == {"violations": 0, "eligible": 2}


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE CONTRACT — what keeps the shipped policy the shipped policy.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_message_policy_is_a_config_decision_and_defaults_to_transfer():
    """The largest behaviour switch in the tool, and it must be a written decision. Which policy
    ships can never be inherited from a function default that an edit could silently change: the
    config names it (``message_policy``, ``"transfer"``), ``calibrate`` reads it, and both the
    shipped policy and the measured floor (``"silent"``) are reachable from the one call site."""
    import importlib
    import inspect

    from rigel.config import CalibrationConfig  # noqa: PLC0415

    # the package re-exports the function under the module's name, so ``import … as`` would bind the
    # FUNCTION and this gate would read one function's body; the module is what reads the switch
    _c = importlib.import_module("rigel.calibration.calibrate")
    assert CalibrationConfig().message_policy == "transfer"
    src = inspect.getsource(_c)
    assert "config.message_policy" in src, (
        "calibrate no longer reads the switch — whichever policy it now hard-codes, the config option is "
        "lying to anyone who sets it."
    )
    assert "TransferPolicy(" in src and "SilentPolicy()" in src, (
        "both arms must be reachable from the one call site; a switch with one arm is not a switch."
    )


def test_solve_chains_parameter_default_is_silent_and_sends_nothing():
    """``SilentPolicy`` is ``solve_chain``'s parameter default (the shipped config installs the
    transfer policy), and it is the MEASURED floor every policy is judged against: win on unstranded
    data, minimal harm on stranded data, never pooled."""
    ctx = _ctx()
    prepared = SilentPolicy().prepare(ctx, None)
    nothing = Received.empty(N, int(ctx.n_grid))
    order = np.arange(N, dtype=np.int64)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    prepared.run_pass(nothing, order, left, np.zeros(N, bool), backward=False)
    prepared.run_pass(nothing, order[::-1], right, np.zeros(N, bool), backward=True)
    assert not nothing.heard.any(), "a silent policy must send nothing"
    assert prepared.solve(nothing, nothing).is_silent


def test_the_backbone_does_not_know_what_a_message_is_about():
    """The structural claim of the split, as a test. The backbone owns the shape of the solve and
    the four assertions; every message-composition choice is a policy. If one of these concepts
    reappears in ``sweep.py``, an operator has leaked back into the backbone and the next reader can
    no longer hold the working system in their head.

    It checks identifiers, from the AST, and not the file's text. Grepping the source would match a
    docstring naming these very words in order to say they are absent, and a test that passes for
    that reason is vacuous.

    ``capture`` has one licensed occurrence and it is not the biology: ``_capture`` is the
    diagnostics hook every instrument passes by keyword. Hybrid capture — the thing the message
    layer argues about — appears nowhere."""
    import ast
    import inspect

    tree = ast.parse(inspect.getsource(SW))
    ident: set[str] = set()
    for region in ast.walk(tree):
        for attr in ("id", "name", "arg", "attr"):
            v = getattr(region, attr, None)
            if isinstance(v, str):
                ident.add(v.lower())
    banned = (
        "SPLICE IN",
        "reframe",
        "pin",
        "enrichment",
        "SPLICE OUT",
        "flank",
        "damp",
        "mismatch",
        "level",
        "lane",
    )
    leaked = {w: sorted(i for i in ident if w in i) for w in banned}
    assert not any(leaked.values()), (
        f"policy concepts leaked into the backbone's identifiers: {leaked}"
    )
    cap = sorted(i for i in ident if "capture" in i)
    assert cap == ["_capture", "sweepcapture"], (
        "the only licensed 'capture' words are the diagnostics hook — the parameter — and its record "
        f"(`blocks.SweepCapture`). Found {cap}"
    )


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE CUBE CHANNEL (the both-stranded locus): a (K, K_t) row per AMBIG slot, final solve only
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_cube_channel_is_checked_per_ambig_slot_and_shape():
    """`_check_message`: a delivered row at an AMBIG slot with ``(K,)`` profiles passes; one at a
    single-strand slot is REFUSED (a cube exists only where both strands are live); a profile that is not
    ``(n_grid,)`` is REFUSED; a non-finite profile is REFUSED by the backbone's assertion."""
    from _psi_reference import cube_rows_of

    ctx = _ctx(free_pos=np.ones(N, bool), free_neg=np.array([i % 2 == 0 for i in range(N)]))
    K = int(ctx.n_grid)

    def rows(slot, prof, k=K):
        return cube_rows_of({slot: (prof, None, 100.0, 50.0, 0.5)}, np.linspace(-10.0, 10.0, k))

    ok = _counts(PsiMessage(cube_rows=rows(0, np.zeros(K))), ctx)
    assert ok["cube_rows_finite"] == {"violations": 0, "eligible": 1}
    with pytest.raises(ValueError, match="not an AMBIG slot"):
        _counts(PsiMessage(cube_rows=rows(1, np.zeros(K))), ctx)
    with pytest.raises(ValueError, match="solve grid"):
        _counts(PsiMessage(cube_rows=rows(0, np.zeros(K + 1), K + 1)), ctx)
    with pytest.raises(AssertionError, match="cube_rows_finite"):
        _counts(PsiMessage(cube_rows=rows(0, np.full(K, np.nan))), ctx)


def test_the_solvers_cube_is_inert_when_absent_and_walls_the_tilt_when_present():
    """`_solve_regions_logodds_all(cube_rows=...)`: ``None`` and an empty table are byte-identical to the
    path without the argument; a wall against low ``f_+`` at one AMBIG slot raises that slot's ``f_pos`` and
    leaves every other slot byte-identical."""
    from _psi_reference import cube_rows_of

    from rigel.calibration.simplex_logodds import CubeRows, _solve_regions_logodds_all

    m, K = 4, 40
    u_pos = np.array([50.0, 50.0, 50.0, 50.0])
    u_neg = np.array([50.0, 50.0, 50.0, 50.0])
    ap = np.ones(m, bool)
    an = np.array([True, True, False, True])
    kw = dict(kappa=0.99, od_g=0.0, od_r=0.0, n_grid=K, L=10.0)
    base = _solve_regions_logodds_all(u_pos, u_neg, ap, an, u_pos + u_neg, np.zeros(m), **kw)
    u = np.linspace(-10.0, 10.0, K)
    for empty in (None, CubeRows.blank(0, u)):
        again = _solve_regions_logodds_all(
            u_pos, u_neg, ap, an, u_pos + u_neg, np.zeros(m), cube_rows=empty, **kw
        )
        assert np.array_equal(again.gdna_frac, base.gdna_frac)
        assert np.array_equal(again.rna_pos_frac, base.rna_pos_frac)
    # "at least 30 % RNA+": a wall below the density that share implies, n = 100, a_r = 100, ρ_ref = 1
    floor = np.where(u < np.log(0.3), -50.0, 0.0)
    wall = cube_rows_of({1: (floor, None, 100.0, 100.0, 1.0)}, u)
    walled = _solve_regions_logodds_all(
        u_pos, u_neg, ap, an, u_pos + u_neg, np.zeros(m), cube_rows=wall, **kw
    )
    assert walled.rna_pos_frac[1] > base.rna_pos_frac[1] + 0.05
    for i in (0, 2, 3):
        assert (
            walled.gdna_frac[i] == base.gdna_frac[i]
            and walled.rna_pos_frac[i] == base.rna_pos_frac[i]
        )


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE STRUCTURAL RULE — a TERMINAL receives nothing, so the chain breaks at it.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_a_terminal_receives_nothing_and_the_kernel_is_never_asked_for_the_hop_into_it():
    """The boundary condition the locus solve stands on: a node marked terminal holds SILENCE from a
    side it has a neighbour on, and the policy's kernel is never called with it as the destination —
    so what a policy WOULD deliver there cannot exist. PERTURBATION: the same kernel with no terminal
    marked delivers its message, which is what proves the mask does the work. What the terminal SENDS
    is untouched: its neighbour still receives from it."""
    ctx = _ctx()
    terminal = [False] * N
    terminal[4] = True
    pol = _Echo()
    held, _br = _passes(pol, ctx, terminal=terminal)
    assert held.has_neighbour[4] and held.silence[4], (
        "a terminal holds SILENCE: a neighbour, nothing present"
    )
    assert (3, 4) not in pol.hops[False], "the kernel was asked for the hop into the terminal"
    assert held.level_gdna.count[5] == 4.0, (
        "the terminal's own sending was blocked; only receiving is"
    )
    assert held.no_neighbour[0], "a terminal rule must not turn an open side into silence"
    loud = _Echo()
    unmasked, _br = _passes(loud, ctx)
    assert unmasked.level_gdna.count[4] == 3.0 and (3, 4) in loud.hops[False]


def test_the_terminal_predicate_is_the_solve_gates_lock_on_a_region():
    """One predicate, two names must not appear: the slots the backbone never delivers into are exactly
    the REGIONS `g1_locked` locks — no admissible RNA strand — read off the source, so a reader cannot
    find a second definition of "terminal" in the file."""
    import inspect
    import re

    src = inspect.getsource(SW)
    assert re.search(r"terminal ?= ?.*is_region & g1_locked\(fp, fn\)", src), "the predicate moved"
    assert re.search(r"_pass\(.*terminal=", src), "the passes are no longer told the terminals"


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE LOCUS BLOCKS — the sweep solved a block at a time is the sweep, for every block size.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def _six(belief):
    return {
        f: np.asarray(getattr(belief, f))
        for f in ("f_g", "f_pos", "f_neg", "var_gdna", "has_composition")
    }


def test_the_block_solve_is_the_chain_solve_for_every_block_size(sweep_inputs):
    """The property the whole decomposition stands on, on the real toy chain through the shipped
    policy: ``solve_chain`` with the chain as one block, one locus per block, and every block size in
    between gives the same belief to the bit and the same ``has_composition`` predicate — and the diagnostic
    capture gathers to the same per-slot arrays. Not vacuous: the toy has several terminals, so
    ``block_slots=1`` makes more than one block."""
    from _transfer_harness import _full_policy
    from rigel.calibration.region_chain import REGION, locus_blocks
    from rigel.calibration.region_geometry import g1_locked

    chain, statics = sweep_inputs["args"][0], sweep_inputs["args"][1]
    terminal = (np.asarray(chain.kind) == REGION) & g1_locked(statics.free_pos, statics.free_neg)
    assert len(locus_blocks(chain, terminal, 1)) > 2, "the toy has fewer than three loci"
    policy = _full_policy(sweep_inputs)[0]
    caps = {}

    kw = dict(sweep_inputs["kw"])
    kw.pop("block_slots", None)  # the capture carries calibrate's; this gate sets its own

    def run(block_slots):
        cap = SweepCapture()
        out = SW.solve_chain(
            *sweep_inputs["args"], **kw, policy=policy, block_slots=block_slots, _capture=cap
        )
        caps[block_slots] = cap
        return _six(out)

    whole = run(None)
    for bs in (1, 2, 5, 17, 10_000):
        got = run(bs)
        for f in whole:
            assert np.array_equal(got[f], whole[f]), f"block_slots={bs}: {f} differs"
        for key in ("f_g", "fg_loc", "tau_lam", "solvable", "mass_global", "count"):
            assert np.array_equal(getattr(caps[bs], key), getattr(caps[None], key)), (
                f"block_slots={bs}: capture {key} differs"
            )
        a, b = caps[bs].lam_rows, caps[None].lam_rows
        assert (a is None) == (b is None) and (a is None or np.array_equal(a, b))
        assert caps[bs].backbone_assertions == caps[None].backbone_assertions
        # what each node HEARD is the same; at a block's first slot — a terminal, which hears nothing
        # by the structural rule — an open side and SILENCE are the same hearing
        for side in ("from_left", "from_right"):
            heard = [getattr(caps[k], side).heard for k in (bs, None)]
            assert np.array_equal(heard[0], heard[1]), (
                f"block_slots={bs}: {side} differs in what was heard"
            )
            for i in np.flatnonzero(terminal):
                assert not any(getattr(caps[k], side).heard[i] for k in (bs, None))


def test_a_block_view_rebases_the_links_and_slices_every_per_slot_array(sweep_inputs):
    """`blocks.block_slice`: a neighbour outside the block is no neighbour; every per-slot array is the chain's
    slice; a 2-D bank keeps its columns; ``n_slots`` follows."""
    chain, statics, geometry, belief, _ra = sweep_inputs["args"]
    n = int(chain.n_slots)
    sl = slice(2, min(9, n))
    c = block_slice(chain, sl)
    assert c.n_slots == sl.stop - sl.start and np.array_equal(c.kind, np.asarray(chain.kind)[sl])
    left = np.asarray(chain.left)[sl] - sl.start
    assert np.array_equal(c.left, np.where((left >= 0) & (left < c.n_slots), left, -1))
    assert c.left[0] == -1 and c.right[-1] == -1
    g = block_slice(geometry, sl)
    assert g.n_slots == c.n_slots and g.unspliced_count.shape == (c.n_slots, 2)
    assert np.array_equal(g.eff_gdna, np.asarray(geometry.eff_gdna)[sl])
    b = block_slice(belief, sl)
    assert np.array_equal(b.f_g, np.asarray(belief.f_g)[sl])
    assert block_slice(statics, sl).boundary_flags.shape == (c.n_slots,)


def test_the_checks_count_only_the_owned_slots():
    """`_check_message(n_owned=...)`: the read-ahead terminal at the end of a block is not counted as
    eligible, and a row array is still required to cover every slot the policy saw."""
    ctx = _ctx()
    K = int(ctx.n_grid)
    c = SW.AssertionCounts()
    SW._check_message(PsiMessage(lam_rows=np.zeros((N, K))), ctx, c, n_owned=N - 1)
    assert c["lam_rows_finite"] == {"violations": 0, "eligible": N - 1}
    assert c["population_at_most_three"]["eligible"] == N - 1
    with pytest.raises(ValueError, match="lam_rows has shape"):
        SW._check_message(PsiMessage(lam_rows=np.zeros((N - 1, K))), ctx, c, n_owned=N - 1)
    merged = SW.AssertionCounts()
    merged.absorb(c)
    merged.absorb(c)
    assert merged["lam_rows_finite"] == {"violations": 0, "eligible": 2 * (N - 1)}


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE MESSAGE MEMO — the message layer is refit-invariant given the grid, so the refit sweeps share it.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def _cache_kw(sweep_inputs):
    kw = dict(sweep_inputs["kw"])
    kw.pop("block_slots", None)
    kw.pop("message_cache", None)
    return kw


def test_a_cache_hit_reproduces_the_uncached_sweep_to_the_bit_and_skips_the_layer(sweep_inputs):
    """Two sweeps on identical inputs through one cache: the second hits every block, never calls the
    policy's `prepare`, and returns the same belief and the same ``has_composition`` as the first and as a
    sweep with no cache at all. Not vacuous: the toy delivers rows, so a stale or empty hit would move
    the numbers."""
    from _transfer_harness import _full_policy

    policy = _full_policy(sweep_inputs)[0]
    kw = _cache_kw(sweep_inputs)
    plain = SW.solve_chain(*sweep_inputs["args"], **kw, policy=policy, block_slots=5)
    cache = MessageCache()
    first = SW.solve_chain(
        *sweep_inputs["args"], **kw, policy=policy, block_slots=5, message_cache=cache
    )
    assert cache.misses > 0 and cache.hits == 0
    n_blocks = cache.misses
    calls = []
    orig = type(policy).prepare

    def spy(self, ctx, library):
        calls.append(ctx.n_slots)
        return orig(self, ctx, library)

    type(policy).prepare = spy
    try:
        second = SW.solve_chain(
            *sweep_inputs["args"], **kw, policy=policy, block_slots=5, message_cache=cache
        )
    finally:
        type(policy).prepare = orig
    assert not calls, "a hit must not prepare the policy"
    assert cache.hits == n_blocks and cache.misses == n_blocks
    for out in (first, second):
        for f in ("f_g", "f_pos", "f_neg", "var_gdna", "has_composition"):
            assert np.array_equal(np.asarray(getattr(out, f)), np.asarray(getattr(plain, f))), f
    assert cache.nbytes > 0


def test_PERTURBATION_the_cache_misses_when_any_input_the_message_layer_reads_changes(sweep_inputs):
    """The key is a digest of EVERY input the message layer reads, so it is safe by construction: a
    changed belief, liveness bit, factory row, observation, library or grid must miss — and a cache that
    hit on any of them would deliver another sweep's messages as this one's."""
    import dataclasses as _dc

    from _transfer_harness import _full_policy

    policy, rows, _g, _w = _full_policy(sweep_inputs)
    kw = dict(
        _cache_kw(sweep_inputs), intron_prior=rows
    )  # the live rows: the layer has something to read
    chain, statics, geometry, belief, ra = sweep_inputs["args"]
    cache = MessageCache()
    SW.solve_chain(chain, statics, geometry, belief, ra, **kw, policy=policy, message_cache=cache)
    base_misses = cache.misses

    def misses_after(pol=policy, **over):
        before = cache.misses
        a = over.pop("args", (chain, statics, geometry, belief, ra))
        SW.solve_chain(*a, **{**kw, **over}, policy=pol, message_cache=cache)
        return cache.misses - before

    assert misses_after() == 0, "identical inputs must hit"
    # the incoming belief (the variance freeze of every own strand profile)
    b2 = _dc.replace(belief, f_g=np.asarray(belief.f_g) * 0.999 + 0.0005)
    assert misses_after(args=(chain, statics, geometry, b2, ra)) == base_misses
    # each perturbation touches ONLY its channel, so the gate isolates that channel of the digest:
    # a constant shift on an intron's live row changes the rows' bits and nothing downstream of them
    # (liveness is a curvature, the claim is max-normalised) …
    from _transfer_harness import _expected_pairs

    intron = int(_expected_pairs(sweep_inputs)[0][0][1])
    r2 = np.asarray(rows).copy()
    r2[intron] += 1e-3
    assert misses_after(intron_prior=r2) >= 1
    # … and an intron's gDNA opportunity is read by the passes alone — by neither the library (which
    # sums intergenic and single-strand-exon opportunities) nor any liveness bit
    eg = np.asarray(geometry.eff_gdna).copy()
    eg[intron] *= 1.5
    g2 = _dc.replace(geometry, eff_gdna=eg)
    assert misses_after(args=(chain, statics, g2, belief, ra)) >= 1
    # the grid: the bracket (the rows keep their shape; a wider bracket re-reads them, so a hit would
    # serve messages laid on another lattice)
    assert misses_after(logodds_window=float(kw["logodds_window"]) + 1.0) == base_misses
    # the library and the policy: another strand model changes every lane's coordinate and every
    # own strand claim — the sweep's own ``rna_sense_frac`` is ψ's input, not the message layer's,
    # so the policy is what carries the strand into the key
    from _transfer_harness import _strand_of
    from rigel.calibration.messages.transfer import TransferPolicy

    kappa, od_g, od_r = _strand_of(sweep_inputs)
    other = TransferPolicy(strand=(kappa * 0.9 + 0.05, od_g, od_r))
    assert misses_after(pol=other) == base_misses
    assert misses_after(rna_sense_frac=float(kw["rna_sense_frac"]) * 0.9 + 0.05) == 0, (
        "the sweep's own strand parameter is ψ's, not the message layer's: a hit is correct"
    )


def test_a_diagnostic_capture_always_runs_the_full_layer(sweep_inputs):
    """An instrument's capture reads the held lists, so with ``_capture`` the layer runs even on a cache
    that would hit — and what it delivers equals the cache's, so the two paths cannot drift."""
    from _transfer_harness import _full_policy

    policy = _full_policy(sweep_inputs)[0]
    kw = _cache_kw(sweep_inputs)
    cache = MessageCache()
    SW.solve_chain(*sweep_inputs["args"], **kw, policy=policy, message_cache=cache)
    cap = SweepCapture()
    hits_before = cache.hits
    out = SW.solve_chain(
        *sweep_inputs["args"], **kw, policy=policy, message_cache=cache, _capture=cap
    )
    assert cache.hits == hits_before, "a captured sweep must not be served from the cache"
    assert cap.from_left is not None and out.has_composition is not None
