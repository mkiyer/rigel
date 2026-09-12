"""The backbone's four assertions, and the contract that keeps the shipped policy the shipped policy.

The four: every node ends the two passes holding one message from each neighbour it has; every
delivered row is one row per slot on the solve grid and finite; a slot's population set has at most
three members; and the write-back touches only solvable slots. TRAPS: perturb-every-gate is the
shape of the whole file — each assertion has a matching perturbation test that constructs a policy
committing exactly that defect and asserts the backbone refuses it, because a gate with no firing
perturbation has not been written yet, it has been typed. Byte-identity against the shipped solver
per condition is not gated here: it needs a real chain and a BAM, and it is
``scripts/design/backbone_parity.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration import sweep as SW
from rigel.calibration.messages import NO_NEIGHBOUR, SILENCE, Message, PsiMessage, StepContext
from rigel.calibration.messages.silent import SilentPolicy


N = 8


def _ctx(*, free_pos=None, free_neg=None, n_grid=60) -> StepContext:
    """A minimal StepContext. Only the fields the assertions read need to be real."""
    ones = np.ones(N)
    fp = np.ones(N, bool) if free_pos is None else np.asarray(free_pos, bool)
    fn = np.zeros(N, bool) if free_neg is None else np.asarray(free_neg, bool)
    return StepContext(
        eff_gdna=ones * 200.0,
        eff_rna=ones * 200.0,
        sj_count=np.zeros((N, 2)),
        sj_count_lo=np.zeros((N, 2)),
        sj_count_hi=np.zeros((N, 2)),
        route_rate_lo=np.zeros((N, 2)),
        route_rate_hi=np.zeros((N, 2)),
        unspliced_count=np.ones((N, 2)) * 5.0,
        n_slot=ones * 10.0,
        spliced_slot=np.zeros(N),
        left=np.arange(-1, N - 1),
        right=np.append(np.arange(1, N), -1),
        is_boundary=np.arange(N) % 2 == 1,
        is_exon_region=np.zeros(N, bool),
        free_pos=fp,
        free_neg=fn,
        exon_pos=np.zeros(N, bool),
        exon_neg=np.zeros(N, bool),
        boundary_flags=np.zeros(N, np.int64),
        own_live=np.zeros(N, bool),
        belief_fg=ones,
        n_grid=n_grid,
        logodds_window=10.0,
    )


def _counts(msg: PsiMessage, ctx: StepContext | None = None):
    c = SW.AssertionCounts()
    SW._check_message(msg, ctx if ctx is not None else _ctx(), c)
    return c


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 1 — THE TWO PHASES: every node holds a message from each neighbour
# it has; a real hop must ARRIVE; a missing neighbour is not silence; the kernel sees indices only.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


class _Echo:
    """A policy whose kernel records every hop and returns a message naming its source — so the pass's
    ORDER and SIDES are observable, and so the solve can be shown the two held lists."""

    name = "echo"

    def __init__(self):
        self.hops = {False: [], True: []}
        self.held = None

    def library(self, view):
        return None

    def prepare(self, ctx, library):
        return self

    def propagate(self, *, backward: bool):
        def receive(s, i):
            self.hops[backward].append((int(s), int(i)))
            return Message(level_gdna=(float(s), 0.0))

        return receive

    def solve(self, from_left, from_right):
        self.held = (list(from_left), list(from_right))
        return PsiMessage.silent()


def test_every_node_holds_a_message_from_each_neighbour_it_has():
    """After the two passes every interior node holds two messages, one
    naming its low neighbour and one its high; the chain's two end nodes hold ONE and ``NO_NEIGHBOUR``
    on the open side — which is not a message and not SILENCE."""
    ctx = _ctx()
    left, right = list(ctx.left), list(ctx.right)
    fl = SW._pass(list(range(ctx.n_slots)), left, _Echo().prepare(ctx, None), backward=False)
    br = SW._pass(list(range(ctx.n_slots))[::-1], right, _Echo().prepare(ctx, None), backward=True)
    for i in range(N):
        if left[i] >= 0:
            assert fl[i].level_gdna[0] == float(left[i]), f"slot {i} holds the wrong low neighbour"
        else:
            assert fl[i] is NO_NEIGHBOUR and fl[i] is not SILENCE
        if right[i] >= 0:
            assert br[i].level_gdna[0] == float(right[i])
        else:
            assert br[i] is NO_NEIGHBOUR
    assert sum(m is None for m in fl) == 1 and sum(m is None for m in br) == 1, "one open side each"


def test_the_passes_run_in_chain_order_and_read_one_side_each():
    """The forward pass visits low→high reading each node's LOW neighbour; the backward pass the
    mirror — so what a source holds from its far side is written before it is asked to send."""
    ctx = _ctx()
    pol = _Echo()
    SW._pass(list(range(ctx.n_slots)), list(ctx.left), pol.prepare(ctx, None), backward=False)
    SW._pass(list(range(ctx.n_slots))[::-1], list(ctx.right), pol.prepare(ctx, None), backward=True)
    assert pol.hops[False] == [(i - 1, i) for i in range(1, N)]
    assert pol.hops[True] == [(i + 1, i) for i in range(N - 2, -1, -1)]


def test_PERTURBATION_a_kernel_that_leaves_a_real_hop_unspoken_is_REFUSED():
    """A hop that carries nothing must still arrive as SILENCE. A kernel returning ``None`` for a node
    that HAS a neighbour is the one thing the pass refuses, because the solve could not then tell
    "nothing to say" from "never spoken to"."""

    class _Mute(_Echo):
        def propagate(self, *, backward: bool):
            return lambda s, i: None

    ctx = _ctx()
    with pytest.raises(AssertionError, match="must still arrive as SILENCE"):
        SW._pass(
            list(range(ctx.n_slots)), list(ctx.left), _Mute().prepare(ctx, None), backward=False
        )


def test_a_policy_that_sends_nothing_leaves_silence_at_every_node_with_a_neighbour():
    """``propagate`` returning no kernel means every node holds SILENCE from that side — delivered,
    distinguishable from the open side of the chain."""
    ctx = _ctx()

    class _Quiet(_Echo):
        def propagate(self, *, backward: bool):
            return None

    fl = SW._pass(
        list(range(ctx.n_slots)), list(ctx.left), _Quiet().prepare(ctx, None), backward=False
    )
    assert fl[0] is NO_NEIGHBOUR and all(m is SILENCE for m in fl[1:])
    assert SILENCE.is_silent and Message(level_gdna=(0.0, 1.0)).is_silent is False


def test_every_lane_of_a_message_survives_the_passes_to_the_solve():
    """The lanes: a kernel that fills every lane — the composition profile and the three level
    claims — hands them to the solve untouched, and a message with any one lane is not silent. The
    backbone carries; it never reads a lane."""
    ctx = _ctx()
    full = Message(
        composition=np.zeros(3),
        level_gdna=(-2.0, 0.1),
        level_rna_pos=(-3.0, 0.2),
        level_rna_neg=(-4.0, 0.3),
    )

    class _Full(_Echo):
        def propagate(self, *, backward: bool):
            return lambda s, i: full

    pol = _Full()
    prepared = pol.prepare(ctx, None)
    fl = SW._pass(list(range(ctx.n_slots)), list(ctx.left), prepared, backward=False)
    br = SW._pass(list(range(ctx.n_slots))[::-1], list(ctx.right), prepared, backward=True)
    prepared.solve(fl, br)
    got_l, got_r = pol.held
    for i in range(N):
        if list(ctx.left)[i] >= 0:
            assert got_l[i] is full
        if list(ctx.right)[i] >= 0:
            assert got_r[i] is full
    assert not full.is_silent and Message().is_silent
    for lane in Message.LANES:
        one = Message(**{lane: np.zeros(2) if lane == "composition" else (0.0, 1.0)})
        assert not one.is_silent, f"a message with only {lane} read as silence"


def test_the_solve_receives_the_two_held_lists_at_the_recipient():
    """Phase 2's inputs are the two held lists indexed AT THE RECIPIENT, straight from the passes."""
    ctx = _ctx()
    pol = _Echo()
    prepared = pol.prepare(ctx, None)
    fl = SW._pass(list(range(ctx.n_slots)), list(ctx.left), prepared, backward=False)
    br = SW._pass(list(range(ctx.n_slots))[::-1], list(ctx.right), prepared, backward=True)
    prepared.solve(fl, br)
    assert pol.held == (fl, br)


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


def test_message_propagation_is_a_config_switch_and_defaults_ON():
    """The largest behaviour switch in the tool, and it must be a written decision. Which policy
    ships can never be inherited from a function default that an edit could silently change: the
    config names it (``message_policy``, ``"transfer"``), ``calibrate`` reads it, and both the
    shipped policy and the measured floor are reachable from the one call site."""
    import inspect

    import rigel.calibration.calibrate as _c  # noqa: PLC0415
    from rigel.config import CalibrationConfig  # noqa: PLC0415

    assert CalibrationConfig().message_propagation is True
    assert CalibrationConfig().message_policy == "transfer"
    src = inspect.getsource(_c)
    assert "config.message_propagation" in src and "config.message_policy" in src, (
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
    prepared = SilentPolicy().prepare(_ctx(), None)
    assert prepared.propagate(backward=False) is None, "a silent policy must send nothing at all"
    assert prepared.propagate(backward=True) is None
    assert prepared.solve([SILENCE] * N, [SILENCE] * N).is_silent


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
    assert cap == ["_capture"], (
        f"the only licensed 'capture' is the diagnostics hook, the parameter. Found {cap}"
    )


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE CUBE CHANNEL (the both-stranded locus): a (K, K_t) row per AMBIG slot, final solve only
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_cube_channel_is_checked_per_ambig_slot_and_shape():
    """`_check_message`: a cube row at an AMBIG slot of shape ``(K, K_t)`` passes; a row at a
    single-strand slot is REFUSED (a cube exists only where both strands are live); a row that is not
    ``(n_grid, ·)`` is REFUSED; a non-finite row is REFUSED by the backbone's assertion."""
    ctx = _ctx(free_pos=np.ones(N, bool), free_neg=np.array([i % 2 == 0 for i in range(N)]))
    K = int(ctx.n_grid)
    ok = _counts(PsiMessage(cube_rows={0: np.zeros((K, 7))}), ctx)
    assert ok["cube_rows_finite"] == {"violations": 0, "eligible": 1}
    with pytest.raises(ValueError, match="not an AMBIG slot"):
        _counts(PsiMessage(cube_rows={1: np.zeros((K, 7))}), ctx)
    with pytest.raises(ValueError, match="expected"):
        _counts(PsiMessage(cube_rows={0: np.zeros((K + 1, 7))}), ctx)
    with pytest.raises(AssertionError, match="cube_rows_finite"):
        _counts(PsiMessage(cube_rows={0: np.full((K, 7), np.nan)}), ctx)


def test_the_solvers_cube_is_inert_when_absent_and_walls_the_tilt_when_present():
    """`_solve_regions_logodds_all(cube_rows=...)`: ``None`` and ``{}`` are byte-identical to the path
    without the argument; a wall against low ``f_+`` at one AMBIG slot raises that slot's ``f_pos`` and
    leaves every other slot byte-identical."""
    from rigel.calibration.simplex_logodds import _solve_regions_logodds_all, _tilt_grid

    m, K, Kt = 4, 40, 24
    u_pos = np.array([50.0, 50.0, 50.0, 50.0])
    u_neg = np.array([50.0, 50.0, 50.0, 50.0])
    ap = np.ones(m, bool)
    an = np.array([True, True, False, True])
    kw = dict(kappa=0.99, od_g=0.0, od_r=0.0, n_grid=K, L=10.0, n_tilt=Kt)
    base = _solve_regions_logodds_all(u_pos, u_neg, ap, an, u_pos + u_neg, np.zeros(m), **kw)
    for empty in (None, {}):
        again = _solve_regions_logodds_all(
            u_pos, u_neg, ap, an, u_pos + u_neg, np.zeros(m), cube_rows=empty, **kw
        )
        assert np.array_equal(again.gdna_frac, base.gdna_frac)
        assert np.array_equal(again.rna_pos_frac, base.rna_pos_frac)
    lam = np.linspace(-10.0, 10.0, K)
    tau = np.sin(_tilt_grid(Kt))
    f_pos = (1.0 - 1.0 / (1.0 + np.exp(-lam)))[:, None] * (1.0 + tau)[None, :] / 2.0
    wall = np.where(f_pos < 0.3, -50.0, 0.0)  # "at least 30 % RNA+"
    walled = _solve_regions_logodds_all(
        u_pos, u_neg, ap, an, u_pos + u_neg, np.zeros(m), cube_rows={1: wall}, **kw
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
    left = list(ctx.left)
    terminal = [False] * N
    terminal[4] = True
    pol = _Echo()
    held = SW._pass(list(range(N)), left, pol.prepare(ctx, None), backward=False, terminal=terminal)
    assert held[4] is SILENCE, "a terminal must hold SILENCE, delivered and empty"
    assert (3, 4) not in pol.hops[False], "the kernel was asked for the hop into the terminal"
    assert held[5].level_gdna[0] == 4.0, "the terminal's own sending was blocked; only receiving is"
    assert held[0] is NO_NEIGHBOUR, "a terminal rule must not turn an open side into silence"
    loud = _Echo()
    unmasked = SW._pass(list(range(N)), left, loud.prepare(ctx, None), backward=False)
    assert unmasked[4].level_gdna[0] == 3.0 and (3, 4) in loud.hops[False]


def test_the_terminal_predicate_is_the_solve_gates_lock_on_a_region():
    """One predicate, two names must not appear: the slots the backbone never delivers into are exactly
    the REGIONS `g1_locked` locks — no admissible RNA strand — read off the source, so a reader cannot
    find a second definition of "terminal" in the file."""
    import inspect
    import re

    src = inspect.getsource(SW)
    assert re.search(r"terminal = .*is_region & g1_locked\(fp, fn\)", src), "the predicate moved"
    assert re.search(r"_pass\(.*terminal=", src), "the passes are no longer told the terminals"


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE LOCUS BLOCKS — the sweep solved a block at a time is the sweep, for every block size.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def _six(belief):
    return {
        f: np.asarray(getattr(belief, f))
        for f in ("f_g", "f_pos", "f_neg", "var_gdna", "var_pos", "var_neg", "informed")
    }


def test_the_block_solve_is_the_chain_solve_for_every_block_size(sweep_inputs):
    """The property the whole decomposition stands on, on the real toy chain through the shipped
    policy: ``solve_chain`` with the chain as one block, one locus per block, and every block size in
    between gives the same belief to the bit and the same ``informed`` predicate — and the diagnostic
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
        cap: dict = {}
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
        for key in (
            "f_g",
            "fg_loc",
            "_tau0_lam",
            "held_composition",
            "solvable",
            "mass_global",
            "count",
        ):
            assert np.array_equal(caps[bs][key], caps[None][key]), (
                f"block_slots={bs}: capture {key} differs"
            )
        a, b = caps[bs]["lam_rows"], caps[None]["lam_rows"]
        assert (a is None) == (b is None) and (a is None or np.array_equal(a, b))
        assert caps[bs]["backbone_assertions"] == caps[None]["backbone_assertions"]
        # what each node HEARD is the same; at a block's first slot — a terminal, which hears nothing
        # by the structural rule — an open side (``None``) and SILENCE are the same hearing
        for side in ("from_left", "from_right"):
            heard = [[m is not None and not m.is_silent for m in caps[k][side]] for k in (bs, None)]
            assert heard[0] == heard[1], f"block_slots={bs}: {side} differs in what was heard"
            for i in np.flatnonzero(terminal):
                assert all(
                    caps[k][side][i] is None or caps[k][side][i].is_silent for k in (bs, None)
                )


def test_a_block_view_rebases_the_links_and_slices_every_per_slot_array(sweep_inputs):
    """`_slots`: a neighbour outside the block is no neighbour; every per-slot array is the chain's
    slice; a 2-D bank keeps its columns; ``n_slots`` follows."""
    chain, statics, geometry, belief, _ra = sweep_inputs["args"]
    n = int(chain.n_slots)
    sl = slice(2, min(9, n))
    c = SW._slots(chain, sl)
    assert c.n_slots == sl.stop - sl.start and np.array_equal(c.kind, np.asarray(chain.kind)[sl])
    left = np.asarray(chain.left)[sl] - sl.start
    assert np.array_equal(c.left, np.where((left >= 0) & (left < c.n_slots), left, -1))
    assert c.left[0] == -1 and c.right[-1] == -1
    g = SW._slots(geometry, sl)
    assert g.n_slots == c.n_slots and g.unspliced_count.shape == (c.n_slots, 2)
    assert np.array_equal(g.eff_gdna, np.asarray(geometry.eff_gdna)[sl])
    b = SW._slots(belief, sl)
    assert np.array_equal(b.f_g, np.asarray(belief.f_g)[sl])
    assert SW._slots(statics, sl).boundary_flags.shape == (c.n_slots,)


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


def _memo_kw(sweep_inputs):
    kw = dict(sweep_inputs["kw"])
    kw.pop("block_slots", None)
    kw.pop("message_memo", None)
    return kw


def test_a_memo_hit_reproduces_the_uncached_sweep_to_the_bit_and_skips_the_layer(sweep_inputs):
    """Two sweeps on identical inputs through one memo: the second hits every block, never calls the
    policy's `prepare`, and returns the same belief and the same ``informed`` as the first and as a
    sweep with no memo at all. Not vacuous: the toy delivers rows, so a stale or empty hit would move
    the numbers."""
    from _transfer_harness import _full_policy

    policy = _full_policy(sweep_inputs)[0]
    kw = _memo_kw(sweep_inputs)
    plain = SW.solve_chain(*sweep_inputs["args"], **kw, policy=policy, block_slots=5)
    memo = SW.MessageMemo()
    first = SW.solve_chain(
        *sweep_inputs["args"], **kw, policy=policy, block_slots=5, message_memo=memo
    )
    assert memo.misses > 0 and memo.hits == 0
    n_blocks = memo.misses
    calls = []
    orig = type(policy).prepare

    def spy(self, ctx, library):
        calls.append(ctx.n_slots)
        return orig(self, ctx, library)

    type(policy).prepare = spy
    try:
        second = SW.solve_chain(
            *sweep_inputs["args"], **kw, policy=policy, block_slots=5, message_memo=memo
        )
    finally:
        type(policy).prepare = orig
    assert not calls, "a hit must not prepare the policy"
    assert memo.hits == n_blocks and memo.misses == n_blocks
    for out in (first, second):
        for f in ("f_g", "f_pos", "f_neg", "var_gdna", "var_pos", "var_neg", "informed"):
            assert np.array_equal(np.asarray(getattr(out, f)), np.asarray(getattr(plain, f))), f
    assert memo.nbytes > 0


def test_PERTURBATION_the_memo_misses_when_any_input_the_message_layer_reads_changes(sweep_inputs):
    """The key is a digest of EVERY input the message layer reads, so it is safe by construction: a
    changed belief, liveness bit, factory row, observation, library or grid must miss — and a memo that
    hit on any of them would deliver another sweep's messages as this one's."""
    import dataclasses as _dc

    from _transfer_harness import _full_policy

    policy, rows, _g, _w = _full_policy(sweep_inputs)
    kw = dict(
        _memo_kw(sweep_inputs), intron_prior=rows
    )  # the live rows: the layer has something to read
    chain, statics, geometry, belief, ra = sweep_inputs["args"]
    memo = SW.MessageMemo()
    SW.solve_chain(chain, statics, geometry, belief, ra, **kw, policy=policy, message_memo=memo)
    base_misses = memo.misses

    def misses_after(pol=policy, **over):
        before = memo.misses
        a = over.pop("args", (chain, statics, geometry, belief, ra))
        SW.solve_chain(*a, **{**kw, **over}, policy=pol, message_memo=memo)
        return memo.misses - before

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
    """An instrument's capture reads the held lists, so with ``_capture`` the layer runs even on a memo
    that would hit — and what it delivers equals the memo's, so the two paths cannot drift."""
    from _transfer_harness import _full_policy

    policy = _full_policy(sweep_inputs)[0]
    kw = _memo_kw(sweep_inputs)
    memo = SW.MessageMemo()
    SW.solve_chain(*sweep_inputs["args"], **kw, policy=policy, message_memo=memo)
    cap: dict = {}
    hits_before = memo.hits
    out = SW.solve_chain(
        *sweep_inputs["args"], **kw, policy=policy, message_memo=memo, _capture=cap
    )
    assert memo.hits == hits_before, "a captured sweep must not be served from the memo"
    assert "from_left" in cap and out.informed is not None
