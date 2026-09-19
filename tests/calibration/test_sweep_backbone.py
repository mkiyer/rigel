"""The backbone's four assertions, and the contract that keeps the shipped policy the shipped policy.

The four: every node ends the two passes holding one message from each neighbour it has; every
delivered row is one row per slot on the solve grid and finite; a slot's population set has at most
three members; and the write-back touches only solvable slots. The kernel counts them per block
(`native/solve_kernel.cpp`) and the backbone judges (`sweep.AssertionCounts`). TRAPS: perturb-every-gate
is the shape of the whole file — each assertion has a matching perturbation test that hands the backbone
exactly that defect and asserts it refuses, because a gate with no firing perturbation has not been
written yet, it has been typed. The pass's order and sides, the two states a node can hold from a side
and the terminal rule are read off the kernel's received tables on hand-built chains through the gates'
harness (`_transfer_harness`). Byte-identity against the shipped solver per condition is not gated here:
it needs a real chain and a BAM, and it is ``scripts/design/rename_identity.py`` with
``scripts/profiling/sweep_replay.py``.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration import sweep as SW
from rigel.calibration.blocks import SweepCapture
from rigel.calibration.message_cache import MessageCache
from rigel.calibration.messages import ChainView
from _transfer_harness import (
    FORWARD,
    Faces,
    LevelLane,
    Prepared,
    Received,
    RowTable,
    _lane_prepared,
    _native_passes,
    norm,
)


N = 8
K = 60


def _ctx(*, free_pos=None, free_neg=None, n_grid=K) -> ChainView:
    """A minimal chain view: a chain of N slots, ``N E N E …``, every side linked. Only the fields the
    assertions and the passes read need to be real."""
    ones = np.ones(N)
    fp = np.ones(N, bool) if free_pos is None else np.asarray(free_pos, bool)
    fn = np.zeros(N, bool) if free_neg is None else np.asarray(free_neg, bool)
    return ChainView(
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
        boundary_flags=np.zeros(N, np.uint16),
        n_grid=n_grid,
        logodds_window=10.0,
    )


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 1 — THE TWO PHASES: every node holds a message from each neighbour it has; a real hop must
# ARRIVE; a missing neighbour is not silence; the passes run in chain order and read one side each.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def _echo_lane(u):
    """A gDNA lane on a chain of EMPTY nodes, each holding its own lower-sided level with a flux witness
    NAMING ITSELF. An empty node forwards what it holds and emits the intersection of its own level and the
    held one with ITS witness, and an empty recipient does not re-price — so after a pass the count a node
    holds names the node it received from, and the pass's order and sides are observable."""
    own = RowTable(N, u.shape[0])
    witness = RowTable(N, 2)
    for i in range(N):
        own[i] = np.where(u < -5.0 + 0.5 * i, -50.0, 0.0)
        witness[i] = np.array([float(i), 1.0])
    return LevelLane(
        "gdna",
        u,
        u,
        0.5,
        np.zeros(N),
        np.ones(N),
        np.ones(N, bool),
        own,
        np.ones((N, 2), bool),
        flux_witness=witness,
    )


def _echo(ctx, terminal=None):
    """The two passes of the echo lane on the chain view's links; ``terminal`` marks the nodes that receive
    nothing. Returns ``(from_left, from_right)``."""
    u = np.linspace(-10.0, 10.0, int(ctx.n_grid))
    prepared = _lane_prepared(_echo_lane(u), ctx.left, ctx.right)
    order = np.arange(N, dtype=np.int64)
    term = np.zeros(N, bool) if terminal is None else np.asarray(terminal, bool)
    tables = []
    for nbr, seq, backward in ((ctx.left, order, False), (ctx.right, order[::-1], True)):
        nbr = np.asarray(nbr, np.int64)
        received = Received.empty(N, int(ctx.n_grid))
        received.has_neighbour[seq] = nbr[seq] >= 0
        prepared.run_pass(received, seq, nbr, term, backward=backward)
        tables.append(received)
    return tuple(tables)


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
    fl, br = _echo(ctx)
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
    """The forward pass visits low→high reading each node's LOW neighbour, the backward pass the mirror,
    so what a source holds from its far side is written before it is asked to send: on a chain of
    FORWARD faces with a distinct own claim at every node, the composition a node holds from the left is
    its low neighbour's claim fused with what THAT neighbour held — the recursion down the chain — and
    from the right the mirror. PERTURBATION: a pass out of order, or one reading the wrong side, composes
    the wrong claims and fails the recursion at the first interior node."""
    ctx = _ctx()
    lam = np.linspace(-10.0, 10.0, K)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    own = RowTable(N, K)
    for i in range(N):
        own[i] = norm(-0.5 * ((lam - (-6.0 + 1.6 * i)) / 0.4) ** 2)
    faces = Faces(lam, left, right)
    faces.kind[left >= 0, 0] = FORWARD
    faces.kind[right >= 0, 1] = FORWARD
    prepared = Prepared.hand_built(own, faces, {}, np.ones(N, bool), np.zeros(N, bool))
    fl, br = _native_passes(prepared, ctx)
    want_l, want_r = [None] * N, [None] * N
    for i in range(1, N):
        s = i - 1
        want_l[i] = norm(own[s] if want_l[s] is None else own[s] + want_l[s])
    for i in range(N - 2, -1, -1):
        s = i + 1
        want_r[i] = norm(own[s] if want_r[s] is None else own[s] + want_r[s])
    for i in range(N):
        assert fl.has_composition[i] == (want_l[i] is not None)
        assert br.has_composition[i] == (want_r[i] is not None)
        if want_l[i] is not None:
            np.testing.assert_allclose(fl.composition[i], want_l[i], atol=1e-12)
        if want_r[i] is not None:
            np.testing.assert_allclose(br.composition[i], want_r[i], atol=1e-12)
    assert not fl.has_level.any() and not br.has_level.any(), "a FORWARD face carries no level"


def test_a_policy_that_sends_nothing_leaves_silence_at_every_node_with_a_neighbour():
    """A pass through tables that carry no rule and no lane writes nothing, so every node holds SILENCE
    from that side — a neighbour and nothing present — distinguishable from the open side of the chain."""
    ctx = _ctx()
    lam = np.linspace(-10.0, 10.0, K)
    quiet = Prepared.hand_built(
        RowTable(N, K), Faces(lam, ctx.left, ctx.right), {}, np.ones(N, bool), np.zeros(N, bool)
    )
    for t in _native_passes(quiet, ctx):
        assert np.array_equal(t.silence, t.has_neighbour) and t.silence.sum() == N - 1
        assert np.array_equal(t.no_neighbour, ~t.has_neighbour) and t.no_neighbour.sum() == 1
        assert not t.heard.any()


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 2 — every delivered ROW is finite (one row per slot on the solve grid is structural now: the
# kernel's rows are the grid's), counted per block by the kernel and judged by the backbone.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════

_NAMES = (
    "population_at_most_three",
    "population_reaches_three",
    "lam_rows_finite",
    "cube_rows_finite",
    "writeback_only_solvable",
)


def _counts(per_block, rows_delivered, cube_delivered, n_owned):
    """`AssertionCounts.of_blocks` on hand-written per-block counts ``{name: [(violations, eligible), …]}``."""
    B = len(n_owned)
    arr = np.zeros((B, len(_NAMES), 2), np.int64)
    for a, name in enumerate(_NAMES):
        for b, pair in enumerate(per_block.get(name, [(0, 0)] * B)):
            arr[b, a] = pair
    delivered = {
        "rows_delivered": np.asarray(rows_delivered, bool),
        "cube_delivered": np.asarray(cube_delivered, bool),
    }
    return SW.AssertionCounts.of_blocks(list(_NAMES), arr, delivered, np.asarray(n_owned))


def test_a_non_finite_delivered_row_is_counted_and_refused():
    """The λ-row channel: finite rows are counted eligible; a non-finite row is COUNTED and, with no
    waiver, raises — and the check is absent from the report when no block delivered rows at all."""
    ok = _counts({"lam_rows_finite": [(0, 5), (0, 3)]}, [True, True], [False, False], [5, 3])
    assert ok["lam_rows_finite"] == {"violations": 0, "eligible": 8}
    assert "cube_rows_finite" not in ok
    with pytest.raises(AssertionError, match="lam_rows_finite"):
        _counts({"lam_rows_finite": [(0, 5), (1, 3)]}, [True, True], [False, False], [5, 3])
    quiet = _counts({}, [False, False], [False, False], [5, 3])
    assert "lam_rows_finite" not in quiet and "cube_rows_finite" not in quiet


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
        counts.note("population_at_most_three", 3, 3)


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
    counts = SW.AssertionCounts()
    with pytest.raises(AssertionError, match="writeback_only_solvable"):
        counts.note("writeback_only_solvable", 1, 2)


def test_a_writeback_confined_to_solvable_is_accepted():
    counts = SW.AssertionCounts()
    counts.note("writeback_only_solvable", 0, 2)
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


def test_solve_chains_parameter_default_is_silent_and_sends_nothing(sweep_inputs):
    """``SilentPolicy`` is ``solve_chain``'s parameter default (the shipped config installs the
    transfer policy), and it is the MEASURED floor every policy is judged against: win on unstranded
    data, minimal harm on stranded data, never pooled. The kernel runs no layer for it: the capture of a
    default sweep names the silent policy and holds no received table and no delivered row."""
    import inspect

    assert inspect.signature(SW.solve_chain).parameters["policy"].default is None
    assert SW._POLICY_KERNEL["silent"] == 0
    cap = SweepCapture()
    out = SW.solve_chain(*sweep_inputs["args"], **sweep_inputs["kw"], _capture=cap)
    assert cap.policy_name == "silent" and cap.from_left is None and cap.lam_rows is None
    assert out.has_composition is not None and "lam_rows_finite" not in cap.backbone_assertions


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
# THE CUBE CHANNEL (the both-stranded locus): a row per delivered AMBIG slot, final solve only
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_cube_channel_is_counted_where_a_cube_was_delivered():
    """The cube rows' finiteness is counted over the delivered rows at a block's OWNED slots — absent
    from the report when no block delivered a cube, a non-finite profile refused — and the kernel's
    solve writes a cube row only at an AMBIG node (`transfer_kernel.h`), so a cube at a single-strand
    slot cannot be built."""
    ok = _counts({"cube_rows_finite": [(0, 2), (0, 0)]}, [True, True], [True, False], [5, 3])
    assert ok["cube_rows_finite"] == {"violations": 0, "eligible": 2}
    with pytest.raises(AssertionError, match="cube_rows_finite"):
        _counts({"cube_rows_finite": [(1, 2), (0, 0)]}, [True, True], [True, False], [5, 3])


def test_the_solvers_cube_is_inert_when_absent_and_walls_the_tilt_when_present():
    """`_solve_regions_logodds_all(cube_rows=...)`: ``None`` and an empty table are byte-identical to the
    path without the argument; a wall against low ``f_+`` at one AMBIG slot raises that slot's ``f_pos`` and
    leaves every other slot byte-identical."""
    from _psi_reference import cube_rows_of

    from rigel.calibration.simplex_logodds import CubeRows, _solve_regions_logodds_all

    m, K_ = 4, 40
    u_pos = np.array([50.0, 50.0, 50.0, 50.0])
    u_neg = np.array([50.0, 50.0, 50.0, 50.0])
    ap = np.ones(m, bool)
    an = np.array([True, True, False, True])
    kw = dict(kappa=0.99, od_g=0.0, od_r=0.0, n_grid=K_, L=10.0)
    base = _solve_regions_logodds_all(u_pos, u_neg, ap, an, u_pos + u_neg, np.zeros(m), **kw)
    u = np.linspace(-10.0, 10.0, K_)
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
    side it has a neighbour on, and the kernel is never asked for the hop into it — so what a policy
    WOULD deliver there cannot exist. PERTURBATION: the same pass with no terminal marked delivers its
    level, which is what proves the mask does the work. What the terminal SENDS is untouched: its
    neighbour still receives from it."""
    ctx = _ctx()
    terminal = np.zeros(N, bool)
    terminal[4] = True
    held, _br = _echo(ctx, terminal=terminal)
    assert held.has_neighbour[4] and held.silence[4], (
        "a terminal holds SILENCE: a neighbour, nothing present"
    )
    assert held.level_gdna.count[5] == 4.0, (
        "the terminal's own sending was blocked; only receiving is"
    )
    assert held.no_neighbour[0], "a terminal rule must not turn an open side into silence"
    unmasked, _br = _echo(ctx)
    assert unmasked.level_gdna.present[4] and unmasked.level_gdna.count[4] == 3.0


def test_the_terminal_predicate_is_the_solve_gates_lock_on_a_region():
    """One predicate, two names must not appear: the slots the backbone never delivers into are exactly
    the REGIONS `g1_locked` locks — no admissible RNA strand — read off the source, so a reader cannot
    find a second definition of "terminal" in the file; and the kernel is told them."""
    import inspect
    import re

    src = inspect.getsource(SW)
    assert re.search(r"terminal ?= ?.*is_region & g1_locked\(fp, fn\)", src), "the predicate moved"
    assert re.search(r"terminal=bits\(terminal\)", src), (
        "the kernel is no longer told the terminals"
    )


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE LOCUS BLOCKS — the sweep solved a block at a time is the sweep, for every block size and every
# thread count.
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
        for key in ("f_g", "fg_loc", "tau_lam", "tau_fac", "solvable", "mass_global", "count"):
            assert np.array_equal(getattr(caps[bs], key), getattr(caps[None], key)), (
                f"block_slots={bs}: capture {key} differs"
            )
        a, b = caps[bs].lam_rows, caps[None].lam_rows
        assert (a is None) == (b is None) and (a is None or np.array_equal(a, b))
        assert caps[bs].backbone_assertions == caps[None].backbone_assertions
        # what each node HEARD is the same; at a block's first slot — a terminal, which hears nothing
        # by the structural rule — an open side and SILENCE are the same hearing
        for side in ("from_left", "from_right"):
            heard = [Received.from_kernel(getattr(caps[k], side)).heard for k in (bs, None)]
            assert np.array_equal(heard[0], heard[1]), (
                f"block_slots={bs}: {side} differs in what was heard"
            )
            for i in np.flatnonzero(terminal):
                assert not any(h[i] for h in heard)


def test_the_block_solve_is_thread_exact_so_the_thread_count_moves_no_number(sweep_inputs):
    """The blocks are pulled one at a time by a pool of threads, each block solved by the same arithmetic
    on its own arena and written to its own slots, so the sweep is BIT-IDENTICAL at every thread count —
    the budget is a resource, not a tunable of the answer. On the toy cut into many blocks, at 1, 2, 3 and
    every core."""
    from _transfer_harness import _full_policy

    policy = _full_policy(sweep_inputs)[0]
    kw = dict(sweep_inputs["kw"])
    kw.pop("block_slots", None)
    kw.pop("n_threads", None)
    serial = _six(
        SW.solve_chain(*sweep_inputs["args"], **kw, policy=policy, block_slots=2, n_threads=1)
    )
    for n_threads in (2, 3, 0):
        got = _six(
            SW.solve_chain(
                *sweep_inputs["args"], **kw, policy=policy, block_slots=2, n_threads=n_threads
            )
        )
        for f in serial:
            assert np.array_equal(got[f], serial[f]), f"{n_threads} threads: {f} moved"


def test_the_checks_count_only_the_owned_slots():
    """The λ-row check's ELIGIBLE set reads as the chain's: where any block delivered rows, a block that
    delivered none holds zero rows — finite rows that were checked — so the published count does not
    depend on how the chain was cut; and the kernel counts each block's OWNED slots only (its read-ahead
    terminal is another block's)."""
    two = _counts(
        {"lam_rows_finite": [(0, 5), (0, 0)], "population_at_most_three": [(0, 5), (0, 3)]},
        [True, False],
        [False, False],
        [5, 3],
    )
    assert two["lam_rows_finite"] == {"violations": 0, "eligible": 8}
    assert two["population_at_most_three"] == {"violations": 0, "eligible": 8}
    merged = SW.AssertionCounts()
    merged.note("lam_rows_finite", 0, 7)
    merged.note("lam_rows_finite", 0, 7)
    assert merged["lam_rows_finite"] == {"violations": 0, "eligible": 14}


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE MESSAGE MEMO — the message layer is refit-invariant given the grid, so the refit sweeps share it.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def _cache_kw(sweep_inputs):
    kw = dict(sweep_inputs["kw"])
    kw.pop("block_slots", None)
    kw.pop("message_cache", None)
    return kw


def test_a_cache_hit_reproduces_the_uncached_sweep_to_the_bit_and_skips_the_layer(sweep_inputs):
    """Two sweeps on identical inputs through one cache: the second is served every block — the kernel
    runs no layer for a served block — and returns the same belief and the same ``has_composition`` as
    the first and as a sweep with no cache at all. Not vacuous: the toy delivers rows, so a stale or empty
    hit would move the numbers."""
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
    second = SW.solve_chain(
        *sweep_inputs["args"], **kw, policy=policy, block_slots=5, message_cache=cache
    )
    assert cache.hits == n_blocks and cache.misses == n_blocks
    for out in (first, second):
        for f in ("f_g", "f_pos", "f_neg", "var_gdna", "has_composition"):
            assert np.array_equal(np.asarray(getattr(out, f)), np.asarray(getattr(plain, f))), f
    assert cache.nbytes > 0


def test_PERTURBATION_the_cache_misses_when_any_input_the_message_layer_reads_changes(sweep_inputs):
    """The key is a digest of EVERY input the message layer reads, so it is safe by construction: a
    changed belief, factory row, observation, library or grid must miss — and a cache that hit on any of
    them would deliver another sweep's messages as this one's."""
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
    """An instrument's capture reads the held tables, so with ``_capture`` the layer runs even on a cache
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


def test_the_factory_rows_enter_the_key_by_the_digest_of_their_inputs():
    """`calibrate.FactoryRows.digest(block)` is a digest of what a block's rows are a pure function of — the
    background's fields, the block's intron mask, counts and opportunities, the grid — never of the rows'
    bytes: two factories on identical inputs digest alike per block; a changed count at one intron changes
    ITS block's digest and no other's; a changed background changes every block's; and a factory whose
    rows are only ever read through the digest hashes 1/K of the bytes. PERTURBATION: a digest that
    skipped the counts fails here."""
    from rigel.calibration.calibrate import FactoryRows
    from rigel.calibration.density_deconv import GdnaBackground

    def factory(count, background):
        f = object.__new__(FactoryRows)
        f.background = background
        f.is_intron = np.array([True, False, True, True, False, True])
        f.count = np.asarray(count, np.float64)
        f.eff = np.array([100.0, 0.0, 250.0, 80.0, 0.0, 300.0])
        f.fg = np.linspace(0.01, 0.99, 7)
        f.shape = (6, 7)
        return f

    bg = GdnaBackground(log_mu_bg=-3.0, alpha=12.0, size=40.5, n_regions=9, informative=True)
    a = factory([5.0, 0.0, 12.0, 3.0, 0.0, 8.0], bg)
    b = factory([5.0, 0.0, 12.0, 3.0, 0.0, 8.0], bg)
    blocks = (SimpleNamespace(start=0, stop=3, end=3), SimpleNamespace(start=3, stop=6, end=6))
    for bl in blocks:
        assert a.digest(bl) == b.digest(bl)
    c = factory([5.0, 0.0, 12.0, 3.0, 0.0, 9.0], bg)  # one intron's count, in the second block
    assert c.digest(blocks[0]) == a.digest(blocks[0]) and c.digest(blocks[1]) != a.digest(blocks[1])
    d = factory([5.0, 0.0, 12.0, 3.0, 0.0, 8.0], GdnaBackground(-3.1, 12.0, 40.5, 9, True))
    assert all(d.digest(bl) != a.digest(bl) for bl in blocks)
    inputs = a.kernel()
    assert inputs[0] == "inputs" and inputs[1].dtype == bool and inputs[-1] is True


def test_the_cache_keys_a_blocks_rows_by_the_factorys_digest(sweep_inputs):
    """The wiring: `solve_chain` hands the kernel the factory's inputs AND asks it for each block's digest,
    and the cache keys on the digest — a factory answering the same digest hits although its rows are
    rebuilt, and one whose digest moves misses (rows given as one array keep digesting by content: the
    perturbation gate above)."""
    from _transfer_harness import _full_policy

    policy, rows, _g, _w = _full_policy(sweep_inputs)
    chain, statics, geometry, belief, ra = sweep_inputs["args"]
    kw = _cache_kw(sweep_inputs)
    kw.pop("intron_prior", None)  # the factory below stands in for the captured rows

    class _Factory:
        def __init__(self, rows, salt):
            self.rows, self.salt, self.asked = np.asarray(rows, np.float64), salt, 0

        def kernel(self):
            return ("rows", self.rows.copy())  # rebuilt: a fresh array each time, the same numbers

        def digest(self, block):
            self.asked += 1
            return f"{self.salt}:{block.start}:{block.stop}".encode()

    cache = MessageCache()
    first = _Factory(rows, "a")
    SW.solve_chain(
        chain,
        statics,
        geometry,
        belief,
        ra,
        **kw,
        intron_prior=first,
        policy=policy,
        message_cache=cache,
    )
    misses = cache.misses
    assert misses > 0 and first.asked == misses
    same = _Factory(rows, "a")
    SW.solve_chain(
        chain,
        statics,
        geometry,
        belief,
        ra,
        **kw,
        intron_prior=same,
        policy=policy,
        message_cache=cache,
    )
    assert cache.misses == misses and cache.hits == misses, "the same digest must hit"
    other = _Factory(rows, "b")
    SW.solve_chain(
        chain,
        statics,
        geometry,
        belief,
        ra,
        **kw,
        intron_prior=other,
        policy=policy,
        message_cache=cache,
    )
    assert cache.misses == 2 * misses, "another digest must miss"
