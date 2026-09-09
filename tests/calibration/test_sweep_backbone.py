"""The BACKBONE's four assertions, and the contract that keeps the shipped policy the shipped policy.

⛔⛔ **TRAPS: perturb-every-gate IS THE WHOLE SHAPE OF THIS FILE.** Writing a gate before the fix is half the discipline;
the other half is breaking the fixed code and watching each gate fire. So every assertion here has a
matching PERTURBATION test that constructs a policy committing exactly that defect and asserts the backbone
refuses it. A gate with no firing perturbation has not been written yet — it has been typed.

⭐ The per-condition byte-identity of the restructure against the shipped solver is NOT gated here, because
it needs a real 70,176-slot chain and a BAM. It is
``scripts/design/backbone_parity.py`` (421,056 output elements and 18,245,830 diagnostic elements, zero
differences when the restructure landed).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration import sweep as SW
from rigel.calibration.messages import NO_NEIGHBOUR, SILENCE, Message, PsiMessage, StepContext
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.simplex_logodds import _logodds_grid


N = 8


def _ctx(*, free_pos=None, free_neg=None, n_grid=60) -> StepContext:
    """A minimal StepContext. Only the fields the assertions read need to be real."""
    ones = np.ones(N)
    fp = np.ones(N, bool) if free_pos is None else np.asarray(free_pos, bool)
    fn = np.zeros(N, bool) if free_neg is None else np.asarray(free_neg, bool)
    lam, grid = _logodds_grid(n_grid, 10.0)
    return StepContext(
        mass=ones * 100.0,
        inv_abundance=ones * 0.5,
        inv_sj_lo=np.zeros((N, 2)),
        inv_sj_hi=np.zeros((N, 2)),
        eff_gdna_global=ones * 200.0,
        eff_rna=ones * 200.0,
        eff_gdna=ones * 200.0,
        eff_sj=np.ones((N, 2)) * 200.0,
        sj_count_lo=np.zeros((N, 2)),
        sj_count_hi=np.zeros((N, 2)),
        sj_count=np.zeros((N, 2)),
        route_rate_lo=np.zeros((N, 2)),
        route_rate_hi=np.zeros((N, 2)),
        route_count_lo=np.zeros((N, 2), dtype=np.int64),
        route_count_hi=np.zeros((N, 2), dtype=np.int64),
        unspliced_count=np.ones((N, 2)) * 5.0,
        n_slot=ones * 10.0,
        spliced_slot=np.zeros(N),
        left=np.arange(-1, N - 1),
        right=np.append(np.arange(1, N), -1),
        is_boundary=np.arange(N) % 2 == 1,
        is_exon_region=np.zeros(N, bool),
        left_interface_certified=np.zeros(N, bool),
        right_interface_certified=np.zeros(N, bool),
        ss_intron_boundary=np.zeros(N, bool),
        free_pos=fp,
        free_neg=fn,
        exon_pos=np.zeros(N, bool),
        exon_neg=np.zeros(N, bool),
        boundary_flags=np.zeros(N, np.int64),
        geometry=None,
        order=list(range(N)),
        left_list=list(range(-1, N - 1)),
        right_list=[*range(1, N), -1],
        own=None,
        belief_fg=ones,
        n_grid=n_grid,
        logodds_window=10.0,
        solve_grid=grid,
    )


def _counts(msg: PsiMessage, ctx: StepContext | None = None):
    c = SW.AssertionCounts()
    SW._check_message(msg, ctx if ctx is not None else _ctx(), c)
    return c


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 1 — THE TWO PHASES (owner ruling 2026-09-04): every node holds a message from each neighbour
# it has; a real hop must ARRIVE; a missing neighbour is not silence; the kernel sees indices only.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


class _Echo:
    """A policy whose kernel records every hop and returns a message naming its source — so the pass's
    ORDER and SIDES are observable, and so the solve can be shown the two held lists."""

    name = "echo"

    def __init__(self):
        self.hops = {False: [], True: []}
        self.held = None

    def prepare(self, ctx):
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
    """⭐⭐ The ruling made executable: after the two passes every interior node holds TWO messages, one
    naming its low neighbour and one its high; the chain's two end nodes hold ONE and ``NO_NEIGHBOUR``
    on the open side — which is not a message and not SILENCE."""
    ctx = _ctx()
    left, right = list(ctx.left), list(ctx.right)
    fl = SW._pass(list(ctx.order), left, _Echo().prepare(ctx), backward=False)
    br = SW._pass(list(ctx.order)[::-1], right, _Echo().prepare(ctx), backward=True)
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
    SW._pass(list(ctx.order), list(ctx.left), pol.prepare(ctx), backward=False)
    SW._pass(list(ctx.order)[::-1], list(ctx.right), pol.prepare(ctx), backward=True)
    assert pol.hops[False] == [(i - 1, i) for i in range(1, N)]
    assert pol.hops[True] == [(i + 1, i) for i in range(N - 2, -1, -1)]


def test_PERTURBATION_a_kernel_that_leaves_a_real_hop_unspoken_is_REFUSED():
    """⛔ A hop that carries nothing must still ARRIVE as SILENCE. A kernel returning ``None`` for a node
    that HAS a neighbour is the one thing the pass refuses, because the solve could not then tell
    "nothing to say" from "never spoken to"."""

    class _Mute(_Echo):
        def propagate(self, *, backward: bool):
            return lambda s, i: None

    ctx = _ctx()
    with pytest.raises(AssertionError, match="ARRIVE as SILENCE"):
        SW._pass(list(ctx.order), list(ctx.left), _Mute().prepare(ctx), backward=False)


def test_a_policy_that_sends_nothing_leaves_silence_at_every_node_with_a_neighbour():
    """``propagate`` returning no kernel means every node holds SILENCE from that side — delivered,
    distinguishable from the open side of the chain."""
    ctx = _ctx()

    class _Quiet(_Echo):
        def propagate(self, *, backward: bool):
            return None

    fl = SW._pass(list(ctx.order), list(ctx.left), _Quiet().prepare(ctx), backward=False)
    assert fl[0] is NO_NEIGHBOUR and all(m is SILENCE for m in fl[1:])
    assert SILENCE.is_silent and Message(level_gdna=(0.0, 1.0)).is_silent is False


def test_every_lane_of_a_message_survives_the_passes_to_the_solve():
    """THE LANES (owner ruling 2026-09-04): a kernel that fills every lane — both composition profiles
    and the three level claims — hands them to the solve untouched, and a message with any one lane is
    not silent. The backbone carries; it never reads a lane."""
    ctx = _ctx()
    full = Message(
        composition=np.zeros(3),
        tilt=np.zeros(2),
        level_gdna=(-2.0, 0.1),
        level_rna_pos=(-3.0, 0.2),
        level_rna_neg=(-4.0, 0.3),
    )

    class _Full(_Echo):
        def propagate(self, *, backward: bool):
            return lambda s, i: full

    pol = _Full()
    prepared = pol.prepare(ctx)
    fl = SW._pass(list(ctx.order), list(ctx.left), prepared, backward=False)
    br = SW._pass(list(ctx.order)[::-1], list(ctx.right), prepared, backward=True)
    prepared.solve(fl, br)
    got_l, got_r = pol.held
    for i in range(N):
        if list(ctx.left)[i] >= 0:
            assert got_l[i] is full
        if list(ctx.right)[i] >= 0:
            assert got_r[i] is full
    assert not full.is_silent and Message().is_silent
    for lane in Message.LANES:
        one = Message(**{lane: np.zeros(2) if lane in ("composition", "tilt") else (0.0, 1.0)})
        assert not one.is_silent, f"a message with only {lane} read as silence"


def test_the_solve_receives_the_two_held_lists_at_the_recipient():
    """Phase 2's inputs are the two held lists indexed AT THE RECIPIENT, straight from the passes."""
    ctx = _ctx()
    pol = _Echo()
    prepared = pol.prepare(ctx)
    fl = SW._pass(list(ctx.order), list(ctx.left), prepared, backward=False)
    br = SW._pass(list(ctx.order)[::-1], list(ctx.right), prepared, backward=True)
    prepared.solve(fl, br)
    assert pol.held == (fl, br)


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 2 — every delivered ROW is one row per slot on the solve grid, and finite.
# (The coordinate and share gates that stood here guarded the retired relay's Gaussian channels —
# TRAPS: off-grid-message-mode — and retired with them on 2026-09-09; a profile on ψ's own grid cannot
# be delivered off-grid nor claim an over-unit share.)
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_lambda_rows_are_checked_for_shape_and_finiteness():
    """The λ-row channel: one row per slot on the solve grid is accepted and counted finite; a row array
    of another shape is REFUSED outright; a non-finite row is COUNTED and, with no waiver, raises."""
    ctx = _ctx()
    K = int(ctx.n_grid)
    ok = _counts(PsiMessage(lam_rows=np.zeros((N, K))), ctx)
    assert ok["flux_rows_finite"] == {"violations": 0, "eligible": N}
    with pytest.raises(ValueError, match="lam_rows has shape"):
        _counts(PsiMessage(lam_rows=np.zeros((N + 1, K))), ctx)
    with pytest.raises(ValueError, match="lam_rows has shape"):
        _counts(PsiMessage(lam_rows=np.zeros(N)), ctx)
    bad = np.zeros((N, K))
    bad[2, 0] = np.nan
    with pytest.raises(AssertionError, match="flux_rows_finite"):
        _counts(PsiMessage(lam_rows=bad), ctx)


def test_a_waiver_is_never_silent():
    """Every waived assertion carries a written reason, so a reader learns the defect rather than the
    exemption. ⛔ An empty reason would be a widened predicate wearing a waiver's clothes."""
    for name, why in SW._KNOWN_VIOLATIONS.items():
        assert len(why) > 80, f"{name}'s waiver does not say what the defect is"


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 3 — |T| <= 3.  AXIOM 0, made executable.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_population_set_is_at_most_three_because_it_is_a_function_of_two_bits():
    """⭐⭐ ``T(slot) = {gDNA} u {RNA+ if free_pos} u {RNA- if free_neg}``, so ``|T| = 1 + free_pos +
    free_neg`` and it is in ``{1, 2, 3}`` for every possible input. **There are THREE populations and there
    is no fourth** — "mature" and "nascent" are not species, and RNA inside an intron is RNA that has not
    spliced at that position. This is structural rather than something to remember, and that is the point."""
    for fp in (True, False):
        for fn in (True, False):
            ctx = _ctx(free_pos=np.full(N, fp), free_neg=np.full(N, fn))
            pop = ctx.population_size()
            assert set(np.unique(pop)) <= {1, 2, 3}
            assert np.all(pop == 1 + int(fp) + int(fn))


def test_PERTURBATION_a_fourth_population_is_REFUSED():
    """⛔ AXIOM 0's tell, executable: a population set with more than three members. A derivation once
    opened with ``{gDNA, nascent+, nascent-, mature+, mature-}`` and every table built on it came out
    wrong."""
    counts = SW.AssertionCounts()
    with pytest.raises(AssertionError, match="population_at_most_three"):
        counts.note("population_at_most_three", np.array([4, 4, 5]) > 3, np.ones(3, bool))


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 4 — the write-back touches only `solvable` slots.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_PERTURBATION_a_writeback_outside_solvable_is_REFUSED():
    """⛔ The silent version of this made an TRAPS: byte-identity-gate identity gate read ``max|delta| = 1.0``: a replay compared
    the solve's raw output against the shipped belief, and the two differ by exactly this mask. Reproducing
    a pipeline stage means reproducing its WRITE-BACK.

    ⭐ A locked slot — one with no admissible RNA strand — is never solved and keeps its signature-binary
    init, because RNA cannot cross a gene boundary so its unspliced mass is purely gDNA."""
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
    """⛔⛔ **THE LARGEST BEHAVIOUR SWITCH IN THE TOOL, AND IT MUST BE A WRITTEN DECISION.** Which policy
    ships can never be inherited from a function default that an edit could silently change: the config
    names it (``message_policy``, ``"transfer"`` since 2026-09-09), ``calibrate`` reads it, and both the
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
    """⭐ ``SilentPolicy`` is ``solve_chain``'s PARAMETER default (the shipped config installs the
    transfer policy), and it is the MEASURED floor every policy is judged against: win on unstranded
    data, minimal harm on stranded data, never pooled."""
    prepared = SilentPolicy().prepare(_ctx())
    assert prepared.propagate(backward=False) is None, "a silent policy must send nothing at all"
    assert prepared.propagate(backward=True) is None
    assert prepared.solve([SILENCE] * N, [SILENCE] * N).is_silent


def test_the_backbone_does_not_know_what_a_message_is_about():
    """⭐⭐⭐ **THE STRUCTURAL CLAIM OF THE WHOLE RESTRUCTURE, as a test.** The backbone owns the shape of the
    solve and the five assertions; every message-composition choice is a policy. If one of these concepts
    reappears in ``sweep.py``, an operator has leaked back into the backbone and the next reader can no
    longer hold the working system in their head.

    ⚠ **It checks IDENTIFIERS, from the AST — not the file's text.** Grepping the source would match the
    module docstring, which names these very words in order to say they are absent; a test that passes for
    that reason is vacuous, and this one failed exactly that way when first written.

    ⚠ **``capture`` has ONE licensed occurrence and it is not the biology**: ``_capture`` is the diagnostics
    hook every instrument passes by keyword, and it carried that name in the shipped solver. Hybrid capture —
    the thing the message layer argues about — appears nowhere."""
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
    )
    leaked = {w: sorted(i for i in ident if w in i) for w in banned}
    assert not any(leaked.values()), (
        f"policy concepts leaked into the backbone's identifiers: {leaked}"
    )
    cap = sorted(i for i in ident if "capture" in i)
    assert cap == ["_capture", "capture"], (
        f"the only licensed 'capture' is the diagnostics hook — the parameter and the StepContext field "
        f"that carries it. Found {cap}"
    )


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE CUBE CHANNEL (the both-stranded locus, 2026-09-08): a (K, K_t) row per AMBIG slot, final solve only
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
