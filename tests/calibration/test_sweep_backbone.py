"""The BACKBONE's five assertions, and the contract that keeps the shipped policy the shipped policy.

⛔⛔ **TRAPS: perturb-every-gate IS THE WHOLE SHAPE OF THIS FILE.** Writing a gate before the fix is half the discipline;
the other half is breaking the fixed code and watching each gate fire. So every assertion here has a
matching PERTURBATION test that constructs a policy committing exactly that defect and asserts the backbone
refuses it. A gate with no firing perturbation has not been written yet — it has been typed.

⭐ The per-condition byte-identity of the restructure against the shipped solver is NOT gated here, because
it needs a real 70,176-slot chain and a BAM. It is
``scripts/design/backbone_parity.py`` (421,056 output elements and 18,245,830 diagnostic elements, zero
differences) and ``ladder_arm_ab.py --arm backbone_relay`` / ``--arm backbone`` on the 36-condition panel.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration import sweep as SW
from rigel.calibration.messages import NeighbourState, PsiMessage, StepContext
from rigel.calibration.messages.relay import RelayPolicy, RelaySwitches
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.simplex_logodds import _logodds_grid, _tilt_grid


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
# ASSERTION 1 — deliver() sees only the two NEIGHBOUR states.  TRAPS: a-message-from-the-destinations-belief, nine recurrences in nine costumes.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_neighbour_state_is_gathered_at_the_source_so_d4_is_not_expressible():
    """⭐⭐ The structural half of TRAPS: a-message-from-the-destinations-belief: what ``deliver`` is handed is already indexed at the source.

    A message computed from the destination's own relayed belief carries zero information and confirms the
    destination — and every one of TRAPS: a-message-from-the-destinations-belief's nine costumes was that. After the gather there is no array in a
    :class:`NeighbourState` that still answers questions about the destination, so none of the nine is
    expressible however the policy is written."""
    relayed = tuple(np.arange(N, dtype=np.float64) + 100.0 * k for k in range(10))
    src = np.array([0, 0, 1, 2, 3, 4, 5, 6])  # slot i's source, clipped at the terminal
    at_src = SW._at_source(relayed, src)
    for a, b in zip(relayed, at_src):
        assert np.array_equal(b, a[src])
        # the destination's own value is NOT what the policy now holds, except where src == dst
        differs = src != np.arange(N)
        assert np.any(b[differs] != a[differs]), (
            "the gather did not actually move anything (TRAPS: could-the-arm-have-fired)"
        )


def test_the_gather_is_exact_so_assertion_one_costs_no_bits():
    """A gather cannot change a double, so making TRAPS: a-message-from-the-destinations-belief structural is free. If it were not, the assertion
    would be trading correctness for identity and could not ship in a byte-identity commit."""
    rng = np.random.default_rng(7)
    a = rng.standard_normal(500) * 1e17
    src = rng.integers(0, 500, 500)
    assert np.array_equal(SW._at_source((a,), src)[0], a[src])


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 2 — every message mode inside its coordinate's own grid.  TRAPS: off-grid-message-mode, 74 % of a zero control's error.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_a_tilt_mode_in_the_angle_coordinate_is_accepted():
    """The control. ``theta = arcsin(tau)`` spans exactly ``[-pi/2, +pi/2]``, and both endpoints are
    legitimate answers — ``tau = +-1`` is "all RNA on one strand"."""
    tilt = _tilt_grid(60)
    m = np.linspace(tilt[0], tilt[-1], N)
    c = _counts(PsiMessage(theta_mode=m, theta_prec=np.ones(N)))
    assert c["mode_in_grid_theta"]["violations"] == 0
    assert c["mode_in_grid_theta"]["eligible"] == N, (
        "TRAPS: could-the-arm-have-fired: the check must have been eligible somewhere"
    )


def test_PERTURBATION_a_log_odds_delivered_into_the_tilt_slot_is_REFUSED():
    """⛔⛔ **TRAPS: off-grid-message-mode ITSELF, as the exact defect that happened.** A raw log-odds ``log(u+/u-)`` was delivered
    into psi's tilt slot; the measured modes were **+-4.6** against a domain of **+-1.5708**, i.e. 2.9x
    outside the whole coordinate. A Gaussian ``-1/2 p (theta - m)^2`` with ``m`` off-grid is MONOTONE across
    every grid point, so the tilt pinned at the boundary, the AMBIG Schur protection that keeps the strand
    term out of ``f_g`` was destroyed, and the strand likelihood explained the residue by calling the mass
    gDNA. **74 % of the zero control's error, one unit error.**"""
    m = np.full(N, 4.6)
    with pytest.raises(AssertionError, match="mode_in_grid_theta"):
        _counts(PsiMessage(theta_mode=m, theta_prec=np.ones(N)))


def test_a_two_ULP_overshoot_at_the_tilt_endpoint_is_NOT_a_defect():
    """⚠ **MEASURED, and it is why the tolerance is the grid's own SPACING rather than nothing.** The
    shipped tilt mode is a convex mean ``(p_a th_a + p_b th_b)/(p_a + p_b)`` of two messages, and when both
    agree on ``tau = +-1`` the division rounds UP: 63 of 4,795 live slots on one real condition overshoot
    ``pi/2`` by exactly 2 ULP. A bare ``m > hi`` reports a correct answer as a defect.

    ⭐ And the spacing separates the two by four orders of magnitude — TRAPS: off-grid-message-mode's real overshoot is **57
    spacings** — so this is not a threshold buying tolerance, it is the coordinate's own resolution."""
    hi = float(_tilt_grid(60)[-1])
    m = np.full(N, np.nextafter(np.nextafter(hi, np.inf), np.inf))
    assert m[0] > hi
    c = _counts(PsiMessage(theta_mode=m, theta_prec=np.ones(N)))
    assert c["mode_in_grid_theta"]["violations"] == 0


def test_an_off_grid_mode_at_ZERO_precision_is_not_eligible():
    """Only a mode carried at POSITIVE precision can pin anything: psi's term is ``-1/2 p (x - m)^2``, so
    at ``p = 0`` the mode is not read at all. ⭐ That is also what makes the shipped lambda channel's
    out-of-domain mode harmless — the emission gate zeroes its precision."""
    c = _counts(PsiMessage(theta_mode=np.full(N, 40.0), theta_prec=np.zeros(N)))
    assert c["mode_in_grid_theta"] == {"violations": 0, "eligible": 0}


def test_PERTURBATION_a_lambda_mode_beyond_the_log_odds_window_is_COUNTED_and_waived():
    """⛔ The shipped policy violates this one and it is WAIVED WITH A REASON, not widened away.

    The lambda mode is ``log(rho_g E_g) - log(rho_R E_r)`` with both arms floored at ``_EPS``, so a message
    carrying ONE component reads ``+-log(1/_EPS)`` ~ ``+-20.7`` against a grid half-width of 10. The
    emission gate zeroes its PRECISION, which makes it harmless but not in-domain. ⭐ The waiver is what
    keeps the count VISIBLE — the alternative, a looser predicate, would make the check vacuous for the
    coordinate error it exists to catch."""
    assert "mode_in_grid_lam" in SW._KNOWN_VIOLATIONS
    c = _counts(PsiMessage(lam_mode=np.full(N, 20.7), lam_prec=np.ones(N)))
    assert c["mode_in_grid_lam"]["violations"] == N, "a waived assertion must still COUNT"


def test_a_waiver_is_never_silent():
    """Every waived assertion carries a written reason, so a reader learns the defect rather than the
    exemption. ⛔ An empty reason would be a widened predicate wearing a waiver's clothes."""
    for name, why in SW._KNOWN_VIOLATIONS.items():
        assert len(why) > 80, f"{name}'s waiver does not say what the defect is"


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 3 — every delivered share in [0, 1].
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_a_share_below_one_is_accepted():
    """The control: a mode of ``log(share)`` is <= 0 for any real share."""
    c = _counts(PsiMessage(gdna_mode=np.full(N, -0.5), gdna_prec=np.ones(N)))
    assert c["share_in_unit_interval"] == {"violations": 0, "eligible": N}


def test_PERTURBATION_an_over_unit_share_is_COUNTED():
    """⛔ The defect this exists for: a certified-RNA density of 10.3 frag/base over 8,300 b of RNA
    opportunity implied **85,477 fragments at a slot holding 5**. A share above 1 is not a large claim, it
    is an impossible one — a component alone accounting for more fragments than the slot observed.

    ⚠ The shipped policy commits it, so it is WAIVED and counted: the shares are ``rho_c E_c / M`` with no
    upper bound, and the certified-RNA arm is a LOWER BOUND used as an equality. The correction is the
    one-sided term, and landing that ALONE is refused by the panel — TRAPS: a-cancelling-defect-pair, it and the missing
    gDNA level channel are a CANCELLING pair and must be priced in ONE arm."""
    assert "share_in_unit_interval" in SW._KNOWN_VIOLATIONS
    c = _counts(PsiMessage(gdna_mode=np.full(N, np.log(85477.0 / 5.0)), gdna_prec=np.ones(N)))
    assert c["share_in_unit_interval"]["violations"] == N


def test_both_rna_strands_are_checked_not_only_the_first():
    """⚠ TRAPS: off-grid-message-mode's sibling failure: a channel absent from the check is a channel no gate can rank. The RNA
    message is a PAIR, and an over-unit claim on the antisense arm is the same defect."""
    ok = np.full(N, -1.0)
    bad = np.full(N, +2.0)
    c = _counts(PsiMessage(rna_mode=(ok, bad), rna_prec=(np.ones(N), np.ones(N))))
    assert c["share_in_unit_interval"]["violations"] == N


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 4 — |T| <= 3.  AXIOM 0, made executable.
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


def test_PERTURBATION_a_message_carrying_a_fourth_component_channel_is_REFUSED():
    """The same axiom on the message side: the packet carries exactly three component channels — gDNA and
    the two RNA strands. A fourth is not a new feature, it is a violated axiom."""
    three = np.zeros(N)
    with pytest.raises(AssertionError, match="message_has_three_components"):
        _counts(PsiMessage(rna_mode=(three, three, three), rna_prec=(three, three, three)))


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# ASSERTION 5 — the write-back touches only `solvable` slots.
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
    """⛔⛔ **THE LARGEST BEHAVIOUR SWITCH IN THE TOOL, AND IT MUST BE A WRITTEN DECISION.** The two
    policies differ by **99.9 %** on the panel total and by **−58 % to +155 %** depending on the stratum, so
    which one ships can never be inherited from a function default that an edit could silently change.

    ⭐ **It is ON as of 2026-08-18 (owner)**, after ~11 days muted. ⚠ The 2026-08-07 mute was a STUDY
    configuration measured on the 36-condition ladder RETIRED on 2026-08-13; what re-opened it is that on
    a slot with NO own evidence the relay is the only thing that can solve it at all — on unstranded
    capture-OFF those exons read ψ's uninformative ½ EXACTLY (mean |error| 0.500) muted, and 0.000087
    with the relay live. ⛔ It is NOT uniformly better — on stranded CONTAMINATED data it is worse
    whole-chain (`g98 ss0.99 capture_off` 1.628×) — and that asymmetry is the debugging target rather
    than a reason to re-mute."""
    import inspect

    import rigel.calibration.calibrate as _c  # noqa: PLC0415
    from rigel.config import CalibrationConfig  # noqa: PLC0415

    assert CalibrationConfig().message_propagation is True
    src = inspect.getsource(_c)
    assert "config.message_propagation" in src, (
        "calibrate no longer reads the switch — whichever policy it now hard-codes, the config option is "
        "lying to anyone who sets it."
    )
    assert "RelayPolicy(" in src and "SilentPolicy()" in src, (
        "both arms must be reachable from the one call site; a switch with one arm is not a switch."
    )


def test_solve_chains_parameter_default_is_silent_and_sends_nothing():
    """⭐ ``SilentPolicy`` is ``solve_chain``'s PARAMETER default (the shipped config installs the
    relay), and it is a MEASURED floor rather than a placeholder: with no
    belief propagation the deliverable is a net improvement on three of the four strata and a large
    regression on exactly one — the stratum where kappa = 1/2 leaves a slot no own composition evidence."""
    relay = SilentPolicy().prepare(_ctx())
    assert relay.scan(backward=False) is None, "a silent policy must relay nothing at all"
    assert relay.scan(backward=True) is None
    nb = NeighbourState(state=(), valid=np.ones(N, bool), src=np.arange(N))
    assert relay.deliver(nb, nb).is_silent


def test_every_head_operator_is_an_independently_named_switch():
    """⭐⭐ One switch per independently-ablatable operator, because the next step prices them ONE AT A
    TIME. TRAPS: all-small-singly-large-jointly is why: removing them as a block is already measured, and when every single
    ablation is small while the joint one is large you go one stage upstream — you do not keep ablating the
    block."""
    sw = RelaySwitches()
    assert sw.off() == (), "the default must be the shipped answer: every switch ON"
    names = sw.names()
    assert len(names) == len(set(names)) >= 15
    assert not any(c.isupper() for n in names for c in n), "snake_case, and no Greek in identifiers"
    for n in names:
        one_off = RelaySwitches(**{n: False})
        assert one_off.off() == (n,)
        assert RelayPolicy(one_off).name.endswith(f"no_{n}"), "an arm must label itself"


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
