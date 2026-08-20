"""The third policy's falsification gates — written before each concept, verified failing first.

The architecture under test is the owner's (2026-08-19, second ruling pass):

* **sources** are the slots' INITIAL BELIEFS (`ctx.own` — g1 slots, the intron factory, the strand
  model on stranded data); the policy invents no source of its own;
* **the message** is the three ABUNDANCES ``{gDNA, RNA+, RNA-}`` with per-component precisions
  (⛔ the word is ABUNDANCE — counts/bp; "level" is retired vocabulary);
* **the sender sends blindly; the RECIPIENT decides** — per hop the destination asks the static table
  whether its transcript POPULATION equals the source's: a REGION's population is the transcripts
  covering it, a BOUNDARY's is the transcripts CROSSING it (one STARTING there is not in it), so
  ``pop(B) = pop(left) ∩ pop(right)`` and a hop is population-equal iff the REGION side gains nothing
  — one bit per (boundary, side), both directions of a pair sharing it;
* populations equal ⇒ the COMPOSITION strategy (rescale by measured totals, §3.5e); different ⇒ the
  ABUNDANCE strategy (cross unscaled, disagreement-damped at a later concept);
* a slot with NO belief RELAYS the running message unmodified; a slot WITH one fuses
  precision-weighted; at solve time the two directional messages combine and BECOME a beliefless
  slot's belief.

⛔ A refused or inadmissible claim loses its VALUE and its PRECISION in one statement
(``TRAPS: zero-the-precision-with-the-value``).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration import sweep as SW
from rigel.calibration.messages import NeighbourState, PsiMessage, StepContext
from rigel.calibration.messages.currency import (
    CurrencyPolicy,
    hop_structure,
    population_equal_from_left,
    splice_in_densities,
    splice_out_unspliced,
)
from rigel.calibration.region_init import RegionInit
from rigel.calibration.simplex_logodds import _log_fg, _logodds_grid
from rigel.calibration.splice_graph import (
    FLAG_ACCEPTOR_NEG,
    FLAG_DONOR_POS,
    FLAG_TES_NEG,
    FLAG_TSS_POS,
)

N = 9


def _own(rho=None, prec=None) -> RegionInit:
    """A RegionInit carrying only what the policy reads: three abundances + three precisions.

    ⛔ Built through the real constructor, never ``__new__`` — a fixture that bypasses ``__init__``
    tests the rebuild, not the subject (this repo has paid for that once already)."""
    z = np.zeros(N)
    r = {c: np.asarray((rho or {}).get(c, z), np.float64) for c in ("g", "p", "n")}
    p = {c: np.asarray((prec or {}).get(c, z), np.float64) for c in ("g", "p", "n")}
    return RegionInit(
        f_g=np.ones(N),
        f_pos=z.copy(),
        f_neg=z.copy(),
        rho_g=r["g"],
        rho_pos=r["p"],
        rho_neg=r["n"],
        prec_g=p["g"],
        prec_pos=p["p"],
        prec_neg=p["n"],
        struct_lock=np.zeros(N, bool),
        tau_lam=z.copy(),
    )


def _ctx(
    *,
    own: RegionInit | None = None,
    mass=None,
    eff=None,
    eff_rna=None,
    inv_abundance=None,
    n_slot=None,
    boundary_flags=None,
    free_pos=None,
    free_neg=None,
    n_grid=60,
) -> StepContext:
    """A minimal chain ``N E N E … N`` (9 slots). Only what the policy reads needs to be real."""
    ones = np.ones(N)
    m = ones * 100.0 if mass is None else np.asarray(mass, np.float64)
    # ⚠ ONE divisor for both components unless a test says otherwise: `eff` sets E_g AND E_r, because
    # a fixture that silently left E_r at another value made the k-form's rescale read a frame gap
    # that the scenario did not intend (found while gating §3.5e's worked numbers).
    e = ones * 200.0 if eff is None else np.asarray(eff, np.float64)
    e_r = e if eff_rna is None else np.asarray(eff_rna, np.float64)
    fp = np.ones(N, bool) if free_pos is None else np.asarray(free_pos, bool)
    fn = np.zeros(N, bool) if free_neg is None else np.asarray(free_neg, bool)
    if boundary_flags is None:
        # a real chain's every BOUNDARY carries at least one structure bit (measured: 0 of 1,043,595
        # human boundaries carry neither); the default is a bare + terminus
        fl = np.where(np.arange(N) % 2 == 1, FLAG_TSS_POS, 0).astype(np.uint16)
    else:
        fl = np.asarray(boundary_flags, np.uint16)
    _, grid = _logodds_grid(n_grid, 10.0)
    return StepContext(
        mass=m,
        # the MODEL-FREE abundance: at M == E the reciprocal-opportunity bank reads M/E; a fixture may
        # override it to state an enrichment ratio directly (that is what the knob reads)
        inv_abundance=m / np.where(e > 0, e, 1.0)
        if inv_abundance is None
        else np.asarray(inv_abundance, np.float64),
        eff_gdna_global=e,
        eff_rna=e_r,
        eff_gdna=e,
        eff_sj=np.ones((N, 2)) * 200.0,
        sj_count=np.zeros((N, 2)),
        unspliced_count=np.stack([m / 2.0, m / 2.0], axis=1),
        # ⚠ the observed COUNTS are what set the knob's noise floor, and they are decoupled from the
        # abundances here on purpose: a scenario may hold the same densities at any depth, and a gate
        # that wants the knob at its COMPOSITION end says so by giving it counts it can trust.
        n_slot=m if n_slot is None else np.asarray(n_slot, np.float64),
        spliced_slot=np.zeros(N),
        left=np.arange(-1, N - 1),
        right=np.append(np.arange(1, N), -1),
        is_boundary=np.arange(N) % 2 == 1,
        is_exon_region=np.arange(N) % 2 == 0,
        free_pos=fp,
        free_neg=fn,
        g1_locked=~fp & ~fn,
        boundary_flags=fl,
        geometry=None,
        order=list(range(N)),
        left_list=list(range(-1, N - 1)),
        right_list=[*range(1, N), -1],
        own=own if own is not None else _own(),
        belief_fg=ones,
        n_grid=n_grid,
        logodds_window=10.0,
        solve_grid=grid,
    )


def _deliver(ctx: StepContext) -> PsiMessage:
    """Run the policy exactly as the backbone would: two scans, gather at source, one deliver."""
    relay = CurrencyPolicy().prepare(ctx)
    states = {}
    for backward, seq, nbr in (
        (False, ctx.order, ctx.left_list),
        (True, ctx.order[::-1], ctx.right_list),
    ):
        kernel = relay.scan(backward=backward)
        assert kernel is not None, "the policy relays abundances, so the scan must not be None"
        step, publish = kernel
        for i in seq:
            s = nbr[i]
            if s >= 0:
                step(s, i)
        states[backward] = publish()
    left = np.asarray(ctx.left, np.int64)
    right = np.asarray(ctx.right, np.int64)
    sl, sr = np.clip(left, 0, N - 1), np.clip(right, 0, N - 1)
    msg = relay.deliver(
        NeighbourState(state=SW._at_source(states[False], sl), valid=left >= 0, src=sl),
        NeighbourState(state=SW._at_source(states[True], sr), valid=right >= 0, src=sr),
    )
    # the message must pass the BACKBONE's own assertions — that is the point of the shared seam
    counts = SW.AssertionCounts()
    SW._check_message(msg, ctx, counts)
    return msg


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# the protocol, and the silent floor
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_policy_satisfies_the_backbone_protocol():
    from rigel.calibration.messages import Policy

    assert isinstance(CurrencyPolicy(), Policy)
    assert CurrencyPolicy().name == "currency"


def test_no_belief_anywhere_delivers_no_claim():
    """With every slot beliefless there is nothing to relay, and 'no claim' must be a PRECISION of
    exactly 0 on every arm — not a weak mode."""
    msg = _deliver(_ctx(own=_own()))
    for p in (msg.gdna_prec, *(msg.rna_prec or (np.zeros(N), np.zeros(N)))):
        assert p is not None and np.all(np.asarray(p) == 0.0)
    assert msg.lam_mode is None and msg.theta_mode is None


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# CONCEPT A — population equality per (boundary, side), the owner's worked example verbatim
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_population_equality_matches_the_owners_terminus_example():
    """TA+ (1000,10000), TB+ (5000,10000), TC− (500,20000); B@5000 is TB's TSS.

    pop(R1 = (1000,5000)) = {TA+, TC−}; pop(B) = {TA+, TC−} (TB STARTS at B and is not in it);
    pop(R2 = (5000,10000)) = {TA+, TB+, TC−}. So the hop R1↔B is population-EQUAL — composition is
    licensed even though B carries a TSS bit — and B↔R2 is NOT (R2 gains TB+)."""
    # chain: … R1(slot 2)  B(slot 3, TSS_POS)  R2(slot 4) …
    fl = np.zeros(N, np.uint16)
    fl[[1, 5, 7]] = FLAG_TES_NEG  # the other boundaries need some bit; irrelevant to the example
    fl[3] = FLAG_TSS_POS
    ctx = _ctx(boundary_flags=fl, free_neg=np.ones(N, bool))
    eq = population_equal_from_left(ctx)
    # the hop INTO slot 3 from its left (R1 → B): populations equal — licensed
    assert eq[3], "R1 -> B must be population-equal: TB+ STARTS at B and is not in B's population"
    # the hop INTO slot 4 from its left (B → R2): R2 gains TB+ — not equal
    assert not eq[4], "B -> R2 must NOT be population-equal: R2 gains TB+"


def test_population_equality_is_a_pair_property():
    """The pair (i, left[i]) IS the pair (left[i], right[left[i]]): both directions of one hop share
    the bit, so the right-hand table is the left-hand one read through ``right``."""
    fl = np.zeros(N, np.uint16)
    fl[[1, 3, 5, 7]] = FLAG_TSS_POS
    eq_l = population_equal_from_left(_ctx(boundary_flags=fl))
    # at a TSS_POS boundary the RIGHT flank gains, so: the hop from the LEFT region into the boundary
    # is equal (eq at the boundary slot), and the hop from the boundary into the RIGHT region is not
    # (eq at the region slot reads the SAME boundary's right-gain through `left`) — one bit, two reads
    assert eq_l[1] and not eq_l[2] and eq_l[3] and not eq_l[4]
    # a slot with no left neighbour has no hop to be equal on
    assert not eq_l[0]


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# CONCEPT B — sources are the initial beliefs; relay-unmodified through beliefless slots; fuse at
# believing ones; a refused arm loses value and precision together
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_a_beliefless_chain_relays_a_claim_unmodified():
    """One believing slot (0); everything else beliefless. The claim must arrive at slot 4 with its
    VALUE and PRECISION exactly as sent — a beliefless slot relays unmodified (owner's lifecycle)."""
    rho = {"g": np.array([0.3, 0, 0, 0, 0, 0, 0, 0, 0])}
    prec = {"g": np.array([0.7, 0, 0, 0, 0, 0, 0, 0, 0])}
    # ⚠ every hop population-equal (no gains): TES_NEG makes the LEFT flank gain − RNA… use a bare
    # DONOR_POS: a splice site changes no population by itself. M = E so the totals are uniform
    # (r = 1) and the implied counts are large (the sampling cap must not bind on a physical claim).
    fl = np.where(np.arange(N) % 2 == 1, FLAG_DONOR_POS, 0).astype(np.uint16)
    msg = _deliver(
        _ctx(own=_own(rho, prec), mass=np.full(N, 1.0e4), eff=np.full(N, 1.0e4), boundary_flags=fl)
    )
    np.testing.assert_allclose(np.asarray(msg.gdna_prec)[4], 0.7, rtol=0, atol=0)
    assert np.asarray(msg.gdna_mode)[4] == pytest.approx(np.log(0.3))


def test_a_believing_slot_dominates_a_weak_message():
    """A strong local belief at slot 2 (precision 100) against a weak incoming claim (precision 0.1):
    what continues past slot 2 is essentially slot 2's own — 'the message will likely dampen or die
    there' (owner's lifecycle)."""
    rho = {"g": np.array([0.09, 0, 0.01, 0, 0, 0, 0, 0, 0])}
    prec = {"g": np.array([0.1, 0, 100.0, 0, 0, 0, 0, 0, 0])}
    fl = np.where(np.arange(N) % 2 == 1, FLAG_DONOR_POS, 0).astype(np.uint16)
    msg = _deliver(
        _ctx(own=_own(rho, prec), mass=np.full(N, 1.0e4), eff=np.full(N, 1.0e4), boundary_flags=fl)
    )
    # delivered at slot 4: the fused running state from the left = (0.1*0.09 + 100*0.01)/100.1
    got = np.exp(np.asarray(msg.gdna_mode)[4])  # share == abundance at M == E
    np.testing.assert_allclose(got, (0.1 * 0.09 + 100.0 * 0.01) / 100.1, rtol=1e-12)
    np.testing.assert_allclose(np.asarray(msg.gdna_prec)[4], 100.1, rtol=0, atol=0)


def test_zero_is_carried_as_a_value_with_its_own_precision():
    """An empty g1 anchor's claim is ``rho_g = 0`` at finite precision. The zero must ARRIVE as the
    lowest expressible share at that precision — never dropped, never floored away
    (``TRAPS: a-ratio-cannot-carry-zero``)."""
    rho = {"g": np.zeros(N)}
    prec = {"g": np.array([0.2026, 0, 0, 0, 0, 0, 0, 0, 0])}
    # ⭐ the sampling cap on a ZERO-count claim is 1/(trigamma(½) − trigamma(1.5)) = 0.25, ABOVE the honest
    # zero-claim precision 1/trigamma(½) ≈ 0.2026 — so the cap never punishes an honest zero
    fl = np.where(np.arange(N) % 2 == 1, FLAG_DONOR_POS, 0).astype(np.uint16)
    msg = _deliver(_ctx(own=_own(rho, prec), boundary_flags=fl))
    lam, _ = _logodds_grid(60, 10.0)
    assert np.asarray(msg.gdna_prec)[4] == pytest.approx(0.2026)
    assert np.asarray(msg.gdna_mode)[4] == pytest.approx(float(_log_fg(lam)[0]))


def test_an_inadmissible_strand_loses_value_and_precision_together():
    """An RNA+ claim relayed into a chain segment whose annotation admits no + strand must arrive as
    NOTHING — value 0 AND precision 0, one statement (``TRAPS: zero-the-precision-with-the-value``)."""
    rho = {"p": np.array([5.0, 0, 0, 0, 0, 0, 0, 0, 0])}
    prec = {"p": np.array([2.0, 0, 0, 0, 0, 0, 0, 0, 0])}
    # ⚠ the claim must CROSS the inadmissible slot (2) and land on an ADMISSIBLE one (4) — otherwise
    # the delivery-time admissibility belt masks a scan that kept the precision (measured: it did)
    fp = np.array([True, True, False, True, True, True, True, True, True])
    fl = np.where(np.arange(N) % 2 == 1, FLAG_TES_NEG, 0).astype(np.uint16)
    msg = _deliver(
        _ctx(own=_own(rho, prec), free_pos=fp, free_neg=np.ones(N, bool), boundary_flags=fl)
    )
    assert msg.rna_prec is not None
    assert np.asarray(msg.rna_prec[0])[4] == 0.0
    assert np.asarray(msg.rna_mode[0])[4] == 0.0


def test_both_directions_fuse_precision_weighted_at_delivery():
    """Believing slots at both ends; the middle slot's delivered claim is the precision-weighted
    combination of the two directional messages (owner: 'the average of the precision-weighted
    messages becomes the belief')."""
    rho = {"g": np.array([0.02, 0, 0, 0, 0, 0, 0, 0, 0.08])}
    prec = {"g": np.array([3.0, 0, 0, 0, 0, 0, 0, 0, 1.0])}
    fl = np.where(np.arange(N) % 2 == 1, FLAG_DONOR_POS, 0).astype(np.uint16)
    msg = _deliver(
        _ctx(own=_own(rho, prec), mass=np.full(N, 1.0e4), eff=np.full(N, 1.0e4), boundary_flags=fl)
    )
    got = np.exp(np.asarray(msg.gdna_mode)[4])
    np.testing.assert_allclose(got, (3.0 * 0.02 + 1.0 * 0.08) / 4.0, rtol=1e-12)
    np.testing.assert_allclose(np.asarray(msg.gdna_prec)[4], 4.0, rtol=0, atol=0)


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# the static structure bits (kept from rung 1 — concept C's inputs)
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_hop_table_reads_term_and_sj_off_the_flags_per_strand():
    kind_is_boundary = np.array([False, True, False, True, False], bool)
    flags = np.array([0, FLAG_TSS_POS, 0, FLAG_DONOR_POS | FLAG_TES_NEG, 0], np.uint16)
    st = hop_structure(kind_is_boundary, flags)
    assert st.term_pos[1] and not st.sj_pos[1] and not st.term_neg[1] and not st.sj_neg[1]
    assert st.sj_pos[3] and not st.term_pos[3] and st.term_neg[3] and not st.sj_neg[3]
    assert not st.term_pos[0] and not st.sj_pos[0]


def test_a_boundary_with_neither_bit_is_REFUSED():
    kind_is_boundary = np.array([False, True, False], bool)
    with pytest.raises(AssertionError):
        hop_structure(kind_is_boundary, np.zeros(3, np.uint16))


def test_sj_plus_term_rules_with_term():
    """Owner: a terminus AT a splice junction is still a terminus — the population still changes."""
    kind_is_boundary = np.array([False, True, False], bool)
    flags = np.array([0, FLAG_TSS_POS | FLAG_ACCEPTOR_NEG, 0], np.uint16)
    st = hop_structure(kind_is_boundary, flags)
    assert st.term_pos[1] and st.sj_neg[1] and not st.sj_pos[1] and not st.term_neg[1]


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# the §3.5e operators — the owner's worked numbers, verbatim (wired at concept C)
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_worked_SPLICE_OUT_numbers():
    out = splice_out_unspliced(src=(10.0, 90.0, 0.0), dst_unspliced=(2.0, 1.0, 0.0), sj=(17.0, 0.0))
    np.testing.assert_allclose(out, (2.0, 1.0, 0.0), rtol=0, atol=1e-12)


def test_the_worked_SPLICE_IN_inverse():
    out = splice_in_densities(src_with_flux=(2.0, 18.0, 0.0), rho_tot_dst=100.0, rho_tot_src=20.0)
    np.testing.assert_allclose(out, (10.0, 90.0, 0.0), rtol=0, atol=1e-12)


def test_the_splice_out_refuses_a_zero_total():
    out = splice_out_unspliced(src=(0.0, 0.0, 0.0), dst_unspliced=(0.0, 0.0, 0.0), sj=(0.0, 0.0))
    assert out is None


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# the config seam — a third value beside the two shipped policies, default unchanged
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_policy_is_selectable_by_config_and_the_default_is_the_relay():
    import inspect

    from rigel.calibration.calibrate import calibrate
    from rigel.config import CalibrationConfig

    assert CalibrationConfig().message_policy == "relay"
    src = inspect.getsource(calibrate)
    for token in ("RelayPolicy()", "SilentPolicy()", "CurrencyPolicy()", "message_policy"):
        assert token in src, f"calibrate no longer wires {token}"


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# CONCEPT C — the COMPOSITION strategy on population-equal hops: belief-free rescale by measured
# totals, the §3.5e flux arithmetic, and the binomial up/down-sampling transfer variance (the sampling variance)
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def _geom(flux_lo=None, flux_hi=None):
    """A fake geometry carrying only the four flank-flux banks the policy reads (per TRANSCRIPT
    strand). ``None`` ⇒ no flux anywhere."""
    from types import SimpleNamespace

    z = np.zeros((N, 2))
    e = np.ones((N, 2))
    return SimpleNamespace(
        sj_count_lo=z if flux_lo is None else np.asarray(flux_lo, np.float64),
        sj_count_hi=z if flux_hi is None else np.asarray(flux_hi, np.float64),
        eff_sj_lo=e,
        eff_sj_hi=e,
    )


def test_a_composition_hop_rescales_by_the_measured_totals():
    """Source region T = 0.5 (M=100/E=200), destination boundary T = 1.0 (M=100/E=100): the claim
    doubles — the rescale is the MEASURED total ratio, no belief anywhere in it (owner ruling)."""
    # the source supplies BOTH live components (the licence's SUPPLY half); one frame, E = 200
    eff = np.full(N, 200.0)
    mass = np.array([100.0, 40.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0])
    rho = {
        "g": np.array([0.1, 0, 0, 0, 0, 0, 0, 0, 0]),
        "p": np.array([0.4, 0, 0, 0, 0, 0, 0, 0, 0]),
    }
    prec = {
        "g": np.array([500.0, 0, 0, 0, 0, 0, 0, 0, 0]),  # over-confident ON PURPOSE — see the cap
        "p": np.array([500.0, 0, 0, 0, 0, 0, 0, 0, 0]),
    }
    fl = np.where(np.arange(N) % 2 == 1, FLAG_DONOR_POS, 0).astype(np.uint16)
    # the two slots' MODEL-FREE totals differ by 2x — a real enrichment — and the counts are deep, so
    # the knob reaches its COMPOSITION end (w -> 1) and the full mass-identity rescale applies
    inv = np.array([0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    ctx = _ctx(
        own=_own(rho, prec),
        mass=mass,
        eff=eff,
        boundary_flags=fl,
        inv_abundance=inv,
        n_slot=np.full(N, 1.0e9),
    )
    object.__setattr__(ctx, "geometry", _geom())
    msg = _deliver(ctx)
    # k = M_dst / S with S = (0.1 + 0.4)*200 = 100 and M_dst = 40 -> k = 0.4; the gDNA claim becomes
    # 0.04, i.e. a share of 0.04*200/40 = 0.2 at the destination
    np.testing.assert_allclose(np.exp(np.asarray(msg.gdna_mode)[1]), 0.2, rtol=1e-8)
    # and the precision is CAPPED at what the source counts support ("still based on 2 fragments"):
    # n_g = 0.1*200 = 20 of n_tot = 100 -> cap = 1/(trigamma(20.5) - trigamma(101.5)), below 500
    from rigel.calibration.messages.variance import count_logvar

    cap = 1.0 / (count_logvar(20.0) - count_logvar(101.0))
    np.testing.assert_allclose(np.asarray(msg.gdna_prec)[1], cap, rtol=1e-6)


def test_upsampling_adds_no_precision():
    """The owner's 122x example: multiplying a 2-fragment claim by a large factor must never RAISE its
    precision — the sampling variance is based on the source counts alone."""
    from rigel.calibration.messages.currency import composition_sampling_logvar

    v_small = composition_sampling_logvar(n_c=2.0, n_tot=10.0)
    v_large = composition_sampling_logvar(n_c=2000.0, n_tot=10000.0)
    assert v_small > v_large > 0.0
    # and a zero-count component gets a WIDE but finite variance (Jeffreys), never a refusal
    v_zero = composition_sampling_logvar(n_c=0.0, n_tot=10.0)
    assert np.isfinite(v_zero) and v_zero > v_small
    # ⛔ the variance is a property of the SHARE, so it must depend on the TOTAL too: 2 fragments out
    # of 2 is a nearly certain share, 2 out of 10 is not. A form that reads only ``n_c`` would pass
    # every assertion above and be wrong (measured — this clause is what fires that perturbation).
    assert composition_sampling_logvar(n_c=2.0, n_tot=2.0) < 0.5 * v_small


def test_the_splice_flux_arithmetic_matches_the_worked_numbers_end_to_end():
    """§3.5e(i) wired: an exon's claim entering an sj boundary is rescaled by the flux-inclusive total
    ratio and the flux is subtracted — the owner's {10, 90, 0} -> {2, 1, 0}, through the actual scan."""
    # slot 0 = the exon (source), slot 1 = the boundary. Totals: T_exon = 100, T_boundary = 3 + 17.
    # Engineer: exon M=100, E_g=1 -> T=100. boundary M=3, E_g=1 -> unspliced 3; flux_lo 17 on +.
    mass = np.array([100.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    eff = np.ones(N)
    rho = {
        "g": np.array([10.0, 0, 0, 0, 0, 0, 0, 0, 0]),
        "p": np.array([90.0, 0, 0, 0, 0, 0, 0, 0, 0]),
    }
    prec = {
        "g": np.array([5.0, 0, 0, 0, 0, 0, 0, 0, 0]),
        "p": np.array([5.0, 0, 0, 0, 0, 0, 0, 0, 0]),
    }
    fl = np.where(np.arange(N) % 2 == 1, FLAG_DONOR_POS, 0).astype(np.uint16)
    flux_lo = np.zeros((N, 2))
    flux_lo[1, 0] = 17.0
    # the owner's example IS an enrichment: the exon's total is 100 and the boundary's 20. State that
    # in the model-free bank and give it deep counts, so the knob reaches w = 1 and the arithmetic is
    # exactly §3.5e's.
    inv = np.full(N, 20.0)
    inv[0] = 100.0
    ctx = _ctx(
        own=_own(rho, prec),
        mass=mass,
        eff=eff,
        boundary_flags=fl,
        inv_abundance=inv,
        n_slot=np.full(N, 1.0e9),
    )
    object.__setattr__(ctx, "geometry", _geom(flux_lo=flux_lo))
    relay = CurrencyPolicy().prepare(ctx)
    step, publish = relay.scan(backward=False)
    step(0, 1)
    vg, vp, vn, pg, pp, pn = publish()
    # ⚠ to the KNOB's own noise floor, not the ULP: ``w = lr²/(lr²+v)`` is 1 − O(v) and ``v > 0`` at
    # any finite depth, so 1e9 counts leave ~1e-9. Demanding exactness here would demand ``v = 0``,
    # i.e. infinite counts — the tolerance IS the derivation.
    np.testing.assert_allclose(vg[1], 2.0, rtol=1e-7)
    np.testing.assert_allclose(vp[1], 1.0, rtol=1e-7)


def test_the_splice_in_flux_joins_the_message_out_of_a_boundary():
    """§3.5e(ii) wired: a boundary's claim entering its exon carries the flux WITH it, then rescales —
    {2, 1, 0} + flux 17 at totals 20 -> 100 becomes {10, 90, 0}."""
    mass = np.array([1.0, 3.0, 100.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    eff = np.ones(N)
    rho = {
        "g": np.array([0, 2.0, 0, 0, 0, 0, 0, 0, 0]),
        "p": np.array([0, 1.0, 0, 0, 0, 0, 0, 0, 0]),
    }
    prec = {
        "g": np.array([0, 5.0, 0, 0, 0, 0, 0, 0, 0]),
        "p": np.array([0, 5.0, 0, 0, 0, 0, 0, 0, 0]),
    }
    fl = np.where(np.arange(N) % 2 == 1, FLAG_DONOR_POS, 0).astype(np.uint16)
    # the boundary's HI face carries the flux (the hop into the RIGHT region is the splice-in)
    flux_hi = np.zeros((N, 2))
    flux_hi[1, 0] = 17.0
    inv = np.full(N, 20.0)
    inv[2] = 100.0
    ctx = _ctx(
        own=_own(rho, prec),
        mass=mass,
        eff=eff,
        boundary_flags=fl,
        inv_abundance=inv,
        n_slot=np.full(N, 1.0e9),
    )
    object.__setattr__(ctx, "geometry", _geom(flux_hi=flux_hi))
    relay = CurrencyPolicy().prepare(ctx)
    step, publish = relay.scan(backward=False)
    step(0, 1)  # beliefless region 0 into the boundary: nothing arrives (prec 0 either side)
    step(1, 2)  # the boundary's own belief splices IN to region 2
    vg, vp, vn, pg, pp, pn = publish()
    np.testing.assert_allclose(vg[2], 10.0, rtol=1e-7)
    np.testing.assert_allclose(vp[2], 90.0, rtol=1e-7)


def test_a_population_unequal_hop_does_not_rescale():
    """The ABUNDANCE strategy: across a terminus the claim crosses UNSCALED whatever the totals say —
    the disagreement becomes a VARIANCE (concept D), never a rescale."""
    eff = np.array([200.0, 100.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0])
    rho = {"g": np.array([0.04, 0, 0, 0, 0, 0, 0, 0, 0])}
    prec = {"g": np.array([50.0, 0, 0, 0, 0, 0, 0, 0, 0])}
    fl = np.where(np.arange(N) % 2 == 1, FLAG_TSS_POS, 0).astype(np.uint16)  # termini everywhere
    ctx = _ctx(own=_own(rho, prec), eff=eff, boundary_flags=fl)
    object.__setattr__(ctx, "geometry", _geom())
    relay = CurrencyPolicy().prepare(ctx)
    step, publish = relay.scan(backward=False)
    step(
        0, 1
    )  # region -> boundary: the LEFT flank gains nothing at a TSS_POS, so this hop IS equal
    step(1, 2)  # boundary -> region: the region GAINS the starting transcript — the unequal hop
    vg = publish()[0]
    np.testing.assert_allclose(vg[2], vg[1], rtol=0, atol=0)  # unscaled across the unequal hop


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# CONCEPT D — THE KNOB. The two strategies are the two ENDS of one continuum, and the data picks the
# point (owner's intuition, 2026-08-20: "there is actually a knob that connects them").
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_enrichment_ratio_is_model_free():
    """⛔ The enrichment ratio must come from the RECIPROCAL-OPPORTUNITY bank, never from
    ``mass / effective_length``: that divisor depends on the composition being solved for, so a total
    abundance built from it swings with the gDNA-vs-RNA length gap (the owner's 0.25-vs-0.33 example).

    Here the two slots hold the SAME mass and the SAME model-free abundance but very different
    effective lengths — a length-model difference alone. A composition-blind ratio reads 1.0; a
    ``mass/E`` ratio would read 4."""
    from rigel.calibration.messages.currency import enrichment_ratio

    eff = np.array([50.0, 200.0, 50.0, 200.0, 50.0, 200.0, 50.0, 200.0, 50.0])
    ctx = _ctx(mass=np.full(N, 100.0), eff=eff, inv_abundance=np.full(N, 2.0))
    r = enrichment_ratio(ctx, backward=False)
    np.testing.assert_allclose(r[1:], 1.0, rtol=0, atol=0)


def test_the_knob_is_zero_when_the_disagreement_is_noise():
    """No enrichment ⇒ the observed log-ratio is sampling noise ⇒ the knob shrinks it to ~0 and the
    claim crosses as an ABUNDANCE. That is the capture-OFF end of the continuum."""
    from rigel.calibration.messages.currency import rescale_weight

    w = rescale_weight(log_ratio=0.0, var_ratio=0.5)
    assert w == 0.0
    w_small = rescale_weight(log_ratio=0.05, var_ratio=0.5)
    assert 0.0 <= w_small < 0.02


def test_the_knob_is_one_when_the_disagreement_dwarfs_the_noise():
    """A 122x enrichment against ordinary counting noise ⇒ the knob goes to ~1 and the claim crosses
    as a COMPOSITION. That is the capture-ON end of the same continuum."""
    from rigel.calibration.messages.currency import rescale_weight

    w = rescale_weight(log_ratio=np.log(122.0), var_ratio=0.05)
    assert w > 0.99


def test_the_knob_is_continuous_and_monotone():
    """⭐ There is no switch: the weight rises monotonically with the disagreement and falls with the
    noise, so every point between "pure abundance" and "pure composition" is reachable and the DATA
    chooses it."""
    from rigel.calibration.messages.currency import rescale_weight

    ws = [rescale_weight(log_ratio=x, var_ratio=0.2) for x in (0.0, 0.1, 0.3, 1.0, 3.0)]
    assert ws == sorted(ws) and ws[0] == 0.0 and ws[-1] > 0.9
    noisier = [rescale_weight(log_ratio=1.0, var_ratio=v) for v in (0.01, 0.1, 1.0, 10.0)]
    # ⛔ STRICTLY decreasing: a form that ignores the noise entirely produces an all-equal list, which
    # a non-strict check accepts (measured — that perturbation passed until this clause was sharpened).
    assert all(a > b for a, b in zip(noisier, noisier[1:])), (
        "more noise must move the knob toward the ABUNDANCE end, strictly"
    )


def test_a_partial_claim_never_reaches_the_composition_end():
    """The SUPPLY half of the licence is not overridden by a large disagreement: without every
    admissible component the mass identity has no right-hand side, so the knob stays at the ABUNDANCE
    end however loud the enrichment is (this is what stopped 73,728 invented gDNA fragments)."""
    inv = np.full(N, 1.0)
    inv[1] = 100.0  # a huge disagreement at the boundary
    rho = {
        "g": np.array([0.2, 0, 0, 0, 0, 0, 0, 0, 0])
    }  # gDNA only — RNA+ is admissible and UNSUPPLIED
    prec = {"g": np.array([9.0, 0, 0, 0, 0, 0, 0, 0, 0])}
    fl = np.where(np.arange(N) % 2 == 1, FLAG_DONOR_POS, 0).astype(np.uint16)
    ctx = _ctx(
        own=_own(rho, prec),
        eff=np.full(N, 100.0),
        mass=np.full(N, 100.0),
        inv_abundance=inv,
        boundary_flags=fl,
    )
    object.__setattr__(ctx, "geometry", _geom())
    relay = CurrencyPolicy().prepare(ctx)
    step, publish = relay.scan(backward=False)
    step(0, 1)
    np.testing.assert_allclose(publish()[0][1], 0.2, rtol=0, atol=0)


def test_a_large_disagreement_damps_an_unsupplied_claim():
    """The other half of the owner's rule: where the claim cannot be rescaled, a large disagreement
    must make it LESS trustworthy — "high disagreement -> high variance -> low precision -> the
    destination ignores the message"."""
    rho = {"g": np.array([0.2, 0, 0, 0, 0, 0, 0, 0, 0])}
    prec = {"g": np.array([9.0, 0, 0, 0, 0, 0, 0, 0, 0])}
    fl = np.where(np.arange(N) % 2 == 1, FLAG_DONOR_POS, 0).astype(np.uint16)

    def _prec_at_1(inv1):
        inv = np.full(N, 1.0)
        inv[1] = inv1
        ctx = _ctx(
            own=_own(rho, prec),
            eff=np.full(N, 100.0),
            mass=np.full(N, 100.0),
            inv_abundance=inv,
            boundary_flags=fl,
        )
        object.__setattr__(ctx, "geometry", _geom())
        relay = CurrencyPolicy().prepare(ctx)
        step, publish = relay.scan(backward=False)
        step(0, 1)
        return publish()[3][1]

    agree = _prec_at_1(1.0)
    disagree = _prec_at_1(100.0)
    assert disagree < 0.5 * agree, "a 100x disagreement did not damp the claim"
