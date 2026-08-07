#!/usr/bin/env python
"""⭐⭐⭐ THE ``eta`` REBUILD, AS A STANDALONE ``node_sweep`` — the structural simplification, whole.

⛔ **A PROTOTYPE.** It is injected in place of `bp_solver.node_sweep` by ``ladder_arm_ab.py --arm eta`` so
the framework can be measured **as a whole** rather than as a sequence of partial edits each fighting the
operators that remain — the obstruction that made the preceding campaign hard to read. When it earns its
place it moves into ``src/`` and the old core is deleted in one commit.

## The derivation, in one paragraph

The transferable invariant is ``eta = lambda - log(E_g/E_r) = log(rho_g/rho_R)``. It crosses any
same-population EDGE as the **IDENTITY**, so the destination's own opportunities do the conversion exactly
and belief-free, and ``sum_c rho_c E_c = M`` holds by construction. A composition may cross iff
``T_s == T_d`` — a two-bit comparison of ``(free_pos, free_neg)`` against the two opportunities. Where the
populations differ, the shared densities cross unchanged and the newly-active population takes the residual
of the destination's OWN count: determined, except for the tilt when both strands are newly live.

The algebra lives in ``tests/calibration/_eta_reference.py`` — ONE home, shared with the gates in
``test_eta_transfer.py``, so the prototype and its falsification cannot drift apart (`TRAPS.md` A11).

## What this deletes, and why each one goes

======================================  ===================================================================
deleted                                 because
======================================  ===================================================================
the **mass pin**, both licence cases    the identity holds by construction, so there is nothing to restore.
(``_lend``/``g1_locked``) + ``_pin_v``  ⭐ An operator with a measured failure in BOTH directions.
the reframe ``r``, ``framed``, and      nothing forms a ratio of totals; the conversion is a known
the flank pair ``rho_lo``/``rho_hi``    geometric constant read off the index.
four ``rho > _EPS`` guards              ``psi0(a+1/2) - log E`` is finite at every ``a >= 0``, so no exact
                                        zero is ever formed and none of them has anything to guard.
``pop`` / ``lend`` / graft-peel as      one ``T`` comparison and one set intersection.
SEPARATE gating concepts
``transfer_logvar``'s non-graft branch  a belief-free exact conversion has no scale noise to charge.
``struct_lock`` as a flag distinct      ``|T| == 1`` IS structural certainty, derived here rather than
from ``g1_locked``                      carried as a second predicate.
======================================  ===================================================================

## ⚠ What it does NOT change, deliberately

* **ψ is untouched.** The sweep hands ``_solve_nodes_logodds_all`` the ``lam_imp``/``theta_imp`` pair it
  already accepts, so the per-node solve, the strand likelihood, the reference and the fitted gDNA prior
  are all the shipped ones. Only the MESSAGES differ, which is what makes the arm readable.
* **The self-solve is the shipped `build_node_init`.** One thing varied (`TRAPS.md` B8).
* **The spliced bank is untouched.** A spliced fragment is certified RNA — a MEASUREMENT of the RNA side,
  needing no deconvolution. The identity here is over the unspliced bank.

⛔ **ACCEPTANCE IS ACCURACY-NEUTRALITY ON THE DELIVERABLE (`abs_err_all_final`), NOT A WIN.** A -37.2 %
pass-0 change was -3.9 % shipped (`TRAPS.md` A15), and nine pass-0 mechanisms from three independent
mathematical directions were priced on the full panel and refused by the zero control. Making the stranded
regression this arm's bar would repeat that campaign's central mistake.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests" / "calibration"))

from _eta_reference import (  # noqa: E402
    GDNA,
    log_opportunity_ratio,
    log_rate_posterior,
    population_set,
    residual_eta_scalar,
    share_from_density,
    tilt_angle,
)
from rigel.calibration.bp_solver import _logodds_grid, _solve_nodes_logodds_all  # noqa: E402
from rigel.calibration.node_chain import NODE  # noqa: E402
from rigel.calibration.node_geometry import NodeBelief, node_global_geometry  # noqa: E402
from rigel.calibration.node_init import build_node_init  # noqa: E402

_EPS = 1.0e-12


def _exon_mask(chain, region_arrays):
    """Per-slot "is an EXON NODE" — only an exon receives a junction's mature flux (an intron carries
    none, and grafting into one was a measured defect in the predecessor)."""
    from rigel.calibration.signature import coarse_type_array

    rt = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    ri = np.clip(np.asarray(chain.obj_idx, dtype=np.int64), 0, rt.shape[0] - 1)
    return (np.asarray(chain.kind) == NODE) & (rt[ri] == 2)

#: incremented every time the prototype's sweep runs — `TRAPS.md` A10, an arm that never fired reads as
#: "no effect" and is indistinguishable from a clean negative result.
FIRED = [0]

#: ⭐ How often a measurement channel's implied fragment count EXCEEDED the destination's own mass, and
#: how much mass sat at those slots. The graft's premise failing is a measurement in its own right and
#: nothing else in the tool reports it; a silent refusal would be `TRAPS.md` A14 with the sign flipped.
REFUSED: dict = {}

#: ⭐ THE ARM'S OWN ATTRIBUTION HANDLE. With the LEVEL channel muted, no ``|N| >= 1`` step can fire and the
#: sweep carries the same-population ``eta`` identity ALONE. That separates the two halves of the
#: derivation — §4's frame-free composition transfer from §2+§5's "a zero count is a positive density,
#: and the destination's residual converts it into a positive ``f_g``" — which is the pairing the zero
#: control is the detector for. Set by ``ladder_arm_ab.py --arm eta_nolevel``.
MUTE_LEVEL = [False]


def _fuse(mode_a, prec_a, mode_b, prec_b):
    """Inverse-variance fusion on one axis. A precision of 0 contributes nothing and cannot move a mode."""
    p = prec_a + prec_b
    if p <= _EPS:
        return mode_a, 0.0
    return (prec_a * mode_a + prec_b * mode_b) / p, p


def eta_node_sweep(
    chain, statics, geometry, belief, region_arrays, *,
    rna_sense_frac, gdna_strand_overdispersion=0.0, rna_strand_overdispersion=0.0,
    n_gdna_obs=0.0, n_rna_obs=0.0, n_grid, logodds_window=10.0, n_tilt=None, n_grid_ss=None,
    gdna_prior=None, intron_prior=None, length_loglik=None, _capture=None,
):
    """A drop-in for `bp_solver.node_sweep`: same signature, same return, different messages."""
    FIRED[0] += 1
    left = np.asarray(chain.left, dtype=np.int64)
    right = np.asarray(chain.right, dtype=np.int64)
    n = int(chain.n_slots)

    # ── 1. THE GEOMETRY AND THE POPULATION SET — both READS, neither a model ─────────────────────────
    mass_global, eff_global = node_global_geometry(geometry)
    M = np.asarray(mass_global, np.float64)
    E_g = np.asarray(eff_global, np.float64)
    E_r = np.asarray(geometry.eff_rna, np.float64)
    CNT = np.asarray(geometry.unspliced_count, np.float64)
    n_slot = CNT.sum(axis=1)
    u_pos, u_neg = CNT[:, 0], CNT[:, 1]
    spliced_slot = np.asarray(geometry.spliced_count, np.float64).sum(axis=1)
    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)

    T = population_set(fp, fn, E_g, E_r)
    geo = log_opportunity_ratio(E_g, E_r)  # THE conversion: lambda = eta + geo
    # ⭐ |T| == 1 with gDNA the only member IS structural composition certainty — derived from two bits
    # and two opportunities, not carried as a separate `struct_lock` flag.
    pure_gdna = (T.sum(axis=1) == 1) & T[:, GDNA]
    solvable = (fp | fn) & (n_slot > 0.0)
    # ⛔⛔ AN EDGE'S |T| = 1 IS NOT A NODE'S |T| = 1, AND §9 IS WRONG TO COLLAPSE THEM. An EDGE's
    # population set is the INTERSECTION of its two endpoints', so ``T(EDGE) = {g}`` says only "nothing
    # but gDNA can cross here CONTIGUOUSLY" — it does not say the mass sitting on the seam is gDNA. At a
    # zero-gDNA library that mass is RNA contamination, and letting the seam emit ``rho_g = M/E_g``
    # declares all of it gDNA: measured **+96,923 %** on the ``g00`` control, unchanged by removing §2,
    # which is how we know this and not §2 is that control's mechanism. It is the phantom-gDNA emitter
    # `strand_evidence`'s docstring names and `TRAPS.md` D4j's cancelling pair.
    # ⭐ So EMISSION is NODE-scoped, and that is derived here rather than imported: a NODE's ``T`` is its
    # own, an EDGE's is an intersection, and only the first licenses a claim about what the mass IS. This
    # reproduces HEAD's node-only ``struct_lock`` — the two predicates differ for a reason.
    is_node_slot = np.asarray(chain.kind) == NODE
    emits_level = pure_gdna & is_node_slot

    _, solve_grid = _logodds_grid(int(n_grid), float(logodds_window))
    global_lp = (
        gdna_prior.logprior(solve_grid, mass_global, eff_global) if gdna_prior is not None else None
    )

    # ── 2. EACH SLOT'S OWN CLAIM ─────────────────────────────────────────────────────────────────────
    ni = build_node_init(
        chain, statics, geometry, kappa=float(rna_sense_frac),
        od_g=gdna_strand_overdispersion, od_r=rna_strand_overdispersion,
        n_gdna_obs=n_gdna_obs, n_rna_obs=n_rna_obs, n_grid=int(n_grid),
        logodds_window=float(logodds_window), n_tilt=n_tilt, n_grid_ss=n_grid_ss,
        belief=belief, global_logprior=global_lp, intron_prior=intron_prior,
        length_loglik=length_loglik,
    )
    f_g_own = np.clip(np.asarray(ni.f_g, np.float64), 1e-12, 1.0 - 1e-12)
    eta_own = (np.log(f_g_own) - np.log1p(-f_g_own)) - geo
    # eta and lambda differ by a CONSTANT, so Var is the same number on both axes: a precision crosses
    # an edge as the identity too, and the conversion has no scale noise to charge.
    prec_own = np.where(pure_gdna, 0.0, np.asarray(ni.tau_lam, np.float64))

    # ── THE LEVEL CLAIM — HEAD's LOCATION, and §2 is deliberately NOT here ───────────────────────────
    # ⛔⛔ §2's Jeffreys LOCATION (``psi0(a+1/2) - log E``) was in this slot and is REFUTED at
    # **+96,299 %** on the zero-gDNA control. Two reasons it comes out, and the second is the general one:
    #
    #  1. It cannot say ZERO. ``exp(psi0(a+1/2))/E > 0`` at every count, so a structurally pure-gDNA slot
    #     over 50 Mb of empty genome asserts a small POSITIVE gDNA density — and at ``g00`` the truth is
    #     exactly 0, so every one of those is a false positive with nothing to cancel it. Measured:
    #     `node_init.rho_g` is an EXACT zero at 60,544 of 70,176 slots, and that is the statement earning
    #     the -98 % at ``g00``. §2 replaces all of them with a floor.
    #  2. It was never needed for what it was credited with. §2's claim that it "dissolves `TRAPS.md`
    #     C0d" is FALSE — C0d is dissolved by §4/§5, which remove the multiplicative transport entirely:
    #     the conversion is an additive constant in log space and the density crossing is the identity,
    #     so nothing divides by a density and zero transports perfectly. The benefit belongs to §4.
    #
    # ⭐ So the level is the slot's OWN observed density, exactly as HEAD forms it — zero at a zero count.
    # A length-aware location remains open as a SEPARATE experiment, and the derived form is the Gamma
    # MODE ``max(a - 1/2, 0)/E`` (exactly 0 at a = 0, length-awareness in the precision) rather than the
    # mean. ⚠ It has a kink at a = 1/2 and Gamma(1/2, E) is badly non-Gaussian, so a (mode, precision)
    # summary needs a Monte-Carlo check before it is trusted.
    rho_g_own = np.asarray(ni.rho_g, np.float64)
    lvl_prec_own = np.where(
        emits_level & (E_g > 0.0) & (not MUTE_LEVEL[0]), np.asarray(ni.prec_g, np.float64), 0.0
    )

    # the tilt: a nuisance parameter, live only where both strands are admissible AND both were counted.
    # ⚠ Also §2-free — a Jeffreys location here would fabricate a tilt out of two zero counts.
    #
    # ⛔⛔⛔ THIS CHANNEL WAS IN THE WRONG COORDINATE AND IT CARRIED 74 % OF THE ``g00`` ERROR.
    # It delivered a raw log-odds ``log(u_pos/u_neg)`` — unbounded, and ±4.6 at a stranded library —
    # into ψ's ``theta`` slot, whose grid is the ANGLE ``arcsin(tau)`` and spans exactly [−π/2, +π/2]
    # (`simplex_logodds._tilt_grid`). A mode 2.9× outside the entire domain makes ψ's Gaussian monotone
    # across the whole grid, so the tilt pins at the boundary ``tau = ±1`` ("all RNA on one strand") and
    # the strand likelihood — no longer free to integrate the nuisance out — explains the counts it then
    # cannot fit by calling the mass gDNA. Measured at ``g00 ss0.99 capture_on``, ψ-boundary ablation:
    # Σ|err| 1,751,145 as-is → 452,326 with this one channel muted → 170,142 with every message muted,
    # against 1,677 for HEAD. ⭐ `tilt_angle_from_counts` carries the derivation, and the precision comes
    # out as exactly the count from two independent routes with no new constant.
    #
    # ⛔⛔ AND THE COORDINATE ALONE IS NOT THE WHOLE FIX — the SOURCE was wrong too, measured. ψ's τ is the
    # tilt of the RNA; a raw COUNT tilt is the tilt of the MASS, and F2 relates them by
    # ``τ_obs = (2κ − 1)(1 − f_g)·τ``, which on an antisense protocol INVERTS THE SIGN (κ = 0.0101 at
    # ss 0.99 ⇒ 2κ − 1 = −0.98). Correcting only the coordinate, still on counts, left ``g00`` at
    # 1,568,412 — no better than the 1,516,976 it started at — because it then asserted the right angle
    # with the WRONG SIGN at a precision of n. ⭐ ``ni.rho_pos``/``ni.rho_neg`` are ``f_pos·M/E_r`` and
    # ``f_neg·M/E_r``, i.e. ψ's OWN coordinates up to a common factor that cancels in the ratio, so the
    # tilt needs neither κ nor a belief about ``f_g``. `tilt_angle` carries both derivations.
    both = T[:, 1] & T[:, 2]
    _theta, _prec_theta = tilt_angle(ni.rho_pos, ni.rho_neg, ni.prec_pos, ni.prec_neg)
    theta_own = np.where(both, _theta, 0.0)
    prec_theta_own = np.where(both, _prec_theta, 0.0)

    # ── THE MEASUREMENT CHANNELS — certified RNA, and a pure-gDNA neighbour's count ──────────────────
    # ⛔ These are NOT deleted by the derivation and must not be dropped with the operators. §9 retires
    # ``graft``/``peel`` as SEPARATE GATING CONCEPTS — the routing question they answered is now the T
    # comparison — but the channels themselves are measurements:
    #   * a SPLICED fragment is certified RNA (gDNA cannot splice), so a junction's mature flux is a
    #     direct measurement of the destination exon's RNA density, needing no deconvolution;
    #   * a structurally pure-gDNA neighbour's count is a direct measurement of the gDNA density.
    # Measured while building this: with both channels absent the per-object error is 2.3x base even
    # though the LIBRARY total is within 4 % — and `bp_solver` records that the RNA measurement factor
    # alone carries 75 % of the posterior precision on the confidently-wrong unstranded exons and is the
    # only thing that lets a zero-gDNA library say "my mass is all RNA".
    SPL = np.asarray(geometry.junction_count, np.float64)
    ESP = np.asarray(geometry.eff_junction, np.float64)
    spl_live = (SPL > 0.0) & (ESP > 0.0)
    spl_rho = np.where(spl_live, SPL / np.where(spl_live, ESP, 1.0), 0.0)
    _, spl_var = log_rate_posterior(SPL, np.maximum(ESP, _EPS))
    spl_prec = np.where(spl_live, 1.0 / spl_var, 0.0)
    spl_rho_l = spl_rho.tolist()
    spl_prec_l = spl_prec.tolist()

    # the two-bit comparison, precomputed once per directed step.
    same_l = np.zeros(n, bool)
    same_r = np.zeros(n, bool)
    vl, vr = left >= 0, right >= 0
    same_l[vl] = np.all(T[left[vl]] == T[vl], axis=1)
    same_r[vr] = np.all(T[right[vr]] == T[vr], axis=1)

    # ── 3. ONE HOP, USED BY BOTH THE SCANS AND THE COMBINE ───────────────────────────────────────────
    # ⭐⭐ `bp_solver`'s `_relay` and `_transport` are DELIBERATE TWINS with a do-not-merge note, because
    # the relay is a sequential scalar scan and the combine is vectorised, and routing the former through
    # the latter's numpy form measured 15.7x per operation. Under this derivation the hop is THREE LINES
    # of scalar arithmetic — a two-bit comparison and, at a population step, one division — so both sides
    # can call the SAME function and the twin-maintenance hazard (`TRAPS.md` B14) disappears with the
    # operators that created it.
    Ml, Egl, Erl = M.tolist(), E_g.tolist(), E_r.tolist()
    both_l = both.tolist()

    # ⛔⛔ THE LEVEL IS THE SOURCE'S OWN MEASUREMENT AND CROSSES EXACTLY ONE STEP — IT IS NOT A RUNNING
    # CHAIN STATE. §5 licenses "the shared gDNA density crosses unchanged" between two ADJACENT objects;
    # fusing it along the chain is a different claim and it is false under capture. Measured while
    # building this: a chain-fused level is dominated by the 1,312 intergenic anchors (huge E_g, hence
    # huge precision), so every exon inherits the OFF-PROBE floor — and the census puts the true gDNA
    # density at an empty anchor's neighbours at **346x** what the anchor itself reports. The shipped mass
    # pin's case (ii) hands an exon its flanking EDGE's OWN enriched measurement for exactly this reason;
    # a one-step level keeps that property without the pin.
    lvl_src = rho_g_own.tolist()  # a DENSITY now, not a log-density
    pl_src = lvl_prec_own.tolist()

    def _hop(s, i, eta_s, pe_s, th_s, pt_s, same):
        """Deliver the state at ``s`` into ``i``'s frame. Returns ``(eta, prec, theta, prec)``.

        At a same-population step the coordinate is delivered VERBATIM — no reframe, no ratio of totals,
        no scale noise to charge. At a population step the shared gDNA density crosses unchanged and the
        destination's OWN observed count closes the residual.
        """
        if same:
            return eta_s, pe_s, (th_s if both_l[i] else 0.0), (pt_s if both_l[i] else 0.0)
        if pl_src[s] > 0.0 and Ml[i] > _EPS and Erl[i] > 0.0:
            return residual_eta_scalar(lvl_src[s], Ml[i], Egl[i], Erl[i]), pl_src[s], 0.0, 0.0
        return 0.0, 0.0, 0.0, 0.0

    def _scan(nbr, same):
        """One directional pass. Index ``i`` holds ``own(i)`` fused with everything UPSTREAM of ``i``."""
        el, pel = eta_own.tolist(), prec_own.tolist()
        tl, ptl = theta_own.tolist(), prec_theta_own.tolist()
        nb, sm = nbr.tolist(), same.tolist()
        for i in range(n):
            s = nb[i]
            if s < 0:
                continue
            de, dp, dt, dtp = _hop(s, i, el[s], pel[s], tl[s], ptl[s], sm[i])
            el[i], pel[i] = _fuse(el[i], pel[i], de, dp)
            tl[i], ptl[i] = _fuse(tl[i], ptl[i], dt, dtp)
        return el, pel, tl, ptl

    fwd = _scan(left, same_l)
    bwd = _scan(right, same_r)

    # ── 4. THE COMBINE ───────────────────────────────────────────────────────────────────────────────
    # ⛔ THE MESSAGE INTO SLOT i IS THE RELAY STATE AT ITS NEIGHBOUR, NEVER AT i ITSELF. Reading the state
    # at ``i`` hands ψ the node's own belief back as an external claim — on top of the strand likelihood
    # that produced it — which triple-counts the own evidence. Measured while building this: the library
    # gDNA fraction collapsed 0.72 -> 0.18 against a truth of 0.79 on the worst condition.
    eta_m, pe_tot = np.zeros(n), np.zeros(n)
    th_m, pt_tot = np.zeros(n), np.zeros(n)
    # the measurement streams: a gDNA density and a certified-RNA density per strand, each an ABSOLUTE
    # claim on the destination's own share, delivered as log f_c rather than as a composition vote.
    rho_mg, prec_mg = np.zeros(n), np.zeros(n)
    rho_mp, prec_mp = np.zeros(n), np.zeros(n)
    rho_mn, prec_mn = np.zeros(n), np.zeros(n)
    li, ri = left.tolist(), right.tolist()
    sml, smr = same_l.tolist(), same_r.tolist()
    is_exon = np.asarray(_exon_mask(chain, region_arrays)).tolist()
    for i in range(n):
        e, p, t, q = 0.0, 0.0, 0.0, 0.0
        mg, pmg, mp, pmp, mn, pmn = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        for s, state, sm in ((li[i], fwd, sml[i]), (ri[i], bwd, smr[i])):
            if s < 0:
                continue
            de, dp, dt, dq = _hop(s, i, state[0][s], state[1][s], state[2][s], state[3][s], sm)
            e, p = _fuse(e, p, de, dp)
            t, q = _fuse(t, q, dt, dq)
            # a pure-gDNA neighbour's own count, as a gDNA DENSITY measurement (frame-free, one step)
            if pl_src[s] > 0.0:
                mg, pmg = _fuse(mg, pmg, lvl_src[s], pl_src[s])
            # certified RNA: only an EXON receives a junction's mature flux, and only from an EDGE
            if is_exon[i] and spl_prec_l[s][0] + spl_prec_l[s][1] > 0.0:
                mp, pmp = _fuse(mp, pmp, spl_rho_l[s][0], spl_prec_l[s][0])
                mn, pmn = _fuse(mn, pmn, spl_rho_l[s][1], spl_prec_l[s][1])
        eta_m[i], pe_tot[i], th_m[i], pt_tot[i] = e, p, t, q
        rho_mg[i], prec_mg[i] = mg, pmg
        rho_mp[i], prec_mp[i] = mp, pmp
        rho_mn[i], prec_mn[i] = mn, pmn

    # A density measurement becomes a claim on this slot's own SHARE using its own (M, E) — no belief.
    # ⛔⛔ AND AN OVER-UNIT SHARE NOW RAISES RATHER THAN BEING CLIPPED (owner, 2026-08-06). This channel
    # bypasses the consistency rule the composition channel obeys, and delivered a "share" far above 1:
    # measured on the g00 census, the median delivered ``log f_g`` at the 6,128 receiving slots was
    # **+22.41**, i.e. an implied f_g of 5.4e9, carried with real precision. ⛔ Clipping that to 1.0 is
    # MAXIMAL gDNA — the worst answer available at a zero-gDNA library — so the clip converted a loud
    # absurdity into a silent, confident error (`TRAPS.md` A6). A prototype must be brittle: the
    # assertion lives in `_eta_reference.share_from_density`, ONE home, gated by `test_eta_transfer`.
    # ⚠ the strand admissibility mask FIRST: the assertion is about the claims actually DELIVERED, and a
    # forbidden strand's mode is never read.
    prec_mp = np.where(np.asarray(fp, bool), prec_mp, 0.0)
    prec_mn = np.where(np.asarray(fn, bool), prec_mn, 0.0)
    # ⛔⛔ AND A SHARE OF ZERO MASS IS UNDEFINED — the assertion below is what made this visible, on its
    # first run, at pass-0: a junction's certified-RNA density of 45.9 was being delivered as a "share"
    # of a slot holding M = 0, i.e. 5.9e304, and the clip turned it into a confident ``f_g = 0``. There
    # is no composition question at a slot with no unspliced mass to split, so the channel has nothing
    # to claim there. A STRUCTURAL gate on ``M > 0``, the same shape as ``_hop``'s own — not a threshold.
    # ⚠ Expected INERT on the deliverable and not credited with anything: ``solvable`` already excludes
    # these slots, so no belief was written back from them either way (`TRAPS.md` A14 — say what an arm
    # could have moved). What it removes is an absurdity travelling with real precision.
    has_mass = M > 0.0
    prec_mg = np.where(has_mass, prec_mg, 0.0)
    prec_mp = np.where(has_mass, prec_mp, 0.0)
    prec_mn = np.where(has_mass, prec_mn, 0.0)

    # ⛔⛔⛔ AND AN OVER-UNIT SHARE IS THE PREMISE BEING CONTRADICTED, SO THE CLAIM IS REFUSED.
    # The assertion found this on its first run, and it is not a boundary case: a junction's certified
    # RNA density of 10.3 frag/base applied to an exon's 8,300 b of UNSPLICED RNA opportunity claims
    # 85,477 fragments at a slot holding 5. The two are different opportunities for different
    # populations, and the claim rests on the graft's premise — ``rho_R(exon) >= the flux crossing its
    # seam`` — used as an equality. `bp_solver` charges that premise TWO variances (`graft_premise_logvar`,
    # `graft_frame_logvar`) and had the mass pin behind them; this rebuild deleted all three and left the
    # claim unbounded.
    # ⭐ When the implied count exceeds what the destination actually holds, the premise is FALSE here, so
    # the message carries no information about the split. Refusing it adds doubt in NEITHER direction —
    # and that is the whole content of §6's graveyard: eleven mechanisms died because they resolved doubt
    # toward gDNA, and saturating this channel at ``share = 1`` is the same rule pointed at RNA.
    # ⚠ Counted, never silent (`TRAPS.md` A14) — a refusal rate is a measurement of how often the graft's
    # premise fails, which nothing else in this tool reports.
    def _refuse_over_unit(rho, eff, prec, key):
        over = (prec > 0.0) & (np.asarray(rho) * np.asarray(eff) > M)
        REFUSED[key] = REFUSED.get(key, 0) + int(over.sum())
        REFUSED[key + "_mass"] = REFUSED.get(key + "_mass", 0.0) + float(M[over].sum())
        return np.where(over, 0.0, prec)

    prec_mg = _refuse_over_unit(rho_mg, E_g, prec_mg, "gdna")
    prec_mp = _refuse_over_unit(rho_mp, E_r, prec_mp, "rna_pos")
    prec_mn = _refuse_over_unit(rho_mn, E_r, prec_mn, "rna_neg")
    mo_g = share_from_density(rho_mg, E_g, M, prec_mg, name="gdna measurement")
    mo_p = share_from_density(rho_mp, E_r, M, prec_mp, name="certified-RNA + measurement")
    mo_n = share_from_density(rho_mn, E_r, M, prec_mn, name="certified-RNA − measurement")
    # convert into each slot's OWN share frame with its OWN known constant — no belief, no ratio.
    # ⚠ clamp to ψ's OWN grid window. `residual_eta_scalar` saturates at 40 nats when the shared claim
    # exceeds the destination's count, and ρ_g = 0 sends it the other way; both are honest statements of
    # "at or beyond the edge", and the grid is the domain on which ψ can represent them.
    lam_msg = np.clip(eta_m + geo, -float(logodds_window), float(logodds_window))
    ok = solvable & np.isfinite(lam_msg) & (pe_tot > _EPS)
    lam_prec = np.where(ok, pe_tot, 0.0)
    ok_t = solvable & np.isfinite(th_m) & (pt_tot > _EPS)

    dc = _solve_nodes_logodds_all(
        u_pos, u_neg, fp, fn, n_slot, spliced_slot,
        kappa=float(rna_sense_frac), od_g=gdna_strand_overdispersion,
        od_r=rna_strand_overdispersion, n_grid=int(n_grid), L=float(logodds_window),
        n_tilt=n_tilt, n_grid_ss=n_grid_ss, global_logprior=global_lp,
        gdna_imp_mode=mo_g, gdna_imp_prec=prec_mg,
        rna_imp_mode=(mo_p, mo_n), rna_imp_prec=(prec_mp, prec_mn),
        lam_imp_mode=np.where(ok, lam_msg, 0.0), lam_imp_prec=lam_prec,
        theta_imp_mode=np.where(ok_t, th_m, 0.0), theta_imp_prec=np.where(ok_t, pt_tot, 0.0),
        lam_logprior=intron_prior, length_loglik=length_loglik,
        fg_ref=np.asarray(belief.f_g, np.float64),
        fpos_ref=np.asarray(belief.f_pos, np.float64),
        fneg_ref=np.asarray(belief.f_neg, np.float64),
    )

    f_g = np.asarray(belief.f_g, np.float64).copy()
    f_pos = np.asarray(belief.f_pos, np.float64).copy()
    f_neg = np.asarray(belief.f_neg, np.float64).copy()
    var_g = np.asarray(belief.var_gdna, np.float64).copy()
    var_pos = np.asarray(belief.var_pos, np.float64).copy()
    var_neg = np.asarray(belief.var_neg, np.float64).copy()
    f_g = np.where(solvable, np.clip(dc.gdna_frac, 0.0, 1.0), f_g)
    f_pos = np.where(solvable, np.clip(dc.rna_pos_frac, 0.0, 1.0), f_pos)
    f_neg = np.where(solvable, np.clip(dc.rna_neg_frac, 0.0, 1.0), f_neg)
    var_g = np.where(solvable, dc.gdna_frac_var, var_g)
    var_pos = np.where(solvable, dc.rna_pos_frac_var, var_pos)
    var_neg = np.where(solvable, dc.rna_neg_frac_var, var_neg)

    if _capture is not None:  # the keys the scoring path reads, plus the rebuild's own state
        # the strand-ONLY belief (no prior, no messages) — the first rung of `solvability_audit`'s
        # ablation ladder. A debug-only solve, exactly as the shipped sweep computes it, so the two arms'
        # ladders are the same measurement.
        fg_strand = _solve_nodes_logodds_all(
            u_pos, u_neg, fp, fn, n_slot, spliced_slot,
            kappa=float(rna_sense_frac), od_g=gdna_strand_overdispersion,
            od_r=rna_strand_overdispersion, n_grid=int(n_grid), L=float(logodds_window),
            n_tilt=n_tilt, n_grid_ss=n_grid_ss, global_logprior=None,
        ).gdna_frac
        _capture.update(
            _tau0_lam=np.asarray(ni.tau_lam, np.float64),
            free_pos=fp, free_neg=fn, var_g=var_g.copy(), intron_prior=intron_prior,
            fg_strand=fg_strand, f_g=f_g.copy(),
            eta=eta_m, prec_eta=pe_tot, lam_msg=lam_msg, geo=geo, n_pop=T.sum(axis=1),
            # ⛔ THE TILT CHANNEL, PUBLISHED. It was the only message stream the dissection could not
            # see, and it turned out to carry 74 % of the ``g00`` error. A channel absent from the
            # capture is a channel no instrument can rank (`TRAPS.md` A14's shape).
            theta_msg=th_m, prec_theta=pt_tot, theta_own=theta_own, prec_theta_own=prec_theta_own,
            pure_gdna=pure_gdna, same_left=same_l, same_right=same_r, solvable=solvable,
            prec_own_lam=prec_own, lvl_prec_own=lvl_prec_own, lvl_own=rho_g_own,
            prec_mg=prec_mg, mo_g=mo_g, rho_mg=rho_mg, prec_mp=prec_mp, prec_mn=prec_mn,
            ni_rho_g=np.asarray(ni.rho_g, np.float64), ni_prec_g=np.asarray(ni.prec_g, np.float64),
            fg_loc=np.asarray(ni.f_g, np.float64), mass_global=mass_global,
            eff_global=eff_global, eff_rna=E_r, solve_grid=solve_grid, global_lp=global_lp,
        )

    return NodeBelief(
        f_pos=f_pos, f_neg=f_neg, f_g=f_g,
        var_pos=var_pos, var_neg=var_neg, var_gdna=var_g,
    )
