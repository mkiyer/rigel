"""THE BACKBONE — two directional scans, one combine, one ψ solve, one write-back, five assertions.

       Gate: ``tests/calibration/test_sweep_backbone.py``

Each slot's unspliced fragment mass is deconvolved into a pie ``(f_pos, f_neg, f_g)`` — sense-RNA /
antisense-RNA / gDNA — over the ``N E N E … N`` chain (`region_chain`), by ONE forward pass and ONE backward
pass. The chain is a forest of linear paths, so that is exact belief propagation, not an iteration.

⭐⭐⭐ **THIS FILE KNOWS NOTHING ABOUT CAPTURE, GRAFTS, REFRAMES, PINS OR ENRICHMENT — those words do not
appear in it, and that is the design rather than a tidiness.** Everything about *what a message says* is a
:mod:`~.messages` policy. What is left here is the shape of the solve and the five invariants no policy may
break:

===================================================  ====================================================
the backbone asserts                                 it would have caught
===================================================  ====================================================
``deliver`` sees only the two NEIGHBOUR states       **TRAPS: a-message-from-the-destinations-belief — nine recurrences in nine costumes**
every message mode lies inside its own grid          **TRAPS: off-grid-message-mode** — the tilt bug, 74 % of a zero control's error
every delivered share is in ``[0, 1]``               the over-unit certified-RNA claim
``|T| <= 3``                                         **AXIOM 0**, made executable
the write-back touches only ``solvable`` slots       the basis mismatch that made a gate read max|Δ| = 1.0
===================================================  ====================================================

⭐⭐ **The assertions live HERE, not in the policy, and that is the entire point.** A future policy can be
as wrong as it likes and still cannot commit any of these — each of which has shipped at least once.

Two gates on this slot in the pipeline
--------------------------------------
Both come from the region SIGNATURE and never from the counts:

* the SOLVE gate (``solvable``) — a slot deconvolves its own split iff it admits >= 1 RNA strand and has
  unspliced mass. A slot with no admissible RNA strand (an intergenic region, or a gene-boundary boundary) is a
  LOCKED all-gDNA object; it is not solved and keeps its signature-binary init, because RNA cannot cross a
  gene boundary so its unspliced mass is purely gDNA.
* the EMISSION gate — which MESSAGES a slot sends. That is a policy question and lives there.
"""

from __future__ import annotations

import numpy as np

from .messages import NeighbourState, PsiMessage, StepContext
from .messages.silent import SilentPolicy
from .region_geometry import (
    RegionBelief,
    RegionGeometry,
    RegionStatics,
    g1_locked,
    region_gdna_geometry,
)
from .region_init import build_region_init
from .signature import coarse_type_array
from .simplex_logodds import (
    _log_fg,
    _logodds_grid,
    _solve_regions_logodds_all,
    _tilt_grid,
)
from .region_chain import BOUNDARY, REGION, RegionChain, RegionDeconv

__all__ = ["AssertionCounts", "chain_boundary_deconv", "chain_region_deconv", "solve_chain"]


# ⛔⛔ ASSERTIONS THE SHIPPED POLICY IS KNOWN TO VIOLATE, with the measurement beside each. They are
# COUNTED and PUBLISHED rather than raised, because this commit is a restructure gated on byte-identity and
# widening an assertion to fit a defect is how a gate becomes vacuous (TRAPS: perturb-every-gate/TRAPS: a-gate-that-reconstructs). Each carries a
# STRICT xfail in the gate file, which is this project's convention for a PROVEN defect whose fix is
# panel-negative on its own — the fix is `ROADMAP.md`'s TRAPS: a-cancelling-defect-pair pair, not a looser predicate here.
#: ``name -> why it is not fatal yet``. Anything NOT in here raises.
_KNOWN_VIOLATIONS: dict[str, str] = {
    # ── the LOG-SHARE grid, low side. All three counts are new: nothing had checked this coordinate. ─────
    "mode_in_grid_gdna": (
        "MEASURED at 15,240 of 50,984 live gDNA modes (29.9 %) on g50 ss0.50 capture_on, pass-0 — so nearly "
        "a THIRD of the level claims the solver delivers sit outside psi's own log-share grid, whose domain "
        "is [log sigma(-L), log sigma(+L)] = [-10.000045, -4.54e-5] and NOT (-inf, 0]. The cause is the "
        "_EPS = 1e-9 floor on the share: log(1e-9) = -20.723, i.e. 10.72 nats below the low end. "
        "⚠ **DO NOT read this as 29.9 % of the error.** TRAPS: off-grid-message-mode's tilt pin was at the WRONG corner; a low-side "
        "share pin is at the corner meaning 'as little of this component as the grid can express', which at "
        "a slot that genuinely holds almost none of it is the RIGHT answer. What is certain is structural "
        "and is TRAPS: off-grid-message-mode's mechanism: the mode has no interior minimum on the grid, so the accumulated precision "
        "buys a CORNER rather than a location. Whether that corner is right is not measured, and the arm "
        "that would measure it is the TRAPS: a-cancelling-defect-pair pair."
    ),
    "mode_in_grid_rna_pos": (
        "MEASURED at 2,311 of 15,446 live modes (15.0 %). Same coordinate and same _EPS floor as "
        "mode_in_grid_gdna — see its entry for why the count is not a share of the error."
    ),
    "mode_in_grid_rna_neg": (
        "MEASURED at 2,180 of 15,629 live modes (14.0 %). Same coordinate and same _EPS floor as "
        "mode_in_grid_gdna, and the antisense arm is checked separately rather than folded into the sense "
        "one because a channel absent from a check is a channel no instrument can rank (TRAPS: off-grid-message-mode's second lesson)."
    ),
    # ── the mass-conservation identity ───────────────────────────────────────────────────────────────────
    "share_sum_at_most_one": (
        "⛔⛔ MEASURED at 31,174 of 50,984 live packets (61.1 %) on g50 ss0.50 capture_on, pass-0 — a "
        "MAJORITY of delivered packets assert that the three components TOGETHER account for more fragments "
        "than the slot observed. This is the identity Sigma_c rho_c E_c = M that the mass pin exists to "
        "restore, and the pin is licensed in only two states, so everywhere else the residual is delivered "
        "rather than fixed. ⭐ It is consistent with a number already in the tree — messages/variance.py "
        "records the over-claim on 52-71 % of regions — but nothing had surfaced it as a checkable invariant, "
        "so nothing could rank it. ⚠ It is a SUPERSET of share_in_unit_interval by construction and the two "
        "are not independent evidence; report both, and never add them."
    ),
    "mode_in_grid_lam": (
        "MEASURED at 180 of 38,993 live lambda modes (0.46 %) on g50 ss0.50 capture_on, pass-0. The mode is "
        "log(rho_g E_g) - log(rho_R E_r) with both arms floored at _EPS, so a message carrying ONE component "
        "reads +-log(1/_EPS) ~ +-20.7 against a grid half-width of 10. The emission gate zeroes its "
        "PRECISION, which is what keeps it harmless — but a mode that would pin the grid if it were ever "
        "carried is not in-domain, and the count is the thing that says how close that is to happening."
    ),
    "share_in_unit_interval": (
        "MEASURED at 3,013 of 15,629 live share claims (19.3 %) on g50 ss0.50 capture_on, pass-0 — so this "
        "is not a tail, it is a fifth of every level claim the solver delivers. The shares are rho_c E_c / M "
        "with no upper bound, and the certified-RNA arm is a LOWER BOUND used as an equality, which is "
        "exactly what the one-sided term corrects. Landing that term ALONE is refused by the panel: TRAPS: a-cancelling-defect-pair, it "
        "and the missing gDNA level channel are a CANCELLING pair and must be priced in ONE arm. "
        "⭐ And the peel is currently SUPPRESSING it — switching the peel off takes the rate to 41.7 %."
    ),
}


class AssertionCounts(dict):
    """How many slots violated each backbone assertion, published into the diagnostics capture.

    ⭐ A count rather than a bool, because TRAPS: could-the-arm-have-fired is the rule: *before believing "the arm changed
    nothing", check it COULD have changed something.* An assertion reporting 0 violations on a substrate
    where the predicate can never fire is not evidence, so the report also carries how many slots were
    ELIGIBLE for each check.
    """

    def note(self, name: str, violated, eligible) -> None:
        n_v = int(np.count_nonzero(violated))
        n_e = int(np.count_nonzero(eligible)) if eligible is not None else int(np.size(violated))
        self[name] = {"violations": n_v, "eligible": n_e}
        if n_v and name not in _KNOWN_VIOLATIONS:
            raise AssertionError(
                f"backbone assertion {name!r} violated at {n_v:,} of {n_e:,} eligible slots. "
                f"The policy may not do this — see rigel.calibration.sweep's module docstring for what "
                f"each assertion catches. If this is a NEW and PROVEN defect whose fix is panel-negative "
                f"alone, add it to _KNOWN_VIOLATIONS with the measurement, never widen the predicate."
            )


def _check_message(msg: PsiMessage, ctx: StepContext, counts: AssertionCounts) -> None:
    """Assertions 2, 3 and 4 — on what the policy actually delivered.

    ⛔ Assertion 1 (TRAPS: a-message-from-the-destinations-belief) is not checked here because it is enforced BY CONSTRUCTION: ``deliver`` is handed
    :class:`~.messages.NeighbourState`, whose relayed arrays are already gathered at the source, so no
    policy can read a neighbour's relayed belief at the destination. A structural impossibility beats a
    check. ⛔ Assertion 5 is checked at the write-back, where its basis lives.
    """
    # ── (4) AXIOM 0, made executable: |T(slot)| = 1 + free_pos + free_neg, and it is <= 3 ALWAYS ────────
    # There are THREE populations and there is no fourth. This is a function of TWO BITS, which is what
    # makes it structural rather than something to remember — and the message packet carries exactly three
    # component channels (gDNA, RNA+, RNA-) for the same reason.
    pop = ctx.population_size()
    counts.note("population_at_most_three", pop > 3, np.ones_like(pop, bool))
    n_rna_channels = 0 if msg.rna_mode is None else len(msg.rna_mode)
    counts.note(
        "message_has_three_components",
        np.array([n_rna_channels not in (0, 2)]),
        np.array([True]),
    )
    # ⛔ TRAPS: could-the-arm-have-fired's ANTI-DEGENERACY CLAUSE, and it is the half that makes the gate mean anything: on a chain
    # where NO slot admits both RNA strands, ``|T| <= 3`` is satisfied by a substrate that never had a
    # three-population slot to test. That is not the axiom holding, it is the check never running — so the
    # eligible set is the slots that actually reach 3, and a substrate with none of them says so.
    counts.note("population_reaches_three", np.array([]), pop >= 3)

    # ── (2) TRAPS: off-grid-message-mode: every mode inside its coordinate's own grid domain ─────────────────────────────────────
    # ⛔ An out-of-range mode is the MOST confident statement a channel can make, not the least: the penalty
    # -1/2 p (x - m)^2 with m off-grid is monotone across every grid point, so it has no interior minimum
    # and precision buys a CORNER rather than a location. That is how a raw log-odds delivered into the tilt
    # slot — 2.9x outside the whole domain — pinned the tilt at ``tau = +-1``, destroyed the AMBIG Schur
    # protection that keeps the strand term out of ``f_g``, and carried 74 % of a zero control's error.
    # Only a mode at POSITIVE precision can pin anything, so that is the eligible set.
    #
    # ⭐⭐ THE TOLERANCE IS THE GRID'S OWN SPACING, AND IT IS DERIVED RATHER THAN TUNED (TRAPS: no-magic-numbers).
    # Within one spacing of the outermost grid point the penalty's minimum is still that boundary CELL,
    # which for the tilt is the legitimate answer "all RNA on one strand" — so a mode there is correct and
    # not a pin. Beyond it the mode asserts something the coordinate cannot represent.
    # ⚠ **And the spacing is load-bearing rather than cosmetic: measured, the shipped tilt mode overshoots
    # ``pi/2`` by 2 ULP (4.44e-16) at 63 of 4,795 live slots**, because it is a convex mean
    # ``(p_a th_a + p_b th_b)/(p_a + p_b)`` of two messages that AGREE on ``tau = +-1`` and the division
    # rounds up. A bare ``m > hi`` reports that as a defect. TRAPS: off-grid-message-mode's real magnitude is **57 spacings**, so
    # this predicate separates them by more than four orders of magnitude and needs no threshold.
    # ⛔⛔ **AND THE LOW SIDE OF THE LOG-SHARE GRID IS THE HALF THAT WAS MISSING.** The gDNA and RNA modes
    # are ``log`` of a share, so their coordinate is psi's own ``log f`` lattice — whose domain is
    # ``[log sigma(-L), log sigma(+L)]`` = **[-10.000045, -4.54e-5]** at L = 10, not ``(-inf, 0]``. The
    # combine floors every share at ``_EPS = 1e-9``, and ``log(1e-9) = -20.723`` is **10.72 nats BELOW the
    # grid's low end** — the same off-grid pin as TRAPS: off-grid-message-mode, pointing the other way, and the first version of this
    # check tested only ``m > 0`` and could not see it. ⭐ The lambda-emission gate exists precisely because
    # of this floor on the lambda channel (its own comment: 90 slots solved to ``f_g = 1.000`` against a mean
    # oracle of 0.814); **no equivalent gate exists for the gDNA or RNA channel**, which is why the count
    # here is worth having.
    lam_grid, _ = _logodds_grid(int(ctx.n_grid), float(ctx.logodds_window))
    share_grid = _log_fg(lam_grid)  # the SHIPPED builder — one home for the domain
    tilt = _tilt_grid(int(ctx.n_grid) if ctx.n_grid else 2)
    for name, mode, prec, grid in (
        ("lam", msg.lam_mode, msg.lam_prec, lam_grid),
        ("theta", msg.theta_mode, msg.theta_prec, tilt),
        ("gdna", msg.gdna_mode, msg.gdna_prec, share_grid),
        (
            "rna_pos",
            None if msg.rna_mode is None else msg.rna_mode[0],
            None if msg.rna_prec is None else msg.rna_prec[0],
            share_grid,
        ),
        (
            "rna_neg",
            None if msg.rna_mode is None else msg.rna_mode[1],
            None if msg.rna_prec is None else msg.rna_prec[1],
            share_grid,
        ),
    ):
        if mode is None or prec is None:
            continue
        m, p = np.asarray(mode, np.float64), np.asarray(prec, np.float64)
        live = p > 0.0
        lo, hi = float(grid[0]), float(grid[-1])
        step = (hi - lo) / max(int(grid.shape[0]) - 1, 1)
        counts.note(f"mode_in_grid_{name}", live & ((m < lo - step) | (m > hi + step)), live)

    # ── (3) every delivered SHARE in [0, 1] ────────────────────────────────────────────────────────────
    # The gDNA and RNA channels' modes are ``log`` of a share of the slot's own observed mass, so a mode
    # above 0 is a claim that a component alone accounts for MORE fragments than the slot holds. That is
    # not a large claim, it is an impossible one — and it is how an implied 85,477 fragments were once
    # delivered to a slot holding 5.
    for name, mode, prec in (
        ("gdna", msg.gdna_mode, msg.gdna_prec),
        (
            "rna_pos",
            None if msg.rna_mode is None else msg.rna_mode[0],
            None if msg.rna_prec is None else msg.rna_prec[0],
        ),
        (
            "rna_neg",
            None if msg.rna_mode is None else msg.rna_mode[1],
            None if msg.rna_prec is None else msg.rna_prec[1],
        ),
    ):
        if mode is None or prec is None:
            continue
        m, p = np.asarray(mode, np.float64), np.asarray(prec, np.float64)
        live = p > 0.0
        counts.note("share_in_unit_interval", live & (m > 0.0), live)

    # ⭐ AND THE SUM, which is the identity the mass pin exists to restore: the three components account for
    # the fragments the slot observed, so ``Sigma_c exp(mode_c) <= 1``. A per-component check passes while
    # three shares of 0.4 each assert 120 % of the slot's mass between them, and that is a different defect
    # from any one of them being over unity.
    if msg.gdna_mode is not None and msg.rna_mode is not None and msg.gdna_prec is not None:
        live_any = np.asarray(msg.gdna_prec, np.float64) > 0.0
        total = np.exp(np.asarray(msg.gdna_mode, np.float64))
        for k, mo in enumerate(msg.rna_mode):
            pk = np.asarray(msg.rna_prec[k], np.float64) if msg.rna_prec is not None else None
            if pk is None:
                continue
            live_any = live_any | (pk > 0.0)
            total = total + np.exp(np.asarray(mo, np.float64))
        counts.note("share_sum_at_most_one", live_any & (total > 1.0), live_any)


def _at_source(state: tuple[np.ndarray, ...], src: np.ndarray) -> tuple[np.ndarray, ...]:
    """Gather every relayed array at the SOURCE slot. ⭐⭐ **THIS IS ASSERTION 1.**

    ``TRAPS.md`` TRAPS: a-message-from-the-destinations-belief has recurred NINE times in nine costumes, and every one of them was a message built
    from the destination's own relayed or fused belief. After this gather the policy holds values FOR THE
    SOURCE and has no way to ask the same array about the destination, so none of the nine is expressible.
    A gather is exact, so this costs nothing in bits — the shipped policy did the same gather one level
    down.
    """
    return tuple(np.asarray(a)[src] for a in state)


def solve_chain(
    chain: RegionChain,
    statics: RegionStatics,
    geometry: RegionGeometry,
    belief: RegionBelief,
    region_arrays,
    *,
    rna_sense_frac: float,
    gdna_strand_overdispersion: float = 0.0,
    rna_strand_overdispersion: float = 0.0,
    n_gdna_obs: float = 0.0,
    n_rna_obs: float = 0.0,
    n_grid: int,
    logodds_window: float = 10.0,
    n_tilt: int | None = None,
    n_grid_ss: int | None = None,
    gdna_prior=None,
    intron_prior=None,
    policy=None,
    _capture: dict | None = None,
) -> RegionBelief:
    """One forward-backward sweep over the chain. Returns the resolved :class:`RegionBelief`.

    ``policy`` is the message-composition policy (:mod:`~.messages`). ⭐ **It defaults to
    :class:`~.messages.silent.SilentPolicy`, which sends nothing** — so a reader of this file plus five
    boundaries holds the whole working system. The shipped answer is
    :class:`~.messages.head.HeadPolicy`, which ``calibrate`` passes explicitly.

    ``gdna_prior=None`` is a first-class PRIOR-FREE solve: ψ then carries the derived reference alone on
    both arms. Prior-free is not reference-free. ⭐ That pass's only job is to be a training substrate for
    the population gDNA hyperprior — it is not the deliverable, and it does not have to answer objects it
    cannot solve.
    """
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    fp, fn = statics.free_pos, statics.free_neg
    f_pos = np.asarray(belief.f_pos, dtype=np.float64).copy()
    f_neg = np.asarray(belief.f_neg, dtype=np.float64).copy()
    f_g = np.asarray(belief.f_g, dtype=np.float64).copy()
    var_pos = np.asarray(belief.var_pos, dtype=np.float64).copy()
    var_neg = np.asarray(belief.var_neg, dtype=np.float64).copy()
    var_g = np.asarray(belief.var_gdna, dtype=np.float64).copy()
    # the INCOMING belief, kept for the diagnostic capture: it is the ``fg_ref`` the final solve freezes
    # its variance at, so a channel-ablation replay must pass the SAME reference to be faithful.
    _fg_init, _fp_init, _fn_init = f_g.copy(), f_pos.copy(), f_neg.copy()

    EG = np.asarray(geometry.eff_gdna, np.float64)
    ER = np.asarray(geometry.eff_rna, np.float64)
    ESP = np.asarray(geometry.eff_junction, np.float64)  # [n, 2] by TRANSCRIPT strand
    SPL = np.asarray(geometry.junction_count, np.float64)  # [n, 2] by TRANSCRIPT strand
    CNT = np.asarray(geometry.unspliced_count, np.float64)  # [n, 2] by GENOME strand
    # the unspliced count is BOTH the density numerator and the Poisson n — one number, not a fractional
    # mass plus a separate integer flux.
    n_slot = CNT.sum(axis=1)
    u_pos, u_neg = CNT[:, 0], CNT[:, 1]
    spliced_slot = np.asarray(geometry.spliced_count, np.float64).sum(axis=1)
    # per-slot "global" gDNA support — the basis the rate prior is fit and projected on.
    mass_global, eff_global = region_gdna_geometry(geometry)

    _, solve_grid = _logodds_grid(int(n_grid), float(logodds_window))
    kappa = float(rna_sense_frac)
    od_g, od_r = gdna_strand_overdispersion, rna_strand_overdispersion

    def _psi(g_arr, msg: PsiMessage, *, fg_ref, fpos_ref, fneg_ref):
        """The per-slot solve (the log-density log-odds backend). Phase A calls it with a silent message;
        the final call passes the combine's four channels.

        ⚠ ``fg_ref`` is the count-zero-information variance freeze: the reference is the incoming belief,
        so the variance — hence the message precision — is evaluated near the truth and not at a flat 1/2.
        It is passed EXPLICITLY rather than closed over, because the write-back rebinds the belief and one
        diagnostic solve below deliberately runs after that.
        """
        return _solve_regions_logodds_all(
            u_pos,
            u_neg,
            fp,
            fn,
            n_slot,
            spliced_slot,
            kappa=kappa,
            od_g=od_g,
            od_r=od_r,
            n_grid=int(n_grid),
            L=float(logodds_window),
            n_tilt=n_tilt,
            n_grid_ss=n_grid_ss,
            global_logprior=g_arr,
            gdna_imp_mode=msg.gdna_mode,
            gdna_imp_prec=msg.gdna_prec,
            rna_imp_mode=msg.rna_mode,
            rna_imp_prec=msg.rna_prec,
            lam_imp_mode=msg.lam_mode,
            lam_imp_prec=msg.lam_prec,
            theta_imp_mode=msg.theta_mode,
            theta_imp_prec=msg.theta_prec,
            # the gDNA intron-factory λ-factor (anchored, per-intron, 0 elsewhere): peels confident gDNA
            # from introns against the intergenic background BEFORE the sweep resolves the pie. Added to ψ,
            # distinct from the gDNA arm; participates in the local solve AND the relay.
            lam_logprior=intron_prior,
            # ⭐ the FRAGMENT-LENGTH λ-factor. It enters the LOCAL solve and the FINAL one, exactly like
            # the intron factory, so a slot's own length evidence both sets its belief and propagates.
            # ``None`` ⇒ byte-identical to the path without it.
            fg_ref=fg_ref,
            fpos_ref=fpos_ref,
            fneg_ref=fneg_ref,
        )

    # ── the SOLVE gate. Structural, from the signature, never from the counts ──────────────────────────
    solvable = (fp | fn) & (n_slot > 0.0)

    # THE gDNA ARM of ψ — the COMPOSITION prior, and ONLY that. A total-density model is an ENRICHMENT
    # model, not a DNA composition prior: letting it vote a slot's f_g is the count-votes-composition
    # regression.
    global_lp = (
        gdna_prior.logprior(solve_grid, mass_global, eff_global) if gdna_prior is not None else None
    )

    # ⭐ Slot ids ARE the genomic visiting order, so the order is ``arange`` and the chain does not store
    # it. The scans are sequential, so iterate as a Python list of ints.
    order_list = list(range(int(chain.n_slots)))
    # per-slot EXON-region flag — the graft's destination class, and a policy input rather than a gate.
    _rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    _ri = np.clip(np.asarray(chain.obj_idx, dtype=np.int64), 0, _rtype.shape[0] - 1)
    is_exon_region = (np.asarray(chain.kind) == REGION) & (_rtype[_ri] == 2)

    # ── (A) the per-slot message-free SELF-SOLVE — the four init sources ──────────────────────────────
    own = build_region_init(
        chain,
        statics,
        geometry,
        kappa=kappa,
        od_g=od_g,
        od_r=od_r,
        n_gdna_obs=n_gdna_obs,
        n_rna_obs=n_rna_obs,
        n_grid=int(n_grid),
        logodds_window=float(logodds_window),
        n_tilt=n_tilt,
        n_grid_ss=n_grid_ss,
        belief=belief,
        global_logprior=global_lp,
        intron_prior=intron_prior,
    )

    ctx = StepContext(
        # observations
        mass=mass_global,
        eff_gdna_global=eff_global,
        eff_rna=ER,
        eff_gdna=EG,
        eff_junction=ESP,
        junction_count=SPL,
        unspliced_count=CNT,
        n_slot=n_slot,
        spliced_slot=spliced_slot,
        # geometry / structure
        left=left,
        right=right,
        is_boundary=np.asarray(chain.kind) != REGION,
        is_exon_region=is_exon_region,
        free_pos=np.asarray(fp, bool),
        free_neg=np.asarray(fn, bool),
        g1_locked=g1_locked(np.asarray(fp, bool), np.asarray(fn, bool)),
        boundary_flags=statics.boundary_flags,
        geometry=geometry,
        order=order_list,
        left_list=left.tolist(),
        right_list=right.tolist(),
        # beliefs — SOURCE-SIDE ONLY (TRAPS: a-message-from-the-destinations-belief)
        own=own,
        belief_fg=f_g,
        # the solve's own scalars
        n_grid=int(n_grid),
        logodds_window=float(logodds_window),
        solve_grid=solve_grid,
        capture=_capture,
    )

    relay = (policy if policy is not None else SilentPolicy()).prepare(ctx)

    # ── (B) the FORWARD scan L→R and (C) the BACKWARD scan R→L ────────────────────────────────────────
    # ⛔ TRAPS: a-comment-quoted-as-a-finding: ONE pass each, accumulating IN PLACE, which on a chain IS the forward half of
    # forward-backward. "In place" is why it cannot be vectorised; it is not an iterative scheme, and a
    # source comment's shorthand for the first crossed into a design doc as if it were the second.
    fwd = _scan(order_list, ctx.left_list, relay, backward=False)
    bwd = _scan(order_list[::-1], ctx.right_list, relay, backward=True)

    # ── (D) the COMBINE — both neighbours into this slot's frame, then ONE batched ψ solve ────────────
    n = int(chain.n_slots)
    sl = np.clip(np.asarray(left, np.int64), 0, n - 1)
    sr = np.clip(np.asarray(right, np.int64), 0, n - 1)
    vl, vr = np.asarray(left, np.int64) >= 0, np.asarray(right, np.int64) >= 0
    msg = relay.deliver(
        NeighbourState(state=_at_source(fwd, sl) if fwd else (), valid=vl, src=sl),
        NeighbourState(state=_at_source(bwd, sr) if bwd else (), valid=vr, src=sr),
    )
    counts = AssertionCounts()
    _check_message(msg, ctx, counts)

    dc_fin = _psi(global_lp, msg, fg_ref=f_g, fpos_ref=f_pos, fneg_ref=f_neg)

    # ── THE WRITE-BACK — only SOLVABLE slots ──────────────────────────────────────────────────────────
    # A locked slot (no admissible RNA strand) or an empty one keeps its signature-binary init. ⛔ The
    # alternative — skip UNIDENTIFIED slots too and defer to the prior — was derived, implemented and
    # EMPIRICALLY REFUTED: it regresses both standalone and with the hyperprior, because the prior resolves
    # an imperfectly-solved slot better than a deferred ``f_g = 1``.
    mg_, mp_, mn_ = dc_fin.gdna_frac, dc_fin.rna_pos_frac, dc_fin.rna_neg_frac
    vg_, vp_, vn_ = dc_fin.gdna_frac_var, dc_fin.rna_pos_frac_var, dc_fin.rna_neg_frac_var
    out_fg = np.where(solvable, np.clip(mg_, 0.0, 1.0), f_g)
    out_fpos = np.where(solvable, np.clip(mp_, 0.0, 1.0), f_pos)
    out_fneg = np.where(solvable, np.clip(mn_, 0.0, 1.0), f_neg)
    out_vg = np.where(solvable, vg_, var_g)
    out_vpos = np.where(solvable, vp_, var_pos)
    out_vneg = np.where(solvable, vn_, var_neg)
    # ── (5) the write-back touched ONLY solvable slots ────────────────────────────────────────────────
    # ⛔ The silent version of this made an TRAPS: byte-identity-gate identity gate read ``max|Δ| = 1.0``: a replay compared the
    # solve's raw output against the shipped belief, and the two differ by exactly this mask. Reproducing a
    # pipeline stage means reproducing its WRITE-BACK.
    untouched = ~np.asarray(solvable, bool)
    counts.note(
        "writeback_only_solvable",
        untouched
        & (
            (out_fg != _fg_init)
            | (out_fpos != _fp_init)
            | (out_fneg != _fn_init)
            | (out_vg != np.asarray(belief.var_gdna, np.float64))
        ),
        untouched,
    )
    f_g, f_pos, f_neg = out_fg, out_fpos, out_fneg
    var_g, var_pos, var_neg = out_vg, out_vpos, out_vneg

    if _capture is not None:  # inert diagnostic hook
        # strand-ONLY local belief (no global prior, no messages) — to split the local error into the
        # strand likelihood vs the global gDNA prior contribution. Same solver, global=None.
        fg_strand = _solve_regions_logodds_all(
            u_pos,
            u_neg,
            fp,
            fn,
            n_slot,
            spliced_slot,
            kappa=kappa,
            od_g=od_g,
            od_r=od_r,
            n_grid=int(n_grid),
            L=float(logodds_window),
            n_tilt=n_tilt,
            n_grid_ss=n_grid_ss,
            global_logprior=None,
        ).gdna_frac
        # the message-free self-solve variances, for the local-error attribution — a debug-only solve, so
        # the production path carries none of it (the self-solve fractions come from ``own``).
        # ⚠ Deliberately AFTER the write-back, so its reference is the OUTGOING belief. That is what the
        # shipped solver did, and an instrument comparing against it depends on the same reference.
        _dc_loc = _psi(global_lp, PsiMessage.silent(), fg_ref=f_g, fpos_ref=f_pos, fneg_ref=f_neg)
        _capture.setdefault("_uni_static", {}).update(
            order=np.asarray(order_list, np.int64),
            n_slot=n_slot.copy(),
            left=np.asarray(left, np.int64),
            right=np.asarray(right, np.int64),
        )
        # ⭐ the RAW per-slot relayed state, both directions. The BACKBONE publishes it because the
        # backbone is the only thing that holds it un-indexed: ``deliver`` receives it gathered at the
        # source, which is assertion 1. The dissect loop reads these keys, so the names are the shipped
        # ones and the values are the arrays exactly as each scan published them.
        _NAMES = ("g", "p", "n", "pg", "pp", "pn", "mg", "mp", "mn", "tau")
        for _tag, _st in (("fwd", fwd), ("bwd", bwd)):
            if _st is not None:
                _capture["_uni_static"].update({f"{_tag}_{k}": v for k, v in zip(_NAMES, _st)})
        _capture.update(
            backbone_assertions=counts,
            fg_loc=own.f_g,
            fg_strand=fg_strand,
            fp_loc=own.f_pos,
            fn_loc=own.f_neg,
            vg_loc=_dc_loc.gdna_frac_var,
            vp_loc=_dc_loc.rna_pos_frac_var,
            f_g=f_g.copy(),
            f_pos=f_pos.copy(),
            f_neg=f_neg.copy(),
            var_g=var_g.copy(),
            solvable=solvable,
            count=CNT,
            spliced=spliced_slot,
            mature=SPL.sum(axis=1),
            free_pos=np.asarray(fp, bool),
            free_neg=np.asarray(fn, bool),
            eff_global=eff_global,
            mass_global=mass_global,
            eff_gdna=EG,
            eff_rna=ER,
            # the full per-slot global prior term on the solve grid, so a diagnostic can replay the solve
            # with message channels ablated (the message help/hurt attribution).
            global_lp=global_lp,
            solve_grid=solve_grid,
            _tau0_lam=own.tau_lam,
            # the incoming belief (the final solve's ``fg_ref``) + the intron-factory λ arm, so an ablation
            # replay reproduces the shipped f_g exactly BEFORE ablating.
            fg_init=_fg_init,
            fpos_init=_fp_init,
            fneg_init=_fn_init,
            intron_prior=intron_prior,
            solvable_mask=solvable,
        )

    return RegionBelief(
        f_pos=f_pos,
        f_neg=f_neg,
        f_g=f_g,
        var_pos=var_pos,
        var_neg=var_neg,
        var_gdna=var_g,
    )


def _scan(seq, nbr, relay, *, backward: bool):  # noqa: D401 — see solve_chain's (B)/(C)
    """ONE directional pass: for each slot in ``seq``, relay from its neighbour of the other kind.

    ⭐ **The whole direction dependence is which neighbour array is read.** ``-1`` is a reference terminal
    and therefore a propagation sink — it is skipped, so a sweep cannot relay across a reference boundary.

    Returns the policy's published state, or ``None`` if the policy relays nothing at all.
    """
    kernel = relay.scan(backward=backward)
    if kernel is None:
        return None
    step, publish = kernel
    for i in seq:
        s = nbr[i]
        if s < 0:
            continue
        step(s, i)
    return publish()


# ──────────────────────────────────────────────────────────────────────────────────────────────────────
# THE OUTPUT CONTRACT — projecting the chain belief back onto the two payload axes.
# Not part of the solve: this is what ``CalibrationResult`` / ``priors`` / ``derive`` consume.
# ──────────────────────────────────────────────────────────────────────────────────────────────────────


def chain_region_deconv(chain: RegionChain, belief: RegionBelief, substrate) -> RegionDeconv:
    """Project the chain belief's REGION slots back onto the REGION axis as a :class:`RegionDeconv` — what
    ``CalibrationResult`` / ``priors`` / ``derive`` consume.

    ⚠ **A region's contained population carries no spliced term any more, and that is structural**: the
    accumulator credits ``region_contained`` only when the fragment used no junction, so a contained
    fragment is unspliced by construction. The predecessor added ``+ mass_spliced`` here; that quantity
    is identically zero on the region axis now, and adding it would be adding a channel that cannot exist.
    """
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    reg = kind == REGION
    count = np.asarray(substrate.region_contained.count, dtype=np.float64).sum(axis=1)
    n = count.shape[0]
    f_g = np.zeros(n)
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    ri = idx[reg]
    f_g[ri] = belief.f_g[reg]
    f_pos[ri] = belief.f_pos[reg]
    f_neg[ri] = belief.f_neg[reg]
    return RegionDeconv(
        gdna_mass=f_g * count,
        rna_mass=(1.0 - f_g) * count,
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
    )


def chain_boundary_deconv(chain: RegionChain, belief: RegionBelief, substrate) -> RegionDeconv:
    """Project the chain belief's BOUNDARY slots onto the CONTIGUOUS-BOUNDARY axis — the crossing flux that
    ``priors`` and ``derive`` consume.

    ⭐ **ONE per-boundary result, not a ``(left, right)`` pair of per-region ones.** The predecessor split
    each boundary's flux onto its two flanking regions and ``priors`` then pooled the two halves straight
    back together — so the split and the re-pool were a no-op, and that exact sum-then-halve pattern is
    what hid a factor of 2 for months. Owner ruling, 2026-07-30:
    ``CalibrationResult``'s per-region ``mass_*_left/right`` become per-boundary arrays.

    The RNA mass is spliced-inclusive: an boundary's certified-RNA crossings (``boundary_spliced``) are RNA
    whatever the unspliced mixture resolves to, since gDNA cannot be spliced.
    """
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    boundary = kind == BOUNDARY
    unspliced = np.asarray(substrate.boundary_unspliced.count, dtype=np.float64).sum(axis=1)
    spliced = np.asarray(substrate.boundary_spliced.count, dtype=np.float64).sum(axis=1)
    n = unspliced.shape[0]
    ei = idx[boundary]
    f_g = np.zeros(n)
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    f_g[ei] = np.asarray(belief.f_g, dtype=np.float64)[boundary]
    # ⭐ THE PER-STRAND RNA SPLIT, PROJECTED ON THIS AXIS TOO. ψ solves the simplex
    # ``(f_g, f_pos, f_neg)`` at EVERY slot — AXIOM 0's `T(slot)`, which is a function of the two
    # `free_*` bits and never of the slot's kind — so an BOUNDARY slot has the same three-way composition a
    # REGION slot does. ⛔ This projection used to emit ``np.zeros(n)`` for both RNA strands, so the
    # crossing axis published a composition that summed to ``f_g`` alone. Nothing consumed it, which is
    # why it survived; a per-transcript prior reading composition per object does.
    f_pos[ei] = np.asarray(belief.f_pos, dtype=np.float64)[boundary]
    f_neg[ei] = np.asarray(belief.f_neg, dtype=np.float64)[boundary]
    return RegionDeconv(
        gdna_mass=f_g * unspliced,
        rna_mass=(1.0 - f_g) * unspliced + spliced,
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
    )
