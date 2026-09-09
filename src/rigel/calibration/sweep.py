"""THE BACKBONE — two directional scans, one combine, one ψ solve, one write-back, five assertions.

       Gate: ``tests/calibration/test_sweep_backbone.py``

Each slot's unspliced fragment mass is deconvolved into a pie ``(f_pos, f_neg, f_g)`` — sense-RNA /
antisense-RNA / gDNA — over the ``N E N E … N`` chain (`region_chain`), by ONE forward pass and ONE backward
pass. The chain is a forest of linear paths, so that is exact belief propagation, not an iteration.

⭐⭐⭐ **THIS FILE KNOWS NOTHING ABOUT CAPTURE, SPLICE IN, REFRAMES, PINS OR ENRICHMENT — those words do not
appear in it, and that is the design rather than a tidiness.** Everything about *what a message says* is a
:mod:`~.messages` policy. What is left here is the shape of the solve and the five invariants no policy may
break:

===================================================  ====================================================
the backbone asserts                                 it would have caught
===================================================  ====================================================
the kernel sees only the two NEIGHBOUR states        **TRAPS: a-message-from-the-destinations-belief — nine recurrences in nine costumes**
every delivered row is one row per slot, finite      a row array off the solve grid, or a NaN reaching ψ
``|T| <= 3``                                         **AXIOM 0**, made executable
the write-back touches only ``solvable`` slots       the basis mismatch that made a gate read max|Δ| = 1.0
===================================================  ====================================================

⭐⭐ **The assertions live HERE, not in the policy, and that is the entire point.** A future policy can be
as wrong as it likes and still cannot commit any of these — each of which has shipped at least once.
(Two more assertions guarded the retired relay's Gaussian channels — a mode inside its coordinate's grid,
**TRAPS: off-grid-message-mode**, and every delivered share in ``[0, 1]``; the channels retired with the
relay on 2026-09-09 and the transfer policy delivers profiles on the solve grid, which cannot commit
either.)

Two gates on this slot in the pipeline
--------------------------------------
Both come from the region SIGNATURE and never from the counts:

* the SOLVE gate (``solvable``) — a slot deconvolves its own split iff it admits >= 1 RNA strand and has
  unspliced mass. A slot with no admissible RNA strand (an intergenic region, or a gene boundary) is a
  LOCKED all-gDNA object; it is not solved and keeps its signature-binary init, because RNA cannot cross a
  gene boundary so its unspliced mass is purely gDNA.
* the EMISSION gate — which MESSAGES a slot sends. That is a policy question and lives there.
"""

from __future__ import annotations

import numpy as np

from .messages import SILENCE, PsiMessage, StepContext
from .messages.silent import SilentPolicy
from .region_geometry import (
    RegionBelief,
    RegionGeometry,
    RegionStatics,
    region_gdna_geometry,
    region_rna_geometry,
)
from .region_init import build_region_init
from .signature import BIT_EXON_NEG, BIT_EXON_POS, coarse_type_array
from .structural_claims import build_structural_claims, interface_masks
from .simplex_logodds import (
    CompositionPriors,
    _logodds_grid,
    _solve_regions_logodds_all,
)
from .region_chain import BOUNDARY, REGION, RegionChain, RegionDeconv

__all__ = ["AssertionCounts", "chain_boundary_deconv", "chain_region_deconv", "solve_chain"]


# ⛔⛔ ASSERTIONS A SHIPPED POLICY IS KNOWN TO VIOLATE, with the measurement beside each. They would be
# COUNTED and PUBLISHED rather than raised, because widening an assertion to fit a defect is how a gate
# becomes vacuous (TRAPS: perturb-every-gate / TRAPS: a-gate-that-reconstructs); each would carry a STRICT
# xfail in the gate file, this project's convention for a PROVEN defect whose fix is panel-negative on its
# own. ⭐ EMPTY since 2026-09-09: every entry guarded the retired relay's Gaussian channels (off-grid share
# modes at the _EPS floor, over-unit shares, the share sum) and retired with them; the transfer policy
# delivers max-normalised profiles on the solve grid. Anything NOT in here raises.
#: ``name -> why it is not fatal yet``.
_KNOWN_VIOLATIONS: dict[str, str] = {}


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
    """The assertions on what the policy actually delivered: the population axiom, and every row
    channel one row per slot on the solve grid and finite.

    ⛔ Assertion 1 (TRAPS: a-message-from-the-destinations-belief) is not checked here because it is enforced BY CONSTRUCTION: the
    propagate kernel is called with two INDICES and builds the message into the destination from the
    source's claim and what the source holds; the backbone writes ``held`` and the policy never reaches
    past its hop. A structural impossibility beats a check. ⛔ Assertion 5 is checked at the write-back,
    where its basis lives.
    """
    # ── (4) AXIOM 0, made executable: |T(slot)| = 1 + free_pos + free_neg, and it is <= 3 ALWAYS ────────
    # There are THREE populations and there is no fourth. This is a function of TWO BITS, which is what
    # makes it structural rather than something to remember — and the message packet carries exactly three
    # component channels (gDNA, RNA+, RNA-) for the same reason.
    pop = ctx.population_size()
    counts.note("population_at_most_three", pop > 3, np.ones_like(pop, bool))
    # ── the λ rows: one row per slot in psi's general evidence currency, or absent ───────────────────
    if msg.lam_rows is not None:
        rows = np.asarray(msg.lam_rows)
        if rows.shape[0] != ctx.n_slots or rows.ndim != 2:
            raise ValueError(
                f"lam_rows has shape {rows.shape}; expected ({ctx.n_slots}, K) — a policy must "
                "deliver one row per slot on the solve grid"
            )
        counts.note("flux_rows_finite", ~np.isfinite(rows).all(axis=1), np.ones(ctx.n_slots, bool))
    # ── the cube channel: a (K, K_t) row per AMBIG slot, or absent ──────────────────────────────────
    if msg.cube_rows is not None:
        amb = np.asarray(ctx.free_pos, bool) & np.asarray(ctx.free_neg, bool)
        for slot, row in msg.cube_rows.items():
            row = np.asarray(row)
            if not (0 <= int(slot) < ctx.n_slots) or not amb[int(slot)]:
                raise ValueError(
                    f"cube_rows carries slot {slot}, which is not an AMBIG slot — a cube exists only "
                    "where both strands are live"
                )
            if row.ndim != 2 or row.shape[0] != int(ctx.n_grid):
                raise ValueError(
                    f"cube_rows[{slot}] has shape {row.shape}; expected ({ctx.n_grid}, K_t) — one row "
                    "over the (λ, θ) cube on the solve grid"
                )
        bad = np.array([not np.isfinite(np.asarray(r)).all() for r in msg.cube_rows.values()], bool)
        counts.note("cube_rows_finite", bad, np.ones(bad.shape[0], bool))
    # ⛔ TRAPS: could-the-arm-have-fired's ANTI-DEGENERACY CLAUSE, and it is the half that makes the gate mean anything: on a chain
    # where NO slot admits both RNA strands, ``|T| <= 3`` is satisfied by a substrate that never had a
    # three-population slot to test. That is not the axiom holding, it is the check never running — so the
    # eligible set is the slots that actually reach 3, and a substrate with none of them says so.
    counts.note("population_reaches_three", np.array([]), pop >= 3)


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
    rna_prior=None,
    intron_prior=None,
    policy=None,
    _capture: dict | None = None,
) -> RegionBelief:
    """One forward-backward sweep over the chain. Returns the resolved :class:`RegionBelief`.

    ``policy`` is the message-composition policy (:mod:`~.messages`). ⭐ **It defaults to
    :class:`~.messages.silent.SilentPolicy`, which sends nothing** — so a reader of this file plus five
    boundaries holds the whole working system. The shipped answer is
    :class:`~.messages.transfer.TransferPolicy`, which ``calibrate`` passes explicitly.

    ``gdna_prior=None`` is a first-class PRIOR-FREE solve: ψ then carries the Jeffreys reference
    measure alone on both arms. Prior-free is not reference-free. ⭐ That pass's only job is to be a training substrate for
    the population gDNA hyperprior — it is not the deliverable, and it does not have to answer objects it
    cannot solve.

    ⛔ ψ carries NO reference location (owner refutation, 2026-08-24): the reference is the
    symmetric Jeffreys measure and asserts nothing; background information enters as the
    ``intron_prior`` λ-factor, a likelihood whose precision scales with counts.
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
    ESP = np.asarray(geometry.eff_sj, np.float64)  # [n, 2] by TRANSCRIPT strand
    SPL = np.asarray(geometry.sj_count, np.float64)  # [n, 2] by TRANSCRIPT strand
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

    def _psi(
        g_arr, msg: PsiMessage, *, fg_ref, fpos_ref, fneg_ref, extra_lam_rows=None, cube_rows=None
    ):
        """The per-slot solve (the log-density log-odds backend). Phase A calls it with a silent message;
        the final call passes the policy's rows (the λ rows and the cube rows).

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
            priors=g_arr,
            # the gDNA intron-factory λ-factor (anchored, per-intron, 0 elsewhere): deconvolves confident gDNA
            # from introns against the intergenic background BEFORE the sweep resolves the pie. Added to ψ,
            # distinct from the gDNA arm; participates in the local solve AND the message layer.
            lam_logprior=(
                intron_prior
                if extra_lam_rows is None
                else (
                    np.asarray(extra_lam_rows, np.float64)
                    if intron_prior is None
                    else intron_prior + np.asarray(extra_lam_rows, np.float64)
                )
            ),
            fg_ref=fg_ref,
            fpos_ref=fpos_ref,
            fneg_ref=fneg_ref,
            # ⭐ the CUBE channel: the RNA level lanes delivered at AMBIG slots, final solve only
            cube_rows=cube_rows,
        )

    # ── the SOLVE gate. Structural, from the signature, never from the counts ──────────────────────────
    solvable = (fp | fn) & (n_slot > 0.0)

    # THE gDNA ARM of ψ — the COMPOSITION prior, and ONLY that. A total-density model is an ENRICHMENT
    # model, not a DNA composition prior: letting it vote a slot's f_g is the count-votes-composition
    # regression.
    # ⭐ ONE construction site for ψ's composition arms. The RNA member stays ``None`` until something
    # fits an RNA landscape; ``None`` there means "that arm takes its derived reference", which is the
    # shipped behaviour and a first-class configuration rather than a gap.
    # ⭐ The RNA arm asks the SAME landscape about the OTHER component: the complementary fraction
    # `1 - f_g` against RNA's own opportunity. `mass_global` is shared because both components split one
    # unspliced population (`region_rna_geometry`).
    _rna_mass, _eff_rna = region_rna_geometry(geometry)
    global_lp = CompositionPriors(
        gdna=gdna_prior.logprior(solve_grid, mass_global, eff_global)
        if gdna_prior is not None
        else None,
        rna=rna_prior.logprior(1.0 - solve_grid, _rna_mass, _eff_rna)
        if rna_prior is not None
        else None,
    )

    # ⭐ Slot ids ARE the genomic visiting order, so the order is ``arange`` and the chain does not store
    # it. The scans are sequential, so iterate as a Python list of ints.
    order_list = list(range(int(chain.n_slots)))
    # per-slot EXON-region flag — the SPLICE IN's destination class, and a policy input rather than a gate.
    _rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    _ri = np.clip(np.asarray(chain.obj_idx, dtype=np.int64), 0, _rtype.shape[0] - 1)
    is_exon_region = (np.asarray(chain.kind) == REGION) & (_rtype[_ri] == 2)
    # the signature's two EXON bits per slot — a strand's opportunity geometry, which the RNA level
    # lanes read to tell a strand's intron from its exon inside an overlapping locus
    _sig = np.asarray(region_arrays.signature).astype(np.int64)[_ri]
    _is_region = np.asarray(chain.kind) == REGION
    exon_pos = _is_region & ((_sig & BIT_EXON_POS) > 0)
    exon_neg = _is_region & ((_sig & BIT_EXON_NEG) > 0)
    # the stage-0 interface masks the message layer consumes — computed by the module that
    # owns the concept and carried here under policy-neutral names (the backbone's vocabulary
    # firewall keeps every message-composition concept out of this file)
    _if_left, _if_right, _ss_b = interface_masks(build_structural_claims(chain, statics))

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
        priors=global_lp,
        intron_prior=intron_prior,
    )

    ctx = StepContext(
        # observations
        mass=mass_global,
        inv_abundance=np.asarray(geometry.inv_abundance, np.float64),
        inv_sj_lo=np.asarray(geometry.inv_sj_lo, np.float64),
        inv_sj_hi=np.asarray(geometry.inv_sj_hi, np.float64),
        eff_gdna_global=eff_global,
        eff_rna=ER,
        eff_gdna=EG,
        eff_sj=ESP,
        sj_count=SPL,
        sj_count_lo=np.asarray(geometry.sj_count_lo, np.float64),
        sj_count_hi=np.asarray(geometry.sj_count_hi, np.float64),
        route_rate_lo=np.asarray(geometry.route_rate_lo, np.float64),
        route_rate_hi=np.asarray(geometry.route_rate_hi, np.float64),
        route_count_lo=np.asarray(geometry.route_count_lo, np.int64),
        route_count_hi=np.asarray(geometry.route_count_hi, np.int64),
        unspliced_count=CNT,
        n_slot=n_slot,
        spliced_slot=spliced_slot,
        # geometry / structure
        left=left,
        right=right,
        is_boundary=np.asarray(chain.kind) != REGION,
        is_exon_region=is_exon_region,
        left_interface_certified=_if_left,
        right_interface_certified=_if_right,
        ss_intron_boundary=_ss_b,
        free_pos=np.asarray(fp, bool),
        free_neg=np.asarray(fn, bool),
        exon_pos=exon_pos,
        exon_neg=exon_neg,
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
        n_tilt=None if n_tilt is None else int(n_tilt),
    )

    prepared = (policy if policy is not None else SilentPolicy()).prepare(ctx)

    # ── PHASE 1, PROPAGATE: (B) the FORWARD pass L→R and (C) the BACKWARD pass R→L ────────────────────
    # ⛔ TRAPS: a-comment-quoted-as-a-finding: ONE pass each, in chain order, which on a chain IS
    # forward-backward. It is not an iterative scheme, and a source comment's shorthand once crossed
    # into a design doc as if it were one. When both passes end every node holds one message from each
    # neighbour it has (owner ruling 2026-09-04): the recipient's kernel wrote it, or SILENCE stands.
    from_left = _pass(order_list, ctx.left_list, prepared, backward=False)
    from_right = _pass(order_list[::-1], ctx.right_list, prepared, backward=True)

    # ── PHASE 2, SOLVE: (D) the policy's half — the two held messages into ψ's channels ──────────────
    msg = prepared.solve(from_left, from_right)
    counts = AssertionCounts()
    _check_message(msg, ctx, counts)

    # ⭐⭐ THE CITIZENSHIP SEAM (owner ruling 2026-08-25): a delivered certified-flux claim joins the
    # λ-factor rows in the FINAL solve only. Phase-A (`build_region_init`) and the own-evidence
    # precision never see it — an imputation may inform the fused answer, never masquerade as the
    # slot's own evidence.
    dc_fin = _psi(
        global_lp,
        msg,
        fg_ref=f_g,
        fpos_ref=f_pos,
        fneg_ref=f_neg,
        extra_lam_rows=msg.lam_rows,
        cube_rows=msg.cube_rows,
    )

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
            priors=None,
        ).gdna_frac
        # the message-free self-solve variances, for the local-error attribution — a debug-only solve, so
        # the production path carries none of it (the self-solve fractions come from ``own``).
        # ⚠ Deliberately AFTER the write-back, so its reference is the OUTGOING belief. That is what the
        # shipped solver did, and an instrument comparing against it depends on the same reference.
        _dc_loc = _psi(global_lp, PsiMessage.silent(), fg_ref=f_g, fpos_ref=f_pos, fneg_ref=f_neg)
        _capture.update(
            backbone_assertions=counts,
            order=np.asarray(order_list, np.int64),
            n_slot=n_slot.copy(),
            left=np.asarray(left, np.int64),
            right=np.asarray(right, np.int64),
            # ⭐ which policy ran, read off the artifact — the witness an instrument's "the arm ran"
            # assertion needs (TRAPS: an-ablation-that-never-ran), never a config flag it did not thread
            policy_name=str(getattr(policy, "name", type(policy).__name__)),
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
            # replay reproduces the shipped f_g exactly BEFORE ablating. ⛔ The final solve's row
            # factor is ``intron_prior`` PLUS the delivered certified-flux rows (`msg.lam_rows`) —
            # a replay passing the bare ``intron_prior`` is unfaithful whenever the stream is live,
            # so both are published and a faithful replay sums them.
            fg_init=_fg_init,
            fpos_init=_fp_init,
            fneg_init=_fn_init,
            intron_prior=intron_prior,
            lam_rows=msg.lam_rows,
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


def _pass(seq, nbr, prepared, *, backward: bool) -> list:
    """ONE directional pass — phase 1 of the two-phase solve: for each node in ``seq``, in chain order,
    the node RECEIVES from its neighbour of the other kind and holds the result.

    ⭐ **The whole direction dependence is which neighbour array is read.** ``-1`` is a reference terminal:
    the node holds ``NO_NEIGHBOUR`` (``None``) there, which is not a message — the chain's two end nodes
    hold one message, every other node two. ⛔ **A real hop must arrive**: a kernel that returns ``None``
    for a node that HAS a neighbour is refused, because the solve could not then tell "nothing to say"
    (:data:`~.messages.SILENCE`) from "never spoken to". A policy that sends nothing at all returns no
    kernel, and every node then holds SILENCE from this side.
    """
    receive = prepared.propagate(backward=backward)
    held: list = [None] * len(seq)
    for i in seq:
        s = nbr[i]
        if s < 0:
            continue
        held[i] = SILENCE if receive is None else receive(s, i)
        if held[i] is None:
            raise AssertionError(
                f"the {'backward' if backward else 'forward'} pass left slot {i} with no message from "
                f"its neighbour {s}: a hop that carries nothing must still ARRIVE as SILENCE (owner "
                "ruling 2026-09-04) — return SILENCE, never None, from a real hop"
            )
    return held


# ──────────────────────────────────────────────────────────────────────────────────────────────────────
# THE OUTPUT CONTRACT — projecting the chain belief back onto the two payload axes.
# Not part of the solve: this is what ``CalibrationResult`` / ``priors`` / ``derive`` consume.
# ──────────────────────────────────────────────────────────────────────────────────────────────────────


def chain_region_deconv(chain: RegionChain, belief: RegionBelief, substrate) -> RegionDeconv:
    """Project the chain belief's REGION slots back onto the REGION axis as a :class:`RegionDeconv` — what
    ``CalibrationResult`` / ``priors`` / ``derive`` consume.

    ⚠ **A region's contained population carries no spliced term any more, and that is structural**: the
    accumulator credits ``region_contained`` only when the fragment used no sj, so a contained
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
