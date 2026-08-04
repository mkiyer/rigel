"""The belief-propagation calibration solver: per-node gDNA-vs-RNA deconvolution by a bidirectional sweep.

Deconvolves each node's unspliced fragment mass into a pie ``(f_pos, f_neg, f_g)`` — sense-RNA /
antisense-RNA / gDNA — over the unified region↔boundary chain (`node_chain`), by a single forward-backward
(L→R then R→L) belief-propagation pass (exact on the chain, a forest of linear paths). Each per-node solve
(`simplex_logodds`, the log-density log-odds solver) reconciles three sources of information: the intrinsic
strand likelihood (the Beta-Binomial tilt — the only signal a count carries), the cross-node imputation
messages, and the population gDNA prior. Its two theoretical legs are the count-zero-information
principle and the composition (enrichment-ratio) message model.


**The UNIFIED (composition) solver — one message mode.** Each message is REFRAMED into the destination's
frame by the enrichment ratio ``r = ρ_tot(dst)/ρ_tot(src)`` (`node_total_density`, lazy + composition-aware,
per-side spliced), the inactive strands are FILTERED at the destination, the mature is ROUTED (graft into an
exon / peel out of an exon, via the boundary's measured spliced), then each component becomes a factor on the
destination's ``f_c`` by the density mode ÷ its own observed mass ``M_dst`` (the pure arithmetic in
`enrichment_frame`). The forward-backward relay carries per-component densities + precisions; the combine
transports both neighbours into the node's frame and runs the ψ solve. The message VARIANCE model is
COMPLETE — laws M1–M11 derived, MC-validated (`scripts/debug/message_variance_mc.py`) and A/B-won. The
density-uniformity NPMLE proxy it replaced is retired.

Module layout. The per-node geometry / belief / statics / init primitives and the pure geometry helpers
(`build_node_geometry`, `build_node_statics`, `init_beliefs`, `node_global_geometry`, `node_total_density`,
`NodeGeometry`/`NodeBelief`/`NodeStatics`) live in the lower `node_geometry` module; the per-node
INITIALIZATION self-solve (the four sources → each node's own ``(density, precision)``) lives in `node_init`
(`build_node_init`). Both are re-exported here for the calibrator's convenience; this module owns:
* `node_sweep` — the single forward-backward unified sweep (build the self-solve → relay → combine → ψ solve).
* `chain_node_deconv` / `chain_edge_deconv` — project the converged belief back onto the node axis and
  the contiguous-edge axis, which is what `CalibrationResult` consumes.
"""

from __future__ import annotations

import math

import numpy as np

from .enrichment_frame import (
    _fmax,
    composition_logvar,
    graft_frame_logvar,
    graft_frame_logvar_scalar,
    graft_premise_logvar,
    mismatch_deflate,
    mismatch_gap,
    peel_continue_share,
    peel_continue_share_scalar,
    peel_share_logvar,
    residual_level,
    residual_level_scalar,
    transfer_logvar,
)
from .node_chain import EDGE, NODE, NodeChain
from .signature import coarse_type_array
from .node_geometry import (
    NodeBelief,
    NodeGeometry,
    NodeStatics,
    build_node_geometry,
    build_node_statics,
    g1_locked,
    init_beliefs,
    node_global_geometry,
    node_total_density,
    terminus_flank_gain,
)
from .node_init import build_node_init, own_composition_logvar
from .simplex_logodds import (
    _logodds_grid,
    _solve_nodes_logodds_all,
)
from .strand_deconv import NodeDeconv

__all__ = [
    "NodeGeometry",
    "build_node_geometry",
    "NodeBelief",
    "NodeStatics",
    "build_node_statics",
    "init_beliefs",
    "node_sweep",
    "node_global_geometry",
    "node_total_density",
    "chain_node_deconv",
    "chain_edge_deconv",
]

_EPS = 1.0e-9

# ⭐ ONE ρ-iteration (straight-line, no loop — see the combine below). The combine's reframe needs the DESTINATION's composition, which is what the messages
# are trying to determine — so a second iteration reframed off the FUSED posterior, i.e. the message from L
# into i was framed using a belief that already contained the message from R into i. That is a BP violation
# (a message may depend on anything except its own destination's other messages), and it was the last of
# three: `_far` (the destination's other NEIGHBOUR) and `_p_nu_own` (the destination's OWN evidence) are
# already deleted. With ONE iteration the frame comes from `_init_belief()`, which is a pure function of the
# chain / substrate / statics and the three fitted library scalars — belief-free — so no posterior feeds back
# and the chain's forward-backward pass is exact BP.
#
# PRICE, measured over all 32 conditions: refit=0 **−0.0002** (4 better / 2 worse / 26 flat — removing it is
# slightly BETTER), refit=1 **+0.0011** (4/4/24), concentrated in unstranded × capON (+0.0104, 0 better /
# 3 worse). Held-fixed-node-set z2 is flat (ALL 8.53 → 8.50; exon-AMBIG 59.4 → 54.3 and boundary 4.27 → 4.09
# better, exon-single 3.57 → 3.84 worse). Fit-substrate mwae +0.0002.
#
# ⛔ That price is also the answer to "does the elegant refactor earn itself?" — NO. Carrying each message as
# a FUNCTION of the destination's state (a pairwise potential ψ(x_s, x_i) evaluated inside the destination's
# own ψ solve) is the principled way to get the composition-aware frame back legally, but it is a rewrite of
# the message representation and it is bidding for ~0.001. Revisit only if the frame becomes load-bearing for
# another reason. A belief-FREE constant frame (`node_total_density`'s pure-gDNA fallback) was also
# prototyped: it costs more (+0.0003 suite, +0.0010 fit-substrate) and is neutral-to-worse on a held-fixed
# z2, so it is not a better trade than simply dropping the iteration.


def node_sweep(
    chain: NodeChain,
    statics: NodeStatics,
    geometry: NodeGeometry,
    belief: NodeBelief,
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
    length_loglik=None,
    _capture: dict | None = None,
):
    """The belief-propagation sweep over the chain — gDNA AND per-strand RNA messages, in COUNT space (no
    ``(M/E)²`` Jacobian) — as a single **FORWARD-BACKWARD** pass.

    The chain is a forest of linear paths, so BP is exact in one forward + one backward pass (vs Gauss-Seidel /
    Jacobi which propagate one hop per pass). A message's precision is the source's own HONEST belief precision
    degraded by the two independent defects a cross-node imputation suffers


        p = 1 / ( Var(log f_c^src) + 1/n_src  +  σ²_transfer  +  b̂² )
                 \\__ strand ___/   \\_count_/    \\_ SCALE _/    \\_ COMPOSITION _/

    — the composition and count precision the source actually earned, the reframe's own scale uncertainty
    ``σ²_transfer = Var(log r)`` (M5; 0 on the matched graft where ``r`` is common-mode), and the
    DerSimonian–Laird estimate of how wrong the imputation PREMISE itself is (``b̂²``, M7 — the message's
    composition against the destination's independent self-solve). All four are computed inside the pass from
    counts and effective lengths: nothing is fitted, there is no precision to refit and no outer fixed-point
    loop. The global prior is ANCHORED (every input fit once before the solve), so the single FB pass is exact.

    BEFORE the pass: the population gDNA hyperprior (:class:`~.gdna_landscape.GdnaLandscape`), fit once on
    the pass-0 result, and passed as ``gdna_prior``. (It replaced the δ-pin :class:`~.npmle.DensityNPMLE` in
    this role in W4; the NPMLE is retired here but still fitted for the Role-A enrichment landscape and QC.)
    ``gdna_prior=None`` is a first-class PRIOR-FREE solve: ψ then
    carries the derived reference alone on both arms (``simplex_logodds._gdna_arm`` / ``_rna_arm``) — prior-free
    is not reference-free.

    The pass (`_unified_solve`): (A) the per-node message-free SELF-SOLVE — each node's own per-component
    ``(density, precision)`` from the four init sources (`node_init.build_node_init`); (B) a FORWARD relay L→R
    fusing each node's own belief with the reframed upstream context (a thin node has own precision 0 and passes
    the context through); (C) a BACKWARD relay R→L; (D) the COMBINE — transport both neighbour messages into the
    node's frame (reframe → route mature → ÷``M_dst``), fuse, and run one batched ψ solve.

    The per-node ψ solve integrates: the strand likelihood, the Jeffreys reference (single-strand nodes), the
    global NB prior (ALL solvable nodes), the gDNA intron-factory λ-factor, and the gDNA + per-strand RNA
    imputation factors from the combine. Each RNA strand flows only where that strand is continuous on both
    endpoints (``free_s``); gDNA is genomically continuous and strand-agnostic. Only G2/G3 nodes with data are
    written; G1 sinks + empty nodes keep their signature-binary init. Returns the resolved
    :class:`NodeBelief`."""
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    fp, fn = statics.free_pos, statics.free_neg
    f_pos = np.asarray(belief.f_pos, dtype=np.float64).copy()
    f_neg = np.asarray(belief.f_neg, dtype=np.float64).copy()
    f_g = np.asarray(belief.f_g, dtype=np.float64).copy()
    # precision state (Phase 1: computed + carried; consumed by the honest message send in Phase 2).
    var_pos = np.asarray(belief.var_pos, dtype=np.float64).copy()
    var_neg = np.asarray(belief.var_neg, dtype=np.float64).copy()
    var_g = np.asarray(belief.var_gdna, dtype=np.float64).copy()
    # the INCOMING belief, kept for the diagnostic capture: it is the ``fg_ref`` the final solve freezes its
    # variance at, so a channel-ablation replay must pass the SAME reference to be faithful (inert in prod).
    _fg_init, _fp_init, _fn_init = f_g.copy(), f_pos.copy(), f_neg.copy()

    # ── ⭐ THE FACES ARE GONE (S5.e). ──────────────────────────────────────────────────────────────
    # These were six ``(left, right)`` tuples indexed by face, because a boundary's two sides lay in
    # differently-sized flanks and therefore had different divisors. A contiguous edge is a 0-bp line
    # with ONE set of numbers, so each is a single array and every ``df``/``sf`` face parameter
    # threaded through `_relay`, `_transport`, `_peel_share` and `_seam_pair` disappears with them.
    # ⚠ What does NOT disappear is DIRECTION: the forward and backward relays are still distinct, and
    # a junction's donor and acceptor are still different lines. A face was a property of the GEOMETRY;
    # a direction is a property of the MESSAGE.
    EG = np.asarray(geometry.eff_gdna, np.float64)
    ER = np.asarray(geometry.eff_rna, np.float64)
    ESP = np.asarray(geometry.eff_junction, np.float64)  # [n, 2] by TRANSCRIPT strand
    SPL = np.asarray(geometry.junction_count, np.float64)  # [n, 2] by TRANSCRIPT strand
    CNT = np.asarray(geometry.unspliced_count, np.float64)  # [n, 2] by GENOME strand
    # the unspliced count is BOTH the density numerator and the Poisson n (S5.e) — one number, not the
    # old fractional ``mass`` plus a separate integer ``flux``.
    n_slot = CNT.sum(axis=1)
    u_pos, u_neg = CNT[:, 0], CNT[:, 1]
    spliced_slot = np.asarray(geometry.spliced_count, np.float64).sum(axis=1)
    # per-slot "global" gDNA support — the basis the pass-0 rate prior is fit + projected on.
    mass_global, eff_global = node_global_geometry(geometry)

    # The per-node solve is the log-density 1-D/2-D log-odds solver (simplex_logodds, O(m·K),
    # genome-scale-tractable). The "solve grid" is the f_g axis the global NB prior is evaluated on (the
    # log-odds σ(λ) lattice).
    _, solve_grid = _logodds_grid(int(n_grid), float(logodds_window))
    kappa = float(rna_sense_frac)
    od_g, od_r = gdna_strand_overdispersion, rna_strand_overdispersion

    def _local_solve(g_arr, gm=None, gp=None, rm=None, rp=None, lam_imp=None, theta_imp=None):
        """The per-node local/final solve (log-density log-odds backend). Returns the :class:`NodeDeconv`
        (the readout ``*_frac``/``*_frac_var`` + the free-coordinate seed ``lam_mean``/``lam_var``/
        ``theta_mean``/``theta_var``). Phase A calls it message-free; phase D passes the FB messages.
        ``lam_imp``/``theta_imp`` are the SINGLE-λ composition message + the θ tilt message (2-tuples of
        ``(mode, prec)``), the rank-1 fix that replaces the two per-component ``gm``/``rm`` messages."""
        return _solve_nodes_logodds_all(
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
            gdna_imp_mode=gm,
            gdna_imp_prec=gp,
            rna_imp_mode=rm,
            rna_imp_prec=rp,
            lam_imp_mode=None if lam_imp is None else lam_imp[0],
            lam_imp_prec=None if lam_imp is None else lam_imp[1],
            theta_imp_mode=None if theta_imp is None else theta_imp[0],
            theta_imp_prec=None if theta_imp is None else theta_imp[1],
            # the gDNA intron factory λ-factor (anchored, per-intron, 0 elsewhere): peels confident gDNA from
            # introns against the intergenic background BEFORE the sweep resolves the pie (design §9.3). Added
            # to ψ, distinct from the gDNA arm; participates in the local solve AND the relay (a confident
            # intron gDNA belief propagates genomically to the flanking exons/boundaries).
            lam_logprior=intron_prior,
            # ⭐ the FRAGMENT-LENGTH λ-factor (`length_likelihood`, P2). It enters the LOCAL solve and the
            # FINAL one, exactly like the intron factory, so a node's own length evidence both sets its
            # belief and propagates through the relay. ``None`` ⇒ byte-identical to the pre-P2 path.
            length_loglik=length_loglik,
            # count-zero-info variance freeze (§2, B1): reference = the incoming belief, so the variance —
            # hence the message precision — is evaluated near the truth, not at a flat ½.
            fg_ref=f_g,
            fpos_ref=f_pos,
            fneg_ref=f_neg,
        )

    # Two distinct gates, both from the region SIGNATURE (never the counts — count-zero-info):
    #   * SOLVE gate (`solvable`): a node deconvolves its own gDNA-vs-RNA split iff it admits ≥1 RNA strand and
    #     has unspliced mass. A G1 node — no admissible RNA strand: an intergenic region, or a gene-boundary seam
    #     (TSS/TES, opposite-strand exon↔exon) — is a LOCKED all-gDNA node ({0,0,1}); it is not solved, keeping
    #     its init (RNA cannot cross a gene boundary, so its unspliced mass is purely gDNA).
    #   * EMISSION gate (per component, in `_scan`): which MESSAGES a node sends. A three-term Boolean over the
    #     components gDNA / +RNA / −RNA, structural and symmetric — defined at the top of the scan loop.
    solvable = (fp | fn) & (n_slot > 0.0)

    # THE gDNA ARM of ψ — the COMPOSITION prior. The NPMLE's two roles are kept SEPARATE
    # this ``gdna_prior`` is the COMPOSITION arm ONLY.
    #   * ``gdna_prior=None`` — the INITIAL prior-free solve: the arm is the derived inert reference (½·log f_c
    #     on both arms, added by simplex_logodds). A total-density NPMLE is an ENRICHMENT model, NOT a DNA
    #     composition prior — letting it vote a node's f_g is the count-votes-composition regression (§0).
    #   * ``gdna_prior`` set — a REFIT with the FITTED gDNA hyperprior (fit on the DECONVOLVED gDNA density,
    #     with ρ_bg pinned as a smooth low-density component — NO clamp/floor). ANCHORED, EXTREMELY WEAK.
    global_lp = (
        gdna_prior.logprior(solve_grid, mass_global, eff_global) if gdna_prior is not None else None
    )

    # (RETIRED: the belief-free NPMLE-PROJECTION σ²_transfer — ``var_proj[dst] + (μ_proj[dst]−μ_proj[src])²``,
    # a density-uniformity proxy that hybrid capture invalidates, and which was identically 0 in pass-0 anyway.
    # σ²_transfer is now the DERIVED ``Var(log r)`` (M5, from ``composition_logvar``) and the composition half of
    # the cliff cost is the DL ``b̂²`` — both computed inside the pass from counts and eff-lengths, so no
    # enrichment prior enters the solver at all. The ENRICHMENT NPMLE itself is still fit in `calibrate` for the
    # QC report + the toy-injection substrate; it simply no longer feeds message precision.)

    # Genomic slot order for the forward/backward scans (within each ref path; left/right break at −1).
    # ⭐ Slot ids ARE the genomic visiting order, so the order is ``arange`` and the chain does not store
    # it. The scans are sequential, so iterate as a Python list of ints (faster than numpy indexing).
    order_list = list(range(int(chain.n_slots)))
    # Per-node EXON-region flag (coarse_type == 2) — the unified relay routes mature into EXON destinations
    # (the graft) and peels it out of EXON sources (`ex_a` in `_unified_solve`).
    _rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    _ri = np.clip(np.asarray(chain.obj_idx, dtype=np.int64), 0, _rtype.shape[0] - 1)
    is_exon_node = ((np.asarray(chain.kind) == NODE) & (_rtype[_ri] == 2)).tolist()

    # ─────────────────────────────────────────────────────────────────────────────────────────────────────
    # THE UNIFIED SOLVER (owner 2026-07-23) — ONE mode: reframe → filter → route → ÷M.
    # ─────────────────────────────────────────────────────────────────────────────────────────────────────
    _uni_msg = (
        None,
        None,
        None,
        None,
        None,
        None,
    )  # the published imputation factors (for the diagnostics)

    # The per-node message-free self-solve — the four init sources (`node_init.build_node_init`): MEASURED /
    # intron-factory / strand-deconvolution / unsolved-default. Its ``rho_*`` / ``prec_*`` seed the relay.
    _ni = build_node_init(
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
        length_loglik=length_loglik,
    )

    def _unified_solve():
        n = f_g.shape[0]
        is_reg_a = np.asarray(chain.kind) == NODE
        is_bnd_a = ~is_reg_a
        ex_a = np.asarray(
            is_exon_node, dtype=bool
        )  # source-is-exon selector for the mature routing
        fp_a, fn_a = np.asarray(fp, bool), np.asarray(fn, bool)
        is_amb = (
            fp_a & fn_a
        )  # AMBIG: both strands live → the tilt θ is a free DOF (its own message)
        M = np.asarray(mass_global, np.float64)
        E_g = np.asarray(eff_global, np.float64)
        E_r = ER  # per-slot RNA-FL eff-length — one number, no faces to sum

        # ── M5 σ²_transfer = Var(log r): the per-node total-density log-variance
        # (`enrichment_frame.composition_logvar`). Var(log ρ_tot) = 1/n + [(1/E_g−1/E_r)/B]²·Var(f_g), with
        # Var(f_g) = (f_g(1−f_g))²/τ_λ CAPPED at f_g(1−f_g) (a fraction's max variance) and 0 for a
        # composition-certain (struct_lock) node. Evaluated at the INPUT belief f_g (consistent with the reframe
        # density rho0). This is the HONEST, prior-free transfer variance that RETIRES the
        # density-uniformity NPMLE proxy (which was identically 0 in pass-0, so pass-0 had NO transfer damping).
        _n_node = n_slot  # the count IS the Poisson n; no separate flux to pool over faces
        _fg0 = np.clip(np.asarray(f_g, np.float64), 0.0, 1.0)
        _fgfr = _fg0 * (1.0 - _fg0)
        _tau = np.asarray(_ni.tau_lam, np.float64)
        _var_fg = np.where(
            np.asarray(_ni.struct_lock, bool),
            0.0,
            np.where(_tau > _EPS, np.minimum(_fgfr * _fgfr / np.maximum(_tau, _EPS), _fgfr), _fgfr),
        )
        logvar_tot = np.asarray(composition_logvar(_fg0, E_g, E_r, _var_fg, _n_node), np.float64)
        # ── THE CROSS-CLIFF PRECISION: σ²_transfer (M5, the SCALE) + b̂² (the DL COMPOSITION MISMATCH) ────────
        # A message's delivered mode error decomposes EXACTLY into two orthogonal defects (MC-validated to
        # machine precision, `scripts/debug/message_variance_mc.py` M7a/M7b):
        #     mo_c − log f_c^dst,true  =  log(s_c^src / s_c^dst,true)  +  log(r̂/r_true)
        #                                 └─ the composition-SHARE mismatch      └─ the reframe's own scale noise
        # so   σ²_c,delivered = v_src,c + σ²_transfer + b_c² .  The two terms are NOT interchangeable:
        #   * σ²_transfer = Var(log r) (M5, `transfer_logvar`) prices the SCALE sampling — 0 on the matched-set
        #     graft (r is common-mode there and cancels in the composition), load-bearing on peel/partial;
        #   * b_c² prices the COMPOSITION drift — the imputation premise ("neighbours share a composition") being
        #     wrong. This is the term a pure enrichment cliff must NOT pay: `(log r)²` (the shipped proxy this
        #     replaces) charged the WHOLE cliff as mismatch, which recovers the stranded arm but over-damps the
        #     extreme-capture arm, where the composition really is preserved across a 1000× enrichment step.
        # b_c is a population (third-source) quantity we do not have prior-free — but the destination has an
        # INDEPENDENT estimate of the same composition: its own message-free self-solve (`_ni`). Treating the
        # message and the own belief as two studies of one quantity, the DerSimonian–Laird between-source
        # estimator recovers b² with NO tuned constant (`_dl` below). Applied in `_transport`.

        # ── ⭐⭐⭐ THE gDNA SCALE RULE: the reframe is a COMPOSITION imputation, and it needs a LICENCE ───
        # Substituting ``rho_c(src) = phi_c(src)·rho_tot(src)`` into the reframe shows what ``r`` actually
        # delivers:
        #       rho_c^msg(dst) = rho_c(src)·rho_tot(dst)/rho_tot(src) = phi_c(src)·rho_tot(dst)
        # — the source's density SHARE applied to the destination's observed total. ⭐ So the reframe is a
        # PURE COMPOSITION IMPUTATION with no level transport in it at all: exact iff the two shares match,
        # and wrong by exactly ``phi_c(src)/phi_c(dst)``. ⛔ When the source carries only gDNA that ratio is
        # ``1/phi_g(dst)`` and the delivered level collapses to ``rho_tot(dst)`` — the destination's OWN
        # total, INDEPENDENTLY of the source's measurement. That is `TRAPS.md` D4 exactly, and it is why the
        # answer at an evidence-free exon sat flat at f_g ≈ 0.91 across a 1000x RNA sweep while its truth
        # spanned 1.00 → 0.001; the measured 28.684 delivered against a true 0.0257 is ``M/E_g`` to the
        # digit, not an approximation.
        #
        # So a message makes TWO claims and they need two scales. A COMPOSITION crosses by ``r``, licensed
        # by the predicate the λ-emission gate already uses — the source must have SUPPLIED both components
        # of the pair. A source carrying one component has no composition to lend, and its authority is a
        # LEVEL. ⭐⭐ For gDNA a level crosses **UNSCALED**, and that is not a patch, it is this project's
        # base assumption: gDNA is uniform along the genome before capture, so its level does not change
        # across an EDGE.
        #
        # ⭐⭐⭐ AND THAT IS THE WHOLE RULE FOR BOTH CAPTURE ARMS, WITH NO BRANCH — because capture is
        # carried by the pure-gDNA population's OWN measurements rather than by a scale factor. The relay's
        # mass pin below re-scales the running level to each object's own observed total, and at a
        # structurally pure-gDNA object that total IS its gDNA density, measured at its own capture stratum.
        # So the level handed to an exon is its flanking EDGE's own enriched measurement, not the
        # off-target floor and not the exon's total. A pooled per-capture-class landscape ratio was built to do this
        # explicitly and measured **byte-identical off capture** and worth 1.2 % of one class on one
        # capture-ON condition; it is deleted, and `node_geometry` records why so it is not rebuilt.
        #
        # ⭐ It also fixes the INTRON-mediated face of the same defect, which was not obvious in advance and
        # is a measurement rather than an argument: a factory-solved intron reports f_g ≈ 1, hence ~zero RNA
        # density, hence zero RNA precision — so it cannot lend a composition either. The toy gate that used
        # to PIN that behaviour (`test_toy_harness`, an exon inheriting its neighbouring intron's
        # composition, 31x on the ladder donor) failed because the dependence is gone: 1.04x.
        # ⛔ What is left unfixed is narrower than it first looked — a source with a real MIXED composition
        # that differs from the destination's, which needs genuine RNA evidence at the source (nascent plus
        # strand or length). That is the premise being wrong rather than absent: M7's DerSimonian-Laird b̂²,
        # which needs the destination to have its own self-solve and so is silent on unstranded data.
        #
        # ⛔ And the residual under capture: the level delivered to an exon is its flanking EDGE's, and a
        # fragment spanning that EDGE lies only PARTLY under the probe while one contained in the exon lies
        # wholly under it, so it is a LOWER bound (for any non-decreasing probe-retention law,
        # off-probe ≤ spanning an EDGE ≤ on-probe). Closing that needs a
        # probe-geometry model the tool does not have. ⛔ The affine-in-overlap extrapolation
        # ``e_g[exon] = 2·e_g[EDGE] − e_g[off-probe]`` is EXACTLY the simulator's own retention law, so fitting
        # it would be scoring against the substrate that generated it. NOT built.

        # ── the per-TRANSCRIPT-STRAND mature (junction) DENSITY at each line ──────────────────────────
        # ⭐ ONE array per strand, not a (left face, right face) pair. The predecessor needed the pair
        # because it could not tell a donor flank from an acceptor flank without guessing from the
        # signature bits; the v8 index states ``(src, dst, strand)`` and `build_node_geometry` has
        # already placed the flux on the lines the junction actually leaves and enters.
        # ⚠ And the ``accept_l``/``accept_r`` test goes with them: a line either carries mature flux or
        # it does not, and ``mature_count == 0`` says so directly.
        def _mature_rho(strand: int) -> np.ndarray:
            c, e = SPL[:, strand], ESP[:, strand]
            live = (c > _EPS) & (e > _EPS)
            return np.where(live, c / np.where(live, e, 1.0), 0.0)

        spl_p = _mature_rho(0)  # + transcript mature density at this line (0 on NODE slots)
        spl_n = _mature_rho(1)

        # ── P1d: the graft's RNA claim is a LOWER BOUND, and the tightest one is the DOMINANT flank ──────
        # A boundary hands the exon next door ``ρ_ν + ρ_μ`` as *the* exon RNA density. Every molecule
        # counted there is in the exon — but the exon may also hold molecules that never touch this seam:
        # ones that reach it by the OTHER flank, or that start/end inside it. So what the graft actually
        # knows is an INEQUALITY,
        #     ρ_R(exon) ≥ ρ_ν(B) + ρ_μ(B)       for every boundary B flanking the exon,
        # and it uses it as an equality. That — not the 100 bp → 2,100 bp extrapolation the term was
        # originally scoped as — is the graft's premise error: measured FLAT in the extrapolation ratio
        # (1.13× over a 6.7× range) and flat in exon length, but 30–100× larger at boundaries carrying a
        # transcript terminus, i.e. exactly where RNA does NOT flow through the seam.
        #
        # Both flanks bound the SAME quantity, so the tightest available bound is the larger flux. This is
        # a structural fact about which molecules cross a seam, with no fitted constant: the alternatives
        # all fail in the direction their algebra predicts — ``min`` is worse than either bound alone
        # (0.514 vs 0.487), ``sum`` over-states by +0.61 nats as a sum of two bounds must, and doubling the
        # arriving flux does nothing (0.503), so it is specifically the OTHER seam's VALUE that is missing.
        _li_a, _ri_a = np.asarray(left, np.int64), np.asarray(right, np.int64)
        _vl_a, _vr_a = _li_a >= 0, _ri_a >= 0
        _sl_a, _sr_a = np.clip(_li_a, 0, n - 1), np.clip(_ri_a, 0, n - 1)

        # ── ⭐⭐⭐ THE POPULATION HALF OF THE COMPOSITION LICENCE: "is the source measuring the same thing
        # I am?" ──────────────────────────────────────────────────────────────────────────────────────────
        # An EDGE counts what spans it CONTIGUOUSLY, so `T(EDGE) = T(left) ∩ T(right)` and a transcript
        # TERMINUS at the EDGE makes one flank's RNA population strictly larger. Where two objects do not
        # measure the same population, the density discrepancy between them is not attributable to
        # hybrid-capture enrichment — enrichment and a population difference are indistinguishable — so the
        # composition may not be imputed across the step, the gDNA LEVEL crosses UNSCALED, and the mass pin
        # is off with it (the identity is only implied under the premise).
        # ⛔ TERMINI ONLY: a DONOR/ACCEPTOR EDGE also changes the population, but there the flux is MEASURED
        # and the graft and the peel exist to route it. `terminus_flank_gain` carries that scope and the
        # genomic-versus-TSS/TES reasoning; deriving it there rather than here is what keeps ONE home for
        # the predicate.
        # THE PAIR ALGEBRA: the chain strictly alternates, so exactly one slot of an adjacent pair is an
        # EDGE, and the pair `(i, left[i])` IS the pair `(left[i], right[left[i]])`. So one array answers
        # the question for the left-hand step of every slot and the other is that array read through
        # `right`. When `i` is the EDGE its left neighbour is its LEFT flank; when `i` is a NODE it is the
        # RIGHT flank of the EDGE on its left.
        _rgain, _lgain = terminus_flank_gain(statics.edge_flags)
        _pop_l_a = np.where(is_bnd_a, ~_lgain, ~_rgain[_sl_a]) & _vl_a
        _pop_r_a = np.where(_vr_a, _pop_l_a[_sr_a], False)

        # ⚠ The two seams must be compared IN THE DESTINATION'S OWN FRAME. Each is lifted by its own
        # enrichment step ``r``, and under capture those steps differ (the two seams sit at different probe
        # positions), so comparing unlifted fluxes reads a capture difference as an abundance difference.
        # Measured on the MODE variant of this idea (replace the claim by the dominant seam's flux): raw is
        # uniformly better off-capture (14/14 conditions, 0 worse) and uniformly worse on (8 worse); lifting
        # first halves the on-capture damage. That MODE fix is NOT landed — under capture the graft already
        # OVER-states (median φ = 2.45, an M8 frame problem), so no bound-tightening can help there — but the
        # same seam pair is what measures the variance below, and it must be framed for that too.
        def _flank_dom(rho, spf):
            """Per-slot: the flux each of its two flanking EDGES sends it, ALREADY lifted into this
            slot's frame, per strand.

            ⭐ The predecessor took four arguments — a ρ per face and a flux per face — and paired them
            crosswise (``spf[1]`` is what a LEFT neighbour presents, ``spf[0]`` what a RIGHT one does).
            With one ρ and one flux per slot the pairing is the identity and only the DIRECTION remains.
            Zero wherever a flank is absent or is not an edge, so on a one-junction exon this degenerates
            exactly to that junction's own lifted flux."""

            def _side(ok, nb):
                r = np.where(
                    ok & (rho[nb] > _EPS) & (rho > _EPS), rho / np.maximum(rho[nb], _EPS), 1.0
                )
                return np.where(ok & is_bnd_a[nb], spf[nb] * r, 0.0)

            return _side(_vl_a, _sl_a), _side(_vr_a, _sr_a)

        # own per-component densities + precisions — the message-free SELF-SOLVE (`node_init.build_node_init`,
        # the four sources). ``rho_*`` are the own densities; ``prec_*`` combine the strand + intron-factory
        # composition evidence ``τ_λ`` (Var(log f_g)=(1−f_g)²/τ, Var(log f_r)=f_g²/τ) with the Poisson count
        # power — 0 at an uninformed node (the honest unstranded statement). Both RNA arms are built from the
        # SAME τ as gDNA (never the local posterior variance, which pools the shared reference into a phantom).
        og, op, on = _ni.rho_g, _ni.rho_pos, _ni.rho_neg
        pg_own, pp_own, pn_own = _ni.prec_g, _ni.prec_pos, _ni.prec_neg

        # ── THREE-STREAM precision seeds (the single-λ combine)
        # §6). The density precision (pg/pp/pn, unchanged) drives the MODE fusion. TWO extra accumulators
        # SEPARATE the RANK-1 composition (→ ONE λ-message) from the INDEPENDENT measurements (→ their own
        # gdna/rna messages), so ψ counts each ONCE — the fix for the two-message double-count:
        #   * τ (COMPOSITION): the Schur λ-precision; 0 at anchors + unstranded non-factory nodes. On unstranded
        #     data τ≈0 ⇒ the λ-message is ~off and the solve rides the measurement channels (preserving the M5
        #     win); on stranded data the composition rides ONE λ-message (no double-count).
        #   * mg/mp/mn (MEASUREMENT): the anchor gDNA count (own, struct_lock only) + the spliced RNA count
        #     (added at the graft). A density-level bound, NOT a composition vote.
        _struct = np.asarray(_ni.struct_lock, bool)
        tau_own = np.asarray(_ni.tau_lam, np.float64)
        mg_own = np.where(_struct, np.asarray(pg_own, np.float64), 0.0)
        mp_own = np.zeros_like(np.asarray(pp_own, np.float64))
        mn_own = np.zeros_like(np.asarray(pn_own, np.float64))

        v_own_g, v_own_r = own_composition_logvar(_ni.f_g, tau_own, _struct)
        # the same three states on the λ axis, for the single-DOF composition stream (Var(λ) = 1/τ_λ — the two
        # per-component arms are perfectly anti-correlated, so their Jacobians sum to 1 and the count cancels).
        v_own_lam = np.where(
            _struct, 0.0, np.where(tau_own > _EPS, 1.0 / np.maximum(tau_own, _EPS), np.inf)
        )

        def _fuse(a, pa, b, pb):  # scalar precision-weighted density fuse
            p = pa + pb
            return ((pa * a + pb * b) / p, p) if p > _EPS else (a, 0.0)

        # ── PIN A RELAYED BELIEF TO A COMPOSITION (÷M_dst, design §2 — applied where it is EARNED) ─────────
        # Three independently-relayed densities are NOT a composition: nothing makes ``Σ_c ρ_c·E_c = M``, and
        # the reframe ``r = ρ_tot(dst)/ρ_tot(src)`` cannot supply it either — ``ρ_tot`` is built from the node's
        # *belief* ``f_g`` while the message carries its own, so each hop multiplies in an uncancelled residual.
        # Geometrically the residual telescopes (exon→bnd 0.431 × bnd→exon 2.298 ≈ 1), but its per-hop spread is
        # enormous (p90 = 209 on bnd→exon, ρ_tot itself swings 500× between an intron and its exon), so the
        # relay is a multiplicative RANDOM WALK: the drift's arithmetic mean ratchets to 6–20× over the measured
        # 9-hop mean adopt-run. That is precisely the defect the ``(λ,θ)`` relay retired in ``_scan`` — "the old
        # three-independent-log-fraction relay violated all three (a boundary relayed ``fbp = 51.9``)" — and it
        # returned here as ``Σ_c f_c`` = 75.6 at introns.
        # The fix is the design's own ÷``M_dst``, applied at EVERY node rather than only at the final combine:
        # a belief about node i is a claim about i's OWN observed mass, so scale it to that mass. An UNINFORMED
        # component (precision 0) is supplied by i's own density in the normalizer — that is what keeps a partial
        # message partial (a seam sending gDNA only still gives ``f_g < 1``, §2) instead of renormalizing it into
        # the shift's "the missing component is absent". A structurally-dead strand has own density 0 and so
        # contributes nothing, correctly. Applied per-message in the combine (`_transport`).
        def _pin_v(g, p, n, pg_, pp_, pn_):
            """The message's densities re-expressed in the destination's MASS FRAME: `Σ_c ρ_c·E_c = M`, an
            IDENTITY under the imputation premise, restored by a common factor `k = M/S`.

            ⚠⚠ **THE RESULT IS A COMPARISON FRAME FOR THE MISMATCH TEST — IT IS NOT THE MESSAGE.** It used to
            replace the delivered densities, and that was a belief-propagation violation ON THE MODE: the
            budget `S` fills every component the message does not supply from the DESTINATION'S OWN density
            (`og/op/on`), so the delivered claim became a function of the destination's own self-solve.
            Measured: with an RNA-only message and `E_g = E_r`, the delivered RNA fraction was **exactly
            `1/(1 + f_g_own)`** (verified to 2.1e-16) — the node confirming itself. On unstranded data
            `f_g_own` is the uninformative ½, so the pin reserved **33.6 %** of the mass budget for gDNA the
            message never claimed and a **zero-gDNA library read back 29.3 % gDNA**; on stranded data, where
            the strand likelihood resolves `f_g_own` to 0.013, the reservation was 1.2 % and the
            false-positive rate 1.4 %. The reservation WAS the false-positive rate. Full derivation:


            Feeding the DL mismatch test is what it is legitimately for. That test is M7's two-study
            random-effects comparison **against the destination's own self-solve**, so destination
            information is the point of it — and it sets a VARIANCE, which can mis-weight a message but
            cannot invent a location. Pinning first is what makes the gap a pure COMPOSITION statement: both
            sides then account for the same total, so the common scale (the reframe residual) is gone from
            `G` and only the share drift is left.

            The partial-claim semantics are load-bearing FOR THAT COMPARISON: a component the claim does not
            SUPPLY (precision 0) contributes the node's OWN density to the budget and does not move, which
            keeps a partial claim partial. Rescaling all three blindly instead regresses capture-OFF 3.6×.

            The weighted alternative (M12 — apportion the correction
            by how badly each component is known, of which this common factor is the `w → 0` limit) was
            derived, implemented and A/B'd, and is NOT used: once λ is fused by its own precision (the
            message packet) it is completely INERT — 0 better / 0 worse / 32 flat — and the per-message
            variant is a net loss."""
            sg = np.where(pg_ > 0.0, g, og)
            sp = np.where(pp_ > 0.0, p, op)
            sn = np.where(pn_ > 0.0, n, on)
            s_ = sg * E_g + (sp + sn) * E_r
            k = np.where((s_ > _EPS) & (M > _EPS), M / np.maximum(s_, _EPS), 1.0)
            return g * k, p * k, n * k

        # ── THE DerSimonian–LAIRD COMPOSITION-MISMATCH DEFLATION (`message_variance_mc.py` M7c/M7d) ─────────
        # The message and the destination's own message-free self-solve are two INDEPENDENT estimators of the
        # same composition. Their observed gap G has, under the null "the shares match", variance
        # v_msg + v_own; anything beyond that is real between-source drift — the DerSimonian–Laird
        # (random-effects, method-of-moments) between-study variance:
        #       b̂² = max(0, G² − v_msg − v_own)       →       p_effective = 1 / (v_msg + b̂²)
        # NO tuned constant: the coefficient is 1 because it is a second moment, and the truncation at 0 is the
        # method's own (a negative variance estimate means "no detectable drift"). At b=0 it is positively biased
        # by 0.4839·(v_msg+v_own) — the OVER-damping direction (the count-zero-info-safe one), and harmless
        # because a message that agrees with the own belief moves the fused mode nowhere.
        # ── the RELAY: accumulate the forward/backward context belief (densities), reframed each hop by the
        # INPUT-belief ρ_tot. Returns each slot's context belief IN ITS OWN FRAME; the combine re-reframes.
        # ⭐ ONE ρ_tot per slot. The predecessor returned a triple (node, left face, right face) because
        # only the ACCEPTOR face carried mature; with the faces gone a line either has mature flux or it
        # does not, so ``rho_with_mature`` already is the per-slot answer.
        rho_node0, rho0 = node_total_density(geometry, np.asarray(f_g, np.float64))

        # ── THE PEEL, BY COMPOSITION: the LEVEL (M11 ⊕ the own belief) and the SHARE (M10) ─────────────────
        # A boundary's unspliced crossing is `gDNA + the RNA that CONTINUES`. The exon's RNA arriving at the
        # seam either continues or splices away, so the honest operator is a SCALING by the continuing share
        #     w = ρ_ν/(ρ_ν+ρ_μ)                                            (M10, `peel_continue_share`)
        # not the subtraction `ρ_R·r − ρ_μ` it replaces. `w` is enrichment-free — capture multiplies the
        # continuing and the splicing channels alike, so `e` cancels identically inside the ratio — and a
        # scaling COMMUTES with the scale error the reframe carries, where a difference AMPLIFIES it by
        # `u = 1/w` (measured 1.77× / 2.39× / 5.01× at u = 2/3/10, M10b). That matters because the exon-face
        # reframe error is irreducible: a boundary samples an `fl_mean` window around a point while an exon
        # samples its whole interior, so with mid-exon probes the two sit at genuinely different capture
        # (0.4–1.3 nats, unchanged even with the ORACLE f_g).
        #
        # THE LEVEL is the whole difficulty, and it has ONE shape: ρ_ν = (1 − f̂_g)·M/E_r. So there is no
        # precedence to choose between — there are THREE INDEPENDENT estimators of the same f̂_g and they FUSE
        # by inverse variance, each contributing exactly the precision it earned and nothing where it has none:
        #   * the node's OWN belief — precision τ_own, i.e. its strand tilt or its own density deconvolution.
        #     Identically 0 on unstranded data and at any node with no factory, so it is INERT there rather
        #     than asserting ψ's uninformative reference as if it were a measurement (that route was measured
        #     at |Δw| = 0.628 on the unstranded arm);
        #   * the node ACROSS the seam, reframed — in practice the factory-solved intron, which is the one
        #     source that knows an RNA-free seam is RNA-free. Also gated by ITS τ_own, and charged the M5
        #     transfer variance for the hop (measured |Δw| = 0.294 on the unstranded arm);
        #   * the MASS identity closed with the MESSAGE's own gDNA claim (`residual_level`, M11) — the generic
        #     density deconvolution with the gDNA prior supplied by the neighbour instead of the intergenic
        #     pool. This is the one that exists EVERYWHERE: at `exon|exon` seams (97 % have no factory within
        #     reach and no strand when unstranded) and at every seam of a low-gDNA library. It is built on the
        #     good channel — gDNA transports across a hop at 0.96–0.99 while the RNA channel is the contaminant.
        # A pure PRECEDENCE over the first two, with a "no evidence ⇒ v_ν = ∞" third arm, silenced the RNA
        # channel outright in exactly the libraries where RNA is the entire signal (`gdna1`/`gdna5`/`none` ×
        # capture-ON, measured 0.0689 → 0.0935 on the worst of them). The fuse has no such arm, because M11 has
        # no such arm — and where two estimators disagree the fuse is what prices the disagreement, not a rank.
        #
        # The level is a λ-axis quantity — a density deconvolution delivers the gDNA-vs-RNA LEVEL and does not
        # assign the tilt — so M11's total is split across the strands by the MESSAGE's own tilt, and each
        # strand's share is formed against that strand's measured spliced flux. `v_μ` uses the spliced COUNT,
        # never the mass: the accumulator deposits fragments fractionally, so at a junction face the median
        # count is 33 against a median mass of 11 and the mass would over-state `v_μ` ~3×.
        # ── ⛔ THE DESTINATION-OWN PLUG-IN IS DELETED — the second BP violation ──────────────────────────
        # `_p_nu_own` fused the DESTINATION's own message-free self-solve into the level that shapes the
        # message arriving AT that destination. In BP a message m_{s→i}(x_i) must not contain φ_i, the
        # destination's own local evidence, because the belief is b_i ∝ φ_i · Π_j m_{j→i} — so φ_i would be
        # counted twice. Measured in the relay: 909/909 forward peel hops carry it. Removing it is FREE
        # (0 better / 0 worse / 32 flat, +0.0000 mwae) and marginally better on trust (boundary-single
        # z2|Q1 4.27 → 4.15), so there is no reason to keep an illegal term that buys nothing.
        # The level now fuses TWO legal estimators: the message's own claim, and M11's `residual_level`
        # (the destination's observed MASS closed against the message's gDNA claim — count-zero-information
        # legal, because the information is the imputed DENSITY and the count only converts it to a share).
        # ⚠ ``v_mu`` uses the spliced COUNT, never a mass — and with S5.e the two are the same array, so
        # the rule is now structural rather than a discipline: `mature_count` IS the junction fragment
        # count. Indexed by TRANSCRIPT strand, matching ``_mu``.
        _mu_s = (spl_p, spl_n)
        _v_mu_s = tuple(
            np.where(c > 0.0, 1.0 / np.maximum(c, _EPS), np.inf) for c in (SPL[:, 0], SPL[:, 1])
        )

        # ── ⭐ THE SEQUENTIAL SCAN'S OPERANDS, AS PYTHON LISTS (`*_l`) ──────────────────────────────────
        # `_relay` reads every one of these ONE ELEMENT AT A TIME — it is a Gauss-Seidel scan and cannot be
        # vectorised — and at genome scale that is ~6 M edge-iterations, ~40 reads each. `.tolist()` is
        # exact on float64 / int64 / bool (identical IEEE-754 doubles, identical ints), so this is
        # BIT-IDENTICAL by construction, and it buys ~3×: `lst[i]` costs a third of `arr[i]`, and it yields
        # a PYTHON float, whose arithmetic is ~3× `np.float64`'s because no intermediate is boxed in a 0-d
        # array. The array forms stay — the vectorised combine (`_transport`, `_peel_share`) needs them.
        (
            og_l,
            op_l,
            on_l,
            pg_own_l,
            pp_own_l,
            pn_own_l,
            mg_own_l,
            mp_own_l,
            mn_own_l,
            tau_own_l,
            M_l,
            E_g_l,
            E_r_l,
            n_node_l,
            logvar_l,
        ) = (
            np.asarray(a, np.float64).tolist()
            for a in (
                og,
                op,
                on,
                pg_own,
                pp_own,
                pn_own,
                mg_own,
                mp_own,
                mn_own,
                tau_own,
                M,
                E_g,
                E_r,
                _n_node,
                logvar_tot,
            )
        )
        ex_l, bnd_l, fp_l, fn_l = (a.tolist() for a in (ex_a, is_bnd_a, fp_a, fn_a))
        # the destination's own composition CERTAINTY — case (ii) of the mass pin's licence, below. Both
        # axes: an `intergenic|exon` EDGE is as structurally pure-gDNA as an intergenic NODE.
        g1_l = g1_locked(fp_a, fn_a).tolist()
        spl_p_l, spl_n_l = spl_p.tolist(), spl_n.tolist()
        SP_l, SN_l = SPL[:, 0].tolist(), SPL[:, 1].tolist()
        mu_l, v_mu_l = [c.tolist() for c in _mu_s], [c.tolist() for c in _v_mu_s]

        # ── ⛔ THE ACROSS-THE-SEAM (`_far`) LEVEL ESTIMATOR IS DELETED — it was a BP VIOLATION ──────────
        # It read the node ACROSS the seam and fused it into the peel share's level. That node is
        # **IDENTICALLY the source of the destination's OTHER message** — verified structurally on 35,421
        # peel edges, every condition — so its evidence reached the destination twice: once through its own
        # message and once inside this one, and the fuse added the two precisions as if independent. That is
        # the one term whose scope was a THIRD node; BP's rule is that a message into `i` may depend on
        # anything except `i`'s message back to its own source, and never on `i`'s OTHER neighbour.
        #
        # Deleting it is the cheapest CORRECT option and it is also a trust win. On a held-fixed node set
        # (HEAD's own confident quartile, so the arms are comparable): boundary-single `z2` **5.58 → 2.98**,
        # boundary-AMBIG 18.89 → 15.84, exon-AMBIG 64.39 → 52.39, ALL 8.98 → 8.28 — because the declared
        # variance there was inflated 2.81× (8.54× in the top decile) by the double count. Price: +0.0007
        # (refit 0) / +0.0009 (refit 1), 0 better / 2 worse / 30 flat.
        #
        # ⚠ The plan's proposed replacement — give the left message the RIGHT MESSAGE's delivered claim as
        # its far level — was prototyped and is **illegal AND worthless**: it does not remove the reuse (the
        # far node is still the other message's source, and its relayed belief carries strictly MORE of that
        # node's data than the raw self-solve did), and it recovers only 0.000117 of the 0.000696 price.
        #
        # ⚠ AND `_far` WAS NOT THE LAST THIRD-NODE DEPENDENCE. When the sweep ran TWO ρ-iterations, the second's
        # reframe faces are built from the destination's FUSED posterior — i.e. from the other message —
        # so the reframe itself still reuses it (measured median |Δlog ρ_face| between iterations 0.0116
        # stranded-capOFF to 0.1242 unstranded-capON, >1 % on 52.7–79.0 % of nodes). Deleting `_far` does
        # NOT make the message BP-legal; it removes the largest and most direct of two violations. The
        # remaining one is recorded in as P4b.
        #
        # The legal way to recover what `_far` supplied, if it is ever wanted back: carry the peel share as
        # a FUNCTION of the destination's state (a proper pairwise potential ψ(x_L, x_k)) instead of a
        # plugged-in point estimate, and let the destination's own ψ solve modulate it. No data reuse at
        # all. Structurally available, unimplemented, unmeasured.

        def _peel_share(tg, tpg, tp, tn):
            """The continuing share ``w`` and ``Var(log w)`` per strand, at every slot, for a message whose
            gDNA claim is ``(tg, tpg)`` and whose RNA claim is ``(tp, tn)``. Returns
            ``((w_p, vw_p), (w_n, vw_n))``; ``Var(log w) = +inf`` (⇒ zero precision, an inert message) only
            where NONE of the three estimators of the level exists.

            ⭐ The ``df`` face selector is gone: the mature flux a peel is measured against is the flux at
            THIS line, one number per transcript strand.

            ⚠ TWIN of `_peel_share_scalar` (the sequential relay's arm) — mirror any change into both."""
            _vg = np.where(np.asarray(tpg, np.float64) > 0.0, 1.0 / np.maximum(tpg, _EPS), np.inf)
            _nu_m, _vlog_m, _vl_m = residual_level(M, _n_node, tg, E_g, E_r, _vg)
            _A = np.asarray(tp, np.float64) + np.asarray(tn, np.float64)
            _a_p = np.where(_A > _EPS, np.asarray(tp, np.float64) / np.maximum(_A, _EPS), 0.0)
            out = []
            for _a, _mu, _vmu in (
                (_a_p, _mu_s[0], _v_mu_s[0]),
                (1.0 - _a_p, _mu_s[1], _v_mu_s[1]),
            ):
                # ── THE FUSE, in LINEAR density space (see `residual_level`'s return contract) ──────────
                # Each estimator is (ρ_i, Var_i); an own/far claim carrying a delta-method log-variance v is
                # the SAME statement as ρ ± ρ√v, so Var_i = ρ_i²·v_i. Linear is the coordinate that lets a
                # confident near-zero claim (the factory-solved intron across an RNA-free seam, measured
                # ρ_ν = 0.0006 ± 0.0012) actually pull the level down; a geometric fuse of positive modes
                # cannot reach zero, which is why "the intron sets the level" kept failing to.
                #
                # M11's level is a λ-axis TOTAL, split onto this strand by the message's own tilt. Mask the ∞
                # BEFORE the product — the standing ``0·inf = nan`` trap, np.where evaluating both branches.
                _nu_ms = _nu_m * _a
                _vl_s = np.where(np.isfinite(_vl_m), _vl_m, 0.0) * _a * _a
                _pm = np.where(
                    np.isfinite(_vl_m) & (_nu_ms > _EPS), 1.0 / np.maximum(_vl_s, _EPS), 0.0
                )
                _pt = _pm
                _live = _pt > _EPS
                _nu = np.where(
                    _live,
                    (_pm * _nu_ms) / np.maximum(_pt, _EPS),
                    0.0,
                )
                # back to the model's currency: the fused level's effective fragment COUNT k = ρ̂²/Var̂, and
                # its log-variance by the same exact rule M11 uses (ψ'(k) → 1/k for a well-determined level,
                # → π²/6 for one earned by a single fragment; never over-confident).
                # M11 already returned this level's log-variance; use it. The round trip it replaces
                # (`k = ρ̂²/V̂` then `ψ'(k)`) is EXACTLY `ψ'(1/v_log)`, and `_peel_share` was discarding
                # the `v_log` it re-derived. That mattered: `residual_level`'s `k ≥ 1` floor is an exact
                # limit of the TRUNCATION — k is the RNA share's effective fragment COUNT, and one
                # fragment is the least evidence a level can carry. Here `k` is a reciprocal variance,
                # not a count, so the same floor became a hard CEILING of `ψ'(1) = π²/6` on `_v_nu`: a
                # level M11 declared 10 nats uncertain was delivered as 1.64, over-stating its
                # confidence 6× on exactly the seams where the level is least determined. Measured:
                # accuracy 0 better / 0 worse / 32 flat at both refits, z2 better on every population.
                _v_nu = np.where(_live, _vlog_m, np.inf)
                _w = np.where(_live, peel_continue_share(_nu, _mu), 0.0)
                if _capture is not None:
                    _capture.setdefault(
                        "_lvl", []
                    ).append(  # inert: the level's provenance, per slot
                        {
                            "nu": np.asarray(_nu).copy(),
                            "v_nu": np.asarray(_v_nu).copy(),
                            "pm": np.asarray(_pm).copy(),
                            "nu_m": np.asarray(_nu_ms).copy(),
                            "mu": np.asarray(_mu).copy(),
                            "w": np.asarray(_w).copy(),
                            "v_g": np.asarray(_vg).copy(),
                            "vl_m": np.where(np.isfinite(_vl_m), _vl_s, np.inf),
                            "phi": np.asarray(tg * E_g / np.maximum(M, _EPS)).copy(),
                        }
                    )
                # a spliced DENSITY with no spliced COUNT cannot be priced (plan §4.7) ⇒ no claim at all.
                _ok = _live & (np.isfinite(_vmu) | ~(_mu > _EPS))
                out.append(
                    (
                        _w,
                        np.where(
                            _ok,
                            peel_share_logvar(
                                1.0 - _w,
                                np.where(_live, _v_nu, 0.0),
                                np.where(np.isfinite(_vmu), _vmu, 0.0),
                            ),
                            np.inf,
                        ),
                    )
                )
            return out

        def _peel_share_scalar(i, tg, tpg, tp, tn):
            """The SCALAR twin of `_peel_share`, for one slot ``i`` — see that docstring for the model.

            ⚠ TWIN: mirror any change into both. It exists because `_relay` is a sequential Gauss-Seidel
            scan (it cannot be vectorised) and calls this once per node per direction, so the array form
            runs ~50 numpy ops on 0-d arrays per call — 0.5–0.7 µs each against ~0.02 µs for the float
            expression. Same arithmetic in the same association order; every `np.where` becomes the branch
            it encodes, so the dead arms (`ζ(2,·)` and `residual_level`'s four `log_ndtr`s at a node with
            no level) are no longer evaluated and discarded. Measured 25× on this path, which is the bulk
            of a genome-scale calibration."""
            _vg = 1.0 / _fmax(tpg, _EPS) if tpg > 0.0 else math.inf
            _nu_m, _vlog_m, _vl_m = residual_level_scalar(
                M_l[i], n_node_l[i], tg, E_g_l[i], E_r_l[i], _vg
            )
            _fin = math.isfinite(_vl_m)
            _A = tp + tn
            _a_p = tp / _A if _A > _EPS else 0.0
            out = []
            for _a, _mu, _vmu in (
                (_a_p, mu_l[0][i], v_mu_l[0][i]),
                (1.0 - _a_p, mu_l[1][i], v_mu_l[1][i]),
            ):
                _nu_ms = _nu_m * _a
                _pm = (
                    1.0 / _fmax((_vl_m if _fin else 0.0) * _a * _a, _EPS)
                    if (_fin and _nu_ms > _EPS)
                    else 0.0
                )
                if not _pm > _EPS:  # no level ⇒ no claim: w = 0 at zero precision, an inert message
                    out.append((0.0, math.inf))
                    continue
                # ⚠ NOT `_nu_ms`: `(a·b)/a ≠ b` in floating point. This is the one-estimator remnant of
                # a three-estimator fuse (`_far` and the destination-own plug-in were deleted); kept
                # verbatim so this commit's delta is attributable to `_v_nu` alone.
                _nu = (_pm * _nu_ms) / _pm
                _v_nu = _vlog_m  # M11's own log-variance — see the twin
                _w = peel_continue_share_scalar(_nu, _mu)
                _wm = 1.0 - _w
                # a spliced DENSITY with no spliced COUNT cannot be priced (plan §4.7) ⇒ no claim at all.
                _ok = math.isfinite(_vmu) or not _mu > _EPS
                out.append(
                    (
                        _w,
                        _wm * _wm * (_v_nu + (_vmu if math.isfinite(_vmu) else 0.0))
                        if _ok
                        else math.inf,
                    )
                )
            return out

        # σ²_transfer per-hop damping: add the edge's transfer log-variance to the message's log-variance
        # (1/p → 1/p + σ²_transfer). This is now the DERIVED, direction-dependent M5 variance (0 on the matched
        # graft, Var(log r) on the peel/plain/anchor — computed per edge from ``logvar_tot`` below), retiring the
        # density-uniformity NPMLE proxy (`_mup`/`_vp`, no longer read).
        def _damp(p, s2t):
            return 1.0 / (1.0 / p + s2t) if p > 0.0 else 0.0

        def _damp_v(p, v):
            """STEP 3: add a log-variance to a precision. 1/p → 1/p + v ; v = ∞ ⇒ 0 (no claim)."""
            return p / (1.0 + p * v) if (p > 0.0 and math.isfinite(v)) else 0.0

        # ── ⚠ `_relay` AND `_transport` ARE DELIBERATE TWINS. DO NOT MERGE THEM. ────────────────────────
        # They implement the SAME per-edge transform in the same order — reframe `r` → detect the graft →
        # the composition-imputation LICENCE and the gDNA arm's own scale `r_g` → σ²_transfer (M5) → graft
        # the source's spliced → damp the three precision streams → graft precision ⊕ M8 ⊕ P1d → peel by
        # composition (M10) → strand filter → pin to M — and the correspondence is line-for-line. Merging
        # them into one polymorphic routine (`k` = a scalar index for the relay, `slice(None)` for the
        # combine) IS structurally possible; `_peel_share` already proves the pattern.
        #
        # ⛔ IT WAS MEASURED AND REJECTED. `_relay` is a SEQUENTIAL Gauss-Seidel scan — it reads
        # `rg[s]`/`pg[s]`/… at a source that an earlier iteration of the same pass may already have written
        # (verified: every one of `rg pg pp pn mg mp mn tau` is both read at `s` and written at `i`), so it
        # cannot be vectorised and its per-node arithmetic is plain Python floats. Routing that through the
        # combine's numpy form costs **15.7×** per operation (0.083 µs → 1.307 µs measured), i.e. **+5.5 s per
        # calibration** on a 74,494-node chain — to save ~60 lines. Bad trade.
        #
        # ⚠ SO THE DUPLICATION IS REAL AND MUST BE MAINTAINED BY HAND. A change to ONE of these must be
        # mirrored in the other. The two also differ deliberately in three edge cases, and those differences
        # are NOT bugs to unify: (1) the relay skips invalid edges with `continue` while the combine masks
        # them to `r = 0`; (2) `_damp` uses the raw `p` where `_dv` uses `max(p, _EPS)`, which differ only for
        # `0 < p < _EPS`; (3) the relay short-circuits the graft block under `if _gr:` while the combine
        # evaluates `graft_frame_logvar(r)` on every edge and masks. Any merge must reconcile all three.
        # The relay's side of the pair is now scalar-native throughout — Python-float operands (the `*_l`
        # block) calling the `*_scalar` twins of the shared primitives — which is what makes it fast, and
        # one more reason the two forms cannot collapse into one.
        def _relay(seq, nbr, rho, pop):
            # every operand here is a Python float or bool — see the `*_l` block above.
            # ⭐ ``dst_face``/``src_face``/``df``/``sf`` are gone: one ρ_tot and one mature flux per slot,
            # so the only thing that still varies between the forward and backward passes is ``nbr``.
            rg, rp, rn = og_l.copy(), op_l.copy(), on_l.copy()
            pg, pp, pn = pg_own_l.copy(), pp_own_l.copy(), pn_own_l.copy()  # full → MODE fusion
            mg, mp, mn = (
                mg_own_l.copy(),
                mp_own_l.copy(),
                mn_own_l.copy(),
            )  # MEASUREMENT (anchor + spliced)
            tau = tau_own_l.copy()  # COMPOSITION (τ_λ) → the λ-message
            for i in seq:
                s = nbr[i]
                if s < 0:
                    continue
                # The reframe is a PAIR quantity: the source's ρ_tot on the face it emits FROM (``src_face[s]`` —
                # WITH the spliced when that face is the acceptor, §6 component-set matching), and the dst's on
                # the face it receives ON. Using the mature-free NODE ρ_tot as the denominator (the earlier
                # form) left an uncancelled ``1+ρ_spl/ρ_unspl`` at every acceptor boundary, and because the relay
                # re-scales the SAME running density each hop that factor COMPOUNDS — measured Σf_c ≈ 71 at
                # introns (a composition must sum to 1) and r up to 10³ into exons.
                rho_src = rho[s]
                r = (
                    (rho[i] / rho_src) if (rho_src > _EPS and rho[i] > _EPS) else 1.0
                )  # no frame ⇒ pass-through
                # GRAFT (boundary → EXON, §6): the boundary's measured mature is a density AT THE SOURCE, so it
                # joins the source's RNA BEFORE the reframe; the peel is measured at the destination and so is
                # applied after. Both use the face that FACES the other endpoint (``sf``/``df``), never the
                # node-pooled sum. Only an EXON receives the graft — an intron carries no mature (`ex_a[i]`,
                # not `is_reg_a[i]`, which grafted the junction's whole mature flux into every flanking intron).
                _gr = ex_l[i] and bnd_l[s]
                # σ²_transfer = Var(log r) (M5): 0 on the matched-set GRAFT (r is common-mode across {g,R} and
                # cancels in the composition — charging it there is a double-count), Var(log r) elsewhere (peel /
                # plain reframe / partial-anchor, where r is load-bearing). The COMPOSITION-mismatch term b̂² is
                # the combine's job (`_transport`): the relay has no destination self-solve to measure a gap
                # against — its running belief is already fused with the messages, so a DL gap here would be
                # feedback, not evidence.
                s2t = 0.0 if _gr else (logvar_l[i] + logvar_l[s])
                gp = spl_p_l[s] if _gr else 0.0
                gn = spl_n_l[s] if _gr else 0.0
                # ⭐⭐ THE gDNA SCALE (the scalar twin of the combine's `r_g` — see the derivation above the
                # relay). ``_lend``: may this source lend a COMPOSITION? Two conditions, and they are
                # different questions about the same step:
                #   * SUPPLY — did the source state both components of the pair from its OWN crossing
                #     population? "Supplied" is a statement about PRECISION, not about a density's value.
                #   * POPULATION — is the source measuring the same RNA population as the destination?
                #     `_pop` is that test for this step (derived where `_pop_l_a` is built).
                # Where either fails, the gDNA LEVEL crosses UNSCALED and the mass pin is off.
                _lend = pop[i] and pg[s] > 0.0 and (pp[s] + pn[s]) > 0.0
                r_g = r if _lend else 1.0
                tg, tp, tn = rg[s] * r_g, (rp[s] + gp) * r, (rn[s] + gn) * r
                tpg, tpp, tpn = (
                    _damp(pg[s], s2t),
                    _damp(pp[s], s2t),
                    _damp(pn[s], s2t),
                )  # full (mode)
                tmg, tmp, tmn = (
                    _damp(mg[s], s2t),
                    _damp(mp[s], s2t),
                    _damp(mn[s], s2t),
                )  # measurement
                ttau = _damp(tau[s], s2t)  # composition
                # The grafted mature is a MEASUREMENT (a junction COUNT), not an imputation, so it carries its
                # own precision and is NOT τ-gated — the source's PREDICTION precision is 0 on unstranded data
                # and would otherwise drop the graft on the floor. It enters BOTH the mode-fusion (full) and the
                # MEASUREMENT stream (never the composition τ — a count is not a composition vote).
                if _gr:
                    # M8 (`graft_frame_logvar`): the grafted spliced density is measured in the DESTINATION
                    # exon's frame (its fragments' blocks lie in the exons), so it has no matched gDNA partner
                    # to cancel ``r`` against and M5's graft-zero does not cover it. Charge the frame step it
                    # is implicitly mis-lifted by. Identically 0 at r = 1.
                    _s2f = s2t + graft_frame_logvar_scalar(r)
                    _sps = SP_l[s]
                    _spc = _sps / (1.0 + _sps * _s2f) if _sps > _EPS else 0.0
                    _sns = SN_l[s]
                    _snc = _sns / (1.0 + _sns * _s2f) if _sns > _EPS else 0.0
                    tpp += _spc
                    tpn += _snc
                    tmp += _spc
                    tmn += _snc
                    _vgp, _vgn = vgp_l[i], vgn_l[i]
                    tpp, tmp = _damp_v(tpp, _vgp), _damp_v(tmp, _vgp)
                    tpn, tmn = _damp_v(tpn, _vgn), _damp_v(tmn, _vgn)

                if (
                    bnd_l[i] and ex_l[s]
                ):  # EXON → boundary: PEEL by COMPOSITION (scale by the share)
                    (_wp, _vwp), (_wn, _vwn) = _peel_share_scalar(i, tg, tpg, tp, tn)
                    tp, tn = tp * _wp, tn * _wn
                    tpp, tmp = _damp_v(tpp, _vwp), _damp_v(tmp, _vwp)
                    tpn, tmn = _damp_v(tpn, _vwn), _damp_v(tmn, _vwn)
                if not fp_l[i]:
                    tp, tpp, tmp = 0.0, 0.0, 0.0
                if not fn_l[i]:
                    tn, tpn, tmn = 0.0, 0.0, 0.0
                # ── PIN THE CONTEXT TO THIS SLOT'S OBSERVED MASS — WHERE NO BELIEF ENTERS THE BUDGET ──────────
                # `Σ_c ρ_c·E_c = M` is an IDENTITY under the imputation premise, not an approximation: a matched
                # reframe delivers ρ_c^msg = a_c·ρ_tot(dst) = ρ_c^dst,true, so the components account for exactly
                # the fragments the slot observed. The pin restores it with a common factor `k = M/S`, and the
                # budget `S` fills every component the context does NOT supply from the destination's own
                # density — which is what keeps a partial claim partial, and also what can make the delivered
                # value a function of the destination's own BELIEF. `TRAPS.md` D4 permits a message to use the
                # destination's CONSTANTS and its OBSERVATIONS, never its beliefs. So the pin is licensed in
                # exactly the two states where no belief reaches `S`:
                #
                #   (i)  `_lend` — the context SUPPLIED the composition, so the premise is granted and the
                #        identity is implied. This is the reframe's own predicate: one licence, one place.
                #   (ii) `g1_l[i]` — the destination is a structurally pure-gDNA object, so there is no
                #        unsupplied component to fill in. `f_g = 1` there is STRUCTURE, not a belief, and
                #        `S = ρ_g·E_g` makes `k = (M/E_g)/L_in`: the pin hands the object its OWN MEASURED
                #        density. ⭐ That is the operator the capture landscape travels on — an
                #        `intergenic|exon` EDGE measures the gDNA density at its own capture stratum, and the
                #        exon behind it has no other way to hear it (a per-capture-class landscape ratio built
                #        to do the same job explicitly measured inert; `node_geometry`'s deleted-landscape note).
                #        ⚠ Gated on `g1_locked`, which spans BOTH axes, and NOT on `node_init`'s node-only
                #        `struct_lock` — the object in question is an EDGE. `g1_locked`'s docstring separates the
                #        two questions; this one is about the destination's own certainty.
                #
                # ⛔ ANYWHERE ELSE IT WAS D4 AT FULL STRENGTH. The closed form is `k = 1/(φ_msg + R_own)`, a
                # saturating map whose fixed point is `(1−R_own)·ρ_tot` with `R_own` the RNA share of the
                # destination's OWN self-solve — EXACTLY ½ at a slot with no composition evidence. So it drove
                # the delivered gDNA FRACTION to ½ regardless of the truth and regardless of what the source
                # measured. Measured on `toy_harness.py --spec nested_exons` (uniform gDNA field, no intron, five
                # evidence-free exons): the delivered level rose 0.071 → 0.102 → 0.097 → 0.162 → 0.192 → 0.215
                # with step distance from the gene ends, 2.8× the truth at the deepest slot, and the gene's
                # mass-weighted |Δf_g| was 0.2618 against 0.0760 licensed. ⚠ It hid because the running product
                # of the rescales TELESCOPES back to 1 at the far end of a gene, so no aggregate, endpoint or
                # conservation check could see it.
                # ⛔ AND IT IS NOT A DELETION: unanchored, the relay's residual compounded to `Σ_c ρ_c E_c / M`
                # p99 = 31–288× and max 519×, and rescaling all three components blindly instead regressed
                # capture-OFF 3.6×. Case (ii) is not a caveat on case (i) either — dropping it delivers the
                # OFF-PROBE floor to every exon under capture (fires `test_gdna_scale_rule`'s capture gate at
                # 20× and 200×, measured).
                # ⚠ THE COMBINE NEEDS NO MIRROR OF THIS, and the DO-NOT-MERGE note above should not send the
                # next reader looking for one: `_transport` leaves the delivered `tg/tp/tn` exactly as measured
                # and uses `_pin_v` only to build the comparison frame for the DL mismatch test.
                if _lend or g1_l[i]:
                    _sg = tg if tpg > 0.0 else og_l[i]
                    _sp = tp if tpp > 0.0 else op_l[i]
                    _sn = tn if tpn > 0.0 else on_l[i]
                    _sv = _sg * E_g_l[i] + (_sp + _sn) * E_r_l[i]
                    if _sv > _EPS and M_l[i] > _EPS:
                        _k = M_l[i] / _sv
                        tg, tp, tn = tg * _k, tp * _k, tn * _k
                rg[i], pg[i] = _fuse(og_l[i], pg_own_l[i], tg, tpg)
                rp[i], pp[i] = _fuse(op_l[i], pp_own_l[i], tp, tpp)
                rn[i], pn[i] = _fuse(on_l[i], pn_own_l[i], tn, tpn)
                # measurement + composition precisions are INDEPENDENT evidence → additive (inverse-variance) fuse
                mg[i] = mg_own_l[i] + tmg
                mp[i] = mp_own_l[i] + tmp
                mn[i] = mn_own_l[i] + tmn
                tau[i] = tau_own_l[i] + ttau
            # back to arrays for the vectorised combine — exact, they are the same doubles
            return tuple(
                np.asarray(a, np.float64) for a in (rg, rp, rn, pg, pp, pn, mg, mp, mn, tau)
            )

        def _seam_pair(rho):
            """Per strand: the graft's premise log-variance — ONE library-level scalar, fitted by method of
            moments from the destination-frame disagreement of exons' flanking seam PAIRS
            (`graft_premise_logvar`) and applied to every graft edge.

            ⚠⚠ **A DEBT, not a model.** The one scalar stands in for a quantity that splits ≥30× on whether
            the boundary carries a transcript TERMINUS — a bit the region map does not have. Re-derive this
            per structural class when TSS/TES land (P1g). See `graft_premise_logvar`."""
            out = []
            for spf, vmu in ((spl_p, 0), (spl_n, 1)):
                fl, fr = _flank_dom(rho, spf)
                # each seam's own noise: its spliced COUNT (never the mass) ⊕ its lift's scale sampling
                # (M5's source leg; the destination's leg is common to both lifts and cancels in ``d``).
                _lv = np.where(np.isfinite(logvar_tot), logvar_tot, 0.0)
                per, pooled = graft_premise_logvar(
                    fl,
                    fr,
                    np.where(_vl_a, _v_mu_s[vmu][_sl_a] + _lv[_sl_a], np.inf),
                    np.where(_vr_a, _v_mu_s[vmu][_sr_a] + _lv[_sr_a], np.inf),
                )
                # ⭐ The POOLED fit is applied to EVERY graft edge; the per-edge value is NOT used. Two
                # reasons, and the statistical one is decisive:
                #
                #  * ``d²`` from ONE pair is a single draw of a χ²₁ scaled by the true variance, so its own
                #    coefficient of variation is √2 — a per-edge "measurement" of a variance is mostly noise,
                #    and it both over- and UNDER-charges around the right mean. The under-charging half is
                #    what does the damage, because it REPLACES the population value on the ~48 % of edges
                #    where it fires. Measured: per-edge+pooled 870,245 confidently-wrong reads / z2 11.12,
                #    pooled alone **762,000 / 8.98**, and exon-single z2 5.68 → **2.46**. Shrinking all the
                #    way to the pooled mean is simply the better estimator.
                #  * it also removes a real BP objection (owner, 2026-07-26): with the per-edge form, the
                #    message from the LEFT seam carried a variance computed from the RIGHT seam's counts, so
                #    a non-adjacent node's data reached the destination twice. Now no message's precision
                #    depends on anything but its own edge and one library-level constant — the same standing
                #    as ``κ`` and both strand overdispersions.
                if (
                    _capture is not None
                ):  # inert: the fitted scalar and the population it was fitted on
                    _ok = (fl > _EPS) & (fr > _EPS)
                    _d = np.log(np.maximum(fl, _EPS)) - np.log(np.maximum(fr, _EPS))
                    _vv = np.where(_vl_a, _v_mu_s[vmu][_sl_a] + _lv[_sl_a], 0.0) + np.where(
                        _vr_a, _v_mu_s[vmu][_sr_a] + _lv[_sr_a], 0.0
                    )
                    _capture.setdefault("_glv", []).append(
                        {
                            "strand": vmu,
                            "omega": pooled,
                            "n_pairs": int(_ok.sum()),
                            "ok": _ok.copy(),
                            "d": _d.copy(),
                            "noise": _vv.copy(),
                            "Ed2": float((_d[_ok] ** 2).mean()) if _ok.any() else 0.0,
                            "Enoise": float(_vv[_ok].mean()) if _ok.any() else 0.0,
                        }
                    )
                out.append(np.full_like(per, pooled))
            return out[0], out[1]

        # the relay runs on the INPUT-belief frame, so its seam pair is formed from that
        vgp_prem, vgn_prem = _seam_pair(rho0)
        vgp_l, vgn_l = vgp_prem.tolist(), vgn_prem.tolist()
        left_l, right_l = left.tolist(), right.tolist()
        rho0_l = rho0.tolist()
        _pop_l_l, _pop_r_l = _pop_l_a.tolist(), _pop_r_a.tolist()
        fwd = _relay(order_list, left_l, rho0_l, _pop_l_l)
        bwd = _relay(order_list[::-1], right_l, rho0_l, _pop_r_l)
        # ── the COMBINE: transport α (from left neighbour) + β (from right neighbour) into the node's frame with
        # the LAZY ρ_tot (two-iteration — the 2nd uses the both-message composition), fuse, ÷M_dst → the ψ solve.
        li, ri, vl, vr, sl, sr = _li_a, _ri_a, _vl_a, _vr_a, _sl_a, _sr_a

        # The VECTORISED twin of `_relay` — see the DO-NOT-MERGE note there, which applies to both.
        def _transport(src, valid, fwd_arrs, rho, pop):
            rg, rp, rn, pg, pp, pn, mg, mp, mn, tau = fwd_arrs
            # A slot with no frame (no mass ⇒ no ρ_tot, §5) cannot reframe: the message passes through at
            # r=1. Falling back to ``rho_src = 1.0`` instead made r the destination's ABSOLUTE density (10³
            # on a short node) — a raw scale masquerading as a ratio. The relay already guards this way.
            framed = valid & (rho[src] > _EPS) & (rho > _EPS)
            r = np.where(framed, rho / np.maximum(rho[src], _EPS), np.where(valid, 1.0, 0.0))
            # GRAFT before the reframe (a density measured AT the source); PEEL after (measured at the dst).
            # The graft only into an EXON — see the relay's twin.
            graft = ex_a & is_bnd_a[src] & valid
            gp = np.where(graft, spl_p[src], 0.0)
            gn = np.where(graft, spl_n[src], 0.0)
            # ⭐⭐ THE gDNA SCALE — see the derivation above the relay, and the relay's own twin. ``lend``
            # asks two things of the step: the λ-emission gate's predicate of the SOURCE — it may lend a
            # composition only if it SUPPLIED both components of the pair — and, of the PAIR, whether the
            # two objects measure the same RNA POPULATION (``pop``; derived where ``_pop_l_a`` is built, and
            # a property of the (EDGE, side) pair rather than of either object). Where either fails, the
            # reframe is a false premise and the gDNA LEVEL crosses UNSCALED — gDNA is uniform before
            # capture, and capture is carried by the pure-gDNA objects' own measurements, not by a scale.
            #
            # ⚠ **A GRAFT edge does NOT license it, and that is a deliberate divergence from the λ gate**,
            # which counts the grafted junction precision as RNA supplied. `TRAPS.md` F9: mature RNA does
            # not cross an intron↔exon EDGE contiguously, so that EDGE's OWN spanning population is
            # gDNA + nascent and the junction flux is a measurement of RNA that lives in the DESTINATION —
            # the routing operator that exists precisely because that component cannot cross by imputation.
            # Using it to license the imputation would be circular. λ is a claim about the pair; the
            # reframe is a claim about the source's own crowding. ⚠ M5's graft-zero rests on the premise
            # this denies; it is a VARIANCE, separately measured and landed, and is left alone.
            #
            # ⚠ **The gDNA conjunct is INERT in practice and is kept for faithfulness, not for effect** —
            # recorded because a perturbation dropping it fires no gate, and that is a fact about the
            # solver rather than a hole. ``pg[src] == 0`` with RNA precision live requires the source's own
            # gDNA density to be 0, and then ``rg[src]·r_g`` is 0 whatever the scale; and on any chain with
            # a structural gDNA anchor (every real one) the relay propagates gDNA precision to every slot,
            # so the state does not arise. It stays because the predicate is a statement about the PAIR.
            lend = pop & (pg[src] > 0.0) & ((pp[src] + pn[src]) > 0.0)
            r_g = np.where(lend, r, np.where(valid, 1.0, 0.0))
            # ⛔ THE SHARE TRANSFER WAS IMPLEMENTED HERE AND REVERTED, 2026-07-27.
            # It delivers the source's own composition share carried onto the destination's scale —
            # `f̂_c = ctx_c·E_c[src]/S_src`, then `t_c = f̂_c·M/E_c` — which is BP-clean, keeps a partial
            # claim partial with the deficit measured AT THE SOURCE, and needs no reframe at all.
            # It hit every target it was aimed at (unstranded × capture-OFF × gDNA-bearing −0.0022/−0.0020,
            # `gdna_none` untouched) **and lost badly everywhere else**: suite +0.0119 (r0) / +0.0121 (r1),
            # capture-ON **+0.0264/+0.0274**, unstranded × capON **+0.0378/+0.0415**, corr −0.033.
            # ⭐ THE LESSON, and it corrects the derivation: **"composition transfers" is a WEAKER premise
            # than "density transfers, reframed".** Under capture an exon and its flanking boundary do NOT
            # share a composition — the boundary sits on the capture slope (0.125× the exon, 2113× the
            # intron at verystrong) — so discarding `r` in favour of a mass ratio throws away the
            # enrichment information the reframe carries. `r` only *looked* inert while `_pin_v` was
            # cancelling it; once the pin no longer rewrites the delivered density, `r` is load-bearing.
            tg, tp, tn = rg[src] * r_g, (rp[src] + gp) * r, (rn[src] + gn) * r
            # σ²_transfer = Var(log r) (M5, the tested pure law): 0 on the matched-set graft (r common-mode ⇒
            # cancels — a double-count otherwise), Var(log r) = logvar_tot[dst]+logvar_tot[src] elsewhere (peel /
            # plain reframe / partial-anchor — r load-bearing). This is the SCALE half of the cliff cost; the
            # COMPOSITION half is the DL b̂² applied at the end of this function. See the relay's twin.
            s2t = transfer_logvar(logvar_tot, logvar_tot[src], graft)

            def _dv(p, s2=s2t):
                return np.where(
                    valid & (p[src] > 0.0), 1.0 / (1.0 / np.maximum(p[src], _EPS) + s2), 0.0
                )

            tpg, tpp, tpn = _dv(pg), _dv(pp), _dv(pn)  # full → mode fusion
            tmg, tmp, tmn = _dv(mg), _dv(mp), _dv(mn)  # measurement (anchor gDNA + spliced RNA)
            ttau = _dv(tau, s2t)  # composition (τ) → the λ-message
            # the graft's MEASUREMENT precision — never τ-gated (see the relay's twin). ``_sp``>0 only on a GRAFT
            # edge, where s2t≡0, so the inf→0 substitution below touches only already-masked entries (a
            # zero-count node has logvar_tot=+inf ⇒ s2t=inf; ``0·inf`` would nan the masked branch np.where evals).
            _sp, _sn = np.where(graft, SPL[:, 0][src], 0.0), np.where(graft, SPL[:, 1][src], 0.0)
            _s2t_spl = np.where(np.isfinite(s2t), s2t, 0.0)
            # M8 — the graft's FRAME-MISLIFT variance (see the relay's twin and `graft_frame_logvar`): the
            # measured spliced already sits in the destination exon's frame, so ``r`` is NOT common-mode for it
            # and M5's graft-zero does not cover it. 0 where r = 1, so it is inert without a capture step.
            _s2t_spl = _s2t_spl + np.where(graft, graft_frame_logvar(r), 0.0)
            _spc = np.where(_sp > _EPS, _sp / (1.0 + _sp * _s2t_spl), 0.0)
            _snc = np.where(_sn > _EPS, _sn / (1.0 + _sn * _s2t_spl), 0.0)
            tpp, tpn = tpp + _spc, tpn + _snc  # into the mode-fusion precision …
            tmp, tmn = (
                tmp + _spc,
                tmn + _snc,
            )  # … and the measurement stream (a count, never composition τ)
            # ⭐ P1d — the graft's PREMISE variance (`graft_premise_logvar`), applied to the WHOLE RNA
            # claim after the spliced arm is folded in, because the premise is about the SUM: measured
            # FLAT in the spliced share w_μ (Var 2.02 → 1.83 across w_μ 0.47 → 1.00, while Var/w_μ²
            # swings 5×), so charging the spliced arm alone would reach only 10–93 % of the delivered
            # confidence while the error contaminates 63–95 % of the delivered density.
            _vgp = np.where(graft, vgp_prem, 0.0)
            _vgn = np.where(graft, vgn_prem, 0.0)
            tpp, tmp = tpp / (1.0 + tpp * _vgp), tmp / (1.0 + tmp * _vgp)
            tpn, tmn = tpn / (1.0 + tpn * _vgn), tmn / (1.0 + tmn * _vgn)

            peel = (
                is_bnd_a & ex_a[src] & valid
            )  # EXON → boundary: PEEL by COMPOSITION (the relay's twin)
            (_wp, _vwp), (_wn, _vwn) = _peel_share(tg, tpg, tp, tn)
            tp = np.where(peel, tp * _wp, tp)
            tn = np.where(peel, tn * _wn, tn)

            def _dv_arr(pr, vv):
                _f = np.isfinite(vv)
                return np.where(
                    peel, np.where(_f, pr / (1.0 + pr * np.where(_f, vv, 0.0)), 0.0), pr
                )

            tpp, tmp = _dv_arr(tpp, _vwp), _dv_arr(tmp, _vwp)
            tpn, tmn = _dv_arr(tpn, _vwn), _dv_arr(tmn, _vwn)
            tp, tpp, tmp = (
                np.where(fp_a, tp, 0.0),
                np.where(fp_a, tpp, 0.0),
                np.where(fp_a, tmp, 0.0),
            )
            tn, tpn, tmn = (
                np.where(fn_a, tn, 0.0),
                np.where(fn_a, tpn, 0.0),
                np.where(fn_a, tmn, 0.0),
            )
            if _capture is not None:  # inert: the PRE-PIN state, for the weighted-rescale prototype
                _capture.setdefault("_pin", []).append(
                    {
                        "src": np.asarray(src).copy(),
                        "valid": np.asarray(valid).copy(),
                        # the gDNA scale rule, per hop: the composition-imputation licence and the two
                        # scales, so an instrument can read the rule off a real run rather than infer it.
                        "lend": np.asarray(lend).copy(),
                        "r": r.copy(),
                        "r_g": r_g.copy(),
                        "tg": tg.copy(),
                        "tp": tp.copy(),
                        "tn": tn.copy(),
                        "tpg": tpg.copy(),
                        "tpp": tpp.copy(),
                        "tpn": tpn.copy(),
                        # the COMMON-mode variance every component shares: the reframe's scale (M5) plus the
                        # source's own count. Everything else in a component's variance is its own.
                        "s2t": np.where(np.isfinite(s2t), s2t, 0.0),
                        "n_src": np.asarray(_n_node)[src].copy(),
                        # the graft: how much of the RNA claim is a MEASUREMENT rather than an imputation
                        "spl_p": _sp.copy(),
                        "spl_n": _sn.copy(),
                        "spl_prec": (_spc + _snc).copy(),
                        "graft": np.asarray(graft).copy(),
                    }
                )
            # ── P1e: the conservation SURPRISE as a DerSimonian–Laird damping term ────────────────────
            # The claim asserts S = Σ_c ρ_c·E_c fragments; the node observed M. That is an IDENTITY, and
            # `_pin_v` restores it by fiat and DISCARDS the residual. Price the residual instead:
            #     δ = log(M/S),  Σ = σ_cm²·11ᵀ + diag(w),  αᵀΣα = Σ_c α_c s_c   (M12's own error model)
            #     b̂²_cons = max(0, δ² − αᵀΣα − 1/n_dst)                          (the DL between-study form)
            # and attribute it by the conditional mean of the error given the observed violation,
            # E[ε | αᵀε = −δ] = −δ·s/(αᵀΣα), i.e. the rank-1 inflation Σ += b̂²·s sᵀ/(αᵀΣα)², whose
            # diagonal is Δv_c = b̂²·(s_c/αᵀΣα)². The scale multiplies b̂² (1.0 = the derived term).
            _p3 = np.stack([tpg, tpp, tpn], axis=-1)
            _sup = _p3 > 0.0
            _mc = np.where(_sup, np.stack([tg, tp, tn], axis=-1), np.stack([og, op, on], axis=-1))
            _mc = _mc * np.stack([E_g, E_r, E_r], axis=-1)
            _S = _mc.sum(axis=-1)
            _okc = valid & (_S > _EPS) & (M > _EPS)
            _al = _mc / np.maximum(_S, _EPS)[..., None]
            _vc = np.where(_sup, 1.0 / np.maximum(_p3, _EPS), 0.0)
            _s2c = (np.where(np.isfinite(s2t), s2t, 0.0) + 1.0 / np.maximum(_n_node[src], _EPS))[
                ..., None
            ]
            _sv = np.where(_sup, _s2c + _al * np.maximum(_vc - _s2c, 0.0), 0.0)
            _aSa = np.sum(_al * _sv, axis=-1)
            _dlt = np.where(_okc, np.log(np.maximum(M, _EPS) / np.maximum(_S, _EPS)), 0.0)
            _den = _aSa + 1.0 / np.maximum(_n_node, _EPS)
            _b2 = np.maximum(_dlt * _dlt - _den, 0.0)
            # ⚠⚠ **PARTLY A DEBT — THIS PRICES A BIAS AS A VARIANCE.** On a large share of its firing mass
            # ``δ`` is systematic (``E[δ]`` ≈ −0.5 to −1.5; bias share 53–77 % on graft × one-component,
            # 98.9–99.2 % at intergenic destinations), and a variance cannot move a mode toward truth. It is
            # landed because it is the only change measured to improve ACCURACY and honest PRECISION together,
            # not because the bias half is derived. It was hoped those strata were inert (intergenic is
            # ``solvable=False``) — **measured and REFUTED: 90–100 % of the damping mass lands on solvable
            # destinations.** The magnitude is also not what works: a flat pooled constant beats the derived
            # ``b̂²`` on 3 of 4 conditions and ``b̂² := δ²`` is identical, while permuting ``b̂²`` FAILS — so
            # ``δ`` selects WHICH message to distrust and the calibration adds nothing (ω_graft's shape again).
            # **When the bias strata are diagnosed, this term must SHRINK.**

            # ── ⭐ THE SCOPE: only the OVER-claim direction is evidence against the MESSAGE ─────────────
            # `S` is a COMPLETE budget: `_pin_v`'s partial-claim semantics fill every component the
            # message does not supply from the node's OWN density. So a shortfall (δ > 0, "your
            # components account for fewer fragments than you observed") can be the node's own density
            # being too low just as easily as the message being wrong — it does not attribute. The
            # OVER-claim direction does: every density is non-negative, so nothing the unsupplied
            # components could be would rescue a budget that already exceeds `M`. Charging only that
            # half is the licensed statement, and it is where the term is a contradiction of a hard
            # observable rather than a disagreement with a guess.
            # `RIGEL_P1E_SCOPE`: `neg` (default) = δ < 0 · `sup` = the strictly harder test, the
            # SUPPLIED components alone over-claim · `all` = unscoped (both directions).
            _b2 = np.where(_dlt < 0.0, _b2, 0.0)
            # ⭐ THE COMMON DIRECTION (the derived law). Because ``αᵀ1 ≡ 1`` — α is a share vector over the
            # same budget S — adding the SCALAR b̂² to every supplied component's log-variance satisfies
            # the constraint ``αᵀΣ'α = αᵀΣα + b̂²`` exactly, and it leaves ``Var(λ)`` **identically
            # unchanged**: a common shift of both arms cannot move the split. That is the whole reason to
            # prefer it. The rank-1 form borrows the CONDITIONAL MEAN's direction and applies it as a
            # variance inflation, which MC shows over-damps λ 5× when the true error is a pure scale
            # error (λ z² 1.00 → 0.21) while still leaving the RNA arm over-confident (z² 2.88).
            # A scalar cannot identify a direction: one observation, one parameter. The rank-1 variant was
            # implemented, MC-refuted and DELETED; do not rebuild it.
            _dg = _dp = _dn = np.where(_sup.any(axis=-1), _b2, 0.0)
            tpg, tmg = tpg / (1.0 + tpg * _dg), tmg / (1.0 + tmg * _dg)
            tpp, tmp = tpp / (1.0 + tpp * _dp), tmp / (1.0 + tmp * _dp)
            tpn, tmn = tpn / (1.0 + tpn * _dn), tmn / (1.0 + tmn * _dn)
            # The mass frame — for the mismatch COMPARISON only. The delivered densities `tg/tp/tn` are left
            # exactly as the source measured and the reframe delivered them: a component's LEVEL is an
            # absolute rate, and re-normalising it against a budget built from the destination's own belief
            # is what made the message self-confirming (`_pin_v`). Note the message
            # packet's other two claims are unaffected either way — `tlam` and `tth` are scale-free, so the
            # pin's common factor cancels from them identically.
            pin_g, pin_p, pin_n = _pin_v(tg, tp, tn, tpg, tpp, tpn)
            # ── the COMPOSITION half of the cliff cost: the DL mismatch deflation, in the PINNED frame.
            # Both sides then account for the same total, so the common scale (the reframe residual) is gone
            # from G and only the share drift is left. Every stream is deflated — measured: the pin recovers
            # only when the composition τ-stream is damped alongside the measurement one.
            # ── THE λ-EMISSION GATE (structural, and PRIOR to any damping question) ────────────────────────
            # A composition message is a claim about the SPLIT, ``λ = log(f_g/f_R)``. A source that carries only
            # ONE component has no such claim to make — λ is not "large" for it, it is UNDEFINED. The canonical
            # case is exactly the geometry ``intergenic | seam | EXON``: RNA cannot cross a gene boundary, so the
            # seam is structurally RNA-free, while the exon it feeds has RNA. The combine builds the λ mode as
            # ``mo_g − mo_R`` with ``mo_R`` floored at ``log(_EPS)``, so a zero-RNA message silently becomes the
            # maximally confident assertion "this node is 100 % gDNA" — a numerical artifact of the floor — while
            # its precision is real, having been accumulated along the relay by nodes that never said that.
            # Measured (`scratchpad/seam_lambda_audit.py`) on unstranded capture-ON: 90 nodes received that
            # artifact and every one of them was solved to f_g = 1.000 against a mean oracle of 0.814.
            # Gating it here, at EMISSION, is the honest fix and it is structural — the same kind of presence
            # test as the ``fp_a``/``fn_a`` strand gates, no threshold. It must NOT be left to the DL term (which
            # only reaches it where the destination has its own evidence, i.e. never on unstranded data) and it
            # is NOT a loss of information: a pure-gDNA source's authority is a DENSITY LEVEL, and that travels
            # on the measurement stream (``tmg``), which is precisely the three-stream separation.
            # The gate asks whether the source SUPPLIED both components of the λ pair. "Supplied" is a
            # statement about PRECISION, not about the density's value: a component carrying zero precision is
            # not supplied however large its density, and a component at zero density with live precision IS
            # supplied — it is the claim "there is none of this here", which is exactly a composition claim.
            # Testing the density conflates the two and silences λ wherever a legitimate zero is emitted.
            ttau = np.where((tpg > 0.0) & ((tpp + tpn) > 0.0), ttau, 0.0)
            g_g, c_g = mismatch_gap(pin_g, og)
            g_p, c_p = mismatch_gap(pin_p, op)
            g_n, c_n = mismatch_gap(pin_n, on)
            tpg = mismatch_deflate(tpg, g_g, c_g, v_own_g)
            tpp = mismatch_deflate(tpp, g_p, c_p, v_own_r)
            tpn = mismatch_deflate(tpn, g_n, c_n, v_own_r)
            tmg = mismatch_deflate(tmg, g_g, c_g, v_own_g)
            tmp = mismatch_deflate(tmp, g_p, c_p, v_own_r)
            tmn = mismatch_deflate(tmn, g_n, c_n, v_own_r)
            # the composition stream is ONE DOF (λ), so its gap is measured on the λ axis — the message's
            # log(f_g/f_R) minus the own belief's, built from the SAME quantities the combine builds
            # ``lam_msg = mo_g − mo_R`` from, so the gap is exactly the error of the claim ψ receives. A
            # contradiction on EITHER arm contradicts the λ claim (a message with no RNA at all is asserting
            # λ = +∞, i.e. "this node is pure gDNA").
            g_R, c_R = mismatch_gap(pin_p + pin_n, op + on)
            _tau_pre = ttau
            ttau = mismatch_deflate(ttau, g_g - g_R, c_g | c_R, v_own_lam)
            if (
                _capture is not None
            ):  # inert: the per-message gaps + the τ-stream kill, for the dissect loop
                _capture.setdefault("_dl", []).append(
                    {
                        "G_g": g_g.copy(),
                        "G_p": g_p.copy(),
                        "G_n": g_n.copy(),
                        "G_lam": (g_g - g_R),
                        "contra": (c_g | c_R).copy(),
                        "tau_pre": _tau_pre.copy(),
                        "tau_post": ttau.copy(),
                    }
                )

            # ── THE MESSAGE PACKET: a LEVEL claim per component, and SEPARATE SPLIT + TILT claims ────────
            # A message must let the destination do three independent things: set each component's LEVEL
            # (``mo_g`` with ``cm_g``, ``mo_p``/``mo_n`` with ``cm_p``/``cm_n``), set the SPLIT (``λ``, with
            # ``τ``), and set the TILT (``θ``). Those are three different claims carried at three different
            # precisions — so the message states the split and the tilt EXPLICITLY. Reading them back off the
            # fused densities instead (as the combine used to) delivers them weighted by the LEVEL precisions
            # rather than their own, which is a mode/precision mismatch: a message with almost no composition
            # evidence but a large, well-counted density would set the split, and one with real composition
            # evidence but little mass would not.
            #
            # Both are scale-free — the node's mass cancels out of ``log(f_g/f_R)`` and out of the tilt — so a
            # message can state them without the destination reconstructing any density, and no rescaling of
            # the level claims can perturb them.
            _tR = tp + tn
            tlam = np.where(
                (tg > _EPS) & (_tR > _EPS),
                np.log(np.maximum(tg * E_g, _EPS)) - np.log(np.maximum(_tR * E_r, _EPS)),
                0.0,
            )
            tth = np.arcsin(
                np.clip(np.where(_tR > _EPS, (tp - tn) / np.maximum(_tR, _EPS), 0.0), -1.0, 1.0)
            )
            return tg, tp, tn, tpg, tpp, tpn, tmg, tmp, tmn, ttau, tlam, tth

        def _fuse_add(
            a, b
        ):  # additive (inverse-variance) fuse of two independent precision streams
            return np.asarray(a, np.float64) + np.asarray(b, np.float64)

        def _fuse_v(a, pa, b, pb):
            p = pa + pb
            return np.where(p > _EPS, (pa * a + pb * b) / np.maximum(p, _EPS), 0.0), p

        # ONE ρ-iteration (see the note at the top of the module), so this is straight-line: the frame
        # `rho0` is already built above, and there is no next iteration to feed. ⭐ The two calls now
        # differ ONLY in which neighbour they read — the face arguments dissolved with the faces.
        ag, ap, an, apg, app, apn, amg, amp, amn, atau, alam, ath = _transport(
            sl, vl, fwd, rho0, _pop_l_a
        )
        bg, bp, bn, bpg, bpp, bpn, bmg, bmp, bmn, btau, blam, bth = _transport(
            sr, vr, bwd, rho0, _pop_r_a
        )
        cg, cpg = _fuse_v(ag, apg, bg, bpg)  # density MODE (full precision-weighted)
        cp, cpp = _fuse_v(ap, app, bp, bpp)
        cn, cpn = _fuse_v(an, apn, bn, bpn)
        cm_g, cm_p, cm_n = _fuse_add(amg, bmg), _fuse_add(amp, bmp), _fuse_add(amn, bmn)
        c_tau = _fuse_add(atau, btau)
        mo_g = np.log(np.maximum(cg * E_g / np.maximum(M, _EPS), _EPS))
        mo_p = np.log(np.maximum(cp * E_r / np.maximum(M, _EPS), _EPS))
        mo_n = np.log(np.maximum(cn * E_r / np.maximum(M, _EPS), _EPS))
        # ── THE THREE-STREAM SINGLE-λ COMBINE (the M6 rank-1 fix) ──
        # (1) COMPOSITION → ONE λ-message on λ = log(f_g/f_R), precision ``c_tau`` (the fused Schur τ) — ψ
        #     counts the composition DOF ONCE, not twice. **Each claim is fused by ITS OWN precision**: the
        #     λ mode is the τ-weighted mean of the two messages' λ, not a ratio read back off the density
        #     fuse (where each component was averaged by its own MODE-FUSION precision). That mismatch
        #     delivered the split at a confidence it was never weighted by — a message with almost no
        #     composition evidence but a large, well-counted density could set it, and one with real
        #     composition evidence but little mass could not. Measured: 0.0889 → 0.0862, 13 better / 3
        #     worse, and it makes the M12 conservation rescale that preceded it completely INERT (0/0/32) —
        #     that operator had been correcting this same mismatch by a longer route.
        cR = cp + cn
        lam_msg = np.where(c_tau > _EPS, (atau * alam + btau * blam) / np.maximum(c_tau, _EPS), 0.0)
        # A λ message exists only where BOTH components of the pair reached this node — the structural
        # presence test. `_transport`'s per-message gate cannot catch every case: a message may carry an
        # RNA DENSITY while contributing zero mode-fusion PRECISION, in which case ``cR`` collapses here.
        c_tau = np.where((cg > _EPS) & (cR > _EPS), c_tau, 0.0)
        # (2) ANCHOR gDNA MEASUREMENT → gdna_imp (mode mo_g, precision ``cm_g``). (3) SPLICED RNA MEASUREMENT
        #     → rna_imp (mode mo_p/mo_n, precision ``cm_p``/``cm_n``). INDEPENDENT of the composition, so
        #     fused separately (an RNA-only spliced measurement constrains f_g via f_R with NO gDNA info).
        # (4) θ TILT (AMBIG) — fused by the MEASURED RNA that carries it, for the same reason as λ.
        _tha, _thb = amp + amn, bmp + bmn
        th_msg = np.where(
            (_tha + _thb) > _EPS,
            (_tha * ath + _thb * bth) / np.maximum(_tha + _thb, _EPS),
            0.0,
        )
        th_prec = np.where(is_amb, cm_p + cm_n, 0.0)
        # ⛔ The RNA MEASUREMENT ψ factor was ablated here behind a flag (now removed). Do not re-try it:
        # it carries 75 % of the posterior precision on the confidently-wrong unstranded × capture-OFF
        # exons, and removing it measured 0.0895 → 0.1033, 4 better / 17 worse — it is also the only thing
        # that lets a zero-gDNA library say "my mass is all RNA" (`gdna_none` 0.1063 → 0.1438).
        dc_fin = _local_solve(
            global_lp,
            mo_g,
            cm_g,
            (mo_p, mo_n),
            (cm_p, cm_n),
            lam_imp=(lam_msg, c_tau),
            theta_imp=(th_msg, th_prec),
        )
        nonlocal _uni_msg
        _uni_msg = (mo_g, cpg, mo_p, cpp, mo_n, cpn)  # publish for the shared diagnostics
        if _capture is not None:  # inert diagnostic: the fused per-component densities + the frames
            _capture.setdefault("_uni", []).append(
                {
                    "cg": cg.copy(),
                    "cp": cp.copy(),
                    "cn": cn.copy(),
                    "pg": cpg.copy(),
                    "pp": cpp.copy(),
                    "pn": cpn.copy(),
                    "ag": ag.copy(),
                    "ap": ap.copy(),
                    "an": an.copy(),
                    "bg": bg.copy(),
                    "bp": bp.copy(),
                    "bn": bn.copy(),
                    "apg": apg.copy(),
                    "app": app.copy(),
                    "bpg": bpg.copy(),
                    "bpp": bpp.copy(),
                    "rho_tot": rho0.copy(),
                    "mo_g": mo_g.copy(),
                    "mo_p": mo_p.copy(),
                    "mo_n": mo_n.copy(),
                    "lam_msg": lam_msg.copy(),
                    "c_tau": c_tau.copy(),
                    "cm_g": cm_g.copy(),
                    "cm_p": cm_p.copy(),
                    "cm_n": cm_n.copy(),
                    "amp": amp.copy(),
                    "bmp": bmp.copy(),
                    "amn": amn.copy(),
                    "bmn": bmn.copy(),
                    "cpg": cpg.copy(),
                    "fg_out": np.clip(np.asarray(dc_fin.gdna_frac, np.float64), 0.0, 1.0),
                }
            )
        if _capture is not None:
            _capture.update(
                _uni_static={
                    "M": M,
                    "E_g": E_g,
                    "E_r": E_r,
                    "spl_p": spl_p,
                    "spl_n": spl_n,
                    "og": og,
                    "op": op,
                    "on": on,
                    "pg_own": pg_own,
                    "pp_own": pp_own,
                    "pn_own": pn_own,
                    "fwd_g": fwd[0],
                    "fwd_p": fwd[1],
                    "fwd_n": fwd[2],
                    "fwd_pg": fwd[3],
                    "fwd_pp": fwd[4],
                    "fwd_pn": fwd[5],
                    "bwd_g": bwd[0],
                    "bwd_p": bwd[1],
                    "bwd_n": bwd[2],
                    "bwd_pg": bwd[3],
                    "bwd_pp": bwd[4],
                    "bwd_pn": bwd[5],
                    "rho_node0": rho_node0,
                    "rho0": rho0,
                    # ── AUDIT_2 instrumentation (invariant scan) ──
                    "order": np.asarray(order_list, np.int64),
                    "logvar_tot": logvar_tot.copy(),
                    "mature_pos": SPL[:, 0].copy(),
                    "mature_neg": SPL[:, 1].copy(),
                    "n_slot": n_slot.copy(),
                    "tau_own": tau_own.copy(),
                    "mg_own": mg_own.copy(),
                    "struct_lock": _struct.copy(),
                    "fwd_mp": fwd[7],
                    "bwd_mp": bwd[7],
                    "fwd_mn": fwd[8],
                    "bwd_mn": bwd[8],
                    "fwd_mg": fwd[6],
                    "bwd_mg": bwd[6],
                    "fwd_tau": fwd[9],
                    "bwd_tau": bwd[9],
                    "is_bnd": is_bnd_a,
                    "is_exon": ex_a,
                    "left": li,
                    "right": ri,
                    "fp": fp_a,
                    "fn": fn_a,
                }
            )
        return dc_fin

    dc_fin = _unified_solve()
    # the final per-component imputation factors, published for the diagnostics (`_uni_msg`, set at the end of
    # `_unified_solve` — `pass0_node_dissect.py`'s ψ channel-ablation replay consumes them).
    mode_g, prec_g, mode_p, prec_p, mode_n, prec_n = _uni_msg
    mg_, mp_, mn_ = dc_fin.gdna_frac, dc_fin.rna_pos_frac, dc_fin.rna_neg_frac
    vg_, vp_, vn_ = dc_fin.gdna_frac_var, dc_fin.rna_pos_frac_var, dc_fin.rna_neg_frac_var
    # write back only SOLVABLE nodes (G1 sinks / empty keep their signature-binary init). The §6B DOF SOLVE-GATE
    # (skip unidentified nodes → keep the f_g=1 init, defer to the prior) was DERIVED, IMPLEMENTED, and
    # EMPIRICALLY REFUTED: it regresses both standalone (refit=0 +0.010) and with the
    # hyperprior (refit=1 +0.025) — the prior resolves an imperfectly-SOLVED node better than a deferred f_g=1.
    f_g = np.where(solvable, np.clip(mg_, 0.0, 1.0), f_g)
    f_pos = np.where(solvable, np.clip(mp_, 0.0, 1.0), f_pos)
    f_neg = np.where(solvable, np.clip(mn_, 0.0, 1.0), f_neg)
    var_g = np.where(solvable, vg_, var_g)
    var_pos = np.where(solvable, vp_, var_pos)
    var_neg = np.where(solvable, vn_, var_neg)

    if _capture is not None:  # inert diagnostic hook
        # strand-ONLY local belief (no global prior, no messages) — to split the LOCAL error into the
        # strand likelihood vs the global gDNA prior contribution. Same log-density solver, global=None.
        fg_strand = _solve_nodes_logodds_all(
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
        # the message-free local self-solve variances (for the local-error attribution) — a debug-only solve
        # so the production path carries none of this (the self-solve fractions come from `_ni`).
        _dc_loc = _local_solve(global_lp)
        _capture.update(
            fg_loc=_ni.f_g,
            fg_strand=fg_strand,
            fp_loc=_ni.f_pos,
            fn_loc=_ni.f_neg,
            vg_loc=_dc_loc.gdna_frac_var,
            vp_loc=_dc_loc.rna_pos_frac_var,
            mode_g=mode_g,
            prec_g=prec_g,
            mode_p=mode_p,
            prec_p=prec_p,
            mode_n=mode_n,
            prec_n=prec_n,
            f_g=f_g.copy(),
            f_pos=f_pos.copy(),
            f_neg=f_neg.copy(),
            var_g=var_g.copy(),
            solvable=solvable,
            # boundary-emission geometry: gDNA emits iff facing unspliced mass>0 (strand-agnostic);
            # RNA iff free_s on both endpoints & (unspliced or spliced facing mass). Capture the faces.
            count=CNT,
            spliced=spliced_slot,
            mature=SPL.sum(axis=1),
            free_pos=np.asarray(fp, bool),
            free_neg=np.asarray(fn, bool),
            # global geometry (μ = clip(ρ·eff_global/mass_global, 0, 1) is the implied prior fraction).
            eff_global=eff_global,
            mass_global=mass_global,
            # per-face geometry for message dissection (logodds diagnostics)
            eff_gdna=EG,
            eff_rna=ER,
            # the full per-node global prior term on the solve grid (strand-free), so a diagnostic can
            # replay _solve_nodes_logodds_all with message channels ablated (message help/hurt attribution).
            global_lp=global_lp,
            solve_grid=solve_grid,
            # DIAGNOSTIC (inert in production): the composition-evidence seed. (The retired NPMLE projection
            # ``_mu_proj``/``_var_proj`` is gone; the LIVE enrichment scale + its variance are ``rho_node0`` and
            # ``logvar_tot`` in ``_uni_static``, from which σ²_transfer = logvar_tot[dst]+logvar_tot[src].)
            _tau0_lam=_ni.tau_lam,
            # the incoming belief (the final solve's ``fg_ref``) + the intron-factory λ arm, so an ablation
            # replay of `_solve_nodes_logodds_all` reproduces the shipped f_g exactly before ablating.
            fg_init=_fg_init,
            fpos_init=_fp_init,
            fneg_init=_fn_init,
            intron_prior=intron_prior,
            solvable_mask=solvable,
        )

    return NodeBelief(
        f_pos=f_pos,
        f_neg=f_neg,
        f_g=f_g,
        var_pos=var_pos,
        var_neg=var_neg,
        var_gdna=var_g,
    )


def chain_node_deconv(chain: NodeChain, belief: NodeBelief, substrate) -> NodeDeconv:
    """Project the chain belief's NODE slots back onto the NODE axis as a :class:`NodeDeconv` — what
    ``CalibrationResult`` / ``priors`` / ``derive`` consume.

    ⚠ **A node's contained population carries no spliced term any more, and that is structural**: the
    accumulator credits ``node_contained`` only when the fragment used no junction, so a contained
    fragment is unspliced by construction. The predecessor added ``+ mass_spliced`` here; that quantity
    is identically zero on the node axis now, and adding it would be adding a channel that cannot exist.
    """
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    reg = kind == NODE
    count = np.asarray(substrate.node_contained.count, dtype=np.float64).sum(axis=1)
    n = count.shape[0]
    f_g = np.zeros(n)
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    ri = idx[reg]
    f_g[ri] = belief.f_g[reg]
    f_pos[ri] = belief.f_pos[reg]
    f_neg[ri] = belief.f_neg[reg]
    return NodeDeconv(
        gdna_mass=f_g * count,
        rna_mass=(1.0 - f_g) * count,
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
    )


def chain_edge_deconv(chain: NodeChain, belief: NodeBelief, substrate) -> NodeDeconv:
    """Project the chain belief's EDGE slots onto the CONTIGUOUS-EDGE axis — the crossing flux that
    ``priors`` and ``derive`` consume.

    ⭐ **ONE per-edge result, not a ``(left, right)`` pair of per-region ones.** The predecessor split
    each edge's flux onto its two flanking regions and ``priors`` then pooled the two halves straight
    back together — so the split and the re-pool were a no-op, and that exact sum-then-halve pattern is
    what hid a factor of 2 for months. Owner ruling, 2026-07-30:
    ``CalibrationResult``'s per-region ``mass_*_left/right`` become per-edge arrays.

    The RNA mass is spliced-inclusive: an edge's certified-RNA crossings (``edge_spliced``) are RNA
    whatever the unspliced mixture resolves to, since gDNA cannot be spliced.
    """
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    edge = kind == EDGE
    unspliced = np.asarray(substrate.edge_unspliced.count, dtype=np.float64).sum(axis=1)
    spliced = np.asarray(substrate.edge_spliced.count, dtype=np.float64).sum(axis=1)
    n = unspliced.shape[0]
    f_g = np.zeros(n)
    f_g[idx[edge]] = np.asarray(belief.f_g, dtype=np.float64)[edge]
    return NodeDeconv(
        gdna_mass=f_g * unspliced,
        rna_mass=(1.0 - f_g) * unspliced + spliced,
        gdna_frac=f_g,
        rna_pos_frac=np.zeros(n),
        rna_neg_frac=np.zeros(n),
    )
