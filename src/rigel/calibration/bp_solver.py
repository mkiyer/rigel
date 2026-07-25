"""The belief-propagation calibration solver: per-node gDNA-vs-RNA deconvolution by a bidirectional sweep.

Deconvolves each node's unspliced fragment mass into a pie ``(f_pos, f_neg, f_g)`` — sense-RNA /
antisense-RNA / gDNA — over the unified region↔boundary chain (`node_chain`), by a single forward-backward
(L→R then R→L) belief-propagation pass (exact on the chain, a forest of linear paths). Each per-node solve
(`simplex_logodds`, the log-density log-odds solver) reconciles three sources of information: the intrinsic
strand likelihood (the Beta-Binomial tilt — the only signal a count carries), the cross-node imputation
messages, and the population gDNA prior. Theory: the count-zero-information principle in
`docs/calibration/CALIBRATION_ARCHITECTURE.md`; the composition (enrichment-ratio) message model in
`docs/calibration/unified_solver_design.md`.

**The UNIFIED (composition) solver — one message mode.** Each message is REFRAMED into the destination's
frame by the enrichment ratio ``r = ρ_tot(dst)/ρ_tot(src)`` (`node_total_density`, lazy + composition-aware,
per-side spliced), the inactive strands are FILTERED at the destination, the mature is ROUTED (graft into an
exon / peel out of an exon, via the boundary's measured spliced), then each component becomes a factor on the
destination's ``f_c`` by the density mode ÷ its own observed mass ``M_dst`` (the pure arithmetic in
`enrichment_frame`). The forward-backward relay carries per-component densities + precisions; the combine
transports both neighbours into the node's frame and runs the ψ solve. ⚠ **The message VARIANCE model
(``σ²_transfer`` / the per-hop damping) is the density-uniformity proxy and is being REDONE**
(`docs/calibration/variance_model_handoff.md` §3-4) — it is deliberately left as-is here.

Module layout. The per-node geometry / belief / statics / init primitives and the pure geometry helpers
(`build_node_geometry`, `build_node_statics`, `init_beliefs`, `node_global_geometry`, `node_total_density`,
`NodeGeometry`/`NodeBelief`/`NodeStatics`) live in the lower `node_geometry` module; the per-node
INITIALIZATION self-solve (the four sources → each node's own ``(density, precision)``) lives in `node_init`
(`build_node_init`). Both are re-exported here for the calibrator's convenience; this module owns:
* `node_sweep` — the single forward-backward unified sweep (build the self-solve → relay → combine → ψ solve).
* `chain_region_deconv` / `chain_boundary_side_deconv` — project the converged belief back to the per-region
  / per-boundary-side masses the `CalibrationResult` consumes.
"""

from __future__ import annotations

import os

import numpy as np

from .enrichment_frame import (
    composition_logvar,
    graft_frame_logvar,
    mismatch_deflate,
    mismatch_gap,
    transfer_logvar,
)
from .node_chain import REGION, NodeChain
from .signature import coarse_type_array
from .node_geometry import (
    NodeBelief,
    NodeGeometry,
    NodeStatics,
    build_node_geometry,
    build_node_statics,
    init_beliefs,
    node_global_geometry,
    node_total_density,
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
    "chain_region_deconv",
    "chain_boundary_side_deconv",
]

_EPS = 1.0e-9

# The number of lazy-ρ_tot combine iterations (`unified_solver_design.md` §10.6): iteration 2 recomputes
# ρ_tot from the both-message composition. A numerical-resolution knob, not a model constant.
_RHO_ITERS = 2


def node_sweep(
    chain: NodeChain,
    statics: NodeStatics,
    geometry: NodeGeometry,
    belief: NodeBelief,
    region_arrays,
    boundary_substrate,
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
    fold_coarse_k: int = 33,
    fold_fine_k: int = 33,
    fold_sigma_coverage: float = 6.0,
    fold_refine_iters: int = 3,
    _capture: dict | None = None,
):
    """The belief-propagation sweep over the chain — gDNA AND per-strand RNA messages, in COUNT space (no
    ``(M/E)²`` Jacobian) — as a single **FORWARD-BACKWARD** pass.

    The chain is a forest of linear paths, so BP is exact in one forward + one backward pass (vs Gauss-Seidel /
    Jacobi which propagate one hop per pass). A message's precision is the source's own HONEST belief precision
    degraded by the two independent defects a cross-node imputation suffers
    (`docs/calibration/message_variance_derivation.md`)::

        p = 1 / ( Var(log f_c^src) + 1/n_src  +  σ²_transfer  +  b̂² )
                 \\__ strand ___/   \\_count_/    \\_ SCALE _/    \\_ COMPOSITION _/

    — the composition and count precision the source actually earned, the reframe's own scale uncertainty
    ``σ²_transfer = Var(log r)`` (M5; 0 on the matched graft where ``r`` is common-mode), and the
    DerSimonian–Laird estimate of how wrong the imputation PREMISE itself is (``b̂²``, M7 — the message's
    composition against the destination's independent self-solve). All four are computed inside the pass from
    counts and effective lengths: nothing is fitted, there is no precision to refit and no outer fixed-point
    loop. The global prior is ANCHORED (every input fit once before the solve), so the single FB pass is exact.

    BEFORE the pass: the pass-0 NPMLE gDNA hyperprior (:class:`~.npmle.DensityNPMLE`), fit once,
    belief-free, and passed as ``gdna_prior``. ``gdna_prior=None`` is a first-class PRIOR-FREE solve: ψ then
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

    # per-face geometry as (left, right) pairs, indexed by face (0=left, 1=right).
    EG = (geometry.eff_gdna_left, geometry.eff_gdna_right)
    ER = (geometry.eff_rna_left, geometry.eff_rna_right)
    ESP = (geometry.eff_spl_left, geometry.eff_spl_right)  # one-sided spliced half-triangle eff-len
    MS = (geometry.mass_left, geometry.mass_right)
    SP = (geometry.spliced_pos_left, geometry.spliced_pos_right)
    SN = (geometry.spliced_neg_left, geometry.spliced_neg_right)
    # per-node "global" gDNA support (region = contained; boundary = both-side crossing over the averaged
    # per-side density length) — the basis the pass-0 rate prior is fit + projected on.
    mass_global, eff_global = node_global_geometry(chain, geometry)

    # The per-node solve is the log-density 1-D/2-D log-odds solver (simplex_logodds, O(m·K),
    # genome-scale-tractable). The "solve grid" is the f_g axis the global NB prior is evaluated on (the
    # log-odds σ(λ) lattice).
    _lam_lo, _fg_lo = _logodds_grid(int(n_grid), float(logodds_window))
    solve_grid = _fg_lo
    kappa = float(rna_sense_frac)
    od_g, od_r = gdna_strand_overdispersion, rna_strand_overdispersion

    def _local_solve(
        g_arr, gm=None, gp=None, rm=None, rp=None, lam_imp=None, theta_imp=None
    ):
        """The per-node local/final solve (log-density log-odds backend). Returns the :class:`NodeDeconv`
        (the readout ``*_frac``/``*_frac_var`` + the free-coordinate seed ``lam_mean``/``lam_var``/
        ``theta_mean``/``theta_var``). Phase A calls it message-free; phase D passes the FB messages.
        ``lam_imp``/``theta_imp`` are the SINGLE-λ composition message + the θ tilt message (2-tuples of
        ``(mode, prec)``), the rank-1 fix that replaces the two per-component ``gm``/``rm`` messages."""
        return _solve_nodes_logodds_all(
            statics.u_pos,
            statics.u_neg,
            fp,
            fn,
            statics.mass_unspliced,
            statics.mass_spliced,
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
    solvable = (fp | fn) & (statics.mass_unspliced > 0.0)

    # THE gDNA ARM of ψ — the COMPOSITION prior. The NPMLE's two roles are kept SEPARATE
    # (docs/calibration/CALIBRATION_MASTER.md §5): this ``gdna_prior`` is the COMPOSITION arm ONLY.
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

    # Genomic node order for the forward/backward scans (within each ref path; left/right break at −1). The
    # scans are sequential per node, so iterate as a Python list of ints (faster than numpy scalar indexing).
    order_list = [int(x) for x in np.asarray(chain.order)]
    # Per-node EXON-region flag (coarse_type == 2) — the unified relay routes mature into EXON destinations
    # (the graft) and peels it out of EXON sources (`ex_a` in `_unified_solve`).
    _rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    _ri = np.clip(np.asarray(chain.ref_idx, dtype=np.int64), 0, _rtype.shape[0] - 1)
    is_exon_node = ((np.asarray(chain.kind) == REGION) & (_rtype[_ri] == 2)).tolist()

    # ─────────────────────────────────────────────────────────────────────────────────────────────────────
    # THE UNIFIED SOLVER (unified_solver_design.md; owner 2026-07-23) — ONE mode: reframe → filter → route → ÷M.
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
    )

    def _unified_solve():
        n = f_g.shape[0]
        is_reg_a = np.asarray(chain.kind) == REGION
        is_bnd_a = ~is_reg_a
        ex_a = np.asarray(
            is_exon_node, dtype=bool
        )  # source-is-exon selector for the mature routing
        fp_a, fn_a = np.asarray(fp, bool), np.asarray(fn, bool)
        is_amb = fp_a & fn_a  # AMBIG: both strands live → the tilt θ is a free DOF (its own message)
        M = np.asarray(mass_global, np.float64)
        E_g = np.asarray(eff_global, np.float64)
        E_r = np.where(is_reg_a, ER[0], ER[0] + ER[1]).astype(
            np.float64
        )  # node-level RNA-FL eff-length

        # ── M5 σ²_transfer = Var(log r): the per-node total-density log-variance (`message_variance_derivation.md`
        # §3; `enrichment_frame.composition_logvar`). Var(log ρ_tot) = 1/n + [(1/E_g−1/E_r)/B]²·Var(f_g), with
        # Var(f_g) = (f_g(1−f_g))²/τ_λ CAPPED at f_g(1−f_g) (a fraction's max variance) and 0 for a
        # composition-certain (struct_lock) node. Evaluated at the INPUT belief f_g (consistent with the reframe
        # densities rho_l0/rho_r0). This is the HONEST, prior-free transfer variance that RETIRES the
        # density-uniformity NPMLE proxy (which was identically 0 in pass-0, so pass-0 had NO transfer damping).
        _n_node = np.where(
            is_reg_a,
            np.asarray(geometry.n_unspl_left, np.float64),
            np.asarray(geometry.n_unspl_left, np.float64)
            + np.asarray(geometry.n_unspl_right, np.float64),
        )
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
        _s2t_off = bool(os.environ.get("RIGEL_S2T_OFF"))  # DIAGNOSTIC ONLY: disable BOTH cliff terms

        # per-strand spliced (mature) DENSITY — one-sided per-face density (0 on regions); the acceptor face.
        def _face_spl(sp):
            return np.where(sp[0] > _EPS, sp[0] / np.maximum(ESP[0], _EPS), 0.0) + np.where(
                sp[1] > _EPS, sp[1] / np.maximum(ESP[1], _EPS), 0.0
            )

        spl_p = _face_spl(SP)  # + mature density at a boundary, both faces pooled (0 on regions)
        spl_n = _face_spl(SN)
        # PER-FACE mature density, indexed [face] — what the routing must use. The node-pooled ``spl_*`` above
        # double-counts on an exon↔exon junction (a donor on one flank, an acceptor on the other): each
        # direction would peel/graft the SUM of both flanks' flux. This is the A3 per-face fix ``_scan`` already
        # carries in ``mature_dilution[df]``; the routing needs the same face selection.
        spl_p_f = tuple(
            np.where(SP[k] > _EPS, SP[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1)
        )
        spl_n_f = tuple(
            np.where(SN[k] > _EPS, SN[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1)
        )
        # ── AN EXON HAS NO RNA CLAIM ACROSS A SPLICE JUNCTION ──────────────────────────────────────────
        # A boundary carrying spliced mass on either face has an ACTIVE junction. Across it, MATURE RNA
        # splices out by definition — what continues unspliced is the nascent read-through, and an exon
        # cannot measure how much of its own RNA that is. So the exon's RNA density is not a claim about what
        # is on the other side, and the honest emission is NONE (density AND precision), not a partial peel.
        #
        # The shipped peel subtracts the MEASURED spliced flux as a proxy for "the mature that leaves". It
        # under-consumes badly. Measured on `nrna_none` (where the oracle RNA at an exon|intron boundary and
        # inside the intron is EXACTLY 0): the relay still delivered ρ_R = 0.53 at the boundary and 0.0041 at
        # the intron, the peel removed only ~16 % of the reframed RNA, and it did not fire at all on 14–17 %
        # of edges (no spliced mass measured on that face). Meanwhile gDNA transported essentially perfectly
        # across the same hops (relayed/oracle 0.96–0.99) — so the RNA channel alone was contaminating the
        # chain, and introns (which the density deconvolution solves accurately) were the visible victims.
        # This is the component-set refusal `enrichment_frame.gdna_fallback_admissible` already describes:
        # where the two ends do not share a component set, emit nothing rather than a wrong frame.
        _junc = ((SP[0] + SN[0] + SP[1] + SN[1]) > _EPS) & is_bnd_a
        accept_l = (
            SP[0] + SN[0]
        ) > _EPS  # the LEFT face carries the spliced (acceptor) ⇒ WITH-spliced ρ_tot
        accept_r = (SP[1] + SN[1]) > _EPS

        # own per-component densities + precisions — the message-free SELF-SOLVE (`node_init.build_node_init`,
        # the four sources). ``rho_*`` are the own densities; ``prec_*`` combine the strand + intron-factory
        # composition evidence ``τ_λ`` (Var(log f_g)=(1−f_g)²/τ, Var(log f_r)=f_g²/τ) with the Poisson count
        # power — 0 at an uninformed node (the honest unstranded statement). Both RNA arms are built from the
        # SAME τ as gDNA (never the local posterior variance, which pools the shared reference into a phantom).
        og, op, on = _ni.rho_g, _ni.rho_pos, _ni.rho_neg
        pg_own, pp_own, pn_own = _ni.prec_g, _ni.prec_pos, _ni.prec_neg

        # ── THREE-STREAM precision seeds (the single-λ combine, message_variance_derivation.md §4 / HANDOFF_4
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

        # ── the OWN-belief composition variances: the DL estimator's "second study" ─────────────────────────
        # Reused from `node_init` so ONE law defines them (the τ_λ FOUNDATION Jacobians). Three states, and the
        # DL term's three regimes fall out of them with no gate and no constant:
        #   * struct_lock (composition CERTAIN) → v_own = 0   ⇒ excess = G² in full. NOTE this regime is INERT at
        #     the combine: `struct_lock = ~solvable & is_region` ⊆ `~solvable`, so such a node's ψ output is
        #     discarded by the write-back gate anyway (and `op+on = 0` there kills the λ gap). What actually
        #     makes an intergenic anchor immovable is its own `pg_own = n` in the relay's fuse, NOT this branch;
        #     it is kept because it is the correct limit, not because it is load-bearing today;
        #   * real evidence τ_own > 0          → v_own finite ⇒ a message that CONFLICTS with a confident own
        #     belief is killed (the stranded arm's fix), one that AGREES is barely touched;
        #   * no evidence τ_own = 0            → v_own = ∞    ⇒ excess ≡ 0: NO mismatch damping. Every AMBIG
        #     node and all of unstranded data sit here — exactly where messages are the ONLY information — so the
        #     M5 unstranded/capture win propagates untouched. (Phase 2's hyperprior is what supplies a finite
        #     v_own there; until then the honest statement is "this node has no opinion to contradict".)
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
        def _pin_v(
            g, p, n, pg_, pp_, pn_
        ):  # per-message pin (the combine): scale a claim to the node's mass
            sg = np.where(pg_ > 0.0, g, og)
            sp = np.where(pp_ > 0.0, p, op)
            sn = np.where(pn_ > 0.0, n, on)
            s = sg * E_g + (sp + sn) * E_r
            k = np.where((s > _EPS) & (M > _EPS), M / np.maximum(s, _EPS), 1.0)
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
        def _rho_faces(fgc):
            """Lazy, composition-aware ρ_tot from the current f_g, split per side (WITH spliced at the acceptor)."""
            ru, rw = node_total_density(chain, geometry, fgc)
            rs = rw - ru  # the one-sided spliced density
            return (
                ru,
                ru + np.where(accept_l, rs, 0.0),
                ru + np.where(accept_r, rs, 0.0),
            )  # node, left-face, right-face

        # ── the RELAY: accumulate the forward/backward context belief (densities), reframed each hop by the
        # INPUT-belief ρ_tot. Returns each node's context belief IN ITS OWN FRAME; the combine re-reframes.
        rho_node0, rho_l0, rho_r0 = _rho_faces(np.asarray(f_g, np.float64))

        # σ²_transfer per-hop damping: add the edge's transfer log-variance to the message's log-variance
        # (1/p → 1/p + σ²_transfer). This is now the DERIVED, direction-dependent M5 variance (0 on the matched
        # graft, Var(log r) on the peel/plain/anchor — computed per edge from ``logvar_tot`` below), retiring the
        # density-uniformity NPMLE proxy (`_mup`/`_vp`, no longer read).
        def _damp(p, s2t):
            return 1.0 / (1.0 / p + s2t) if p > 0.0 else 0.0

        def _relay(seq, nbr, dst_face, src_face, df, sf):
            rg, rp, rn = og.copy(), op.copy(), on.copy()
            pg, pp, pn = pg_own.copy(), pp_own.copy(), pn_own.copy()  # full → MODE fusion
            mg, mp, mn = mg_own.copy(), mp_own.copy(), mn_own.copy()  # MEASUREMENT (anchor gDNA + spliced RNA)
            tau = tau_own.copy()  # COMPOSITION (τ_λ) → the λ-message
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
                rho_src = src_face[s]
                r = (
                    (dst_face[i] / rho_src) if (rho_src > _EPS and dst_face[i] > _EPS) else 1.0
                )  # no frame ⇒ pass-through
                # GRAFT (boundary → EXON, §6): the boundary's measured mature is a density AT THE SOURCE, so it
                # joins the source's RNA BEFORE the reframe; the peel is measured at the destination and so is
                # applied after. Both use the face that FACES the other endpoint (``sf``/``df``), never the
                # node-pooled sum. Only an EXON receives the graft — an intron carries no mature (`ex_a[i]`,
                # not `is_reg_a[i]`, which grafted the junction's whole mature flux into every flanking intron).
                _gr = ex_a[i] and is_bnd_a[s]
                # σ²_transfer = Var(log r) (M5): 0 on the matched-set GRAFT (r is common-mode across {g,R} and
                # cancels in the composition — charging it there is a double-count), Var(log r) elsewhere (peel /
                # plain reframe / partial-anchor, where r is load-bearing). The COMPOSITION-mismatch term b̂² is
                # the combine's job (`_transport`): the relay has no destination self-solve to measure a gap
                # against — its running belief is already fused with the messages, so a DL gap here would be
                # feedback, not evidence.
                s2t = 0.0 if (_s2t_off or _gr) else (logvar_tot[i] + logvar_tot[s])
                gp = spl_p_f[sf][s] if _gr else 0.0
                gn = spl_n_f[sf][s] if _gr else 0.0
                tg, tp, tn = rg[s] * r, (rp[s] + gp) * r, (rn[s] + gn) * r
                tpg, tpp, tpn = _damp(pg[s], s2t), _damp(pp[s], s2t), _damp(pn[s], s2t)  # full (mode)
                tmg, tmp, tmn = _damp(mg[s], s2t), _damp(mp[s], s2t), _damp(mn[s], s2t)  # measurement
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
                    _s2f = s2t + (0.0 if _s2t_off else float(graft_frame_logvar(r)))
                    _spc = SP[sf][s] / (1.0 + SP[sf][s] * _s2f) if SP[sf][s] > _EPS else 0.0
                    _snc = SN[sf][s] / (1.0 + SN[sf][s] * _s2f) if SN[sf][s] > _EPS else 0.0
                    tpp += _spc
                    tpn += _snc
                    tmp += _spc
                    tmn += _snc
                if is_bnd_a[i] and ex_a[s]:  # EXON → boundary: PEEL the departing mature
                    if _junc[i]:  # across a junction the exon has NO RNA claim — see `_junc` above
                        tp, tn, tpp, tpn, tmp, tmn = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                    else:
                        tp, tn = max(tp - spl_p_f[df][i], 0.0), max(tn - spl_n_f[df][i], 0.0)
                if not fp_a[i]:
                    tp, tpp, tmp = 0.0, 0.0, 0.0
                if not fn_a[i]:
                    tn, tpn, tmn = 0.0, 0.0, 0.0
                # ── ANCHOR THE CONTEXT TO THIS NODE'S OBSERVED MASS (the scalar twin of `_pin_v`) ───────────
                # `Σ_c ρ_c·E_c = M` is an IDENTITY under the imputation premise, not an approximation: a matched
                # reframe delivers ρ_c^msg = a_c·ρ_tot(dst) = ρ_c^dst,true, so the components account for exactly
                # the fragments the node observed. Its violation therefore MEASURES the premise error — against a
                # hard observable, with no prior. Enforcing it here is what `_pin_v`'s own docstring already
                # prescribes ("applied at EVERY node rather than only at the final combine"); until now the relay
                # ran unanchored and the residual COMPOUNDED — a multiplicative random walk reaching
                # `Σ_c ρ_c E_c / M` p99 = 31–288× and max 519× (52–71 % of nodes over 1).
                # Note what this does NOT do: one hop's violation from composition error alone is bounded by the
                # eff-length ratio (k = [Σ a_c^dst E_c]/[Σ a_c^src E_c] ⇒ ×1.04 on a contained region, ×1.50 at a
                # boundary crossing), and plain reframes measure exactly that (median ×1.05–1.11). The heavy tail
                # is the ROUTING: graft p90 ×11.6–84.4, peel median ×1.31–1.58 — they add/subtract an ABSOLUTE
                # measured density into a RELATIVE claim. So this is an anti-drift anchor, not a composition test.
                # The `_pin_v` semantics are load-bearing: a component the context does not supply is filled from
                # the node's OWN density, so a PARTIAL claim stays partial. Rescaling all three blindly instead
                # regresses capture-OFF 3.6× (derivation: `scratchpad/derive_2_relay_pin.py`).
                _sg = tg if tpg > 0.0 else og[i]
                _sp = tp if tpp > 0.0 else op[i]
                _sn = tn if tpn > 0.0 else on[i]
                _sv = _sg * E_g[i] + (_sp + _sn) * E_r[i]
                if _sv > _EPS and M[i] > _EPS:
                    _k = M[i] / _sv
                    tg, tp, tn = tg * _k, tp * _k, tn * _k
                rg[i], pg[i] = _fuse(og[i], pg_own[i], tg, tpg)
                rp[i], pp[i] = _fuse(op[i], pp_own[i], tp, tpp)
                rn[i], pn[i] = _fuse(on[i], pn_own[i], tn, tpn)
                # measurement + composition precisions are INDEPENDENT evidence → additive (inverse-variance) fuse
                mg[i] = mg_own[i] + tmg
                mp[i] = mp_own[i] + tmp
                mn[i] = mn_own[i] + tmn
                tau[i] = tau_own[i] + ttau
            return rg, rp, rn, pg, pp, pn, mg, mp, mn, tau

        # dst faces its source on its LEFT (face 0); the source faces the dst on its RIGHT (face 1) — and mirrored
        # for the backward pass. The face pair selects BOTH the per-side ρ_tot and the per-face mature flux.
        fwd = _relay(order_list, left, rho_l0, rho_r0, 0, 1)
        bwd = _relay(order_list[::-1], right, rho_r0, rho_l0, 1, 0)

        # ── the COMBINE: transport α (from left neighbour) + β (from right neighbour) into the node's frame with
        # the LAZY ρ_tot (two-iteration — the 2nd uses the both-message composition), fuse, ÷M_dst → the ψ solve.
        li, ri = np.asarray(left), np.asarray(right)
        vl, vr = li >= 0, ri >= 0
        sl, sr = np.clip(li, 0, n - 1), np.clip(ri, 0, n - 1)

        def _transport(src, valid, df, sf, fwd_arrs, dst_face_v, src_face_v):
            rg, rp, rn, pg, pp, pn, mg, mp, mn, tau = fwd_arrs
            # A node with no frame (no mass ⇒ no ρ_tot, §5) cannot reframe: the message passes through at r=1.
            # Falling back to ``rho_src = 1.0`` instead made r the destination's ABSOLUTE density (10³ on a
            # short node) — a raw scale masquerading as a ratio. The relay already guards this way.
            framed = valid & (src_face_v[src] > _EPS) & (dst_face_v > _EPS)
            r = np.where(
                framed, dst_face_v / np.maximum(src_face_v[src], _EPS), np.where(valid, 1.0, 0.0)
            )
            # GRAFT before the reframe (a density measured AT the source); PEEL after (measured at the dst).
            # Both per-FACE, and the graft only into an EXON — see the relay's twin.
            graft = ex_a & is_bnd_a[src] & valid
            gp = np.where(graft, spl_p_f[sf][src], 0.0)
            gn = np.where(graft, spl_n_f[sf][src], 0.0)
            tg, tp, tn = rg[src] * r, (rp[src] + gp) * r, (rn[src] + gn) * r
            # σ²_transfer = Var(log r) (M5, the tested pure law): 0 on the matched-set graft (r common-mode ⇒
            # cancels — a double-count otherwise), Var(log r) = logvar_tot[dst]+logvar_tot[src] elsewhere (peel /
            # plain reframe / partial-anchor — r load-bearing). This is the SCALE half of the cliff cost; the
            # COMPOSITION half is the DL b̂² applied at the end of this function. See the relay's twin.
            s2t = (
                np.zeros_like(r) if _s2t_off else transfer_logvar(logvar_tot, logvar_tot[src], graft)
            )

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
            _sp, _sn = np.where(graft, SP[sf][src], 0.0), np.where(graft, SN[sf][src], 0.0)
            _s2t_spl = np.where(np.isfinite(s2t), s2t, 0.0)
            # M8 — the graft's FRAME-MISLIFT variance (see the relay's twin and `graft_frame_logvar`): the
            # measured spliced already sits in the destination exon's frame, so ``r`` is NOT common-mode for it
            # and M5's graft-zero does not cover it. 0 where r = 1, so it is inert without a capture step.
            _s2t_spl = _s2t_spl + (
                np.zeros_like(r) if _s2t_off else np.where(graft, graft_frame_logvar(r), 0.0)
            )
            _spc = np.where(_sp > _EPS, _sp / (1.0 + _sp * _s2t_spl), 0.0)
            _snc = np.where(_sn > _EPS, _sn / (1.0 + _sn * _s2t_spl), 0.0)
            tpp, tpn = tpp + _spc, tpn + _snc  # into the mode-fusion precision …
            tmp, tmn = tmp + _spc, tmn + _snc  # … and the measurement stream (a count, never composition τ)
            peel = is_bnd_a & ex_a[src] & valid  # EXON → boundary: PEEL the departing mature
            _kill = peel & _junc  # no RNA claim across a junction (see the relay's twin)
            tp = np.where(peel, np.maximum(tp - spl_p_f[df], 0.0), tp)
            tn = np.where(peel, np.maximum(tn - spl_n_f[df], 0.0), tn)
            tp, tn = np.where(_kill, 0.0, tp), np.where(_kill, 0.0, tn)
            tpp, tpn = np.where(_kill, 0.0, tpp), np.where(_kill, 0.0, tpn)
            tmp, tmn = np.where(_kill, 0.0, tmp), np.where(_kill, 0.0, tmn)
            tp, tpp, tmp = np.where(fp_a, tp, 0.0), np.where(fp_a, tpp, 0.0), np.where(fp_a, tmp, 0.0)
            tn, tpn, tmn = np.where(fn_a, tn, 0.0), np.where(fn_a, tpn, 0.0), np.where(fn_a, tmn, 0.0)
            tg, tp, tn = _pin_v(
                tg, tp, tn, tpg, tpp, tpn
            )  # the message is a claim about THIS node's mass
            # ── the COMPOSITION half of the cliff cost: the DL mismatch deflation, on the PINNED densities.
            # Pinning first is what makes the gap a pure COMPOSITION statement: `_pin_v` has already rescaled the
            # message to this node's own mass, so the common scale (the reframe residual) is gone from G and only
            # the share drift is left. Every stream is deflated — the anchor recovers only when the composition
            # τ-stream is damped alongside the measurement one (ablation-confirmed, HANDOFF_4 §6).
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
            ttau = np.where((tg > _EPS) & ((tp + tn) > _EPS), ttau, 0.0)
# (intron λ gate removed — see if EXP-E alone suffices)
            if not _s2t_off:
                g_g, c_g = mismatch_gap(tg, og)
                g_p, c_p = mismatch_gap(tp, op)
                g_n, c_n = mismatch_gap(tn, on)
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
                g_R, c_R = mismatch_gap(tp + tn, op + on)
                _tau_pre = ttau
                ttau = mismatch_deflate(ttau, g_g - g_R, c_g | c_R, v_own_lam)
                if _capture is not None:  # inert: the per-message gaps + the τ-stream kill, for the dissect loop
                    _capture.setdefault("_dl", []).append(
                        {
                            "df": df,
                            "G_g": g_g.copy(),
                            "G_p": g_p.copy(),
                            "G_n": g_n.copy(),
                            "G_lam": (g_g - g_R),
                            "contra": (c_g | c_R).copy(),
                            "tau_pre": _tau_pre.copy(),
                            "tau_post": ttau.copy(),
                        }
                    )
            return tg, tp, tn, tpg, tpp, tpn, tmg, tmp, tmn, ttau

        def _fuse_add(a, b):  # additive (inverse-variance) fuse of two independent precision streams
            return np.asarray(a, np.float64) + np.asarray(b, np.float64)

        def _fuse_v(a, pa, b, pb):
            p = pa + pb
            return np.where(p > _EPS, (pa * a + pb * b) / np.maximum(p, _EPS), 0.0), p

        dc_fin = None
        f_cur = np.asarray(f_g, np.float64).copy()
        for _ in range(_RHO_ITERS):
            _, rho_lf, rho_rf = _rho_faces(f_cur)
            ag, ap, an, apg, app, apn, amg, amp, amn, atau = _transport(
                sl, vl, 0, 1, fwd, rho_lf, rho_rf
            )  # left msg: dst face 0, src face 1
            bg, bp, bn, bpg, bpp, bpn, bmg, bmp, bmn, btau = _transport(
                sr, vr, 1, 0, bwd, rho_rf, rho_lf
            )  # right msg: dst face 1, src face 0
            cg, cpg = _fuse_v(ag, apg, bg, bpg)  # density MODE (full precision-weighted)
            cp, cpp = _fuse_v(ap, app, bp, bpp)
            cn, cpn = _fuse_v(an, apn, bn, bpn)
            # measurement + composition precisions: additive (independent left + right evidence)
            cm_g, cm_p, cm_n = _fuse_add(amg, bmg), _fuse_add(amp, bmp), _fuse_add(amn, bmn)
            c_tau = _fuse_add(atau, btau)
            mo_g = np.log(np.maximum(cg * E_g / np.maximum(M, _EPS), _EPS))
            mo_p = np.log(np.maximum(cp * E_r / np.maximum(M, _EPS), _EPS))
            mo_n = np.log(np.maximum(cn * E_r / np.maximum(M, _EPS), _EPS))
            # ── THE THREE-STREAM SINGLE-λ COMBINE (the M6 rank-1 fix, message_variance_derivation.md §4) ──
            # (1) COMPOSITION → ONE λ-message on λ = log(f_g/f_R): mode ``mo_g − mo_R`` (from the density relay),
            #     precision ``c_tau`` (the fused Schur τ) — ψ counts the composition DOF ONCE, not twice.
            cR = cp + cn
            mo_R = np.log(np.maximum(cR * E_r / np.maximum(M, _EPS), _EPS))
            lam_msg = mo_g - mo_R
            # …and the same structural gate on the FUSED pair, which is what this mode is actually built from.
            # `_transport`'s per-message gate cannot catch every case: a message may carry an RNA DENSITY while
            # contributing zero mode-fusion PRECISION, in which case `cR` collapses to 0 here and `mo_R` hits
            # the floor again. A λ message exists only when BOTH components of the pair reached this node.
            c_tau = np.where((cg > _EPS) & (cR > _EPS), c_tau, 0.0)
            # (2) ANCHOR gDNA MEASUREMENT → gdna_imp (mode mo_g, precision ``cm_g``). (3) SPLICED RNA MEASUREMENT
            #     → rna_imp (mode mo_p/mo_n, precision ``cm_p``/``cm_n``). INDEPENDENT of the composition, so
            #     fused separately (an RNA-only spliced measurement constrains f_g via f_R with NO gDNA info).
            # θ TILT (AMBIG): τ_tilt = (c_p − c_n)/(c_p + c_n) → arcsin; precision the summed measured RNA.
            tau_tilt = np.clip(np.where(cR > _EPS, (cp - cn) / np.maximum(cR, _EPS), 0.0), -1.0, 1.0)
            th_msg = np.arcsin(tau_tilt)
            th_prec = np.where(is_amb, cm_p + cm_n, 0.0)
            dc_fin = _local_solve(
                global_lp, mo_g, cm_g, (mo_p, mo_n), (cm_p, cm_n),
                lam_imp=(lam_msg, c_tau), theta_imp=(th_msg, th_prec),
            )
            # ⚠ KNOWN DEFECT, deliberately NOT fixed here (it is orthogonal to this term, and the fix REGRESSES
            # the gated arm — decide it on its own evidence, not as a rider). ψ's output at an UNSOLVABLE node
            # (no free RNA strand: an intergenic region, a TSS/TES seam) is discarded by the write-back gate,
            # but it still feeds the NEXT iteration's reframe here — and the solver returns 0 for a node it
            # never solved, so every gDNA anchor is re-framed as if it were 100 % RNA (ρ_tot off by E_g/E_r,
            # up to 1.8×, on every edge incident to an anchor). The one-line fix is
            # ``np.where(solvable, ..., f_g)``; measured A/B (arm ``dl4``): aggregate 0.0969→0.0972 (refit=0),
            # 0.0828→0.0832 (refit=1), unstranded-capON 0.1702→0.1740 — i.e. the CORRECT arithmetic scores
            # slightly worse, which means something downstream is compensating for it. Investigate the pair.
            f_cur = np.clip(np.asarray(dc_fin.gdna_frac, np.float64), 0.0, 1.0)
            nonlocal _uni_msg
            _uni_msg = (mo_g, cpg, mo_p, cpp, mo_n, cpn)  # publish for the shared diagnostics
            if (
                _capture is not None
            ):  # inert diagnostic: the fused per-component densities + the frames
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
                        "rho_lf": rho_lf.copy(),
                        "rho_rf": rho_rf.copy(),
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
                        "fg_out": f_cur.copy(),
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
                    "rho_l0": rho_l0,
                    "rho_r0": rho_r0,
                    # ── AUDIT_2 instrumentation (invariant scan) ──
                    "order": np.asarray(order_list, np.int64),
                    "logvar_tot": logvar_tot.copy(),
                    "SP_l": SP[0].copy(), "SP_r": SP[1].copy(),
                    "SN_l": SN[0].copy(), "SN_r": SN[1].copy(),
                    "n_unspl_l": np.asarray(geometry.n_unspl_left, np.float64),
                    "n_unspl_r": np.asarray(geometry.n_unspl_right, np.float64),
                    "spl_n_pos_l": np.asarray(geometry.spliced_n_pos_left, np.float64),
                    "spl_n_pos_r": np.asarray(geometry.spliced_n_pos_right, np.float64),
                    "spl_n_neg_l": np.asarray(geometry.spliced_n_neg_left, np.float64),
                    "spl_n_neg_r": np.asarray(geometry.spliced_n_neg_right, np.float64),
                    "tau_own": tau_own.copy(), "mg_own": mg_own.copy(),
                    "struct_lock": _struct.copy(),
                    "fwd_mp": fwd[7], "bwd_mp": bwd[7],
                    "fwd_mn": fwd[8], "bwd_mn": bwd[8],
                    "fwd_mg": fwd[6], "bwd_mg": bwd[6],
                    "fwd_tau": fwd[9], "bwd_tau": bwd[9],
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
    # EMPIRICALLY REFUTED (solve_gate_design.md): it regresses both standalone (refit=0 +0.010) and with the
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
            statics.u_pos,
            statics.u_neg,
            fp,
            fn,
            statics.mass_unspliced,
            statics.mass_spliced,
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
            vn_loc=_dc_loc.rna_neg_frac_var,
            a_fwd=None,
            b_bwd=None,
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
            mass_l=MS[0],
            mass_r=MS[1],
            spl_l=SP[0] + SN[0],
            spl_r=SP[1] + SN[1],
            free_pos=np.asarray(fp, bool),
            free_neg=np.asarray(fn, bool),
            # global geometry (μ = clip(ρ·eff_global/mass_global, 0, 1) is the implied prior fraction).
            eff_global=eff_global,
            mass_global=mass_global,
            # per-face geometry for message dissection (logodds diagnostics)
            eff_gdna_l=EG[0],
            eff_gdna_r=EG[1],
            eff_rna_l=ER[0],
            eff_rna_r=ER[1],
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


def chain_region_deconv(chain: NodeChain, belief: NodeBelief, substrate) -> NodeDeconv:
    """Project the chain belief's REGION nodes back to a region-keyed :class:`NodeDeconv` — the transitional
    region projection the existing ``CalibrationResult`` / ``priors`` / ``derive`` consume (the per-node
    first-class schema rewire is P4). gDNA / RNA masses from each region's solved ``f_g`` over its contained
    unspliced (+ the always-RNA contained spliced) mass."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    reg = kind == REGION
    mass_u = np.asarray(substrate.contained.mass_unspliced, dtype=np.float64)
    mass_s = np.asarray(substrate.contained.mass_spliced, dtype=np.float64)
    R = mass_u.shape[0]
    f_g = np.zeros(R)
    f_pos = np.zeros(R)
    f_neg = np.zeros(R)
    ri = idx[reg]
    f_g[ri] = belief.f_g[reg]
    f_pos[ri] = belief.f_pos[reg]
    f_neg[ri] = belief.f_neg[reg]
    return NodeDeconv(
        gdna_mass=f_g * mass_u,
        rna_mass=(1.0 - f_g) * mass_u + mass_s,
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
    )


def chain_boundary_side_deconv(chain: NodeChain, belief: NodeBelief, substrate):
    """Project the chain belief's BOUNDARY ``f_g`` onto each region's two SIDE views — the boundary-flux that
    ``priors``' pooled-seam gDNA eff-len + ``derive`` consume.

    Region ``r``'s left/right boundary IS its left/right chain neighbour; that boundary's solved ``f_g`` splits
    ``r``'s side crossing mass (``substrate.left[r]`` / ``substrate.right[r]`` — the boundary flux already
    projected onto the region by the D1 side-attribution) into gDNA / RNA (the RNA spliced-inclusive, matching
    the contained projection). One boundary pie applied to its side mass. Returns ``(left, right)`` region-keyed
    :class:`NodeDeconv`."""
    kind = np.asarray(chain.kind)
    reg_nodes = np.asarray(chain.order)[kind == REGION]
    ridx = np.asarray(chain.ref_idx, dtype=np.int64)[reg_nodes]
    R = int(ridx.max()) + 1 if ridx.size else 0
    fg = np.asarray(belief.f_g, dtype=np.float64)
    left_fg = np.zeros(R)
    right_fg = np.zeros(R)
    left_fg[ridx] = fg[np.asarray(chain.left)[reg_nodes]]  # f_g of r's left-flank boundary node
    right_fg[ridx] = fg[np.asarray(chain.right)[reg_nodes]]

    def _side(side_fg, view):
        m_u = np.asarray(view.mass_unspliced, dtype=np.float64)
        m_s = np.asarray(view.mass_spliced, dtype=np.float64)
        return NodeDeconv(
            gdna_mass=side_fg * m_u,
            rna_mass=(1.0 - side_fg) * m_u + m_s,
            gdna_frac=side_fg,
            rna_pos_frac=np.zeros(R),
            rna_neg_frac=np.zeros(R),
        )

    return _side(left_fg, substrate.left), _side(right_fg, substrate.right)
