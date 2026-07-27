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

import math

import numpy as np
from scipy.special import zeta

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

# ⭐ ONE ρ-iteration. The combine's reframe needs the DESTINATION's composition, which is what the messages
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
_RHO_ITERS = 1


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
    _, solve_grid = _logodds_grid(int(n_grid), float(logodds_window))
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
    # (docs/calibration/archive/CALIBRATION_MASTER.md §5): this ``gdna_prior`` is the COMPOSITION arm ONLY.
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
        accept_l = (
            SP[0] + SN[0]
        ) > _EPS  # the LEFT face carries the spliced (acceptor) ⇒ WITH-spliced ρ_tot
        accept_r = (SP[1] + SN[1]) > _EPS

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

# ⚠ The two seams must be compared IN THE DESTINATION'S OWN FRAME. Each is lifted by its own
        # enrichment step ``r``, and under capture those steps differ (the two seams sit at different probe
        # positions), so comparing unlifted fluxes reads a capture difference as an abundance difference.
        # Measured on the MODE variant of this idea (replace the claim by the dominant seam's flux): raw is
        # uniformly better off-capture (14/14 conditions, 0 worse) and uniformly worse on (8 worse); lifting
        # first halves the on-capture damage. That MODE fix is NOT landed — under capture the graft already
        # OVER-states (median φ = 2.45, an M8 frame problem), so no bound-tightening can help there — but the
        # same seam pair is what measures the variance below, and it must be framed for that too.
        def _flank_dom(rho_lface, rho_rface, spf):
            """Per-node: the flux each of its two flanking BOUNDARIES sends it, ALREADY lifted into this
            node's frame, per strand.

            ``spf[1]`` is a node's right face — what a LEFT neighbour presents to us — and ``spf[0]`` its
            left face, what a RIGHT neighbour presents; the matching frame steps are the same pairs the
            relay and the combine form. Zero wherever a flank is absent or is not a boundary, so on a
            one-junction exon this degenerates exactly to that junction's own lifted flux."""

            def _side(ok, nb, src_face, dst_face, face):
                r = np.where(
                    ok & (src_face[nb] > _EPS) & (dst_face > _EPS),
                    dst_face / np.maximum(src_face[nb], _EPS),
                    1.0,
                )
                return np.where(ok & is_bnd_a[nb], spf[face][nb] * r, 0.0)

            return (
                _side(_vl_a, _sl_a, rho_rface, rho_lface, 1),
                _side(_vr_a, _sr_a, rho_lface, rho_rface, 0),
            )


        # P1d is ON by default; `RIGEL_GLV_OFF=1` ablates it (bit-identical to the pre-P1d path).
        vgp_prem = vgn_prem = None  # per-ρ-iteration, in the destination's frame; set before each transport

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
            """Scale a claim to the node's observed mass: `Σ_c ρ_c·E_c = M`, an IDENTITY under the imputation
            premise, restored by a common factor `k = M/S`.

            The partial-claim semantics are load-bearing: a component the claim does not SUPPLY (precision 0)
            contributes the node's OWN density to the mass budget and does not move, which is what keeps a
            partial claim partial — a message carrying gDNA only still delivers `f_g < 1`. Rescaling all
            three blindly instead regresses capture-OFF 3.6×.

            The weighted alternative (M12 — apportion the correction
            by how badly each component is known, of which this common factor is the `w → 0` limit) was
            derived, implemented and A/B'd, and is NOT used: once λ is fused by its own precision (the
            message packet) it is completely INERT — 0 better / 0 worse / 32 flat — and the per-message
            variant is a net loss. Recorded in `weighted_rescale_design.md` §9 and `PASS0_FINISH_PLAN.md`."""
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
        _mu_f = ((spl_p_f[0], spl_n_f[0]), (spl_p_f[1], spl_n_f[1]))
        _v_mu_f = tuple(
            tuple(np.where(c > 0.0, 1.0 / np.maximum(c, _EPS), np.inf) for c in cs)
            for cs in (
                (
                    np.asarray(geometry.spliced_n_pos_left, np.float64),
                    np.asarray(geometry.spliced_n_neg_left, np.float64),
                ),
                (
                    np.asarray(geometry.spliced_n_pos_right, np.float64),
                    np.asarray(geometry.spliced_n_neg_right, np.float64),
                ),
            )
        )

        # ── ⭐ THE SEQUENTIAL SCAN'S OPERANDS, AS PYTHON LISTS (`*_l`) ──────────────────────────────────
        # `_relay` reads every one of these ONE ELEMENT AT A TIME — it is a Gauss-Seidel scan and cannot be
        # vectorised — and at genome scale that is ~6 M edge-iterations, ~40 reads each. `.tolist()` is
        # exact on float64 / int64 / bool (identical IEEE-754 doubles, identical ints), so this is
        # BIT-IDENTICAL by construction, and it buys ~3×: `lst[i]` costs a third of `arr[i]`, and it yields
        # a PYTHON float, whose arithmetic is ~3× `np.float64`'s because no intermediate is boxed in a 0-d
        # array. The array forms stay — the vectorised combine (`_transport`, `_peel_share`) needs them.
        (og_l, op_l, on_l, pg_own_l, pp_own_l, pn_own_l, mg_own_l, mp_own_l, mn_own_l, tau_own_l,
         M_l, E_g_l, E_r_l, n_node_l, logvar_l) = (
            np.asarray(a, np.float64).tolist()
            for a in (og, op, on, pg_own, pp_own, pn_own, mg_own, mp_own, mn_own, tau_own,
                      M, E_g, E_r, _n_node, logvar_tot)
        )
        ex_l, bnd_l, fp_l, fn_l = (a.tolist() for a in (ex_a, is_bnd_a, fp_a, fn_a))
        spl_p_l, spl_n_l, SP_l, SN_l = (
            [f.tolist() for f in t] for t in (spl_p_f, spl_n_f, SP, SN)
        )
        mu_l, v_mu_l = ([[c.tolist() for c in face] for face in t] for t in (_mu_f, _v_mu_f))

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
        # ⚠ AND `_far` WAS NOT THE LAST THIRD-NODE DEPENDENCE. `_RHO_ITERS = 2`, and the second iteration's
        # reframe faces are built from the destination's FUSED posterior — i.e. from the other message —
        # so the reframe itself still reuses it (measured median |Δlog ρ_face| between iterations 0.0116
        # stranded-capOFF to 0.1242 unstranded-capON, >1 % on 52.7–79.0 % of nodes). Deleting `_far` does
        # NOT make the message BP-legal; it removes the largest and most direct of two violations. The
        # remaining one is recorded in `PASS0_FINISH_PLAN.md` as P4b.
        #
        # The legal way to recover what `_far` supplied, if it is ever wanted back: carry the peel share as
        # a FUNCTION of the destination's state (a proper pairwise potential ψ(x_L, x_k)) instead of a
        # plugged-in point estimate, and let the destination's own ψ solve modulate it. No data reuse at
        # all. Structurally available, unimplemented, unmeasured.

        def _peel_share(df, tg, tpg, tp, tn):
            """The continuing share ``w`` and ``Var(log w)`` per strand, on face ``df`` of every node, for a
            message whose gDNA claim is ``(tg, tpg)`` and whose RNA claim is ``(tp, tn)``. Returns
            ``((w_p, vw_p), (w_n, vw_n))``; ``Var(log w) = +inf`` (⇒ zero precision, an inert message) only
            where NONE of the three estimators of the level exists.

            ⚠ TWIN of `_peel_share_scalar` (the sequential relay's arm) — mirror any change into both."""
            _vg = np.where(np.asarray(tpg, np.float64) > 0.0, 1.0 / np.maximum(tpg, _EPS), np.inf)
            _nu_m, _, _vl_m = residual_level(M, _n_node, tg, E_g, E_r, _vg)
            _A = np.asarray(tp, np.float64) + np.asarray(tn, np.float64)
            _a_p = np.where(_A > _EPS, np.asarray(tp, np.float64) / np.maximum(_A, _EPS), 0.0)
            out = []
            for _a, _mu, _vmu in (
                (_a_p, _mu_f[df][0], _v_mu_f[df][0]),
                (1.0 - _a_p, _mu_f[df][1], _v_mu_f[df][1]),
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
                _kk = np.maximum(_nu * _nu * _pt, 1.0)
                _v_nu = np.where(_live, zeta(2.0, _kk), np.inf)
                _w = np.where(_live, peel_continue_share(_nu, _mu), 0.0)
                if _capture is not None:
                    _capture.setdefault("_lvl", []).append(  # inert: the level's provenance, per face
                        {"df": df, "nu": np.asarray(_nu).copy(), "v_nu": np.asarray(_v_nu).copy(),
                         "pm": np.asarray(_pm).copy(),
                         "nu_m": np.asarray(_nu_ms).copy(),
                         "mu": np.asarray(_mu).copy(), "w": np.asarray(_w).copy(),
                         "v_g": np.asarray(_vg).copy(),
                         "vl_m": np.where(np.isfinite(_vl_m), _vl_s, np.inf),
                         "phi": np.asarray(tg * E_g / np.maximum(M, _EPS)).copy()}
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

        def _peel_share_scalar(i, df, tg, tpg, tp, tn):
            """The SCALAR twin of `_peel_share`, for one node ``i`` — see that docstring for the model.

            ⚠ TWIN: mirror any change into both. It exists because `_relay` is a sequential Gauss-Seidel
            scan (it cannot be vectorised) and calls this once per node per direction, so the array form
            runs ~50 numpy ops on 0-d arrays per call — 0.5–0.7 µs each against ~0.02 µs for the float
            expression. Same arithmetic in the same association order; every `np.where` becomes the branch
            it encodes, so the dead arms (`ζ(2,·)` and `residual_level`'s four `log_ndtr`s at a node with
            no level) are no longer evaluated and discarded. Measured 25× on this path, which is the bulk
            of a genome-scale calibration."""
            _vg = 1.0 / _fmax(tpg, _EPS) if tpg > 0.0 else math.inf
            _nu_m, _, _vl_m = residual_level_scalar(M_l[i], n_node_l[i], tg, E_g_l[i], E_r_l[i], _vg)
            _fin = math.isfinite(_vl_m)
            _A = tp + tn
            _a_p = tp / _A if _A > _EPS else 0.0
            out = []
            for _a, _mu, _vmu in (
                (_a_p, mu_l[df][0][i], v_mu_l[df][0][i]),
                (1.0 - _a_p, mu_l[df][1][i], v_mu_l[df][1][i]),
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
                _nu = (_pm * _nu_ms) / _pm  # ⚠ NOT `_nu_ms` — `(a·b)/a ≠ b` in floating point
                _v_nu = float(zeta(2.0, _fmax(_nu * _nu * _pm, 1.0)))
                _w = peel_continue_share_scalar(_nu, _mu)
                _wm = 1.0 - _w
                # a spliced DENSITY with no spliced COUNT cannot be priced (plan §4.7) ⇒ no claim at all.
                _ok = math.isfinite(_vmu) or not _mu > _EPS
                out.append(
                    (_w, _wm * _wm * (_v_nu + (_vmu if math.isfinite(_vmu) else 0.0)) if _ok else math.inf)
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
        # They implement the SAME eight-step per-edge transform in the same order — reframe `r` → detect the
        # graft → σ²_transfer (M5) → graft the source's spliced → damp the three precision streams → graft
        # precision ⊕ M8 ⊕ P1d → peel by composition (M10) → strand filter → pin to M — and the correspondence
        # is line-for-line. Merging them into one polymorphic routine (`k` = a scalar index for the relay,
        # `slice(None)` for the combine) IS structurally possible; `_peel_share` already proves the pattern.
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
        def _relay(seq, nbr, dst_face, src_face, df, sf):
            # every operand here is a Python float or bool — see the `*_l` block above
            rg, rp, rn = og_l.copy(), op_l.copy(), on_l.copy()
            pg, pp, pn = pg_own_l.copy(), pp_own_l.copy(), pn_own_l.copy()  # full → MODE fusion
            mg, mp, mn = mg_own_l.copy(), mp_own_l.copy(), mn_own_l.copy()  # MEASUREMENT (anchor + spliced)
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
                rho_src = src_face[s]
                r = (
                    (dst_face[i] / rho_src) if (rho_src > _EPS and dst_face[i] > _EPS) else 1.0
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
                gp = spl_p_l[sf][s] if _gr else 0.0
                gn = spl_n_l[sf][s] if _gr else 0.0
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
                    _s2f = s2t + graft_frame_logvar_scalar(r)
                    _sps = SP_l[sf][s]
                    _spc = _sps / (1.0 + _sps * _s2f) if _sps > _EPS else 0.0
                    _sns = SN_l[sf][s]
                    _snc = _sns / (1.0 + _sns * _s2f) if _sns > _EPS else 0.0
                    tpp += _spc
                    tpn += _snc
                    tmp += _spc
                    tmn += _snc
                    _vgp, _vgn = vgp_l[i], vgn_l[i]
                    tpp, tmp = _damp_v(tpp, _vgp), _damp_v(tmp, _vgp)
                    tpn, tmn = _damp_v(tpn, _vgn), _damp_v(tmn, _vgn)

                if bnd_l[i] and ex_l[s]:  # EXON → boundary: PEEL by COMPOSITION (scale by the share)
                    (_wp, _vwp), (_wn, _vwn) = _peel_share_scalar(i, df, tg, tpg, tp, tn)
                    tp, tn = tp * _wp, tn * _wn
                    tpp, tmp = _damp_v(tpp, _vwp), _damp_v(tmp, _vwp)
                    tpn, tmn = _damp_v(tpn, _vwn), _damp_v(tmn, _vwn)
                if not fp_l[i]:
                    tp, tpp, tmp = 0.0, 0.0, 0.0
                if not fn_l[i]:
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
                np.asarray(a, np.float64)
                for a in (rg, rp, rn, pg, pp, pn, mg, mp, mn, tau)
            )

        # dst faces its source on its LEFT (face 0); the source faces the dst on its RIGHT (face 1) — and mirrored
        # for the backward pass. The face pair selects BOTH the per-side ρ_tot and the per-face mature flux.
        def _seam_pair(rho_lface, rho_rface):
            """Per strand: the graft's premise log-variance — ONE library-level scalar, fitted by method of
            moments from the destination-frame disagreement of exons' flanking seam PAIRS
            (`graft_premise_logvar`) and applied to every graft edge.

            ⚠⚠ **A DEBT, not a model.** The one scalar stands in for a quantity that splits ≥30× on whether
            the boundary carries a transcript TERMINUS — a bit the region map does not have. Re-derive this
            per structural class when TSS/TES land (P1g). See `graft_premise_logvar`."""
            out = []
            for spf, vmu in ((spl_p_f, 0), (spl_n_f, 1)):
                fl, fr = _flank_dom(rho_lface, rho_rface, spf)
                # each seam's own noise: its spliced COUNT (never the mass) ⊕ its lift's scale sampling
                # (M5's source leg; the destination's leg is common to both lifts and cancels in ``d``).
                _lv = np.where(np.isfinite(logvar_tot), logvar_tot, 0.0)
                per, pooled = graft_premise_logvar(
                    fl,
                    fr,
                    np.where(_vl_a, _v_mu_f[1][vmu][_sl_a] + _lv[_sl_a], np.inf),
                    np.where(_vr_a, _v_mu_f[0][vmu][_sr_a] + _lv[_sr_a], np.inf),
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
                if _capture is not None:  # inert: the fitted scalar and the population it was fitted on
                    _ok = (fl > _EPS) & (fr > _EPS)
                    _d = np.log(np.maximum(fl, _EPS)) - np.log(np.maximum(fr, _EPS))
                    _vv = np.where(_vl_a, _v_mu_f[1][vmu][_sl_a] + _lv[_sl_a], 0.0) + np.where(
                        _vr_a, _v_mu_f[0][vmu][_sr_a] + _lv[_sr_a], 0.0
                    )
                    _capture.setdefault("_glv", []).append(
                        {"strand": vmu, "omega": pooled, "n_pairs": int(_ok.sum()),
                         "ok": _ok.copy(), "d": _d.copy(), "noise": _vv.copy(),
                         "Ed2": float((_d[_ok] ** 2).mean()) if _ok.any() else 0.0,
                         "Enoise": float(_vv[_ok].mean()) if _ok.any() else 0.0}
                    )
                out.append(np.full_like(per, pooled))
            return out[0], out[1]
        # the relay runs on the INPUT-belief faces, so its seam pair is formed from those
        vgp_prem, vgn_prem = _seam_pair(rho_l0, rho_r0)
        vgp_l, vgn_l = vgp_prem.tolist(), vgn_prem.tolist()
        left_l, right_l = left.tolist(), right.tolist()
        rho_l0_l, rho_r0_l = rho_l0.tolist(), rho_r0.tolist()
        fwd = _relay(order_list, left_l, rho_l0_l, rho_r0_l, 0, 1)
        bwd = _relay(order_list[::-1], right_l, rho_r0_l, rho_l0_l, 1, 0)
        # ── the COMBINE: transport α (from left neighbour) + β (from right neighbour) into the node's frame with
        # the LAZY ρ_tot (two-iteration — the 2nd uses the both-message composition), fuse, ÷M_dst → the ψ solve.
        li, ri, vl, vr, sl, sr = _li_a, _ri_a, _vl_a, _vr_a, _sl_a, _sr_a

        # The VECTORISED twin of `_relay` — see the DO-NOT-MERGE note there, which applies to both.
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
                transfer_logvar(logvar_tot, logvar_tot[src], graft)
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
            _s2t_spl = _s2t_spl + np.where(graft, graft_frame_logvar(r), 0.0)
            _spc = np.where(_sp > _EPS, _sp / (1.0 + _sp * _s2t_spl), 0.0)
            _snc = np.where(_sn > _EPS, _sn / (1.0 + _sn * _s2t_spl), 0.0)
            tpp, tpn = tpp + _spc, tpn + _snc  # into the mode-fusion precision …
            tmp, tmn = tmp + _spc, tmn + _snc  # … and the measurement stream (a count, never composition τ)
            # ⭐ P1d — the graft's PREMISE variance (`graft_premise_logvar`), applied to the WHOLE RNA
            # claim after the spliced arm is folded in, because the premise is about the SUM: measured
            # FLAT in the spliced share w_μ (Var 2.02 → 1.83 across w_μ 0.47 → 1.00, while Var/w_μ²
            # swings 5×), so charging the spliced arm alone would reach only 10–93 % of the delivered
            # confidence while the error contaminates 63–95 % of the delivered density.
            _vgp = np.where(graft, vgp_prem, 0.0)
            _vgn = np.where(graft, vgn_prem, 0.0)
            tpp, tmp = tpp / (1.0 + tpp * _vgp), tmp / (1.0 + tmp * _vgp)
            tpn, tmn = tpn / (1.0 + tpn * _vgn), tmn / (1.0 + tmn * _vgn)

            peel = is_bnd_a & ex_a[src] & valid  # EXON → boundary: PEEL by COMPOSITION (the relay's twin)
            (_wp, _vwp), (_wn, _vwn) = _peel_share(df, tg, tpg, tp, tn)
            tp = np.where(peel, tp * _wp, tp)
            tn = np.where(peel, tn * _wn, tn)

            def _dv_arr(pr, vv):
                _f = np.isfinite(vv)
                return np.where(peel, np.where(_f, pr / (1.0 + pr * np.where(_f, vv, 0.0)), 0.0), pr)

            tpp, tmp = _dv_arr(tpp, _vwp), _dv_arr(tmp, _vwp)
            tpn, tmn = _dv_arr(tpn, _vwn), _dv_arr(tmn, _vwn)
            tp, tpp, tmp = np.where(fp_a, tp, 0.0), np.where(fp_a, tpp, 0.0), np.where(fp_a, tmp, 0.0)
            tn, tpn, tmn = np.where(fn_a, tn, 0.0), np.where(fn_a, tpn, 0.0), np.where(fn_a, tmn, 0.0)
            if _capture is not None:  # inert: the PRE-PIN state, for the weighted-rescale prototype
                _capture.setdefault("_pin", []).append(
                    {
                        "df": df, "src": np.asarray(src).copy(), "valid": np.asarray(valid).copy(),
                        "tg": tg.copy(), "tp": tp.copy(), "tn": tn.copy(),
                        "tpg": tpg.copy(), "tpp": tpp.copy(), "tpn": tpn.copy(),
                        # the COMMON-mode variance every component shares: the reframe's scale (M5) plus the
                        # source's own count. Everything else in a component's variance is its own.
                        "s2t": np.where(np.isfinite(s2t), s2t, 0.0),
                        "n_src": np.asarray(_n_node)[src].copy(),
                        # the graft: how much of the RNA claim is a MEASUREMENT rather than an imputation
                        "spl_p": _sp.copy(), "spl_n": _sn.copy(), "spl_prec": (_spc + _snc).copy(),
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
            # **When the bias strata are diagnosed, this term must SHRINK.** See `variance_ledger.md` §6.

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
            # implemented, MC-refuted and DELETED — see `variance_ledger.md` §6; do not rebuild it.
            _dg = _dp = _dn = np.where(_sup.any(axis=-1), _b2, 0.0)
            tpg, tmg = tpg / (1.0 + tpg * _dg), tmg / (1.0 + tmg * _dg)
            tpp, tmp = tpp / (1.0 + tpp * _dp), tmp / (1.0 + tmp * _dp)
            tpn, tmn = tpn / (1.0 + tpn * _dn), tmn / (1.0 + tmn * _dn)
            tg, tp, tn = _pin_v(tg, tp, tn, tpg, tpp, tpn)  # a claim about THIS node's mass
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
            # The gate asks whether the source SUPPLIED both components of the λ pair. "Supplied" is a
            # statement about PRECISION, not about the density's value: a component carrying zero precision is
            # not supplied however large its density, and a component at zero density with live precision IS
            # supplied — it is the claim "there is none of this here", which is exactly a composition claim.
            # Testing the density conflates the two and silences λ wherever a legitimate zero is emitted.
            ttau = np.where((tpg > 0.0) & ((tpp + tpn) > 0.0), ttau, 0.0)
# (intron λ gate removed — see if EXP-E alone suffices)
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

        def _fuse_add(a, b):  # additive (inverse-variance) fuse of two independent precision streams
            return np.asarray(a, np.float64) + np.asarray(b, np.float64)

        def _fuse_v(a, pa, b, pb):
            p = pa + pb
            return np.where(p > _EPS, (pa * a + pb * b) / np.maximum(p, _EPS), 0.0), p

        # ONE ρ-iteration (see `_RHO_ITERS`), so this is straight-line: the frames are the INPUT-belief
        # faces `rho_l0`/`rho_r0` already built above, and there is no next iteration to feed.
        ag, ap, an, apg, app, apn, amg, amp, amn, atau, alam, ath = _transport(
            sl, vl, 0, 1, fwd, rho_l0, rho_r0
        )  # left msg: dst face 0, src face 1
        bg, bp, bn, bpg, bpp, bpn, bmg, bmp, bmn, btau, blam, bth = _transport(
            sr, vr, 1, 0, bwd, rho_r0, rho_l0
        )  # right msg: dst face 1, src face 0
        cg, cpg = _fuse_v(ag, apg, bg, bpg)  # density MODE (full precision-weighted)
        cp, cpp = _fuse_v(ap, app, bp, bpp)
        cn, cpn = _fuse_v(an, apn, bn, bpn)
        cm_g, cm_p, cm_n = _fuse_add(amg, bmg), _fuse_add(amp, bmp), _fuse_add(amn, bmn)
        c_tau = _fuse_add(atau, btau)
        mo_g = np.log(np.maximum(cg * E_g / np.maximum(M, _EPS), _EPS))
        mo_p = np.log(np.maximum(cp * E_r / np.maximum(M, _EPS), _EPS))
        mo_n = np.log(np.maximum(cn * E_r / np.maximum(M, _EPS), _EPS))
        # ── THE THREE-STREAM SINGLE-λ COMBINE (the M6 rank-1 fix, message_variance_derivation.md §4) ──
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
        # DIAGNOSTIC ONLY (P1, `docs/calibration/PASS0_FINISH_PLAN.md`): ablate the RNA MEASUREMENT ψ
        # factor. It is the channel that carries 75 % of the posterior precision on the confidently-wrong
        # unstranded × capture-OFF exons — and ablating it is NOT the fix: 0.0895 → 0.1033, 4 better / 17
        # worse, because it is also the only thing that lets a zero-gDNA library say "my mass is all RNA"
        # (`gdna_none` 0.1063 → 0.1438). Kept for developing the real fix.
        dc_fin = _local_solve(
            global_lp, mo_g, cm_g, (mo_p, mo_n), (cm_p, cm_n),
            lam_imp=(lam_msg, c_tau), theta_imp=(th_msg, th_prec),
        )
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
                    "rho_lf": rho_l0.copy(),
                    "rho_rf": rho_r0.copy(),
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
    # EMPIRICALLY REFUTED (docs/calibration/archive/solve_gate_design.md): it regresses both standalone (refit=0 +0.010) and with the
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
