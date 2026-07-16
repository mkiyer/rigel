"""The belief-propagation calibration solver: per-node gDNA-vs-RNA deconvolution by a bidirectional sweep.

Deconvolves each node's unspliced fragment mass into a pie ``(f_pos, f_neg, f_g)`` — sense-RNA /
antisense-RNA / gDNA — over the unified region↔boundary chain (`node_chain`), by a single forward-backward
(L→R then R→L) belief-propagation pass (exact on the chain, a forest of linear paths). Each per-node solve
(`simplex_logodds`, the log-density log-odds solver) reconciles three sources of information: the intrinsic
strand likelihood (the Beta-Binomial tilt — the only signal a count carries), the cross-node imputation
messages, and the population gDNA prior. Theory: the count-zero-information principle in
`docs/calibration/CALIBRATION_ARCHITECTURE.md`; the ψ reference in
`docs/calibration/reference_prior_derivation.md`; the message-precision model — the source's own honest
belief precision ``Var(log ρ_c^src) = Var(log f_c^src) + 1/n_c^src`` — in the `_scan` docstring below.

The gDNA population prior is the pass-0 count-space **`GdnaRatePrior`** (`gdna_rate_prior`), fit ONCE before
the sweep on all nodes' total unspliced density and projected onto each node's ψ (`GdnaRatePrior.logprior`) —
extremely weak, so strand + messages dominate. It replaces the retired seed/floor/global + the density KDE.

Module layout. The per-node geometry / belief / statics / init primitives (`build_node_geometry`,
`build_node_statics`, `init_beliefs`, `NodeGeometry`/`NodeBelief`/`NodeStatics`) live in the
lower `node_geometry` module and are re-exported here for the calibrator's convenience; this module owns:
* `node_global_geometry` — the per-node ``(mass, eff)`` support the rate prior is fit on and projected onto.
* `node_sweep` — the single forward-backward sweep. Message precision is the source's own HONEST belief
  precision ``pr = n_src/(n_src·vb_src + 1)`` = ``1/(Var(log f_c^src) + 1/n_src)`` — strand (composition) and
  count (sampling), per component, from the RUNNING belief. It saturates at ``1/vb_src``, so a source can
  never send more confidence than its own belief supports; there is no var~mean fixed point / no outer loop.
  The adjacent-pair overdispersion ``σ²_transfer`` is a PRIOR and is currently ZERO — it returns later as an
  NPMLE projection, not a total-density fit (`adjacent_disagreement_variance`, retained but no longer in the
  solve, is that retired total-density estimator).
* `chain_region_deconv` / `chain_boundary_side_deconv` — project the converged belief back to the per-region
  / per-boundary-side masses the `CalibrationResult` consumes.
"""

from __future__ import annotations

import math

import numpy as np

from .node_chain import BOUNDARY, REGION, NodeChain
from .node_geometry import (
    NodeBelief,
    NodeGeometry,
    NodeStatics,
    build_node_geometry,
    build_node_statics,
    init_beliefs,
)
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
    "chain_region_deconv",
    "chain_boundary_side_deconv",
]

_EPS = 1.0e-9


def _log_sigmoid(x):
    """``log σ(x)`` — stable (never forms ``1−σ``). Vectorized. ``log σ(−x) = _log_sigmoid(−x)``."""
    return -np.logaddexp(0.0, -x)


def _fold_lambda(mu, var, factors, *, L, coarse_k, fine_k, sigma_cov, refine):
    """The coherent-relay fold: Expectation-Propagation moment-match of a Gaussian belief ``N(λ; mu, var)``
    against gDNA / RNA-total message factors, returning ``(μ', σ²')`` (the DOF-pie relay,
    ``docs/calibration/dof_pie_relay_implementation_plan.md`` §5.3/§6.1; prototype ``dof_pie_relay_check.py`` C8).

    ``factors`` = list of ``(is_gdna, mode, prec)``: a factor ``exp(−½·prec·(log f_c(λ) − mode)²)`` on
    ``log f_g = log σ(λ)`` (``is_gdna=True``) or ``log f_r = log σ(−λ)`` (``is_gdna=False``). The message factor
    is not log-concave, so the peak is located by a COARSE grid ``argmax`` (no fragile Newton), then a FINE grid
    re-centered on the peak + its curvature resolves the moments (``refine`` self-correcting iterations). All
    resolution knobs are :class:`CalibrationConfig` fields (numerical, not model). Empty ``factors`` ⇒ unchanged."""
    if not factors:
        return mu, var
    var = max(var, _EPS)
    sig_l = math.sqrt(var)

    def _logpost(lam):
        lp = -0.5 * (lam - mu) ** 2 / var
        for is_g, a, w in factors:
            lc = _log_sigmoid(lam) if is_g else _log_sigmoid(-lam)
            lp = lp - 0.5 * w * (lc - a) ** 2
        return lp

    # coarse bracket: the prior's 6σ window UNION the message targets (each factor's implied λ), clamped to ±L
    lo, hi = mu - sigma_cov * sig_l, mu + sigma_cov * sig_l
    for is_g, a, _w in factors:
        fc = min(math.exp(min(a, 0.0)), 1.0 - 1e-9)  # e^a ∈ (0,1]; a>0 (implied f>1) saturates
        lam_t = math.log(fc / (1.0 - fc))
        lam_t = lam_t if is_g else -lam_t
        lo, hi = min(lo, lam_t), max(hi, lam_t)
    lo, hi = max(-L, lo - 1.0), min(L, hi + 1.0)
    lam_c = np.linspace(lo, hi, coarse_k)
    psi_c = _logpost(lam_c)
    j = int(np.clip(np.argmax(psi_c), 1, coarse_k - 2))
    h = (hi - lo) / (coarse_k - 1)
    curv = -(psi_c[j + 1] - 2.0 * psi_c[j] + psi_c[j - 1]) / (h * h)  # −ψ″ > 0 near a max
    center, sig_hat = float(lam_c[j]), 1.0 / math.sqrt(max(curv, 1e-9))
    m1, v1 = center, sig_hat * sig_hat
    for _ in range(refine):  # re-center + re-width on the moments (converges fast; captures a skewed width)
        half = max(sigma_cov * sig_hat, 1.5 * h)
        flo, fhi = max(-L, center - half), min(L, center + half)
        lam_f = np.linspace(flo, fhi, fine_k)
        lp = _logpost(lam_f)
        lp -= lp.max()
        p = np.exp(lp)
        p /= np.trapezoid(p, lam_f)
        m1 = float(np.trapezoid(lam_f * p, lam_f))
        v1 = max(float(np.trapezoid(lam_f * lam_f * p, lam_f)) - m1 * m1, _EPS)
        center, sig_hat = m1, math.sqrt(v1)
    return m1, v1


def _poisson_moment_var(resid, ns, nd) -> float:
    """The §2 Poisson disagreement-variance moment estimator (the total-density σ²_imp).
    ``Var(d) = σ²_imp + 1/n_i + 1/n_j`` ⇒ subtract each pair's known Poisson
    sampling and inverse-variance-weight-average (``w = nᵢnⱼ/(nᵢ+nⱼ)``, harmonic — 0 at zero count, no
    threshold). ``resid`` must already be oriented to a single sign convention (so the systematic frame
    offset is one mode); it is median-centered here. ``ns``/``nd`` are symmetric in ``w`` and ``v`` so their
    orientation is irrelevant. Returns ``max(·, 0)``; a degenerate (no pairs) input returns 1.0."""
    resid = np.asarray(resid, float)
    ns = np.asarray(ns, float)
    nd = np.asarray(nd, float)
    ok = np.isfinite(resid) & (ns > _EPS) & (nd > _EPS)
    resid, ns, nd = resid[ok], ns[ok], nd[ok]
    if resid.size < 2:
        return 1.0
    dc = resid - np.median(resid)  # remove the systematic frame offset
    w = (ns * nd) / (ns + nd)
    den = float(np.sum(w))
    if den <= _EPS:
        return 1.0
    return float(max(np.sum(w * (dc * dc - (1.0 / ns + 1.0 / nd))) / den, 0.0))


def _adjacent_edges(chain: NodeChain):
    """Adjacent boundary↔region edges as (src, dst, src_is_boundary), each undirected edge once (i→right[i]).
    src presents its RIGHT face to dst; dst its LEFT face."""
    kind = np.asarray(chain.kind)
    right = np.asarray(chain.right)
    src = np.arange(kind.shape[0])
    dst = right
    m = dst >= 0
    src, dst = src[m], dst[m]
    return src, dst, kind[src] == BOUNDARY


def _adjacent_log_density_residuals(chain: NodeChain, geometry: NodeGeometry):
    """Oriented adjacent boundary↔region log gDNA-density residuals + the two source/dst gDNA counts, feeding
    the σ²_imp moment estimator. ``ρ = mass / eff_gdna`` oriented boundary→region so the systematic frame
    offset is one median-removable mode (the naive total-density frame, RNA included, evaluated pre-solve)."""
    ML = np.asarray(geometry.mass_left)
    MR = np.asarray(geometry.mass_right)
    EGL = np.asarray(geometry.eff_gdna_left)
    EGR = np.asarray(geometry.eff_gdna_right)
    src, dst, s_bnd = _adjacent_edges(chain)
    n_i, e_i, n_j, e_j = MR[src], EGR[src], ML[dst], EGL[dst]
    ok = (n_i > _EPS) & (n_j > _EPS) & (e_i > _EPS) & (e_j > _EPS)
    n_i, e_i, n_j, e_j, s_bnd = n_i[ok], e_i[ok], n_j[ok], e_j[ok], s_bnd[ok]
    lr_i = np.log(n_i / e_i)
    lr_j = np.log(n_j / e_j)
    resid = np.where(s_bnd, lr_i - lr_j, lr_j - lr_i)  # orient boundary→region (one mode)
    return resid, n_i, n_j


def adjacent_disagreement_variance(chain: NodeChain, geometry: NodeGeometry) -> float:
    """v1 total-density σ²_imp (``disagreement_shrinkage_prior_design_v2.md`` §2): the mean adjacent-node
    naive-gDNA log-density disagreement variance (Poisson sampling removed) — the population-average message
    process variance."""
    resid, n_i, n_j = _adjacent_log_density_residuals(chain, geometry)
    return _poisson_moment_var(resid, n_i, n_j)


def node_global_geometry(chain: NodeChain, geometry: NodeGeometry):
    """Per-node 'global' gDNA support ``(mass, eff)``: a REGION uses its contained mass over its contained
    gDNA eff-length; a BOUNDARY uses its both-side crossing mass over the SUMMED per-side density length
    ``E_l + E_r``. This is the basis the pass-0 gDNA-rate prior (`GdnaRatePrior`) is fit on and projected
    onto — shared by :func:`node_sweep` and ``calibrate`` so the fit and the projection use one definition.

    The boundary sum is ``E_l + E_r`` (not the old ``½(E_l+E_r)``) because ``eff_gdna_*`` is now the true
    per-face DENSITY length ``E[min(ℓ,R)]/2`` (`effective_length.boundary_side_eff_length`). The old ½ here
    was silently cancelling the ½ that was *missing* from the face length — which is why this frame read the
    correct ρ while every per-face MESSAGE read ρ/2. Both frames are now the same one."""
    is_reg = np.asarray(chain.kind) == REGION
    egl = np.asarray(geometry.eff_gdna_left, dtype=np.float64)
    egr = np.asarray(geometry.eff_gdna_right, dtype=np.float64)
    msl = np.asarray(geometry.mass_left, dtype=np.float64)
    msr = np.asarray(geometry.mass_right, dtype=np.float64)
    mass = np.where(is_reg, msl, msl + msr)
    eff = np.where(is_reg, egl, egl + egr)
    return mass, eff


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
    n_grid: int,
    logodds_window: float = 10.0,
    n_tilt: int | None = None,
    n_grid_ss: int | None = None,
    gdna_prior=None,
    fold_coarse_k: int = 33,
    fold_fine_k: int = 33,
    fold_sigma_coverage: float = 6.0,
    fold_refine_iters: int = 3,
    _capture: dict | None = None,
):
    """The belief-propagation sweep over the chain — gDNA AND per-strand RNA messages, in COUNT space (no
    ``(M/E)²`` Jacobian) — as a single **FORWARD-BACKWARD** pass.

    The chain is a forest of linear paths, so BP is exact in one forward + one backward pass (vs Gauss-Seidel /
    Jacobi which propagate one hop per pass). The message precision is the source's own HONEST belief precision
    ``pr = n_src/(n_src·vb_src + 1)`` = ``1/(Var(log f_c^src) + 1/n_src)`` — strand (composition) and count
    (sampling), and nothing else. The adjacent-pair overdispersion ``σ²_transfer`` is a PRIOR and is currently
    **zero**; see the ``_scan`` docstring for the derivation. There is no precision to refit and no outer
    fixed-point loop. The global prior is ANCHORED (every input fit once before the solve), so the single FB
    pass is exact.

    BEFORE the pass: the pass-0 NPMLE gDNA-rate prior (:class:`~.gdna_rate_prior.GdnaRatePrior`), fit once,
    belief-free, and passed as ``gdna_prior``. ``gdna_prior=None`` is a first-class PRIOR-FREE solve: ψ then
    carries the derived reference alone on both arms (``simplex_logodds._gdna_arm`` / ``_rna_arm``) — prior-free
    is not reference-free.

    The pass: (A) a batched message-free LOCAL solve → per-component (fraction, precision) — also the
    disagreement anchor ``lf*_loc`` / ``v*_loc``; (B) a FORWARD scan L→R accumulating the left-context message α
    and the forward belief (local ⊗ incoming message — NOT the reverse, true tree BP — so a thin node passes the
    upstream α on: the relay); (C) a BACKWARD scan R→L → β; (D) combine α⊗β and one batched FINAL solve.

    The per-node ψ solve integrates: the strand likelihood, the Jeffreys reference (single-strand nodes), the
    count-space global NB prior (ALL solvable nodes — the fork is dissolved), and the gDNA + per-strand RNA
    imputation messages (mature spliced RNA is carried by the boundary message terms; the dead node-local
    spliced floor was removed in the 2026-07 cleanup). The emission gate
    is a three-term Boolean over the components (gDNA / +RNA / −RNA): each RNA strand flows only where that strand
    is continuous on both endpoints (``free_s``); gDNA is genomically continuous and strand-agnostic, so it flows
    wherever there is unspliced facing mass — including across TSS/TES and other G1 seams (a locked all-gDNA node
    is a confident emitter). Only G2/G3 nodes with data are written; G1 sinks + empty nodes keep their init; a
    node with no neighbour message is solved on its own strand+global evidence (the intrinsic solve folds into the
    final batched solve at prec=0). Returns the resolved :class:`NodeBelief`."""
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    fp, fn = statics.free_pos, statics.free_neg
    mrp, mrn = statics.mrna_active_pos, statics.mrna_active_neg  # mature-crossing gate (see `_scan`)
    f_pos = np.asarray(belief.f_pos, dtype=np.float64).copy()
    f_neg = np.asarray(belief.f_neg, dtype=np.float64).copy()
    f_g = np.asarray(belief.f_g, dtype=np.float64).copy()
    # precision state (Phase 1: computed + carried; consumed by the honest message send in Phase 2).
    var_pos = np.asarray(belief.var_pos, dtype=np.float64).copy()
    var_neg = np.asarray(belief.var_neg, dtype=np.float64).copy()
    var_g = np.asarray(belief.var_gdna, dtype=np.float64).copy()

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

    def _local_solve(g_arr, gm=None, gp=None, rm=None, rp=None):
        """The per-node local/final solve (log-density log-odds backend). Returns the :class:`NodeDeconv`
        (the readout ``*_frac``/``*_frac_var`` + the free-coordinate seed ``lam_mean``/``lam_var``/
        ``theta_mean``/``theta_var``). Phase A calls it message-free; phase D passes the FB messages."""
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

    # The pass-0 gDNA-rate prior (`GdnaRatePrior`, fit ONCE before the sweep on ALL nodes' total unspliced
    # density) projected onto the f_g solve grid → the (n_nodes, K) additive term = log P(f_g·M/E) + the RNA
    # Jeffreys. ANCHORED (constant within the pass). It is EXTREMELY WEAK (n_eff≈0.15 pseudo-obs, never >1),
    # so it can never override a node's strand likelihood or the boundary messages — it supplies the
    # gDNA-vs-RNA population shape + the zero anchor, then lets strand/messages peel RNA out of the f_g=1
    # start. Replaces the retired seed/floor/global (`_gdna_seed_estimate`/`_floor_estimate`/`_global_logprior`)
    # and the density KDE (`_kde_logprior`) with this single count-space model (docs/calibration/
    # npmle_struggles.md §8-9). ``gdna_prior=None`` ⇒ no prior (strand + messages alone; a graceful degenerate).
    global_lp = (
        gdna_prior.logprior(solve_grid, mass_global, eff_global) if gdna_prior is not None else None
    )

    # Genomic node order for the forward/backward scans (within each ref path; left/right break at −1). The
    # scans are sequential per node, so iterate as a Python list of ints (faster than numpy scalar indexing).
    order_list = [int(x) for x in np.asarray(chain.order)]
    n_nodes = f_g.shape[0]

    # FORWARD-BACKWARD solve — ONE exact pass on the chain (a forest of linear paths) given the message-free
    # local beliefs + the ANCHORED global prior. (A) message-free LOCAL beliefs (one batched solve) →
    # per-component (fraction, precision); (B) FORWARD scan L→R accumulating the left-context message α + the
    # forward belief; (C) BACKWARD scan R→L → β; (D) combine α⊗β and one batched FINAL solve. Each node sees the
    # WHOLE-chain context (relay) in one pass — a thin node's forward belief is dominated by, and passes on, the
    # upstream α. ``true`` tree BP: the forward belief excludes β (no double-counting). The message precision is
    # the belief-free Poisson disagreement-variance (in `_scan`) so there is no var~mean fixed point ⇒ no outer loop.

    # (A) LOCAL message-free beliefs (backend-dispatched). The relay now carries the FREE-COORDINATE belief
    # (μ_λ, σ²_λ [, μ_θ, σ²_θ]) — the coherent (λ,θ) relay (docs/calibration/dof_pie_relay_derivation.md); the
    # fraction readout (`*_loc`) is retained for the locked reseed + the diagnostic capture.
    dc_loc = _local_solve(global_lp)
    fg_loc, fp_loc, fn_loc = dc_loc.gdna_frac, dc_loc.rna_pos_frac, dc_loc.rna_neg_frac
    vg_loc, vp_loc, vn_loc = dc_loc.gdna_frac_var, dc_loc.rna_pos_frac_var, dc_loc.rna_neg_frac_var
    lam_loc, lvar_loc = dc_loc.lam_mean.copy(), dc_loc.lam_var.copy()
    thm_loc, thv_loc = dc_loc.theta_mean.copy(), dc_loc.theta_var.copy()

    # G1-EMISSION FIX. A G1-locked seam (intergenic / TSS / TES / opposite-strand exon↔exon) is neither
    # single-strand nor AMBIG, so `_local_solve` SKIPS it and returns f_g=0. The sweep seeds each node's
    # running EMISSION belief (`fbg` in `_scan`) from that local solve, so a locked all-gDNA seam emitted
    # n_src = f_g·mass = 0 — SILENT, contradicting "a locked all-gDNA node is a confident emitter" (a
    # high-count gene-boundary gDNA crossing sent nothing to the adjacent exon). Reseed the locked nodes'
    # running belief from their LOCKED init so they emit their structural gDNA. Their local precision is
    # already 1/ε (vg_loc=0 for skipped nodes) so they stay locked and are not swayed by incoming messages;
    # their own FINAL value is unchanged (the solvable write-back keeps the init for non-solvable nodes) —
    # this only restores their OUTGOING messages. Empty nodes (mass 0) are reseeded but never emit (sm>0).
    #
    # Measured (post wire-fix, both suites): quick_3to1_5mb all 0.1167→0.0980 (−16%), worst 0.6703→0.5533;
    # ambig_dense_10mb all 0.1536→0.1519, worst 0.7002→0.6863. The gdna_none false-positive guard PASSES —
    # zero-gDNA error and over-call are UNCHANGED (0.0088 / 0.0543; 0.12 M / 1.29→1.27 M): un-muting does not
    # manufacture gDNA, because in a zero-gDNA library these seams carry ~no mass and the `sm > _EPS` gate
    # keeps them silent. This — not the σ²_imp precision cap — is why enriched exons never heard "gDNA here".
    locked = ~np.asarray(solvable, bool)
    fg_loc = np.where(locked, f_g, fg_loc)
    fp_loc = np.where(locked, f_pos, fp_loc)
    fn_loc = np.where(locked, f_neg, fn_loc)
    # locked nodes carry their structural {0,0,1} lock as a coordinate: μ_λ = logit(f_g)=+L (all-gDNA), σ²_λ=0
    # (unmovable). Their θ is irrelevant (no RNA emitted). This is the emission seed only; the FINAL write-back
    # keeps the init for non-solvable nodes.
    _Lw = float(logodds_window)
    lam_locked = np.clip(
        np.log(np.clip(f_g, _EPS, 1.0 - _EPS) / np.clip(1.0 - f_g, _EPS, 1.0 - _EPS)), -_Lw, _Lw
    )
    lam_loc = np.where(locked, lam_locked, lam_loc)
    lvar_loc = np.where(locked, 0.0, lvar_loc)

    def _scan(seq, nbr, sf, df):
        """Sequential scan — the coherent ``(λ,θ)`` relay (docs/calibration/dof_pie_relay_derivation.md;
        implementation plan §3-§6). Each node folds its ONE incoming message onto its running FREE-COORDINATE
        belief and emits per-component density messages derived from that belief. Returns the per-node message
        ``(mode, prec)`` per component — the return format is UNCHANGED, so ``_comb`` + the final solve consume
        it as before.

        **Why coordinates.** The running belief is stored as ``(μ_λ, σ²_λ)`` — ``λ = logit(f_g)``, the
        gDNA-vs-RNA-total log-odds — so the relayed pie ``f_g = σ(μ_λ)``, ``f_r = 1−f_g`` (split into ``f_pos``/
        ``f_neg`` by the tilt) is a COMPOSITION by construction: it sums to 1, every ``f_c ∈ [0,1]``, and
        ``n_c = f_c·M ≤ M``. The old three-independent-log-fraction relay violated all three (a boundary relayed
        ``fbp = 51.9``). The tilt ``θ`` is a NUISANCE and is held at the local seed in this pass (v1) — the
        per-strand messages still inform the AMBIG tilt at the FINAL solve; only the running-belief tilt is not
        relayed onward.

        **Message precision — the count-term (docs §5.2).** ``pr = 1/(Var(log f_c^src) + 1/M_src)`` written
        ``M_src/(M_src·Var + 1)`` so ``M_src=0 ⇒ pr=0``. ``M_src = sm`` is the source facing UNSPLICED mass; the
        spliced MEASUREMENT ``n_mat`` rides in the CONTENT but NOT the precision (item 2 owns the imputation
        precision; the mature-measurement channel is priority #3). ``σ²_transfer = 0`` (returns with the NPMLE).

        **The fold.** ``_fold_lambda`` — a two-stage EP moment-match of the running Gaussian against the gDNA
        (on ``log f_g``) and RNA-total (on ``log f_r``) factors: the two ``λ``-messages in tension on one axis.

        **The mature-crossing gate is now EXPLICIT** (``send_s = mrna_active[dst] or not mrna_active[src]``):
        the count-term no longer silences a gated edge via a zero sub-count, so the gate is spelled out."""
        mu_lam, var_lam = lam_loc.copy(), lvar_loc.copy()  # running λ belief (starts at the local seed)
        mu_th, var_th = thm_loc.copy(), thv_loc.copy()  # tilt θ: seeded, NOT relayed (v1 — a nuisance)
        amg, apg = np.zeros(n_nodes), np.zeros(n_nodes)  # gDNA message (mode, prec)
        amp, app = np.zeros(n_nodes), np.zeros(n_nodes)  # RNA-pos
        amn, apn = np.zeros(n_nodes), np.zeros(n_nodes)  # RNA-neg
        EGs, EGd, ERs, ERd = EG[sf], EG[df], ER[sf], ER[df]
        MSs, MSd, SPs, SNs = MS[sf], MS[df], SP[sf], SN[sf]
        ESPs = ESP[sf]  # source-face spliced eff-len (for the mature-RNA MEASUREMENT message)
        SPd, SNd, ESPd = SP[df], SN[df], ESP[df]  # dest-face spliced — the mature ABSORBED at a junction
        for i in seq:
            lsrc = nbr[i]
            if lsrc < 0:
                continue
            md = MSd[i] if MSd[i] > _EPS else _EPS
            egd = EGd[i] if EGd[i] > _EPS else _EPS
            erd = ERd[i] if ERd[i] > _EPS else _EPS
            sm = MSs[lsrc]  # source facing UNSPLICED mass = the count-term M_src
            # source COHERENT fractions + log-fraction variances from its (λ,θ) belief (delta-method).
            ls, vls = mu_lam[lsrc], var_lam[lsrc]
            sin_t = math.sin(mu_th[lsrc])
            cos_t = math.cos(mu_th[lsrc])
            fg_s = 1.0 / (1.0 + math.exp(-ls))
            fr_s = 1.0 - fg_s
            fp_s = fr_s * (1.0 + sin_t) / 2.0
            fn_s = fr_s * (1.0 - sin_t) / 2.0
            v_logfg = fr_s * fr_s * vls  # Var(log f_g)  = (1−f_g)²·σ²_λ
            v_logfr = fg_s * fg_s * vls  # Var(log f_r)  =  f_g²·σ²_λ   (RNA-total)
            vts = var_th[lsrc]
            v_logfp = v_logfr + (cos_t / max(1.0 + sin_t, _EPS)) ** 2 * vts  # +θ term (0 for single-strand)
            v_logfn = v_logfr + (cos_t / max(1.0 - sin_t, _EPS)) ** 2 * vts
            # STRUCTURAL emission gates + the EXPLICIT mature-crossing gate.
            emit_g = sm > _EPS
            send_p = mrp[i] or not mrp[lsrc]
            send_n = mrn[i] or not mrn[lsrc]
            emit_p = fp[lsrc] and fp[i] and (sm > _EPS or SPs[lsrc] > _EPS) and send_p
            emit_n = fn[lsrc] and fn[i] and (sm > _EPS or SNs[lsrc] > _EPS) and send_n
            lam_factors = []
            # ---- gDNA density message (a factor on log f_g) ----
            if emit_g:
                eg = EGs[lsrc] if EGs[lsrc] > _EPS else _EPS
                rho = fg_s * sm / eg  # source gDNA DENSITY (f_g^src ≤ 1 ⇒ n_src ≤ sm by construction)
                mo = math.log(max(rho, 1.0 / egd) / (md / egd))  # → dst log-f_g frame
                pr = sm / (sm * v_logfg + 1.0)  # COUNT-TERM: 1/(Var(log f_g) + 1/M_src)
                amg[i], apg[i] = mo, pr
                if pr > 0.0:
                    lam_factors.append((True, mo, pr))
            # ---- RNA per-strand densities: emission (per strand) + the RNA-TOTAL λ-factor ----
            # The imputed density is the NASCENT unspliced RNA (`f_s·sm/E_r`, gated) + the junction spliced
            # MEASUREMENT (`n_mat/E_spl`, B→exon; rides in CONTENT, NOT precision — #3) − the dst absorption.
            rho_r = 0.0
            if emit_p:
                er = ERs[lsrc] if ERs[lsrc] > _EPS else _EPS
                esp = ESPs[lsrc] if ESPs[lsrc] > _EPS else _EPS
                rho_pos = (
                    fp_s * sm / er
                    + SPs[lsrc] / esp
                    - SPd[i] / (ESPd[i] if ESPd[i] > _EPS else _EPS)
                )
                mo = math.log(max(rho_pos, 1.0 / erd) / (md / erd))  # → dst log-f_pos frame
                pr = sm / (sm * v_logfp + 1.0)  # IMPUTATION precision (unspliced M_src; n_mat excluded)
                amp[i], app[i] = mo, pr
                if rho_pos > 0.0:
                    rho_r += rho_pos
            if emit_n:
                er = ERs[lsrc] if ERs[lsrc] > _EPS else _EPS
                esp = ESPs[lsrc] if ESPs[lsrc] > _EPS else _EPS
                rho_neg = (
                    fn_s * sm / er
                    + SNs[lsrc] / esp
                    - SNd[i] / (ESPd[i] if ESPd[i] > _EPS else _EPS)
                )
                mo = math.log(max(rho_neg, 1.0 / erd) / (md / erd))  # → dst log-f_neg frame
                pr = sm / (sm * v_logfn + 1.0)
                amn[i], apn[i] = mo, pr
                if rho_neg > 0.0:
                    rho_r += rho_neg
            # RNA-TOTAL density message (a factor on log f_r) — the second λ-factor; gDNA + RNA-total in
            # tension on the SAME axis λ (this is what makes the pie coherent — docs §2.2).
            if rho_r > _EPS:
                pr_r = sm / (sm * v_logfr + 1.0)
                if pr_r > 0.0:
                    lam_factors.append((False, math.log(max(rho_r, 1.0 / erd) / (md / erd)), pr_r))
            # ---- FOLD the λ-messages onto the dst's running belief (two-stage EP moment-match) ----
            if lam_factors:
                mu_lam[i], var_lam[i] = _fold_lambda(
                    mu_lam[i],
                    var_lam[i],
                    lam_factors,
                    L=_Lw,
                    coarse_k=fold_coarse_k,
                    fine_k=fold_fine_k,
                    sigma_cov=fold_sigma_coverage,
                    refine=fold_refine_iters,
                )
        if _capture is not None:
            # DIAGNOSTIC — the coherent relay pie: f_g=σ(μ_λ), f_pos/f_neg from (f_r, θ). Sums to 1, each ≤1
            # BY CONSTRUCTION (the S2 invariant). Forward scan appends [0], backward [1].
            _fg = 1.0 / (1.0 + np.exp(-mu_lam))
            _fr = 1.0 - _fg
            _st = np.sin(mu_th)
            _capture.setdefault("_relay_pie", []).append(
                (_fg, _fr * (1.0 + _st) / 2.0, _fr * (1.0 - _st) / 2.0)
            )
        return amg, apg, amp, app, amn, apn

    a = _scan(order_list, left, 1, 0)  # forward (α: left context)
    b = _scan(order_list[::-1], right, 0, 1)  # backward (β: right context)

    # (D) combine α⊗β (precision-weighted product) per component → one batched FINAL solve.
    def _comb(am_a, ap_a, am_b, ap_b):
        pc = ap_a + ap_b
        return np.where(pc > _EPS, (ap_a * am_a + ap_b * am_b) / np.maximum(pc, _EPS), 0.0), pc

    mode_g, prec_g = _comb(a[0], a[1], b[0], b[1])
    mode_p, prec_p = _comb(a[2], a[3], b[2], b[3])
    mode_n, prec_n = _comb(a[4], a[5], b[4], b[5])
    # (D) FINAL solve with the FB messages (backend-dispatched). The final solve is UNCHANGED — it consumes the
    # per-component messages `_comb` produced (now COHERENT); the coordinate relay only changed how they are built.
    dc_fin = _local_solve(global_lp, mode_g, prec_g, (mode_p, mode_n), (prec_p, prec_n))
    mg_, mp_, mn_ = dc_fin.gdna_frac, dc_fin.rna_pos_frac, dc_fin.rna_neg_frac
    vg_, vp_, vn_ = dc_fin.gdna_frac_var, dc_fin.rna_pos_frac_var, dc_fin.rna_neg_frac_var
    # write back only SOLVABLE nodes (G1 sinks / empty keep their signature-binary init).
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
        _capture.update(
            fg_loc=fg_loc,
            fg_strand=fg_strand,
            fp_loc=fp_loc,
            fn_loc=fn_loc,
            vg_loc=vg_loc,
            vp_loc=vp_loc,
            vn_loc=vn_loc,
            a_fwd=a,
            b_bwd=b,
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
