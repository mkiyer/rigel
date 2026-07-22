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

The gDNA population prior is the pass-0 count-space **`DensityNPMLE`** (`npmle`), fit ONCE before
the sweep on all nodes' total unspliced density and projected onto each node's ψ (`DensityNPMLE.logprior`) —
extremely weak, so strand + messages dominate. It replaces the retired seed/floor/global + the density KDE.

Module layout. The per-node geometry / belief / statics / init primitives (`build_node_geometry`,
`build_node_statics`, `init_beliefs`, `NodeGeometry`/`NodeBelief`/`NodeStatics`) live in the
lower `node_geometry` module and are re-exported here for the calibrator's convenience; this module owns:
* `node_global_geometry` — the per-node ``(mass, eff)`` support the rate prior is fit on and projected onto.
* `node_sweep` — the single forward-backward sweep. Message precision is the source's own HONEST belief
  precision ``pr = n_src/(n_src·vb_src + 1)`` = ``1/(Var(log f_c^src) + 1/n_src)`` — strand (composition) and
  count (sampling), per component, from the RUNNING belief — plus the belief-free ``σ²_transfer`` (the NPMLE
  enrichment-crossing projection ``F1``, `DensityNPMLE.project`), added to the message variance. The retired
  total-density disagreement estimator is gone from production (relocated to `scripts/debug/`).
* `chain_region_deconv` / `chain_boundary_side_deconv` — project the converged belief back to the per-region
  / per-boundary-side masses the `CalibrationResult` consumes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .node_chain import REGION, NodeChain
from .signature import coarse_type_array
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

# THE SPLICED ABSORPTION (production): the boundary subtracts ITS OWN spliced-RNA density (``SPd[i]/ESPd[i]``,
# the dst-face spliced) from the incoming message's RNA density and adds the DST exon's spliced (``SPs/ESPs``). A
# splice junction directly measures a pure-RNA spliced population; that mass is accounted for BEFORE the residual
# unspliced RNA is imputed onto the boundary's crossing. The two densities are disjoint independent measurements,
# so the residual can go negative — which just means the boundary's spliced already absorbed all the incoming RNA
# (and more) ⇒ ~zero unspliced RNA left (handled by the honest clamp below).
#
# THE HONEST CLAMP (message_absorption_fix.md Step 0): a residual unspliced-RNA density that the spliced
# absorption drives at/below the count floor ``1/erd`` is a genuine but IMPRECISE zero — a short boundary cannot
# measure a confident zero — so it is sent with the ONE-COUNT precision of that floor (``n_eff = rho_floor·erd =
# 1``), NOT the source's full count ``sm`` (see the ``n_eff`` sites in ``_scan``). This stops the saturated
# subtraction from being laundered into a CONFIDENT "no RNA" (the confident-FP seed); the honest weak zero lets
# the Phase-2 DNA hyperprior pin the node later.
#
# THE τ-PRECISION CORE (docs/calibration/message_precision_derivation.md §3). A message's COMPOSITION precision
# is sourced from a REFERENCE-FREE evidence quantity ``τ`` (I_strand + I_struct + relayed pr) instead of the
# running BELIEF variance ``σ²_λ`` — which pools the shared Beta(½,½) reference and manufactures confidence on
# composition-vacuous unstranded chains (the 35× phantom cascade). A vacuous source (unstranded, no structural
# lock) has ``τ = 0 ⇒ pr = 0``; a strand anchor / structural lock propagates honestly. The strand seed I_strand
# carries the DERIVED overdispersion noise floor (see ``_scan``) so a κ≈½ sampling whisper cannot seed phantom
# precision. Measured: confidently-wrong mass (the NPMLE-corrupting quantity) 10.4% → 1.4% vs the pooled-belief
# precision. I_struct is composition-certain ONLY for true intergenic REGION locks, NOT G1 boundary SEAMS
# (TSS/TES, opposite-strand exon↔exon) — a seam is locked to gDNA by structure but sits between RNA-carrying
# exons, so its crossing mass is RNA-contaminated and making it certain turns it into a phantom-gDNA emitter that
# compounds along the chain (measured: τ→12 in gdna_none); a true intergenic region carries ~0 mass in a
# zero-gDNA library, so it is safe.
#
# THE CLIFF-CROSSING MESSAGE MODE (docs/calibration/cliff_message_derivation.md). Two regimes, gated per edge by
# whether either endpoint is an EXON region (``is_exon_node`` / ``use_shift`` in ``node_sweep``):
#   * CLEAN edges (intron/intergenic ↔ boundary, no mature): the eff-length-frame LOG-ODDS SHIFT — impute the
#     dst ``f_c`` as the source COMPOSITION (per-component imputed mass ``M_c = ρ_c^src·E_c^dst`` normalized by
#     ``ΣM``, ≡ ``λ_dst = λ_src + log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)``). The capture enrichment cancels
#     identically (CLIFF-INVARIANT); the shift is nonzero because gDNA and RNA have different FL distributions.
#     The intron factory makes the source accurate here, so the confident shift is safe.
#   * EXON ↔ boundary edges: the DENSITY mode ``log(ρ_c^src·E_c^dst / md)`` (÷ the dst's OBSERVED total), which
#     keeps the gDNA f_g DECOUPLED from the error-prone mature removal + anchored to real data — robust where
#     the (often unstranded) exon source is NOT reliable. The ``±`` spliced absorption rides in ``rho_pos/neg``.
# §9 records the REJECTED unifiers (the composition shift on all edges, and the gDNA-only / all-component +logR
# "hybrid enrichment-corrected density"): all regress on unstranded exons — that is an identifiability floor, not
# a mode defect. PRECISION is the lever, not the mode. Do not re-attempt a composition-conserving mode on exons.



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


def node_global_geometry(chain: NodeChain, geometry: NodeGeometry):
    """Per-node 'global' gDNA support ``(mass, eff)``: a REGION uses its contained mass over its contained
    gDNA eff-length; a BOUNDARY uses its both-side crossing mass over the SUMMED per-side density length
    ``E_l + E_r``. This is the basis the enrichment NPMLE (`DensityNPMLE`) is fit on and projected
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


def _boundary_spliced_mass_increment(chain, geometry, eff):
    """The mature-INCLUSIVE mass increment for boundaries (0 on regions): fold the boundary's spliced (mature)
    density ``ρ_spl = Σ_side spl_side/E_spl_side`` (per-frame, NOT raw mass — the composition-conditional trap)
    into an equivalent unspliced-frame mass ``ρ_spl·eff``. Used ONLY on the exon↔boundary σ²_transfer edge, where
    the exon flank's contained density already carries within-exon mature (`bp_solver` directional σ²_T)."""
    is_bnd = np.asarray(chain.kind) != REGION
    spl_l = np.asarray(geometry.spliced_pos_left, np.float64) + np.asarray(geometry.spliced_neg_left, np.float64)
    spl_r = np.asarray(geometry.spliced_pos_right, np.float64) + np.asarray(geometry.spliced_neg_right, np.float64)
    rho_spl = spl_l / np.maximum(geometry.eff_spl_left, _EPS) + spl_r / np.maximum(geometry.eff_spl_right, _EPS)
    return np.where(is_bnd, rho_spl * eff, 0.0)


@dataclass(frozen=True)
class StrandEvidence:
    """The per-node reference-free **composition evidence** — the identifiability model in one place
    (`docs/calibration/message_system_derivation.md` §6B). It seeds the message precision (the τ cavity) AND is
    the natural home for the pass-0 solvability gate (§6B: a node skips its solve iff a free axis has zero total
    precision — Phase C). Two evidence channels seed it here; the running τ accumulates relayed message precision
    on top during the sweep:

    * ``tau0_lam`` / ``tau0_th`` — **I_strand**: the differential-κ strand Fisher info on the λ (gDNA-vs-RNA) and
      θ (tilt) axes. Zero on unstranded data (κ=½) by the DERIVED noise-floor deadband
      ``disc = 4·max(0, (κ−½)² − σ²_d)`` — this, not any gate, is what kills the phantom.
    * ``struct_lock`` — **I_struct**: signature composition-certainty (intergenic REGION nodes today; a lock ⇒
      ``ev=0`` ⇒ high message precision governed only by the honest count/transfer terms).
    """

    tau0_lam: np.ndarray  # I_strand on the λ axis (float64[m])
    tau0_th: np.ndarray  # I_strand on the θ axis (float64[m])
    struct_lock: np.ndarray  # I_struct — composition-certain nodes (bool[m])


def _compile_strand_evidence(
    u_pos,
    u_neg,
    fg_loc,
    *,
    kappa: float,
    od_g: float,
    od_r: float,
    n_gdna_obs: float,
    n_rna_obs: float,
    is_region,
    locked,
) -> StrandEvidence:
    """Compile the reference-free strand/structure evidence seeds (`StrandEvidence`), evaluated at the
    message-free local ``fg_loc``. Pure: no cross-node coupling, unit-testable.

    ``I_strand(λ) = N·(2κ−1)²·[f_g(1−f_g)]² / (4 p(1−p))``, ``p = κ + f_g(½−κ)`` — the strand Fisher info,
    IDENTICALLY 0 at κ=½ (unstranded). The count enters as the OVERDISPERSED effective count
    ``N_eff = N/(1+(N−1)ω)`` (the honest power saturates at ~1/ω, not the raw depth), and the discriminability
    ``disc`` carries the DERIVED noise floor ``σ²_d = ¼·(1/N_rna+ω_r) + ¼·(1/N_gdna+ω_g)`` — a κ within √σ²_d of
    ½ is not composition signal (the deadband that kills the phantom; `message_precision_derivation.md`). The
    ``1/N_gdna`` term gates a gDNA-free library (N_gdna=0 ⇒ σ²_d→∞ ⇒ disc=0). ``I_struct`` is composition-certain
    ONLY for true intergenic REGION nodes, never G1 boundary SEAMS (whose crossing mass is RNA-contaminated)."""
    _n_raw = np.asarray(u_pos, dtype=np.float64) + np.asarray(u_neg, dtype=np.float64)
    _n_str = _n_raw / (1.0 + np.maximum(_n_raw - 1.0, 0.0) * od_r)
    _fgl = np.clip(np.asarray(fg_loc, dtype=np.float64), _EPS, 1.0 - _EPS)
    _pmix = np.clip(kappa + _fgl * (0.5 - kappa), _EPS, 1.0 - _EPS)
    _sig2_d = 0.25 * (1.0 / max(float(n_rna_obs), _EPS) + od_r) + 0.25 * (
        1.0 / max(float(n_gdna_obs), _EPS) + od_g
    )
    _disc = 4.0 * max(0.0, (kappa - 0.5) ** 2 - _sig2_d)
    _i_strand = _n_str * _disc * (_fgl * (1.0 - _fgl)) ** 2 / (4.0 * _pmix * (1.0 - _pmix))
    struct_lock = np.asarray(locked, dtype=bool) & np.asarray(is_region, dtype=bool)
    return StrandEvidence(
        tau0_lam=_i_strand.copy(), tau0_th=_i_strand.copy(), struct_lock=struct_lock
    )


# ── The two message MODES (docs/calibration/message_propagation_arithmetic.md) ─────────────────────────────
# A message conveys, per component, a target **log-fraction** ``log f_c`` in the dst frame. There are exactly two
# ways to set it, and (Stage-1 verified, ``mode_verify.py``) the solver computes each bit-for-bit:
#   • the composition SHIFT (§3/§4a): normalize the imputed per-component masses by their IMPUTED total ⇒ the
#     capture scale cancels (cliff-invariant). Correct on every gate-EQUAL edge.
#   • the DENSITY mode (§8): ÷ the dst's OBSERVED total ``md`` ⇒ the enrichment is baked in (fails across a cliff).
#     Retained on exon / unequal-gate edges pending the gate-equality flip (message_mode_implementation_plan.md).
# Both carry the SAME one-fragment resolution floor — a node can never be more certain than its one-fragment
# opportunity (``comp_fl = 1/M_dst`` for the shift; ``max(ρ_c, 1/E_c)`` for the density). NOT an arbitrary epsilon;
# it is also the domain guard that keeps ``log`` finite when the mature subtraction drives ``ρ_c ≤ 0``.


def _mode_shift(mass_c: float, den: float, comp_fl: float) -> float:
    """Composition-SHIFT mode ``log(max(M_c/ΣM, comp_fl))`` — cliff-invariant (§4a)."""
    return math.log(max(mass_c / den, comp_fl))


def _mode_density(rho_c: float, eff_c: float, md: float) -> float:
    """DENSITY mode ``log(max(ρ_c, 1/E_c)·E_c / md)`` — observed-total anchor (§8)."""
    return math.log(max(rho_c, 1.0 / eff_c) / (md / eff_c))


def _pred_precision(count: float, v_log: float, s2t: float) -> float:
    """The deconvolution-PREDICTION precision ``1 / (Var(log f) + 1/count + σ²_transfer)``
    (emission_and_precision_derivation.md §1). It is **zero** when the composition is *unseen*
    (``v_log = ∞`` — a no-evidence source) or there is no count — so a composition-vacuous source emits
    ``pr→0`` and is ignored, with no emission gate and no nan (``count·∞`` is never formed). The spliced
    MEASUREMENT precision is a separate, independent channel and does NOT pass through here (§4)."""
    if count > 0.0 and math.isfinite(v_log):
        return count / (count * (v_log + s2t) + 1.0)
    return 0.0


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
    enrichment_prior=None,
    intron_prior=None,
    fold_coarse_k: int = 33,
    fold_fine_k: int = 33,
    fold_sigma_coverage: float = 6.0,
    fold_refine_iters: int = 3,
    transfer_variance: bool = True,
    _capture: dict | None = None,
):
    """The belief-propagation sweep over the chain — gDNA AND per-strand RNA messages, in COUNT space (no
    ``(M/E)²`` Jacobian) — as a single **FORWARD-BACKWARD** pass.

    The chain is a forest of linear paths, so BP is exact in one forward + one backward pass (vs Gauss-Seidel /
    Jacobi which propagate one hop per pass). The message precision is the source's own HONEST belief precision
    plus the belief-free transfer variance: ``pr = 1/(Var(log f_c^src) + 1/n_src + σ²_transfer)`` — strand
    (composition), count (sampling), and the enrichment-crossing damping ``σ²_transfer`` from the NPMLE
    projection (fit once, belief-free — ``docs/calibration/npmle_projection_variance_design.md``; off when
    ``transfer_variance=False`` or there is no prior). There is no precision to refit and no outer
    fixed-point loop. The global prior is ANCHORED (every input fit once before the solve), so the single FB
    pass is exact.

    BEFORE the pass: the pass-0 NPMLE gDNA hyperprior (:class:`~.npmle.DensityNPMLE`), fit once,
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
    # integer unspliced flux per face — the Poisson n for message PRECISION (mass is the density numerator
    # only; it is fractional AND split across nodes by the accumulator, so 1/mass is not a counting variance).
    MSN = (geometry.n_unspl_left, geometry.n_unspl_right)
    SP = (geometry.spliced_pos_left, geometry.spliced_pos_right)
    SN = (geometry.spliced_neg_left, geometry.spliced_neg_right)
    # integer spliced COUNT (flux) — MASS is the density numerator, COUNT is what a Poisson VARIANCE needs
    # (`Var(log ρ_m)=1/n`, not `1/mass`; node_geometry §spliced_n_*). Consumed by the mature-dilution variance.
    SPN = (geometry.spliced_n_pos_left, geometry.spliced_n_pos_right)
    SNN = (geometry.spliced_n_neg_left, geometry.spliced_n_neg_right)
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

    # The belief-free PROJECTION message transfer variance (transfer_variance_formal_derivation.md): each node's
    # total density projected onto the NPMLE ENRICHMENT landscape → (mu_proj, var_proj). The per-edge
    # σ²_transfer = var_proj[dst] + (mu_proj[dst] − mu_proj[src])² damps a message across a capture-enrichment
    # crossing (mode gap²) and floors at h². This is Role A — message PRECISION, not a composition claim — so it
    # runs on the ENRICHMENT NPMLE (``enrichment_prior``), INDEPENDENTLY of the composition arm above. Fit once,
    # belief-free, ANCHORED. Falls back to ``gdna_prior`` for callers that pass one prior (backward compat).
    # σ²_transfer = 0 with no prior or when disabled (``transfer_variance=False``, for ablation).
    proj_prior = enrichment_prior if enrichment_prior is not None else gdna_prior
    if transfer_variance and proj_prior is not None:
        mu_proj, var_proj = proj_prior.project(mass_global, eff_global)
    else:
        mu_proj = np.zeros_like(mass_global)
        var_proj = np.zeros_like(mass_global)
    # DIRECTIONAL spliced-density σ²_transfer: a mature-INCLUSIVE boundary projection, selected per-edge (in
    # ``_scan``) ONLY on the exon↔boundary edge — where the exon endpoint's contained density already carries
    # within-exon mature (unspliced, no junction). Without it the exon-vs-boundary mode gap is inflated by the
    # mature-inclusion asymmetry (the boundary crossing is mature-free) and σ²_transfer GAGS the reliable dense
    # exon edge. The intron↔boundary edge stays mature-FREE (the intron has no mature) — hence per-EDGE, not
    # per-node. Fit-basis note: the enrichment NPMLE is fit on the bare density; this mature-inclusive projection
    # only re-reads the exon-facing boundary onto that landscape (a higher density → a mode nearer the exon's).
    if transfer_variance and proj_prior is not None:
        _mass_mat = mass_global + _boundary_spliced_mass_increment(chain, geometry, eff_global)
        mu_proj_mat, var_proj_mat = proj_prior.project(_mass_mat, eff_global)
    else:
        mu_proj_mat, var_proj_mat = mu_proj, var_proj

    # Genomic node order for the forward/backward scans (within each ref path; left/right break at −1). The
    # scans are sequential per node, so iterate as a Python list of ints (faster than numpy scalar indexing).
    order_list = [int(x) for x in np.asarray(chain.order)]
    n_nodes = f_g.shape[0]
    # Per-node EXON-region flag (coarse_type == 2). The cliff-crossing log-odds shift
    # (`cliff_message_derivation.md`) applies on CLEAN transitions only — a message whose region endpoint is an
    # intron/intergenic (no mature), where the intron factory makes the source accurate. An EXON endpoint keeps
    # the DENSITY mode (observed-md anchor), which is the correct permanent design: a composition-conserving
    # mode on exon edges over-trusts the (often unstranded) exon source and regresses (§9 — the exon floor is
    # identifiability, not the mode). ``is_exon_node`` gates the shift per edge.
    _rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    _ri = np.clip(np.asarray(chain.ref_idx, dtype=np.int64), 0, _rtype.shape[0] - 1)
    is_exon_node = ((np.asarray(chain.kind) == REGION) & (_rtype[_ri] == 2)).tolist()
    is_intron_node = ((np.asarray(chain.kind) == REGION) & (_rtype[_ri] == 1)).tolist()  # INTRON == 1
    _is_bnd = (np.asarray(chain.kind) != REGION).tolist()  # boundary-dst edges: geo-mean crossing + density frame

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
    thm_loc = dc_loc.theta_mean.copy()  # tilt θ mode (seeded, not relayed); its variance is not used in v1

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

    # The §14 GEO-MEAN crossing decode (`rho_g_cross = sqrt(rho_g^exon * rho_g^intron)`) was RETIRED here
    # (roadmap R3). It was an unweighted, pre-scan stand-in for the MISSING exon->boundary imputation and it
    # bypassed the precision machinery entirely. Both exon<->boundary directions now carry real,
    # precision-bearing messages (shift +/- c_b), and the alpha-beta integrator -- already a
    # precision-weighted GEOMETRIC mean -- reconciles a boundary's two flanks properly. Removal measured
    # INERT: byte-identical on the gdna_none phantom guard (9 scenarios), the grounded full-transcript toy,
    # and the calibration/native/golden suites, because its only remaining destinations were structurally
    # pinned seams.
    _is_reg_arr = np.asarray(chain.kind) == REGION

    # --- MATURE-DILUTION per boundary (docs/calibration/exon_boundary_mature_dilution_plan.md) ---
    # A mature fragment is SPLICED, so it cannot cross an intron-exon junction contiguously: the boundary's
    # unspliced crossing is mature-FREE (ρ_bg+ν) while the exon carries mature too (ρ_bg+ν+μ). The boundary's
    # TOTAL (unspliced + spliced) therefore has the SAME component set as the exon, so the capture enrichment
    # e(x) CANCELS and the composition transfer across the cliff is exactly
    #     f_g^boundary = f_g^exon · (D_B + S_B)/D_B     ⇒  an ADDITIVE  c_b = log1p(S_B/D_B)  on the log-f mode
    #     (+c_b exon→boundary: remove mature;  −c_b boundary→exon: restore it).
    # ZERO free constants; the correction uses only the boundary's OWN measured spliced/unspliced split, so it
    # never extrapolates within-exon mature via eff-length ratios (the step that failed before, §9).
    # SUMMED-EFF frame — validated on the toy (~2% @20k frags); the per-side density sum over-corrects ~2×
    # because spliced mass lands on ONE face only.
    _D_B = (MS[0] + MS[1]) / np.maximum(EG[0] + EG[1], _EPS)  # unspliced crossing density
    _spl_mass = SP[0] + SN[0] + SP[1] + SN[1]  # spliced (mature) MASS
    _S_B = _spl_mass / np.maximum(ESP[0] + ESP[1], _EPS)  # spliced (mature) density
    _r_mat = np.where(_D_B > _EPS, _S_B / np.maximum(_D_B, _EPS), 0.0)
    mature_dilution = np.where(~_is_reg_arr, np.log1p(_r_mat), 0.0)
    # Counting-noise variance of c_b, propagated through c=log(1+r):  Var(c) = [r/(1+r)]²·Var(log r), and
    # Var(log r) = 1/n_s + 1/n_d for two independent Poisson channels. The variance MUST use the INTEGER
    # COUNTS, not the fractional mass — mass sums per-fragment shares, so 1/mass mis-states the Poisson
    # variance (node_geometry §spliced_n_*: `Var(log ρ_m)=1/n`, Kish n_eff ≥ mass). Gamma posterior (n+1)
    # keeps it finite at n=0. NOTE: this models SAMPLING noise only — a probe-ATTENUATED junction is a
    # systematic BIAS in the mode, which no variance term can represent (see the plan doc §5 risk 1).
    _spl_n = SPN[0] + SNN[0] + SPN[1] + SNN[1]  # integer spliced flux
    _unspl_n = np.asarray(statics.u_pos, np.float64) + np.asarray(statics.u_neg, np.float64)
    _var_mat = np.where(
        ~_is_reg_arr,
        (_r_mat / (1.0 + _r_mat)) ** 2 * (1.0 / (_spl_n + 1.0) + 1.0 / (_unspl_n + 1.0)),
        0.0,
    )

    # ---- the reference-free composition evidence (StrandEvidence) ----
    # I_strand (the differential-κ deadband Fisher info — 0 on unstranded, the real phantom fix) + I_struct
    # (signature composition-certainty), evaluated at the message-free local f_g. The running τ (below)
    # accumulates relayed message precision on top during the sweep. Derivation + the solvability model this
    # object is the home for: `docs/calibration/message_system_derivation.md` §6B.
    _ev = _compile_strand_evidence(
        statics.u_pos,
        statics.u_neg,
        fg_loc,
        kappa=kappa,
        od_g=od_g,
        od_r=od_r,
        n_gdna_obs=n_gdna_obs,
        n_rna_obs=n_rna_obs,
        is_region=(np.asarray(chain.kind) == REGION),
        locked=locked,
    )
    tau0_lam, tau0_th, struct_lock = _ev.tau0_lam, _ev.tau0_th, _ev.struct_lock

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
        precision; the mature-measurement channel is priority #3). ``σ²_transfer`` (``s2t``) is the belief-free
        NPMLE-projection enrichment-crossing damping ``var_proj[dst] + (mu_proj[dst] − mu_proj[src])²``, shared
        by every component at pass-0 (capture is nucleic-acid-agnostic).

        **The fold.** ``_fold_lambda`` — a two-stage EP moment-match of the running Gaussian against the gDNA
        (on ``log f_g``) and RNA-total (on ``log f_r``) factors: the two ``λ``-messages in tension on one axis.

        **Emission is gated by the STRUCTURAL per-strand ``free_s`` continuity only** (the mature-crossing gate
        was dismantled, 2026-07-16): each RNA strand flows wherever it is continuous on both endpoints. Mature RNA
        at a junction is reconciled in DENSITY space by the spliced source/sink absorption (``±SPs``/``−absorb``
        in ``rho_pos``/``rho_neg``), not by a separate gate. A full node-local nascent/mature split
        (``ρ_nascent = ρ_RNA − ρ_mature``) is not modelled here — the per-locus EM separates nascent from mature
        downstream."""
        mu_lam, var_lam = lam_loc.copy(), lvar_loc.copy()  # running λ belief (starts at the local seed)
        mu_th = thm_loc.copy()  # tilt θ: seeded, NOT relayed (v1 — a nuisance)
        # running reference-free evidence precision τ; accumulates the relayed pr (the cavity along this
        # sweep's direction).
        tau_lam = tau0_lam.copy()
        tau_th = tau0_th  # seed-only (θ not relayed) — read, not mutated
        amg, apg = np.zeros(n_nodes), np.zeros(n_nodes)  # gDNA message (mode, prec)
        amp, app = np.zeros(n_nodes), np.zeros(n_nodes)  # RNA-pos
        amn, apn = np.zeros(n_nodes), np.zeros(n_nodes)  # RNA-neg
        # DIAGNOSTIC (inert unless _capture): per-edge precision-term decomposition for the honest-precision
        # audit — the source count `sm`, the source composition variances `v_logf*`/`vls`, the transfer var
        # `s2t`, the source `fg_s`. Reveals WHICH term lets a message run to high precision.
        _pt = {k: np.zeros(n_nodes) for k in ("sm", "vlfg", "vlfp", "s2t", "fgs", "vls")}
        EGs, EGd, ERs, ERd = EG[sf], EG[df], ER[sf], ER[df]
        MSs, MSd, SPs, SNs = MS[sf], MS[df], SP[sf], SN[sf]
        MSNs = MSN[sf]  # source-face integer unspliced flux — the Poisson n for the message precision
        ESPs = ESP[sf]  # source-face spliced eff-len (for the mature-RNA MEASUREMENT message)
        SPd, SNd, ESPd = SP[df], SN[df], ESP[df]  # dest-face spliced — the mature ABSORBED at a junction
        for i in seq:
            lsrc = nbr[i]
            if lsrc < 0:
                continue
            md = MSd[i] if MSd[i] > _EPS else _EPS
            egd = EGd[i] if EGd[i] > _EPS else _EPS
            erd = ERd[i] if ERd[i] > _EPS else _EPS
            sm = MSs[lsrc]  # source facing UNSPLICED MASS — the DENSITY numerator (ρ = mass/eff)
            smn = MSNs[lsrc]  # source facing integer unspliced COUNT — the Poisson n for the PRECISION
            # belief-free message transfer variance for THIS edge (src=lsrc → dst=i): the enrichment-crossing
            # damping (F1). Shared by every component (capture is nucleic-acid-agnostic ⇒ the crossing is a
            # discontinuity in gDNA AND RNA at pass-0). 0 when the projection is off.
            # DIRECTIONAL spliced-density σ²_transfer: use the mature-INCLUSIVE boundary projection on an
            # exon↔boundary edge only (a boundary endpoint whose OTHER endpoint is an exon region);
            # intron↔boundary stays mature-free. mu_proj_mat == mu_proj when there is no projection ⇒ no-op there.
            _dst_mat = _is_bnd[i] and is_exon_node[lsrc]
            _src_mat = _is_bnd[lsrc] and is_exon_node[i]
            _mu_d = mu_proj_mat[i] if _dst_mat else mu_proj[i]
            _mu_s = mu_proj_mat[lsrc] if _src_mat else mu_proj[lsrc]
            _var_d = var_proj_mat[i] if _dst_mat else var_proj[i]
            s2t = _var_d + (_mu_d - _mu_s) ** 2
            # source COHERENT fractions from its (λ,θ) belief (the MODE; the belief is still the point
            # estimate). ``vls`` kept for the diagnostic capture.
            ls, vls = mu_lam[lsrc], var_lam[lsrc]
            sin_t = math.sin(mu_th[lsrc])
            cos_t = math.cos(mu_th[lsrc])
            fg_s = 1.0 / (1.0 + math.exp(-ls))
            fr_s = 1.0 - fg_s
            fp_s = fr_s * (1.0 + sin_t) / 2.0
            fn_s = fr_s * (1.0 - sin_t) / 2.0
            # PRECISION SOURCE (the log-fraction variances): the REFERENCE-FREE evidence variance ``1/τ`` — a
            # structural lock is composition-CERTAIN (``ev=0`` ⇒ high pr, governed only by the honest
            # 1/M_src+s2t); a source with no composition evidence (``τ=0``, not locked) is gated OUT
            # (``lam_ev``/``th_ev`` False ⇒ pr=0), the phantom collapse.
            # The composition evidence VARIANCE on λ=logit(f_g) (reference-free — never the belief variance,
            # which pools the shared reference into a phantom; message_precision_derivation.md §2). Three states,
            # nothing in between (emission_and_precision_derivation.md §2):
            #   structural lock  → composition CERTAIN → variance 0
            #   real evidence τ  → Var(log f_c) = jac²·(1/τ_λ) → finite
            #   NO evidence τ=0  → composition UNSEEN → variance ∞ → the message carries zero composition precision
            # Set ∞ DIRECTLY (not fr_s²·∞, which would be 0·∞=nan at the λ-window edge) — coordinate-clean and safe.
            lock_s = bool(struct_lock[lsrc])
            if lock_s:
                v_logfg = v_logfr = 0.0
            elif tau_lam[lsrc] > _EPS:
                _ev_lam = 1.0 / tau_lam[lsrc]
                v_logfg = fr_s * fr_s * _ev_lam  # Var(log f_g) = (1−f_g)²·(1/τ_λ)
                v_logfr = fg_s * fg_s * _ev_lam  # Var(log f_r) =  f_g²·(1/τ_λ)     (RNA-total)
            else:
                v_logfg = v_logfr = math.inf
            ev_th = (1.0 / tau_th[lsrc]) if tau_th[lsrc] > _EPS else 0.0  # θ tilt evidence (single-strand lock ⇒ 0)
            v_logfp = v_logfr + (cos_t / max(1.0 + sin_t, _EPS)) ** 2 * ev_th  # +θ term (0 for single-strand)
            v_logfn = v_logfr + (cos_t / max(1.0 - sin_t, _EPS)) ** 2 * ev_th
            if _capture is not None:
                _pt["sm"][i], _pt["vlfg"][i], _pt["vlfp"][i] = sm, v_logfg, v_logfp
                _pt["s2t"][i], _pt["fgs"][i], _pt["vls"][i] = s2t, fg_s, vls
                _pt.setdefault("tau_src", np.zeros(n_nodes))[i] = tau_lam[lsrc]
                _pt.setdefault("lock_src", np.zeros(n_nodes))[i] = 1.0 if struct_lock[lsrc] else 0.0
            # STRUCTURAL emission only — the source ALWAYS emits (emission_and_precision_derivation.md §3). The
            # only gates are structural: gDNA flows genomically; each RNA strand flows where that strand is
            # CONTINUOUS (`free_s`) on both endpoints (facing mass — unspliced `sm` or spliced `SPs/SNs`). There
            # is NO evidence gate: a composition-vacuous source (τ=0) carries ``v_log*=∞`` above ⇒ its PREDICTION
            # precision is 0 (``_pred_precision``) ⇒ the recipient ignores it. This RETIRES the τ-gag bug class —
            # a gate keyed on τ once silenced BOTH channels, hiding the spliced MEASUREMENT (52% of junctions on
            # unstranded data); now the measurement is a separate channel that no gate can suppress (§4).
            emit_g = sm > _EPS
            emit_p = fp[lsrc] and fp[i] and (sm > _EPS or SPs[lsrc] > _EPS)
            emit_n = fn[lsrc] and fn[i] and (sm > _EPS or SNs[lsrc] > _EPS)
            # THE MODE PREDICATE (message_propagation_arithmetic.md §4/§7). The eff-length-frame log-odds SHIFT on
            # CLEAN edges (neither endpoint an exon region — the cliff is invariant there); the DENSITY mode
            # (observed-md anchor; the ``±SPs/−absorb`` mature reconciliation rides in ``rho_pos``/``rho_neg``) on
            # exon edges, where a composition-conserving shift over-trusts the unstranded exon (§9).
            # INTRON REGION → BOUNDARY also uses the shift (owner 2026-07-22): an intron and its IE/EI boundary
            # have IDENTICAL active components (gDNA + nascent RNA, mature-free), so composition is invariant
            # across the enrichment cliff and the fraction scales safely — the initial design, restored after the
            # geo-mean §14 refactor scoped every boundary-dst edge onto the density frame. The EXON source →
            # boundary direction (mature-contaminated) still needs its own derivation, so it keeps the geo-mean.
            _ex_s, _ex_d = is_exon_node[lsrc], is_exon_node[i]
            use_shift = (not _ex_s) and (not _ex_d) and (not _is_bnd[i] or is_intron_node[lsrc])
            # EXON↔BOUNDARY: the mature-dilution edges. The gDNA factor rides the composition SHIFT plus the
            # additive ±c_b; the RNA factors keep their existing treatment (the exon's unspliced RNA is
            # mature-contaminated and stays suppressed, §14 change-1), so `use_shift` (which governs the RNA
            # modes) is left untouched and only the gDNA mode switches via ``use_shift_g``.
            _c_mat = _v_mat = 0.0
            if _is_bnd[i] and _ex_s:  # exon → boundary: REMOVE the exon's mature
                _c_mat, _v_mat = mature_dilution[i], _var_mat[i]
            elif _ex_d and _is_bnd[lsrc]:  # boundary → exon: RESTORE the mature
                _c_mat, _v_mat = -mature_dilution[lsrc], _var_mat[lsrc]
            use_shift_g = use_shift or (_c_mat != 0.0) or (_is_bnd[i] and _ex_s) or (_ex_d and _is_bnd[lsrc])
            # ---- source per-component DENSITIES (shared by the mode AND the precision path) ----
            # ``rho_g`` = source gDNA density; ``rho_pos``/``rho_neg`` = the per-strand RNA density the source
            # imputes, carrying the FULL mature accounting (§5, BOTH directions): ``− absorb`` removes the
            # SOURCE exon's contained mature (``SPd/ESPd``, nonzero only when the source-facing flank is an
            # exon) and ``+ SPs/esp`` adds the DST exon's mature (``SPs/ESPs``, nonzero only when the
            # dst-facing flank is an exon). SPd/SPs route to the exon flank only ⇒ at most one fires per hop,
            # so a boundary/intron end contributes no mature. These densities are what the PRECISION path
            # (``n_eff``, ``rho_r``) consumes, unchanged.
            _eg = EGs[lsrc] if EGs[lsrc] > _EPS else _EPS
            _er = ERs[lsrc] if ERs[lsrc] > _EPS else _EPS
            _esp = ESPs[lsrc] if ESPs[lsrc] > _EPS else _EPS
            _espd = ESPd[i] if ESPd[i] > _EPS else _EPS
            rho_g = fg_s * sm / _eg
            absorb_p = SPd[i] / _espd  # source exon mature density (removed from the imputed RNA)
            absorb_n = SNd[i] / _espd
            rho_pos = fp_s * sm / _er + SPs[lsrc] / _esp - absorb_p  # +SPs/esp = the DST exon's mature (added)
            rho_neg = fn_s * sm / _er + SNs[lsrc] / _esp - absorb_n
            if _is_bnd[i]:
                if use_shift:
                    # INTRON → BOUNDARY: composition invariance. Transfer the intron's OWN gDNA + nascent
                    # (identical active components, mature-free) across the enrichment cliff via the shift below —
                    # NO geo-mean crossing (that decodes the boundary from BOTH flanks) and NO mature
                    # reconciliation (the intron carries no mature; the boundary's spliced is not the intron's to
                    # absorb). ``rho_g`` stays the intron's own gDNA density; drop the ``±SPs/−absorb`` terms.
                    rho_pos = fp_s * sm / _er
                    rho_neg = fn_s * sm / _er
                elif _ex_s:
                    # EXON → boundary: the MATURE-DILUTION path. ``rho_g`` stays the exon's OWN gDNA density
                    # (NO geo-mean override — that unweighted hack stood in for this missing imputation); the
                    # +c_b term converts the exon's mature-contaminated composition into the boundary's
                    # mature-free crossing composition. RNA stays suppressed to the spliced MEASUREMENT.
                    rho_pos = SPs[lsrc] / _esp - absorb_p
                    rho_neg = SNs[lsrc] / _esp - absorb_n
                # else: INTERGENIC → boundary (a TSS/TES seam) keeps the source's own densities on the
                # density mode. A seam is a distinct case — the intergenic flank carries gDNA only (no
                # nascent) — so neither composition invariance nor the mature-dilution rule applies
                # unmodified. OUT OF SCOPE pending its own derivation (roadmap R6).
            # ---- the cliff-crossing LOG-ODDS SHIFT (cliff_message_derivation.md §3) ----
            # Impute the dst's f_c as the source COMPOSITION: the per-component imputed MASS ``M_c = ρ_c^src ·
            # E_c^dst``, normalized by the IMPUTED total ``ΣM``. Equivalent to the log-odds shift
            # ``λ_dst = λ_src + log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)``: the capture enrichment ``e(x)``
            # cancels identically (CLIFF-INVARIANT), and the shift is nonzero because gDNA and RNA have different
            # FL distributions ⇒ different per-component eff-lengths in the contained (region) vs crossing
            # (boundary) frames. This REPLACES the retired ``log(ρ_c^src/ρ_total^dst)`` density mode, which
            # divided an enriched-frame numerator by the dst's single-component OBSERVED total ``md/E`` and so
            # failed across the ~10²–10³× cliff (isolated intron→boundary |Δf_g| 0.65 → 0.17). MC-validated
            # across gaussian/gamma/bimodal/uniform FL pairs (scripts/debug/cliff_message_mc.py).
            if use_shift_g:
                # ΣM must count only what PHYSICALLY crosses to the dst: gDNA always (genomic); RNA strand s
                # only where s is structurally continuous on BOTH endpoints (``fp/fn`` — the transcript-
                # structure gate). Without this, a gene-end / TES / opposite-strand seam (where the RNA does
                # NOT cross, so the boundary's unspliced crossing is PURE gDNA — oracle f_g=1) has its gDNA
                # composition wrongly DILUTED by the non-crossing exon RNA, crushing f_g toward 0 (verified:
                # density imputes ≈1 there, ungated composition ≈0.03; 646:14 edges worse). The density mode
                # is immune because it divides by the dst's OBSERVED total (genuinely ≈pure gDNA at a seam).
                cont_p = bool(fp[lsrc]) and bool(fp[i])
                cont_n = bool(fn[lsrc]) and bool(fn[i])
                Mg = rho_g * egd  # → dst gDNA mass (always crosses)
                Mp = (max(rho_pos, 0.0) * erd) if cont_p else 0.0  # → dst +RNA mass (0 if it can't cross)
                Mn = (max(rho_neg, 0.0) * erd) if cont_n else 0.0
                _den = Mg + Mp + Mn
                _den = _den if _den > _EPS else _EPS  # gDNA + crossing-RNA imputed mass
                comp_fl = 1.0 / md  # one-fragment floor ≡ the retired density floor log(1/M_dst) = −log M_dst
            lam_factors = []
            # ---- gDNA message (a factor on log f_g) ----
            if emit_g:
                # gDNA rides the composition SHIFT on every clean edge AND on the exon↔boundary edges, where the
                # additive ±c_b converts between the exon's mature-INCLUSIVE and the boundary's mature-FREE
                # composition (enrichment cancels). ``_v_mat`` carries c_b's counting noise into the precision,
                # so a count-starved junction self-limits instead of asserting an uncorrected mode.
                mo = _mode_shift(Mg, _den, comp_fl) if use_shift_g else _mode_density(rho_g, egd, md)
                mo += _c_mat
                pr = _pred_precision(smn, v_logfg + _v_mat, s2t)  # 1/(Var(log f_g)+Var(c_b)+1/n_src+σ²_T)
                amg[i], apg[i] = mo, pr
                if pr > 0.0:
                    lam_factors.append((True, mo, pr))
            # ---- RNA per-strand messages (emission per strand) + the RNA-TOTAL λ-factor ----
            # The RNA a splice-junction boundary imputes folds TWO sources into one message: the deconvolved
            # UNSPLICED RNA (a PREDICTION at count-zero-info-degraded precision ``pr``) and the SPLICED fragments
            # (a direct MEASUREMENT of mature RNA — already in the ``rho_*`` density via ``±`` above).
            # R2-BLOCKED (measured 2026-07-22): undamping this MEASUREMENT precision is correct in principle
            # (a junction COUNT is not an imputation, so an enrichment cliff should not attenuate it) but it
            # REGRESSED the gdna_none phantom guard +49%. Cause: ``pr += S`` attaches the MEASUREMENT's
            # confidence to the PREDICTION's mode ``mo``, which the mature absorption can drive to a clamped
            # ~zero — undamped, that laundered a weak "no RNA" into a CONFIDENT one (exactly what the clamp
            # note below forbids). σ²_transfer was silently holding this unsound merge together. R2 is
            # therefore blocked on item E (prediction⊕measurement MERGE — the measurement needs its OWN mode).
            # The honest clamp on the PREDICTION side stays: when the
            # absorption saturates the residual to the count floor, that is a genuine but WEAK ~zero — backed by
            # the floor's ONE count, not the source's full count — never laundered into a confident "no RNA".
            rho_r = 0.0
            if emit_p:
                # HYBRID keeps RNA on the PURE density mode (decoupled from the gDNA cliff correction —
                # density's robustness: the gDNA f_g must NOT be hostage to the error-prone mature removal).
                mo = _mode_shift(Mp, _den, comp_fl) if use_shift else _mode_density(rho_pos, erd, md)
                n_eff = (max(rho_pos, 0.0) * erd) if (rho_pos <= 1.0 / erd) else smn  # honest clamp
                pr = _pred_precision(n_eff, v_logfp, s2t)  # PREDICTION channel (0 on a composition-vacuous source)
                if SPs[lsrc] > _EPS:  # MEASUREMENT channel (independent evidence — always credited)
                    pr += SPs[lsrc] / (1.0 + SPs[lsrc] * s2t)  # damped: see R2-BLOCKED note above
                amp[i], app[i] = mo, pr
                if rho_pos > 0.0:
                    rho_r += rho_pos
            if emit_n:
                mo = _mode_shift(Mn, _den, comp_fl) if use_shift else _mode_density(rho_neg, erd, md)
                n_eff = (max(rho_neg, 0.0) * erd) if (rho_neg <= 1.0 / erd) else smn  # honest clamp
                pr = _pred_precision(n_eff, v_logfn, s2t)  # PREDICTION channel (0 on a composition-vacuous source)
                if SNs[lsrc] > _EPS:  # MEASUREMENT channel (independent evidence — always credited)
                    pr += SNs[lsrc] / (1.0 + SNs[lsrc] * s2t)  # damped: see R2-BLOCKED note above
                amn[i], apn[i] = mo, pr
                if rho_neg > 0.0:
                    rho_r += rho_neg
            # RNA-TOTAL message (a factor on log f_r) — the second λ-factor; gDNA + RNA-total in tension on the
            # SAME axis λ (this is what makes the pie coherent — docs §2.2).
            if rho_r > _EPS:
                # PREDICTION channel (τ-gated) + spliced MEASUREMENT channel (independent — always)
                pr_r = _pred_precision(smn, v_logfr, s2t)  # PREDICTION channel (0 on a composition-vacuous source)
                s_spl = SPs[lsrc] + SNs[lsrc]  # total spliced (motif-stranded ⇒ one term nonzero)
                if s_spl > _EPS:
                    pr_r += s_spl / (1.0 + s_spl * s2t)  # damped: see R2-BLOCKED note above
                if pr_r > 0.0:
                    mo_r = _mode_shift(Mp + Mn, _den, comp_fl) if use_shift else _mode_density(rho_r, erd, md)
                    lam_factors.append((False, mo_r, pr_r))
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
                # The τ cavity: the dst's evidence precision grows by the message it absorbed, so it relays
                # that (real) evidence onward — but NEVER the reference (τ starts at I_strand/I_struct, never at
                # the reference precision a=0.1153). CRITICAL FRAME CONVERSION: the message precision ``f[2]``
                # is in the dst's LOG-FRACTION frame (log f_g for a gDNA factor, log f_r for an RNA-total
                # factor); τ lives in the λ=logit(f_g) frame. Convert by the Jacobian² — ``(d log f_g/dλ)² =
                # (1−f_g)²`` for gDNA, ``(d log f_r/dλ)² = f_g²`` for RNA-total — BEFORE accumulating. Without
                # this, τ over-counts by 1/jac² and self-bootstraps: a tiny seed amplifies 1/(1−f_g)²≈4× per hop
                # and saturates at 1/s2t (toy: `scratchpad/tau_toy.py`; real: τ→12 in gdna_none). This
                # conversion makes the relayed λ-evidence ≈ the source's own (s2t-damped), never amplified.
                _fgd = 1.0 / (1.0 + math.exp(-mu_lam[i]))  # dst f_g after the fold
                _jd_g, _jd_r = (1.0 - _fgd) ** 2, _fgd * _fgd
                tau_lam[i] += math.fsum(f[2] * (_jd_g if f[0] else _jd_r) for f in lam_factors)
            if _capture is not None:  # inert (production _capture is None) — Stage-1 mode-verify per-edge record
                _capture.setdefault("_edge_modes", []).append(
                    {
                        "src": int(lsrc), "dst": int(i), "use_shift": bool(use_shift),
                        "ref_s": int(chain.ref_idx[lsrc]), "ref_d": int(chain.ref_idx[i]),
                        "fp_s": bool(fp[lsrc]), "fp_d": bool(fp[i]),
                        "fn_s": bool(fn[lsrc]), "fn_d": bool(fn[i]),
                        "exon_s": bool(is_exon_node[lsrc]), "exon_d": bool(is_exon_node[i]),
                        "kind_s": int(chain.kind[lsrc]), "kind_d": int(chain.kind[i]),
                        "rho_g": float(rho_g), "rho_pos": float(rho_pos), "rho_neg": float(rho_neg),
                        "egs": float(_eg), "ers": float(_er), "egd": float(egd), "erd": float(erd),
                        "md": float(md), "sm": float(sm),
                        # DIAGNOSTIC (inert): the decomposed RNA density terms for the boundary RNA-imputation trace
                        "rna_src": float((fp_s + fn_s) * sm / _er),   # source RNA density (nascent + within-exon mature)
                        "mat_abs": float((SPd[i] + SNd[i]) / _espd),  # boundary junction-mature density (absorbed)
                        "mat_add": float((SPs[lsrc] + SNs[lsrc]) / _esp),  # source-face mature added
                        "esp_s": float(_esp), "esp_d": float(_espd), "fr_s": float(fr_s),
                        "mode_g": float(amg[i]), "mode_p": float(amp[i]), "mode_n": float(amn[i]),
                        "prec_g": float(apg[i]), "prec_p": float(app[i]), "prec_n": float(apn[i]),
                    }
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
            _capture.setdefault("_prec_terms", []).append(
                {**{k: v.copy() for k, v in _pt.items()}, "pr_g": apg.copy(), "pr_p": app.copy()}
            )
            # running-belief (λ) trajectory at the END of this scan — shows how the relay sharpens σ²_λ.
            _capture.setdefault("_relay_lam", []).append((mu_lam.copy(), var_lam.copy()))
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
            # DIAGNOSTIC (inert in production): the composition-evidence seed + the bare enrichment projection
            # (mu_proj, var_proj) that σ²_transfer = var_proj[dst]+(mu_proj[dst]−mu_proj[src])² is built from, so a
            # diagnostic can recompute per-flank σ²_transfer + the DOF verdict offline (boundary-rule tooling).
            _tau0_lam=tau0_lam,
            _mu_proj=mu_proj,
            _var_proj=var_proj,
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
