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
import os
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
    "node_total_density",
    "chain_region_deconv",
    "chain_boundary_side_deconv",
]

_EPS = 1.0e-9

# The variance of λ = logit(f_g) under the node's OWN reference Beta(½,½) (`simplex_logodds._JEFFREYS_REF`) —
# the honest relay seed for a node that observed nothing, so its posterior IS its reference. DERIVED, not
# tuned: ``Var(logit f) = ψ′(a) + ψ′(b) = 2·ψ′(½) = 2·(π²/2) = π²`` (and ``E[λ] = ψ(½) − ψ(½) = 0``, which the
# skipped local solve already returns). Used only by the B1a data-free seed below.
_REF_LAM_VAR = math.pi * math.pi

# B1b (message_layer_derivation.md §12.4/§12.7) — the gDNA DENSITY relay. **LANDED as the default.**
# Measured on all 32 ambig_dense_10mb scenarios vs the oracle: mwae 0.1396 → 0.1361, corr 0.671 → 0.688,
# 16 better / 6 worse / 10 flat, ``gdna_none`` 3,818,168 → 3,704,635 (−3.0 %), and the STRANDED arm neutral
# (0.0209 → 0.0214, 10 of 12 scenarios bit-identical). It also removes the ship blocker: corr(precision,|err|)
# on boundaries flips from +0.09…+0.67 (confidently wrong) to negative, with their message precision dropping
# 0.62–0.90 → 0.02–0.06 — honestly weak instead of loudly wrong.
# ``RIGEL_B1B=0`` restores the pre-B1b path for A/B; "2" additionally makes ``struct_lock`` STRUCTURAL, which
# was measured to be a wash (0.1359 vs 0.1361) because the fusion subsumes it.
_B1B_MODE = os.environ.get("RIGEL_B1B", "1")
_B1B = _B1B_MODE != "0"
# N4a — the per-strand RNA relay at IDENTITY transport. DEFAULT OFF: its premise was measured FALSE.
# `free_s` only says a strand is PRESENT on both flanks, not that it is the SAME transcript, and RNA is
# transcript-specific rather than a genomic background. Measured RNA-density step across the very boundaries
# N4a fires on (stranded scenarios, |log rho_R(right) - log rho_R(left)|): median **2.314** (a 10x step),
# p90 7.53 (~1800x) -- versus gDNA's ~0.02-0.09 (§6.5 E1). Battery: 0.1361 -> 0.1372, 0 better / 4 worse.
# The machinery is kept for N4b, where the transport is MEASURED (the routing split) rather than assumed 1.
_N4A = os.environ.get("RIGEL_N4A", "0") != "0"

# N4b-var (§13.6f/g) — the ROUTING transport on exon→junction edges, WITH the variance the peel needs.
# Under §13.0 there is one RNA; at a junction it splits by ROUTE. The exon sends its whole RNA density and the
# receiving boundary subtracts the departing share IT measured:  rho_cont = rho_R - rho_dep.
# The subtraction was measured to be right in the MEDIAN (0.914 of the oracle-required peel) and
# catastrophically short in a TAIL (p10 0.306; 4.9 % below 0.1) — a PRECISION failure, not a bias. Applied as
# a point estimate it over-states an intron's RNA by ~7x and crushes its f_g (0.97 -> 0.38 on the worst node).
# The delta method on a DIFFERENCE supplies the missing term, with no free parameter:
#     Var(log rho_cont) = (rho_R/rho_cont)^2 * Var(log rho_R) + (rho_dep/rho_cont)^2 * (1/n_spliced)
# The weights DIVERGE as rho_cont -> 0, so a near-total peel yields a near-worthless residual — which is
# exactly the honest statement, and is what stops the intron collapse.
_N4B = os.environ.get("RIGEL_N4B", "0") != "0"
_B1B_STRUCT_LOCK = _B1B_MODE == "2"

# E2 (enrichment_ratio_generalization.md; owner 2026-07-23) — the ENRICHMENT-RATIO REFRAME on the gDNA relay.
# A message carries a DENSITY, and a density is only comparable between two nodes after scaling by their
# enrichment ratio. B1b relays the gDNA density but transports the MODE UNSCALED (only σ²_transfer damps the
# precision), so a density measured in a depleted intron is blended into an enriched exon still claiming the
# intron's value — the "transport across a scale change" §13.6k localized as the real failure. E2 scales the
# incoming running density into the DESTINATION's frame by the belief-free ratio ρframe(dst)/ρframe(src),
# ρframe = M/E_g (validated median-0 across direct edges AND tiny-region relays:
# scripts/debug/enrichment_relay_validation.py). A data-free node has no frame and passes the carried frame
# through (owner: tiny regions relay; the ratio is formed between the framed nodes across the gap). This is a
# MODE fix only — σ²_transfer is untouched here (it is a proxy for this very correction, §B.1 of
# message_arithmetic_reconciliation.md, and its neutralization is a separate, owner-driven step).
_E2 = os.environ.get("RIGEL_E2", "0") != "0"

# UNIFIED (unified_solver_design.md; owner 2026-07-23) — the ONE-MODE solver that retires the shift/density
# split. Every message is REFRAMED into the node's frame by the enrichment ratio r = ρ_tot(node)/ρ_tot(src)
# (ρ_tot lazy + composition-aware, per-side spliced), FILTERED (kill inactive strands at the destination), the
# mature ROUTED (absorb/graft via the boundary's measured spliced), then turned into a per-component factor by
# the density mode ÷ M_dst — which (§2 theorem, test_enrichment_frame.py) equals the shift on a matched message
# and handles a partial source where the shift asserts f_g=1. The forward-backward α/β skeleton + the ψ solve
# stay; only the message construction changes. ρ_tot uses the both-message composition via a two-iteration
# combine. PRECISION is unchanged here (σ²_transfer + 1/(Var+1/n)); the variance of r and retiring σ²_transfer
# are the NEXT task (R2/R3). Default off ⇒ the current ``_scan`` path, byte-identical.
_UNIFIED = os.environ.get("RIGEL_UNIFIED", "0") != "0"
_UNIFIED_RHO_ITERS = int(os.environ.get("RIGEL_UNIFIED_RHO_ITERS", "2"))  # two-iteration lazy ρ_tot (owner)
_UNIFIED_ROUTE = os.environ.get("RIGEL_UNIFIED_ROUTE", "1") != "0"  # mature peel/graft (diagnostic gate)
_UNIFIED_S2T = os.environ.get("RIGEL_UNIFIED_S2T", "1")  # R3 ablation: "1" keep, "flat" drop cliff, "0" retire
# Pin the RELAY's running belief to a composition? OFF: the relay carries the absolute DENSITY currency
# (design §10.5) and only the COMBINE pins (`_pin_v` per message). Pinning per hop cancels the reframe
# identically — every component is scaled by the same ``r``, so ``k = M/Σρ_c·E_c`` divides it straight back
# out — and substitutes the destination's own ZERO-precision local default for uninformed components, which
# on unstranded data is fg_loc ≈ ½ at every node: a self-fulfilling wash that overwrites the anchor.
_UNIFIED_RELAY_PIN = os.environ.get("RIGEL_UNIFIED_RELAY_PIN", "0") != "0"

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


def node_total_density(chain: NodeChain, geometry: NodeGeometry, f_g):
    """The LAZY, composition-aware node total density (`unified_solver_design.md` §2, owner 2026-07-23):
    the SUM of component densities, each in its OWN FL frame, from the current belief ``f_g``::

        ρ_unspliced = f_g · (M/E_g^gDNA)  +  (1−f_g) · (M/E_r^RNA)      gDNA-FL for gDNA, RNA-FL for RNA
        ρ_spliced   = spliced_mass / E_spl^RNA                          one-sided, boundary only

    Returns ``(rho_unspliced, rho_with_spliced)`` per node. ``rho_with_spliced`` adds ρ_spliced and is the
    total density used to form the enrichment ratio toward the exon/acceptor side (mature-bearing); the
    mature-free ``rho_unspliced`` is used toward the intron side (§6). This is NEVER a pure-gDNA precompute —
    ``f_g`` is the best current composition (self-solve + messages + measured spliced); gDNA-FL alone
    (``f_g = 1``) is only the fallback where composition is genuinely unknown, and the bounding lemma (§2)
    bounds *that* fallback, not this. Mass/eff are the node-level (both-face-pooled) quantities of
    :func:`node_global_geometry`; the RNA eff-length is the RNA-FL twin, summed the same way."""
    is_reg = np.asarray(chain.kind) == REGION
    mass, eff_g = node_global_geometry(chain, geometry)
    erl = np.asarray(geometry.eff_rna_left, dtype=np.float64)
    err_ = np.asarray(geometry.eff_rna_right, dtype=np.float64)
    eff_r = np.where(is_reg, erl, erl + err_)
    fg = np.clip(np.asarray(f_g, dtype=np.float64), 0.0, 1.0)
    rho_unspl = mass * (fg / np.maximum(eff_g, _EPS) + (1.0 - fg) / np.maximum(eff_r, _EPS))
    # one-sided spliced (mature) DENSITY: spliced mass lands on ONE face, so divide by THAT face's E_spl (a
    # summed-eff divisor would under-state it ~2×). Sum the per-face densities (only the acceptor face is nonzero).
    spl_l = np.asarray(geometry.spliced_pos_left, np.float64) + np.asarray(geometry.spliced_neg_left, np.float64)
    spl_r = np.asarray(geometry.spliced_pos_right, np.float64) + np.asarray(geometry.spliced_neg_right, np.float64)
    espl_l = np.maximum(np.asarray(geometry.eff_spl_left, np.float64), _EPS)
    espl_r = np.maximum(np.asarray(geometry.eff_spl_right, np.float64), _EPS)
    rho_spliced = np.where(spl_l > _EPS, spl_l / espl_l, 0.0) + np.where(spl_r > _EPS, spl_r / espl_r, 0.0)
    return rho_unspl, rho_unspl + rho_spliced


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


def _lambda_factor_precision(lam_logprior, lam_grid):
    """``I_factory`` — the composition evidence a λ-factor carries, read off its own CURVATURE.

    The gDNA intron factory (`gdna_intron_factory.py`) deconvolves an intron against the intergenic background
    NegBinom and tabulates ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)`` over the σ(λ) solve grid. That curve is a
    genuine, reference-FREE likelihood on ``λ`` — external information (``ρ_bg``) about this node's composition
    — so its precision belongs in ``τ_λ`` alongside ``I_strand`` and ``I_struct``. This is the precision
    `gdna_intron_factory_design.md` §4 derives: *"read straight off the curvature of the tabulated log P(f_g) on
    the solve grid ... by grid moment-matching, so the belief is (mode, precision) with no extra work."*

    Why it is registered here and not as a ``struct_lock``: certainty is binary and unbounded, and extending
    ``struct_lock`` was tried and REVERTED (§E.8 — a lock asserts a floor at raw counting precision and pinned a
    vacuous exon at f_g = 0.9616). The curvature is the honest alternative — ``Var(g) = μ + μ²/α_eff`` makes a
    low-count / high-overdispersion intron peel IMPRECISE and a dense one confident, so the factory self-limits
    exactly where §E.8's certainty did not. A FLAT row (every non-intron node, or an uninformative background)
    returns τ = 0 — no evidence, unchanged behaviour.

    ``τ_λ = 1/Var_λ`` under the normalized factor. Returns ``(m,)``; all-zero when ``lam_logprior`` is None."""
    if lam_logprior is None:
        return None
    lp = np.asarray(lam_logprior, dtype=np.float64)
    lam = np.asarray(lam_grid, dtype=np.float64)
    live = np.ptp(lp, axis=1) > _EPS  # a flat factor carries NO information (τ=0), never the grid's own width
    tau = np.zeros(lp.shape[0], dtype=np.float64)
    if not bool(live.any()):
        return tau
    w = np.exp(lp[live] - np.max(lp[live], axis=1, keepdims=True))
    w /= np.maximum(w.sum(axis=1, keepdims=True), _EPS)
    mu = w @ lam
    var = w @ (lam * lam) - mu * mu
    tau[live] = np.where(var > _EPS, 1.0 / np.maximum(var, _EPS), 0.0)
    return tau


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
    # I_struct — composition CERTAINTY. TRIED AND REVERTED (2026-07-22, Part E.8): making this STRUCTURAL
    # (``struct_lock = locked``, dropping ``& is_region``) is what the gDNA-factory design calls for — a locked
    # node admits NO RNA strand, so its unspliced fragments are necessarily gDNA, and that argument is identical
    # for an intergenic REGION and a G1 BOUNDARY seam. Measured: it un-silences the factory exactly as intended
    # (``lock_src`` 0→1, ``Var(log f_g)`` inf→0 at every seam) and the phantom guard IMPROVES
    # (3,821,731 → 3,814,453). It was reverted because it breaks the phantom RED LINE on the prior-free path:
    # with no enrichment prior ``σ²_transfer = 0``, so a terminal G1 lock asserts its density floor at raw
    # counting precision (n = 47–124) and pins a vacuous exon at f_g = 0.9616
    # (``test_tau_gag_fix_deconvolution_prediction_stays_gated``). The missing piece is NOT the certainty — it
    # is a sound statement of the FLOOR: a source's gDNA DENSITY bounds a destination's gDNA FRACTION only
    # where enrichment is comparable, and ``_mode_density`` across a cliff yields "fractions" of 186–324
    # (measured, boundary→intron). σ²_transfer is what currently makes that floor safe — which REFRAMES R4:
    # the cliff term is not merely a legacy proxy, it is the mechanism holding the density floor together.
    # The old restriction, kept until the floor is derived:
    # B1b mode 2: composition certainty is a property of the COMPONENT SET, not the node kind — a TSS/TES seam
    # admits no RNA strand either, so it is gDNA-certain exactly as an intergenic region is. Rides with the
    # density relay (alone it amplifies the M6 mode defect; §12.7).
    struct_lock = np.asarray(locked, dtype=bool)
    if not _B1B_STRUCT_LOCK:
        struct_lock = struct_lock & np.asarray(is_region, dtype=bool)
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
    """Composition-SHIFT mode ``log(max(M_c/ΣM, comp_fl))`` — cliff-invariant (§4a).

    ⚠ MEASURED DEFECT (`message_layer_derivation.md` §9.1 M5): ``comp_fl = 1/md`` is a "one-fragment floor"
    computed from a fractional MASS, not an integer count, so it EXCEEDS 1 whenever the destination face
    carries less than one fragment — asserting a gDNA *fraction* above 1 (up to 1e9 on an empty face) on
    27–67 % of live edges. The identity is exact: ``frac(exp(mode) > 1) == frac(md < 1)`` per edge family.
    Capping at ``log 1`` was A/B'd and is **inert** (§9.1 M7) — the damage is done by modes at ~1.0, not
    above it — so this is a symptom of the frame error, not the lever. It disappears with the per-face frame.
    """
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


def _routing_precision(rho_r, rho_dep, v_logf, n_src, n_spl, s2t):
    """Precision of the RNA that CONTINUES past a junction: ``rho_cont = rho_r - rho_dep`` (§13.6f).

    A difference of two noisy densities. Its LOG-variance is the delta method with the shares SQUARED, and
    those shares diverge as ``rho_cont -> 0`` — so a near-total peel is honestly near-worthless rather than
    confidently zero, which is precisely the intron case that collapsed without this term::

        Var(log rho_cont) = (rho_r/rho_cont)^2 * [Var(log f_r) + 1/n_src]
                          + (rho_dep/rho_cont)^2 * [1/n_spliced]

    Returns 0 (no information) when nothing continues or either input is unusable — never a confident zero."""
    rho_cont = rho_r - rho_dep
    if rho_cont <= _EPS or rho_r <= _EPS or not math.isfinite(v_logf):
        return 0.0
    w_r, w_d = rho_r / rho_cont, rho_dep / rho_cont
    v_r = v_logf + (1.0 / n_src if n_src > 0.0 else math.inf)
    if not math.isfinite(v_r):
        return 0.0
    v = w_r * w_r * v_r + w_d * w_d * (1.0 / (n_spl + 1.0)) + s2t
    return 1.0 / v if v > 0.0 else 0.0


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
    # the INCOMING belief, kept for the diagnostic capture: it is the ``fg_ref`` the final solve freezes its
    # variance at, so a channel-ablation replay must pass the SAME reference to be faithful (inert in prod).
    _fg_init, _fp_init, _fn_init = f_g.copy(), f_pos.copy(), f_neg.copy()

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
    # ── B1a (message_layer_derivation.md §12): ``locked`` covers TWO DIFFERENT STATES and must not treat them
    # alike when seeding the relay coordinate.
    #
    #   STRUCTURAL gDNA  ``~(fp|fn)``               — no admissible RNA strand (intergenic, a no-RNA-crossing
    #                                                 seam), so the unspliced mass is NECESSARILY gDNA.
    #                                                 μ_λ = logit(f_g) = +L at σ²_λ = 0 is CORRECT: this is a
    #                                                 composition CERTAINTY, and it is what the G1-emission
    #                                                 fix above exists to propagate.
    #   DATA-FREE        ``(fp|fn) & mass == 0``    — admits RNA but observed NOTHING. ``init_beliefs`` already
    #                                                 records this honestly as ``var_gdna = inf``; pinning it
    #                                                 at +L/σ²=0 overrides "no information" with "certainly
    #                                                 pure gDNA". Measured: **790 nodes/scenario, 23.3 % of the
    #                                                 chain** (463 region, 327 boundary).
    #
    # These nodes are SILENT today (``emit_g = sm > _EPS`` and they carry no mass), so this split is expected
    # to be byte-identical — it is the safety PREREQUISITE for B1b, which restores the relay through them. Were
    # B1b landed first, all 790 would begin broadcasting "certainly pure gDNA" at full confidence.
    #
    # A data-free node's honest λ seed is its own REFERENCE, not a lock: with no likelihood the posterior IS
    # Beta(½,½) (`_JEFFREYS_REF`), whose mean and variance in λ = logit(f_g) are DERIVED, not chosen —
    #     E[λ] = ψ(½) − ψ(½) = 0     Var[λ] = ψ′(½) + ψ′(½) = 2·(π²/2) = π²
    # and the skipped local solve already returns λ = 0, so only the variance needs writing.
    _Lw = float(logodds_window)
    lam_locked = np.clip(
        np.log(np.clip(f_g, _EPS, 1.0 - _EPS) / np.clip(1.0 - f_g, _EPS, 1.0 - _EPS)), -_Lw, _Lw
    )
    struct_gdna = ~(np.asarray(fp, bool) | np.asarray(fn, bool))  # composition-CERTAIN by structure
    no_data = locked & ~struct_gdna  # admits RNA, observed nothing → no information
    lam_loc = np.where(struct_gdna, lam_locked, lam_loc)
    lvar_loc = np.where(struct_gdna, 0.0, np.where(no_data, _REF_LAM_VAR, lvar_loc))

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
    # ── A3 (2026-07-23): the DENOMINATOR. Spliced mass lands on ONE face, so it must be divided by THAT
    # face's opportunity — not by the sum of both faces', which counts an opportunity that can never receive
    # mass and so under-states S_B by ~2×. Exact closed-form check (`scripts/debug/cb_denominator_check.py`,
    # ρ_g=0.5 ρ_ν=0.3 ρ_μ=1.0 ⇒ r_true=1.25):
    #       SUMMED-EFF  (previous)  r = 0.714   c_b = 0.5390     ← under-peels by 33 % in log space
    #       PER-FACE    (this)      r = 1.250   c_b = 0.8109     ← EXACT against the truth
    # The old "the per-side density sum over-corrects ~2×" note was a real observation of a DIFFERENT bug:
    # ``eff_spl`` was the half-triangle, which A1/A2 showed is up to 199× low on a continuing face (1.88× at
    # R=100). Summing the denominators cancelled part of that. Two compensating errors; A2 fixed the frame, so
    # the correct per-face denominator can now be used.
    # ⚠ RESIDUAL, NOT fixed here: ``_D_B`` divides a mass containing NASCENT (RNA FL) by the **gDNA**-FL
    # opportunity, so it is not ρ_g+ρ_ν — measured ``D_B/(ρ_g+ρ_ν) = 0.875``, i.e. c_b still over-peels ~14 %.
    # Removing that needs the composition, which is what the solve produces, so it belongs in the restructure.
    #
    # ── DIRECTIONAL (2026-07-23). "Is a splice junction" and "mature crosses contiguously" are INDEPENDENT
    # properties of a boundary, and this term used to conflate them by summing BOTH faces' spliced flux.
    # An exon↔exon boundary can be a splice junction AND carry contiguous mature at the same time:
    #     TA+ (1000,2000)(10000,11000)   TB+ (1000,2000)(9000,11000)
    # The signature changes at 10000, so there is an exon↔exon boundary there; it is ALSO TA's splice
    # ACCEPTOR (2000→10000), while TB's exon runs straight through 9000→11000. So mature BOTH departs
    # (TA, into the spliced channel) and continues unspliced (TB) at the same point.
    #
    # The peel must remove only the mature that TERMINATES on the SOURCE's side. The accumulator already
    # routes spliced mass ONE-SIDED to its exon flank (`node_geometry`, junction_strand × exon bit), so the
    # source-facing face's flux IS that quantity: in the example the acceptor's exon flank is the RIGHT, so
    # right-exon→boundary peels TA's mature and left-exon→boundary peels NOTHING (its face carries zero
    # spliced, and TB's mature demonstrably continues). Summing both faces peeled the full junction flux from
    # BOTH directions.
    #
    # Measured on the battery: 78 boundaries/scenario are exon↔exon AND a splice junction — 4.6 % of
    # boundaries but **19.1 % of boundary mass** (a further 18.2 % is plain exon↔exon, where mature also
    # crosses). Of the 78, 38 carry spliced on ONE face (the flank with zero flux was being falsely peeled)
    # and 40 on BOTH (a donor one side, an acceptor the other — each direction was peeled by the SUM).
    #
    # On a CLASSIC exon↔intron junction this is byte-identical: spliced lands only on the exon flank, and
    # ``c_b`` is applied only on an exon-source edge, so the source-facing face already holds the whole flux.
    _D_B = (MS[0] + MS[1]) / np.maximum(EG[0] + EG[1], _EPS)  # unspliced crossing density (node-level)
    _esp_face = (np.maximum(ESP[0], _EPS), np.maximum(ESP[1], _EPS))  # A3: per-FACE opportunity
    _unspl_n = np.asarray(statics.u_pos, np.float64) + np.asarray(statics.u_neg, np.float64)
    # Counting-noise variance of c_b, propagated through c=log(1+r):  Var(c) = [r/(1+r)]²·Var(log r), and
    # Var(log r) = 1/n_s + 1/n_d for two independent Poisson channels. The variance MUST use the INTEGER
    # COUNTS, not the fractional mass — mass sums per-fragment shares, so 1/mass mis-states the Poisson
    # variance (node_geometry §spliced_n_*: `Var(log ρ_m)=1/n`, Kish n_eff ≥ mass). Gamma posterior (n+1)
    # keeps it finite at n=0. NOTE: this models SAMPLING noise only — a probe-ATTENUATED junction is a
    # systematic BIAS in the mode, which no variance term can represent (see the plan doc §5 risk 1).
    # Indexed by FACE (0=left, 1=right); ``_scan`` selects the face that faces the message SOURCE (``df``).
    mature_dilution: list[np.ndarray] = []
    _var_mat: list[np.ndarray] = []
    for _f in (0, 1):
        _s_face = (SP[_f] + SN[_f]) / _esp_face[_f]  # spliced (mature) density terminating on THIS face
        _r = np.where(_D_B > _EPS, _s_face / np.maximum(_D_B, _EPS), 0.0)
        _n_face = SPN[_f] + SNN[_f]  # integer spliced flux on THIS face
        mature_dilution.append(np.where(~_is_reg_arr, np.log1p(_r), 0.0))
        _var_mat.append(
            np.where(
                ~_is_reg_arr,
                (_r / (1.0 + _r)) ** 2 * (1.0 / (_n_face + 1.0) + 1.0 / (_unspl_n + 1.0)),
                0.0,
            )
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
    # ── I_factory: the gDNA intron factory's own curvature, registered as composition evidence ────────────
    # The factory is the architecture's THIRD information source (the population baseline ρ_bg) made node-local,
    # and it is reference-free — so a factory-solved intron has genuinely learned its λ and must be allowed to
    # SAY so. Without this it learns in silence: `_local_solve` gets `lam_logprior` and the intron's own f_g goes
    # 0.504 → 0.975 (oracle 1.000), but τ_λ stays 0, so Var(log f_g) = ∞, `_own_pg` = 0, and the message it emits
    # carries no precision at all. Measured on gdna300_ss_0.50_nrna_none_capture_on: turning the factory on moved
    # the intron's own belief +93 % while `prec_g` at its neighbours was BIT-IDENTICAL, and the correction
    # attenuated to +16 % at the flanking boundary and +0.6 % at the exon — the factory learning, and the relay
    # never hearing it. τ is the channel that carries "this node knows its composition", which is exactly what
    # the factory establishes (`gdna_intron_factory_design.md` §4).
    _tau_fac = _lambda_factor_precision(intron_prior, _lam_lo)
    if _tau_fac is not None:
        tau0_lam = tau0_lam + _tau_fac

    # ── B1b: the per-node OWN gDNA-density estimate (message-free) — the relay's own-term ──────────────────
    # The message content is today ``rho_g = f_g^src · sm / E_g``, built from the SOURCE's own facing mass. A
    # data-free node then has no density to send at all (§12.1 D1), and a LOW-count node throws away what it
    # just heard and substitutes its own noisy estimate — measured at the intergenic→seam→exon hop, where a
    # 4,114-count message (precision 8.21) is replaced by an 11-count one (precision 0.41) and, under strong
    # capture, by nothing at all (§12.7).
    #
    # Both are the same defect. The fix is one precision-weighted fusion of (own estimate, relayed estimate),
    # in the frame-free DENSITY currency, carried along the scan:
    #     own:  ρ_g = f_g^loc · M_node / E_g^node ,  precision 1/(Var(log f_g^loc) + 1/n_node)
    #     run(i) = fuse( own(i) , transported run(src) )
    # A node with no mass has own precision 0 and therefore ADOPTS the relay exactly — pass-through with no
    # special case. An intergenic region's own term dominates ⇒ it SOURCES (and overwrites, i.e. sinks). A seam
    # blends. The owner's source/sink/relay taxonomy is not three rules; it is this one fusion (§12.7).
    #
    # NODE-level (both faces pooled), not per-face: a boundary's two faces observe the SAME fragments — the
    # split is mass bookkeeping — so pooling is the unbiased, maximally-precise estimate (`node_global_geometry`).
    _is_reg_b1 = np.asarray(chain.kind) == REGION
    _n_node = np.where(
        _is_reg_b1,
        np.asarray(geometry.n_unspl_left, np.float64),
        np.asarray(geometry.n_unspl_left, np.float64) + np.asarray(geometry.n_unspl_right, np.float64),
    )
    # message-free composition variance on log f_g, mirroring the in-scan three-state rule but seeded from
    # tau0 (the node's OWN evidence) rather than the running, message-accumulated tau.
    _fr_loc = 1.0 - np.clip(fg_loc, 0.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        _v_logfg_loc = np.where(
            struct_lock, 0.0, np.where(tau0_lam > _EPS, _fr_loc * _fr_loc / np.maximum(tau0_lam, _EPS), np.inf)
        )
    _own_rho_g = np.where(
        (mass_global > _EPS) & (eff_global > _EPS), fg_loc * mass_global / np.maximum(eff_global, _EPS), 0.0
    )
    _own_lg = np.log(np.maximum(_own_rho_g, _EPS))
    # NB: substitute a finite value BEFORE the product — ``np.where`` evaluates both branches, so leaving the
    # ∞ in would form ``0·∞ = nan`` at every count-free node. Same hazard ``_pred_precision`` avoids by
    # branching rather than masking.
    # ── E2: the per-node ENRICHMENT FRAME = the NPMLE enrichment projection ``mu_proj`` (the same quantity
    # σ²_transfer is built from — ``s2t = var_proj[dst] + (mu_proj[dst] − mu_proj[src])²``). A relayed gDNA
    # density is scaled into a destination's frame by ``mu_proj(dst) − (frame it arrived in)``: σ²_transfer
    # *penalizes* that enrichment step; E2 *corrects* by it (the same measured Δ), which is what lets the
    # penalty be retired later (message_arithmetic_reconciliation.md §B.1). Chosen over the raw ``M/E_g``
    # density because the projection SNAPS each node onto the learned landscape modes, washing out the
    # per-node composition/frame scatter that ``M/E_g`` carries (the region-vs-boundary offset): at capture-off
    # the landscape is a single mode, so every node projects to ≈ the same ``mu_proj`` ⇒ reframe ≈ 0 and no
    # spurious scaling where there is no enrichment. A node with no mass has no frame and passes the carried
    # frame through. Inert unless ``_E2``; and ``mu_proj ≡ 0`` when no enrichment prior ⇒ reframe ≡ 0 (inert).
    _has_frame = mass_global > _EPS
    _log_rho_frame = np.where(_has_frame, np.asarray(mu_proj, np.float64), 0.0)
    # ── The gDNA WALL (owner's SINK model, §12.9). An intergenic region is a LOCUS WALL: it SOURCES its own
    # measured gDNA and SINKS whatever arrives, so a claim from one locus never propagates into the next.
    # This is a deliberate STRUCTURAL choice, not a consequence of continuity (gDNA *is* genomically
    # continuous). Three justifications: (a) ρ_g uniformity is only LOCAL — CNV and mappability break it at
    # range; (b) an intergenic region measures ρ_g directly (f_g ≡ 1, no solve), so it has nothing to learn
    # from inside the locus; (c) it bounds error propagation to one locus.
    #
    # Letting precision decide instead is NOT safe, and that is measured: 11.3 % of intergenic-sourced gDNA
    # messages already differ from the region's own measurement at capture-off and **60.9 % at capture-on**.
    # The magnitude is negligible today only because the intergenic count (n≈4,114) out-shouts the incoming —
    # but that count falls to **n≈6 under very strong capture** (§12.7), where it would lose. The wall is
    # currently absent, merely out-shouted.
    #
    # Conditioned on HAVING a measurement: a mass-free intergenic region has nothing to defend and must stay
    # transparent, or the wall re-severs the chain B1b just reconnected (279 dead structural-gDNA nodes/scenario).
    _gdna_wall = (
        _is_reg_b1 & (np.asarray(_rtype[_ri], np.int64) == 0) & (mass_global > _EPS)
    ).tolist()
    # ── N4a (§13.6): the per-strand RNA relay, IDENTITY transport only ────────────────────────────────────
    # Same fusion as the gDNA relay, but RNA is NOT continuous everywhere: at a splice junction the routing
    # splits it (§13.0 — one RNA, two routes), so the relay is restricted here to edges whose boundary
    # endpoint carries NO spliced channel. There the whole RNA flux continues contiguously and the transport
    # is the identity, exactly as for gDNA. Junction edges (the routing transport) are N4b and are untouched.
    _ern_node = np.where(  # node-level RNA-FL opportunity — the twin of ``eff_global``
        _is_reg_b1,
        np.asarray(geometry.eff_rna_left, np.float64),
        np.asarray(geometry.eff_rna_left, np.float64) + np.asarray(geometry.eff_rna_right, np.float64),
    )
    _spl_node = (
        np.asarray(geometry.spliced_pos_left, np.float64)
        + np.asarray(geometry.spliced_neg_left, np.float64)
        + np.asarray(geometry.spliced_pos_right, np.float64)
        + np.asarray(geometry.spliced_neg_right, np.float64)
    )
    _no_spl_bnd = ((~_is_reg_b1) & (_spl_node <= _EPS)).tolist()  # boundary with no routing split

    def _own_rna(frac_loc, var_loc, live):
        """The message-free own-density estimate for one RNA strand, node-level (both faces pooled)."""
        rho = np.where(
            (mass_global > _EPS) & (_ern_node > _EPS) & live,
            np.asarray(frac_loc, np.float64) * mass_global / np.maximum(_ern_node, _EPS),
            0.0,
        )
        v = np.asarray(var_loc, np.float64)
        v_ok = np.isfinite(v)
        v_fin = np.where(v_ok, v, 0.0)
        pr = np.where(
            (_n_node > 0.0) & v_ok & (rho > _EPS),
            _n_node / (_n_node * np.maximum(v_fin, 0.0) + 1.0),
            0.0,
        )
        return np.log(np.maximum(rho, _EPS)), pr

    _own_lp, _own_pp = _own_rna(fp_loc, vp_loc, np.asarray(fp, bool))
    _own_ln, _own_pn = _own_rna(fn_loc, vn_loc, np.asarray(fn, bool))

    _v_fin = np.where(np.isfinite(_v_logfg_loc), _v_logfg_loc, 0.0)
    _own_pg = np.where(
        (_n_node > 0.0) & np.isfinite(_v_logfg_loc) & (_own_rho_g > _EPS),
        _n_node / (_n_node * np.maximum(_v_fin, 0.0) + 1.0),
        0.0,
    )

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
        # B1b: the running fused gDNA log-density + its precision, seeded from each node's OWN estimate.
        # Because the scan visits genomically, ``rl_g[lsrc]`` already carries everything upstream of the
        # source when the source is read — i.e. it IS the cavity for the message to ``i`` (it excludes i's
        # side), so no message is ever fed back into its own origin.
        rl_g, rp_g = _own_lg.copy(), _own_pg.copy()
        # E2: the frame each running gDNA density is currently expressed in (``log ρframe`` of the last FRAMED
        # node it passed through). A framed node's own seed is in its own frame; a data-free node inherits the
        # carried frame from upstream (set in the fusion below). Reset per scan, like ``rl_g``.
        rlfr_g = np.where(_has_frame, _log_rho_frame, np.nan).astype(np.float64)
        # N4a: the per-strand RNA running state (identity-transport edges only — see the fuse below)
        rl_p, rp_p = _own_lp.copy(), _own_pp.copy()
        rl_n, rp_n = _own_ln.copy(), _own_pn.copy()
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
            # B1b: a node emits gDNA iff its RUNNING density carries precision — which it does when it has its
            # own count OR when it inherited a claim from upstream. This is what reconnects the 565 severed
            # relays/scenario (§12.2); the old ``sm > _EPS`` keyed emission to the node's own mass alone.
            emit_g = (rp_g[lsrc] > 0.0) if _B1B else (sm > _EPS)
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
            # EXON↔BOUNDARY: the mature reconciliation. The two directions are NOT mirror images, because the
            # two endpoints know different things (message_arithmetic_reconciliation.md Part E):
            #   * exon → boundary  = PEEL. The dst's unspliced crossing is mature-FREE, but the exon cannot
            #     decompose its own lumped RNA (a region carries NO spliced channel — `spliced_* = 0` for every
            #     REGION node, node_geometry). Only the DST measures the mature share, as its own S_B/D_B ⇒ the
            #     peel must be an additive SHIFT ``+c_b`` in the dst's frame. No eff-length extrapolation.
            #   * boundary → exon  = GRAFT. The dst HAS mature, and the src MEASURES it (S_B). So the mature
            #     enters as an ordinary per-component density (``+SPs/_esp`` in ``rho_pos``/``rho_neg`` below),
            #     transported by the shared rule ``M_c = ρ_c·E_c^dst``. It is already inside the shift's
            #     normalizer ``_den`` ⇒ applying ``−c_b`` on top would restore the SAME mature TWICE.
            # Measured on the grounded toy (boundary→exon, capture ON): lump-only f_g 0.502, c_b-only 0.528,
            # BOTH (the old code) 0.4615 — strictly below either single accounting, i.e. a systematic
            # gDNA UNDER-call at exons (oracle 0.754). PEEL uses the shift; GRAFT uses a density term; never
            # both on one edge. The graft is also what makes item E well-posed: only a separated ρ_μ/ρ_ν can
            # carry provenance-distinct precision.
            # P0b TRIED AND REJECTED (measured, 2026-07-22). The additive ``+c_b`` moves the gDNA mode ALONE,
            # so the gDNA and RNA λ-factors assert incompatible pies (Part E.5 defect 2). The coherent
            # alternative — peel in the RNA NUMERATOR (``rho_pos = fp_s·sm/_er − absorb_p``) so the correction
            # flows through the SHARED normalizer ``_den`` and f_g + f_r ≡ 1 by construction, with ``c_b``
            # dropped — is WORSE on every axis: gdna_none 3,321,740 → 3,524,635; boundary corr(prec,|err|)
            # −0.205 → +0.033 at capture off (the honesty inversion partly returns); AMBIG-region corr
            # unchanged (0.078 → 0.072). Likely cause: ``absorb_*`` is a density in the SPLICED eff-length
            # frame (``spliced_side_eff_length``) while ``rho_pos`` is in the CONTAINED RNA frame, and that
            # spliced frame is known-wrong on a face where the exon continues (effective_length.py). Do not
            # retry the numerator peel before that frame is fixed.
            _c_mat = _v_mat = 0.0
            _n4b_edge = False
            if _is_bnd[i] and _ex_s and not _N4B:  # exon → boundary: PEEL the exon's mature
                # DIRECTIONAL: the face of the destination boundary that FACES this source (``df``). Only the
                # mature terminating on the source's own side is peeled — an exon↔exon boundary that is also
                # a splice junction carries contiguous mature through it, which must NOT be peeled.
                _c_mat, _v_mat = mature_dilution[df][i], _var_mat[df][i]
            # GATE-EQUALITY (Part E.7) — DERIVED, prototyped, reverted with the factory unlock it pairs with.
            # The composition SHIFT transfers the source's whole pie, so it is valid only where the source can
            # supply EVERY component the destination admits. Where a structural gate zeroes a strand the dst
            # DOES admit — a TSS/TES seam, an opposite-strand exon↔exon junction — the source is BLIND to real
            # destination RNA and the shift asserts "the dst is pure gDNA" (measured: f_g = 1.0000 emitted into
            # an exon whose oracle is 0.754). The predicate is
            #     gate_unequal = (fp[dst] and not cont_p) or (fn[dst] and not cont_n)
            # and it is what would DEMOTE the density mode to its correct role (the one-sided gDNA floor at a
            # seam) rather than retiring it. Blocked on the same missing piece as the factory: a density floor
            # is only sound where enrichment is comparable (see `_compile_strand_evidence`).
            use_shift_g = use_shift or (_is_bnd[i] and _ex_s) or (_ex_d and _is_bnd[lsrc])
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
            _b1b_pg = 0.0
            if _B1B:
                # B1b (§12.4): the SOURCE's FUSED running density is the message content, and its precision
                # replaces the per-face ``_pred_precision(smn, v_logfg, ·)`` because the composition + counting
                # terms already live inside ``rp_g`` (see ``_own_pg``). σ²_transfer and Var(c_b) are added per
                # edge below.
                # E2: REFRAME the source's running density into the DESTINATION's enrichment frame before it
                # is used (as the message content ``rho_g`` AND fused below). A density measured at ``src``'s
                # enrichment ``e(src)`` predicts ``dst``'s as ``ρ(src)·e(dst)/e(src) = ρ(src)·ρframe(dst)/
                # ρframe(src)``. ``rlfr_g[lsrc]`` is the frame the incoming density is in (its own if framed,
                # the last framed node's if it relayed through data-free nodes). A dst with no frame does not
                # reframe — it passes the density through in the frame it arrived in (``_reframe = 0``). Inert
                # unless ``_E2`` (``_reframe = 0`` ⇒ byte-identical to B1b).
                _reframe = 0.0
                if _E2 and _has_frame[i] and math.isfinite(rlfr_g[lsrc]):
                    _reframe = _log_rho_frame[i] - rlfr_g[lsrc]
                _rl_src = rl_g[lsrc] + _reframe
                if rp_g[lsrc] > 0.0:
                    # Message content frame. On a SHIFT edge the mode normalizes the imputed masses
                    # (``Mg/ΣM``), and that cancellation requires EVERY component in the source's frame — so
                    # reframing gDNA alone would break it; keep ``rho_g`` in the source frame (``rl_g[lsrc]``,
                    # which the fusion holds in lsrc's own frame). On a DENSITY edge ``rho_g`` is divided by the
                    # dst's OBSERVED ``md``, which carries ``e(dst)``, so the source density must be reframed
                    # into the dst frame (``_rl_src``). The fusion below always uses the reframed ``_rl_src``.
                    rho_g = math.exp(rl_g[lsrc] if (_E2 and use_shift_g) else _rl_src)
                    _b1b_pg = rp_g[lsrc]
                else:
                    rho_g = 0.0
                # RELAY: fuse the transported claim into the DESTINATION's running state, so a node with no
                # mass (own precision 0) passes it on unchanged and a low-count node blends rather than
                # overwrites. σ²_transfer therefore ACCUMULATES per hop. The intergenic WALL (§12.9) is the one
                # node that refuses the fusion: it sources its own measured gDNA and sinks what arrives.
                _mp = 1.0 / (1.0 / _b1b_pg + s2t) if _b1b_pg > 0.0 else 0.0
                if _mp > 0.0 and not _gdna_wall[i]:
                    if rp_g[i] > 0.0:
                        _tp = rp_g[i] + _mp
                        rl_g[i] = (rp_g[i] * rl_g[i] + _mp * _rl_src) / _tp
                        rp_g[i] = _tp
                    else:
                        rl_g[i], rp_g[i] = _rl_src, _mp
                    # E2: rl_g[i] is now expressed in i's own frame if i is framed; a data-free node inherits
                    # the frame the density arrived in, so a framed receiver downstream reframes across the gap.
                    if _E2:
                        rlfr_g[i] = _log_rho_frame[i] if _has_frame[i] else rlfr_g[lsrc]
            if _N4A:
                # N4a: IDENTITY-transport RNA relay. Only where the edge's boundary endpoint has no spliced
                # channel — there the whole RNA flux continues contiguously, so transport is the identity and
                # the fusion is exactly the gDNA one. A junction edge splits the flux by ROUTE (§13.0) and is
                # left to N4b; its densities below are unchanged.
                _b = i if _is_bnd[i] else lsrc
                if _no_spl_bnd[_b]:
                    for _live, _rl, _rp in ((fp, rl_p, rp_p), (fn, rl_n, rp_n)):
                        if not (_live[lsrc] and _live[i]):
                            continue  # SINK: the strand is not continuous across this edge (§12.9)
                        _mpr = 1.0 / (1.0 / _rp[lsrc] + s2t) if _rp[lsrc] > 0.0 else 0.0
                        if _mpr <= 0.0:
                            continue
                        if _rp[i] > 0.0:
                            _tpr = _rp[i] + _mpr
                            _rl[i] = (_rp[i] * _rl[i] + _mpr * _rl[lsrc]) / _tpr
                            _rp[i] = _tpr
                        else:
                            _rl[i], _rp[i] = _rl[lsrc], _mpr
            absorb_p = SPd[i] / _espd  # source exon mature density (removed from the imputed RNA)
            absorb_n = SNd[i] / _espd
            rho_pos = fp_s * sm / _er + SPs[lsrc] / _esp - absorb_p  # +SPs/esp = the DST exon's mature (added)
            rho_neg = fn_s * sm / _er + SNs[lsrc] / _esp - absorb_n
            _n4a = _N4A and _no_spl_bnd[i if _is_bnd[i] else lsrc]
            if _n4a:  # N4a: the FUSED per-strand density replaces the source's own facing estimate
                if fp[lsrc] and fp[i] and rp_p[lsrc] > 0.0:
                    rho_pos = math.exp(rl_p[lsrc])
                if fn[lsrc] and fn[i] and rp_n[lsrc] > 0.0:
                    rho_neg = math.exp(rl_n[lsrc])
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
                    # RNA stays suppressed to the spliced MEASUREMENT. TWO alternatives were implemented and
                    # MEASURED on all 32 scenarios, and BOTH are net-negative — the flag is RETIRED rather than
                    # left as a knob (message_layer_derivation.md §11):
                    #   P0  transport the exon's own RNA, peel via the additive +c_b   mwae 0.1416 (base 0.1396)
                    #   P0b transport it, peel in the RNA NUMERATOR, c_b dropped       mwae 0.1443
                    # Both help every LOW-gDNA scenario and hurt every HIGH-gDNA and every STRANDED one: the
                    # exon's within-exon MATURE arrives at the boundary counted as nascent, biasing f_g down.
                    # A3 (the per-face c_b denominator) closed 55 % of P0's gap; the residual is the mixed-FL
                    # frame error in D_B, which cannot be removed without the composition — Stage C, not a patch.
                    if _N4B:  # ROUTING: send the exon's whole RNA; the boundary subtracts what departs
                        rho_pos = fp_s * sm / _er - absorb_p
                        rho_neg = fn_s * sm / _er - absorb_n
                        _n4b_edge = True
                    else:
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
                if _B1B:  # precision already carries Var(log f_g)+1/n via the fused running density
                    pr = 1.0 / (1.0 / _b1b_pg + _v_mat + s2t) if _b1b_pg > 0.0 else 0.0
                else:
                    pr = _pred_precision(smn, v_logfg + _v_mat, s2t)  # 1/(Var(log f_g)+Var(c_b)+1/n+σ²_T)
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
                if _n4b_edge:  # N4b-var: variance of a DIFFERENCE of two noisy densities (§13.6f)
                    pr = _routing_precision(fp_s * sm / _er, absorb_p, v_logfp, smn, SPN[df][i], s2t)
                else:
                    n_eff = (max(rho_pos, 0.0) * erd) if (rho_pos <= 1.0 / erd) else smn  # honest clamp
                    pr = _pred_precision(n_eff, v_logfp, s2t)  # PREDICTION channel (0 if vacuous source)
                if SPs[lsrc] > _EPS:  # MEASUREMENT channel (independent evidence — always credited)
                    pr += SPs[lsrc] / (1.0 + SPs[lsrc] * s2t)  # damped: see R2-BLOCKED note above
                amp[i], app[i] = mo, pr
                if rho_pos > 0.0:
                    rho_r += rho_pos
            if emit_n:
                mo = _mode_shift(Mn, _den, comp_fl) if use_shift else _mode_density(rho_neg, erd, md)
                if _n4b_edge:
                    pr = _routing_precision(fn_s * sm / _er, absorb_n, v_logfn, smn, SNN[df][i], s2t)
                else:
                    n_eff = (max(rho_neg, 0.0) * erd) if (rho_neg <= 1.0 / erd) else smn  # honest clamp
                    pr = _pred_precision(n_eff, v_logfn, s2t)  # PREDICTION channel (0 if vacuous source)
                if SNs[lsrc] > _EPS:  # MEASUREMENT channel (independent evidence — always credited)
                    pr += SNs[lsrc] / (1.0 + SNs[lsrc] * s2t)  # damped: see R2-BLOCKED note above
                amn[i], apn[i] = mo, pr
                if rho_neg > 0.0:
                    rho_r += rho_neg
            # RNA-TOTAL message (a factor on log f_r) — the second λ-factor; gDNA + RNA-total in tension on the
            # SAME axis λ (this is what makes the pie coherent — docs §2.2).
            if rho_r > _EPS:
                # PREDICTION channel (τ-gated) + spliced MEASUREMENT channel (independent — always)
                if _n4b_edge:
                    # N4b-var: the RNA-TOTAL factor is the one in tension with gDNA on the λ axis, so it is
                    # the one that sets f_g — it needs the SAME difference-variance as the per-strand factors.
                    # Fixing only the per-strand channels left this factor asserting the routing residual at
                    # undiminished confidence, and the stranded regression was unmoved (0.0261 -> 0.0260).
                    pr_r = _routing_precision(
                        fr_s * sm / _er, absorb_p + absorb_n, v_logfr, smn,
                        SPN[df][i] + SNN[df][i], s2t,
                    )
                else:
                    pr_r = _pred_precision(smn, v_logfr, s2t)  # 0 on a composition-vacuous source
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

    # ─────────────────────────────────────────────────────────────────────────────────────────────────────
    # THE UNIFIED SOLVER (unified_solver_design.md; owner 2026-07-23) — ONE mode: reframe → filter → route → ÷M.
    # ─────────────────────────────────────────────────────────────────────────────────────────────────────
    _uni_msg = (None, None, None, None, None, None)  # the unified path's published imputation factors

    def _unified_solve():
        n = f_g.shape[0]
        is_reg_a = np.asarray(chain.kind) == REGION
        is_bnd_a = ~is_reg_a
        ex_a = np.asarray(is_exon_node, dtype=bool)  # source-is-exon selector for the mature routing
        fp_a, fn_a = np.asarray(fp, bool), np.asarray(fn, bool)
        M = np.asarray(mass_global, np.float64)
        E_g = np.asarray(eff_global, np.float64)
        E_r = np.where(is_reg_a, ER[0], ER[0] + ER[1]).astype(np.float64)  # node-level RNA-FL eff-length

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
        spl_p_f = tuple(np.where(SP[k] > _EPS, SP[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
        spl_n_f = tuple(np.where(SN[k] > _EPS, SN[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
        accept_l = (SP[0] + SN[0]) > _EPS  # the LEFT face carries the spliced (acceptor) ⇒ WITH-spliced ρ_tot
        accept_r = (SP[1] + SN[1]) > _EPS

        # own per-component densities + precisions (the message-free SELF-SOLVE — reused from the B1b/N4a seeds)
        og = np.maximum(np.asarray(_own_rho_g, np.float64), 0.0)
        op = np.where(np.asarray(_own_pp, np.float64) > 0.0, np.exp(np.asarray(_own_lp, np.float64)), 0.0)
        on = np.where(np.asarray(_own_pn, np.float64) > 0.0, np.exp(np.asarray(_own_ln, np.float64)), 0.0)
        pg_own = np.asarray(_own_pg, np.float64)
        # ── COHERENT own precision on BOTH arms (the unified relay's own components are independent, so an
        # incoherent pair is not self-correcting the way ``_scan``'s single λ axis is).
        # ``_own_pg`` is built from the REFERENCE-FREE evidence ``tau0_lam``; ``_own_pp``/``_own_pn`` are built
        # from ``vp_loc``/``vn_loc`` — the local solve's POSTERIOR variance, which pools the shared reference
        # (the phantom `message_precision_derivation.md` §2 forbids as a message precision). On unstranded data
        # that asymmetry is total: τ=0 ⇒ Var(log f_g)=∞ ⇒ the gDNA arm is silent at EVERY non-locked node
        # (measured p_own_g = 0.00), while the RNA arm still asserts "f_r = ½" at precision 0.35. Since ``_pin``
        # normalizes the pie to the node's mass, that one-sided claim pushes f_g down everywhere — the residual
        # wash. Rebuild the RNA precision from the SAME τ, with ``_scan``'s own Jacobian: Var(log f_r) = f_g²/τ_λ
        # (vs Var(log f_g) = (1−f_g)²/τ_λ). τ=0 ⇒ BOTH arms silent, which is the honest unstranded statement.
        _fg_l = np.clip(fg_loc, 0.0, 1.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            _v_logfr_loc = np.where(
                struct_lock, 0.0, np.where(tau0_lam > _EPS, _fg_l * _fg_l / np.maximum(tau0_lam, _EPS), np.inf)
            )
        _vr_fin = np.where(np.isfinite(_v_logfr_loc), _v_logfr_loc, 0.0)

        def _own_rna_prec(rho):
            return np.where(
                (_n_node > 0.0) & np.isfinite(_v_logfr_loc) & (rho > _EPS),
                _n_node / (_n_node * np.maximum(_vr_fin, 0.0) + 1.0),
                0.0,
            )

        pp_own, pn_own = _own_rna_prec(op), _own_rna_prec(on)

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
        # contributes nothing, correctly.
        def _pin(i, g, p, n, pg_, pp_, pn_):
            sg = g if pg_ > 0.0 else og[i]
            sp = p if pp_ > 0.0 else op[i]
            sn = n if pn_ > 0.0 else on[i]
            s = sg * E_g[i] + (sp + sn) * E_r[i]
            if s <= _EPS or M[i] <= _EPS:
                return g, p, n
            k = M[i] / s
            return g * k, p * k, n * k

        def _pin_v(g, p, n, pg_, pp_, pn_):  # the vectorized twin (the combine's per-message pin)
            sg = np.where(pg_ > 0.0, g, og)
            sp = np.where(pp_ > 0.0, p, op)
            sn = np.where(pn_ > 0.0, n, on)
            s = sg * E_g + (sp + sn) * E_r
            k = np.where((s > _EPS) & (M > _EPS), M / np.maximum(s, _EPS), 1.0)
            return g * k, p * k, n * k

        def _rho_faces(fgc):
            """Lazy, composition-aware ρ_tot from the current f_g, split per side (WITH spliced at the acceptor)."""
            ru, rw = node_total_density(chain, geometry, fgc)
            rs = rw - ru  # the one-sided spliced density
            return ru, ru + np.where(accept_l, rs, 0.0), ru + np.where(accept_r, rs, 0.0)  # node, left-face, right-face

        # ── the RELAY: accumulate the forward/backward context belief (densities), reframed each hop by the
        # INPUT-belief ρ_tot. Returns each node's context belief IN ITS OWN FRAME; the combine re-reframes.
        rho_node0, rho_l0, rho_r0 = _rho_faces(np.asarray(f_g, np.float64))

        _mup, _vp = np.asarray(mu_proj, np.float64), np.asarray(var_proj, np.float64)
        if _UNIFIED_S2T == "0":  # R3 ablation: retire σ²_transfer entirely (mode-correct ⇒ no cliff penalty)
            _mup, _vp = np.zeros_like(_mup), np.zeros_like(_vp)
        elif _UNIFIED_S2T == "flat":  # R3 ablation: keep the floor var_proj, drop the cliff-height double-count
            _mup = np.zeros_like(_mup)

        def _damp(p, s2t):  # σ²_transfer per-hop damping (PRECISION unchanged this task — R2/R3 replace it)
            return 1.0 / (1.0 / p + s2t) if p > 0.0 else 0.0

        def _relay(seq, nbr, dst_face, src_face, df, sf):
            rg, rp, rn = og.copy(), op.copy(), on.copy()
            pg, pp, pn = pg_own.copy(), pp_own.copy(), pn_own.copy()
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
                r = (dst_face[i] / rho_src) if (rho_src > _EPS and dst_face[i] > _EPS) else 1.0  # no frame ⇒ pass-through
                s2t = _vp[i] + (_mup[i] - _mup[s]) ** 2  # enrichment-crossing damping (bounds relay confidence)
                # GRAFT (boundary → EXON, §6): the boundary's measured mature is a density AT THE SOURCE, so it
                # joins the source's RNA BEFORE the reframe; the peel is measured at the destination and so is
                # applied after. Both use the face that FACES the other endpoint (``sf``/``df``), never the
                # node-pooled sum. Only an EXON receives the graft — an intron carries no mature (`ex_a[i]`,
                # not `is_reg_a[i]`, which grafted the junction's whole mature flux into every flanking intron).
                _gr = _UNIFIED_ROUTE and ex_a[i] and is_bnd_a[s]
                gp = spl_p_f[sf][s] if _gr else 0.0
                gn = spl_n_f[sf][s] if _gr else 0.0
                tg, tp, tn = rg[s] * r, (rp[s] + gp) * r, (rn[s] + gn) * r
                tpg, tpp, tpn = _damp(pg[s], s2t), _damp(pp[s], s2t), _damp(pn[s], s2t)
                # The grafted mature is a MEASUREMENT (a junction COUNT), not an imputation, so it carries its
                # own precision and is NOT τ-gated — the source's PREDICTION precision is 0 on unstranded data
                # and would otherwise drop the graft on the floor. Gating both channels on τ is the τ-gag bug
                # `_scan` retired (it silenced 52 % of spliced junctions unstranded); keep them separate.
                if _gr:
                    tpp += SP[sf][s] / (1.0 + SP[sf][s] * s2t) if SP[sf][s] > _EPS else 0.0
                    tpn += SN[sf][s] / (1.0 + SN[sf][s] * s2t) if SN[sf][s] > _EPS else 0.0
                if _UNIFIED_ROUTE and is_bnd_a[i] and ex_a[s]:  # EXON → boundary: PEEL the departing mature
                    tp, tn = max(tp - spl_p_f[df][i], 0.0), max(tn - spl_n_f[df][i], 0.0)
                if not fp_a[i]:
                    tp, tpp = 0.0, 0.0
                if not fn_a[i]:
                    tn, tpn = 0.0, 0.0
                rg[i], pg[i] = _fuse(og[i], pg_own[i], tg, tpg)
                rp[i], pp[i] = _fuse(op[i], pp_own[i], tp, tpp)
                rn[i], pn[i] = _fuse(on[i], pn_own[i], tn, tpn)
                if _UNIFIED_RELAY_PIN:
                    rg[i], rp[i], rn[i] = _pin(i, rg[i], rp[i], rn[i], pg[i], pp[i], pn[i])
            return rg, rp, rn, pg, pp, pn

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
            rg, rp, rn, pg, pp, pn = fwd_arrs
            # A node with no frame (no mass ⇒ no ρ_tot, §5) cannot reframe: the message passes through at r=1.
            # Falling back to ``rho_src = 1.0`` instead made r the destination's ABSOLUTE density (10³ on a
            # short node) — a raw scale masquerading as a ratio. The relay already guards this way.
            framed = valid & (src_face_v[src] > _EPS) & (dst_face_v > _EPS)
            r = np.where(framed, dst_face_v / np.maximum(src_face_v[src], _EPS), np.where(valid, 1.0, 0.0))
            # GRAFT before the reframe (a density measured AT the source); PEEL after (measured at the dst).
            # Both per-FACE, and the graft only into an EXON — see the relay's twin.
            graft = _UNIFIED_ROUTE & ex_a & is_bnd_a[src] & valid
            gp = np.where(graft, spl_p_f[sf][src], 0.0)
            gn = np.where(graft, spl_n_f[sf][src], 0.0)
            tg, tp, tn = rg[src] * r, (rp[src] + gp) * r, (rn[src] + gn) * r
            # σ²_transfer per-hop damping (vectorized) — same as the relay; PRECISION unchanged this task.
            s2t = _vp + (_mup - _mup[src]) ** 2

            def _dv(p):
                return np.where(valid & (p[src] > 0.0), 1.0 / (1.0 / np.maximum(p[src], _EPS) + s2t), 0.0)

            tpg, tpp, tpn = _dv(pg), _dv(pp), _dv(pn)
            # the graft's MEASUREMENT precision — never τ-gated (see the relay's twin)
            _sp, _sn = np.where(graft, SP[sf][src], 0.0), np.where(graft, SN[sf][src], 0.0)
            tpp = tpp + np.where(_sp > _EPS, _sp / (1.0 + _sp * s2t), 0.0)
            tpn = tpn + np.where(_sn > _EPS, _sn / (1.0 + _sn * s2t), 0.0)
            if _UNIFIED_ROUTE:
                peel = is_bnd_a & ex_a[src] & valid  # EXON → boundary
                tp = np.where(peel, np.maximum(tp - spl_p_f[df], 0.0), tp)
                tn = np.where(peel, np.maximum(tn - spl_n_f[df], 0.0), tn)
            tp, tpp = np.where(fp_a, tp, 0.0), np.where(fp_a, tpp, 0.0)
            tn, tpn = np.where(fn_a, tn, 0.0), np.where(fn_a, tpn, 0.0)
            tg, tp, tn = _pin_v(tg, tp, tn, tpg, tpp, tpn)  # the message is a claim about THIS node's mass
            return tg, tp, tn, tpg, tpp, tpn

        def _fuse_v(a, pa, b, pb):
            p = pa + pb
            return np.where(p > _EPS, (pa * a + pb * b) / np.maximum(p, _EPS), 0.0), p

        dc_fin = None
        f_cur = np.asarray(f_g, np.float64).copy()
        for _ in range(max(1, _UNIFIED_RHO_ITERS)):
            _, rho_lf, rho_rf = _rho_faces(f_cur)
            ag, ap, an, apg, app, apn = _transport(sl, vl, 0, 1, fwd, rho_lf, rho_rf)  # left msg: dst face 0, src face 1
            bg, bp, bn, bpg, bpp, bpn = _transport(sr, vr, 1, 0, bwd, rho_rf, rho_lf)  # right msg: dst face 1, src face 0
            cg, cpg = _fuse_v(ag, apg, bg, bpg)
            cp, cpp = _fuse_v(ap, app, bp, bpp)
            cn, cpn = _fuse_v(an, apn, bn, bpn)
            mo_g = np.log(np.maximum(cg * E_g / np.maximum(M, _EPS), _EPS))
            mo_p = np.log(np.maximum(cp * E_r / np.maximum(M, _EPS), _EPS))
            mo_n = np.log(np.maximum(cn * E_r / np.maximum(M, _EPS), _EPS))
            dc_fin = _local_solve(global_lp, mo_g, cpg, (mo_p, mo_n), (cpp, cpn))
            f_cur = np.clip(np.asarray(dc_fin.gdna_frac, np.float64), 0.0, 1.0)
            nonlocal _uni_msg
            _uni_msg = (mo_g, cpg, mo_p, cpp, mo_n, cpn)  # publish for the shared diagnostics
            if _capture is not None:  # inert diagnostic: the fused per-component densities + the frames
                _capture.setdefault("_uni", []).append(
                    {
                        "cg": cg.copy(), "cp": cp.copy(), "cn": cn.copy(),
                        "pg": cpg.copy(), "pp": cpp.copy(), "pn": cpn.copy(),
                        "ag": ag.copy(), "ap": ap.copy(), "an": an.copy(),
                        "bg": bg.copy(), "bp": bp.copy(), "bn": bn.copy(),
                        "apg": apg.copy(), "app": app.copy(), "bpg": bpg.copy(), "bpp": bpp.copy(),
                        "rho_lf": rho_lf.copy(), "rho_rf": rho_rf.copy(),
                        "mo_g": mo_g.copy(), "mo_p": mo_p.copy(), "mo_n": mo_n.copy(),
                        "fg_out": f_cur.copy(),
                    }
                )
        if _capture is not None:
            _capture.update(
                _uni_static={
                    "M": M, "E_g": E_g, "E_r": E_r, "spl_p": spl_p, "spl_n": spl_n,
                    "og": og, "op": op, "on": on, "pg_own": pg_own, "pp_own": pp_own, "pn_own": pn_own,
                    "fwd_g": fwd[0], "fwd_p": fwd[1], "fwd_n": fwd[2],
                    "fwd_pg": fwd[3], "fwd_pp": fwd[4], "fwd_pn": fwd[5],
                    "bwd_g": bwd[0], "bwd_p": bwd[1], "bwd_n": bwd[2],
                    "bwd_pg": bwd[3], "bwd_pp": bwd[4], "bwd_pn": bwd[5],
                    "rho_node0": rho_node0, "rho_l0": rho_l0, "rho_r0": rho_r0,
                    "is_bnd": is_bnd_a, "is_exon": ex_a, "left": li, "right": ri,
                    "fp": fp_a, "fn": fn_a,
                }
            )
        return dc_fin

    if _UNIFIED:
        dc_fin = _unified_solve()
        a = b = None  # the α/β raw-scan tuples have no unified twin (the relay carries densities, not modes)
        # The unified path's final per-component imputation factors, published into the SAME capture slots the
        # default path uses (`_uni_msg`, set at the end of `_unified_solve`). This is what lets the shared
        # diagnostics — `pass0_node_dissect.py`'s ψ channel-ablation replay above all — run on BOTH solvers.
        mode_g, prec_g, mode_p, prec_p, mode_n, prec_n = _uni_msg
    else:
        a = _scan(order_list, left, 1, 0)  # forward (α: left context)
        b = _scan(order_list[::-1], right, 0, 1)  # backward (β: right context)

        # (D) combine α⊗β (precision-weighted product) per component → one batched FINAL solve.
        def _comb(am_a, ap_a, am_b, ap_b):
            pc = ap_a + ap_b
            return np.where(pc > _EPS, (ap_a * am_a + ap_b * am_b) / np.maximum(pc, _EPS), 0.0), pc

        mode_g, prec_g = _comb(a[0], a[1], b[0], b[1])
        mode_p, prec_p = _comb(a[2], a[3], b[2], b[3])
        mode_n, prec_n = _comb(a[4], a[5], b[4], b[5])
        # (D) FINAL solve with the FB messages (backend-dispatched). The final solve is UNCHANGED — it consumes
        # the per-component messages `_comb` produced (now COHERENT); the relay only changed how they are built.
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
