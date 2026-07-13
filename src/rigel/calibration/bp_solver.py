"""The belief-propagation calibration solver: per-node gDNA-vs-RNA deconvolution by a bidirectional sweep.

Deconvolves each node's unspliced fragment mass into a pie ``(f_pos, f_neg, f_g)`` — sense-RNA /
antisense-RNA / gDNA — over the unified region↔boundary chain (`node_chain`), by a single forward-backward
(L→R then R→L) belief-propagation pass (exact on the chain, a forest of linear paths). Each per-node solve
(`simplex_logodds`, the log-density log-odds solver) reconciles three sources of information: the intrinsic
strand likelihood (the Beta-Binomial tilt — the only signal a count carries), the cross-node imputation
messages, and the population gDNA prior. Theory: the count-zero-information principle in
`docs/calibration/CALIBRATION_ARCHITECTURE.md`; the message-precision model (belief-free Poisson
disagreement-variance ``σ²_msg = σ²_imp + 1/n_src``) in `docs/calibration/disagreement_shrinkage_prior_design_v2.md`.

Module layout. The per-node geometry / belief / density / statics / init primitives (`build_node_geometry`,
`node_densities`, `build_node_statics`, `init_beliefs`, `NodeGeometry`/`NodeBelief`/`NodeStatics`) now live in
the lower `node_geometry` module and are re-exported here for the calibrator's convenience; this module owns
the global gDNA prior + the sweep:
* `_gdna_seed_estimate` — the ANCHORED global gDNA prior: the population baseline
  ρ_global + seed between-node spread σ²_g (the seed-based non-circular firewall). Fit ONCE before the solve
  (NOT a per-message reliability — the message precision is the separately-fit belief-free σ²_imp).
* `node_sweep` — the single forward-backward sweep. Message precision is the **belief-free Poisson
  disagreement-variance** ``pr = n_src/(n_src·σ²_imp + 1)`` (``σ²_imp`` = the empirical adjacent-node total-density
  imputation floor, `adjacent_disagreement_variance`, fit once) — bounded so a confident source can never send an
  overconfident message, and there is no var~mean fixed point / no outer loop.
* `chain_region_deconv` / `chain_boundary_side_deconv` — project the converged belief back to the per-region
  / per-boundary-side masses the `CalibrationResult` consumes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .node_chain import BOUNDARY, REGION, NodeChain
from .node_geometry import (
    NodeBelief,
    NodeDensities,
    NodeGeometry,
    NodeStatics,
    _node_region_type,
    build_node_geometry,
    build_node_statics,
    init_beliefs,
    node_densities,
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
    "NodeDensities",
    "node_densities",
    "NodeStatics",
    "build_node_statics",
    "init_beliefs",
    "node_sweep",
    "chain_region_deconv",
    "chain_boundary_side_deconv",
]

_EPS = 1.0e-9
_MSG_PSEUDOCOUNT = 1.0  # one pseudo-observation: a density from a finite count can never have zero
#                         sampling variance, so a message precision can never escape to ∞ (var stabilizer).
# ── The pass-1 global gDNA prior is STABILITY-ONLY: it gives low/zero-gDNA nodes a finite baseline so the
#    solve stays sane, capped at one pseudo-observation (`_GLOBAL_STAB_PREC`) so it can NEVER drive a
#    solution. A node resolves from its STRAND likelihood + the messages (strong strand pins the tilt;
#    unstranded relies on the messages) and, in PASS 2, the trained Phase-2 gDNA-density KDE — the only
#    trained prior (it supersedes the earlier enrichment-transfer ê(z), now removed). ──
_GLOBAL_STAB_PREC = 1.0  # one pseudo-observation — the global can never override a node's own data.


def _strand_discriminability(kappa: float) -> float:
    """Strand-discriminability weight ``(2κ−1)²`` ∈ [0,1]: 0 at κ=½ (unstranded — no strand information),
    1 at κ=1 (fully stranded). Down-weights strand-derived terms where the strand is uninformative."""
    return float((2.0 * float(kappa) - 1.0) ** 2)


def _gdna_seed_estimate(
    chain, statics, geometry, region_arrays, boundary_substrate, f_g_init, kappa
):
    """The honest, NON-CIRCULAR population gDNA prior, fit ONCE on gDNA-clean seed nodes (§4.3).

    Returns ``(rho_global: float, sigma2_g: _LogLinearVarMean, var_mean: float)`` — the exposure-pooled gDNA
    rate + the deterministic gDNA between-node spread σ²_g(μ) + its rate-estimate variance. Inputs are
    belief-independent (structural ``M/E`` + the strand-ONLY init
    ``f_g``), so this is computed once before the sweep and never refit (breaks the per-pass circularity).

    Seeds (per-node gDNA density + weight):
      * **structural** (always; the only path for UNSTRANDED data) — intergenic & intron regions and
        intergenic-exon & intron-exon boundary crossings (exon-facing side). Density = observed ``M/E``
        (gDNA-clean by structure under the **nascent-sparse** assumption). Weight 1.
      * **strand** (single-strand nodes not already structural — mainly exons; reach the capture-enriched
        range). Density = strand-deconv ``f_g_init·M/E``. Weight ``(2κ−1)²`` (the strand discriminability) so
        the seeds fade to 0 as κ→½ — stranded data gets the extra coverage, unstranded falls back to the
        structural seeds automatically (no hard threshold)."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    is_reg = kind == REGION
    is_bnd = kind == BOUNDARY
    node_rtype, rtype = _node_region_type(chain, region_arrays)
    R = rtype.shape[0]
    EGl, EGr = geometry.eff_gdna_left, geometry.eff_gdna_right
    Ml, Mr = geometry.mass_left, geometry.mass_right

    # boundary flank types → clean (intron/intergenic)-exon boundary + which side is the exon.
    blr = np.asarray(boundary_substrate.left_region, dtype=np.int64)
    brr = np.asarray(boundary_substrate.right_region, dtype=np.int64)
    Bn = blr.shape[0]
    bi_ = np.clip(idx, 0, Bn - 1)
    lt = np.where(
        (blr[bi_] >= 0) & is_bnd, rtype[np.clip(blr[bi_], 0, R - 1)], -1
    )  # left flank type
    rt = np.where(
        (brr[bi_] >= 0) & is_bnd, rtype[np.clip(brr[bi_], 0, R - 1)], -1
    )  # right flank type
    left_clean = (lt == 0) | (lt == 1)
    right_clean = (rt == 0) | (rt == 1)
    clean_exon_bnd = is_bnd & ((left_clean & (rt == 2)) | (right_clean & (lt == 2)))
    exon_on_right = clean_exon_bnd & (rt == 2)

    # representative (mass, gDNA-eff): region = contained; clean-exon boundary = exon-facing side.
    mass = np.where(is_reg, Ml, np.where(exon_on_right, Mr, Ml))
    eff = np.maximum(np.where(is_reg, EGl, np.where(exon_on_right, EGr, EGl)), _EPS)
    rho_obs = mass / eff  # observed total density on the representative face (M/E)

    # Structural (M/E) seeds: intergenic regions (pure gDNA) + clean (intron/intergenic)-exon boundary
    # crossings (the gDNA-clean, capture-OBSERVABLE enriched-edge signal). Intron REGION interiors are NOT
    # structural — their contained mass carries the nascent RNA (so M/E is not gDNA-clean) and they are
    # DEPLETED under capture; a single-strand intron is instead strand-deconvolved below, and the gDNA-clean
    # intron signal the user intends is the intron-exon BOUNDARY crossing (already structural).
    struct_seed = (is_reg & (node_rtype == 0)) | clean_exon_bnd
    strand_seed = np.asarray(statics.strand_obs, bool) & ~struct_seed & (mass > 0.0)
    w_strand = _strand_discriminability(kappa)  # ∈[0,1]; →0 at κ→½
    dens = np.where(
        struct_seed, rho_obs, np.where(strand_seed, np.asarray(f_g_init, float) * rho_obs, 0.0)
    )
    seed_w = np.where(struct_seed, 1.0, np.where(strand_seed, w_strand, 0.0))
    # Structural seeds are kept even at ZERO count (an intergenic region with 0 gDNA fragments over a large
    # exposure E is the STRONGEST evidence that gDNA is scarce — it drives ρ_global→0). Strand seeds need
    # counts to deconvolve, so they alone require mass>0.
    is_seed = struct_seed | strand_seed

    # exposure-pooled gDNA rate with ONE pseudocount of TOTAL gDNA (a=1; the global, not per-node — the
    # Poisson–Gamma posterior Gamma(a+G, E_tot)): G = Σ(w·gcount) (gcount = dens·E = mass for structural,
    # f_g·mass for strand), E_tot = Σ(w·E). ρ_global=(a+G)/E_tot (→1/E_tot≈0⁺ when G=0: gDNA scarce, never
    # zero); var_mean=(a+G)/E_tot² is the rate-estimate variance (1/CV²=a+G floors at the 1 pseudocount).
    sw = seed_w * is_seed
    G = float(np.sum(sw * dens * eff))
    E_tot = max(float(np.sum(sw * eff)), _EPS)
    rho_global = (1.0 + G) / E_tot
    # var_mean = the global rate-estimate variance, log-density form (D4, the delta method = this design's own
    # 1/count principle): Var(log ρ̂) ≈ Var(ρ̂)/ρ̂² = [(1+G)/E_tot²]/[(1+G)/E_tot]² = 1/(1+G) — the inverse
    # pooled effective count (NOT the opaque trigamma; transparent + consistent with pois_log=1/(count+1)).
    # Gentle N→~0.5 at zero-gDNA (G=0).
    var_mean = 1.0 / (1.0 + G)
    return rho_global, _fit_seed_varmean(chain, dens, eff, is_seed, seed_w), var_mean


@dataclass(frozen=True, slots=True)
class _LogLinearVarMean:
    """Deterministic ``σ²_g(μ) = max(a + b·log μ, 0)`` — the closed-form WLS replacement for the bistable
    ``MonotoneVarMean`` P-spline (whose GCV-λ ``argmin`` was the calibration cross-process nondeterminism
    root cause; `calibrate_cross_process_nondeterminism.md`). Fit on the seed-edge Poisson-corrected excess
    ``(raw − offset)`` against the edge-midpoint log-density; every quantity is a continuous algebraic function
    of the data (weighted sums, one ratio, a ``max``) so a machine-ε input change moves ``σ²_g`` by machine-ε —
    NO discrete selection (no ``argmin`` / GCV / IRLS / active-set). It captures the load-bearing monotone μ
    dependence (σ²_g≈0 at low ρ_global ⇒ strong prior suppresses no-gDNA exon FP; large at high ρ_global under
    capture ⇒ weak prior spares enriched exons) that a single scalar cannot. Design:
    `docs/calibration/M2_loglinear_sigma2g_design.md`."""

    a: float
    b: float

    @classmethod
    def fit(cls, means, raws, offs, wts) -> "_LogLinearVarMean":
        m = np.asarray(means, dtype=np.float64)
        y = np.asarray(raws, dtype=np.float64) - np.asarray(
            offs, dtype=np.float64
        )  # Poisson-corrected excess
        w = np.asarray(wts, dtype=np.float64)
        ok = (m > 0.0) & (w > 0.0) & np.isfinite(y)
        m, y, w = m[ok], y[ok], w[ok]
        if m.shape[0] < 2:
            return cls(0.0, 0.0)  # <2 edges → σ²_g=0 (the strong-prior / no-gDNA regime)
        x = np.log(m)
        sw = float(np.sum(w))
        xbar = float(np.sum(w * x) / sw)
        ybar = float(np.sum(w * y) / sw)
        sxx = float(np.sum(w * (x - xbar) ** 2))
        if (
            sxx <= _EPS
        ):  # no density spread (all seeds ~one level) → flat law: b=0, a=weighted-mean excess
            return cls(max(ybar, 0.0), 0.0)
        b = float(np.sum(w * (x - xbar) * (y - ybar)) / sxx)
        a = ybar - b * xbar
        return cls(a, b)

    def predict(self, x) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        return np.maximum(self.a + self.b * np.log(np.maximum(x, _EPS)), 0.0)


def _fit_seed_varmean(chain, dens, eff, is_seed, seed_w):
    """σ²_g(μ) on adjacent SEED-SEED edges from the per-node seed gDNA density (the gDNA between-node spread),
    in the log-density form (twin of :func:`_edge_varmean`). Edge weight = ``min`` of the two endpoints' seed
    weights (the weaker endpoint limits reliability). The fit is the deterministic closed-form log-linear law
    (:class:`_LogLinearVarMean`) — NOT the retired bistable P-spline."""
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    is_seed = np.asarray(is_seed, bool)
    means, raws, offs, wts = [], [], [], []
    for nbr in (left, right):
        idx = np.where((nbr >= 0) & is_seed)[0]
        if not idx.size:
            continue
        s = nbr[idx]
        keep = is_seed[s]
        idx, s = idx[keep], s[keep]
        dr, sr, de, se = dens[idx], dens[s], eff[idx], eff[s]
        ok = (dr > 0.0) & (sr > 0.0)
        dr, sr, de, se = dr[ok], sr[ok], de[ok], se[ok]
        means.append(0.5 * (dr + sr))
        raws.append((np.log(dr) - np.log(sr)) ** 2)
        offs.append(1.0 / (dr * de + 1.0) + 1.0 / (sr * se + 1.0))
        wts.append(np.minimum(seed_w[idx][ok], seed_w[s][ok]))
    cat = lambda p: np.concatenate(p) if p else np.zeros(0)  # noqa: E731
    return _LogLinearVarMean.fit(cat(means), cat(raws), cat(offs), cat(wts))


def _floor_estimate(chain, geometry, region_arrays, f_g_init, kappa):
    """The CONSERVATIVE gDNA background FLOOR — the population background gDNA density pooled over ALL
    intergenic + ALL intron REGIONS (NOT boundaries, NOT exons). Applied by :func:`_global_logprior`: it pins a
    floor-level intron to ``f_g≈1`` and gives an intron whose density SUBSTANTIALLY EXCEEDS the floor its
    nascent (``target = ρ_floor·E/M``) — the nascent-from-self-evidence principle, no gate. Enriched exons are
    excluded (capture-enriched on-target; their gDNA density is handled by the Phase-2 density KDE).

    Nascent RNA is sparse, rare, locus-dependent and unpredictable (there is no library-level nascent rate), so
    the floor is CONSERVATIVE — assume nascent absent and set a HIGH bar, catching only BLATANT nascent and
    conservatively leaving subtle nascent as gDNA (low sensitivity, high specificity; the intended behavior on
    real data, NOT an error). Each region's density is a continuous **strand-weighted** blend
    ``dens_g = (M/E)·gdna_frac``, ``gdna_frac = (1−w)·1 + w·f_g_init`` with ``w = (2κ−1)²`` — the
    strand-specificity weight (0 unstranded, 1 fully stranded), NOT a binary gate:
      * UNSTRANDED (κ=½ ⇒ w=0): ``gdna_frac=1`` ⇒ full observed density ``M/E`` — where nascent is undetectable
        we conservatively count all mass as gDNA background (the high, assume-absent floor).
      * STRANDED (κ→1 ⇒ w→1): ``gdna_frac=f_g_init`` ⇒ down-weight a strand-identified RNA intron OUT of the
        gDNA floor, so the floor DEFERS to the strand signal and stranded nascent is still recovered.
    (Validated 2026-07-08: dropping the weighting — raw pooled density — is identical on unstranded data but
    fights the strand on stranded data, regressing stranded-nascent recovery ~0.4% on the AMBIG suite. The
    strand weighting is kept.) Intergenic regions are locked ``f_g_init=1`` ⇒ their full ``M/E`` always counts.

    Returns ``(rho_floor, s2_floor, var_mean_floor, floor_mask, s2_bg)``: the exposure-pooled floor rate, the
    between-region log-density SPREAD (the floor tightness — biological/CNV excess over the per-region Poisson
    floor), the rate-estimate variance ``1/(1+G)``, the per-chain-node intergenic+intron mask, and the pooled
    background log-density spread ``s2_bg`` (the far-tail threshold for the intron density likelihood).
    """
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    node_rtype, _ = _node_region_type(chain, region_arrays)
    floor_mask = is_reg & ((node_rtype == 0) | (node_rtype == 1))  # intergenic + intron REGIONS
    EGl = np.maximum(
        np.asarray(geometry.eff_gdna_left, np.float64), _EPS
    )  # region: contained gDNA eff-len
    Ml = np.asarray(geometry.mass_left, np.float64)  # region: contained mass
    # strand-WEIGHTED gDNA density: down-weight known-RNA introns where the strand is informative (w→1);
    # falls back to the full observed density where it is not (w→0). Continuous, not a binary gate.
    w_str = _strand_discriminability(kappa)
    gdna_frac = (1.0 - w_str) + w_str * np.asarray(
        f_g_init, np.float64
    )  # ∈ (0,1]; =1 unstranded/intergenic
    dens_g = (Ml / EGl) * gdna_frac
    eff = EGl[floor_mask]
    g_dens = dens_g[floor_mask]
    # Exposure-pooled floor rate (Poisson–Gamma, 1 pseudocount on the TOTAL gDNA): zero-gDNA depleted
    # regions are KEPT (0 gDNA over a large E is the strongest evidence gDNA is scarce → ρ_floor→0⁺).
    G = float(np.sum(g_dens * eff))
    E_tot = max(float(np.sum(eff)), _EPS)
    rho_floor = (1.0 + G) / E_tot
    var_mean_floor = 1.0 / (1.0 + G)

    # Between-region SPREAD of log gDNA-density over the POPULATED floor regions (eff-weighted population
    # variance minus the per-region log-Poisson floor → the excess/biological spread, ≥0). Tight for a
    # coherent depleted population (a confident floor); naturally widens on real data (GC/mappability).
    def _logdens_spread(dvals, evals):
        """Eff-weighted between-node log-density spread minus the per-node log-Poisson floor (the excess/
        biological spread, ≥0); 0.0 for <2 populated nodes."""
        p = dvals > 0.0
        if int(np.sum(p)) < 2:
            return 0.0
        lr_ = np.log(dvals[p])
        w_ = evals[p]
        sw_ = float(np.sum(w_))
        mu_ = float(np.sum(w_ * lr_) / sw_)
        s2_raw_ = float(np.sum(w_ * (lr_ - mu_) ** 2) / sw_)
        pois_ = float(np.mean(1.0 / (dvals[p] * evals[p] + 1.0)))
        return max(s2_raw_ - pois_, 0.0)

    s2_floor = _logdens_spread(g_dens, eff)
    # Background log-density spread for the intron density likelihood = the intergenic + intron floor spread.
    # It must capture the gDNA density variation RELEVANT TO INTRONS (GC / length / mappability / capture
    # off-target) — intergenic-only is too tight, so ordinary intron gDNA variation would read as far-tail
    # "excess" (phantom nascent). Sparse nascent barely inflates it, and the strand weighting down-weights
    # stranded-RNA introns out of it. So ``s2_bg = s2_floor``.
    s2_bg = s2_floor
    return rho_floor, s2_floor, var_mean_floor, floor_mask, s2_bg


def _global_logprior(
    fgg,
    mass_global,
    eff_global,
    rho_global,
    sigma2_g,
    var_mean,
    *,
    floor_mask=None,
    rho_floor=None,
    s2_floor_total=None,
    s2_bg=None,
):
    """Precompute the count-space global as a per-node BINOMIAL pseudo-count on f_g, ``(n_nodes, P)`` (§4):

        glob(f_g) = α·log f_g + (N − α)·log(1 − f_g),   α = N·μ,   μ = clip(ρ·E/M, 0, 1),   N = ρ² / σ²_prior

    The genome-wide baseline uses ``ρ = ρ_global`` and ``σ²_prior = var_mean + σ²_g(ρ_global)``; the
    depleted intergenic/intron floor override (``floor_mask``/``rho_floor``) pins that population to its
    confident floor.

    ``N`` is the **M-INDEPENDENT** population confidence (``= 1/CV²`` of the rate) — a hyperprior is naturally
    imprecise, so it can NEVER overrule a node's own (strand) evidence; the MEAN ``μ`` keeps the honest
    density→fraction ``E/M`` conversion. ``σ²_prior`` = the rate-estimate variance ``var_mean = (1+G)/E_tot²``
    (which carries the 1-pseudocount floor, so ``N → 1`` at uniform zero-gDNA — a gentle zero-baseline, never
    ``0/0``) PLUS the between-node spread ``σ²_g`` (large under capture ⇒ ``N → ρ²/σ²_g ≈`` small, imprecise;
    small under uniform present gDNA ⇒ ``N`` large, confident). Two-sided (mode μ); applied to ALL solvable
    nodes (the strand_obs fork is dissolved)."""
    s2_between = max(float(sigma2_g.predict(np.array([max(rho_global, _EPS)]))[0]), 0.0)
    s2_flat = var_mean + s2_between
    # ── LOG-DENSITY global: a Gaussian on log f_g (D-plan §1.4). var_mean = 1/(1+G) (D4); the
    #    M-independent precision is N_log = 1/Var(log ρ) directly — NO ρ² factor (s2_flat is ALREADY a
    #    log-variance, not a density variance). target = log(implied gDNA fraction). ──
    eff = np.maximum(eff_global, _EPS)
    mass = np.maximum(mass_global, _EPS)
    n_glob = 1.0 / max(s2_flat, _EPS)
    target = np.log(np.clip(rho_global * eff / mass, _EPS, 1.0))
    n_node = np.full(target.shape, n_glob, dtype=np.float64)
    # DEPLETED-REGION FLOOR (`_floor_estimate`): intergenic + intron REGION nodes are a coherent depleted
    # population, so override them with the CONFIDENT floor (ρ_floor at the tight floor spread). target =
    # ρ_floor·E/M ⇒ a floor-level intron pins to f_g≈1; an intron with density EXCESS over the floor gets
    # nascent (low f_g) — the nascent-from-self-evidence principle, no gate. Exons/boundaries keep the
    # genome-wide baseline. On real data the floor and a node's own strand evidence agree (same physical
    # density); they conflict only in the documented all-RNA-intron case the floor assumption excludes.
    if floor_mask is not None and rho_floor is not None and s2_floor_total is not None:
        fm = np.asarray(floor_mask, bool)
        target[fm] = np.log(np.clip(float(rho_floor) * eff[fm] / mass[fm], _EPS, 1.0))
        n_node[fm] = 1.0 / max(float(s2_floor_total), _EPS)
    # STABILITY-ONLY: cap the WHOLE global (flat + floor) at one pseudo-observation so it cannot drive or
    # drag any node — the single-strand solve is carried by strand + messages, the global only keeps
    # low/zero-gDNA nodes finite. (The target VALUES — floor for depleted, ρ_global elsewhere — stay; only
    # their weight is capped.)
    n_node = np.minimum(n_node, _GLOBAL_STAB_PREC)
    # PHASE 2c — the intron DENSITY LIKELIHOOD (`three_component_mature_nascent_design.md` §10.3-10.4). For the
    # depleted floor nodes (intergenic + pure introns — off-target, so the intergenic background is a clean
    # control), the density-vs-background term is a LIKELIHOOD (the node's count vs the population), NOT the
    # capped hyperprior, so its precision is DATA-DRIVEN `1/(σ²_bg + 1/N)` (N = the node's unspliced mass) and
    # is applied AFTER the stability cap (it overrides the clamp for these nodes only). It is BOUNDED — `≤ N`
    # (Poisson) and `≤ 1/σ²_bg` (finite empirical spread) — never infinite, and cannot collapse (σ²_bg, N are
    # fixed); the outgoing MESSAGE precision stays separately ceilinged by σ²_imp in `_scan`, so no
    # overconfident message can leave a node regardless of this local certainty. (`node_sweep` always supplies
    # `floor_mask`/`s2_bg`, so this density likelihood is active for floor nodes; the guards keep the no-floor
    # path — used by the genome-wide-baseline-only callers, e.g. scripts/debug/node_error_attribution.py — valid.)
    log_fg = np.log(np.maximum(np.asarray(fgg, np.float64), _EPS))  # (K,)
    term = -0.5 * n_node[:, None] * (log_fg[None, :] - target[:, None]) ** 2
    if floor_mask is not None and s2_bg is not None:
        fm = np.asarray(floor_mask, bool)
        n_node[fm] = 1.0 / (float(s2_bg) + 1.0 / mass[fm])
        # The density term is a background Gaussian on log ρ_g whose MODE is f_g = ρ_bg/ρ. Alone, that pins any
        # intron ABOVE the background mean below 1 — reading ordinary gDNA density variation as nascent. The
        # scale-free RNA parsimony (Jeffreys ``−log(1−f_g)``, the same term `_kde_logprior` carries) corrects
        # this: gDNA is the residual after a typical-magnitude RNA, so the BULK (within σ²_bg) pins to f_g≈1 and
        # only a genuine far-tail density EXCESS becomes nascent (`three_component_mature_nascent_design.md`
        # §10.3). The floor nodes use THIS density likelihood in place of the capped floor, so the term is
        # rebuilt for them (mode + data-driven precision + parsimony).
        jeff = -np.log1p(-np.minimum(np.asarray(fgg, np.float64), 1.0 - _EPS))  # (K,)
        term[fm, :] = (
            -0.5 * n_node[fm, None] * (log_fg[None, :] - target[fm, None]) ** 2 + jeff[None, :]
        )
    return term


# Lattice density for the tabulate-and-interpolate KDE evaluation in _kde_logprior: this many lattice
# nodes per KDE bandwidth. Linear-interpolation error scales as ~(1/PTS_PER_BW)^2, so 16 keeps it far
# below the kernel's own smoothing — a numerical-accuracy knob, not a model parameter.
_KDE_LATTICE_PTS_PER_BW = 16


def _kde_logprior(fgg, mass_global, eff_global, gdna_prior):
    """The GENERATIVE two-density prior term ``(n_nodes, K)`` on the f_g solve grid (design:
    ``ambig_boundary_spliced_deconvolution.md``; derived by the density-prior-integration workflow).

    The node's total density ``d = M/E`` splits into gDNA density ``ρ_g = f_g·d`` and RNA density
    ``ρ_r = (1−f_g)·d``. Two independent density priors:

      * gDNA: the empirical population KDE ``P(log ρ_g)`` — evaluated with **real Gaussian tails**
        (:meth:`GdnaDensityPrior.logpdf_kernel`, NOT the clamped interpolation, whose constant tail lets a
        high-density node drift to ``f_g≈0.5`` → false-positive gDNA).
      * RNA: a scale-free **Jeffreys** prior ``p(ρ_r) ∝ 1/ρ_r`` (RNA spans >10⁴× — no informative scale).
        Its Jacobian into the f_g coordinate is exactly ``1/(1−f_g)`` ⇒ the ``−log(1−f_g)`` term below.

    The Jeffreys term is the crux: without it gDNA is priored but RNA is free, so the cheapest explanation of
    any flat-strand node is "dump mass into free RNA, park gDNA at the tall depleted KDE mode" — the cliff and
    the false-positives. With it, lowering f_g raises ρ_r (penalised), so gDNA is the residual after a
    typical-magnitude RNA. Both terms are on the f_g grid only (the Jeffreys is node-independent, O(K)).
    NO tuned constants (the Jeffreys exponent is 1, the KDE coordinate is native-log)."""
    eff = np.maximum(np.asarray(eff_global, np.float64), _EPS)
    mass = np.maximum(np.asarray(mass_global, np.float64), _EPS)
    log_me = np.log(mass) - np.log(eff)  # (m,) = log(M/E)
    fg = np.minimum(np.maximum(np.asarray(fgg, np.float64), _EPS), 1.0 - _EPS)  # (K,)
    log_rho = np.log(fg)[None, :] + log_me[:, None]  # (m,K) = log ρ_g at each grid point
    # log ρ_g lies on a bounded 1-D interval and logpdf_kernel is a smooth 1-D function, so evaluate the
    # EXACT kernel (real quadratic tails) on a dense lattice spanning the query range and linearly
    # interpolate the m·K points off it — O(L·n_train + m·K) instead of O(m·K·n_train), which at genome
    # scale (m·K ~ 10^8) is the difference between milliseconds and a multi-TiB OOM. The lattice spans the
    # FULL query range, so no point is extrapolated: the real-tail behaviour that makes logpdf_kernel (not
    # the clamped logpdf) the correct choice is preserved to interpolation accuracy (≪ bandwidth).
    flat = log_rho.ravel()
    lo = float(np.min(flat))
    hi = float(np.max(flat))
    h_bw = max(float(gdna_prior.bandwidth), _EPS)
    # lattice step = bandwidth / PTS_PER_BW ⇒ interp error ~(1/PTS_PER_BW)^2, far below the kernel smoothing.
    n_lat = int(np.clip(np.ceil((hi - lo) / h_bw * _KDE_LATTICE_PTS_PER_BW) + 1.0, 256, 65536))
    if hi - lo <= _EPS or flat.size <= n_lat:
        # Small query set (or degenerate range): the exact kernel is no costlier than building the lattice,
        # so evaluate it directly — keeps small/golden cases BIT-EXACT. Tabulation only engages when the
        # query set is large enough that interpolation genuinely saves work (the genome-scale m·K ~ 10^8).
        kde_term = gdna_prior.logpdf_kernel(flat).reshape(log_rho.shape)
    else:
        lattice = np.linspace(lo, hi, n_lat)
        tab = gdna_prior.logpdf_kernel(lattice)  # (n_lat,) exact kernel — cheap
        kde_term = np.interp(flat, lattice, tab).reshape(log_rho.shape)
    # MIXTURE BRIDGE (Fix 1; ε = ``gdna_prior.mixture_bridge`` from CalibrationConfig; ε=0 ⇒ bit-exact
    # legacy KDE). The KDE is estimated from clean (unimodal) REGION nodes, so it has a deep VALLEY between
    # the depleted and enriched modes — but a node's current-belief gDNA density is generically a spatial
    # MIXTURE of enriched (in-probe) and depleted (off-target) positions (a capture boundary, a sparse-probe
    # region — a whole density BAND, not just boundaries), which lands in that valley by construction. The
    # valley penalty (~10² nats) makes such a node collapse to f_g≈0 (park ρ_g at the tall depleted mode,
    # dump the rest into free RNA), emitting a pathologic RNA message that then crushes its neighbours. Mixing
    # in a uniform "any-mixture" bridge over the observed gDNA-density support floors the valley (no collapse)
    # while leaving the KDE's real Gaussian tails OUTSIDE the support intact (the high-ρ_g false-positive
    # suppression is unchanged, since the bridge is bounded to [blo,bhi]≈[depleted,enriched]). The KDE input
    # is the current-belief gDNA density (refined by strand/spliced/messages), NEVER total density. See
    # docs/calibration/boundary_kde_valley_collapse_and_simplex_precision.md.
    eps = float(getattr(gdna_prior, "mixture_bridge", 0.0) or 0.0)
    if eps > 0.0:
        tx = np.asarray(gdna_prior.train_x, np.float64)
        trim = float(
            getattr(gdna_prior, "bridge_trim_pct", 0.5)
        )  # % support trim (config.calib_kde_bridge_trim_pct)
        blo, bhi = (
            (float(np.percentile(tx, trim)), float(np.percentile(tx, 100.0 - trim)))
            if tx.size
            else (0.0, 0.0)
        )
        if bhi > blo:
            uni = np.where((log_rho >= blo) & (log_rho <= bhi), -math.log(bhi - blo), -np.inf)
            kde_term = np.logaddexp(math.log1p(-eps) + kde_term, math.log(eps) + uni)
    jeffreys = -np.log1p(-fg)  # (K,) RNA Jeffreys 1/(1−f_g)
    return kde_term + jeffreys[None, :]


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
    disagreement_sigma2: float,
    emit_locked: bool = False,
    _capture: dict | None = None,
):
    """The belief-propagation sweep over the chain — gDNA AND per-strand RNA messages, in COUNT space (no
    ``(M/E)²`` Jacobian) — as a single **FORWARD-BACKWARD** pass.

    The chain is a forest of linear paths, so BP is exact in one forward + one backward pass (vs Gauss-Seidel /
    Jacobi which propagate one hop per pass). The message precision is the belief-free Poisson
    disagreement-variance ``pr = n_src/(n_src·σ²_imp + 1)`` (``σ²_imp`` fit once by
    :func:`adjacent_disagreement_variance`), NOT a fitted ``σ²_bio(μ)`` var~mean curve — so there is no precision
    to refit and no outer fixed-point loop. The global prior is ANCHORED (every input fit once before the
    solve), so the single FB pass is exact.

    BEFORE the pass: the non-circular population gDNA prior on gDNA-clean seeds (:func:`_gdna_seed_estimate`) —
    exposure-pooled ``ρ_global`` + the seed gDNA spread ``σ²_g`` (the GLOBAL prior's between-node spread, NOT a
    message reliability), the conservative intergenic+intron floor, and (pass 2) the trained Phase-2 gDNA-density
    KDE; all anchored (fit once before the solve).

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
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
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
    SP = (geometry.spliced_pos_left, geometry.spliced_pos_right)
    SN = (geometry.spliced_neg_left, geometry.spliced_neg_right)
    # per-node "global" gDNA support: region = its contained support; boundary = the AVERAGED per-side density
    # length ½(E_l+E_r) over the total crossing mass.
    eff_global = np.where(is_reg, EG[0], 0.5 * (EG[0] + EG[1]))
    mass_global = np.where(is_reg, MS[0], MS[0] + MS[1])

    # The per-node solve is the log-density 1-D/2-D log-odds solver (simplex_logodds, O(m·K),
    # genome-scale-tractable). The "solve grid" is the f_g axis the global NB prior is evaluated on (the
    # log-odds σ(λ) lattice).
    _lam_lo, _fg_lo = _logodds_grid(int(n_grid), float(logodds_window))
    solve_grid = _fg_lo
    kappa = float(rna_sense_frac)
    od_g, od_r = gdna_strand_overdispersion, rna_strand_overdispersion
    _zero_spl = np.zeros_like(np.asarray(statics.mass_unspliced, dtype=np.float64))

    def _local_solve(g_arr, gm=None, gp=None, rm=None, rp=None):
        """The per-node local/final solve (log-density log-odds backend). Returns
        ``(f_g median, f_pos mean, f_neg mean, var_g, var_pos, var_neg)``. Phase A calls it message-free;
        phase D passes the FB messages (fraction-space)."""
        dc = _solve_nodes_logodds_all(
            statics.u_pos,
            statics.u_neg,
            fp,
            fn,
            statics.strand_obs,
            statics.mass_unspliced,
            _zero_spl,
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
        )
        return (
            dc.gdna_frac,
            dc.rna_pos_frac,
            dc.rna_neg_frac,
            dc.gdna_frac_var,
            dc.rna_pos_frac_var,
            dc.rna_neg_frac_var,
        )

    # Two distinct gates, both from the region SIGNATURE (never the counts — count-zero-info):
    #   * SOLVE gate (`solvable`): a node deconvolves its own gDNA-vs-RNA split iff it admits ≥1 RNA strand and
    #     has unspliced mass. A G1 node — no admissible RNA strand: an intergenic region, or a gene-boundary seam
    #     (TSS/TES, opposite-strand exon↔exon) — is a LOCKED all-gDNA node ({0,0,1}); it is not solved, keeping
    #     its init (RNA cannot cross a gene boundary, so its unspliced mass is purely gDNA).
    #   * EMISSION gate (per component, in `_scan`): which MESSAGES a node sends. A three-term Boolean over the
    #     components gDNA / +RNA / −RNA, structural and symmetric — defined at the top of the scan loop.
    solvable = (fp | fn) & (statics.mass_unspliced > 0.0)

    # The NON-CIRCULAR population gDNA prior on gDNA-clean seeds (§4.3), ANCHORED — computed ONCE, never refit:
    # the exposure-pooled rate ρ_global (the global baseline) + var_mean (its variance) + the seed gDNA spread
    # σ²_g (``gdna_vm`` — the GLOBAL prior's between-node spread, NOT a message reliability). Strand seeds use the
    # strand-ONLY init f_g. (Message reliability is the belief-free σ²_imp in `_scan` — no var~mean message fit.)
    rho_global, gdna_vm, var_mean = _gdna_seed_estimate(
        chain, statics, geometry, region_arrays, boundary_substrate, f_g, kappa
    )
    # The DEPLETED gDNA FLOOR from intergenic + intron REGIONS (the user's empirical-prior directive): a
    # confident floor that pins depleted intron nodes (the nascent-hallucination fix), letting an intron
    # with density excess over the floor carry nascent. Exons/boundaries keep the genome-wide global.
    rho_floor, s2_floor, var_mean_floor, floor_mask, s2_bg = _floor_estimate(
        chain, geometry, region_arrays, f_g, kappa
    )
    # The global gDNA prior on f_g, an (n_nodes, K) log-term on the solve grid, applied to ALL nodes. It is
    # ALWAYS the weak, exposure-pooled stability floor (`_global_logprior`: M-INDEPENDENT, capped at one
    # pseudo-observation, so it never overrules a node's own strand evidence; the depleted-floor override).
    # PASS 2 ADDS the trained Phase-2 mixture ON TOP (`_kde_logprior`). This layering is the design's
    # "floor + mixture" (§3.4): the pooled floor is the always-present "gDNA is scarce" baseline — its
    # downward pull is what resolves a balanced AMBIG node (whose two-strand RNA is strand-degenerate with
    # gDNA) toward ~0 in a gDNA-poor context — while the KDE supplies the enriched-mode structure the floor
    # lacks. Replacing the floor with the KDE (an earlier attempt) lost the downward anchor, because the
    # per-node KDE's RNA-residual mode sits ABOVE the exposure-pooled ρ_global (dissection: gDNA=0 AMBIG
    # went 0.02→0.38). ANCHORED — every input is fit once, so the prior is CONSTANT within a pass.
    # PASS 1 (gdna_prior=None): the weak stability floor only (strand + messages carry the single-strand solve;
    # the floor keeps low/zero-gDNA nodes finite and anchors the depleted population that trains the KDE).
    # The floor nodes (intergenic + introns) get the DATA-DRIVEN intron density likelihood (`s2_bg`, uncapped
    # in `_global_logprior`) — the SOLE bleed-stopper now that the mature message split is retired.
    _s2_bg = s2_bg
    global_lp = _global_logprior(
        solve_grid,
        mass_global,
        eff_global,
        rho_global,
        gdna_vm,
        var_mean,
        floor_mask=floor_mask,
        rho_floor=rho_floor,
        s2_floor_total=var_mean_floor + s2_floor,
        s2_bg=_s2_bg,
    )
    # PASS 2 (gdna_prior set): ADD the generative two-density prior — the empirical gDNA-density KDE (real
    # tails) × the Jeffreys RNA prior 1/(1−f_g) (`_kde_logprior`). This is the density-prior INTEGRATION: the
    # strand landscape (ψ, added in the solve) × the population landscape. The Jeffreys term removes the
    # gDNA-priored/RNA-free asymmetry that caused the boundary cliff + the false-positive gDNA; the real KDE
    # tails stop a high-density node drifting to f_g≈0.5. Derivation + validation:
    # ambig_boundary_spliced_deconvolution.md. Applied to ALL solvable nodes (self-scaling: a confident strand
    # dominates ψ, an AMBIG/thin node leans on the population).
    if gdna_prior is not None:
        kde_lp = _kde_logprior(solve_grid, mass_global, eff_global, gdna_prior)
        if _s2_bg is not None:
            # PHASE 2c — node-type-specific density prior (no double-count). The FLOOR nodes (intergenic +
            # introns, off-target) use the intergenic-background density LIKELIHOOD (Gaussian + Jeffreys, built
            # in `_global_logprior` above), so they are EXCLUDED from the capture-aware KDE here — the KDE is
            # the density model for on-target ENRICHED exons/boundaries. This gives exactly ONE density prior
            # and ONE Jeffreys RNA-parsimony per node; without the exclusion, floor nodes would carry the
            # Jeffreys twice (once from each), doubling the parsimony and over-pinning introns to gDNA.
            kde_lp = kde_lp.copy()
            kde_lp[np.asarray(floor_mask, bool), :] = 0.0
        global_lp = global_lp + kde_lp

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

    # (A) LOCAL message-free beliefs (backend-dispatched).
    fg_loc, fp_loc, fn_loc, vg_loc, vp_loc, vn_loc = _local_solve(global_lp)
    pg_loc = 1.0 / np.maximum(
        vg_loc, _EPS
    )  # local precision (var floored: a sharp belief ⇒ large finite)
    pp_loc = 1.0 / np.maximum(vp_loc, _EPS)
    pn_loc = 1.0 / np.maximum(vn_loc, _EPS)

    # G1-EMISSION FIX (`emit_locked`): a G1-locked seam (intergenic / TSS / TES / opposite-strand exon↔exon)
    # is neither single-strand nor AMBIG, so `_local_solve` SKIPS it and returns f_g=0. The sweep seeds each
    # node's running EMISSION belief (`fbg` in `_scan`) from that local solve, so a locked all-gDNA seam emits
    # n_src = f_g·mass = 0 — it is SILENT, contradicting this module's "a locked all-gDNA node is a confident
    # emitter" (a high-count gene-boundary gDNA crossing sends nothing to the adjacent exon). Reseed the locked
    # nodes' running belief from their LOCKED init belief so they emit their structural gDNA. Their local
    # precision is already 1/ε (vg_loc=0 for skipped nodes) so they stay locked and are not swayed by incoming
    # messages. Their own FINAL value is unchanged — the solvable write-back keeps the init for non-solvable
    # nodes — this only restores their OUTGOING messages. Empty nodes (mass 0) are also reseeded but never emit
    # (the sm>0 gate). Default off (bit-identical) pending all-24 validation, incl. gdna_none false-positive.
    if emit_locked:
        locked = ~np.asarray(solvable, bool)
        fg_loc = np.where(locked, f_g, fg_loc)
        fp_loc = np.where(locked, f_pos, fp_loc)
        fn_loc = np.where(locked, f_neg, fn_loc)

    # Belief-free Poisson message precision (`disagreement_shrinkage_prior_design_v2.md`): σ²_msg = σ²_imp +
    # 1/n_src ⇒ pr = n_src/(n_src·σ²_imp + 1). ONE total-density scalar for every channel (gDNA + both RNA
    # strands); σ²_imp is the empirical adjacent-node imputation floor (`adjacent_disagreement_variance`).
    sig_imp = float(disagreement_sigma2)

    def _scan(seq, nbr, sf, df):
        """Sequential scan: project the running belief from each node's ``nbr`` (src face ``sf`` → dst face
        ``df``) into the dst, then fold it into the dst's running belief (local ⊗ incoming message; NOT the
        reverse message — true tree BP). Returns the per-node message (mode, prec) per component. O(m).

        BELIEF-FREE Poisson disagreement-variance precision (``disagreement_shrinkage_prior_design_v2.md``):
        ``σ²_msg = σ²_imp + 1/n_src`` ⇒ ``pr = n_src/(n_src·σ²_imp + 1)`` — denom ≥ 1, ``pr=0`` exactly at
        ``n_src=0`` (no message), no clamp. ``σ²_imp`` (``sig_imp``) is the empirical adjacent-node total-density
        imputation floor; a confident source cannot shrink it to 0, so no message is ever overconfident. Same
        scalar for all three components (gDNA / ±RNA). Replaces the retired ``σ²_bio(μ)`` var~mean curve."""
        fbg, fbp, fbn = (
            fg_loc.copy(),
            fp_loc.copy(),
            fn_loc.copy(),
        )  # running belief (starts at local)
        vbg, vbp, vbn = vg_loc.copy(), vp_loc.copy(), vn_loc.copy()
        amg, apg = np.zeros(n_nodes), np.zeros(n_nodes)  # gDNA message (mode, prec)
        amp, app = np.zeros(n_nodes), np.zeros(n_nodes)  # RNA-pos
        amn, apn = np.zeros(n_nodes), np.zeros(n_nodes)  # RNA-neg
        EGs, EGd, ERs, ERd = EG[sf], EG[df], ER[sf], ER[df]
        MSs, MSd, SPs, SNs = MS[sf], MS[df], SP[sf], SN[sf]
        ESPs = ESP[sf]  # source-face spliced eff-len (for the mature-RNA MEASUREMENT message)
        SPd, SNd, ESPd = (
            SP[df],
            SN[df],
            ESP[df],
        )  # DEST-face spliced — the mature ABSORBED at a junction
        # (subtracted from an exon→boundary message so only NASCENT crosses into the intron side).
        # The running belief combines in LOG-fraction space (the message is a Gaussian on log f_c).
        # Precompute the local log-fractions (constant across the scan) for the combine.
        lfg_loc = np.log(np.maximum(fg_loc, _EPS))
        lfp_loc = np.log(np.maximum(fp_loc, _EPS))
        lfn_loc = np.log(np.maximum(fn_loc, _EPS))
        for i in seq:
            lsrc = nbr[i]
            if lsrc < 0:
                continue
            md = MSd[i] if MSd[i] > _EPS else _EPS
            egd = EGd[i] if EGd[i] > _EPS else _EPS
            erd = ERd[i] if ERd[i] > _EPS else _EPS
            sm = MSs[lsrc]  # source facing unspliced mass
            # STRUCTURAL emission gate — one Boolean per component; src AND dst must admit it. gDNA is
            # genomically universal (admitted everywhere) ⇒ it gates on facing mass alone; each RNA strand
            # transmits only where THAT strand is continuous across the edge (free_pos / free_neg on both).
            emit_g = sm > _EPS
            emit_p = fp[lsrc] and fp[i] and (sm > _EPS or SPs[lsrc] > _EPS)
            emit_n = fn[lsrc] and fn[i] and (sm > _EPS or SNs[lsrc] > _EPS)
            # gDNA — a G1 seam (intergenic / TSS / TES) is a locked, confident all-gDNA emitter.
            if emit_g:
                eg = EGs[lsrc] if EGs[lsrc] > _EPS else _EPS
                n_src = fbg[lsrc] * sm  # source gDNA COUNT (deconvolved)
                rho = n_src / eg  # source gDNA DENSITY ρ_g_src — the MESSAGE currency
                # DENSITY message (NO fractions in the wire): the content is the source gDNA density
                # ρ_g_src; the RECEIVER re-expresses it in its OWN log-f_g solve frame via its gDNA
                # density base M_dst/E_gdna_dst (= md/egd), flooring ρ at the dst min-observable density
                # 1/egd. (Fractions are not comparable across nodes — only densities are.)
                mo = math.log(max(rho, 1.0 / egd) / (md / egd))
                # Poisson disagreement-variance: σ²_msg = σ²_imp + 1/n_src ⇒ pr = n_src/(n_src·σ²_imp + 1) —
                # denom ≥ 1, pr=0 exactly at n_src=0 (no message; "zero density is not a measurement"), no clamp.
                pr = n_src / (n_src * sig_imp + 1.0)
                amg[i], apg[i] = mo, pr
                pt = pg_loc[i] + pr
                fbg[i] = math.exp((pg_loc[i] * lfg_loc[i] + pr * mo) / pt)
                vbg[i] = 1.0 / pt
            # RNA messages (± symmetric) — the imputed density is the total unspliced RNA (nascent + exon-body
            # mature) crossing contiguously, plus the junction spliced source/absorption. MATURE crosses only at
            # a junction boundary via the one-sided spliced terms (the geometry gates it — `spliced_efflen…`):
            #   ρ = src unspliced RNA ``fbp·sm/E_r``
            #     + src-face mature ``SPs/E_spl``  — ADDED at a B→exon source (MEASURES the exon's mature)
            #     − dst-face mature ``SPd/E_spl``  — ABSORBED at an exon→B dest (spliced sink; leaves nascent).
            # At most one of SPs/SPd is nonzero per edge, and both are 0 on intron↔boundary edges, so introns
            # carry pure nascent with no gate. Precision: the belief-free σ²_imp + 1/n_src (`sig_imp`).
            if emit_p:
                # LUMPED RNA⁺ message (the 2-component {RNA,gDNA} model): total unspliced RNA + the junction
                # spliced MEASUREMENT/absorption, NO mature-vs-nascent decomposition. RNA flows into introns; the
                # intron-density likelihood (`_global_logprior`, `_s2_bg`) is the sole bleed-stopper.
                er = ERs[lsrc] if ERs[lsrc] > _EPS else _EPS
                esp = ESPs[lsrc] if ESPs[lsrc] > _EPS else _EPS
                n_nasc = (
                    fbp[lsrc] * sm
                )  # source total unspliced RNA count (nascent + exon-body mature)
                n_mat = SPs[
                    lsrc
                ]  # source-face spliced (>0 only B→exon): MEASURES the exon's mature
                rho_mat_dst = SPd[i] / (
                    ESPd[i] if ESPd[i] > _EPS else _EPS
                )  # dst-face mature absorbed (exon→B)
                rho = (
                    n_nasc / er + n_mat / esp - rho_mat_dst
                )  # total-RNA density (+ MEASUREMENT into an exon)
                mo = math.log(
                    max(rho, 1.0 / erd) / (md / erd)
                )  # → dst log-f_pos frame (floored at min-observable)
                n_src = n_nasc + n_mat  # source RNA⁺ count (Poisson sampling)
                pr = n_src / (n_src * sig_imp + 1.0)
                amp[i], app[i] = mo, pr
                pt = pp_loc[i] + pr
                fbp[i] = math.exp((pp_loc[i] * lfp_loc[i] + pr * mo) / pt)
                vbp[i] = 1.0 / pt
            # RNA-neg — symmetric (mature on the −strand junction motif; same 3-term nascent message).
            if emit_n:
                # LUMPED RNA⁻ message — mirror of RNA⁺ (mature on the −strand junction motif; spliced source SNs).
                er = ERs[lsrc] if ERs[lsrc] > _EPS else _EPS
                esp = ESPs[lsrc] if ESPs[lsrc] > _EPS else _EPS
                n_nasc = fbn[lsrc] * sm
                n_mat = SNs[lsrc]
                rho_mat_dst = SNd[i] / (
                    ESPd[i] if ESPd[i] > _EPS else _EPS
                )  # dst-face mature absorbed (exon→B)
                rho = n_nasc / er + n_mat / esp - rho_mat_dst
                mo = math.log(max(rho, 1.0 / erd) / (md / erd))  # → dst log-f_neg frame
                n_src = n_nasc + n_mat  # source RNA⁻ count (Poisson sampling)
                pr = n_src / (n_src * sig_imp + 1.0)
                amn[i], apn[i] = mo, pr
                pt = pn_loc[i] + pr
                fbn[i] = math.exp((pn_loc[i] * lfn_loc[i] + pr * mo) / pt)
                vbn[i] = 1.0 / pt
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
    # (D) FINAL solve with the FB messages (backend-dispatched).
    mg_, mp_, mn_, vg_, vp_, vn_ = _local_solve(
        global_lp, mode_g, prec_g, (mode_p, mode_n), (prec_p, prec_n)
    )
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
            statics.strand_obs,
            statics.mass_unspliced,
            _zero_spl,
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
            rho_global=rho_global,
            gdna_vm=gdna_vm,
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
            rho_floor=rho_floor,
            floor_mask=floor_mask,
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
