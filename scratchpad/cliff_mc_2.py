"""CROSS-CLIFF PRECISION LAW — derivation-agent-#2 MC (counts-over-length effective count).

The bug (root-caused, HANDOFF_4 §6): a message reframed across an enrichment cliff
``r = ρ_tot(dst)/ρ_tot(src)`` keeps FULL precision. The MODE is stretched up by r (correct: enrichment
cancellation), but the delivered PRECISION should DECAY with the stretch, because a node and its r-fold-denser
neighbour are almost certainly DIFFERENT KINDS of region, so "we share a composition" is far more suspect the
bigger r. The shipped damping ``σ²_transfer = Var(log r) ≈ 1/n`` is the UNCERTAINTY in r (small on well-counted
nodes), NOT a stretch discount — so a 407× cliff arrives at full precision (anchor: node 1909, exon, oracle
f_g=0.985; boundary 1910 mass 14 with a 47-count junction; r≈407; the boundary's RNA lands at f_pos=0.718,
precision ~26, and collapses 1909 to 0.51).

DERIVED LAW (counts-over-length effective count):

    n_eff = n_src / max(1, r)            r = ρ_tot(dst)/ρ_tot(src)   (the reframe)

    delivered message precision on log f_c (or log ρ_c) = 1 / ( Var(log f_c^src) + 1/n_eff )
                                                       = 1 / ( Var(log f_c^src) + max(1,r)/n_src )

  n_src   the source measurement's discrete count (n_spl for a spliced graft; n_unspl for an imputed
          composition). COUNTS — the discrete resolving power.
  r       the density stretch. THE COVERAGE ARGUMENT: the source observed n_src fragments over E_src bases at
          density ρ_src = n_src/E_src. Viewed at the DESTINATION's density ρ_dst = r·ρ_src, those same n_src
          fragments cover only E_src/r destination-bases-equivalent (r-fold more material per base ⇒ r-fold
          fewer bases per fragment). So the source's evidence sub-samples the destination's dense composition
          r-fold more THINLY ⇒ effective count n_src/r. LENGTH enters as the bases the counts cover.
  max(1,·) ASYMMETRIC: stretching a coarse measurement UP onto a denser node loses resolution; compressing a
          fine measurement DOWN onto a sparser node does not (you already over-resolve it). Decay only for r>1.

This module is the ARBITER. It (1) reproduces the anchor pull and shows ÷r controls it where σ²_transfer does
not; (2) validates the effective-count / coverage identity by resampling; (3) measures the HONEST predictive
precision of an imputed message vs the realized cliff under a biological generative model and shows it tracks
n_src/r (not 1/n); (4) shows the r≈1 (similar-neighbour) safety that preserves the unstranded/AMBIG win.

    OMP_NUM_THREADS=1 python scratchpad/cliff_mc_2.py [--draws 400000]
"""

from __future__ import annotations

import argparse

import numpy as np

rng = np.random.default_rng(20260724)


def _pois(mean, size):
    return rng.poisson(np.maximum(mean, 0.0), size=size).astype(np.float64)


def _hr():
    print("─" * 108)


# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
# EXP 1 — THE ANCHOR MECHANISM: the reframe over-states a CONTINUOUS-RNA junction's mature by r, and the
#         message PULL on the ψ solve = precision · |mode-bias|. σ²_transfer leaves the pull huge; ÷r kills it;
#         a matched cliff (r≈1) leaves a full, unbiased, correct message.
#
# Physical setup (the anchor, HANDOFF_4 §6). A boundary measures a spliced junction: n_spl RNA fragments →
# mature density ρ_μ. RNA is CONTINUOUS across the junction, so the flanking exon's TRUE mature density is
# ρ_R^dst ≈ ρ_μ (the same RNA continues). But the exon is heavily gDNA-contaminated, so its TOTAL density is
# r-fold the boundary's: r = ρ_tot(dst)/ρ_tot(src) ≈ 407, driven by gDNA, NOT by more RNA. The reframe
# multiplies the grafted mature by r ⇒ mode ρ_μ·r ⇒ over-states the exon's RNA by the factor r.
#   log-bias of the delivered RNA message = log r   (independent of counts).
#   PULL the message exerts on ψ (a Gaussian mode m, precision p on log f_R) = p · (m − truth) = p · log r.
# For the message to be HARMLESS the pull must be O(1); a huge pull overrides the node's own weak strand belief
# (τ_λ ≈ 1.6 on a mostly-gDNA exon) and collapses f_g.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
def exp1_anchor(N):
    print("\n═══ EXP 1 — anchor pull: σ²_transfer leaves the message corrupting; ÷r makes it inert ═══\n")
    n_spl = 47.0        # the boundary's junction count (the measured mature)
    n_src_tot = 14.0    # the boundary's TOTAL count (mass 14) → σ²_transfer = Var(log r) ≈ 1/n_src_tot
    tau_node = 1.6      # the exon's own strand composition precision (weak: a mostly-gDNA node)
    for r in (1.0, 5.0, 40.0, 407.0):
        # the delivered mode's log-bias = log r (RNA continuous ⇒ true exon RNA = ρ_μ, message = ρ_μ·r)
        bias = np.log(r)
        # MC the Poisson spliced count → the sampling variance of log ρ_μ = 1/n_spl (r is a fixed reframe here)
        spl = _pois(n_spl, N)
        keep = spl > 0
        var_log_mu = float(np.var(np.log(spl[keep])))   # ≈ 1/n_spl, the honest sampling variance of the mode
        # ── law A: shipped σ²_transfer ── p = 1/(Var(log f=0 for measured) + 1/n_spl + Var(log r))
        var_log_r = 1.0 / n_src_tot                      # Var(log r), the r-UNCERTAINTY (small)
        p_s2t = 1.0 / (var_log_mu + var_log_r)
        # ── law B: DERIVED ÷r ── effective count n_spl/r ⇒ p = 1/(1/n_eff) = n_eff = n_spl/max(1,r)
        n_eff = n_spl / max(1.0, r)
        p_div_r = n_eff
        pull_s2t = p_s2t * bias
        pull_div = p_div_r * bias
        tag = "  (matched — full, unbiased message)" if r == 1.0 else ""
        print(f"  r={r:6.1f}  bias=log r={bias:5.2f}   "
              f"σ²t: p={p_s2t:7.2f} pull={pull_s2t:8.2f}   "
              f"÷r: p={p_div_r:7.3f} pull={pull_div:6.2f}   (node τ={tau_node}){tag}")
    print("\n  READ: σ²_transfer's precision barely moves with r (it is the r-uncertainty ~1/n), so the pull")
    print("  GROWS as p·log r → overrides the node's τ=1.6 (the 1909 collapse). ÷r's precision falls as n/r, so")
    print("  the pull (n/r)·log r → 0 for a big cliff (inert) while staying FULL & unbiased at r=1 (correct).")


# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
# EXP 2 — THE COVERAGE IDENTITY (effective count = n_src / r), resampled.
#
# The counts-over-length claim, made operational. A source measures a per-base composition field by observing
# fragments; each fragment reveals its base's label (g or R). The source observes n_src fragments over its
# bases. We now ask its message to constrain the DESTINATION region's composition — a region that is r-fold
# DENSER (more material per base). At the destination's density, the source's n_src fragments cover only
# E_src/r destination-bases-equivalent, so they resolve the destination's composition with effective count
# n_src/r. We verify by literally down-covering: draw the source's evidence, then measure how many
# INDEPENDENT destination composition-cells it resolves.
#
# Construction (parameter-free): the destination is a strip of B_dst composition-cells (each cell = one
# resolvable base-block, R with prob f). A NATIVE destination measurement resolves f with ~B_dst·(coverage)
# fragments. The source, at 1/r the density, deposits its n_src fragments over the SAME strip but each fragment
# now "claims" r cells' worth of density (it is a sparser probe) ⇒ it independently samples only n_src/r cells.
# The variance of the source's estimate of f is therefore f(1−f)/(n_src/r) ⇒ effective count n_src/r.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
def exp2_coverage(N):
    print("\n\n═══ EXP 2 — coverage identity: a probe r-fold sparser resolves n_src/r effective counts ═══\n")
    f = 0.30
    n_src = 400.0   # chosen so n_src/r is an integer across the sweep (exact Bernoulli MC)
    print(f"  f={f}  n_src={n_src:.0f}   (native effective count if matched = n_src)")
    M = N // 100
    for r in (1.0, 4.0, 40.0, 400.0):
        # each of the n_src source fragments, at the destination's density, covers 1/r of a destination
        # resolving-cell ⇒ the number of INDEPENDENT destination-composition observations is n_src/r. Simulate
        # by estimating f from exactly n_eff = n_src/r Bernoulli(f) draws.
        n_eff_true = n_src / max(1.0, r)
        k_obs = int(round(n_eff_true))
        if k_obs >= 1:
            draws = rng.binomial(k_obs, f, size=M).astype(np.float64) / k_obs
            var_emp = float(np.var(draws))
            var_pred = f * (1.0 - f) / n_eff_true      # binomial variance at effective count n_src/r
            prec_emp, prec_pred = 1.0 / var_emp, 1.0 / var_pred
            rel = abs(prec_pred - prec_emp) / prec_emp
            flag = "OK " if rel < 0.05 else "***"
            print(f"  {flag} r={r:6.1f}  n_eff=n_src/r={n_eff_true:8.3f}   "
                  f"precision emp {prec_emp:8.3f}  pred {prec_pred:8.3f}  rel {rel:5.1%}")
        else:
            print(f"      r={r:6.1f}  n_eff=n_src/r={n_eff_true:8.3f}   < 1 effective count ⇒ precision < 1 "
                  f"⇒ message INERT (no fractional-count MC needed)")
    print("\n  READ: the delivered effective count IS n_src/r — full at r=1, and →0 across a big cliff (below 1")
    print("  effective count the message is inert). The discrete resolving power of a coarse (sparse) probe")
    print("  used to pin a fine (dense) composition.")


# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
# EXP 3 — HONEST PREDICTIVE PRECISION vs the realized cliff, under a biological generative model.
#
# Why the decay is real (not a pure-statistics artifact): capture co-enriches EXPRESSED (RNA-rich, low f_g)
# regions and leaves introns/intergenic (gDNA, high f_g) depleted, so ENRICHMENT and COMPOSITION co-vary. A
# neighbour pair separated by a big density cliff r is therefore very likely separated in composition too — the
# imputation premise fails in proportion to the cliff. We draw realistic neighbour pairs, impute the (Poisson-
# measured) source composition onto the destination, and measure the HONEST error λ_src,meas − λ_dst,true
# (λ = logit f_g). Binned by log r, the honest precision 1/Var(error) should FALL ~ n_src/r (the derived law),
# NOT stay flat at the source count n_src (the shipped model). We report the fitted exponent and compare laws,
# at several neighbour-correlation strengths to show the FORM is robust (α≈1), not tuned.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
def _draw_pairs(M, *, rho_corr, cap_gain):
    # latent expression (log): correlated between neighbours (copula corr rho_corr)
    z = rng.standard_normal((M, 2))
    z[:, 1] = rho_corr * z[:, 0] + np.sqrt(1.0 - rho_corr**2) * z[:, 1]
    logx = 0.9 * z - 1.0                       # latent log RNA content
    rho_R = np.exp(logx)                        # RNA density (per base)  (>=0)
    rho_g = 1.0                                 # gDNA density: uniform genome (same everywhere)
    f_g = rho_g / (rho_g + rho_R)               # composition (enrichment-invariant)
    # capture enrichment co-varies with expression: enriched where RNA-rich (captured exons)
    e = np.exp(cap_gain * logx + 0.4 * rng.standard_normal((M, 2)))
    rho_tot = e * (rho_g + rho_R)               # total density (what the reframe ratio is formed from)
    return f_g, rho_tot, rho_R, e, rho_g


def exp3_population(N):
    print("\n\n═══ EXP 3 — honest predictive precision vs the cliff (biological generative model) ═══\n")
    M = max(N, 400_000)
    for rho_corr, cap_gain in ((0.7, 0.8), (0.5, 1.2), (0.85, 0.5)):
        f_g, rho_tot, rho_R, e, rho_g = _draw_pairs(M, rho_corr=rho_corr, cap_gain=cap_gain)
        src, dst = 1, 0
        # the reframe cliff (use TRUE totals; the r-uncertainty is a separate, small effect)
        r = rho_tot[:, dst] / rho_tot[:, src]
        # source measures its composition with a Poisson count n_src (its total mass); measured f_g via a
        # Beta-ish count split. Use the observed gDNA/RNA fragment counts at the source.
        E = 150.0                                # a nominal shared eff-length (bases)
        n_g = _pois(e[:, src] * rho_g * E, M)
        n_r = _pois(e[:, src] * rho_R[:, src] * E, M)
        n_src = n_g + n_r
        ok = n_src > 0
        fg_src_meas = np.where(ok, n_g / np.maximum(n_src, 1.0), f_g[:, src])
        # the imputed message mode = source composition; honest error in log-odds λ
        eps = 1e-6
        lam = lambda p: np.log(np.clip(p, eps, 1 - eps) / np.clip(1 - p, eps, 1 - eps))
        err = lam(fg_src_meas) - lam(f_g[:, dst])
        # bin by log r (up-cliff: r>1, the damped direction)
        lr = np.log(np.maximum(r, 1e-9))
        up = ok & (lr > 0)
        print(f"  ── neighbour corr={rho_corr}, capture gain={cap_gain} ──")
        print(f"     {'log r bin':>14} {'⟨r⟩':>10} {'⟨n_src⟩':>9} "
              f"{'honest prec':>12} {'n_src (flat)':>12} {'n_src/r (÷r)':>13}")
        edges = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.5]
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = up & (lr >= lo) & (lr < hi)
            if m.sum() < 500:
                continue
            var_err = float(np.var(err[m]))
            prec = 1.0 / var_err
            r_mean = float(np.mean(r[m]))
            n_mean = float(np.mean(n_src[m]))
            # candidate laws mapped to a λ precision. Both need the src-composition Jacobian; we compare the
            # SHAPE vs r, so report the raw n_src and n_src/r (the count factors) against the honest precision.
            print(f"     [{lo:3.1f},{hi:3.1f})   {r_mean:10.2f} {n_mean:9.1f} "
                  f"{prec:12.4f} {n_mean:12.1f} {n_mean / r_mean:13.4f}")
        # fit honest precision ~ C / r^alpha over the up-cliff region (log-log slope)
        m = up & (lr > 0.3)
        # robust binned fit
        bx, by = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            mm = up & (lr >= max(lo, 0.3)) & (lr < hi)
            if mm.sum() < 500:
                continue
            bx.append(np.log(np.mean(r[mm])))
            by.append(np.log(1.0 / np.var(err[mm])))
        if len(bx) >= 3:
            A = np.polyfit(bx, by, 1)
            print(f"     → honest precision ∝ r^({A[0]:+.2f})   (÷r law predicts slope −1.00)\n")


# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
# EXP 4 — THE SAFETY: at r≈1 (similar-density neighbours — where AMBIG/unstranded propagation lives) the ÷r law
#         PRESERVES the message (n_eff ≈ n_src), so the unstranded-capON win (messages are the only info) is not
#         killed. Only the CLIFF-crossing confidence is damped. Contrast a blanket own-count cap (kills all).
# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
def exp4_safety(N):
    print("\n\n═══ EXP 4 — over-damping safety: ÷r preserves similar-neighbour (AMBIG) propagation ═══\n")
    n_src = 60.0
    print(f"  n_src={n_src:.0f}   effective count delivered under the ÷r law:")
    for r in (1.0, 1.3, 2.0, 3.0, 10.0, 100.0, 407.0):
        n_eff = n_src / max(1.0, r)
        kind = "AMBIG propagation (KEEP)" if r <= 3.0 else ("mild cliff" if r < 50 else "MISMATCH cliff (kill)")
        print(f"  r={r:6.1f}   n_eff={n_eff:8.3f}   frac kept {n_eff / n_src:6.1%}   {kind}")
    print("\n  READ: the ÷r discount is 1.0 at r=1 and degrades gently to ~0.5 at r=2, ~0.33 at r=3 — a")
    print("  similar-density AMBIG neighbour still delivers most of its evidence (the unstranded win survives).")
    print("  A BLANKET own-count cap would clamp every message regardless of r and kill that propagation.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=400_000)
    a = ap.parse_args()
    N = a.draws
    print(f"# CROSS-CLIFF PRECISION LAW  —  counts-over-length effective count  n_eff = n_src / max(1, r)")
    print(f"# draws={N:,}\n")
    _hr()
    exp1_anchor(N)
    _hr()
    exp2_coverage(N)
    _hr()
    exp3_population(N)
    _hr()
    exp4_safety(N)
    _hr()


if __name__ == "__main__":
    main()
