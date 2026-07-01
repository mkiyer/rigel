"""Study the Phase-2 KDE under REALISTIC zero-inflation + class imbalance.

Real transcriptomes annotate >500k transcripts, of which a small minority are expressed in any sample →
a vast majority of exons/introns/boundaries carry ZERO counts. This harness builds a zero-inflated toy
(many genes, few expressed) and dissects:
  * the substrate composition (zero vs nonzero, by node kind) — the imbalance;
  * where the zero-count nodes land (their per-node 1/E density is E-dependent → do the zeros SMEAR?);
  * the pooled floor ρ_floor = (a+ΣG)/(ΣE) vs the per-node zeros;
  * the KDE under different weightings: PRECISION (current) vs UNIT vs CLASS-BALANCED.

Usage: python -m scripts.debug.zero_inflation_study [--n-genes N] [--frac-expr F] [--gdna G]
"""
from __future__ import annotations

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from scripts.debug.toy_prod import run, SCRATCH  # noqa: E402
from rigel.calibration.gdna_density_prior import (  # noqa: E402
    build_training_substrate, GdnaDensityPrior, TrainingSubstrate, KIND_NAMES,
)

_EPS = 1e-9


def _genes(n_genes, frac_expr, seed=13):
    """n_genes single-strand genes; the first ceil(frac·n) are expressed, the rest unexpressed (abundance 0).
    Varied exon sizes so E-lengths span a realistic range (the zero-count 1/E smear test)."""
    rng = np.random.RandomState(seed)  # noqa: NPY002  (deterministic layout only, not in the solve path)
    genes, x, INTRON, GAP = [], 5000, 2000, 4000
    n_expr = max(1, int(np.ceil(frac_expr * n_genes)))
    sizes_pool = [150, 300, 600, 1200, 2500, 5000]
    for k in range(n_genes):
        s = x
        ex = []
        for _ in range(3):
            sz = int(sizes_pool[rng.randint(len(sizes_pool))])
            ex.append((s, s + sz)); s = s + sz + INTRON
        genes.append((f"G{k}", "+" if k % 2 == 0 else "-", ex, 100.0 if k < n_expr else 0.0))
        x = s + GAP
    return genes, x + GAP, n_expr


def _refit(sub, weight, bandwidth):
    """Refit the KDE with a replaced weight vector (same densities/std)."""
    s2 = TrainingSubstrate(log_rho=sub.log_rho, weight=weight, node_kind=sub.node_kind,
                           node_index=sub.node_index, log_rho_std=sub.log_rho_std)
    return GdnaDensityPrior.fit(s2, bandwidth=bandwidth)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-genes", type=int, default=80)
    ap.add_argument("--frac-expr", type=float, default=0.1)
    ap.add_argument("--gdna", type=float, default=0.3)
    ap.add_argument("--capture", action="store_true", default=True)
    ap.add_argument("--seed", type=int, default=13)
    args = ap.parse_args()

    genes, glen, n_expr = _genes(args.n_genes, args.frac_expr, args.seed)
    dbg: dict = {}
    _rdf, _solved, truth, tg, tr, sig = run(
        "zeroinfl", genes, kappa=0.99, n_rna=8000 * n_expr, gdna_fraction=args.gdna,
        capture=True, capture_strength=20.0, genome_length=glen, seed=args.seed, _debug=dbg)

    sub = build_training_substrate(
        dbg["chain"], dbg["belief"], dbg["geometry"], dbg["statics"],
        dbg["region_arrays"], dbg["boundary_substrate"], min_eff_length=dbg["gdna_fl_mean"])
    geom = dbg["geometry"]; belief = dbg["belief"]
    fg = np.asarray(belief.f_g); Ml = np.asarray(geom.mass_left)
    gcount_all = fg * Ml  # deconvolved gDNA count per node

    # per-training-node gDNA count (via node_index)
    gcount = gcount_all[sub.node_index]
    is_zero = gcount < 1.0

    print(f"\n=== zero-inflation study: {args.n_genes} genes, {n_expr} expressed "
          f"({100 * args.frac_expr:.0f}%), gdna={args.gdna} ===")
    print(f"substrate n={sub.n}  |  zero-count (gcount<1): {int(is_zero.sum())} "
          f"({100 * is_zero.mean():.0f}%)   nonzero: {int((~is_zero).sum())}")
    for code, name in KIND_NAMES.items():
        m = sub.node_kind == code
        if not m.any():
            continue
        z = m & is_zero
        print(f"  {name:>10}: n={int(m.sum()):4d}  zero={int(z.sum()):4d}  "
              f"nonzero={int((m & ~is_zero).sum()):3d}  "
              f"logρ[{sub.log_rho[m].min():+.1f},{sub.log_rho[m].max():+.1f}]  Σw={sub.weight[m].sum():.0f}")

    # the pooled floor vs the per-node zeros
    EGl = np.maximum(np.asarray(geom.eff_gdna_left), _EPS)[sub.node_index]
    G_zero = gcount[is_zero]; E_zero = EGl[is_zero]
    if is_zero.any():
        rho_floor = (1.0 + G_zero.sum()) / max(E_zero.sum(), _EPS)
        print(f"  POOLED floor over zero nodes: log ρ_floor={np.log(rho_floor):+.2f}  "
              f"(ΣG={G_zero.sum():.0f}, ΣE={E_zero.sum():.0f})")
        print(f"  per-node zero logρ SPREAD: [{sub.log_rho[is_zero].min():+.2f}, "
              f"{sub.log_rho[is_zero].max():+.2f}]  (E-dependent 1/E smear)")

    # KDE under different weightings (fixed bandwidth for comparability)
    h = 0.5
    pr_prec = _refit(sub, sub.weight, h)
    pr_unit = _refit(sub, np.ones_like(sub.weight), h)
    # class-balanced: each kind contributes equal total weight
    wcb = np.ones_like(sub.weight)
    for code in np.unique(sub.node_kind):
        m = sub.node_kind == code
        wcb[m] = 1.0 / max(int(m.sum()), 1)
    pr_cb = _refit(sub, wcb, h)
    print(f"  modes  precision={[round(m[0],1) for m in pr_prec.modes[:4]]}")
    print(f"         unit     ={[round(m[0],1) for m in pr_unit.modes[:4]]}")
    print(f"         classbal ={[round(m[0],1) for m in pr_cb.modes[:4]]}")

    # oracle
    reff = np.asarray(dbg["region_eff_len"]); orc = np.where((tg > 0) & (reff > 0), tg / np.maximum(reff, 1e-9), np.nan)
    orc_log = np.log(orc[np.isfinite(orc) & (orc > 0)])

    fig, ax = plt.subplots(figsize=(11, 6))
    if orc_log.size:
        ax.hist(orc_log, bins=40, density=True, alpha=0.15, color="gray", label="oracle region gDNA density")
    for pr, lab, c in ((pr_prec, "precision-weighted", "tab:red"),
                       (pr_unit, "unit weight", "tab:blue"),
                       (pr_cb, "class-balanced", "tab:green")):
        ax.plot(pr.x_grid, np.exp(pr.logP_grid), lw=2, color=c, label=f"{lab} (h={h})")
    ax.plot(sub.log_rho[is_zero], np.full(int(is_zero.sum()), -0.02), "|", color="navy", ms=8, alpha=0.4,
            label=f"zero-count nodes (n={int(is_zero.sum())})")
    ax.plot(sub.log_rho[~is_zero], np.full(int((~is_zero).sum()), -0.05), "|", color="darkred", ms=8, alpha=0.6,
            label=f"nonzero nodes (n={int((~is_zero).sum())})")
    ax.set_xlabel("log ρ_g"); ax.set_ylabel("density")
    ax.set_title(f"Zero-inflation: {args.n_genes} genes {int(100*args.frac_expr)}% expressed, gdna={args.gdna} "
                 f"— {int(100*is_zero.mean())}% zero-count nodes")
    smax = max(np.exp(pr_unit.logP_grid).max(), np.exp(pr_cb.logP_grid).max())
    ax.set_ylim(bottom=-0.08, top=2.0 * float(smax)); ax.legend(fontsize=8)
    fig.tight_layout()
    out = SCRATCH / f"zeroinfl_{args.n_genes}g_{int(100*args.frac_expr)}pct_g{args.gdna}.png"
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
