"""P2.2 — the Phase-2 gDNA-density-prior plotting framework.

Run a production-faithful toy through the real ``calibrate()`` (via ``toy_prod``, using its ``_debug``
hook to grab the solved chain internals), build the Phase-2 training substrate, fit ``P(log ρ_g)`` with
each bandwidth estimator, and PLOT the fit against the oracle. This is how we *see* the mixture prior — the
number/separation of modes vs bandwidth, whether the fitted modes land on the true gDNA-density levels, and
how the estimators behave — BEFORE any bandwidth is chosen for production or wired into the solve.

Usage:
    python -m scripts.debug.plot_gdna_prior                 # the default bimodal capture scenario
    python -m scripts.debug.plot_gdna_prior --kappa 0.5     # unstranded (structural training nodes only)
    python -m scripts.debug.plot_gdna_prior --no-capture    # unimodal (uniform gDNA)
Outputs a PNG under the scratchpad and prints the substrate + mode summary.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from scripts.debug.toy_prod import run, SCRATCH  # noqa: E402
from rigel.calibration.gdna_density_prior import (  # noqa: E402
    build_training_substrate, GdnaDensityPrior, KIND_NAMES,
)
from rigel.calibration.bp_solver import REGION  # noqa: E402

_KIND_COLOR = {0: "tab:blue", 1: "tab:green", 2: "tab:red", 3: "tab:orange"}


def _default_genes(n_genes: int):
    """n single-strand genes, alternating expressed / unexpressed, multi-exon (300/1200/5000 bp), well
    spaced (intergenic deserts + introns → the depleted training nodes). Expressed exons + their captured gDNA →
    the enriched training nodes; unexpressed exons → pure-gDNA enriched training nodes."""
    genes, x, INTRON, GAP, SIZES = [], 5000, 3000, 6000, (300, 1200, 5000)
    for k in range(n_genes):
        s = x
        ex = []
        for sz in SIZES:
            ex.append((s, s + sz)); s = s + sz + INTRON
        genes.append((f"G{k}", "+" if k % 2 == 0 else "-", ex, 100.0 if k % 2 == 0 else 0.0))
        x = s + GAP
    return genes, x + GAP


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", type=float, default=0.99)
    ap.add_argument("--n-genes", type=int, default=10)
    ap.add_argument("--gdna", type=float, default=0.3)
    ap.add_argument("--no-capture", action="store_true")
    ap.add_argument("--boundaries", action="store_true", help="include clean-exon boundary crossings (off by default)")
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()

    genes, glen = _default_genes(args.n_genes)
    dbg: dict = {}
    rdf, solved, truth, tg, tr, sig = run(
        "gdna_prior", genes, kappa=args.kappa, n_rna=12000 * args.n_genes, gdna_fraction=args.gdna,
        capture=not args.no_capture, capture_strength=20.0, genome_length=glen, seed=args.seed, _debug=dbg,
    )
    kappa_fit = float(dbg["rna_sense_frac"])

    sub = build_training_substrate(
        dbg["chain"], dbg["belief"], dbg["geometry"], dbg["statics"],
        dbg["region_arrays"], dbg["boundary_substrate"],
        include_boundaries=args.boundaries,
    )
    print(f"scenario: kappa_fit={kappa_fit:.3f} capture={not args.no_capture} gdna={args.gdna} "
          f"n_genes={args.n_genes}")
    print(f"substrate: n={sub.n} training nodes  n_eff={sub.n_eff:.1f}")
    for code, name in KIND_NAMES.items():
        m = sub.node_kind == code
        if m.any():
            print(f"  {name:11s} n={int(m.sum()):4d}  log_rho[{sub.log_rho[m].min():+.2f},"
                  f"{sub.log_rho[m].max():+.2f}]  Σw={sub.weight[m].sum():.1f}")

    # oracle per-region gDNA density = true gDNA count / gDNA eff-len (log). The reference the fit must match.
    reff = np.asarray(dbg["region_eff_len"], dtype=np.float64)
    oracle_rho = np.where((tg > 0) & (reff > 0), tg / np.maximum(reff, 1e-9), np.nan)
    oracle_log = np.log(oracle_rho[np.isfinite(oracle_rho) & (oracle_rho > 0)])

    priors = {}
    for bw in ("silverman", "lscv", 0.5):
        try:
            priors[str(bw)] = GdnaDensityPrior.fit(sub, bandwidth=bw)
        except Exception as e:  # noqa: BLE001
            print(f"  [{bw}] fit failed: {e}")
    for bw, pr in priors.items():
        modes = ", ".join(f"{x:+.2f}(logρ)" for x, _ in pr.modes[:4])
        print(f"  [{bw}] h={pr.bandwidth:.3f}  modes: {modes or '(none)'}")

    # ---- plot ----
    fig, ax = plt.subplots(figsize=(10, 6))
    for code, name in KIND_NAMES.items():
        m = sub.node_kind == code
        if m.any():
            ax.plot(sub.log_rho[m], np.full(int(m.sum()), -0.02), "|", color=_KIND_COLOR[code],
                    markersize=12, alpha=0.6, label=f"training: {name} (n={int(m.sum())})")
    if oracle_log.size:
        ax.hist(oracle_log, bins=30, density=True, alpha=0.18, color="gray",
                label="oracle region gDNA density")
    for bw, pr in priors.items():
        ax.plot(pr.x_grid, np.exp(pr.logP_grid), lw=2, label=f"P̂(logρ_g) [{bw}, h={pr.bandwidth:.2f}]")
        for x, lp in pr.modes[:4]:
            ax.plot([x], [np.exp(lp)], "v", color="black", markersize=7)
    ax.set_xlabel("log ρ_g  (gDNA density)")
    ax.set_ylabel("density")
    ax.set_title(f"Phase-2 gDNA density prior — κ_fit={kappa_fit:.2f}, "
                 f"capture={'on' if not args.no_capture else 'off'}, gdna={args.gdna}")
    ax.legend(fontsize=8, loc="upper right")
    # cap y at ~2× the smooth (fixed-h) curve so the tight auto-bandwidth spikes don't squash the structure
    smooth = priors.get("0.5")
    ytop = 2.0 * float(np.exp(smooth.logP_grid).max()) if smooth is not None else None
    ax.set_ylim(bottom=-0.05, top=ytop)
    fig.tight_layout()
    out = SCRATCH / f"gdna_prior_k{args.kappa}_cap{0 if args.no_capture else 1}_g{args.gdna}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=110)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
