"""Zero-inflated NPMLE prototype — count-space P(ρ) vs the density-space KDE, on real cached data.

Reviewer's Option 1 (docs/calibration/gdna_prior_zero_handling.md §4): model the latent gDNA RATE as a
discrete mixture over a fixed log-rate grid, fit by Poisson EM on the confident set's (count, E) —
INCLUDING count=0. k_i ~ Poisson(ρ_j·E_i); the lowest grid point is the zero atom (its weight = π). No
bandwidth (⇒ no comb-smearing), zero native (count=0 drives weight to ρ→0), exposure native (per-node E
is the Poisson offset — for us E = the gDNA-FL eff-length, ~fragment length for boundaries).

We fit P(ρ) on the confident (≈pure-gDNA) set; the query projection (f_g deconvolution) is unchanged —
this only replaces WHAT the density prior is fit as. Plots P(ρ) vs the density KDE (fit on nonzero) to
show whether the fake boundary "enriched mode" evaporates and the zero atom captures the count-0 majority.

    python scripts/debug/npmle_prototype.py --cache C.pkl[,...] --out DIR
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import logsumexp

from rigel.calibration.gdna_density_prior import GdnaDensityPrior, TrainingSubstrate

sys.path.insert(0, str(Path(__file__).parent))
from calib_cache import load  # noqa: E402
from pass0_kde_zero import confident  # noqa: E402

_EPS = 1.0e-12


def npmle_fit(count, eff, log_rho_grid, n_iter=120, chunk=40000):
    """Poisson-EM NPMLE: weights w_j over the log-rate grid s.t. k_i ~ Σ_j w_j Poisson(ρ_j·E_i).
    Chunked over nodes (n×g matrix); the per-node lgamma(k+1) cancels in the E-step normalization."""
    rho = np.exp(log_rho_grid)  # (g,)
    g = rho.shape[0]
    n = count.shape[0]
    logw = np.full(g, -np.log(g))
    for _ in range(n_iter):
        acc = np.full(
            g, -np.inf
        )  # log Σ_i r_ij, accumulated in log space via logsumexp over chunks
        tot = 0.0
        for s in range(0, n, chunk):
            k = count[s : s + chunk][:, None]  # (c,1)
            E = eff[s : s + chunk][:, None]
            lam = rho[None, :] * E  # (c,g) expected counts
            logpois = k * np.log(rho[None, :] * E) - lam  # drop lgamma(k+1) (cancels in E-step)
            lp = logw[None, :] + logpois  # (c,g)
            r = np.exp(lp - logsumexp(lp, axis=1, keepdims=True))  # responsibilities (c,g)
            csum = r.sum(axis=0)  # (g,)
            acc = np.logaddexp(acc, np.log(np.maximum(csum, _EPS)))
            tot += r.shape[0]
        logw = acc - np.log(tot)  # M-step (normalized), in log space
        logw -= logsumexp(logw)
    return np.exp(logw)  # (g,) weights over the grid


def run(cache_path, outdir: Path):
    cond = Path(cache_path).stem
    count, eff, is_bnd = confident(load(cache_path))
    # deterministic subsample for the prototype (keeps the count-0 proportion; fit is stable at this n).
    rng = np.random.default_rng(0)
    if count.shape[0] > 100_000:
        idx = rng.choice(count.shape[0], 100_000, replace=False)
        count, eff, is_bnd = count[idx], eff[idx], is_bnd[idx]

    dens_nz = count[count > 0] / eff[count > 0]
    lo = np.log(0.1 / float(eff.max()))  # zero atom: expected count ~0.1 over the largest E ⇒ ρ→0
    hi = np.log(dens_nz.max() * 3.0)
    log_rho = np.linspace(lo, hi, 200)
    w = npmle_fit(count.astype(float), eff.astype(float), log_rho)
    pi_zero_data = float((count == 0).mean())
    pi_atom = float(w[:3].sum())  # weight in the lowest 3 grid points (the ρ→0 atom region)

    # density KDE (nonzero only — what the current substrate sees) for contrast
    x_nz = np.log(np.maximum(dens_nz, 1e-9))
    sub = TrainingSubstrate(
        log_rho=x_nz,
        weight=np.ones_like(x_nz),
        node_kind=np.zeros_like(x_nz, np.int64),
        node_index=np.arange(x_nz.shape[0]),
        log_rho_std=np.sqrt(1.0 / (count[count > 0] + 1.0)),
    )
    kde = GdnaDensityPrior.fit(sub, bandwidth="silverman")

    fig, ax = plt.subplots(1, 1, figsize=(11, 5.5))
    lr10 = log_rho / np.log(10)
    ax.bar(
        lr10,
        w,
        width=(lr10[1] - lr10[0]) * 0.9,
        alpha=0.55,
        color="tab:red",
        label=f"NPMLE P(ρ)  [zero-atom π={pi_atom:.2f}]",
    )
    xs = np.linspace(x_nz.min(), x_nz.max(), 400)
    kde_pdf = np.exp(kde.logpdf_kernel(xs))
    kde_pdf = kde_pdf / kde_pdf.sum() * (kde_pdf.size / 400.0)  # rough scale to overlay shape
    ax.plot(
        xs / np.log(10),
        kde_pdf / kde_pdf.max() * w.max(),
        "k-",
        lw=1.8,
        label="density KDE (nonzero, scaled) — fake enriched mode?",
    )
    ax.set_title(
        f"{cond}: count-space NPMLE P(ρ) vs density KDE   (count-0 in data = {pi_zero_data:.2f})"
    )
    ax.set_xlabel("log10 gDNA rate ρ")
    ax.set_ylabel("weight")
    ax.legend(fontsize=9)
    fig.tight_layout()
    p = outdir / f"npmle_{cond}.png"
    fig.savefig(p, dpi=110)
    plt.close(fig)
    print(
        f"{cond:10s} n={count.shape[0]} count0_data={pi_zero_data:.3f} NPMLE_zero_atom={pi_atom:.3f} "
        f"-> {p}",
        flush=True,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    for c in args.cache.split(","):
        run(c, outdir)


if __name__ == "__main__":
    main()
