"""How WEAK is the pass-0 all-nodes NPMLE prior? — measure its per-node strength in the solver's units.

The pass-0 prior must be EXTREMELY weak so the strand likelihood (when present) and the boundary messages
dominate and can peel RNA out of the f_g=1 start. The solve works in log-odds λ (f_g=σ(λ), λ∈[−L,L], L=10);
the prior enters as an additive term t(λ)=logP(σ(λ)·d) on each node's ψ. Its STRENGTH is how much it
concentrates λ:

  * SD_λ of the prior-alone posterior p(λ)∝exp(t(λ)). A DIFFUSE prior over [−L,L] has SD_λ = 2L/√12 ≈ 5.77
    (maximally weak — says nothing). A strong prior pins λ to SD_λ≈1.
  * pseudo-observations n_eff ≈ prior curvature / (strand curvature of ONE balanced BB observation = ¼).
    Production caps the global prior at n_eff≈1 (`_GLOBAL_STAB_PREC`) so a node's own count always wins;
    we want the NPMLE prior comparably weak (n_eff ≲ 1), else it over-fixes gDNA and RNA can't be peeled back.

Fits the ALL-nodes total-density prior (kernel), then plots SD_λ and n_eff vs node density + histograms.

    python scripts/debug/npmle_variance.py --beliefs DIR/COND.beliefs.npz[,...] --out DIR
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from npmle_fusion import fusion_fit, node_obs, smooth_logprior  # noqa: E402

_L = 10.0
_DIFFUSE_SD = 2.0 * _L / np.sqrt(12.0)  # SD_λ of a uniform prior over [−L,L] — the maximally-weak reference


def _prior_strength(logP, log_rho, dens, lam):
    """Per-node prior-alone posterior over λ: SD_λ and curvature-based pseudo-obs. dens=(n,), lam=(K,)."""
    fg = 1.0 / (1.0 + np.exp(-lam))  # σ(λ)  (K,)
    lrho = np.log(fg)[None, :] + np.log(np.maximum(dens, 1e-30))[:, None]  # (n,K) log ρ_g
    t = np.interp(lrho.ravel(), log_rho, logP, left=logP[0], right=logP[-1]).reshape(lrho.shape)  # (n,K)
    w = np.exp(t - t.max(axis=1, keepdims=True))
    w /= w.sum(axis=1, keepdims=True)
    m = (w * lam[None, :]).sum(1)
    var = (w * (lam[None, :] - m[:, None]) ** 2).sum(1)
    sd = np.sqrt(np.maximum(var, 1e-12))
    # pseudo-observations: prior λ-precision (1/var) in units of one balanced BB observation's curvature (¼).
    n_eff = (1.0 / np.maximum(var, 1e-12)) / 0.25
    return sd, m, n_eff


def run(npz_path, outdir: Path):
    cond = Path(npz_path).stem.replace(".beliefs", "")
    bc = dict(np.load(npz_path))
    bc["e_max"] = float(bc["e_max"])
    gh, E, sg, fit = node_obs(bc, 2, "total")
    g2, E2 = gh[fit], E[fit]
    dens = g2 / E2  # total density (f_g=1); count-0 → 0 (handled by the 1e-30 floor → pushed to λ edge)
    ldp = np.log(dens[dens > 0])
    ln10 = np.log(10.0)
    lo = float(np.percentile(ldp, 0.5)) - ln10
    hi = float(np.percentile(ldp, 99.5)) + 0.5 * ln10
    log_rho = np.linspace(lo, hi, 300)
    P, _ = fusion_fit(g2, E2, sg[fit], log_rho, method="kernel")
    logP = smooth_logprior(P, log_rho, h_dex=0.03)

    lam = np.linspace(-_L, _L, 400)
    # sample nonzero-density nodes for the scatter (all for the histogram stats).
    sd, mlam, n_eff = _prior_strength(logP, log_rho, dens, lam)
    nz = dens > 0
    ld10 = np.log10(dens[nz])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    hb = axes[0].hexbin(ld10, sd[nz], gridsize=60, cmap="viridis", bins="log")
    axes[0].axhline(_DIFFUSE_SD, color="red", ls="--", lw=1.5, label=f"diffuse (weakest) SD_λ={_DIFFUSE_SD:.2f}")
    axes[0].axhline(1.0, color="orange", ls=":", lw=1.2, label="strong (SD_λ≈1)")
    axes[0].set_title(f"{cond}: prior strength vs density — high SD_λ = WEAK")
    axes[0].set_xlabel("log10 total density")
    axes[0].set_ylabel("prior-alone SD_λ (log-odds)")
    axes[0].set_ylim(0, _DIFFUSE_SD * 1.1)
    axes[0].legend(fontsize=8, loc="lower left")
    fig.colorbar(hb, ax=axes[0], label="log node count")
    # histogram of pseudo-obs (n_eff): want ≲1 (weaker than a node's own count).
    axes[1].hist(np.clip(n_eff[nz], 0, 10), bins=60, color="tab:purple", alpha=0.7)
    axes[1].axvline(1.0, color="red", ls="--", lw=1.5, label="1 pseudo-obs (production cap)")
    med = float(np.median(n_eff[nz]))
    p90 = float(np.percentile(n_eff[nz], 90))
    axes[1].set_title(f"{cond}: prior pseudo-obs n_eff  (median={med:.2f}, p90={p90:.2f})")
    axes[1].set_xlabel("prior pseudo-observations n_eff (clip 10)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    out = outdir / f"variance_{cond}.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(
        f"{cond:10s} SD_λ[median={np.median(sd[nz]):.2f} (diffuse={_DIFFUSE_SD:.2f})]  "
        f"n_eff[median={med:.2f} p90={p90:.2f} frac>1={float((n_eff[nz] > 1).mean()):.2f}]  -> {out}",
        flush=True,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--beliefs", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    for p in args.beliefs.split(","):
        run(p, outdir)


if __name__ == "__main__":
    main()
