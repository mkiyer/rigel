"""NPMLE DIAGNOSTICS — the two faces of the pass-0 NPMLE (mean + variance), and how STRONG it is as a prior
relative to the belief-propagation messages.

Per condition, two panels:

  (A) the NPMLE's two faces on ONE density axis:
        * MEAN     — P(log ρ), the fitted enrichment landscape (left y-axis);
        * VARIANCE — var_proj(observed density) = the projection variance, whose floor is h² and whose spikes
          mark mode ambiguity (right y-axis). This is what the message precision reads.

  (B) STRENGTH, as pseudo-observations (log y): three per-node distributions —
        * PRIOR n_eff       — how many pseudo-obs the NPMLE prior pins f_g with (want ≲ a read; production ~0.15);
        * MESSAGE precision — the combined fwd+bwd gDNA message precision a node RECEIVES (from `_capture`);
        * projection CAP    — 1/var_proj, the MAX precision any message through this node may carry (the
          count-zero-information ceiling).
     If PRIOR n_eff ≪ MESSAGE precision, the messages dominate the prior (as intended); the CAP is the ceiling
     the σ²_transfer damping imposes on the messages.

    OMP_NUM_THREADS=1 python npmle_strength_diag.py [--conditions a,b] [--out DIR]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.insert(0, str(Path("/Users/mkiyer/proj/rigel/scripts/debug")))
from selfsolve_diag import _scan_and_truth  # noqa: E402
from npmle_movie_frame import solve_belief  # runs the real solve; stashes prior/cap on .last  # noqa: E402

from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12
_LN10 = np.log(10.0)


def prior_neff(prior, mass, eff, L=10.0, K=400):
    """Per-node prior strength in pseudo-observations: 1/Var(λ) of the prior-alone f_g posterior, in units of
    one balanced BB observation's curvature (¼) — the same measure `npmle_variance._prior_strength` uses."""
    lam = np.linspace(-L, L, K)
    fg = 1.0 / (1.0 + np.exp(-lam))
    t = prior.logprior(fg, mass, eff)
    w = np.exp(t - t.max(axis=1, keepdims=True))
    w /= w.sum(axis=1, keepdims=True)
    m = (w * lam[None, :]).sum(1)
    var = (w * (lam[None, :] - m[:, None]) ** 2).sum(1)
    return (1.0 / np.maximum(var, 1e-9)) / 0.25


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    outdir = Path(a.out) if a.out else Path(os.environ.get("RIGEL_SCRATCH", "/tmp"))
    outdir.mkdir(parents=True, exist_ok=True)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    conds = a.conditions.split(",") if a.conditions else sorted(p.stem for p in cache.glob("*.pkl"))

    n = len(conds)
    fig, axes = plt.subplots(n, 2, figsize=(14, 3.0 * n))
    axes = np.atleast_2d(axes)
    for row, c in enumerate(conds):
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        solve_belief(inp, index, cfg)  # runs the solve; stashes prior/cap
        st = solve_belief.last
        prior, capd, mass, eff = st["prior"], st["cap"], st["mass_g"], st["eff_g"]
        e = np.maximum(eff, _EPS)
        live = e > 1e-9 * 1.001
        h = float(prior.bandwidth) * _LN10

        # (A) MEAN P(logρ) + VARIANCE var_proj(density)
        axA = axes[row, 0]
        x = prior.log_rho / _LN10
        P = np.exp(prior.logP - prior.logP.max())
        axA.fill_between(x, 0, P / P.max(), color="0.82", label="MEAN P(log ρ) (scaled)")
        og = np.linspace(prior.log_rho[0] + 0.3, prior.log_rho[-1] - 0.3, 200)
        _mu, vp = prior.project(np.exp(og) * 1000.0, np.full_like(og, 1000.0))  # density = exp(og)
        axV = axA.twinx()
        axV.plot(og / _LN10, vp, "k-", lw=1.8, label="var_proj (right axis)")
        axV.axhline(h * h, color="tab:green", ls=":", lw=1.2, label=f"floor h²={h * h:.2f}")
        axA.set_title(f"{c.replace('gdna_gdna300_', '')}: NPMLE mean + var_proj", fontsize=8)
        axA.set_xlabel("log₁₀ ρ", fontsize=8)
        axA.set_ylabel("P(log ρ)", fontsize=8)
        axV.set_ylabel("var_proj (nats²)", fontsize=8)
        axA.legend(fontsize=6, loc="upper left")
        axV.legend(fontsize=6, loc="upper right")

        # (B) strength distributions
        axB = axes[row, 1]
        neff = prior_neff(prior, mass[live], eff[live])
        _mu2, vpn = prior.project(mass[live], eff[live])
        cap = 1.0 / np.maximum(vpn, 1e-9)
        prec = np.asarray(capd["prec_g"], float)[live]
        prec = prec[prec > _EPS]
        bins = np.logspace(-2, 3, 50)
        axB.hist(np.clip(neff, 1e-2, 1e3), bins=bins, alpha=0.55, color="tab:purple", label=f"prior n_eff (med {np.median(neff):.2f})")
        axB.hist(np.clip(prec, 1e-2, 1e3), bins=bins, alpha=0.55, color="tab:orange", label=f"msg precision (med {np.median(prec):.2f})")
        axB.hist(np.clip(cap, 1e-2, 1e3), bins=bins, alpha=0.4, color="tab:blue", label=f"proj CAP 1/var (med {np.median(cap):.2f})")
        axB.set_xscale("log")
        axB.set_title("strength (pseudo-obs): prior vs message vs cap", fontsize=8)
        axB.set_xlabel("pseudo-observations", fontsize=8)
        axB.legend(fontsize=6)
        print(f"{c:52s} prior_neff med={np.median(neff):.2f}  msg med={np.median(prec):.2f}  "
              f"cap med={np.median(cap):.2f}", flush=True)
    fig.suptitle("NPMLE diagnostics — mean+variance (left) and prior vs message strength (right)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    out = outdir / "npmle_strength_diag.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
