"""NPMLE SMOOTHNESS DIAGNOSTIC — why is the deconvolved-gDNA NPMLE fit BUMPY where the oracle is smooth?

Hypothesis: the bumps are sparse FALSE POSITIVES — nodes the unstranded first pass over-called (true f_g low,
solved f_g high), scattered across ρ. The oracle is smooth because the TRUE densities cluster tightly.

Three panels on one low-DNA scenario:
  A — the raw data the NPMLE fits: histogram of the deconvolved density ĝ/E per node, split into truly-gDNA vs
      FALSE-POSITIVE (over-called) nodes, with the TRUE G/E overlaid. Shows where the scatter/bumps come from.
  B — the fit vs bandwidth h ∈ {0.10, 0.15, 0.25, 0.40} against the oracle — does more smoothing remove the
      bounce, and at what cost to the modes?
  C — does precision-weighting absorb it? var_gdna (the belief width the NPMLE deconvolves) vs the deconvolved
      density, colored by false-positive — are the bumpy nodes UNCERTAIN (so the NPMLE should down-weight them)?

    OMP_NUM_THREADS=1 python npmle_smoothness_diag.py [--condition C] [--out DIR]
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
from flagship_interrogate import _oracle_per_node, _solve  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.npmle import DensityNPMLE  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12
_LN10 = np.log(10.0)


def _curve(pr, ax, **kw):
    x = np.asarray(pr.log_rho) / _LN10
    p = np.exp(np.asarray(pr.logP) - np.asarray(pr.logP).max())
    p = p / p.sum()
    ax.plot(x, p, **kw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--condition", default="gdna_gdna5_ss_0.50_nrna_present_capture_on")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    inp = _scan_and_truth(suite, a.condition, index, cfg, work, cache)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

    dbg = _solve(inp, ra, 0)
    chain = dbg["chain"]
    cap = dbg["capture"]
    belief = dbg["belief"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fg = np.asarray(belief.f_g, float)
    var_g = np.asarray(belief.var_gdna, float)
    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)
    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    isr = np.asarray(chain.kind) == REGION
    cls = np.where(fp & fn, "AMBIG", np.where(fp | fn, "single", "gonly"))
    cls = np.where(~isr, "bndry", cls)
    live = (eff > 1e-9 * 1.001) & (mass > _EPS)

    sel = live & isr & np.isin(cls, ["single", "gonly"])  # the SELECTED training set (no AMBIG/boundary)
    g_hat = fg * mass
    d = np.log10(np.maximum(g_hat[sel] / np.maximum(eff[sel], _EPS), 1e-12))  # deconvolved density
    tot = (G + R)[sel]
    fo = np.where(tot > _EPS, G[sel] / np.maximum(tot, _EPS), 0.0)
    vg = var_g[sel]
    td = np.log10(np.maximum(G[sel] / np.maximum(eff[sel], _EPS), 1e-12))  # TRUE gDNA density
    # false positive: truly ~RNA (fo<0.2) but the solve put real gDNA there (f_g>0.25)
    fp_mask = (fo < 0.2) & (fg[sel] > 0.25)
    true_g = fo > 0.5

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))
    # A — the raw data
    ax = axes[0]
    bins = np.linspace(-5, 2.5, 60)
    ax.hist(d[true_g], bins=bins, color="tab:green", alpha=0.6, label=f"truly-gDNA nodes (n={int(true_g.sum())})")
    ax.hist(d[fp_mask], bins=bins, color="tab:red", alpha=0.6, label=f"FALSE POSITIVE (over-called, n={int(fp_mask.sum())})")
    ax.hist(td[G[sel] > 0], bins=bins, histtype="step", color="k", lw=1.6, label="TRUE gDNA density G/E")
    ax.set_title("A. the data the NPMLE fits: deconvolved ĝ/E per node")
    ax.set_xlabel("log10 density"); ax.legend(fontsize=7)
    # B — bandwidth sweep vs oracle
    ax = axes[1]
    oracle = DensityNPMLE.fit(G[live & isr], eff[live & isr])
    _curve(oracle, ax, color="k", lw=2.4, label="ORACLE (true G)")
    for h, c in ((0.10, "tab:cyan"), (0.15, "tab:blue"), (0.25, "tab:orange"), (0.40, "tab:red")):
        pr = DensityNPMLE.fit(g_hat[sel], eff[sel], var_g=vg, bandwidth=h)
        _curve(pr, ax, color=c, lw=1.3, label=f"fit h={h}")
    ax.set_xlim(-5, 2.5)
    ax.set_title("B. fit vs bandwidth (precision-weighted) vs oracle")
    ax.set_xlabel("log10 ρ_g"); ax.legend(fontsize=7)
    # C — does precision absorb it? var vs density
    ax = axes[2]
    ax.scatter(d[true_g], vg[true_g], s=8, c="tab:green", alpha=0.4, label="truly-gDNA")
    ax.scatter(d[fp_mask], vg[fp_mask], s=12, c="tab:red", alpha=0.5, label="false positive")
    ax.axhline(float(np.median(vg)), color="0.5", ls="--", lw=0.8, label=f"median var={np.median(vg):.2f}")
    ax.set_title("C. is the scatter UNCERTAIN? var_gdna vs deconvolved density")
    ax.set_xlabel("log10 ĝ/E"); ax.set_ylabel("var_gdna (belief width)"); ax.legend(fontsize=7)

    fig.suptitle(f"NPMLE smoothness — {a.condition}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    outdir = Path(a.out) if a.out else Path(os.environ.get("RIGEL_SCRATCH", "/tmp"))
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "npmle_smoothness_diag.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    # numbers
    print(f"=== {a.condition} ===")
    print(f"selected training nodes n={int(sel.sum())}  truly-gDNA={int(true_g.sum())}  "
          f"FALSE-POSITIVE(over-called)={int(fp_mask.sum())} ({fp_mask.mean()*100:.0f}%)")
    print(f"false-positive median var_gdna={np.median(vg[fp_mask]):.2f} vs truly-gDNA median var={np.median(vg[true_g]):.2f} "
          f"vs all median={np.median(vg):.2f}")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
