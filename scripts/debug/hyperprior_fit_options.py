"""gDNA HYPERPRIOR FIT OPTIONS — after the Phase-1a first pass, fit the deconvolved-gDNA NPMLE under different
training-node choices and compare each to the ORACLE NPMLE (fit on the TRUE per-node gDNA densities).

Per scenario, five densities `P(log ρ_g)` are drawn on one panel:
  * ORACLE  (bold black) — DensityNPMLE fit on the TRUE gDNA count G per REGION node (the target the hyperprior
    should recover);
  * the 2×2 fit options on the SOLVED gDNA `ĝ = f_g·M` — {precision-weighted by var_gdna | flat} × {boundaries
    excluded | included}, AMBIG always excluded (the hyperprior's target ⇒ non-circular).

A distance-to-oracle (L1 between normalized densities on a common grid) is printed + annotated: lower = the fit
option better recovers the true gDNA-density distribution. Run across the DNA gradient + controls.

    OMP_NUM_THREADS=1 python hyperprior_fit_options.py [--out DIR]
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

from rigel.calibration.background_reference import measure_background  # noqa: E402
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.effective_length import region_eff_length  # noqa: E402
from rigel.calibration.npmle import DensityNPMLE  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12
_LN10 = np.log(10.0)

CONDS = [
    ("gdna_gdna300_ss_0.99_nrna_present_capture_on", "CONTROL stranded ss.99 gdna300 capON"),
    ("gdna_gdna300_ss_0.50_nrna_present_capture_off", "CONTROL capOFF gdna300 unstr"),
    ("gdna_none_ss_0.50_nrna_present_capture_on", "ZERO DNA unstr capON"),
    ("gdna_gdna1_ss_0.50_nrna_present_capture_on", "LOW 1% unstr capON"),
    ("gdna_gdna5_ss_0.50_nrna_present_capture_on", "LOW 5% unstr capON"),
    ("gdna_gdna100_ss_0.50_nrna_present_capture_on", "HIGH 100% unstr capON"),
    ("gdna_gdna300_ss_0.50_nrna_present_capture_on", "VERY HIGH 300% unstr capON (flagship)"),
]

_GRID = np.linspace(-8.0, 3.0, 400)  # common log10 ρ grid for the distance metric (wide: covers the ρ_bg floor)


def _true_bg(G, eff, mask):
    """The TRUE aggregate background from the oracle gDNA counts over the pooled (intergenic) regions — so the
    oracle is itself the ideal HYBRID (true aggregate + true individual), a like-for-like target."""
    from rigel.calibration.background_reference import BackgroundReference

    sg, se, nr = float(G[mask].sum()), float(eff[mask].sum()), int(mask.sum())
    lr = float(np.log(sg) - np.log(se)) if sg > 1e-12 else -np.inf
    return BackgroundReference(log_rho_bg=lr, sigma_bg=(1 / np.sqrt(sg) if sg > 1e-12 else np.inf),
                               n_counts=sg, eff_total=se, n_regions=nr)


def _dens_on_grid(pr):
    """Normalized density of a DensityNPMLE interpolated onto the common log10-ρ grid."""
    x = np.asarray(pr.log_rho) / _LN10
    p = np.exp(np.asarray(pr.logP) - np.asarray(pr.logP).max())
    d = np.interp(_GRID, x, p, left=0.0, right=0.0)
    s = d.sum()
    return d / s if s > 0 else d


def _fit(g, eff, mask, var=None, background=None):
    if mask.sum() < 5:
        return None
    return DensityNPMLE.fit(
        g[mask], eff[mask], var_g=(None if var is None else var[mask]),
        background=background, bandwidth=0.15,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    ra = RegionArrays.from_index(index)

    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    axes = axes.ravel()
    for ax, (cond, title) in zip(axes, CONDS):
        inp = _scan_and_truth(suite, cond, index, cfg, work, cache)
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
        g_hat = fg * mass

        ridx = np.asarray(chain.ref_idx, np.int64)
        sig_arr = index.nodes_df["signature"].to_numpy()
        intergenic = isr & (ridx < sig_arr.shape[0]) & (sig_arr[np.clip(ridx, 0, sig_arr.shape[0] - 1)] == 0)
        region = live & isr
        sel = region & np.isin(cls, ["single", "gonly"])  # no AMBIG, no boundary
        selb = sel | (live & (cls == "bndry"))  # + boundaries
        sel_bg = sel & ~intergenic  # +bg: pooled (intergenic) regions move to the aggregate cell, not individual
        true_gfrac = G[region].sum() / max((G[region].sum() + R[region].sum()), 1.0)
        # the aggregate background scalar (intergenic-only) — injected as a genome-length Poisson cell
        bg = measure_background(CalibrationSubstrate.from_payload(inp["payload"], ra), ra,
                                region_eff_length(ra.region_size_bp, np.asarray(inp["gdna_fl_pmf"])))

        fits = {
            "ORACLE (hybrid, true)": (_fit(G, eff, region & ~intergenic, background=_true_bg(G, eff, intergenic & region)), "k", 2.4, "-"),
            "prec, −bnd (SELECTED)": (_fit(g_hat, eff, sel, var_g), "tab:blue", 1.8, "-"),
            "flat, −bnd": (_fit(g_hat, eff, sel), "tab:cyan", 1.3, "--"),
            "prec, +bnd": (_fit(g_hat, eff, selb, var_g), "tab:red", 1.5, "-"),
            "prec, −bnd, +bg (HYBRID)": (_fit(g_hat, eff, sel_bg, var_g, background=bg), "tab:green", 2.0, "-"),
        }
        oracle_d = _dens_on_grid(fits["ORACLE (hybrid, true)"][0])
        labels = []
        for name, (pr, c, lw, ls) in fits.items():
            if pr is None:
                continue
            d = _dens_on_grid(pr)
            dist = float(np.abs(d - oracle_d).sum()) if name != "ORACLE (hybrid, true)" else 0.0
            lab = name if name == "ORACLE (hybrid, true)" else f"{name}  L1={dist:.2f}"
            ax.plot(_GRID, d, color=c, lw=lw, ls=ls, label=lab)
            labels.append((name, dist))
        ax.set_title(f"{title}\n[true gDNAfrac={true_gfrac:.2f}]", fontsize=9)
        ax.set_xlabel("log10 ρ_g", fontsize=8)
        ax.set_xlim(-8, 3)
        ax.legend(fontsize=6.5, loc="upper left")
        ax.tick_params(labelsize=7)
        print(f"{title[:44]:44s}  " + "  ".join(f"{n.split()[0]}:{d:.2f}" for n, d in labels if n != "ORACLE (hybrid, true)"),
              flush=True)

    axes[-1].axis("off")
    fig.suptitle(
        "gDNA hyperprior NPMLE — fit options vs the ORACLE (true-G) NPMLE, across the DNA gradient. "
        "L1 = distance to oracle (lower=better). AMBIG always excluded.",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    outdir = Path(a.out) if a.out else Path(os.environ.get("RIGEL_SCRATCH", "/tmp"))
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "hyperprior_fit_options.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
