"""FRAME-BY-FRAME gDNA-NPMLE movie: pass-0 (initial prior-free solve) → pass-1 (gDNA-hyperprior REFIT), with
the ORACLE gDNA NPMLE superimposed. "What happens to the deconvolved-gDNA density after the refit?"

Now that the Phase-2 refit is wired into `calibrate` (`calib_refit_iters`), a single `_solve(inp, ra, 1)`
exposes both the initial belief (`belief_pass0`) and the refit belief (`belief`). Per condition we fit the
SELECTED gDNA NPMLE (region, single/gonly, precision-weighted) on each, plus the oracle on the TRUE gDNA at the
same nodes, and overlay. L1-to-oracle before (pass0) vs after (refit) = the refit's payoff.

    OMP_NUM_THREADS=1 python scripts/debug/refit_gdna_movie.py [--out DIR]
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
import pandas as pd
from scipy.stats import wasserstein_distance
import sys

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402
from flagship_interrogate import _oracle_per_node, _solve  # noqa: E402
from hyperprior_fit_options import _GRID, _dens_on_grid, _fit  # noqa: E402
from hyperprior_suite_vs_oracle import _factors  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--out", default="docs/calibration/figures")
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path("/tmp/rigel_selfsolve"); cache = suite / "_selfsolve_cache"
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    conds = sorted(d.name for d in suite.iterdir()
                   if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_"))
    outdir = Path(a.out); outdir.mkdir(parents=True, exist_ok=True)

    n = len(conds); ncol = 4; nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.2 * nrow))
    axes = np.asarray(axes).ravel()
    rows = []
    for ax, cond in zip(axes, conds):
        inp = _scan_and_truth(suite, cond, index, cfg, work, cache)
        dbg = _solve(inp, ra, 1)  # pass-0 + 1 refit
        chain, cap = dbg["chain"], dbg["capture"]
        b0, b1 = dbg["belief_pass0"], dbg["belief"]
        Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
        G = Gp + Gn
        mass = np.asarray(cap["mass_global"], float)
        eff = np.asarray(cap["eff_global"], float)
        fp, fn = np.asarray(cap["free_pos"], bool), np.asarray(cap["free_neg"], bool)
        isr = np.asarray(chain.kind) == REGION
        live = (eff > 1e-9) & (mass > _EPS)
        sel = live & isr & ((fp ^ fn) | (~fp & ~fn))  # SELECTED: single or gonly, no AMBIG/boundary
        f0 = _fit(np.asarray(b0.f_g, float) * mass, eff, sel, np.asarray(b0.var_gdna, float))
        f1 = _fit(np.asarray(b1.f_g, float) * mass, eff, sel, np.asarray(b1.var_gdna, float))
        orc = _fit(G, eff, sel)
        f = _factors(cond)
        gfrac = G[live & isr].sum() / max((G[live & isr].sum() + (Rp + Rn)[live & isr].sum()), 1.0)
        od = _dens_on_grid(orc) if orc is not None else None
        # EMD (Wasserstein-1) on the log-ρ axis = the log-displacement of mass; NON-saturating (L1 saturates at
        # 2.0 on disjoint support = the degenerate zero-gDNA case). This is the metric fix.
        emd = lambda pr: (float(wasserstein_distance(_GRID, _GRID, _dens_on_grid(pr), od)) if (pr is not None and od is not None) else np.nan)  # noqa: E731
        l1 = lambda pr: (float(np.abs(_dens_on_grid(pr) - od).sum()) if (pr is not None and od is not None) else np.nan)  # noqa: E731
        l1_0, l1_1 = l1(f0), l1(f1)
        emd_0, emd_1 = emd(f0), emd(f1)
        rows.append(dict(condition=cond.replace("gdna_", "").replace("_nrna", ""), gdna=f["gdna"], ss=f["ss"],
                         cap=f["capture"], gfrac=float(gfrac), emd_pass0=emd_0, emd_refit=emd_1,
                         d_emd=emd_1 - emd_0, l1_pass0=l1_0, l1_refit=l1_1, d=l1_1 - l1_0))
        if orc is not None:
            ax.plot(_GRID, _dens_on_grid(orc), "k", lw=2.2, label="ORACLE")
        if f0 is not None:
            ax.plot(_GRID, _dens_on_grid(f0), "tab:blue", lw=1.5, ls=":", label=f"pass-0 L1={l1_0:.2f}")
        if f1 is not None:
            ax.plot(_GRID, _dens_on_grid(f1), "tab:red", lw=1.7, label=f"REFIT L1={l1_1:.2f}")
        ax.set_title(f"{f['gdna']} ss{f['ss']} cap-{f['capture']} [gf={gfrac:.2f}]", fontsize=7.5)
        ax.set_xlim(-8, 3); ax.tick_params(labelsize=6); ax.legend(fontsize=5.5, loc="upper left")
        print(f"{cond:52s} EMD pass0={emd_0:.3f} -> refit={emd_1:.3f} ({emd_1-emd_0:+.3f})", flush=True)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("gDNA-NPMLE frame-by-frame: pass-0 (dotted blue) → REFIT (red) vs ORACLE (black). L1=dist to oracle.", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(outdir / "refit_gdna_movie.png", dpi=110)

    df = pd.DataFrame(rows).sort_values("d_emd")
    df.to_csv(outdir / "refit_gdna_movie.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print("\n=== REFIT vs PASS-0 — EMD (Wasserstein, log-displacement; the non-degenerate metric) & L1 ===")
    print("    (Δ = refit − pass0; NEGATIVE = refit BETTER)")
    print(df[["condition", "gdna", "ss", "cap", "gfrac", "emd_pass0", "emd_refit", "d_emd", "d"]].to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\nBy EMD — improved={int((df.d_emd<-0.02).sum())}  worse={int((df.d_emd>0.02).sum())}  flat={int((df.d_emd.abs()<=0.02).sum())}"
          f"  (L1 hid the sparse-condition regressions the floor bug causes)")
    print(f"fig -> {outdir/'refit_gdna_movie.png'}")


if __name__ == "__main__":
    main()
