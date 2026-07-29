"""gDNA HYPERPRIOR vs ORACLE — across the FULL ambig suite (all conditions), the step-1 landscape.

For every condition: run the pass-0 INITIAL solve (production `calibrate`, refit_iters=0), fit the gDNA
hyperprior NPMLE on the DECONVOLVED gDNA (`ĝ = f_g·M`, precision-weighted by `var_gdna`, region nodes only —
boundaries + AMBIG excluded ⇒ non-circular), and compare to the ORACLE NPMLE (fit on the TRUE per-node gDNA).
Two fitted variants are drawn/scored vs oracle: SELECTED (`prec,−bnd`) and HYBRID (`+bg` aggregate background).

Outputs:
  * a plot GRID (one panel/condition): SELECTED + HYBRID + ORACLE densities `P(log10 ρ_g)`.
  * a RANKED table (worst L1 first) → the conditions whose hyperprior most deviates from oracle = the
    diagnostic targets for step 2.

    OMP_NUM_THREADS=1 python scripts/debug/hyperprior_suite_vs_oracle.py [--suite DIR] [--out DIR]
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
import sys

sys.path.insert(0, str(Path("/Users/mkiyer/proj/rigel/scripts/debug")))
from selfsolve_diag import _scan_and_truth  # noqa: E402
from flagship_interrogate import _oracle_per_node, _solve  # noqa: E402
from hyperprior_fit_options import _GRID, _dens_on_grid, _fit, _true_bg  # noqa: E402

from rigel.calibration.background_reference import measure_background  # noqa: E402
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.effective_length import region_eff_length  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12


def _factors(cond: str) -> dict:
    t = cond.split("_")
    g = t[1] if len(t) > 1 else "?"
    get = lambda k: (t[t.index(k) + 1] if k in t and t.index(k) + 1 < len(t) else "?")  # noqa: E731
    return dict(gdna=g, ss=get("ss"), nrna=get("nrna"), capture=get("capture"))


def _eval_condition(suite, cond, index, cfg, work, cache, ra):
    inp = _scan_and_truth(suite, cond, index, cfg, work, cache)
    dbg = _solve(inp, ra, 0)  # pass-0 INITIAL solve (prior-free), production toggles
    chain, cap, belief = dbg["chain"], dbg["capture"], dbg["belief"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fg = np.asarray(belief.f_g, float)
    var_g = np.asarray(belief.var_gdna, float)
    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)
    fp, fn = np.asarray(cap["free_pos"], bool), np.asarray(cap["free_neg"], bool)
    isr = np.asarray(chain.kind) == REGION
    cls = np.where(fp & fn, "AMBIG", np.where(fp | fn, "single", "gonly"))
    cls = np.where(~isr, "bndry", cls)
    live = (eff > 1e-9 * 1.001) & (mass > _EPS)
    g_hat = fg * mass
    ridx = np.asarray(chain.ref_idx, np.int64)
    sig_arr = index.nodes_df["signature"].to_numpy()
    intergenic = isr & (ridx < sig_arr.shape[0]) & (sig_arr[np.clip(ridx, 0, sig_arr.shape[0] - 1)] == 0)
    region = live & isr
    sel = region & np.isin(cls, ["single", "gonly"])   # no AMBIG, no boundary (non-circular)
    sel_bg = sel & ~intergenic                          # pooled intergenic → aggregate cell
    true_gfrac = G[region].sum() / max((G[region].sum() + R[region].sum()), 1.0)
    bg = measure_background(
        CalibrationSubstrate.from_payload(inp["payload"], ra), ra,
        region_eff_length(ra.region_size_bp, np.asarray(inp["gdna_fl_pmf"])),
    )
    oracle = _fit(G, eff, region & ~intergenic, background=_true_bg(G, eff, intergenic & region))
    selected = _fit(g_hat, eff, sel, var_g)
    hybrid = _fit(g_hat, eff, sel_bg, var_g, background=bg)
    od = _dens_on_grid(oracle) if oracle is not None else None
    l1 = lambda pr: (float(np.abs(_dens_on_grid(pr) - od).sum()) if (pr is not None and od is not None) else np.nan)  # noqa: E731
    return dict(
        oracle=oracle, selected=selected, hybrid=hybrid,
        l1_selected=l1(selected), l1_hybrid=l1(hybrid),
        true_gfrac=float(true_gfrac), n_train=int(sel.sum()),
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
    conds = sorted(d.name for d in suite.iterdir()
                   if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_"))
    outdir = Path(a.out) if a.out else Path(os.environ.get("RIGEL_SCRATCH", "/tmp"))
    outdir.mkdir(parents=True, exist_ok=True)

    n = len(conds)
    ncol = 4
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.2 * nrow))
    axes = np.asarray(axes).ravel()
    rows = []
    for ax, cond in zip(axes, conds):
        r = _eval_condition(suite, cond, index, cfg, work, cache, ra)
        f = _factors(cond)
        rows.append(dict(condition=cond, **f, **{k: r[k] for k in ("l1_selected", "l1_hybrid", "true_gfrac", "n_train")}))
        if r["oracle"] is not None:
            ax.plot(_GRID, _dens_on_grid(r["oracle"]), "k", lw=2.2, label="ORACLE")
        if r["selected"] is not None:
            ax.plot(_GRID, _dens_on_grid(r["selected"]), "tab:blue", lw=1.6, label=f"SELECTED L1={r['l1_selected']:.2f}")
        if r["hybrid"] is not None:
            ax.plot(_GRID, _dens_on_grid(r["hybrid"]), "tab:green", lw=1.6, ls="--", label=f"HYBRID L1={r['l1_hybrid']:.2f}")
        ax.set_title(f"{f['gdna']} ss{f['ss']} {f['nrna']} cap-{f['capture']}\n[gDNAfrac={r['true_gfrac']:.2f}]", fontsize=7.5)
        ax.set_xlim(-8, 3); ax.tick_params(labelsize=6); ax.legend(fontsize=5.5, loc="upper left")
        print(f"{cond:52s} L1_sel={r['l1_selected']:.3f} L1_hyb={r['l1_hybrid']:.3f} gfrac={r['true_gfrac']:.2f}", flush=True)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("gDNA hyperprior NPMLE vs ORACLE across the ambig suite — P(log10 ρ_g). L1=distance to oracle (lower=better).", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    figpath = outdir / "hyperprior_suite_vs_oracle.png"
    fig.savefig(figpath, dpi=110)

    df = pd.DataFrame(rows).sort_values("l1_selected", ascending=False)
    tsv = outdir / "hyperprior_suite_vs_oracle.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    pd.set_option("display.width", 200)
    print("\n=== RANKED by L1_selected (worst first = diagnostic targets) ===")
    print(df[["condition", "gdna", "ss", "nrna", "capture", "true_gfrac", "l1_selected", "l1_hybrid"]].to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\nfig  -> {figpath}\ntable-> {tsv}")


if __name__ == "__main__":
    main()
