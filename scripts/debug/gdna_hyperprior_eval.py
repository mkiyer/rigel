"""Role-B gDNA hyperprior — the STAGE verification harness (the ruler; docs/CARRY_FORWARD.md §6).

Runs the REAL calibrate end-to-end on the cached ambig_dense_10mb suite and, for the current code, emits:
  (a) docs/figures/<stage>_fit_vs_oracle.png — the fitted gDNA prior P(log10 rho_g) vs the ORACLE landscape,
      per scenario, with node 1055's true / pass-0 / accessible-max marked and the none_* pole shown.
  (b) a refit-vs-oracle metric table on stdout: per-(ss,cap,node-type) mass-weighted f_g L1-to-oracle, plus
      POLE 1 (node 1055, must not collapse), POLE 2 (none_* captured exons, must stay low), INVARIANT
      (ss0.99 stranded exons, must not move), and a JSON dump for A/B across stages.

Usage:  OMP_NUM_THREADS=1 python scripts/debug/gdna_hyperprior_eval.py --stage baseline
A/B a code change: run with --stage before, change code, run with --stage after; diff the JSONs / plots."""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402
from flagship_interrogate import _oracle_per_node  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12
LN10 = np.log(10.0)
GRID = np.linspace(-5.0, 2.5, 320)  # log10 rho_g, for the landscape plots
H = 0.15  # decades, the oracle-render bandwidth
TYPE = {0: "intergenic", 1: "intron", 2: "exon"}

ap = argparse.ArgumentParser()
ap.add_argument("--stage", default="baseline", help="label for figure/JSON filenames")
ap.add_argument("--outdir", default="/Users/mkiyer/proj/rigel/docs/figures")
ap.add_argument("--additive", action="store_true", help="use the additive Role-B prior (D1)")
ap.add_argument("--enrich", action="store_true", help="enrichment-condition the Role-B prior (Stage 1)")
a = ap.parse_args()

suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index"))
ra = RegionArrays.from_index(index)
cfg = PipelineConfig()
work = Path("/tmp/rigel_selfsolve")
cache = suite / "_selfsolve_cache"
conds = sorted(
    d.name for d in suite.iterdir()
    if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_")
)
PANELS = [
    "gdna_gdna100_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_none_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.99_nrna_none_capture_on",
]


def run(inp, iters):
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=iters, gdna_prior_additive=a.additive,
                             gdna_prior_enrichment_condition=a.enrich)
    dbg = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    return dbg


def render_unit(a_cells):
    z = (GRID[:, None] - a_cells[None, :]) / H
    d = (np.exp(-0.5 * z * z) / (H * np.sqrt(2 * np.pi))).sum(axis=1)
    return d / max(float(d.sum()), _EPS)


def prior_on_grid(hp):
    """Fitted DensityNPMLE.logP (natural-log rho axis) -> normalized density on the log10 GRID.
    Outside the fitted grid the density is UNDEFINED (the solver clamps, but rendering that as a flat
    plateau is misleading) -> return NaN so the plot shows the prior only where it is actually fit."""
    if hp is None:
        return None
    lr10 = np.asarray(hp.log_rho, float) / LN10
    lp = np.interp(GRID, lr10, np.asarray(hp.logP, float), left=hp.logP[0], right=hp.logP[-1])
    d = np.exp(lp - lp.max())
    d = d / max(float(d.sum()), _EPS)
    d[(GRID < lr10.min()) | (GRID > lr10.max())] = np.nan  # don't draw the interp clamp
    return d


buck = defaultdict(lambda: np.zeros(3))  # (ss,cap,type) -> [sum|refit-orc|*w, w, unused]
poles = {}
panels_data = {}
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, work, cache)
    d0 = run(inp, 0)
    d1 = run(inp, 1)
    chain = d0["chain"]
    cap = d0["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G = Gp + Gn
    Rt = Rp + Rn
    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)
    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    isr = np.asarray(chain.kind) == REGION
    idx = np.asarray(chain.ref_idx, np.int64)
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    ntype = np.where(isr, rtype[np.clip(idx, 0, len(rtype) - 1)], -1)
    single = fp ^ fn
    live = (eff > 1e-9) & (mass > _EPS)
    fo = np.divide(G, np.maximum(G + Rt, _EPS))
    f0 = np.asarray(d0["belief"].f_g, float)
    f1 = np.asarray(d1["belief"].f_g, float)
    v1 = np.asarray(d1["belief"].var_gdna, float)

    ss = "0.99" if "0.99" in cond else "0.50"
    capt = "capON" if "capture_on" in cond else "capOFF"
    for t in (0, 1, 2):
        m = live & isr & (ntype == t)
        if m.any():
            w = mass[m]
            buck[(ss, capt, TYPE[t])] += [float(np.sum(np.abs(f1[m] - fo[m]) * w)), float(np.sum(w)), 0.0]

    # POLE 1: node 1055 (max-G captured single exon, true f_g>0.6) — must not collapse
    crush = live & isr & (ntype == 2) & single & (fo > 0.6)
    if crush.any() and "gdna100" in cond and ss == "0.50" and capt == "capON":
        i = int(np.argmax(G * crush))
        poles["pole1_node1055"] = dict(cond=cond, node=i, true=float(fo[i]), pass0=float(f0[i]),
                                       refit=float(f1[i]), var=float(v1[i]))
    # POLE 2: none_* captured single exons — must stay low
    ex = live & isr & (ntype == 2) & single
    if "gdna_none" in cond and capt == "capON" and ex.any():
        poles.setdefault("pole2_none", {})[cond] = float(np.mean(f1[ex]))
    # INVARIANT: ss0.99 captured exons L1 already in buck; record node-1055 equivalent
    if crush.any() and "gdna100" in cond and ss == "0.99" and capt == "capON":
        i = int(np.argmax(G * crush))
        poles["invariant_node1055_ss099"] = dict(true=float(fo[i]), refit=float(f1[i]))

    if cond in PANELS:
        dO = render_unit(np.log10(np.maximum(G[live & isr], 1.0)) - np.log10(eff[live & isr]))
        dP = prior_on_grid(d1.get("gdna_hyperprior"))
        marks = None
        if crush.any():
            i = int(np.argmax(G * crush))
            rt = mass[i] / eff[i]
            marks = dict(true=float(np.log10(max(fo[i] * rt, 1e-9))),
                         p0=float(np.log10(max(f0[i] * rt, 1e-9))), acc=float(np.log10(rt)),
                         fo=float(fo[i]), f1=float(f1[i]))
        panels_data[cond] = (dO, dP, marks)

# ---- landscape figure ----
outdir = Path(a.outdir)
outdir.mkdir(parents=True, exist_ok=True)
fig, axes = plt.subplots(2, 2, figsize=(13, 8))
axes = axes.ravel()
for ax, cond in zip(axes, PANELS):
    dO, dP, marks = panels_data[cond]
    ax.fill_between(GRID, dO, color="0.82", label="ORACLE landscape")
    ax.plot(GRID, dO, "k", lw=1.4)
    if dP is not None:
        ax.plot(GRID, dP, "tab:red", lw=1.8, label="fitted gDNA prior")
    if marks is not None:
        ax.axvline(marks["true"], color="tab:green", lw=2.0, label=f"node TRUE ({marks['fo']:.2f})")
        ax.axvline(marks["p0"], color="tab:orange", lw=1.4, ls="--", label="node PASS-0")
        ax.axvline(marks["acc"], color="0.5", lw=1.0, ls=":", label="accessible max")
        ax.set_title(f"{cond.replace('gdna_', '')}\n refit f_g={marks['f1']:.3f} (true {marks['fo']:.2f})",
                     fontsize=8)
    else:
        ax.set_title(cond.replace("gdna_", ""), fontsize=8)
    ax.set_xlim(-5, 2.5)
    ax.set_xlabel("log10 gDNA density", fontsize=8)
    ax.legend(fontsize=6, loc="upper right")
    ax.tick_params(labelsize=7)
fig.suptitle(f"[{a.stage}] fitted gDNA prior (red) vs ORACLE landscape (gray) + refit f_g at the poles",
             fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.97))
figpath = outdir / f"{a.stage}_fit_vs_oracle.png"
fig.savefig(figpath, dpi=120)

# ---- metric table + JSON ----
print(f"\n===== STAGE: {a.stage} — refit-vs-oracle =====")
print(f"{'ss':>4s} {'cap':>6s} {'type':>11s} | {'f_g L1 (mass-wtd)':>18s}")
table = {}
for k in sorted(buck):
    s = buck[k]
    l1 = s[0] / max(s[1], _EPS)
    table["/".join(k)] = round(l1, 4)
    print(f"{k[0]:>4s} {k[1]:>6s} {k[2]:>11s} | {l1:18.4f}")

p1 = poles.get("pole1_node1055", {})
p2 = poles.get("pole2_none", {})
inv = poles.get("invariant_node1055_ss099", {})
print("\n--- POLES / INVARIANT ---")
if p1:
    verdict = "CRUSH" if p1["refit"] < 0.5 * p1["pass0"] else "held"
    print(f"POLE 1 node {p1['node']}: true={p1['true']:.3f} pass0={p1['pass0']:.3f} "
          f"refit={p1['refit']:.4f} (var {p1['var']:.2f})  -> {verdict}")
if p2:
    print("POLE 2 none_* captured-exon mean f_g (want ~0): "
          + "  ".join(f"{c.replace('gdna_none_', '')}={v:.3f}" for c, v in sorted(p2.items())))
if inv:
    print(f"INVARIANT ss0.99 node-1055-eq: true={inv['true']:.3f} refit={inv['refit']:.3f} "
          f"(want refit~=true)")

jsonpath = outdir / f"{a.stage}_metrics.json"
jsonpath.write_text(json.dumps({"stage": a.stage, "table": table, "poles": poles}, indent=2))
print(f"\nfig  -> {figpath}\njson -> {jsonpath}")
