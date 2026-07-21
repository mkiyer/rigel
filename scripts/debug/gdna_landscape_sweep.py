"""GOAL-1 sweep engine: find the DNA-landscape recipe whose shape is closest to the ORACLE and ROBUST across
the full scenario spectrum (stranded/unstranded x capture on/off x nascent x DNA 0..high). Pass-0 runs ONCE per
scenario (cached); the recipe sweep (subset x fit x smoothing x ρ_bg) is cheap on top.

Dimensions explored:
  FIT      kde (point-estimate, floors at 1/E) | poisson (additive, zero-native)
  SMOOTH   post-convolve the density by a Gaussian of h_smooth decades (fixes the sparse-node spikiness)
  SUBSET   single | single|gonly | +zero-count-struct | variance percentiles | structural-only(invariants)
  ρ_bg     add the pooled background Σg/ΣE over structural nodes as a sharp depleted anchor (zero-DNA)

Metric = EMD(fit, oracle) in log10 ρ_g. Reports MEAN and MAX (worst group) per recipe, ranked by MAX (robust).
"""
from __future__ import annotations

import dataclasses
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.special import gammaln  # noqa: E402
from scipy.stats import wasserstein_distance  # noqa: E402

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
GRID = np.linspace(-5.0, 2.5, 260)
_lnat = GRID * LN10
_dg = GRID[1] - GRID[0]


def smooth_matrix(h_dec):
    if h_dec <= 0:
        return None
    d = (GRID[:, None] - GRID[None, :]) / h_dec
    k = np.exp(-0.5 * d * d)
    return k / k.sum(0, keepdims=True)


def norm(d):
    return d / max(float(d.sum()), _EPS)


def poisson_density(g, E, h_smooth=0.0, bg=None):
    """Additive Poisson (zero-native) + optional pooled-background anchor + optional post-smoothing."""
    parts = []
    if g.size:
        lam = np.exp(_lnat)[None, :] * E[:, None]
        ll = g[:, None] * np.log(np.maximum(lam, _EPS)) - lam - gammaln(g[:, None] + 1.0)
        ll -= ll.max(1, keepdims=True)
        pn = np.exp(ll)
        pn /= np.maximum(pn.sum(1, keepdims=True), _EPS)
        parts.append(pn.sum(0))
    if bg is not None:  # aggregate Σg over ΣE as ONE Poisson obs, weighted by n_bg (pooled precision)
        gsum, esum, w = bg
        lam = np.exp(_lnat) * esum
        ll = gsum * np.log(np.maximum(lam, _EPS)) - lam - gammaln(gsum + 1.0)
        pnb = np.exp(ll - ll.max())
        parts.append(w * pnb / max(pnb.sum(), _EPS))
    if not parts:
        return None
    d = np.sum(parts, 0)
    S = smooth_matrix(h_smooth)
    if S is not None:
        d = S @ d
    return norm(d)


def kde_density(g, E, h=0.15, bg=None):
    if not g.size and bg is None:
        return None
    a = np.log10(np.maximum(g, 1.0)) - np.log10(E) if g.size else np.array([])
    d = np.zeros_like(GRID)
    if a.size:
        z = (GRID[:, None] - a[None, :]) / h
        d = (np.exp(-0.5 * z * z) / (h * np.sqrt(2 * np.pi))).sum(1)
    if bg is not None:
        gsum, esum, w = bg
        ab = np.log10(max(gsum, 1.0)) - np.log10(esum)
        d = d + w * np.exp(-0.5 * ((GRID - ab) / h) ** 2) / (h * np.sqrt(2 * np.pi))
    return norm(d)


# ---- pass-0 once per scenario, cache node arrays ----
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index"))
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
cfg = PipelineConfig()
work = Path("/tmp/rigel_selfsolve")
cache = suite / "_selfsolve_cache"
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_"))
rtype_all = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)

scen = []
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, work, cache)
    dbg = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]),
              dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    chain = dbg["chain"]
    cap = dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    isr = np.asarray(chain.kind) == REGION
    idx = np.asarray(chain.ref_idx, np.int64)
    ntype = np.where(isr, rtype_all[np.clip(idx, 0, len(rtype_all) - 1)], -1)
    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)
    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    f0 = np.asarray(dbg["belief"].f_g, float)
    var = np.asarray(dbg["belief"].var_gdna, float)
    dbkt = next((k for k in ("none", "gdna1000", "gdna300", "gdna100", "gdna5", "gdna1") if k in cond), "?")
    scen.append(dict(
        cond=cond, G=Gp + Gn, mass=mass, eff=eff, g_hat=f0 * mass, var=var, ntype=ntype, isr=isr,
        single=fp ^ fn, gonly=~fp & ~fn, ambig=fp & fn,
        group=("ON" if "capture_on" in cond else "OFF", dbkt),
    ))


# ---- subset builders ----
def base(s):
    return (s["eff"] > 1e-9) & s["isr"] & ~s["ambig"] & (s["mass"] > _EPS)


def struct_zero(s):
    return (s["eff"] > 1e-9) & s["isr"] & ((s["ntype"] == 0) | (s["ntype"] == 1)) & (s["mass"] <= _EPS)


def vpct(s, sel, p):
    b = base(s)
    thr = np.percentile(s["var"][b], p) if b.any() else np.inf
    return sel & (s["var"] <= thr)


SUBSETS = {
    "single": lambda s: base(s) & s["single"],
    "sgonly": lambda s: base(s) & (s["single"] | s["gonly"]),
    "sgz": lambda s: (base(s) & (s["single"] | s["gonly"])) | struct_zero(s),
    "sgz_v50": lambda s: vpct(s, (base(s) & (s["single"] | s["gonly"])), 50) | struct_zero(s),
    "sgz_v35": lambda s: vpct(s, (base(s) & (s["single"] | s["gonly"])), 35) | struct_zero(s),
    "sgz_v25": lambda s: vpct(s, (base(s) & (s["single"] | s["gonly"])), 25) | struct_zero(s),
    "sgz_v15": lambda s: vpct(s, (base(s) & (s["single"] | s["gonly"])), 15) | struct_zero(s),
    "sgz_v10": lambda s: vpct(s, (base(s) & (s["single"] | s["gonly"])), 10) | struct_zero(s),
    "struct_only": lambda s: (base(s) & s["gonly"]) | struct_zero(s),
}


def rho_bg(s):
    """Pooled background Σg/ΣE over structural (intergenic/intron/gonly) nodes, incl zero-count."""
    m = (s["eff"] > 1e-9) & s["isr"] & ((s["ntype"] == 0) | (s["ntype"] == 1) | s["gonly"])
    return (float(s["g_hat"][m].sum()), float(s["eff"][m].sum()), 3.0) if m.any() else None


# ---- recipes: (subset, fit, kwargs, use_rho_bg) ----
def R(sub, fit, sm=0.0, bg=False):
    return dict(sub=sub, fit=fit, sm=sm, bg=bg)


RECIPES = {
    "single/kde": R("single", "kde"),
    "sgonly/pois": R("sgonly", "poisson"),
    "sgz/pois": R("sgz", "poisson"),
    "sgz/pois.sm2": R("sgz", "poisson", sm=0.2),
    "sgz/pois.sm3": R("sgz", "poisson", sm=0.3),
    "sgz_v50/pois.sm2": R("sgz_v50", "poisson", sm=0.2),
    "sgz_v35/pois.sm2": R("sgz_v35", "poisson", sm=0.2),
    "sgz_v25/pois.sm2": R("sgz_v25", "poisson", sm=0.2),
    "sgz_v15/pois.sm2": R("sgz_v15", "poisson", sm=0.2),
    "sgz_v10/pois.sm2": R("sgz_v10", "poisson", sm=0.2),
    "sgz_v25/pois.sm0": R("sgz_v25", "poisson", sm=0.0),
    "sgz_v25/pois.sm3": R("sgz_v25", "poisson", sm=0.3),
    "sgz_v25/pois.sm2.bg": R("sgz_v25", "poisson", sm=0.2, bg=True),
    "sgz_v25/kde": R("sgz_v25", "kde"),
    "struct/pois.sm2.bg": R("struct_only", "poisson", sm=0.2, bg=True),
}


def fit_recipe(s, rec):
    sel = SUBSETS[rec["sub"]](s)
    bg = rho_bg(s) if rec["bg"] else None
    if rec["fit"] == "poisson":
        return poisson_density(s["g_hat"][sel], s["eff"][sel], h_smooth=rec["sm"], bg=bg)
    return kde_density(s["g_hat"][sel], s["eff"][sel], h=(0.15 if rec["sm"] == 0 else rec["sm"]), bg=bg)


# ---- sweep ----
res = defaultdict(dict)  # recipe -> group -> [emd]
per_group = defaultdict(lambda: defaultdict(list))
for s in scen:
    om = (s["eff"] > 1e-9) & s["isr"]
    oracle = poisson_density(s["G"][om], s["eff"][om])
    for rname, rec in RECIPES.items():
        fit = fit_recipe(s, rec)
        e = float(wasserstein_distance(GRID, GRID, fit, oracle)) if fit is not None else np.nan
        per_group[rname][s["group"]].append(e)

print(f"{'recipe':24s} {'MEAN':>6s} {'MAX':>6s}  {'worst-group':>14s}   (per-group means)")
rows = []
groups = sorted({g for s in scen for g in [s["group"]]})
for rname in RECIPES:
    gm = {g: np.nanmean(per_group[rname][g]) for g in groups if per_group[rname][g]}
    mean = np.nanmean(list(gm.values()))
    worst = max(gm, key=gm.get)
    rows.append((rname, mean, gm[worst], worst, gm))
for rname, mean, mx, worst, gm in sorted(rows, key=lambda r: r[1]):
    tag = f"{worst[0]}/{worst[1].replace('gdna', '')}"
    print(f"{rname:24s} {mean:6.3f} {mx:6.3f}  {tag:>14s}   "
          + " ".join(f"{g[0][-2:]}{g[1].replace('gdna', '')[:3]}:{gm[g]:.2f}" for g in groups))

print("\nRANKED by MEAN:")
for rname, mean, mx, worst, gm in sorted(rows, key=lambda r: r[1])[:5]:
    print(f"  {rname:24s} mean={mean:.3f} max={mx:.3f}")
print("RANKED by MAX (robustness):")
for rname, mean, mx, worst, gm in sorted(rows, key=lambda r: r[2])[:5]:
    print(f"  {rname:24s} max={mx:.3f} mean={mean:.3f}")

# ---- shape plot: the winner recipe vs oracle across the spectrum ----
WINNER = "sgz_v25/pois.sm2"
PANELS = ["gdna_gdna100_ss_0.50_nrna_none_capture_on", "gdna_none_ss_0.50_nrna_none_capture_on",
          "gdna_gdna300_ss_0.99_nrna_none_capture_off", "gdna_gdna1_ss_0.50_nrna_present_capture_on",
          "gdna_none_ss_0.99_nrna_none_capture_off", "gdna_gdna100_ss_0.99_nrna_present_capture_on"]
fig, axes = plt.subplots(2, 3, figsize=(16, 8))
axes = axes.ravel()
by_cond = {s["cond"]: s for s in scen}
for ax, cond in zip(axes, PANELS):
    s = by_cond[cond]
    om = (s["eff"] > 1e-9) & s["isr"]
    oracle = poisson_density(s["G"][om], s["eff"][om])
    fit = fit_recipe(s, RECIPES[WINNER])
    e = float(wasserstein_distance(GRID, GRID, fit, oracle))
    ax.fill_between(GRID, oracle, color="0.82", label="ORACLE")
    ax.plot(GRID, oracle, "k", lw=1.3)
    ax.plot(GRID, fit, "tab:red", lw=1.7, label=f"{WINNER}\nEMD={e:.2f}")
    ax.set_title(cond.replace("gdna_", ""), fontsize=8)
    ax.set_xlim(-5, 2.5)
    ax.set_xlabel("log10 gDNA density", fontsize=8)
    ax.legend(fontsize=7, loc="upper right")
    ax.tick_params(labelsize=7)
fig.suptitle(f"GOAL 1 winner: {WINNER} (confident-25% + zero-count struct, Poisson, smooth) vs ORACLE",
             fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
figp = Path("/Users/mkiyer/proj/rigel/docs/figures/landscape_winner.png")
fig.savefig(figp, dpi=120)
print(f"\nfig -> {figp}")
