"""GOAL 1 explorer: fit the DNA-density landscape as close to the ORACLE as possible, across ALL scenarios,
sweeping TRAINING SUBSETS x FIT METHODS. Metric = EMD(fitted landscape, oracle landscape) in log10 ρ_g space.
The projection (Goal 2) is fixed; an accurate landscape ⇒ accurate projection.

Fit methods:
  kde     — point-estimate unit-mass Gaussian KDE at log10(max(g_hat,1)/E); floors at 1/E (the resolution wall)
            ⇒ CANNOT represent DNA below 1/E (the zero-DNA problem).
  poisson — ADDITIVE Poisson-lognormal: each node deposits its own Poisson likelihood P(g_hat|ρE), unit mass,
            summed (no EM competition). ZERO-NATIVE: g_hat=0 ⇒ e^{−ρE} decays toward ρ→0 ⇒ represents the
            depleted-toward-zero mode. This is the zero-DNA instrumentation.

Training subsets (all EXCLUDE ambiguous both-strand nodes + boundaries):
  single       — single-strand region nodes (fp^fn): inherently solvable (one transcript, one strand).
  single_gonly — single-strand + structural gDNA (|gonly): the v0.7.1 substrate (adds the intergenic/off-target).
  lowvar       — low pass-0 variance nodes (var_gdna ≤ per-scenario median): confidently solved.
  single_lowvar— single-strand AND low variance.
  sgonly_lowvar— (single|gonly) AND low variance.

Oracle landscape = additive Poisson on the TRUE gDNA counts G (zero-native). Prints an EMD table grouped by
capture x DNA-level, and overlays 4 key scenarios."""
from __future__ import annotations

import argparse
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
import dataclasses  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12
LN10 = np.log(10.0)
GRID = np.linspace(-5.0, 2.5, 260)  # log10 ρ_g
H = 0.15  # decades, KDE bandwidth
_lrho_nat = GRID * LN10  # natural-log grid for the Poisson

ap = argparse.ArgumentParser()
ap.add_argument("--outdir", default="/Users/mkiyer/proj/rigel/docs/figures")
a = ap.parse_args()
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index"))
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
cfg = PipelineConfig()
work = Path("/tmp/rigel_selfsolve")
cache = suite / "_selfsolve_cache"
conds = sorted(
    d.name for d in suite.iterdir()
    if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_")
)
PANELS = [
    "gdna_gdna100_ss_0.50_nrna_none_capture_on",
    "gdna_none_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna1_ss_0.50_nrna_present_capture_on",
]


def fit_kde(g_hat, eff):
    """Point-estimate unit-mass KDE, floored at 1/E."""
    if g_hat.size == 0:
        return None
    a_c = np.log10(np.maximum(g_hat, 1.0)) - np.log10(eff)
    z = (GRID[:, None] - a_c[None, :]) / H
    d = (np.exp(-0.5 * z * z) / (H * np.sqrt(2 * np.pi))).sum(1)
    return d / max(float(d.sum()), _EPS)


def fit_poisson(g_hat, eff):
    """Additive Poisson: each node's normalized P(g_hat|ρE) on the grid, summed (unit mass). Zero-native."""
    if g_hat.size == 0:
        return None
    lam = np.exp(_lrho_nat)[None, :] * eff[:, None]  # (n,G) = ρ·E
    ll = g_hat[:, None] * np.log(np.maximum(lam, _EPS)) - lam - gammaln(g_hat[:, None] + 1.0)  # (n,G)
    ll -= ll.max(1, keepdims=True)
    pn = np.exp(ll)
    pn /= np.maximum(pn.sum(1, keepdims=True), _EPS)  # per-node normalized
    d = pn.sum(0)
    return d / max(float(d.sum()), _EPS)


METHODS = {"kde": fit_kde, "poisson": fit_poisson}


def dna_bucket(cond):
    for k in ("none", "gdna1000", "gdna300", "gdna100", "gdna5", "gdna1"):
        if k in cond:
            return k
    return "?"


emd = defaultdict(list)  # (subset, method, group) -> [emd,...]
panel_curves = {}
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, work, cache)
    dbg = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]),
              dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    chain = dbg["chain"]
    cap = dbg["capture"]
    b0 = dbg["belief"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G = Gp + Gn
    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)
    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    isr = np.asarray(chain.kind) == REGION
    idx = np.asarray(chain.ref_idx, np.int64)
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    ntype = np.where(isr, rtype[np.clip(idx, 0, len(rtype) - 1)], -1)
    single = fp ^ fn
    gonly = ~fp & ~fn
    ambig = fp & fn
    live = (eff > 1e-9) & (mass > _EPS)
    f0 = np.asarray(b0.f_g, float)
    var = np.asarray(b0.var_gdna, float)
    g_hat = f0 * mass
    base = live & isr & ~ambig  # never ambiguous
    vmed = np.median(var[base]) if base.any() else 0.0
    lowv = var <= vmed
    # ZERO-DNA instrumentation: include zero-count STRUCTURAL nodes (intergenic/intron, eff>0, mass≈0). Their
    # g_hat=0 ⇒ Poisson(0|ρE) decays toward ρ→0 ⇒ the "representation of what zero DNA looks like" (pooled
    # subtle global DNA). These are exactly the nodes the mass>0 `live` filter drops.
    haveE = eff > 1e-9
    struct_zero = haveE & isr & ((ntype == 0) | (ntype == 1)) & (mass <= _EPS)
    SUBSETS = {
        "single": base & single,
        "single_gonly": base & (single | gonly),
        "lowvar": base & lowv,
        "sgonly_lowvar": base & (single | gonly) & lowv,
        "sgonly_zero": (base & (single | gonly)) | struct_zero,
        "sgonly_lowvar_zero": (base & (single | gonly) & lowv) | struct_zero,
    }
    # ORACLE landscape: additive Poisson on the TRUE counts G over all live region nodes (zero-native)
    om = live & isr
    oracle = fit_poisson(G[om], eff[om])
    grp = ("capON" if "capture_on" in cond else "capOFF", dna_bucket(cond))
    for sname, sel in SUBSETS.items():
        for mname, fn_ in METHODS.items():
            fit = fn_(g_hat[sel], eff[sel])
            e = float(wasserstein_distance(GRID, GRID, fit, oracle)) if fit is not None else np.nan
            emd[(sname, mname, grp)].append(e)
    if cond in PANELS:
        sg = base & (single | gonly)
        sgz = sg | struct_zero
        panel_curves[cond] = dict(
            oracle=oracle,
            kde=fit_kde(g_hat[sg], eff[sg]),
            poisson=fit_poisson(g_hat[sg], eff[sg]),
            poisson_zero=fit_poisson(g_hat[sgz], eff[sgz]),
        )

# ---- table: mean EMD per (subset, method), grouped ----
groups = sorted({g for (_, _, g) in emd})
subsets = ["single", "single_gonly", "lowvar", "sgonly_lowvar", "sgonly_zero", "sgonly_lowvar_zero"]
print(f"\n{'EMD-to-oracle (lower=better)':30s}", end="")
for g in groups:
    print(f"{g[0][-2:]+'/'+g[1].replace('gdna',''):>11s}", end="")
print(f"{'MEAN':>8s}")
best = {}
for sname in subsets:
    for mname in METHODS:
        row = []
        for g in groups:
            vals = emd.get((sname, mname, g), [])
            row.append(np.nanmean(vals) if vals else np.nan)
        m = np.nanmean(row)
        best[(sname, mname)] = m
        print(f"{sname+'/'+mname:30s}", end="")
        for v in row:
            print(f"{v:11.2f}", end="")
        print(f"{m:8.2f}")

print("\nRANKED by mean EMD:")
for (s, m), v in sorted(best.items(), key=lambda kv: kv[1])[:6]:
    print(f"  {s}/{m:10s}  mean EMD={v:.3f}")

# ---- panel plot ----
outdir = Path(a.outdir)
fig, axes = plt.subplots(2, 2, figsize=(13, 8))
axes = axes.ravel()
for ax, cond in zip(axes, PANELS):
    pc = panel_curves[cond]
    ax.fill_between(GRID, pc["oracle"], color="0.82", label="ORACLE (Poisson on G)")
    ax.plot(GRID, pc["oracle"], "k", lw=1.3)
    if pc["kde"] is not None:
        ax.plot(GRID, pc["kde"], "tab:blue", lw=1.4, ls="--", label="KDE (single|gonly)")
    if pc["poisson"] is not None:
        ax.plot(GRID, pc["poisson"], "tab:red", lw=1.6, label="Poisson (single|gonly)")
    if pc["poisson_zero"] is not None:
        ax.plot(GRID, pc["poisson_zero"], "tab:green", lw=1.5, label="Poisson (+ zero-count struct)")
    ax.set_title(cond.replace("gdna_", ""), fontsize=8)
    ax.set_xlim(-5, 2.5)
    ax.set_xlabel("log10 gDNA density", fontsize=8)
    ax.legend(fontsize=6, loc="upper right")
    ax.tick_params(labelsize=7)
fig.suptitle("GOAL 1: fitted DNA landscape vs ORACLE (Poisson-on-G). KDE floors at 1/E; Poisson is zero-native.",
             fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.97))
figp = outdir / "landscape_explore.png"
fig.savefig(figp, dpi=120)
print(f"\nfig -> {figp}")
