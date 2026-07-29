"""MECHANISTIC DISSECTION of the node-1055 f_g crush. Answers the owner's question directly:
'if the prior has a mode at the enriched level and the node's density is there, why does f_g go to 0 not 5?'

Strategy: hold node 1055's REAL data fixed; swap ONLY the prior SHAPE; watch the f_g read-out.
  Part A  reproduce the real end-to-end crush + dump node facts + the REAL fitted prior's modes.
  Part B  ISOLATED delta-pin (flat/unstranded strand => psi = logP_g(node) + logP_r_ref), 4 prior shapes:
            (a) the REAL fitted prior            (b) synthetic DEPLETED-only mode
            (c) synthetic ENRICHED mode at TRUTH (d) synthetic BIMODAL (depleted + enriched)
          For each: logP_g sampled across f_g, and the median / mean / argmax read-out.
  Part C  FULL real refit with the synthetic ENRICHED prior swapped in (monkeypatch the fit) — the faithful
          confirmation: does node 1055 then land on the enriched mode?

If (c)/(C) pin node 1055 to ~0.9 => NO delta-pin/solver bug; the crush is the FITTED PRIOR lacking an
enriched mode (a self-fulfilling under-call). If they still go to 0 => a real solver/read-out bug."""
from __future__ import annotations

import dataclasses
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402
from flagship_interrogate import _oracle_per_node  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402

calmod = sys.modules["rigel.calibration.calibrate"]  # the MODULE (import-as would bind the function)
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.calibration.simplex_logodds import (  # noqa: E402
    _DEFAULT_L, _gdna_arm, _logodds_grid, _posterior_median_fg, _rna_arm,
)
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12
LN10 = np.log(10.0)
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index"))
ra = RegionArrays.from_index(index)
cfg = PipelineConfig()
work = Path("/tmp/rigel_selfsolve")
cache = suite / "_selfsolve_cache"
COND = "gdna_gdna100_ss_0.50_nrna_none_capture_on"
inp = _scan_and_truth(suite, COND, index, cfg, work, cache)


def run(iters):
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=iters)
    dbg = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    return dbg


# ---------- Part A: real crush + node facts + real prior ----------
d0 = run(0)
d1 = run(1)
chain = d0["chain"]
cap = d0["capture"]
hp = d1["gdna_hyperprior"]
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
f0 = np.asarray(d0["belief"].f_g, float)
f1 = np.asarray(d1["belief"].f_g, float)
fo = np.divide(G, np.maximum(G + Rt, _EPS))
sel = isr & (ntype == 2) & single & (eff > 1e-9) & (mass > _EPS) & (fo > 0.6)
i = int(np.argmax(G * sel))  # full-array index of the max-gDNA crush-pole node
M, E = float(mass[i]), float(eff[i])
rho_tot = M / E
log_me = np.log(M) - np.log(E)  # natural log(M/E)
true_lrg10 = np.log10(max(fo[i] * rho_tot, 1e-12))

print("=" * 92)
print(f"NODE {i}  ({COND})   — captured single-strand exon, UNSTRANDED library (ss0.50)")
print("=" * 92)
print(f"  mass M = {M:.0f}   eff E = {E:.0f}   ->  total density rho_tot = M/E = {rho_tot:.2f} "
      f"(log10 {np.log10(rho_tot):+.2f})")
print(f"  ORACLE: G={G[i]:.0f} gDNA, R={Rt[i]:.0f} RNA  ->  true f_g = {fo[i]:.3f}   "
      f"(true rho_g log10 {true_lrg10:+.2f})")
print(f"  pass-0 f_g = {f0[i]:.3f}  (var {np.asarray(d0['belief'].var_gdna, float)[i]:.2f}) "
      f"->  REFIT f_g = {f1[i]:.4f}   <<< the crush")
lr10 = hp.log_rho / LN10
top = np.argsort(-np.asarray(hp.weights, float))[:6]
print("  REAL fitted prior modes (log10 rho_g, weight): "
      + "  ".join(f"[{lr10[j]:+.2f},w{np.asarray(hp.weights, float)[j]:.2f}]"
                  for j in sorted(top, key=lambda k: lr10[k])))
print(f"  fitted prior grid spans log10 rho_g in [{lr10.min():+.2f}, {lr10.max():+.2f}]; "
      f"node needs mass near {true_lrg10:+.2f} to pin the truth "
      f"(accessible max at f_g=1 is {np.log10(rho_tot):+.2f}).")

# ---------- Part B: isolated delta-pin, 4 prior shapes ----------
n_grid = int(cfg.calibration.sweep_n_grid_single_strand)
lam, fg = _logodds_grid(n_grid, _DEFAULT_L)


def _lse(a):
    m = np.max(a)
    return m + np.log(np.sum(np.exp(a - m)))


def gauss_logP(grid, centers_log10, weights, h_dec=0.15):
    """A normalized mixture log-density on the natural-log grid, modes at centers_log10 (decades)."""
    h = h_dec * LN10
    dens = np.zeros_like(grid)
    for c, w in zip(centers_log10, weights):
        z = (grid - c * LN10) / h
        dens += w * np.exp(-0.5 * z * z) / (h * np.sqrt(2 * np.pi))
    dens /= max(dens.sum(), _EPS)
    return np.log(np.maximum(dens, 1e-30))


# a fresh wide grid so the enriched mode is representable (the real prior's grid may not reach it)
gridN = np.linspace(-6 * LN10, 3 * LN10, 700)


def dpin_isolated(logP, log_rho, label):
    """psi = logP_g(node)  +  RNA Jeffreys reference (flat/unstranded strand => strand term is constant)."""
    lrg = np.log(fg) + log_me  # (K,) natural log rho_g(f)
    lp_g = np.interp(lrg, log_rho, logP, left=logP[0], right=logP[-1])  # the gDNA arm
    psi = _gdna_arm(lam, lp_g[None, :])[0] + _rna_arm(lam)[0]  # + flat strand (const)
    post = np.exp(psi - _lse(psi))
    med = float(_posterior_median_fg(post[None, :], fg)[0])
    mean = float(post @ fg)
    amax = float(fg[int(np.argmax(psi))])
    # sample logP_g across f_g so the owner/reviewer can read the pin directly
    samp = {}
    for fq in (0.001, 0.01, 0.1, 0.375, 0.5, 0.9):
        samp[fq] = float(np.interp(np.log(fq) + log_me, log_rho, logP, left=logP[0], right=logP[-1]))
    return med, mean, amax, samp


priors = [
    ("(a) REAL fitted prior", hp.logP, hp.log_rho),
    ("(b) synth DEPLETED-only @ -2.0", gauss_logP(gridN, [-2.0], [1.0]), gridN),
    (f"(c) synth ENRICHED @ truth {true_lrg10:+.2f}", gauss_logP(gridN, [true_lrg10], [1.0]), gridN),
    (f"(d) synth BIMODAL (-2.0 & {true_lrg10:+.2f})", gauss_logP(gridN, [-2.0, true_lrg10], [0.5, 0.5]), gridN),
]
print("\n" + "-" * 92)
print("PART B — ISOLATED delta-pin (flat unstranded strand; psi = logP_g(node) + RNA ref). "
      "read-out of f_g:")
print("-" * 92)
print(f"  {'prior shape':46s} {'median':>7s} {'mean':>7s} {'argmax':>7s}   logP_g at f_g=[.001 .01 .1 .375 .5 .9]")
for label, logP, log_rho in priors:
    med, mean, amax, samp = dpin_isolated(logP, log_rho, label)
    s = " ".join(f"{samp[q]:+5.1f}" for q in (0.001, 0.01, 0.1, 0.375, 0.5, 0.9))
    print(f"  {label:46s} {med:7.3f} {mean:7.3f} {amax:7.3f}   [{s}]")

# ---------- Part C: FULL real refit with the synthetic ENRICHED prior swapped in ----------
print("\n" + "-" * 92)
print("PART C — FULL real refit, but the fitted prior REPLACED by the synthetic ENRICHED mode "
      f"@ {true_lrg10:+.2f} (all else identical). Faithful confirmation:")
print("-" * 92)
enr_logP = gauss_logP(gridN, [true_lrg10], [1.0])
synth = dataclasses.replace(hp, log_rho=gridN.copy(), logP=enr_logP.copy(),
                            weights=np.exp(enr_logP - enr_logP.max()))
_orig_fit = calmod._fit_gdna_hyperprior


def _patched_fit(*a, **k):
    return synth


calmod._fit_gdna_hyperprior = _patched_fit
try:
    dC = run(1)
finally:
    calmod._fit_gdna_hyperprior = _orig_fit
fC = float(np.asarray(dC["belief"].f_g, float)[i])
print(f"  node {i}: true f_g={fo[i]:.3f} | real-prior refit={f1[i]:.4f} | "
      f"ENRICHED-prior refit={fC:.4f}   "
      f"({'PINS the enriched mode -> NO solver bug' if fC > 0.5 else 'STILL crushes -> SOLVER BUG'})")
print("=" * 92)
