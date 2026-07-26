"""RC2 (B) — EXACT-REPLAY per-channel ψ ablation for single-strand (full-rank) nodes.

`scripts/debug/pass0_node_dissect.py` no longer replays the shipped solve faithfully: `_uni_msg` publishes
the MODE-FUSION precisions (`cpg/cpp/cpn`) while ψ is actually handed the MEASUREMENT precisions
(`cm_g/cm_p/cm_n`), and the λ / θ messages are not replayed at all (measured fidelity 6.0e-01). This script
rebuilds ψ term-by-term from `_uni[-1]` (the arrays the shipped final `_local_solve` really received) and
validates against the shipped `f_g` before ablating anything.

ψ (single-strand, 1-D on λ) =
      strand BB  +  gDNA arm (½·log f_g, or the fitted global prior)  +  RNA arm (½·log(1−f_g))
    + intron-factory λ prior  +  gdna_imp  +  rna_imp₊  +  rna_imp₋  +  λ-imp

    OMP_NUM_THREADS=1 python scratchpad/rc2_b_ablate.py [COND] [--nodes 2197,...] [--top 6]
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.simplex import _mixture_strand_loglik  # noqa: E402
from rigel.calibration.simplex_logodds import (  # noqa: E402
    _log1m_fg,
    _log_fg,
    _logodds_grid,
    _posterior_median_fg,
    _regrid_global,
    _lse,
)
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
_EPS = 1e-9

ap = argparse.ArgumentParser()
ap.add_argument("cond", nargs="?", default="gdna_gdna300_ss_0.99_nrna_none_capture_on")
ap.add_argument("--nodes", default=None)
ap.add_argument("--top", type=int, default=6)
ap.add_argument("--refit", type=int, default=0)
a = ap.parse_args()

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, a.cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
dbg: dict = {}
cc = dataclasses.replace(cfg.calibration, calib_refit_iters=a.refit)
calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                 np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
uni, us = cap["_uni"][-1], cap["_uni_static"]
prs = dbg["calibration_priors"]

Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
G, R = Gp + Gn, Rp + Rn
fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
fg_ship = np.asarray(cap["f_g"])
self_fg = np.asarray(cap["fg_loc"])
mass = np.asarray(cap["mass_global"])
tau0 = np.asarray(cap["_tau0_lam"])
rt, _ = _node_region_type(chain, ra)
kind = np.asarray(chain.kind)
cls = np.array([CLS[int(rt[i])] if kind[i] == REGION else "boundary" for i in range(fg_ship.size)])
fp_ = np.asarray(st.free_pos, bool)
fn_ = np.asarray(st.free_neg, bool)
solvable = np.asarray(cap["solvable_mask"], bool)
single = fp_ ^ fn_
is_amb = fp_ & fn_

# ── the arrays the shipped final `_local_solve` actually received ─────────────────────────────────────────
MO_G, CM_G = np.asarray(uni["mo_g"]), np.asarray(uni["cm_g"])
MO_P, CM_P = np.asarray(uni["mo_p"]), np.asarray(uni["cm_p"])
MO_N, CM_N = np.asarray(uni["mo_n"]), np.asarray(uni["cm_n"])
LAM_M, LAM_P = np.asarray(uni["lam_msg"]), np.asarray(uni["c_tau"])
cp, cn = np.asarray(uni["cp"]), np.asarray(uni["cn"])
cR = cp + cn
tau_tilt = np.clip(np.where(cR > _EPS, (cp - cn) / np.maximum(cR, _EPS), 0.0), -1.0, 1.0)
TH_M = np.arcsin(tau_tilt)
TH_P = np.where(is_amb, CM_P + CM_N, 0.0)
GLP, IP = cap["global_lp"], cap["intron_prior"]
FGR, FPR, FNR = cap["fg_init"], cap["fpos_init"], cap["fneg_init"]
KAP = float(dbg["rna_sense_frac"])
OD_G, OD_R = float(prs.gdna_strand_overdispersion), float(prs.rna_strand_overdispersion)

K_SS = int(cc.sweep_n_grid_single_strand or cc.sweep_n_grid)
L = float(cc.sweep_logodds_window)
lam, fgrid = _logodds_grid(K_SS, L)
LOG_FG, LOG_FA = _log_fg(lam), _log1m_fg(lam)
GLP_SS = _regrid_global(GLP, cc.sweep_n_grid, K_SS, L) if GLP is not None else None
IP_SS = _regrid_global(IP, cc.sweep_n_grid, K_SS, L) if IP is not None else None

CHANNELS = ("strand", "ref_g", "ref_r", "factory", "gdna_imp", "rna_imp_p", "rna_imp_n", "lam_imp")


def psi_terms(i: int) -> dict[str, np.ndarray]:
    """Every ψ term of node ``i`` on the fine λ grid, as (K,) arrays. Single-strand only."""
    up, un = float(st.u_pos[i]), float(st.u_neg[i])
    n = up + un
    pos_live = bool(fp_[i]) and not bool(fn_[i])
    neg_live = bool(fn_[i]) and not bool(fp_[i])
    f_act = 1.0 - fgrid
    f_pos = f_act if pos_live else np.zeros_like(f_act)
    f_neg = f_act if neg_live else np.zeros_like(f_act)
    t = {}
    t["strand"] = _mixture_strand_loglik(
        np.array([up])[:, None], np.array([n])[:, None], fgrid[None, :],
        f_pos[None, :], f_neg[None, :], KAP, OD_G, OD_R,
        np.array([FGR[i]])[:, None], np.array([FPR[i]])[:, None], np.array([FNR[i]])[:, None],
    )[0]
    t["ref_g"] = 0.5 * LOG_FG if GLP_SS is None else np.asarray(GLP_SS[i], np.float64)
    t["ref_r"] = 0.5 * LOG_FA
    t["factory"] = np.zeros(K_SS) if IP_SS is None else np.asarray(IP_SS[i], np.float64)
    t["gdna_imp"] = -0.5 * CM_G[i] * (LOG_FG - MO_G[i]) ** 2
    t["rna_imp_p"] = -0.5 * CM_P[i] * (LOG_FA - MO_P[i]) ** 2
    t["rna_imp_n"] = -0.5 * CM_N[i] * (LOG_FA - MO_N[i]) ** 2
    t["lam_imp"] = -0.5 * LAM_P[i] * (lam - LAM_M[i]) ** 2
    return t


def readout(psi: np.ndarray) -> float:
    post = np.exp(psi - _lse(psi[None, :], axis=1, keepdims=True))
    return float(_posterior_median_fg(post, fgrid)[0])


def solve_with(i: int, drop=()) -> float:
    t = psi_terms(i)
    psi = np.zeros(K_SS)
    for k, v in t.items():
        if k not in drop:
            psi = psi + v
    return readout(psi)


# ── select target nodes ───────────────────────────────────────────────────────────────────────────────────
err = np.abs(fg_ship - fo)
serr = np.abs(self_fg - fo)
ok = np.isfinite(fo) & (mass > _EPS) & solvable & single & (tau0 > 0.0)
if a.nodes:
    targets = [int(x) for x in a.nodes.split(",")]
else:
    hurt = ok & (err > serr + 0.02) & (cls == "exon")
    targets = np.argsort(np.where(hurt, -(err * mass), 0.0))[: a.top].astype(int).tolist()

# ── replay fidelity over ALL single-strand nodes ──────────────────────────────────────────────────────────
sub = np.flatnonzero(ok)
rep = np.array([solve_with(int(i)) for i in sub])
fid = np.max(np.abs(rep - fg_ship[sub]))
print(f"# cond={a.cond}  refit={a.refit}   kappa={KAP:.4f} od_g={OD_G:.3g} od_r={OD_R:.3g}")
print(f"# REPLAY FIDELITY over {sub.size:,} full-rank single-strand nodes: "
      f"max|replay-shipped| = {fid:.3e}  {'OK' if fid < 1e-9 else '*** MISMATCH ***'}")

MSGS = ("gdna_imp", "rna_imp_p", "rna_imp_n", "lam_imp")
INTRINSIC = ("strand", "ref_g", "ref_r", "factory")


def solve_keep(i: int, keep_msgs=()) -> float:
    return solve_with(i, drop=tuple(m for m in MSGS if m not in keep_msgs))


print(f"\n{'=' * 132}\nPER-CHANNEL ABLATION  (each column = that ONE channel REMOVED; 'ALL' = shipped)\n{'=' * 132}")
hdr = f"{'node':>6}{'cls':>9}{'mass':>11}{'oracle':>8}{'self':>8}{'ALL':>8}"
hdr += "".join(f"{'-' + c:>11}" for c in CHANNELS)
print(hdr)
for i in targets:
    i = int(i)
    base = solve_with(i)
    row = f"{i:>6}{cls[i]:>9}{mass[i]:>11,.0f}{fo[i]:>8.3f}{self_fg[i]:>8.3f}{base:>8.3f}"
    row += "".join(f"{solve_with(i, drop=(c,)):>11.3f}" for c in CHANNELS)
    print(row)

print(f"\n{'=' * 132}\nKEEP-ONLY  (intrinsic ψ = strand+ref_g+ref_r+factory, then ONE message channel added)"
      f"\n{'=' * 132}")
print(f"{'node':>6}{'oracle':>8}{'ALL':>8}{'intrinsic':>11}{'+gdna':>9}{'+rna+':>9}{'+rna-':>9}"
      f"{'+lam':>9}{'   -[rna-,lam]':>15}{'  -[rna+,rna-,lam]':>19}")
for i in targets:
    i = int(i)
    print(f"{i:>6}{fo[i]:>8.3f}{solve_with(i):>8.3f}{solve_keep(i):>11.3f}"
          f"{solve_keep(i, ('gdna_imp',)):>9.3f}{solve_keep(i, ('rna_imp_p',)):>9.3f}"
          f"{solve_keep(i, ('rna_imp_n',)):>9.3f}{solve_keep(i, ('lam_imp',)):>9.3f}"
          f"{solve_with(i, drop=('rna_imp_n', 'lam_imp')):>15.3f}"
          f"{solve_with(i, drop=('rna_imp_p', 'rna_imp_n', 'lam_imp')):>19.3f}")

# ── the energy budget: each term's value at the oracle f_g vs at the shipped f_g ──────────────────────────
print(f"\n{'=' * 130}\nENERGY BUDGET — ψ_term(oracle f_g) − ψ_term(shipped f_g)  [nats; NEGATIVE = the term "
      f"OPPOSES the oracle]\n{'=' * 130}")
for i in targets:
    i = int(i)
    t = psi_terms(i)
    jo = int(np.argmin(np.abs(fgrid - fo[i])))
    js = int(np.argmin(np.abs(fgrid - fg_ship[i])))
    tot = 0.0
    parts = []
    for c in CHANNELS:
        d_ = float(t[c][jo] - t[c][js])
        tot += d_
        parts.append(f"{c}={d_:+.1f}")
    print(f"  node {i:<6} [{cls[i]}] oracle={fo[i]:.3f}(grid {fgrid[jo]:.3f}) shipped={fg_ship[i]:.3f}"
          f"(grid {fgrid[js]:.3f})   TOTAL={tot:+.1f}")
    print("        " + "  ".join(parts))
    print(f"        MSGS: gdna mo={np.exp(MO_G[i]):.4f} p={CM_G[i]:.4g} | rna+ mo={np.exp(MO_P[i]):.4f} "
          f"p={CM_P[i]:.4g} | rna- mo={np.exp(MO_N[i]):.4f} p={CM_N[i]:.4g} | "
          f"lam f_g={1 / (1 + np.exp(-LAM_M[i])):.4f} p={LAM_P[i]:.4g} | tau_own={tau0[i]:.4g}")
    print(f"        ORACLE gDNA={G[i]:,.0f} RNA+={Rp[i]:,.0f} RNA-={Rn[i]:,.0f}  "
          f"u_pos={st.u_pos[i]:,.0f} u_neg={st.u_neg[i]:,.0f} free=({int(fp_[i])},{int(fn_[i])})")

# ── AGGREGATE: is ψ's Jeffreys RNA reference what holds the f_g→1 vertex bin off the truth? ───────────────
print(f"\n{'=' * 132}\nAGGREGATE over ALL full-rank single-strand nodes, by ORACLE bin "
      f"(mass-weighted mean solved f_g)\n{'=' * 132}")
VAR = {
    "ALL(shipped)": lambda i: solve_with(i),
    "no ref_r": lambda i: solve_with(i, drop=("ref_r",)),
    "no ref_g+ref_r": lambda i: solve_with(i, drop=("ref_g", "ref_r")),
    "intrinsic only": lambda i: solve_keep(i),
    "intrinsic-ref_r": lambda i: solve_with(i, drop=("ref_r",) + MSGS),
    "no rna_imp": lambda i: solve_with(i, drop=("rna_imp_p", "rna_imp_n")),
    "no lam_imp": lambda i: solve_with(i, drop=("lam_imp",)),
    "no rna+lam": lambda i: solve_with(i, drop=("rna_imp_p", "rna_imp_n", "lam_imp")),
    "no gdna_imp": lambda i: solve_with(i, drop=("gdna_imp",)),
}
res = {k: np.array([f(int(i)) for i in sub]) for k, f in VAR.items()}
BINS = [(0.0, 0.30), (0.30, 0.60), (0.60, 0.90), (0.90, 0.99), (0.99, 1.01)]
print(f"{'oracle bin':<13}{'n':>6}{'mass':>12}{'oracle':>8}{'self':>8}" + "".join(f"{k:>16}" for k in VAR))
for lo, hi in BINS:
    m = (fo[sub] >= lo) & (fo[sub] < hi)
    if not m.any():
        continue
    w = mass[sub][m]
    row = "".join(f"{np.average(res[k][m], weights=w):>16.3f}" for k in VAR)
    print(f"[{lo:.2f},{hi:.2f})  {int(m.sum()):>6}{w.sum():>12,.0f}"
          f"{np.average(fo[sub][m], weights=w):>8.3f}{np.average(self_fg[sub][m], weights=w):>8.3f}{row}")
print(f"\n{'mwae':<13}{'':>6}{'':>12}{'':>8}{np.average(np.abs(self_fg[sub] - fo[sub]), weights=mass[sub]):>8.4f}"
      + "".join(f"{np.average(np.abs(res[k] - fo[sub]), weights=mass[sub]):>16.4f}" for k in VAR))
