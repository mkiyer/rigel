"""Does a boundary RECEIVE the information it needs to solve? (owner's thesis: it does, but the arithmetic
produces a bad solve.) Decisive test, per unstranded boundary POSITION across the 20 ss0.50 scenarios:

    corr_solved  = corr(solved f_g, oracle)            — the current solve quality (~0.04, coin-flip)
    corr_msg     = corr(exp(mode_g_combined), oracle)  — what the COMBINED incoming gDNA message implies
    corr_msgL/R  = corr(exp(mode_g_left/right), oracle) — per-flank (fwd scan = left msg, bwd = right msg)

The message mode is log(f_g^dst) (both _mode_shift and _mode_density return log f_g in the dst frame), so
exp(mode_g) is the f_g the message tells the boundary to adopt. On unstranded data the final solve is
message-dominated (strand≈0, prior weak), so:
  * corr_msg ≫ corr_solved  ⇒ the info IS there, the FOLD/COMBINE loses it (owner right → fold bug).
  * corr_msg ≈ corr_solved ≈ 0  ⇒ the message CONTENT is wrong (mode/source bug) or genuinely uninformative.
A per-flank split shows whether ONE flank (e.g. the enriched exon) carries the answer while the combine dilutes it."""
from __future__ import annotations
import dataclasses
import importlib
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth
from flagship_interrogate import _oracle_per_node
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.bp_solver import REGION
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-9
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
_SS = sys.argv[1] if len(sys.argv) > 1 else "0.50"
_CAP = sys.argv[2] if len(sys.argv) > 2 else ""  # "" all | "capture_off" | "capture_on" | "capture_verystrong"
def _capmatch(name):
    if not _CAP:
        return True
    if _CAP == "capture_on":  # exact (not verystrong)
        return "capture_on" in name and "verystrong" not in name
    return _CAP in name
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_")
               and _SS in d.name and _capmatch(d.name))


def _impl(mode, prec):
    """implied f_g from a message = exp(mode) clipped; NaN where no message (prec≈0) so it's excluded."""
    f = np.exp(np.clip(mode, -30, 0.0))  # mode = log f_g ≤ 0
    return np.where(prec > _EPS, np.clip(f, 0.0, 1.0), np.nan)


acc = {k: [] for k in ("nid", "fo", "fg", "mass", "msg", "msgL", "msgR", "pg", "spl", "cliff")}
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache")
    dbg = {}; cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    chain, cap = dbg["chain"], dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    fg = np.asarray(cap["f_g"]); mass = np.asarray(cap["mass_global"]); kind = np.asarray(chain.kind)
    mode_g, prec_g = np.asarray(cap["mode_g"]), np.asarray(cap["prec_g"])  # COMBINED gDNA message
    a_fwd, b_bwd = cap["a_fwd"], cap["b_bwd"]
    msg = _impl(mode_g, prec_g)
    msgL = _impl(np.asarray(a_fwd[0]), np.asarray(a_fwd[1]))  # left-flank gDNA message
    msgR = _impl(np.asarray(b_bwd[0]), np.asarray(b_bwd[1]))  # right-flank gDNA message
    spl = np.asarray(cap["spl_l"]) + np.asarray(cap["spl_r"])
    mu_proj = np.log(np.maximum(np.asarray(cap["_uni_static"]["rho_node0"]), 1e-12)); left, right = np.asarray(chain.left), np.asarray(chain.right)
    cliff = np.array([abs((mu_proj[left[i]] if left[i] >= 0 else mu_proj[i]) -
                          (mu_proj[right[i]] if right[i] >= 0 else mu_proj[i])) for i in range(len(fg))])
    bnd = (kind != REGION) & np.isfinite(fo) & (mass > _EPS)
    for i in np.where(bnd)[0]:
        acc["nid"].append(int(i)); acc["fo"].append(fo[i]); acc["fg"].append(fg[i]); acc["mass"].append(mass[i])
        acc["msg"].append(msg[i]); acc["msgL"].append(msgL[i]); acc["msgR"].append(msgR[i])
        acc["pg"].append(prec_g[i]); acc["spl"].append(spl[i]); acc["cliff"].append(cliff[i])
A = {k: np.array(v) for k, v in acc.items()}


def _pos_corr(nid, x, y, minobs=6):
    """per-position corr(x,y) across scenarios (NaN-aware), returned per unique nid (NaN if too few / no var)."""
    out = {}
    for u in np.unique(nid):
        m = nid == u
        xs, ys = x[m], y[m]
        ok = np.isfinite(xs) & np.isfinite(ys)
        if ok.sum() < minobs or xs[ok].std() < 1e-6 or ys[ok].std() < 1e-6:
            out[u] = np.nan
        else:
            out[u] = np.corrcoef(xs[ok], ys[ok])[0, 1]
    return out


_MINOBS = 4 if _CAP else 6  # a single capture regime has ~5 gDNA levels
cs = _pos_corr(A["nid"], A["fg"], A["fo"], _MINOBS)      # solved vs oracle
cm = _pos_corr(A["nid"], A["msg"], A["fo"], _MINOBS)     # combined message vs oracle
cl = _pos_corr(A["nid"], A["msgL"], A["fo"], _MINOBS)    # left-flank message vs oracle
cr = _pos_corr(A["nid"], A["msgR"], A["fo"], _MINOBS)    # right-flank message vs oracle
uniq = np.unique(A["nid"])
massp = {u: A["mass"][A["nid"] == u].mean() for u in uniq}


def _mean(d, keys):
    v = np.array([d[u] for u in keys]); v = v[np.isfinite(v)]
    return float(v.mean()) if len(v) else float("nan"), len(v)


print(f"[ss_{_SS}] {len(conds)} scenarios, {len(uniq)} boundary positions\n")
print("Per-position corr with oracle (hold boundary fixed, vary gDNA level):")
for name, d in (("solved f_g", cs), ("COMBINED msg (exp mode_g)", cm),
                ("left-flank msg", cl), ("right-flank msg", cr)):
    mval, n = _mean(d, uniq)
    print(f"  {name:>28} : mean corr_i = {mval:>7.3f}  (n_pos={n})")

# the crux: on positions where BOTH solved and message corr are defined, does the message beat the solve?
best = {u: np.nanmax([cl.get(u, np.nan), cr.get(u, np.nan)]) for u in uniq}  # best single flank
both_def = [u for u in uniq if np.isfinite(cs[u]) and np.isfinite(cm[u])]
mb_msg = np.array([cm[u] for u in both_def]); mb_solv = np.array([cs[u] for u in both_def])
print(f"\nOn {len(both_def)} positions where both are defined:")
print(f"  mean corr(message) = {mb_msg.mean():.3f}   mean corr(solved) = {mb_solv.mean():.3f}   "
      f"Δ = {mb_msg.mean() - mb_solv.mean():+.3f}")
print(f"  positions where message beats solve by >0.1: {np.mean(mb_msg - mb_solv > 0.1):.1%}")
bestv = np.array([best[u] for u in both_def]); bestv = bestv[np.isfinite(bestv)]
print(f"  mean BEST-single-flank corr = {bestv.mean():.3f}  (does one flank carry more than the combine?)")

# split the message-carries-info test by spliced presence (the CUT-6 axis)
spl_pos = {u: (A["spl"][A["nid"] == u] > _EPS).mean() > 0.5 for u in uniq}
for lab, sel in (("spliced boundaries", True), ("no-spliced boundaries", False)):
    ks = [u for u in both_def if spl_pos[u] == sel]
    if len(ks) < 5:
        continue
    mm = np.array([cm[u] for u in ks]); ss = np.array([cs[u] for u in ks])
    print(f"  [{lab:>22}] corr(msg)={mm.mean():.3f}  corr(solved)={ss.mean():.3f}  Δ={mm.mean()-ss.mean():+.3f}  (n={len(ks)})")
