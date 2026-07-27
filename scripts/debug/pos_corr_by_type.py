"""VALIDATE the honest per-position metric against ALL node types. The pooled corr(f_g, oracle) in
`docs/calibration/archive/solve_gate_design.md` says single-strand + AMBIG REGIONS resolve (0.63/0.69) while boundaries are coin-flips
(0.13). But `boundary_dissect.py` showed the pooled corr is Simpson-confounded by between-structure signal — and
the honest per-position corr (hold a node FIXED, vary the gDNA level across the 20 unstranded scenarios, ask if
its solve tracks its oracle) collapses EVERY boundary stratum to ~0.

The decisive question this script answers: does the per-position metric ALSO collapse for regions? If regions
RESOLVE per-position (high corr_i) but boundaries do NOT, the metric is validated and the finding is clean —
boundaries are unidentifiable on unstranded data, regions are not. If regions ALSO collapse per-position, then
the "DOF valid for regions" conclusion rests on the same confound and must be revisited."""
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
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
_SS = sys.argv[1] if len(sys.argv) > 1 else "0.50"  # "0.50" unstranded (default) | "0.99" stranded
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_") and _SS in d.name)
print(f"[strand-specificity filter: ss_{_SS}]")

# per (nid, scenario): fg, fo, mass, type, dof — accumulate then compute per-position corr across scenarios.
acc = {"nid": [], "fg": [], "fo": [], "mass": [], "typ": [], "dof": []}
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache")
    dbg = {}; cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    chain, cap = dbg["chain"], dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    fg = np.asarray(cap["f_g"]); mass = np.asarray(cap["mass_global"])
    fp = np.asarray(cap["free_pos"], bool); fn = np.asarray(cap["free_neg"], bool)
    pg, pp, pn = np.asarray(cap["prec_g"]), np.asarray(cap["prec_p"]), np.asarray(cap["prec_n"])
    t0 = np.asarray(cap["_tau0_lam"]); kind = np.asarray(chain.kind)
    lam_id = (t0 > _EPS) | (pg > _EPS) | (pp + pn > _EPS)
    th_id = (t0 > _EPS) | (pp > _EPS) | (pn > _EPS)
    nfree = fp.astype(int) + fn.astype(int)
    dof_solv = np.where(nfree >= 2, lam_id & th_id, np.where(nfree == 1, lam_id, False))
    typ = np.where(kind != REGION, 3, np.where(nfree >= 2, 2, np.where(nfree == 1, 1, 0)))  # 0=intgc 1=single 2=AMBIG 3=bnd
    ok = np.isfinite(fo) & (mass > _EPS)
    for i in np.where(ok)[0]:
        acc["nid"].append(int(i)); acc["fg"].append(fg[i]); acc["fo"].append(fo[i])
        acc["mass"].append(mass[i]); acc["typ"].append(int(typ[i])); acc["dof"].append(bool(dof_solv[i]))

A = {k: np.array(v) for k, v in acc.items()}
_MINOBS = 6
_TN = {0: "intergenic", 1: "single-reg", 2: "AMBIG-reg", 3: "boundary"}

# per-position corr across scenarios
uniq = np.unique(A["nid"])
rows = []  # (nid, corr_i, mass, typ, dof)
for u in uniq:
    m = A["nid"] == u
    fog, fgg = A["fo"][m], A["fg"][m]
    if m.sum() < _MINOBS or fog.std() < 1e-6 or fgg.std() < 1e-6:
        continue
    ci = np.corrcoef(fgg, fog)[0, 1]
    rows.append((u, ci, A["mass"][m].mean(), int(np.round(np.median(A["typ"][m]))),
                 bool(A["dof"][m].mean() > 0.5)))
corr = np.array([r[1] for r in rows]); pmass = np.array([r[2] for r in rows])
ptyp = np.array([r[3] for r in rows]); pdof = np.array([r[4] for r in rows])
print(f"pooled {len(conds)} unstranded scenarios; {len(rows)} positions with ≥{_MINOBS} informative scenarios\n")
print("PER-POSITION corr_i (hold node fixed, vary gDNA level) — the honest resolving-power metric")
print(f"{'type':>12} {'DOF':>11} | {'n_pos':>6} | {'mean corr_i':>11} | {'mass-wtd':>9} | {'frac>0.5':>8} | note")
for t in (1, 2, 3, 0):
    for d in (True, False):
        m = (ptyp == t) & (pdof == d)
        if m.sum() < 5:
            continue
        c = corr[m]
        mc, wc = float(np.mean(c)), float(np.average(c, weights=pmass[m]))
        fpos = float(np.mean(c > 0.5))
        note = "RESOLVES" if mc > 0.35 else ("coin-flip" if mc < 0.20 else "weak")
        print(f"{_TN[t]:>12} {('SOLVABLE' if d else 'UNSOLV'):>11} | {int(m.sum()):>6} | {mc:>11.3f} | "
              f"{wc:>9.3f} | {fpos:>8.2f} | {note}")

# also the coarse pooled corr per type (to show the Simpson gap side-by-side)
print("\nPOOLED corr (all node-observations, between-structure-confounded — for contrast)")
print(f"{'type':>12} {'DOF':>11} | {'n_obs':>6} | {'pooled corr':>11}")
for t in (1, 2, 3, 0):
    for d in (True, False):
        m = (A["typ"] == t) & (A["dof"] == d)
        if m.sum() < 5 or A["fo"][m].std() < 1e-9 or A["fg"][m].std() < 1e-9:
            continue
        c = np.corrcoef(A["fg"][m], A["fo"][m])[0, 1]
        print(f"{_TN[t]:>12} {('SOLVABLE' if d else 'UNSOLV'):>11} | {int(m.sum()):>6} | {c:>11.3f}")
