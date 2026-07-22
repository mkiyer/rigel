"""Node-type-by-node-type test of the DOF solvability criterion: for the nodes it would mark UNSOLVABLE, is the
current forced-solve f_g MEANINGFUL (correlated with the oracle) or a COIN-FLIP (uncorrelated)? If coin-flip, we
were right to withhold them; if correlated, we are wrongly judging them unsolvable. Pools the UNSTRANDED (ss0.50)
ambig_dense_10mb scenarios. Error is NOT the metric here (unsolved f_g is an arbitrary default) — CORRELATION is."""
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
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_") and "0.50" in d.name)

rows = {}  # (type, solvable) -> lists
def push(k, fg, fo, m):
    rows.setdefault(k, ([], [], []))
    rows[k][0].append(fg); rows[k][1].append(fo); rows[k][2].append(m)

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
    typ = np.where(kind != REGION, "boundary", np.where(nfree >= 2, "AMBIG-reg", np.where(nfree == 1, "single-reg", "intergenic")))
    ok = np.isfinite(fo) & (mass > _EPS)
    for i in np.where(ok)[0]:
        push((typ[i], bool(dof_solv[i])), fg[i], fo[i], mass[i])

print(f"pooled {len(conds)} unstranded scenarios\n")
print(f"{'node type':>12} {'DOF':>10} | {'n':>5} {'mass%':>6} | {'corr(fg,oracle)':>15} | {'mean|fg-or|':>11} | note")
tot_mass = sum(sum(v[2]) for v in rows.values())
for typ in ("single-reg", "AMBIG-reg", "boundary", "intergenic"):
    for solv in (True, False):
        k = (typ, solv)
        if k not in rows:
            continue
        fg = np.array(rows[k][0]); fo = np.array(rows[k][1]); m = np.array(rows[k][2])
        corr = np.corrcoef(fg, fo)[0, 1] if len(fg) > 2 and fg.std() > 1e-9 and fo.std() > 1e-9 else float("nan")
        err = float(np.average(np.abs(fg - fo), weights=m))
        note = "COIN-FLIP → withhold" if (not np.isnan(corr) and abs(corr) < 0.15) else ("meaningful" if not np.isnan(corr) and corr > 0.3 else "")
        print(f"{typ:>12} {('SOLVABLE' if solv else 'UNSOLVABLE'):>10} | {len(fg):>5} {100*m.sum()/tot_mass:>5.1f}% | {corr:>15.3f} | {err:>11.3f} | {note}")
