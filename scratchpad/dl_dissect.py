"""DL residual dissection — where is the stranded arm's remaining error mass, and is DL protecting there?

The DL composition-mismatch term damps a message ONLY where the destination has its own composition evidence
(τ_own > 0 ⇒ v_own finite). Where τ_own = 0 (AMBIG boundaries; every node on unstranded data) v_own = ∞ and
DL disables itself BY DESIGN — that is what preserves the M5 unstranded/capture win. So the question this
script answers is: is the residual stranded error concentrated on the DL-PROTECTED nodes (⇒ the law is too
weak / mis-specified) or on the UNPROTECTED τ_own=0 nodes (⇒ it is the known Phase-2 gap: AMBIG has no own
opinion to contradict a wrong loud neighbour, and the hyperprior is what supplies one)?

    OMP_NUM_THREADS=1 python scratchpad/dl_dissect.py [COND] [--top 20]
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
_EPS = 1e-9

ap = argparse.ArgumentParser()
ap.add_argument("cond", nargs="?", default="gdna_gdna300_ss_0.99_nrna_present_capture_on")
ap.add_argument("--top", type=int, default=20)
ap.add_argument("--refit", type=int, default=0)
ap.add_argument("--only", choices=["protected", "unprotected"], help="restrict the top-N listing")
a = ap.parse_args()

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, a.cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
dbg: dict = {}
cc = dataclasses.replace(cfg.calibration, calib_refit_iters=a.refit)
calmod.calibrate(
    inp["payload"], ra, inp["strand_model"],
    np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
)
chain, cap = dbg["chain"], dbg["capture"]
uni = cap["_uni"][-1]
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
G, R = Gp + Gn, Rp + Rn
fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
fg = np.asarray(cap["f_g"])
mass = np.asarray(cap["mass_global"])
rt, _ = _node_region_type(chain, ra)
kind = np.asarray(chain.kind)
fp = np.asarray(dbg["statics"].free_pos, bool)
fn = np.asarray(dbg["statics"].free_neg, bool)
tau0 = np.asarray(cap["_tau0_lam"])
fgloc = np.asarray(cap["fg_loc"])
left, right = np.asarray(chain.left), np.asarray(chain.right)
ok = np.isfinite(fo) & (mass > _EPS)
err = np.abs(fg - fo)
emass = err * mass

print(f"# {a.cond}   refit={a.refit}")
print(f"# nodes with mass+oracle: {int(ok.sum()):,}   mwae={np.average(err[ok], weights=mass[ok]):.4f}")

# ── attribution: is the error where DL protects (τ_own>0) or where it disables itself (τ_own=0)? ──
tot = emass[ok].sum()
print(f"\n{'stratum':<34}{'nodes':>8}{'mass':>14}{'err-mass':>12}{'share':>8}{'mwae':>9}")


def _row(label, m):
    m = m & ok
    if not m.any():
        return
    print(f"{label:<34}{int(m.sum()):>8}{mass[m].sum():>14,.0f}{emass[m].sum():>12,.0f}"
          f"{emass[m].sum() / tot:>7.1%}{np.average(err[m], weights=mass[m]):>9.4f}")


amb = fp & fn
single = fp ^ fn
print("── by DL protection state (v_own from τ_own) ──")
_row("PROTECTED   τ_own>0", tau0 > _EPS)
_row("UNPROTECTED τ_own=0 (v_own=∞)", tau0 <= _EPS)
print("── by strand DOF ──")
_row("AMBIG (both strands free)", amb)
_row("single-strand", single)
_row("no free strand (locked)", ~(fp | fn))
print("── τ_own=0 split ──")
_row("τ=0 & AMBIG", (tau0 <= _EPS) & amb)
_row("τ=0 & single-strand", (tau0 <= _EPS) & single)
_row("τ=0 & locked", (tau0 <= _EPS) & ~(fp | fn))
print("── by node class ──")
for k, name in CLS.items():
    _row(name, (rt == k) if k >= 0 else (kind != REGION))

# ── the local (message-free) self-solve error, for the same strata: how much does the message ADD? ──
print("\n── message harm: |self-solve − oracle| vs |solved − oracle| (mass-weighted) ──")
eloc = np.abs(fgloc - fo)
print(f"{'stratum':<34}{'self-solve':>12}{'solved':>10}{'Δ (msg)':>10}")
for label, m in (
    ("ALL", ok),
    ("PROTECTED τ_own>0", ok & (tau0 > _EPS)),
    ("UNPROTECTED τ_own=0", ok & (tau0 <= _EPS)),
    ("τ=0 & AMBIG", ok & (tau0 <= _EPS) & amb),
    ("τ=0 & locked", ok & (tau0 <= _EPS) & ~(fp | fn)),
):
    if not m.any():
        continue
    el = np.average(eloc[m], weights=mass[m])
    es = np.average(err[m], weights=mass[m])
    print(f"{label:<34}{el:>12.4f}{es:>10.4f}{es - el:>+10.4f}")

# ── the worst nodes by error MASS, with their message provenance ──
sel_top = ok
if a.only == "protected":
    sel_top = ok & (tau0 > _EPS)
elif a.only == "unprotected":
    sel_top = ok & (tau0 <= _EPS)
idx = np.argsort(np.where(sel_top, emass, -1.0))[::-1][: a.top]
print(f"\n── top {a.top} nodes by error MASS ──")
for i in idx:
    i = int(i)
    cl = CLS[int(rt[i])] if kind[i] == REGION else "boundary"
    dof = "AMBIG" if amb[i] else ("single" if single[i] else "locked")
    print(f"\nnode {i} [{cl}/{dof}] mass={mass[i]:,.0f} err-mass={emass[i]:,.0f}  "
          f"oracle={fo[i]:.3f} solved={fg[i]:.3f} self-solve={fgloc[i]:.3f} τ_own={tau0[i]:.4g}")
    print(f"   λ-msg f_g_eq={1 / (1 + np.exp(-uni['lam_msg'][i])):.3f} prec={uni['c_tau'][i]:.4g}   "
          f"gDNA-meas mode={np.exp(uni['mo_g'][i]):.3f} prec={uni['cm_g'][i]:.4g}   "
          f"RNA-meas +={np.exp(uni['mo_p'][i]):.3f}/{uni['cm_p'][i]:.4g} −={np.exp(uni['mo_n'][i]):.3f}/{uni['cm_n'][i]:.4g}")
    for tag, j in (("L", int(left[i])), ("R", int(right[i]))):
        if j >= 0:
            clj = CLS[int(rt[j])] if kind[j] == REGION else "boundary"
            print(f"   nbr {tag} {j} [{clj}] oracle={fo[j]:.3f} solved={fg[j]:.3f} "
                  f"self={fgloc[j]:.3f} τ={tau0[j]:.4g} mass={mass[j]:,.0f}")

if os.environ.get("DL_SUMMARY_ONLY"):
    pass
