"""TARGET step 2 — partition the exon error by WHICH CHANNELS ARRIVED.

Step 1: at kappa = 1/2 every exon self-solves to ~0.5 (the uninformed default), so an exon's answer is 100 %
IMPORTED. The worst-node table showed at least three distinct modes:
  (a) node 3041: ALL channel precisions 0 -- no message arrived at all, stuck at 0.542, 22,806 error reads;
  (b) node 1075: cm_g = 0 but cm_p = 2.62 -- an RNA-ONLY message, dragged to f_g = 0.305 (oracle 0.822);
  (c) node 3111/1909/1167: a CONFIDENT but WRONG gDNA message (cm_g 17-27, mo_g 0.65-0.80).

Those are different bugs. This measures how much error mass each carries, and how far each exon sits from the
nearest node that has intrinsic information (an intergenic struct_lock anchor or a factory-solved intron) --
because at kappa = 1/2 those are the ONLY two sources in the entire solve.

    OMP_NUM_THREADS=1 python scratchpad/t2_strata.py
"""

from __future__ import annotations

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
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.50_nrna_none_capture_on"

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
dbg: dict = {}
calmod.calibrate(
    inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
    np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
)
chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
uni, us = cap["_uni"][-1], cap["_uni_static"]
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
T = Gp + Gn + Rp + Rn
fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
mass = np.asarray(cap["mass_global"])
ship = np.asarray(cap["f_g"])
self_fg = np.asarray(cap["fg_loc"])
solvable = np.asarray(cap["solvable"], bool)
tau0 = np.asarray(cap["_tau0_lam"])
struct = np.asarray(us["struct_lock"], bool)
rt, _ = _node_region_type(chain, ra)
kind = np.asarray(chain.kind)
is_exon = (kind == REGION) & (rt == 2)
ok = np.isfinite(fo) & (mass > _EPS) & solvable
err = np.abs(ship - fo) * mass
E = err[ok].sum()

cg, cpn = uni["cm_g"], uni["cm_p"] + uni["cm_n"]
ct = uni["c_tau"]
has_g, has_r, has_l = cg > _EPS, cpn > _EPS, ct > _EPS
anymsg = has_g | has_r | has_l

# ── distance (in chain hops) to the nearest node with INTRINSIC information ────────────────────────────────
# At kappa = 1/2 there are exactly two: a struct_lock intergenic anchor, and a factory-solved intron.
src = struct | (tau0 > _EPS)
left, right = np.asarray(chain.left), np.asarray(chain.right)
INF = 10**6
dist = np.where(src, 0, INF).astype(np.int64)
order = np.asarray(chain.order)
for _ in range(2):  # forward then backward sweeps until stable (chain is a forest of paths)
    for i in order:
        li = left[i]
        if li >= 0 and dist[li] + 1 < dist[i]:
            dist[i] = dist[li] + 1
    for i in order[::-1]:
        ri = right[i]
        if ri >= 0 and dist[ri] + 1 < dist[i]:
            dist[i] = dist[ri] + 1

print(f"{COND}\n  ERR = {E:,.0f} reads over {int(ok.sum()):,} solvable nodes; "
      f"{int((src & ok).sum())} of them carry INTRINSIC information "
      f"(struct_lock anchors {int((struct & ok).sum())}, factory introns {int(((tau0 > _EPS) & ok).sum())})")

print(f"\n  CHANNEL ARRIVAL — which messages reached the node (EXONS only)")
print(f"  {'stratum':<26}{'n':>6}{'reads':>12}{'ERR':>11}{'share':>7}{'mwae':>8}{'self':>8}{'signed':>9}")
strata = {
    "both gDNA + RNA": has_g & has_r,
    "gDNA only": has_g & ~has_r,
    "RNA only": ~has_g & has_r,
    "neither (lam only)": ~has_g & ~has_r & has_l,
    "NOTHING arrived": ~anymsg,
}
for lab, m0 in strata.items():
    m = ok & is_exon & m0
    if m.sum() < 1:
        continue
    w = mass[m]
    print(f"  {lab:<26}{int(m.sum()):>6}{w.sum():>12,.0f}{err[m].sum():>11,.0f}"
          f"{err[m].sum() / E:>7.1%}{np.average(np.abs(ship - fo)[m], weights=w):>8.4f}"
          f"{np.average(np.abs(self_fg - fo)[m], weights=w):>8.4f}"
          f"{np.average((ship - fo)[m], weights=w):>+9.4f}")

print(f"\n  DISTANCE to the nearest intrinsic-information node (EXONS only)")
print(f"  {'hops':<8}{'n':>6}{'reads':>12}{'ERR':>11}{'share':>7}{'mwae':>8}{'signed':>9}{'cm_g':>9}")
for lo, hi, lab in ((0, 1, "0-1"), (2, 3, "2-3"), (4, 5, "4-5"), (6, 9, "6-9"),
                    (10, 10**5, "10+"), (INF, INF, "unreachable")):
    m = ok & is_exon & (dist >= lo) & (dist <= hi)
    if m.sum() < 1:
        continue
    w = mass[m]
    print(f"  {lab:<8}{int(m.sum()):>6}{w.sum():>12,.0f}{err[m].sum():>11,.0f}"
          f"{err[m].sum() / E:>7.1%}{np.average(np.abs(ship - fo)[m], weights=w):>8.4f}"
          f"{np.average((ship - fo)[m], weights=w):>+9.4f}{np.average(cg[m], weights=w):>9.2f}")

# ── is the gDNA message ACCURATE where it arrives? (the (c) mode) ─────────────────────────────────────────
print(f"\n  gDNA MESSAGE ACCURACY on exons where it arrived (mo_g is a claimed f_g)")
print(f"  {'cm_g band':<14}{'n':>6}{'reads':>12}{'ERR':>11}{'share':>7}"
      f"{'|mo_g-orc|':>11}{'bias':>9}{'z2':>8}")
mg = np.clip(np.exp(uni["mo_g"]), 0.0, 1.0)
for lo, hi in ((0, 1), (1, 5), (5, 20), (20, 1e9)):
    m = ok & is_exon & has_g & (cg >= lo) & (cg < hi)
    if m.sum() < 3:
        continue
    w = mass[m]
    dm = np.log(np.maximum(mg[m], 1e-9)) - np.log(np.maximum(fo[m], 1e-9))
    z2 = np.average(dm * dm, weights=w) / np.average(1.0 / np.maximum(cg[m], _EPS), weights=w)
    print(f"  {f'{lo:g}-{hi:g}':<14}{int(m.sum()):>6}{w.sum():>12,.0f}{err[m].sum():>11,.0f}"
          f"{err[m].sum() / E:>7.1%}{np.average(np.abs(mg - fo)[m], weights=w):>11.4f}"
          f"{np.average((mg - fo)[m], weights=w):>+9.4f}{z2:>8.1f}")
