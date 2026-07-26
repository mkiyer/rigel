"""TARGET step 3 — the MEASUREMENT-STREAM ASYMMETRY hypothesis.

Step 2, on exons, by which channels arrived:
    both gDNA+RNA   414 nodes  943 k err (69.4 %)  signed -0.1262
    gDNA only       144 nodes   74 k err ( 5.4 %)  signed -0.0053   <-- essentially UNBIASED
    RNA only         90 nodes  200 k err (14.7 %)  signed -0.2070   <-- catastrophically RNA-biased
    NOTHING           6 nodes   33 k err ( 2.4 %)  stuck at the 0.5 default

and the gDNA message is well calibrated where it arrives (z2 = 0.4-1.5, bias -> -0.02 at high precision).
So the gDNA channel is sound and the imbalance is what hurts.

HYPOTHESIS. In `bp_solver`, the gDNA MEASUREMENT stream is seeded ONLY at struct_lock nodes
(`mg_own = where(struct_lock, pg_own, 0)`), while the RNA measurement stream is seeded at every spliced GRAFT.
But in THIS scenario the intrinsic information is 572 factory-solved introns and (checked below) very few
usable struct_lock anchors -- and a factory intron's confident gDNA density is NOT admitted to the
measurement stream. So an exon deep inside a gene sees a spliced RNA measurement with no gDNA counterpart.

    OMP_NUM_THREADS=1 python scratchpad/t3_streams.py
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
ship, self_fg = np.asarray(cap["f_g"]), np.asarray(cap["fg_loc"])
solvable = np.asarray(cap["solvable"], bool)
tau0 = np.asarray(cap["_tau0_lam"])
struct = np.asarray(us["struct_lock"], bool)
mg_own = np.asarray(us["mg_own"])
rt, _ = _node_region_type(chain, ra)
kind = np.asarray(chain.kind)
CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
cls = np.array([CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(len(mass))])
ok = np.isfinite(fo) & (mass > _EPS) & solvable
err = np.abs(ship - fo) * mass
E = err[ok].sum()
is_exon = cls == "exon"

print(f"{COND}\n{'=' * 110}")
print("WHO SEEDS THE TWO MEASUREMENT STREAMS?")
print(f"  struct_lock nodes (seed mg_own):      {int(struct.sum()):>6}   "
      f"of which mg_own>0: {int((mg_own > 0).sum()):>6}   total mass {mass[struct].sum():>12,.0f}")
print(f"  factory introns (tau_lam>0):          {int((tau0 > _EPS).sum()):>6}   "
      f"mg_own>0 among them: {int(((tau0 > _EPS) & (mg_own > 0)).sum()):>6}   "
      f"total mass {mass[tau0 > _EPS].sum():>12,.0f}")
sp = np.asarray(us['SP_l']) + np.asarray(us['SP_r']) + np.asarray(us['SN_l']) + np.asarray(us['SN_r'])
print(f"  boundaries with spliced (seed RNA):   {int((sp > _EPS).sum()):>6}   "
      f"total spliced {sp.sum():>12,.0f}")
print(f"  by class, struct_lock count: " + "  ".join(
    f"{c}={int((struct & (cls == c)).sum())}" for c in ("intergenic", "intron", "exon", "boundary")))

# where do the exons' gDNA vs RNA measurement precisions actually come from?
cg, cpn = uni["cm_g"], uni["cm_p"] + uni["cm_n"]
print(f"\nEXON MEASUREMENT PRECISION, mass-weighted (cm_g = gDNA count stream, cm_p+n = spliced RNA stream)")
w = mass[ok & is_exon]
print(f"  exons: n={int((ok & is_exon).sum())}  mean cm_g={np.average(cg[ok & is_exon], weights=w):.2f}  "
      f"mean cm_p+n={np.average(cpn[ok & is_exon], weights=w):.2f}   "
      f"ratio RNA/gDNA={np.average(cpn[ok & is_exon], weights=w) / max(np.average(cg[ok & is_exon], weights=w), _EPS):.2f}")

# ── the RNA-only exons: is there a factory intron in the neighbourhood whose gDNA is being ignored? ───────
left, right = np.asarray(chain.left), np.asarray(chain.right)


def nbrs(i, hops=2):
    """node ids within `hops` chain steps."""
    seen, frontier = {i}, [i]
    for _ in range(hops):
        nxt = []
        for j in frontier:
            for k in (left[j], right[j]):
                if k >= 0 and k not in seen:
                    seen.add(int(k))
                    nxt.append(int(k))
        frontier = nxt
    return seen - {i}


rna_only = ok & is_exon & (cg <= _EPS) & (cpn > _EPS)
gd_only = ok & is_exon & (cg > _EPS) & (cpn <= _EPS)
none_ar = ok & is_exon & (cg <= _EPS) & (cpn <= _EPS) & (uni["c_tau"] <= _EPS)
for lab, m in (("RNA-only", rna_only), ("gDNA-only", gd_only), ("NOTHING", none_ar)):
    idx = np.flatnonzero(m)
    if idx.size == 0:
        continue
    near_fac = np.array([any(tau0[j] > _EPS for j in nbrs(int(i), 2)) for i in idx])
    near_str = np.array([any(struct[j] for j in nbrs(int(i), 3)) for i in idx])
    wm = mass[idx]
    print(f"\n  {lab}: n={idx.size}  err={err[idx].sum():,.0f} ({err[idx].sum() / E:.1%})  "
          f"has a FACTORY intron within 2 hops: {np.average(near_fac, weights=wm):.0%}   "
          f"a struct_lock within 3 hops: {np.average(near_str, weights=wm):.0%}")

print(f"\n  WORST 'NOTHING arrived' and 'RNA-only' exons")
print(f"  {'node':>6}{'stratum':>10}{'reads':>10}{'orc':>6}{'self':>6}{'ship':>6}{'ERR':>9}"
      f"{'cm_g':>7}{'cm_p+n':>8}{'c_tau':>7}|{'nbr classes (2 hops)':<44}")
cand = np.flatnonzero((rna_only | none_ar))
for i in cand[np.argsort(-err[cand])][:12]:
    nb = sorted(nbrs(int(i), 2))
    tag = "NOTHING" if none_ar[i] else "RNA-only"
    desc = ", ".join(f"{cls[j][:4]}{'*' if tau0[j] > _EPS else ''}{'L' if struct[j] else ''}" for j in nb)
    print(f"  {i:>6}{tag:>10}{mass[i]:>10,.0f}{fo[i]:>6.3f}{self_fg[i]:>6.3f}{ship[i]:>6.3f}"
          f"{err[i]:>9,.0f}{cg[i]:>7.2f}{cpn[i]:>8.2f}{uni['c_tau'][i]:>7.2f}|{desc[:44]:<44}")
print("\n  (* = factory intron, L = struct_lock)")
