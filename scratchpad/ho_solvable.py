"""TASK 5 — "SOLVABLE PRIOR-FREE": what predicates exist today, do they agree, and do they separate
accurate from inaccurate nodes?

Candidates measured, per node, over the whole suite:
  cap["solvable"]        = (free_pos | free_neg) & mass_unspliced > 0      (bp_solver's WRITE-BACK gate)
  struct_lock            = ~solvable & is_region                           (node_init source 1: MEASURED)
  tau_own  (_tau0_lam)   = I_strand·[single-strand] + I_density(factory)   (node_init sources 2+3)
  single / ambig / gonly = free_pos ^ free_neg / & / ~&~                   (statics, the DOF class)

and the four node_init SOURCES are reconstructed from them.

    OMP_NUM_THREADS=1 python scratchpad/ho_solvable.py [--conds N]
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
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9

ap = argparse.ArgumentParser()
ap.add_argument("--refit", type=int, default=0)
a = ap.parse_args()

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

COL: dict[str, list] = {k: [] for k in (
    "cond", "mass", "err", "eloc", "var", "tau", "solv", "lock", "isr", "amb", "single",
    "gonly", "rt", "strand",
)}
for i, cond in enumerate(conds):
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=a.refit), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    fg = np.asarray(cap["f_g"], float)
    isr = np.asarray(chain.kind) == REGION
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    rt, _ = _node_region_type(chain, ra)
    COL["cond"] += [cond] * fg.shape[0]
    COL["mass"] += np.asarray(cap["mass_global"], float).tolist()
    COL["err"] += np.abs(fg - fo).tolist()
    COL["eloc"] += np.abs(np.asarray(cap["fg_loc"], float) - fo).tolist()
    COL["var"] += np.asarray(cap["var_g"], float).tolist()
    COL["tau"] += np.asarray(cap["_tau0_lam"], float).tolist()
    COL["solv"] += np.asarray(cap["solvable"], bool).tolist()
    COL["lock"] += np.asarray(cap["_uni_static"]["struct_lock"], bool).tolist()
    COL["isr"] += isr.tolist()
    COL["amb"] += (fp & fn).tolist()
    COL["single"] += (fp ^ fn).tolist()
    COL["gonly"] += (~fp & ~fn).tolist()
    COL["rt"] += np.where(isr, rt, -1).astype(np.int64).tolist()
    COL["strand"] += ["ss_0.99" if "ss_0.99" in cond else "ss_0.50"] * fg.shape[0]
    print(f"  [{i + 1:>2}/{len(conds)}] {cond}", flush=True)

d = {k: np.asarray(v) for k, v in COL.items()}
mass, err, eloc, var, tau = d["mass"], d["err"], d["eloc"], d["var"], d["tau"]
ok = np.isfinite(err) & (mass > _EPS)
solv, lock, isr = d["solv"], d["lock"], d["isr"]
amb, single, gonly, rt = d["amb"], d["single"], d["gonly"], d["rt"]
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}


def _w(m, x=err):
    m = m & ok
    if not m.any():
        return float("nan")
    return float(np.average(x[m], weights=mass[m]))


def _z2(m):
    m = m & ok & np.isfinite(var)
    if not m.any():
        return float("nan")
    den = float(np.sum(mass[m] * var[m]))
    return float(np.sum(mass[m] * err[m] ** 2)) / den if den > 0 else float("nan")


print(f"\n{'=' * 108}\nPREDICATE AGREEMENT (whole suite, {int(ok.sum()):,} nodes with mass+oracle)\n{'=' * 108}")
print(f"  {'predicate':<34}{'nodes':>10}{'mass%':>9}{'mwae':>9}{'self mwae':>11}{'z2':>9}")
PREDS = [
    ("ALL", np.ones(ok.shape, bool)),
    ("solvable (write-back gate)", solv),
    ("  ~solvable (locked/empty)", ~solv),
    ("struct_lock (MEASURED)", lock),
    ("tau_own > 0", tau > _EPS),
    ("  tau_own == 0", tau <= _EPS),
    ("tau_own>0 OR struct_lock", (tau > _EPS) | lock),
    ("  ...and its complement", ~((tau > _EPS) | lock)),
    ("single-strand (DOF)", single),
    ("AMBIG (DOF)", amb),
    ("gonly (DOF)", gonly),
]
tot = float(mass[ok].sum())
for lab, m in PREDS:
    mm = m & ok
    print(f"  {lab:<34}{int(mm.sum()):>10,}{mass[mm].sum() / tot:>8.1%}"
          f"{_w(mm):>9.4f}{_w(mm, eloc):>11.4f}{_z2(mm):>9.2f}")

print(f"\n{'=' * 108}\nDO THEY AGREE? cross-tab of tau_own>0 against the DOF class and struct_lock\n{'=' * 108}")
print(f"  {'':<26}{'tau>0 n':>10}{'tau>0 mass%':>13}{'tau=0 n':>10}{'tau=0 mass%':>13}"
      f"{'mwae tau>0':>12}{'mwae tau=0':>12}")
for lab, m in (("single-strand", single), ("AMBIG", amb), ("gonly", gonly),
               ("struct_lock", lock), ("region", isr), ("boundary", ~isr)):
    a1, a0 = m & (tau > _EPS) & ok, m & (tau <= _EPS) & ok
    print(f"  {lab:<26}{int(a1.sum()):>10,}{mass[a1].sum() / tot:>12.1%}"
          f"{int(a0.sum()):>10,}{mass[a0].sum() / tot:>12.1%}{_w(a1):>12.4f}{_w(a0):>12.4f}")

print(f"\n{'=' * 108}\nTHE FOUR node_init SOURCES — does the source predict accuracy?\n{'=' * 108}")
SRC = [
    ("1 MEASURED (struct_lock)", lock),
    ("2+3 tau_own>0, intron", (tau > _EPS) & ~lock & (rt == 1)),
    ("2+3 tau_own>0, exon", (tau > _EPS) & ~lock & (rt == 2)),
    ("2+3 tau_own>0, boundary", (tau > _EPS) & ~lock & ~isr),
    ("4 UNSOLVED tau=0, exon", (tau <= _EPS) & ~lock & (rt == 2)),
    ("4 UNSOLVED tau=0, intron", (tau <= _EPS) & ~lock & (rt == 1)),
    ("4 UNSOLVED tau=0, boundary", (tau <= _EPS) & ~lock & ~isr),
]
print(f"  {'source':<32}{'nodes':>10}{'mass%':>9}{'mwae':>9}{'self mwae':>11}{'z2':>9}"
      f"{'  |  stranded / unstranded mwae':<34}")
for lab, m in SRC:
    mm = m & ok
    if not mm.any():
        continue
    s99 = mm & (d["strand"] == "ss_0.99")
    s50 = mm & (d["strand"] == "ss_0.50")
    print(f"  {lab:<32}{int(mm.sum()):>10,}{mass[mm].sum() / tot:>8.1%}{_w(mm):>9.4f}"
          f"{_w(mm, eloc):>11.4f}{_z2(mm):>9.2f}  |  {_w(s99):>10.4f} /{_w(s50):>10.4f}")

print(f"\n{'=' * 108}\nSEPARATION POWER — is tau_own>0 a usable 'this node can be solved prior-free' test?\n{'=' * 108}")
for lab, sub in (("ALL", ok), ("stranded ss_0.99", ok & (d["strand"] == "ss_0.99")),
                 ("unstranded ss_0.50", ok & (d["strand"] == "ss_0.50"))):
    a1, a0 = sub & (tau > _EPS), sub & (tau <= _EPS)
    ml = f"  {lab:<20} tau>0 mwae={_w(a1):.4f} (mass {mass[a1].sum() / tot:>5.1%})"
    ml += f"   tau=0 mwae={_w(a0):.4f} (mass {mass[a0].sum() / tot:>5.1%})"
    ml += f"   ratio {_w(a0) / max(_w(a1), 1e-9):>6.1f}x"
    print(ml)
    # error MASS share carried by each side
    em = err * mass
    print(f"  {'':<20} error-mass share: tau>0 {em[a1].sum() / em[sub].sum():>5.1%}"
          f"   tau=0 {em[a0].sum() / em[sub].sum():>5.1%}")
np.savez_compressed("/tmp/ho_solvable.npz", **d)
print("\nwrote /tmp/ho_solvable.npz")
