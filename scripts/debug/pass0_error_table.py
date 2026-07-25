"""PASS-0 STATE OF PLAY — the full suite in READS, with the TRUST view.

`pass0_oracle_bench.py` reports mwae (a rate); this reports the **absolute error in fragments** — how many
reads pass-0 assigns to the wrong side of the gDNA/RNA split — because that is what the downstream EM and the
hyperprior fit actually inherit. It also reports the question that decides whether an error MATTERS:

    is the solver's own stated precision HONEST about it?

Pass-0's deliverable is not a low score. It is a **substrate of trustworthy nodes**: if the nodes that are
wrong also say they are uncertain, the hyperprior fit can down-weight them and the calibration still
succeeds. Error on nodes that are CONFIDENTLY wrong is the kind that poisons everything downstream. So every
row splits the error mass by the node's own posterior ``Var(log f_g)``, at the suite-wide QUARTILES (a
data-defined split, not a threshold): an honest solver puts far less than 25 % of its error in the
most-confident quartile.

    OMP_NUM_THREADS=1 python scripts/debug/pass0_error_table.py [--refit 0] [--out /tmp/pass0_state.npz]
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
ap.add_argument("--out", default="/tmp/pass0_state.npz")
a = ap.parse_args()

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

C = {k: [] for k in ("cond", "mass", "err", "amb", "var", "cls", "self")}
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
    mass = np.asarray(cap["mass_global"])
    ok = np.isfinite(fo) & (mass > _EPS) & np.asarray(cap["solvable"], bool)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
    cls = np.array([CLS[int(rt[j])] if kind[j] == REGION else "boundary" for j in range(len(mass))])
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    C["cond"] += [cond] * int(ok.sum())
    C["mass"] += mass[ok].tolist()
    C["err"] += (np.abs(np.asarray(cap["f_g"]) - fo) * mass)[ok].tolist()
    C["amb"] += (fp & fn)[ok].tolist()
    C["var"] += np.asarray(cap["var_g"])[ok].tolist()
    C["cls"] += cls[ok].tolist()
    C["self"] += (np.abs(np.asarray(cap["fg_loc"]) - fo) * mass)[ok].tolist()
    print(f"  [{i + 1:>2}/{len(conds)}] {cond}", flush=True)

d = {k: np.asarray(v) for k, v in C.items()}
np.savez_compressed(a.out, **d)
cond, mass, err, amb, var = d["cond"], d["mass"], d["err"], d["amb"].astype(bool), d["var"]
sel = d["self"]

# the TRUST split — suite-wide quartiles of the node's OWN stated Var(log f_g) (data-defined, no threshold)
finite = np.isfinite(var)
q1, q3 = np.quantile(var[finite], [0.25, 0.75])
conf = finite & (var <= q1)  # most-confident quartile
unsure = (~finite) | (var >= q3)  # least-confident quartile (non-finite = no opinion at all)

print(f"\n{'=' * 132}\nPASS-0 STATE OF PLAY  (refit={a.refit})   "
      f"error = Σ mass·|f_g − oracle|, i.e. FRAGMENTS on the wrong side of the gDNA/RNA split\n{'=' * 132}")
print(f"{'scenario':<48}{'reads':>12}{'ERR reads':>11}{'mwae':>8}{'selfERR':>10}|"
      f"{'single':>10}{'AMBIG':>10}|{'errQ1conf':>10}{'errQ4unsure':>12}")
order = sorted(set(cond), key=lambda c: -err[cond == c].sum())
for c in order:
    m = cond == c
    e = err[m].sum()
    print(f"{c[5:]:<48}{mass[m].sum():>12,.0f}{e:>11,.0f}{e / mass[m].sum():>8.4f}"
          f"{sel[m].sum():>10,.0f}|{err[m & ~amb].sum():>10,.0f}{err[m & amb].sum():>10,.0f}|"
          f"{err[m & conf].sum() / max(e, _EPS):>10.1%}{err[m & unsure].sum() / max(e, _EPS):>12.1%}")
E = err.sum()
print(f"{'-' * 132}\n{'TOTAL':<48}{mass.sum():>12,.0f}{E:>11,.0f}{E / mass.sum():>8.4f}{sel.sum():>10,.0f}|"
      f"{err[~amb].sum():>10,.0f}{err[amb].sum():>10,.0f}|"
      f"{err[conf].sum() / E:>10.1%}{err[unsure].sum() / E:>12.1%}")

print(f"\n{'axis rollup':<34}{'reads':>13}{'ERR reads':>12}{'mwae':>8}{'share of ERR':>14}")
AX = {"capture off": lambda c: "capture_off" in c, "capture on": lambda c: "capture_on" in c,
      "capture verystrong": lambda c: "verystrong" in c,
      "stranded ss_0.99": lambda c: "ss_0.99" in c, "unstranded ss_0.50": lambda c: "ss_0.50" in c,
      "nRNA present": lambda c: "nrna_present" in c, "nRNA none": lambda c: "nrna_none" in c,
      "  unstranded × capON": lambda c: "ss_0.50" in c and "capture_on" in c,
      "  unstranded × verystrong": lambda c: "ss_0.50" in c and "verystrong" in c,
      "  stranded × capON": lambda c: "ss_0.99" in c and "capture_on" in c}
for lab, f in AX.items():
    m = np.array([f(c) for c in cond])
    if not m.any():
        continue
    print(f"{lab:<34}{mass[m].sum():>13,.0f}{err[m].sum():>12,.0f}"
          f"{err[m].sum() / mass[m].sum():>8.4f}{err[m].sum() / E:>14.1%}")

print(f"\n{'node class':<20}{'reads':>13}{'ERR reads':>12}{'mwae':>8}{'share of ERR':>14}"
      f"{'errQ1conf':>11}")
for c in ("exon", "boundary", "intron", "intergenic"):
    for lab, m2 in ((" single", ~amb), (" AMBIG", amb)):
        m = (d["cls"] == c) & m2
        if not m.any():
            continue
        print(f"{c + lab:<20}{mass[m].sum():>13,.0f}{err[m].sum():>12,.0f}"
              f"{err[m].sum() / max(mass[m].sum(), _EPS):>8.4f}{err[m].sum() / E:>14.1%}"
              f"{err[m & conf].sum() / max(err[m].sum(), _EPS):>11.1%}")
print(f"\n  TRUST: Var(log f_g) quartiles q1={q1:.4g} q3={q3:.4g}; "
      f"{(~finite).sum():,} nodes have no finite variance. An HONEST solver keeps the most-confident "
      f"quartile's share of error well under 25 %.")
