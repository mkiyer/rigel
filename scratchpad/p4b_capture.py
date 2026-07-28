"""P4b — price the `_RHO_ITERS = 2` BP violation. CAPTURE pass.

Runs the full 32-condition suite for one (arm, refit) and dumps every per-node quantity the report needs:
mass, oracle f_g, solved f_g, stated Var(log f_g), AMBIG flag, class, solvable, and the hyperprior's own
fit-substrate membership (`_fit_gdna_hyperprior`'s `live & isr & (single|gonly)`), plus — on the 2-iteration
arm — the per-node |Δ log ρ_face| between the two lazy-ρ iterations.

    RIGEL_RHO_ITERS=1 python scratchpad/p4b_capture.py --arm rho1 --refit 0 --out /tmp/p4b_rho1_r0.npz
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
ap.add_argument("--arm", required=True)
ap.add_argument("--refit", type=int, default=0)
ap.add_argument("--out", required=True)
a = ap.parse_args()

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

KEYS = ("cond", "nid", "mass", "fo", "fg", "var", "amb", "cls", "solv", "fit", "dlr", "selferr")
C: dict = {k: [] for k in KEYS}
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
    eff = np.asarray(cap["eff_global"])
    fg = np.asarray(cap["f_g"])
    kind = np.asarray(chain.kind)
    rt, _ = _node_region_type(chain, ra)
    CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
    cls = np.array([CLS[int(rt[j])] if kind[j] == REGION else "boundary" for j in range(len(mass))])
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    # the hyperprior's OWN training-substrate selector (`calibrate._fit_gdna_hyperprior`, additive=False,
    # background=None ⇒ no intergenic drop): live & REGION & (single | gonly)
    single = fp ^ fn
    gonly = ~fp & ~fn
    live = (eff > 1.0e-9) & (mass > 1.0e-12)
    fit_sub = live & (kind == REGION) & (single | gonly)
    # |Δ log ρ_face| between the two lazy-ρ iterations (only meaningful on the 2-iteration arm)
    uni = cap.get("_uni") or []
    if len(uni) >= 2:
        d0 = np.abs(np.log(np.maximum(uni[1]["rho_lf"], _EPS)) - np.log(np.maximum(uni[0]["rho_lf"], _EPS)))
        d1 = np.abs(np.log(np.maximum(uni[1]["rho_rf"], _EPS)) - np.log(np.maximum(uni[0]["rho_rf"], _EPS)))
        dlr = np.maximum(d0, d1)
    else:
        dlr = np.zeros_like(mass)
    ok = np.isfinite(fo) & (mass > _EPS)
    n = int(ok.sum())
    C["cond"] += [cond] * n
    C["nid"] += np.nonzero(ok)[0].tolist()
    C["mass"] += mass[ok].tolist()
    C["fo"] += fo[ok].tolist()
    C["fg"] += fg[ok].tolist()
    C["var"] += np.asarray(cap["var_g"])[ok].tolist()
    C["amb"] += (fp & fn)[ok].tolist()
    C["cls"] += cls[ok].tolist()
    C["solv"] += np.asarray(cap["solvable"], bool)[ok].tolist()
    C["fit"] += fit_sub[ok].tolist()
    C["dlr"] += dlr[ok].tolist()
    C["selferr"] += np.abs(np.asarray(cap["fg_loc"]) - fo)[ok].tolist()
    print(f"  [{i + 1:>2}/{len(conds)}] {a.arm} r{a.refit} {cond}", flush=True)

np.savez_compressed(a.out, **{k: np.asarray(v) for k, v in C.items()})
print(f"wrote {a.out}")
