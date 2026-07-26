"""SUITE-WIDE pass-0 error interrogation — assumption-free: find where the error IS, then look at it.

Runs every scenario in the suite, records one row per node, and writes a single table for interrogation.
The organizing question is NOT "is the answer wrong" but **"do the MESSAGES make it worse than the node's own
message-free self-solve, and where"** — because a solvable node that self-solves well and is then dragged off
by a message is a MODEL BUG (like the intergenic-seam λ artifact), whereas a node that self-solves badly and
stays bad is an information problem (the hyperprior's job).

Scope note (owner, 2026-07-25): pass-0 has no population prior, so AMBIG / τ_own=0 nodes are NOT full rank and
cannot be solved here — they must end up either near-zero precision or unsolved. They are reported separately
and are NOT the target. **The target is SOLVABLE nodes that solve WRONG.**

    OMP_NUM_THREADS=1 python scratchpad/suite_dissect.py [--refit 0] [--out /tmp/suite_nodes.npz]
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
ap.add_argument("--out", default="/tmp/suite_nodes.npz")
ap.add_argument("--conds", default=None, help="comma-separated subset")
a = ap.parse_args()

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
conds = (
    a.conds.split(",")
    if a.conds
    else sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
)

cols: dict[str, list] = {k: [] for k in (
    "cond", "node", "cls", "dof", "solvable", "mass", "oracle", "self", "solved", "tau_own",
    "c_tau", "lam_fg", "cm_g", "mo_g", "cm_p", "mo_p", "cm_n", "mo_n",
    "nl_cls", "nl_oracle", "nl_mass", "nr_cls", "nr_oracle", "nr_mass", "n_unspl",
)}
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}

for ci, cond in enumerate(conds):
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=a.refit), _debug=dbg,
    )
    chain, cap = dbg["chain"], dbg["capture"]
    uni, us = cap["_uni"][-1], cap["_uni_static"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    fg, mass = np.asarray(cap["f_g"]), np.asarray(cap["mass_global"])
    self_fg = np.asarray(cap["fg_loc"])
    tau = np.asarray(cap["_tau0_lam"])
    fp = np.asarray(dbg["statics"].free_pos, bool)
    fn = np.asarray(dbg["statics"].free_neg, bool)
    solvable = np.asarray(cap["solvable"], bool)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    n = fg.shape[0]
    cls = np.array([CLS[int(rt[i])] if kind[i] == REGION else "boundary" for i in range(n)])
    dof = np.where(fp & fn, "ambig", np.where(fp ^ fn, "single", "locked"))
    keep = np.flatnonzero(np.isfinite(fo) & (mass > _EPS))
    nu = np.asarray(us["n_unspl_l"], np.float64) + np.asarray(us["n_unspl_r"], np.float64)

    def _nb(idx, arr, default):
        return np.array([arr[j] if j >= 0 else default for j in idx])

    li, ri = left[keep], right[keep]
    cols["cond"] += [cond] * keep.size
    cols["node"] += keep.tolist()
    cols["cls"] += cls[keep].tolist()
    cols["dof"] += dof[keep].tolist()
    for k, v in (("solvable", solvable), ("mass", mass), ("oracle", fo), ("self", self_fg),
                 ("solved", fg), ("tau_own", tau), ("c_tau", uni["c_tau"]),
                 ("cm_g", uni["cm_g"]), ("cm_p", uni["cm_p"]), ("cm_n", uni["cm_n"]), ("n_unspl", nu)):
        cols[k] += np.asarray(v)[keep].tolist()
    cols["lam_fg"] += (1.0 / (1.0 + np.exp(-np.asarray(uni["lam_msg"])[keep]))).tolist()
    for k, v in (("mo_g", "mo_g"), ("mo_p", "mo_p"), ("mo_n", "mo_n")):
        cols[k] += np.exp(np.asarray(uni[v])[keep]).tolist()
    cols["nl_cls"] += _nb(li, cls, "-").tolist()
    cols["nr_cls"] += _nb(ri, cls, "-").tolist()
    cols["nl_oracle"] += _nb(li, fo, np.nan).tolist()
    cols["nr_oracle"] += _nb(ri, fo, np.nan).tolist()
    cols["nl_mass"] += _nb(li, mass, 0.0).tolist()
    cols["nr_mass"] += _nb(ri, mass, 0.0).tolist()
    print(f"  [{ci + 1:>2}/{len(conds)}] {cond}", flush=True)

out = {k: np.asarray(v) for k, v in cols.items()}
np.savez_compressed(a.out, **out)
err = np.abs(out["solved"] - out["oracle"])
serr = np.abs(out["self"] - out["oracle"])
em = err * out["mass"]
print(f"\nwrote {err.size:,} node rows -> {a.out}")
print(f"total error mass {em.sum():,.0f}   suite mwae {np.average(err, weights=out['mass']):.4f}")

# ── the headline split: where does the error mass live, and do MESSAGES help or hurt there? ──
print(f"\n{'stratum':<30}{'nodes':>8}{'err-mass':>12}{'share':>8}{'self-mwae':>11}{'solved-mwae':>13}{'Δ msg':>9}")


def row(lab, m):
    if not m.any():
        return
    w = out["mass"][m]
    print(f"{lab:<30}{int(m.sum()):>8}{em[m].sum():>12,.0f}{em[m].sum() / em.sum():>7.1%}"
          f"{np.average(serr[m], weights=w):>11.4f}{np.average(err[m], weights=w):>13.4f}"
          f"{np.average(err[m], weights=w) - np.average(serr[m], weights=w):>+9.4f}")


solv, amb = out["solvable"].astype(bool), out["dof"] == "ambig"
row("ALL", np.ones_like(solv))
row("SOLVABLE", solv)
row("  solvable & single-strand", solv & ~amb)
row("  solvable & AMBIG", solv & amb)
row("NOT solvable (locked)", ~solv)
print()
for c in ("exon", "intron", "boundary", "intergenic"):
    row(f"SOLVABLE single-strand {c}", solv & ~amb & (out["cls"] == c))
