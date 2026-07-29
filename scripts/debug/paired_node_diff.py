"""PAIRED PER-NODE DIFF between two solver arms — which exact nodes cause a regression, and why.

The aggregate mwae delta is an average; it cannot say whether damage is CONCENTRATED at one topology (a bug)
or SPREAD thinly (a systematic small bias -- the honest contention between messages and strand). This runs the
same scenarios under two arms, joins per node, and ranks by mass-weighted error INCREASE.

    OMP_NUM_THREADS=1 python scripts/debug/paired_node_diff.py --arm base
    RIGEL_ROUTE_TEST=1 OMP_NUM_THREADS=1 python scripts/debug/paired_node_diff.py --arm route
    OMP_NUM_THREADS=1 python scripts/debug/paired_node_diff.py --report

``--set`` selects the scenario family (default: stranded ss_0.99).
"""
from __future__ import annotations

import argparse
import dataclasses
import importlib
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION, node_global_geometry  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

_EPS = 1e-9
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
OUT = Path("/tmp/paired_node_diff")
OUT.mkdir(exist_ok=True)

ap = argparse.ArgumentParser()
ap.add_argument("--arm")
ap.add_argument("--report", action="store_true")
ap.add_argument("--set", default="ss_0.99")
args = ap.parse_args()

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
CONDS = sorted(
    d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists() and args.set in d.name
)

if args.arm:
    store = {}
    for cond in CONDS:
        inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
        dbg: dict = {}
        cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
        calmod.calibrate(
            inp["payload"], ra, inp["strand_model"],
            np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
        )
        chain, cap, st, g = dbg["chain"], dbg["capture"], dbg["statics"], dbg["geometry"]
        Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
        G, R = Gp + Gn, Rp + Rn
        mass, _e = node_global_geometry(chain, g)
        rt, _ = _node_region_type(chain, ra)
        spl = (
            np.asarray(g.spliced_pos_left, float) + np.asarray(g.spliced_neg_left, float)
            + np.asarray(g.spliced_pos_right, float) + np.asarray(g.spliced_neg_right, float)
        )
        left, right = np.asarray(chain.left), np.asarray(chain.right)
        # is this node ADJACENT to a splice junction (i.e. inside the routing change's blast radius)?
        adj_j = np.zeros(mass.shape[0], bool)
        for nb in (left, right):
            ok = nb >= 0
            adj_j |= np.where(ok, spl[np.clip(nb, 0, None)] > _EPS, False)
        store[cond] = dict(
            fg=np.asarray(cap["f_g"]),
            fo=np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan),
            mass=mass, kind=np.asarray(chain.kind), rt=rt, spl=spl, adj_j=adj_j,
            tau0=np.asarray(cap["_tau0_lam"]),
            pg=np.asarray(cap["prec_g"]), pp=np.asarray(cap["prec_p"]), pn=np.asarray(cap["prec_n"]),
            fp=np.asarray(st.free_pos, bool), fn=np.asarray(st.free_neg, bool),
        )
        print(f"  {args.arm:>6} {cond}", flush=True)
    (OUT / f"{args.arm}.pkl").write_bytes(pickle.dumps(store))
    print(f"wrote {OUT / f'{args.arm}.pkl'}")

if args.report:
    A = pickle.loads((OUT / "base.pkl").read_bytes())
    B = pickle.loads((OUT / "route.pkl").read_bytes())
    tot_mass = 0.0
    dnum = 0.0
    rows = []
    for cond in A:
        a, b = A[cond], B[cond]
        ok = np.isfinite(a["fo"]) & (a["mass"] > _EPS)
        ea = np.abs(a["fg"] - a["fo"])
        eb = np.abs(b["fg"] - b["fo"])
        d = (eb - ea) * a["mass"]
        tot_mass += a["mass"][ok].sum()
        dnum += d[ok].sum()
        for i in np.where(ok)[0]:
            rows.append(
                (d[i], a["mass"][i], ea[i], eb[i], a["fg"][i], b["fg"][i], a["fo"][i],
                 int(a["kind"][i]), int(a["rt"][i]), bool(a["adj_j"][i]), a["tau0"][i],
                 a["spl"][i], int(a["fp"][i]) + int(a["fn"][i]), cond)
            )
    print(f"\nPAIRED DIFF over {len(A)} scenarios  (base -> route)")
    print(f"   mass-weighted mwae delta = {dnum/tot_mass:+.5f}\n")

    d = np.array([r[0] for r in rows])
    order = np.argsort(-d)
    # concentration: how much of the total damage comes from the worst k% of nodes?
    pos = d[d > 0].sum()
    dsort = np.sort(d)[::-1]
    cum = np.cumsum(dsort)
    for frac in (0.001, 0.01, 0.05, 0.10):
        k = max(1, int(frac * d.size))
        print(f"   worst {frac*100:>5.1f}% of nodes ({k:>5}) carry {100*cum[k-1]/max(pos,1e-12):>5.1f}% of the total damage")
    print(f"   nodes harmed: {(d>0).sum():>6}   helped: {(d<0).sum():>6}   unchanged: {(np.abs(d)<1e-12).sum():>6}")

    def klass(r):
        kind, rt, adj, nstr = r[7], r[8], r[9], r[12]
        base = "boundary" if kind != REGION else {0: "intergenic", 1: "intron", 2: "exon"}[rt]
        if kind != REGION:
            base += " (junction)" if r[11] > _EPS else " (seam)"
        elif adj:
            base += " adj-junction"
        return f"{base} [{nstr}-strand]"

    agg: dict[str, list] = {}
    for r in rows:
        agg.setdefault(klass(r), []).append(r)
    print(f"\n   {'class':<34} | {'n':>6} {'mass%':>6} | {'damage%':>8} | {'|err| base->route':>19} | {'mean Δf_g':>9}")
    print("   " + "-" * 96)
    for k in sorted(agg, key=lambda k: -sum(x[0] for x in agg[k])):
        v = agg[k]
        dm = sum(x[0] for x in v)
        m = np.array([x[1] for x in v])
        ea = np.average([x[2] for x in v], weights=m)
        eb = np.average([x[3] for x in v], weights=m)
        dfg = np.average([x[5] - x[4] for x in v], weights=m)
        print(f"   {k:<34} | {len(v):>6} {100*m.sum()/tot_mass:>5.1f}% | {100*dm/max(pos,1e-12):>7.1f}% |"
              f" {ea:>8.4f} -> {eb:>8.4f} | {dfg:>+9.4f}")

    print("\n   TOP 12 individual nodes by mass-weighted error increase:")
    print(f"   {'class':<34} {'mass':>9} {'oracle':>7} {'f_g base':>9} {'f_g route':>10} {'Δ(m·err)':>10}")
    for j in order[:12]:
        r = rows[j]
        print(f"   {klass(r):<34} {r[1]:>9.1f} {r[6]:>7.3f} {r[4]:>9.3f} {r[5]:>10.3f} {r[0]:>10.2f}")
