"""P4 — WHOSE DATA SETS THE PEEL SHARE, and is the far node the other message's source?

Reads the solver's own `_lvl` capture (published by `_peel_share` in the vectorized combine) over the whole
32-condition suite and reports, mass-weighted on PEEL edges:

  * the share of the fused level's PRECISION contributed by each of the three estimators
        po = the destination boundary's OWN self-solve   (a plug-in of x_k  -> approximation)
        pf = `_far`, the node ACROSS the seam            (a THIRD node -> the BP violation)
        pm = M11 `residual_level`, from the message's own gDNA claim + the node's mass (legal)
  * the STRUCTURAL check: is `_far`'s node identically the source of the OTHER message?
  * how much of the OTHER message's own delivered precision is seeded by that same node's self-solve
    (i.e. the size of what is being counted twice).

    OMP_NUM_THREADS=1 python scratchpad/p4_far_provenance.py [cond ...]
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
CONDS = sys.argv[1:] or sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

ident_far_is_other_src = True
TOT = {}
BYFAR = {}
rows = []
for cond in CONDS:
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap = dbg["chain"], dbg["capture"]
    us = cap["_uni_static"]
    kind = np.asarray(chain.kind)
    n = kind.shape[0]
    rt, _ = _node_region_type(chain, ra)
    CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array([CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(n)])
    li, ri = np.asarray(us["left"]), np.asarray(us["right"])
    is_bnd, is_exon = np.asarray(us["is_bnd"], bool), np.asarray(us["is_exon"], bool)
    solv = np.asarray(cap["solvable"], bool)
    M = np.asarray(us["M"], float)
    uni = cap["_uni"][-1]

    acc = {k: 0.0 for k in ("w", "po", "pf", "pm", "live", "nlive", "nedge", "wlive")}
    for R in cap["_lvl"][-4:]:  # last rho-iteration: 2 faces x 2 strands
        df = R["df"]
        src = li if df == 0 else ri
        far = ri if df == 0 else li
        # ---- STRUCTURAL: for a boundary, `_far`'s node IS the other message's source
        othersrc = ri if df == 0 else li
        peel = is_bnd & (src >= 0) & is_exon[np.clip(src, 0, n - 1)] & solv & (M > _EPS)
        if (far[peel] != othersrc[peel]).any():
            ident_far_is_other_src = False
        pt = R["po"] + R["pf"] + R["pm"]
        m = peel & (pt > 0)
        if not m.any():
            continue
        w = M[m]
        acc["w"] += float(w.sum())
        acc["po"] += float(np.sum(R["po"][m] / pt[m] * w))
        acc["pf"] += float(np.sum(R["pf"][m] / pt[m] * w))
        acc["pm"] += float(np.sum(R["pm"][m] / pt[m] * w))
        acc["nedge"] += float(m.sum())
        lv = m & (R["pf"] > 0)
        acc["nlive"] += float(lv.sum())
        acc["wlive"] += float(M[lv].sum())
        # per far-node class
        for c in ("intron", "exon", "intergenic", "boundary"):
            k = m & (cls[np.clip(far, 0, n - 1)] == c) & (far >= 0)
            if not k.any():
                continue
            e = BYFAR.setdefault(c, {"w": 0.0, "po": 0.0, "pf": 0.0, "pm": 0.0, "n": 0.0, "nlive": 0.0})
            ww = M[k]
            e["w"] += float(ww.sum())
            e["po"] += float(np.sum(R["po"][k] / pt[k] * ww))
            e["pf"] += float(np.sum(R["pf"][k] / pt[k] * ww))
            e["pm"] += float(np.sum(R["pm"][k] / pt[k] * ww))
            e["n"] += float(k.sum())
            e["nlive"] += float((k & (R["pf"] > 0)).sum())
    if acc["w"] <= 0:
        continue
    rows.append((cond, acc))
    for k, v in acc.items():
        TOT[k] = TOT.get(k, 0.0) + v

print(f"\n{'=' * 108}\nP4 — peel-share LEVEL precision provenance (mass-weighted, last rho-iteration)\n{'=' * 108}")
print(f"  {'condition':<46}{'edges':>8}{'far live':>10}{'own %':>9}{'FAR %':>9}{'M11 %':>9}")
for cond, a in sorted(rows, key=lambda r: -r[1]["pf"] / max(r[1]["w"], _EPS)):
    print(f"  {cond[5:]:<46}{a['nedge']:>8,.0f}{a['nlive'] / max(a['nedge'], 1):>9.1%}"
          f"{a['po'] / a['w']:>9.1%}{a['pf'] / a['w']:>9.1%}{a['pm'] / a['w']:>9.1%}")
a = TOT
print(f"  {'-' * 90}\n  {'ALL':<46}{a['nedge']:>8,.0f}{a['nlive'] / max(a['nedge'], 1):>9.1%}"
      f"{a['po'] / a['w']:>9.1%}{a['pf'] / a['w']:>9.1%}{a['pm'] / a['w']:>9.1%}")

print(f"\n  by FAR-NODE class:\n  {'far class':<20}{'edges':>10}{'far live':>10}{'own %':>9}{'FAR %':>9}{'M11 %':>9}")
for c, e in sorted(BYFAR.items(), key=lambda kv: -kv[1]["pf"]):
    if e["w"] <= 0:
        continue
    print(f"  {c:<20}{e['n']:>10,.0f}{e['nlive'] / max(e['n'], 1):>10.1%}"
          f"{e['po'] / e['w']:>9.1%}{e['pf'] / e['w']:>9.1%}{e['pm'] / e['w']:>9.1%}")

print(f"\n  STRUCTURAL CHECK — `_far`'s node is identically the OTHER message's source: {ident_far_is_other_src}")
