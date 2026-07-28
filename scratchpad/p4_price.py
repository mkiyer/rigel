"""P4 — WHERE the deletion's price comes from, and whether it is legally recoverable.

On the two conditions that pay for removing `_far`, look at the boundaries whose far node is an intron and
ask why the intron's OWN message cannot do the job `_far` was doing: compare the two messages' RNA claims and
the precisions they arrive with.

    OMP_NUM_THREADS=1 python scratchpad/p4_price.py [cond ...]
"""

from __future__ import annotations

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
_EPS = 1e-9

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
]


def run(cond, p4):
    os.environ["RIGEL_P4"] = p4
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    os.environ.pop("RIGEL_P4", None)
    return inp, dbg


for cond in CONDS:
    inp, dbg = run(cond, "")
    _, dbg2 = run(cond, "m11")
    chain, cap, cap2 = dbg["chain"], dbg["capture"], dbg2["capture"]
    us, uni = cap["_uni_static"], cap["_uni"][-1]
    kind = np.asarray(chain.kind)
    n = kind.shape[0]
    rt, _ = _node_region_type(chain, ra)
    CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array([CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(n)])
    li, ri = np.asarray(us["left"]), np.asarray(us["right"])
    is_bnd, is_exon = np.asarray(us["is_bnd"], bool), np.asarray(us["is_exon"], bool)
    solv = np.asarray(cap["solvable"], bool)
    M = np.asarray(us["M"], float)
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    fh, fm = np.asarray(cap["f_g"]), np.asarray(cap2["f_g"])
    # boundaries fed by an EXON on one side and an INTRON on the other — `_far`'s working population
    lc, rc = cls[np.clip(li, 0, n - 1)], cls[np.clip(ri, 0, n - 1)]
    tgt = (
        is_bnd & solv & (M > _EPS) & np.isfinite(fo) & (li >= 0) & (ri >= 0)
        & (((lc == "exon") & (rc == "intron")) | ((lc == "intron") & (rc == "exon")))
    )
    exl = tgt & (lc == "exon")  # the exon message arrives from the left
    aR, bR = uni["ap"] + uni["an"], uni["bp"] + uni["bn"]
    # exon-side / intron-side claim + its mode-fusion precision (pos arm; unstranded ⇒ representative)
    exR = np.where(exl, aR, bR)
    inR = np.where(exl, bR, aR)
    exP = np.where(exl, uni["app"], uni["bpp"])
    inP = np.where(exl, uni["bpp"], uni["app"])
    w = M[tgt]
    g = lambda x: float(np.average(x[tgt], weights=w))  # noqa: E731
    print(f"\n{'=' * 100}\n{cond[5:]}   exon|B|intron boundaries: {int(tgt.sum())}  (mass {w.sum():,.0f})")
    print(f"{'=' * 100}")
    print(f"  oracle f_g {g(fo):.3f}   HEAD f_g {g(fh):.3f} (|err| {g(np.abs(fh - fo)):.3f})"
          f"   M11-only f_g {g(fm):.3f} (|err| {g(np.abs(fm - fo)):.3f})")
    print(f"  EXON-side msg   RNA {g(exR):>9.4f}   mode-prec {g(exP):>10.4f}")
    print(f"  INTRON-side msg RNA {g(inR):>9.4f}   mode-prec {g(inP):>10.4f}"
          f"   → intron share of the RNA fuse weight: {g(inP) / max(g(inP) + g(exP), _EPS):.1%}")
    print(f"  fused RNA {g(uni['cp'] + uni['cn']):>9.4f}     "
          f"(a linear-space precision-weighted mean of the two, so a confident NEAR-ZERO claim can only")
    print("   pull the level down if its precision is large — which is the arm `_far` was supplying.)")
