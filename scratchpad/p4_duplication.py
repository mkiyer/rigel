"""P4 — IS `_far` A DUPLICATE? Compare the far level fused into message L's peel share against the RNA
claim the far node's OWN message delivers to the same boundary.

If the two numbers are the same quantity, `_far` adds no information — it only makes the destination count
the far node's evidence twice, and the per-component inverse-variance fuse then adds the two precisions as
if they were independent.

`cap["_lvl"]` holds 4 records per rho-iteration, appended (df=0,pos), (df=0,neg), (df=1,pos), (df=1,neg).

    OMP_NUM_THREADS=1 python scratchpad/p4_duplication.py [cond ...]
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
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna_none_ss_0.50_nrna_present_capture_off",
]

print(f"\n{'=' * 116}\nP4 — the far LEVEL vs the SAME NODE's own delivered message (peel edges, far = intron)"
      f"\n{'=' * 116}")
print(f"  {'condition':<44}{'edges':>7}{'nu_far':>9}{'msgRNA':>9}{'ratio':>8}"
      f"{'corr(log)':>10}{'med|dlog|':>10}{'p90':>8}")
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
    rec = cap["_lvl"][-4:]
    nuf, msg, wts = [], [], []
    for df, (rp, rn) in ((0, (rec[0], rec[1])), (1, (rec[2], rec[3]))):
        assert rp["df"] == df and rn["df"] == df
        src = li if df == 0 else ri
        far = ri if df == 0 else li
        peel = is_bnd & (src >= 0) & is_exon[np.clip(src, 0, n - 1)] & solv & (M > _EPS)
        peel = peel & (far >= 0) & (cls[np.clip(far, 0, n - 1)] == "intron") & (rp["pf"] + rn["pf"] > 0)
        if not peel.any():
            continue
        other = (uni["bp"] + uni["bn"]) if df == 0 else (uni["ap"] + uni["an"])
        nuf.append((rp["nu_f"] + rn["nu_f"])[peel])
        msg.append(other[peel])
        wts.append(M[peel])
    if not nuf:
        continue
    a, b, w = np.concatenate(nuf), np.concatenate(msg), np.concatenate(wts)
    k = (a > _EPS) & (b > _EPS)
    la, lb = np.log(a[k]), np.log(b[k])
    cc = float(np.corrcoef(la, lb)[0, 1]) if k.sum() > 3 else float("nan")
    ma, mb = np.average(a[k], weights=w[k]), np.average(b[k], weights=w[k])
    print(f"  {cond[5:]:<44}{int(k.sum()):>7}{ma:>9.3f}{mb:>9.3f}{ma / max(mb, _EPS):>8.2f}"
          f"{cc:>10.3f}{float(np.median(np.abs(la - lb))):>10.3f}"
          f"{float(np.quantile(np.abs(la - lb), 0.9)):>8.3f}")
print("\n  nu_far = the level `_far` puts into the peel share of the message arriving from the OTHER side")
print("  msgRNA = the total RNA density that same far node's OWN message delivers to the same boundary")
