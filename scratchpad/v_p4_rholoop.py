"""ADVERSARIAL ADDITION to the P4 audit: `_far` is NOT the only route by which the RIGHT neighbour's data
reaches the LEFT message. `_RHO_ITERS = 2` and `f_cur = dc_fin.gdna_frac` (the FUSED belief) — so iteration
2's reframe ratio r = dst_face[k]/src_face[src] is built from node k's OWN posterior, which already contains
m_{R->k}. Measure how live that loop is.
"""
from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

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
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
]
print(f"  {'condition':<44}{'nodes':>8}{'med|dlog rho_lf|':>18}{'p90':>9}{'frac>1%':>9}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    u = dbg["capture"]["_uni"]
    a, b = u[0]["rho_lf"], u[1]["rho_lf"]
    k = (a > _EPS) & (b > _EPS)
    dl = np.abs(np.log(b[k]) - np.log(a[k]))
    print(f"  {cond[5:]:<44}{int(k.sum()):>8}{np.median(dl):>18.4f}"
          f"{np.quantile(dl, 0.9):>9.4f}{np.mean(dl > 0.01):>9.1%}")
print("\n  rho_lf is the LEFT message's dst_face frame. Iteration 1 uses the INPUT belief; iteration 2 uses")
print("  the belief FUSED FROM BOTH MESSAGES. Any non-zero column means the loop is live and the left")
print("  message's reframe depends on the right neighbour's data — a third-node dependence `_far`'s")
print("  deletion does NOT remove.")
