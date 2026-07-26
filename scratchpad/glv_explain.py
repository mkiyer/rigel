"""What exactly is P1d's fitted scalar? Dump the fit, its population, and its effect.

    OMP_NUM_THREADS=1 python scratchpad/glv_explain.py [cond ...]
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
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "gdna_none_ss_0.50_nrna_present_capture_off",
]
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

print("THE FIT:  omega = max(0, E[d^2] - E[noise]) / 2      d = log(flux_leftseam / flux_rightseam)")
print("one scalar per STRAND per library; the first fit shown is the relay's (the input-belief frame).\n")
hdr = (f"{'condition':<40}{'str':>4}{'pairs':>7}{'E[d^2]':>9}{'E[noise]':>10}{'OMEGA':>9}"
       f"{'sd(d)':>8}{'med|d|':>8}{'p90|d|':>8}{'exon%':>7}{'intr%':>7}")
print(hdr)
print("-" * len(hdr))
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
    rt, _ = _node_region_type(chain, ra)
    is_reg = np.asarray(chain.kind) == REGION
    fits = cap.get("_glv", [])
    for f in fits[:2]:  # the relay's two strands
        ok, d = f["ok"], f["d"]
        cls = rt[ok]
        n = max(int(ok.sum()), 1)
        ex = float(np.mean(cls == 2)) * 100 if ok.any() else 0.0   # 0=intergenic 1=intron 2=exon
        it = float(np.mean(cls == 1)) * 100 if ok.any() else 0.0
        print(f"{cond[5:]:<40}{'+-'[f['strand']]:>4}{f['n_pairs']:>7}{f['Ed2']:>9.4f}"
              f"{f['Enoise']:>10.4f}{f['omega']:>9.4f}"
              f"{np.std(d[ok]) if ok.any() else 0:>8.3f}"
              f"{np.median(np.abs(d[ok])) if ok.any() else 0:>8.3f}"
              f"{np.percentile(np.abs(d[ok]), 90) if ok.any() else 0:>8.3f}"
              f"{ex:>7.1f}{it:>7.1f}"
              f"{(np.std(d[ok]**2)/np.sqrt(n)/2 if ok.any() else 0):>9.3f}")
    # how many fits happen per sweep, and how stable they are
    om = [f["omega"] for f in fits if f["strand"] == 0]
    print(f"{'':<40}{'':>4}  fits/sweep(+strand) = {len(om)}   values = "
          f"{['%.4f' % v for v in om]}   regions with 2 live seams / all regions = "
          f"{int(fits[0]['ok'].sum())}/{int(is_reg.sum())}")
