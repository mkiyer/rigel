"""Check finding #8 on the DOC'S OWN quantity.

`weighted_rescale_design.md` §3.3 bounds a composition-only conservation violation by "x1.04 contained,
x1.50 at a boundary crossing".  Those two numbers are `composition_logvar`'s own docstring examples for
    coeff = (1/E_g - 1/E_r) / B            (exp(0.036) = 1.037,  exp(0.40) = 1.49)
i.e. d log rho_tot / d f_g at the node.  The report instead measures |log(E_g/E_r)| (the TOTAL variation of
log B over f_g in [0,1]) and gets 0.357 / 0.566-0.692.  Both are legitimate; this measures the DOC's coeff
directly on the suite so the two can be compared on the same footing.
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
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
]
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

print(f"{'cond':<34}{'node set':>10}{'n':>6}{'med coeff':>11}{'p90 coeff':>11}"
      f"{'med |log Eg/Er|':>17}{'exp(med coeff)':>16}")
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
    cap = dbg["capture"]
    us = cap["_uni_static"]
    M, E_g, E_r, is_bnd = us["M"], us["E_g"], us["E_r"], us["is_bnd"]
    fg = np.clip(np.asarray(cap["f_g"], float), 0.0, 1.0)
    ig, ir = 1.0 / np.maximum(E_g, _EPS), 1.0 / np.maximum(E_r, _EPS)
    B = fg * ig + (1.0 - fg) * ir
    coeff = np.abs((ig - ir) / np.maximum(B, _EPS))
    lever = np.abs(np.log(np.maximum(E_g, _EPS) / np.maximum(E_r, _EPS)))
    for nm, m in (("region", ~is_bnd & (M > _EPS)), ("boundary", is_bnd & (M > _EPS))):
        print(f"{cond[5:]:<34}{nm:>10}{int(m.sum()):>6}{np.median(coeff[m]):>11.3f}"
              f"{np.percentile(coeff[m], 90):>11.3f}{np.median(lever[m]):>17.3f}"
              f"{np.exp(np.median(coeff[m])):>16.3f}")
