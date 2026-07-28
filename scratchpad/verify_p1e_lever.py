"""Finding #8, checked at the right granularity.

`p1e_1_delta.py` prints one `bound` per (cond, EDGE CLASS) where the edge classes are
    plain = valid & ~graft & ~peel
which mixes REGION and BOUNDARY destinations in one row.  The report quotes its 0.357 as
"median |log(E_g/E_r)| ... 0.357 on region nodes (x1.43)".  This splits the same statistic by destination
kind, exactly as `p1e_2_audit.py` splits the delta census, and also reports the DOC's own quantity
coeff = |(1/E_g - 1/E_r)/B| on the same populations.
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

print(f"{'cond':<34}{'edge':>7}{'dst':>5}{'n':>7}{'med|d|':>9}{'med lever':>11}"
      f"{'med coeff':>11}{'x lever':>9}{'x coeff':>9}{'>lever':>8}")
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
    M, E_g, E_r, is_bnd, is_exon = us["M"], us["E_g"], us["E_r"], us["is_bnd"], us["is_exon"]
    og, op, on = us["og"], us["op"], us["on"]
    fg = np.clip(np.asarray(cap["f_g"], float), 0.0, 1.0)
    ig, ir = 1.0 / np.maximum(E_g, _EPS), 1.0 / np.maximum(E_r, _EPS)
    B = fg * ig + (1.0 - fg) * ir
    coeff = np.abs((ig - ir) / np.maximum(B, _EPS))
    lever = np.abs(np.log(np.maximum(E_g, _EPS) / np.maximum(E_r, _EPS)))
    for ent in cap["_pin"][-2:]:
        src, valid = ent["src"], ent["valid"]
        p3 = np.stack([ent["tpg"], ent["tpp"], ent["tpn"]], axis=1)
        sup = p3 > 0.0
        d3 = np.stack([ent["tg"], ent["tp"], ent["tn"]], axis=1)
        o3 = np.stack([og, op, on], axis=1)
        E3 = np.stack([E_g, E_r, E_r], axis=1)
        S = (np.where(sup, d3, o3) * E3).sum(axis=1)
        ok = valid & (S > _EPS) & (M > _EPS)
        delta = np.log(np.maximum(M, _EPS) / np.maximum(S, _EPS))
        graft = ent["graft"]
        peel = is_bnd & is_exon[src] & valid
        plain = valid & ~graft & ~peel
        for nm, msk in (("plain", plain), ("graft", graft), ("peel", peel)):
            for dk, dm in (("reg", ~is_bnd), ("bnd", is_bnd)):
                m = msk & dm & ok
                if m.sum() < 20:
                    continue
                print(f"{cond[5:]:<34}{nm:>7}{dk:>5}{int(m.sum()):>7}"
                      f"{np.median(np.abs(delta[m])):>9.3f}{np.median(lever[m]):>11.3f}"
                      f"{np.median(coeff[m]):>11.3f}{np.exp(np.median(lever[m])):>9.3f}"
                      f"{np.exp(np.median(coeff[m])):>9.3f}"
                      f"{np.mean(np.abs(delta[m]) > lever[m]):>8.1%}")
