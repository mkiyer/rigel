"""LEVEL PROVENANCE — which of the three estimators actually sets the seam's RNA level, and how wrong is it.

Reads the solver's own `_lvl` capture (published from `_peel_share` in the vectorized combine, so these are
exactly the numbers the shipped answer was built from — no offline re-derivation) and scores the fused level
and each component against the ORACLE RNA density at the boundary.

    OMP_NUM_THREADS=1 python scratchpad/w4_provenance.py [cond ...]
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
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna1_ss_0.50_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

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
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION
    n = kind.shape[0]

    def pool(k):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    Ru = pool("mat_uns_pos") + pool("nas_uns_pos") + pool("mat_uns_neg") + pool("nas_uns_neg")
    G = pool("gdna_pos") + pool("gdna_neg")
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    nu_true = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), 0.0)
    rt, _ = _node_region_type(chain, ra)
    CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array([CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(n)])
    li, ri = us["left"], us["right"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    solv = np.asarray(cap["solvable"], bool)

    # the LAST rho-iteration's records: 2 directions x 2 strands = 4, appended in order
    recs = cap["_lvl"][-4:]
    print(f"\n{'=' * 122}\n{cond[5:]}   (oracle f_g at boundaries "
          f"{np.average(np.where(G + Ru > _EPS, G / np.maximum(G + Ru, _EPS), 0)[is_bnd & solv], weights=M[is_bnd & solv]):.3f})"
          f"\n{'=' * 122}")
    print(f"  {'stratum':<26}{'n':>6}{'nu_true':>9}{'nu_fused':>10}{'nu_own':>9}{'nu_far':>9}{'nu_mass':>9}"
          f"{'rho_mu':>9}{'w':>7}{'w_true':>8}{'phi':>7}{'sd_g':>7}{'sd_m':>7}{'sd_f':>7}"
          f"   weight own/far/mass")
    for lab, msel in (
        ("ALL peel edges", None),
        ("  far = intron", "intron"),
        ("  far = exon", "exon"),
    ):
        num = {k: 0.0 for k in ("t", "f", "o", "fa", "m", "mu", "w", "wt", "phi", "sg", "sm", "sf")}
        wo = wf = wm = wsum = 0.0
        cnt = 0
        for R in recs:
            df = R["df"]
            src = li if df == 0 else ri
            far = ri if df == 0 else li
            s = np.clip(src, 0, n - 1)
            peel = is_bnd & is_exon[s] & (src >= 0) & solv & (M > _EPS)
            if msel is not None:
                cf = np.array([cls[j] if j >= 0 else "edge" for j in far])
                peel = peel & (cf == msel)
            peel = peel & (R["po"] + R["pf"] + R["pm"] > 0)
            if not peel.any():
                continue
            w = M[peel]
            cnt += int(peel.sum())
            wsum += float(w.sum())
            num["t"] += float(np.sum(nu_true[peel] * w))
            num["f"] += float(np.sum(R["nu"][peel] * w))
            num["o"] += float(np.sum(R["nu_o"][peel] * w))
            num["fa"] += float(np.sum(R["nu_f"][peel] * w))
            num["m"] += float(np.sum(R["nu_m"][peel] * w))
            num["mu"] += float(np.sum(R["mu"][peel] * w))
            num["w"] += float(np.sum(R["w"][peel] * w))
            wt = np.where(nu_true[peel] + R["mu"][peel] > _EPS,
                          nu_true[peel] / np.maximum(nu_true[peel] + R["mu"][peel], _EPS), 1.0)
            num["wt"] += float(np.sum(wt * w))
            num["phi"] += float(np.sum(np.minimum(R["phi"][peel], 3.0) * w))
            num["sg"] += float(np.sum(np.sqrt(np.minimum(R["v_g"][peel], 100.0)) * w))
            num["sm"] += float(np.sum(np.sqrt(np.minimum(R["vl_m"][peel], 1e9)) * w))
            _vf = np.where(R["pf"][peel] > 0, 1.0 / np.maximum(R["pf"][peel], 1e-30), 0.0)
            num["sf"] += float(np.sum(np.sqrt(np.minimum(_vf, 1e9)) * w))
            pt = R["po"][peel] + R["pf"][peel] + R["pm"][peel]
            wo += float(np.sum(R["po"][peel] / pt * w))
            wf += float(np.sum(R["pf"][peel] / pt * w))
            wm += float(np.sum(R["pm"][peel] / pt * w))
        if wsum <= 0:
            continue
        g = lambda k: num[k] / wsum  # noqa: E731
        print(f"  {lab:<26}{cnt:>6}{g('t'):>9.3f}{g('f'):>10.3f}{g('o'):>9.3f}{g('fa'):>9.3f}"
              f"{g('m'):>9.3f}{g('mu'):>9.3f}{g('w'):>7.3f}{g('wt'):>8.3f}"
              f"{g('phi'):>7.3f}{g('sg'):>7.3f}{g('sm'):>7.3f}{g('sf'):>7.3f}"
              f"   {wo / wsum:.0%} / {wf / wsum:.0%} / {wm / wsum:.0%}")
