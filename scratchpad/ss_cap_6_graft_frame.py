"""SINGLE-STRAND x CAPTURE, step 6 — the graft's REQUIRED frame factor, and what predicts it.

Step 5: the delivered exon claim is  rho_g^src*E_g : (rho_nu + rho_mu)*E_r  (r and the pin cancel, fid 1e-15),
and rho_mu/rho_nu = 400..37000, so the GRAFT is the entire RNA claim. Off-capture it lands within 0.5 nats;
on-capture it overshoots by 1.2-3.3 nats of RNA.

Solve for the factor c that each edge would need on rho_mu to deliver the ORACLE composition:

    rho_nu + c*rho_mu = rho_g * E_g * (1-fo_g) / (E_r * fo_g)

and regress log c on the candidate predictors: log r (the S5.1 "graft outside the reframe" fix predicts
log c = -log r), and 0 (the shipped model predicts log c = 0).

    OMP_NUM_THREADS=1 python scratchpad/ss_cap_6_graft_frame.py
"""

from __future__ import annotations

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
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def q(x, w, p):
    o = np.argsort(x)
    cw = np.cumsum(w[o]) / w[o].sum()
    return x[o][np.searchsorted(cw, p)]


print(f"{'condition':<44}{'n':>6}{'med log c':>11}{'p25':>8}{'p75':>8}|{'med log r':>11}"
      f"{'corr':>7}{'slope':>7}|{'|logc|':>8}{'|logc+logr|':>12}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st, geo = dbg["chain"], dbg["capture"], dbg["statics"], dbg["geometry"]
    us = cap["_uni_static"]
    uni = cap["_uni"][-1]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    T = G + R
    fo_g = np.where(T > _EPS, G / np.maximum(T, _EPS), np.nan)
    E_g, E_r, M = us["E_g"], us["E_r"], us["M"]
    fp, fn, li, ri = us["fp"], us["fn"], us["left"], us["right"]
    n = len(M)
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    SP, SN = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    splf = {k: (np.where(SP[k] > _EPS, SP[k] / np.maximum(ESP[k], _EPS), 0.0),
                np.where(SN[k] > _EPS, SN[k] / np.maximum(ESP[k], _EPS), 0.0)) for k in (0, 1)}
    rt, _ = _node_region_type(chain, ra)
    isex = (np.asarray(chain.kind) == REGION) & (rt == 2)
    mass = np.asarray(cap["mass_global"])
    rho_l0, rho_r0, rho_node = us["rho_l0"], us["rho_r0"], us["rho_node0"]

    LC, LR, W = [], [], []
    for src_i, sf, relay, face in ((li, 1, ("fwd_g", "fwd_p", "fwd_n"), rho_r0),
                                   (ri, 0, ("bwd_g", "bwd_p", "bwd_n"), rho_l0)):
        s = np.clip(src_i, 0, n - 1)
        rg, rp, rn = (us[k][s] for k in relay)
        nu = np.where(fp, rp, 0.0) + np.where(fn, rn, 0.0)
        mu = np.where(fp, splf[sf][0][s], 0.0) + np.where(fn, splf[sf][1][s], 0.0)
        want = rg * E_g * (1.0 - fo_g) / np.maximum(E_r * fo_g, _EPS)
        c = (want - nu) / np.maximum(mu, _EPS)
        r = rho_node / np.maximum(face[s], _EPS)
        m = (isex & (src_i >= 0) & us["is_bnd"][s] & np.isfinite(fo_g) & (fo_g > 1e-6)
             & (fo_g < 1 - 1e-6) & (mass > _EPS) & (rg > _EPS) & (mu > _EPS) & (c > _EPS)
             & (face[s] > _EPS) & (rho_node > _EPS))
        LC.append(np.log(c[m]))
        LR.append(np.log(r[m]))
        W.append(mass[m])
    lc, lr, w = np.concatenate(LC), np.concatenate(LR), np.concatenate(W)
    if lc.size < 10:
        continue
    cw = np.cov(np.stack([lc, lr]), aweights=w)
    corr = cw[0, 1] / np.sqrt(cw[0, 0] * cw[1, 1])
    print(f"{cond[5:]:<44}{lc.size:>6}{q(lc, w, .5):>11.3f}{q(lc, w, .25):>8.3f}{q(lc, w, .75):>8.3f}|"
          f"{q(lr, w, .5):>11.3f}{corr:>7.3f}{cw[0, 1] / cw[1, 1]:>7.3f}|"
          f"{np.average(np.abs(lc), weights=w):>8.3f}{np.average(np.abs(lc + lr), weights=w):>12.3f}")
