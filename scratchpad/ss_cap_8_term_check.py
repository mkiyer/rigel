"""SINGLE-STRAND x CAPTURE, step 8 — is w_mu^2*(log r)^2 the right size for the graft's frame mis-lift?

DERIVATION. On a GRAFT edge the delivered share is  rho_g^src*E_g : (rho_nu^src + rho_mu)*E_r, and r cancels
(step 5, fid 1e-15). M5 sets sigma^2_transfer = 0 there because "r is common-mode across {g,R}". That is true
for the IMPUTED rho_nu, which travels with rho_g from the same source -- but rho_mu is an ABSOLUTE spliced
measurement whose molecules are anchored in the flanking EXONS, so it is already in the DESTINATION's frame
(step 7: rho_R(exon)/rho_spl(bnd) = 1.02-1.86, capture-invariant). It has no matched gDNA partner to cancel r
against -- exactly M5's PEEL/partial case, where sigma^2_transfer is load-bearing.

So the grafted share carries the full frame step as a systematic error, of magnitude |log r_true|, and the
model's own estimate of that step is log r. By M2's delta method on the sum the share weight is w_mu^2:

    Var(log rho_R^msg) += w_mu^2 * (log r)^2 ,      w_mu = rho_mu / (rho_nu + rho_mu)

This checks the SIZE against the realized error: z2 = E[(mode-truth)^2]/E[v] per message, before and after.

    OMP_NUM_THREADS=1 python scratchpad/ss_cap_8_term_check.py
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
    "gdna_none_ss_0.50_nrna_none_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

print("GRAFT edges into exons. z2 = E[(delivered lambda - true lambda)^2] / E[v]  (1.0 = honest).")
print(f"{'condition':<44}{'n':>5}{'bias':>8}{'|bias|':>8}|{'v_now':>8}{'z2_now':>8}|"
      f"{'v_new':>8}{'z2_new':>8}|{'w_mu':>7}{'med|logr|':>10}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, geo = dbg["chain"], dbg["capture"], dbg["geometry"]
    us = cap["_uni_static"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    T = G + R
    fo_g = np.where(T > _EPS, G / np.maximum(T, _EPS), np.nan)
    lam_true = np.log(np.maximum(fo_g, _EPS)) - np.log(np.maximum(1 - fo_g, _EPS))
    E_g, E_r = us["E_g"], us["E_r"]
    fp, fn, li, ri = us["fp"], us["fn"], us["left"], us["right"]
    n = len(E_g)
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    SP, SN = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    splf = {k: (np.where(SP[k] > _EPS, SP[k] / np.maximum(ESP[k], _EPS), 0.0),
                np.where(SN[k] > _EPS, SN[k] / np.maximum(ESP[k], _EPS), 0.0)) for k in (0, 1)}
    rt, _ = _node_region_type(chain, ra)
    isex = (np.asarray(chain.kind) == REGION) & (rt == 2)
    mass = np.asarray(cap["mass_global"])
    rho_l0, rho_r0, rho_node = us["rho_l0"], us["rho_r0"], us["rho_node0"]

    A = {k: [] for k in ("d", "vn", "vw", "wm", "lr", "w")}
    for src_i, sf, relay, mrelay, face in (
        (li, 1, ("fwd_g", "fwd_p", "fwd_n"), ("fwd_mp", "fwd_mn"), rho_r0),
        (ri, 0, ("bwd_g", "bwd_p", "bwd_n"), ("bwd_mp", "bwd_mn"), rho_l0),
    ):
        s = np.clip(src_i, 0, n - 1)
        rg, rp, rn = (us[k][s] for k in relay)
        nu = np.where(fp, rp, 0.0) + np.where(fn, rn, 0.0)
        mu = np.where(fp, splf[sf][0][s], 0.0) + np.where(fn, splf[sf][1][s], 0.0)
        spl_cnt = np.where(fp, SP[sf][s], 0.0) + np.where(fn, SN[sf][s], 0.0)
        lam_del = np.log(np.maximum(rg * E_g, _EPS)) - np.log(np.maximum((nu + mu) * E_r, _EPS))
        r = rho_node / np.maximum(face[s], _EPS)
        w_mu = mu / np.maximum(nu + mu, _EPS)
        # the graft's shipped RNA precision is the spliced COUNT (s2t = 0 on a graft edge)
        v_now = 1.0 / np.maximum(spl_cnt, _EPS)
        v_new = v_now + w_mu * w_mu * np.log(np.maximum(r, _EPS)) ** 2
        m = (isex & (src_i >= 0) & us["is_bnd"][s] & np.isfinite(lam_true) & (mass > _EPS)
             & (rg > _EPS) & (mu > _EPS) & (spl_cnt > _EPS) & (face[s] > _EPS) & (rho_node > _EPS))
        A["d"].append((lam_del - lam_true)[m]); A["vn"].append(v_now[m]); A["vw"].append(v_new[m])
        A["wm"].append(w_mu[m]); A["lr"].append(np.abs(np.log(np.maximum(r, _EPS)))[m])
        A["w"].append(mass[m])
    a = {k: np.concatenate(v) for k, v in A.items()}
    w = a["w"]
    if w.size < 5:
        continue
    d2 = np.average(a["d"] ** 2, weights=w)
    print(f"{cond[5:]:<44}{w.size:>5}{np.average(a['d'], weights=w):>+8.3f}"
          f"{np.average(np.abs(a['d']), weights=w):>8.3f}|"
          f"{np.average(a['vn'], weights=w):>8.4f}{d2 / np.average(a['vn'], weights=w):>8.1f}|"
          f"{np.average(a['vw'], weights=w):>8.4f}{d2 / np.average(a['vw'], weights=w):>8.1f}|"
          f"{np.average(a['wm'], weights=w):>7.3f}{np.median(a['lr']):>10.3f}")
