"""SINGLE-STRAND x CAPTURE, step 5 — WHERE does the exon RNA over-claim bias come from?

Step 4: off-capture the exon RNA message is perfectly calibrated (z2 = 1.0, bias -0.04); on-capture it is
biased +0.5..+1.4 nats of RNA over-claim at an essentially correct precision (z2 = 5-7).

The delivered SHARE is  f_g : f_R = rho_g^src * E_g(dst) : (rho_nu^src + rho_mu) * E_r(dst)  -- the reframe r
and the pin factor k multiply all components alike and CANCEL. So the bias lives in the SOURCE's composition
plus the GRAFT. This splits it:

    lambda_delivered   = log(rho_g^src * E_g) - log((rho_nu^src + rho_mu) * E_r)
    lambda_nograft     = log(rho_g^src * E_g) - log( rho_nu^src              * E_r)
    d_graft            = lambda_delivered - lambda_nograft   = -log(1 + rho_mu/rho_nu)

    OMP_NUM_THREADS=1 python scratchpad/ss_cap_5_bias_source.py
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
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def lg(x):
    return np.log(np.maximum(np.asarray(x, np.float64), _EPS))


print(f"{'condition':<42}{'side':<6}{'n':>6}{'bias_del':>10}{'bias_nog':>10}{'d_graft':>9}"
      f"{'mu/nu':>9}{'r_bnd>ex':>10}{'srcRNAsh':>10}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st, geo = dbg["chain"], dbg["capture"], dbg["statics"], dbg["geometry"]
    uni, us = cap["_uni"][-1], cap["_uni_static"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    T = G + R
    fo_g = np.where(T > _EPS, G / np.maximum(T, _EPS), np.nan)
    lam_true = lg(fo_g) - lg(1.0 - fo_g)

    E_g, E_r, M = us["E_g"], us["E_r"], us["M"]
    fp, fn = us["fp"], us["fn"]
    li, ri = us["left"], us["right"]
    n = len(M)
    # per-face mature density  spl_*_f[face] = SP[face]/ESP[face]
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    SP = (us["SP_l"], us["SP_r"])
    SN = (us["SN_l"], us["SN_r"])
    splf = {k: (np.where(SP[k] > _EPS, SP[k] / np.maximum(ESP[k], _EPS), 0.0),
                np.where(SN[k] > _EPS, SN[k] / np.maximum(ESP[k], _EPS), 0.0)) for k in (0, 1)}
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    isex = (kind == REGION) & (rt == 2)
    mass = np.asarray(cap["mass_global"])

    for side, (src, valid_src, sf, relay, pub) in {
        "left": (np.clip(li, 0, n - 1), li >= 0, 1, ("fwd_g", "fwd_p", "fwd_n"), ("ag", "ap", "an")),
        "right": (np.clip(ri, 0, n - 1), ri >= 0, 0, ("bwd_g", "bwd_p", "bwd_n"), ("bg", "bp", "bn")),
    }.items():
        rg, rp, rn = (us[k][src] for k in relay)
        gp = np.where(us["is_bnd"][src], splf[sf][0][src], 0.0)
        gn = np.where(us["is_bnd"][src], splf[sf][1][src], 0.0)
        nu = np.where(fp, rp, 0.0) + np.where(fn, rn, 0.0)
        mu = np.where(fp, gp, 0.0) + np.where(fn, gn, 0.0)
        lam_del = lg(rg * E_g) - lg((nu + mu) * E_r)
        lam_nog = lg(rg * E_g) - lg(nu * E_r)
        # FIDELITY: the published post-transport densities must reproduce lam_del (r and k cancel)
        pg_, pp_, pn_ = (uni[k] for k in pub)
        lam_pub = lg(pg_ * E_g) - lg((pp_ + pn_) * E_r)
        m = (isex & valid_src & us["is_bnd"][src] & np.isfinite(lam_true) & (mass > _EPS)
             & (rg > _EPS) & (nu > _EPS) & (pg_ > _EPS) & ((pp_ + pn_) > _EPS))
        if m.sum() < 5:
            continue
        w = mass[m]
        fid = np.max(np.abs(lam_del[m] - lam_pub[m]))
        r_be = np.where(us["rho_l0"][src] > _EPS, us["rho_node0"] / np.maximum(us["rho_r0"][src], _EPS), 1.0)
        print(f"{cond[5:] if side == 'left' else '':<42}{side:<6}{int(m.sum()):>6}"
              f"{np.average((lam_del - lam_true)[m], weights=w):>+10.3f}"
              f"{np.average((lam_nog - lam_true)[m], weights=w):>+10.3f}"
              f"{np.average((lam_del - lam_nog)[m], weights=w):>+9.3f}"
              f"{np.average((mu / np.maximum(nu, _EPS))[m], weights=w):>9.2f}"
              f"{np.median(r_be[m]):>10.2f}{np.average((mu + nu)[m] * E_r[m] / np.maximum((mu + nu)[m] * E_r[m] + rg[m] * E_g[m], _EPS), weights=w):>10.3f}"
              f"   fid={fid:.1e}")
    print()
