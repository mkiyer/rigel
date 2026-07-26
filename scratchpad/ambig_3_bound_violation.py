"""AMBIG study, step 3 — do incoming messages VIOLATE the node's own structural bound?

The strand mixture is exactly  p = 1/2 + (kappa-1/2)*d  with d = f_+ - f_-  (gDNA is 50/50, so f_g drops out
of the mean entirely -- that is the "strand constrains the tilt, not f_g" result). So d is DIRECTLY OBSERVED:

    d_hat = (p_obs - 1/2) / (kappa - 1/2)        p_obs = u_pos / (u_pos + u_neg)

and the simplex gives |d| <= f_+ + f_- = 1 - f_g, i.e. the HARD structural bound

    f_g  <=  B = 1 - |d_hat|          (equality iff the minority strand carries no RNA)

This is a bound, not an estimate -- reading f_g = B would be a PRIOR ("one strand is silent"). But f_g > B is
IMPOSSIBLE, up to the sampling noise of d_hat. `node_init` nevertheless registers tau_own = 0 for every AMBIG
node ⇒ v_own = inf ⇒ M7's DL term is INERT there, so a message asserting f_g above the bound is accepted at
full precision. This measures how often that happens and what it costs.

    OMP_NUM_THREADS=1 python scratchpad/ambig_3_bound_violation.py
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

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def mw(x, w):
    return np.average(x, weights=w) if np.sum(w) > 0 else np.nan


print("AMBIG nodes.  B = 1-|d_hat| is a HARD upper bound on f_g.  sd_B from the strand sampling noise.")
print(f"{'condition':<42}{'n':>5}{'B-err':>8}{'self':>8}{'ship':>8}|"
      f"{'msg>B':>8}{'mass>B':>8}{'z_viol':>8}{'errViol':>9}{'errOK':>8}|{'shipB':>8}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    res = calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    uni = cap["_uni"][-1]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    T = G + R
    fo = np.where(T > _EPS, G / np.maximum(T, _EPS), np.nan)
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    amb = fp & fn
    up = np.asarray(st.u_pos, float)
    un = np.asarray(st.u_neg, float)
    n = up + un
    kap = float(res.rna_sense_frac)
    od_r = float(res.rna_strand_overdispersion)
    n_eff = n / (1.0 + np.maximum(n - 1.0, 0.0) * od_r)
    p_obs = np.where(n > 0, up / np.maximum(n, _EPS), 0.5)
    denom = kap - 0.5
    d_hat = (p_obs - 0.5) / denom if abs(denom) > _EPS else np.zeros_like(p_obs)
    B = np.clip(1.0 - np.abs(d_hat), _EPS, 1.0)
    # sampling sd of the bound: Var(d_hat) = p(1-p)/(N_eff (kappa-1/2)^2)
    var_d = np.where(n_eff > 0, p_obs * (1 - p_obs) / np.maximum(n_eff, _EPS), np.inf) / max(
        denom * denom, _EPS
    )
    sd_B = np.sqrt(np.maximum(var_d, 0.0))
    msg_fg = 1.0 / (1.0 + np.exp(-uni["lam_msg"]))  # the delivered lambda message as a fraction
    ship = np.asarray(cap["f_g"])
    self_fg = np.asarray(cap["fg_loc"])
    mass = np.asarray(cap["mass_global"])
    ok = (np.isfinite(fo) & (mass > _EPS) & np.asarray(cap["solvable"], bool) & amb
          & (n > 0) & (uni["c_tau"] > _EPS) & np.isfinite(sd_B))
    if ok.sum() < 5:
        print(f"{cond[5:]:<42}{'(no scoreable AMBIG lambda messages)':>60}")
        continue
    w = mass[ok]
    viol = msg_fg[ok] > B[ok]
    z = (msg_fg[ok] - B[ok]) / np.maximum(sd_B[ok], _EPS)
    err = np.abs(ship - fo)[ok]
    print(f"{cond[5:]:<42}{int(ok.sum()):>5}{mw(np.abs(B - fo)[ok], w):>8.4f}"
          f"{mw(np.abs(self_fg - fo)[ok], w):>8.4f}{mw(err, w):>8.4f}|"
          f"{viol.mean():>8.1%}{mw(viol.astype(float), w):>8.1%}"
          f"{mw(np.where(viol, z, 0.0), w):>8.1f}"
          f"{mw(err[viol], w[viol]) if viol.any() else np.nan:>9.4f}"
          f"{mw(err[~viol], w[~viol]) if (~viol).any() else np.nan:>8.4f}|"
          f"{mw(np.abs(np.minimum(ship, B) - fo)[ok], w):>8.4f}")
