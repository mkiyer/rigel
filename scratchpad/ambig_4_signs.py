"""AMBIG study, step 4 — SIGNS: which way do messages push, and where does the truth sit in [0, B]?

Step 3 established B = 1-|d_hat| is a hard structural bound and an excellent predictor (err 0.010-0.021 under
capture), but that message>B is NOT the damage marker. So the damage must be directional. This measures the
signed errors of the self-solve, the delivered lambda message, and the shipped answer, plus where the truth
sits inside the node's own admissible interval [0, B].

    OMP_NUM_THREADS=1 python scratchpad/ambig_4_signs.py
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
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def mw(x, w):
    return np.average(x, weights=w) if np.sum(w) > 0 else np.nan


print("AMBIG nodes, mass-weighted SIGNED errors (positive = too much gDNA).  frac = (oracle f_g)/B.")
print(f"{'condition':<42}{'kap':>6}{'oracle':>8}{'B':>7}{'self':>8}{'ship':>8}|"
      f"{'sSelf':>8}{'sMsg':>8}{'sShip':>8}|{'oracle/B':>9}{'self/B':>8}{'c_tau':>8}{'cm_g':>8}")
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
    up, un = np.asarray(st.u_pos, float), np.asarray(st.u_neg, float)
    n = up + un
    kap = float(res.rna_sense_frac)
    p_obs = np.where(n > 0, up / np.maximum(n, _EPS), 0.5)
    den = kap - 0.5
    d_hat = (p_obs - 0.5) / den if abs(den) > _EPS else np.zeros_like(p_obs)
    B = np.clip(1.0 - np.abs(d_hat), _EPS, 1.0)
    msg_fg = 1.0 / (1.0 + np.exp(-uni["lam_msg"]))
    ship, self_fg = np.asarray(cap["f_g"]), np.asarray(cap["fg_loc"])
    mass = np.asarray(cap["mass_global"])
    ok = np.isfinite(fo) & (mass > _EPS) & np.asarray(cap["solvable"], bool) & amb & (n > 0)
    w = mass[ok]
    live = ok & (uni["c_tau"] > _EPS)
    wl = mass[live]
    print(f"{cond[5:]:<42}{kap:>6.2f}{mw(fo[ok], w):>8.4f}{mw(B[ok], w):>7.4f}"
          f"{mw(np.abs(self_fg - fo)[ok], w):>8.4f}{mw(np.abs(ship - fo)[ok], w):>8.4f}|"
          f"{mw((self_fg - fo)[ok], w):>+8.4f}"
          f"{mw((msg_fg - fo)[live], wl) if live.any() else np.nan:>+8.4f}"
          f"{mw((ship - fo)[ok], w):>+8.4f}|"
          f"{mw((fo / np.maximum(B, _EPS))[ok], w):>9.3f}"
          f"{mw((self_fg / np.maximum(B, _EPS))[ok], w):>8.3f}"
          f"{mw(uni['c_tau'][ok], w):>8.2f}{mw(uni['cm_g'][ok], w):>8.2f}")
