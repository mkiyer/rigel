"""AMBIG study, step 5 — the CEILING: what is solving AMBIG actually worth, and what would deliver it?

Steps 1-4 showed: (i) the simplex bound B = 1-|d_hat| is a hard structural fact and, under capture, an almost
unbiased f_g estimate (err 0.010-0.021) while the solver sits 0.09-0.14 BELOW it; (ii) B - f_g = 2*min(f_+,f_-)
is exactly TWICE THE MINORITY STRAND'S RNA -- so "the truth sits at B" is a statement about the minority
strand, i.e. a population fact, not a prior-free deduction; (iii) damping the lambda message at stranded AMBIG
nodes (either by a DOF gate or by a bound-derived v_own) buys only ~0.0006 aggregate.

This bounds the remaining prize: suite mwae with the AMBIG nodes' answer replaced by (a) the bound B, and
(b) the ORACLE. The gap between shipped and (b) is what a correct AMBIG solve is worth; the gap between
(a) and (b) is what is left after the bound, i.e. what the hyperprior must actually supply.

    OMP_NUM_THREADS=1 python scratchpad/ambig_5_ceiling.py
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
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def mw(x, w):
    return np.average(x, weights=w) if np.sum(w) > 0 else np.nan


print("Suite mwae with the AMBIG answer replaced.  'AMBIG share' = AMBIG fraction of total err-mass.")
print(f"{'condition':<46}{'shipped':>9}{'AMB:=B':>9}{'AMB:=orc':>10}{'AMBshare':>10}"
      f"|{'ambShip':>9}{'ambB':>8}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    res = calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
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
    B = np.clip(1.0 - np.abs(d_hat), 0.0, 1.0)
    ship = np.asarray(cap["f_g"])
    mass = np.asarray(cap["mass_global"])
    ok = np.isfinite(fo) & (mass > _EPS)
    w, o = mass[ok], fo[ok]
    e_ship = np.abs(ship[ok] - o)
    a = amb[ok]
    e_B = np.where(a, np.abs(B[ok] - o), e_ship)
    e_or = np.where(a, 0.0, e_ship)
    print(f"{cond[5:]:<46}{mw(e_ship, w):>9.4f}{mw(e_B, w):>9.4f}{mw(e_or, w):>10.4f}"
          f"{(e_ship * w)[a].sum() / (e_ship * w).sum():>10.1%}"
          f"|{mw(e_ship[a], w[a]):>9.4f}{mw(np.abs(B[ok] - o)[a], w[a]):>8.4f}")
