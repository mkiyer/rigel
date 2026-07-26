"""AMBIG study, step 6 — validate the PRIOR-FREE design before building it.

Step 5: snapping AMBIG to the simplex bound B = 1-|d_hat| captures 85-91 % of the whole AMBIG prize on
stranded capture-ON, and is catastrophic on unstranded (d_hat is pure noise at kappa = 1/2).

But B is an UPPER bound, and using it as an estimate asserts the minority strand is silent -- a population
claim. The exact identity is

    f_g  =  1 - |d|  -  2*min(f_+, f_-)                    (d = f_+ - f_-, so 1-|d| = f_g + 2*min)

so the missing piece is the MINORITY STRAND'S RNA FRACTION -- and the solver already imputes that, per strand,
from that strand's own neighbours (`mo_p` / `mo_n`). That makes the constrained estimator PRIOR-FREE:

    f_g_hat  =  clip( 1 - |d_hat| - 2*min(f_+^imp, f_-^imp),  0,  1 )

This scores that estimator offline against the shipped answer, the bound, and the oracle -- with no solver
change -- so the design is validated (or refuted) before any implementation.

    OMP_NUM_THREADS=1 python scratchpad/ambig_6_design.py
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
    "gdna_gdna100_ss_0.99_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "gdna_none_ss_0.99_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def mw(x, w):
    return np.average(x, weights=w) if np.sum(w) > 0 else np.nan


print("AMBIG answer replaced offline.  'cons' = clip(1-|d| - 2*min(f+imp, f-imp)).  suite mwae / AMBIG mwae.")
print(f"{'condition':<44}{'ship':>8}{'B':>8}{'cons':>8}{'oracle':>8} |"
      f"{'ambShip':>9}{'ambB':>8}{'ambCons':>9}{'minorIm':>9}{'minorOr':>9}")
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
    minor_or = np.where(T > _EPS, np.minimum(Rp, Rn) / np.maximum(T, _EPS), 0.0)
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    amb = fp & fn
    up, un = np.asarray(st.u_pos, float), np.asarray(st.u_neg, float)
    n = up + un
    kap = float(res.rna_sense_frac)
    p_obs = np.where(n > 0, up / np.maximum(n, _EPS), 0.5)
    den = kap - 0.5
    d_hat = (p_obs - 0.5) / den if abs(den) > _EPS else np.zeros_like(p_obs)
    B = np.clip(1.0 - np.abs(d_hat), 0.0, 1.0)
    # the per-strand imputed RNA fractions the combine already delivers to psi
    f_imp_p = np.clip(np.exp(uni["mo_p"]), 0.0, 1.0)
    f_imp_n = np.clip(np.exp(uni["mo_n"]), 0.0, 1.0)
    minor_im = np.minimum(f_imp_p, f_imp_n)
    cons = np.clip(B - 2.0 * minor_im, 0.0, 1.0)
    ship = np.asarray(cap["f_g"])
    mass = np.asarray(cap["mass_global"])
    ok = np.isfinite(fo) & (mass > _EPS)
    w, o, a = mass[ok], fo[ok], amb[ok]
    e_ship = np.abs(ship[ok] - o)
    def suite(v):
        return mw(np.where(a, np.abs(v[ok] - o), e_ship), w)
    print(f"{cond[5:]:<44}{mw(e_ship, w):>8.4f}{suite(B):>8.4f}{suite(cons):>8.4f}"
          f"{mw(np.where(a, 0.0, e_ship), w):>8.4f} |{mw(e_ship[a], w[a]):>9.4f}"
          f"{mw(np.abs(B[ok] - o)[a], w[a]):>8.4f}{mw(np.abs(cons[ok] - o)[a], w[a]):>9.4f}"
          f"{mw(minor_im[ok][a], w[a]):>9.4f}{mw(minor_or[ok][a], w[a]):>9.4f}")
