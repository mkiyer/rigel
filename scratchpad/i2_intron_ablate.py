"""INTRON STUDY step 2 — WHICH message channel drags introns off their factory self-solve?

Step 1: messages make introns worse in 16/16 strata (Δmsg +0.0002 .. +0.0275). The factory self-solve is good
(oracle 1.0000, self 0.9821) and the solve degrades it. Introns never RECEIVE a spliced graft (the graft is
gated on `ex_a`, exons only), so the suspects are the lambda (composition) message, the gDNA measurement
message, and the tilt message.

Also re-measures the factory's calibration WITHOUT the vertex artifact -- step 1's z2 of 20-26 came from nodes
whose oracle f_g is exactly 1.0 or 0.0, where logit is clipped. Restricted to 0.02 <= f_g <= 0.98 the honest
number is much smaller, and that matters: it decides whether the NB precision is over-confident or simply
weak.

    OMP_NUM_THREADS=1 python scratchpad/i2_intron_ablate.py
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
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "gdna_gdna5_ss_0.50_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
cc = cfg.calibration
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

print("INTRONS only. Error in READS. 'nomsg' = the pure factory self-solve.")
print(f"{'condition':<44}{'reads':>10}{'shipped':>9}{'-lam':>9}{'-gdna':>9}{'-rna':>9}"
      f"{'-tilt':>9}{'nomsg':>9}|{'z2 fac*':>9}{'n*':>5}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    res = calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cc, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    uni = cap["_uni"][-1]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    amb = fp & fn
    mass = np.asarray(cap["mass_global"])
    ship = np.asarray(cap["f_g"])
    solvable = np.asarray(cap["solvable"], bool)
    cp, cn = uni["cp"], uni["cn"]
    cR = cp + cn
    th_msg = np.arcsin(np.clip(np.where(cR > _EPS, (cp - cn) / np.maximum(cR, _EPS), 0.0), -1, 1))
    th_prec = np.where(amb, uni["cm_p"] + uni["cm_n"], 0.0)
    Z = np.zeros_like(uni["cm_g"])

    def psi(*, lam=True, gdna=True, rna=True, tilt=True):
        return _solve_nodes_logodds_all(
            st.u_pos, st.u_neg, fp, fn, st.mass_unspliced, st.mass_spliced,
            kappa=float(res.rna_sense_frac), od_g=float(res.gdna_strand_overdispersion),
            od_r=float(res.rna_strand_overdispersion),
            n_grid=int(cc.sweep_n_grid), L=float(cc.sweep_logodds_window),
            n_tilt=cc.sweep_n_tilt, n_grid_ss=cc.sweep_n_grid_single_strand,
            global_logprior=cap["global_lp"],
            gdna_imp_mode=uni["mo_g"], gdna_imp_prec=uni["cm_g"] if gdna else Z,
            rna_imp_mode=(uni["mo_p"], uni["mo_n"]),
            rna_imp_prec=(uni["cm_p"], uni["cm_n"]) if rna else (Z, Z),
            lam_logprior=cap["intron_prior"],
            lam_imp_mode=uni["lam_msg"], lam_imp_prec=uni["c_tau"] if lam else Z,
            theta_imp_mode=th_msg, theta_imp_prec=th_prec if tilt else Z,
            fg_ref=cap["fg_init"], fpos_ref=cap["fpos_init"], fneg_ref=cap["fneg_init"],
        ).gdna_frac

    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    is_intron = (kind == REGION) & (rt == 1)
    ok = np.isfinite(fo) & (mass > _EPS) & solvable & is_intron

    def er(v):
        return (np.abs(np.where(solvable, v, ship) - fo) * mass)[ok].sum()

    row = [er(np.where(solvable, psi(), ship)), er(psi(lam=False)), er(psi(gdna=False)),
           er(psi(rna=False)), er(psi(tilt=False)),
           er(psi(lam=False, gdna=False, rna=False, tilt=False))]
    # factory calibration, vertex-free
    tau = np.asarray(cap["_tau0_lam"])
    self_fg = np.asarray(cap["fg_loc"])
    interior = ok & (fo > 0.02) & (fo < 0.98) & (tau > _EPS)
    if interior.sum() > 5:
        w = mass[interior]
        lg = lambda x: np.log(np.clip(x, 1e-6, 1 - 1e-6) / (1 - np.clip(x, 1e-6, 1 - 1e-6)))  # noqa: E731
        z2 = (np.average((lg(self_fg[interior]) - lg(fo[interior])) ** 2, weights=w)
              / np.average(1.0 / tau[interior], weights=w))
    else:
        z2 = np.nan
    print(f"{cond[5:]:<44}{mass[ok].sum():>10,.0f}" + "".join(f"{x:>9,.0f}" for x in row)
          + f"|{z2:>9.2f}{int(interior.sum()):>5}")
