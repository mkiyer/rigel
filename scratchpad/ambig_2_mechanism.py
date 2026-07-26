"""AMBIG study, step 2 — the mechanism: does psi extract the simplex bound, and is DL inert?

Step 1: on stranded data the bound f_g <= 1-|d| is TIGHT on 96-97 % of AMBIG mass (slack 0.008-0.014), so the
strand tilt identifies f_g -- yet `tau_own = 0` on every AMBIG node (node_init gates the strand lambda-term to
single-strand: the Schur complement of an INTERIOR rank-1 Fisher). Messages then TRIPLE the error.

Two things to establish:
  (a) does the message-free self-solve track 1-|d_obs|?  (i.e. psi's 2-D cube already sees the bound)
  (b) is the damage the MESSAGES, and is the DL term inert because v_own = 1/tau_own = inf?

Uses a bit-exact psi replay (max|delta| = 0), now covering AMBIG nodes too.

    OMP_NUM_THREADS=1 python scratchpad/ambig_2_mechanism.py
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
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
cc = cfg.calibration
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def mw(x, w):
    return np.average(x, weights=w) if np.sum(w) > 0 else np.nan


print(f"{'condition':<44}{'stratum':<10}{'REPLAY':>9}{'ship':>8}{'self':>8}"
      f"{'-lam':>8}{'-gdna':>8}{'-rna':>8}{'-tilt':>8}{'nomsg':>8}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    res = calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cc, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    uni = cap["_uni"][-1]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    T = G + R
    fo = np.where(T > _EPS, G / np.maximum(T, _EPS), np.nan)
    d_abs = np.where(T > _EPS, np.abs(Rp - Rn) / np.maximum(T, _EPS), np.nan)
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    amb = fp & fn
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

    solvable = np.asarray(cap["solvable"], bool)
    ship = np.asarray(cap["f_g"])
    mass = np.asarray(cap["mass_global"])
    self_fg = np.asarray(cap["fg_loc"])
    base = np.where(solvable, psi(), ship)
    fid = np.max(np.abs(base - ship))
    arms = {"-lam": psi(lam=False), "-gdna": psi(gdna=False), "-rna": psi(rna=False),
            "-tilt": psi(tilt=False), "nomsg": psi(lam=False, gdna=False, rna=False, tilt=False)}
    ok = np.isfinite(fo) & (mass > _EPS) & solvable
    for lab, m in (("AMBIG", ok & amb), ("single", ok & ~amb)):
        w, o = mass[m], fo[m]
        row = [mw(np.abs(base[m] - o), w), mw(np.abs(ship[m] - o), w), mw(np.abs(self_fg[m] - o), w)]
        row += [mw(np.abs(np.where(solvable, v, ship)[m] - o), w) for v in arms.values()]
        print(f"{cond[5:] if lab == 'AMBIG' else '':<44}{lab:<10}" + "".join(f"{x:>8.4f}" for x in row))
    # (a) does the self-solve track the simplex bound 1-|d| ?
    a = ok & amb
    w = mass[a]
    print(f"{'':<44}{'':<10}  self vs oracle {mw(np.abs(self_fg - fo)[a], w):.4f} | "
          f"self vs (1-|d|) {mw(np.abs(self_fg - (1 - d_abs))[a], w):.4f} | "
          f"(1-|d|) vs oracle {mw(np.abs(1 - d_abs - fo)[a], w):.4f} | replay fid {fid:.1e}")
