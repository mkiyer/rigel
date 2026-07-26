"""TARGET step 6 — the CEILING: how much of the residual is addressable by a better MODE?

Step 4 showed the gDNA channel's error is the transport ratio, and that `r = rho_tot(dst)/rho_tot(src)` is
structurally not the capture step whenever composition changes (capture scales both components alike, so
rho_tot ratio = capture x content-ratio). Step 5 refuted the left/right gap as a mismatch statistic.

So: how much of the remaining error could ANY prior-free improvement to the message MODES buy? This uses the
bit-exact psi replay and substitutes ORACLE modes into one channel at a time, keeping every precision as
shipped. It is an upper bound on mode work, and the gap it leaves is what only the hyperprior can supply.

    OMP_NUM_THREADS=1 python scratchpad/t6_ceiling.py
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
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
cc = cfg.calibration
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

print(f"{'condition':<44}{'shipped':>11}{'orc mo_g':>11}{'orc mo_R':>11}{'orc BOTH':>11}"
      f"{'x10 prec':>11}{'nomsg':>11}")
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
    fop = np.where(T > _EPS, Rp / np.maximum(T, _EPS), 0.0)
    fon = np.where(T > _EPS, Rn / np.maximum(T, _EPS), 0.0)
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
    L = np.log

    def psi(mo_g=None, mo_p=None, mo_n=None, lam=None, pscale=1.0, off=False):
        return _solve_nodes_logodds_all(
            st.u_pos, st.u_neg, fp, fn, st.mass_unspliced, st.mass_spliced,
            kappa=float(res.rna_sense_frac), od_g=float(res.gdna_strand_overdispersion),
            od_r=float(res.rna_strand_overdispersion),
            n_grid=int(cc.sweep_n_grid), L=float(cc.sweep_logodds_window),
            n_tilt=cc.sweep_n_tilt, n_grid_ss=cc.sweep_n_grid_single_strand,
            global_logprior=cap["global_lp"],
            gdna_imp_mode=uni["mo_g"] if mo_g is None else mo_g,
            gdna_imp_prec=Z if off else uni["cm_g"] * pscale,
            rna_imp_mode=(uni["mo_p"] if mo_p is None else mo_p, uni["mo_n"] if mo_n is None else mo_n),
            rna_imp_prec=(Z, Z) if off else (uni["cm_p"] * pscale, uni["cm_n"] * pscale),
            lam_logprior=cap["intron_prior"],
            lam_imp_mode=uni["lam_msg"] if lam is None else lam,
            lam_imp_prec=Z if off else uni["c_tau"] * pscale,
            theta_imp_mode=th_msg, theta_imp_prec=Z if off else th_prec * pscale,
            fg_ref=cap["fg_init"], fpos_ref=cap["fpos_init"], fneg_ref=cap["fneg_init"],
        ).gdna_frac

    ok = np.isfinite(fo) & (mass > _EPS) & solvable

    def errof(v):
        return (np.abs(np.where(solvable, v, ship) - fo) * mass)[ok].sum()

    og = L(np.maximum(fo, 1e-9))
    op, on = L(np.maximum(fop, 1e-9)), L(np.maximum(fon, 1e-9))
    olam = og - L(np.maximum(fop + fon, 1e-9))
    row = [errof(np.where(solvable, psi(), ship)), errof(psi(mo_g=og, lam=olam)),
           errof(psi(mo_p=op, mo_n=on)), errof(psi(mo_g=og, mo_p=op, mo_n=on, lam=olam)),
           errof(psi(pscale=10.0)), errof(psi(off=True))]
    print(f"{cond[5:]:<44}" + "".join(f"{x:>11,.0f}" for x in row))
