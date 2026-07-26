"""SINGLE-STRAND x CAPTURE, step 3 — an EXACT psi replay + per-channel ABLATION.

The suite table says the exon degradation under capture is a MESSAGE defect (self flat, solved 4-19x worse).
This replays the shipped psi solve bit-for-bit from the captured channels, verifies fidelity, then ablates
each channel to attribute the damage:

    lam   : the single-lambda COMPOSITION message (lam_msg, c_tau)
    gdna  : the gDNA MEASUREMENT message          (mo_g,  cm_g)
    rna   : the spliced RNA MEASUREMENT messages  (mo_p/mo_n, cm_p/cm_n)

    OMP_NUM_THREADS=1 python scratchpad/ss_cap_3_replay.py [COND ...]
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
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.99_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
cc = cfg.calibration
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def run(cond):
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
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    fp = np.asarray(st.free_pos, bool)
    fn = np.asarray(st.free_neg, bool)
    is_amb = fp & fn
    cp, cn = uni["cp"], uni["cn"]
    cR = cp + cn
    tau_tilt = np.clip(np.where(cR > _EPS, (cp - cn) / np.maximum(cR, _EPS), 0.0), -1.0, 1.0)
    th_msg = np.arcsin(tau_tilt)
    th_prec = np.where(is_amb, uni["cm_p"] + uni["cm_n"], 0.0)

    def psi(*, lam=True, gdna=True, rna=True):
        return _solve_nodes_logodds_all(
            st.u_pos, st.u_neg, fp, fn, st.mass_unspliced, st.mass_spliced,
            kappa=float(res.rna_sense_frac),
            od_g=float(res.gdna_strand_overdispersion),
            od_r=float(res.rna_strand_overdispersion),
            n_grid=int(cc.sweep_n_grid), L=float(cc.sweep_logodds_window),
            n_tilt=cc.sweep_n_tilt, n_grid_ss=cc.sweep_n_grid_single_strand,
            global_logprior=cap["global_lp"],
            gdna_imp_mode=uni["mo_g"], gdna_imp_prec=uni["cm_g"] if gdna else np.zeros_like(uni["cm_g"]),
            rna_imp_mode=(uni["mo_p"], uni["mo_n"]),
            rna_imp_prec=(uni["cm_p"], uni["cm_n"]) if rna
            else (np.zeros_like(uni["cm_p"]), np.zeros_like(uni["cm_n"])),
            lam_logprior=cap["intron_prior"],
            lam_imp_mode=uni["lam_msg"],
            lam_imp_prec=uni["c_tau"] if lam else np.zeros_like(uni["c_tau"]),
            theta_imp_mode=th_msg, theta_imp_prec=th_prec,
            fg_ref=cap["fg_init"], fpos_ref=cap["fpos_init"], fneg_ref=cap["fneg_init"],
        ).gdna_frac

    solvable = np.asarray(cap["solvable"], bool)
    shipped = np.asarray(cap["f_g"])
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
    cls = np.array([CLS[int(rt[i])] if kind[i] == REGION else "boundary" for i in range(len(shipped))])
    mass = np.asarray(cap["mass_global"])
    ok = np.isfinite(fo) & (mass > _EPS) & ~is_amb
    return dict(fo=fo, mass=mass, ok=ok, cls=cls, solvable=solvable, shipped=shipped,
                self_fg=np.asarray(cap["fg_loc"]), psi=psi, cap=cap, uni=uni)


def mw(x, w):
    return np.average(x, weights=w) if np.sum(w) > 0 else np.nan


print(f"{'condition':<44}{'stratum':<10}{'REPLAY':>9}{'shipped':>9}{'self':>9}"
      f"{'-lam':>9}{'-gdna':>9}{'-rna':>9}{'nomsg':>9}")
for cond in CONDS:
    d = run(cond)
    base = np.where(d["solvable"], d["psi"](), d["shipped"])
    fid = np.max(np.abs(base - d["shipped"]))
    arms = {
        "-lam": d["psi"](lam=False), "-gdna": d["psi"](gdna=False), "-rna": d["psi"](rna=False),
        "nomsg": d["psi"](lam=False, gdna=False, rna=False),
    }
    for strat in ("ALL", "exon", "boundary", "intron"):
        m = d["ok"] & d["solvable"] & (np.ones_like(d["ok"]) if strat == "ALL" else (d["cls"] == strat))
        if not m.any():
            continue
        w, o = d["mass"][m], d["fo"][m]
        row = [mw(np.abs(base[m] - o), w), mw(np.abs(d["shipped"][m] - o), w),
               mw(np.abs(d["self_fg"][m] - o), w)]
        row += [mw(np.abs(np.where(d["solvable"], v, d["shipped"])[m] - o), w) for v in arms.values()]
        lab = cond[5:] if strat == "ALL" else ""
        print(f"{lab:<44}{strat:<10}" + "".join(f"{x:>9.4f}" for x in row))
    print(f"{'':<44}{'(replay fidelity max|Δ| = ' + f'{fid:.2e})':<10}")
