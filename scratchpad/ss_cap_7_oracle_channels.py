"""SINGLE-STRAND x CAPTURE, step 7 — what actually moves under capture? ORACLE channel densities only.

Step 6 killed both candidate frame factors (c=1 is exact off-capture; c=1/r is wrong off AND on; log r does
not predict log c on-capture). c is the ratio the graft must supply:

    c = [R_true(exon)/E_r(exon)] / [S(bnd)/E_spl(bnd)]

Both E's are capture-independent, so c's 4.5x move must be in R_true/S -- the exon's UNSPLICED oracle RNA
mass against the flanking junction's SPLICED mass. This measures the oracle channel densities directly, with
no solver in the loop, off vs on capture.

    OMP_NUM_THREADS=1 python scratchpad/ss_cap_7_oracle_channels.py
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
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def q(x, w, p):
    o = np.argsort(x)
    return x[o][np.searchsorted(np.cumsum(w[o]) / w[o].sum(), p)]


print("Per boundary->exon edge, ORACLE only (medians, mass-weighted). rho = mass/eff_len.")
print(f"{'condition':<44}{'n':>5}{'rhoR_ex':>10}{'rhoSpl_b':>10}{'log c':>8}|"
      f"{'rhoG_ex':>10}{'rhoG_b':>10}{'log capstep':>12}|{'R_ex/S_b':>10}{'Er/Espl':>9}")
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
    E_g, E_r = us["E_g"], us["E_r"]
    li, ri = us["left"], us["right"]
    n = len(E_g)
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    SPm = (us["SP_l"] + us["SN_l"], us["SP_r"] + us["SN_r"])
    rt, _ = _node_region_type(chain, ra)
    isex = (np.asarray(chain.kind) == REGION) & (rt == 2)
    mass = np.asarray(cap["mass_global"])

    A = {k: [] for k in ("rR", "rS", "lc", "rGe", "rGb", "cs", "RS", "EE", "w")}
    for src_i, sf in ((li, 1), (ri, 0)):
        s = np.clip(src_i, 0, n - 1)
        rho_R = R / np.maximum(E_r, _EPS)                      # exon unspliced RNA density
        rho_S = SPm[sf][s] / np.maximum(ESP[sf][s], _EPS)      # junction spliced density (source face)
        rho_Ge = G / np.maximum(E_g, _EPS)                     # exon gDNA density
        rho_Gb = G[s] / np.maximum(E_g[s], _EPS)               # boundary gDNA density (= capture step ref)
        m = (isex & (src_i >= 0) & us["is_bnd"][s] & (mass > _EPS) & (R > _EPS)
             & (SPm[sf][s] > _EPS) & (G > _EPS) & (G[s] > _EPS) & (ESP[sf][s] > _EPS))
        A["rR"].append(rho_R[m]); A["rS"].append(rho_S[m])
        A["lc"].append(np.log(rho_R[m] / rho_S[m]))
        A["rGe"].append(rho_Ge[m]); A["rGb"].append(rho_Gb[m])
        A["cs"].append(np.log(rho_Ge[m] / rho_Gb[m]))
        A["RS"].append(R[m] / SPm[sf][s][m])
        A["EE"].append(E_r[m] / ESP[sf][s][m])
        A["w"].append(mass[m])
    a = {k: np.concatenate(v) for k, v in A.items()}
    w = a["w"]
    print(f"{cond[5:]:<44}{w.size:>5}{q(a['rR'], w, .5):>10.4f}{q(a['rS'], w, .5):>10.4f}"
          f"{q(a['lc'], w, .5):>8.3f}|{q(a['rGe'], w, .5):>10.4f}{q(a['rGb'], w, .5):>10.4f}"
          f"{q(a['cs'], w, .5):>12.3f}|{q(a['RS'], w, .5):>10.4f}{q(a['EE'], w, .5):>9.3f}")
