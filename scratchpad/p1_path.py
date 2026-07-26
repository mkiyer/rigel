"""PATH STUDY step 1 — the exon -> (exon|intron boundary) -> intron message path, hop by hop.

OWNER: the messages themselves are bad. They travel through exon/intron boundaries, are contaminated with the
exon's RNA, and over-estimate RNA -- which then poisons introns, which were accurate. Gating messages AT the
intron is a band-aid on the symptom.

The claim is sharply testable in `nrna_none`: at an exon->intron seam MATURE RNA SPLICES OUT and there is no
nascent RNA, so ~ZERO RNA continues past the seam. The oracle says the boundary's unspliced RNA is ~0 and the
intron's is exactly 0. The PEEL (`tp - spl_*_f[df]`) exists precisely to strip the departing mature. So:

    if the relay's RNA density arriving at the boundary (post-peel) is >> the oracle's, the peel is
    UNDER-CONSUMING and the exon's RNA is leaking across a seam it should not cross.

This traces the density at each hop against the oracle: at the EXON, at the exon|intron BOUNDARY, and at the
INTRON.

    OMP_NUM_THREADS=1 python scratchpad/p1_path.py
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
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def q(x, w, p=0.5):
    o = np.argsort(x)
    return x[o][np.searchsorted(np.cumsum(w[o]) / max(w[o].sum(), _EPS), p)]


for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st, geo = dbg["chain"], dbg["capture"], dbg["statics"], dbg["geometry"]
    us = cap["_uni_static"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    E_g, E_r, M = us["E_g"], us["E_r"], us["M"]
    n = len(M)
    mass = np.asarray(cap["mass_global"])
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array([CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(n)])
    li, ri = us["left"], us["right"]
    rho_R_true = np.where(E_r > _EPS, R / np.maximum(E_r, _EPS), np.nan)
    rho_g_true = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)
    # the relay's running RNA density (forward from the left, backward from the right)
    fwd_R = us["fwd_p"] + us["fwd_n"]
    bwd_R = us["bwd_p"] + us["bwd_n"]
    fwd_g, bwd_g = us["fwd_g"], us["bwd_g"]

    # classify boundaries by their flanking pair
    is_bnd = kind != REGION
    pair = np.array(["-"] * n, dtype=object)
    for i in np.flatnonzero(is_bnd):
        a = cls[li[i]] if li[i] >= 0 else "edge"
        b = cls[ri[i]] if ri[i] >= 0 else "edge"
        pair[i] = " | ".join(sorted([str(a), str(b)]))
    ex_int = is_bnd & (pair == "exon | intron")

    print(f"\n{'=' * 112}\n{cond[5:]}\n{'=' * 112}")
    print(f"  exon|intron boundaries: {int(ex_int.sum())}   "
          f"oracle RNA FRACTION there (median, mass-wtd): "
          f"{q(np.where(G + R > _EPS, R / np.maximum(G + R, _EPS), 0.0)[ex_int], mass[ex_int]):.4f}")
    print(f"  introns:                {int((cls == 'intron').sum())}   "
          f"oracle RNA fraction: "
          f"{q(np.where(G + R > _EPS, R / np.maximum(G + R, _EPS), 0.0)[cls == 'intron'], mass[cls == 'intron']):.4f}")

    print(f"\n  {'node kind':<26}{'n':>5}{'oracle rho_R':>14}{'relayed rho_R':>15}"
          f"{'ratio':>10}{'oracle rho_g':>14}{'relayed rho_g':>15}{'ratio':>9}")
    for lab, m, rr, rg in (
        ("EXON (source)", (cls == "exon") & (mass > _EPS), fwd_R, fwd_g),
        ("exon|intron BOUNDARY", ex_int & (mass > _EPS), fwd_R, fwd_g),
        ("INTRON", (cls == "intron") & (mass > _EPS), fwd_R, fwd_g),
    ):
        mm = m & np.isfinite(rho_R_true) & np.isfinite(rho_g_true)
        if mm.sum() < 3:
            continue
        w = mass[mm]
        oR, cR_ = rho_R_true[mm], rr[mm]
        oG, cG = rho_g_true[mm], rg[mm]
        print(f"  {lab:<26}{int(mm.sum()):>5}{q(oR, w):>14.4f}{q(cR_, w):>15.4f}"
              f"{q(cR_ / np.maximum(oR, _EPS), w):>10.1f}{q(oG, w):>14.4f}{q(cG, w):>15.4f}"
              f"{q(cG / np.maximum(oG, _EPS), w):>9.2f}")

    # ── the PEEL: on exon -> boundary edges, how much RNA does the relay carry across vs what should? ──
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    SP, SN = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    print(f"\n  THE PEEL on exon -> exon|intron boundary edges "
          f"(what CONTINUES past the seam; oracle says ~0 in nrna_none)")
    print(f"  {'dir':<8}{'n':>5}{'rho_R(exon)':>13}{'reframed':>11}{'peel amt':>11}"
          f"{'post-peel':>11}{'oracle@bnd':>12}{'over-claim':>12}{'peel==0':>9}")
    for lab, src_i, sf, df, face_s, face_d in (
        ("L->R", li, 1, 0, us["rho_r0"], us["rho_l0"]),
        ("R->L", ri, 0, 1, us["rho_l0"], us["rho_r0"]),
    ):
        s = np.clip(src_i, 0, n - 1)
        m = ex_int & (src_i >= 0) & (cls[s] == "exon") & (mass > _EPS) & np.isfinite(rho_R_true)
        if m.sum() < 3:
            continue
        rho_src = (us["fwd_p"] + us["fwd_n"])[s] if lab == "L->R" else (us["bwd_p"] + us["bwd_n"])[s]
        r = np.where(face_s[s] > _EPS, face_d / np.maximum(face_s[s], _EPS), 1.0)
        reframed = rho_src * r
        peel = (np.where(SP[df] > _EPS, SP[df] / np.maximum(ESP[df], _EPS), 0.0)
                + np.where(SN[df] > _EPS, SN[df] / np.maximum(ESP[df], _EPS), 0.0))
        post = np.maximum(reframed - peel, 0.0)
        w = mass[m]
        print(f"  {lab:<8}{int(m.sum()):>5}{q(rho_src[m], w):>13.4f}{q(reframed[m], w):>11.4f}"
              f"{q(peel[m], w):>11.4f}{q(post[m], w):>11.4f}{q(rho_R_true[m], w):>12.4f}"
              f"{q((post / np.maximum(rho_R_true, _EPS))[m], w):>12.1f}"
              f"{np.average((peel[m] <= _EPS).astype(float), weights=w):>9.0%}")
