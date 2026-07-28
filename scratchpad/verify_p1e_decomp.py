"""ADVERSARIAL check of P1e's §0 SEMANTICS and of the null's practical size.

(1) THE EXACT DECOMPOSITION.  With Lambda_src = sum_c rho_c^relay E_c^src / M_src (the source's own
    conservation, measured ~1), F = M*B + sigma_spl and beta_c = rho_c^relay E_c^src / (Lambda_src M_src)
    the source's MASS shares, k_c = E_c^dst/E_c^src:

        delta = log(M_dst/S)
              = log(B_src/B_dst) + log((1+q_src)/(1+q_dst)) - log(sum_c beta_c k_c) - log(Lambda_src)

    The report reads delta as "the message's composition disagreeing with the destination's belief".
    But the MESSAGE's composition beta enters ONLY through log(sum_c beta_c k_c) -- and that term is
    IDENTICALLY ZERO whenever k_g = k_R.  The first term depends on NO message quantity at all: it is the
    two NODES' current BELIEF compositions.  This measures which term carries delta.

(2) HOW BIG IS THE NULL ERROR THE REPORT CALLS "THE EXPENSIVE ONE"?  Compare
        b2_correct = max(0, d^2 - sd2 - (q/(1+q))^2/n_dst)   vs   b2_sketch = max(0, d^2 - sd2 - 1/n_dst)
    on the real messages: fraction that change firing state, and the median relative change in b2.

    OMP_NUM_THREADS=1 python scratchpad/verify_p1e_decomp.py
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

print(f"{'cond':<34}{'edge':>7}{'n':>6}{'sd(d)':>8}{'sd(Brat)':>10}{'sd(bk)':>8}"
      f"{'sd(qrat)':>10}{'sd(Lam)':>9}{'corr(d,Brat)':>14}{'corr(d,-bk)':>13}")
print("-" * 121)
ROWS2 = []
for cond in CONDS:
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    cap = dbg["capture"]
    us, uni = cap["_uni_static"], cap["_uni"][-1]
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    og, op, on = us["og"], us["op"], us["on"]
    n_node = np.where(~is_bnd, us["n_unspl_l"], us["n_unspl_l"] + us["n_unspl_r"])
    f_cur = np.clip(np.asarray(cap["_uni"][-2]["fg_out"], float), 0.0, 1.0)
    B = f_cur / np.maximum(E_g, _EPS) + (1.0 - f_cur) / np.maximum(E_r, _EPS)
    # the relayed claim per node, as it entered the combine (rho_l0 / rho_r0 are the relay outputs)
    for ent in cap["_pin"][-2:]:
        df, src, valid = ent["df"], ent["src"], ent["valid"]
        face_d = uni["rho_lf"] if df == 0 else uni["rho_rf"]
        face_s = (uni["rho_rf"] if df == 0 else uni["rho_lf"])[src]   # the SOURCE's outgoing face
        q_d = np.maximum(face_d / np.maximum(M * B, _EPS) - 1.0, 0.0)
        q_s = np.maximum(face_s / np.maximum((M * B)[src], _EPS) - 1.0, 0.0)
        tg, tp, tn = ent["tg"], ent["tp"], ent["tn"]
        p3 = np.stack([ent["tpg"], ent["tpp"], ent["tpn"]], axis=1)
        sup = p3 > 0.0
        d3 = np.stack([tg, tp, tn], axis=1)
        o3 = np.stack([og, op, on], axis=1)
        E3d = np.stack([E_g, E_r, E_r], axis=1)
        E3s = np.stack([E_g[src], E_r[src], E_r[src]], axis=1)
        m3 = np.where(sup, d3, o3) * E3d
        S = m3.sum(axis=1)
        ok = valid & (S > _EPS) & (M > _EPS) & (M[src] > _EPS)
        delta = np.log(np.maximum(M, _EPS) / np.maximum(S, _EPS))
        # the message's own MASS shares expressed in the SOURCE frame (undo r, which is a common factor)
        ms = np.where(sup, d3, o3) * E3s
        Ssrc = ms.sum(axis=1)
        beta = ms / np.maximum(Ssrc, _EPS)[:, None]
        lam_src = Ssrc / np.maximum(M[src], _EPS)
        k = E3d / np.maximum(E3s, _EPS)
        bk = np.log(np.maximum((beta * k).sum(axis=1), _EPS))
        Brat = np.log(np.maximum(B[src], _EPS) / np.maximum(B, _EPS))
        qrat = np.log((1.0 + q_s) / (1.0 + q_d))
        resid = delta - (Brat + qrat - bk - np.log(np.maximum(lam_src, _EPS)))
        graft = ent["graft"]
        peel = is_bnd & is_exon[src] & valid
        plain = valid & ~graft & ~peel
        for name, msk in (("plain", plain), ("graft", graft), ("peel", peel)):
            m = msk & ok & np.isfinite(delta) & np.isfinite(Brat) & np.isfinite(bk)
            if m.sum() < 20:
                continue

            def _c(a, b, m=m):
                a, b = a[m], b[m]
                if a.std() < 1e-12 or b.std() < 1e-12:
                    return float("nan")
                return float(np.corrcoef(a, b)[0, 1])
            print(f"{cond[5:]:<34}{name:>7}{int(m.sum()):>6}{delta[m].std():>8.3f}"
                  f"{Brat[m].std():>10.3f}{bk[m].std():>8.3f}{qrat[m].std():>10.3f}"
                  f"{np.log(np.maximum(lam_src[m], _EPS)).std():>9.3f}"
                  f"{_c(delta, Brat):>14.3f}{_c(delta, -bk):>13.3f}")
            if name == "plain":
                print(f"{'':<34}{'':>7}      max|residual of the identity| = "
                      f"{float(np.nanmax(np.abs(resid[m]))):.3e}")
        # (2) the null comparison
        scm = np.maximum(ent["s2t"], 0.0) + 1.0 / np.maximum(ent["n_src"], _EPS)
        v3 = np.where(sup, 1.0 / np.maximum(p3, _EPS), np.inf)
        alpha = m3 / np.maximum(S, _EPS)[:, None]
        w3 = np.maximum(np.where(sup, v3, 0.0) - scm[:, None], 0.0)
        sd2 = scm + (alpha * alpha * w3).sum(axis=1)
        nu = (q_d / (1.0 + q_d)) ** 2 / np.maximum(n_node, _EPS)
        b_ok = np.maximum(0.0, delta**2 - sd2 - nu)
        b_sk = np.maximum(0.0, delta**2 - sd2 - 1.0 / np.maximum(n_node, _EPS))
        mm = ok
        f_ok, f_sk = b_ok[mm] > 0, b_sk[mm] > 0
        fire_change = float(np.mean(f_ok != f_sk))
        both = f_ok & f_sk
        rel = np.median((b_ok[mm][both] - b_sk[mm][both]) / np.maximum(b_ok[mm][both], _EPS)) \
            if both.any() else float("nan")
        ROWS2.append((cond[5:], df, int(mm.sum()), float(np.median(1.0 / np.maximum(n_node[mm], _EPS))),
                      float(np.median(sd2[mm])), float(np.median(b_ok[mm][f_ok])) if f_ok.any() else 0.0,
                      fire_change, float(rel)))

print("\n(2) the '+1/n_dst' null error, sized on the real messages")
print(f"{'cond':<34}{'df':>3}{'n':>7}{'med 1/n_dst':>13}{'med sd2':>10}{'med b2(fire)':>14}"
      f"{'fire flips':>12}{'med rel db2':>13}")
for r in ROWS2:
    print(f"{r[0]:<34}{r[1]:>3}{r[2]:>7}{r[3]:>13.4f}{r[4]:>10.4f}{r[5]:>14.4f}"
          f"{r[6]:>11.2%}{r[7]:>13.2%}")
