"""P1e step 1 — MEASURE the conservation surprise on the solver's real PRE-PIN messages.

`_capture["_pin"]` publishes exactly the state `_pin_v` is about to erase: the per-message per-component
densities (tg, tp, tn), their mode-fusion precisions (tpg, tpp, tpn), the common-mode variance legs
(s2t = M5's Var(log r), n_src), and the graft mask. This script forms, per message:

    S = tg*E_g + (tp+tn)*E_r        the fragment count the claim ASSERTS
    delta = log(M/S)                the conservation residual `_pin_v` discards
    alpha_c = m_c/S                 the claimed mass shares (sum 1)
    Sigma = s_cm^2 11^T + diag(w)   w_c = max(0, v_c - s_cm^2),  s_cm^2 = s2t + 1/n_src
    sd2 = alpha^T Sigma alpha = s_cm^2 + sum_c alpha_c^2 w_c
    z2  = delta^2 / sd2

and reports, per EDGE CLASS (plain / graft / peel) and destination class:
  * the size of delta against the EFF-LENGTH BOUND |log(E_g/E_r)| -- the ceiling a composition error alone
    can produce on a matched reframe (the ROADMAP's x1.04 / x1.50);
  * z2, with and without the +1/n_dst the sketch proposes;
  * the DIRECTION split s = Sigma alpha: how much of it is COMMON (s_cm^2) vs DIFFERENTIAL (alpha_c w_c),
    which is what decides whether P1e's rank-1 inflation reaches the lambda channel at all;
  * whether delta predicts the delivered message's composition error against the ORACLE.

    OMP_NUM_THREADS=1 python scratchpad/p1e_1_delta.py [cond ...]
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
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "gdna_none_ss_0.50_nrna_present_capture_off",
    "gdna_gdna1_ss_0.50_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def solve(cond):
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    return inp, dbg


ROWS = []
for cond in CONDS:
    inp, dbg = solve(cond)
    chain, cap = dbg["chain"], dbg["capture"]
    us = cap["_uni_static"]
    uni_last = cap["_uni"][-1]
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    n = M.shape[0]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    li, ri = us["left"], us["right"]
    og, op, on = us["og"], us["op"], us["on"]
    n_node = np.where(~is_bnd, us["n_unspl_l"], us["n_unspl_l"] + us["n_unspl_r"])
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    CLS = {0: "intergenic", 1: "intron", 2: "exon"}
    dcls = np.array(
        [CLS.get(int(rt[j]), "?") if kind[j] == REGION else "boundary" for j in range(n)]
    )

    # ── the exact M-cancellation check: for a REGION destination the face density is M*B(f_g), so S/M is
    # M-free.  B is recomputed from f_cur (the belief the last rho-iteration reframed on) with no mass.
    f_cur = np.clip(np.asarray(cap["_uni"][-2]["fg_out"], float), 0.0, 1.0) if len(cap["_uni"]) > 1 \
        else None
    if f_cur is not None:
        B = f_cur / np.maximum(E_g, _EPS) + (1.0 - f_cur) / np.maximum(E_r, _EPS)
        reg = ~is_bnd & (M > _EPS)
        rel = np.abs(uni_last["rho_lf"][reg] / np.maximum(M[reg], _EPS) - B[reg]) / np.maximum(B[reg], _EPS)
        print(f"[{cond[5:]}]  M-cancellation: max |rho_face/M - B(f_g)| / B over REGION nodes "
              f"= {float(np.nanmax(rel)):.3e}   (0 => S is proportional to M => delta is M-free)")

    for ent in cap["_pin"][-2:]:  # the LAST rho-iteration's two messages
        df, src, valid = ent["df"], ent["src"], ent["valid"]
        tg, tp, tn = ent["tg"], ent["tp"], ent["tn"]
        tpg, tpp, tpn = ent["tpg"], ent["tpp"], ent["tpn"]
        s2t, nsrc, graft = ent["s2t"], ent["n_src"], ent["graft"]
        peel = is_bnd & is_exon[src] & valid
        plain = valid & ~graft & ~peel

        p3 = np.stack([tpg, tpp, tpn], axis=1)
        d3 = np.stack([tg, tp, tn], axis=1)
        o3 = np.stack([og, op, on], axis=1)
        E3 = np.stack([E_g, E_r, E_r], axis=1)
        sup = p3 > 0.0
        m3 = np.where(sup, d3, o3) * E3          # `_pin_v`'s partial-claim mass budget
        S = m3.sum(axis=1)
        ok = valid & (S > _EPS) & (M > _EPS)
        alpha = m3 / np.maximum(S, _EPS)[:, None]
        v3 = np.where(sup, 1.0 / np.maximum(p3, _EPS), np.inf)
        scm = np.maximum(s2t, 0.0) + 1.0 / np.maximum(nsrc, _EPS)
        w3 = np.maximum(np.where(sup, v3, 0.0) - scm[:, None], 0.0)
        s_c = np.where(sup, scm[:, None] + alpha * w3, 0.0)
        sd2 = scm + (alpha * alpha * w3).sum(axis=1)     # alpha^T Sigma alpha (sum alpha = 1)
        delta = np.log(np.maximum(M, _EPS) / np.maximum(S, _EPS))
        bound = np.abs(np.log(np.maximum(E_g, _EPS) / np.maximum(E_r, _EPS)))
        z2 = delta * delta / np.maximum(sd2, _EPS)
        z2n = delta * delta / np.maximum(sd2 + 1.0 / np.maximum(n_node, _EPS), _EPS)
        # the delivered composition and its oracle error
        num = np.where(sup[:, 0], tg, og) * E_g
        den = num + (np.where(sup[:, 1], tp, op) + np.where(sup[:, 2], tn, on)) * E_r
        fmsg = np.where(den > _EPS, num / np.maximum(den, _EPS), np.nan)
        # the DIRECTION: common vs differential part of s
        com = scm
        dif = np.abs(s_c[:, 0] - np.where(alpha[:, 1] + alpha[:, 2] > 0,
                                          (alpha[:, 1] * s_c[:, 1] + alpha[:, 2] * s_c[:, 2])
                                          / np.maximum(alpha[:, 1] + alpha[:, 2], _EPS), 0.0))
        for name, msk in (("plain", plain), ("graft", graft), ("peel", peel)):
            m = msk & ok & (M > _EPS)
            if m.sum() < 5:
                continue
            ROWS.append(dict(
                cond=cond[5:], edge=name, n=int(m.sum()), mass=float(M[m].sum()),
                med_absd=float(np.median(np.abs(delta[m]))),
                p90_absd=float(np.percentile(np.abs(delta[m]), 90)),
                med_ratio=float(np.exp(np.median(-delta[m]))),
                bound=float(np.median(bound[m])),
                over=float(np.mean(np.abs(delta[m]) > bound[m])),
                z2=float(np.mean(z2[m])), z2med=float(np.median(z2[m])),
                z2n=float(np.mean(z2n[m])),
                fr_com=float(np.median(com[m] / np.maximum(sd2[m], _EPS))),
                dif_rel=float(np.median(dif[m] / np.maximum(s_c[m, 0], _EPS))),
                err=float(np.nanmean(np.abs(fmsg[m] - fo[m]))),
            ))
        # delta-vs-error, exon destinations only, this message
        m = ok & (dcls == "exon") & np.isfinite(fo) & np.isfinite(fmsg)
        if m.sum() > 50:
            q = np.quantile(np.abs(delta[m]), [0.25, 0.5, 0.75])
            bins = np.digitize(np.abs(delta[m]), q)
            e = [float(np.mean(np.abs(fmsg[m][bins == b] - fo[m][bins == b]))) for b in range(4)]
            zz = [float(np.median(z2[m][bins == b])) for b in range(4)]
            print(f"    df={df} EXON dst  |delta| quartiles {np.round(q, 3)}  "
                  f"-> mean |f_g^msg - oracle| {np.round(e, 4)}   z2 med {np.round(zz, 2)}")

print()
hdr = ("cond", "edge", "n", "med|d|", "p90|d|", "S/M", "bound", ">bnd", "z2", "z2med",
       "z2+1/n", "com/sd2", "dif/s_g", "|f-orc|")
print(f"{hdr[0]:<42}{hdr[1]:>7}{hdr[2]:>7}" + "".join(f"{h:>9}" for h in hdr[3:]))
print("-" * (42 + 14 + 9 * 11))
for r in ROWS:
    print(f"{r['cond']:<42}{r['edge']:>7}{r['n']:>7}{r['med_absd']:>9.3f}{r['p90_absd']:>9.3f}"
          f"{1.0 / r['med_ratio']:>9.3f}{r['bound']:>9.3f}{r['over']:>9.1%}"
          f"{r['z2']:>9.1f}{r['z2med']:>9.2f}{r['z2n']:>9.1f}{r['fr_com']:>9.2f}"
          f"{r['dif_rel']:>9.2f}{r['err']:>9.3f}")
