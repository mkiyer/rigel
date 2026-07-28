"""ADVERSARIAL CONTROL for P1e's §4.2-leg-4 offline effect measurement.

`p1e_3_channel.py` (B) shows D:shipped 0.3389 -> D:common 0.2044 on capture-OFF gDNA-rich exons.  The fused
composition only moves when the LEFT and RIGHT message are damped by DIFFERENT amounts, so the claim is that
the DERIVED b2_cons is a valid per-message trust statistic.  The null hypotheses this script tests:

  const : replace b2 by ONE constant (the mean of the firing b2) on the SAME firing subset.  Removes all
          per-message information, keeps the exposure and the magnitude.  Should be ~inert if the mechanism
          really is differential.
  perm  : randomly PERMUTE the per-message b2 among the ok messages of the same edge class.  Keeps the exact
          marginal distribution of b2 and the amount of differential damping, destroys its association with
          delta.  If `perm` recovers most of the gain, b2 is not carrying information -- noise in the
          left/right weights is.
  d2    : b2 := delta^2 (no null subtraction) -- does the -sd2-nu subtraction matter at all?
  graft : damp every GRAFT edge by the pooled mean b2 and nothing else (a structural, delta-free rule).

Same offline re-fuse and the same score as `p1e_3_channel.py` (B).

    OMP_NUM_THREADS=1 python scratchpad/verify_p1e_control.py [cond ...]
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
rng = np.random.default_rng(7)
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "gdna_gdna1_ss_0.50_nrna_present_capture_on",
    "gdna_none_ss_0.50_nrna_present_capture_off",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def msg_state(ent, us, uni, M, E_g, E_r, is_bnd, is_exon, n_node, B):
    df, src, valid = ent["df"], ent["src"], ent["valid"]
    og, op, on = us["og"], us["op"], us["on"]
    face = uni["rho_lf"] if df == 0 else uni["rho_rf"]
    q = np.maximum(face / np.maximum(M * B, _EPS) - 1.0, 0.0)
    cM = q / (1.0 + q)
    tg, tp, tn = ent["tg"], ent["tp"], ent["tn"]
    p3 = np.stack([ent["tpg"], ent["tpp"], ent["tpn"]], axis=1)
    d3 = np.stack([tg, tp, tn], axis=1)
    o3 = np.stack([og, op, on], axis=1)
    E3 = np.stack([E_g, E_r, E_r], axis=1)
    sup = p3 > 0.0
    m3 = np.where(sup, d3, o3) * E3
    S = m3.sum(axis=1)
    ok = valid & (S > _EPS) & (M > _EPS)
    alpha = m3 / np.maximum(S, _EPS)[:, None]
    v3 = np.where(sup, 1.0 / np.maximum(p3, _EPS), np.inf)
    scm = np.maximum(ent["s2t"], 0.0) + 1.0 / np.maximum(ent["n_src"], _EPS)
    w3 = np.maximum(np.where(sup, v3, 0.0) - scm[:, None], 0.0)
    sd2 = scm + (alpha * alpha * w3).sum(axis=1)
    delta = np.log(np.maximum(M, _EPS) / np.maximum(S, _EPS))
    nu = cM * cM / np.maximum(n_node, _EPS)
    b2 = np.where(ok, np.maximum(0.0, delta * delta - sd2 - nu), 0.0)
    d2 = np.where(ok, delta * delta, 0.0)
    graft = ent["graft"]
    peel = is_bnd & is_exon[src] & valid
    return dict(df=df, src=src, ok=ok, b2=b2, d2=d2, p3=p3, d3=d3, sup=sup,
                graft=graft, peel=peel, plain=valid & ~graft & ~peel)


def cls_key(X):
    k = np.zeros(X["ok"].shape, np.int8)
    k[X["graft"]] = 1
    k[X["peel"]] = 2
    return k


print("mass-weighted |f_g - oracle| of the composition the exon node is HANDED (offline re-fuse)")
hdr = ("condition", "n", "shipped", "common", "const", "perm", "d2", "graft")
print(f"{hdr[0]:<44}{hdr[1]:>5}" + "".join(f"{h:>10}" for h in hdr[2:]))
print("-" * (49 + 10 * 6))
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
    chain, cap = dbg["chain"], dbg["capture"]
    us, uni = cap["_uni_static"], cap["_uni"][-1]
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    n = M.shape[0]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    n_node = np.where(~is_bnd, us["n_unspl_l"], us["n_unspl_l"] + us["n_unspl_r"])
    f_cur = np.clip(np.asarray(cap["_uni"][-2]["fg_out"], float), 0.0, 1.0)
    B = f_cur / np.maximum(E_g, _EPS) + (1.0 - f_cur) / np.maximum(E_r, _EPS)
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    CLS = {0: "intergenic", 1: "intron", 2: "exon"}
    dcls = np.array([CLS.get(int(rt[j]), "?") if kind[j] == REGION else "boundary" for j in range(n)])

    A = msg_state(cap["_pin"][-2], us, uni, M, E_g, E_r, is_bnd, is_exon, n_node, B)
    Bm = msg_state(cap["_pin"][-1], us, uni, M, E_g, E_r, is_bnd, is_exon, n_node, B)
    sel = ((dcls == "exon") & np.isfinite(fo) & (M > _EPS)
           & np.asarray(cap["solvable"], bool) & A["ok"] & Bm["ok"])

    # the pooled mean of the FIRING b2, and its per-edge-class means, over both messages
    allb = np.concatenate([A["b2"][A["ok"]], Bm["b2"][Bm["ok"]]])
    mfire = float(allb[allb > 0].mean()) if (allb > 0).any() else 0.0
    kA, kB = cls_key(A), cls_key(Bm)

    def bvec(X, k, arm):
        if arm == "common":
            return X["b2"]
        if arm == "d2":
            return X["d2"]
        if arm == "const":
            return np.where(X["b2"] > 0, mfire, 0.0)
        if arm == "graft":
            return np.where(X["graft"] & X["ok"], mfire, 0.0)
        if arm == "perm":  # permute within (edge class) among ok messages
            out = np.zeros_like(X["b2"])
            for c in (0, 1, 2):
                m = X["ok"] & (k == c)
                idx = np.flatnonzero(m)
                if idx.size:
                    out[idx] = X["b2"][rng.permutation(idx)]
            return out
        return np.zeros_like(X["b2"])

    row = {"cond": cond[5:], "n": int(sel.sum())}
    for arm in ("shipped", "common", "const", "perm", "d2", "graft"):
        pg, dens = [], []
        for X, k in ((A, kA), (Bm, kB)):
            p = X["p3"].copy()
            if arm != "shipped":
                v = np.where(X["sup"], 1.0 / np.maximum(p, _EPS), np.inf)
                bb = np.repeat(bvec(X, k, arm)[:, None], 3, axis=1)
                p = np.where(X["sup"], 1.0 / np.maximum(v + bb, _EPS), 0.0)
            pg.append(p)
            dens.append(X["d3"])
        pa, pb = pg
        cd = np.where((pa + pb) > _EPS, (pa * dens[0] + pb * dens[1]) / np.maximum(pa + pb, _EPS), 0.0)
        gm_, rm_ = cd[:, 0] * E_g, (cd[:, 1] + cd[:, 2]) * E_r
        f = np.where((gm_ + rm_) > _EPS, gm_ / np.maximum(gm_ + rm_, _EPS), np.nan)
        m2 = sel & np.isfinite(f)
        row[arm] = float(np.average(np.abs(f[m2] - fo[m2]), weights=M[m2]))
    print(f"{row['cond']:<44}{row['n']:>5}" + "".join(
        f"{row[a]:>10.4f}" for a in ("shipped", "common", "const", "perm", "d2", "graft")))
