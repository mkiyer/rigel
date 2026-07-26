"""RC2 (A) — WHERE does the gDNA-message mode bias come from?  Exact re-implementation of `_transport`,
staged.

`bp_solver._transport` is rebuilt here verbatim from the diagnostic capture and VERIFIED against the shipped
`ag/bg/cg/mo_g` before any attribution. Then the message's implied f_g is read out at each stage of its
construction, so the bias `mo_g − oracle` can be split into:

    S0  the SOURCE neighbour's OWN message-free self-solve f_g        (source-belief error)
    S1  the relay-accumulated composition share AT the source          (+ the per-hop `_fuse` drift)
    S2  after GRAFT / PEEL routing + the strand filter                 (+ routing)
    S3  after the REFRAME by r                                         (+ scale; share-invariant by algebra)
    S4  after `_pin_v`                                                 (+ the own-density substitution)
    S5  after the left/right per-component mode fusion  == mo_g        (+ the asymmetric combine)

    OMP_NUM_THREADS=1 python scratchpad/rc2_a_decomp.py [COND ...]
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.enrichment_frame import (  # noqa: E402
    mismatch_deflate,
    mismatch_gap,
    transfer_logvar,
)
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.node_init import own_composition_logvar  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
_EPS = 1e-9
BINS = [(0.0, 0.30), (0.30, 0.60), (0.60, 0.90), (0.90, 0.99), (0.99, 1.01)]

ap = argparse.ArgumentParser()
ap.add_argument("conds", nargs="*", default=["gdna_gdna300_ss_0.99_nrna_present_capture_on",
                                             "gdna_gdna300_ss_0.99_nrna_none_capture_on",
                                             "gdna_gdna100_ss_0.99_nrna_present_capture_on"])
a = ap.parse_args()

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
ACC: dict[str, list] = {}


def run(cond: str):
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]),
                     dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    chain, cap, st, geo = dbg["chain"], dbg["capture"], dbg["statics"], dbg["geometry"]
    uni, us = cap["_uni"][-1], cap["_uni_static"]
    n = np.asarray(cap["f_g"]).shape[0]

    # ── the statics `_unified_solve` builds ──
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    og, op, on = us["og"], us["op"], us["on"]
    is_bnd, ex_a = us["is_bnd"].astype(bool), us["is_exon"].astype(bool)
    is_reg = ~is_bnd
    fp_a, fn_a = us["fp"].astype(bool), us["fn"].astype(bool)
    li, ri = us["left"], us["right"]
    logvar_tot = us["logvar_tot"]
    tau_own, struct = us["tau_own"], us["struct_lock"].astype(bool)
    ni_fg = np.asarray(cap["fg_loc"])
    v_own_g, v_own_r = own_composition_logvar(ni_fg, tau_own, struct)
    v_own_lam = np.where(struct, 0.0, np.where(tau_own > _EPS, 1.0 / np.maximum(tau_own, _EPS), np.inf))

    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    SP = (us["SP_l"], us["SP_r"])
    SN = (us["SN_l"], us["SN_r"])
    spl_p_f = tuple(np.where(SP[k] > _EPS, SP[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
    spl_n_f = tuple(np.where(SN[k] > _EPS, SN[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))

    vl, vr = li >= 0, ri >= 0
    sl, sr = np.clip(li, 0, n - 1), np.clip(ri, 0, n - 1)

    def pin_v(g, p, nn, pg_, pp_, pn_):
        sg = np.where(pg_ > 0.0, g, og)
        sp = np.where(pp_ > 0.0, p, op)
        sn = np.where(pn_ > 0.0, nn, on)
        s = sg * E_g + (sp + sn) * E_r
        k = np.where((s > _EPS) & (M > _EPS), M / np.maximum(s, _EPS), 1.0)
        return g * k, p * k, nn * k, k

    def transport(src, valid, df, sf, arrs, dst_face_v, src_face_v):
        """Verbatim `bp_solver._transport`, plus the per-stage implied-f_g readouts."""
        rg, rp, rn, pg, pp, pn, mg, mp, mn, tau = arrs
        framed = valid & (src_face_v[src] > _EPS) & (dst_face_v > _EPS)
        r = np.where(framed, dst_face_v / np.maximum(src_face_v[src], _EPS), np.where(valid, 1.0, 0.0))
        graft = ex_a & is_bnd[src] & valid
        gp = np.where(graft, spl_p_f[sf][src], 0.0)
        gn = np.where(graft, spl_n_f[sf][src], 0.0)
        # ── S1: the relay-accumulated share AT THE SOURCE (source geometry) ──
        s1 = _share(rg[src], rp[src] + rn[src], E_g[src], E_r[src])
        # ── S2: after graft (pre-reframe, source frame) ──
        s2 = _share(rg[src], rp[src] + gp + rn[src] + gn, E_g[src], E_r[src])
        tg, tp, tn = rg[src] * r, (rp[src] + gp) * r, (rn[src] + gn) * r
        s2t = transfer_logvar(logvar_tot, logvar_tot[src], graft)

        def _dv(p_, s2_=s2t):
            return np.where(valid & (p_[src] > 0.0), 1.0 / (1.0 / np.maximum(p_[src], _EPS) + s2_), 0.0)

        tpg, tpp, tpn = _dv(pg), _dv(pp), _dv(pn)
        tmg, tmp, tmn = _dv(mg), _dv(mp), _dv(mn)
        ttau = _dv(tau, s2t)
        _sp, _sn = np.where(graft, SP[sf][src], 0.0), np.where(graft, SN[sf][src], 0.0)
        _s2t_spl = np.where(np.isfinite(s2t), s2t, 0.0)
        _spc = np.where(_sp > _EPS, _sp / (1.0 + _sp * _s2t_spl), 0.0)
        _snc = np.where(_sn > _EPS, _sn / (1.0 + _sn * _s2t_spl), 0.0)
        tpp, tpn = tpp + _spc, tpn + _snc
        tmp, tmn = tmp + _spc, tmn + _snc
        peel = is_bnd & ex_a[src] & valid
        tp = np.where(peel, np.maximum(tp - spl_p_f[df], 0.0), tp)
        tn = np.where(peel, np.maximum(tn - spl_n_f[df], 0.0), tn)
        tp, tpp, tmp = np.where(fp_a, tp, 0.0), np.where(fp_a, tpp, 0.0), np.where(fp_a, tmp, 0.0)
        tn, tpn, tmn = np.where(fn_a, tn, 0.0), np.where(fn_a, tpn, 0.0), np.where(fn_a, tmn, 0.0)
        # ── S3: after peel/filter + reframe, DESTINATION geometry, pre-pin ──
        s3 = _share(tg, tp + tn, E_g, E_r)
        pre = (tg.copy(), tp.copy(), tn.copy())
        tg, tp, tn, kpin = pin_v(tg, tp, tn, tpg, tpp, tpn)
        # ── S4: the per-message implied f_g actually delivered (post-pin) ──
        s4 = tg * E_g / np.maximum(M, _EPS)
        g_g, c_g = mismatch_gap(tg, og)
        g_p, c_p = mismatch_gap(tp, op)
        g_n, c_n = mismatch_gap(tn, on)
        tpg = mismatch_deflate(tpg, g_g, c_g, v_own_g)
        tpp = mismatch_deflate(tpp, g_p, c_p, v_own_r)
        tpn = mismatch_deflate(tpn, g_n, c_n, v_own_r)
        tmg = mismatch_deflate(tmg, g_g, c_g, v_own_g)
        tmp = mismatch_deflate(tmp, g_p, c_p, v_own_r)
        tmn = mismatch_deflate(tmn, g_n, c_n, v_own_r)
        ttau = np.where((tg > _EPS) & ((tp + tn) > _EPS), ttau, 0.0)
        g_R, c_R = mismatch_gap(tp + tn, op + on)
        ttau = mismatch_deflate(ttau, g_g - g_R, c_g | c_R, v_own_lam)
        return dict(tg=tg, tp=tp, tn=tn, tpg=tpg, tpp=tpp, tpn=tpn, tmg=tmg, tmp=tmp, tmn=tmn,
                    ttau=ttau, s1=s1, s2=s2, s3=s3, s4=s4, r=r, kpin=kpin, graft=graft, peel=peel,
                    pre_g=pre[0], pre_p=pre[1], pre_n=pre[2], valid=valid, src=src)

    def _share(g, rr, eg, er):
        num = g * eg
        den = num + rr * er
        return np.where(den > _EPS, num / np.maximum(den, _EPS), np.nan)

    fwd = (us["fwd_g"], us["fwd_p"], us["fwd_n"], us["fwd_pg"], us["fwd_pp"], us["fwd_pn"],
           us["fwd_mg"], us["fwd_mp"], us["fwd_mn"], us["fwd_tau"])
    bwd = (us["bwd_g"], us["bwd_p"], us["bwd_n"], us["bwd_pg"], us["bwd_pp"], us["bwd_pn"],
           us["bwd_mg"], us["bwd_mp"], us["bwd_mn"], us["bwd_tau"])
    A = transport(sl, vl, 0, 1, fwd, uni["rho_lf"], uni["rho_rf"])
    B = transport(sr, vr, 1, 0, bwd, uni["rho_rf"], uni["rho_lf"])

    # ── FIDELITY ──
    cpg_ = A["tpg"] + B["tpg"]
    cg = np.where(cpg_ > _EPS, (A["tpg"] * A["tg"] + B["tpg"] * B["tg"]) / np.maximum(cpg_, _EPS), 0.0)
    mo_g = np.log(np.maximum(cg * E_g / np.maximum(M, _EPS), _EPS))
    f1 = np.max(np.abs(A["tg"] - uni["ag"]))
    f2 = np.max(np.abs(B["tg"] - uni["bg"]))
    f3 = np.max(np.abs(mo_g - uni["mo_g"]))
    f4 = np.max(np.abs(A["tmg"] + B["tmg"] - uni["cm_g"]))

    # ── the oracle + the per-node frame ──
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    cls = np.array([CLS[int(rt[i])] if kind[i] == REGION else "boundary" for i in range(n)])
    mass = np.asarray(cap["mass_global"])
    print(f"\n### {cond}")
    print(f"  transport fidelity  |ag|={f1:.2e} |bg|={f2:.2e} |mo_g|={f3:.2e} |cm_g|={f4:.2e}"
          f"   {'OK' if max(f1, f2, f3, f4) < 1e-9 else '*** MISMATCH ***'}")

    ACC.setdefault("rows", []).append(dict(
        cond=cond, fo=fo, cls=cls, mass=mass, mo_g=np.exp(mo_g), tau_own=tau_own, ni_fg=ni_fg,
        cm_g=uni["cm_g"], A=A, B=B, sl=sl, sr=sr, vl=vl, vr=vr, fg=np.asarray(cap["f_g"]),
        E_g=E_g, E_r=E_r, M=M, og=og, op=op, on=on, cg=cg, cpg=cpg_,
        cp=uni["cp"], cn=uni["cn"], cpp=uni["pp"], cpn=uni["pn"],
    ))


for c in a.conds:
    run(c)

# ────────────────────────────────────────────────────────────────────────────────────────────────────────
# THE STAGED ATTRIBUTION, pooled over the conditions, on FULL-RANK EXONS with a live gDNA message.
# For each side (left/right) we track the message's implied f_g through S0..S4, then the combine → mo_g.
# ────────────────────────────────────────────────────────────────────────────────────────────────────────
cat: dict[str, list] = {}


def push(**kw):
    for k, v in kw.items():
        cat.setdefault(k, []).append(np.asarray(v, np.float64))


for R_ in ACC["rows"]:
    fo, mass, tau_own, cls = R_["fo"], R_["mass"], R_["tau_own"], R_["cls"]
    for side, T, s_idx, val in (("L", R_["A"], R_["sl"], R_["vl"]), ("R", R_["B"], R_["sr"], R_["vr"])):
        m = (np.isfinite(fo) & (mass > 0) & (tau_own > 0) & (cls == "exon") & val
             & (T["tpg"] > 0) & (R_["cm_g"] > 0))
        if not m.any():
            continue
        s = s_idx[m]
        push(oracle=fo[m], mass=mass[m], own=R_["ni_fg"][m], src_oracle=fo[s], src_own=R_["ni_fg"][s],
             s1=T["s1"][m], s2=T["s2"][m], s3=T["s3"][m], s4=T["s4"][m], mo_g=R_["mo_g"][m],
             r=T["r"][m], kpin=T["kpin"][m], graft=T["graft"][m].astype(float),
             peel=T["peel"][m].astype(float),
             # was the message PARTIAL (an RNA arm not supplied ⇒ the pin substitutes the dst's own)?
             part_p=(T["tpp"][m] <= 0).astype(float), part_n=(T["tpn"][m] <= 0).astype(float),
             src_is_bnd=np.ones(int(m.sum())))
D = {k: np.concatenate(v) for k, v in cat.items()}
w = D["mass"]
print(f"\n{'=' * 150}\nSTAGED ATTRIBUTION of the per-side gDNA message, full-rank EXONS with a live gDNA "
      f"message  (n={w.size:,} messages)\n{'=' * 150}")
print(f"{'oracle bin':<13}{'n':>6}{'mass':>12}{'oracle':>8}{'SRC oracle':>11}{'S0 src own':>11}"
      f"{'S1 relay':>10}{'S2 graft':>10}{'S3 reframe':>11}{'S4 pin':>9}{'S5 mo_g':>9}"
      f"{'  |  dS1':>9}{'dS2':>8}{'dS3':>8}{'dS4':>8}{'dS5':>8}")
for lo, hi in BINS:
    m = (D["oracle"] >= lo) & (D["oracle"] < hi)
    if not m.any():
        continue
    ww = w[m]

    def av(k):
        return np.average(D[k][m], weights=ww)

    s0, s1, s2, s3, s4, s5 = av("src_own"), av("s1"), av("s2"), av("s3"), av("s4"), av("mo_g")
    print(f"[{lo:.2f},{hi:.2f})  {int(m.sum()):>6}{ww.sum():>12,.0f}{av('oracle'):>8.3f}"
          f"{av('src_oracle'):>11.3f}{s0:>11.3f}{s1:>10.3f}{s2:>10.3f}{s3:>11.3f}{s4:>9.3f}{s5:>9.3f}"
          f"{s1 - s0:>+12.3f}{s2 - s1:>+8.3f}{s3 - s2:>+8.3f}{s4 - s3:>+8.3f}{s5 - s4:>+8.3f}")
print(f"\n{'oracle bin':<13}{'graft%':>8}{'peel%':>8}{'partial +':>10}{'partial -':>10}"
      f"{'median r':>10}{'p90 r':>10}{'median kpin':>13}{'p10/p90 kpin':>18}")
for lo, hi in BINS:
    m = (D["oracle"] >= lo) & (D["oracle"] < hi)
    if not m.any():
        continue
    print(f"[{lo:.2f},{hi:.2f})  {D['graft'][m].mean():>8.1%}{D['peel'][m].mean():>8.1%}"
          f"{D['part_p'][m].mean():>10.1%}{D['part_n'][m].mean():>10.1%}"
          f"{np.median(D['r'][m]):>10.3f}{np.quantile(D['r'][m], 0.9):>10.3f}"
          f"{np.median(D['kpin'][m]):>13.3f}"
          f"{np.quantile(D['kpin'][m], 0.1):>9.3f}/{np.quantile(D['kpin'][m], 0.9):<9.3f}")

# ────────────────────────────────────────────────────────────────────────────────────────────────────────
# WHAT the two dominant stages are made of, in ABSOLUTE RNA-fraction terms at the destination.
#   relay_rna  = (rp+rn)[src]·r·E_r/M   — the imputed RNA the source relayed
#   graft_rna  = (gp+gn)·r·E_r/M        — the boundary's MEASURED mature spliced density, reframed
#   deliv_rna  = (tp+tn)·E_r/M          — what the message finally claims (post-peel/filter/pin)
# against the destination's TRUE RNA fraction 1−oracle.
# ────────────────────────────────────────────────────────────────────────────────────────────────────────
cat2: dict[str, list] = {}
for R_ in ACC["rows"]:
    fo, mass, tau_own, cls = R_["fo"], R_["mass"], R_["tau_own"], R_["cls"]
    E_g, E_r, M, og = R_["E_g"], R_["E_r"], R_["M"], R_["og"]
    for T, s_idx, val in ((R_["A"], R_["sl"], R_["vl"]), (R_["B"], R_["sr"], R_["vr"])):
        m = (np.isfinite(fo) & (mass > 0) & (tau_own > 0) & (cls == "exon") & val
             & (T["tpg"] > 0) & (R_["cm_g"] > 0))
        if not m.any():
            continue
        s = s_idx[m]
        for k, v in dict(
            oracle=fo[m], mass=mass[m], src_oracle=fo[s], src_own=R_["ni_fg"][s],
            s1=T["s1"][m],
            relay_rna=(T["pre_p"][m] + T["pre_n"][m]) * E_r[m] / np.maximum(M[m], _EPS),
            deliv_rna=(T["tp"][m] + T["tn"][m]) * E_r[m] / np.maximum(M[m], _EPS),
            deliv_g=T["tg"][m] * E_g[m] / np.maximum(M[m], _EPS),
            own_g=og[m] * E_g[m] / np.maximum(M[m], _EPS),
            true_rna=1.0 - fo[m],
        ).items():
            cat2.setdefault(k, []).append(np.asarray(v, np.float64))
D2 = {k: np.concatenate(v) for k, v in cat2.items()}
w2 = D2["mass"]
print(f"\n{'=' * 150}\nABSOLUTE RNA CLAIM of the message vs the destination's TRUTH (full-rank exons)\n{'=' * 150}")
print(f"{'oracle bin':<13}{'n':>6}{'TRUE rna frac':>15}{'relay rna':>11}{'delivered rna':>15}"
      f"{'delivered g':>13}{'dst own g':>11}{'  |  SRC oracle':>16}{'SRC own':>9}{'S1 relay share':>16}"
      f"{'  relay err vs SRC':>19}{'own err vs SRC':>16}")
for lo, hi in BINS:
    m = (D2["oracle"] >= lo) & (D2["oracle"] < hi)
    if not m.any():
        continue
    ww = w2[m]

    def av2(k):
        return np.average(D2[k][m], weights=ww)

    print(f"[{lo:.2f},{hi:.2f})  {int(m.sum()):>6}{av2('true_rna'):>15.4f}{av2('relay_rna'):>11.3f}"
          f"{av2('deliv_rna'):>15.3f}{av2('deliv_g'):>13.3f}{av2('own_g'):>11.3f}"
          f"{av2('src_oracle'):>19.3f}{av2('src_own'):>9.3f}{av2('s1'):>16.3f}"
          f"{av2('s1') - av2('src_oracle'):>+19.3f}{av2('src_own') - av2('src_oracle'):>+16.3f}")

# ── per-node stage walk for the anchor nodes ──
print(f"\n{'=' * 150}\nPER-NODE STAGE WALK (left/right message separately)\n{'=' * 150}")
WANT = {"gdna_gdna300_ss_0.99_nrna_none_capture_on": [2197, 1909, 3237, 2173, 1523, 2683]}
for R_ in ACC["rows"]:
    ids = WANT.get(R_["cond"])
    if not ids:
        continue
    print(f"\n# {R_['cond']}")
    for i in ids:
        fo, ni = R_["fo"], R_["ni_fg"]
        print(f"  node {i} [{R_['cls'][i]}] mass={R_['mass'][i]:,.0f} oracle={fo[i]:.4f} "
              f"own={ni[i]:.4f} SOLVED={R_['fg'][i]:.4f}  mo_g={R_['mo_g'][i]:.4f} "
              f"cm_g={R_['cm_g'][i]:.4g}  E_g={R_['E_g'][i]:,.0f} E_r={R_['E_r'][i]:,.0f} M={R_['M'][i]:,.0f}")
        for side, T, sx, vv in (("L", R_["A"], R_["sl"], R_["vl"]), ("R", R_["B"], R_["sr"], R_["vr"])):
            if not vv[i]:
                continue
            s = int(sx[i])
            rr = (T["pre_p"][i] + T["pre_n"][i]) * R_["E_r"][i] / max(R_["M"][i], _EPS)
            dr = (T["tp"][i] + T["tn"][i]) * R_["E_r"][i] / max(R_["M"][i], _EPS)
            print(f"     {side}: src={s}[{R_['cls'][s]}] src_oracle={fo[s]:.4f} src_own={ni[s]:.4f} "
                  f"| S1={T['s1'][i]:.4f} S2={T['s2'][i]:.4f} S3={T['s3'][i]:.4f} S4={T['s4'][i]:.4f} "
                  f"| r={T['r'][i]:.4g} kpin={T['kpin'][i]:.4g} graft={bool(T['graft'][i])} "
                  f"| rna: relayed={rr:.4f} delivered={dr:.4f} (TRUE {1 - fo[i]:.4f}) "
                  f"| prec g={T['tpg'][i]:.4g} p={T['tpp'][i]:.4g} n={T['tpn'][i]:.4g} "
                  f"mg={T['tmg'][i]:.4g} mp={T['tmp'][i]:.4g} mn={T['tmn'][i]:.4g}")
