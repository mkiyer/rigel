"""P-2 / PASS-0 — THE gDNA OVER-CALL, TRACED TO THE NODE AND DECOMPOSED TO ITS ORIGIN.

Part F measured that the gDNA density delivered into EXON nodes at capture-OFF is **5.4× too large**
(0.732 decades), that the source is an intron or boundary on **100 %** of those edges, and that dropping the
reframe `r` from the gDNA component takes the error to 0.170 dec. This is the forensic version: which nodes,
what they look like, what the raw counts are, and where the number actually comes from.

⭐ **THE DECOMPOSITION IS EXACT, NOT A MODEL.** The delivered gDNA density is `tg = rg_src · r`, so

    log10( tg / ρ_g^true(dst) )  =  log10( rg_src / ρ_g^true(src) )   ← the SOURCE's own error
                                 +  log10( r )                        ← the REFRAME
                                 +  log10( ρ_g^true(src) / ρ_g^true(dst) )   ← TRUE spatial difference

Every term is separately observable. gDNA is uniform in these libraries, so the third term is ~0 by
construction, which means the whole error must live in the first two — and the point of the table is to say
WHICH.

Run: OMP_NUM_THREADS=1 python scratchpad/p2r_g_trace.py [--cond ...] [--top 6]
"""
from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-12
CLS = ("intergenic", "intron", "exon", "boundary")
BITS = [(0x8, "intron+"), (0x4, "intron-"), (0x2, "exon+"), (0x1, "exon-")]


def sig_str(s):
    s = int(s)
    on = [n for b, n in BITS if s & b]
    return f"0x{s:X}[{','.join(on) if on else 'intergenic'}]"


def node_loc(i, chain, ix, ra):
    """(label, ref, start, end) for a chain node — region interval or boundary position."""
    k = int(np.asarray(chain.kind)[i])
    r = int(np.asarray(chain.ref_idx)[i])
    if k == REGION:
        row = ix.region_df.iloc[r]
        return (f"region {r}", row["ref_name"], int(row["start"]), int(row["end"]),
                int(row["signature"]))
    row = ix.boundary_df.iloc[r]
    return (f"boundary {r}", row["ref_name"], int(row["position"]), int(row["position"]), -1)


def transcripts_over(ix, ref, start, end):
    t = ix.t_df
    m = (t["ref"] == ref) & (t["start"] < max(end, start + 1)) & (t["end"] > start)
    return t[m]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cond", default="gdna_gdna100_ss_0.50_nrna_none_capture_off")
    ap.add_argument("--top", type=int, default=6)
    a = ap.parse_args()

    ix = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(ix.region_df, ix.ref_name_to_id)
    inp = _scan_and_truth(SUITE, a.cond, ix, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]),
              dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    cap, chain, st = dbg["capture"], dbg["chain"], dbg["statics"]
    S = cap["_uni_static"]
    E_g, E_r, M = (np.asarray(S[k], float) for k in ("E_g", "E_r", "M"))
    og, op, on = (np.asarray(S[k], float) for k in ("og", "op", "on"))
    rl, rr = np.asarray(S["rho_l0"], float), np.asarray(S["rho_r0"], float)
    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)
    fg = np.asarray(cap["f_g"], float)
    rt, _ = _node_region_type(chain, ra)
    cls = np.where(np.asarray(chain.kind) != REGION, 3, rt)
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    rho_true = G / np.maximum(E_g, _EPS)                      # the TRUE gDNA density at every node
    up, un = np.asarray(st.u_pos, float), np.asarray(st.u_neg, float)
    tau0 = np.asarray(cap["_tau0_lam"], float)
    MG, PG = np.asarray(cap["mode_g"], float), np.asarray(cap["prec_g"], float)
    MP, PP = np.asarray(cap["mode_p"], float), np.asarray(cap["prec_p"], float)
    MN, PN = np.asarray(cap["mode_n"], float), np.asarray(cap["prec_n"], float)

    # the two per-message states: (capture record, dst face, src face)
    msgs = [(cap["_pin"][-2], rl, rr, "LEFT"), (cap["_pin"][-1], rr, rl, "RIGHT")]

    print(f"# condition = {a.cond}")
    print(f"# {int((cls == 2).sum())} exon nodes, {int((cls == 1).sum())} intron, "
          f"{int((cls == 3).sum())} boundary, {int((cls == 0).sum())} intergenic\n")

    # ── rank exon destinations by how much gDNA over-call mass they carry ────────────────────────────
    rows = []
    for rec, dstf, srcf, tag in msgs:
        src = np.asarray(rec["src"], np.int64)
        valid = np.asarray(rec["valid"], bool)
        tg, tpg = np.asarray(rec["tg"], float), np.asarray(rec["tpg"], float)
        framed = valid & (srcf[src] > _EPS) & (dstf > _EPS)
        r = np.where(framed, dstf / np.maximum(srcf[src], _EPS), 1.0)
        sel = (valid & (tpg > 0) & (tg > _EPS) & (cls == 2) & (rho_true > _EPS) & (mass > _EPS))
        for i in np.nonzero(sel)[0]:
            i = int(i)
            e = np.log10(tg[i] / rho_true[i])
            rows.append(dict(i=i, tag=tag, src=int(src[i]), r=float(r[i]), tg=float(tg[i]),
                             tpg=float(tpg[i]), err=float(e), w=float(mass[i] * abs(e))))
    rows.sort(key=lambda d: -d["w"])

    # ── THE DECOMPOSITION, over all such edges ───────────────────────────────────────────────────────
    print("=" * 108)
    print("THE EXACT DECOMPOSITION of the delivered gDNA error, over every exon destination with a live")
    print("gDNA message.  log10(tg/rho_true_dst) = [source's own error] + [log10 r] + [true spatial diff]")
    print("=" * 108)
    tot = np.array([r_["w"] for r_ in rows])
    def col(f):
        v = np.array([f(r_) for r_ in rows])
        return np.average(v, weights=np.maximum(tot, _EPS)), np.median(v)
    src_err = []
    reframe = []
    spatial = []
    for r_ in rows:
        i, s = r_["i"], r_["src"]
        rg_src = r_["tg"] / max(r_["r"], _EPS)          # the relayed context gDNA density AT THE SOURCE
        se = np.log10(max(rg_src, _EPS) / max(rho_true[s], _EPS)) if rho_true[s] > _EPS else np.nan
        src_err.append(se)
        reframe.append(np.log10(max(r_["r"], _EPS)))
        spatial.append(np.log10(max(rho_true[s], _EPS) / max(rho_true[i], _EPS))
                       if rho_true[s] > _EPS else np.nan)
    src_err, reframe, spatial = np.array(src_err), np.array(reframe), np.array(spatial)
    tot_e = np.array([r_["err"] for r_ in rows])
    ok = np.isfinite(src_err) & np.isfinite(spatial)
    w = np.maximum(tot, _EPS)[ok]
    print(f"\n{'term':44s} {'mass-wt mean':>13s} {'median':>9s} {'|mean|':>9s}")
    for nm, v in (("TOTAL  log10(delivered / true)", tot_e),
                  ("  (1) the SOURCE's own gDNA error", src_err),
                  ("  (2) the REFRAME   log10(r)", reframe),
                  ("  (3) TRUE spatial difference", spatial)):
        vv = v[ok]
        print(f"{nm:44s} {np.average(vv, weights=w):13.3f} {np.median(vv):9.3f} "
              f"{np.average(np.abs(vv), weights=w):9.3f}")
    resid = tot_e[ok] - (src_err[ok] + reframe[ok] + spatial[ok])
    print(f"{'  identity residual (must be ~0)':44s} {np.abs(resid).max():13.2e}")
    print(f"\n  n edges = {ok.sum()}   (gDNA is uniform in this library, so term (3) SHOULD be ~0)")

    # ── the worst nodes, in full ─────────────────────────────────────────────────────────────────────
    seen = set()
    shown = 0
    for r_ in rows:
        if shown >= a.top:
            break
        if r_["i"] in seen:
            continue
        seen.add(r_["i"])
        shown += 1
        i, s = r_["i"], r_["src"]
        lbl, ref, x0, x1, sg = node_loc(i, chain, ix, ra)
        slbl, sref, sx0, sx1, ssg = node_loc(s, chain, ix, ra)
        print("\n" + "=" * 108)
        print(f"NODE {i}  [{CLS[int(cls[i])]}]  {lbl}  {ref}:{x0:,}-{x1:,}  len={x1 - x0:,}  "
              f"sig={sig_str(sg)}")
        print("=" * 108)

        print("\n-- TRANSCRIPT STRUCTURE over this interval " + "-" * 60)
        tt = transcripts_over(ix, ref, x0, x1)
        if tt.empty:
            print("   (no transcript overlaps — intergenic)")
        for _, t in tt.iterrows():
            ex = ix.get_exon_intervals(int(t["t_index"]))
            tags = []
            if ex is not None:
                for (a0, b0) in ex:
                    if b0 > x0 and a0 < x1:
                        tags.append(f"EXON {a0:,}-{b0:,}")
                for k in range(len(ex) - 1):
                    a0, b0 = int(ex[k][1]), int(ex[k + 1][0])
                    if b0 > x0 and a0 < x1:
                        tags.append(f"intron {a0:,}-{b0:,}")
            print(f"   {t['t_id']:10s} {t['g_name']:12s} strand={'+' if t['strand'] == 1 else '-'} "
                  f"span={t['start']:,}-{t['end']:,} n_exons={t['n_exons']} "
                  f"nrna={bool(t['is_nrna'])}")
            print(f"      exons: {[(int(p), int(q)) for p, q in ex] if ex is not None else None}")
            print(f"      overlapping this region: {'; '.join(tags) if tags else 'span only'}")

        print("\n-- RAW COUNTS AND DENSITIES AT THIS NODE " + "-" * 62)
        print(f"   unspliced  u_pos={up[i]:,.1f}  u_neg={un[i]:,.1f}   mass M={M[i]:,.1f}   "
              f"eff E_g={E_g[i]:,.1f}  E_r={E_r[i]:,.1f}")
        print(f"   ORACLE     gDNA={G[i]:,.1f}  RNA={R[i]:,.1f}   true f_g={fo[i]:.4f}   "
              f"true rho_g={rho_true[i]:.6f}   total rho={M[i] / max(eff[i], _EPS):.6f}")
        print(f"   OWN solve  og={og[i]:.6f} op={op[i]:.6f} on={on[i]:.6f}  "
              f"f_g_own={og[i] / max(og[i] + op[i] + on[i], _EPS):.4f}  tau0_lam={tau0[i]:.4g} "
              f"({'NO' if tau0[i] <= 0 else 'has'} intrinsic strand evidence)")
        print(f"   SOLVED     f_g={fg[i]:.4f}   vs oracle {fo[i]:.4f}   "
              f"|err|={abs(fg[i] - fo[i]):.4f}   err mass={mass[i] * abs(fg[i] - fo[i]):,.0f}")

        print("\n-- THE INCOMING MESSAGES " + "-" * 78)
        for rec, dstf, srcf, tag in msgs:
            sv = int(np.asarray(rec["src"])[i])
            if not bool(np.asarray(rec["valid"], bool)[i]):
                print(f"   {tag:5s}  (no message)")
                continue
            tg = float(np.asarray(rec["tg"], float)[i])
            tp = float(np.asarray(rec["tp"], float)[i])
            tpg_ = float(np.asarray(rec["tpg"], float)[i])
            tpp_ = float(np.asarray(rec["tpp"], float)[i])
            rv = (float(dstf[i]) / max(float(srcf[sv]), _EPS)
                  if (srcf[sv] > _EPS and dstf[i] > _EPS) else 1.0)
            sl2, sref2, sa, sb, ssg2 = node_loc(sv, chain, ix, ra)
            rg_src = tg / max(rv, _EPS)
            print(f"   {tag:5s}  from node {sv} [{CLS[int(cls[sv])]}] {sl2} {sref2}:{sa:,}-{sb:,} "
                  f"sig={sig_str(ssg2)}")
            print(f"          SOURCE truth: gDNA={G[sv]:,.1f} RNA={R[sv]:,.1f} f_g={fo[sv]:.4f} "
                  f"rho_g={rho_true[sv]:.6f}   SOURCE solved f_g={fg[sv]:.4f}")
            print(f"          reframe r = rho_face(dst)/rho_face(src) = {float(dstf[i]):.6f}/"
                  f"{float(srcf[sv]):.6f} = {rv:.3f}")
            print(f"          relayed gDNA at source rg_src={rg_src:.6f}  ->  x r  ->  "
                  f"DELIVERED tg={tg:.6f}   (prec {tpg_:.4g})")
            print(f"          delivered RNA tp={tp:.6f} (prec {tpp_:.4g})")
            if rho_true[sv] > _EPS:
                print(f"          ORIGIN:  source's own error {np.log10(max(rg_src, _EPS) / rho_true[sv]):+.3f}"
                      f"  +  reframe {np.log10(max(rv, _EPS)):+.3f}"
                      f"  +  true spatial {np.log10(rho_true[sv] / max(rho_true[i], _EPS)):+.3f}"
                      f"  =  {np.log10(tg / max(rho_true[i], _EPS)):+.3f} dec "
                      f"({tg / max(rho_true[i], _EPS):.1f}x)")

        print("\n-- WHAT psi RECEIVES " + "-" * 82)
        print(f"   gDNA level  mode f_g={np.exp(MG[i]):.4f} (prec {PG[i]:.4g})     "
              f"[oracle f_g={fo[i]:.4f}]")
        print(f"   RNA+ level  mode f_R={np.exp(MP[i]):.4f} (prec {PP[i]:.4g})     "
              f"RNA- mode={np.exp(MN[i]):.4f} (prec {PN[i]:.4g})")
        print(f"   => solved f_g = {fg[i]:.4f}")


if __name__ == "__main__":
    main()
