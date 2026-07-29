"""UNIFIED-SOLVER MESSAGE AUDIT — the fused per-component densities vs the TRUTH, per node class.

The unified solver (`RIGEL_UNIFIED=1`) washes ``f_g`` toward ~0.3 regardless of truth
(`docs/calibration/archive/UNIFIED_SOLVER_HANDOFF.md` §4). This probe answers the §5 question: **is the message
wrong, or does the ψ solve mis-weight a right message?** It dumps, per node, the fused message's implied
composition ``(f_g, f_p, f_n)^msg = (cg·E_g, cp·E_r, cn·E_r)/M`` alongside the oracle's, plus the own-only
seed, the left/right half-messages, the reframe ratios, and the precisions — so each candidate cause of
`§5` can be read off directly instead of inferred.

    RIGEL_UNIFIED=1 OMP_NUM_THREADS=1 python scripts/debug/unified_message_audit.py [cond] [--npz out.npz]
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-12


def run(cond: str):
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_index(index)
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    chain = dbg["chain"]
    cap = dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    cls = np.where(kind != 0, 3, rt)  # 0 intergenic, 1 intron, 2 exon, 3 boundary
    return dict(cap=cap, chain=chain, cls=cls, Gp=Gp, Gn=Gn, Rp=Rp, Rn=Rn, inp=inp,
                mass=np.asarray(cap["mass_global"]))


_NAMES = {0: "intergenic", 1: "intron", 2: "exon", 3: "boundary"}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cond = args[0] if args else "gdna_gdna300_ss_0.50_nrna_present_capture_on"
    d = run(cond)
    cap, cls = d["cap"], d["cls"]
    st = cap["_uni_static"]
    iters = cap["_uni"]
    M, E_g, E_r = st["M"], st["E_g"], st["E_r"]
    G = d["Gp"] + d["Gn"]
    Rp, Rn = d["Rp"], d["Rn"]
    tot = G + Rp + Rn
    ok = (tot > 1e-9) & (M > 1e-9)
    # TRUE composition (of the unspliced mass) — the target the message factor should reproduce
    tf_g, tf_p, tf_n = G / np.maximum(tot, _EPS), Rp / np.maximum(tot, _EPS), Rn / np.maximum(tot, _EPS)

    print(f"# cond={cond}  nodes={M.size}  iters={len(iters)}")
    print(f"# mass identity: median |M/(G+R) - 1| = {np.median(np.abs(M[ok] / tot[ok] - 1.0)):.4f}")

    def imp(c, E):
        return c * E / np.maximum(M, _EPS)

    for it, rec in enumerate(iters):
        mg, mp, mn = imp(rec["cg"], E_g), imp(rec["cp"], E_r), imp(rec["cn"], E_r)
        print(f"\n=== iter {it + 1}: fused-message implied composition vs TRUTH (mass-weighted) ===")
        print(f"{'class':<11}{'n':>6}{'msg_fg':>9}{'true_fg':>9}{'msg_fr':>9}{'true_fr':>9}"
              f"{'SUMf':>8}{'solved':>8}{'pr_g':>9}{'pr_r':>9}")
        for c in (0, 1, 2, 3):
            m = ok & (cls == c)
            if not m.sum():
                continue
            w = M[m]

            def A(x):
                return float(np.average(x[m], weights=w))

            print(f"{_NAMES[c]:<11}{int(m.sum()):>6}{A(mg):>9.3f}{A(tf_g):>9.3f}"
                  f"{A(mp + mn):>9.3f}{A(tf_p + tf_n):>9.3f}{A(mg + mp + mn):>8.3f}"
                  f"{A(rec['fg_out']):>8.3f}{A(rec['pg']):>9.2f}{A(rec['pp'] + rec['pn']):>9.2f}")

    # ---- own-only (no message) implied composition — is the SEED already wrong? ----
    og, op, on = st["og"], st["op"], st["on"]
    print("\n=== OWN seed (message-free) implied composition vs TRUTH ===")
    print(f"{'class':<11}{'own_fg':>9}{'true_fg':>9}{'own_fr':>9}{'true_fr':>9}{'SUMf':>8}"
          f"{'p_own_g':>9}{'p_own_r':>9}")
    for c in (0, 1, 2, 3):
        m = ok & (cls == c)
        if not m.sum():
            continue
        w = M[m]

        def A(x):
            return float(np.average(x[m], weights=w))

        print(f"{_NAMES[c]:<11}{A(imp(og, E_g)):>9.3f}{A(tf_g):>9.3f}{A(imp(op + on, E_r)):>9.3f}"
              f"{A(tf_p + tf_n):>9.3f}{A(imp(og, E_g) + imp(op + on, E_r)):>8.3f}"
              f"{A(st['pg_own']):>9.2f}{A(st['pp_own'] + st['pn_own']):>9.2f}")

    # ---- the relay belief (fwd/bwd, in each node's OWN frame) — where does RNA accumulate? ----
    print("\n=== RELAY belief (fwd, own frame) implied composition vs TRUTH ===")
    print(f"{'class':<11}{'fwd_fg':>9}{'fwd_fr':>9}{'bwd_fg':>9}{'bwd_fr':>9}{'true_fg':>9}"
          f"{'p_fwd_g':>9}{'p_fwd_r':>9}")
    for c in (0, 1, 2, 3):
        m = ok & (cls == c)
        if not m.sum():
            continue
        w = M[m]

        def A(x):
            return float(np.average(x[m], weights=w))

        print(f"{_NAMES[c]:<11}{A(imp(st['fwd_g'], E_g)):>9.3f}{A(imp(st['fwd_p'] + st['fwd_n'], E_r)):>9.3f}"
              f"{A(imp(st['bwd_g'], E_g)):>9.3f}{A(imp(st['bwd_p'] + st['bwd_n'], E_r)):>9.3f}"
              f"{A(tf_g):>9.3f}{A(st['fwd_pg']):>9.2f}{A(st['fwd_pp'] + st['fwd_pn']):>9.2f}")

    # ---- the reframe ratio r, per edge type ----
    li, ri_ = st["left"], st["right"]
    rec = iters[-1]
    rho_lf, rho_rf, rho0 = rec["rho_lf"], rec["rho_rf"], st["rho_node0"]
    vl, vr = li >= 0, ri_ >= 0
    r_l = np.where(vl & (rho_rf[np.clip(li, 0, M.size - 1)] > _EPS),
                   rho_lf / np.maximum(rho_rf[np.clip(li, 0, M.size - 1)], _EPS), np.nan)
    r_r = np.where(vr & (rho_lf[np.clip(ri_, 0, M.size - 1)] > _EPS),
                   rho_rf / np.maximum(rho_lf[np.clip(ri_, 0, M.size - 1)], _EPS), np.nan)
    print("\n=== reframe ratio r = rho_tot(dst,face)/rho_tot(src,face), by DST class ===")
    print(f"{'class':<11}{'median r_L':>12}{'median r_R':>12}{'geomean r_L':>13}{'p10':>9}{'p90':>9}")
    for c in (0, 1, 2, 3):
        m = ok & (cls == c)
        a = r_l[m & np.isfinite(r_l)]
        b = r_r[m & np.isfinite(r_r)]
        if a.size == 0:
            continue
        print(f"{_NAMES[c]:<11}{np.median(a):>12.3f}{(np.median(b) if b.size else np.nan):>12.3f}"
              f"{np.exp(np.mean(np.log(np.maximum(a, 1e-9)))):>13.3f}"
              f"{np.percentile(a, 10):>9.3f}{np.percentile(a, 90):>9.3f}")
    print(f"\n# rho_tot(node0) vs M/E_g : median ratio "
          f"{np.median((rho0[ok] / np.maximum(M[ok] / E_g[ok], _EPS))):.3f}")
    print(f"# E_g/E_r median = {np.median(E_g[ok] / np.maximum(E_r[ok], _EPS)):.3f}")

    if "--npz" in sys.argv:
        out = sys.argv[sys.argv.index("--npz") + 1]
        np.savez(out, cls=cls, M=M, E_g=E_g, E_r=E_r, tf_g=tf_g, tf_p=tf_p, tf_n=tf_n,
                 ok=ok, **{k: v for k, v in st.items() if isinstance(v, np.ndarray)},
                 **{f"it{len(iters)}_{k}": v for k, v in iters[-1].items()})
        print(f"\n[wrote {out}]")


if __name__ == "__main__":
    main()
