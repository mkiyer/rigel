"""PASS-0 NODE DISSECTION — WHY does a node solve wrong? Channel-ablation replay of the ψ solve.

For the worst nodes of one scenario (ranked by **error mass** ``mass·|f_g−oracle|``), this replays
`simplex_logodds._solve_nodes_logodds_all` with the ψ channels switched on one at a time, using the SAME
inputs the shipped solve used (including ``fg_ref`` — the incoming belief the variance is frozen at, and the
intron-factory λ arm). The replay is validated against the shipped ``f_g`` before anything is ablated, so an
attribution here is a measurement, not a reconstruction.

The channels (`CALIBRATION_ARCHITECTURE.md` — the three legitimate information sources):

    S      strand Beta-Binomial likelihood ALONE          (the only INTRINSIC gDNA/RNA signal)
    S+P    + the population gDNA prior / reference        (global_logprior + the intron factory)
    S+P+G  + the gDNA imputation message
    S+P+R  + the per-strand RNA imputation messages (no gDNA message)
    ALL    everything = the shipped solve

Reading it: whichever step moves ``f_g`` AWAY from the oracle is the channel at fault, and the per-node dump
(oracle mature/nascent/gDNA split, message modes + precisions, geometry, neighbours) says why.

    OMP_NUM_THREADS=1 python scripts/debug/pass0_node_dissect.py [cond] [--top 12] [--cls exon]
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
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-9
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
CLASSES = ("intergenic", "intron", "exon", "boundary")


def _pools(inp, chain, keys):
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION
    a = sum(np.asarray(inp["region_pools"][k], float) for k in keys)
    b = sum(np.asarray(inp["boundary_pools"][k], float) for k in keys)
    return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])


def load(cond, factory: bool | None = None):
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_index(index)
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    # PASS-0 ONLY (no refit). ``factory=None`` ⇒ follow the shipped config default.
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0,
                             **({} if factory is None else {"intron_factory": bool(factory)}))
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    return index, ra, inp, dbg, cc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cond", nargs="?", default="gdna_gdna300_ss_0.50_nrna_present_capture_on")
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--cls", default=None, choices=list(CLASSES))
    ap.add_argument("--factory", default=None, choices=["0", "1"],
                    help="override the gDNA intron factory (default: follow the shipped config)")
    a = ap.parse_args()

    index, ra, inp, dbg, cc = load(a.cond, factory=None if a.factory is None else a.factory == "1")
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    mat = _pools(inp, chain, ("mat_uns_pos", "mat_uns_neg"))
    nas = _pools(inp, chain, ("nas_uns_pos", "nas_uns_neg"))
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    fg = np.asarray(cap["f_g"])
    mass = np.asarray(cap["mass_global"])
    eff = np.asarray(cap["eff_global"])
    rt, _ = _node_region_type(chain, ra)
    cls = np.where(np.asarray(chain.kind) != 0, 3, rt)
    ok = np.isfinite(fo) & (mass > _EPS) & np.asarray(cap["solvable_mask"], bool)

    # ── the channel-ablation replay ────────────────────────────────────────────────────────────────────
    base_kw = dict(
        kappa=float(dbg["rna_sense_frac"]), od_g=0.0, od_r=0.0,
        n_grid=cc.sweep_n_grid, L=cc.sweep_logodds_window,
        n_tilt=cc.sweep_n_tilt, n_grid_ss=cc.sweep_n_grid_single_strand,
        fg_ref=cap["fg_init"], fpos_ref=cap["fpos_init"], fneg_ref=cap["fneg_init"],
    )
    # the strand overdispersions the sweep actually used
    prs = dbg["calibration_priors"]
    base_kw["od_g"] = float(prs.gdna_strand_overdispersion)
    base_kw["od_r"] = float(prs.rna_strand_overdispersion)
    pos_args = (st.u_pos, st.u_neg, st.free_pos, st.free_neg, st.mass_unspliced, st.mass_spliced)
    GLP, IP = cap["global_lp"], cap["intron_prior"]
    MG, PG = cap["mode_g"], cap["prec_g"]
    MP, PP, MN, PN = cap["mode_p"], cap["prec_p"], cap["mode_n"], cap["prec_n"]

    def solve(prior=False, gdna=False, rna=False):
        kw = dict(base_kw)
        if prior:
            kw.update(global_logprior=GLP, lam_logprior=IP)
        if gdna:
            kw.update(gdna_imp_mode=MG, gdna_imp_prec=PG)
        if rna:
            kw.update(rna_imp_mode=(MP, MN), rna_imp_prec=(PP, PN))
        return np.asarray(_solve_nodes_logodds_all(*pos_args, **kw).gdna_frac)

    A = {
        "S": solve(),
        "S+P": solve(prior=True),
        "S+P+G": solve(prior=True, gdna=True),
        "S+P+R": solve(prior=True, rna=True),
        "ALL": solve(prior=True, gdna=True, rna=True),
    }
    fid = float(np.max(np.abs(A["ALL"][ok] - fg[ok])))
    print(f"# cond={a.cond}")
    print(f"# REPLAY FIDELITY  max|replay(ALL) - shipped f_g| = {fid:.2e}"
          f"   {'OK' if fid < 1e-6 else '*** MISMATCH — attribution below is unreliable ***'}")

    # ── aggregate attribution, by class ───────────────────────────────────────────────────────────────
    print("\n=== CHANNEL ATTRIBUTION — mass-weighted mean f_g at each stage (oracle in the last column) ===")
    print(f"{'class':<11}{'n':>6}{'S':>9}{'S+P':>9}{'S+P+G':>9}{'S+P+R':>9}{'ALL':>9}{'oracle':>9}")
    for ci, cn in enumerate(CLASSES):
        m = ok & (cls == ci)
        if not m.sum():
            continue
        w = mass[m]
        row = "".join(f"{np.average(A[k][m], weights=w):>9.3f}" for k in A)
        print(f"{cn:<11}{int(m.sum()):>6}{row}{np.average(fo[m], weights=w):>9.3f}")
    print(f"\n{'class':<11}{'n':>6}{'|S|':>9}{'|S+P|':>9}{'|S+P+G|':>9}{'|S+P+R|':>9}{'|ALL|':>9}  (mwae)")
    for ci, cn in enumerate(CLASSES):
        m = ok & (cls == ci)
        if not m.sum():
            continue
        w = mass[m]
        row = "".join(f"{np.average(np.abs(A[k][m] - fo[m]), weights=w):>9.4f}" for k in A)
        print(f"{cn:<11}{int(m.sum()):>6}{row}")

    # ── the worst nodes ───────────────────────────────────────────────────────────────────────────────
    sel = ok & ((cls == CLASSES.index(a.cls)) if a.cls else True)
    errmass = mass * np.abs(fg - fo)
    order = np.argsort(np.where(sel, -errmass, 0.0))[: a.top]
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    tau0 = np.asarray(cap["_tau0_lam"])
    up, un = np.asarray(st.u_pos, float), np.asarray(st.u_neg, float)
    fp_, fn_ = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    spl = np.asarray(cap["spl_l"]) + np.asarray(cap["spl_r"])

    print(f"\n\n{'='*118}\nWORST {a.top} NODES BY ERROR MASS"
          f"{' (class=' + a.cls + ')' if a.cls else ''}\n{'='*118}")
    for i in order:
        i = int(i)
        if errmass[i] <= 0:
            continue
        cn = CLASSES[int(cls[i])]
        tot = G[i] + R[i]
        print(f"\n── node {i}  [{cn}]  errmass={errmass[i]:,.0f}  mass={mass[i]:,.1f}  eff={eff[i]:,.1f}")
        print(f"   ORACLE   f_g={fo[i]:.3f}   gDNA={G[i]:,.1f}  RNA={R[i]:,.1f}"
              f"  (mature={mat[i]:,.1f}  nascent={nas[i]:,.1f})"
              f"   RNA+/-={Rp[i]:,.0f}/{Rn[i]:,.0f}  gDNA+/-={Gp[i]:,.0f}/{Gn[i]:,.0f}")
        print(f"   SOLVE    S={A['S'][i]:.3f}  S+P={A['S+P'][i]:.3f}  S+P+G={A['S+P+G'][i]:.3f}"
              f"  S+P+R={A['S+P+R'][i]:.3f}  ALL={A['ALL'][i]:.3f}   |err|={abs(fg[i]-fo[i]):.3f}")
        print(f"   COUNTS   u_pos={up[i]:.0f} u_neg={un[i]:.0f}  free=({int(fp_[i])},{int(fn_[i])})"
              f"  spliced={spl[i]:,.1f}  tau0_lam={tau0[i]:.4g}"
              f"  (tau0=0 ⇒ NO intrinsic strand evidence)")
        pri = None
        if GLP is not None:
            g = np.asarray(cap["solve_grid"])
            pri = float(g[int(np.argmax(np.asarray(GLP)[i]))])
        print(f"   MSG      gDNA mode={np.exp(MG[i]):.3f} prec={PG[i]:.3g}"
              f"   RNA+ mode={np.exp(MP[i]):.3f} prec={PP[i]:.3g}"
              f"   RNA- mode={np.exp(MN[i]):.3f} prec={PN[i]:.3g}"
              + (f"   prior argmax f_g={pri:.3f}" if pri is not None else "   prior=None"))
        nb = [("L", int(left[i])), ("R", int(right[i]))]
        for tag, j in nb:
            if j < 0:
                continue
            print(f"   nbr {tag}   node {j} [{CLASSES[int(cls[j])]}]  oracle={fo[j]:.3f}"
                  f"  solved={fg[j]:.3f}  mass={mass[j]:,.1f}  rho={mass[j]/max(eff[j],_EPS):.4f}")
        _ = tot


if __name__ == "__main__":
    main()
