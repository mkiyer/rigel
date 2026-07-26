"""BOUNDARY STUDY step 1 — what a boundary OBSERVES, what it must deconvolve, and what reaches it.

OWNER'S MODEL. An exon region cannot know more than its gDNA-vs-RNA split (RNA is ONE species; only its
direction is unknown). The BOUNDARY is the crux: it directly measures how much RNA is SPLICED (leaves for
another exon) and how much is UNSPLICED (continues). That split is observed. The hard part is deconvolving the
UNSPLICED crossing mass  U = gDNA + RNA_unspliced,  and three things can do it:
   (1) strand-specific data (the tilt),
   (2) the intron gDNA factory next door (gDNA density),
   (3) the exon's own gDNA estimate, via the exon -> boundary message.

This measures, per boundary class (flanking pair x junction present/absent):
   - the ORACLE composition of the crossing mass: gDNA / unspliced-RNA, and the spliced RNA beside it;
   - how well f_g on the unspliced mass is solved (self vs solved);
   - which of the three sources actually arrives (strand tau, factory-intron neighbour, exon gDNA message);
   - the ENRICHMENT ordering exon vs boundary vs intron (oracle gDNA densities = capture);
   - COVERAGE of the regime the owner flags as untested: unspliced RNA density > spliced RNA density.

    OMP_NUM_THREADS=1 python scratchpad/b1_boundary.py
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
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def mw(x, w):
    return np.average(x, weights=w) if np.sum(w) > 0 else np.nan


for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st, geo = dbg["chain"], dbg["capture"], dbg["statics"], dbg["geometry"]
    us = cap["_uni_static"]
    uni = cap["_uni"][-1]
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION

    def pool(k):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    G = pool("gdna_pos") + pool("gdna_neg")
    R_u = pool("mat_uns_pos") + pool("nas_uns_pos") + pool("mat_uns_neg") + pool("nas_uns_neg")
    R_s = pool("mat_spl") + pool("nas_spl")
    U = G + R_u
    fo = np.where(U > _EPS, G / np.maximum(U, _EPS), np.nan)
    mass = np.asarray(cap["mass_global"])
    ship, self_fg = np.asarray(cap["f_g"]), np.asarray(cap["fg_loc"])
    tau0 = np.asarray(cap["_tau0_lam"])
    E_g = us["E_g"]
    li, ri = us["left"], us["right"]
    n = len(mass)
    rt, _ = _node_region_type(chain, ra)
    CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array([CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(n)])
    ms = np.asarray(st.mass_spliced, float)
    is_bnd = ~isr
    pair = np.array(["-"] * n, dtype=object)
    for i in np.flatnonzero(is_bnd):
        a = cls[li[i]] if li[i] >= 0 else "edge"
        b = cls[ri[i]] if ri[i] >= 0 else "edge"
        pair[i] = " | ".join(sorted([str(a), str(b)]))
    ok = is_bnd & np.isfinite(fo) & (mass > _EPS) & np.asarray(cap["solvable"], bool)
    err = np.abs(ship - fo) * mass
    E = err[ok].sum()
    # which sources reach this boundary?
    nb_fac = np.zeros(n, bool)
    for i in np.flatnonzero(is_bnd):
        for k in (li[i], ri[i]):
            if k >= 0 and tau0[k] > _EPS:
                nb_fac[i] = True

    print(f"\n{'=' * 126}\n{cond[5:]}   boundary ERR {E:,.0f} reads\n{'=' * 126}")
    print(f"  {'flank pair':<26}{'junc':>5}{'n':>5}{'reads':>10}{'ERR':>9}{'share':>7}"
          f"|{'orc f_g':>8}{'self':>7}{'ship':>7}|{'R_s/U':>7}{'R_u/R_s':>8}"
          f"|{'tau>0':>6}{'facNbr':>7}{'cm_g':>7}{'c_tau':>7}")
    rows = []
    for p in sorted(set(pair[ok])):
        for jl, jm in (("yes", ms > _EPS), ("no", ms <= _EPS)):
            m = ok & (pair == p) & jm
            if m.sum() < 3:
                continue
            rows.append((err[m].sum(), p, jl, m))
    for e, p, jl, m in sorted(rows, reverse=True)[:10]:
        w = mass[m]
        print(f"  {p:<26}{jl:>5}{int(m.sum()):>5}{w.sum():>10,.0f}{e:>9,.0f}{e / E:>7.1%}"
              f"|{mw(fo[m], w):>8.4f}{mw(np.abs(self_fg - fo)[m], w):>7.4f}{mw(np.abs(ship - fo)[m], w):>7.4f}"
              f"|{mw((R_s / np.maximum(U, _EPS))[m], w):>7.2f}"
              f"{mw((R_u / np.maximum(R_s, _EPS))[m], w):>8.3f}"
              f"|{np.mean(tau0[m] > _EPS):>6.0%}{np.mean(nb_fac[m]):>7.0%}"
              f"{mw(uni['cm_g'][m], w):>7.2f}{mw(uni['c_tau'][m], w):>7.2f}")

    # ── the ENRICHMENT ordering (oracle gDNA density = capture, since gDNA content is uniform) ──
    rho_g = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)
    ex_b, in_b = [], []
    for i in np.flatnonzero(is_bnd & np.isfinite(rho_g)):
        for k in (li[i], ri[i]):
            if k >= 0 and np.isfinite(rho_g[k]) and rho_g[k] > 0 and rho_g[i] > 0:
                (ex_b if cls[k] == "exon" else in_b if cls[k] == "intron" else []).append(rho_g[i] / rho_g[k])
    ex_b, in_b = np.asarray(ex_b), np.asarray(in_b)
    if ex_b.size and in_b.size:
        print(f"\n  ENRICHMENT (oracle gDNA density ratios): boundary/EXON median {np.median(ex_b):.3f} "
              f"(bnd>exon on {np.mean(ex_b > 1):.0%})   boundary/INTRON median {np.median(in_b):.3f} "
              f"(bnd<intron on {np.mean(in_b < 1):.0%})")
    # ── COVERAGE of the untested regime: unspliced RNA density > spliced RNA density at a boundary ──
    ESP = np.asarray(geo.eff_spl_left, float) + np.asarray(geo.eff_spl_right, float)
    E_r = us["E_r"]
    d_u = np.where(E_r > _EPS, R_u / np.maximum(E_r, _EPS), np.nan)
    d_s = np.where(ESP > _EPS, R_s / np.maximum(ESP, _EPS), np.nan)
    jm = ok & (ms > _EPS) & np.isfinite(d_u) & np.isfinite(d_s) & (d_s > 0)
    if jm.sum() > 5:
        rr = (d_u / d_s)[jm]
        print(f"  COVERAGE  unspliced/spliced RNA DENSITY at junction boundaries: median {np.median(rr):.4f}, "
              f"p90 {np.quantile(rr, 0.9):.4f}, max {rr.max():.3f}  —  fraction > 1: "
              f"{np.mean(rr > 1):.1%}  (the owner's untested regime)")
