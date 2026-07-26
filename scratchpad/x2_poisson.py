"""P1d — THE CENTRAL MEASUREMENT: decompose Var(log phi) into EVERY Poisson term plus the residual.

x1_graft.py measured the graft share

    phi = [ rho_nu(bnd) + rho_mu(bnd) ] / [ rho_R(exon) * step ]

and variance_ledger.md §5.2 subtracted only ONE Poisson term (1/n_spliced) from it. But phi is built from
FIVE independently-counted objects:

    numerator    Ru(bnd)       the boundary's TRUE unspliced RNA count (oracle pool)   -> rho_nu
                 S(bnd,face)   the boundary's per-face spliced count (SP+SN)           -> rho_mu
    denominator  (Ru+Rs)(exon) the exon's TRUE RNA count (oracle pool)                  -> rho_R
    step         G(bnd), G(exon)  the oracle gDNA counts forming rho_g(b)/rho_g(i)

The exact Poisson log-variance by the delta method, with the EXACT trigamma psi'(n) (NOT 1/n — at n=1 the
two differ by 64%, and the low-count bin is exactly where the open question lives):

    v_num   = w_nu^2 * psi'(n_unspl_R_bnd) + w_mu^2 * psi'(n_spl)     (M2's share rule: the numerator is a SUM)
    v_den   = psi'(n_R_exon)
    v_step  = psi'(n_g_bnd) + psi'(n_g_exon)
    v_poisson_total = v_num + v_den + v_step

    RESIDUAL = Var(log phi) - E[v_poisson_total]        <- the PREMISE error, if anything survives

    OMP_NUM_THREADS=1 python scratchpad/x2_poisson.py
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np
from scipy.special import polygamma

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
_TINY = 1e-12

CONDS = [
    ("gdna_gdna300_ss_0.99_nrna_present_capture_off", True),
    ("gdna_gdna300_ss_0.50_nrna_present_capture_off", True),
    ("gdna_gdna100_ss_0.50_nrna_present_capture_off", True),
    ("gdna_gdna300_ss_0.50_nrna_none_capture_off", True),
    ("gdna_gdna300_ss_0.99_nrna_present_capture_on", False),
    ("gdna_gdna100_ss_0.50_nrna_present_capture_verystrong", False),
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def psi1(n):
    """EXACT trigamma. psi'(n) ~ 1/n^2 as n->0+, so guard only at n<=0."""
    return polygamma(1, np.maximum(np.asarray(n, float), _TINY))


def psi1_floor1(n):
    """Sensitivity variant: the project's own structural ONE-FRAGMENT floor (bp_solver `_kk = max(.,1)`),
    which caps a count's log-variance at psi'(1) = pi^2/6 = 1.6449."""
    return polygamma(1, np.maximum(np.asarray(n, float), 1.0))


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# PER-EDGE HARVEST
# ────────────────────────────────────────────────────────────────────────────────────────────────────
FIELDS = [
    "face", "n_spl", "n_unspl_R_bnd", "n_R_exon", "n_g_bnd", "n_g_exon",
    "E_r_bnd", "E_r_exon", "E_g_bnd", "E_g_exon", "eff_spl",
    "len_exon", "len_bnd_win", "rho_nu", "rho_mu", "rho_R_ex", "step", "phi",
    "n_spl_pos", "n_spl_neg", "n_Rs_exon", "n_Ru_exon",
]


def harvest(cond: str) -> dict:
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, geo = dbg["chain"], dbg["capture"], dbg["geometry"]
    us = cap["_uni_static"]
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION

    def pool(k):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    G = pool("gdna_pos") + pool("gdna_neg")
    Ru = pool("mat_uns_pos") + pool("nas_uns_pos") + pool("mat_uns_neg") + pool("nas_uns_neg")
    Rs = pool("mat_spl") + pool("nas_spl")
    E_g, E_r = np.asarray(us["E_g"], float), np.asarray(us["E_r"], float)
    li, ri = us["left"], us["right"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    SPf, SNf = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    _ = _node_region_type(chain, ra)

    # genomic bp length per node (regions only; boundaries have no extent -> nan)
    rsz = np.asarray(ra.region_size_bp, float)
    len_node = np.where(isr, rsz[np.clip(idx, 0, rsz.shape[0] - 1)], np.nan)

    rho_g = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)
    rho_R_ex = np.where(E_r > _EPS, (Ru + Rs) / np.maximum(E_r, _EPS), np.nan)
    rho_nu_b = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)

    rec = {f: [] for f in FIELDS}
    for face, nbr in ((1, li), (0, ri)):  # exon i, its neighbour boundary on that side
        for i in np.flatnonzero(is_exon):
            b = nbr[i]
            if b < 0 or not is_bnd[b]:
                continue
            sp, sn = float(SPf[face][b]), float(SNf[face][b])
            espl = float(ESP[face][b])
            mu = (sp + sn) / max(espl, _EPS)
            if not (mu > _EPS) or not np.isfinite(rho_R_ex[i]) or rho_R_ex[i] <= _EPS:
                continue
            if not (rho_g[b] > _EPS and rho_g[i] > _EPS):
                continue
            step = rho_g[b] / rho_g[i]
            phi = (rho_nu_b[b] + mu) / (rho_R_ex[i] * step)
            if not (np.isfinite(phi) and phi > 0):
                continue
            rec["face"].append(face)
            rec["n_spl"].append(sp + sn)
            rec["n_spl_pos"].append(sp)
            rec["n_spl_neg"].append(sn)
            rec["n_unspl_R_bnd"].append(float(Ru[b]))
            rec["n_R_exon"].append(float(Ru[i] + Rs[i]))
            rec["n_Ru_exon"].append(float(Ru[i]))
            rec["n_Rs_exon"].append(float(Rs[i]))
            rec["n_g_bnd"].append(float(G[b]))
            rec["n_g_exon"].append(float(G[i]))
            rec["E_r_bnd"].append(float(E_r[b]))
            rec["E_r_exon"].append(float(E_r[i]))
            rec["E_g_bnd"].append(float(E_g[b]))
            rec["E_g_exon"].append(float(E_g[i]))
            rec["eff_spl"].append(espl)
            rec["len_exon"].append(float(len_node[i]))
            rec["len_bnd_win"].append(float(E_r[b]))
            rec["rho_nu"].append(float(rho_nu_b[b]))
            rec["rho_mu"].append(float(mu))
            rec["rho_R_ex"].append(float(rho_R_ex[i]))
            rec["step"].append(float(step))
            rec["phi"].append(float(phi))
    out = {k: np.asarray(v, float) for k, v in rec.items()}
    # ── the exact Poisson log-variance decomposition ────────────────────────────────────────────────
    tot = out["rho_nu"] + out["rho_mu"]
    w_nu = np.where(tot > _EPS, out["rho_nu"] / np.maximum(tot, _EPS), 0.0)
    w_mu = 1.0 - w_nu
    out["w_nu"], out["w_mu"] = w_nu, w_mu
    for tag, f in (("", psi1), ("_f1", psi1_floor1)):
        v_num = np.where(w_nu > 0, w_nu**2 * f(out["n_unspl_R_bnd"]), 0.0) + w_mu**2 * f(out["n_spl"])
        v_den = f(out["n_R_exon"])
        v_step = f(out["n_g_bnd"]) + f(out["n_g_exon"])
        out["v_num" + tag], out["v_den" + tag], out["v_step" + tag] = v_num, v_den, v_step
        out["v_tot" + tag] = v_num + v_den + v_step
    # ── the EXACT variant ───────────────────────────────────────────────────────────────────────────
    # A SUM's variance must be combined in LINEAR space. For a Gamma(n,1) posterior on a count,
    # Var(N) = n, so Var(rho_i) = n_i/E_i^2 and
    #     Var(log(rho_nu+rho_mu)) = [Var(rho_nu)+Var(rho_mu)] / (rho_nu+rho_mu)^2 = w_nu^2/n_nu + w_mu^2/n_mu
    # i.e. for a SUM the per-component factor is 1/n, NOT psi'(n) — psi' is the LOG-space variance and is
    # correct only for a quantity entering log phi on its OWN (the denominator, and each leg of the step).
    # This matters exactly where the open question lives: w^2*psi'(n) -> w^2/n^2 as n->0, a finite spurious
    # floor, while the exact w^2/n -> 0 (w vanishes linearly in n).
    ve_num = np.where(w_nu > 0, w_nu**2 / np.maximum(out["n_unspl_R_bnd"], _TINY), 0.0) + \
        w_mu**2 / np.maximum(out["n_spl"], _TINY)
    out["v_num_ex"] = ve_num
    out["v_den_ex"] = psi1(out["n_R_exon"])
    out["v_step_ex"] = psi1(out["n_g_bnd"]) + psi1(out["n_g_exon"])
    out["v_tot_ex"] = out["v_num_ex"] + out["v_den_ex"] + out["v_step_ex"]
    out["log_phi"] = np.log(out["phi"])
    return out


DATA = {}
for cond, _off in CONDS:
    DATA[cond] = harvest(cond)
    print(f"harvested {cond:<52} edges={DATA[cond]['phi'].size}", flush=True)
print()


def cat(conds, keys=None):
    keys = keys or list(DATA[conds[0]].keys())
    return {k: np.concatenate([DATA[c][k] for c in conds]) for k in keys}


# ────────────────────────────────────────────────────────────────────────────────────────────────────
def hdr(title):
    print("\n" + "=" * 128)
    print(title)
    print("=" * 128)


# ── TABLE 1 ─────────────────────────────────────────────────────────────────────────────────────────
hdr("(1) PER-CONDITION — Var(log phi) vs the FULL Poisson budget   [trigamma psi'(n), no floor]")
print(f"{'condition':<50}{'n':>6}{'E[log phi]':>11}{'Var(logphi)':>12}"
      f"{'E[v_num]':>10}{'E[v_den]':>10}{'E[v_step]':>11}{'E[v_pois]':>11}{'RESIDUAL':>11}{'%expl':>8}")
print("-" * 128)
for cond, _off in CONDS:
    d = DATA[cond]
    lp = d["log_phi"]
    V = float(np.var(lp))
    vp = float(np.mean(d["v_tot"]))
    print(f"{cond[5:]:<50}{lp.size:>6}{np.mean(lp):>11.3f}{V:>12.4f}"
          f"{np.mean(d['v_num']):>10.4f}{np.mean(d['v_den']):>10.4f}{np.mean(d['v_step']):>11.4f}"
          f"{vp:>11.4f}{V - vp:>11.4f}{100 * vp / max(V, _EPS):>7.1f}%")
for lbl, tg in (("the ONE-FRAGMENT floor psi'(max(n,1)) (project precedent, bp_solver `_kk`)", "_f1"),
                ("*** THE EXACT FORM *** (sum combined in LINEAR space: w^2/n; psi' on den+step)", "_ex")):
    print(f"\n  same, with {lbl}:")
    print(f"{'condition':<50}{'n':>6}{'':>11}{'Var(logphi)':>12}"
          f"{'E[v_num]':>10}{'E[v_den]':>10}{'E[v_step]':>11}{'E[v_pois]':>11}{'RESIDUAL':>11}{'%expl':>8}")
    print("-" * 128)
    for cond, _off in CONDS:
        d = DATA[cond]
        V = float(np.var(d["log_phi"]))
        vp = float(np.mean(d["v_tot" + tg]))
        print(f"{cond[5:]:<50}{d['phi'].size:>6}{'':>11}{V:>12.4f}"
              f"{np.mean(d['v_num' + tg]):>10.4f}{np.mean(d['v_den' + tg]):>10.4f}"
              f"{np.mean(d['v_step' + tg]):>11.4f}{vp:>11.4f}{V - vp:>11.4f}"
              f"{100 * vp / max(V, _EPS):>7.1f}%")

OFF = [c for c, o in CONDS if o]
P = cat(OFF)


def shape_table(d, key, edges, label, tag="_ex", ledger=False):
    print(f"\n{label}")
    extra = f"{'1/n_spl':>10}{'ledger exc':>11}" if ledger else ""
    print(f"{'bin':>14}{'n':>6}{'med '+key:>12}{'E[logphi]':>11}{'Var(logphi)':>12}"
          f"{'E[v_num]':>10}{'E[v_den]':>10}{'E[v_step]':>11}{'E[v_pois]':>11}{'RESIDUAL':>11}{'%expl':>8}"
          + extra)
    print("-" * (116 + len(extra)))
    x = d[key]
    lo = [-np.inf] + list(edges)
    hi = list(edges) + [np.inf]
    for a, b in zip(lo, hi):
        m = (x >= a) & (x < b)
        if m.sum() < 3:
            continue
        lp = d["log_phi"][m]
        V = float(np.var(lp))
        vp = float(np.mean(d["v_tot" + tag][m]))
        nm = f"<{b:g}" if a == -np.inf else (f">{a:g}" if b == np.inf else f"{a:g}-{b:g}")
        ex = ""
        if ledger:
            inv = float(np.mean(1.0 / np.maximum(d["n_spl"][m], _TINY)))
            ex = f"{inv:>10.4f}{V - inv:>11.4f}"
        print(f"{nm:>14}{m.sum():>6}{np.median(x[m]):>12.1f}{np.mean(lp):>11.3f}{V:>12.4f}"
              f"{np.mean(d['v_num'+tag][m]):>10.4f}{np.mean(d['v_den'+tag][m]):>10.4f}"
              f"{np.mean(d['v_step'+tag][m]):>11.4f}{vp:>11.4f}{V - vp:>11.4f}"
              f"{100 * vp / max(V, _EPS):>7.1f}%" + ex)


# ── TABLE 2 ─────────────────────────────────────────────────────────────────────────────────────────
hdr("(2) THE SHAPE TABLE REDONE — binned by JUNCTION SPLICED COUNT, pooled capture-OFF "
    f"({len(OFF)} conditions, n={P['phi'].size})")
shape_table(P, "n_spl", [30, 100, 300, 1000],
            "  *** THE EXACT FORM *** (last two cols reproduce ledger 5.2's Var vs 1/n and its 'excess')",
            tag="_ex", ledger=True)
shape_table(P, "n_spl", [30, 100, 300, 1000], "  [delta-method w^2*psi'(n) on the sum, no floor]", tag="")
shape_table(P, "n_spl", [30, 100, 300, 1000], "  [delta-method w^2*psi'(n), one-fragment floor]",
            tag="_f1")

# ── TABLE 3 ─────────────────────────────────────────────────────────────────────────────────────────
hdr("(3) THE SAME SHAPE, RE-BINNED — does the structure live in the junction count or elsewhere?")
shape_table(P, "n_R_exon", [30, 100, 300, 1000], "  binned by the EXON's RNA count n_R(exon)")
shape_table(P, "n_unspl_R_bnd", [30, 100, 300, 1000],
            "  binned by the BOUNDARY's unspliced RNA count n_unspl_R(bnd)")
shape_table(P, "n_g_bnd", [30, 100, 300, 1000], "  binned by the BOUNDARY's gDNA count n_g(bnd)")
P2 = dict(P)
P2["spl_share"] = P["rho_mu"] / np.maximum(P["rho_mu"] + P["rho_nu"], _EPS)
P2["extrap"] = P["len_exon"] / np.maximum(P["eff_spl"], _EPS)
shape_table(P2, "spl_share", [0.2, 0.4, 0.6, 0.8],
            "  binned by the SPLICED SHARE of the numerator w_mu (the M2 weight)")
shape_table(P2, "extrap", [4, 8, 16, 32],
            "  binned by the EXTRAPOLATION RATIO len(exon)/eff_spl (the '12-21x' of ledger §3)")
shape_table(P2, "len_exon", [400, 700, 1200], "  binned by the EXON's genomic LENGTH (bp)")
shape_table(P2, "rho_mu", [0.3, 1.5, 5.0], "  binned by the junction's spliced DENSITY rho_mu")
shape_table(P2, "rho_R_ex", [0.4, 2.0, 7.5], "  binned by the exon's RNA DENSITY rho_R(exon)")

# ── TABLE 4 ─────────────────────────────────────────────────────────────────────────────────────────
hdr("(4) MEDIAN-BASED ROBUST VIEW of table (2) — matched statistics only "
    "(ledger §5.1: a median vs an sd overturned a conclusion once already)")
print("  robust sd  = 1.4826 * MAD(log phi)  (the exact normal-consistency constant)")
print("  med|logphi| is comparable to med sqrt(v_pois) ONLY as a scale, not as a variance.\n")
def robust_table(d, key, edges, label):
    print(f"\n{label}")
    print(f"{'bin':>14}{'n':>6}{'med|logphi|':>13}{'med sqrt(vp)':>14}{'ratio':>8}"
          f"{'robust sd':>11}{'robust var':>12}{'med v_pois':>12}{'rob RESID':>11}"
          f"{'mean RESID':>12}{'rob/mean':>10}")
    print("-" * 126)
    x = d[key]
    lo = [-np.inf] + list(edges)
    hi = list(edges) + [np.inf]
    for a, b in zip(lo, hi):
        m = (x >= a) & (x < b)
        if m.sum() < 3:
            continue
        lp = d["log_phi"][m]
        vp = d["v_tot_ex"][m]
        med_abs = float(np.median(np.abs(lp)))
        med_sq = float(np.median(np.sqrt(vp)))
        rsd = 1.4826 * float(np.median(np.abs(lp - np.median(lp))))
        rres = rsd**2 - float(np.median(vp))
        mres = float(np.var(lp) - np.mean(vp))
        nm = f"<{b:g}" if a == -np.inf else (f">{a:g}" if b == np.inf else f"{a:g}-{b:g}")
        print(f"{nm:>14}{m.sum():>6}{med_abs:>13.4f}{med_sq:>14.4f}"
              f"{med_abs / max(med_sq, _EPS):>8.2f}{rsd:>11.4f}{rsd**2:>12.4f}"
              f"{np.median(vp):>12.4f}{rres:>11.4f}{mres:>12.4f}{rres / max(mres, _EPS):>10.2f}")


robust_table(P, "n_spl", [30, 100, 300, 1000], "  pooled capture-OFF, binned by junction spliced count")

# ── EXTRAS ──────────────────────────────────────────────────────────────────────────────────────────
hdr("(5) THE LOW-COUNT REGIME (where trigamma and 1/n diverge) and the WELL-COUNTED VIEW")
five = ["n_spl", "n_unspl_R_bnd", "n_R_exon", "n_g_bnd", "n_g_exon"]
print(f"{'condition':<50}" + "".join(f"{'frac ' + k + '<10':>22}" for k in ("n_g_bnd", "n_unspl_R"))
      + f"{'frac ANY<10':>13}{'frac ANY<1':>12}{'frac all>=30':>14}")
print("-" * 128)
for cond, _off in CONDS:
    d = DATA[cond]
    fg = float(np.mean(d["n_g_bnd"] < 10))
    fu = float(np.mean(d["n_unspl_R_bnd"] < 10))
    anyl = np.zeros(d["phi"].size, bool)
    any1 = np.zeros(d["phi"].size, bool)
    ok30 = np.ones(d["phi"].size, bool)
    for k in five:
        anyl |= d[k] < 10
        any1 |= d[k] < 1
        ok30 &= d[k] >= 30
    print(f"{cond[5:]:<50}{fg:>22.3f}{fu:>22.3f}{np.mean(anyl):>13.3f}"
          f"{np.mean(any1):>12.3f}{np.mean(ok30):>14.3f}")

print("\n  WELL-COUNTED ONLY (all five counts >= 30) — the cleanest view of the PREMISE error:")
print(f"{'condition':<50}{'n':>6}{'E[logphi]':>11}{'Var(logphi)':>12}"
      f"{'E[v_num]':>10}{'E[v_den]':>10}{'E[v_step]':>11}{'E[v_pois]':>11}{'RESIDUAL':>11}{'%expl':>8}")
print("-" * 128)
for cond, _off in CONDS:
    d = DATA[cond]
    ok = np.ones(d["phi"].size, bool)
    for k in five:
        ok &= d[k] >= 30
    if ok.sum() < 5:
        print(f"{cond[5:]:<50}{ok.sum():>6}   (too few)")
        continue
    lp = d["log_phi"][ok]
    V = float(np.var(lp))
    vp = float(np.mean(d["v_tot_ex"][ok]))
    print(f"{cond[5:]:<50}{ok.sum():>6}{np.mean(lp):>11.3f}{V:>12.4f}"
          f"{np.mean(d['v_num_ex'][ok]):>10.4f}{np.mean(d['v_den_ex'][ok]):>10.4f}"
          f"{np.mean(d['v_step_ex'][ok]):>11.4f}{vp:>11.4f}{V - vp:>11.4f}"
          f"{100 * vp / max(V, _EPS):>7.1f}%")

okP = np.ones(P["phi"].size, bool)
for k in five:
    okP &= P[k] >= 30
Pw = {k: v[okP] for k, v in P.items()}
print(f"\n  pooled capture-OFF, well-counted: n={okP.sum()}  E[log phi]={np.mean(Pw['log_phi']):.3f}  "
      f"Var={np.var(Pw['log_phi']):.4f}  E[v_pois]={np.mean(Pw['v_tot_ex']):.4f}  "
      f"RESID={np.var(Pw['log_phi']) - np.mean(Pw['v_tot_ex']):.4f}")
shape_table(Pw, "n_spl", [100, 300, 1000], "  well-counted, binned by junction spliced count")
robust_table(Pw, "n_spl", [100, 300, 1000], "  well-counted, MEDIAN view of the same bins")
print("\n  per-condition ROBUST vs MEAN, all edges (do the two views agree?):")
print(f"{'condition':<50}{'mean Var':>10}{'robust var':>12}{'E[v_pois]':>11}{'med v_pois':>12}"
      f"{'mean RESID':>12}{'rob RESID':>11}")
print("-" * 128)
for cond, _off in CONDS:
    d = DATA[cond]
    lp, vp = d["log_phi"], d["v_tot_ex"]
    rsd = 1.4826 * float(np.median(np.abs(lp - np.median(lp))))
    print(f"{cond[5:]:<50}{np.var(lp):>10.4f}{rsd**2:>12.4f}{np.mean(vp):>11.4f}"
          f"{np.median(vp):>12.4f}{np.var(lp) - np.mean(vp):>12.4f}{rsd**2 - np.median(vp):>11.4f}")

# ── the shared-fragment correlation caveat, quantified ───────────────────────────────────────────────
hdr("(6) INDEPENDENCE AUDIT — the law-of-total-variance subtraction assumes the per-edge Poisson error "
    "is independent of the premise error")
print("  a) DO the five counts share fragments? The exon's oracle SPLICED pool n_Rs(exon) is the test:")
print("     a spliced fragment crosses a junction, so the accumulator books it as BOUNDARY mass, never as")
print("     contained-region mass. If n_Rs(exon) == 0 identically, all five counts are drawn from DISJOINT")
print("     fragment sets and their Poisson errors are mutually independent GIVEN the true densities.")
print(f"{'condition':<50}{'med n_spl/n_R_exon':>21}{'mean n_Rs/n_R exon':>20}{'corr(logn_spl,logn_Rex)':>26}")
print("-" * 128)
for cond, _off in CONDS:
    d = DATA[cond]
    r1 = np.median(d["n_spl"] / np.maximum(d["n_R_exon"], _EPS))
    r2 = np.mean(d["n_Rs_exon"] / np.maximum(d["n_R_exon"], _EPS))
    cc = np.corrcoef(np.log(np.maximum(d["n_spl"], _TINY)),
                     np.log(np.maximum(d["n_R_exon"], _TINY)))[0, 1]
    print(f"{cond[5:]:<50}{r1:>21.3f}{r2:>20.3f}{cc:>26.3f}")
print("\n  b) is the per-edge Poisson budget CORRELATED with |log phi| (it must not be, for the")
print("     law-of-total-variance subtraction to be a clean decomposition)?")
print(f"{'condition':<50}{'corr(v_pois,|logphi|)':>24}{'corr(log v_pois, log|logphi|)':>32}")
print("-" * 128)
for cond, _off in CONDS:
    d = DATA[cond]
    a, b = d["v_tot_ex"], np.abs(d["log_phi"])
    ok = b > 0
    print(f"{cond[5:]:<50}{np.corrcoef(a, b)[0, 1]:>24.3f}"
          f"{np.corrcoef(np.log(a[ok]), np.log(b[ok]))[0, 1]:>32.3f}")

hdr("(7) TAIL CONCENTRATION + BOOTSTRAP SE — how trustworthy is the mean-based RESIDUAL?")
rng = np.random.default_rng(0)
print(f"{'condition':<50}{'n':>6}{'top1% of Var':>14}{'top5% of Var':>14}"
      f"{'RESIDUAL':>11}{'boot SE':>10}{'boot 95% CI':>22}")
print("-" * 128)
for cond, _off in CONDS:
    d = DATA[cond]
    lp, vp = d["log_phi"], d["v_tot_ex"]
    dev = (lp - lp.mean()) ** 2
    o = np.sort(dev)[::-1]
    n1, n5 = max(1, int(0.01 * o.size)), max(1, int(0.05 * o.size))
    bs = np.array([
        (lambda s: float(np.var(lp[s]) - np.mean(vp[s])))(rng.integers(0, lp.size, lp.size))
        for _ in range(2000)
    ])
    print(f"{cond[5:]:<50}{lp.size:>6}{o[:n1].sum() / o.sum():>14.3f}"
          f"{o[:n5].sum() / o.sum():>14.3f}{np.var(lp) - np.mean(vp):>11.4f}{bs.std():>10.4f}"
          f"{f'[{np.percentile(bs, 2.5):.3f}, {np.percentile(bs, 97.5):.3f}]':>22}")
print("\n  pooled capture-OFF (n=%d): RESIDUAL=%.4f" % (P["phi"].size,
                                                        np.var(P["log_phi"]) - np.mean(P["v_tot_ex"])))
_lp, _vp = P["log_phi"], P["v_tot_ex"]
_bs = np.array([(lambda s: float(np.var(_lp[s]) - np.mean(_vp[s])))(rng.integers(0, _lp.size, _lp.size))
                for _ in range(2000)])
print(f"    bootstrap SE = {_bs.std():.4f}   95% CI = [{np.percentile(_bs, 2.5):.3f}, "
      f"{np.percentile(_bs, 97.5):.3f}]")
print("  the four capture-OFF per-condition residuals: "
      + ", ".join(f"{np.var(DATA[c]['log_phi']) - np.mean(DATA[c]['v_tot_ex']):.4f}" for c in OFF))
_r = np.array([np.var(DATA[c]["log_phi"]) - np.mean(DATA[c]["v_tot_ex"]) for c in OFF])
print(f"    -> mean {_r.mean():.4f}  sd {_r.std(ddof=1):.4f}  ({100 * _r.std(ddof=1) / _r.mean():.1f}% "
      f"spread across gDNA depth x strand-specificity x nascent present/absent)")

np.savez("/Users/mkiyer/proj/rigel/scratchpad/x2_poisson_edges.npz",
         **{f"{c}__{k}": v for c, _o in CONDS for k, v in DATA[c].items()})
print("\nper-edge records -> scratchpad/x2_poisson_edges.npz")
