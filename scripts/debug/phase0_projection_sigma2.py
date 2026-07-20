"""PHASE 0 — does the belief-free NPMLE projection σ²_transfer reproduce the oracle enrichment-crossing strata?

The Phase-0 gate (docs/calibration/npmle_projection_variance_design.md §6): fit P_total belief-free, compute
the projection σ²_transfer for every edge, and check it reproduces the ORACLE gDNA transfer-variance strata
(dep-dep small, enr-enr moderate, crossing large) WITHOUT ever touching a solved belief. Settles the two form
decisions (§5.1 h, §5.2 the σ² formula) against the oracle, not by tuning.

Per node i, project the belief-free total log-density d_i = log(M_i/E_i) onto P_total = Σ_j w_j N(μ_j, h²):
    r_ij ∝ w_j·exp(−½((d_i−μ_j)/h)²)      μ_proj_i = Σ_j r_ij μ_j     Var_proj_i = Σ_j r_ij(μ_j−μ_proj_i)² + h²
Per edge (src→dst), three candidate σ²_transfer formulas:
    F0 gap-only : (μ_proj_dst − μ_proj_src)² + 2h²
    F1 predictive: Var_proj_dst + (μ_proj_dst − μ_proj_src)²         (the design's default)
    F2 symmetric : Var_proj_src + Var_proj_dst + (μ_proj_dst − μ_proj_src)²

ORACLE = the true gDNA log-rate difference Δ = log ρ_g,dst − log ρ_g,src; per stratum its Poisson-removed
variance (bp_solver._poisson_moment_var) is the truth the projection must reproduce. Strata by the BELIEF-FREE
total-density regime (median split) — the same label the projection sees. Also a calibration curve: bin edges
by predicted σ², plot realized Var(Δ_true) — ideally on y=x.

    OMP_NUM_THREADS=1 python phase0_projection_sigma2.py [--h 0.25] [--out DIR]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

sys.path.insert(0, str(Path("/Users/mkiyer/proj/rigel/scripts/debug")))
from selfsolve_diag import _scan_and_truth  # noqa: E402
from npmle_projection_variance import build, fit_weights  # noqa: E402

from _disagreement_variance import _adjacent_edges, _poisson_moment_var
from rigel.calibration.node_chain import build_node_chain
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12
_LN10 = np.log(10.0)


def project_all(d, log_rho, w, h):
    """Vectorized projection of node log-densities d (n,) onto the mixture -> (mu_proj, var_proj)."""
    z = (d[:, None] - log_rho[None, :]) / h  # (n,G)
    lr = np.log(np.maximum(w, _EPS))[None, :] - 0.5 * z * z
    lr -= lr.max(axis=1, keepdims=True)
    r = np.exp(lr)
    r /= r.sum(axis=1, keepdims=True)
    mu = (r * log_rho[None, :]).sum(1)
    var = (r * (log_rho[None, :] - mu[:, None]) ** 2).sum(1) + h * h
    return mu, var


def measure(inp, index, h):
    M, G, R, eff = build(inp, index)  # live nodes only
    e = np.maximum(eff, _EPS)
    d_tot = np.log(np.maximum(M, 1.0) / e)  # belief-free total log-density (opportunity floor)
    lg = np.log(np.maximum(G, 1.0) / e)  # oracle gDNA log-rate

    # rebuild the edge list on the SAME live-node indexing build() used
    ra_pl = inp["payload"]
    chain = build_node_chain(ra_pl.ref_region_offsets, ra_pl.ref_boundary_offsets)
    src_full, dst_full, _ = _adjacent_edges(chain)
    # build() kept nodes with eff>1e-9*1.001; reconstruct that mask to remap edge indices
    from rigel.calibration.bp_solver import node_global_geometry
    from rigel.calibration.node_geometry import build_node_geometry
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(ra_pl, ra)
    bsub = BoundarySubstrate.from_payload(ra_pl)
    geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])
    _m, eff_all = node_global_geometry(chain, geom)
    live = np.asarray(eff_all, float) > 1e-9 * 1.001
    remap = -np.ones(live.shape[0], np.int64)
    remap[live] = np.arange(int(live.sum()))
    ok = live[src_full] & live[dst_full]
    src, dst = remap[src_full[ok]], remap[dst_full[ok]]

    # fit P_total belief-free
    dens_all = np.concatenate([M / e, G / e, R / e])
    dpos = dens_all[dens_all > 0]
    lo = float(np.log(np.percentile(dpos, 0.2))) - 1.0 * _LN10
    hi = float(np.log(dpos.max())) + 0.4 * _LN10
    log_rho = np.linspace(lo, hi, 320)
    w, _dens = fit_weights(M, e, log_rho, h)
    mu_p, var_p = project_all(d_tot, log_rho, w, h * _LN10)

    gap2 = (mu_p[dst] - mu_p[src]) ** 2
    F0 = gap2 + 2.0 * (h * _LN10) ** 2
    F1 = var_p[dst] + gap2
    F2 = var_p[src] + var_p[dst] + gap2

    # belief-free regime (median split on total density) -> strata
    nid = np.unique(np.concatenate([src, dst]))
    reg = (d_tot > np.median(d_tot[nid])).astype(int)
    s, dd = reg[src], reg[dst]
    strata = {"dep-dep": (s == 0) & (dd == 0), "enr-enr": (s == 1) & (dd == 1), "MIXED": s != dd}

    Dgd = lg[dst] - lg[src]  # oracle gDNA log-rate difference
    rows = []
    for name, m in strata.items():
        if int(m.sum()) < 2:
            continue
        orc_raw = float(np.var(Dgd[m]))  # RAW disagreement variance (includes Poisson — apples-to-apples w/ F*)
        orc_pois = _poisson_moment_var(Dgd[m], np.maximum(G[src[m]], _EPS), np.maximum(G[dst[m]], _EPS))
        rows.append(dict(stratum=name, n=int(m.sum()), orc_raw=orc_raw, orc_pois=orc_pois,
                         F0=float(F0[m].mean()), F1=float(F1[m].mean()), F2=float(F2[m].mean())))
    return pd.DataFrame(rows), (F1, np.abs(Dgd) ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--hs", default="0.15,0.25,0.35")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    outdir = Path(a.out) if a.out else Path(os.environ.get("RIGEL_SCRATCH", "/tmp"))
    outdir.mkdir(parents=True, exist_ok=True)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    conds = sorted(p.stem for p in cache.glob("*.pkl"))
    hs = [float(x) for x in a.hs.split(",")]

    calib_pts = []
    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        cap = "capON" if "capture_on" in c else "capOFF"
        nas = "+nas" if "nrna_present" in c else "    "
        print(f"\n{'=' * 92}\n{c}  [{cap} {nas}]")
        for h in hs:
            df, cp = measure(inp, index, h)
            if h == 0.25:
                calib_pts.append((c, cp))
            print(f"  h={h}:")
            print(df.to_string(index=False, float_format=lambda x: f"{x:,.3g}",
                               columns=["stratum", "n", "orc_raw", "orc_pois", "F0", "F1", "F2"]))

    # calibration curve at h=0.25 (predicted F1 vs realized Var of the true gDNA Δ), pooled
    fig, ax = plt.subplots(figsize=(6.5, 6))
    allF1 = np.concatenate([cp[0] for _, cp in calib_pts])
    allD2 = np.concatenate([cp[1] for _, cp in calib_pts])
    q = np.quantile(allF1, np.linspace(0, 1, 13))
    q = np.unique(q)
    xm, ym = [], []
    for lo, hi in zip(q[:-1], q[1:]):
        m = (allF1 >= lo) & (allF1 < hi)
        if m.sum() >= 20:
            xm.append(float(np.median(allF1[m])))
            ym.append(float(np.mean(allD2[m])))  # realized E[Δ²] (raw; includes Poisson — upper bound)
    ax.loglog(xm, ym, "o-", label="realized E[Δ²_true] per predicted-σ² bin")
    lim = [min(xm + ym), max(xm + ym)]
    ax.loglog(lim, lim, "k--", alpha=0.5, label="y = x (perfect calibration)")
    ax.set_xlabel("predicted σ²_transfer (F1, projection)")
    ax.set_ylabel("realized E[(Δ log ρ_g,true)²]")
    ax.set_title("Phase 0 calibration — projection σ² vs oracle realized variance (h=0.25, pooled)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = outdir / "phase0_projection_calibration.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
