"""A/B — the NPMLE projection σ²_transfer wired into the real solve vs the shipped σ²=0.

Reproduces the production `calibrate` sweep (fit DensityNPMLE on total density → init_beliefs → node_sweep)
with `transfer_variance` toggled, and scores the deconvolved per-REGION f_g against the oracle truth
(mass-weighted |Δf_g|). This is the "run a real solve and look at the results" step — Phase 2 of
docs/calibration/npmle_projection_variance_design.md.

    OMP_NUM_THREADS=1 python sigma2_transfer_ab.py [--conditions a,b]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import sys

sys.path.insert(0, str(Path("/Users/mkiyer/proj/rigel/scripts/debug")))
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION, node_global_geometry, node_sweep
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from rigel.calibration.npmle import DensityNPMLE
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate,
    fit_rna_strand_from_substrate,
    overdispersion_for_beta,
)
from rigel.calibration.node_chain import build_node_chain
from rigel.calibration.node_geometry import build_node_geometry, build_node_statics, init_beliefs
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12


def solve(inp, index, cfg, transfer_variance):
    cc = cfg.calibration
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    pl = inp["payload"]
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    gfl, rfl = np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"])
    reg_el = region_eff_length(ra.region_size_bp, gfl)
    bnd_el = boundary_side_eff_length(gfl, ra.region_size_bp)
    fl_mean = boundary_eff_length(gfl)
    kappa = float(fit_strand_balance(inp["strand_model"]).rna_sense_frac)
    nd = node_gdna_density(sub, ra, reg_el, fl_mean)
    od_g = fit_gdna_strand_from_substrate(
        sub, ra, nd, bnd_el, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(cc.gdna_strand_prior_alpha_beta),
        prior_weight=cc.gdna_strand_prior_weight,
    ).gdna_strand_overdispersion
    od_r = fit_rna_strand_from_substrate(
        sub, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(cc.rna_strand_prior_alpha_beta),
        prior_weight=cc.rna_strand_prior_weight,
    ).rna_strand_overdispersion

    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    geom = build_node_geometry(chain, sub, bsub, ra, gfl, rfl)
    st = build_node_statics(chain, sub, bsub, ra)
    b0 = init_beliefs(
        chain, sub, bsub, ra, rna_sense_frac=kappa, gdna_strand_overdispersion=od_g,
        rna_strand_overdispersion=od_r, n_grid=cc.sweep_n_grid,
        n_grid_ss=cc.sweep_n_grid_single_strand, logodds_window=cc.sweep_logodds_window, statics=st,
    )
    mass_g, eff_g = node_global_geometry(chain, geom)
    prior = DensityNPMLE.fit(mass_g, eff_g, bandwidth=cc.npmle_bandwidth)
    belief = node_sweep(
        chain, st, geom, b0, ra, bsub, rna_sense_frac=kappa, gdna_strand_overdispersion=od_g,
        rna_strand_overdispersion=od_r, n_grid=cc.sweep_n_grid, logodds_window=cc.sweep_logodds_window,
        n_tilt=cc.sweep_n_tilt, n_grid_ss=cc.sweep_n_grid_single_strand, gdna_prior=prior,
        transfer_variance=transfer_variance,
    )
    # per-REGION f_g vs oracle true_fg (unspliced mass basis)
    is_reg = np.asarray(chain.kind) == REGION
    ridx = np.asarray(chain.ref_idx, np.int64)[is_reg]
    fg = np.asarray(belief.f_g, float)[is_reg]
    p = inp["region_pools"]
    G = np.asarray(p["gdna_pos"], float) + np.asarray(p["gdna_neg"], float)
    R = sum(np.asarray(p[k], float) for k in ("mat_uns_pos", "mat_uns_neg", "nas_uns_pos", "nas_uns_neg"))
    G, R = G[ridx], R[ridx]
    tot = G + R
    true_fg = np.where(tot > _EPS, G / np.maximum(tot, _EPS), np.nan)
    ok = np.isfinite(true_fg) & (tot > _EPS)
    w = tot[ok]
    mwae = float((w * np.abs(fg[ok] - true_fg[ok])).sum() / max(w.sum(), 1.0))
    # signed gDNA-mass error (net leak): cal gDNA mass − true gDNA mass
    cal_g = fg[ok] * tot[ok]
    net = float((cal_g - G[ok]).sum())
    return mwae, net, float(G.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    conds = a.conditions.split(",") if a.conditions else sorted(p.stem for p in cache.glob("*.pkl"))
    rows = []
    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        mwae_off, net_off, gtot = solve(inp, index, cfg, transfer_variance=False)
        mwae_on, net_on, _ = solve(inp, index, cfg, transfer_variance=True)
        cap = "capON" if "capture_on" in c else "capOFF"
        nas = "+nas" if "nrna_present" in c else "    "
        rows.append(dict(cond=c.replace("gdna_gdna300_", ""), cap=cap, nas=nas,
                         mwae_OFF=mwae_off, mwae_ON=mwae_on, d_mwae=mwae_on - mwae_off,
                         net_OFF=net_off, net_ON=net_on, true_g=gtot))
    df = pd.DataFrame(rows)
    print(df.to_string(index=False, float_format=lambda x: f"{x:,.4g}"))
    print(f"\nMEAN mwae:  OFF={df.mwae_OFF.mean():.4f}  ON={df.mwae_ON.mean():.4f}  "
          f"Δ={df.mwae_ON.mean() - df.mwae_OFF.mean():+.4f}")


if __name__ == "__main__":
    main()
