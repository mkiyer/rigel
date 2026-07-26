"""AMBIG study, step 7 — the OWNER'S CORRECTION: kappa is measured, has a distribution, and (kappa-1/2)
in a denominator is degenerate.

The plug-in estimator d_hat = (p_obs - 1/2)/(kappa - 1/2) divides by a SMALL, UNCERTAIN quantity. Two errors
compound:
  (a) p_obs has its own sampling spread even for pure gDNA -- 1000 gDNA counts split 506/494, not 500/500
      (beta-binomial). A tiny |p_obs - 1/2| becomes a large |d_hat| when divided by a tiny (kappa - 1/2).
  (b) kappa itself is FITTED from spliced counts, with posterior Beta(n_same+1, n_opp+1). So
      Var(kappa) = kappa(1-kappa)/(n_obs+3) -- which is ALREADY computed as
      `StrandBalance.rna_strand_overdispersion` and currently labelled "QC-only, NOT fed into the deconv".

This measures, per condition: kappa_hat, sd(kappa), the z-score |kappa-1/2|/sd(kappa) (is the library
stranded AT ALL?), and what the plug-in d_hat does. It also reports the NON-DEGENERATE form of the same
constraint,  (1 - f_g)*|kappa - 1/2| >= |p_obs - 1/2|,  which never divides.

    OMP_NUM_THREADS=1 python scratchpad/ambig_7_kappa_degeneracy.py
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

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.strand_balance import fit_strand_balance  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna_none_ss_0.50_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

print("kappa is FITTED: posterior Beta(n_same+1, n_opp+1) ⇒ Var(kappa) = kappa(1-kappa)/(n_obs+3).")
print(f"{'condition':<42}{'kappa':>9}{'sd(kap)':>10}{'|k-.5|/sd':>11}{'n_obs':>9}"
      f"|{'med|d_hat|':>11}{'%|d|>1':>8}{'medB':>7}{'orcFg':>7}|{'A=|k-.5|':>9}{'med|e|':>8}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    sb = fit_strand_balance(inp["strand_model"])
    dbg: dict = {}
    res = calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    cap, st = dbg["capture"], dbg["statics"]
    chain = dbg["chain"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    T = G + R
    fo = np.where(T > _EPS, G / np.maximum(T, _EPS), np.nan)
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    amb = fp & fn
    up, un = np.asarray(st.u_pos, float), np.asarray(st.u_neg, float)
    n = up + un
    mass = np.asarray(cap["mass_global"])
    kap = float(res.rna_sense_frac)
    # kappa's OWN posterior sd, from the same Beta the mean came from
    n_obs = float(sb.n_observations)
    sd_kap = np.sqrt(kap * (1.0 - kap) / (n_obs + 3.0))
    z_lib = abs(kap - 0.5) / max(sd_kap, _EPS)
    p_obs = np.where(n > 0, up / np.maximum(n, _EPS), 0.5)
    e = p_obs - 0.5  # the observed strand EXCESS -- the non-degenerate statistic
    den = kap - 0.5
    d_hat = np.where(abs(den) > _EPS, e / (den if abs(den) > _EPS else 1.0), 0.0)
    ok = np.isfinite(fo) & (mass > _EPS) & amb & (n > 0)
    w = mass[ok]

    def q(x, p=0.5):
        o = np.argsort(x)
        return x[o][np.searchsorted(np.cumsum(w[o]) / w[o].sum(), p)]

    B = np.clip(1.0 - np.abs(d_hat), 0.0, 1.0)
    print(f"{cond[5:]:<42}{kap:>9.5f}{sd_kap:>10.5f}{z_lib:>11.1f}{n_obs:>9,.0f}"
          f"|{q(np.abs(d_hat[ok])):>11.3f}{np.average((np.abs(d_hat[ok]) > 1).astype(float), weights=w):>8.1%}"
          f"{q(B[ok]):>7.3f}{np.average(fo[ok], weights=w):>7.3f}"
          f"|{abs(den):>9.5f}{q(np.abs(e[ok])):>8.5f}")
