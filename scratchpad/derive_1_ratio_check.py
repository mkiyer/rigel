"""DERIVATION STEP 1 — is the mass-bound violation caused by a MIS-SPECIFIED r, or by scale-blindness?

This decides what to derive next (`docs/calibration/density_composition_reconciliation.md` §3, and the
orchestrator's flag that §4.2 may partly subsume §3). Two hypotheses:

  H1  THE RATIO IS WRONG.  The single r = rho_tot(dst)/rho_tot(src) conflates channels that have different
      capture exposure (the boundary is a SLOPE, not a cliff -- §4.2), so it mis-estimates the enrichment step
      and manufactures density. If so, the fix is a correct per-channel ratio, and a sigma^2_pin term would be
      pricing an error we could simply remove.

  H2  THE PREMISE IS WRONG.  r estimates capture correctly, but RNA does not follow capture -- RNA density is
      set by TRANSCRIPT ABUNDANCE, which is per-gene and not spatially smooth. If so, the reframe is doing its
      job and the failure is the imputation premise, which is what sigma^2_pin / b-hat^2 exist to price.

THE ORACLE GIVES US A GROUND TRUTH FOR r. Genomic DNA is uniformly present along the genome, so ANY difference
in true gDNA density between two nodes IS the capture efficiency difference:

        r_true^gDNA(dst<-src)  =  [G(dst)/E_g(dst)] / [G(src)/E_g(src)]      <- the true enrichment step
        r_true^RNA (dst<-src)  =  [R(dst)/E_r(dst)] / [R(src)/E_r(src)]      <- what RNA ACTUALLY does

Then:
  * corr(log r_model, log r_true^gDNA) high and slope ~1  =>  the reframe estimates capture WELL  => H2
  * r_true^RNA tracking r_true^gDNA                       =>  the imputation premise HOLDS
  * r_true^RNA scattering around r_true^gDNA              =>  the premise FAILS, by that much

    OMP_NUM_THREADS=1 python scratchpad/derive_1_ratio_check.py
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
_EPS = 1e-12
CLS = {0: "interg", 1: "intron", 2: "exon", -1: "bound"}

CONDS = [
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def _stats(x):
    x = x[np.isfinite(x)]
    return (np.median(x), np.percentile(x, 25), np.percentile(x, 75), x.size) if x.size else (np.nan,) * 3 + (0,)


for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap = dbg["chain"], dbg["capture"]
    us = cap["_uni_static"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    Eg, Er = us["E_g"], us["E_r"]
    rho_l0, rho_r0 = us["rho_l0"], us["rho_r0"]
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    kind = np.asarray(chain.kind)
    rt, _ = _node_region_type(chain, ra)
    cls = np.array([CLS[int(rt[i])] if kind[i] == REGION else "bound" for i in range(G.size)])

    rho_g_true = np.where(Eg > _EPS, G / np.maximum(Eg, _EPS), np.nan)
    rho_r_true = np.where(Er > _EPS, R / np.maximum(Er, _EPS), np.nan)

    lm, lg, lr, tag = [], [], [], []
    for i in range(G.size):
        for s, dstf, srcf in ((int(left[i]), rho_l0, rho_r0), (int(right[i]), rho_r0, rho_l0)):
            if s < 0:
                continue
            if not (dstf[i] > _EPS and srcf[s] > _EPS):
                continue
            if not (rho_g_true[i] > _EPS and rho_g_true[s] > _EPS):
                continue
            lm.append(np.log(dstf[i] / srcf[s]))
            lg.append(np.log(rho_g_true[i] / rho_g_true[s]))
            lr.append(
                np.log(rho_r_true[i] / rho_r_true[s])
                if (rho_r_true[i] > _EPS and rho_r_true[s] > _EPS)
                else np.nan
            )
            tag.append(f"{cls[s]}->{cls[i]}")
    lm, lg, lr, tag = np.array(lm), np.array(lg), np.array(lr), np.array(tag)

    print(f"\n{'=' * 96}\n{cond}\n{'=' * 96}")
    ok = np.isfinite(lm) & np.isfinite(lg)
    sl = float(np.polyfit(lg[ok], lm[ok], 1)[0])
    print(f"H1 — does the MODEL's r estimate the true capture step?   n={int(ok.sum()):,}")
    print(f"     corr(log r_model, log r_true^gDNA) = {np.corrcoef(lm[ok], lg[ok])[0, 1]:.3f}   slope = {sl:.3f}"
          f"   median |log r_model - log r_true^gDNA| = {np.median(np.abs(lm[ok] - lg[ok])):.3f}"
          f"  (= x{np.exp(np.median(np.abs(lm[ok] - lg[ok]))):.2f})")
    ok2 = ok & np.isfinite(lr)
    print(f"H2 — does RNA FOLLOW capture?                             n={int(ok2.sum()):,}")
    print(f"     corr(log r_true^RNA, log r_true^gDNA) = {np.corrcoef(lr[ok2], lg[ok2])[0, 1]:.3f}"
          f"   median |log r_true^RNA - log r_true^gDNA| = {np.median(np.abs(lr[ok2] - lg[ok2])):.3f}"
          f"  (= x{np.exp(np.median(np.abs(lr[ok2] - lg[ok2]))):.2f})")
    print(f"\n     {'edge type':<18}{'n':>6}{'|logr_mdl-logr_gDNA|':>22}{'|logr_RNA-logr_gDNA|':>22}")
    for t in ("exon->bound", "bound->exon", "intron->bound", "bound->intron", "bound->interg", "interg->bound"):
        m = ok & (tag == t)
        m2 = ok2 & (tag == t)
        if m.sum() < 5:
            continue
        a = np.median(np.abs(lm[m] - lg[m]))
        b = np.median(np.abs(lr[m2] - lg[m2])) if m2.sum() >= 5 else np.nan
        print(f"     {t:<18}{int(m.sum()):>6}{a:>16.3f} (x{np.exp(a):.1f}){b:>16.3f} (x{np.exp(b):.1f})")
