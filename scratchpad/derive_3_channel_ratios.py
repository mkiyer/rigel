"""DERIVATION STEP 3 (§4.2) — measure the THREE per-channel enrichment ratios at a boundary face.

OWNER'S MODEL. A boundary does not sit on a vertical cliff; it sits on an enrichment/depletion FACE, and the
slope depends on probe placement. Its observed mass has TWO physically distinct channels:

  * UNSPLICED crossing  — fragments straddling the exon<->intron seam. Part of the footprint lies in the exon
                          (baited) and part in the intron (not baited) => an INTERMEDIATE capture exposure.
  * SPLICED crossing    — junction fragments, both blocks inside exons => the EXON's exposure, and HIGHER
                          still if a probe spans the junction (then spliced can exceed unspliced).

So we may assume neither `rho(boundary) == rho(exon)` on the spliced side, nor `rho(boundary) == rho(intron)`
on the unspliced side. The solver has ONE ratio, r = rho_tot(dst face)/rho_tot(src face), and `_rho_faces`
folds the one-sided spliced density INTO that face total — mixing the two channels into a single number.

HOW WE MEASURE THE TRUTH. Genomic DNA is uniformly present along the genome, so ANY ratio of true gDNA
densities IS a ratio of capture efficiencies:

    e(B^unspliced)/e(X) = [G(B)/E_g(B)] / [G(X)/E_g(X)]        <- boundary unspliced vs the EXON
    e(B^unspliced)/e(I) = [G(B)/E_g(B)] / [G(I)/E_g(I)]        <- boundary unspliced vs the INTRON
    e(B^spliced)/e(X)   = [S(B)/E_spl(B)] / [R(X)/E_r(X)]      <- the spliced channel vs the exon's RNA

The first two are pure capture measurements (gDNA has no expression component). The third is the one already
measured in the reconciliation doc (0.97 / 0.87 / 0.54 as capture strengthens).

    OMP_NUM_THREADS=1 python scratchpad/derive_3_channel_ratios.py
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


def q(x):
    x = np.asarray(x, np.float64)
    x = x[np.isfinite(x) & (x > 0)]
    return (np.median(x), np.percentile(x, 25), np.percentile(x, 75), x.size) if x.size else (np.nan,) * 3 + (0,)


for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, geom = dbg["chain"], dbg["geometry"]
    us = dbg["capture"]["_uni_static"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    Eg, Er = us["E_g"], us["E_r"]
    rho_l0, rho_r0 = us["rho_l0"], us["rho_r0"]
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    kind = np.asarray(chain.kind)
    rt, _ = _node_region_type(chain, ra)
    cls = np.array([CLS[int(rt[i])] if kind[i] == REGION else "bound" for i in range(G.size)])
    ESPl, ESPr = np.asarray(geom.eff_spl_left), np.asarray(geom.eff_spl_right)
    SPl, SPr = np.asarray(geom.spliced_pos_left), np.asarray(geom.spliced_pos_right)
    SNl, SNr = np.asarray(geom.spliced_neg_left), np.asarray(geom.spliced_neg_right)

    rg_true = np.where(Eg > _EPS, G / np.maximum(Eg, _EPS), np.nan)  # gDNA density = pure capture probe
    rr_true = np.where(Er > _EPS, R / np.maximum(Er, _EPS), np.nan)

    u_vs_exon, u_vs_intron, s_vs_exon, mdl_vs_true_exon, mdl_vs_true_intron = [], [], [], [], []
    for b in range(G.size):
        if kind[b] == REGION:
            continue
        lo, ro = int(left[b]), int(right[b])
        if lo < 0 or ro < 0:
            continue
        flank = {cls[lo]: lo, cls[ro]: ro}
        if not (rg_true[b] > _EPS):
            continue
        # the boundary's UNSPLICED channel vs each flank, in pure gDNA (= pure capture)
        if "exon" in flank and rg_true[flank["exon"]] > _EPS:
            x = flank["exon"]
            u_vs_exon.append(rg_true[b] / rg_true[x])
            # what the MODEL's face ratio says for that same edge (boundary face -> exon)
            fb = rho_r0[b] if x > b else rho_l0[b]
            fx = rho_l0[x] if x > b else rho_r0[x]
            if fb > _EPS and fx > _EPS:
                mdl_vs_true_exon.append((fb / fx) / (rg_true[b] / rg_true[x]))
            # the SPLICED channel vs the exon's RNA
            S = (SPr[b] + SNr[b]) if x > b else (SPl[b] + SNl[b])
            E = ESPr[b] if x > b else ESPl[b]
            if S > 20 and E > _EPS and rr_true[x] > _EPS and R[x] > 50:
                s_vs_exon.append((S / E) / rr_true[x])
        if "intron" in flank and rg_true[flank["intron"]] > _EPS:
            i2 = flank["intron"]
            u_vs_intron.append(rg_true[b] / rg_true[i2])
            fb = rho_r0[b] if i2 > b else rho_l0[b]
            fi = rho_l0[i2] if i2 > b else rho_r0[i2]
            if fb > _EPS and fi > _EPS:
                mdl_vs_true_intron.append((fb / fi) / (rg_true[b] / rg_true[i2]))

    print(f"\n{'=' * 100}\n{cond}\n{'=' * 100}")
    print(f"  {'TRUE per-channel enrichment (oracle gDNA = pure capture)':<52}{'median':>9}{'p25':>9}{'p75':>9}{'n':>7}")
    for lab, v in (
        ("e(boundary UNSPLICED) / e(EXON)", u_vs_exon),
        ("e(boundary UNSPLICED) / e(INTRON)", u_vs_intron),
        ("e(boundary SPLICED)  / e(EXON)   [RNA]", s_vs_exon),
    ):
        m, a, b_, n = q(v)
        print(f"  {lab:<52}{m:>9.3f}{a:>9.3f}{b_:>9.3f}{n:>7}")
    print(f"\n  {'MODEL face ratio / TRUE ratio (1.0 = the model has it right)':<52}{'median':>9}{'p25':>9}{'p75':>9}{'n':>7}")
    for lab, v in (
        ("boundary<->EXON edge", mdl_vs_true_exon),
        ("boundary<->INTRON edge", mdl_vs_true_intron),
    ):
        m, a, b_, n = q(v)
        print(f"  {lab:<52}{m:>9.3f}{a:>9.3f}{b_:>9.3f}{n:>7}")
