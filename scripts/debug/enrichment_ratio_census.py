"""E0 — ENRICHMENT-RATIO CENSUS. Is delta_e measurable, and with WHICH estimator?

delta_e = log[rho_tot(dst)/rho_tot(src)] is the edge's log-enrichment step. It is confounded by any
composition or content difference, so it is measurable exactly where the two endpoints carry the SAME ACTIVE
COMPONENT SET:

  A  intergenic region  <->  intergenic-exon/intron SEAM boundary     {g}      <-> {g}        MEASURABLE
  B  intron region      <->  intron-exon JUNCTION boundary            {g,R}    <-> {g,R}      MEASURABLE
                             ... boundary UNSPLICED crossing only (exclude spliced)
  C  exon region        <->  intron-exon JUNCTION boundary            {g,R,Rdep}<->{g,R,Rdep} MEASURABLE
                             ... boundary TOTAL, unspliced + SPLICED
  D  SEAM boundary      <->  exon region                              {g}      <-> {g,R,Rdep} NOT measurable

TWO ESTIMATORS, and they are NOT the same quantity (adversarial review, 2026-07-23):

  GDNA-FRAME   rho = M / E_g                      belief-free, composition-free. What the ORIGINAL census
                                                  measured (via `node_global_geometry`, whose eff is
                                                  eff_gdna_*). Cheap and robust, but it is NOT a total density.
  BLENDED      rho = M*[f_g/E_g + (1-f_g)/E_r]    the TOTAL density the enrichment framework actually needs
                                                  (docs/calibration/archive/enrichment_ratio_generalization.md §1). Composition-
                                                  dependent and taken from a SOLVE. The two coincide only at
                                                  f_g = 1.

An earlier draft claimed the framework's estimator was "already validated"; it was not -- the validation
covered the gDNA-frame estimator on the 20 UNSTRANDED scenarios only. This script closes that gap: BOTH
estimators, ALL 32 conditions, stranded and unstranded reported separately.

THE GATE: at capture OFF there is no enrichment, so a MEASURABLE family (A/B/C) must centre on ~0 under the
BLENDED estimator. If it does not, the framework's premise is wrong and the implementation must not proceed.

Usage:  OMP_NUM_THREADS=1 python scripts/debug/enrichment_ratio_census.py
"""
from __future__ import annotations

import dataclasses
import importlib
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

_EPS = 1e-12
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
INTERGENIC, INTRON, EXON = 0, 1, 2

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
CONDS = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

# (strandedness, capture, family, estimator) -> [delta]
acc: dict[tuple, list] = defaultdict(list)
n_scen = 0

for cond in CONDS:
    strand = "STRANDED  " if "ss_0.99" in cond else "unstranded"
    capst = "ON " if "capture_on" in cond else ("VSTR" if "verystrong" in cond else "OFF")
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
    )
    chain, geom, cap = dbg["chain"], dbg["geometry"], dbg["capture"]
    n_scen += 1
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    rtype, _ = _node_region_type(chain, ra)

    # node-level mass and per-component effective lengths (both faces pooled for a boundary)
    ml, mr = np.asarray(geom.mass_left, float), np.asarray(geom.mass_right, float)
    mass = np.where(is_reg, ml, ml + mr)
    egl, egr = np.asarray(geom.eff_gdna_left, float), np.asarray(geom.eff_gdna_right, float)
    erl, err_ = np.asarray(geom.eff_rna_left, float), np.asarray(geom.eff_rna_right, float)
    E_g = np.where(is_reg, egl, egl + egr)
    E_r = np.where(is_reg, erl, erl + err_)
    f_g = np.clip(np.asarray(cap["f_g"], float), 0.0, 1.0)

    # spliced DENSITY, per face (A3: spliced lands on ONE face, divide by THAT face's opportunity)
    spl_l = np.asarray(geom.spliced_pos_left, float) + np.asarray(geom.spliced_neg_left, float)
    spl_r = np.asarray(geom.spliced_pos_right, float) + np.asarray(geom.spliced_neg_right, float)
    esp_l = np.maximum(np.asarray(geom.eff_spl_left, float), _EPS)
    esp_r = np.maximum(np.asarray(geom.eff_spl_right, float), _EPS)
    rho_spl = np.where(spl_l > _EPS, spl_l / esp_l, 0.0) + np.where(spl_r > _EPS, spl_r / esp_r, 0.0)

    ok = (mass > _EPS) & (E_g > _EPS) & (E_r > _EPS)
    # ESTIMATOR 1 — gDNA-frame (what the original census measured): belief-free
    gdna_unspl = np.where(ok, mass / np.maximum(E_g, _EPS), np.nan)
    # ESTIMATOR 2 — BLENDED total density (what the framework needs): composition-dependent
    blend_unspl = np.where(ok, mass * (f_g / np.maximum(E_g, _EPS) + (1.0 - f_g) / np.maximum(E_r, _EPS)), np.nan)

    EST = {
        "gDNA-frame M/E_g": (gdna_unspl, gdna_unspl + rho_spl),
        "BLENDED rho_tot": (blend_unspl, blend_unspl + rho_spl),
    }

    for b in np.where(~is_reg)[0]:
        lf, rg = left[b], right[b]
        if lf < 0 or rg < 0:
            continue
        tl, tr = rtype[lf], rtype[rg]
        for reg, other in ((lf, tr), (rg, tl)):
            t = rtype[reg]
            if t == INTERGENIC:
                fam, use_total = "A intgc-reg <-> seam-bnd", False
            elif t == INTRON and other == EXON:
                fam, use_total = "B intron-reg <-> junc (unspl)", False
            elif t == EXON and other == INTRON:
                fam, use_total = "C exon-reg <-> junc (TOTAL)", True
            elif t == EXON and other == INTERGENIC:
                fam, use_total = "D exon-reg <-> seam  [UNMEAS]", True
            else:
                continue
            for est_name, (unspl, total) in EST.items():
                arr = total if use_total else unspl
                rs, rb = arr[reg], arr[b]
                if np.isfinite(rs) and np.isfinite(rb) and rs > _EPS and rb > _EPS:
                    acc[(strand, capst, fam, est_name)].append(float(np.log(rb) - np.log(rs)))

FAMS = [
    "A intgc-reg <-> seam-bnd",
    "B intron-reg <-> junc (unspl)",
    "C exon-reg <-> junc (TOTAL)",
    "D exon-reg <-> seam  [UNMEAS]",
]
ESTS = ["gDNA-frame M/E_g", "BLENDED rho_tot"]

print(f"\nE0 CENSUS — {n_scen} conditions (20 unstranded + 12 stranded), BOTH estimators.")
print("delta = log rho(boundary) - log rho(region). A MEASURABLE family must centre on ~0 at capture OFF.\n")
for strand in ("unstranded", "STRANDED  "):
    for capst in ("OFF", "ON ", "VSTR"):
        rows = [(f, e) for f in FAMS for e in ESTS if acc.get((strand, capst, f, e))]
        if not rows:
            continue
        print(f"  --- {strand.strip()} / capture {capst.strip()} ---")
        print(f"      {'family':<30} | {'estimator':<18} | {'n':>6} {'median':>9} {'mean':>9} {'sd':>8}")
        for f, e in rows:
            v = np.asarray(acc[(strand, capst, f, e)])
            v = v[np.isfinite(v)]
            if v.size < 4:
                continue
            print(f"      {f:<30} | {e:<18} | {v.size:>6} {np.median(v):>9.3f} {v.mean():>9.3f} {v.std():>8.3f}")
        print()

print("=" * 100)
print("GATE — capture OFF, MEASURABLE families (A/B/C) must centre on ~0 under the BLENDED estimator:")
fail = []
for strand in ("unstranded", "STRANDED  "):
    for f in FAMS[:3]:
        v = np.asarray(acc.get((strand, "OFF", f, "BLENDED rho_tot"), []))
        v = v[np.isfinite(v)]
        if v.size < 4:
            continue
        med = float(np.median(v))
        ok = abs(med) < 0.25  # a factor of ~1.28 -- the bounding lemma allows up to ~1.5x from composition
        print(f"   {strand.strip():<11} {f:<30} median {med:>+7.3f}   {'PASS' if ok else 'FAIL'}")
        if not ok:
            fail.append((strand.strip(), f, med))
    v = np.asarray(acc.get((strand, "OFF", FAMS[3], "BLENDED rho_tot"), []))
    v = v[np.isfinite(v)]
    if v.size >= 4:
        print(f"   {strand.strip():<11} {FAMS[3]:<30} median {float(np.median(v)):>+7.3f}   (control: must be FAR from 0)")
print()
print("VERDICT:", "GATE PASSED — the blended estimator is unbiased on the measurable families."
      if not fail else f"GATE FAILED on {len(fail)}: {fail}")
