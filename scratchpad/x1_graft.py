"""P1d — the EXTRAPOLATION premise at the graft, measured against the oracle.

THE GRAFT'S CLAIM. A boundary hands the exon next door its RNA: the RNA that CROSSES the seam unspliced
(rho_nu) plus the RNA that SPLICES at the junction (rho_mu, a direct measurement). The premise is that those
two together are the exon's whole RNA density:

    rho_R(exon)  ==  [ rho_nu(bnd) + rho_mu(bnd) ]  /  (the capture step between them)

That is exact for an internal exon of a single isoform: everything in the exon entered through this seam. It
is NOT exact when some of the exon's RNA never passes this seam -- unspliced/nascent RNA that has not spliced
yet, transcripts that start or end inside the exon, or isoforms that reach the exon by a different junction.

So define the GRAFT SHARE, the exact mirror of the peel's continuing share w:

    phi  =  [ rho_nu(bnd) + rho_mu(bnd) ]  /  rho_R(exon)        in a COMMON frame,  phi <= 1

phi = 1 => the premise holds and the graft needs no extrapolation term at all, only its count and frame
variance. phi < 1 => the graft SYSTEMATICALLY UNDER-STATES the exon's RNA by 1/phi, which is a BIAS, and a
variance term cannot fix a bias.

THE FRAME. Both sides are put in a common frame using the ORACLE capture step. gDNA content is uniform along
the genome, so the oracle gDNA density ratio IS the capture ratio -- no model estimate enters.

    OMP_NUM_THREADS=1 python scratchpad/x1_graft.py [cond ...]
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
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
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

print(f"{'condition':<44}{'edges':>7}{'phi p25':>9}{'phi MED':>9}{'phi p75':>9}"
      f"{'  1/phi bias (nats)':>20}{'sd(log phi)':>13}{'2-junction |gap|':>18}")
print("-" * 130)
for cond in CONDS:
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
    n = kind.shape[0]

    def pool(k):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    G = pool("gdna_pos") + pool("gdna_neg")
    Ru = pool("mat_uns_pos") + pool("nas_uns_pos") + pool("mat_uns_neg") + pool("nas_uns_neg")
    Rs = pool("mat_spl") + pool("nas_spl")
    E_g, E_r = us["E_g"], us["E_r"]
    li, ri = us["left"], us["right"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    SPf, SNf = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    rt, _ = _node_region_type(chain, ra)

    # ORACLE densities. gDNA is uniform along the genome, so its density ratio IS the capture step.
    rho_g = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)
    # the exon's RNA density: unspliced RNA it contains PLUS the spliced RNA whose blocks lie in it
    rho_R_ex = np.where(E_r > _EPS, (Ru + Rs) / np.maximum(E_r, _EPS), np.nan)
    rho_nu_b = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)

    phis, gaps = [], []
    for face, nbr in ((1, li), (0, ri)):  # exon i, its neighbour boundary on that side
        for i in np.flatnonzero(is_exon):
            b = nbr[i]
            if b < 0 or not is_bnd[b]:
                continue
            mu = (SPf[face][b] + SNf[face][b]) / max(ESP[face][b], _EPS)
            if not (mu > _EPS) or not np.isfinite(rho_R_ex[i]) or rho_R_ex[i] <= _EPS:
                continue
            if not (rho_g[b] > _EPS and rho_g[i] > _EPS):
                continue
            step = rho_g[b] / rho_g[i]  # oracle capture step boundary→exon
            phis.append((rho_nu_b[b] + mu) / (rho_R_ex[i] * step))
    # the candidate prior-free estimator: do the exon's TWO flanking junctions agree?
    for i in np.flatnonzero(is_exon):
        lb, rb = li[i], ri[i]
        if lb < 0 or rb < 0:
            continue
        ml = (SPf[1][lb] + SNf[1][lb]) / max(ESP[1][lb], _EPS)
        mr = (SPf[0][rb] + SNf[0][rb]) / max(ESP[0][rb], _EPS)
        if ml > _EPS and mr > _EPS:
            gaps.append(abs(np.log(ml / mr)))
    phis = np.asarray([p for p in phis if np.isfinite(p) and p > 0])
    gaps = np.asarray(gaps)
    if phis.size < 5:
        print(f"{cond[5:]:<44}{phis.size:>7}   (too few)")
        continue
    lp = np.log(phis)
    print(f"{cond[5:]:<44}{phis.size:>7}"
          f"{np.percentile(phis, 25):>9.3f}{np.median(phis):>9.3f}{np.percentile(phis, 75):>9.3f}"
          f"{-np.median(lp):>20.3f}{np.std(lp):>13.3f}"
          f"{(np.median(gaps) if gaps.size else np.nan):>18.3f}")
