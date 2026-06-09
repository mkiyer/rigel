#!/usr/bin/env python3
"""TEMP DEBUG / DRY-RUN — prototype the local boundary-anchored density imputation.

Replaces the undamped global sweep with a LOCAL, directional inward imputation and compares two
density estimators against oracle truth, per region:
  A1  side_gdna_mass / boundary_side_eff_len            (user's stated form)
  A2  side_gdna_count / fl_mean                         (count-rate form)
Plus the run-fill: regions with no observable boundary side (run interiors) inherit the nearest
anchored neighbour. Observable regions use their own contained density.

Goal: see whether the local estimate recovers the true per-region gDNA density (it should be ~10-30×
the broken global), and which normalization is closest / unbiased. Surfaces implementation issues.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pysam

from rigel.config import BamScanConfig, CalibrationConfig
from rigel.index import TranscriptIndex
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.effective_length import region_eff_length, boundary_eff_length, boundary_side_eff_length
from rigel.calibration.density_model import node_gdna_density, count_observable_masks
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.signature import TS_NEG, TS_NONE, TS_AMBIG, coarse_type_array
from rigel.pipeline import scan_and_buffer, _check_region_payload_alignment
from rigel.splice import SpliceType
from rigel.sim.read_name import parse_origin

_TYPE = {0: "interg", 1: "intron", 2: "EXON"}


def _clean_frac(sense, total, kappa):
    sfrac = np.where(total > 0, sense / np.maximum(total, 1e-9), 0.5)
    d = 0.5 - kappa
    if abs(d) <= 1e-6:
        return np.ones_like(sfrac)
    return np.clip((sfrac - kappa) / d, 0.0, 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True, type=Path)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--start", required=True, type=int)
    ap.add_argument("--end", required=True, type=int)
    args = ap.parse_args()

    index = TranscriptIndex.load(args.index)
    _s, sm, fla, _b, pl = scan_and_buffer(str(args.bam), index, BamScanConfig(sj_strand_tag="auto"))
    sm.finalize()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    _check_region_payload_alignment(ra, pl)
    flm = build_fl_models(global_counts=fla.global_model.counts,
                          rna_counts=fla.category_models[SpliceType.SPLICED_ANNOT].counts,
                          gdna_counts=gdna_fl_mass(pl), max_size=fla.max_size)
    gpmf = flm.gdna_pmf
    cal = calibrate(pl, ra, sm, gpmf, CalibrationConfig())
    sub = CalibrationSubstrate.from_payload(pl, ra)
    reg_eff = region_eff_length(ra.region_size_bp, gpmf)
    bside = boundary_side_eff_length(gpmf, ra.region_size_bp)
    fl_mean = boundary_eff_length(gpmf)
    kappa = cal.rna_sense_frac
    ts = np.asarray(ra.strand_class)
    ctype = coarse_type_array(np.asarray(ra.signature))
    reg_obs, bnd_obs = count_observable_masks(np.asarray(ra.signature), np.asarray(ra.ref_id))
    n = ra.n_regions

    def side_sense(view):
        return np.where(ts == TS_NEG, view.n_unspliced_neg, view.n_unspliced_pos).astype(np.float64)

    # per-region cleaned contained gDNA (for observable regions)
    c = sub.contained
    cont_tot = (c.n_unspliced_pos + c.n_unspliced_neg).astype(np.float64)
    cont_gf = np.where(ts == TS_AMBIG, 1.0, _clean_frac(side_sense(c), cont_tot, kappa))
    cont_gf = np.where(ts == TS_NONE, 1.0, cont_gf)
    cont_gdna_mass = cont_gf * c.mass_unspliced
    cont_gdna_cnt = cont_gf * cont_tot

    # per-side cleaned gDNA (left side = right of left boundary; right side = left of right boundary)
    lgf = np.where(ts == TS_AMBIG, 1.0, _clean_frac(side_sense(sub.left), sub.left.n_unspliced.astype(float), kappa))
    rgf = np.where(ts == TS_AMBIG, 1.0, _clean_frac(side_sense(sub.right), sub.right.n_unspliced.astype(float), kappa))
    left_g_mass = lgf * sub.left.mass_unspliced
    left_g_cnt = lgf * sub.left.n_unspliced.astype(float)
    right_g_mass = rgf * sub.right.mass_unspliced
    right_g_cnt = rgf * sub.right.n_unspliced.astype(float)

    # boundary-observability per region SIDE: left side usable if boundary (r-1,r) observable; right side if (r,r+1)
    left_bnd_obs = np.zeros(n, bool)
    right_bnd_obs = np.zeros(n, bool)
    rid = np.asarray(ra.ref_id)
    left_bnd_obs[1:] = bnd_obs[:-1] & (rid[1:] == rid[:-1])
    right_bnd_obs[:-1] = bnd_obs[:-1] & (rid[:-1] == rid[1:])

    # impute density per region (A1 mass-based, A2 count-based)
    dens_A1 = np.full(n, np.nan)
    dens_A2 = np.full(n, np.nan)
    for r in range(n):
        if reg_obs[r]:
            dens_A1[r] = cont_gdna_mass[r] / max(reg_eff[r], 1e-9)
            dens_A2[r] = cont_gdna_cnt[r] / max(reg_eff[r], 1e-9)
            continue
        m = []
        cnt = []
        if left_bnd_obs[r]:
            m.append(left_g_mass[r] / max(bside[r], 1e-9))
            cnt.append(left_g_cnt[r] / max(fl_mean, 1e-9))
        if right_bnd_obs[r]:
            m.append(right_g_mass[r] / max(bside[r], 1e-9))
            cnt.append(right_g_cnt[r] / max(fl_mean, 1e-9))
        if m:
            dens_A1[r] = float(np.mean(m))
            dens_A2[r] = float(np.mean(cnt))
    # run-fill: regions still NaN (no observable side) inherit nearest anchored neighbour (both dirs, avg)
    for arr in (dens_A1, dens_A2):
        # L->R carry then R->L carry, average the two where interior
        carryL = arr.copy()
        for r in range(1, n):
            if np.isnan(carryL[r]) and rid[r] == rid[r-1]:
                carryL[r] = carryL[r-1]
        carryR = arr.copy()
        for r in range(n - 2, -1, -1):
            if np.isnan(carryR[r]) and rid[r] == rid[r+1]:
                carryR[r] = carryR[r+1]
        both = np.where(np.isnan(arr), np.nanmean(np.vstack([carryL, carryR]), axis=0), arr)
        arr[:] = both

    # oracle true contained gDNA density per region (template span)
    rref = index.ref_name_to_id[args.ref]
    win = np.where((rid == rref) & (ra.end > args.start) & (ra.start < args.end))[0]
    ws, we = ra.start[win], ra.end[win]
    tg = np.zeros(len(win))
    with pysam.AlignmentFile(str(args.bam), "rb") as b:
        for r in b:
            if r.is_read2 or r.is_secondary or r.is_supplementary or r.reference_name != args.ref:
                continue
            fl = abs(r.template_length) or ((r.reference_end or r.reference_start) - r.reference_start)
            s = min(r.reference_start, r.next_reference_start) if r.template_length else r.reference_start
            e = s + fl
            if e <= args.start or s >= args.end or parse_origin(r.query_name).kind != "gdna":
                continue
            j = np.searchsorted(we, s, side="right")
            if 0 <= j < len(win) and ws[j] <= s and e <= we[j]:
                tg[j] += 1

    nd = node_gdna_density(sub, ra, reg_eff, fl_mean, gdna_frac=cont_gf)
    print(f"\n=== IMPUTE DRY-RUN {args.ref}:{args.start:,}-{args.end:,}  (kappa={kappa:.3f}, fl_mean={fl_mean:.0f}) ===")
    print(f"{'reg':>4} {'type':>6} {'obs':>3} {'Lobs':>4} {'Robs':>4} | {'true_dens':>9} {'CURRENT':>8} | {'A1(mass)':>8} {'A2(cnt)':>8}")
    for k, i in enumerate(win):
        td = tg[k] / max(reg_eff[i], 1e-9)
        print(f"{i:>4} {_TYPE[ctype[i]]:>6} {'Y' if reg_obs[i] else 'n':>3} "
              f"{'Y' if left_bnd_obs[i] else 'n':>4} {'Y' if right_bnd_obs[i] else 'n':>4} | "
              f"{td:>9.3f} {nd.density[i]:>8.3f} | {dens_A1[i]:>8.3f} {dens_A2[i]:>8.3f}")


if __name__ == "__main__":
    main()
