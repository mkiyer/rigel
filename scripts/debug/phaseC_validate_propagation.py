"""Validate src/rigel/calibration/propagation.py against the oracle on locus 21 (both nrna conditions).

Calls the production module (propagate_regions) and compares per-region gdna_frac to oracle truth, to
confirm it reproduces the v2 prototype (AMBIG mean ~0.056 on nrna_rnd) before wiring into calibrate.

Usage:  python scripts/debug/phaseC_validate_propagation.py [condition]
"""
import dataclasses
import sys

import numpy as np
import pysam

from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import (
    boundary_eff_length, boundary_side_eff_length, region_eff_length,
)
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate, fit_rna_strand_from_substrate, overdispersion_for_beta,
)
from rigel.calibration.propagation import propagate_regions
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"


def run(cond):
    bam = f"{SUITE}/{cond}/sim_oracle.bam"
    idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
    cfg = PipelineConfig()
    ccfg = cfg.calibration
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    reg_el = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    bnd_el = boundary_side_eff_length(fl.gdna_pmf, ra.region_size_bp)
    fl_mean = boundary_eff_length(fl.gdna_pmf)
    e_rna = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    rna_fl_mean = boundary_eff_length(fl.rna_pmf)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    nd = node_gdna_density(sub, ra, reg_el, fl_mean, need_count_variance=False)
    od_g = fit_gdna_strand_from_substrate(sub, ra, nd, bnd_el, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.gdna_strand_prior_alpha_beta),
        prior_weight=ccfg.gdna_strand_prior_weight).gdna_strand_overdispersion
    od_r = fit_rna_strand_from_substrate(sub, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.rna_strand_prior_alpha_beta),
        prior_weight=ccfg.rna_strand_prior_weight).rna_strand_overdispersion

    result = propagate_regions(
        sub, ra, rna_region_eff_len=e_rna, rna_fl_mean=rna_fl_mean, rna_sense_frac=kappa,
        gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r,
        count_gdna_frac=nd.count_gdna_frac, n_grid=ccfg.n_grid,
    )
    fg = np.asarray(result.gdna_frac)

    # oracle
    df = idx.region_df
    starts = {r: g["start"].to_numpy() for r, g in df.groupby("ref_name")}
    ids = {r: g["region_id"].to_numpy() for r, g in df.groupby("ref_name")}
    og = np.zeros(len(df))
    om = np.zeros(len(df))
    with pysam.AlignmentFile(bam, "rb") as b:
        for r in b.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1:
                continue
            if r.cigartuples and any(op == 3 for op, _ in r.cigartuples):
                continue
            ref = r.reference_name
            if ref not in starts:
                continue
            i = int(np.searchsorted(starts[ref], r.reference_start, side="right") - 1)
            if i < 0:
                continue
            rid = int(ids[ref][i])
            (og if parse_origin(r.query_name).kind == "gdna" else om)[rid] += 1
    o_fg = np.where(og + om > 0, og / np.maximum(og + om, 1e-9), np.nan)
    return ra, idx, fg, o_fg, sub


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_rnd_capture_on"
    ra, idx, fg, o_fg, sub = run(cond)
    ts = np.asarray(ra.strand_class)
    df = idx.region_df
    s = df["start"].to_numpy()
    e = df["end"].to_numpy()
    refn = df["ref_name"].to_numpy()
    rid = df["region_id"].to_numpy()
    U = (sub.contained.n_unspliced_pos + sub.contained.n_unspliced_neg).astype(float)
    sel = np.flatnonzero((refn == "chr_syn") & (e > 964416) & (s < 1004165) & (ts == TS_AMBIG) & (U > 200))
    print(f"=== {cond}: propagate_regions vs oracle (locus 21 AMBIG) ===")
    print(f"{'rid':>5}{'U':>9}{'oracle_fg':>11}{'prop_fg':>9}")
    for i in sel:
        print(f"{rid[i]:>5}{U[i]:>9.0f}{o_fg[i]:>11.2f}{fg[i]:>9.2f}")
    err = np.abs(fg[sel] - o_fg[sel])
    print(f"\n  AMBIG mean|prop_fg - oracle_fg| = {np.nanmean(err):.3f}  (v2 prototype: ~0.056 on nrna_rnd)")


if __name__ == "__main__":
    main()
