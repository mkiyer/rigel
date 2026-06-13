"""Dissect the flagship through the PROPAGATION calibration: per-region/side error vs oracle, traced.

Runs the propagation deconvolution (propagate -> regions, left, right) on a flagship condition and joins
each region to its ORACLE contained-unspliced origin truth. Surfaces:
  - per-strand-class deficit (oracle gDNA - calib gDNA);
  - the worst regions by gDNA mass error, with the full propagation internals (nascent fields, seedable,
    mature, residual) so each error traces to its source;
  - the boundary-side gDNA totals.

Usage:  python scripts/debug/propagation_flagship_dissect.py [condition] [n_worst]
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
from rigel.calibration.propagation import _compute_fields, _solve_regions, _solve_sides
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import (
    BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS, TS_AMBIG, TS_NEG, TS_NONE, TS_POS,
)
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
SC = {TS_POS: "POS", TS_NEG: "NEG", TS_NONE: "NONE", TS_AMBIG: "AMBIG"}
_SIG = ((BIT_EXON_POS, "E+"), (BIT_EXON_NEG, "E-"), (BIT_INTRON_POS, "I+"), (BIT_INTRON_NEG, "I-"))


def decode(s):
    return "|".join(t for b, t in _SIG if s & b) or "."


def oracle_contained(bam, starts, ids, R):
    g, m, n = np.zeros(R), np.zeros(R), np.zeros(R)
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
            k = parse_origin(r.query_name).kind
            (g if k == "gdna" else n if k == "nrna" else m)[rid] += 1.0
    return g, m, n


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
    n_worst = int(sys.argv[2]) if len(sys.argv) > 2 else 20
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

    fields = _compute_fields(sub, ra, e_rna, rna_fl_mean, kappa, od_g, od_r, ccfg.n_grid)
    regions = _solve_regions(sub, ra, fields, e_rna, kappa, od_g, od_r, nd.count_gdna_frac, ccfg.n_grid)
    left, right = _solve_sides(sub, ra, fields, rna_fl_mean, kappa, od_g, od_r, nd.count_gdna_frac,
                               ccfg.n_grid)

    df = idx.region_df
    starts = {r: g["start"].to_numpy() for r, g in df.groupby("ref_name")}
    ids = {r: g["region_id"].to_numpy() for r, g in df.groupby("ref_name")}
    R = len(df)
    og, om, on = oracle_contained(bam, starts, ids, R)
    o_tot = og + om + on
    o_gfrac = np.where(o_tot > 0, og / np.maximum(o_tot, 1e-9), np.nan)

    ts = np.asarray(ra.strand_class)
    sig = np.asarray(ra.signature).astype(np.int64)
    calib_g = np.asarray(regions.gdna_mass)
    fg = np.asarray(regions.gdna_frac)
    rid = df["region_id"].to_numpy()
    s = df["start"].to_numpy()
    e = df["end"].to_numpy()
    refn = df["ref_name"].to_numpy()
    nasc_pos = np.where(fields.seedable_pos, fields.nasc_density_pos * e_rna, 0.0)
    nasc_neg = np.where(fields.seedable_neg, fields.nasc_density_neg * e_rna, 0.0)

    print(f"\n===== {cond}: PROPAGATION calibration dissection =====")
    print(f"kappa={kappa:.4f} od_g={od_g:.4f} od_r={od_r:.4f}")
    print(f"TOTALS contained-unspl: oracle gDNA={og.sum():.0f} mRNA={om.sum():.0f} nRNA={on.sum():.0f}  "
          f"calib gDNA={calib_g.sum():.0f}  deficit={og.sum()-calib_g.sum():+.0f}")
    lg = np.asarray(left.gdna_mass).sum() + np.asarray(right.gdna_mass).sum()
    print(f"  side gDNA total = {lg:.0f}")
    deficit = og - calib_g
    print("  deficit by region class (oracle_gdna - calib_gdna):")
    for cv, nm in ((TS_AMBIG, "AMBIG"), (TS_POS, "POS"), (TS_NEG, "NEG"), (TS_NONE, "NONE")):
        m = ts == cv
        tot = deficit.sum()
        print(f"    {nm:>5}: n={int(m.sum()):>5} oracle_gdna={og[m].sum():>9.0f} calib={calib_g[m].sum():>9.0f}"
              f" deficit={deficit[m].sum():>+9.0f} ({100*deficit[m].sum()/max(tot,1):>4.0f}%)")

    order = np.argsort(-np.abs(deficit))
    print(f"\n--- {n_worst} regions, largest |gDNA deficit| (propagation) ---")
    print(f"{'ref:start':>16}{'sig':>6}{'cls':>6}{'o_g':>7}{'o_m':>7}{'o_n':>7}{'ogfr':>6}{'U':>7}"
          f"{'nascP':>6}{'nascN':>6}{'sdP':>4}{'sdN':>4}{'M+':>6}{'M-':>6}{'fg':>5}{'cgd':>7}{'defic':>7}")
    for i in order[:n_worst]:
        c = sub.contained
        U = float(c.n_unspliced_pos[i] + c.n_unspliced_neg[i])
        print(f"{f'{refn[i]}:{s[i]}':>16}{decode(int(sig[i])):>6}{SC[ts[i]]:>6}{og[i]:>7.0f}{om[i]:>7.0f}"
              f"{on[i]:>7.0f}{o_gfrac[i]:>6.2f}{U:>7.0f}{nasc_pos[i]:>6.0f}{nasc_neg[i]:>6.0f}"
              f"{int(fields.seedable_pos[i]):>4}{int(fields.seedable_neg[i]):>4}"
              f"{fields.mature_pos[i]:>6.0f}{fields.mature_neg[i]:>6.0f}{fg[i]:>5.2f}{calib_g[i]:>7.0f}"
              f"{deficit[i]:>+7.0f}")


if __name__ == "__main__":
    main()
