"""Phase-3a flagship forensic debug: gdna_gdna300_ss_0.99_nrna_none_capture_on.

Replicates calibrate's REGION path step-by-step (the A/B showed the regression is the region
deconvolution, not the splice or sides), capturing every per-region intermediate, and joins each region
to its ORACLE truth (contained-unspliced fragment origins). Surfaces the regions where the calibration
most under-calls gDNA (the gDNA→RNA leak), with the full cleaning calculation so we can see PRECISELY
why a true-gDNA region is read as RNA.

Usage:  python scripts/debug/phase3_flagship_debug.py [condition] [n_worst]
"""
import dataclasses
import sys

import numpy as np
import pysam

from rigel.calibration.effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate,
    fit_rna_strand_from_substrate,
    overdispersion_for_beta,
)
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG, TS_NEG, TS_NONE, TS_POS
from rigel.calibration.splice_junction import region_splice_gdna_frac
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.strand_deconv import cleaned_gdna_count, deconv_regions, strand_deconvolve
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
SC = {TS_POS: "POS", TS_NEG: "NEG", TS_NONE: "NONE", TS_AMBIG: "AMBIG"}


def oracle_contained_unspliced(bam_path, starts, ids, R):
    """Per region: (gdna, mrna, nrna) counts of CONTAINED-UNSPLICED read1 fragments, from origins."""
    g, m, n = np.zeros(R), np.zeros(R), np.zeros(R)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1:
                continue
            if r.cigartuples and any(op == 3 for op, _ in r.cigartuples):
                continue  # spliced — not contained-unspliced
            ref = r.reference_name
            if ref not in starts:
                continue
            i = int(np.searchsorted(starts[ref], r.reference_start, side="right") - 1)
            if i < 0:
                continue
            rid = int(ids[ref][i])
            kind = parse_origin(r.query_name).kind
            (g if kind == "gdna" else n if kind == "nrna" else m)[rid] += 1.0
    return g, m, n


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
    n_worst = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    bam = f"{SUITE}/{cond}/sim_oracle.bam"
    idx = TranscriptIndex.load(f"{SUITE}/rigel_index")

    cfg = PipelineConfig()
    ccfg = cfg.calibration
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    st, sm, flm, buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size,
    )

    # ---- replicate calibrate's region path, capturing intermediates --------------------------
    sub = CalibrationSubstrate.from_payload(pl, ra)
    region_eff_len = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    boundary_eff_len = boundary_side_eff_length(fl.gdna_pmf, ra.region_size_bp)
    fl_mean = boundary_eff_length(fl.gdna_pmf)
    rna_fl_mean = boundary_eff_length(fl.rna_pmf)
    region_eff_len_rna = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)

    nd_raw = node_gdna_density(sub, ra, region_eff_len, fl_mean, need_count_variance=False)
    gdna_strand = fit_gdna_strand_from_substrate(
        sub, ra, nd_raw, boundary_eff_len, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.gdna_strand_prior_alpha_beta),
        prior_weight=ccfg.gdna_strand_prior_weight,
    )
    od_g = gdna_strand.gdna_strand_overdispersion
    rna_strand = fit_rna_strand_from_substrate(
        sub, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.rna_strand_prior_alpha_beta),
        prior_weight=ccfg.rna_strand_prior_weight,
    )
    od_r = rna_strand.rna_strand_overdispersion
    i0 = ccfg.gdna_strand_info_scale

    cont_split, left_split, right_split = strand_deconvolve(
        sub, ra, rna_sense_frac=kappa, gdna_strand_overdispersion=od_g,
        rna_strand_overdispersion=od_r, deconv_quantile=ccfg.gdna_deconv_quantile, n_grid=ccfg.n_grid,
    )

    def rawc(v):
        return v.n_unspliced_pos.astype(np.float64) + v.n_unspliced_neg.astype(np.float64)

    raw_c, raw_l, raw_r = rawc(sub.contained), rawc(sub.left), rawc(sub.right)
    cl_l = cleaned_gdna_count(left_split, raw_l, i0)
    cl_r = cleaned_gdna_count(right_split, raw_r, i0)
    # g_count: RAW contained own-density + CLEANED boundary imputation (the count-only spatial estimate).
    nd = node_gdna_density(
        sub, ra, region_eff_len, fl_mean, need_count_variance=False, gdna_counts=(raw_c, cl_l, cl_r),
    )
    region_count_frac, n_up = region_splice_gdna_frac(
        sub, ra, nd.count_gdna_frac, eff_gdna=fl_mean, eff_rna=rna_fl_mean,
        eff_gdna_region=region_eff_len, eff_rna_region=region_eff_len_rna,
        left_gdna_unspl=cl_l, right_gdna_unspl=cl_r,
    )
    # THE COMBINE: g = w·g_strand + (1−w)·g_count, w=I/(I+I₀), via deconv_regions (per-node weight).
    region_node_density = dataclasses.replace(nd, count_gdna_frac=region_count_frac)
    regions = deconv_regions(
        sub, ra, region_node_density, rna_sense_frac=kappa, gdna_strand_overdispersion=od_g,
        rna_strand_overdispersion=od_r, deconv_quantile=ccfg.gdna_deconv_quantile,
        n_grid=ccfg.n_grid, info_scale=i0,
    )
    region_frac = np.asarray(regions.gdna_frac)  # the combined fraction
    calib_gdna = np.asarray(regions.gdna_mass)

    # ---- oracle truth per region -------------------------------------------------------------
    df = idx.region_df
    starts = {ref: gg["start"].to_numpy() for ref, gg in df.groupby("ref_name")}
    ids = {ref: gg["region_id"].to_numpy() for ref, gg in df.groupby("ref_name")}
    R = len(df)
    og, om, on = oracle_contained_unspliced(bam, starts, ids, R)
    o_tot = og + om + on
    o_gfrac = np.where(o_tot > 0, og / np.maximum(o_tot, 1e-9), np.nan)

    ts = np.asarray(ra.strand_class)
    sig = np.asarray(ra.signature)
    ref_name = df["ref_name"].to_numpy()
    start = df["start"].to_numpy()
    end = df["end"].to_numpy()
    w = np.where(cont_split.info > 0, cont_split.info / (cont_split.info + i0), 0.0)
    reg_obs, _ = nd.region_count_observable, nd.boundary_count_observable

    print(f"\n===== {cond} =====")
    print(f"kappa={kappa:.4f}  od_gdna={od_g:.4f}  od_rna={od_r:.4f}  i0={i0}  splice_upgraded={n_up}")
    print(f"TOTALS: oracle contained-unspl gDNA={og.sum():.0f} mRNA={om.sum():.0f}  "
          f"calib gDNA mass={calib_gdna.sum():.0f}  deficit={og.sum()-calib_gdna.sum():.0f}")

    deficit = og - calib_gdna
    print("  deficit by region class (oracle_gdna - calib_gdna):")
    for cls_val, cls_name in [(TS_AMBIG, "AMBIG"), (TS_POS, "POS"), (TS_NEG, "NEG"), (TS_NONE, "NONE")]:
        m = ts == cls_val
        tot = deficit.sum()
        print(f"    {cls_name:>5}: n={int(m.sum()):>5}  oracle_gdna={og[m].sum():>9.0f}  "
              f"calib_gdna={calib_gdna[m].sum():>9.0f}  deficit={deficit[m].sum():>+9.0f}  "
              f"({100*deficit[m].sum()/tot:>4.0f}% of total)")

    # worst regions by gDNA mass under-call (oracle gDNA - calib gDNA); optional class filter (argv[3])
    cls_filter = sys.argv[3] if len(sys.argv) > 3 else None
    order = [i for i in np.argsort(-deficit) if cls_filter is None or SC[ts[i]] == cls_filter]
    print(f"\n--- {n_worst} {cls_filter or 'all'} regions, largest gDNA UNDER-call (oracle_gdna - calib_gdna) ---")
    print("  (combine: reg_gf = w·gstr + (1−w)·gcount)")
    print(f"{'ref:start-end':>20}{'sig':>4}{'cls':>5}{'robs':>5}{'o_gdna':>7}{'o_mrna':>7}{'o_gfr':>6}"
          f"{'rawc':>7}{'gstr':>6}{'info':>8}{'w':>6}{'gcount':>7}{'reg_gf':>7}{'cgdna':>7}{'defic':>7}")
    for i in order[:n_worst]:
        if deficit[i] <= 0:
            break
        print(f"{f'{ref_name[i]}:{start[i]}-{end[i]}':>20}{int(sig[i]):>4}{SC[ts[i]]:>5}"
              f"{int(reg_obs[i]):>5}{og[i]:>7.0f}{om[i]:>7.0f}{o_gfrac[i]:>6.2f}"
              f"{raw_c[i]:>7.0f}{cont_split.gdna_frac[i]:>6.2f}{cont_split.info[i]:>8.0f}{w[i]:>6.2f}"
              f"{region_count_frac[i]:>7.2f}{region_frac[i]:>7.2f}{calib_gdna[i]:>7.0f}{deficit[i]:>7.0f}")


if __name__ == "__main__":
    main()
