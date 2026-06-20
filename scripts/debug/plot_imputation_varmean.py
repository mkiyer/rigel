"""Plot the gDNA + RNA imputation var~mean fits (Step 2 diagnostic).

Runs the REAL calibrate() (Step 2b: gDNA imputation shipped) on a scenario rich in single-strand multi-exon
genes (the fit set) + overlapping opposite-strand pairs (AMBIG nodes — where the imputation prior is APPLIED
but never TRAINED), with gDNA + nascent RNA present. It then reconstructs, at the converged f_g, the two
Poisson-offset imputation var~mean fits exactly as the production gDNA path (imputation.gdna_imputation_prior)
and the (reverted) per-strand RNA path build them, and plots:

  • the fit points (μ = the queried density, raw_var = the single-predictor residual (ρ_dst−ρ_src)²),
  • the learned σ²_bio(μ) curve (the EXCESS over the Poisson floor — what the prior precision uses),
  • the per-point Poisson offset V_p = ρ/L (the computed sampling floor that is subtracted),
  • a rug of the AMBIG region densities (where the RNA prior is applied — to expose train/apply mismatch).

The instructive contrast: σ²_bio is fit on SMOOTH within-transcript / boundary→region pairs (low dispersion),
so applying it at AMBIG (a different, higher-dispersion context) yields an OVER-CONFIDENT prior — the 2c
failure (complex-loci AMBIG +44%).

Usage:  python scripts/debug/plot_imputation_varmean.py [out.png]   (default /tmp/imputation_varmean.png)
"""
from __future__ import annotations

import dataclasses
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rigel.calibration.calibrate import calibrate
from rigel.calibration.density_model import node_gdna_density
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
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.run_fill import same_ref_left_right
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    TS_AMBIG,
    TS_NEG,
    TS_POS,
)
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.strand_deconv import cleaned_gdna_count, deconv_sides, strand_deconvolve
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.variance_model import MonotoneVarMean, pair_imputation_points
from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
from rigel.splice import SpliceType

_EPS = 1.0e-9


def build_scenario(work_dir):
    """Multi-exon alternating +/− genes (single-strand fit set) + overlapping pairs (AMBIG) + gDNA + nRNA."""
    win = 14000
    n_single = 10
    glen = (n_single + 8) * win
    sc = Scenario("imp_varmean", genome_length=glen, seed=11, work_dir=str(work_dir))
    rng = np.random.default_rng(3)
    for gi in range(n_single):
        base = (gi + 1) * win
        strand = "+" if gi % 2 == 0 else "-"
        exons = [(base + 1000, base + 2000), (base + 4000, base + 5000),
                 (base + 7000, base + 8000), (base + 10000, base + 11000)]
        sc.add_gene(f"g{gi}", strand, [{"t_id": f"g{gi}_t", "exons": exons,
                                        "abundance": int(rng.integers(120, 260))}])
    # overlapping opposite-strand pairs → AMBIG nodes
    for k in range(3):
        base = (n_single + 1 + k) * win
        sc.add_gene(f"a{k}+", "+", [{"t_id": f"a{k}p", "exons": [(base + 1000, base + 1500),
                                     (base + 4000, base + 6000)], "abundance": 140}])
        sc.add_gene(f"a{k}-", "-", [{"t_id": f"a{k}n", "exons": [(base + 5000, base + 7000),
                                     (base + 10000, base + 10500)], "abundance": 140}])
    res = sc.build_oracle(
        n_fragments=max(80000, (n_single + 6) * 4000),
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=100, strand_specificity=0.99, seed=11),
        gdna_config=GDNAConfig(abundance=120, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000),
        nrna_abundance=25.0,
    )
    return sc, res


def main(out_png="/tmp/imputation_varmean.png"):
    with tempfile.TemporaryDirectory() as td:
        sc, res = build_scenario(Path(td))
        idx, bam = res.index, str(res.bam_path)
        cfg = PipelineConfig()
        ccfg = cfg.calibration
        scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
        _st, sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
        ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
        fl = build_fl_models(global_counts=flm.global_model.counts,
                             rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                             gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
        result = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, ccfg)
        sub = CalibrationSubstrate.from_payload(pl, ra)
        sc.cleanup()

    # ---- geometry + converged f_g ----------------------------------------------------
    eff_len = np.maximum(region_eff_length(ra.region_size_bp, fl.gdna_pmf), _EPS)
    boundary_eff_len = boundary_side_eff_length(fl.gdna_pmf, ra.region_size_bp)
    fl_mean = boundary_eff_length(fl.gdna_pmf)
    eff_rna = np.maximum(region_eff_length(ra.region_size_bp, fl.rna_pmf), _EPS)
    mass_unspl = np.asarray(sub.contained.mass_unspliced, dtype=np.float64)
    mass_u = np.maximum(mass_unspl, _EPS)
    g_mass = np.asarray(result.mass_gdna_contained, dtype=np.float64)
    f_g = np.where(mass_unspl > _EPS, g_mass / mass_u, 0.0)
    sig = np.asarray(ra.signature).astype(np.int64)
    ts = np.asarray(ra.strand_class)
    kappa = float(result.rna_sense_frac)
    od_g, od_r = result.gdna_strand_overdispersion, result.rna_strand_overdispersion

    # ---- replicate calibrate's deconv_sides anchors (the gDNA predictors) -------------
    nd_raw = node_gdna_density(sub, ra, eff_len, fl_mean)
    gd = fit_gdna_strand_from_substrate(sub, ra, nd_raw, boundary_eff_len, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.gdna_strand_prior_alpha_beta),
        prior_weight=ccfg.gdna_strand_prior_weight).gdna_strand_overdispersion
    rd_ = fit_rna_strand_from_substrate(sub, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.rna_strand_prior_alpha_beta),
        prior_weight=ccfg.rna_strand_prior_weight).rna_strand_overdispersion
    _, lsplit, rsplit = strand_deconvolve(sub, ra, rna_sense_frac=kappa,
        gdna_strand_overdispersion=gd, rna_strand_overdispersion=rd_,
        deconv_quantile=ccfg.gdna_deconv_quantile, n_grid=ccfg.n_grid)
    i0 = ccfg.gdna_strand_info_scale

    def _raw(v):
        return v.n_unspliced_pos.astype(np.float64) + v.n_unspliced_neg.astype(np.float64)

    cl = cleaned_gdna_count(lsplit, _raw(sub.left), i0)
    cr = cleaned_gdna_count(rsplit, _raw(sub.right), i0)
    node_density = node_gdna_density(sub, ra, eff_len, fl_mean,
                                     gdna_counts=(_raw(sub.contained), cl, cr))
    left, right = deconv_sides(sub, ra, node_density, boundary_eff_len, rna_sense_frac=kappa,
        gdna_strand_overdispersion=gd, rna_strand_overdispersion=rd_,
        deconv_quantile=ccfg.gdna_deconv_quantile, n_grid=ccfg.n_grid, info_scale=i0)

    # ---- gDNA imputation fit inputs (mirror imputation.gdna_imputation_prior) ----------
    bco = np.asarray(node_density.boundary_count_observable, dtype=bool)
    obs = np.asarray(node_density.region_count_observable, dtype=bool)
    ls_same, rs_same = same_ref_left_right(np.asarray(ra.ref_id))
    n = ra.n_regions
    left_obs = np.zeros(n, bool); right_obs = np.zeros(n, bool)
    if n > 1:
        left_obs[1:] = bco[:-1] & ls_same[1:]
        right_obs[:-1] = bco[:-1] & rs_same[:-1]
    inv_side = 1.0 / np.maximum(boundary_eff_len, _EPS)
    g_left = np.asarray(left.gdna_mass, float); g_right = np.asarray(right.gdna_mass, float)
    d_left = g_left * inv_side; d_right = g_right * inv_side
    rho_g = g_mass / eff_len
    g_means, g_raw, g_off = pair_imputation_points(
        rho_g, d_left, d_right,
        region_eligible=~obs, left_ok=left_obs & (g_left > 0.0), right_ok=right_obs & (g_right > 0.0),
        region_var_samp=rho_g / eff_len, left_var_samp=d_left * inv_side, right_var_samp=d_right * inv_side,
    )
    g_curve = MonotoneVarMean.fit_offset(g_means, g_raw, g_off) if g_means.size else None

    # ---- RNA imputation fit inputs (mirror the reverted per-strand 2c builder) ---------
    rna_density = (1.0 - f_g) * mass_unspl / eff_rna
    var_samp_r = rna_density / eff_rna
    has = {TS_POS: (sig & (BIT_EXON_POS | BIT_INTRON_POS)) != 0,
           TS_NEG: (sig & (BIT_EXON_NEG | BIT_INTRON_NEG)) != 0}

    def _l(a):
        o = np.full(n, np.nan); o[1:] = a[:-1]; return o

    def _r(a):
        o = np.full(n, np.nan); o[:-1] = a[1:]; return o

    rho_l, rho_r = _l(rna_density), _r(rna_density)
    vs_l, vs_r = _l(var_samp_r), _r(var_samp_r)
    r_m, r_raw, r_off = [], [], []
    for s in (TS_POS, TS_NEG):
        pred = ts == s
        pl_ = np.zeros(n, bool); pr_ = np.zeros(n, bool)
        if n > 1:
            pl_[1:] = pred[:-1] & ls_same[1:]; pr_[:-1] = pred[1:] & rs_same[:-1]
        lok = pl_ & has[s] & (np.nan_to_num(rho_l) > 0.0)
        rok = pr_ & has[s] & (np.nan_to_num(rho_r) > 0.0)
        m, rw, off = pair_imputation_points(rna_density, rho_l, rho_r,
            region_eligible=(ts == s), left_ok=lok, right_ok=rok,
            region_var_samp=var_samp_r, left_var_samp=vs_l, right_var_samp=vs_r)
        r_m.append(m); r_raw.append(rw); r_off.append(off)
    rm = np.concatenate(r_m); rr = np.concatenate(r_raw)
    ro = np.concatenate([o for o in r_off if o is not None]) if any(o is not None for o in r_off) else None
    r_curve = MonotoneVarMean.fit_offset(rm, rr, ro) if (ro is not None and rm.size) else None
    ambig_rna_dens = rna_density[(ts == TS_AMBIG) & (mass_unspl > _EPS)]

    # ---- plot --------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for ax, (title, means, raw, off, curve, rug) in zip(axes, [
        ("gDNA imputation var~mean (boundary-side → exon/AMBIG region)", g_means, g_raw, g_off, g_curve, None),
        ("RNA imputation var~mean (single-strand neighbour → region, ± pooled)", rm, rr, ro, r_curve,
         ambig_rna_dens),
    ]):
        m = np.asarray(means); rw = np.asarray(raw); of = np.asarray(off)
        ok = (m > _EPS) & (rw > 0.0)
        ax.scatter(m[ok], rw[ok], s=10, alpha=0.35, color="0.5",
                   label="raw residual² = σ²_bio + V_p (fit points)")
        ofk = (m > _EPS) & (of > 0.0)
        ax.scatter(m[ofk], of[ofk], s=8, alpha=0.35, color="tab:orange",
                   label="Poisson offset V_p = ρ/L (subtracted floor)")
        if curve is not None and m[ok].size:
            grid = np.logspace(np.log10(m[ok].min()), np.log10(m[ok].max()), 200)
            ax.plot(grid, np.maximum(curve.predict(grid), 1e-12), color="tab:blue", lw=2.2,
                    label="learned σ²_bio(μ) (the prior precision uses 1/(σ²_bio+ρ_src/L_src))")
        if rug is not None and rug.size:
            ax.scatter(rug, np.full(rug.size, ax.get_ylim()[0] if False else 1e-9), marker="|",
                       color="tab:red", s=200, label=f"AMBIG region densities (apply, N={rug.size})")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("μ  (queried density)"); ax.set_ylabel("variance")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, loc="upper left"); ax.grid(True, which="both", alpha=0.15)
    fig.suptitle("Imputation var~mean fits — σ²_bio (learned excess) vs the Poisson floor V_p", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    print(f"wrote {out_png}")
    print(f"  gDNA fit: {g_means.size} points, μ-range [{g_means.min():.4g}, {g_means.max():.4g}]"
          if g_means.size else "  gDNA fit: 0 points")
    print(f"  RNA  fit: {rm.size} points, μ-range [{rm.min():.4g}, {rm.max():.4g}]"
          if rm.size else "  RNA fit: 0 points")
    print(f"  AMBIG apply nodes: {ambig_rna_dens.size}, RNA-density range "
          f"[{ambig_rna_dens.min():.4g}, {ambig_rna_dens.max():.4g}]" if ambig_rna_dens.size else
          "  AMBIG apply nodes: 0")
    if g_curve is not None and r_curve is not None and ambig_rna_dens.size:
        print(f"  σ²_bio_gDNA at median μ: {float(g_curve.predict(np.array([np.median(g_means)]))[0]):.4g}")
        print(f"  σ²_bio_RNA  at median fit μ: {float(r_curve.predict(np.array([np.median(rm)]))[0]):.4g}")
        print(f"  σ²_bio_RNA  at AMBIG μ (extrapolated apply): "
              f"{float(r_curve.predict(np.array([np.median(ambig_rna_dens)]))[0]):.4g}")
    return out_png


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/tmp/imputation_varmean.png")
