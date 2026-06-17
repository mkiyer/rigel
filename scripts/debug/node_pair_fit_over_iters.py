"""Diagnostic — the gDNA (and RNA) node-PAIR `var~mean` IMPUTATION fit OVER calibrate iterations.

Replicates the production calibrate() iterative loop (src/rigel/calibration/calibrate.py:268-350)
VERBATIM — importing the SAME internals, NOT modifying production code — on a capture-ON, ss=0.99,
+nascent scenario (the data-rich case; reuses the `capon` builder from phase1_m2_calibrate_internals).

At EACH PASS it captures the gDNA node-PAIR IMPUTATION `MonotoneVarMean` (the eff-len-FIXED model: one
fit point per region<->boundary adjacency, mean = region gDNA density, raw_var = (d_region - d_side)²,
side density = side gDNA MASS / E[min(ℓ,L_side)]). `.to_dataframe()` → the individual (mean, var) fit
POINTS + the sampled monotone spline CURVE. Same for the RNA node-pair fit (standalone/inert in the
solve, but we fit it per pass on the current gDNA estimate).

Outputs (PNG ≥150 dpi + per-pass feather/TSV fallbacks) to ~/Downloads/rigel_runs/node_pair_fit_plots/:
  1. gdna_imputation_fit_by_pass.png       — small-multiples log-log scatter(points)+line(curve) per pass
  2. gdna_imputation_error_convergence.png — per-pass median/p90 raw_var + mean|f_g-prev| over passes
  3. rna_imputation_fit_by_pass.png        — same as (1) for the RNA node-pair fit

    python scripts/debug/node_pair_fit_over_iters.py
"""
import dataclasses
import os
import pathlib
import tempfile

import numpy as np
import pandas as pd

from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
from rigel.sim.capture import CaptureConfig
from rigel.sim.capture.design import write_random_capture_probes
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.config import PipelineConfig
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.splice import SpliceType

# the SAME production internals calibrate() uses (verbatim replication; no prod edits)
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.density_model import node_gdna_density, count_observable_masks
from rigel.calibration.effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate,
    fit_rna_strand_from_substrate,
    overdispersion_for_beta,
)
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.strand_deconv import cleaned_gdna_count, deconv_sides, strand_deconvolve
from rigel.calibration.splice_junction import region_splice_gdna_frac
from rigel.calibration.variance_model import (
    direct_points,
    fit_direct_varmean,
    fit_pair_imputation_varmean,
    fit_pair_imputation_rna_varmean,
)
from rigel.calibration.simplex_sweep import deconv_regions_sweep
from rigel.calibration.run_fill import same_ref_left_right

OUT = os.path.expanduser("~/Downloads/rigel_runs/node_pair_fit_plots")
SCENARIO = "capon (capture-ON, ss=0.99, +nascent)"


# --- scenario builder (copied VERBATIM from scripts/debug/phase1_m2_calibrate_internals.py; that
#     module's top-level imports reference renamed-away symbols in this eff-len-fixed tree, so we inline
#     its `capon` builder rather than import it). ---
WIN = 12000
N_GENES = 24


def _build_genes(sc):
    rng = np.random.default_rng(7)
    for gi in range(N_GENES):
        base = (gi + 1) * WIN
        strand = "+" if gi % 2 == 0 else "-"
        e1 = (base + 1000, base + 2000)
        e2 = (base + 4000, base + 5500)
        e3 = (base + 8000, base + 9000)
        ab = int(rng.integers(60, 200))
        sc.add_gene(f"g{gi}", strand, [{"t_id": f"g{gi}_t", "exons": [e1, e2, e3], "abundance": ab}])
    sbase = (N_GENES + 2) * WIN
    sc.add_gene("seedP", "+", [{"t_id": "SP", "exons": [(sbase, sbase + 500),
                (sbase + 1500, sbase + 2000), (sbase + 3000, sbase + 3500)], "abundance": 180}])
    sc.add_gene("seedM", "-", [{"t_id": "SM", "exons": [(sbase + 5000, sbase + 5500),
                (sbase + 6500, sbase + 7000), (sbase + 8000, sbase + 8500)], "abundance": 180}])


def build_capon(work):
    """The capon (capture-ON, ss=0.99, +nascent) scenario. Returns (sc, pl, ra, sm, gpmf, rpmf, cfg, R)."""
    glen = (N_GENES + 4) * WIN
    sc = Scenario("node_pair_capon", genome_length=glen, seed=13, work_dir=work)
    _build_genes(sc)
    sim_cfg = ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                            read_length=100, strand_specificity=0.99, seed=13)
    transcripts = sc.annotation.get_transcripts()
    probes = work + "/probes.tsv"
    write_random_capture_probes(transcripts, pathlib.Path(probes),
                                capture_fraction=0.6, probe_length=120, probe_density=1.0, seed=7)
    capture_cfg = CaptureConfig(probes=probes, off_target_weight=1.0, binding_per_base=10.0,
                                gdna_split_penalty=0.2)
    gdna_cfg = GDNAConfig(abundance=120.0, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)
    res = sc.build_oracle(n_fragments=max(120000, N_GENES * 4000), sim_config=sim_cfg,
                          gdna_config=gdna_cfg, capture_config=capture_cfg, nrna_abundance=25.0)
    idx = res.index
    bam = str(res.bam_path)
    scan = dataclasses.replace(PipelineConfig().scan, sj_strand_tag=_native_detect_sj_tag(bam))
    st, sm, flm, buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    cfg = PipelineConfig().calibration
    return sc, pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg, len(idx.region_df)


def run(pl, ra, sm, gdna_fl_pmf, rna_fl_pmf, cfg):
    """VERBATIM replication of calibrate() setup + loop, capturing per-pass node-pair fits."""
    substrate = CalibrationSubstrate.from_payload(pl, ra)

    # --- setup (calibrate.py:99-215) ---
    region_eff_len = region_eff_length(ra.region_size_bp, gdna_fl_pmf)
    boundary_eff_len = boundary_side_eff_length(gdna_fl_pmf, ra.region_size_bp)
    fl_mean = boundary_eff_length(gdna_fl_pmf)
    rna_fl_mean = boundary_eff_length(rna_fl_pmf)
    region_eff_len_rna = region_eff_length(ra.region_size_bp, rna_fl_pmf)

    rna_sense_frac = float(fit_strand_balance(sm).rna_sense_frac)
    need_count_variance = float(cfg.gdna_deconv_quantile) != 0.5
    node_density_raw = node_gdna_density(
        substrate, ra, region_eff_len, fl_mean, need_count_variance=need_count_variance
    )
    gdna_strand = fit_gdna_strand_from_substrate(
        substrate, ra, node_density_raw, boundary_eff_len, rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(cfg.gdna_strand_prior_alpha_beta),
        prior_weight=cfg.gdna_strand_prior_weight)
    gdna_od = gdna_strand.gdna_strand_overdispersion
    rna_strand = fit_rna_strand_from_substrate(
        substrate, rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(cfg.rna_strand_prior_alpha_beta),
        prior_weight=cfg.rna_strand_prior_weight)
    rna_od = rna_strand.rna_strand_overdispersion

    _, left_split, right_split = strand_deconvolve(
        substrate, ra, rna_sense_frac=rna_sense_frac, gdna_strand_overdispersion=gdna_od,
        rna_strand_overdispersion=rna_od, deconv_quantile=cfg.gdna_deconv_quantile, n_grid=cfg.n_grid)
    i0 = cfg.gdna_strand_info_scale

    def _raw(view):
        return view.n_unspliced_pos.astype(np.float64) + view.n_unspliced_neg.astype(np.float64)

    cleaned_left = cleaned_gdna_count(left_split, _raw(substrate.left), i0)
    cleaned_right = cleaned_gdna_count(right_split, _raw(substrate.right), i0)

    node_density = node_gdna_density(
        substrate, ra, region_eff_len, fl_mean, need_count_variance=need_count_variance,
        gdna_counts=(_raw(substrate.contained), cleaned_left, cleaned_right))
    region_count_frac, _ = region_splice_gdna_frac(
        substrate, ra, node_density.count_gdna_frac, eff_gdna=fl_mean, eff_rna=rna_fl_mean,
        eff_gdna_region=region_eff_len, eff_rna_region=region_eff_len_rna,
        left_gdna_unspl=cleaned_left, right_gdna_unspl=cleaned_right)

    left, right = deconv_sides(
        substrate, ra, node_density, boundary_eff_len, rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_od, rna_strand_overdispersion=rna_od,
        deconv_quantile=cfg.gdna_deconv_quantile, n_grid=cfg.n_grid, info_scale=i0)

    # --- loop precomputation (calibrate.py:238-265) ---
    c = substrate.contained
    u_tot = c.n_unspliced_pos.astype(np.float64) + c.n_unspliced_neg.astype(np.float64)
    eff_len = np.maximum(np.asarray(region_eff_len, dtype=np.float64), 1e-9)
    mass_u = np.maximum(np.asarray(c.mass_unspliced, dtype=np.float64), 1e-9)
    geom2 = (eff_len / mass_u) ** 2
    obs = np.asarray(node_density.region_count_observable, dtype=bool)
    gdna_left = np.asarray(left.gdna_mass, dtype=np.float64)
    gdna_right = np.asarray(right.gdna_mass, dtype=np.float64)
    bco = np.asarray(node_density.boundary_count_observable, dtype=bool)
    ls_same, rs_same = same_ref_left_right(np.asarray(ra.ref_id))
    left_obs = np.zeros(ra.n_regions, dtype=bool)
    right_obs = np.zeros(ra.n_regions, dtype=bool)
    if ra.n_regions > 1:
        left_obs[1:] = bco[:-1] & ls_same[1:]
        right_obs[:-1] = bco[:-1] & rs_same[:-1]
    region_eligible_g = ~obs
    inv_side_len_g = 1.0 / np.maximum(boundary_eff_len, 1e-9)
    gdna_c = u_tot.copy()  # Pass-0 all-gDNA init
    prev_fg = None

    gdna_fits = []   # per-pass MonotoneVarMean (gDNA imputation)
    rna_fits = []    # per-pass MonotoneVarMean (RNA imputation) or None
    summaries = []   # per-pass dict

    for _pass in range(int(cfg.sweep_max_passes)):
        rho_global = (
            float(gdna_c[obs].sum() / max(float(eff_len[obs].sum()), 1e-9)) if obs.any() else 0.0
        )
        direct = fit_direct_varmean(
            direct_points(substrate, ra, region_eff_len, boundary_eff_len,
                          gdna_views=(gdna_c, gdna_left, gdna_right)))
        region_density_g = gdna_c / eff_len
        # ---- THE gDNA NODE-PAIR IMPUTATION FIT (capture this) ----
        imputation = fit_pair_imputation_varmean(
            region_density_g, gdna_left * inv_side_len_g, gdna_right * inv_side_len_g,
            region_eligible=region_eligible_g, left_ok=left_obs & (gdna_left > 0.0),
            right_ok=right_obs & (gdna_right > 0.0), ref_id=ra.ref_id)
        gdna_fits.append(imputation)

        # ---- the RNA node-pair imputation fit on the current gDNA estimate (standalone/inert) ----
        fg_now = np.where(mass_u > 1e-12, gdna_c / mass_u, 0.0)
        fg_now = np.clip(fg_now, 0.0, 1.0)
        rna_boundary_side_eff_len = boundary_side_eff_length(rna_fl_pmf, ra.region_size_bp)
        try:
            rna_fit = fit_pair_imputation_rna_varmean(
                substrate, ra, region_eff_len_rna, rna_boundary_side_eff_len,
                gdna_frac=fg_now, left_gdna_frac=left_split.gdna_frac,
                right_gdna_frac=right_split.gdna_frac,
                cleaned_left=cleaned_left, cleaned_right=cleaned_right)
        except Exception as e:  # noqa: BLE001
            print(f"  [pass {_pass}] RNA fit raised {type(e).__name__}: {e}")
            rna_fit = None
        rna_fits.append(rna_fit)

        mu = gdna_c / eff_len
        var_d = np.where(obs, direct.predict(mu), imputation.predict(mu))
        cf = np.clip(region_count_frac, 0.0, 1.0)
        sig2_frac = np.maximum(np.minimum(var_d * geom2, cf * (1.0 - cf)), 1e-12)
        active_mu = mu[u_tot > 0.0]
        sigma_d_global = (1.4826 * float(np.median(np.abs(active_mu - np.median(active_mu))))
                          if active_mu.size else rho_global)
        sig2_glob = np.maximum(sigma_d_global ** 2 * geom2, 1e-12)
        tau_count = np.minimum(1.0 / sig2_frac, mass_u)
        tau_global = np.clip(1.0 / sig2_glob, 1.0, mass_u)
        regions = deconv_regions_sweep(
            substrate, ra, rna_sense_frac=rna_sense_frac, gdna_strand_overdispersion=gdna_od,
            rna_strand_overdispersion=rna_od, count_gdna_frac=region_count_frac,
            count_precision=tau_count, n_grid=cfg.sweep_n_grid, rho_global=rho_global,
            region_eff_len=region_eff_len, info_scale=i0, global_tau=tau_global)

        fg = np.asarray(regions.gdna_frac, dtype=np.float64)
        delta = float(np.mean(np.abs(fg - prev_fg))) if prev_fg is not None else float("nan")

        # per-pass residual summary on the fit POINTS
        gpts = imputation.fit_mean, imputation.fit_var
        n_g = int(gpts[0].size)
        med_g = float(np.median(gpts[1])) if n_g else float("nan")
        p90_g = float(np.percentile(gpts[1], 90)) if n_g else float("nan")
        rmu = imputation.fit_mean
        summaries.append({
            "pass": _pass, "rho_global": rho_global,
            "n_points_gdna": n_g, "median_raw_var_gdna": med_g, "p90_raw_var_gdna": p90_g,
            "edf_gdna": float(imputation.edf), "lam_gdna": float(imputation.lam),
            "mu_lo_gdna": float(rmu.min()) if n_g else float("nan"),
            "mu_hi_gdna": float(rmu.max()) if n_g else float("nan"),
            "fg_delta": delta, "mean_fg": float(np.mean(fg)),
            "n_points_rna": int(rna_fit.fit_mean.size) if rna_fit is not None else 0,
            "median_raw_var_rna": (float(np.median(rna_fit.fit_var))
                                   if rna_fit is not None and rna_fit.fit_var.size else float("nan")),
            "edf_rna": float(rna_fit.edf) if rna_fit is not None else float("nan"),
        })
        print(f"  PASS {_pass}: rho_global={rho_global:.5g}  gDNA pts={n_g}  "
              f"med_raw_var={med_g:.4g}  p90={p90_g:.4g}  edf={imputation.edf:.2f}  "
              f"mu∈[{summaries[-1]['mu_lo_gdna']:.4g},{summaries[-1]['mu_hi_gdna']:.4g}]  "
              f"|f_g-prev|={delta:.3e}  RNA pts={summaries[-1]['n_points_rna']}")

        if prev_fg is not None and delta < cfg.sweep_convergence_delta:
            print(f"    -> CONVERGED (delta < {cfg.sweep_convergence_delta}) at pass {_pass}; "
                  "production loop breaks here")
            break
        prev_fg = fg
        gdna_c = np.asarray(regions.gdna_mass, dtype=np.float64)

    return gdna_fits, rna_fits, pd.DataFrame(summaries)


def _plot_fits_by_pass(fits, title, out_png, kind_label):
    """Small-multiples grid: one panel per pass, log-log scatter(points)+line(curve)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    have = [(i, f) for i, f in enumerate(fits) if f is not None and f.fit_mean.size > 0]
    if not have:
        print(f"  [{kind_label}] no passes with >=1 fit point — writing a note panel.")
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, f"{kind_label}: no fit points in any pass", ha="center", va="center")
        ax.set_axis_off()
        fig.suptitle(title)
        fig.savefig(out_png, dpi=160, bbox_inches="tight")
        plt.close(fig)
        return 0

    n = len(have)
    ncol = min(n, 3)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 4.2 * nrow), squeeze=False)
    for ax in axes.flat:
        ax.set_visible(False)
    for panel, (pi, f) in enumerate(have):
        ax = axes.flat[panel]
        ax.set_visible(True)
        df = f.to_dataframe()
        pts = df[df["kind"] == "point"]
        crv = df[df["kind"] == "curve"]
        ax.scatter(pts["mean"], pts["var"], s=14, alpha=0.45, color="C0",
                   edgecolors="none", label="fit points")
        if not crv.empty:
            ax.plot(crv["mean"], crv["var"], color="C3", lw=2.0, label="monotone spline")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("node gDNA density  μ = mass / E[max(0,L−ℓ)]" if "gDNA" in kind_label
                      else "node RNA density  μ")
        ax.set_ylabel("imputation residual var  raw_var = (d_region − d_side)²")
        med = float(np.median(pts["var"])) if len(pts) else float("nan")
        ax.set_title(f"pass {pi}   n={len(pts)}   median raw_var={med:.3g}\n"
                     f"edf={f.edf:.2f}  lam={f.lam:.3g}", fontsize=9)
        ax.legend(fontsize=7, loc="best")
        ax.grid(True, which="both", ls=":", alpha=0.3)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_png, dpi=160)
    plt.close(fig)
    return n


def _plot_convergence(summary, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax1 = plt.subplots(figsize=(8, 5))
    p = summary["pass"].to_numpy()
    ax1.plot(p, summary["median_raw_var_gdna"], "o-", color="C0",
             label="gDNA median raw_var (fit points)")
    ax1.plot(p, summary["p90_raw_var_gdna"], "s--", color="C0", alpha=0.6,
             label="gDNA p90 raw_var")
    if summary["median_raw_var_rna"].notna().any():
        ax1.plot(p, summary["median_raw_var_rna"], "^-", color="C2",
                 label="RNA median raw_var")
    ax1.set_yscale("log")
    ax1.set_xlabel("calibrate pass (iteration)")
    ax1.set_ylabel("imputation residual variance  raw_var  (log)")
    ax1.set_xticks(p)
    ax1.grid(True, which="both", ls=":", alpha=0.3)

    ax2 = ax1.twinx()
    fgd = summary["fg_delta"].to_numpy()
    # pass-0 delta is NaN (no prev); fg_delta may be exactly 0 (the solve converged). Only log-scale
    # ax2 if there is at least one positive value to show; otherwise keep it linear and annotate.
    fgd_pos = np.where(np.isfinite(fgd) & (fgd > 0), fgd, np.nan)
    any_pos = bool(np.isfinite(fgd_pos).any())
    ax2.plot(p, fgd_pos, "d-", color="C3", label="mean|f_g − prev_f_g|")
    if any_pos:
        ax2.set_yscale("log")
        ax2.set_ylabel("solve convergence  mean|f_g − prev_f_g|  (log)", color="C3")
    else:
        ax2.set_ylim(-0.05, 1.0)
        ax2.set_ylabel("solve convergence  mean|f_g − prev_f_g|", color="C3")
    y0 = ax2.get_ylim()[0]
    for xi, raw in zip(p, fgd):
        if np.isfinite(raw):
            ax2.annotate(f"Δf_g={raw:.1e}", (xi, y0), color="C3", fontsize=8,
                         ha="center", va="bottom")
    ax2.tick_params(axis="y", labelcolor="C3")

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=8, loc="best")
    ax1.set_title(f"Node-pair imputation error vs iteration — {SCENARIO}", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def _dump_fit_tables(fits, prefix):
    for i, f in enumerate(fits):
        if f is None:
            continue
        df = f.to_dataframe()
        df.insert(0, "pass", i)
        df.to_csv(os.path.join(OUT, f"{prefix}_pass{i}.tsv"), sep="\t", index=False)
        try:
            df.to_feather(os.path.join(OUT, f"{prefix}_pass{i}.feather"))
        except Exception:
            pass


def main():
    os.makedirs(OUT, exist_ok=True)
    with tempfile.TemporaryDirectory() as work:
        sc, pl, ra, sm, gpmf, rpmf, cfg, R = build_capon(work)
        print(f"=== {SCENARIO}  (R={R} regions) ===")
        try:
            gdna_fits, rna_fits, summary = run(pl, ra, sm, gpmf, rpmf, cfg)
        finally:
            try:
                sc.cleanup()
            except Exception:
                pass

    summary.to_csv(os.path.join(OUT, "per_pass_summary.tsv"), sep="\t", index=False)
    try:
        summary.to_feather(os.path.join(OUT, "per_pass_summary.feather"))
    except Exception:
        pass
    _dump_fit_tables(gdna_fits, "gdna_fit")
    _dump_fit_tables(rna_fits, "rna_fit")

    g_png = os.path.join(OUT, "gdna_imputation_fit_by_pass.png")
    c_png = os.path.join(OUT, "gdna_imputation_error_convergence.png")
    r_png = os.path.join(OUT, "rna_imputation_fit_by_pass.png")
    n_g = _plot_fits_by_pass(
        gdna_fits, f"gDNA node-pair imputation var~mean by pass — {SCENARIO}", g_png, "gDNA")
    _plot_convergence(summary, c_png)
    n_r = _plot_fits_by_pass(
        rna_fits, f"RNA node-pair imputation var~mean by pass — {SCENARIO}", r_png, "RNA")

    readme = os.path.join(OUT, "README.txt")
    with open(readme, "w") as fh:
        fh.write(_readme_text(summary, n_g, n_r))

    print("\n=== SAVED ===")
    for pth in (g_png, c_png, r_png, readme,
                os.path.join(OUT, "per_pass_summary.tsv")):
        print("  " + pth)
    print("\n=== per-pass summary ===")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(summary.to_string(index=False))


def _readme_text(summary, n_g, n_r):
    s0 = summary.iloc[0]
    sf = summary.iloc[-1]
    mrv = summary["median_raw_var_gdna"].to_numpy()
    decreasing = bool(len(mrv) >= 2 and np.all(np.diff(mrv) <= 1e-12))
    verdict = ("YES — median gDNA imputation residual decreases monotonically across passes"
               if decreasing else
               "median gDNA imputation residual is NOT strictly monotone (see per_pass_summary.tsv)")
    return f"""Node-pair var~mean IMPUTATION fit over calibrate iterations
Scenario: {SCENARIO}
Generated by scripts/debug/node_pair_fit_over_iters.py (replicates calibrate.py:268-350 verbatim;
production code NOT modified).

THE MODEL (eff-len-FIXED node-pair imputation, CALIBRATION_PLAN_v5 §3):
  One fit point per region<->boundary-side ADJACENCY (one imputed-region dest + one observable boundary
  side as the single predictor). In density space (log-log):
    mean    = region gDNA density  = region gDNA mass / E[max(0, L − ℓ)]   (the queried axis)
    raw_var = (d_region − d_side)²  (the FULL single-predictor imputation residual, dof=1 Jensen)
  The side density = side gDNA MASS / E[min(ℓ, L_side)]  (the per-side DENSITY length, the eff-len fix;
  NOT mass / E[ℓ], which under-states short-flank density and fabricates a spurious residual). A region
  with both flanks eligible contributes TWO points; one eligible flank, ONE point. MonotoneVarMean =
  monotone-increasing P-spline (GCV-λ, IRLS-robust) in log-log, power-law fallback when <{18} points.

FILES:
  gdna_imputation_fit_by_pass.png       — small-multiples (one panel per pass): log-log scatter of the
                                          individual fit POINTS (mean=node density, y=raw_var) + the
                                          monotone spline CURVE. Title per panel: pass index, # points,
                                          median raw_var, edf, lam. {n_g} panel(s).
  gdna_imputation_error_convergence.png — headline convergence view: per-pass median & p90 raw_var
                                          (left log axis) + the solve's mean|f_g − prev_f_g| (right log
                                          axis) vs pass number.
  rna_imputation_fit_by_pass.png        — same as the first plot for the RNA node-pair fit
                                          (fit_pair_imputation_rna_varmean; standalone/inert in the solve,
                                          fit each pass on the current gDNA estimate). {n_r} panel(s).
  per_pass_summary.tsv / .feather       — per-pass scalars (rho_global, n_points, median/p90 raw_var,
                                          edf, lam, mu range, f_g delta, mean f_g; gDNA + RNA).
  gdna_fit_pass{{N}}.tsv / .feather       — per-pass to_dataframe(): the fit POINTS + the sampled CURVE.
  rna_fit_pass{{N}}.tsv / .feather        — same for the RNA fit.

HEADLINE (gDNA imputation):
  median raw_var  pass {int(s0['pass'])} -> pass {int(sf['pass'])}:  {s0['median_raw_var_gdna']:.4g} -> {sf['median_raw_var_gdna']:.4g}
  # points        pass {int(s0['pass'])} -> pass {int(sf['pass'])}:  {int(s0['n_points_gdna'])} -> {int(sf['n_points_gdna'])}
  fit mu-range (final pass):  [{sf['mu_lo_gdna']:.4g}, {sf['mu_hi_gdna']:.4g}]
  rho_global      pass {int(s0['pass'])} -> pass {int(sf['pass'])}:  {s0['rho_global']:.4g} -> {sf['rho_global']:.4g}

DOES THE ERROR DECREASE OVER ITERATIONS?  {verdict}
"""


if __name__ == "__main__":
    main()
