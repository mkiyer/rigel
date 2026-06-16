"""Phase-1 Measurement 2 — calibrate() per-pass internals, instrumented.

Replicates the calibrate() iterative loop (src/rigel/calibration/calibrate.py:90-297) line-for-line,
importing the SAME internals, and dumps per pass:

  (1) rho_global + whether log(rho_global) is inside the DIRECT fit range [x_lo,x_hi].
  (2) p10/p50/p90 over nodes of geom2, var_d, sig2_frac, tau_count; PLUS the realized count-prior
      influence: fraction of nodes where tau_count >~ I, tau_count << I, tau_count ~ 0
      (I = N*(2k-1)^2 = strand Fisher info).
  (3) convergence trace: per-pass mean|f_g - prev_f_g|.
  (4) wall-time: var~mean refit (varmean_points + fit_direct + fit_imputation) vs deconv_regions_sweep;
      node count R.
  (5) triplet / pair mu-coverage: count + mu min/median/max of (i) all-three-observable triplets and
      (ii) genomically-adjacent BOTH-count-observable region pairs, vs the mu-range of the EXON
      (non-count-observable) nodes.

Builds THREE scenarios (capture-ON+nascent, capture-OFF+nascent, zero-DNA) — multi-exon genes, ss=0.99.

    python scripts/debug/phase1_m2_calibrate_internals.py [scenario=all|capon|capoff|zerodna]
"""
import dataclasses
import sys
import tempfile
import time

import numpy as np

from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
from rigel.sim.capture import CaptureConfig
from rigel.sim.capture.design import write_random_capture_probes
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.config import PipelineConfig
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.splice import SpliceType

# the SAME internals calibrate() uses
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
    fit_direct_varmean,
    fit_imputation_varmean,
    varmean_points,
)
from rigel.calibration.simplex_sweep import deconv_regions_sweep
from rigel.calibration.run_fill import same_ref_left_right


# --------------------------------------------------------------------------------------------------
# Scenario construction — a few dozen multi-exon genes on both strands, with strand-training seeds.
# --------------------------------------------------------------------------------------------------
WIN = 12000
N_GENES = 24  # a few dozen


def _build_genes(sc):
    """Multi-exon genes alternating strand, each in its own window; plus +/- strand-training seeds."""
    rng = np.random.default_rng(7)
    for gi in range(N_GENES):
        base = (gi + 1) * WIN
        strand = "+" if gi % 2 == 0 else "-"
        # 3-exon gene with introns (gives count-observable introns + boundaries)
        e1 = (base + 1000, base + 2000)
        e2 = (base + 4000, base + 5500)
        e3 = (base + 8000, base + 9000)
        ab = int(rng.integers(60, 200))
        sc.add_gene(f"g{gi}", strand, [{"t_id": f"g{gi}_t", "exons": [e1, e2, e3], "abundance": ab}])
    # strand-training seeds (multi-exon, both strands), in the trailing windows
    sbase = (N_GENES + 2) * WIN
    sc.add_gene("seedP", "+", [{"t_id": "SP", "exons": [(sbase, sbase + 500),
                (sbase + 1500, sbase + 2000), (sbase + 3000, sbase + 3500)], "abundance": 180}])
    sc.add_gene("seedM", "-", [{"t_id": "SM", "exons": [(sbase + 5000, sbase + 5500),
                (sbase + 6500, sbase + 7000), (sbase + 8000, sbase + 8500)], "abundance": 180}])


def build_scenario(kind, work):
    """kind in {capon, capoff, zerodna}. Returns (pl, ra, sm, gdna_pmf, rna_pmf, cfg, R)."""
    glen = (N_GENES + 4) * WIN
    sc = Scenario(f"phase1_{kind}", genome_length=glen, seed=13, work_dir=work)
    _build_genes(sc)

    sim_cfg = ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                            read_length=100, strand_specificity=0.99, seed=13)

    capture_cfg = None
    if kind == "capon":
        # design probes on ~half the transcripts → hybrid-capture enrichment
        transcripts = sc.annotation.get_transcripts()
        probes = work + "/probes.tsv"
        write_random_capture_probes(transcripts, __import__("pathlib").Path(probes),
                                    capture_fraction=0.6, probe_length=120, probe_density=1.0, seed=7)
        capture_cfg = CaptureConfig(probes=probes, off_target_weight=1.0, binding_per_base=10.0,
                                    gdna_split_penalty=0.2)

    if kind == "zerodna":
        gdna_cfg = None
        nrna = 0.0
        gdna_present = False
    else:
        gdna_cfg = GDNAConfig(abundance=120.0, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)
        nrna = 25.0
        gdna_present = True

    res = sc.build_oracle(n_fragments=max(120000, N_GENES * 4000), sim_config=sim_cfg,
                          gdna_config=gdna_cfg, capture_config=capture_cfg, nrna_abundance=nrna)
    idx = res.index
    bam = str(res.bam_path)
    scan = dataclasses.replace(PipelineConfig().scan, sj_strand_tag=_native_detect_sj_tag(bam))
    st, sm, flm, buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    cfg = PipelineConfig().calibration
    R = len(idx.region_df)
    return sc, pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg, R, gdna_present


def _pct(a):
    a = np.asarray(a, dtype=np.float64)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return (np.nan, np.nan, np.nan)
    return tuple(np.percentile(a, [10, 50, 90]))


def _rng_str(a):
    a = np.asarray(a, dtype=np.float64)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return "  (none)"
    return f"min={a.min():.4g} med={np.median(a):.4g} max={a.max():.4g}"


def instrument(kind, pl, ra, sm, gdna_fl_pmf, rna_fl_pmf, cfg, R, gdna_present):
    print(f"\n{'='*92}\n=== M2 SCENARIO: {kind}  (R={R} regions, gDNA_present={gdna_present}) ===\n{'='*92}")

    substrate = CalibrationSubstrate.from_payload(pl, ra)
    region_eff_len = region_eff_length(ra.region_size_bp, gdna_fl_pmf)
    boundary_eff_len = boundary_side_eff_length(gdna_fl_pmf, ra.region_size_bp)
    fl_mean = boundary_eff_length(gdna_fl_pmf)
    rna_fl_mean = boundary_eff_length(rna_fl_pmf)
    region_eff_len_rna = region_eff_length(ra.region_size_bp, rna_fl_pmf)

    balance = fit_strand_balance(sm)
    rna_sense_frac = float(balance.rna_sense_frac)
    print(f"  rna_sense_frac (kappa) = {rna_sense_frac:.4f}")

    need_count_variance = float(cfg.gdna_deconv_quantile) != 0.5
    node_density_raw = node_gdna_density(substrate, ra, region_eff_len, fl_mean,
                                         need_count_variance=need_count_variance)
    gdna_strand = fit_gdna_strand_from_substrate(
        substrate, ra, node_density_raw, boundary_eff_len, rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(cfg.gdna_strand_prior_alpha_beta),
        prior_weight=cfg.gdna_strand_prior_weight)
    gdna_strand_overdispersion = gdna_strand.gdna_strand_overdispersion
    rna_strand = fit_rna_strand_from_substrate(
        substrate, rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(cfg.rna_strand_prior_alpha_beta),
        prior_weight=cfg.rna_strand_prior_weight)
    rna_strand_overdispersion = rna_strand.rna_strand_overdispersion
    print(f"  gdna_strand_overdispersion={gdna_strand_overdispersion:.4g} "
          f"rna_strand_overdispersion={rna_strand_overdispersion:.4g}")

    _, left_split, right_split = strand_deconvolve(
        substrate, ra, rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        deconv_quantile=cfg.gdna_deconv_quantile, n_grid=cfg.n_grid)
    i0 = cfg.gdna_strand_info_scale

    def _raw_count(view):
        return view.n_unspliced_pos.astype(np.float64) + view.n_unspliced_neg.astype(np.float64)

    cleaned_left = cleaned_gdna_count(left_split, _raw_count(substrate.left), i0)
    cleaned_right = cleaned_gdna_count(right_split, _raw_count(substrate.right), i0)

    node_density = node_gdna_density(
        substrate, ra, region_eff_len, fl_mean, need_count_variance=need_count_variance,
        gdna_counts=(_raw_count(substrate.contained), cleaned_left, cleaned_right))
    region_count_frac, n_splice_upgraded = region_splice_gdna_frac(
        substrate, ra, node_density.count_gdna_frac, eff_gdna=fl_mean, eff_rna=rna_fl_mean,
        eff_gdna_region=region_eff_len, eff_rna_region=region_eff_len_rna,
        left_gdna_unspl=cleaned_left, right_gdna_unspl=cleaned_right)

    left, right = deconv_sides(
        substrate, ra, node_density, boundary_eff_len, rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        deconv_quantile=cfg.gdna_deconv_quantile, n_grid=cfg.n_grid, info_scale=i0)

    c = substrate.contained
    u_tot = c.n_unspliced_pos.astype(np.float64) + c.n_unspliced_neg.astype(np.float64)
    U = u_tot
    eff_len = np.maximum(np.asarray(region_eff_len, dtype=np.float64), 1e-9)
    mass_u = np.maximum(np.asarray(c.mass_unspliced, dtype=np.float64), 1e-9)
    geom2 = (eff_len / mass_u) ** 2
    obs = np.asarray(node_density.region_count_observable, dtype=bool)
    gdna_left = np.asarray(left.gdna_mass, dtype=np.float64)
    gdna_right = np.asarray(right.gdna_mass, dtype=np.float64)
    gdna_c = u_tot.copy()

    # ---- (5) triplet / pair mu-coverage — geometry is pass-independent; compute once ----
    sig = np.asarray(ra.signature)
    ref_id = np.asarray(ra.ref_id)
    rco, bco = count_observable_masks(sig, ref_id)
    ls, rs = same_ref_left_right(ref_id)
    la = np.zeros(R, bool); rb = np.zeros(R, bool)
    if R > 1:
        la[1:] = bco[:-1] & ls[1:]
        rb[:-1] = bco[:-1] & rs[:-1]
    triplet = la & rb & (rco & (eff_len > 1e-9))   # region + BOTH boundaries count-observable
    exon_node = ~obs                                # non-count-observable (exonic / AMBIG) nodes
    # adjacent BOTH-count-observable region pairs (the Phase-7 s_edge support)
    adj_pair_left = np.zeros(R, bool)
    if R > 1:
        same = ref_id[:-1] == ref_id[1:]
        both = rco[:-1] & rco[1:] & same
        adj_pair_left[:-1] = both   # region r is left member of an observable pair

    # strand Fisher info (pass-independent — depends on counts + kappa)
    kappa = rna_sense_frac
    I_fisher = U * (2.0 * kappa - 1.0) ** 2

    # ---- the iterative loop, instrumented ----
    prev_fg = None
    fg_history = []
    print(f"\n  -- PER-PASS (max_passes={cfg.sweep_max_passes}, K={cfg.sweep_n_grid}) --")
    for _pass in range(int(cfg.sweep_max_passes)):
        t0 = time.perf_counter()
        rho_global = (float(gdna_c[obs].sum() / max(float(eff_len[obs].sum()), 1e-9))
                      if obs.any() else 0.0)
        pts = varmean_points(substrate, ra, region_eff_len, fl_mean,
                             gdna_views=(gdna_c, gdna_left, gdna_right))
        direct = fit_direct_varmean(pts)
        imputation = fit_imputation_varmean(pts)
        t_varmean = time.perf_counter() - t0

        mu = gdna_c / eff_len
        var_d = np.where(obs, direct.predict(mu), imputation.predict(mu))
        sig2_frac = np.maximum(var_d * geom2, 1e-12)
        active_mu = mu[u_tot > 0.0]
        sigma_d_global = (1.4826 * float(np.median(np.abs(active_mu - np.median(active_mu))))
                          if active_mu.size else rho_global)
        sig2_glob = np.maximum(sigma_d_global ** 2 * geom2, 1e-12)
        tau_count = np.minimum(1.0 / sig2_frac, mass_u)
        tau_global = np.clip(1.0 / sig2_glob, 1.0, mass_u)

        # (1) rho_global in DIRECT fit range?
        lrho = np.log(max(rho_global, 1e-12))
        in_range = direct.x_lo <= lrho <= direct.x_hi

        t1 = time.perf_counter()
        regions = deconv_regions_sweep(
            substrate, ra, rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            count_gdna_frac=region_count_frac, count_precision=tau_count,
            n_grid=cfg.sweep_n_grid, rho_global=rho_global, region_eff_len=region_eff_len,
            info_scale=i0, global_tau=tau_global)
        t_sweep = time.perf_counter() - t1

        fg = np.asarray(regions.gdna_frac, dtype=np.float64)
        # (3) convergence
        if prev_fg is not None:
            delta = float(np.mean(np.abs(fg - prev_fg)))
        else:
            delta = float("nan")
        fg_history.append(fg.copy())

        # (2) realized count-prior influence vs strand Fisher info.
        # NOTE: the sweep down-weights tau_count internally by (1-w), w=I/(I+I0). Report BOTH the raw
        # tau_count (what the var~mean emits) and the EFFECTIVE precision the sweep actually applies.
        from rigel.calibration.signature import TS_NONE
        info = U * (2.0 * kappa - 1.0) ** 2
        # the sweep's internal (1-w) count down-weight, w=I/(I+I0); w=0 at TS_NONE (intergenic)
        w_strand = np.where(np.asarray(ra.strand_class) == TS_NONE, 0.0, info / (info + float(i0)))
        eff_tau_count = (1.0 - w_strand) * tau_count
        # classification on ACTIVE (U>0) nodes
        act = U > 0.0
        Ic = I_fisher[act]
        tc = tau_count[act]
        etc = eff_tau_count[act]
        n_act = int(act.sum())
        # count is a live lever: tau_count >= I (within 2x); strand dominates: tau_count << I (<0.1*I);
        # count inert: tau_count ~ 0 (< 1e-6 absolute or eff_tau ~ 0)
        live = (tc >= 0.5 * Ic) & (tc > 1e-6)
        strand_dom = (tc < 0.1 * np.maximum(Ic, 1e-12)) & (Ic > 1e-9)
        inert = tc < 1e-6
        # also: how many nodes have I==0 (strand uninformative -> count is the only magnitude lever)
        strand_silent = Ic <= 1e-12

        print(f"\n  PASS {_pass}: rho_global={rho_global:.5g}  "
              f"log(rho) in DIRECT[{direct.x_lo:.3g},{direct.x_hi:.3g}]? {in_range}  "
              f"(direct.edf={direct.edf:.2f} lam={direct.lam:.3g}; "
              f"imput.edf={imputation.edf:.2f})")
        print(f"    convergence mean|f_g-prev| = {delta if not np.isnan(delta) else float('nan'):.3e}"
              + ("  (pass 0 — no prev)" if np.isnan(delta) else ""))
        print(f"    wall: var~mean refit = {t_varmean*1e3:7.1f} ms   sweep = {t_sweep*1e3:7.1f} ms")
        p = _pct(geom2[act]); print(f"    geom2     p10/p50/p90 = {p[0]:.4g} / {p[1]:.4g} / {p[2]:.4g}")
        p = _pct(var_d[act]); print(f"    var_d     p10/p50/p90 = {p[0]:.4g} / {p[1]:.4g} / {p[2]:.4g}")
        p = _pct(sig2_frac[act]); print(f"    sig2_frac p10/p50/p90 = {p[0]:.4g} / {p[1]:.4g} / {p[2]:.4g}")
        p = _pct(tau_count[act]); print(f"    tau_count p10/p50/p90 = {p[0]:.4g} / {p[1]:.4g} / {p[2]:.4g}")
        p = _pct(eff_tau_count[act]); print(f"    eff_tau_c p10/p50/p90 = {p[0]:.4g} / {p[1]:.4g} / {p[2]:.4g}  (after (1-w) sweep down-weight)")
        print(f"    count-influence over {n_act} active nodes:")
        print(f"      tau_count >~ I (count LIVE lever): {100*live.mean():5.1f}%")
        print(f"      tau_count << I (strand dominates): {100*strand_dom.mean():5.1f}%")
        print(f"      tau_count ~ 0  (count inert)      : {100*inert.mean():5.1f}%")
        print(f"      strand silent I~0 (count is only): {100*strand_silent.mean():5.1f}%")

        if not np.isnan(delta) and delta < 1e-3:
            print(f"    -> CONVERGED (delta < 1e-3) at pass {_pass}; loop would break here")
            prev_fg = fg
            gdna_c = np.asarray(regions.gdna_mass, dtype=np.float64)
            break
        prev_fg = fg
        gdna_c = np.asarray(regions.gdna_mass, dtype=np.float64)

    # (3) full convergence summary
    print("\n  -- convergence trace (mean|f_g - prev_f_g|) --")
    deltas = [float(np.mean(np.abs(fg_history[i] - fg_history[i-1]))) for i in range(1, len(fg_history))]
    print("    " + " -> ".join(f"{d:.3e}" for d in deltas) if deltas else "    (single pass)")
    if deltas:
        mono = all(deltas[i] <= deltas[i-1] + 1e-12 for i in range(1, len(deltas)))
        print(f"    monotone-decreasing? {mono}   below 1e-3 within {cfg.sweep_max_passes} passes? "
              f"{any(d < 1e-3 for d in deltas)}")

    # (5) mu-coverage — at the converged gDNA estimate
    mu = gdna_c / eff_len
    print("\n  -- (5) mu-coverage (gDNA density mu = gdna_c/eff_len, at converged estimate) --")
    print(f"    all-three-observable TRIPLETS: n={int(triplet.sum())}  mu {_rng_str(mu[triplet])}")
    print(f"    adjacent BOTH-observable PAIRS: n={int(adj_pair_left.sum())}  "
          f"mu {_rng_str(mu[adj_pair_left & (U>0)])}")
    print(f"    EXON (non-count-observable) nodes: n={int((exon_node & (U>0)).sum())}  "
          f"mu {_rng_str(mu[exon_node & (U>0)])}")
    # does the training reach the exon means?
    if triplet.sum() and (exon_node & (U > 0)).sum():
        tmax = np.nanmax(mu[triplet]) if triplet.sum() else np.nan
        ex = mu[exon_node & (U > 0)]
        frac_exon_above = float(np.mean(ex > tmax)) if ex.size else float("nan")
        print(f"    fraction of EXON nodes with mu ABOVE the max triplet mu (= extrapolation region): "
              f"{100*frac_exon_above:.1f}%")
    print(f"    n_splice_upgraded = {n_splice_upgraded}")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    kinds = ["capon", "capoff", "zerodna"] if which == "all" else [which]
    for kind in kinds:
        with tempfile.TemporaryDirectory() as work:
            try:
                sc, pl, ra, sm, gpmf, rpmf, cfg, R, gpres = build_scenario(kind, work)
                instrument(kind, pl, ra, sm, gpmf, rpmf, cfg, R, gpres)
            finally:
                try:
                    sc.cleanup()
                except Exception:
                    pass


if __name__ == "__main__":
    main()
