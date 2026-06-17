"""Monotone-increasing P-spline (SCAM) var~mean fitter + the DIRECT / node-PAIR IMPUTATION builders.

Pins: monotonicity by construction, power-law recovery, robustness to an outlier, flat extrapolation
outside the fit range, the too-few-points fallback, the diagnostics dataframe, the DIRECT own-count
extractor, and the unified node-PAIR imputation reliability (CALIBRATION_PLAN_v5 §3).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.variance_model import (
    MonotoneVarMean,
    direct_points,
    fit_pair_imputation_rna_varmean,
    fit_pair_imputation_varmean,
)


def _powerlaw(n=400, a=0.01, b=2.0, seed=0):
    rng = np.random.default_rng(seed)
    mean = np.exp(rng.uniform(np.log(0.01), np.log(10.0), n))
    var = a * mean**b * np.exp(rng.normal(0.0, 0.3, n))  # log-normal noise around the law
    return mean, var


def test_fit_is_monotone_on_a_grid():
    mean, var = _powerlaw()
    fit = MonotoneVarMean.fit(mean, var)
    grid = np.logspace(-2, 1, 200)
    pred = fit.predict(grid)
    assert np.all(np.diff(pred) >= -1e-9)  # non-decreasing everywhere


def test_recovers_power_law_trend():
    mean, var = _powerlaw(a=0.02, b=2.0, n=600)
    fit = MonotoneVarMean.fit(mean, var)
    xs = np.array([0.05, 0.5, 5.0])
    pred = fit.predict(xs)
    truth = 0.02 * xs**2.0
    # within a factor of ~2 of the planted trend across 2 decades
    assert np.all(pred > truth / 2.5) and np.all(pred < truth * 2.5)


def test_robust_to_a_wild_outlier():
    mean, var = _powerlaw(n=400, seed=1)
    var = var.copy()
    var[mean.argmax()] = 1e6  # a wild high-variance outlier at the top mean
    robust = MonotoneVarMean.fit(mean, var, robust_iters=2)
    naive = MonotoneVarMean.fit(mean, var, robust_iters=0)
    top = float(mean.max())
    # robust fit at the top is far less dragged toward the 1e6 outlier than the naive fit
    assert robust.predict(np.array([top]))[0] < 0.5 * naive.predict(np.array([top]))[0]


def test_predict_clips_to_fit_range():
    mean, var = _powerlaw()
    fit = MonotoneVarMean.fit(mean, var)
    lo_in = fit.predict(np.array([float(np.exp(fit.x_lo))]))[0]
    lo_out = fit.predict(np.array([1e-9]))[0]  # far below the fit range
    hi_in = fit.predict(np.array([float(np.exp(fit.x_hi))]))[0]
    hi_out = fit.predict(np.array([1e9]))[0]  # far above
    assert np.isclose(lo_out, lo_in, rtol=1e-6)  # flat extrapolation, not a runaway
    assert np.isclose(hi_out, hi_in, rtol=1e-6)


def test_too_few_points_falls_back_monotone():
    mean = np.array([0.1, 1.0, 5.0])
    var = np.array([1e-3, 0.1, 2.0])
    fit = MonotoneVarMean.fit(mean, var)  # < k points ⇒ power-law fallback
    pred = fit.predict(np.logspace(-1, 1, 50))
    assert np.all(np.diff(pred) >= -1e-9)


def test_to_dataframe_has_points_and_curve():
    mean, var = _powerlaw(n=200)
    fit = MonotoneVarMean.fit(mean, var)
    df = fit.to_dataframe()
    assert set(df["kind"].unique()) == {"point", "curve"}
    assert (df["var"] > 0).all()


def test_direct_points_extracts_count_observable_only():
    # DirectPoints is the count-observable subset of the triplet extraction: one point per
    # count-observable region with >= 2 measurements. An exonic (not count-observable) region never
    # contributes a DIRECT point even if it has flanking observable crossings. Side densities use the
    # per-side density length E[min(ℓ,L_side)] (passed as boundary_side_eff_len), not fl_mean.
    from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS

    class _Ra:
        # 3 same-ref intron-only regions (all count-observable) flanked so internal boundaries are
        # observable (no shared exon bit) → the middle region gets contained + both flanks (k=3).
        signature = np.array([BIT_INTRON_POS, BIT_INTRON_POS, BIT_INTRON_POS, BIT_EXON_POS])
        ref_id = np.zeros(4, dtype=np.int64)

    class _C:
        n_unspliced_pos = np.array([40.0, 60.0, 90.0, 99.0])
        n_unspliced_neg = np.zeros(4)

    class _S:
        n_unspliced_pos = np.array([5.0, 6.0, 7.0, 0.0])
        n_unspliced_neg = np.zeros(4)

    class _Sub:
        contained = _C()
        left = _S()
        right = _S()

    eff = np.array([10.0, 10.0, 30.0, 10.0])
    side_len = np.array([2.0, 2.0, 2.0, 2.0])  # per-side density length E[min(ℓ,L_side)] (here ≈ fl_mean)
    pts = direct_points(_Sub(), _Ra(), eff, side_len)
    # the exonic region (index 3) is excluded; only count-observable regions with >=2 measurements appear.
    # kcount is 2 or 3 (own + observable flanks); dof = kcount - 1.
    assert pts.mean.size >= 1
    assert np.all((pts.kcount >= 2.0) & (pts.kcount <= 3.0))
    # all means/raw_vars finite & positive
    assert np.all(np.isfinite(pts.mean)) and np.all(pts.mean > 0.0)
    assert np.all(np.isfinite(pts.raw_var)) and np.all(pts.raw_var > 0.0)


def test_pair_imputation_builder_densifies_and_is_monotone():
    # The node-PAIR builder: one point per (observable side → eligible region) adjacency, at
    # mean = region density, raw_var = (d_region − d_side)². A both-flanks region → 2 points; a
    # one-flank region → 1 point (the densification the both-sides triplet missed).
    rng = np.random.default_rng(7)
    n = 300
    rd = np.exp(rng.uniform(np.log(0.05), np.log(8.0), n))
    # planted: each side density is the region density plus log-normal multiplicative error.
    ld = rd * np.exp(rng.normal(0.0, 0.4, n))
    rrd = rd * np.exp(rng.normal(0.0, 0.4, n))
    ref = np.zeros(n, dtype=np.int64)
    elig = rng.random(n) < 0.7
    lok = rng.random(n) < 0.8
    rok = rng.random(n) < 0.8
    fit = fit_pair_imputation_varmean(
        rd, ld, rrd, region_eligible=elig, left_ok=lok, right_ok=rok, ref_id=ref
    )
    # point count == number of (eligible region, observable flank) adjacencies (both-flanks → 2)
    expected = int((elig & lok).sum() + (elig & rok).sum())
    assert fit.fit_mean.size == expected
    # monotone over the fit range
    g = np.logspace(np.log10(np.exp(fit.x_lo)), np.log10(np.exp(fit.x_hi)), 100)
    assert np.all(np.diff(fit.predict(g)) >= -1e-9)
    # the queried axis is the REGION density (means are drawn from rd[eligible & flank-ok])
    assert fit.fit_mean.min() >= rd[elig].min() * (1 - 1e-9)


def test_pair_imputation_single_flank_contributes():
    # A region eligible with only ONE observable flank still contributes exactly one point.
    rd = np.array([1.0, 2.0, 3.0])
    ld = np.array([0.5, 0.0, 0.0])
    rrd = np.array([0.0, 0.0, 0.0])
    elig = np.array([True, True, True])
    lok = np.array([True, False, False])  # only region 0 has a left flank
    rok = np.array([False, False, False])
    ref = np.zeros(3, dtype=np.int64)
    fit = fit_pair_imputation_varmean(
        rd, ld, rrd, region_eligible=elig, left_ok=lok, right_ok=rok, ref_id=ref
    )
    assert fit.fit_mean.size == 1
    assert np.isclose(fit.fit_mean[0], 1.0)  # mean = region density of region 0


def test_jensen_offset_values():
    # The Jensen df-offset Δ = log(dof/2) − ψ(dof/2): positive, decreasing in dof, → 0 as dof → ∞.
    from rigel.calibration.variance_model import _jensen_offset

    off = _jensen_offset(np.array([1.0, 2.0, 1e6]))
    assert off[0] == pytest.approx(1.2703628454614782, rel=1e-9)  # dof=1 (k=2 disagreement)
    assert off[1] == pytest.approx(0.5772156649015329, rel=1e-9)  # dof=2 (Euler–Mascheroni)
    assert off[2] == pytest.approx(0.0, abs=1e-5)  # dof → ∞ (no correction)
    assert off[0] > off[1] > off[2] >= 0.0


def test_jensen_offset_inflates_recovered_variance():
    # With dof passed, the fit target log(var) is shifted UP by Δ_k>0 ⇒ the recovered variance is
    # uniformly larger than the un-corrected fit (removes the small-dof over-confidence).
    mean, var = _powerlaw(n=400, seed=3)
    dof = np.full(mean.shape[0], 1.0)  # every point a k=2 disagreement
    corrected = MonotoneVarMean.fit(mean, var, dof=dof)
    plain = MonotoneVarMean.fit(mean, var)  # dof=None ⇒ back-compat, no offset
    g = np.logspace(-1.5, 0.5, 40)
    pc, pp = corrected.predict(g), plain.predict(g)
    assert np.all(pc > pp)  # inflated everywhere
    # the inflation ≈ exp(Δ_1) = exp(1.2704) ≈ 3.56× (the verified k=2 factor)
    ratio = np.median(pc / pp)
    assert ratio == pytest.approx(np.exp(1.2703628454614782), rel=0.15)


# --- Phase A: the RNA imputation var~mean reliability fit -----------------------------------------


def _multi_exon_single_strand_substrate(work_dir, *, gdna_abundance=120.0, nrna_abundance=25.0):
    """A real multi-exon single-strand scenario → (substrate, ra, region_eff_len_rna,
    rna_boundary_side_eff_len, fg, left_split, right_split, cleaned_left, cleaned_right) for the RNA-builder
    test.

    4-exon genes (alternating +/−) give INTERNAL exon regions flanked by introns on both sides, so both
    adjacent boundary sides carry same-strand spliced junction crossings — the both-sides-eligible set the
    RNA imputation fit trains on. Runs the calibrator to convergence and recomputes the per-side strand
    split + cleaned counts exactly as calibrate() does.
    """
    import dataclasses

    from rigel.calibration.calibrate import calibrate
    from rigel.calibration.effective_length import (
        boundary_side_eff_length,
        region_eff_length,
    )
    from rigel.calibration.fl import build_fl_models, gdna_fl_mass
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.strand_balance import fit_strand_balance
    from rigel.calibration.strand_deconv import cleaned_gdna_count, strand_deconvolve
    from rigel.calibration.substrate import CalibrationSubstrate
    from rigel.config import PipelineConfig
    from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
    from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
    from rigel.splice import SpliceType

    win = 14000
    n_genes = 12
    glen = (n_genes + 2) * win
    sc = Scenario("rna_varmean", genome_length=glen, seed=11, work_dir=str(work_dir))
    rng = np.random.default_rng(3)
    for gi in range(n_genes):
        base = (gi + 1) * win
        strand = "+" if gi % 2 == 0 else "-"
        exons = [
            (base + 1000, base + 2000),
            (base + 4000, base + 5000),
            (base + 7000, base + 8000),
            (base + 10000, base + 11000),
        ]
        ab = int(rng.integers(120, 260))
        sc.add_gene(f"g{gi}", strand, [{"t_id": f"g{gi}_t", "exons": exons, "abundance": ab}])
    gd = (
        GDNAConfig(abundance=gdna_abundance, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)
        if gdna_abundance > 0
        else None
    )
    res = sc.build_oracle(
        n_fragments=max(60000, n_genes * 4000),
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=100, strand_specificity=0.99, seed=11),
        gdna_config=gd, nrna_abundance=float(nrna_abundance),
    )
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
    mass_unspl = np.asarray(sub.contained.mass_unspliced, dtype=np.float64)
    fg = np.where(mass_unspl > 1e-12,
                  np.asarray(result.mass_gdna_contained, dtype=np.float64) / np.maximum(mass_unspl, 1e-12),
                  0.0)
    region_eff_len_rna = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    rna_boundary_side_eff_len = boundary_side_eff_length(fl.rna_pmf, ra.region_size_bp)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    _, lsplit, rsplit = strand_deconvolve(
        sub, ra, rna_sense_frac=kappa,
        gdna_strand_overdispersion=result.gdna_strand_overdispersion,
        rna_strand_overdispersion=result.rna_strand_overdispersion,
        deconv_quantile=ccfg.gdna_deconv_quantile, n_grid=ccfg.n_grid,
    )

    def _raw(v):
        return v.n_unspliced_pos.astype(np.float64) + v.n_unspliced_neg.astype(np.float64)

    i0 = ccfg.gdna_strand_info_scale
    cl = cleaned_gdna_count(lsplit, _raw(sub.left), i0)
    cr = cleaned_gdna_count(rsplit, _raw(sub.right), i0)
    sc.cleanup()
    return sub, ra, region_eff_len_rna, rna_boundary_side_eff_len, fg, lsplit, rsplit, cl, cr


def test_rna_imputation_varmean_real_fit(tmp_path):
    """On a multi-exon single-strand scenario the RNA imputation fit is real (>0 pts), monotone, finite,
    and SPANS its eligible RNA-density range (no extrapolation — the 2a fit-and-query-on-the-same-axis)."""
    (sub, ra, rel_rna, rna_side_len, fg, lsplit, rsplit, cl, cr) = _multi_exon_single_strand_substrate(
        tmp_path
    )
    fit = fit_pair_imputation_rna_varmean(
        sub, ra, rel_rna, rna_side_len,
        gdna_frac=fg, left_gdna_frac=lsplit.gdna_frac, right_gdna_frac=rsplit.gdna_frac,
        cleaned_left=cl, cleaned_right=cr,
    )
    # a real fit with > 0 points
    assert fit.fit_mean.size > 0
    # monotone over the fit range
    grid = np.logspace(np.log10(np.exp(fit.x_lo)), np.log10(np.exp(fit.x_hi)), 100)
    pred = fit.predict(grid)
    assert np.all(np.diff(pred) >= -1e-9)
    # finite, positive predictions
    assert np.all(np.isfinite(pred)) and np.all(pred > 0.0)
    # the fit spans its own data (predict at the fit points is flat-clipped, never a runaway)
    assert np.all(np.isfinite(fit.predict(fit.fit_mean)))
    assert np.exp(fit.x_lo) <= float(fit.fit_mean.min()) * (1 + 1e-9)
    assert np.exp(fit.x_hi) >= float(fit.fit_mean.max()) * (1 - 1e-9)


def test_rna_imputation_varmean_nan_gdna_frac_fallback():
    """The builder must not crash on the NaN per-side gdna_frac / cleaned-count fallback path, and the
    cleaned-count fallback (side_unspliced − cleaned) must be honored where gdna_frac is NaN."""
    from rigel.calibration.substrate import CalibrationSubstrate
    from tests.calibration._synthetic import make_synthetic_payload

    payload, ra = make_synthetic_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)
    r = ra.n_regions
    # No region in this tiny 3-region payload is both-sides-eligible (no internal exon flanked by
    # spliced-carrying observable boundaries), so the fit is empty but well-defined — the builder must
    # not crash on the all-NaN per-side gdna_frac with the cleaned-count fallback.
    nan = np.full(r, np.nan)
    raw_left = sub.left.n_unspliced_pos.astype(float) + sub.left.n_unspliced_neg.astype(float)
    raw_right = sub.right.n_unspliced_pos.astype(float) + sub.right.n_unspliced_neg.astype(float)
    fit = fit_pair_imputation_rna_varmean(
        sub, ra, np.full(r, 50.0), np.full(r, 50.0),  # per-side RNA density length array
        gdna_frac=np.zeros(r), left_gdna_frac=nan, right_gdna_frac=nan,
        cleaned_left=raw_left, cleaned_right=raw_right,  # cleaned == raw ⇒ RNA-removed = 0
    )
    # an empty / thin fit is still a valid MonotoneVarMean (power-law fallback); predict is finite.
    pred = fit.predict(np.array([1.0, 5.0]))
    assert np.all(np.isfinite(pred))
