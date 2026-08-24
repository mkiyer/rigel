"""Gates for the RNA-anchored evidence factor (`calibration.rna_anchor`).

The factor anchors the RNA SIDE of the unspliced count — certified splice flux at complete-flank
exons, the adjacent intron's excess-over-background (nascent) rate at eligible boundaries — so no
gDNA enrichment estimate appears anywhere and hybrid capture cannot mis-scale it. Written BEFORE
the module (the falsification-first rule); each gate was watched firing against a deliberately
broken build before the real one landed.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration import rna_anchor as RA
from rigel.calibration.simplex_logodds import _logodds_grid

K = 61
WINDOW = 10.0
_, FG = _logodds_grid(K, WINDOW)


# ── the pure exon-anchor rows ────────────────────────────────────────────────────────────────────


def _exon_rows(unspl, flux, opp, eff_r, pad=True):
    """Probe slots, optionally padded with a well-behaved calibration population so the left-tail
    dispersion estimator is live (it needs a real population; a handful of probes cannot feed it —
    that limit has its own gate below)."""
    unspl = list(np.asarray(unspl, np.float64))
    flux = list(np.asarray(flux, np.float64))
    opp = list(np.asarray(opp, np.float64))
    eff_r = list(np.asarray(eff_r, np.float64))
    n_probe = len(unspl)
    if pad:
        rng = np.random.default_rng(11)
        for _ in range(60):
            f, o, er = 300.0, 200.0, 400.0
            mu = (f + 0.5) / o * er
            unspl.append(float(mu * np.exp(rng.normal(0.0, 0.1))))
            flux.append(f)
            opp.append(o)
            eff_r.append(er)
    rows = RA.exon_anchor_rows(
        unspl=np.asarray(unspl, np.float64),
        flux=np.asarray(flux, np.float64),
        flux_opportunity=np.asarray(opp, np.float64),
        rna_opportunity=np.asarray(eff_r, np.float64),
        fg_grid=FG,
    )
    return rows[:n_probe]


def _argbest_fg(row):
    return float(FG[int(np.argmax(row))])


def test_zero_flux_delivers_gdna_and_high_flux_refutes_it():
    """The two ends of the mechanism: no accounted RNA ⇒ the unspliced mass is gDNA; flux that
    fully accounts the count ⇒ it is RNA. Scored on where each factor row PEAKS."""
    rows = _exon_rows(
        unspl=[1000.0, 1000.0],
        flux=[0.0, 500.0],
        opp=[200.0, 200.0],
        eff_r=[400.0, 400.0],  # exon 1: expected RNA = (500.5/200)*400 ≈ 1001 ≈ the whole count
    )
    assert _argbest_fg(rows[0]) > 0.95, "zero flux must deliver: the count is gDNA"
    assert _argbest_fg(rows[1]) < 0.05, "flux accounting the whole count must refute gDNA"


def test_the_zero_flux_anchor_is_jeffreys_not_a_wall():
    """flux = 0 takes the Jeffreys-Poisson posterior mean (½/opportunity), never a hard zero: the
    v0 prototype floored the mean at ~0 and confidently asserted gDNA at zero-observed-flux exons
    with real RNA — the measured zero-control failure this gate pins. A LARGE opportunity earns a
    tighter zero than a small one."""
    small = _exon_rows([100.0], [0.0], [50.0], [400.0])[0]
    large = _exon_rows([100.0], [0.0], [5000.0], [400.0])[0]
    # both deliver, but the small-opportunity anchor must be measurably softer on the RNA end
    assert large[0] < small[0] - 1.0, (
        "an anchor with 100x the opportunity behind its zero must punish the RNA end harder "
        f"(small {small[0]:.2f}, large {large[0]:.2f})"
    )


def test_pair_dispersion_widens_the_factor():
    """The left-tail dispersion estimator: a population whose counts scatter widely below their
    predictions must yield WIDER factors than a tight population — the owner's the-anchor-is-an-
    imputation requirement, whose absence was the v0 zero-control failure."""
    rng = np.random.default_rng(3)
    n = 200
    flux = np.full(n, 300.0)
    opp = np.full(n, 200.0)
    eff_r = np.full(n, 400.0)
    mu = (flux + 0.5) / opp * eff_r
    tight = mu * np.exp(rng.normal(0.0, 0.05, n))
    wide = mu * np.exp(rng.normal(0.0, 1.0, n))
    r_tight = _exon_rows(tight, flux, opp, eff_r, pad=False)
    r_wide = _exon_rows(wide, flux, opp, eff_r, pad=False)
    # spread of the log-factor over the grid = the factor's strength; wide data ⇒ weaker factor
    s_tight = float(np.median(np.ptp(r_tight, axis=1)))
    s_wide = float(np.median(np.ptp(r_wide, axis=1)))
    assert s_wide < 0.5 * s_tight, (
        f"wide pairs must weaken the factor ({s_wide:.2f} vs {s_tight:.2f})"
    )


def test_too_few_pairs_makes_the_factor_near_flat():
    """Below the estimator's minimum population the pair dispersion is unknowable, and an unknown
    dispersion must read as a WEAK factor, never a confident one. ⚠ The probes deliberately
    STRADDLE their prediction (three sit below it) so a weakened population guard would actually
    fit a variance from them — a fixture with no negative residuals cannot see that defect
    (TRAPS: could-the-arm-have-fired; found by the perturbation sweep reading INERT)."""
    mu = (300.0 + 0.5) / 200.0 * 400.0  # = 601: the slots' shared prediction
    obs = [mu * f for f in (0.95, 0.98, 0.99, 1.4, 1.6)]
    rows = _exon_rows(obs, [300.0] * 5, [200.0] * 5, [400.0] * 5, pad=False)
    assert float(np.ptp(rows)) < 1.0, (
        "with too few pairs to estimate dispersion the factor must be weak"
    )


def test_flank_disagreement_measures_transport_dispersion_gdna_free():
    """The PRIMARY dispersion estimator: two complete flanks of one exon estimate the same RNA, so
    their disagreement measures the transport dispersion with no gDNA anywhere in the computation
    — it cannot starve on contaminated data the way the left-tail fallback does (at high gDNA no
    observation ever falls below its RNA prediction). Tight flanks ⇒ small V; scattered ⇒ large;
    and the fluxes' own counting noise is subtracted, so pure-Poisson disagreement reads ≈ 0."""
    rng = np.random.default_rng(7)
    n = 100
    opp = np.full(n, 200.0)
    base = rng.uniform(200.0, 600.0, n)
    tight = RA.flank_disagreement_log_variance(
        base, opp, base * np.exp(rng.normal(0, 0.02, n)), opp
    )
    wide = RA.flank_disagreement_log_variance(base, opp, base * np.exp(rng.normal(0, 0.8, n)), opp)
    assert tight is not None and wide is not None
    assert wide > 10.0 * max(tight, 1e-6), (
        f"scattered flanks must read larger (t {tight}, w {wide})"
    )
    few = RA.flank_disagreement_log_variance(base[:5], opp[:5], base[:5], opp[:5])
    assert few is None, "below the population minimum the estimator must refuse"


def test_flank_disagreement_recovers_a_known_dispersion_absolutely():
    """SCALE gate, two regimes (a ratio-only gate let two scale defects pass as inert — found by
    the perturbation sweep — and these absolute checks make them fire):

    * LARGE flux, known log-sd: counting noise is negligible, so recovery within ±50 % pins the
      pair-halving (each flank carries HALF the two-flank disagreement);
    * SMALL flux, ZERO true dispersion: the disagreement is pure counting noise, so the estimate
      must read ≈ 0 — which is exactly what the counting-noise subtraction achieves."""
    rng = np.random.default_rng(19)
    n = 4000
    opp = np.full(n, 200.0)

    s_true = 0.3
    fa = rng.poisson(500.0 * np.exp(rng.normal(0, s_true, n))).astype(float)
    fb = rng.poisson(500.0 * np.exp(rng.normal(0, s_true, n))).astype(float)
    V = RA.flank_disagreement_log_variance(fa, opp, fb, opp)
    assert V is not None
    assert 0.5 * s_true**2 < V < 1.5 * s_true**2, (
        f"injected {s_true**2:.3f}, recovered {V:.3f} — outside the 50% band"
    )

    fa0 = rng.poisson(8.0, n).astype(float)
    fb0 = rng.poisson(8.0, n).astype(float)
    V0 = RA.flank_disagreement_log_variance(fa0, opp, fb0, opp)
    assert V0 is not None and V0 < 0.03, (
        f"pure counting noise must read ~0 dispersion, got {V0:.3f}"
    )


def test_the_fitted_center_recenters_the_anchor():
    """A population whose observations sit uniformly 25 % ABOVE their predictions with tight
    scatter and NO gDNA must be read as a transport center, not as gDNA: the factor's peak must
    sit at the RNA end. The predecessor pinned the center at the raw prediction and converted the
    measured ~10-25 % transport bias into phantom gDNA at every RNA-rich exon of a clean library
    (whole-library 17k → 235k at the zero control) — this gate pins the repair."""
    rng = np.random.default_rng(23)
    n = 400
    flux = np.full(n, 300.0)
    opp = np.full(n, 200.0)
    eff_r = np.full(n, 400.0)
    mu = (flux + 0.5) / opp * eff_r
    # scatter 0.3 STRADDLES the +25% offset, so the lower-quantile fit is live — with tighter
    # scatter no residual is negative, the fit correctly refuses, and the factor goes flat, which
    # made the first version of this gate pass vacuously with the center pinned at zero
    # (TRAPS: could-the-arm-have-fired, caught by the perturbation sweep reading MISSED).
    obs = mu * 1.25 * np.exp(rng.normal(0.0, 0.3, n))
    rows = RA.exon_anchor_rows(
        unspl=obs, flux=flux, flux_opportunity=opp, rna_opportunity=eff_r, fg_grid=FG
    )
    assert float(np.median(np.ptp(rows, axis=1))) > 1.0, "the factor must be LIVE for this gate"
    med_peak = float(np.median(FG[np.argmax(rows, axis=1)]))
    assert med_peak < 0.05, (
        f"a uniform +25% transport offset must be absorbed by the center, not read as gDNA "
        f"(median peak f_g = {med_peak:.3f})"
    )


def test_a_gdna_contaminated_tail_refuses_the_center_fit():
    """The fit's self-consistency guard: when gDNA inflates most residuals (a `g98`-like
    population), the observed negative fraction falls far below what the fitted log-normal
    predicts, and the fit must REFUSE rather than return a spuriously positive center — measured,
    an unguarded center on `g98`'s tail gave back most of the high-gDNA win."""
    rng = np.random.default_rng(29)
    n = 2000
    mu = np.full(n, 600.0)
    # transport scatter 0.2, plus gDNA inflating 97% of slots by 2-50x — only ~3% stay clean
    clean = rng.random(n) < 0.03
    obs = mu * np.exp(rng.normal(0.0, 0.2, n))
    obs[~clean] *= rng.uniform(2.0, 50.0, int((~clean).sum()))
    assert RA.left_fit_center_spread(obs, mu) is None, (
        "a 97%-contaminated residual population must refuse the fit"
    )
    # and the same scatter with NO contamination must be accepted, with the center near zero
    fit = RA.left_fit_center_spread(mu * np.exp(rng.normal(0.0, 0.2, n)), mu)
    assert fit is not None and abs(fit[0]) < 0.1


# ── the pure boundary-anchor rows ────────────────────────────────────────────────────────────────


def _boundary_rows(unspl_b, unspl_i, eff_g_i, eff_r_b, rho_bg):
    """Probe slots padded with a well-behaved boundary population (nascent-rich introns whose
    boundaries scatter tightly around the prediction) so the left-tail estimator is live."""
    unspl_b = list(np.asarray(unspl_b, np.float64))
    unspl_i = list(np.asarray(unspl_i, np.float64))
    eff_g_i = list(np.asarray(eff_g_i, np.float64))
    eff_r_b = list(np.asarray(eff_r_b, np.float64))
    n_probe = len(unspl_b)
    rng = np.random.default_rng(13)
    for _ in range(60):
        ci, eg, er = 9000.0, 10000.0, 600.0
        mu = ((ci + 0.5) / eg - rho_bg) * er
        unspl_b.append(float(mu * np.exp(rng.normal(0.0, 0.1))))
        unspl_i.append(ci)
        eff_g_i.append(eg)
        eff_r_b.append(er)
    rows = RA.boundary_anchor_rows(
        unspl_boundary=np.asarray(unspl_b, np.float64),
        unspl_intron=np.asarray(unspl_i, np.float64),
        gdna_opportunity_intron=np.asarray(eff_g_i, np.float64),
        rna_opportunity_boundary=np.asarray(eff_r_b, np.float64),
        background_rate=float(rho_bg),
        fg_grid=FG,
    )
    return rows[:n_probe]


def test_a_clean_intron_delivers_the_boundary_and_a_nascent_one_yields():
    """The boundary mechanism: an intron at the background rate predicts zero nascent ⇒ the
    boundary crossing is gDNA; an intron far above background licenses RNA at the boundary."""
    rows = _boundary_rows(
        unspl_b=[500.0, 500.0],
        unspl_i=[1000.0, 9000.0],
        eff_g_i=[10000.0, 10000.0],
        eff_r_b=[600.0, 600.0],
        rho_bg=0.1,  # intron 0 sits AT background; intron 1 is 9x above it
    )
    assert _argbest_fg(rows[0]) > 0.95, "an at-background intron must deliver the boundary as gDNA"
    assert _argbest_fg(rows[1]) < 0.6, "a nascent-rich intron must license RNA at the boundary"


def test_a_shallow_anchor_speaks_more_weakly_than_a_deep_one():
    """The subtraction-noise propagation (the v1→v2 fix): the SAME nascent rate estimated from a
    10x shallower intron must yield a measurably weaker factor — the anchor's depth, not just its
    point estimate, sets how hard the factor may push. (An earlier draft demanded near-FLATNESS
    when the excess sits inside its noise; that was the wrong ask — a 500-count crossing against a
    nascent band of 0–4 fragments is a justified gDNA call at either end of the band. What must
    scale is confidence, and this pins it.)"""
    rows = _boundary_rows(
        unspl_b=[500.0, 500.0],
        unspl_i=[2000.0, 20.0],  # both 2x background, backed by 100x different intron depths
        eff_g_i=[10000.0, 100.0],
        eff_r_b=[600.0, 600.0],
        rho_bg=0.1,
    )
    deep, shallow = float(np.ptp(rows[0])), float(np.ptp(rows[1]))
    assert shallow < 0.5 * deep, (
        f"a 100x shallower anchor must weaken the factor (deep {deep:.1f}, shallow {shallow:.1f})"
    )


# ── assembly on the chain ────────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def anchor_toy(tmp_path_factory):
    """The smallest real substrate: the prior_vs_oracle toy (contaminated, nascent in the introns),
    scanned once, with the chain pieces and full `calibrate` kwargs built the production way."""
    import dataclasses
    import importlib.util
    import sys
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "tpo_for_anchor", Path(__file__).parent / "test_prior_vs_oracle.py"
    )
    m = importlib.util.module_from_spec(spec)
    sys.modules["tpo_for_anchor"] = m
    spec.loader.exec_module(m)
    toy = m.toy.__wrapped__(tmp_path_factory)

    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.region_chain import build_region_chain
    from rigel.calibration.region_geometry import build_region_geometry, build_region_statics
    from rigel.calibration.sj_opportunity import crossing_probability_from_index
    from rigel.calibration.splice_graph import (
        build_boundary_flags_array,
        build_sj_geometry_arrays,
    )
    from rigel.calibration.structural_claims import build_structural_claims
    from rigel.calibration.substrate import CalibrationSubstrate
    from rigel.config import PipelineConfig
    from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer

    index = toy.index
    scan = dataclasses.replace(
        PipelineConfig().scan, sj_strand_tag=_native_detect_sj_tag(str(toy.bam_path))
    )
    _stats, strand_model, _buf, payload = scan_and_buffer(str(toy.bam_path), index, scan)
    ra = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
    fl = build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, int(payload.max_length)),
        gdna_opportunity=gdna_opportunity_from_index(index, int(payload.max_length)),
    )
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    sj = build_sj_geometry_arrays(index)
    bflags = build_boundary_flags_array(index)
    geometry = build_region_geometry(chain, substrate, ra, sj, fl.gdna_pmf, fl.rna_pmf, None)
    statics = build_region_statics(chain, ra, bflags)
    claims = build_structural_claims(chain, statics)
    kwargs = dict(
        region_arrays=ra,
        strand_model=strand_model,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        sj=sj,
        boundary_flags=bflags,
    )
    return dict(
        payload=payload,
        chain=chain,
        statics=statics,
        geometry=geometry,
        region_arrays=ra,
        claims=claims,
        calibrate_kwargs=kwargs,
    )


def test_the_factor_touches_only_the_two_claimed_populations(anchor_toy, monkeypatch):
    """No leakage: rows are exactly zero at every slot outside complete-flank exons and eligible
    ss-intron boundaries — the factor may not whisper at intergenic space, introns, or anywhere
    the claim set does not license.

    ⚠ On this toy the dispersion estimator refuses — every slot is gDNA-contaminated ABOVE its RNA
    prediction, so there are no negative residuals to fit — and the factor is legitimately ALL-FLAT
    (the conservative unknown-dispersion limit, asserted first). The leakage check then pins the
    estimator to a known variance so the factor goes live (the estimator has its own unit gates;
    TRAPS: could-the-arm-have-fired — a flat factor cannot leak, so flatness proves nothing)."""
    t = anchor_toy
    flat = RA.build_rna_anchor_factor(
        t["chain"],
        t["statics"],
        t["geometry"],
        t["region_arrays"],
        n_grid=K,
        logodds_window=WINDOW,
    )
    assert flat is not None and not np.any(np.ptp(flat, axis=1) > 0.0), (
        "below the pair minimum the factor must be flat everywhere"
    )
    monkeypatch.setattr(RA, "left_fit_center_spread", lambda observed, predicted: (0.0, 0.01))
    out = RA.build_rna_anchor_factor(
        t["chain"],
        t["statics"],
        t["geometry"],
        t["region_arrays"],
        n_grid=K,
        logodds_window=WINDOW,
    )
    touched = np.ptp(out, axis=1) > 0.0
    assert np.any(touched), (
        "with the minimum lowered the toy's eligible slots must go live "
        "(TRAPS: could-the-arm-have-fired)"
    )
    eligible = RA.eligible_slots(t["chain"], t["statics"], t["geometry"], t["region_arrays"])
    assert not np.any(touched & ~eligible), (
        "the factor moved a slot outside its claimed populations"
    )


def test_the_flag_gates_the_builder_both_ways(anchor_toy, monkeypatch):
    """`rna_anchor = False` never reaches the builder; `True` always does — asserted by making the
    builder explosive, because an output comparison can pass while a mechanism runs and cancels
    (TRAPS: could-the-arm-have-fired)."""
    import sys

    import rigel.calibration.calibrate  # noqa: F401

    calibrate_mod = sys.modules["rigel.calibration.calibrate"]
    from rigel.config import CalibrationConfig

    class Boom(Exception):
        pass

    def boom(*a, **kw):
        raise Boom

    monkeypatch.setattr(calibrate_mod, "build_rna_anchor_factor", boom)
    t = anchor_toy
    calibrate_mod.calibrate(
        payload=t["payload"],
        config=CalibrationConfig(rna_anchor=False, message_propagation=False),
        **t["calibrate_kwargs"],
    )  # flag OFF: the explosive builder must never be reached
    with pytest.raises(Boom):
        calibrate_mod.calibrate(
            payload=t["payload"],
            config=CalibrationConfig(rna_anchor=True, message_propagation=False),
            **t["calibrate_kwargs"],
        )


def test_the_flag_exists_and_defaults_on():
    from rigel.config import CalibrationConfig

    assert CalibrationConfig().rna_anchor is True
