"""The abundance landscape: the pre-pass-0 total-density field and the mode census read off it.

`fit_abundance_landscape` fits the population estimator (`landscape.fit_landscape`) on the wall-exact
measured totals from `total_abundance` and then reports every local maximum with its basin, the
depleted mode (the basin containing the pooled intergenic anchor rate), the enriched mode (the
largest-mass basin strictly above it), the span `R`, and a per-region enriched-basin responsibility
`w_i`. The fixtures are synthetic Poisson populations whose truth is stated by hand, so every
assertion here is absolute rather than a recorded output. What the gates mostly hold is that no
significance threshold exists anywhere in the census: basins must PARTITION the density
(Σ basin_mass = 1), a unimodal fit must report NO enriched mode with `span_R` exactly 1 and `w ≡ 0`,
and the anchor-consistency verdict must use the depleted mode's own fitted width as its tolerance —
the density's own statement of its resolution, never a chosen number.
"""

from __future__ import annotations


import numpy as np
import pytest

from rigel.calibration.abundance_landscape import (
    fit_abundance_landscape,
)
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS
from rigel.calibration.total_abundance import RegionWallMask

import pandas as pd


def parts(
    counts,
    lengths,
    signatures,
    *,
    start_exact=None,
    end_exact=None,
    seed_starts=None,
):
    """A hand-built (substrate-like, region_arrays, wall_mask) triple.

    ``counts`` are the per-region START totals we want the landscape to see; the substrate carries
    them split over two strand columns (never evenly, so a single-column read cannot pass) and the
    END bank mirrors the same totals with a different split — the ledger closes by construction.
    """
    counts = np.asarray(counts, dtype=np.int64)
    lengths = np.asarray(lengths, dtype=np.int64)
    n = counts.shape[0]
    bounds = np.concatenate([[0], np.cumsum(lengths)])
    frame = pd.DataFrame(
        {
            "region_id": np.arange(n, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * n, dtype="string"),
            "start": bounds[:-1],
            "end": bounds[1:],
            "length": lengths,
            "signature": np.asarray(signatures, dtype=np.uint8),
        }
    )
    ra = RegionArrays.from_frame(frame, {"chr1": 0})
    lo = counts // 3
    s_bank = np.stack([counts - lo, lo], axis=1)
    e_bank = np.stack([lo, counts - lo], axis=1)

    class _Sub:
        region_start_count = s_bank
        region_end_count = e_bank

    se = np.ones(n, dtype=bool) if start_exact is None else np.asarray(start_exact, dtype=bool)
    ee = np.ones(n, dtype=bool) if end_exact is None else np.asarray(end_exact, dtype=bool)
    big = np.full(n, 1.0e6)
    mask = RegionWallMask(
        n_regions=n,
        d_low=np.where(ee, big, 0.0),
        d_high=np.where(se, big, 0.0),
        start_exact=se,
        end_exact=ee,
        w_max=500,
    )
    return _Sub(), ra, mask


def bimodal_parts(rng=None):
    """Two Poisson populations 2.5 decades apart, plus intergenic anchors AT the depleted level.

    600 intergenic/intron regions of 10 kb at rho = 0.01 (expected count 100) and 220 exon regions of
    1 kb at rho = 3.16 (expected count 3162) — decisively separated, deterministically drawn.
    """
    rng = rng or np.random.default_rng(7)
    n_lo, n_hi = 600, 220
    len_lo, len_hi = 10_000, 1_000
    rho_lo, rho_hi = 0.01, 10.0 ** (np.log10(0.01) + 2.5)
    c_lo = rng.poisson(rho_lo * len_lo, n_lo)
    c_hi = rng.poisson(rho_hi * len_hi, n_hi)
    counts = np.concatenate([c_lo, c_hi])
    lengths = np.concatenate([np.full(n_lo, len_lo), np.full(n_hi, len_hi)])
    sig = np.concatenate(
        [
            np.zeros(n_lo // 2, np.uint8),  # intergenic — the anchor pool
            np.full(n_lo - n_lo // 2, BIT_INTRON_POS, np.uint8),
            np.full(n_hi, BIT_EXON_POS, np.uint8),
        ]
    )
    return counts, lengths, sig, rho_lo, rho_hi


# ---------------------------------------------------------------------------
# the bimodal fixture — the capture shape
# ---------------------------------------------------------------------------


def test_a_bimodal_field_yields_TWO_modes_at_the_right_places():
    counts, lengths, sig, rho_lo, rho_hi = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    assert al is not None
    assert al.enriched is not None
    # located within the fit's own resolution: the mode's fitted width, not a chosen tolerance
    assert abs(al.depleted.log_rho - np.log(rho_lo)) <= max(al.depleted.width, 0.1)
    assert abs(al.enriched.log_rho - np.log(rho_hi)) <= max(al.enriched.width, 0.1)
    assert al.span_R == pytest.approx(rho_hi / rho_lo, rel=0.5)


def test_the_basins_PARTITION_the_density():
    counts, lengths, sig, *_ = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    assert sum(m.basin_mass for m in al.modes) == pytest.approx(1.0, rel=0, abs=1e-9)
    # basins tile the grid: each starts where the previous ended
    for a, b in zip(al.modes, al.modes[1:], strict=False):
        assert a.hi == pytest.approx(b.lo, rel=0, abs=1e-12)


def test_w_separates_the_two_populations():
    counts, lengths, sig, *_ = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    n_lo = 600
    w = al.w_slot
    assert np.nanmax(w[:n_lo]) < 0.5, "a depleted region reads enriched"
    assert np.nanmin(w[n_lo:]) > 0.5, "an enriched region reads depleted"


def test_the_anchor_agrees_with_the_depleted_mode_and_the_flag_can_FAIL():
    """Getting this fixture right depends on the census's own robustness: any COHERENT anchor pool,
    large or tiny, drags a local maximum along with it through its own kernels, and "the basin
    containing the anchor" then follows it — so a coherent shift can NEVER flip the flag. What the
    flag guards is an INCOHERENT anchor population: anchors whose pooled rate is unrepresentative of
    any of them, landing in a basin whose peak is far away. That is also the honest real-data failure
    — an intergenic pool contaminated in both directions — so the fixture states it directly."""
    counts, lengths, sig, *_ = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    assert al.anchor_consistent
    assert al.anchor_gap_nats <= max(al.depleted.width, 0.1)

    # heterogeneous anchors: 3 at rho 0.001 and 3 at rho 0.5 against an 800-region bulk at 0.02.
    # The POOLED anchor rate (~0.25, log -1.38) sits where no peak is: gap 2.52 nats vs a depleted
    # width of 0.24 — the flag must flip, and the gap must be reported.
    rng = np.random.default_rng(3)
    n_bulk = 800
    lengths2 = np.full(n_bulk + 6, 5_000)
    counts2 = np.concatenate(
        [
            rng.poisson(0.001 * 5_000, 3),
            rng.poisson(0.5 * 5_000, 3),
            rng.poisson(0.02 * 5_000, n_bulk),
        ]
    )
    sig2 = np.concatenate([np.zeros(6, np.uint8), np.full(n_bulk, BIT_INTRON_POS, np.uint8)])
    sub2, ra2, mask2 = parts(counts2, lengths2, sig2)
    al2 = fit_abundance_landscape(sub2, ra2, mask2)
    assert not al2.anchor_consistent
    assert al2.anchor_gap_nats > al2.depleted.width


# ---------------------------------------------------------------------------
# the unimodal fixture — the capture-OFF shape
# ---------------------------------------------------------------------------


def unimodal_parts():
    rng = np.random.default_rng(11)
    n = 800
    lengths = np.full(n, 5_000)
    counts = rng.poisson(0.02 * 5_000, n)
    sig = np.concatenate([np.zeros(n // 2, np.uint8), np.full(n - n // 2, BIT_EXON_POS, np.uint8)])
    return counts, lengths, sig


def test_a_uniform_field_is_UNIMODAL_with_span_exactly_one_and_w_zero():
    counts, lengths, sig = unimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    assert al.enriched is None
    assert al.span_R == 1.0
    live = ~np.isnan(al.w_slot)
    assert live.any()
    assert np.all(al.w_slot[live] == 0.0)


# ---------------------------------------------------------------------------
# the mask, the selection, and the degenerate inputs
# ---------------------------------------------------------------------------


def test_a_NOT_model_free_region_contributes_NOTHING_and_reads_NaN():
    """A double-walled region must neither train the fit nor receive a responsibility. The perturbing
    half: hand the excluded region an enormous count — if it leaks into the fit, the census moves."""
    counts, lengths, sig = unimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    base = fit_abundance_landscape(sub, ra, mask)

    counts2 = counts.copy()
    counts2[0] = 10_000_000  # would manufacture an enriched mode if admitted
    excl = np.ones(counts.shape[0], dtype=bool)
    excl[0] = False
    sub2, ra2, mask2 = parts(counts2, lengths, sig, start_exact=excl, end_exact=excl)
    al = fit_abundance_landscape(sub2, ra2, mask2)
    assert np.isnan(al.w_slot[0])
    assert al.enriched is None, "an excluded region's count leaked into the fit"
    assert al.n_train == base.n_train - 1


def test_the_fit_is_DETERMINISTIC():
    counts, lengths, sig, *_ = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    a = fit_abundance_landscape(sub, ra, mask)
    b = fit_abundance_landscape(sub, ra, mask)
    np.testing.assert_array_equal(a.landscape.logP, b.landscape.logP)
    np.testing.assert_array_equal(a.w_slot, b.w_slot)


def test_too_few_training_regions_returns_None():
    sub, ra, mask = parts([5], [1000], [0])
    assert fit_abundance_landscape(sub, ra, mask) is None


def test_a_toy_with_NO_anchor_pool_falls_back_to_the_largest_basin():
    """All-exon annotation (a toy): no intergenic region exists, so the anchor rate is undefined and
    the depleted mode falls back to the largest-mass basin — reported, not asserted consistent."""
    counts, lengths, sig = unimodal_parts()
    sig = np.full_like(sig, BIT_EXON_POS)
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    assert al is not None
    assert np.isnan(al.anchor_log_rho)
    assert al.depleted.basin_mass == max(m.basin_mass for m in al.modes)


def test_the_START_and_END_banks_are_both_consumed_where_both_walls_clear():
    """With both sides exact the training count is the side-selected pair's pooled rate — (S+E)/2 per
    region. A fit reading S alone shifts every centre by the fixture's uneven strand split; pin the
    pooled count through n_train's companion: the depleted mode of a constant field must sit at
    Σ(S+E)/Σ(2ℓ), which the fixture states by hand."""
    n = 400
    lengths = np.full(n, 10_000)
    counts = np.full(n, 200)  # exact constant field: rho = 0.02
    sig = np.zeros(n, np.uint8)
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    assert abs(al.depleted.log_rho - np.log(0.02)) <= max(al.depleted.width, 0.05)


# ---------------------------------------------------------------------------
# gates added because the perturbation pass found holes — each names the perturbation it catches
# ---------------------------------------------------------------------------


def test_the_fit_is_EXACTLY_fit_landscape_on_the_selected_pair_with_var_zero():
    """Pins every argument of the underlying fit at once: the side-selected (counts, exposure) pair,
    mass ≡ count (the total's own ceiling), var ≡ 0 (a direct measurement has no deconvolution
    ambiguity) and the zero-count anchor rule. PERTURBATION: a var-perturbation that fabricates
    ambiguity moves NO other gate here — the fixtures are too well separated for a uniform
    down-weighting to move a mode — and byte-equality with an independent call is immune to that."""
    from rigel.calibration.landscape import fit_landscape
    from rigel.calibration.total_abundance import region_counts_and_exposure

    counts, lengths, sig, *_ = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    c, e, free = region_counts_and_exposure(sub, ra, mask)
    sel = free & (e > 0)
    ref = fit_landscape(c[sel], c[sel], e[sel], np.zeros(int(sel.sum())), anchor=(c[sel] == 0.0))
    np.testing.assert_array_equal(al.landscape.logP, ref.logP)
    np.testing.assert_array_equal(al.landscape.log_rho, ref.log_rho)


def test_w_matches_an_INDEPENDENT_posterior_recomputation():
    """The documented formula, re-implemented with scipy's own Poisson pmf (no shared helper): the
    region's kernel TIMES the fitted density, normalised on the grid, integrated over the enriched
    basin — asserted TIGHT, so any formula drift fires.

    PERTURBATION: dropping the landscape factor and using the kernel alone cannot be made to fail
    this oracle, and that is a property of the census rather than a hole. For every TRAINED region
    the two agree to floating point on every fixture that can be built, valley straddlers included,
    because the census partitions at density MINIMA and a training region's own kernel raises the
    density at its centre — so the cuts avoid it and its kernel mass stays inside one basin, and even
    two isolated wide-kernel regions form their own micro-basin rather than straddle. The landscape
    factor in `w` is therefore a formula commitment that matters for a FUTURE non-training query, not
    a behavioural difference on the training population."""
    from scipy.stats import poisson

    from rigel.calibration.total_abundance import region_counts_and_exposure

    counts, lengths, sig, *_ = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    c, e, free = region_counts_and_exposure(sub, ra, mask)
    sel = free & (e > 0)
    lam = np.exp(al.landscape.log_rho)[None, :] * e[sel][:, None]
    kern = poisson.pmf(np.round(c[sel])[:, None], lam)
    kern /= np.maximum(kern.sum(axis=1, keepdims=True), 1e-300)
    post = kern * np.exp(al.landscape.logP)[None, :]
    post /= np.maximum(post.sum(axis=1, keepdims=True), 1e-300)
    basin = (al.landscape.log_rho >= al.enriched.lo) & (al.landscape.log_rho <= al.enriched.hi)
    w_ref = post[:, basin].sum(axis=1)
    np.testing.assert_allclose(al.w_slot[sel], w_ref, rtol=0, atol=1e-9)


def test_the_anchor_picks_depleted_and_enriched_stays_ABOVE_it_with_three_modes():
    """THREE modes, anchors at the MIDDLE one, the LARGEST basin at the bottom: a depleted-by-mass
    rule picks the bottom (wrong), and an enriched-anywhere rule picks the bottom too (span < 1).
    PERTURBATION: on a two-mode fixture the anchor basin IS the largest and everything above the
    depleted basin IS the enriched mode, so neither of those two perturbations can fire there — this
    third mode is what gives them something to break."""
    rng = np.random.default_rng(5)
    n_low, n_mid, n_hi = 900, 200, 120
    lengths = np.concatenate([np.full(n_low, 20_000), np.full(n_mid, 10_000), np.full(n_hi, 1_000)])
    counts = np.concatenate(
        [
            rng.poisson(0.0001 * 20_000, n_low),  # the mass-dominant silent bulk
            rng.poisson(0.02 * 10_000, n_mid),  # the anchored depleted level
            rng.poisson(3.0 * 1_000, n_hi),  # the enriched minority
        ]
    )
    sig = np.concatenate(
        [
            np.full(n_low, BIT_INTRON_POS, np.uint8),
            np.zeros(n_mid, np.uint8),  # intergenic — the anchors, at the MIDDLE mode
            np.full(n_hi, BIT_EXON_POS, np.uint8),
        ]
    )
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    assert al.anchor_consistent
    assert abs(al.depleted.log_rho - np.log(0.02)) <= max(al.depleted.width, 0.1)
    assert al.enriched is not None
    assert al.enriched.log_rho > al.depleted.log_rho
    assert al.span_R > 1.0
    assert al.span_R == pytest.approx(3.0 / 0.02, rel=0.5)


# ---------------------------------------------------------------------------
# the calibrate() integration — no wall inputs fits nothing, the wall inputs fit, an injection is kept
# ---------------------------------------------------------------------------


def _calibrate_parts():
    import sys as _sys

    _sys.path.insert(0, "tests/calibration")
    from _synthetic import (
        make_gdna_fl_pmf,
        make_strand_models,
        make_synthetic_payload,
        make_synthetic_sj,
    )

    payload, ra = make_synthetic_payload()
    import dataclasses

    hist = np.zeros(201, dtype=np.uint32)
    hist[60] = 36
    payload = dataclasses.replace(payload, deposited_lengths=hist)
    return payload, ra, make_strand_models(0.95, 40), make_gdna_fl_pmf(), make_synthetic_sj()


def test_without_the_wall_inputs_the_landscape_is_SKIPPED_LOUDLY_and_nothing_raises(caplog):
    """The landscape is fit whenever the wall inputs (``mature_walls``, ``boundary_reach``) are passed,
    so a caller that supplies none gets no panel rather than an exception — the landscape feeds only the
    QC report's density panel, and a toy or a unit fixture that never wanted one runs without it.

    Nothing fits on unmasked totals either way: the choice is between refusing and not fitting, and
    the mask's wall bias is excluded under both. The skip is not a silent fallback — it is logged at
    WARNING and the object is ``None`` rather than a quietly different estimate: this object is read by
    the report and the debug bundle and by nothing in the solve."""
    import logging

    from rigel.calibration import calibrate
    from rigel.config import CalibrationConfig

    payload, ra, sm, pmf, sj = _calibrate_parts()
    d: dict = {}
    with caplog.at_level(logging.WARNING):
        calibrate(
            payload=payload,
            region_arrays=ra,
            strand_model=sm,
            gdna_fl_pmf=pmf,
            rna_fl_pmf=pmf,
            config=CalibrationConfig(),
            sj=sj,
            _debug=d,
        )
    assert d["abundance_landscape"] is None
    assert any(
        "wall inputs" in r.message % r.args if r.args else "wall inputs" in r.message
        for r in caplog.records
    ), "the skip must be LOUD — no silent no-op"


def test_the_flag_FITS_a_landscape_and_the_default_fits_NOTHING():
    """There is no config flag: passing the wall inputs is what fits a landscape, and omitting them
    fits nothing — both in the debug capture and in the priors bundle."""
    from rigel.calibration import calibrate
    from rigel.calibration.splice_graph import MatureWallDistances
    from rigel.config import CalibrationConfig

    payload, ra, sm, pmf, sj = _calibrate_parts()
    walls = MatureWallDistances(
        d_low=np.zeros((3, 2)), d_high=np.zeros((3, 2)), covered=np.zeros((3, 2), dtype=bool)
    )
    reach = (np.full((2, 2), 5000.0), np.full((2, 2), 5000.0))

    d_off: dict = {}
    calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=sm,
        gdna_fl_pmf=pmf,
        rna_fl_pmf=pmf,
        config=CalibrationConfig(),
        sj=sj,
        _debug=d_off,
    )
    assert "abundance_landscape" not in d_off or d_off.get("abundance_landscape") is None
    assert d_off["calibration_priors"].abundance_landscape is None

    d_on: dict = {}
    calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=sm,
        gdna_fl_pmf=pmf,
        rna_fl_pmf=pmf,
        config=CalibrationConfig(),
        sj=sj,
        mature_walls=walls,
        boundary_reach=reach,
        _debug=d_on,
    )
    al = d_on["abundance_landscape"]
    assert al is not None
    assert al.n_train >= 2
    assert d_on["calibration_priors"].abundance_landscape is al


def test_an_INJECTED_landscape_is_taken_verbatim_and_nothing_is_refit():
    """The toy path: a population-fitted landscape rides the priors bundle into a toy calibrate. The
    injected object must come back IDENTICALLY (`is`), and no wall inputs are needed because nothing
    is fit."""
    from rigel.calibration import calibrate
    from rigel.calibration.calibrate import InjectedCalibrationPriors
    from rigel.config import CalibrationConfig

    payload, ra, sm, pmf, sj = _calibrate_parts()
    counts, lengths, sig, *_ = bimodal_parts()
    sub2, ra2, mask2 = parts(counts, lengths, sig)
    donor = fit_abundance_landscape(sub2, ra2, mask2)
    assert donor is not None

    d: dict = {}
    calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=sm,
        gdna_fl_pmf=pmf,
        rna_fl_pmf=pmf,
        config=CalibrationConfig(),
        sj=sj,
        injected_priors=InjectedCalibrationPriors(abundance_landscape=donor),
        _debug=d,
    )
    assert d["abundance_landscape"] is donor


# ---------------------------------------------------------------------------
# the sweep handle and the exported basin-split rule
# ---------------------------------------------------------------------------


def test_knn_scale_passes_through_and_the_default_is_BIT_IDENTICAL():
    """`knn_scale` exists so a bandwidth sweep can be priced without restating the fit. The default
    must be bit-identical to omitting it (the shipped fit is not moved by the handle existing), and a
    different scale must MOVE the curve (an inert handle would let a sweep silently score one arm
    five times)."""
    counts, lengths, sig, *_ = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    base = fit_abundance_landscape(sub, ra, mask)
    explicit = fit_abundance_landscape(sub, ra, mask, knn_scale=0.5)
    assert np.array_equal(base.landscape.logP, explicit.landscape.logP)
    assert np.array_equal(base.landscape.log_rho, explicit.landscape.log_rho)

    wide = fit_abundance_landscape(sub, ra, mask, knn_scale=4.0)
    assert not np.array_equal(base.landscape.logP, wide.landscape.logP)


def test_split_basins_is_the_shipped_rule_importable_on_its_own():
    """The depleted/enriched selection is one rule with one home: `split_basins` must reproduce what
    `fit_abundance_landscape` publishes (same objects, not similar ones), pick the anchor-containing
    basin as depleted, and fall back to the largest-mass basin when the anchor is NaN — so an
    instrument scoring ANOTHER estimator's curve can apply the identical rule instead of restating
    it."""
    from rigel.calibration.abundance_landscape import AbundanceMode, split_basins

    counts, lengths, sig, *_ = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    dep, enr = split_basins(al.modes, al.anchor_log_rho)
    assert dep is al.depleted and enr is al.enriched

    lo = AbundanceMode(log_rho=-4.0, basin_mass=0.6, width=0.3, lo=-6.0, hi=-2.0)
    hi = AbundanceMode(log_rho=1.0, basin_mass=0.4, width=0.3, lo=-2.0, hi=3.0)
    dep2, enr2 = split_basins((lo, hi), -4.5)
    assert dep2 is lo and enr2 is hi
    # anchor NaN -> the LARGEST basin is depleted; nothing sits above it here, so enriched is None
    dep3, enr3 = split_basins((lo, hi), float("nan"))
    assert dep3 is lo
    assert enr3 is hi  # hi sits strictly above lo's basin, so it is still the enriched candidate
    # an anchor INSIDE the upper basin makes IT depleted and leaves nothing above
    dep4, enr4 = split_basins((lo, hi), 0.5)
    assert dep4 is hi and enr4 is None
    # TWO basins above with unequal masses: enriched must be the LARGER one. PERTURBATION: the
    # self-consistency assertions above cannot see a rule change that applies to both sides at once —
    # a min/max flip passes every one of them, and only this case catches it.
    mid = AbundanceMode(log_rho=0.0, basin_mass=0.30, width=0.2, lo=-2.0, hi=1.5)
    top = AbundanceMode(log_rho=2.5, basin_mass=0.10, width=0.2, lo=1.5, hi=4.0)
    base = AbundanceMode(log_rho=-4.0, basin_mass=0.60, width=0.3, lo=-6.0, hi=-2.0)
    dep5, enr5 = split_basins((base, mid, top), -4.5)
    assert dep5 is base and enr5 is mid


# ── the located enriched mode: the ruler's reference ────────────────────────────────────────────────


def _hand_landscape(peaks, sds, masses, members, walls=(), spread=0.1, lo=-16.0, hi=3.0, n=400):
    """A DensityLandscape rendered by hand — a mixture of Gaussians in natural-log density — with its
    kernels stated outright: per peak, ``members`` LOCATED kernels spread evenly across ±``spread`` nat
    of it, and ``walls`` as ``(centre, n)`` pairs of location-free kernels (anchors' resolution walls) —
    so a test controls where each basin sits, how much it holds, and how many located kernels stand
    behind it."""
    from rigel.calibration.landscape import DensityLandscape

    x = np.linspace(lo, hi, n)
    p = np.zeros_like(x)
    for mu, sd, m in zip(peaks, sds, masses, strict=True):
        p += m * np.exp(-0.5 * ((x - mu) / sd) ** 2) / sd
    p /= p.sum()
    centre = [
        np.linspace(mu - spread, mu + spread, k) for mu, k in zip(peaks, members, strict=True)
    ]
    located = [np.ones(k, dtype=bool) for k in members]
    for mu, k in walls:
        centre.append(np.linspace(mu - 0.1, mu + 0.1, k))
        located.append(np.zeros(k, dtype=bool))
    return DensityLandscape(
        log_rho=x,
        logP=np.log(p + 1e-300),
        n_train=int(sum(members) + sum(k for _, k in walls)),
        centre=np.concatenate(centre),
        located=np.concatenate(located),
    )


def test_a_unimodal_landscape_names_NO_enriched_mode():
    """Capture-OFF: one basin, nothing above the depleted mode, no reference, no contraction."""
    from rigel.calibration.abundance_landscape import located_enriched_mode

    assert located_enriched_mode(_hand_landscape([-3.0], [0.3], [1000], [400])) is None


def test_a_bimodal_landscape_names_the_LOCATED_upper_mode_at_its_peak_with_its_members():
    """Capture-ON: the basin above the depleted one holding the most located kernels, resolved by them
    (60 members against k = √300 ≈ 17); the reference is its peak and the members are published."""
    from rigel.calibration.abundance_landscape import located_enriched_mode

    ref = located_enriched_mode(_hand_landscape([-7.0, 0.5], [0.4, 0.25], [700, 300], [240, 60]))
    assert ref is not None and abs(ref.mode.log_rho - 0.5) < 0.05
    assert ref.n_members == 60


def test_a_basin_of_fewer_than_k_located_kernels_is_not_a_mode_however_narrow():
    """The blank contig's shadow transcription, a few false-positive fragments on slivers of support, or a
    sparse library's tail of a dozen measured exons: a basin whose located members number at most
    k = √n_located has no k-th neighbour inside itself and is no mode, whatever the density's cut made
    of it. PERTURBATION: the same landscape with the cluster grown past k IS a mode."""
    from rigel.calibration.abundance_landscape import located_enriched_mode

    # 400 located kernels in all ⇒ k = 20; 20 members do not resolve the basin, 21 do
    assert (
        located_enriched_mode(_hand_landscape([-14.0, -5.0], [0.5, 0.2], [998, 2], [380, 20]))
        is None
    )
    grown = located_enriched_mode(_hand_landscape([-14.0, -5.0], [0.5, 0.2], [998, 2], [379, 21]))
    assert grown is not None and abs(grown.mode.log_rho - (-5.0)) < 0.05 and grown.n_members == 21
    # and the members are read at the POPULATION's k, not their own √21: the same 21 strewn across
    # ±1.5 nat put every member's 20th neighbour ~3 nat away, and that is no location either
    strewn = _hand_landscape([-14.0, -5.0], [0.5, 2.0], [998, 2], [379, 21], spread=1.5)
    assert located_enriched_mode(strewn) is None


def test_anchors_walls_are_not_members_and_cannot_locate_a_basin():
    """`ISSUES: the-ruler-reference-on-sparse-real-libraries`, the reader's half: a basin above the bulk
    packed with 20,000 anchors' resolution walls around twelve located kernels names no reference — the
    walls are not members, and twelve is below k = √n_located. PERTURBATION: admitting every centre as a
    member, walls included, reads the walls as a located mode of 20,012 members."""
    from rigel.calibration.abundance_landscape import located_enriched_mode

    ls = _hand_landscape([-7.0, -1.5], [0.4, 0.2], [700, 300], [1000, 12], walls=[(-1.5, 20_000)])
    assert located_enriched_mode(ls) is None
    # and with the same twelve grown into a population the walls change nothing either way
    ls2 = _hand_landscape([-7.0, -1.5], [0.4, 0.2], [700, 300], [1000, 200], walls=[(-1.5, 20_000)])
    ref = located_enriched_mode(ls2)
    assert ref is not None and ref.n_members == 200


def test_the_enriched_mode_is_chosen_by_LOCATED_MEMBERS_above_the_depleted_one():
    """Two basins above the depleted one: the reference is the one holding the located population, not
    the one holding the most rendered mass — here the upper basin carries more density (walls render
    mass too) but thirty located kernels against the middle basin's 250, so the middle one is the
    reference. PERTURBATION: chosen by rendered mass, the upper basin is the candidate, and with thirty
    members against k = 31 it is no mode at all. A basin with no located member is nothing."""
    from rigel.calibration.abundance_landscape import located_enriched_mode

    ls = _hand_landscape([-7.0, -1.0, 0.8], [0.4, 0.3, 0.1], [700, 100, 250], [700, 250, 30])
    ref = located_enriched_mode(ls)
    assert ref is not None and abs(ref.mode.log_rho - (-1.0)) < 0.05 and ref.n_members == 250
    empty = _hand_landscape([-7.0, 0.5], [0.4, 0.25], [700, 300], [700, 0])
    assert located_enriched_mode(empty) is None


def _walls_around_a_tail(seed=0, n_long=20_000, n_dep=200, n_short=600, n_tail=12):
    """A sparse human-like library, refitted through the loop's own E-step: the depleted bulk is 20,000
    long anchors (walls at 10^-5) with 200 located intron pieces at the depleted level; 600 SHORT anchors —
    20–60 bp intergenic slivers that sequenced nothing — whose resolution walls 1/E sit between 10^-1.8
    and 10^-1.3; and twelve located exon kernels near 10^-1.3, a tail and not a population."""
    from rigel.calibration.landscape import fit_landscape

    rng = np.random.default_rng(seed)
    eff = np.concatenate(
        [
            np.full(n_long, 1e5),
            np.full(n_dep, 1e5),
            rng.uniform(20.0, 60.0, n_short),
            np.full(n_tail, 300.0),
        ]
    )
    count = np.concatenate(
        [
            np.zeros(n_long),
            rng.poisson(3e-5 * 1e5, n_dep) + 1.0,
            np.zeros(n_short),
            rng.poisson(0.05 * 300.0, n_tail) + 1.0,
        ]
    )
    mass = np.maximum(count * 2.0, 1.0)
    var = np.where(count > 0, 0.1, np.inf)
    anchor = count <= 0.0
    ls = None
    for _ in range(
        3
    ):  # calib_refit_iters: the E-step moves the walls' mass to the bulk, not their centres
        ls = fit_landscape(count, mass, eff, var, anchor=anchor, prev=ls)
    return ls


def test_a_fitted_basin_of_walls_around_a_handful_of_located_kernels_is_NOT_a_mode():
    """LBX0190's defect through the estimator itself (the falsification test, verified failing on the
    shipped reader, which read a located mode of 612 members here): the walls' centres pack the basin
    above the bulk so densely that every rendered width there sits at the grid step. A wall is not a
    location: the basin's located members are twelve, fewer than k = √212, and it names no reference."""
    from rigel.calibration.abundance_landscape import located_enriched_mode

    ls = _walls_around_a_tail()
    assert ls is not None
    assert int(ls.located.sum()) == 212
    assert located_enriched_mode(ls) is None
