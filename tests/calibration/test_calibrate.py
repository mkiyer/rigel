"""calibrate(): the acyclic single-pass calibrator — schema + structural invariants.

These pin the *mechanics* (a valid result, conservation on each axis, bounded masses, sane exposure,
κ_rna provenance, the confidence knob). Converged *biology* (paralog rescue, exon→RNA) needs realistic
data and is covered by the scenario suite.

⭐ **Every expected number here is arithmetic on ``_synthetic.make_synthetic_payload``'s banks**, not a
recorded output — so a change in the solver moves the composition but never the totals, and a change in
the AXES fails immediately. The fixture's three axes are deliberately three different lengths
(3 regions / 2 boundaries / 1 sj).
"""

from __future__ import annotations

import numpy as np
import pytest
from _synthetic import (
    make_gdna_fl_pmf,
    make_strand_models,
    make_synthetic_sj,
    make_synthetic_payload,
)

from rigel.calibration import calibrate
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig

# region_contained_count summed over the two genome-strand columns: [10+2, 1+20, 7+8]
REGION_TOTAL = np.array([12.0, 21.0, 15.0])
# boundary_unspliced + boundary_spliced: [4+1, 2+3] + [0+0, 6+0]
BOUNDARY_TOTAL = np.array([5.0, 11.0])
BOUNDARY_SPLICED = np.array([0.0, 6.0])
SJ_FLUX = np.array([13.0])  # sj_count [9, 4]


_UNSET = object()


def _run(config=None, sj=_UNSET):
    payload, ra = make_synthetic_payload()
    return calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=make_strand_models(0.95, 40),
        gdna_fl_pmf=make_gdna_fl_pmf(),
        rna_fl_pmf=make_gdna_fl_pmf(),  # a valid pmf; these tests pin mechanics, not the splice fraction
        config=config or CalibrationConfig(),
        sj=make_synthetic_sj() if sj is _UNSET else sj,
    )


def test_returns_valid_result_on_all_three_axes():
    result = _run()
    assert isinstance(result, CalibrationResult)
    assert (result.n_regions, result.n_boundaries, result.n_sj) == (3, 2, 1)


def test_mass_conserved_per_region():
    """A region's gDNA + RNA equals its contained count. There is no spliced term to add: the
    accumulator credits ``region_contained`` only when the fragment used no sj."""
    result = _run()
    np.testing.assert_allclose(result.mass_gdna_region + result.mass_rna_region, REGION_TOTAL)


def test_mass_conserved_per_boundary_INCLUDING_the_spliced_crossings():
    """A boundary's gDNA + RNA equals unspliced + spliced. ``mass_rna_boundary`` is spliced-INCLUSIVE — a
    certified-RNA crossing is RNA whatever the unspliced mixture resolves to, since gDNA cannot splice
    — so dropping it here would lose 6 of boundary 1's 11 fragments."""
    result = _run()
    np.testing.assert_allclose(result.mass_gdna_boundary + result.mass_rna_boundary, BOUNDARY_TOTAL)
    np.testing.assert_allclose(result.mass_rna_spliced_boundary, BOUNDARY_SPLICED)
    assert np.all(result.mass_rna_boundary >= result.mass_rna_spliced_boundary - 1e-9)


def test_sj_flux_is_exported_VERBATIM_and_never_deconvolved():
    """⭐ The third axis (owner ruling, 2026-07-30). A sj boundary is pure mature RNA by
    construction, so there is nothing to split: the result carries ``sj_count`` summed over the
    genome-strand columns, exactly.

    ⚠ It is two orders of magnitude away from ``mass_rna_spliced_boundary`` on real data at the same boundary
    — 13 vs 0/6 even in this toy — which is why folding the two into one "mature" number names
    nothing.
    """
    result = _run()
    np.testing.assert_array_equal(result.count_rna_sj, SJ_FLUX)


def test_the_conserved_sj_mass_recovers_the_ACCUMULATORS_OWN_sj_mass_BANK():
    """⭐⭐ ``sj_conserved_mass`` is ``sj_count × (sj_mass / sj_count)``, so it must come back as
    ``sj_mass`` itself — the bank the scanner wrote, not an approximation of it.

    ⛔ **This is the gate that says the published quantity is the CONSERVED one.** The fixture's
    sj carries 13 incidences and ``sj_mass`` 1.3, a **10×** gap, so the two cannot be confused by
    coincidence — which they could on real data at a boundary where every fragment used one sj. On
    ``g00 ss0.99 capture_off`` the same round trip agrees with the bank to 9.1e-13.
    """
    payload, _ = make_synthetic_payload()
    result = _run()
    # ⭐ SUMMED over the strand columns: `sj_mass` went per strand on 2026-08-13 and `substrate` folds
    # it, because the incidence→fragment conversion has no strand in it. This gate therefore now pins
    # the FOLD as well as the conversion — and the fixture's columns are unequal (0.9 / 0.4), so a fold
    # that took one column or their mean cannot pass.
    np.testing.assert_allclose(
        result.sj_conserved_mass, payload.sj_mass.sum(axis=1), rtol=0, atol=1e-12
    )
    assert payload.sj_mass[0, 0] != payload.sj_mass[0, 1], "the fixture cannot separate the fold rules"
    # ⚠ Could this have failed? The incidence is 10× the mass here, so passing by accident is not
    # available (`TRAPS: could-the-arm-have-fired`).
    assert not np.allclose(payload.sj_mass.sum(axis=1), SJ_FLUX)


def test_masses_bounded_by_their_own_totals():
    result = _run()
    for g, tot in (
        (result.mass_gdna_region, REGION_TOTAL),
        (result.mass_gdna_boundary, BOUNDARY_TOTAL),
    ):
        assert np.all(g >= -1e-9)
        assert np.all(g <= tot + 1e-9)


def test_an_intergenic_region_is_ALL_gDNA():
    """Region 2 carries no exon or intron bit, so no RNA can be contained in it — a structural lock,
    not an inference. ⚠ ``mass_rna_region[2] == 0`` exactly; a floored or smoothed answer here would be
    manufacturing RNA where the annotation says none exists."""
    result = _run()
    assert result.mass_rna_region[2] == 0.0
    assert result.mass_gdna_region[2] == REGION_TOTAL[2]


# --- the two geometric supports ---------------------------------------------------------------


def test_the_supports_are_the_TWO_FRAMES_of_one_formula_family():
    """⭐ Arithmetic, not a recorded number. Every region is 100 bp and the gDNA pmf is a delta at 50:

        contained  E_f[(100 − 50 + 1)+]                    = 51   — starts that FIT inside
        crossing   E_f[min(w−1, R_lo, R_hi, ...)] at R = ∞  = 49   — offsets that SPAN the boundary

    The two differ by 2 and neither is the region length, which is the whole point: ``region_size_bp``
    ignores the fit-inside constraint and ``mean_FL`` ignores the ``−1``.
    """
    result = _run()
    np.testing.assert_allclose(result.gdna_region_eff_len, [51.0, 51.0, 51.0])
    np.testing.assert_allclose(result.gdna_boundary_eff_len, [49.0, 49.0])


def test_gdna_density_global_is_a_ratio_of_SUMS_over_both_axes():
    """Σ gDNA mass / Σ gDNA support, pooled across regions AND boundaries — never a mean of per-object
    ratios, which is a different number whenever the supports differ."""
    result = _run()
    expected = (result.mass_gdna_region.sum() + result.mass_gdna_boundary.sum()) / (
        result.gdna_region_eff_len.sum() + result.gdna_boundary_eff_len.sum()
    )
    assert result.gdna_density_global == pytest.approx(expected)


# --- library scalars --------------------------------------------------------------------------


def test_gdna_strand_overdispersion_populated():
    assert 0.0 <= _run().gdna_strand_overdispersion < 1.0


def test_rna_strand_overdispersion_populated():
    # clamped to the Beta(2,2) ceiling (0.2)
    assert 0.0 <= _run().rna_strand_overdispersion <= 0.2


def test_kappa_matches_strand_balance():
    # κ_rna is the posterior-predictive strand fit; the calibrator passes it through.
    from rigel.calibration.strand_balance import fit_strand_balance

    sb = fit_strand_balance(make_strand_models(0.95, 40))
    assert _run().rna_sense_frac == sb.rna_sense_frac


def test_density_and_supports_sane():
    result = _run()
    assert np.isfinite(result.gdna_density_global) and result.gdna_density_global >= 0.0
    for arr in (result.gdna_region_eff_len, result.gdna_boundary_eff_len):
        assert np.all(np.isfinite(arr)) and np.all(arr >= 0.0)
    assert 0.0 <= result.rna_sense_frac <= 1.0


# --- the sj axis is CHECKED, not assumed ------------------------------------------------


def test_a_mismatched_sj_axis_is_REFUSED():
    """⛔ A sj axis built against a different graph would place every splice on the wrong boundary,
    and nothing downstream would fault on it — the shape is plausible either way. Refuse at the door
    (476,719 of 476,732 real fragments once vanished inside a deposit
    while every golden test passed)."""
    with pytest.raises(ValueError, match="sj axis"):
        _run(sj=None)  # the payload declares one sj; an empty axis is not it


def test_the_sj_axis_length_is_what_is_checked_not_its_content():
    from rigel.calibration.splice_graph import SpliceJunctionGeometry
    from rigel.types import Strand

    two = SpliceJunctionGeometry(
        src_region=np.array([0, 0], dtype=np.int64),
        dst_region=np.array([2, 2], dtype=np.int64),
        strand=np.full(2, int(Strand.POS), dtype=np.int8),
        reach_lo=np.full(2, 100.0),
        reach_hi=np.full(2, 100.0),
    )
    with pytest.raises(ValueError, match="2 boundaries but the payload has 1"):
        _run(sj=two)


# --- the intron factory -----------------------------------------------------------------------


def test_intron_factory_flag_runs_and_conserves_mass():
    import dataclasses

    result = _run(dataclasses.replace(CalibrationConfig(), intron_factory=True))
    assert isinstance(result, CalibrationResult)
    np.testing.assert_allclose(result.mass_gdna_region + result.mass_rna_region, REGION_TOTAL)


def test_intron_factory_noop_without_introns():
    # Correct scoping: with no INTRON regions (the synthetic is +exon/−exon/intergenic) the factory is a
    # graceful no-op — byte-identical to the flag-off calibration.
    import dataclasses

    off = _run(CalibrationConfig())
    on = _run(dataclasses.replace(CalibrationConfig(), intron_factory=True))
    np.testing.assert_array_equal(off.mass_gdna_region, on.mass_gdna_region)
    np.testing.assert_array_equal(off.mass_rna_region, on.mass_rna_region)
