"""``assemble_priors`` — the acyclic ``CalibrationResult`` as a per-locus EM prior, and the locus
projection it stands on.

A region owns the fragments contained in it; a boundary owns the fragments that cross it. A locus
collects both — its regions by genomic overlap, its boundaries by touching those regions, so a
locus of ``k`` contiguous regions carries ``k + 1`` boundaries including its two outer ones — and
no boundary's mass is ever folded into a region's total. The second half of this file gates that
projection; the first half gates what ``assemble_priors`` builds on it, whose bedrock invariant is
that under a uniform gDNA field every object's ``min(m/ρ_ref, S)`` returns its own effective
support ``S``, so ``gdna_eff_len == span == ΣS`` exactly and an unenriched library contracts
nothing. Dividing by the genomic ``region_size_bp`` instead fabricates a contraction.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.priors import (
    _boundary_locus_shares,
    _region_locus_shares,
    assemble_priors,
    contended_boundaries,
)
from rigel.calibration.region_arrays import RegionArrays, boundary_region_indices
from rigel.calibration.result import CalibrationResult
from rigel.calibration.signature import BIT_EXON_POS
from rigel.config import CalibrationConfig
from rigel.locus import Locus, MultiLocus


def _result(
    *,
    region_g,
    region_r,
    region_eff,
    boundary_g=None,
    boundary_r=None,
    boundary_eff=None,
    boundary_spliced=None,
    mass_per_crossing=None,
    gdna_density_global=0.01,
    gdna_reference_density=None,
    gdna_reference_members=0,
    rna_region_eff=None,
    rna_boundary_eff=None,
    efficiency=None,
    efficiency_boundary=None,
) -> CalibrationResult:
    """Build a result on the three axes. One reference with ``n`` regions owns exactly ``n − 1`` boundaries.

    ``boundary_*`` default to zeros, so a caller that cares only about contained mass writes only the region
    arrays — but the boundary axis is still the RIGHT LENGTH, because a boundary axis inconsistent with its
    own region axis is a mis-shaped fixture, not a "no boundaries" one.

    The RNA supports default to the gDNA ones, so a test about projection, conservation or
    re-keying — not about the length tilt — keeps both components on one support and its ``g:r``
    ratio is unaffected by it. The tilt is exercised where it belongs, in `test_prior_units.py`, by
    giving the two components genuinely different opportunities.

    ``efficiency`` / ``efficiency_boundary`` are the per-object capture efficiencies (1 everywhere by
    default, as a field with no reference must carry).
    """
    ng = np.asarray(region_g, dtype=np.float64)
    n = ng.shape[0]
    ne = max(n - 1, 0)
    ez = np.zeros(ne, dtype=np.float64)
    region_eff_arr = np.asarray(region_eff, dtype=np.float64)
    boundary_eff_arr = (
        ez.copy() if boundary_eff is None else np.asarray(boundary_eff, dtype=np.float64)
    )
    return CalibrationResult(
        count_gdna_region=ng,
        count_rna_region=np.asarray(region_r, dtype=np.float64),
        count_gdna_boundary=ez.copy()
        if boundary_g is None
        else np.asarray(boundary_g, dtype=np.float64),
        count_rna_boundary=ez.copy()
        if boundary_r is None
        else np.asarray(boundary_r, dtype=np.float64),
        count_rna_spliced_boundary=(
            ez.copy()
            if boundary_spliced is None
            else np.asarray(boundary_spliced, dtype=np.float64)
        ),
        # Geometry, not a split. 1.0 is the identity — a boundary whose flanks both exceed every
        # fragment length, where one crossing IS one fragment. A test exercising K-inflation overrides it.
        boundary_mass_per_crossing=(
            np.ones_like(ez)
            if mass_per_crossing is None
            else np.asarray(mass_per_crossing, dtype=np.float64)
        ),
        count_rna_sj=np.zeros(0, dtype=np.float64),
        boundary_spliced_mass_per_crossing=np.ones_like(ez),
        sj_mass_per_crossing=np.ones(0, dtype=np.float64),
        gdna_region_eff_len=region_eff_arr,
        gdna_boundary_eff_len=boundary_eff_arr,
        rna_region_eff_len=(
            region_eff_arr
            if rna_region_eff is None
            else np.asarray(rna_region_eff, dtype=np.float64)
        ),
        rna_boundary_eff_len=(
            boundary_eff_arr
            if rna_boundary_eff is None
            else np.asarray(rna_boundary_eff, dtype=np.float64)
        ),
        gdna_frac_region=np.zeros_like(ng),
        rna_pos_frac_region=np.zeros_like(ng),
        rna_neg_frac_region=np.zeros_like(ng),
        gdna_frac_boundary=ez.copy(),
        rna_pos_frac_boundary=ez.copy(),
        rna_neg_frac_boundary=ez.copy(),
        gdna_density_global=gdna_density_global,
        gdna_reference_density=gdna_reference_density,
        gdna_reference_members=0 if gdna_reference_density is None else 1,
        gdna_capture_efficiency_region=(
            np.ones(n) if efficiency is None else np.asarray(efficiency, dtype=np.float64)
        ),
        gdna_capture_efficiency_boundary=(
            np.ones(ne)
            if efficiency_boundary is None
            else np.asarray(efficiency_boundary, dtype=np.float64)
        ),
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n,
        n_boundaries=ne,
        n_sj=0,
        config=CalibrationConfig(),
    )


def _uniform_field(region_eff, boundary_eff, rho) -> CalibrationResult:
    """A genuinely UNIFORM gDNA field: every object's mass is ``ρ × its own effective support``.

    That is the accumulator's deposition law stated directly — ``ρ·E_f[(L−w+1)+]`` contained,
    ``ρ·E_f[w−1]`` crossing — with no faces, no halving and nothing to cancel. The factor-1 invariant
    then says ``gdna_eff_len == span == Σ S`` exactly, and ``G/eff_len`` recovers the true ρ.
    """
    region_eff = np.asarray(region_eff, dtype=np.float64)
    boundary_eff = np.asarray(boundary_eff, dtype=np.float64)
    return _result(
        region_g=rho * region_eff,
        region_r=np.zeros_like(region_eff),
        region_eff=region_eff,
        boundary_g=rho * boundary_eff,
        boundary_eff=boundary_eff,
        gdna_density_global=rho,
    )


def _regions(starts, ends, signature=None) -> RegionArrays:
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    n = starts.shape[0]
    return RegionArrays(
        ref_id=np.zeros(n, dtype=np.int32),
        start=starts,
        end=ends,
        signature=(
            np.zeros(n, dtype=np.uint8) if signature is None else np.asarray(signature, np.uint8)
        ),
        strand_class=np.zeros(n, dtype=np.int8),
        region_size_bp=(ends - starts).astype(np.float64),
        ref_offsets=np.array([0, n], dtype=np.int32),
        n_refs=1,
    )


def _ml(locus_id, blocks) -> MultiLocus:
    loci = tuple(Locus(ref=str(rid), ref_id=rid, start=s, end=e) for rid, s, e in blocks)
    return MultiLocus(
        multi_locus_id=locus_id,
        transcript_indices=np.array([], dtype=np.int32),
        unit_indices=np.array([], dtype=np.int32),
        gdna_span=sum(e - s for _, s, e in blocks),
        loci=loci,
    )


# --- the bedrock invariant: factor = 1 under uniform gDNA -----------------------------------------


def test_factor_one_under_uniform_gdna():
    # THE correctness criterion. A uniform (unenriched) gDNA field carries no reference, every efficiency
    # is 1, and gdna_eff_len = span = Σ S EXACTLY: region_eff=[120,200,80] (region 1 is SHORT),
    # boundary_eff=[120,120] at q = 1, ρ=0.02 over 3 same-ref regions ⇒ span = 400 + 240 = 640; and the
    # gDNA per-position rate G/eff_len recovers the true ρ.
    region_eff = [120.0, 200.0, 80.0]
    boundary_eff = [120.0, 120.0]
    rho = 0.02
    span = sum(region_eff) + sum(boundary_eff)  # 640
    cal = _uniform_field(region_eff, boundary_eff, rho)
    ra = _regions([0, 120, 320], [120, 320, 400])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 400)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [span], rtol=1e-9)
    # The prior is a CONSERVED FRAGMENT COUNT read out of the bank: the contained mass ρ·Σregion_eff
    # (8.0) plus the crossing mass ρ·Σboundary_eff (4.8) rescaled by q, which this fixture sets to the
    # identity 1.0 — flanks exceeding every fragment length, where one crossing IS one fragment.
    # The wrong answer it must not equal is `ρ · span_bp` = 8.0, a density rule that reaches fragment
    # units by dividing the mass by its own opportunity and re-integrating: it drops the 4.8 of
    # crossing fragments entirely, because a 0-bp boundary contributes no genomic span to integrate.
    np.testing.assert_allclose(
        priors.gdna_prior_count, [rho * (sum(region_eff) + sum(boundary_eff))], rtol=1e-9
    )
    assert not np.isclose(priors.gdna_prior_count[0], rho * 400.0)  # not the ρ·span_bp density rule


def test_factor_one_holds_for_any_density():
    # The length does not depend on ρ at all: a 50000× denser uniform library reads the same span.
    region_eff = [300.0, 150.0]
    boundary_eff = [200.0]
    span = 650.0
    ra = _regions([0, 300], [300, 450])
    for rho in (1e-4, 0.01, 0.5, 5.0):
        priors = assemble_priors(
            _uniform_field(region_eff, boundary_eff, rho), ra, [_ml(0, [(0, 0, 450)])]
        )
        np.testing.assert_allclose(priors.gdna_eff_len, [span], rtol=1e-9)


def test_eff_len_uses_effective_support_not_genomic_size():
    # PROOF the length is the EFFECTIVE support gdna_region_eff_len, NOT the genomic region_size_bp:
    # genomic sizes 100+100+100 = 300, supports 120+200+80 = 400 plus the crossings 300 ⇒ 700. The count
    # is on these objects, so the length is too; a base count would read 300 and something else.
    region_eff = [120.0, 200.0, 80.0]
    boundary_eff = [150.0, 150.0]
    span = 700.0
    cal = _uniform_field(region_eff, boundary_eff, 0.03)
    ra = _regions([0, 100, 200], [100, 200, 300])  # genomic sizes 100,100,100 (Σ=300) ≠ region_eff
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [span], rtol=1e-9)
    assert not np.isclose(priors.gdna_eff_len[0], 300.0 + sum(boundary_eff))


def test_a_boundary_enters_the_length_at_its_crossing_support_CONVERTED_BY_q():
    # The count converts a boundary's mass by q, the conserved mass per crossing; the length converts the
    # boundary's support by the SAME q, so a start is counted once however many boundaries its fragment
    # crosses. Three 40-bp pieces, fragments of ~181 bp, q = ½ on both boundaries: every crossing start spans
    # a piece and is in both supports, so the length is 30 + ½·360 = 210 and the density 4.2/210 = 0.02,
    # the field's (region 0.2/10, boundary 3.6/180). The unconverted 30 + 360 = 390 read the same locus at
    # half its density; it was pinned here until 2026-09-20 on a capture-OFF transcript-number measurement.
    region_eff = [10.0, 10.0, 10.0]
    boundary_eff = [180.0, 180.0]
    cal = _result(
        region_g=[0.2, 0.2, 0.2],
        region_r=[0.0, 0.0, 0.0],
        region_eff=region_eff,
        boundary_g=[3.6, 3.6],
        boundary_eff=boundary_eff,
        mass_per_crossing=[0.5, 0.5],
    )
    ra = _regions([0, 40, 80], [40, 80, 120])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 120)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [30.0 + 0.5 * 360.0], rtol=1e-9)
    np.testing.assert_allclose(priors.gdna_prior_count / priors.gdna_eff_len, [0.02], rtol=1e-9)
    assert not np.isclose(priors.gdna_eff_len[0], 30.0 + 360.0)


@pytest.mark.parametrize(
    "lengths",
    [(60,), (50, 60, 70)],
    ids=["one-length", "three-lengths"],
)
def test_the_length_counts_each_crossing_start_once_as_the_count_counts_each_fragment_once(lengths):
    """THE CONSISTENCY THE TWO LEDGERS MUST SHARE, by brute-force enumeration against the deposit rule
    (the executable specification is ``tests/native/_accumulator_reference.py``: a crossing fragment
    deposits a count of +1 at EVERY boundary it crosses and a mass summing to 1 across them).

    The pseudocount converts a boundary's incidence count to fragments by ``q = mass / count``. The
    length must convert the boundary's incidence SUPPORT by the same ``q``, so that a start position is
    counted once however many boundaries its fragment crosses. Three 40-bp pieces inside a long reference
    with fragments longer than a piece: every crossing start that spans a piece would otherwise be
    counted at both of its boundaries, and the locus's gDNA would read at HALF the field's density.

    Under a uniform field of one fragment per start position the pseudocount is the number of distinct
    starts overlapping the locus, and so must the length be: their ratio is the field's density, 1.
    """
    from rigel.calibration.effective_length import (
        UNBOUNDED_REACH,
        contained_eff_length,
        crossing_eff_length,
    )

    starts = np.array([0, 400, 440, 480, 520])
    ends = np.array([400, 440, 480, 520, 1000])
    ra = _regions(starts, ends)
    lo, hi = boundary_region_indices(ra.ref_id)
    bpos = ends[lo]  # the boundary's coordinate: the end of its lower flank
    n_r, n_e = len(starts), len(bpos)
    pmf = np.zeros(max(lengths) + 1)
    pmf[list(lengths)] = 1.0 / len(lengths)
    S_r = contained_eff_length(ends - starts, pmf)
    S_e = crossing_eff_length(pmf, np.full(n_e, UNBOUNDED_REACH), np.full(n_e, UNBOUNDED_REACH))

    # the enumeration: one fragment per (start, length) at the pmf's weight, deposited by the rule
    count_r, count_e, mass_e = np.zeros(n_r), np.zeros(n_e), np.zeros(n_e)
    distinct_starts_in_locus = 0.0
    for w in lengths:
        wt = 1.0 / len(lengths)
        for s in range(0, 1000 - w + 1):
            e = s + w
            if s < 520 and e > 400:
                distinct_starts_in_locus += wt
            crossed = [j for j in range(n_e) if s < bpos[j] < e]
            if not crossed:
                r = int(np.searchsorted(ends, s, side="right"))
                assert starts[r] <= s and e <= ends[r]
                count_r[r] += wt
            else:
                for j in crossed:
                    count_e[j] += wt
                    mass_e[j] += wt / len(crossed)
    # every object reads its own density, 1 — the enumeration and the module's geometry agree exactly
    np.testing.assert_allclose(count_r, S_r, rtol=1e-12)
    np.testing.assert_allclose(count_e, S_e, rtol=1e-12)
    q = mass_e / count_e

    cal = _result(
        region_g=count_r,
        region_r=np.zeros(n_r),
        region_eff=S_r,
        boundary_g=count_e,
        boundary_eff=S_e,
        mass_per_crossing=q,
    )
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 400, 520)])])
    # the pseudocount already counts each fragment once
    np.testing.assert_allclose(priors.gdna_prior_count, [distinct_starts_in_locus], rtol=1e-12)
    # ...and so must the length: the density under a uniform field is the field's, 1
    np.testing.assert_allclose(priors.gdna_eff_len, [distinct_starts_in_locus], rtol=1e-12)
    np.testing.assert_allclose(priors.gdna_prior_count / priors.gdna_eff_len, [1.0], rtol=1e-12)
    # PERTURBATION: the incidence total counts a spanning start at both of its boundaries
    incidences = float(S_r[1:4].sum() + S_e.sum())
    assert incidences > distinct_starts_in_locus * 1.2
    assert not np.isclose(priors.gdna_eff_len[0], incidences)


def test_every_OBJECT_has_the_same_density_under_a_uniform_field():
    """The precondition for the per-object ``min()`` factor-1 identity, asserted on the objects
    themselves rather than on a folded total.

    Summing each boundary's mass into a flank region before dividing would check the fold's density
    and never a boundary's own. Regions and boundaries are peers, so each axis is checked on its own.
    """
    region_eff = np.array([120.0, 200.0, 80.0])
    boundary_eff = np.array([120.0, 120.0])
    rho = 0.02
    cal = _uniform_field(region_eff, boundary_eff, rho)
    np.testing.assert_allclose(cal.count_gdna_region / cal.gdna_region_eff_len, rho, rtol=1e-9)
    np.testing.assert_allclose(cal.count_gdna_boundary / cal.gdna_boundary_eff_len, rho, rtol=1e-9)


# --- mass / projection (independent of the support choice) ----------------------------------------


def test_single_locus_projects_both_components():
    # The priors are CONSERVED FRAGMENT COUNTS. No crossing mass here, so both are the contained
    # mass alone — one deposit per contained fragment, nothing to convert:
    #   gDNA: Σm = 4.5      (a ρ·span density rule would give 4.5/750 · 450 = 2.7)
    #   RNA : Σm = 12.0     (and 7.2)
    # The g:r RATIO is 0.375 either way, because a common divisor cancels from a ratio. That is
    # exactly why the ratio is NOT what discriminates the two rules — the totals are.
    cal = _result(
        region_g=[1.0, 2.0, 1.5],
        region_r=[3.0, 4.0, 5.0],
        region_eff=[100.0, 200.0, 150.0],
        boundary_eff=[150.0, 150.0],
    )
    ra = _regions([0, 100, 300], [100, 300, 450])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 450)])])
    np.testing.assert_allclose(priors.rna_prior_count, [12.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [4.5])
    np.testing.assert_allclose(
        priors.gdna_prior_count[0] / priors.rna_prior_count[0], 4.5 / 12.0, rtol=1e-9
    )
    span = (100.0 + 200.0 + 150.0) + 2 * 150.0
    assert 0.0 < priors.gdna_eff_len[0] <= span + 1e-9


def test_gdna_mass_conservation_regions_plus_boundaries():
    # CONSERVATION OF MASS. total gDNA = Σ region mass + Σ boundary mass, every object counted exactly once.
    #   regions = [2,3,1] (Σ=6); boundaries = [2,3] (Σ=5)  ⇒  total gDNA = 11.
    # Two axes summed once each: there are no per-face half-arrays and no terminal slots to zero by
    # hand, so there is no half of a boundary that can be counted twice or not at all.
    cal = _result(
        region_g=[2.0, 3.0, 1.0],
        region_r=[0.0, 0.0, 0.0],
        region_eff=[100.0, 100.0, 100.0],
        boundary_g=[2.0, 3.0],
        boundary_eff=[50.0, 50.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    # OBJECT conservation is what this test is named for: the locus covers every region, so it collects
    # every region AND every boundary, and the prior is their total — nothing dropped, nothing double-counted.
    # And here the prior equals the raw sum, 11.0, which is NOT evidence that it *is* a raw sum.
    # This fixture's q is the identity 1.0, so incidence and fragment coincide by construction and this
    # test CANNOT tell the two rules apart. The discrimination lives in the q ≠ 1 test below and in
    # `test_prior_units.py`; asserting 11.0 here would otherwise read as a ruling that it is a raw sum.
    np.testing.assert_allclose(priors.gdna_prior_count, [11.0])
    np.testing.assert_allclose(
        cal.boundary_mass_per_crossing, 1.0
    )  # ...the reason it coincides, pinned
    np.testing.assert_allclose(
        priors.gdna_prior_count.sum(), cal.count_gdna_region.sum() + cal.count_gdna_boundary.sum()
    )
    assert contended_boundaries(ra, [_ml(0, [(0, 0, 300)])], 1).size == 0  # nothing double-claimed


def test_the_crossing_mass_is_rescaled_by_the_conserved_share():
    """The one test in this file that separates a conserved count from a raw incidence sum.

    Every other fixture here leaves ``boundary_mass_per_crossing`` at the identity 1.0, where one crossing
    IS one fragment and the two rules coincide — so they all pass under either. This one sets ``q`` to
    ``[0.5, 0.25]``: a fragment crossing boundary 0 deposited on 2 objects on average, boundary 1 on 4.

        gDNA = 3 (contained, one deposit each) + 4·0.5 + 8·0.25 = 7.0     raw sum would be 15.0
        RNA  = 6                               + 4·0.5 + 4·0.25 = 9.0     raw sum would be 14.0

    The CONTAINED term must not be rescaled — a contained fragment touches exactly one region and is
    already a count. Rescaling it too would give 3·? and is the other wrong answer this pins out.
    """
    cal = _result(
        region_g=[1.0, 1.0, 1.0],
        region_r=[2.0, 2.0, 2.0],
        region_eff=[100.0, 100.0, 100.0],
        boundary_g=[4.0, 8.0],
        boundary_r=[4.0, 4.0],
        boundary_eff=[50.0, 50.0],
        mass_per_crossing=[0.5, 0.25],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.gdna_prior_count, [7.0])
    np.testing.assert_allclose(priors.rna_prior_count, [9.0])
    assert not np.isclose(priors.gdna_prior_count[0], 15.0)  # not the raw incidence sum
    assert not np.isclose(priors.rna_prior_count[0], 14.0)


def test_spliced_mass_withheld_from_rna_prior():
    # A spliced fragment has no gDNA candidate in the EM (gDNA does not splice) → it is guaranteed-RNA
    # and assigned directly, so it must NOT load rna_prior_count. Region RNA [3,4,5] (Σ=12) plus boundary RNA
    # [4,4] of which [1,3] is spliced ⇒ RNA mass = 12 + (4−1) + (4−3) = 16 (NOT 20).
    # q is the identity here, so the conserved fragment count is that 16 unchanged; gDNA is its
    # contained 4.5 (no crossing mass). The WITHHOLDING is what this test pins: without it the RNA
    # mass would be 20 and the prior 20.0.
    cal = _result(
        region_g=[1.0, 2.0, 1.5],
        region_r=[3.0, 4.0, 5.0],
        region_eff=[100.0, 200.0, 150.0],
        boundary_r=[4.0, 4.0],
        boundary_spliced=[1.0, 3.0],
        boundary_eff=[150.0, 150.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.rna_prior_count, [16.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [4.5])
    assert priors.rna_prior_count[0] < 20.0  # the spliced mass really is withheld


def test_the_sj_flux_does_NOT_enter_the_rna_prior():
    """A sj fragment is certified RNA in exactly the sense a spliced crossing is withheld for — it
    has no gDNA candidate in the EM — so counting it would load the RNA side of a split that
    arbitrates only unspliced fragments. A locus whose RNA is fully spliced should get a near-zero
    ``rna_prior_count``: its unspliced fragments really are gDNA or nascent.

    The result carries the flux for QC (`test_calibrate`); ``assemble_priors`` must ignore it, and
    that is a deliberate asymmetry rather than an oversight.
    """
    base = _result(
        region_g=[1.0, 1.0], region_r=[2.0, 2.0], region_eff=[100.0, 100.0], boundary_eff=[50.0]
    )
    import dataclasses

    loud = dataclasses.replace(
        base,
        count_rna_sj=np.array([10_000.0]),
        sj_mass_per_crossing=np.ones(1),
        n_sj=1,
    )
    ra = _regions([0, 100], [100, 200])
    ml = [_ml(0, [(0, 0, 200)])]
    quiet_priors = assemble_priors(base, ra, ml)
    loud_priors = assemble_priors(loud, ra, ml)
    np.testing.assert_array_equal(quiet_priors.rna_prior_count, loud_priors.rna_prior_count)
    np.testing.assert_array_equal(quiet_priors.gdna_prior_count, loud_priors.gdna_prior_count)


def test_region_split_between_two_loci():
    # One region [0,100) straddling two adjacent loci ([0,50), [50,100)) → overlap shares 0.5/0.5;
    # each locus gets half of every projected quantity. Single region ⇒ no boundary ⇒ span = region_eff.
    cal = _result(region_g=[5.0], region_r=[10.0], region_eff=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 50)]), _ml(1, [(0, 50, 100)])])
    np.testing.assert_allclose(priors.rna_prior_count, [5.0, 5.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [2.5, 2.5])
    np.testing.assert_allclose(priors.gdna_eff_len, [50.0, 50.0])


def test_intergenic_region_dropped():
    # Region 1 ([200,300)) overlaps no locus → its mass is dropped, not allocated.
    cal = _result(region_g=[1.0, 50.0], region_r=[5.0, 99.0], region_eff=[100.0, 100.0])
    ra = _regions([0, 200], [100, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])
    np.testing.assert_allclose(priors.rna_prior_count, [5.0])  # the intergenic 99 is gone
    np.testing.assert_allclose(priors.gdna_prior_count, [1.0])
    assert priors.gdna_eff_len[0] > 0.0


def test_a_locus_keeps_the_outer_boundary_against_its_INTERGENIC_flank():
    """A locus's far-left outer boundary has an intergenic left flank — a region the projection drops.

    A fragment crossing that boundary overlaps the locus, so it is one of its EM candidates and its
    mass must load the locus's prior. Folding a boundary's mass into one flank region loses this
    boundary into the dropped intergenic flank, and then needs an explicit intergenic re-key to get
    it back; there is nothing to re-key when the boundary is its own object touching region 1.

    Were the outer boundary dropped the locus would see 3 rather than 10 — a large under-count with
    no shape error anywhere to give it away.
    """
    # region 0 intergenic, regions 1-2 exonic ⇒ boundary 0 is the far-LEFT outer boundary, boundary 1 is interior.
    cal = _result(
        region_g=[0.0, 0.0, 0.0],
        region_r=[0.0, 0.0, 0.0],
        region_eff=[100.0, 100.0, 100.0],
        boundary_g=[7.0, 3.0],
        boundary_eff=[50.0, 50.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300], signature=[0, BIT_EXON_POS, BIT_EXON_POS])
    ml = [_ml(0, [(0, 100, 300)])]  # the locus is regions 1-2 only
    priors = assemble_priors(cal, ra, ml)
    np.testing.assert_allclose(priors.gdna_prior_count, [10.0])  # 7 + 3, nothing lost
    assert contended_boundaries(ra, ml, 1).size == 0


# --- Laplace shrinkage toward the (effective) span ------------------------------------------------


def test_evidence_free_region_gives_zero_gdna_prior():
    # Acyclic: no observed gDNA ⇒ zero gDNA pseudocount. With G=0 the Laplace-smoothed IPR is
    # (0+1)²/(1/span) = span exactly, so the eff-len is the uniform effective span (single region ⇒
    # span = region_eff = 100) — never a tiny length.
    cal = _result(region_g=[0.0], region_r=[0.0], region_eff=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])
    np.testing.assert_allclose(priors.rna_prior_count, [0.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [0.0])
    np.testing.assert_allclose(priors.gdna_eff_len, [100.0])  # effective-span fallback


def _six_region_ra():
    return _regions([0, 100, 200, 300, 400, 500], [100, 200, 300, 400, 500, 600])


def _stray_on_a_dead_boundary_cal(stray: float) -> CalibrationResult:
    """7 regions at a gDNA density of 1.0, the reference the result states (ρ_ref = 1.0), 6 boundaries;
    region 0 is intergenic. Boundary 0 has support 0 and carries ``stray`` mass."""
    return _result(
        region_g=[100.0] * 7,
        region_r=[0.0] * 7,
        region_eff=[100.0] * 7,
        boundary_g=[stray, 5.0, 50.0, 50.0, 50.0, 50.0],
        boundary_eff=[0.0, 50.0, 50.0, 50.0, 50.0, 50.0],
        gdna_density_global=1.0,
        gdna_reference_density=1.0,
        gdna_reference_members=1,
        efficiency=[1.0, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0],
    )


def test_stray_mass_on_a_zero_opportunity_boundary_never_reaches_the_eff_len():
    """The length reads the efficiencies and the reach, never a boundary's mass: stray mass on a
    zero-opportunity boundary moves the length not at all, while the prior still counts it — its other
    half is `test_prior_units.test_mass_on_a_zero_opportunity_object_STILL_COUNTS_because_a_count_has_no_divisor`.
    ``mass > 0`` with ``support == 0`` is an ordinary configuration: ``contained_eff_length`` is exactly
    0 wherever an object is shorter than that component's shortest fragment, a fifth of the regions on a
    real chromosome, and the solver can still put mass there because ``f_g`` is an inference."""
    ra = _regions(
        list(range(0, 700, 100)), list(range(100, 800, 100)), signature=[0] + [BIT_EXON_POS] * 6
    )
    ml = [_ml(0, [(0, 100, 700)])]  # the locus is regions 1-6; region 0 is intergenic and dropped
    quiet = assemble_priors(_stray_on_a_dead_boundary_cal(0.0), ra, ml).gdna_eff_len[0]
    for stray in (20.0, 5000.0):
        loud = assemble_priors(_stray_on_a_dead_boundary_cal(stray), ra, ml)
        np.testing.assert_allclose(loud.gdna_eff_len, [quiet], rtol=1e-12)
        # non-vacuity: the locus is genuinely contracted (region 1 at efficiency 0.2)
        assert quiet < 6 * 100.0 + 5 * 50.0
        # ...and the stray mass is NOT silently discarded everywhere — the prior still counts it.
        assert loud.gdna_prior_count[0] == pytest.approx(
            assemble_priors(_stray_on_a_dead_boundary_cal(0.0), ra, ml).gdna_prior_count[0] + stray
        )


def test_empty_multiloci_returns_empty():
    cal = _result(region_g=[1.0], region_r=[1.0], region_eff=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [])
    assert priors.rna_prior_count.shape == (0,)
    assert priors.gdna_prior_count.shape == (0,)
    assert priors.gdna_eff_len.shape == (0,)


def test_region_count_mismatch_raises():
    cal = _result(region_g=[1.0, 2.0], region_r=[1.0, 2.0], region_eff=[10.0, 20.0])
    ra = _regions([0], [100])  # 1 region vs calibration's 2
    with pytest.raises(ValueError, match="regions"):
        assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])


def test_gdna_eff_len_is_the_span_when_every_object_sits_at_the_reference():
    # Priors-side factor-1 WITH a reference: 6 regions at the reference (efficiency 1 everywhere) read
    # the span, 850, exactly — the same as with no reference on the result.
    import dataclasses

    ra = _six_region_ra()
    rho = 2.0
    cal = dataclasses.replace(
        _uniform_field(np.full(6, 100.0), np.full(5, 50.0), rho),
        gdna_reference_density=rho,
        gdna_reference_members=1,
    )
    span = 6 * 100.0 + 5 * 50.0  # 850
    np.testing.assert_allclose(
        assemble_priors(cal, ra, [_ml(0, [(0, 0, 600)])]).gdna_eff_len[0], span, rtol=1e-9
    )


def test_the_length_is_the_counts_objects_at_their_efficiencies():
    """The sum stated: the locus is region 0 (support 100 at efficiency 0.1) and its ONE boundary
    (support 50 at q = 1 and efficiency 0.4); the length is 10 + 20. Neither object's MASS enters —
    the masses are the count's business, and the length is consistent with the count by reading the
    same objects, never by re-deriving a density from them."""
    ra = _six_region_ra()
    mg = np.full(6, 100.0)
    mg[0] = 10.0
    eg = np.full(5, 50.0)
    eg[0] = 250.0
    cal = _result(
        region_g=mg,
        region_r=np.zeros(6),
        region_eff=np.full(6, 100.0),
        boundary_g=eg,
        boundary_eff=np.full(5, 50.0),
        gdna_density_global=1.0,
        gdna_reference_density=1.0,
        gdna_reference_members=1,
        efficiency=[0.1, 1.0, 1.0, 1.0, 1.0, 1.0],
        efficiency_boundary=[0.4, 1.0, 1.0, 1.0, 1.0],
    )
    eff = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])]).gdna_eff_len[0]
    np.testing.assert_allclose(eff, 100.0 * 0.1 + 50.0 * 0.4, rtol=1e-12)
    # PERTURBATION: the mass-based plug-in rule would read min(10, 100) + min(250, 50) = 60
    assert not np.isclose(eff, 60.0)
    # ...and a rule over the locus's bases at the regions' efficiencies alone would read 10
    assert not np.isclose(eff, 10.0)


def test_a_locus_yield_below_one_base_is_not_floored():
    """A 100-base support at efficiency 0.001 beside a 50-base crossing support at 0.001 is a yield of
    0.15 fragments per unit abundance and reads 0.15; at efficiency 0 it reads 0, which the EM takes as
    "cannot emit". The 1 bp floor was a geometric guard that under capture clamped every depleted locus
    shorter than a kilobase to one yield, breaking the E-step's invariance to a common thinning."""
    ra = _six_region_ra()
    for c, expect in ((0.001, 0.15), (0.0, 0.0)):
        cal = _result(
            region_g=np.full(6, 100.0),
            region_r=np.zeros(6),
            region_eff=np.full(6, 100.0),
            boundary_g=np.full(5, 50.0),
            boundary_eff=np.full(5, 50.0),
            gdna_density_global=1.0,
            gdna_reference_density=1.0,
            gdna_reference_members=1,
            efficiency=[c, 1.0, 1.0, 1.0, 1.0, 1.0],
            efficiency_boundary=[c, 1.0, 1.0, 1.0, 1.0],
        )
        eff = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])]).gdna_eff_len[0]
        np.testing.assert_allclose(eff, expect, rtol=1e-12, atol=0.0)


# ── the locus projection underneath: a locus collects REGIONS and BOUNDARIES alike ────────────


def _regions_from_bounds(bounds, signature=None, ref_id=None) -> RegionArrays:
    """Regions tiling one (or more) references from the given region_bound positions."""
    bounds = np.asarray(bounds, dtype=np.int64)
    starts, ends = bounds[:-1], bounds[1:]
    n = starts.shape[0]
    rid = np.zeros(n, dtype=np.int32) if ref_id is None else np.asarray(ref_id, np.int32)
    n_refs = int(rid.max()) + 1 if n else 1
    offsets = np.searchsorted(rid, np.arange(n_refs + 1)).astype(np.int32)
    return RegionArrays(
        ref_id=rid,
        start=starts,
        end=ends,
        signature=(
            np.full(n, BIT_EXON_POS, dtype=np.uint8)
            if signature is None
            else np.asarray(signature, np.uint8)
        ),
        strand_class=np.zeros(n, dtype=np.int8),
        region_size_bp=(ends - starts).astype(np.float64),
        ref_offsets=offsets,
        n_refs=n_refs,
    )


# --- the rule -------------------------------------------------------------------------------------


def test_a_locus_of_k_regions_carries_k_plus_1_boundaries():
    """The rule, stated as a count. 5 regions; the locus is the middle 3, flanked by intergenic.

    Its boundaries are the 4 touching those 3 regions: the two interior ones and the two outer ones.
    A form that kept only interior boundaries would give 2, and a left-keying form with no re-key
    would keep only one outer boundary and give 3.
    """
    ra = _regions_from_bounds([0, 100, 200, 300, 400, 500])
    ml = [_ml(0, [(0, 100, 400)])]  # regions 1,2,3
    e, lid, w = _boundary_locus_shares(ra, ml, 1)
    assert sorted(e.tolist()) == [0, 1, 2, 3], "expected the 4 boundaries touching regions 1-3"
    np.testing.assert_allclose(w, 1.0)
    assert set(lid.tolist()) == {0}


def test_an_boundary_between_two_intergenic_regions_belongs_to_no_locus():
    """The complement, and it is what makes the rule a rule rather than "keep everything".

    A fragment crossing a boundary with no locus region on either side overlaps no transcript, so it is a
    candidate nowhere and must load no prior.
    """
    ra = _regions_from_bounds([0, 100, 200, 300, 400])
    ml = [_ml(0, [(0, 300, 400)])]  # region 3 only
    e, lid, _w = _boundary_locus_shares(ra, ml, 1)
    assert sorted(e.tolist()) == [2], "only the boundary touching region 3 may be kept"
    assert set(lid.tolist()) == {0}


def test_a_region_touching_ONE_locus_gives_it_everything():
    """The share is an allocation, not an overlap fraction — easy to misread, and load-bearing.

    ``_region_locus_shares`` normalises across the loci a region touches, so a region overlapping
    exactly one locus contributes all its mass there however small the overlap. That is what makes
    the projection conserve: a region's mass is never partially discarded, only distributed.
    """
    ra = _regions_from_bounds([0, 100, 200])
    ml = [_ml(0, [(0, 0, 160)])]  # region 1 = [100,200) overlaps by only 60 bp
    r_idx, _lid, r_w = _region_locus_shares(ra, ml, 1)
    got = dict(zip(r_idx.tolist(), r_w.tolist()))
    assert got[0] == pytest.approx(1.0)
    assert got[1] == pytest.approx(1.0), "a 60 % overlap with the ONLY locus still allocates 100 %"


def test_an_boundary_takes_the_MAX_share_of_its_two_flanks():
    """A region genuinely split BETWEEN two loci has fractional shares; its boundaries take the larger.

    ``max``, not a sum and not a mean: if a region is part of a locus then its two boundaries are
    part of that locus, so a boundary inherits the stronger of its two flanks' memberships.

    region 0 = [0,100) lies wholly in locus 0; region 1 = [100,200) is split 50/50 between loci 0 and 1.
    The single boundary between them is therefore ``max(1.0, 0.5) = 1.0`` in locus 0 and
    ``max(0.0, 0.5) = 0.5`` in locus 1.
    """
    ra = _regions_from_bounds([0, 100, 200])
    ml = [_ml(0, [(0, 0, 150)]), _ml(1, [(0, 150, 200)])]
    r_idx, r_lid, r_w = _region_locus_shares(ra, ml, 2)
    got = dict(zip(zip(r_idx.tolist(), r_lid.tolist()), r_w.tolist()))
    assert got[(0, 0)] == pytest.approx(1.0)
    assert got[(1, 0)] == pytest.approx(0.5) and got[(1, 1)] == pytest.approx(0.5)
    e, lid, w = _boundary_locus_shares(ra, ml, 2)
    assert e.tolist() == [0, 0]
    assert dict(zip(lid.tolist(), w.tolist())) == {0: pytest.approx(1.0), 1: pytest.approx(0.5)}


def test_the_region_locus_overlap_is_traversed_ONCE_per_assembly(monkeypatch):
    """`_region_locus_shares` says "computed exactly once" and on the pipeline's path it now is: the
    assembler needs the region projection itself, so it hands the SAME triples to the boundary
    projection, which used to traverse for them again. On the deep library that second traversal is
    1.7 s of every run.

    Counted through a spy rather than asserted about the output, because the output was never wrong —
    what was wrong was paying for it twice. PERTURBATION: dropping ``region_shares`` at the call site
    makes this two.
    """
    import rigel.calibration.priors as P

    calls = {"n": 0}
    real = P._region_locus_shares

    def counting(*a, **k):
        calls["n"] += 1
        return real(*a, **k)

    monkeypatch.setattr(P, "_region_locus_shares", counting)
    cal = _result(
        region_g=[1.0, 2.0, 1.5],
        region_r=[3.0, 4.0, 5.0],
        region_eff=[100.0, 200.0, 150.0],
        boundary_eff=[150.0, 150.0],
    )
    ra = _regions([0, 100, 300], [100, 300, 450])
    P.assemble_priors(cal, ra, [_ml(0, [(0, 0, 450)])])
    assert calls["n"] == 1, f"the region-to-locus overlap was traversed {calls['n']} times"


def test_the_boundary_shares_are_the_same_whether_the_triples_are_handed_in_or_not():
    """The new argument is a hand-down, not a second rule: the same triples in gives the same projection
    out, to the bit, as computing them inside."""
    ra = _regions_from_bounds([0, 100, 200])
    ml = [_ml(0, [(0, 0, 150)]), _ml(1, [(0, 150, 200)])]
    inside = _boundary_locus_shares(ra, ml, 2)
    handed = _boundary_locus_shares(ra, ml, 2, region_shares=_region_locus_shares(ra, ml, 2))
    for a, b in zip(inside, handed):
        assert np.array_equal(a, b)


def test_a_contended_boundary_carries_no_mass():
    """The claim the rule rests on, as a measurement rather than an assumption.

    Two adjacent regions in different multi-loci would give their shared boundary to both, so its
    shares sum to 2. That is unreachable for a boundary that carries mass: any fragment crossing it
    overlaps transcripts in both loci, is a candidate in both, and the union-find would already have
    merged them into one multi-locus.

    So the assertion is not "it cannot happen" but "where it happens, the mass is zero" — and the
    projection reports such boundaries rather than silently renormalising them.
    """
    ra = _regions_from_bounds([0, 100, 200])
    ml = [_ml(0, [(0, 0, 100)]), _ml(1, [(0, 100, 200)])]  # adjacent, no intergenic between
    e, lid, w = _boundary_locus_shares(ra, ml, 2)
    assert e.tolist() == [0, 0], "the contended boundary reaches both loci"
    assert sorted(lid.tolist()) == [0, 1]
    assert float(w.sum()) == pytest.approx(2.0), "shares sum above 1 — the reportable configuration"


def test_two_references_do_not_share_an_boundary():
    """A boundary exists only between two regions of the SAME reference; the axis must not straddle refs."""
    ra = _regions_from_bounds([0, 100, 200, 0, 100, 200][:4], ref_id=[0, 0, 1])
    lo, hi = boundary_region_indices(np.asarray(ra.ref_id))
    assert lo.tolist() == [0], "one boundary, inside reference 0 only"
    assert hi.tolist() == [1]


def test_empty_loci_returns_empty():
    ra = _regions_from_bounds([0, 100, 200])
    e, lid, w = _boundary_locus_shares(ra, [], 0)
    assert e.size == lid.size == w.size == 0


def test_the_region_projection_is_unchanged_by_the_refactor():
    """``_project_regions_to_loci`` must keep its exact behaviour: the region half does not change
    when the boundary half does. It is expressed through ``_region_locus_shares``, so this pins that
    the shared helper introduced no drift — shares normalise across the loci a region touches, and a
    region touching none is dropped.
    """
    from rigel.calibration.priors import _project_regions_to_loci

    ra = _regions_from_bounds([0, 100, 200, 300])
    ml = [_ml(0, [(0, 0, 50)]), _ml(1, [(0, 50, 100)])]  # region 0 split 50/50, regions 1-2 outside
    out = _project_regions_to_loci(ra, ml, 2, {"m": np.array([10.0, 99.0, 99.0])})
    np.testing.assert_allclose(out["m"], [5.0, 5.0])
