"""``assemble_priors``: the acyclic ``CalibrationResult`` as a per-locus EM prior, and the locus projection it stands
on.

A region owns the fragments contained in it; a boundary owns the fragments that cross it. A locus collects both (its
regions by genomic overlap, its boundaries by touching those regions, so ``k`` contiguous regions carry ``k + 1``
boundaries including the two outer ones), and no boundary's mass is folded into a region's total. The second half of the
module tests that projection, the first half what ``assemble_priors`` builds on it. The gDNA length is gDNA's share of
each object (a region's contained support ``S_r``, a boundary's conserved share ``M_e``) at the object's own capture
efficiency. With no landscape every efficiency is 1, so ``gdna_eff_len == Σ S_r + Σ M_e`` exactly and an unenriched
library contracts nothing; the genomic ``region_size_bp`` in place of ``S`` would fabricate a contraction.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.priors import (
    boundary_locus_shares,
    region_locus_shares,
    assemble_priors,
)
from rigel.calibration.region_arrays import RegionArrays, boundary_region_indices
from rigel.calibration.result import CalibrationResult
from rigel.calibration.signature import BIT_EXON_POS
from rigel.locus import Locus, MultiLocus


def _boundary_shares(region_arrays, multi_loci, n_loci):
    """`boundary_locus_shares` on the region shares it is handed in production."""
    return boundary_locus_shares(
        region_arrays, region_locus_shares(region_arrays, multi_loci, n_loci), n_loci
    )


def contended_boundaries(region_arrays, multi_loci, n_loci) -> np.ndarray:
    """``int64[]``: boundaries whose locus shares sum above 1, i.e. reached by two multi-loci at once."""
    e, _lid, w = _boundary_shares(region_arrays, multi_loci, n_loci)
    if e.size == 0:
        return np.zeros(0, np.int64)
    total = np.zeros(int(e.max()) + 1, dtype=np.float64)
    np.add.at(total, e, w)
    return np.flatnonzero(total > 1.0 + 1e-9).astype(np.int64)


def _result(
    *,
    region_g,
    region_r,
    region_eff,
    boundary_g=None,
    boundary_r=None,
    boundary_eff=None,
    boundary_share=None,
    boundary_spliced=None,
    mass_per_crossing=None,
    gdna_density_global=0.01,
    gdna_bulk_density=None,
    efficiency=None,
    efficiency_boundary=None,
) -> CalibrationResult:
    """Build a result on the three axes. One reference with ``n`` regions owns exactly ``n − 1`` boundaries.

    ``boundary_*`` default to zeros, so a caller that cares only about contained mass writes only the region arrays, but
    the boundary axis still has the right length. ``efficiency`` / ``efficiency_boundary`` are the per-object capture
    efficiencies (1 by default, as a library with no landscape carries). ``boundary_share`` is gDNA's conserved share at
    each boundary; it defaults to the crossing support, which it equals where both flanks exceed every fragment length
    (a crossing start then crosses one boundary, ``q`` at its identity 1).
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
        # Geometry, not a split. 1.0 is the identity: a boundary whose flanks both exceed every fragment length, where
        # one crossing is one fragment. A test exercising K-inflation overrides it.
        boundary_mass_per_crossing=(
            np.ones_like(ez)
            if mass_per_crossing is None
            else np.asarray(mass_per_crossing, dtype=np.float64)
        ),
        count_rna_sj=np.zeros(0, dtype=np.float64),
        boundary_spliced_mass_per_crossing=np.ones_like(ez),
        sj_mass_per_crossing=np.ones(0, dtype=np.float64),
        sj_src_region=np.zeros(0, dtype=np.int64),
        sj_dst_region=np.zeros(0, dtype=np.int64),
        gdna_region_eff_len=region_eff_arr,
        gdna_boundary_eff_len=boundary_eff_arr,
        gdna_boundary_conserved_len=(
            boundary_eff_arr
            if boundary_share is None
            else np.asarray(boundary_share, dtype=np.float64)
        ),
        gdna_frac_region=np.zeros_like(ng),
        rna_pos_frac_region=np.zeros_like(ng),
        rna_neg_frac_region=np.zeros_like(ng),
        gdna_frac_boundary=ez.copy(),
        rna_pos_frac_boundary=ez.copy(),
        rna_neg_frac_boundary=ez.copy(),
        gdna_density_global=gdna_density_global,
        gdna_bulk_density=gdna_bulk_density,
        gdna_capture_weight=1.0,
        gdna_capture_saturation=1.0,
        gdna_capture_efficiency_region=(
            np.ones(n) if efficiency is None else np.asarray(efficiency, dtype=np.float64)
        ),
        gdna_capture_efficiency_boundary=(
            np.ones(ne)
            if efficiency_boundary is None
            else np.asarray(efficiency_boundary, dtype=np.float64)
        ),
        rna_sense_frac=0.9,
        n_regions=n,
        n_boundaries=ne,
        n_sj=0,
    )


def _uniform_field(region_eff, boundary_eff, rho) -> CalibrationResult:
    """A uniform gDNA field: every object's mass is ``ρ × its own effective support``, the accumulator's deposition law
    stated directly (``ρ·E_f[(L−w+1)+]`` contained, ``ρ·E_f[w−1]`` crossing). Then ``gdna_eff_len == span == Σ S``
    exactly, and ``G/eff_len`` recovers the true ρ.
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
        span_bp=sum(e - s for _, s, e in blocks),
        loci=loci,
    )


# --- the bedrock invariant: factor = 1 under uniform gDNA -----------------------------------------


def test_factor_one_under_uniform_gdna():
    # A uniform (unenriched) gDNA field is one level, every efficiency is 1, and gdna_eff_len = span = Σ S exactly:
    # region_eff=[120,200,80] (region 1 is short), boundary shares [120,120] (q = 1), ρ=0.02 over 3 same-ref regions ⇒
    # span = 400 + 240 = 640; and the gDNA per-position rate G/eff_len recovers the true ρ.
    region_eff = [120.0, 200.0, 80.0]
    boundary_eff = [120.0, 120.0]
    rho = 0.02
    span = sum(region_eff) + sum(boundary_eff)  # 640
    cal = _uniform_field(region_eff, boundary_eff, rho)
    ra = _regions([0, 120, 320], [120, 320, 400])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 400)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [span], rtol=1e-9)
    # The count is a conserved fragment count read out of the bank: the contained mass ρ·Σregion_eff (8.0) plus the
    # crossing mass ρ·Σboundary_eff (4.8) rescaled by q, here the identity 1.0. It must not equal `ρ · span_bp` = 8.0, a
    # density rule that would drop the 4.8 crossing fragments because a 0-bp boundary has no genomic span to integrate.
    np.testing.assert_allclose(
        priors.gdna_count, [rho * (sum(region_eff) + sum(boundary_eff))], rtol=1e-9
    )
    assert not np.isclose(priors.gdna_count[0], rho * 400.0)  # not the ρ·span_bp density rule


def test_factor_one_holds_for_any_density():
    # The length does not depend on ρ: a 50000× denser uniform library reads the same span.
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
    # The length is the effective support gdna_region_eff_len, not the genomic region_size_bp: genomic sizes 100+100+100
    # = 300, supports 120+200+80 = 400 plus the crossings 300 ⇒ 700. The count is on these objects, so the length is
    # too.
    region_eff = [120.0, 200.0, 80.0]
    boundary_eff = [150.0, 150.0]
    span = 700.0
    cal = _uniform_field(region_eff, boundary_eff, 0.03)
    ra = _regions([0, 100, 200], [100, 200, 300])  # genomic sizes 100,100,100 (Σ=300) ≠ region_eff
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [span], rtol=1e-9)
    assert not np.isclose(priors.gdna_eff_len[0], 300.0 + sum(boundary_eff))


def test_a_boundary_enters_the_length_at_gdnas_conserved_share_never_the_counts_q():
    """The count converts a boundary's mass by ``q``, the accumulator's mass per crossing pooled over gDNA and RNA; the
    length reads gDNA's own conserved share at the boundary, gDNA's starts split as the deposit rule splits them. Three
    40-bp pieces with fragments longer than a piece: each boundary's share is ``½·40 + ½·40 = 40``, so the length is 30
    + 80 = 110, and it stays 110 when ``q`` moves (RNA crossing the same boundaries moves the pooled ``q``, never gDNA's
    geometry) while the count follows ``q``.
    """
    cal = _result(
        region_g=[0.2, 0.2, 0.2],
        region_r=[0.0, 0.0, 0.0],
        region_eff=[10.0, 10.0, 10.0],
        boundary_g=[3.6, 3.6],
        boundary_eff=[180.0, 180.0],
        boundary_share=[40.0, 40.0],
        mass_per_crossing=[0.5, 0.5],
    )
    ra = _regions([0, 40, 80], [40, 80, 120])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 120)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [110.0], rtol=1e-12)
    import dataclasses

    skewed = assemble_priors(
        dataclasses.replace(cal, boundary_mass_per_crossing=np.array([0.8, 0.8])),
        ra,
        [_ml(0, [(0, 0, 120)])],
    )
    np.testing.assert_allclose(skewed.gdna_eff_len, [110.0], rtol=1e-12)
    assert skewed.gdna_count[0] > priors.gdna_count[0] + 1.0
    assert not np.isclose(priors.gdna_eff_len[0], 30.0 + 0.5 * 360.0)


@pytest.mark.parametrize(
    "lengths",
    [(60,), (50, 60, 70)],
    ids=["one-length", "three-lengths"],
)
def test_the_length_counts_each_crossing_start_once_as_the_count_each_fragment(lengths):
    """The two ledgers count consistently, checked by brute-force enumeration against the deposit rule
    (``tests/native/_accumulator_reference.py``: a crossing fragment deposits a count of +1 at every boundary it crosses
    and a mass summing to 1 across them).

    The count converts a boundary's incidence count to fragments by ``q = mass / count``. The length reads each
    boundary's conserved share (`effective_length.conserved_cut_shares`), which splits a crossing start's unit over the
    boundaries its fragment crosses, so a start is counted once however many boundaries it crosses. With three 40-bp
    pieces and fragments longer than a piece, an incidence total would count a spanning start at both its boundaries and
    read the locus's gDNA at half the field's density. Under a uniform field of one fragment per start position, the
    count and the length both equal the number of distinct starts overlapping the locus, so their ratio is 1.
    """
    from rigel.calibration.effective_length import (
        UNBOUNDED_REACH,
        conserved_cut_shares,
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
    length = (ends - starts).astype(np.float64)
    M_e = np.add(
        *conserved_cut_shares(pmf, length[lo], length[hi], UNBOUNDED_REACH, UNBOUNDED_REACH)
    )

    cal = _result(
        region_g=count_r,
        region_r=np.zeros(n_r),
        region_eff=S_r,
        boundary_g=count_e,
        boundary_eff=S_e,
        boundary_share=M_e,
        mass_per_crossing=q,
    )
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 400, 520)])])
    # the count counts each fragment once
    np.testing.assert_allclose(priors.gdna_count, [distinct_starts_in_locus], rtol=1e-12)
    # and so must the length: the density under a uniform field is the field's, 1
    np.testing.assert_allclose(priors.gdna_eff_len, [distinct_starts_in_locus], rtol=1e-12)
    np.testing.assert_allclose(priors.gdna_count / priors.gdna_eff_len, [1.0], rtol=1e-12)
    # an incidence total would count a spanning start at both of its boundaries
    incidences = float(S_r[1:4].sum() + S_e.sum())
    assert incidences > distinct_starts_in_locus * 1.2
    assert not np.isclose(priors.gdna_eff_len[0], incidences)


def test_every_object_has_the_same_density_under_a_uniform_field():
    """Each object has the field's density under a uniform field, the precondition for ``G / gdna_eff_len`` recovering
    ρ. Regions and boundaries are peers, so each axis is checked on its own rather than through a folded total.
    """
    region_eff = np.array([120.0, 200.0, 80.0])
    boundary_eff = np.array([120.0, 120.0])
    rho = 0.02
    cal = _uniform_field(region_eff, boundary_eff, rho)
    np.testing.assert_allclose(cal.count_gdna_region / cal.gdna_region_eff_len, rho, rtol=1e-9)
    np.testing.assert_allclose(cal.count_gdna_boundary / cal.gdna_boundary_eff_len, rho, rtol=1e-9)


# --- mass / projection (independent of the support choice) ----------------------------------------


def test_single_locus_projects_the_gdna_count():
    # The count is a conserved fragment count. No crossing mass here, so it is the contained mass alone, one deposit per
    # contained fragment: Σm = 4.5, where a ρ·span density rule would give 4.5/750 · 450 = 2.7. The RNA mass beside it
    # (Σm = 12.0) is not read.
    cal = _result(
        region_g=[1.0, 2.0, 1.5],
        region_r=[3.0, 4.0, 5.0],
        region_eff=[100.0, 200.0, 150.0],
        boundary_eff=[150.0, 150.0],
    )
    ra = _regions([0, 100, 300], [100, 300, 450])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 450)])])
    np.testing.assert_allclose(priors.gdna_count, [4.5])
    span = (100.0 + 200.0 + 150.0) + 2 * 150.0
    assert 0.0 < priors.gdna_eff_len[0] <= span + 1e-9


def test_gdna_mass_conservation_regions_plus_boundaries():
    # Conservation of mass: total gDNA = Σ region mass + Σ boundary mass, every object counted once.
    #   regions = [2,3,1] (Σ=6); boundaries = [2,3] (Σ=5)  ⇒  total gDNA = 11.
    cal = _result(
        region_g=[2.0, 3.0, 1.0],
        region_r=[0.0, 0.0, 0.0],
        region_eff=[100.0, 100.0, 100.0],
        boundary_g=[2.0, 3.0],
        boundary_eff=[50.0, 50.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    # The locus covers every region, so it collects every region and every boundary, and the count is their total. Here
    # the count equals the raw sum, 11.0, only because this fixture's q is the identity 1.0; the q ≠ 1 test below and
    # `test_prior_units.py` tell the two rules apart.
    np.testing.assert_allclose(priors.gdna_count, [11.0])
    np.testing.assert_allclose(cal.boundary_mass_per_crossing, 1.0)  # the reason it coincides
    np.testing.assert_allclose(
        priors.gdna_count.sum(), cal.count_gdna_region.sum() + cal.count_gdna_boundary.sum()
    )
    assert contended_boundaries(ra, [_ml(0, [(0, 0, 300)])], 1).size == 0  # nothing double-claimed


def test_the_crossing_mass_is_rescaled_by_the_conserved_share():
    """The one test in this module that separates a conserved count from a raw incidence sum: every other fixture leaves
    ``boundary_mass_per_crossing`` at 1.0, where the two coincide. Here ``q`` is ``[0.5, 0.25]`` (a fragment crossing
    boundary 0 deposited on 2 objects on average, boundary 1 on 4)::

        gDNA = 3 (contained, one deposit each) + 4·0.5 + 8·0.25 = 7.0     raw sum would be 15.0

    The contained term is not rescaled: a contained fragment touches exactly one region and is already a count.
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
    np.testing.assert_allclose(priors.gdna_count, [7.0])
    assert not np.isclose(priors.gdna_count[0], 15.0)  # not the raw incidence sum


def test_the_rna_mass_spliced_or_not_never_enters_the_gdna_count():
    # Calibration's RNA count does not reach the EM (its RNA side is its own count of the locus's fragments less this
    # one), so none of the RNA mass here, region RNA [3,4,5] beside boundary RNA [4,4] of which [1,3] is spliced, may
    # load the gDNA count. q is the identity and there is no gDNA crossing mass, so the count is gDNA's contained 4.5.
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
    np.testing.assert_allclose(priors.gdna_count, [4.5])


def test_the_sj_flux_does_not_enter_the_gdna_count():
    """A sj fragment is certified RNA with no gDNA candidate in the EM, so it has no place in calibration's count of the
    locus's gDNA fragments. Its conserved mass is in the count of every fragment, which the gDNA count is a share of.
    """
    base = _result(
        region_g=[1.0, 1.0], region_r=[2.0, 2.0], region_eff=[100.0, 100.0], boundary_eff=[50.0]
    )
    import dataclasses

    loud = dataclasses.replace(
        base,
        count_rna_sj=np.array([10_000.0]),
        sj_mass_per_crossing=np.ones(1),
        sj_src_region=np.array([0], dtype=np.int64),
        sj_dst_region=np.array([1], dtype=np.int64),
        n_sj=1,
    )
    ra = _regions([0, 100], [100, 200])
    ml = [_ml(0, [(0, 0, 200)])]
    quiet_priors = assemble_priors(base, ra, ml)
    loud_priors = assemble_priors(loud, ra, ml)
    np.testing.assert_array_equal(quiet_priors.gdna_count, loud_priors.gdna_count)
    np.testing.assert_allclose(loud_priors.fragment_count - quiet_priors.fragment_count, [10_000.0])


# --- the fragment count: what the gDNA count is a share of ---------------------------------------


def _three_axis_result() -> CalibrationResult:
    """Three regions, two boundaries and one sj joining regions 0 and 2, with every bank occupied and the two
    conversions different, so no term can stand in for another:

    regions      gDNA [10, 0, 20], RNA [5, 1, 5]                                   -> 41
    boundaries   unspliced gDNA [4, 6] and RNA [4, 2] at q = [0.5, 0.8]            -> 4 + 6.4
                 spliced [5, 1] at their own q = [0.25, 0.4]                       -> 1.25 + 0.4
    sj           7 incidences at 0.5 a crossing                                    -> 3.5

    Every fragment: 56.55. gDNA: 30 + 4·0.5 + 6·0.8 = 36.8.
    """
    import dataclasses

    base = _result(
        region_g=[10.0, 0.0, 20.0],
        region_r=[5.0, 1.0, 5.0],
        region_eff=[100.0, 100.0, 100.0],
        boundary_g=[4.0, 6.0],
        boundary_r=[9.0, 3.0],  # spliced-inclusive
        boundary_spliced=[5.0, 1.0],
        boundary_eff=[50.0, 50.0],
        mass_per_crossing=[0.5, 0.8],
    )
    return dataclasses.replace(
        base,
        boundary_spliced_mass_per_crossing=np.array([0.25, 0.4]),
        count_rna_sj=np.array([7.0]),
        sj_mass_per_crossing=np.array([0.5]),
        sj_src_region=np.array([0], dtype=np.int64),
        sj_dst_region=np.array([2], dtype=np.int64),
        n_sj=1,
    )


def test_the_fragment_count_holds_every_fragment_on_the_locus_once():
    """One locus over every object: its fragment count is the library's, gDNA and RNA, contained, crossing unspliced,
    crossing spliced and jumping, each at its own conversion. The gDNA count is the gDNA part of the same sum."""
    cal = _three_axis_result()
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.fragment_count, [56.55], rtol=1e-12)
    np.testing.assert_allclose(priors.gdna_count, [36.8], rtol=1e-12)
    np.testing.assert_allclose(
        priors.fragment_count, [cal.library_gdna_fragments + cal.library_rna_fragments], rtol=1e-12
    )


def test_the_fragment_count_takes_the_gdna_counts_own_shares():
    """Region 1 is split between two loci. Each locus's two counts take the same share of it and of the boundaries
    that touch it, so the pair of fragment counts adds every region's mass once and each locus's gDNA count stays a
    part of its own fragment count."""
    cal = _result(
        region_g=[10.0, 8.0, 20.0],
        region_r=[5.0, 4.0, 5.0],
        region_eff=[100.0, 100.0, 100.0],
        boundary_eff=[50.0, 50.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 150)]), _ml(1, [(0, 150, 300)])])
    np.testing.assert_allclose(priors.gdna_count, [10.0 + 4.0, 4.0 + 20.0])
    np.testing.assert_allclose(priors.fragment_count, [15.0 + 6.0, 6.0 + 25.0])
    np.testing.assert_allclose(priors.fragment_count.sum(), cal.region_conserved_mass.sum())


def test_an_sj_belongs_to_the_locus_of_the_exon_regions_it_joins():
    """Two loci on one reference with an intergenic region between them, and one sj inside the second. Its conserved
    mass goes to that locus alone; an sj whose regions are in no locus goes nowhere."""
    import dataclasses

    base = _result(
        region_g=np.zeros(5),
        region_r=np.ones(5),
        region_eff=np.full(5, 100.0),
        boundary_eff=np.full(4, 50.0),
    )
    ra = _regions([0, 100, 200, 300, 400], [100, 200, 300, 400, 500])
    loci = [_ml(0, [(0, 0, 100)]), _ml(1, [(0, 200, 500)])]
    quiet = assemble_priors(base, ra, loci).fragment_count
    inside = dataclasses.replace(
        base,
        count_rna_sj=np.array([12.0]),
        sj_mass_per_crossing=np.array([0.5]),
        sj_src_region=np.array([2], dtype=np.int64),
        sj_dst_region=np.array([4], dtype=np.int64),
        n_sj=1,
    )
    np.testing.assert_allclose(assemble_priors(inside, ra, loci).fragment_count - quiet, [0.0, 6.0])
    outside = dataclasses.replace(
        inside,
        sj_src_region=np.array([1], dtype=np.int64),
        sj_dst_region=np.array([1], dtype=np.int64),
    )
    np.testing.assert_allclose(
        assemble_priors(outside, ra, loci).fragment_count - quiet, [0.0, 0.0]
    )


def test_the_gdna_count_never_exceeds_the_fragment_count_in_floating_point():
    """Every term of the gDNA count is a part of a term of the fragment count, summed in the same order, so the
    inequality is exact, with no tolerance, on any masses and any loci; and a locus that is all gDNA has the two
    equal to the bit, a share of exactly one."""
    import dataclasses

    rng = np.random.default_rng(20261114)
    for trial in range(200):
        n = int(rng.integers(3, 12))
        bounds = np.concatenate([[0], np.cumsum(rng.integers(20, 400, n))])
        ra = _regions(bounds[:-1], bounds[1:])
        pure_gdna = trial % 4 == 0
        spliced = np.zeros(n - 1) if pure_gdna else rng.uniform(0.0, 50.0, n - 1)
        base = _result(
            region_g=rng.uniform(0.0, 1e4, n),
            region_r=np.zeros(n) if pure_gdna else rng.uniform(0.0, 1e3, n),
            region_eff=np.full(n, 100.0),
            boundary_g=rng.uniform(0.0, 1e3, n - 1),
            boundary_r=spliced + (0.0 if pure_gdna else rng.uniform(0.0, 300.0, n - 1)),
            boundary_spliced=spliced,
            boundary_eff=np.full(n - 1, 50.0),
            mass_per_crossing=rng.uniform(0.05, 1.0, n - 1),
        )
        j = 0 if pure_gdna else int(rng.integers(0, 6))
        src = np.sort(rng.integers(0, n - 1, j))
        cal = dataclasses.replace(
            base,
            boundary_spliced_mass_per_crossing=rng.uniform(0.05, 1.0, n - 1),
            count_rna_sj=rng.uniform(0.0, 500.0, j),
            sj_mass_per_crossing=rng.uniform(0.05, 1.0, j),
            sj_src_region=src.astype(np.int64),
            sj_dst_region=np.minimum(src + 1 + rng.integers(0, 3, j), n - 1).astype(np.int64),
            n_sj=j,
        )
        cut = int(bounds[int(rng.integers(1, n))])
        loci = [_ml(0, [(0, 0, cut)]), _ml(1, [(0, cut, int(bounds[-1]))])]
        priors = assemble_priors(cal, ra, loci)
        assert np.all(priors.gdna_count <= priors.fragment_count), (trial, priors)
        if pure_gdna:
            np.testing.assert_array_equal(priors.gdna_count, priors.fragment_count)


def test_region_split_between_two_loci():
    # One region [0,100) straddling two adjacent loci ([0,50), [50,100)) gives overlap shares 0.5/0.5, so each locus
    # gets half of every projected quantity. One region means no boundary, so span = region_eff.
    cal = _result(region_g=[5.0], region_r=[10.0], region_eff=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 50)]), _ml(1, [(0, 50, 100)])])
    np.testing.assert_allclose(priors.gdna_count, [2.5, 2.5])
    np.testing.assert_allclose(priors.gdna_eff_len, [50.0, 50.0])


def test_intergenic_region_dropped():
    # Region 1 ([200,300)) overlaps no locus, so its mass is dropped, not allocated.
    cal = _result(region_g=[1.0, 50.0], region_r=[5.0, 99.0], region_eff=[100.0, 100.0])
    ra = _regions([0, 200], [100, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])
    np.testing.assert_allclose(priors.gdna_count, [1.0])  # the intergenic 50 is gone
    assert priors.gdna_eff_len[0] > 0.0


def test_a_locus_keeps_the_outer_boundary_against_its_intergenic_flank():
    """A locus's far-left outer boundary has an intergenic left flank, a region the projection drops. A fragment
    crossing that boundary overlaps the locus, so it is one of its EM candidates and its mass must load the locus's
    count; the boundary is its own object touching region 1, so it is kept without any re-keying. Dropping it would give
    3 rather than 10, an under-count with no shape error to show it.
    """
    # region 0 intergenic, regions 1-2 exonic, so boundary 0 is the far-left outer boundary and boundary 1 is interior.
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
    np.testing.assert_allclose(priors.gdna_count, [10.0])  # 7 + 3, nothing lost
    assert contended_boundaries(ra, ml, 1).size == 0


# --- an evidence-free locus -----------------------------------------------------------------------


def test_evidence_free_region_gives_zero_gdna_count():
    # No observed gDNA means a zero gDNA count. The length reads supports and efficiencies, never the mass, so it is
    # still the effective span (one region at efficiency 1 ⇒ span = region_eff = 100).
    cal = _result(region_g=[0.0], region_r=[0.0], region_eff=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])
    np.testing.assert_allclose(priors.gdna_count, [0.0])
    np.testing.assert_allclose(priors.gdna_eff_len, [100.0])  # the span: the length reads no mass


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
        gdna_bulk_density=1.0,
        efficiency=[1.0, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0],
    )


def test_stray_mass_on_a_zero_opportunity_boundary_never_reaches_the_eff_len():
    """The length reads the efficiencies and the reach, never a boundary's mass, so stray mass on a zero-opportunity
    boundary does not move it, while the count still counts it, because a count has no divisor and the mass is fragments
    the accumulator really deposited. Both halves are asserted. ``mass > 0`` with ``support == 0`` is ordinary:
    ``contained_eff_length`` is 0 wherever an object is shorter than that component's shortest fragment, as for a fifth
    of the regions on a real chromosome, and the solver can still put mass there because ``f_g`` is an inference.
    """
    ra = _regions(
        list(range(0, 700, 100)), list(range(100, 800, 100)), signature=[0] + [BIT_EXON_POS] * 6
    )
    ml = [_ml(0, [(0, 100, 700)])]  # the locus is regions 1-6; region 0 is intergenic and dropped
    quiet = assemble_priors(_stray_on_a_dead_boundary_cal(0.0), ra, ml).gdna_eff_len[0]
    for stray in (20.0, 5000.0):
        loud = assemble_priors(_stray_on_a_dead_boundary_cal(stray), ra, ml)
        np.testing.assert_allclose(loud.gdna_eff_len, [quiet], rtol=1e-12)
        # the locus is contracted (region 1 at efficiency 0.2)
        assert quiet < 6 * 100.0 + 5 * 50.0
        # and the stray mass is not discarded: the count still counts it.
        assert loud.gdna_count[0] == pytest.approx(
            assemble_priors(_stray_on_a_dead_boundary_cal(0.0), ra, ml).gdna_count[0] + stray
        )


def test_empty_multiloci_returns_empty():
    cal = _result(region_g=[1.0], region_r=[1.0], region_eff=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [])
    assert priors.gdna_count.shape == (0,)
    assert priors.gdna_eff_len.shape == (0,)


def test_region_count_mismatch_raises():
    cal = _result(region_g=[1.0, 2.0], region_r=[1.0, 2.0], region_eff=[10.0, 20.0])
    ra = _regions([0], [100])  # 1 region vs calibration's 2
    with pytest.raises(ValueError, match="regions"):
        assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])


def test_gdna_eff_len_is_the_plain_support_when_every_object_sits_at_the_bulk():
    # Factor 1 with a landscape: 6 regions at the bulk (efficiency 1 everywhere) read the plain support, 850, exactly,
    # as with no landscape on the result.
    import dataclasses

    ra = _six_region_ra()
    rho = 2.0
    cal = dataclasses.replace(
        _uniform_field(np.full(6, 100.0), np.full(5, 50.0), rho),
        gdna_bulk_density=rho,
    )
    span = 6 * 100.0 + 5 * 50.0  # 850
    np.testing.assert_allclose(
        assemble_priors(cal, ra, [_ml(0, [(0, 0, 600)])]).gdna_eff_len[0], span, rtol=1e-9
    )


def test_the_length_is_the_counts_objects_at_their_efficiencies():
    """The locus is region 0 (support 100 at efficiency 0.1) and its one boundary (support 50 at q = 1 and efficiency
    0.4), so the length is 10 + 20. Neither object's mass enters: the length agrees with the count by reading the same
    objects, not by re-deriving a density from them.
    """
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
        gdna_bulk_density=1.0,
        efficiency=[0.1, 1.0, 1.0, 1.0, 1.0, 1.0],
        efficiency_boundary=[0.4, 1.0, 1.0, 1.0, 1.0],
    )
    eff = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])]).gdna_eff_len[0]
    np.testing.assert_allclose(eff, 100.0 * 0.1 + 50.0 * 0.4, rtol=1e-12)
    # a mass-based plug-in rule would read min(10, 100) + min(250, 50) = 60
    assert not np.isclose(eff, 60.0)
    # and a rule over the locus's bases at the regions' efficiencies alone would read 10
    assert not np.isclose(eff, 10.0)


def test_a_locus_yield_below_one_base_is_not_floored():
    """A 100-base support at efficiency 0.001 beside a 50-base crossing support at 0.001 is a yield of 0.15 fragments
    per unit abundance and reads 0.15; at efficiency 0 it reads 0, which the EM takes as "cannot emit". A 1 bp floor
    would clamp every depleted locus shorter than a kilobase under capture to one yield, breaking the E-step's
    invariance to a common thinning."""
    ra = _six_region_ra()
    for c, expect in ((0.001, 0.15), (0.0, 0.0)):
        cal = _result(
            region_g=np.full(6, 100.0),
            region_r=np.zeros(6),
            region_eff=np.full(6, 100.0),
            boundary_g=np.full(5, 50.0),
            boundary_eff=np.full(5, 50.0),
            gdna_density_global=1.0,
            gdna_bulk_density=1.0,
            efficiency=[c, 1.0, 1.0, 1.0, 1.0, 1.0],
            efficiency_boundary=[c, 1.0, 1.0, 1.0, 1.0],
        )
        eff = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])]).gdna_eff_len[0]
        np.testing.assert_allclose(eff, expect, rtol=1e-12, atol=0.0)


def test_an_enriched_locus_reads_above_its_plain_support_and_nothing_clips_it():
    """In the bulk unit an enriched object contributes its share at its efficiency, far above 1, so a captured locus's
    yield exceeds its plain support, and nothing may clip it there."""
    ra = _six_region_ra()
    cal = _result(
        region_g=np.full(6, 100.0),
        region_r=np.zeros(6),
        region_eff=np.full(6, 100.0),
        boundary_g=np.zeros(5),
        boundary_eff=np.full(5, 50.0),
        gdna_density_global=1.0,
        gdna_bulk_density=1e-3,
        efficiency=np.full(6, 850.0),
        efficiency_boundary=np.full(5, 850.0),
    )
    plain = 6 * 100.0 + 5 * 50.0
    eff = assemble_priors(cal, ra, [_ml(0, [(0, 0, 600)])]).gdna_eff_len[0]
    np.testing.assert_allclose(eff, 850.0 * plain, rtol=1e-12)


def test_a_common_factor_on_every_efficiency_scales_the_length_and_not_the_count():
    """The unit is free: the EM reads every component's length in the efficiencies' unit, so one factor on every
    efficiency scales the gDNA length by exactly that factor, and the count, which has no divisor, not at all."""
    ra = _six_region_ra()
    rng = np.random.default_rng(11)
    kw = dict(
        region_g=rng.uniform(0.0, 50.0, 6),
        region_r=np.zeros(6),
        region_eff=rng.uniform(40.0, 160.0, 6),
        boundary_g=rng.uniform(0.0, 20.0, 5),
        boundary_eff=rng.uniform(20.0, 80.0, 5),
        gdna_density_global=1.0,
        gdna_bulk_density=1e-3,
    )
    c_r, c_b = rng.uniform(0.5, 900.0, 6), rng.uniform(0.5, 900.0, 5)
    loci = [_ml(0, [(0, 0, 300)]), _ml(1, [(0, 300, 600)])]
    one = assemble_priors(_result(**kw, efficiency=c_r, efficiency_boundary=c_b), ra, loci)
    scaled = assemble_priors(
        _result(**kw, efficiency=37.5 * c_r, efficiency_boundary=37.5 * c_b), ra, loci
    )
    np.testing.assert_allclose(scaled.gdna_eff_len, 37.5 * one.gdna_eff_len, rtol=1e-12)
    np.testing.assert_array_equal(scaled.gdna_count, one.gdna_count)


# ── the locus projection underneath: a locus collects regions and boundaries alike ────────────


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
        region_size_bp=(ends - starts).astype(np.float64),
        ref_offsets=offsets,
        n_refs=n_refs,
    )


# --- the rule -------------------------------------------------------------------------------------


def test_a_locus_of_k_regions_carries_k_plus_1_boundaries():
    """5 regions; the locus is the middle 3, flanked by intergenic. Its boundaries are the 4 touching those 3 regions:
    the two interior and the two outer ones. Keeping only interior boundaries would give 2, and keying each boundary to
    its left region would keep one outer boundary and give 3.
    """
    ra = _regions_from_bounds([0, 100, 200, 300, 400, 500])
    ml = [_ml(0, [(0, 100, 400)])]  # regions 1,2,3
    e, lid, w = _boundary_shares(ra, ml, 1)
    assert sorted(e.tolist()) == [0, 1, 2, 3], "expected the 4 boundaries touching regions 1-3"
    np.testing.assert_allclose(w, 1.0)
    assert set(lid.tolist()) == {0}


def test_an_boundary_between_two_intergenic_regions_belongs_to_no_locus():
    """A fragment crossing a boundary with no locus region on either side overlaps no transcript, so it is a candidate
    nowhere and loads no prior."""
    ra = _regions_from_bounds([0, 100, 200, 300, 400])
    ml = [_ml(0, [(0, 300, 400)])]  # region 3 only
    e, lid, _w = _boundary_shares(ra, ml, 1)
    assert sorted(e.tolist()) == [2], "only the boundary touching region 3 may be kept"
    assert set(lid.tolist()) == {0}


def test_a_region_touching_one_locus_gives_it_everything():
    """The share is an allocation, not an overlap fraction. ``region_locus_shares`` normalises across the loci a region
    touches, so a region overlapping one locus gives it all its mass however small the overlap; a region's mass is only
    distributed, never partly discarded.
    """
    ra = _regions_from_bounds([0, 100, 200])
    ml = [_ml(0, [(0, 0, 160)])]  # region 1 = [100,200) overlaps by only 60 bp
    r_idx, _lid, r_w = region_locus_shares(ra, ml, 1)
    got = dict(zip(r_idx.tolist(), r_w.tolist()))
    assert got[0] == pytest.approx(1.0)
    assert got[1] == pytest.approx(1.0), "a 60 % overlap with the ONLY locus still allocates 100 %"


def test_a_boundary_takes_the_max_share_of_its_two_flanks():
    """A region split between two loci has fractional shares, and each of its boundaries takes the larger share of its
    two flanks: if a region is part of a locus, so are its boundaries. Region 0 = [0,100) lies wholly in locus 0; region
    1 = [100,200) is split 50/50 between loci 0 and 1. The boundary between them is ``max(1.0, 0.5) = 1.0`` in locus 0
    and ``max(0.0, 0.5) = 0.5`` in locus 1.
    """
    ra = _regions_from_bounds([0, 100, 200])
    ml = [_ml(0, [(0, 0, 150)]), _ml(1, [(0, 150, 200)])]
    r_idx, r_lid, r_w = region_locus_shares(ra, ml, 2)
    got = dict(zip(zip(r_idx.tolist(), r_lid.tolist()), r_w.tolist()))
    assert got[(0, 0)] == pytest.approx(1.0)
    assert got[(1, 0)] == pytest.approx(0.5) and got[(1, 1)] == pytest.approx(0.5)
    e, lid, w = _boundary_shares(ra, ml, 2)
    assert e.tolist() == [0, 0]
    assert dict(zip(lid.tolist(), w.tolist())) == {0: pytest.approx(1.0), 1: pytest.approx(0.5)}


def test_the_region_locus_overlap_is_traversed_once_per_assembly(monkeypatch):
    """The assembler needs the region projection itself, so it hands the same triples to the boundary projection instead
    of traversing again. Counted through a spy, since a second traversal gives the same output and only costs time.
    """
    import rigel.calibration.priors as P

    calls = {"n": 0}
    real = P.region_locus_shares

    def counting(*a, **k):
        calls["n"] += 1
        return real(*a, **k)

    monkeypatch.setattr(P, "region_locus_shares", counting)
    cal = _result(
        region_g=[1.0, 2.0, 1.5],
        region_r=[3.0, 4.0, 5.0],
        region_eff=[100.0, 200.0, 150.0],
        boundary_eff=[150.0, 150.0],
    )
    ra = _regions([0, 100, 300], [100, 300, 450])
    P.assemble_priors(cal, ra, [_ml(0, [(0, 0, 450)])])
    assert calls["n"] == 1, f"the region-to-locus overlap was traversed {calls['n']} times"


def test_a_contended_boundary_carries_no_mass():
    """Two adjacent regions in different multi-loci would give their shared boundary to both, so its shares sum to 2.
    That cannot happen to a boundary that carries mass: a fragment crossing it overlaps transcripts in both loci, and
    the union-find would already have merged them. So where it happens the mass is zero, and the projection reports such
    boundaries rather than renormalising them.
    """
    ra = _regions_from_bounds([0, 100, 200])
    ml = [_ml(0, [(0, 0, 100)]), _ml(1, [(0, 100, 200)])]  # adjacent, no intergenic between
    e, lid, w = _boundary_shares(ra, ml, 2)
    assert e.tolist() == [0, 0], "the contended boundary reaches both loci"
    assert sorted(lid.tolist()) == [0, 1]
    assert float(w.sum()) == pytest.approx(2.0), "shares sum above 1 — the reportable configuration"


def test_two_references_do_not_share_an_boundary():
    """A boundary exists only between two regions of the same reference; the axis must not straddle refs."""
    ra = _regions_from_bounds([0, 100, 200, 0, 100, 200][:4], ref_id=[0, 0, 1])
    lo, hi = boundary_region_indices(np.asarray(ra.ref_id))
    assert lo.tolist() == [0], "one boundary, inside reference 0 only"
    assert hi.tolist() == [1]


def test_empty_loci_returns_empty():
    ra = _regions_from_bounds([0, 100, 200])
    e, lid, w = _boundary_shares(ra, [], 0)
    assert e.size == lid.size == w.size == 0


def test_a_region_split_between_loci_is_shared_and_one_outside_every_locus_dropped():
    """Shares normalise across the loci a region touches, and a region touching none emits nothing."""
    ra = _regions_from_bounds([0, 100, 200, 300])
    ml = [_ml(0, [(0, 0, 50)]), _ml(1, [(0, 50, 100)])]  # region 0 split 50/50, regions 1-2 outside
    r_idx, lid, w = region_locus_shares(ra, ml, 2)
    assert r_idx.tolist() == [0, 0] and sorted(lid.tolist()) == [0, 1]
    np.testing.assert_allclose(w, [0.5, 0.5])
