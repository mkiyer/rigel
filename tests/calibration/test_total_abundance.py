"""Gates for the MEASURED TOTAL — the composition-free per-slot abundance and its wall mask.

The quantity under test: per slot, the TOTAL fragment density (every component pooled), formed with no
composition model anywhere in it. At a BOUNDARY that is the shipped exact banks plus the certified
spliced arm at its incidence divisor; at a REGION it is the START/END banks over ``ell``, side-selected
by the wall rule — a side is exact iff its template distance clears ``w_max - 1``, distances taken at
the component minimum. Written BEFORE the implementation and verified failing, per the falsification
discipline; the brute-force oracles here share no helper with the implementation.

Conventions borrowed from ``test_region_geometry.py``: absolute expected values (rtol=0), every fixture
bank distinct so a wrong-bank read cannot pass, and refusal gates for degenerate inputs.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.effective_length import UNBOUNDED_REACH, crossing_eff_length, fl_mean
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain
from rigel.calibration.region_geometry import build_region_geometry
from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS
from rigel.calibration.splice_graph import (
    MatureWallDistances,
    mature_wall_distances_kernel,
)
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.total_abundance import (
    RegionWallMask,
    build_region_wall_mask,
    build_total_abundance,
    w_max_from_deposited_lengths,
)
from rigel.types import Strand

from _synthetic import delta_pmf, make_synthetic_payload, make_synthetic_sj

GDNA_PMF = delta_pmf(50, 200)
RNA_PMF = delta_pmf(80, 200)  # mu_r - 1 == 79, the certified-spliced incidence divisor


def regions_from_bounds(bounds, signatures, ref_name="chr1"):
    """A real :class:`RegionArrays` over one reference's ``bounds`` — no SimpleNamespace stand-in,
    because the mask builder reads ``start``/``end``/``ref_offsets`` and a stand-in would let an
    attribute drift pass silently."""
    bounds = np.asarray(bounds, dtype=np.int64)
    n = bounds.size - 1
    frame = pd.DataFrame(
        {
            "region_id": np.arange(n, dtype=np.int64),
            "ref_name": pd.array([ref_name] * n, dtype="string"),
            "start": bounds[:-1],
            "end": bounds[1:],
            "length": bounds[1:] - bounds[:-1],
            "signature": np.asarray(signatures, dtype=np.uint8),
        }
    )
    return RegionArrays.from_frame(frame, {ref_name: 0})


def mature_none(n_regions):
    """The no-mature-template-anywhere distances (`covered` all False)."""
    z = np.zeros((n_regions, 2), dtype=np.float64)
    return MatureWallDistances(
        d_low=z.copy(), d_high=z.copy(), covered=np.zeros((n_regions, 2), dtype=bool)
    )


def hand_mask(n_regions, start_exact, end_exact, w_max=201):
    """A :class:`RegionWallMask` stated by hand, for assembly gates that are not about the mask."""
    se = np.asarray(start_exact, dtype=bool)
    ee = np.asarray(end_exact, dtype=bool)
    big = np.full(n_regions, float(w_max), dtype=np.float64)
    return RegionWallMask(
        n_regions=n_regions,
        d_low=np.where(ee, big, 0.0),
        d_high=np.where(se, big, 0.0),
        start_exact=se,
        end_exact=ee,
        w_max=w_max,
    )


# ---------------------------------------------------------------------------
# w_max — read from the support end of deposited_lengths, never chosen
# ---------------------------------------------------------------------------


def test_w_max_is_READ_from_the_support_end_of_deposited_lengths():
    hist = np.zeros(501, dtype=np.uint32)
    hist[50] = 3
    assert w_max_from_deposited_lengths(hist) == 50
    hist[500] = (
        1  # one fragment at the top bin must move the answer — the support end, not a quantile
    )
    assert w_max_from_deposited_lengths(hist) == 500


def test_w_max_REFUSES_an_empty_histogram_AND_SAYS_WHY():
    """The refusal must NAME the quantity. Dropping the check leaves numpy's own "zero-size array to
    reduction operation" — still a ValueError, so an unmatched ``raises`` would pass on an ungated
    implementation (measured: that perturbation did not fire until this match was added)."""
    with pytest.raises(ValueError, match="deposited_lengths"):
        w_max_from_deposited_lengths(np.zeros(300, dtype=np.uint32))


# ---------------------------------------------------------------------------
# The mature wall distances — SPLICED-template bases past each wall, MAX over isoforms
# ---------------------------------------------------------------------------
#
# One reference, bounds 0/100/200/300/400/500 (five 100 bp regions).
#   T0 (POS): exons [0,200) + [300,400)  -> exonic total 300, intron [200,300)
#   T1 (POS): exon  [0,100)              -> exonic total 100 (the short isoform)
#   T2 (NEG): exon  [300,500)            -> exonic total 200
# Hand enumeration (independent of the kernel):
#   R0: T0 d_low 0,   d_high 200 ; T1 d_low 0, d_high 0    -> MAX collapse: (0, 200)
#   R1: T0 d_low 100, d_high 100                            -> (100, 100)
#   R2: intron — no mature template covers it
#   R3: T0 d_low 200, d_high 0   (POS); T2 d_low 0, d_high 100 (NEG)
#   R4: T2 d_low 100, d_high 0   (NEG only)


@pytest.fixture
def kernel_parts():
    ra = regions_from_bounds([0, 100, 200, 300, 400, 500], [0] * 5)
    t_index = np.array([0, 0, 1, 2], dtype=np.int64)
    ref_id = np.zeros(4, dtype=np.int64)
    ex_start = np.array([0, 300, 0, 300], dtype=np.int64)
    ex_end = np.array([200, 400, 100, 500], dtype=np.int64)
    strand_of = np.array([int(Strand.POS), int(Strand.POS), int(Strand.NEG)], dtype=np.int64)
    return ra, t_index, ref_id, ex_start, ex_end, strand_of


def test_mature_distances_are_SPLICED_template_bases_never_genomic(kernel_parts):
    """R1's high wall sits 200 GENOMIC bases from T0's end but only 100 SPLICED bases — the intron
    does not exist on the mature template. A genomic implementation reads 200 here and fails."""
    ra, t, r, a, b, s = kernel_parts
    m = mature_wall_distances_kernel(t, r, a, b, s, ra)
    assert m.d_high[1, 0] == 100.0
    assert m.d_low[1, 0] == 100.0


def test_mature_distances_MAX_collapse_over_isoforms(kernel_parts):
    """R0's high wall: T0 continues 200 spliced bases, T1 ends flush (0). The collapse is the
    MAXIMUM — the wall binds only if it binds for every covering template — so the answer is 200.
    A MIN collapse reads 0 and fails."""
    ra, t, r, a, b, s = kernel_parts
    m = mature_wall_distances_kernel(t, r, a, b, s, ra)
    assert m.d_high[0, 0] == 200.0
    assert m.d_low[0, 0] == 0.0


def test_a_region_inside_an_intron_is_NOT_mature_covered(kernel_parts):
    ra, t, r, a, b, s = kernel_parts
    m = mature_wall_distances_kernel(t, r, a, b, s, ra)
    assert not m.covered[2, 0] and not m.covered[2, 1]
    assert m.covered[0, 0] and m.covered[3, 0] and m.covered[3, 1] and m.covered[4, 1]


def test_the_strand_columns_are_INDEPENDENT(kernel_parts):
    """R3 is covered by T0 (POS, flush high wall) and T2 (NEG, 100 to go). One array per strand;
    a kernel that pools strands cannot get both columns right."""
    ra, t, r, a, b, s = kernel_parts
    m = mature_wall_distances_kernel(t, r, a, b, s, ra)
    assert m.d_high[3, 0] == 0.0 and m.d_low[3, 0] == 200.0
    assert m.d_high[3, 1] == 100.0 and m.d_low[3, 1] == 0.0
    assert not m.covered[4, 0] and m.covered[4, 1]
    assert m.d_high[4, 1] == 0.0 and m.d_low[4, 1] == 100.0


def test_a_template_ending_flush_at_the_wall_is_covered_at_distance_ZERO(kernel_parts):
    """Flush is a DISTANCE of zero, not an absence: T0's last exon ends exactly at R3's high wall.
    Zero must survive into the mask (where it binds) rather than being read as no-template."""
    ra, t, r, a, b, s = kernel_parts
    m = mature_wall_distances_kernel(t, r, a, b, s, ra)
    assert m.covered[3, 0]
    assert m.d_high[3, 0] == 0.0


def test_the_kernel_REFUSES_an_exon_whose_endpoint_is_not_a_region_bound():
    """Region bounds sit at every exon endpoint on a real index; an exon that violates that would make
    every distance silently wrong, so it is refused, never absorbed."""
    ra = regions_from_bounds([0, 100, 200], [0, 0])
    with pytest.raises(ValueError):
        mature_wall_distances_kernel(
            np.array([0], dtype=np.int64),
            np.zeros(1, dtype=np.int64),
            np.array([0], dtype=np.int64),
            np.array([150], dtype=np.int64),  # not a region bound
            np.array([int(Strand.POS)], dtype=np.int64),
            ra,
        )


# ---------------------------------------------------------------------------
# The wall mask — a side is exact iff its component-minimum distance clears w_max - 1
# ---------------------------------------------------------------------------
#
# One reference, bounds 0/1000/1100/2000/3000 (R0..R3), w_max = 201 so the bar is 200.
# gDNA distances (contig): R0 (0, 2000) · R1 (1000, 1900) · R2 (1100, 1000) · R3 (2000, 0).


def gdna_only_parts(signatures=(0, 0, 0, 0)):
    ra = regions_from_bounds([0, 1000, 1100, 2000, 3000], signatures)
    n_boundaries = 3
    far = np.full((n_boundaries, 2), 5000.0)
    return ra, far.copy(), far.copy()


def test_the_gDNA_contig_wall_binds_at_the_OUTERMOST_regions():
    """With no RNA anywhere the only template is the chromosome: the first region's low wall and the
    last region's high wall are the contig edges (distance 0); every interior side clears 200."""
    ra, lo, hi = gdna_only_parts()
    mask = build_region_wall_mask(ra, mature_none(4), lo, hi, w_max=201)
    np.testing.assert_array_equal(mask.start_exact, [True, True, True, False])
    np.testing.assert_array_equal(mask.end_exact, [False, True, True, True])
    assert mask.d_low[0] == 0.0 and mask.d_high[3] == 0.0
    assert mask.d_low[2] == 1100.0 and mask.d_high[2] == 1000.0


def test_a_side_is_exact_IFF_the_minimum_clears_w_max_minus_1():
    """Bracketed at the bar itself: a distance of exactly ``w_max - 1`` is exact and one below is
    not — an off-by-one in the rule fails one of the two."""
    ra, lo, hi = gdna_only_parts(signatures=(0, BIT_EXON_POS, 0, 0))
    mature = mature_none(4)
    mature.covered[1, 0] = True
    mature.d_high[1, 0] = 200.0  # exactly w_max - 1
    mature.d_low[1, 0] = 5000.0
    mask = build_region_wall_mask(ra, mature, lo, hi, w_max=201)
    assert mask.start_exact[1]
    mask = build_region_wall_mask(ra, mature, lo, hi, w_max=202)  # bar rises to 201
    assert not mask.start_exact[1]


def test_a_short_mature_template_makes_its_own_region_INEXACT():
    """The component minimum takes the mature arm where one covers: a 50-base spliced remainder is
    below the bar, so the START side is not exact and the mask reports the binding distance."""
    short = mature_none(4)
    short.covered[1, 0] = True
    short.d_high[1, 0] = 50.0
    short.d_low[1, 0] = 5000.0
    ra, lo, hi = gdna_only_parts(signatures=(0, BIT_EXON_POS, 0, 0))
    mask = build_region_wall_mask(ra, short, lo, hi, w_max=201)
    assert not mask.start_exact[1]
    assert mask.d_high[1] == 50.0


def test_a_mature_covered_region_MUST_carry_the_matching_exon_bit():
    """The ``free_s ∧ mrna_active_s`` licence is REDUNDANT with ``covered`` on a consistent index —
    an exon covering a region puts that strand's exon bit in its signature — so it is ASSERTED here
    rather than carried as an ungated term. A covered flag at an intergenic region means the
    annotation and the distances disagree, and every wall verdict on that slot would be arbitrary."""
    inconsistent = mature_none(4)
    inconsistent.covered[1, 0] = True
    inconsistent.d_high[1, 0] = 50.0
    ra, lo, hi = gdna_only_parts(signatures=(0, 0, 0, 0))
    with pytest.raises(ValueError, match="exon bit"):
        build_region_wall_mask(ra, inconsistent, lo, hi, w_max=201)


def test_the_nascent_licence_keeps_a_ZERO_reach_from_binding_where_NO_transcript_exists():
    """⭐ A reach of 0 is an ANSWER — no template of that strand crosses there — and 40.6 % / 42.9 %
    of contiguous boundaries have none. At a region the annotation admits no RNA at, that zero must
    NOT enter the minimum: the only template there is the contig, so the slot stays exact. At an
    intron region, where the annotation does admit nascent RNA, the same zero binds. Dropping the
    licence collapses the two."""
    ra, lo, hi = gdna_only_parts(signatures=(0, 0, 0, 0))
    lo[:] = 0.0
    hi[:] = 0.0  # no template crosses any boundary of this reference
    mask = build_region_wall_mask(ra, mature_none(4), lo, hi, w_max=201)
    np.testing.assert_array_equal(mask.start_exact, [True, True, True, False])
    np.testing.assert_array_equal(mask.end_exact, [False, True, True, True])

    ra_intron, lo, hi = gdna_only_parts(signatures=(0, BIT_INTRON_POS, 0, 0))
    lo[:] = 0.0
    hi[:] = 0.0
    mask = build_region_wall_mask(ra_intron, mature_none(4), lo, hi, w_max=201)
    assert not mask.start_exact[1] and not mask.end_exact[1]


def test_the_nascent_arm_reads_the_GENOMIC_reach_where_no_exon_covers():
    """An intron region has nascent RNA and no mature template: its RNA distance is the contiguous
    boundary reach at the wall. Short reach binds; long reach clears."""
    ra, lo, hi = gdna_only_parts(signatures=(0, BIT_INTRON_POS, 0, 0))
    hi[1, 0] = 100.0  # boundary 1 is R1's high wall; POS column
    mask = build_region_wall_mask(ra, mature_none(4), lo, hi, w_max=201)
    assert not mask.start_exact[1]
    assert mask.d_high[1] == 100.0
    hi[1, 0] = 300.0
    mask = build_region_wall_mask(ra, mature_none(4), lo, hi, w_max=201)
    assert mask.start_exact[1]


def test_the_mature_wall_binds_BELOW_the_nascent_reach():
    """At an exon region both forms exist and the minimum is over both: a long genomic reach must not
    rescue a short spliced template. This is the measured mature-differential wall as a gate."""
    ra, lo, hi = gdna_only_parts(signatures=(0, BIT_EXON_POS, 0, 0))
    # nascent reach long on purpose (5000); mature spliced template short
    mature = mature_none(4)
    mature.covered[1, 0] = True
    mature.d_high[1, 0] = 120.0
    mature.d_low[1, 0] = 5000.0
    mask = build_region_wall_mask(ra, mature, lo, hi, w_max=201)
    assert not mask.start_exact[1]
    assert mask.d_high[1] == 120.0


def test_a_DOUBLE_WALLED_region_has_no_exact_side():
    """A template shorter than ``2 * (w_max - 1)`` can leave both sides binding; the mask must say so
    on both, and the assembly below turns that into a not-model-free slot rather than a number."""
    ra = regions_from_bounds([0, 100, 273], [0, 0], ref_name="ERCC-1")
    n_boundaries = 1
    lo = np.zeros((n_boundaries, 2))
    hi = np.zeros((n_boundaries, 2))
    mask = build_region_wall_mask(ra, mature_none(2), lo, hi, w_max=500)
    assert not mask.start_exact.any() and not mask.end_exact.any()


# ---------------------------------------------------------------------------
# The assembly — REGION side selection, BOUNDARY banks, one total per slot
# ---------------------------------------------------------------------------
#
# The synthetic payload's chain alternates REGION/BOUNDARY over 5 slots (three regions, two
# boundaries between them). For the REGION gates the END bank
# is replaced so no region has S == E (the shipped fixture's strand sums tie at 11/12/13 on both banks,
# which would hide a side swap): S sums (11, 12, 13), E sums (13, 11, 12), both totals still 36 so the
# ledger holds.


@pytest.fixture
def assembly_parts():
    payload, region_arrays = make_synthetic_payload()
    payload = dataclasses.replace(
        payload,
        region_end_count=np.array([[9, 4], [6, 5], [7, 5]], dtype=np.uint32),
    )
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    geometry = build_region_geometry(
        chain, substrate, region_arrays, make_synthetic_sj(), GDNA_PMF, RNA_PMF
    )
    return payload, region_arrays, substrate, chain, geometry


def slot_of(chain, kind, obj):
    k = np.asarray(chain.kind)
    o = np.asarray(chain.obj_idx)
    (idx,) = np.nonzero((k == kind) & (o == obj))
    assert idx.size == 1
    return int(idx[0])


def test_a_REGION_total_is_the_EXACT_SIDE_over_length_and_the_average_where_both_clear(
    assembly_parts,
):
    """Absolute values, stated by hand: S sums (11, 12, 13), E sums (13, 11, 12), ell = 100.
    R0 start-only -> 0.11; R1 end-only -> 0.11; R2 both -> (13 + 12) / 200 = 0.125. Distinct S and E
    per region, so using the wrong side (or ignoring the mask) cannot pass."""
    payload, ra, substrate, chain, geometry = assembly_parts
    mask = hand_mask(3, start_exact=[True, False, True], end_exact=[False, True, True])
    t = build_total_abundance(chain, substrate, ra, geometry, mask, RNA_PMF)
    assert t.total[slot_of(chain, REGION, 0)] == pytest.approx(11 / 100, rel=0, abs=1e-15)
    assert t.total[slot_of(chain, REGION, 1)] == pytest.approx(11 / 100, rel=0, abs=1e-15)
    assert t.total[slot_of(chain, REGION, 2)] == pytest.approx(25 / 200, rel=0, abs=1e-15)
    assert t.start_used[slot_of(chain, REGION, 0)] and not t.end_used[slot_of(chain, REGION, 0)]
    assert not t.start_used[slot_of(chain, REGION, 1)] and t.end_used[slot_of(chain, REGION, 1)]
    assert t.start_used[slot_of(chain, REGION, 2)] and t.end_used[slot_of(chain, REGION, 2)]


def test_the_START_and_END_banks_are_summed_over_BOTH_strand_columns(assembly_parts):
    """R0's start bank is [6, 5]: a single-column read gives 0.06 or 0.05, never 0.11."""
    payload, ra, substrate, chain, geometry = assembly_parts
    mask = hand_mask(3, start_exact=[True, True, True], end_exact=[False, False, False])
    t = build_total_abundance(chain, substrate, ra, geometry, mask, RNA_PMF)
    assert t.total[slot_of(chain, REGION, 0)] == pytest.approx(11 / 100, rel=0, abs=1e-15)


def test_a_double_walled_REGION_reads_NaN_and_is_NOT_model_free(assembly_parts):
    payload, ra, substrate, chain, geometry = assembly_parts
    mask = hand_mask(3, start_exact=[False, True, True], end_exact=[False, True, True])
    t = build_total_abundance(chain, substrate, ra, geometry, mask, RNA_PMF)
    s = slot_of(chain, REGION, 0)
    assert np.isnan(t.total[s])
    assert not t.model_free[s]
    assert not t.start_used[s] and not t.end_used[s]
    assert t.model_free[slot_of(chain, REGION, 1)]


def test_a_BOUNDARY_total_is_the_exact_banks_PLUS_certified_spliced_at_its_incidence_divisor(
    assembly_parts,
):
    """Stated by hand off the fixture's distinct banks. Boundary inv_length_sum is 0.20 on both;
    the sj (POS, region 0 -> region 2) leaves at the first boundary and enters at the second, each
    with inv_length_sum 1.3; the certified spliced count is 0 at the first and 6 at the second, and
    its divisor is mu_r - 1 = 79::

        first boundary   0.20 + 1.3 + 0/79 = 1.5
        second boundary  0.20 + 1.3 + 6/79

    Every term has its own value, so a missing or double-counted arm cannot pass."""
    payload, ra, substrate, chain, geometry = assembly_parts
    mask = hand_mask(3, start_exact=[True] * 3, end_exact=[True] * 3)
    t = build_total_abundance(chain, substrate, ra, geometry, mask, RNA_PMF)
    lo_b, hi_b = slot_of(chain, BOUNDARY, 0), slot_of(chain, BOUNDARY, 1)
    assert t.total[lo_b] == pytest.approx(0.20 + 1.3, rel=0, abs=1e-12)
    assert t.total[hi_b] == pytest.approx(0.20 + 1.3 + 6.0 / 79.0, rel=0, abs=1e-12)
    assert t.model_free[lo_b] and t.model_free[hi_b]


def test_the_certified_spliced_arm_uses_the_COUNT_never_the_mass(assembly_parts):
    """The fixture's spliced mass at the second boundary is 3.0 where its count is 6: an implementation reading the
    mass lands at 0.20 + 1.3 + 3/79 and fails the absolute assertion above. Pinned separately so the
    failure names the defect."""
    payload, ra, substrate, chain, geometry = assembly_parts
    mask = hand_mask(3, start_exact=[True] * 3, end_exact=[True] * 3)
    t = build_total_abundance(chain, substrate, ra, geometry, mask, RNA_PMF)
    hi_b = slot_of(chain, BOUNDARY, 1)
    wrong = 0.20 + 1.3 + 3.0 / 79.0
    assert abs(t.total[hi_b] - wrong) > 1e-6


def test_the_spliced_divisor_IS_the_unbounded_crossing_eff_length():
    """The certified-divisor rule: the incidence divisor is ``mu_r - 1`` exactly, which is the
    unbounded crossing effective length — one identity, no second derivation. Bracketed with a
    two-point pmf so a mean computed off the wrong support fails."""
    pmf = np.zeros(201)
    pmf[60] = 0.25
    pmf[100] = 0.75
    assert crossing_eff_length(pmf, UNBOUNDED_REACH, UNBOUNDED_REACH) == pytest.approx(
        fl_mean(pmf) - 1.0, rel=0, abs=1e-9
    )
    assert fl_mean(pmf) - 1.0 == pytest.approx(89.0, rel=0, abs=1e-12)


def test_the_builder_REFUSES_a_payload_whose_ledger_does_not_close(assembly_parts):
    """sum(S) == sum(E) is the accumulator's ledger; a payload that violates it is corrupted input
    and the builder must refuse it rather than average two different populations."""
    payload, ra, substrate, chain, geometry = assembly_parts
    bad = dataclasses.replace(
        payload, region_end_count=np.array([[9, 4], [6, 5], [7, 6]], dtype=np.uint32)
    )
    bad_substrate = CalibrationSubstrate.from_payload(bad, ra)
    mask = hand_mask(3, start_exact=[True] * 3, end_exact=[True] * 3)
    with pytest.raises(ValueError):
        build_total_abundance(chain, bad_substrate, ra, geometry, mask, RNA_PMF)


# ---------------------------------------------------------------------------
# The CONSUMER SWAP — `background_abundance`, and the rule that the PAIR changes, not the estimator
# ---------------------------------------------------------------------------


def _bg_parts():
    """A payload whose START banks and CONTAINED banks disagree, so an estimator reading the wrong
    pair cannot pass by coincidence. Region 2 is the intergenic one (signature 0)."""
    payload, ra = make_synthetic_payload()
    payload = dataclasses.replace(
        payload, region_end_count=np.array([[9, 4], [6, 5], [7, 5]], dtype=np.uint32)
    )
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    return payload, ra, substrate


def test_the_background_fit_takes_the_PAIR_and_keeps_its_own_ESTIMATOR():
    """`fit_intron_background` pools `Gamma(Σcounts + ½, Σexposure)`. Handing it the measured pair must
    move the ANSWER while leaving the estimator's own arithmetic intact — so the swapped fit equals
    the shipped fit computed on the swapped pair, exactly."""
    from rigel.calibration.density_deconv import fit_gdna_background, fit_intron_background

    payload, ra, substrate = _bg_parts()
    eff = np.array([500.0, 500.0, 500.0])
    shipped = fit_intron_background(substrate, ra, eff, include_introns=False)

    counts = np.array([11.0, 12.0, 13.0])  # the START sums
    exposure = np.array([100.0, 100.0, 100.0])
    swapped = fit_intron_background(
        substrate, ra, eff, include_introns=False, counts_exposure=(counts, exposure)
    )
    # region 2 is the only intergenic one, so the pool is exactly {2} either way
    expected = fit_gdna_background(counts[2:3], exposure[2:3])
    assert swapped.log_mu_bg == pytest.approx(expected.log_mu_bg, rel=0, abs=1e-15)
    assert swapped.log_mu_bg != shipped.log_mu_bg


def test_a_ZERO_exposure_drops_a_region_from_the_pool_with_no_second_predicate():
    """A double-walled region comes back with exposure 0, and the pool predicate already excludes a
    zero support — so the not-model-free population excludes itself. If it did NOT, the region would
    contribute counts over no exposure and the pooled rate would be infinite."""
    from rigel.calibration.density_deconv import fit_intron_background

    payload, ra, substrate = _bg_parts()
    eff = np.array([500.0, 500.0, 500.0])
    live = fit_intron_background(
        substrate, ra, eff, counts_exposure=(np.array([1.0, 1.0, 7.0]), np.array([10.0, 10.0, 10.0]))
    )
    dropped = fit_intron_background(
        substrate, ra, eff, counts_exposure=(np.array([1.0, 1.0, 7.0]), np.array([10.0, 10.0, 0.0]))
    )
    assert live.informative and live.n_regions == 1
    assert not dropped.informative and dropped.n_regions == 0


def test_the_config_flag_REFUSES_an_unknown_value():
    from rigel.config import CalibrationConfig

    with pytest.raises(ValueError, match="background_abundance"):
        CalibrationConfig(background_abundance="start_bank")


def test_measured_total_REFUSES_without_the_wall_inputs():
    """⛔ It must refuse, not fall back: a background rate that silently changed estimator because an
    argument was missing is worse than either estimator."""
    from _synthetic import make_gdna_fl_pmf, make_strand_models

    from rigel.calibration import calibrate
    from rigel.config import CalibrationConfig

    payload, ra = make_synthetic_payload()
    with pytest.raises(ValueError, match="mature_walls"):
        calibrate(
            payload=payload,
            region_arrays=ra,
            strand_model=make_strand_models(0.95, 40),
            gdna_fl_pmf=make_gdna_fl_pmf(),
            rna_fl_pmf=make_gdna_fl_pmf(),
            config=CalibrationConfig(background_abundance="measured_total"),
            sj=make_synthetic_sj(),
        )


def test_the_shipped_default_is_BIT_IDENTICAL_and_the_flag_is_NOT_INERT():
    """⭐ Both halves matter. The default must reproduce the tree before the flag existed, and the
    non-default must MOVE something — an arm that cannot move a number has not been tested
    (`TRAPS: an-ablation-that-never-ran`)."""
    from _synthetic import make_gdna_fl_pmf, make_strand_models

    from rigel.calibration import calibrate
    from rigel.calibration.splice_graph import MatureWallDistances
    from rigel.config import CalibrationConfig

    payload, ra = make_synthetic_payload()
    payload = dataclasses.replace(payload, deposited_lengths=np.zeros(201, dtype=np.uint32))
    hist = np.zeros(201, dtype=np.uint32)
    hist[60] = 36
    payload = dataclasses.replace(payload, deposited_lengths=hist)

    def run(source, **kw):
        return calibrate(
            payload=payload,
            region_arrays=ra,
            strand_model=make_strand_models(0.95, 40),
            gdna_fl_pmf=make_gdna_fl_pmf(),
            rna_fl_pmf=make_gdna_fl_pmf(),
            config=CalibrationConfig(background_abundance=source),
            sj=make_synthetic_sj(),
            **kw,
        )

    n_r = 3
    z = np.zeros((n_r, 2))
    walls = MatureWallDistances(
        d_low=z.copy(), d_high=z.copy(), covered=np.zeros((n_r, 2), dtype=bool)
    )
    reach = (np.full((2, 2), 5000.0), np.full((2, 2), 5000.0))
    # ⚠ Asserted on the FITTED BACKGROUND rather than on the deliverable, and the reason is honest:
    # this fixture has three regions and ONE intergenic one, so the pooled floor cannot move a
    # deconvolution no matter which pair it is fitted from. What the flag controls at this level is
    # WHICH PAIR the background was fitted from, and that is exactly what is checked here; whether the
    # deliverable moves is a PANEL question (`total_abundance_audit.py` arm ⓕ prices the estimator, and
    # the panel prices the deliverable). A unit test that asserted the deliverable would either be
    # vacuous or would be pinning the fixture's arithmetic.
    d_base: dict = {}
    d_swap: dict = {}
    base = run("contained", _debug=d_base)
    swapped = run("measured_total", mature_walls=walls, boundary_reach=reach, _debug=d_swap)
    assert np.isfinite(base.gdna_density_global) and np.isfinite(swapped.gdna_density_global)
    # ⭐ Asserted on `intron_background` — the intergenic gDNA background that REACHES ψ (via the intron
    # λ-factor). ⚠ It used to assert on the aggregate `background` field, which was deleted 2026-08-21
    # as a second implementation of this same pool that no caller consumed; pointing the gate at the
    # LIVE consumer is what it should have done from the start.
    bg_base = d_base["calibration_priors"].intron_background
    bg_swap = d_swap["calibration_priors"].intron_background
    assert bg_base is not None and bg_swap is not None, "the intron background was never fitted"
    assert bg_base.log_mu_bg != bg_swap.log_mu_bg, (
        "background_abundance='measured_total' fitted the SAME background as 'contained' — the arm "
        "never ran"
    )
