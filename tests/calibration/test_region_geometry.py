"""The chain's geometry: its shape, its per-slot arrays, its structural flags and its numbers.

A reference with ``k`` regions lays out ``2k - 1`` chain slots running REGION, BOUNDARY, …, REGION,
so a BOUNDARY always has a region on both sides and no slot is a reference terminal. This file
gates that shape, the ``RegionArrays`` the slots are built from with the two-way mapping between a
region and the contiguous boundary on its right, the splice graph's structural flags landing on the
payload's own boundary axis, and ``RegionGeometry`` itself — one set of numbers per slot, with no
face concept in it. Two strand axes meet here and are asserted apart: ``count`` is keyed by genome
strand, ``sj_count`` by the transcript strand each junction carries. The brute-force oracles share
no helper with the implementation and the flag tests match by genomic coordinate rather than by
index, so neither can be fooled by the arithmetic it is checking.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.region_arrays import (
    RegionArrays,
    boundary_region_indices,
    region_right_boundary,
)
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain
from rigel.calibration.region_geometry import (
    RegionGeometry,
    build_region_geometry,
    region_gdna_geometry,
)
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_POS,
    TS_NEG,
    TS_NONE,
    TS_POS,
)
from rigel.calibration.splice_graph import (
    FLAG_ACCEPTOR_POS,
    FLAG_DONOR_POS,
    FLAG_TES_POS,
    FLAG_TSS_POS,
    SpliceJunctionGeometry,
    build_boundary_flags_array,
    build_region_partition_arrays,
    is_splice_site,
    is_terminus,
)
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.types import Strand

from _index_builder import build_test_index
from _synthetic import make_synthetic_payload


# ---------------------------------------------------------------------------
# brute force — integer start positions, counted in a Python loop. No shared helper.
# ---------------------------------------------------------------------------


def brute_contained(region_len: int, pmf: np.ndarray) -> float:
    """``E_f[#{s : [s, s+w) lies inside [0, region_len)}]`` by enumeration."""
    total = 0.0
    for w, p in enumerate(pmf):
        if p <= 0.0:
            continue
        starts = sum(1 for s in range(0, region_len + 1) if s >= 0 and s + w <= region_len)
        total += p * starts
    return total


def brute_crossing(pmf: np.ndarray, reach_lo: float, reach_hi: float) -> float:
    """``E_f[#{a : 1 <= a <= w-1, a <= R_lo, w-a <= R_hi}]`` by enumeration.

    ``a`` is how many of the molecule's bases lie to the LEFT of the 0-bp boundary. The molecule must have
    at least one base on each side and must fit in what remains of its own template either way.
    """
    total = 0.0
    for w, p in enumerate(pmf):
        if p <= 0.0:
            continue
        starts = sum(1 for a in range(1, w) if a <= reach_lo and (w - a) <= reach_hi)
        total += p * starts
    return total


def spike_pmf(mean: int, max_size: int = 200) -> np.ndarray:
    pmf = np.zeros(max_size + 1, dtype=np.float64)
    pmf[mean] = 1.0
    return pmf


def two_point_pmf(a: int, b: int, wa: float = 0.5, max_size: int = 200) -> np.ndarray:
    pmf = np.zeros(max_size + 1, dtype=np.float64)
    pmf[a] = wa
    pmf[b] = 1.0 - wa
    return pmf


# ---------------------------------------------------------------------------
# the fixture: the 3-region / 2-boundary / 1-sj synthetic payload
# ---------------------------------------------------------------------------

GDNA_PMF = spike_pmf(50)
RNA_PMF = spike_pmf(80)


@pytest.fixture
def parts():
    payload, region_arrays = make_synthetic_payload()
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    # the fixture's one sj: regions are [0,100) [100,200) [200,300), so an intron running from the
    # end of region 0 to the start of region 2 has its DONOR at boundary 0 and its ACCEPTOR at boundary 1.
    sj = SpliceJunctionGeometry(
        src_region=np.array([0], dtype=np.int64),
        dst_region=np.array([2], dtype=np.int64),
        strand=np.array([Strand.POS], dtype=np.int8),
        reach_lo=np.array([1000.0]),
        reach_hi=np.array([1000.0]),
    )
    return payload, region_arrays, substrate, chain, sj


@pytest.fixture
def geometry(parts):
    payload, region_arrays, substrate, chain, sj = parts
    return build_region_geometry(chain, substrate, region_arrays, sj, GDNA_PMF, RNA_PMF)


# ---------------------------------------------------------------------------
# 1. the dissolution itself
# ---------------------------------------------------------------------------


def test_NO_FIELD_NAMES_A_FACE():
    """No field may name a face. A surviving ``_left``/``_right`` would mean the boundary's two
    sides still carry different numbers, which is the model that put a large share of "spliced"
    mass at positions with no annotated splice site."""
    fields = set(RegionGeometry.__dataclass_fields__)
    for name in fields:
        assert not name.endswith("_left"), f"{name} is a face"
        assert not name.endswith("_right"), f"{name} is a face"
    # and the specific names, so a rename cannot smuggle the concept back
    for dead in (
        "mass_left",
        "mass_right",
        "n_unspl_left",
        "n_unspl_right",
        "eff_gdna_left",
        "eff_gdna_right",
        "eff_rna_left",
        "eff_rna_right",
        "eff_spl_left",
        "eff_spl_right",
        "spliced_pos_left",
        "spliced_neg_right",
        "spliced_n_pos_left",
        "spliced_n_neg_right",
    ):
        assert dead not in fields


def test_ONE_SET_OF_NUMBERS_PER_SLOT(geometry, parts):
    _, _, _, chain, _ = parts
    n = chain.n_slots
    assert geometry.n_slots == n
    assert geometry.unspliced_count.shape == (n, 2)
    assert geometry.eff_gdna.shape == (n,)
    assert geometry.eff_rna.shape == (n,)
    assert geometry.sj_count.shape == (n, 2)
    assert geometry.eff_sj.shape == (n, 2)


def test_the_chain_is_N_E_N_E_N_and_the_geometry_is_addressed_by_SLOT(parts, geometry):
    """3 regions and 2 boundaries interleave into 5 slots. A geometry keyed by region id or by
    boundary id instead would have the wrong length, which is the axis mix-up that drops nearly
    every fragment while every golden test stays green."""
    _, _, _, chain, _ = parts
    assert chain.n_slots == 5
    np.testing.assert_array_equal(chain.kind, [REGION, BOUNDARY, REGION, BOUNDARY, REGION])
    assert geometry.unspliced_count.shape[0] == 5


# ---------------------------------------------------------------------------
# 2. the routing — which population lands on which slot kind
# ---------------------------------------------------------------------------


def test_a_REGION_slot_carries_region_contained_and_an_BOUNDARY_slot_carries_boundary_unspliced(
    geometry, parts
):
    """The two populations live on axes that are off by one per reference, and the fixture gives every
    bank a distinct total so a consumer reading the wrong one cannot pass by coincidence."""
    payload, _, _, chain, _ = parts
    region_slots = np.flatnonzero(np.asarray(chain.kind) == REGION)
    boundary_slots = np.flatnonzero(np.asarray(chain.kind) == BOUNDARY)
    np.testing.assert_array_equal(
        geometry.unspliced_count[region_slots], payload.region_contained_count.astype(np.float64)
    )
    np.testing.assert_array_equal(
        geometry.unspliced_count[boundary_slots],
        payload.boundary_unspliced_count.astype(np.float64),
    )


def test_the_count_columns_are_GENOME_STRAND_unpermuted(geometry, parts):
    """POS is column 0. Storing some banks by genome strand and others by sense is how a large
    share of the spliced deposits land in the opposite column from their unspliced neighbours at
    the same boundary."""
    payload, _, _, chain, _ = parts
    first_region = int(np.flatnonzero(np.asarray(chain.kind) == REGION)[0])
    assert geometry.unspliced_count[first_region, 0] == float(payload.region_contained_count[0, 0])
    assert geometry.unspliced_count[first_region, 1] == float(payload.region_contained_count[0, 1])
    assert payload.region_contained_count[0, 0] != payload.region_contained_count[0, 1]


def test_the_count_IS_the_flux_so_mass_and_n_no_longer_diverge(geometry, parts):
    """A fractional mass (the density numerator) and an integer count (the Poisson power) are two
    numbers only if one fragment's mass is split across objects. The accumulator deposits ``+1`` on
    every object the fragment touched, so there is one number and ``Var(log rho) = 1/n`` is honest
    against it. Assert it is exactly integral."""
    np.testing.assert_array_equal(geometry.unspliced_count, np.rint(geometry.unspliced_count))


def test_a_REGION_carries_NO_MATURE_FLUX_because_contained_is_unspliced_by_construction(
    geometry, parts
):
    """Structural, from the specification: ``region_contained`` is credited only when ``not sj_ids``
    (`_accumulator_reference.deposit`). So a region's contained population can never contain a spliced
    fragment, and any mature flux appearing on a REGION slot is a routing bug."""
    _, _, _, chain, _ = parts
    region_slots = np.asarray(chain.kind) == REGION
    assert np.all(geometry.sj_count[region_slots] == 0.0)


# ---------------------------------------------------------------------------
# 3. the divisors — enumerated, not restated
# ---------------------------------------------------------------------------


def test_a_REGION_divisor_is_the_CONTAINED_placements_count(geometry, parts):
    _, region_arrays, _, chain, _ = parts
    region_slots = np.flatnonzero(np.asarray(chain.kind) == REGION)
    for slot, obj in zip(region_slots, np.asarray(chain.obj_idx)[region_slots]):
        length = int(region_arrays.region_size_bp[obj])
        assert geometry.eff_gdna[slot] == pytest.approx(brute_contained(length, GDNA_PMF))
        assert geometry.eff_rna[slot] == pytest.approx(brute_contained(length, RNA_PMF))


def test_a_CONTIGUOUS_BOUNDARY_divisor_is_the_UNBOUNDED_crossing_count__the_A7_RULING(
    geometry, parts
):
    """At a contiguous boundary both components pass ``UNBOUNDED_REACH``, so both divisors collapse
    to ``mu - 1`` exactly. gDNA is unbounded by physics — its template is the chromosome, so
    ``taper_g = 1`` — and RNA is unbounded by ruling, which defers the RNA taper to a change that
    can be A/B'd on its own.

    This test pins that ruling: turning the RNA taper on must break it, which is the point.
    """
    _, _, _, chain, _ = parts
    boundary_slots = np.flatnonzero(np.asarray(chain.kind) == BOUNDARY)
    mu_g = float(np.dot(np.arange(GDNA_PMF.size), GDNA_PMF))
    mu_r = float(np.dot(np.arange(RNA_PMF.size), RNA_PMF))
    for slot in boundary_slots:
        assert geometry.eff_gdna[slot] == pytest.approx(mu_g - 1.0)
        assert geometry.eff_rna[slot] == pytest.approx(mu_r - 1.0)
        # and against the enumerator, so this is not just restating `fl_mean - 1`
        assert geometry.eff_gdna[slot] == pytest.approx(brute_crossing(GDNA_PMF, 1e12, 1e12))


def test_the_SJ_divisor_uses_its_REAL_EXONIC_REACH__the_other_half_of_A7(parts):
    """The other half of the same ruling: a sj boundary is used only by a molecule that spliced
    across it, so what remains either side is exonic and the reach is real. Leaving the sj divisor
    unbounded ships a divisor wrong by up to 4x at a first exon.

    Here the reach binds: 30 bases of exon either side of an 80 bp molecule.
    """
    payload, region_arrays, substrate, chain, _ = parts
    tight = SpliceJunctionGeometry(
        src_region=np.array([0], dtype=np.int64),
        dst_region=np.array([2], dtype=np.int64),
        strand=np.array([Strand.POS], dtype=np.int8),
        reach_lo=np.array([30.0]),
        reach_hi=np.array([30.0]),
    )
    g = build_region_geometry(chain, substrate, region_arrays, tight, GDNA_PMF, RNA_PMF)
    expected = brute_crossing(RNA_PMF, 30.0, 30.0)
    mu_r_minus_1 = float(np.dot(np.arange(RNA_PMF.size), RNA_PMF)) - 1.0
    assert expected < mu_r_minus_1, (
        "the fixture must make the reach BIND, or the test proves nothing"
    )
    boundary_slots = np.flatnonzero(np.asarray(chain.kind) == BOUNDARY)
    live = [s for s in boundary_slots if g.sj_count[s].sum() > 0]
    assert live, "the sj must reach some boundary slot"
    for slot in live:
        assert g.eff_sj[slot, 0] == pytest.approx(expected)


def test_a_reach_of_ZERO_gives_ZERO_opportunity_and_is_not_a_sentinel(parts):
    """`splice_graph`: a reach of 0 is meaningful — no strand-s molecule can occupy that side, so the
    opportunity is genuinely zero. It must NOT be read as 'unset' and replaced by the mean."""
    _, region_arrays, substrate, chain, _ = parts
    dead = SpliceJunctionGeometry(
        src_region=np.array([0], dtype=np.int64),
        dst_region=np.array([2], dtype=np.int64),
        strand=np.array([Strand.POS], dtype=np.int8),
        reach_lo=np.array([0.0]),
        reach_hi=np.array([1000.0]),
    )
    g = build_region_geometry(chain, substrate, region_arrays, dead, GDNA_PMF, RNA_PMF)
    boundary_slots = np.flatnonzero(np.asarray(chain.kind) == BOUNDARY)
    assert np.all(g.eff_sj[boundary_slots, 0] == 0.0)


def test_a_divisor_of_ZERO_is_NOT_FLOORED(parts):
    """`effective_length`'s own contract: an object with no opportunity must return 0, and the
    caller must treat 0 as "no evidence". Flooring every divisor to ``_EPS`` instead produces
    densities of ~1e9 on the substantial minority of fine-partition regions where the contained
    effective length collapses to exactly 0.

    A 100 bp region cannot contain a 150 bp fragment: the opportunity is exactly zero.
    """
    _, region_arrays, substrate, chain, sj = parts
    huge = spike_pmf(150)
    g = build_region_geometry(chain, substrate, region_arrays, sj, huge, huge)
    region_slots = np.flatnonzero(np.asarray(chain.kind) == REGION)
    assert np.all(g.eff_gdna[region_slots] == 0.0)
    assert np.all(g.eff_rna[region_slots] == 0.0)


def test_the_divisors_differ_between_the_two_COMPONENTS(geometry, parts):
    """gDNA and RNA have different length distributions and therefore different opportunity at the same
    object. A single shared divisor would make the 2x2 deconvolution singular."""
    _, _, _, chain, _ = parts
    region_slots = np.flatnonzero(np.asarray(chain.kind) == REGION)
    assert not np.allclose(geometry.eff_gdna[region_slots], geometry.eff_rna[region_slots])


def test_inv_abundance_is_FILLED_from_the_substrate_banks_and_the_values_are_ABSOLUTE(parts):
    """The fill gate. Every other occurrence of ``inv_abundance`` in the suite constructs it
    synthetically, so the one line that fills it from the substrate needs its own gate. The
    fixture's banks are distinct per object (region contained ``[0.24, 0.42, 0.30]``, boundary
    unspliced ``[0.20, 0.20]``), so a fill that reads the wrong bank or the wrong axis cannot pass.

    Absolute values, hard-coded, never flatness or a reconstruction: an accuracy assertion catches
    every planted off-by-one where a flatness assertion misses the pure scale errors.

    A REGION slot carries the CONTAINED bank, whose expectation is ``rho·P(w<=ell)`` — a truncated
    density shape (TRAPS: a-cancellation-is-conditional-on-its-support). A truncation-free "start"
    bank was A/B'd through this fill and refused on the panel
    (ISSUES: the-truncation-free-region-bank).
    """
    payload, region_arrays, substrate, chain, sj = parts
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    r = kind == REGION

    g = build_region_geometry(chain, substrate, region_arrays, sj, GDNA_PMF, RNA_PMF)
    np.testing.assert_allclose(
        g.inv_abundance[r],
        np.array([0.24, 0.42, 0.30])[obj[r]],
        rtol=0,
        atol=0,
    )
    b = kind == BOUNDARY
    np.testing.assert_allclose(
        g.inv_abundance[b],
        np.array([0.20, 0.20])[np.clip(obj[b], 0, 1)],
        rtol=0,
        atol=0,
    )


def test_a_MIXED_pmf_is_not_collapsed_to_its_mean(parts):
    """The divisor is ``E_f[placements]``, not ``placements(E_f[w])`` — the two differ whenever the
    opportunity is non-linear in ``w``, which is every contained frame. Enumerated on a two-point pmf.
    """
    _, region_arrays, substrate, chain, sj = parts
    pmf = two_point_pmf(20, 150)  # mean 85; 150 does not fit in a 100 bp region, 20 does
    g = build_region_geometry(chain, substrate, region_arrays, sj, pmf, pmf)
    region_slots = np.flatnonzero(np.asarray(chain.kind) == REGION)
    slot = int(region_slots[0])
    length = int(region_arrays.region_size_bp[int(np.asarray(chain.obj_idx)[slot])])
    assert g.eff_gdna[slot] == pytest.approx(brute_contained(length, pmf))
    assert g.eff_gdna[slot] != pytest.approx(brute_contained(length, spike_pmf(85)))


# ---------------------------------------------------------------------------
# 4. the sj -> boundary incidence, and the transcript-strand keying
# ---------------------------------------------------------------------------


def test_a_sj_deposits_on_BOTH_the_donor_and_the_acceptor_boundary(geometry, parts):
    """A sj ``(src, dst)`` has its DONOR at the boundary to the right of ``src`` and its ACCEPTOR at
    the boundary to the left of ``dst``. Molecules leave the template at the first and arrive at the
    second, so both boundaries genuinely saw the flux — and the index states both explicitly rather
    than leaving it to be guessed from the exon bits.

    The fixture's sj runs region 0 -> region 2, so it lands on boundary 0 (donor) and boundary 1 (acceptor)
    — i.e. on BOTH boundary slots.
    """
    payload, _, _, chain, _ = parts
    boundary_slots = np.flatnonzero(np.asarray(chain.kind) == BOUNDARY)
    flux = float(payload.sj_count[0].sum())
    assert flux > 0
    for slot in boundary_slots:
        assert geometry.sj_count[slot, 0] == pytest.approx(flux)


def test_the_mature_flux_is_keyed_by_the_SJ_OWN_STRAND_not_the_align_column(parts):
    """The 'sense is derived, never stored' rule, made executable. The accumulator's ``sj_count``
    columns are the genome strand the read aligned to; the sj's strand is a property of the
    ANNOTATION. A ``-`` sj's flux belongs to the ``-`` transcript column however its reads
    aligned — and the fixture's sj has 9 POS-aligned and 4 NEG-aligned reads, so a column-wise
    copy would put 4 fragments in the wrong place.
    """
    payload, region_arrays, substrate, chain, _ = parts
    neg = SpliceJunctionGeometry(
        src_region=np.array([0], dtype=np.int64),
        dst_region=np.array([2], dtype=np.int64),
        strand=np.array([Strand.NEG], dtype=np.int8),
        reach_lo=np.array([1000.0]),
        reach_hi=np.array([1000.0]),
    )
    g = build_region_geometry(chain, substrate, region_arrays, neg, GDNA_PMF, RNA_PMF)
    boundary_slots = np.flatnonzero(np.asarray(chain.kind) == BOUNDARY)
    total = float(payload.sj_count[0].sum())
    assert payload.sj_count[0, 0] != payload.sj_count[0, 1], (
        "the fixture must make the axes distinct"
    )
    for slot in boundary_slots:
        assert g.sj_count[slot, 1] == pytest.approx(total)  # the - transcript column
        assert g.sj_count[slot, 0] == 0.0


def test_several_sj_on_one_boundary_POOL_their_counts_AND_their_divisors(parts):
    """Two sj sharing a donor boundary are two estimates of one rate, so the pooled statement is
    ``sum(count) / sum(E)`` — the ratio of sums, never the mean of ratios.
    ``rho_bg = sum(g)/sum(E)``). Averaging the divisors instead would mis-weight the deeper sj.
    """
    payload, region_arrays, substrate, chain, _ = parts
    two = SpliceJunctionGeometry(
        src_region=np.array([0, 0], dtype=np.int64),
        dst_region=np.array([2, 2], dtype=np.int64),
        strand=np.array([Strand.POS, Strand.POS], dtype=np.int8),
        reach_lo=np.array([1000.0, 30.0]),
        reach_hi=np.array([1000.0, 30.0]),
    )
    # the fixture payload has one sj row; give the second its own
    import dataclasses

    payload2 = dataclasses.replace(
        payload,
        sj_count=np.array([[9, 4], [5, 1]], dtype=np.uint32),
        # ONE column per sj, not the count's (n_sj, 2): `sj_inv_length_sum` is 1-D in the
        # executable specification (`tests/native/_accumulator_reference.py`), which wins.
        sj_inv_length_sum=np.zeros(2, dtype=np.uint64),
        sj_mass=np.zeros(2, dtype=np.uint64),
        ref_sj_offsets=np.array([0, 2], dtype=np.int64),
    )
    sub2 = CalibrationSubstrate.from_payload(payload2, region_arrays)
    g = build_region_geometry(
        chain,
        substrate=sub2,
        region_arrays=region_arrays,
        sj=two,
        gdna_fl_pmf=GDNA_PMF,
        rna_fl_pmf=RNA_PMF,
    )
    slot = int(np.flatnonzero(np.asarray(chain.kind) == BOUNDARY)[0])
    assert g.sj_count[slot, 0] == pytest.approx(9 + 4 + 5 + 1)
    assert g.eff_sj[slot, 0] == pytest.approx(
        brute_crossing(RNA_PMF, 1e12, 1e12) + brute_crossing(RNA_PMF, 30.0, 30.0)
    )


def two_reference_parts(payload):
    """chr1 with 3 regions / 2 boundaries, then chr2 with 2 regions / 1 boundary.

    The second reference is the point. Slot ids run ``N E N E N`` per reference, so within the
    FIRST reference region ``i`` always sits at slot ``2i`` — which means a geometry that assumes that
    layout instead of reading the chain's own adjacency is indistinguishable from a correct one on any
    single-reference fixture. chr2's region 3 sits at slot 5, not slot 6. This helper is what makes
    that assumption observable: the perturbation ``slot_of_region = 2 * arange`` passes every other
    test in this file.
    """
    import dataclasses

    import pandas as pd

    from rigel.calibration.region_arrays import RegionArrays

    p2 = dataclasses.replace(
        payload,
        ref_region_offsets=np.array([0, 3, 5], dtype=np.int64),
        ref_boundary_offsets=np.array([0, 2, 3], dtype=np.int64),
        ref_sj_offsets=np.array([0, 0, 1], dtype=np.int64),
        n_refs=2,
        region_contained_count=np.vstack(
            [payload.region_contained_count, np.array([[1, 1], [2, 2]], np.uint32)]
        ),
        region_contained_inv_opportunity_sum=np.concatenate(
            [payload.region_contained_inv_opportunity_sum, np.zeros(2, np.uint64)]
        ),
        region_start_count=np.concatenate(
            [payload.region_start_count, np.zeros((2, 2), np.uint32)]
        ),
        region_end_count=np.concatenate([payload.region_end_count, np.zeros((2, 2), np.uint32)]),
        region_span_count=np.concatenate([payload.region_span_count, np.zeros((2, 2), np.uint32)]),
        boundary_unspliced_count=np.vstack(
            [payload.boundary_unspliced_count, np.array([[3, 3]], np.uint32)]
        ),
        boundary_unspliced_inv_length_sum=np.concatenate(
            [payload.boundary_unspliced_inv_length_sum, np.zeros(1, np.uint64)]
        ),
        boundary_spliced_count=np.vstack(
            [payload.boundary_spliced_count, np.zeros((1, 2), np.uint32)]
        ),
        region_bounds=np.array([0, 100, 200, 300, 0, 100, 200], dtype=np.int64),
        ref_region_bound_offsets=np.array([0, 4, 7], dtype=np.int64),
    )
    df = pd.DataFrame(
        {
            "region_id": np.arange(5, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * 3 + ["chr2"] * 2, dtype="string"),
            "start": np.array([0, 100, 200, 0, 100], dtype=np.int64),
            "end": np.array([100, 200, 300, 100, 200], dtype=np.int64),
            "length": np.array([100, 100, 100, 100, 100], dtype=np.int64),
            "signature": np.zeros(5, dtype=np.uint8),
        }
    )
    ra2 = RegionArrays.from_frame(df, {"chr1": 0, "chr2": 1})
    sub2 = CalibrationSubstrate.from_payload(p2, ra2)
    chain2 = build_region_chain(p2.ref_region_offsets, p2.ref_boundary_offsets)
    return p2, ra2, sub2, chain2


def _boundary_slot_of(chain, boundary_obj_id: int) -> int:
    boundary_slots = np.flatnonzero(np.asarray(chain.kind) == BOUNDARY)
    obj = np.asarray(chain.obj_idx)[boundary_slots]
    return int(boundary_slots[obj == boundary_obj_id][0])


def test_a_sj_never_lands_on_ANOTHER_REFERENCE(parts):
    """Getting the incidence wrong shifts a sj onto a neighbouring reference's boundaries — invisible
    in aggregate, and exactly the class of defect the substrate's per-reference offset check exists
    for."""
    payload, _, _, _, _ = parts
    p2, ra2, sub2, chain2 = two_reference_parts(payload)
    j = SpliceJunctionGeometry(
        src_region=np.array([0], dtype=np.int64),  # chr1's first region
        dst_region=np.array([2], dtype=np.int64),
        strand=np.array([Strand.POS], dtype=np.int8),
        reach_lo=np.array([1000.0]),
        reach_hi=np.array([1000.0]),
    )
    g = build_region_geometry(chain2, sub2, ra2, j, GDNA_PMF, RNA_PMF)
    assert g.sj_count[_boundary_slot_of(chain2, 2)].sum() == 0.0, "a chr1 sj reached chr2"


def test_a_sj_on_a_LATER_REFERENCE_lands_on_that_references_own_boundary(parts):
    """The layout-assumption killer. Region ``i`` sits at slot ``2i`` only within the FIRST
    reference; chr2's region 3 sits at slot 5. A geometry that computes the slot arithmetically instead
    of reading ``chain.left``/``chain.right`` passes every single-reference test and then puts chr2's
    mature flux on the wrong slot entirely — here on a REGION, which cannot carry mature at all.

    PERTURBATION: ``slot_of_region = 2 * arange(n_regions)`` passes every other test in this file
    and fails only here.
    """
    payload, _, _, _, _ = parts
    p2, ra2, sub2, chain2 = two_reference_parts(payload)
    j = SpliceJunctionGeometry(
        src_region=np.array([3], dtype=np.int64),  # chr2's FIRST region
        dst_region=np.array([4], dtype=np.int64),  # chr2's second region
        strand=np.array([Strand.POS], dtype=np.int8),
        reach_lo=np.array([1000.0]),
        reach_hi=np.array([1000.0]),
    )
    g = build_region_geometry(chain2, sub2, ra2, j, GDNA_PMF, RNA_PMF)
    flux = float(p2.sj_count[0].sum())
    chr2_boundary = _boundary_slot_of(chain2, 2)
    # donor (right of region 3) and acceptor (left of region 4) are the SAME boundary, so it takes the flux twice
    assert g.sj_count[chr2_boundary, 0] == pytest.approx(2.0 * flux)
    assert np.all(g.sj_count[np.asarray(chain2.kind) == REGION] == 0.0)
    for e in (0, 1):  # chr1's two boundaries saw nothing
        assert g.sj_count[_boundary_slot_of(chain2, e)].sum() == 0.0


# ---------------------------------------------------------------------------
# 5. the derived per-region quantities
# ---------------------------------------------------------------------------


def test_region_gdna_geometry_no_longer_SUMS_TWO_FACES(geometry, parts):
    """Returning ``mass_l + mass_r`` over ``E_l + E_r`` at a boundary needs a ``1/2`` on each side to
    cancel a missing ``1/2``. With one set of numbers per boundary there is nothing to sum and no
    halves to cancel, which is what this pins."""
    _, _, _, chain, _ = parts
    mass, eff = region_gdna_geometry(geometry)
    np.testing.assert_array_equal(mass, geometry.unspliced_count.sum(axis=1))
    np.testing.assert_array_equal(eff, geometry.eff_gdna)


def test_spliced_count_and_sj_count_are_DIFFERENT_POPULATIONS(geometry, parts):
    """The naming hazard this schema exists to remove. Both populations are certified RNA and both
    could be called "mature", so one word would cover both and distinguish neither. They are different
    molecules:

    * ``spliced_count`` — crossed this boundary CONTIGUOUSLY, having spliced somewhere else;
    * ``sj_count`` — never crossed it, it JUMPED from here.

    The fixture makes them differ on both axes at once. At boundary 0 the sj's donor sits on the boundary
    but nothing crossed it contiguously (13 vs 0); at boundary 1 both are live and unequal (13 vs 6). A
    consumer that read one for the other would be off by the whole gene's mature output at a donor
    boundary — two orders of magnitude at a real donor.
    """
    payload, _, _, chain, _ = parts
    boundary_slots = np.flatnonzero(np.asarray(chain.kind) == BOUNDARY)
    spliced = geometry.spliced_count[boundary_slots].sum(axis=1)
    sj = geometry.sj_count[boundary_slots].sum(axis=1)
    np.testing.assert_array_equal(spliced, payload.boundary_spliced_count.sum(axis=1))
    assert sj[0] == pytest.approx(float(payload.sj_count[0].sum()))
    # boundary 0: sj flux with NO contiguous spliced crossing — they cannot be the same array
    assert spliced[0] == 0.0 and sj[0] > 0.0
    # boundary 1: both live, and unequal
    assert spliced[1] > 0.0 and sj[1] > 0.0 and spliced[1] != sj[1]


def test_the_word_MATURE_names_no_field(geometry):
    """It fits both spliced populations, so it cannot be the name of either — one word on two
    concepts names neither. The fields carry the accumulator's own three bank names."""
    fields = set(RegionGeometry.__dataclass_fields__)
    assert not any("mature" in f for f in fields), fields
    assert {"unspliced_count", "spliced_count", "sj_count"} <= fields


# ── the route-summed certified rates ───────────────────────────────────────────────────────────


def test_route_rates_are_the_sum_of_per_route_rates(geometry, parts):
    """`route_rate_lo/hi` must equal the independently derived per-route sum Σ flux_J / A_J over
    each face's junctions — the disjoint-routes law at the observation's source. Re-derived here
    from the sj axis by a separate implementation (TRAPS: self-checking-validator)."""
    _payload, _ra, substrate, chain, sj = parts
    from rigel.calibration.effective_length import crossing_eff_length
    from rigel.types import Strand as _S

    n = int(chain.n_slots)
    want_lo = np.zeros((n, 2))
    want_hi = np.zeros((n, 2))
    if sj.n_sj:
        flux = np.asarray(substrate.sj.count, np.float64).sum(axis=1)
        # the same pmf the fixture built the geometry with
        eff = crossing_eff_length(RNA_PMF, sj.reach_lo, sj.reach_hi)
        slot_of_region = np.zeros(int(chain.n_regions_total), np.int64)
        is_region = np.asarray(chain.kind) == REGION
        slot_of_region[np.asarray(chain.obj_idx)[is_region]] = np.flatnonzero(is_region)
        donor = np.asarray(chain.right)[slot_of_region[np.asarray(sj.src_region, np.int64)]]
        acceptor = np.asarray(chain.left)[slot_of_region[np.asarray(sj.dst_region, np.int64)]]
        col = np.where(np.asarray(sj.strand) == np.int8(_S.POS), 0, 1)
        for j in range(int(sj.n_sj)):
            if eff[j] <= 0:
                continue
            want_lo[donor[j], col[j]] += flux[j] / eff[j]
            want_hi[acceptor[j], col[j]] += flux[j] / eff[j]
    np.testing.assert_allclose(np.asarray(geometry.route_rate_lo), want_lo, rtol=1e-12)
    np.testing.assert_allclose(np.asarray(geometry.route_rate_hi), want_hi, rtol=1e-12)


def test_route_rate_dominates_the_pooled_ratio(geometry):
    """The sum of per-route rates is ≥ the pooled
    ratio-of-sums at every face (equality iff every route agrees), so the pooled form's k-route
    under-read cannot survive. Vacuity-guarded: the fixture must expose live faces."""
    rr = np.asarray(geometry.route_rate_lo) + np.asarray(geometry.route_rate_hi)
    jc = np.asarray(geometry.sj_count_lo) + np.asarray(geometry.sj_count_hi)
    ej = np.asarray(geometry.eff_sj_lo) + np.asarray(geometry.eff_sj_hi)
    live = ej > 0
    assert live.any(), "the fixture must carry sj flux"
    pooled = np.where(live, jc / np.where(live, ej, 1.0), 0.0)
    assert np.all(rr[live] >= pooled[live] - 1e-12)


# ── RegionArrays' own geometry, and the region ↔ contiguous-boundary index mapping ─────────────


def _scrambled_region_df():
    # Two refs, rows deliberately out of (ref_id, start) order.
    rows = [
        ("chr2", 150, 400, BIT_EXON_NEG),
        ("chr1", 100, 200, BIT_EXON_POS),
        ("chr2", 0, 150, 0),
        ("chr1", 0, 100, 0),
        ("chr1", 200, 300, BIT_INTRON_POS),
    ]
    return pd.DataFrame(
        {
            "region_id": np.arange(len(rows), dtype=np.int64),
            "ref_name": pd.array([r[0] for r in rows], dtype="string"),
            "start": np.array([r[1] for r in rows], dtype=np.int64),
            "end": np.array([r[2] for r in rows], dtype=np.int64),
            "length": np.array([r[2] - r[1] for r in rows], dtype=np.int64),
            "signature": np.array([r[3] for r in rows], dtype=np.uint8),
        }
    )


def test_from_region_df_csr_ordering():
    df = _scrambled_region_df()
    ra = RegionArrays.from_frame(df, {"chr1": 0, "chr2": 1})

    # Sorted by (ref_id, start).
    np.testing.assert_array_equal(ra.ref_id, [0, 0, 0, 1, 1])
    np.testing.assert_array_equal(ra.start, [0, 100, 200, 0, 150])
    np.testing.assert_array_equal(ra.end, [100, 200, 300, 150, 400])
    np.testing.assert_array_equal(ra.ref_offsets, [0, 3, 5])
    np.testing.assert_array_equal(ra.region_size_bp, [100, 100, 100, 150, 250])
    np.testing.assert_array_equal(ra.strand_class, [TS_NONE, TS_POS, TS_POS, TS_NONE, TS_NEG])
    assert ra.n_regions == 5
    assert ra.n_refs == 2


def test_from_region_df_requires_signature():
    df = _scrambled_region_df().drop(columns=["signature"])
    try:
        RegionArrays.from_frame(df, {"chr1": 0, "chr2": 1})
    except ValueError as exc:
        assert "signature" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("expected ValueError for missing signature column")


# ---------------------------------------------------------------------------
# The region ↔ contiguous-boundary mapping.
# ---------------------------------------------------------------------------
#
# Topology deliberately MULTI-REFERENCE with three different shapes, because the whole class of
# defect here is "the first reference happens to work". ref0 = 3 regions / 2 boundaries, ref1 = 2 regions /
# 1 boundary, ref2 = 1 region / 0 boundaries (a legal reference that owns no boundary at all). E == N − n_refs
# only counts non-empty refs: 6 regions, 3 refs, 3 boundaries.
REF_ID = np.array([0, 0, 0, 1, 1, 2], dtype=np.int32)
REF_REGION_OFFSETS = np.array([0, 3, 5, 6], dtype=np.int64)
REF_BOUNDARY_OFFSETS = np.array([0, 2, 3, 3], dtype=np.int64)


def test_region_right_boundary_is_minus_one_at_every_reference_end():
    # region 2 is chr0's last, region 4 chr1's last, region 5 chr2's only — none owns a boundary to its right.
    np.testing.assert_array_equal(region_right_boundary(REF_ID), [0, 1, -1, 2, -1, -1])


def test_boundary_region_indices_are_the_adjacent_pair():
    lo, hi = boundary_region_indices(REF_ID)
    np.testing.assert_array_equal(lo, [0, 1, 3])
    np.testing.assert_array_equal(hi, [1, 2, 4])
    # A boundary NEVER straddles two references — the invariant a k+1 axis cannot state.
    assert np.all(REF_ID[lo] == REF_ID[hi])


def test_a_single_region_reference_owns_no_boundary():
    lo, _ = boundary_region_indices(np.array([0, 1, 2], dtype=np.int32))
    assert lo.size == 0
    np.testing.assert_array_equal(
        region_right_boundary(np.array([0, 1, 2], dtype=np.int32)), [-1, -1, -1]
    )


def test_the_two_directions_round_trip():
    lo, hi = boundary_region_indices(REF_ID)
    right = region_right_boundary(REF_ID)
    np.testing.assert_array_equal(
        right[lo], np.arange(lo.size)
    )  # boundary → its left region → itself
    has_boundary = right >= 0
    np.testing.assert_array_equal(lo[right[has_boundary]], np.flatnonzero(has_boundary))


def test_boundary_numbering_matches_the_chain_built_from_the_payload_offsets():
    """The gate that matters: re-derive the SAME numbering by a DIFFERENT algorithm.

    ``boundary_region_indices`` counts adjacent same-reference region pairs; ``build_region_chain`` walks the
    payload's two CSR offset arrays and lays out ``N E N E … N`` slot by slot. They must agree, or the
    calibration result's per-boundary arrays are keyed to a different axis than the payload's — the
    exact class of defect that drops nearly every real fragment while every golden test stays
    green. A validator that called the builder's own helper would prove nothing here
    (TRAPS: self-checking-validator).
    """
    chain = build_region_chain(REF_REGION_OFFSETS, REF_BOUNDARY_OFFSETS)
    is_boundary = np.asarray(chain.kind) == BOUNDARY
    obj = np.asarray(chain.obj_idx)
    left_slot = np.asarray(chain.left)[is_boundary]
    right_slot = np.asarray(chain.right)[is_boundary]

    # the chain's own answer: boundary obj_idx e sits between these two region obj_idx values
    chain_edge_id = obj[is_boundary]
    chain_lo = obj[left_slot]
    chain_hi = obj[right_slot]
    order = np.argsort(chain_edge_id)

    lo, hi = boundary_region_indices(REF_ID)
    np.testing.assert_array_equal(chain_edge_id[order], np.arange(lo.size))
    np.testing.assert_array_equal(chain_lo[order], lo)
    np.testing.assert_array_equal(chain_hi[order], hi)


def test_ref_id_must_be_grouped():
    """A scrambled ``ref_id`` cannot produce a valid boundary axis, and must say so rather than
    silently manufacturing boundaries that straddle references."""
    with pytest.raises(ValueError, match="grouped"):
        region_right_boundary(np.array([0, 1, 0], dtype=np.int32))


# ── the chain's shape: a reference alternates REGION/BOUNDARY and owns no terminal slot ────────


def _chain(regions_per_ref):
    """Build from per-reference region counts, deriving the boundary offsets the way the payload does."""
    region_offsets = np.concatenate([[0], np.cumsum(regions_per_ref)]).astype(np.int64)
    boundaries_per_ref = np.maximum(np.asarray(regions_per_ref) - 1, 0)
    boundary_offsets = np.concatenate([[0], np.cumsum(boundaries_per_ref)]).astype(np.int64)
    return build_region_chain(region_offsets, boundary_offsets)


def test_a_reference_alternates_REGION_BOUNDARY_and_ENDS_ON_A_REGION():
    chain = _chain([4])
    assert chain.n_regions_total == 4 and chain.n_boundaries_total == 3
    assert list(chain.kind) == [REGION, BOUNDARY, REGION, BOUNDARY, REGION, BOUNDARY, REGION]
    assert list(chain.obj_idx) == [0, 0, 1, 1, 2, 2, 3]


def test_the_chain_length_is_2k_minus_1_not_2k_plus_1():
    """The shape in one assertion: `k` regions plus `k - 1` boundaries is `2k - 1` slots, not the
    `2k + 1` a `k + 1` boundary axis would give."""
    for k in (1, 2, 5, 50):
        assert _chain([k]).n_slots == 2 * k - 1


def test_a_SINGLE_REGION_reference_has_no_boundaries_and_is_still_a_slot():
    """1 bp and single-region references are legal — human has many thousands of length-1 regions —
    and a reference with one region has zero interior boundaries, not one."""
    chain = _chain([1])
    assert chain.n_slots == 1
    assert list(chain.kind) == [REGION]
    assert chain.left[0] == -1 and chain.right[0] == -1


def test_a_reference_with_NO_regions_contributes_NOTHING():
    chain = _chain([3, 0, 2])
    assert chain.n_slots == (2 * 3 - 1) + 0 + (2 * 2 - 1)
    assert chain.n_regions_total == 5 and chain.n_boundaries_total == 3


def test_references_do_not_BLEED_into_each_other():
    """The last region of one reference must not be adjacent to the first of the next. A chain that
    wrapped would carry density across a chromosome boundary, which is exactly the class of error that
    survives every aggregate check."""
    chain = _chain([3, 2])
    last_of_first = 2 * 3 - 2  # the final REGION slot of reference 0
    assert chain.kind[last_of_first] == REGION
    assert chain.right[last_of_first] == -1
    assert chain.left[last_of_first + 1] == -1


def test_adjacency_is_symmetric_and_alternates_TYPE():
    chain = _chain([4, 3])
    for slot in range(chain.n_slots):
        right = int(chain.right[slot])
        if right >= 0:
            assert int(chain.left[right]) == slot, "left/right must be inverses"
            assert chain.kind[right] != chain.kind[slot], "the chain is bipartite"


def test_an_BOUNDARY_always_has_a_region_on_BOTH_sides():
    """A boundary sits BETWEEN two regions, so it can never be a terminal. This is the invariant a
    ``k + 1`` shape cannot state, because its terminal boundaries have exactly one flank."""
    chain = _chain([5, 2, 1])
    for slot in np.flatnonzero(np.asarray(chain.kind) == BOUNDARY):
        assert chain.left[slot] >= 0 and chain.right[slot] >= 0


def test_boundary_e_sits_between_region_e_and_region_e_plus_one():
    """The endpoints are IMPLICIT — that is the design's word — so this arithmetic is the contract."""
    chain = _chain([4])
    for slot in np.flatnonzero(np.asarray(chain.kind) == BOUNDARY):
        e = int(chain.obj_idx[slot])
        assert int(chain.obj_idx[chain.left[slot]]) == e
        assert int(chain.obj_idx[chain.right[slot]]) == e + 1


def test_obj_idx_is_a_bijection_onto_each_axis():
    """Every region and every boundary appears exactly once; nothing is visited twice or skipped."""
    chain = _chain([4, 1, 3])
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx)
    assert sorted(idx[kind == REGION]) == list(range(chain.n_regions_total))
    assert sorted(idx[kind == BOUNDARY]) == list(range(chain.n_boundaries_total))


def test_an_INCONSISTENT_boundary_count_is_REFUSED_with_the_arithmetic_named():
    """Both offset arrays come from ONE payload, so a mismatch is an accumulator inconsistency rather
    than a stale index — and the error must say so, or the reader rebuilds the index for nothing."""
    with pytest.raises(ValueError, match="k regions has exactly k-1"):
        build_region_chain(np.array([0, 4], np.int64), np.array([0, 4], np.int64))


def test_the_predecessors_TERMINAL_SLOT_SHAPE_is_GONE():
    """The gate is on the shape, not on the names: ``REGION`` and ``BOUNDARY`` are live constants,
    so banning the words would ban the vocabulary rather than the defect.

    A ``B R B R … R B`` layout gives a reference ``k + 1`` boundary slots, the two outermost
    carrying no data and existing only so every region has an object on each side. The shipped
    layout is ``R B R … B R``: ``2k − 1`` slots, starting and ending with a REGION, and a BOUNDARY
    that always has a region on both sides — an invariant the other shape cannot state. That is
    what must not regress, and it is what this pins.
    """
    chain = _chain([4])
    assert chain.n_slots == 2 * 4 - 1, "k regions must give 2k-1 slots, never k + (k+1)"
    assert chain.kind[0] == REGION and chain.kind[-1] == REGION, (
        "the chain starts and ends with a REGION"
    )
    # the invariant a k+1 axis cannot state: every BOUNDARY has a region on BOTH sides. The SHAPE is
    # what is pinned here, not the spelling of the constant.
    for i in range(chain.n_slots):
        if chain.kind[i] == BOUNDARY:
            assert chain.left[i] >= 0 and chain.right[i] >= 0


# ── the splice graph's structural flags, on the payload's own contiguous-boundary axis ─────────


#: t0 splices [400,700); t1 ENDS at 700, where t0's intron ends — so position 700 is an ACCEPTOR
#: *and* a TES, the case the two independent flag bits exist for. t2 on chr2 keeps the
#: per-reference offsets honest.
GTF = """\
chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t701\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t651\t700\t.\t+\t.\tgene_id "g2"; transcript_id "t1";
chr2\ttest\texon\t301\t500\t.\t+\t.\tgene_id "g3"; transcript_id "t2";
chr2\ttest\texon\t801\t1000\t.\t+\t.\tgene_id "g3"; transcript_id "t2";
"""

REFS = {"chr1": 1500, "chr2": 1500}


@pytest.fixture(scope="module")
def index(tmp_path_factory):
    return build_test_index(tmp_path_factory, GTF, name="boundary_flags", refs=REFS)


def _boundary_positions(index) -> np.ndarray:
    """Genomic position of every contiguous boundary, in boundary order.

    A reference contributing ``c`` bounds owns ``c − 1`` regions and ``c − 2`` interior boundaries,
    and boundary ``e`` sits at bound ``e + 1`` — so the interior bounds, per reference, are the
    boundary coordinates.
    """
    positions, region_bound_offsets, _types = build_region_partition_arrays(index)
    out = []
    for f in range(len(region_bound_offsets) - 1):
        lo, hi = int(region_bound_offsets[f]), int(region_bound_offsets[f + 1])
        if hi - lo >= 2:
            out.append(positions[lo + 1 : hi - 1])
    return np.concatenate(out) if out else np.zeros(0, np.int64)


def test_there_is_EXACTLY_ONE_ENTRY_PER_BOUNDARY_and_no_padding(index):
    """A reference with ``k`` regions contributes ``k − 1`` entries, never ``k + 1``: two data-free
    terminals exist only so every region has an object on each side, and a boundary is defined by
    having a region on both sides."""
    flags = build_boundary_flags_array(index)
    assert flags.dtype == np.uint16
    n_regions = len(index.regions_df)
    n_refs_with_regions = index.regions_df["ref_name"].nunique()
    assert flags.shape[0] == n_regions - n_refs_with_regions
    assert flags.shape == _boundary_positions(index).shape


def test_each_entry_carries_the_flags_of_its_own_GENOMIC_POSITION(index):
    """THE mapping assertion, checked position by position rather than index by index — matching by
    COORDINATE cannot be fooled by an off-by-one in the index arithmetic, which is the failure this
    module exists to catch."""
    flags = build_boundary_flags_array(index)
    positions = _boundary_positions(index)
    regions, boundaries = index.regions_df, index.edges_df

    contiguous = boundaries[boundaries["kind"] == 0]
    end_of_src = regions["end"].to_numpy(np.int64)[contiguous["src"].to_numpy(np.int64)]
    want = dict(zip(end_of_src.tolist(), contiguous["flags"].to_numpy(np.uint16).tolist()))
    got = {int(p): int(f) for p, f in zip(positions, flags)}
    assert got == want
    assert any(want.values()), "the fixture must set SOME flag, or this asserts nothing"


def test_the_terminus_and_splice_site_predicates(index):
    """Position 700 is BOTH a TES (t1 ends) and an ACCEPTOR (t0's intron ends) — the case the 4-bit
    signature cannot represent, and the reason the two predicates are independent bits."""
    flags = build_boundary_flags_array(index)
    at = {int(p): int(f) for p, f in zip(_boundary_positions(index), flags)}

    assert at[700] & FLAG_TES_POS and at[700] & FLAG_ACCEPTOR_POS
    assert at[400] & FLAG_DONOR_POS  # t0's intron starts
    assert at[200] & FLAG_TSS_POS  # t0's 5' end

    arr = np.array([at[700], at[400], at[200]], dtype=np.uint16)
    np.testing.assert_array_equal(is_terminus(arr, Strand.POS), [True, False, True])
    np.testing.assert_array_equal(is_splice_site(arr, Strand.POS), [True, True, False])
    # the fixture is + only, so nothing may leak onto the − strand
    assert not is_terminus(arr, Strand.NEG).any()
    assert not is_splice_site(arr, Strand.NEG).any()


def test_flags_align_with_a_REALLY_SCANNED_payload(tmp_path):
    """The end-to-end statement: the array indexes the payload the C++ accumulator actually
    produced, not a re-derivation of the arithmetic that produced it.

    ``RegionStatics`` then places each boundary's bits on the right chain SLOT — asserted through the chain
    rather than assumed, because that is the second place an off-by-one could hide.
    """
    from rigel.calibration.region_chain import BOUNDARY, build_region_chain
    from rigel.calibration.region_geometry import build_region_statics
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.config import BamScanConfig
    from rigel.pipeline import scan_and_buffer
    from rigel.sim import ReadSimConfig, Scenario

    sc = Scenario("boundary_flags", genome_length=5000, seed=11, work_dir=tmp_path / "ef")
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 50}])
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(2500, 2700), (3000, 3200)], "abundance": 30}])
    result = sc.build_oracle(
        n_fragments=200,
        sim_config=ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=11,
        ),
    )
    index = result.index
    _stats, _sm, _buf, payload = scan_and_buffer(
        str(result.bam_path), index, BamScanConfig(sj_strand_tag="auto")
    )

    flags = build_boundary_flags_array(index)
    assert flags.shape[0] == int(payload.ref_boundary_offsets[-1]), (
        "the flags array does not address the payload's contiguous-boundary axis"
    )

    ra = RegionArrays.from_index(index)
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    statics = build_region_statics(chain, ra, flags)
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, np.int64)
    assert not statics.boundary_flags[kind != BOUNDARY].any(), (
        "a REGION slot carries boundary flags"
    )
    np.testing.assert_array_equal(
        statics.boundary_flags[kind == BOUNDARY], flags[idx[kind == BOUNDARY]]
    )
    assert statics.boundary_flags.any(), "the scenario must set SOME flag, or this asserts nothing"
    sc.cleanup()


def test_a_wrong_length_is_refused(index):
    """A silently mis-sized array would shift every flag by one boundary — invisible in aggregate, and
    exactly what the solver must not inherit. A ``k + 1`` boundary axis is the specific wrong length
    most likely to be handed in, so that is the one used here."""
    from rigel.calibration.region_chain import build_region_chain
    from rigel.calibration.region_geometry import build_region_statics
    from rigel.calibration.region_arrays import RegionArrays

    ra = RegionArrays.from_index(index)
    rno = np.asarray(ra.ref_offsets, np.int64)
    reo = np.zeros_like(rno)
    np.cumsum(np.maximum(np.diff(rno) - 1, 0), out=reo[1:])
    chain = build_region_chain(rno, reo)

    class _Sub:
        pass

    old_shape = rno + np.arange(rno.shape[0], dtype=np.int64)  # a k+1 boundary axis
    with pytest.raises(ValueError, match="one per contiguous boundary"):
        build_region_statics(chain, ra, np.zeros(int(old_shape[-1]), dtype=np.uint16))


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE LOCUS — `locus_blocks`: the chain cut at its message terminals, merged up to a block size.
# Pure topology; which slot is a terminal is the caller's predicate.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def _two_reference_chain():
    """Two references: 7 regions (13 slots) then 3 regions (5 slots). Terminals at REGION slots 0, 6,
    12 (the end of reference 0) and 13 (the start of reference 1) — so the fixture has a terminal that
    is also a reference start, a locus at each end, and one interior locus of 5 slots."""
    from rigel.calibration.region_chain import REGION, build_region_chain

    chain = build_region_chain(np.array([0, 7, 10]), np.array([0, 6, 8]))
    assert chain.n_slots == 18
    terminal = np.zeros(chain.n_slots, bool)
    for s in (0, 6, 12, 13):
        assert chain.kind[s] == REGION
        terminal[s] = True
    return chain, terminal


def test_locus_blocks_own_every_slot_once_in_chain_order_and_cut_only_at_terminals_or_starts():
    from rigel.calibration.region_chain import locus_blocks

    chain, terminal = _two_reference_chain()
    left = np.asarray(chain.left)
    for block_slots in (None, 1, 2, 5, 6, 7, 100):
        blocks = locus_blocks(chain, terminal, block_slots)
        assert blocks[0].start == 0 and blocks[-1].stop == chain.n_slots
        for a, b in zip(blocks, blocks[1:]):
            assert a.stop == b.start, f"blocks do not abut at {a} / {b}"
        for b in blocks:
            assert b.start < b.stop, "an empty block"
            assert terminal[b.start] or left[b.start] < 0, f"{b} starts inside a locus"
            # the slot read beyond the owned range is a terminal linked to the block's last slot
            assert b.end in (b.stop, b.stop + 1)
            if b.end == b.stop + 1:
                assert terminal[b.stop] and left[b.stop] >= 0
            else:
                assert b.stop == chain.n_slots or left[b.stop] < 0


def test_locus_blocks_merge_loci_up_to_the_block_size_and_never_split_one():
    from rigel.calibration.region_chain import LocusBlock, locus_blocks

    chain, terminal = _two_reference_chain()
    # the loci: [0,6) [6,12) [12,13) [13,18) — cuts at 0, 6, 12 (terminals) and 13 (terminal + start)
    one_per_locus = locus_blocks(chain, terminal, 1)
    assert [(b.start, b.stop) for b in one_per_locus] == [(0, 6), (6, 12), (12, 13), (13, 18)]
    # a 6-slot locus does not fit a 5-slot block: it becomes a block of its own length
    assert all(b.stop - b.start <= 6 for b in locus_blocks(chain, terminal, 5))
    assert [(b.start, b.stop) for b in locus_blocks(chain, terminal, 5)] == one_per_locus_ranges(
        one_per_locus
    )
    # 7 slots: [0,6) cannot take [6,12); [6,12) cannot take [12,13)? it can (7 slots): [6,13)
    assert [(b.start, b.stop) for b in locus_blocks(chain, terminal, 7)] == [
        (0, 6),
        (6, 13),
        (13, 18),
    ]
    # one block of the whole chain, reading nothing beyond it
    assert locus_blocks(chain, terminal, None) == [LocusBlock(0, 18, 18)]
    # the read-ahead: [0,6) reads slot 6 (a terminal linked to slot 5); [6,12) reads 12; [12,13)
    # reads NOT 13, which is a reference start — no message crosses a reference boundary
    assert [(b.stop, b.end) for b in one_per_locus] == [(6, 7), (12, 13), (13, 13), (18, 18)]


def one_per_locus_ranges(blocks):
    return [(b.start, b.stop) for b in blocks]


def test_locus_blocks_refuse_a_wrong_shaped_predicate_and_a_zero_block_size():
    from rigel.calibration.region_chain import locus_blocks

    chain, terminal = _two_reference_chain()
    with pytest.raises(ValueError, match="one flag per slot"):
        locus_blocks(chain, terminal[:-1])
    with pytest.raises(ValueError, match="block_slots"):
        locus_blocks(chain, terminal, 0)
