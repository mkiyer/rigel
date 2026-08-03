"""NodeGeometry — ONE set of numbers per chain slot. Every face concept dissolves.

     (S5.e)
    Divisors: ``rigel.calibration.effective_length`` (S5.c) — the one placements formula

⛔ **What this replaces.** ``build_node_geometry`` was 213 lines producing **18 per-face arrays** — a
``_left``/``_right`` pair for the unspliced mass, its integer flux, three effective lengths and four
spliced channels. The pairs existed because a boundary's two sides lay in differently-sized flanks and
therefore had **different divisors**. A contiguous edge is a 0-bp line: one set of numbers, seen
identically from both sides. So the pairs go, and with them the junction-strand routing, the exon-bit
flank gating and the ``_continues``/``_eff_spl_face`` far-boundary machinery — all of which existed to
*guess* which flank a spliced deposit belonged to, because the old accumulator attributed a splice to
the node's edge rather than to the junction's own coordinate. The v8
index states ``(src, dst, strand)`` explicitly, so the guess is replaced by index arithmetic.

⭐ **TWO STRAND AXES, AND THEY ARE NOT THE SAME AXIS.** ``count`` is keyed by **genome** strand — where
the read aligned, the accumulator's one convention. ``junction_count`` is keyed by **transcript** strand,
derived from each junction's own annotated strand exactly as the design prescribes ("sense derived from
a junction's own strand"). Putting them in one array under one name is the two-conventions-in-one-schema
defect the redesign exists to remove, so they are asserted apart here.

⚠ **The brute-force oracles below share no helper with the implementation**.
trap 1). They enumerate integer start positions in an explicit Python loop.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.node_chain import EDGE, NODE, build_node_chain
from rigel.calibration.node_geometry import (
    NodeGeometry,
    build_node_geometry,
    node_global_geometry,
    node_total_density,
)
from rigel.calibration.splice_graph import JunctionGeometry
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.types import Strand

from _synthetic import make_synthetic_payload


# ---------------------------------------------------------------------------
# brute force — integer start positions, counted in a Python loop. No shared helper.
# ---------------------------------------------------------------------------


def brute_contained(node_len: int, pmf: np.ndarray) -> float:
    """``E_f[#{s : [s, s+w) lies inside [0, node_len)}]`` by enumeration."""
    total = 0.0
    for w, p in enumerate(pmf):
        if p <= 0.0:
            continue
        starts = sum(1 for s in range(0, node_len + 1) if s >= 0 and s + w <= node_len)
        total += p * starts
    return total


def brute_crossing(pmf: np.ndarray, reach_lo: float, reach_hi: float) -> float:
    """``E_f[#{a : 1 <= a <= w-1, a <= R_lo, w-a <= R_hi}]`` by enumeration.

    ``a`` is how many of the molecule's bases lie to the LEFT of the 0-bp line. The molecule must have
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
# the fixture: the 3-node / 2-edge / 1-junction synthetic payload
# ---------------------------------------------------------------------------

GDNA_PMF = spike_pmf(50)
RNA_PMF = spike_pmf(80)


@pytest.fixture
def parts():
    payload, region_arrays = make_synthetic_payload()
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    chain = build_node_chain(payload.ref_node_offsets, payload.ref_edge_offsets)
    # the fixture's one junction: nodes are [0,100) [100,200) [200,300), so an intron running from the
    # end of node 0 to the start of node 2 has its DONOR at line 0 and its ACCEPTOR at line 1.
    junctions = JunctionGeometry(
        src_node=np.array([0], dtype=np.int64),
        dst_node=np.array([2], dtype=np.int64),
        strand=np.array([Strand.POS], dtype=np.int8),
        reach_lo=np.array([1000.0]),
        reach_hi=np.array([1000.0]),
    )
    return payload, region_arrays, substrate, chain, junctions


@pytest.fixture
def geometry(parts):
    payload, region_arrays, substrate, chain, junctions = parts
    return build_node_geometry(chain, substrate, region_arrays, junctions, GDNA_PMF, RNA_PMF)


# ---------------------------------------------------------------------------
# 1. the dissolution itself
# ---------------------------------------------------------------------------


def test_NO_FIELD_NAMES_A_FACE():
    """⛔ The 18 per-face arrays are the thing S5.e deletes. A surviving ``_left``/``_right`` would mean
    the boundary's two sides still carry different numbers — which is the model that produced the 18 %
    of 'spliced' mass sitting at positions with no annotated splice site."""
    fields = set(NodeGeometry.__dataclass_fields__)
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
    assert geometry.junction_count.shape == (n, 2)
    assert geometry.eff_junction.shape == (n, 2)


def test_the_chain_is_N_E_N_E_N_and_the_geometry_is_addressed_by_SLOT(parts, geometry):
    """3 nodes and 2 edges interleave into 5 slots. A geometry keyed by node id or by edge id instead
    would have the wrong length, which is the axis mix-up that once dropped 476,719 of 476,732
    fragments while every golden test passed."""
    _, _, _, chain, _ = parts
    assert chain.n_slots == 5
    np.testing.assert_array_equal(chain.kind, [NODE, EDGE, NODE, EDGE, NODE])
    assert geometry.unspliced_count.shape[0] == 5


# ---------------------------------------------------------------------------
# 2. the routing — which population lands on which slot kind
# ---------------------------------------------------------------------------


def test_a_NODE_slot_carries_node_contained_and_an_EDGE_slot_carries_edge_unspliced(
    geometry, parts
):
    """The two populations live on axes that are off by one per reference, and the fixture gives every
    bank a distinct total so a consumer reading the wrong one cannot pass by coincidence."""
    payload, _, _, chain, _ = parts
    node_slots = np.flatnonzero(np.asarray(chain.kind) == NODE)
    edge_slots = np.flatnonzero(np.asarray(chain.kind) == EDGE)
    np.testing.assert_array_equal(
        geometry.unspliced_count[node_slots], payload.node_contained_count.astype(np.float64)
    )
    np.testing.assert_array_equal(
        geometry.unspliced_count[edge_slots], payload.edge_unspliced_count.astype(np.float64)
    )


def test_the_count_columns_are_GENOME_STRAND_unpermuted(geometry, parts):
    """⭐ POS is column 0. The predecessor stored some banks by genome strand and others by sense, which
    is how 40–44 % of ``node_spanning`` deposits landed in the opposite column."""
    payload, _, _, chain, _ = parts
    first_node = int(np.flatnonzero(np.asarray(chain.kind) == NODE)[0])
    assert geometry.unspliced_count[first_node, 0] == float(payload.node_contained_count[0, 0])
    assert geometry.unspliced_count[first_node, 1] == float(payload.node_contained_count[0, 1])
    assert payload.node_contained_count[0, 0] != payload.node_contained_count[0, 1]


def test_the_count_IS_the_flux_so_mass_and_n_no_longer_diverge(geometry, parts):
    """⭐ The predecessor carried ``mass_*`` (fractional, the density numerator) AND ``n_unspl_*``
    (integer, the Poisson power) because the old accumulator split one fragment's mass across objects.
    The new one deposits ``+1`` on every object the fragment touched, so there is ONE number and
    ``Var(log rho) = 1/n`` is honest against it. Assert it is exactly integral."""
    np.testing.assert_array_equal(geometry.unspliced_count, np.rint(geometry.unspliced_count))


def test_a_NODE_carries_NO_MATURE_FLUX_because_contained_is_unspliced_by_construction(
    geometry, parts
):
    """⛔ Structural, from the specification: ``node_contained`` is credited only when ``not sj_ids``
    (`_accumulator_reference.deposit`). So a node's contained population can never contain a spliced
    fragment, and any mature flux appearing on a NODE slot is a routing bug."""
    _, _, _, chain, _ = parts
    node_slots = np.asarray(chain.kind) == NODE
    assert np.all(geometry.junction_count[node_slots] == 0.0)


# ---------------------------------------------------------------------------
# 3. the divisors — enumerated, not restated
# ---------------------------------------------------------------------------


def test_a_NODE_divisor_is_the_CONTAINED_placements_count(geometry, parts):
    _, region_arrays, _, chain, _ = parts
    node_slots = np.flatnonzero(np.asarray(chain.kind) == NODE)
    for slot, obj in zip(node_slots, np.asarray(chain.obj_idx)[node_slots]):
        length = int(region_arrays.region_size_bp[obj])
        assert geometry.eff_gdna[slot] == pytest.approx(brute_contained(length, GDNA_PMF))
        assert geometry.eff_rna[slot] == pytest.approx(brute_contained(length, RNA_PMF))


def test_a_CONTIGUOUS_EDGE_divisor_is_the_UNBOUNDED_crossing_count__the_A7_RULING(geometry, parts):
    """⭐ **A7, ruled 2026-07-30**: at a contiguous edge BOTH components pass ``UNBOUNDED_REACH``, so
    both divisors collapse to ``mu - 1`` exactly. gDNA is unbounded by physics (its template is the
    chromosome, ``taper_g = 1``); RNA is unbounded **by the ruling**, which keeps S5.e varying exactly
    one thing and defers the taper to S5.g where it can be A/B'd against S5.f's first baseline.

    ⚠ This test PINS the ruling. Turning the RNA taper on later must break it — that is the point.
    """
    _, _, _, chain, _ = parts
    edge_slots = np.flatnonzero(np.asarray(chain.kind) == EDGE)
    mu_g = float(np.dot(np.arange(GDNA_PMF.size), GDNA_PMF))
    mu_r = float(np.dot(np.arange(RNA_PMF.size), RNA_PMF))
    for slot in edge_slots:
        assert geometry.eff_gdna[slot] == pytest.approx(mu_g - 1.0)
        assert geometry.eff_rna[slot] == pytest.approx(mu_r - 1.0)
        # and against the enumerator, so this is not just restating `fl_mean - 1`
        assert geometry.eff_gdna[slot] == pytest.approx(brute_crossing(GDNA_PMF, 1e12, 1e12))


def test_the_JUNCTION_divisor_uses_its_REAL_EXONIC_REACH__the_other_half_of_A7(parts):
    """⭐ **A7, ruled**: a junction edge is used only by a molecule that spliced across it, so what
    remains either side is exonic and the reach is real. A junction is a BRAND-NEW population — the
    predecessor had no junction divisor at all — so wiring it regresses nothing, while leaving it
    unbounded would ship a divisor wrong by up to 4x at a first exon (199.0 at
    R=550 against 50.0 at R=50).

    Here the reach BINDS: 30 bases of exon either side of an 80 bp molecule.
    """
    payload, region_arrays, substrate, chain, _ = parts
    tight = JunctionGeometry(
        src_node=np.array([0], dtype=np.int64),
        dst_node=np.array([2], dtype=np.int64),
        strand=np.array([Strand.POS], dtype=np.int8),
        reach_lo=np.array([30.0]),
        reach_hi=np.array([30.0]),
    )
    g = build_node_geometry(chain, substrate, region_arrays, tight, GDNA_PMF, RNA_PMF)
    expected = brute_crossing(RNA_PMF, 30.0, 30.0)
    mu_r_minus_1 = float(np.dot(np.arange(RNA_PMF.size), RNA_PMF)) - 1.0
    assert expected < mu_r_minus_1, (
        "the fixture must make the reach BIND, or the test proves nothing"
    )
    edge_slots = np.flatnonzero(np.asarray(chain.kind) == EDGE)
    live = [s for s in edge_slots if g.junction_count[s].sum() > 0]
    assert live, "the junction must reach some edge slot"
    for slot in live:
        assert g.eff_junction[slot, 0] == pytest.approx(expected)


def test_a_reach_of_ZERO_gives_ZERO_opportunity_and_is_not_a_sentinel(parts):
    """`splice_graph`: a reach of 0 is meaningful — no strand-s molecule can occupy that side, so the
    opportunity is genuinely zero. It must NOT be read as 'unset' and replaced by the mean."""
    _, region_arrays, substrate, chain, _ = parts
    dead = JunctionGeometry(
        src_node=np.array([0], dtype=np.int64),
        dst_node=np.array([2], dtype=np.int64),
        strand=np.array([Strand.POS], dtype=np.int8),
        reach_lo=np.array([0.0]),
        reach_hi=np.array([1000.0]),
    )
    g = build_node_geometry(chain, substrate, region_arrays, dead, GDNA_PMF, RNA_PMF)
    edge_slots = np.flatnonzero(np.asarray(chain.kind) == EDGE)
    assert np.all(g.eff_junction[edge_slots, 0] == 0.0)


def test_a_divisor_of_ZERO_is_NOT_FLOORED(parts):
    """⛔ and `effective_length`'s own contract: an object with no
    opportunity must return 0 and the CALLER must treat 0 as 'no evidence'. The predecessor floored
    every divisor to ``_EPS``, which is the defect that produced densities of ~1e9 on the 12.4 % of
    fine-partition nodes where the contained effective length collapses to exactly 0.

    A 100 bp node cannot contain a 150 bp fragment: the opportunity is exactly zero.
    """
    _, region_arrays, substrate, chain, junctions = parts
    huge = spike_pmf(150)
    g = build_node_geometry(chain, substrate, region_arrays, junctions, huge, huge)
    node_slots = np.flatnonzero(np.asarray(chain.kind) == NODE)
    assert np.all(g.eff_gdna[node_slots] == 0.0)
    assert np.all(g.eff_rna[node_slots] == 0.0)


def test_the_divisors_differ_between_the_two_COMPONENTS(geometry, parts):
    """gDNA and RNA have different length distributions and therefore different opportunity at the same
    object. A single shared divisor would make the 2x2 deconvolution singular."""
    _, _, _, chain, _ = parts
    node_slots = np.flatnonzero(np.asarray(chain.kind) == NODE)
    assert not np.allclose(geometry.eff_gdna[node_slots], geometry.eff_rna[node_slots])


def test_a_MIXED_pmf_is_not_collapsed_to_its_mean(parts):
    """The divisor is ``E_f[placements]``, not ``placements(E_f[w])`` — the two differ whenever the
    opportunity is non-linear in ``w``, which is every contained frame. Enumerated on a two-point pmf.
    """
    _, region_arrays, substrate, chain, junctions = parts
    pmf = two_point_pmf(20, 150)  # mean 85; 150 does not fit in a 100 bp node, 20 does
    g = build_node_geometry(chain, substrate, region_arrays, junctions, pmf, pmf)
    node_slots = np.flatnonzero(np.asarray(chain.kind) == NODE)
    slot = int(node_slots[0])
    length = int(region_arrays.region_size_bp[int(np.asarray(chain.obj_idx)[slot])])
    assert g.eff_gdna[slot] == pytest.approx(brute_contained(length, pmf))
    assert g.eff_gdna[slot] != pytest.approx(brute_contained(length, spike_pmf(85)))


# ---------------------------------------------------------------------------
# 4. the junction -> edge incidence, and the transcript-strand keying
# ---------------------------------------------------------------------------


def test_a_junction_deposits_on_BOTH_the_donor_and_the_acceptor_line(geometry, parts):
    """A junction ``(src, dst)`` has its DONOR at the line to the right of ``src`` and its ACCEPTOR at
    the line to the left of ``dst``. Molecules leave the template at the first and arrive at the
    second, so both lines genuinely saw the flux — and the index states both explicitly, which is what
    replaces the old exon-bit guess.

    The fixture's junction runs node 0 -> node 2, so it lands on edge 0 (donor) and edge 1 (acceptor)
    — i.e. on BOTH edge slots.
    """
    payload, _, _, chain, _ = parts
    edge_slots = np.flatnonzero(np.asarray(chain.kind) == EDGE)
    flux = float(payload.sj_count[0].sum())
    assert flux > 0
    for slot in edge_slots:
        assert geometry.junction_count[slot, 0] == pytest.approx(flux)


def test_the_mature_flux_is_keyed_by_the_JUNCTIONS_OWN_STRAND_not_the_align_column(parts):
    """⭐ The 'sense is derived, never stored' rule, made executable. The accumulator's ``sj_count``
    columns are the **genome** strand the read aligned to; the junction's strand is a property of the
    ANNOTATION. A ``-`` junction's flux belongs to the ``-`` transcript column however its reads
    aligned — and the fixture's junction has 9 POS-aligned and 4 NEG-aligned reads, so a column-wise
    copy would put 4 fragments in the wrong place.
    """
    payload, region_arrays, substrate, chain, _ = parts
    neg = JunctionGeometry(
        src_node=np.array([0], dtype=np.int64),
        dst_node=np.array([2], dtype=np.int64),
        strand=np.array([Strand.NEG], dtype=np.int8),
        reach_lo=np.array([1000.0]),
        reach_hi=np.array([1000.0]),
    )
    g = build_node_geometry(chain, substrate, region_arrays, neg, GDNA_PMF, RNA_PMF)
    edge_slots = np.flatnonzero(np.asarray(chain.kind) == EDGE)
    total = float(payload.sj_count[0].sum())
    assert payload.sj_count[0, 0] != payload.sj_count[0, 1], (
        "the fixture must make the axes distinct"
    )
    for slot in edge_slots:
        assert g.junction_count[slot, 1] == pytest.approx(total)  # the - transcript column
        assert g.junction_count[slot, 0] == 0.0


def test_several_junctions_on_one_line_POOL_their_counts_AND_their_divisors(parts):
    """Two junctions sharing a donor line are two estimates of one rate, so the pooled statement is
    ``sum(count) / sum(E)`` — the ratio of sums, never the mean of ratios.
    ``rho_bg = sum(g)/sum(E)``). Averaging the divisors instead would mis-weight the deeper junction.
    """
    payload, region_arrays, substrate, chain, _ = parts
    two = JunctionGeometry(
        src_node=np.array([0, 0], dtype=np.int64),
        dst_node=np.array([2, 2], dtype=np.int64),
        strand=np.array([Strand.POS, Strand.POS], dtype=np.int8),
        reach_lo=np.array([1000.0, 30.0]),
        reach_hi=np.array([1000.0, 30.0]),
    )
    # the fixture payload has one junction row; give the second its own
    import dataclasses

    payload2 = dataclasses.replace(
        payload,
        sj_count=np.array([[9, 4], [5, 1]], dtype=np.uint32),
        sj_inv_length_sum=np.zeros((2, 2), dtype=np.uint64),
        sj_length_sum=np.zeros((2, 2), dtype=np.uint64),
        ref_sj_offsets=np.array([0, 2], dtype=np.int64),
    )
    sub2 = CalibrationSubstrate.from_payload(payload2, region_arrays)
    g = build_node_geometry(
        chain,
        substrate=sub2,
        region_arrays=region_arrays,
        junctions=two,
        gdna_fl_pmf=GDNA_PMF,
        rna_fl_pmf=RNA_PMF,
    )
    slot = int(np.flatnonzero(np.asarray(chain.kind) == EDGE)[0])
    assert g.junction_count[slot, 0] == pytest.approx(9 + 4 + 5 + 1)
    assert g.eff_junction[slot, 0] == pytest.approx(
        brute_crossing(RNA_PMF, 1e12, 1e12) + brute_crossing(RNA_PMF, 30.0, 30.0)
    )


def two_reference_parts(payload):
    """chr1 with 3 nodes / 2 edges, then chr2 with 2 nodes / 1 edge.

    ⭐ **The second reference is the point.** Slot ids run ``N E N E N`` per reference, so within the
    FIRST reference node ``i`` always sits at slot ``2i`` — which means a geometry that assumes that
    layout instead of reading the chain's own adjacency is indistinguishable from a correct one on any
    single-reference fixture. chr2's node 3 sits at slot 5, not slot 6. This helper is what makes that
    assumption observable; it was added after a perturbation (`slot_of_node = 2 * arange`) passed all
    22 tests.
    """
    import dataclasses

    import pandas as pd

    from rigel.calibration.region_arrays import RegionArrays

    p2 = dataclasses.replace(
        payload,
        ref_node_offsets=np.array([0, 3, 5], dtype=np.int64),
        ref_edge_offsets=np.array([0, 2, 3], dtype=np.int64),
        ref_sj_offsets=np.array([0, 0, 1], dtype=np.int64),
        n_refs=2,
        node_contained_count=np.vstack(
            [payload.node_contained_count, np.array([[1, 1], [2, 2]], np.uint32)]
        ),
        node_contained_inv_length_sum=np.vstack(
            [payload.node_contained_inv_length_sum, np.zeros((2, 2), np.uint64)]
        ),
        node_contained_length_sum=np.vstack(
            [payload.node_contained_length_sum, np.zeros((2, 2), np.uint64)]
        ),
        node_spanning_count=np.vstack([payload.node_spanning_count, np.zeros((2, 2), np.uint32)]),
        node_spanning_inv_length_sum=np.vstack(
            [payload.node_spanning_inv_length_sum, np.zeros((2, 2), np.uint64)]
        ),
        node_spanning_length_sum=np.vstack(
            [payload.node_spanning_length_sum, np.zeros((2, 2), np.uint64)]
        ),
        node_start_count=np.concatenate([payload.node_start_count, np.zeros(2, np.uint32)]),
        edge_unspliced_count=np.vstack(
            [payload.edge_unspliced_count, np.array([[3, 3]], np.uint32)]
        ),
        edge_unspliced_inv_length_sum=np.vstack(
            [payload.edge_unspliced_inv_length_sum, np.zeros((1, 2), np.uint64)]
        ),
        edge_unspliced_length_sum=np.vstack(
            [payload.edge_unspliced_length_sum, np.zeros((1, 2), np.uint64)]
        ),
        edge_spliced_count=np.vstack([payload.edge_spliced_count, np.zeros((1, 2), np.uint32)]),
        edge_spliced_inv_length_sum=np.vstack(
            [payload.edge_spliced_inv_length_sum, np.zeros((1, 2), np.uint64)]
        ),
        edge_spliced_length_sum=np.vstack(
            [payload.edge_spliced_length_sum, np.zeros((1, 2), np.uint64)]
        ),
        cut_positions=np.array([0, 100, 200, 300, 0, 100, 200], dtype=np.int64),
        ref_cut_offsets=np.array([0, 4, 7], dtype=np.int64),
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
    chain2 = build_node_chain(p2.ref_node_offsets, p2.ref_edge_offsets)
    return p2, ra2, sub2, chain2


def _edge_slot_of(chain, edge_obj_id: int) -> int:
    edge_slots = np.flatnonzero(np.asarray(chain.kind) == EDGE)
    obj = np.asarray(chain.obj_idx)[edge_slots]
    return int(edge_slots[obj == edge_obj_id][0])


def test_a_junction_never_lands_on_ANOTHER_REFERENCE(parts):
    """Getting the incidence wrong shifts a junction onto a neighbouring reference's lines — invisible
    in aggregate, and exactly the class of defect the substrate's per-reference offset check exists
    for."""
    payload, _, _, _, _ = parts
    p2, ra2, sub2, chain2 = two_reference_parts(payload)
    j = JunctionGeometry(
        src_node=np.array([0], dtype=np.int64),  # chr1's first node
        dst_node=np.array([2], dtype=np.int64),
        strand=np.array([Strand.POS], dtype=np.int8),
        reach_lo=np.array([1000.0]),
        reach_hi=np.array([1000.0]),
    )
    g = build_node_geometry(chain2, sub2, ra2, j, GDNA_PMF, RNA_PMF)
    assert g.junction_count[_edge_slot_of(chain2, 2)].sum() == 0.0, "a chr1 junction reached chr2"


def test_a_junction_on_a_LATER_REFERENCE_lands_on_that_references_own_line(parts):
    """⛔ **The layout-assumption killer.** Node ``i`` sits at slot ``2i`` only within the FIRST
    reference; chr2's node 3 sits at slot 5. A geometry that computes the slot arithmetically instead
    of reading ``chain.left``/``chain.right`` passes every single-reference test and then puts chr2's
    mature flux on the wrong slot entirely — here on a NODE, which cannot carry mature at all.

    Found by perturbation: ``slot_of_node = 2 * arange(n_nodes)`` passed all 22 tests before this one
    existed.
    """
    payload, _, _, _, _ = parts
    p2, ra2, sub2, chain2 = two_reference_parts(payload)
    j = JunctionGeometry(
        src_node=np.array([3], dtype=np.int64),  # chr2's FIRST node
        dst_node=np.array([4], dtype=np.int64),  # chr2's second node
        strand=np.array([Strand.POS], dtype=np.int8),
        reach_lo=np.array([1000.0]),
        reach_hi=np.array([1000.0]),
    )
    g = build_node_geometry(chain2, sub2, ra2, j, GDNA_PMF, RNA_PMF)
    flux = float(p2.sj_count[0].sum())
    chr2_edge = _edge_slot_of(chain2, 2)
    # donor (right of node 3) and acceptor (left of node 4) are the SAME line, so it takes the flux twice
    assert g.junction_count[chr2_edge, 0] == pytest.approx(2.0 * flux)
    assert np.all(g.junction_count[np.asarray(chain2.kind) == NODE] == 0.0)
    for e in (0, 1):  # chr1's two lines saw nothing
        assert g.junction_count[_edge_slot_of(chain2, e)].sum() == 0.0


# ---------------------------------------------------------------------------
# 5. the derived per-node quantities
# ---------------------------------------------------------------------------


def test_node_global_geometry_no_longer_SUMS_TWO_FACES(geometry, parts):
    """The predecessor returned ``mass_l + mass_r`` over ``E_l + E_r`` at a boundary, and the docstring
    carried a long note about a ``1/2`` that was cancelling a missing ``1/2``. With one set of numbers
    there is nothing to sum and the note goes with it."""
    _, _, _, chain, _ = parts
    mass, eff = node_global_geometry(geometry)
    np.testing.assert_array_equal(mass, geometry.unspliced_count.sum(axis=1))
    np.testing.assert_array_equal(eff, geometry.eff_gdna)


def test_node_total_density_is_the_SUM_of_component_densities_each_in_its_own_frame(
    geometry, parts
):
    """`unified_solver_design` §2: rho = f_g*(M/E_g) + (1-f_g)*(M/E_r), each component in its OWN
    length frame, plus the one-sided mature density. Never one shared divisor."""
    _, _, _, chain, _ = parts
    f_g = np.full(chain.n_slots, 0.25)
    rho_u, rho_w = node_total_density(geometry, f_g)
    mass = geometry.unspliced_count.sum(axis=1)
    expected = mass * (
        0.25 / np.where(geometry.eff_gdna > 0, geometry.eff_gdna, np.inf)
        + 0.75 / np.where(geometry.eff_rna > 0, geometry.eff_rna, np.inf)
    )
    np.testing.assert_allclose(rho_u, expected)
    assert np.all(rho_w >= rho_u)


def test_a_zero_opportunity_object_emits_ZERO_DENSITY_not_infinity(parts):
    """⛔ The other half of trap 23. With the divisor no longer floored, the consumer must divide by
    zero *safely* — a node that cannot hold the component has no density, not an infinite one."""
    _, region_arrays, substrate, chain, junctions = parts
    huge = spike_pmf(150)
    g = build_node_geometry(chain, substrate, region_arrays, junctions, huge, huge)
    rho_u, rho_w = node_total_density(g, np.full(chain.n_slots, 0.5))
    assert np.all(np.isfinite(rho_u))
    assert np.all(np.isfinite(rho_w))
    node_slots = np.asarray(chain.kind) == NODE
    assert np.all(rho_u[node_slots] == 0.0)


def test_mature_density_pools_the_two_TRANSCRIPT_strands(geometry, parts):
    """rho_mature = sum over transcript strands of count/E, each strand in its own frame — and a strand
    with no flux contributes nothing rather than a 0/0."""
    _, _, _, chain, _ = parts
    rho_u, rho_w = node_total_density(geometry, np.zeros(chain.n_slots))
    expected = np.zeros(chain.n_slots)
    for s in (0, 1):
        c, e = geometry.junction_count[:, s], geometry.eff_junction[:, s]
        expected += np.where((c > 0) & (e > 0), c / np.where(e > 0, e, 1.0), 0.0)
    np.testing.assert_allclose(rho_w - rho_u, expected)


def test_spliced_count_and_junction_count_are_DIFFERENT_POPULATIONS(geometry, parts):
    """⛔ **The naming hazard this schema exists to remove.** Both populations are certified RNA and both
    could be called "mature", so one word would cover both and distinguish neither. They are different
    molecules:

    * ``spliced_count`` — crossed this line CONTIGUOUSLY, having spliced somewhere else;
    * ``junction_count`` — never crossed it, it JUMPED from here.

    The fixture makes them differ on both axes at once. At edge 0 the junction's donor sits on the line
    but nothing crossed it contiguously (13 vs 0); at edge 1 both are live and unequal (13 vs 6). A
    consumer that read one for the other would be off by the whole gene's mature output at a donor seam
    — which is the 2-versus-251 case in the design log.
    """
    payload, _, _, chain, _ = parts
    edge_slots = np.flatnonzero(np.asarray(chain.kind) == EDGE)
    spliced = geometry.spliced_count[edge_slots].sum(axis=1)
    junction = geometry.junction_count[edge_slots].sum(axis=1)
    np.testing.assert_array_equal(spliced, payload.edge_spliced_count.sum(axis=1))
    assert junction[0] == pytest.approx(float(payload.sj_count[0].sum()))
    # edge 0: junction flux with NO contiguous spliced crossing — they cannot be the same array
    assert spliced[0] == 0.0 and junction[0] > 0.0
    # edge 1: both live, and unequal
    assert spliced[1] > 0.0 and junction[1] > 0.0 and spliced[1] != junction[1]


def test_the_word_MATURE_names_no_field(geometry):
    """It fits both spliced populations, so it cannot be the name of either.
    trap 27 — one word on two concepts). The fields carry the accumulator's own three bank names."""
    fields = set(NodeGeometry.__dataclass_fields__)
    assert not any("mature" in f for f in fields), fields
    assert {"unspliced_count", "spliced_count", "junction_count"} <= fields
