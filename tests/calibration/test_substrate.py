"""CalibrationSubstrate: one type, the payload's own axes, and what those axes conserve.

This module tests the substrate the solver is handed (the populations each on its own axis, the derived per-crossing
quantities) and the alignment guard that refuses a geometry whose object count or per-reference offsets have drifted
from the payload's; a population read against the wrong axis would drop nearly every fragment inside ``deposit()`` with
every golden test green. The last block runs on a real scan, so the payload, the index and the geometry must agree with
each other, and its invariants are re-derived from the fragment buffer and the index, sources independent of the
accumulator.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.errors import CalibrationSubstrateError
from rigel.calibration.region_arrays import RegionArrays, boundary_region_indices
from rigel.calibration.substrate import CalibrationSubstrate, PopulationView
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario

from _synthetic import make_synthetic_payload


@pytest.fixture
def substrate():
    payload, region_arrays = make_synthetic_payload()
    return CalibrationSubstrate.from_payload(payload, region_arrays), payload, region_arrays


# ---------------------------------------------------------------------------
# the four populations, on their own axes
# ---------------------------------------------------------------------------


def test_every_population_is_present_on_the_right_axis(substrate):
    """Regions, contiguous boundaries and sj are three axes off by one per reference; a population read against the
    wrong axis would drop nearly every fragment with every golden test green, so the shapes are asserted.
    """
    sub, payload, _ = substrate
    # The populations do not all carry the same channels: a channel is stored where a named consumer reads it. ``None``
    # means "this population does not measure that", which differs from "it measured zero", and the view keeps them
    # distinguishable.
    for view, n, has_mass in (
        (sub.region_contained, payload.n_regions, False),
        (sub.boundary_unspliced, payload.n_boundaries, True),
        (sub.boundary_spliced, payload.n_boundaries, True),
        (sub.sj, payload.n_sj, True),
    ):
        assert view.count.shape == (n, 2)
        if has_mass:
            # One column: a mass carries no strand axis, while ``count`` keeps two. ``sj_mass`` arrives from the payload
            # with two columns and is folded at this boundary, which this shape assertion pins for the sj row.
            assert view.mass is not None and view.mass.shape == (n,)
        else:
            assert view.mass is None, (
                f"{view.name}.mass must be None, not zeros — a zero array cannot be told "
                f"apart from a real measurement of nothing"
            )
    assert payload.n_regions != payload.n_boundaries, "the fixture must not let an axis mix-up pass"


def test_the_columns_are_genome_strand_and_nothing_is_re_oriented(substrate):
    """One convention. Sense/antisense is transcript-relative, derived by the consumer from a sj's own strand, and never
    stored, so no field here may be named for it."""
    sub, payload, _ = substrate
    np.testing.assert_array_equal(sub.region_contained.count, payload.region_contained_count)
    for name in dir(sub):
        assert "sense" not in name, (
            f"{name} names a transcript-relative concept the schema does not store"
        )


def test_no_population_is_a_view_of_another(substrate):
    """Two populations holding the same numbers under two names cannot be told apart downstream, so nothing here may
    alias another population's array."""
    sub, _, _ = substrate
    banks = [sub.region_contained, sub.boundary_unspliced, sub.boundary_spliced]
    totals = [int(b.count.sum()) for b in banks]
    assert len(set(totals)) == len(totals), "the fixture gives every bank a distinct total"


# ---------------------------------------------------------------------------
# derived quantities
# ---------------------------------------------------------------------------


def test_the_sj_mass_arrives_per_strand_and_is_folded_here(substrate):
    """``sj_mass`` arrives per strand for artifact detection, and the strand axis stops here:
    :attr:`PopulationView.mass` is strand-agnostic (the mass turns an object-incidence total into a fragment count, a
    question with no strand in it). The fixture's two columns are unequal (0.9 / 0.4), so a fold that took one column,
    the max or the mean cannot pass.
    """
    sub, payload, _ = substrate
    assert payload.sj_mass.ndim == 2, "the payload bank is per strand"
    assert payload.sj_mass.shape[1] == 2
    assert payload.sj_mass[0, 0] != payload.sj_mass[0, 1], (
        "the fixture cannot separate the fold rules"
    )
    assert sub.sj.mass.ndim == 1, "PopulationView.mass is strand-agnostic"
    np.testing.assert_allclose(sub.sj.mass, payload.sj_mass.sum(axis=1))


def test_mass_per_crossing_at_zero_count_is_the_identity_not_zero():
    """An object the accumulator never saw must emit nothing, never a floored value: a "no data" default of 100 % gDNA
    would seed false gDNA into neighbouring exons. For a mass the null answer is 1.0 rather than NaN, because this
    factor rescales whatever mass the deconvolution placed at the boundary: 0 would delete it, 1.0 leaves it alone.
    """
    view = PopulationView(
        name="sj",
        count=np.zeros((2, 2), np.int64),
        mass=np.zeros(2, np.float64),
    )
    np.testing.assert_array_equal(view.mass_per_crossing, np.ones(2))


def test_total_count_sums_the_two_strands(substrate):
    sub, payload, _ = substrate
    np.testing.assert_array_equal(
        sub.region_contained.total_count, payload.region_contained_count.sum(axis=1)
    )


# ---------------------------------------------------------------------------
# the alignment check: a geometry must describe the payload it is handed with
# ---------------------------------------------------------------------------


def test_a_geometry_with_the_wrong_object_count_is_refused():
    """A geometry/payload mismatch otherwise surfaces as a shape error deep in the solver, pointing
    nowhere near its cause."""
    payload, region_arrays = make_synthetic_payload()
    import dataclasses

    trimmed = dataclasses.replace(
        region_arrays,
        start=region_arrays.start[:2],
        end=region_arrays.end[:2],
        ref_id=region_arrays.ref_id[:2],
        signature=region_arrays.signature[:2],
        region_size_bp=region_arrays.region_size_bp[:2],
    )
    with pytest.raises(CalibrationSubstrateError, match="2 objects but the payload has 3"):
        CalibrationSubstrate.from_payload(payload, trimmed)


def test_a_geometry_with_drifted_per_reference_offsets_is_refused():
    """The count can match while the per-reference slicing does not, which would drop nearly every real fragment inside
    ``deposit()``, so matching totals are not sufficient."""
    import dataclasses

    payload, region_arrays = make_synthetic_payload()
    drifted = dataclasses.replace(region_arrays, ref_offsets=np.array([0, 2], dtype=np.int32))
    with pytest.raises(CalibrationSubstrateError, match="per-reference offsets"):
        CalibrationSubstrate.from_payload(payload, drifted)


# ── the region geometry ↔ index/payload alignment guard ───────────────────────────────────────


def _region_arrays(mini_index):
    """The live geometry: the region partition the scanner deposits into."""
    return RegionArrays.from_index(mini_index)


def test_region_arrays_align_with_index(mini_index):
    ra = _region_arrays(mini_index)

    # R_total == len(regions_df).
    assert ra.n_regions == len(mini_index.regions_df)

    # Per-ref offsets equal the grouped region counts in ref_names order.
    counts = [int((mini_index.regions_df["ref_name"] == ref).sum()) for ref in mini_index.ref_names]
    expected = np.concatenate([[0], np.cumsum(counts)]).astype(np.int32)
    np.testing.assert_array_equal(ra.ref_offsets, expected)


def _payload(n_regions, ref_region_offsets):
    """The two fields the alignment guard reads. Both are on the region axis: the payload's boundary and sj axes are
    sized from it (``E = N − n_refs``), so a region-axis mismatch must be caught at the door."""
    return SimpleNamespace(
        n_regions=n_regions, ref_region_offsets=np.asarray(ref_region_offsets, np.int64)
    )


def test_alignment_guard_accepts_matching_payload(mini_index):
    ra = _region_arrays(mini_index)
    # Should not raise.
    CalibrationSubstrate.check_alignment(_payload(ra.n_regions, ra.ref_offsets), ra)


def test_alignment_guard_rejects_count_mismatch(mini_index):
    ra = _region_arrays(mini_index)
    with pytest.raises(CalibrationSubstrateError):
        CalibrationSubstrate.check_alignment(_payload(ra.n_regions + 1, ra.ref_offsets), ra)


def test_alignment_guard_rejects_offset_mismatch(mini_index):
    ra = _region_arrays(mini_index)
    bad = ra.ref_offsets.astype(np.int64).copy()
    bad[1] += 1  # perturb an offset so it no longer matches the geometry
    with pytest.raises(CalibrationSubstrateError):
        CalibrationSubstrate.check_alignment(_payload(ra.n_regions, bad), ra)


# ── conservation on a real scanned payload: the three axes, end to end ────────────────────────


SEED = 4321


@pytest.fixture(scope="module")
def scanned(tmp_path_factory):
    work = tmp_path_factory.mktemp("subcons")
    sc = Scenario("subcons", genome_length=6000, seed=SEED, work_dir=work / "subcons")
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(300, 600), (900, 1200)], "abundance": 60}])
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(3000, 3300), (3700, 4000)], "abundance": 40}])
    # gDNA is needed here: without it every fragment is mature RNA, which either sits inside one exon region or jumps,
    # so ``boundary_unspliced``, the population calibration deconvolves, would be identically zero and a conservation
    # test would pass over an empty bank.
    from rigel.sim.reads import GDNAConfig

    result = sc.build_oracle(
        n_fragments=300,
        sim_config=ReadSimConfig(frag_mean=200, frag_std=30, seed=SEED),
        gdna_config=GDNAConfig(abundance=400.0, frag_mean=250, frag_std=60),
    )
    config = PipelineConfig(em=EMConfig(seed=SEED), scan=BamScanConfig(sj_strand_tag="auto"))
    _, _, buffer, payload = scan_and_buffer(str(result.bam_path), result.index, config.scan)
    ra = RegionArrays.from_frame(result.index.regions_df, result.index.ref_name_to_id)
    yield payload, ra, buffer, result.index
    sc.cleanup()


def test_one_region_start_per_accepted_fragment(scanned):
    """``region_start_count`` is incremented once, at the region holding the fragment's first base, for every fragment
    the accumulator accepts, so its total is the accepted count. It is checked against the payload's QC tally, a
    separate counter written at the same point in the deposit.
    """
    payload, _ra, _buffer, _index = scanned
    assert int(np.asarray(payload.region_start_count).sum()) == int(payload.qc.deposited)
    assert int(payload.qc.deposited) > 0, "a scan that deposited nothing proves nothing"


def test_every_buffered_fragment_is_accepted_or_dropped_for_a_named_reason(scanned):
    """The accepted count plus every drop reason must account for the whole buffer, an independent source filled by the
    scanner's own fragment grouping."""
    payload, _ra, buffer, _index = scanned
    qc = payload.qc
    dropped = (
        qc.dropped_too_long
        + qc.dropped_empty
        + qc.dropped_strand_undefined
        + qc.deferred_undetermined_gap
    )
    assert qc.deposited + dropped == buffer.total_fragments


def test_the_boundary_axis_is_n_minus_the_nonempty_references(scanned):
    """``E = N − (references that own at least one region)``, re-derived from the index rather than the payload's own
    offsets. A reference with one region owns no boundary; a reference with none owns nothing; neither is a special
    case."""
    payload, ra, _buffer, index = scanned
    regions_per_ref = np.diff(np.asarray(ra.ref_offsets, dtype=np.int64))
    expected_boundaries = int(np.maximum(regions_per_ref - 1, 0).sum())
    assert int(payload.n_boundaries) == expected_boundaries
    assert int(payload.n_regions) == len(index.regions_df)


def test_contained_deposits_never_exceed_the_accepted_fragments(scanned):
    """A fragment lies wholly inside at most one region, so the contained bank cannot out-count the
    fragments. It is an inequality, not an identity: a fragment that crosses any boundary is contained
    in nothing and contributes only to the boundary banks."""
    payload, ra, _buffer, _index = scanned
    sub = CalibrationSubstrate.from_payload(payload, ra)
    contained = int(np.asarray(sub.region_contained.count).sum())
    assert 0 < contained <= int(np.asarray(payload.region_start_count).sum())


def test_no_boundary_straddles_a_reference(scanned):
    """Every contiguous boundary joins two regions of the same reference, an invariant a ``k + 1`` axis could not state."""
    _payload, ra, _buffer, _index = scanned
    lo, hi = boundary_region_indices(np.asarray(ra.ref_id))
    ref = np.asarray(ra.ref_id)
    np.testing.assert_array_equal(ref[lo], ref[hi])
    np.testing.assert_array_equal(hi, lo + 1)


def test_the_two_boundary_banks_are_disjoint_populations(scanned):
    """``boundary_unspliced`` and ``boundary_spliced`` are different molecules at the same boundary (crossed
    contiguously having spliced nowhere, versus having spliced elsewhere), so a fragment lands in exactly one. Checked
    on real data as a non-degeneracy check: the two banks cannot be the same array."""
    payload, ra, _buffer, _index = scanned
    sub = CalibrationSubstrate.from_payload(payload, ra)
    unspliced = np.asarray(sub.boundary_unspliced.count)
    spliced = np.asarray(sub.boundary_spliced.count)
    assert unspliced.sum() > 0, "the fixture must exercise the unspliced bank"
    assert not np.array_equal(unspliced, spliced)


def test_the_sj_axis_matches_the_payloads_own(scanned):
    """The substrate's sj population is exactly ``n_sj`` rows: the third axis, which a consumer must not size from
    ``n_regions`` or ``n_boundaries``."""
    payload, ra, _buffer, _index = scanned
    sub = CalibrationSubstrate.from_payload(payload, ra)
    assert sub.n_sj == int(payload.n_sj)
    assert np.asarray(sub.sj.count).shape == (int(payload.n_sj), 2)
    assert np.asarray(sub.sj.count).sum() > 0, "the fixture must exercise the sj axis"
