"""``CalibrationResult.__post_init__`` — the intrinsic invariants of the three-axis schema.

``n_regions``, ``n_boundaries`` and ``n_sj`` are independent axes (``E = N − n_refs``, and ``J`` is
unrelated to either), so every fixture here uses three different lengths on purpose: with one
length for all three a fixture cannot tell an axis mix-up from a correct result. The gates hold
each array's shape against its own axis, the dtype rule, the value bounds on the three-way
composition, the deliberate absence of a closure assertion, the per-face names the schema must not
grow back, and the conserved sj mass — a derived property rather than a stored field, so an arm
that replaces the incidence array cannot leave a stale mass beside it.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig

N_REGIONS, N_BOUNDARIES, N_SJ = 4, 3, 2


def _valid_kwargs() -> dict:
    region = np.ones(N_REGIONS, dtype=np.float64)
    boundary = np.ones(N_BOUNDARIES, dtype=np.float64)
    return dict(
        count_gdna_region=np.zeros(N_REGIONS),
        count_rna_region=region.copy(),
        count_gdna_boundary=np.zeros(N_BOUNDARIES),
        count_rna_boundary=boundary.copy(),
        count_rna_spliced_boundary=np.zeros(N_BOUNDARIES),
        # Geometry, not a split: the mean conserved fragment-mass one crossing carries. 1.0 is the
        # identity — a boundary whose flanks both exceed every fragment length, where an incidence IS
        # a fragment — so a fixture that does not exercise K-inflation states it explicitly.
        boundary_mass_per_crossing=np.ones(N_BOUNDARIES),
        # These two must not be equal, and a `ones`/`ones` pair would make them so. `count_rna_sj`
        # is an incidence count and `sj_mass_per_crossing` converts it to the conserved mass; with
        # the conversion at the identity no gate in this file could tell the two quantities apart, so
        # the fixture could not fail the one confusion the pair exists to prevent
        # (`TRAPS: could-the-arm-have-fired` — a fixture is an arm). Conserved mass = [2.0, 1.5].
        count_rna_sj=np.array([4.0, 6.0]),
        boundary_spliced_mass_per_crossing=boundary.copy(),
        sj_mass_per_crossing=np.array([0.5, 0.25]),
        gdna_region_eff_len=region.copy(),
        gdna_boundary_eff_len=boundary.copy(),
        rna_region_eff_len=region.copy(),
        rna_boundary_eff_len=boundary.copy(),
        # The three-way composition ψ solves, per object. Not renormalised — it fails to close on a
        # substantial minority of both axes on real data, so a fixture that pretends otherwise would
        # be asserting something the shipped solver does not produce.
        gdna_frac_region=np.zeros(N_REGIONS),
        rna_pos_frac_region=region.copy(),
        rna_neg_frac_region=np.zeros(N_REGIONS),
        gdna_frac_boundary=np.zeros(N_BOUNDARIES),
        rna_pos_frac_boundary=boundary.copy(),
        rna_neg_frac_boundary=np.zeros(N_BOUNDARIES),
        gdna_density_global=1e-3,
        gdna_reference_density=None,
        gdna_reference_members=0,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=N_REGIONS,
        n_boundaries=N_BOUNDARIES,
        n_sj=N_SJ,
        config=CalibrationConfig(),
    )


def test_valid_result_constructs():
    CalibrationResult(**_valid_kwargs())


def test_zero_gdna_library_constructs():
    # Graceful zero-gDNA: gdna_density_global = 0 and all gDNA mass = 0 must be valid, not a failure.
    kw = _valid_kwargs()
    kw["gdna_density_global"] = 0.0
    kw["count_gdna_region"] = np.zeros(N_REGIONS)
    kw["count_gdna_boundary"] = np.zeros(N_BOUNDARIES)
    CalibrationResult(**kw)


def test_a_library_with_no_sj_constructs():
    """``J = 0`` is legal — a single-exon-only reference has no sj boundary at all — and must not
    be confused with "no sj flux"."""
    kw = _valid_kwargs()
    kw["count_rna_sj"] = np.zeros(0)
    kw["sj_mass_per_crossing"] = np.zeros(0)
    kw["n_sj"] = 0
    assert CalibrationResult(**kw).count_rna_sj.shape == (0,)


# ---------------------------------------------------------------------------
# Each array is pinned to ITS OWN axis.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "field,n_expected",
    [
        ("count_gdna_region", N_REGIONS),
        ("count_rna_region", N_REGIONS),
        ("gdna_region_eff_len", N_REGIONS),
        ("rna_region_eff_len", N_REGIONS),
        ("count_gdna_boundary", N_BOUNDARIES),
        ("count_rna_boundary", N_BOUNDARIES),
        ("count_rna_spliced_boundary", N_BOUNDARIES),
        ("gdna_boundary_eff_len", N_BOUNDARIES),
        ("rna_boundary_eff_len", N_BOUNDARIES),
        ("count_rna_sj", N_SJ),
    ],
)
def test_every_array_is_pinned_to_its_own_axis(field, n_expected):
    """The one defect this gate exists for: an array keyed to the wrong axis. With E = N − n_refs
    the two lengths differ by only a handful genome-wide, so a mis-keyed array is a plausible shape
    and a shape check is the only thing that catches it before the numbers go silently wrong."""
    for wrong in (N_REGIONS, N_BOUNDARIES, N_SJ):
        if wrong == n_expected:
            continue
        kw = _valid_kwargs()
        kw[field] = np.ones(wrong, dtype=np.float64)
        with pytest.raises(ValueError, match=f"expected \\({n_expected},\\)"):
            CalibrationResult(**kw)


def test_the_error_names_the_axis_it_expected():
    kw = _valid_kwargs()
    kw["count_gdna_boundary"] = np.ones(N_REGIONS)
    with pytest.raises(ValueError, match="count_gdna_boundary"):
        CalibrationResult(**kw)


# ---------------------------------------------------------------------------
# Value invariants.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "field",
    [
        "count_gdna_region",
        "gdna_region_eff_len",
        "count_gdna_boundary",
        "gdna_boundary_eff_len",
        "rna_region_eff_len",
        "rna_boundary_eff_len",
    ],
)
def test_rejects_negative(field):
    kw = _valid_kwargs()
    arr = np.asarray(kw[field], dtype=np.float64).copy()
    arr[0] = -1.0
    kw[field] = arr
    with pytest.raises(ValueError, match="non-negative"):
        CalibrationResult(**kw)


def test_rejects_non_finite_array():
    kw = _valid_kwargs()
    kw["gdna_region_eff_len"] = np.array([np.inf, 1.0, 1.0, 1.0])
    with pytest.raises(ValueError, match="non-finite"):
        CalibrationResult(**kw)


def test_accepts_an_integer_count_array():
    """The accumulator's primary per-object observable is an integer count, so an exact integer
    array is a better input here than a float one. ``count_rna_sj`` is the sharpest case: it is
    the sj flux verbatim, never deconvolved, so it arrives integral."""
    for dtype in (np.int64, np.int32, np.uint32, np.uint64):
        kw = _valid_kwargs()
        kw["count_rna_sj"] = np.array([1, 1], dtype=dtype)
        assert CalibrationResult(**kw).count_rna_sj.dtype == dtype


def test_still_rejects_a_narrower_float():
    """Integers are exact; float32 is not. Admitting it would silently mix precisions through
    arithmetic that is float64 everywhere else, which is what the dtype gate is actually for."""
    kw = _valid_kwargs()
    kw["count_rna_region"] = np.ones(N_REGIONS, dtype=np.float32)
    with pytest.raises(ValueError, match="float64 or an integer count"):
        CalibrationResult(**kw)


def test_still_rejects_a_negative_integer_count():
    kw = _valid_kwargs()
    kw["count_rna_region"] = np.array([1, -1, 1, 1], dtype=np.int64)
    with pytest.raises(ValueError, match="non-negative"):
        CalibrationResult(**kw)


@pytest.mark.parametrize(
    "field,value",
    [
        ("gdna_density_global", -1.0),
        ("gdna_density_global", np.inf),
        ("rna_sense_frac", 1.5),
        ("rna_sense_frac", -0.1),
        ("gdna_strand_overdispersion", 1.0),
        ("rna_strand_overdispersion", -0.1),
    ],
)
def test_rejects_bad_scalars(field, value):
    kw = _valid_kwargs()
    kw[field] = value
    with pytest.raises(ValueError):
        CalibrationResult(**kw)


@pytest.mark.parametrize("field", ["n_regions", "n_boundaries", "n_sj"])
def test_rejects_negative_axis_length(field):
    kw = _valid_kwargs()
    kw[field] = -1
    with pytest.raises(ValueError):
        CalibrationResult(**kw)


# ---------------------------------------------------------------------------
# What the schema must NOT carry any more.
# ---------------------------------------------------------------------------


def test_the_per_face_fields_are_gone():
    """``gdna_boundary_len`` has no successor, and neither do the ``left``/``right`` mass pairs.

    ``gdna_boundary_len`` was ``E[min(ℓ,L)]/2``: a per-face divisor, halved because a boundary was
    treated as two sides that were then summed back together. A boundary is one object, so its
    divisor is the per-boundary ``gdna_boundary_eff_len`` — one number at a 0-bp boundary, with no ½
    in it. Anything naming the old fields is reading a convention the schema does not have.
    """
    # What is banned is the per-FACE convention, not a vocabulary: `gdna_region_eff_len` and
    # `n_regions` are live field names. Everything listed below names the convention — a divisor or a
    # mass belonging to one side of a boundary — rather than a word that fell out of fashion.
    fields = set(CalibrationResult.__dataclass_fields__)
    assert not fields & {
        "gdna_boundary_len",
        "mass_gdna_left",
        "mass_gdna_right",
        "mass_rna_left",
        "mass_rna_right",
        "mass_gdna_contained",
        "mass_rna_contained",
        "mass_rna_spliced",
    }


def test_count_rna_spliced_has_no_region_twin():
    """Structural, not an omission: the accumulator credits ``region_contained`` only when the fragment
    used NO sj, so a region's contained population cannot hold a spliced molecule. A
    ``count_rna_spliced_region`` field would be a channel that cannot exist."""
    assert "count_rna_spliced_boundary" in CalibrationResult.__dataclass_fields__
    assert "count_rna_spliced_region" not in CalibrationResult.__dataclass_fields__


# ── the three-way composition ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "name",
    [
        "gdna_frac_region",
        "rna_pos_frac_region",
        "rna_neg_frac_region",
        "gdna_frac_boundary",
        "rna_pos_frac_boundary",
        "rna_neg_frac_boundary",
    ],
)
def test_a_composition_component_above_one_is_REFUSED(name):
    """Each of the three is a fraction of its object's unspliced population, so it is bounded by 1.
    Non-negativity and finiteness come from the shared axis check; this is the upper bound, the half
    a "sum of shares" schema cannot get from the axis check alone."""
    kw = _valid_kwargs()
    arr = np.asarray(kw[name], dtype=np.float64).copy()
    arr[0] = 1.5
    kw[name] = arr
    with pytest.raises(ValueError, match="must not exceed 1"):
        CalibrationResult(**kw)


def test_a_composition_that_does_NOT_close_is_ACCEPTED_and_that_is_deliberate():
    """This pins a decision, not a behaviour. ``f_g + f_pos + f_neg`` fails to reach 1 on a
    substantial minority of both axes on real data, because ``sweep``'s write-back clips the three
    posterior means independently and an unsolvable slot keeps an init instead.

    ⛔ The schema therefore does not assert closure, and this test exists so that nobody adds the
    assertion without first fixing ψ, and nobody "repairs" the symptom by renormalising the arrays
    at publication — which would make a short-closing object indistinguishable from a solved one.
    """
    kw = _valid_kwargs()
    kw["gdna_frac_region"] = np.full(N_REGIONS, 0.25)
    kw["rna_pos_frac_region"] = np.full(N_REGIONS, 0.25)
    kw["rna_neg_frac_region"] = np.full(N_REGIONS, 0.25)  # sums to 0.75, not 1
    res = CalibrationResult(**kw)
    total = res.gdna_frac_region + res.rna_pos_frac_region + res.rna_neg_frac_region
    assert np.allclose(total, 0.75), "the composition was silently renormalised at publication"


# ---------------------------------------------------------------------------
# The conserved sj mass — the third axis's incidence→fragment conversion.
# ---------------------------------------------------------------------------


def test_the_conserved_sj_mass_is_the_incidence_TIMES_its_own_conversion():
    """The arithmetic, on a fixture where the two are far enough apart that the gate can fail."""
    res = CalibrationResult(**_valid_kwargs())
    np.testing.assert_allclose(res.sj_conserved_mass, [2.0, 1.5])


def test_the_sj_INCIDENCE_is_NOT_the_sj_MASS():
    """The confusion this property exists for: ``count_rna_sj`` is named like a mass and is an
    incidence count — a fragment deposits ``+1`` on every sj it uses, so on real data it runs well
    over twice the conserved mass. This fires the moment anyone "simplifies" the property to return
    the incidence array, which is the specific edit that looks correct.
    """
    res = CalibrationResult(**_valid_kwargs())
    incidence = np.asarray(res.count_rna_sj, dtype=np.float64)
    assert not np.allclose(res.sj_conserved_mass, incidence), (
        "sj_conserved_mass returned the INCIDENCE count — the conversion was dropped"
    )
    # And the direction is fixed, not merely different: an incidence over-counts, never under-counts,
    # because a fragment using K sj books K of them and one unit of mass.
    assert np.all(incidence >= res.sj_conserved_mass - 1e-12)


def test_a_sj_NOTHING_crossed_has_ZERO_conserved_mass_not_the_identity():
    """``mass_per_crossing`` is deliberately 1.0 where nothing crossed — the identity, so a
    deconvolution's mass at an unobserved boundary is rescaled by 1 rather than deleted. Multiplying
    it by the zero incidence is what turns that identity back into the ``0`` that is correct here.

    This is the gate that fires if the property ever grows a ``where(count > 0, …)`` branch falling
    back to the conversion factor, which would publish one unit of RNA mass at a sj no fragment ever
    used — a false positive on the axis that is certified RNA by construction. Zero-count sj are a
    large minority of the axis on real data, not a corner case, so that fallback would invent mass
    at scale.
    """
    kw = _valid_kwargs()
    kw["count_rna_sj"] = np.array([0.0, 6.0])
    kw["sj_mass_per_crossing"] = np.array([1.0, 0.25])
    res = CalibrationResult(**kw)
    assert res.sj_conserved_mass[0] == 0.0
    np.testing.assert_allclose(res.sj_conserved_mass[1], 1.5)


def test_library_rna_fragments_READS_the_property_so_there_is_ONE_home():
    """``library_rna_fragments`` must read the property rather than re-spell the conversion. Two
    spellings of one conversion is how a caller reading the property and a caller reading the
    library count come to disagree, so this pins that moving one moves the other.
    """
    kw = _valid_kwargs()
    before = CalibrationResult(**kw)
    kw["sj_mass_per_crossing"] = np.array([0.5, 0.25]) * 3.0
    after = CalibrationResult(**kw)
    delta_property = float(after.sj_conserved_mass.sum() - before.sj_conserved_mass.sum())
    delta_library = after.library_rna_fragments - before.library_rna_fragments
    assert delta_property > 0.0, (
        "the perturbation did not move the property — the arm could not fire"
    )
    np.testing.assert_allclose(delta_library, delta_property)


def test_the_conserved_mass_survives_the_oracle_arms_dataclass_replace():
    """Why it is a property and not a stored field. ``count_rna_sj`` is in
    ``prior_vs_oracle.OVERRIDE_FIELDS``: an arm swaps it with ``dataclasses.replace``. A stored array
    would survive that swap and go on describing the array it replaced —
    ``TRAPS: a-hash-that-misses-its-artifact`` in dataclass form.
    """
    res = CalibrationResult(**_valid_kwargs())
    swapped = dataclasses.replace(res, count_rna_sj=np.array([40.0, 60.0]))
    np.testing.assert_allclose(swapped.sj_conserved_mass, [20.0, 15.0])
    assert not np.allclose(swapped.sj_conserved_mass, res.sj_conserved_mass)


def test_the_reference_density_is_None_or_positive_and_finite():
    """`gdna_reference_density` is the fully-captured gDNA level or ``None`` (no enriched mode): a zero,
    a negative or a non-finite reference is refused, since the ruler divides by it."""
    kw = _valid_kwargs()
    assert CalibrationResult(**kw).gdna_reference_density is None
    kw["gdna_reference_density"] = 0.37
    kw["gdna_reference_members"] = 12
    assert CalibrationResult(**kw).gdna_reference_density == 0.37
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        kw["gdna_reference_density"] = bad
        with pytest.raises(ValueError):
            CalibrationResult(**kw)


def test_the_reference_members_count_the_kernels_behind_a_reference_and_are_ZERO_without_one():
    """`gdna_reference_members` is the regime: the located kernels the enriched mode rests on. It is
    positive exactly when a reference is present — a reference from no kernel and a member count with
    no reference are both refused, so a consumer reading the pair cannot see a half-published state."""
    kw = _valid_kwargs()
    assert CalibrationResult(**kw).gdna_reference_members == 0
    kw["gdna_reference_density"] = 0.37
    kw["gdna_reference_members"] = 12
    assert CalibrationResult(**kw).gdna_reference_members == 12
    kw["gdna_reference_members"] = 0
    with pytest.raises(ValueError, match="gdna_reference_members"):
        CalibrationResult(**kw)
    kw["gdna_reference_density"] = None
    kw["gdna_reference_members"] = 3
    with pytest.raises(ValueError, match="gdna_reference_members"):
        CalibrationResult(**kw)
