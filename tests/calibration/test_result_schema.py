"""CalibrationResult.__post_init__ intrinsic invariants — the THREE-AXIS schema (S5.f).

⭐ **Three axes, and they are different lengths on purpose.** ``n_nodes``, ``n_edges`` and
``n_junctions`` are independent (``E = N − n_refs``, ``J`` unrelated to either), so every fixture here
uses three DIFFERENT lengths. A fixture that used one length for all three could not tell an
axis mix-up from a correct result, which is exactly how the predecessor's per-region ``left``/``right``
pair survived being pooled straight back together.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig

N_NODES, N_EDGES, N_JUNCTIONS = 4, 3, 2


def _valid_kwargs() -> dict:
    node = np.ones(N_NODES, dtype=np.float64)
    edge = np.ones(N_EDGES, dtype=np.float64)
    return dict(
        mass_gdna_node=np.zeros(N_NODES),
        mass_rna_node=node.copy(),
        mass_gdna_edge=np.zeros(N_EDGES),
        mass_rna_edge=edge.copy(),
        mass_rna_spliced_edge=np.zeros(N_EDGES),
        # ⭐ GEOMETRY, not a split: the mean conserved fragment-mass one crossing carries. 1.0 is the
        # identity — a line whose flanks both exceed every fragment length, where an incidence IS
        # a fragment — so a fixture that does not exercise K-inflation states it explicitly.
        edge_mass_per_crossing=np.ones(N_EDGES),
        # ⛔ These two must NOT be equal, and a `ones`/`ones` pair is what they were. `mass_rna_junction`
        # is an INCIDENCE count and `junction_mass_per_crossing` converts it to the conserved mass; with
        # the conversion at the identity no gate in this file could tell the two quantities apart, so
        # the fixture could not fail the one confusion the pair exists to prevent
        # (`TRAPS: could-the-arm-have-fired` — a fixture is an arm). Conserved mass = [2.0, 1.5].
        mass_rna_junction=np.array([4.0, 6.0]),
        edge_spliced_mass_per_crossing=edge.copy(),
        junction_mass_per_crossing=np.array([0.5, 0.25]),
        gdna_node_eff_len=node.copy(),
        gdna_edge_eff_len=edge.copy(),
        rna_node_eff_len=node.copy(),
        rna_edge_eff_len=edge.copy(),
        # ⭐ The three-way composition ψ solves, per object. ⛔ NOT renormalised — see the field
        # docstring: it fails to close on ~25 % of both axes on real data, so a fixture that pretends
        # otherwise would be asserting something the shipped solver does not produce.
        gdna_frac_node=np.zeros(N_NODES),
        rna_pos_frac_node=node.copy(),
        rna_neg_frac_node=np.zeros(N_NODES),
        gdna_frac_edge=np.zeros(N_EDGES),
        rna_pos_frac_edge=edge.copy(),
        rna_neg_frac_edge=np.zeros(N_EDGES),
        gdna_density_global=1e-3,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_nodes=N_NODES,
        n_edges=N_EDGES,
        n_junctions=N_JUNCTIONS,
        config=CalibrationConfig(),
    )


def test_valid_result_constructs():
    CalibrationResult(**_valid_kwargs())


def test_zero_gdna_library_constructs():
    # Graceful zero-gDNA: gdna_density_global = 0 and all gDNA mass = 0 must be valid, not a failure.
    kw = _valid_kwargs()
    kw["gdna_density_global"] = 0.0
    kw["mass_gdna_node"] = np.zeros(N_NODES)
    kw["mass_gdna_edge"] = np.zeros(N_EDGES)
    CalibrationResult(**kw)


def test_a_library_with_no_junctions_constructs():
    """``J = 0`` is legal — a single-exon-only reference has no junction edge at all — and must not
    be confused with "no junction flux"."""
    kw = _valid_kwargs()
    kw["mass_rna_junction"] = np.zeros(0)
    kw["junction_mass_per_crossing"] = np.zeros(0)
    kw["n_junctions"] = 0
    assert CalibrationResult(**kw).mass_rna_junction.shape == (0,)


# ---------------------------------------------------------------------------
# Each array is pinned to ITS OWN axis.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "field,n_expected",
    [
        ("mass_gdna_node", N_NODES),
        ("mass_rna_node", N_NODES),
        ("gdna_node_eff_len", N_NODES),
        ("rna_node_eff_len", N_NODES),
        ("mass_gdna_edge", N_EDGES),
        ("mass_rna_edge", N_EDGES),
        ("mass_rna_spliced_edge", N_EDGES),
        ("gdna_edge_eff_len", N_EDGES),
        ("rna_edge_eff_len", N_EDGES),
        ("mass_rna_junction", N_JUNCTIONS),
    ],
)
def test_every_array_is_pinned_to_its_own_axis(field, n_expected):
    """⛔ The one defect this gate exists for: an array keyed to the WRONG axis. With E = N − n_refs
    the two lengths differ by only a handful genome-wide, so a mis-keyed array is a plausible shape
    and a shape check is the only thing that catches it before the numbers go silently wrong."""
    for wrong in (N_NODES, N_EDGES, N_JUNCTIONS):
        if wrong == n_expected:
            continue
        kw = _valid_kwargs()
        kw[field] = np.ones(wrong, dtype=np.float64)
        with pytest.raises(ValueError, match=f"expected \\({n_expected},\\)"):
            CalibrationResult(**kw)


def test_the_error_names_the_axis_it_expected():
    kw = _valid_kwargs()
    kw["mass_gdna_edge"] = np.ones(N_NODES)
    with pytest.raises(ValueError, match="mass_gdna_edge"):
        CalibrationResult(**kw)


# ---------------------------------------------------------------------------
# Value invariants.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "field",
    [
        "mass_gdna_node",
        "gdna_node_eff_len",
        "mass_gdna_edge",
        "gdna_edge_eff_len",
        "rna_node_eff_len",
        "rna_edge_eff_len",
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
    kw["gdna_node_eff_len"] = np.array([np.inf, 1.0, 1.0, 1.0])
    with pytest.raises(ValueError, match="non-finite"):
        CalibrationResult(**kw)


def test_accepts_an_integer_count_array():
    """⚠ The accumulator's primary per-object observable is an integer COUNT, so an exact integer
    array is a BETTER input here than a float one. ``mass_rna_junction`` is the sharpest case: it is
    the junction flux verbatim, never deconvolved, so it arrives integral."""
    for dtype in (np.int64, np.int32, np.uint32, np.uint64):
        kw = _valid_kwargs()
        kw["mass_rna_junction"] = np.array([1, 1], dtype=dtype)
        assert CalibrationResult(**kw).mass_rna_junction.dtype == dtype


def test_still_rejects_a_narrower_float():
    """Integers are exact; float32 is not. Admitting it would silently mix precisions through
    arithmetic that is float64 everywhere else, which is what the dtype gate is actually for."""
    kw = _valid_kwargs()
    kw["mass_rna_node"] = np.ones(N_NODES, dtype=np.float32)
    with pytest.raises(ValueError, match="float64 or an integer count"):
        CalibrationResult(**kw)


def test_still_rejects_a_negative_integer_count():
    kw = _valid_kwargs()
    kw["mass_rna_node"] = np.array([1, -1, 1, 1], dtype=np.int64)
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


@pytest.mark.parametrize("field", ["n_nodes", "n_edges", "n_junctions"])
def test_rejects_negative_axis_length(field):
    kw = _valid_kwargs()
    kw[field] = -1
    with pytest.raises(ValueError):
        CalibrationResult(**kw)


# ---------------------------------------------------------------------------
# What the schema must NOT carry any more.
# ---------------------------------------------------------------------------


def test_the_per_face_fields_are_gone():
    """⛔ ``gdna_boundary_len`` has NO successor, and neither do the ``left``/``right`` mass pairs.

    ``gdna_boundary_len`` was ``E[min(ℓ,L)]/2`` — a per-FACE divisor, halved because a boundary had two
    sides that were then summed back together. S5.c deleted the quantity and S5.e deleted the faces;
    its replacement is the per-edge ``gdna_edge_eff_len``, ONE number at a 0-bp line with no ½ in it.
    Anything still naming the old fields is reading a convention that no longer exists

    """
    fields = set(CalibrationResult.__dataclass_fields__)
    assert not fields & {
        "gdna_boundary_len",
        "gdna_region_eff_len",
        "mass_gdna_left",
        "mass_gdna_right",
        "mass_rna_left",
        "mass_rna_right",
        "mass_gdna_contained",
        "mass_rna_contained",
        "mass_rna_spliced",
        "n_regions",
    }


def test_mass_rna_spliced_has_no_node_twin():
    """⚠ Structural, not an omission: the accumulator credits ``node_contained`` only when the fragment
    used NO junction, so a node's contained population cannot hold a spliced molecule. A
    ``mass_rna_spliced_node`` field would be a channel that cannot exist."""
    assert "mass_rna_spliced_edge" in CalibrationResult.__dataclass_fields__
    assert "mass_rna_spliced_node" not in CalibrationResult.__dataclass_fields__


# ── the three-way composition ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "name",
    [
        "gdna_frac_node",
        "rna_pos_frac_node",
        "rna_neg_frac_node",
        "gdna_frac_edge",
        "rna_pos_frac_edge",
        "rna_neg_frac_edge",
    ],
)
def test_a_composition_component_above_one_is_REFUSED(name):
    """⭐ Each of the three is a FRACTION of its object's unspliced population, so it is bounded by 1.
    ⚠ Non-negativity and finiteness come from the shared axis check; this is the upper bound, which is
    the half a "sum of shares" schema cannot get from the axis check alone."""
    kw = _valid_kwargs()
    arr = np.asarray(kw[name], dtype=np.float64).copy()
    arr[0] = 1.5
    kw[name] = arr
    with pytest.raises(ValueError, match="must not exceed 1"):
        CalibrationResult(**kw)


def test_a_composition_that_does_NOT_close_is_ACCEPTED_and_that_is_deliberate():
    """⛔⛔ **PINNING A DECISION, NOT A BEHAVIOUR.** ``f_g + f_pos + f_neg`` fails to reach 1 on ~25 % of
    both axes on real data — measured 74.72 % / 77.24 % closing, the rest with median 0.978 and a p5 of
    0.869 — because ``sweep``'s write-back clips the three posterior means INDEPENDENTLY and an
    unsolvable slot keeps an init instead.

    ⭐ The schema therefore does NOT assert closure, and this test exists so that nobody adds the
    assertion without first fixing ψ, and nobody "repairs" the symptom by renormalising the arrays at
    publication — which would make a 15 %-short object indistinguishable from a solved one.
    """
    kw = _valid_kwargs()
    kw["gdna_frac_node"] = np.full(N_NODES, 0.25)
    kw["rna_pos_frac_node"] = np.full(N_NODES, 0.25)
    kw["rna_neg_frac_node"] = np.full(N_NODES, 0.25)  # sums to 0.75, not 1
    res = CalibrationResult(**kw)
    total = res.gdna_frac_node + res.rna_pos_frac_node + res.rna_neg_frac_node
    assert np.allclose(total, 0.75), "the composition was silently renormalised at publication"


# ---------------------------------------------------------------------------
# ⭐⭐ The CONSERVED junction mass — the third axis's incidence→fragment conversion.
# ---------------------------------------------------------------------------


def test_the_conserved_junction_mass_is_the_incidence_TIMES_its_own_conversion():
    """⭐ The arithmetic, on a fixture where the two are 2–4× apart so the gate can fail."""
    res = CalibrationResult(**_valid_kwargs())
    np.testing.assert_allclose(res.junction_conserved_mass, [2.0, 1.5])


def test_the_junction_INCIDENCE_is_NOT_the_junction_MASS():
    """⛔⛔ **THE TRAP THIS PROPERTY EXISTS FOR.** ``mass_rna_junction`` is named like a mass and is an
    incidence count — a fragment deposits ``+1`` on every junction it uses, measured **2.0719×** per
    unit of conserved mass on ``g00 ss0.99 capture_off``. This fires the moment anyone "simplifies" the
    property to return the incidence array, which is the specific edit that looks correct.
    """
    res = CalibrationResult(**_valid_kwargs())
    incidence = np.asarray(res.mass_rna_junction, dtype=np.float64)
    assert not np.allclose(res.junction_conserved_mass, incidence), (
        "junction_conserved_mass returned the INCIDENCE count — the conversion was dropped"
    )
    # ⭐ And the direction is fixed, not merely different: an incidence over-counts, never under-counts,
    # because a fragment using K junctions books K of them and one unit of mass.
    assert np.all(incidence >= res.junction_conserved_mass - 1e-12)


def test_a_junction_NOTHING_crossed_has_ZERO_conserved_mass_not_the_identity():
    """⛔ ``mass_per_crossing`` is deliberately **1.0** where nothing crossed — the identity, so a
    deconvolution's mass at an unobserved line is rescaled by 1 rather than deleted. Multiplying it by
    the zero incidence is what turns that identity back into the ``0`` that is correct here.

    ⚠ This is the gate that fires if the property ever grows a ``where(count > 0, …)`` branch that
    falls back to the conversion factor, which would publish **1.0 units of RNA mass at a junction no
    fragment ever used** — a false positive on the axis that is certified RNA by construction. ⛔ Not a
    corner case: **4,636 of 13,482** junctions are zero-count on ``g00 ss0.99 capture_off``, so that
    fallback would invent mass on a third of the axis.
    """
    kw = _valid_kwargs()
    kw["mass_rna_junction"] = np.array([0.0, 6.0])
    kw["junction_mass_per_crossing"] = np.array([1.0, 0.25])
    res = CalibrationResult(**kw)
    assert res.junction_conserved_mass[0] == 0.0
    np.testing.assert_allclose(res.junction_conserved_mass[1], 1.5)


def test_library_rna_fragments_READS_the_property_so_there_is_ONE_home():
    """⛔ The junction term used to be spelled out a second time inside ``library_rna_fragments``. Two
    spellings of one conversion is how a caller reading the property and a caller reading the library
    count come to disagree, so this pins that moving one moves the other.
    """
    kw = _valid_kwargs()
    before = CalibrationResult(**kw)
    kw["junction_mass_per_crossing"] = np.array([0.5, 0.25]) * 3.0
    after = CalibrationResult(**kw)
    delta_property = float(
        after.junction_conserved_mass.sum() - before.junction_conserved_mass.sum()
    )
    delta_library = after.library_rna_fragments - before.library_rna_fragments
    assert delta_property > 0.0, (
        "the perturbation did not move the property — the arm could not fire"
    )
    np.testing.assert_allclose(delta_library, delta_property)


def test_the_conserved_mass_survives_the_oracle_arms_dataclass_replace():
    """⭐⭐ **WHY IT IS A PROPERTY AND NOT A STORED FIELD.** ``mass_rna_junction`` is in
    ``prior_vs_oracle.OVERRIDE_FIELDS``: an arm swaps it with ``dataclasses.replace``. A stored array
    would survive that swap and go on describing the array it replaced —
    ``TRAPS: a-hash-that-misses-its-artifact`` in dataclass form.
    """
    res = CalibrationResult(**_valid_kwargs())
    swapped = dataclasses.replace(res, mass_rna_junction=np.array([40.0, 60.0]))
    np.testing.assert_allclose(swapped.junction_conserved_mass, [20.0, 15.0])
    assert not np.allclose(swapped.junction_conserved_mass, res.junction_conserved_mass)
