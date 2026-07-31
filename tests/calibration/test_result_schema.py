"""CalibrationResult.__post_init__ intrinsic invariants — the THREE-AXIS schema (S5.f).

⭐ **Three axes, and they are different lengths on purpose.** ``n_nodes``, ``n_edges`` and
``n_junctions`` are independent (``E = N − n_refs``, ``J`` unrelated to either), so every fixture here
uses three DIFFERENT lengths. A fixture that used one length for all three could not tell an
axis mix-up from a correct result, which is exactly how the predecessor's per-region ``left``/``right``
pair survived being pooled straight back together.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig

N_NODES, N_EDGES, N_JUNCTIONS = 4, 3, 2


def _valid_kwargs() -> dict:
    node = np.ones(N_NODES, dtype=np.float64)
    edge = np.ones(N_EDGES, dtype=np.float64)
    junction = np.ones(N_JUNCTIONS, dtype=np.float64)
    return dict(
        mass_gdna_node=np.zeros(N_NODES),
        mass_rna_node=node.copy(),
        mass_gdna_edge=np.zeros(N_EDGES),
        mass_rna_edge=edge.copy(),
        mass_rna_spliced_edge=np.zeros(N_EDGES),
        mass_rna_junction=junction.copy(),
        gdna_node_eff_len=node.copy(),
        gdna_edge_eff_len=edge.copy(),
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
        ("mass_gdna_edge", N_EDGES),
        ("mass_rna_edge", N_EDGES),
        ("mass_rna_spliced_edge", N_EDGES),
        ("gdna_edge_eff_len", N_EDGES),
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
    ["mass_gdna_node", "gdna_node_eff_len", "mass_gdna_edge", "gdna_edge_eff_len"],
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
    (`CARRY_FORWARD.md` §3 trap 2).
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
