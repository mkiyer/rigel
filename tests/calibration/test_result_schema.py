"""CalibrationResult.__post_init__ intrinsic invariants (acyclic schema)."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig


def _valid_kwargs(n_regions: int = 2) -> dict:
    z = np.zeros(n_regions, dtype=np.float64)
    o = np.ones(n_regions, dtype=np.float64)
    return dict(
        mass_gdna_contained=z.copy(),
        mass_rna_contained=o.copy(),
        mass_gdna_left=z.copy(),
        mass_rna_left=z.copy(),
        mass_gdna_right=z.copy(),
        mass_rna_right=z.copy(),
        mass_rna_spliced=z.copy(),
        gdna_boundary_len=o.copy(),
        gdna_region_eff_len=o.copy(),
        gdna_density_global=1e-3,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n_regions,
        config=CalibrationConfig(),
    )


def test_valid_result_constructs():
    CalibrationResult(**_valid_kwargs())


def test_zero_gdna_library_constructs():
    # Graceful zero-gDNA (decision F): gdna_density_global = 0 and all gDNA mass = 0 must be valid.
    kw = _valid_kwargs()
    kw["gdna_density_global"] = 0.0
    for k in (
        "mass_gdna_contained",
        "mass_gdna_left",
        "mass_gdna_right",
    ):
        kw[k] = np.zeros(2, dtype=np.float64)
    CalibrationResult(**kw)


@pytest.mark.parametrize(
    "field,value",
    [
        ("mass_gdna_contained", np.array([-1.0, 0.0])),  # negative mass
        ("gdna_region_eff_len", np.array([1.0, -1.0])),  # negative effective support
    ],
)
def test_rejects_negative_region_arrays(field, value):
    kw = _valid_kwargs()
    kw[field] = value
    with pytest.raises(ValueError):
        CalibrationResult(**kw)


def test_rejects_non_finite_array():
    kw = _valid_kwargs()
    kw["gdna_region_eff_len"] = np.array([np.inf, 1.0])
    with pytest.raises(ValueError):
        CalibrationResult(**kw)


def test_accepts_an_integer_count_array():
    """⚠ The gate was float64-ONLY, and accumulator v5 makes the primary per-object observable an
    integer COUNT (spec §6: no floats in the data model). An exact integer array is a better input
    here than a float one, so it is admitted — and rejecting it would block every W5 consumer arm
    at its last step."""
    for dtype in (np.int64, np.int32, np.uint32, np.uint64):
        kw = _valid_kwargs()
        kw["mass_rna_contained"] = np.array([1, 1], dtype=dtype)
        assert CalibrationResult(**kw).mass_rna_contained.dtype == dtype


def test_still_rejects_a_narrower_float():
    """Integers are exact; float32 is not. Admitting it would silently mix precisions through
    arithmetic that is float64 everywhere else, which is what the dtype gate is actually for."""
    kw = _valid_kwargs()
    kw["mass_rna_contained"] = np.array([1.0, 1.0], dtype=np.float32)
    with pytest.raises(ValueError, match="float64 or an integer count"):
        CalibrationResult(**kw)


def test_still_rejects_a_negative_integer_count():
    kw = _valid_kwargs()
    kw["mass_rna_contained"] = np.array([1, -1], dtype=np.int64)
    with pytest.raises(ValueError, match="non-negative"):
        CalibrationResult(**kw)


def test_rejects_length_mismatch():
    kw = _valid_kwargs()
    kw["mass_gdna_contained"] = np.zeros(3, dtype=np.float64)  # n_regions is 2
    with pytest.raises(ValueError):
        CalibrationResult(**kw)


@pytest.mark.parametrize(
    "field,value",
    [
        ("gdna_density_global", -1.0),
        ("gdna_density_global", np.inf),
        ("rna_sense_frac", 1.5),
        ("rna_sense_frac", -0.1),
    ],
)
def test_rejects_bad_scalars(field, value):
    kw = _valid_kwargs()
    kw[field] = value
    with pytest.raises(ValueError):
        CalibrationResult(**kw)
