"""rigel._accumulator — Python wrapper around the native fractional Accumulator.

The native class (``rigel._bam_impl.Accumulator``) exposes the raw
storage as 2-D numpy views: ``regions_contained`` (N,4 uint32),
``boundaries_mass_left/right`` (N+1, 4 float32),
``boundaries_flux_left/right`` (N+1, 4 uint32 — per-side flux).

The spec tests in ``tests/native/test_accumulator_spec.py`` expect
dataclass-style indexed access (``acc.regions[i].contained[ch]``,
``acc.boundaries[i].mass_left[ch]``), mirroring the Python reference
in ``tests/native/_accumulator_reference.py``. This module supplies
that thin façade.

The wrapper holds a reference to the native instance for the lifetime
of the row views (the numpy arrays are zero-copy views over native
storage).
"""

from __future__ import annotations

import numpy as np

from ._bam_impl import Accumulator as _NativeAccumulator


class _RegionRow:
    """View over one region's channel counts (uint32[4])."""

    __slots__ = ("contained",)

    def __init__(self, contained_row: np.ndarray) -> None:
        self.contained = contained_row


class _BoundaryRow:
    """View over one boundary's per-channel mass and per-side flux."""

    __slots__ = ("mass_left", "mass_right", "flux_left", "flux_right")

    def __init__(
        self,
        mass_left_row: np.ndarray,
        mass_right_row: np.ndarray,
        flux_left_row: np.ndarray,
        flux_right_row: np.ndarray,
    ) -> None:
        self.mass_left = mass_left_row
        self.mass_right = mass_right_row
        self.flux_left = flux_left_row
        self.flux_right = flux_right_row


class _RegionList:
    """Sequence-like view over Accumulator's regions."""

    __slots__ = ("_owner", "_contained")

    def __init__(self, owner: "Accumulator") -> None:
        self._owner = owner
        self._contained = owner._native.regions_contained

    def __len__(self) -> int:
        return self._contained.shape[0]

    def __getitem__(self, i: int) -> _RegionRow:
        return _RegionRow(self._contained[i])

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


class _BoundaryList:
    """Sequence-like view over Accumulator's boundaries."""

    __slots__ = ("_owner", "_ml", "_mr", "_fl_l", "_fl_r")

    def __init__(self, owner: "Accumulator") -> None:
        self._owner = owner
        self._ml = owner._native.boundaries_mass_left
        self._mr = owner._native.boundaries_mass_right
        self._fl_l = owner._native.boundaries_flux_left
        self._fl_r = owner._native.boundaries_flux_right

    def __len__(self) -> int:
        return self._ml.shape[0]

    def __getitem__(self, i: int) -> _BoundaryRow:
        return _BoundaryRow(self._ml[i], self._mr[i], self._fl_l[i], self._fl_r[i])

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


class Accumulator:
    """Fractional accumulator (one reference).

    Mirrors the API of ``tests.native._accumulator_reference.Accumulator``
    over the native C++ implementation in ``rigel._bam_impl``.

    Construct with the sorted boundary-position array (int64, length N+1):

        acc = Accumulator(
            boundary_positions=np.array([0, 100, 200, 400], dtype=np.int64),
        )

    Then deposit fragments:

        acc.deposit(blocks=[(50, 100), (100, 180)], spliced=False, primary=True)
    """

    __slots__ = ("_native", "regions", "boundaries")

    def __init__(self, boundary_positions: np.ndarray) -> None:
        positions = np.ascontiguousarray(boundary_positions, dtype=np.int64)
        self._native = _NativeAccumulator(positions)
        self.regions = _RegionList(self)
        self.boundaries = _BoundaryList(self)

    @property
    def n_regions(self) -> int:
        return int(self._native.n_regions)

    @property
    def n_boundaries(self) -> int:
        return int(self._native.n_boundaries)

    @property
    def boundary_positions(self) -> np.ndarray:
        return np.asarray(self._native.boundary_positions)

    def region_of_pos(self, pos: int) -> int:
        return int(self._native.region_of_pos(int(pos)))

    def left_boundary_of(self, region_idx: int) -> int:
        return int(region_idx)

    def right_boundary_of(self, region_idx: int) -> int:
        return int(region_idx) + 1

    def deposit(
        self,
        blocks: list[tuple[int, int]],
        spliced: bool,
        primary: bool,
    ) -> None:
        if not blocks:
            return
        starts = np.empty(len(blocks), dtype=np.int64)
        ends = np.empty(len(blocks), dtype=np.int64)
        for i, (s, e) in enumerate(blocks):
            starts[i] = s
            ends[i] = e
        self._native.deposit(starts, ends, bool(spliced), bool(primary))

    def total_mass_deposited(self) -> float:
        total = 0.0
        total += float(self._native.regions_contained.sum())
        total += float(self._native.boundaries_mass_left.sum())
        total += float(self._native.boundaries_mass_right.sum())
        return total


__all__ = ["Accumulator"]
