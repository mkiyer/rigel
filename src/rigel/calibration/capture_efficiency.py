"""The capture efficiency of every region and boundary: the object's own gDNA density read against the bulk, under the
gDNA landscape as a two-group prior.

gDNA is one template at a uniform rate before capture, so an object's gDNA density after capture measures how strongly
capture enriched it. The bulk is the landscape's largest basin (``landscape.split_basins``), since gDNA off the probes
always outnumbers gDNA on them, and it is the unit: an efficiency is a relative yield, not a probability, and a factor
common to every efficiency cancels in the EM.

- Prior (Efron's two groups): the bulk basin collapses onto its peak ``x_b`` with the basin's whole mass ``π₀``, since
  off the probes there is one template and the basin's width is measurement; elsewhere the prior is the population's
  own rendered density (``DensityLandscape.logP_free``). Nothing is located, counted or gated, and a one-basin
  landscape reads exactly 1 everywhere.
- Likelihood: the object's deconvolved gDNA mass ``k`` on its gDNA support ``S``, the Poisson
  ``ℓ(x) = k·log(S·eˣ) − S·eˣ``.
- Readout: ``log e = E[min(x, s(x))] − x_b`` over the posterior, so the bulk reads 1. The mean is of the log density,
  so a lone prior level of mass ``μ`` at ``Δ`` nats above the bulk lifts an object with no evidence by ``e^{μΔ}``
  rather than a linear mean's ``1 + μ(e^Δ − 1)``. ``s(x)`` is the peak of the enriched basin holding ``x`` (``+∞``
  outside every enriched basin), so capture saturates per level.
- Saturation: ``exp(x_sat − x_b)`` for the most massive enriched basin's peak ``x_sat``, or 1 when there is none; the
  junction price saturates at it (``capture_eff_length``).

Each object reads the fragments the deposit rule gives it: a region its contained gDNA count on its contained support,
a boundary its crossing count on its crossing support, never apportioned onto the regions beside it. An object with
no gDNA support (``S = 0``) reads 1. Objects are read in row tiles (``simplex_logodds._block_rows``) with every sum
along one object's row, so the tile height changes no bit. Layer 5; ``calibrate`` publishes the reading on the result.

Gate: ``tests/calibration/test_capture_efficiency.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import CONSTANTS
from .landscape import DensityLandscape, census, split_basins
from .simplex_logodds import _block_rows

__all__ = ["CaptureReading", "capture_efficiencies"]

_TINY = np.finfo(np.float64).tiny


@dataclass(frozen=True, slots=True)
class CaptureReading:
    """The capture reader's answer for one library."""

    #: float64[n_regions] — each region's efficiency against the bulk, from its contained count.
    efficiency_region: np.ndarray
    #: float64[n_boundaries] — each boundary's efficiency against the bulk, from its crossing count.
    efficiency_boundary: np.ndarray
    #: The bulk basin's peak density, gDNA fragments per bp: the unit the efficiencies are read in.
    bulk_density: float
    #: The most massive enriched basin's peak over the bulk, ``> 1``; exactly 1 when nothing lies above the bulk.
    saturation: float


@dataclass(frozen=True, slots=True)
class _TwoGroups:
    """The two-group prior on the landscape's grid, and the level each grid cell reads at."""

    x: np.ndarray
    #: the bulk peak's grid index: the census reads its peaks off this grid
    ib: int
    #: ``log π₀``, the bulk basin's mass; ``-inf`` when the population puts none there
    log_pi0: float
    #: the alternative's log mass per cell: ``-inf`` inside the bulk basin and wherever the population is empty
    log_alt: np.ndarray
    #: ``min(x, s(x))``: the level each cell reads at, saturated at its enriched basin's peak
    level: np.ndarray
    saturation: float


def _two_groups(landscape: DensityLandscape) -> _TwoGroups:
    x = np.asarray(landscape.log_rho, dtype=np.float64)
    free = np.asarray(landscape.logP_free, dtype=np.float64)
    if x.ndim != 1 or x.size < 2 or not np.all(np.isfinite(x)) or free.shape != x.shape:
        raise ValueError(
            "capture_efficiencies: the landscape needs a finite grid of two cells or more, with its own "
            "density on the same grid."
        )
    if np.any(np.isnan(free)) or np.any(free == np.inf) or not np.any(np.isfinite(free)):
        raise ValueError(
            "capture_efficiencies: the landscape's own density is empty or not finite."
        )
    bulk, above = split_basins(census(landscape))
    ib = int(np.argmin(np.abs(x - bulk.log_rho)))
    in_bulk = (x >= bulk.lo) & (x <= bulk.hi)
    p = np.exp(free - np.max(free))
    p = p / p.sum()
    pi0 = float(p[in_bulk].sum())
    alt = np.where(in_bulk, 0.0, p)
    log_alt = np.where(alt > 0.0, np.log(np.maximum(alt, _TINY)), -np.inf)
    cap = np.full(x.shape, np.inf)
    for m in above:
        cap[(x >= m.lo) & (x <= m.hi)] = m.log_rho
    saturation = (
        float(np.exp(max(above, key=lambda m: m.basin_mass).log_rho - x[ib])) if above else 1.0
    )
    return _TwoGroups(
        x=x,
        ib=ib,
        log_pi0=float(np.log(pi0)) if pi0 > 0.0 else -np.inf,
        log_alt=log_alt,
        level=np.minimum(x, cap),
        saturation=saturation,
    )


def _read(prior: _TwoGroups, k: np.ndarray, S: np.ndarray) -> np.ndarray:
    """``exp(E[min(x, s(x))] − x_b)`` per object: the two-group posterior of its own count ``k`` on its support
    ``S``, the bulk as the one level ``x_b`` with its basin's mass. ``S = 0``: 1."""
    x = prior.x
    xb = float(x[prior.ib])
    rate = np.exp(x)
    n = k.shape[0]
    out = np.ones(n)
    rows = _block_rows(x.size, CONSTANTS.resources.capture_tile_cell_bytes)
    for r0 in range(0, n, rows):
        sl = slice(r0, min(r0 + rows, n))
        lam = rate[None, :] * S[sl, None]
        ll = np.where(
            S[sl, None] > 0.0,
            k[sl, None] * np.log(np.maximum(lam, CONSTANTS.calibration.rate_log_floor)) - lam,
            0.0,
        )
        ll -= ll.max(1, keepdims=True)
        a_null = prior.log_pi0 + ll[:, prior.ib]
        a_alt = prior.log_alt[None, :] + ll
        top = np.maximum(a_null, a_alt.max(1))
        top = np.where(np.isfinite(top), top, 0.0)
        w_null = np.exp(a_null - top)
        w_alt = np.exp(a_alt - top[:, None])
        den = np.maximum(w_null + w_alt.sum(1), _TINY)
        # the level moment as a sum along the object's own row, not a matrix–vector product: BLAS need not sum a row
        # in one order at every tile height, and the height is a resource setting
        mean_level = (w_null * xb + (w_alt * prior.level).sum(1)) / den
        out[sl] = np.where(S[sl] > 0.0, np.exp(mean_level - xb), 1.0)
    return out


def capture_efficiencies(
    landscape: DensityLandscape,
    count_region: np.ndarray,
    eff_region: np.ndarray,
    count_boundary: np.ndarray,
    eff_boundary: np.ndarray,
) -> CaptureReading:
    """Every region's efficiency from its contained count on its contained support, every boundary's from its
    crossing count on its crossing support, with the bulk density and the saturation level they are read against.

    The counts are the deconvolved gDNA masses, the supports the objects' own gDNA opportunities. Both must be
    finite; a negative value (a deconvolution's rounding) reads as 0.
    """
    prior = _two_groups(landscape)

    def clean(values, name: str) -> np.ndarray:
        v = np.asarray(values, dtype=np.float64)
        if v.ndim != 1 or not np.all(np.isfinite(v)):
            raise ValueError(f"capture_efficiencies: {name} must be a finite 1-D array.")
        return np.maximum(v, 0.0)

    k_r, s_r = clean(count_region, "count_region"), clean(eff_region, "eff_region")
    k_b, s_b = clean(count_boundary, "count_boundary"), clean(eff_boundary, "eff_boundary")
    if k_r.shape != s_r.shape or k_b.shape != s_b.shape:
        raise ValueError(
            "capture_efficiencies: every count needs its own support, on the same axis."
        )
    return CaptureReading(
        efficiency_region=_read(prior, k_r, s_r),
        efficiency_boundary=_read(prior, k_b, s_b),
        bulk_density=float(np.exp(prior.x[prior.ib])),
        saturation=prior.saturation,
    )
