"""The capture efficiency of every region and boundary — the object's own gDNA density read against the bulk, with
the fitted gDNA landscape as a two-group prior and the object's OWN count as the likelihood.

gDNA is one template at a uniform rate before capture, so after capture an object's gDNA density measures how
strongly capture enriched it. Off the probes the template is not enriched at all: that is the BULK. The landscape
``P(log ρ)`` is the population's own statement of where gDNA densities sit, and its largest basin
(`landscape.census`, `landscape.split_basins`) is the bulk, since gDNA off the probes always outnumbers gDNA on them.

THE PRIOR — Efron's two groups. The bulk basin collapses onto its peak ``x_b`` and carries the basin's whole mass
``π₀``: off the probes there is one uniform template, and the basin's width is measurement, not spectrum. Everywhere
else the prior is the population's own rendered density as it is (``DensityLandscape.logP_free``): the enriched
spectrum, at whatever levels and in whatever quantity it sits. No mode is located, no member is counted, no width is
tested and nothing is gated. A library whose landscape is one basin reads exactly 1 everywhere, because the prior
is then a point.

THE LIKELIHOOD. The object's own count: its deconvolved gDNA mass ``k`` on its gDNA support ``S``, as the Poisson
``ℓ(x) = k·log(S·eˣ) − S·eˣ`` (a mass, so the Poisson continues through the gamma function).

THE READOUT. ``log e = E[min(x, s(x))] − x_b`` over the posterior, so the bulk reads exactly 1:

- the posterior mean of the LOG density. An object with no evidence of its own reads the prior, and a lone level of
  prior mass ``μ`` at ``Δ`` nats above the bulk lifts it by ``e^{μΔ}``, where a linear mean would lift it by
  ``1 + μ(e^Δ − 1)``;
- CAPTURE SATURATES PER LEVEL: ``s(x)`` is the peak of the enriched basin holding ``x`` (``+∞`` outside every
  enriched basin), so a density above its own basin's peak reads that peak.

THE SATURATION LEVEL. ``saturation = exp(x_sat − x_b)``, with ``x_sat`` the peak of the most massive enriched basin,
and 1 when there is none: the fully captured level against the bulk. The junction price saturates at it
(`capture_eff_length`).

THE UNIT IS THE BULK. An efficiency is a relative yield, not a probability: an enriched object reads far above 1.
Every consumer prices every component's opportunity at these efficiencies, so a factor common to all of them
cancels in the EM (`capture_eff_length`, `priors`).

EACH OBJECT READS ITS OWN EVIDENCE — the fragments the deposit rule gives it. A region holds the fragments contained
in it: its gDNA count on its contained support. A boundary holds the fragments that cross it: its gDNA count on its
crossing support. That is the split the capture-contracted length prices (`capture_eff_length`, `priors`): a
component's contained share of a region at the region's efficiency, and its conserved share at a boundary at the
boundary's. The crossings are read once, by the boundary that holds them, and never apportioned onto the regions
beside it as well. An object with no gDNA support (``S = 0``) has no likelihood and reads 1, the null: no gDNA
fragment could sit there to witness its capture.

The objects are read in row tiles (`simplex_logodds._block_rows`), so a million objects never exist against the
grid at once; every sum runs along one object's own row, so the tile's height moves no bit.

Layer 5: a posterior under the density prior. The consumer is `calibrate`, which publishes the efficiencies, the
bulk density and the saturation level on the result. Gate: ``tests/calibration/test_capture_efficiency.py``.
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
        # the level moment as a sum over the object's own row, not a matrix–vector product: BLAS need not sum a
        # row in one order at every tile height (a one-row tile takes another kernel), and the height is a
        # resource setting (``tests/calibration/test_capture_efficiency.py``)
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
