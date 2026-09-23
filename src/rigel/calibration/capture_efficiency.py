"""The capture efficiency of every region and boundary — the posterior mean of the object's clipped gDNA
density against the fully captured level, under the fitted gDNA landscape, from the object's OWN count.

gDNA is one template at a uniform rate before capture, so its density after capture at an object is that
object's capture efficiency up to the one unit ``ρ_ref`` — the located enriched mode of the landscape
(`abundance_landscape.located_enriched_mode`), the fully captured level. An object's efficiency is

    c̃_o = E[ min(ρ_o / ρ_ref, 1) | k_o, S_o ],

the expectation under the landscape ``P(log ρ)`` as the prior (the population's own statement of where
gDNA densities sit) and the Poisson counting rule as the likelihood. At high depth it is the plug-in
``min(k/S/ρ_ref, 1)``; at low depth the population's mixture weighted by the object's own likelihood; an
object with no evidence reads the population's clipped mean. No constant enters and no floor: the
multimapper floor ``C/(C+1)`` this replaced was a 30–40× bias on every unprobed transcript
(`ISSUES: ruler-multimapper-floor-caps-the-correction`), protecting against a plug-in's exact zero that a
posterior mean never produces.

EACH OBJECT READS ITS OWN EVIDENCE — the fragments the deposit rule gives it. A region holds the fragments
contained in it: its gDNA count ``k_r`` on its contained support ``S_r``. A boundary holds the fragments
that cross it: its gDNA count ``k_e`` on its crossing support ``S_e``. That is the split the
capture-contracted length prices (`capture_eff_length`, `priors`): a component's contained share of a
region at the region's efficiency and its conserved share at a boundary at the boundary's, so the
crossings are read once, by the boundary that holds them, and never apportioned onto the regions beside it
as well. A region too short to contain a fragment has no contained share to price and no count to read;
its efficiency is the population's, and no length multiplies it.

Layer 5: a posterior under the density prior. The consumer is `calibrate`, which publishes the
efficiencies on the result for `capture_eff_length` and `priors`. Gate:
``tests/calibration/test_capture_efficiency.py``.
"""

from __future__ import annotations

import numpy as np

from .landscape import DensityLandscape
from .simplex_logodds import _block_rows

__all__ = ["capture_efficiencies"]

_EPS = 1e-300


def _posterior(
    landscape: DensityLandscape, rho_ref: float, k: np.ndarray, S: np.ndarray
) -> np.ndarray:
    """``E[min(ρ/ρ_ref, 1)]`` per object on the landscape's grid: the prior times the Poisson likelihood of
    the object's count ``k`` on its support ``S`` (``S = 0``: no likelihood, the prior alone). Row tiles, as
    ψ tiles (`simplex_logodds._block_rows`): a million objects never exist against the grid at once."""
    rho = np.exp(np.asarray(landscape.log_rho, dtype=np.float64))
    lp = np.asarray(landscape.logP, dtype=np.float64)
    clipped = np.minimum(rho / rho_ref, 1.0)
    n = k.shape[0]
    out = np.empty(n)
    rows = _block_rows(rho.size, 8)
    for r0 in range(0, n, rows):
        sl = slice(r0, min(r0 + rows, n))
        lam = rho[None, :] * S[sl, None]
        ll = np.where(S[sl, None] > 0.0, k[sl, None] * np.log(np.maximum(lam, _EPS)) - lam, 0.0)
        lw = ll + lp[None, :]
        lw -= lw.max(1, keepdims=True)
        w = np.exp(lw)
        w /= w.sum(1, keepdims=True)
        out[sl] = w @ clipped
    return out


def capture_efficiencies(
    landscape: DensityLandscape,
    rho_ref: float,
    count_region: np.ndarray,
    eff_region: np.ndarray,
    count_boundary: np.ndarray,
    eff_boundary: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """``(efficiency_region, efficiency_boundary)``: every region's ``c̃_r`` from its contained count on its
    contained support, every boundary's ``c̃_e`` from its crossing count on its crossing support.

    The counts are the deconvolved gDNA masses (a mass, so the Poisson enters through the gamma function as
    a continuation), the supports the objects' own gDNA opportunities.
    """

    def clean(values):
        return np.maximum(np.asarray(values, dtype=np.float64), 0.0)

    return (
        _posterior(landscape, rho_ref, clean(count_region), clean(eff_region)),
        _posterior(landscape, rho_ref, clean(count_boundary), clean(eff_boundary)),
    )
