"""Per-region gDNA exposure — closed-form Gamma posterior (doc 03 §4, G4).

After the three E-step views are allocated, the per-region total gDNA mass is the
**D1 side-attributed** sum — contained plus each boundary side, each already
attributed to its single region by PR 2.5's per-side flux/mass. There is **no ½
half-split** (the doc 01 §7 / doc 03 §4 ½ is superseded by D1; PR04a §I.3):

    M_g_tot = M_g_contained + M_g_left + M_g_right

The Gamma-Poisson conjugacy (doc 01 §4.2) then gives a closed-form posterior on
the per-region exposure ω (multiplicative on ρ_0)::

    α_post = 1/φ + M_g_tot
    β_post = 1/φ + ρ_0 · L_eff
    ω             = α_post / β_post          (posterior mean)
    Var(log ω)    = 1 / α_post               (delta method)

Closed form, O(R), no iteration. Each boundary side carries its own region's
exposure: the left/right sides of a boundary belong to different regions and
fold into those regions' ω (III.4).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class Exposure:
    """Per-region exposure posterior + the gDNA mass it was computed from."""

    omega: np.ndarray  # float64[R] — posterior mean E[ω | data], > 0
    log_omega_var: np.ndarray  # float64[R] — delta-method Var(log ω), > 0
    m_g_tot: np.ndarray  # float64[R] — D1-aggregated total gDNA mass


def exposure_posterior(
    m_g_contained: np.ndarray,
    m_g_left: np.ndarray,
    m_g_right: np.ndarray,
    *,
    rho_0: float,
    L_eff: np.ndarray,
    exposure_dispersion: float,
) -> Exposure:
    """Closed-form Gamma exposure posterior with D1 side-attribution (no ½).

    ``exposure_dispersion`` is the variance of the exposure prior
    ``ω ~ Gamma(s, s)`` with shape/rate ``s = 1/exposure_dispersion`` (was ``φ``).
    """
    inv_disp = 1.0 / exposure_dispersion
    m_g_tot = m_g_contained + m_g_left + m_g_right
    alpha_post = inv_disp + m_g_tot
    beta_post = inv_disp + rho_0 * np.asarray(L_eff, dtype=np.float64)
    return Exposure(
        omega=alpha_post / beta_post,
        log_omega_var=1.0 / alpha_post,
        m_g_tot=m_g_tot,
    )


__all__ = ["Exposure", "exposure_posterior"]
