"""rigel.calibration.density_loco — Locoregional NB EB shrinkage.

A single closed-form formula: shrink a per-locus NB density estimate
toward the corresponding global density using the type-specific NB
overdispersion :math:`\\kappa`.  Per
``docs/calibration/calibration_v6_plan.md`` §2.6.1:

.. math::

    \\hat\\rho_{\\mathrm{loco}} =
        \\frac{N_{\\mathrm{loco}} + \\kappa \\cdot \\hat\\rho_{\\mathrm{global}}}
             {L_{\\mathrm{eff,loco}} + \\kappa}

Limits:

* :math:`\\kappa \\to 0` ⇒ ``N / L_eff`` (pure local).
* :math:`\\kappa \\to \\infty` ⇒ ``ρ_global`` (pure global).
* ``L_eff == 0`` ⇒ ``ρ_global`` (degenerate locus, no support).

Single source of truth for the locoregional shrinkage; never inline
this formula at call sites.
"""

from __future__ import annotations


__all__ = ["shrink_to_loco"]


def shrink_to_loco(
    n_loco: float,
    leff_loco: float,
    rho_global: float,
    kappa: float,
) -> float:
    """Closed-form NB EB shrinkage.

    Parameters
    ----------
    n_loco
        Per-locus count of fragments contributing to this density type.
    leff_loco
        Per-locus effective length (bp) for this density type.
    rho_global
        Global density (fragments / bp) for this type.
    kappa
        NB overdispersion for this density type.

    Returns
    -------
    float
        Shrunk per-locus density.  Always finite if inputs are finite.
    """
    if not (n_loco >= 0.0 and leff_loco >= 0.0 and rho_global >= 0.0 and kappa >= 0.0):
        raise ValueError(
            f"shrink_to_loco: all inputs must be finite and non-negative; got "
            f"n_loco={n_loco!r}, leff_loco={leff_loco!r}, "
            f"rho_global={rho_global!r}, kappa={kappa!r}."
        )
    denom = leff_loco + kappa
    if denom <= 0.0:
        return rho_global
    return (n_loco + kappa * rho_global) / denom
