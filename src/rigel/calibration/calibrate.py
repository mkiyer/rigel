"""calibrate() — the joint fractional-accumulator + calibration-v6 orchestrator.

Consumes the accumulator payload, the region geometry, and the trained strand
model and returns a :class:`~rigel.calibration.result.CalibrationResult` with
the recovered library hyperparameters and per-region deconvoluted mass.

**PR 4a — single E-step pass.** The body seeds the global gDNA density ``ρ_0``
(:mod:`rigel.calibration.density`), fits the RNA strand-balance model
(:mod:`rigel.calibration.strand_balance`, PR 3), then runs **one** E-step
(:mod:`rigel.calibration.estep`) over the three substrate views and the
closed-form exposure posterior (:mod:`rigel.calibration.exposure`). All library
hyperparameters are **fixed** at their initial values; the M-step and the outer
loop that iterates to convergence are PR 5.

The single pass uses ``ω = 1``, ``π_g_prior = 0.5``, and ``μ_d = 0`` (no prior
RNA mass — doc 03 §3.1), so the count log-Bayes-factor is identically 0 and the
allocation is strand-driven + spliced-deterministic. The count channel is fully
built (and unit-tested with a non-zero ``μ_d``); it engages once PR 5's outer
loop feeds ``M_d`` back. Hence ``n_iterations = 1`` and ``converged = False``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .density import estimate_gdna_density
from .estep import estep_view
from .exposure import exposure_posterior
from .result import CalibrationResult
from .strand_balance import fit_strand_balance
from .substrate import CalibrationSubstrate

if TYPE_CHECKING:
    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays

# Doc 03 §7 initialization values (iteration-0). PR 5's M-step owns these
# definitively; here they seed the single pass / the outer-loop warm start.
# (ρ_0 is owned by density.py; ρ_r_bb / κ_rna by the strand-balance model.)
_RHO_D_BB_INIT = 0.01
_EPS_S_INIT = 1.0e-3

# Uninformative gDNA-mixing prior for the single pass (doc 03 §7). PR 5's outer
# loop replaces it with the data-driven §5.6 update.
_PI_G_PRIOR_INIT = 0.5


def initial_hyperparameters(
    substrate: CalibrationSubstrate,
    config: "CalibrationConfig",
) -> tuple[float, float, float]:
    """Doc 03 §7 data-driven initial hyperparameters ``(φ, ρ_d_bb, ε_s)``.

    ``φ`` is the NB moment estimator, floored to stay valid with little/no data.
    Reused as the PR 5 outer-loop warm start. ``ρ_0`` comes from
    :func:`rigel.calibration.density.estimate_gdna_density`; ``κ_rna`` / ``ρ_r_bb``
    from the strand-balance model — neither is computed here.
    """
    n_u = substrate.contained.n_unspliced.astype(np.float64)
    total_n = float(n_u.sum())

    if total_n > 0.0 and n_u.size > 1:
        mean = float(n_u.mean())
        phi = float(n_u.var()) / mean**2 - 1.0 / mean if mean > 0.0 else config.phi_floor
    else:
        phi = config.phi_floor
    phi = max(phi, config.phi_floor)

    return phi, _RHO_D_BB_INIT, _EPS_S_INIT


def calibrate(
    *,
    payload: "AccumulatorPayload",
    region_arrays: "RegionArrays",
    strand_model: "StrandModels",
    config: "CalibrationConfig",
) -> CalibrationResult:
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    r = substrate.n_regions

    # --- Library hyperparameters (fixed this PR) ---
    density = estimate_gdna_density(substrate, region_arrays)  # ρ_0 seed (§I.1 / III.2)
    rho_0 = density.rho_0
    phi, rho_d_bb, eps_s = initial_hyperparameters(substrate, config)
    strand = fit_strand_balance(substrate, strand_model)  # κ_rna, ρ_r_bb (PR 3)

    # --- Single E-step pass over the three views (ω ≡ 1, μ_d ≡ 0) ---
    L_eff = substrate.L_eff
    ts_class = substrate.ts_class
    pass_kwargs = dict(
        omega=np.ones(r, dtype=np.float64),
        rho_0=rho_0,
        L_eff=L_eff,
        phi=phi,
        kappa_rna=strand.kappa_rna,
        rho_r_bb=strand.rho_r_bb,
        rho_d_bb=rho_d_bb,
        eps_s=eps_s,
        pi_g_prior=np.full(r, _PI_G_PRIOR_INIT, dtype=np.float64),
        m_d_unspl_prev=np.zeros(r, dtype=np.float64),
    )
    contained = estep_view(substrate.contained, ts_class, **pass_kwargs)
    left = estep_view(substrate.left, ts_class, **pass_kwargs)
    right = estep_view(substrate.right, ts_class, **pass_kwargs)

    # --- Exposure posterior (G4): D1 aggregation, no ½ ---
    exposure = exposure_posterior(
        contained.m_g, left.m_g, right.m_g, rho_0=rho_0, L_eff=L_eff, phi=phi
    )

    # One pass: the mass-change diagnostic is the change from the zero-gDNA
    # initialization, max_r |M_g_tot − 0| / (0 + 1) = max_r M_g_tot.
    mass_change = float(exposure.m_g_tot.max()) if r > 0 else 0.0

    return CalibrationResult(
        mass_g_contained=contained.m_g,
        mass_d_contained=contained.m_d,
        mass_g_left=left.m_g,
        mass_d_left=left.m_d,
        mass_g_right=right.m_g,
        mass_d_right=right.m_d,
        omega=exposure.omega,
        log_omega_var=exposure.log_omega_var,
        rho_0=rho_0,
        phi=phi,
        rho_d_bb=rho_d_bb,
        kappa_rna=strand.kappa_rna,
        rho_r_bb=strand.rho_r_bb,
        eps_s=eps_s,
        n_iterations=1,
        converged=False,
        mass_change_history=np.array([mass_change], dtype=np.float64),
        n_regions=r,
        config=config,
    )


__all__ = ["calibrate", "initial_hyperparameters"]
