"""calibrate() — the joint fractional-accumulator + calibration-v6 orchestrator.

Consumes the accumulator payload, the region geometry, and the trained strand
model and returns a :class:`~rigel.calibration.result.CalibrationResult` with
the recovered library hyperparameters and per-region deconvoluted mass.

**Placeholder during PR 2.** The body builds the real
:class:`~rigel.calibration.substrate.CalibrationSubstrate` and returns the
trivial "iteration-0" result — **no gDNA inferred**: all contained/boundary
mass attributed to RNA, unit exposure ``ω = 1`` (uniform; there is no data to
train an exposure on), hyperparameters at their doc 03 §7 initial values. The
real deconvolution (E-step/exposure → M-step/outer loop) replaces this body in
PR 2–PR 5; see ``docs/acc_caljointmodel/00_implementation_plan.md``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .result import CalibrationResult
from .substrate import CalibrationSubstrate

if TYPE_CHECKING:
    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays

# Doc 03 §7 initialization table (iteration-0 dispersions). PR 5's M-step owns
# these definitively; here they seed the placeholder / the outer-loop warm start.
_RHO_D_BB_INIT = 0.01
_RHO_R_BB_INIT = 0.01
_EPS_S_INIT = 1.0e-3


def initial_hyperparameters(
    substrate: CalibrationSubstrate,
    config: "CalibrationConfig",
) -> tuple[float, float, float, float, float]:
    """Doc 03 §7 data-driven initial hyperparameters (rho_0, phi, ρ_d, ρ_r, ε_s).

    ``rho_0`` is the half-and-half density assumption; ``phi`` is the NB moment
    estimator. Both are floored to stay valid when there is little/no data.
    Reused as the PR 5 outer-loop warm start.
    """
    n_u = substrate.contained.n_unspliced.astype(np.float64)
    total_L = float(substrate.L_eff.sum())
    total_n = float(n_u.sum())

    if total_L > 0.0:
        # Half of the unspliced fragments assumed gDNA; floored at "one
        # fragment over the genome" so rho_0 stays positive with no data.
        rho_0 = max(0.5 * total_n / total_L, 1.0 / total_L)
    else:
        rho_0 = 1.0

    if total_n > 0.0 and n_u.size > 1:
        mean = float(n_u.mean())
        phi = float(n_u.var()) / mean**2 - 1.0 / mean if mean > 0.0 else config.phi_floor
    else:
        phi = config.phi_floor
    phi = max(phi, config.phi_floor)

    return rho_0, phi, _RHO_D_BB_INIT, _RHO_R_BB_INIT, _EPS_S_INIT


def calibrate(
    *,
    payload: "AccumulatorPayload",
    region_arrays: "RegionArrays",
    strand_model: "StrandModels",
    config: "CalibrationConfig",
) -> CalibrationResult:
    # strand_model is unused by the placeholder; PR 3/PR 4 consume it for the
    # strand-balance model and the E-step. It stays in the signature so the
    # live call site is exercised end-to-end.
    del strand_model

    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    r = substrate.n_regions
    rho_0, phi, rho_d_bb, rho_r_bb, eps_s = initial_hyperparameters(substrate, config)

    # No gDNA inferred: all mass → RNA, exposure uniform at the prior mean 1.0.
    # log_omega_var = 1/(1/phi + M_g_tot) = phi since M_g_tot ≡ 0 (doc 04 §6.4).
    return CalibrationResult(
        mass_g_contained=np.zeros(r, dtype=np.float64),
        mass_d_contained=substrate.contained.mass_unspliced + substrate.contained.mass_spliced,
        mass_g_left=np.zeros(r, dtype=np.float64),
        mass_d_left=substrate.left.mass_unspliced + substrate.left.mass_spliced,
        mass_g_right=np.zeros(r, dtype=np.float64),
        mass_d_right=substrate.right.mass_unspliced + substrate.right.mass_spliced,
        omega=np.ones(r, dtype=np.float64),
        log_omega_var=np.full(r, phi, dtype=np.float64),
        rho_0=rho_0,
        phi=phi,
        rho_d_bb=rho_d_bb,
        rho_r_bb=rho_r_bb,
        eps_s=eps_s,
        n_iterations=0,
        converged=True,
        mass_change_history=np.empty(0, dtype=np.float64),
        n_regions=r,
        config=config,
    )


__all__ = ["calibrate", "initial_hyperparameters"]
