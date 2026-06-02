"""calibrate() — the joint fractional-accumulator + calibration-v6 orchestrator.

Consumes the accumulator payload, the region geometry, and the trained strand
model and returns a :class:`~rigel.calibration.result.CalibrationResult` with
the recovered library hyperparameters and per-region deconvoluted mass.

**PR 5 — the working calibrator (M-step + outer loop, doc 03 §1).** Seeds ``ρ_0``
(:mod:`rigel.calibration.density`) and the RNA strand model (κ_rna, ρ_r_bb fixed
from the spliced channel — :mod:`rigel.calibration.strand_balance`, PR 3), then
iterates: **E-step** (:mod:`rigel.calibration.estep`) over the three views — the
**count channel is live** (μ_d from the previous iteration; zero on iter 1, so
iteration 1 reproduces the PR 4a single pass) and uses the **FL-corrected
per-view gDNA effective length** (contained ``E_f[max(0,L−ℓ)]``, boundary
``μ_FL``; :mod:`rigel.calibration.effective_length`, PR 4c) — then the **exposure
posterior** (G4, physical length), the **AMBIG sweep**
(:mod:`rigel.calibration.sweep`, D7 leg 2), the **M-step**
(:mod:`rigel.calibration.mstep` — fits ρ_0, ε_s, φ, ρ_d_bb), and the π_g-prior
update. The loop runs until the max-region relative mass change drops below
``mass_rel_tol`` (``converged=True``) or ``max_outer_iterations`` is hit. The
exposure ``ω`` keeps physical length (4a's total-mass form is exact); the
FL-corrected lengths enter only the per-view count channel and the sweep.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from .density import estimate_gdna_density
from .effective_length import boundary_eff_length, region_eff_length
from .estep import estep_view
from .exposure import exposure_posterior
from .mstep import (
    decodable_node_masks,
    fit_rho_d_bb,
    update_exposure_dispersion,
    update_pi_g_prior,
    update_rho_0,
)
from .result import CalibrationResult
from .strand_balance import fit_strand_balance
from .substrate import CalibrationSubstrate
from .sweep import sweep_ambig_exposure

if TYPE_CHECKING:
    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays

logger = logging.getLogger(__name__)

# Doc 03 §7 initialization values (iteration-0). PR 5's M-step owns these
# definitively; here they seed the single pass / the outer-loop warm start.
# (ρ_0 is owned by density.py; ρ_r_bb / κ_rna by the strand-balance model.)
_RHO_D_BB_INIT = 0.01

# Uninformative gDNA-mixing prior for the single pass (doc 03 §7). PR 5's outer
# loop replaces it with the data-driven §5.6 update.
_PI_G_PRIOR_INIT = 0.5


def initial_hyperparameters(
    substrate: CalibrationSubstrate,
    config: "CalibrationConfig",
) -> tuple[float, float]:
    """Doc 03 §7 data-driven initial hyperparameters ``(exposure_dispersion, ρ_d_bb)``.

    The ``exposure_dispersion`` seed is the NB moment estimator on the contained
    unspliced counts, floored to stay valid with little/no data — just a warm
    start; the outer-loop EM M-step (:func:`update_exposure_dispersion`) takes
    over. ``ρ_0`` comes from :func:`rigel.calibration.density.estimate_gdna_density`;
    ``κ_rna`` / ``ρ_r_bb`` from the strand-balance model — neither is computed here.
    """
    n_u = substrate.contained.n_unspliced.astype(np.float64)
    total_n = float(n_u.sum())
    floor = config.exposure_dispersion_floor

    if total_n > 0.0 and n_u.size > 1:
        mean = float(n_u.mean())
        exposure_dispersion = float(n_u.var()) / mean**2 - 1.0 / mean if mean > 0.0 else floor
    else:
        exposure_dispersion = floor
    exposure_dispersion = max(exposure_dispersion, floor)

    return exposure_dispersion, _RHO_D_BB_INIT


def calibrate(
    *,
    payload: "AccumulatorPayload",
    region_arrays: "RegionArrays",
    strand_model: "StrandModels",
    gdna_fl_pmf: np.ndarray,
    config: "CalibrationConfig",
) -> CalibrationResult:
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    r = substrate.n_regions
    ts_class = substrate.ts_class
    l_phys = substrate.L_eff  # exposure ω uses physical length (PR 4a, exact)

    # --- Initial hyperparameters ---
    rho_0 = estimate_gdna_density(substrate, region_arrays).rho_0  # ρ_0 seed (PR 4c §III.2)
    exposure_dispersion, rho_d_bb = initial_hyperparameters(substrate, config)
    strand = fit_strand_balance(substrate, strand_model)  # κ_rna, ρ_r_bb fixed (PR 3 / III.1)
    kappa_rna, rho_r_bb = strand.kappa_rna, strand.rho_r_bb

    # FL-corrected per-view gDNA exposures for the count channel (PR 4c / PR05):
    # contained = E_f[max(0, L−ℓ)], boundary = μ_FL. (ω above keeps physical L.)
    region_eff_len = region_eff_length(l_phys, gdna_fl_pmf)
    mu_fl = boundary_eff_length(gdna_fl_pmf)
    boundary_eff = np.full(r, mu_fl, dtype=np.float64)

    # Per-node decodability (PR 8): ρ_0 + exposure dispersion are fit from decodable
    # nodes only; AMBIG regions are imputed, never used to set the global density.
    dec_region, dec_left_bnd, dec_right_bnd = decodable_node_masks(
        ts_class, region_arrays.ref_id
    )

    # --- Outer EM loop (doc 03 §1) ---
    omega = np.ones(r, dtype=np.float64)
    log_omega_var = np.full(r, exposure_dispersion, dtype=np.float64)
    pi_g_prior = np.full(r, _PI_G_PRIOR_INIT, dtype=np.float64)
    m_d_cont = np.zeros(r, dtype=np.float64)
    m_d_left = np.zeros(r, dtype=np.float64)
    m_d_right = np.zeros(r, dtype=np.float64)
    n_u_cont = substrate.contained.n_unspliced
    m_g_tot_prev: "np.ndarray | None" = None
    mass_changes: list[float] = []
    converged = False
    n_iter = 0
    contained = left = right = None

    for it in range(1, config.max_outer_iterations + 1):
        n_iter = it
        shared = dict(
            omega=omega,
            rho_0=rho_0,
            exposure_dispersion=exposure_dispersion,
            kappa_rna=kappa_rna,
            rho_r_bb=rho_r_bb,
            rho_d_bb=rho_d_bb,
            pi_g_prior=pi_g_prior,
        )
        # E-step: count channel now live (μ_d from the previous iteration; zero on
        # iter 1). Per-view FL-corrected L_eff (contained vs boundary).
        contained = estep_view(
            substrate.contained, ts_class, L_eff=region_eff_len, m_d_unspl_prev=m_d_cont, **shared
        )
        left = estep_view(
            substrate.left, ts_class, L_eff=boundary_eff, m_d_unspl_prev=m_d_left, **shared
        )
        right = estep_view(
            substrate.right, ts_class, L_eff=boundary_eff, m_d_unspl_prev=m_d_right, **shared
        )

        # Exposure (G4, physical L) + AMBIG sweep (PR 4b).
        exposure = exposure_posterior(
            contained.m_g,
            left.m_g,
            right.m_g,
            rho_0=rho_0,
            L_eff=l_phys,
            exposure_dispersion=exposure_dispersion,
        )
        swept = sweep_ambig_exposure(
            substrate,
            region_arrays,
            alloc_contained=contained,
            alloc_left=left,
            alloc_right=right,
            region_eff_len=region_eff_len,
            mu_fl=mu_fl,
            rho_0=rho_0,
            exposure_dispersion=exposure_dispersion,
            base_omega=exposure.omega,
            base_log_omega_var=exposure.log_omega_var,
            dec_region=dec_region,
        )
        omega = swept.omega
        log_omega_var = swept.log_omega_var

        # --- M-step (doc 03 §5): fit ρ_0, exposure_dispersion, ρ_d_bb (κ_rna/ρ_r_bb fixed) ---
        # PR 8: ρ_0 and the dispersion are fit from DECODABLE nodes only — an AMBIG
        # region's undecodable contained gDNA neither sets the global density nor the
        # dispersion; its decodable boundaries still contribute.
        rho_0 = update_rho_0(
            contained.m_g,
            left.m_g,
            right.m_g,
            omega,
            l_phys,
            dec_region=dec_region,
            dec_left_bnd=dec_left_bnd,
            dec_right_bnd=dec_right_bnd,
        )
        # Exposure dispersion: the proper EM M-step from the (pre-sweep) Gamma
        # exposure posteriors — NOT a count-NB fit. This is what breaks the φ
        # limit cycle (docs/acc_caljointmodel/calibration_oscillation_diagnosis.md).
        # Decodable regions only (PR 8): AMBIG pre-sweep ω is imputed-over garbage.
        exposure_dispersion = update_exposure_dispersion(
            exposure.omega[dec_region],
            exposure.log_omega_var[dec_region],
            floor=config.exposure_dispersion_floor,
        )
        k_plus_g = np.maximum(
            contained.k_sense.astype(np.float64) - kappa_rna * contained.m_d_unspl, 0.0
        )
        rho_d_bb = fit_rho_d_bb(k_plus_g, contained.m_g_unspl)
        pi_g_prior = update_pi_g_prior(omega, rho_0, region_eff_len, n_u_cont)

        # Carry the unspliced RNA mass into the next iteration's count channel.
        m_d_cont, m_d_left, m_d_right = contained.m_d_unspl, left.m_d_unspl, right.m_d_unspl

        # --- Convergence (doc 03 §9): max-region relative mass change ---
        m_g_tot = exposure.m_g_tot
        if r == 0:
            delta = 0.0
        elif m_g_tot_prev is None:
            delta = float(m_g_tot.max())  # change from the zero-gDNA init
        else:
            delta = float(np.max(np.abs(m_g_tot - m_g_tot_prev) / (m_g_tot_prev + 1.0)))
        mass_changes.append(delta)
        m_g_tot_prev = m_g_tot.copy()
        logger.debug(
            "[CAL iter %d] rho_0=%.4g exp_disp=%.4g rho_d_bb=%.4g delta=%.4g",
            it,
            rho_0,
            exposure_dispersion,
            rho_d_bb,
            delta,
        )
        if it >= 2 and delta < config.mass_rel_tol:
            converged = True
            break

    if not converged and r > 0:
        logger.warning(
            "[CAL] not converged in %d iterations (last mass change=%.4g > tol=%.4g)",
            n_iter,
            mass_changes[-1] if mass_changes else float("nan"),
            config.mass_rel_tol,
        )

    return CalibrationResult(
        mass_g_contained=contained.m_g,
        mass_d_contained=contained.m_d,
        mass_g_left=left.m_g,
        mass_d_left=left.m_d,
        mass_g_right=right.m_g,
        mass_d_right=right.m_d,
        omega=omega,
        log_omega_var=log_omega_var,
        rho_0=rho_0,
        exposure_dispersion=exposure_dispersion,
        rho_d_bb=rho_d_bb,
        kappa_rna=kappa_rna,
        rho_r_bb=rho_r_bb,
        n_iterations=n_iter,
        converged=converged,
        mass_change_history=np.array(mass_changes, dtype=np.float64),
        n_regions=r,
        config=config,
    )


__all__ = ["calibrate", "initial_hyperparameters"]
