"""Orchestrator: turn region counts + per-fragment lengths into a CalibrationResult.

v4: two-component mappability-aware EM mixture deconvolution.

See ``docs/calibration/calibration_v4_em_mixture_plan.md``.
"""

from __future__ import annotations

import logging
import math

import numpy as np
import pandas as pd

from ..frag_length_model import FragmentLengthModel
from ._em import run_em
from ._fl_model import build_gdna_fl_model
from ._result import CalibrationResult
from ._stats import compute_region_stats


logger = logging.getLogger(__name__)


def calibrate_gdna(
    region_counts: pd.DataFrame,
    fl_table: pd.DataFrame,
    region_df: pd.DataFrame,
    strand_specificity: float,
    *,
    mean_frag_len: float,
    intergenic_fl_model: FragmentLengthModel | None = None,
    fl_prior_ess: float = 1000.0,
    strand_specificity_noise_floor: float = 0.001,
    strand_specificity_ci_epsilon: float = 0.0,
    strand_llr_mode: str = "binomial",
    max_iterations: int = 50,
    convergence_tol: float = 1e-4,
    diagnostics: bool = False,
) -> CalibrationResult:
    """Calibrate gDNA contamination via the v4 mixture EM.

    Parameters
    ----------
    strand_specificity_noise_floor
        Biological floor on ``1 − ss`` (structural antisense rate under
        RNA).  Default 0.001.  See ``CalibrationConfig``.
    strand_specificity_ci_epsilon
        Measurement-uncertainty floor on ``1 − ss`` from the strand
        trainer (e.g. ``StrandModels.strand_specificity_ci_epsilon``).
        The effective floor used in the strand LLR is
        ``max(noise_floor, ci_epsilon)``.  Default 0.0 (rely on
        noise_floor alone).
    strand_llr_mode
        ``"binomial"`` (default) or ``"betabinom"``.  The latter uses
        a shared Beta-Binomial dispersion parameter κ, re-estimated
        inside the EM loop via the exact mixture marginal log-likelihood.
        See ``docs/calibration/calibration_v3_vs_v4_comparison.md``.
    """
    if len(region_counts) != len(region_df):
        raise ValueError("region_counts and region_df length mismatch")

    stats = compute_region_stats(region_counts, region_df)
    mappable_bp = stats["mappable_bp"]
    n_total = stats["n_total"].astype(np.float64)

    eps_floor_bio = float(strand_specificity_noise_floor)
    eps_floor_ci = float(strand_specificity_ci_epsilon)
    eps_floor = max(eps_floor_bio, eps_floor_ci)
    dominant = "ε_CI" if eps_floor_ci > eps_floor_bio else "ε_bio"
    logger.info(
        "Calibration strand LLR floor: ss_est=%.6f  ε_bio=%.4f  ε_CI=%.4f  "
        "→ ε_eff=%.4f (%s dominates)  max_per_antisense_llr=%.2f nats",
        float(strand_specificity), eps_floor_bio, eps_floor_ci, eps_floor,
        dominant, math.log(0.5 / max(eps_floor, 1e-12)),
    )

    fit = run_em(
        stats,
        fl_table,
        float(strand_specificity),
        mean_frag_len=float(mean_frag_len),
        strand_noise_floor=eps_floor,
        strand_llr_mode=str(strand_llr_mode),
        max_iterations=max_iterations,
        convergence_tol=convergence_tol,
        diagnostics=diagnostics,
    )

    region_e_gdna = fit.lam_G * mappable_bp.astype(np.float64)
    region_e_gdna = np.minimum(region_e_gdna, n_total)
    region_e_gdna = np.where(np.isfinite(region_e_gdna), region_e_gdna, 0.0)
    region_e_gdna = np.maximum(region_e_gdna, 0.0)

    gdna_fl_model = build_gdna_fl_model(
        fl_table,
        stats,
        region_weight=fit.gamma,
        strand_specificity=float(strand_specificity),
        intergenic_fl_model=intergenic_fl_model,
        fl_prior_ess=float(fl_prior_ess),
    )

    total_e = float(region_e_gdna.sum())
    total_n = float(n_total.sum())
    logger.info(
        "Calibration v4: lam_G=%.3e, pi=%.3f, pi_soft=%.3f, iters=%d%s, "
        "E[gDNA]=%.0f, frac=%.3f",
        fit.lam_G, fit.pi, fit.pi_soft, fit.n_iter,
        " (converged)" if fit.converged else " (NOT converged)",
        total_e, total_e / max(total_n, 1.0),
    )

    return CalibrationResult(
        region_e_gdna=region_e_gdna,
        region_n_total=n_total,
        gdna_fl_model=gdna_fl_model,
        lambda_gdna=float(fit.lam_G),
        strand_specificity=float(strand_specificity),
        region_gamma=fit.gamma,
        mu_R=float(fit.mu_R),
        sigma_R=float(fit.sigma_R),
        mixing_pi=float(fit.pi),
        mixing_pi_soft=float(fit.pi_soft),
        strand_used=fit.strand_used,
        strand_z=fit.strand_z,
        strand_llr_mode=fit.strand_llr_mode,
        kappa=fit.kappa,
        em_n_iter=fit.n_iter,
        em_converged=fit.converged,
        n_eligible=fit.n_eligible,
        n_soft=fit.n_soft,
        n_spliced_hard=fit.n_spliced_hard,
    )
