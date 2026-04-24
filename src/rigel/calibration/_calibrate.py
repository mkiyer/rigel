"""Orchestrator: turn region counts + per-fragment lengths into a CalibrationResult.

v5: two-component mappability-aware EM mixture deconvolution with
    Beta-Binomial(κ_G) vs ρ-integrated Beta-Binomial(κ_R) strand channel.

See ``docs/calibration/strand_channel_theory_and_redesign.md``.
"""

from __future__ import annotations

import logging

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
    kappa_gdna_min: float = 3.0,
    max_iterations: int = 50,
    convergence_tol: float = 1e-4,
    diagnostics: bool = False,
    capture_e: tuple[np.ndarray, np.ndarray] | None = None,
) -> CalibrationResult:
    """Calibrate gDNA contamination via the v5 mixture EM.

    Parameters
    ----------
    capture_e
        Optional ``(E_on, E_off)`` tuple of per-region mappable-bp
        attributions produced by
        :func:`rigel.calibration._annotate.annotate_capture_class`.
        When provided, the gDNA density channel is fit as a composite
        Poisson ``λ_on·E_on + λ_off·E_off`` via a nested EM; when
        ``None`` a single global density model is used (default for
        whole-genome / total-RNA libraries).
    """
    if len(region_counts) != len(region_df):
        raise ValueError("region_counts and region_df length mismatch")

    stats = compute_region_stats(region_counts, region_df)
    mappable_bp = stats["mappable_bp"]
    n_total = stats["n_total"].astype(np.float64)

    fit = run_em(
        stats,
        fl_table,
        float(strand_specificity),
        mean_frag_len=float(mean_frag_len),
        max_iterations=max_iterations,
        convergence_tol=convergence_tol,
        kappa_gdna_min=float(kappa_gdna_min),
        diagnostics=diagnostics,
        capture_e=capture_e,
    )

    # Per-region expected gDNA count under the composite-Poisson
    # model: rate_G_i = λ_on·E_on_i + λ_off·E_off_i.  In non-capture
    # mode ``capture_e`` is None and we fall back to the single-rate
    # ``λ_G · mappable_bp`` product (bit-identical to pre-composite).
    if fit.capture_class_mode and capture_e is not None:
        e_on, e_off = capture_e
        e_on = np.asarray(e_on, dtype=np.float64)
        e_off = np.asarray(e_off, dtype=np.float64)
        region_e_gdna = fit.lam_G_on * e_on + fit.lam_G_off * e_off
    else:
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
    if fit.capture_class_mode:
        enrichment = (fit.lam_G_on or 0.0) / max(fit.lam_G_off or 0.0, 1e-12)
        logger.info(
            "Calibration v5 (capture-class, composite-Poisson): "
            "λ_on=%.3e, λ_off=%.3e (%.1fx), λ_eff=%.3e, π=%.3f, π_soft=%.3f, "
            "κ_G=%.2f, κ_R=%.2f, iters=%d%s, E[gDNA]=%.0f, frac=%.3f",
            fit.lam_G_on, fit.lam_G_off, enrichment, fit.lam_G,
            fit.pi, fit.pi_soft, fit.kappa_G, fit.kappa_R, fit.n_iter,
            " (converged)" if fit.converged else " (NOT converged)",
            total_e, total_e / max(total_n, 1.0),
        )
    else:
        logger.info(
            "Calibration v5: lam_G=%.3e, pi=%.3f, pi_soft=%.3f, "
            "κ_G=%.2f, κ_R=%.2f, iters=%d%s, E[gDNA]=%.0f, frac=%.3f",
            fit.lam_G, fit.pi, fit.pi_soft, fit.kappa_G, fit.kappa_R, fit.n_iter,
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
        region_gamma_strand=fit.region_gamma_strand,
        mu_R=float(fit.mu_R),
        sigma_R=float(fit.sigma_R),
        mixing_pi=float(fit.pi),
        mixing_pi_soft=float(fit.pi_soft),
        kappa_G=float(fit.kappa_G),
        kappa_R=float(fit.kappa_R),
        em_n_iter=fit.n_iter,
        em_converged=fit.converged,
        n_eligible=fit.n_eligible,
        n_soft=fit.n_soft,
        n_spliced_hard=fit.n_spliced_hard,
        lam_G_on=float(fit.lam_G_on) if fit.capture_class_mode else None,
        lam_G_off=float(fit.lam_G_off) if fit.capture_class_mode else None,
        capture_class_mode=bool(fit.capture_class_mode),
    )
