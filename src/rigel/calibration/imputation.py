"""The cross-node imputation layer (CALIBRATION_ARCHITECTURE.md §1.2 + imputation_variance_model.md).

Reintroduced in Step 2 (the crude ``q_rna`` odds-propagation was removed in Step 1). The imputation gives a
node a PRIOR on its composition from a higher-precision neighbour:

    mean      = the neighbour's density (the IDENTITY — adjacent nodes of one field share a rate; the
                factor-1-under-uniform eff-len makes the slope 1 by construction, so there is no fitted slope)
    precision = 1 / (σ²_bio(μ) + ρ_src/L_src)
                = 1 / (the LEARNED biological dispersion  +  the predictor's COMPUTED Poisson sampling noise)

``σ²_bio`` is the Poisson-offset ``var~mean`` (:func:`variance_model.fit_pair_imputation_varmean`, which
subtracts the computed sampling floor and learns only the excess). ``ρ_src/L_src`` is the **second** §0
count→precision channel — a high-count predictor is a tight prior, a low-count one diffuse — never a
composition vote. The density→fraction conversion uses the definitional ``(eff/mass)²`` jacobian.

This module builds the **gDNA** prior (genomic propagation from the observable boundary-side gDNA crossings).
The transcript-structure-aware RNA prior is the sibling builder (Step 2c).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .variance_model import fit_pair_imputation_varmean

__all__ = ["ImputationPrior", "gdna_imputation_prior"]

_EPS = 1.0e-9


@dataclass(frozen=True, slots=True)
class ImputationPrior:
    """A per-node Gaussian prior on a component fraction: pull ``f`` toward ``frac`` at ``precision``.

    ``precision == 0`` ⇒ no prior at that node (it falls to the strand likelihood + the global foundation).
    Fraction-space (the sweep solves fractions); the density→fraction jacobian is already applied.
    """

    frac: np.ndarray  # μ_imp per node (the imputed mean, density→fraction-converted), 0 where absent
    precision: np.ndarray  # τ_imp per node (fraction-space), 0 where absent


def gdna_imputation_prior(
    region_density,
    left_density,
    right_density,
    *,
    region_eff_len,
    side_eff_len,
    mass_u,
    region_eligible,
    left_ok,
    right_ok,
    ref_id=None,
) -> ImputationPrior:
    """The gDNA imputation prior on ``f_g`` from the observable boundary-side gDNA crossings.

    All inputs are per-region (length R). ``region_density`` = the region's CURRENT gDNA density
    ``ρ_g = gdna_mass/eff_len`` (the iterating estimate); ``left/right_density`` = the flanking boundary
    sides' (fixed, deconv'd) gDNA crossing densities; ``side_eff_len`` = the per-side density length
    ``L_side``; ``region_eligible`` = the imputed dests (the non-count-observable exon/AMBIG regions);
    ``*_ok`` = that flank exists, is observable, carries gDNA mass.

    Fits ``σ²_bio(μ)`` (Poisson-offset) on the ``(observable side → eligible region)`` pairs, then for each
    eligible region combines its flanks precision-weighted: ``μ = side density`` (identity), ``τ =
    1/(σ²_bio(ρ_g) + ρ_side/L_side)``. Converts to the ``f_g`` prior via the ``(eff/mass)²`` jacobian.
    """
    rd = np.asarray(region_density, dtype=np.float64)
    ld = np.asarray(left_density, dtype=np.float64)
    rrd = np.asarray(right_density, dtype=np.float64)
    le = np.maximum(np.asarray(region_eff_len, dtype=np.float64), _EPS)
    se = np.maximum(np.asarray(side_eff_len, dtype=np.float64), _EPS)
    mu_mass = np.maximum(np.asarray(mass_u, dtype=np.float64), _EPS)
    elig = np.asarray(region_eligible, dtype=bool)
    lok = np.asarray(left_ok, dtype=bool)
    rok = np.asarray(right_ok, dtype=bool)

    # Per-node Poisson sampling variance ρ/L (= C/L²): region (the query node) + each side (the predictor).
    region_var_samp = rd / le
    left_var_samp = ld / se
    right_var_samp = rrd / se

    curve = fit_pair_imputation_varmean(
        rd, ld, rrd,
        region_eligible=elig, left_ok=lok, right_ok=rok,
        region_var_samp=region_var_samp, left_var_samp=left_var_samp, right_var_samp=right_var_samp,
        ref_id=ref_id,
    )

    # σ²_bio queried at the region's own current density (the fit-and-query-on-the-same-axis contract).
    sig2_bio = np.maximum(curve.predict(rd), 0.0)

    # Per-flank precision τ = 1/(σ²_bio + predictor sampling noise); mean = the side density (identity).
    tau_l = np.where(lok, 1.0 / np.maximum(sig2_bio + left_var_samp, _EPS), 0.0)
    tau_r = np.where(rok, 1.0 / np.maximum(sig2_bio + right_var_samp, _EPS), 0.0)
    tau_d = tau_l + tau_r  # precision-weighted combine of the (≤2) flanks
    mu_d = np.where(tau_d > 0.0, (tau_l * ld + tau_r * rrd) / np.maximum(tau_d, _EPS), 0.0)

    # Only eligible regions with an observable flank carry a prior.
    active = elig & (tau_d > 0.0)
    jac = le / mu_mass  # df/dρ — the density→fraction jacobian
    frac = np.where(active, np.clip(mu_d * jac, 0.0, 1.0), 0.0)
    precision = np.where(active, tau_d / np.maximum(jac, _EPS) ** 2, 0.0)  # τ_f = τ_ρ·(mass/eff)²
    return ImputationPrior(frac=frac.astype(np.float64), precision=precision.astype(np.float64))
