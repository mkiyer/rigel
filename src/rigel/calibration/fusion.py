"""Closed-form Fractional Mixture Allocation (FMA) fusion engine (PR 07 Phase 3).

Fuses per-region gDNA evidence from two independent channels:

* **Density** — count-per-opportunity geometry, exposed as a Negative-Binomial
  predictive mean/variance per region (``DensityEvidence``).
* **Strand** — orientation evidence conditional on the observed unspliced
  mass (``StrandGdnaEvidence``).

The fusion is a single vectorized closed-form pass: no grid integration, no
per-region loop, no exact/approximate branching. See
``docs/newcalib/pr07_strand_density_integration_v2.md`` for the derivation.

Properties (smooth in every input):

* ``π_r > 0, κ → 0``: ``D̂_r → T_anti + T_sense · 0.5π/(1−0.5π)``.
* ``π_r → 0, κ → 0``: ``D̂_r → T_anti`` (FMA does not assume RNA/gDNA
  symmetry; the textbook ``2·T_anti`` estimator requires non-trivial ``π_r``).
* ``κ → 0.5`` (unstranded): ``D̂_r = T_r · π_r = μ_density``.
* ``π_r → 1``: ``D̂_r → T_r``.
* Structurally absent strand (``applicable=False``): ``D̂_r = T_r · π_r``.
* ``T_r == 0``: ``D̂_r == 0`` with no NaN/inf anywhere in the dataclass.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .density_model import DensityEvidence
from .strand_evidence import StrandGdnaEvidence

__all__ = [
    "FUSED_STRAND_STRUCTURAL_ABSENT",
    "FUSED_STRAND_NEAR_UNSTRANDED",
    "FUSED_DENSITY_LOW_OPPORTUNITY",
    "FUSED_LOW_INFORMATION",
    "FUSED_CLIPPED_TO_TOTAL",
    "FUSED_DEGENERATE",
    "FusedRegionGdnaEvidence",
    "fuse_density_and_strand",
]


# Flag bits. uint16.
FUSED_STRAND_STRUCTURAL_ABSENT: int = 1 << 0
FUSED_STRAND_NEAR_UNSTRANDED:   int = 1 << 1
FUSED_DENSITY_LOW_OPPORTUNITY:  int = 1 << 2
FUSED_LOW_INFORMATION:          int = 1 << 3
FUSED_CLIPPED_TO_TOTAL:         int = 1 << 4
FUSED_DEGENERATE:               int = 1 << 5

_VAR_FLOOR: float = 1.0e-9
_LOW_INFO_FLOOR: float = 1.0e-12
# Used to flag near-unstranded regions for diagnostics only; the math is
# continuous and does not branch on this threshold.
_NEAR_UNSTRANDED_THR: float = 0.05


@dataclass(frozen=True, slots=True)
class FusedRegionGdnaEvidence:
    """Per-region FMA-fused gDNA mass posterior."""

    total_mass: np.ndarray              # float64[R]
    gdna_mass: np.ndarray               # float64[R] — D̂_r, clipped to [0, T_r]
    rna_mass: np.ndarray                # float64[R] — T_r − D̂_r
    variance: np.ndarray                # float64[R]
    pi_prior: np.ndarray                # float64[R]
    q_sense: np.ndarray                 # float64[R]
    q_anti: np.ndarray                  # float64[R]
    density_information: np.ndarray     # float64[R]
    strand_information: np.ndarray      # float64[R]
    density_weight: np.ndarray          # float64[R]
    strand_weight: np.ndarray           # float64[R]
    n_contained: np.ndarray             # uint32[R]
    n_left: np.ndarray                  # uint32[R]
    n_right: np.ndarray                 # uint32[R]
    flags: np.ndarray                   # uint16[R]


def fuse_density_and_strand(
    *,
    density_evidence: DensityEvidence,
    strand_evidence: StrandGdnaEvidence,
    n_contained: np.ndarray,
    n_left: np.ndarray,
    n_right: np.ndarray,
) -> FusedRegionGdnaEvidence:
    """Vectorized FMA fusion of density and strand channels.

    Parameters mirror the plan's contract: density's NB predictive provides
    ``μ_density`` and ``Var[density]``; strand provides ``T_sense``, ``T_anti``
    and the continuous gain factor ``g_lib · g_region · applicable``. The
    per-compartment uint32 support vectors are forwarded into the output for
    downstream consumers (background / exposure refit).
    """
    T = np.asarray(strand_evidence.n_total, dtype=np.float64)
    Ts = np.asarray(strand_evidence.n_sense, dtype=np.float64)
    Ta = np.asarray(strand_evidence.n_anti, dtype=np.float64)
    R = int(T.size)

    if density_evidence.mean_unbounded is None:
        raise ValueError("fuse_density_and_strand: DensityEvidence.mean_unbounded is required")
    if density_evidence.information is None:
        raise ValueError("fuse_density_and_strand: DensityEvidence.information is required")
    mu_d = np.asarray(density_evidence.mean_unbounded, dtype=np.float64)
    info_d = np.asarray(density_evidence.information, dtype=np.float64)

    # 1. Density-channel prior, safe against T_r == 0.
    pi = np.zeros(R, dtype=np.float64)
    np.divide(mu_d, T, out=pi, where=T > 0.0)
    np.clip(pi, 0.0, 1.0, out=pi)

    # 2. Effective per-region crosstalk. Continuous; structural applicability
    #    only — no count-driven branches.
    kappa_lib = float(min(strand_evidence.p_r1_sense, 1.0 - strand_evidence.p_r1_sense))
    g_lib = float(strand_evidence.global_info_scale)
    g_region = np.asarray(strand_evidence.region_info_gain, dtype=np.float64)
    applicable = np.asarray(strand_evidence.applicable, dtype=bool).astype(np.float64)
    g = g_lib * g_region * applicable
    kappa_eff = 0.5 - (0.5 - kappa_lib) * g  # κ_lib when g==1; 0.5 when g==0

    # 3. Allocation weights. At the degenerate corner ``denom == 0`` (which
    #    happens only at π=0 ∧ κ_eff=0 for denom_a, or π=0 ∧ κ_eff=1 for
    #    denom_s) the limiting value of ``q`` along the channel-only path
    #    is 1 — see plan §1.1. Setting q=1 at the sentinel reproduces the
    #    ``D̂ → T_anti`` strict limit and makes the function smoothly
    #    extend through the indeterminate point.
    denom_s = 0.5 * pi + (1.0 - kappa_eff) * (1.0 - pi)
    denom_a = 0.5 * pi + kappa_eff * (1.0 - pi)
    q_sense = np.where(
        denom_s > 0.0,
        np.divide(0.5 * pi, np.maximum(denom_s, _VAR_FLOOR)),
        1.0,
    )
    q_anti = np.where(
        denom_a > 0.0,
        np.divide(0.5 * pi, np.maximum(denom_a, _VAR_FLOOR)),
        1.0,
    )

    # 4. Mass allocation by linearity of expectation.
    D_hat = Ts * q_sense + Ta * q_anti
    clipped = (D_hat < 0.0) | (D_hat > T)
    np.clip(D_hat, 0.0, T, out=D_hat)
    R_hat = T - D_hat

    # 5. Posterior variance via Bernoulli on each strand at its mixture
    #    probability, scaled from physical fragment count N to fractional
    #    mass T by (T/N)^2.
    N_total = np.asarray(strand_evidence.support_total, dtype=np.float64)
    p_sense = np.full(R, 0.5, dtype=np.float64)
    np.divide(Ts, T, out=p_sense, where=T > 0.0)
    N_s = N_total * p_sense
    N_a = N_total * (1.0 - p_sense)
    scale = np.zeros(R, dtype=np.float64)
    np.divide(T * T, N_total * N_total, out=scale, where=N_total > 0.0)
    var_strand_raw = (N_s * q_sense * (1.0 - q_sense) + N_a * q_anti * (1.0 - q_anti)) * scale

    # Precisions: density from NB-predictive; strand from binomial mixture,
    # scaled by the per-region effective gain ``g`` so structurally absent
    # or low-info regions contribute zero strand information.
    tau_d = info_d
    info_s = np.zeros(R, dtype=np.float64)
    np.divide(1.0, var_strand_raw, out=info_s, where=var_strand_raw > _VAR_FLOOR)
    tau_s = info_s * g
    tau = tau_d + tau_s
    variance = np.where(tau > _LOW_INFO_FLOOR, 1.0 / np.maximum(tau, _VAR_FLOOR), np.inf)
    density_weight = np.where(tau > _LOW_INFO_FLOOR, tau_d / np.maximum(tau, _VAR_FLOOR), 0.0)
    strand_weight = np.where(tau > _LOW_INFO_FLOOR, tau_s / np.maximum(tau, _VAR_FLOOR), 0.0)

    # Finiteness invariant: every division above was guarded.
    assert np.all(np.isfinite(D_hat)), "fuse_density_and_strand produced non-finite gdna_mass"
    assert np.all(np.isfinite(q_sense))
    assert np.all(np.isfinite(q_anti))

    # Flags (diagnostic only; no consumer branches on them).
    flags = np.zeros(R, dtype=np.uint16)
    flags[np.asarray(strand_evidence.structural_absent, dtype=bool)] |= FUSED_STRAND_STRUCTURAL_ABSENT
    flags[abs(0.5 - kappa_eff) < _NEAR_UNSTRANDED_THR] |= FUSED_STRAND_NEAR_UNSTRANDED
    if density_evidence.applicable is not None:
        flags[~np.asarray(density_evidence.applicable, dtype=bool)] |= FUSED_DENSITY_LOW_OPPORTUNITY
    flags[tau <= _LOW_INFO_FLOOR] |= FUSED_LOW_INFORMATION
    flags[clipped] |= FUSED_CLIPPED_TO_TOTAL
    flags[T <= 0.0] |= FUSED_DEGENERATE

    return FusedRegionGdnaEvidence(
        total_mass=T,
        gdna_mass=D_hat,
        rna_mass=R_hat,
        variance=variance,
        pi_prior=pi,
        q_sense=q_sense,
        q_anti=q_anti,
        density_information=tau_d.astype(np.float64, copy=False),
        strand_information=tau_s,
        density_weight=density_weight,
        strand_weight=strand_weight,
        n_contained=np.asarray(n_contained, dtype=np.uint32, order="C"),
        n_left=np.asarray(n_left, dtype=np.uint32, order="C"),
        n_right=np.asarray(n_right, dtype=np.uint32, order="C"),
        flags=flags,
    )
