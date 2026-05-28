"""Calibration E/M iteration for the two-state expression region model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np
from scipy.special import logsumexp

from ._arrays import RegionArrays
from .background_model import BackgroundModel
from .boundaries import BoundaryTable
from .boundary_model import (
    BoundaryLocalPosterior,
    build_boundary_local_posterior,
    predict_contained_gdna_from_excess,
)
from .boundary_sweep import BoundarySweepResult, run_boundary_sweep
from .density_observation import DensityObservation
from .exposure import RegionExposure, estimate_region_exposure
from .latent_states import (
    N_STATES,
    STATE_EXPRESSED,
    STATE_UNEXPRESSED,
    build_logbf_expression,
    build_logbf_strand,
    build_state_log_prior,
    build_state_log_tensor,
    normalize_state_log_tensor,
)
from .strand_deconv import RegionGdnaChannelEstimate

__all__ = [
    "FLAG_STATE_AMBIGUOUS",
    "FLAG_STRAND_UNINFORMATIVE",
    "FLAG_BOUNDARY_SPARSE",
    "FLAG_EXPRESSED_UNCERTAIN",
    "FLAG_EXACT_STRAND_POSTERIOR",
    "FLAG_M_IMPUTED_FROM_BACKGROUND",
    "FLAG_M_CLIPPED_TO_TOTAL",
    "METHOD_STRAND",
    "METHOD_BOUNDARY",
    "METHOD_BACKGROUND_FALLBACK",
    "RegionUnsplicedMass",
    "BackgroundDensity",
    "RegionCalibration",
    "CalibrationStepResult",
    "build_region_unspliced_mass",
    "estimate_background_density",
    "calibration_e_step",
    "calibration_m_step",
    "run_calibration_iteration",
]

FLAG_STATE_AMBIGUOUS: int = 1 << 0
FLAG_STRAND_UNINFORMATIVE: int = 1 << 1
FLAG_BOUNDARY_SPARSE: int = 1 << 2
FLAG_EXPRESSED_UNCERTAIN: int = 1 << 3
FLAG_EXACT_STRAND_POSTERIOR: int = 1 << 4

# PR 03: per-region flags on RegionUnsplicedMass.
FLAG_M_IMPUTED_FROM_BACKGROUND: int = 1 << 0
FLAG_M_CLIPPED_TO_TOTAL: int = 1 << 1

# PR 03: three-tier M_r method constants.
METHOD_STRAND: int = 1
METHOD_BOUNDARY: int = 2
METHOD_BACKGROUND_FALLBACK: int = 3
_METHOD_VALUES: frozenset[int] = frozenset(
    {METHOD_STRAND, METHOD_BOUNDARY, METHOD_BACKGROUND_FALLBACK}
)

_EPS: float = 1.0e-12


@dataclass(frozen=True, slots=True)
class RegionUnsplicedMass:
    """Mass-conserving unspliced decomposition per region.

    All primary mass tensors are ``float64`` so that
    ``gdna_mass + rna_mass == total_mass`` is an exact identity (no
    tolerance). The three-tier ``M_r`` hierarchy (``METHOD_STRAND`` /
    ``METHOD_BOUNDARY`` / ``METHOD_BACKGROUND_FALLBACK``) writes
    ``method[r]`` to record which tier owned the estimate.
    """

    # Primary mass tensors (float64 for exact conservation).
    total_mass: np.ndarray         # float64[R]  T_r
    gdna_mass: np.ndarray          # float64[R]  M_r
    rna_mass: np.ndarray           # float64[R]  R_r = T_r - M_r

    # Geometric denominator.
    region_size_bp: np.ndarray     # float64[R]  end - start in bp

    # Effective sample size (from PR 02a).
    unspliced_counts: np.ndarray   # uint64[R]   physical unspliced fragment support

    # Provenance + quality.
    method: np.ndarray             # uint8[R]    METHOD_STRAND / METHOD_BOUNDARY /
                                   #             METHOD_BACKGROUND_FALLBACK
    precision: np.ndarray          # float64[R]  estimator reliability (M-step weight)
    flags: np.ndarray              # uint16[R]   diagnostic bits

    def __post_init__(self) -> None:
        total = np.ascontiguousarray(np.asarray(self.total_mass, dtype=np.float64))
        if total.ndim != 1:
            raise ValueError(f"total_mass must be 1D; got shape {total.shape}.")
        region_count = int(total.shape[0])
        object.__setattr__(self, "total_mass", total)

        for name in ("gdna_mass", "rna_mass", "region_size_bp", "precision"):
            values = np.ascontiguousarray(np.asarray(getattr(self, name), dtype=np.float64))
            if values.shape != (region_count,):
                raise ValueError(
                    f"{name} must have shape ({region_count},); got {values.shape}."
                )
            object.__setattr__(self, name, values)

        counts = np.ascontiguousarray(np.asarray(self.unspliced_counts, dtype=np.uint64))
        if counts.shape != (region_count,):
            raise ValueError(
                f"unspliced_counts must have shape ({region_count},); got {counts.shape}."
            )
        object.__setattr__(self, "unspliced_counts", counts)

        method = np.ascontiguousarray(np.asarray(self.method, dtype=np.uint8))
        if method.shape != (region_count,):
            raise ValueError(f"method must have shape ({region_count},); got {method.shape}.")
        if region_count > 0:
            valid = np.isin(method, np.array(sorted(_METHOD_VALUES), dtype=np.uint8))
            if not bool(valid.all()):
                bad = np.unique(method[~valid]).tolist()
                raise ValueError(
                    f"method must be one of {sorted(_METHOD_VALUES)}; got values {bad}."
                )
        object.__setattr__(self, "method", method)

        flags = np.ascontiguousarray(np.asarray(self.flags, dtype=np.uint16))
        if flags.shape != (region_count,):
            raise ValueError(f"flags must have shape ({region_count},); got {flags.shape}.")
        object.__setattr__(self, "flags", flags)

        if np.any(total < 0.0):
            raise ValueError("total_mass must be non-negative.")
        if np.any(self.gdna_mass < 0.0):
            raise ValueError("gdna_mass must be non-negative.")
        if np.any(self.gdna_mass > total):
            raise ValueError("gdna_mass must not exceed total_mass per region.")
        if np.any(self.rna_mass < 0.0):
            raise ValueError("rna_mass must be non-negative.")
        if np.any(self.rna_mass > total):
            raise ValueError("rna_mass must not exceed total_mass per region.")
        # Mass conservation must be exact in float64; the construction code is
        # responsible for deriving rna_mass = total_mass - gdna_mass.
        if not np.array_equal(self.gdna_mass + self.rna_mass, total):
            raise ValueError("gdna_mass + rna_mass must equal total_mass exactly.")
        if np.any(self.region_size_bp <= 0.0):
            raise ValueError("region_size_bp must be strictly positive.")
        if not np.all(np.isfinite(self.precision)) or np.any(self.precision < 0.0):
            raise ValueError("precision must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class BackgroundDensity:
    """Library-wide gDNA fragment-mass density (PR 03).

    Carries the robust geometric-mean view (``rho0_mean``, ``log_dispersion``)
    and the Gamma conjugacy view (``alpha0``, ``beta0``). The two are kept
    consistent by construction: ``rho0_mean == alpha0 / beta0``.
    """

    rho0_mean: float            # gDNA fragment mass per bp
    alpha0: float               # Gamma shape (alpha_floor + Sum w*gdna_mass after rescale)
    beta0: float                # Gamma rate  (beta_floor  + Sum w*region_size_bp after rescale)
    log_dispersion: float       # weighted sigma-equivalent of log(rho_hat_r), in nats
    n_effective_regions: float  # Sum of weights over pool
    n_regions_in_pool: int      # raw pool size
    method_histogram: tuple     # (n_tier1, n_tier2, n_tier3_excluded)
    fit_status: str             # see Section 6.6 of PR03 plan

    _VALID_STATUSES: ClassVar[tuple[str, ...]] = (
        "converged",
        "iterating",
        "fallback_bootstrap",
        "prior_only",
    )

    def __post_init__(self) -> None:
        rho0 = float(self.rho0_mean)
        alpha0 = float(self.alpha0)
        beta0 = float(self.beta0)
        log_disp = float(self.log_dispersion)
        n_eff = float(self.n_effective_regions)
        n_pool = int(self.n_regions_in_pool)
        if not np.isfinite(rho0) or rho0 <= 0.0:
            raise ValueError(f"rho0_mean must be finite and > 0; got {self.rho0_mean!r}.")
        if not np.isfinite(alpha0) or alpha0 <= 0.0:
            raise ValueError(f"alpha0 must be finite and > 0; got {self.alpha0!r}.")
        if not np.isfinite(beta0) or beta0 <= 0.0:
            raise ValueError(f"beta0 must be finite and > 0; got {self.beta0!r}.")
        if not np.isfinite(log_disp) or log_disp < 0.0:
            raise ValueError(
                f"log_dispersion must be finite and >= 0; got {self.log_dispersion!r}."
            )
        if not np.isfinite(n_eff) or n_eff < 0.0:
            raise ValueError(
                f"n_effective_regions must be finite and >= 0; got {self.n_effective_regions!r}."
            )
        if n_pool < 0:
            raise ValueError(
                f"n_regions_in_pool must be >= 0; got {self.n_regions_in_pool!r}."
            )
        hist = tuple(int(x) for x in self.method_histogram)
        if len(hist) != 3 or any(x < 0 for x in hist):
            raise ValueError(
                "method_histogram must be a 3-tuple of non-negative ints; "
                f"got {self.method_histogram!r}."
            )
        status = str(self.fit_status)
        if status not in self._VALID_STATUSES:
            raise ValueError(
                f"fit_status must be one of {self._VALID_STATUSES}; got {self.fit_status!r}."
            )
        object.__setattr__(self, "rho0_mean", rho0)
        object.__setattr__(self, "alpha0", alpha0)
        object.__setattr__(self, "beta0", beta0)
        object.__setattr__(self, "log_dispersion", log_disp)
        object.__setattr__(self, "n_effective_regions", n_eff)
        object.__setattr__(self, "n_regions_in_pool", n_pool)
        object.__setattr__(self, "method_histogram", hist)
        object.__setattr__(self, "fit_status", status)

    @classmethod
    def from_bootstrap(cls, model: BackgroundModel) -> "BackgroundDensity":
        """Seed the first iteration's BackgroundDensity from a bootstrap fit."""
        return cls(
            rho0_mean=float(model.rho_off_mean),
            alpha0=float(model.rho_off_alpha),
            beta0=float(model.rho_off_beta),
            log_dispersion=float(np.log(10.0)),
            n_effective_regions=0.0,
            n_regions_in_pool=0,
            method_histogram=(0, 0, 0),
            fit_status="fallback_bootstrap",
        )


@dataclass(frozen=True, slots=True)
class RegionCalibration:
    """Region-level two-state calibration output for downstream wiring."""

    p_states: np.ndarray
    mu_gdna: np.ndarray
    upper_gdna: np.ndarray
    rna_lower: np.ndarray
    region_unspliced_mass: RegionUnsplicedMass
    background_density: BackgroundDensity
    region_exposure: RegionExposure
    kappa_d: float | None
    n_passes: int
    converged: bool
    flags: np.ndarray
    pass_diagnostics: tuple[dict[str, object], ...]
    background_model: BackgroundModel | None = None
    boundary_local: BoundaryLocalPosterior | None = None
    boundary_sweep: BoundarySweepResult | None = None

    def __post_init__(self) -> None:
        p_states = np.asarray(self.p_states, dtype=np.float32)
        if p_states.ndim != 2 or p_states.shape[1] != N_STATES:
            raise ValueError(f"p_states must have shape (R, {N_STATES}); got {p_states.shape}.")
        region_count = int(p_states.shape[0])
        row_sums = p_states.sum(axis=1)
        if not np.allclose(row_sums, 1.0, rtol=1.0e-5, atol=1.0e-5):
            raise ValueError("p_states rows must sum to 1 within tolerance.")
        object.__setattr__(self, "p_states", p_states)

        for field_name in ("mu_gdna", "upper_gdna", "rna_lower"):
            values = _as_float32_vector(field_name, getattr(self, field_name), region_count)
            object.__setattr__(self, field_name, values)

        if np.any(self.mu_gdna < 0.0):
            raise ValueError("mu_gdna must be non-negative.")
        if np.any(self.upper_gdna + 1.0e-5 < self.mu_gdna):
            raise ValueError("upper_gdna must be >= mu_gdna within tolerance.")
        if np.any(self.rna_lower < 0.0):
            raise ValueError("rna_lower must be non-negative.")

        if not isinstance(self.region_unspliced_mass, RegionUnsplicedMass):
            raise TypeError(
                "region_unspliced_mass must be a RegionUnsplicedMass; "
                f"got {type(self.region_unspliced_mass).__name__}."
            )
        if self.region_unspliced_mass.total_mass.shape != (region_count,):
            raise ValueError(
                "region_unspliced_mass arrays must match p_states region count; "
                f"got {self.region_unspliced_mass.total_mass.shape}, "
                f"expected {(region_count,)}."
            )
        if not isinstance(self.background_density, BackgroundDensity):
            raise TypeError(
                "background_density must be a BackgroundDensity; "
                f"got {type(self.background_density).__name__}."
            )
        if not isinstance(self.region_exposure, RegionExposure):
            raise TypeError(
                "region_exposure must be a RegionExposure; "
                f"got {type(self.region_exposure).__name__}."
            )
        if self.region_exposure.omega.shape != (region_count,):
            raise ValueError(
                "region_exposure arrays must match p_states region count; "
                f"got {self.region_exposure.omega.shape}, expected {(region_count,)}."
            )

        flags = np.asarray(self.flags, dtype=np.uint16)
        if flags.shape != (region_count,):
            raise ValueError(f"flags must have shape ({region_count},); got {flags.shape}.")
        object.__setattr__(self, "flags", flags)

        if self.kappa_d is not None and (
            not np.isfinite(self.kappa_d) or float(self.kappa_d) <= 0.0
        ):
            raise ValueError(f"kappa_d must be None or finite and positive; got {self.kappa_d!r}.")
        if int(self.n_passes) < 1:
            raise ValueError(f"n_passes must be >= 1; got {self.n_passes!r}.")
        diagnostics = tuple(dict(item) for item in self.pass_diagnostics)
        object.__setattr__(self, "pass_diagnostics", diagnostics)

    @property
    def rho_off(self) -> float:
        """Library-wide gDNA density (mean of ``background_density``)."""
        return float(self.background_density.rho0_mean)

    @property
    def p_unexpressed(self) -> np.ndarray:
        return self.p_states[:, STATE_UNEXPRESSED]

    @property
    def p_expressed(self) -> np.ndarray:
        return self.p_states[:, STATE_EXPRESSED]


@dataclass(frozen=True, slots=True)
class CalibrationStepResult:
    """Single E-step output before scalar M-step refits."""

    p_states: np.ndarray
    mu_gdna: np.ndarray
    upper_gdna: np.ndarray
    rna_lower: np.ndarray
    region_unspliced_mass: RegionUnsplicedMass
    region_exposure: RegionExposure
    flags: np.ndarray
    local_posterior: BoundaryLocalPosterior
    sweep: BoundarySweepResult
    log_tensor: np.ndarray
    sum_log_evidence: float

    @property
    def p_unexpressed(self) -> np.ndarray:
        return self.p_states[:, STATE_UNEXPRESSED]

    @property
    def p_expressed(self) -> np.ndarray:
        return self.p_states[:, STATE_EXPRESSED]



def _as_float32_vector(name: str, values: np.ndarray, region_count: int) -> np.ndarray:
    array = np.asarray(values, dtype=np.float32)
    if array.shape != (region_count,):
        raise ValueError(f"{name} must have shape ({region_count},); got {array.shape}.")
    return array


def _as_float64_vector(name: str, values: np.ndarray, region_count: int) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (region_count,):
        raise ValueError(f"{name} must have shape ({region_count},); got {array.shape}.")
    return array


def _validate_region_inputs(
    region_arrays: RegionArrays,
    observation: DensityObservation,
    boundaries: BoundaryTable,
    background: BackgroundModel,
    strand_channels: RegionGdnaChannelEstimate | None,
) -> int:
    region_count = int(np.asarray(observation.contained_leff).shape[0])
    if np.asarray(region_arrays.signature).shape != (region_count,):
        raise ValueError(
            "region_arrays.signature must match observation region count; "
            f"got {np.asarray(region_arrays.signature).shape}, expected {(region_count,)}."
        )
    if int(boundaries.ref_region_offsets[-1]) != region_count:
        raise ValueError(
            "boundaries and observation disagree on region count; "
            f"got {int(boundaries.ref_region_offsets[-1])} and {region_count}."
        )
    if np.asarray(background.seed_mask).shape != (region_count,):
        raise ValueError(
            "background.seed_mask must match observation region count; "
            f"got {np.asarray(background.seed_mask).shape}, expected {(region_count,)}."
        )
    if strand_channels is not None:
        _as_float64_vector("strand_channels.contained_mean", strand_channels.contained_mean, region_count)
    return region_count


def _validate_local_posterior(local_posterior: BoundaryLocalPosterior, region_count: int) -> None:
    _as_float64_vector("local_posterior.alpha_excess", local_posterior.alpha_excess, region_count)
    _as_float64_vector("local_posterior.beta_excess", local_posterior.beta_excess, region_count)


# ---------------------------------------------------------------------------
# PR 03: Three-tier M_r hierarchy
# ---------------------------------------------------------------------------

_STRAND_RELIABILITY_EPS: float = 1.0e-6


def build_region_unspliced_mass(
    observation: DensityObservation,
    *,
    region_size_bp: np.ndarray,
    unspliced_counts: np.ndarray,
    strand_channels: RegionGdnaChannelEstimate | None,
    local_posterior: BoundaryLocalPosterior,
    sweep: BoundarySweepResult,
    background_density: BackgroundDensity,
    force_zero_gdna: bool = False,
) -> RegionUnsplicedMass:
    """Three-tier M_r hierarchy (Section 4 of PR03 plan v3).

    For each region the first tier whose condition fires owns ``M_r``:

    - Tier 1 (METHOD_STRAND): any strand-channel reliability >= 1e-6.
    - Tier 2 (METHOD_BOUNDARY): boundary excess evidence
      (``alpha_excess > 0 or beta_excess > 0``).
    - Tier 3 (METHOD_BACKGROUND_FALLBACK): ``M_r = clip(rho0_mean * region_size_bp, 0, T_r)``.

    All primary mass tensors are float64; ``rna_mass = total_mass - gdna_mass`` is
    derived after clipping so conservation is exact.
    """
    total = _as_float64_vector(
        "observation.observed_compatible_count",
        observation.observed_compatible_count,
        int(np.asarray(observation.observed_compatible_count).shape[0]),
    )
    region_count = int(total.shape[0])
    region_bp = _as_float64_vector("region_size_bp", region_size_bp, region_count)
    counts = np.ascontiguousarray(np.asarray(unspliced_counts, dtype=np.uint64))
    if counts.shape != (region_count,):
        raise ValueError(
            f"unspliced_counts must have shape ({region_count},); got {counts.shape}."
        )
    _validate_local_posterior(local_posterior, region_count)
    for name in ("mu_sweep", "alpha_excess", "beta_excess"):
        src = sweep if name == "mu_sweep" else local_posterior
        _as_float64_vector(f"sweep/local.{name}", getattr(src, name), region_count)

    gdna = np.zeros(region_count, dtype=np.float64)
    precision = np.zeros(region_count, dtype=np.float64)
    method = np.full(region_count, METHOD_BACKGROUND_FALLBACK, dtype=np.uint8)
    flags = np.zeros(region_count, dtype=np.uint16)
    tier1_gdna_raw = np.zeros(region_count, dtype=np.float64)

    if bool(force_zero_gdna):
        flags |= np.uint16(FLAG_M_IMPUTED_FROM_BACKGROUND)
        return RegionUnsplicedMass(
            total_mass=total,
            gdna_mass=gdna,
            rna_mass=total - gdna,
            region_size_bp=region_bp,
            unspliced_counts=counts,
            method=method,
            precision=precision,
            flags=flags,
        )

    # ---- Tier 1: strand deconvolution ----
    tier1_mask = np.zeros(region_count, dtype=bool)
    if strand_channels is not None:
        contained_mean = _as_float64_vector(
            "strand_channels.contained_mean",
            strand_channels.contained_mean,
            region_count,
        )
        left_mean = _as_float64_vector(
            "strand_channels.boundary_left_mean",
            strand_channels.boundary_left_mean,
            region_count,
        )
        right_mean = _as_float64_vector(
            "strand_channels.boundary_right_mean",
            strand_channels.boundary_right_mean,
            region_count,
        )
        contained_rel = np.clip(
            _as_float64_vector(
                "strand_channels.contained_reliability",
                strand_channels.contained_reliability,
                region_count,
            ),
            0.0,
            1.0,
        )
        left_rel = np.clip(
            _as_float64_vector(
                "strand_channels.boundary_left_reliability",
                strand_channels.boundary_left_reliability,
                region_count,
            ),
            0.0,
            1.0,
        )
        right_rel = np.clip(
            _as_float64_vector(
                "strand_channels.boundary_right_reliability",
                strand_channels.boundary_right_reliability,
                region_count,
            ),
            0.0,
            1.0,
        )
        contained_prec = _as_float64_vector(
            "strand_channels.contained_precision",
            strand_channels.contained_precision,
            region_count,
        )
        left_prec = _as_float64_vector(
            "strand_channels.boundary_left_precision",
            strand_channels.boundary_left_precision,
            region_count,
        )
        right_prec = _as_float64_vector(
            "strand_channels.boundary_right_precision",
            strand_channels.boundary_right_precision,
            region_count,
        )

        tier1_mask = (
            (contained_rel >= _STRAND_RELIABILITY_EPS)
            | (left_rel >= _STRAND_RELIABILITY_EPS)
            | (right_rel >= _STRAND_RELIABILITY_EPS)
        )

        tier1_gdna = (
            contained_mean * contained_rel
            + left_mean * left_rel
            + right_mean * right_rel
        )
        tier1_gdna_raw = tier1_gdna
        tier1_prec = np.maximum.reduce(
            [
                contained_prec * contained_rel,
                left_prec * left_rel,
                right_prec * right_rel,
            ]
        )
        strand_flags = np.asarray(strand_channels.flags, dtype=np.uint16)

        gdna = np.where(tier1_mask, np.clip(tier1_gdna, 0.0, total), gdna)
        precision = np.where(tier1_mask, tier1_prec, precision)
        method = np.where(tier1_mask, np.uint8(METHOD_STRAND), method)
        flags = np.where(tier1_mask, strand_flags, flags).astype(np.uint16, copy=False)

    # ---- Tier 2: boundary sweep imputation (for regions Tier 1 did not own) ----
    alpha_excess = np.asarray(local_posterior.alpha_excess, dtype=np.float64)
    beta_excess = np.asarray(local_posterior.beta_excess, dtype=np.float64)
    mu_sweep = np.asarray(sweep.mu_sweep, dtype=np.float64)
    tier2_candidate = (~tier1_mask) & ((alpha_excess > 0.0) | (beta_excess > 0.0))

    tier2_gdna = np.clip(mu_sweep, 0.0, total)
    # PR03 plan referenced `sweep.contained_gdna_precision`, which is not a field
    # on BoundarySweepResult. Use boundary-excess ESS (alpha + beta) as a
    # proportional precision proxy: this is the underlying sample-size weight of
    # the sweep estimator at that region.
    tier2_prec = alpha_excess + beta_excess

    gdna = np.where(tier2_candidate, tier2_gdna, gdna)
    precision = np.where(tier2_candidate, tier2_prec, precision)
    method = np.where(tier2_candidate, np.uint8(METHOD_BOUNDARY), method)

    # ---- Tier 3: background fallback ----
    tier3_mask = ~(tier1_mask | tier2_candidate)
    rho0 = float(background_density.rho0_mean)
    tier3_gdna_raw = rho0 * region_bp
    tier3_gdna = np.clip(tier3_gdna_raw, 0.0, total)
    gdna = np.where(tier3_mask, tier3_gdna, gdna)
    # precision already 0.0 for Tier 3 rows.
    flags = np.where(
        tier3_mask,
        (flags | np.uint16(FLAG_M_IMPUTED_FROM_BACKGROUND)).astype(np.uint16),
        flags,
    ).astype(np.uint16, copy=False)

    # Detect clipping. Strand-tier raw value is tier1_gdna; boundary-tier raw is
    # mu_sweep; background-tier raw is rho0 * region_bp. We compare each region's
    # raw value against total to flag clipping.
    raw_value = np.zeros(region_count, dtype=np.float64)
    raw_value = np.where(tier1_mask, tier1_gdna_raw, raw_value)
    raw_value = np.where(tier2_candidate, mu_sweep, raw_value)
    raw_value = np.where(tier3_mask, tier3_gdna_raw, raw_value)
    clipped = (raw_value > total) | (raw_value < 0.0)
    flags = np.where(
        clipped,
        (flags | np.uint16(FLAG_M_CLIPPED_TO_TOTAL)).astype(np.uint16),
        flags,
    ).astype(np.uint16, copy=False)

    # Final sanitisation: replace any non-finite gdna with 0.0 before deriving rna.
    gdna = np.nan_to_num(gdna, nan=0.0, posinf=0.0, neginf=0.0)
    gdna = np.clip(gdna, 0.0, total)
    rna = total - gdna  # exact in float64 because gdna <= total by clipping.

    precision = np.nan_to_num(precision, nan=0.0, posinf=0.0, neginf=0.0)
    precision = np.maximum(precision, 0.0)

    return RegionUnsplicedMass(
        total_mass=total,
        gdna_mass=gdna,
        rna_mass=rna,
        region_size_bp=region_bp,
        unspliced_counts=counts,
        method=method,
        precision=precision,
        flags=flags,
    )


# ---------------------------------------------------------------------------
# PR 03: Robust geomean + Huber rho0 estimator (Section 6.3-6.7)
# ---------------------------------------------------------------------------

_LOG_DISPERSION_FLOOR_DEFAULT: float = float(np.log(1.1))


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    """Weighted median (lower-median convention) over finite-weight entries."""
    mask = weights > 0.0
    if not np.any(mask):
        return float("nan")
    v = values[mask]
    w = weights[mask]
    order = np.argsort(v, kind="mergesort")
    v = v[order]
    w = w[order]
    cum = np.cumsum(w)
    half = 0.5 * cum[-1]
    idx = int(np.searchsorted(cum, half, side="left"))
    if idx >= v.shape[0]:
        idx = v.shape[0] - 1
    return float(v[idx])


def estimate_background_density(
    region_unspliced_mass: RegionUnsplicedMass,
    p_unexpressed: np.ndarray,
    *,
    previous: BackgroundDensity,
    damping: float = 0.5,
    huber_k: float = 1.5,
    alpha_floor: float = 1.0,
    beta_floor: float = 1.0,
    rho_floor: float = 1.0e-12,
    convergence_log_tol: float = 0.01,
    log_dispersion_floor: float = _LOG_DISPERSION_FLOOR_DEFAULT,
) -> BackgroundDensity:
    """Robust geometric-mean + Huber refit of the off-target background density.

    Implements Section 6.1-6.7 of the PR 03 plan.
    On an empty pool the previous ``BackgroundDensity`` is carried unchanged
    (``fit_status="fallback_bootstrap"``).
    """
    if not np.isfinite(damping) or not 0.0 <= float(damping) <= 1.0:
        raise ValueError(f"damping must be finite and in [0, 1]; got {damping!r}.")
    if huber_k <= 0.0 or not np.isfinite(huber_k):
        raise ValueError(f"huber_k must be a positive finite float; got {huber_k!r}.")
    if alpha_floor <= 0.0 or beta_floor <= 0.0:
        raise ValueError(
            f"alpha_floor/beta_floor must be positive; got {alpha_floor!r}, {beta_floor!r}."
        )
    if rho_floor <= 0.0:
        raise ValueError(f"rho_floor must be positive; got {rho_floor!r}.")
    if convergence_log_tol < 0.0:
        raise ValueError(
            f"convergence_log_tol must be non-negative; got {convergence_log_tol!r}."
        )
    if log_dispersion_floor <= 0.0:
        raise ValueError(
            f"log_dispersion_floor must be positive; got {log_dispersion_floor!r}."
        )

    gdna = np.asarray(region_unspliced_mass.gdna_mass, dtype=np.float64)
    region_bp = np.asarray(region_unspliced_mass.region_size_bp, dtype=np.float64)
    counts = np.asarray(region_unspliced_mass.unspliced_counts, dtype=np.uint64)
    method = np.asarray(region_unspliced_mass.method, dtype=np.uint8)
    precision = np.asarray(region_unspliced_mass.precision, dtype=np.float64)

    region_count = int(gdna.shape[0])
    p_unx = np.asarray(p_unexpressed, dtype=np.float64)
    if p_unx.shape != (region_count,):
        raise ValueError(
            f"p_unexpressed must have shape ({region_count},); got {p_unx.shape}."
        )

    n_strand = int(np.count_nonzero(method == METHOD_STRAND))
    n_boundary = int(np.count_nonzero(method == METHOD_BOUNDARY))
    n_background = int(np.count_nonzero(method == METHOD_BACKGROUND_FALLBACK))
    method_histogram = (n_strand, n_boundary, n_background)

    pool_mask = (
        ((method == METHOD_STRAND) | (method == METHOD_BOUNDARY))
        & (counts >= np.uint64(1))
        & (region_bp >= 1.0)
    )
    n_in_pool = int(np.count_nonzero(pool_mask))

    prev_rho0 = max(float(previous.rho0_mean), rho_floor)

    if n_in_pool == 0:
        return BackgroundDensity(
            rho0_mean=prev_rho0,
            alpha0=float(previous.alpha0),
            beta0=float(previous.beta0),
            log_dispersion=float(previous.log_dispersion),
            n_effective_regions=0.0,
            n_regions_in_pool=0,
            method_histogram=method_histogram,
            fit_status="fallback_bootstrap",
        )

    # ---- Section 6.3.1: Bayesian-shrunk per-region density ----
    pseudo_mass = alpha_floor
    pseudo_size = alpha_floor / prev_rho0
    rho_hat = (gdna + pseudo_mass) / (region_bp + pseudo_size)
    rho_hat = np.maximum(rho_hat, rho_floor)
    log_rho_hat = np.log(rho_hat)

    weights = precision * counts.astype(np.float64) * p_unx
    weights = np.where(pool_mask, weights, 0.0)
    weights = np.nan_to_num(weights, nan=0.0, posinf=0.0, neginf=0.0)
    weights = np.maximum(weights, 0.0)

    w_sum = float(weights.sum())
    if w_sum <= 0.0:
        return BackgroundDensity(
            rho0_mean=prev_rho0,
            alpha0=float(previous.alpha0),
            beta0=float(previous.beta0),
            log_dispersion=float(previous.log_dispersion),
            n_effective_regions=0.0,
            n_regions_in_pool=n_in_pool,
            method_histogram=method_histogram,
            fit_status="fallback_bootstrap",
        )

    # ---- Section 6.3.2: robust center ----
    log_center_raw = float(np.sum(weights * log_rho_hat) / w_sum)
    abs_dev = np.abs(log_rho_hat - log_center_raw)
    log_mad = _weighted_median(abs_dev, weights)
    if not np.isfinite(log_mad):
        log_mad = 0.0
    log_dispersion = max(1.4826 * log_mad, log_dispersion_floor)

    lo = log_center_raw - huber_k * log_dispersion
    hi = log_center_raw + huber_k * log_dispersion
    clipped = np.clip(log_rho_hat, lo, hi)
    log_center = float(np.sum(weights * clipped) / w_sum)

    log_prev = float(np.log(prev_rho0))
    log_center_damped = (1.0 - damping) * log_prev + damping * log_center
    rho0_mean_next = float(np.exp(log_center_damped))
    rho0_mean_next = max(rho0_mean_next, rho_floor)

    # ---- Section 6.4: Gamma conjugacy + rescale to match rho0_mean_next ----
    alpha_hat = alpha_floor + float(np.sum(weights * gdna))
    beta_hat = beta_floor + float(np.sum(weights * region_bp))
    alpha_damped = (1.0 - damping) * float(previous.alpha0) + damping * alpha_hat
    beta_damped = (1.0 - damping) * float(previous.beta0) + damping * beta_hat
    if alpha_damped <= 0.0 or beta_damped <= 0.0:
        alpha0 = float(previous.alpha0)
        beta0 = float(previous.beta0)
    else:
        raw_mean = alpha_damped / beta_damped
        if raw_mean <= 0.0 or not np.isfinite(raw_mean):
            alpha0 = alpha_damped
            beta0 = beta_damped
        else:
            scale = rho0_mean_next / raw_mean
            sqrt_scale = float(np.sqrt(scale))
            alpha0 = alpha_damped * sqrt_scale
            beta0 = beta_damped / sqrt_scale

    # ---- Section 6.6: fit_status ----
    log_jump = abs(float(np.log(rho0_mean_next)) - log_prev)
    if previous.fit_status == "fallback_bootstrap" and previous.n_regions_in_pool == 0:
        fit_status = "prior_only"
    elif log_jump < convergence_log_tol:
        fit_status = "converged"
    else:
        fit_status = "iterating"

    return BackgroundDensity(
        rho0_mean=rho0_mean_next,
        alpha0=float(alpha0),
        beta0=float(beta0),
        log_dispersion=float(log_dispersion),
        n_effective_regions=float(w_sum),
        n_regions_in_pool=n_in_pool,
        method_histogram=method_histogram,
        fit_status=fit_status,
    )


def _derive_region_flags(
    p_states: np.ndarray,
    local_posterior: BoundaryLocalPosterior,
    strand_channels: RegionGdnaChannelEstimate | None,
) -> np.ndarray:
    region_count = int(p_states.shape[0])
    flags = np.zeros(region_count, dtype=np.uint16)
    flags[np.max(p_states, axis=1) < 0.6] |= FLAG_STATE_AMBIGUOUS

    p_expressed = p_states[:, STATE_EXPRESSED]
    uncertain_expression = (p_expressed > 0.3) & (p_expressed < 0.7)
    flags[uncertain_expression] |= FLAG_EXPRESSED_UNCERTAIN

    sparse_boundary = (np.asarray(local_posterior.alpha_excess) <= 0.0) | (
        np.asarray(local_posterior.beta_excess) <= 0.0
    )
    flags[sparse_boundary] |= FLAG_BOUNDARY_SPARSE

    if strand_channels is None or not np.isfinite(strand_channels.kappa_d):
        flags |= FLAG_STRAND_UNINFORMATIVE
    return flags


def calibration_e_step(
    region_arrays: RegionArrays,
    observation: DensityObservation,
    boundaries: BoundaryTable,
    background: BackgroundModel,
    local_posterior: BoundaryLocalPosterior | None = None,
    strand_channels: RegionGdnaChannelEstimate | None = None,
    *,
    pass_index: int = 0,
    confidence: float = 0.95,
    background_boost: float = 1.0,
    transfer_weight: np.ndarray | None = None,
    unspliced_counts: np.ndarray,
    background_density: BackgroundDensity,
    previous_region_exposure: RegionExposure | None = None,
    force_zero_gdna_mass: bool = False,
) -> CalibrationStepResult:
    """Run one two-state expression calibration E-step."""
    region_count = _validate_region_inputs(
        region_arrays,
        observation,
        boundaries,
        background,
        strand_channels,
    )
    if local_posterior is None:
        local_posterior = build_boundary_local_posterior(
            observation,
            boundaries,
            background,
            strand_channels=strand_channels,
            confidence=confidence,
        )
    else:
        _validate_local_posterior(local_posterior, region_count)

    sweep = run_boundary_sweep(
        local_posterior,
        boundaries,
        observation,
        background,
        transfer_weight=transfer_weight,
        strand_channels=strand_channels,
        confidence=confidence,
    )
    state_log_prior = build_state_log_prior(
        region_arrays,
        background,
        pass_index=pass_index,
        background_boost=background_boost,
    )
    logbf_expression = build_logbf_expression(observation, strand_channels)
    logbf_strand = build_logbf_strand(strand_channels, region_count)
    log_tensor = build_state_log_tensor(
        state_log_prior,
        logbf_expression,
        logbf_strand,
    )
    p_states = normalize_state_log_tensor(log_tensor)

    contained_leff = _as_float64_vector(
        "observation.contained_leff", observation.contained_leff, region_count
    )
    off_target_mu = np.maximum(float(background.rho_off_mean), 0.0) * np.maximum(
        contained_leff,
        0.0,
    )
    upper_off_target = predict_contained_gdna_from_excess(
        background,
        contained_leff,
        np.zeros(region_count, dtype=np.float64),
        np.zeros(region_count, dtype=np.float64),
        confidence=confidence,
    )[1].astype(np.float64, copy=False)
    mu_gdna = off_target_mu.astype(np.float32)
    upper_gdna = upper_off_target.astype(np.float32)
    upper_gdna = np.maximum(upper_gdna, mu_gdna)

    if strand_channels is None:
        rna_lower = np.zeros(region_count, dtype=np.float32)
    else:
        rna_lower = np.maximum(
            _as_float64_vector(
                "strand_channels.contained_rna_lower",
                strand_channels.contained_rna_lower,
                region_count,
            ),
            0.0,
        ).astype(np.float32)

    prior_mass = build_region_unspliced_mass(
        observation,
        region_size_bp=region_arrays.region_size_bp,
        unspliced_counts=unspliced_counts,
        strand_channels=strand_channels,
        local_posterior=local_posterior,
        sweep=sweep,
        background_density=background_density,
        force_zero_gdna=force_zero_gdna_mass,
    )

    region_exposure = estimate_region_exposure(
        prior_mass,
        background_density,
        p_states[:, STATE_UNEXPRESSED].astype(np.float64),
        previous=previous_region_exposure,
    )
    flags = _derive_region_flags(p_states, local_posterior, strand_channels)
    sum_log_evidence = float(np.sum(logsumexp(log_tensor, axis=1), dtype=np.float64))

    return CalibrationStepResult(
        p_states=p_states,
        mu_gdna=mu_gdna,
        upper_gdna=upper_gdna,
        rna_lower=rna_lower,
        region_unspliced_mass=prior_mass,
        region_exposure=region_exposure,
        flags=flags,
        local_posterior=local_posterior,
        sweep=sweep,
        log_tensor=log_tensor,
        sum_log_evidence=sum_log_evidence,
    )


def calibration_m_step(
    observation: DensityObservation,
    background: BackgroundModel,
    p_states: np.ndarray,
    strand_channels: RegionGdnaChannelEstimate | None = None,
    *,
    damping: float = 0.5,
    alpha_floor: float = 1.0,
    beta_floor: float = 1.0,
) -> tuple[BackgroundModel, float | None]:
    """Refit scalar calibration parameters from state posteriors."""
    if not np.isfinite(damping) or not 0.0 <= float(damping) <= 1.0:
        raise ValueError(f"damping must be finite and in [0, 1]; got {damping!r}.")
    if not np.isfinite(alpha_floor) or float(alpha_floor) <= 0.0:
        raise ValueError(f"alpha_floor must be finite and > 0; got {alpha_floor!r}.")
    if not np.isfinite(beta_floor) or float(beta_floor) <= 0.0:
        raise ValueError(f"beta_floor must be finite and > 0; got {beta_floor!r}.")

    states = np.asarray(p_states, dtype=np.float64)
    if states.ndim != 2 or states.shape[1] != N_STATES:
        raise ValueError(f"p_states must have shape (R, {N_STATES}); got {states.shape}.")
    region_count = int(states.shape[0])
    contained_leff = _as_float64_vector(
        "observation.contained_leff", observation.contained_leff, region_count
    )
    if strand_channels is None:
        gdna_count = _as_float64_vector(
            "observation.contained_count", observation.contained_count, region_count
        )
        kappa_d = None
    else:
        gdna_count = _as_float64_vector(
            "strand_channels.contained_mean", strand_channels.contained_mean, region_count
        )
        kappa_d = float(strand_channels.kappa_d)

    seed_mask = np.asarray(background.seed_mask, dtype=bool)
    if seed_mask.shape != (region_count,):
        raise ValueError(
            "background.seed_mask must match region count; "
            f"got {seed_mask.shape}, expected {(region_count,)}."
        )
    unexpressed_weight = states[:, STATE_UNEXPRESSED] * seed_mask.astype(np.float64)
    weighted_fragments = float(
        np.sum(unexpressed_weight * gdna_count, dtype=np.float64)
    )
    weighted_eff_length = float(
        np.sum(unexpressed_weight * contained_leff, dtype=np.float64)
    )
    alpha_hat = float(alpha_floor) + weighted_fragments
    beta_hat = float(beta_floor) + weighted_eff_length
    eta = float(damping)
    alpha_next = (1.0 - eta) * float(background.rho_off_alpha) + eta * alpha_hat
    beta_next = (1.0 - eta) * float(background.rho_off_beta) + eta * beta_hat
    beta_next = max(beta_next, np.finfo(np.float64).tiny)
    rho_next = alpha_next / beta_next

    n_seed_regions = int(np.count_nonzero((unexpressed_weight > 0.0) & seed_mask))
    if n_seed_regions == 0 or weighted_eff_length <= 0.0:
        fit_status = "prior_only"
    elif background.fit_status == "ok":
        fit_status = "ok"
    else:
        fit_status = "sparse"

    next_background = BackgroundModel(
        rho_off_alpha=float(alpha_next),
        rho_off_beta=float(beta_next),
        rho_off_mean=float(rho_next),
        seed_mask=seed_mask.copy(),
        top_t_exclusion_mask=np.asarray(background.top_t_exclusion_mask, dtype=bool).copy(),
        n_seed_regions=n_seed_regions,
        n_fragments=weighted_fragments,
        eff_length=weighted_eff_length,
        fit_status=fit_status,
        flags=np.asarray(background.flags, dtype=np.uint16).copy(),
    )

    return next_background, kappa_d


def _relative_change(current_value: float, previous_value: float | None) -> float:
    if previous_value is None:
        return float("inf")
    denominator = max(abs(float(previous_value)), _EPS)
    return abs(float(current_value) - float(previous_value)) / denominator


def run_calibration_iteration(
    region_arrays: RegionArrays,
    observation: DensityObservation,
    boundaries: BoundaryTable,
    background: BackgroundModel,
    strand_channels: RegionGdnaChannelEstimate | None = None,
    *,
    local_posterior: BoundaryLocalPosterior | None = None,
    transfer_weight: np.ndarray | None = None,
    max_calibration_passes: int = 5,
    p_tol: float = 0.01,
    rho_tol: float = 0.02,
    damping: float = 0.5,
    confidence: float = 0.95,
    background_boost: float = 1.0,
    unspliced_counts: np.ndarray,
    force_zero_gdna_mass: bool = False,
) -> RegionCalibration:
    """Run the two-state expression calibration loop."""
    if int(max_calibration_passes) < 1:
        raise ValueError(
            "max_calibration_passes must be >= 1; "
            f"got {max_calibration_passes!r}."
        )
    if not np.isfinite(p_tol) or float(p_tol) < 0.0:
        raise ValueError(f"p_tol must be finite and >= 0; got {p_tol!r}.")
    if not np.isfinite(rho_tol) or float(rho_tol) < 0.0:
        raise ValueError(f"rho_tol must be finite and >= 0; got {rho_tol!r}.")

    current_background = background
    current_kappa = float(strand_channels.kappa_d) if strand_channels is not None else None
    # Bootstrap the BackgroundDensity from the incoming BackgroundModel and
    # refit it in lock-step with the BackgroundModel each pass.
    current_density: BackgroundDensity = BackgroundDensity.from_bootstrap(background)
    previous_region_exposure: RegionExposure | None = None
    previous_p_states: np.ndarray | None = None
    previous_rho: float | None = None
    diagnostics: list[dict[str, object]] = []
    last_step: CalibrationStepResult | None = None
    converged = False

    for pass_index in range(int(max_calibration_passes)):
        step = calibration_e_step(
            region_arrays,
            observation,
            boundaries,
            current_background,
            local_posterior=local_posterior,
            strand_channels=strand_channels,
            pass_index=pass_index,
            confidence=confidence,
            background_boost=background_boost,
            transfer_weight=transfer_weight,
            unspliced_counts=unspliced_counts,
            background_density=current_density,
            previous_region_exposure=previous_region_exposure,
            force_zero_gdna_mass=force_zero_gdna_mass,
        )
        if previous_p_states is None:
            max_state_shift = float("inf")
        else:
            max_state_shift = float(np.max(np.abs(step.p_states - previous_p_states)))
        rho_shift = _relative_change(current_background.rho_off_mean, previous_rho)
        if previous_region_exposure is None or previous_region_exposure.tau2_method not in {
            "moment",
            "moment_damped",
        }:
            exposure_tau2_shift = float("inf")
        else:
            exposure_tau2_shift = _relative_change(
                step.region_exposure.tau2,
                previous_region_exposure.tau2,
            )
        converged = bool(
            previous_p_states is not None
            and max_state_shift < float(p_tol)
            and rho_shift < float(rho_tol)
        )
        diagnostics.append(
            {
                "pass_index": int(pass_index),
                "max_state_shift": float(max_state_shift),
                "rho_off": float(current_background.rho_off_mean),
                "relative_rho_shift": float(rho_shift),
                "kappa_d": current_kappa,
                "sum_log_evidence": float(step.sum_log_evidence),
                "exposure_tau2": float(step.region_exposure.tau2),
                "exposure_tau2_hat": float(step.region_exposure.tau2_hat),
                "exposure_relative_tau2_shift": float(exposure_tau2_shift),
                "exposure_tau2_method": str(step.region_exposure.tau2_method),
                "exposure_tau2_pool_size": int(step.region_exposure.tau2_pool_size),
                "converged": bool(converged),
                "n_regions_expressed": int(np.count_nonzero(step.p_expressed > 0.5)),
                "n_regions_unexpressed": int(np.count_nonzero(step.p_unexpressed > 0.5)),
            }
        )
        last_step = step
        if converged or pass_index == int(max_calibration_passes) - 1:
            break

        previous_p_states = step.p_states.copy()
        previous_region_exposure = step.region_exposure
        previous_rho = float(current_background.rho_off_mean)
        current_background, current_kappa = calibration_m_step(
            observation,
            current_background,
            step.p_states,
            strand_channels,
            damping=damping,
        )
        current_density = estimate_background_density(
            step.region_unspliced_mass,
            step.p_unexpressed.astype(np.float64),
            previous=current_density,
            damping=damping,
        )

    if last_step is None:  # pragma: no cover - max_calibration_passes guard prevents this.
        raise RuntimeError("run_calibration_iteration: no calibration pass was executed.")

    return RegionCalibration(
        p_states=last_step.p_states,
        mu_gdna=last_step.mu_gdna,
        upper_gdna=last_step.upper_gdna,
        rna_lower=last_step.rna_lower,
        region_unspliced_mass=last_step.region_unspliced_mass,
        background_density=current_density,
        region_exposure=last_step.region_exposure,
        kappa_d=current_kappa,
        n_passes=len(diagnostics),
        converged=bool(converged),
        flags=last_step.flags,
        pass_diagnostics=tuple(diagnostics),
        background_model=current_background,
        boundary_local=last_step.local_posterior,
        boundary_sweep=last_step.sweep,
    )
