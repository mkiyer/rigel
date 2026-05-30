"""CalibrationResult — the calibration output schema (doc 04 §5).

Per-region deconvoluted mass (G2 contained, G3 boundary), the per-region
exposure posterior (G4), the five fitted library hyperparameters, and
convergence diagnostics. ``__post_init__`` enforces the *intrinsic*
invariants — those checkable from the result alone (shapes, non-negativity,
parameter ranges, monotone mass-change). Mass conservation against the raw
fragment counts is verified by the calibrator / tests, since it requires the
substrate (which the result does not carry).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import CalibrationConfig
from .errors import CalibrationConvergenceError

# Tolerance for the "mass-change is non-increasing" EM-monotonicity check.
_MONOTONE_ATOL = 1.0e-9


def _check_region_array(arr: np.ndarray, name: str, n_regions: int, *, strictly_positive: bool):
    """Validate a per-region float64 array: shape, dtype, finite, sign."""
    if not isinstance(arr, np.ndarray):
        raise ValueError(f"CalibrationResult.{name} must be a numpy array.")
    if arr.dtype != np.float64:
        raise ValueError(f"CalibrationResult.{name} must be float64; got {arr.dtype}.")
    if arr.shape != (n_regions,):
        raise ValueError(
            f"CalibrationResult.{name} has shape {arr.shape}; expected ({n_regions},)."
        )
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"CalibrationResult.{name} contains non-finite values.")
    if strictly_positive:
        if np.any(arr <= 0.0):
            raise ValueError(f"CalibrationResult.{name} must be strictly positive.")
    elif np.any(arr < 0.0):
        raise ValueError(f"CalibrationResult.{name} must be non-negative.")


@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """Per-region deconvoluted mass + exposure posteriors + library hyperparameters.

    Each output traces to one of G2–G4.
    """

    # --- G2: contained-fragment mass (float64[R]) ---
    mass_g_contained: np.ndarray
    mass_d_contained: np.ndarray

    # --- G3: boundary mass attributed to each region (float64[R]) ---
    mass_g_left: np.ndarray
    mass_d_left: np.ndarray
    mass_g_right: np.ndarray
    mass_d_right: np.ndarray

    # --- G4: per-region exposure posterior (float64[R]) ---
    omega: np.ndarray  # E[omega | data]   (Gamma posterior mean)
    log_omega_var: np.ndarray  # Var(log omega | data)  (delta method)

    # --- library hyperparameters ---
    rho_0: float  # > 0
    phi: float  # > 0  (NB dispersion = Gamma exposure prior shape)
    rho_d_bb: float  # in (0, 1), gDNA strand BB dispersion (kappa_d = 0.5 fixed)
    kappa_rna: float  # in (0, 1), RNA sense mean (from StrandModel; PR 3)
    rho_r_bb: float  # in (0, 1), RNA strand BB dispersion (fit by PR 3)
    eps_s: float  # in (0, 1), gDNA splice-artifact rate

    # --- convergence diagnostics ---
    n_iterations: int
    converged: bool
    mass_change_history: np.ndarray  # float64[n_iterations]

    # --- provenance ---
    n_regions: int
    config: CalibrationConfig

    def __post_init__(self) -> None:
        n = self.n_regions
        if n < 0:
            raise ValueError(f"CalibrationResult.n_regions must be >= 0; got {n}.")

        for name in (
            "mass_g_contained",
            "mass_d_contained",
            "mass_g_left",
            "mass_d_left",
            "mass_g_right",
            "mass_d_right",
        ):
            _check_region_array(getattr(self, name), name, n, strictly_positive=False)
        for name in ("omega", "log_omega_var"):
            _check_region_array(getattr(self, name), name, n, strictly_positive=True)

        if not self.rho_0 > 0.0:
            raise ValueError(f"CalibrationResult.rho_0 must be > 0; got {self.rho_0}.")
        if not self.phi > 0.0:
            raise ValueError(f"CalibrationResult.phi must be > 0; got {self.phi}.")
        for name in ("rho_d_bb", "kappa_rna", "rho_r_bb", "eps_s"):
            value = float(getattr(self, name))
            if not 0.0 < value < 1.0:
                raise ValueError(f"CalibrationResult.{name} must be in (0, 1); got {value}.")

        if self.n_iterations < 0:
            raise ValueError(
                f"CalibrationResult.n_iterations must be >= 0; got {self.n_iterations}."
            )
        if not isinstance(self.converged, (bool, np.bool_)):
            raise ValueError("CalibrationResult.converged must be a bool.")

        hist = self.mass_change_history
        if not isinstance(hist, np.ndarray) or hist.dtype != np.float64:
            raise ValueError("CalibrationResult.mass_change_history must be a float64 array.")
        if hist.shape != (self.n_iterations,):
            raise ValueError(
                f"CalibrationResult.mass_change_history has shape {hist.shape}; "
                f"expected ({self.n_iterations},)."
            )
        if hist.size > 1:
            # EM theory: the mass-change diagnostic must not increase between
            # iterations. A real increase indicates an implementation bug.
            increases = np.diff(hist)
            atol = _MONOTONE_ATOL * (1.0 + np.abs(hist[:-1]))
            if np.any(increases > atol):
                bad = int(np.flatnonzero(increases > atol)[0])
                raise CalibrationConvergenceError(
                    "CalibrationResult.mass_change_history increased between iterations "
                    f"{bad} and {bad + 1} ({hist[bad]!r} -> {hist[bad + 1]!r}); EM "
                    "monotonicity violated."
                )


__all__ = ["CalibrationResult"]
