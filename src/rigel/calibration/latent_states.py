"""Two-state expression calibration tensor for the calibration path."""

from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

from ._arrays import RegionArrays
from .background_model import BackgroundModel
from .density_observation import DensityObservation
from .fractional_evidence import is_intergenic, is_intron_only
from .strand_deconv import RegionGdnaChannelEstimate

__all__ = [
    "STATE_UNEXPRESSED",
    "STATE_EXPRESSED",
    "N_STATES",
    "STATE_IS_EXPRESSED",
    "STATE_NAMES",
    "build_state_log_prior",
    "build_logbf_expression",
    "build_logbf_strand",
    "build_state_log_tensor",
    "normalize_state_log_tensor",
    "build_state_tensor",
]

STATE_UNEXPRESSED: int = 0
STATE_EXPRESSED: int = 1
N_STATES: int = 2

STATE_IS_EXPRESSED: np.ndarray = np.array([False, True], dtype=bool)
STATE_NAMES: tuple[str, ...] = ("unexpressed", "expressed")

_MIN_LOG: float = -20.0


def _validate_region_vector(name: str, values: np.ndarray, region_count: int) -> np.ndarray:
    array = np.asarray(values)
    if array.shape != (region_count,):
        raise ValueError(f"{name} must have shape ({region_count},); got {array.shape}.")
    return array


def _observation_region_count(observation: DensityObservation) -> int:
    return int(np.asarray(observation.contained_count).shape[0])


def _validate_background(background: BackgroundModel, region_count: int) -> None:
    seed_mask = np.asarray(background.seed_mask)
    if seed_mask.shape != (region_count,):
        raise ValueError(
            "background.seed_mask must match region count; "
            f"got {seed_mask.shape}, expected {(region_count,)}."
        )
    if not np.isfinite(background.rho_off_mean) or float(background.rho_off_mean) < 0.0:
        raise ValueError(
            "background.rho_off_mean must be finite and non-negative; "
            f"got {background.rho_off_mean!r}."
        )


def build_state_log_prior(
    region_arrays: RegionArrays,
    background: BackgroundModel,
    *,
    pass_index: int = 0,
    background_boost: float = 1.0,
) -> np.ndarray:
    """Return annotation-derived expression-state log priors with an early background boost."""
    signatures = np.asarray(region_arrays.signature)
    region_count = int(signatures.shape[0])
    _validate_background(background, region_count)
    if int(pass_index) < 0:
        raise ValueError(f"pass_index must be >= 0; got {pass_index!r}.")
    boost = float(background_boost)
    if not np.isfinite(boost) or boost < 0.0:
        raise ValueError(f"background_boost must be finite and >= 0; got {background_boost!r}.")

    state_log_prior = np.zeros((region_count, N_STATES), dtype=np.float64)

    intergenic = is_intergenic(signatures)
    intron_only = is_intron_only(signatures)

    state_log_prior[intergenic, STATE_UNEXPRESSED] = 1.0
    state_log_prior[intergenic, STATE_EXPRESSED] = -2.0

    state_log_prior[intron_only, STATE_UNEXPRESSED] = 0.75
    state_log_prior[intron_only, STATE_EXPRESSED] = -1.0

    if int(pass_index) < 2 and boost > 0.0:
        seed_mask = np.asarray(background.seed_mask, dtype=bool)
        state_log_prior[seed_mask, STATE_UNEXPRESSED] += boost

    return state_log_prior


def build_logbf_expression(
    observation: DensityObservation,
    strand_channels: RegionGdnaChannelEstimate | None = None,
    *,
    spliced_log_penalty: float = -10.0,
    rna_lower_log_penalty: float = -5.0,
    rna_lower_threshold: float = 0.1,
) -> np.ndarray:
    """Return expression-vs-not-expression log likelihood factors."""
    region_count = _observation_region_count(observation)
    spliced = _validate_region_vector("observation.spliced_count", observation.spliced_count, region_count)
    if not np.isfinite(spliced_log_penalty) or float(spliced_log_penalty) > 0.0:
        raise ValueError(
            "spliced_log_penalty must be finite and <= 0; "
            f"got {spliced_log_penalty!r}."
        )
    if not np.isfinite(rna_lower_log_penalty) or float(rna_lower_log_penalty) > 0.0:
        raise ValueError(
            "rna_lower_log_penalty must be finite and <= 0; "
            f"got {rna_lower_log_penalty!r}."
        )
    if not np.isfinite(rna_lower_threshold) or float(rna_lower_threshold) < 0.0:
        raise ValueError(
            "rna_lower_threshold must be finite and >= 0; "
            f"got {rna_lower_threshold!r}."
        )

    logbf = np.zeros((region_count, 2), dtype=np.float64)
    logbf[np.asarray(spliced, dtype=np.float64) > 0.0, 0] += float(spliced_log_penalty)

    if strand_channels is not None:
        rna_lower = _validate_region_vector(
            "strand_channels.contained_rna_lower",
            strand_channels.contained_rna_lower,
            region_count,
        )
        logbf[np.asarray(rna_lower, dtype=np.float64) > float(rna_lower_threshold), 0] += float(
            rna_lower_log_penalty
        )

    return logbf


def build_logbf_strand(
    strand_channels: RegionGdnaChannelEstimate | None,
    region_count: int,
    *,
    max_abs_logbf: float = 3.0,
) -> np.ndarray:
    """Return a conservative state-specific strand summary factor.

    The exact small-count strand tensor needs raw folded counts, which are not
    carried by `RegionGdnaChannelEstimate`. This summary term uses RNA lower
    bounds conservatively until the exact wrapper is wired in.
    """
    count = int(region_count)
    if count < 0:
        raise ValueError(f"region_count must be >= 0; got {region_count!r}.")
    logbf = np.zeros((count, N_STATES), dtype=np.float64)
    if strand_channels is None:
        return logbf
    if not np.isfinite(max_abs_logbf) or float(max_abs_logbf) < 0.0:
        raise ValueError(f"max_abs_logbf must be finite and >= 0; got {max_abs_logbf!r}.")

    rna_lower = _validate_region_vector(
        "strand_channels.contained_rna_lower", strand_channels.contained_rna_lower, count
    ).astype(np.float64, copy=False)
    precision = _validate_region_vector(
        "strand_channels.contained_precision", strand_channels.contained_precision, count
    ).astype(np.float64, copy=False)
    precision = np.clip(np.nan_to_num(precision, nan=0.0, posinf=1.0, neginf=0.0), 0.0, 1.0)
    evidence = np.log1p(np.maximum(rna_lower, 0.0)) * (0.5 + 0.5 * precision)
    evidence = np.minimum(evidence, float(max_abs_logbf))
    logbf[:, STATE_UNEXPRESSED] -= evidence
    return logbf


def build_state_log_tensor(
    state_log_prior: np.ndarray,
    logbf_expression: np.ndarray,
    logbf_strand: np.ndarray,
) -> np.ndarray:
    """Assemble the unnormalized two-state expression log tensor."""
    prior = np.asarray(state_log_prior, dtype=np.float64)
    if prior.ndim != 2 or prior.shape[1] != N_STATES:
        raise ValueError(f"state_log_prior must have shape (R, {N_STATES}); got {prior.shape}.")
    region_count = int(prior.shape[0])
    expression = np.asarray(logbf_expression, dtype=np.float64)
    strand = np.asarray(logbf_strand, dtype=np.float64)
    expected_state_shape = (region_count, N_STATES)
    if expression.shape != expected_state_shape:
        raise ValueError(f"logbf_expression must have shape {expected_state_shape}; got {expression.shape}.")
    if strand.shape != expected_state_shape:
        raise ValueError(f"logbf_strand must have shape {expected_state_shape}; got {strand.shape}.")

    log_tensor = prior.copy() + expression
    log_tensor += strand
    return np.nan_to_num(log_tensor, nan=_MIN_LOG, posinf=0.0, neginf=_MIN_LOG)


def normalize_state_log_tensor(log_tensor: np.ndarray) -> np.ndarray:
    """Normalize a state log tensor to row-stochastic state probabilities."""
    tensor = np.asarray(log_tensor, dtype=np.float64)
    if tensor.ndim != 2 or tensor.shape[1] != N_STATES:
        raise ValueError(f"log_tensor must have shape (R, {N_STATES}); got {tensor.shape}.")
    log_norm = logsumexp(tensor, axis=1, keepdims=True)
    probabilities = np.exp(tensor - log_norm)
    probabilities = np.nan_to_num(probabilities, nan=0.0, posinf=1.0, neginf=0.0)
    row_sum = probabilities.sum(axis=1, keepdims=True)
    probabilities = np.divide(
        probabilities,
        row_sum,
        out=np.full_like(probabilities, 1.0 / float(N_STATES)),
        where=row_sum > 0.0,
    )
    return probabilities.astype(np.float32)


def build_state_tensor(
    state_log_prior: np.ndarray,
    logbf_expression: np.ndarray,
    logbf_strand: np.ndarray,
) -> np.ndarray:
    """Assemble and normalize the two-state posterior tensor."""
    log_tensor = build_state_log_tensor(
        state_log_prior,
        logbf_expression,
        logbf_strand,
    )
    return normalize_state_log_tensor(log_tensor)
