"""Off-target gDNA background model for the v6 calibration path."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .density_observation import DensityObservation
from .strand_deconv import RegionGdnaChannelEstimate


__all__ = [
    "BACKGROUND_MIN_EFF_LENGTH",
    "FLAG_BACKGROUND_CANDIDATE",
    "FLAG_BACKGROUND_SEED",
    "FLAG_TOP_T_EXCLUDED",
    "FLAG_STRAND_RNA_EXCLUDED",
    "FLAG_LOW_OPPORTUNITY",
    "BackgroundModel",
    "fit_background_model",
]


BACKGROUND_MIN_EFF_LENGTH: float = 1.0

FLAG_BACKGROUND_CANDIDATE: int = 1 << 0
FLAG_BACKGROUND_SEED: int = 1 << 1
FLAG_TOP_T_EXCLUDED: int = 1 << 2
FLAG_STRAND_RNA_EXCLUDED: int = 1 << 3
FLAG_LOW_OPPORTUNITY: int = 1 << 4


@dataclass(frozen=True, slots=True)
class BackgroundModel:
    """Gamma posterior for the off-target gDNA density ``rho_off``."""

    rho_off_alpha: float
    rho_off_beta: float
    rho_off_mean: float
    seed_mask: np.ndarray  # bool[R]
    top_t_exclusion_mask: np.ndarray  # bool[R]
    n_seed_regions: int
    n_fragments: float
    eff_length: float
    fit_status: str  # "ok" | "sparse" | "prior_only"
    flags: np.ndarray  # uint16[R]


def _validate_fraction(name: str, value: float) -> float:
    out = float(value)
    if not np.isfinite(out) or out < 0.0 or out >= 1.0:
        raise ValueError(f"{name} must be finite and in [0, 1); got {value!r}.")
    return out


def _top_t_mask(density: np.ndarray, candidate: np.ndarray, fraction: float) -> np.ndarray:
    mask = np.zeros(candidate.shape, dtype=bool)
    idx = np.where(candidate)[0]
    if idx.size == 0 or fraction <= 0.0:
        return mask
    n_exclude = int(np.ceil(float(idx.size) * fraction))
    if n_exclude <= 0:
        return mask
    n_exclude = min(n_exclude, int(idx.size))
    local_density = density[idx]
    order = np.argpartition(local_density, idx.size - n_exclude)[idx.size - n_exclude :]
    mask[idx[order]] = True
    return mask


def fit_background_model(
    observation: DensityObservation,
    strand_channels: RegionGdnaChannelEstimate | None = None,
    *,
    top_t_fraction: float = 0.01,
    min_seed_regions: int = 20,
    alpha_floor: float = 1.0,
    beta_floor: float = 1.0,
    min_eff_length: float = BACKGROUND_MIN_EFF_LENGTH,
) -> BackgroundModel:
    """Fit the off-target gDNA density from background-like regions.

    The model is intentionally simple for Phase I: a weak Gamma prior updated
    by contained gDNA-compatible counts over background seed regions.  When
    strand deconvolution is available, the count is the contained gDNA mean and
    any region with contained or boundary RNA lower-bound evidence is excluded.
    """
    fraction = _validate_fraction("top_t_fraction", top_t_fraction)
    if int(min_seed_regions) < 0:
        raise ValueError(f"min_seed_regions must be >= 0; got {min_seed_regions!r}.")
    if not np.isfinite(alpha_floor) or float(alpha_floor) <= 0.0:
        raise ValueError(f"alpha_floor must be finite and > 0; got {alpha_floor!r}.")
    if not np.isfinite(beta_floor) or float(beta_floor) <= 0.0:
        raise ValueError(f"beta_floor must be finite and > 0; got {beta_floor!r}.")
    if not np.isfinite(min_eff_length) or float(min_eff_length) < 0.0:
        raise ValueError(f"min_eff_length must be finite and >= 0; got {min_eff_length!r}.")

    contained_leff = np.asarray(observation.contained_leff, dtype=np.float64)
    region_count = int(contained_leff.shape[0])
    if strand_channels is None:
        gdna_count = np.asarray(observation.contained_count, dtype=np.float64)
        strand_rna_evidence = np.zeros(region_count, dtype=bool)
    else:
        gdna_count = np.asarray(strand_channels.contained_mean, dtype=np.float64)
        if gdna_count.shape != (region_count,):
            raise ValueError(
                "fit_background_model: strand_channels.contained_mean has shape "
                f"{gdna_count.shape}; expected {(region_count,)}."
            )
        strand_rna_evidence = (
            (np.asarray(strand_channels.contained_rna_lower) > 0.0)
            | (np.asarray(strand_channels.boundary_left_rna_lower) > 0.0)
            | (np.asarray(strand_channels.boundary_right_rna_lower) > 0.0)
        )

    low_opportunity = contained_leff <= float(min_eff_length)
    candidate = (
        np.asarray(observation.is_anchor, dtype=bool)
        & (np.asarray(observation.spliced_count) <= 0.0)
        & ~low_opportunity
    )
    density = np.divide(
        gdna_count,
        np.maximum(contained_leff, 1.0e-12),
        out=np.zeros(region_count, dtype=np.float64),
        where=contained_leff > 0.0,
    )
    top_t_exclusion = _top_t_mask(density, candidate, fraction)
    seed_mask = candidate & ~top_t_exclusion & ~strand_rna_evidence

    flags = np.zeros(region_count, dtype=np.uint16)
    flags[candidate] |= FLAG_BACKGROUND_CANDIDATE
    flags[top_t_exclusion] |= FLAG_TOP_T_EXCLUDED
    flags[strand_rna_evidence] |= FLAG_STRAND_RNA_EXCLUDED
    flags[low_opportunity] |= FLAG_LOW_OPPORTUNITY
    flags[seed_mask] |= FLAG_BACKGROUND_SEED

    n_seed_regions = int(np.count_nonzero(seed_mask))
    n_fragments = float(np.sum(gdna_count[seed_mask], dtype=np.float64))
    eff_length = float(np.sum(contained_leff[seed_mask], dtype=np.float64))
    alpha = float(alpha_floor) + n_fragments
    beta = float(beta_floor) + eff_length
    rho = alpha / beta

    if n_seed_regions == 0 or eff_length <= 0.0:
        fit_status = "prior_only"
    elif n_seed_regions < int(min_seed_regions):
        fit_status = "sparse"
    else:
        fit_status = "ok"

    return BackgroundModel(
        rho_off_alpha=float(alpha),
        rho_off_beta=float(beta),
        rho_off_mean=float(rho),
        seed_mask=seed_mask.astype(bool, copy=False),
        top_t_exclusion_mask=top_t_exclusion.astype(bool, copy=False),
        n_seed_regions=n_seed_regions,
        n_fragments=n_fragments,
        eff_length=eff_length,
        fit_status=fit_status,
        flags=flags,
    )