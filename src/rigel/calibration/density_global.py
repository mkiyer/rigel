"""Global gDNA density estimates for the calibration path."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from ..frag_length_model import FragmentLengthModel
from ._arrays import PayloadArrays, RegionArrays
from .fractional_evidence import is_intergenic, is_intron_only
from .strand_summary import STRAND_CONTRAST_NUMERICAL_FLOOR


__all__ = [
    "DensityType",
    "GlobalGdnaDensity",
    "GlobalDensityTable",
    "StrandBalanceEstimate",
    "compute_global_densities",
    "estimate_strand_balance",
    "l_eff_contained",
    "STRAND_CONTRAST_NUMERICAL_FLOOR",
]


DensityType = Literal["INTERGENIC", "INTRON", "EXON-INTRON", "EXON-CONTAINED"]

STRAND_KAPPA_DEFAULT: float = 1.0e6
STRAND_KAPPA_MIN: float = 1.0e-3
STRAND_KAPPA_MAX: float = 1.0e6
MIN_REGIONS_FOR_STRAND_MOM: int = 2


def l_eff_contained(
    spans_bp: np.ndarray,
    gdna_fl: FragmentLengthModel,
) -> np.ndarray:
    """FL-PMF-weighted contained effective length, in bp."""
    return gdna_fl.compute_all_transcript_eff_lens(
        np.asarray(spans_bp, dtype=np.int64),
        min_value=0.0,
    )


@dataclass(frozen=True, slots=True)
class GlobalGdnaDensity:
    """One density estimate for one global calibration category."""

    type: DensityType
    rho: float
    n_fragments: float
    eff_length_bp: float
    n_regions_used: int
    n_fragments_estimated: float | None = None
    n_rows_eligible: int | None = None
    strand_active: bool = False
    rho_uncorrected: float | None = None
    strand_correction_fragments: float = 0.0
    n_strand_informative_regions: int = 0
    strand_informative_exposure_fraction: float = 0.0
    n_uninf_observed: int = 0

    def __post_init__(self) -> None:
        if self.n_fragments_estimated is None:
            object.__setattr__(self, "n_fragments_estimated", float(self.n_fragments))
        if self.n_rows_eligible is None:
            object.__setattr__(self, "n_rows_eligible", int(self.n_regions_used))
        if self.rho_uncorrected is None:
            object.__setattr__(self, "rho_uncorrected", float(self.rho))

    def to_summary_dict(self) -> dict[str, object]:
        return {
            "type": self.type,
            "rho": float(self.rho),
            "n_fragments": float(self.n_fragments),
            "eff_length_bp": float(self.eff_length_bp),
            "n_regions_used": int(self.n_regions_used),
            "n_fragments_estimated": float(self.n_fragments_estimated),
            "n_rows_eligible": int(self.n_rows_eligible),
            "rho_uncorrected": float(self.rho_uncorrected),
        }


@dataclass(frozen=True, slots=True)
class StrandBalanceEstimate:
    """Symmetric beta-binomial strand overdispersion estimate for gDNA."""

    kappa: float
    n_regions: int
    n_fragments: float
    n_pos: float
    n_neg: float
    observed_pos_fraction: float
    residual_sum: float
    binomial_variance_sum: float
    max_overdispersed_variance_sum: float
    overdispersion_factor: float
    fallback_used: bool = False
    fallback_reason: str = ""

    @property
    def alpha(self) -> float:
        return float(self.kappa / 2.0)

    @property
    def beta(self) -> float:
        return float(self.kappa / 2.0)

    def to_summary_dict(self) -> dict[str, object]:
        return {
            "model": "symmetric_beta_binomial",
            "mean_fixed": 0.5,
            "kappa": float(self.kappa),
            "alpha": self.alpha,
            "beta": self.beta,
            "n_regions": int(self.n_regions),
            "n_fragments": float(self.n_fragments),
            "n_pos": float(self.n_pos),
            "n_neg": float(self.n_neg),
            "observed_pos_fraction": float(self.observed_pos_fraction),
            "residual_sum": float(self.residual_sum),
            "binomial_variance_sum": float(self.binomial_variance_sum),
            "max_overdispersed_variance_sum": float(self.max_overdispersed_variance_sum),
            "overdispersion_factor": float(self.overdispersion_factor),
            "fallback_used": bool(self.fallback_used),
            "fallback_reason": str(self.fallback_reason),
        }


@dataclass(frozen=True, slots=True)
class GlobalDensityTable:
    """The complete global-density block of the calibration result."""

    intergenic: GlobalGdnaDensity
    intron: GlobalGdnaDensity
    exon_intron: GlobalGdnaDensity
    gdna_fl: FragmentLengthModel
    exon_contained: GlobalGdnaDensity | None = None
    strand_balance: StrandBalanceEstimate | None = None

    @property
    def gdna_fl_mean(self) -> float:
        return float(self.gdna_fl.mean)

    @property
    def contained_n_fragments(self) -> float:
        return float(self.intergenic.n_fragments + self.intron.n_fragments)

    @property
    def contained_eff_length_bp(self) -> float:
        return float(self.intergenic.eff_length_bp + self.intron.eff_length_bp)

    @property
    def contained_rho(self) -> float:
        denom = self.contained_eff_length_bp
        return float(self.contained_n_fragments / denom) if denom > 0.0 else 0.0

    def to_summary_dict(self) -> dict[str, object]:
        return {
            "gdna_fl_mean": self.gdna_fl_mean,
            "contained_global": {
                "rho": self.contained_rho,
                "n_fragments": self.contained_n_fragments,
                "eff_length_bp": self.contained_eff_length_bp,
                "n_regions_used": int(self.intergenic.n_regions_used + self.intron.n_regions_used),
            },
            "intergenic": self.intergenic.to_summary_dict(),
            "intron": self.intron.to_summary_dict(),
            "exon_intron": self.exon_intron.to_summary_dict(),
            "strand_balance": (
                self.strand_balance.to_summary_dict() if self.strand_balance is not None else None
            ),
        }


def compute_global_densities(
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    *,
    gdna_fl: FragmentLengthModel,
) -> GlobalDensityTable:
    """Estimate initial global gDNA densities from contained pure regions."""
    intergenic_mask = is_intergenic(region_arrays.signature)
    intron_mask = is_intron_only(region_arrays.signature)

    intergenic = _contained_density_for_mask(
        "INTERGENIC",
        intergenic_mask,
        region_arrays,
        payload_arrays,
        gdna_fl,
    )
    intron = _contained_density_for_mask(
        "INTRON",
        intron_mask,
        region_arrays,
        payload_arrays,
        gdna_fl,
    )

    return GlobalDensityTable(
        intergenic=intergenic,
        intron=intron,
        exon_intron=_empty_density("EXON-INTRON"),
        gdna_fl=gdna_fl,
        strand_balance=estimate_strand_balance(
            payload_arrays.contained_unspliced_pos,
            payload_arrays.contained_unspliced_neg,
            intergenic_mask | intron_mask,
        ),
    )


def estimate_strand_balance(
    pos_counts: np.ndarray,
    neg_counts: np.ndarray,
    region_mask: np.ndarray,
) -> StrandBalanceEstimate:
    """Estimate beta-binomial strand overdispersion with fixed mean 0.5."""
    pos = np.asarray(pos_counts, dtype=np.float64)
    neg = np.asarray(neg_counts, dtype=np.float64)
    if pos.shape != neg.shape:
        raise ValueError(
            f"estimate_strand_balance: pos.shape ({pos.shape}) != neg.shape ({neg.shape})."
        )

    mask = np.asarray(region_mask, dtype=bool)
    if mask.shape != pos.shape:
        raise ValueError(
            f"estimate_strand_balance: region_mask.shape ({mask.shape}) != pos.shape ({pos.shape})."
        )

    total = pos + neg
    eligible = mask & (total > 0.0)
    n_regions = int(eligible.sum())
    if n_regions < MIN_REGIONS_FOR_STRAND_MOM:
        return _strand_fallback(
            pos[eligible],
            neg[eligible],
            n_regions,
            "n_regions < MIN_REGIONS_FOR_STRAND_MOM",
        )

    pos_e = pos[eligible]
    neg_e = neg[eligible]
    total_e = total[eligible]
    residual = pos_e - 0.5 * total_e
    residual_sum = float(np.sum(residual * residual))
    binomial_variance_sum = float(0.25 * np.sum(total_e))
    max_variance_sum = float(0.25 * np.sum(total_e * total_e))
    n_fragments = float(total_e.sum())
    n_pos = float(pos_e.sum())
    n_neg = float(neg_e.sum())
    observed_pos_fraction = float(n_pos / n_fragments) if n_fragments > 0.0 else 0.5
    overdispersion_factor = (
        float(residual_sum / binomial_variance_sum) if binomial_variance_sum > 0.0 else 0.0
    )

    if binomial_variance_sum <= 0.0:
        return _strand_fallback(pos_e, neg_e, n_regions, "no positive strand exposure")
    if residual_sum <= binomial_variance_sum:
        return StrandBalanceEstimate(
            kappa=STRAND_KAPPA_MAX,
            n_regions=n_regions,
            n_fragments=n_fragments,
            n_pos=n_pos,
            n_neg=n_neg,
            observed_pos_fraction=observed_pos_fraction,
            residual_sum=residual_sum,
            binomial_variance_sum=binomial_variance_sum,
            max_overdispersed_variance_sum=max_variance_sum,
            overdispersion_factor=overdispersion_factor,
            fallback_used=True,
            fallback_reason="residual variance <= binomial expectation",
        )

    numerator = max_variance_sum - residual_sum
    denominator = residual_sum - binomial_variance_sum
    if denominator <= 0.0 or not np.isfinite(denominator):
        kappa = STRAND_KAPPA_MAX
        fallback_used = True
        fallback_reason = "invalid MoM denominator"
    else:
        kappa = numerator / denominator
        fallback_used = False
        fallback_reason = ""
    if not np.isfinite(kappa):
        kappa = STRAND_KAPPA_MAX
        fallback_used = True
        fallback_reason = "non-finite MoM estimate"
    kappa = float(np.clip(kappa, STRAND_KAPPA_MIN, STRAND_KAPPA_MAX))

    return StrandBalanceEstimate(
        kappa=kappa,
        n_regions=n_regions,
        n_fragments=n_fragments,
        n_pos=n_pos,
        n_neg=n_neg,
        observed_pos_fraction=observed_pos_fraction,
        residual_sum=residual_sum,
        binomial_variance_sum=binomial_variance_sum,
        max_overdispersed_variance_sum=max_variance_sum,
        overdispersion_factor=overdispersion_factor,
        fallback_used=fallback_used,
        fallback_reason=fallback_reason,
    )


def _contained_density_for_mask(
    density_type: DensityType,
    region_mask: np.ndarray,
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    gdna_fl: FragmentLengthModel,
) -> GlobalGdnaDensity:
    spans = region_arrays.end - region_arrays.start
    candidate = np.asarray(region_mask, dtype=bool)
    eff_all = l_eff_contained(spans[candidate], gdna_fl)
    counts_all = np.asarray(
        payload_arrays.contained_unspliced_total[candidate],
        dtype=np.float64,
    )
    eligible = eff_all > 0.0
    if not np.any(eligible):
        return _empty_density(density_type, n_rows_eligible=int(candidate.sum()))

    counts = counts_all[eligible]
    eff = eff_all[eligible]
    n_fragments = float(counts.sum())
    eff_total = float(eff.sum())
    rho = float(n_fragments / eff_total) if eff_total > 0.0 else 0.0
    return GlobalGdnaDensity(
        type=density_type,
        rho=rho,
        n_fragments=n_fragments,
        eff_length_bp=eff_total,
        n_regions_used=int(eligible.sum()),
        n_rows_eligible=int(candidate.sum()),
    )


def _empty_density(
    density_type: DensityType,
    *,
    n_rows_eligible: int = 0,
) -> GlobalGdnaDensity:
    return GlobalGdnaDensity(
        type=density_type,
        rho=0.0,
        n_fragments=0.0,
        eff_length_bp=0.0,
        n_regions_used=0,
        n_rows_eligible=n_rows_eligible,
    )


def _strand_fallback(
    pos: np.ndarray,
    neg: np.ndarray,
    n_regions: int,
    reason: str,
) -> StrandBalanceEstimate:
    n_pos = float(np.asarray(pos, dtype=np.float64).sum())
    n_neg = float(np.asarray(neg, dtype=np.float64).sum())
    n_fragments = n_pos + n_neg
    return StrandBalanceEstimate(
        kappa=STRAND_KAPPA_DEFAULT,
        n_regions=int(n_regions),
        n_fragments=n_fragments,
        n_pos=n_pos,
        n_neg=n_neg,
        observed_pos_fraction=float(n_pos / n_fragments) if n_fragments > 0.0 else 0.5,
        residual_sum=0.0,
        binomial_variance_sum=0.0,
        max_overdispersed_variance_sum=0.0,
        overdispersion_factor=0.0,
        fallback_used=True,
        fallback_reason=reason,
    )
