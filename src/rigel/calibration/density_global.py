"""rigel.calibration.density_global — Global gDNA density estimates.

Reduces the M3 ``CalibrationScanPayload`` plus the M2 region table to
exactly three scalar gDNA densities (fragments / bp), one for each
fragment-mask category that produces a usable gDNA estimate:

* ``INTERGENIC``   — every-block-intergenic fragments.
* ``INTRON``       — every-block-intronic fragments.
* ``EXON-INTRON``  — capture-aware boundary flux on internal exon edges.

Each density is paired with an NB overdispersion estimate (:class:`KappaEstimate`)
that the M6 locoregional EB shrinkage consumes.

See ``docs/calibration/calibration_v6_plan.md`` §2.6 and
``docs/calibration/m4_implementation_plan.md`` §3.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from ..frag_length_model import FragmentLengthModel
from ._exposure import boundary_crossing_exposure
from ._kappa import KappaEstimate, estimate_kappa
from ._orient import ORIENT_OPP, ORIENT_SAME, ORIENT_UNINF, StrandSummary
from .regions import RegionStrand, RegionType
from .scan_payload import (
    MASK_INTERGENIC as _MASK_INTERGENIC,
    CalibrationScanPayload,
)


# (mask constants are imported from scan_payload — single source of truth)

DensityType = Literal["INTERGENIC", "INTRON", "EXON-INTRON"]

_ERR_SUFFIX = " Rebuild the index or rerun the scan."
STRAND_CONTRAST_NUMERICAL_FLOOR = 1.0e-3


def l_eff_contained(
    spans_bp: np.ndarray,
    gdna_fl: FragmentLengthModel,
) -> np.ndarray:
    """FL-PMF-weighted containment effective length.

    For each region span :math:`L`, returns
    :math:`\\sum_{\\ell=1}^{L} h(\\ell)\\,(L - \\ell + 1)`,
    the expected number of distinct fragment start positions whose
    fragment is fully contained in a region of that span, under the
    gDNA fragment-length distribution :math:`h`.

    This is the MLE denominator for the Poisson rate :math:`\\rho`
    that pairs with the C++ scanner's *contained* numerator (every
    aligned block in mask-bit positions). Reuses the salmon-style
    eCDF cache already implemented in :class:`FragmentLengthModel` for
    transcript effective lengths.
    """
    return gdna_fl.compute_all_transcript_eff_lens(
        np.asarray(spans_bp, dtype=np.int64),
        min_value=0.0,
    )


@dataclass(frozen=True, slots=True)
class GlobalGdnaDensity:
    """One density estimate for one mask category."""

    type: DensityType
    rho: float  # fragments / bp
    n_fragments: int  # raw observed numerator
    eff_length_bp: float  # denominator used for rho
    n_regions_used: int  # rows that contributed to the denominator
    kappa: KappaEstimate  # NB overdispersion + fallback diagnostic
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


@dataclass(frozen=True, slots=True)
class GlobalDensityTable:
    """The complete global-density block of the calibration result.

    Carries the finalized gDNA FL model that was used to compute the
    densities; downstream consumers (per-locus prior, full-genome gDNA
    projection) reuse it for FL-PMF-weighted effective lengths so the
    numerator (contained fragments) and denominator (FL-weighted span)
    measure the same geometric event.
    """

    intergenic: GlobalGdnaDensity
    intron: GlobalGdnaDensity
    exon_intron: GlobalGdnaDensity
    gdna_fl: FragmentLengthModel  # the model used to weight L_eff

    @property
    def gdna_fl_mean(self) -> float:
        """Mean of the gDNA FL distribution; used for capture-window geometry
        in the EXON-INTRON branch and as a summary statistic."""
        return float(self.gdna_fl.mean)

    def to_summary_dict(self) -> dict[str, dict[str, float | int | bool | str]]:
        """Flatten to a JSON-serializable nested dict for ``summary.json``."""
        return {
            "gdna_fl_mean": self.gdna_fl_mean,  # type: ignore[dict-item]
            "INTERGENIC": _density_to_dict(self.intergenic),
            "INTRON": _density_to_dict(self.intron),
            "EXON-INTRON": _density_to_dict(self.exon_intron),
        }


def _density_to_dict(d: GlobalGdnaDensity) -> dict[str, float | int | bool | str]:
    return {
        "rho": d.rho,
        "n_fragments": d.n_fragments,
        "n_fragments_estimated": float(d.n_fragments_estimated),
        "eff_length_bp": d.eff_length_bp,
        "n_regions_used": d.n_regions_used,
        "n_rows_eligible": int(d.n_rows_eligible),
        "strand_active": bool(d.strand_active),
        "rho_uncorrected": float(d.rho_uncorrected),
        "strand_correction_fragments": float(d.strand_correction_fragments),
        "n_strand_informative_regions": int(d.n_strand_informative_regions),
        "strand_informative_exposure_fraction": float(d.strand_informative_exposure_fraction),
        "n_uninf_observed": int(d.n_uninf_observed),
        "kappa": d.kappa.value,
        "kappa_n_regions": d.kappa.n_regions,
        "kappa_fallback": d.kappa.fallback_used,
        "kappa_reason": d.kappa.fallback_reason,
    }


@dataclass(frozen=True, slots=True)
class _ChannelRows:
    """Per-row observations and exposure for one global density channel."""

    type_label: DensityType
    n_observed: np.ndarray
    n_same: np.ndarray
    n_opp: np.ndarray
    n_uninf: np.ndarray
    exposure: np.ndarray
    strand_identifiable: np.ndarray


def compute_global_densities(
    region_df: pd.DataFrame,
    payload: CalibrationScanPayload,
    *,
    gdna_fl: FragmentLengthModel,
    splicing_anchor_tolerance: int = 0,
    strand_summary: StrandSummary | None = None,
) -> GlobalDensityTable:
    """Compute the three global gDNA densities + per-type $\\kappa$.

    Pure function: no I/O, no scanner state, no mutation of inputs.

    The denominator for INTERGENIC and INTRON is the FL-PMF-weighted
    containment effective length (:func:`l_eff_contained`), matching
    the *contained* fragment numerator emitted by the C++ scanner.
    The EXON-INTRON denominator is the per-side boundary-crossing
    exposure :math:`B_\\text{cross}(K) = \\sum_\\ell h(\\ell)\\,
    \\max(\\ell - 2q(K) + 1, 0)`, with :math:`q(K) = \\max(K, 1)`,
    summed over eligible boundary sides; see
    :func:`rigel.calibration._exposure.boundary_crossing_exposure`.
    ``K`` here must equal the boundary tolerance the scanner used
    when populating ``payload.u_left`` / ``payload.u_right``; the
    numerator/denominator must agree.
    """
    n_regions = len(region_df)
    if payload.per_region_counts.shape[0] != n_regions:
        raise ValueError(
            f"compute_global_densities: region_df has {n_regions} rows but "
            f"payload.per_region_counts has {payload.per_region_counts.shape[0]} "
            f"rows.{_ERR_SUFFIX}"
        )
    if payload.u_left.shape != (n_regions,) or payload.u_right.shape != (n_regions,):
        raise ValueError(
            f"compute_global_densities: payload boundary-flux shapes "
            f"{payload.u_left.shape}/{payload.u_right.shape} != ({n_regions},)."
            f"{_ERR_SUFFIX}"
        )
    if not (gdna_fl.mean > 0.0):
        raise ValueError(
            f"compute_global_densities: gdna_fl.mean must be > 0 "
            f"(got {gdna_fl.mean!r}); the scanner-trained gDNA FL model is "
            f"required upstream of the global density pass."
        )
    if strand_summary is None:
        strand_summary = StrandSummary.uninformative()

    types = region_df["type"].to_numpy()
    strands = _region_strands(region_df, n_regions)
    spans = region_df["end"].to_numpy().astype(np.int64) - region_df["start"].to_numpy().astype(
        np.int64
    )
    leff = l_eff_contained(spans, gdna_fl)
    prc = payload.per_region_counts  # (R, 8) int64

    intergenic = _compute_density(
        _channel_intergenic(
            region_mask=(types == int(RegionType.INTERGENIC)),
            counts_col=prc[:, _MASK_INTERGENIC],
            leff=leff,
        ),
        strand_summary=strand_summary,
    )
    intron = _compute_density(
        _channel_intron(
            region_mask=(types == int(RegionType.INTRON)),
            counts_by_orient=payload.intron_counts_by_orient,
            leff=leff,
            strands=strands,
        ),
        strand_summary=strand_summary,
    )
    exon_intron = _compute_density(
        _channel_exon_intron(
            region_mask=(types == int(RegionType.EXON)),
            bf_left=region_df["boundary_flux_left"].to_numpy(),
            bf_right=region_df["boundary_flux_right"].to_numpy(),
            u_left=payload.u_left,
            u_right=payload.u_right,
            u_left_by_orient=payload.u_left_by_orient,
            u_right_by_orient=payload.u_right_by_orient,
            b_cross=boundary_crossing_exposure(
                gdna_fl, splicing_anchor_tolerance=splicing_anchor_tolerance
            ),
            strands=strands,
        ),
        strand_summary=strand_summary,
    )

    return GlobalDensityTable(
        intergenic=intergenic,
        intron=intron,
        exon_intron=exon_intron,
        gdna_fl=gdna_fl,
    )


def _region_strands(region_df: pd.DataFrame, n_regions: int) -> np.ndarray:
    if "strand" not in region_df.columns:
        return np.zeros(n_regions, dtype=np.uint8)
    return region_df["strand"].to_numpy(np.uint8, copy=False)


def _strand_identifiable_rows(strands: np.ndarray) -> np.ndarray:
    return (strands == int(RegionStrand.POS)) | (strands == int(RegionStrand.NEG))


def _channel_intergenic(
    *,
    region_mask: np.ndarray,
    counts_col: np.ndarray,
    leff: np.ndarray,
) -> _ChannelRows:
    sub_counts = counts_col[region_mask].astype(np.int64, copy=False)
    return _ChannelRows(
        type_label="INTERGENIC",
        n_observed=sub_counts,
        n_same=np.zeros_like(sub_counts),
        n_opp=np.zeros_like(sub_counts),
        n_uninf=sub_counts,
        exposure=leff[region_mask].astype(np.float64, copy=False),
        strand_identifiable=np.zeros(sub_counts.shape, dtype=bool),
    )


def _channel_intron(
    *,
    region_mask: np.ndarray,
    counts_by_orient: np.ndarray,
    leff: np.ndarray,
    strands: np.ndarray,
) -> _ChannelRows:
    same = counts_by_orient[:, ORIENT_SAME]
    opp = counts_by_orient[:, ORIENT_OPP]
    uninf = counts_by_orient[:, ORIENT_UNINF]
    observed = same + opp + uninf
    return _ChannelRows(
        type_label="INTRON",
        n_observed=observed[region_mask].astype(np.int64, copy=False),
        n_same=same[region_mask].astype(np.int64, copy=False),
        n_opp=opp[region_mask].astype(np.int64, copy=False),
        n_uninf=uninf[region_mask].astype(np.int64, copy=False),
        exposure=leff[region_mask].astype(np.float64, copy=False),
        strand_identifiable=_strand_identifiable_rows(strands)[region_mask],
    )


def _channel_exon_intron(
    *,
    region_mask: np.ndarray,  # (R,) bool: type == EXON
    bf_left: np.ndarray,  # (R,) bool: boundary_flux_left flag
    bf_right: np.ndarray,  # (R,) bool: boundary_flux_right flag
    u_left: np.ndarray,  # (R,) int64: per-region left-edge flux
    u_right: np.ndarray,  # (R,) int64: per-region right-edge flux
    u_left_by_orient: np.ndarray,
    u_right_by_orient: np.ndarray,
    b_cross: float,  # FL-PMF boundary-crossing exposure per side
    strands: np.ndarray,
) -> _ChannelRows:
    """Capture-aware EXON boundary-flux density.

    Per-region count: ``u_L * 1_L + u_R * 1_R``.
    Per-region effective length: ``(1_L + 1_R) * B_cross``, where
    :math:`B_\\text{cross} = \\sum_\\ell h(\\ell)\\,\\max(\\ell - 1, 0)` is the
    expected number of *strict* boundary-crossing start positions for a
    single fragment under the gDNA fragment-length distribution. This
    replaces the historical ``mean_FL`` factor, which was off by a
    1 bp shift (a fragment of length :math:`\\ell` has :math:`\\ell - 1`
    strict boundary-crossing positions, not :math:`\\ell`).

    The eligibility filter (``1_L | 1_R`` on EXON rows) ensures that
    terminal exon boundaries (TSS/TES) — where capture probes do not
    tile — never enter either side of the ratio.
    """
    left_enabled = region_mask & bf_left.astype(bool)
    right_enabled = region_mask & bf_right.astype(bool)
    sides = left_enabled.astype(np.int64) + right_enabled.astype(np.int64)
    eligible = sides > 0

    per_region_counts = (
        u_left * left_enabled.astype(np.int64) + u_right * right_enabled.astype(np.int64)
    ).astype(np.int64)
    same = (
        u_left_by_orient[:, ORIENT_SAME] * left_enabled.astype(np.int64)
        + u_right_by_orient[:, ORIENT_SAME] * right_enabled.astype(np.int64)
    ).astype(np.int64)
    opp = (
        u_left_by_orient[:, ORIENT_OPP] * left_enabled.astype(np.int64)
        + u_right_by_orient[:, ORIENT_OPP] * right_enabled.astype(np.int64)
    ).astype(np.int64)
    uninf = (
        u_left_by_orient[:, ORIENT_UNINF] * left_enabled.astype(np.int64)
        + u_right_by_orient[:, ORIENT_UNINF] * right_enabled.astype(np.int64)
    ).astype(np.int64)
    per_region_leff = sides.astype(np.float64) * b_cross

    return _ChannelRows(
        type_label="EXON-INTRON",
        n_observed=per_region_counts[eligible].astype(np.int64, copy=False),
        n_same=same[eligible].astype(np.int64, copy=False),
        n_opp=opp[eligible].astype(np.int64, copy=False),
        n_uninf=uninf[eligible].astype(np.int64, copy=False),
        exposure=per_region_leff[eligible].astype(np.float64, copy=False),
        strand_identifiable=_strand_identifiable_rows(strands)[eligible],
    )


def _gdna_count_moment(
    n_same: np.ndarray,
    n_opp: np.ndarray,
    *,
    signed_strand_contrast: float,
) -> np.ndarray:
    """Unbiased per-row gDNA count moments from same/opposite counts."""
    total = n_same.astype(np.float64, copy=False) + n_opp.astype(np.float64, copy=False)
    difference = n_same.astype(np.float64, copy=False) - n_opp.astype(np.float64, copy=False)
    return total - difference / signed_strand_contrast


def _compute_density(
    rows: _ChannelRows,
    *,
    strand_summary: StrandSummary,
) -> GlobalGdnaDensity:
    exposure = rows.exposure.astype(np.float64, copy=False)
    n_observed = int(rows.n_observed.sum())
    total_exposure = float(exposure.sum())
    rho_uncorrected = n_observed / total_exposure if total_exposure > 0.0 else 0.0

    positive_exposure = exposure > 0.0
    informative = rows.strand_identifiable & positive_exposure
    informative_exposure = float(exposure[informative].sum())
    signed_strand_contrast = strand_summary.signed_strand_contrast
    effective_min_contrast = max(
        STRAND_CONTRAST_NUMERICAL_FLOOR,
        strand_summary.signed_strand_contrast_margin(confidence=0.99),
    )
    has_usable_contrast = (
        abs(signed_strand_contrast) >= effective_min_contrast
        and abs(signed_strand_contrast) > 0.0
    )
    strand_active = bool(has_usable_contrast and informative_exposure > 0.0)

    if strand_active:
        corrected_numerator_raw = float(
            _gdna_count_moment(
                rows.n_same[informative],
                rows.n_opp[informative],
                signed_strand_contrast=signed_strand_contrast,
            ).sum()
        )
        corrected_numerator = max(corrected_numerator_raw, 0.0)
        rho = corrected_numerator / informative_exposure
        n_fragments_estimated = rho * total_exposure
        eff_length_bp = informative_exposure
        n_regions_used = int(np.count_nonzero(informative))
    else:
        rho = rho_uncorrected
        n_fragments_estimated = float(n_observed)
        eff_length_bp = total_exposure
        n_regions_used = int(np.count_nonzero(positive_exposure))

    kappa = estimate_kappa(rows.n_observed, exposure, rho_uncorrected)
    exposure_fraction = informative_exposure / total_exposure if total_exposure > 0.0 else 0.0

    return GlobalGdnaDensity(
        type=rows.type_label,
        rho=rho,
        n_fragments=n_observed,
        eff_length_bp=eff_length_bp,
        n_regions_used=n_regions_used,
        kappa=kappa,
        n_fragments_estimated=n_fragments_estimated,
        n_rows_eligible=int(rows.n_observed.size),
        strand_active=strand_active,
        rho_uncorrected=rho_uncorrected,
        strand_correction_fragments=n_fragments_estimated - float(n_observed),
        n_strand_informative_regions=int(np.count_nonzero(informative)),
        strand_informative_exposure_fraction=exposure_fraction,
        n_uninf_observed=int(rows.n_uninf.sum()),
    )


# ---------------------------------------------------------------------------
# Global gDNA rate estimate (projected over entire genome)
# ---------------------------------------------------------------------------


def estimate_global_gdna_fragments(
    global_densities: GlobalDensityTable,
    region_df: pd.DataFrame,
) -> float:
    """Predict total gDNA fragments across the entire genome.

    Projects each density (intergenic, intron, exon-intron) onto the
    FL-PMF-weighted containment effective length of its region class,
    then sums. The exon-intron density (boundary flux) is projected
    onto ALL exon regions, estimating the exonic gDNA that is invisible
    to direct classification (looks like mRNA).

    Parameters
    ----------
    global_densities : GlobalDensityTable
        The three global density estimates from calibration. Carries
        the gDNA FL model used to weight the projection.
    region_df : pd.DataFrame
        The full region partition from the index (columns: start, end, type).

    Returns
    -------
    float
        Predicted total gDNA fragments across the entire genome.
    """
    types = region_df["type"].to_numpy()
    spans = region_df["end"].to_numpy().astype(np.int64) - region_df["start"].to_numpy().astype(
        np.int64
    )
    leff = l_eff_contained(spans, global_densities.gdna_fl)

    ig_mask = types == int(RegionType.INTERGENIC)
    n_gdna_ig = global_densities.intergenic.rho * float(leff[ig_mask].sum())

    in_mask = types == int(RegionType.INTRON)
    n_gdna_in = global_densities.intron.rho * float(leff[in_mask].sum())

    ex_mask = types == int(RegionType.EXON)
    n_gdna_ex = global_densities.exon_intron.rho * float(leff[ex_mask].sum())

    return n_gdna_ig + n_gdna_in + n_gdna_ex
