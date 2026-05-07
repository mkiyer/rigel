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
from ._kappa import KappaEstimate, estimate_kappa
from .regions import RegionType
from .scan_payload import (
    MASK_INTERGENIC as _MASK_INTERGENIC,
    MASK_INTRON as _MASK_INTRON,
    CalibrationScanPayload,
)


# (mask constants are imported from scan_payload — single source of truth)

DensityType = Literal["INTERGENIC", "INTRON", "EXON-INTRON"]

_ERR_SUFFIX = " Rebuild the index or rerun the scan."


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
    rho: float                 # fragments / bp
    n_fragments: int           # numerator
    eff_length_bp: float       # denominator
    n_regions_used: int        # rows that contributed to the denominator
    kappa: KappaEstimate       # NB overdispersion + fallback diagnostic


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
    gdna_fl: FragmentLengthModel    # the model used to weight L_eff

    @property
    def gdna_fl_mean(self) -> float:
        """Mean of the gDNA FL distribution; used for capture-window geometry
        in the EXON-INTRON branch and as a summary statistic."""
        return float(self.gdna_fl.mean)

    def to_summary_dict(self) -> dict[str, dict[str, float | int | bool | str]]:
        """Flatten to a JSON-serializable nested dict for ``summary.json``."""
        return {
            "gdna_fl_mean": self.gdna_fl_mean,  # type: ignore[dict-item]
            "INTERGENIC":   _density_to_dict(self.intergenic),
            "INTRON":       _density_to_dict(self.intron),
            "EXON-INTRON":  _density_to_dict(self.exon_intron),
        }


def _density_to_dict(d: GlobalGdnaDensity) -> dict[str, float | int | bool | str]:
    return {
        "rho":             d.rho,
        "n_fragments":     d.n_fragments,
        "eff_length_bp":   d.eff_length_bp,
        "n_regions_used":  d.n_regions_used,
        "kappa":           d.kappa.value,
        "kappa_n_regions": d.kappa.n_regions,
        "kappa_fallback":  d.kappa.fallback_used,
        "kappa_reason":    d.kappa.fallback_reason,
    }


def compute_global_densities(
    region_df: pd.DataFrame,
    payload: CalibrationScanPayload,
    *,
    gdna_fl: FragmentLengthModel,
) -> GlobalDensityTable:
    """Compute the three global gDNA densities + per-type $\\kappa$.

    Pure function: no I/O, no scanner state, no mutation of inputs.

    The denominator for INTERGENIC and INTRON is the FL-PMF-weighted
    containment effective length (:func:`l_eff_contained`), matching
    the *contained* fragment numerator emitted by the C++ scanner.
    The EXON-INTRON denominator is a capture-window geometry that
    uses the FL mean (:attr:`FragmentLengthModel.mean`).
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

    types = region_df["type"].to_numpy()
    spans = (
        region_df["end"].to_numpy().astype(np.int64)
        - region_df["start"].to_numpy().astype(np.int64)
    )
    leff = l_eff_contained(spans, gdna_fl)
    prc = payload.per_region_counts          # (R, 8) int64

    intergenic = _density_simple(
        type_label="INTERGENIC",
        region_mask=(types == int(RegionType.INTERGENIC)),
        counts_col=prc[:, _MASK_INTERGENIC],
        leff=leff,
    )
    intron = _density_simple(
        type_label="INTRON",
        region_mask=(types == int(RegionType.INTRON)),
        counts_col=prc[:, _MASK_INTRON],
        leff=leff,
    )
    exon_intron = _density_exon_intron(
        region_mask=(types == int(RegionType.EXON)),
        bf_left=region_df["boundary_flux_left"].to_numpy(),
        bf_right=region_df["boundary_flux_right"].to_numpy(),
        u_left=payload.u_left,
        u_right=payload.u_right,
        gdna_fl_mean=float(gdna_fl.mean),
    )

    return GlobalDensityTable(
        intergenic=intergenic,
        intron=intron,
        exon_intron=exon_intron,
        gdna_fl=gdna_fl,
    )


def _density_simple(
    *,
    type_label: DensityType,
    region_mask: np.ndarray,    # (R,) bool: which rows are this type
    counts_col: np.ndarray,     # (R,) int64: per-region count for this mask
    leff: np.ndarray,           # (R,) float64: per-region L_eff
) -> GlobalGdnaDensity:
    """INTERGENIC / INTRON: one count column, one denominator column."""
    sub_counts = counts_col[region_mask]
    sub_leff = leff[region_mask]
    n_frag = int(sub_counts.sum())
    eff_len = float(sub_leff.sum())
    rho = n_frag / eff_len if eff_len > 0.0 else 0.0
    n_used = int(np.count_nonzero(sub_leff > 0.0))
    kappa = estimate_kappa(sub_counts, sub_leff, rho)
    return GlobalGdnaDensity(
        type=type_label,
        rho=rho,
        n_fragments=n_frag,
        eff_length_bp=eff_len,
        n_regions_used=n_used,
        kappa=kappa,
    )


def _density_exon_intron(
    *,
    region_mask: np.ndarray,    # (R,) bool: type == EXON
    bf_left: np.ndarray,        # (R,) bool: boundary_flux_left flag
    bf_right: np.ndarray,       # (R,) bool: boundary_flux_right flag
    u_left: np.ndarray,         # (R,) int64: per-region left-edge flux
    u_right: np.ndarray,        # (R,) int64: per-region right-edge flux
    gdna_fl_mean: float,
) -> GlobalGdnaDensity:
    """Capture-aware EXON boundary-flux density.

    Per-region count: ``u_L * 1_L + u_R * 1_R``.
    Per-region effective length: ``(1_L + 1_R) * gdna_fl_mean``.

    The eligibility filter (``1_L | 1_R`` on EXON rows) ensures that
    terminal exon boundaries (TSS/TES) — where capture probes do not
    tile — never enter either side of the ratio.
    """
    sides = bf_left.astype(np.int64) + bf_right.astype(np.int64)  # 0, 1, or 2
    eligible = region_mask & (sides > 0)

    per_region_counts = (
        u_left * bf_left.astype(np.int64) +
        u_right * bf_right.astype(np.int64)
    ).astype(np.int64)
    per_region_leff = (sides.astype(np.float64) * gdna_fl_mean)

    sub_counts = per_region_counts[eligible]
    sub_leff = per_region_leff[eligible]
    n_frag = int(sub_counts.sum())
    eff_len = float(sub_leff.sum())
    rho = n_frag / eff_len if eff_len > 0.0 else 0.0
    n_used = int(np.count_nonzero(eligible))
    kappa = estimate_kappa(sub_counts, sub_leff, rho)
    return GlobalGdnaDensity(
        type="EXON-INTRON",
        rho=rho,
        n_fragments=n_frag,
        eff_length_bp=eff_len,
        n_regions_used=n_used,
        kappa=kappa,
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
    spans = (
        region_df["end"].to_numpy().astype(np.int64)
        - region_df["start"].to_numpy().astype(np.int64)
    )
    leff = l_eff_contained(spans, global_densities.gdna_fl)

    ig_mask = types == int(RegionType.INTERGENIC)
    n_gdna_ig = global_densities.intergenic.rho * float(leff[ig_mask].sum())

    in_mask = types == int(RegionType.INTRON)
    n_gdna_in = global_densities.intron.rho * float(leff[in_mask].sum())

    ex_mask = types == int(RegionType.EXON)
    n_gdna_ex = global_densities.exon_intron.rho * float(leff[ex_mask].sum())

    return n_gdna_ig + n_gdna_in + n_gdna_ex
