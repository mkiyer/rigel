"""rigel.calibration — Analytical gDNA–RNA deconvolution.

Separates genomic DNA (gDNA) contamination from Total RNA using two
complementary signals, blended by strand specificity (SS):

1. **Strand decomposition** (dominant when SS → 1.0) — For regions
   with unambiguous gene-strand annotation, analytically decompose
   unspliced reads into gDNA and RNA using the known SS.  gDNA is
   strand-symmetric (50/50); RNA follows SS.  Direct formula, no
   iteration.

2. **Density baseline** (dominant when SS → 0.5) — Estimate a global
   gDNA background density λ_G from the low-density tail of the
   unspliced density distribution.  Cap per-region gDNA counts at
   λ_G × L_r.

Blending weight: w = (2·SS − 1)².  At SS=1.0 → pure strand.
At SS=0.5 → pure density.

Calibration does NOT distinguish mature from nascent RNA — that is the
downstream per-locus EM's job via exon/intron geometry.

Outputs
-------
* ``region_e_gdna`` — per-region E[N_gDNA], continuous physical counts.
* ``region_n_total`` — per-region total fragment counts.
* ``gdna_fl_model`` — gDNA fragment-length distribution.
* ``lambda_gdna`` — global gDNA background density (frags/bp).

Strand constants
----------------
Fragment strands use the integer convention from ``constants.h``
and ``types.py``:

    STRAND_POS = 1  (+ strand / forward)
    STRAND_NEG = 2  (− strand / reverse)

Gene-strand context from region metadata:

    gene_strand == +1  → only plus-strand transcripts
    gene_strand == -1  → only minus-strand transcripts
    gene_strand ==  0  → both strands or intergenic (no decomposition)

Under R1-antisense convention (dUTP / TruSeq Stranded):

    gene_strand == +1: antisense unspliced = n_pos (R1 on + = antisense)
    gene_strand == -1: antisense unspliced = n_neg (R1 on − = antisense)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .frag_length_model import FragmentLengthModel
from .types import STRAND_POS, STRAND_NEG

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Small constant to avoid division by zero.
# ---------------------------------------------------------------------------
_EPS: float = 1e-12


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CalibrationResult:
    """Continuous gDNA deconvolution results.

    For every genomic region, provides the expected number of gDNA
    fragments.  No binary classification, no probabilities — physical
    fragment count estimates.

    Downstream usage:

    * ``region_e_gdna`` — per-region E[N_gDNA] used by
      ``locus.compute_locus_gdna_priors()`` for per-locus Dirichlet
      priors.
    * ``region_n_total`` — per-region total fragment counts for RNA
      budget computation.
    * ``gdna_fl_model`` — fragment-length distribution for gDNA
      scoring in ``scoring.FragmentScorer``.
    * ``lambda_gdna`` — global gDNA background density.
    """

    #: Per-region expected gDNA fragment count.
    region_e_gdna: np.ndarray  # (n_regions,) float64 ≥ 0

    #: Per-region total fragment count (unspliced + spliced).
    region_n_total: np.ndarray  # (n_regions,) float64

    #: gDNA fragment-length model for downstream scoring.
    gdna_fl_model: FragmentLengthModel | None

    #: Global gDNA background density (frags / bp).
    lambda_gdna: float

    #: Library strand specificity (echoed for downstream reference).
    strand_specificity: float

    def to_summary_dict(self) -> dict:
        """Return a concise JSON-serializable summary."""
        total_e_gdna = float(self.region_e_gdna.sum())
        total_n = float(self.region_n_total.sum())
        d: dict = {
            "lambda_gdna": round(self.lambda_gdna, 8),
            "total_expected_gdna": round(total_e_gdna, 1),
            "gdna_fraction": round(total_e_gdna / max(total_n, 1.0), 4),
            "strand_specificity": round(self.strand_specificity, 4),
            "gdna_fl_mean": None,
            "gdna_fl_observations": 0,
        }
        if self.gdna_fl_model is not None:
            d["gdna_fl_mean"] = round(self.gdna_fl_model.mean, 2)
            d["gdna_fl_observations"] = self.gdna_fl_model.n_observations
        return d


# ---------------------------------------------------------------------------
# Per-region summary statistics
# ---------------------------------------------------------------------------


def compute_region_stats(
    region_counts: pd.DataFrame,
    region_df: pd.DataFrame,
) -> dict[str, np.ndarray]:
    """Compute per-region summary statistics.

    Parameters
    ----------
    region_counts : pd.DataFrame
        Columns: ``n_unspliced_pos``, ``n_unspliced_neg``,
        ``n_spliced_pos``, ``n_spliced_neg`` (float32).
    region_df : pd.DataFrame
        Region metadata with ``tx_pos``, ``tx_neg``, ``length`` columns.

    Returns
    -------
    dict with keys: n_pos, n_neg, n_unspliced, n_spliced, n_total,
    strand_ratio, splice_rate, density, gene_strand, region_length,
    tx_pos, tx_neg, exon_pos, exon_neg, ref.
    """
    n_pos = region_counts["n_unspliced_pos"].values.astype(np.float64)
    n_neg = region_counts["n_unspliced_neg"].values.astype(np.float64)
    n_spliced = (
        region_counts["n_spliced_pos"].values + region_counts["n_spliced_neg"].values
    ).astype(np.float64)

    n_unspliced = n_pos + n_neg
    n_total = n_unspliced + n_spliced

    strand_ratio = np.full(len(n_pos), np.nan)
    has_unspliced = n_unspliced > 0
    strand_ratio[has_unspliced] = n_pos[has_unspliced] / n_unspliced[has_unspliced]

    splice_rate = np.zeros(len(n_pos), dtype=np.float64)
    has_total = n_total > 0
    splice_rate[has_total] = n_spliced[has_total] / n_total[has_total]

    region_length = region_df["length"].values.astype(np.float64)
    density = np.zeros(len(n_pos), dtype=np.float64)
    has_length = region_length > 0
    density[has_length] = n_total[has_length] / region_length[has_length]

    # Gene strand context: +1 only-sense, -1 only-antisense, 0 ambiguous.
    tx_pos = region_df["tx_pos"].values.astype(bool)
    tx_neg = region_df["tx_neg"].values.astype(bool)
    gene_strand = np.zeros(len(region_df), dtype=np.int8)
    gene_strand[tx_pos & ~tx_neg] = 1
    gene_strand[~tx_pos & tx_neg] = -1

    # Annotation flags.
    exon_pos = (
        region_df["exon_pos"].values.astype(bool)
        if "exon_pos" in region_df.columns
        else np.zeros(len(region_df), dtype=bool)
    )
    exon_neg = (
        region_df["exon_neg"].values.astype(bool)
        if "exon_neg" in region_df.columns
        else np.zeros(len(region_df), dtype=bool)
    )
    ref = (
        region_df["ref"].values
        if "ref" in region_df.columns
        else np.full(len(region_df), "unknown", dtype=object)
    )

    return {
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_unspliced": n_unspliced,
        "n_spliced": n_spliced,
        "n_total": n_total,
        "strand_ratio": strand_ratio,
        "splice_rate": splice_rate,
        "density": density,
        "gene_strand": gene_strand,
        "region_length": region_length,
        "tx_pos": tx_pos,
        "tx_neg": tx_neg,
        "exon_pos": exon_pos,
        "exon_neg": exon_neg,
        "ref": ref,
    }


# ---------------------------------------------------------------------------
# Empirical strand histogram
# ---------------------------------------------------------------------------


def compute_sense_fraction(
    stats: dict[str, np.ndarray],
) -> np.ndarray:
    """Convert per-region strand_ratio to gene-strand-normalised sense fraction.

    R1-antisense convention (dUTP / TruSeq Stranded):

    * ``gene_strand == +1``: RNA reads land on − strand, so
      ``sense_frac = 1 − strand_ratio``.
    * ``gene_strand == -1``: RNA reads land on + strand, so
      ``sense_frac = strand_ratio``.
    * ``gene_strand == 0``: ambiguous — set to NaN (no strand signal).

    For gDNA, sense_frac is ~0.5 regardless of gene_strand.
    For RNA, sense_frac ≈ SS (close to 1.0 for stranded libraries).

    Returns an array of shape ``(n_regions,)`` with NaN where strand
    signal is unavailable (gene_strand == 0 or n_unspliced < 2).
    """
    strand_ratio = stats["strand_ratio"]
    gene_strand = stats["gene_strand"]
    n_unspliced = stats["n_unspliced"]

    sf = np.full(len(strand_ratio), np.nan, dtype=np.float64)

    valid = np.isfinite(strand_ratio) & (n_unspliced >= 2)
    plus = valid & (gene_strand == 1)
    minus = valid & (gene_strand == -1)

    sf[plus] = 1.0 - strand_ratio[plus]
    sf[minus] = strand_ratio[minus]

    return sf


# ---------------------------------------------------------------------------
# Length-weighted percentile
# ---------------------------------------------------------------------------


def _length_weighted_percentile(
    values: np.ndarray,
    weights: np.ndarray,
    percentile: float,
) -> float:
    """Compute a weighted percentile using linear interpolation.

    Parameters
    ----------
    values : array of float
        The values to compute the percentile over.
    weights : array of float (> 0)
        Weights for each value (e.g. region lengths).
    percentile : float
        Target percentile in [0, 100].

    Returns
    -------
    float
        The weighted percentile value.
    """
    mask = weights > 0
    v = values[mask]
    w = weights[mask]
    if len(v) == 0:
        return 0.0

    order = np.argsort(v)
    v = v[order]
    w = w[order]

    cumw = np.cumsum(w)
    total = cumw[-1]
    if total <= 0:
        return 0.0

    # Normalised cumulative weight, centred at midpoint of each bin
    frac = (cumw - 0.5 * w) / total * 100.0
    target = float(np.clip(percentile, frac[0], frac[-1]))
    return float(np.interp(target, frac, v))


# ---------------------------------------------------------------------------
# gDNA fragment-length model builder
# ---------------------------------------------------------------------------


def _build_gdna_fl_model(
    fl_table: pd.DataFrame,
    stats: dict[str, np.ndarray],
    eligible: np.ndarray,
    strand_specificity: float,
    density_percentile: float,
    intergenic_fl_model: FragmentLengthModel | None,
    max_fl: int = 1000,
    min_ess: int = 50,
) -> FragmentLengthModel | None:
    """Build gDNA fragment-length model from high-confidence gDNA fragments.

    Two source pathways, blended by strand specificity weight:

    Stranded (SS > 0.5):
        Antisense unspliced fragments in regions with unambiguous gene
        strand.  At SS=0.95, ~91% of these are gDNA.

    Density-based (SS ≈ 0.5):
        Fragments from low-density regions (bottom Kth percentile of
        unspliced density).

    Falls back to ``intergenic_fl_model`` if insufficient observations.
    """
    if fl_table is None or len(fl_table) == 0:
        return intergenic_fl_model

    gene_strand = stats["gene_strand"]
    has_strand = gene_strand != 0
    region_length = stats["region_length"]
    n_unspliced = stats["n_unspliced"]
    has_frag_strand = "frag_strand" in fl_table.columns

    # Region-level data for density-based selection
    d_unspliced = np.zeros_like(region_length)
    pos = (region_length > 0) & eligible
    d_unspliced[pos] = n_unspliced[pos] / region_length[pos]

    density_threshold = _length_weighted_percentile(
        d_unspliced[eligible],
        region_length[eligible],
        density_percentile,
    )

    low_density_regions = eligible & (d_unspliced <= density_threshold)

    # Collect FL observations from each pathway
    rids = fl_table["region_id"].values
    flens = fl_table["frag_len"].values

    # Stranded pathway: antisense unspliced fragments
    strand_mask = np.zeros(len(fl_table), dtype=bool)
    if has_frag_strand and strand_specificity > 0.5:
        fstrands = fl_table["frag_strand"].values
        for rid_val in range(len(gene_strand)):
            if not eligible[rid_val] or not has_strand[rid_val]:
                continue
            gs = gene_strand[rid_val]
            rid_match = rids == rid_val
            if gs == 1:
                # gene on +: antisense = R1 on + strand = STRAND_POS
                strand_mask |= rid_match & (fstrands == STRAND_POS)
            elif gs == -1:
                # gene on −: antisense = R1 on − strand = STRAND_NEG
                strand_mask |= rid_match & (fstrands == STRAND_NEG)

    # Density pathway: fragments from low-density regions
    density_mask = np.isin(rids, np.where(low_density_regions)[0])

    # Blend: use stranded if available, supplement with density
    w = (2.0 * strand_specificity - 1.0) ** 2
    n_strand = int(strand_mask.sum())
    n_density = int(density_mask.sum())

    if n_strand >= min_ess and w > 0:
        # Sufficient stranded observations
        selected = strand_mask
    elif n_strand > 0 and n_density > 0:
        # Combine both sources
        selected = strand_mask | density_mask
    elif n_density >= min_ess:
        selected = density_mask
    else:
        logger.info(
            f"[CAL-FL] Insufficient FL observations "
            f"(strand={n_strand}, density={n_density}); "
            f"using intergenic fallback"
        )
        return intergenic_fl_model

    selected_fl = flens[selected]
    valid = (selected_fl > 0) & (selected_fl <= max_fl)
    if valid.sum() < min_ess:
        return intergenic_fl_model

    model = FragmentLengthModel(max_size=max_fl)
    hist = np.bincount(selected_fl[valid].astype(np.intp), minlength=max_fl + 1)
    model.counts[: max_fl + 1] = hist[: max_fl + 1].astype(np.float64)
    model._total_weight = float(hist.sum())
    model.finalize()
    return model


# ---------------------------------------------------------------------------
# Main calibration entry point
# ---------------------------------------------------------------------------


def calibrate_gdna(
    region_counts: pd.DataFrame,
    fl_table: pd.DataFrame,
    region_df: pd.DataFrame,
    strand_specificity: float,
    density_percentile: float = 10.0,
    intergenic_fl_model: FragmentLengthModel | None = None,
    **kwargs,
) -> CalibrationResult:
    """Analytical gDNA–RNA deconvolution.

    Separates gDNA from Total RNA using strand decomposition (primary)
    and density baseline (secondary), blended by strand specificity.

    Parameters
    ----------
    region_counts : pd.DataFrame
        Per-region fragment counts (4 columns: unspliced_pos/neg,
        spliced_pos/neg).
    fl_table : pd.DataFrame
        Fragment-length observations with columns ``region_id``,
        ``frag_len``, and optionally ``frag_strand``.
    region_df : pd.DataFrame
        Region metadata (``tx_pos``, ``tx_neg``, ``length``, etc.).
    strand_specificity : float
        Library strand specificity ∈ [0.5, 1.0].
    density_percentile : float
        Percentile of unspliced density for λ_G estimation (unstranded).
    intergenic_fl_model : FragmentLengthModel or None
        Fallback FL model when insufficient gDNA observations.
    **kwargs
        Ignored (accepts old parameters for API compatibility during
        transition).

    Returns
    -------
    CalibrationResult
    """
    if kwargs:
        ignored = ", ".join(sorted(kwargs.keys()))
        logger.debug(f"[CAL] Ignoring deprecated parameters: {ignored}")

    SS = float(np.clip(strand_specificity, 0.5, 1.0))
    stats = compute_region_stats(region_counts, region_df)
    n_regions = len(stats["n_total"])

    # -- Eligibility --
    eligible = (stats["n_total"] > 0) & (stats["region_length"] > 0)
    n_eligible = int(eligible.sum())

    # -- Early exit: no data --
    if n_eligible == 0:
        logger.warning("[CAL] No eligible regions — returning zero gDNA")
        return CalibrationResult(
            region_e_gdna=np.zeros(n_regions, dtype=np.float64),
            region_n_total=stats["n_total"].astype(np.float64),
            gdna_fl_model=intergenic_fl_model,
            lambda_gdna=0.0,
            strand_specificity=SS,
        )

    n_unspliced = stats["n_unspliced"]
    n_pos = stats["n_pos"]
    n_neg = stats["n_neg"]
    gene_strand = stats["gene_strand"]
    region_length = stats["region_length"]
    has_strand = gene_strand != 0

    # ==================================================================
    # Strand pathway: E[N_gDNA]^strand per region
    # ==================================================================
    #
    # For regions with unambiguous gene strand, decompose unspliced
    # reads into gDNA and RNA using the strand balance.
    #
    # Under R1-antisense convention:
    #   gene_strand == +1: antisense unspliced = n_pos
    #   gene_strand == -1: antisense unspliced = n_neg
    #
    # Formula: E[N_gDNA] = max(0, (N_anti - N_unspliced·(1-SS)) / (SS-0.5))

    n_anti = np.where(
        gene_strand == 1,
        n_pos,
        np.where(gene_strand == -1, n_neg, 0.0),
    )

    denom = max(SS - 0.5, _EPS)
    e_gdna_strand = np.maximum(
        0.0,
        (n_anti - n_unspliced * (1.0 - SS)) / denom,
    )
    # Zero out regions without strand annotation or ineligible
    e_gdna_strand[~has_strand | ~eligible] = 0.0
    # Cap at observed unspliced count
    e_gdna_strand = np.minimum(e_gdna_strand, n_unspliced)

    # -- Strand-derived λ_G --
    strand_mask = has_strand & eligible
    strand_e_total = float(e_gdna_strand[strand_mask].sum())
    strand_l_total = float(region_length[strand_mask].sum())
    lambda_g_strand = strand_e_total / max(strand_l_total, 1.0)

    # ==================================================================
    # Density pathway: E[N_gDNA]^density per region
    # ==================================================================
    #
    # gDNA produces a roughly uniform background density λ_G.
    # Estimate λ_G from the low-density tail of unspliced density.
    #
    # Requires a minimum number of eligble regions to form a meaningful
    # density distribution.  With too few regions, the percentile is
    # not informative (e.g. single expressed gene → its own density
    # becomes the "background").
    _MIN_DENSITY_REGIONS = 20

    d_unspliced = np.zeros(n_regions, dtype=np.float64)
    pos = (region_length > 0) & eligible
    d_unspliced[pos] = n_unspliced[pos] / region_length[pos]

    if n_eligible >= _MIN_DENSITY_REGIONS:
        lambda_g_density = _length_weighted_percentile(
            d_unspliced[eligible],
            region_length[eligible],
            density_percentile,
        )
    else:
        lambda_g_density = 0.0
        logger.info(
            f"[CAL] Density pathway disabled: only {n_eligible} eligible "
            f"regions (need {_MIN_DENSITY_REGIONS})"
        )

    e_gdna_density = np.minimum(
        n_unspliced,
        lambda_g_density * region_length,
    )
    e_gdna_density[~eligible] = 0.0

    # ==================================================================
    # Blending
    # ==================================================================
    #
    # w = (2·SS − 1)²: smoothly transitions from pure density (SS=0.5)
    # to pure strand (SS=1.0).

    w = (2.0 * SS - 1.0) ** 2

    # Blended λ_G
    lambda_gdna = w * lambda_g_strand + (1.0 - w) * lambda_g_density

    # Per-region blend (regions without strand info use density only)
    w_per_region = np.where(has_strand, w, 0.0)
    e_gdna = w_per_region * e_gdna_strand + (1.0 - w_per_region) * e_gdna_density

    # ==================================================================
    # gDNA fragment-length model
    # ==================================================================

    gdna_fl = _build_gdna_fl_model(
        fl_table,
        stats,
        eligible,
        SS,
        density_percentile,
        intergenic_fl_model,
    )

    # -- Logging --
    total_e_gdna = float(e_gdna.sum())
    total_unspliced = float(n_unspliced[eligible].sum())
    logger.info(
        f"[CAL] Strand pathway: λ_G={lambda_g_strand:.2e} (from {int(strand_mask.sum())} regions)"
    )
    logger.info(
        f"[CAL] Density pathway: λ_G={lambda_g_density:.2e} (percentile={density_percentile})"
    )
    logger.info(
        f"[CAL] Blended: λ_G={lambda_gdna:.2e}, w_strand={w:.3f}, "
        f"E[gDNA]={total_e_gdna:.0f} / {total_unspliced:.0f} unspliced"
    )

    return CalibrationResult(
        region_e_gdna=e_gdna,
        region_n_total=stats["n_total"].astype(np.float64),
        gdna_fl_model=gdna_fl,
        lambda_gdna=lambda_gdna,
        strand_specificity=SS,
    )
