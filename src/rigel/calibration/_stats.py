"""Per-region summary statistics (shared helpers, no behavioural change)."""

from __future__ import annotations

import numpy as np
import pandas as pd


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
        Region metadata.  Must have ``length`` and ``tx_pos`` /
        ``tx_neg`` columns.  ``mappable_effective_length`` (float32)
        is read if present and otherwise defaults to ``length`` (the
        ``--no-mappability`` index path supplies this column directly).

    Returns
    -------
    dict with keys: n_pos, n_neg, n_unspliced, n_spliced, n_total,
    strand_ratio, splice_rate, density, gene_strand, region_length,
    mappable_bp, tx_pos, tx_neg, exon_pos, exon_neg, ref.

    Notes
    -----
    ``mappable_bp`` is the Poisson-rate denominator for the region:
    the number of uniquely-mappable bases available to generate
    fragments.  Falls back to raw ``length`` when the index does not
    provide ``mappable_effective_length``.
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

    # Mappable bp: Poisson-rate denominator per region.  Uses the
    # index-provided mappable effective length when available; older
    # indexes without that column fall back to raw span.
    if "mappable_effective_length" in region_df.columns:
        mappable_bp = region_df["mappable_effective_length"].values.astype(np.float64)
    else:
        mappable_bp = region_length.copy()

    density = np.zeros(len(n_pos), dtype=np.float64)
    has_length = region_length > 0
    density[has_length] = n_total[has_length] / region_length[has_length]

    tx_pos = region_df["tx_pos"].values.astype(bool)
    tx_neg = region_df["tx_neg"].values.astype(bool)
    gene_strand = np.zeros(len(region_df), dtype=np.int8)
    gene_strand[tx_pos & ~tx_neg] = 1
    gene_strand[~tx_pos & tx_neg] = -1

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
        "mappable_bp": mappable_bp,
        "tx_pos": tx_pos,
        "tx_neg": tx_neg,
        "exon_pos": exon_pos,
        "exon_neg": exon_neg,
        "ref": ref,
    }


def compute_sense_fraction(
    stats: dict[str, np.ndarray],
) -> np.ndarray:
    """Convert per-region strand_ratio to gene-strand-normalised sense fraction.

    R1-antisense convention (dUTP / TruSeq Stranded):

    * ``gene_strand == +1``: RNA reads land on − strand, so
      ``sense_frac = 1 − strand_ratio``.
    * ``gene_strand == -1``: RNA reads land on + strand, so
      ``sense_frac = strand_ratio``.
    * ``gene_strand == 0``: ambiguous — set to NaN.

    For gDNA, sense_frac is ~0.5 regardless of gene_strand.
    For RNA, sense_frac ≈ SS.
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
