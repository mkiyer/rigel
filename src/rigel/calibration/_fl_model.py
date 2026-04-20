"""gDNA fragment-length model — Tukey-weighted, prior-blended (no cliffs)."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from ..frag_length_model import FragmentLengthModel


logger = logging.getLogger(__name__)


def build_gdna_fl_model(
    fl_table: pd.DataFrame,
    stats: dict[str, np.ndarray],
    region_weight: np.ndarray,
    *,
    strand_specificity: float,
    intergenic_fl_model: FragmentLengthModel | None,
    fl_prior_ess: float,
    max_size: int = 1000,
) -> FragmentLengthModel:
    """Build the gDNA fragment-length model.

    Each fragment receives a soft weight ``w = w_density[region_id] · α``
    where ``α`` is 1 for unstranded data, the antisense indicator for
    stranded data, or 1 for unstranded tx_strand=0 regions.  The
    weighted histogram is blended with ``intergenic_fl_model.counts`` as
    a Dirichlet-Multinomial prior with effective sample size
    ``fl_prior_ess`` — a Bayesian replacement for the v3 ``min_ess``
    cliff: when no high-weight gDNA-like fragments are present the
    model equals the intergenic prior, and as evidence accumulates it
    smoothly transitions to the data.

    Always returns a finalised model — never ``None``.  The "no
    information" fallback is the uniform distribution implied by
    Laplace smoothing on an empty histogram.
    """
    n_regions = len(stats["tx_strand"])
    region_weight = np.asarray(region_weight, dtype=np.float64)
    if region_weight.shape[0] != n_regions:
        raise ValueError("region_weight length must match stats arrays")

    model = FragmentLengthModel(max_size=max_size)

    if fl_table is not None and len(fl_table) > 0:
        rid = fl_table["region_id"].to_numpy(dtype=np.intp)
        flen = fl_table["frag_len"].to_numpy(dtype=np.intp)

        # Per-fragment region weight (soft Tukey biweight).
        keep_idx = (rid >= 0) & (rid < n_regions)
        rid = rid[keep_idx]
        flen = flen[keep_idx]
        w = region_weight[rid]

        # Stranded selection: keep antisense fragments in stranded regions
        # at full weight; in unstranded regions (tx_strand == 0) keep all
        # fragments at full weight.  No cliff — this is just a deterministic
        # selector consistent with the strand pathway's likelihood.
        if "frag_strand" in fl_table.columns and float(strand_specificity) > 0.5:
            fstrand = fl_table["frag_strand"].to_numpy(dtype=np.int8)[keep_idx]
            gs = stats["tx_strand"][rid]
            # R1-antisense convention: antisense fragments are POS in +1
            # genes and NEG in -1 genes (matches `_stats.compute_sense_fraction`).
            anti = np.zeros_like(w)
            anti[(gs == 1) & (fstrand == 1)] = 1.0
            anti[(gs == -1) & (fstrand == -1)] = 1.0
            anti[gs == 0] = 1.0
            w = w * anti

        # Bin into the histogram via weighted bincount.
        in_range = (flen >= 0) & (flen <= max_size) & (w > 0)
        if in_range.any():
            sub_flen = flen[in_range]
            sub_w = w[in_range]
            counts = np.bincount(
                sub_flen, weights=sub_w, minlength=max_size + 1
            ).astype(np.float64)
            model.counts += counts[: max_size + 1]
            model._total_weight += float(sub_w.sum())

    have_prior = intergenic_fl_model is not None and intergenic_fl_model.total_weight > 0

    if have_prior:
        model.finalize(
            prior_counts=intergenic_fl_model.counts,
            prior_ess=float(fl_prior_ess),
        )
    else:
        model.finalize()

    logger.info(
        "Built gDNA FL model: total_weight=%.1f, mean=%.1f, prior_ess=%.0f%s",
        model.total_weight,
        model.mean,
        float(fl_prior_ess) if have_prior else 0.0,
        " (intergenic prior blended)" if have_prior else " (no prior)",
    )
    return model
