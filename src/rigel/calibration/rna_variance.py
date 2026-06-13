"""Phase-2a — per-strand RNA ``var~mean`` from splice-junction boundary pairs (CALIBRATION_PLAN §4, Q2).

The RNA mirror of the gDNA ``density_variance_curve``, with the rev-3 correction: **only the two flanking
splice-junction boundaries** are used (they own pure mature-RNA spliced counts), **never the contained
region** (a gDNA+RNA mixture — unusable as an RNA measurement until a first deconvolution pass).

A region that is the **exon side** of an eligible splice junction on **both** boundaries, on the same
strand, with nonzero spliced flux, gives two measurements of the same mature-RNA density (`d_L`, `d_R` =
spliced flux × ``1/fl_mean_rna``). Their disagreement ``¼(d_L−d_R)²`` at mean ``½(d_L+d_R)`` is one
``var~mean`` observation; the robust log-log LOESS (reused from ``density_model.density_variance_curve``)
fits ``σ²_RNA(μ)``. It is **motif-stranded**, so it works at ``ss=0.5`` (the splice strand comes from the
junction motif, not the unspliced read orientation). The fitted curve feeds the propagation process
variance ``Q_RNA`` (the odds-coupling decay) in Phase 2b.

Deviations from the plan: none — this is the plan's Q2, 2-boundary-only. Open: pooled vs per-strand fit
(here per-strand, simplest); the transcript-set dedup (count each boundary-set once) is deferred; where a
strand has ``< _LOESS_MIN_FIT`` pairs the curve is ``NaN`` and the caller falls back to a constant ``Q_RNA``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .density_model import density_variance_curve
from .mature_density import mature_density
from .signature import TS_NEG, TS_POS
from .splice_junction import splice_junction_eligibility

__all__ = ["RnaVariance", "rna_spliced_variance"]


@dataclass(frozen=True, slots=True)
class RnaVariance:
    """Per-region per-strand RNA ``var~mean`` result + the fit-point diagnostics (for plotting)."""

    sigma2_pos: np.ndarray  # float64[R] — σ²_RNA at each region's + mature density (NaN where undefined)
    sigma2_neg: np.ndarray  # float64[R] — σ²_RNA at each region's − mature density
    fit_mu: np.ndarray  # pooled fit-point means ½(d_L+d_R)
    fit_var: np.ndarray  # pooled fit-point raw variances ¼(d_L−d_R)²
    fit_strand: np.ndarray  # +1 / −1 per fit point


def _boundary_spliced_densities(substrate, region_arrays, strand: int, inv_fl: float):
    """Per-region ``(d_left, d_right)`` spliced densities where the region is the EXON side of an eligible
    same-strand splice junction on that boundary (else ``NaN``). Mirrors :func:`mature_density`'s anchoring."""
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    ref_id = np.asarray(region_arrays.ref_id)
    r = sig.shape[0]
    left_spl = (substrate.left.n_spliced_sense + substrate.left.n_spliced_antisense).astype(np.float64)
    right_spl = (substrate.right.n_spliced_sense + substrate.right.n_spliced_antisense).astype(np.float64)
    d_left = np.full(r, np.nan, dtype=np.float64)
    d_right = np.full(r, np.nan, dtype=np.float64)
    for i in range(r):
        if i > 0 and ref_id[i] == ref_id[i - 1]:
            el = splice_junction_eligibility(int(sig[i - 1]), int(sig[i]))
            if el is not None and el.exon_side == "R" and strand in el.strands:
                d_left[i] = left_spl[i] * inv_fl
        if i < r - 1 and ref_id[i] == ref_id[i + 1]:
            el = splice_junction_eligibility(int(sig[i]), int(sig[i + 1]))
            if el is not None and el.exon_side == "L" and strand in el.strands:
                d_right[i] = right_spl[i] * inv_fl
    return d_left, d_right


def rna_spliced_variance(
    substrate, region_arrays, rna_region_eff_len: np.ndarray, fl_mean_rna: float
) -> RnaVariance:
    """Per-strand RNA ``var~mean`` σ²_RNA(μ) from splice-junction boundary pairs (see module docstring).

    ``rna_region_eff_len`` = ``region_eff_length(region_size_bp, rna_pmf)``; ``fl_mean_rna`` =
    ``boundary_eff_length(rna_pmf)`` (the RNA crossing eff-length) — same as :func:`mature_density`.
    """
    md = mature_density(substrate, region_arrays, rna_region_eff_len, fl_mean_rna)
    inv_fl = 1.0 / fl_mean_rna if fl_mean_rna > 0.0 else 0.0
    out: dict[int, np.ndarray] = {}
    fit_mu, fit_var, fit_strand = [], [], []
    for strand, m_s, sv in ((TS_POS, md.m_pos, 1), (TS_NEG, md.m_neg, -1)):
        d_left, d_right = _boundary_spliced_densities(substrate, region_arrays, strand, inv_fl)
        left_ok = np.isfinite(d_left) & (d_left > 0.0)
        right_ok = np.isfinite(d_right) & (d_right > 0.0)
        # 2-boundary only (contained=None): the curve evaluated at each region's mature density m_s.
        out[strand] = density_variance_curve(
            np.asarray(m_s, dtype=np.float64), d_left=d_left, d_right=d_right,
            left_ok=left_ok, right_ok=right_ok,
        )
        pair = left_ok & right_ok
        mu = 0.5 * (d_left[pair] + d_right[pair])
        fit_mu.append(mu)
        fit_var.append(0.25 * (d_left[pair] - d_right[pair]) ** 2)
        fit_strand.append(np.full(mu.shape, sv, dtype=np.int64))
    return RnaVariance(
        sigma2_pos=out[TS_POS], sigma2_neg=out[TS_NEG],
        fit_mu=np.concatenate(fit_mu), fit_var=np.concatenate(fit_var),
        fit_strand=np.concatenate(fit_strand),
    )
