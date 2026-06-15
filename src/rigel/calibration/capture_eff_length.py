"""Capture-aware EM effective lengths — extend the gDNA IPR contraction to the mRNA/nRNA components.

**The bug it fixes.** The hybrid-capture inverse-participation-ratio (IPR) effective-length contraction
was applied to the **gDNA component only** (``priors.assemble_priors``); mRNA and nRNA components used the
capture-blind full FL-marginal length. Under capture the gDNA component is artificially concentrated and
out-competes nRNA (and mature) for the capture-enriched reads — the nascent→gDNA siphon (nRNA effective
length measured 5–12× the gDNA component's on the siphon loci).

**The fix.** Contract *every* transcript's EM effective length by the **same per-region gDNA enrichment**,
over the transcript's own region set. gDNA is source-uniform within a locus, so its per-region density is
the probe-enrichment pattern — carrying none of the ~10⁴-fold expression dynamic range that would poison a
coverage-based readout. The region set is the regions the transcript's exons overlap: **exon regions for
mRNA** (introns are gaps), the **full span for nRNA** (its single full-span exon covers introns too).

Contracted in **FL-marginal units** (scale the existing ``effective_lengths`` by the enrichment IPR over
its regions), so a uniform enrichment reduces *exactly* to the input length ⇒ capture-off is bit-identical;
only the captured case contracts. No new readout, no new constant — it reuses the calibration's per-region
gDNA mass and the same IPR shape as the gDNA component. See
``docs/calibration/capture_effective_length_design.md``.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ..types import IntervalType
from .priors import _transport_boundary_flux

if TYPE_CHECKING:
    from ..index import TranscriptIndex
    from .region_arrays import RegionArrays
    from .result import CalibrationResult

__all__ = ["transcript_capture_eff_lengths"]


def _exon_region_incidence(
    index: "TranscriptIndex", region_arrays: "RegionArrays"
) -> tuple[np.ndarray, np.ndarray]:
    """``(inc_t, inc_r)`` — per-transcript region membership over its exon set.

    mRNA → the regions its **exon** intervals overlap; nRNA synthetics (absent from the exon-interval
    table) → the regions its **full span** overlaps (from ``t_df``). Region boundaries do NOT always
    coincide with exon edges (only ~90% do — the partition is annotation-WIDE), so an exon/span ``[a, b)``
    is mapped to every region it **overlaps**: the first region whose ``end > a`` (the one containing ``a``)
    through the last region whose ``start < b``. Reading the lower bound from region *starts* alone
    (``searchsorted(starts, a, side='left')``) skips the region that contains ``a`` whenever ``a`` falls in
    a region's interior — dropping fully-contained exons/spans entirely — so ``lo`` is read from the region
    *ends* (matching the accumulator reference's ``searchsorted(side='right')`` containment idiom).
    Annotation-only (sample-independent).
    """
    starts = np.asarray(region_arrays.start, dtype=np.int64)
    ends = np.asarray(region_arrays.end, dtype=np.int64)
    ref_off = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
    name_to_id = index.ref_name_to_id
    inc_t_parts: list[np.ndarray] = []
    inc_r_parts: list[np.ndarray] = []
    seen: set[int] = set()

    def _add(t: int, ref_name: object, a: int, b: int) -> None:
        rid = name_to_id.get(str(ref_name))
        if rid is None:
            return
        lo0, hi0 = int(ref_off[rid]), int(ref_off[rid + 1])
        # overlap of region [start_r, end_r) with [a, b): first region with end_r > a (contains/after a)
        # through the last with start_r < b. (searchsorted(starts, a) would skip the region containing a
        # when a is interior to it, dropping fully-contained exons/spans.)
        lo = lo0 + int(np.searchsorted(ends[lo0:hi0], a, side="right"))
        hi = lo0 + int(np.searchsorted(starts[lo0:hi0], b, side="left"))
        if hi > lo:
            inc_r_parts.append(np.arange(lo, hi, dtype=np.int64))
            inc_t_parts.append(np.full(hi - lo, int(t), dtype=np.int64))

    iv = pd.read_feather(os.path.join(index.index_dir, "intervals.feather"))
    ex = iv[(iv["interval_type"] == int(IntervalType.EXON)) & (iv["t_index"] >= 0)]
    for t, ref_name, a, b in zip(ex["t_index"], ex["ref"], ex["start"], ex["end"], strict=True):
        _add(int(t), ref_name, int(a), int(b))
        seen.add(int(t))

    tdf = index.t_df
    if tdf is not None and "is_synthetic" in tdf.columns:
        syn = tdf[tdf["is_synthetic"].to_numpy(dtype=bool)]
        for t, ref_name, a, b in zip(
            syn["t_index"], syn["ref"], syn["start"], syn["end"], strict=True
        ):
            if int(t) in seen:
                continue
            _add(int(t), ref_name, int(a), int(b))

    if not inc_t_parts:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    return np.concatenate(inc_t_parts), np.concatenate(inc_r_parts)


def transcript_capture_eff_lengths(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    index: "TranscriptIndex",
    fl_eff_lengths: np.ndarray,
) -> np.ndarray:
    """Capture-contract each transcript's EM effective length by the gDNA-enrichment IPR over its regions.

    The contraction is the **same Laplace-smoothed IPR the gDNA component uses**, over the transcript's
    region set, expressed as a factor on the FL-marginal length::

        eff_ipr_t = (G_t + 1)² / [ Σ(m_r²/L_r) + (2·G_t + 1)/span_t ],  capped at span_t   (G_t = Σ m_r)
        eff_em_t  = fl_eff_lengths_t · (eff_ipr_t / span_t)                                 (capped at fl)

    where ``m_r`` = transported per-region gDNA mass, ``L_r`` = ``gdna_geom_len``, ``span_t`` = Σ L_r over
    the transcript's regions. A single O(incidence) pass (``np.add.at``) does every transcript at once.
    Properties (all inherited from the gDNA IPR — no new logic):

    * **uniform gDNA** (capture off) ⇒ ``eff_ipr_t = span_t`` ⇒ factor 1 ⇒ ``eff_em_t = fl`` — bit-identical;
    * **sparse gDNA** ⇒ the Laplace ``+1`` keeps ``eff_ipr_t ≈ span_t`` ⇒ ≈ no contraction (prevents
      spurious over-contraction on thin evidence);
    * **concentrated gDNA** (capture) ⇒ contracts toward the probed regions, matching the gDNA component;
    * **no gDNA over the transcript** ⇒ factor 1 ⇒ ``fl``.

    ``fl_eff_lengths`` is the existing FL-marginal ``effective_lengths`` (the contraction never exceeds it).
    """
    fl = np.asarray(fl_eff_lengths, dtype=np.float64)
    length = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    gdna_region = _transport_boundary_flux(
        calibration.mass_gdna_contained,
        calibration.mass_gdna_left,
        calibration.mass_gdna_right,
        length,
        calibration.gdna_boundary_len,
        np.asarray(region_arrays.ref_id),
    )
    geom = np.maximum(np.asarray(calibration.gdna_geom_len, dtype=np.float64), 1e-9)
    m_sq_over_l = gdna_region**2 / geom

    inc_t, inc_r = _exon_region_incidence(index, region_arrays)
    n_t = fl.shape[0]
    g_sum = np.zeros(n_t)  # G_t  = Σ m_r        (component gDNA mass)
    supp = np.zeros(n_t)  # Σ m_r²/L_r
    span = np.zeros(n_t)  # span_t = Σ L_r
    n_reg = np.zeros(n_t)  # number of regions in the footprint
    if inc_t.size:
        np.add.at(g_sum, inc_t, gdna_region[inc_r])
        np.add.at(supp, inc_t, m_sq_over_l[inc_r])
        np.add.at(span, inc_t, geom[inc_r])
        np.add.at(n_reg, inc_t, 1.0)

    with np.errstate(divide="ignore", invalid="ignore"):
        safe_span = np.maximum(span, 1e-9)
        eff_ipr = np.minimum((g_sum + 1.0) ** 2 / (supp + (2.0 * g_sum + 1.0) / safe_span), span)
        factor_raw = np.where(span > 0.0, eff_ipr / safe_span, 1.0)  # ∈ (0,1]; 1 ⇒ no contraction
        # Evidence-weighted SHRINKAGE of the contraction toward 1 (no contraction). The IPR estimates an
        # effective footprint from the per-region gDNA counts; that estimate is only trustworthy when there
        # are enough fragments to resolve the distribution. Empirical-Bayes weight w = G/(G + n_reg) — one
        # pseudo-fragment per footprint region (the canonical Dirichlet prior over the regions, NO tunable
        # constant): w→1 when the gDNA is abundant relative to the footprint (real capture concentration ⇒
        # full contraction) and w→0 when the gDNA is sparse — a few false-positive fragments in a zero-DNA
        # library, or genuinely thin evidence — so the contraction smoothly SHRINKS back to no-contraction
        # rather than over-contracting on noise. Degrades to ≈ the old factor wherever gDNA is real (G ≫ n_reg).
        w = np.where(n_reg > 0.0, g_sum / (g_sum + n_reg), 0.0)
        factor = w * factor_raw + (1.0 - w)
    return np.minimum(fl * factor, fl)
