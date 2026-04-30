"""SRD v1 calibration orchestrator (``calibrate_gdna``).

Single buffer walk → categorize → pool assembly → 1-D mixture fit →
EB-smoothed FL models. Library-agnostic.
"""

from __future__ import annotations

import logging

import numpy as np

from ..buffer import FragmentBuffer
from ..frag_length_model import FragmentLengthModels
from ..index import TranscriptIndex
from ..splice import SPLICE_UNSPLICED
from ._categorize import (
    FragmentCategory,
    FragmentStrand,
    N_CATEGORIES,
    N_FRAGMENT_STRANDS,
    categorize_chunk,
)
from ._fl_empirical_bayes import build_gdna_fl, build_global_fl, build_rna_fl
from ._fl_mixture import fit_fl_mixture
from ._result import CalibrationResult, GdnaFlQuality


logger = logging.getLogger(__name__)


# Quality thresholds.
_PI_MIN_GOOD = 0.02


def calibrate_gdna(
    buffer: FragmentBuffer,
    index: TranscriptIndex,
    frag_length_models: FragmentLengthModels,
    strand_specificity: float,
    *,
    read1_sense: bool = True,
    exon_fit_tolerance_bp: int = 0,
    fl_prior_ess: float = 500.0,
    max_iter: int = 1000,
    tol: float = 1e-4,
) -> CalibrationResult:
    """Run SRD v1 calibration over a finalized fragment buffer.

    Parameters
    ----------
    buffer
        Finalized :class:`FragmentBuffer`.
    index
        Loaded :class:`TranscriptIndex`.
    frag_length_models
        Trained :class:`FragmentLengthModels`. Provides ``global_model``
        (EB anchor), ``rna_model`` (spliced histogram), and ``max_size``.
    strand_specificity
        Library SS estimate from the strand model. Passed through to the
        result for diagnostics; the orchestrator does not branch on it.
        The pool is library-agnostic by construction.

    Other Parameters
    ----------------
    exon_fit_tolerance_bp
        Tolerance for the exon-fit geometric test (Pass 0).
    fl_prior_ess
        Dirichlet pseudo-observation count for shrinking ``RNA_FL`` and
        ``gDNA_FL`` toward ``global_FL``. Small pools naturally collapse
        to the global FL prior; no separate ``min_pool_size`` gate is
        needed.
    """
    # Per-transcript strand lookup, needed for the per-T EXON_CONTAINED
    # check inside `categorize_chunk`.  After this snapshot the full
    # index object is unused (the schema-change comment from Phase 3
    # still applies to everything else).
    t_strand_arr = np.asarray(index.t_to_strand_arr)
    del index
    max_size = int(frag_length_models.max_size)
    n_bins = max_size + 1

    # ---- Pass 0: walk chunks, categorize geometrically ------------------
    # Categorization itself filters to (num_hits == 1) AND (splice_type ==
    # UNSPLICED).  We separately tally multimappers and spliced fragments
    # for diagnostics.
    cat_chunks: list[np.ndarray] = []
    strand_chunks: list[np.ndarray] = []
    len_chunks: list[np.ndarray] = []
    n_multimap_excluded = 0
    n_spliced = 0

    for chunk in buffer.iter_chunks():
        unique_mask = chunk.num_hits == 1
        n_multimap_excluded += int((~unique_mask).sum())
        # Count spliced unique-mapper fragments for the QC warning.
        n_spliced += int(
            (unique_mask & (chunk.splice_type != np.uint8(SPLICE_UNSPLICED))).sum()
        )
        if not unique_mask.any():
            continue

        cc = categorize_chunk(
            chunk,
            t_strand_arr=t_strand_arr,
            exon_fit_tolerance_bp=exon_fit_tolerance_bp,
            read1_sense=read1_sense,
        )
        if not cc.keep.any():
            continue
        cat_chunks.append(cc.category[cc.keep])
        strand_chunks.append(cc.strand[cc.keep])
        len_chunks.append(cc.frag_length[cc.keep])

    global_counts = np.asarray(frag_length_models.global_model.counts, dtype=np.float64)
    spliced_counts = np.asarray(frag_length_models.rna_model.counts, dtype=np.float64)

    # ---- Empty buffer / all-multimapper short-circuit ------------------
    if not cat_chunks:
        return _fallback_result(
            global_counts=global_counts,
            spliced_counts=spliced_counts,
            max_size=max_size,
            strand_specificity=strand_specificity,
            n_multimap_excluded=n_multimap_excluded,
            exon_fit_tolerance_bp=exon_fit_tolerance_bp,
            fl_prior_ess=fl_prior_ess,
            extra={"warning": "empty buffer (no unique fragments)"},
        )

    category = np.concatenate(cat_chunks)
    strand = np.concatenate(strand_chunks)
    frag_length = np.concatenate(len_chunks)

    # 2-D (category, strand) tally; categorize_chunk only labels
    # unique-mapper UNSPLICED fragments with known read strand, so this
    # matches the pre-filter.
    flat = (
        category.astype(np.intp) * N_FRAGMENT_STRANDS + strand.astype(np.intp)
    )
    cat_counts = np.bincount(
        flat, minlength=N_CATEGORIES * N_FRAGMENT_STRANDS
    ).astype(np.int64).reshape((N_CATEGORIES, N_FRAGMENT_STRANDS))

    # ---- Pass 1: pool assembly -----------------------------------------
    # Pool = INTERGENIC ∪ INTRONIC ∪ EXON_INCOMPATIBLE.
    # Library-agnostic; no SS dependence.  EXON_CONTAINED is treated as
    # mature mRNA and held out -- at imperfect SS it is dominated by RNA
    # from highly-expressed loci and would dilute the gDNA_FL signal we
    # are trying to recover.  SRD v2 intentionally does not subdivide on
    # the antisense sub-label here; that is the v2.1 follow-up
    # (SS-gated antisense-exonic pool inclusion).
    # SRD v3 Phase 1: strand sub-label is now in transcript frame
    # (SENSE / ANTISENSE / AMBIG); the pool membership rule itself is
    # still strand-agnostic.  Per-bucket strand asymmetry is surfaced
    # via the diagnostic counters below; the upcoming v3 Phase 2 joint
    # mixture will consume the strand axis directly.
    pool_mask = (
        (category == int(FragmentCategory.INTERGENIC))
        | (category == int(FragmentCategory.INTRONIC))
        | (category == int(FragmentCategory.EXON_INCOMPATIBLE))
    )
    pool_lens = frag_length[pool_mask]
    n_pool_categorized = int(pool_mask.sum())

    # SRD v2 Phase 1: drop out-of-range fragments instead of silently
    # clipping into bin 0 / bin max_size. The previous `np.clip(pool_lens,
    # 0, max_size)` mapped per-candidate `-1` sentinels into bin 0 and
    # saturated long fragments into bin `max_size`, producing a
    # non-physical gDNA_FL with mode=0. With `frag_length` now sourced
    # from `genomic_footprint` (always >= 0), the only out-of-range case
    # is fragments longer than `max_size`.
    #
    # POLICY (SRD v2 Phase 7): fragments with `genomic_footprint >
    # max_frag_length` are biologically anomalous (true chimeras,
    # unannotated splicing, read-through into intergenic, annotation
    # errors). We do not attempt to root-cause classify them — the hard
    # decision is to exclude them from FL training entirely. The dropped
    # count is reported via `n_pool_dropped_out_of_range` for visibility.
    in_range = (pool_lens >= 0) & (pool_lens <= max_size)
    n_pool_dropped_out_of_range = int(pool_lens.size - in_range.sum())
    pool_hist = np.bincount(
        pool_lens[in_range].astype(np.intp),
        minlength=n_bins,
    ).astype(np.float64)

    # SRD v2: zero-candidate (intergenic) fragments flow through the
    # buffer as INTERGENIC and are already counted in `pool_lens`.
    #
    # `n_pool` is the denominator that actually fed the 1-D mixture EM
    # (i.e., in-range pool only). OOR fragments are tracked separately
    # via `n_pool_dropped_out_of_range`.
    n_pool = n_pool_categorized - n_pool_dropped_out_of_range

    # ---- Diagnostic: intronic strand asymmetry (nRNA pollution surface).
    # SRD v3 Phase 1: now expressed in transcript frame, so SENSE >>
    # ANTISENSE in a stranded library is a direct nRNA-contamination
    # signal (where v2 saw a meaningless POS vs NEG asymmetry tied to
    # which genomic strand the local genes happened to live on).
    n_pool_intronic_strand_sense = int(
        cat_counts[int(FragmentCategory.INTRONIC), int(FragmentStrand.SENSE)]
    )
    n_pool_intronic_strand_antisense = int(
        cat_counts[int(FragmentCategory.INTRONIC), int(FragmentStrand.ANTISENSE)]
    )
    n_pool_intronic_strand_ambig = int(
        cat_counts[int(FragmentCategory.INTRONIC), int(FragmentStrand.AMBIG)]
    )

    # ---- Pass 2: 1-D mixture EM (always run; EB handles small pools) ---
    rna_probs = _spliced_probs(spliced_counts, global_counts, n_bins)

    quality: GdnaFlQuality
    extra: dict[str, object] = {}

    fit = fit_fl_mixture(pool_hist, rna_probs, max_iter=max_iter, tol=tol)
    gdna_counts: np.ndarray | None = fit.gdna_counts
    pi = float(fit.pi)
    mixture_converged = bool(fit.converged)
    mixture_iterations = int(fit.n_iter)

    if not mixture_converged:
        # Hard failure: degenerate input. Fall back to global FL shape.
        quality = "fallback"
        gdna_counts = None
        extra["fallback_reason"] = "mixture did not converge"
    elif pi < _PI_MIN_GOOD:
        # Converged but no meaningful gDNA signal.
        quality = "weak"
    else:
        quality = "good"

    # QC-fail: near-zero spliced fragments → loud warning.
    if n_spliced < 100:
        extra["warning_low_spliced"] = (
            f"only {n_spliced} spliced fragments; RNA_FL collapses to global FL"
        )
        logger.warning(
            "[CAL-SRD] only %d spliced fragments — RNA_FL collapses to global FL",
            n_spliced,
        )

    # ---- Pass 3: build the three FL models -----------------------------
    global_fl = build_global_fl(global_counts, max_size)
    rna_fl = build_rna_fl(
        spliced_counts, global_counts, max_size, prior_ess=fl_prior_ess
    )
    gdna_fl = build_gdna_fl(
        gdna_counts, global_counts, max_size, prior_ess=fl_prior_ess
    )

    logger.info(
        "[CAL-SRD] quality=%s SS=%.3f n_pool=%d "
        "pi=%.4f conv=%s n_iter=%d categories=%s n_multimap_excluded=%d "
        "n_spliced=%d n_pool_dropped_out_of_range=%d "
        "intronic_strand=(sense=%d,antisense=%d,ambig=%d)",
        quality,
        float(strand_specificity),
        n_pool,
        pi,
        mixture_converged,
        mixture_iterations,
        cat_counts.tolist(),
        n_multimap_excluded,
        n_spliced,
        n_pool_dropped_out_of_range,
        n_pool_intronic_strand_sense,
        n_pool_intronic_strand_antisense,
        n_pool_intronic_strand_ambig,
    )

    return CalibrationResult(
        gdna_fl_model=gdna_fl,
        rna_fl_model=rna_fl,
        global_fl_model=global_fl,
        gdna_fl_quality=quality,
        strand_specificity=float(strand_specificity),
        category_counts=cat_counts,
        n_multimap_excluded=n_multimap_excluded,
        n_spliced=n_spliced,
        n_pool=n_pool,
        pi_pool=pi,
        mixture_converged=mixture_converged,
        mixture_iterations=mixture_iterations,
        n_pool_intronic_strand_sense=n_pool_intronic_strand_sense,
        n_pool_intronic_strand_antisense=n_pool_intronic_strand_antisense,
        n_pool_intronic_strand_ambig=n_pool_intronic_strand_ambig,
        n_pool_dropped_out_of_range=n_pool_dropped_out_of_range,
        exon_fit_tolerance_bp=int(exon_fit_tolerance_bp),
        fl_prior_ess=float(fl_prior_ess),
        extra=extra,
    )


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _spliced_probs(
    spliced_counts: np.ndarray, global_counts: np.ndarray, n_bins: int
) -> np.ndarray:
    """Probability vector from spliced histogram, falling back to global FL."""
    counts = spliced_counts if float(spliced_counts.sum()) > 0 else global_counts
    counts = np.asarray(counts, dtype=np.float64)
    if counts.size != n_bins:
        adj = np.zeros(n_bins, dtype=np.float64)
        m = min(counts.size, n_bins)
        adj[:m] = counts[:m]
        counts = adj
    s = float(counts.sum())
    if s <= 0:
        return np.full(n_bins, 1.0 / n_bins, dtype=np.float64)
    return counts / s


def _fallback_result(
    *,
    global_counts: np.ndarray,
    spliced_counts: np.ndarray,
    max_size: int,
    strand_specificity: float,
    n_multimap_excluded: int,
    exon_fit_tolerance_bp: int,
    fl_prior_ess: float,
    extra: dict[str, object],
) -> CalibrationResult:
    """Build a result where ``gDNA_FL = global_FL``."""
    global_fl = build_global_fl(global_counts, max_size)
    rna_fl = build_rna_fl(
        spliced_counts, global_counts, max_size, prior_ess=fl_prior_ess
    )
    gdna_fl = build_gdna_fl(
        None, global_counts, max_size, prior_ess=fl_prior_ess
    )
    return CalibrationResult(
        gdna_fl_model=gdna_fl,
        rna_fl_model=rna_fl,
        global_fl_model=global_fl,
        gdna_fl_quality="fallback",
        strand_specificity=float(strand_specificity),
        category_counts=np.zeros((N_CATEGORIES, N_FRAGMENT_STRANDS), dtype=np.int64),
        n_multimap_excluded=int(n_multimap_excluded),
        n_spliced=0,
        n_pool=0,
        pi_pool=0.0,
        mixture_converged=False,
        mixture_iterations=0,
        n_pool_intronic_strand_sense=0,
        n_pool_intronic_strand_antisense=0,
        n_pool_intronic_strand_ambig=0,
        n_pool_dropped_out_of_range=0,
        exon_fit_tolerance_bp=int(exon_fit_tolerance_bp),
        fl_prior_ess=float(fl_prior_ess),
        extra=dict(extra),
    )
