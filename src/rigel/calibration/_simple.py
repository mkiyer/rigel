"""SRD v1 calibration orchestrator (the new ``calibrate_gdna``).

Single buffer walk → categorize → pool assembly → 1-D mixture fit →
EB-smoothed FL models. Library-agnostic, no regional priors, no
density modelling, no overlap clusters.

This is the Phase 2 staging entry point. It runs parallel to v5's
``_calibrate.calibrate_gdna``; Phase 3 wires it into ``pipeline.py``.
"""

from __future__ import annotations

import logging

import numpy as np

from ..buffer import FragmentBuffer
from ..frag_length_model import FragmentLengthModels
from ..index import TranscriptIndex
from ._categorize import (
    FragmentCategory,
    N_CATEGORIES,
    build_t_to_ref_id,
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
    exon_fit_tolerance_bp: int = 5,
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
    t_to_strand = np.asarray(index.t_to_strand_arr, dtype=np.int8)
    t_to_ref_id = build_t_to_ref_id(index.t_df)
    max_size = int(frag_length_models.max_size)
    n_bins = max_size + 1

    # ---- Pass 0: walk chunks, apply unique-map filter, categorize -------
    cat_chunks: list[np.ndarray] = []
    len_chunks: list[np.ndarray] = []
    n_multimap_excluded = 0

    for chunk in buffer.iter_chunks():
        unique_mask = chunk.num_hits == 1
        n_multimap_excluded += int((~unique_mask).sum())
        if not unique_mask.any():
            continue

        cc = categorize_chunk(
            chunk,
            t_to_strand=t_to_strand,
            t_to_ref_id=t_to_ref_id,
            exon_fit_tolerance_bp=exon_fit_tolerance_bp,
        )
        cat_chunks.append(cc.category[unique_mask])
        len_chunks.append(cc.frag_length[unique_mask])

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
    frag_length = np.concatenate(len_chunks)
    cat_counts = np.bincount(category, minlength=N_CATEGORIES).astype(np.int64)

    # ---- Pass 1: pool assembly -----------------------------------------
    # Pool = EXON_INCOMPATIBLE ∪ INTRONIC ∪ INTERGENIC.
    # Library-agnostic; no SS dependence. UNSPLICED_ANTISENSE_EXONIC is
    # tracked as a diagnostic (n_antisense_exonic) but excluded from the
    # pool: at imperfect SS it is dominated by RNA from highly-expressed
    # loci and would dilute the gDNA_FL signal we are trying to recover.
    pool_mask = (
        (category == int(FragmentCategory.EXON_INCOMPATIBLE))
        | (category == int(FragmentCategory.INTRONIC))
        | (category == int(FragmentCategory.INTERGENIC))
    )
    pool_lens = frag_length[pool_mask]
    n_pool_categorized = int(pool_mask.sum())
    pool_hist = np.bincount(
        np.clip(pool_lens, 0, max_size).astype(np.intp),
        minlength=n_bins,
    ).astype(np.float64)

    # Pool B contribution from the C++ scanner's intergenic FL accumulator.
    # Unique-mapper unspliced fragments that resolve to ZERO transcript
    # candidates are dropped before reaching the buffer (resolve_context.h:1043),
    # so they never reach `categorize_chunk`. The scanner sets them aside
    # in `frag_length_models.intergenic` (see bam_scanner.cpp:1340-1344 and
    # pipeline._replay_fraglen_observations). They are pure-by-construction
    # gDNA/extra-genic candidates and belong in Pool B.
    intergenic_counts = np.zeros(n_bins, dtype=np.float64)
    n_intergenic_unique = 0
    intergenic_model = getattr(frag_length_models, "intergenic", None)
    if intergenic_model is not None:
        ic = np.asarray(intergenic_model.counts, dtype=np.float64)
        if ic.size != n_bins:
            adj = np.zeros(n_bins, dtype=np.float64)
            m = min(ic.size, n_bins)
            adj[:m] = ic[:m]
            ic = adj
        intergenic_counts = ic
        n_intergenic_unique = int(ic.sum())
        pool_hist = pool_hist + ic

    n_pool = n_pool_categorized + n_intergenic_unique

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
    n_spliced = int(cat_counts[int(FragmentCategory.SPLICED)])
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
        "[CAL-SRD] quality=%s SS=%.3f n_pool=%d (categorized=%d + intergenic=%d) "
        "pi=%.4f conv=%s n_iter=%d categories=%s n_multimap_excluded=%d",
        quality,
        float(strand_specificity),
        n_pool,
        n_pool_categorized,
        n_intergenic_unique,
        pi,
        mixture_converged,
        mixture_iterations,
        cat_counts.tolist(),
        n_multimap_excluded,
    )

    return CalibrationResult(
        gdna_fl_model=gdna_fl,
        rna_fl_model=rna_fl,
        global_fl_model=global_fl,
        gdna_fl_quality=quality,
        strand_specificity=float(strand_specificity),
        category_counts=cat_counts,
        n_multimap_excluded=n_multimap_excluded,
        n_pool=n_pool,
        pi_pool=pi,
        mixture_converged=mixture_converged,
        mixture_iterations=mixture_iterations,
        n_intergenic_unique=n_intergenic_unique,
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
        category_counts=np.zeros(N_CATEGORIES, dtype=np.int64),
        n_multimap_excluded=int(n_multimap_excluded),
        n_pool=0,
        pi_pool=0.0,
        mixture_converged=False,
        mixture_iterations=0,
        n_intergenic_unique=0,
        exon_fit_tolerance_bp=int(exon_fit_tolerance_bp),
        fl_prior_ess=float(fl_prior_ess),
        extra=dict(extra),
    )
