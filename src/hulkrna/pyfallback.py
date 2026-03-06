"""hulkrna.pyfallback — Python fallback implementations for debugging.

Contains Python-only code paths that duplicate the C++ hot paths.
These exist for validation and debugging (``HULK_FORCE_PYTHON=1``)
and will be removed once C++ is established as the gold standard.

Dispatch points in the main modules delegate here when the Python
fallback path is requested.  Do NOT import from this module directly
in production code.
"""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

import numpy as np

from ._em_impl import (
    EM_LOG_EPSILON as _EM_LOG_EPSILON,
    run_locus_em_native as _run_locus_em_native,
)
from .annotate import (
    POOL_CODE_CHIMERIC,
    POOL_CODE_GDNA,
    POOL_CODE_MRNA,
    POOL_CODE_NRNA,
)
from .buffer import (
    FRAG_AMBIG_SAME_STRAND,
    FRAG_CHIMERIC,
    FRAG_MULTIMAPPER,
    FRAG_UNAMBIG,
)
from .splice import SPLICE_ANNOT, SPLICE_UNSPLICED

if TYPE_CHECKING:
    from .buffer import FragmentBuffer
    from .estimator import (
        AbundanceEstimator,
        LocusEMInput,
        ScoredFragments,
    )

logger = logging.getLogger(__name__)

# Convergence tolerance — mirrors the default in estimator.py but is
# defined here to avoid importing from estimator at module level (which
# would create a circular-import risk).
_EM_CONVERGENCE_DELTA = 1e-6


# ======================================================================
# Scan fallback
# ======================================================================


def scan_python(
    router,
    buffer: FragmentBuffer,
    log_every: int,
) -> ScoredFragments:
    """Per-fragment Python scan path (fallback when C++ unavailable).

    Reimplements ``FragmentRouter._scan_native`` in pure Python.
    Moved here from ``scan.py`` to centralise all fallback code.

    Parameters
    ----------
    router : FragmentRouter
        The ``FragmentRouter`` instance whose state (accumulators,
        estimator, stats, strand_models, etc.) is read and mutated.
    buffer : FragmentBuffer
        Fragment buffer to scan.
    log_every : int
        Log progress every *log_every* fragments.

    Returns
    -------
    ScoredFragments
    """
    # Late import to avoid circular dependency at module level.
    from .estimator import ScoredFragments

    estimator = router.estimator
    stats = router.stats
    index = router.index
    t_to_g = router.ctx.t_to_g
    annotations = router.annotations

    n_processed = 0
    for chunk in buffer.iter_chunks():
        frag_classes = chunk.fragment_classes

        for i in range(chunk.size):
            fc = frag_classes[i]

            # --- Chimeric ---
            if fc == FRAG_CHIMERIC:
                if annotations is not None:
                    chimeric_fid = int(chunk.frag_id[i])
                    annotations.add(
                        frag_id=chimeric_fid,
                        best_tid=-1,
                        best_gid=-1,
                        pool=POOL_CODE_CHIMERIC,
                        posterior=0.0,
                        frag_class=FRAG_CHIMERIC,
                        n_candidates=0,
                        splice_type=int(chunk[i].splice_type),
                    )
                n_processed += 1
                if n_processed % log_every == 0:
                    logger.debug(
                        f"  Scan: {n_processed:,} / "
                        f"{buffer.total_fragments:,}"
                    )
                continue

            bf = chunk[i]
            is_spliced_annot = bf.splice_type == SPLICE_ANNOT

            # --- Pre-EM strand/intronic accumulation ---
            if fc == FRAG_UNAMBIG or fc == FRAG_AMBIG_SAME_STRAND:
                first_t = int(bf.t_inds[0])
                t_strand = int(index.t_to_strand_arr[first_t])

                is_anti = estimator.is_antisense(
                    bf.exon_strand, t_strand,
                    router.strand_models.exonic_spliced,
                )

                is_unspliced = bf.splice_type == SPLICE_UNSPLICED

                n_cand = len(bf.t_inds)
                weight = 1.0 / n_cand
                for k, t_idx in enumerate(bf.t_inds):
                    t_idx_int = int(t_idx)

                    has_unambig_intron = (
                        bf.unambig_intron_bp is not None
                        and bf.unambig_intron_bp[k] > 0
                    )

                    if not has_unambig_intron:
                        if is_anti:
                            estimator.transcript_exonic_antisense[
                                t_idx_int
                            ] += weight
                        else:
                            estimator.transcript_exonic_sense[
                                t_idx_int
                            ] += weight

                    if is_unspliced:
                        if is_anti:
                            estimator.transcript_unspliced_antisense[
                                t_idx_int
                            ] += weight
                        else:
                            estimator.transcript_unspliced_sense[
                                t_idx_int
                            ] += weight

                    if has_unambig_intron:
                        if is_anti:
                            estimator.transcript_intronic_antisense[
                                t_idx_int
                            ] += weight
                        else:
                            estimator.transcript_intronic_sense[
                                t_idx_int
                            ] += weight

            # --- Deterministic unique: SPLICED_ANNOT + FRAG_UNAMBIG ---
            if fc == FRAG_UNAMBIG and is_spliced_annot:
                estimator.assign_unambig(bf, index, router.strand_models)
                stats.deterministic_unambig_units += 1

                if annotations is not None:
                    det_tid = int(next(iter(bf.t_inds)))
                    det_gid = int(t_to_g[det_tid])
                    annotations.add(
                        frag_id=int(chunk.frag_id[i]),
                        best_tid=det_tid,
                        best_gid=det_gid,
                        pool=POOL_CODE_MRNA,
                        posterior=1.0,
                        frag_class=FRAG_UNAMBIG,
                        n_candidates=1,
                        splice_type=int(bf.splice_type),
                    )

            elif fc == FRAG_MULTIMAPPER:
                fid = int(chunk.frag_id[i])
                if fid != router.mm_fid:
                    if router.mm_pending:
                        router._flush_mm_group()
                        stats.em_routed_multimapper_units += 1
                    router.mm_fid = fid
                    router.mm_pending.clear()
                    router.mm_pending.append(bf)
                else:
                    router.mm_pending.append(bf)

            else:
                router._add_single_fragment(bf, chunk, i, fc)

            n_processed += 1
            if n_processed % log_every == 0:
                logger.debug(
                    f"  Scan: {n_processed:,} / "
                    f"{buffer.total_fragments:,}"
                )

    # Flush last multimapper group
    if router.mm_pending:
        router._flush_mm_group()
        stats.em_routed_multimapper_units += 1

    n_units = len(router.offsets) - 1
    n_candidates = len(router.t_indices_list)

    def _to_np(arr, dtype):
        if len(arr) == 0:
            return np.empty(0, dtype=dtype)
        return np.frombuffer(arr, dtype=dtype).copy()

    return ScoredFragments(
        offsets=_to_np(router.offsets, np.int64),
        t_indices=_to_np(router.t_indices_list, np.int32),
        log_liks=_to_np(router.log_liks_list, np.float64),
        count_cols=_to_np(router.count_cols_list, np.uint8),
        coverage_weights=_to_np(
            router.coverage_weights_list, np.float64,
        ),
        tx_starts=_to_np(router.tx_starts_list, np.int32),
        tx_ends=_to_np(router.tx_ends_list, np.int32),
        locus_t_indices=_to_np(router.locus_t_list, np.int32),
        locus_count_cols=_to_np(router.locus_ct_list, np.uint8),
        is_spliced=_to_np(
            router.is_spliced_list, np.int8,
        ).astype(bool),
        gdna_log_liks=_to_np(router.gdna_ll_list, np.float64),
        genomic_footprints=_to_np(
            router.genomic_footprints_list, np.int32,
        ),
        frag_ids=_to_np(router.frag_id_list, np.int64),
        frag_class=_to_np(router.frag_class_list, np.int8),
        splice_type=_to_np(router.splice_type_list, np.uint8),
        n_units=n_units,
        n_candidates=n_candidates,
        nrna_base_index=router.ctx.nrna_base,
    )


# ======================================================================
# Per-locus EM solver fallback
# ======================================================================


def run_locus_em(
    estimator: AbundanceEstimator,
    locus_em: LocusEMInput,
    *,
    em_iterations: int = 1000,
    em_convergence_delta: float = _EM_CONVERGENCE_DELTA,
) -> tuple[np.ndarray, np.ndarray]:
    """Run EM for a single locus sub-problem.

    Delegates to the per-locus C++ solver (``run_locus_em_native``).
    This wrapper is part of the fallback path because the *batch*
    C++ path (``batch_locus_em``) replaces the entire Python per-locus
    loop.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (theta, alpha) — converged parameters for this locus.
    """
    theta, alpha_out, _em_totals, _nrna_frac_out = _run_locus_em_native(
        locus_em.offsets,
        locus_em.t_indices,
        locus_em.log_liks,
        locus_em.coverage_weights,
        locus_em.tx_starts,
        locus_em.tx_ends,
        locus_em.bias_profiles,
        locus_em.unambig_totals,
        locus_em.effective_lengths,
        (locus_em.prior > 0).astype(np.float64),
        locus_em.n_components,
        estimator.em_config.prior_alpha,
        estimator.em_config.prior_gamma,
        em_iterations,
        em_convergence_delta,
        estimator.em_config.mode == "vbem",
        (
            estimator.em_config.prune_threshold
            if estimator.em_config.prune_threshold is not None
            else -1.0
        ),
        locus_em.n_transcripts,
        locus_em.nrna_frac_alpha,
        locus_em.nrna_frac_beta,
    )
    return np.asarray(theta), np.asarray(alpha_out)


# ======================================================================
# Per-locus posterior assignment fallback
# ======================================================================


def assign_locus_ambiguous(
    estimator: AbundanceEstimator,
    locus_em: LocusEMInput,
    theta: np.ndarray,
    *,
    confidence_threshold: float = 0.95,
    unit_annotations: list | None = None,
) -> dict[str, float]:
    """Assign ambiguous units within one locus, scatter to global arrays.

    Each unit is assigned fractionally across all components (mRNA,
    nRNA, gDNA) using the converged posterior probabilities.

    Returns
    -------
    dict[str, float]
        Per-pool fragment counts: ``{"mrna": ..., "nrna": ..., "gdna": ...}``.
    """
    offsets = locus_em.offsets
    t_indices = locus_em.t_indices
    log_liks = locus_em.log_liks
    count_cols_arr = locus_em.count_cols
    n_units = len(offsets) - 1

    if n_units == 0 or len(t_indices) == 0:
        return {"mrna": 0.0, "nrna": 0.0, "gdna": 0.0}

    n_t = locus_em.n_transcripts
    gdna_idx = 2 * n_t  # the single gDNA component index
    local_to_global_t = locus_em.local_to_global_t
    seg_lengths = np.diff(offsets).astype(np.intp)

    # Compute posteriors from converged theta
    eff_len = locus_em.effective_lengths
    log_weights = np.log(theta + _EM_LOG_EPSILON) - np.log(eff_len)
    log_posteriors = log_liks + log_weights[t_indices]

    offsets_int = offsets[:-1].astype(np.intp)
    seg_max = np.maximum.reduceat(log_posteriors, offsets_int)
    log_posteriors -= np.repeat(seg_max, seg_lengths)
    posteriors = np.exp(log_posteriors)
    seg_sum = np.add.reduceat(posteriors, offsets_int)
    # Guard: segments where every candidate is -inf → zero posteriors
    bad_seg = (seg_sum == 0) | ~np.isfinite(seg_sum)
    seg_sum[bad_seg] = 1.0
    posteriors /= np.repeat(seg_sum, seg_lengths)
    if bad_seg.any():
        bad_mask = np.repeat(bad_seg, seg_lengths)
        posteriors[bad_mask] = 0.0

    # --- Classify all candidates at once ---
    mrna_mask = t_indices < n_t
    nrna_mask = (t_indices >= n_t) & (t_indices < gdna_idx)
    gdna_mask = t_indices == gdna_idx

    # --- Vectorized mRNA assignment ---
    mrna_count = 0.0
    if mrna_mask.any():
        mrna_p = posteriors[mrna_mask]
        mrna_count = float(mrna_p.sum())
        mrna_local_t = t_indices[mrna_mask]
        mrna_global_t = local_to_global_t[mrna_local_t]
        mrna_cols = count_cols_arr[mrna_mask]

        np.add.at(
            estimator.em_counts,
            (mrna_global_t, mrna_cols),
            mrna_p,
        )

        # High-confidence: per-unit max of mRNA posteriors
        mrna_p_for_max = np.where(mrna_mask, posteriors, 0.0)
        unit_max_mrna = np.maximum.reduceat(
            mrna_p_for_max, offsets_int,
        )
        conf_units = unit_max_mrna >= confidence_threshold
        if conf_units.any():
            conf_expanded = np.repeat(conf_units, seg_lengths)
            hc_mask = mrna_mask & conf_expanded
            if hc_mask.any():
                hc_p = posteriors[hc_mask]
                hc_local_t = t_indices[hc_mask]
                hc_global_t = local_to_global_t[hc_local_t]
                hc_cols = count_cols_arr[hc_mask]
                np.add.at(
                    estimator.em_high_conf_counts,
                    (hc_global_t, hc_cols),
                    hc_p,
                )

        # Confidence tracking
        if estimator._em_posterior_sum is None:
            estimator._em_posterior_sum = np.zeros(
                estimator.num_transcripts, dtype=np.float64,
            )
            estimator._em_n_assigned = np.zeros(
                estimator.num_transcripts, dtype=np.float64,
            )
        np.add.at(
            estimator._em_posterior_sum,
            mrna_global_t, mrna_p * mrna_p,
        )
        np.add.at(estimator._em_n_assigned, mrna_global_t, mrna_p)

    # --- Vectorized nRNA assignment ---
    nrna_count = 0.0
    if nrna_mask.any():
        nrna_p = posteriors[nrna_mask]
        nrna_count = float(nrna_p.sum())
        nrna_local_t = t_indices[nrna_mask] - n_t
        nrna_global_t = local_to_global_t[nrna_local_t]
        np.add.at(estimator.nrna_em_counts, nrna_global_t, nrna_p)

    # --- Vectorized gDNA assignment ---
    gdna_count = 0.0
    if gdna_mask.any():
        gdna_p = posteriors[gdna_mask]
        gdna_count = float(gdna_p.sum())

        # Per-unit gDNA sum for locus attribution
        gdna_p_full = np.where(gdna_mask, posteriors, 0.0)
        gdna_per_unit = np.add.reduceat(
            gdna_p_full, offsets_int,
        )
        # Scatter to gdna_locus_counts using per-unit locus_t/locus_ct
        locus_t_arr = locus_em.locus_t_indices
        locus_ct_arr = locus_em.locus_count_cols
        valid_gdna_units = (gdna_per_unit > 0) & (locus_t_arr >= 0)
        if valid_gdna_units.any():
            np.add.at(
                estimator.gdna_locus_counts,
                (
                    locus_t_arr[valid_gdna_units],
                    locus_ct_arr[valid_gdna_units],
                ),
                gdna_per_unit[valid_gdna_units],
            )

    # --- Per-unit annotation extraction (opt-in) ---
    if unit_annotations is not None:
        global_unit_indices = locus_em.locus.unit_indices
        t_to_g = estimator._t_to_g
        for u in range(n_units):
            s = int(offsets[u])
            e = int(offsets[u + 1])
            if e <= s:
                continue
            seg_posteriors = posteriors[s:e]
            seg_t_indices = t_indices[s:e]
            best_local = int(np.argmax(seg_posteriors))
            best_posterior = float(seg_posteriors[best_local])
            best_comp = int(seg_t_indices[best_local])

            if best_comp < n_t:
                g_tid = int(local_to_global_t[best_comp])
                g_gid = int(t_to_g[g_tid])
                pool_code = POOL_CODE_MRNA
            elif best_comp < gdna_idx:
                g_tid = int(local_to_global_t[best_comp - n_t])
                g_gid = int(t_to_g[g_tid])
                pool_code = POOL_CODE_NRNA
            else:
                g_tid = -1
                g_gid = -1
                pool_code = POOL_CODE_GDNA
            unit_annotations.append((
                int(global_unit_indices[u]),
                g_tid,
                g_gid,
                pool_code,
                best_posterior,
                e - s,  # n_candidates
            ))

    return {"mrna": mrna_count, "nrna": nrna_count, "gdna": gdna_count}


# ======================================================================
# Per-locus EM loop fallback (pipeline.py Python branch)
# ======================================================================


def run_per_locus_em_loop(
    *,
    estimator: AbundanceEstimator,
    loci: list,
    em_data: ScoredFragments,
    index,
    geometry,
    gdna_inits: np.ndarray,
    em_config,
    annotations,
    build_locus_meta,
) -> float:
    """Python per-locus EM loop — fallback for the batch C++ path.

    Iterates over every locus, calls ``build_locus_em_data`` →
    ``run_locus_em`` → ``assign_locus_ambiguous``, and maps unit
    annotations back to fragment IDs.

    Returns
    -------
    float
        Total gDNA EM count across all loci.
    """
    from .locus import build_locus_em_data

    total_gdna_em = 0.0
    _unit_ann_list: list | None = [] if annotations is not None else None
    _LARGE_LOCUS_LOG_THRESHOLD = 100

    # Pre-extract DataFrame columns for build_locus_em_data
    _locus_cache: dict = {
        "t_starts": index.t_df["start"].values,
        "t_ends": index.t_df["end"].values,
        "t_lengths": index.t_df["length"].values,
        "local_map": np.full(
            index.num_transcripts * 2 + 1, -1, dtype=np.int32,
        ),
    }

    for i, locus in enumerate(loci):
        if len(locus.unit_indices) == 0:
            continue

        locus_em = build_locus_em_data(
            locus, em_data, estimator, index, geometry.mean_frag,
            gdna_init=gdna_inits[i],
            _cache=_locus_cache,
        )
        theta, _alpha = run_locus_em(
            estimator,
            locus_em,
            em_iterations=em_config.iterations,
            em_convergence_delta=em_config.convergence_delta,
        )
        pool_counts = assign_locus_ambiguous(
            estimator,
            locus_em,
            theta,
            confidence_threshold=em_config.confidence_threshold,
            unit_annotations=_unit_ann_list,
        )
        gdna_assigned = pool_counts["gdna"]
        total_gdna_em += gdna_assigned

        estimator.locus_results.append(build_locus_meta(
            locus,
            mrna=pool_counts["mrna"],
            nrna=pool_counts["nrna"],
            gdna=gdna_assigned,
            gdna_init=gdna_inits[i],
        ))

        if len(locus.transcript_indices) > _LARGE_LOCUS_LOG_THRESHOLD:
            logger.debug(
                f"  Locus {locus.locus_id}: "
                f"{len(locus.transcript_indices)} transcripts, "
                f"{len(locus.unit_indices)} units, "
                f"gDNA={gdna_assigned:.0f}"
            )

    # --- Map EM unit annotations to frag_ids ---
    if _unit_ann_list:
        for (
            global_unit_idx, g_tid, g_gid,
            pool_code, posterior, n_cand,
        ) in _unit_ann_list:
            fid = int(em_data.frag_ids[global_unit_idx])
            fc = int(em_data.frag_class[global_unit_idx])
            st = int(em_data.splice_type[global_unit_idx])
            annotations.add(
                frag_id=fid,
                best_tid=g_tid,
                best_gid=g_gid,
                pool=pool_code,
                posterior=posterior,
                frag_class=fc,
                n_candidates=n_cand,
                splice_type=st,
            )

    return total_gdna_em
