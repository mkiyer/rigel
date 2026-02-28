"""hulkrna.scan — CSR builder / buffer scanner.

Converts a ``FragmentBuffer`` into ``ScanData`` for the EM solver.
The ``EmDataBuilder`` class replaces the old 620-line
``_scan_and_build_em_data`` closure-based function.

**Bug fixes vs. the old triplicated code:**

1. **MM nRNA gate was too strict** — the old ``_flush_mm_group`` used
   ``any_spliced`` (ANNOT or UNANNOT) to gate the entire nRNA pool.
   The single-fragment path only skipped nRNA for ``SPLICED_ANNOT``.
   Now both paths use ``is_annot_spliced`` (ANNOT only), so
   ``SPLICED_UNANNOT`` multimapper groups correctly get nRNA candidates.

2. **MM splice_type preferred first-found over most-informative** —
   the old code took the splice_type of the first spliced hit it found,
   potentially choosing UNANNOT when an ANNOT hit existed later.  Now
   we always prefer ``SPLICED_ANNOT > SPLICED_UNANNOT > UNSPLICED``.

3. **Dead parameters removed** — ``_finalize_unit`` had unused ``best_t``
   / ``best_ct`` parameters.  Eliminated.
"""

import math

import numpy as np

from .buffer import (
    FragmentBuffer,
    FRAG_UNIQUE,
    FRAG_ISOFORM_AMBIG,
    FRAG_MULTIMAPPER,
    FRAG_CHIMERIC,
)
from .estimator import AbundanceEstimator, ScanData
from .index import HulkIndex
from .scoring import (
    ScoringContext,
    LOG_SAFE_FLOOR,
    LOG_HALF,
    SPLICE_UNSPLICED,
    SPLICE_UNANNOT,
    SPLICE_ANNOT,
    STRAND_NEG,
    frag_len_log_lik,
    frag_len_log_lik_batch,
    genomic_to_transcript_pos,
    genomic_to_transcript_pos_bisect,
    compute_coverage_weight,
    compute_coverage_weight_batch,
)
from .stats import PipelineStats


class EmDataBuilder:
    """Builds the global CSR (ScanData) from a FragmentBuffer.

    State that was previously captured via closures is now explicit
    instance attributes.  Every method is a candidate for Cython / C
    acceleration.

    Parameters
    ----------
    ctx : ScoringContext
        Pre-computed scoring parameters.
    counter : AbundanceEstimator
        Accumulates unique counts during scan.
    stats : PipelineStats
        Pipeline statistics accumulator.
    index : HulkIndex
        Reference index.
    annotations : AnnotationTable or None
        If provided, deterministic-unique and chimeric fragments are
        annotated immediately during the scan.
    """

    def __init__(
        self,
        ctx: ScoringContext,
        counter: AbundanceEstimator,
        stats: PipelineStats,
        index: HulkIndex,
        strand_models,
        annotations=None,
    ):
        self.ctx = ctx
        self.counter = counter
        self.stats = stats
        self.index = index
        self.strand_models = strand_models
        self.annotations = annotations

        # CSR accumulation arrays (array.array for zero-copy → numpy)
        from array import array as _array
        self.offsets = _array('q', [0])          # int64
        self.t_indices_list = _array('i')         # int32
        self.log_liks_list = _array('d')          # float64
        self.count_cols_list = _array('B')        # uint8
        self.coverage_weights_list = _array('d')  # float64
        self.tx_starts_list = _array('i')         # int32
        self.tx_ends_list = _array('i')           # int32

        # Per-unit metadata arrays
        self.locus_t_list = _array('i')           # int32
        self.locus_ct_list = _array('B')          # uint8
        self.is_spliced_list = _array('b')        # int8 (→ bool at finalize)
        self.gdna_ll_list = _array('d')           # float64
        self.genomic_footprints_list = _array('i')  # int32
        self.frag_id_list = _array('q')           # int64
        self.frag_class_list = _array('b')        # int8
        self.splice_type_list = _array('B')       # uint8

        # Multimapper group state
        self.mm_pending: list = []
        self.mm_fid: int = -1

    # ------------------------------------------------------------------
    # Winner-take-all mRNA scoring (shared by single-fragment and MM)
    # ------------------------------------------------------------------

    def _score_wta_mrna(
        self,
        t_inds,
        exon_bp,
        frag_lengths,
        exon_strand: int,
        splice_type: int,
        nm: int,
        read_length: int,
        genomic_start: int = -1,
    ) -> dict[int, tuple[int, float, int, float, int, int]]:
        """Score mRNA candidates with winner-take-all gating.

        Returns dict mapping t_idx → (overhang, log_lik, count_col,
        coverage_weight, tx_start, tx_end).  Only WTA winners (minimum
        overhang) are included.

        Vectorized implementation: replaces per-candidate Python loops
        with NumPy array operations for overhangs, strand log-probs,
        fragment-length log-liks, and coverage weights.
        """
        ctx = self.ctx
        n_cand = len(t_inds)
        if n_cand == 0:
            return {}

        fl = read_length if read_length > 0 else 1
        log_nm = nm * ctx.mismatch_log_penalty if nm > 0 else 0.0

        # --- Vectorized arrays ---
        t_inds_arr = np.asarray(t_inds, dtype=np.int32)

        # Exon BP and validity
        if exon_bp is not None:
            ebp_arr = np.asarray(exon_bp, dtype=np.int32)
        else:
            ebp_arr = np.full(n_cand, fl, dtype=np.int32)
        valid = ebp_arr > 0
        if not np.any(valid):
            return {}

        # Overhangs (vectorized)
        oh_arr = np.maximum(fl - ebp_arr, 0)

        # Fragment lengths
        if frag_lengths is not None:
            fl_arr = np.asarray(frag_lengths, dtype=np.int32)
        else:
            fl_arr = np.full(n_cand, -1, dtype=np.int32)

        # Fragment-length log-likelihoods (vectorized LUT)
        log_fl_arr = frag_len_log_lik_batch(ctx, fl_arr)

        # Strand log-likelihoods (vectorized)
        t_strands = ctx.t_strand_arr[t_inds_arr]
        if exon_strand == 1 or exon_strand == 2:
            same_mask = (exon_strand == t_strands)
            log_strand = np.where(
                same_mask, ctx.log_p_sense, ctx.log_p_antisense,
            )
            anti_arr = np.where(same_mask, ctx.anti_flag, not ctx.anti_flag)
        else:
            log_strand = np.full(n_cand, LOG_HALF)
            anti_arr = np.zeros(n_cand, dtype=bool)

        # Combined log-likelihood (vectorized)
        log_lik_arr = (
            log_strand + log_fl_arr
            + oh_arr * ctx.overhang_log_penalty + log_nm
        )

        # Count columns (vectorized)
        ct_arr = splice_type * 2 + anti_arr.astype(np.int32)

        # --- Coverage weights and transcript positions ---
        cov_wt_arr = np.ones(n_cand, dtype=np.float64)
        tx_s_arr = np.zeros(n_cand, dtype=np.int32)
        t_lens = ctx.t_length_arr[t_inds_arr]
        tx_e_arr = np.where(fl_arr > 0, fl_arr, t_lens).astype(np.int32)

        if genomic_start >= 0:
            need_pos = valid & (fl_arr > 0)
            need_idx = np.where(need_pos)[0]
            if len(need_idx) > 0:
                t_exon_data = ctx._t_exon_data
                # Compute transcript positions using bisect on pre-computed data
                n_need = len(need_idx)
                gs_starts = np.empty(n_need, dtype=np.int32)
                gs_ends = np.empty(n_need, dtype=np.int32)
                gs_t_lens = np.empty(n_need, dtype=np.int32)
                cov_valid_mask = np.ones(n_need, dtype=bool)

                for j in range(n_need):
                    k = need_idx[j]
                    t_idx_int = int(t_inds_arr[k])
                    t_len = int(t_lens[k])
                    flen_k = int(fl_arr[k])
                    exon_data = t_exon_data.get(t_idx_int) if t_exon_data else None
                    if exon_data is not None and t_len > 0:
                        starts, ends, cumsum_before = exon_data
                        t_strand = int(t_strands[k])
                        tx_s = genomic_to_transcript_pos_bisect(
                            genomic_start, starts, ends,
                            cumsum_before, t_strand, t_len,
                        )
                        gs_starts[j] = tx_s
                        gs_ends[j] = tx_s + flen_k
                        gs_t_lens[j] = t_len
                    else:
                        gs_starts[j] = 0
                        gs_ends[j] = flen_k
                        gs_t_lens[j] = max(flen_k, 1)
                        cov_valid_mask[j] = False

                # Batch coverage weight (vectorized)
                cov_ends_clipped = np.minimum(gs_ends, gs_t_lens)
                cov_wts = compute_coverage_weight_batch(
                    gs_starts, cov_ends_clipped, gs_t_lens,
                )
                # Only apply for candidates with valid exon data
                cov_wts[~cov_valid_mask] = 1.0

                tx_s_arr[need_idx] = gs_starts
                tx_e_arr[need_idx] = gs_ends  # unclipped for bias model
                cov_wt_arr[need_idx] = cov_wts

        # --- WTA gating: keep candidates with minimum overhang ---
        # Exclude invalid candidates from the min computation
        oh_for_wta = oh_arr.copy()
        oh_for_wta[~valid] = np.iinfo(np.int32).max
        min_oh = int(oh_for_wta.min())
        winners_mask = valid & (oh_arr == min_oh)

        result: dict[int, tuple[int, float, int, float, int, int]] = {}
        for k in np.where(winners_mask)[0]:
            result[int(t_inds_arr[k])] = (
                int(oh_arr[k]),
                float(log_lik_arr[k]),
                int(ct_arr[k]),
                float(cov_wt_arr[k]),
                int(tx_s_arr[k]),
                int(tx_e_arr[k]),
            )
        return result

    # ------------------------------------------------------------------
    # Winner-take-all nRNA scoring (shared by single-fragment and MM)
    # ------------------------------------------------------------------

    def _score_wta_nrna(
        self,
        t_inds,
        exon_bp,
        intron_bp,
        exon_strand: int,
        nm: int,
        read_length: int,
        genomic_footprint: int,
        genomic_start: int = -1,
    ) -> dict[int, tuple[int, float, float, int, int]]:
        """Score nRNA candidates with winner-take-all gating.

        Returns dict mapping t_idx → (overhang, nrna_log_lik, cov_wt,
        tx_start, tx_end).  Only WTA winners (minimum overhang) are
        included.

        Vectorized implementation: replaces per-candidate Python loops
        with NumPy array operations.
        """
        ctx = self.ctx
        n_cand = len(t_inds)
        if n_cand == 0:
            return {}

        fl = read_length if read_length > 0 else 1
        log_fl = frag_len_log_lik(ctx, genomic_footprint)
        log_nm = nm * ctx.mismatch_log_penalty if nm > 0 else 0.0

        # --- Vectorized arrays ---
        t_inds_arr = np.asarray(t_inds, dtype=np.int32)

        # Single-exon transcript filter (vectorized)
        t_spans = ctx.t_span_arr[t_inds_arr]
        t_exonics = ctx.t_length_arr[t_inds_arr]
        multi_exon = t_spans > t_exonics

        # Genic overlap
        if exon_bp is not None:
            ebp_arr = np.asarray(exon_bp, dtype=np.int32)
        else:
            ebp_arr = np.full(n_cand, fl, dtype=np.int32)
        if intron_bp is not None:
            ibp_arr = np.asarray(intron_bp, dtype=np.int32)
        else:
            ibp_arr = np.zeros(n_cand, dtype=np.int32)

        span_bp = ebp_arr + ibp_arr
        valid = multi_exon & (span_bp > 0)
        if not np.any(valid):
            return {}

        # Overhangs (vectorized)
        oh_arr = np.maximum(fl - span_bp, 0)

        # Strand log-likelihoods (vectorized)
        t_strands = ctx.t_strand_arr[t_inds_arr]
        if exon_strand == 1 or exon_strand == 2:
            same_mask = (exon_strand == t_strands)
            log_strand = np.where(
                same_mask, ctx.log_p_sense, ctx.log_p_antisense,
            )
        else:
            log_strand = np.full(n_cand, LOG_HALF)

        # Combined nRNA log-likelihood (vectorized)
        nrna_ll_arr = (
            log_strand + log_fl + oh_arr * ctx.overhang_log_penalty + log_nm
        )

        # --- Coverage weights (nRNA uses genomic span) ---
        cov_wt_arr = np.ones(n_cand, dtype=np.float64)
        nrna_tx_s_arr = np.zeros(n_cand, dtype=np.int32)
        nrna_gfp = genomic_footprint if genomic_footprint > 0 else 0
        nrna_tx_e_arr = np.where(
            nrna_gfp > 0,
            np.int32(nrna_gfp),
            t_spans,
        ).astype(np.int32)

        if genomic_start >= 0 and genomic_footprint > 0:
            need_idx = np.where(valid)[0]
            if len(need_idx) > 0:
                t_span_vals = t_spans[need_idx]
                tx_starts_g = ctx.t_start_arr[t_inds_arr[need_idx]]
                frag_pos = genomic_start - tx_starts_g
                frag_end_pos = frag_pos + genomic_footprint

                # Clamped for coverage weight (vectorized)
                frag_pos_c = np.clip(frag_pos, 0, t_span_vals)
                frag_end_c = np.clip(frag_end_pos, 0, t_span_vals)
                cov_valid = frag_end_c > frag_pos_c

                if np.any(cov_valid):
                    cov_sub = need_idx[cov_valid]
                    cov_wt_arr[cov_sub] = compute_coverage_weight_batch(
                        frag_pos_c[cov_valid],
                        frag_end_c[cov_valid],
                        t_span_vals[cov_valid],
                    )

                # Unclipped for bias model
                nrna_tx_s_arr[need_idx] = np.maximum(0, frag_pos)
                nrna_tx_e_arr[need_idx] = (
                    np.maximum(0, frag_pos) + genomic_footprint
                )

        # --- WTA gating ---
        oh_for_wta = oh_arr.copy()
        oh_for_wta[~valid] = np.iinfo(np.int32).max
        min_oh = int(oh_for_wta.min())
        winners_mask = valid & (oh_arr == min_oh)

        result: dict[int, tuple[int, float, float, int, int]] = {}
        for k in np.where(winners_mask)[0]:
            result[int(t_inds_arr[k])] = (
                int(oh_arr[k]),
                float(nrna_ll_arr[k]),
                float(cov_wt_arr[k]),
                int(nrna_tx_s_arr[k]),
                int(nrna_tx_e_arr[k]),
            )
        return result

    # ------------------------------------------------------------------
    # Emit candidates to CSR
    # ------------------------------------------------------------------

    def _emit_mrna(
        self, winners: dict[int, tuple[int, float, int, float, int, int]],
    ) -> tuple[float, int, int]:
        """Append mRNA winner candidates to CSR lists.

        Returns (best_log_lik, best_t_idx, best_count_col).
        """
        best_ll = -np.inf
        best_t = -1
        best_ct = 0
        for t_idx, (oh, log_lik, ct, cov_wt, tx_s, tx_e) in winners.items():
            self.t_indices_list.append(t_idx)
            self.log_liks_list.append(log_lik)
            self.count_cols_list.append(ct)
            self.coverage_weights_list.append(cov_wt)
            self.tx_starts_list.append(tx_s)
            self.tx_ends_list.append(tx_e)
            if log_lik > best_ll:
                best_ll = log_lik
                best_t = t_idx
                best_ct = ct
        return best_ll, best_t, best_ct

    def _emit_nrna(self, winners: dict[int, tuple[int, float, float, int, int]]) -> None:
        """Append nRNA winner candidates to CSR lists."""
        nrna_base = self.ctx.nrna_base
        for t_idx, (oh, nrna_ll, cov_wt, tx_s, tx_e) in winners.items():
            self.t_indices_list.append(nrna_base + t_idx)
            self.log_liks_list.append(nrna_ll)
            self.count_cols_list.append(0)
            self.coverage_weights_list.append(cov_wt)
            self.tx_starts_list.append(tx_s)
            self.tx_ends_list.append(tx_e)

    # ------------------------------------------------------------------
    # Per-unit metadata
    # ------------------------------------------------------------------

    def _finalize_unit_metadata(self, bf) -> None:
        """Record is_spliced and gDNA log-likelihood for this EM unit."""
        ctx = self.ctx
        st = bf.splice_type
        is_spl = (st == SPLICE_ANNOT or st == SPLICE_UNANNOT)
        self.is_spliced_list.append(is_spl)
        self.genomic_footprints_list.append(int(bf.genomic_footprint))

        if not is_spl:
            es = bf.exon_strand
            sp = ctx.gdna_splice_penalties.get(st, 1.0)
            if es == 1 or es == 2:
                p_ig = ctx.ig_p if es == 1 else (1.0 - ctx.ig_p)
            else:
                p_ig = 0.5
            log_s = math.log(max(p_ig, LOG_SAFE_FLOOR))
            log_i = frag_len_log_lik(ctx, bf.genomic_footprint)
            self.gdna_ll_list.append(
                log_s + log_i + math.log(max(sp, LOG_SAFE_FLOOR))
            )
        else:
            self.gdna_ll_list.append(-np.inf)

    # ------------------------------------------------------------------
    # Single-fragment EM unit
    # ------------------------------------------------------------------

    def _add_single_fragment(self, bf, chunk, i, fc) -> None:
        """Process a non-chimeric, non-multimapper EM unit."""
        is_spliced_annot = bf.splice_type == SPLICE_ANNOT

        # --- mRNA candidates (WTA) ---
        mrna_winners = self._score_wta_mrna(
            bf.t_inds, bf.exon_bp, bf.frag_lengths,
            bf.exon_strand, bf.splice_type, bf.nm, bf.read_length,
            bf.genomic_start,
        )
        best_ll, best_t, best_ct = self._emit_mrna(mrna_winners)

        # --- nRNA candidates (independent pool, skip only SPLICED_ANNOT) ---
        if not is_spliced_annot:
            nrna_winners = self._score_wta_nrna(
                bf.t_inds, bf.exon_bp, bf.intron_bp,
                bf.exon_strand, bf.nm, bf.read_length,
                bf.genomic_footprint, bf.genomic_start,
            )
            self._emit_nrna(nrna_winners)

        if len(self.t_indices_list) > self.offsets[-1]:
            self.offsets.append(len(self.t_indices_list))
            self.locus_t_list.append(best_t)
            self.locus_ct_list.append(best_ct)
            self._finalize_unit_metadata(bf)
            self.frag_id_list.append(int(chunk.frag_id[i]))
            self.frag_class_list.append(fc)
            self.splice_type_list.append(int(bf.splice_type))
        else:
            self.stats.n_gated_out += 1

        if fc == FRAG_UNIQUE:
            self.stats.em_routed_unique_units += 1
        elif fc == FRAG_ISOFORM_AMBIG:
            self.stats.em_routed_isoform_ambig_units += 1
        else:
            self.stats.em_routed_gene_ambig_units += 1

    # ------------------------------------------------------------------
    # Multimapper group flush
    # ------------------------------------------------------------------

    def _flush_mm_group(self) -> None:
        """Flush accumulated multimapper alignments as one EM unit.

        Uses winner-take-all gating with independent mRNA and nRNA pools.

        **Bug fixes:**
        1. nRNA pool gate now uses is_annot_spliced (ANNOT only),
           consistent with the single-fragment path.  Previously
           SPLICED_UNANNOT MM groups incorrectly lost all nRNA candidates.
        2. splice_type prefers ANNOT > UNANNOT > UNSPLICED (most
           informative), not first-found.
        """
        ctx = self.ctx
        mm_pending = self.mm_pending

        # --- Determine splice status ---
        # BUG FIX #1: Use is_annot_spliced (ANNOT only) for nRNA gate,
        # matching the single-fragment path.  The is_any_spliced flag
        # (ANNOT or UNANNOT) is still used for the is_spliced metadata
        # and gDNA exclusion, which is correct.
        is_any_spliced = False
        is_annot_spliced = False
        for bf in mm_pending:
            st = bf.splice_type
            if st == SPLICE_ANNOT:
                is_annot_spliced = True
                is_any_spliced = True
                break  # ANNOT is highest priority, no need to continue
            elif st == SPLICE_UNANNOT:
                is_any_spliced = True

        # ----- mRNA pool (WTA across all hits) -----
        merged_mrna: dict[int, tuple[int, float, int, float, int, int]] = {}
        for bf in mm_pending:
            winners = self._score_wta_mrna(
                bf.t_inds, bf.exon_bp, bf.frag_lengths,
                bf.exon_strand, bf.splice_type, bf.nm, bf.read_length,
                bf.genomic_start,
            )
            # Merge: keep best per transcript (lower oh wins, then higher ll)
            for t_idx, (oh, ll, ct, cov_wt, tx_s, tx_e) in winners.items():
                prev = merged_mrna.get(t_idx)
                if prev is None or oh < prev[0] or (
                    oh == prev[0] and ll > prev[1]
                ):
                    merged_mrna[t_idx] = (oh, ll, ct, cov_wt, tx_s, tx_e)

        # Global WTA across all transcripts
        if merged_mrna:
            min_oh = min(v[0] for v in merged_mrna.values())
            final_mrna = {
                t: v for t, v in merged_mrna.items() if v[0] == min_oh
            }
        else:
            final_mrna = {}

        best_ll, best_t, best_ct = self._emit_mrna(final_mrna)

        # ----- nRNA pool (independent, skip only SPLICED_ANNOT) -----
        # BUG FIX #1: gate on is_annot_spliced, not is_any_spliced
        if not is_annot_spliced:
            merged_nrna: dict[int, tuple[int, float, float, int, int]] = {}
            for bf in mm_pending:
                if bf.splice_type == SPLICE_ANNOT:
                    continue
                winners = self._score_wta_nrna(
                    bf.t_inds, bf.exon_bp, bf.intron_bp,
                    bf.exon_strand, bf.nm, bf.read_length,
                    bf.genomic_footprint, bf.genomic_start,
                )
                for t_idx, (oh, nrna_ll, cov_wt, tx_s, tx_e) in winners.items():
                    prev = merged_nrna.get(t_idx)
                    if prev is None or oh < prev[0] or (
                        oh == prev[0] and nrna_ll > prev[1]
                    ):
                        merged_nrna[t_idx] = (oh, nrna_ll, cov_wt, tx_s, tx_e)

            if merged_nrna:
                min_oh = min(v[0] for v in merged_nrna.values())
                final_nrna = {
                    t: v for t, v in merged_nrna.items()
                    if v[0] == min_oh
                }
            else:
                final_nrna = {}

            self._emit_nrna(final_nrna)

        # --- Record unit if any candidates were admitted ---
        if len(self.t_indices_list) > self.offsets[-1]:
            self.offsets.append(len(self.t_indices_list))
            self.locus_t_list.append(best_t)
            self.locus_ct_list.append(best_ct)
            self.frag_id_list.append(self.mm_fid)
            self.frag_class_list.append(FRAG_MULTIMAPPER)

            # BUG FIX #2: prefer most-informative splice type
            # ANNOT > UNANNOT > UNSPLICED
            mm_stype = SPLICE_UNSPLICED
            for bf in mm_pending:
                st = int(bf.splice_type)
                if st == SPLICE_ANNOT:
                    mm_stype = SPLICE_ANNOT
                    break
                elif st == SPLICE_UNANNOT:
                    mm_stype = SPLICE_UNANNOT
            self.splice_type_list.append(mm_stype)

            # Per-unit metadata
            self.is_spliced_list.append(is_any_spliced)

            if not is_any_spliced:
                # gDNA log-lik: best across unspliced MM hits
                best_gdna_ll = -np.inf
                best_footprint = 0
                for bf in mm_pending:
                    st = bf.splice_type
                    if st == SPLICE_ANNOT or st == SPLICE_UNANNOT:
                        continue
                    es = bf.exon_strand
                    sp = ctx.gdna_splice_penalties.get(st, 1.0)
                    if es == 1 or es == 2:
                        p_ig = ctx.ig_p if es == 1 else (1.0 - ctx.ig_p)
                    else:
                        p_ig = 0.5
                    log_s = math.log(max(p_ig, LOG_SAFE_FLOOR))
                    log_i = frag_len_log_lik(ctx, bf.genomic_footprint)
                    gdna_ll = (
                        log_s + log_i
                        + math.log(max(sp, LOG_SAFE_FLOOR))
                    )
                    if gdna_ll > best_gdna_ll:
                        best_gdna_ll = gdna_ll
                        best_footprint = int(bf.genomic_footprint)
                self.gdna_ll_list.append(best_gdna_ll)
                self.genomic_footprints_list.append(best_footprint)
            else:
                self.gdna_ll_list.append(-np.inf)
                # Use first hit's footprint as representative
                self.genomic_footprints_list.append(
                    int(mm_pending[0].genomic_footprint)
                )
        else:
            self.stats.n_gated_out += 1

    # ------------------------------------------------------------------
    # Main scan driver
    # ------------------------------------------------------------------

    def scan(self, buffer: FragmentBuffer, log_every: int) -> ScanData:
        """Single pass over buffer.  Returns packed ScanData.

        Deterministic-unique (SPLICED_ANNOT + FRAG_UNIQUE) fragments
        are counted directly via ``counter.assign_unique``.  All other
        exonic fragments build EM units.  Chimeric fragments are
        recorded in the annotation table (if active) and skipped.
        """
        import logging
        logger = logging.getLogger(__name__)

        counter = self.counter
        stats = self.stats
        index = self.index
        t_to_g = self.ctx.t_to_g
        annotations = self.annotations

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
                            pool=4,  # chimeric
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
                if fc == FRAG_UNIQUE or fc == FRAG_ISOFORM_AMBIG:
                    g_idx = int(t_to_g[int(bf.t_inds[0])])
                    gene_strand = int(index.g_to_strand_arr[g_idx])

                    is_anti = counter.is_antisense(
                        bf.exon_strand, gene_strand,
                        self.strand_models.exonic_spliced,
                    )

                    is_unspliced = (bf.splice_type == SPLICE_UNSPLICED)
                    if is_unspliced:
                        if is_anti:
                            counter.gene_antisense_all[g_idx] += 1.0
                        else:
                            counter.gene_sense_all[g_idx] += 1.0

                    # Transcript-level intronic (for nRNA init)
                    n_cand = len(bf.t_inds)
                    weight = 1.0 / n_cand
                    for k, t_idx in enumerate(bf.t_inds):
                        t_idx_int = int(t_idx)
                        has_intron = (
                            bf.unambig_intron_bp is not None
                            and bf.unambig_intron_bp[k] > 0
                        )
                        if has_intron:
                            if is_anti:
                                counter.transcript_intronic_antisense[
                                    t_idx_int
                                ] += weight
                            else:
                                counter.transcript_intronic_sense[
                                    t_idx_int
                                ] += weight

                # --- Deterministic unique: SPLICED_ANNOT + FRAG_UNIQUE ---
                if fc == FRAG_UNIQUE and is_spliced_annot:
                    counter.assign_unique(bf, index, self.strand_models)
                    stats.deterministic_unique_units += 1

                    if annotations is not None:
                        det_tid = int(next(iter(bf.t_inds)))
                        det_gid = int(t_to_g[det_tid])
                        annotations.add(
                            frag_id=int(chunk.frag_id[i]),
                            best_tid=det_tid,
                            best_gid=det_gid,
                            pool=0,  # mRNA
                            posterior=1.0,
                            frag_class=FRAG_UNIQUE,
                            n_candidates=1,
                            splice_type=int(bf.splice_type),
                        )

                elif fc == FRAG_MULTIMAPPER:
                    fid = int(chunk.frag_id[i])
                    if fid != self.mm_fid:
                        if self.mm_pending:
                            self._flush_mm_group()
                            stats.em_routed_multimapper_units += 1
                        self.mm_fid = fid
                        self.mm_pending.clear()
                        self.mm_pending.append(bf)
                    else:
                        self.mm_pending.append(bf)

                else:
                    self._add_single_fragment(bf, chunk, i, fc)

                n_processed += 1
                if n_processed % log_every == 0:
                    logger.debug(
                        f"  Scan: {n_processed:,} / "
                        f"{buffer.total_fragments:,}"
                    )

        # Flush last multimapper group
        if self.mm_pending:
            self._flush_mm_group()
            stats.em_routed_multimapper_units += 1

        n_units = len(self.offsets) - 1
        n_candidates = len(self.t_indices_list)

        # Zero-copy conversion: array.array → numpy via buffer protocol
        def _to_np(arr, dtype):
            if len(arr) == 0:
                return np.empty(0, dtype=dtype)
            return np.frombuffer(arr, dtype=dtype).copy()

        return ScanData(
            offsets=_to_np(self.offsets, np.int64),
            t_indices=_to_np(self.t_indices_list, np.int32),
            log_liks=_to_np(self.log_liks_list, np.float64),
            count_cols=_to_np(self.count_cols_list, np.uint8),
            coverage_weights=_to_np(
                self.coverage_weights_list, np.float64,
            ),
            tx_starts=_to_np(self.tx_starts_list, np.int32),
            tx_ends=_to_np(self.tx_ends_list, np.int32),
            locus_t_indices=_to_np(self.locus_t_list, np.int32),
            locus_count_cols=_to_np(self.locus_ct_list, np.uint8),
            is_spliced=_to_np(self.is_spliced_list, np.int8).astype(bool),
            gdna_log_liks=_to_np(self.gdna_ll_list, np.float64),
            genomic_footprints=_to_np(
                self.genomic_footprints_list, np.int32,
            ),
            frag_ids=_to_np(self.frag_id_list, np.int64),
            frag_class=_to_np(self.frag_class_list, np.int8),
            splice_type=_to_np(self.splice_type_list, np.uint8),
            n_units=n_units,
            n_candidates=n_candidates,
            nrna_base_index=self.ctx.nrna_base,
        )
