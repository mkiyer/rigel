"""hulkrna.scan — CSR builder / buffer scanner.

Converts a ``FragmentBuffer`` into ``ScoredFragments`` for the EM solver.
The ``FragmentRouter`` class builds per-fragment scoring data, routes
unambiguous fragments to deterministic assignment, and accumulates
ambiguous fragments into CSR-formatted equivalence-class data for
the locus-level EM.
"""

import math

import numpy as np

from .buffer import (
    FragmentBuffer,
    FRAG_UNAMBIG,
    FRAG_AMBIG_SAME_STRAND,
    FRAG_MULTIMAPPER,
    FRAG_CHIMERIC,
)
from .estimator import AbundanceEstimator, ScoredFragments
from .index import TranscriptIndex
from .scoring import (
    FragmentScorer,
    LOG_SAFE_FLOOR,
    LOG_HALF,
    frag_len_log_lik,
    genomic_to_transcript_pos_bisect,
    compute_fragment_weight,
)
from .types import STRAND_POS, STRAND_NEG
from .splice import SPLICE_UNSPLICED, SPLICE_UNANNOT, SPLICE_ANNOT
from .frag_length_model import _TAIL_DECAY_LP
from .stats import PipelineStats
from .annotate import POOL_CODE_MRNA, POOL_CODE_CHIMERIC

# Default splice penalty for splice types not found in the penalty map.
_DEFAULT_SPLICE_PENALTY = 1.0

# Pre-computed negative infinity (avoids repeated np.inf attribute lookup).
_NEG_INF = float('-inf')


class FragmentRouter:
    """Builds the global CSR (ScoredFragments) from a FragmentBuffer.

    State that was previously captured via closures is now explicit
    instance attributes.  Every method is a candidate for Cython / C
    acceleration.

    Parameters
    ----------
    ctx : FragmentScorer
        Pre-computed scoring parameters.
    estimator : AbundanceEstimator
        Accumulates unambig counts during scan.
    stats : PipelineStats
        Pipeline statistics accumulator.
    index : TranscriptIndex
        Reference index.
    annotations : AnnotationTable or None
        If provided, deterministic-unambig and chimeric fragments are
        annotated immediately during the scan.
    """

    def __init__(
        self,
        ctx: FragmentScorer,
        estimator: AbundanceEstimator,
        stats: PipelineStats,
        index: TranscriptIndex,
        strand_models,
        annotations=None,
    ):
        self.ctx = ctx
        self.estimator = estimator
        self.stats = stats
        self.index = index
        self.strand_models = strand_models
        self.annotations = annotations

        # CSR accumulation arrays (array.array for zero-copy → numpy)
        from array import array as _array
        self.offsets = _array('q', [0])          # int64
        self.t_indices_list = _array('i')         # int32
        self.log_liks_list = _array('f')          # float32
        self.count_cols_list = _array('B')        # uint8
        self.coverage_weights_list = _array('f')  # float32
        self.tx_starts_list = _array('i')         # int32
        self.tx_ends_list = _array('i')           # int32

        # Per-unit metadata arrays
        self.locus_t_list = _array('i')           # int32
        self.locus_ct_list = _array('B')          # uint8
        self.is_spliced_list = _array('b')        # int8 (→ bool at finalize)
        self.gdna_ll_list = _array('f')           # float32
        self.genomic_footprints_list = _array('i')  # int32
        self.frag_id_list = _array('i')           # int32
        self.frag_class_list = _array('b')        # int8
        self.splice_type_list = _array('B')       # uint8

        # Multimapper group state
        self.mm_pending: list = []
        self.mm_fid: int = -1

        # Native scoring context (set by FragmentScorer.from_models)
        self._native_ctx = getattr(ctx, '_native_ctx', None)

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

        Dispatches to C++ NativeFragmentScorer when available.
        """
        nc = self._native_ctx
        if nc is not None:
            return nc.score_wta_mrna(
                t_inds, exon_bp, frag_lengths,
                exon_strand, splice_type, nm, read_length, genomic_start,
            )
        ctx = self.ctx
        n_cand = len(t_inds)
        if n_cand == 0:
            return {}

        rl = read_length if read_length > 0 else 1
        log_nm = nm * ctx.mismatch_log_penalty if nm > 0 else 0.0
        oh_log_pen = ctx.overhang_log_penalty

        # Pre-fetch context arrays for tight loop
        t_strand_arr = ctx.t_strand_arr
        t_length_arr = ctx.t_length_arr
        log_p_sense = ctx.log_p_sense
        log_p_antisense = ctx.log_p_antisense
        r1_antisense = ctx.r1_antisense
        fl_log_prob = ctx.fl_log_prob
        fl_max_size = ctx.fl_max_size
        fl_tail_base = ctx.fl_tail_base
        _tail_decay = _TAIL_DECAY_LP
        t_exon_data = ctx._t_exon_data
        has_strand = (exon_strand == 1 or exon_strand == 2)
        has_genomic = genomic_start >= 0

        # --- Single pass: score + find min overhang ---
        min_oh = 0x7FFFFFFF  # sentinel
        scored: list[tuple[int, int, float, int, float, int, int]] = []

        for k in range(n_cand):
            t_idx = int(t_inds[k])
            ebp = int(exon_bp[k]) if exon_bp is not None else rl
            if ebp <= 0:
                continue

            oh = rl - ebp
            if oh < 0:
                oh = 0

            # Fragment length
            flen = int(frag_lengths[k]) if frag_lengths is not None else -1

            # Fragment-length log-likelihood (scalar LUT)
            if flen <= 0 or fl_log_prob is None:
                log_fl = 0.0
            elif flen <= fl_max_size:
                log_fl = float(fl_log_prob[flen])
            else:
                log_fl = fl_tail_base + (flen - fl_max_size) * _tail_decay

            # Strand log-likelihood
            if has_strand:
                t_strand = int(t_strand_arr[t_idx])
                same = (exon_strand == t_strand)
                log_strand = log_p_sense if same else log_p_antisense
                is_anti = r1_antisense if same else (not r1_antisense)
            else:
                log_strand = LOG_HALF
                is_anti = False

            log_lik = log_strand + log_fl + oh * oh_log_pen + log_nm
            ct = splice_type * 2 + int(is_anti)

            # Coverage weight + transcript position
            t_len = int(t_length_arr[t_idx])
            tx_s = 0
            tx_e = flen if flen > 0 else t_len
            cov_wt = 1.0

            if has_genomic and flen > 0:
                exon_data = t_exon_data.get(t_idx) if t_exon_data else None
                if exon_data is not None and t_len > 0:
                    starts, ends, cumsum_before = exon_data
                    t_strand_v = int(t_strand_arr[t_idx])
                    tx_s = genomic_to_transcript_pos_bisect(
                        genomic_start, starts, ends,
                        cumsum_before, t_strand_v, t_len,
                    )
                    tx_e = tx_s + flen
                    # Clipped end for coverage weight
                    cov_end = tx_e if tx_e < t_len else t_len
                    cov_wt = compute_fragment_weight(tx_s, cov_end, t_len)
                else:
                    tx_s = 0
                    tx_e = flen

            scored.append((t_idx, oh, log_lik, ct, cov_wt, tx_s, tx_e))
            if oh < min_oh:
                min_oh = oh

        if not scored:
            return {}

        # --- WTA gating: keep only winners with minimum overhang ---
        result: dict[int, tuple[int, float, int, float, int, int]] = {}
        for t_idx, oh, log_lik, ct, cov_wt, tx_s, tx_e in scored:
            if oh == min_oh:
                result[t_idx] = (oh, log_lik, ct, cov_wt, tx_s, tx_e)
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

        Dispatches to C++ NativeFragmentScorer when available.
        """
        nc = self._native_ctx
        if nc is not None:
            return nc.score_wta_nrna(
                t_inds, exon_bp, intron_bp,
                exon_strand, nm, read_length,
                genomic_footprint, genomic_start,
            )
        ctx = self.ctx
        n_cand = len(t_inds)
        if n_cand == 0:
            return {}

        rl = read_length if read_length > 0 else 1
        log_fl = frag_len_log_lik(ctx, genomic_footprint)
        log_nm = nm * ctx.mismatch_log_penalty if nm > 0 else 0.0
        oh_log_pen = ctx.overhang_log_penalty

        # Pre-fetch context arrays
        t_strand_arr = ctx.t_strand_arr
        t_span_arr = ctx.t_span_arr
        t_length_arr = ctx.t_length_arr
        t_start_arr = ctx.t_start_arr
        log_p_sense = ctx.log_p_sense
        log_p_antisense = ctx.log_p_antisense
        has_strand = (exon_strand == 1 or exon_strand == 2)
        has_genomic = genomic_start >= 0 and genomic_footprint > 0
        nrna_gfp = genomic_footprint if genomic_footprint > 0 else 0

        min_oh = 0x7FFFFFFF
        scored: list[tuple[int, int, float, float, int, int]] = []
        # scored[j] = (t_idx, oh, nrna_ll, cov_wt, tx_s, tx_e)

        for k in range(n_cand):
            t_idx = int(t_inds[k])

            # Single-exon transcript filter
            t_span = int(t_span_arr[t_idx])
            t_exonic = int(t_length_arr[t_idx])
            if t_span <= t_exonic:
                continue  # single-exon → skip

            # Genic overlap
            ebp = int(exon_bp[k]) if exon_bp is not None else rl
            ibp = int(intron_bp[k]) if intron_bp is not None else 0
            span_bp = ebp + ibp
            if span_bp <= 0:
                continue

            oh = rl - span_bp
            if oh < 0:
                oh = 0

            # Strand log-likelihood
            if has_strand:
                t_strand = int(t_strand_arr[t_idx])
                log_strand = log_p_sense if (exon_strand == t_strand) else log_p_antisense
            else:
                log_strand = LOG_HALF

            nrna_ll = log_strand + log_fl + oh * oh_log_pen + log_nm

            # Coverage weight (nRNA uses genomic span)
            tx_s = 0
            tx_e = nrna_gfp if nrna_gfp > 0 else t_span
            cov_wt = 1.0

            if has_genomic:
                tx_start_g = int(t_start_arr[t_idx])
                frag_pos = genomic_start - tx_start_g
                frag_end_pos = frag_pos + genomic_footprint

                # Clamped for fragment weight (ternary instead of min/max)
                frag_pos_c = frag_pos if frag_pos < t_span else t_span
                if frag_pos_c < 0:
                    frag_pos_c = 0
                frag_end_c = frag_end_pos if frag_end_pos < t_span else t_span
                if frag_end_c < 0:
                    frag_end_c = 0
                if frag_end_c > frag_pos_c and t_span > 0:
                    cov_wt = compute_fragment_weight(
                        frag_pos_c, frag_end_c, t_span,
                    )

                # Unclipped for bias model (ternary instead of max)
                tx_s = frag_pos if frag_pos > 0 else 0
                tx_e = tx_s + genomic_footprint

            scored.append((t_idx, oh, nrna_ll, cov_wt, tx_s, tx_e))
            if oh < min_oh:
                min_oh = oh

        if not scored:
            return {}

        # --- WTA gating ---
        result: dict[int, tuple[int, float, float, int, int]] = {}
        for t_idx, oh, nrna_ll, cov_wt, tx_s, tx_e in scored:
            if oh == min_oh:
                result[t_idx] = (oh, nrna_ll, cov_wt, tx_s, tx_e)
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
        if not winners:
            return _NEG_INF, -1, 0
        # Direct extend — avoids intermediate list allocation
        _ti = self.t_indices_list
        _ll_l = self.log_liks_list
        _ct_l = self.count_cols_list
        _cw_l = self.coverage_weights_list
        _ts_l = self.tx_starts_list
        _te_l = self.tx_ends_list
        best_ll = _NEG_INF
        best_t = -1
        best_ct = 0
        for t_idx, (oh, log_lik, ct, cov_wt, tx_s, tx_e) in winners.items():
            _ti.append(t_idx)
            _ll_l.append(log_lik)
            _ct_l.append(ct)
            _cw_l.append(cov_wt)
            _ts_l.append(tx_s)
            _te_l.append(tx_e)
            if log_lik > best_ll:
                best_ll = log_lik
                best_t = t_idx
                best_ct = ct
        return best_ll, best_t, best_ct

    def _emit_nrna(self, winners: dict[int, tuple[int, float, float, int, int]]) -> None:
        """Append nRNA winner candidates to CSR lists."""
        if not winners:
            return
        nrna_base = self.ctx.nrna_base
        _ti = self.t_indices_list
        _ll_l = self.log_liks_list
        _ct_l = self.count_cols_list
        _cw_l = self.coverage_weights_list
        _ts_l = self.tx_starts_list
        _te_l = self.tx_ends_list
        for t_idx, (oh, nrna_ll, cov_wt, tx_s, tx_e) in winners.items():
            _ti.append(nrna_base + t_idx)
            _ll_l.append(nrna_ll)
            _ct_l.append(0)
            _cw_l.append(cov_wt)
            _ts_l.append(tx_s)
            _te_l.append(tx_e)

    # ------------------------------------------------------------------
    # Per-unit metadata
    # ------------------------------------------------------------------

    def _gdna_log_lik(self, bf) -> float:
        """Compute gDNA log-likelihood for a single unspliced fragment.

        gDNA has no strand bias: strand probability is always 0.5.
        """
        ctx = self.ctx
        sp = ctx.gdna_splice_penalties.get(
            bf.splice_type, _DEFAULT_SPLICE_PENALTY,
        )
        log_s = LOG_HALF
        log_i = frag_len_log_lik(ctx, bf.genomic_footprint)
        return log_s + log_i + math.log(max(sp, LOG_SAFE_FLOOR))

    def _finalize_unit_metadata(self, bf) -> None:
        """Record is_spliced and gDNA log-likelihood for this EM unit."""
        st = bf.splice_type
        is_spl = (st == SPLICE_ANNOT or st == SPLICE_UNANNOT)
        self.is_spliced_list.append(is_spl)
        self.genomic_footprints_list.append(int(bf.genomic_footprint))

        if not is_spl:
            self.gdna_ll_list.append(self._gdna_log_lik(bf))
        else:
            self.gdna_ll_list.append(_NEG_INF)

    # ------------------------------------------------------------------
    # Single-fragment EM unit
    # ------------------------------------------------------------------

    def _add_single_fragment(self, bf, chunk, i, fc) -> None:
        """Process a non-chimeric, non-multimapper EM unit."""
        nc = self._native_ctx
        if nc is not None:
            # Fused score+emit: single C++ call, bytes output
            r = nc.score_emit_fragment(
                bf.t_inds, bf.exon_bp, bf.frag_lengths, bf.intron_bp,
                bf.exon_strand, bf.splice_type, bf.nm, bf.read_length,
                bf.genomic_footprint, bf.genomic_start,
            )
            best_ll, best_t, best_ct = r[0], r[1], r[2]
            self.t_indices_list.frombytes(r[3])
            self.log_liks_list.frombytes(r[4])
            self.count_cols_list.frombytes(r[5])
            self.coverage_weights_list.frombytes(r[6])
            self.tx_starts_list.frombytes(r[7])
            self.tx_ends_list.frombytes(r[8])
        else:
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

        if fc == FRAG_UNAMBIG:
            self.stats.em_routed_unambig_units += 1
        elif fc == FRAG_AMBIG_SAME_STRAND:
            self.stats.em_routed_ambig_same_strand_units += 1
        else:
            self.stats.em_routed_ambig_opp_strand_units += 1

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
                best_gdna_ll = _NEG_INF
                best_footprint = 0
                for bf in mm_pending:
                    st = bf.splice_type
                    if st == SPLICE_ANNOT or st == SPLICE_UNANNOT:
                        continue
                    gdna_ll = self._gdna_log_lik(bf)
                    if gdna_ll > best_gdna_ll:
                        best_gdna_ll = gdna_ll
                        best_footprint = int(bf.genomic_footprint)
                self.gdna_ll_list.append(best_gdna_ll)
                self.genomic_footprints_list.append(best_footprint)
            else:
                self.gdna_ll_list.append(_NEG_INF)
                # Use first hit's footprint as representative
                self.genomic_footprints_list.append(
                    int(mm_pending[0].genomic_footprint)
                )
        else:
            self.stats.n_gated_out += 1

    # ------------------------------------------------------------------
    # Main scan driver
    # ------------------------------------------------------------------

    def scan(self, buffer: FragmentBuffer, log_every: int) -> ScoredFragments:
        """Single pass over buffer.  Returns packed ScoredFragments.

        Deterministic-unique (SPLICED_ANNOT + FRAG_UNAMBIG) fragments
        are assigned directly via ``estimator.assign_unambig``.  All other
        exonic fragments build EM units.  Chimeric fragments are
        recorded in the annotation table (if active) and skipped.

        Uses the C++ chunk-level fast path when a NativeFragmentScorer is
        available and the buffer yields proper ``_FinalizedChunk`` objects
        (with columnar NumPy arrays).  Falls back to the per-fragment
        Python path otherwise.
        """
        nc = self._native_ctx
        if nc is not None:
            # Peek at the first chunk to see if it has columnar arrays.
            # Test mocks may use a simple _Chunk that only supports __getitem__.
            chunks = list(buffer.iter_chunks())
            if chunks and hasattr(chunks[0], 't_offsets'):
                return self._scan_native(chunks, buffer, log_every)
        return self._scan_python(buffer, log_every)

    def _scan_native(
        self, chunks, buffer: FragmentBuffer, log_every: int,
    ) -> ScoredFragments:
        """Chunk-level C++ scan path (fast).

        Uses the fused two-pass C++ scorer: pass 1 counts EM units
        and candidates (accumulating unambig counts), pass 2 fills
        pre-allocated numpy arrays at exact sizes.  Multimapper
        fragments remain on the Python path and are concatenated
        after the C++ pass.
        """
        import logging
        logger = logging.getLogger(__name__)

        nc = self._native_ctx
        estimator = self.estimator
        stats = self.stats
        index = self.index
        t_to_g = self.ctx.t_to_g
        annotations = self.annotations

        from .scoring import LOG_SAFE_FLOOR
        _gdna_sp = self.ctx.gdna_splice_penalties.get(
            SPLICE_UNSPLICED, _DEFAULT_SPLICE_PENALTY,
        )
        gdna_log_sp = math.log(max(_gdna_sp, LOG_SAFE_FLOOR))

        t_strand_arr = np.ascontiguousarray(
            index.t_to_strand_arr, dtype=np.int8)

        # ---- Prepare chunk arrays as tuples for C++ ----
        chunk_arrays = []
        for chunk in chunks:
            chunk_arrays.append((
                np.ascontiguousarray(
                    chunk.t_offsets, dtype=np.int64),
                np.ascontiguousarray(
                    chunk.t_indices, dtype=np.int32),
                np.ascontiguousarray(
                    chunk.frag_lengths, dtype=np.int32),
                np.ascontiguousarray(
                    chunk.exon_bp, dtype=np.int32),
                np.ascontiguousarray(
                    chunk.intron_bp, dtype=np.int32),
                np.ascontiguousarray(
                    chunk.unambig_intron_bp, dtype=np.int32),
                np.ascontiguousarray(
                    chunk.splice_type, dtype=np.uint8),
                np.ascontiguousarray(
                    chunk.exon_strand, dtype=np.uint8),
                np.ascontiguousarray(
                    chunk.fragment_classes, dtype=np.uint8),
                np.ascontiguousarray(
                    chunk.frag_id, dtype=np.int64),
                np.ascontiguousarray(
                    chunk.read_length, dtype=np.uint32),
                np.ascontiguousarray(
                    chunk.genomic_footprint, dtype=np.int32),
                np.ascontiguousarray(
                    chunk.genomic_start, dtype=np.int32),
                np.ascontiguousarray(
                    chunk.nm, dtype=np.uint16),
            ))

        # ---- Fused two-pass C++ scoring ----
        result = nc.fused_score_buffer(
            chunk_arrays,
            t_strand_arr,
            estimator.transcript_exonic_sense,
            estimator.transcript_exonic_antisense,
            estimator.transcript_unspliced_sense,
            estimator.transcript_unspliced_antisense,
            estimator.transcript_intronic_sense,
            estimator.transcript_intronic_antisense,
            estimator.unambig_counts,
            gdna_log_sp,
        )

        (
            cpp_offsets, cpp_ti, cpp_ll, cpp_ct, cpp_cw,
            cpp_ts, cpp_te,
            cpp_locus_t, cpp_locus_ct, cpp_is_spliced,
            cpp_gdna_ll, cpp_gfp, cpp_fid, cpp_fclass,
            cpp_stype,
            det_tids, det_fids,
            n_det, n_em_u, n_em_as, n_em_ao, n_gated, n_chim,
        ) = result

        stats.deterministic_unambig_units += int(n_det)
        stats.em_routed_unambig_units += int(n_em_u)
        stats.em_routed_ambig_same_strand_units += int(n_em_as)
        stats.em_routed_ambig_opp_strand_units += int(n_em_ao)
        stats.n_gated_out += int(n_gated)

        # ---- Annotations (lightweight per-chunk pass) ----
        n_processed = 0
        for ci, chunk in enumerate(chunks):
            frag_classes = chunk_arrays[ci][8]  # pre-computed

            if annotations is not None:
                chimeric_indices = np.where(
                    frag_classes == FRAG_CHIMERIC)[0]
                for idx in chimeric_indices:
                    annotations.add(
                        frag_id=int(chunk.frag_id[idx]),
                        best_tid=-1,
                        best_gid=-1,
                        pool=POOL_CODE_CHIMERIC,
                        posterior=0.0,
                        frag_class=FRAG_CHIMERIC,
                        n_candidates=0,
                        splice_type=int(chunk.splice_type[idx]),
                    )

            # Multimapper handling (Python path)
            mm_indices = np.where(
                frag_classes == FRAG_MULTIMAPPER)[0]
            for idx in mm_indices:
                bf = chunk[int(idx)]
                fid = int(chunk.frag_id[idx])
                if fid != self.mm_fid:
                    if self.mm_pending:
                        self._flush_mm_group()
                        stats.em_routed_multimapper_units += 1
                    self.mm_fid = fid
                    self.mm_pending.clear()
                    self.mm_pending.append(bf)
                else:
                    self.mm_pending.append(bf)

            n_processed += chunk.size
            if n_processed % log_every < chunk.size:
                logger.debug(
                    f"  Scan: {n_processed:,} / "
                    f"{buffer.total_fragments:,}"
                )

        # Det-unambig annotations
        if annotations is not None and len(det_tids) > 0:
            for j in range(len(det_tids)):
                tid = int(det_tids[j])
                gid = int(t_to_g[tid])
                annotations.add(
                    frag_id=int(det_fids[j]),
                    best_tid=tid,
                    best_gid=gid,
                    pool=POOL_CODE_MRNA,
                    posterior=1.0,
                    frag_class=FRAG_UNAMBIG,
                    n_candidates=1,
                    splice_type=int(SPLICE_ANNOT),
                )

        # Flush last multimapper group
        if self.mm_pending:
            self._flush_mm_group()
            stats.em_routed_multimapper_units += 1

        # Free chunk_arrays — no longer needed
        del chunk_arrays

        # ---- Merge C++ and MM results ----
        mm_n_units = len(self.offsets) - 1
        mm_n_cands = len(self.t_indices_list)

        def _to_np(arr, dtype):
            if len(arr) == 0:
                return np.empty(0, dtype=dtype)
            return np.frombuffer(arr, dtype=dtype).copy()

        if mm_n_units > 0:
            # Merge: concatenate C++ arrays with MM arrays
            mm_off = _to_np(self.offsets, np.int64)
            offsets = np.concatenate([
                cpp_offsets,
                mm_off[1:] + cpp_offsets[-1],
            ])
            t_indices = np.concatenate([
                cpp_ti,
                _to_np(self.t_indices_list, np.int32),
            ])
            log_liks = np.concatenate([
                cpp_ll,
                _to_np(self.log_liks_list, np.float32),
            ])
            count_cols = np.concatenate([
                cpp_ct,
                _to_np(self.count_cols_list, np.uint8),
            ])
            coverage_weights = np.concatenate([
                cpp_cw,
                _to_np(self.coverage_weights_list, np.float32),
            ])
            tx_starts = np.concatenate([
                cpp_ts,
                _to_np(self.tx_starts_list, np.int32),
            ])
            tx_ends = np.concatenate([
                cpp_te,
                _to_np(self.tx_ends_list, np.int32),
            ])
            locus_t_indices = np.concatenate([
                cpp_locus_t,
                _to_np(self.locus_t_list, np.int32),
            ])
            locus_count_cols = np.concatenate([
                cpp_locus_ct,
                _to_np(self.locus_ct_list, np.uint8),
            ])
            is_spliced = np.concatenate([
                cpp_is_spliced,
                _to_np(self.is_spliced_list, np.int8),
            ]).astype(bool)
            gdna_log_liks = np.concatenate([
                cpp_gdna_ll,
                _to_np(self.gdna_ll_list, np.float32),
            ])
            genomic_footprints = np.concatenate([
                cpp_gfp,
                _to_np(self.genomic_footprints_list, np.int32),
            ])
            frag_ids = np.concatenate([
                cpp_fid,
                _to_np(self.frag_id_list, np.int32),
            ])
            frag_class = np.concatenate([
                cpp_fclass,
                _to_np(self.frag_class_list, np.int8),
            ])
            splice_type = np.concatenate([
                cpp_stype,
                _to_np(self.splice_type_list, np.uint8),
            ])
            n_units = int(len(cpp_offsets) - 1) + mm_n_units
            n_candidates = int(len(cpp_ti)) + mm_n_cands
        else:
            # No multimappers — use C++ arrays directly
            offsets = cpp_offsets
            t_indices = cpp_ti
            log_liks = cpp_ll
            count_cols = cpp_ct
            coverage_weights = cpp_cw
            tx_starts = cpp_ts
            tx_ends = cpp_te
            locus_t_indices = cpp_locus_t
            locus_count_cols = cpp_locus_ct
            is_spliced = np.asarray(
                cpp_is_spliced, dtype=np.int8).astype(bool)
            gdna_log_liks = cpp_gdna_ll
            genomic_footprints = cpp_gfp
            frag_ids = cpp_fid
            frag_class = cpp_fclass
            splice_type = cpp_stype
            n_units = int(len(cpp_offsets) - 1)
            n_candidates = int(len(cpp_ti))

        return ScoredFragments(
            offsets=offsets,
            t_indices=t_indices,
            log_liks=log_liks,
            count_cols=count_cols,
            coverage_weights=coverage_weights,
            tx_starts=tx_starts,
            tx_ends=tx_ends,
            locus_t_indices=locus_t_indices,
            locus_count_cols=locus_count_cols,
            is_spliced=is_spliced,
            gdna_log_liks=gdna_log_liks,
            genomic_footprints=genomic_footprints,
            frag_ids=frag_ids,
            frag_class=frag_class,
            splice_type=splice_type,
            n_units=n_units,
            n_candidates=n_candidates,
            nrna_base_index=self.ctx.nrna_base,
        )

    def _scan_python(
        self, buffer: FragmentBuffer, log_every: int,
    ) -> ScoredFragments:
        """Per-fragment Python scan path (fallback when native unavailable).

        Implementation lives in :mod:`hulkrna.pyfallback` to keep the
        main module lean.  Will be removed once C++ is the sole path.
        """
        from .pyfallback import scan_python as _scan_python_fb
        return _scan_python_fb(self, buffer, log_every)
