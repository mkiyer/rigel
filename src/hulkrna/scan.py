"""hulkrna.scan — CSR builder / buffer scanner.

Converts a ``FragmentBuffer`` into ``ScoredFragments`` for the EM solver.
The ``FragmentRouter`` class builds per-fragment scoring data, routes
unique fragments to deterministic assignment, and accumulates
ambiguous fragments into CSR-formatted equivalence-class data for
the locus-level EM.
"""

import math

import numpy as np

from .buffer import (
    FragmentBuffer,
    FRAG_UNIQUE,
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
        Accumulates unique counts during scan.
    stats : PipelineStats
        Pipeline statistics accumulator.
    index : TranscriptIndex
        Reference index.
    annotations : AnnotationTable or None
        If provided, deterministic-unique and chimeric fragments are
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
        anti_flag = ctx.anti_flag
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
                is_anti = anti_flag if same else (not anti_flag)
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
        """Compute gDNA log-likelihood for a single unspliced fragment."""
        ctx = self.ctx
        sp = ctx.gdna_splice_penalties.get(
            bf.splice_type, _DEFAULT_SPLICE_PENALTY,
        )
        es = bf.exon_strand
        if es == STRAND_POS or es == STRAND_NEG:
            p_ig = ctx.ig_p if es == STRAND_POS else (1.0 - ctx.ig_p)
        else:
            p_ig = 0.5
        log_s = math.log(max(p_ig, LOG_SAFE_FLOOR))
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

        if fc == FRAG_UNIQUE:
            self.stats.em_routed_unique_units += 1
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

        Deterministic-unique (SPLICED_ANNOT + FRAG_UNIQUE) fragments
        are assigned directly via ``estimator.assign_unique``.  All other
        exonic fragments build EM units.  Chimeric fragments are
        recorded in the annotation table (if active) and skipped.
        """
        import logging
        logger = logging.getLogger(__name__)

        estimator = self.estimator
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
                if fc == FRAG_UNIQUE or fc == FRAG_AMBIG_SAME_STRAND:
                    first_t = int(bf.t_inds[0])
                    t_strand = int(index.t_to_strand_arr[first_t])

                    is_anti = estimator.is_antisense(
                        bf.exon_strand, t_strand,
                        self.strand_models.exonic_spliced,
                    )

                    is_unspliced = (bf.splice_type == SPLICE_UNSPLICED)

                    # Transcript-level unspliced sense/antisense (for gDNA init)
                    # and intronic sense/antisense (for nRNA init).
                    n_cand = len(bf.t_inds)
                    weight = 1.0 / n_cand
                    for k, t_idx in enumerate(bf.t_inds):
                        t_idx_int = int(t_idx)

                        if is_unspliced:
                            if is_anti:
                                estimator.transcript_unspliced_antisense[
                                    t_idx_int
                                ] += weight
                            else:
                                estimator.transcript_unspliced_sense[
                                    t_idx_int
                                ] += weight

                        has_unambig_intron = (
                            bf.unambig_intron_bp is not None
                            and bf.unambig_intron_bp[k] > 0
                        )
                        if has_unambig_intron:
                            if is_anti:
                                estimator.transcript_intronic_antisense[
                                    t_idx_int
                                ] += weight
                            else:
                                estimator.transcript_intronic_sense[
                                    t_idx_int
                                ] += weight

                # --- Deterministic unique: SPLICED_ANNOT + FRAG_UNIQUE ---
                if fc == FRAG_UNIQUE and is_spliced_annot:
                    estimator.assign_unique(bf, index, self.strand_models)
                    stats.deterministic_unique_units += 1

                    if annotations is not None:
                        det_tid = int(next(iter(bf.t_inds)))
                        det_gid = int(t_to_g[det_tid])
                        annotations.add(
                            frag_id=int(chunk.frag_id[i]),
                            best_tid=det_tid,
                            best_gid=det_gid,
                            pool=POOL_CODE_MRNA,
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

        return ScoredFragments(
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
