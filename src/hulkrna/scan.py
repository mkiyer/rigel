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
    frag_len_log_lik,
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

        # CSR accumulation lists
        self.offsets: list[int] = [0]
        self.t_indices_list: list[int] = []
        self.log_liks_list: list[float] = []
        self.count_cols_list: list[int] = []

        # Per-unit metadata
        self.locus_t_list: list[int] = []
        self.locus_ct_list: list[int] = []
        self.is_spliced_list: list[bool] = []
        self.gdna_ll_list: list[float] = []
        self.frag_id_list: list[int] = []
        self.frag_class_list: list[int] = []
        self.splice_type_list: list[int] = []

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
    ) -> dict[int, tuple[int, float, int]]:
        """Score mRNA candidates with winner-take-all gating.

        Returns dict mapping t_idx → (overhang, log_lik, count_col).
        Only WTA winners (minimum overhang) are included.
        """
        ctx = self.ctx
        fl = read_length if read_length > 0 else 1
        log_nm = nm * ctx.mismatch_log_penalty if nm > 0 else 0.0

        # Strand log-probabilities
        if exon_strand == 1 or exon_strand == 2:
            log_same = ctx.log_p_sense
            log_diff = ctx.log_p_antisense
        else:
            log_same = log_diff = LOG_HALF

        # Pass 1 — find minimum overhang and collect per-transcript best
        best_per_t: dict[int, tuple[int, float, int]] = {}
        for k in range(len(t_inds)):
            t_idx_int = int(t_inds[k])
            ebp = int(exon_bp[k]) if exon_bp is not None else fl
            if ebp <= 0:
                continue
            oh = max(fl - ebp, 0)

            # Per-candidate fragment-length log-probability
            log_fl = frag_len_log_lik(
                ctx, int(frag_lengths[k]) if frag_lengths is not None else -1
            )

            # Inlined scoring
            t_strand = int(ctx.t_strand_arr[t_idx_int])
            if exon_strand == t_strand:
                log_strand = log_same
                anti = ctx.anti_flag
            else:
                log_strand = log_diff
                anti = not ctx.anti_flag

            log_lik = log_strand + log_fl + oh * ctx.overhang_log_penalty + log_nm
            ct = splice_type * 2 + int(anti)

            prev = best_per_t.get(t_idx_int)
            if prev is None or oh < prev[0] or (
                oh == prev[0] and log_lik > prev[1]
            ):
                best_per_t[t_idx_int] = (oh, log_lik, ct)

        if not best_per_t:
            return {}

        # Pass 2 — keep only WTA winners
        min_oh = min(v[0] for v in best_per_t.values())
        return {
            t: v for t, v in best_per_t.items() if v[0] == min_oh
        }

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
    ) -> dict[int, tuple[int, float]]:
        """Score nRNA candidates with winner-take-all gating.

        Returns dict mapping t_idx → (overhang, nrna_log_lik).
        Only WTA winners (minimum overhang) are included.
        """
        ctx = self.ctx
        fl = read_length if read_length > 0 else 1
        log_fl = frag_len_log_lik(ctx, genomic_footprint)
        log_nm = nm * ctx.mismatch_log_penalty if nm > 0 else 0.0

        if exon_strand == 1 or exon_strand == 2:
            log_same = ctx.log_p_sense
            log_diff = ctx.log_p_antisense
        else:
            log_same = log_diff = LOG_HALF

        best_per_t: dict[int, tuple[int, float]] = {}
        for k in range(len(t_inds)):
            t_idx_int = int(t_inds[k])
            ebp = int(exon_bp[k]) if exon_bp is not None else fl
            ibp = int(intron_bp[k]) if intron_bp is not None else 0
            span_bp = ebp + ibp
            if span_bp <= 0:
                continue
            oh = max(fl - span_bp, 0)

            g_idx = int(ctx.t_to_g[t_idx_int])
            gene_strand = int(ctx.g_strand_arr[g_idx])
            if exon_strand == gene_strand:
                log_strand = log_same
            else:
                log_strand = log_diff

            nrna_ll = log_strand + log_fl + oh * ctx.overhang_log_penalty + log_nm

            prev = best_per_t.get(t_idx_int)
            if prev is None or oh < prev[0] or (
                oh == prev[0] and nrna_ll > prev[1]
            ):
                best_per_t[t_idx_int] = (oh, nrna_ll)

        if not best_per_t:
            return {}

        min_oh = min(v[0] for v in best_per_t.values())
        return {
            t: v for t, v in best_per_t.items() if v[0] == min_oh
        }

    # ------------------------------------------------------------------
    # Emit candidates to CSR
    # ------------------------------------------------------------------

    def _emit_mrna(
        self, winners: dict[int, tuple[int, float, int]],
    ) -> tuple[float, int, int]:
        """Append mRNA winner candidates to CSR lists.

        Returns (best_log_lik, best_t_idx, best_count_col).
        """
        best_ll = -np.inf
        best_t = -1
        best_ct = 0
        for t_idx, (oh, log_lik, ct) in winners.items():
            self.t_indices_list.append(t_idx)
            self.log_liks_list.append(log_lik)
            self.count_cols_list.append(ct)
            if log_lik > best_ll:
                best_ll = log_lik
                best_t = t_idx
                best_ct = ct
        return best_ll, best_t, best_ct

    def _emit_nrna(self, winners: dict[int, tuple[int, float]]) -> None:
        """Append nRNA winner candidates to CSR lists."""
        nrna_base = self.ctx.nrna_base
        for t_idx, (oh, nrna_ll) in winners.items():
            self.t_indices_list.append(nrna_base + t_idx)
            self.log_liks_list.append(nrna_ll)
            self.count_cols_list.append(0)

    # ------------------------------------------------------------------
    # Per-unit metadata
    # ------------------------------------------------------------------

    def _finalize_unit_metadata(self, bf) -> None:
        """Record is_spliced and gDNA log-likelihood for this EM unit."""
        ctx = self.ctx
        st = bf.splice_type
        is_spl = (st == SPLICE_ANNOT or st == SPLICE_UNANNOT)
        self.is_spliced_list.append(is_spl)

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
        )
        best_ll, best_t, best_ct = self._emit_mrna(mrna_winners)

        # --- nRNA candidates (independent pool, skip only SPLICED_ANNOT) ---
        if not is_spliced_annot:
            nrna_winners = self._score_wta_nrna(
                bf.t_inds, bf.exon_bp, bf.intron_bp,
                bf.exon_strand, bf.nm, bf.read_length,
                bf.genomic_footprint,
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
        merged_mrna: dict[int, tuple[int, float, int]] = {}
        for bf in mm_pending:
            winners = self._score_wta_mrna(
                bf.t_inds, bf.exon_bp, bf.frag_lengths,
                bf.exon_strand, bf.splice_type, bf.nm, bf.read_length,
            )
            # Merge: keep best per transcript (lower oh wins, then higher ll)
            for t_idx, (oh, ll, ct) in winners.items():
                prev = merged_mrna.get(t_idx)
                if prev is None or oh < prev[0] or (
                    oh == prev[0] and ll > prev[1]
                ):
                    merged_mrna[t_idx] = (oh, ll, ct)

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
            merged_nrna: dict[int, tuple[int, float]] = {}
            for bf in mm_pending:
                if bf.splice_type == SPLICE_ANNOT:
                    continue
                winners = self._score_wta_nrna(
                    bf.t_inds, bf.exon_bp, bf.intron_bp,
                    bf.exon_strand, bf.nm, bf.read_length,
                    bf.genomic_footprint,
                )
                for t_idx, (oh, nrna_ll) in winners.items():
                    prev = merged_nrna.get(t_idx)
                    if prev is None or oh < prev[0] or (
                        oh == prev[0] and nrna_ll > prev[1]
                    ):
                        merged_nrna[t_idx] = (oh, nrna_ll)

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
                self.gdna_ll_list.append(best_gdna_ll)
            else:
                self.gdna_ll_list.append(-np.inf)
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

        return ScanData(
            offsets=np.array(self.offsets, dtype=np.int64),
            t_indices=np.array(self.t_indices_list, dtype=np.int32),
            log_liks=np.array(self.log_liks_list, dtype=np.float64),
            count_cols=np.array(self.count_cols_list, dtype=np.uint8),
            locus_t_indices=np.array(self.locus_t_list, dtype=np.int32),
            locus_count_cols=np.array(self.locus_ct_list, dtype=np.uint8),
            is_spliced=np.array(self.is_spliced_list, dtype=bool),
            gdna_log_liks=np.array(self.gdna_ll_list, dtype=np.float64),
            frag_ids=np.array(self.frag_id_list, dtype=np.int64),
            frag_class=np.array(self.frag_class_list, dtype=np.int8),
            splice_type=np.array(self.splice_type_list, dtype=np.uint8),
            n_units=n_units,
            n_candidates=n_candidates,
            nrna_base_index=self.ctx.nrna_base,
        )
