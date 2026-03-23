/**
 * scoring.cpp — nanobind C++ extension for hot-path scoring functions.
 *
 * This module provides performance-critical kernels called from the
 * Python scoring loop.  Functions are exposed as ``rigel._scoring_impl``.
 *
 * Contents:
 *   - compute_fragment_weight(frag_start, frag_end, transcript_length)
 *   - NativeFragmentScorer class:
 *       .fused_score_buffer(...)     — two-pass bulk scoring (production hot path)
 *
 * Build:
 *   Part of the rigel scikit-build-core build — see CMakeLists.txt.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "constants.h"

namespace nb = nanobind;

// Use constants from constants.h via the rigel namespace.
using rigel::LOG_HALF;
using rigel::TAIL_DECAY_LP;
using rigel::STRAND_NEG;
using rigel::SPLICE_UNSPLICED;
using rigel::SPLICE_SPLICED_UNANNOT;
using rigel::SPLICE_SPLICED_ANNOT;
using rigel::FRAG_UNAMBIG;
using rigel::FRAG_AMBIG_SAME_STRAND;
using rigel::FRAG_AMBIG_OPP_STRAND;
using rigel::FRAG_MULTIMAPPER;
using rigel::FRAG_CHIMERIC;
using rigel::POOL_CODE_MRNA;
using rigel::POOL_CODE_NRNA;
using rigel::POOL_CODE_GDNA;
using rigel::POOL_CODE_INTERGENIC;
using rigel::POOL_CODE_CHIMERIC;

// Maximum candidates handled via stack-allocated buffer before falling
// back to heap allocation.  Avoids malloc in the common case where a
// fragment overlaps a moderate number of transcripts.
static constexpr int SCORED_STACK_CAPACITY = 64;

// ----------------------------------------------------------------
// Array type aliases
// ----------------------------------------------------------------

using i32_1d  = nb::ndarray<const int32_t,  nb::ndim<1>, nb::c_contig>;
using i8_1d   = nb::ndarray<const int8_t,  nb::ndim<1>, nb::c_contig>;
using f64_1d  = nb::ndarray<const double,  nb::ndim<1>, nb::c_contig>;
using i64_1d  = nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig>;
using u8_1d   = nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig>;
using u16_1d  = nb::ndarray<const uint16_t, nb::ndim<1>, nb::c_contig>;
using u32_1d  = nb::ndarray<const uint32_t, nb::ndim<1>, nb::c_contig>;
using f64_mut = nb::ndarray<double, nb::ndim<1>, nb::c_contig>;
using f64_2d_mut = nb::ndarray<double, nb::ndim<2>, nb::c_contig>;

// ----------------------------------------------------------------
// compute_fragment_weight
// ----------------------------------------------------------------
//
// Inverse coverage-capacity weight for a fragment on a transcript.
//
// Uses the trapezoid coverage model: under uniform random sampling
// the expected coverage at transcript-relative position x on a
// transcript of spliced length L for a fragment of length f is
//   min(x, w, L - x)   where  w = min(f, L/2).
//
// A fragment in the plateau region gets weight 1.0.  A fragment near
// a transcript edge gets weight > 1.0, reflecting that it is less
// likely to be generated yet was still observed — stronger evidence
// for that transcript.

static double compute_fragment_weight(int32_t frag_start,
                                      int32_t frag_end,
                                      int32_t transcript_length) {
    int32_t f_len = frag_end - frag_start;
    if (f_len <= 0 || transcript_length <= 0)
        return 1.0;

    double L = static_cast<double>(transcript_length);
    double w = std::min(static_cast<double>(f_len), L / 2.0);
    double Lw = L - w;

    double fs = static_cast<double>(frag_start);
    double fe = static_cast<double>(frag_end);

    double total_area = 0.0;

    // Left ramp: c(x) = x  for x in [0, w)
    {
        double s = std::clamp(fs, 0.0, w);
        double e = std::clamp(fe, 0.0, w);
        if (s < e)
            total_area += (e - s) * (s + e) * 0.5;
    }

    // Plateau: c(x) = w  for x in [w, L - w)
    {
        double s = std::clamp(fs, w, Lw);
        double e = std::clamp(fe, w, Lw);
        if (s < e)
            total_area += (e - s) * w;
    }

    // Right ramp: c(x) = L - x  for x in [L - w, L]
    {
        double s = std::clamp(fs, Lw, L);
        double e = std::clamp(fe, Lw, L);
        if (s < e)
            total_area += (e - s) * ((L - s) + (L - e)) * 0.5;
    }

    double area_per_base = total_area / static_cast<double>(f_len);
    if (area_per_base < 1.0)
        area_per_base = 1.0;

    return w / area_per_base;
}

// ================================================================
// NativeFragmentScorer
// ================================================================
//
// Holds all pre-computed scoring parameters and index arrays so that
// scoring and candidate pruning execute entirely in C++ with
// zero Python round-trips.  Constructed once per pipeline run from
// the Python FragmentScorer.

class NativeFragmentScorer {

    // --- Scalar scoring parameters ---
    double log_p_sense_;
    double log_p_antisense_;
    bool   r1_antisense_;
    double oh_log_pen_;
    double mm_log_pen_;
    int32_t fl_max_size_;
    double  fl_tail_base_;
    bool    has_fl_lut_;
    int32_t n_transcripts_;

    // --- Copied index arrays ---
    std::vector<int8_t>  t_strand_;   // int8[n_transcripts]
    std::vector<int32_t> t_length_;   // spliced exonic length
    std::vector<int32_t> t_span_;     // genomic span (incl introns)
    std::vector<int32_t> t_start_;    // genomic start coordinate
    std::vector<uint8_t> t_is_nrna_;  // uint8[n_transcripts] — 1=nRNA, 0=mRNA

    // --- Pool-separated likelihood pruning ---
    double max_ll_delta_;             // Δ = -log(ε), pre-computed

    // --- Fragment-length LUT — RNA model (copied from numpy) ---
    std::vector<double> fl_log_prob_;

    // --- Fragment-length LUT — gDNA model (copied from numpy) ---
    std::vector<double> gdna_fl_log_prob_;
    int32_t gdna_fl_max_size_;
    double  gdna_fl_tail_base_;
    bool    has_gdna_fl_lut_;

    // --- Exon data in CSR layout for genomic→transcript mapping ---
    std::vector<int32_t> exon_offsets_;  // [n_transcripts + 1]
    std::vector<int32_t> exon_starts_;   // flattened
    std::vector<int32_t> exon_ends_;     // flattened
    std::vector<int32_t> exon_cumsum_;   // flattened cumsum_before

    // --- Internal helpers (inlined in the hot loop) ---

    inline double frag_len_log_lik(int32_t flen) const {
        if (flen <= 0 || !has_fl_lut_) return 0.0;
        if (flen <= fl_max_size_)
            return fl_log_prob_[static_cast<size_t>(flen)];
        return fl_tail_base_ + (flen - fl_max_size_) * TAIL_DECAY_LP;
    }

    inline double gdna_frag_len_log_lik(int32_t flen) const {
        if (flen <= 0 || !has_gdna_fl_lut_) return 0.0;
        if (flen <= gdna_fl_max_size_)
            return gdna_fl_log_prob_[static_cast<size_t>(flen)];
        return gdna_fl_tail_base_ + (flen - gdna_fl_max_size_) * TAIL_DECAY_LP;
    }

    inline int32_t genomic_to_tx_pos(int32_t genomic_pos,
                                     int32_t t_idx) const {
        int32_t begin   = exon_offsets_[t_idx];
        int32_t end     = exon_offsets_[t_idx + 1];
        int32_t n_exons = end - begin;
        if (n_exons <= 0) return 0;

        const int32_t* starts = exon_starts_.data() + begin;
        const int32_t* ends   = exon_ends_.data()   + begin;
        const int32_t* cumsum = exon_cumsum_.data()  + begin;

        // bisect_right(starts, genomic_pos) - 1
        int ei = static_cast<int>(
            std::upper_bound(starts, starts + n_exons, genomic_pos)
            - starts
        ) - 1;

        int32_t offset;
        if (ei < 0) {
            offset = 0;
        } else if (genomic_pos >= ends[ei]) {
            // In intron after exon ei (or past last exon)
            offset = cumsum[ei] + (ends[ei] - starts[ei]);
        } else {
            // Inside exon ei
            offset = cumsum[ei] + (genomic_pos - starts[ei]);
        }

        int32_t t_len = t_length_[t_idx];
        if (offset < 0) offset = 0;
        else if (offset > t_len) offset = t_len;

        if (t_strand_[t_idx] == STRAND_NEG)
            offset = t_len - offset;
        return offset;
    }

public:

    // ---------------------------------------------------------------
    // Constructor — copies all data from Python into C++ owned storage
    // ---------------------------------------------------------------

    NativeFragmentScorer(
        double log_p_sense,
        double log_p_antisense,
        bool   r1_antisense,
        double oh_log_pen,
        double mm_log_pen,
        nb::object fl_log_prob_obj,
        int32_t fl_max_size,
        double  fl_tail_base,
        nb::object gdna_fl_log_prob_obj,
        int32_t gdna_fl_max_size,
        double  gdna_fl_tail_base,
        i8_1d   t_strand_arr,
        i32_1d  t_length_arr,
        i32_1d  t_span_arr,
        i32_1d  t_start_arr,
        nb::object t_exon_data_obj,
        u8_1d   t_is_nrna_arr,
        double  pruning_max_ll_delta)
      : log_p_sense_(log_p_sense),
        log_p_antisense_(log_p_antisense),
        r1_antisense_(r1_antisense),
        oh_log_pen_(oh_log_pen),
        mm_log_pen_(mm_log_pen),
        fl_max_size_(fl_max_size),
        fl_tail_base_(fl_tail_base),
        has_fl_lut_(false),
        gdna_fl_max_size_(gdna_fl_max_size),
        gdna_fl_tail_base_(gdna_fl_tail_base),
        has_gdna_fl_lut_(false),
        max_ll_delta_(pruning_max_ll_delta)
    {
        // Copy index arrays
        n_transcripts_ = static_cast<int32_t>(t_strand_arr.shape(0));
        {
            const auto* p = t_strand_arr.data();
            t_strand_.assign(p, p + n_transcripts_);
        }
        {
            const auto* p = t_length_arr.data();
            t_length_.assign(p, p + n_transcripts_);
        }
        {
            const auto* p = t_span_arr.data();
            t_span_.assign(p, p + n_transcripts_);
        }
        {
            const auto* p = t_start_arr.data();
            t_start_.assign(p, p + n_transcripts_);
        }
        {
            const auto* p = t_is_nrna_arr.data();
            t_is_nrna_.assign(p, p + n_transcripts_);
        }

        // Copy fragment-length LUT (RNA model)
        if (!fl_log_prob_obj.is_none()) {
            auto fl_arr = nb::cast<f64_1d>(fl_log_prob_obj);
            const double* p = fl_arr.data();
            int32_t n = static_cast<int32_t>(fl_arr.shape(0));
            fl_log_prob_.assign(p, p + n);
            has_fl_lut_ = true;
        }

        // Copy fragment-length LUT (gDNA model)
        if (!gdna_fl_log_prob_obj.is_none()) {
            auto fl_arr = nb::cast<f64_1d>(gdna_fl_log_prob_obj);
            const double* p = fl_arr.data();
            int32_t n = static_cast<int32_t>(fl_arr.shape(0));
            gdna_fl_log_prob_.assign(p, p + n);
            has_gdna_fl_lut_ = true;
        }

        // Build exon CSR from Python dict
        exon_offsets_.assign(n_transcripts_ + 1, 0);
        if (!t_exon_data_obj.is_none()) {
            nb::dict exon_dict = nb::cast<nb::dict>(t_exon_data_obj);

            // First pass: count exons per transcript
            for (auto item : exon_dict) {
                int32_t t_idx = nb::cast<int32_t>(item.first);
                nb::tuple tup = nb::cast<nb::tuple>(item.second);
                nb::object starts_obj = tup[0];
                exon_offsets_[t_idx + 1] =
                    static_cast<int32_t>(nb::len(starts_obj));
            }

            // Prefix sum
            for (int32_t i = 0; i < n_transcripts_; ++i)
                exon_offsets_[i + 1] += exon_offsets_[i];

            int32_t total = exon_offsets_[n_transcripts_];
            exon_starts_.resize(total);
            exon_ends_.resize(total);
            exon_cumsum_.resize(total);

            // Second pass: fill data
            for (auto item : exon_dict) {
                int32_t t_idx = nb::cast<int32_t>(item.first);
                nb::tuple tup = nb::cast<nb::tuple>(item.second);
                nb::tuple starts_t = nb::cast<nb::tuple>(tup[0]);
                nb::tuple ends_t   = nb::cast<nb::tuple>(tup[1]);
                nb::tuple cum_t    = nb::cast<nb::tuple>(tup[2]);
                int32_t off = exon_offsets_[t_idx];
                int32_t n   = static_cast<int32_t>(nb::len(starts_t));
                for (int32_t j = 0; j < n; ++j) {
                    exon_starts_[off + j] = nb::cast<int32_t>(starts_t[j]);
                    exon_ends_[off + j]   = nb::cast<int32_t>(ends_t[j]);
                    exon_cumsum_[off + j]  = nb::cast<int32_t>(cum_t[j]);
                }
            }
        }
    }

    // ---------------------------------------------------------------
    // fused_score_buffer — two-pass scoring of entire buffer in C++
    // ---------------------------------------------------------------
    //
    // Replaces the Python accumulation loop in FragmentRouter._scan_native.
    // Pass 1: count EM units + candidates (accumulate unambig counts).
    // Pass 2: fill pre-allocated arrays at exact sizes.
    // Returns capsule-backed numpy arrays (zero-copy to Python).
    //
    // Float64 output for log_liks, coverage_weights, gdna_log_liks.

private:

    // Cached raw pointers for one prepared chunk.
    struct ChunkPtrs {
        const int64_t*  t_off;
        const int32_t*  t_ind;
        const int32_t*  f_len;
        const int32_t*  e_bp;
        const int32_t*  i_bp;
        const uint8_t*  s_type;
        const uint8_t*  e_str;
        const uint8_t*  fc;
        const int64_t*  f_id;
        const uint32_t* r_len;
        const int32_t*  g_fp;
        const int32_t*  g_sta;
        const uint16_t* nm;
        int N;
    };

    // Pre-allocated output arrays + write cursors.
    struct FillState {
        // CSR candidate arrays
        int32_t*  ti;
        double*   ll;
        uint8_t*  ct;
        double*   cw;
        int32_t*  ts;
        int32_t*  te;
        // CSR offsets (n_units + 1)
        int64_t*  offsets;
        // Per-unit metadata
        int32_t*  locus_t;
        uint8_t*  locus_ct;
        int8_t*   is_spliced;
        double*   gdna_ll;
        int32_t*  gfp;
        int32_t*  fid;
        int8_t*   fclass;
        uint8_t*  stype;
        // Det-unambig
        int32_t*  det_ti;
        int64_t*  det_fid;
        // Cursors
        int64_t unit_cur;
        int64_t cand_cur;
        int64_t det_cur;

        // Pending multimapper group state (persists across chunks)
        int64_t mm_fid;     // frag_id of current pending group (-1 = none)
        // Members stored as (chunk_index, row_index) pairs
        std::vector<std::pair<int, int>> mm_members;
        // Pointer to all chunk data (set before scoring begins)
        const std::vector<ChunkPtrs>* all_cps;

        FillState()
            : ti(nullptr), ll(nullptr), ct(nullptr), cw(nullptr),
              ts(nullptr), te(nullptr), offsets(nullptr),
              locus_t(nullptr), locus_ct(nullptr), is_spliced(nullptr),
              gdna_ll(nullptr), gfp(nullptr), fid(nullptr),
              fclass(nullptr), stype(nullptr),
              det_ti(nullptr), det_fid(nullptr),
              unit_cur(0), cand_cur(0), det_cur(0),
              mm_fid(-1), all_cps(nullptr) {}
    };

    // Merged candidate entry for MM cross-alignment merge.
    struct MergedMrna {
        int32_t overhang;
        double  log_lik;
        int32_t count_col;
        double  coverage_wt;
        int32_t tx_start, tx_end;
    };
    // ---------------------------------------------------------------
    // flush_mm_group — score and emit one multimapper group
    // ---------------------------------------------------------------
    //
    // Members are identified by (chunk_idx, row_idx) pairs stored in
    // st.mm_members.  Uses pool-separated likelihood pruning (same as
    // the non-MM path), but merges across all alignments in the group
    // before applying the pruning gate.
    //
    // Mirrors the Python _flush_mm_group logic exactly.

    template<bool FillMode>
    void flush_mm_group(
        FillState& st,
        double gdna_log_sp,
        int64_t& stat_mm,
        int64_t& stat_gated) const
    {
        static constexpr double NEG_INF =
            -std::numeric_limits<double>::infinity();

        const auto& members = st.mm_members;
        const auto& cps = *st.all_cps;
        if (members.empty()) return;

        // --- Determine splice status across all hits ---
        bool is_any_spliced = false;
        for (auto& [ci, ri] : members) {
            int stype = cps[ci].s_type[ri];
            if (stype == SPLICE_SPLICED_ANNOT ||
                stype == SPLICE_SPLICED_UNANNOT) {
                is_any_spliced = true;
                break;
            }
        }

        // --- mRNA pool: score each hit, merge across alignments ---
        std::unordered_map<int32_t, MergedMrna> merged_mrna;

        for (auto& [ci, ri] : members) {
            const auto& cp = cps[ci];
            int64_t start = cp.t_off[ri];
            int64_t end   = cp.t_off[ri + 1];
            int n_cand    = static_cast<int>(end - start);
            if (n_cand <= 0) continue;

            int stype    = cp.s_type[ri];
            int exon_str = cp.e_str[ri];
            int rl = cp.r_len[ri] > 0
                   ? static_cast<int>(cp.r_len[ri]) : 1;
            int nm = cp.nm[ri];
            double log_nm = nm > 0 ? nm * mm_log_pen_ : 0.0;
            bool has_strand  = (exon_str == 1 || exon_str == 2);
            int genomic_start_val = cp.g_sta[ri];
            bool has_genomic = (genomic_start_val >= 0);

            for (int64_t k = start; k < end; ++k) {
                int32_t t_idx = cp.t_ind[k];
                int32_t ebp   = cp.e_bp[k];
                if (ebp <= 0) continue;

                int32_t oh = rl - ebp;
                if (oh < 0) oh = 0;

                int32_t flen = cp.f_len[k];
                double log_fl = frag_len_log_lik(flen);

                double log_strand;
                bool is_anti;
                if (has_strand) {
                    bool same = (exon_str ==
                        static_cast<int>(t_strand_[t_idx]));
                    log_strand = same ? log_p_sense_
                                     : log_p_antisense_;
                    is_anti = same ? r1_antisense_
                                  : !r1_antisense_;
                } else {
                    log_strand = LOG_HALF;
                    is_anti = false;
                }

                double log_lik = log_strand + log_fl
                               + oh * oh_log_pen_ + log_nm;
                int32_t ct = stype * 2 + (is_anti ? 1 : 0);

                // Coverage weight + transcript position
                int32_t t_len = t_length_[t_idx];
                int32_t tx_s = 0;
                int32_t tx_e = flen > 0 ? flen : t_len;
                double cov_wt = 1.0;

                if (has_genomic && flen > 0) {
                    int32_t n_exons =
                        exon_offsets_[t_idx + 1]
                        - exon_offsets_[t_idx];
                    if (n_exons > 0 && t_len > 0) {
                        tx_s = genomic_to_tx_pos(
                            genomic_start_val, t_idx);
                        tx_e = tx_s + flen;
                        int32_t cov_end =
                            tx_e < t_len ? tx_e : t_len;
                        cov_wt = compute_fragment_weight(
                            tx_s, cov_end, t_len);
                    } else {
                        tx_s = 0;
                        tx_e = flen;
                    }
                }

                // Merge: keep best per transcript
                // (lower overhang wins; ties broken by higher log_lik)
                auto it = merged_mrna.find(t_idx);
                if (it == merged_mrna.end()) {
                    merged_mrna[t_idx] = {oh, log_lik, ct,
                                          cov_wt, tx_s, tx_e};
                } else {
                    auto& prev = it->second;
                    if (oh < prev.overhang ||
                        (oh == prev.overhang &&
                         log_lik > prev.log_lik))
                    {
                        prev = {oh, log_lik, ct,
                                cov_wt, tx_s, tx_e};
                    }
                }
            }
        }

        // Pool-separated likelihood pruning (multimapper)
        double mrna_best = NEG_INF;
        double nrna_best = NEG_INF;
        for (auto& [t_idx, mc] : merged_mrna) {
            if (t_is_nrna_[t_idx])
                nrna_best = std::max(nrna_best, mc.log_lik);
            else
                mrna_best = std::max(mrna_best, mc.log_lik);
        }

        int64_t emit_start = st.cand_cur;
        double best_ll = NEG_INF;
        int32_t best_t = -1;
        int32_t best_ct = 0;

        for (auto& [t_idx, mc] : merged_mrna) {
            double pool_best = t_is_nrna_[t_idx]
                             ? nrna_best : mrna_best;
            if (pool_best - mc.log_lik <= max_ll_delta_) {
                if constexpr (FillMode) {
                    int64_t c = st.cand_cur;
                    st.ti[c] = t_idx;
                    st.ll[c] = mc.log_lik;
                    st.ct[c] = static_cast<uint8_t>(mc.count_col);
                    st.cw[c] = mc.coverage_wt;
                    st.ts[c] = mc.tx_start;
                    st.te[c] = mc.tx_end;
                }
                ++st.cand_cur;
                if (mc.log_lik > best_ll) {
                    best_ll = mc.log_lik;
                    best_t  = t_idx;
                    best_ct = mc.count_col;
                }
            }
        }

        // --- Finalize this MM EM unit ---
        int64_t new_cands = st.cand_cur - emit_start;
        if (new_cands > 0) {
            if constexpr (FillMode) {
                st.offsets[st.unit_cur + 1] = st.cand_cur;
                st.locus_t[st.unit_cur] = best_t;
                st.locus_ct[st.unit_cur] =
                    static_cast<uint8_t>(best_ct);
                st.fid[st.unit_cur] =
                    static_cast<int32_t>(st.mm_fid);
                st.fclass[st.unit_cur] =
                    static_cast<int8_t>(FRAG_MULTIMAPPER);

                // Splice type: ANNOT > UNANNOT > UNSPLICED
                int mm_stype = SPLICE_UNSPLICED;
                for (auto& [ci, ri] : members) {
                    int stype = cps[ci].s_type[ri];
                    if (stype == SPLICE_SPLICED_ANNOT) {
                        mm_stype = SPLICE_SPLICED_ANNOT;
                        break;
                    } else if (stype == SPLICE_SPLICED_UNANNOT) {
                        mm_stype = SPLICE_SPLICED_UNANNOT;
                    }
                }
                st.stype[st.unit_cur] =
                    static_cast<uint8_t>(mm_stype);
                st.is_spliced[st.unit_cur] =
                    is_any_spliced ? 1 : 0;

                if (!is_any_spliced) {
                    // gDNA log-lik: best across unspliced hits
                    double best_gdna_ll = NEG_INF;
                    int32_t best_fp = 0;
                    for (auto& [ci, ri] : members) {
                        int stype = cps[ci].s_type[ri];
                        if (stype == SPLICE_SPLICED_ANNOT ||
                            stype == SPLICE_SPLICED_UNANNOT)
                            continue;
                        int32_t gfp_val = cps[ci].g_fp[ri];
                        int nm_val = cps[ci].nm[ri];
                        double hit_log_nm = nm_val > 0
                            ? nm_val * mm_log_pen_ : 0.0;
                        double gdna_fl =
                            gdna_frag_len_log_lik(gfp_val);
                        double gdna_ll_val =
                            gdna_fl + gdna_log_sp + LOG_HALF + hit_log_nm;
                        if (gdna_ll_val > best_gdna_ll) {
                            best_gdna_ll = gdna_ll_val;
                            best_fp = gfp_val;
                        }
                    }
                    st.gdna_ll[st.unit_cur] = best_gdna_ll;
                    st.gfp[st.unit_cur] = best_fp;
                } else {
                    st.gdna_ll[st.unit_cur] = NEG_INF;
                    // Use first hit's footprint
                    auto& [ci0, ri0] = members[0];
                    st.gfp[st.unit_cur] = cps[ci0].g_fp[ri0];
                }
            }
            ++st.unit_cur;
            ++stat_mm;
        } else {
            st.cand_cur = emit_start;
            ++stat_gated;
        }
    }

    // Templated inner loop — shared by count and fill passes.
    // FillMode=false: count only, update accumulators (out is unused).
    // FillMode=true:  write into pre-allocated arrays, skip accumulators.
    template<bool FillMode>
    void score_chunk_impl(
        const ChunkPtrs& cp,
        int chunk_idx,
        const int8_t* t_str,
        double* ua_counts, int ua_ncols,
        double gdna_log_sp,
        FillState& st,
        int64_t& stat_det, int64_t& stat_em_u,
        int64_t& stat_em_as, int64_t& stat_em_ao,
        int64_t& stat_gated, int64_t& stat_chim,
        int64_t& stat_mm) const
    {
        static constexpr double NEG_INF =
            -std::numeric_limits<double>::infinity();

        const int N            = cp.N;
        const int64_t*  t_off  = cp.t_off;
        const int32_t*  t_ind  = cp.t_ind;
        const int32_t*  f_len  = cp.f_len;
        const int32_t*  e_bp   = cp.e_bp;
        const int32_t*  i_bp   = cp.i_bp;
        const uint8_t*  s_type = cp.s_type;
        const uint8_t*  e_str  = cp.e_str;
        const uint8_t*  fc     = cp.fc;
        const int64_t*  f_id   = cp.f_id;
        const uint32_t* r_len  = cp.r_len;
        const int32_t*  g_fp   = cp.g_fp;
        const int32_t*  g_sta  = cp.g_sta;
        const uint16_t* nm_arr = cp.nm;

        for (int i = 0; i < N; ++i) {
            int fclass = fc[i];

            if (fclass == FRAG_CHIMERIC) { ++stat_chim; continue; }

            // ---- Multimapper accumulation ----
            if (fclass == FRAG_MULTIMAPPER) {
                int64_t fid_val = f_id[i];
                if (fid_val != st.mm_fid) {
                    // New group — flush previous if any
                    if (!st.mm_members.empty()) {
                        flush_mm_group<FillMode>(
                            st, gdna_log_sp,
                            stat_mm, stat_gated);
                    }
                    st.mm_fid = fid_val;
                    st.mm_members.clear();
                }
                st.mm_members.emplace_back(chunk_idx, i);
                continue;
            }

            // Non-MM fragment: flush any pending MM group
            if (!st.mm_members.empty()) {
                flush_mm_group<FillMode>(
                    st, gdna_log_sp,
                    stat_mm, stat_gated);
                st.mm_members.clear();
                st.mm_fid = -1;
            }

            int stype      = s_type[i];
            int exon_str   = e_str[i];
            int64_t start  = t_off[i];
            int64_t end    = t_off[i + 1];
            int n_cand     = static_cast<int>(end - start);

            // ---- Pre-EM strand accumulation (count pass only) ----
            if constexpr (!FillMode) {
                if (fclass == FRAG_UNAMBIG ||
                    fclass == FRAG_AMBIG_SAME_STRAND)
                {
                    // Det-unambig: SPLICE_ANNOT + FRAG_UNAMBIG
                    if (fclass == FRAG_UNAMBIG &&
                        stype == SPLICE_SPLICED_ANNOT) {
                        int32_t t_idx = t_ind[start];
                        int t_sv =
                            static_cast<int>(t_str[t_idx]);
                        bool anti;
                        if (exon_str == 1 || exon_str == 2) {
                            bool same = (exon_str == t_sv);
                            double lp = same ? log_p_sense_
                                            : log_p_antisense_;
                            anti = std::exp(lp) < 0.5;
                        } else {
                            anti = false;
                        }
                        int col =
                            stype * 2 + (anti ? 1 : 0);
                        if (col < ua_ncols)
                            ua_counts[t_idx * ua_ncols + col]
                                += 1.0;
                        ++stat_det;
                        ++st.det_cur;
                        continue;
                    }
                }
            } else {
                // Fill pass: skip accum, still skip det-unambig
                if (fclass == FRAG_UNAMBIG ||
                    fclass == FRAG_AMBIG_SAME_STRAND)
                {
                    if (fclass == FRAG_UNAMBIG &&
                        stype == SPLICE_SPLICED_ANNOT)
                    {
                        // Record det-unambig for annotations
                        st.det_ti[st.det_cur] = t_ind[start];
                        st.det_fid[st.det_cur] =
                            f_id[i];
                        ++st.det_cur;
                        ++stat_det;
                        continue;
                    }
                }
            }

            // ---- EM-routed scoring ----
            int rl = r_len[i] > 0 ? static_cast<int>(r_len[i]) : 1;
            int nm = nm_arr[i];
            double log_nm = nm > 0 ? nm * mm_log_pen_ : 0.0;
            bool has_strand  = (exon_str == 1 || exon_str == 2);
            int genomic_start = g_sta[i];
            int genomic_footprint = g_fp[i];
            bool has_genomic = (genomic_start >= 0);

            int64_t emit_start = st.cand_cur;
            double best_ll = NEG_INF;
            int32_t best_t = -1;
            int32_t best_ct = 0;

            if (n_cand > 0) {
                // ===== mRNA scoring =====
                struct MrnaScored {
                    int32_t t_idx, oh, ct, tx_s, tx_e;
                    double  log_lik, cov_wt;
                };
                MrnaScored m_stack[SCORED_STACK_CAPACITY];
                std::vector<MrnaScored> m_heap;
                MrnaScored* m_scored = m_stack;
                if (n_cand > SCORED_STACK_CAPACITY) {
                    m_heap.resize(n_cand);
                    m_scored = m_heap.data();
                }
                int m_n = 0;

                for (int64_t k = start; k < end; ++k) {
                    int32_t t_idx = t_ind[k];
                    int32_t ebp   = e_bp[k];
                    if (ebp <= 0) continue;

                    int32_t oh = rl - ebp;
                    if (oh < 0) oh = 0;

                    int32_t flen = f_len[k];
                    double log_fl = frag_len_log_lik(flen);

                    double log_strand;
                    bool is_anti;
                    if (has_strand) {
                        bool same = (exon_str ==
                            static_cast<int>(t_strand_[t_idx]));
                        log_strand = same ? log_p_sense_
                                         : log_p_antisense_;
                        is_anti = same ? r1_antisense_
                                      : !r1_antisense_;
                    } else {
                        log_strand = LOG_HALF;
                        is_anti = false;
                    }

                    double log_lik = log_strand + log_fl
                                   + oh * oh_log_pen_ + log_nm;
                    int32_t ct = stype * 2 + (is_anti ? 1 : 0);

                    int32_t t_len = t_length_[t_idx];
                    int32_t tx_s = 0;
                    int32_t tx_e = flen > 0 ? flen : t_len;
                    double cov_wt = 1.0;

                    if (has_genomic && flen > 0) {
                        int32_t n_exons =
                            exon_offsets_[t_idx + 1]
                            - exon_offsets_[t_idx];
                        if (n_exons > 0 && t_len > 0) {
                            tx_s = genomic_to_tx_pos(
                                genomic_start, t_idx);
                            tx_e = tx_s + flen;
                            int32_t cov_end =
                                tx_e < t_len ? tx_e : t_len;
                            cov_wt = compute_fragment_weight(
                                tx_s, cov_end, t_len);
                        } else {
                            tx_s = 0;
                            tx_e = flen;
                        }
                    }

                    m_scored[m_n++] = {t_idx, oh, ct, tx_s,
                                       tx_e, log_lik, cov_wt};
                }

                // Pool-separated likelihood pruning:
                // Find best log_lik within each pool (mRNA vs nRNA)
                double mrna_best = NEG_INF;
                double nrna_best = NEG_INF;
                for (int j = 0; j < m_n; ++j) {
                    if (t_is_nrna_[m_scored[j].t_idx])
                        nrna_best = std::max(nrna_best,
                                             m_scored[j].log_lik);
                    else
                        mrna_best = std::max(mrna_best,
                                             m_scored[j].log_lik);
                }

                // Emit candidates within Δ of their pool's best
                for (int j = 0; j < m_n; ++j) {
                    auto& s = m_scored[j];
                    double pool_best = t_is_nrna_[s.t_idx]
                                     ? nrna_best : mrna_best;
                    if (pool_best - s.log_lik <= max_ll_delta_) {
                        if constexpr (FillMode) {
                            int64_t c = st.cand_cur;
                            st.ti[c] = s.t_idx;
                            st.ll[c] = s.log_lik;
                            st.ct[c] = static_cast<uint8_t>(
                                s.ct);
                            st.cw[c] = s.cov_wt;
                            st.ts[c] = s.tx_s;
                            st.te[c] = s.tx_e;
                        }
                        ++st.cand_cur;
                        if (s.log_lik > best_ll) {
                            best_ll = s.log_lik;
                            best_t  = s.t_idx;
                            best_ct = s.ct;
                        }
                    }
                }
            }

            // ---- Finalize this EM unit ----
            int64_t new_cands = st.cand_cur - emit_start;
            if (new_cands > 0) {
                if constexpr (FillMode) {
                    st.offsets[st.unit_cur + 1] = st.cand_cur;
                    st.locus_t[st.unit_cur] = best_t;
                    st.locus_ct[st.unit_cur] =
                        static_cast<uint8_t>(best_ct);
                    st.fid[st.unit_cur] =
                        static_cast<int32_t>(f_id[i]);
                    st.fclass[st.unit_cur] =
                        static_cast<int8_t>(fclass);
                    st.stype[st.unit_cur] =
                        static_cast<uint8_t>(stype);

                    bool is_spl =
                        (stype == SPLICE_SPLICED_ANNOT
                         || stype == SPLICE_SPLICED_UNANNOT);
                    st.is_spliced[st.unit_cur] =
                        is_spl ? 1 : 0;
                    st.gfp[st.unit_cur] = genomic_footprint;

                    if (!is_spl) {
                        double gdna_fl =
                            gdna_frag_len_log_lik(
                                genomic_footprint);
                        st.gdna_ll[st.unit_cur] =
                                gdna_fl + gdna_log_sp + LOG_HALF + log_nm;
                    } else {
                        st.gdna_ll[st.unit_cur] = NEG_INF;
                    }
                }
                ++st.unit_cur;

                if (fclass == FRAG_UNAMBIG)
                    ++stat_em_u;
                else if (fclass == FRAG_AMBIG_SAME_STRAND)
                    ++stat_em_as;
                else
                    ++stat_em_ao;
            } else {
                // Revert cand_cur (no candidates emitted)
                st.cand_cur = emit_start;
                ++stat_gated;
            }
        }
    }

public:

    nb::tuple fused_score_buffer(
        nb::list chunk_arrays,
        i8_1d   t_to_strand_arr,
        f64_2d_mut unambig_counts,
        double  gdna_log_splice_pen_unspliced)
    {
        const int8_t* t_str = t_to_strand_arr.data();
        double* ua_c   = unambig_counts.data();
        int ua_ncols   =
            static_cast<int>(unambig_counts.shape(1));

        size_t n_chunks = chunk_arrays.size();

        // ---- Extract and cache chunk pointers ----
        // The Python list keeps all tuples alive; tuples keep
        // the numpy arrays alive, so raw pointers stay valid.
        std::vector<ChunkPtrs> cps(n_chunks);
        for (size_t ci = 0; ci < n_chunks; ++ci) {
            nb::tuple c =
                nb::borrow<nb::tuple>(chunk_arrays[ci]);
            auto& cp  = cps[ci];
            cp.t_off  = nb::cast<i64_1d>(c[0]).data();
            cp.t_ind  = nb::cast<i32_1d>(c[1]).data();
            cp.f_len  = nb::cast<i32_1d>(c[2]).data();
            cp.e_bp   = nb::cast<i32_1d>(c[3]).data();
            cp.i_bp   = nb::cast<i32_1d>(c[4]).data();
            cp.s_type = nb::cast<u8_1d>(c[5]).data();
            cp.e_str  = nb::cast<u8_1d>(c[6]).data();
            cp.fc     = nb::cast<u8_1d>(c[7]).data();
            cp.f_id   = nb::cast<i64_1d>(c[8]).data();
            cp.r_len  = nb::cast<u32_1d>(c[9]).data();
            cp.g_fp   = nb::cast<i32_1d>(c[10]).data();
            cp.g_sta  = nb::cast<i32_1d>(c[11]).data();
            cp.nm     = nb::cast<u16_1d>(c[12]).data();
            cp.N      = static_cast<int>(
                nb::cast<u8_1d>(c[5]).shape(0));
        }

        // ============ PASS 1: COUNT ============
        FillState count_st{};
        count_st.all_cps = &cps;
        int64_t stat_det = 0, stat_em_u = 0;
        int64_t stat_em_as = 0, stat_em_ao = 0;
        int64_t stat_gated = 0, stat_chim = 0;
        int64_t stat_mm = 0;

        for (size_t ci = 0; ci < n_chunks; ++ci) {
            score_chunk_impl<false>(
                cps[ci],
                static_cast<int>(ci),
                t_str,
                ua_c, ua_ncols,
                gdna_log_splice_pen_unspliced,
                count_st,
                stat_det, stat_em_u, stat_em_as,
                stat_em_ao, stat_gated, stat_chim,
                stat_mm);
        }
        // Flush final pending MM group from count pass
        if (!count_st.mm_members.empty()) {
            flush_mm_group<false>(
                count_st,
                gdna_log_splice_pen_unspliced,
                stat_mm, stat_gated);
            count_st.mm_members.clear();
            count_st.mm_fid = -1;
        }

        int64_t total_units = count_st.unit_cur;
        int64_t total_cands = count_st.cand_cur;
        int64_t total_det   = count_st.det_cur;

        // ============ ALLOCATE ============
        auto* v_offsets = new std::vector<int64_t>(
            total_units + 1, 0);
        auto* v_ti  = new std::vector<int32_t>(total_cands);
        auto* v_ll  = new std::vector<double>(total_cands);
        auto* v_ct  = new std::vector<uint8_t>(total_cands);
        auto* v_cw  = new std::vector<double>(total_cands);
        auto* v_ts  = new std::vector<int32_t>(total_cands);
        auto* v_te  = new std::vector<int32_t>(total_cands);
        auto* v_lt  = new std::vector<int32_t>(total_units);
        auto* v_lct = new std::vector<uint8_t>(total_units);
        auto* v_isp = new std::vector<int8_t>(total_units);
        auto* v_gll = new std::vector<double>(total_units);
        auto* v_gfp = new std::vector<int32_t>(total_units);
        auto* v_fid = new std::vector<int32_t>(total_units);
        auto* v_fc  = new std::vector<int8_t>(total_units);
        auto* v_st  = new std::vector<uint8_t>(total_units);
        auto* v_dti = new std::vector<int32_t>(total_det);
        auto* v_dfid = new std::vector<int64_t>(total_det);

        // ============ PASS 2: FILL ============
        FillState fill_st{};
        fill_st.offsets   = v_offsets->data();
        fill_st.ti        = v_ti->data();
        fill_st.ll        = v_ll->data();
        fill_st.ct        = v_ct->data();
        fill_st.cw        = v_cw->data();
        fill_st.ts        = v_ts->data();
        fill_st.te        = v_te->data();
        fill_st.locus_t   = v_lt->data();
        fill_st.locus_ct  = v_lct->data();
        fill_st.is_spliced = v_isp->data();
        fill_st.gdna_ll   = v_gll->data();
        fill_st.gfp       = v_gfp->data();
        fill_st.fid       = v_fid->data();
        fill_st.fclass    = v_fc->data();
        fill_st.stype     = v_st->data();
        fill_st.det_ti    = v_dti->data();
        fill_st.det_fid   = v_dfid->data();

        // Reset stats for fill pass (should match)
        fill_st.all_cps = &cps;
        int64_t f_det = 0, f_em_u = 0, f_em_as = 0;
        int64_t f_em_ao = 0, f_gated = 0, f_chim = 0;
        int64_t f_mm = 0;

        for (size_t ci = 0; ci < n_chunks; ++ci) {
            score_chunk_impl<true>(
                cps[ci],
                static_cast<int>(ci),
                t_str,
                ua_c, ua_ncols,
                gdna_log_splice_pen_unspliced,
                fill_st,
                f_det, f_em_u, f_em_as,
                f_em_ao, f_gated, f_chim,
                f_mm);
        }
        // Flush final pending MM group from fill pass
        if (!fill_st.mm_members.empty()) {
            flush_mm_group<true>(
                fill_st,
                gdna_log_splice_pen_unspliced,
                f_mm, f_gated);
            fill_st.mm_members.clear();
            fill_st.mm_fid = -1;
        }

        // ============ RETURN NUMPY CAPSULES ============
        auto mk_i64 = [](std::vector<int64_t>* v)
            -> nb::object {
            size_t n = v->size();
            nb::capsule del(v, [](void* p) noexcept {
                delete static_cast<
                    std::vector<int64_t>*>(p);
            });
            return nb::ndarray<nb::numpy, int64_t,
                               nb::ndim<1>>(
                v->data(), {n}, del).cast();
        };
        auto mk_i32 = [](std::vector<int32_t>* v)
            -> nb::object {
            size_t n = v->size();
            nb::capsule del(v, [](void* p) noexcept {
                delete static_cast<
                    std::vector<int32_t>*>(p);
            });
            return nb::ndarray<nb::numpy, int32_t,
                               nb::ndim<1>>(
                v->data(), {n}, del).cast();
        };
        auto mk_f64 = [](std::vector<double>* v)
            -> nb::object {
            size_t n = v->size();
            nb::capsule del(v, [](void* p) noexcept {
                delete static_cast<
                    std::vector<double>*>(p);
            });
            return nb::ndarray<nb::numpy, double,
                               nb::ndim<1>>(
                v->data(), {n}, del).cast();
        };
        auto mk_u8 = [](std::vector<uint8_t>* v)
            -> nb::object {
            size_t n = v->size();
            nb::capsule del(v, [](void* p) noexcept {
                delete static_cast<
                    std::vector<uint8_t>*>(p);
            });
            return nb::ndarray<nb::numpy, uint8_t,
                               nb::ndim<1>>(
                v->data(), {n}, del).cast();
        };
        auto mk_i8 = [](std::vector<int8_t>* v)
            -> nb::object {
            size_t n = v->size();
            nb::capsule del(v, [](void* p) noexcept {
                delete static_cast<
                    std::vector<int8_t>*>(p);
            });
            return nb::ndarray<nb::numpy, int8_t,
                               nb::ndim<1>>(
                v->data(), {n}, del).cast();
        };

        return nb::make_tuple(
            // CSR arrays
            mk_i64(v_offsets),
            mk_i32(v_ti),
            mk_f64(v_ll),
            mk_u8(v_ct),
            mk_f64(v_cw),
            mk_i32(v_ts),
            mk_i32(v_te),
            // Per-unit metadata
            mk_i32(v_lt),
            mk_u8(v_lct),
            mk_i8(v_isp),
            mk_f64(v_gll),
            mk_i32(v_gfp),
            mk_i32(v_fid),
            mk_i8(v_fc),
            mk_u8(v_st),
            // Det-unambig
            mk_i32(v_dti),
            mk_i64(v_dfid),
            // Stats
            stat_det,
            stat_em_u,
            stat_em_as,
            stat_em_ao,
            stat_gated,
            stat_chim,
            stat_mm
        );
    }

};

// ----------------------------------------------------------------
// nanobind module definition
// ----------------------------------------------------------------

NB_MODULE(_scoring_impl, m) {
    m.doc() = "C++ hot-path scoring kernels for rigel (nanobind)";

    m.def("compute_fragment_weight", &compute_fragment_weight,
          nb::arg("frag_start"), nb::arg("frag_end"),
          nb::arg("transcript_length"),
          "Inverse coverage-capacity weight for a fragment on a transcript.\n\n"
          "Uses the trapezoid coverage model.  Plateau → 1.0; edge → > 1.0.");

    nb::class_<NativeFragmentScorer>(m, "NativeFragmentScorer")
        .def(nb::init<
                 double, double, bool, double, double,
                 nb::object, int32_t, double,
                 nb::object, int32_t, double,
                 i8_1d, i32_1d, i32_1d, i32_1d,
                 nb::object,
                 u8_1d, double>(),
             nb::arg("log_p_sense"),
             nb::arg("log_p_antisense"),
             nb::arg("r1_antisense"),
             nb::arg("overhang_log_penalty"),
             nb::arg("mismatch_log_penalty"),
             nb::arg("fl_log_prob").none(),
             nb::arg("fl_max_size"),
             nb::arg("fl_tail_base"),
             nb::arg("gdna_fl_log_prob").none(),
             nb::arg("gdna_fl_max_size"),
             nb::arg("gdna_fl_tail_base"),
             nb::arg("t_strand_arr"),
             nb::arg("t_length_arr"),
             nb::arg("t_span_arr"),
             nb::arg("t_start_arr"),
             nb::arg("t_exon_data").none(),
             nb::arg("t_is_nrna_arr"),
             nb::arg("pruning_max_ll_delta"))
        .def("fused_score_buffer",
             &NativeFragmentScorer::fused_score_buffer,
             nb::arg("chunk_arrays"),
             nb::arg("t_to_strand_arr"),
             nb::arg("unambig_counts"),
             nb::arg("gdna_log_splice_pen_unspliced"));

    // Export scoring constants for Python-side parity tests
    m.attr("LOG_HALF")      = rigel::LOG_HALF;
    m.attr("TAIL_DECAY_LP") = rigel::TAIL_DECAY_LP;

    // Fragment classification constants (mirrors rigel.buffer)
    m.attr("FRAG_UNAMBIG")           = rigel::FRAG_UNAMBIG;
    m.attr("FRAG_AMBIG_SAME_STRAND") = rigel::FRAG_AMBIG_SAME_STRAND;
    m.attr("FRAG_AMBIG_OPP_STRAND")  = rigel::FRAG_AMBIG_OPP_STRAND;
    m.attr("FRAG_MULTIMAPPER")       = rigel::FRAG_MULTIMAPPER;
    m.attr("FRAG_CHIMERIC")          = rigel::FRAG_CHIMERIC;

    // Pool code constants (mirrors rigel.annotate)
    m.attr("POOL_CODE_MRNA")       = rigel::POOL_CODE_MRNA;
    m.attr("POOL_CODE_NRNA")       = rigel::POOL_CODE_NRNA;
    m.attr("POOL_CODE_GDNA")       = rigel::POOL_CODE_GDNA;
    m.attr("POOL_CODE_INTERGENIC") = rigel::POOL_CODE_INTERGENIC;
    m.attr("POOL_CODE_CHIMERIC")   = rigel::POOL_CODE_CHIMERIC;

    // Internal tuning constant
    m.attr("SCORED_STACK_CAPACITY") = SCORED_STACK_CAPACITY;
}
