/**
 * scoring.cpp — nanobind C++ extension for hot-path scoring functions.
 *
 * This module provides performance-critical kernels called from the
 * Python scoring loop.  Functions are exposed as ``rigel._scoring_impl``.
 *
 * Contents:
 *   - compute_fragment_weight(frag_start, frag_end, transcript_length)
 *   - NativeFragmentScorer class:
 *       .score_wta_mrna(...)         — mRNA WTA scoring kernel
 *       .score_wta_nrna(...)         — nRNA WTA scoring kernel
 *       .score_emit_fragment(...)    — fused score+emit (bytes output)
 *
 * Build:
 *   Part of the rigel scikit-build-core build — see CMakeLists.txt.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "constants.h"

namespace nb = nanobind;

// Use constants from constants.h via the rigel namespace.
using rigel::LOG_HALF;
using rigel::TAIL_DECAY_LP;
using rigel::STRAND_NEG;

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
// score_wta_mrna / score_wta_nrna can execute entirely in C++ with
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
    int32_t nrna_base_;
    int32_t n_transcripts_;

    // --- Copied index arrays ---
    std::vector<int8_t>  t_strand_;   // int8[n_transcripts]
    std::vector<int32_t> t_length_;   // spliced exonic length
    std::vector<int32_t> t_span_;     // genomic span (incl introns)
    std::vector<int32_t> t_start_;    // genomic start coordinate

    // --- Fragment-length LUT (copied from numpy) ---
    std::vector<double> fl_log_prob_;

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
        i8_1d   t_strand_arr,
        i32_1d  t_length_arr,
        i32_1d  t_span_arr,
        i32_1d  t_start_arr,
        int32_t nrna_base,
        nb::object t_exon_data_obj)
      : log_p_sense_(log_p_sense),
        log_p_antisense_(log_p_antisense),
        r1_antisense_(r1_antisense),
        oh_log_pen_(oh_log_pen),
        mm_log_pen_(mm_log_pen),
        fl_max_size_(fl_max_size),
        fl_tail_base_(fl_tail_base),
        has_fl_lut_(false),
        nrna_base_(nrna_base)
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

        // Copy fragment-length LUT
        if (!fl_log_prob_obj.is_none()) {
            auto fl_arr = nb::cast<f64_1d>(fl_log_prob_obj);
            const double* p = fl_arr.data();
            int32_t n = static_cast<int32_t>(fl_arr.shape(0));
            fl_log_prob_.assign(p, p + n);
            has_fl_lut_ = true;
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
    // score_wta_mrna — mRNA winner-take-all scoring kernel
    // ---------------------------------------------------------------

    nb::dict score_wta_mrna(
        i32_1d t_inds,
        i32_1d exon_bp,
        i32_1d frag_lengths,
        int exon_strand,
        int splice_type,
        int nm,
        int read_length,
        int genomic_start)
    {
        int n_cand = static_cast<int>(t_inds.shape(0));
        if (n_cand == 0) return nb::dict();

        const int32_t* ti = t_inds.data();
        const int32_t* eb = exon_bp.data();
        const int32_t* fl = frag_lengths.data();

        int rl = read_length > 0 ? read_length : 1;
        double log_nm = nm > 0 ? nm * mm_log_pen_ : 0.0;
        bool has_strand  = (exon_strand == 1 || exon_strand == 2);
        bool has_genomic = (genomic_start >= 0);

        // Scratch storage (stack for typical small candidate sets)
        struct Scored {
            int32_t t_idx, oh, ct, tx_s, tx_e;
            double  log_lik, cov_wt;
        };
        Scored stack_buf[64];
        std::vector<Scored> heap_buf;
        Scored* scored = stack_buf;
        if (n_cand > 64) {
            heap_buf.resize(n_cand);
            scored = heap_buf.data();
        }

        int n_scored = 0;
        int32_t min_oh = 0x7FFFFFFF;

        for (int k = 0; k < n_cand; ++k) {
            int32_t t_idx = ti[k];
            int32_t ebp   = eb[k];
            if (ebp <= 0) continue;

            int32_t oh = rl - ebp;
            if (oh < 0) oh = 0;

            // Fragment-length log-likelihood (scalar LUT)
            int32_t flen = fl[k];
            double log_fl = frag_len_log_lik(flen);

            // Strand log-likelihood
            double log_strand;
            bool is_anti;
            if (has_strand) {
                bool same = (exon_strand ==
                             static_cast<int>(t_strand_[t_idx]));
                log_strand = same ? log_p_sense_ : log_p_antisense_;
                is_anti    = same ? r1_antisense_ : !r1_antisense_;
            } else {
                log_strand = LOG_HALF;
                is_anti = false;
            }

            double log_lik = log_strand + log_fl
                           + oh * oh_log_pen_ + log_nm;
            int32_t ct = splice_type * 2 + (is_anti ? 1 : 0);

            // Coverage weight + transcript position
            int32_t t_len = t_length_[t_idx];
            int32_t tx_s = 0;
            int32_t tx_e = flen > 0 ? flen : t_len;
            double cov_wt = 1.0;

            if (has_genomic && flen > 0) {
                int32_t n_exons = exon_offsets_[t_idx + 1]
                                - exon_offsets_[t_idx];
                if (n_exons > 0 && t_len > 0) {
                    tx_s = genomic_to_tx_pos(genomic_start, t_idx);
                    tx_e = tx_s + flen;
                    int32_t cov_end = tx_e < t_len ? tx_e : t_len;
                    cov_wt = compute_fragment_weight(tx_s, cov_end, t_len);
                } else {
                    tx_s = 0;
                    tx_e = flen;
                }
            }

            scored[n_scored++] = {t_idx, oh, ct, tx_s, tx_e,
                                  log_lik, cov_wt};
            if (oh < min_oh) min_oh = oh;
        }

        if (n_scored == 0) return nb::dict();

        // WTA gating: keep only minimum-overhang winners
        nb::dict result;
        for (int j = 0; j < n_scored; ++j) {
            auto& s = scored[j];
            if (s.oh == min_oh) {
                result[nb::int_(s.t_idx)] = nb::make_tuple(
                    s.oh, s.log_lik, s.ct, s.cov_wt, s.tx_s, s.tx_e);
            }
        }
        return result;
    }

    // ---------------------------------------------------------------
    // score_wta_nrna — nRNA winner-take-all scoring kernel
    // ---------------------------------------------------------------

    nb::dict score_wta_nrna(
        i32_1d t_inds,
        i32_1d exon_bp,
        i32_1d intron_bp,
        int exon_strand,
        int nm,
        int read_length,
        int genomic_footprint,
        int genomic_start)
    {
        int n_cand = static_cast<int>(t_inds.shape(0));
        if (n_cand == 0) return nb::dict();

        const int32_t* ti = t_inds.data();
        const int32_t* eb = exon_bp.data();
        const int32_t* ib = intron_bp.data();

        int rl = read_length > 0 ? read_length : 1;
        double log_fl = frag_len_log_lik(genomic_footprint);
        double log_nm = nm > 0 ? nm * mm_log_pen_ : 0.0;
        bool has_strand  = (exon_strand == 1 || exon_strand == 2);
        bool has_genomic = (genomic_start >= 0 && genomic_footprint > 0);
        int32_t nrna_gfp = genomic_footprint > 0 ? genomic_footprint : 0;

        struct Scored {
            int32_t t_idx, oh, tx_s, tx_e;
            double  nrna_ll, cov_wt;
        };
        Scored stack_buf[64];
        std::vector<Scored> heap_buf;
        Scored* scored = stack_buf;
        if (n_cand > 64) {
            heap_buf.resize(n_cand);
            scored = heap_buf.data();
        }

        int n_scored = 0;
        int32_t min_oh = 0x7FFFFFFF;

        for (int k = 0; k < n_cand; ++k) {
            int32_t t_idx = ti[k];

            // Single-exon transcript filter
            int32_t t_span   = t_span_[t_idx];
            int32_t t_exonic = t_length_[t_idx];
            if (t_span <= t_exonic) continue;

            int32_t ebp = eb[k];
            int32_t ibp = ib[k];
            int32_t span_bp = ebp + ibp;
            if (span_bp <= 0) continue;

            int32_t oh = rl - span_bp;
            if (oh < 0) oh = 0;

            // Strand log-likelihood
            double log_strand;
            if (has_strand) {
                bool same = (exon_strand ==
                             static_cast<int>(t_strand_[t_idx]));
                log_strand = same ? log_p_sense_ : log_p_antisense_;
            } else {
                log_strand = LOG_HALF;
            }

            double nrna_ll = log_strand + log_fl
                           + oh * oh_log_pen_ + log_nm;

            // Coverage weight (nRNA uses genomic span, not tx coords)
            int32_t tx_s = 0;
            int32_t tx_e = nrna_gfp > 0 ? nrna_gfp : t_span;
            double cov_wt = 1.0;

            if (has_genomic) {
                int32_t tx_start_g   = t_start_[t_idx];
                int32_t frag_pos     = genomic_start - tx_start_g;
                int32_t frag_end_pos = frag_pos + genomic_footprint;

                // Clamped for fragment weight
                int32_t fpc = frag_pos < t_span ? frag_pos : t_span;
                if (fpc < 0) fpc = 0;
                int32_t fec = frag_end_pos < t_span
                            ? frag_end_pos : t_span;
                if (fec < 0) fec = 0;
                if (fec > fpc && t_span > 0) {
                    cov_wt = compute_fragment_weight(fpc, fec, t_span);
                }

                // Unclipped for bias model
                tx_s = frag_pos > 0 ? frag_pos : 0;
                tx_e = tx_s + genomic_footprint;
            }

            scored[n_scored++] = {t_idx, oh, tx_s, tx_e,
                                  nrna_ll, cov_wt};
            if (oh < min_oh) min_oh = oh;
        }

        if (n_scored == 0) return nb::dict();

        // WTA gating
        nb::dict result;
        for (int j = 0; j < n_scored; ++j) {
            auto& s = scored[j];
            if (s.oh == min_oh) {
                result[nb::int_(s.t_idx)] = nb::make_tuple(
                    s.oh, s.nrna_ll, s.cov_wt, s.tx_s, s.tx_e);
            }
        }
        return result;
    }

    // ---------------------------------------------------------------
    // score_emit_fragment — fused score + emit for single fragments
    // ---------------------------------------------------------------
    //
    // Combines mRNA and nRNA WTA scoring with CSR emission into a
    // single C++ call.  Returns (best_ll, best_t, best_ct,
    //   ti_bytes, ll_bytes, ct_bytes, cw_bytes, ts_bytes, te_bytes)
    // where each *_bytes is a raw bytes buffer suitable for
    // array.array.frombytes().

    nb::tuple score_emit_fragment(
        i32_1d t_inds,
        i32_1d exon_bp,
        i32_1d frag_lengths,
        i32_1d intron_bp,
        int exon_strand,
        int splice_type,
        int nm,
        int read_length,
        int genomic_footprint,
        int genomic_start)
    {
        static constexpr int SPLICE_ANNOT_VAL = 2;
        static constexpr double NEG_INF = -std::numeric_limits<double>::infinity();

        int n_cand = static_cast<int>(t_inds.shape(0));

        // Result vectors (small — typically 5-20 winners)
        std::vector<int32_t> out_ti, out_ts, out_te;
        std::vector<double>  out_ll, out_cw;
        std::vector<uint8_t> out_ct;

        double best_ll = NEG_INF;
        int32_t best_t = -1, best_ct = 0;

        if (n_cand == 0)
            goto pack_output;

        {
            const int32_t* ti = t_inds.data();
            const int32_t* eb = exon_bp.data();
            const int32_t* fl = frag_lengths.data();
            const int32_t* ib = intron_bp.data();

            int rl = read_length > 0 ? read_length : 1;
            double log_nm = nm > 0 ? nm * mm_log_pen_ : 0.0;
            bool has_strand  = (exon_strand == 1 || exon_strand == 2);
            bool has_genomic = (genomic_start >= 0);

            // ========== mRNA scoring ==========
            struct MrnaScored {
                int32_t t_idx, oh, ct, tx_s, tx_e;
                double  log_lik, cov_wt;
            };
            MrnaScored m_stack[64];
            std::vector<MrnaScored> m_heap;
            MrnaScored* m_scored = m_stack;
            if (n_cand > 64) {
                m_heap.resize(n_cand);
                m_scored = m_heap.data();
            }
            int m_n = 0;
            int32_t m_min_oh = 0x7FFFFFFF;

            for (int k = 0; k < n_cand; ++k) {
                int32_t t_idx = ti[k];
                int32_t ebp   = eb[k];
                if (ebp <= 0) continue;

                int32_t oh = rl - ebp;
                if (oh < 0) oh = 0;

                int32_t flen = fl[k];
                double log_fl = frag_len_log_lik(flen);

                double log_strand;
                bool is_anti;
                if (has_strand) {
                    bool same = (exon_strand ==
                                 static_cast<int>(t_strand_[t_idx]));
                    log_strand = same ? log_p_sense_ : log_p_antisense_;
                    is_anti    = same ? r1_antisense_ : !r1_antisense_;
                } else {
                    log_strand = LOG_HALF;
                    is_anti = false;
                }

                double log_lik = log_strand + log_fl
                               + oh * oh_log_pen_ + log_nm;
                int32_t ct = splice_type * 2 + (is_anti ? 1 : 0);

                int32_t t_len = t_length_[t_idx];
                int32_t tx_s = 0;
                int32_t tx_e = flen > 0 ? flen : t_len;
                double cov_wt = 1.0;

                if (has_genomic && flen > 0) {
                    int32_t n_exons = exon_offsets_[t_idx + 1]
                                    - exon_offsets_[t_idx];
                    if (n_exons > 0 && t_len > 0) {
                        tx_s = genomic_to_tx_pos(genomic_start, t_idx);
                        tx_e = tx_s + flen;
                        int32_t cov_end = tx_e < t_len ? tx_e : t_len;
                        cov_wt = compute_fragment_weight(
                            tx_s, cov_end, t_len);
                    } else {
                        tx_s = 0;
                        tx_e = flen;
                    }
                }

                m_scored[m_n++] = {t_idx, oh, ct, tx_s, tx_e,
                                   log_lik, cov_wt};
                if (oh < m_min_oh) m_min_oh = oh;
            }

            // mRNA WTA winners → output vectors
            for (int j = 0; j < m_n; ++j) {
                auto& s = m_scored[j];
                if (s.oh == m_min_oh) {
                    out_ti.push_back(s.t_idx);
                    out_ll.push_back(s.log_lik);
                    out_ct.push_back(static_cast<uint8_t>(s.ct));
                    out_cw.push_back(s.cov_wt);
                    out_ts.push_back(s.tx_s);
                    out_te.push_back(s.tx_e);
                    if (s.log_lik > best_ll) {
                        best_ll = s.log_lik;
                        best_t  = s.t_idx;
                        best_ct = s.ct;
                    }
                }
            }

            // ========== nRNA scoring (skip if SPLICE_ANNOT) ==========
            if (splice_type != SPLICE_ANNOT_VAL) {
                struct NrnaScored {
                    int32_t t_idx, oh, tx_s, tx_e;
                    double  nrna_ll, cov_wt;
                };
                NrnaScored n_stack[64];
                std::vector<NrnaScored> n_heap;
                NrnaScored* n_scored = n_stack;
                if (n_cand > 64) {
                    n_heap.resize(n_cand);
                    n_scored = n_heap.data();
                }
                int n_n = 0;
                int32_t n_min_oh = 0x7FFFFFFF;

                double n_log_fl = frag_len_log_lik(genomic_footprint);
                bool n_has_genomic = (genomic_start >= 0
                                      && genomic_footprint > 0);
                int32_t nrna_gfp = genomic_footprint > 0
                                 ? genomic_footprint : 0;

                for (int k = 0; k < n_cand; ++k) {
                    int32_t t_idx = ti[k];
                    int32_t t_span   = t_span_[t_idx];
                    int32_t t_exonic = t_length_[t_idx];
                    if (t_span <= t_exonic) continue;

                    int32_t ebp = eb[k];
                    int32_t ibp = ib[k];
                    int32_t span_bp = ebp + ibp;
                    if (span_bp <= 0) continue;

                    int32_t oh = rl - span_bp;
                    if (oh < 0) oh = 0;

                    double log_strand;
                    if (has_strand) {
                        bool same = (exon_strand ==
                                     static_cast<int>(t_strand_[t_idx]));
                        log_strand = same ? log_p_sense_
                                         : log_p_antisense_;
                    } else {
                        log_strand = LOG_HALF;
                    }

                    double nrna_ll = log_strand + n_log_fl
                                   + oh * oh_log_pen_ + log_nm;

                    int32_t tx_s = 0;
                    int32_t tx_e = nrna_gfp > 0 ? nrna_gfp : t_span;
                    double cov_wt = 1.0;

                    if (n_has_genomic) {
                        int32_t tx_start_g = t_start_[t_idx];
                        int32_t frag_pos     = genomic_start - tx_start_g;
                        int32_t frag_end_pos = frag_pos + genomic_footprint;
                        int32_t fpc = frag_pos < t_span
                                    ? frag_pos : t_span;
                        if (fpc < 0) fpc = 0;
                        int32_t fec = frag_end_pos < t_span
                                    ? frag_end_pos : t_span;
                        if (fec < 0) fec = 0;
                        if (fec > fpc && t_span > 0)
                            cov_wt = compute_fragment_weight(
                                fpc, fec, t_span);
                        tx_s = frag_pos > 0 ? frag_pos : 0;
                        tx_e = tx_s + genomic_footprint;
                    }

                    n_scored[n_n++] = {t_idx, oh, tx_s, tx_e,
                                       nrna_ll, cov_wt};
                    if (oh < n_min_oh) n_min_oh = oh;
                }

                // nRNA WTA winners → output vectors (with nrna_base)
                for (int j = 0; j < n_n; ++j) {
                    auto& s = n_scored[j];
                    if (s.oh == n_min_oh) {
                        out_ti.push_back(nrna_base_ + s.t_idx);
                        out_ll.push_back(s.nrna_ll);
                        out_ct.push_back(0);  // nRNA count_col = 0
                        out_cw.push_back(s.cov_wt);
                        out_ts.push_back(s.tx_s);
                        out_te.push_back(s.tx_e);
                    }
                }
            }
        }

    pack_output:
        // Pack vectors as raw bytes for array.array.frombytes()
        auto to_bytes = [](const void* data, size_t nbytes) -> nb::bytes {
            if (nbytes == 0)
                return nb::bytes(static_cast<const char*>(nullptr), 0);
            return nb::bytes(
                reinterpret_cast<const char*>(data), nbytes);
        };

        return nb::make_tuple(
            best_ll, best_t, best_ct,
            to_bytes(out_ti.data(), out_ti.size() * sizeof(int32_t)),
            to_bytes(out_ll.data(), out_ll.size() * sizeof(double)),
            to_bytes(out_ct.data(), out_ct.size() * sizeof(uint8_t)),
            to_bytes(out_cw.data(), out_cw.size() * sizeof(double)),
            to_bytes(out_ts.data(), out_ts.size() * sizeof(int32_t)),
            to_bytes(out_te.data(), out_te.size() * sizeof(int32_t)));
    }

    // ---------------------------------------------------------------
    // scan_chunk — process entire chunk in C++ (Phase 1 optimization)
    // ---------------------------------------------------------------
    //
    // Processes all non-multimapper, non-chimeric fragments in one call.
    // Replaces per-fragment Python dispatch, BufferedFragment creation,
    // pre-EM accumulation, assign_unambig, and score_emit_fragment calls.
    //
    // Returns a tuple of numpy arrays + scalar stats.

    nb::tuple scan_chunk(
        // Chunk columnar arrays
        i64_1d  chunk_t_offsets,      // int64[N+1]
        i32_1d  chunk_t_indices,      // int32[M]
        i32_1d  chunk_frag_lengths,   // int32[M]
        i32_1d  chunk_exon_bp,        // int32[M]
        i32_1d  chunk_intron_bp,      // int32[M]
        i32_1d  chunk_unambig_intron, // int32[M]
        u8_1d   chunk_splice_type,    // uint8[N]
        u8_1d   chunk_exon_strand,    // uint8[N]
        u8_1d   chunk_frag_classes,   // uint8[N]
        i64_1d  chunk_frag_id,        // int64[N]
        u32_1d  chunk_read_length,    // uint32[N]
        i32_1d  chunk_genomic_fp,     // int32[N]
        i32_1d  chunk_genomic_start,  // int32[N]
        u16_1d  chunk_nm,             // uint16[N]
        // Index arrays for strand classification
        i8_1d   t_to_strand_arr,      // int8[T]
        // Pre-EM accumulator arrays (modified in-place)
        f64_mut exonic_sense,
        f64_mut exonic_antisense,
        f64_mut unspliced_sense,
        f64_mut unspliced_antisense,
        f64_mut intronic_sense,
        f64_mut intronic_antisense,
        f64_2d_mut unambig_counts,    // float64[T, 6]
        // gDNA splice penalty for unspliced (log already applied outside)
        double  gdna_log_splice_pen_unspliced)
    {
        static constexpr int FRAG_UNAMBIG            = 0;
        static constexpr int FRAG_AMBIG_SAME_STRAND  = 1;
        // static constexpr int FRAG_AMBIG_OPP_STRAND = 2;
        static constexpr int FRAG_MULTIMAPPER        = 3;
        static constexpr int FRAG_CHIMERIC           = 4;
        static constexpr int SPLICE_UNSPLICED_VAL    = 0;
        // static constexpr int SPLICE_UNANNOT_VAL   = 1;
        static constexpr int SPLICE_ANNOT_VAL        = 2;
        static constexpr double NEG_INF =
            -std::numeric_limits<double>::infinity();

        int N = static_cast<int>(chunk_splice_type.shape(0));

        const int64_t*  t_off  = chunk_t_offsets.data();
        const int32_t*  t_ind  = chunk_t_indices.data();
        const int32_t*  f_len  = chunk_frag_lengths.data();
        const int32_t*  e_bp   = chunk_exon_bp.data();
        const int32_t*  i_bp   = chunk_intron_bp.data();
        const int32_t*  ui_bp  = chunk_unambig_intron.data();
        const uint8_t*  s_type = chunk_splice_type.data();
        const uint8_t*  e_str  = chunk_exon_strand.data();
        const uint8_t*  fc     = chunk_frag_classes.data();
        const int64_t*  f_id   = chunk_frag_id.data();
        const uint32_t* r_len  = chunk_read_length.data();
        const int32_t*  g_fp   = chunk_genomic_fp.data();
        const int32_t*  g_sta  = chunk_genomic_start.data();
        const uint16_t* nm_arr = chunk_nm.data();
        const int8_t*   t_str  = t_to_strand_arr.data();

        double* acc_exonic_s     = exonic_sense.data();
        double* acc_exonic_a     = exonic_antisense.data();
        double* acc_unspl_s      = unspliced_sense.data();
        double* acc_unspl_a      = unspliced_antisense.data();
        double* acc_intron_s     = intronic_sense.data();
        double* acc_intron_a     = intronic_antisense.data();
        double* ua_counts        = unambig_counts.data();
        int ua_ncols = static_cast<int>(unambig_counts.shape(1));

        // CSR output vectors (will be returned as numpy arrays)
        std::vector<int64_t>  out_offsets;
        std::vector<int32_t>  out_ti, out_ts, out_te;
        std::vector<double>   out_ll, out_cw;
        std::vector<uint8_t>  out_ct;

        // Per-unit metadata
        std::vector<int32_t>  out_locus_t;
        std::vector<uint8_t>  out_locus_ct;
        std::vector<int8_t>   out_is_spliced;
        std::vector<double>   out_gdna_ll;
        std::vector<int32_t>  out_gfp;
        std::vector<int64_t>  out_fid;
        std::vector<int8_t>   out_fclass;
        std::vector<uint8_t>  out_stype;

        // Det-unambig tracking (for annotation support)
        std::vector<int32_t>  det_t_indices;
        std::vector<int64_t>  det_frag_ids;

        // Stats
        int64_t n_det_unambig  = 0;
        int64_t n_em_unambig   = 0;
        int64_t n_em_ambig_ss  = 0;
        int64_t n_em_ambig_os  = 0;
        int64_t n_gated_out    = 0;
        int64_t n_chimeric     = 0;

        // Reserve reasonable capacity (typical ~50% of N goes to EM)
        out_offsets.reserve(N / 2 + 1);
        int64_t csr_len = 0;  // tracks current CSR offset

        for (int i = 0; i < N; ++i) {
            int fclass = fc[i];

            // Skip chimeric
            if (fclass == FRAG_CHIMERIC) {
                ++n_chimeric;
                continue;
            }
            // Skip multimapper (handled in Python)
            if (fclass == FRAG_MULTIMAPPER)
                continue;

            int stype      = s_type[i];
            int exon_str   = e_str[i];
            int64_t start  = t_off[i];
            int64_t end    = t_off[i + 1];
            int n_cand     = static_cast<int>(end - start);

            // ---- Pre-EM strand accumulation ----
            // For FRAG_UNAMBIG and FRAG_AMBIG_SAME_STRAND
            if (fclass == FRAG_UNAMBIG ||
                fclass == FRAG_AMBIG_SAME_STRAND)
            {
                int32_t first_t = t_ind[start];
                int first_t_strand = static_cast<int>(t_str[first_t]);

                // is_antisense: same logic as Python estimator.is_antisense
                bool is_anti;
                if (exon_str == 1 || exon_str == 2) {
                    bool same = (exon_str == first_t_strand);
                    double log_p = same ? log_p_sense_ : log_p_antisense_;
                    is_anti = std::exp(log_p) < 0.5;
                } else {
                    is_anti = false;
                }

                bool is_unspliced = (stype == SPLICE_UNSPLICED_VAL);
                double weight = n_cand > 0 ? 1.0 / n_cand : 0.0;

                for (int64_t k = start; k < end; ++k) {
                    int32_t t_idx = t_ind[k];
                    bool has_ui = (ui_bp[k] > 0);

                    // Exonic accumulation (exclude intronic fragments)
                    if (!has_ui) {
                        if (is_anti)
                            acc_exonic_a[t_idx] += weight;
                        else
                            acc_exonic_s[t_idx] += weight;
                    }

                    // Unspliced accumulation
                    if (is_unspliced) {
                        if (is_anti)
                            acc_unspl_a[t_idx] += weight;
                        else
                            acc_unspl_s[t_idx] += weight;
                    }

                    // Intronic accumulation
                    if (has_ui) {
                        if (is_anti)
                            acc_intron_a[t_idx] += weight;
                        else
                            acc_intron_s[t_idx] += weight;
                    }
                }

                // ---- Deterministic unique: SPLICE_ANNOT + FRAG_UNAMBIG ----
                if (fclass == FRAG_UNAMBIG && stype == SPLICE_ANNOT_VAL) {
                    int32_t t_idx = t_ind[start];
                    int t_strand_val = static_cast<int>(t_str[t_idx]);
                    bool anti;
                    if (exon_str == 1 || exon_str == 2) {
                        bool same = (exon_str == t_strand_val);
                        double log_p = same ? log_p_sense_
                                           : log_p_antisense_;
                        anti = std::exp(log_p) < 0.5;
                    } else {
                        anti = false;
                    }
                    int col = stype * 2 + (anti ? 1 : 0);
                    if (col < ua_ncols)
                        ua_counts[t_idx * ua_ncols + col] += 1.0;
                    ++n_det_unambig;

                    // Track for annotation support
                    det_t_indices.push_back(t_idx);
                    det_frag_ids.push_back(f_id[i]);

                    continue;  // Det-unambig DOES NOT go to EM
                }
            }

            // ---- EM-routed scoring (non-MM, non-chimeric, non-det) ----
            // Reuses the same logic as score_emit_fragment, inlined here.

            int rl = r_len[i] > 0 ? static_cast<int>(r_len[i]) : 1;
            int nm = nm_arr[i];
            double log_nm = nm > 0 ? nm * mm_log_pen_ : 0.0;
            bool has_strand  = (exon_str == 1 || exon_str == 2);
            int genomic_start = g_sta[i];
            int genomic_footprint = g_fp[i];
            bool has_genomic = (genomic_start >= 0);

            // Result vectors for this fragment
            size_t emit_start = out_ti.size();
            double best_ll = NEG_INF;
            int32_t best_t = -1;
            int32_t best_ct = 0;

            if (n_cand > 0) {
                // ========== mRNA scoring ==========
                struct MrnaScored {
                    int32_t t_idx, oh, ct, tx_s, tx_e;
                    double  log_lik, cov_wt;
                };
                MrnaScored m_stack[64];
                std::vector<MrnaScored> m_heap;
                MrnaScored* m_scored = m_stack;
                if (n_cand > 64) {
                    m_heap.resize(n_cand);
                    m_scored = m_heap.data();
                }
                int m_n = 0;
                int32_t m_min_oh = 0x7FFFFFFF;

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
                        is_anti = same ? r1_antisense_ : !r1_antisense_;
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
                        int32_t n_exons = exon_offsets_[t_idx + 1]
                                        - exon_offsets_[t_idx];
                        if (n_exons > 0 && t_len > 0) {
                            tx_s = genomic_to_tx_pos(genomic_start, t_idx);
                            tx_e = tx_s + flen;
                            int32_t cov_end = tx_e < t_len ? tx_e : t_len;
                            cov_wt = compute_fragment_weight(
                                tx_s, cov_end, t_len);
                        } else {
                            tx_s = 0;
                            tx_e = flen;
                        }
                    }

                    m_scored[m_n++] = {t_idx, oh, ct, tx_s, tx_e,
                                       log_lik, cov_wt};
                    if (oh < m_min_oh) m_min_oh = oh;
                }

                // mRNA WTA winners → output
                for (int j = 0; j < m_n; ++j) {
                    auto& s = m_scored[j];
                    if (s.oh == m_min_oh) {
                        out_ti.push_back(s.t_idx);
                        out_ll.push_back(s.log_lik);
                        out_ct.push_back(static_cast<uint8_t>(s.ct));
                        out_cw.push_back(s.cov_wt);
                        out_ts.push_back(s.tx_s);
                        out_te.push_back(s.tx_e);
                        if (s.log_lik > best_ll) {
                            best_ll = s.log_lik;
                            best_t  = s.t_idx;
                            best_ct = s.ct;
                        }
                    }
                }

                // ========== nRNA scoring (skip if SPLICE_ANNOT) ==========
                if (stype != SPLICE_ANNOT_VAL) {
                    double n_log_fl = frag_len_log_lik(genomic_footprint);
                    bool n_has_genomic = (genomic_start >= 0
                                          && genomic_footprint > 0);
                    int32_t nrna_gfp = genomic_footprint > 0
                                     ? genomic_footprint : 0;

                    struct NrnaScored {
                        int32_t t_idx, oh, tx_s, tx_e;
                        double  nrna_ll, cov_wt;
                    };
                    NrnaScored n_stack[64];
                    std::vector<NrnaScored> n_heap;
                    NrnaScored* n_scored = n_stack;
                    if (n_cand > 64) {
                        n_heap.resize(n_cand);
                        n_scored = n_heap.data();
                    }
                    int n_n = 0;
                    int32_t n_min_oh = 0x7FFFFFFF;

                    for (int64_t k = start; k < end; ++k) {
                        int32_t t_idx = t_ind[k];
                        int32_t t_span   = t_span_[t_idx];
                        int32_t t_exonic = t_length_[t_idx];
                        if (t_span <= t_exonic) continue;

                        int32_t ebp = e_bp[k];
                        int32_t ibp = i_bp[k];
                        int32_t span_bp = ebp + ibp;
                        if (span_bp <= 0) continue;

                        int32_t oh = rl - span_bp;
                        if (oh < 0) oh = 0;

                        double log_strand;
                        if (has_strand) {
                            bool same = (exon_str ==
                                         static_cast<int>(t_strand_[t_idx]));
                            log_strand = same ? log_p_sense_
                                             : log_p_antisense_;
                        } else {
                            log_strand = LOG_HALF;
                        }

                        double nrna_ll = log_strand + n_log_fl
                                       + oh * oh_log_pen_ + log_nm;

                        int32_t tx_s = 0;
                        int32_t tx_e = nrna_gfp > 0 ? nrna_gfp : t_span;
                        double cov_wt = 1.0;

                        if (n_has_genomic) {
                            int32_t tx_start_g = t_start_[t_idx];
                            int32_t frag_pos = genomic_start - tx_start_g;
                            int32_t frag_end_pos =
                                frag_pos + genomic_footprint;
                            int32_t fpc = frag_pos < t_span
                                        ? frag_pos : t_span;
                            if (fpc < 0) fpc = 0;
                            int32_t fec = frag_end_pos < t_span
                                        ? frag_end_pos : t_span;
                            if (fec < 0) fec = 0;
                            if (fec > fpc && t_span > 0)
                                cov_wt = compute_fragment_weight(
                                    fpc, fec, t_span);
                            tx_s = frag_pos > 0 ? frag_pos : 0;
                            tx_e = tx_s + genomic_footprint;
                        }

                        n_scored[n_n++] = {t_idx, oh, tx_s, tx_e,
                                           nrna_ll, cov_wt};
                        if (oh < n_min_oh) n_min_oh = oh;
                    }

                    // nRNA WTA winners → output (with nrna_base offset)
                    for (int j = 0; j < n_n; ++j) {
                        auto& s = n_scored[j];
                        if (s.oh == n_min_oh) {
                            out_ti.push_back(nrna_base_ + s.t_idx);
                            out_ll.push_back(s.nrna_ll);
                            out_ct.push_back(0);  // nRNA count_col = 0
                            out_cw.push_back(s.cov_wt);
                            out_ts.push_back(s.tx_s);
                            out_te.push_back(s.tx_e);
                        }
                    }
                }
            }

            // ---- Finalize this EM unit ----
            size_t new_cands = out_ti.size() - emit_start;
            if (new_cands > 0) {
                csr_len += static_cast<int64_t>(new_cands);
                out_offsets.push_back(csr_len);
                out_locus_t.push_back(best_t);
                out_locus_ct.push_back(static_cast<uint8_t>(best_ct));
                out_fid.push_back(f_id[i]);
                out_fclass.push_back(static_cast<int8_t>(fclass));
                out_stype.push_back(static_cast<uint8_t>(stype));

                // is_spliced and gDNA log-likelihood
                bool is_spl = (stype == SPLICE_ANNOT_VAL || stype == 1);
                out_is_spliced.push_back(is_spl ? 1 : 0);
                out_gfp.push_back(genomic_footprint);

                if (!is_spl) {
                    // gDNA log-lik for unspliced fragments
                    double gdna_fl = frag_len_log_lik(genomic_footprint);
                    out_gdna_ll.push_back(
                        LOG_HALF + gdna_fl + gdna_log_splice_pen_unspliced);
                } else {
                    out_gdna_ll.push_back(NEG_INF);
                }

                // Stats
                if (fclass == FRAG_UNAMBIG)
                    ++n_em_unambig;
                else if (fclass == FRAG_AMBIG_SAME_STRAND)
                    ++n_em_ambig_ss;
                else
                    ++n_em_ambig_os;
            } else {
                ++n_gated_out;
            }
        }

        // ---- Pack output as numpy arrays ----
        // Helper: wrap std::vector<T> as a numpy array (capsule-based
        // ownership transfer).
        auto to_np_i32 = [](std::vector<int32_t>&& v) -> nb::object {
            size_t n = v.size();
            if (n == 0) {
                auto* empty = new std::vector<int32_t>();
                nb::capsule del(empty, [](void* p) noexcept {
                    delete static_cast<std::vector<int32_t>*>(p);
                });
                return nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(
                    empty->data(), {0}, del).cast();
            }
            auto* owner = new std::vector<int32_t>(std::move(v));
            nb::capsule del(owner, [](void* p) noexcept {
                delete static_cast<std::vector<int32_t>*>(p);
            });
            return nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(
                owner->data(), {n}, del).cast();
        };

        auto to_np_i64 = [](std::vector<int64_t>&& v) -> nb::object {
            size_t n = v.size();
            if (n == 0) {
                auto* empty = new std::vector<int64_t>();
                nb::capsule del(empty, [](void* p) noexcept {
                    delete static_cast<std::vector<int64_t>*>(p);
                });
                return nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(
                    empty->data(), {0}, del).cast();
            }
            auto* owner = new std::vector<int64_t>(std::move(v));
            nb::capsule del(owner, [](void* p) noexcept {
                delete static_cast<std::vector<int64_t>*>(p);
            });
            return nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(
                owner->data(), {n}, del).cast();
        };

        auto to_np_f64 = [](std::vector<double>&& v) -> nb::object {
            size_t n = v.size();
            if (n == 0) {
                auto* empty = new std::vector<double>();
                nb::capsule del(empty, [](void* p) noexcept {
                    delete static_cast<std::vector<double>*>(p);
                });
                return nb::ndarray<nb::numpy, double, nb::ndim<1>>(
                    empty->data(), {0}, del).cast();
            }
            auto* owner = new std::vector<double>(std::move(v));
            nb::capsule del(owner, [](void* p) noexcept {
                delete static_cast<std::vector<double>*>(p);
            });
            return nb::ndarray<nb::numpy, double, nb::ndim<1>>(
                owner->data(), {n}, del).cast();
        };

        auto to_np_u8 = [](std::vector<uint8_t>&& v) -> nb::object {
            size_t n = v.size();
            if (n == 0) {
                auto* empty = new std::vector<uint8_t>();
                nb::capsule del(empty, [](void* p) noexcept {
                    delete static_cast<std::vector<uint8_t>*>(p);
                });
                return nb::ndarray<nb::numpy, uint8_t, nb::ndim<1>>(
                    empty->data(), {0}, del).cast();
            }
            auto* owner = new std::vector<uint8_t>(std::move(v));
            nb::capsule del(owner, [](void* p) noexcept {
                delete static_cast<std::vector<uint8_t>*>(p);
            });
            return nb::ndarray<nb::numpy, uint8_t, nb::ndim<1>>(
                owner->data(), {n}, del).cast();
        };

        auto to_np_i8 = [](std::vector<int8_t>&& v) -> nb::object {
            size_t n = v.size();
            if (n == 0) {
                auto* empty = new std::vector<int8_t>();
                nb::capsule del(empty, [](void* p) noexcept {
                    delete static_cast<std::vector<int8_t>*>(p);
                });
                return nb::ndarray<nb::numpy, int8_t, nb::ndim<1>>(
                    empty->data(), {0}, del).cast();
            }
            auto* owner = new std::vector<int8_t>(std::move(v));
            nb::capsule del(owner, [](void* p) noexcept {
                delete static_cast<std::vector<int8_t>*>(p);
            });
            return nb::ndarray<nb::numpy, int8_t, nb::ndim<1>>(
                owner->data(), {n}, del).cast();
        };

        return nb::make_tuple(
            // CSR arrays
            to_np_i64(std::move(out_offsets)),
            to_np_i32(std::move(out_ti)),
            to_np_f64(std::move(out_ll)),
            to_np_u8(std::move(out_ct)),
            to_np_f64(std::move(out_cw)),
            to_np_i32(std::move(out_ts)),
            to_np_i32(std::move(out_te)),
            // Per-unit metadata
            to_np_i32(std::move(out_locus_t)),
            to_np_u8(std::move(out_locus_ct)),
            to_np_i8(std::move(out_is_spliced)),
            to_np_f64(std::move(out_gdna_ll)),
            to_np_i32(std::move(out_gfp)),
            to_np_i64(std::move(out_fid)),
            to_np_i8(std::move(out_fclass)),
            to_np_u8(std::move(out_stype)),
            // Det-unambig data (for annotations)
            to_np_i32(std::move(det_t_indices)),
            to_np_i64(std::move(det_frag_ids)),
            // Stats
            n_det_unambig,
            n_em_unambig,
            n_em_ambig_ss,
            n_em_ambig_os,
            n_gated_out,
            n_chimeric
        );
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
        const int32_t*  ui_bp;
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
    };

    // Templated inner loop — shared by count and fill passes.
    // FillMode=false: count only, update accumulators (out is unused).
    // FillMode=true:  write into pre-allocated arrays, skip accumulators.
    template<bool FillMode>
    void score_chunk_impl(
        const ChunkPtrs& cp,
        const int8_t* t_str,
        double* acc_es, double* acc_ea,
        double* acc_us, double* acc_ua,
        double* acc_is, double* acc_ia,
        double* ua_counts, int ua_ncols,
        double gdna_log_sp,
        FillState& st,
        int64_t& stat_det, int64_t& stat_em_u,
        int64_t& stat_em_as, int64_t& stat_em_ao,
        int64_t& stat_gated, int64_t& stat_chim) const
    {
        static constexpr int FRAG_UNAMBIG            = 0;
        static constexpr int FRAG_AMBIG_SAME_STRAND  = 1;
        static constexpr int FRAG_MULTIMAPPER        = 3;
        static constexpr int FRAG_CHIMERIC           = 4;
        static constexpr int SPLICE_UNSPLICED_VAL    = 0;
        static constexpr int SPLICE_ANNOT_VAL        = 2;
        static constexpr double NEG_INF =
            -std::numeric_limits<double>::infinity();

        const int N            = cp.N;
        const int64_t*  t_off  = cp.t_off;
        const int32_t*  t_ind  = cp.t_ind;
        const int32_t*  f_len  = cp.f_len;
        const int32_t*  e_bp   = cp.e_bp;
        const int32_t*  i_bp   = cp.i_bp;
        const int32_t*  ui_bp  = cp.ui_bp;
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
            if (fclass == FRAG_MULTIMAPPER) continue;

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
                    int32_t first_t = t_ind[start];
                    int first_t_strand =
                        static_cast<int>(t_str[first_t]);

                    bool is_anti;
                    if (exon_str == 1 || exon_str == 2) {
                        bool same = (exon_str == first_t_strand);
                        double log_p = same ? log_p_sense_
                                           : log_p_antisense_;
                        is_anti = std::exp(log_p) < 0.5;
                    } else {
                        is_anti = false;
                    }

                    bool is_unspliced = (stype == SPLICE_UNSPLICED_VAL);
                    double weight =
                        n_cand > 0 ? 1.0 / n_cand : 0.0;

                    for (int64_t k = start; k < end; ++k) {
                        int32_t t_idx = t_ind[k];
                        bool has_ui = (ui_bp[k] > 0);

                        if (!has_ui) {
                            if (is_anti) acc_ea[t_idx] += weight;
                            else         acc_es[t_idx] += weight;
                        }
                        if (is_unspliced) {
                            if (is_anti) acc_ua[t_idx] += weight;
                            else         acc_us[t_idx] += weight;
                        }
                        if (has_ui) {
                            if (is_anti) acc_ia[t_idx] += weight;
                            else         acc_is[t_idx] += weight;
                        }
                    }

                    // Det-unambig: SPLICE_ANNOT + FRAG_UNAMBIG
                    if (fclass == FRAG_UNAMBIG &&
                        stype == SPLICE_ANNOT_VAL) {
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
                        stype == SPLICE_ANNOT_VAL)
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
                MrnaScored m_stack[64];
                std::vector<MrnaScored> m_heap;
                MrnaScored* m_scored = m_stack;
                if (n_cand > 64) {
                    m_heap.resize(n_cand);
                    m_scored = m_heap.data();
                }
                int m_n = 0;
                int32_t m_min_oh = 0x7FFFFFFF;

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
                    if (oh < m_min_oh) m_min_oh = oh;
                }

                // mRNA WTA winners → output
                for (int j = 0; j < m_n; ++j) {
                    auto& s = m_scored[j];
                    if (s.oh == m_min_oh) {
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

                // ===== nRNA scoring (skip if SPLICE_ANNOT) =====
                if (stype != SPLICE_ANNOT_VAL) {
                    double n_log_fl =
                        frag_len_log_lik(genomic_footprint);
                    bool n_has_genomic =
                        (genomic_start >= 0
                         && genomic_footprint > 0);
                    int32_t nrna_gfp =
                        genomic_footprint > 0
                        ? genomic_footprint : 0;

                    struct NrnaScored {
                        int32_t t_idx, oh, tx_s, tx_e;
                        double  nrna_ll, cov_wt;
                    };
                    NrnaScored n_stack[64];
                    std::vector<NrnaScored> n_heap;
                    NrnaScored* n_scored = n_stack;
                    if (n_cand > 64) {
                        n_heap.resize(n_cand);
                        n_scored = n_heap.data();
                    }
                    int n_n = 0;
                    int32_t n_min_oh = 0x7FFFFFFF;

                    for (int64_t k = start; k < end; ++k) {
                        int32_t t_idx = t_ind[k];
                        int32_t t_span   = t_span_[t_idx];
                        int32_t t_exonic = t_length_[t_idx];
                        if (t_span <= t_exonic) continue;

                        int32_t ebp = e_bp[k];
                        int32_t ibp = i_bp[k];
                        int32_t span_bp = ebp + ibp;
                        if (span_bp <= 0) continue;

                        int32_t oh = rl - span_bp;
                        if (oh < 0) oh = 0;

                        double log_strand;
                        if (has_strand) {
                            bool same = (exon_str ==
                                static_cast<int>(
                                    t_strand_[t_idx]));
                            log_strand = same
                                ? log_p_sense_
                                : log_p_antisense_;
                        } else {
                            log_strand = LOG_HALF;
                        }

                        double nrna_ll =
                            log_strand + n_log_fl
                            + oh * oh_log_pen_ + log_nm;

                        int32_t tx_s = 0;
                        int32_t tx_e = nrna_gfp > 0
                            ? nrna_gfp : t_span;
                        double cov_wt = 1.0;

                        if (n_has_genomic) {
                            int32_t tx_start_g =
                                t_start_[t_idx];
                            int32_t frag_pos =
                                genomic_start - tx_start_g;
                            int32_t frag_end_pos =
                                frag_pos + genomic_footprint;
                            int32_t fpc =
                                frag_pos < t_span
                                ? frag_pos : t_span;
                            if (fpc < 0) fpc = 0;
                            int32_t fec =
                                frag_end_pos < t_span
                                ? frag_end_pos : t_span;
                            if (fec < 0) fec = 0;
                            if (fec > fpc && t_span > 0)
                                cov_wt =
                                    compute_fragment_weight(
                                        fpc, fec, t_span);
                            tx_s = frag_pos > 0
                                ? frag_pos : 0;
                            tx_e = tx_s + genomic_footprint;
                        }

                        n_scored[n_n++] = {t_idx, oh, tx_s,
                                           tx_e, nrna_ll,
                                           cov_wt};
                        if (oh < n_min_oh) n_min_oh = oh;
                    }

                    // nRNA WTA winners → output
                    for (int j = 0; j < n_n; ++j) {
                        auto& s = n_scored[j];
                        if (s.oh == n_min_oh) {
                            if constexpr (FillMode) {
                                int64_t c = st.cand_cur;
                                st.ti[c] =
                                    nrna_base_ + s.t_idx;
                                st.ll[c] = s.nrna_ll;
                                st.ct[c] = 0;
                                st.cw[c] = s.cov_wt;
                                st.ts[c] = s.tx_s;
                                st.te[c] = s.tx_e;
                            }
                            ++st.cand_cur;
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
                        (stype == SPLICE_ANNOT_VAL
                         || stype == 1);
                    st.is_spliced[st.unit_cur] =
                        is_spl ? 1 : 0;
                    st.gfp[st.unit_cur] = genomic_footprint;

                    if (!is_spl) {
                        double gdna_fl =
                            frag_len_log_lik(
                                genomic_footprint);
                        st.gdna_ll[st.unit_cur] =
                                LOG_HALF + gdna_fl
                                + gdna_log_sp;
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
        f64_mut exonic_sense,
        f64_mut exonic_antisense,
        f64_mut unspliced_sense,
        f64_mut unspliced_antisense,
        f64_mut intronic_sense,
        f64_mut intronic_antisense,
        f64_2d_mut unambig_counts,
        double  gdna_log_splice_pen_unspliced)
    {
        const int8_t* t_str = t_to_strand_arr.data();
        double* acc_es = exonic_sense.data();
        double* acc_ea = exonic_antisense.data();
        double* acc_us = unspliced_sense.data();
        double* acc_ua = unspliced_antisense.data();
        double* acc_is = intronic_sense.data();
        double* acc_ia = intronic_antisense.data();
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
            cp.ui_bp  = nb::cast<i32_1d>(c[5]).data();
            cp.s_type = nb::cast<u8_1d>(c[6]).data();
            cp.e_str  = nb::cast<u8_1d>(c[7]).data();
            cp.fc     = nb::cast<u8_1d>(c[8]).data();
            cp.f_id   = nb::cast<i64_1d>(c[9]).data();
            cp.r_len  = nb::cast<u32_1d>(c[10]).data();
            cp.g_fp   = nb::cast<i32_1d>(c[11]).data();
            cp.g_sta  = nb::cast<i32_1d>(c[12]).data();
            cp.nm     = nb::cast<u16_1d>(c[13]).data();
            cp.N      = static_cast<int>(
                nb::cast<u8_1d>(c[6]).shape(0));
        }

        // ============ PASS 1: COUNT ============
        FillState count_st{};
        int64_t stat_det = 0, stat_em_u = 0;
        int64_t stat_em_as = 0, stat_em_ao = 0;
        int64_t stat_gated = 0, stat_chim = 0;

        for (size_t ci = 0; ci < n_chunks; ++ci) {
            score_chunk_impl<false>(
                cps[ci], t_str,
                acc_es, acc_ea, acc_us, acc_ua,
                acc_is, acc_ia, ua_c, ua_ncols,
                gdna_log_splice_pen_unspliced,
                count_st,
                stat_det, stat_em_u, stat_em_as,
                stat_em_ao, stat_gated, stat_chim);
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
        int64_t f_det = 0, f_em_u = 0, f_em_as = 0;
        int64_t f_em_ao = 0, f_gated = 0, f_chim = 0;

        for (size_t ci = 0; ci < n_chunks; ++ci) {
            score_chunk_impl<true>(
                cps[ci], t_str,
                acc_es, acc_ea, acc_us, acc_ua,
                acc_is, acc_ia, ua_c, ua_ncols,
                gdna_log_splice_pen_unspliced,
                fill_st,
                f_det, f_em_u, f_em_as,
                f_em_ao, f_gated, f_chim);
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
            stat_chim
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
                 i8_1d, i32_1d, i32_1d, i32_1d,
                 int32_t, nb::object>(),
             nb::arg("log_p_sense"),
             nb::arg("log_p_antisense"),
             nb::arg("r1_antisense"),
             nb::arg("overhang_log_penalty"),
             nb::arg("mismatch_log_penalty"),
             nb::arg("fl_log_prob").none(),
             nb::arg("fl_max_size"),
             nb::arg("fl_tail_base"),
             nb::arg("t_strand_arr"),
             nb::arg("t_length_arr"),
             nb::arg("t_span_arr"),
             nb::arg("t_start_arr"),
             nb::arg("nrna_base"),
             nb::arg("t_exon_data").none())
        .def("score_wta_mrna",
             &NativeFragmentScorer::score_wta_mrna,
             nb::arg("t_inds"), nb::arg("exon_bp"),
             nb::arg("frag_lengths"),
             nb::arg("exon_strand"), nb::arg("splice_type"),
             nb::arg("nm"), nb::arg("read_length"),
             nb::arg("genomic_start"))
        .def("score_wta_nrna",
             &NativeFragmentScorer::score_wta_nrna,
             nb::arg("t_inds"), nb::arg("exon_bp"),
             nb::arg("intron_bp"),
             nb::arg("exon_strand"), nb::arg("nm"),
             nb::arg("read_length"), nb::arg("genomic_footprint"),
             nb::arg("genomic_start"))
        .def("score_emit_fragment",
             &NativeFragmentScorer::score_emit_fragment,
             nb::arg("t_inds"), nb::arg("exon_bp"),
             nb::arg("frag_lengths"), nb::arg("intron_bp"),
             nb::arg("exon_strand"), nb::arg("splice_type"),
             nb::arg("nm"), nb::arg("read_length"),
             nb::arg("genomic_footprint"), nb::arg("genomic_start"))
        .def("scan_chunk",
             &NativeFragmentScorer::scan_chunk,
             nb::arg("chunk_t_offsets"),
             nb::arg("chunk_t_indices"),
             nb::arg("chunk_frag_lengths"),
             nb::arg("chunk_exon_bp"),
             nb::arg("chunk_intron_bp"),
             nb::arg("chunk_unambig_intron"),
             nb::arg("chunk_splice_type"),
             nb::arg("chunk_exon_strand"),
             nb::arg("chunk_frag_classes"),
             nb::arg("chunk_frag_id"),
             nb::arg("chunk_read_length"),
             nb::arg("chunk_genomic_fp"),
             nb::arg("chunk_genomic_start"),
             nb::arg("chunk_nm"),
             nb::arg("t_to_strand_arr"),
             nb::arg("exonic_sense"),
             nb::arg("exonic_antisense"),
             nb::arg("unspliced_sense"),
             nb::arg("unspliced_antisense"),
             nb::arg("intronic_sense"),
             nb::arg("intronic_antisense"),
             nb::arg("unambig_counts"),
             nb::arg("gdna_log_splice_pen_unspliced"))
        .def("fused_score_buffer",
             &NativeFragmentScorer::fused_score_buffer,
             nb::arg("chunk_arrays"),
             nb::arg("t_to_strand_arr"),
             nb::arg("exonic_sense"),
             nb::arg("exonic_antisense"),
             nb::arg("unspliced_sense"),
             nb::arg("unspliced_antisense"),
             nb::arg("intronic_sense"),
             nb::arg("intronic_antisense"),
             nb::arg("unambig_counts"),
             nb::arg("gdna_log_splice_pen_unspliced"));

    // Export scoring constants for Python-side parity tests
    m.attr("LOG_HALF")      = rigel::LOG_HALF;
    m.attr("TAIL_DECAY_LP") = rigel::TAIL_DECAY_LP;
}
