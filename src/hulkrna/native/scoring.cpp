/**
 * scoring.cpp — nanobind C++ extension for hot-path scoring functions.
 *
 * This module provides performance-critical kernels called from the
 * Python scoring loop.  Functions are exposed as ``hulkrna._scoring_impl``.
 *
 * Contents:
 *   - compute_fragment_weight(frag_start, frag_end, transcript_length)
 *   - NativeFragmentScorer class:
 *       .score_wta_mrna(...)         — mRNA WTA scoring kernel
 *       .score_wta_nrna(...)         — nRNA WTA scoring kernel
 *       .score_emit_fragment(...)    — fused score+emit (bytes output)
 *
 * Build:
 *   Part of the hulkrna scikit-build-core build — see CMakeLists.txt.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

// ----------------------------------------------------------------
// Constants
// ----------------------------------------------------------------

static constexpr double LOG_HALF = -0.6931471805599453;        // log(0.5)
static constexpr double TAIL_DECAY_LP = -0.01005033585350145;  // log(0.99)
static constexpr int STRAND_NEG = 2;

// ----------------------------------------------------------------
// Array type aliases
// ----------------------------------------------------------------

using i32_1d = nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig>;
using i8_1d  = nb::ndarray<const int8_t,  nb::ndim<1>, nb::c_contig>;
using f64_1d = nb::ndarray<const double,  nb::ndim<1>, nb::c_contig>;

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
    bool   anti_flag_;
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
        bool   anti_flag,
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
        anti_flag_(anti_flag),
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
                is_anti    = same ? anti_flag_ : !anti_flag_;
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
                    is_anti    = same ? anti_flag_ : !anti_flag_;
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
};

// ----------------------------------------------------------------
// nanobind module definition
// ----------------------------------------------------------------

NB_MODULE(_scoring_impl, m) {
    m.doc() = "C++ hot-path scoring kernels for hulkrna (nanobind)";

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
             nb::arg("anti_flag"),
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
             nb::arg("genomic_footprint"), nb::arg("genomic_start"));
}
