/**
 * scoring.cpp — nanobind C++ extension for hot-path scoring functions.
 *
 * This module provides performance-critical kernels called from the
 * Python scoring loop.  Functions are exposed as ``rigel._scoring_impl``.
 *
 * Contents:
 *   - compute_fragment_weight (the coverage weight, used by the kernels below)
 *   - NativeFragmentScorer class (scoring parameters + per-chunk kernels)
 *   - StreamingScorer class (stateful per-chunk scoring for streaming memory)
 *
 * Build:
 *   Part of the rigel scikit-build-core build — see CMakeLists.txt.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "constants.h"
#include "ndarray_util.h"

namespace nb = nanobind;

// Use constants from constants.h via the rigel namespace.
using rigel::LOG_HALF;
using rigel::TAIL_DECAY_LP;
using rigel::STRAND_NEG;
using rigel::SPLICE_UNSPLICED;
using rigel::SPLICE_SPLICED_UNANNOT;
using rigel::SPLICE_SPLICED_ANNOT;
using rigel::SPLICE_IMPLICIT;
using rigel::SPLICE_ARTIFACT;
using rigel::vec_to_ndarray;
using rigel::FRAG_UNAMBIG;
using rigel::FRAG_AMBIG_SAME_STRAND;
using rigel::FRAG_AMBIG_OPP_STRAND;
using rigel::FRAG_MULTIMAPPER;
using rigel::FRAG_CHIMERIC;

// Maximum candidates handled via stack-allocated buffer before falling
// back to heap allocation.  Avoids malloc in the common case where a
// fragment overlaps a moderate number of transcripts.
static constexpr int SCORED_STACK_CAPACITY = 64;

// ----------------------------------------------------------------
// Online log-sum-exp accumulator (numerically stable)
// ----------------------------------------------------------------
//
// Maintains running (max_v, sum_v) such that the logsumexp of the
// stream of inputs x_1,..,x_n equals  max_v + log(sum_v).
// Uses an explicit has_v flag instead of an infinity sentinel because
// this extension is compiled with fast-math.
static inline void lse_update(
    double& max_v,
    double& sum_v,
    bool& has_v,
    double x)
{
    if (!has_v) {
        max_v = x;
        sum_v = 1.0;
        has_v = true;
    } else if (x > max_v) {
        sum_v = sum_v * std::exp(max_v - x) + 1.0;
        max_v = x;
    } else {
        sum_v += std::exp(x - max_v);
    }
}

// ----------------------------------------------------------------
// Array type aliases
// ----------------------------------------------------------------

using i32_1d  = nb::ndarray<const int32_t,  nb::ndim<1>, nb::c_contig>;
using i8_1d   = nb::ndarray<const int8_t,  nb::ndim<1>, nb::c_contig>;
using f64_1d  = nb::ndarray<const double,  nb::ndim<1>, nb::c_contig>;
using i64_1d  = nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig>;
using u8_1d   = nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig>;
using u16_1d  = nb::ndarray<const uint16_t, nb::ndim<1>, nb::c_contig>;
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
// a transcript boundary gets weight > 1.0, reflecting that it is less
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

    /// Can gDNA have produced this alignment? The ONE rule, for a single fragment and for every hit of a
    /// multimapper alike, before any likelihood is read.
    ///   UNSPLICED  yes.
    ///   ARTIFACT   yes: every junction it carried was rejected by the blacklist, so it is unspliced.
    ///   IMPLICIT   yes: its mate gap is unobserved, so gDNA reads it unbroken at the footprint and a
    ///              transcript at its own length, and the fragment length tells them apart.
    ///   SPLICED_ANNOT / SPLICED_UNANNOT  no: a sequenced junction that survived the blacklist.
    /// Whether gDNA then COMPETES is the likelihood's: its term is the unspliced term at the genomic
    /// footprint whatever the label, and it is pruned exactly as an RNA candidate is (`gdna_competes`).
    inline bool gdna_can_explain(int stype) const {
        return stype == SPLICE_UNSPLICED || stype == SPLICE_ARTIFACT || stype == SPLICE_IMPLICIT;
    }

    /// The gDNA candidate is kept when its likelihood is within the pruning bound of the unit's best RNA
    /// candidate — `max_ll_delta_ = -log(pruning_min_posterior)`, the bound every RNA candidate meets within
    /// its pool — and is otherwise no candidate at all. A footprint the gDNA length law cannot produce is
    /// pruned by this, with no length rule of its own.
    /// ⛔ The callers carry the outcome as a bool, never as an infinity tested later: this target is built
    /// with -ffast-math (CMakeLists.txt), under which `std::isfinite(-inf)` may return true.
    inline bool gdna_competes(double gdna_ll, double best_rna_ll) const {
        return best_rna_ll - gdna_ll <= max_ll_delta_;
    }

    /// The fragment's coverage weight on transcript ``t_idx`` (``compute_fragment_weight``), read on the interval it
    /// occupies in the transcript's OWN 5'->3' coordinates, cut to the transcript. Its start is measured as the
    /// resolver measures the fragment's length (``tx_frag_length``): inside an exon at the exon's offset; in an
    /// intron, or before the first exon, the overhang's length before the next exon's start; past the last exon
    /// that far past the end. So it occupies ``[fwd, fwd + flen)`` forward, and on a minus-strand transcript, whose
    /// 5'->3' runs against the genome, the mirror ``[t_len - fwd - flen, t_len - fwd)``. A fragment and its mirror
    /// image on the other strand weigh the same. Held to an independent evaluation of the coverage model on both
    /// strands in ``tests/test_pipeline_routing.py``.
    inline double coverage_weight(int32_t genomic_start, int32_t flen, int32_t t_idx) const {
        const int32_t begin   = exon_offsets_[t_idx];
        const int32_t n_exons = exon_offsets_[t_idx + 1] - begin;
        const int32_t t_len   = t_length_[t_idx];
        if (n_exons <= 0 || t_len <= 0) return 1.0;

        const int32_t* starts = exon_starts_.data() + begin;
        const int32_t* ends   = exon_ends_.data()   + begin;
        const int32_t* cumsum = exon_cumsum_.data()  + begin;

        // the last exon starting at or before the fragment's start (-1: before the first exon)
        const int ei = static_cast<int>(
            std::upper_bound(starts, starts + n_exons, genomic_start) - starts) - 1;
        int32_t fwd;
        if (ei >= 0 && genomic_start < ends[ei]) {
            fwd = cumsum[ei] + (genomic_start - starts[ei]);
        } else if (ei + 1 < n_exons) {
            fwd = cumsum[ei + 1] - (starts[ei + 1] - genomic_start);
        } else {
            fwd = t_len + (genomic_start - ends[ei]);
        }

        int32_t lo, hi;
        if (t_strand_[t_idx] == STRAND_NEG) {
            hi = t_len - fwd;
            lo = hi - flen;
        } else {
            lo = fwd;
            hi = fwd + flen;
        }
        return compute_fragment_weight(std::max<int32_t>(lo, 0), std::min(hi, t_len), t_len);
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
        i32_1d  exon_offsets_arr,
        i32_1d  exon_starts_arr,
        i32_1d  exon_ends_arr,
        i32_1d  exon_cumsum_arr,
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

        // Copy pre-built exon CSR arrays directly
        {
            const auto* p = exon_offsets_arr.data();
            int32_t n = static_cast<int32_t>(exon_offsets_arr.shape(0));
            exon_offsets_.assign(p, p + n);
        }
        {
            const auto* p = exon_starts_arr.data();
            int32_t n = static_cast<int32_t>(exon_starts_arr.shape(0));
            exon_starts_.assign(p, p + n);
        }
        {
            const auto* p = exon_ends_arr.data();
            int32_t n = static_cast<int32_t>(exon_ends_arr.shape(0));
            exon_ends_.assign(p, p + n);
        }
        {
            const auto* p = exon_cumsum_arr.data();
            int32_t n = static_cast<int32_t>(exon_cumsum_arr.shape(0));
            exon_cumsum_.assign(p, p + n);
        }
    }

    friend class StreamingScorer;

private:

    // Cached raw pointers for one prepared chunk.
    struct ChunkPtrs {
        const int32_t*  t_off;
        const int32_t*  t_ind;
        const int32_t*  f_len;
        const uint16_t* e_bp;
        const uint8_t*  s_type;
        const uint8_t*  e_str;
        const uint8_t*  fc;
        const int64_t*  f_id;
        const uint16_t* r_len;
        const int32_t*  g_fp;
        const int32_t*  g_sta;
        const uint16_t* nm;
        int N;
    };

    // Merged candidate entry for MM cross-alignment merge.
    struct MergedMrna {
        int32_t overhang;
        double  log_lik;
        int32_t count_col;
        double  coverage_wt;
    };

    // Growing output vectors + write cursors for single-pass scoring.
    struct FillState {
        // CSR candidate arrays (growing)
        std::vector<int32_t>*  v_ti;
        std::vector<float>*    v_ll;
        std::vector<uint8_t>*  v_ct;
        std::vector<float>*    v_cw;
        // CSR offsets (growing; starts with {0})
        std::vector<int64_t>*  v_offsets;
        // Per-unit metadata (growing)
        std::vector<int32_t>*  v_locus_t;
        std::vector<uint8_t>*  v_locus_ct;
        std::vector<int8_t>*   v_is_spliced;
        std::vector<float>*    v_gdna_ll;
        std::vector<int64_t>*  v_gmid;  // genomic midpoint per unit; INT64_MIN sentinel
        std::vector<int64_t>*  v_fid;
        std::vector<int8_t>*   v_fclass;
        std::vector<uint8_t>*  v_stype;
        // Det-unambig (growing)
        std::vector<int32_t>*  v_det_ti;
        std::vector<int64_t>*  v_det_fid;
        // Chimeric accumulation (growing)
        std::vector<int64_t>*  v_chim_fid;
        std::vector<uint8_t>*  v_chim_stype;
        // Unambig counts accumulation (single-pass)
        double* ua_counts;
        int ua_ncols;
        const int8_t* t_str;
        // Cursors
        int64_t cand_cur;
        int64_t unit_cur;
        int64_t det_cur;

        // Eager MM group accumulators (persists across chunks).
        // Candidates are scored and merged eagerly as alignments
        // arrive.  flush_mm_group() only prunes + emits.
        int64_t mm_fid;
        int     mm_n_members;
        std::unordered_map<int32_t, MergedMrna> mm_merged;
        int     mm_best_stype;

        // Online logsumexp of per-hit gDNA log-liks. The emitted
        // gdna_log_lik at flush time is
        //   lse_max + log(lse_sum) - log(nh_gdna).
        double  mm_gdna_lse_max;
        double  mm_gdna_lse_sum;
        bool    mm_gdna_lse_has;
        int32_t mm_nh_gdna;
        // First valid unspliced (genomic-midpoint, gDNA-bearing) hit
        // observed in the current MM group. INT64_MIN means none.
        int64_t mm_gmid_first;

        FillState()
            : v_ti(nullptr), v_ll(nullptr), v_ct(nullptr),
                            v_cw(nullptr),
              v_offsets(nullptr), v_locus_t(nullptr),
              v_locus_ct(nullptr), v_is_spliced(nullptr),
              v_gdna_ll(nullptr), v_gmid(nullptr), v_fid(nullptr),
              v_fclass(nullptr), v_stype(nullptr),
              v_det_ti(nullptr), v_det_fid(nullptr),
              v_chim_fid(nullptr), v_chim_stype(nullptr),
              ua_counts(nullptr), ua_ncols(0), t_str(nullptr),
              cand_cur(0), unit_cur(0), det_cur(0),
              mm_fid(-1), mm_n_members(0),
              mm_best_stype(SPLICE_UNSPLICED),
              mm_gdna_lse_max(0.0),
              mm_gdna_lse_sum(0.0),
              mm_gdna_lse_has(false),
              mm_nh_gdna(0),
              mm_gmid_first(std::numeric_limits<int64_t>::min()) {}

        void reset_mm_group() {
            mm_merged.clear();  // reuses hash bucket allocation
            mm_fid = -1;
            mm_n_members = 0;
            mm_best_stype = SPLICE_UNSPLICED;
            mm_gdna_lse_max = 0.0;
            mm_gdna_lse_sum = 0.0;
            mm_gdna_lse_has = false;
            mm_nh_gdna = 0;
            mm_gmid_first = std::numeric_limits<int64_t>::min();
        }
    };

    // ---------------------------------------------------------------
    // score_mm_alignment — eagerly score one MM alignment
    // ---------------------------------------------------------------
    //
    // Called during score_chunk_impl() for each MM alignment.
    // Reads from the current chunk's live arrays and merges scored
    // candidates into st.mm_merged.  Also updates the per-group
    // metadata accumulators (splice status, gDNA, footprint).

    void score_mm_alignment(
        const ChunkPtrs& cp, int row,
        FillState& st, double gdna_log_sp) const
    {
        int64_t start = cp.t_off[row];
        int64_t end   = cp.t_off[row + 1];
        int n_cand    = static_cast<int>(end - start);

        int stype    = cp.s_type[row];
        int exon_str = cp.e_str[row];
        int rl = cp.r_len[row] > 0
               ? static_cast<int>(cp.r_len[row]) : 1;
        int nm = cp.nm[row];
        double log_nm = nm > 0 ? nm * mm_log_pen_ : 0.0;
        bool has_strand  = (exon_str == 1 || exon_str == 2);
        int genomic_start_val = cp.g_sta[row];
        bool has_genomic = (genomic_start_val >= 0);

        // Update per-group metadata accumulators (the reported splice type only)
        if (stype == SPLICE_SPLICED_ANNOT) {
            st.mm_best_stype = SPLICE_SPLICED_ANNOT;
        } else if (stype == SPLICE_SPLICED_UNANNOT &&
                   st.mm_best_stype != SPLICE_SPLICED_ANNOT) {
            st.mm_best_stype = SPLICE_SPLICED_UNANNOT;
        }

        // gDNA per-hit contribution, from every hit gDNA can explain (`gdna_can_explain`). A spliced hit
        // no longer removes gDNA from the group: its other hits are places gDNA could have come from.
        // Effective-length normalization is component-level in the EM; the scorer emits log h_G(ell_f)
        // plus non-length score terms. At flush time:
        //    gdna_log_lik = lse_max + log(lse_sum) - log(nh_gdna)
        // the mean over the hits gDNA can explain (nh_gdna of them), not over NH.
        int32_t gfp_val = cp.g_fp[row];
        if (gdna_can_explain(stype) && n_cand > 0)
        {
            double gdna_fl    = gdna_frag_len_log_lik(gfp_val);
            double hit_log_ll = gdna_fl + gdna_log_sp + LOG_HALF + log_nm;
            lse_update(st.mm_gdna_lse_max,
                       st.mm_gdna_lse_sum,
                       st.mm_gdna_lse_has,
                       hit_log_ll);
            ++st.mm_nh_gdna;

            // Track the first valid (has_genomic + has gDNA hypothesis)
            // midpoint observed in this MM group for regional weighting.
            if (has_genomic &&
                st.mm_gmid_first ==
                    std::numeric_limits<int64_t>::min()) {
                st.mm_gmid_first =
                    static_cast<int64_t>(genomic_start_val)
                    + static_cast<int64_t>(gfp_val) / 2;
            }
        }

        // Track member count for NH bookkeeping
        ++st.mm_n_members;

        // Score and merge mRNA candidates
        if (n_cand <= 0) return;

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

            double cov_wt = (has_genomic && flen > 0)
                ? coverage_weight(genomic_start_val, flen, t_idx) : 1.0;

            auto it = st.mm_merged.find(t_idx);
            if (it == st.mm_merged.end()) {
                st.mm_merged[t_idx] = {oh, log_lik, ct,
                                       cov_wt};
            } else {
                auto& prev = it->second;
                if (oh < prev.overhang ||
                    (oh == prev.overhang &&
                     log_lik > prev.log_lik))
                {
                    prev = {oh, log_lik, ct,
                            cov_wt};
                }
            }
        }
    }

    // ---------------------------------------------------------------
    // flush_mm_group — prune and emit one multimapper group
    // ---------------------------------------------------------------
    //
    // All scoring has already been done eagerly by score_mm_alignment.
    // This method only prunes mm_merged and emits to the CSR vectors,
    // then resets the accumulators for the next group.

    void flush_mm_group(
        FillState& st,
        double gdna_log_sp,
        int64_t& stat_mm,
        int64_t& stat_gated) const
    {
        static constexpr double NEG_INF =
            -std::numeric_limits<double>::infinity();

        if (st.mm_n_members == 0) return;

        // Pool-separated likelihood pruning
        double mrna_best = NEG_INF;
        double nrna_best = NEG_INF;
        for (auto& [t_idx, mc] : st.mm_merged) {
            if (t_is_nrna_[t_idx])
                nrna_best = std::max(nrna_best, mc.log_lik);
            else
                mrna_best = std::max(mrna_best, mc.log_lik);
        }

        int64_t emit_start = st.cand_cur;
        double best_ll = NEG_INF;
        int32_t best_t = -1;
        int32_t best_ct = 0;

        for (auto& [t_idx, mc] : st.mm_merged) {
            double pool_best = t_is_nrna_[t_idx]
                             ? nrna_best : mrna_best;
            if (pool_best - mc.log_lik <= max_ll_delta_) {
                st.v_ti->push_back(t_idx);
                st.v_ll->push_back(static_cast<float>(mc.log_lik));
                st.v_ct->push_back(
                    static_cast<uint8_t>(mc.count_col));
                st.v_cw->push_back(static_cast<float>(mc.coverage_wt));
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
            st.v_offsets->push_back(st.cand_cur);
            st.v_locus_t->push_back(best_t);
            st.v_locus_ct->push_back(
                static_cast<uint8_t>(best_ct));
            st.v_fid->push_back(st.mm_fid);
            st.v_fclass->push_back(
                static_cast<int8_t>(FRAG_MULTIMAPPER));

            st.v_stype->push_back(
                static_cast<uint8_t>(st.mm_best_stype));

            // One per-unit gDNA log-lik over the group's gDNA-explainable hits, kept when it competes
            // with the group's best RNA candidate. The EM applies the per-locus gDNA effective length.
            double final_gdna_ll = 0.0;
            bool has_gdna = false;
            if (st.mm_nh_gdna > 0 &&
                st.mm_gdna_lse_has &&
                st.mm_gdna_lse_sum > 0.0) {
                final_gdna_ll = st.mm_gdna_lse_max
                              + std::log(st.mm_gdna_lse_sum)
                              - std::log((double)st.mm_nh_gdna);
                has_gdna = gdna_competes(final_gdna_ll, best_ll);
            }
            // "spliced" = gDNA is no candidate: the group is certified RNA.
            st.v_is_spliced->push_back(has_gdna ? 0 : 1);
            if (has_gdna) {
                st.v_gdna_ll->push_back(static_cast<float>(final_gdna_ll));
                st.v_gmid->push_back(st.mm_gmid_first);
            } else {
                st.v_gdna_ll->push_back(static_cast<float>(NEG_INF));
                st.v_gmid->push_back(
                    std::numeric_limits<int64_t>::min());
            }

            ++st.unit_cur;
            ++stat_mm;
        } else {
            ++stat_gated;
        }

        st.reset_mm_group();
    }

    // Single-pass inner loop — scores fragments and pushes results
    // to growing output vectors.  Also accumulates ua_counts for
    // det-unambig fragments (merged from the former count pass).
    void score_chunk_impl(
        const ChunkPtrs& cp,
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
        const int32_t*  t_off  = cp.t_off;
        const int32_t*  t_ind  = cp.t_ind;
        const int32_t*  f_len  = cp.f_len;
        const uint16_t* e_bp   = cp.e_bp;
        const uint8_t*  s_type = cp.s_type;
        const uint8_t*  e_str  = cp.e_str;
        const uint8_t*  fc     = cp.fc;
        const int64_t*  f_id   = cp.f_id;
        const uint16_t* r_len  = cp.r_len;
        const int32_t*  g_fp   = cp.g_fp;
        const int32_t*  g_sta  = cp.g_sta;
        const uint16_t* nm_arr = cp.nm;

        for (int i = 0; i < N; ++i) {
            int fclass = fc[i];

            if (fclass == FRAG_CHIMERIC) {
                // Accumulate chimeric info for annotation
                st.v_chim_fid->push_back(f_id[i]);
                st.v_chim_stype->push_back(s_type[i]);
                ++stat_chim;
                continue;
            }

            // ---- Multimapper: eager scoring ----
            if (fclass == FRAG_MULTIMAPPER) {
                int64_t fid_val = f_id[i];
                if (fid_val != st.mm_fid) {
                    if (st.mm_n_members > 0) {
                        flush_mm_group(
                            st, gdna_log_sp,
                            stat_mm, stat_gated);
                    }
                    st.mm_fid = fid_val;
                }
                score_mm_alignment(cp, i, st, gdna_log_sp);
                continue;
            }

            // Non-MM fragment: flush any pending MM group
            if (st.mm_n_members > 0) {
                flush_mm_group(
                    st, gdna_log_sp,
                    stat_mm, stat_gated);
            }

            int stype      = s_type[i];
            int exon_str   = e_str[i];
            int64_t start  = t_off[i];
            int64_t end    = t_off[i + 1];
            int n_cand     = static_cast<int>(end - start);

            // zero-candidate fragments are present in the
            // buffer for calibration only (intergenic / dropped
            // by merge).  Skip silently — they must not contribute
            // to stat_gated, otherwise accountability double-counts
            // (the legacy n_intergenic_unspliced/_spliced counter
            // already represents these fragments).
            if (n_cand <= 0) continue;

            // ---- Det-unambig: both ua_counts + det arrays ----
            if (fclass == FRAG_UNAMBIG ||
                fclass == FRAG_AMBIG_SAME_STRAND)
            {
                if (fclass == FRAG_UNAMBIG &&
                    stype == SPLICE_SPLICED_ANNOT) {
                    int32_t t_idx = t_ind[start];

                    // ua_counts accumulation
                    {
                        int t_sv =
                            static_cast<int>(st.t_str[t_idx]);
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
                        if (col < st.ua_ncols)
                            st.ua_counts[
                                t_idx * st.ua_ncols + col]
                                += 1.0;
                    }

                    // Det-unambig output arrays
                    st.v_det_ti->push_back(t_idx);
                    st.v_det_fid->push_back(f_id[i]);
                    ++st.det_cur;
                    ++stat_det;
                    continue;
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
                    int32_t t_idx, oh, ct;
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

                    double cov_wt = (has_genomic && flen > 0)
                        ? coverage_weight(genomic_start, flen, t_idx) : 1.0;

                    m_scored[m_n++] = {t_idx, oh, ct, log_lik, cov_wt};
                }

                // Pool-separated likelihood pruning
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

                for (int j = 0; j < m_n; ++j) {
                    auto& s = m_scored[j];
                    double pool_best = t_is_nrna_[s.t_idx]
                                     ? nrna_best : mrna_best;
                    if (pool_best - s.log_lik <= max_ll_delta_) {
                        st.v_ti->push_back(s.t_idx);
                        st.v_ll->push_back(static_cast<float>(s.log_lik));
                        st.v_ct->push_back(
                            static_cast<uint8_t>(s.ct));
                        st.v_cw->push_back(static_cast<float>(s.cov_wt));
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
                st.v_offsets->push_back(st.cand_cur);
                st.v_locus_t->push_back(best_t);
                st.v_locus_ct->push_back(
                    static_cast<uint8_t>(best_ct));
                st.v_fid->push_back(f_id[i]);
                st.v_fclass->push_back(
                    static_cast<int8_t>(fclass));
                st.v_stype->push_back(
                    static_cast<uint8_t>(stype));

                // Non-MM path (NH=1 or ambig-same-strand): single-term gDNA
                // expression, kept when it competes with the unit's best RNA
                // candidate. EM applies per-locus gDNA effective length.
                double gdna_ll = 0.0;
                bool has_gdna = false;
                if (gdna_can_explain(stype) && best_t >= 0) {
                    gdna_ll = gdna_frag_len_log_lik(genomic_footprint) + gdna_log_sp + LOG_HALF + log_nm;
                    has_gdna = gdna_competes(gdna_ll, best_ll);
                }
                // "spliced" = gDNA is no candidate: certified RNA.
                st.v_is_spliced->push_back(
                    has_gdna ? 0 : 1);

                if (has_gdna) {
                    st.v_gdna_ll->push_back(static_cast<float>(gdna_ll));
                    if (has_genomic) {
                        st.v_gmid->push_back(
                            static_cast<int64_t>(genomic_start)
                            + static_cast<int64_t>(genomic_footprint) / 2);
                    } else {
                        st.v_gmid->push_back(
                            std::numeric_limits<int64_t>::min());
                    }
                } else {
                    st.v_gdna_ll->push_back(static_cast<float>(NEG_INF));
                    st.v_gmid->push_back(
                        std::numeric_limits<int64_t>::min());
                }

                ++st.unit_cur;

                if (fclass == FRAG_UNAMBIG)
                    ++stat_em_u;
                else if (fclass == FRAG_AMBIG_SAME_STRAND)
                    ++stat_em_as;
                else
                    ++stat_em_ao;
            } else {
                ++stat_gated;
            }
        }
    }

    // (fused_score_buffer removed — use StreamingScorer instead)
};

// ================================================================
// StreamingScorer — stateful per-chunk scoring for streaming memory
// ================================================================
//
// Processes one chunk at a time via score_chunk(), allowing Python
// to free each chunk's arrays immediately after the call returns.
// Multimapper groups are scored eagerly (no cross-chunk references).
// Call finish() after all chunks to flush the final MM group and
// retrieve the CSR result arrays.

class StreamingScorer {
    const NativeFragmentScorer& scorer_;

    // Growing output vectors (heap-allocated, capsule-transferred)
    std::vector<int64_t>*  v_offsets_;
    std::vector<int32_t>*  v_ti_;
    std::vector<float>*    v_ll_;
    std::vector<uint8_t>*  v_ct_;
    std::vector<float>*    v_cw_;
    std::vector<int32_t>*  v_lt_;
    std::vector<uint8_t>*  v_lct_;
    std::vector<int8_t>*   v_isp_;
    std::vector<float>*    v_gll_;
    std::vector<int64_t>*  v_gmid_;
    std::vector<int64_t>*  v_fid_;
    std::vector<int8_t>*   v_fc_;
    std::vector<uint8_t>*  v_st_;
    std::vector<int32_t>*  v_dti_;
    std::vector<int64_t>*  v_dfid_;
    std::vector<int64_t>*  v_chim_fid_;
    std::vector<uint8_t>*  v_chim_stype_;

    NativeFragmentScorer::FillState st_;
    double gdna_log_sp_;

    // Statistics
    int64_t stat_det_, stat_em_u_, stat_em_as_, stat_em_ao_;
    int64_t stat_gated_, stat_chim_, stat_mm_;

    bool finished_;

public:
    StreamingScorer(
        NativeFragmentScorer& scorer,
        i8_1d   t_to_strand_arr,
        f64_2d_mut unambig_counts,
        double  gdna_log_splice_pen_unspliced)
      : scorer_(scorer),
        gdna_log_sp_(gdna_log_splice_pen_unspliced),
        stat_det_(0), stat_em_u_(0), stat_em_as_(0),
        stat_em_ao_(0), stat_gated_(0), stat_chim_(0),
        stat_mm_(0), finished_(false)
    {
        // Allocate growing vectors
        v_offsets_ = new std::vector<int64_t>();
        v_offsets_->push_back(0);

        v_ti_  = new std::vector<int32_t>();
        v_ll_  = new std::vector<float>();
        v_ct_  = new std::vector<uint8_t>();
        v_cw_  = new std::vector<float>();
        v_lt_  = new std::vector<int32_t>();
        v_lct_ = new std::vector<uint8_t>();
        v_isp_ = new std::vector<int8_t>();
        v_gll_ = new std::vector<float>();
        v_gmid_ = new std::vector<int64_t>();
        v_fid_ = new std::vector<int64_t>();
        v_fc_  = new std::vector<int8_t>();
        v_st_  = new std::vector<uint8_t>();
        v_dti_ = new std::vector<int32_t>();
        v_dfid_ = new std::vector<int64_t>();
        v_chim_fid_   = new std::vector<int64_t>();
        v_chim_stype_ = new std::vector<uint8_t>();

        // Wire up FillState
        st_ = NativeFragmentScorer::FillState{};
        st_.v_offsets    = v_offsets_;
        st_.v_ti         = v_ti_;
        st_.v_ll         = v_ll_;
        st_.v_ct         = v_ct_;
        st_.v_cw         = v_cw_;
        st_.v_locus_t    = v_lt_;
        st_.v_locus_ct   = v_lct_;
        st_.v_is_spliced = v_isp_;
        st_.v_gdna_ll    = v_gll_;
        st_.v_gmid       = v_gmid_;
        st_.v_fid        = v_fid_;
        st_.v_fclass     = v_fc_;
        st_.v_stype      = v_st_;
        st_.v_det_ti     = v_dti_;
        st_.v_det_fid    = v_dfid_;
        st_.v_chim_fid   = v_chim_fid_;
        st_.v_chim_stype = v_chim_stype_;
        st_.ua_counts    = unambig_counts.data();
        st_.ua_ncols     =
            static_cast<int>(unambig_counts.shape(1));
        st_.t_str        = t_to_strand_arr.data();
    }

    ~StreamingScorer() {
        // If finish() was never called, clean up owned vectors.
        // After finish() these are nullptr (ownership transferred).
        delete v_offsets_;
        delete v_ti_;
        delete v_ll_;
        delete v_ct_;
        delete v_cw_;
        delete v_lt_;
        delete v_lct_;
        delete v_isp_;
        delete v_gll_;
        delete v_gmid_;
        delete v_fid_;
        delete v_fc_;
        delete v_st_;
        delete v_dti_;
        delete v_dfid_;
        delete v_chim_fid_;
        delete v_chim_stype_;
    }

    // Non-copyable
    StreamingScorer(const StreamingScorer&) = delete;
    StreamingScorer& operator=(const StreamingScorer&) = delete;

    void score_chunk(nb::tuple chunk_arrays) {
        if (finished_)
            throw std::runtime_error(
                "StreamingScorer: score_chunk called after finish");

        NativeFragmentScorer::ChunkPtrs cp{};
        cp.t_off  = nb::cast<i32_1d>(chunk_arrays[0]).data();
        cp.t_ind  = nb::cast<i32_1d>(chunk_arrays[1]).data();
        cp.f_len  = nb::cast<i32_1d>(chunk_arrays[2]).data();
        cp.e_bp   = nb::cast<u16_1d>(chunk_arrays[3]).data();
        cp.s_type = nb::cast<u8_1d>(chunk_arrays[4]).data();
        cp.e_str  = nb::cast<u8_1d>(chunk_arrays[5]).data();
        cp.fc     = nb::cast<u8_1d>(chunk_arrays[6]).data();
        cp.f_id   = nb::cast<i64_1d>(chunk_arrays[7]).data();
        cp.r_len  = nb::cast<u16_1d>(chunk_arrays[8]).data();
        cp.g_fp   = nb::cast<i32_1d>(chunk_arrays[9]).data();
        cp.g_sta  = nb::cast<i32_1d>(chunk_arrays[10]).data();
        cp.nm     = nb::cast<u16_1d>(chunk_arrays[11]).data();
        cp.N      = static_cast<int>(
            nb::cast<u8_1d>(chunk_arrays[4]).shape(0));

        scorer_.score_chunk_impl(
            cp, gdna_log_sp_, st_,
            stat_det_, stat_em_u_, stat_em_as_,
            stat_em_ao_, stat_gated_, stat_chim_,
            stat_mm_);
    }

    nb::tuple finish() {
        if (finished_)
            throw std::runtime_error(
                "StreamingScorer: finish called twice");
        finished_ = true;

        // Flush final pending MM group
        if (st_.mm_n_members > 0) {
            scorer_.flush_mm_group(
                st_, gdna_log_sp_,
                stat_mm_, stat_gated_);
        }

        // Build capsule-backed numpy arrays (zero-copy to Python)
        // Transfer ownership to capsules — set pointers to null
        // so destructor doesn't double-free.
        auto result = nb::make_tuple(
            vec_to_ndarray(v_offsets_),
            vec_to_ndarray(v_ti_),
            vec_to_ndarray(v_ll_),
            vec_to_ndarray(v_ct_),
            vec_to_ndarray(v_cw_),
            vec_to_ndarray(v_lt_),
            vec_to_ndarray(v_lct_),
            vec_to_ndarray(v_isp_),
            vec_to_ndarray(v_gll_),
            vec_to_ndarray(v_fid_),
            vec_to_ndarray(v_fc_),
            vec_to_ndarray(v_st_),
            vec_to_ndarray(v_dti_),
            vec_to_ndarray(v_dfid_),
            vec_to_ndarray(v_chim_fid_),
            vec_to_ndarray(v_chim_stype_),
            stat_det_,
            stat_em_u_,
            stat_em_as_,
            stat_em_ao_,
            stat_gated_,
            stat_chim_,
            stat_mm_,
            vec_to_ndarray(v_gmid_)
        );

        // Null out pointers — capsules now own the memory
        v_offsets_ = nullptr;
        v_ti_ = nullptr;  v_ll_ = nullptr;
        v_ct_ = nullptr;  v_cw_ = nullptr;
        v_lt_ = nullptr;  v_lct_ = nullptr;
        v_isp_ = nullptr; v_gll_ = nullptr;
        v_gmid_ = nullptr;
        v_fid_ = nullptr;
        v_fc_ = nullptr;  v_st_ = nullptr;
        v_dti_ = nullptr; v_dfid_ = nullptr;
        v_chim_fid_ = nullptr;
        v_chim_stype_ = nullptr;

        return result;
    }
};

// ----------------------------------------------------------------
// nanobind module definition
// ----------------------------------------------------------------

NB_MODULE(_scoring_impl, m) {
    m.doc() = "C++ hot-path scoring kernels for rigel (nanobind)";

    nb::class_<NativeFragmentScorer>(m, "NativeFragmentScorer")
        .def(nb::init<
                 double, double, bool, double, double,
                 nb::object, int32_t, double,
                 nb::object, int32_t, double,
                 i8_1d, i32_1d,
                 i32_1d, i32_1d, i32_1d, i32_1d,
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
             nb::arg("exon_offsets"),
             nb::arg("exon_starts"),
             nb::arg("exon_ends"),
             nb::arg("exon_cumsum"),
             nb::arg("t_is_nrna_arr"),
             nb::arg("pruning_max_ll_delta"))
;

    nb::class_<StreamingScorer>(m, "StreamingScorer")
        .def(nb::init<
                 NativeFragmentScorer&,
                 i8_1d,
                 f64_2d_mut,
                 double>(),
             nb::arg("scorer"),
             nb::arg("t_to_strand_arr"),
             nb::arg("unambig_counts"),
             nb::arg("gdna_log_splice_pen_unspliced"))
        .def("score_chunk",
             &StreamingScorer::score_chunk,
             nb::arg("chunk_arrays"))
        .def("finish",
             &StreamingScorer::finish);

    // Export scoring constants for Python-side parity tests
    m.attr("LOG_HALF")      = rigel::LOG_HALF;
    m.attr("TAIL_DECAY_LP") = rigel::TAIL_DECAY_LP;
}
