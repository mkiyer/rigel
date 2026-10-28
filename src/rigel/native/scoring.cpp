/**
 * scoring.cpp — nanobind C++ extension for hot-path scoring functions.
 *
 * This module provides the performance-critical scoring kernels. Python
 * builds one NativeFragmentScorer per run (``rigel.scoring``) and feeds the
 * buffer's chunks to a StreamingScorer (``rigel.routing.FragmentRouter``); the
 * per-fragment loop runs here, in C++.  Functions are exposed as
 * ``rigel._scoring_impl``.
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
using rigel::SPLICE_SPLICED_ANNOT;
using rigel::SPLICE_IMPLICIT;
using rigel::SPLICE_ARTIFACT;
using rigel::vec_to_ndarray;
using rigel::FRAG_UNAMBIG;
using rigel::FRAG_AMBIG_SAME_STRAND;
using rigel::FRAG_MULTIMAPPER;

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
using f32_1d  = nb::ndarray<const float,   nb::ndim<1>, nb::c_contig>;
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
    // `--overhang-alpha 0` is a hard gate: log 0 = -inf, and 0 * -inf is NaN, so the gate is carried as a flag and
    // the penalty as 0; a candidate with any overhang is skipped, which is what the -inf would have done
    bool oh_gate_;
    double oh_log_pen_;
    // `--mismatch-alpha 0` is a hard gate too, carried by the penalty on the EM path: nm * -inf is -inf
    // (nm > 0), which prunes an RNA candidate and keeps a gDNA term from competing; a multimapper's sum over its
    // places skips a -inf term, so the order of its hits cannot matter. The deterministic path, which scores
    // nothing, reads the flag
    bool mm_gate_;
    double mm_log_pen_;
    int32_t fl_max_size_;
    double  fl_tail_base_;
    bool    has_fl_lut_;
    int32_t n_transcripts_;

    // --- Copied index arrays ---
    std::vector<int8_t>  t_strand_;   // int8[n_transcripts]
    std::vector<int32_t> t_length_;   // spliced exonic length
    std::vector<uint8_t> t_is_synthetic_;  // uint8[n_transcripts] — 1 = a synthetic span, 0 = an annotated transcript
    std::vector<int32_t> t_ref_;           // int32[n_transcripts] — its reference: a genic place's, for gDNA's places

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
    ///   SPLICED_ANNOT  no: a sequenced junction that survived the blacklist.
    /// Whether gDNA then COMPETES is the likelihood's: its term is the unspliced term at the genomic
    /// footprint whatever the label, and it is pruned exactly as an RNA candidate is (`gdna_competes`).
    inline bool gdna_can_explain(int stype) const {
        return stype == SPLICE_UNSPLICED || stype == SPLICE_ARTIFACT || stype == SPLICE_IMPLICIT;
    }

    /// The locus's gDNA candidate is kept when its likelihood is within the pruning bound of the unit's best RNA
    /// candidate — `max_ll_delta_ = -log(pruning_min_posterior)`, the bound every RNA candidate meets within
    /// its pool — and is otherwise no candidate at all. A footprint the gDNA length law cannot produce is
    /// pruned by this, with no length rule of its own. The intergenic term is never pruned (`flush_mm_group`).
    inline bool gdna_competes(double gdna_ll, double best_rna_ll) const {
        return best_rna_ll - gdna_ll <= max_ll_delta_;
    }

    /// The fragment's coverage weight on transcript ``t_idx`` (``compute_fragment_weight``), read on the interval it
    /// occupies in the transcript's OWN 5'->3' coordinates, cut to the transcript. Its start is projected onto the
    /// transcript: inside an exon at the exon's offset; in an intron, or before the first exon, the overhang's length
    /// before the next exon's start; past the last exon that far past the end. It extends ``flen``, the length of the
    /// path the transcript implies. So it occupies ``[fwd, fwd + flen)`` forward, and on a minus-strand transcript,
    /// whose 5'->3' runs against the genome, the mirror ``[t_len - fwd - flen, t_len - fwd)``. A fragment and its
    /// mirror image on the other strand weigh the same. Held to an independent evaluation of the coverage model on
    /// both strands in ``tests/test_pipeline_routing.py``.
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
        u8_1d   t_is_synthetic_arr,
        i32_1d  t_ref_arr,
        double  pruning_max_ll_delta)
      : log_p_sense_(log_p_sense),
        log_p_antisense_(log_p_antisense),
        r1_antisense_(r1_antisense),
        oh_gate_(std::isinf(oh_log_pen) && oh_log_pen < 0.0),
        oh_log_pen_(std::isinf(oh_log_pen) && oh_log_pen < 0.0 ? 0.0 : oh_log_pen),
        mm_gate_(std::isinf(mm_log_pen) && mm_log_pen < 0.0),
        mm_log_pen_(mm_log_pen),
        fl_max_size_(fl_max_size),
        fl_tail_base_(fl_tail_base),
        has_fl_lut_(false),
        max_ll_delta_(pruning_max_ll_delta),
        gdna_fl_max_size_(gdna_fl_max_size),
        gdna_fl_tail_base_(gdna_fl_tail_base),
        has_gdna_fl_lut_(false)
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
            const auto* p = t_is_synthetic_arr.data();
            t_is_synthetic_.assign(p, p + n_transcripts_);
        }
        {
            const auto* p = t_ref_arr.data();
            t_ref_.assign(p, p + n_transcripts_);
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
        // The log capture weight of each row's place, applied to a multimapper's every candidate at that place;
        // nullptr reads 0 (a unique fragment's is common to its candidates and never applied).
        const float*    lw;
        int N;
    };

    // A component's term at one hit of a multimapper, keyed by the component — a transcript (annotated or a synthetic
    // span alike), or for the locus's gDNA the hit's reference — with the hit's genomic extent `[start, end)`, which
    // tells the component's places apart (`merge_places`).
    struct PlaceTerm {
        int32_t key;
        int32_t start, end;
        double  log_lik;
        int32_t count_col;
        double  coverage_wt;
    };

    // One key's term over all of a multimapper's hits, and what the unit records of its best hit.
    struct MergedCandidate {
        int32_t key;
        double  log_lik;
        int32_t count_col;
        double  coverage_wt;
    };

    // Each key's term is the SUM over its places, RNA's and gDNA's alike: hits whose extents overlap are alternatives
    // at one place and count once, at the best of them; distinct places add. A transcript's unit records the count
    // column and coverage weight of its best hit, the first in input order among equals. `terms` is sorted in place.
    static void merge_places(std::vector<PlaceTerm>& terms, std::vector<MergedCandidate>& out) {
        static constexpr double NEG_INF = -std::numeric_limits<double>::infinity();
        out.clear();
        std::stable_sort(terms.begin(), terms.end(), [](const PlaceTerm& a, const PlaceTerm& b) {
            return a.key != b.key ? a.key < b.key : a.start < b.start;
        });
        for (size_t i = 0; i < terms.size();) {
            const int32_t key = terms[i].key;
            const PlaceTerm* best = &terms[i];
            double lse_max = 0.0, lse_sum = 0.0;
            bool lse_has = false;
            while (i < terms.size() && terms[i].key == key) {
                // one place: the hits chained by overlapping extents, read at the best of them
                double place = terms[i].log_lik;
                int32_t end = terms[i].end;
                if (terms[i].log_lik > best->log_lik) best = &terms[i];
                for (++i; i < terms.size() && terms[i].key == key && terms[i].start < end; ++i) {
                    place = std::max(place, terms[i].log_lik);
                    end = std::max(end, terms[i].end);
                    if (terms[i].log_lik > best->log_lik) best = &terms[i];
                }
                if (place > NEG_INF) lse_update(lse_max, lse_sum, lse_has, place);
            }
            out.push_back({key, lse_has ? lse_max + std::log(lse_sum) : NEG_INF, best->count_col, best->coverage_wt});
        }
    }

    // The log of the sum over every place of every key: gDNA's term, its keys the references its places are on.
    // `per_key` receives each key's own sum.
    static double sum_over_places(std::vector<PlaceTerm>& terms, std::vector<MergedCandidate>& per_key) {
        merge_places(terms, per_key);
        double lse_max = 0.0, lse_sum = 0.0;
        bool lse_has = false;
        for (const auto& mc : per_key)
            if (mc.log_lik > -std::numeric_limits<double>::infinity()) lse_update(lse_max, lse_sum, lse_has, mc.log_lik);
        return lse_has ? lse_max + std::log(lse_sum) : -std::numeric_limits<double>::infinity();
    }

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
        std::vector<int8_t>*   v_is_spliced;
        std::vector<float>*    v_gdna_ll;
        std::vector<float>*    v_ig_ll;
        std::vector<int64_t>*  v_fid;
        std::vector<int8_t>*   v_fclass;
        std::vector<uint8_t>*  v_stype;
        // Det-unambig (growing)
        std::vector<int32_t>*  v_det_ti;
        std::vector<int64_t>*  v_det_fid;
        // The fragments a hard gate dropped, and their classes: no pool holds them, and the annotated BAM says so
        std::vector<int64_t>*  v_gated_fid;
        std::vector<int8_t>*   v_gated_fclass;
        // The multimappers whose only places the gates left are intergenic: intergenic gDNA, in no locus
        std::vector<int64_t>*  v_intergenic_fid;
        // Unambig counts accumulation (single-pass)
        double* ua_counts;
        int ua_ncols;
        const int8_t* t_str;
        // Candidate write cursor
        int64_t cand_cur;

        // Eager MM group accumulators (persists across chunks).
        // Candidates are scored eagerly as alignments arrive;
        // flush_mm_group() merges them over places, prunes and emits.
        int64_t mm_fid;
        int     mm_n_members;
        std::vector<PlaceTerm> mm_terms;
        std::vector<MergedCandidate> mm_merged;
        int     mm_best_stype;

        // gDNA's terms at the places it can explain: at a place with a transcript candidate the locus's, merged over
        // places at the flush; at a place with none the intergenic term's, an online log-sum-exp (no reference is
        // known there, so each is its own place).
        std::vector<PlaceTerm> mm_gdna_terms;
        std::vector<MergedCandidate> mm_gdna_merged;
        double  mm_ig_lse_max;
        double  mm_ig_lse_sum;
        bool    mm_ig_lse_has;

        FillState()
            : v_ti(nullptr), v_ll(nullptr), v_ct(nullptr),
                            v_cw(nullptr),
              v_offsets(nullptr), v_is_spliced(nullptr),
              v_gdna_ll(nullptr), v_ig_ll(nullptr), v_fid(nullptr),
              v_fclass(nullptr), v_stype(nullptr),
              v_det_ti(nullptr), v_det_fid(nullptr),
              v_gated_fid(nullptr), v_gated_fclass(nullptr), v_intergenic_fid(nullptr),
              ua_counts(nullptr), ua_ncols(0), t_str(nullptr),
              cand_cur(0),
              mm_fid(-1), mm_n_members(0),
              mm_best_stype(SPLICE_UNSPLICED),
              mm_ig_lse_max(0.0),
              mm_ig_lse_sum(0.0),
              mm_ig_lse_has(false) {}

        void reset_mm_group() {
            mm_terms.clear();  // all four keep their allocation
            mm_merged.clear();
            mm_gdna_terms.clear();
            mm_gdna_merged.clear();
            mm_fid = -1;
            mm_n_members = 0;
            mm_best_stype = SPLICE_UNSPLICED;
            mm_ig_lse_max = 0.0;
            mm_ig_lse_sum = 0.0;
            mm_ig_lse_has = false;
        }
    };

    // ---------------------------------------------------------------
    // score_mm_alignment — eagerly score one MM alignment
    // ---------------------------------------------------------------
    //
    // Called during score_chunk_impl() for each MM alignment.
    // Reads from the current chunk's live arrays and appends its scored
    // candidates to st.mm_terms.  Also updates the per-group
    // splice type and gDNA accumulators.

    void score_mm_alignment(
        const ChunkPtrs& cp, int row,
        FillState& st) const
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
        // The place's capture weight, on every candidate at it, RNA's and gDNA's alike.
        const double log_w = cp.lw ? static_cast<double>(cp.lw[row]) : 0.0;
        // The hit's genomic extent, which tells a component's places apart.
        int32_t gfp_val = cp.g_fp[row];
        const int32_t hit_start = genomic_start_val;
        const int32_t hit_end = genomic_start_val + std::max(gfp_val, 1);

        // Update per-group metadata accumulators (the reported splice type only)
        if (stype == SPLICE_SPLICED_ANNOT) st.mm_best_stype = SPLICE_SPLICED_ANNOT;

        // gDNA's term at this place, from every hit gDNA can explain (`gdna_can_explain`), SUMMED over the
        // places: the molecule came from one of them. A place with a transcript candidate is the locus's, on its
        // candidates' reference; a place with none is intergenic, in no locus, and its sum is the intergenic term.
        // A term at -inf (the mismatch gate, or a length the law cannot produce) adds nothing, so it is skipped and
        // the hits' order cannot matter. The EM applies each component's rate and length.
        if (gdna_can_explain(stype))
        {
            double hit_log_ll = gdna_frag_len_log_lik(gfp_val) + LOG_HALF + log_nm + log_w;
            if (hit_log_ll > -std::numeric_limits<double>::infinity()) {
                if (n_cand > 0)
                    st.mm_gdna_terms.push_back({t_ref_[cp.t_ind[start]], hit_start, hit_end, hit_log_ll, 0, 0.0});
                else
                    lse_update(st.mm_ig_lse_max, st.mm_ig_lse_sum, st.mm_ig_lse_has, hit_log_ll);
            }
        }

        // Count the group's alignments — read only as the group's non-empty test, before a flush
        ++st.mm_n_members;

        // Score every transcript candidate, annotated or the synthetic span
        if (n_cand <= 0) return;

        for (int64_t k = start; k < end; ++k) {
            int32_t t_idx = cp.t_ind[k];
            int32_t ebp   = cp.e_bp[k];
            if (ebp <= 0) continue;

            int32_t oh = rl - ebp;
            if (oh < 0) oh = 0;
            if (oh_gate_ && oh > 0) continue;

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
                           + oh * oh_log_pen_ + log_nm + log_w;
            int32_t ct = stype * 2 + (is_anti ? 1 : 0);

            double cov_wt = (has_genomic && flen > 0)
                ? coverage_weight(genomic_start_val, flen, t_idx) : 1.0;

            st.mm_terms.push_back({t_idx, hit_start, hit_end, log_lik, ct, cov_wt});
        }
    }

    // ---------------------------------------------------------------
    // flush_mm_group — prune and emit one multimapper group
    // ---------------------------------------------------------------
    //
    // All scoring has already been done eagerly by score_mm_alignment. This method merges each transcript's terms
    // over its places (`merge_places`), prunes them, decides the group's outcome once every place has been read, and
    // resets the accumulators for the next group:
    //   - an EM unit, when an RNA candidate survives: its transcripts hold its locus;
    //   - intergenic gDNA, when none does but an intergenic place does: no locus holds it, and it is booked as a
    //     multimapper none of whose hits resolves is (its gDNA at a genic place has no locus without a candidate);
    //   - gated, when no place survives.

    void flush_mm_group(
        FillState& st,
        int64_t& stat_mm,
        int64_t& stat_gated) const
    {
        static constexpr double NEG_INF =
            -std::numeric_limits<double>::infinity();

        if (st.mm_n_members == 0) return;

        merge_places(st.mm_terms, st.mm_merged);

        // Pool-separated likelihood pruning: annotated transcripts (`annotated_best`) and synthetic spans
        // (`synthetic_best`, by `t_is_synthetic_`) are each pruned against their own pool's best
        double annotated_best = NEG_INF;
        double synthetic_best = NEG_INF;
        for (const auto& mc : st.mm_merged) {
            if (t_is_synthetic_[mc.key])
                synthetic_best = std::max(synthetic_best, mc.log_lik);
            else
                annotated_best = std::max(annotated_best, mc.log_lik);
        }

        int64_t emit_start = st.cand_cur;
        double best_ll = NEG_INF;

        for (const auto& mc : st.mm_merged) {
            const int32_t t_idx = mc.key;
            double pool_best = t_is_synthetic_[t_idx]
                             ? synthetic_best : annotated_best;
            if (pool_best - mc.log_lik <= max_ll_delta_) {
                st.v_ti->push_back(t_idx);
                st.v_ll->push_back(static_cast<float>(mc.log_lik));
                st.v_ct->push_back(
                    static_cast<uint8_t>(mc.count_col));
                st.v_cw->push_back(static_cast<float>(mc.coverage_wt));
                ++st.cand_cur;
                if (mc.log_lik > best_ll) best_ll = mc.log_lik;
            }
        }

        // The group's two gDNA readings, each the sum over its places: the locus's, and the intergenic term no locus
        // carries. -inf where no place of the kind survived its gates.
        const double gdna_ll = sum_over_places(st.mm_gdna_terms, st.mm_gdna_merged);
        const double ig_ll = st.mm_ig_lse_has ? st.mm_ig_lse_max + std::log(st.mm_ig_lse_sum) : NEG_INF;

        int64_t new_cands = st.cand_cur - emit_start;
        if (new_cands > 0) {
            st.v_offsets->push_back(st.cand_cur);
            st.v_fid->push_back(st.mm_fid);
            st.v_fclass->push_back(
                static_cast<int8_t>(FRAG_MULTIMAPPER));

            st.v_stype->push_back(
                static_cast<uint8_t>(st.mm_best_stype));

            // The locus's gDNA is kept when it competes with the group's best RNA candidate. The intergenic term is
            // kept whenever a place survived: it competes at the intergenic rate against the locus's rates and
            // lengths, none of which are in these scores, so no difference of scores bounds its posterior.
            const bool has_gdna = gdna_competes(gdna_ll, best_ll);
            const bool has_ig = ig_ll > NEG_INF;
            // "spliced" = gDNA is no candidate anywhere: the group is certified RNA.
            st.v_is_spliced->push_back((has_gdna || has_ig) ? 0 : 1);
            st.v_gdna_ll->push_back(static_cast<float>(has_gdna ? gdna_ll : NEG_INF));
            st.v_ig_ll->push_back(static_cast<float>(ig_ll));

            ++stat_mm;
        } else if (ig_ll > NEG_INF) {
            st.v_intergenic_fid->push_back(st.mm_fid);
        } else {
            ++stat_gated;
            st.v_gated_fid->push_back(st.mm_fid);
            st.v_gated_fclass->push_back(static_cast<int8_t>(FRAG_MULTIMAPPER));
        }

        st.reset_mm_group();
    }

    // Single-pass inner loop — scores fragments and pushes results
    // to growing output vectors.  Also accumulates ua_counts for
    // det-unambig fragments.
    void score_chunk_impl(
        const ChunkPtrs& cp,
        FillState& st,
        int64_t& stat_det, int64_t& stat_em_u,
        int64_t& stat_em_as, int64_t& stat_em_ao,
        int64_t& stat_gated,
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

            // ---- Multimapper: eager scoring ----
            if (fclass == FRAG_MULTIMAPPER) {
                int64_t fid_val = f_id[i];
                if (fid_val != st.mm_fid) {
                    if (st.mm_n_members > 0) {
                        flush_mm_group(
                            st,
                            stat_mm, stat_gated);
                    }
                    st.mm_fid = fid_val;
                }
                score_mm_alignment(cp, i, st);
                continue;
            }

            // Non-MM fragment: flush any pending MM group
            if (st.mm_n_members > 0) {
                flush_mm_group(
                    st,
                    stat_mm, stat_gated);
            }

            int stype      = s_type[i];
            int exon_str   = e_str[i];
            int64_t start  = t_off[i];
            int64_t end    = t_off[i + 1];
            int n_cand     = static_cast<int>(end - start);

            // A fragment with no candidate is skipped and not counted in stat_gated.
            if (n_cand <= 0) continue;

            // ---- Det-unambig: both ua_counts + det arrays ----
            if (fclass == FRAG_UNAMBIG &&
                stype == SPLICE_SPLICED_ANNOT) {
                int32_t t_idx = t_ind[start];

                // The hard gates hold on this path too: at alpha 0 a mismatch, or an overhang past the one
                // candidate's exons, excludes it, so the fragment is gated exactly as an EM-routed one would be
                const int rl_det = r_len[i] > 0 ? static_cast<int>(r_len[i]) : 1;
                if ((mm_gate_ && nm_arr[i] > 0) || (oh_gate_ && rl_det - e_bp[start] > 0)) {
                    ++stat_gated;
                    st.v_gated_fid->push_back(f_id[i]);
                    st.v_gated_fclass->push_back(static_cast<int8_t>(fclass));
                    continue;
                }

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
                ++stat_det;
                continue;
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

            // ===== candidate scoring: annotated transcripts and synthetic spans alike =====
            struct ScoredCandidate {
                int32_t t_idx, oh, ct;
                double  log_lik, cov_wt;
            };
            ScoredCandidate scored_stack[SCORED_STACK_CAPACITY];
            std::vector<ScoredCandidate> scored_heap;
            ScoredCandidate* scored = scored_stack;
            if (n_cand > SCORED_STACK_CAPACITY) {
                scored_heap.resize(n_cand);
                scored = scored_heap.data();
            }
            int n_scored = 0;

            for (int64_t k = start; k < end; ++k) {
                int32_t t_idx = t_ind[k];
                int32_t ebp   = e_bp[k];
                if (ebp <= 0) continue;

                int32_t oh = rl - ebp;
                if (oh < 0) oh = 0;
                if (oh_gate_ && oh > 0) continue;

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

                scored[n_scored++] = {t_idx, oh, ct, log_lik, cov_wt};
            }

            // Pool-separated likelihood pruning: annotated transcripts (`annotated_best`) and synthetic spans
            // (`synthetic_best`, by `t_is_synthetic_`) are each pruned against their own pool's best
            double annotated_best = NEG_INF;
            double synthetic_best = NEG_INF;
            for (int j = 0; j < n_scored; ++j) {
                if (t_is_synthetic_[scored[j].t_idx])
                    synthetic_best = std::max(synthetic_best, scored[j].log_lik);
                else
                    annotated_best = std::max(annotated_best, scored[j].log_lik);
            }

            for (int j = 0; j < n_scored; ++j) {
                auto& s = scored[j];
                double pool_best = t_is_synthetic_[s.t_idx]
                                 ? synthetic_best : annotated_best;
                if (pool_best - s.log_lik <= max_ll_delta_) {
                    st.v_ti->push_back(s.t_idx);
                    st.v_ll->push_back(static_cast<float>(s.log_lik));
                    st.v_ct->push_back(
                        static_cast<uint8_t>(s.ct));
                    st.v_cw->push_back(static_cast<float>(s.cov_wt));
                    ++st.cand_cur;
                    if (s.log_lik > best_ll) best_ll = s.log_lik;
                }
            }

            // ---- Finalize this EM unit ----
            int64_t new_cands = st.cand_cur - emit_start;
            if (new_cands > 0) {
                st.v_offsets->push_back(st.cand_cur);
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
                if (gdna_can_explain(stype)) {
                    gdna_ll = gdna_frag_len_log_lik(genomic_footprint) + LOG_HALF + log_nm;
                    has_gdna = gdna_competes(gdna_ll, best_ll);
                }
                // "spliced" = gDNA is no candidate: certified RNA.
                st.v_is_spliced->push_back(
                    has_gdna ? 0 : 1);

                if (has_gdna) {
                    st.v_gdna_ll->push_back(static_cast<float>(gdna_ll));
                } else {
                    st.v_gdna_ll->push_back(static_cast<float>(NEG_INF));
                }
                // One place: no intergenic reading beside a resolved hit.
                st.v_ig_ll->push_back(static_cast<float>(NEG_INF));

                if (fclass == FRAG_UNAMBIG)
                    ++stat_em_u;
                else if (fclass == FRAG_AMBIG_SAME_STRAND)
                    ++stat_em_as;
                else
                    ++stat_em_ao;
            } else {
                ++stat_gated;
                st.v_gated_fid->push_back(f_id[i]);
                st.v_gated_fclass->push_back(static_cast<int8_t>(fclass));
            }
        }
    }
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
    std::vector<int8_t>*   v_isp_;
    std::vector<float>*    v_gll_;
    std::vector<float>*    v_igll_;
    std::vector<int64_t>*  v_fid_;
    std::vector<int8_t>*   v_fc_;
    std::vector<uint8_t>*  v_st_;
    std::vector<int32_t>*  v_dti_;
    std::vector<int64_t>*  v_dfid_;
    std::vector<int64_t>*  v_gfid_;
    std::vector<int8_t>*   v_gfc_;
    std::vector<int64_t>*  v_igfid_;

    NativeFragmentScorer::FillState st_;

    // Statistics
    int64_t stat_det_, stat_em_u_, stat_em_as_, stat_em_ao_;
    int64_t stat_gated_, stat_mm_;

    bool finished_;

public:
    StreamingScorer(
        NativeFragmentScorer& scorer,
        i8_1d   t_to_strand_arr,
        f64_2d_mut unambig_counts)
      : scorer_(scorer),
        stat_det_(0), stat_em_u_(0), stat_em_as_(0),
        stat_em_ao_(0), stat_gated_(0),
        stat_mm_(0), finished_(false)
    {
        // Allocate growing vectors
        v_offsets_ = new std::vector<int64_t>();
        v_offsets_->push_back(0);

        v_ti_  = new std::vector<int32_t>();
        v_ll_  = new std::vector<float>();
        v_ct_  = new std::vector<uint8_t>();
        v_cw_  = new std::vector<float>();
        v_isp_ = new std::vector<int8_t>();
        v_gll_ = new std::vector<float>();
        v_igll_ = new std::vector<float>();
        v_fid_ = new std::vector<int64_t>();
        v_fc_  = new std::vector<int8_t>();
        v_st_  = new std::vector<uint8_t>();
        v_dti_ = new std::vector<int32_t>();
        v_dfid_ = new std::vector<int64_t>();
        v_gfid_ = new std::vector<int64_t>();
        v_gfc_ = new std::vector<int8_t>();
        v_igfid_ = new std::vector<int64_t>();

        // Wire up FillState
        st_ = NativeFragmentScorer::FillState{};
        st_.v_offsets    = v_offsets_;
        st_.v_ti         = v_ti_;
        st_.v_ll         = v_ll_;
        st_.v_ct         = v_ct_;
        st_.v_cw         = v_cw_;
        st_.v_is_spliced = v_isp_;
        st_.v_gdna_ll    = v_gll_;
        st_.v_ig_ll      = v_igll_;
        st_.v_fid        = v_fid_;
        st_.v_fclass     = v_fc_;
        st_.v_stype      = v_st_;
        st_.v_det_ti     = v_dti_;
        st_.v_det_fid    = v_dfid_;
        st_.v_gated_fid    = v_gfid_;
        st_.v_gated_fclass = v_gfc_;
        st_.v_intergenic_fid = v_igfid_;
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
        delete v_isp_;
        delete v_gll_;
        delete v_igll_;
        delete v_fid_;
        delete v_fc_;
        delete v_st_;
        delete v_dti_;
        delete v_dfid_;
        delete v_gfid_;
        delete v_gfc_;
        delete v_igfid_;
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
        cp.lw     = nb::len(chunk_arrays) > 12 ? nb::cast<f32_1d>(chunk_arrays[12]).data() : nullptr;
        cp.N      = static_cast<int>(
            nb::cast<u8_1d>(chunk_arrays[4]).shape(0));

        scorer_.score_chunk_impl(
            cp, st_,
            stat_det_, stat_em_u_, stat_em_as_,
            stat_em_ao_, stat_gated_,
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
                st_,
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
            vec_to_ndarray(v_isp_),
            vec_to_ndarray(v_gll_),
            vec_to_ndarray(v_fid_),
            vec_to_ndarray(v_fc_),
            vec_to_ndarray(v_st_),
            vec_to_ndarray(v_dti_),
            vec_to_ndarray(v_dfid_),
            vec_to_ndarray(v_gfid_),
            vec_to_ndarray(v_gfc_),
            stat_det_,
            stat_em_u_,
            stat_em_as_,
            stat_em_ao_,
            stat_gated_,
            stat_mm_,
            vec_to_ndarray(v_igll_),
            vec_to_ndarray(v_igfid_)
        );

        // Null out pointers — capsules now own the memory
        v_offsets_ = nullptr;
        v_ti_ = nullptr;  v_ll_ = nullptr;
        v_ct_ = nullptr;  v_cw_ = nullptr;
        v_isp_ = nullptr; v_gll_ = nullptr; v_igll_ = nullptr;
        v_fid_ = nullptr;
        v_fc_ = nullptr;  v_st_ = nullptr;
        v_dti_ = nullptr; v_dfid_ = nullptr;
        v_gfid_ = nullptr; v_gfc_ = nullptr; v_igfid_ = nullptr;

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
                 u8_1d, i32_1d, double>(),
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
             nb::arg("t_is_synthetic_arr"),
             nb::arg("t_ref_arr"),
             nb::arg("pruning_max_ll_delta"))
;

    nb::class_<StreamingScorer>(m, "StreamingScorer")
        .def(nb::init<
                 NativeFragmentScorer&,
                 i8_1d,
                 f64_2d_mut>(),
             nb::arg("scorer"),
             nb::arg("t_to_strand_arr"),
             nb::arg("unambig_counts"))
        .def("score_chunk",
             &StreamingScorer::score_chunk,
             nb::arg("chunk_arrays"))
        .def("finish",
             &StreamingScorer::finish);

    // Export scoring constants for Python-side parity tests
    m.attr("LOG_HALF")      = rigel::LOG_HALF;
    m.attr("TAIL_DECAY_LP") = rigel::TAIL_DECAY_LP;

    // The fragment classes this scorer routes on, the one definition rigel.buffer reads
    m.attr("FRAG_UNAMBIG")           = rigel::FRAG_UNAMBIG;
    m.attr("FRAG_AMBIG_SAME_STRAND") = rigel::FRAG_AMBIG_SAME_STRAND;
    m.attr("FRAG_AMBIG_OPP_STRAND")  = rigel::FRAG_AMBIG_OPP_STRAND;
    m.attr("FRAG_MULTIMAPPER")       = rigel::FRAG_MULTIMAPPER;
}
