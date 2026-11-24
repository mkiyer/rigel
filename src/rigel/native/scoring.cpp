// scoring.cpp — the fragment scorer, the native module `_scoring_impl`.
//
// Python builds one NativeFragmentScorer per run (`rigel.scoring.build_scorer`) and feeds the buffer's chunks to a
// StreamingScorer (`rigel.routing.FragmentRouter`), which scores every fragment's candidates and emits the EM's CSR
// units.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "constants.h"
#include "ndarray_util.h"
#include "testing/bind_testing.h"

namespace nb = nanobind;

using rigel::LOG_HALF;
using rigel::TAIL_DECAY_LP;
using rigel::STRAND_POS;
using rigel::STRAND_NEG;
using rigel::SPLICE_UNSPLICED;
using rigel::SPLICE_SPLICED_ANNOT;
using rigel::SPLICE_SPLICED_IMPLICIT;
using rigel::vec_to_ndarray;
using rigel::FRAG_UNAMBIG;
using rigel::FRAG_AMBIG_SAME_STRAND;
using rigel::FRAG_MULTIMAPPER;
using rigel::SCORED_STACK_CAPACITY;
using rigel::Array1;

using rigel::NEG_INF;

// An online log-sum-exp: keeps (max, sum) so that logsumexp(x_1..x_n) = max + log(sum); -inf over no terms.
struct LogSumExp {
    double max = 0.0, sum = 0.0;
    bool has = false;

    void add(double x) {
        if (!has) {
            max = x;
            sum = 1.0;
            has = true;
        } else if (x > max) {
            sum = sum * std::exp(max - x) + 1.0;
            max = x;
        } else {
            sum += std::exp(x - max);
        }
    }
    double value() const { return has ? max + std::log(sum) : NEG_INF; }
};

using i32_1d  = Array1<int32_t>;
using i8_1d   = Array1<int8_t>;
using f64_1d  = Array1<double>;
using u8_1d   = Array1<uint8_t>;
using f64_2d_out = nb::ndarray<double, nb::ndim<2>, nb::c_contig>;

// A chunk's array, by name. It must be of the exact dtype: its data is read after the cast that checked it.
template <typename T>
static const T* chunk_array(nb::dict arrays, const char* name) {
    return nb::cast<Array1<T>>(arrays[name], false).data();
}

// An array's elements, copied.
template <typename T>
static std::vector<T> copy_of(const Array1<T>& a) {
    return std::vector<T>(a.data(), a.data() + a.shape(0));
}

// Inverse coverage-capacity weight of a fragment on a transcript, under the trapezoid coverage model: with uniform
// sampling, the expected coverage at position x of a transcript of spliced length L, for a fragment of length f, is
//   min(x, w, L - x),   w = min(f, L/2).
// A fragment on the plateau weighs 1; one near an end weighs more, being less likely to be generated yet observed.
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

// A fragment-length log-probability table: its entries up to max_size, then a fixed decay per base beyond it. An
// absent table reads 0, as does a length of 0 or less.
struct LengthTable {
    std::vector<double> log_prob;
    int32_t max_size = 0;
    double  tail_base = 0.0;
    bool    present = false;

    LengthTable(nb::object table, int32_t max_size_, double tail_base_)
        : max_size(max_size_), tail_base(tail_base_), present(!table.is_none()) {
        if (present) log_prob = copy_of(nb::cast<f64_1d>(table));
    }

    inline double log_lik(int32_t flen) const {
        if (flen <= 0 || !present) return 0.0;
        if (flen <= max_size)
            return log_prob[static_cast<size_t>(flen)];
        return tail_base + (flen - max_size) * TAIL_DECAY_LP;
    }
};

// The scoring parameters and index arrays, copied once per run, so scoring and pruning make no Python round-trip.
class NativeFragmentScorer {

    double log_p_sense_;
    double log_p_antisense_;
    bool   r1_antisense_;
    // An overhang alpha of 0 is a hard gate. Its log is -inf and 0 * -inf is NaN, so the gate is a flag and the
    // penalty 0: a candidate with any overhang is skipped.
    bool oh_gate_;
    double oh_log_pen_;
    // Log penalties per NM mismatch and per soft-clipped base (`rigel.scoring.sequence_log_penalties`). Read only on a
    // multimapper's hits, which they tell apart; a single hit's sequence is common to its candidates and cancels. In a
    // library with no mismatches a mismatch is -inf and the hit drops out of the sum over places.
    double mismatch_log_pen_;
    double clip_log_pen_;
    // Fragment-length log-probability tables, RNA's and gDNA's
    LengthTable rna_fl_, gdna_fl_;

    // Per transcript
    std::vector<int8_t>  t_strand_;
    std::vector<int32_t> t_length_;        // spliced exonic length
    std::vector<uint8_t> t_is_synthetic_;  // 1 = a synthetic span, 0 = an annotated transcript
    std::vector<int32_t> t_ref_;           // its reference, which keys gDNA's places

    double max_ll_delta_;             // the pruning bound, -log(pruning_min_posterior)

    // Exons in CSR layout, for mapping a genomic position onto a transcript
    std::vector<int32_t> exon_offsets_;  // [n_transcripts + 1]
    std::vector<int32_t> exon_starts_;
    std::vector<int32_t> exon_ends_;
    std::vector<int32_t> exon_cumsum_;   // exonic bases before each exon

    /// Whether gDNA can have produced an alignment, for a single fragment and every hit of a multimapper alike:
    ///   UNSPLICED      yes.
    ///   IMPLICIT       yes: its mate gap is unobserved, so gDNA reads it unbroken at the footprint and a transcript
    ///                  at its own length, and the fragment length tells them apart.
    ///   SPLICED_ANNOT  no: it carries a sequenced junction.
    /// Whether gDNA then competes is the likelihood's to decide (`gdna_competes`).
    inline bool gdna_can_explain(int stype) const {
        return stype == SPLICE_UNSPLICED || stype == SPLICE_SPLICED_IMPLICIT;
    }

    /// The locus's gDNA candidate is kept when its likelihood is within `max_ll_delta_` of the best RNA candidate at
    /// its place, the bound every RNA candidate meets within its pool there (a unique fragment has one place). A
    /// footprint the gDNA length law cannot produce is pruned here, with no length rule of its own. The intergenic
    /// term is never pruned (`flush_mm_group`).
    inline bool gdna_competes(double gdna_ll, double best_rna_ll) const {
        return best_rna_ll - gdna_ll <= max_ll_delta_;
    }

    /// The fragment's coverage weight on transcript ``t_idx`` (``compute_fragment_weight``), on the interval it
    /// occupies in the transcript's own 5'->3' coordinates, cut to the transcript. Its start projects inside an exon to
    /// the exon's offset, in an intron or before the first exon to that many bases before the next exon's start, and
    /// past the last exon to that many bases past the end. It spans ``flen``, the path length the transcript implies:
    /// ``[fwd, fwd + flen)`` on a plus-strand transcript and the mirror ``[t_len - fwd - flen, t_len - fwd)`` on a
    /// minus-strand one, so a fragment and its mirror image weigh the same.
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

    // Copies every input into storage the scorer owns.
    NativeFragmentScorer(
        double log_p_sense,
        double log_p_antisense,
        bool   r1_antisense,
        double oh_log_pen,
        double mismatch_log_pen,
        double clip_log_pen,
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
        mismatch_log_pen_(mismatch_log_pen),
        clip_log_pen_(clip_log_pen),
        rna_fl_(fl_log_prob_obj, fl_max_size, fl_tail_base),
        gdna_fl_(gdna_fl_log_prob_obj, gdna_fl_max_size, gdna_fl_tail_base),
        t_strand_(copy_of(t_strand_arr)),
        t_length_(copy_of(t_length_arr)),
        t_is_synthetic_(copy_of(t_is_synthetic_arr)),
        t_ref_(copy_of(t_ref_arr)),
        max_ll_delta_(pruning_max_ll_delta),
        exon_offsets_(copy_of(exon_offsets_arr)),
        exon_starts_(copy_of(exon_starts_arr)),
        exon_ends_(copy_of(exon_ends_arr)),
        exon_cumsum_(copy_of(exon_cumsum_arr)) {}

    friend class StreamingScorer;

private:

    // Raw pointers into one chunk's arrays.
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
        const uint16_t* clip;
        // Each row's log place weight, applied to every multimapper candidate at that place; nullptr reads 0. A unique
        // fragment's is common to its candidates and never applied.
        const float*    lw;
        int N;
    };

    // A component's term at one hit of a multimapper. The key is a transcript (annotated or synthetic span), or for the
    // locus's gDNA the hit's reference; the genomic extent `[start, end)` tells the component's places apart.
    struct PlaceTerm {
        int32_t key;
        int32_t start, end;
        double  log_lik;
        int32_t count_col;
        double  coverage_wt;
    };

    // A candidate's term (a transcript's, or for gDNA a reference's) and what the unit records of it: one hit's, or a
    // multimapper's summed over its places with the count column and weight of its best hit.
    struct Candidate {
        int32_t key;
        double  log_lik;
        int32_t count_col;
        double  coverage_wt;
    };

    // Each key's term is the sum over its places, for RNA and gDNA alike: hits whose extents overlap are alternatives
    // at one place and count once, at the best of them; distinct places add. The unit records the count column and
    // coverage weight of the best hit, the first in input order among equals. Sorts `terms` in place.
    static void merge_places(std::vector<PlaceTerm>& terms, std::vector<Candidate>& out) {
        out.clear();
        std::stable_sort(terms.begin(), terms.end(), [](const PlaceTerm& a, const PlaceTerm& b) {
            return a.key != b.key ? a.key < b.key : a.start < b.start;
        });
        for (size_t i = 0; i < terms.size();) {
            const int32_t key = terms[i].key;
            const PlaceTerm* best = &terms[i];
            LogSumExp places;
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
                if (place > NEG_INF) places.add(place);
            }
            out.push_back({key, places.value(), best->count_col, best->coverage_wt});
        }
    }

    // gDNA's term: the log of the sum over every place of every key, the keys being references. `per_key` receives
    // each key's own sum.
    static double sum_over_places(std::vector<PlaceTerm>& terms, std::vector<Candidate>& per_key) {
        merge_places(terms, per_key);
        LogSumExp all;
        for (const auto& mc : per_key)
            if (mc.log_lik > NEG_INF) all.add(mc.log_lik);
        return all.value();
    }

    // How the scored fragments were routed.
    struct ScoreStats {
        int64_t deterministic = 0;                       // counted directly into the unambiguous counts
        int64_t em_unambig = 0, em_ambig_same = 0, em_ambig_opp = 0;   // EM units, by fragment class
        int64_t gated = 0;                               // no candidate survived a gate
        int64_t multimapper = 0;                         // multimapper EM units
    };

    // A multimapper's place: the group's hits whose extents overlap on one reference, chained by extent as
    // `merge_places` chains one key's hits. The best RNA term there in each pool.
    struct PlaceBest {
        double annotated = NEG_INF;
        double synthetic = NEG_INF;
        double rna() const { return std::max(annotated, synthetic); }
    };

    // One term's extent on its reference, for finding the places: `term` indexes the RNA terms, then the gDNA terms.
    struct TermSpan {
        int32_t ref, start, end;
        uint32_t term;
    };

    // The open multimapper group, which persists across chunks: alignments are scored as they arrive, and
    // flush_mm_group() prunes them within each place, merges them over places and emits.
    struct MultimapperGroup {
        int64_t fid = -1;
        int     n_members = 0;
        int     best_stype = SPLICE_UNSPLICED;
        std::vector<PlaceTerm> terms;
        std::vector<Candidate> merged;
        // gDNA's terms at the places it can explain. A place with a transcript candidate adds to the locus's term,
        // merged over places at the flush; a place with none adds to the intergenic term, an online log-sum-exp in
        // which each hit is its own place.
        std::vector<PlaceTerm> gdna_terms;
        std::vector<Candidate> gdna_merged;
        LogSumExp intergenic;
        // The flush's scratch: each term's place, each place's best, and the terms kept within their place
        std::vector<TermSpan> spans;
        std::vector<int32_t> place_of;
        std::vector<PlaceBest> places;
        std::vector<PlaceTerm> kept, gdna_kept;

        void reset() {
            terms.clear();  // every vector keeps its allocation
            merged.clear();
            gdna_terms.clear();
            gdna_merged.clear();
            fid = -1;
            n_members = 0;
            best_stype = SPLICE_UNSPLICED;
            intergenic = LogSumExp{};
        }
    };

    // The output of one scoring pass, grown chunk by chunk.
    struct FillState {
        // CSR candidate arrays, and their offsets from {0}
        std::vector<int32_t> ti;
        std::vector<float>   ll;
        std::vector<uint8_t> ct;
        std::vector<float>   cw;
        std::vector<int64_t> offsets{0};
        int64_t cand_cur = 0;
        // Per unit
        std::vector<int8_t>  is_spliced;
        std::vector<float>   gdna_ll, ig_ll;
        std::vector<int64_t> fid;
        std::vector<int8_t>  fclass;
        std::vector<uint8_t> stype;
        // Deterministic unambiguous assignments
        std::vector<int32_t> det_ti;
        std::vector<int64_t> det_fid;
        // The fragments a hard gate dropped, and their classes: no pool holds them, and the annotated BAM says so
        std::vector<int64_t> gated_fid;
        std::vector<int8_t>  gated_fclass;
        // Multimappers whose only places left by the gates are intergenic: intergenic gDNA, in no locus
        std::vector<int64_t> intergenic_fid;
        // Unambiguous counts, [n_transcripts, ua_ncols], the caller's array
        double* ua_counts = nullptr;
        int ua_ncols = 0;
        ScoreStats stats;
        MultimapperGroup mm;
    };

    // The count column of a fragment of splice type `stype` and exon strand `exon_str` on transcript `t_idx`: two per
    // splice type, sense then antisense. A fragment is antisense when its read-1 orientation disagrees with the
    // library's, and sense when it has no exon strand.
    inline int32_t count_col(int stype, int exon_str, int32_t t_idx) const {
        const bool stranded = exon_str == STRAND_POS || exon_str == STRAND_NEG;
        const bool antisense = stranded && ((exon_str == t_strand_[t_idx]) == r1_antisense_);
        return stype * 2 + (antisense ? 1 : 0);
    }

    // The term of row `row`'s candidate entry `k`: strand, the fragment length the transcript implies, and the
    // overhang past its exons. A term the whole hit shares is the caller's. False when the entry is no candidate: no
    // exonic base, or an overhang under the overhang gate.
    inline bool candidate_term(const ChunkPtrs& cp, int row, int64_t k, Candidate& c) const {
        const int32_t t_idx = cp.t_ind[k];
        const int32_t ebp = cp.e_bp[k];
        if (ebp <= 0) return false;
        const int rl = cp.r_len[row] > 0 ? static_cast<int>(cp.r_len[row]) : 1;
        const int32_t oh = std::max(rl - ebp, 0);
        if (oh_gate_ && oh > 0) return false;

        const int exon_str = cp.e_str[row];
        double log_strand = LOG_HALF;
        if (exon_str == STRAND_POS || exon_str == STRAND_NEG)
            log_strand = exon_str == t_strand_[t_idx] ? log_p_sense_ : log_p_antisense_;
        const int32_t flen = cp.f_len[k];
        const int32_t genomic_start = cp.g_sta[row];
        c.key = t_idx;
        c.log_lik = log_strand + rna_fl_.log_lik(flen) + oh * oh_log_pen_;
        c.count_col = count_col(cp.s_type[row], exon_str, t_idx);
        c.coverage_wt = (genomic_start >= 0 && flen > 0) ? coverage_weight(genomic_start, flen, t_idx) : 1.0;
        return true;
    }

    // gDNA's term at a hit, before any term the hit shares: the gDNA length law at the genomic footprint, and ½ for a
    // strand it does not prefer.
    inline double gdna_term(int32_t footprint) const {
        return gdna_fl_.log_lik(footprint) + LOG_HALF;
    }

    // Scores one multimapper alignment from the current chunk: appends its candidates to st.mm.terms and updates the
    // group's splice type and gDNA terms.
    void score_mm_alignment(
        const ChunkPtrs& cp, int row,
        FillState& st) const
    {
        int64_t start = cp.t_off[row];
        int64_t end   = cp.t_off[row + 1];
        int n_cand    = static_cast<int>(end - start);

        int stype    = cp.s_type[row];
        // The hit's sequence and its place weight apply to every candidate at it, RNA and gDNA alike.
        const int nm = cp.nm[row], clip = cp.clip[row];
        const double log_seq = (nm > 0 ? nm * mismatch_log_pen_ : 0.0) + (clip > 0 ? clip * clip_log_pen_ : 0.0);
        const double log_w = cp.lw ? static_cast<double>(cp.lw[row]) : 0.0;
        // The hit's genomic extent, which tells a component's places apart.
        int32_t gfp_val = cp.g_fp[row];
        const int32_t hit_start = cp.g_sta[row];
        const int32_t hit_end = hit_start + std::max(gfp_val, 1);

        // The group's reported splice type
        if (stype == SPLICE_SPLICED_ANNOT) st.mm.best_stype = SPLICE_SPLICED_ANNOT;

        // gDNA's term is summed over every place it can explain, since the molecule came from one of them. A place
        // with a transcript candidate is the locus's, on the candidates' reference; a place with none is intergenic.
        // A -inf term adds nothing and is skipped, so the hits' order cannot matter. The EM applies each component's
        // rate and length.
        if (gdna_can_explain(stype))
        {
            double hit_log_ll = gdna_term(gfp_val) + log_seq + log_w;
            if (hit_log_ll > NEG_INF) {
                if (n_cand > 0)
                    st.mm.gdna_terms.push_back({t_ref_[cp.t_ind[start]], hit_start, hit_end, hit_log_ll, 0, 0.0});
                else
                    st.mm.intergenic.add(hit_log_ll);
            }
        }

        // Read only as the group's non-empty test
        ++st.mm.n_members;

        Candidate c;
        for (int64_t k = start; k < end; ++k)
            if (candidate_term(cp, row, k, c))
                st.mm.terms.push_back({c.key, hit_start, hit_end, c.log_lik + log_seq + log_w, c.count_col,
                                       c.coverage_wt});
    }

    // Writes one candidate to the CSR.
    static void emit(const Candidate& c, FillState& st) {
        st.ti.push_back(c.key);
        st.ll.push_back(static_cast<float>(c.log_lik));
        st.ct.push_back(static_cast<uint8_t>(c.count_col));
        st.cw.push_back(static_cast<float>(c.coverage_wt));
        ++st.cand_cur;
    }

    // Writes a one-place unit's candidates to the CSR, pruned by pool: annotated transcripts and synthetic spans are
    // each kept within `max_ll_delta_` of their own pool's best. Returns the best kept term, -inf when none is kept.
    double emit_pruned(const Candidate* c, int n, FillState& st) const {
        double annotated_best = NEG_INF;
        double synthetic_best = NEG_INF;
        for (int j = 0; j < n; ++j) {
            double& pool_best = t_is_synthetic_[c[j].key] ? synthetic_best : annotated_best;
            pool_best = std::max(pool_best, c[j].log_lik);
        }
        double best_ll = NEG_INF;
        for (int j = 0; j < n; ++j) {
            const double pool_best = t_is_synthetic_[c[j].key] ? synthetic_best : annotated_best;
            if (pool_best - c[j].log_lik <= max_ll_delta_) {
                emit(c[j], st);
                best_ll = std::max(best_ll, c[j].log_lik);
            }
        }
        return best_ll;
    }

    // A multimapper's terms pruned within each place, never across places: at a place, each RNA term is kept within
    // `max_ll_delta_` of its pool's best there, and the gDNA term within it of the place's best RNA, or always where
    // the place has no RNA. Every place keeps its best, so none is dropped whole, and which place the molecule came
    // from is the EM's to decide. Fills `mm.kept` and `mm.gdna_kept`.
    void keep_within_places(MultimapperGroup& mm) const {
        const auto n_rna = static_cast<uint32_t>(mm.terms.size());
        mm.spans.clear();
        for (uint32_t i = 0; i < n_rna; ++i)
            mm.spans.push_back({t_ref_[mm.terms[i].key], mm.terms[i].start, mm.terms[i].end, i});
        for (uint32_t i = 0; i < mm.gdna_terms.size(); ++i)
            mm.spans.push_back({mm.gdna_terms[i].key, mm.gdna_terms[i].start, mm.gdna_terms[i].end, n_rna + i});
        std::sort(mm.spans.begin(), mm.spans.end(), [](const TermSpan& a, const TermSpan& b) {
            return a.ref != b.ref ? a.ref < b.ref : a.start < b.start;
        });
        mm.place_of.assign(mm.spans.size(), -1);
        mm.places.clear();
        int32_t ref = 0, end = 0;
        for (const TermSpan& s : mm.spans) {
            if (mm.places.empty() || s.ref != ref || s.start >= end) {
                mm.places.emplace_back();
                ref = s.ref;
                end = s.end;
            } else {
                end = std::max(end, s.end);
            }
            mm.place_of[s.term] = static_cast<int32_t>(mm.places.size()) - 1;
        }
        for (uint32_t i = 0; i < n_rna; ++i) {
            PlaceBest& best = mm.places[mm.place_of[i]];
            double& pool = t_is_synthetic_[mm.terms[i].key] ? best.synthetic : best.annotated;
            pool = std::max(pool, mm.terms[i].log_lik);
        }
        mm.kept.clear();
        for (uint32_t i = 0; i < n_rna; ++i) {
            const PlaceBest& best = mm.places[mm.place_of[i]];
            const double pool = t_is_synthetic_[mm.terms[i].key] ? best.synthetic : best.annotated;
            if (pool - mm.terms[i].log_lik <= max_ll_delta_) mm.kept.push_back(mm.terms[i]);
        }
        // A place with no RNA term reads -inf, which every gDNA term competes with.
        mm.gdna_kept.clear();
        for (uint32_t i = 0; i < mm.gdna_terms.size(); ++i)
            if (gdna_competes(mm.gdna_terms[i].log_lik, mm.places[mm.place_of[n_rna + i]].rna()))
                mm.gdna_kept.push_back(mm.gdna_terms[i]);
    }

    // Closes the unit whose candidates were just written: its CSR offset, identity, class and splice type, and its two
    // gDNA terms, -inf where gDNA is no candidate. "Spliced" means gDNA is no candidate anywhere: certified RNA.
    static void close_unit(FillState& st, int64_t fid, int fclass, int stype, bool has_gdna, double gdna_ll,
                           double ig_ll) {
        st.offsets.push_back(st.cand_cur);
        st.fid.push_back(fid);
        st.fclass.push_back(static_cast<int8_t>(fclass));
        st.stype.push_back(static_cast<uint8_t>(stype));
        st.is_spliced.push_back((has_gdna || ig_ll > NEG_INF) ? 0 : 1);
        st.gdna_ll.push_back(static_cast<float>(has_gdna ? gdna_ll : NEG_INF));
        st.ig_ll.push_back(static_cast<float>(ig_ll));
    }

    // A fragment no pool holds: a hard gate dropped every candidate.
    static void gate(FillState& st, int64_t fid, int fclass) {
        ++st.stats.gated;
        st.gated_fid.push_back(fid);
        st.gated_fclass.push_back(static_cast<int8_t>(fclass));
    }

    // Closes the open multimapper group: prunes its terms within each place (`keep_within_places`), merges each
    // transcript's kept terms over its places, decides the group's outcome and resets the accumulators. The outcome is
    //   - an EM unit, when an RNA candidate survives;
    //   - intergenic gDNA, when none does but an intergenic place does: no locus holds it (gDNA at a genic place has
    //     no locus without a candidate);
    //   - gated, when no place survives.
    void flush_mm_group(FillState& st) const {
        MultimapperGroup& mm = st.mm;
        if (mm.n_members == 0) return;
        keep_within_places(mm);
        merge_places(mm.kept, mm.merged);
        const int64_t emit_start = st.cand_cur;
        for (const Candidate& c : mm.merged) emit(c, st);
        // The group's two gDNA terms, each summed over its places: the locus's, over the places that kept it, and the
        // intergenic one; -inf where no place of the kind survived.
        const double gdna_ll = sum_over_places(mm.gdna_kept, mm.gdna_merged);
        const double ig_ll = mm.intergenic.value();
        if (st.cand_cur > emit_start) {
            // The intergenic term is kept whenever a place survived: its rate and the locus's rates and lengths are not
            // in these scores, so no difference of scores bounds its posterior.
            close_unit(st, mm.fid, FRAG_MULTIMAPPER, mm.best_stype, gdna_ll > NEG_INF, gdna_ll, ig_ll);
            ++st.stats.multimapper;
        } else if (ig_ll > NEG_INF) {
            st.intergenic_fid.push_back(mm.fid);
        } else {
            gate(st, mm.fid, FRAG_MULTIMAPPER);
        }
        mm.reset();
    }

    // A deterministic unambiguous fragment (one spliced annotated candidate), counted directly. The overhang gate
    // holds here too, so the fragment is gated exactly as an EM-routed one would be.
    void count_deterministic(const ChunkPtrs& cp, int i, FillState& st) const {
        const int64_t start = cp.t_off[i];
        const int32_t t_idx = cp.t_ind[start];
        const int rl = cp.r_len[i] > 0 ? static_cast<int>(cp.r_len[i]) : 1;
        if (oh_gate_ && rl - cp.e_bp[start] > 0) {
            gate(st, cp.f_id[i], cp.fc[i]);
            return;
        }
        const int32_t col = count_col(cp.s_type[i], cp.e_str[i], t_idx);
        if (col < st.ua_ncols) st.ua_counts[t_idx * st.ua_ncols + col] += 1.0;
        st.det_ti.push_back(t_idx);
        st.det_fid.push_back(cp.f_id[i]);
        ++st.stats.deterministic;
    }

    // An EM-routed fragment: one hit, so no sequence term (it would be common to every candidate). Every candidate is
    // scored, the pruned survivors written as the unit, and its one gDNA term kept when it competes with the unit's
    // best RNA candidate; the EM applies the locus's gDNA effective length.
    void score_em_unit(const ChunkPtrs& cp, int i, FillState& st) const {
        const int64_t start = cp.t_off[i], end = cp.t_off[i + 1];
        const int n_cand = static_cast<int>(end - start);
        Candidate scored_stack[SCORED_STACK_CAPACITY];
        std::vector<Candidate> scored_heap;
        Candidate* scored = scored_stack;
        if (n_cand > SCORED_STACK_CAPACITY) {
            scored_heap.resize(n_cand);
            scored = scored_heap.data();
        }
        int n_scored = 0;
        Candidate c;
        for (int64_t k = start; k < end; ++k)
            if (candidate_term(cp, i, k, c)) scored[n_scored++] = c;
        const int64_t emit_start = st.cand_cur;
        const double best_ll = emit_pruned(scored, n_scored, st);
        const int fclass = cp.fc[i], stype = cp.s_type[i];
        if (st.cand_cur == emit_start) {
            gate(st, cp.f_id[i], fclass);
            return;
        }
        double gdna_ll = 0.0;
        bool has_gdna = false;
        if (gdna_can_explain(stype)) {
            gdna_ll = gdna_term(cp.g_fp[i]);
            has_gdna = gdna_competes(gdna_ll, best_ll);
        }
        // One place: no intergenic reading beside a resolved hit.
        close_unit(st, cp.f_id[i], fclass, stype, has_gdna, gdna_ll, NEG_INF);
        if (fclass == FRAG_UNAMBIG) ++st.stats.em_unambig;
        else if (fclass == FRAG_AMBIG_SAME_STRAND) ++st.stats.em_ambig_same;
        else ++st.stats.em_ambig_opp;
    }

    // Scores one chunk's fragments into the output. A multimapper's alignments join its open group, which the next
    // fragment of another identity closes; a fragment with no candidate is skipped and not counted as gated.
    void score_chunk_impl(const ChunkPtrs& cp, FillState& st) const {
        for (int i = 0; i < cp.N; ++i) {
            if (cp.fc[i] == FRAG_MULTIMAPPER) {
                if (cp.f_id[i] != st.mm.fid) {
                    flush_mm_group(st);
                    st.mm.fid = cp.f_id[i];
                }
                score_mm_alignment(cp, i, st);
                continue;
            }
            flush_mm_group(st);
            if (cp.t_off[i + 1] - cp.t_off[i] <= 0) continue;
            if (cp.fc[i] == FRAG_UNAMBIG && cp.s_type[i] == SPLICE_SPLICED_ANNOT) count_deterministic(cp, i, st);
            else score_em_unit(cp, i, st);
        }
    }
};

// Scores one chunk per score_chunk() call, so Python can free each chunk's arrays as soon as it returns; a multimapper
// group holds no reference into an earlier chunk. finish() closes the last group and returns the CSR arrays. The
// scorer and the unambiguous counts are kept alive by the binding for as long as this object lives.
class StreamingScorer {
    const NativeFragmentScorer& scorer_;
    NativeFragmentScorer::FillState st_;
    bool finished_ = false;

public:
    StreamingScorer(NativeFragmentScorer& scorer, f64_2d_out unambig_counts) : scorer_(scorer) {
        st_.ua_counts = unambig_counts.data();
        st_.ua_ncols = static_cast<int>(unambig_counts.shape(1));
    }

    // One chunk's arrays, by name (`FinalizedChunk.to_scoring_arrays`); `log_weight`, each row's place weight, is
    // optional and reads 0.
    void score_chunk(nb::dict arrays) {
        if (finished_) throw std::runtime_error("StreamingScorer: score_chunk called after finish");
        NativeFragmentScorer::ChunkPtrs cp{};
        cp.t_off  = chunk_array<int32_t>(arrays, "t_offsets");
        cp.t_ind  = chunk_array<int32_t>(arrays, "t_indices");
        cp.f_len  = chunk_array<int32_t>(arrays, "frag_lengths");
        cp.e_bp   = chunk_array<uint16_t>(arrays, "exon_bp");
        cp.s_type = chunk_array<uint8_t>(arrays, "splice_type");
        cp.e_str  = chunk_array<uint8_t>(arrays, "align_strand");
        cp.fc     = chunk_array<uint8_t>(arrays, "fragment_classes");
        cp.f_id   = chunk_array<int64_t>(arrays, "frag_id");
        cp.r_len  = chunk_array<uint16_t>(arrays, "read_length");
        cp.g_fp   = chunk_array<int32_t>(arrays, "genomic_footprint");
        cp.g_sta  = chunk_array<int32_t>(arrays, "genomic_start");
        cp.nm     = chunk_array<uint16_t>(arrays, "nm");
        cp.clip   = chunk_array<uint16_t>(arrays, "clip");
        cp.lw     = arrays.contains("log_weight") ? chunk_array<float>(arrays, "log_weight") : nullptr;
        cp.N      = static_cast<int>(nb::cast<u8_1d>(arrays["splice_type"], false).shape(0));
        scorer_.score_chunk_impl(cp, st_);
    }

    // The CSR arrays, the per-unit arrays and the routing counts, handed to numpy without a copy.
    nb::tuple finish() {
        if (finished_) throw std::runtime_error("StreamingScorer: finish called twice");
        finished_ = true;
        scorer_.flush_mm_group(st_);
        NativeFragmentScorer::FillState& s = st_;
        const auto& n = s.stats;
        return nb::make_tuple(
            vec_to_ndarray(std::move(s.offsets)), vec_to_ndarray(std::move(s.ti)), vec_to_ndarray(std::move(s.ll)),
            vec_to_ndarray(std::move(s.ct)), vec_to_ndarray(std::move(s.cw)), vec_to_ndarray(std::move(s.is_spliced)),
            vec_to_ndarray(std::move(s.gdna_ll)), vec_to_ndarray(std::move(s.fid)), vec_to_ndarray(std::move(s.fclass)),
            vec_to_ndarray(std::move(s.stype)), vec_to_ndarray(std::move(s.det_ti)),
            vec_to_ndarray(std::move(s.det_fid)), vec_to_ndarray(std::move(s.gated_fid)),
            vec_to_ndarray(std::move(s.gated_fclass)), n.deterministic, n.em_unambig, n.em_ambig_same, n.em_ambig_opp,
            n.gated, n.multimapper, vec_to_ndarray(std::move(s.ig_ll)), vec_to_ndarray(std::move(s.intergenic_fid)));
    }
};

NB_MODULE(_scoring_impl, m) {
    m.doc() = "The fragment scorer: per-fragment candidate likelihoods, streamed chunk by chunk into the EM's units.";

    nb::class_<NativeFragmentScorer>(m, "NativeFragmentScorer")
        .def(nb::init<
                 double, double, bool, double, double, double,
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
             nb::arg("clip_log_penalty"),
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
             nb::arg("pruning_max_ll_delta"));

    nb::class_<StreamingScorer>(m, "StreamingScorer")
        .def(nb::init<NativeFragmentScorer&, f64_2d_out>(), nb::arg("scorer"), nb::arg("unambig_counts").noconvert(),
             nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>())
        .def("score_chunk", &StreamingScorer::score_chunk, nb::arg("arrays"))
        .def("finish", &StreamingScorer::finish);

    // The fragment classes this scorer routes on, the one definition rigel.buffer reads
    m.attr("FRAG_UNAMBIG")           = rigel::FRAG_UNAMBIG;
    m.attr("FRAG_AMBIG_SAME_STRAND") = rigel::FRAG_AMBIG_SAME_STRAND;
    m.attr("FRAG_AMBIG_OPP_STRAND")  = rigel::FRAG_AMBIG_OPP_STRAND;
    m.attr("FRAG_MULTIMAPPER")       = rigel::FRAG_MULTIMAPPER;

    rigel::bind_testing(m.def_submodule("testing", "The scorer's constants, bound for the parity tests."));
}
