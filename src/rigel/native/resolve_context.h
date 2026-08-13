/**
 * resolve_context.h — FragmentResolver, ResolvedFragment, FragmentAccumulator.
 *
 * Extracted from resolve.cpp so that bam_scanner.cpp can call
 * _resolve_core() directly in C++ and accumulate results into
 * FragmentAccumulator without crossing the Python boundary.
 *
 * Depends on constants.h and cgranges.h.
 */

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

extern "C" {
#include "cgranges.h"
}

#include "constants.h"
#include "ndarray_util.h"

namespace nb = nanobind;

namespace rigel {

static inline uint16_t checked_u16_buffer_value(
    int32_t value,
    const char* column,
    bool zero_if_missing = false)
{
    if (value < 0 && zero_if_missing) return 0;
    if (value < 0 || value > static_cast<int32_t>(std::numeric_limits<uint16_t>::max())) {
        throw std::runtime_error(
            std::string("Fragment buffer column '") + column
            + "' cannot store value " + std::to_string(value)
            + " as uint16");
    }
    return static_cast<uint16_t>(value);
}

// ================================================================
// ResolvedFragment — C++ result object exposed to Python
// ================================================================

class ResolvedFragment {
public:
    std::vector<int32_t> t_inds;
    int32_t ambig_strand = 0;
    int32_t splice_type = 0;
    int32_t align_strand = 0;
    int32_t sj_strand = 0;
    int32_t genomic_footprint = -1;
    int32_t genomic_start = -1;
    int32_t merge_criteria = MC_EMPTY;
    int32_t read_length = 0;
    int32_t chimera_type = CHIMERA_NONE;
    int32_t chimera_gap = -1;
    int32_t num_hits = 1;
    int32_t nm = 0;

    // Leftmost ANNOTATED junction (ref, start, end); -1 when the fragment
    // crosses none.  With `sj_strand` this is the per-junction SJ strand
    // table's key — see RawResolveResult in constants.h for why it is the
    // first ANNOTATED junction and why only one junction is credited.
    int32_t sj_key_ref = -1;
    int32_t sj_key_start = -1;
    int32_t sj_key_end = -1;

    // Per-transcript parallel arrays (same order as t_inds)
    std::vector<int32_t> frag_lengths;      // -1 if missing
    std::vector<int32_t> exon_bp;
    std::vector<int32_t> intron_bp;

    // SRD v2: strand-aware collapsed overlap counts
    int32_t exon_bp_pos = 0;
    int32_t exon_bp_neg = 0;
    int32_t tx_bp_pos = 0;
    int32_t tx_bp_neg = 0;

    // The annotated introns found inside this fragment's UNSEQUENCED gaps, and whether the candidate
    // transcripts disagreed about them. Carried out of `RawResolveResult` unchanged.
    //
    // ⚠ INSTRUMENTATION, and it is load-bearing for one gate that is otherwise unobservable. The
    // deposit adapter unions these with the observed CIGAR-N introns and the accumulator normalises
    // the union, so an intron the detector re-derived from an observed splice would be merged away
    // and leave `L` untouched — i.e. "it was never emitted" and "it was emitted and absorbed" are
    // indistinguishable downstream. The rule is precisely about the case where they
    // are NOT the same thing (a NEAR match shortens `L`), so the emission itself has to be visible.
    std::vector<IntronBlock> gap_introns;
    std::vector<int32_t>     gap_intron_offsets;
    int32_t                  n_gap_hypotheses = 0;

    // --- Properties for Python model training ---

    bool get_is_chimeric() const {
        return chimera_type != CHIMERA_NONE;
    }

    bool get_is_same_strand() const {
        return !ambig_strand;
    }

    bool get_is_strand_qualified() const {
        return splice_type == SPLICE_SPLICED_ANNOT
            && !ambig_strand
            && (align_strand == STRAND_POS || align_strand == STRAND_NEG)
            && (sj_strand == STRAND_POS || sj_strand == STRAND_NEG)
            && chimera_type == CHIMERA_NONE;
    }

    int32_t get_first_t_ind() const {
        return t_inds.empty() ? -1 : t_inds[0];
    }

    bool get_has_frag_lengths() const {
        for (auto fl : frag_lengths)
            if (fl > 0) return true;
        return false;
    }

    int32_t get_unique_frag_length() const {
        int32_t val = -1;
        for (auto fl : frag_lengths) {
            if (fl > 0) {
                if (val == -1) val = fl;
                else if (fl != val) return -1;
            }
        }
        return val;
    }

    // ⛔ `get_unique_frag_length_mrna` was DELETED by TRAPS: pure-and-length-censored.
    // It was definition **B**: a TRANSCRIPT-SPACE fragment length, gated on every non-nRNA candidate
    // agreeing, silently discarding the 4.6 % that did not. Summed into one array with definition
    // **A** — a GENOMIC footprint over a disjoint subset — and called the library's unconditional
    // fragment-length distribution, which it was neither unconditional nor one quantity.
    //
    // The tool now has ONE definition: the accumulator's `L`, the total length of the fragment's own
    // path (span minus region_bound introns, mate gap included), proven exhaustively in TRAPS: two-divisors-opposite-sign and binned for
    // every deposited fragment by TRAPS: a-purity-filter-is-a-length-filter's `deposited_lengths`.
    // Return t_inds as a Python frozenset for compatibility
    nb::object get_t_inds() const {
        nb::set s;
        for (int32_t t : t_inds) s.add(nb::cast(t));
        return nb::frozenset(s);
    }

    // Return frag_lengths as a Python dict for compatibility
    nb::dict get_frag_lengths() const {
        nb::dict d;
        for (size_t i = 0; i < t_inds.size(); i++) {
            if (frag_lengths[i] > 0)
                d[nb::cast(t_inds[i])] = nb::cast(frag_lengths[i]);
        }
        return d;
    }

    // Every hypothesis as a list of ``[(ref_id, start, end, strand), ...]``, in enumeration order.
    // ⭐ An EMPTY inner list is the unspliced (genomic) hypothesis.
    nb::object get_gap_hypotheses() const {
        nb::list out;
        for (int32_t h = 0; h < n_gap_hypotheses; ++h) {
            nb::list path;
            for (int32_t i = gap_intron_offsets[h]; i < gap_intron_offsets[h + 1]; ++i) {
                path.append(nb::make_tuple(gap_introns[i].ref_id, gap_introns[i].start,
                                           gap_introns[i].end, gap_introns[i].strand));
            }
            out.append(path);
        }
        return out;
    }

    // Return overlap_bp as a Python dict for compatibility
    nb::dict get_overlap_bp() const {
        nb::dict d;
        for (size_t i = 0; i < t_inds.size(); i++) {
            d[nb::cast(t_inds[i])] = nb::make_tuple(
                exon_bp[i], intron_bp[i]);
        }
        return d;
    }

    // Build from RawResolveResult (used by resolve_fragment and bam_scanner)
    static ResolvedFragment from_core(RawResolveResult& cr) {
        ResolvedFragment r;
        size_t n_t = cr.t_inds.size();
        r.t_inds = std::move(cr.t_inds);
        r.ambig_strand = cr.ambig_strand;
        r.splice_type = cr.splice_type;
        r.align_strand = cr.align_strand;
        r.sj_strand = cr.sj_strand;
        r.genomic_footprint = cr.genomic_footprint;
        r.genomic_start = cr.genomic_start;
        r.merge_criteria = cr.merge_criteria;
        r.read_length = cr.read_length;
        r.chimera_type = cr.chimera_type;
        r.chimera_gap = cr.chimera_gap;
        r.exon_bp = std::move(cr.t_exon_bp);
        r.intron_bp = std::move(cr.t_intron_bp);
        r.exon_bp_pos = cr.exon_bp_pos;
        r.exon_bp_neg = cr.exon_bp_neg;
        r.tx_bp_pos = cr.tx_bp_pos;
        r.tx_bp_neg = cr.tx_bp_neg;
        r.sj_key_ref = cr.sj_key_ref;
        r.sj_key_start = cr.sj_key_start;
        r.sj_key_end = cr.sj_key_end;
        // ⛔ COPIED, NOT MOVED, and the difference is a bug. `bam_scanner.cpp` builds the
        // ResolvedFragment first and calls `deposit_to_accumulator(frag, cr)` afterwards, so moving
        // this out of `cr` would hand the deposit adapter an empty intron list and silently stop
        // cutting every gap intron — the exact defect this field exists to make visible.
        r.gap_introns = cr.gap_introns;
        r.gap_intron_offsets = cr.gap_intron_offsets;
        r.n_gap_hypotheses = cr.n_gap_hypotheses();
        if (cr.frag_lengths.size() == n_t) {
            r.frag_lengths = std::move(cr.frag_lengths);
        } else {
            r.frag_lengths.assign(n_t, -1);
        }
        return r;
    }
};

// ================================================================
// FragmentAccumulator — C++ columnar buffer replacing _AccumulatorChunk
// ================================================================

class FragmentAccumulator {
public:
    std::vector<uint8_t>  splice_type_;
    std::vector<uint8_t>  align_strand_;
    std::vector<uint8_t>  sj_strand_;
    std::vector<uint8_t>  ambig_strand_;
    std::vector<uint16_t> num_hits_;
    std::vector<uint8_t>  merge_criteria_;
    std::vector<uint8_t>  chimera_type_;
    std::vector<int32_t>  t_indices_;
    std::vector<int32_t>  t_offsets_;
    std::vector<int32_t>  frag_lengths_;
    std::vector<uint16_t> exon_bp_;
    std::vector<int64_t>  frag_id_;
    std::vector<uint16_t> read_length_;
    std::vector<int32_t>  genomic_footprint_;
    std::vector<int32_t>  genomic_start_;
    std::vector<uint16_t> nm_;
    int32_t size_ = 0;

    FragmentAccumulator() {
        t_offsets_.push_back(0);
    }

    /// Pre-allocate internal vectors to avoid repeated reallocation.
    /// @param n_fragments  Expected number of fragments.
    /// @param n_candidates Expected total CSR entries (sum of per-fragment
    ///                     transcript candidate counts).
    void reserve(int64_t n_fragments, int64_t n_candidates) {
        splice_type_.reserve(n_fragments);
        align_strand_.reserve(n_fragments);
        sj_strand_.reserve(n_fragments);
        ambig_strand_.reserve(n_fragments);
        num_hits_.reserve(n_fragments);
        merge_criteria_.reserve(n_fragments);
        chimera_type_.reserve(n_fragments);
        frag_id_.reserve(n_fragments);
        read_length_.reserve(n_fragments);
        genomic_footprint_.reserve(n_fragments);
        genomic_start_.reserve(n_fragments);
        nm_.reserve(n_fragments);
        t_indices_.reserve(n_candidates);
        t_offsets_.reserve(n_fragments + 1);
        frag_lengths_.reserve(n_candidates);
        exon_bp_.reserve(n_candidates);
    }

    void append(const ResolvedFragment& r, int64_t frag_id) {
        splice_type_.push_back(static_cast<uint8_t>(r.splice_type));
        align_strand_.push_back(static_cast<uint8_t>(r.align_strand));
        sj_strand_.push_back(static_cast<uint8_t>(r.sj_strand));
        ambig_strand_.push_back(static_cast<uint8_t>(r.ambig_strand));
        num_hits_.push_back(static_cast<uint16_t>(r.num_hits));
        merge_criteria_.push_back(static_cast<uint8_t>(r.merge_criteria));
        chimera_type_.push_back(static_cast<uint8_t>(r.chimera_type));

        // CSR data — parallel arrays to t_indices
        t_indices_.insert(t_indices_.end(),
                          r.t_inds.begin(), r.t_inds.end());
        frag_lengths_.insert(frag_lengths_.end(),
                             r.frag_lengths.begin(), r.frag_lengths.end());
        for (int32_t bp : r.exon_bp) {
            exon_bp_.push_back(
                checked_u16_buffer_value(bp, "exon_bp"));
        }
        t_offsets_.push_back(static_cast<int32_t>(t_indices_.size()));

        frag_id_.push_back(frag_id);
        read_length_.push_back(
            checked_u16_buffer_value(r.read_length, "read_length"));
        genomic_footprint_.push_back(r.genomic_footprint);
        genomic_start_.push_back(r.genomic_start);
        nm_.push_back(static_cast<uint16_t>(r.nm));
        size_++;
    }

    int32_t get_size() const { return size_; }

    /// Finalize accumulator → dict of raw bytes (for numpy frombuffer).
    /// The transcript strand array is kept for API compatibility.
    nb::dict finalize(const std::vector<int32_t>& t_strand_arr) {
        (void)t_strand_arr;

        auto to_bytes = [](const void* data, size_t nbytes) -> nb::bytes {
            if (nbytes == 0)
                return nb::bytes(static_cast<const char*>(nullptr), 0);
            return nb::bytes(
                reinterpret_cast<const char*>(data), nbytes);
        };

        nb::dict result;
        result["splice_type"] = to_bytes(
            splice_type_.data(), splice_type_.size());
        result["align_strand"] = to_bytes(
            align_strand_.data(), align_strand_.size());
        result["sj_strand"] = to_bytes(
            sj_strand_.data(), sj_strand_.size());
        result["num_hits"] = to_bytes(
            num_hits_.data(), num_hits_.size() * sizeof(uint16_t));
        result["merge_criteria"] = to_bytes(
            merge_criteria_.data(), merge_criteria_.size());
        result["chimera_type"] = to_bytes(
            chimera_type_.data(), chimera_type_.size());
        result["t_indices"] = to_bytes(
            t_indices_.data(), t_indices_.size() * sizeof(int32_t));
        result["t_offsets"] = to_bytes(
            t_offsets_.data(), t_offsets_.size() * sizeof(int32_t));
        result["frag_lengths"] = to_bytes(
            frag_lengths_.data(), frag_lengths_.size() * sizeof(int32_t));
        result["exon_bp"] = to_bytes(
            exon_bp_.data(), exon_bp_.size() * sizeof(uint16_t));
        result["ambig_strand"] = to_bytes(ambig_strand_.data(), ambig_strand_.size());
        result["frag_id"] = to_bytes(
            frag_id_.data(), frag_id_.size() * sizeof(int64_t));
        result["read_length"] = to_bytes(
            read_length_.data(), read_length_.size() * sizeof(uint16_t));
        result["genomic_footprint"] = to_bytes(
            genomic_footprint_.data(),
            genomic_footprint_.size() * sizeof(int32_t));
        result["genomic_start"] = to_bytes(
            genomic_start_.data(),
            genomic_start_.size() * sizeof(int32_t));
        result["nm"] = to_bytes(
            nm_.data(), nm_.size() * sizeof(uint16_t));
        result["size"] = nb::cast(size_);

        return result;
    }

    /// Zero-copy finalize: moves internal vectors to heap-allocated
    /// storage and returns capsule-backed numpy arrays.  The accumulator
    /// is consumed (left empty) after this call.
    nb::dict finalize_zero_copy(const std::vector<int32_t>& t_strand_arr) {
        int32_t n = size_;
        (void)t_strand_arr;

        nb::dict result;
        result["splice_type"]       = vec_to_ndarray(std::move(splice_type_));
        result["align_strand"]       = vec_to_ndarray(std::move(align_strand_));
        result["sj_strand"]         = vec_to_ndarray(std::move(sj_strand_));
        result["num_hits"]          = vec_to_ndarray(std::move(num_hits_));
        result["merge_criteria"]    = vec_to_ndarray(std::move(merge_criteria_));
        result["chimera_type"]      = vec_to_ndarray(std::move(chimera_type_));
        result["t_indices"]         = vec_to_ndarray(std::move(t_indices_));
        result["t_offsets"]         = vec_to_ndarray(std::move(t_offsets_));
        result["frag_lengths"]      = vec_to_ndarray(std::move(frag_lengths_));
        result["exon_bp"]           = vec_to_ndarray(std::move(exon_bp_));
        result["ambig_strand"]      = vec_to_ndarray(std::move(ambig_strand_));
        result["frag_id"]           = vec_to_ndarray(std::move(frag_id_));
        result["read_length"]       = vec_to_ndarray(std::move(read_length_));
        result["genomic_footprint"] = vec_to_ndarray(std::move(genomic_footprint_));
        result["genomic_start"]     = vec_to_ndarray(std::move(genomic_start_));
        result["nm"]                = vec_to_ndarray(std::move(nm_));
        result["size"]              = nb::cast(n);

        size_ = 0;
        return result;
    }
};

// ================================================================
// ResolverScratch — per-thread mutable scratch buffers
// ================================================================
// Extracted from FragmentResolver so that multiple threads can call
// _resolve_core() and compute_frag_lengths() concurrently, each with
// its own scratch. The single-threaded path uses FragmentResolver's
// internal scratch_ member.

struct GapBlock {
    int32_t ref_id;
    int32_t start;
    int32_t end;
};

struct ResolverScratch {
    // Per-transcript bp counters (sparse-cleaned)
    std::vector<int32_t> t_exon_bp;
    std::vector<int32_t> t_transcript_bp;
    std::vector<uint8_t> t_dirty;
    std::vector<int32_t> dirty_indices;

    // cgranges query buffers (reusable per-call)
    int64_t* buf = nullptr;
    int64_t  buf_cap = 0;

    // Reusable scratch for SPLICED_IMPLICIT detection
    // (per-fragment block sort + PE-gap collection)
    std::vector<ExonBlock> implicit_blocks;
    std::vector<GapBlock> implicit_gaps;
    std::vector<IntronBlock> hypothesis_path;  // one candidate's path, while enumerating

    // Reusable sorted-vector set-operation buffers for _resolve_core.
    std::vector<std::vector<int32_t>> exon_t_sets;
    std::vector<std::vector<int32_t>> transcript_t_sets;
    std::vector<std::vector<int32_t>> sj_t_sets;
    std::vector<int32_t> tmp_a;
    std::vector<int32_t> tmp_b;
    std::vector<int32_t> tmp_out;
    std::vector<int32_t> tmp_union;
    std::vector<int32_t> all_overlap_t;
    std::vector<int32_t> nrna_t;

    ResolverScratch() = default;

    explicit ResolverScratch(int32_t n_transcripts)
        : t_exon_bp(n_transcripts, 0),
          t_transcript_bp(n_transcripts, 0),
          t_dirty(n_transcripts, 0)
    {
        dirty_indices.reserve(512);
    }

    ~ResolverScratch() {
        free(buf);
    }

    // Non-copyable, movable
    ResolverScratch(const ResolverScratch&) = delete;
    ResolverScratch& operator=(const ResolverScratch&) = delete;

    ResolverScratch(ResolverScratch&& o) noexcept
        : t_exon_bp(std::move(o.t_exon_bp)),
          t_transcript_bp(std::move(o.t_transcript_bp)),
          t_dirty(std::move(o.t_dirty)),
          dirty_indices(std::move(o.dirty_indices)),
          buf(o.buf),
          buf_cap(o.buf_cap),
          implicit_blocks(std::move(o.implicit_blocks)),
          implicit_gaps(std::move(o.implicit_gaps)),
          exon_t_sets(std::move(o.exon_t_sets)),
          transcript_t_sets(std::move(o.transcript_t_sets)),
          sj_t_sets(std::move(o.sj_t_sets)),
          tmp_a(std::move(o.tmp_a)),
          tmp_b(std::move(o.tmp_b)),
          tmp_out(std::move(o.tmp_out)),
          tmp_union(std::move(o.tmp_union)),
          all_overlap_t(std::move(o.all_overlap_t)),
          nrna_t(std::move(o.nrna_t))
    {
        o.buf = nullptr; o.buf_cap = 0;
    }

    ResolverScratch& operator=(ResolverScratch&& o) noexcept {
        if (this != &o) {
            free(buf);
            t_exon_bp = std::move(o.t_exon_bp);
            t_transcript_bp = std::move(o.t_transcript_bp);
            t_dirty = std::move(o.t_dirty);
            dirty_indices = std::move(o.dirty_indices);
            implicit_blocks = std::move(o.implicit_blocks);
            implicit_gaps = std::move(o.implicit_gaps);
            exon_t_sets = std::move(o.exon_t_sets);
            transcript_t_sets = std::move(o.transcript_t_sets);
            sj_t_sets = std::move(o.sj_t_sets);
            tmp_a = std::move(o.tmp_a);
            tmp_b = std::move(o.tmp_b);
            tmp_out = std::move(o.tmp_out);
            tmp_union = std::move(o.tmp_union);
            all_overlap_t = std::move(o.all_overlap_t);
            nrna_t = std::move(o.nrna_t);
            buf = o.buf; buf_cap = o.buf_cap;
            o.buf = nullptr; o.buf_cap = 0;
        }
        return *this;
    }

    void mark_dirty(int32_t t_idx) {
        if (!t_dirty[t_idx]) {
            t_dirty[t_idx] = 1;
            dirty_indices.push_back(t_idx);
        }
    }

    void clean() {
        for (int32_t t : dirty_indices) {
            t_exon_bp[t] = 0;
            t_transcript_bp[t] = 0;
            t_dirty[t] = 0;
        }
        dirty_indices.clear();
    }

    void reset_per_fragment() noexcept {
        for (int32_t t : dirty_indices) {
            t_exon_bp[t] = 0;
            t_transcript_bp[t] = 0;
            t_dirty[t] = 0;
        }
        dirty_indices.clear();

        for (auto& v : exon_t_sets) v.clear();
        for (auto& v : transcript_t_sets) v.clear();
        sj_t_sets.clear();
        tmp_a.clear();
        tmp_b.clear();
        tmp_out.clear();
        tmp_union.clear();
        all_overlap_t.clear();
        nrna_t.clear();
        implicit_blocks.clear();
        implicit_gaps.clear();
    }
};

// ================================================================
// FragmentResolver — holds all index data for C++ resolution
// ================================================================

class FragmentResolver {
public:
    // --- Overlap index (main cgranges) ---
    cgranges_t* cr_ = nullptr;
    std::vector<int8_t>  iv_type_;
    std::vector<int32_t> t_set_data_;
    std::vector<int32_t> t_set_offsets_;

    // --- SJ exact-match map ---
    using SJMap = std::unordered_map<SJKey, std::pair<int32_t, int32_t>, SJKeyHash>;
    SJMap sj_map_;
    std::vector<int32_t> sj_map_data_;

    // --- SJ artifact blacklist (from alignable) ---
    // Strand-agnostic: keyed on (ref_id, start, end).  A CIGAR splice
    // junction is rejected at BAM-scan time when either its left or
    // right anchor is <= the stored maximum.
    using SJBlacklistMap = std::unordered_map<
        SJBlacklistKey, SJBlacklistEntry, SJBlacklistKeyHash>;
    SJBlacklistMap sj_blacklist_;



    // --- Transcript metadata ---
    std::vector<int32_t> t_to_g_arr_;
    int32_t n_transcripts_ = 0;

    // --- Transcript strand metadata (direct lookup, no gene indirection) ---
    std::vector<int32_t> t_strand_arr_;

    // --- Per-transcript nRNA status (1 = nRNA synthetic, 0 = annotated) ---
    std::vector<uint8_t> t_is_nrna_;

    // --- Per-transcript nRNA parent index (-1 if none, else index of
    //     the parent nRNA entity for derive-on-demand resolution).
    //     Synthetic nRNAs are not in cgranges; they are surfaced into
    //     each fragment's candidate set during _resolve_core by
    //     mapping each real-tx hit to its nrna_parent_ entry.
    std::vector<int32_t> nrna_parent_;

    // --- Per-transcript exon structure (CSR layout) for FL computation ---
    std::vector<int32_t> exon_offsets_;    // [n_transcripts + 1]
    std::vector<int32_t> exon_starts_;     // [total_exons] — genomic start coords
    std::vector<int32_t> exon_ends_;       // [total_exons] — genomic end coords
    std::vector<int32_t> exon_cumsum_;     // [total_exons] — cumulative spliced bp before each exon
    std::vector<int32_t> t_length_;        // [n_transcripts] — spliced transcript lengths

    // --- SPLICED_IMPLICIT discriminant tolerance (bp) ---
    int32_t splicing_anchor_tolerance_ = 0;

    // --- Reference name <-> ID ---
    std::unordered_map<std::string, int32_t> ref_to_id_;
    std::vector<std::string> id_to_ref_;

    // --- Per-call scratch (mutable, NOT thread-safe) ---
    // Single-threaded callers use this. Multi-threaded callers pass
    // their own ResolverScratch to the overloaded methods.
    ResolverScratch scratch_;

    // ----------------------------------------------------------------
    FragmentResolver() = default;

    ~FragmentResolver() {
        if (cr_) cr_destroy(cr_);
    }

    // Non-copyable
    FragmentResolver(const FragmentResolver&) = delete;
    FragmentResolver& operator=(const FragmentResolver&) = delete;

    // ----------------------------------------------------------------
    // Ref-ID helpers
    // ----------------------------------------------------------------

    int32_t get_or_create_ref_id(const std::string& ref) {
        auto it = ref_to_id_.find(ref);
        if (it != ref_to_id_.end()) return it->second;
        int32_t id = static_cast<int32_t>(id_to_ref_.size());
        ref_to_id_[ref] = id;
        id_to_ref_.push_back(ref);
        return id;
    }

    // Seed the reference-ID space in a CANONICAL order (index.ref_names) so the
    // resolver's ref ids match index.ref_name_to_id — the same space used by
    // t_df, region_df, RegionArrays and the calibration accumulator partition.
    // Must be called BEFORE build_overlap_index: otherwise ref ids are assigned
    // by first-seen interval order (gene/transcript order), which scrambles them
    // relative to ref_names and silently mis-routes the accumulator deposit on
    // multi-reference genomes (a single-ref synthetic index hides the bug).
    void set_ref_names(const std::vector<std::string>& ref_names) {
        ref_to_id_.clear();
        id_to_ref_.clear();
        ref_to_id_.reserve(ref_names.size());
        id_to_ref_.reserve(ref_names.size());
        for (const auto& name : ref_names) get_or_create_ref_id(name);
    }

    // ----------------------------------------------------------------
    // Build methods — called once at index load time
    // ----------------------------------------------------------------

    void build_overlap_index(
        const std::vector<std::string>& refs,
        const std::vector<int32_t>& starts,
        const std::vector<int32_t>& ends,
        const std::vector<int32_t>& iv_types,
        const std::vector<int32_t>& tset_data,
        const std::vector<int32_t>& tset_offsets)
    {
        size_t n = refs.size();
        for (const auto& ref : refs) get_or_create_ref_id(ref);

        cr_ = cr_init();
        for (size_t i = 0; i < n; i++)
            cr_add(cr_, refs[i].c_str(), starts[i], ends[i],
                   static_cast<int32_t>(i));
        cr_index(cr_);

        iv_type_.resize(n);
        for (size_t i = 0; i < n; i++)
            iv_type_[i] = static_cast<int8_t>(iv_types[i]);

        t_set_data_ = tset_data;
        t_set_offsets_ = tset_offsets;
    }

    void build_sj_map(
        const std::vector<std::string>& refs,
        const std::vector<int32_t>& starts,
        const std::vector<int32_t>& ends,
        const std::vector<int32_t>& strands,
        const std::vector<int32_t>& tset_data,
        const std::vector<int32_t>& tset_offsets)
    {
        sj_map_data_ = tset_data;
        size_t n = refs.size();
        for (size_t i = 0; i < n; i++) {
            int32_t ref_id = get_or_create_ref_id(refs[i]);
            SJKey key{ref_id, starts[i], ends[i], strands[i]};
            sj_map_[key] = {tset_offsets[i],
                            tset_offsets[i + 1] - tset_offsets[i]};
        }
    }

    /// Build the splice-junction artifact blacklist (from alignable).
    ///
    /// Strand-agnostic: a single entry per (ref, start, end).  At BAM
    /// scan time ``sj_blacklist_lookup`` is used to test a CIGAR-derived
    /// junction; rejection uses the OR rule (either anchor <= max).
    void build_sj_blacklist_map(
        const std::vector<std::string>& refs,
        const std::vector<int32_t>& starts,
        const std::vector<int32_t>& ends,
        const std::vector<int32_t>& max_anchor_left,
        const std::vector<int32_t>& max_anchor_right)
    {
        size_t n = refs.size();
        sj_blacklist_.clear();
        sj_blacklist_.reserve(n);
        for (size_t i = 0; i < n; i++) {
            int32_t ref_id = get_or_create_ref_id(refs[i]);
            SJBlacklistKey key{ref_id, starts[i], ends[i]};
            sj_blacklist_[key] = {max_anchor_left[i], max_anchor_right[i]};
        }
    }

    bool has_sj_blacklist() const { return !sj_blacklist_.empty(); }

    /// Return blacklist entry for a junction, or nullptr if not present.
    inline const SJBlacklistEntry* sj_blacklist_lookup(
        int32_t ref_id, int32_t start, int32_t end) const
    {
        auto it = sj_blacklist_.find(SJBlacklistKey{ref_id, start, end});
        return (it == sj_blacklist_.end()) ? nullptr : &it->second;
    }

    void set_metadata(const std::vector<int32_t>& t_to_g,
                      int32_t n_transcripts) {
        t_to_g_arr_ = t_to_g;
        n_transcripts_ = n_transcripts;

        // Initialize internal scratch buffers
        scratch_ = ResolverScratch(n_transcripts_);
    }

    /// Set gene strand mapping (for BAM scanner model training)

    /// Set per-transcript strand array (direct lookup, no gene indirection).
    void set_transcript_strands(const std::vector<int32_t>& t_strand) {
        t_strand_arr_ = t_strand;
    }

    /// Set per-transcript nRNA status (uint8, 1 = nRNA synthetic).
    /// ⚠ On a NON-SYNTHETIC row `is_nrna` means "single-exon, so mature == nascent" — NOT
    /// "manufactured span". Read `nrna_mask()`'s own comment before using it as a realness filter.
    void set_nrna_status(const std::vector<uint8_t>& t_is_nrna) {
        t_is_nrna_ = t_is_nrna;
    }

    /// Set per-transcript nRNA parent index (int32; -1 = no parent).
    /// Used by _resolve_core to derive synthetic-nRNA candidates from
    /// real-transcript hits without bloating the cgranges index.
    void set_nrna_parent_index(const std::vector<int32_t>& nrna_parent) {
        nrna_parent_ = nrna_parent;
    }

    const uint8_t* nrna_mask() const {
        return t_is_nrna_.empty() ? nullptr : t_is_nrna_.data();
    }

    /// Build per-transcript exon CSR index for FL computation.
    void build_exon_index(
        const std::vector<int32_t>& offsets,
        const std::vector<int32_t>& starts,
        const std::vector<int32_t>& ends,
        const std::vector<int32_t>& cumsum,
        const std::vector<int32_t>& lengths)
    {
        exon_offsets_ = offsets;
        exon_starts_ = starts;
        exon_ends_ = ends;
        exon_cumsum_ = cumsum;
        t_length_ = lengths;
    }

    bool has_exon_index() const {
        return !exon_offsets_.empty();
    }

    // ----------------------------------------------------------------
    // SPLICED_IMPLICIT detection
    // ----------------------------------------------------------------

    /// Set the splicing-anchor tolerance K (bp) used for one-sided slack
    /// in the SPLICED_IMPLICIT per-intron whole-containment discriminant.
    /// K must be >= 0; default 0 (strict containment).
    void set_splicing_anchor_tolerance(int32_t K) {
        if (K < 0) K = 0;
        splicing_anchor_tolerance_ = K;
    }

    int32_t splicing_anchor_tolerance() const {
        return splicing_anchor_tolerance_;
    }

    /// Does transcript ``t`` contradict the alignment? True when one of its introns overlaps the
    /// INTERIOR of an aligned block -- the read has bases there, so the molecule cannot have spliced them
    /// out, and ``t`` is not a candidate explanation for anything.
    ///
    /// ⛔ THE CANDIDATE SET IS NOT ALREADY THIS. `cr.t_inds` comes from `merge_sets`, which falls back to
    /// a UNION when the intersection is empty, so a transcript the reads contradict can be in it. Without
    /// this test the enumeration would emit hypotheses no molecule has.
    ///
    /// ⚠ Abutting is fine: an intron may START where a block ends. And the overlap must exceed the
    /// SPLICING ANCHOR TOLERANCE K before it counts -- the same slack the gap match uses, because it is
    /// the same physical thing: an aligner routinely runs a read a base or two past a donor. Without
    /// this the predicate rejects exactly the fragments K exists to accept.
    inline bool transcript_contradicts_blocks(int32_t t,
                                              const std::vector<ExonBlock>& blocks,
                                              int32_t K) const
    {
        if (t < 0 || t + 1 >= static_cast<int32_t>(exon_offsets_.size())) return false;
        const int32_t begin = exon_offsets_[t];
        const int32_t n_introns = exon_offsets_[t + 1] - begin - 1;
        if (n_introns <= 0) return false;
        const int32_t* intron_starts = exon_ends_.data() + begin;
        const int32_t* intron_ends   = exon_starts_.data() + begin + 1;

        for (const ExonBlock& block : blocks) {
            // First intron that could end inside this block, then walk while it still starts before it.
            int32_t i = static_cast<int32_t>(
                std::upper_bound(intron_ends, intron_ends + n_introns, block.start) - intron_ends);
            for (; i < n_introns && intron_starts[i] < block.end; ++i) {
                const int32_t overlap = std::min(intron_ends[i], block.end) -
                                        std::max(intron_starts[i], block.start);
                if (overlap > K) return true;   // the read has bases the transcript splices out
            }
        }
        return false;
    }

    /// ⭐ EVERY explanation of this fragment's UNSEQUENCED gaps, grouped so that identical paths are ONE
    /// hypothesis. The accumulator arbitrates between them; this only enumerates.
    ///
    /// A hypothesis is a PATH THROUGH THE ANNOTATION, not "an intron": a gap may hold several introns and
    /// several whole exons that no read ever touched. Each compatible transcript determines exactly one
    /// path -- its own introns lying inside the gaps -- so the set is finite and small.
    ///
    /// ⛔ THE GAPS ARE ONLY THE UNSEQUENCED ONES. The walk over consecutive blocks emits every hole, and a
    /// read `50M200N50M` has two blocks, so a sequenced intron arrives indistinguishable from a mate gap.
    /// Searching it again is not merely redundant: the ±K anchor tolerance lets a DIFFERENT nearby intron
    /// answer for it, the two normalise into one wider interval, and L comes out too SHORT. Gaps matching
    /// an observed intron by EXACT (start, end) equality are dropped -- never by overlap, which would also
    /// drop a genuine gap intron abutting an observed one.
    ///
    /// ⚠ The path is compared WHOLE, across every gap, never per gap. A per-gap union can take gap A's
    /// intron from T1 and gap B's from T2 and emit a path no single molecule has; grouping by the whole
    /// path makes that unrepresentable rather than merely avoided.
    void enumerate_gap_hypotheses(const std::vector<ExonBlock>& exons,
                                  const std::vector<IntronBlock>& observed_introns,
                                  const std::vector<int32_t>& candidate_t,
                                  bool certified_rna,
                                  ResolverScratch& scratch,
                                  RawResolveResult& cr) const
    {
        cr.gap_introns.clear();
        cr.gap_sj_strand.clear();
        cr.gap_supporting.clear();
        cr.gap_intron_offsets.assign(1, 0);
        cr.gap_supporting_offsets.assign(1, 0);

        auto& gaps = scratch.implicit_gaps;
        gaps.clear();
        if (exons.size() >= 2 && has_exon_index()) {
            auto exon_block_less = [](const ExonBlock& a, const ExonBlock& b) {
                if (a.ref_id != b.ref_id) return a.ref_id < b.ref_id;
                if (a.start != b.start) return a.start < b.start;
                if (a.end != b.end) return a.end < b.end;
                return a.strand < b.strand;
            };
            const std::vector<ExonBlock>* ordered = &exons;
            if (!std::is_sorted(exons.begin(), exons.end(), exon_block_less)) {
                auto& blocks = scratch.implicit_blocks;
                blocks.assign(exons.begin(), exons.end());
                std::sort(blocks.begin(), blocks.end(), exon_block_less);
                ordered = &blocks;
            }
            int32_t cur_ref = (*ordered)[0].ref_id;
            int32_t cur_end = (*ordered)[0].end;
            for (size_t k = 1; k < ordered->size(); ++k) {
                const auto& block = (*ordered)[k];
                if (block.ref_id != cur_ref) {
                    cur_ref = block.ref_id;
                    cur_end = block.end;
                } else if (block.start > cur_end) {
                    gaps.push_back({cur_ref, cur_end, block.start});
                    cur_end = block.end;
                } else {
                    cur_end = std::max(cur_end, block.end);
                }
            }
            gaps.erase(std::remove_if(gaps.begin(), gaps.end(),
                                      [&observed_introns](const GapBlock& gap) {
                                          for (const auto& observed : observed_introns) {
                                              if (observed.ref_id == gap.ref_id &&
                                                  observed.start == gap.start &&
                                                  observed.end == gap.end) {
                                                  return true;
                                              }
                                          }
                                          return false;
                                      }),
                       gaps.end());
        }

        // ⭐ No unsequenced gap ⇒ one hypothesis, the unspliced one, and nothing to arbitrate. The
        // degenerate case is the general case: `deposit` still receives a set, of size one.
        const int32_t K = splicing_anchor_tolerance_;
        if (gaps.empty()) {
            emit_unspliced_hypothesis(cr);
            return;
        }

        auto& path = scratch.hypothesis_path;
        bool any_candidate_implies_nothing = false;

        for (int32_t t : candidate_t) {
            if (transcript_contradicts_blocks(t, exons, K)) continue;
            path.clear();
            const int32_t strand = (t >= 0 && t < static_cast<int32_t>(t_strand_arr_.size()))
                                       ? t_strand_arr_[t]
                                       : STRAND_NONE;
            for (const GapBlock& gap : gaps)
                collect_transcript_introns_in_gap(t, gap, K, strand, path);
            if (path.empty()) {
                any_candidate_implies_nothing = true;
                continue;
            }
            add_hypothesis(cr, path, strand, t);
        }

        // ⭐ The UNSPLICED hypothesis. Present when some compatible transcript implies nothing in the
        // gaps (a retained-intron isoform), and ALWAYS when no annotated junction was sequenced -- then
        // the molecule may be gDNA or nascent, and the gap is real template.
        if (any_candidate_implies_nothing || !certified_rna || cr.n_gap_hypotheses() == 0)
            emit_unspliced_hypothesis(cr);
    }

private:
    /// Append every intron of ``t`` lying inside ``gap`` (±K), in genomic order.
    ///
    /// ⛔ EVERY one, not the first. Returning the first and stopping is what made a mate gap spanning two
    /// annotated introns keep only one region_bound -- measured at 98.5 % of the fragment-length tail that survived
    /// the earlier single-region_bound form -- and it also broke the ambiguity test, because two transcripts differing only in their
    /// SECOND intron read as agreeing.
    inline void collect_transcript_introns_in_gap(int32_t t,
                                                  const GapBlock& gap,
                                                  int32_t K,
                                                  int32_t strand,
                                                  std::vector<IntronBlock>& out) const
    {
        if (t < 0 || t + 1 >= static_cast<int32_t>(exon_offsets_.size())) return;
        const int32_t begin = exon_offsets_[t];
        const int32_t n_introns = exon_offsets_[t + 1] - begin - 1;
        if (n_introns <= 0 || gap.end <= gap.start) return;

        const int32_t* intron_starts = exon_ends_.data() + begin;
        const int32_t* intron_ends   = exon_starts_.data() + begin + 1;
        const int64_t  max_end       = static_cast<int64_t>(gap.end) + K;

        int32_t i = static_cast<int32_t>(
            std::lower_bound(intron_starts, intron_starts + n_introns,
                             static_cast<int64_t>(gap.start) - K) - intron_starts);
        const int32_t by_end = static_cast<int32_t>(
            std::upper_bound(intron_ends, intron_ends + n_introns, gap.start) - intron_ends);
        for (i = std::max(i, by_end); i < n_introns; ++i) {
            if (intron_starts[i] >= gap.end) break;                       // no overlap with the gap
            if (static_cast<int64_t>(intron_ends[i]) > max_end) break;    // protrudes past the slack
            if (intron_ends[i] > intron_starts[i])
                out.push_back({gap.ref_id, intron_starts[i], intron_ends[i], strand});
        }
    }

    /// Record ``path`` as a hypothesis, or credit ``t`` to the identical one already there.
    static void add_hypothesis(RawResolveResult& cr,
                               const std::vector<IntronBlock>& path,
                               int32_t strand,
                               int32_t t)
    {
        for (int32_t h = 0; h < cr.n_gap_hypotheses(); ++h) {
            const int32_t lo = cr.gap_intron_offsets[h], hi = cr.gap_intron_offsets[h + 1];
            if (hi - lo != static_cast<int32_t>(path.size())) continue;
            bool same = true;
            for (int32_t i = 0; i < hi - lo && same; ++i) {
                same = cr.gap_introns[lo + i].start == path[i].start &&
                       cr.gap_introns[lo + i].end == path[i].end;
            }
            if (!same) continue;
            // ⛔ ONE PATH, BUT ITS SUPPORTERS MAY DISAGREE ABOUT THE STRAND. Grouping on coordinates is
            // right -- two transcripts implying the same introns imply the same molecule -- but the
            // hypothesis carries ONE `sj_strand`, and keeping the first supporter's would make the answer
            // depend on the order the transcripts happen to sit in the GTF. AMBIGUOUS is what this state
            // is called everywhere else (the fragment-level `sj_strand` uses it for exactly this:
            // contradictory evidence, not missing evidence), and `deposit` already refuses to credit a
            // junction on it. ⚠ Idempotent: once AMBIGUOUS, any further disagreement keeps it there.
            //
            // ⚠ Unreachable on human data -- 0 of 404,168 junction coordinates are annotated on both
            // strands, and the index warns that it is biologically impossible. Fixed because the
            // alternative is a silent answer that flips when two GTF boundaries are swapped.
            if (cr.gap_sj_strand[h] != strand) cr.gap_sj_strand[h] = STRAND_AMBIGUOUS;
            // ⚠ Inserted in place so the supporting lists stay contiguous per hypothesis.
            cr.gap_supporting.insert(cr.gap_supporting.begin() + cr.gap_supporting_offsets[h + 1], t);
            for (int32_t k = h + 1; k < static_cast<int32_t>(cr.gap_supporting_offsets.size()); ++k)
                ++cr.gap_supporting_offsets[k];
            return;
        }
        cr.gap_introns.insert(cr.gap_introns.end(), path.begin(), path.end());
        cr.gap_intron_offsets.push_back(static_cast<int32_t>(cr.gap_introns.size()));
        cr.gap_sj_strand.push_back(strand);
        cr.gap_supporting.push_back(t);
        cr.gap_supporting_offsets.push_back(static_cast<int32_t>(cr.gap_supporting.size()));
    }

    /// The empty path: region_bound nothing. ⚠ Idempotent, because two routes can both call for it.
    static void emit_unspliced_hypothesis(RawResolveResult& cr) {
        for (int32_t h = 0; h < cr.n_gap_hypotheses(); ++h) {
            if (cr.gap_intron_offsets[h] == cr.gap_intron_offsets[h + 1]) return;
        }
        cr.gap_intron_offsets.push_back(static_cast<int32_t>(cr.gap_introns.size()));
        cr.gap_sj_strand.push_back(STRAND_NONE);
        cr.gap_supporting_offsets.push_back(static_cast<int32_t>(cr.gap_supporting.size()));
    }

public:

    /// Map a genomic position to transcript-space offset for FL computation.
    ///
    /// Designed for projecting fragment ENDPOINTS (gstart, gend) only.
    /// When a position falls outside an exon, the overhang distance is added
    /// rather than snapped/clipped, because fragment endpoints represent the
    /// physical extent of the molecule and those bases must be counted in FL.
    ///
    /// - Strand-agnostic (no strand flip): always returns forward-strand
    ///   spliced coordinate.
    /// - NOT clamped: returns negative values (before first exon), values
    ///   > t_len (past last exon), or tx-position + intronic-overhang
    ///   (endpoint in internal intron).
    inline int32_t genomic_to_tx_pos(int32_t genomic_pos, int32_t t_idx) const {
        int32_t begin   = exon_offsets_[t_idx];
        int32_t end     = exon_offsets_[t_idx + 1];
        int32_t n_exons = end - begin;
        if (n_exons <= 0) return 0;

        const int32_t* starts = exon_starts_.data() + begin;
        const int32_t* ends   = exon_ends_.data()   + begin;
        const int32_t* cumsum = exon_cumsum_.data()  + begin;

        // bisect_right(starts, genomic_pos) - 1
        int ei = static_cast<int>(
            std::upper_bound(starts, starts + n_exons, genomic_pos) - starts) - 1;

        if (ei < 0) {
            // Before first exon: negative offset
            return genomic_pos - starts[0];
        }
        if (genomic_pos >= ends[ei]) {
            // Past end of exon[ei] — either in intron or past last exon.
            // Always add the overhang: endpoint bases must be counted in FL.
            return cumsum[ei] + (ends[ei] - starts[ei]) + (genomic_pos - ends[ei]);
        }
        // Inside exon ei
        return cumsum[ei] + (genomic_pos - starts[ei]);
    }

    nb::dict get_ref_to_id() const {
        nb::dict d;
        for (const auto& [k, v] : ref_to_id_)
            d[nb::cast(k)] = nb::cast(v);
        return d;
    }

    // ----------------------------------------------------------------
    // Fragment-length computation via transcript-space projection
    // ----------------------------------------------------------------

    /// Thread-safe fragment-length projection.  Results are aligned to
    /// t_inds: frag_lengths[i] is the length for t_inds[i], or -1.
    void compute_frag_lengths_aligned(
        const std::vector<ExonBlock>& exons,
        const std::vector<IntronBlock>& /*introns*/,
        const std::vector<int32_t>& t_inds,
        std::vector<int32_t>& frag_lengths,
        ResolverScratch& /*scratch*/) const
    {
        frag_lengths.assign(t_inds.size(), -1);
        if (exons.empty() || t_inds.empty()) return;

        // Single exon block: FL is trivially block length (same for all candidates)
        if (exons.size() == 1) {
            int32_t fl = exons[0].end - exons[0].start;
            if (fl > 0) std::fill(frag_lengths.begin(), frag_lengths.end(), fl);
            return;
        }

        // Find outermost genomic positions across all alignment blocks
        int32_t gstart = exons[0].start;
        int32_t gend = exons[0].end;
        for (const auto& e : exons) {
            gstart = std::min(gstart, e.start);
            gend = std::max(gend, e.end);
        }

        for (size_t i = 0; i < t_inds.size(); i++) {
            int32_t t = t_inds[i];
            int32_t tx_s = genomic_to_tx_pos(gstart, t);
            int32_t tx_e = genomic_to_tx_pos(gend, t);
            int32_t fl = std::abs(tx_e - tx_s);
            if (fl > 0) frag_lengths[i] = fl;
        }
    }

    /// Backward-compatible wrapper using internal scratch.
    std::vector<int32_t> compute_frag_lengths(
        const std::vector<ExonBlock>& exons,
        const std::vector<IntronBlock>& introns,
        const std::vector<int32_t>& t_inds)
    {
        std::vector<int32_t> frag_lengths;
        compute_frag_lengths_aligned(exons, introns, t_inds,
                                     frag_lengths, scratch_);
        return frag_lengths;
    }

    // ----------------------------------------------------------------
    // SJ map lookup helper
    // ----------------------------------------------------------------

    void sj_lookup_into(int32_t ref_id, int32_t start,
                        int32_t end, int32_t strand,
                        std::vector<int32_t>& out) const {
        out.clear();
        SJKey key{ref_id, start, end, strand};
        auto it = sj_map_.find(key);
        if (it != sj_map_.end()) {
            auto [off, cnt] = it->second;
            out.insert(out.end(), sj_map_data_.begin() + off,
                       sj_map_data_.begin() + off + cnt);
            sort_unique(out);
            return;
        }
        // Strand-agnostic fallback
        if (strand != STRAND_POS && strand != STRAND_NEG) {
            SJKey pk{ref_id, start, end, STRAND_POS};
            SJKey nk{ref_id, start, end, STRAND_NEG};
            auto pi = sj_map_.find(pk);
            auto ni = sj_map_.find(nk);
            if (pi != sj_map_.end() || ni != sj_map_.end()) {
                if (pi != sj_map_.end()) {
                    auto [off, cnt] = pi->second;
                    for (int32_t k = 0; k < cnt; k++)
                        out.push_back(sj_map_data_[off + k]);
                }
                if (ni != sj_map_.end()) {
                    auto [off, cnt] = ni->second;
                    for (int32_t k = 0; k < cnt; k++)
                        out.push_back(sj_map_data_[off + k]);
                }
                sort_unique(out);
                return;
            }
        }
    }

    std::vector<int32_t> sj_lookup(int32_t ref_id, int32_t start,
                                   int32_t end, int32_t strand) const {
        std::vector<int32_t> out;
        sj_lookup_into(ref_id, start, end, strand, out);
        return out;
    }

    // ----------------------------------------------------------------
    // Scratch buffer management (delegates to scratch_)
    // ----------------------------------------------------------------

    void mark_dirty(int32_t t_idx) {
        scratch_.mark_dirty(t_idx);
    }

    void clean_scratch() {
        scratch_.clean();
    }

    // ================================================================
    // _resolve_core — shared resolution logic behind resolve_fragment()
    // ================================================================

    /// Thread-safe overload: uses caller-supplied scratch buffers.
    bool _resolve_core(
        const std::vector<ExonBlock>& exons,
        const std::vector<IntronBlock>& introns,
        int32_t genomic_footprint,
        RawResolveResult& cr,
        ResolverScratch& scratch)
    {
        int n_exons = static_cast<int>(exons.size());
        if (n_exons == 0) return false;
        scratch.reset_per_fragment();

        cr.genomic_footprint = genomic_footprint;

        // --- Interchromosomal chimera detection ---
        cr.chimera_type = CHIMERA_NONE;
        cr.chimera_gap = -1;
        scratch.tmp_a.clear();
        scratch.tmp_a.reserve(static_cast<size_t>(n_exons));
        for (const auto& e : exons) scratch.tmp_a.push_back(e.ref_id);
        sort_unique(scratch.tmp_a);
        bool is_interchromosomal = (scratch.tmp_a.size() > 1);
        if (is_interchromosomal) cr.chimera_type = CHIMERA_TRANS;

        // --- Query each exon block ---
        int n_introns = static_cast<int>(introns.size());
        auto& exon_t_sets = scratch.exon_t_sets;
        auto& transcript_t_sets = scratch.transcript_t_sets;
        exon_t_sets.resize(static_cast<size_t>(n_exons));
        transcript_t_sets.resize(static_cast<size_t>(n_exons));
        cr.align_strand = STRAND_NONE;
        cr.read_length = 0;

        for (int bi = 0; bi < n_exons; bi++) {
            const auto& eb = exons[bi];
            cr.align_strand |= eb.strand;
            int32_t bstart = eb.start;
            int32_t bend = eb.end;
            cr.read_length += bend - bstart;

            if (eb.ref_id < 0 ||
                eb.ref_id >= static_cast<int32_t>(id_to_ref_.size())) {
                continue;
            }

            auto& block_exon_t = exon_t_sets[static_cast<size_t>(bi)];
            auto& block_transcript_t = transcript_t_sets[static_cast<size_t>(bi)];
            block_exon_t.clear();
            block_transcript_t.clear();

            const char* ref_str = id_to_ref_[eb.ref_id].c_str();
            int64_t n = cr_overlap(cr_, ref_str, bstart, bend,
                                   &scratch.buf, &scratch.buf_cap);

            for (int64_t hi = 0; hi < n; hi++) {
                int64_t idx = scratch.buf[hi];
                int32_t h_start = cr_start(cr_, idx);
                int32_t h_end   = cr_end(cr_, idx);
                int32_t label   = cr_label(cr_, idx);
                int8_t  itype   = iv_type_[label];

                int32_t off = t_set_offsets_[label];
                int32_t cnt = t_set_offsets_[label + 1] - off;

                if (itype == ITYPE_EXON) {
                    int32_t clo = std::max(bstart, h_start);
                    int32_t chi = std::min(bend, h_end);
                    int32_t bp = (chi > clo) ? (chi - clo) : 0;
                    bool any_pos = false, any_neg = false;
                    for (int32_t k = 0; k < cnt; k++) {
                        int32_t ti = t_set_data_[off + k];
                        block_exon_t.push_back(ti);
                        if (bp > 0) {
                            scratch.mark_dirty(ti);
                            scratch.t_exon_bp[ti] += bp;
                            int32_t s = (ti >= 0 && ti < static_cast<int32_t>(t_strand_arr_.size()))
                                        ? t_strand_arr_[ti] : STRAND_NONE;
                            if (s == STRAND_POS) any_pos = true;
                            else if (s == STRAND_NEG) any_neg = true;
                        }
                    }
                    if (bp > 0) {
                        if (any_pos) cr.exon_bp_pos += bp;
                        if (any_neg) cr.exon_bp_neg += bp;
                    }
                } else if (itype == ITYPE_TRANSCRIPT) {
                    int32_t clo = std::max(bstart, h_start);
                    int32_t chi = std::min(bend, h_end);
                    int32_t bp = (chi > clo) ? (chi - clo) : 0;
                    bool any_pos = false, any_neg = false;
                    for (int32_t k = 0; k < cnt; k++) {
                        int32_t ti = t_set_data_[off + k];
                        block_transcript_t.push_back(ti);
                        if (bp > 0) {
                            scratch.mark_dirty(ti);
                            scratch.t_transcript_bp[ti] += bp;
                            int32_t s = (ti >= 0 && ti < static_cast<int32_t>(t_strand_arr_.size()))
                                        ? t_strand_arr_[ti] : STRAND_NONE;
                            if (s == STRAND_POS) any_pos = true;
                            else if (s == STRAND_NEG) any_neg = true;
                        }
                    }
                    if (bp > 0) {
                        if (any_pos) cr.tx_bp_pos += bp;
                        if (any_neg) cr.tx_bp_neg += bp;
                    }
                }
            }

            sort_unique(block_exon_t);
            sort_unique(block_transcript_t);
        }

        // --- Intrachromosomal chimera detection ---
        if (!is_interchromosomal) {
            auto cr_res = detect_chimera(exons, exon_t_sets);
            if (cr_res.type != CHIMERA_NONE) {
                cr.chimera_type = cr_res.type;
                cr.chimera_gap = cr_res.gap;
            }
        }

        // --- Derive nRNA candidates from real-tx hits ---------------
        // Synthetic nRNAs are intentionally absent from the cgranges
        // index (see _gen_transcript_intervals).  Each real-tx hit
        // carries a one-to-one link to its parent nRNA entity via
        // nrna_parent_; collect the unique parents the fragment
        // overlaps, accumulate fragment-bp inside each parent's span,
        // and inject them into every per-block exon_t_set so that
        // merge_sets() naturally picks them up.  Run AFTER chimera
        // detection (nRNAs would otherwise mask cis-chimeras since
        // they cover both blocks by construction) and AFTER the
        // strand-aware exon_bp_pos/_neg accumulation (nRNAs must not
        // pollute the calibration overlap counts).
        if (!nrna_parent_.empty()) {
            auto add_nrna = [&](int32_t ti) {
                if (ti < 0 ||
                    ti >= static_cast<int32_t>(nrna_parent_.size())) return;
                int32_t n = nrna_parent_[ti];
                if (n >= 0 && n != ti) scratch.nrna_t.push_back(n);
            };
            for (const auto& s : exon_t_sets)
                for (int32_t ti : s) add_nrna(ti);
            for (const auto& s : transcript_t_sets)
                for (int32_t ti : s) add_nrna(ti);

            sort_unique(scratch.nrna_t);
            for (int32_t n : scratch.nrna_t) {
                // Synthetic nRNAs have exactly one exon by construction;
                // annotated nascent-equivs are handled via cgranges and
                // never appear here (their nrna_parent_[ti] == ti was
                // filtered out by the n != ti guard above).
                int32_t ns = exon_starts_[exon_offsets_[n]];
                int32_t ne = exon_ends_  [exon_offsets_[n]];
                int32_t bp_total = 0;
                for (const auto& eb : exons) {
                    int32_t lo = std::max(eb.start, ns);
                    int32_t hi = std::min(eb.end,   ne);
                    if (hi > lo) bp_total += (hi - lo);
                }
                if (bp_total <= 0) continue;
                scratch.mark_dirty(n);
                scratch.t_exon_bp[n] += bp_total;
                scratch.t_transcript_bp[n] += bp_total;
                // Insert into every block's exon_t_set so merge_sets'
                // intersection logic surfaces the nRNA across blocks.
                // Sets are sorted; use lower_bound to keep them sorted
                // and unique.
                for (auto& bset : exon_t_sets) {
                    auto it = std::lower_bound(bset.begin(), bset.end(), n);
                    if (it == bset.end() || *it != n)
                        bset.insert(it, n);
                }
            }
        }

        // --- SJ matching ---
        auto& sj_t_sets_vec = scratch.sj_t_sets;
        bool has_annotated_sj = false;
        bool has_unannotated_sj = false;
        cr.sj_strand = STRAND_NONE;
        // Reset alongside sj_strand: `cr` is an out-param taken by reference, so a
        // caller that reuses one across fragments would otherwise carry the previous
        // fragment's junction key forward. Every call site builds a fresh `cr` today.
        cr.sj_key_ref = cr.sj_key_start = cr.sj_key_end = -1;

        for (int ii = 0; ii < n_introns; ii++) {
            sj_t_sets_vec.emplace_back();
            auto& sj_t = sj_t_sets_vec.back();
            sj_lookup_into(introns[ii].ref_id, introns[ii].start,
                           introns[ii].end, introns[ii].strand, sj_t);
            if (!sj_t.empty()) {
                has_annotated_sj = true;
                cr.sj_strand |= introns[ii].strand;
                if (cr.sj_key_start < 0) {
                    // First ANNOTATED junction (introns arrive sorted by
                    // ref/start/end) — the per-junction strand table's key.
                    cr.sj_key_ref   = introns[ii].ref_id;
                    cr.sj_key_start = introns[ii].start;
                    cr.sj_key_end   = introns[ii].end;
                }
            } else {
                sj_t_sets_vec.pop_back();
                has_unannotated_sj = true;
            }
        }

        // --- Resolution: merge and classify ---
        bool any_exon = false, any_transcript = false;
        for (const auto& s : exon_t_sets)
            if (!s.empty()) { any_exon = true; break; }
        for (const auto& s : transcript_t_sets)
            if (!s.empty()) { any_transcript = true; break; }

        if (has_annotated_sj && any_exon) {
            auto exon_merge = merge_sets(exon_t_sets);
            auto sj_merge   = merge_sets(sj_t_sets_vec);

            if (!sj_merge.is_empty()) {
                if (!exon_merge.is_empty()) {
                    auto isect = vec_intersect(exon_merge.t_inds,
                                               sj_merge.t_inds);
                    cr.t_inds = isect.empty() ? std::move(sj_merge.t_inds)
                                              : std::move(isect);
                } else {
                    cr.t_inds = std::move(sj_merge.t_inds);
                }
                cr.merge_criteria = sj_merge.criteria;
            } else {
                std::vector<std::vector<int32_t>> all_sets;
                all_sets.insert(all_sets.end(),
                                exon_t_sets.begin(), exon_t_sets.end());
                all_sets.insert(all_sets.end(),
                                sj_t_sets_vec.begin(), sj_t_sets_vec.end());
                auto mr = merge_sets(all_sets);
                cr.t_inds = std::move(mr.t_inds);
                cr.merge_criteria = mr.criteria;
            }
            cr.splice_type = SPLICE_SPLICED_ANNOT;

        } else if (any_exon || any_transcript) {
            std::vector<std::vector<int32_t>> parts;
            int32_t best_criteria = MC_UNION;

            if (any_exon) {
                auto em = merge_sets(exon_t_sets);
                if (!em.is_empty()) {
                    parts.push_back(std::move(em.t_inds));
                    best_criteria = em.criteria;
                }
            }
            if (any_transcript) {
                auto tm = merge_sets(transcript_t_sets);
                if (!tm.is_empty())
                    parts.push_back(std::move(tm.t_inds));
            }

            if (!parts.empty()) {
                cr.t_inds = std::move(parts[0]);
                for (size_t pi = 1; pi < parts.size(); pi++) {
                    scratch.tmp_union.clear();
                    scratch.tmp_union.reserve(cr.t_inds.size() + parts[pi].size());
                    std::set_union(cr.t_inds.begin(), cr.t_inds.end(),
                                   parts[pi].begin(), parts[pi].end(),
                                   std::back_inserter(scratch.tmp_union));
                    cr.t_inds.assign(scratch.tmp_union.begin(),
                                     scratch.tmp_union.end());
                }
                cr.merge_criteria =
                    (parts.size() == 1) ? best_criteria : MC_UNION;
            }

            cr.splice_type = has_unannotated_sj ? SPLICE_SPLICED_UNANNOT
                                                : SPLICE_UNSPLICED;

        } else {
            // No exon and no transcript overlap -> truly intergenic.
            // SRD v2: do NOT early-exit; let the fragment flow through
            // the buffer with empty t_inds so calibration can categorize
            // it as INTERGENIC.
            cr.splice_type = has_unannotated_sj ? SPLICE_SPLICED_UNANNOT
                                                : SPLICE_UNSPLICED;
        }

        // SRD v2: SPLICE_ARTIFACT promotion.  When the alignment had
        // CIGAR splice junctions that were rejected by the artifact
        // blacklist (set by caller in cr.n_sj_blacklisted before this
        // call), promote any non-spliced classification so the
        // fragment is held out of calibration.
        if (cr.n_sj_blacklisted > 0 &&
            cr.splice_type != SPLICE_SPLICED_ANNOT &&
            cr.splice_type != SPLICE_SPLICED_UNANNOT) {
            cr.splice_type = SPLICE_ARTIFACT;
        }

        // --- Overlap profiles (parallel to t_inds) ---
        // Skip when t_inds is empty (intergenic / unresolved).
        auto& all_overlap_t = scratch.all_overlap_t;
        all_overlap_t.clear();
        for (const auto& s : exon_t_sets)
            all_overlap_t.insert(all_overlap_t.end(), s.begin(), s.end());
        for (const auto& s : transcript_t_sets)
            all_overlap_t.insert(all_overlap_t.end(), s.begin(), s.end());
        sort_unique(all_overlap_t);

        bool any_overlap = !all_overlap_t.empty() && cr.read_length > 0;
        cr.t_exon_bp.reserve(cr.t_inds.size());
        cr.t_intron_bp.reserve(cr.t_inds.size());
        for (int32_t t : cr.t_inds) {
            if (any_overlap && std::binary_search(
                    all_overlap_t.begin(), all_overlap_t.end(), t)) {
                cr.t_exon_bp.push_back(scratch.t_exon_bp[t]);
                cr.t_intron_bp.push_back(
                    std::max(scratch.t_transcript_bp[t] - scratch.t_exon_bp[t], 0));
            } else if (!any_overlap) {
                cr.t_exon_bp.push_back(0);
                cr.t_intron_bp.push_back(cr.read_length);
            } else {
                cr.t_exon_bp.push_back(cr.read_length);
                cr.t_intron_bp.push_back(0);
            }
        }

        // --- ambig_strand ---
        int32_t first_strand = -999;
        bool mixed = false;
        for (int32_t t : cr.t_inds) {
            if (t >= 0 && t < static_cast<int32_t>(t_strand_arr_.size())) {
                int32_t s = t_strand_arr_[t];
                if (first_strand == -999) {
                    first_strand = s;
                } else if (s != first_strand) {
                    mixed = true;
                    break;
                }
            }
        }
        cr.ambig_strand = mixed ? 1 : 0;

        // --- Fragment lengths ---
        if (cr.chimera_type == CHIMERA_NONE) {
            compute_frag_lengths_aligned(exons, introns, cr.t_inds,
                                         cr.frag_lengths, scratch);
        }

        // --- gap-hypothesis enumeration ---
        //
        // ⭐ EVERY fragment, whatever its splice type. It used to run only on fragments already
        // classified SPLICE_UNSPLICED, so one carrying an observed CIGAR-N splice never had its mate gap
        // examined and kept that intron inside L. ⚠ UNSPLICED never meant "one aligned block": an
        // unspliced paired-end fragment already has two blocks and a mate gap, and that case always
        // worked. The missed population is SPLICED fragments that ALSO have a gap intron -- long by
        // construction, and so exactly the tail the measurement found.
        //
        // ⛔ THE ACCUMULATOR ARBITRATES. This only enumerates; whether the fragment deposits or is held
        // for the second pass is decided where the outcome is reported.
        //
        // ⚠ `chimera_type == CHIMERA_NONE` is KEPT: a chimeric fragment's blocks are not one molecule, so
        // its gaps are not introns. Such a fragment gets the unspliced hypothesis alone.
        if (cr.chimera_type == CHIMERA_NONE) {
            enumerate_gap_hypotheses(exons, introns, cr.t_inds,
                                     cr.splice_type == SPLICE_SPLICED_ANNOT, scratch, cr);
        } else {
            cr.gap_introns.clear();
            cr.gap_sj_strand.assign(1, STRAND_NONE);
            cr.gap_intron_offsets.assign(2, 0);
            cr.gap_supporting.clear();
            cr.gap_supporting_offsets.assign(2, 0);
        }

        // ⛔ THE SPLICE_IMPLICIT PROMOTION STAYS UNSPLICED-ONLY, and it is now purely DESCRIPTIVE.
        // `splice_type` is the scanner's census of what it SAW; it feeds scoring, the buffer, the strand
        // training and `rigel report`. Re-labelling an observed SPLICED_ANNOT fragment would silently move
        // mass between reported categories. ⚠ The fragment may still be DEFERRED -- the two axes are
        // independent, which is why the umbrella census is its own axis and not a splice type.
        if (cr.splice_type == SPLICE_UNSPLICED && cr.n_gap_hypotheses() > 1) {
            cr.splice_type = SPLICE_IMPLICIT;
        }

        // --- Genomic start ---
        cr.genomic_start = exons[0].start;
        for (const auto& e : exons)
            cr.genomic_start = std::min(cr.genomic_start, e.start);

        scratch.clean();
        return true;
    }

    /// Backward-compatible wrapper using internal scratch.
    bool _resolve_core(
        const std::vector<ExonBlock>& exons,
        const std::vector<IntronBlock>& introns,
        int32_t genomic_footprint,
        RawResolveResult& cr)
    {
        return _resolve_core(exons, introns, genomic_footprint, cr, scratch_);
    }

    // ================================================================
    // resolve_fragment — accepts Python Fragment, returns ResolvedFragment
    // ================================================================

    nb::object resolve_fragment(nb::object frag) {
        nb::object exons_obj = frag.attr("exons");
        size_t n_exons = nb::len(exons_obj);
        if (n_exons == 0) return nb::none();

        std::vector<ExonBlock> exons(n_exons);
        for (size_t i = 0; i < n_exons; i++) {
            nb::object e = exons_obj[nb::int_(i)];
            std::string ref = nb::cast<std::string>(e[nb::int_(0)]);
            auto it = ref_to_id_.find(ref);
            exons[i] = {
                (it != ref_to_id_.end()) ? it->second : -1,
                nb::cast<int32_t>(e[nb::int_(1)]),
                nb::cast<int32_t>(e[nb::int_(2)]),
                nb::cast<int32_t>(e[nb::int_(3)])
            };
        }

        nb::object introns_obj = frag.attr("introns");
        size_t n_introns = nb::len(introns_obj);
        std::vector<IntronBlock> introns(n_introns);
        for (size_t i = 0; i < n_introns; i++) {
            nb::object intr = introns_obj[nb::int_(i)];
            std::string ref = nb::cast<std::string>(intr[nb::int_(0)]);
            auto it = ref_to_id_.find(ref);
            introns[i] = {
                (it != ref_to_id_.end()) ? it->second : -1,
                nb::cast<int32_t>(intr[nb::int_(1)]),
                nb::cast<int32_t>(intr[nb::int_(2)]),
                nb::cast<int32_t>(intr[nb::int_(3)])
            };
        }

        int32_t genomic_footprint = nb::cast<int32_t>(
            frag.attr("genomic_footprint"));

        RawResolveResult cr;
        if (!_resolve_core(exons, introns, genomic_footprint, cr))
            return nb::none();

        return nb::cast(ResolvedFragment::from_core(cr));
    }
};

}  // namespace rigel
