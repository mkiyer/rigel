/**
 * resolve_context.h — FragmentResolver, ResolvedFragment, FragmentAccumulator.
 *
 * Shared by resolve.cpp (the Python binding) and bam_scanner.cpp, which
 * calls _resolve_core() directly in C++ and accumulates results into
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
#include <type_traits>
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

static inline uint16_t checked_u16_buffer_value(int32_t value, const char* column)
{
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
    int32_t num_hits = 1;
    int32_t nm = 0;

    // Leftmost ANNOTATED sj (ref, start, end); -1 when the fragment
    // crosses none.  With `sj_strand` this is the per-sj strand
    // table's key — see RawResolveResult in constants.h for why it is the
    // first ANNOTATED sj and why only one sj is credited.
    int32_t sj_key_ref = -1;
    int32_t sj_key_start = -1;
    int32_t sj_key_end = -1;

    // Per-transcript parallel arrays (same order as t_inds)
    std::vector<int32_t> frag_lengths;      // -1 if missing
    std::vector<int32_t> exon_bp;
    std::vector<int32_t> intron_bp;

    // The annotated introns found inside this fragment's UNSEQUENCED gaps, one run per gap hypothesis
    // (`gap_intron_offsets`; an empty run is the unspliced hypothesis). Carried out of `RawResolveResult`
    // unchanged.
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

    // Return t_inds as a Python frozenset
    nb::object get_t_inds() const {
        nb::set s;
        for (int32_t t : t_inds) s.add(nb::cast(t));
        return nb::frozenset(s);
    }

    // Return frag_lengths as a Python dict {t_ind: length}, positive lengths only
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

    // Return overlap_bp as a Python dict {t_ind: (exon_bp, intron_bp)}
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
        r.exon_bp = std::move(cr.t_exon_bp);
        r.intron_bp = std::move(cr.t_intron_bp);
        r.sj_key_ref = cr.sj_key_ref;
        r.sj_key_start = cr.sj_key_start;
        r.sj_key_end = cr.sj_key_end;
        // ⛔ COPIED, NOT MOVED, and the difference is a bug. `bam_scanner.cpp` builds the
        // ResolvedFragment first and calls `collect_placement(frag, cr)` afterwards, so moving
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
// FragmentAccumulator — C++ columnar buffer
// ================================================================

class FragmentAccumulator {
public:
    std::vector<uint8_t>  splice_type_;
    std::vector<uint8_t>  align_strand_;
    std::vector<uint8_t>  ambig_strand_;
    std::vector<uint16_t> num_hits_;
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
        ambig_strand_.reserve(n_fragments);
        num_hits_.reserve(n_fragments);
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
        ambig_strand_.push_back(static_cast<uint8_t>(r.ambig_strand));
        num_hits_.push_back(static_cast<uint16_t>(r.num_hits));
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

    /// Finalize: moves internal vectors to heap-allocated storage and
    /// returns capsule-backed numpy arrays.  The accumulator is consumed
    /// (left empty) after this call.
    nb::dict finalize() {
        int32_t n = size_;

        nb::dict result;
        result["splice_type"]       = vec_to_ndarray(std::move(splice_type_));
        result["align_strand"]       = vec_to_ndarray(std::move(align_strand_));
        result["num_hits"]          = vec_to_ndarray(std::move(num_hits_));
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
// Kept apart from FragmentResolver so that multiple threads can call
// _resolve_core() concurrently, each with its own scratch. The single-threaded path uses FragmentResolver's
// internal scratch_ member.

struct GapBlock {
    int32_t ref_id;
    int32_t start;
    int32_t end;
};

/// One mate of a candidate pair, as `FragmentResolver::validate_pair` reads it: its aligned blocks, split at each
/// CIGAR N and ascending, so the holes between consecutive blocks are its observed skips; and the soft clips at the
/// alignment's two ends (the read's 5' clip is the leading one on a forward mate, the trailing one on a reverse mate).
struct PairMate {
    const std::vector<std::pair<int32_t, int32_t>>* blocks;
    int32_t clip_lead;
    int32_t clip_trail;
};

/// The cgranges query buffer and its capacity, owned. `cr_overlap` grows it with `realloc` through the
/// addresses of `ptr` and `cap`, so it is released with `free`. A move hands both over and leaves the source
/// empty, as a fresh buffer is: a source keeping its `cap` with no `ptr` would be written through null by
/// its next query, because cgranges reallocates only once the count reaches `cap`.
struct OverlapQueryBuffer {
    int64_t* ptr = nullptr;
    int64_t  cap = 0;

    OverlapQueryBuffer() = default;
    OverlapQueryBuffer(const OverlapQueryBuffer&) = delete;
    OverlapQueryBuffer& operator=(const OverlapQueryBuffer&) = delete;
    OverlapQueryBuffer(OverlapQueryBuffer&& o) noexcept
        : ptr(std::exchange(o.ptr, nullptr)), cap(std::exchange(o.cap, 0)) {}
    OverlapQueryBuffer& operator=(OverlapQueryBuffer&& o) noexcept {
        if (this != &o) {
            free(ptr);
            ptr = std::exchange(o.ptr, nullptr);
            cap = std::exchange(o.cap, 0);
        }
        return *this;
    }
    ~OverlapQueryBuffer() { free(ptr); }
};

struct ResolverScratch {
    // Per-transcript bp counters (sparse-cleaned)
    std::vector<int32_t> t_exon_bp;
    std::vector<int32_t> t_transcript_bp;
    std::vector<uint8_t> t_dirty;
    std::vector<int32_t> dirty_indices;
    // dirty_indices' initial capacity, in transcripts one fragment touches. A capacity hint only; it grows
    // past it.
    static constexpr size_t DIRTY_INDICES_RESERVE = 512;

    // cgranges query buffer (reusable per-call)
    OverlapQueryBuffer buf;

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
    std::vector<int32_t> tmp_union;
    std::vector<int32_t> all_overlap_t;
    std::vector<int32_t> nrna_t;

    // Reusable scratch for the pair validator (`validate_pair`)
    std::vector<std::pair<int32_t, int32_t>> pair_skips;
    std::vector<std::pair<int32_t, int32_t>> pair_aligned;
    std::vector<std::pair<int32_t, int32_t>> pair_covered;
    std::vector<GapBlock> pair_stretches;
    std::vector<IntronBlock> pair_cut;
    std::vector<IntronBlock> hit_junctions;  // every junction a hit observes (the scan's `observed_junctions`)
    std::vector<int32_t> pair_candidates;
    std::vector<int32_t> pair_tmp;

    ResolverScratch() = default;

    explicit ResolverScratch(int32_t n_transcripts)
        : t_exon_bp(n_transcripts, 0),
          t_transcript_bp(n_transcripts, 0),
          t_dirty(n_transcripts, 0)
    {
        dirty_indices.reserve(DIRTY_INDICES_RESERVE);
    }

    // Non-copyable. ⚠ The moves are the members' own, defaulted, so a member added later cannot be left out
    // of them, as a hand-written member list once left out `hypothesis_path`.
    ResolverScratch(const ResolverScratch&) = delete;
    ResolverScratch& operator=(const ResolverScratch&) = delete;
    ResolverScratch(ResolverScratch&&) = default;
    ResolverScratch& operator=(ResolverScratch&&) = default;

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
        tmp_union.clear();
        all_overlap_t.clear();
        nrna_t.clear();
        implicit_blocks.clear();
        implicit_gaps.clear();
    }
};

// The defaulted moves take their exception specification from the members, so this holds them to the
// `noexcept` the hand-written moves promised: a member whose move can throw fails here, at compile time.
static_assert(std::is_nothrow_move_constructible_v<ResolverScratch> &&
                  std::is_nothrow_move_assignable_v<ResolverScratch>,
              "ResolverScratch: every member must move without throwing");

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
    // sj is rejected at BAM-scan time when either its left or
    // right anchor is <= the stored maximum.
    using SJBlacklistMap = std::unordered_map<
        SJBlacklistKey, SJBlacklistEntry, SJBlacklistKeyHash>;
    SJBlacklistMap sj_blacklist_;



    // --- Transcript metadata ---
    int32_t n_transcripts_ = 0;

    // --- Transcript strand metadata (direct lookup, no gene indirection) ---
    std::vector<int32_t> t_strand_arr_;

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

    // --- Intron-overhang tolerance K (bp), the gap-hypothesis enumeration's slack ---
    int32_t intron_overhang_tolerance_ = 0;
    //: the library's fragment-length limit; 0 ⇒ no molecule is compatible, so the compatibility
    //: rescue in `detect_chimera` cannot fire and every disjoint candidate stays a chimera.
    int32_t max_fragment_length_ = 0;

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
    /// sj; rejection uses the OR rule (either anchor <= max).
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

    /// Return blacklist entry for a sj, or nullptr if not present.
    inline const SJBlacklistEntry* sj_blacklist_lookup(
        int32_t ref_id, int32_t start, int32_t end) const
    {
        auto it = sj_blacklist_.find(SJBlacklistKey{ref_id, start, end});
        return (it == sj_blacklist_.end()) ? nullptr : &it->second;
    }

    void set_metadata(int32_t n_transcripts) {
        n_transcripts_ = n_transcripts;

        // Initialize internal scratch buffers
        scratch_ = ResolverScratch(n_transcripts_);
    }

    /// Set per-transcript strand array (direct lookup, no gene indirection).
    void set_transcript_strands(const std::vector<int32_t>& t_strand) {
        t_strand_arr_ = t_strand;
    }

    /// Set per-transcript nRNA parent index (int32; -1 = no parent).
    /// Used by _resolve_core to derive synthetic-nRNA candidates from
    /// real-transcript hits without bloating the cgranges index.
    void set_nrna_parent_index(const std::vector<int32_t>& nrna_parent) {
        nrna_parent_ = nrna_parent;
    }

    /// Build per-transcript exon CSR index for FL computation.
    void build_exon_index(
        const std::vector<int32_t>& offsets,
        const std::vector<int32_t>& starts,
        const std::vector<int32_t>& ends,
        const std::vector<int32_t>& cumsum)
    {
        exon_offsets_ = offsets;
        exon_starts_ = starts;
        exon_ends_ = ends;
        exon_cumsum_ = cumsum;
    }

    bool has_exon_index() const {
        return !exon_offsets_.empty();
    }

    // ----------------------------------------------------------------
    // Chimera detection
    // ----------------------------------------------------------------

    /// The library's fragment-length limit (`BamScanConfig.max_frag_length`), used by
    /// `detect_chimera` to tell an ordinary genomic molecule from a rearrangement.
    void set_max_fragment_length(int32_t L) {
        if (L < 0) throw std::invalid_argument(
            "set_max_fragment_length: must be >= 0, got " + std::to_string(L));
        max_fragment_length_ = L;
    }

    // ----------------------------------------------------------------
    // SPLICED_IMPLICIT detection
    // ----------------------------------------------------------------

    /// Set the intron-overhang tolerance K (bp): how many bases an alignment may run into an annotated
    /// intron, and an annotated intron overhang a mate gap, and still be read as spliced there. It is the
    /// slack in both tests of the gap-hypothesis enumeration, which every non-chimeric fragment with an
    /// unsequenced gap goes through: an annotated intron in a gap may overlap up to K aligned bases on
    /// either side and is still cut whole (`collect_transcript_introns_in_gap`), and a transcript is
    /// contradicted only by an intron overlapping more than K bases of an aligned block
    /// (`transcript_contradicts_blocks`). Negative K is clamped to 0; the resolver's own default, before
    /// this is called, is 0.
    void set_intron_overhang_tolerance(int32_t K) {
        if (K < 0) K = 0;
        intron_overhang_tolerance_ = K;
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
    /// INTRON-OVERHANG TOLERANCE K before it counts -- the same slack the gap match uses, because it is
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

    // ----------------------------------------------------------------
    // The pair validator: two mates are one molecule, or why not
    // ----------------------------------------------------------------

    /// Every annotated transcript whose span overlaps [start, end) on ``ref_id``, sorted and unique. The synthetic
    /// nascent spans are not in the overlap index, and have no introns to offer.
    void transcripts_overlapping(int32_t ref_id, int32_t start, int32_t end,
                                 ResolverScratch& scratch, std::vector<int32_t>& out) const
    {
        out.clear();
        if (!cr_ || ref_id < 0 || ref_id >= static_cast<int32_t>(id_to_ref_.size())) return;
        const int64_t n = cr_overlap(cr_, id_to_ref_[ref_id].c_str(), start, end,
                                     &scratch.buf.ptr, &scratch.buf.cap);
        for (int64_t k = 0; k < n; ++k) {
            const int32_t label = cr_label(cr_, scratch.buf.ptr[k]);
            if (iv_type_[label] != ITYPE_TRANSCRIPT) continue;
            for (int32_t i = t_set_offsets_[label]; i < t_set_offsets_[label + 1]; ++i)
                out.push_back(t_set_data_[i]);
        }
        sort_unique(out);
    }

    /// Whether the forward mate ``f`` and the reverse mate ``r``, on reference ``ref_id``, can be one molecule:
    /// `HIT_USABLE`, or `HIT_OUTWARD` or `HIT_LENGTH` (the caller has checked the reference and the strands).
    ///
    /// The forward mate's 5' soft clip is walked back along the genome and the reverse mate's forward, giving the
    /// molecule [a, b); the mates must face inward, with both mates' aligned bases inside it. Its length counts the
    /// bases of [a, b) with every observed skip of either mate removed, and must lie in [1, max_fragment_length].
    /// Failing that, one annotated transcript that holds every junction the hit observes (``junctions``, matched
    /// as `sj_lookup_into` matches them), or, with none observed, any transcript over the molecule, may witness a
    /// shorter path: its other introns in the stretches neither mate aligned are cut whole by the K rule of
    /// `collect_transcript_introns_in_gap`, and it offers no path when one of them overlaps an aligned block by more
    /// than K. One transcript supplies both the observed junctions and the implied introns, so a junction one
    /// transcript holds and an intron only another holds never make a path together (the targeted review's
    /// hybrid). Measured against per-base truth on 47.9 M STAR templates, this rejects no correct pair.
    int32_t validate_pair(int32_t ref_id, const PairMate& f, const PairMate& r,
                          const std::vector<IntronBlock>& junctions, ResolverScratch& scratch) const
    {
        const int32_t K = intron_overhang_tolerance_;
        const int64_t max_len = max_fragment_length_;
        const int32_t a = f.blocks->front().first - f.clip_lead;
        const int32_t b = r.blocks->back().second + r.clip_trail;
        if (!(a < b && f.blocks->back().second <= b && r.blocks->front().first >= a)) return HIT_OUTWARD;

        auto& skips = scratch.pair_skips;  // each mate's holes between consecutive blocks, the two merged
        skips.clear();
        for (const PairMate* m : {&f, &r})
            for (size_t i = 1; i < m->blocks->size(); ++i)
                skips.push_back({(*m->blocks)[i - 1].second, (*m->blocks)[i].first});
        std::sort(skips.begin(), skips.end());
        skips.erase(std::unique(skips.begin(), skips.end()), skips.end());
        int64_t removed = 0;  // the union of the skips inside [a, b)
        int32_t run_s = 0, run_e = 0;
        bool open = false;
        for (auto [s, e] : skips) {
            s = std::max(s, a);
            e = std::min(e, b);
            if (e <= s) continue;
            if (open && s <= run_e) {
                run_e = std::max(run_e, e);
            } else {
                if (open) removed += run_e - run_s;
                run_s = s;
                run_e = e;
                open = true;
            }
        }
        if (open) removed += run_e - run_s;
        const int64_t length = static_cast<int64_t>(b) - a - removed;
        if (length >= 1 && length <= max_len) return HIT_USABLE;

        auto& cands = scratch.pair_candidates;  // the witnesses
        cands.clear();
        if (junctions.empty()) {
            transcripts_overlapping(ref_id, a, b, scratch, cands);
        } else {
            for (size_t j = 0; j < junctions.size(); ++j) {
                sj_lookup_into(junctions[j].ref_id, junctions[j].start, junctions[j].end,
                               junctions[j].strand, scratch.pair_tmp);
                if (j == 0) {
                    cands.swap(scratch.pair_tmp);
                } else {
                    auto end = std::set_intersection(cands.begin(), cands.end(), scratch.pair_tmp.begin(),
                                                     scratch.pair_tmp.end(), cands.begin());
                    cands.erase(end, cands.end());
                }
                if (cands.empty()) return HIT_LENGTH;
            }
        }

        auto& aligned = scratch.pair_aligned;  // both mates' aligned blocks
        aligned.assign(f.blocks->begin(), f.blocks->end());
        aligned.insert(aligned.end(), r.blocks->begin(), r.blocks->end());
        std::sort(aligned.begin(), aligned.end());
        auto& covered = scratch.pair_covered;  // aligned or skipped
        covered.assign(aligned.begin(), aligned.end());
        covered.insert(covered.end(), skips.begin(), skips.end());
        std::sort(covered.begin(), covered.end());
        auto& stretches = scratch.pair_stretches;  // what neither mate aligned or skipped
        stretches.clear();
        int32_t cur = a;
        for (const auto& [s, e] : covered) {
            if (s > cur) stretches.push_back({ref_id, cur, std::min(s, b)});
            cur = std::max(cur, e);
        }
        if (cur < b) stretches.push_back({ref_id, cur, b});

        for (int32_t t : cands) {
            if (t < 0 || t + 1 >= static_cast<int32_t>(exon_offsets_.size())) continue;
            const int32_t begin = exon_offsets_[t];
            const int32_t n_introns = exon_offsets_[t + 1] - begin - 1;
            const int32_t* intron_starts = exon_ends_.data() + begin;
            const int32_t* intron_ends = exon_starts_.data() + begin + 1;
            bool path = true;
            for (int32_t i = 0; i < n_introns && path; ++i) {
                const int32_t i0 = intron_starts[i], i1 = intron_ends[i];
                if (i1 <= a || i0 >= b) continue;
                if (std::binary_search(skips.begin(), skips.end(), std::make_pair(i0, i1))) continue;  // observed
                for (const auto& [x, y] : aligned)
                    if (std::min(i1, y) - std::max(i0, x) > K) { path = false; break; }
            }
            if (!path) continue;
            auto& cut = scratch.pair_cut;
            cut.clear();
            for (const GapBlock& g : stretches) collect_transcript_introns_in_gap(t, g, K, STRAND_NONE, cut);
            std::sort(cut.begin(), cut.end(), [](const IntronBlock& x, const IntronBlock& y) {
                return x.start != y.start ? x.start < y.start : x.end < y.end;
            });
            int64_t extra = 0;
            for (size_t i = 0; i < cut.size(); ++i)
                if (i == 0 || cut[i].start != cut[i - 1].start || cut[i].end != cut[i - 1].end)
                    extra += cut[i].end - cut[i].start;
            if (length - extra >= 1 && length - extra <= max_len) return HIT_USABLE;
        }
        return HIT_LENGTH;
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
    /// Searching it again is not merely redundant: the ±K intron-overhang tolerance lets a DIFFERENT
    /// nearby intron answer for it, the two normalise into one wider interval, and L comes out too SHORT.
    /// Gaps matching an observed intron by EXACT (start, end) equality are dropped -- never by overlap,
    /// which would also drop a genuine gap intron abutting an observed one.
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
        cr.gap_intron_offsets.assign(1, 0);

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
        const int32_t K = intron_overhang_tolerance_;
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
            add_hypothesis(cr, path, strand);
        }

        // ⭐ The UNSPLICED hypothesis. Present when some compatible transcript implies nothing in the
        // gaps (a retained-intron isoform), and ALWAYS when no annotated sj was sequenced -- then
        // the molecule may be gDNA or RNA that has not spliced there, and the gap is real template.
        if (any_candidate_implies_nothing || !certified_rna || cr.n_gap_hypotheses() == 0)
            emit_unspliced_hypothesis(cr);
    }

private:
    /// Append every intron of ``t`` lying inside ``gap`` (±K), in genomic order.
    ///
    /// ⛔ EVERY one, not the first. Stopping at the first would cut only one of two annotated introns a mate
    /// gap spans, and two transcripts differing only in their SECOND intron would read as agreeing.
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

    /// Record ``path`` as a hypothesis, or fold its strand into the identical one already there.
    static void add_hypothesis(RawResolveResult& cr,
                               const std::vector<IntronBlock>& path,
                               int32_t strand)
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
            // sj on it. ⚠ Idempotent: once AMBIGUOUS, any further disagreement keeps it there.
            //
            // ⚠ An sj annotated on both strands is biologically impossible, and the index warns about
            // one. Handled because the alternative is a silent answer that flips when two GTF
            // lines are swapped.
            if (cr.gap_sj_strand[h] != strand) cr.gap_sj_strand[h] = STRAND_AMBIGUOUS;
            return;
        }
        cr.gap_introns.insert(cr.gap_introns.end(), path.begin(), path.end());
        cr.gap_intron_offsets.push_back(static_cast<int32_t>(cr.gap_introns.size()));
        cr.gap_sj_strand.push_back(strand);
    }

    /// The empty path: cut nothing. ⚠ Idempotent, because two routes can both call for it.
    static void emit_unspliced_hypothesis(RawResolveResult& cr) {
        for (int32_t h = 0; h < cr.n_gap_hypotheses(); ++h) {
            if (cr.gap_intron_offsets[h] == cr.gap_intron_offsets[h + 1]) return;
        }
        cr.gap_intron_offsets.push_back(static_cast<int32_t>(cr.gap_introns.size()));
        cr.gap_sj_strand.push_back(STRAND_NONE);
    }

public:

    /// The fragment's length on transcript ``t_idx``: its genomic span ``[gstart, gend)`` less every
    /// intron of the transcript it JUMPS, one lying strictly inside the span (``gstart < intron start``
    /// and ``intron end < gend``). An intron holding an endpoint is not jumped: the read's bases in it are
    /// the molecule's own (overhang) and count, so a start in an intron is measured back from the NEXT
    /// exon's start and an end forward from the PREVIOUS exon's end. This is the accumulator's L for the
    /// path the transcript implies — the span with the path's introns cut out.
    ///
    /// Strand-agnostic, and at least 1 whenever ``gend > gstart``: a jumped intron has a base of the
    /// fragment on each side. Held to a base-by-base count in ``tests/test_transcript_space_fl.py``.
    ///
    /// ⚠ An intron strictly inside the span is jumped whether or not a read covered it: only the two
    /// endpoints are read, so a read-through of a whole intron is not modelled.
    inline int32_t tx_frag_length(int32_t gstart, int32_t gend, int32_t t_idx) const {
        const int32_t begin   = exon_offsets_[t_idx];
        const int32_t n_exons = exon_offsets_[t_idx + 1] - begin;
        const int32_t span    = gend - gstart;
        if (n_exons < 2) return span;

        const int32_t* starts = exon_starts_.data() + begin;
        const int32_t* ends   = exon_ends_.data()   + begin;
        const int32_t* cumsum = exon_cumsum_.data()  + begin;

        // Intron k is [ends[k], starts[k+1]); it is jumped iff ends[k] > gstart and starts[k+1] < gend,
        // so the jumped introns are one contiguous run k = lo .. hi.
        const int lo = static_cast<int>(std::upper_bound(ends, ends + n_exons, gstart) - ends);
        const int hi = static_cast<int>(std::lower_bound(starts, starts + n_exons, gend) - starts) - 2;
        if (lo > hi) return span;
        // Their total: the genomic distance from exon lo's start to exon hi+1's start, less the exonic
        // bases between those two starts.
        const int32_t jumped = (starts[hi + 1] - starts[lo]) - (cumsum[hi + 1] - cumsum[lo]);
        return span - jumped;
    }

    // ----------------------------------------------------------------
    // Fragment-length computation
    // ----------------------------------------------------------------

    /// Thread-safe: the fragment's length on each candidate (``tx_frag_length``). Results are aligned
    /// to t_inds: frag_lengths[i] is the length for t_inds[i], or -1.
    void compute_frag_lengths(
        const std::vector<ExonBlock>& exons,
        const std::vector<int32_t>& t_inds,
        std::vector<int32_t>& frag_lengths) const
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
            int32_t fl = tx_frag_length(gstart, gend, t_inds[i]);
            if (fl > 0) frag_lengths[i] = fl;
        }
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

    // ================================================================
    // _resolve_core — shared resolution logic behind resolve_fragment()
    // ================================================================

    /// Thread-safe overload: uses caller-supplied scratch buffers. ``exons`` must be non-empty.
    void _resolve_core(
        const std::vector<ExonBlock>& exons,
        const std::vector<IntronBlock>& introns,
        int32_t genomic_footprint,
        RawResolveResult& cr,
        ResolverScratch& scratch)
    {
        int n_exons = static_cast<int>(exons.size());
        scratch.reset_per_fragment();

        cr.genomic_footprint = genomic_footprint;

        // --- Interchromosomal chimera detection ---
        cr.chimera_type = CHIMERA_NONE;
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
                                   &scratch.buf.ptr, &scratch.buf.cap);

            for (int64_t hi = 0; hi < n; hi++) {
                int64_t idx = scratch.buf.ptr[hi];
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
                    for (int32_t k = 0; k < cnt; k++) {
                        int32_t ti = t_set_data_[off + k];
                        block_exon_t.push_back(ti);
                        if (bp > 0) {
                            scratch.mark_dirty(ti);
                            scratch.t_exon_bp[ti] += bp;
                        }
                    }
                } else if (itype == ITYPE_TRANSCRIPT) {
                    int32_t clo = std::max(bstart, h_start);
                    int32_t chi = std::min(bend, h_end);
                    int32_t bp = (chi > clo) ? (chi - clo) : 0;
                    for (int32_t k = 0; k < cnt; k++) {
                        int32_t ti = t_set_data_[off + k];
                        block_transcript_t.push_back(ti);
                        if (bp > 0) {
                            scratch.mark_dirty(ti);
                            scratch.t_transcript_bp[ti] += bp;
                        }
                    }
                }
            }

            sort_unique(block_exon_t);
            sort_unique(block_transcript_t);
        }

        // --- Intrachromosomal chimera detection ---
        if (!is_interchromosomal) {
            cr.chimera_type = detect_chimera(exons, exon_t_sets, max_fragment_length_);
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
        // they cover both blocks by construction).
        if (!nrna_parent_.empty()) {
            auto add_nrna = [&](int32_t ti) {
                if (ti < 0 ||
                    ti >= static_cast<int32_t>(nrna_parent_.size())) return;
                int32_t n = nrna_parent_[ti];
                if (n >= 0) scratch.nrna_t.push_back(n);
            };
            for (const auto& s : exon_t_sets)
                for (int32_t ti : s) add_nrna(ti);
            for (const auto& s : transcript_t_sets)
                for (int32_t ti : s) add_nrna(ti);

            sort_unique(scratch.nrna_t);
            for (int32_t n : scratch.nrna_t) {
                // Synthetic nRNAs have exactly one exon by construction;
                // annotated nascent-equivs are handled via cgranges and
                // never appear here (nrna_parent_ holds only synthetic
                // parents; a non-synthetic one is installed as -1).
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
        // fragment's sj key forward. Every call site builds a fresh `cr` today.
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
                    // First ANNOTATED sj (introns arrive sorted by
                    // ref/start/end) — the per-sj strand table's key.
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
            // Never empty: `sj_t_sets_vec` keeps only non-empty lookups, and `has_annotated_sj` means one.
            auto sj_merge   = merge_sets(sj_t_sets_vec);

            if (!exon_merge.is_empty()) {
                auto isect = vec_intersect(exon_merge.t_inds,
                                           sj_merge.t_inds);
                cr.t_inds = isect.empty() ? std::move(sj_merge.t_inds)
                                          : std::move(isect);
            } else {
                cr.t_inds = std::move(sj_merge.t_inds);
            }
            cr.merge_criteria = sj_merge.criteria;
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
            // do NOT early-exit; let the fragment flow through
            // the buffer with empty t_inds so calibration can categorize
            // it as INTERGENIC.
            cr.splice_type = has_unannotated_sj ? SPLICE_SPLICED_UNANNOT
                                                : SPLICE_UNSPLICED;
        }

        // SPLICE_ARTIFACT promotion.  When the alignment had
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
        bool seen = false;
        int32_t first_strand = 0;
        bool mixed = false;
        for (int32_t t : cr.t_inds) {
            if (t >= 0 && t < static_cast<int32_t>(t_strand_arr_.size())) {
                int32_t s = t_strand_arr_[t];
                if (!seen) {
                    first_strand = s;
                    seen = true;
                } else if (s != first_strand) {
                    mixed = true;
                    break;
                }
            }
        }
        cr.ambig_strand = mixed ? 1 : 0;

        // --- Fragment lengths ---
        if (cr.chimera_type == CHIMERA_NONE) {
            compute_frag_lengths(exons, cr.t_inds, cr.frag_lengths);
        }

        // --- gap-hypothesis enumeration ---
        //
        // ⭐ EVERY fragment, whatever its splice type: a SPLICED fragment can also carry an intron in its
        // mate gap.
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
        }

        // ⛔ THE SPLICE_IMPLICIT PROMOTION IS UNSPLICED-ONLY.
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
    }

    /// The single-threaded path, on the resolver's own scratch.
    void _resolve_core(
        const std::vector<ExonBlock>& exons,
        const std::vector<IntronBlock>& introns,
        int32_t genomic_footprint,
        RawResolveResult& cr)
    {
        _resolve_core(exons, introns, genomic_footprint, cr, scratch_);
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
        _resolve_core(exons, introns, genomic_footprint, cr);
        return nb::cast(ResolvedFragment::from_core(cr));
    }
};

}  // namespace rigel
