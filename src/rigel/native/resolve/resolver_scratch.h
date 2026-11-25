// resolver_scratch.h — the resolver's per-thread working memory. Kept apart from FragmentResolver, which is read-only
// once built, so every scan worker resolves concurrently on its own scratch.
#pragma once

#include <cstdint>
#include <cstdlib>
#include <type_traits>
#include <utility>
#include <vector>

#include "constants.h"
#include "fragment_blocks.h"

namespace rigel {

/// An unsequenced stretch of a fragment on one reference: between two of its blocks, or between its mates.
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

/// The cgranges query buffer and its capacity, owned. `cr_overlap` grows it with `realloc` through the addresses of
/// `ptr` and `cap`, so it is released with `free`. A move hands both over and leaves the source empty, as a fresh
/// buffer is: a source keeping its `cap` with no `ptr` would be written through null by its next query, because
/// cgranges reallocates only once the count reaches `cap`.
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
    // Each transcript's aligned bases on its exons, cleaned sparsely: only the touched transcripts are reset
    std::vector<int32_t> t_exon_bp;
    std::vector<uint8_t> t_touched;
    std::vector<int32_t> touched;

    OverlapQueryBuffer buf;

    // One fragment's resolution: each block's exon-interval and transcript-interval sets, the junctions' holders, the
    // candidates a block overlaps, and the synthetic spans its candidates imply
    std::vector<std::vector<int32_t>> exon_t_sets;
    std::vector<std::vector<int32_t>> transcript_t_sets;
    std::vector<int32_t> sj_holders;
    std::vector<int32_t> all_overlap_t;
    std::vector<int32_t> nascent_t;

    // The gap-hypothesis enumeration: the blocks in order, the unsequenced gaps, one candidate's path and its
    // normalised introns and segments (`path_length`), and the blocks as the hold test reads them
    std::vector<ExonBlock> implicit_blocks;
    std::vector<GapBlock> implicit_gaps;
    std::vector<IntronBlock> hypothesis_path;
    std::vector<std::pair<int64_t, int64_t>> path_introns;
    std::vector<std::pair<int64_t, int64_t>> path_segments;
    std::vector<std::pair<int32_t, int32_t>> held_blocks;

    // The pair rule (`validate_pair`)
    std::vector<std::pair<int32_t, int32_t>> pair_skips;
    std::vector<std::pair<int32_t, int32_t>> pair_aligned;
    std::vector<std::pair<int32_t, int32_t>> pair_covered;
    std::vector<GapBlock> pair_stretches;
    std::vector<IntronBlock> pair_cut;
    std::vector<int32_t> pair_candidates;
    std::vector<int32_t> pair_tmp;

    ResolverScratch() = default;

    explicit ResolverScratch(int32_t n_transcripts) : t_exon_bp(n_transcripts, 0), t_touched(n_transcripts, 0) {
        touched.reserve(TOUCHED_TRANSCRIPTS_RESERVE);
    }

    // Non-copyable; the moves are defaulted, so every member moves, including one added later.
    ResolverScratch(const ResolverScratch&) = delete;
    ResolverScratch& operator=(const ResolverScratch&) = delete;
    ResolverScratch(ResolverScratch&&) = default;
    ResolverScratch& operator=(ResolverScratch&&) = default;

    void add_exon_bp(int32_t t, int32_t bp) {
        if (!t_touched[t]) {
            t_touched[t] = 1;
            touched.push_back(t);
        }
        t_exon_bp[t] += bp;
    }

    void reset_per_fragment() noexcept {
        for (int32_t t : touched) {
            t_exon_bp[t] = 0;
            t_touched[t] = 0;
        }
        touched.clear();
        for (auto& v : exon_t_sets) v.clear();
        for (auto& v : transcript_t_sets) v.clear();
        sj_holders.clear();
        all_overlap_t.clear();
        nascent_t.clear();
        implicit_blocks.clear();
        implicit_gaps.clear();
    }
};

// The defaulted moves take their exception specification from the members: a member whose move can throw fails here,
// at compile time.
static_assert(std::is_nothrow_move_constructible_v<ResolverScratch> &&
                  std::is_nothrow_move_assignable_v<ResolverScratch>,
              "ResolverScratch: every member must move without throwing");

}  // namespace rigel
