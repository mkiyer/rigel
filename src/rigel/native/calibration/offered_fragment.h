// offered_fragment.h — one fragment as the scan offers it to the accumulator: its extent, its observed introns, its two
// strands and every explanation of its unsequenced gaps.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "constants.h"
#include "fragment_blocks.h"

namespace rigel::accumulator {

/// One hypothesis about what a fragment's unsequenced gaps contain. Each candidate transcript implies one answer (its
/// own introns inside the gaps); transcripts implying the same introns are one hypothesis. The empty hypothesis is the
/// unspliced one, the gap read as real template (gDNA, or RNA unspliced there), so no separate gDNA flag or
/// synthetic-span hypothesis is needed. `introns` are the implied ones only; introns the CIGAR states are cut under
/// every hypothesis and live on `OfferedFragment`.
struct GapHypothesis {
    const IntronBlock*  introns;        // implied; empty => the unspliced (genomic) hypothesis
    std::size_t         n_introns;
    std::int32_t        sj_strand;      // inferred; an observed motif always wins

    bool is_unspliced() const noexcept { return n_introns == 0; }
};

/// The unspliced hypothesis: cut nothing in the gaps, so the path is the extent less the observed introns.
inline constexpr GapHypothesis kUnsplicedHypothesis{nullptr, 0, STRAND_NONE};

/// One fragment offered to the accumulator, with every explanation of its unsequenced gaps; the conservation
/// identity is `deposited + deferred + dropped_* == offered`. `[start, end)` is the full extent, mate gap included.
/// `observed_introns` may overlap or repeat (mates can disagree about an intron's end); `deposit` normalises them, so L
/// is always the total of the path's segments. `hypotheses` is never empty: with no gap it is the unspliced one.
struct OfferedFragment {
    std::int64_t         start;
    std::int64_t         end;
    const IntronBlock*   observed_introns;   // CIGAR-N: cut under every hypothesis
    std::size_t          n_observed_introns;
    std::int32_t         align_strand;
    std::int32_t         sj_strand;
    const GapHypothesis* hypotheses;
    std::size_t          n_hypotheses;
    /// The hit's sequence over both mates: its edits (NM) and soft-clipped bases. Never deposited; the second
    /// pass's draw prices a multimapper's placements by them, as the EM prices its hits.
    std::int32_t         nm = 0;
    std::int32_t         clip = 0;
};

/// ``fragment``, an empty hypothesis set answered as the unspliced one (the specification's default
/// `UNSPLICED_ONLY`).
inline OfferedFragment with_unspliced_default(const OfferedFragment& fragment) noexcept {
    OfferedFragment out = fragment;
    if (out.n_hypotheses == 0) {
        out.hypotheses = &kUnsplicedHypothesis;
        out.n_hypotheses = 1;
    }
    return out;
}

/// One `OfferedFragment` with the storage its spans point into, filled one way whatever it is rebuilt from: a held
/// bank (`rebuild_offered`) or a fragment marshalled out of Python. Not copyable, since the fragment points into this
/// object's vectors; a move keeps it valid.
class OfferedStorage {
public:
    OfferedStorage() = default;
    OfferedStorage(const OfferedStorage&) = delete;
    OfferedStorage& operator=(const OfferedStorage&) = delete;
    OfferedStorage(OfferedStorage&&) = default;
    OfferedStorage& operator=(OfferedStorage&&) = default;

    /// Forget the last fragment, keeping the capacity.
    void clear() noexcept {
        observed_.clear();
        implied_.clear();
        spans_.clear();
        first_implied_.clear();
    }
    void add_observed(std::int32_t start, std::int32_t end) { observed_.push_back({0, start, end, 0}); }
    /// Open the next hypothesis; `add_implied` fills it.
    void add_hypothesis(std::int32_t sj_strand) {
        first_implied_.push_back(implied_.size());
        spans_.push_back({nullptr, 0, sj_strand});
    }
    void add_implied(std::int32_t start, std::int32_t end) {
        implied_.push_back({0, start, end, 0});
        ++spans_.back().n_introns;
    }
    /// The fragment over everything added since `clear`. Each span records its offset while `implied_` grows and
    /// becomes a pointer only here, so a reallocation during the fill invalidates nothing.
    const OfferedFragment& finish(std::int64_t start, std::int64_t end, std::int32_t align_strand,
                                  std::int32_t sj_strand, std::int32_t nm, std::int32_t clip) {
        for (std::size_t s = 0; s < spans_.size(); ++s) spans_[s].introns = implied_.data() + first_implied_[s];
        offered_ = {start,        end,           observed_.data(), observed_.size(), align_strand,
                    sj_strand,    spans_.data(), spans_.size(),    nm,               clip};
        return offered_;
    }
    const OfferedFragment& offered() const noexcept { return offered_; }

private:
    std::vector<IntronBlock>   observed_;
    std::vector<IntronBlock>   implied_;
    std::vector<GapHypothesis> spans_;
    std::vector<std::size_t>   first_implied_;
    OfferedFragment            offered_{};
};

/// A held bank's flat arrays, read in place: the CSR `DeferredFragments` keeps and the specification's
/// `Tally.deferred_arrays()` names, every array int64 and each intron a (start, end) pair. `nm` and `clip` may be
/// null, read as 0.
struct HeldArrays {
    const std::int64_t* start;
    const std::int64_t* end;
    const std::int64_t* align_strand;
    const std::int64_t* sj_strand;
    const std::int64_t* observed_intron_offsets;
    const std::int64_t* observed_introns;
    const std::int64_t* hypothesis_sj_strand;
    const std::int64_t* hypothesis_intron_offsets;
    const std::int64_t* hypothesis_introns;
    const std::int64_t* nm;
    const std::int64_t* clip;
};

/// Placement `p` of a held bank with its hypotheses `[g0, g1)` (ids on the bank's flat hypothesis axis), rebuilt
/// into `storage`. The bank's coordinates fit int32: the scan wrote them from int32 blocks, and a bank handed in from
/// Python is checked where it is bound.
inline const OfferedFragment& rebuild_offered(const HeldArrays& bank, std::size_t p, std::size_t g0, std::size_t g1,
                                              OfferedStorage& storage) {
    const auto coordinate = [](std::int64_t x) { return static_cast<std::int32_t>(x); };
    storage.clear();
    for (std::int64_t i = bank.observed_intron_offsets[p]; i < bank.observed_intron_offsets[p + 1]; ++i)
        storage.add_observed(coordinate(bank.observed_introns[2 * i]), coordinate(bank.observed_introns[2 * i + 1]));
    for (std::size_t g = g0; g < g1; ++g) {
        storage.add_hypothesis(static_cast<std::int32_t>(bank.hypothesis_sj_strand[g]));
        for (std::int64_t i = bank.hypothesis_intron_offsets[g]; i < bank.hypothesis_intron_offsets[g + 1]; ++i)
            storage.add_implied(coordinate(bank.hypothesis_introns[2 * i]),
                                coordinate(bank.hypothesis_introns[2 * i + 1]));
    }
    return storage.finish(bank.start[p], bank.end[p], static_cast<std::int32_t>(bank.align_strand[p]),
                          static_cast<std::int32_t>(bank.sj_strand[p]),
                          bank.nm != nullptr ? static_cast<std::int32_t>(bank.nm[p]) : 0,
                          bank.clip != nullptr ? static_cast<std::int32_t>(bank.clip[p]) : 0);
}

}  // namespace rigel::accumulator
