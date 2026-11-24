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

/// One rebuilt `OfferedFragment` with the storage its spans point into (`DeferredFragments::offered_of`).
struct OfferedStorage {
    std::vector<IntronBlock>    observed;
    std::vector<IntronBlock>    implied;
    std::vector<GapHypothesis>  spans;
    OfferedFragment             offered{};
};

}  // namespace rigel::accumulator
