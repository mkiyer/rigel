// deferred_fragments.h — the bank of fragments held whole for the second pass, as flat CSR arrays in a canonical order.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "calibration/offered_fragment.h"

namespace rigel::accumulator {

/// Fragments held whole for the second pass: more than one surviving explanation of the gap, or of which placement
/// the molecule took. The fragment is stored, not its consequences; the drain replays it through
/// `Accumulator::deposit`, so there is one deposit path.
///
/// Records hold placements, placements hold observed introns and hypotheses, hypotheses hold implied introns: four
/// cumulative offset arrays starting at 0 (an empty queue is `{0}`). The append order depends on the worker count,
/// so the export sorts by content (`canonicalise`). Every array is int64, strands included, because the
/// specification emits one dtype and the parity gate compares dtypes.
struct DeferredFragments {
    /// Records -> placements. A record is one fragment: a gap-held unique mapper holds one placement, a
    /// multimapping fragment one per admitted alignment, each on its own reference.
    std::vector<std::int64_t> placement_offsets{0};
    std::vector<std::int64_t> ref;                         // which reference: the drain needs the region_bound axis
    std::vector<std::int64_t> start, end;                  // the clipped extent
    std::vector<std::int64_t> align_strand, sj_strand;
    std::vector<std::int64_t> observed_intron_offsets{0};  // in pairs, into observed_introns
    std::vector<std::int64_t> observed_introns;            // flat (start, end)
    std::vector<std::int64_t> hypothesis_offsets{0};       // into the per-hypothesis arrays below
    std::vector<std::int64_t> hypothesis_sj_strand;
    std::vector<std::int64_t> hypothesis_intron_offsets{0};
    std::vector<std::int64_t> hypothesis_introns;
    std::vector<std::int64_t> nm, clip;                    // the hit's sequence, per placement; never deposited

    std::size_t n_records()    const noexcept { return placement_offsets.size() - 1; }
    std::size_t n_placements() const noexcept { return start.size(); }
    std::size_t n_hypotheses(std::size_t p) const noexcept {
        return static_cast<std::size_t>(hypothesis_offsets[p + 1] - hypothesis_offsets[p]);
    }

    /// Append one placement, clipped to `start`/`end`, to the record being built; `close_record` ends it. Every
    /// hypothesis is kept, even one the length filter removed: the second pass re-scores with a length law the
    /// first pass lacked.
    void append_placement(const OfferedFragment& fragment, std::int64_t ref_id, std::int64_t start,
                          std::int64_t end);
    void close_record() { placement_offsets.push_back(static_cast<std::int64_t>(start.size())); }

    /// One fragment held on its gap: one placement, one record.
    void append(const OfferedFragment& fragment, std::int64_t ref_id, std::int64_t clip_start,
                std::int64_t clip_end) {
        append_placement(fragment, ref_id, clip_start, clip_end);
        close_record();
    }

    /// Copy placement `p` of `from` onto the record being built here: every per-placement array, for the canonical
    /// sort and the set's offer.
    void copy_placement(const DeferredFragments& from, std::size_t p);

    /// Rebuild placement `p`'s offered fragment with its hypotheses `[h0, h1)` (indices local to the placement)
    /// into `storage`, which owns every span the result points into.
    const OfferedFragment& offered_of(std::size_t p, std::size_t h0, std::size_t h1,
                                      OfferedStorage& storage) const;

    /// The deposit-content key between two placements of this bank: the specification's
    /// `(ref, start, end, align_strand, sj_strand, observed_introns, hypotheses)`, compared as Python compares
    /// tuples (a prefix sorts first). Returns -1, 0 or 1; placements equal here deposit identically.
    int compare_placements(std::size_t a, std::size_t b) const noexcept;

    /// Forget every record, keeping the capacity (the set's per-fragment scratch).
    void clear();

    /// Concatenate `other`, shifting its offsets. The order is not canonical until `canonicalise` runs.
    void merge_from(const DeferredFragments& other);

    /// Reorder the records into the specification's canonical order: a record's key is the sequence of its
    /// placements' `compare_placements` keys, a prefix sorting first, then their (nm, clip). Records that tie on
    /// both are identical, so the order is total. Idempotent: the export calls it on banks already sorted.
    void canonicalise();

private:
    /// The arrays the offset arrays point into, and the offset arrays themselves (each starting at 0), of a bank.
    template <class Bank>
    static auto value_arrays(Bank& b) noexcept {
        return std::array{&b.ref, &b.start, &b.end, &b.align_strand, &b.sj_strand, &b.observed_introns,
                          &b.hypothesis_sj_strand, &b.hypothesis_introns, &b.nm, &b.clip};
    }
    template <class Bank>
    static auto offset_arrays(Bank& b) noexcept {
        return std::array{&b.placement_offsets, &b.observed_intron_offsets, &b.hypothesis_offsets,
                          &b.hypothesis_intron_offsets};
    }
};

}  // namespace rigel::accumulator
