// accumulator_set.h — one Accumulator per reference over the index's flat partition, and the offer: a fragment with
// every placement it aligned to, deposited, held or rejected once.
#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include "calibration/accumulator.h"

namespace rigel::accumulator {

/// One alignment of a fragment on one reference, with every explanation of its unsequenced gaps: what the scanner
/// offers per hit. A unique mapper offers one; a multimapping fragment one per resolved hit.
struct Placement {
    std::int32_t    ref_id;
    OfferedFragment offered;
};

/// Reusable scratch for `AccumulatorSet::offer`: allocated once per worker, capacity kept across fragments.
struct OfferScratch {
    DeferredFragments        admitted;   // one record per admitted placement, canonicalised, then de-duplicated
    std::vector<std::size_t> kept;       // the de-duplicated placements, in canonical order
    OfferedStorage           storage;    // the materialised placement being arbitrated or deposited
    DepositScratch           deposit;    // the accumulator's own scratch
    std::vector<std::pair<std::size_t, std::size_t>> union_survivors;  // (kept index, hypothesis) within the limit
};

/// `region_bounds` is the concatenated, reference-major region_bound array; reference f owns
/// region_bounds[ref_region_bound_offsets[f] .. ref_region_bound_offsets[f+1]). A reference with fewer than two
/// region_bounds owns no regions and no boundaries, which is legal.
class AccumulatorSet {
public:
    AccumulatorSet(const std::int64_t* region_bounds, std::size_t n_positions,
                   const std::int64_t* ref_region_bound_offsets, std::size_t n_refs, const std::uint8_t* region_types,
                   std::size_t n_region_types, int max_length);

    std::size_t n_refs() const noexcept { return accs_.size(); }

    /// Install the sj CSR for every reference at once, from the flat arrays the index emits (`build_sj_arrays`),
    /// keyed by the flat region_bound index. Reference `f` owns region_bounds `[c0, c1)`; the CSR is sorted by
    /// flat low end and references are region_bound-major, so its sj are the contiguous slot range
    /// `[offsets[c0], offsets[c1])`, and per reference:
    ///
    ///     offsets        -> offsets[c0 .. c1]        - offsets[c0]   (length n_region_bounds + 1)
    ///     boundary_right -> boundary_right[j0 .. j1] - c0            (a ref-local region_bound index)
    ///
    /// A reference's sj ids are `slot - j0`, so the payload's sj axis is the flat slot order concatenated in
    /// reference order, and `edges_df.edge_row` stays a join key that never crosses the ABI. The per-reference
    /// values are int32, which holds at human scale (about 1e6 region_bounds and 4e5 sj); one that does not fit
    /// throws.
    void set_sj(const std::int64_t* offsets,         // n_region_bounds_total + 1, over the flat region_bound axis
                std::size_t n_offsets,
                const std::int64_t* boundary_right,  // n_sj_total, flat region_bound indices
                const std::int8_t* sj_strand, std::size_t n_sj, const std::int64_t* ref_region_bound_offsets);

    Accumulator&       at(std::int32_t ref_id);
    const Accumulator& at(std::int32_t ref_id) const;

    /// `Accumulator::placement_masses` on reference `ref_id`, its objects on the flat axes: regions, boundaries and
    /// sj in reference order, as the payload holds them.
    void placement_masses(std::int32_t ref_id, const OfferedFragment& fragment, DepositScratch& scratch,
                          std::vector<PlaceMass>& out) const;

    /// Offer one fragment with every placement it aligned to; deposit it, hold it or reject it. The specification
    /// is `_accumulator_reference.Accumulator.offer`:
    ///
    ///   * one placement: `at(ref).deposit`, bit for bit, counters included;
    ///   * the strand check and the clip exclude placements one at a time; with none left the fragment is rejected
    ///     once, `dropped_strand_undefined` if every exclusion was the strand, else `dropped_empty`;
    ///   * identical admitted placements are one; with one left, `deposit` at it with all its hypotheses;
    ///   * with several left, the survivors are the (placement, hypothesis) pairs within the length limit, all of
    ///     them if none is; one survivor deposits alone, leaving the umbrella census unchanged; two or more hold the
    ///     fragment whole in this set's bank, `kDeferredPlacements`.
    DepositOutcome offer(const Placement* placements, std::size_t n, OfferScratch& scratch);

    /// The QC counters, the umbrella census and the held-fragment bank over every reference and the set's own
    /// records: what the payload exports. The bank is in canonical order, so it is non-const.
    DepositCounters   total_counters() const;
    GapCensus         total_gap_census() const;
    DeferredFragments total_deferred();

    void merge_from(const AccumulatorSet& other);

private:
    /// The admitted placements, one record each, clipped, into `scratch.admitted`; the fragment's rejection,
    /// counted, when none is admitted.
    std::optional<DepositOutcome> admit_placements(const Placement* placements, std::size_t n, OfferScratch& scratch);
    /// Every (kept placement, hypothesis) pair within the length limit, or every pair when none is, into
    /// `scratch.union_survivors`.
    void arbitrate_union(OfferScratch& scratch) const;

    std::vector<Accumulator> accs_;
    DeferredFragments        deferred_;  // the multi-placement records: a record spans references
    DepositCounters          counters_;  // the fragments the set rejected whole and the ones it held
    // Each reference's first region, boundary and sj on the flat axes.
    std::vector<std::int64_t> region_base_, boundary_base_, sj_base_;
};

}  // namespace rigel::accumulator
