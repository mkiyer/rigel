// accumulator_set.cpp — the per-reference accumulators over the flat partition, and the multi-placement offer.
#include "calibration/accumulator_set.h"

#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace rigel::accumulator {

namespace {

/// ``value`` as the int32 a per-reference sj CSR stores; one that does not fit throws rather than wrapping.
std::int32_t to_int32(std::int64_t value, const char* what) {
    if (value < std::numeric_limits<std::int32_t>::min() || value > std::numeric_limits<std::int32_t>::max()) {
        throw std::invalid_argument(std::string("accumulator set: set_sj's ") + what + " " + std::to_string(value) +
                                    " does not fit the per-reference int32 CSR");
    }
    return static_cast<std::int32_t>(value);
}

}  // namespace

AccumulatorSet::AccumulatorSet(const std::int64_t* region_bounds, std::size_t n_positions,
                               const std::int64_t* ref_region_bound_offsets, std::size_t n_refs,
                               const std::uint8_t* region_types, std::size_t n_region_types, int max_length) {
    if (ref_region_bound_offsets == nullptr ||
        static_cast<std::size_t>(ref_region_bound_offsets[n_refs]) != n_positions) {
        throw std::invalid_argument(
            "accumulator set: ref_region_bound_offsets must end at n_positions = " + std::to_string(n_positions));
    }
    accs_.reserve(n_refs);

    // Reference f's first region and boundary on the flat axes: the counts the earlier references own.
    std::size_t region_base = 0, boundary_base = 0;
    for (std::size_t f = 0; f < n_refs; ++f) {
        const std::size_t lo = static_cast<std::size_t>(ref_region_bound_offsets[f]);
        const std::size_t hi = static_cast<std::size_t>(ref_region_bound_offsets[f + 1]);
        std::vector<std::int64_t> ref_region_bounds(region_bounds + lo, region_bounds + hi);
        const std::size_t n_regions = hi - lo >= 2 ? hi - lo - 1 : 0;

        std::vector<std::uint8_t> types;
        if (region_types != nullptr && n_region_types > 0 && n_regions > 0) {
            if (region_base + n_regions > n_region_types) {
                throw std::invalid_argument(
                    "accumulator set: region_types has " + std::to_string(n_region_types) +
                    " entries but reference " + std::to_string(f) + " needs " +
                    std::to_string(region_base + n_regions));
            }
            types.assign(region_types + region_base, region_types + region_base + n_regions);
        }
        accs_.emplace_back(std::move(ref_region_bounds), std::move(types), max_length,
                           static_cast<std::int32_t>(f));
        region_base_.push_back(static_cast<std::int64_t>(region_base));
        boundary_base_.push_back(static_cast<std::int64_t>(boundary_base));
        sj_base_.push_back(0);
        region_base += n_regions;
        boundary_base += hi - lo >= 2 ? hi - lo - 2 : 0;
    }
}

void AccumulatorSet::set_sj(const std::int64_t* offsets, std::size_t n_offsets, const std::int64_t* boundary_right,
                            const std::int8_t* sj_strand, std::size_t n_sj,
                            const std::int64_t* ref_region_bound_offsets) {
    if (offsets == nullptr || ref_region_bound_offsets == nullptr) {
        throw std::invalid_argument(
            "accumulator set: set_sj needs both the CSR offsets and ref_region_bound_offsets");
    }
    const std::size_t n_region_bounds = static_cast<std::size_t>(ref_region_bound_offsets[accs_.size()]);
    if (n_offsets != n_region_bounds + 1) {
        throw std::invalid_argument(
            "accumulator set: sj CSR offsets must have length n_region_bounds + 1 = " +
            std::to_string(n_region_bounds + 1) + ", got " + std::to_string(n_offsets));
    }
    if (static_cast<std::size_t>(offsets[n_region_bounds]) != n_sj) {
        throw std::invalid_argument(
            "accumulator set: the sj CSR ends at " + std::to_string(offsets[n_region_bounds]) +
            " but " + std::to_string(n_sj) + " sj were given");
    }

    std::vector<std::int32_t> ref_offsets, ref_boundary_right;
    std::vector<std::int8_t>  ref_strand;
    for (std::size_t f = 0; f < accs_.size(); ++f) {
        const std::int64_t c0 = ref_region_bound_offsets[f];
        const std::int64_t c1 = ref_region_bound_offsets[f + 1];
        if (c1 <= c0) continue;  // a reference with no region_bounds owns no regions and no sj
        const std::int64_t j0 = offsets[c0];
        const std::int64_t j1 = offsets[c1];
        sj_base_[f] = j0;

        ref_offsets.clear();
        ref_offsets.reserve(static_cast<std::size_t>(c1 - c0) + 1);
        for (std::int64_t c = c0; c <= c1; ++c) ref_offsets.push_back(to_int32(offsets[c] - j0, "sj offset"));
        ref_boundary_right.clear();
        ref_strand.clear();
        ref_boundary_right.reserve(static_cast<std::size_t>(j1 - j0));
        ref_strand.reserve(static_cast<std::size_t>(j1 - j0));
        for (std::int64_t k = j0; k < j1; ++k) {
            ref_boundary_right.push_back(to_int32(boundary_right[k] - c0, "high-end region_bound"));
            ref_strand.push_back(sj_strand[k]);
        }
        accs_[f].set_sj(ref_offsets, ref_boundary_right, ref_strand);
    }
}

const Accumulator& AccumulatorSet::at(std::int32_t ref_id) const {
    if (ref_id < 0 || static_cast<std::size_t>(ref_id) >= accs_.size()) {
        throw std::out_of_range("accumulator set: ref_id " + std::to_string(ref_id) + " of " +
                                std::to_string(accs_.size()));
    }
    return accs_[static_cast<std::size_t>(ref_id)];
}

Accumulator& AccumulatorSet::at(std::int32_t ref_id) {
    return const_cast<Accumulator&>(std::as_const(*this).at(ref_id));
}

void AccumulatorSet::placement_masses(std::int32_t ref_id, const OfferedFragment& fragment, DepositScratch& scratch,
                                      std::vector<PlaceMass>& out) const {
    at(ref_id).placement_masses(fragment, scratch, out);
    const std::size_t f = static_cast<std::size_t>(ref_id);
    for (auto& m : out) {
        m.index += m.kind == kPlaceRegion ? region_base_[f] : m.kind == kPlaceBoundary ? boundary_base_[f] : sj_base_[f];
    }
}

void AccumulatorSet::merge_from(const AccumulatorSet& other) {
    if (accs_.size() != other.accs_.size()) {
        throw std::invalid_argument(
            "accumulator set: merge_from requires the same reference count (this has " +
            std::to_string(accs_.size()) + ", other has " + std::to_string(other.accs_.size()) + ")");
    }
    for (std::size_t f = 0; f < accs_.size(); ++f) accs_[f].merge_from(other.accs_[f]);
    deferred_.merge_from(other.deferred_);
    counters_.merge_from(other.counters_);
}

DepositCounters AccumulatorSet::total_counters() const {
    DepositCounters total = counters_;
    for (const Accumulator& a : accs_) total.merge_from(a.counters());
    return total;
}

GapCensus AccumulatorSet::total_gap_census() const {
    GapCensus total;
    for (const Accumulator& a : accs_) total.merge_from(a.gap_census());
    return total;
}

DeferredFragments AccumulatorSet::total_deferred() {
    // Each reference's bank is canonical on its own, but the set's records span references, so the union is
    // re-canonicalised.
    DeferredFragments all;
    for (Accumulator& a : accs_) all.merge_from(a.deferred_canonical());
    deferred_.canonicalise();
    all.merge_from(deferred_);
    all.canonicalise();
    return all;
}

// ============================================================================
// the offer
// ============================================================================

DepositOutcome AccumulatorSet::offer(const Placement* placements, std::size_t n, OfferScratch& sc) {
    if (n == 0) throw std::invalid_argument("accumulator set: offer needs at least one placement");
    // One placement is the unique mapper's path, bit for bit.
    if (n == 1) return at(placements[0].ref_id).deposit(placements[0].offered, sc.deposit);

    if (const std::optional<DepositOutcome> rejected = admit_placements(placements, n, sc)) return *rejected;
    // Canonical order (one record per placement), identical placements collapsed to one.
    sc.admitted.canonicalise();
    sc.kept.clear();
    for (std::size_t p = 0; p < sc.admitted.n_placements(); ++p)
        if (p == 0 || sc.admitted.compare_placements(p - 1, p) != 0) sc.kept.push_back(p);
    const auto accumulator_of = [&](std::size_t p) -> Accumulator& {
        return accs_[static_cast<std::size_t>(sc.admitted.ref[p])];
    };

    if (sc.kept.size() == 1) {
        // One placement left: the fragment's one open question, if any, is its gap -- `deposit`'s.
        const std::size_t p = sc.kept[0];
        return accumulator_of(p).deposit(sc.admitted.offered_of(p, 0, sc.admitted.n_hypotheses(p), sc.storage),
                                         sc.deposit);
    }

    arbitrate_union(sc);
    if (sc.union_survivors.size() == 1) {
        // One pair left: deposited at its placement with that hypothesis alone, and the umbrella census restored,
        // since its question has no single answer for this fragment.
        const auto [k, h] = sc.union_survivors[0];
        const std::size_t p = sc.kept[k];
        Accumulator& a = accumulator_of(p);
        const GapCensus before = a.gap_census();
        const DepositOutcome outcome = a.deposit(sc.admitted.offered_of(p, h, h + 1, sc.storage), sc.deposit);
        a.restore_gap_census(before);
        return outcome;
    }
    // Held whole: every kept placement and hypothesis, in canonical order, as one record.
    for (const std::size_t p : sc.kept) deferred_.copy_placement(sc.admitted, p);
    deferred_.close_record();
    ++counters_.deferred_multiple_placements;
    return DepositOutcome::kDeferredPlacements;
}

std::optional<DepositOutcome> AccumulatorSet::admit_placements(const Placement* placements, std::size_t n,
                                                               OfferScratch& sc) {
    std::size_t strand_rejected = 0, empty_rejected = 0;
    sc.admitted.clear();
    for (std::size_t i = 0; i < n; ++i) {
        const Placement& pl = placements[i];
        std::int64_t s = 0, e = 0;
        const int column = at(pl.ref_id).admit(pl.offered, &s, &e);  // an unknown reference is refused, as for one
        if (column == Accumulator::ADMIT_NO_STRAND) { ++strand_rejected; continue; }
        if (column == Accumulator::ADMIT_EMPTY) { ++empty_rejected; continue; }
        sc.admitted.append(with_unspliced_default(pl.offered), pl.ref_id, s, e);  // one record per placement, clipped
    }
    if (sc.admitted.n_records() > 0) return std::nullopt;
    // Every placement excluded: rejected once on this set's counters, as strand-undefined only when the strand was
    // every exclusion.
    if (strand_rejected > 0 && empty_rejected == 0) {
        ++counters_.dropped_strand_undefined;
        return DepositOutcome::kStrandUndefined;
    }
    ++counters_.dropped_empty;
    return DepositOutcome::kEmpty;
}

void AccumulatorSet::arbitrate_union(OfferScratch& sc) const {
    sc.union_survivors.clear();
    for (std::size_t k = 0; k < sc.kept.size(); ++k) {
        const std::size_t p = sc.kept[k];
        const OfferedFragment& whole = sc.admitted.offered_of(p, 0, sc.admitted.n_hypotheses(p), sc.storage);
        accs_[static_cast<std::size_t>(sc.admitted.ref[p])].arbitrate(
            whole, sc.admitted.start[p], sc.admitted.end[p], sc.deposit, /*all_stand_if_empty=*/false, nullptr);
        for (const ScoredHypothesis& sh : sc.deposit.survivors) sc.union_survivors.emplace_back(k, sh.index);
    }
    if (!sc.union_survivors.empty()) return;
    // The filter would empty the union: every pair stands, and the drain's kTooLong counts the chosen one.
    for (std::size_t k = 0; k < sc.kept.size(); ++k)
        for (std::size_t h = 0; h < sc.admitted.n_hypotheses(sc.kept[k]); ++h) sc.union_survivors.emplace_back(k, h);
}

}  // namespace rigel::accumulator
