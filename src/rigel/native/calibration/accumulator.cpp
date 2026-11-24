/**
 * accumulator.cpp — the deposit rule. Specification: tests/native/_accumulator_reference.py, byte for byte; where
 * the two disagree, the Python file is right. Every ordering decision in `deposit` is load-bearing, and the
 * comments say why: read them before reordering anything.
 */
#include "calibration/accumulator.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

#include "path_length.h"

namespace rigel::accumulator {

void DepositCounters::merge_from(const DepositCounters& other) noexcept {
    deposited                    += other.deposited;
    dropped_too_long             += other.dropped_too_long;
    dropped_empty                += other.dropped_empty;
    dropped_strand_undefined     += other.dropped_strand_undefined;
    deferred_undetermined_gap    += other.deferred_undetermined_gap;
    deferred_multiple_placements += other.deferred_multiple_placements;
}

void GapCensus::merge_from(const GapCensus& other) noexcept {
    resolved_spliced       += other.resolved_spliced;
    deferred_rna_or_gdna   += other.deferred_rna_or_gdna;
    deferred_which_introns += other.deferred_which_introns;
    deferred_both          += other.deferred_both;
}

// ============================================================================
// construction
// ============================================================================

Accumulator::Accumulator(std::vector<std::int64_t> region_bounds, std::vector<std::uint8_t> region_types,
                         int max_length, std::int32_t ref_id)
    : region_bounds_(std::move(region_bounds)), ref_id_(ref_id), max_length_(max_length) {
    if (ref_id_ < 0) {
        throw std::invalid_argument(
            "accumulator: ref_id must be >= 0, got " + std::to_string(ref_id_) +
            ". It is stamped into every deferred record and the second pass replays those onto that "
            "reference's region_bound axis, so there is no such thing as a fragment held for no reference.");
    }
    if (max_length_ < 1) {
        throw std::invalid_argument(
            "accumulator: max_length must be >= 1, got " + std::to_string(max_length_) +
            ". It is the fragment-length limit applied to L as well as the pool-histogram width, so at 0 "
            "every real fragment is dropped as too long and the whole tally is silently empty.");
    }
    for (std::size_t i = 1; i < region_bounds_.size(); ++i) {
        if (region_bounds_[i] <= region_bounds_[i - 1]) {
            throw std::invalid_argument(
                "accumulator: region_bounds must strictly increase, but region_bounds[" + std::to_string(i) + "] = " +
                std::to_string(region_bounds_[i]) + " <= region_bounds[" + std::to_string(i - 1) + "] = " +
                std::to_string(region_bounds_[i - 1]));
        }
    }

    // c region_bounds own c-1 regions and c-2 boundaries; fewer than two own neither, which is legal.
    const std::size_t n_regions = region_bounds_.size() >= 2 ? region_bounds_.size() - 1 : 0;
    const std::size_t n_boundaries = region_bounds_.size() >= 2 ? region_bounds_.size() - 2 : 0;
    regions_.assign(n_regions, Region{});
    boundaries_.assign(n_boundaries, Boundary{});
    region_start_count_.assign(n_regions * kNStrandColumns, 0u);

    // Required whenever this reference owns a region: the specification always has one type per region, and an
    // empty array would silently emit no length pools.
    if (n_regions > 0 && region_types.empty()) {
        throw std::invalid_argument(
            "accumulator: region_types is empty but this reference has " + std::to_string(n_regions) +
            " regions. It types the fragment-length pools, so an empty array would silently emit no pools "
            "at all rather than fail.");
    }
    if (!region_types.empty()) {
        if (region_types.size() != n_regions) {
            throw std::invalid_argument(
                "accumulator: region_types has " + std::to_string(region_types.size()) +
                " entries but this reference has " + std::to_string(n_regions) + " regions");
        }
        region_types_ = std::move(region_types);
    }
    pool_lengths_.assign(kNFragmentPools * (static_cast<std::size_t>(max_length_) + 1), 0);
    deposited_lengths_.assign(static_cast<std::size_t>(max_length_) + 1, 0u);
}

void Accumulator::set_sj(std::vector<std::int32_t> offsets, std::vector<std::int32_t> boundary_right,
                         std::vector<std::int8_t> sj_strand) {
    if (!offsets.empty() && offsets.size() != region_bounds_.size() + 1) {
        throw std::invalid_argument(
            "accumulator: sj CSR offsets must have length n_region_bounds + 1 = " +
            std::to_string(region_bounds_.size() + 1) + ", got " + std::to_string(offsets.size()));
    }
    if (boundary_right.size() != sj_strand.size()) {
        throw std::invalid_argument(
            "accumulator: sj boundary_right has " + std::to_string(boundary_right.size()) +
            " entries but sj_strand has " + std::to_string(sj_strand.size()));
    }
    if (!offsets.empty() && static_cast<std::size_t>(offsets.back()) != boundary_right.size()) {
        throw std::invalid_argument(
            "accumulator: sj CSR ends at " + std::to_string(offsets.back()) +
            " but there are " + std::to_string(boundary_right.size()) + " sj");
    }
    sj_offsets_      = std::move(offsets);
    sj_boundary_right_ = std::move(boundary_right);
    sj_strand_       = std::move(sj_strand);
    sj_.assign(sj_boundary_right_.size(), SpliceJunction{});
}

// ============================================================================
// locating things on the region_bound axis
// ============================================================================

std::int64_t Accumulator::region_of_pos(std::int64_t position) const noexcept {
    if (regions_.empty()) return -1;
    const auto it = std::upper_bound(region_bounds_.begin(), region_bounds_.end(), position);
    const std::int64_t region = (it - region_bounds_.begin()) - 1;
    return std::min<std::int64_t>(std::max<std::int64_t>(region, 0),
                                  static_cast<std::int64_t>(regions_.size()) - 1);
}

std::int64_t Accumulator::exact_region_bound(std::int64_t position) const noexcept {
    const auto it = std::lower_bound(region_bounds_.begin(), region_bounds_.end(), position);
    if (it == region_bounds_.end() || *it != position) return -1;
    return it - region_bounds_.begin();
}

std::int64_t Accumulator::sj_edge_id(std::int64_t intron_start, std::int64_t intron_end,
                                     std::int32_t sj_strand) const noexcept {
    if (sj_offsets_.empty()) return -1;
    // Both ends of every annotated intron are region_bounds, so an intron whose ends are not is unannotated
    // without consulting the table.
    const std::int64_t donor = exact_region_bound(intron_start);
    if (donor < 0) return -1;
    const std::int64_t acceptor = exact_region_bound(intron_end);
    if (acceptor < 0) return -1;

    // The strand filter applies only to a definite strand. NONE means the aligner wrote no motif tag (STAR writes
    // XS, minimap2 ts, some neither), so demanding one would drop every annotated sj of such a BAM. AMBIGUOUS comes
    // only from `sj_edge_ids` and is treated as NONE: the first annotated sj with these ends answers.
    const bool definite = (sj_strand == STRAND_POS || sj_strand == STRAND_NEG);
    const std::int32_t lo = sj_offsets_[static_cast<std::size_t>(donor)];
    const std::int32_t hi = sj_offsets_[static_cast<std::size_t>(donor) + 1];
    for (std::int32_t k = lo; k < hi; ++k) {  // a donor has only a few sj
        if (sj_boundary_right_[static_cast<std::size_t>(k)] != acceptor) continue;
        if (definite && sj_strand_[static_cast<std::size_t>(k)] != sj_strand) continue;
        return k;
    }
    return -1;
}

void Accumulator::sj_edge_ids(const std::int64_t* starts, const std::int64_t* ends, const std::int32_t* sj_strand,
                              std::size_t n, std::int64_t* out) const noexcept {
    for (std::size_t i = 0; i < n; ++i) out[i] = sj_edge_id(starts[i], ends[i], sj_strand[i]);
}

std::pair<std::int64_t, std::int64_t> Accumulator::crossed_bounds(std::int64_t seg_start,
                                                                  std::int64_t seg_end) const noexcept {
    const std::int64_t first =
        std::upper_bound(region_bounds_.begin(), region_bounds_.end(), seg_start) - region_bounds_.begin();
    const std::int64_t last =
        std::lower_bound(region_bounds_.begin(), region_bounds_.end(), seg_end) - region_bounds_.begin();
    return {first, last};
}

bool Accumulator::clip(const OfferedFragment& fragment, std::int64_t& start, std::int64_t& end) const noexcept {
    if (region_bounds_.size() < 2) return false;
    start = std::max(fragment.start, region_bounds_.front());
    end = std::min(fragment.end, region_bounds_.back());
    return end > start;
}

// ============================================================================
// the deposit
// ============================================================================

std::int64_t Accumulator::hypothesis_length(const OfferedFragment& fragment, const GapHypothesis& hypothesis,
                                            std::int64_t start, std::int64_t end, DepositScratch& scratch) const {
    // The one definition of L, `path_length`, which the resolver also scores candidates by. Observed and implied
    // introns are disjoint (the enumeration never searches a gap the CIGAR explained), so this is their union.
    return path_length(fragment.observed_introns, fragment.n_observed_introns, hypothesis.introns,
                       hypothesis.n_introns, start, end, scratch.introns, scratch.segments);
}

DepositOutcome Accumulator::deposit(const OfferedFragment& fragment, DepositScratch& scratch) {
    // The strand check comes first: a fragment with no genome strand has no column, and it wins over the deferral
    // because every fragment counts once and the second pass resolves which path, not which strand.
    std::int64_t start = 0, end = 0;
    const int column = admit(fragment, &start, &end);
    if (column == -1) {
        ++counters_.dropped_strand_undefined;
        return DepositOutcome::kStrandUndefined;
    }
    if (column == -2) {
        ++counters_.dropped_empty;
        return DepositOutcome::kEmpty;
    }

    // ── arbitration: which hypotheses survive, and is exactly one left? ───────────────────────────
    const OfferedFragment arbitrated = with_unspliced_default(fragment);
    bool any_spliced_hypothesis = false;
    arbitrate(arbitrated, start, end, scratch, /*all_stand_if_empty=*/true, &any_spliced_hypothesis);
    const auto& survivors = scratch.survivors;
    if (any_spliced_hypothesis) record_gap_resolution(arbitrated, survivors);

    if (survivors.size() > 1) {
        deferred_.append(arbitrated, ref_id_, start, end);
        ++counters_.deferred_undetermined_gap;
        return DepositOutcome::kDeferred;
    }
    return deposit_survivor(arbitrated, arbitrated.hypotheses[survivors.front().index], column, start, end,
                            scratch);
}

int Accumulator::admit(const OfferedFragment& fragment, std::int64_t* start, std::int64_t* end) const noexcept {
    const int column = strand_column(fragment.align_strand);
    if (column < 0) return -1;
    // Clip to the reference. L is the clipped length, and the clip precedes arbitration because a hypothesis is
    // filtered on its L.
    return clip(fragment, *start, *end) ? column : -2;
}

void Accumulator::arbitrate(const OfferedFragment& arbitrated, std::int64_t start, std::int64_t end,
                            DepositScratch& scratch, bool all_stand_if_empty, bool* any_spliced) const {
    // A hypothesis implying L > `max_length_` is not a molecule this library contains; for the unspliced one, whose
    // L is the genomic span, this says a fragment spanning more than the limit must be RNA. With
    // `all_stand_if_empty`, a filter that would empty the set leaves every hypothesis standing, and kTooLong then
    // counts the fragment (at once, or at the drain if held). The set's offer applies that clause to the union.
    auto& survivors = scratch.survivors;
    survivors.clear();
    bool any = false;
    for (std::size_t h = 0; h < arbitrated.n_hypotheses; ++h) {
        survivors.push_back({h, hypothesis_length(arbitrated, arbitrated.hypotheses[h], start, end, scratch)});
        any |= !arbitrated.hypotheses[h].is_unspliced();
    }
    const auto too_long = [this](const ScoredHypothesis& s) { return s.length > max_length_; };
    if (!all_stand_if_empty || !std::all_of(survivors.begin(), survivors.end(), too_long))
        survivors.erase(std::remove_if(survivors.begin(), survivors.end(), too_long), survivors.end());
    if (any_spliced != nullptr) *any_spliced = any;
}

DepositOutcome Accumulator::deposit_survivor(const OfferedFragment& fragment, const GapHypothesis& chosen,
                                             const int column, const std::int64_t start, const std::int64_t end,
                                             DepositScratch& scratch) {
    // Re-normalised rather than cached per hypothesis: one code path from a hypothesis to its introns.
    const std::int64_t length = hypothesis_length(fragment, chosen, start, end, scratch);
    if (length <= 0) {
        ++counters_.dropped_empty;
        return DepositOutcome::kEmpty;
    }
    if (length > max_length_) {
        ++counters_.dropped_too_long;
        return DepositOutcome::kTooLong;
    }

    // An implied strand is used only when no motif was sequenced on this fragment: an observed motif is evidence,
    // an implied strand only an inference from the annotation. Which annotated sj the path uses is resolved before
    // the crossings, because `spliced` chooses the bank the crossings land in.
    const std::int32_t sj_strand = fragment.sj_strand == STRAND_NONE ? chosen.sj_strand : fragment.sj_strand;
    const bool spliced = resolve_annotated_sj(sj_strand, scratch);
    const std::size_t col = static_cast<std::size_t>(column);

    // The path's own first and last covered base, not the fragment's extent: introns [(150,480)] over [150,500)
    // put the whole path in [480,500), and the extent would credit a region the path never touches.
    const std::pair<std::int64_t, std::int64_t> covered = covered_bases(scratch.segments);
    const std::int64_t first_region = region_of_pos(covered.first);
    const std::int64_t last_region  = region_of_pos(covered.second);
    region_start_count_[static_cast<std::size_t>(first_region) * kNStrandColumns + col] += 1u;
    // The deposited-length histogram moves with the start count and the deposited counter, so all three describe
    // one set of fragments by construction. `length` is already clipped and within the limit.
    deposited_lengths_[static_cast<std::size_t>(length)] += 1u;
    ++counters_.deposited;

    // 0 at L == 1: a length-1 molecule cannot cross a 0-bp boundary, and 1/(L-1) would divide by zero; an L == 1
    // path on an annotated sj then books a count against density 0, which is correct.
    const double inv_boundary = length >= 2 ? 1.0 / static_cast<double>(length - 1) : 0.0;
    const std::int64_t sole_boundary = deposit_crossings(spliced, col, length, inv_boundary, scratch);
    for (const std::int32_t id : scratch.sj_ids) {
        SpliceJunction& sj = sj_[static_cast<std::size_t>(id)];
        sj.count[col] += 1u;
        sj.inv_length_sum += inv_boundary;
    }

    // ── contained: the whole path lies inside one region ────────────────────────────────────────────
    // Not merely "crossed no boundary": an unannotated intron can swallow every boundary between two blocks,
    // leaving a fragment that crosses nothing yet straddles two regions. It deposits on no object but still
    // increments region_start_count, so the loss is visible.
    std::int64_t contained_region = -1;
    if (!spliced && first_region == last_region) {
        contained_region = first_region;
        regions_[static_cast<std::size_t>(contained_region)].contained_count[col] += 1u;
    }

    const std::int64_t pool = fragment_pool(spliced, contained_region, sole_boundary);
    if (pool >= 0) {
        pool_lengths_[static_cast<std::size_t>(pool) * (static_cast<std::size_t>(max_length_) + 1) +
                      static_cast<std::size_t>(length)] += 1;
    }
    return DepositOutcome::kDeposited;
}

std::int64_t Accumulator::deposit_crossings(bool spliced, std::size_t col, std::int64_t length, double inv_boundary,
                                            const DepositScratch& scratch) {
    // A boundary is crossed iff it lies strictly inside a segment, so per segment the crossed boundaries are a
    // contiguous index range.
    std::int64_t n_crossed = 0;
    std::int64_t sole_boundary = -1;
    for (std::size_t segment = 0; segment < scratch.segments.size(); ++segment) {
        const std::int64_t seg_start = scratch.segments[segment].first;
        const std::int64_t seg_end = scratch.segments[segment].second;
        if (seg_end <= seg_start) continue;  // an empty first or last segment: no base to book, no boundary inside it
        const std::pair<std::int64_t, std::int64_t> crossed = crossed_bounds(seg_start, seg_end);
        const std::int64_t first = crossed.first, last = crossed.second;
        for (std::int64_t bound = first; bound < last; ++bound) {
            // A spliced crossing books no length moment: it is certified RNA, which nothing deconvolves.
            Boundary& boundary = boundaries_[static_cast<std::size_t>(bound - 1)];
            if (spliced) {
                boundary.spliced_count[col] += 1u;
            } else {
                boundary.unspliced_count[col] += 1u;
                boundary.unspliced_inv_length_sum += inv_boundary;
            }
        }
        for_each_mass_share(segment, seg_start, seg_end, first, last, length, scratch.sj_id_at_gap,
                            [&](std::int8_t kind, std::int64_t id, double share) {
                                if (kind == kPlaceBoundary) {
                                    Boundary& boundary = boundaries_[static_cast<std::size_t>(id)];
                                    if (spliced) boundary.spliced_mass += share;
                                    else         boundary.unspliced_mass += share;
                                } else {
                                    // `count`'s column, so `mass[c] / count[c]` is a per-strand mean.
                                    sj_[static_cast<std::size_t>(id)].mass[col] += share;
                                }
                            });
        if (last > first) {
            sole_boundary = (n_crossed == 0 && last - first == 1) ? first : -1;
            n_crossed += last - first;
        }
    }
    return sole_boundary;
}

bool Accumulator::resolve_annotated_sj(std::int32_t sj_strand, DepositScratch& scratch) const {
    // Per intron position, because the conserved mass needs to know which of a segment's ends is a sj, and a filtered
    // list loses the alignment between intron `i` and the gap after segment `i`.
    scratch.sj_ids.clear();
    scratch.sj_id_at_gap.assign(scratch.introns.size(), -1);
    // The motif tag is read once per record, so AMBIGUOUS means the mates disagreed about one molecule:
    // contradictory evidence, so no splice is trusted and every intron stays -1.
    if (sj_strand == STRAND_AMBIGUOUS) return false;
    for (std::size_t gap = 0; gap < scratch.introns.size(); ++gap) {
        const std::int64_t id = sj_edge_id(scratch.introns[gap].first, scratch.introns[gap].second, sj_strand);
        if (id >= 0) {
            scratch.sj_id_at_gap[gap] = static_cast<std::int32_t>(id);
            scratch.sj_ids.push_back(static_cast<std::int32_t>(id));
        }
    }
    return !scratch.sj_ids.empty();
}

template <class Visit>
void Accumulator::for_each_mass_share(std::size_t segment, std::int64_t seg_start, std::int64_t seg_end,
                                      std::int64_t first, std::int64_t last, std::int64_t length,
                                      const std::vector<std::int32_t>& sj_id_at_gap, Visit&& visit) const {
    const std::int32_t left_sj = segment > 0 ? sj_id_at_gap[segment - 1] : -1;
    const std::int32_t right_sj = segment < sj_id_at_gap.size() ? sj_id_at_gap[segment] : -1;
    const std::int64_t n_slices = last - first + 1;
    for (std::int64_t i = 0; i < n_slices; ++i) {
        const std::int64_t lo = (i == 0) ? seg_start : region_bounds_[static_cast<std::size_t>(first + i - 1)];
        const std::int64_t hi = (i == n_slices - 1) ? seg_end : region_bounds_[static_cast<std::size_t>(first + i)];
        const std::int64_t left_bound  = (i > 0) ? first + i - 1 : -1;
        const std::int64_t right_bound = (i < n_slices - 1) ? first + i : -1;
        // The segment's own ends are sj bounds and its interior ends boundaries, so each side of a slice has at
        // most one bound.
        const std::int32_t left_jid  = (left_bound  < 0) ? left_sj  : -1;
        const std::int32_t right_jid = (right_bound < 0) ? right_sj : -1;
        const std::int64_t n_bounds = (left_bound >= 0 ? 1 : 0) + (right_bound >= 0 ? 1 : 0)
                                    + (left_jid  >= 0 ? 1 : 0) + (right_jid  >= 0 ? 1 : 0);
        if (n_bounds == 0) continue;
        // A coverage fraction in (0, 1], deposited directly as double.
        const double share = static_cast<double>(hi - lo)
                           / (static_cast<double>(length) * static_cast<double>(n_bounds));
        for (const std::int64_t bound : {left_bound, right_bound})
            if (bound >= 0) visit(kPlaceBoundary, bound - 1, share);
        for (const std::int32_t jid : {left_jid, right_jid})
            if (jid >= 0) visit(kPlaceSj, jid, share);
    }
}

void Accumulator::record_gap_resolution(const OfferedFragment& fragment,
                                        const std::vector<ScoredHypothesis>& survivors) noexcept {
    // The one site, so the census's four classes cannot drift apart. The caller has established that at least one
    // hypothesis is spliced; counting any other fragment would make the umbrella's denominator the whole library.
    std::int64_t n_spliced = 0;
    bool unspliced_survives = false;
    for (const ScoredHypothesis& s : survivors) {
        if (fragment.hypotheses[s.index].is_unspliced()) unspliced_survives = true;
        else ++n_spliced;
    }
    if (survivors.size() == 1) {
        // Necessarily spliced: the unspliced path is the longest, so it is never the last one left.
        ++gap_census_.resolved_spliced;
    } else if (!unspliced_survives) {
        ++gap_census_.deferred_which_introns;   // certified RNA; only the structure is open
    } else if (n_spliced == 1) {
        ++gap_census_.deferred_rna_or_gdna;     // one bit: was anything spliced at all?
    } else {
        ++gap_census_.deferred_both;
    }
}

std::int64_t Accumulator::fragment_pool(bool spliced, std::int64_t contained_region,
                                        std::int64_t sole_boundary) const noexcept {
    // Priority: a splice on the one surviving path is certified RNA; a contained fragment is typed by its
    // region; a single-boundary crossing is a "splash" read typed by its two flanks. Anything else (an exonic
    // contained fragment, a multi-boundary crossing) is a mixture and enters nothing.
    if (spliced) return static_cast<std::int64_t>(FragmentPool::kRnaSpliced);
    if (contained_region >= 0) {
        switch (region_types_[static_cast<std::size_t>(contained_region)]) {
            case kTypeIntergenic: return static_cast<std::int64_t>(FragmentPool::kDnaIntergenic);
            case kTypeIntron:     return static_cast<std::int64_t>(FragmentPool::kDnaIntronic);
            default:              return -1;  // exonic is a gDNA/RNA mixture, absent by design
        }
    }
    if (sole_boundary >= 1) {
        const std::uint8_t left  = region_types_[static_cast<std::size_t>(sole_boundary) - 1];
        const std::uint8_t right = region_types_[static_cast<std::size_t>(sole_boundary)];
        const std::uint8_t lo    = std::min(left, right);
        const std::uint8_t hi    = std::max(left, right);
        if (lo == kTypeIntron && hi == kTypeExon) {
            return static_cast<std::int64_t>(FragmentPool::kDnaIntronExon);
        }
        if (lo == kTypeIntergenic && hi == kTypeExon) {
            return static_cast<std::int64_t>(FragmentPool::kDnaIntergenicExon);
        }
    }
    return -1;
}

std::int64_t Accumulator::length_under(const OfferedFragment& fragment, std::size_t hypothesis_index,
                                       DepositScratch& scratch) const {
    // The clip `deposit` applies, in the same order: L is measured after clipping to the reference.
    std::int64_t start = 0, end = 0;
    if (hypothesis_index >= fragment.n_hypotheses || !clip(fragment, start, end)) return 0;
    return hypothesis_length(fragment, fragment.hypotheses[hypothesis_index], start, end, scratch);
}

void Accumulator::placement_masses(const OfferedFragment& fragment, DepositScratch& scratch,
                                   std::vector<PlaceMass>& out) const {
    // `deposit_survivor`'s clip, sj resolution, crossing and containment on the observed path, writing nothing.
    out.clear();
    std::int64_t start = 0, end = 0;
    if (!clip(fragment, start, end)) return;
    const std::int64_t length = hypothesis_length(fragment, kUnsplicedHypothesis, start, end, scratch);
    if (length <= 0) return;
    const bool spliced = resolve_annotated_sj(fragment.sj_strand, scratch);

    for (std::size_t segment = 0; segment < scratch.segments.size(); ++segment) {
        const std::int64_t seg_start = scratch.segments[segment].first;
        const std::int64_t seg_end = scratch.segments[segment].second;
        if (seg_end <= seg_start) continue;  // an empty first or last segment: no base to book, no boundary inside it
        const std::pair<std::int64_t, std::int64_t> crossed = crossed_bounds(seg_start, seg_end);
        for_each_mass_share(segment, seg_start, seg_end, crossed.first, crossed.second, length, scratch.sj_id_at_gap,
                            [&out](std::int8_t kind, std::int64_t id, double share) {
                                out.push_back({kind, id, share});
                            });
    }
    const std::pair<std::int64_t, std::int64_t> covered = covered_bases(scratch.segments);
    const std::int64_t first_region = region_of_pos(covered.first);
    const std::int64_t last_region  = region_of_pos(covered.second);
    if (!spliced && first_region == last_region) out.push_back({kPlaceRegion, first_region, 1.0});

    // One entry per object: a boundary bounds two slices of one segment.
    std::sort(out.begin(), out.end(), [](const PlaceMass& a, const PlaceMass& b) {
        return a.kind != b.kind ? a.kind < b.kind : a.index < b.index;
    });
    std::size_t w = 0;
    for (std::size_t r = 0; r < out.size(); ++r) {
        if (w > 0 && out[w - 1].kind == out[r].kind && out[w - 1].index == out[r].index) out[w - 1].mass += out[r].mass;
        else out[w++] = out[r];
    }
    out.resize(w);
}

const DeferredFragments& Accumulator::deferred_canonical() {
    deferred_.canonicalise();
    return deferred_;
}

// ============================================================================
// the per-worker merge
// ============================================================================

void Accumulator::check_mergeable(const Accumulator& other) const {
    // Self-merge first: every other check passes against itself, and the held-record bank would append to itself
    // forever.
    if (&other == this) {
        throw std::invalid_argument("accumulator: an accumulator cannot merge itself");
    }
    // `max_length` is the width of both length histograms, and the merge walks `this`'s width, so a narrower
    // `other` would be read past its end.
    if (max_length_ != other.max_length_) {
        throw std::invalid_argument(
            "accumulator: merge_from requires the same max_length (this has " + std::to_string(max_length_) +
            ", other has " + std::to_string(other.max_length_) +
            "); it is the width of the pool and deposited-length histograms");
    }
    // The reference id and the region_bound positions must both match: two references can share an axis, and a
    // size check alone would pass two different axes.
    if (ref_id_ != other.ref_id_) {
        throw std::invalid_argument(
            "accumulator: merge_from requires the same reference (this is ref " +
            std::to_string(ref_id_) + ", other is ref " + std::to_string(other.ref_id_) + ")");
    }
    if (region_bounds_ != other.region_bounds_) {
        throw std::invalid_argument(
            "accumulator: merge_from requires identical region_bound positions (this has " +
            std::to_string(region_bounds_.size()) + ", other has " + std::to_string(other.region_bounds_.size()) + ")");
    }
    if (sj_.size() != other.sj_.size()) {
        throw std::invalid_argument(
            "accumulator: merge_from requires the same sj bank (this has " +
            std::to_string(sj_.size()) + ", other has " +
            std::to_string(other.sj_.size()) + ")");
    }
}

void Accumulator::merge_from(const Accumulator& other) {
    // Every precondition is checked before anything is merged, so a refusal leaves `this` unchanged.
    check_mergeable(other);
    gap_census_.merge_from(other.gap_census_);
    deferred_.merge_from(other.deferred_);

    for (std::size_t i = 0; i < regions_.size(); ++i)
        for (std::size_t c = 0; c < kNStrandColumns; ++c)
            regions_[i].contained_count[c] += other.regions_[i].contained_count[c];
    for (std::size_t i = 0; i < region_start_count_.size(); ++i) region_start_count_[i] += other.region_start_count_[i];
    for (std::size_t i = 0; i < deposited_lengths_.size(); ++i) deposited_lengths_[i] += other.deposited_lengths_[i];
    for (std::size_t i = 0; i < boundaries_.size(); ++i) {
        for (std::size_t c = 0; c < kNStrandColumns; ++c) {
            boundaries_[i].unspliced_count[c] += other.boundaries_[i].unspliced_count[c];
            boundaries_[i].spliced_count[c]   += other.boundaries_[i].spliced_count[c];
        }
        // Outside the column loop: one value per boundary, not per strand.
        boundaries_[i].unspliced_inv_length_sum += other.boundaries_[i].unspliced_inv_length_sum;
        boundaries_[i].unspliced_mass += other.boundaries_[i].unspliced_mass;
        boundaries_[i].spliced_mass   += other.boundaries_[i].spliced_mass;
    }
    for (std::size_t i = 0; i < sj_.size(); ++i) {
        for (std::size_t c = 0; c < kNStrandColumns; ++c) {
            sj_[i].count[c] += other.sj_[i].count[c];
            // Inside the column loop: the sj mass is the one mass with a strand (`SpliceJunction::mass`).
            sj_[i].mass[c]  += other.sj_[i].mass[c];
        }
        sj_[i].inv_length_sum += other.sj_[i].inv_length_sum;
    }
    for (std::size_t i = 0; i < pool_lengths_.size(); ++i) pool_lengths_[i] += other.pool_lengths_[i];
    counters_.merge_from(other.counters_);
}

}  // namespace rigel::accumulator
