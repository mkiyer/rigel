/**
 * accumulator.h — the per-fragment tally built during the single-pass BAM scan. Specification:
 * tests/native/_accumulator_reference.py, reproduced byte for byte; where the two disagree, the Python file is right.
 *
 * One Accumulator holds one reference, described by its sorted region_bound positions. `c` region_bounds give
 * `c - 1` regions and `c - 2` interior boundaries; a boundary is the 0-bp junction between two adjacent regions:
 *
 *       region_bounds    0        100       200       600        c = 4
 *       regions   [  n0  ][   n1   ][   n2   ]            c - 1 = 3
 *       boundaries            boundary 1    boundary 2               c - 2 = 2
 *
 * An sj is an annotated link from an intron's low end to its high end. A fragment is a path: its aligned blocks
 * joined across the mate gap, broken by introns. Regions count the fragments contained in them, boundaries the
 * fragments crossing them. Each object stores a count, the reciprocal-opportunity sum and, on boundaries and sj, the
 * conserved mass.
 * Counts are integers and fractions are double: integer banks merge bit-identically across workers; double banks
 * are re-associated by the merge and agree to rounding.
 *
 * With A(w) the number of admissible start positions (w - 1 across a 0-bp boundary) and the deposit 1/A(w):
 *
 *       E[count]      = rho * E_f[A(w)]
 *       E[sum 1/A]    = rho * P(A > 0)     = rho at a boundary, where P(w >= 2) = 1
 *
 * Two independent strands: `align_strand`, the genomic strand a read aligned to, selects the column; `sj_strand`,
 * a junction's strand from its motif (GT..AG is +, CT..AC is -), only resolves an intron against the annotation.
 * Sense vs antisense is derived from the two and never stored.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "calibration/deferred_fragments.h"
#include "calibration/offered_fragment.h"
#include "constants.h"

namespace rigel::accumulator {

// ============================================================================
// the strand column
// ============================================================================

//: Every array has two columns, one per genome strand. `Strand` has four values (POS|NEG == AMBIGUOUS, NONE is no
//: strand), so only POS and NEG name a column; a fragment with neither is rejected, never filed under one.
inline constexpr std::size_t kNStrandColumns = 2;

/// Array column for `align_strand`, or -1 if it names no single genome strand. The -1 must stay: mapping every
/// non-POS strand to column 1 would book an undefined strand into the minus column.
inline int strand_column(std::int32_t align_strand) noexcept {
    if (align_strand == STRAND_POS) return 0;
    if (align_strand == STRAND_NEG) return 1;
    return -1;
}

// ============================================================================
// what each object stores
// ============================================================================

/// A region: an interval, holding the fragments contained in it, in two strand columns. No spliced fragment is
/// contained: both ends of an annotated intron are region_bounds.
struct Region {
    std::uint32_t contained_count[kNStrandColumns];
};
static_assert(sizeof(Region) == 8, "Region must be 8 bytes with no padding");

/// A contiguous boundary: the 0-bp boundary between two adjacent regions. `spliced` means the fragment used an
/// annotated sj somewhere, not that this boundary is one; gDNA cannot splice, so a spliced crossing is certified RNA.
struct Boundary {
    std::uint32_t unspliced_count[kNStrandColumns];
    std::uint32_t spliced_count[kNStrandColumns];
    /// One value: the aligned strand says nothing about gDNA vs RNA (the counts keep two for the strand model).
    double unspliced_inv_length_sum;
    /// The conserved mass: where `unspliced_count` adds 1 on every boundary a fragment crosses, the mass adds up
    /// to one per fragment across them. One value: turning incidences into a fragment count has no strand in it.
    double unspliced_mass;
    /// The same rule for spliced fragments: a per-boundary certified-RNA term, not a count of spliced fragments.
    double spliced_mass;
};
static_assert(sizeof(Boundary) == 40, "Boundary must be 40 bytes with no padding");

/// An sj: one exact jump from an intron's low end to its high end. Spliced by construction, so it holds no unspliced
/// fragment, and not a genomic position, so it has no structural flags.
struct SpliceJunction {
    std::uint32_t count[kNStrandColumns];
    /// Read by the second pass to score a held fragment's sj evidence.
    double inv_length_sum;
    /// The conserved mass's third axis: a block's end is an sj wherever its intron resolved to an annotated
    /// sj, and each slice of a block shares its `slice_len / L` equally among the objects bounding it
    /// (`Accumulator::deposit`). Gate: tests/native/test_conserved_mass.py. Per strand, in `count`'s column, so
    /// `mass[c] / count[c]` is a per-strand mean: an artifactual junction accumulates on both strands as gDNA does,
    /// and the count alone cannot tell a sj used by many short fragments from one used by few long ones.
    double mass[kNStrandColumns];
};
static_assert(sizeof(SpliceJunction) == 32, "SpliceJunction must be 32 bytes with no padding");

// ============================================================================
// the fragment-length pools
// ============================================================================

//: Five pools, each defined by structure. Only kRnaSpliced is certified (gDNA cannot splice); the four DNA pools
//: are gDNA-dominated, not pure, and `calibration.fl` deconvolves them by a two-pool contrast. Exonic contained
//: fragments and multi-boundary crossings have no pool: they are gDNA/RNA mixtures by structure.
enum class FragmentPool : std::uint8_t {
    kDnaIntergenic     = 0,  // contained in an intergenic region
    kDnaIntronic       = 1,  // contained in an intronic region
    kDnaIntronExon     = 2,  // crossing exactly one boundary, flanks {intron, exon} -- a "splash" read
    kDnaIntergenicExon = 3,  // crossing exactly one boundary, flanks {intergenic, exon}
    kRnaSpliced        = 4,  // using an annotated sj on its one surviving path, splice sequenced or implied
};
inline constexpr std::size_t kNFragmentPools = 5;

//: Coarse region types, as `signature.coarse_type_array` emits them.
inline constexpr std::uint8_t kTypeIntergenic = 0;
inline constexpr std::uint8_t kTypeIntron     = 1;
inline constexpr std::uint8_t kTypeExon       = 2;

// ============================================================================
// the deposit
// ============================================================================

/// Why a fragment did or did not deposit. Every rejection is counted.
enum class DepositOutcome : std::uint8_t {
    kDeposited       = 0,
    kTooLong         = 1,  // L above the fragment-length limit
    kEmpty           = 2,  // no path left after clipping to the reference
    kStrandUndefined = 3,  // align_strand is NONE or AMBIGUOUS, so it names no column
    kDeferred        = 4,  // >1 surviving hypothesis: the gap is undetermined; held whole for the second pass
    kDeferredPlacements = 5,  // >1 placement admits a surviving (placement, hypothesis): held whole, every
                              // placement and hypothesis, for the second pass
};

/// The QC counter this outcome increments, which is also the specification's key for it.
inline const char* outcome_key(DepositOutcome outcome) noexcept {
    switch (outcome) {
        case DepositOutcome::kDeposited:       return "deposited";
        case DepositOutcome::kTooLong:         return "dropped_too_long";
        case DepositOutcome::kEmpty:           return "dropped_empty";
        case DepositOutcome::kStrandUndefined: return "dropped_strand_undefined";
        case DepositOutcome::kDeferred:        return "deferred_undetermined_gap";
        case DepositOutcome::kDeferredPlacements: return "deferred_multiple_placements";
    }
    return "";
}

/// One object a placement deposits on, and its share of the fragment's unit (`Accumulator::placement_masses`):
/// a region containing the placement takes the whole unit; the boundaries and junctions a crossing placement
/// crosses share it by its slices.
inline constexpr std::int8_t kPlaceRegion   = 0;
inline constexpr std::int8_t kPlaceBoundary = 1;
inline constexpr std::int8_t kPlaceSj       = 2;
struct PlaceMass {
    std::int8_t  kind;    // kPlaceRegion, kPlaceBoundary or kPlaceSj
    std::int64_t index;   // ref-local from an `Accumulator`, flat from an `AccumulatorSet`
    double       mass;
};

/// The umbrella census: every fragment whose enumeration produced at least one spliced hypothesis, partitioned by
/// how the gap was resolved. Exhaustive and exclusive, so `sum(GapCensus)` is the umbrella and the three
/// deferred_* sum to `deferred_undetermined_gap`.
///
/// Its own axis, not a `splice_type` (a SPLICED_ANNOT fragment can have an intron in its mate gap too), and it
/// classifies the arbitration, not the deposit: a resolved_* fragment can still be dropped as too long. There is no
/// `resolved_unspliced`: L_spliced <= L_unspliced, so under `L <= max_length` every spliced path survives whenever
/// the unspliced one does.
struct GapCensus {
    std::int64_t resolved_spliced        = 0;  // one survivor, and it necessarily cuts something
    std::int64_t deferred_rna_or_gdna    = 0;  // unspliced vs one spliced path: was anything spliced?
    std::int64_t deferred_which_introns  = 0;  // >= 2 spliced paths, none unspliced: certified RNA
    std::int64_t deferred_both           = 0;  // both questions at once

    void merge_from(const GapCensus& other) noexcept;
};

struct ScoredHypothesis {
    std::size_t  index;     // into OfferedFragment::hypotheses
    std::int64_t length;    // L under this hypothesis
};

/// Reusable scratch so `deposit` allocates nothing on the per-fragment path. One instance per worker; the vectors
/// keep their capacity across fragments.
struct DepositScratch {
    std::vector<std::pair<std::int64_t, std::int64_t>> introns;   // normalised: sorted, disjoint, clipped
    std::vector<std::pair<std::int64_t, std::int64_t>> segments;  // the path, introns cut out
    std::vector<std::int32_t>                         sj_ids;     // annotated splice junctions used
    /// The same resolution per intron position, -1 where unannotated, so the conserved mass can tell which of a
    /// block's ends is a sj: block `i` is bounded by `sj_id_at_gap[i-1]` and `sj_id_at_gap[i]`.
    std::vector<std::int32_t>                         sj_id_at_gap;
    std::vector<ScoredHypothesis>                     survivors;  // arbitration, per fragment
};

/// The QC denominators: every conservation statement downstream must name what it excluded.
struct DepositCounters {
    std::int64_t deposited              = 0;
    std::int64_t dropped_too_long       = 0;
    std::int64_t dropped_empty          = 0;
    std::int64_t dropped_strand_undefined = 0;
    //: >1 surviving hypothesis, so the gap is undetermined. Not dropped: the fragment is held whole in
    //: `DeferredFragments`, and deposited + deferred + dropped_* == offered.
    std::int64_t deferred_undetermined_gap = 0;
    //: >1 placement admits a surviving pair, so which placement is undetermined; held whole in the set's bank.
    //: deposited + deferred_undetermined_gap + deferred_multiple_placements + dropped_* == offered, counting
    //: fragments (one offer per read name), never hits.
    std::int64_t deferred_multiple_placements = 0;

    void merge_from(const DepositCounters& other) noexcept;
};

// ============================================================================
// Accumulator — one reference
// ============================================================================

class Accumulator {
public:
    /// `region_bounds` is this reference's strictly increasing region_bound positions, moved in; 0 or 1 of them is
    /// a reference with no regions, which deposits nothing. `region_types` is empty or one coarse type per region
    /// and types the length pools. `max_length` is the limit on L and the width of the pool histograms, and must
    /// be >= 1. `ref_id` is which reference this accumulator is, with no default: it is stamped into every deferred
    /// record, and the second pass replays those onto that reference's region_bound axis.
    explicit Accumulator(std::vector<std::int64_t> region_bounds, std::vector<std::uint8_t> region_types,
                         int max_length, std::int32_t ref_id);

    /// Install this reference's splice junctions as a CSR keyed by the low end's region_bound index, which the deposit
    /// already computes while locating the boundaries its path crosses.
    ///
    /// The sj id is the slot: `sj_boundary_right[k]` and bank entry `k` are the same k (not a row of
    /// `edges.feather`). Slot order is part of the contract: the caller sorts on (low-end region_bound, high-end
    /// region_bound, sj_strand), matching `Partition.from_region_bounds` in the specification.
    void set_sj(std::vector<std::int32_t> offsets,         // size n_region_bounds + 1
                std::vector<std::int32_t> boundary_right,  // the high end's region_bound index, not a coordinate
                std::vector<std::int8_t>  sj_strand);      // the sj's annotated strand

    std::size_t n_regions()    const noexcept { return regions_.size(); }
    std::size_t n_boundaries() const noexcept { return boundaries_.size(); }
    std::size_t n_sj()         const noexcept { return sj_.size(); }

    const Region*         regions_data()    const noexcept { return regions_.data(); }
    const Boundary*       boundaries_data() const noexcept { return boundaries_.data(); }
    const SpliceJunction* sj_data()         const noexcept { return sj_.data(); }

    /// One uint32 per region and strand column, counting fragments whose first covered base lies in it. Invariant:
    /// `sum(region_start_count) == deposited`, a number the scanner knows independently.
    const std::uint32_t* region_start_count_data() const noexcept { return region_start_count_.data(); }

    /// Length histograms, pool-major: pool p occupies [p*(max_length+1), (p+1)*(max_length+1)), binned at L.
    const std::int64_t* pool_lengths_data() const noexcept { return pool_lengths_.data(); }
    std::size_t         pool_lengths_size() const noexcept { return pool_lengths_.size(); }

    /// Every deposited fragment binned at its L, unconditioned: the anchor for shrinking the conditioned pools, in
    /// their frame. Invariant: `sum(deposited_lengths) == sum(region_start_count) == deposited`.
    const std::uint32_t* deposited_lengths_data() const noexcept { return deposited_lengths_.data(); }
    std::size_t          deposited_lengths_size() const noexcept { return deposited_lengths_.size(); }
    int                  max_length()             const noexcept { return max_length_; }

    const DepositCounters& counters() const noexcept { return counters_; }
    const GapCensus&       gap_census() const noexcept { return gap_census_; }

    /// The deferred queue in its canonical order. Non-const because it sorts; there is no accessor for the append
    /// order, which depends on the worker count.
    const DeferredFragments& deferred_canonical();

    /// Index of the region containing `position`, clamped into [0, n_regions - 1]. `deposit` has already clipped
    /// the path to this reference; the clamp matches the spec's
    /// `min(max(searchsorted(region_bounds, p, 'right') - 1, 0), n_region_bounds - 2)`.
    std::int64_t region_of_pos(std::int64_t position) const noexcept;

    /// Deposit one fragment. Allocates nothing: `scratch` is reused across calls.
    DepositOutcome deposit(const OfferedFragment& fragment, DepositScratch& scratch);

    /// `admit`'s refusals: the fragment has no genome strand, or it clips to nothing.
    static constexpr int ADMIT_NO_STRAND = -1;
    static constexpr int ADMIT_EMPTY = -2;

    /// The strand check and the clip, counting nothing: the strand column with the clipped extent written to
    /// `start`/`end`, or ADMIT_NO_STRAND when the fragment has no genome strand, ADMIT_EMPTY when it clips to nothing.
    /// The first two steps of `deposit`, exposed so the set's offer admits each placement without counting per
    /// placement.
    int admit(const OfferedFragment& fragment, std::int64_t* start, std::int64_t* end) const noexcept;

    /// Every hypothesis's L, and the survivors of the one filter `L <= max_length` into `scratch.survivors` (the
    /// third step of `deposit`). With `all_stand_if_empty`, every hypothesis survives when the filter would empty
    /// the set, as in `deposit`; without it the caller sees the empty set (the set's offer applies that clause
    /// over every placement at once).
    void arbitrate(const OfferedFragment& fragment, std::int64_t start, std::int64_t end,
                   DepositScratch& scratch, bool all_stand_if_empty, bool* any_spliced) const;

    /// Put the umbrella census back: the set's offer deposits a multi-placement fragment's one survivor through
    /// `deposit`, whose census question has no single answer for such a fragment.
    void restore_gap_census(const GapCensus& census) noexcept { gap_census_ = census; }

    /// `L` under one hypothesis, depositing nothing: what the second pass scores against, exposed so the tool has
    /// one definition of fragment length. Returns the clipped `L`, or 0 when the fragment clips away entirely.
    std::int64_t length_under(const OfferedFragment& fragment, std::size_t hypothesis_index,
                              DepositScratch& scratch) const;

    /// The deposit rule's mass at every object the placement's observed path (observed introns cut, no implied
    /// one) deposits on, depositing nothing: what `deposit_survivor` would credit on that path, one entry per
    /// object, the masses summing to 1. Empty where the deposit would credit no object.
    void placement_masses(const OfferedFragment& fragment, DepositScratch& scratch,
                          std::vector<PlaceMass>& out) const;

    /// The annotated sj id for each of `n` introns, or -1: `sj_edge_id` over arrays, exposed for the second pass
    /// so the lookup has one definition. Ids are slots in this reference's CSR, as `set_sj` defines them.
    void sj_edge_ids(const std::int64_t* starts, const std::int64_t* ends, const std::int32_t* sj_strand,
                     std::size_t n, std::int64_t* out) const noexcept;

    /// Element-wise sum of `other` into this accumulator. Requires a different accumulator with the same
    /// `max_length`, reference, region_bound positions and sj count, all checked before anything is merged: a
    /// refusal throws and leaves this accumulator unchanged.
    void merge_from(const Accumulator& other);

private:
    /// The fragment's extent clipped to this reference, into `start`/`end`: false when nothing is left.
    bool clip(const OfferedFragment& fragment, std::int64_t& start, std::int64_t& end) const noexcept;

    /// L for one hypothesis: its implied introns unioned with the observed ones, normalised and clipped into
    /// [start, end). Leaves the normalised introns in `scratch.introns` and the path in `scratch.segments`. The one
    /// definition of L, for the winner and the length filter alike.
    std::int64_t hypothesis_length(const OfferedFragment& fragment, const GapHypothesis& hypothesis,
                                   std::int64_t start, std::int64_t end, DepositScratch& scratch) const;

    /// The umbrella census, at its one site. Called only when at least one hypothesis is spliced.
    void record_gap_resolution(const OfferedFragment& fragment,
                               const std::vector<ScoredHypothesis>& survivors) noexcept;

    /// Write out the one surviving hypothesis: L, the pools and every channel of every object the path touches.
    /// The last step of `deposit`.
    DepositOutcome deposit_survivor(const OfferedFragment& fragment, const GapHypothesis& chosen, int column,
                                    std::int64_t start, std::int64_t end, DepositScratch& scratch);

    /// The crossings and the conserved mass of the path in `scratch`, segment by segment, into the banks; returns the
    /// boundary it crosses when it crosses exactly one, else -1.
    std::int64_t deposit_crossings(bool spliced, std::size_t column, std::int64_t length, double inv_boundary,
                                   const DepositScratch& scratch);

    /// Resolve each of the path's introns (`scratch.introns`) to an annotated sj under `sj_strand`: per intron into
    /// `scratch.sj_id_at_gap` (-1 where unannotated), the hits into `scratch.sj_ids`. True when the path uses one, so
    /// it is spliced.
    bool resolve_annotated_sj(std::int32_t sj_strand, DepositScratch& scratch) const;

    /// One segment of a path that holds a base, `[start, end)`, the `index`-th on the path (the gap after it is
    /// intron `index`), with the region_bounds strictly inside it as the index range `[first, last)`: the boundaries
    /// it crosses are first - 1 .. last - 2.
    struct PathSegment {
        std::size_t  index;
        std::int64_t start, end, first, last;
    };

    /// Each segment of the path in `scratch.segments` that holds a base, as `visit(segment)`. A boundary is crossed
    /// iff it lies strictly inside a segment, so per segment the crossed boundaries are a contiguous index range.
    template <class Visit>
    void for_each_segment(const DepositScratch& scratch, Visit&& visit) const;

    /// Every share of one segment's conserved mass, as `visit(kind, id, share)`: the crossed region_bounds split the
    /// segment into slices, and each slice's `slice_len / length` is shared equally among the objects bounding it, so
    /// the shares over a fragment sum to 1. The segment's own ends are bounds wherever its intron resolved to an
    /// annotated sj, its interior ends boundaries. Coverage-weighted rather than 1/K: both conserve, but only this one
    /// says where the fragment sat. An unbounded segment yields nothing: the contained case, or an unannotated
    /// intron's.
    template <class Visit>
    void for_each_mass_share(const PathSegment& segment, std::int64_t length,
                             const std::vector<std::int32_t>& sj_id_at_gap, Visit&& visit) const;

    /// The one length pool this fragment belongs to, or -1 for none. A fragment arrives here only when exactly one
    /// hypothesis survived, so its L is determined, and an implied splice counts as a sequenced one: barring
    /// inferred lengths would bar fragments whose mates sit far apart, filtering the pool by length.
    std::int64_t fragment_pool(bool spliced, std::int64_t contained_region, std::int64_t sole_boundary) const noexcept;

    /// The annotated sj id for one intron, or -1 if it is not an annotated sj.
    std::int64_t sj_edge_id(std::int64_t intron_start, std::int64_t intron_end, std::int32_t sj_strand) const noexcept;

    /// The flat region_bound index of `position`, or -1 if it is not a region_bound on this reference.
    std::int64_t exact_region_bound(std::int64_t position) const noexcept;

    /// Throws unless `other` can merge into this accumulator (`merge_from`).
    void check_mergeable(const Accumulator& other) const;

    std::vector<std::int64_t>   region_bounds_;       // n_region_bounds, strictly increasing
    std::vector<Region>         regions_;             // n_region_bounds - 1
    std::vector<Boundary>       boundaries_;          // n_region_bounds - 2, the interior boundaries
    std::vector<SpliceJunction> sj_;                  // one per annotated sj on this reference
    // The start-region bank, flat [n_regions * kNStrandColumns] and separate so Region stays 8 bytes: the path's
    // first covered base, by align strand. Its opportunity is the region length for every fragment length, blind
    // to a wall only at the template's downstream end.
    std::vector<std::uint32_t>  region_start_count_;

    std::vector<std::int32_t>   sj_offsets_;          // n_region_bounds + 1, CSR over the low end's region_bound index
    std::vector<std::int32_t>   sj_boundary_right_;   // n_sj
    std::vector<std::int8_t>    sj_strand_;           // n_sj, the annotated strand

    std::vector<std::uint8_t>   region_types_;        // n_regions
    std::int32_t                ref_id_ = 0;          // stamped into every deferred record
    int                         max_length_ = 0;
    std::vector<std::int64_t>   pool_lengths_;        // kNFragmentPools * (max_length + 1)
    std::vector<std::uint32_t>  deposited_lengths_;   // max_length + 1 -- unconditional given deposit
    DepositCounters             counters_;
    GapCensus                   gap_census_;
    DeferredFragments           deferred_;
};

}  // namespace rigel::accumulator
