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
 * A sj boundary is an annotated donor->acceptor link. A fragment is a path: its aligned blocks joined across the
 * mate gap, broken by introns. Regions count the fragments contained in them, boundaries the fragments crossing
 * them. Each object stores a count, the reciprocal-opportunity sum and, on boundaries and sj, the conserved mass.
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

#include "../constants.h"

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

/// A sj boundary: one exact donor->acceptor jump. Spliced by construction, so it holds no unspliced fragment, and
/// not a genomic position, so it has no structural flags.
struct SpliceJunction {
    std::uint32_t count[kNStrandColumns];
    /// Read by the second pass to score a held fragment's sj evidence.
    double inv_length_sum;
    /// The conserved mass's third axis: a block's end is a sj boundary wherever its intron resolved to an annotated
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

/// One hypothesis about what a fragment's unsequenced gaps contain. Each candidate transcript implies one answer
/// (its own introns inside the gaps); transcripts implying the same introns are one hypothesis. The empty
/// hypothesis is the unspliced one, the gap read as real template (gDNA, or RNA unspliced there), so no separate
/// gDNA flag or synthetic-span hypothesis is needed. `introns` are the implied ones only; introns the CIGAR states
/// are cut under every hypothesis and live on `OfferedFragment`.
struct GapHypothesis {
    const IntronBlock*  introns;        // implied; empty => the unspliced (genomic) hypothesis
    std::size_t         n_introns;
    std::int32_t        sj_strand;      // inferred; an observed motif always wins

    bool is_unspliced() const noexcept { return n_introns == 0; }
};

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

/// One fragment offered to the accumulator, with every explanation of its unsequenced gaps; the conservation
/// identity is `deposited + deferred + dropped_* == offered`. `[start, end)` is the full extent, mate gap included.
/// `observed_introns` may overlap or repeat (mates can disagree about an acceptor); `deposit` normalises them, so L
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

/// One rebuilt `OfferedFragment` with the storage its spans point into: the inverse of
/// `DeferredFragments::append_placement`, used by the set's offer.
struct OfferedStorage {
    std::vector<IntronBlock>    observed;
    std::vector<IntronBlock>    implied;
    std::vector<GapHypothesis>  spans;
    OfferedFragment             offered{};
};

/// Fragments held whole for the second pass: more than one surviving explanation of the gap, or of which
/// placement the molecule took. The fragment is stored, not its consequences; the drain replays it through
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

    /// Append one placement, clipped to `start`/`end`, to the record being built; `close_record` ends it. Every
    /// hypothesis is kept, even one the length filter removed: the second pass re-scores with a length law the
    /// first pass lacked.
    void append_placement(const OfferedFragment& fragment,
                          std::int64_t ref_id,
                          std::int64_t start,
                          std::int64_t end);
    void close_record() { placement_offsets.push_back(static_cast<std::int64_t>(start.size())); }

    /// One fragment held on its gap: one placement, one record.
    void append(const OfferedFragment& fragment,
                std::int64_t ref_id,
                std::int64_t clip_start,
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
    std::vector<std::int32_t>                         sj_ids;     // annotated sj boundaries used
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
    explicit Accumulator(std::vector<std::int64_t> region_bounds,
                         std::vector<std::uint8_t> region_types,
                         int max_length,
                         std::int32_t ref_id);

    /// Install this reference's sj boundaries as a CSR keyed by the donor region_bound index, which the deposit
    /// already computes while locating the boundaries its path crosses.
    ///
    /// The sj-boundary id is the slot: `sj_boundary_right[k]` and bank entry `k` are the same k (not a row of
    /// `edges.feather`). Slot order is part of the contract: the caller sorts on (donor region_bound, acceptor
    /// region_bound, sj_strand), matching `Partition.from_region_bounds` in the specification.
    void set_sj(std::vector<std::int32_t> offsets,       // size n_region_bounds + 1
                       std::vector<std::int32_t> boundary_right,  // acceptor region_bound index, not a coordinate
                       std::vector<std::int8_t>  sj_strand);    // the sj's annotated strand

    std::size_t n_regions()    const noexcept { return regions_.size(); }
    std::size_t n_boundaries()    const noexcept { return boundaries_.size(); }
    std::size_t n_sj() const noexcept { return sj_.size(); }

    Region*               regions_data()      noexcept { return regions_.data(); }
    const Region*         regions_data() const noexcept { return regions_.data(); }
    Boundary*     boundaries_data()      noexcept { return boundaries_.data(); }
    const Boundary* boundaries_data() const noexcept { return boundaries_.data(); }
    SpliceJunction*       sj_data()      noexcept { return sj_.data(); }
    const SpliceJunction* sj_data() const noexcept { return sj_.data(); }

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
    int                 max_length()        const noexcept { return max_length_; }

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

    /// The strand check and the clip, counting nothing: the strand column with the clipped extent written to
    /// `start`/`end`, or -1 when the fragment has no genome strand, -2 when it clips to nothing. The first two
    /// steps of `deposit`, exposed so the set's offer admits each placement without counting per placement.
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
    std::int64_t length_under(const OfferedFragment& fragment,
                              std::size_t hypothesis_index,
                              DepositScratch& scratch) const;

    /// The deposit rule's mass at every object the placement's observed path (observed introns cut, no implied
    /// one) deposits on, depositing nothing: what `deposit_survivor` would credit on that path, one entry per
    /// object, the masses summing to 1. Empty where the deposit would credit no object.
    void placement_masses(const OfferedFragment& fragment,
                          DepositScratch& scratch,
                          std::vector<PlaceMass>& out) const;

    /// The annotated sj id for each of `n` introns, or -1: `sj_edge_id` over arrays, exposed for the second pass
    /// so the lookup has one definition. Ids are slots in this reference's CSR, as `set_sj` defines them.
    void sj_edge_ids(const std::int64_t* starts,
                     const std::int64_t* ends,
                     const std::int32_t* sj_strand,
                     std::size_t n,
                     std::int64_t* out) const noexcept;

    /// Element-wise sum of `other` into this accumulator. Requires a different accumulator with the same
    /// `max_length`, reference, region_bound positions and sj count, all checked before anything is merged: a
    /// refusal throws and leaves this accumulator unchanged.
    void merge_from(const Accumulator& other);

private:
    /// The one length pool this fragment belongs to, or -1 for none. A fragment arrives here only when exactly one
    /// hypothesis survived, so its L is determined, and an implied splice counts as a sequenced one: barring
    /// inferred lengths would bar fragments whose mates sit far apart, filtering the pool by length.
    std::int64_t fragment_pool(bool spliced,
                               std::int64_t contained_region,
                               std::int64_t sole_boundary) const noexcept;

    /// L for one hypothesis: its implied introns unioned with the observed ones, normalised and clipped into
    /// [start, end). Leaves the normalised list in `scratch.introns`. The one definition of L, for the winner and
    /// the length filter alike.
    std::int64_t hypothesis_length(const OfferedFragment& fragment,
                                   const GapHypothesis& hypothesis,
                                   std::int64_t start,
                                   std::int64_t end,
                                   DepositScratch& scratch) const;

    /// The umbrella census, at its one site. Called only when at least one hypothesis is spliced.
    void record_gap_resolution(const OfferedFragment& fragment,
                               const std::vector<ScoredHypothesis>& survivors) noexcept;

    /// Write out the one surviving hypothesis: L, the pools and every channel of every object the path touches.
    /// The last step of `deposit`.
    DepositOutcome deposit_survivor(const OfferedFragment& fragment, const GapHypothesis& chosen,
                                    int column, std::int64_t start, std::int64_t end,
                                    DepositScratch& scratch);

    /// The annotated sj-boundary id for one intron, or -1 if it is not an annotated sj.
    std::int64_t sj_edge_id(std::int64_t intron_start,
                            std::int64_t intron_end,
                            std::int32_t sj_strand) const noexcept;

    /// The flat region_bound index of `position`, or -1 if it is not a region_bound on this reference.
    std::int64_t exact_region_bound(std::int64_t position) const noexcept;

    std::vector<std::int64_t>  region_bounds_;              // n_region_bounds, strictly increasing
    std::vector<Region>          regions_;             // n_region_bounds - 1
    std::vector<Boundary> boundaries_;            // n_region_bounds - 2, the interior boundaries
    std::vector<SpliceJunction>  sj_;         // one per annotated sj on this reference
    // The start-region bank, flat [n_regions * kNStrandColumns] and separate so Region stays 8 bytes: the path's
    // first covered base, by align strand. Its opportunity is the region length for every fragment length, blind
    // to a wall only at the template's downstream end.
    std::vector<std::uint32_t> region_start_count_;

    std::vector<std::int32_t>  sj_offsets_;        // n_region_bounds + 1, CSR over the donor region_bound index
    std::vector<std::int32_t>  sj_boundary_right_;   // n_sj
    std::vector<std::int8_t>   sj_strand_;         // n_sj, the annotated strand

    std::vector<std::uint8_t>  region_types_;        // n_regions
    std::int32_t               ref_id_ = 0;        // stamped into every deferred record
    int                        max_length_ = 0;
    std::vector<std::int64_t>  pool_lengths_;      // kNFragmentPools * (max_length + 1)
    std::vector<std::uint32_t> deposited_lengths_; // max_length + 1 -- unconditional given deposit
    DepositCounters            counters_;
    GapCensus                  gap_census_;
    DeferredFragments          deferred_;
};

// ============================================================================
// AccumulatorSet — one Accumulator per reference over a flat partition
// ============================================================================
//
// `region_bounds` is the concatenated, reference-major region_bound array; reference f owns
// region_bounds[ref_region_bound_offsets[f] .. ref_region_bound_offsets[f+1]). A reference with fewer than two
// region_bounds owns no regions and no boundaries, which is legal.

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

class AccumulatorSet {
public:
    AccumulatorSet(const std::int64_t* region_bounds,
                   std::size_t n_positions,
                   const std::int64_t* ref_region_bound_offsets,
                   std::size_t n_refs,
                   const std::uint8_t* region_types,
                   std::size_t n_region_types,
                   int max_length);

    std::size_t n_refs() const noexcept { return accs_.size(); }

    /// Install the sj CSR for every reference at once, from the flat arrays the index emits (`build_sj_arrays`),
    /// keyed by the flat region_bound index. Reference `f` owns region_bounds `[c0, c1)`; the CSR is sorted by
    /// flat donor and references are region_bound-major, so its sj are the contiguous slot range
    /// `[offsets[c0], offsets[c1])`, and per reference:
    ///
    ///     offsets        -> offsets[c0 .. c1]        - offsets[c0]   (length n_region_bounds + 1)
    ///     boundary_right -> boundary_right[j0 .. j1] - c0            (a ref-local region_bound index)
    ///
    /// A reference's sj ids are `slot - j0`, so the payload's sj axis is the flat slot order concatenated in
    /// reference order, and `edges_df.edge_row` stays a join key that never crosses the ABI. The narrowing to int32
    /// holds at human scale (about 1e6 region_bounds and 4e5 sj).
    void set_sj(const std::int64_t* offsets,       // n_region_bounds_total + 1, over the flat region_bound axis
                       std::size_t n_offsets,
                       const std::int64_t* boundary_right,  // n_sj_total, flat region_bound indices
                       const std::int8_t* sj_strand,
                       std::size_t n_sj,
                       const std::int64_t* ref_region_bound_offsets);

    Accumulator&       at(std::int32_t ref_id);
    const Accumulator& at(std::int32_t ref_id) const;

    /// `Accumulator::placement_masses` on reference `ref_id`, its objects on the flat axes: regions, boundaries and
    /// sj in reference order, as the payload holds them.
    void placement_masses(std::int32_t ref_id,
                          const OfferedFragment& fragment,
                          DepositScratch& scratch,
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

    /// This set's own bank of multi-placement records, in canonical order. The export merges it with every
    /// reference's gap bank and canonicalises the union.
    const DeferredFragments& deferred_canonical();

    /// The set's own counters: the fragments it rejected whole and the ones it held on their placements.
    const DepositCounters& counters() const noexcept { return counters_; }

    void merge_from(const AccumulatorSet& other);

private:
    std::vector<Accumulator> accs_;
    DeferredFragments        deferred_;
    DepositCounters          counters_;
    // Each reference's first region, boundary and sj on the flat axes.
    std::vector<std::int64_t> region_base_, boundary_base_, sj_base_;
};

}  // namespace rigel::accumulator
