/**
 * accumulator.h — the per-fragment tally built during the single-pass BAM scan.
 *
 *     SPEC:   tests/native/_accumulator_reference.py  — this file must reproduce it BYTE FOR BYTE.
 *             Where the two disagree, the Python file wins.
 *     Design: docs/ACCUMULATOR_DESIGN.md          Plan: docs/IMPLEMENTATION_PLAN.md §3
 *
 * THE MODEL
 *   The genome is a graph. One Accumulator holds ONE reference, described by its sorted CUT positions.
 *   A reference contributing `c` cuts owns `c - 1` NODES and `c - 2` interior LINES, and a line is a
 *   0-bp CONTIGUOUS EDGE between two adjacent nodes:
 *
 *       cuts    0        100       200       600        c = 4
 *       nodes   [  n0  ][   n1   ][   n2   ]            c - 1 = 3
 *       lines            line 1    line 2               c - 2 = 2
 *
 *   A JUNCTION EDGE is a directed donor->acceptor link taken from the annotation. A fragment is a
 *   PATH: its aligned blocks joined across the mate gap, broken by introns.
 *
 *   Nodes count fragments CONTAINED (the whole path fits inside one node) and SPANNING (one segment
 *   covers the node whole); edges count fragments CROSSING. Every object stores a uint32 count and a
 *   three integer sums: count, inv_length_sum (fixed point) and length_sum.
 *
 * WHY THREE SUMS AND NOT ONE
 *   With `placements` the number of admissible start positions -- L at a node, L-1 at a 0-bp line:
 *
 *       E[count]   = rho * E[placements]
 *       E[inv_length_sum] = rho * E[placements * (1/placements)] = rho   <- at an EDGE, exactly
 *       E[length_sum]     = rho * E[placements * L]
 *
 * `inv_length_sum` is deliberately NOT called `density`: it is one at an edge, exactly, and is NOT one
 * at a node, where the opportunity (node - L + 1)+ does not cancel. `length_sum` is the second tilt, and
 * without it the (count, inv_length_sum) pair has determinant mu_g - mu_r and so carries NO information
 * about the gDNA/RNA split whenever the two components share a mean length. See
 * docs/NODE_DENSITY_DERIVATION.md.
 *
 *   The opportunity factor cancels identically at an edge for ANY length distribution, which is why no
 *   divisor and no length model appear there. It does not cancel at a node.
 *
 * TWO STRANDS, AND THEY ARE INDEPENDENT
 *   align_strand   the genomic strand the read ALIGNED to. Every read has one. Selects the column.
 *   sj_strand      a splice junction's strand, from its genomic MOTIF (GT..AG is +, its reverse
 *                  complement CT..AC is -). Spliced reads only. Resolves an intron against the
 *                  annotation, and nothing else.
 *   Comparing them yields sense vs antisense, which is DERIVED and never stored. The old 4-channel axis
 *   collapsed that comparison into one bool named `primary`; that concept is gone.
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

//: Every array has exactly two columns, and they ARE the two genome strands.
//: `Strand` has FOUR values -- OR semantics make POS|NEG == AMBIGUOUS, and NONE means no strand at all
//: -- so only POS and NEG name a column. A fragment carrying neither is REJECTED, never filed under one
//: of them: that would credit it to the wrong strand rather than lose it, which is the error class this
//: single convention exists to delete.
inline constexpr std::size_t kNStrandColumns = 2;

/// Array column for `align_strand`, or -1 if it names no single genome strand.
///
/// ⚠ The -1 is the whole point. `align_strand == STRAND_POS ? 0 : 1` compiles, is shorter, and silently
/// books an undefined strand into the minus column.
inline int strand_column(std::int32_t align_strand) noexcept {
    if (align_strand == STRAND_POS) return 0;
    if (align_strand == STRAND_NEG) return 1;
    return -1;
}

// ============================================================================
// the fixed-point reciprocal-opportunity sum
// ============================================================================

//: Densities accumulate as round(kInvLengthScale / placements) in uint64. Integer addition is associative,
//: so the per-worker merge is bit-identical at any thread count -- which float accumulation is not, and
//: that nondeterminism propagated to a ~2.6 % difference in the calibration output. The scale is 2^32
//: because it holds the quantisation error below float32's own epsilon while leaving ample headroom under
//: the uint64 ceiling at realistic depth.
inline constexpr std::uint64_t kInvLengthScale = 1ull << 32;

/// round(kInvLengthScale / placements), rounding halves AWAY FROM ZERO.
///
/// ⚠ The rounding mode is part of the contract -- byte-identity with the Python reference is undefined
/// without it. `placements` must be positive; every caller guards it (a length-1 molecule cannot cross a
/// 0-bp line, so its edge quantum is 0 rather than a division).
inline std::uint64_t inv_length_quantum(std::int64_t placements) noexcept {
    const auto p = static_cast<std::uint64_t>(placements);
    return (2 * kInvLengthScale + p) / (2 * p);
}

// ============================================================================
// what each object stores
// ============================================================================

/// A node: an interval. Two disjoint populations, each two columns.
struct Node {
    std::uint32_t contained_count[kNStrandColumns];
    std::uint32_t spanning_count[kNStrandColumns];
    std::uint64_t contained_inv_length_sum[kNStrandColumns];
    std::uint64_t spanning_inv_length_sum[kNStrandColumns];
    std::uint64_t contained_length_sum[kNStrandColumns];
    std::uint64_t spanning_length_sum[kNStrandColumns];
};
static_assert(sizeof(Node) == 80, "Node must be 80 bytes with no padding");

/// A contiguous edge: the 0-bp line between two adjacent nodes. `spliced` means the FRAGMENT used an
/// annotated junction somewhere -- not that this line is one. gDNA cannot be spliced, so a spliced
/// crossing is a certified RNA crossing.
struct ContiguousEdge {
    std::uint32_t unspliced_count[kNStrandColumns];
    std::uint32_t spliced_count[kNStrandColumns];
    std::uint64_t unspliced_inv_length_sum[kNStrandColumns];
    std::uint64_t spliced_inv_length_sum[kNStrandColumns];
    std::uint64_t unspliced_length_sum[kNStrandColumns];
    std::uint64_t spliced_length_sum[kNStrandColumns];
};
static_assert(sizeof(ContiguousEdge) == 80, "ContiguousEdge must be 80 bytes with no padding");

/// A junction edge: one exact donor->acceptor jump. Spliced by construction, so there is no unspliced
/// population; and it is not a genomic position, so it carries no structural flags.
struct JunctionEdge {
    std::uint32_t count[kNStrandColumns];
    std::uint64_t inv_length_sum[kNStrandColumns];
    std::uint64_t length_sum[kNStrandColumns];
};
static_assert(sizeof(JunctionEdge) == 40, "JunctionEdge must be 40 bytes with no padding");

// ============================================================================
// the fragment-length pools
// ============================================================================

//: Five pools, each PURE BY CONSTRUCTION. Purity removes the circularity: a length model is fitted from
//: a population known to be one component, so nothing is estimated from the fragments it will explain.
//:
//: There is deliberately NO pool for an exonic contained fragment or a multi-line crossing -- those are
//: gDNA/RNA mixtures, and an impure pool is worse than a missing one.
enum class FragmentPool : std::uint8_t {
    kDnaIntergenic     = 0,  // contained in an intergenic node
    kDnaIntronic       = 1,  // contained in an intronic node
    kDnaIntronExon     = 2,  // crossing exactly one line, flanks {intron, exon} -- a "splash" read
    kDnaIntergenicExon = 3,  // crossing exactly one line, flanks {intergenic, exon}
    kRnaSpliced        = 4,  // using an annotated junction, splice OBSERVED
};
inline constexpr std::size_t kNFragmentPools = 5;

//: Coarse node types, as `signature.coarse_type_array` emits them.
inline constexpr std::uint8_t kTypeIntergenic = 0;
inline constexpr std::uint8_t kTypeIntron     = 1;
inline constexpr std::uint8_t kTypeExon       = 2;

// ============================================================================
// the deposit
// ============================================================================

/// Why a fragment did or did not deposit. Every rejection is counted, never silent.
enum class DepositOutcome : std::uint8_t {
    kDeposited       = 0,
    kTooLong         = 1,  // L above the fragment-length limit
    kEmpty           = 2,  // no path left after clipping to the reference
    kStrandUndefined = 3,  // align_strand is NONE or AMBIGUOUS, so it names no column
    kAmbiguousPath   = 4,  // the candidates imply >1 intron set, so L is undetermined (design §9.1)
};

/// The QC counter this outcome increments — and the specification's own key for it, so the two cannot
/// drift apart. Kept beside the enum for the same reason.
inline const char* outcome_key(DepositOutcome outcome) noexcept {
    switch (outcome) {
        case DepositOutcome::kDeposited:       return "deposited";
        case DepositOutcome::kTooLong:         return "dropped_too_long";
        case DepositOutcome::kEmpty:           return "dropped_empty";
        case DepositOutcome::kStrandUndefined: return "dropped_strand_undefined";
        case DepositOutcome::kAmbiguousPath:   return "dropped_ambiguous_path";
    }
    return "";
}

/// One fragment, as the scanner has it. `[start, end)` is the full genomic extent -- leftmost block
/// start to rightmost block end, MATE GAP INCLUDED, because the gap is part of the molecule.
///
/// ⚠ `introns` need not be sorted, disjoint, or de-duplicated; `deposit` normalises them. That is
/// deliberate: a real BAM produces overlapping introns when the mates disagree about an acceptor, and
/// normalising inside is what lets L be DEFINED as the total of the path's segments rather than computed
/// by a second, independent formula that disagrees with it.
struct FragmentPath {
    std::int64_t       start;
    std::int64_t       end;
    const IntronBlock* introns;
    std::size_t        n_introns;
    std::int32_t       align_strand;
    std::int32_t       sj_strand;

    /// SPLICE_IMPLICIT: the splice was never sequenced, so it is barred from the pure-RNA pool. This is a
    /// flag on a DEPOSIT -- such a fragment overlaps an annotated intron and matches in every other way.
    bool               sj_implicit;

    /// ⛔ The candidate transcripts imply DIFFERENT INTRON SETS, so `L`, both quanta, the pool bin and the
    /// set of crossed lines are all undetermined. Deposits on NOTHING and is counted; the second pass
    /// resolves it with the fragment length and the strand. Design §9.1.
    ///
    /// ⚠ Not the same thing as `sj_implicit`, and the accumulator cannot decide it -- only the caller has
    /// the candidate list. It is an outcome here rather than a `return` in the scanner so that the loss is
    /// COUNTED instead of vanishing.
    bool               path_ambiguous;
};

/// Reusable scratch so `deposit` allocates nothing on the per-fragment path.
///
/// ⭐ Measured on the shipped accumulator: the one per-fragment `std::vector` cost 22.8 ns, 18 % of the
/// deposit -- and it is invisible to any profiler that samples by function, because the time is
/// attributed to `malloc`. One instance per worker; the vectors keep their capacity across fragments.
struct DepositScratch {
    std::vector<std::pair<std::int64_t, std::int64_t>> introns;   // normalised: sorted, disjoint, clipped
    std::vector<std::pair<std::int64_t, std::int64_t>> segments;  // the path, introns cut out
    std::vector<std::int32_t>                         sj_ids;     // annotated junction edges used
};

/// The QC denominators. Not optional and not derivable afterwards: every conservation statement
/// downstream must be able to name what it excluded.
struct DepositCounters {
    std::int64_t deposited              = 0;
    std::int64_t dropped_too_long       = 0;
    std::int64_t dropped_empty          = 0;
    std::int64_t dropped_strand_undefined = 0;
    std::int64_t dropped_ambiguous_path = 0;  // >1 implied intron set; the second pass owes these
    std::int64_t unannotated_introns    = 0;  // observed introns with no annotated junction
    std::int64_t contradictory_sj_strand = 0;  // the mates' motif tags disagreed; no splice trusted
    std::int64_t sj_implicit_fragments  = 0;  // SPLICE_IMPLICIT: the splice was not sequenced
    std::int64_t introns_absorbed       = 0;  // overlapping or abutting introns merged away

    void merge_from(const DepositCounters& other) noexcept;
};

// ============================================================================
// Accumulator — one reference
// ============================================================================

class Accumulator {
public:
    /// `cuts` is this reference's sorted, strictly increasing cut positions and is moved in. Length may
    /// be 0 or 1, which is a reference with no nodes: legal, and it deposits nothing.
    ///
    /// `node_types` is either empty or one coarse type per node; it types the length pools. `max_length`
    /// is the fragment-length limit applied to L and the width of the pool histograms, and must be >= 1
    /// -- at 0 every real fragment would be dropped as too long and the whole tally would be silently
    /// empty.
    explicit Accumulator(std::vector<std::int64_t> cuts,
                         std::vector<std::uint8_t> node_types,
                         int max_length);

    /// Install this reference's junction edges as a CSR keyed by DONOR CUT INDEX -- the index the
    /// deposit already computes while locating the lines its path crosses.
    ///
    /// ⚠ The junction-edge id IS the slot: `sj_acceptor_cut[k]` and the bank entry `k` are the same k.
    /// There is no indirection to a row in `edges.feather`; using that row as a bank index writes past
    /// the end of a 404,168-entry array, because the highest such row is 1,447,755.
    ///
    /// ⚠ Slot ORDER is part of the contract, because the id is the rank: the caller must sort on
    /// (donor cut, acceptor cut, sj_strand), matching `Partition.from_cuts` in the Python spec.
    void set_junctions(std::vector<std::int32_t> offsets,       // size n_cuts + 1
                       std::vector<std::int32_t> acceptor_cut,  // acceptor CUT INDEX, not a coordinate
                       std::vector<std::int8_t>  sj_strand);    // the junction's ANNOTATED strand

    std::size_t n_nodes()    const noexcept { return nodes_.size(); }
    std::size_t n_edges()    const noexcept { return edges_.size(); }
    std::size_t n_junctions() const noexcept { return junctions_.size(); }
    std::size_t n_cuts()     const noexcept { return cuts_.size(); }

    const std::int64_t* cuts_data() const noexcept { return cuts_.data(); }
    Node*               nodes_data()      noexcept { return nodes_.data(); }
    const Node*         nodes_data() const noexcept { return nodes_.data(); }
    ContiguousEdge*     edges_data()      noexcept { return edges_.data(); }
    const ContiguousEdge* edges_data() const noexcept { return edges_.data(); }
    JunctionEdge*       junctions_data()      noexcept { return junctions_.data(); }
    const JunctionEdge* junctions_data() const noexcept { return junctions_.data(); }

    /// One uint32 per node counting fragments whose FIRST COVERED BASE lies in it.
    ///
    /// ⭐ This is the accumulator's one real invariant: `sum(node_start_count) == deposited`, checkable
    /// against a number the scanner knows independently. The three "conservation identities" it replaced
    /// were tautologies -- each right-hand side could only be evaluated by re-running the deposit, so a
    /// deliberately broken replay satisfied all three while 91 % of the crossings were junk.
    const std::uint32_t* node_start_count_data() const noexcept { return node_start_count_.data(); }

    /// Length histograms, pool-major: pool p occupies [p*(max_length+1), (p+1)*(max_length+1)), binned
    /// at L. Empty when this reference has no node types.
    const std::int64_t* pool_lengths_data() const noexcept { return pool_lengths_.data(); }
    std::size_t         pool_lengths_size() const noexcept { return pool_lengths_.size(); }
    int                 max_length()        const noexcept { return max_length_; }

    const DepositCounters& counters() const noexcept { return counters_; }

    /// Index of the node containing `position`, clamped into [0, n_nodes - 1].
    ///
    /// ⚠ Clamped, not -1 on miss: `deposit` has already clipped the path into this reference, so a
    /// position outside cannot arrive, and clamping keeps this byte-compatible with the spec's
    /// `min(max(searchsorted(cuts, p, 'right') - 1, 0), n_cuts - 2)`.
    std::int64_t node_of_pos(std::int64_t position) const noexcept;

    /// Deposit one fragment. Allocates nothing: `scratch` is reused across calls.
    DepositOutcome deposit(const FragmentPath& path, DepositScratch& scratch);

    /// Element-wise sum of `other` into this accumulator. Requires identical cut positions.
    void merge_from(const Accumulator& other);

private:
    /// The one length pool this fragment belongs to, or -1 for none.
    std::int64_t fragment_pool(bool spliced,
                               bool sj_implicit,
                               std::int64_t contained_node,
                               std::int64_t sole_line) const noexcept;

    /// The annotated junction-edge id for one intron, or -1 if it is not an annotated junction.
    std::int64_t sj_edge_id(std::int64_t intron_start,
                            std::int64_t intron_end,
                            std::int32_t sj_strand) const noexcept;

    /// The flat cut index of `position`, or -1 if it is not a cut on this reference.
    std::int64_t exact_cut(std::int64_t position) const noexcept;

    std::vector<std::int64_t>  cuts_;              // n_cuts, strictly increasing
    std::vector<Node>          nodes_;             // n_cuts - 1
    std::vector<ContiguousEdge> edges_;            // n_cuts - 2, the interior lines
    std::vector<JunctionEdge>  junctions_;         // one per annotated junction on this reference
    std::vector<std::uint32_t> node_start_count_;  // n_nodes -- its own array, so Node stays 48 B

    std::vector<std::int32_t>  sj_offsets_;        // n_cuts + 1, CSR over the donor cut index
    std::vector<std::int32_t>  sj_acceptor_cut_;   // n_junctions
    std::vector<std::int8_t>   sj_strand_;         // n_junctions, the ANNOTATED strand

    std::vector<std::uint8_t>  node_types_;        // n_nodes, or empty (no pools)
    int                        max_length_ = 0;
    std::vector<std::int64_t>  pool_lengths_;      // kNFragmentPools * (max_length + 1), or empty
    DepositCounters            counters_;
};

// ============================================================================
// AccumulatorSet — one Accumulator per reference over a flat partition
// ============================================================================
//
// `cut_positions` is the concatenated, reference-major cut array; reference f owns
// cut_positions[ref_cut_offsets[f] .. ref_cut_offsets[f+1]). A reference with fewer than 2 cuts owns no
// nodes and no edges, which is legal.
//
class AccumulatorSet {
public:
    AccumulatorSet(const std::int64_t* cut_positions,
                   std::size_t n_positions,
                   const std::int64_t* ref_cut_offsets,
                   std::size_t n_refs,
                   const std::uint8_t* node_types,
                   std::size_t n_node_types,
                   int max_length);

    std::size_t n_refs() const noexcept { return accs_.size(); }

    /// Install the junction CSR for every reference at once, from the FLAT arrays the index emits
    /// (`build_junction_edge_arrays`), keyed by the flat cut index.
    ///
    /// The slicing is pure arithmetic and lives here so it exists once. Reference `f` owns cuts
    /// `[c0, c1)`, and because the CSR is sorted by flat donor cut while references are cut-major, its
    /// junctions are the contiguous SLOT range `[offsets[c0], offsets[c1])`. So per reference:
    ///
    ///     offsets      -> offsets[c0 .. c1]      - offsets[c0]      (length n_cuts + 1)
    ///     acceptor_cut -> acceptor_cut[j0 .. j1] - c0               (a ref-local cut index)
    ///
    /// ⚠ Two consequences, both load-bearing. A reference's junction-edge ids are `slot - j0`, so the
    /// payload's junction axis is exactly the flat slot order concatenated in reference order — which is
    /// what lets `edges_df.edge_row` stay a join key that never crosses the ABI. And the narrowing to
    /// int32 is safe by census: 1.04 M cuts and 404,168 junctions at human scale.
    void set_junctions(const std::int64_t* offsets,       // n_cuts_total + 1, over the FLAT cut axis
                       std::size_t n_offsets,
                       const std::int64_t* acceptor_cut,  // n_junctions_total, FLAT cut indices
                       const std::int8_t* sj_strand,
                       std::size_t n_junctions,
                       const std::int64_t* ref_cut_offsets);

    Accumulator&       at(std::int32_t ref_id);
    const Accumulator& at(std::int32_t ref_id) const;

    void merge_from(const AccumulatorSet& other);

private:
    std::vector<Accumulator> accs_;
};

}  // namespace rigel::accumulator
