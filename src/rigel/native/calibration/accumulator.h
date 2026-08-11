/**
 * accumulator.h — the per-fragment tally built during the single-pass BAM scan.
 *
 *     SPEC:   tests/native/_accumulator_reference.py  — this file must reproduce it BYTE FOR BYTE.
 *             Where the two disagree, the Python file wins.
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
 *   Nodes count fragments CONTAINED (the whole path fits inside one node); edges count fragments
 *   CROSSING. Each population stores only the channels something READS, and they differ: count,
 *   inv_length_sum (fixed point), length_sum, and -- on the contiguous edges -- the conserved mass.
 *
 * WHY MORE THAN ONE SUM
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
// ONE NUMERIC CONVENTION
// ============================================================================

//: ⭐⭐⭐ A COUNT IS AN INTEGER. A FRACTION IS double. There is no fixed point in the tally and no
//: scale constant to decode.
//:
//: ⛔ The predecessor accumulated every fraction as round(2^32 / placements) in uint64, because integer
//: addition is associative and therefore bit-identical across worker counts. The argument was sound; the
//: price it quoted was not. The ~2.6 % it cited was measured on a **float32** accumulator (~3.7e-7 per
//: cell). double is ~1e-16 -- 3.4e9x finer -- reaching the deliverable at ~1e-11, five orders below
//: EMConfig.convergence_delta = 1e-6.
//:
//: ⭐⭐ And the fixed point was LESS ACCURATE, measured against exact rational arithmetic on the
//: reciprocal-opportunity theorem (each length contributes exactly one density unit):
//:
//:      node_len 151    fixed 7.0e-10    double 5.8e-15      120,000x better
//:      node_len 400    fixed 1.7e-08    double 1.0e-13      170,000x better
//:      node_len 1000   fixed 2.0e-07    double 2.8e-13      714,000x better
//:
//: ⚠ The exactness the old gates asserted was a property of their FIXTURES: 1/2 + 1/3 + 1/6 lands back
//: on 2^32 because two rounding errors cancel, while 1/3 + 1/3 + 1/3 is one quantum short -- and double
//: is exact on both.
//:
//: ⛔ What is genuinely given up is bit-identity across worker counts, since float addition is not
//: associative. Owner ruling 2026-08-10: one convention, and this is it.

// ============================================================================
// what each object stores
// ============================================================================

/// A node: an interval. ONE population — the fragments contained inside it — in two strand columns.
///
/// ⚠ A `spanning` population (one segment covering the node whole) was removed on evidence: it reached
/// no evidence-starved node the node's own bounding EDGES did not already reach off capture, and 141
/// nodes / 822 fragments (0.008 %) under it. Its mass is not lost — a spanning fragment crosses both of
/// the node's lines and is deposited there.
/// ⛔ Consequence, and it is structural: no spliced fragment touches the node axis at all, because a
/// spliced fragment can never be `contained` (both endpoints of an annotated intron are cuts).
struct Node {
    std::uint32_t contained_count[kNStrandColumns];
    /// ⭐ ONE value, while `contained_count` keeps two. The length moments are strand-AGNOSTIC -- which
    /// strand a read aligned to says nothing about whether the molecule was gDNA or RNA -- and every
    /// consumer summed the two columns before using them. The COUNTS keep both because the strand model
    /// is a Beta-Binomial over them, per strand.
    double contained_inv_opportunity_sum;
    std::uint64_t contained_length_sum;
};
static_assert(sizeof(Node) == 24, "Node must be 24 bytes with no padding");

/// A contiguous edge: the 0-bp line between two adjacent nodes. `spliced` means the FRAGMENT used an
/// annotated junction somewhere -- not that this line is one. gDNA cannot be spliced, so a spliced
/// crossing is a certified RNA crossing.
struct ContiguousEdge {
    std::uint32_t unspliced_count[kNStrandColumns];
    std::uint32_t spliced_count[kNStrandColumns];
    /// ⭐ ONE value each -- strand-agnostic, see `Node`.
    double unspliced_inv_length_sum;
    std::uint64_t unspliced_length_sum;
    /// ⭐⭐ THE CONSERVED MASS, fixed point. A COUNT and a MASS are two different deposits and one
    /// number cannot be both: `unspliced_count` is `+1` on every line a fragment crosses, so a fragment
    /// books `max(K, 1)` of them; this sums to ONE per fragment, across all the lines it crosses.
    ///
    /// ⛔ ONE VALUE, NOT TWO, WHILE EVERY BANK ABOVE IS PER STRAND — deliberate. `strand_deconv` reads
    /// the counts per column; nothing reads a mass per strand, because the mass exists to turn an
    /// object-incidence total into a fragment count and that question has no strand in it.
    double unspliced_mass;
    /// ⭐ The same rule, routed by the same `spliced` flag — so `mass` is not the one channel that
    /// ignores the split. ⛔ A PARTIAL, never a conservation ledger: a spliced fragment's blocks with no
    /// interior line deposit nothing (their accounting is on the junction axis), so this sums to
    /// `crossed_block_len / L`. It is a per-LINE certified-RNA term, commensurate with the unspliced
    /// mass at the same line, and is NOT "the number of spliced fragments here".
    double spliced_mass;
};
static_assert(sizeof(ContiguousEdge) == 48, "ContiguousEdge must be 48 bytes with no padding");

/// A junction edge: one exact donor->acceptor jump. Spliced by construction, so there is no unspliced
/// population; and it is not a genomic position, so it carries no structural flags.
struct JunctionEdge {
    std::uint32_t count[kNStrandColumns];
    /// ⭐ LIVE: `second_pass` scores a held fragment's junction evidence with it. `length_sum` was
    /// removed — nothing read it, and `pool_lengths`' RNA_SPLICED row already carries that
    /// population's length distribution.
    double inv_length_sum;
    /// ⭐⭐⭐ THE CONSERVED MASS'S THIRD AXIS. A spliced fragment's block that contains no interior
    /// line deposits on neither edge bank, and is not `contained` either -- its path spans a junction,
    /// so it lies in no single node. Such a fragment existed on the incidence axis and on no conserved
    /// one, which is why a library fragment count was not computable. Measured on the origin-split
    /// oracle at ladder g50 capture_off: 1,222,375 of 4,830,713 RNA fragments (25.3 %) are in that
    /// population, against 0 of 4,997,761 gDNA fragments, because gDNA cannot splice.
    ///
    /// ⛔ The rule ADDS a boundary class; it does not re-apportion an existing one. A block that
    /// crossed a line is untouched, so `unspliced_mass` and `spliced_mass` are byte-identical to what
    /// they were. Spec: `_accumulator_reference.py`; gates: `tests/native/test_conserved_mass.py`.
    double mass;
};
static_assert(sizeof(JunctionEdge) == 24, "JunctionEdge must be 24 bytes with no padding");

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
    kDeferred        = 4,  // >1 surviving hypothesis: the gap is undetermined, so the fragment
                           // is held WHOLE for the second pass
};

/// The QC counter this outcome increments — and the specification's own key for it, so the two cannot
/// drift apart. Kept beside the enum for the same reason.
inline const char* outcome_key(DepositOutcome outcome) noexcept {
    switch (outcome) {
        case DepositOutcome::kDeposited:       return "deposited";
        case DepositOutcome::kTooLong:         return "dropped_too_long";
        case DepositOutcome::kEmpty:           return "dropped_empty";
        case DepositOutcome::kStrandUndefined: return "dropped_strand_undefined";
        case DepositOutcome::kDeferred:        return "deferred_undetermined_gap";
    }
    return "";
}

/// ONE hypothesis about what a fragment's UNSEQUENCED gaps contain.
///
/// A mate gap may hold no intron, one, or several, and which it is cannot be observed -- the bases are
/// not there. Each candidate transcript determines exactly one answer (its own introns lying inside the
/// gaps), so the hypotheses are finite and small, and two transcripts implying the same introns are ONE
/// hypothesis.
///
/// ⭐ THE EMPTY HYPOTHESIS IS THE UNSPLICED ONE, and it is the genomic explanation: cutting nothing means
/// the gap is real template, i.e. the molecule is gDNA -- or nascent RNA, which is the same unspliced
/// span. That is why the accumulator needs no separate "could this be gDNA?" flag, and why the nascent
/// shadow transcript is not a candidate: it IS this hypothesis.
///
/// ⚠ `introns` are the IMPLIED ones only. Introns the CIGAR actually stated are cut under EVERY
/// hypothesis and live on `OfferedFragment` instead, because they are not in doubt.
struct GapHypothesis {
    const IntronBlock*  introns;        // implied; empty => the unspliced (genomic) hypothesis
    std::size_t         n_introns;
    std::int32_t        sj_strand;      // ⚠ an INFERENCE; an observed motif always wins
    const std::int32_t* supporting_t;   // candidate transcripts implying this path -- the second pass
    std::size_t         n_supporting;   //   weights hypotheses by their abundance; pass one never reads

    bool is_unspliced() const noexcept { return n_introns == 0; }
};

/// One fragment offered to the accumulator, with every explanation of its unsequenced gaps.
///
/// ⭐ Named for the population the conservation identity counts:
/// `deposited + deferred + dropped_* == offered`.
///
/// `[start, end)` is the full genomic extent -- leftmost block start to rightmost block end, MATE GAP
/// INCLUDED, because the gap is part of the molecule.
///
/// ⚠ `observed_introns` need not be sorted, disjoint, or de-duplicated; `deposit` normalises them. That
/// is deliberate: a real BAM produces overlapping introns when the mates disagree about an acceptor, and
/// normalising inside is what lets L be DEFINED as the total of the path's segments rather than computed
/// by a second, independent formula that disagrees with it.
///
/// ⚠ `hypotheses` must never be empty -- a fragment with no unsequenced gap still has ONE hypothesis,
/// the unspliced one. The degenerate case is the general case, not a branch.
struct OfferedFragment {
    std::int64_t         start;
    std::int64_t         end;
    const IntronBlock*   observed_introns;   // CIGAR-N: cut under EVERY hypothesis
    std::size_t          n_observed_introns;
    std::int32_t         align_strand;
    std::int32_t         sj_strand;
    const GapHypothesis* hypotheses;
    std::size_t          n_hypotheses;
};

/// ⭐ The umbrella census (owner ruling, 2026-08-01): every fragment whose enumeration produced at least
/// one non-unspliced hypothesis, partitioned by how the gap was RESOLVED. Exhaustive and mutually
/// exclusive, so `sum(GapCensus) == the umbrella` and the three deferred_* == `deferred_undetermined_gap`.
///
/// ⛔ Its own axis, NOT a `splice_type`. The umbrella cuts ACROSS the splice census: a certified-RNA
/// SPLICED_ANNOT fragment with an intron in its mate gap needs resolving exactly as much as an UNSPLICED
/// one does, so putting these on `splice_type` would need two labels per fragment.
///
/// ⚠ These classify the ARBITRATION, not the deposit: a resolved_* fragment can still be rejected
/// afterwards as TOO_LONG, which is a different question with its own counter.
///
/// ⛔ THERE IS NO `resolved_unspliced`, AND IT IS NOT AN OMISSION. The field existed and no fragment could
/// enter it: a spliced hypothesis CUTS bases the unspliced one keeps, so L_spliced <= L_unspliced always,
/// and the one filter is `L <= max_length`. If the unspliced path survives the filter then every spliced
/// path survives it too, so the survivor set can never be exactly {unspliced} while a spliced path was
/// offered -- which is the condition for being in this census at all. The ORDERING is pinned directly by
/// `test_gap_hypothesis_arbitration.test_the_GENOMIC_hypothesis_is_ALWAYS_the_LONGEST`.
struct GapCensus {
    std::int64_t resolved_spliced        = 0;  // one survivor, and it necessarily cuts something
    std::int64_t deferred_rna_or_gdna    = 0;  // unspliced vs ONE spliced path: was anything spliced?
    std::int64_t deferred_which_introns  = 0;  // >= 2 spliced paths, none unspliced: certified RNA
    std::int64_t deferred_both           = 0;  // both questions at once

    void merge_from(const GapCensus& other) noexcept;
};

/// ⭐ Fragments whose gap has more than one surviving explanation, held WHOLE for the second pass
/// ( calls this the side buffer).
///
/// The FRAGMENT is stored, never its consequences. Object ids are large, derived, and would have to be
/// kept consistent with the partition; the fragment is small and replays exactly. The drain re-enters
/// `Accumulator::deposit` with the chosen hypothesis, so there is no second deposit path, no duplicated
/// crossing logic, and byte-identity with the specification is preserved for free.
///
/// Two nested variable-length levels -- fragments hold hypotheses, hypotheses hold introns -- so there
/// are two offset arrays. Offsets are cumulative and start at 0, so an empty queue is `{0}`, never `{}`.
///
/// ⛔ ORDER IS OBSERVABLE HERE AND NOWHERE ELSE. Every other bank is a sum of integers and integer
/// addition is associative, so a per-worker merge is exact whatever order the chunks arrived in. This is
/// a LIST. Concatenating per-worker queues gives a different byte sequence at 1, 2, 4 and 8 workers with
/// identical contents -- so the EXPORT sorts on the record's own content, exactly as
/// `Tally.deferred_arrays()` does in the specification. Two records that tie are identical records.
///
/// ⚠ EVERY ARRAY IS int64, INCLUDING THE TWO STRAND COLUMNS. They are int32 everywhere else in the
/// scanner, but the specification's flattening emits one dtype and the parity gate compares dtypes -- so
/// widening here costs 8 bytes per deferred record and removes a conversion that would otherwise have to
/// happen at the export, where getting it wrong compares equal by value.
struct DeferredFragments {
    std::vector<std::int64_t> ref;                         // which reference: the drain needs the cut axis
    std::vector<std::int64_t> start, end;                  // the CLIPPED extent
    std::vector<std::int64_t> align_strand, sj_strand;
    std::vector<std::int64_t> observed_intron_offsets{0};  // in PAIRS, into observed_introns
    std::vector<std::int64_t> observed_introns;            // flat (start, end)
    std::vector<std::int64_t> hypothesis_offsets{0};       // into the per-hypothesis arrays below
    std::vector<std::int64_t> hypothesis_sj_strand;
    std::vector<std::int64_t> hypothesis_intron_offsets{0};
    std::vector<std::int64_t> hypothesis_introns;
    std::vector<std::int64_t> hypothesis_t_offsets{0};
    std::vector<std::int64_t> hypothesis_t;

    std::size_t size() const noexcept { return start.size(); }

    /// Append one fragment with every hypothesis it arrived with. `start`/`end` are the CLIPPED extent,
    /// because that is what the drain must replay.
    ///
    /// ⚠ EVERY hypothesis, including any the length filter removed. The record is what was OFFERED: the
    /// second pass re-scores from scratch with a fragment-length distribution the first pass did not have,
    /// so pre-pruning here would decide with the weaker evidence and hide the decision.
    void append(const OfferedFragment& fragment,
                std::int64_t ref_id,
                std::int64_t start,
                std::int64_t end);

    /// Concatenate `other`, shifting its offsets. ⚠ Order is not canonical until `canonicalise` runs.
    void merge_from(const DeferredFragments& other);

    /// ⭐ Reorder the records into the ONE canonical order, which is the specification's:
    ///
    ///     (ref, start, end, align_strand, sj_strand, observed_introns, hypotheses)
    ///
    /// compared exactly as Python compares those tuples -- element-wise, and a PREFIX sorts BEFORE the
    /// longer sequence it is a prefix of. Two records that tie on this key are identical records, so their
    /// relative order cannot be observed and the sort needs no tie-break.
    ///
    /// ⚠ Idempotent, and it must be: the export calls it and a merged accumulator may already be sorted.
    void canonicalise();
};

/// Reusable scratch so `deposit` allocates nothing on the per-fragment path.
///
/// ⭐ Measured on the shipped accumulator: the one per-fragment `std::vector` cost 22.8 ns, 18 % of the
/// deposit -- and it is invisible to any profiler that samples by function, because the time is
/// attributed to `malloc`. One instance per worker; the vectors keep their capacity across fragments.
struct ScoredHypothesis {
    std::size_t  index;     // into OfferedFragment::hypotheses
    std::int64_t length;    // L under this hypothesis
    std::int64_t absorbed;  // introns normalise merged away while computing it
};

struct DepositScratch {
    std::vector<std::pair<std::int64_t, std::int64_t>> introns;   // normalised: sorted, disjoint, clipped
    std::vector<std::pair<std::int64_t, std::int64_t>> segments;  // the path, introns cut out
    std::vector<std::int32_t>                         sj_ids;     // annotated junction edges used
    /// ⭐ The same resolution kept PER INTRON POSITION, -1 where unannotated. `sj_ids` is filtered, so
    /// it cannot say which of a block's two ends is a junction — and the conserved mass needs exactly
    /// that. One entry per intron, so block `i` is bounded by `sj_id_at_gap[i-1]` and `sj_id_at_gap[i]`.
    std::vector<std::int32_t>                         sj_id_at_gap;
    std::vector<ScoredHypothesis>                     survivors;  // arbitration, per fragment
};

/// The QC denominators. Not optional and not derivable afterwards: every conservation statement
/// downstream must be able to name what it excluded.
struct DepositCounters {
    std::int64_t deposited              = 0;
    std::int64_t dropped_too_long       = 0;
    std::int64_t dropped_empty          = 0;
    std::int64_t dropped_strand_undefined = 0;
    //: ⭐ >1 surviving hypothesis, so the gap is undetermined. ⚠ NOT dropped -- the fragment is held
    //: WHOLE in `DeferredFragments` and the identity is deposited + deferred + dropped_* == offered.
    //: A name saying `dropped` for a population that is kept is how a recoverable loss gets read as
    //: a permanent one.
    std::int64_t deferred_undetermined_gap = 0;
    std::int64_t unannotated_introns    = 0;  // observed introns with no annotated junction
    std::int64_t contradictory_sj_strand = 0;  // the mates' motif tags disagreed; no splice trusted
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
    ///
    /// ⚠ `ref_id` is which reference this accumulator IS, and it is required rather than defaulted. An
    /// Accumulator is described by its cut positions alone and has no other way to know: it is stamped
    /// into every DEFERRED record, and the second pass replays those through `deposit` -- onto the wrong
    /// cut axis if the stamp is wrong, which is the failure mode `merge_from`'s cut comparison exists for.
    explicit Accumulator(std::vector<std::int64_t> cuts,
                         std::vector<std::uint8_t> node_types,
                         int max_length,
                         std::int32_t ref_id);

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

    /// ⭐ TRAPS: a-purity-filter-is-a-length-filter — EVERY deposited fragment, binned at its own L, with NO purity condition.
    ///
    /// The five pure pools above are deliberately CONDITIONED (an impure
    /// pool is worse than a missing one), so they cannot serve as the unconditional anchor an
    /// empirical-Bayes shrinkage needs -- which is why that anchor was taken from the SCANNER, which
    /// measures length by two other rules over another population.
    /// This row removes that reason: anchor and pools become one measurement of one quantity.
    ///
    /// It is "unconditional GIVEN DEPOSIT" and the name says so: it excludes what the accumulator
    /// rejects (too long / ambiguous path / strand-undefined / empty), each of which is counted in
    /// `DepositCounters`. That is exactly the population the pools are drawn from.
    ///
    /// Invariant, the same externally-checkable form as node_start_count's and a DIFFERENT statement:
    /// `sum(deposited_lengths) == sum(node_start_count) == deposited`.
    const std::uint32_t* deposited_lengths_data() const noexcept { return deposited_lengths_.data(); }
    std::size_t          deposited_lengths_size() const noexcept { return deposited_lengths_.size(); }
    int                 max_length()        const noexcept { return max_length_; }

    const DepositCounters& counters() const noexcept { return counters_; }
    const GapCensus&       gap_census() const noexcept { return gap_census_; }
    std::int32_t           ref_id()     const noexcept { return ref_id_; }

    /// The deferred queue in its ONE canonical order. ⛔ Non-const because it sorts: there is deliberately
    /// no accessor that hands out the append order, because a caller who exported that would produce a
    /// different byte sequence at every worker count and the determinism gate would fail on a difference
    /// that means nothing.
    const DeferredFragments& deferred_canonical();

    /// Index of the node containing `position`, clamped into [0, n_nodes - 1].
    ///
    /// ⚠ Clamped, not -1 on miss: `deposit` has already clipped the path into this reference, so a
    /// position outside cannot arrive, and clamping keeps this byte-compatible with the spec's
    /// `min(max(searchsorted(cuts, p, 'right') - 1, 0), n_cuts - 2)`.
    std::int64_t node_of_pos(std::int64_t position) const noexcept;

    /// Deposit one fragment. Allocates nothing: `scratch` is reused across calls.
    DepositOutcome deposit(const OfferedFragment& fragment, DepositScratch& scratch);

    /// ⭐ `L` under ONE hypothesis, without depositing anything — what the SECOND PASS scores against.
    ///
    /// ⛔ Exposed rather than reimplemented. The tool has ONE
    /// definition of fragment length, and the scorer needs a length per *counterfactual* hypothesis. A
    /// Python reimplementation would be a second definition of exactly the quantity that audit existed to
    /// unify — and it would be the one the drain then disagreed with.
    ///
    /// Returns the clipped `L`, or 0 when the fragment clips away entirely.
    std::int64_t length_under(const OfferedFragment& fragment,
                              std::size_t hypothesis_index,
                              DepositScratch& scratch) const;

    /// Element-wise sum of `other` into this accumulator. Requires identical cut positions.
    void merge_from(const Accumulator& other);

private:
    /// The one length pool this fragment belongs to, or -1 for none.
    ///
    /// ⭐ DETERMINACY, NOT PROVENANCE. There used to be an `sj_implicit` argument barring a fragment
    /// whose splice was inferred rather than sequenced. It is gone: a fragment reaches this line only
    /// when exactly ONE hypothesis survived, so its L is not in doubt however it was arrived at.
    /// Measured before deleting it -- the pool reads +0.67 % mean / +2.40 % sd against truth under
    /// determinacy and -9.58 % / -22.46 % under provenance, because barring inferred lengths
    /// preferentially bars fragments whose mates sit far apart. A purity filter on a length pool is a
    /// length filter.
    std::int64_t fragment_pool(bool spliced,
                               std::int64_t contained_node,
                               std::int64_t sole_line) const noexcept;

    /// L for one hypothesis: its implied introns UNIONED with the observed ones, normalised and clipped
    /// into [start, end). Leaves the normalised list in `scratch.introns`. ⭐ ONE definition of L and one
    /// code path to it, whether the hypothesis wins or is only being scored for the length filter.
    std::int64_t hypothesis_length(const OfferedFragment& fragment,
                                   const GapHypothesis& hypothesis,
                                   std::int64_t start,
                                   std::int64_t end,
                                   DepositScratch& scratch,
                                   std::int64_t* absorbed) const;

    /// ⭐ The umbrella census, at its ONE site. Called only when the fragment had a question to answer,
    /// i.e. at least one hypothesis was not the unspliced one.
    void record_gap_resolution(const OfferedFragment& fragment,
                               const std::vector<ScoredHypothesis>& survivors) noexcept;

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
    std::int32_t               ref_id_ = 0;        // stamped into every deferred record
    int                        max_length_ = 0;
    std::vector<std::int64_t>  pool_lengths_;      // kNFragmentPools * (max_length + 1), or empty
    std::vector<std::uint32_t> deposited_lengths_; // max_length + 1 -- unconditional given deposit
    DepositCounters            counters_;
    GapCensus                  gap_census_;
    DeferredFragments          deferred_;
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
