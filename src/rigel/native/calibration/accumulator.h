/**
 * accumulator.h — the per-fragment tally built during the single-pass BAM scan.
 *
 *     SPEC:   tests/native/_accumulator_reference.py  — this file must reproduce it BYTE FOR BYTE.
 *             Where the two disagree, the Python file wins.
 *
 * THE MODEL
 *   The genome is a graph. One Accumulator holds ONE reference, described by its sorted REGION_BOUND positions.
 *   A reference contributing `c` region_bounds owns `c - 1` REGIONS and `c - 2` interior BOUNDARIES, and a boundary is a
 *   0-bp CONTIGUOUS BOUNDARY between two adjacent regions:
 *
 *       region_bounds    0        100       200       600        c = 4
 *       regions   [  n0  ][   n1   ][   n2   ]            c - 1 = 3
 *       boundaries            boundary 1    boundary 2               c - 2 = 2
 *
 *   A SJ BOUNDARY is a directed donor->acceptor link taken from the annotation. A fragment is a
 *   PATH: its aligned blocks joined across the mate gap, broken by introns.
 *
 *   Regions count fragments CONTAINED (the whole path fits inside one region); boundaries count fragments
 *   CROSSING. Each population stores only the channels something READS, and they differ: count,
 *   inv_length_sum (fixed point), length_sum, and -- on the contiguous boundaries -- the conserved mass.
 *
 * WHY MORE THAN ONE SUM
 *   With `placements` the number of admissible start positions -- L at a region, L-1 at a 0-bp boundary:
 *
 *       E[count]   = rho * E[placements]
 *       E[inv_length_sum] = rho * E[placements * (1/placements)] = rho   <- at an BOUNDARY, exactly
 *       E[length_sum]     = rho * E[placements * L]
 *
 * `inv_length_sum` is deliberately NOT called `density`: it is one at an boundary, exactly, and is NOT one
 * at a region, where the opportunity (region - L + 1)+ does not cancel. `length_sum` is the second tilt, and
 * without it the (count, inv_length_sum) pair has determinant mu_g - mu_r and so carries NO information
 * about the gDNA/RNA split whenever the two components share a mean length. See
 *
 *   The opportunity factor cancels identically at an boundary for ANY length distribution, which is why no
 *   divisor and no length model appear there. It does not cancel at a region.
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
//:      region_len 151    fixed 7.0e-10    double 5.8e-15      120,000x better
//:      region_len 400    fixed 1.7e-08    double 1.0e-13      170,000x better
//:      region_len 1000   fixed 2.0e-07    double 2.8e-13      714,000x better
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

/// A region: an interval. ONE population — the fragments contained inside it — in two strand columns.
///
/// ⚠ A `spanning` population (one segment covering the region whole) was removed on evidence: it reached
/// no evidence-starved region the region's own bounding BOUNDARIES did not already reach off capture, and 141
/// regions / 822 fragments (0.008 %) under it. Its mass is not lost — a spanning fragment crosses both of
/// the region's boundaries and is deposited there.
/// ⛔ Consequence, and it is structural: no spliced fragment touches the region axis at all, because a
/// spliced fragment can never be `contained` (both endpoints of an annotated intron are region_bounds).
struct Region {
    std::uint32_t contained_count[kNStrandColumns];
    /// ⭐ ONE value, while `contained_count` keeps two. The length moments are strand-AGNOSTIC -- which
    /// strand a read aligned to says nothing about whether the molecule was gDNA or RNA -- and every
    /// consumer summed the two columns before using them. The COUNTS keep both because the strand model
    /// is a Beta-Binomial over them, per strand.
    double contained_inv_opportunity_sum;
};
static_assert(sizeof(Region) == 16, "Region must be 16 bytes with no padding");

/// A contiguous boundary: the 0-bp boundary between two adjacent regions. `spliced` means the FRAGMENT used an
/// annotated sj somewhere -- not that this boundary is one. gDNA cannot be spliced, so a spliced
/// crossing is a certified RNA crossing.
struct Boundary {
    std::uint32_t unspliced_count[kNStrandColumns];
    std::uint32_t spliced_count[kNStrandColumns];
    /// ⭐ ONE value -- strand-agnostic, see `Region`.
    double unspliced_inv_length_sum;
    /// ⭐⭐ THE CONSERVED MASS, fixed point. A COUNT and a MASS are two different deposits and one
    /// number cannot be both: `unspliced_count` is `+1` on every boundary a fragment crosses, so a fragment
    /// books `max(K, 1)` of them; this sums to ONE per fragment, across all the boundaries it crosses.
    ///
    /// ⛔ ONE VALUE, NOT TWO, AND THE RULING STANDS **HERE** WHILE IT WAS REVERSED ON THE SJ AXIS
    /// (2026-08-13) — the premise that changed is specific to sj and does not reach this bank.
    /// `strand_deconv` reads the counts per column; nothing reads a BOUNDARY's mass per strand, because at a
    /// boundary the mass exists to turn an object-incidence total into a fragment count and that question
    /// has no strand in it. ⚠ `one-thing-varied`: widening this too would have been a second change with
    /// no named consumer. See `SpliceJunction::mass`.
    double unspliced_mass;
    /// ⭐ The same rule, routed by the same `spliced` flag — so `mass` is not the one channel that
    /// ignores the split. ⛔ A PARTIAL, never a conservation ledger: a spliced fragment's blocks with no
    /// interior boundary deposit nothing (their accounting is on the sj axis), so this sums to
    /// `crossed_block_len / L`. It is a per-BOUNDARY certified-RNA term, commensurate with the unspliced
    /// mass at the same boundary, and is NOT "the number of spliced fragments here".
    double spliced_mass;
};
static_assert(sizeof(Boundary) == 40, "Boundary must be 40 bytes with no padding");

/// A sj boundary: one exact donor->acceptor jump. Spliced by construction, so there is no unspliced
/// population; and it is not a genomic position, so it carries no structural flags.
struct SpliceJunction {
    std::uint32_t count[kNStrandColumns];
    /// ⭐ LIVE: `second_pass` scores a held fragment's sj evidence with it. `length_sum` was
    /// removed — nothing read it, and `pool_lengths`' RNA_SPLICED row already carries that
    /// population's length distribution.
    double inv_length_sum;
    /// ⭐⭐⭐ THE CONSERVED MASS'S THIRD AXIS. A spliced fragment's block that contains no interior
    /// boundary deposits on neither boundary bank, and is not `contained` either -- its path spans a sj,
    /// so it lies in no single region. Such a fragment existed on the incidence axis and on no conserved
    /// one, which is why a library fragment count was not computable. Measured on the origin-split
    /// oracle at ladder g50 capture_off: 1,222,375 of 4,830,713 RNA fragments (25.3 %) are in that
    /// population, against 0 of 4,997,761 gDNA fragments, because gDNA cannot splice.
    ///
    /// ⛔ The rule ADDS a boundary class; it does not re-apportion an existing one. A block that
    /// crossed a boundary is untouched, so `unspliced_mass` and `spliced_mass` are byte-identical to what
    /// they were. Spec: `_accumulator_reference.py`; gates: `tests/native/test_conserved_mass.py`.
    ///
    /// ⭐⭐⭐ **TWO VALUES, AND THIS REVERSES `Boundary::unspliced_mass`'s ONE-VALUE RULING ON
    /// THIS AXIS ONLY (owner, 2026-08-12). THE REVERSAL IS ADMISSIBLE BECAUSE THE PREMISE CHANGED, AND
    /// THE PREMISE IS RECORDED HERE SO IT IS NOT RE-LITIGATED IN EITHER DIRECTION.**
    /// The ruling was *"nothing reads a mass per strand"*. That is now false for sj and only for
    /// sj: an ARTIFACTUAL splice junction accumulates SYMMETRICALLY on both strands, exactly as
    /// gDNA does, so the strand model the tool already has can detect one — but only if it is given a
    /// per-strand observable, and the COUNT is not enough because a count cannot separate a sj
    /// used by many short fragments from one used by few long ones.
    /// ⚠ The second reason is structural: without this bank, artifact filtering needs TWO passes over
    /// the BAM (tally, filter, re-accumulate the mass), which is the one thing the single-pass
    /// architecture exists to avoid.
    ///
    /// ⛔ **The column is `col` — the SAME genome-strand column `count` is deposited at**, so
    /// `mass[c] / count[c]` is a per-strand mean and not a ratio of two different populations.
    /// ⚠ Summed over columns this is byte-comparable to the single accumulator it replaces, but NOT
    /// bit-identical: float addition is not associative and the deposit order per column differs.
    /// Agreement is ~1e-15 relative, which is the convention this file already documents.
    double mass[kNStrandColumns];
};
static_assert(sizeof(SpliceJunction) == 32, "SpliceJunction must be 32 bytes with no padding");

// ============================================================================
// the fragment-length pools
// ============================================================================

//: Five pools, each PURE BY CONSTRUCTION. Purity removes the circularity: a length model is fitted from
//: a population known to be one component, so nothing is estimated from the fragments it will explain.
//:
//: There is deliberately NO pool for an exonic contained fragment or a multi-boundary crossing -- those are
//: gDNA/RNA mixtures, and an impure pool is worse than a missing one.
enum class FragmentPool : std::uint8_t {
    kDnaIntergenic     = 0,  // contained in an intergenic region
    kDnaIntronic       = 1,  // contained in an intronic region
    kDnaIntronExon     = 2,  // crossing exactly one boundary, flanks {intron, exon} -- a "splash" read
    kDnaIntergenicExon = 3,  // crossing exactly one boundary, flanks {intergenic, exon}
    kRnaSpliced        = 4,  // using an annotated sj, splice OBSERVED
};
inline constexpr std::size_t kNFragmentPools = 5;

//: Coarse region types, as `signature.coarse_type_array` emits them.
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
/// ⚠ `introns` are the IMPLIED ones only. Introns the CIGAR actually stated are region_bound under EVERY
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
    const IntronBlock*   observed_introns;   // CIGAR-N: region_bound under EVERY hypothesis
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
/// ⛔ Its own axis, NOT a `splice_type`. The umbrella region_bounds ACROSS the splice census: a certified-RNA
/// SPLICED_ANNOT fragment with an intron in its mate gap needs resolving exactly as much as an UNSPLICED
/// one does, so putting these on `splice_type` would need two labels per fragment.
///
/// ⚠ These classify the ARBITRATION, not the deposit: a resolved_* fragment can still be rejected
/// afterwards as TOO_LONG, which is a different question with its own counter.
///
/// ⛔ THERE IS NO `resolved_unspliced`, AND IT IS NOT AN OMISSION. The field existed and no fragment could
/// enter it: a spliced hypothesis REGION_BOUNDS bases the unspliced one keeps, so L_spliced <= L_unspliced always,
/// and the one filter is `L <= max_length`. If the unspliced path survives the filter then every spliced
/// path survives it too, so the survivor set can never be exactly {unspliced} while a spliced path was
/// offered -- which is the condition for being in this census at all. The ORDERING is pinned directly by
/// `test_gap_hypothesis_arbitration.test_the_GENOMIC_hypothesis_is_ALWAYS_the_LONGEST`.
struct GapCensus {
    std::int64_t resolved_spliced        = 0;  // one survivor, and it necessarily region_bounds something
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
    std::vector<std::int64_t> ref;                         // which reference: the drain needs the region_bound axis
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
    std::vector<std::pair<std::int64_t, std::int64_t>> segments;  // the path, introns region_bound out
    std::vector<std::int32_t>                         sj_ids;     // annotated sj boundaries used
    /// ⭐ The same resolution kept PER INTRON POSITION, -1 where unannotated. `sj_ids` is filtered, so
    /// it cannot say which of a block's two ends is a sj — and the conserved mass needs exactly
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
    std::int64_t unannotated_introns    = 0;  // observed introns with no annotated sj
    std::int64_t contradictory_sj_strand = 0;  // the mates' motif tags disagreed; no splice trusted
    std::int64_t introns_absorbed       = 0;  // overlapping or abutting introns merged away

    void merge_from(const DepositCounters& other) noexcept;
};

// ============================================================================
// Accumulator — one reference
// ============================================================================

class Accumulator {
public:
    /// `region_bounds` is this reference's sorted, strictly increasing region_bound positions and is moved in. Length may
    /// be 0 or 1, which is a reference with no regions: legal, and it deposits nothing.
    ///
    /// `region_types` is either empty or one coarse type per region; it types the length pools. `max_length`
    /// is the fragment-length limit applied to L and the width of the pool histograms, and must be >= 1
    /// -- at 0 every real fragment would be dropped as too long and the whole tally would be silently
    /// empty.
    ///
    /// ⚠ `ref_id` is which reference this accumulator IS, and it is required rather than defaulted. An
    /// Accumulator is described by its region_bound positions alone and has no other way to know: it is stamped
    /// into every DEFERRED record, and the second pass replays those through `deposit` -- onto the wrong
    /// region_bound axis if the stamp is wrong, which is the failure mode `merge_from`'s region_bound comparison exists for.
    explicit Accumulator(std::vector<std::int64_t> region_bounds,
                         std::vector<std::uint8_t> region_types,
                         int max_length,
                         std::int32_t ref_id);

    /// Install this reference's sj boundaries as a CSR keyed by DONOR REGION_BOUND INDEX -- the index the
    /// deposit already computes while locating the boundaries its path crosses.
    ///
    /// ⚠ The sj-boundary id IS the slot: `sj_boundary_right[k]` and the bank entry `k` are the same k.
    /// There is no indirection to a row in `edges.feather`; using that row as a bank index writes past
    /// the end of a 404,168-entry array, because the highest such row is 1,447,755.
    ///
    /// ⚠ Slot ORDER is part of the contract, because the id is the rank: the caller must sort on
    /// (donor region_bound, acceptor region_bound, sj_strand), matching `Partition.from_region_bounds` in the Python spec.
    void set_sj(std::vector<std::int32_t> offsets,       // size n_region_bounds + 1
                       std::vector<std::int32_t> boundary_right,  // acceptor REGION_BOUND INDEX, not a coordinate
                       std::vector<std::int8_t>  sj_strand);    // the sj's ANNOTATED strand

    std::size_t n_regions()    const noexcept { return regions_.size(); }
    std::size_t n_boundaries()    const noexcept { return boundaries_.size(); }
    std::size_t n_sj() const noexcept { return sj_.size(); }
    std::size_t n_region_bounds()     const noexcept { return region_bounds_.size(); }

    const std::int64_t* region_bounds_data() const noexcept { return region_bounds_.data(); }
    Region*               regions_data()      noexcept { return regions_.data(); }
    const Region*         regions_data() const noexcept { return regions_.data(); }
    Boundary*     boundaries_data()      noexcept { return boundaries_.data(); }
    const Boundary* boundaries_data() const noexcept { return boundaries_.data(); }
    SpliceJunction*       sj_data()      noexcept { return sj_.data(); }
    const SpliceJunction* sj_data() const noexcept { return sj_.data(); }

    /// One uint32 per region counting fragments whose FIRST COVERED BASE lies in it.
    ///
    /// ⭐ This is the accumulator's one real invariant: `sum(region_start_count) == deposited`, checkable
    /// against a number the scanner knows independently. The three "conservation identities" it replaced
    /// were tautologies -- each right-hand side could only be evaluated by re-running the deposit, so a
    /// deliberately broken replay satisfied all three while 91 % of the crossings were junk.
    const std::uint32_t* region_start_count_data() const noexcept { return region_start_count_.data(); }

    /// Length histograms, pool-major: pool p occupies [p*(max_length+1), (p+1)*(max_length+1)), binned
    /// at L. Empty when this reference has no region types.
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
    /// Invariant, the same externally-checkable form as region_start_count's and a DIFFERENT statement:
    /// `sum(deposited_lengths) == sum(region_start_count) == deposited`.
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

    /// Index of the region containing `position`, clamped into [0, n_regions - 1].
    ///
    /// ⚠ Clamped, not -1 on miss: `deposit` has already clipped the path into this reference, so a
    /// position outside cannot arrive, and clamping keeps this byte-compatible with the spec's
    /// `min(max(searchsorted(region_bounds, p, 'right') - 1, 0), n_region_bounds - 2)`.
    std::int64_t region_of_pos(std::int64_t position) const noexcept;

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

    /// Element-wise sum of `other` into this accumulator. Requires identical region_bound positions.
    void merge_from(const Accumulator& other);

private:
    /// The one length pool this fragment belongs to, or -1 for none.
    ///
    /// ⭐ DETERMINACY, NOT PROVENANCE. There used to be an `sj_implicit` argument barring a fragment
    /// whose splice was inferred rather than sequenced. It is gone: a fragment reaches this boundary only
    /// when exactly ONE hypothesis survived, so its L is not in doubt however it was arrived at.
    /// Measured before deleting it -- the pool reads +0.67 % mean / +2.40 % sd against truth under
    /// determinacy and -9.58 % / -22.46 % under provenance, because barring inferred lengths
    /// preferentially bars fragments whose mates sit far apart. A purity filter on a length pool is a
    /// length filter.
    std::int64_t fragment_pool(bool spliced,
                               std::int64_t contained_region,
                               std::int64_t sole_boundary) const noexcept;

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
    std::vector<std::uint32_t> region_start_count_;  // n_regions -- its own array, so Region stays 48 B

    std::vector<std::int32_t>  sj_offsets_;        // n_region_bounds + 1, CSR over the donor region_bound index
    std::vector<std::int32_t>  sj_boundary_right_;   // n_sj
    std::vector<std::int8_t>   sj_strand_;         // n_sj, the ANNOTATED strand

    std::vector<std::uint8_t>  region_types_;        // n_regions, or empty (no pools)
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
// `region_bounds` is the concatenated, reference-major region_bound array; reference f owns
// region_bounds[ref_region_bound_offsets[f] .. ref_region_bound_offsets[f+1]). A reference with fewer than 2 region_bounds owns no
// regions and no boundaries, which is legal.
//
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

    /// Install the sj CSR for every reference at once, from the FLAT arrays the index emits
    /// (`build_sj_arrays`), keyed by the flat region_bound index.
    ///
    /// The slicing is pure arithmetic and lives here so it exists once. Reference `f` owns region_bounds
    /// `[c0, c1)`, and because the CSR is sorted by flat donor region_bound while references are region_bound-major, its
    /// sj are the contiguous SLOT range `[offsets[c0], offsets[c1])`. So per reference:
    ///
    ///     offsets      -> offsets[c0 .. c1]      - offsets[c0]      (length n_region_bounds + 1)
    ///     boundary_right -> boundary_right[j0 .. j1] - c0               (a ref-local region_bound index)
    ///
    /// ⚠ Two consequences, both load-bearing. A reference's sj-boundary ids are `slot - j0`, so the
    /// payload's sj axis is exactly the flat slot order concatenated in reference order — which is
    /// what lets `edges_df.edge_row` stay a join key that never crosses the ABI. And the narrowing to
    /// int32 is safe by census: 1.04 M region_bounds and 404,168 sj at human scale.
    void set_sj(const std::int64_t* offsets,       // n_region_bounds_total + 1, over the FLAT region_bound axis
                       std::size_t n_offsets,
                       const std::int64_t* boundary_right,  // n_sj_total, FLAT region_bound indices
                       const std::int8_t* sj_strand,
                       std::size_t n_sj,
                       const std::int64_t* ref_region_bound_offsets);

    Accumulator&       at(std::int32_t ref_id);
    const Accumulator& at(std::int32_t ref_id) const;

    void merge_from(const AccumulatorSet& other);

private:
    std::vector<Accumulator> accs_;
};

}  // namespace rigel::accumulator
