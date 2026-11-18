/**
 * constants.h — Shared constants for rigel C++ extensions.
 *
 * The enum-mirror constants must match the Python values each block names; the fragment classes are
 * defined here alone and Python reads them from the native module.
 */

#pragma once

#include <cstddef>
#include <cstdint>

namespace rigel {

// ================================================================
// Enum-mirror constants
// ================================================================

// IntervalType (rigel.types.IntervalType)
static constexpr int8_t ITYPE_EXON           = 0;
static constexpr int8_t ITYPE_TRANSCRIPT      = 1;

// SpliceType (rigel.splice.SpliceType). One spliced type: every junction a fragment observes is held by one
// annotated transcript (the scan's annotated-path rule), so no fragment is spliced across an unannotated junction.
static constexpr int32_t SPLICE_UNSPLICED       = 0;
static constexpr int32_t SPLICE_SPLICED_ANNOT   = 1;
// The two derived classes, neither a CIGAR splice:
static constexpr int32_t SPLICE_IMPLICIT        = 2;  // PE gap spans an annotated intron
static constexpr int32_t SPLICE_ARTIFACT        = 3;  // CIGAR sj rejected by blacklist
// Must equal `len(rigel.splice.SpliceType)`: it sizes the scanner's per-fragment splice census.
static constexpr size_t NUM_SPLICE_TYPES        = 4;

// Strand (rigel.types.Strand)
static constexpr int32_t STRAND_NONE      = 0;
static constexpr int32_t STRAND_POS       = 1;
static constexpr int32_t STRAND_NEG       = 2;
static constexpr int32_t STRAND_AMBIGUOUS = 3;  // POS | NEG

// FragmentClass: the fragment classification codes, defined here alone. `_scoring_impl` routes on them and exports
// them, and rigel.buffer reads them from it.
static constexpr int32_t FRAG_UNAMBIG           = 0;
static constexpr int32_t FRAG_AMBIG_SAME_STRAND = 1;
static constexpr int32_t FRAG_AMBIG_OPP_STRAND  = 2;
static constexpr int32_t FRAG_MULTIMAPPER       = 3;

// The per-fragment outcome bitfield, written to the BAM as the ZF:i tag. src/rigel/annotate.py holds the full schema
// and composes every value but the three the C++ stamp sites write, so only those values' bits are named here.
static constexpr int32_t AF_RESOLVED_BIT         = 0x01;
static constexpr int32_t AF_GDNA_BIT             = 0x04;
static constexpr int32_t AF_INTERGENIC_BIT       = 0x20;
static constexpr int32_t AF_MULTIMAPPER_DROP_BIT = 0x80;

// The composed values the C++ stamp sites write.
static constexpr int32_t AF_UNRESOLVED        = 0x00;  // withheld: the template has no usable hit
static constexpr int32_t AF_GDNA_INTERGENIC   = AF_RESOLVED_BIT | AF_GDNA_BIT | AF_INTERGENIC_BIT;  // 0x25
static constexpr int32_t AF_MULTIMAPPER_DROP  = AF_MULTIMAPPER_DROP_BIT;                            // 0x80

// Why a hit is not usable (`hit_verdict`); HIT_USABLE when it is. Its two mates are not one molecule
// (`FragmentResolver::validate_pair`):
static constexpr int32_t HIT_USABLE       = 0;
static constexpr int32_t HIT_REFERENCE    = 1;  // the mates lie on two references
static constexpr int32_t HIT_STRAND       = 2;  // on one strand
static constexpr int32_t HIT_OUTWARD      = 3;  // their 5' ends facing away from each other
static constexpr int32_t HIT_LENGTH       = 4;  // no path gives a fragment length in [1, max_fragment_length]
static constexpr int32_t HIT_MISSING_MATE = 5;  // one mate only: no fragment length at all
// Or they are, and no one annotated transcript explains the hit: holds every junction it observes and overlaps both
// mates (the annotated-path rule, `FragmentResolver::pair_holders`):
static constexpr int32_t HIT_CONTRADICTORY  = 6;  // two junctions overlap unequal, or their motifs name both strands
static constexpr int32_t HIT_NOVEL_JUNCTION = 7;  // a junction no transcript holds
static constexpr int32_t HIT_UNSUPPORTED    = 8;  // each junction held, but no one holder of them all over both mates
// Or it lies where Rigel can place nothing:
static constexpr int32_t HIT_OFF_INDEX      = 9;   // both mates on one reference the index does not hold
static constexpr int32_t HIT_CHIMERIC       = 10;  // a mate aligned in pieces: its record carries SA, or a 5' hard clip
static constexpr int32_t N_HIT_VERDICTS     = 11;

/// A molecule the annotated-path rule withheld: its template, left with no usable hit, is an unannotated splice.
inline bool is_unannotated_splice(int32_t verdict) {
    return verdict == HIT_CONTRADICTORY || verdict == HIT_NOVEL_JUNCTION || verdict == HIT_UNSUPPORTED;
}

// Why a template is withheld: it has candidate hits and none is usable (`judge_template`). The label is the annotated
// BAM's ZC and the stats export's name.
static constexpr int32_t WITHHELD_NONE               = 0;
static constexpr int32_t WITHHELD_UNPAIRED           = 1;  // none of its hits is a molecule
static constexpr int32_t WITHHELD_UNANNOTATED_SPLICE = 2;  // the annotated-path rule withheld one of its hits
static constexpr int32_t WITHHELD_OFF_INDEX          = 3;  // one of its hits lies on a reference outside the index
static constexpr int32_t WITHHELD_CHIMERIC           = 4;  // one of its hits has a mate aligned in pieces
static constexpr int32_t N_WITHHELD                  = 5;

// The rest of a read pair's outcomes, numbered after the WITHHELD_* ones (`pair_digest`)
static constexpr int32_t PAIR_EXCLUDED_BY_OPTION  = N_WITHHELD;
static constexpr int32_t PAIR_PLACED_UNIQUE       = N_WITHHELD + 1;
static constexpr int32_t PAIR_PLACED_MULTIMAPPING = N_WITHHELD + 2;

// ================================================================
// Scoring constants (shared by scoring.cpp and the Python side)
// ================================================================
static constexpr double LOG_HALF       = -0.6931471805599453;       // log(0.5): a strand with no preference
static constexpr double TAIL_DECAY_LP  = -0.01005033585350145;      // log(0.99): per-bp penalty past the length table
// Candidates scored in a stack buffer before falling back to the heap, so the common fragment needs no malloc.
static constexpr int SCORED_STACK_CAPACITY = 64;

// ================================================================
// The locus EM (_em_impl)
// ================================================================

// The EM's log floor, log(θ + EM_LOG_EPSILON), and the smallest total it divides by.
static constexpr double EM_LOG_EPSILON = 1e-300;
// One SQUAREM iteration runs three EM steps (two plain, one stabilising), so `iterations / 3` SQUAREM iterations
// spend the budget `EMConfig.iterations` names.
static constexpr int SQUAREM_BUDGET_DIVISOR = 3;
// VBEM's floor: the smallest alpha a component keeps. digamma(1e-300) ≈ −1e300, so after max-subtraction the E-step
// gives a component at the floor exactly zero responsibility and the evidence-proportional prior gives it nothing: it
// stays dead unless it holds deterministic fragments. Only the EM's own step may put a component there, never
// SQUAREM's extrapolation (`backtracked_squarem_step`).
static constexpr double VBEM_SQUAREM_PRIOR_FLOOR = EM_LOG_EPSILON;

// The assignment codes are labels, not quantities. `_em_impl` exports them for rigel.estimator, and
// `batch_locus_em_partitioned` refuses any other value.
static constexpr int ASSIGN_FRACTIONAL = 0;  // every candidate keeps its posterior share
static constexpr int ASSIGN_SAMPLE     = 1;  // count first, then one component per unit

// A unit's winner beside a transcript index: none, the locus's own gDNA, or the intergenic term.
static constexpr int32_t EM_WINNER_NONE       = -1;
static constexpr int32_t EM_WINNER_GDNA       = -2;
static constexpr int32_t EM_WINNER_INTERGENIC = -3;

// The columns of a unit's class masses, which the EM returns when asked (`batch_locus_em_partitioned`'s
// `t_is_synthetic`): its posterior on the annotated transcripts, on the synthetic (nascent) spans, on the locus's gDNA
// component and on the intergenic term. `_em_impl` exports their names in this order.
static constexpr int EM_CLASS_ANNOTATED  = 0;
static constexpr int EM_CLASS_NASCENT    = 1;
static constexpr int EM_CLASS_GDNA       = 2;
static constexpr int EM_CLASS_INTERGENIC = 3;
static constexpr int N_EM_CLASSES        = 4;

// ================================================================
// The calibration solve (_solve_impl)
// ================================================================

// The message policies `solve_blocks` switches on; the silent policy runs no message layer.
static constexpr int POLICY_SILENT   = 0;
static constexpr int POLICY_TRANSFER = 1;
// A slot holds its own composition evidence when its strand evidence exceeds this guard; the strand evidence clips
// its fractions to it.
static constexpr double OWN_EVIDENCE_EPS = 1.0e-9;

// ================================================================
// The scan's capacities: bounds on memory and hand-off, not model parameters
// ================================================================

// Input-queue capacity per worker: the read-name batches the reader may hold queued before it blocks.
static constexpr int INPUT_BATCHES_PER_WORKER = 4;
// Output-queue capacity per worker: the finished chunks the workers may hold queued before they block on the main
// thread.
static constexpr int OUTPUT_CHUNKS_PER_WORKER = 2;
// The transcripts one fragment touches, the resolver's initial capacity for their list: a hint, it grows past it.
static constexpr size_t TOUCHED_TRANSCRIPTS_RESERVE = 512;

}  // namespace rigel
