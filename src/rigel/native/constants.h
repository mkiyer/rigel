/**
 * constants.h — Shared constants and helper types for rigel C++ extensions.
 *
 * The enum-mirror constants must match the Python values each block names; the fragment classes are
 * defined here alone and Python reads them from the native module.
 */

#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <set>
#include <string>
#include <utility>
#include <vector>

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
// the two derived classes — neither is a CIGAR splice:
static constexpr int32_t SPLICE_IMPLICIT        = 2;  // PE gap spans an annotated intron
static constexpr int32_t SPLICE_ARTIFACT        = 3;  // CIGAR sj rejected by blacklist
// ⚠ Must equal `len(rigel.splice.SpliceType)`. It sizes the scanner's per-fragment splice census,
// so a category added on one side and not the other is caught by `test_every_splice_type_is_censused`
// rather than reading zero through the stats dict's `.get(key, 0)`.
static constexpr size_t NUM_SPLICE_TYPES        = 4;

// MergeOutcome (named for the gates in tests/_resolution_reference.py)
static constexpr int32_t MC_INTERSECTION          = 0;
static constexpr int32_t MC_INTERSECTION_NONEMPTY = 1;
static constexpr int32_t MC_UNION                 = 2;
static constexpr int32_t MC_EMPTY                 = 3;

// Strand (rigel.types.Strand)
static constexpr int32_t STRAND_NONE      = 0;
static constexpr int32_t STRAND_POS       = 1;
static constexpr int32_t STRAND_NEG       = 2;
static constexpr int32_t STRAND_AMBIGUOUS = 3;  // POS | NEG

// FragmentClass — the fragment classification codes, defined here alone: `_scoring_impl`, which routes on
// them, exports them, and rigel.buffer reads them from it (pinned by tests/test_buffer.py).
static constexpr int32_t FRAG_UNAMBIG           = 0;
static constexpr int32_t FRAG_AMBIG_SAME_STRAND = 1;
static constexpr int32_t FRAG_AMBIG_OPP_STRAND  = 2;
static constexpr int32_t FRAG_MULTIMAPPER       = 3;

// Assignment flags bitfield (written to BAM as ZF:i tag).
// ZF is the unified per-fragment outcome bitfield.  See
// src/rigel/annotate.py for the full schema + invariants; it composes
// every value but the three the C++ stamp sites write, so only those
// three values' bits are named here.
//
// Primitive bits
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

// ================================================================
// Scoring constants (shared by scoring.cpp and Python side)
// ================================================================
static constexpr double LOG_HALF       = -0.6931471805599453;       // log(0.5)
static constexpr double TAIL_DECAY_LP  = -0.01005033585350145;      // log(0.99)

// ================================================================
// Geometric helper types shared across modules
// ================================================================

struct ExonBlock {
    int32_t ref_id;
    int32_t start;
    int32_t end;
    int32_t strand;
};

struct IntronBlock {
    int32_t ref_id;
    int32_t start;
    int32_t end;
    int32_t strand;
};

// Merge result: sorted vector of transcript indices + criteria
struct MergeResult {
    std::vector<int32_t> t_inds;
    int32_t criteria;
    bool is_empty() const { return t_inds.empty(); }
};

// Fold one int32 key field into a running hash (the boost hash_combine step). 0x9e3779b9 is 2^32/φ, the
// golden-ratio increment that spreads consecutive values across the bits.
inline void hash_combine(size_t& h, int32_t v) {
    h ^= std::hash<int32_t>()(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
}

// SJ exact-match map key
struct SJKey {
    int32_t ref_id;
    int32_t start;
    int32_t end;
    int32_t strand;

    bool operator==(const SJKey& o) const {
        return ref_id == o.ref_id && start == o.start &&
               end == o.end && strand == o.strand;
    }
};

struct SJKeyHash {
    size_t operator()(const SJKey& k) const {
        size_t h = std::hash<int32_t>()(k.ref_id);
        hash_combine(h, k.start);
        hash_combine(h, k.end);
        hash_combine(h, k.strand);
        return h;
    }
};

// Splice-junction artifact blacklist key.
// Strand-agnostic: alignable records a single entry per (ref, start, end).
struct SJBlacklistKey {
    int32_t ref_id;
    int32_t start;
    int32_t end;

    bool operator==(const SJBlacklistKey& o) const {
        return ref_id == o.ref_id && start == o.start && end == o.end;
    }
};

struct SJBlacklistKeyHash {
    size_t operator()(const SJBlacklistKey& k) const {
        size_t h = std::hash<int32_t>()(k.ref_id);
        hash_combine(h, k.start);
        hash_combine(h, k.end);
        return h;
    }
};

struct SJBlacklistEntry {
    int32_t max_anchor_left;
    int32_t max_anchor_right;
};

// ================================================================
// RawResolveResult — internal result from _resolve_core()
// ================================================================

struct RawResolveResult {
    std::vector<int32_t> t_inds;
    int32_t ambig_strand = 0;
    int32_t splice_type = 0;
    int32_t align_strand = 0;
    int32_t sj_strand = 0;
    int32_t genomic_footprint = -1;
    int32_t genomic_start = -1;
    int32_t merge_criteria = MC_EMPTY;
    int32_t read_length = 0;
    // Parallel arrays to t_inds
    std::vector<int32_t> frag_lengths;
    std::vector<int32_t> t_exon_bp;
    std::vector<int32_t> t_intron_bp;

    // Per-fragment count of CIGAR splice junctions rejected by the
    // alignment-time blacklist.  Caller sets this BEFORE calling
    // _resolve_core so the resolver can promote SPLICE_UNSPLICED
    // to SPLICE_ARTIFACT.
    int32_t n_sj_blacklisted = 0;

    // The per-sj strand table's key ---
    // Coordinates of the LEFTMOST ANNOTATED CIGAR-N sj this fragment
    // crosses; -1 when it crosses none.  `sj_strand` above is the fragment's
    // motif strand (one XS/ts tag per fragment), which completes the key.
    //
    // ⚠ Crediting the FIRST annotated sj is a CHOICE, recorded here so it
    // is not mistaken for a derivation.  A qualified fragment carries ONE sense
    // bit: `sj_strand` is read from the BAM XS/ts tag, which is per RECORD, and
    // mates that disagree OR to STRAND_AMBIGUOUS, which the qualification
    // rejects.  So all K sj of a fragment necessarily agree.
    //
    // Two measured reasons for credit-one:
    //   (a) the 2×2 marginal identity REQUIRES one row per fragment;
    //   (b) the K >= 2 stratum is a materially cleaner population (minor rate
    //       0.35x at K=2, 0.14x at K=3, 0 at K>=4), so crediting all K shifts
    //       kappa by 21-34%.
    //
    // ⚠ A cost of the deterministic leftmost pick, recorded so it is not
    // rediscovered as a bug: 3-5% of sj never receive an observation (leftmost
    // and rightmost agree, so this is concentration on fewer sj, not a 5'/3'
    // bias). Determinism wins anyway: this repo has goldens.
    //
    // ⚠ ANNOTATED, not merely leftmost: `sj_strand` is the OR of the ANNOTATED
    // introns' strands only, so an unannotated intron may carry a different (or
    // absent) tag.  Keying on an annotated sj is what makes the 2×2
    // marginal identity hold unconditionally.
    int32_t sj_key_ref = -1;
    int32_t sj_key_start = -1;
    int32_t sj_key_end = -1;

    // ⭐ EVERY explanation of this fragment's UNSEQUENCED gaps, flat. Hypothesis h owns introns
    // [gap_intron_offsets[h], gap_intron_offsets[h+1]).
    //
    // A mate gap may hold no intron, one, or several, and which it is CANNOT be observed -- the bases are
    // not there. Each compatible transcript determines exactly one answer (its own introns lying inside
    // the gaps), so the set is finite and small, and two transcripts implying the same introns are ONE
    // hypothesis. The accumulator arbitrates.
    //
    // ⚠ An EMPTY hypothesis (no introns) is the UNSPLICED one, and it is the genomic explanation: cutting
    // nothing means the gap is real template. It is present whenever some compatible transcript implies
    // nothing there, or whenever the fragment carries no annotated CIGAR-N sj and so could be gDNA.
    // ⛔ It is never added as the synthetic span candidate -- the synthetic span IS this hypothesis.
    std::vector<IntronBlock> gap_introns;
    std::vector<int32_t>     gap_intron_offsets;      // n_hypotheses + 1, always starts at 0
    std::vector<int32_t>     gap_sj_strand;           // n_hypotheses

    int32_t n_gap_hypotheses() const {
        return gap_intron_offsets.empty()
                   ? 0
                   : static_cast<int32_t>(gap_intron_offsets.size()) - 1;
    }
};

// ================================================================
// Set operations
// ================================================================

inline void sort_unique(std::vector<int32_t>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

// Sorted-vector intersection.
inline std::vector<int32_t> vec_intersect(const std::vector<int32_t>& a,
                                          const std::vector<int32_t>& b) {
    std::vector<int32_t> out;
    out.reserve(std::min(a.size(), b.size()));
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                          std::back_inserter(out));
    return out;
}

// Progressive set merging.
inline MergeResult merge_sets(const std::vector<std::vector<int32_t>>& sets) {
    if (sets.empty()) return {{}, MC_EMPTY};

    // Separate non-empty sets
    std::vector<const std::vector<int32_t>*> non_empty;
    for (const auto& s : sets)
        if (!s.empty()) non_empty.push_back(&s);
    if (non_empty.empty()) return {{}, MC_EMPTY};

    bool has_empty_set = (non_empty.size() < sets.size());

    // 1. Intersection of ALL sets (only meaningful if none are empty)
    if (!has_empty_set) {
        std::vector<int32_t> isect = *non_empty[0];
        for (size_t i = 1; i < non_empty.size() && !isect.empty(); i++)
            isect = vec_intersect(isect, *non_empty[i]);
        if (!isect.empty())
            return {std::move(isect), MC_INTERSECTION};
    } else {
        // 2. Intersection of non-empty sets only
        if (non_empty.size() == 1) {
            return {*non_empty[0], MC_INTERSECTION_NONEMPTY};
        }
        std::vector<int32_t> isect = *non_empty[0];
        for (size_t i = 1; i < non_empty.size() && !isect.empty(); i++)
            isect = vec_intersect(isect, *non_empty[i]);
        if (!isect.empty())
            return {std::move(isect), MC_INTERSECTION_NONEMPTY};
    }

    // 3. Union of all sets
    std::vector<int32_t> sorted_union;
    std::vector<int32_t> tmp;
    for (const auto* s : non_empty) {
        if (sorted_union.empty()) {
            sorted_union = *s;
            continue;
        }
        tmp.clear();
        tmp.reserve(sorted_union.size() + s->size());
        std::set_union(sorted_union.begin(), sorted_union.end(),
                       s->begin(), s->end(),
                       std::back_inserter(tmp));
        sorted_union.swap(tmp);
    }
    if (!sorted_union.empty()) return {std::move(sorted_union), MC_UNION};
    return {{}, MC_EMPTY};
}

}  // namespace rigel
