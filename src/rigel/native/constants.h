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
#include <iterator>
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
// The two derived classes, neither a CIGAR splice:
static constexpr int32_t SPLICE_IMPLICIT        = 2;  // PE gap spans an annotated intron
static constexpr int32_t SPLICE_ARTIFACT        = 3;  // CIGAR sj rejected by blacklist
// Must equal `len(rigel.splice.SpliceType)`: it sizes the scanner's per-fragment splice census.
static constexpr size_t NUM_SPLICE_TYPES        = 4;

// MergeOutcome: the resolver's merge levels (mirrored in tests/_resolution_reference.py)
static constexpr int32_t MC_INTERSECTION          = 0;
static constexpr int32_t MC_INTERSECTION_NONEMPTY = 1;
static constexpr int32_t MC_UNION                 = 2;
static constexpr int32_t MC_EMPTY                 = 3;

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

// ================================================================
// Scoring constants (shared by scoring.cpp and the Python side)
// ================================================================
static constexpr double LOG_HALF       = -0.6931471805599453;       // log(0.5): a strand with no preference
static constexpr double TAIL_DECAY_LP  = -0.01005033585350145;      // log(0.99): per-bp penalty past the length table

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

    // CIGAR junctions rejected by the alignment-time blacklist. The caller sets it before _resolve_core, so the
    // resolver can promote SPLICE_UNSPLICED to SPLICE_ARTIFACT.
    int32_t n_sj_blacklisted = 0;

    // The per-junction strand table's key: the leftmost annotated CIGAR-N junction this fragment crosses (-1 when
    // none), completed by `sj_strand`, the fragment's one motif strand. A qualified fragment carries one strand bit
    // (mates whose XS/ts tags disagree OR to STRAND_AMBIGUOUS and are rejected), so all its junctions agree, and
    // crediting one junction per fragment keeps the 2×2 marginal identity. The pick is the leftmost so it is
    // deterministic, and it must be an annotated junction because `sj_strand` is the OR of annotated introns only.
    int32_t sj_key_ref = -1;
    int32_t sj_key_start = -1;
    int32_t sj_key_end = -1;

    // Every explanation of this fragment's unsequenced gaps, flat: hypothesis h owns introns
    // [gap_intron_offsets[h], gap_intron_offsets[h+1]) (`FragmentResolver::enumerate_gap_hypotheses`). The empty
    // hypothesis is the unspliced, genomic one; the synthetic span is this hypothesis, never a separate one.
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
