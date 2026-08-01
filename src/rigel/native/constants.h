/**
 * constants.h — Shared constants and helper types for rigel C++ extensions.
 *
 * Must match the Python IntEnum values exactly (rigel.types,
 * rigel.categories).  Included by resolve.cpp, bam_scanner.cpp, and
 * any future native modules.
 */

#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
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

// SpliceType (rigel.splice.SpliceType)
static constexpr int32_t SPLICE_UNSPLICED       = 0;
static constexpr int32_t SPLICE_SPLICED_UNANNOT = 1;
static constexpr int32_t SPLICE_SPLICED_ANNOT   = 2;
// SRD v2 additions:
static constexpr int32_t SPLICE_IMPLICIT        = 3;  // PE gap spans an annotated intron
static constexpr int32_t SPLICE_ARTIFACT        = 4;  // CIGAR junction rejected by blacklist
// ⚠ Must equal `len(rigel.splice.SpliceType)`. It sizes the scanner's per-fragment splice census,
// so a category added on one side and not the other is caught by `test_every_splice_type_is_censused`
// rather than reading zero through the stats dict's `.get(key, 0)`.
static constexpr size_t NUM_SPLICE_TYPES        = 5;

// MergeOutcome (rigel.types.MergeOutcome)
static constexpr int32_t MC_INTERSECTION          = 0;
static constexpr int32_t MC_INTERSECTION_NONEMPTY = 1;
static constexpr int32_t MC_UNION                 = 2;
static constexpr int32_t MC_EMPTY                 = 3;

// ChimeraType (rigel.types.ChimeraType)
static constexpr int32_t CHIMERA_NONE           = 0;
static constexpr int32_t CHIMERA_TRANS           = 1;
static constexpr int32_t CHIMERA_CIS_STRAND_SAME = 2;
static constexpr int32_t CHIMERA_CIS_STRAND_DIFF = 3;

// Strand (rigel.types.Strand)
static constexpr int32_t STRAND_NONE      = 0;
static constexpr int32_t STRAND_POS       = 1;
static constexpr int32_t STRAND_NEG       = 2;
static constexpr int32_t STRAND_AMBIGUOUS = 3;  // POS | NEG

// FragmentClass (rigel.buffer — fragment classification codes)
static constexpr int32_t FRAG_UNAMBIG           = 0;
static constexpr int32_t FRAG_AMBIG_SAME_STRAND = 1;
static constexpr int32_t FRAG_AMBIG_OPP_STRAND  = 2;
static constexpr int32_t FRAG_MULTIMAPPER       = 3;
static constexpr int32_t FRAG_CHIMERIC          = 4;

// Assignment flags bitfield (written to BAM as ZF:i tag).
// ZF is the unified per-fragment outcome bitfield.  See
// src/rigel/annotate.py for the full schema + invariants.  C++ stamp
// sites reference only the composed values below.
//
// Primitive bits
static constexpr int32_t AF_RESOLVED_BIT         = 0x01;
static constexpr int32_t AF_MRNA_BIT             = 0x02;
static constexpr int32_t AF_GDNA_BIT             = 0x04;
static constexpr int32_t AF_NRNA_BIT             = 0x08;
static constexpr int32_t AF_SYNTHETIC_BIT        = 0x10;
static constexpr int32_t AF_INTERGENIC_BIT       = 0x20;
static constexpr int32_t AF_CHIMERIC_BIT         = 0x40;
static constexpr int32_t AF_MULTIMAPPER_DROP_BIT = 0x80;

// Canonical composed values (the only legitimate ZF outputs).
static constexpr int32_t AF_UNRESOLVED        = 0x00;
static constexpr int32_t AF_MRNA              = AF_RESOLVED_BIT | AF_MRNA_BIT;               // 0x03
static constexpr int32_t AF_NRNA              = AF_RESOLVED_BIT | AF_NRNA_BIT;               // 0x09
static constexpr int32_t AF_NRNA_SYNTH        = AF_NRNA | AF_SYNTHETIC_BIT;                  // 0x19
static constexpr int32_t AF_GDNA_EM           = AF_RESOLVED_BIT | AF_GDNA_BIT;               // 0x05
static constexpr int32_t AF_GDNA_INTERGENIC   = AF_GDNA_EM | AF_INTERGENIC_BIT;              // 0x25
static constexpr int32_t AF_CHIMERIC          = AF_CHIMERIC_BIT;                             // 0x40
static constexpr int32_t AF_MULTIMAPPER_DROP  = AF_MULTIMAPPER_DROP_BIT;                     // 0x80

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

// Chimera detection result
struct ChimeraResult {
    int32_t type;
    int32_t gap;
};

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
        h ^= std::hash<int32_t>()(k.start)  + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int32_t>()(k.end)    + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int32_t>()(k.strand) + 0x9e3779b9 + (h << 6) + (h >> 2);
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
        h ^= std::hash<int32_t>()(k.start) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int32_t>()(k.end)   + 0x9e3779b9 + (h << 6) + (h >> 2);
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
    int32_t chimera_type = CHIMERA_NONE;
    int32_t chimera_gap = -1;
    // Parallel arrays to t_inds
    std::vector<int32_t> frag_lengths;
    std::vector<int32_t> t_exon_bp;
    std::vector<int32_t> t_intron_bp;

    // --- SRD v2: strand-aware collapsed overlap counts ---
    // bp of fragment overlapping ANY (+/-)-strand transcript's exon / span.
    int32_t exon_bp_pos = 0;
    int32_t exon_bp_neg = 0;
    int32_t tx_bp_pos = 0;
    int32_t tx_bp_neg = 0;
    // Per-fragment count of CIGAR splice junctions rejected by the
    // alignment-time blacklist.  Caller sets this BEFORE calling
    // _resolve_core so the resolver can promote SPLICE_UNSPLICED
    // to SPLICE_ARTIFACT.
    int32_t n_sj_blacklisted = 0;

    // --- The per-junction SJ strand table's key (docs/CARRY_FORWARD.md) ---
    // Coordinates of the LEFTMOST ANNOTATED CIGAR-N junction this fragment
    // crosses; -1 when it crosses none.  `sj_strand` above is the fragment's
    // motif strand (one XS/ts tag per fragment), which completes the key.
    //
    // ⚠ Crediting the FIRST annotated junction is a CHOICE, recorded here so it
    // is not mistaken for a derivation.  A qualified fragment carries ONE sense
    // bit: `sj_strand` is read from the BAM XS/ts tag, which is per RECORD, and
    // mates that disagree OR to STRAND_AMBIGUOUS, which the qualification
    // rejects.  So all K junctions of a fragment necessarily agree.
    //
    // Two measured reasons for credit-one (2026-07-28):
    //   (a) the 2×2 marginal identity REQUIRES one row per fragment;
    //   (b) the K >= 2 stratum is a materially cleaner population (minor rate
    //       0.35x at K=2, 0.14x at K=3, 0 at K>=4), so crediting all K shifts
    //       kappa by 21-34% — and the fit feeds kappa in as the MoM node mean.
    // ⛔ NOT because credit-all "inflates" the dispersion: measured null bias is
    // -7.3e-5 +/- 3.0e-4 on LBX0190, i.e. none. That reasoning was refuted.
    // ⛔ A 1/K split is provably BIASED (4-12 sigma): Var(sum w_i X_i) = pq*sum w^2
    // but the estimator subtracts pq*sum w. Do not use it.
    //
    // ⚠ Two costs of the deterministic leftmost pick, recorded so they are not
    // rediscovered as bugs: od_mom is 0.85-0.99x the random-pick value, and 3-5%
    // of junctions never receive an observation (leftmost and rightmost agree, so
    // this is concentration on fewer junctions, not a 5'/3' bias). Determinism
    // wins anyway — random-pick jitter is 3.8% CV and this repo has goldens.
    //
    // ⚠ ANNOTATED, not merely leftmost: `sj_strand` is the OR of the ANNOTATED
    // introns' strands only, so an unannotated intron may carry a different (or
    // absent) tag.  Keying on an annotated junction is what makes the 2×2
    // marginal identity of §2.1 hold unconditionally.
    int32_t sj_key_ref = -1;
    int32_t sj_key_start = -1;
    int32_t sj_key_end = -1;

    // Implicit-splice introns: annotated introns found wholly inside a
    // paired-end mate gap (SPLICE_IMPLICIT). One per matched gap, carrying the
    // matched transcript's strand. The accumulator cuts these out of the
    // fragment span and orients the spliced channel by their strand (the splice
    // motif itself was not sequenced). Empty unless splice_type == SPLICE_IMPLICIT.
    std::vector<IntronBlock> implicit_introns;

    // ⛔ The candidate transcripts do NOT agree on which introns the mate gaps contain, so the implied
    // intron set -- and therefore L, the density quanta, the length pool and the set of lines the
    // fragment crosses -- is undetermined. `implicit_introns` above is then only ONE candidate's answer
    // and must not be tallied: it is the accumulator's `path_ambiguous`, deferred to the second pass,
    // which can separate the candidates by fragment length and strand. Design §9.1.
    //
    // ⚠ Meaningless unless splice_type == SPLICE_IMPLICIT.
    bool implicit_ambiguous = false;
};

// ================================================================
// Set operations
// ================================================================

inline void sort_unique(std::vector<int32_t>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

// Check whether two sorted vectors share at least one element.
inline bool has_intersection(const std::vector<int32_t>& a,
                             const std::vector<int32_t>& b) {
    auto ai = a.begin(), bi = b.begin();
    while (ai != a.end() && bi != b.end()) {
        if (*ai < *bi) ++ai;
        else if (*ai > *bi) ++bi;
        else return true;
    }
    return false;
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

// Progressive set merging identical to Python merge_sets_with_criteria().
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

// Intrachromosomal chimera detection
inline ChimeraResult detect_chimera(
    const std::vector<ExonBlock>& exons,
    const std::vector<std::vector<int32_t>>& exon_t_sets)
{
    // Filter to blocks with non-empty transcript sets
    std::vector<int> item_idx;
    for (size_t i = 0; i < exon_t_sets.size(); i++)
        if (!exon_t_sets[i].empty()) item_idx.push_back(static_cast<int>(i));

    int n = static_cast<int>(item_idx.size());
    if (n <= 1) return {CHIMERA_NONE, -1};

    // Union-find
    std::vector<int> parent(n);
    std::iota(parent.begin(), parent.end(), 0);

    auto find = [&](int x) -> int {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra != rb) parent[ra] = rb;
    };

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (has_intersection(exon_t_sets[item_idx[i]],
                                 exon_t_sets[item_idx[j]]))
                unite(i, j);

    // Group into connected components with vector-backed root lookup.
    std::vector<int> roots(n);
    for (int i = 0; i < n; i++) roots[i] = find(i);
    std::vector<int> unique_roots = roots;
    std::sort(unique_roots.begin(), unique_roots.end());
    unique_roots.erase(std::unique(unique_roots.begin(), unique_roots.end()),
                       unique_roots.end());
    if (unique_roots.size() <= 1) return {CHIMERA_NONE, -1};

    std::vector<std::vector<int>> components(unique_roots.size());
    for (int i = 0; i < n; i++) {
        auto it = std::lower_bound(unique_roots.begin(), unique_roots.end(), roots[i]);
        components[static_cast<size_t>(it - unique_roots.begin())].push_back(i);
    }

    // Strand characterisation
    std::vector<int32_t> unique_strands;
    unique_strands.reserve(components.size());
    for (const auto& members : components) {
        int32_t strand = STRAND_NONE;
        for (int idx : members) strand |= exons[item_idx[idx]].strand;
        unique_strands.push_back(strand);
    }
    sort_unique(unique_strands);
    int32_t chimera_type = (unique_strands.size() == 1)
        ? CHIMERA_CIS_STRAND_SAME : CHIMERA_CIS_STRAND_DIFF;

    // Minimum gap between components
    int32_t min_gap = std::numeric_limits<int32_t>::max();
    for (size_t ci = 0; ci < components.size(); ci++) {
        for (size_t cj = ci + 1; cj < components.size(); cj++) {
            for (int bi : components[ci]) {
                for (int bj : components[cj]) {
                    const auto& blki = exons[item_idx[bi]];
                    const auto& blkj = exons[item_idx[bj]];
                    int32_t gap;
                    if (blki.end <= blkj.start) gap = blkj.start - blki.end;
                    else if (blkj.end <= blki.start) gap = blki.start - blkj.end;
                    else gap = 0;
                    min_gap = std::min(min_gap, gap);
                }
            }
        }
    }
    return {chimera_type, min_gap};
}

}  // namespace rigel
