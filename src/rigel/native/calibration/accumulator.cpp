/**
 * accumulator.cpp — the deposit rule.
 *
 *     SPEC: tests/native/_accumulator_reference.py — byte for byte. Where the two disagree, it wins.
 *
 * ⚠ Every ordering decision in `deposit` is load-bearing and each one below was a real bug first. Read
 * the comments before reordering anything.
 */
#include "accumulator.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace rigel::accumulator {

void DepositCounters::merge_from(const DepositCounters& other) noexcept {
    deposited                += other.deposited;
    dropped_too_long         += other.dropped_too_long;
    dropped_empty            += other.dropped_empty;
    dropped_strand_undefined += other.dropped_strand_undefined;
    deferred_undetermined_gap += other.deferred_undetermined_gap;
    unannotated_introns      += other.unannotated_introns;
    contradictory_sj_strand  += other.contradictory_sj_strand;
    introns_absorbed         += other.introns_absorbed;
}

void GapCensus::merge_from(const GapCensus& other) noexcept {
    resolved_spliced       += other.resolved_spliced;
    deferred_rna_or_gdna   += other.deferred_rna_or_gdna;
    deferred_which_introns += other.deferred_which_introns;
    deferred_both          += other.deferred_both;
}

void DeferredFragments::append(const OfferedFragment& fragment,
                               std::int64_t ref_id,
                               std::int64_t clipped_start,
                               std::int64_t clipped_end) {
    ref.push_back(ref_id);
    start.push_back(clipped_start);
    end.push_back(clipped_end);
    align_strand.push_back(fragment.align_strand);
    sj_strand.push_back(fragment.sj_strand);
    for (std::size_t i = 0; i < fragment.n_observed_introns; ++i) {
        observed_introns.push_back(fragment.observed_introns[i].start);
        observed_introns.push_back(fragment.observed_introns[i].end);
    }
    observed_intron_offsets.push_back(static_cast<std::int64_t>(observed_introns.size() / 2));
    for (std::size_t h = 0; h < fragment.n_hypotheses; ++h) {
        const GapHypothesis& hypothesis = fragment.hypotheses[h];
        hypothesis_sj_strand.push_back(hypothesis.sj_strand);
        for (std::size_t i = 0; i < hypothesis.n_introns; ++i) {
            hypothesis_introns.push_back(hypothesis.introns[i].start);
            hypothesis_introns.push_back(hypothesis.introns[i].end);
        }
        hypothesis_intron_offsets.push_back(static_cast<std::int64_t>(hypothesis_introns.size() / 2));
        for (std::size_t i = 0; i < hypothesis.n_supporting; ++i)
            hypothesis_t.push_back(hypothesis.supporting_t[i]);
        hypothesis_t_offsets.push_back(static_cast<std::int64_t>(hypothesis_t.size()));
    }
    hypothesis_offsets.push_back(static_cast<std::int64_t>(hypothesis_sj_strand.size()));
}

void DeferredFragments::merge_from(const DeferredFragments& other) {
    // ⚠ CSR concatenation: the values append, and the offsets append SHIFTED by this queue's current
    // extent. Every offset array starts at 0, so `other`'s leading 0 is skipped rather than added.
    const auto shift_append = [](std::vector<std::int64_t>& into,
                                 const std::vector<std::int64_t>& from) {
        const std::int64_t base = into.back();
        for (std::size_t i = 1; i < from.size(); ++i) into.push_back(base + from[i]);
    };
    ref.insert(ref.end(), other.ref.begin(), other.ref.end());
    start.insert(start.end(), other.start.begin(), other.start.end());
    end.insert(end.end(), other.end.begin(), other.end.end());
    align_strand.insert(align_strand.end(), other.align_strand.begin(), other.align_strand.end());
    sj_strand.insert(sj_strand.end(), other.sj_strand.begin(), other.sj_strand.end());
    shift_append(observed_intron_offsets, other.observed_intron_offsets);
    observed_introns.insert(observed_introns.end(),
                            other.observed_introns.begin(), other.observed_introns.end());
    shift_append(hypothesis_offsets, other.hypothesis_offsets);
    hypothesis_sj_strand.insert(hypothesis_sj_strand.end(),
                                other.hypothesis_sj_strand.begin(), other.hypothesis_sj_strand.end());
    shift_append(hypothesis_intron_offsets, other.hypothesis_intron_offsets);
    hypothesis_introns.insert(hypothesis_introns.end(),
                              other.hypothesis_introns.begin(), other.hypothesis_introns.end());
    shift_append(hypothesis_t_offsets, other.hypothesis_t_offsets);
    hypothesis_t.insert(hypothesis_t.end(), other.hypothesis_t.begin(), other.hypothesis_t.end());
}

namespace {

/// Python's sequence comparison, exactly: element-wise, and a PREFIX sorts BEFORE the longer sequence it
/// is a prefix of. Returns -1, 0 or 1.
///
/// ⚠ A flat run of (start, end) PAIRS compares identically to Python's tuple-of-pairs, because a pair has
/// fixed length 2 -- so there is no need to compare pairwise and the same helper serves both the intron
/// runs and the supporting-transcript runs.
inline int lex_compare(const std::int64_t* a, std::size_t na,
                       const std::int64_t* b, std::size_t nb) noexcept
{
    const std::size_t common = std::min(na, nb);
    for (std::size_t i = 0; i < common; ++i) {
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    }
    if (na == nb) return 0;
    return na < nb ? -1 : 1;
}

inline int scalar_compare(std::int64_t a, std::int64_t b) noexcept {
    if (a == b) return 0;
    return a < b ? -1 : 1;
}

}  // namespace

void DeferredFragments::canonicalise() {
    const std::size_t n = size();
    if (n < 2) return;

    // ⭐ The specification's key, in the specification's order:
    //     (ref, start, end, align_strand, sj_strand, observed_introns, hypotheses)
    // where `hypotheses` is the sequence of (introns, sj_strand, supporting_t) triples. Compared as Python
    // compares those tuples -- see `lex_compare` for the prefix rule, which is the part a hand-rolled
    // comparator gets wrong: a one-intron path is NOT equal to a two-intron path that starts with it.
    const auto record_less = [this](std::size_t a, std::size_t b) {
        for (const std::vector<std::int64_t>* column : {&ref, &start, &end, &align_strand, &sj_strand}) {
            const int c = scalar_compare((*column)[a], (*column)[b]);
            if (c != 0) return c < 0;
        }
        const auto observed_run = [this](std::size_t i) {
            const std::size_t lo = static_cast<std::size_t>(observed_intron_offsets[i]);
            const std::size_t hi = static_cast<std::size_t>(observed_intron_offsets[i + 1]);
            return std::pair{observed_introns.data() + 2 * lo, 2 * (hi - lo)};
        };
        const auto [a_obs, a_n_obs] = observed_run(a);
        const auto [b_obs, b_n_obs] = observed_run(b);
        const int c_obs = lex_compare(a_obs, a_n_obs, b_obs, b_n_obs);
        if (c_obs != 0) return c_obs < 0;

        const std::size_t a_h0 = static_cast<std::size_t>(hypothesis_offsets[a]);
        const std::size_t a_h1 = static_cast<std::size_t>(hypothesis_offsets[a + 1]);
        const std::size_t b_h0 = static_cast<std::size_t>(hypothesis_offsets[b]);
        const std::size_t b_h1 = static_cast<std::size_t>(hypothesis_offsets[b + 1]);
        const std::size_t common = std::min(a_h1 - a_h0, b_h1 - b_h0);
        for (std::size_t k = 0; k < common; ++k) {
            const std::size_t ga = a_h0 + k, gb = b_h0 + k;
            const auto intron_run = [this](std::size_t g) {
                const std::size_t lo = static_cast<std::size_t>(hypothesis_intron_offsets[g]);
                const std::size_t hi = static_cast<std::size_t>(hypothesis_intron_offsets[g + 1]);
                return std::pair{hypothesis_introns.data() + 2 * lo, 2 * (hi - lo)};
            };
            const auto [a_in, a_n_in] = intron_run(ga);
            const auto [b_in, b_n_in] = intron_run(gb);
            int c = lex_compare(a_in, a_n_in, b_in, b_n_in);
            if (c != 0) return c < 0;
            c = scalar_compare(hypothesis_sj_strand[ga], hypothesis_sj_strand[gb]);
            if (c != 0) return c < 0;
            const std::size_t a_t0 = static_cast<std::size_t>(hypothesis_t_offsets[ga]);
            const std::size_t a_t1 = static_cast<std::size_t>(hypothesis_t_offsets[ga + 1]);
            const std::size_t b_t0 = static_cast<std::size_t>(hypothesis_t_offsets[gb]);
            const std::size_t b_t1 = static_cast<std::size_t>(hypothesis_t_offsets[gb + 1]);
            c = lex_compare(hypothesis_t.data() + a_t0, a_t1 - a_t0,
                            hypothesis_t.data() + b_t0, b_t1 - b_t0);
            if (c != 0) return c < 0;
        }
        return (a_h1 - a_h0) < (b_h1 - b_h0);
    };

    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    // ⚠ Not a stable sort, and it does not need to be: two records that tie on the whole key are
    // byte-identical records, so no permutation of them is observable in the exported arrays.
    std::sort(order.begin(), order.end(), record_less);

    bool identity = true;
    for (std::size_t i = 0; i < n && identity; ++i) identity = (order[i] == i);
    if (identity) return;  // keeps a second export free, which is what makes this idempotent in practice

    DeferredFragments out;
    out.ref.reserve(n);
    out.start.reserve(n);
    out.end.reserve(n);
    out.align_strand.reserve(n);
    out.sj_strand.reserve(n);
    out.observed_intron_offsets.reserve(n + 1);
    out.observed_introns.reserve(observed_introns.size());
    out.hypothesis_offsets.reserve(n + 1);
    out.hypothesis_sj_strand.reserve(hypothesis_sj_strand.size());
    out.hypothesis_intron_offsets.reserve(hypothesis_sj_strand.size() + 1);
    out.hypothesis_introns.reserve(hypothesis_introns.size());
    out.hypothesis_t_offsets.reserve(hypothesis_sj_strand.size() + 1);
    out.hypothesis_t.reserve(hypothesis_t.size());

    for (const std::size_t i : order) {
        out.ref.push_back(ref[i]);
        out.start.push_back(start[i]);
        out.end.push_back(end[i]);
        out.align_strand.push_back(align_strand[i]);
        out.sj_strand.push_back(sj_strand[i]);
        for (std::size_t p = static_cast<std::size_t>(observed_intron_offsets[i]);
             p < static_cast<std::size_t>(observed_intron_offsets[i + 1]); ++p) {
            out.observed_introns.push_back(observed_introns[2 * p]);
            out.observed_introns.push_back(observed_introns[2 * p + 1]);
        }
        out.observed_intron_offsets.push_back(
            static_cast<std::int64_t>(out.observed_introns.size() / 2));
        for (std::size_t g = static_cast<std::size_t>(hypothesis_offsets[i]);
             g < static_cast<std::size_t>(hypothesis_offsets[i + 1]); ++g) {
            out.hypothesis_sj_strand.push_back(hypothesis_sj_strand[g]);
            for (std::size_t p = static_cast<std::size_t>(hypothesis_intron_offsets[g]);
                 p < static_cast<std::size_t>(hypothesis_intron_offsets[g + 1]); ++p) {
                out.hypothesis_introns.push_back(hypothesis_introns[2 * p]);
                out.hypothesis_introns.push_back(hypothesis_introns[2 * p + 1]);
            }
            out.hypothesis_intron_offsets.push_back(
                static_cast<std::int64_t>(out.hypothesis_introns.size() / 2));
            for (std::size_t t = static_cast<std::size_t>(hypothesis_t_offsets[g]);
                 t < static_cast<std::size_t>(hypothesis_t_offsets[g + 1]); ++t) {
                out.hypothesis_t.push_back(hypothesis_t[t]);
            }
            out.hypothesis_t_offsets.push_back(static_cast<std::int64_t>(out.hypothesis_t.size()));
        }
        out.hypothesis_offsets.push_back(static_cast<std::int64_t>(out.hypothesis_sj_strand.size()));
    }
    *this = std::move(out);
}

// ============================================================================
// construction
// ============================================================================

Accumulator::Accumulator(std::vector<std::int64_t> cuts,
                         std::vector<std::uint8_t> region_types,
                         int max_length,
                         std::int32_t ref_id)
    : cuts_(std::move(cuts)), ref_id_(ref_id), max_length_(max_length)
{
    if (ref_id_ < 0) {
        throw std::invalid_argument(
            "accumulator: ref_id must be >= 0, got " + std::to_string(ref_id_) +
            ". It is stamped into every deferred record and the second pass replays those onto that "
            "reference's cut axis, so there is no such thing as a fragment held for no reference.");
    }
    if (max_length_ < 1) {
        throw std::invalid_argument(
            "accumulator: max_length must be >= 1, got " + std::to_string(max_length_) +
            ". It is the fragment-length limit applied to L as well as the pool-histogram width, so at 0 "
            "every real fragment is dropped as too long and the whole tally is silently empty.");
    }
    for (std::size_t i = 1; i < cuts_.size(); ++i) {
        if (cuts_[i] <= cuts_[i - 1]) {
            throw std::invalid_argument(
                "accumulator: cuts must strictly increase, but cuts[" + std::to_string(i) + "] = " +
                std::to_string(cuts_[i]) + " <= cuts[" + std::to_string(i - 1) + "] = " +
                std::to_string(cuts_[i - 1]));
        }
    }

    // A reference contributing c cuts owns c-1 regions and c-2 interior boundaries; one contributing fewer than
    // two cuts owns neither, which is legal and deposits nothing.
    const std::size_t n_regions = cuts_.size() >= 2 ? cuts_.size() - 1 : 0;
    const std::size_t n_boundaries = cuts_.size() >= 2 ? cuts_.size() - 2 : 0;
    regions_.assign(n_regions, Region{});
    boundaries_.assign(n_boundaries, Boundary{});
    region_start_count_.assign(n_regions, 0u);

    // ⚠ REQUIRED whenever this reference owns a region, and it throws rather than quietly skipping the
    // length pools. The shipped accumulator disabled pooling silently on an empty type array, which is a
    // whole output going missing with nothing to notice it by; and the specification has no such state at
    // all -- `Partition.from_cuts` always materialises a type per region -- so a C++ mode the spec cannot
    // express could only ever be a way to disagree with it.
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
    // ⭐ TRAPS: a-purity-filter-is-a-length-filter — allocated ALWAYS, unlike pool_lengths_ which is empty when a reference has no region types.
    // The unconditional histogram does not depend on region typing -- a fragment has a length whether or
    // not its region can be classified -- and an anchor that silently vanished on an untyped reference
    // would be exactly the kind of conditioning this row exists to remove.
    deposited_lengths_.assign(static_cast<std::size_t>(max_length_) + 1, 0u);
}

void Accumulator::set_junctions(std::vector<std::int32_t> offsets,
                                std::vector<std::int32_t> boundary_right,
                                std::vector<std::int8_t>  sj_strand)
{
    if (!offsets.empty() && offsets.size() != cuts_.size() + 1) {
        throw std::invalid_argument(
            "accumulator: junction CSR offsets must have length n_cuts + 1 = " +
            std::to_string(cuts_.size() + 1) + ", got " + std::to_string(offsets.size()));
    }
    if (boundary_right.size() != sj_strand.size()) {
        throw std::invalid_argument(
            "accumulator: junction boundary_right has " + std::to_string(boundary_right.size()) +
            " entries but sj_strand has " + std::to_string(sj_strand.size()));
    }
    if (!offsets.empty() && static_cast<std::size_t>(offsets.back()) != boundary_right.size()) {
        throw std::invalid_argument(
            "accumulator: junction CSR ends at " + std::to_string(offsets.back()) +
            " but there are " + std::to_string(boundary_right.size()) + " junctions");
    }
    sj_offsets_      = std::move(offsets);
    sj_boundary_right_ = std::move(boundary_right);
    sj_strand_       = std::move(sj_strand);
    junctions_.assign(sj_boundary_right_.size(), JunctionEdge{});
}

// ============================================================================
// locating things on the cut axis
// ============================================================================

std::int64_t Accumulator::region_of_pos(std::int64_t position) const noexcept {
    if (regions_.empty()) return -1;
    const auto it = std::upper_bound(cuts_.begin(), cuts_.end(), position);
    const std::int64_t region = (it - cuts_.begin()) - 1;
    return std::min<std::int64_t>(std::max<std::int64_t>(region, 0),
                                  static_cast<std::int64_t>(regions_.size()) - 1);
}

std::int64_t Accumulator::exact_cut(std::int64_t position) const noexcept {
    const auto it = std::lower_bound(cuts_.begin(), cuts_.end(), position);
    if (it == cuts_.end() || *it != position) return -1;
    return it - cuts_.begin();
}

std::int64_t Accumulator::sj_edge_id(std::int64_t intron_start,
                                     std::int64_t intron_end,
                                     std::int32_t sj_strand) const noexcept
{
    if (sj_offsets_.empty()) return -1;
    // ⭐ Every annotated intron has BOTH endpoints as partition cuts (measured: 404,168 of 404,168), so
    // "is this intron annotated?" reduces to the binary search the deposit already performs. If the start
    // is not a cut the table is never consulted -- and 70.4 % of cuts are not a donor at all.
    const std::int64_t donor = exact_cut(intron_start);
    if (donor < 0) return -1;
    const std::int64_t acceptor = exact_cut(intron_end);
    if (acceptor < 0) return -1;

    // ⚠ The strand filter applies only when the observed strand is DEFINITE. NONE means the aligner wrote
    // no motif tag at all (STAR writes XS, minimap2 ts, some write neither), so on such a BAM every
    // spliced fragment arrives with NONE -- demanding a strand there would delete 100 % of that aligner's
    // annotated junctions. AMBIGUOUS never reaches here; `deposit` rejects the whole fragment's splices.
    const bool definite = (sj_strand == STRAND_POS || sj_strand == STRAND_NEG);
    const std::int32_t lo = sj_offsets_[static_cast<std::size_t>(donor)];
    const std::int32_t hi = sj_offsets_[static_cast<std::size_t>(donor) + 1];
    for (std::int32_t k = lo; k < hi; ++k) {  // one to three iterations at human scale
        if (sj_boundary_right_[static_cast<std::size_t>(k)] != acceptor) continue;
        if (definite && sj_strand_[static_cast<std::size_t>(k)] != sj_strand) continue;
        return k;
    }
    return -1;
}

// ============================================================================
// the deposit
// ============================================================================

namespace {

/// The introns as a sorted, DISJOINT set clipped to [start, end); returns how many were absorbed.
///
/// ⚠ Overlapping AND ABUTTING introns both merge. Overlapping ones are contradictory observations of one
/// molecule -- a real BAM produces them when the mates disagree about an acceptor. Abutting ones imply a
/// zero-length exon, which is physically impossible, so no single molecule can legitimately use both.
///
/// ⚠ Sort the RAW pairs and clip inside the loop, in that order, to match the spec exactly.
std::int64_t normalise_introns(const IntronBlock* observed,
                               std::size_t n_observed,
                               const IntronBlock* implied,
                               std::size_t n_implied,
                               std::int64_t start,
                               std::int64_t end,
                               std::vector<std::pair<std::int64_t, std::int64_t>>& out)
{
    // ⚠ Both ranges are emplaced BEFORE the sort, so the union needs no second scratch vector and the
    // caller never has to think about which list came first.
    out.clear();
    out.reserve(n_observed + n_implied);
    for (std::size_t i = 0; i < n_observed; ++i) {
        out.emplace_back(static_cast<std::int64_t>(observed[i].start),
                         static_cast<std::int64_t>(observed[i].end));
    }
    for (std::size_t i = 0; i < n_implied; ++i) {
        out.emplace_back(static_cast<std::int64_t>(implied[i].start),
                         static_cast<std::int64_t>(implied[i].end));
    }
    std::sort(out.begin(), out.end());

    std::int64_t absorbed = 0;
    std::size_t kept = 0;
    for (std::size_t i = 0; i < out.size(); ++i) {
        const std::int64_t intron_start = std::max(out[i].first, start);
        const std::int64_t intron_end   = std::min(out[i].second, end);
        if (intron_end <= intron_start) {  // zero-length, or entirely outside the fragment
            ++absorbed;
            continue;
        }
        if (kept > 0 && intron_start <= out[kept - 1].second) {
            out[kept - 1].second = std::max(out[kept - 1].second, intron_end);
            ++absorbed;
            continue;
        }
        out[kept++] = {intron_start, intron_end};  // kept <= i, so this never clobbers an unread slot
    }
    out.resize(kept);
    return absorbed;
}

/// The path's contiguous genomic segments: [start, end) with the (normalised) introns cut out.
void build_segments(std::int64_t start,
                    std::int64_t end,
                    const std::vector<std::pair<std::int64_t, std::int64_t>>& introns,
                    std::vector<std::pair<std::int64_t, std::int64_t>>& out)
{
    out.clear();
    std::int64_t cursor = start;
    for (const auto& [intron_start, intron_end] : introns) {
        if (intron_start > cursor) out.emplace_back(cursor, intron_start);
        cursor = intron_end;
    }
    if (end > cursor) out.emplace_back(cursor, end);
}

}  // namespace

std::int64_t Accumulator::hypothesis_length(const OfferedFragment& fragment,
                                           const GapHypothesis& hypothesis,
                                           std::int64_t start,
                                           std::int64_t end,
                                           DepositScratch& scratch,
                                           std::int64_t* absorbed) const
{
    // ⚠ ONE definition of L: the total length of the path's segments. Deriving it any other way invites
    // two formulas for one quantity, and the obvious second formula -- (end - start) - sum(intron) --
    // disagrees by up to 1.5x on overlapping introns and goes NEGATIVE on a wide overlap.
    //
    // The observed introns and the implied ones are UNIONED. They are disjoint by construction -- the
    // enumeration never searches a gap the CIGAR already explained -- so this is a union and not a merge,
    // and `normalise_introns` sorts and clips it either way.
    *absorbed = normalise_introns(fragment.observed_introns, fragment.n_observed_introns,
                                  hypothesis.introns, hypothesis.n_introns,
                                  start, end, scratch.introns);
    build_segments(start, end, scratch.introns, scratch.segments);
    std::int64_t length = 0;
    for (const auto& [a, b] : scratch.segments) length += b - a;
    return length;
}

DepositOutcome Accumulator::deposit(const OfferedFragment& fragment, DepositScratch& scratch) {
    // ⚠ FIRST, before any geometry: the strand is a property of the fragment alone, and one with no
    // single genome strand has no column in any bank. Booking it anyway would credit it to the WRONG
    // strand; rejecting it here is what lets the loss be counted instead of vanishing.
    //
    // ⚠ And it must win over the deferral, because every fragment counts exactly ONCE and a fragment can
    // be both. The queue sizes the population the SECOND PASS CAN RECOVER, and a fragment with no genome
    // strand is not recoverable -- that pass resolves which PATH, not which strand the read aligned to.
    const int column = strand_column(fragment.align_strand);
    if (column < 0) {
        ++counters_.dropped_strand_undefined;
        return DepositOutcome::kStrandUndefined;
    }
    if (cuts_.size() < 2) {
        ++counters_.dropped_empty;
        return DepositOutcome::kEmpty;
    }

    // Clip to the reference. L is the CLIPPED length, so the placement count stays consistent -- and the
    // clip must precede arbitration, because a hypothesis is filtered on its L.
    const std::int64_t start = std::max(fragment.start, cuts_.front());
    const std::int64_t end   = std::min(fragment.end,   cuts_.back());
    if (end <= start) {
        ++counters_.dropped_empty;
        return DepositOutcome::kEmpty;
    }

    // ── arbitration: which hypotheses survive, and is exactly one left? ───────────────────────────
    //
    // ⭐ Short-read chemistry does not sequence molecules past `max_length_` -- the same statement that
    // makes kTooLong a rejection -- so a hypothesis implying a longer L is not a molecule this library
    // contains. Applied to the UNSPLICED hypothesis this is exactly "a fragment whose genomic span
    // exceeds the limit must be RNA": that hypothesis's L IS the span. There is no second rule.
    // ⚠ Unless the filter would empty the set, in which case the survivors stand and the ordinary
    // kTooLong rejection counts them, as it did before any of this.
    // ⛔⛔ EVERY FRAGMENT HAS AT LEAST ONE HYPOTHESIS -- "cut nothing beyond what was sequenced" -- and
    // the executable specification makes that its DEFAULT (`UNSPLICED_ONLY`), not an option:
    // "the degenerate case is the general case, not a branch". A caller offering an EMPTY set is
    // asking the same question, so answer it the same way rather than crashing.
    //
    // ⛔ Without this, an empty set walked straight past the `survivors.size() > 1` deferral and
    // dereferenced `survivors.front()` on an EMPTY vector, indexing a NULL `hypotheses` -- a hard
    // segfault, and one nothing could reach from the scanner because the scanner always offers the
    // genomic path. It cost a session to find precisely because it was unreachable from production.
    static const GapHypothesis kUnsplicedOnly{nullptr, 0, STRAND_NONE, nullptr, 0};
    OfferedFragment offered_or_default = fragment;
    if (fragment.n_hypotheses == 0) {
        offered_or_default.hypotheses = &kUnsplicedOnly;
        offered_or_default.n_hypotheses = 1;
    }
    const OfferedFragment& arbitrated = offered_or_default;

    auto& survivors = scratch.survivors;
    survivors.clear();
    bool any_spliced_hypothesis = false;
    for (std::size_t h = 0; h < arbitrated.n_hypotheses; ++h) {
        std::int64_t absorbed = 0;
        const std::int64_t candidate_length =
            hypothesis_length(arbitrated, arbitrated.hypotheses[h], start, end, scratch, &absorbed);
        any_spliced_hypothesis |= !arbitrated.hypotheses[h].is_unspliced();
        survivors.push_back({h, candidate_length, absorbed});
    }
    const std::size_t n_offered = survivors.size();
    survivors.erase(std::remove_if(survivors.begin(), survivors.end(),
                                   [this](const ScoredHypothesis& s) {
                                       return s.length > max_length_;
                                   }),
                    survivors.end());
    if (survivors.empty()) {
        for (std::size_t h = 0; h < n_offered; ++h) {
            std::int64_t absorbed = 0;
            const std::int64_t candidate_length =
                hypothesis_length(arbitrated, arbitrated.hypotheses[h], start, end, scratch, &absorbed);
            survivors.push_back({h, candidate_length, absorbed});
        }
    }
    if (any_spliced_hypothesis) record_gap_resolution(arbitrated, survivors);

    if (survivors.size() > 1) {
        deferred_.append(arbitrated, ref_id_, start, end);
        ++counters_.deferred_undetermined_gap;
        return DepositOutcome::kDeferred;
    }

    // The single survivor. ⚠ Re-normalised rather than cached per hypothesis: one extra normalise on the
    // winner is cheaper than carrying a normalised list for every candidate, and it keeps ONE code path
    // from a hypothesis to its introns.
    const GapHypothesis& chosen = arbitrated.hypotheses[survivors.front().index];
    std::int64_t absorbed = 0;
    const std::int64_t length =
        hypothesis_length(fragment, chosen, start, end, scratch, &absorbed);

    // ⭐ An implied strand is used ONLY when no motif was sequenced anywhere on this fragment. An observed
    // motif is evidence; an implied strand is an inference from the annotation, and mixing an inference
    // into an observation is how `primary` went wrong.
    const std::int32_t sj_strand =
        fragment.sj_strand == STRAND_NONE ? chosen.sj_strand : fragment.sj_strand;

    if (length <= 0) {
        ++counters_.dropped_empty;
        return DepositOutcome::kEmpty;
    }
    if (length > max_length_) {
        ++counters_.dropped_too_long;
        return DepositOutcome::kTooLong;
    }
    counters_.introns_absorbed += absorbed;

    // ── which annotated junctions does this path use? this also picks the boundary bank ────────────────
    // ⚠ Resolved BEFORE the crossing loop, because `spliced` chooses which bank the crossings land in.
    // ⭐ Resolved PER INTRON POSITION into `sj_id_at_gap`, with -1 where the annotation has no such
    // junction, because the conserved mass needs to know WHICH of a block's two ends is a junction
    // boundary. A filtered list cannot answer that -- dropping the unannotated entries destroys the
    // alignment between intron `i` and the gap between blocks `i` and `i+1`.
    auto& sj_ids      = scratch.sj_ids;
    auto& sj_id_at_gap = scratch.sj_id_at_gap;
    sj_ids.clear();
    sj_id_at_gap.assign(scratch.introns.size(), -1);
    if (sj_strand == STRAND_AMBIGUOUS) {
        // The motif tag is read once per RECORD, so AMBIGUOUS means the mates DISAGREED about one
        // molecule: contradictory evidence, not missing evidence. Trust no splice, and count it on its own
        // denominator -- folding it into `unannotated_introns` would poison the one metric whose job is
        // measuring annotation coverage.
        ++counters_.contradictory_sj_strand;
    } else {
        for (std::size_t gap = 0; gap < scratch.introns.size(); ++gap) {
            const auto& [intron_start, intron_end] = scratch.introns[gap];
            const std::int64_t id = sj_edge_id(intron_start, intron_end, sj_strand);
            if (id >= 0) {
                sj_id_at_gap[gap] = static_cast<std::int32_t>(id);
                sj_ids.push_back(static_cast<std::int32_t>(id));
            }
        }
        counters_.unannotated_introns +=
            static_cast<std::int64_t>(scratch.introns.size()) - static_cast<std::int64_t>(sj_ids.size());
    }
    const bool spliced = !sj_ids.empty();

    // ⚠ The path's own first and last COVERED base, not the fragment's extent. With a leading intron the
    // molecule does not begin at `start`: introns [(150,480)] over [150,500) has its whole path in
    // [480,500), a different region, and using the extent would credit a region it never touches.
    const std::int64_t first_base = scratch.segments.front().first;
    const std::int64_t last_base  = scratch.segments.back().second - 1;

    const std::int64_t first_region = region_of_pos(first_base);
    region_start_count_[static_cast<std::size_t>(first_region)] += 1u;
    // ⭐ TRAPS: a-purity-filter-is-a-length-filter — incremented HERE -- beside the start count and the DEPOSITED counter -- so all three
    // describe one population by construction rather than by agreement. `length` is already clipped to
    // the reference and gated by the length limit above.
    deposited_lengths_[static_cast<std::size_t>(length)] += 1u;
    ++counters_.deposited;

    // ── crossings, per contiguous SEGMENT of the path ─────────────────────────────────────────────
    // A boundary is crossed iff it lies strictly inside a segment, so per segment the crossed boundaries are a
    // contiguous index range and no container is needed. A region is SPANNED iff ONE segment crosses both of
    // its boundaries -- not merely "both boundaries crossed", which would count a region the fragment JUMPS OVER,
    // whose two boundaries are touched by the two flanking segments from opposite sides.
    //
    // ⚠ 0 at L == 1: a length-1 molecule cannot cross a 0-bp boundary, and 1/(L-1) would divide by zero. Its
    // residue is the schema's only count/density co-support violation -- an L == 1 path on an annotated
    // junction books a count against density 0, which is correct.
    const double inv_boundary = length >= 2 ? 1.0 / static_cast<double>(length - 1) : 0.0;
    const std::size_t   col          = static_cast<std::size_t>(column);

    std::int64_t n_crossed = 0;
    std::int64_t sole_boundary = -1;
    for (std::size_t block = 0; block < scratch.segments.size(); ++block) {
        const auto& [seg_start, seg_end] = scratch.segments[block];
        const std::int64_t first =
            std::upper_bound(cuts_.begin(), cuts_.end(), seg_start) - cuts_.begin();
        const std::int64_t last =
            std::lower_bound(cuts_.begin(), cuts_.end(), seg_end) - cuts_.begin();

        for (std::int64_t boundary_idx = first; boundary_idx < last; ++boundary_idx) {
            // ⭐ TWO channels on the spliced bank, FOUR on the unspliced one, and the asymmetry is
            // the design: a spliced crossing is certified RNA, nothing deconvolves it, so its length
            // moments have no consumer and are not stored.
            Boundary& boundary = boundaries_[static_cast<std::size_t>(boundary_idx - 1)];
            if (spliced) {
                boundary.spliced_count[col] += 1u;
            } else {
                boundary.unspliced_count[col] += 1u;
                boundary.unspliced_inv_length_sum += inv_boundary;
            }
        }
        // ── ⭐⭐⭐ THE CONSERVED MASS, per SLICE, over ONE BOUNDARY SET ────────────────────────────
        // The crossed cuts split this block into `last - first + 1` slices. Each slice's
        // `slice_len / length` is shared EQUALLY between the objects that bound it, so every bounded
        // slice disposes of exactly its own bases:  Sum over the fragment = Sum slice_len / length = 1.
        //
        // ⭐⭐ A JUNCTION IS A BOUNDARY EXACTLY LIKE A BOUNDARY, and that is the whole rule. A block's
        // interior boundaries are the boundaries it crosses; its two ENDS are boundaries too whenever the
        // intron there resolved to an annotated junction. So a fragment's 1.0 is shared across every
        // object it crosses -- boundaries and junctions together -- rather than boundaries first and junctions
        // only with whatever is left over.
        //
        // ⛔ The predecessor gave a boundary-crossing block's bases entirely to boundaries, so a junction whose
        // two flanking blocks both crossed a boundary received NOTHING while `sj_count` credited it. That
        // still conserved -- the total was 1.0 -- but it is not a sharing.
        //
        // ⭐ Coverage-weighted, NOT `1/K`. Both conserve; only this one says WHERE the fragment sat, and
        // only this one is expressible per base -- which is how the two are told apart at all.
        // ⚠ An UNSPLICED path has no junction boundaries, so this reduces to the previous rule exactly
        // and `unspliced_mass` is byte-identical. A single block with no boundary and no annotated junction
        // is bounded by nothing and deposits nothing: for a one-block path that is the CONTAINED case,
        // already whole in `contained_count`; for a multi-block one it is an unannotated intron's block.
        {
            const std::int32_t left_junction =
                block > 0 ? sj_id_at_gap[block - 1] : -1;
            const std::int32_t right_junction =
                block < sj_id_at_gap.size() ? sj_id_at_gap[block] : -1;
            const std::int64_t n_slices = last - first + 1;
            for (std::int64_t i = 0; i < n_slices; ++i) {
                const std::int64_t lo = (i == 0) ? seg_start : cuts_[static_cast<std::size_t>(first + i - 1)];
                const std::int64_t hi = (i == n_slices - 1) ? seg_end : cuts_[static_cast<std::size_t>(first + i)];
                const std::int64_t left_boundary  = (i > 0) ? first + i - 1 : -1;
                const std::int64_t right_boundary = (i < n_slices - 1) ? first + i : -1;
                // The block's own ends are junction boundaries; its interior ends are boundaries. A slice
                // therefore has at most one boundary of each kind on each side, never both.
                const std::int32_t left_jid  = (left_boundary  < 0) ? left_junction  : -1;
                const std::int32_t right_jid = (right_boundary < 0) ? right_junction : -1;
                const std::int64_t n_bounds = (left_boundary >= 0 ? 1 : 0) + (right_boundary >= 0 ? 1 : 0)
                                            + (left_jid  >= 0 ? 1 : 0) + (right_jid  >= 0 ? 1 : 0);
                if (n_bounds == 0) continue;
                // ⭐ float64, deposited directly: the share is a coverage fraction in (0, 1] and needs
                // no fixed point. Conservation is 2.1e6x tighter than the 2^-32 grid it replaced.
                const double share = static_cast<double>(hi - lo)
                                   / (static_cast<double>(length) * static_cast<double>(n_bounds));
                for (const std::int64_t boundary_idx : {left_boundary, right_boundary}) {
                    if (boundary_idx < 0) continue;
                    Boundary& boundary = boundaries_[static_cast<std::size_t>(boundary_idx - 1)];
                    if (spliced) boundary.spliced_mass += share;
                    else         boundary.unspliced_mass += share;
                }
                for (const std::int32_t jid : {left_jid, right_jid}) {
                    // ⭐ `col` — the SAME genome-strand column `count` is deposited at below, so
                    // `mass[c] / count[c]` is a per-strand mean rather than a ratio across populations.
                    if (jid >= 0) junctions_[static_cast<std::size_t>(jid)].mass[col] += share;
                }
            }
        }
        if (last > first) {
            sole_boundary = (n_crossed == 0 && last - first == 1) ? first : -1;
            n_crossed += last - first;
        }
    }

    for (const std::int32_t id : sj_ids) {
        JunctionEdge& junction = junctions_[static_cast<std::size_t>(id)];
        junction.count[col] += 1u;
        junction.inv_length_sum += inv_boundary;
    }

    // ── contained: the WHOLE path lies inside ONE region ────────────────────────────────────────────
    // ⚠ Not merely "crossed no boundary". An unannotated intron can swallow every boundary between two blocks,
    // leaving a fragment that crosses nothing yet straddles two regions. Such a fragment deposits on NO
    // object but still increments region_start_count, so the loss is visible rather than silent.
    std::int64_t contained_region = -1;
    if (sj_ids.empty() && first_region == region_of_pos(last_base)) {
        contained_region = first_region;
        Region& region = regions_[static_cast<std::size_t>(contained_region)];
        // ⭐⭐⭐ THE RECIPROCAL-OPPORTUNITY DEPOSIT -- what makes this channel a DENSITY. A length-`w`
        // fragment contained in a region of length `ell` had `ell - w + 1` admissible start positions, so
        // `1/(ell - w + 1)` cancels the opportunity identically and E[SUM] = rho for ANY length
        // distribution. `1/L` does NOT cancel it: measured, that channel read 25.67 density units for
        // short fragments and 1.60 for long ones at the same true density. An BOUNDARY is the `ell -> 0`
        // limit, which is why `1/(L-1)` is right there and was wrong here.
        // ⚠ `A >= 1` is structural: the fragment IS contained, so `w <= ell`.
        const std::int64_t region_len =
            cuts_[static_cast<std::size_t>(contained_region) + 1] -
            cuts_[static_cast<std::size_t>(contained_region)];
        region.contained_count[col] += 1u;
        region.contained_inv_opportunity_sum += 1.0 / static_cast<double>(region_len - length + 1);
    }

    if (!pool_lengths_.empty()) {
        const std::int64_t pool = fragment_pool(spliced, contained_region, sole_boundary);
        if (pool >= 0) {
            pool_lengths_[static_cast<std::size_t>(pool) * (static_cast<std::size_t>(max_length_) + 1) +
                          static_cast<std::size_t>(length)] += 1;
        }
    }
    return DepositOutcome::kDeposited;
}

void Accumulator::record_gap_resolution(const OfferedFragment& fragment,
                                        const std::vector<ScoredHypothesis>& survivors) noexcept
{
    // ⭐ The ONE site, so the subclasses cannot drift out of sync. The caller has already established
    // that at least one hypothesis was not the unspliced one -- a fragment that never had a question to
    // answer is not in this population, and counting it would quietly make the umbrella's denominator
    // the whole library.
    std::int64_t n_spliced = 0;
    bool unspliced_survives = false;
    for (const ScoredHypothesis& s : survivors) {
        if (fragment.hypotheses[s.index].is_unspliced()) unspliced_survives = true;
        else ++n_spliced;
    }
    if (survivors.size() == 1) {
        // ⛔ The sole survivor is necessarily the SPLICED one, and there is no branch for the other case
        // because it cannot occur: the unspliced path is always the longest, so it can never be the last
        // one left. See `GapCensus`.
        ++gap_census_.resolved_spliced;
    } else if (!unspliced_survives) {
        ++gap_census_.deferred_which_introns;   // certified RNA; only the structure is open
    } else if (n_spliced == 1) {
        ++gap_census_.deferred_rna_or_gdna;     // one bit: was anything spliced at all?
    } else {
        ++gap_census_.deferred_both;
    }
}

std::int64_t Accumulator::fragment_pool(bool spliced,
                                        std::int64_t contained_region,
                                        std::int64_t sole_boundary) const noexcept
{
    // Priority, so that every pool stays pure: an OBSERVED splice is unambiguously RNA; a contained
    // fragment is typed by its region; a single-boundary crossing is a "splash" read typed by its two flanks.
    // Anything else -- an exonic contained fragment, a multi-boundary crossing -- is a mixture and enters
    // nothing.
    //
    // ⭐ DETERMINACY, NOT PROVENANCE: a fragment reaches here only when exactly ONE hypothesis survived,
    // so its L is not in doubt however it was arrived at. See the declaration for the measurement that
    // deleted the old `sj_implicit` bar.
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

std::int64_t Accumulator::length_under(const OfferedFragment& fragment,
                                       std::size_t hypothesis_index,
                                       DepositScratch& scratch) const
{
    // ⚠ The SAME clip `deposit` applies, and in the same order: L is measured after clipping to the
    // reference, so a scorer that skipped it would score a length the deposit can never produce.
    if (cuts_.size() < 2 || hypothesis_index >= fragment.n_hypotheses) return 0;
    const std::int64_t start = std::max(fragment.start, cuts_.front());
    const std::int64_t end   = std::min(fragment.end,   cuts_.back());
    if (end <= start) return 0;
    std::int64_t absorbed = 0;
    return hypothesis_length(fragment, fragment.hypotheses[hypothesis_index], start, end, scratch,
                             &absorbed);
}

const DeferredFragments& Accumulator::deferred_canonical() {
    deferred_.canonicalise();
    return deferred_;
}

// ============================================================================
// the per-worker merge
// ============================================================================

void Accumulator::merge_from(const Accumulator& other) {
    // ⚠ The cut arrays must match element-wise. A ref-id mismatch here once silently dropped 476,719 of
    // 476,732 fragments while every golden test passed, so this compares positions rather than sizes.
    //
    // ⭐ And now it can compare the ref id itself, which is what that bug was actually about — two
    // references with coincidentally identical cut axes would have passed the position check.
    if (ref_id_ != other.ref_id_) {
        throw std::invalid_argument(
            "accumulator: merge_from requires the same reference (this is ref " +
            std::to_string(ref_id_) + ", other is ref " + std::to_string(other.ref_id_) + ")");
    }
    if (cuts_ != other.cuts_) {
        throw std::invalid_argument(
            "accumulator: merge_from requires identical cut positions (this has " +
            std::to_string(cuts_.size()) + ", other has " + std::to_string(other.cuts_.size()) + ")");
    }
    if (junctions_.size() != other.junctions_.size()) {
        throw std::invalid_argument(
            "accumulator: merge_from requires the same junction bank (this has " +
            std::to_string(junctions_.size()) + ", other has " +
            std::to_string(other.junctions_.size()) + ")");
    }

    gap_census_.merge_from(other.gap_census_);
    deferred_.merge_from(other.deferred_);

    // Integer addition is associative, so the result is identical at any worker count, on any machine.
    for (std::size_t i = 0; i < regions_.size(); ++i) {
        for (std::size_t c = 0; c < kNStrandColumns; ++c) {
            regions_[i].contained_count[c]   += other.regions_[i].contained_count[c];
        }
        // ⚠ Outside the column loop — these have ONE value per region, not one per strand.
        regions_[i].contained_inv_opportunity_sum += other.regions_[i].contained_inv_opportunity_sum;
        region_start_count_[i] += other.region_start_count_[i];
    }
    for (std::size_t i = 0; i < deposited_lengths_.size(); ++i) {
        deposited_lengths_[i] += other.deposited_lengths_[i];
    }
    for (std::size_t i = 0; i < boundaries_.size(); ++i) {
        for (std::size_t c = 0; c < kNStrandColumns; ++c) {
            boundaries_[i].unspliced_count[c]   += other.boundaries_[i].unspliced_count[c];
            boundaries_[i].spliced_count[c]     += other.boundaries_[i].spliced_count[c];
        }
        // ⚠ Outside the column loop — these have ONE value per boundary, not one per strand.
        boundaries_[i].unspliced_inv_length_sum += other.boundaries_[i].unspliced_inv_length_sum;
        boundaries_[i].unspliced_mass += other.boundaries_[i].unspliced_mass;
        boundaries_[i].spliced_mass   += other.boundaries_[i].spliced_mass;
    }
    for (std::size_t i = 0; i < junctions_.size(); ++i) {
        for (std::size_t c = 0; c < kNStrandColumns; ++c) {
            junctions_[i].count[c]   += other.junctions_[i].count[c];
            // ⭐ INSIDE the column loop, unlike every other mass in this file: the junction mass is the
            // one that carries a strand. See `JunctionEdge::mass` for the premise that changed.
            junctions_[i].mass[c]    += other.junctions_[i].mass[c];
        }
        junctions_[i].inv_length_sum += other.junctions_[i].inv_length_sum;
    }
    // ⚠ Throws rather than skipping. A size mismatch means the two were built with different `max_length`,
    // and silently not merging would lose one side's pools entirely — the same class of defect as the
    // shipped accumulator disabling its pools on an empty type array.
    if (pool_lengths_.size() != other.pool_lengths_.size()) {
        throw std::invalid_argument(
            "accumulator: merge_from requires the same pool-histogram width (this has " +
            std::to_string(pool_lengths_.size()) + " bins, other has " +
            std::to_string(other.pool_lengths_.size()) + "); they were built with different max_length");
    }
    for (std::size_t i = 0; i < pool_lengths_.size(); ++i) {
        pool_lengths_[i] += other.pool_lengths_[i];
    }
    counters_.merge_from(other.counters_);
}

// ============================================================================
// AccumulatorSet
// ============================================================================

AccumulatorSet::AccumulatorSet(const std::int64_t* cut_positions,
                               std::size_t n_positions,
                               const std::int64_t* ref_cut_offsets,
                               std::size_t n_refs,
                               const std::uint8_t* region_types,
                               std::size_t n_region_types,
                               int max_length)
{
    if (ref_cut_offsets == nullptr || static_cast<std::size_t>(ref_cut_offsets[n_refs]) != n_positions) {
        throw std::invalid_argument(
            "accumulator set: ref_cut_offsets must end at n_positions = " + std::to_string(n_positions));
    }
    accs_.reserve(n_refs);

    // A reference contributing c cuts owns c-1 regions, so the region offset for reference f is
    // ref_cut_offsets[f] - (number of earlier references that own any region).
    std::size_t region_base = 0;
    for (std::size_t f = 0; f < n_refs; ++f) {
        const std::size_t lo = static_cast<std::size_t>(ref_cut_offsets[f]);
        const std::size_t hi = static_cast<std::size_t>(ref_cut_offsets[f + 1]);
        std::vector<std::int64_t> cuts(cut_positions + lo, cut_positions + hi);
        const std::size_t n_regions = hi - lo >= 2 ? hi - lo - 1 : 0;

        std::vector<std::uint8_t> types;
        if (region_types != nullptr && n_region_types > 0 && n_regions > 0) {
            if (region_base + n_regions > n_region_types) {
                throw std::invalid_argument(
                    "accumulator set: region_types has " + std::to_string(n_region_types) +
                    " entries but reference " + std::to_string(f) + " needs " +
                    std::to_string(region_base + n_regions));
            }
            types.assign(region_types + region_base, region_types + region_base + n_regions);
        }
        accs_.emplace_back(std::move(cuts), std::move(types), max_length,
                           static_cast<std::int32_t>(f));
        region_base += n_regions;
    }
}

void AccumulatorSet::set_junctions(const std::int64_t* offsets,
                                   std::size_t n_offsets,
                                   const std::int64_t* boundary_right,
                                   const std::int8_t* sj_strand,
                                   std::size_t n_junctions,
                                   const std::int64_t* ref_cut_offsets)
{
    if (offsets == nullptr || ref_cut_offsets == nullptr) {
        throw std::invalid_argument(
            "accumulator set: set_junctions needs both the CSR offsets and ref_cut_offsets");
    }
    const std::size_t n_cuts = static_cast<std::size_t>(ref_cut_offsets[accs_.size()]);
    if (n_offsets != n_cuts + 1) {
        throw std::invalid_argument(
            "accumulator set: junction CSR offsets must have length n_cuts + 1 = " +
            std::to_string(n_cuts + 1) + ", got " + std::to_string(n_offsets));
    }
    if (static_cast<std::size_t>(offsets[n_cuts]) != n_junctions) {
        throw std::invalid_argument(
            "accumulator set: the junction CSR ends at " + std::to_string(offsets[n_cuts]) +
            " but " + std::to_string(n_junctions) + " junctions were given");
    }

    std::vector<std::int32_t> ref_offsets, ref_boundary_right;
    std::vector<std::int8_t>  ref_strand;
    for (std::size_t f = 0; f < accs_.size(); ++f) {
        const std::int64_t c0 = ref_cut_offsets[f];
        const std::int64_t c1 = ref_cut_offsets[f + 1];
        if (c1 <= c0) continue;  // a reference with no cuts owns no regions and no junctions
        const std::int64_t j0 = offsets[c0];
        const std::int64_t j1 = offsets[c1];

        ref_offsets.clear();
        ref_offsets.reserve(static_cast<std::size_t>(c1 - c0) + 1);
        for (std::int64_t c = c0; c <= c1; ++c) {
            ref_offsets.push_back(static_cast<std::int32_t>(offsets[c] - j0));
        }
        ref_boundary_right.clear();
        ref_strand.clear();
        ref_boundary_right.reserve(static_cast<std::size_t>(j1 - j0));
        ref_strand.reserve(static_cast<std::size_t>(j1 - j0));
        for (std::int64_t k = j0; k < j1; ++k) {
            ref_boundary_right.push_back(static_cast<std::int32_t>(boundary_right[k] - c0));
            ref_strand.push_back(sj_strand[k]);
        }
        accs_[f].set_junctions(ref_offsets, ref_boundary_right, ref_strand);
    }
}

Accumulator& AccumulatorSet::at(std::int32_t ref_id) {
    if (ref_id < 0 || static_cast<std::size_t>(ref_id) >= accs_.size()) {
        throw std::out_of_range("accumulator set: ref_id " + std::to_string(ref_id) + " of " +
                                std::to_string(accs_.size()));
    }
    return accs_[static_cast<std::size_t>(ref_id)];
}

const Accumulator& AccumulatorSet::at(std::int32_t ref_id) const {
    if (ref_id < 0 || static_cast<std::size_t>(ref_id) >= accs_.size()) {
        throw std::out_of_range("accumulator set: ref_id " + std::to_string(ref_id) + " of " +
                                std::to_string(accs_.size()));
    }
    return accs_[static_cast<std::size_t>(ref_id)];
}

void AccumulatorSet::merge_from(const AccumulatorSet& other) {
    if (accs_.size() != other.accs_.size()) {
        throw std::invalid_argument(
            "accumulator set: merge_from requires the same reference count (this has " +
            std::to_string(accs_.size()) + ", other has " + std::to_string(other.accs_.size()) + ")");
    }
    for (std::size_t f = 0; f < accs_.size(); ++f) accs_[f].merge_from(other.accs_[f]);
}

}  // namespace rigel::accumulator
