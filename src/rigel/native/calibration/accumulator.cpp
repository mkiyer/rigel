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
    dropped_ambiguous_path   += other.dropped_ambiguous_path;
    unannotated_introns      += other.unannotated_introns;
    contradictory_sj_strand  += other.contradictory_sj_strand;
    sj_implicit_fragments    += other.sj_implicit_fragments;
    introns_absorbed         += other.introns_absorbed;
}

// ============================================================================
// construction
// ============================================================================

Accumulator::Accumulator(std::vector<std::int64_t> cuts,
                         std::vector<std::uint8_t> node_types,
                         int max_length)
    : cuts_(std::move(cuts)), max_length_(max_length)
{
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

    // A reference contributing c cuts owns c-1 nodes and c-2 interior lines; one contributing fewer than
    // two cuts owns neither, which is legal and deposits nothing.
    const std::size_t n_nodes = cuts_.size() >= 2 ? cuts_.size() - 1 : 0;
    const std::size_t n_edges = cuts_.size() >= 2 ? cuts_.size() - 2 : 0;
    nodes_.assign(n_nodes, Node{});
    edges_.assign(n_edges, ContiguousEdge{});
    node_start_count_.assign(n_nodes, 0u);

    // ⚠ REQUIRED whenever this reference owns a node, and it throws rather than quietly skipping the
    // length pools. The shipped accumulator disabled pooling silently on an empty type array, which is a
    // whole output going missing with nothing to notice it by; and the specification has no such state at
    // all -- `Partition.from_cuts` always materialises a type per node -- so a C++ mode the spec cannot
    // express could only ever be a way to disagree with it.
    if (n_nodes > 0 && node_types.empty()) {
        throw std::invalid_argument(
            "accumulator: node_types is empty but this reference has " + std::to_string(n_nodes) +
            " nodes. It types the fragment-length pools, so an empty array would silently emit no pools "
            "at all rather than fail.");
    }
    if (!node_types.empty()) {
        if (node_types.size() != n_nodes) {
            throw std::invalid_argument(
                "accumulator: node_types has " + std::to_string(node_types.size()) +
                " entries but this reference has " + std::to_string(n_nodes) + " nodes");
        }
        node_types_ = std::move(node_types);
    }
    pool_lengths_.assign(kNFragmentPools * (static_cast<std::size_t>(max_length_) + 1), 0);
}

void Accumulator::set_junctions(std::vector<std::int32_t> offsets,
                                std::vector<std::int32_t> acceptor_cut,
                                std::vector<std::int8_t>  sj_strand)
{
    if (!offsets.empty() && offsets.size() != cuts_.size() + 1) {
        throw std::invalid_argument(
            "accumulator: junction CSR offsets must have length n_cuts + 1 = " +
            std::to_string(cuts_.size() + 1) + ", got " + std::to_string(offsets.size()));
    }
    if (acceptor_cut.size() != sj_strand.size()) {
        throw std::invalid_argument(
            "accumulator: junction acceptor_cut has " + std::to_string(acceptor_cut.size()) +
            " entries but sj_strand has " + std::to_string(sj_strand.size()));
    }
    if (!offsets.empty() && static_cast<std::size_t>(offsets.back()) != acceptor_cut.size()) {
        throw std::invalid_argument(
            "accumulator: junction CSR ends at " + std::to_string(offsets.back()) +
            " but there are " + std::to_string(acceptor_cut.size()) + " junctions");
    }
    sj_offsets_      = std::move(offsets);
    sj_acceptor_cut_ = std::move(acceptor_cut);
    sj_strand_       = std::move(sj_strand);
    junctions_.assign(sj_acceptor_cut_.size(), JunctionEdge{});
}

// ============================================================================
// locating things on the cut axis
// ============================================================================

std::int64_t Accumulator::node_of_pos(std::int64_t position) const noexcept {
    if (nodes_.empty()) return -1;
    const auto it = std::upper_bound(cuts_.begin(), cuts_.end(), position);
    const std::int64_t node = (it - cuts_.begin()) - 1;
    return std::min<std::int64_t>(std::max<std::int64_t>(node, 0),
                                  static_cast<std::int64_t>(nodes_.size()) - 1);
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
        if (sj_acceptor_cut_[static_cast<std::size_t>(k)] != acceptor) continue;
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
std::int64_t normalise_introns(const IntronBlock* introns,
                               std::size_t n_introns,
                               std::int64_t start,
                               std::int64_t end,
                               std::vector<std::pair<std::int64_t, std::int64_t>>& out)
{
    out.clear();
    out.reserve(n_introns);
    for (std::size_t i = 0; i < n_introns; ++i) {
        out.emplace_back(static_cast<std::int64_t>(introns[i].start),
                         static_cast<std::int64_t>(introns[i].end));
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

DepositOutcome Accumulator::deposit(const FragmentPath& path, DepositScratch& scratch) {
    // ⚠ FIRST, before any geometry: the strand is a property of the fragment alone, and one with no
    // single genome strand has no column in any bank. Booking it anyway would credit it to the WRONG
    // strand; rejecting it here is what lets the loss be counted instead of vanishing.
    const int column = strand_column(path.align_strand);
    if (column < 0) {
        ++counters_.dropped_strand_undefined;
        return DepositOutcome::kStrandUndefined;
    }
    // ⚠ SECOND, and the order against the strand is part of the contract, because every fragment must
    // count exactly ONCE and a fragment can be both. It is filed under the strand, because that is which
    // denominator stays honest: `dropped_ambiguous_path` sizes the population the SECOND PASS CAN RECOVER,
    // and a fragment with no genome strand is not recoverable -- the second pass resolves which
    // transcript, not which strand the read aligned to.
    if (path.path_ambiguous) {
        ++counters_.dropped_ambiguous_path;
        return DepositOutcome::kAmbiguousPath;
    }
    if (cuts_.size() < 2) {
        ++counters_.dropped_empty;
        return DepositOutcome::kEmpty;
    }

    // Clip to the reference. L is the CLIPPED length, so the placement count stays consistent.
    const std::int64_t start = std::max(path.start, cuts_.front());
    const std::int64_t end   = std::min(path.end,   cuts_.back());
    if (end <= start) {
        ++counters_.dropped_empty;
        return DepositOutcome::kEmpty;
    }

    // ⚠ ONE definition of L: the total length of the path's segments. Deriving it any other way invites
    // two formulas for one quantity, and the obvious second formula -- (end - start) - sum(intron) --
    // disagrees by up to 1.5x on overlapping introns and goes NEGATIVE on a wide overlap.
    const std::int64_t absorbed =
        normalise_introns(path.introns, path.n_introns, start, end, scratch.introns);
    build_segments(start, end, scratch.introns, scratch.segments);
    std::int64_t length = 0;
    for (const auto& [a, b] : scratch.segments) length += b - a;

    if (length <= 0) {
        ++counters_.dropped_empty;
        return DepositOutcome::kEmpty;
    }
    if (length > max_length_) {
        ++counters_.dropped_too_long;
        return DepositOutcome::kTooLong;
    }
    counters_.introns_absorbed += absorbed;

    // ── which annotated junctions does this path use? this also picks the edge bank ────────────────
    // ⚠ Resolved BEFORE the crossing loop, because `spliced` chooses which bank the crossings land in.
    auto& sj_ids = scratch.sj_ids;
    sj_ids.clear();
    if (path.sj_strand == STRAND_AMBIGUOUS) {
        // The motif tag is read once per RECORD, so AMBIGUOUS means the mates DISAGREED about one
        // molecule: contradictory evidence, not missing evidence. Trust no splice, and count it on its own
        // denominator -- folding it into `unannotated_introns` would poison the one metric whose job is
        // measuring annotation coverage.
        ++counters_.contradictory_sj_strand;
    } else {
        for (const auto& [intron_start, intron_end] : scratch.introns) {
            const std::int64_t id = sj_edge_id(intron_start, intron_end, path.sj_strand);
            if (id >= 0) sj_ids.push_back(static_cast<std::int32_t>(id));
        }
        counters_.unannotated_introns +=
            static_cast<std::int64_t>(scratch.introns.size()) - static_cast<std::int64_t>(sj_ids.size());
    }
    const bool spliced = !sj_ids.empty();

    // ⚠ The path's own first and last COVERED base, not the fragment's extent. With a leading intron the
    // molecule does not begin at `start`: introns [(150,480)] over [150,500) has its whole path in
    // [480,500), a different node, and using the extent would credit a node it never touches.
    const std::int64_t first_base = scratch.segments.front().first;
    const std::int64_t last_base  = scratch.segments.back().second - 1;

    const std::int64_t first_node = node_of_pos(first_base);
    node_start_count_[static_cast<std::size_t>(first_node)] += 1u;
    ++counters_.deposited;
    if (path.sj_implicit) ++counters_.sj_implicit_fragments;

    // ── crossings, per contiguous SEGMENT of the path ─────────────────────────────────────────────
    // A line is crossed iff it lies strictly inside a segment, so per segment the crossed lines are a
    // contiguous index range and no container is needed. A node is SPANNED iff ONE segment crosses both of
    // its lines -- not merely "both lines crossed", which would count a node the fragment JUMPS OVER,
    // whose two lines are touched by the two flanking segments from opposite sides.
    //
    // ⚠ quantum_edge is 0 at L == 1: a length-1 molecule cannot cross a 0-bp line, and `inv_length_quantum`
    // would divide by zero. Its residue is the schema's only count/density co-support violation -- an
    // L == 1 path on an annotated junction books a count against density 0, which is correct.
    const std::uint64_t quantum_edge = length >= 2 ? inv_length_quantum(length - 1) : 0;
    const std::uint64_t quantum_node = inv_length_quantum(length);
    const std::size_t   col          = static_cast<std::size_t>(column);

    std::int64_t n_crossed = 0;
    std::int64_t sole_line = -1;
    for (const auto& [seg_start, seg_end] : scratch.segments) {
        const std::int64_t first =
            std::upper_bound(cuts_.begin(), cuts_.end(), seg_start) - cuts_.begin();
        const std::int64_t last =
            std::lower_bound(cuts_.begin(), cuts_.end(), seg_end) - cuts_.begin();

        for (std::int64_t line = first; line < last; ++line) {
            ContiguousEdge& edge = edges_[static_cast<std::size_t>(line - 1)];
            if (spliced) {
                edge.spliced_count[col] += 1u;
                edge.spliced_inv_length_sum[col] += quantum_edge;
                edge.spliced_length_sum[col] += static_cast<std::uint64_t>(length);
            } else {
                edge.unspliced_count[col] += 1u;
                edge.unspliced_inv_length_sum[col] += quantum_edge;
                edge.unspliced_length_sum[col] += static_cast<std::uint64_t>(length);
            }
        }
        for (std::int64_t line = first; line + 1 < last; ++line) {  // the node between two crossed lines
            Node& node = nodes_[static_cast<std::size_t>(line)];
            node.spanning_count[col] += 1u;
            node.spanning_inv_length_sum[col] += quantum_node;
            node.spanning_length_sum[col] += static_cast<std::uint64_t>(length);
        }
        if (last > first) {
            sole_line = (n_crossed == 0 && last - first == 1) ? first : -1;
            n_crossed += last - first;
        }
    }

    for (const std::int32_t id : sj_ids) {
        JunctionEdge& junction = junctions_[static_cast<std::size_t>(id)];
        junction.count[col] += 1u;
        junction.inv_length_sum[col] += quantum_edge;
        junction.length_sum[col] += static_cast<std::uint64_t>(length);
    }

    // ── contained: the WHOLE path lies inside ONE node ────────────────────────────────────────────
    // ⚠ Not merely "crossed no line". An unannotated intron can swallow every line between two blocks,
    // leaving a fragment that crosses nothing yet straddles two nodes. Such a fragment deposits on NO
    // object but still increments node_start_count, so the loss is visible rather than silent.
    std::int64_t contained_node = -1;
    if (sj_ids.empty() && first_node == node_of_pos(last_base)) {
        contained_node = first_node;
        Node& node = nodes_[static_cast<std::size_t>(contained_node)];
        node.contained_count[col] += 1u;
        node.contained_inv_length_sum[col] += quantum_node;
        node.contained_length_sum[col] += static_cast<std::uint64_t>(length);
    }

    if (!pool_lengths_.empty()) {
        const std::int64_t pool = fragment_pool(spliced, path.sj_implicit, contained_node, sole_line);
        if (pool >= 0) {
            pool_lengths_[static_cast<std::size_t>(pool) * (static_cast<std::size_t>(max_length_) + 1) +
                          static_cast<std::size_t>(length)] += 1;
        }
    }
    return DepositOutcome::kDeposited;
}

std::int64_t Accumulator::fragment_pool(bool spliced,
                                        bool sj_implicit,
                                        std::int64_t contained_node,
                                        std::int64_t sole_line) const noexcept
{
    // Priority, so that every pool stays pure: an OBSERVED splice is unambiguously RNA; a contained
    // fragment is typed by its node; a single-line crossing is a "splash" read typed by its two flanks.
    // Anything else -- an exonic contained fragment, a multi-line crossing -- is a mixture and enters
    // nothing. An IMPLICIT splice enters nothing either: it was never sequenced, so certifying it as RNA
    // would make the pure-RNA pool depend on the very length model it is used to fit.
    if (spliced) {
        return sj_implicit ? -1 : static_cast<std::int64_t>(FragmentPool::kRnaSpliced);
    }
    if (contained_node >= 0) {
        switch (node_types_[static_cast<std::size_t>(contained_node)]) {
            case kTypeIntergenic: return static_cast<std::int64_t>(FragmentPool::kDnaIntergenic);
            case kTypeIntron:     return static_cast<std::int64_t>(FragmentPool::kDnaIntronic);
            default:              return -1;  // exonic is a gDNA/RNA mixture, absent by design
        }
    }
    if (sole_line >= 1) {
        const std::uint8_t left  = node_types_[static_cast<std::size_t>(sole_line) - 1];
        const std::uint8_t right = node_types_[static_cast<std::size_t>(sole_line)];
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

// ============================================================================
// the per-worker merge
// ============================================================================

void Accumulator::merge_from(const Accumulator& other) {
    // ⚠ The cut arrays must match element-wise. A ref-id mismatch here once silently dropped 476,719 of
    // 476,732 fragments while every golden test passed, so this compares positions rather than sizes.
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

    // Integer addition is associative, so the result is identical at any worker count, on any machine.
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        for (std::size_t c = 0; c < kNStrandColumns; ++c) {
            nodes_[i].contained_count[c]   += other.nodes_[i].contained_count[c];
            nodes_[i].spanning_count[c]    += other.nodes_[i].spanning_count[c];
            nodes_[i].contained_inv_length_sum[c] += other.nodes_[i].contained_inv_length_sum[c];
            nodes_[i].spanning_inv_length_sum[c]  += other.nodes_[i].spanning_inv_length_sum[c];
            nodes_[i].contained_length_sum[c] += other.nodes_[i].contained_length_sum[c];
            nodes_[i].spanning_length_sum[c]  += other.nodes_[i].spanning_length_sum[c];
        }
        node_start_count_[i] += other.node_start_count_[i];
    }
    for (std::size_t i = 0; i < edges_.size(); ++i) {
        for (std::size_t c = 0; c < kNStrandColumns; ++c) {
            edges_[i].unspliced_count[c]   += other.edges_[i].unspliced_count[c];
            edges_[i].spliced_count[c]     += other.edges_[i].spliced_count[c];
            edges_[i].unspliced_inv_length_sum[c] += other.edges_[i].unspliced_inv_length_sum[c];
            edges_[i].spliced_inv_length_sum[c]   += other.edges_[i].spliced_inv_length_sum[c];
            edges_[i].unspliced_length_sum[c] += other.edges_[i].unspliced_length_sum[c];
            edges_[i].spliced_length_sum[c]   += other.edges_[i].spliced_length_sum[c];
        }
    }
    for (std::size_t i = 0; i < junctions_.size(); ++i) {
        for (std::size_t c = 0; c < kNStrandColumns; ++c) {
            junctions_[i].count[c]   += other.junctions_[i].count[c];
            junctions_[i].inv_length_sum[c] += other.junctions_[i].inv_length_sum[c];
            junctions_[i].length_sum[c] += other.junctions_[i].length_sum[c];
        }
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
                               const std::uint8_t* node_types,
                               std::size_t n_node_types,
                               int max_length)
{
    if (ref_cut_offsets == nullptr || static_cast<std::size_t>(ref_cut_offsets[n_refs]) != n_positions) {
        throw std::invalid_argument(
            "accumulator set: ref_cut_offsets must end at n_positions = " + std::to_string(n_positions));
    }
    accs_.reserve(n_refs);

    // A reference contributing c cuts owns c-1 nodes, so the node offset for reference f is
    // ref_cut_offsets[f] - (number of earlier references that own any node).
    std::size_t node_base = 0;
    for (std::size_t f = 0; f < n_refs; ++f) {
        const std::size_t lo = static_cast<std::size_t>(ref_cut_offsets[f]);
        const std::size_t hi = static_cast<std::size_t>(ref_cut_offsets[f + 1]);
        std::vector<std::int64_t> cuts(cut_positions + lo, cut_positions + hi);
        const std::size_t n_nodes = hi - lo >= 2 ? hi - lo - 1 : 0;

        std::vector<std::uint8_t> types;
        if (node_types != nullptr && n_node_types > 0 && n_nodes > 0) {
            if (node_base + n_nodes > n_node_types) {
                throw std::invalid_argument(
                    "accumulator set: node_types has " + std::to_string(n_node_types) +
                    " entries but reference " + std::to_string(f) + " needs " +
                    std::to_string(node_base + n_nodes));
            }
            types.assign(node_types + node_base, node_types + node_base + n_nodes);
        }
        accs_.emplace_back(std::move(cuts), std::move(types), max_length);
        node_base += n_nodes;
    }
}

void AccumulatorSet::set_junctions(const std::int64_t* offsets,
                                   std::size_t n_offsets,
                                   const std::int64_t* acceptor_cut,
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

    std::vector<std::int32_t> ref_offsets, ref_acceptor;
    std::vector<std::int8_t>  ref_strand;
    for (std::size_t f = 0; f < accs_.size(); ++f) {
        const std::int64_t c0 = ref_cut_offsets[f];
        const std::int64_t c1 = ref_cut_offsets[f + 1];
        if (c1 <= c0) continue;  // a reference with no cuts owns no nodes and no junctions
        const std::int64_t j0 = offsets[c0];
        const std::int64_t j1 = offsets[c1];

        ref_offsets.clear();
        ref_offsets.reserve(static_cast<std::size_t>(c1 - c0) + 1);
        for (std::int64_t c = c0; c <= c1; ++c) {
            ref_offsets.push_back(static_cast<std::int32_t>(offsets[c] - j0));
        }
        ref_acceptor.clear();
        ref_strand.clear();
        ref_acceptor.reserve(static_cast<std::size_t>(j1 - j0));
        ref_strand.reserve(static_cast<std::size_t>(j1 - j0));
        for (std::int64_t k = j0; k < j1; ++k) {
            ref_acceptor.push_back(static_cast<std::int32_t>(acceptor_cut[k] - c0));
            ref_strand.push_back(sj_strand[k]);
        }
        accs_[f].set_junctions(ref_offsets, ref_acceptor, ref_strand);
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
