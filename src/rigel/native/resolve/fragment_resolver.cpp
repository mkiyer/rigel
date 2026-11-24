// fragment_resolver.cpp — FragmentResolver's index, its pair rule, its resolution and its gap-hypothesis enumeration.
#include "resolve/fragment_resolver.h"

#include <algorithm>
#include <cstdlib>
#include <stdexcept>

#include "constants.h"
#include "path_length.h"
#include "resolve/sorted_sets.h"

namespace rigel {

// ================================================================
// The index
// ================================================================

void FragmentResolver::set_ref_names(const std::vector<std::string>& ref_names) {
    ref_to_id_.clear();
    id_to_ref_.clear();
    ref_to_id_.reserve(ref_names.size());
    id_to_ref_.reserve(ref_names.size());
    for (const auto& name : ref_names)
        if (ref_to_id_.emplace(name, static_cast<int32_t>(id_to_ref_.size())).second) id_to_ref_.push_back(name);
}

int32_t FragmentResolver::installed_ref_id(const std::string& ref, const char* builder) const {
    const auto it = ref_to_id_.find(ref);
    if (it == ref_to_id_.end())
        throw std::invalid_argument(std::string(builder) + ": reference '" + ref +
                                    "' is not one set_ref_names installed");
    return it->second;
}

void FragmentResolver::build_overlap_index(const std::vector<std::string>& refs, const std::vector<int32_t>& starts,
                                           const std::vector<int32_t>& ends, const std::vector<int32_t>& iv_types,
                                           const std::vector<int32_t>& tset_data,
                                           const std::vector<int32_t>& tset_offsets) {
    const size_t n = refs.size();
    for (const auto& ref : refs) installed_ref_id(ref, "build_overlap_index");

    cr_ = cr_init();
    for (size_t i = 0; i < n; i++) cr_add(cr_, refs[i].c_str(), starts[i], ends[i], static_cast<int32_t>(i));
    cr_index(cr_);

    iv_type_.resize(n);
    for (size_t i = 0; i < n; i++) iv_type_[i] = static_cast<int8_t>(iv_types[i]);
    t_set_data_ = tset_data;
    t_set_offsets_ = tset_offsets;
}

void FragmentResolver::build_sj_map(const std::vector<std::string>& refs, const std::vector<int32_t>& starts,
                                    const std::vector<int32_t>& ends, const std::vector<int32_t>& strands,
                                    const std::vector<int32_t>& tset_data, const std::vector<int32_t>& tset_offsets) {
    sj_map_data_ = tset_data;
    for (size_t i = 0; i < refs.size(); i++) {
        const SJKey key{installed_ref_id(refs[i], "build_sj_map"), starts[i], ends[i], strands[i]};
        sj_map_[key] = {tset_offsets[i], tset_offsets[i + 1] - tset_offsets[i]};
    }
}

void FragmentResolver::build_exon_index(const std::vector<int32_t>& offsets, const std::vector<int32_t>& starts,
                                        const std::vector<int32_t>& ends) {
    exon_offsets_ = offsets;
    exon_starts_ = starts;
    exon_ends_ = ends;
}

void FragmentResolver::set_max_fragment_length(int32_t L) {
    if (L < 0) throw std::invalid_argument("set_max_fragment_length: must be >= 0, got " + std::to_string(L));
    max_fragment_length_ = L;
}

int32_t FragmentResolver::ref_id(const char* name) const {
    const auto it = ref_to_id_.find(name);
    return it == ref_to_id_.end() ? -1 : it->second;
}

int32_t FragmentResolver::transcript_strand(int32_t t) const {
    return (t >= 0 && t < static_cast<int32_t>(t_strand_arr_.size())) ? t_strand_arr_[t] : STRAND_NONE;
}

// ================================================================
// A transcript's introns
// ================================================================

FragmentResolver::TranscriptIntrons FragmentResolver::introns_of(int32_t t) const {
    if (!in_exon_index(t)) return {nullptr, nullptr, 0};
    const int32_t begin = exon_offsets_[t];
    const int32_t n = exon_offsets_[t + 1] - begin - 1;
    if (n <= 0) return {nullptr, nullptr, 0};
    return {exon_ends_.data() + begin, exon_starts_.data() + begin + 1, n};
}

// The candidate set needs this test because `merge_sets` falls back to a union, which can hold a contradicted
// transcript; the pair rule's witness and the gap enumeration ask it too. Abutting is fine (an intron may start where
// a block ends), and K is the slack the gap match uses: an aligner routinely runs a read a base or two past a splice
// site.
bool FragmentResolver::transcript_contradicts_blocks(int32_t t,
                                                     const std::vector<std::pair<int32_t, int32_t>>& blocks,
                                                     int32_t K) const {
    const TranscriptIntrons in = introns_of(t);
    for (const auto& [start, end] : blocks) {
        // The first intron that could end inside this block, then on while one still starts before its end.
        int32_t i = static_cast<int32_t>(std::upper_bound(in.ends, in.ends + in.n, start) - in.ends);
        for (; i < in.n && in.starts[i] < end; ++i)
            if (std::min(in.ends[i], end) - std::max(in.starts[i], start) > K) return true;
    }
    return false;
}

// A mate gap can span two annotated introns, and two transcripts differing only in the second must stay apart, so
// every intron in the gap is kept, in order.
void FragmentResolver::collect_transcript_introns_in_gap(int32_t t, const GapBlock& gap, int32_t K, int32_t strand,
                                                         std::vector<IntronBlock>& out) const {
    const TranscriptIntrons in = introns_of(t);
    if (in.n == 0 || gap.end <= gap.start) return;
    const int64_t max_end = static_cast<int64_t>(gap.end) + K;
    int32_t i = static_cast<int32_t>(
        std::lower_bound(in.starts, in.starts + in.n, static_cast<int64_t>(gap.start) - K) - in.starts);
    const int32_t by_end = static_cast<int32_t>(std::upper_bound(in.ends, in.ends + in.n, gap.start) - in.ends);
    for (i = std::max(i, by_end); i < in.n; ++i) {
        if (in.starts[i] >= gap.end) break;                      // no overlap with the gap
        if (static_cast<int64_t>(in.ends[i]) > max_end) break;  // protrudes past the slack
        if (in.ends[i] > in.starts[i]) out.push_back({gap.ref_id, in.starts[i], in.ends[i], strand});
    }
}

// Overlap, not containment: a mate that runs past an annotated end still overlaps, and the EM prices its extra bases as
// overhang.
bool FragmentResolver::overlaps_both(int32_t t, int32_t f_lo, int32_t f_hi, int32_t r_lo, int32_t r_hi) const {
    if (!in_exon_index(t)) return false;
    const int32_t first = exon_offsets_[t], last = exon_offsets_[t + 1] - 1;
    if (last < first) return false;
    const int32_t t_lo = exon_starts_[first], t_hi = exon_ends_[last];
    return t_lo < f_hi && t_hi > f_lo && t_lo < r_hi && t_hi > r_lo;
}

void FragmentResolver::retain_overlapping_both(std::vector<int32_t>& ts, int32_t f_lo, int32_t f_hi, int32_t r_lo,
                                               int32_t r_hi) const {
    ts.erase(std::remove_if(ts.begin(), ts.end(),
                            [&](int32_t t) { return !overlaps_both(t, f_lo, f_hi, r_lo, r_hi); }),
             ts.end());
}

// The synthetic nascent spans are not in the overlap index, and have no introns to offer.
void FragmentResolver::transcripts_overlapping(int32_t ref_id, int32_t start, int32_t end, ResolverScratch& scratch,
                                               std::vector<int32_t>& out) const {
    out.clear();
    if (!cr_ || ref_id < 0 || ref_id >= static_cast<int32_t>(id_to_ref_.size())) return;
    const int64_t n = cr_overlap(cr_, id_to_ref_[ref_id].c_str(), start, end, &scratch.buf.ptr, &scratch.buf.cap);
    for (int64_t k = 0; k < n; ++k) {
        const int32_t label = cr_label(cr_, scratch.buf.ptr[k]);
        if (iv_type_[label] != ITYPE_TRANSCRIPT) continue;
        for (int32_t i = t_set_offsets_[label]; i < t_set_offsets_[label + 1]; ++i) out.push_back(t_set_data_[i]);
    }
    sort_unique(out);
}

// ================================================================
// The annotated junctions
// ================================================================

void FragmentResolver::sj_lookup_into(int32_t ref_id, int32_t start, int32_t end, int32_t strand,
                                      std::vector<int32_t>& out) const {
    out.clear();
    const auto append = [&](const std::pair<int32_t, int32_t>& slice) {
        const auto [off, cnt] = slice;
        out.insert(out.end(), sj_map_data_.begin() + off, sj_map_data_.begin() + off + cnt);
    };
    const auto it = sj_map_.find(SJKey{ref_id, start, end, strand});
    if (it != sj_map_.end()) {
        append(it->second);
    } else if (strand != STRAND_POS && strand != STRAND_NEG) {
        // No definite strand: the junction on either strand
        const auto pos = sj_map_.find(SJKey{ref_id, start, end, STRAND_POS});
        const auto neg = sj_map_.find(SJKey{ref_id, start, end, STRAND_NEG});
        if (pos != sj_map_.end()) append(pos->second);
        if (neg != sj_map_.end()) append(neg->second);
    }
    sort_unique(out);
}

/// The transcripts that hold every junction in ``junctions`` (a hit's, sorted and unique, at least one), each matched
/// as `sj_lookup_into` matches it, into ``holders``: `HIT_USABLE` when there are any. Otherwise the reason there are
/// none, the first that applies: `HIT_CONTRADICTORY` when no one molecule can carry the junctions as observed (two
/// overlap without being equal, or their motifs name both strands), `HIT_NOVEL_JUNCTION` when one is held by no
/// transcript, and `HIT_UNSUPPORTED` when each is held by some transcript and none holds them all.
int32_t FragmentResolver::junction_holders(const std::vector<IntronBlock>& junctions, ResolverScratch& scratch,
                                           std::vector<int32_t>& holders) const {
    holders.clear();
    int32_t motifs = STRAND_NONE;
    for (size_t j = 0; j < junctions.size(); ++j) {
        const IntronBlock& x = junctions[j];
        motifs |= x.strand;
        if (j == 0) continue;
        // Sorted, so two that overlap unequal leave such a pair adjacent.
        const IntronBlock& w = junctions[j - 1];
        if (x.ref_id == w.ref_id && x.start < w.end && (x.start != w.start || x.end != w.end))
            return HIT_CONTRADICTORY;
    }
    if (motifs == STRAND_AMBIGUOUS) return HIT_CONTRADICTORY;
    for (size_t j = 0; j < junctions.size(); ++j) {
        const IntronBlock& x = junctions[j];
        sj_lookup_into(x.ref_id, x.start, x.end, x.strand, scratch.pair_tmp);
        if (scratch.pair_tmp.empty()) {
            holders.clear();
            return HIT_NOVEL_JUNCTION;
        }
        if (j == 0) {
            holders.swap(scratch.pair_tmp);
        } else {
            holders.erase(std::set_intersection(holders.begin(), holders.end(), scratch.pair_tmp.begin(),
                                                scratch.pair_tmp.end(), holders.begin()),
                          holders.end());
        }
    }
    return holders.empty() ? HIT_UNSUPPORTED : HIT_USABLE;
}

int32_t FragmentResolver::pair_holders(const std::vector<IntronBlock>& junctions, int32_t f_lo, int32_t f_hi,
                                       int32_t r_lo, int32_t r_hi, ResolverScratch& scratch,
                                       std::vector<int32_t>& holders) const {
    const int32_t why = junction_holders(junctions, scratch, holders);
    if (why != HIT_USABLE) return why;
    retain_overlapping_both(holders, f_lo, f_hi, r_lo, r_hi);
    return holders.empty() ? HIT_UNSUPPORTED : HIT_USABLE;
}

// ================================================================
// The pair rule
// ================================================================

namespace {

/// The length of [lo, hi) less the union of the ``skips`` inside it.
int64_t length_less_skips(int32_t lo, int32_t hi, const std::vector<std::pair<int32_t, int32_t>>& skips) {
    int64_t removed = 0;
    int32_t run_s = 0, run_e = 0;
    bool open = false;
    for (auto [s, e] : skips) {
        s = std::max(s, lo);
        e = std::min(e, hi);
        if (e <= s) continue;
        if (open && s <= run_e) {
            run_e = std::max(run_e, e);
        } else {
            if (open) removed += run_e - run_s;
            run_s = s;
            run_e = e;
            open = true;
        }
    }
    if (open) removed += run_e - run_s;
    return static_cast<int64_t>(hi) - lo - removed;
}

/// The stretches of [lo, hi) on ``ref_id`` that no interval of ``covered`` (sorted) covers.
void uncovered_stretches(int32_t ref_id, int32_t lo, int32_t hi, const std::vector<std::pair<int32_t, int32_t>>& covered,
                         std::vector<GapBlock>& out) {
    out.clear();
    int32_t cur = lo;
    for (const auto& [s, e] : covered) {
        if (s > cur) out.push_back({ref_id, cur, std::min(s, hi)});
        cur = std::max(cur, e);
    }
    if (cur < hi) out.push_back({ref_id, cur, hi});
}

/// The total length of the distinct introns in ``cut``, which it sorts.
int64_t distinct_intron_length(std::vector<IntronBlock>& cut) {
    std::sort(cut.begin(), cut.end(), [](const IntronBlock& x, const IntronBlock& y) {
        return x.start != y.start ? x.start < y.start : x.end < y.end;
    });
    int64_t total = 0;
    for (size_t i = 0; i < cut.size(); ++i)
        if (i == 0 || cut[i].start != cut[i - 1].start || cut[i].end != cut[i - 1].end)
            total += cut[i].end - cut[i].start;
    return total;
}

}  // namespace

int32_t FragmentResolver::validate_pair(int32_t ref_id, const PairMate& f, const PairMate& r,
                                        const std::vector<IntronBlock>& junctions, ResolverScratch& scratch) const {
    const int32_t K = intron_overhang_tolerance_;
    const int64_t max_len = max_fragment_length_;
    const auto in_range = [max_len](int64_t length) { return length >= 1 && length <= max_len; };
    const int32_t a = f.blocks->front().first - f.clip_lead;
    const int32_t b = r.blocks->back().second + r.clip_trail;
    if (a >= b) return HIT_OUTWARD;
    const int32_t lo = std::min(a, r.blocks->front().first);
    const int32_t hi = std::max(b, f.blocks->back().second);

    auto& skips = scratch.pair_skips;  // each mate's holes between consecutive blocks, the two merged
    skips.clear();
    for (const PairMate* m : {&f, &r})
        for (size_t i = 1; i < m->blocks->size(); ++i)
            skips.push_back({(*m->blocks)[i - 1].second, (*m->blocks)[i].first});
    std::sort(skips.begin(), skips.end());
    skips.erase(std::unique(skips.begin(), skips.end()), skips.end());
    const int64_t length = length_less_skips(lo, hi, skips);
    if (in_range(length)) return HIT_USABLE;

    // The witnesses: the transcripts over the forward mate that also overlap the reverse one, holding every observed
    // junction when there is one.
    const int32_t f_lo = f.blocks->front().first, f_hi = f.blocks->back().second;
    const int32_t r_lo = r.blocks->front().first, r_hi = r.blocks->back().second;
    auto& cands = scratch.pair_candidates;
    if (junctions.empty()) {
        transcripts_overlapping(ref_id, f_lo, f_hi, scratch, cands);
        retain_overlapping_both(cands, f_lo, f_hi, r_lo, r_hi);
    } else if (pair_holders(junctions, f_lo, f_hi, r_lo, r_hi, scratch, cands) != HIT_USABLE) {
        return HIT_LENGTH;
    }

    auto& aligned = scratch.pair_aligned;  // both mates' aligned blocks
    aligned.assign(f.blocks->begin(), f.blocks->end());
    aligned.insert(aligned.end(), r.blocks->begin(), r.blocks->end());
    std::sort(aligned.begin(), aligned.end());
    auto& covered = scratch.pair_covered;  // aligned or skipped
    covered.assign(aligned.begin(), aligned.end());
    covered.insert(covered.end(), skips.begin(), skips.end());
    std::sort(covered.begin(), covered.end());
    auto& stretches = scratch.pair_stretches;  // what neither mate aligned or skipped
    uncovered_stretches(ref_id, lo, hi, covered, stretches);

    for (int32_t t : cands) {
        if (!in_exon_index(t) || transcript_contradicts_blocks(t, aligned, K)) continue;
        auto& cut = scratch.pair_cut;
        cut.clear();
        for (const GapBlock& g : stretches) collect_transcript_introns_in_gap(t, g, K, STRAND_NONE, cut);
        if (in_range(length - distinct_intron_length(cut))) return HIT_USABLE;
    }
    return HIT_LENGTH;
}

// ================================================================
// Resolution
// ================================================================

void FragmentResolver::resolve(const std::vector<ExonBlock>& exons, const std::vector<IntronBlock>& introns,
                               int32_t genomic_footprint, ResolvedFragment& out, ResolverScratch& scratch) const {
    if (exons.empty()) throw std::invalid_argument("resolve: a fragment needs at least one aligned block");
    scratch.reset_per_fragment();
    out.genomic_footprint = genomic_footprint;

    query_blocks(exons, out, scratch);
    add_nascent_spans(exons, scratch);
    read_junctions(introns, out, scratch);
    choose_candidates(!introns.empty(), out, scratch);
    candidate_exon_bp(out, scratch);
    out.ambig_strand = candidates_on_both_strands(out.t_inds) ? 1 : 0;

    // Every fragment, whatever its splice type: a spliced fragment can also carry an intron in its mate gap. This only
    // enumerates; whether the fragment deposits or waits for the second pass is decided where the outcome is reported.
    enumerate_gap_hypotheses(exons, introns, out.splice_type == SPLICE_SPLICED_ANNOT, out, scratch);
    // SPLICE_SPLICED_IMPLICIT promotes unspliced fragments only. `splice_type` is the scanner's census of what it saw
    // and feeds scoring, the buffer, strand training and `rigel report`, so relabelling an observed SPLICED_ANNOT
    // fragment would move mass between reported categories. Deferral is a separate axis.
    if (out.splice_type == SPLICE_UNSPLICED && out.n_gap_hypotheses() > 1) out.splice_type = SPLICE_SPLICED_IMPLICIT;

    out.genomic_start = fragment_extent(exons).first;
}

/// Each block's exon-interval and transcript-interval transcripts, sorted and unique, and every transcript's aligned
/// bases on its exons; the fragment's aligned strand and length.
void FragmentResolver::query_blocks(const std::vector<ExonBlock>& exons, ResolvedFragment& out,
                                    ResolverScratch& scratch) const {
    scratch.exon_t_sets.resize(exons.size());
    scratch.transcript_t_sets.resize(exons.size());
    out.align_strand = STRAND_NONE;
    out.read_length = 0;
    for (size_t bi = 0; bi < exons.size(); bi++) {
        const ExonBlock& eb = exons[bi];
        out.align_strand |= eb.strand;
        out.read_length += eb.end - eb.start;
        if (eb.ref_id < 0 || eb.ref_id >= static_cast<int32_t>(id_to_ref_.size())) continue;

        auto& block_exon_t = scratch.exon_t_sets[bi];
        auto& block_transcript_t = scratch.transcript_t_sets[bi];
        const int64_t n = cr_overlap(cr_, id_to_ref_[eb.ref_id].c_str(), eb.start, eb.end, &scratch.buf.ptr,
                                     &scratch.buf.cap);
        for (int64_t k = 0; k < n; k++) {
            const int64_t idx = scratch.buf.ptr[k];
            const int32_t label = cr_label(cr_, idx);
            const int32_t* ts = t_set_data_.data() + t_set_offsets_[label];
            const int32_t* ts_end = t_set_data_.data() + t_set_offsets_[label + 1];
            if (iv_type_[label] == ITYPE_EXON) {
                const int32_t bp =
                    std::max(std::min(eb.end, cr_end(cr_, idx)) - std::max(eb.start, cr_start(cr_, idx)), 0);
                for (const int32_t* t = ts; t < ts_end; ++t) {
                    block_exon_t.push_back(*t);
                    if (bp > 0) scratch.add_exon_bp(*t, bp);
                }
            } else if (iv_type_[label] == ITYPE_TRANSCRIPT) {
                block_transcript_t.insert(block_transcript_t.end(), ts, ts_end);
            }
        }
        sort_unique(block_exon_t);
        sort_unique(block_transcript_t);
    }
}

/// The synthetic nascent spans of every annotated transcript the blocks overlap: each span the fragment has bases in
/// joins every block's exon set (so the merge's intersection surfaces it across blocks), with those bases as exonic.
void FragmentResolver::add_nascent_spans(const std::vector<ExonBlock>& exons, ResolverScratch& scratch) const {
    if (nrna_parent_.empty()) return;
    const auto add_parent = [&](int32_t t) {
        if (t < 0 || t >= static_cast<int32_t>(nrna_parent_.size())) return;
        if (nrna_parent_[t] >= 0) scratch.nrna_t.push_back(nrna_parent_[t]);
    };
    for (const auto& s : scratch.exon_t_sets)
        for (int32_t t : s) add_parent(t);
    for (const auto& s : scratch.transcript_t_sets)
        for (int32_t t : s) add_parent(t);
    sort_unique(scratch.nrna_t);

    for (int32_t span : scratch.nrna_t) {
        // A synthetic span has one exon by construction.
        const int32_t ns = exon_starts_[exon_offsets_[span]];
        const int32_t ne = exon_ends_[exon_offsets_[span]];
        int32_t bp = 0;
        for (const auto& eb : exons) {
            const int32_t lo = std::max(eb.start, ns), hi = std::min(eb.end, ne);
            if (hi > lo) bp += hi - lo;
        }
        if (bp <= 0) continue;
        scratch.add_exon_bp(span, bp);
        for (auto& block_set : scratch.exon_t_sets) {
            const auto it = std::lower_bound(block_set.begin(), block_set.end(), span);
            if (it == block_set.end() || *it != span) block_set.insert(it, span);
        }
    }
}

/// The observed junctions' holders, the fragment's motif strand and its strand-table key. The scan's hit verdict
/// guarantees one annotated transcript holds every junction, so a spliced fragment without holders is a broken
/// precondition, refused.
void FragmentResolver::read_junctions(const std::vector<IntronBlock>& introns, ResolvedFragment& out,
                                      ResolverScratch& scratch) const {
    scratch.sj_holders.clear();
    out.sj_strand = STRAND_NONE;
    out.sj_key_ref = out.sj_key_start = out.sj_key_end = -1;
    if (introns.empty()) return;
    if (junction_holders(introns, scratch, scratch.sj_holders) != HIT_USABLE)
        throw std::invalid_argument(
            "resolve: a spliced fragment needs one annotated transcript holding every junction it observes (the "
            "scan's hit verdict withholds any other)");
    for (const auto& x : introns) out.sj_strand |= x.strand;
    out.sj_key_ref = introns[0].ref_id;
    out.sj_key_start = introns[0].start;
    out.sj_key_end = introns[0].end;
}

/// The candidates and the splice type. A spliced fragment's candidates are the junctions' holders its blocks overlap,
/// or the holders themselves; an unspliced one's the merged exon sets united with the merged transcript sets. With no
/// overlap at all the fragment is intergenic: no candidate, and it still flows through the buffer, so calibration can
/// count it.
void FragmentResolver::choose_candidates(bool spliced, ResolvedFragment& out, ResolverScratch& scratch) const {
    if (spliced) {
        const std::vector<int32_t> by_blocks = merge_sets(scratch.exon_t_sets);
        std::vector<int32_t> isect =
            by_blocks.empty() ? std::vector<int32_t>{} : sorted_intersection(by_blocks, scratch.sj_holders);
        out.t_inds = isect.empty() ? scratch.sj_holders : std::move(isect);
        out.splice_type = SPLICE_SPLICED_ANNOT;
    } else {
        out.t_inds = sorted_union(merge_sets(scratch.exon_t_sets), merge_sets(scratch.transcript_t_sets));
        out.splice_type = SPLICE_UNSPLICED;
    }
}

/// Each candidate's aligned bases on its exons: its tally when a block overlaps it, the whole read when it overlaps
/// nothing the blocks did, and 0 when the blocks overlap nothing at all.
void FragmentResolver::candidate_exon_bp(ResolvedFragment& out, ResolverScratch& scratch) const {
    auto& all_overlap_t = scratch.all_overlap_t;
    for (const auto& s : scratch.exon_t_sets) all_overlap_t.insert(all_overlap_t.end(), s.begin(), s.end());
    for (const auto& s : scratch.transcript_t_sets) all_overlap_t.insert(all_overlap_t.end(), s.begin(), s.end());
    sort_unique(all_overlap_t);

    const bool any_overlap = !all_overlap_t.empty() && out.read_length > 0;
    out.exon_bp.clear();
    out.exon_bp.reserve(out.t_inds.size());
    for (int32_t t : out.t_inds) {
        if (!any_overlap)
            out.exon_bp.push_back(0);
        else if (std::binary_search(all_overlap_t.begin(), all_overlap_t.end(), t))
            out.exon_bp.push_back(scratch.t_exon_bp[t]);
        else
            out.exon_bp.push_back(out.read_length);
    }
}

bool FragmentResolver::candidates_on_both_strands(const std::vector<int32_t>& candidates) const {
    bool seen = false;
    int32_t first = 0;
    for (int32_t t : candidates) {
        if (t < 0 || t >= static_cast<int32_t>(t_strand_arr_.size())) continue;
        if (!seen) {
            first = t_strand_arr_[t];
            seen = true;
        } else if (t_strand_arr_[t] != first) {
            return true;
        }
    }
    return false;
}

// ================================================================
// The gap hypotheses
// ================================================================

namespace {

/// Record ``path`` as a hypothesis, or fold its strand into the identical one already there.
void add_hypothesis(ResolvedFragment& out, const std::vector<IntronBlock>& path, int32_t strand) {
    for (int32_t h = 0; h < out.n_gap_hypotheses(); ++h) {
        const int32_t lo = out.gap_intron_offsets[h], hi = out.gap_intron_offsets[h + 1];
        if (hi - lo != static_cast<int32_t>(path.size())) continue;
        bool same = true;
        for (int32_t i = 0; i < hi - lo && same; ++i)
            same = out.gap_introns[lo + i].start == path[i].start && out.gap_introns[lo + i].end == path[i].end;
        if (!same) continue;
        // One path, but its supporters may disagree about the strand. The hypothesis carries one `sj_strand`, and
        // keeping the first supporter's would depend on GTF order, so disagreement makes it STRAND_AMBIGUOUS
        // (contradictory evidence, as for the fragment-level `sj_strand`), which `deposit` does not credit; once
        // ambiguous it stays so. A junction annotated on both strands is biologically impossible and the index warns
        // about it.
        if (out.gap_sj_strand[h] != strand) out.gap_sj_strand[h] = STRAND_AMBIGUOUS;
        return;
    }
    out.gap_introns.insert(out.gap_introns.end(), path.begin(), path.end());
    out.gap_intron_offsets.push_back(static_cast<int32_t>(out.gap_introns.size()));
    out.gap_sj_strand.push_back(strand);
}

/// The empty path: cut nothing. Idempotent, because two routes can both call for it.
void emit_unspliced_hypothesis(ResolvedFragment& out) {
    for (int32_t h = 0; h < out.n_gap_hypotheses(); ++h)
        if (out.gap_intron_offsets[h] == out.gap_intron_offsets[h + 1]) return;
    out.gap_intron_offsets.push_back(static_cast<int32_t>(out.gap_introns.size()));
    out.gap_sj_strand.push_back(STRAND_NONE);
}

}  // namespace

/// The holes between consecutive blocks on one reference, less every one an observed junction explains exactly. The
/// walk also yields each observed CIGAR-N intron (`50M200N50M` is two blocks), and searching it again would let a nearby
/// intron within ±K answer for it and shorten L. Those holes are dropped by exact (start, end) equality, never by
/// overlap, which would also drop a genuine gap intron abutting an observed one.
void FragmentResolver::unsequenced_gaps(const std::vector<ExonBlock>& exons,
                                        const std::vector<IntronBlock>& observed_introns,
                                        ResolverScratch& scratch) const {
    auto& gaps = scratch.implicit_gaps;
    gaps.clear();
    if (exons.size() < 2 || exon_offsets_.empty()) return;
    const std::vector<ExonBlock>* ordered = &exons;
    if (!std::is_sorted(exons.begin(), exons.end(), exon_block_less)) {
        scratch.implicit_blocks.assign(exons.begin(), exons.end());
        std::sort(scratch.implicit_blocks.begin(), scratch.implicit_blocks.end(), exon_block_less);
        ordered = &scratch.implicit_blocks;
    }
    int32_t cur_ref = (*ordered)[0].ref_id;
    int32_t cur_end = (*ordered)[0].end;
    for (size_t k = 1; k < ordered->size(); ++k) {
        const ExonBlock& block = (*ordered)[k];
        if (block.ref_id != cur_ref) {
            cur_ref = block.ref_id;
            cur_end = block.end;
        } else if (block.start > cur_end) {
            gaps.push_back({cur_ref, cur_end, block.start});
            cur_end = block.end;
        } else {
            cur_end = std::max(cur_end, block.end);
        }
    }
    const auto observed = [&observed_introns](const GapBlock& gap) {
        for (const auto& x : observed_introns)
            if (x.ref_id == gap.ref_id && x.start == gap.start && x.end == gap.end) return true;
        return false;
    };
    gaps.erase(std::remove_if(gaps.begin(), gaps.end(), observed), gaps.end());
}

/// Every explanation of this fragment's unsequenced gaps, identical paths grouped into one hypothesis; the accumulator
/// arbitrates between them, this only enumerates. A hypothesis is a path through the annotation: a gap may hold several
/// introns and whole exons no read touched. Each compatible candidate implies one path, its own introns inside the
/// gaps, and paths are compared whole across every gap, so no hypothesis mixes one transcript's intron in one gap with
/// another's in the next.
///
/// Each candidate's length comes from the same path: ``out.frag_lengths[i]`` is `path_length` of the fragment's extent
/// less the observed junctions and candidate i's gap introns (-1 when the path has no base), the length the
/// accumulator deposits. A candidate the reads contradict (`transcript_contradicts_blocks`) joins no hypothesis but
/// keeps a length: its molecule with the covered intron retained, which the EM scores beside the covered bases'
/// overhang penalty.
void FragmentResolver::enumerate_gap_hypotheses(const std::vector<ExonBlock>& exons,
                                                const std::vector<IntronBlock>& observed_introns,
                                                bool certified_rna, ResolvedFragment& out,
                                                ResolverScratch& scratch) const {
    out.gap_introns.clear();
    out.gap_sj_strand.clear();
    out.gap_intron_offsets.assign(1, 0);

    const std::pair<int32_t, int32_t> span = fragment_extent(exons);
    const auto length_of = [&](const std::vector<IntronBlock>& path) -> int32_t {
        const int64_t L = path_length(observed_introns.data(), observed_introns.size(), path.data(), path.size(),
                                      span.first, span.second, scratch.path_introns, scratch.path_segments);
        return L > 0 ? static_cast<int32_t>(L) : -1;
    };

    unsequenced_gaps(exons, observed_introns, scratch);
    const auto& gaps = scratch.implicit_gaps;
    const std::vector<int32_t>& candidates = out.t_inds;
    auto& path = scratch.hypothesis_path;
    path.clear();
    // No unsequenced gap: one hypothesis, the unspliced one, and one length on every candidate, the extent less the
    // observed junctions.
    if (gaps.empty()) {
        out.frag_lengths.assign(candidates.size(), length_of(path));
        emit_unspliced_hypothesis(out);
        return;
    }

    const int32_t K = intron_overhang_tolerance_;
    auto& blocks = scratch.held_blocks;  // as `transcript_contradicts_blocks` reads them
    blocks.clear();
    for (const auto& e : exons) blocks.push_back({e.start, e.end});

    bool any_candidate_implies_nothing = false;
    out.frag_lengths.assign(candidates.size(), -1);
    for (size_t i = 0; i < candidates.size(); ++i) {
        const int32_t t = candidates[i];
        const int32_t strand = transcript_strand(t);
        path.clear();
        for (const GapBlock& gap : gaps) collect_transcript_introns_in_gap(t, gap, K, strand, path);
        out.frag_lengths[i] = length_of(path);
        if (transcript_contradicts_blocks(t, blocks, K)) continue;
        if (path.empty()) {
            any_candidate_implies_nothing = true;
            continue;
        }
        add_hypothesis(out, path, strand);
    }

    // The unspliced hypothesis: present when some compatible transcript implies nothing in the gaps (a retained-intron
    // isoform), when no other hypothesis survives, and always when no annotated junction was sequenced, since then the
    // molecule may be gDNA or RNA that has not spliced there.
    if (any_candidate_implies_nothing || !certified_rna || out.n_gap_hypotheses() == 0) emit_unspliced_hypothesis(out);
}

}  // namespace rigel
