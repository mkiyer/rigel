// deferred_fragments.cpp — the held-fragment bank: appending, rebuilding, comparing and canonically ordering records.
#include "calibration/deferred_fragments.h"

#include <algorithm>
#include <utility>

namespace rigel::accumulator {

void DeferredFragments::append_placement(const OfferedFragment& fragment, std::int64_t ref_id,
                                         std::int64_t clipped_start, std::int64_t clipped_end) {
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
    }
    hypothesis_offsets.push_back(static_cast<std::int64_t>(hypothesis_sj_strand.size()));
    nm.push_back(fragment.nm);
    clip.push_back(fragment.clip);
}

void DeferredFragments::copy_placement(const DeferredFragments& from, std::size_t i) {
    ref.push_back(from.ref[i]);
    start.push_back(from.start[i]);
    end.push_back(from.end[i]);
    align_strand.push_back(from.align_strand[i]);
    sj_strand.push_back(from.sj_strand[i]);
    for (std::size_t p = static_cast<std::size_t>(from.observed_intron_offsets[i]);
         p < static_cast<std::size_t>(from.observed_intron_offsets[i + 1]); ++p) {
        observed_introns.push_back(from.observed_introns[2 * p]);
        observed_introns.push_back(from.observed_introns[2 * p + 1]);
    }
    observed_intron_offsets.push_back(static_cast<std::int64_t>(observed_introns.size() / 2));
    for (std::size_t g = static_cast<std::size_t>(from.hypothesis_offsets[i]);
         g < static_cast<std::size_t>(from.hypothesis_offsets[i + 1]); ++g) {
        hypothesis_sj_strand.push_back(from.hypothesis_sj_strand[g]);
        for (std::size_t p = static_cast<std::size_t>(from.hypothesis_intron_offsets[g]);
             p < static_cast<std::size_t>(from.hypothesis_intron_offsets[g + 1]); ++p) {
            hypothesis_introns.push_back(from.hypothesis_introns[2 * p]);
            hypothesis_introns.push_back(from.hypothesis_introns[2 * p + 1]);
        }
        hypothesis_intron_offsets.push_back(static_cast<std::int64_t>(hypothesis_introns.size() / 2));
    }
    hypothesis_offsets.push_back(static_cast<std::int64_t>(hypothesis_sj_strand.size()));
    nm.push_back(from.nm[i]);
    clip.push_back(from.clip[i]);
}

const OfferedFragment& DeferredFragments::offered_of(std::size_t p, std::size_t h0, std::size_t h1,
                                                     OfferedStorage& storage) const {
    const HeldArrays bank{start.data(),
                          end.data(),
                          align_strand.data(),
                          sj_strand.data(),
                          observed_intron_offsets.data(),
                          observed_introns.data(),
                          hypothesis_sj_strand.data(),
                          hypothesis_intron_offsets.data(),
                          hypothesis_introns.data(),
                          nm.data(),
                          clip.data()};
    const std::size_t g = static_cast<std::size_t>(hypothesis_offsets[p]);
    return rebuild_offered(bank, p, g + h0, g + h1, storage);
}

namespace {

/// Python's sequence comparison (a prefix sorts first), returning -1, 0 or 1. A flat run of (start, end) pairs
/// compares as Python's tuple of pairs, so one helper serves the observed and the implied intron runs.
int lex_compare(const std::int64_t* a, std::size_t na, const std::int64_t* b, std::size_t nb) noexcept {
    const std::size_t common = std::min(na, nb);
    for (std::size_t i = 0; i < common; ++i)
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    if (na == nb) return 0;
    return na < nb ? -1 : 1;
}

int scalar_compare(std::int64_t a, std::int64_t b) noexcept {
    if (a == b) return 0;
    return a < b ? -1 : 1;
}

/// Run ``i`` of a flat (start, end) array under its pair offsets: its first value and its value count.
std::pair<const std::int64_t*, std::size_t> intron_run(const std::vector<std::int64_t>& offsets,
                                                       const std::vector<std::int64_t>& introns, std::size_t i) {
    const std::size_t lo = static_cast<std::size_t>(offsets[i]);
    const std::size_t hi = static_cast<std::size_t>(offsets[i + 1]);
    return {introns.data() + 2 * lo, 2 * (hi - lo)};
}

}  // namespace

int DeferredFragments::compare_placements(std::size_t a, std::size_t b) const noexcept {
    for (const std::vector<std::int64_t>* column : {&ref, &start, &end, &align_strand, &sj_strand}) {
        const int c = scalar_compare((*column)[a], (*column)[b]);
        if (c != 0) return c;
    }
    const auto [a_obs, a_n_obs] = intron_run(observed_intron_offsets, observed_introns, a);
    const auto [b_obs, b_n_obs] = intron_run(observed_intron_offsets, observed_introns, b);
    const int c_obs = lex_compare(a_obs, a_n_obs, b_obs, b_n_obs);
    if (c_obs != 0) return c_obs;

    const std::size_t a_h0 = static_cast<std::size_t>(hypothesis_offsets[a]);
    const std::size_t b_h0 = static_cast<std::size_t>(hypothesis_offsets[b]);
    const std::size_t common = std::min(n_hypotheses(a), n_hypotheses(b));
    for (std::size_t k = 0; k < common; ++k) {
        const std::size_t ga = a_h0 + k, gb = b_h0 + k;
        const auto [a_in, a_n_in] = intron_run(hypothesis_intron_offsets, hypothesis_introns, ga);
        const auto [b_in, b_n_in] = intron_run(hypothesis_intron_offsets, hypothesis_introns, gb);
        int c = lex_compare(a_in, a_n_in, b_in, b_n_in);
        if (c != 0) return c;
        c = scalar_compare(hypothesis_sj_strand[ga], hypothesis_sj_strand[gb]);
        if (c != 0) return c;
    }
    return scalar_compare(static_cast<std::int64_t>(n_hypotheses(a)), static_cast<std::int64_t>(n_hypotheses(b)));
}

void DeferredFragments::clear() {
    for (std::vector<std::int64_t>* v : value_arrays(*this)) v->clear();
    for (std::vector<std::int64_t>* v : offset_arrays(*this)) {
        v->clear();
        v->push_back(0);
    }
}

void DeferredFragments::merge_from(const DeferredFragments& other) {
    // CSR concatenation: offsets append shifted by this bank's extent, skipping `other`'s leading 0.
    const auto into_offsets = offset_arrays(*this);
    const auto from_offsets = offset_arrays(other);
    for (std::size_t k = 0; k < into_offsets.size(); ++k) {
        std::vector<std::int64_t>& into = *into_offsets[k];
        const std::vector<std::int64_t>& from = *from_offsets[k];
        const std::int64_t base = into.back();
        for (std::size_t i = 1; i < from.size(); ++i) into.push_back(base + from[i]);
    }
    const auto into_values = value_arrays(*this);
    const auto from_values = value_arrays(other);
    for (std::size_t k = 0; k < into_values.size(); ++k)
        into_values[k]->insert(into_values[k]->end(), from_values[k]->begin(), from_values[k]->end());
}

void DeferredFragments::canonicalise() {
    const std::size_t n = n_records();
    if (n < 2) return;

    // The specification's key: the sequence of the placements' `compare_placements` keys, a prefix sorting first.
    const auto record_less = [this](std::size_t a, std::size_t b) {
        const std::size_t a0 = static_cast<std::size_t>(placement_offsets[a]);
        const std::size_t a1 = static_cast<std::size_t>(placement_offsets[a + 1]);
        const std::size_t b0 = static_cast<std::size_t>(placement_offsets[b]);
        const std::size_t b1 = static_cast<std::size_t>(placement_offsets[b + 1]);
        const std::size_t common = std::min(a1 - a0, b1 - b0);
        for (std::size_t k = 0; k < common; ++k) {
            const int c = compare_placements(a0 + k, b0 + k);
            if (c != 0) return c < 0;
        }
        if (a1 - a0 != b1 - b0) return (a1 - a0) < (b1 - b0);
        // (nm, clip) breaks only a content tie, so it never reorders records the content orders (the draw's stream
        // follows that order); it makes the order total, and a set's offer keeps the fewest-mismatch copy.
        for (std::size_t k = 0; k < common; ++k) {
            if (nm[a0 + k] != nm[b0 + k]) return nm[a0 + k] < nm[b0 + k];
            if (clip[a0 + k] != clip[b0 + k]) return clip[a0 + k] < clip[b0 + k];
        }
        return false;
    };

    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    // Not stable, and need not be: records that tie on the whole key are byte-identical.
    std::sort(order.begin(), order.end(), record_less);

    bool identity = true;
    for (std::size_t i = 0; i < n && identity; ++i) identity = (order[i] == i);
    if (identity) return;  // already canonical: a repeat export copies nothing

    DeferredFragments out;
    const auto out_arrays = value_arrays(out);
    const auto in_arrays = value_arrays(*this);
    for (std::size_t k = 0; k < out_arrays.size(); ++k) out_arrays[k]->reserve(in_arrays[k]->size());
    const auto out_offsets = offset_arrays(out);
    const auto in_offsets = offset_arrays(*this);
    for (std::size_t k = 0; k < out_offsets.size(); ++k) out_offsets[k]->reserve(in_offsets[k]->size());
    for (const std::size_t r : order) {
        for (std::size_t p = static_cast<std::size_t>(placement_offsets[r]);
             p < static_cast<std::size_t>(placement_offsets[r + 1]); ++p) {
            out.copy_placement(*this, p);
        }
        out.close_record();
    }
    *this = std::move(out);
}

}  // namespace rigel::accumulator
