// scan_result.cpp — what the scan returns: its stats, the strand observations, and the calibration payload.
#include <algorithm>
#include <string>

#include "ndarray_util.h"
#include "scan/accumulator_python.h"
#include "scan/bam_scanner.h"
#include "scan/labels.h"

namespace nb = nanobind;

namespace rigel {

using accumulator::Accumulator;
using accumulator::AccumulatorSet;
using accumulator::kNFragmentPools;
using accumulator::kNStrandColumns;

nb::dict BamScanStats::to_dict() const {
    nb::dict out;
    out["total"] = total;
    out["qc_fail"] = qc_fail;
    out["unmapped"] = unmapped;
    out["secondary"] = secondary;
    out["supplementary"] = supplementary;
    out["duplicate"] = duplicate;
    out["proper_pair"] = proper_pair;
    pairs.export_to(out);
    out["n_fragments"] = n_fragments;
    out["n_intergenic"] = n_intergenic;
    out["n_with_exon"] = n_with_exon;
    out["n_with_annotated_sj"] = n_with_annotated_sj;
    out["n_same_strand"] = n_same_strand;
    out["n_ambig_strand"] = n_ambig_strand;
    out["n_strand_trained"] = n_strand_trained;
    out["n_strand_skipped_no_sj"] = n_strand_skipped_no_sj;
    out["n_strand_skipped_ambig_strand"] = n_strand_skipped_ambig_strand;
    out["n_strand_skipped_ambiguous"] = n_strand_skipped_ambiguous;
    out["n_multimapper_groups"] = n_multimapper_groups;
    out["n_multimapper_alignments"] = n_multimapper_alignments;
    out["n_sj_observed"] = records.n_sj_observed;
    out["n_sj_blacklisted"] = records.n_sj_blacklisted;
    out["n_sj_tag_tagged"] = records.sj_tags.tagged;
    out["n_sj_tag_untagged"] = records.sj_tags.untagged;
    out["n_sj_tag_conflicting"] = records.sj_tags.conflicting;
    out["n_sj_tag_invalid"] = records.sj_tags.invalid;
    out["n_nm_tag_missing"] = records.nm_tag.missing;
    out["n_nm_tag_invalid"] = records.nm_tag.invalid;
    out["n_nm_conflicting"] = n_nm_conflicting;
    out["n_unique_unspliced_mismatches"] = n_unique_unspliced_mismatches;
    out["n_unique_unspliced_bases"] = n_unique_unspliced_bases;
    out["n_nh_tag_missing"] = records.nh_tag.missing;
    out["n_nh_tag_invalid"] = records.nh_tag.invalid;
    for (int32_t v = HIT_USABLE + 1; v < N_HIT_VERDICTS; v++)
        out[("n_hits_rejected_" + std::string(hit_verdict_label(v))).c_str()] = hits_rejected[v];
    out["n_nh_above_reported"] = n_nh_above_reported;
    out["n_read_pairs_lacking_mate"] = n_read_pairs_lacking_mate;
    for (size_t i = 0; i < splice_census.size(); i++)
        out[("n_census_" + std::string(splice_type_label(static_cast<int32_t>(i)))).c_str()] = splice_census[i];
    return out;
}

namespace {

/// The strand observations: the all-exonic diagnostic's 2×2 cell by cell, and the per-junction strand table as six
/// parallel arrays sorted by (ref, start, end, motif strand), since the table is an unordered map merged across
/// workers and no downstream number may depend on scheduling or hash order.
nb::dict strand_dict(const StrandObservations& obs) {
    nb::dict out;
    out["exonic_pos_pos"] = obs.exonic.pos_pos;
    out["exonic_pos_neg"] = obs.exonic.pos_neg;
    out["exonic_neg_pos"] = obs.exonic.neg_pos;
    out["exonic_neg_neg"] = obs.exonic.neg_neg;

    std::vector<SJKey> keys;
    keys.reserve(obs.sj_strand_table.size());
    for (const auto& [key, counts] : obs.sj_strand_table) keys.push_back(key);
    std::sort(keys.begin(), keys.end());
    const size_t n = keys.size();
    std::vector<int32_t> ref(n);
    std::vector<int64_t> start(n), end(n);
    std::vector<int8_t> motif(n);
    std::vector<uint64_t> sense(n), antisense(n);
    for (size_t i = 0; i < n; i++) {
        const SJStrandCounts& c = obs.sj_strand_table.at(keys[i]);
        ref[i] = keys[i].ref_id;
        start[i] = keys[i].start;  // int32 upstream; widened once, here
        end[i] = keys[i].end;
        motif[i] = static_cast<int8_t>(keys[i].strand);
        sense[i] = c.n_sense;
        antisense[i] = c.n_antisense;
    }
    out["sj_ref_id"] = vec_to_ndarray(std::move(ref));
    out["sj_start"] = vec_to_ndarray(std::move(start));
    out["sj_end"] = vec_to_ndarray(std::move(end));
    out["sj_motif_strand"] = vec_to_ndarray(std::move(motif));
    out["sj_n_sense"] = vec_to_ndarray(std::move(sense));
    out["sj_n_antisense"] = vec_to_ndarray(std::move(antisense));
    return out;
}

/// Each reference's first object on one flat axis, and the total: n_refs + 1 cumulative offsets.
template <class Count>
std::vector<int64_t> ref_offsets(const AccumulatorSet& set, Count count) {
    std::vector<int64_t> offsets(set.n_refs() + 1, 0);
    for (size_t f = 0; f < set.n_refs(); ++f)
        offsets[f + 1] = offsets[f] + static_cast<int64_t>(count(set.at(static_cast<int32_t>(f))));
    return offsets;
}

/// One bank of every reference, concatenated in reference order on its flat axis: `rows(a)` rows of ``cols`` values
/// from each accumulator, `value(a, i, c)` the value at row i, column c.
template <class T, class Rows, class Value>
nb::object concatenated(const AccumulatorSet& set, size_t cols, Rows rows, Value value) {
    std::vector<T> out;
    for (size_t f = 0; f < set.n_refs(); ++f) {
        const Accumulator& a = set.at(static_cast<int32_t>(f));
        for (size_t i = 0; i < rows(a); ++i)
            for (size_t c = 0; c < cols; ++c) out.push_back(value(a, i, c));
    }
    return vec_to_ndarray(std::move(out));
}

/// A library-wide histogram: every reference's ``n`` bins summed. A size mismatch throws, since skipping it would drop
/// that reference's bins silently; the deposited-length anchor's sum must equal `deposited`.
template <class T, class Bins, class Size>
nb::object summed(const AccumulatorSet& set, size_t n, Bins bins, Size size, const char* what) {
    std::vector<T> total(n, 0);
    for (size_t f = 0; f < set.n_refs(); ++f) {
        const Accumulator& a = set.at(static_cast<int32_t>(f));
        if (size(a) != n) {
            throw std::runtime_error("build_result: reference " + std::to_string(f) + " has " +
                                     std::to_string(size(a)) + " " + what + " bins but the payload expects " +
                                     std::to_string(n));
        }
        for (size_t i = 0; i < n; ++i) total[i] += bins(a)[i];
    }
    return vec_to_ndarray(std::move(total));
}

/// Every bank of the tally under the specification's `Tally` field names, each a fresh owning array, emitted flat (the
/// two-column ones are reshaped in Python), so Python's buffers are independent of the accumulators' struct layout.
void put_banks(nb::dict& cal, const AccumulatorSet& set, size_t pool_row) {
    using A = Accumulator;
    constexpr size_t kCols = kNStrandColumns;
    const auto regions = [](const A& a) { return a.n_regions(); };
    const auto boundaries = [](const A& a) { return a.n_boundaries(); };
    const auto sj = [](const A& a) { return a.n_sj(); };

    cal["region_contained_count"] = concatenated<uint32_t>(set, kCols, regions, [](const A& a, size_t i, size_t c) {
        return a.regions_data()[i].contained_count[c];
    });
    cal["region_start_count"] = concatenated<uint32_t>(set, kCols, regions, [](const A& a, size_t i, size_t c) {
        return a.region_start_count_data()[i * kCols + c];
    });
    cal["boundary_unspliced_count"] = concatenated<uint32_t>(set, kCols, boundaries, [](const A& a, size_t i, size_t c) {
        return a.boundaries_data()[i].unspliced_count[c];
    });
    // One value per boundary: the length moment and the conserved masses carry no strand
    cal["boundary_unspliced_inv_length_sum"] =
        concatenated<double>(set, 1, boundaries, [](const A& a, size_t i, size_t) {
            return a.boundaries_data()[i].unspliced_inv_length_sum;
        });
    cal["boundary_spliced_count"] = concatenated<uint32_t>(set, kCols, boundaries, [](const A& a, size_t i, size_t c) {
        return a.boundaries_data()[i].spliced_count[c];
    });
    cal["boundary_unspliced_mass"] = concatenated<double>(set, 1, boundaries, [](const A& a, size_t i, size_t) {
        return a.boundaries_data()[i].unspliced_mass;
    });
    cal["boundary_spliced_mass"] = concatenated<double>(set, 1, boundaries, [](const A& a, size_t i, size_t) {
        return a.boundaries_data()[i].spliced_mass;
    });
    cal["sj_count"] = concatenated<uint32_t>(set, kCols, sj, [](const A& a, size_t i, size_t c) {
        return a.sj_data()[i].count[c];
    });
    cal["sj_inv_length_sum"] = concatenated<double>(set, 1, sj, [](const A& a, size_t i, size_t) {
        return a.sj_data()[i].inv_length_sum;
    });
    // Two columns, unlike every other mass: the junction mass carries a strand (`SpliceJunction::mass`)
    cal["sj_mass"] = concatenated<double>(set, kCols, sj, [](const A& a, size_t i, size_t c) {
        return a.sj_data()[i].mass[c];
    });
    // The pools and the unconditional length histogram are library-wide, so summed over references, not concatenated
    cal["pool_lengths"] = summed<int64_t>(
        set, kNFragmentPools * pool_row, [](const A& a) { return a.pool_lengths_data(); },
        [](const A& a) { return a.pool_lengths_size(); }, "pool");
    cal["deposited_lengths"] = summed<uint32_t>(
        set, pool_row, [](const A& a) { return a.deposited_lengths_data(); },
        [](const A& a) { return a.deposited_lengths_size(); }, "deposited-length");
}

}  // namespace

nb::dict BamScanner::build_result() {
    AccumulatorSet& set = *acc_set_;
    nb::dict cal;
    // Echo the partition back, so a consumer can locate every object without reloading the index.
    cal["region_bounds"] = vec_to_ndarray(std::vector<int64_t>(partition_.region_bounds));
    cal["ref_region_bound_offsets"] = vec_to_ndarray(std::vector<int64_t>(partition_.ref_region_bound_offsets));
    cal["ref_region_offsets"] = vec_to_ndarray(ref_offsets(set, [](const Accumulator& a) { return a.n_regions(); }));
    cal["ref_boundary_offsets"] =
        vec_to_ndarray(ref_offsets(set, [](const Accumulator& a) { return a.n_boundaries(); }));
    cal["ref_sj_offsets"] = vec_to_ndarray(ref_offsets(set, [](const Accumulator& a) { return a.n_sj(); }));
    put_banks(cal, set, static_cast<size_t>(partition_.max_length) + 1);
    // The QC denominators, the side buffer and the umbrella census go through the exporters the bindings use, so the
    // payload and the bound accumulators share their keys.
    cal["qc"] = qc_dict(set.total_counters());
    cal["deferred"] = deferred_dict(set.total_deferred());
    cal["gap_resolution"] = gap_census_dict(set.total_gap_census());
    cal["n_strand_columns"] = static_cast<int>(kNStrandColumns);
    cal["max_length"] = partition_.max_length;
    cal["n_refs"] = static_cast<int>(set.n_refs());
    // The library's mismatch rate travels with the tally, so a cached scan drains at its own rate.
    cal["unique_unspliced_mismatches"] = stats_.n_unique_unspliced_mismatches;
    cal["unique_unspliced_bases"] = stats_.n_unique_unspliced_bases;

    nb::dict result;
    result["stats"] = stats_.to_dict();
    result["strand_observations"] = strand_dict(strand_obs_);
    result["calibration"] = cal;
    return result;
}

}  // namespace rigel
