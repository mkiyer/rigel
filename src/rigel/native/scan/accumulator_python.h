// accumulator_python.h — the accumulator's values to and from Python, shared by the scan's payload and the bindings:
// a fragment marshalled from Python objects, and the three non-array banks as dicts keyed as the specification's
// `Tally` is.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <nanobind/nanobind.h>

#include "calibration/accumulator_set.h"
#include "ndarray_util.h"

namespace rigel {

/// One `OfferedFragment` marshalled out of Python, with the storage its spans point into. It reads attributes off
/// whatever it is handed, so the parity gate passes the specification's own `GapHypothesis` objects. Not copyable,
/// since `offered` points into this object's vectors; a move keeps it valid.
struct MarshalledFragment {
    std::vector<IntronBlock> observed;
    std::vector<IntronBlock> implied;
    std::vector<accumulator::GapHypothesis> spans;
    accumulator::OfferedFragment offered{};

    MarshalledFragment(int64_t start, int64_t end, nanobind::iterable introns, int32_t align_strand,
                       int32_t sj_strand, nanobind::object hypotheses, int32_t nm = 0, int32_t clip = 0) {
        namespace nb = nanobind;
        const auto intron_of = [](nb::handle item) {
            auto pair = nb::cast<nb::tuple>(item);
            return IntronBlock{0, nb::cast<int32_t>(pair[0]), nb::cast<int32_t>(pair[1]), 0};
        };
        for (nb::handle item : introns) observed.push_back(intron_of(item));
        // One pass over `hypotheses`, which may be a one-shot iterable (a generator). While `implied` grows, each span
        // records its offset into it, turned into a pointer only after the last hypothesis is read, so a reallocation
        // during the fill invalidates nothing.
        std::vector<size_t> first_intron;
        for (nb::handle h : hypotheses) {
            first_intron.push_back(implied.size());
            for (nb::handle item : nb::getattr(h, "introns")) implied.push_back(intron_of(item));
            const int32_t strand = nb::cast<int32_t>(nb::getattr(h, "sj_strand"));
            spans.push_back({nullptr, implied.size() - first_intron.back(), strand});
        }
        for (size_t s = 0; s < spans.size(); ++s) spans[s].introns = implied.data() + first_intron[s];

        offered.start = start;
        offered.end = end;
        offered.observed_introns = observed.data();
        offered.n_observed_introns = observed.size();
        offered.align_strand = align_strand;
        offered.sj_strand = sj_strand;
        offered.hypotheses = spans.data();
        offered.n_hypotheses = spans.size();
        offered.nm = nm;
        offered.clip = clip;
    }

    MarshalledFragment(const MarshalledFragment&) = delete;
    MarshalledFragment& operator=(const MarshalledFragment&) = delete;
    MarshalledFragment(MarshalledFragment&&) = default;
    MarshalledFragment& operator=(MarshalledFragment&&) = default;
};

/// The scratch every bound accumulator call on this thread shares; the scan's workers keep their own.
inline accumulator::DepositScratch& binding_scratch() {
    static thread_local accumulator::DepositScratch scratch;
    return scratch;
}

/// The deposit counters: the QC denominators, naming what each conservation statement downstream excludes; none is
/// derivable after the fact.
inline nanobind::dict qc_dict(const accumulator::DepositCounters& c) {
    nanobind::dict out;
    out["deposited"]                    = c.deposited;
    out["dropped_too_long"]             = c.dropped_too_long;
    out["dropped_empty"]                = c.dropped_empty;
    out["dropped_strand_undefined"]     = c.dropped_strand_undefined;
    out["deferred_undetermined_gap"]    = c.deferred_undetermined_gap;
    out["deferred_multiple_placements"] = c.deferred_multiple_placements;
    return out;
}

/// The umbrella census, keyed by the specification's `GapResolution` values. There is no `gap_resolved_unspliced` key:
/// the unspliced hypothesis is always the longest, so it is never the sole survivor (`GapCensus`).
inline nanobind::dict gap_census_dict(const accumulator::GapCensus& census) {
    nanobind::dict out;
    out["gap_resolved_spliced"]       = census.resolved_spliced;
    out["gap_deferred_rna_or_gdna"]   = census.deferred_rna_or_gdna;
    out["gap_deferred_which_introns"] = census.deferred_which_introns;
    out["gap_deferred_both"]          = census.deferred_both;
    return out;
}

/// The deferred queue as the CSR `Tally.deferred_arrays()` specifies: same keys, same dtype (int64 throughout), same
/// canonical order. Every array is copied out, since a later canonicalise may reseat the vectors under a numpy view.
inline nanobind::dict deferred_dict(const accumulator::DeferredFragments& deferred) {
    nanobind::dict out;
    const auto put = [&out](const char* name, const std::vector<int64_t>& v) {
        out[name] = vec_to_ndarray(std::vector<int64_t>(v));
    };
    put("placement_offsets", deferred.placement_offsets);
    put("ref", deferred.ref);
    put("start", deferred.start);
    put("end", deferred.end);
    put("align_strand", deferred.align_strand);
    put("sj_strand", deferred.sj_strand);
    put("observed_intron_offsets", deferred.observed_intron_offsets);
    put("observed_introns", deferred.observed_introns);
    put("hypothesis_offsets", deferred.hypothesis_offsets);
    put("hypothesis_sj_strand", deferred.hypothesis_sj_strand);
    put("hypothesis_intron_offsets", deferred.hypothesis_intron_offsets);
    put("hypothesis_introns", deferred.hypothesis_introns);
    put("nm", deferred.nm);
    put("clip", deferred.clip);
    return out;
}

}  // namespace rigel
