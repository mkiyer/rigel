// testing/bam.cpp — rigel._bam_impl.testing: the read-outs the parity gates take from an accumulator and the
// accumulator set the multimapper offer is held to, and the scan's splice labels.
#include <cstdint>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include "scan/accumulator_python.h"
#include "scan/labels.h"
#include "testing/bind_testing.h"

namespace nb = nanobind;

namespace rigel {
namespace {

using accumulator::Accumulator;
using accumulator::AccumulatorSet;

template <class T>
using InArray = nb::ndarray<const T, nb::ndim<1>, nb::c_contig>;

/// The production Accumulator with what only the gates read: its counters, its census, its held bank in canonical
/// order, its region lookup, and the per-worker merge.
struct InspectedAccumulator : Accumulator {
    using Accumulator::Accumulator;
};

void bind_inspected_accumulator(nb::module_& m) {
    nb::class_<InspectedAccumulator, Accumulator>(m, "Accumulator")
        .def("__init__",
             [](InspectedAccumulator* self, InArray<int64_t> region_bounds, InArray<uint8_t> region_types,
                int max_length, int32_t ref) {
                 new (self) InspectedAccumulator(
                     std::vector<int64_t>(region_bounds.data(), region_bounds.data() + region_bounds.shape(0)),
                     std::vector<uint8_t>(region_types.data(), region_types.data() + region_types.shape(0)),
                     max_length, ref);
             },
             nb::arg("region_bounds"), nb::arg("region_types"), nb::arg("max_length"), nb::arg("ref"))
        .def_prop_ro("qc", [](const InspectedAccumulator& a) { return qc_dict(a.counters()); })
        // Keyed by the specification's `GapResolution` values, so the gate compares dicts directly.
        .def_prop_ro("gap_resolution", [](const InspectedAccumulator& a) { return gap_census_dict(a.gap_census()); })
        // Copied, not viewed, because `canonicalise` may replace the vectors on the next call and a view would dangle.
        .def_prop_ro("deferred", [](InspectedAccumulator& a) { return deferred_dict(a.deferred_canonical()); })
        .def("region_of_pos", [](const InspectedAccumulator& a, int64_t pos) { return a.region_of_pos(pos); },
             nb::arg("pos"))
        // The parallel scan's per-worker merge, so the determinism contract can be tested directly: shard one corpus K
        // ways, merge, and compare with the unsharded run.
        .def("merge_from",
             [](InspectedAccumulator& a, const InspectedAccumulator& other) { a.merge_from(other); },
             nb::arg("other"));
}

/// The accumulator set, built from the flat partition as the scanner builds its own: `at(ref)` hands out each
/// reference's accumulator, and `deferred` / `qc` / `gap_resolution` are the union the payload exports.
void bind_accumulator_set(nb::module_& m) {
    nb::class_<AccumulatorSet>(m, "AccumulatorSet")
        .def("__init__",
             [](AccumulatorSet* self, InArray<int64_t> region_bounds, InArray<int64_t> ref_region_bound_offsets,
                InArray<uint8_t> region_types, int max_length) {
                 new (self) AccumulatorSet(region_bounds.data(), region_bounds.shape(0),
                                           ref_region_bound_offsets.data(), ref_region_bound_offsets.shape(0) - 1,
                                           region_types.data(), region_types.shape(0), max_length);
             },
             nb::arg("region_bounds"), nb::arg("ref_region_bound_offsets"), nb::arg("region_types"),
             nb::arg("max_length"))
        .def("set_sj",
             [](AccumulatorSet& s, InArray<int64_t> offsets, InArray<int64_t> boundary_right,
                InArray<int8_t> sj_strand, InArray<int64_t> ref_region_bound_offsets) {
                 s.set_sj(offsets.data(), offsets.shape(0), boundary_right.data(), sj_strand.data(),
                          sj_strand.shape(0), ref_region_bound_offsets.data());
             },
             nb::arg("offsets"), nb::arg("boundary_right"), nb::arg("sj_strand"), nb::arg("ref_region_bound_offsets"))
        .def_prop_ro("n_refs", &AccumulatorSet::n_refs)
        .def("at", [](AccumulatorSet& s, int32_t ref) -> Accumulator& { return s.at(ref); }, nb::arg("ref"),
             nb::rv_policy::reference_internal)
        // One fragment with every placement it aligned to; each placement reads attributes off whatever it is handed
        // (`ref`, `start`, `end`, `align_strand`, `sj_strand`, `observed_introns`, `hypotheses`, `nm`, `clip`), so the
        // gate passes the specification's own `Placement` objects. Returns the QC key the offer incremented.
        .def("offer",
             [](AccumulatorSet& s, nb::iterable placements) {
                 static accumulator::OfferScratch offer_scratch;
                 // One pass over `placements`, as over each placement's hypotheses, so a one-shot iterable offers the
                 // same fragment as a list. The `Placement`s are taken after the loop, once `marshalled` has stopped
                 // growing.
                 std::vector<MarshalledFragment> marshalled;
                 std::vector<int32_t> refs;
                 for (nb::handle p : placements) {
                     marshalled.emplace_back(nb::cast<int64_t>(nb::getattr(p, "start")),
                                             nb::cast<int64_t>(nb::getattr(p, "end")),
                                             nb::cast<nb::iterable>(nb::getattr(p, "observed_introns")),
                                             nb::cast<int32_t>(nb::getattr(p, "align_strand")),
                                             nb::cast<int32_t>(nb::getattr(p, "sj_strand")),
                                             nb::getattr(p, "hypotheses"),
                                             nb::cast<int32_t>(nb::getattr(p, "nm")),
                                             nb::cast<int32_t>(nb::getattr(p, "clip")));
                     refs.push_back(nb::cast<int32_t>(nb::getattr(p, "ref")));
                 }
                 std::vector<accumulator::Placement> pls;
                 pls.reserve(marshalled.size());
                 for (size_t i = 0; i < marshalled.size(); ++i) pls.push_back({refs[i], marshalled[i].offered});
                 return std::string(accumulator::outcome_key(s.offer(pls.data(), pls.size(), offer_scratch)));
             },
             nb::arg("placements"))
        .def_prop_ro("deferred", [](AccumulatorSet& s) { return deferred_dict(s.total_deferred()); })
        .def_prop_ro("qc", [](const AccumulatorSet& s) { return qc_dict(s.total_counters()); })
        .def_prop_ro("gap_resolution", [](const AccumulatorSet& s) { return gap_census_dict(s.total_gap_census()); });
}

}  // namespace

void bind_testing(nb::module_ m) {
    bind_inspected_accumulator(m);
    bind_accumulator_set(m);
    // The label a `SpliceType` code takes in the census keys and the ZS tag
    m.def("splice_type_label", [](int32_t code) { return std::string(splice_type_label(code)); }, nb::arg("code"));
}

}  // namespace rigel
