// resolve_module.cpp — the rigel._resolve_impl module: the FragmentResolver the index builds and the scan reads.
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "resolve/fragment_resolver.h"
#include "testing/bind_testing.h"

namespace nb = nanobind;
using rigel::FragmentResolver;

NB_MODULE(_resolve_impl, m) {
    m.doc() = "C++ fragment resolution kernel for rigel (nanobind).";

    nb::class_<FragmentResolver>(m, "FragmentResolver")
        .def(nb::init<>(), "Create an empty resolve context.")
        .def("set_ref_names", &FragmentResolver::set_ref_names, nb::arg("ref_names"),
             "Install the ref-name -> ID space in canonical index.ref_names order.\n\n"
             "Call it first, so the resolver's ref ids match index.ref_name_to_id\n"
             "(the space used by t_df, region_df and the calibration accumulator\n"
             "partition). The builders refuse a reference it did not install.")
        .def("build_overlap_index", &FragmentResolver::build_overlap_index, nb::arg("refs"), nb::arg("starts"),
             nb::arg("ends"), nb::arg("iv_types"), nb::arg("tset_data"), nb::arg("tset_offsets"),
             "Build the main overlap interval index from collapsed data.\n\n"
             "Raises ValueError for an interval on a reference set_ref_names did not install.")
        .def("build_sj_map", &FragmentResolver::build_sj_map, nb::arg("refs"), nb::arg("starts"), nb::arg("ends"),
             nb::arg("strands"), nb::arg("tset_data"), nb::arg("tset_offsets"),
             "Build the splice-junction exact-match lookup map.\n\n"
             "Raises ValueError for a junction on a reference set_ref_names did not install.")
        .def("set_metadata", &FragmentResolver::set_metadata, nb::arg("n_transcripts"),
             "Set the transcript count, which sizes each scan worker's scratch.")
        .def("set_transcript_strands", &FragmentResolver::set_transcript_strands, nb::arg("t_strand"),
             "Set per-transcript strand array (direct lookup, no gene indirection).")
        .def("set_nascent_parent_index", &FragmentResolver::set_nascent_parent_index, nb::arg("nascent_parent"),
             "Set each transcript's synthetic nascent span (int32; -1 = none).\n"
             "Nascent-span candidates are derived from annotated-transcript hits\n"
             "during resolution; synthetics are not in cgranges.")
        .def("build_exon_index", &FragmentResolver::build_exon_index, nb::arg("offsets"), nb::arg("starts"),
             nb::arg("ends"),
             "Build the per-transcript exon CSR: each candidate's introns, for its path and its length.");

    rigel::bind_testing(m.def_submodule("testing", "The resolver's test seams; production never imports them."));
}
