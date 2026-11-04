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
             "Seed the ref-name -> ID space in canonical index.ref_names order.\n\n"
             "Call it before build_overlap_index so the resolver's ref ids\n"
             "match index.ref_name_to_id (the space used by t_df, region_df and the\n"
             "calibration accumulator partition).")
        .def("build_overlap_index", &FragmentResolver::build_overlap_index, nb::arg("refs"), nb::arg("starts"),
             nb::arg("ends"), nb::arg("iv_types"), nb::arg("tset_data"), nb::arg("tset_offsets"),
             "Build the main overlap interval index from collapsed data.")
        .def("build_sj_map", &FragmentResolver::build_sj_map, nb::arg("refs"), nb::arg("starts"), nb::arg("ends"),
             nb::arg("strands"), nb::arg("tset_data"), nb::arg("tset_offsets"),
             "Build the splice-junction exact-match lookup map.")
        .def("build_sj_blacklist_map", &FragmentResolver::build_sj_blacklist_map, nb::arg("refs"),
             nb::arg("starts"), nb::arg("ends"), nb::arg("max_anchor_left"), nb::arg("max_anchor_right"),
             "Build the splice-junction artifact blacklist map.\n\n"
             "Entries are keyed by (ref, start, end) without strand.\n"
             "A CIGAR sj is rejected when either its left or right\n"
             "anchor is <= the blacklist maximum for that sj.")
        .def("set_metadata", &FragmentResolver::set_metadata, nb::arg("n_transcripts"),
             "Set the transcript count, which sizes each scan worker's scratch.")
        .def("set_transcript_strands", &FragmentResolver::set_transcript_strands, nb::arg("t_strand"),
             "Set per-transcript strand array (direct lookup, no gene indirection).")
        .def("set_nrna_parent_index", &FragmentResolver::set_nrna_parent_index, nb::arg("nrna_parent"),
             "Set each transcript's synthetic nascent span (int32; -1 = none).\n"
             "Nascent-span candidates are derived from annotated-transcript hits\n"
             "during resolution; synthetics are not in cgranges.")
        .def("build_exon_index", &FragmentResolver::build_exon_index, nb::arg("offsets"), nb::arg("starts"),
             nb::arg("ends"),
             "Build the per-transcript exon CSR: each candidate's introns, for its path and its length.")
        .def("set_max_fragment_length", &FragmentResolver::set_max_fragment_length,
             nb::arg("max_fragment_length"),
             "The library's fragment-length limit: no pair whose molecule is longer, along the genome and "
             "along every witness transcript, is one.")
        .def("set_intron_overhang_tolerance", &FragmentResolver::set_intron_overhang_tolerance, nb::arg("K"),
             "Set the intron-overhang tolerance K (bp): how many bases an alignment\n"
             "may run into an annotated intron, and an annotated intron overhang a\n"
             "mate gap, and still be read as spliced there. In the gap-hypothesis\n"
             "enumeration an annotated intron in an unsequenced gap may overlap up\n"
             "to K aligned bases on either side and is still cut whole, and a\n"
             "transcript is contradicted only by an intron overlapping more than K\n"
             "bases of an aligned block. Negative K is clamped to 0; the resolver's\n"
             "own default, before this is called, is 0.");

    rigel::bind_testing(m.def_submodule("testing", "The resolver's test seams; production never imports them."));
}
