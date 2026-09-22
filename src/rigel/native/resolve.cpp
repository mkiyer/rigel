/**
 * resolve.cpp — nanobind C++ extension for fragment resolution.
 *
 * Fragment resolution: cgranges overlap queries, set merging, chimera
 * detection and fragment-length computation, calling cgranges directly
 * without round-tripping to Python.
 *
 * The core logic (FragmentResolver, ResolvedFragment, FragmentAccumulator,
 * constants, helper types) lives in shared headers so that
 * bam_scanner.cpp can call _resolve_core() directly:
 *   - native/constants.h
 *   - native/resolve_context.h
 *
 * Module: rigel._resolve_impl
 *
 * Build:
 *   Part of the rigel scikit-build-core build — see CMakeLists.txt.
 */

#include "resolve_context.h"

using namespace rigel;

// ================================================================
// nanobind module definition
// ================================================================

NB_MODULE(_resolve_impl, m) {
    m.doc() = "C++ fragment resolution kernel for rigel (nanobind).\n\n"
              "Ports resolve_fragment() to C++ with direct cgranges queries.";

    // --- ResolvedFragment ---
    nb::class_<ResolvedFragment>(m, "ResolvedFragment")
        .def_rw("num_hits", &ResolvedFragment::num_hits)
        .def_rw("nm", &ResolvedFragment::nm)
        .def_ro("splice_type", &ResolvedFragment::splice_type)
        .def_ro("align_strand", &ResolvedFragment::align_strand)
        .def_ro("sj_strand", &ResolvedFragment::sj_strand)
        .def_ro("ambig_strand", &ResolvedFragment::ambig_strand)
        .def_ro("chimera_type", &ResolvedFragment::chimera_type)
        .def_ro("merge_criteria", &ResolvedFragment::merge_criteria)
        .def_ro("read_length", &ResolvedFragment::read_length)
        .def_ro("genomic_footprint", &ResolvedFragment::genomic_footprint)
        .def_ro("genomic_start", &ResolvedFragment::genomic_start)
        .def_ro("exon_bp_pos", &ResolvedFragment::exon_bp_pos)
        .def_ro("exon_bp_neg", &ResolvedFragment::exon_bp_neg)
        .def_ro("tx_bp_pos", &ResolvedFragment::tx_bp_pos)
        .def_ro("tx_bp_neg", &ResolvedFragment::tx_bp_neg)
        .def_prop_ro("is_same_strand", &ResolvedFragment::get_is_same_strand)
        .def_prop_ro("is_strand_qualified",
                     &ResolvedFragment::get_is_strand_qualified)
        .def_prop_ro("t_inds", &ResolvedFragment::get_t_inds)
        .def_prop_ro("frag_lengths", &ResolvedFragment::get_frag_lengths)
        .def_prop_ro("overlap_bp", &ResolvedFragment::get_overlap_bp)
        .def_ro("n_gap_hypotheses", &ResolvedFragment::n_gap_hypotheses)
        .def_prop_ro("gap_hypotheses",
                     &ResolvedFragment::get_gap_hypotheses)
        ;

    // --- FragmentAccumulator ---
    nb::class_<FragmentAccumulator>(m, "FragmentAccumulator")
        .def(nb::init<>(), "Create an empty native accumulator.")
        .def("append", &FragmentAccumulator::append,
             nb::arg("result"), nb::arg("frag_id"),
             "Append a ResolvedFragment to the accumulator.")
        .def_prop_ro("size", &FragmentAccumulator::get_size,
                     "Number of fragments in the accumulator.")
        .def("finalize", &FragmentAccumulator::finalize,
             "Finalize to a dict of numpy arrays, leaving the accumulator empty.")
        ;

    // --- FragmentResolver ---
    nb::class_<FragmentResolver>(m, "FragmentResolver")
        .def(nb::init<>(), "Create an empty resolve context.")
        .def("set_ref_names", &FragmentResolver::set_ref_names,
             nb::arg("ref_names"),
             "Seed the ref-name -> ID space in canonical index.ref_names order.\n\n"
             "Must be called BEFORE build_overlap_index so the resolver's ref ids\n"
             "match index.ref_name_to_id (the space used by t_df, region_df and the\n"
             "calibration accumulator partition).")
        .def("build_overlap_index", &FragmentResolver::build_overlap_index,
             nb::arg("refs"), nb::arg("starts"), nb::arg("ends"),
             nb::arg("iv_types"), nb::arg("tset_data"),
             nb::arg("tset_offsets"),
             "Build the main overlap interval index from collapsed data.")
        .def("build_sj_map", &FragmentResolver::build_sj_map,
             nb::arg("refs"), nb::arg("starts"), nb::arg("ends"),
             nb::arg("strands"), nb::arg("tset_data"),
             nb::arg("tset_offsets"),
             "Build the splice-junction exact-match lookup map.")
        .def("build_sj_blacklist_map",
             &FragmentResolver::build_sj_blacklist_map,
             nb::arg("refs"), nb::arg("starts"), nb::arg("ends"),
             nb::arg("max_anchor_left"), nb::arg("max_anchor_right"),
             "Build the splice-junction artifact blacklist map.\n\n"
             "SpliceJunctions are keyed by (ref, start, end) without strand.\n"
             "A CIGAR sj is rejected when EITHER its left or right\n"
             "anchor is <= the blacklist maximum for that sj.")
        .def("set_metadata", &FragmentResolver::set_metadata,
             nb::arg("n_transcripts"),
             "Set the transcript count and allocate scratch buffers.")
        .def("resolve_fragment", &FragmentResolver::resolve_fragment,
             nb::arg("frag"),
             "Resolve a Fragment object to its compatible transcript set.\n\n"
             "Returns a ResolvedFragment or None for intergenic fragments.")
        .def("set_transcript_strands", &FragmentResolver::set_transcript_strands,
             nb::arg("t_strand"),
             "Set per-transcript strand array (direct lookup, no gene indirection).")
        .def("set_nrna_parent_index",
             &FragmentResolver::set_nrna_parent_index,
             nb::arg("nrna_parent"),
             "Set per-transcript nRNA parent index (int32; -1 = none).\n"
             "Synthetic nRNA candidates are derived from real-tx hits\n"
             "during _resolve_core; synthetics are not in cgranges.")
        .def("build_exon_index", &FragmentResolver::build_exon_index,
             nb::arg("offsets"), nb::arg("starts"), nb::arg("ends"),
             nb::arg("cumsum"),
             "Build per-transcript exon CSR index for FL computation.")
        .def("set_max_fragment_length",
             &FragmentResolver::set_max_fragment_length,
             nb::arg("max_fragment_length"),
             "The library's fragment-length limit. A same-strand, same-reference pair whose implied "
             "fragment length is within it is an ordinary genomic molecule, not a chimera.")
        .def("set_splicing_anchor_tolerance",
             &FragmentResolver::set_splicing_anchor_tolerance,
             nb::arg("K"),
             "Set splicing-anchor tolerance K (bp, >= 0) used by the\n"
             "SPLICED_IMPLICIT per-intron whole-containment discriminant.")
        ;

    // --- The SpliceType constants, mirrored for the parity gate against rigel.splice ---
    m.attr("SPLICE_UNSPLICED")       = rigel::SPLICE_UNSPLICED;
    m.attr("SPLICE_SPLICED_UNANNOT") = rigel::SPLICE_SPLICED_UNANNOT;
    m.attr("SPLICE_SPLICED_ANNOT")   = rigel::SPLICE_SPLICED_ANNOT;
    m.attr("SPLICE_IMPLICIT")        = rigel::SPLICE_IMPLICIT;
    m.attr("SPLICE_ARTIFACT")        = rigel::SPLICE_ARTIFACT;
}
