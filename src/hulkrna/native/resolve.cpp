/**
 * resolve.cpp — nanobind C++ extension for fragment resolution.
 *
 * Ports the entire resolve_fragment() logic from resolution.py to C++,
 * including cgranges overlap queries, set merging, chimera detection,
 * and fragment-length computation.  Calls cgranges directly from C++
 * without round-tripping to Python.
 *
 * The core logic (ResolveContext, ResolvedResult, NativeAccumulator,
 * constants, helper types) now lives in shared headers so that
 * bam_scanner.cpp can call _resolve_core() directly:
 *   - native/constants.h
 *   - native/resolve_context.h
 *
 * Module: hulkrna._resolve_impl
 *
 * Build:
 *   Part of the hulkrna scikit-build-core build — see CMakeLists.txt.
 */

#include "resolve_context.h"

using namespace hulk;

// ================================================================
// nanobind module definition
// ================================================================

NB_MODULE(_resolve_impl, m) {
    m.doc() = "C++ fragment resolution kernel for hulkrna (nanobind).\n\n"
              "Ports resolve_fragment() to C++ with direct cgranges queries.";

    // --- ResolvedResult ---
    nb::class_<ResolvedResult>(m, "ResolvedResult")
        .def_rw("num_hits", &ResolvedResult::num_hits)
        .def_rw("nm", &ResolvedResult::nm)
        .def_ro("splice_type", &ResolvedResult::splice_type)
        .def_ro("exon_strand", &ResolvedResult::exon_strand)
        .def_ro("sj_strand", &ResolvedResult::sj_strand)
        .def_ro("ambig_strand", &ResolvedResult::ambig_strand)
        .def_ro("chimera_type", &ResolvedResult::chimera_type)
        .def_ro("chimera_gap", &ResolvedResult::chimera_gap)
        .def_ro("merge_criteria", &ResolvedResult::merge_criteria)
        .def_ro("read_length", &ResolvedResult::read_length)
        .def_ro("genomic_footprint", &ResolvedResult::genomic_footprint)
        .def_ro("genomic_start", &ResolvedResult::genomic_start)
        .def_prop_ro("is_chimeric", &ResolvedResult::get_is_chimeric)
        .def_prop_ro("is_same_strand", &ResolvedResult::get_is_same_strand)
        .def_prop_ro("is_strand_qualified",
                     &ResolvedResult::get_is_strand_qualified)
        .def_prop_ro("first_t_ind", &ResolvedResult::get_first_t_ind)
        .def_prop_ro("has_frag_lengths",
                     &ResolvedResult::get_has_frag_lengths)
        .def_prop_ro("unique_frag_length",
                     &ResolvedResult::get_unique_frag_length)
        .def_prop_ro("t_inds", &ResolvedResult::get_t_inds)
        .def_prop_ro("frag_lengths", &ResolvedResult::get_frag_lengths)
        .def_prop_ro("overlap_bp", &ResolvedResult::get_overlap_bp)
        ;

    // --- NativeAccumulator ---
    nb::class_<NativeAccumulator>(m, "NativeAccumulator")
        .def(nb::init<>(), "Create an empty native accumulator.")
        .def("append", &NativeAccumulator::append,
             nb::arg("result"), nb::arg("frag_id"),
             "Append a ResolvedResult to the accumulator.")
        .def_prop_ro("size", &NativeAccumulator::get_size,
                     "Number of fragments in the accumulator.")
        .def("finalize", &NativeAccumulator::finalize,
             nb::arg("t_strand_arr"),
             "Finalize to a dict of raw bytes for numpy conversion.")
        ;

    // --- ResolveContext ---
    nb::class_<ResolveContext>(m, "ResolveContext")
        .def(nb::init<>(), "Create an empty resolve context.")
        .def("build_overlap_index", &ResolveContext::build_overlap_index,
             nb::arg("refs"), nb::arg("starts"), nb::arg("ends"),
             nb::arg("iv_types"), nb::arg("tset_data"),
             nb::arg("tset_offsets"),
             "Build the main overlap interval index from collapsed data.")
        .def("build_sj_map", &ResolveContext::build_sj_map,
             nb::arg("refs"), nb::arg("starts"), nb::arg("ends"),
             nb::arg("strands"), nb::arg("tset_data"),
             nb::arg("tset_offsets"),
             "Build the splice-junction exact-match lookup map.")
        .def("build_sj_gap_index", &ResolveContext::build_sj_gap_index,
             nb::arg("refs"), nb::arg("starts"), nb::arg("ends"),
             nb::arg("t_indices"), nb::arg("strands"),
             "Build the SJ gap cgranges index for fragment-length computation.")
        .def("set_metadata", &ResolveContext::set_metadata,
             nb::arg("t_to_g"), nb::arg("n_transcripts"),
             "Set transcript-to-gene mapping and allocate scratch buffers.")
        .def("get_ref_to_id", &ResolveContext::get_ref_to_id,
             "Return the ref-name → integer-ID mapping as a Python dict.")
        .def("resolve", &ResolveContext::resolve,
             nb::arg("exon_ref_ids"), nb::arg("exon_starts"),
             nb::arg("exon_ends"), nb::arg("exon_strands"),
             nb::arg("intron_ref_ids"), nb::arg("intron_starts"),
             nb::arg("intron_ends"), nb::arg("intron_strands"),
             nb::arg("genomic_footprint"),
             "Resolve a fragment to its compatible transcript set.\n\n"
             "Returns a 13-element tuple or None for intergenic fragments.")
        .def("resolve_fragment", &ResolveContext::resolve_fragment,
             nb::arg("frag"),
             "Resolve a Fragment object directly.\n\n"
             "Returns a ResolvedResult or None for intergenic fragments.\n"
             "Eliminates Python marshaling overhead of the resolve() path.")
        .def("set_gene_strands", &ResolveContext::set_gene_strands,
             nb::arg("g_to_strand"),
             "Set gene strand mapping for BAM scanner model training.")
        .def("set_transcript_strands", &ResolveContext::set_transcript_strands,
             nb::arg("t_strand"),
             "Set per-transcript strand array (direct lookup, no gene indirection).")
        ;

    // --- Expose C++ enum constants as module-level attributes ---
    // These mirror the Python IntEnum values in hulkrna.types / hulkrna.splice.
    // Single authoritative source: constants.h

    // Strand
    m.attr("STRAND_NONE")      = hulk::STRAND_NONE;
    m.attr("STRAND_POS")       = hulk::STRAND_POS;
    m.attr("STRAND_NEG")       = hulk::STRAND_NEG;
    m.attr("STRAND_AMBIGUOUS") = hulk::STRAND_AMBIGUOUS;

    // SpliceType
    m.attr("SPLICE_UNSPLICED")       = hulk::SPLICE_UNSPLICED;
    m.attr("SPLICE_SPLICED_UNANNOT") = hulk::SPLICE_SPLICED_UNANNOT;
    m.attr("SPLICE_SPLICED_ANNOT")   = hulk::SPLICE_SPLICED_ANNOT;

    // MergeCriteria
    m.attr("MC_INTERSECTION")          = hulk::MC_INTERSECTION;
    m.attr("MC_INTERSECTION_NONEMPTY") = hulk::MC_INTERSECTION_NONEMPTY;
    m.attr("MC_UNION")                 = hulk::MC_UNION;
    m.attr("MC_EMPTY")                 = hulk::MC_EMPTY;

    // ChimeraType
    m.attr("CHIMERA_NONE")           = hulk::CHIMERA_NONE;
    m.attr("CHIMERA_TRANS")          = hulk::CHIMERA_TRANS;
    m.attr("CHIMERA_CIS_STRAND_SAME") = hulk::CHIMERA_CIS_STRAND_SAME;
    m.attr("CHIMERA_CIS_STRAND_DIFF") = hulk::CHIMERA_CIS_STRAND_DIFF;

    // IntervalType
    m.attr("ITYPE_EXON")           = hulk::ITYPE_EXON;
    m.attr("ITYPE_INTRON")         = static_cast<int8_t>(1);  // TRANSCRIPT slot
    m.attr("ITYPE_UNAMBIG_INTRON") = hulk::ITYPE_UNAMBIG_INTRON;
}
