// testing/resolve.cpp — rigel._resolve_impl.testing: a Python fragment resolved as the scan resolves a hit, the row it
// yields, the chunk rows gather into, and the splice codes the parity gate holds to rigel.splice.
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include "resolve/fragment_resolver.h"
#include "testing/bind_testing.h"

namespace nb = nanobind;

namespace rigel {
namespace {

/// A fragment's ``(ref, start, end, strand)`` blocks, the reference by name.
template <class Block>
std::vector<Block> read_blocks(const FragmentResolver& resolver, nb::object blocks) {
    std::vector<Block> out;
    const size_t n = nb::len(blocks);
    for (size_t i = 0; i < n; i++) {
        nb::object b = blocks[nb::int_(i)];
        out.push_back({resolver.ref_id(nb::cast<std::string>(b[nb::int_(0)]).c_str()), nb::cast<int32_t>(b[nb::int_(1)]),
                       nb::cast<int32_t>(b[nb::int_(2)]), nb::cast<int32_t>(b[nb::int_(3)])});
    }
    return out;
}

/// A fragment-like object (``exons``, ``introns``, ``genomic_footprint``) resolved as the scan resolves a hit, its
/// junctions sorted as the scan hands them over; None when it has no exon block.
nb::object resolve_fragment(const FragmentResolver& resolver, nb::object frag) {
    std::vector<ExonBlock> exons = read_blocks<ExonBlock>(resolver, frag.attr("exons"));
    if (exons.empty()) return nb::none();
    std::vector<IntronBlock> introns = read_blocks<IntronBlock>(resolver, frag.attr("introns"));
    sort_unique_junctions(introns);
    ResolverScratch scratch(resolver.n_transcripts());
    ResolvedFragment out;
    resolver.resolve(exons, introns, nb::cast<int32_t>(frag.attr("genomic_footprint")), out, scratch);
    return nb::cast(std::move(out));
}

/// Per candidate, ``value[i]`` as ``{t_ind: value}``, the entries ``keep`` admits.
template <class Keep>
nb::dict by_candidate(const ResolvedFragment& r, const std::vector<int32_t>& value, Keep keep) {
    nb::dict d;
    for (size_t i = 0; i < r.t_inds.size(); i++)
        if (keep(value[i])) d[nb::cast(r.t_inds[i])] = nb::cast(value[i]);
    return d;
}

}  // namespace

void bind_testing(nb::module_ m) {
    nb::class_<ResolvedFragment>(m, "ResolvedFragment")
        .def_rw("num_hits", &ResolvedFragment::num_hits)
        .def_rw("nm", &ResolvedFragment::nm)
        .def_ro("splice_type", &ResolvedFragment::splice_type)
        .def_ro("align_strand", &ResolvedFragment::align_strand)
        .def_ro("ambig_strand", &ResolvedFragment::ambig_strand)
        .def_ro("read_length", &ResolvedFragment::read_length)
        .def_ro("genomic_footprint", &ResolvedFragment::genomic_footprint)
        .def_ro("genomic_start", &ResolvedFragment::genomic_start)
        .def_prop_ro("is_same_strand", &ResolvedFragment::is_same_strand)
        .def_prop_ro("t_inds",
                     [](const ResolvedFragment& r) {
                         nb::set s;
                         for (int32_t t : r.t_inds) s.add(nb::cast(t));
                         return nb::frozenset(s);
                     })
        .def_prop_ro("frag_lengths",
                     [](const ResolvedFragment& r) {
                         return by_candidate(r, r.frag_lengths, [](int32_t length) { return length > 0; });
                     })
        .def_prop_ro("exon_bp",
                     [](const ResolvedFragment& r) { return by_candidate(r, r.exon_bp, [](int32_t) { return true; }); })
        .def_prop_ro("n_gap_hypotheses", &ResolvedFragment::n_gap_hypotheses)
        // Every hypothesis as a list of ``(ref_id, start, end, strand)``, in enumeration order; an empty list is the
        // unspliced (genomic) hypothesis.
        .def_prop_ro("gap_hypotheses", [](const ResolvedFragment& r) {
            nb::list out;
            for (int32_t h = 0; h < r.n_gap_hypotheses(); ++h) {
                nb::list path;
                for (int32_t i = r.gap_intron_offsets[h]; i < r.gap_intron_offsets[h + 1]; ++i) {
                    const IntronBlock& x = r.gap_introns[i];
                    path.append(nb::make_tuple(x.ref_id, x.start, x.end, x.strand));
                }
                out.append(path);
            }
            return out;
        });

    nb::class_<FragmentAccumulator>(m, "FragmentAccumulator")
        .def(nb::init<>(), "Create an empty native accumulator.")
        .def("append", &FragmentAccumulator::append, nb::arg("result"), nb::arg("frag_id"),
             "Append a ResolvedFragment to the accumulator.")
        .def_prop_ro("size", &FragmentAccumulator::size, "Number of fragments in the accumulator.")
        .def("finalize", &FragmentAccumulator::finalize,
             "Finalize to a dict of numpy arrays, leaving the accumulator empty.");

    m.def("resolve_fragment", &resolve_fragment, nb::arg("resolver"), nb::arg("frag"),
          "Resolve a Fragment object to its candidate transcript set (t_inds), as the scan resolves a hit.\n"
          "A candidate need not be compatible with the alignment: the merge falls\n"
          "back to a union when the per-block sets share no transcript.\n\n"
          "Returns a ResolvedFragment, or None when the fragment has no exon blocks.");

    // The SpliceType codes, mirrored for the parity gate against rigel.splice
    m.attr("SPLICE_UNSPLICED") = SPLICE_UNSPLICED;
    m.attr("SPLICE_SPLICED_ANNOT") = SPLICE_SPLICED_ANNOT;
    m.attr("SPLICE_IMPLICIT") = SPLICE_IMPLICIT;
    // The Strand, IntervalType and ZF codes, mirrored for the parity gates against rigel.types and rigel.annotate
    m.attr("STRAND_NONE") = STRAND_NONE;
    m.attr("STRAND_POS") = STRAND_POS;
    m.attr("STRAND_NEG") = STRAND_NEG;
    m.attr("STRAND_AMBIGUOUS") = STRAND_AMBIGUOUS;
    m.attr("ITYPE_EXON") = ITYPE_EXON;
    m.attr("ITYPE_TRANSCRIPT") = ITYPE_TRANSCRIPT;
    m.attr("AF_RESOLVED_BIT") = AF_RESOLVED_BIT;
    m.attr("AF_GDNA_BIT") = AF_GDNA_BIT;
    m.attr("AF_INTERGENIC_BIT") = AF_INTERGENIC_BIT;
    m.attr("AF_MULTIMAPPER_DROP_BIT") = AF_MULTIMAPPER_DROP_BIT;
    m.attr("AF_UNRESOLVED") = AF_UNRESOLVED;
    m.attr("AF_GDNA_INTERGENIC") = AF_GDNA_INTERGENIC;
    m.attr("AF_MULTIMAPPER_DROP") = AF_MULTIMAPPER_DROP;
}

}  // namespace rigel
