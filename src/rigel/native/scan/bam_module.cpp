// bam_module.cpp — the rigel._bam_impl module (built against htslib): the scan, the annotated-BAM pass, and the
// per-reference accumulator the second pass and the scan cache rebuild.
#include <climits>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "scan/accumulator_python.h"
#include "scan/annotation_writer.h"
#include "scan/bam_scanner.h"
#include "testing/bind_testing.h"

namespace nb = nanobind;
using namespace rigel;
using accumulator::Accumulator;
using accumulator::kNFragmentPools;
using accumulator::kNStrandColumns;

namespace {

/// A read-only numpy view of one field of a struct bank, owned by ``owner``: ``rows`` rows ``row_stride`` elements
/// apart (a struct's size in the field's elements), ``cols`` values each.
template <class T>
nb::object table_view(nb::handle owner, const T* first, size_t rows, size_t cols, int64_t row_stride) {
    return nb::ndarray<nb::numpy, const T, nb::ndim<2>>(first, {rows, cols}, owner, {row_stride, int64_t{1}}).cast();
}

template <class T>
nb::object column_view(nb::handle owner, const T* first, size_t rows, int64_t row_stride) {
    return nb::ndarray<nb::numpy, const T, nb::ndim<1>>(first, {rows}, owner, {row_stride}).cast();
}

/// The elements of a ``T`` field between one row of a ``Row`` bank and the next.
template <class Row, class T>
constexpr int64_t stride_of() {
    return sizeof(Row) / sizeof(T);
}

/// The held bank's arrays (`DeferredFragments`), read in place: placement `p` with a run of its hypotheses becomes
/// the `OfferedFragment` the per-record bindings marshal from the same record, so a batch over the bank answers as
/// those calls would, one Python call per reference instead of one per placement and hypothesis.
struct HeldBank {
    Array1<int64_t> start, end, align_strand, sj_strand, observed_intron_offsets, observed_introns,
        hypothesis_offsets, hypothesis_sj_strand, hypothesis_intron_offsets, hypothesis_introns;
    accumulator::OfferedStorage storage;

    HeldBank(Array1<int64_t> start, Array1<int64_t> end, Array1<int64_t> align_strand, Array1<int64_t> sj_strand,
             Array1<int64_t> observed_intron_offsets, Array1<int64_t> observed_introns,
             Array1<int64_t> hypothesis_offsets, Array1<int64_t> hypothesis_sj_strand,
             Array1<int64_t> hypothesis_intron_offsets, Array1<int64_t> hypothesis_introns)
        : start(start), end(end), align_strand(align_strand), sj_strand(sj_strand),
          observed_intron_offsets(observed_intron_offsets), observed_introns(observed_introns),
          hypothesis_offsets(hypothesis_offsets), hypothesis_sj_strand(hypothesis_sj_strand),
          hypothesis_intron_offsets(hypothesis_intron_offsets), hypothesis_introns(hypothesis_introns) {
        const size_t n = start.shape(0);
        if (end.shape(0) != n || align_strand.shape(0) != n || sj_strand.shape(0) != n ||
            observed_intron_offsets.shape(0) != n + 1 || hypothesis_offsets.shape(0) != n + 1 ||
            hypothesis_intron_offsets.shape(0) != hypothesis_sj_strand.shape(0) + 1)
            throw std::invalid_argument("HeldBank: the arrays do not describe one bank");
        // every coordinate as the per-record bindings cast it, refused where it would not fit
        for (const Array1<int64_t>* introns : {&observed_introns, &hypothesis_introns})
            for (size_t i = 0; i < introns->shape(0); ++i)
                if (introns->data()[i] < INT32_MIN || introns->data()[i] > INT32_MAX)
                    throw std::overflow_error("held bank: an intron coordinate exceeds int32");
    }

    int64_t n_placements() const { return static_cast<int64_t>(start.shape(0)); }

    /// The hypotheses of placement `p`, as a run `[h0, h1)` of the flat hypothesis axis.
    std::pair<int64_t, int64_t> hypotheses_of(int64_t p) const {
        if (p < 0 || p >= n_placements()) throw std::out_of_range("held bank: placement out of range");
        return {hypothesis_offsets.data()[p], hypothesis_offsets.data()[p + 1]};
    }

    /// Placement `p` offered with hypotheses `[h0, h1)`; the pointers stay valid until the next call.
    const accumulator::OfferedFragment& offer(int64_t p, int64_t h0, int64_t h1) {
        const accumulator::HeldArrays bank{start.data(),
                                           end.data(),
                                           align_strand.data(),
                                           sj_strand.data(),
                                           observed_intron_offsets.data(),
                                           observed_introns.data(),
                                           hypothesis_sj_strand.data(),
                                           hypothesis_intron_offsets.data(),
                                           hypothesis_introns.data(),
                                           nullptr,
                                           nullptr};
        return accumulator::rebuild_offered(bank, static_cast<size_t>(p), static_cast<size_t>(h0),
                                            static_cast<size_t>(h1), storage);
    }
};

/// The Accumulator's banks as numpy views, under the specification's `Tally` field names: what the second pass adds
/// up after its drain, and what the scan cache's key hashes.
void bind_accumulator_banks(nb::class_<Accumulator>& c) {
    using accumulator::Boundary;
    using accumulator::Region;
    using accumulator::SpliceJunction;
    const auto acc = [](nb::handle h) -> const Accumulator& { return nb::cast<const Accumulator&>(h); };
    c.def_prop_ro("region_contained_count", [acc](nb::handle h) {
         const Accumulator& a = acc(h);
         return table_view(h, &a.regions_data()[0].contained_count[0], a.n_regions(), kNStrandColumns,
                           stride_of<Region, uint32_t>());
     })
        .def_prop_ro("region_start_count", [acc](nb::handle h) {
            const Accumulator& a = acc(h);
            return table_view(h, a.region_start_count_data(), a.n_regions(), kNStrandColumns,
                              static_cast<int64_t>(kNStrandColumns));
        })
        .def_prop_ro("boundary_unspliced_count", [acc](nb::handle h) {
            const Accumulator& a = acc(h);
            return table_view(h, &a.boundaries_data()[0].unspliced_count[0], a.n_boundaries(), kNStrandColumns,
                              stride_of<Boundary, uint32_t>());
        })
        .def_prop_ro("boundary_spliced_count", [acc](nb::handle h) {
            const Accumulator& a = acc(h);
            return table_view(h, &a.boundaries_data()[0].spliced_count[0], a.n_boundaries(), kNStrandColumns,
                              stride_of<Boundary, uint32_t>());
        })
        .def_prop_ro("boundary_unspliced_inv_length_sum", [acc](nb::handle h) {
            const Accumulator& a = acc(h);
            return column_view(h, &a.boundaries_data()[0].unspliced_inv_length_sum, a.n_boundaries(),
                               stride_of<Boundary, double>());
        })
        // The conserved masses: one value per boundary rather than per strand
        .def_prop_ro("boundary_unspliced_mass", [acc](nb::handle h) {
            const Accumulator& a = acc(h);
            return column_view(h, &a.boundaries_data()[0].unspliced_mass, a.n_boundaries(),
                               stride_of<Boundary, double>());
        })
        .def_prop_ro("boundary_spliced_mass", [acc](nb::handle h) {
            const Accumulator& a = acc(h);
            return column_view(h, &a.boundaries_data()[0].spliced_mass, a.n_boundaries(),
                               stride_of<Boundary, double>());
        })
        .def_prop_ro("sj_count", [acc](nb::handle h) {
            const Accumulator& a = acc(h);
            return table_view(h, &a.sj_data()[0].count[0], a.n_sj(), kNStrandColumns,
                              stride_of<SpliceJunction, uint32_t>());
        })
        .def_prop_ro("sj_inv_length_sum", [acc](nb::handle h) {
            const Accumulator& a = acc(h);
            return column_view(h, &a.sj_data()[0].inv_length_sum, a.n_sj(), stride_of<SpliceJunction, double>());
        })
        // Two columns, unlike every other mass: the junction mass carries a strand, in `sj_count`'s columns
        .def_prop_ro("sj_mass", [acc](nb::handle h) {
            const Accumulator& a = acc(h);
            return table_view(h, &a.sj_data()[0].mass[0], a.n_sj(), kNStrandColumns,
                              stride_of<SpliceJunction, double>());
        })
        .def_prop_ro("pool_lengths", [acc](nb::handle h) {
            const Accumulator& a = acc(h);
            const size_t row = static_cast<size_t>(a.max_length()) + 1;
            return table_view(h, a.pool_lengths_data(), kNFragmentPools, row, static_cast<int64_t>(row));
        })
        .def_prop_ro("deposited_lengths", [acc](nb::handle h) {
            const Accumulator& a = acc(h);
            return column_view(h, a.deposited_lengths_data(), static_cast<size_t>(a.max_length()) + 1, int64_t{1});
        });
}

/// One reference's accumulator, rebuilt by the second pass (one per reference from the payload, re-depositing each
/// drained fragment) and by the scan cache's key. The scan itself deposits in C++.
void bind_accumulator(nb::module_& m) {
    nb::class_<Accumulator> c(m, "Accumulator");
    c.def("__init__",
          [](Accumulator* self, Array1<int64_t> region_bounds, Array1<uint8_t> region_types, int max_length,
             int32_t ref) {
              new (self) Accumulator(
                  std::vector<int64_t>(region_bounds.data(), region_bounds.data() + region_bounds.shape(0)),
                  std::vector<uint8_t>(region_types.data(), region_types.data() + region_types.shape(0)),
                  max_length, ref);
          },
          nb::arg("region_bounds"), nb::arg("region_types"), nb::arg("max_length"), nb::arg("ref"),
          "One reference's sorted region_bound positions, one coarse type per region, the\n"
          "fragment-length limit (which is also the pool-histogram width), and WHICH\n"
          "reference this is — stamped into every deferred record, because the second\n"
          "pass replays those onto that reference's region_bound axis.")
        .def("set_sj",
             [](Accumulator& a, Array1<int32_t> offsets, Array1<int32_t> boundary_right, Array1<int8_t> sj_strand) {
                 a.set_sj(std::vector<int32_t>(offsets.data(), offsets.data() + offsets.shape(0)),
                          std::vector<int32_t>(boundary_right.data(), boundary_right.data() + boundary_right.shape(0)),
                          std::vector<int8_t>(sj_strand.data(), sj_strand.data() + sj_strand.shape(0)));
             },
             nb::arg("offsets"), nb::arg("boundary_right"), nb::arg("sj_strand"),
             "The sj CSR for THIS reference, keyed by the ref-local low-end region_bound\n"
             "index. The sj id is the slot; slot order is a contract.");
    bind_accumulator_banks(c);

    // Returns the QC key the deposit incremented, the specification's `DepositOutcome.value`. `introns` are `(start,
    // end)` pairs on this reference: the scan's adapter restricts a fragment's introns to one reference and
    // de-duplicates them, and `deposit` normalises by coordinate, never reading `IntronBlock::ref_id`.
    c.def("deposit",
          [](Accumulator& a, int64_t start, int64_t end, nb::iterable introns, int32_t align_strand, int32_t sj_strand,
             nb::object hypotheses, int32_t nm, int32_t clip) {
              accumulator::OfferedStorage storage;
              const accumulator::OfferedFragment& frag =
                  marshal_fragment(storage, start, end, introns, align_strand, sj_strand, hypotheses, nm, clip);
              return std::string(accumulator::outcome_key(a.deposit(frag, binding_scratch())));
          },
          nb::arg("start"), nb::arg("end"), nb::arg("observed_introns") = nb::tuple(),
          nb::arg("align_strand") = STRAND_POS, nb::arg("sj_strand") = STRAND_NONE, nb::arg("hypotheses"),
          nb::arg("nm") = 0, nb::arg("clip") = 0)
        // L under one hypothesis, without depositing: what the second pass scores against. Exposed rather than
        // reimplemented in Python, so the scorer and the drain share the deposit's one definition of length.
        .def("length_under",
             [](const Accumulator& a, int64_t start, int64_t end, nb::iterable introns, int32_t align_strand,
                int32_t sj_strand, nb::object hypotheses, size_t hypothesis_index) {
                 accumulator::OfferedStorage storage;
                 const accumulator::OfferedFragment& frag =
                     marshal_fragment(storage, start, end, introns, align_strand, sj_strand, hypotheses);
                 return a.length_under(frag, hypothesis_index, binding_scratch());
             },
             nb::arg("start"), nb::arg("end"), nb::arg("observed_introns") = nb::tuple(),
             nb::arg("align_strand") = STRAND_POS, nb::arg("sj_strand") = STRAND_NONE,
             nb::arg("hypotheses") = nb::tuple(), nb::arg("hypothesis_index") = 0)
        // The annotated junction id for many introns at once: `sj_edge_id`, the deposit's own lookup, per element, so
        // the second pass makes one call per reference. The ids are ref-local slots, the space `set_sj` pins to
        // `ref_sj_offsets`, so a caller on the payload's flat axis adds that base. -1 where the intron is not
        // annotated.
        .def("sj_edge_ids",
             [](const Accumulator& a, Array1<int64_t> starts, Array1<int64_t> ends, Array1<int32_t> sj_strand) {
                 const size_t n = starts.shape(0);
                 if (ends.shape(0) != n || sj_strand.shape(0) != n)
                     throw std::invalid_argument("sj_edge_ids: starts, ends and sj_strand must be aligned per intron");
                 std::vector<int64_t> out(n);
                 a.sj_edge_ids(starts.data(), ends.data(), sj_strand.data(), n, out.data());
                 return vec_to_ndarray(std::move(out));
             },
             nb::arg("starts"), nb::arg("ends"), nb::arg("sj_strand"),
             "The annotated sj id per intron, or -1 — `sj_edge_id` over arrays. Ids are\n"
             "ref-local slots, as `set_sj` defines them.")
        // `(kinds, ids, masses)`, ids ref-local, one entry per object the placement's observed path deposits on: the
        // second pass reads a held placement's path through the deposit's own walk.
        .def("placement_masses",
             [](const Accumulator& a, int64_t start, int64_t end, nb::iterable introns, int32_t align_strand,
                int32_t sj_strand) {
                 accumulator::OfferedStorage storage;
                 const accumulator::OfferedFragment& frag =
                     marshal_fragment(storage, start, end, introns, align_strand, sj_strand, nb::tuple());
                 std::vector<accumulator::PlaceMass> out;
                 a.placement_masses(frag, binding_scratch(), out);
                 std::vector<int8_t> kinds;
                 std::vector<int64_t> ids;
                 std::vector<double> masses;
                 for (const auto& x : out) {
                     kinds.push_back(x.kind);
                     ids.push_back(x.index);
                     masses.push_back(x.mass);
                 }
                 return nb::make_tuple(vec_to_ndarray(std::move(kinds)), vec_to_ndarray(std::move(ids)),
                                       vec_to_ndarray(std::move(masses)));
             },
             nb::arg("start"), nb::arg("end"), nb::arg("observed_introns") = nb::tuple(),
             nb::arg("align_strand") = STRAND_POS, nb::arg("sj_strand") = STRAND_NONE)
        // The second pass's reads of a batch of this reference's held placements, off the held bank in place:
        // `length_under` at every hypothesis of each placement, in the flat axis's order, and each placement's
        // observed path (`placement_masses`' kinds and ref-local ids, CSR over the batch). Returns
        // `(lengths, path_offsets, path_kinds, path_ids)`.
        .def("held_reads",
             [](const Accumulator& a, Array1<int64_t> placements, HeldBank& bank) {
                 accumulator::DepositScratch& scratch = binding_scratch();
                 std::vector<int64_t> lengths, path_offsets{0}, path_ids;
                 std::vector<int8_t> path_kinds;
                 std::vector<accumulator::PlaceMass> masses;
                 for (size_t i = 0; i < placements.shape(0); ++i) {
                     const int64_t p = placements.data()[i];
                     const auto [h0, h1] = bank.hypotheses_of(p);
                     const accumulator::OfferedFragment all = bank.offer(p, h0, h1);
                     for (int64_t k = 0; k < h1 - h0; ++k)
                         lengths.push_back(a.length_under(all, static_cast<size_t>(k), scratch));
                     accumulator::OfferedFragment bare = all;  // the path binding offers no hypothesis
                     bare.n_hypotheses = 0;
                     a.placement_masses(bare, scratch, masses);
                     for (const auto& x : masses) {
                         path_kinds.push_back(x.kind);
                         path_ids.push_back(x.index);
                     }
                     path_offsets.push_back(static_cast<int64_t>(path_ids.size()));
                 }
                 return nb::make_tuple(vec_to_ndarray(std::move(lengths)), vec_to_ndarray(std::move(path_offsets)),
                                       vec_to_ndarray(std::move(path_kinds)), vec_to_ndarray(std::move(path_ids)));
             },
             nb::arg("placements"), nb::arg("bank"))
        // The drain's replay of a batch of this reference's held records, in the order given: placement
        // `placements[i]` deposited with hypothesis `hypotheses[i]` (flat ids) alone, as `deposit` would. Returns
        // the count of each outcome by its key.
        .def("deposit_held",
             [](Accumulator& a, Array1<int64_t> placements, Array1<int64_t> hypotheses, HeldBank& bank) {
                 if (hypotheses.shape(0) != placements.shape(0))
                     throw std::invalid_argument("deposit_held: one hypothesis per placement");
                 accumulator::DepositScratch& scratch = binding_scratch();
                 nb::dict counts;
                 for (size_t i = 0; i < placements.shape(0); ++i) {
                     const int64_t p = placements.data()[i];
                     const int64_t h = hypotheses.data()[i];
                     const auto [h0, h1] = bank.hypotheses_of(p);
                     if (h < h0 || h >= h1) throw std::out_of_range("deposit_held: a hypothesis of another placement");
                     const char* key = accumulator::outcome_key(a.deposit(bank.offer(p, h, h + 1), scratch));
                     counts[key] = (counts.contains(key) ? nb::cast<int64_t>(counts[key]) : int64_t{0}) + 1;
                 }
                 return counts;
             },
             nb::arg("placements"), nb::arg("hypotheses"), nb::arg("bank"));
}

/// The held bank as the batch reads take it: built once from `DeferredFragments`' arrays, which it keeps alive.
void bind_held_bank(nb::module_& m) {
    nb::class_<HeldBank>(m, "HeldBank")
        .def(nb::init<Array1<int64_t>, Array1<int64_t>, Array1<int64_t>, Array1<int64_t>, Array1<int64_t>,
                      Array1<int64_t>, Array1<int64_t>, Array1<int64_t>, Array1<int64_t>, Array1<int64_t>>(),
             nb::arg("start"), nb::arg("end"), nb::arg("align_strand"), nb::arg("sj_strand"),
             nb::arg("observed_intron_offsets"), nb::arg("observed_introns"), nb::arg("hypothesis_offsets"),
             nb::arg("hypothesis_sj_strand"), nb::arg("hypothesis_intron_offsets"), nb::arg("hypothesis_introns"));
}

void bind_scanner(nb::module_& m) {
    nb::exception<UnpairedRecordError>(m, "UnpairedRecordError", PyExc_ValueError);
    nb::class_<BamScanner>(m, "BamScanner")
        .def(nb::init<const FragmentResolver&, const std::string&, bool, bool, int32_t, int32_t>(), nb::arg("ctx"),
             nb::arg("sj_tag_spec"), nb::arg("skip_duplicates"), nb::arg("include_multimap"),
             nb::arg("max_fragment_length"), nb::arg("intron_overhang_tolerance"),
             "Create a BAM scanner.\n\n"
             "Parameters\n"
             "----------\n"
             "ctx : FragmentResolver\n"
             "    The resolve context from index building.\n"
             "sj_tag_spec : str\n"
             "    SJ strand tag specification: 'auto' (XS and ts read on every record, a disagreement\n"
             "    ambiguous), 'XS', 'ts', 'XS,ts', 'ts,XS' (the first present wins), or 'none'. Any other\n"
             "    spec raises ValueError.\n"
             "skip_duplicates : bool\n"
             "    Skip duplicate-flagged reads.\n"
             "include_multimap : bool\n"
             "    Include multimapping reads.\n"
             "max_fragment_length, intron_overhang_tolerance : int\n"
             "    This scan's own pair-rule length limit (>= 1) and intron-overhang tolerance (>= 0).\n",
             nb::keep_alive<1, 2>())  // the scanner keeps its resolver alive
        .def("scan", &BamScanner::scan, nb::arg("bam_path"), nb::arg("chunk_callback"), nb::arg("chunk_size"),
             nb::arg("n_workers"), nb::arg("n_decomp_threads"), nb::arg("qname_batch_size"),
             "Scan BAM file with streaming chunk output.\n\n"
             "Reads the BAM in a background thread, resolves fragments in\n"
             "worker threads, and streams finalized chunks to a Python\n"
             "callback on the main thread.\n\n"
             "Parameters\n"
             "----------\n"
             "bam_path : str\n"
             "    Path to name-sorted BAM file.\n"
             "chunk_callback : callable\n"
             "    Called with a dict of capsule-backed numpy arrays for\n"
             "    each finalized chunk.\n"
             "chunk_size : int\n"
             "    Target number of fragments per chunk.\n"
             "n_workers : int\n"
             "    Number of worker threads.\n"
             "n_decomp_threads : int\n"
             "    Number of htslib BGZF decompression threads.\n"
             "qname_batch_size : int\n"
             "    Read-name groups per input queue item.\n\n"
             "Returns\n"
             "-------\n"
             "dict\n"
             "    Dict with keys: 'stats', 'strand_observations',\n"
             "    'calibration'.\n")
        .def("set_regions", &BamScanner::set_regions, nb::arg("region_bounds"), nb::arg("ref_region_bound_offsets"),
             nb::arg("n_refs"), nb::arg("region_types"), nb::arg("max_length"),
             "Install the accumulator's region partition. Call set_sj next.\n\n"
             "region_bounds : int64[n_region_bounds_total]\n"
             "    Flat sorted region_bound positions for all references. A reference\n"
             "    contributing c region_bounds owns c-1 regions and c-2 interior boundaries.\n"
             "ref_region_bound_offsets : int64[n_refs + 1]\n"
             "    Per-ref offsets into region_bounds; ref f spans\n"
             "    region_bounds[ref_region_bound_offsets[f]:ref_region_bound_offsets[f+1]].\n"
             "n_refs : int\n"
             "    Number of references (matches FragmentResolver).\n"
             "region_types : uint8[n_regions_total]\n"
             "    Coarse region type (0=intergenic, 1=intron, 2=exon), ref-major.\n"
             "    It types the fragment-length pools.\n"
             "max_length : int\n"
             "    The fragment-length limit applied to L, and the width of the\n"
             "    pool histograms. Must be >= 1.\n")
        .def("set_sj", &BamScanner::set_sj, nb::arg("offsets"), nb::arg("boundary_right"), nb::arg("sj_strand"),
             "Install the annotated splice junctions, as build_sj_arrays\n"
             "emits them. A SECOND call because set_regions refuses to run twice.\n\n"
             "offsets : int64[n_region_bounds_total + 1]\n"
             "    CSR over the flat region_bound axis, keyed by the LOW END's region_bound index.\n"
             "boundary_right : int64[n_sj]\n"
             "    Flat region_bound index of each sj's high end. The sj id\n"
             "    IS the slot here; edges_df.edge_row is a join key and is not\n"
             "    passed. Slot order is a contract: sort on\n"
             "    (low-end region_bound, high-end region_bound, sj_strand).\n"
             "sj_strand : int8[n_sj]\n"
             "    Each sj's ANNOTATED strand.\n\n"
             "Not optional: scan() refuses to run without it, because a missing\n"
             "table makes every observed intron read as unannotated and empties\n"
             "the splice junctions silently. Pass empty arrays to declare none.\n");
}

void bind_annotation_writer(nb::module_& m) {
    using Writer = BamAnnotationWriter;
    nb::class_<Writer>(m, "BamAnnotationWriter")
        .def(nb::init<const FragmentResolver&, const std::string&, bool, bool, int32_t, int32_t>(), nb::arg("ctx"),
             nb::arg("sj_tag_spec"), nb::arg("skip_duplicates"), nb::arg("include_multimap"),
             nb::arg("max_fragment_length"), nb::arg("intron_overhang_tolerance"),
             "Create a BAM annotation writer.\n\n"
             "Parameters\n"
             "----------\n"
             "ctx : FragmentResolver\n"
             "    The resolve context from index building.\n"
             "sj_tag_spec : str\n"
             "    SJ strand tag specification (must match pass 1); the specs BamScanner accepts.\n"
             "skip_duplicates : bool\n"
             "    Skip duplicate-flagged reads (must match pass 1).\n"
             "include_multimap : bool\n"
             "    Include multimapping reads (must match pass 1).\n"
             "max_fragment_length, intron_overhang_tolerance : int\n"
             "    Pass 1's pair-rule length limit and intron-overhang tolerance.\n",
             nb::keep_alive<1, 2>())
        .def("write",
             [](const Writer& w, const std::string& bam_path, const std::string& output_path,
                Array1<int64_t> frag_ids, Array1<int32_t> best_tid, Array1<int32_t> best_gid,
                Array1<uint8_t> tx_flags, Array1<float> posterior,
                Array1<int8_t> frag_class, Array1<int32_t> n_candidates,
                Array1<int32_t> locus_id, const std::vector<std::string>& t_ids,
                const std::vector<std::string>& g_ids, const std::vector<std::string>& g_names) {
                 const size_t n = frag_ids.shape(0);
                 for (const size_t len : {best_tid.shape(0), best_gid.shape(0), tx_flags.shape(0), posterior.shape(0),
                                          frag_class.shape(0), n_candidates.shape(0), locus_id.shape(0)})
                     if (len != n) throw std::invalid_argument("write: the annotation arrays are not one length");
                 const Writer::Annotations annotations{frag_ids,   best_tid,     best_gid, tx_flags, posterior,
                                                       frag_class, n_candidates, locus_id, static_cast<int64_t>(n)};
                 return w.write(bam_path, output_path, annotations, t_ids, g_ids, g_names);
             },
             nb::arg("bam_path"), nb::arg("output_path"), nb::arg("ann_frag_ids"), nb::arg("ann_best_tid"),
             nb::arg("ann_best_gid"), nb::arg("ann_tx_flags"), nb::arg("ann_posterior"), nb::arg("ann_frag_class"),
             nb::arg("ann_n_candidates"), nb::arg("ann_locus_id"), nb::arg("t_ids"),
             nb::arg("g_ids"), nb::arg("g_names"),
             "Stamp annotation tags and write BAM.\n\n"
             "Parameters\n"
             "----------\n"
             "bam_path : str\n"
             "    Path to name-sorted input BAM.\n"
             "output_path : str\n"
             "    Path for output annotated BAM.\n"
             "ann_frag_ids .. ann_locus_id : ndarray\n"
             "    Annotation arrays, one entry per fragment, all one length.\n"
             "t_ids, g_ids, g_names : list[str]\n"
             "    Transcript / gene ID / gene name lookup arrays.\n\n"
             "Returns\n"
             "-------\n"
             "dict\n"
             "    The read pairs and records as this pass judged them, counted as the scan counts them, the\n"
             "    digest of its read pairs, and the records written.\n");
}

}  // namespace

NB_MODULE(_bam_impl, m) {
    m.doc() = "C++ BAM scanner and annotation writer for rigel using htslib.\n\n"
              "Provides BamScanner (pass 1: BAM → resolve → buffer) and\n"
              "BamAnnotationWriter (pass 2: stamp tags → write BAM).";

    // The kinds of object a placement's mass lands on (the buffer's place masses, which the place weight reads)
    m.attr("PLACE_REGION") = static_cast<int>(accumulator::kPlaceRegion);
    m.attr("PLACE_BOUNDARY") = static_cast<int>(accumulator::kPlaceBoundary);
    m.attr("PLACE_SJ") = static_cast<int>(accumulator::kPlaceSj);

    bind_held_bank(m);
    bind_accumulator(m);
    bind_scanner(m);
    bind_annotation_writer(m);
    bind_testing(
        m.def_submodule("testing", "The scan's and the accumulator's test seams; production never imports them."));
}
