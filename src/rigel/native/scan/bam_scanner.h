// bam_scanner.h — the single-pass BAM scan (rigel._bam_impl.BamScanner). A reader thread reads a name-collated BAM
// read pair by read pair; worker threads judge each template, resolve its usable hits into the buffer's chunks and offer
// its placements to the calibration tally; the main thread hands each finished chunk to Python. The stats, the strand
// observations and the calibration payload are returned when the whole BAM is read.
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "calibration/accumulator_set.h"
#include "constants.h"
#include "fragment_blocks.h"
#include "ndarray_util.h"
#include "resolve/fragment_resolver.h"
#include "scan/bam_record.h"
#include "scan/template_hits.h"
#include "thread_queue.h"

namespace rigel {

/// Per-junction sense / antisense counts; uint64, so no count can overflow.
struct SJStrandCounts {
    uint64_t n_sense = 0;
    uint64_t n_antisense = 0;
};

/// A 2×2 of fragment counts, alignment strand × reference strand — the cells of Python's `StrandModel`.
struct StrandContingency {
    int64_t pos_pos = 0;  // align POS, reference POS
    int64_t pos_neg = 0;  // align POS, reference NEG
    int64_t neg_pos = 0;  // align NEG, reference POS
    int64_t neg_neg = 0;  // align NEG, reference NEG

    /// One fragment; the caller has qualified both strands to POS or NEG.
    void add(int32_t align_strand, int32_t ref_strand) {
        const bool ref_pos = ref_strand == STRAND_POS;
        if (align_strand == STRAND_POS) ++(ref_pos ? pos_pos : pos_neg);
        else                            ++(ref_pos ? neg_pos : neg_neg);
    }

    void merge_from(const StrandContingency& o) {
        pos_pos += o.pos_pos;
        pos_neg += o.pos_neg;
        neg_pos += o.neg_pos;
        neg_neg += o.neg_neg;
    }
};

/// The strand models' training observations, over unique-mapper fragments.
struct StrandObservations {
    // The 2×2 of (align_strand, transcript_strand) over fragments whose candidate transcripts share one strand
    StrandContingency exonic;

    // The per-junction strand table, keyed by (ref, start, end, motif strand): each strand-qualified fragment credits
    // its leftmost annotated junction with one sense (align_strand == sj_strand) or antisense observation. Its marginal
    // over motif strand is the spliced strand model's 2×2 (Python's `exonic_spliced`); the table also keeps the
    // dispersion across junctions.
    std::unordered_map<SJKey, SJStrandCounts, SJKeyHash> sj_strand_table;

    /// Add ``o``'s observations, emptying its table.
    void merge_from(StrandObservations& o);
};

struct BamScanStats {
    // Records: each flag counter counts every record carrying the flag; `proper_pair` the primary records flagged as a
    // proper pair
    int64_t total = 0;
    int64_t qc_fail = 0;
    int64_t unmapped = 0;
    int64_t secondary = 0;
    int64_t supplementary = 0;
    int64_t duplicate = 0;
    int64_t proper_pair = 0;

    PairLedger pairs;      // the read pairs and their records as judged, filtered records among them
    RecordCounts records;  // what the records' tags and junctions held

    // Resolution-level
    int64_t n_fragments = 0;
    int64_t n_intergenic = 0;  // templates no hit of which resolves: unspliced, gDNA by construction
    int64_t n_with_exon = 0;
    int64_t n_with_annotated_sj = 0;
    int64_t n_same_strand = 0;
    int64_t n_ambig_strand = 0;

    // Strand model training
    int64_t n_strand_trained = 0;
    int64_t n_strand_skipped_no_sj = 0;
    int64_t n_strand_skipped_ambig_strand = 0;
    int64_t n_strand_skipped_ambiguous = 0;

    // Multimapper
    int64_t n_multimapper_groups = 0;
    int64_t n_multimapper_alignments = 0;

    // Alignments whose records disagree on NM (each keeps the larger)
    int64_t n_nm_conflicting = 0;
    // The library's mismatch rate, read on its unique unspliced hits: their NM, and the read bases they align. A
    // junction's alignment is the aligner's choice of where to cut, and a multimapper's hits are the EM's question.
    int64_t n_unique_unspliced_mismatches = 0;
    int64_t n_unique_unspliced_bases = 0;

    // Hits that are not usable, by `HIT_*` reason
    std::array<int64_t, N_HIT_VERDICTS> hits_rejected{};
    // Templates whose NH exceeds the records reported for either mate: the aligner truncated its output
    int64_t n_nh_above_reported = 0;
    // Read pairs a record of which calls its mate mapped while the run holds no record of that mate: the mate is
    // elsewhere in the BAM (one not grouped by name) or was removed
    int64_t n_read_pairs_lacking_mate = 0;

    // The splice census, `rigel report`'s splice breakdown: one observation per fragment the scanner offers the
    // accumulator, under the resolver's classification of its first hit with a candidate, else of its first hit.
    // Indexed by `SpliceType`, so the increment, the merge and the export cannot miss a category; the names live in
    // `splice_type_label`.
    std::array<int64_t, NUM_SPLICE_TYPES> splice_census{};

    /// One record's flags: each counted on every record carrying it, a proper pair on the primary records only.
    void count_flags(uint16_t flag);
    void merge_from(const BamScanStats& o);
    /// Every counter under its stats-export key.
    nanobind::dict to_dict() const;
};

/// The accumulator's partition and junction CSR as `set_regions` and `set_sj` installed them. Each scan worker
/// deposits into its own set built from it, merged into the scanner's when the BAM is read.
struct ScanPartition {
    std::vector<int64_t> region_bounds;             // flat, reference-major
    std::vector<int64_t> ref_region_bound_offsets;  // n_refs + 1
    std::vector<uint8_t> region_types;              // coarse region type, ref-major, one per region
    int max_length = 0;                             // the fragment-length limit AND the pool-histogram width
    std::vector<int64_t> sj_offsets;                // n_region_bounds + 1, CSR over the flat region_bound axis
    std::vector<int64_t> sj_boundary_right;         // flat region_bound index of each sj's high end
    std::vector<int8_t> sj_strand;                  // each sj's annotated strand
    bool has_sj = false;

    std::unique_ptr<accumulator::AccumulatorSet> make_set() const;
};

/// A record the scan cannot read as paired-end input. Python receives it as `UnpairedRecordError`, a ValueError.
struct UnpairedRecordError : std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

class BamScanner {
public:
    BamScanner(const FragmentResolver& resolver, const std::string& sj_tag_spec, bool skip_duplicates,
               bool include_multimap, int32_t max_fragment_length, int32_t intron_overhang_tolerance);

    /// Install the region partition the accumulator deposits onto, before scan(); `set_sj` follows, a separate call
    /// because this one refuses to run twice. ``ref_region_bound_offsets`` has n_refs + 1 entries; a reference with
    /// offsets[f+1] == offsets[f] has no partition and accepts no deposits. ``max_length`` is the fragment-length limit
    /// applied to L, and the pool-histogram width.
    void set_regions(Array1<int64_t> region_bounds, Array1<int64_t> ref_region_bound_offsets, int32_t n_refs,
                     Array1<uint8_t> region_types, int max_length);

    /// Install the annotated junctions, the flat CSR `build_sj_arrays` emits, keyed by the flat region_bound index.
    /// Required even for an annotation with none (pass empty arrays): scan() refuses to run without it, because a
    /// missing table raises no error elsewhere and every observed intron would read as unannotated.
    void set_sj(Array1<int64_t> offsets, Array1<int64_t> boundary_right, Array1<int8_t> sj_strand);

    /// Read the BAM at ``bam_path``, streaming each chunk of ``chunk_size`` fragments to ``chunk_callback``; the
    /// stats, the strand observations and the calibration payload when the whole BAM is read.
    nanobind::dict scan(const std::string& bam_path, nanobind::callable chunk_callback, int64_t chunk_size,
                        int n_workers, int n_decomp_threads, int qname_batch_size);

private:
    class Worker;
    struct QnameBatch;

    /// scan() needs both the partition and the junction table.
    void require_partition() const;
    /// The reader thread's work: the BAM's read pairs, each with an eligible record a group for the workers.
    void read_bam(const std::string& bam_path, int n_decomp_threads, int qname_batch_size,
                  BoundedQueue<QnameBatch>& input);
    /// Every worker's stats, strand observations and accumulator set, into the scanner's.
    void merge_workers(std::vector<std::unique_ptr<Worker>>& workers);
    nanobind::dict build_result();

    const FragmentResolver* resolver_;
    SJTagMode sj_tag_mode_;
    bool skip_duplicates_;
    bool include_multimap_;
    ResolverPolicy policy_;  // this scan's own settings, never the shared resolver's

    ScanPartition partition_;
    std::unique_ptr<accumulator::AccumulatorSet> acc_set_;  // the workers' sets merge here
    BamScanStats stats_;
    StrandObservations strand_obs_;
};

}  // namespace rigel
