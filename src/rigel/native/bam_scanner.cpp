/**
 * bam_scanner.cpp — C++ BAM scanner using htslib for rigel.
 *
 * A single C++ scan that:
 *
 *   1. Reads BAM records via htslib
 *   2. Groups records by query name (name-sorted BAM)
 *   3. Parses CIGAR → exon blocks + splice junctions
 *   4. Assembles each hit's mates into one fragment (AssembledFragment)
 *   5. Calls FragmentResolver::_resolve_core() directly
 *   6. Appends results to FragmentAccumulator
 *   7. Offers each read name's placements to the calibration tally at once (AccumulatorSet::offer)
 *   8. Collects model training observations as arrays
 *   9. Streams FragmentAccumulator chunks to Python through `chunk_callback` as they fill, and returns
 *      the stats, the strand observations and the calibration payload once the whole BAM is read
 *
 * Module: rigel._bam_impl
 *
 * Build:
 *   Part of the rigel scikit-build-core build — see CMakeLists.txt.
 *   Requires htslib (linked via -lhts).
 */

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "resolve_context.h"
#include "thread_queue.h"
#include "ndarray_util.h"
#include "calibration/accumulator.h"

namespace nb = nanobind;
using namespace rigel;

// ================================================================
// CIGAR op codes (BAM_C* from htslib/sam.h)
// ================================================================

static constexpr int CIG_MATCH    = BAM_CMATCH;     // 0 M
static constexpr int CIG_DEL      = BAM_CDEL;       // 2 D
static constexpr int CIG_REF_SKIP = BAM_CREF_SKIP;  // 3 N
static constexpr int CIG_EQUAL    = BAM_CEQUAL;     // 7 =
static constexpr int CIG_DIFF     = BAM_CDIFF;      // 8 X

// Does this CIGAR op advance the reference position?
static inline bool cigar_advances_ref(int op) {
    return op == CIG_MATCH || op == CIG_DEL ||
           op == CIG_EQUAL || op == CIG_DIFF;
}

// ================================================================
// Streaming-scan capacities: bounds on memory and hand-off, not model parameters
// ================================================================

// Input-queue capacity per worker: the read-name batches the reader may hold queued before it blocks.
static constexpr int INPUT_BATCHES_PER_WORKER = 4;
// Output-queue capacity per worker: the finished chunks the workers may hold queued before they block on
// the main thread.
static constexpr int OUTPUT_CHUNKS_PER_WORKER = 2;

// Pre-size a worker's chunk: chunk_size fragment rows and 3/2 candidate-transcript rows per fragment. A
// capacity hint only; the columns grow past it.
static inline void reserve_chunk(FragmentAccumulator& acc, int64_t chunk_size) {
    acc.reserve(chunk_size, chunk_size * 3 / 2);
}

// ================================================================
// SJ strand tag configuration
// ================================================================

enum class SJTagMode : uint8_t {
    NONE = 0,     // No SJ strand tag
    XS_ONLY = 1,  // Check XS tag only
    TS_ONLY = 2,  // Check ts tag only (minimap2)
    XS_TS = 3,    // Check XS first, then ts
    TS_XS = 4,    // Check ts first, then XS
    AUTO = 5,     // Read both on every record; a disagreement is AMBIGUOUS
};

// Parse SJ tag configuration string → SJTagMode enum. Any other spec is refused.
static SJTagMode parse_sj_tag_spec(const std::string& spec) {
    if (spec == "auto")  return SJTagMode::AUTO;
    if (spec == "XS")    return SJTagMode::XS_ONLY;
    if (spec == "ts")    return SJTagMode::TS_ONLY;
    if (spec == "XS,ts") return SJTagMode::XS_TS;
    if (spec == "ts,XS") return SJTagMode::TS_XS;
    if (spec == "none" || spec.empty()) return SJTagMode::NONE;
    throw std::invalid_argument(
        "sj_tag_spec must be 'auto', 'XS', 'ts', 'XS,ts', 'ts,XS' or 'none'; got '" + spec + "'");
}

// What the junction-strand tag said on each spliced record of the scan (the annotation pass reads the same tags
// and counts nothing). A tag holding neither '+' nor '-' reads as absent; under AUTO a record whose XS and ts both
// give a strand and disagree is AMBIGUOUS.
struct SJTagCounts {
    int64_t tagged = 0;       // spliced records given a strand
    int64_t untagged = 0;     // spliced records no tag gave a strand
    int64_t conflicting = 0;  // spliced records whose XS and ts disagree (AUTO only)
    int64_t invalid = 0;      // XS or ts tags holding neither '+' nor '-'

    void merge(const SJTagCounts& o) {
        tagged += o.tagged;
        untagged += o.untagged;
        conflicting += o.conflicting;
        invalid += o.invalid;
    }
};

// An integer tag on each record of the scan (NM, NH): absent, or held as anything but an integer of 0..2^31-1 (read
// as absent). The annotation pass reads the same tags and counts nothing.
struct IntTagCounts {
    int64_t missing = 0;
    int64_t invalid = 0;

    void merge(const IntTagCounts& o) {
        missing += o.missing;
        invalid += o.invalid;
    }
};

// ================================================================
// Fragment class / pool / splice-type label helpers
// ================================================================

static const char* frag_class_label(int code) {
    switch (code) {
        case FRAG_UNAMBIG:           return "unambig";
        case FRAG_AMBIG_SAME_STRAND: return "ambig_same_strand";
        case FRAG_AMBIG_OPP_STRAND:  return "ambig_opp_strand";
        case FRAG_MULTIMAPPER:       return "multimapper";
        // Chimeric / intergenic fragments never entered EM; ZC reports the
        // *input-ambiguity* axis only and is blanked for these cases.  The
        // outcome is carried entirely by ZF (AF_CHIMERIC / AF_GDNA_INTERGENIC).
        case FRAG_CHIMERIC:          return ".";
        case -1:                     return ".";
        default:                     return "unknown";
    }
}

static const char* splice_type_label(int code) {
    switch (code) {
        case SPLICE_UNSPLICED:       return "unspliced";
        case SPLICE_SPLICED_ANNOT:   return "spliced_annot";
        case SPLICE_IMPLICIT:        return "spliced_implicit";
        case SPLICE_ARTIFACT:        return "splice_artifact";
        default:                     return "unknown";
    }
}

// Why a hit is not usable (`HIT_*`), as the stats export names it.
static const char* hit_verdict_label(int32_t verdict) {
    switch (verdict) {
        case HIT_REFERENCE:      return "reference";
        case HIT_STRAND:         return "strand";
        case HIT_OUTWARD:        return "outward";
        case HIT_LENGTH:         return "length";
        case HIT_MISSING_MATE:   return "missing_mate";
        case HIT_CONTRADICTORY:  return "contradictory";
        case HIT_NOVEL_JUNCTION: return "novel_junction";
        case HIT_UNSUPPORTED:    return "unsupported";
        case HIT_OFF_INDEX:      return "off_index";
        case HIT_CHIMERIC:       return "chimeric";
        default:                 return "unknown";
    }
}

// Why a template is withheld: it has candidate hits and none is usable (`judge_hits`). The label is the annotated
// BAM's ZC and the stats export's name.
static constexpr int32_t WITHHELD_NONE               = 0;
static constexpr int32_t WITHHELD_UNPAIRED           = 1;  // none of its hits is a molecule
static constexpr int32_t WITHHELD_UNANNOTATED_SPLICE = 2;  // the annotated-path rule withheld one of its hits
static constexpr int32_t WITHHELD_OFF_INDEX          = 3;  // one of its hits lies on a reference outside the index
static constexpr int32_t WITHHELD_CHIMERIC           = 4;  // one of its hits has a mate aligned in pieces
static constexpr int32_t N_WITHHELD                  = 5;

static const char* withheld_label(int32_t outcome) {
    switch (outcome) {
        case WITHHELD_UNPAIRED:           return "unpaired";
        case WITHHELD_UNANNOTATED_SPLICE: return "unannotated_splice";
        case WITHHELD_OFF_INDEX:          return "off_index";
        case WITHHELD_CHIMERIC:           return "chimeric";
        default:                          return "unknown";
    }
}

// ================================================================
// A read error is never the end of the file
// ================================================================

// htslib's sam_read1 returns -1 at the end of the file and less than -1 on an error. A loop that stops on either
// would quantify a truncated or corrupt BAM, or a CRAM whose reference cannot be found, on the records read before
// the error, with no sign of it; each read loop keeps the return and throws this instead.
static std::runtime_error read_error(const std::string& path, int ret, int64_t n_records) {
    return std::runtime_error("Failed to read " + path + " after " + std::to_string(n_records) +
                              " records (htslib error " + std::to_string(ret) +
                              "): a truncated or corrupt BAM, or a CRAM whose reference cannot be found");
}

// ================================================================
// Build htslib tid → internal ref_id mapping
// ================================================================

static void build_tid_mapping(
    const bam_hdr_t* hdr,
    const std::unordered_map<std::string, int32_t>& ref_to_id,
    std::vector<int32_t>& tid_to_ref_id)
{
    int32_t n_targets = hdr->n_targets;
    tid_to_ref_id.resize(n_targets, -1);
    for (int32_t i = 0; i < n_targets; i++) {
        auto it = ref_to_id.find(hdr->target_name[i]);
        if (it != ref_to_id.end()) {
            tid_to_ref_id[i] = it->second;
        }
    }
}

// ================================================================
// Lightweight record for one BAM alignment in a qname group
// ================================================================

/// CIGAR-derived splice junction with per-read anchor lengths.
///
/// The anchor on each side of the sj is the number of
/// reference-advancing CIGAR bases between that side and the nearest
/// clip (or the end of the alignment).  Both anchors are required to
/// test against the artifact blacklist.
struct SJCigarEntry {
    int32_t start;
    int32_t end;
    int32_t strand;
    int32_t anchor_left;
    int32_t anchor_right;
};

struct ParsedAlignment {
    int32_t ref_id;         // the index's reference id; -1 for a reference outside the index
    int32_t tid = -1;       // the BAM's reference id
    int32_t ref_start;      // pos
    uint16_t flag;
    // Per-record count of CIGAR sj dropped by the splice-artifact
    // blacklist.  Populated by filter_blacklisted_sjs in Pass 1 (summed
    // per fragment for the resolver's SPLICE_ARTIFACT promotion; the
    // n_sj_blacklisted stat is accumulated through parse_bam_record's
    // pointer argument, not from here) and Pass 2 (for the ZB:i
    // annotated-BAM tag).
    // Fits into the padding after `flag` — zero effective memory growth.
    uint16_t n_sj_blacklisted = 0;
    int32_t nm;             // NM tag (edit distance), 0 if absent
    int32_t nh;             // NH tag, 1 if absent

    // CIGAR-parsed exon blocks and splice junctions
    std::vector<std::pair<int32_t, int32_t>> exons;  // (start, end)
    std::vector<SJCigarEntry> sjs;

    // The soft clips at the CIGAR's two ends, which the pair validator walks from the read's 5' end
    int32_t clip_lead = 0;
    int32_t clip_trail = 0;

    // The CIGAR as the alignment's identity (`canonical_cigar`)
    std::vector<uint32_t> cigar;

    // Whether the record carries SA, naming the other pieces of a chimeric alignment; and, on a location's
    // representative, whether the location is one such piece (`template_hits`)
    bool has_sa = false;
    bool chimeric = false;

    bool is_read1()         const { return flag & BAM_FREAD1; }
    bool is_secondary()     const { return flag & BAM_FSECONDARY; }
    bool is_supplementary() const { return flag & BAM_FSUPPLEMENTARY; }
    bool is_reverse()       const { return flag & BAM_FREVERSE; }
    bool is_proper_pair()   const { return flag & BAM_FPROPER_PAIR; }
};

// ================================================================
// QnameGroup — all BAM alignments for one query name, pre-parsed
// ================================================================
//
// The producer thread reads raw htslib bam1_t records, immediately
// parses each into a ParsedAlignment, and packages them here.
// Workers receive fully parsed value types — no raw bam1_t* pointers
// cross the thread boundary.

struct QnameGroup {
    std::vector<ParsedAlignment> records;
    int64_t frag_id = 0;
};

struct QnameBatch {
    std::vector<QnameGroup> groups;
};

// ================================================================
// Model training observation collectors
// ================================================================

/// Per-sj sense / antisense counts, sized as counts: uint64 removes the
/// overflow question permanently at a cost of ~2.4 MB over a human sj set.
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

struct StrandObservations {
    // exonic: the 2×2 of (align_strand, transcript_strand) over unique-mapper, non-chimeric fragments
    // whose candidate transcripts share one strand
    StrandContingency exonic;

    // ⭐ The per-sj strand table.
    // A sj is uniquely specified by (ref, start, end, motif strand); each
    // strand-qualified fragment credits its leftmost ANNOTATED sj with one
    // sense (align_strand == sj_strand) or antisense observation.
    //
    // ⭐ This is a strict REFINEMENT of the spliced strand model's 2×2 (Python's
    // `exonic_spliced`), which is exactly its marginal: summing n_sense over
    // motif-POS sj gives pos_pos, n_antisense over motif-POS gives neg_pos, and
    // symmetrically for motif-NEG.  Populated from the SAME
    // `get_is_strand_qualified()` branch, one observation per fragment, so the
    // identity is unconditional.  It exists because a dispersion ACROSS sj cannot
    // be recovered from the 2×2.
    std::unordered_map<SJKey, SJStrandCounts, SJKeyHash> sj_strand_table;
};

// ================================================================
// Stats counters
// ================================================================

struct BamScanStats {
    // BAM-level
    int64_t total = 0;
    int64_t qc_fail = 0;
    int64_t unmapped = 0;
    int64_t secondary = 0;
    int64_t supplementary = 0;
    int64_t duplicate = 0;
    int64_t n_read_names = 0;
    int64_t unique = 0;
    int64_t multimapping = 0;
    int64_t proper_pair = 0;
    int64_t improper_pair = 0;
    int64_t mate_unmapped = 0;

    // Resolution-level
    int64_t n_fragments = 0;
    int64_t n_chimeric = 0;
    int64_t n_chimeric_trans = 0;
    int64_t n_chimeric_cis_strand_same = 0;
    int64_t n_chimeric_cis_strand_diff = 0;
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

    // Splice-junction artifact blacklist (per read, across all alignments)
    int64_t n_sj_observed = 0;
    int64_t n_sj_blacklisted = 0;

    // The junction-strand tag on each spliced record; NM and NH on each record; and the alignments whose records
    // disagree on NM (each keeps the larger)
    SJTagCounts sj_tags;
    IntTagCounts nm_tag, nh_tag;
    int64_t n_nm_conflicting = 0;

    // Hits that are not usable, by `HIT_*` reason, and the templates left with none, by `WITHHELD_*` outcome
    std::array<int64_t, N_HIT_VERDICTS> hits_rejected{};
    std::array<int64_t, N_WITHHELD> withheld{};
    // Templates whose NH exceeds the records reported for either mate: the aligner truncated its output
    int64_t n_nh_above_reported = 0;

    // ── the per-fragment splice census ────────────────────────────────────────────────────────────
    //
    // ⭐ ONE observation per HIT the scanner offers the accumulator (a multimapper's hits one by one), filed
    // under the resolver's own classification. This IS `rigel report`'s splice breakdown, which used to be
    // read off the fragment-length CATEGORY MODELS — so it counted only the fragments that also yielded a
    // length observation, a population nobody had stated. The scanner classifies, so the scanner counts.
    //
    // ⚠ An ARRAY and not five named fields, deliberately: `census[st]++` is one statement, the merge
    // is one loop, the export is one loop, and a category added to `SpliceType` cannot slip past any
    // of the three. The five names exist once, in `splice_type_label`.
    std::array<int64_t, NUM_SPLICE_TYPES> splice_census{};

    // Censused, then never offered to `deposit()`: the deposit adapter could not express the
    // fragment as ONE molecule on ONE region_bound axis (chiefly blocks on more than one reference), or no
    // accumulator was installed. ⭐ It exists to close the books — without it those `return`s are a
    // silent loss, which is the pattern this accumulator was rewritten to delete. With it:
    //
    //   Σ census − census[SPLICE_ARTIFACT] == qc.deposited + qc.deferred_undetermined_gap
    //       + qc.deferred_multiple_placements + Σ qc.dropped_* + n_deposit_not_offered
    //
    // ⚠ Exact while no read name offers more than one placement. The census and this counter count
    // assembled hits and the qc counts one offer per read name, so a multimapper offering k placements
    // adds k on the left and 1 on the right.
    //
    // ⚠ Artifacts are subtracted rather than counted here because holding them out is a POSITIVE
    // scanner decision about the data, not a failure to represent it — and `census[SPLICE_ARTIFACT]`
    // already reports it exactly.
    int64_t n_deposit_not_offered = 0;

    /// Add all counters from another stats struct (for merge).
    void merge_from(const BamScanStats& o) {
        total += o.total;
        qc_fail += o.qc_fail;
        unmapped += o.unmapped;
        secondary += o.secondary;
        supplementary += o.supplementary;
        duplicate += o.duplicate;
        n_read_names += o.n_read_names;
        unique += o.unique;
        multimapping += o.multimapping;
        proper_pair += o.proper_pair;
        improper_pair += o.improper_pair;
        mate_unmapped += o.mate_unmapped;
        n_fragments += o.n_fragments;
        n_chimeric += o.n_chimeric;
        n_chimeric_trans += o.n_chimeric_trans;
        n_chimeric_cis_strand_same += o.n_chimeric_cis_strand_same;
        n_chimeric_cis_strand_diff += o.n_chimeric_cis_strand_diff;
        n_intergenic += o.n_intergenic;
        n_with_exon += o.n_with_exon;
        n_with_annotated_sj += o.n_with_annotated_sj;
        n_same_strand += o.n_same_strand;
        n_ambig_strand += o.n_ambig_strand;
        n_strand_trained += o.n_strand_trained;
        n_strand_skipped_no_sj += o.n_strand_skipped_no_sj;
        n_strand_skipped_ambig_strand += o.n_strand_skipped_ambig_strand;
        n_strand_skipped_ambiguous += o.n_strand_skipped_ambiguous;
        n_multimapper_groups += o.n_multimapper_groups;
        n_multimapper_alignments += o.n_multimapper_alignments;
        n_sj_observed += o.n_sj_observed;
        n_sj_blacklisted += o.n_sj_blacklisted;
        sj_tags.merge(o.sj_tags);
        nm_tag.merge(o.nm_tag);
        nh_tag.merge(o.nh_tag);
        n_nm_conflicting += o.n_nm_conflicting;
        for (int32_t v = 0; v < N_HIT_VERDICTS; v++) hits_rejected[v] += o.hits_rejected[v];
        for (int32_t w = 0; w < N_WITHHELD; w++) withheld[w] += o.withheld[w];
        n_nh_above_reported += o.n_nh_above_reported;
        for (size_t i = 0; i < splice_census.size(); i++)
            splice_census[i] += o.splice_census[i];
        n_deposit_not_offered += o.n_deposit_not_offered;
    }
};

// ================================================================
// AssembledFragment — one hit's mates, merged into exon blocks and introns
// ================================================================

struct AssembledFragment {
    std::vector<ExonBlock> exons;
    std::vector<IntronBlock> introns;
    int32_t nm = 0;

    int32_t genomic_footprint() const {
        if (exons.empty()) return -1;
        int32_t min_start = exons[0].start;
        int32_t max_end = exons[0].end;
        for (size_t i = 1; i < exons.size(); i++) {
            min_start = std::min(min_start, exons[i].start);
            max_end = std::max(max_end, exons[i].end);
        }
        return max_end - min_start;
    }

    bool has_introns() const { return !introns.empty(); }
};

// ================================================================
// WorkerState — per-thread mutable state for parallel scan
// ================================================================

struct WorkerState {
    ResolverScratch scratch;
    FragmentAccumulator accumulator;
    BamScanStats stats;
    StrandObservations strand_obs;
    // Per-worker accumulator set (one Accumulator per ref). Owned
    // here so worker threads never touch the shared template. Null when no
    // region partition was provided.
    std::unique_ptr<rigel::accumulator::AccumulatorSet> acc_set;

    // Reusable scratch so the deposit path allocates no per-fragment memory. ⭐ Measured on the shipped
    // accumulator: one per-fragment std::vector cost 22.8 ns, 18 % of the deposit, and was invisible to a
    // sampling profiler because the time is attributed to `malloc`. Capacity amortizes to the high-water
    // mark and is never released.
    //
    // `offer_scratch` is the accumulator set's (its `deposit` member holds the accumulator's normalised
    // introns, path segments, sj ids) and is opaque here. The placement arenas below are the ADAPTER's: a
    // fragment's introns restricted to the reference of the hit being collected and de-duplicated, which is
    // what `OfferedFragment::observed_introns` points at.

    // ⭐ ONE OFFER PER FRAGMENT. Every non-chimeric hit of a fragment is collected as a placement — its
    // observed introns and its gap hypotheses' implied introns COPIED into these arenas, because the
    // resolver's result dies with the hit — and after the last hit the whole set is offered to the
    // accumulator set at once (`AccumulatorSet::offer`), which is `deposit` for a unique mapper and the
    // multimapper rule for the rest. Spans point into the arenas, so they are resolved from offsets only
    // once nothing more will grow.
    struct PendingPlacement {
        int32_t ref_id;
        int64_t start, end;
        int32_t align_strand, sj_strand;
        std::size_t obs0, obs1, hyp0, hyp1;
    };
    struct PendingHypothesis {
        std::size_t imp0, imp1;
        int32_t sj_strand;
    };
    std::vector<PendingPlacement>  pending;
    std::vector<PendingHypothesis> pending_hypotheses;
    std::vector<IntronBlock>       arena_observed, arena_implied;
    std::vector<rigel::accumulator::GapHypothesis> arena_spans;
    std::vector<rigel::accumulator::Placement>     placements;
    rigel::accumulator::OfferScratch               offer_scratch;

    explicit WorkerState(int32_t n_transcripts)
        : scratch(n_transcripts) {}
};

// Merge observations
static void merge_strand_obs(StrandObservations& dst, StrandObservations& src) {
    dst.exonic.merge_from(src.exonic);
    for (const auto& [key, counts] : src.sj_strand_table) {
        auto& dst_counts = dst.sj_strand_table[key];
        dst_counts.n_sense     += counts.n_sense;
        dst_counts.n_antisense += counts.n_antisense;
    }
    src.sj_strand_table.clear();
}

// ================================================================
// CIGAR parsing
// ================================================================
//
// In addition to emitting exon blocks and splice junctions, we compute
// per-sj anchor lengths (left, right).  An "anchor" is the count
// of reference-advancing CIGAR bases on one side of the N op, bounded
// by clips/ends.  Insertions never contribute; deletions do (they
// occupy reference positions in the aligner's evidence for the anchor).
//
// Example: CIGAR 10M 500N 15M 1000N 125M
//   sj 1 (500N):  anchor_left = 10,  anchor_right = 15 + 125 = 140
//   sj 2 (1000N): anchor_left = 10 + 15 = 25, anchor_right = 125

/// The soft clips at the CIGAR's first and last operation.
static void read_soft_clips(const bam1_t* b, int32_t& lead, int32_t& trail) {
    lead = trail = 0;
    const int32_t n_cigar = b->core.n_cigar;
    if (n_cigar == 0) return;
    const uint32_t* cigar = bam_get_cigar(b);
    if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) lead = bam_cigar_oplen(cigar[0]);
    if (bam_cigar_op(cigar[n_cigar - 1]) == BAM_CSOFT_CLIP) trail = bam_cigar_oplen(cigar[n_cigar - 1]);
}

static void parse_cigar(
    const bam1_t* b,
    std::vector<std::pair<int32_t, int32_t>>& exons,
    std::vector<SJCigarEntry>& sjs)
{
    exons.clear();
    sjs.clear();

    int32_t n_cigar = b->core.n_cigar;
    if (n_cigar == 0) return;

    const uint32_t* cigar = bam_get_cigar(b);
    int32_t pos = b->core.pos;     // 0-based genomic position
    int32_t start = pos;           // current exon start
    int32_t ref_advanced = 0;      // ref-advancing bases seen so far

    // First pass: build exons, record each sj's left anchor and
    // its CIGAR intron length so the right anchor can be back-filled.
    for (int32_t i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        if (op == CIG_REF_SKIP) {
            if (pos > start) exons.push_back({start, pos});
            SJCigarEntry sj{};
            sj.start = pos;
            sj.end = pos + len;
            sj.strand = STRAND_NONE;
            sj.anchor_left = ref_advanced;
            sj.anchor_right = 0;  // patched below
            sjs.push_back(sj);
            start = pos + len;
            pos = start;
        } else if (cigar_advances_ref(op)) {
            pos += len;
            ref_advanced += len;
        }
        // INS, SOFT_CLIP, HARD_CLIP, PAD: do not advance the reference
        // and therefore do not contribute to anchor lengths.
    }
    if (pos > start) exons.push_back({start, pos});

    // Back-fill right anchors: right = total - left.  ``ref_advanced`` is
    // the grand total of ref-advancing bases across the whole CIGAR.
    for (auto& sj : sjs) {
        sj.anchor_right = ref_advanced - sj.anchor_left;
    }
}

// ================================================================
// Splice-junction artifact filtering
// ================================================================
//
// Drop each CIGAR splice junction whose (left, right) anchors sit
// inside the blacklist's observed artifact envelope.  Rejection rule:
// either anchor <= max → sj is indistinguishable from known
// gDNA-derived artifacts and is treated as unspliced.
//
// Returns the number of sj dropped (for scanner statistics).
static inline int32_t filter_blacklisted_sjs(
    std::vector<SJCigarEntry>& sjs,
    const FragmentResolver& resolver,
    int32_t ref_id)
{
    if (sjs.empty() || !resolver.has_sj_blacklist() || ref_id < 0) return 0;
    int32_t n_before = static_cast<int32_t>(sjs.size());
    sjs.erase(std::remove_if(sjs.begin(), sjs.end(),
        [&](const SJCigarEntry& sj) {
            const auto* hit = resolver.sj_blacklist_lookup(
                ref_id, sj.start, sj.end);
            if (!hit) return false;
            return sj.anchor_left  <= hit->max_anchor_left
                || sj.anchor_right <= hit->max_anchor_right;
        }), sjs.end());
    return n_before - static_cast<int32_t>(sjs.size());
}

// ================================================================
// Read SJ strand tag from a BAM record
// ================================================================

static int32_t read_sj_strand(const bam1_t* b, SJTagMode mode, SJTagCounts* counts) {
    if (mode == SJTagMode::NONE) {
        if (counts) ++counts->untagged;
        return STRAND_NONE;
    }

    auto try_tag = [&](const char* tag, bool is_ts) -> int32_t {  // a strand, or -1 when absent or unreadable
        uint8_t* aux = bam_aux_get(b, tag);
        if (!aux) return -1;  // tag not present
        // Dispatch on the htslib type code (*aux) to call the correct
        // accessor.  Type 'A' = single char, type 'Z'/'H' = string
        // (take first character).  Any other type or value reads as absent, counted.
        char val = 0;
        if (*aux == 'A') {
            val = bam_aux2A(aux);
        } else if (*aux == 'Z' || *aux == 'H') {
            const char* s = bam_aux2Z(aux);
            if (s && s[0]) val = s[0];
        }
        int32_t strand;
        if (val == '+') strand = STRAND_POS;
        else if (val == '-') strand = STRAND_NEG;
        else {
            if (counts) ++counts->invalid;
            return -1;
        }

        // minimap2's 'ts' tag is alignment-relative: flip for reverse reads
        if (is_ts && (b->core.flag & BAM_FREVERSE)) {
            strand = (strand == STRAND_POS) ? STRAND_NEG : STRAND_POS;
        }
        return strand;
    };

    int32_t result = -1;

    switch (mode) {
        case SJTagMode::XS_ONLY:
            result = try_tag("XS", false);
            break;
        case SJTagMode::TS_ONLY:
            result = try_tag("ts", true);
            break;
        case SJTagMode::XS_TS:
            result = try_tag("XS", false);
            if (result < 0) result = try_tag("ts", true);
            break;
        case SJTagMode::TS_XS:
            result = try_tag("ts", true);
            if (result < 0) result = try_tag("XS", false);
            break;
        case SJTagMode::AUTO: {
            // Both tags, every record: neither name proves which is right, so a disagreement is AMBIGUOUS
            const int32_t xs = try_tag("XS", false);
            const int32_t ts = try_tag("ts", true);
            if (xs >= 0 && ts >= 0 && xs != ts) {
                if (counts) ++counts->conflicting;
                return STRAND_AMBIGUOUS;
            }
            result = (xs >= 0) ? xs : ts;
            break;
        }
        default:
            break;
    }

    if (result < 0) {
        if (counts) ++counts->untagged;
        return STRAND_NONE;
    }
    if (counts) ++counts->tagged;
    return result;
}

/// A record's CIGAR as its alignment's identity: M, = and X read as one operation, adjacent runs of it merged (`50=50X`
/// and `100M` are one alignment); every other operation, soft and hard clips included, as written.
static void canonical_cigar(const bam1_t* b, std::vector<uint32_t>& out)
{
    out.clear();
    const uint32_t* cigar = bam_get_cigar(b);
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        const uint32_t len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CEQUAL || op == BAM_CDIFF) op = BAM_CMATCH;
        if (op == BAM_CMATCH && !out.empty() && bam_cigar_op(out.back()) == BAM_CMATCH)
            out.back() = bam_cigar_gen(bam_cigar_oplen(out.back()) + len, BAM_CMATCH);
        else
            out.push_back(bam_cigar_gen(len, op));
    }
}

/// A record's integer tag ``tag``: its value, or ``absent`` when the record lacks it or holds it as anything but an
/// integer of 0..2^31-1, each counted in ``counts``.
static int32_t read_int_tag(const bam1_t* b, const char* tag, int32_t absent, IntTagCounts* counts)
{
    const uint8_t* aux = bam_aux_get(b, tag);
    if (!aux) {
        if (counts) ++counts->missing;
        return absent;
    }
    int64_t v = -1;
    switch (*aux) {
        case 'c': case 'C': case 's': case 'S': case 'i': case 'I': v = bam_aux2i(aux); break;
        default: break;
    }
    if (v < 0 || v > std::numeric_limits<int32_t>::max()) {
        if (counts) ++counts->invalid;
        return absent;
    }
    return static_cast<int32_t>(v);
}

// ================================================================
// Parse a raw htslib bam1_t into a ParsedAlignment value type
// ================================================================
// Extracts all fields the pipeline needs (scalars, tags, CIGAR-derived
// exon blocks and splice junctions) so that downstream code never
// touches bam1_t directly. ``counts`` takes the per-record counts (the
// scan's); the annotation pass, which reads the same records again, passes null.

static ParsedAlignment parse_bam_record(
    const bam1_t* b,
    const std::vector<int32_t>& tid_to_ref_id,
    SJTagMode sj_tag_mode,
    const FragmentResolver* resolver,
    BamScanStats* counts)
{
    ParsedAlignment rec;
    rec.flag = b->core.flag;
    rec.tid = b->core.tid;
    rec.ref_id = (rec.tid >= 0 && rec.tid < static_cast<int32_t>(tid_to_ref_id.size())) ? tid_to_ref_id[rec.tid] : -1;
    rec.ref_start = b->core.pos;
    canonical_cigar(b, rec.cigar);
    rec.has_sa = bam_aux_get(b, "SA") != nullptr;

    // BAM aux tags: a missing NM reads 0, a disclosed approximation (counted, and warned of); a missing NH reads 1
    rec.nm = read_int_tag(b, "NM", 0, counts ? &counts->nm_tag : nullptr);
    rec.nh = read_int_tag(b, "NH", 1, counts ? &counts->nh_tag : nullptr);

    // CIGAR → exon blocks + splice junctions
    parse_cigar(b, rec.exons, rec.sjs);
    read_soft_clips(b, rec.clip_lead, rec.clip_trail);
    if (!rec.sjs.empty()) {
        const int32_t sj_strand = read_sj_strand(b, sj_tag_mode, counts ? &counts->sj_tags : nullptr);
        for (auto& sj : rec.sjs) sj.strand = sj_strand;
    }

    // Splice-junction artifact filtering.  Per-read, immediately after
    // CIGAR parsing.  Downstream fragment assembly unions across mates,
    // so this is the permissive rule: a sj survives as long as at
    // least one read supports it with anchors above the blacklist envelope.
    if (counts) counts->n_sj_observed += static_cast<int64_t>(rec.sjs.size());
    if (resolver && resolver->has_sj_blacklist() && !rec.sjs.empty()) {
        const int32_t n_dropped = filter_blacklisted_sjs(rec.sjs, *resolver, rec.ref_id);
        rec.n_sj_blacklisted = static_cast<uint16_t>(n_dropped);
        if (counts) counts->n_sj_blacklisted += n_dropped;
    }

    return rec;
}

// ================================================================
// A hit's junctions, and its fragment
// ================================================================

/// Every junction a hit observes: the kept junctions of both mates' records, ascending by reference, coordinates and
/// motif strand, each once. The pair validator and the fragment the resolver reads take this one set.
static void observed_junctions(
    const std::vector<ParsedAlignment*>& r1_reads,
    const std::vector<ParsedAlignment*>& r2_reads,
    std::vector<IntronBlock>& out)
{
    out.clear();
    for (const auto* reads : {&r1_reads, &r2_reads})
        for (const auto* rec : *reads)
            for (const auto& sj : rec->sjs) out.push_back({rec->ref_id, sj.start, sj.end, sj.strand});
    const auto key = [](const IntronBlock& x) { return std::tie(x.ref_id, x.start, x.end, x.strand); };
    std::sort(out.begin(), out.end(), [&](const IntronBlock& x, const IntronBlock& y) { return key(x) < key(y); });
    out.erase(std::unique(out.begin(), out.end(),
                          [&](const IntronBlock& x, const IntronBlock& y) { return key(x) == key(y); }),
              out.end());
}

/// A hit's records as one fragment: their aligned blocks merged per reference and strand, R2's strand flipped, so a
/// block carries the fragment's strand; the hit's junctions (`observed_junctions`); and NM summed.
static AssembledFragment build_fragment(
    const std::vector<ParsedAlignment*>& r1_reads,
    const std::vector<ParsedAlignment*>& r2_reads)
{
    AssembledFragment frag;
    // Key: (ref_id << 32) | strand → the blocks on it
    std::unordered_map<int64_t, std::vector<std::pair<int32_t, int32_t>>> exon_dict;
    for (const bool r2 : {false, true}) {
        for (const auto* rec : r2 ? r2_reads : r1_reads) {
            const int32_t strand = (rec->is_reverse() != r2) ? STRAND_NEG : STRAND_POS;
            auto& blocks = exon_dict[(static_cast<int64_t>(rec->ref_id) << 32) | static_cast<uint32_t>(strand)];
            blocks.insert(blocks.end(), rec->exons.begin(), rec->exons.end());
            frag.nm += rec->nm;
        }
    }

    // Merge overlapping/adjacent exon blocks within each (ref, strand) group
    std::vector<ExonBlock> merged_exons;
    for (auto& [key, intervals] : exon_dict) {
        if (intervals.empty()) continue;
        int32_t ref_id = static_cast<int32_t>(key >> 32);
        int32_t strand = static_cast<int32_t>(key & 0xFFFFFFFF);

        std::sort(intervals.begin(), intervals.end());
        int32_t cur_start = intervals[0].first;
        int32_t cur_end = intervals[0].second;

        for (size_t i = 1; i < intervals.size(); i++) {
            if (intervals[i].first <= cur_end) {
                cur_end = std::max(cur_end, intervals[i].second);
            } else {
                merged_exons.push_back({ref_id, cur_start, cur_end, strand});
                cur_start = intervals[i].first;
                cur_end = intervals[i].second;
            }
        }
        merged_exons.push_back({ref_id, cur_start, cur_end, strand});
    }

    // Sort exon blocks for deterministic ordering
    std::sort(merged_exons.begin(), merged_exons.end(),
              [](const ExonBlock& a, const ExonBlock& b) {
                  if (a.ref_id != b.ref_id) return a.ref_id < b.ref_id;
                  if (a.start != b.start) return a.start < b.start;
                  if (a.end != b.end) return a.end < b.end;
                  return a.strand < b.strand;
              });

    frag.exons = std::move(merged_exons);
    observed_junctions(r1_reads, r2_reads, frag.introns);
    return frag;
}

// ================================================================
// A template's hits: every R1 alignment with every R2 alignment
// ================================================================

/// A template's candidate hits, each its R1 records and its R2 records; and, per record, the record it is written with
/// in the annotated BAM when it carries no hit of its own, an alias of an alignment another record represents (-1
/// otherwise).
struct TemplateHits {
    std::vector<std::pair<std::vector<ParsedAlignment*>, std::vector<ParsedAlignment*>>> hits;
    std::vector<int32_t> alias_of;
    int32_t n_nm_conflicting = 0;  // locations whose records disagree on NM
};

/// Whether two records carry one alignment: one mate, one BAM reference, one strand, one leftmost aligned position and
/// one CIGAR (`canonical_cigar`). STAR repeats a mate's record in every pair alignment that shares it.
static bool same_alignment(const ParsedAlignment& x, const ParsedAlignment& y)
{
    return x.is_read1() == y.is_read1() && x.tid == y.tid && x.is_reverse() == y.is_reverse() &&
           x.ref_start == y.ref_start && x.cigar == y.cigar;
}

/// Whether a record's read is hard-clipped at its 5' end: the CIGAR's first operation on a forward record, its last on
/// a reverse one. Those bases are not in the record, so the molecule's end is not either.
static bool five_prime_hard_clipped(const ParsedAlignment& r)
{
    if (r.cigar.empty()) return false;
    return bam_cigar_op(r.is_reverse() ? r.cigar.back() : r.cigar.front()) == BAM_CHARD_CLIP;
}

/// A template's candidate hits, found by geometry alone, for any aligner: the aligner's own pairing (STAR's HI) is not
/// read. The records are merged into mate locations, one per alignment (`same_alignment`): a primary record represents
/// its location, the others are its aliases. A location is chimeric, one piece of an alignment in pieces, when a record
/// of it carries SA or its 5' end is hard-clipped, or when that mate has a supplementary record and none of the mate's
/// records carries SA: which alignment the piece belongs to is then unknown, so every location of the mate is. A
/// supplementary record is no location and joins no hit; it is written with the same mate's primary location (an
/// alias), or its first location when none is primary. A record with no aligned base is a location too, which the
/// verdict passes through unjudged. Every R1 location with every R2 location on one BAM reference with opposite strands
/// is a candidate for the hit verdict (`hit_verdict`).
/// When one mate has no location, each location of the other is a one-mate hit (`HIT_MISSING_MATE`); when both have
/// locations and no pair is a candidate, the first location of each is offered, so the verdict counts the template's
/// reason once. On the STAR panel this finds the correct pair wherever a template's places hold one, where STAR's HI
/// groups miss it in 0.31 % of the templates with several places.
static TemplateHits template_hits(std::vector<ParsedAlignment>& records)
{
    struct Location {
        std::vector<int32_t> members;  // record indices, in record order
        int32_t rep = -1;
    };
    TemplateHits out;
    std::vector<Location> locs[2];
    std::vector<int32_t> supplementary[2];
    out.alias_of.assign(records.size(), -1);
    for (int32_t k = 0; k < static_cast<int32_t>(records.size()); ++k) {
        const ParsedAlignment& r = records[k];
        const int m = r.is_read1() ? 0 : 1;
        if (r.is_supplementary()) {
            supplementary[m].push_back(k);
            continue;
        }
        auto it = std::find_if(locs[m].begin(), locs[m].end(),
                               [&](const Location& l) { return same_alignment(records[l.members.front()], r); });
        if (it == locs[m].end()) locs[m].push_back({{k}, -1});
        else it->members.push_back(k);
    }
    for (auto& mate : locs) {
        for (auto& l : mate) {
            l.rep = l.members.front();
            for (int32_t k : l.members) {
                if (!records[k].is_secondary()) {
                    l.rep = k;
                    break;
                }
            }
            for (int32_t k : l.members)
                if (k != l.rep) out.alias_of[k] = l.rep;
            // Records of one alignment that disagree on NM are an input defect: counted, and the larger NM kept.
            int32_t lo = records[l.rep].nm, hi = lo;
            for (int32_t k : l.members) {
                lo = std::min(lo, records[k].nm);
                hi = std::max(hi, records[k].nm);
            }
            if (lo != hi) {
                ++out.n_nm_conflicting;
                records[l.rep].nm = hi;
            }
        }
    }
    for (int m = 0; m < 2; ++m) {
        bool sa = false;
        for (int32_t k : supplementary[m]) sa |= records[k].has_sa;
        for (const auto& l : locs[m])
            for (int32_t k : l.members) sa |= records[k].has_sa;
        const bool unnamed_chain = !supplementary[m].empty() && !sa;
        for (const auto& l : locs[m]) {
            ParsedAlignment& rep = records[l.rep];
            rep.chimeric = unnamed_chain || five_prime_hard_clipped(rep);
            for (int32_t k : l.members) rep.chimeric |= records[k].has_sa;
        }
        if (locs[m].empty()) continue;
        auto primary = std::find_if(locs[m].begin(), locs[m].end(),
                                    [&](const Location& l) { return !records[l.rep].is_secondary(); });
        const int32_t owner = (primary == locs[m].end() ? locs[m].front() : *primary).rep;
        for (int32_t k : supplementary[m]) out.alias_of[k] = owner;
    }
    auto reads = [&](const Location& l) { return std::vector<ParsedAlignment*>{&records[l.rep]}; };
    if (locs[0].empty() || locs[1].empty()) {
        for (const auto& l : locs[0]) out.hits.push_back({reads(l), {}});
        for (const auto& l : locs[1]) out.hits.push_back({{}, reads(l)});
        return out;
    }
    for (const auto& l1 : locs[0]) {
        const ParsedAlignment& a = records[l1.rep];
        for (const auto& l2 : locs[1]) {
            const ParsedAlignment& b = records[l2.rep];
            if (a.tid == b.tid && a.is_reverse() != b.is_reverse()) out.hits.push_back({reads(l1), reads(l2)});
        }
    }
    if (out.hits.empty()) out.hits.push_back({reads(locs[0].front()), reads(locs[1].front())});
    return out;
}

// ================================================================
// The hit verdict: two mates that can be one molecule, explained by one annotated transcript when spliced
// ================================================================

/// `HIT_USABLE` when a hit's two mates can be one molecule (`FragmentResolver::validate_pair`), or when the hit is not
/// a pair this rule judges: a mate with no aligned base, which stays as it was for the change that owns it. Otherwise
/// the reason they cannot; a hit with one mate gives no fragment length at all (`HIT_MISSING_MATE`), one on a
/// reference the index does not hold has no partition, transcript or length along one (`HIT_OFF_INDEX`), and one with
/// a mate aligned in pieces is a chimeric alignment, which Rigel does not quantify (`HIT_CHIMERIC`). Mates on two
/// references are two references, wherever they lie: references are told apart by the BAM's ids, as every reference
/// outside the index shares the index id -1. A transcript that witnesses a path must overlap both mates and hold every
/// junction the hit observes (``junctions``). ``f`` and ``r`` are set to the forward and reverse mates when the rule
/// judged them. Each side of a hit is one mate location's record (`template_hits`).
static int32_t pair_verdict(
    const std::pair<std::vector<ParsedAlignment*>, std::vector<ParsedAlignment*>>& hit,
    const std::vector<IntronBlock>& junctions,
    const FragmentResolver& ctx,
    ResolverScratch& scratch,
    const ParsedAlignment*& f,
    const ParsedAlignment*& r)
{
    f = r = nullptr;
    if (hit.first.empty() || hit.second.empty()) return HIT_MISSING_MATE;
    const ParsedAlignment* m1 = hit.first.front();
    const ParsedAlignment* m2 = hit.second.front();
    if (m1->ref_id < 0 || m2->ref_id < 0) return m1->tid != m2->tid ? HIT_REFERENCE : HIT_OFF_INDEX;
    if (m1->chimeric || m2->chimeric) return HIT_CHIMERIC;
    if (m1->exons.empty() || m2->exons.empty()) return HIT_USABLE;
    if (m1->ref_id != m2->ref_id) return HIT_REFERENCE;
    if (m1->is_reverse() == m2->is_reverse()) return HIT_STRAND;
    f = m1->is_reverse() ? m2 : m1;
    r = m1->is_reverse() ? m1 : m2;
    const PairMate pf{&f->exons, f->clip_lead, f->clip_trail};
    const PairMate pr{&r->exons, r->clip_lead, r->clip_trail};
    return ctx.validate_pair(f->ref_id, pf, pr, junctions, scratch);
}

/// Why a hit is not usable, or `HIT_USABLE`. Its two mates must be able to be one molecule (`pair_verdict`), and when
/// it is spliced, one annotated transcript must explain it: hold every junction it observes (`observed_junctions`) and
/// overlap both mates (`FragmentResolver::pair_holders`): the annotated-path rule. A hit the pair rule passes through
/// unjudged is held to its junctions alone. The molecule is judged first, so a hit the rule withholds always has a
/// fragment length.
static int32_t hit_verdict(
    const std::pair<std::vector<ParsedAlignment*>, std::vector<ParsedAlignment*>>& hit,
    const FragmentResolver& ctx,
    ResolverScratch& scratch)
{
    auto& junctions = scratch.hit_junctions;
    observed_junctions(hit.first, hit.second, junctions);
    const ParsedAlignment* f = nullptr;
    const ParsedAlignment* r = nullptr;
    const int32_t pair = pair_verdict(hit, junctions, ctx, scratch, f, r);
    if (pair != HIT_USABLE || junctions.empty()) return pair;
    if (!f) return ctx.junction_holders(junctions, scratch, scratch.pair_candidates);
    return ctx.pair_holders(junctions, f->exons.front().first, f->exons.back().second, r->exons.front().first,
                            r->exons.back().second, scratch, scratch.pair_candidates);
}

/// Every hit's verdict into ``verdicts`` (`hit_verdict`: HIT_USABLE or the reason), and why the template is withheld
/// when it has hits and none is usable: off the index when one of them lies on a reference the index does not hold,
/// else chimeric when one of them has a mate aligned in pieces, else an unannotated splice when the annotated-path
/// rule withheld one of them, else unpaired (WITHHELD_NONE otherwise). The scan and the annotation pass both judge a
/// template here, so the BAM the annotation pass writes says what the scan counted.
static int32_t judge_hits(const TemplateHits& th, const FragmentResolver& ctx, ResolverScratch& scratch,
                          std::vector<int32_t>& verdicts)
{
    verdicts.clear();
    bool usable = false, off_index = false, chimeric = false, unannotated_splice = false;
    for (const auto& hit : th.hits) {
        const int32_t v = hit_verdict(hit, ctx, scratch);
        verdicts.push_back(v);
        usable |= v == HIT_USABLE;
        off_index |= v == HIT_OFF_INDEX;
        chimeric |= v == HIT_CHIMERIC;
        unannotated_splice |= is_unannotated_splice(v);
    }
    if (usable || th.hits.empty()) return WITHHELD_NONE;
    if (off_index) return WITHHELD_OFF_INDEX;
    if (chimeric) return WITHHELD_CHIMERIC;
    return unannotated_splice ? WITHHELD_UNANNOTATED_SPLICE : WITHHELD_UNPAIRED;
}

/// A template's NH: the first record's above 1, else 1 (a missing NH reads 1).
static int32_t template_nh(const std::vector<ParsedAlignment>& records)
{
    for (const auto& r : records)
        if (r.nh > 1) return r.nh;
    return 1;
}

/// Whether a template is a multimapper: an NH above 1, or a secondary record (for an aligner that writes no NH).
static bool is_multimapper(const std::vector<ParsedAlignment>& records, int32_t nh)
{
    return nh > 1 ||
           std::any_of(records.begin(), records.end(), [](const ParsedAlignment& r) { return r.is_secondary(); });
}

/// One `OfferedFragment` marshalled out of Python, with the storage its spans point into.
///
/// ⭐ Shared by `Accumulator.deposit`, `Accumulator.length_under` and `AccumulatorSet.offer` so they cannot
/// disagree about how a hypothesis set crosses the ABI. ⚠ Reads attributes off whatever it is handed, so the
/// parity gate can pass the SAME `GapHypothesis` objects to the specification and to the binding; a binding
/// with its own tuple convention would be a second representation to keep in step.
///
/// ⛔ Not copyable: `offered` points into this object's own vectors, so a copy's `offered` would point into
/// the original's. A move hands the vectors' buffers over unchanged, so the moved-to `offered` stays valid.
struct MarshalledFragment {
    std::vector<IntronBlock> observed;
    std::vector<IntronBlock> implied;
    std::vector<rigel::accumulator::GapHypothesis> spans;
    rigel::accumulator::OfferedFragment offered{};

    MarshalledFragment(int64_t start, int64_t end, nb::iterable introns,
                       int32_t align_strand, int32_t sj_strand, nb::object hypotheses) {
        for (nb::handle item : introns) {
            auto pair = nb::cast<nb::tuple>(item);
            observed.push_back({0, nb::cast<int32_t>(pair[0]), nb::cast<int32_t>(pair[1]), 0});
        }
        // ⚠ ONE pass over `hypotheses`: a second pass over a one-shot iterable (a generator) sees nothing,
        // and an empty set deposits as the unspliced fragment alone. While `implied` grows, each span holds
        // an OFFSET into it, turned into a pointer only after the last hypothesis is read, so a
        // reallocation during the fill invalidates nothing. Pinned by
        // `test_a_GENERATOR_of_hypotheses_is_the_same_fragment_as_the_tuple`.
        for (nb::handle h : hypotheses) {
            const std::size_t i0 = implied.size();
            for (nb::handle item : nb::getattr(h, "introns")) {
                auto pair = nb::cast<nb::tuple>(item);
                implied.push_back({0, nb::cast<int32_t>(pair[0]), nb::cast<int32_t>(pair[1]), 0});
            }
            spans.push_back({nullptr, implied.size() - i0,
                             nb::cast<int32_t>(nb::getattr(h, "sj_strand"))});
            spans.back().introns = reinterpret_cast<const IntronBlock*>(i0);
        }
        for (auto& span : spans) {  // offsets -> pointers, now that nothing more will grow
            span.introns = implied.data() + reinterpret_cast<std::size_t>(span.introns);
        }

        offered.start = start;
        offered.end = end;
        offered.observed_introns = observed.data();
        offered.n_observed_introns = observed.size();
        offered.align_strand = align_strand;
        offered.sj_strand = sj_strand;
        offered.hypotheses = spans.data();
        offered.n_hypotheses = spans.size();
    }

    MarshalledFragment(const MarshalledFragment&) = delete;
    MarshalledFragment& operator=(const MarshalledFragment&) = delete;
    MarshalledFragment(MarshalledFragment&&) = default;
    MarshalledFragment& operator=(MarshalledFragment&&) = default;
};

// ================================================================
// the three non-array banks, as Python dicts
// ================================================================
//
// ⭐ Written ONCE and used by the bound `Accumulator` and `AccumulatorSet` (the parity surface) and by
// `build_result` (the payload). The specification's key strings live here and nowhere else on this side of
// the ABI, so the exports cannot disagree about a name — which is the whole reason the payload, the
// reference and the parity gate share one vocabulary.

/// The deposit counters: the QC denominators. Every conservation statement downstream has to be able to
/// name what it excluded, and none of these is derivable after the fact.
static nb::dict qc_dict(const rigel::accumulator::DepositCounters& c) {
    nb::dict out;
    out["deposited"]                    = c.deposited;
    out["dropped_too_long"]             = c.dropped_too_long;
    out["dropped_empty"]                = c.dropped_empty;
    out["dropped_strand_undefined"]     = c.dropped_strand_undefined;
    out["deferred_undetermined_gap"]    = c.deferred_undetermined_gap;
    out["deferred_multiple_placements"] = c.deferred_multiple_placements;
    return out;
}

/// The umbrella census, keyed by the specification's `GapResolution` values.
///
/// ⛔ There is no `gap_resolved_unspliced` key. See `GapCensus`: the unspliced hypothesis is always the
/// longest, so it can never be the sole survivor, and the class it would name is unreachable.
static nb::dict gap_census_dict(const rigel::accumulator::GapCensus& census) {
    nb::dict out;
    out["gap_resolved_spliced"]       = census.resolved_spliced;
    out["gap_deferred_rna_or_gdna"]   = census.deferred_rna_or_gdna;
    out["gap_deferred_which_introns"] = census.deferred_which_introns;
    out["gap_deferred_both"]          = census.deferred_both;
    return out;
}

/// The deferred queue as the CSR `Tally.deferred_arrays()` specifies — same keys, same dtype (int64
/// throughout), same canonical order.
///
/// ⚠ Every array is COPIED out. The caller has just canonicalised, which may have reseated the underlying
/// vectors, and a later export would reseat them again — a numpy view over them would dangle and read as
/// plausible coordinates. The queue is 1–3.5 % of a library, so the copy is bounded and one-off.
static nb::dict deferred_dict(const rigel::accumulator::DeferredFragments& deferred) {
    nb::dict out;
    const auto put = [&out](const char* name, const std::vector<int64_t>& v) {
        out[name] = rigel::vec_to_ndarray(std::vector<int64_t>(v));
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
    return out;
}

// ================================================================
// BamScanner — main scanning class
// ================================================================

class BamScanner {
public:
    // Set by the one constructor, from the caller's scan config.
    FragmentResolver* ctx_;
    SJTagMode sj_tag_mode_;
    bool skip_duplicates_;
    bool include_multimap_;

    // Mapping from htslib tid → internal ref_id (FragmentResolver's numbering)
    std::vector<int32_t> tid_to_ref_id_;

    // Results
    BamScanStats stats_;
    StrandObservations strand_obs_;

    // The accumulator's region partition, installed by set_regions() and the sj CSR by
    // set_sj(). One Accumulator per BAM reference; deposits go into a per-worker copy during
    // scan() and are merged back here.
    //
    // ⚠ These are MEMBERS because the per-worker sets are built from them inside scan(). Anything
    // installed after that point is invisible to the workers, which is a silent half-empty tally.
    std::vector<int64_t> region_bounds_;
    std::vector<int64_t> ref_region_bound_offsets_;
    std::vector<uint8_t> region_types_;      // coarse region type, ref-major, one per region
    int                  max_length_ = 0;  // the fragment-length limit AND the pool-histogram width
    std::vector<int64_t> sj_offsets_;       // n_region_bounds + 1, CSR over the flat region_bound axis
    std::vector<int64_t> sj_boundary_right_;  // flat region_bound index of each sj's high end
    std::vector<int8_t>  sj_strand_;        // each sj's ANNOTATED strand
    bool                 sj_set_ = false;
    std::unique_ptr<rigel::accumulator::AccumulatorSet> acc_set_;

    BamScanner(FragmentResolver& ctx,
               const std::string& sj_tag_spec,
               bool skip_duplicates,
               bool include_multimap)
        : ctx_(&ctx),
          sj_tag_mode_(parse_sj_tag_spec(sj_tag_spec)),
          skip_duplicates_(skip_duplicates),
          include_multimap_(include_multimap)
    {}

    // ----------------------------------------------------------------
    // The accumulator's partition
    // ----------------------------------------------------------------
    //
    // Install the region partition the accumulator deposits onto. Must be called BEFORE scan(); without it
    // the calibration payload comes back empty. `set_sj` is a SECOND call, because this one
    // refuses to run twice.
    //
    // Inputs:
    //   region_bounds:   flat int64[n_region_bounds_total], the concatenated sorted region_bound positions of every
    //                    reference. A reference contributing c region_bounds owns c-1 regions and c-2 boundaries.
    //   ref_region_bound_offsets: int64[n_refs + 1] offsets into region_bounds. A reference with
    //                    offsets[f+1] == offsets[f] has no partition and accepts no deposits.
    //   n_refs:          number of references (must match ctx_->ref_to_id_).
    //   region_types:      uint8[n_regions_total], the coarse region type, ref-major; it types the length pools.
    //   max_length:      the fragment-length limit applied to L, and the pool-histogram width.
    void set_regions(
        nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> region_bounds,
        nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> ref_region_bound_offsets,
        int32_t n_refs,
        nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig> region_types,
        int max_length)
    {
        if (n_refs < 0) {
            throw std::invalid_argument(
                "set_regions: n_refs must be non-negative");
        }
        const std::size_t n_off = ref_region_bound_offsets.shape(0);
        if (n_off != static_cast<std::size_t>(n_refs) + 1) {
            throw std::invalid_argument(
                "set_regions: ref_region_bound_offsets must have length n_refs + 1");
        }
        if (acc_set_) {
            throw std::runtime_error(
                "set_regions: regions already set on this BamScanner; "
                "create a new instance");
        }
        // ⚠ Validated HERE, at the boundary, because `BamScanConfig.max_frag_length` has no CLI flag and
        // no validation anywhere else. It gates L *and* sizes the pool histograms, so at 0 every real
        // fragment is dropped as too long and the whole tally is silently empty.
        if (max_length < 1) {
            throw std::invalid_argument(
                "set_regions: max_length must be >= 1, got " + std::to_string(max_length) +
                ". It is the fragment-length limit applied to L as well as the pool-histogram width, so "
                "at 0 every real fragment is dropped as too long and the tally is silently empty.");
        }
        region_bounds_.assign(
            region_bounds.data(),
            region_bounds.data() + region_bounds.shape(0));
        ref_region_bound_offsets_.assign(
            ref_region_bound_offsets.data(),
            ref_region_bound_offsets.data() + n_off);
        region_types_.assign(
            region_types.data(),
            region_types.data() + region_types.shape(0));
        max_length_ = max_length;
        acc_set_ = make_accumulator_set();
    }

    // ----------------------------------------------------------------
    // The accumulator's sj boundaries
    // ----------------------------------------------------------------
    //
    // A SECOND method rather than two more arguments to set_regions, because that one throws if called
    // twice. Takes the flat CSR `build_sj_arrays` emits, keyed by the flat region_bound index;
    // AccumulatorSet slices it per reference.
    //
    // ⚠ Required even when there are none: pass empty arrays to say "this annotation has no sj".
    // scan() refuses to run without it, because a missing table is not an error anywhere — every observed
    // intron simply reads as unannotated, and 404,168 annotated sj become zero silently.
    void set_sj(
        nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> offsets,
        nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> boundary_right,
        nb::ndarray<const int8_t, nb::ndim<1>, nb::c_contig> sj_strand)
    {
        if (!acc_set_) {
            throw std::runtime_error(
                "set_sj: call set_regions first — the sj CSR is keyed by region_bound index, so "
                "there is nothing to key it against until the partition is installed");
        }
        if (boundary_right.shape(0) != sj_strand.shape(0)) {
            throw std::invalid_argument(
                "set_sj: boundary_right has " + std::to_string(boundary_right.shape(0)) +
                " entries but sj_strand has " + std::to_string(sj_strand.shape(0)));
        }
        sj_offsets_.assign(offsets.data(), offsets.data() + offsets.shape(0));
        sj_boundary_right_.assign(
            boundary_right.data(), boundary_right.data() + boundary_right.shape(0));
        sj_strand_.assign(sj_strand.data(), sj_strand.data() + sj_strand.shape(0));
        sj_set_ = true;
        install_sj(*acc_set_);
    }

    // Build a per-worker set over the same partition, locally writable so workers never contend.
    std::unique_ptr<rigel::accumulator::AccumulatorSet> make_accumulator_set() const {
        auto set = std::make_unique<rigel::accumulator::AccumulatorSet>(
            region_bounds_.data(),
            region_bounds_.size(),
            ref_region_bound_offsets_.data(),
            ref_region_bound_offsets_.empty() ? 0 : ref_region_bound_offsets_.size() - 1,
            region_types_.empty() ? nullptr : region_types_.data(),
            region_types_.size(),
            max_length_);
        install_sj(*set);
        return set;
    }

    void install_sj(rigel::accumulator::AccumulatorSet& set) const {
        if (!sj_set_) return;
        set.set_sj(sj_offsets_.data(),
                          sj_offsets_.size(),
                          sj_boundary_right_.data(),
                          sj_strand_.data(),
                          sj_boundary_right_.size(),
                          ref_region_bound_offsets_.data());
    }

    // ----------------------------------------------------------------
    // Main scan entry point — streaming chunk architecture
    // ----------------------------------------------------------------
    //
    // Three-thread architecture:
    //   Reader thread (background): BAM I/O → input_queue
    //   Worker threads (background): input_queue → accumulate → output_queue
    //   Main thread: output_queue → finalize → chunk_callback
    //
    // The main thread holds the GIL while calling the Python callback
    // and releases it while waiting on the output queue.

    nb::dict scan(const std::string& bam_path,
                  nb::callable chunk_callback,
                  int64_t chunk_size,
                  int n_workers,
                  int n_decomp_threads,
                  int qname_batch_size)
    {
        if (n_workers < 1) n_workers = 1;
        if (chunk_size < 1) chunk_size = 1;
        if (qname_batch_size < 1) qname_batch_size = 1;

        // ⚠ Explicit rather than defaulted, because the failure is invisible: with no sj table every
        // observed intron reads as unannotated, so the spliced banks and every sj boundary stay at zero
        // and the tally still looks well-formed. Pass empty arrays to declare an annotation with none.
        if (acc_set_ && !sj_set_) {
            throw std::runtime_error(
                "scan: set_regions was called but set_sj was not. The sj table is not "
                "optional — without it every observed intron reads as unannotated, so the sj boundaries "
                "and the spliced banks are silently empty. Call set_sj with empty arrays if this "
                "annotation genuinely has no sj.");
        }

        // Two queues: input (SPMC) and output (MPSC)
        BoundedQueue<QnameBatch> input_queue(
            static_cast<size_t>(n_workers * INPUT_BATCHES_PER_WORKER));
        BoundedQueue<FragmentAccumulator> output_queue(
            static_cast<size_t>(n_workers * OUTPUT_CHUNKS_PER_WORKER));

        // Per-worker state (local to scan — not a class member)
        std::vector<std::unique_ptr<WorkerState>> worker_states;
        worker_states.reserve(n_workers);
        for (int i = 0; i < n_workers; i++) {
            int32_t n_transcripts = ctx_->n_transcripts_;
            auto ws = std::make_unique<WorkerState>(n_transcripts);
            reserve_chunk(ws->accumulator, chunk_size);
            // Per-worker accumulator: the same partition and the same sj CSR as the shared
            // template, locally writable so workers don't contend.
            if (acc_set_) {
                ws->acc_set = make_accumulator_set();
            }
            worker_states.push_back(std::move(ws));
        }

        // Capture read-only config
        bool include_multimap = include_multimap_;
        FragmentResolver* ctx = ctx_;

        // ---- Launch worker threads ----
        // Workers pop from input_queue, process groups, accumulate, and
        // push full chunks to output_queue. A worker that throws keeps the
        // first worker exception and aborts both queues.
        std::exception_ptr worker_exception;
        std::mutex worker_exception_mutex;
        std::vector<std::thread> workers;
        workers.reserve(n_workers);
        for (int i = 0; i < n_workers; i++) {
            workers.emplace_back([&input_queue, &output_queue,
                                  &worker_states, &worker_exception,
                                  &worker_exception_mutex, ctx, i,
                                  include_multimap, chunk_size]()
            {
                try {
                    WorkerState& ws = *worker_states[i];
                    QnameBatch batch;
                    bool output_aborted = false;
                    while (input_queue.pop(batch)) {
                        for (auto& group : batch.groups) {
                            process_qname_group_threaded(
                                group, *ctx, ws, include_multimap);
                            // Emit a chunk when accumulator reaches threshold
                            if (ws.accumulator.get_size() >= chunk_size) {
                                if (!output_queue.push(std::move(ws.accumulator))) {
                                    output_aborted = true;
                                    break;  // aborted
                                }
                                ws.accumulator = FragmentAccumulator();
                                reserve_chunk(ws.accumulator, chunk_size);
                            }
                        }
                        batch.groups.clear();
                        if (output_aborted) break;
                    }
                    // Flush remaining fragments
                    if (!output_aborted && ws.accumulator.get_size() > 0) {
                        output_queue.push(std::move(ws.accumulator));
                    }
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lock(worker_exception_mutex);
                        if (!worker_exception) {
                            worker_exception = std::current_exception();
                        }
                    }
                    input_queue.abort();
                    output_queue.abort();
                }
            });
        }

        // ---- Launch reader thread (BAM I/O → input_queue) ----
        std::exception_ptr reader_exception;
        std::thread reader_thread([&]() {
            try {
                htsFile* fp = hts_open(bam_path.c_str(), "rb");
                if (!fp) {
                    throw std::runtime_error(
                        "Failed to open BAM file: " + bam_path);
                }
                if (n_decomp_threads > 0) {
                    hts_set_threads(fp, n_decomp_threads);
                }
                bam_hdr_t* hdr = sam_hdr_read(fp);
                if (!hdr) {
                    hts_close(fp);
                    throw std::runtime_error(
                        "Failed to read BAM header: " + bam_path);
                }
                build_tid_mapping(hdr, ctx->ref_to_id_, tid_to_ref_id_);

                bam1_t* b = bam_init1();
                if (!b) {
                    bam_hdr_destroy(hdr);
                    hts_close(fp);
                    throw std::runtime_error("Failed to allocate bam1_t");
                }

                QnameGroup current_group;
                QnameBatch current_batch;
                current_batch.groups.reserve(
                    static_cast<size_t>(qname_batch_size));
                std::string current_qname;
                int64_t frag_id = 0;
                bool input_aborted = false;

                auto push_group = [&](QnameGroup&& group) -> bool {
                    current_batch.groups.push_back(std::move(group));
                    if (static_cast<int>(current_batch.groups.size()) < qname_batch_size) {
                        return true;
                    }
                    if (!input_queue.push(std::move(current_batch))) {
                        return false;
                    }
                    current_batch = QnameBatch{};
                    current_batch.groups.reserve(
                        static_cast<size_t>(qname_batch_size));
                    return true;
                };

                auto flush_batch = [&]() -> bool {
                    if (current_batch.groups.empty()) return true;
                    return input_queue.push(std::move(current_batch));
                };

                int read_ret = 0;
                while ((read_ret = sam_read1(fp, hdr, b)) >= 0) {
                    stats_.total++;

                    uint16_t flag = b->core.flag;
                    if (flag & BAM_FQCFAIL) { stats_.qc_fail++; continue; }
                    if (flag & BAM_FUNMAP) { stats_.unmapped++; continue; }
                    if (flag & BAM_FDUP) {
                        stats_.duplicate++;
                        if (skip_duplicates_) continue;
                    }
                    if (!(flag & BAM_FPAIRED)) {
                        bam_destroy1(b);
                        bam_hdr_destroy(hdr);
                        hts_close(fp);
                        throw std::runtime_error(
                            "Input BAM must be paired-end, "
                            "but found unpaired read");
                    }
                    if (flag & BAM_FSECONDARY) stats_.secondary++;
                    if (flag & BAM_FSUPPLEMENTARY) stats_.supplementary++;

                    const char* qname = bam_get_qname(b);
                    if (!current_group.records.empty() &&
                        current_qname != qname) {
                        current_group.frag_id = frag_id++;
                        if (!push_group(std::move(current_group))) {
                            input_aborted = true;
                            break;  // aborted
                        }
                        current_group = QnameGroup{};
                    }
                    current_qname = qname;
                    current_group.records.push_back(
                        parse_bam_record(b, tid_to_ref_id_, sj_tag_mode_, ctx_, &stats_));
                }

                if (read_ret < -1) {
                    bam_destroy1(b);
                    bam_hdr_destroy(hdr);
                    hts_close(fp);
                    throw read_error(bam_path, read_ret, stats_.total);
                }

                // Flush last group
                if (!input_aborted && !current_group.records.empty()) {
                    current_group.frag_id = frag_id++;
                    input_aborted = !push_group(std::move(current_group));
                }
                if (!input_aborted) {
                    flush_batch();
                }

                bam_destroy1(b);
                bam_hdr_destroy(hdr);
                hts_close(fp);

                // Signal workers: no more input
                input_queue.close();

                // Wait for workers to finish, then close output queue
                for (auto& w : workers) w.join();

                // Merge worker results (stats, strand, accumulator)
                for (auto& ws_ptr : worker_states) {
                    WorkerState& ws = *ws_ptr;
                    stats_.merge_from(ws.stats);
                    merge_strand_obs(strand_obs_, ws.strand_obs);
                    if (acc_set_ && ws.acc_set) {
                        acc_set_->merge_from(*ws.acc_set);
                    }
                }

                output_queue.close();
            } catch (...) {
                reader_exception = std::current_exception();
                input_queue.abort();
                // Still need to join workers before closing output
                for (auto& w : workers) {
                    if (w.joinable()) w.join();
                }
                output_queue.abort();
            }
        });

        // ---- Main thread: drain output queue with GIL cycling ----
        std::exception_ptr main_exception;
        try {
            FragmentAccumulator acc;
            // Release GIL while waiting for chunks
            {
                nb::gil_scoped_release release;
                while (output_queue.pop(acc)) {
                    // Acquire GIL for finalize + Python callback
                    nb::gil_scoped_acquire acquire;
                    nb::dict chunk = acc.finalize();
                    chunk_callback(chunk);
                    // GIL released again at loop top
                }
            }
        } catch (...) {
            main_exception = std::current_exception();
            // Abort both queues to unblock reader/workers
            input_queue.abort();
            output_queue.abort();
        }

        // Wait for reader thread to finish
        reader_thread.join();

        // Propagate exceptions: reader, then worker, then main thread. The reader's comes first because a BAM
        // it cannot read is the root cause of whatever followed. An abort makes no other thread throw (an
        // aborted push or pop returns false), so no stored exception is an echo of another's.
        if (reader_exception) {
            std::rethrow_exception(reader_exception);
        }
        if (worker_exception) {
            std::rethrow_exception(worker_exception);
        }
        if (main_exception) {
            std::rethrow_exception(main_exception);
        }

        return build_result();
    }

private:

    // ----------------------------------------------------------------
    // Process one qname group (worker thread, pre-parsed records)
    // ----------------------------------------------------------------
    // Static so it doesn't capture 'this' (shares no mutable BamScanner
    // state — all mutable state is in WorkerState).

    static void process_qname_group_threaded(
        QnameGroup& group,
        FragmentResolver& ctx,
        WorkerState& ws,
        bool include_multimap)
    {
        if (group.records.empty()) return;

        auto& records = group.records;
        int64_t frag_id = group.frag_id;
        BamScanStats& stats = ws.stats;
        FragmentAccumulator& accumulator = ws.accumulator;
        StrandObservations& strand_obs = ws.strand_obs;
        ResolverScratch& scratch = ws.scratch;

        // ── the deposit adapter ───────────────────────────────────────────────────────────────────────
        //
        // Turn one assembled fragment into the accumulator's `OfferedFragment` and deposit it. The
        // accumulator owns the whole deposit rule; this function's only job is to say what the fragment IS
        // — its extent on one reference, the introns cut out of it, and the two independent strands.
        //
        // ⭐ THE TWO STRANDS ARE INDEPENDENT.
        //   align_strand  where the read ALIGNED. Every read has one. It selects the array column.
        //   sj_strand     the splice junction's strand. Spliced reads only. It resolves an intron against
        //                 the annotation, and nothing else.
        // Sense/antisense is DERIVED by a consumer from the fragment strand and the sj's own strand.
        //
        // ⚠ No strand gate here. The deposit rejects an undefined `align_strand` itself and COUNTS it.
        const auto collect_placement =
            [&ws](const AssembledFragment& f, const RawResolveResult& cr) {
                const int32_t st = cr.splice_type;

                // ── the splice census ─────────────────────────────────────────────────────────────
                //
                // ⭐ FIRST, and before any gate, because this is the count of what the scanner SAW:
                // one observation per hit reaching this adapter.
                //
                // ⚠ It is deliberately NOT gated on `ws.acc_set`: what the scanner classified does
                // not depend on whether an accumulator happened to be installed.
                if (st >= 0 && static_cast<size_t>(st) < ws.stats.splice_census.size())
                    ws.stats.splice_census[st]++;

                // Hold artifact splices (blacklisted CIGAR-N) out of the accumulator entirely — no
                // deposit, no length pool. Their true span is unrecoverable: a blacklisted sj may be
                // a real-but-rejected sj OR a wholly incorrect alignment, and the span would be
                // derived from that suspect alignment, so any reconstruction injects a false assumption.
                //
                // ⚠ Counted by the census above and NOT by `n_deposit_not_offered`: this is a
                // positive decision about the data, not a failure to represent it.
                if (st == SPLICE_ARTIFACT) return;

                if (!ws.acc_set || f.exons.empty()) { ws.stats.n_deposit_not_offered++; return; }
                const int32_t ref_id = f.exons.front().ref_id;
                if (ref_id < 0 ||
                    static_cast<std::size_t>(ref_id) >= ws.acc_set->n_refs()) {
                    ws.stats.n_deposit_not_offered++;
                    return;
                }

                // ── the introns the CIGAR actually stated ─────────────────────────────────────
                //
                // ⭐ These are cut under EVERY hypothesis, so they are NOT part of the hypothesis set.
                // ⚠ The union with a hypothesis's implied introns happens inside `deposit`, where L is
                // defined — the adapter used to do it here, and doing it in two places is how one
                // quantity ends up with two formulas.
                //
                // ⛔ Restricting to `ref_id` is DEFINITIONAL, not defensive: `OfferedFragment` means the
                // introns on the reference being deposited, and `Accumulator::deposit` normalises them by
                // coordinate alone — it never looks at `ref_id` — so an intron from another reference would
                // be cut out of this reference's path. The exon-block check below cannot catch one: a record
                // whose CIGAR has an N and no M, D, = or X op adds an intron and no exon block.
                //
                // ⚠ The de-duplication is NOT to prevent double-crediting a sj -- `deposit` normalises a
                // duplicate away either way. `parse_bam_record` reads XS/ts once per RECORD and
                // `observed_junctions` keys a hit's junctions on (ref, start, end, strand), so a pair where one
                // mate carries the tag yields the same intron twice; dropping the second keeps it out of a
                // held record, which stores its observed introns as offered and sorts on them. The junctions
                // are ordered, so consecutive comparison is enough; ⛔ they must NOT be re-keyed or merged in
                // `observed_junctions`, which the resolver shares: a merged junction may be one no transcript
                // holds, and the resolver refuses a fragment with one.
                auto& observed = ws.arena_observed;
                const std::size_t obs0 = observed.size();
                for (const auto& intron : f.introns) {
                    if (intron.ref_id != ref_id) continue;
                    if (observed.size() > obs0 && observed.back().start == intron.start &&
                        observed.back().end == intron.end) {
                        continue;
                    }
                    observed.push_back(intron);
                }

                // ── every explanation of the unsequenced gaps ─────────────────────────────────────
                //
                // ⭐ The resolver enumerated them; the ACCUMULATOR arbitrates. This only re-presents the
                // flat CSR as the span-of-spans the deposit interface takes. There is always at least
                // one — the unspliced hypothesis — so `deposit` never receives an empty set.
                const std::size_t hyp0 = ws.pending_hypotheses.size();
                for (std::int32_t h = 0; h < cr.n_gap_hypotheses(); ++h) {
                    const std::int32_t i0 = cr.gap_intron_offsets[h], i1 = cr.gap_intron_offsets[h + 1];
                    const std::size_t imp0 = ws.arena_implied.size();
                    ws.arena_implied.insert(ws.arena_implied.end(), cr.gap_introns.begin() + i0,
                                            cr.gap_introns.begin() + i1);
                    ws.pending_hypotheses.push_back({imp0, ws.arena_implied.size(), cr.gap_sj_strand[h]});
                }

                // The molecule's extent on this reference: leftmost block start to rightmost block end,
                // MATE GAP INCLUDED, because the gap is part of the molecule and must count toward L.
                //
                // ⛔ A fragment with blocks on MORE THAN ONE REFERENCE deposits nothing. It is not one
                // molecule, and an `OfferedFragment` cannot express it — it carries one
                // extent on one region_bound axis.
                std::int64_t start = 0, end = 0;
                bool any = false;
                for (const auto& block : f.exons) {
                    if (block.ref_id != ref_id) {
                        ws.stats.n_deposit_not_offered++;
                        ws.arena_observed.resize(obs0);
                        ws.pending_hypotheses.resize(hyp0);
                        return;
                    }
                    if (!any) {
                        start = block.start;
                        end = block.end;
                        any = true;
                    } else {
                        start = std::min<std::int64_t>(start, block.start);
                        end = std::max<std::int64_t>(end, block.end);
                    }
                }
                if (!any) {
                    ws.stats.n_deposit_not_offered++;
                    ws.arena_observed.resize(obs0);
                    ws.pending_hypotheses.resize(hyp0);
                    return;
                }

                // ⚠ The OBSERVED motif strand, straight from `cr`. A hypothesis carries the strand its
                // supporting transcripts imply, and `deposit` falls back to that only when nothing was
                // sequenced. An observed motif is evidence; an implied strand is an inference from the
                // annotation.
                ws.pending.push_back({ref_id, start, end, cr.align_strand, cr.sj_strand,
                                      obs0, observed.size(), hyp0, ws.pending_hypotheses.size()});
            };

        // ── the offer: every placement of this fragment, at once ──────────────────────────────────────
        const auto offer_placements = [&ws]() {
            if (!ws.acc_set || ws.pending.empty()) return;
            // Offsets → pointers, now that nothing more will grow.
            ws.arena_spans.clear();
            ws.arena_spans.reserve(ws.pending_hypotheses.size());
            for (const auto& h : ws.pending_hypotheses) {
                ws.arena_spans.push_back({ws.arena_implied.data() + h.imp0, h.imp1 - h.imp0, h.sj_strand});
            }
            ws.placements.clear();
            ws.placements.reserve(ws.pending.size());
            for (const auto& p : ws.pending) {
                rigel::accumulator::Placement pl;
                pl.ref_id = p.ref_id;
                pl.offered.start = p.start;
                pl.offered.end = p.end;
                pl.offered.observed_introns = ws.arena_observed.data() + p.obs0;
                pl.offered.n_observed_introns = p.obs1 - p.obs0;
                pl.offered.align_strand = p.align_strand;
                pl.offered.sj_strand = p.sj_strand;
                pl.offered.hypotheses = ws.arena_spans.data() + p.hyp0;
                pl.offered.n_hypotheses = p.hyp1 - p.hyp0;
                ws.placements.push_back(pl);
            }
            ws.acc_set->offer(ws.placements.data(), ws.placements.size(), ws.offer_scratch);
        };

        // Per-worker state refs
        stats.n_read_names++;
        ws.pending.clear();
        ws.pending_hypotheses.clear();
        ws.arena_observed.clear();
        ws.arena_implied.clear();

        const int32_t nh = template_nh(records);
        if (is_multimapper(records, nh)) {
            stats.multimapping++;
            if (!include_multimap) return;
        } else {
            stats.unique++;
        }

        auto hit_result = template_hits(records);
        stats.n_nm_conflicting += hit_result.n_nm_conflicting;

        // NH counts a template's pair alignments, and an aligner writes a mate's record in each; a mate reported in
        // fewer records means the aligner truncated its output (counted). The mate counters, once per template, read
        // its primary records.
        {
            int32_t n_r1 = 0, n_r2 = 0;
            const ParsedAlignment* p1 = nullptr;
            const ParsedAlignment* p2 = nullptr;
            for (const auto& r : records) {
                if (r.is_supplementary()) continue;
                ++(r.is_read1() ? n_r1 : n_r2);
                if (r.is_secondary()) continue;
                if (r.is_read1() && !p1) p1 = &r;
                else if (!r.is_read1() && !p2) p2 = &r;
            }
            if (nh > std::max(n_r1, n_r2)) stats.n_nh_above_reported++;
            if (!p1 || !p2) stats.mate_unmapped++;
            else if (p1->is_proper_pair()) stats.proper_pair++;
            else stats.improper_pair++;
        }

        auto& all_hits = hit_result.hits;

        // ⭐ A hit is usable when its two mates can be one molecule and one annotated transcript holds every junction
        // it observes (`hit_verdict`). A hit that is not is dropped, counted by its reason; a template left with none
        // is withheld, counted by its outcome (`judge_hits`), and is no fragment.
        auto& verdicts = scratch.hit_verdicts;
        const int32_t withheld = judge_hits(hit_result, ctx, scratch, verdicts);
        for (int32_t v : verdicts)
            if (v != HIT_USABLE) stats.hits_rejected[v]++;
        if (withheld != WITHHELD_NONE) {
            stats.withheld[withheld]++;
            return;
        }
        size_t kept = 0;
        for (size_t h = 0; h < all_hits.size(); ++h) {
            if (verdicts[h] != HIT_USABLE) continue;
            if (kept != h) all_hits[kept] = std::move(all_hits[h]);
            ++kept;
        }
        all_hits.resize(kept);

        // ⭐ The template's multiplicity is its valid hits, not the aligner's NH: a hit dropped above is no place the
        // template could be from.
        int32_t num_hits = static_cast<int32_t>(all_hits.size());
        bool is_unique_mapper = (num_hits == 1);
        int32_t n_buffered_mm = 0;

        // Track per-molecule intergenic status across all hits.
        // For multimappers: count intergenic only if NO hit resolves,
        // and count it exactly once (not per-hit).
        bool any_hit_resolved = false;
        bool any_non_chimeric_resolved = false;
        bool any_unresolved = false;

        // Track chimera per-fragment (not per-hit).
        // A fragment is chimeric only if it resolves but ALL
        // resolved hits are chimeric.
        bool any_hit_chimeric = false;
        int32_t worst_chimera_type = CHIMERA_NONE;

        // Count ONE fragment per physical molecule, not per hit.
        stats.n_fragments++;

        for (const auto& [r1_reads, r2_reads] : all_hits) {
            AssembledFragment frag = build_fragment(r1_reads, r2_reads);

            if (frag.exons.empty()) continue;

            // aggregate per-record blacklist counts into a
            // per-fragment total so the resolver can promote
            // SPLICE_UNSPLICED -> SPLICE_ARTIFACT.
            int32_t frag_n_sj_blacklisted = 0;
            for (const auto* r : r1_reads)
                frag_n_sj_blacklisted += r->n_sj_blacklisted;
            for (const auto* r : r2_reads)
                frag_n_sj_blacklisted += r->n_sj_blacklisted;

            RawResolveResult cr;
            cr.n_sj_blacklisted = frag_n_sj_blacklisted;
            ctx._resolve_core(
                frag.exons, frag.introns,
                frag.genomic_footprint(), cr, scratch);

            // Empty t_inds: an intergenic fragment. It is unspliced: a spliced hit the scan keeps has a transcript
            // holding its junctions (`hit_verdict`). Counted once all hits are seen, to count the template once.
            if (cr.t_inds.empty()) {
                any_unresolved = true;

                // unique-mapper zero-candidate fragments are
                // appended to the buffer so that calibration's
                // geometric categorization can label them INTERGENIC.
                //
                // We do NOT set any_hit_resolved here -- the
                // n_intergenic telemetry counter
                // (driven by !any_hit_resolved below) must keep firing
                // so the PipelineStats accountability sum stays
                // consistent.  The C++ scorer skips empty-t_inds
                // fragments silently (n_cand <= 0 early skip) so they
                // don't increment stat_gated either.
                if (is_unique_mapper) {
                    ResolvedFragment ig_result = ResolvedFragment::from_core(cr);
                    ig_result.num_hits = num_hits;
                    ig_result.nm = frag.nm;
                    accumulator.append(ig_result, frag_id);
                }

                // An intergenic hit is a placement for calibration whatever the fragment's other hits:
                // gDNA sits anywhere, and a repeat's gDNA often sits exactly here (the EM never sees a
                // multimapper's intergenic hits, and that is unchanged).
                collect_placement(frag, cr);

                continue;
            }

            any_hit_resolved = true;

            ResolvedFragment result = ResolvedFragment::from_core(cr);
            result.num_hits = num_hits;
            result.nm = frag.nm;

            if (result.get_is_chimeric()) {
                any_hit_chimeric = true;
                if (result.chimera_type > worst_chimera_type)
                    worst_chimera_type = result.chimera_type;
                accumulator.append(result, frag_id);
                continue;
            }

            any_non_chimeric_resolved = true;

            // For multimappers, count exon/splice/strand stats only
            // on the first resolved non-chimeric hit to avoid inflation.
            bool count_stats = is_unique_mapper || (n_buffered_mm == 0);

            if (count_stats) {
                stats.n_with_exon++;
                if (result.splice_type == SPLICE_SPLICED_ANNOT) {
                    stats.n_with_annotated_sj++;
                }

                if (result.get_is_same_strand()) {
                    stats.n_same_strand++;
                } else {
                    stats.n_ambig_strand++;
                }
            }

            if (is_unique_mapper) {
                if (result.get_is_strand_qualified()) {
                    // ⭐ ONE observation per strand-qualified fragment, keyed by
                    // SJ. `sj_strand` is qualified to POS/NEG here and is the
                    // motif strand of every annotated sj the fragment crosses,
                    // so `sense` is exactly the align==sj bit the 2×2 is built from
                    // and `sj_key_*` is necessarily set.
                    //
                    // The table's marginal IS the 2x2 the strand fit needs, so no
                    // per-fragment label vector is pushed alongside it.
                    // `stats.n_strand_trained` below counts the fragments, and Python
                    // asserts it against the table's total depth — the invariant that
                    // one fragment credits one sj.
                    SJKey sj_key{result.sj_key_ref, result.sj_key_start,
                                 result.sj_key_end, result.sj_strand};
                    auto& counts = strand_obs.sj_strand_table[sj_key];
                    if (result.align_strand == result.sj_strand) counts.n_sense++;
                    else                                         counts.n_antisense++;
                    stats.n_strand_trained++;
                } else if (result.splice_type != SPLICE_SPLICED_ANNOT) {
                    stats.n_strand_skipped_no_sj++;
                } else if (!result.get_is_same_strand()) {
                    stats.n_strand_skipped_ambig_strand++;
                } else {
                    stats.n_strand_skipped_ambiguous++;
                }

                if (result.get_is_same_strand() &&
                    (result.align_strand == STRAND_POS ||
                     result.align_strand == STRAND_NEG)) {
                    int32_t t_idx = result.get_first_t_ind();
                    if (t_idx >= 0 &&
                        t_idx < static_cast<int32_t>(ctx.t_strand_arr_.size())) {
                        int32_t t_strand = ctx.t_strand_arr_[t_idx];
                        if (t_strand == STRAND_POS ||
                            t_strand == STRAND_NEG) {
                            strand_obs.exonic.add(result.align_strand, t_strand);
                        }
                    }
                }

            }

            // Every resolved non-chimeric hit is a placement for calibration; the fragment is offered
            // whole after its last hit (`offer_placements`).
            collect_placement(frag, cr);

            accumulator.append(result, frag_id);

            if (!is_unique_mapper) {
                n_buffered_mm++;
            }
        }

        // ⭐ THE OFFER: one per fragment, every placement at once.
        offer_placements();

        // Chimera counting: per-fragment, not per-hit.
        // A fragment is chimeric only if it resolves but every
        // resolved hit is chimeric (no non-chimeric alternative).
        if (any_hit_chimeric && !any_non_chimeric_resolved) {
            stats.n_chimeric++;
            if (worst_chimera_type == CHIMERA_TRANS)
                stats.n_chimeric_trans++;
            else if (worst_chimera_type == CHIMERA_CIS_STRAND_SAME)
                stats.n_chimeric_cis_strand_same++;
            else if (worst_chimera_type == CHIMERA_CIS_STRAND_DIFF)
                stats.n_chimeric_cis_strand_diff++;
        }

        // Intergenic counting: only if NO hit resolved to transcripts.
        // For a multimapper, if any hit resolves, the fragment goes
        // through the EM which handles gDNA allocation — don't also
        // count it as intergenic (that would double-count the molecule).
        if (!any_hit_resolved && any_unresolved) stats.n_intergenic++;

        if (n_buffered_mm > 0) {
            stats.n_multimapper_groups++;
            stats.n_multimapper_alignments += n_buffered_mm;
        }
    }

    // ----------------------------------------------------------------
    // Build Python result dict
    // ----------------------------------------------------------------

    nb::dict build_result() {
        nb::dict result;

        // Stats
        nb::dict stats_dict;
        stats_dict["total"] = stats_.total;
        stats_dict["qc_fail"] = stats_.qc_fail;
        stats_dict["unmapped"] = stats_.unmapped;
        stats_dict["secondary"] = stats_.secondary;
        stats_dict["supplementary"] = stats_.supplementary;
        stats_dict["duplicate"] = stats_.duplicate;
        stats_dict["n_read_names"] = stats_.n_read_names;
        stats_dict["unique"] = stats_.unique;
        stats_dict["multimapping"] = stats_.multimapping;
        stats_dict["proper_pair"] = stats_.proper_pair;
        stats_dict["improper_pair"] = stats_.improper_pair;
        stats_dict["mate_unmapped"] = stats_.mate_unmapped;
        stats_dict["n_fragments"] = stats_.n_fragments;
        stats_dict["n_chimeric"] = stats_.n_chimeric;
        stats_dict["n_chimeric_trans"] = stats_.n_chimeric_trans;
        stats_dict["n_chimeric_cis_strand_same"] = stats_.n_chimeric_cis_strand_same;
        stats_dict["n_chimeric_cis_strand_diff"] = stats_.n_chimeric_cis_strand_diff;
        stats_dict["n_intergenic"] = stats_.n_intergenic;
        stats_dict["n_with_exon"] = stats_.n_with_exon;
        stats_dict["n_with_annotated_sj"] = stats_.n_with_annotated_sj;
        stats_dict["n_same_strand"] = stats_.n_same_strand;
        stats_dict["n_ambig_strand"] = stats_.n_ambig_strand;
        stats_dict["n_strand_trained"] = stats_.n_strand_trained;
        stats_dict["n_strand_skipped_no_sj"] = stats_.n_strand_skipped_no_sj;
        stats_dict["n_strand_skipped_ambig_strand"] = stats_.n_strand_skipped_ambig_strand;
        stats_dict["n_strand_skipped_ambiguous"] = stats_.n_strand_skipped_ambiguous;
        stats_dict["n_multimapper_groups"] = stats_.n_multimapper_groups;
        stats_dict["n_multimapper_alignments"] = stats_.n_multimapper_alignments;
        stats_dict["n_sj_observed"] = stats_.n_sj_observed;
        stats_dict["n_sj_blacklisted"] = stats_.n_sj_blacklisted;
        stats_dict["n_sj_tag_tagged"] = stats_.sj_tags.tagged;
        stats_dict["n_sj_tag_untagged"] = stats_.sj_tags.untagged;
        stats_dict["n_sj_tag_conflicting"] = stats_.sj_tags.conflicting;
        stats_dict["n_sj_tag_invalid"] = stats_.sj_tags.invalid;
        stats_dict["n_nm_tag_missing"] = stats_.nm_tag.missing;
        stats_dict["n_nm_tag_invalid"] = stats_.nm_tag.invalid;
        stats_dict["n_nm_conflicting"] = stats_.n_nm_conflicting;
        stats_dict["n_nh_tag_missing"] = stats_.nh_tag.missing;
        stats_dict["n_nh_tag_invalid"] = stats_.nh_tag.invalid;
        for (int32_t v = HIT_USABLE + 1; v < N_HIT_VERDICTS; v++)
            stats_dict[("n_hits_rejected_" + std::string(hit_verdict_label(v))).c_str()] = stats_.hits_rejected[v];
        for (int32_t w = WITHHELD_NONE + 1; w < N_WITHHELD; w++)
            stats_dict[("n_withheld_" + std::string(withheld_label(w))).c_str()] = stats_.withheld[w];
        stats_dict["n_nh_above_reported"] = stats_.n_nh_above_reported;
        // ⭐ The splice census, keyed off the ONE name table. `splice_type_label`'s strings are the
        // `SpliceType` member names lower-cased, so `rigel.splice.census_field` derives the very same
        // key on the Python side without a second table to drift from this one.
        for (size_t i = 0; i < stats_.splice_census.size(); i++) {
            stats_dict[("n_census_" + std::string(splice_type_label(static_cast<int>(i)))).c_str()] =
                stats_.splice_census[i];
        }
        stats_dict["n_deposit_not_offered"] = stats_.n_deposit_not_offered;
        result["stats"] = stats_dict;

        // Strand observations
        nb::dict strand_dict;
        // The all-exonic diagnostic's 2×2, cell by cell.
        strand_dict["exonic_pos_pos"]       = strand_obs_.exonic.pos_pos;
        strand_dict["exonic_pos_neg"]       = strand_obs_.exonic.pos_neg;
        strand_dict["exonic_neg_pos"]       = strand_obs_.exonic.neg_pos;
        strand_dict["exonic_neg_neg"]       = strand_obs_.exonic.neg_neg;

        // The per-sj strand table, as six parallel arrays.  Sorted by
        // (ref, start, end, motif strand) because the source is an unordered_map
        // merged across workers, and every downstream number must not depend on
        // thread scheduling or hash order.
        {
            std::vector<SJKey> keys;
            keys.reserve(strand_obs_.sj_strand_table.size());
            for (const auto& [key, counts] : strand_obs_.sj_strand_table) keys.push_back(key);
            std::sort(keys.begin(), keys.end(), [](const SJKey& a, const SJKey& b) {
                if (a.ref_id != b.ref_id) return a.ref_id < b.ref_id;
                if (a.start  != b.start)  return a.start  < b.start;
                if (a.end    != b.end)    return a.end    < b.end;
                return a.strand < b.strand;
            });
            size_t n_sj = keys.size();
            std::vector<int32_t>  sj_ref(n_sj);
            std::vector<int64_t>  sj_start(n_sj), sj_end(n_sj);
            std::vector<int8_t>   sj_motif(n_sj);
            std::vector<uint64_t> sj_sense(n_sj), sj_antisense(n_sj);
            for (size_t i = 0; i < n_sj; i++) {
                const auto& c = strand_obs_.sj_strand_table.at(keys[i]);
                sj_ref[i]       = keys[i].ref_id;
                sj_start[i]     = keys[i].start;   // int32 upstream; widened once, here
                sj_end[i]       = keys[i].end;
                sj_motif[i]     = static_cast<int8_t>(keys[i].strand);
                sj_sense[i]     = c.n_sense;
                sj_antisense[i] = c.n_antisense;
            }
            strand_dict["sj_ref_id"]       = vec_to_ndarray(std::move(sj_ref));
            strand_dict["sj_start"]        = vec_to_ndarray(std::move(sj_start));
            strand_dict["sj_end"]          = vec_to_ndarray(std::move(sj_end));
            strand_dict["sj_motif_strand"] = vec_to_ndarray(std::move(sj_motif));
            strand_dict["sj_n_sense"]      = vec_to_ndarray(std::move(sj_sense));
            strand_dict["sj_n_antisense"]  = vec_to_ndarray(std::move(sj_antisense));
        }
        result["strand_observations"] = strand_dict;

        // ── the accumulator payload ───────────────────────────────────────────────────────────────────
        //
        // ⚠ The keys here are the FIELD NAMES OF THE SPECIFICATION'S `Tally`, character for character
        // (`tests/native/_accumulator_reference.py`). That is deliberate: the payload, the reference and
        // the parity gate then all say one name per quantity, and no mapping table can drift.
        //
        // ⚠ NOT zero-copy, and the comment that used to claim it was is deleted rather than carried
        // forward. Each array is a fresh owning vector concatenated across references, so Python gets one
        // contiguous buffer and does not depend on the AccumulatorSet's internal layout. Arrays are emitted
        // FLAT; the two-column ones are reshaped Python-side.
        //
        // ⚠ And the boundary and sj copies cannot be a memcpy: the accumulator stores AoS — a 40 B `Boundary`
        // interleaves two uint32 count pairs with three float64 channels, a 32 B `SpliceJunction` one count pair
        // with three float64s — while the payload is SoA, so their transpose is a strided read by construction.
        // (The 8 B `Region` is one count pair, already the payload's layout.) That is the price of hot structs
        // with no padding, which is the right trade — this loop is O(partition) once per scan, the deposit is
        // O(fragments).
        if (acc_set_) {
            using rigel::accumulator::Accumulator;
            using rigel::accumulator::Boundary;
            using rigel::accumulator::SpliceJunction;
            using rigel::accumulator::Region;
            using rigel::accumulator::kNFragmentPools;
            using rigel::accumulator::kNStrandColumns;

            const std::size_t n_refs = acc_set_->n_refs();

            std::vector<int64_t> ref_region_offsets(n_refs + 1, 0);
            std::vector<int64_t> ref_boundary_offsets(n_refs + 1, 0);
            std::vector<int64_t> ref_sj_offsets(n_refs + 1, 0);
            for (std::size_t f = 0; f < n_refs; ++f) {
                const Accumulator& a = acc_set_->at(static_cast<int32_t>(f));
                ref_region_offsets[f + 1] =
                    ref_region_offsets[f] + static_cast<int64_t>(a.n_regions());
                ref_boundary_offsets[f + 1] =
                    ref_boundary_offsets[f] + static_cast<int64_t>(a.n_boundaries());
                ref_sj_offsets[f + 1] =
                    ref_sj_offsets[f] + static_cast<int64_t>(a.n_sj());
            }
            const auto n_regions = static_cast<std::size_t>(ref_region_offsets.back());
            const auto n_boundaries = static_cast<std::size_t>(ref_boundary_offsets.back());
            const auto n_sj    = static_cast<std::size_t>(ref_sj_offsets.back());

            std::vector<uint32_t> region_contained_count(n_regions * kNStrandColumns, 0u);
            std::vector<uint32_t> region_start_count(n_regions * kNStrandColumns, 0u);
            std::vector<uint32_t> boundary_unspliced_count(n_boundaries * kNStrandColumns, 0u);
            std::vector<double> boundary_unspliced_inv_length_sum(n_boundaries, 0.0);
            std::vector<uint32_t> boundary_spliced_count(n_boundaries * kNStrandColumns, 0u);
            // ⚠ NOT multiplied by kNStrandColumns — the conserved mass has ONE value per boundary.
            std::vector<double> boundary_unspliced_mass(n_boundaries, 0.0);
            std::vector<double> boundary_spliced_mass(n_boundaries, 0.0);
            std::vector<uint32_t> sj_count(n_sj * kNStrandColumns, 0u);
            std::vector<double> sj_inv_length_sum(n_sj, 0.0);
            // MULTIPLIED by kNStrandColumns, unlike every other mass here: the sj mass is the
            // one that carries a strand (premise in `SpliceJunction::mass`).
            std::vector<double> sj_mass(n_sj * kNStrandColumns, 0.0);

            const std::size_t pool_row = static_cast<std::size_t>(max_length_) + 1;
            std::vector<int64_t> pool_lengths(kNFragmentPools * pool_row, 0);
            // ⭐ The unconditional length histogram, summed over references exactly like the pools.
            std::vector<uint32_t> deposited_lengths(pool_row, 0u);

            rigel::accumulator::DepositCounters qc;
            rigel::accumulator::GapCensus      gap_census;
            // ⭐ The side buffer, concatenated in REFERENCE ORDER — which is already the one canonical
            // order, and needs no second sort. Each reference's bank is canonicalised on its own, and
            // every record in it carries that reference's id, so `(ref, start, end, …)` ascends across the
            // concatenation by construction. `merge_from` is exactly a CSR concatenation.
            rigel::accumulator::DeferredFragments deferred;

            for (std::size_t f = 0; f < n_refs; ++f) {
                Accumulator& a = acc_set_->at(static_cast<int32_t>(f));
                deferred.merge_from(a.deferred_canonical());
                gap_census.merge_from(a.gap_census());
                const Region* regions = a.regions_data();
                const Boundary* boundaries = a.boundaries_data();
                const SpliceJunction* sj = a.sj_data();
                const uint32_t* starts = a.region_start_count_data();
                const auto region_base = static_cast<std::size_t>(ref_region_offsets[f]);
                const auto boundary_base = static_cast<std::size_t>(ref_boundary_offsets[f]);
                const auto sj_base   = static_cast<std::size_t>(ref_sj_offsets[f]);

                for (std::size_t i = 0; i < a.n_regions(); ++i) {
                    for (std::size_t c = 0; c < kNStrandColumns; ++c) {
                        region_start_count[(region_base + i) * kNStrandColumns + c] = starts[i * kNStrandColumns + c];
                    }
                    for (std::size_t c = 0; c < kNStrandColumns; ++c) {
                        const std::size_t o = (region_base + i) * kNStrandColumns + c;
                        region_contained_count[o]   = regions[i].contained_count[c];
                    }
                }
                for (std::size_t i = 0; i < a.n_boundaries(); ++i) {
                    for (std::size_t c = 0; c < kNStrandColumns; ++c) {
                        const std::size_t o = (boundary_base + i) * kNStrandColumns + c;
                        boundary_unspliced_count[o]   = boundaries[i].unspliced_count[c];
                        boundary_spliced_count[o]     = boundaries[i].spliced_count[c];
                    }
                    // ⚠ Outside the column loop, and keyed WITHOUT kNStrandColumns: one value per boundary.
                    boundary_unspliced_inv_length_sum[boundary_base + i] = boundaries[i].unspliced_inv_length_sum;
                    boundary_unspliced_mass[boundary_base + i] = boundaries[i].unspliced_mass;
                    boundary_spliced_mass[boundary_base + i]   = boundaries[i].spliced_mass;
                }
                for (std::size_t i = 0; i < a.n_sj(); ++i) {
                    for (std::size_t c = 0; c < kNStrandColumns; ++c) {
                        const std::size_t o = (sj_base + i) * kNStrandColumns + c;
                        sj_count[o]   = sj[i].count[c];
                        sj_mass[o]    = sj[i].mass[c];
                    }
                    sj_inv_length_sum[sj_base + i] = sj[i].inv_length_sum;
                }
                // The pools are library-wide, so the per-reference histograms are SUMMED, not concatenated.
                // ⚠ A size mismatch throws instead of being skipped: skipping would drop that reference's
                // pools with nothing to notice it by, which is the failure mode this rework has already
                // removed twice.
                if (a.pool_lengths_size() != pool_lengths.size()) {
                    throw std::runtime_error(
                        "build_result: reference " + std::to_string(f) + " has " +
                        std::to_string(a.pool_lengths_size()) + " pool bins but the payload expects " +
                        std::to_string(pool_lengths.size()));
                }
                const int64_t* pools = a.pool_lengths_data();
                for (std::size_t i = 0; i < pool_lengths.size(); ++i) pool_lengths[i] += pools[i];
                // ⭐ The unconditional length histogram takes the same size guard for the same reason: a silently skipped reference would leave the
                // anchor short by that reference's fragments with nothing to notice it by, and the
                // sum(deposited_lengths) == deposited invariant is what would fire -- but only if the
                // array is the right length in the first place.
                if (a.deposited_lengths_size() != deposited_lengths.size()) {
                    throw std::runtime_error(
                        "build_result: reference " + std::to_string(f) + " has " +
                        std::to_string(a.deposited_lengths_size()) +
                        " deposited-length bins but the payload expects " +
                        std::to_string(deposited_lengths.size()));
                }
                const uint32_t* dep = a.deposited_lengths_data();
                for (std::size_t i = 0; i < deposited_lengths.size(); ++i) deposited_lengths[i] += dep[i];
                qc.merge_from(a.counters());
            }
            // ⭐ The multi-placement records live on the SET (a record spans references), so the bank is
            // the union of every reference's gap bank and the set's own, in ONE canonical order: the
            // reference-order concatenation above is canonical only among one-placement records.
            deferred.merge_from(acc_set_->deferred_canonical());
            deferred.canonicalise();
            qc.merge_from(acc_set_->counters());

            nb::dict cal;
            // Echo the partition back, so a consumer can locate every object without reloading the index.
            cal["region_bounds"]   = vec_to_ndarray(std::vector<int64_t>(region_bounds_));
            cal["ref_region_bound_offsets"] = vec_to_ndarray(std::vector<int64_t>(ref_region_bound_offsets_));
            cal["ref_region_offsets"] = vec_to_ndarray(std::move(ref_region_offsets));
            cal["ref_boundary_offsets"] = vec_to_ndarray(std::move(ref_boundary_offsets));
            cal["ref_sj_offsets"]   = vec_to_ndarray(std::move(ref_sj_offsets));

            cal["region_contained_count"]   = vec_to_ndarray(std::move(region_contained_count));
            cal["region_start_count"]       = vec_to_ndarray(std::move(region_start_count));
            cal["boundary_unspliced_count"]   = vec_to_ndarray(std::move(boundary_unspliced_count));
            cal["boundary_unspliced_inv_length_sum"] = vec_to_ndarray(std::move(boundary_unspliced_inv_length_sum));
            cal["boundary_spliced_count"]     = vec_to_ndarray(std::move(boundary_spliced_count));
            cal["boundary_unspliced_mass"]    = vec_to_ndarray(std::move(boundary_unspliced_mass));
            cal["boundary_spliced_mass"]      = vec_to_ndarray(std::move(boundary_spliced_mass));
            cal["sj_count"]               = vec_to_ndarray(std::move(sj_count));
            cal["sj_inv_length_sum"]             = vec_to_ndarray(std::move(sj_inv_length_sum));
            cal["sj_mass"]                = vec_to_ndarray(std::move(sj_mass));
            cal["pool_lengths"]           = vec_to_ndarray(std::move(pool_lengths));
            cal["deposited_lengths"]      = vec_to_ndarray(std::move(deposited_lengths));

            // ⭐ The QC denominators, the side buffer and the umbrella census, through the SAME three
            // exporters the parity surface uses, so the payload and the bound Accumulator cannot disagree
            // about a key.
            cal["qc"]             = qc_dict(qc);
            cal["deferred"]       = deferred_dict(deferred);
            cal["gap_resolution"] = gap_census_dict(gap_census);

            cal["n_strand_columns"] = static_cast<int>(kNStrandColumns);
            cal["n_fragment_pools"] = static_cast<int>(kNFragmentPools);
            cal["max_length"]       = max_length_;
            cal["n_refs"]           = static_cast<int>(n_refs);

            result["calibration"] = cal;
        } else {
            result["calibration"] = nb::none();
        }

        return result;
    }
};

// ================================================================
// BamAnnotationWriter — the annotated-BAM pass's tag stamping
// ================================================================
//
// Reuses the same BAM reading, hit grouping, fragment building, and
// multimapper pairing infrastructure as BamScanner, but instead of
// resolving + buffering, it:
//   1. Looks up per-fragment annotation by frag_id
//   2. Stamps Z* BAM tags on every record
//   3. Writes tagged records to an output BAM
//
// The annotation table is passed as sliced NumPy arrays (one entry
// per annotated fragment).  Unannotated frag_ids (intergenic) are
// stamped with intergenic defaults.

// ndarray aliases (const, 1-D, C-contiguous)
using ci64_1d = nb::ndarray<const int64_t,  nb::ndim<1>, nb::c_contig>;
using ci32_1d = nb::ndarray<const int32_t,  nb::ndim<1>, nb::c_contig>;
using cu8_1d  = nb::ndarray<const uint8_t,  nb::ndim<1>, nb::c_contig>;
using cf32_1d = nb::ndarray<const float,    nb::ndim<1>, nb::c_contig>;
using ci8_1d  = nb::ndarray<const int8_t,   nb::ndim<1>, nb::c_contig>;
using ci16_1d = nb::ndarray<const int16_t,  nb::ndim<1>, nb::c_contig>;

// Rigel owns the tags it stamps: an aligner's tag of the same name (HISAT2 writes ZS:i, its second-best score) is
// replaced, never kept — htslib will not update a tag in place across types — and a failed write throws.
static void set_tag_int(bam1_t* r, const char* tag, int64_t v) {
    if (uint8_t* old = bam_aux_get(r, tag)) bam_aux_del(r, old);
    if (bam_aux_update_int(r, tag, v) < 0)
        throw std::runtime_error(std::string("Failed to write the ") + tag + " tag");
}

static void set_tag_float(bam1_t* r, const char* tag, float v) {
    if (uint8_t* old = bam_aux_get(r, tag)) bam_aux_del(r, old);
    if (bam_aux_update_float(r, tag, v) < 0)
        throw std::runtime_error(std::string("Failed to write the ") + tag + " tag");
}

static void set_tag_str(bam1_t* r, const char* tag, const char* v) {
    if (uint8_t* old = bam_aux_get(r, tag)) bam_aux_del(r, old);
    if (bam_aux_update_str(r, tag, static_cast<int>(strlen(v) + 1), v) < 0)
        throw std::runtime_error(std::string("Failed to write the ") + tag + " tag");
}

class BamAnnotationWriter {
public:
    // Set by the one constructor, from the caller's scan config (pass 1's).
    FragmentResolver* ctx_;
    SJTagMode sj_tag_mode_;
    bool skip_duplicates_;
    bool include_multimap_;
    std::vector<int32_t> tid_to_ref_id_;

    BamAnnotationWriter(FragmentResolver& ctx,
                        const std::string& sj_tag_spec,
                        bool skip_duplicates,
                        bool include_multimap)
        : ctx_(&ctx),
          sj_tag_mode_(parse_sj_tag_spec(sj_tag_spec)),
          skip_duplicates_(skip_duplicates),
          include_multimap_(include_multimap)
    {}

    // ----------------------------------------------------------------
    // Main write entry point
    // ----------------------------------------------------------------

    nb::dict write(
        const std::string& bam_path,
        const std::string& output_path,
        // Annotation arrays (sliced to n_annotations, indexed by row)
        ci64_1d ann_frag_ids,
        ci32_1d ann_best_tid,
        ci32_1d ann_best_gid,
        cu8_1d  ann_tx_flags,
        cf32_1d ann_posterior,
        ci8_1d  ann_frag_class,
        ci16_1d ann_n_candidates,
        cu8_1d  ann_splice_type,
        ci32_1d ann_locus_id,
        int64_t ann_size,
        // String ID lookup
        const std::vector<std::string>& t_ids,
        const std::vector<std::string>& g_ids,
        const std::vector<std::string>& g_names)
    {
        // Build frag_id → row lookup
        std::unordered_map<int64_t, int64_t> frag_to_row;
        frag_to_row.reserve(static_cast<size_t>(ann_size));
        const int64_t* fid_ptr = ann_frag_ids.data();
        for (int64_t i = 0; i < ann_size; i++) {
            frag_to_row[fid_ptr[i]] = i;
        }

        // Raw data pointers
        const int32_t* tid_ptr  = ann_best_tid.data();
        const int32_t* gid_ptr  = ann_best_gid.data();
        const uint8_t* flags_ptr = ann_tx_flags.data();
        const float*   post_ptr = ann_posterior.data();
        const int8_t*  fc_ptr   = ann_frag_class.data();
        const int16_t* nc_ptr   = ann_n_candidates.data();
        const uint8_t* st_ptr   = ann_splice_type.data();
        const int32_t* lid_ptr  = ann_locus_id.data();

        // Open input BAM
        htsFile* fp = hts_open(bam_path.c_str(), "rb");
        if (!fp)
            throw std::runtime_error("Failed to open BAM: " + bam_path);

        bam_hdr_t* hdr = sam_hdr_read(fp);
        if (!hdr) {
            hts_close(fp);
            throw std::runtime_error("Failed to read header: " + bam_path);
        }

        // Build tid mapping
        build_tid_mapping(hdr, ctx_->ref_to_id_, tid_to_ref_id_);

        // Open output BAM
        htsFile* out = hts_open(output_path.c_str(), "wb");
        if (!out) {
            bam_hdr_destroy(hdr);
            hts_close(fp);
            throw std::runtime_error("Failed to open output: " + output_path);
        }
        if (sam_hdr_write(out, hdr) < 0) {
            hts_close(out);
            bam_hdr_destroy(hdr);
            hts_close(fp);
            throw std::runtime_error("Failed to write header: " + output_path);
        }

        bam1_t* b = bam_init1();

        // Qname grouping — collect raw records + ParsedAlignments
        std::vector<bam1_t*> raw_group;       // cloned bam1_t*
        std::vector<ParsedAlignment> light_group;  // parsed lightweight records
        std::string current_qname;
        int64_t frag_id = 0;

        // Summary counters
        int64_t n_read_groups = 0;
        int64_t n_annotated = 0;
        int64_t n_intergenic = 0;
        int64_t n_records_written = 0;
        int64_t n_filtered_passthrough = 0;

        // Scratch buffer for resolve calls (single-threaded, reused)
        ResolverScratch scratch(ctx_->n_transcripts_);

        auto process_and_write = [&]() {
            if (light_group.empty()) return;

            n_read_groups++;

            if (!include_multimap_ && is_multimapper(light_group, template_nh(light_group))) {
                // Write records without annotation (match pass 1 skip)
                // Pass 1 assigns frag_id per qname group before worker-side
                // multimapper filtering, so we must still advance frag_id
                // here to keep pass-2 lookup synchronized.
                for (size_t i = 0; i < raw_group.size(); i++) {
                    set_tag_int(
                        raw_group[i], "ZB",
                        static_cast<int>(light_group[i].n_sj_blacklisted));
                }
                stamp_and_write_hit(raw_group, out, hdr, ".", ".", ".",
                                    -1, -1,
                                    AF_MULTIMAPPER_DROP, 0.0f,
                                    ".", 1, 0, "unknown", -1);
                n_intergenic++;
                n_records_written += static_cast<int64_t>(raw_group.size());
                frag_id++;
                return;
            }

            // The template's candidate hits, as the scan formed them
            auto hit_result = template_hits(light_group);
            auto& all_hits = hit_result.hits;

            // The scan's verdicts, replayed (`judge_hits`): a rejected hit's records are written as a secondary
            // hit's, and a template left with no hit is written unresolved, its class the outcome's label.
            auto& verdicts = scratch.hit_verdicts;
            const int32_t withheld = judge_hits(hit_result, *ctx_, scratch, verdicts);
            if (withheld != WITHHELD_NONE) {
                for (size_t i = 0; i < raw_group.size(); i++) {
                    set_tag_int(
                        raw_group[i], "ZB",
                        static_cast<int>(light_group[i].n_sj_blacklisted));
                }
                stamp_and_write_hit(raw_group, out, hdr, ".", ".", ".",
                                    -1, -1,
                                    AF_UNRESOLVED, 0.0f,
                                    withheld_label(withheld), 0, 0, "unknown", -1);
                n_records_written += static_cast<int64_t>(raw_group.size());
                frag_id++;
                return;
            }

            const int32_t num_hits = static_cast<int32_t>(std::count(verdicts.begin(), verdicts.end(), HIT_USABLE));
            const int first_valid_hit = static_cast<int>(
                std::find(verdicts.begin(), verdicts.end(), HIT_USABLE) - verdicts.begin());

            // Look up annotation for this frag_id
            auto ann_it = frag_to_row.find(frag_id);
            bool has_ann = (ann_it != frag_to_row.end());

            // Each input record is written EXACTLY once, with the first hit that holds it: one mate location joins
            // several hits. An alias of an alignment is written with the record that represents it
            // (`TemplateHits::alias_of`).
            std::vector<bool> rec_written(raw_group.size(), false);
            std::vector<std::vector<int32_t>> aliases(raw_group.size());
            for (size_t i = 0; i < hit_result.alias_of.size(); ++i)
                if (hit_result.alias_of[i] >= 0) aliases[hit_result.alias_of[i]].push_back(static_cast<int32_t>(i));

            for (int hit_idx = 0;
                 hit_idx < static_cast<int>(all_hits.size());
                 hit_idx++)
            {
                const auto& [r1_reads, r2_reads] = all_hits[hit_idx];

                // Gather bam1_t* for all records in this hit, and their aliases, skipping any already emitted by an
                // earlier hit. Stamp the per-record ZB:i tag (splice-artifact count) on each bam1_t as we select it:
                // a record-local property of the CIGAR, not a hit-level value.
                std::vector<bam1_t*> hit_raws;
                auto take = [&](ptrdiff_t idx) {
                    if (rec_written[idx]) return;
                    rec_written[idx] = true;
                    set_tag_int(raw_group[idx], "ZB", static_cast<int>(light_group[idx].n_sj_blacklisted));
                    hit_raws.push_back(raw_group[idx]);
                };
                for (const auto* reads : {&r1_reads, &r2_reads}) {
                    for (const auto* lr : *reads) {
                        const ptrdiff_t idx = lr - &light_group[0];
                        take(idx);
                        for (int32_t a : aliases[idx]) take(a);
                    }
                }
                if (hit_raws.empty()) continue;

                if (has_ann) {
                    int64_t row = ann_it->second;
                    int32_t best_tid_val = tid_ptr[row];
                    int32_t best_gid_val = gid_ptr[row];
                    uint8_t flags_val    = flags_ptr[row];
                    float   post_val     = post_ptr[row];
                    int8_t  fc_val       = fc_ptr[row];
                    int16_t nc_val       = nc_ptr[row];
                    uint8_t st_val       = st_ptr[row];

                    const char* t_id_str = (best_tid_val >= 0 &&
                        best_tid_val < static_cast<int32_t>(t_ids.size()))
                        ? t_ids[best_tid_val].c_str() : ".";
                    const char* g_id_str = (best_gid_val >= 0 &&
                        best_gid_val < static_cast<int32_t>(g_ids.size()))
                        ? g_ids[best_gid_val].c_str() : ".";
                    const char* g_name_str = (best_gid_val >= 0 &&
                        best_gid_val < static_cast<int32_t>(g_names.size()))
                        ? g_names[best_gid_val].c_str() : ".";

                    // Determine primary hit
                    bool is_primary = false;
                    if (verdicts[hit_idx] != HIT_USABLE) {
                        is_primary = false;
                    } else if (num_hits == 1) {
                        is_primary = true;
                    } else if (best_tid_val >= 0) {
                        // Resolve this hit to see if assigned transcript
                        // is among its candidates
                        AssembledFragment frag = build_fragment(
                            r1_reads, r2_reads);
                        if (!frag.exons.empty()) {
                            RawResolveResult cr;
                            ctx_->_resolve_core(
                                frag.exons, frag.introns,
                                frag.genomic_footprint(), cr, scratch);
                            for (int32_t ti : cr.t_inds) {
                                if (ti == best_tid_val) {
                                    is_primary = true;
                                    break;
                                }
                            }
                        }
                    }

                    int32_t lid_val = lid_ptr[row];

                    stamp_and_write_hit(
                        hit_raws, out, hdr,
                        t_id_str, g_id_str, g_name_str,
                        static_cast<int>(best_tid_val),
                        static_cast<int>(best_gid_val),
                        static_cast<int>(flags_val), post_val,
                        frag_class_label(fc_val),
                        is_primary ? 1 : 0,
                        static_cast<int>(nc_val),
                        splice_type_label(st_val),
                        static_cast<int>(lid_val));

                    n_annotated++;
                } else {
                    // No EM annotation for this frag_id: the deterministic
                    // intergenic fast-path. It is gDNA by construction (a spliced
                    // hit the scan keeps has a transcript holding its junctions),
                    // and ZF carries the intergenic+gDNA bits, so a single AND-mask
                    // on ZF captures total gDNA. ZW=1.0 because there is no
                    // competing hypothesis.
                    stamp_and_write_hit(
                        hit_raws, out, hdr,
                        ".", ".", ".",
                        -1, -1,
                        AF_GDNA_INTERGENIC, 1.0f,
                        ".",
                        (hit_idx == first_valid_hit) ? 1 : 0,
                        0, "unknown", -1);
                    n_intergenic++;
                }

                n_records_written +=
                    static_cast<int64_t>(hit_raws.size());
            }

            // Defensive end-of-group sweep: any raw record not yet
            // emitted via a hit is written as a single intergenic
            // pass-through so output record count == input record count.
            std::vector<bam1_t*> orphan;
            for (size_t i = 0; i < raw_group.size(); i++) {
                if (!rec_written[i]) {
                    set_tag_int(
                        raw_group[i], "ZB",
                        static_cast<int>(light_group[i].n_sj_blacklisted));
                    orphan.push_back(raw_group[i]);
                }
            }
            if (!orphan.empty()) {
                stamp_and_write_hit(
                    orphan, out, hdr,
                    ".", ".", ".",
                    -1, -1,
                    AF_GDNA_INTERGENIC, 1.0f,
                    ".", 0, 0, "unknown", -1);
                n_intergenic++;
                n_records_written +=
                    static_cast<int64_t>(orphan.size());
            }

            frag_id++;
        };

        // ---- Main BAM read loop ----
        int read_ret = 0;
        while ((read_ret = sam_read1(fp, hdr, b)) >= 0) {
            uint16_t flag = b->core.flag;

            // Filtered reads: write through as-is (no annotation tags)
            // to preserve every record from the original BAM.  These
            // records are also skipped by Pass 1, so they never receive
            // a frag_id and must NOT advance the frag_id counter.
            bool filtered = (flag & BAM_FQCFAIL)
                         || (flag & BAM_FUNMAP)
                         || ((flag & BAM_FDUP) && skip_duplicates_)
                         || !(flag & BAM_FPAIRED);
            if (filtered) {
                // A filtered record that opens a new read name: write the held group first, or this read name
                // would be split around it and the output would no longer be collated
                if (!light_group.empty() && current_qname != bam_get_qname(b)) {
                    process_and_write();
                    for (auto* r : raw_group) bam_destroy1(r);
                    raw_group.clear();
                    light_group.clear();
                }
                if (sam_write1(out, hdr, b) < 0)
                    throw std::runtime_error(
                        "Failed to write filtered BAM record");
                n_filtered_passthrough++;
                n_records_written++;
                continue;
            }

            const char* qname = bam_get_qname(b);

            if (!light_group.empty() && current_qname != qname) {
                process_and_write();
                // Free cloned records
                for (auto* r : raw_group) bam_destroy1(r);
                raw_group.clear();
                light_group.clear();
            }

            current_qname = qname;

            // Clone raw record for output
            raw_group.push_back(bam_dup1(b));

            // The scan's own parser, the blacklist included, so the annotated BAM sees the evidence the scan saw;
            // its drop count is the per-record ZB:i tag stamped below.
            light_group.push_back(
                parse_bam_record(b, tid_to_ref_id_, sj_tag_mode_, ctx_, nullptr));
        }

        if (read_ret < -1) {
            for (auto* r : raw_group) bam_destroy1(r);
            bam_destroy1(b);
            bam_hdr_destroy(hdr);
            hts_close(fp);
            hts_close(out);
            throw read_error(bam_path, read_ret, n_records_written);
        }

        // Process last group
        if (!light_group.empty()) {
            process_and_write();
            for (auto* r : raw_group) bam_destroy1(r);
        }

        // Cleanup
        bam_destroy1(b);
        bam_hdr_destroy(hdr);
        hts_close(fp);
        hts_close(out);

        // Build summary dict
        nb::dict summary;
        summary["n_read_groups"] = n_read_groups;
        summary["n_annotated"] = n_annotated;
        summary["n_intergenic"] = n_intergenic;
        summary["n_records_written"] = n_records_written;
        summary["n_filtered_passthrough"] = n_filtered_passthrough;
        return summary;
    }

private:
    // Stamp BAM tags on all records in a hit and write to output.
    static void stamp_and_write_hit(
        const std::vector<bam1_t*>& records,
        htsFile* out,
        const bam_hdr_t* hdr,
        const char* zt,    // transcript ID
        const char* zg,    // gene ID
        const char* zr,    // gene name / symbol
        int zi,            // transcript index (-1 = unassigned)
        int zj,            // gene index (-1 = unassigned)
        int zf,            // assignment flags bitfield
        float zw,          // posterior
        const char* zc,    // fragment class label
        int zh,            // primary hit flag
        int zn,            // n_candidates
        const char* zs,    // splice type label
        int zl)            // locus_id (-1 = no locus)
    {
        for (bam1_t* r : records) {
            set_tag_str(r, "ZT", zt);
            set_tag_str(r, "ZG", zg);
            set_tag_str(r, "ZR", zr);
            set_tag_int(r, "ZI", zi);
            set_tag_int(r, "ZJ", zj);
            set_tag_int(r, "ZF", zf);
            set_tag_float(r, "ZW", zw);
            set_tag_str(r, "ZC", zc);
            set_tag_int(r, "ZH", zh);
            set_tag_int(r, "ZN", zn);
            set_tag_str(r, "ZS", zs);
            set_tag_int(r, "ZL", zl);

            if (sam_write1(out, hdr, r) < 0) {
                throw std::runtime_error("Failed to write BAM record");
            }
        }
    }
};

// ================================================================
// nanobind module definition
// ================================================================

NB_MODULE(_bam_impl, m) {
    m.doc() = "C++ BAM scanner and annotation writer for rigel using htslib.\n\n"
              "Provides BamScanner (pass 1: BAM → resolve → buffer) and\n"
              "BamAnnotationWriter (pass 2: stamp tags → write BAM).";

    // ----------------------------------------------------------------
    // Accumulator — ONE REFERENCE, bound as the parity surface
    // ----------------------------------------------------------------
    //
    // ⚠ THIS BINDING EXISTS TO BE COMPARED. `tests/native/test_accumulator_native_parity.py` drives it and
    // `tests/native/_accumulator_reference.py` with identical fragments and requires every array, dtype and
    // QC counter to match. The scan path does not go through here — it calls `deposit` in C++ — so the
    // property names below are chosen to be the SPECIFICATION'S `Tally` FIELD NAMES, character for
    // character. The gate reads them off `dataclasses.fields(Tally)`, so a name that drifts fails loudly
    // instead of quietly dropping out of the comparison.
    //
    // ⚠ Every view into a struct bank needs EXPLICIT STRIDES. The storage is AoS with mixed widths — a 40 B
    // `Boundary` holds two uint32 pairs then three float64s, a 32 B `SpliceJunction` one uint32 pair then three
    // float64s — so the row stride is `sizeof(struct)/sizeof(element)`, not the column count. The 8 B `Region`
    // is one uint32 pair, so there the two coincide; the START bank is a flat array and needs none.
    {
        using rigel::accumulator::Accumulator;
        using rigel::accumulator::Boundary;
        using rigel::accumulator::DepositScratch;
        using rigel::accumulator::GapHypothesis;
        using rigel::accumulator::OfferedFragment;
        using rigel::accumulator::SpliceJunction;
        using rigel::accumulator::Region;
        using rigel::accumulator::kNFragmentPools;
        using rigel::accumulator::kNStrandColumns;

        // One scratch per thread, shared by every bound `Accumulator`. The scan path uses
        // `WorkerState::offer_scratch.deposit` instead.
        static thread_local DepositScratch binding_scratch;

        nb::class_<Accumulator>(m, "Accumulator")
            .def("__init__",
                 [](Accumulator* self,
                    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> region_bounds,
                    nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig> region_types,
                    int max_length,
                    int32_t ref) {
                     std::vector<int64_t> c(region_bounds.data(), region_bounds.data() + region_bounds.shape(0));
                     std::vector<uint8_t> t(region_types.data(),
                                            region_types.data() + region_types.shape(0));
                     new (self) Accumulator(std::move(c), std::move(t), max_length, ref);
                 },
                 nb::arg("region_bounds"),
                 nb::arg("region_types"),
                 nb::arg("max_length"),
                 nb::arg("ref"),
                 "One reference's sorted region_bound positions, one coarse type per region, the\n"
                 "fragment-length limit (which is also the pool-histogram width), and WHICH\n"
                 "reference this is — stamped into every deferred record, because the second\n"
                 "pass replays those onto that reference's region_bound axis.")
            .def("set_sj",
                 [](Accumulator& a,
                    nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig> offsets,
                    nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig> boundary_right,
                    nb::ndarray<const int8_t, nb::ndim<1>, nb::c_contig> sj_strand) {
                     a.set_sj(
                         std::vector<int32_t>(offsets.data(), offsets.data() + offsets.shape(0)),
                         std::vector<int32_t>(boundary_right.data(),
                                              boundary_right.data() + boundary_right.shape(0)),
                         std::vector<int8_t>(sj_strand.data(),
                                             sj_strand.data() + sj_strand.shape(0)));
                 },
                 nb::arg("offsets"),
                 nb::arg("boundary_right"),
                 nb::arg("sj_strand"),
                 "The sj CSR for THIS reference, keyed by the ref-local donor region_bound\n"
                 "index. The sj-boundary id is the slot; slot order is a contract.")
            .def_prop_ro("n_regions",     [](const Accumulator& a) { return a.n_regions(); })
            .def_prop_ro("n_boundaries",     [](const Accumulator& a) { return a.n_boundaries(); })
            .def_prop_ro("n_sj", [](const Accumulator& a) { return a.n_sj(); })
            .def_prop_ro("max_length",  [](const Accumulator& a) { return a.max_length(); })

            // ── regions ────────────────────────────────────────────────────────────────────────────────
            .def_prop_ro("region_contained_count", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Region) / sizeof(uint32_t);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<2>>(
                    &a.regions_data()[0].contained_count[0], {a.n_regions(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("region_start_count", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<2>>(
                    a.region_start_count_data(), {a.n_regions(), kNStrandColumns}, h).cast();
            })

            // ── contiguous boundaries ─────────────────────────────────────────────────────────────────────
            .def_prop_ro("boundary_unspliced_count", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Boundary) / sizeof(uint32_t);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<2>>(
                    &a.boundaries_data()[0].unspliced_count[0], {a.n_boundaries(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("boundary_spliced_count", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Boundary) / sizeof(uint32_t);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<2>>(
                    &a.boundaries_data()[0].spliced_count[0], {a.n_boundaries(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("boundary_unspliced_inv_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Boundary) / sizeof(double);
                return nb::ndarray<nb::numpy, const double, nb::ndim<1>>(
                    &a.boundaries_data()[0].unspliced_inv_length_sum, {a.n_boundaries()}, h, {row}).cast();
            })
            // ⭐ The conserved mass — ndim<1>, because it has one value per boundary and not one per
            // strand. The stride is still the struct row, so this is a view into the same array.
            .def_prop_ro("boundary_unspliced_mass", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Boundary) / sizeof(double);
                return nb::ndarray<nb::numpy, const double, nb::ndim<1>>(
                    &a.boundaries_data()[0].unspliced_mass, {a.n_boundaries()}, h, {row}).cast();
            })
            .def_prop_ro("boundary_spliced_mass", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Boundary) / sizeof(double);
                return nb::ndarray<nb::numpy, const double, nb::ndim<1>>(
                    &a.boundaries_data()[0].spliced_mass, {a.n_boundaries()}, h, {row}).cast();
            })

            // ── sj boundaries ───────────────────────────────────────────────────────────────────────
            .def_prop_ro("sj_count", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(SpliceJunction) / sizeof(uint32_t);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<2>>(
                    &a.sj_data()[0].count[0], {a.n_sj(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("sj_inv_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(SpliceJunction) / sizeof(double);
                return nb::ndarray<nb::numpy, const double, nb::ndim<1>>(
                    &a.sj_data()[0].inv_length_sum, {a.n_sj()}, h, {row}).cast();
            })
            // ⭐ ndim<2>, unlike every other mass here — the sj mass carries a strand, and its
            // columns are `sj_count`'s columns. See `SpliceJunction::mass`.
            .def_prop_ro("sj_mass", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(SpliceJunction) / sizeof(double);
                return nb::ndarray<nb::numpy, const double, nb::ndim<2>>(
                    &a.sj_data()[0].mass[0], {a.n_sj(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })

            // ── the length pools, and the denominators ───────────────────────────────────────────────
            .def_prop_ro("pool_lengths", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                const std::size_t row = static_cast<std::size_t>(a.max_length()) + 1;
                return nb::ndarray<nb::numpy, const int64_t, nb::ndim<2>>(
                    a.pool_lengths_data(), {kNFragmentPools, row}, h).cast();
            })
            .def_prop_ro("deposited_lengths", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<1>>(
                    a.deposited_lengths_data(),
                    {static_cast<std::size_t>(a.max_length()) + 1}, h).cast();
            })
            .def_prop_ro("qc", [](const Accumulator& a) {
                return qc_dict(a.counters());
            })

            // ⭐ The umbrella census. The keys are the specification's `GapResolution` VALUES, character
            // for character, so the parity gate compares dicts and no mapping table exists to drift.
            .def_prop_ro("gap_resolution", [](const Accumulator& a) {
                return gap_census_dict(a.gap_census());
            })

            // ⭐ The deferred queue, flattened exactly as `Tally.deferred_arrays()` specifies and in the
            // SAME canonical order — the one bank whose order is observable, so the export sorts it.
            //
            // ⚠ Copies rather than viewing the C++ vectors. `canonicalise` may replace them wholesale on
            // the next call, and a numpy view over a vector that is about to be reseated is a dangling
            // pointer that reads as plausible numbers.
            .def_prop_ro("deferred", [](Accumulator& a) {
                return deferred_dict(a.deferred_canonical());
            })

            .def("region_of_pos",
                 [](const Accumulator& a, int64_t pos) { return a.region_of_pos(pos); },
                 nb::arg("pos"))

            // Returns the QC KEY the deposit incremented, which is the specification's
            // `DepositOutcome.value` — so the parity gate compares outcomes without an enum mapping, and
            // the string it compares is the name of the counter that moved.
            //
            // ⚠ `introns` here are `(start, end)` pairs on THIS reference. The scan-path adapter is what
            // restricts a fragment's introns to one reference and de-duplicates them; `deposit` itself
            // normalises by coordinate and never reads `IntronBlock::ref_id`.
            .def("deposit",
                 [](Accumulator& a,
                    int64_t start,
                    int64_t end,
                    nb::iterable introns,
                    int32_t align_strand,
                    int32_t sj_strand,
                    nb::object hypotheses) {
                     MarshalledFragment m(start, end, introns, align_strand, sj_strand, hypotheses);
                     return std::string(
                         rigel::accumulator::outcome_key(a.deposit(m.offered, binding_scratch)));
                 },
                 nb::arg("start"),
                 nb::arg("end"),
                 nb::arg("observed_introns") = nb::tuple(),
                 nb::arg("align_strand") = STRAND_POS,
                 nb::arg("sj_strand") = STRAND_NONE,
                 nb::arg("hypotheses"))

            // ⭐ L under ONE hypothesis, without depositing — what the second pass scores against.
            // ⛔ Exposed rather than reimplemented in Python: the pool rules left the tool with ONE definition of
            // fragment length, and a scorer that computed its own would be a second definition of exactly
            // the quantity that audit unified — and the one the drain would then disagree with.
            .def("length_under",
                 [](const Accumulator& a,
                    int64_t start,
                    int64_t end,
                    nb::iterable introns,
                    int32_t align_strand,
                    int32_t sj_strand,
                    nb::object hypotheses,
                    std::size_t hypothesis_index) {
                     MarshalledFragment m(start, end, introns, align_strand, sj_strand, hypotheses);
                     return a.length_under(m.offered, hypothesis_index, binding_scratch);
                 },
                 nb::arg("start"),
                 nb::arg("end"),
                 nb::arg("observed_introns") = nb::tuple(),
                 nb::arg("align_strand") = STRAND_POS,
                 nb::arg("sj_strand") = STRAND_NONE,
                 nb::arg("hypotheses") = nb::tuple(),
                 nb::arg("hypothesis_index") = 0)

            // ⭐ The annotated sj id for MANY introns at once — `sj_edge_id`, the deposit's own lookup,
            // applied per element. The second pass asks it once per hypothesis intron, which as a numpy
            // call each was 0.9 M calls on a deep library; per reference it is one call. The ids are
            // REF-LOCAL slots, the same space `set_sj` above pins to `ref_sj_offsets`, so a caller holding
            // the payload's flat axis adds that base. -1 where the intron is not annotated.
            .def("sj_edge_ids",
                 [](const Accumulator& a,
                    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> starts,
                    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> ends,
                    nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig> sj_strand) {
                     const std::size_t n = starts.shape(0);
                     if (ends.shape(0) != n || sj_strand.shape(0) != n)
                         throw std::invalid_argument(
                             "sj_edge_ids: starts, ends and sj_strand must be aligned per intron");
                     std::vector<int64_t> out(n);
                     a.sj_edge_ids(starts.data(), ends.data(), sj_strand.data(), n, out.data());
                     return rigel::vec_to_ndarray(std::move(out));
                 },
                 nb::arg("starts"),
                 nb::arg("ends"),
                 nb::arg("sj_strand"),
                 "The annotated sj id per intron, or -1 — `sj_edge_id` over arrays. Ids are\n"
                 "ref-local slots, as `set_sj` defines them.")

            // Element-wise sum of `other` into this accumulator — the per-worker merge the parallel scan
            // performs internally, exposed so the DETERMINISM contract can be tested directly: shard one
            // fragment corpus K ways, merge, and compare with the unsharded run. The count banks are integers
            // and merge bit for bit; the five float64 banks (the inverse-length sums and the conserved masses)
            // agree to within the float64 representation, because the merge re-associates their sums.
            .def("merge_from",
                 [](Accumulator& a, const Accumulator& other) { a.merge_from(other); },
                 nb::arg("other"));
    }


        // ⭐ The accumulator SET — the parity surface for `offer`, the multimapper rule. Built from the
        // flat partition exactly as the scanner builds its own; `at(ref)` hands out each reference's
        // accumulator, and `deferred` / `qc` / `gap_resolution` are the union the payload exports.
        nb::class_<rigel::accumulator::AccumulatorSet>(m, "AccumulatorSet")
            .def("__init__",
                 [](rigel::accumulator::AccumulatorSet* self,
                    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> region_bounds,
                    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> ref_region_bound_offsets,
                    nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig> region_types,
                    int max_length) {
                     new (self) rigel::accumulator::AccumulatorSet(
                         region_bounds.data(), region_bounds.shape(0), ref_region_bound_offsets.data(),
                         ref_region_bound_offsets.shape(0) - 1, region_types.data(),
                         region_types.shape(0), max_length);
                 },
                 nb::arg("region_bounds"), nb::arg("ref_region_bound_offsets"), nb::arg("region_types"),
                 nb::arg("max_length"))
            .def("set_sj",
                 [](rigel::accumulator::AccumulatorSet& s,
                    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> offsets,
                    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> boundary_right,
                    nb::ndarray<const int8_t, nb::ndim<1>, nb::c_contig> sj_strand,
                    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> ref_region_bound_offsets) {
                     s.set_sj(offsets.data(), offsets.shape(0), boundary_right.data(), sj_strand.data(),
                              sj_strand.shape(0), ref_region_bound_offsets.data());
                 },
                 nb::arg("offsets"), nb::arg("boundary_right"), nb::arg("sj_strand"),
                 nb::arg("ref_region_bound_offsets"))
            .def_prop_ro("n_refs", [](const rigel::accumulator::AccumulatorSet& s) { return s.n_refs(); })
            .def("at",
                 [](rigel::accumulator::AccumulatorSet& s, int32_t ref) -> rigel::accumulator::Accumulator& {
                     return s.at(ref);
                 },
                 nb::arg("ref"), nb::rv_policy::reference_internal)
            // One fragment with every placement it aligned to; each placement reads attributes off
            // whatever it is handed (`ref`, `start`, `end`, `align_strand`, `sj_strand`,
            // `observed_introns`, `hypotheses`), so the parity gate passes the specification's own
            // `Placement` objects. Returns the QC key the offer incremented.
            .def("offer",
                 [](rigel::accumulator::AccumulatorSet& s, nb::iterable placements) {
                     static rigel::accumulator::OfferScratch offer_scratch;
                     // ⚠ ONE pass over `placements`, as over each placement's hypotheses, so a one-shot
                     // iterable offers the same fragment as a list. The `Placement`s are taken after the
                     // loop, once `marshalled` has stopped growing.
                     std::vector<MarshalledFragment> marshalled;
                     std::vector<int32_t> refs;
                     for (nb::handle p : placements) {
                         marshalled.emplace_back(nb::cast<int64_t>(nb::getattr(p, "start")),
                                                 nb::cast<int64_t>(nb::getattr(p, "end")),
                                                 nb::cast<nb::iterable>(nb::getattr(p, "observed_introns")),
                                                 nb::cast<int32_t>(nb::getattr(p, "align_strand")),
                                                 nb::cast<int32_t>(nb::getattr(p, "sj_strand")),
                                                 nb::getattr(p, "hypotheses"));
                         refs.push_back(nb::cast<int32_t>(nb::getattr(p, "ref")));
                     }
                     std::vector<rigel::accumulator::Placement> pls;
                     pls.reserve(marshalled.size());
                     for (std::size_t i = 0; i < marshalled.size(); ++i) pls.push_back({refs[i], marshalled[i].offered});
                     return std::string(rigel::accumulator::outcome_key(s.offer(pls.data(), pls.size(), offer_scratch)));
                 },
                 nb::arg("placements"))
            .def_prop_ro("deferred", [](rigel::accumulator::AccumulatorSet& s) {
                rigel::accumulator::DeferredFragments all;
                for (std::size_t f = 0; f < s.n_refs(); ++f)
                    all.merge_from(s.at(static_cast<int32_t>(f)).deferred_canonical());
                all.merge_from(s.deferred_canonical());
                all.canonicalise();
                return deferred_dict(all);
            })
            .def_prop_ro("qc", [](const rigel::accumulator::AccumulatorSet& s) {
                rigel::accumulator::DepositCounters c = s.counters();
                for (std::size_t f = 0; f < s.n_refs(); ++f) c.merge_from(s.at(static_cast<int32_t>(f)).counters());
                return qc_dict(c);
            })
            .def_prop_ro("gap_resolution", [](const rigel::accumulator::AccumulatorSet& s) {
                rigel::accumulator::GapCensus g;
                for (std::size_t f = 0; f < s.n_refs(); ++f) g.merge_from(s.at(static_cast<int32_t>(f)).gap_census());
                return gap_census_dict(g);
            });


    nb::class_<BamScanner>(m, "BamScanner")
        .def(nb::init<FragmentResolver&, const std::string&, bool, bool>(),
             nb::arg("ctx"),
             nb::arg("sj_tag_spec"),
             nb::arg("skip_duplicates"),
             nb::arg("include_multimap"),
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
             "    Include multimapping reads.\n",
             nb::keep_alive<1, 2>()  // BamScanner must keep FragmentResolver alive
        )
        .def("scan", &BamScanner::scan,
             nb::arg("bam_path"),
             nb::arg("chunk_callback"),
             nb::arg("chunk_size"),
             nb::arg("n_workers"),
             nb::arg("n_decomp_threads"),
             nb::arg("qname_batch_size"),
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
          .def("set_regions",
                 [](BamScanner& self,
                     nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> region_bounds,
                     nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> ref_region_bound_offsets,
                     int32_t n_refs,
                     nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig> region_types,
                     int max_length) {
                      self.set_regions(region_bounds, ref_region_bound_offsets, n_refs,
                                       region_types, max_length);
                 },
                 nb::arg("region_bounds"),
                 nb::arg("ref_region_bound_offsets"),
                 nb::arg("n_refs"),
                 nb::arg("region_types"),
                 nb::arg("max_length"),
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
          .def("set_sj",
                 [](BamScanner& self,
                     nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> offsets,
                     nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> boundary_right,
                     nb::ndarray<const int8_t, nb::ndim<1>, nb::c_contig> sj_strand) {
                      self.set_sj(offsets, boundary_right, sj_strand);
                 },
                 nb::arg("offsets"),
                 nb::arg("boundary_right"),
                 nb::arg("sj_strand"),
                 "Install the annotated sj boundaries, as build_sj_arrays\n"
                 "emits them. A SECOND call because set_regions refuses to run twice.\n\n"
                 "offsets : int64[n_region_bounds_total + 1]\n"
                 "    CSR over the flat region_bound axis, keyed by the DONOR region_bound index.\n"
                 "boundary_right : int64[n_sj]\n"
                 "    Flat region_bound index of each sj's high end. The sj-boundary id\n"
                 "    IS the slot here; edges_df.edge_row is a join key and is not\n"
                 "    passed. Slot order is a contract: sort on\n"
                 "    (donor region_bound, acceptor region_bound, sj_strand).\n"
                 "sj_strand : int8[n_sj]\n"
                 "    Each sj's ANNOTATED strand.\n\n"
                 "Not optional: scan() refuses to run without it, because a missing\n"
                 "table makes every observed intron read as unannotated and empties\n"
                 "the sj boundaries silently. Pass empty arrays to declare none.\n")
        ;

    nb::class_<BamAnnotationWriter>(m, "BamAnnotationWriter")
        .def(nb::init<FragmentResolver&, const std::string&, bool, bool>(),
             nb::arg("ctx"),
             nb::arg("sj_tag_spec"),
             nb::arg("skip_duplicates"),
             nb::arg("include_multimap"),
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
             "    Include multimapping reads (must match pass 1).\n",
             nb::keep_alive<1, 2>()
        )
        .def("write", &BamAnnotationWriter::write,
             nb::arg("bam_path"),
             nb::arg("output_path"),
             nb::arg("ann_frag_ids"),
             nb::arg("ann_best_tid"),
             nb::arg("ann_best_gid"),
             nb::arg("ann_tx_flags"),
             nb::arg("ann_posterior"),
             nb::arg("ann_frag_class"),
             nb::arg("ann_n_candidates"),
             nb::arg("ann_splice_type"),
             nb::arg("ann_locus_id"),
             nb::arg("ann_size"),
             nb::arg("t_ids"),
             nb::arg("g_ids"),
             nb::arg("g_names"),
             "Stamp annotation tags and write BAM.\n\n"
             "Parameters\n"
             "----------\n"
             "bam_path : str\n"
             "    Path to name-sorted input BAM.\n"
             "output_path : str\n"
             "    Path for output annotated BAM.\n"
             "ann_frag_ids .. ann_splice_type : ndarray\n"
             "    Annotation arrays (sliced to ann_size).\n"
             "ann_size : int\n"
             "    Number of annotation entries.\n"
             "t_ids, g_ids, g_names : list[str]\n"
             "    Transcript / gene ID / gene name lookup arrays.\n\n"
             "Returns\n"
             "-------\n"
             "dict\n"
             "    Summary: n_read_groups, n_annotated, n_intergenic, etc.\n")
        ;

}
