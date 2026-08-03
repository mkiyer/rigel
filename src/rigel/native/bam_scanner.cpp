/**
 * bam_scanner.cpp — C++ BAM scanner using htslib for rigel.
 *
 * Replaces the Python BAM parsing hot path (bam.py parse_bam_file +
 * fragment.py Fragment.from_reads + resolution.py resolve_fragment +
 * buffer.append) with a single C++ scan that:
 *
 *   1. Reads BAM records via htslib (no pysam overhead)
 *   2. Groups records by query name (name-sorted BAM)
 *   3. Parses CIGAR → exon blocks + splice junctions
 *   4. Builds Fragment-equivalent structures in C++
 *   5. Calls FragmentResolver::_resolve_core() directly
 *   6. Appends results to FragmentAccumulator
 *   7. Collects model training observations as arrays
 *   8. Returns everything to Python in one call
 *
 * Module: rigel._bam_impl
 *
 * Build:
 *   Part of the rigel scikit-build-core build — see CMakeLists.txt.
 *   Requires htslib (linked via -lhts).
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <sys/stat.h>

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
static constexpr int CIG_INS      = BAM_CINS;       // 1 I
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
// SJ strand tag configuration
// ================================================================

enum class SJTagMode : uint8_t {
    NONE = 0,     // No SJ strand tag
    XS_ONLY = 1,  // Check XS tag only
    TS_ONLY = 2,  // Check ts tag only (minimap2)
    XS_TS = 3,    // Check XS first, then ts
    TS_XS = 4,    // Check ts first, then XS
};

// Parse SJ tag configuration string → SJTagMode enum.
static SJTagMode parse_sj_tag_spec(const std::string& spec) {
    if (spec == "XS")    return SJTagMode::XS_ONLY;
    if (spec == "ts")    return SJTagMode::TS_ONLY;
    if (spec == "XS,ts") return SJTagMode::XS_TS;
    if (spec == "ts,XS") return SJTagMode::TS_XS;
    if (spec == "none" || spec.empty()) return SJTagMode::NONE;
    return SJTagMode::XS_TS;  // default
}

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

// ================================================================
// ZF assignment-flag constants are defined in constants.h.  This file
// references them through the unqualified names via the using-decls
// below so the stamp sites read cleanly.
// ================================================================

using rigel::AF_UNRESOLVED;
using rigel::AF_MRNA;
using rigel::AF_NRNA;
using rigel::AF_NRNA_SYNTH;
using rigel::AF_GDNA_EM;
using rigel::AF_GDNA_INTERGENIC;
using rigel::AF_CHIMERIC;
using rigel::AF_MULTIMAPPER_DROP;

static const char* splice_type_label(int code) {
    switch (code) {
        case SPLICE_UNSPLICED:       return "unspliced";
        case SPLICE_SPLICED_UNANNOT: return "spliced_unannot";
        case SPLICE_SPLICED_ANNOT:   return "spliced_annot";
        case SPLICE_IMPLICIT:        return "spliced_implicit";
        case SPLICE_ARTIFACT:        return "splice_artifact";
        default:                     return "unknown";
    }
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
/// The anchor on each side of the junction is the number of
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
    int32_t ref_id;         // tid
    int32_t ref_start;      // pos
    int32_t mate_ref_id;    // mtid
    int32_t mate_ref_start; // mpos
    uint16_t flag;
    // Per-record count of CIGAR junctions dropped by the splice-artifact
    // blacklist.  Populated by filter_blacklisted_sjs in Pass 1 (for
    // scanner stats parity) and Pass 2 (for the ZB:i annotated-BAM tag).
    // Fits into the padding after `flag` — zero effective memory growth.
    uint16_t n_sj_blacklisted = 0;
    int32_t nm;             // NM tag (edit distance), 0 if absent
    int32_t nh;             // NH tag, 1 if absent
    int32_t hi;             // HI tag, -1 if absent
    int32_t sj_strand;      // from XS/ts tag, STRAND_NONE if absent

    // CIGAR-parsed exon blocks and splice junctions
    std::vector<std::pair<int32_t, int32_t>> exons;  // (start, end)
    std::vector<SJCigarEntry> sjs;

    bool is_read1()         const { return flag & BAM_FREAD1; }
    bool is_read2()         const { return flag & BAM_FREAD2; }
    bool is_secondary()     const { return flag & BAM_FSECONDARY; }
    bool is_supplementary() const { return flag & BAM_FSUPPLEMENTARY; }
    bool is_reverse()       const { return flag & BAM_FREVERSE; }
    bool is_proper_pair()   const { return flag & BAM_FPROPER_PAIR; }
    bool is_mate_unmapped() const { return flag & BAM_FMUNMAP; }
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

/// Per-junction sense / antisense counts.  ⚠ These ARE counts, unlike the
/// int8 label vectors below (whose count is the vector's size()), so they are
/// sized as counts: uint64 removes the overflow question permanently at a cost
/// of ~2.4 MB over a human junction set.
struct SJStrandCounts {
    uint64_t n_sense = 0;
    uint64_t n_antisense = 0;
};

struct StrandObservations {
    // exonic: (align_strand, transcript_strand) for unique-gene unambiguous strand
    std::vector<int8_t> exonic_obs;
    std::vector<int8_t> exonic_truth;

    // ⭐ The per-junction SJ strand table.
    // A junction is uniquely specified by (ref, start, end, motif strand); each
    // strand-qualified fragment credits its leftmost ANNOTATED junction with one
    // sense (align_strand == sj_strand) or antisense observation.
    //
    // ⭐ This is a strict REFINEMENT of the exonic_spliced 2×2 above, which is
    // exactly its marginal: summing n_sense over motif-POS junctions gives
    // pos_pos, n_antisense over motif-POS gives neg_pos, and symmetrically for
    // motif-NEG.  Populated from the SAME `get_is_strand_qualified()` branch, one
    // observation per fragment, so the identity is unconditional.  It exists
    // because a dispersion ACROSS junctions cannot be recovered from the 2×2.
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
    int64_t n_intergenic_unspliced = 0;
    int64_t n_intergenic_spliced = 0;
    int64_t n_with_exon = 0;
    int64_t n_with_annotated_sj = 0;
    int64_t n_with_unannotated_sj = 0;
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

    // ── the per-fragment splice census ────────────────────────────────────────────────────────────
    //
    // ⭐ ONE observation per fragment the scanner offers the accumulator, filed under the resolver's
    // own classification. This IS `rigel report`'s splice breakdown, which used to be read off the
    // fragment-length CATEGORY MODELS — so it counted only the fragments that also yielded a length
    // observation, a population nobody had stated. The scanner classifies, so the scanner counts.
    //
    // ⚠ An ARRAY and not five named fields, deliberately: `census[st]++` is one statement, the merge
    // is one loop, the export is one loop, and a category added to `SpliceType` cannot slip past any
    // of the three. The five names exist once, in `splice_type_label`.
    std::array<int64_t, NUM_SPLICE_TYPES> splice_census{};

    // Censused, then never offered to `deposit()`: the deposit adapter could not express the
    // fragment as ONE molecule on ONE cut axis (chiefly blocks on more than one reference), or no
    // accumulator was installed. ⭐ It exists to close the books — without it those `return`s are a
    // silent loss, which is the pattern this accumulator was rewritten to delete. With it:
    //
    //   Σ census − census[SPLICE_ARTIFACT] == qc.deposited + Σ qc.dropped_* + n_deposit_not_offered
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
        n_intergenic_unspliced += o.n_intergenic_unspliced;
        n_intergenic_spliced += o.n_intergenic_spliced;
        n_with_exon += o.n_with_exon;
        n_with_annotated_sj += o.n_with_annotated_sj;
        n_with_unannotated_sj += o.n_with_unannotated_sj;
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
        for (size_t i = 0; i < splice_census.size(); i++)
            splice_census[i] += o.splice_census[i];
        n_deposit_not_offered += o.n_deposit_not_offered;
    }
};

// ================================================================
// Fragment-equivalent structure (built in C++, no Python object)
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
    // Per-worker fractional accumulator (one Accumulator per ref). Owned
    // here so worker threads never touch the shared template. Null when no
    // region partition was provided.
    std::unique_ptr<rigel::accumulator::AccumulatorSet> acc_set;

    // Reusable scratch so the deposit path allocates no per-fragment memory. ⭐ Measured on the shipped
    // accumulator: one per-fragment std::vector cost 22.8 ns, 18 % of the deposit, and was invisible to a
    // sampling profiler because the time is attributed to `malloc`. Capacity amortizes to the high-water
    // mark and is never released.
    //
    // Two buffers, because they belong to two different owners. `deposit_scratch` is the accumulator's
    // (normalised introns, path segments, junction ids) and is opaque here. `deposit_introns` is the
    // ADAPTER's: the fragment's introns restricted to the reference being deposited and de-duplicated,
    // which is what `FragmentPath::introns` points at.
    //
    // ⛔ That restriction is mandatory, not tidiness. `Accumulator::deposit` normalises introns by
    // coordinate alone — it never looks at `IntronBlock::ref_id` — so an intron from another reference
    // would be cut out of this reference's path. `fragment_genomic_spans` filtered by ref explicitly; the
    // filter has to survive its deletion. Multi-reference fragments DO arrive here: the intergenic call
    // site is not chimera-gated, and `detect_chimera` returns CHIMERA_NONE when the blocks carry empty
    // transcript sets.
    rigel::accumulator::DepositScratch deposit_scratch;
    std::vector<IntronBlock>           deposit_introns;
    std::vector<rigel::accumulator::GapHypothesis> gap_hypotheses;

    WorkerState(int32_t n_transcripts, int64_t /*unused*/)
        : scratch(n_transcripts) {}
};

// Merge observations
static void merge_strand_obs(StrandObservations& dst, StrandObservations& src) {
    dst.exonic_obs.insert(dst.exonic_obs.end(),
                          src.exonic_obs.begin(), src.exonic_obs.end());
    dst.exonic_truth.insert(dst.exonic_truth.end(),
                            src.exonic_truth.begin(), src.exonic_truth.end());
    for (const auto& [key, counts] : src.sj_strand_table) {
        auto& dst_counts = dst.sj_strand_table[key];
        dst_counts.n_sense     += counts.n_sense;
        dst_counts.n_antisense += counts.n_antisense;
    }
    src.sj_strand_table.clear();
}

// ================================================================
// CIGAR parsing — ports parse_read() from bam.py
// ================================================================
//
// In addition to emitting exon blocks and splice junctions, we compute
// per-junction anchor lengths (left, right).  An "anchor" is the count
// of reference-advancing CIGAR bases on one side of the N op, bounded
// by clips/ends.  Insertions never contribute; deletions do (they
// occupy reference positions in the aligner's evidence for the anchor).
//
// Example: CIGAR 10M 500N 15M 1000N 125M
//   junction 1 (500N):  anchor_left = 10,  anchor_right = 15 + 125 = 140
//   junction 2 (1000N): anchor_left = 10 + 15 = 25, anchor_right = 125

static void parse_cigar(
    const bam1_t* b,
    int32_t ref_id,
    int32_t sj_strand,
    std::vector<std::pair<int32_t, int32_t>>& exons,
    std::vector<SJCigarEntry>& sjs)
{
    (void)ref_id;  // ref_id carried by the caller; kept for API stability
    exons.clear();
    sjs.clear();

    int32_t n_cigar = b->core.n_cigar;
    if (n_cigar == 0) return;

    const uint32_t* cigar = bam_get_cigar(b);
    int32_t pos = b->core.pos;     // 0-based genomic position
    int32_t start = pos;           // current exon start
    int32_t ref_advanced = 0;      // ref-advancing bases seen so far

    // First pass: build exons, record each junction's left anchor and
    // its CIGAR intron length so the right anchor can be back-filled.
    const size_t sj_base = sjs.size();
    for (int32_t i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        if (op == CIG_REF_SKIP) {
            if (pos > start) exons.push_back({start, pos});
            SJCigarEntry sj{};
            sj.start = pos;
            sj.end = pos + len;
            sj.strand = sj_strand;
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
    for (size_t k = sj_base; k < sjs.size(); k++) {
        sjs[k].anchor_right = ref_advanced - sjs[k].anchor_left;
    }
}

// ================================================================
// Splice-junction artifact filtering
// ================================================================
//
// Drop each CIGAR splice junction whose (left, right) anchors sit
// inside the blacklist's observed artifact envelope.  Rejection rule:
// either anchor <= max → junction is indistinguishable from known
// gDNA-derived artifacts and is treated as unspliced.
//
// Returns the number of junctions dropped (for scanner statistics).
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

static int32_t read_sj_strand(const bam1_t* b, SJTagMode mode) {
    if (mode == SJTagMode::NONE) return STRAND_NONE;

    auto try_tag = [&](const char* tag, bool is_ts) -> int32_t {
        uint8_t* aux = bam_aux_get(b, tag);
        if (!aux) return -1;  // tag not present
        // Dispatch on the htslib type code (*aux) to call the correct
        // accessor.  Type 'A' = single char, type 'Z'/'H' = string
        // (take first character).  Any other type → val 0 → STRAND_NONE.
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
        else return STRAND_NONE;

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
        default:
            break;
    }

    return (result >= 0) ? result : STRAND_NONE;
}

// ================================================================
// Parse a raw htslib bam1_t into a ParsedAlignment value type
// ================================================================
// Extracts all fields the pipeline needs (scalars, tags, CIGAR-derived
// exon blocks and splice junctions) so that downstream code never
// touches bam1_t directly.  Used by both single-threaded scan() and
// the producer thread in scan().

static ParsedAlignment parse_bam_record(
    const bam1_t* b,
    const std::vector<int32_t>& tid_to_ref_id,
    SJTagMode sj_tag_mode,
    const FragmentResolver* resolver,
    int64_t* n_sj_observed,
    int64_t* n_sj_blacklisted)
{
    ParsedAlignment rec;
    rec.flag = b->core.flag;
    rec.ref_id = b->core.tid;
    rec.ref_start = b->core.pos;
    rec.mate_ref_id = b->core.mtid;
    rec.mate_ref_start = b->core.mpos;

    // BAM aux tags
    rec.nm = 0;
    uint8_t* nm_aux = bam_aux_get(b, "NM");
    if (nm_aux) rec.nm = bam_aux2i(nm_aux);

    rec.nh = 1;
    uint8_t* nh_aux = bam_aux_get(b, "NH");
    if (nh_aux) rec.nh = bam_aux2i(nh_aux);

    rec.hi = -1;
    uint8_t* hi_aux = bam_aux_get(b, "HI");
    if (hi_aux) rec.hi = bam_aux2i(hi_aux);

    // CIGAR → exon blocks + splice junctions
    int32_t mapped_ref_id = (rec.ref_id >= 0 &&
        rec.ref_id < static_cast<int32_t>(tid_to_ref_id.size()))
        ? tid_to_ref_id[rec.ref_id] : -1;
    parse_cigar(b, mapped_ref_id, STRAND_NONE, rec.exons, rec.sjs);
    rec.sj_strand = STRAND_NONE;
    if (!rec.sjs.empty()) {
        rec.sj_strand = read_sj_strand(b, sj_tag_mode);
        for (auto& sj : rec.sjs) sj.strand = rec.sj_strand;
    }
    rec.ref_id = mapped_ref_id;

    // Splice-junction artifact filtering.  Per-read, immediately after
    // CIGAR parsing.  Downstream fragment assembly unions across mates,
    // so this is the permissive rule: a junction survives as long as at
    // least one read supports it with anchors above the blacklist envelope.
    if (resolver && resolver->has_sj_blacklist() && !rec.sjs.empty()) {
        if (n_sj_observed) *n_sj_observed += static_cast<int64_t>(rec.sjs.size());
        int32_t n_dropped = filter_blacklisted_sjs(rec.sjs, *resolver, mapped_ref_id);
        rec.n_sj_blacklisted = static_cast<uint16_t>(n_dropped);
        if (n_sj_blacklisted) *n_sj_blacklisted += n_dropped;
    } else if (n_sj_observed) {
        *n_sj_observed += static_cast<int64_t>(rec.sjs.size());
    }

    return rec;
}

// ================================================================
// Build AssembledFragment from a hit (r1_records, r2_records)
// ================================================================

static AssembledFragment build_fragment(
    const std::vector<ParsedAlignment*>& r1_reads,
    const std::vector<ParsedAlignment*>& r2_reads)
{
    // Key: (ref_id, strand) → list of exon intervals
    // We use a simple map with int64_t key = (ref_id << 32) | strand
    std::unordered_map<int64_t, std::vector<std::pair<int32_t, int32_t>>> exon_dict;
    std::set<std::tuple<int32_t, int32_t, int32_t, int32_t>> intron_set;
    // (ref_id, start, end, strand)

    auto make_key = [](int32_t ref_id, int32_t strand) -> int64_t {
        return (static_cast<int64_t>(ref_id) << 32) |
               (static_cast<uint32_t>(strand));
    };

    int32_t nm_total = 0;

    // Process R1 reads
    for (const auto* rec : r1_reads) {
        int32_t ref_strand = rec->is_reverse() ? STRAND_NEG : STRAND_POS;
        int64_t key = make_key(rec->ref_id, ref_strand);
        for (const auto& [s, e] : rec->exons) {
            exon_dict[key].push_back({s, e});
        }
        for (const auto& sj : rec->sjs) {
            intron_set.insert({rec->ref_id, sj.start, sj.end, sj.strand});
        }
        nm_total += rec->nm;
    }

    // Process R2 reads — flip strand (R2 strand flip convention)
    for (const auto* rec : r2_reads) {
        int32_t ref_strand = rec->is_reverse() ? STRAND_NEG : STRAND_POS;
        // Flip: POS→NEG, NEG→POS
        ref_strand = (ref_strand == STRAND_POS) ? STRAND_NEG : STRAND_POS;
        int64_t key = make_key(rec->ref_id, ref_strand);
        for (const auto& [s, e] : rec->exons) {
            exon_dict[key].push_back({s, e});
        }
        for (const auto& sj : rec->sjs) {
            intron_set.insert({rec->ref_id, sj.start, sj.end, sj.strand});
        }
        nm_total += rec->nm;
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

    // Sort intron blocks
    std::vector<IntronBlock> sorted_introns;
    sorted_introns.reserve(intron_set.size());
    for (const auto& [ref, s, e, st] : intron_set) {
        sorted_introns.push_back({ref, s, e, st});
    }
    // Already sorted since std::set is ordered

    AssembledFragment frag;
    frag.exons = std::move(merged_exons);
    frag.introns = std::move(sorted_introns);
    frag.nm = nm_total;
    return frag;
}

// ================================================================
// Hit grouping — ports _group_records_by_hit from bam.py
// ================================================================

struct AlignmentGroup {
    // Hits: each is (r1_reads, r2_reads)
    std::vector<std::pair<std::vector<ParsedAlignment*>,
                          std::vector<ParsedAlignment*>>> hits;
    // Secondary locations for transcript-aware pairing
    std::vector<std::vector<ParsedAlignment*>> sec_r1_locs;
    std::vector<std::vector<ParsedAlignment*>> sec_r2_locs;
};

static AlignmentGroup group_records_by_hit(
    std::vector<ParsedAlignment>& usable)
{
    AlignmentGroup result;

    // Detect HI tags
    bool has_hi = false;
    for (const auto& r : usable) {
        if (r.hi >= 0) { has_hi = true; break; }
    }

    if (has_hi) {
        // Group by HI value
        std::unordered_map<int32_t,
            std::pair<std::vector<ParsedAlignment*>,
                      std::vector<ParsedAlignment*>>> groups;

        for (auto& r : usable) {
            int32_t hi = (r.hi >= 0) ? r.hi : 0;
            auto& [r1s, r2s] = groups[hi];
            if (r.is_read1())
                r1s.push_back(&r);
            else if (r.is_read2())
                r2s.push_back(&r);
        }

        // Sort by HI for determinism
        std::vector<int32_t> hi_keys;
        hi_keys.reserve(groups.size());
        for (auto& [hi, _] : groups) hi_keys.push_back(hi);
        std::sort(hi_keys.begin(), hi_keys.end());

        for (int32_t hi : hi_keys) {
            result.hits.push_back(std::move(groups[hi]));
        }
        return result;
    }

    // No HI tags — separate primary from secondary
    std::vector<ParsedAlignment*> primary_r1, primary_r2;
    std::vector<std::pair<std::vector<ParsedAlignment*>,
                          std::vector<ParsedAlignment*>>> singleton_hits;

    for (auto& r : usable) {
        if (r.is_supplementary()) {
            if (r.is_read1())
                primary_r1.push_back(&r);
            else
                primary_r2.push_back(&r);
        } else if (r.is_secondary()) {
            if (r.is_mate_unmapped()) {
                if (r.is_read1())
                    singleton_hits.push_back({{&r}, {}});
                else
                    singleton_hits.push_back({{}, {&r}});
            } else if (r.is_read1()) {
                result.sec_r1_locs.push_back({&r});
            } else {
                result.sec_r2_locs.push_back({&r});
            }
        } else {
            // Primary
            if (r.is_read1())
                primary_r1.push_back(&r);
            else
                primary_r2.push_back(&r);
        }
    }

    bool has_secondary = !result.sec_r1_locs.empty() ||
                         !result.sec_r2_locs.empty();

    if (has_secondary) {
        // Multimapper without HI: include primary as additional locations
        if (!primary_r1.empty())
            result.sec_r1_locs.insert(result.sec_r1_locs.begin(), primary_r1);
        if (!primary_r2.empty())
            result.sec_r2_locs.insert(result.sec_r2_locs.begin(), primary_r2);
        result.hits = std::move(singleton_hits);
        return result;
    }

    // Unique mapper: primary pair as sole hit
    if (!primary_r1.empty() || !primary_r2.empty()) {
        result.hits.push_back({std::move(primary_r1),
                               std::move(primary_r2)});
    }
    result.hits.insert(result.hits.end(),
                       singleton_hits.begin(), singleton_hits.end());
    return result;
}

// ================================================================
// Multimapper pairing — ports pair_multimapper_reads from resolution.py
// ================================================================
// Thread-safe: uses caller-supplied ResolverScratch for all resolve calls.

static std::vector<std::pair<std::vector<ParsedAlignment*>,
                              std::vector<ParsedAlignment*>>>
pair_multimapper_reads(
    std::vector<std::vector<ParsedAlignment*>>& sec_r1_locs,
    std::vector<std::vector<ParsedAlignment*>>& sec_r2_locs,
    FragmentResolver& ctx,
    ResolverScratch& scratch)
{
    using Hit = std::pair<std::vector<ParsedAlignment*>,
                          std::vector<ParsedAlignment*>>;
    std::vector<Hit> paired;

    if (sec_r1_locs.empty() && sec_r2_locs.empty())
        return paired;

    struct ResolvedLoc {
        std::vector<ParsedAlignment*>* reads;
        std::vector<int32_t> t_inds;
        int32_t ref_id;
        int32_t ref_start;
    };

    auto resolve_location = [&](std::vector<ParsedAlignment*>& reads,
                                bool is_r2) -> std::vector<int32_t>
    {
        AssembledFragment frag = is_r2
            ? build_fragment({}, reads)
            : build_fragment(reads, {});
        std::vector<int32_t> t_inds;
        if (!frag.exons.empty()) {
            RawResolveResult cr;
            if (ctx._resolve_core(frag.exons, frag.introns,
                                   frag.genomic_footprint(), cr, scratch)) {
                t_inds = std::move(cr.t_inds);
            }
        }
        return t_inds;
    };

    std::vector<ResolvedLoc> r1_resolved;
    r1_resolved.reserve(sec_r1_locs.size());
    for (auto& r1_reads : sec_r1_locs) {
        auto t_inds = resolve_location(r1_reads, false);
        int32_t ref_id = r1_reads.empty() ? -1 : r1_reads[0]->ref_id;
        int32_t ref_start = r1_reads.empty() ? -1 : r1_reads[0]->ref_start;
        r1_resolved.push_back({&r1_reads, std::move(t_inds), ref_id, ref_start});
    }

    std::vector<ResolvedLoc> r2_resolved;
    r2_resolved.reserve(sec_r2_locs.size());
    for (auto& r2_reads : sec_r2_locs) {
        auto t_inds = resolve_location(r2_reads, true);
        int32_t ref_id = r2_reads.empty() ? -1 : r2_reads[0]->ref_id;
        int32_t ref_start = r2_reads.empty() ? -1 : r2_reads[0]->ref_start;
        r2_resolved.push_back({&r2_reads, std::move(t_inds), ref_id, ref_start});
    }

    // STRICT — pair by transcript-set intersection
    std::unordered_set<int> r1_paired, r2_paired;

    for (int i = 0; i < static_cast<int>(r1_resolved.size()); i++) {
        if (r1_resolved[i].t_inds.empty()) continue;
        for (int j = 0; j < static_cast<int>(r2_resolved.size()); j++) {
            if (r2_resolved[j].t_inds.empty()) continue;
            if (has_intersection(r1_resolved[i].t_inds,
                                 r2_resolved[j].t_inds)) {
                paired.push_back({*r1_resolved[i].reads,
                                  *r2_resolved[j].reads});
                r1_paired.insert(i);
                r2_paired.insert(j);
            }
        }
    }

    // FALLBACK — same-reference closest distance
    std::vector<int> unmatched_r1, unmatched_r2;
    for (int i = 0; i < static_cast<int>(r1_resolved.size()); i++)
        if (!r1_paired.count(i)) unmatched_r1.push_back(i);
    for (int j = 0; j < static_cast<int>(r2_resolved.size()); j++)
        if (!r2_paired.count(j)) unmatched_r2.push_back(j);

    if (!unmatched_r1.empty() && !unmatched_r2.empty()) {
        std::vector<std::tuple<int32_t, int, int>> candidates;
        for (int i : unmatched_r1) {
            for (int j : unmatched_r2) {
                if (r1_resolved[i].ref_id == r2_resolved[j].ref_id &&
                    r1_resolved[i].ref_id >= 0) {
                    int32_t dist = std::abs(r1_resolved[i].ref_start -
                                            r2_resolved[j].ref_start);
                    candidates.push_back({dist, i, j});
                }
            }
        }
        std::sort(candidates.begin(), candidates.end());
        for (const auto& [dist, i, j] : candidates) {
            if (!r1_paired.count(i) && !r2_paired.count(j)) {
                paired.push_back({*r1_resolved[i].reads,
                                  *r2_resolved[j].reads});
                r1_paired.insert(i);
                r2_paired.insert(j);
            }
        }
    }

    // CROSS-PAIR remaining unmatched R1 × R2 combinatorially.
    // This produces all possible pairings.  Some may be chimeric
    // (e.g. R1 chr12 × R2 chrX) — those are detected downstream
    // by _resolve_core and skipped during scoring.  The non-chimeric
    // combinations (e.g. R1 chr12 × R2 chr12) are valid PE hits
    // that preserve fragment-length information and give the EM
    // better resolution than singletons would.
    std::vector<int> final_r1, final_r2;
    for (int i = 0; i < static_cast<int>(r1_resolved.size()); i++)
        if (!r1_paired.count(i)) final_r1.push_back(i);
    for (int j = 0; j < static_cast<int>(r2_resolved.size()); j++)
        if (!r2_paired.count(j)) final_r2.push_back(j);

    if (!final_r1.empty() && !final_r2.empty()) {
        for (int i : final_r1) {
            for (int j : final_r2) {
                paired.push_back({*r1_resolved[i].reads,
                                  *r2_resolved[j].reads});
                r1_paired.insert(i);
                r2_paired.insert(j);
            }
        }
    }

    // SINGLETONS — only if one side has reads but the other doesn't.
    for (int i = 0; i < static_cast<int>(r1_resolved.size()); i++)
        if (!r1_paired.count(i))
            paired.push_back({*r1_resolved[i].reads, {}});
    for (int j = 0; j < static_cast<int>(r2_resolved.size()); j++)
        if (!r2_paired.count(j))
            paired.push_back({{}, *r2_resolved[j].reads});

    return paired;
}

/// One `OfferedFragment` marshalled out of Python, with the storage its spans point into.
///
/// ⭐ Shared by `Accumulator.deposit` and `Accumulator.length_under` so the two cannot disagree about how
/// a hypothesis set crosses the ABI. ⚠ Reads attributes off whatever it is handed, so the parity gate can
/// pass the SAME `GapHypothesis` objects to the specification and to the binding; a binding with its own
/// tuple convention would be a second representation to keep in step.
struct MarshalledFragment {
    std::vector<IntronBlock> observed;
    std::vector<IntronBlock> implied;
    std::vector<int32_t>     supporting;
    std::vector<rigel::accumulator::GapHypothesis> spans;
    rigel::accumulator::OfferedFragment offered{};

    MarshalledFragment(int64_t start, int64_t end, nb::iterable introns,
                       int32_t align_strand, int32_t sj_strand, nb::object hypotheses) {
        for (nb::handle item : introns) {
            auto pair = nb::cast<nb::tuple>(item);
            observed.push_back({0, nb::cast<int32_t>(pair[0]), nb::cast<int32_t>(pair[1]), 0});
        }
        // ⛔ Reserved up front: the GapHypothesis spans below point INTO these vectors, so a reallocation
        // while filling them would dangle every pointer already handed out.
        std::size_t n_implied = 0, n_supporting = 0;
        for (nb::handle h : hypotheses) {
            n_implied += nb::len(nb::getattr(h, "introns"));
            n_supporting += nb::len(nb::getattr(h, "supporting_t_inds"));
        }
        implied.reserve(n_implied);
        supporting.reserve(n_supporting);

        for (nb::handle h : hypotheses) {
            const std::size_t i0 = implied.size(), t0 = supporting.size();
            for (nb::handle item : nb::getattr(h, "introns")) {
                auto pair = nb::cast<nb::tuple>(item);
                implied.push_back({0, nb::cast<int32_t>(pair[0]), nb::cast<int32_t>(pair[1]), 0});
            }
            for (nb::handle t : nb::getattr(h, "supporting_t_inds"))
                supporting.push_back(nb::cast<int32_t>(t));
            spans.push_back({nullptr, implied.size() - i0,
                             nb::cast<int32_t>(nb::getattr(h, "sj_strand")),
                             nullptr, supporting.size() - t0});
            spans.back().introns = reinterpret_cast<const IntronBlock*>(i0);
            spans.back().supporting_t = reinterpret_cast<const int32_t*>(t0);
        }
        for (auto& span : spans) {  // offsets -> pointers, now that nothing more will grow
            span.introns = implied.data() + reinterpret_cast<std::size_t>(span.introns);
            span.supporting_t = supporting.data() + reinterpret_cast<std::size_t>(span.supporting_t);
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
};

// ================================================================
// the two non-array banks, as Python dicts
// ================================================================
//
// ⭐ Written ONCE and used by both the bound `Accumulator` (the parity surface) and `build_result` (the
// payload). The specification's key strings live here and nowhere else on this side of the ABI, so the two
// exports cannot disagree about a name — which is the whole reason the payload, the reference and the
// parity gate share one vocabulary.

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
    put("hypothesis_t_offsets", deferred.hypothesis_t_offsets);
    put("hypothesis_t", deferred.hypothesis_t);
    return out;
}

// ================================================================
// BamScanner — main scanning class
// ================================================================

class BamScanner {
public:
    FragmentResolver* ctx_;
    SJTagMode sj_tag_mode_ = SJTagMode::XS_ONLY;
    bool skip_duplicates_ = true;
    bool include_multimap_ = false;

    // Mapping from htslib tid → internal ref_id (FragmentResolver's numbering)
    std::vector<int32_t> tid_to_ref_id_;

    // Results
    BamScanStats stats_;
    StrandObservations strand_obs_;

    // The accumulator's node partition, installed by set_regions() and the junction CSR by
    // set_junctions(). One Accumulator per BAM reference; deposits go into a per-worker copy during
    // scan() and are merged back here.
    //
    // ⚠ These are MEMBERS because the per-worker sets are built from them inside scan(). Anything
    // installed after that point is invisible to the workers, which is a silent half-empty tally.
    std::vector<int64_t> cut_positions_;
    std::vector<int64_t> ref_cut_offsets_;
    std::vector<uint8_t> node_types_;      // coarse node type, ref-major, one per node
    int                  max_length_ = 0;  // the fragment-length limit AND the pool-histogram width
    std::vector<int64_t> sj_offsets_;       // n_cuts + 1, CSR over the flat cut axis
    std::vector<int64_t> sj_acceptor_cut_;  // flat cut index of each junction's high end
    std::vector<int8_t>  sj_strand_;        // each junction's ANNOTATED strand
    bool                 junctions_set_ = false;
    std::unique_ptr<rigel::accumulator::AccumulatorSet> acc_set_;

    BamScanner(FragmentResolver& ctx,
               const std::string& sj_tag_spec,
               bool skip_duplicates,
               bool include_multimap)
        : ctx_(&ctx),
          skip_duplicates_(skip_duplicates),
          include_multimap_(include_multimap)
    {
        sj_tag_mode_ = parse_sj_tag_spec(sj_tag_spec);
    }

    // ----------------------------------------------------------------
    // The accumulator's partition
    // ----------------------------------------------------------------
    //
    // Install the node partition the accumulator deposits onto. Must be called BEFORE scan(); without it
    // the calibration payload comes back empty. `set_junctions` is a SECOND call, because this one
    // refuses to run twice.
    //
    // Inputs:
    //   cut_positions:   flat int64[n_cuts_total], the concatenated sorted cut positions of every
    //                    reference. A reference contributing c cuts owns c-1 nodes and c-2 lines.
    //   ref_cut_offsets: int64[n_refs + 1] offsets into cut_positions. A reference with
    //                    offsets[f+1] == offsets[f] has no partition and accepts no deposits.
    //   n_refs:          number of references (must match ctx_->ref_to_id_).
    //   node_types:      uint8[n_nodes_total], the coarse node type, ref-major; it types the length pools.
    //   max_length:      the fragment-length limit applied to L, and the pool-histogram width.
    void set_regions(
        nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> cut_positions,
        nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> ref_cut_offsets,
        int32_t n_refs,
        nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig> node_types,
        int max_length)
    {
        if (n_refs < 0) {
            throw std::invalid_argument(
                "set_regions: n_refs must be non-negative");
        }
        const std::size_t n_off = ref_cut_offsets.shape(0);
        if (n_off != static_cast<std::size_t>(n_refs) + 1) {
            throw std::invalid_argument(
                "set_regions: ref_cut_offsets must have length n_refs + 1");
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
        cut_positions_.assign(
            cut_positions.data(),
            cut_positions.data() + cut_positions.shape(0));
        ref_cut_offsets_.assign(
            ref_cut_offsets.data(),
            ref_cut_offsets.data() + n_off);
        node_types_.assign(
            node_types.data(),
            node_types.data() + node_types.shape(0));
        max_length_ = max_length;
        acc_set_ = make_accumulator_set();
    }

    // ----------------------------------------------------------------
    // The accumulator's junction edges
    // ----------------------------------------------------------------
    //
    // A SECOND method rather than two more arguments to set_regions, because that one throws if called
    // twice. Takes the flat CSR `build_junction_edge_arrays` emits, keyed by the flat cut index;
    // AccumulatorSet slices it per reference.
    //
    // ⚠ Required even when there are none: pass empty arrays to say "this annotation has no junctions".
    // scan() refuses to run without it, because a missing table is not an error anywhere — every observed
    // intron simply reads as unannotated, and 404,168 annotated junctions become zero silently.
    void set_junctions(
        nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> offsets,
        nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> acceptor_cut,
        nb::ndarray<const int8_t, nb::ndim<1>, nb::c_contig> sj_strand)
    {
        if (!acc_set_) {
            throw std::runtime_error(
                "set_junctions: call set_regions first — the junction CSR is keyed by cut index, so "
                "there is nothing to key it against until the partition is installed");
        }
        if (acceptor_cut.shape(0) != sj_strand.shape(0)) {
            throw std::invalid_argument(
                "set_junctions: acceptor_cut has " + std::to_string(acceptor_cut.shape(0)) +
                " entries but sj_strand has " + std::to_string(sj_strand.shape(0)));
        }
        sj_offsets_.assign(offsets.data(), offsets.data() + offsets.shape(0));
        sj_acceptor_cut_.assign(
            acceptor_cut.data(), acceptor_cut.data() + acceptor_cut.shape(0));
        sj_strand_.assign(sj_strand.data(), sj_strand.data() + sj_strand.shape(0));
        junctions_set_ = true;
        install_junctions(*acc_set_);
    }

    // Build a per-worker set over the same partition, locally writable so workers never contend.
    std::unique_ptr<rigel::accumulator::AccumulatorSet> make_accumulator_set() const {
        auto set = std::make_unique<rigel::accumulator::AccumulatorSet>(
            cut_positions_.data(),
            cut_positions_.size(),
            ref_cut_offsets_.data(),
            ref_cut_offsets_.empty() ? 0 : ref_cut_offsets_.size() - 1,
            node_types_.empty() ? nullptr : node_types_.data(),
            node_types_.size(),
            max_length_);
        install_junctions(*set);
        return set;
    }

    void install_junctions(rigel::accumulator::AccumulatorSet& set) const {
        if (!junctions_set_) return;
        set.set_junctions(sj_offsets_.data(),
                          sj_offsets_.size(),
                          sj_acceptor_cut_.data(),
                          sj_strand_.data(),
                          sj_acceptor_cut_.size(),
                          ref_cut_offsets_.data());
    }

    // ----------------------------------------------------------------
    // Main scan entry point — streaming chunk architecture
    // ----------------------------------------------------------------
    //
    // Three-thread architecture:
    //   Reader thread (background): BAM I/O → input_queue
    //   Worker threads (background): input_queue → accumulate → output_queue
    //   Main thread: output_queue → finalize_zero_copy → chunk_callback
    //
    // The main thread holds the GIL while calling the Python callback
    // and releases it while waiting on the output queue.

    nb::dict scan(const std::string& bam_path,
                  nb::callable chunk_callback,
                  const std::vector<int32_t>& t_strand_arr,
                  int64_t chunk_size = 1000000,
                  int n_workers = 1,
                  int n_decomp_threads = 2,
                  int qname_batch_size = 512)
    {
        if (n_workers < 1) n_workers = 1;
        if (chunk_size < 1) chunk_size = 1;
        if (qname_batch_size < 1) qname_batch_size = 1;

        // ⚠ Explicit rather than defaulted, because the failure is invisible: with no junction table every
        // observed intron reads as unannotated, so the spliced banks and every junction edge stay at zero
        // and the tally still looks well-formed. Pass empty arrays to declare an annotation with none.
        if (acc_set_ && !junctions_set_) {
            throw std::runtime_error(
                "scan: set_regions was called but set_junctions was not. The junction table is not "
                "optional — without it every observed intron reads as unannotated, so the junction edges "
                "and the spliced banks are silently empty. Call set_junctions with empty arrays if this "
                "annotation genuinely has no junctions.");
        }

        // Two queues: input (SPMC) and output (MPSC)
        BoundedQueue<QnameBatch> input_queue(
            static_cast<size_t>(n_workers * 4));
        BoundedQueue<FragmentAccumulator> output_queue(
            static_cast<size_t>(n_workers * 2));

        // Per-worker state (local to scan — not a class member)
        std::vector<std::unique_ptr<WorkerState>> worker_states;
        worker_states.reserve(n_workers);
        for (int i = 0; i < n_workers; i++) {
            int32_t n_transcripts = ctx_->n_transcripts_;
            auto ws = std::make_unique<WorkerState>(n_transcripts, 0);
            // Pre-allocate accumulator for chunk_size
            ws->accumulator.reserve(chunk_size, chunk_size * 3 / 2);
            // Per-worker accumulator: the same partition and the same junction CSR as the shared
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
        // push full chunks to output_queue.
        std::vector<std::thread> workers;
        workers.reserve(n_workers);
        for (int i = 0; i < n_workers; i++) {
            workers.emplace_back([&input_queue, &output_queue,
                                  &worker_states, ctx, i,
                                  include_multimap, chunk_size]()
            {
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
                            ws.accumulator.reserve(chunk_size, chunk_size * 3 / 2);
                        }
                    }
                    batch.groups.clear();
                    if (output_aborted) break;
                }
                // Flush remaining fragments
                if (!output_aborted && ws.accumulator.get_size() > 0) {
                    output_queue.push(std::move(ws.accumulator));
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

                while (sam_read1(fp, hdr, b) >= 0) {
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
                        parse_bam_record(b, tid_to_ref_id_, sj_tag_mode_,
                                          ctx_,
                                          &stats_.n_sj_observed,
                                          &stats_.n_sj_blacklisted));
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
                    // Acquire GIL for finalize_zero_copy + Python callback
                    nb::gil_scoped_acquire acquire;
                    nb::dict chunk = acc.finalize_zero_copy(t_strand_arr);
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

        // Propagate exceptions: prefer reader exception (root cause)
        if (reader_exception) {
            std::rethrow_exception(reader_exception);
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
        // Turn one assembled fragment into the accumulator's `FragmentPath` and deposit it. The
        // accumulator owns the whole deposit rule; this function's only job is to say what the fragment IS
        // — its extent on one reference, the introns cut out of it, and the two independent strands.
        //
        // ⭐ THE TWO STRANDS ARE INDEPENDENT, and collapsing them is the bug this rewrite deletes.
        //   align_strand  where the read ALIGNED. Every read has one. It selects the array column.
        //   sj_strand     the splice junction's strand. Spliced reads only. It resolves an intron against
        //                 the annotation, and nothing else.
        // The shipped code compared them into a third concept, a bool named `primary`, and used it to pick
        // a channel labelled *sense* — which is how a dUTP first-strand library ended up with 0.6 % of its
        // spliced fragments in that column. `primary` is deleted, not renamed, and nothing replaces it:
        // sense/antisense is DERIVED by a consumer from the fragment strand and the junction's own strand.
        //
        // ⚠ No strand gate here either. The deposit rejects an undefined `align_strand` itself and COUNTS
        // it, which is the point — the old `align_ok`/`motif_ok` gate returned early, so the loss vanished.
        const auto deposit_to_accumulator =
            [&ws](const AssembledFragment& f, const RawResolveResult& cr) {
                const int32_t st = cr.splice_type;

                // ── the splice census ─────────────────────────────────────────────────────────────
                //
                // ⭐ FIRST, and before any gate, because this is the count of what the scanner SAW.
                // One observation per fragment reaching this adapter — unique mapper, resolved,
                // non-chimeric — which is exactly the population the accumulator is offered, so the
                // report's splice breakdown and its fragment-length histograms describe the same
                // fragments. That was never true of the category models this replaces.
                //
                // ⚠ It is deliberately NOT gated on `ws.acc_set`: what the scanner classified does
                // not depend on whether an accumulator happened to be installed.
                if (st >= 0 && static_cast<size_t>(st) < ws.stats.splice_census.size())
                    ws.stats.splice_census[st]++;

                // Hold artifact splices (blacklisted CIGAR-N) out of the accumulator entirely — no
                // deposit, no length pool. Their true span is unrecoverable: a blacklisted junction may be
                // a real-but-rejected junction OR a wholly incorrect alignment, and the span would be
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
                // coordinate alone — it never looks at `ref_id`.
                //
                // ⚠ The de-duplication is for the QC denominator, NOT to prevent double-crediting a
                // junction: `parse_bam_record` reads XS/ts once per RECORD and `build_fragment` keys
                // `intron_set` on (ref, start, end, strand), so a pair where one mate carries the tag
                // yields the same intron twice. Each duplicate would increment `introns_absorbed`, on
                // every spliced pair from such an aligner, and that counter is reported. `intron_set` is
                // ordered, so consecutive comparison is enough; ⛔ it must NOT be re-keyed or merged in
                // `build_fragment`, which shares it with the resolver (measured: merging demotes
                // SPLICED_ANNOT to SPLICED_UNANNOT and widens `t_inds`).
                auto& observed = ws.deposit_introns;
                observed.clear();
                for (const auto& intron : f.introns) {
                    if (intron.ref_id != ref_id) continue;
                    if (!observed.empty() && observed.back().start == intron.start &&
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
                auto& hypotheses = ws.gap_hypotheses;
                hypotheses.clear();
                for (std::int32_t h = 0; h < cr.n_gap_hypotheses(); ++h) {
                    const std::int32_t i0 = cr.gap_intron_offsets[h];
                    const std::int32_t t0 = cr.gap_supporting_offsets[h];
                    hypotheses.push_back({
                        cr.gap_introns.data() + i0,
                        static_cast<std::size_t>(cr.gap_intron_offsets[h + 1] - i0),
                        cr.gap_sj_strand[h],
                        cr.gap_supporting.data() + t0,
                        static_cast<std::size_t>(cr.gap_supporting_offsets[h + 1] - t0),
                    });
                }

                // The molecule's extent on this reference: leftmost block start to rightmost block end,
                // MATE GAP INCLUDED, because the gap is part of the molecule and must count toward L.
                //
                // ⛔ A fragment with blocks on MORE THAN ONE REFERENCE deposits nothing. It is not one
                // molecule (design §3.3), and an `OfferedFragment` cannot express it — it carries one
                // extent on one cut axis. The shipped code had no such check on the intergenic path: it
                // computed a span per reference and deposited ALL of them onto `exons.front().ref_id`, so
                // chr7 coordinates landed on chr1's cut axis. `ws.span_ref` recorded which reference each
                // span belonged to and **nothing ever read it**, which is how that survived.
                //
                // ⚠ Deliberately narrow: this tests multi-reference, NOT `cr.chimera_type`. That field is
                // also set for single-reference *cis* chimeras, which the intergenic path deposits today,
                // and stopping those is a change to WHAT COUNTS AS A FRAGMENT — its own arm with its own
                // before/after measurement, per the plan's TODO.
                std::int64_t start = 0, end = 0;
                bool any = false;
                for (const auto& block : f.exons) {
                    if (block.ref_id != ref_id) { ws.stats.n_deposit_not_offered++; return; }
                    if (!any) {
                        start = block.start;
                        end = block.end;
                        any = true;
                    } else {
                        start = std::min<std::int64_t>(start, block.start);
                        end = std::max<std::int64_t>(end, block.end);
                    }
                }
                if (!any) { ws.stats.n_deposit_not_offered++; return; }

                rigel::accumulator::OfferedFragment offered;
                offered.start = start;
                offered.end = end;
                offered.observed_introns = observed.data();
                offered.n_observed_introns = observed.size();
                offered.align_strand = cr.align_strand;
                // ⚠ The OBSERVED motif strand, straight from `cr`. A hypothesis carries the strand its
                // supporting transcripts imply, and `deposit` falls back to that only when nothing was
                // sequenced. An observed motif is evidence; an implied strand is an inference from the
                // annotation, and mixing an inference into an observation is how `primary` went wrong.
                offered.sj_strand = cr.sj_strand;
                offered.hypotheses = hypotheses.data();
                offered.n_hypotheses = hypotheses.size();

                ws.acc_set->at(ref_id).deposit(offered, ws.deposit_scratch);
            };

        // Per-worker state refs
        stats.n_read_names++;

        int32_t nh = 1;
        for (const auto& r : records) {
            if (r.nh > 1) { nh = r.nh; break; }
        }

        bool has_secondary = false;
        for (const auto& r : records) {
            if (r.is_secondary()) { has_secondary = true; break; }
        }
        bool is_multimap = (nh > 1) || has_secondary;

        if (is_multimap) {
            stats.multimapping++;
            if (!include_multimap) return;
        } else {
            stats.unique++;
        }

        auto hit_result = group_records_by_hit(records);

        // Track hit-level mate stats
        for (const auto& [r1s, r2s] : hit_result.hits) {
            if (r1s.empty() || r2s.empty()) {
                stats.mate_unmapped++;
            } else {
                bool found_proper = false;
                for (const auto* r : r1s) {
                    if (!r->is_supplementary() && !r->is_secondary()) {
                        found_proper = r->is_proper_pair();
                        break;
                    }
                }
                if (found_proper) stats.proper_pair++;
                else stats.improper_pair++;
            }
        }

        // Multimapper secondary pairing using per-worker scratch
        std::vector<std::pair<std::vector<ParsedAlignment*>,
                              std::vector<ParsedAlignment*>>> secondary_pairs;
        if (!hit_result.sec_r1_locs.empty() ||
            !hit_result.sec_r2_locs.empty()) {
            // Thread-safe multimapper pairing using per-worker scratch
            secondary_pairs = pair_multimapper_reads(
                hit_result.sec_r1_locs,
                hit_result.sec_r2_locs,
                ctx, scratch);
        }

        auto& all_hits = hit_result.hits;
        all_hits.insert(all_hits.end(),
                        secondary_pairs.begin(), secondary_pairs.end());

        int32_t num_hits = std::max(nh, static_cast<int32_t>(all_hits.size()));
        bool is_unique_mapper = (num_hits == 1);
        int32_t n_buffered_mm = 0;

        // Track per-molecule intergenic status across all hits.
        // For multimappers: count intergenic only if NO hit resolves,
        // and count it exactly once (not per-hit).
        bool any_hit_resolved = false;
        bool any_non_chimeric_resolved = false;
        bool any_unresolved_spliced = false;
        bool any_unresolved_unspliced = false;

        // Track chimera per-fragment (not per-hit).
        // A fragment is chimeric only if it resolves but ALL
        // resolved hits are chimeric.
        bool any_hit_chimeric = false;
        int32_t worst_chimera_type = CHIMERA_NONE;

        // Count ONE fragment per physical molecule (not per hit).
        // (Phase A burndown 2026-05-29: per-fragment calibration
        // observation capture removed; rebuilt against the new
        // fractional-accumulator binding in Phase B — see
        stats.n_fragments++;

        for (const auto& [r1_reads, r2_reads] : all_hits) {
            AssembledFragment frag = build_fragment(r1_reads, r2_reads);

            if (frag.exons.empty()) continue;

            // SRD v2: aggregate per-record blacklist counts into a
            // per-fragment total so the resolver can promote
            // SPLICE_UNSPLICED -> SPLICE_ARTIFACT.
            int32_t frag_n_sj_blacklisted = 0;
            for (const auto* r : r1_reads)
                frag_n_sj_blacklisted += r->n_sj_blacklisted;
            for (const auto* r : r2_reads)
                frag_n_sj_blacklisted += r->n_sj_blacklisted;

            RawResolveResult cr;
            cr.n_sj_blacklisted = frag_n_sj_blacklisted;
            bool resolved = ctx._resolve_core(
                frag.exons, frag.introns,
                frag.genomic_footprint(), cr, scratch);

            // SRD v2: _resolve_core now returns true even when t_inds
            // is empty (truly intergenic).  Treat empty-t_inds as the
            // legacy "unresolved" path for stats but ALSO append to
            // the buffer (for unique mappers) so calibration can
            // categorize the fragment as INTERGENIC.
            if (!resolved || cr.t_inds.empty()) {
                // Defer intergenic counting until after all hits
                // are processed to avoid multi-counting.
                if (frag.has_introns()) {
                    any_unresolved_spliced = true;
                } else {
                    any_unresolved_unspliced = true;
                }

                // SRD v2: unique-mapper zero-candidate fragments are
                // appended to the buffer so that calibration's
                // geometric categorization can label them INTERGENIC.
                //
                // We do NOT set any_hit_resolved here -- the
                // n_intergenic_unspliced/_spliced telemetry counter
                // (driven by !any_hit_resolved below) must keep firing
                // so the PipelineStats accountability sum stays
                // consistent.  The C++ scorer skips empty-t_inds
                // fragments silently (n_cand <= 0 early skip) so they
                // don't increment stat_gated either.
                if (resolved && is_unique_mapper) {
                    ResolvedFragment ig_result = ResolvedFragment::from_core(cr);
                    ig_result.num_hits = num_hits;
                    ig_result.nm = frag.nm;
                    accumulator.append(ig_result, frag_id);

                    // Phase 0: deposit the intergenic fragment's genomic-DNA
                    // mass into the calibration accumulator. Without this, an
                    // intergenic region's contained mass is identically zero
                    // and the count-clue density loses its baseline signal.
                    // (Intergenic ⇒ no candidate transcripts ⇒ no implicit
                    // introns; cr drives the channel/spans uniformly.)
                    deposit_to_accumulator(frag, cr);
                }

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
                } else if (result.splice_type == SPLICE_SPLICED_UNANNOT) {
                    stats.n_with_unannotated_sj++;
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
                    // JUNCTION. `sj_strand` is qualified to POS/NEG here and is the
                    // motif strand of every annotated junction the fragment crosses,
                    // so `sense` is exactly the align==sj bit the 2×2 is built from
                    // and `sj_key_*` is necessarily set.
                    //
                    // Two parallel int8 label vectors used to be pushed here as well,
                    // one pair per fragment, and the Python 2×2 was counted from them.
                    // The table's marginal IS that 2×2 (verified exact on 32 synthetic
                    // conditions and 4 real libraries), so they were pure duplication
                    // and were deleted 2026-07-28. `stats.n_strand_trained` below still
                    // counts the fragments, and Python asserts it against the table's
                    // total depth — the invariant that one fragment credits one junction.
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
                            strand_obs.exonic_obs.push_back(static_cast<int8_t>(result.align_strand));
                            strand_obs.exonic_truth.push_back(static_cast<int8_t>(t_strand));
                        }
                    }
                }

                // Fractional accumulator deposit (resolved unique-mapper,
                // non-chimeric). See deposit_to_accumulator above.
                deposit_to_accumulator(frag, cr);
            }

            // Region accumulation for resolved (non-chimeric) fragments.
            // Every resolved hit of a multimapper accumulates with
            // 1/num_hits weight; summed across all hits of a single
            // molecule the total is 1.0, matching the 1/NH crediting
            // convention of the E_i (mappable_effective_length)
            // denominator.  No count_stats gate: that would bias the
            // numerator high by a factor of NH relative to E_i.

            accumulator.append(result, frag_id);

            if (!is_unique_mapper) {
                n_buffered_mm++;
            }
        }

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
        if (!any_hit_resolved &&
            (any_unresolved_spliced || any_unresolved_unspliced)) {
            if (any_unresolved_spliced) {
                stats.n_intergenic_spliced++;
            } else {
                stats.n_intergenic_unspliced++;
            }
        }

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
        stats_dict["n_intergenic_unspliced"] = stats_.n_intergenic_unspliced;
        stats_dict["n_intergenic_spliced"] = stats_.n_intergenic_spliced;
        stats_dict["n_with_exon"] = stats_.n_with_exon;
        stats_dict["n_with_annotated_sj"] = stats_.n_with_annotated_sj;
        stats_dict["n_with_unannotated_sj"] = stats_.n_with_unannotated_sj;
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
        strand_dict["exonic_obs"]           = vec_to_ndarray(std::move(strand_obs_.exonic_obs));
        strand_dict["exonic_truth"]         = vec_to_ndarray(std::move(strand_obs_.exonic_truth));

        // The per-junction SJ strand table, as six parallel arrays.  Sorted by
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
        // ⚠ And this copy cannot be a memcpy: the accumulator stores AoS (a 48 B `Node` interleaving four
        // channel pairs) while the payload is SoA, so the transpose is a strided read by construction. That
        // is the price of the hot struct being 48 B with no padding, which is the right trade — this loop is
        // O(partition) once per scan, the deposit is O(fragments).
        if (acc_set_) {
            using rigel::accumulator::Accumulator;
            using rigel::accumulator::ContiguousEdge;
            using rigel::accumulator::JunctionEdge;
            using rigel::accumulator::Node;
            using rigel::accumulator::kNFragmentPools;
            using rigel::accumulator::kNStrandColumns;

            const std::size_t n_refs = acc_set_->n_refs();

            std::vector<int64_t> ref_node_offsets(n_refs + 1, 0);
            std::vector<int64_t> ref_edge_offsets(n_refs + 1, 0);
            std::vector<int64_t> ref_sj_offsets(n_refs + 1, 0);
            for (std::size_t f = 0; f < n_refs; ++f) {
                const Accumulator& a = acc_set_->at(static_cast<int32_t>(f));
                ref_node_offsets[f + 1] =
                    ref_node_offsets[f] + static_cast<int64_t>(a.n_nodes());
                ref_edge_offsets[f + 1] =
                    ref_edge_offsets[f] + static_cast<int64_t>(a.n_edges());
                ref_sj_offsets[f + 1] =
                    ref_sj_offsets[f] + static_cast<int64_t>(a.n_junctions());
            }
            const auto n_nodes = static_cast<std::size_t>(ref_node_offsets.back());
            const auto n_edges = static_cast<std::size_t>(ref_edge_offsets.back());
            const auto n_sj    = static_cast<std::size_t>(ref_sj_offsets.back());

            std::vector<uint32_t> node_contained_count(n_nodes * kNStrandColumns, 0u);
            std::vector<uint64_t> node_contained_inv_length_sum(n_nodes * kNStrandColumns, 0u);
            std::vector<uint64_t> node_contained_length_sum(n_nodes * kNStrandColumns, 0u);
            std::vector<uint32_t> node_spanning_count(n_nodes * kNStrandColumns, 0u);
            std::vector<uint64_t> node_spanning_inv_length_sum(n_nodes * kNStrandColumns, 0u);
            std::vector<uint64_t> node_spanning_length_sum(n_nodes * kNStrandColumns, 0u);
            std::vector<uint32_t> node_start_count(n_nodes, 0u);
            std::vector<uint32_t> edge_unspliced_count(n_edges * kNStrandColumns, 0u);
            std::vector<uint64_t> edge_unspliced_inv_length_sum(n_edges * kNStrandColumns, 0u);
            std::vector<uint64_t> edge_unspliced_length_sum(n_edges * kNStrandColumns, 0u);
            std::vector<uint32_t> edge_spliced_count(n_edges * kNStrandColumns, 0u);
            std::vector<uint64_t> edge_spliced_inv_length_sum(n_edges * kNStrandColumns, 0u);
            std::vector<uint64_t> edge_spliced_length_sum(n_edges * kNStrandColumns, 0u);
            std::vector<uint32_t> sj_count(n_sj * kNStrandColumns, 0u);
            std::vector<uint64_t> sj_inv_length_sum(n_sj * kNStrandColumns, 0u);
            std::vector<uint64_t> sj_length_sum(n_sj * kNStrandColumns, 0u);

            const std::size_t pool_row = static_cast<std::size_t>(max_length_) + 1;
            std::vector<int64_t> pool_lengths(kNFragmentPools * pool_row, 0);
            // ⭐ C1: the unconditional histogram, summed over references exactly like the pools.
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
                const Node* nodes = a.nodes_data();
                const ContiguousEdge* edges = a.edges_data();
                const JunctionEdge* junctions = a.junctions_data();
                const uint32_t* starts = a.node_start_count_data();
                const auto node_base = static_cast<std::size_t>(ref_node_offsets[f]);
                const auto edge_base = static_cast<std::size_t>(ref_edge_offsets[f]);
                const auto sj_base   = static_cast<std::size_t>(ref_sj_offsets[f]);

                for (std::size_t i = 0; i < a.n_nodes(); ++i) {
                    node_start_count[node_base + i] = starts[i];
                    for (std::size_t c = 0; c < kNStrandColumns; ++c) {
                        const std::size_t o = (node_base + i) * kNStrandColumns + c;
                        node_contained_count[o]   = nodes[i].contained_count[c];
                        node_contained_inv_length_sum[o] = nodes[i].contained_inv_length_sum[c];
                        node_contained_length_sum[o] = nodes[i].contained_length_sum[c];
                        node_spanning_count[o]    = nodes[i].spanning_count[c];
                        node_spanning_inv_length_sum[o]  = nodes[i].spanning_inv_length_sum[c];
                        node_spanning_length_sum[o] = nodes[i].spanning_length_sum[c];
                    }
                }
                for (std::size_t i = 0; i < a.n_edges(); ++i) {
                    for (std::size_t c = 0; c < kNStrandColumns; ++c) {
                        const std::size_t o = (edge_base + i) * kNStrandColumns + c;
                        edge_unspliced_count[o]   = edges[i].unspliced_count[c];
                        edge_unspliced_inv_length_sum[o] = edges[i].unspliced_inv_length_sum[c];
                        edge_unspliced_length_sum[o] = edges[i].unspliced_length_sum[c];
                        edge_spliced_count[o]     = edges[i].spliced_count[c];
                        edge_spliced_inv_length_sum[o]   = edges[i].spliced_inv_length_sum[c];
                        edge_spliced_length_sum[o] = edges[i].spliced_length_sum[c];
                    }
                }
                for (std::size_t i = 0; i < a.n_junctions(); ++i) {
                    for (std::size_t c = 0; c < kNStrandColumns; ++c) {
                        const std::size_t o = (sj_base + i) * kNStrandColumns + c;
                        sj_count[o]   = junctions[i].count[c];
                        sj_inv_length_sum[o] = junctions[i].inv_length_sum[c];
                        sj_length_sum[o] = junctions[i].length_sum[c];
                    }
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
                // ⭐ C1. Same size guard, same reason: a silently skipped reference would leave the
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

            nb::dict cal;
            // Echo the partition back, so a consumer can locate every object without reloading the index.
            cal["cut_positions"]   = vec_to_ndarray(std::vector<int64_t>(cut_positions_));
            cal["ref_cut_offsets"] = vec_to_ndarray(std::vector<int64_t>(ref_cut_offsets_));
            cal["ref_node_offsets"] = vec_to_ndarray(std::move(ref_node_offsets));
            cal["ref_edge_offsets"] = vec_to_ndarray(std::move(ref_edge_offsets));
            cal["ref_sj_offsets"]   = vec_to_ndarray(std::move(ref_sj_offsets));

            cal["node_contained_count"]   = vec_to_ndarray(std::move(node_contained_count));
            cal["node_contained_length_sum"] = vec_to_ndarray(std::move(node_contained_length_sum));
            cal["node_contained_inv_length_sum"] = vec_to_ndarray(std::move(node_contained_inv_length_sum));
            cal["node_spanning_count"]    = vec_to_ndarray(std::move(node_spanning_count));
            cal["node_spanning_length_sum"] = vec_to_ndarray(std::move(node_spanning_length_sum));
            cal["node_spanning_inv_length_sum"]  = vec_to_ndarray(std::move(node_spanning_inv_length_sum));
            cal["node_start_count"]       = vec_to_ndarray(std::move(node_start_count));
            cal["edge_unspliced_count"]   = vec_to_ndarray(std::move(edge_unspliced_count));
            cal["edge_unspliced_length_sum"] = vec_to_ndarray(std::move(edge_unspliced_length_sum));
            cal["edge_unspliced_inv_length_sum"] = vec_to_ndarray(std::move(edge_unspliced_inv_length_sum));
            cal["edge_spliced_count"]     = vec_to_ndarray(std::move(edge_spliced_count));
            cal["edge_spliced_length_sum"] = vec_to_ndarray(std::move(edge_spliced_length_sum));
            cal["edge_spliced_inv_length_sum"]   = vec_to_ndarray(std::move(edge_spliced_inv_length_sum));
            cal["sj_count"]               = vec_to_ndarray(std::move(sj_count));
            cal["sj_length_sum"] = vec_to_ndarray(std::move(sj_length_sum));
            cal["sj_inv_length_sum"]             = vec_to_ndarray(std::move(sj_inv_length_sum));
            cal["pool_lengths"]           = vec_to_ndarray(std::move(pool_lengths));
            cal["deposited_lengths"]      = vec_to_ndarray(std::move(deposited_lengths));

            // The QC denominators (design §10.3). Every conservation statement downstream has to be able
            // to name what it excluded, and none of these is derivable after the fact.
            nb::dict qc_dict;
            qc_dict["deposited"]                = qc.deposited;
            qc_dict["dropped_too_long"]         = qc.dropped_too_long;
            qc_dict["dropped_empty"]            = qc.dropped_empty;
            qc_dict["dropped_strand_undefined"] = qc.dropped_strand_undefined;
            qc_dict["deferred_undetermined_gap"] = qc.deferred_undetermined_gap;
            qc_dict["unannotated_introns"]      = qc.unannotated_introns;
            qc_dict["contradictory_sj_strand"]  = qc.contradictory_sj_strand;
            qc_dict["introns_absorbed"]         = qc.introns_absorbed;
            cal["qc"] = qc_dict;

            // ⭐ The side buffer and the umbrella census, through the SAME two exporters the parity
            // surface uses, so the payload and the bound Accumulator cannot disagree about a key.
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
// BamAnnotationWriter — second-pass BAM tag stamping
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

class BamAnnotationWriter {
public:
    FragmentResolver* ctx_;
    SJTagMode sj_tag_mode_ = SJTagMode::XS_ONLY;
    bool skip_duplicates_ = true;
    bool include_multimap_ = false;
    std::vector<int32_t> tid_to_ref_id_;

    BamAnnotationWriter(FragmentResolver& ctx,
                        const std::string& sj_tag_spec,
                        bool skip_duplicates,
                        bool include_multimap)
        : ctx_(&ctx),
          skip_duplicates_(skip_duplicates),
          include_multimap_(include_multimap)
    {
        sj_tag_mode_ = parse_sj_tag_spec(sj_tag_spec);
    }

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
        int64_t n_chimeric = 0;
        int64_t n_records_written = 0;
        int64_t n_filtered_passthrough = 0;

        // Scratch buffer for resolve calls (single-threaded, reused)
        ResolverScratch scratch(ctx_->n_transcripts_);

        auto process_and_write = [&]() {
            if (light_group.empty()) return;

            n_read_groups++;

            // Determine NH
            int32_t nh = 1;
            for (const auto& r : light_group) {
                if (r.nh > 1) { nh = r.nh; break; }
            }

            // Detect multimap
            bool has_secondary = false;
            for (const auto& r : light_group) {
                if (r.is_secondary()) { has_secondary = true; break; }
            }
            bool is_multimap = (nh > 1) || has_secondary;

            if (is_multimap && !include_multimap_) {
                // Write records without annotation (match pass 1 skip)
                // Pass 1 assigns frag_id per qname group before worker-side
                // multimapper filtering, so we must still advance frag_id
                // here to keep pass-2 lookup synchronized.
                for (size_t i = 0; i < raw_group.size(); i++) {
                    bam_aux_update_int(
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

            // Group into hits + secondary locations
            auto hit_result = group_records_by_hit(light_group);

            // Multimapper secondary pairing
            std::vector<std::pair<std::vector<ParsedAlignment*>,
                                  std::vector<ParsedAlignment*>>> secondary_pairs;
            if (!hit_result.sec_r1_locs.empty() ||
                !hit_result.sec_r2_locs.empty()) {
                secondary_pairs = pair_multimapper_reads(
                    hit_result.sec_r1_locs,
                    hit_result.sec_r2_locs,
                    *ctx_, scratch);
            }

            // Combine all hits
            auto& all_hits = hit_result.hits;
            all_hits.insert(all_hits.end(),
                            secondary_pairs.begin(),
                            secondary_pairs.end());

            int32_t num_hits = std::max(nh,
                static_cast<int32_t>(all_hits.size()));

            // Look up annotation for this frag_id
            auto ann_it = frag_to_row.find(frag_id);
            bool has_ann = (ann_it != frag_to_row.end());

            // Deduplicate raw records across hits: the cross-pair
            // combinatorial block in pair_multimapper_reads intentionally
            // produces overlapping hits for EM scoring, so a single raw
            // record can appear in multiple hits.  For the annotated BAM
            // output we must write each input record EXACTLY once.
            // First-containing-hit wins.
            std::vector<bool> rec_written(raw_group.size(), false);

            for (int hit_idx = 0;
                 hit_idx < static_cast<int>(all_hits.size());
                 hit_idx++)
            {
                const auto& [r1_reads, r2_reads] = all_hits[hit_idx];

                // Gather bam1_t* for all records in this hit,
                // skipping any already emitted by an earlier hit.
                // Stamp the per-record ZB:i tag (splice-artifact count)
                // on each bam1_t as we select it \u2014 this is a record-local
                // property of the CIGAR, not a hit-level value, so it
                // lives outside the uniform stamp_and_write_hit helper.
                std::vector<bam1_t*> hit_raws;
                for (const auto* lr : r1_reads) {
                    ptrdiff_t idx = lr - &light_group[0];
                    if (rec_written[idx]) continue;
                    rec_written[idx] = true;
                    bam_aux_update_int(
                        raw_group[idx], "ZB",
                        static_cast<int>(light_group[idx].n_sj_blacklisted));
                    hit_raws.push_back(raw_group[idx]);
                }
                for (const auto* lr : r2_reads) {
                    ptrdiff_t idx = lr - &light_group[0];
                    if (rec_written[idx]) continue;
                    rec_written[idx] = true;
                    bam_aux_update_int(
                        raw_group[idx], "ZB",
                        static_cast<int>(light_group[idx].n_sj_blacklisted));
                    hit_raws.push_back(raw_group[idx]);
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
                    if (num_hits == 1) {
                        is_primary = true;
                    } else if (best_tid_val >= 0) {
                        // Resolve this hit to see if assigned transcript
                        // is among its candidates
                        AssembledFragment frag = build_fragment(
                            r1_reads, r2_reads);
                        if (!frag.exons.empty()) {
                            RawResolveResult cr;
                            if (ctx_->_resolve_core(
                                    frag.exons, frag.introns,
                                    frag.genomic_footprint(), cr, scratch)) {
                                for (int32_t ti : cr.t_inds) {
                                    if (ti == best_tid_val) {
                                        is_primary = true;
                                        break;
                                    }
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
                    // No EM annotation for this frag_id: deterministic
                    // intergenic gDNA fast-path.  Stamp ZF with the
                    // intergenic+gDNA outcome bits so a single AND-mask on
                    // ZF captures total gDNA.  ZW=1.0 because there is no
                    // competing hypothesis.
                    stamp_and_write_hit(
                        hit_raws, out, hdr,
                        ".", ".", ".",
                        -1, -1,
                        AF_GDNA_INTERGENIC, 1.0f,
                        ".",
                        (hit_idx == 0) ? 1 : 0,
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
                    bam_aux_update_int(
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
        while (sam_read1(fp, hdr, b) >= 0) {
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

            // Build ParsedAlignment (same as BamScanner)
            ParsedAlignment rec;
            rec.ref_id = b->core.tid;
            rec.ref_start = b->core.pos;
            rec.mate_ref_id = b->core.mtid;
            rec.mate_ref_start = b->core.mpos;
            rec.flag = flag;

            rec.nm = 0;
            uint8_t* nm_aux = bam_aux_get(b, "NM");
            if (nm_aux) rec.nm = bam_aux2i(nm_aux);

            rec.nh = 1;
            uint8_t* nh_aux = bam_aux_get(b, "NH");
            if (nh_aux) rec.nh = bam_aux2i(nh_aux);

            rec.hi = -1;
            uint8_t* hi_aux = bam_aux_get(b, "HI");
            if (hi_aux) rec.hi = bam_aux2i(hi_aux);

            int32_t mapped_ref_id = (rec.ref_id >= 0 &&
                rec.ref_id < static_cast<int32_t>(tid_to_ref_id_.size()))
                ? tid_to_ref_id_[rec.ref_id] : -1;
            parse_cigar(b, mapped_ref_id, STRAND_NONE, rec.exons, rec.sjs);
            rec.sj_strand = STRAND_NONE;
            if (!rec.sjs.empty()) {
                rec.sj_strand = read_sj_strand(b, sj_tag_mode_);
                for (auto& sj : rec.sjs) sj.strand = rec.sj_strand;
            }
            rec.ref_id = mapped_ref_id;

            // Apply the same splice-junction blacklist filter that the
            // scanner uses so the annotated BAM sees a consistent view
            // of the evidence.  Capture the drop count for the per-record
            // ZB:i tag stamped below.
            if (ctx_ && ctx_->has_sj_blacklist()) {
                int32_t n_dropped =
                    filter_blacklisted_sjs(rec.sjs, *ctx_, mapped_ref_id);
                rec.n_sj_blacklisted = static_cast<uint16_t>(n_dropped);
            }

            light_group.push_back(std::move(rec));
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
        summary["n_chimeric"] = n_chimeric;
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
            bam_aux_update_str(r, "ZT",
                static_cast<int>(strlen(zt) + 1), zt);
            bam_aux_update_str(r, "ZG",
                static_cast<int>(strlen(zg) + 1), zg);
            bam_aux_update_str(r, "ZR",
                static_cast<int>(strlen(zr) + 1), zr);
            bam_aux_update_int(r, "ZI", zi);
            bam_aux_update_int(r, "ZJ", zj);
            bam_aux_update_int(r, "ZF", zf);
            bam_aux_update_float(r, "ZW", zw);
            bam_aux_update_str(r, "ZC",
                static_cast<int>(strlen(zc) + 1), zc);
            bam_aux_update_int(r, "ZH", zh);
            bam_aux_update_int(r, "ZN", zn);
            bam_aux_update_str(r, "ZS",
                static_cast<int>(strlen(zs) + 1), zs);
            bam_aux_update_int(r, "ZL", zl);

            if (sam_write1(out, hdr, r) < 0) {
                throw std::runtime_error("Failed to write BAM record");
            }
        }
    }
};

// ================================================================
// SJ strand tag auto-detection (ports detect_sj_strand_tag)
// ================================================================

static std::string detect_sj_strand_tag_native(
    const std::string& bam_path,
    int max_spliced_reads = 1000)
{
    htsFile* fp = hts_open(bam_path.c_str(), "rb");
    if (!fp) throw std::runtime_error("Failed to open BAM: " + bam_path);

    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        throw std::runtime_error("Failed to read header: " + bam_path);
    }

    bam1_t* b = bam_init1();
    bool found_xs = false;
    bool found_ts = false;
    int n_spliced = 0;

    while (sam_read1(fp, hdr, b) >= 0 && n_spliced < max_spliced_reads) {
        if (b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL)) continue;
        if (b->core.n_cigar == 0) continue;

        // Check for CIGAR N (ref skip)
        const uint32_t* cigar = bam_get_cigar(b);
        bool has_splice = false;
        for (int i = 0; i < static_cast<int>(b->core.n_cigar); i++) {
            if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
                has_splice = true;
                break;
            }
        }
        if (!has_splice) continue;

        n_spliced++;
        if (!found_xs && bam_aux_get(b, "XS")) found_xs = true;
        if (!found_ts && bam_aux_get(b, "ts")) found_ts = true;

        if (found_xs && found_ts) break;
    }

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(fp);

    // Return detected tags in priority order (XS first)
    if (found_xs && found_ts) return "XS,ts";
    if (found_xs) return "XS";
    if (found_ts) return "ts";
    return "none";
}

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
    // ⚠ Every two-column view needs EXPLICIT STRIDES. The storage is AoS with mixed widths — a 48 B `Node`
    // holds two uint32 pairs then two uint64 pairs — so the row stride is `sizeof(struct)/sizeof(element)`,
    // not the column count. The binding this replaced passed no strides at all and worked only because the
    // old `Region` was exactly four contiguous uint32.
    {
        using rigel::accumulator::Accumulator;
        using rigel::accumulator::ContiguousEdge;
        using rigel::accumulator::DepositScratch;
        using rigel::accumulator::GapHypothesis;
        using rigel::accumulator::OfferedFragment;
        using rigel::accumulator::JunctionEdge;
        using rigel::accumulator::Node;
        using rigel::accumulator::kNFragmentPools;
        using rigel::accumulator::kNStrandColumns;

        // One scratch per bound instance. The class is single-threaded from Python and the scan path uses
        // its own per-worker scratch, so this is only here to keep the signature allocation-free.
        static thread_local DepositScratch binding_scratch;

        nb::class_<Accumulator>(m, "Accumulator")
            .def("__init__",
                 [](Accumulator* self,
                    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> cuts,
                    nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig> node_types,
                    int max_length,
                    int32_t ref) {
                     std::vector<int64_t> c(cuts.data(), cuts.data() + cuts.shape(0));
                     std::vector<uint8_t> t(node_types.data(),
                                            node_types.data() + node_types.shape(0));
                     new (self) Accumulator(std::move(c), std::move(t), max_length, ref);
                 },
                 nb::arg("cuts"),
                 nb::arg("node_types"),
                 nb::arg("max_length"),
                 nb::arg("ref"),
                 "One reference's sorted cut positions, one coarse type per node, the\n"
                 "fragment-length limit (which is also the pool-histogram width), and WHICH\n"
                 "reference this is — stamped into every deferred record, because the second\n"
                 "pass replays those onto that reference's cut axis.")
            .def("set_junctions",
                 [](Accumulator& a,
                    nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig> offsets,
                    nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig> acceptor_cut,
                    nb::ndarray<const int8_t, nb::ndim<1>, nb::c_contig> sj_strand) {
                     a.set_junctions(
                         std::vector<int32_t>(offsets.data(), offsets.data() + offsets.shape(0)),
                         std::vector<int32_t>(acceptor_cut.data(),
                                              acceptor_cut.data() + acceptor_cut.shape(0)),
                         std::vector<int8_t>(sj_strand.data(),
                                             sj_strand.data() + sj_strand.shape(0)));
                 },
                 nb::arg("offsets"),
                 nb::arg("acceptor_cut"),
                 nb::arg("sj_strand"),
                 "The junction CSR for THIS reference, keyed by the ref-local donor cut\n"
                 "index. The junction-edge id is the slot; slot order is a contract.")
            .def_prop_ro("n_nodes",     [](const Accumulator& a) { return a.n_nodes(); })
            .def_prop_ro("n_edges",     [](const Accumulator& a) { return a.n_edges(); })
            .def_prop_ro("n_junctions", [](const Accumulator& a) { return a.n_junctions(); })
            .def_prop_ro("n_cuts",      [](const Accumulator& a) { return a.n_cuts(); })
            .def_prop_ro("max_length",  [](const Accumulator& a) { return a.max_length(); })

            // ── nodes ────────────────────────────────────────────────────────────────────────────────
            .def_prop_ro("node_contained_count", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Node) / sizeof(uint32_t);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<2>>(
                    &a.nodes_data()[0].contained_count[0], {a.n_nodes(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("node_spanning_count", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Node) / sizeof(uint32_t);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<2>>(
                    &a.nodes_data()[0].spanning_count[0], {a.n_nodes(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("node_contained_inv_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Node) / sizeof(uint64_t);
                return nb::ndarray<nb::numpy, const uint64_t, nb::ndim<2>>(
                    &a.nodes_data()[0].contained_inv_length_sum[0], {a.n_nodes(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("node_contained_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Node) / sizeof(uint64_t);
                return nb::ndarray<nb::numpy, const uint64_t, nb::ndim<2>>(
                    &a.nodes_data()[0].contained_length_sum[0], {a.n_nodes(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("node_spanning_inv_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Node) / sizeof(uint64_t);
                return nb::ndarray<nb::numpy, const uint64_t, nb::ndim<2>>(
                    &a.nodes_data()[0].spanning_inv_length_sum[0], {a.n_nodes(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("node_spanning_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(Node) / sizeof(uint64_t);
                return nb::ndarray<nb::numpy, const uint64_t, nb::ndim<2>>(
                    &a.nodes_data()[0].spanning_length_sum[0], {a.n_nodes(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("node_start_count", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<1>>(
                    a.node_start_count_data(), {a.n_nodes()}, h).cast();
            })

            // ── contiguous edges ─────────────────────────────────────────────────────────────────────
            .def_prop_ro("edge_unspliced_count", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(ContiguousEdge) / sizeof(uint32_t);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<2>>(
                    &a.edges_data()[0].unspliced_count[0], {a.n_edges(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("edge_spliced_count", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(ContiguousEdge) / sizeof(uint32_t);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<2>>(
                    &a.edges_data()[0].spliced_count[0], {a.n_edges(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("edge_unspliced_inv_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(ContiguousEdge) / sizeof(uint64_t);
                return nb::ndarray<nb::numpy, const uint64_t, nb::ndim<2>>(
                    &a.edges_data()[0].unspliced_inv_length_sum[0], {a.n_edges(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("edge_unspliced_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(ContiguousEdge) / sizeof(uint64_t);
                return nb::ndarray<nb::numpy, const uint64_t, nb::ndim<2>>(
                    &a.edges_data()[0].unspliced_length_sum[0], {a.n_edges(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("edge_spliced_inv_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(ContiguousEdge) / sizeof(uint64_t);
                return nb::ndarray<nb::numpy, const uint64_t, nb::ndim<2>>(
                    &a.edges_data()[0].spliced_inv_length_sum[0], {a.n_edges(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("edge_spliced_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(ContiguousEdge) / sizeof(uint64_t);
                return nb::ndarray<nb::numpy, const uint64_t, nb::ndim<2>>(
                    &a.edges_data()[0].spliced_length_sum[0], {a.n_edges(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })

            // ── junction edges ───────────────────────────────────────────────────────────────────────
            .def_prop_ro("sj_count", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(JunctionEdge) / sizeof(uint32_t);
                return nb::ndarray<nb::numpy, const uint32_t, nb::ndim<2>>(
                    &a.junctions_data()[0].count[0], {a.n_junctions(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("sj_inv_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(JunctionEdge) / sizeof(uint64_t);
                return nb::ndarray<nb::numpy, const uint64_t, nb::ndim<2>>(
                    &a.junctions_data()[0].inv_length_sum[0], {a.n_junctions(), kNStrandColumns}, h,
                    {row, int64_t{1}}).cast();
            })
            .def_prop_ro("sj_length_sum", [](nb::handle h) {
                auto& a = nb::cast<Accumulator&>(h);
                constexpr int64_t row = sizeof(JunctionEdge) / sizeof(uint64_t);
                return nb::ndarray<nb::numpy, const uint64_t, nb::ndim<2>>(
                    &a.junctions_data()[0].length_sum[0], {a.n_junctions(), kNStrandColumns}, h,
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
                const auto& c = a.counters();
                nb::dict qc;
                qc["deposited"]                = c.deposited;
                qc["dropped_too_long"]         = c.dropped_too_long;
                qc["dropped_empty"]            = c.dropped_empty;
                qc["dropped_strand_undefined"] = c.dropped_strand_undefined;
                qc["deferred_undetermined_gap"] = c.deferred_undetermined_gap;
                qc["unannotated_introns"]      = c.unannotated_introns;
                qc["contradictory_sj_strand"]  = c.contradictory_sj_strand;
                qc["introns_absorbed"]         = c.introns_absorbed;
                return qc;
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

            .def("node_of_pos",
                 [](const Accumulator& a, int64_t pos) { return a.node_of_pos(pos); },
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
            // ⛔ Exposed rather than reimplemented in Python: C0/C2 left the tool with ONE definition of
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

            // Element-wise sum of `other` into this accumulator — the per-worker merge the parallel scan
            // performs internally, exposed so the DETERMINISM contract can be tested directly: shard one
            // fragment corpus K ways, merge, and require bit-identity with the unsharded run. Every channel
            // is an integer now, so that identity is exact rather than approximate.
            .def("merge_from",
                 [](Accumulator& a, const Accumulator& other) { a.merge_from(other); },
                 nb::arg("other"));
    }


    nb::class_<BamScanner>(m, "BamScanner")
        .def(nb::init<FragmentResolver&, const std::string&, bool, bool>(),
             nb::arg("ctx"),
             nb::arg("sj_tag_spec"),
             nb::arg("skip_duplicates") = true,
             nb::arg("include_multimap") = false,
             "Create a BAM scanner.\n\n"
             "Parameters\n"
             "----------\n"
             "ctx : FragmentResolver\n"
             "    The resolve context from index building.\n"
             "sj_tag_spec : str\n"
             "    SJ strand tag specification: 'XS', 'ts', 'XS,ts', 'ts,XS', or 'none'.\n"
             "skip_duplicates : bool\n"
             "    Skip duplicate-flagged reads (default True).\n"
             "include_multimap : bool\n"
             "    Include multimapping reads (default False).\n",
             nb::keep_alive<1, 2>()  // BamScanner must keep FragmentResolver alive
        )
        .def("scan", &BamScanner::scan,
             nb::arg("bam_path"),
             nb::arg("chunk_callback"),
             nb::arg("t_strand_arr"),
             nb::arg("chunk_size") = 1000000,
             nb::arg("n_workers") = 1,
             nb::arg("n_decomp_threads") = 2,
             nb::arg("qname_batch_size") = 512,
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
             "t_strand_arr : list[int]\n"
             "    Per-transcript strand array (int32).\n"
             "chunk_size : int\n"
             "    Target number of fragments per chunk (default 1000000).\n"
             "n_workers : int\n"
             "    Number of worker threads (default 1).\n"
             "n_decomp_threads : int\n"
             "    Number of htslib BGZF decompression threads (default 2).\n"
             "qname_batch_size : int\n"
             "    Read-name groups per input queue item (default 512).\n\n"
             "Returns\n"
             "-------\n"
             "dict\n"
             "    Dict with keys: 'stats', 'strand_observations',\n"
             "    'calibration'.\n")
          .def("set_regions",
                 [](BamScanner& self,
                     nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> cut_positions,
                     nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> ref_cut_offsets,
                     int32_t n_refs,
                     nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig> node_types,
                     int max_length) {
                      self.set_regions(cut_positions, ref_cut_offsets, n_refs,
                                       node_types, max_length);
                 },
                 nb::arg("cut_positions"),
                 nb::arg("ref_cut_offsets"),
                 nb::arg("n_refs"),
                 nb::arg("node_types"),
                 nb::arg("max_length"),
                 "Install the accumulator's node partition. Call set_junctions next.\n\n"
                 "cut_positions : int64[n_cuts_total]\n"
                 "    Flat sorted cut positions for all references. A reference\n"
                 "    contributing c cuts owns c-1 nodes and c-2 interior lines.\n"
                 "ref_cut_offsets : int64[n_refs + 1]\n"
                 "    Per-ref offsets into cut_positions; ref f spans\n"
                 "    cut_positions[ref_cut_offsets[f]:ref_cut_offsets[f+1]].\n"
                 "n_refs : int\n"
                 "    Number of references (matches FragmentResolver).\n"
                 "node_types : uint8[n_nodes_total]\n"
                 "    Coarse node type (0=intergenic, 1=intron, 2=exon), ref-major.\n"
                 "    It types the fragment-length pools.\n"
                 "max_length : int\n"
                 "    The fragment-length limit applied to L, and the width of the\n"
                 "    pool histograms. Must be >= 1.\n")
          .def("set_junctions",
                 [](BamScanner& self,
                     nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> offsets,
                     nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig> acceptor_cut,
                     nb::ndarray<const int8_t, nb::ndim<1>, nb::c_contig> sj_strand) {
                      self.set_junctions(offsets, acceptor_cut, sj_strand);
                 },
                 nb::arg("offsets"),
                 nb::arg("acceptor_cut"),
                 nb::arg("sj_strand"),
                 "Install the annotated junction edges, as build_junction_edge_arrays\n"
                 "emits them. A SECOND call because set_regions refuses to run twice.\n\n"
                 "offsets : int64[n_cuts_total + 1]\n"
                 "    CSR over the flat cut axis, keyed by the DONOR cut index.\n"
                 "acceptor_cut : int64[n_junctions]\n"
                 "    Flat cut index of each junction's high end. The junction-edge id\n"
                 "    IS the slot here; edges_df.edge_row is a join key and is not\n"
                 "    passed. Slot order is a contract: sort on\n"
                 "    (donor cut, acceptor cut, sj_strand).\n"
                 "sj_strand : int8[n_junctions]\n"
                 "    Each junction's ANNOTATED strand.\n\n"
                 "Not optional: scan() refuses to run without it, because a missing\n"
                 "table makes every observed intron read as unannotated and empties\n"
                 "the junction edges silently. Pass empty arrays to declare none.\n")
        ;

    nb::class_<BamAnnotationWriter>(m, "BamAnnotationWriter")
        .def(nb::init<FragmentResolver&, const std::string&, bool, bool>(),
             nb::arg("ctx"),
             nb::arg("sj_tag_spec"),
             nb::arg("skip_duplicates") = true,
             nb::arg("include_multimap") = false,
             "Create a BAM annotation writer.\n\n"
             "Parameters\n"
             "----------\n"
             "ctx : FragmentResolver\n"
             "    The resolve context from index building.\n"
             "sj_tag_spec : str\n"
             "    SJ strand tag specification (must match pass 1).\n"
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

    m.def("detect_sj_strand_tag", &detect_sj_strand_tag_native,
          nb::arg("bam_path"),
          nb::arg("max_spliced_reads") = 1000,
          "Auto-detect SJ strand tags in a BAM file.\n\n"
          "Returns a tag specification string: 'XS', 'ts', 'XS,ts', or 'none'.");
}
