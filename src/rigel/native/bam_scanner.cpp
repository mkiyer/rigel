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
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <sys/stat.h>

#include "resolve_context.h"
#include "thread_queue.h"

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
        case FRAG_CHIMERIC:          return "chimeric";
        case -1:                     return "intergenic";
        default:                     return "unknown";
    }
}

static const char* splice_type_label(int code) {
    switch (code) {
        case SPLICE_UNSPLICED:       return "unspliced";
        case SPLICE_SPLICED_UNANNOT: return "spliced_unannot";
        case SPLICE_SPLICED_ANNOT:   return "spliced_annot";
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

// ================================================================
// Model training observation collectors
// ================================================================

struct StrandObservations {
    // exonic_spliced: (exon_strand, sj_strand) for strand-qualified
    std::vector<int8_t> exonic_spliced_obs;
    std::vector<int8_t> exonic_spliced_truth;

    // exonic: (exon_strand, gene_strand) for unique-gene unambiguous strand
    std::vector<int8_t> exonic_obs;
    std::vector<int8_t> exonic_truth;

    // intergenic: (intergenic_strand, POS) for unique-mapper intergenic
    std::vector<int8_t> intergenic_obs;
    std::vector<int8_t> intergenic_truth;
};

struct FragLenObservations {
    // Resolved fragments with unambiguous fragment length
    std::vector<int32_t> lengths;
    std::vector<int8_t>  splice_types;  // SpliceType enum: {0, 1, 2}

    // Intergenic fragment lengths (splice_type = -1 meaning None)
    std::vector<int32_t> intergenic_lengths;
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

    // Fragment length model training
    int64_t n_frag_length_unambiguous = 0;
    int64_t n_frag_length_ambiguous = 0;
    int64_t n_frag_length_intergenic = 0;

    // Multimapper
    int64_t n_multimapper_groups = 0;
    int64_t n_multimapper_alignments = 0;

    // Splice-junction artifact blacklist (per read, across all alignments)
    int64_t n_sj_observed = 0;
    int64_t n_sj_blacklisted = 0;

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
        n_frag_length_unambiguous += o.n_frag_length_unambiguous;
        n_frag_length_ambiguous += o.n_frag_length_ambiguous;
        n_frag_length_intergenic += o.n_frag_length_intergenic;
        n_multimapper_groups += o.n_multimapper_groups;
        n_multimapper_alignments += o.n_multimapper_alignments;
        n_sj_observed += o.n_sj_observed;
        n_sj_blacklisted += o.n_sj_blacklisted;
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
    FragLenObservations fraglen_obs;
    RegionAccumulator region_acc;

    explicit WorkerState(int32_t n_transcripts)
        : scratch(n_transcripts) {}
};

// Merge observations
static void merge_strand_obs(StrandObservations& dst, StrandObservations& src) {
    dst.exonic_spliced_obs.insert(dst.exonic_spliced_obs.end(),
                                  src.exonic_spliced_obs.begin(),
                                  src.exonic_spliced_obs.end());
    dst.exonic_spliced_truth.insert(dst.exonic_spliced_truth.end(),
                                    src.exonic_spliced_truth.begin(),
                                    src.exonic_spliced_truth.end());
    dst.exonic_obs.insert(dst.exonic_obs.end(),
                          src.exonic_obs.begin(), src.exonic_obs.end());
    dst.exonic_truth.insert(dst.exonic_truth.end(),
                            src.exonic_truth.begin(), src.exonic_truth.end());
    dst.intergenic_obs.insert(dst.intergenic_obs.end(),
                              src.intergenic_obs.begin(),
                              src.intergenic_obs.end());
    dst.intergenic_truth.insert(dst.intergenic_truth.end(),
                                src.intergenic_truth.begin(),
                                src.intergenic_truth.end());
}

static void merge_fraglen_obs(FragLenObservations& dst, FragLenObservations& src) {
    dst.lengths.insert(dst.lengths.end(),
                       src.lengths.begin(), src.lengths.end());
    dst.splice_types.insert(dst.splice_types.end(),
                            src.splice_types.begin(), src.splice_types.end());
    dst.intergenic_lengths.insert(dst.intergenic_lengths.end(),
                                  src.intergenic_lengths.begin(),
                                  src.intergenic_lengths.end());
}

static void merge_region_acc(RegionAccumulator& dst, RegionAccumulator& src) {
    dst.merge_from(src);
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

    rec.sj_strand = read_sj_strand(b, sj_tag_mode);

    // CIGAR → exon blocks + splice junctions
    int32_t mapped_ref_id = (rec.ref_id >= 0 &&
        rec.ref_id < static_cast<int32_t>(tid_to_ref_id.size()))
        ? tid_to_ref_id[rec.ref_id] : -1;
    parse_cigar(b, mapped_ref_id, rec.sj_strand,
                rec.exons, rec.sjs);
    rec.ref_id = mapped_ref_id;

    // Splice-junction artifact filtering.  Per-read, immediately after
    // CIGAR parsing.  Downstream fragment assembly unions across mates,
    // so this is the permissive rule: a junction survives as long as at
    // least one read supports it with anchors above the blacklist envelope.
    if (resolver && resolver->has_sj_blacklist() && !rec.sjs.empty()) {
        if (n_sj_observed) *n_sj_observed += static_cast<int64_t>(rec.sjs.size());
        int32_t n_dropped = filter_blacklisted_sjs(rec.sjs, *resolver, mapped_ref_id);
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
    FragLenObservations fraglen_obs_;
    RegionAccumulator region_acc_;

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
                  int n_decomp_threads = 2)
    {
        if (n_workers < 1) n_workers = 1;
        if (chunk_size < 1) chunk_size = 1;

        // Two queues: input (SPMC) and output (MPSC)
        BoundedQueue<QnameGroup> input_queue(
            static_cast<size_t>(n_workers * 2));
        BoundedQueue<FragmentAccumulator> output_queue(
            static_cast<size_t>(n_workers * 2));

        // Per-worker state (local to scan — not a class member)
        std::vector<std::unique_ptr<WorkerState>> worker_states;
        worker_states.reserve(n_workers);
        for (int i = 0; i < n_workers; i++) {
            int32_t n_transcripts = ctx_->n_transcripts_;
            auto ws = std::make_unique<WorkerState>(n_transcripts);
            if (ctx_->has_region_index()) {
                ws->region_acc.init(
                    ctx_->n_regions(),
                    ctx_->region_cr_,
                    &ctx_->id_to_ref_);
            }
            // Pre-allocate accumulator for chunk_size
            ws->accumulator.reserve(chunk_size, chunk_size * 3 / 2);
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
                QnameGroup group;
                while (input_queue.pop(group)) {
                    process_qname_group_threaded(
                        group, *ctx, ws, include_multimap);
                    // Emit a chunk when accumulator reaches threshold
                    if (ws.accumulator.get_size() >= chunk_size) {
                        if (!output_queue.push(std::move(ws.accumulator))) {
                            break;  // aborted
                        }
                        ws.accumulator = FragmentAccumulator();
                        ws.accumulator.reserve(chunk_size, chunk_size * 3 / 2);
                    }
                }
                // Flush remaining fragments
                if (ws.accumulator.get_size() > 0) {
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
                std::string current_qname;
                int64_t frag_id = 0;

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
                        if (!input_queue.push(std::move(current_group))) {
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
                if (!current_group.records.empty()) {
                    current_group.frag_id = frag_id++;
                    input_queue.push(std::move(current_group));
                }

                bam_destroy1(b);
                bam_hdr_destroy(hdr);
                hts_close(fp);

                // Signal workers: no more input
                input_queue.close();

                // Wait for workers to finish, then close output queue
                for (auto& w : workers) w.join();

                // Merge worker results (stats, strand, fraglen, regions)
                for (auto& ws_ptr : worker_states) {
                    WorkerState& ws = *ws_ptr;
                    stats_.merge_from(ws.stats);
                    merge_strand_obs(strand_obs_, ws.strand_obs);
                    merge_fraglen_obs(fraglen_obs_, ws.fraglen_obs);
                    if (ws.region_acc.enabled()) {
                        merge_region_acc(region_acc_, ws.region_acc);
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
        FragLenObservations& fraglen_obs = ws.fraglen_obs;
        ResolverScratch& scratch = ws.scratch;
        RegionAccumulator& region_acc = ws.region_acc;

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
        stats.n_fragments++;

        for (const auto& [r1_reads, r2_reads] : all_hits) {
            AssembledFragment frag = build_fragment(r1_reads, r2_reads);

            if (frag.exons.empty()) continue;

            RawResolveResult cr;
            bool resolved = ctx._resolve_core(
                frag.exons, frag.introns,
                frag.genomic_footprint(), cr, scratch);

            if (!resolved) {
                // Defer intergenic counting until after all hits
                // are processed to avoid multi-counting.
                if (frag.has_introns()) {
                    any_unresolved_spliced = true;
                } else {
                    any_unresolved_unspliced = true;
                }

                if (is_unique_mapper && !frag.has_introns()) {
                    int32_t flen = frag.genomic_footprint();
                    if (flen > 0) {
                        fraglen_obs.intergenic_lengths.push_back(flen);
                        stats.n_frag_length_intergenic++;
                    }
                }

                if (is_unique_mapper) {
                    int32_t ig_strand = STRAND_NONE;
                    for (const auto& eb : frag.exons) {
                        ig_strand |= eb.strand;
                    }
                    if (ig_strand == STRAND_POS || ig_strand == STRAND_NEG) {
                        strand_obs.intergenic_obs.push_back(static_cast<int8_t>(ig_strand));
                        strand_obs.intergenic_truth.push_back(static_cast<int8_t>(STRAND_POS));
                    }
                }

                // Region accumulation for intergenic fragments.
                // Every hit accumulates with 1/num_hits weight; the
                // sum across hits of one molecule is 1.0, matching
                // the 1/NH crediting of E_i (mappable_effective_length).
                if (region_acc.enabled()) {
                    region_acc.accumulate(
                        frag.exons, frag.has_introns(), is_unique_mapper,
                        num_hits);
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
                    strand_obs.exonic_spliced_obs.push_back(static_cast<int8_t>(result.exon_strand));
                    strand_obs.exonic_spliced_truth.push_back(static_cast<int8_t>(result.sj_strand));
                    stats.n_strand_trained++;
                } else if (result.splice_type != SPLICE_SPLICED_ANNOT) {
                    stats.n_strand_skipped_no_sj++;
                } else if (!result.get_is_same_strand()) {
                    stats.n_strand_skipped_ambig_strand++;
                } else {
                    stats.n_strand_skipped_ambiguous++;
                }

                if (result.get_is_same_strand() &&
                    (result.exon_strand == STRAND_POS ||
                     result.exon_strand == STRAND_NEG)) {
                    int32_t t_idx = result.get_first_t_ind();
                    if (t_idx >= 0 &&
                        t_idx < static_cast<int32_t>(ctx.t_strand_arr_.size())) {
                        int32_t t_strand = ctx.t_strand_arr_[t_idx];
                        if (t_strand == STRAND_POS ||
                            t_strand == STRAND_NEG) {
                            strand_obs.exonic_obs.push_back(static_cast<int8_t>(result.exon_strand));
                            strand_obs.exonic_truth.push_back(static_cast<int8_t>(t_strand));
                        }
                    }
                }

                int32_t ufl = result.get_unique_frag_length_mrna(
                    ctx.nrna_mask());
                if (ufl > 0) {
                    fraglen_obs.lengths.push_back(ufl);
                    fraglen_obs.splice_types.push_back(static_cast<int8_t>(result.splice_type));
                    stats.n_frag_length_unambiguous++;
                } else {
                    stats.n_frag_length_ambiguous++;
                }
            }

            // Region accumulation for resolved (non-chimeric) fragments.
            // Every resolved hit of a multimapper accumulates with
            // 1/num_hits weight; summed across all hits of a single
            // molecule the total is 1.0, matching the 1/NH crediting
            // convention of the E_i (mappable_effective_length)
            // denominator.  No count_stats gate: that would bias the
            // numerator high by a factor of NH relative to E_i.
            if (region_acc.enabled()) {
                bool frag_spliced = (result.splice_type == SPLICE_SPLICED_ANNOT
                                  || result.splice_type == SPLICE_SPLICED_UNANNOT);
                region_acc.accumulate(
                    frag.exons, frag_spliced, is_unique_mapper,
                    num_hits);
            }

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
        stats_dict["n_frag_length_unambiguous"] = stats_.n_frag_length_unambiguous;
        stats_dict["n_frag_length_ambiguous"] = stats_.n_frag_length_ambiguous;
        stats_dict["n_frag_length_intergenic"] = stats_.n_frag_length_intergenic;
        stats_dict["n_multimapper_groups"] = stats_.n_multimapper_groups;
        stats_dict["n_multimapper_alignments"] = stats_.n_multimapper_alignments;
        stats_dict["n_sj_observed"] = stats_.n_sj_observed;
        stats_dict["n_sj_blacklisted"] = stats_.n_sj_blacklisted;
        result["stats"] = stats_dict;

        // Strand observations
        nb::dict strand_dict;
        strand_dict["exonic_spliced_obs"]   = vec_to_ndarray(std::move(strand_obs_.exonic_spliced_obs));
        strand_dict["exonic_spliced_truth"] = vec_to_ndarray(std::move(strand_obs_.exonic_spliced_truth));
        strand_dict["exonic_obs"]           = vec_to_ndarray(std::move(strand_obs_.exonic_obs));
        strand_dict["exonic_truth"]         = vec_to_ndarray(std::move(strand_obs_.exonic_truth));
        strand_dict["intergenic_obs"]       = vec_to_ndarray(std::move(strand_obs_.intergenic_obs));
        strand_dict["intergenic_truth"]     = vec_to_ndarray(std::move(strand_obs_.intergenic_truth));
        result["strand_observations"] = strand_dict;

        // Fragment length observations
        nb::dict fraglen_dict;
        fraglen_dict["lengths"]            = vec_to_ndarray(std::move(fraglen_obs_.lengths));
        fraglen_dict["splice_types"]       = vec_to_ndarray(std::move(fraglen_obs_.splice_types));
        fraglen_dict["intergenic_lengths"] = vec_to_ndarray(std::move(fraglen_obs_.intergenic_lengths));
        result["frag_length_observations"] = fraglen_dict;

        // Region evidence (from gDNA calibration region accumulator)
        if (region_acc_.n_regions > 0 && !region_acc_.counts.empty()) {
            nb::dict region_dict;
            region_dict["counts"]          = vec_to_ndarray(std::move(region_acc_.counts));
            region_dict["n_regions"]       = region_acc_.n_regions;
            region_dict["fl_region_ids"]   = vec_to_ndarray(std::move(region_acc_.fl_region_ids));
            region_dict["fl_frag_lens"]    = vec_to_ndarray(std::move(region_acc_.fl_frag_lens));
            region_dict["fl_frag_strands"] = vec_to_ndarray(std::move(region_acc_.fl_frag_strands));
            result["region_evidence"] = region_dict;
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
                stamp_and_write_hit(raw_group, out, hdr, ".", ".", ".",
                                    -1, -1,
                                    0, 0.0f,
                                    "intergenic", 1, 0, "unknown", -1);
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

            for (int hit_idx = 0;
                 hit_idx < static_cast<int>(all_hits.size());
                 hit_idx++)
            {
                const auto& [r1_reads, r2_reads] = all_hits[hit_idx];

                // Gather bam1_t* for all records in this hit
                std::vector<bam1_t*> hit_raws;
                for (const auto* lr : r1_reads) {
                    ptrdiff_t idx = lr - &light_group[0];
                    hit_raws.push_back(raw_group[idx]);
                }
                for (const auto* lr : r2_reads) {
                    ptrdiff_t idx = lr - &light_group[0];
                    hit_raws.push_back(raw_group[idx]);
                }

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
                    // Intergenic / no annotation
                    stamp_and_write_hit(
                        hit_raws, out, hdr,
                        ".", ".", ".",
                        -1, -1,
                        0, 0.0f,
                        "intergenic",
                        (hit_idx == 0) ? 1 : 0,
                        0, "unknown", -1);
                    n_intergenic++;
                }

                n_records_written +=
                    static_cast<int64_t>(hit_raws.size());
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

            rec.sj_strand = read_sj_strand(b, sj_tag_mode_);

            int32_t mapped_ref_id = (rec.ref_id >= 0 &&
                rec.ref_id < static_cast<int32_t>(tid_to_ref_id_.size()))
                ? tid_to_ref_id_[rec.ref_id] : -1;
            parse_cigar(b, mapped_ref_id, rec.sj_strand,
                        rec.exons, rec.sjs);
            rec.ref_id = mapped_ref_id;

            // Apply the same splice-junction blacklist filter that the
            // scanner uses so the annotated BAM sees a consistent view
            // of the evidence.
            if (ctx_ && ctx_->has_sj_blacklist()) {
                filter_blacklisted_sjs(rec.sjs, *ctx_, mapped_ref_id);
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
             "    Number of htslib BGZF decompression threads (default 2).\n\n"
             "Returns\n"
             "-------\n"
             "dict\n"
             "    Dict with keys: 'stats', 'strand_observations',\n"
             "    'frag_length_observations'.\n")
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
