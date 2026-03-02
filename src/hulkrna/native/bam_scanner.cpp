/**
 * bam_scanner.cpp — C++ BAM scanner using htslib for hulkrna.
 *
 * Replaces the Python BAM parsing hot path (bam.py parse_bam_file +
 * fragment.py Fragment.from_reads + resolution.py resolve_fragment +
 * buffer.append) with a single C++ scan that:
 *
 *   1. Reads BAM records via htslib (no pysam overhead)
 *   2. Groups records by query name (name-sorted BAM)
 *   3. Parses CIGAR → exon blocks + splice junctions
 *   4. Builds Fragment-equivalent structures in C++
 *   5. Calls ResolveContext::_resolve_core() directly
 *   6. Appends results to NativeAccumulator
 *   7. Collects model training observations as arrays
 *   8. Returns everything to Python in one call
 *
 * Module: hulkrna._bam_impl
 *
 * Build:
 *   Part of the hulkrna scikit-build-core build — see CMakeLists.txt.
 *   Requires htslib (linked via -lhts).
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "resolve_context.h"

namespace nb = nanobind;
using namespace hulk;

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

// ================================================================
// Lightweight record for one BAM alignment in a qname group
// ================================================================

struct LightRecord {
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
    std::vector<std::pair<int32_t, int32_t>> exons;              // (start, end)
    std::vector<std::tuple<int32_t, int32_t, int32_t>> sjs;      // (start, end, strand)

    bool is_read1()         const { return flag & BAM_FREAD1; }
    bool is_read2()         const { return flag & BAM_FREAD2; }
    bool is_secondary()     const { return flag & BAM_FSECONDARY; }
    bool is_supplementary() const { return flag & BAM_FSUPPLEMENTARY; }
    bool is_reverse()       const { return flag & BAM_FREVERSE; }
    bool is_proper_pair()   const { return flag & BAM_FPROPER_PAIR; }
    bool is_mate_unmapped() const { return flag & BAM_FMUNMAP; }
};

// ================================================================
// Model training observation collectors
// ================================================================

struct StrandObservations {
    // exonic_spliced: (exon_strand, sj_strand) for strand-qualified
    std::vector<int32_t> exonic_spliced_obs;
    std::vector<int32_t> exonic_spliced_truth;

    // exonic: (exon_strand, gene_strand) for unique-gene unambiguous strand
    std::vector<int32_t> exonic_obs;
    std::vector<int32_t> exonic_truth;

    // intergenic: (intergenic_strand, POS) for unique-mapper intergenic
    std::vector<int32_t> intergenic_obs;
    std::vector<int32_t> intergenic_truth;
};

struct FragLenObservations {
    // Resolved fragments with unambiguous fragment length
    std::vector<int32_t> lengths;
    std::vector<int32_t> splice_types;  // SpliceType enum value

    // Intergenic fragment lengths (splice_type = -1 meaning None)
    std::vector<int32_t> intergenic_lengths;
};

// ================================================================
// Stats counters
// ================================================================

struct ScanStats {
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
    int64_t n_unique_gene = 0;
    int64_t n_multi_gene = 0;

    // Strand model training
    int64_t n_strand_trained = 0;
    int64_t n_strand_skipped_no_sj = 0;
    int64_t n_strand_skipped_multi_gene = 0;
    int64_t n_strand_skipped_ambiguous = 0;

    // Fragment length model training
    int64_t n_frag_length_unambiguous = 0;
    int64_t n_frag_length_ambiguous = 0;
    int64_t n_frag_length_intergenic = 0;

    // Multimapper
    int64_t n_multimapper_groups = 0;
    int64_t n_multimapper_alignments = 0;
};

// ================================================================
// Fragment-equivalent structure (built in C++, no Python object)
// ================================================================

struct NativeFragment {
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
// CIGAR parsing — ports parse_read() from bam.py
// ================================================================

static void parse_cigar(
    const bam1_t* b,
    int32_t ref_id,
    int32_t sj_strand,
    std::vector<std::pair<int32_t, int32_t>>& exons,
    std::vector<std::tuple<int32_t, int32_t, int32_t>>& sjs)
{
    exons.clear();
    sjs.clear();

    int32_t n_cigar = b->core.n_cigar;
    if (n_cigar == 0) return;

    const uint32_t* cigar = bam_get_cigar(b);
    int32_t pos = b->core.pos;  // 0-based
    int32_t start = pos;

    for (int32_t i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        if (op == CIG_REF_SKIP) {
            // End current exon, record splice junction
            if (pos > start) {
                exons.push_back({start, pos});
            }
            sjs.push_back({pos, pos + len, sj_strand});
            start = pos + len;
            pos = start;
        } else if (cigar_advances_ref(op)) {
            pos += len;
        }
        // INS, SOFT_CLIP, HARD_CLIP, PAD: don't advance ref
    }

    // Final exon
    if (pos > start) {
        exons.push_back({start, pos});
    }
}

// ================================================================
// Read SJ strand tag from a BAM record
// ================================================================

static int32_t read_sj_strand(const bam1_t* b, SJTagMode mode) {
    if (mode == SJTagMode::NONE) return STRAND_NONE;

    auto try_tag = [&](const char* tag, bool is_ts) -> int32_t {
        uint8_t* aux = bam_aux_get(b, tag);
        if (!aux) return -1;  // tag not present
        char val = bam_aux2A(aux);
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
// Build NativeFragment from a hit (r1_records, r2_records)
// ================================================================

static NativeFragment build_fragment(
    const std::vector<LightRecord*>& r1_reads,
    const std::vector<LightRecord*>& r2_reads)
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
        for (const auto& [s, e, st] : rec->sjs) {
            intron_set.insert({rec->ref_id, s, e, st});
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
        for (const auto& [s, e, st] : rec->sjs) {
            intron_set.insert({rec->ref_id, s, e, st});
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

    NativeFragment frag;
    frag.exons = std::move(merged_exons);
    frag.introns = std::move(sorted_introns);
    frag.nm = nm_total;
    return frag;
}

// ================================================================
// Hit grouping — ports _group_records_by_hit from bam.py
// ================================================================

struct HitGroupResult {
    // Hits: each is (r1_reads, r2_reads)
    std::vector<std::pair<std::vector<LightRecord*>,
                          std::vector<LightRecord*>>> hits;
    // Secondary locations for transcript-aware pairing
    std::vector<std::vector<LightRecord*>> sec_r1_locs;
    std::vector<std::vector<LightRecord*>> sec_r2_locs;
};

static HitGroupResult group_records_by_hit(
    std::vector<LightRecord>& usable)
{
    HitGroupResult result;

    // Detect HI tags
    bool has_hi = false;
    for (const auto& r : usable) {
        if (r.hi >= 0) { has_hi = true; break; }
    }

    if (has_hi) {
        // Group by HI value
        std::unordered_map<int32_t,
            std::pair<std::vector<LightRecord*>,
                      std::vector<LightRecord*>>> groups;

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
    std::vector<LightRecord*> primary_r1, primary_r2;
    std::vector<std::pair<std::vector<LightRecord*>,
                          std::vector<LightRecord*>>> singleton_hits;

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

struct ResolvedLocation {
    std::vector<LightRecord*> reads;
    std::vector<int32_t> t_inds;  // sorted transcript indices
    int32_t ref_id;
    int32_t ref_start;
};

static std::vector<std::pair<std::vector<LightRecord*>,
                              std::vector<LightRecord*>>>
pair_multimapper_reads_native(
    std::vector<std::vector<LightRecord*>>& sec_r1_locs,
    std::vector<std::vector<LightRecord*>>& sec_r2_locs,
    ResolveContext& ctx)
{
    using Hit = std::pair<std::vector<LightRecord*>,
                          std::vector<LightRecord*>>;
    std::vector<Hit> paired;

    if (sec_r1_locs.empty() && sec_r2_locs.empty())
        return paired;

    // Resolve each R1 location
    std::vector<ResolvedLocation> r1_resolved;
    r1_resolved.reserve(sec_r1_locs.size());
    for (auto& r1_reads : sec_r1_locs) {
        NativeFragment frag = build_fragment(r1_reads, {});
        std::vector<int32_t> t_inds;

        if (!frag.exons.empty()) {
            CoreResult cr;
            if (ctx._resolve_core(frag.exons, frag.introns,
                                   frag.genomic_footprint(), cr)) {
                t_inds = std::move(cr.t_inds);
            }
        }

        int32_t ref_id = r1_reads.empty() ? -1 : r1_reads[0]->ref_id;
        int32_t ref_start = r1_reads.empty() ? -1 : r1_reads[0]->ref_start;
        r1_resolved.push_back({r1_reads, std::move(t_inds), ref_id, ref_start});
    }

    // Resolve each R2 location
    std::vector<ResolvedLocation> r2_resolved;
    r2_resolved.reserve(sec_r2_locs.size());
    for (auto& r2_reads : sec_r2_locs) {
        NativeFragment frag = build_fragment({}, r2_reads);
        std::vector<int32_t> t_inds;

        if (!frag.exons.empty()) {
            CoreResult cr;
            if (ctx._resolve_core(frag.exons, frag.introns,
                                   frag.genomic_footprint(), cr)) {
                t_inds = std::move(cr.t_inds);
            }
        }

        int32_t ref_id = r2_reads.empty() ? -1 : r2_reads[0]->ref_id;
        int32_t ref_start = r2_reads.empty() ? -1 : r2_reads[0]->ref_start;
        r2_resolved.push_back({r2_reads, std::move(t_inds), ref_id, ref_start});
    }

    // STRICT — pair by transcript-set intersection
    std::unordered_set<int> r1_paired, r2_paired;

    for (int i = 0; i < static_cast<int>(r1_resolved.size()); i++) {
        if (r1_resolved[i].t_inds.empty()) continue;
        for (int j = 0; j < static_cast<int>(r2_resolved.size()); j++) {
            if (r2_resolved[j].t_inds.empty()) continue;
            if (has_intersection(r1_resolved[i].t_inds,
                                 r2_resolved[j].t_inds)) {
                paired.push_back({r1_resolved[i].reads,
                                  r2_resolved[j].reads});
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
                paired.push_back({r1_resolved[i].reads,
                                  r2_resolved[j].reads});
                r1_paired.insert(i);
                r2_paired.insert(j);
            }
        }
    }

    // CROSS-PAIR remaining
    std::vector<int> final_r1, final_r2;
    for (int i = 0; i < static_cast<int>(r1_resolved.size()); i++)
        if (!r1_paired.count(i)) final_r1.push_back(i);
    for (int j = 0; j < static_cast<int>(r2_resolved.size()); j++)
        if (!r2_paired.count(j)) final_r2.push_back(j);

    if (!final_r1.empty() && !final_r2.empty()) {
        for (int i : final_r1) {
            for (int j : final_r2) {
                paired.push_back({r1_resolved[i].reads,
                                  r2_resolved[j].reads});
                r1_paired.insert(i);
                r2_paired.insert(j);
            }
        }
    }

    // SINGLETONS
    for (int i = 0; i < static_cast<int>(r1_resolved.size()); i++)
        if (!r1_paired.count(i))
            paired.push_back({r1_resolved[i].reads, {}});
    for (int j = 0; j < static_cast<int>(r2_resolved.size()); j++)
        if (!r2_paired.count(j))
            paired.push_back({{}, r2_resolved[j].reads});

    return paired;
}

// ================================================================
// BamScanner — main scanning class
// ================================================================

class BamScanner {
public:
    ResolveContext* ctx_;
    SJTagMode sj_tag_mode_ = SJTagMode::XS_ONLY;
    bool skip_duplicates_ = true;
    bool include_multimap_ = false;

    // Mapping from htslib tid → internal ref_id (ResolveContext's numbering)
    std::vector<int32_t> tid_to_ref_id_;

    // Results
    ScanStats stats_;
    NativeAccumulator accumulator_;
    StrandObservations strand_obs_;
    FragLenObservations fraglen_obs_;

    BamScanner(ResolveContext& ctx,
               const std::string& sj_tag_spec,
               bool skip_duplicates,
               bool include_multimap)
        : ctx_(&ctx),
          skip_duplicates_(skip_duplicates),
          include_multimap_(include_multimap)
    {
        // Parse SJ tag spec
        if (sj_tag_spec == "XS")
            sj_tag_mode_ = SJTagMode::XS_ONLY;
        else if (sj_tag_spec == "ts")
            sj_tag_mode_ = SJTagMode::TS_ONLY;
        else if (sj_tag_spec == "XS,ts")
            sj_tag_mode_ = SJTagMode::XS_TS;
        else if (sj_tag_spec == "ts,XS")
            sj_tag_mode_ = SJTagMode::TS_XS;
        else if (sj_tag_spec == "none" || sj_tag_spec.empty())
            sj_tag_mode_ = SJTagMode::NONE;
        else
            sj_tag_mode_ = SJTagMode::XS_TS;  // default
    }

    // ----------------------------------------------------------------
    // Main scan entry point
    // ----------------------------------------------------------------

    nb::dict scan(const std::string& bam_path) {
        // Open BAM file
        htsFile* fp = hts_open(bam_path.c_str(), "rb");
        if (!fp) {
            throw std::runtime_error("Failed to open BAM file: " + bam_path);
        }

        bam_hdr_t* hdr = sam_hdr_read(fp);
        if (!hdr) {
            hts_close(fp);
            throw std::runtime_error("Failed to read BAM header: " + bam_path);
        }

        // Build tid → ref_id mapping
        build_tid_mapping(hdr);

        bam1_t* b = bam_init1();
        if (!b) {
            bam_hdr_destroy(hdr);
            hts_close(fp);
            throw std::runtime_error("Failed to allocate bam1_t");
        }

        // Query-name grouping: collect records for same qname
        std::vector<LightRecord> current_group;
        std::string current_qname;
        int64_t frag_id = 0;

        while (sam_read1(fp, hdr, b) >= 0) {
            stats_.total++;

            // Apply early filters
            uint16_t flag = b->core.flag;

            if (flag & BAM_FQCFAIL) {
                stats_.qc_fail++;
                continue;
            }
            if (flag & BAM_FUNMAP) {
                stats_.unmapped++;
                continue;
            }
            if (flag & BAM_FDUP) {
                stats_.duplicate++;
                if (skip_duplicates_) continue;
            }

            // Enforce paired-end
            if (!(flag & BAM_FPAIRED)) {
                bam_destroy1(b);
                bam_hdr_destroy(hdr);
                hts_close(fp);
                throw std::runtime_error(
                    "Input BAM must be paired-end, but found unpaired read");
            }

            // Count secondary/supplementary
            if (flag & BAM_FSECONDARY) stats_.secondary++;
            if (flag & BAM_FSUPPLEMENTARY) stats_.supplementary++;

            // Get query name
            const char* qname = bam_get_qname(b);

            // Check if we've moved to a new query name group
            if (!current_group.empty() && current_qname != qname) {
                process_qname_group(current_group, frag_id);
                current_group.clear();
                frag_id++;
            }

            current_qname = qname;

            // Build LightRecord from htslib bam1_t
            LightRecord rec;
            rec.ref_id = b->core.tid;
            rec.ref_start = b->core.pos;
            rec.mate_ref_id = b->core.mtid;
            rec.mate_ref_start = b->core.mpos;
            rec.flag = flag;

            // Read tags
            rec.nm = 0;
            uint8_t* nm_aux = bam_aux_get(b, "NM");
            if (nm_aux) rec.nm = bam_aux2i(nm_aux);

            rec.nh = 1;
            uint8_t* nh_aux = bam_aux_get(b, "NH");
            if (nh_aux) rec.nh = bam_aux2i(nh_aux);

            rec.hi = -1;
            uint8_t* hi_aux = bam_aux_get(b, "HI");
            if (hi_aux) rec.hi = bam_aux2i(hi_aux);

            // Read SJ strand tag
            rec.sj_strand = read_sj_strand(b, sj_tag_mode_);

            // Parse CIGAR
            int32_t mapped_ref_id = (rec.ref_id >= 0 &&
                rec.ref_id < static_cast<int32_t>(tid_to_ref_id_.size()))
                ? tid_to_ref_id_[rec.ref_id] : -1;
            parse_cigar(b, mapped_ref_id, rec.sj_strand,
                        rec.exons, rec.sjs);

            // Map htslib ref_id to internal ref_id for the record
            rec.ref_id = mapped_ref_id;

            current_group.push_back(std::move(rec));
        }

        // Process last group
        if (!current_group.empty()) {
            process_qname_group(current_group, frag_id);
        }

        // Clean up htslib
        bam_destroy1(b);
        bam_hdr_destroy(hdr);
        hts_close(fp);

        // Build result dict
        return build_result();
    }

private:
    // ----------------------------------------------------------------
    // Build tid → internal ref_id mapping
    // ----------------------------------------------------------------

    void build_tid_mapping(const bam_hdr_t* hdr) {
        int32_t n_targets = hdr->n_targets;
        tid_to_ref_id_.resize(n_targets, -1);

        for (int32_t i = 0; i < n_targets; i++) {
            const char* name = hdr->target_name[i];
            auto it = ctx_->ref_to_id_.find(name);
            if (it != ctx_->ref_to_id_.end()) {
                tid_to_ref_id_[i] = it->second;
            }
            // Refs not in the index get -1 (intergenic)
        }
    }

    // ----------------------------------------------------------------
    // Process one query-name group
    // ----------------------------------------------------------------

    void process_qname_group(
        std::vector<LightRecord>& records,
        int64_t frag_id)
    {
        if (records.empty()) return;

        stats_.n_read_names++;

        // Determine NH from first record (they should all agree)
        int32_t nh = 1;
        for (const auto& r : records) {
            if (r.nh > 1) { nh = r.nh; break; }
        }

        // Detect multimap
        bool has_secondary = false;
        for (const auto& r : records) {
            if (r.is_secondary()) { has_secondary = true; break; }
        }
        bool is_multimap = (nh > 1) || has_secondary;

        if (is_multimap) {
            stats_.multimapping++;
            if (!include_multimap_) return;
        } else {
            stats_.unique++;
        }

        // Group into hits + secondary locations
        auto hit_result = group_records_by_hit(records);

        // Track hit-level mate stats (primary hits only)
        for (const auto& [r1s, r2s] : hit_result.hits) {
            if (r1s.empty() || r2s.empty()) {
                stats_.mate_unmapped++;
            } else {
                // Check proper-pair on primary R1
                bool found_proper = false;
                for (const auto* r : r1s) {
                    if (!r->is_supplementary() && !r->is_secondary()) {
                        found_proper = r->is_proper_pair();
                        break;
                    }
                }
                if (found_proper) stats_.proper_pair++;
                else stats_.improper_pair++;
            }
        }

        // Handle multimapper secondary pairing
        std::vector<std::pair<std::vector<LightRecord*>,
                              std::vector<LightRecord*>>> secondary_pairs;
        if (!hit_result.sec_r1_locs.empty() ||
            !hit_result.sec_r2_locs.empty()) {
            secondary_pairs = pair_multimapper_reads_native(
                hit_result.sec_r1_locs,
                hit_result.sec_r2_locs,
                *ctx_);
        }

        // Combine all hits
        auto& all_hits = hit_result.hits;
        all_hits.insert(all_hits.end(),
                        secondary_pairs.begin(), secondary_pairs.end());

        int32_t num_hits = std::max(nh, static_cast<int32_t>(all_hits.size()));
        bool is_unique_mapper = (num_hits == 1);
        int32_t n_buffered_mm = 0;

        for (const auto& [r1_reads, r2_reads] : all_hits) {
            NativeFragment frag = build_fragment(r1_reads, r2_reads);
            stats_.n_fragments++;

            if (frag.exons.empty()) continue;

            // Resolve fragment
            CoreResult cr;
            bool resolved = ctx_->_resolve_core(
                frag.exons, frag.introns,
                frag.genomic_footprint(), cr);

            if (!resolved) {
                // Intergenic
                if (frag.has_introns()) {
                    stats_.n_intergenic_spliced++;
                } else {
                    stats_.n_intergenic_unspliced++;
                }

                // Train intergenic fragment length (unique mappers, unspliced)
                if (is_unique_mapper && !frag.has_introns()) {
                    int32_t flen = frag.genomic_footprint();
                    if (flen > 0) {
                        fraglen_obs_.intergenic_lengths.push_back(flen);
                        stats_.n_frag_length_intergenic++;
                    }
                }

                // Train intergenic strand model
                if (is_unique_mapper) {
                    int32_t ig_strand = STRAND_NONE;
                    for (const auto& eb : frag.exons) {
                        ig_strand |= eb.strand;
                    }
                    if (ig_strand == STRAND_POS || ig_strand == STRAND_NEG) {
                        strand_obs_.intergenic_obs.push_back(ig_strand);
                        strand_obs_.intergenic_truth.push_back(STRAND_POS);
                    }
                }
                continue;
            }

            // Build ResolvedResult for accumulator
            ResolvedResult result = ResolvedResult::from_core(cr);
            result.num_hits = num_hits;
            result.nm = frag.nm;

            // --- Handle chimeric fragments ---
            if (result.get_is_chimeric()) {
                stats_.n_chimeric++;
                if (result.chimera_type == CHIMERA_TRANS)
                    stats_.n_chimeric_trans++;
                else if (result.chimera_type == CHIMERA_CIS_STRAND_SAME)
                    stats_.n_chimeric_cis_strand_same++;
                else if (result.chimera_type == CHIMERA_CIS_STRAND_DIFF)
                    stats_.n_chimeric_cis_strand_diff++;
                accumulator_.append(result, frag_id);
                continue;
            }

            // --- Stats: overlap type ---
            stats_.n_with_exon++;
            if (result.splice_type == SPLICE_SPLICED_ANNOT) {
                stats_.n_with_annotated_sj++;
            } else if (result.splice_type == SPLICE_SPLICED_UNANNOT) {
                stats_.n_with_unannotated_sj++;
            }

            if (result.get_is_unique_gene()) {
                stats_.n_unique_gene++;
            } else {
                stats_.n_multi_gene++;
            }

            if (is_unique_mapper) {
                // --- Train strand models ---
                if (result.get_is_strand_qualified()) {
                    strand_obs_.exonic_spliced_obs.push_back(result.exon_strand);
                    strand_obs_.exonic_spliced_truth.push_back(result.sj_strand);
                    stats_.n_strand_trained++;
                } else if (result.splice_type != SPLICE_SPLICED_ANNOT) {
                    stats_.n_strand_skipped_no_sj++;
                } else if (!result.get_is_unique_gene()) {
                    stats_.n_strand_skipped_multi_gene++;
                } else {
                    stats_.n_strand_skipped_ambiguous++;
                }

                // Train exonic fallback model
                if (result.get_is_unique_gene() &&
                    (result.exon_strand == STRAND_POS ||
                     result.exon_strand == STRAND_NEG)) {
                    int32_t t_idx = result.get_first_t_ind();
                    if (t_idx >= 0 &&
                        t_idx < ctx_->n_transcripts_) {
                        int32_t gene_idx = ctx_->t_to_g_arr_[t_idx];
                        if (gene_idx >= 0 && gene_idx <
                            static_cast<int32_t>(ctx_->g_to_strand_arr_.size())) {
                            int32_t gene_strand = ctx_->g_to_strand_arr_[gene_idx];
                            if (gene_strand == STRAND_POS ||
                                gene_strand == STRAND_NEG) {
                                strand_obs_.exonic_obs.push_back(result.exon_strand);
                                strand_obs_.exonic_truth.push_back(gene_strand);
                            }
                        }
                    }
                }

                // Train fragment length model
                int32_t ufl = result.get_unique_frag_length();
                if (ufl > 0) {
                    fraglen_obs_.lengths.push_back(ufl);
                    fraglen_obs_.splice_types.push_back(result.splice_type);
                    stats_.n_frag_length_unambiguous++;
                } else {
                    stats_.n_frag_length_ambiguous++;
                }
            }

            // Buffer resolved fragment
            accumulator_.append(result, frag_id);

            if (!is_unique_mapper) {
                n_buffered_mm++;
            }
        }

        if (n_buffered_mm > 0) {
            stats_.n_multimapper_groups++;
            stats_.n_multimapper_alignments += n_buffered_mm;
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
        stats_dict["n_unique_gene"] = stats_.n_unique_gene;
        stats_dict["n_multi_gene"] = stats_.n_multi_gene;
        stats_dict["n_strand_trained"] = stats_.n_strand_trained;
        stats_dict["n_strand_skipped_no_sj"] = stats_.n_strand_skipped_no_sj;
        stats_dict["n_strand_skipped_multi_gene"] = stats_.n_strand_skipped_multi_gene;
        stats_dict["n_strand_skipped_ambiguous"] = stats_.n_strand_skipped_ambiguous;
        stats_dict["n_frag_length_unambiguous"] = stats_.n_frag_length_unambiguous;
        stats_dict["n_frag_length_ambiguous"] = stats_.n_frag_length_ambiguous;
        stats_dict["n_frag_length_intergenic"] = stats_.n_frag_length_intergenic;
        stats_dict["n_multimapper_groups"] = stats_.n_multimapper_groups;
        stats_dict["n_multimapper_alignments"] = stats_.n_multimapper_alignments;
        result["stats"] = stats_dict;

        // Strand observations
        nb::dict strand_dict;
        strand_dict["exonic_spliced_obs"] = std::move(strand_obs_.exonic_spliced_obs);
        strand_dict["exonic_spliced_truth"] = std::move(strand_obs_.exonic_spliced_truth);
        strand_dict["exonic_obs"] = std::move(strand_obs_.exonic_obs);
        strand_dict["exonic_truth"] = std::move(strand_obs_.exonic_truth);
        strand_dict["intergenic_obs"] = std::move(strand_obs_.intergenic_obs);
        strand_dict["intergenic_truth"] = std::move(strand_obs_.intergenic_truth);
        result["strand_observations"] = strand_dict;

        // Fragment length observations
        nb::dict fraglen_dict;
        fraglen_dict["lengths"] = std::move(fraglen_obs_.lengths);
        fraglen_dict["splice_types"] = std::move(fraglen_obs_.splice_types);
        fraglen_dict["intergenic_lengths"] = std::move(fraglen_obs_.intergenic_lengths);
        result["frag_length_observations"] = fraglen_dict;

        // Accumulator size
        result["accumulator_size"] = accumulator_.get_size();

        return result;
    }

public:
    // Expose accumulator for finalize
    NativeAccumulator& get_accumulator() { return accumulator_; }

    nb::dict finalize_accumulator(const std::vector<int64_t>& t_to_g_arr) {
        return accumulator_.finalize(t_to_g_arr);
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
    m.doc() = "C++ BAM scanner for hulkrna using htslib (nanobind).\n\n"
              "Replaces the Python BAM parsing hot path with a single\n"
              "C++ scan: BAM → resolve → buffer in one call.";

    nb::class_<BamScanner>(m, "BamScanner")
        .def(nb::init<ResolveContext&, const std::string&, bool, bool>(),
             nb::arg("ctx"),
             nb::arg("sj_tag_spec"),
             nb::arg("skip_duplicates") = true,
             nb::arg("include_multimap") = false,
             "Create a BAM scanner.\n\n"
             "Parameters\n"
             "----------\n"
             "ctx : ResolveContext\n"
             "    The resolve context from index building.\n"
             "sj_tag_spec : str\n"
             "    SJ strand tag specification: 'XS', 'ts', 'XS,ts', 'ts,XS', or 'none'.\n"
             "skip_duplicates : bool\n"
             "    Skip duplicate-flagged reads (default True).\n"
             "include_multimap : bool\n"
             "    Include multimapping reads (default False).\n",
             nb::keep_alive<1, 2>()  // BamScanner must keep ResolveContext alive
        )
        .def("scan", &BamScanner::scan,
             nb::arg("bam_path"),
             "Scan BAM file and return results dict.\n\n"
             "Returns a dict with keys:\n"
             "  'stats' — dict of all counters\n"
             "  'strand_observations' — dict of strand training arrays\n"
             "  'frag_length_observations' — dict of frag length arrays\n"
             "  'accumulator_size' — number of resolved fragments\n")
        .def("finalize_accumulator", &BamScanner::finalize_accumulator,
             nb::arg("t_to_g_arr"),
             "Finalize the internal accumulator to raw bytes dict.\n\n"
             "Parameters\n"
             "----------\n"
             "t_to_g_arr : list[int]\n"
             "    Transcript-to-gene mapping array (int64).\n")
        ;

    m.def("detect_sj_strand_tag", &detect_sj_strand_tag_native,
          nb::arg("bam_path"),
          nb::arg("max_spliced_reads") = 1000,
          "Auto-detect SJ strand tags in a BAM file.\n\n"
          "Returns a tag specification string: 'XS', 'ts', 'XS,ts', or 'none'.");
}
