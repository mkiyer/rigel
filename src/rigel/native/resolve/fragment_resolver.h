// fragment_resolver.h — FragmentResolver: the annotation as the scan reads it. Built once from the index; then, per
// hit, the pair rule (are the two mates one molecule?) and the resolution (which transcripts could the fragment be
// from, what are its lengths on them, and what could its unsequenced gaps hold?).
#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C" {
#include "cgranges.h"
}

#include "fragment_blocks.h"
#include "resolve/resolved_fragment.h"
#include "resolve/resolver_scratch.h"

namespace rigel {

class FragmentResolver {
public:
    FragmentResolver() = default;
    ~FragmentResolver() {
        if (cr_) cr_destroy(cr_);
    }
    FragmentResolver(const FragmentResolver&) = delete;
    FragmentResolver& operator=(const FragmentResolver&) = delete;

    // ── the index, built once at load ──────────────────────────────────────────────────────────────────────────

    /// The reference-ID space, in canonical order (index.ref_names), so the resolver's ref ids match
    /// index.ref_name_to_id, the space t_df, region_df, RegionArrays and the calibration accumulator partition use.
    /// Call it first: the builders below name only references it installed, and a BAM reference it did not install
    /// maps to -1, outside the index.
    void set_ref_names(const std::vector<std::string>& ref_names);

    /// The overlap index over the collapsed intervals: each interval's type (`ITYPE_*`) and its transcripts, CSR.
    /// An interval on a reference `set_ref_names` did not install is an index defect, and throws.
    void build_overlap_index(const std::vector<std::string>& refs, const std::vector<int32_t>& starts,
                             const std::vector<int32_t>& ends, const std::vector<int32_t>& iv_types,
                             const std::vector<int32_t>& tset_data, const std::vector<int32_t>& tset_offsets);

    /// The annotated junctions, keyed by (ref, start, end, strand), each with the transcripts holding it, CSR. A
    /// junction on a reference `set_ref_names` did not install is an index defect, and throws.
    void build_sj_map(const std::vector<std::string>& refs, const std::vector<int32_t>& starts,
                      const std::vector<int32_t>& ends, const std::vector<int32_t>& strands,
                      const std::vector<int32_t>& tset_data, const std::vector<int32_t>& tset_offsets);

    /// Each transcript's exons, CSR: its introns, for its path through a fragment's gaps and its length.
    void build_exon_index(const std::vector<int32_t>& offsets, const std::vector<int32_t>& starts,
                          const std::vector<int32_t>& ends);

    /// The transcript count, which sizes each scan worker's scratch.
    void set_metadata(int32_t n_transcripts) { n_transcripts_ = n_transcripts; }

    /// Each transcript's strand, looked up directly.
    void set_transcript_strands(const std::vector<int32_t>& t_strand) { t_strand_arr_ = t_strand; }

    /// Each transcript's synthetic nascent span (-1 = none). The spans stay out of the overlap index: `resolve` adds
    /// each annotated candidate's span to the fragment's candidates.
    void set_nrna_parent_index(const std::vector<int32_t>& nrna_parent) { nrna_parent_ = nrna_parent; }

    /// The library's fragment-length limit (`BamScanConfig.max_frag_length`): no pair whose molecule is longer, along
    /// the genome and along every witness transcript, is one (`validate_pair`).
    void set_max_fragment_length(int32_t L);

    /// The intron-overhang tolerance K (bp): how many bases an aligned read end may lie inside an annotated intron, at
    /// the edge of an unsequenced gap or at a fragment's outer end. An annotated intron in an unaligned gap may overlap
    /// up to K aligned bases and is still cut whole, and only an intron holding more than K aligned bases contradicts a
    /// transcript. 0 is the default; a negative K is refused.
    void set_intron_overhang_tolerance(int32_t K) {
        if (K < 0) throw std::invalid_argument("set_intron_overhang_tolerance: K must be >= 0, got " + std::to_string(K));
        intron_overhang_tolerance_ = K;
    }

    // ── what the scan reads ────────────────────────────────────────────────────────────────────────────────────

    /// The index's id of a reference name, or -1 for a reference outside the index.
    int32_t ref_id(const char* name) const;
    int32_t n_transcripts() const { return n_transcripts_; }
    /// Transcript t's strand, STRAND_NONE for an index the strand table does not hold.
    int32_t transcript_strand(int32_t t) const;

    // ── the pair rule: two mates are one molecule, or why not ─────────────────────────────────────────────────

    /// Whether the forward mate ``f`` and the reverse mate ``r`` on ``ref_id`` can be one molecule: `HIT_USABLE`,
    /// `HIT_OUTWARD` or `HIT_LENGTH` (the caller has checked the reference and the strands).
    ///
    /// Walking the forward mate's 5' soft clip back and the reverse mate's forward gives the 5' ends a and b, which
    /// must face each other (a < b). The molecule is [lo, hi), from the lower of a and the reverse mate's first aligned
    /// base to the higher of b and the forward mate's last (a 5'-trimmed read starts inside its molecule, and its mate
    /// may run past it). Its length, [lo, hi) less every observed skip of either mate, must lie in
    /// [1, max_fragment_length]. Failing that, one annotated transcript that overlaps both mates and holds every
    /// observed junction (`pair_holders`) may witness a shorter path: its other introns in the stretches neither mate
    /// aligned are cut by the K rule, and it offers none when it contradicts the blocks. One transcript supplies both
    /// the observed junctions and the implied introns, so introns from two transcripts never combine into one path.
    int32_t validate_pair(int32_t ref_id, const PairMate& f, const PairMate& r,
                          const std::vector<IntronBlock>& junctions, ResolverScratch& scratch) const;

    /// The transcripts that explain a pair, into ``holders``: each holds every junction in ``junctions``
    /// (`junction_holders`) and overlaps both mates' aligned spans, [f_lo, f_hi) and [r_lo, r_hi). `HIT_USABLE` when
    /// there are any; otherwise `junction_holders`' reason, or `HIT_UNSUPPORTED` when every holder misses a mate.
    int32_t pair_holders(const std::vector<IntronBlock>& junctions, int32_t f_lo, int32_t f_hi, int32_t r_lo,
                         int32_t r_hi, ResolverScratch& scratch, std::vector<int32_t>& holders) const;

    // ── resolution ─────────────────────────────────────────────────────────────────────────────────────────────

    /// Resolve one fragment into ``out``: its candidate transcripts (a candidate need not be compatible with the
    /// alignment: the merge falls back to a union when the per-block sets share no transcript), each one's length and
    /// exonic bases, its splice type and strands, and every explanation of its unsequenced gaps. ``exons`` is
    /// non-empty; ``introns`` are its observed junctions as `sort_unique_junctions` leaves them, every one held by one
    /// annotated transcript (the scan's annotated-path rule). Sets every field of ``out`` but the hit's own
    /// (num_hits, nm, clip).
    void resolve(const std::vector<ExonBlock>& exons, const std::vector<IntronBlock>& introns,
                 int32_t genomic_footprint, ResolvedFragment& out, ResolverScratch& scratch) const;

private:
    /// Transcript t's introns: `n` (start, end) pairs, the gaps between its consecutive exons. None for a transcript
    /// the exon index does not hold, or one with a single exon.
    struct TranscriptIntrons {
        const int32_t* starts;
        const int32_t* ends;
        int32_t n;
    };

    /// The id `set_ref_names` gave a reference a builder names; throws for one it did not install.
    int32_t installed_ref_id(const std::string& ref, const char* builder) const;
    bool in_exon_index(int32_t t) const { return t >= 0 && t + 1 < static_cast<int32_t>(exon_offsets_.size()); }
    TranscriptIntrons introns_of(int32_t t) const;

    /// Does transcript ``t`` contradict the alignment? True when one of its introns overlaps the interior of an aligned
    /// block by more than K bases: the read has bases there, so the molecule cannot have spliced them out.
    bool transcript_contradicts_blocks(int32_t t, const std::vector<std::pair<int32_t, int32_t>>& blocks,
                                       int32_t K) const;
    /// Append every intron of ``t`` lying inside ``gap`` (±K), in genomic order, with ``strand``.
    void collect_transcript_introns_in_gap(int32_t t, const GapBlock& gap, int32_t K, int32_t strand,
                                           std::vector<IntronBlock>& out) const;
    /// Whether transcript ``t``'s span overlaps both [f_lo, f_hi) and [r_lo, r_hi).
    bool overlaps_both(int32_t t, int32_t f_lo, int32_t f_hi, int32_t r_lo, int32_t r_hi) const;
    void retain_overlapping_both(std::vector<int32_t>& ts, int32_t f_lo, int32_t f_hi, int32_t r_lo,
                                 int32_t r_hi) const;
    /// Every annotated transcript whose span overlaps [start, end) on ``ref_id``, sorted and unique.
    void transcripts_overlapping(int32_t ref_id, int32_t start, int32_t end, ResolverScratch& scratch,
                                 std::vector<int32_t>& out) const;

    /// The transcripts holding one junction, sorted and unique; a junction with no definite strand matches either.
    void sj_lookup_into(int32_t ref_id, int32_t start, int32_t end, int32_t strand, std::vector<int32_t>& out) const;
    int32_t junction_holders(const std::vector<IntronBlock>& junctions, ResolverScratch& scratch,
                             std::vector<int32_t>& holders) const;

    // The steps of `resolve`, in order
    void query_blocks(const std::vector<ExonBlock>& exons, ResolvedFragment& out, ResolverScratch& scratch) const;
    void add_nascent_spans(const std::vector<ExonBlock>& exons, ResolverScratch& scratch) const;
    void read_junctions(const std::vector<IntronBlock>& introns, ResolvedFragment& out,
                        ResolverScratch& scratch) const;
    void choose_candidates(bool spliced, ResolvedFragment& out, ResolverScratch& scratch) const;
    void candidate_exon_bp(ResolvedFragment& out, ResolverScratch& scratch) const;
    bool candidates_on_both_strands(const std::vector<int32_t>& candidates) const;
    void enumerate_gap_hypotheses(const std::vector<ExonBlock>& exons, const std::vector<IntronBlock>& observed_introns,
                                  bool certified_rna, ResolvedFragment& out, ResolverScratch& scratch) const;
    void unsequenced_gaps(const std::vector<ExonBlock>& exons, const std::vector<IntronBlock>& observed_introns,
                          ResolverScratch& scratch) const;

    // The overlap index: cgranges over the collapsed intervals, each with its type and transcripts (CSR)
    cgranges_t* cr_ = nullptr;
    std::vector<int8_t> iv_type_;
    std::vector<int32_t> t_set_data_;
    std::vector<int32_t> t_set_offsets_;

    // The annotated junctions: each key's slice (offset, count) of sj_map_data_
    std::unordered_map<SJKey, std::pair<int32_t, int32_t>, SJKeyHash> sj_map_;
    std::vector<int32_t> sj_map_data_;

    int32_t n_transcripts_ = 0;
    std::vector<int32_t> t_strand_arr_;
    std::vector<int32_t> nrna_parent_;

    // Each transcript's exons, CSR: [n_transcripts + 1] offsets into the genomic starts and ends
    std::vector<int32_t> exon_offsets_;
    std::vector<int32_t> exon_starts_;
    std::vector<int32_t> exon_ends_;

    int32_t intron_overhang_tolerance_ = 0;
    int32_t max_fragment_length_ = 0;

    std::unordered_map<std::string, int32_t> ref_to_id_;
    std::vector<std::string> id_to_ref_;
};

}  // namespace rigel
