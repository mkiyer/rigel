// bam_scanner.cpp — the scan: the partition it deposits onto, the reader, the workers and the per-template work.
#include "scan/bam_scanner.h"

#include <algorithm>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>

namespace nb = nanobind;

namespace rigel {

using accumulator::AccumulatorSet;

// ================================================================
// The stats and the strand observations
// ================================================================

void BamScanStats::count_flags(uint16_t flag) {
    total++;
    qc_fail += (flag & BAM_FQCFAIL) != 0;
    unmapped += (flag & BAM_FUNMAP) != 0;
    duplicate += (flag & BAM_FDUP) != 0;
    secondary += (flag & BAM_FSECONDARY) != 0;
    supplementary += (flag & BAM_FSUPPLEMENTARY) != 0;
    proper_pair += (flag & (BAM_FPROPER_PAIR | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) == BAM_FPROPER_PAIR;
}

void BamScanStats::merge_from(const BamScanStats& o) {
    total += o.total;
    qc_fail += o.qc_fail;
    unmapped += o.unmapped;
    secondary += o.secondary;
    supplementary += o.supplementary;
    duplicate += o.duplicate;
    proper_pair += o.proper_pair;
    pairs.merge_from(o.pairs);
    records.merge_from(o.records);
    n_fragments += o.n_fragments;
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
    n_nm_conflicting += o.n_nm_conflicting;
    n_unique_unspliced_mismatches += o.n_unique_unspliced_mismatches;
    n_unique_unspliced_bases += o.n_unique_unspliced_bases;
    for (int32_t v = 0; v < N_HIT_VERDICTS; v++) hits_rejected[v] += o.hits_rejected[v];
    n_nh_above_reported += o.n_nh_above_reported;
    n_read_pairs_lacking_mate += o.n_read_pairs_lacking_mate;
    for (size_t i = 0; i < splice_census.size(); i++) splice_census[i] += o.splice_census[i];
}

void StrandObservations::merge_from(StrandObservations& o) {
    exonic.merge_from(o.exonic);
    for (const auto& [key, counts] : o.sj_strand_table) {
        auto& into = sj_strand_table[key];
        into.n_sense += counts.n_sense;
        into.n_antisense += counts.n_antisense;
    }
    o.sj_strand_table.clear();
}

// ================================================================
// The partition
// ================================================================

std::unique_ptr<AccumulatorSet> ScanPartition::make_set() const {
    auto set = std::make_unique<AccumulatorSet>(
        region_bounds.data(), region_bounds.size(), ref_region_bound_offsets.data(),
        ref_region_bound_offsets.empty() ? 0 : ref_region_bound_offsets.size() - 1,
        region_types.empty() ? nullptr : region_types.data(), region_types.size(), max_length);
    if (has_sj) {
        set->set_sj(sj_offsets.data(), sj_offsets.size(), sj_boundary_right.data(), sj_strand.data(),
                    sj_boundary_right.size(), ref_region_bound_offsets.data());
    }
    return set;
}

BamScanner::BamScanner(const FragmentResolver& resolver, const std::string& sj_tag_spec, bool skip_duplicates,
                       bool include_multimap)
    : resolver_(&resolver),
      sj_tag_mode_(parse_sj_tag_spec(sj_tag_spec)),
      skip_duplicates_(skip_duplicates),
      include_multimap_(include_multimap) {}

void BamScanner::set_regions(Array1<int64_t> region_bounds, Array1<int64_t> ref_region_bound_offsets,
                             int32_t n_refs, Array1<uint8_t> region_types, int max_length) {
    if (n_refs < 0) throw std::invalid_argument("set_regions: n_refs must be non-negative");
    const size_t n_off = ref_region_bound_offsets.shape(0);
    if (n_off != static_cast<size_t>(n_refs) + 1) {
        throw std::invalid_argument("set_regions: ref_region_bound_offsets must have length n_refs + 1");
    }
    if (acc_set_) {
        throw std::runtime_error("set_regions: regions already set on this BamScanner; create a new instance");
    }
    // Validated at the boundary as `BamScanConfig.max_frag_length` is: it gates L and sizes the pool histograms, so
    // at 0 every fragment would drop as too long.
    if (max_length < 1) {
        throw std::invalid_argument(
            "set_regions: max_length must be >= 1, got " + std::to_string(max_length) +
            ". It is the fragment-length limit applied to L as well as the pool-histogram width, so "
            "at 0 every real fragment is dropped as too long and the tally is silently empty.");
    }
    partition_.region_bounds.assign(region_bounds.data(), region_bounds.data() + region_bounds.shape(0));
    partition_.ref_region_bound_offsets.assign(ref_region_bound_offsets.data(), ref_region_bound_offsets.data() + n_off);
    partition_.region_types.assign(region_types.data(), region_types.data() + region_types.shape(0));
    partition_.max_length = max_length;
    acc_set_ = partition_.make_set();
}

void BamScanner::set_sj(Array1<int64_t> offsets, Array1<int64_t> boundary_right, Array1<int8_t> sj_strand) {
    if (!acc_set_) {
        throw std::runtime_error(
            "set_sj: call set_regions first — the sj CSR is keyed by region_bound index, so "
            "there is nothing to key it against until the partition is installed");
    }
    if (boundary_right.shape(0) != sj_strand.shape(0)) {
        throw std::invalid_argument("set_sj: boundary_right has " + std::to_string(boundary_right.shape(0)) +
                                    " entries but sj_strand has " + std::to_string(sj_strand.shape(0)));
    }
    partition_.sj_offsets.assign(offsets.data(), offsets.data() + offsets.shape(0));
    partition_.sj_boundary_right.assign(boundary_right.data(), boundary_right.data() + boundary_right.shape(0));
    partition_.sj_strand.assign(sj_strand.data(), sj_strand.data() + sj_strand.shape(0));
    partition_.has_sj = true;
    acc_set_->set_sj(partition_.sj_offsets.data(), partition_.sj_offsets.size(), partition_.sj_boundary_right.data(),
                     partition_.sj_strand.data(), partition_.sj_boundary_right.size(),
                     partition_.ref_region_bound_offsets.data());
}

// ================================================================
// The worker: one template at a time
// ================================================================

namespace {

/// A read pair's eligible records, parsed by the reader, with its frag_id and name hash: no raw record crosses a
/// thread.
struct QnameGroup {
    std::vector<ParsedAlignment> records;
    int64_t frag_id = 0;
    size_t name_hash = 0;  // `std::hash` of the read name, for the digest both passes fold
};

/// Pre-size a worker's chunk: chunk_size fragment rows and 3/2 candidate-transcript rows per fragment. A capacity hint
/// only; the columns grow past it.
void reserve_chunk(FragmentChunk& chunk, int64_t chunk_size) { chunk.reserve(chunk_size, chunk_size * 3 / 2); }

/// The introns a hit's reads stated, appended to ``out`` each once: a pair where only one mate carries the strand tag
/// yields the same intron twice under the junctions' (ref, start, end, strand) key, and the accumulator's held records
/// keep each once. They must not be merged by the junctions' key instead, which the resolver shares: a merged junction
/// may be one no transcript holds.
void append_observed_introns(const std::vector<IntronBlock>& introns, std::vector<IntronBlock>& out) {
    const size_t first = out.size();
    for (const auto& intron : introns) {
        if (out.size() > first && out.back().start == intron.start && out.back().end == intron.end) continue;
        out.push_back(intron);
    }
}

}  // namespace

struct BamScanner::QnameBatch {
    std::vector<QnameGroup> groups;
};

/// One worker's state and its work on one template. It shares nothing mutable with the scanner or the other workers:
/// its stats, strand observations and accumulator set are merged once the BAM is read.
class BamScanner::Worker {
public:
    Worker(const FragmentResolver& resolver, const ScanPartition& partition, bool include_multimap, int64_t chunk_size)
        : resolver_(resolver),
          include_multimap_(include_multimap),
          chunk_size_(chunk_size),
          scratch_(resolver.n_transcripts()),
          acc_set_(partition.make_set()) {
        reserve_chunk(chunk_, chunk_size_);
    }

    /// Take batches until the input closes, pushing each full chunk to ``output``; the last, partial chunk too unless
    /// the output was aborted.
    void run(BoundedQueue<QnameBatch>& input, BoundedQueue<FragmentChunk>& output) {
        QnameBatch batch;
        while (input.pop(batch)) {
            for (QnameGroup& group : batch.groups) {
                scan_template(group);
                if (chunk_.size() < chunk_size_) continue;
                if (!output.push(std::move(chunk_))) return;
                chunk_ = FragmentChunk();
                reserve_chunk(chunk_, chunk_size_);
            }
            batch.groups.clear();
        }
        if (chunk_.size() > 0) output.push(std::move(chunk_));
    }

    // What the worker found, merged into the scanner's once the BAM is read
    BamScanStats stats;
    StrandObservations strand_obs;
    const AccumulatorSet& accumulators() const { return *acc_set_; }

private:
    /// One placement waiting for the offer: offsets into the arenas, resolved to pointers once nothing more grows.
    struct PendingPlacement {
        int32_t ref_id;
        int64_t start, end;
        int32_t align_strand, sj_strand;
        size_t obs0, obs1, hyp0, hyp1;
        int32_t nm, clip;  // the hit's sequence: the second pass's draw prices a multimapper's placements by it
    };
    struct PendingHypothesis {
        size_t imp0, imp1;
        int32_t sj_strand;
    };

    void scan_template(QnameGroup& group);
    void place_template(TemplateJudgment& judgment, int64_t frag_id);
    void count_resolved(const ResolvedFragment& r, bool unique_mapper, bool first_resolved_hit);
    void train_strand(const ResolvedFragment& r);
    void collect_placement(const AssembledFragment& f, const ResolvedFragment& r);
    void offer_placements();
    void add_place_masses(const AssembledFragment& f, const ResolvedFragment& r);

    const FragmentResolver& resolver_;
    const bool include_multimap_;
    const int64_t chunk_size_;
    ResolverScratch scratch_;
    JudgeScratch judge_;
    FragmentChunk chunk_;
    std::unique_ptr<AccumulatorSet> acc_set_;

    // The deposit path allocates no per-fragment memory (capacity is never released). Each resolved hit is collected
    // as a pending placement, its observed and implied introns copied into the arenas because the resolver's result
    // dies with the hit; after the last hit the whole set is offered at once (`AccumulatorSet::offer`).
    std::vector<PendingPlacement> pending_;
    std::vector<PendingHypothesis> pending_hypotheses_;
    std::vector<IntronBlock> arena_observed_, arena_implied_;
    std::vector<accumulator::GapHypothesis> arena_spans_;
    std::vector<accumulator::Placement> placements_;
    accumulator::OfferScratch offer_scratch_;

    // A multimapper's genic place's deposit masses, for its capture weight (`AccumulatorSet::placement_masses`)
    accumulator::DepositScratch mass_scratch_;
    std::vector<IntronBlock> mass_observed_;
    std::vector<accumulator::PlaceMass> masses_;
};

void BamScanner::Worker::scan_template(QnameGroup& group) {
    auto& records = group.records;
    if (records.empty()) return;

    // NH counts a template's pair alignments, and an aligner writes a mate's record in each; a mate reported in fewer
    // records means the aligner truncated its output (counted). A template none of whose primary or secondary records
    // is eligible reports no alignment to compare NH with.
    int32_t n_r1 = 0, n_r2 = 0;
    for (const auto& r : records)
        if (!r.is_supplementary()) ++(r.is_read1() ? n_r1 : n_r2);
    const int32_t reported = std::max(n_r1, n_r2);
    if (reported > 0 && template_nh(records) > reported) stats.n_nh_above_reported++;

    TemplateJudgment judgment = judge_template(records, group.frag_id, group.name_hash, include_multimap_, resolver_,
                                               scratch_, judge_, stats.pairs);
    stats.n_nm_conflicting += judgment.hits.n_nm_conflicting;
    for (int32_t v : judge_.verdicts)
        if (v != HIT_USABLE) stats.hits_rejected[v]++;
    if (judgment.placed()) place_template(judgment, group.frag_id);
}

/// A placed template's usable hits, each resolved into the chunk and collected as a placement, then offered together.
/// Its multiplicity is its usable hits, not the aligner's NH: a hit that is not usable is no place the template could
/// be from.
void BamScanner::Worker::place_template(TemplateJudgment& judgment, int64_t frag_id) {
    std::vector<Hit>& hits = judgment.hits.hits;
    size_t kept = 0;
    for (size_t h = 0; h < hits.size(); ++h)
        if (judge_.verdicts[h] == HIT_USABLE) hits[kept++] = hits[h];
    hits.resize(kept);

    pending_.clear();
    pending_hypotheses_.clear();
    arena_observed_.clear();
    arena_implied_.clear();

    const int32_t num_hits = judgment.n_usable;
    const bool unique_mapper = num_hits == 1;
    stats.n_fragments++;  // one fragment per molecule, not per hit
    // The splice census counts fragments, as the accumulator's offer does: a template once, as its first hit with a
    // candidate (the hit `count_resolved` reads), else as its first hit.
    int32_t census_type = -1, first_type = -1;
    bool any_hit_resolved = false, any_unresolved = false;
    int32_t n_buffered_mm = 0;
    // A multimapper's intergenic hits wait for the last hit: the EM receives them only beside a hit that resolves (a
    // place its gDNA can be at outside every locus); a template none of whose hits resolves stays certified gDNA and
    // reaches no unit.
    std::vector<ResolvedFragment> held_intergenic;

    for (const Hit& hit : hits) {
        AssembledFragment frag;
        ResolvedFragment r;
        resolve_hit(hit, resolver_, scratch_, frag, r);
        r.num_hits = num_hits;
        r.nm = frag.nm;
        r.clip = frag.clip;
        if (first_type < 0) first_type = r.splice_type;
        if (census_type < 0 && !r.t_inds.empty()) census_type = r.splice_type;
        if (unique_mapper && frag.introns.empty()) {
            stats.n_unique_unspliced_mismatches += frag.nm;
            stats.n_unique_unspliced_bases += frag.aligned;
        }
        // Every hit is a placement for calibration, an intergenic one whatever the fragment's other hits: gDNA sits
        // anywhere, and a repeat's gDNA often sits exactly here.
        collect_placement(frag, r);

        // No candidate: an intergenic fragment, unspliced (a spliced hit the scan keeps has a transcript holding its
        // junctions). A unique mapper's is appended, so calibration's geometric categorisation labels it INTERGENIC;
        // the C++ scorer skips a fragment with no candidates, so it does not count as gated either.
        if (r.t_inds.empty()) {
            any_unresolved = true;
            if (unique_mapper) chunk_.append(r, frag_id);
            else held_intergenic.push_back(std::move(r));
            continue;
        }
        any_hit_resolved = true;
        count_resolved(r, unique_mapper, unique_mapper || n_buffered_mm == 0);
        // A multimapper's genic place carries its deposit masses, which its capture weight reads; a unique fragment's
        // would be common to its candidates, and an intergenic place weighs 1.
        if (!unique_mapper) add_place_masses(frag, r);
        chunk_.append(r, frag_id);
        if (!unique_mapper) n_buffered_mm++;
    }
    if (any_hit_resolved) {
        for (const ResolvedFragment& ig : held_intergenic) {
            chunk_.append(ig, frag_id);
            n_buffered_mm++;
        }
    }
    offer_placements();
    if (census_type < 0) census_type = first_type;
    if (census_type >= 0 && static_cast<size_t>(census_type) < stats.splice_census.size())
        stats.splice_census[census_type]++;

    // Intergenic only when no hit resolved, counted once: a multimapper with a resolved hit goes through the EM, which
    // allocates its gDNA, so counting it here would count the molecule twice.
    if (!any_hit_resolved && any_unresolved) stats.n_intergenic++;
    if (n_buffered_mm > 0) {
        stats.n_multimapper_groups++;
        stats.n_multimapper_alignments += n_buffered_mm;
    }
}

/// A resolved hit's exon, splice and strand counts, on a multimapper's first resolved hit only; a unique mapper also
/// trains the strand models.
void BamScanner::Worker::count_resolved(const ResolvedFragment& r, bool unique_mapper, bool first_resolved_hit) {
    if (first_resolved_hit) {
        stats.n_with_exon++;
        if (r.splice_type == SPLICE_SPLICED_ANNOT) stats.n_with_annotated_sj++;
        ++(r.is_same_strand() ? stats.n_same_strand : stats.n_ambig_strand);
    }
    if (unique_mapper) train_strand(r);
}

void BamScanner::Worker::train_strand(const ResolvedFragment& r) {
    if (r.is_strand_qualified()) {
        // One observation per strand-qualified fragment, keyed by its junction (`sj_strand` is POS or NEG here, so
        // `sj_key_*` is set). Python checks `n_strand_trained` against the table's total depth.
        auto& counts = strand_obs.sj_strand_table[SJKey{r.sj_key_ref, r.sj_key_start, r.sj_key_end, r.sj_strand}];
        ++(r.align_strand == r.sj_strand ? counts.n_sense : counts.n_antisense);
        stats.n_strand_trained++;
    } else if (r.splice_type != SPLICE_SPLICED_ANNOT) {
        stats.n_strand_skipped_no_sj++;
    } else if (!r.is_same_strand()) {
        stats.n_strand_skipped_ambig_strand++;
    } else {
        stats.n_strand_skipped_ambiguous++;
    }
    if (r.is_same_strand() && (r.align_strand == STRAND_POS || r.align_strand == STRAND_NEG)) {
        const int32_t t_strand = resolver_.transcript_strand(r.t_inds.front());
        if (t_strand == STRAND_POS || t_strand == STRAND_NEG) strand_obs.exonic.add(r.align_strand, t_strand);
    }
}

// ── the deposit adapter ──────────────────────────────────────────────────────────────────────────────────────────
//
// One hit becomes one `Placement` of the accumulator's `OfferedFragment`; the accumulator owns the deposit rule. The
// two strands are independent: `align_strand` (every read has one) selects the array column, and `sj_strand` (spliced
// reads only) resolves an intron against the annotation. There is no strand gate here: the deposit rejects and counts
// an undefined `align_strand` itself.

void BamScanner::Worker::collect_placement(const AssembledFragment& f, const ResolvedFragment& r) {
    // The introns the CIGAR stated are cut under every hypothesis, so are no part of the hypothesis set: `deposit`
    // unions them with each hypothesis's implied introns. All lie on the hit's reference, as `Accumulator::deposit`
    // assumes when it normalises by coordinate alone.
    const size_t obs0 = arena_observed_.size();
    append_observed_introns(f.introns, arena_observed_);

    // Every explanation of the unsequenced gaps: the resolver enumerates them and the accumulator arbitrates; this
    // re-presents the flat CSR as the span-of-spans `deposit` takes. There is always at least one, the unspliced
    // hypothesis, so `deposit` never receives an empty set.
    const size_t hyp0 = pending_hypotheses_.size();
    for (int32_t h = 0; h < r.n_gap_hypotheses(); ++h) {
        const int32_t i0 = r.gap_intron_offsets[h], i1 = r.gap_intron_offsets[h + 1];
        const size_t imp0 = arena_implied_.size();
        arena_implied_.insert(arena_implied_.end(), r.gap_introns.begin() + i0, r.gap_introns.begin() + i1);
        pending_hypotheses_.push_back({imp0, arena_implied_.size(), r.gap_sj_strand[h]});
    }

    // The molecule's extent: leftmost block start to rightmost block end, the mate gap included, because the gap is
    // part of the molecule and counts toward L. The observed motif strand: `deposit` falls back to a hypothesis's
    // implied strand only when nothing was sequenced, since an observed motif is evidence and an implied one an
    // inference. A usable hit lies on a reference the index, and so the accumulator, holds.
    const std::pair<int32_t, int32_t> extent = fragment_extent(f.exons);
    pending_.push_back({f.exons.front().ref_id, extent.first, extent.second, r.align_strand, r.sj_strand, obs0,
                        arena_observed_.size(), hyp0, pending_hypotheses_.size(), f.nm, f.clip});
}

/// The offer: every placement of the template at once.
void BamScanner::Worker::offer_placements() {
    if (pending_.empty()) return;
    // Offsets to pointers, now that nothing more will grow.
    arena_spans_.clear();
    arena_spans_.reserve(pending_hypotheses_.size());
    for (const auto& h : pending_hypotheses_)
        arena_spans_.push_back({arena_implied_.data() + h.imp0, h.imp1 - h.imp0, h.sj_strand});
    placements_.clear();
    placements_.reserve(pending_.size());
    for (const auto& p : pending_) {
        accumulator::Placement pl;
        pl.ref_id = p.ref_id;
        pl.offered.start = p.start;
        pl.offered.end = p.end;
        pl.offered.observed_introns = arena_observed_.data() + p.obs0;
        pl.offered.n_observed_introns = p.obs1 - p.obs0;
        pl.offered.align_strand = p.align_strand;
        pl.offered.sj_strand = p.sj_strand;
        pl.offered.hypotheses = arena_spans_.data() + p.hyp0;
        pl.offered.n_hypotheses = p.hyp1 - p.hyp0;
        pl.offered.nm = p.nm;
        pl.offered.clip = p.clip;
        placements_.push_back(pl);
    }
    acc_set_->offer(placements_.data(), placements_.size(), offer_scratch_);
}

/// A place's deposit masses, which its capture weight reads: the deposit rule's mass at every object the hit's observed
/// path (the extent `collect_placement` offers, with the CIGAR's introns cut) deposits on, on the flat axes, into the
/// chunk's next row.
void BamScanner::Worker::add_place_masses(const AssembledFragment& f, const ResolvedFragment& r) {
    mass_observed_.clear();
    append_observed_introns(f.introns, mass_observed_);
    const std::pair<int32_t, int32_t> extent = fragment_extent(f.exons);
    const accumulator::OfferedFragment off{extent.first,          extent.second, mass_observed_.data(),
                                           mass_observed_.size(), r.align_strand, r.sj_strand, nullptr, 0};
    acc_set_->placement_masses(f.exons.front().ref_id, off, mass_scratch_, masses_);
    for (const auto& m : masses_) chunk_.add_place_mass(m.kind, m.index, m.mass);
}

// ================================================================
// The scan: a reader thread, worker threads, and the main thread draining chunks to Python
// ================================================================

namespace {

/// Whether a read pair holds every mate its records call mapped, read from every record's flag, filtered records
/// included: a record names its own segment (read 1 or read 2) and, unless its mate is flagged unmapped, expects the
/// other one in the run.
class MateCheck {
public:
    void add(uint16_t flag) {
        const uint8_t segment = (flag & BAM_FREAD1) ? kRead1 : (flag & BAM_FREAD2) ? kRead2 : 0;
        seen_ |= segment;
        if (segment != 0 && (flag & BAM_FPAIRED) && !(flag & BAM_FMUNMAP)) expected_ |= segment ^ (kRead1 | kRead2);
    }
    /// Whether the read pair just read lacks a mate one of its records calls mapped; ready for the next one.
    bool close() {
        const bool lacking = (expected_ & ~seen_) != 0;
        seen_ = expected_ = 0;
        return lacking;
    }

private:
    static constexpr uint8_t kRead1 = 1, kRead2 = 2;
    uint8_t seen_ = 0, expected_ = 0;
};

/// The first exception any of several threads throws.
class FirstError {
public:
    void keep(std::exception_ptr e) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!error_) error_ = e;
    }
    std::exception_ptr get() const { return error_; }

private:
    std::mutex mutex_;
    std::exception_ptr error_;
};

}  // namespace

void BamScanner::read_bam(const std::string& bam_path, int n_decomp_threads, int qname_batch_size,
                          BoundedQueue<QnameBatch>& input) {
    HtsFile fp(hts_open(bam_path.c_str(), "rb"));
    if (!fp) throw std::runtime_error("Failed to open BAM file: " + bam_path);
    if (n_decomp_threads > 0) hts_set_threads(fp.get(), n_decomp_threads);
    BamHeader hdr(sam_hdr_read(fp.get()));
    if (!hdr) throw std::runtime_error("Failed to read BAM header: " + bam_path);
    const RecordParser parser(hdr.get(), *resolver_, sj_tag_mode_);
    BamRecord b = new_bam_record();

    const auto batch_size = static_cast<size_t>(qname_batch_size);
    QnameBatch batch;
    batch.groups.reserve(batch_size);
    // Hand the current batch to the workers; false when the queue was aborted.
    const auto send_batch = [&]() -> bool {
        if (batch.groups.empty()) return true;
        if (!input.push(std::move(batch))) return false;
        batch = QnameBatch{};
        batch.groups.reserve(batch_size);
        return true;
    };

    ReadPairRuns runs;
    MateCheck mates;
    QnameGroup group;
    int64_t frag_id = 0;
    // Close the read pair just read: one with an eligible record is a group, given the next frag_id; one without is
    // excluded. False when the queue was aborted.
    const auto close_pair = [&]() -> bool {
        stats_.n_read_pairs_lacking_mate += mates.close();
        if (group.records.empty()) {
            stats_.pairs.n_excluded++;
            return true;
        }
        group.frag_id = frag_id++;
        batch.groups.push_back(std::move(group));
        group = QnameGroup{};
        return batch.groups.size() < batch_size || send_batch();
    };

    int read_ret = 0;
    while ((read_ret = sam_read1(fp.get(), hdr.get(), b.get())) >= 0) {
        const uint16_t flag = b->core.flag;
        stats_.count_flags(flag);
        const bool pair_open = runs.open();
        if (runs.opens(b.get())) {
            if (pair_open && !close_pair()) return;
            stats_.pairs.n_read_pairs++;
        }
        mates.add(flag);
        if (is_filtered(flag, skip_duplicates_)) {
            stats_.pairs.filtered++;
            continue;
        }
        if (!(flag & BAM_FPAIRED)) {
            throw UnpairedRecordError("record " + std::string(bam_get_qname(b.get())) +
                                      " is not paired (its FLAG lacks 0x1): Rigel reads paired-end libraries only");
        }
        if (group.records.empty()) group.name_hash = runs.name_hash();
        group.records.push_back(parser.parse(b.get(), &stats_.records));
    }
    if (read_ret < -1) throw read_error(bam_path, read_ret, stats_.total);
    if (runs.open() && !close_pair()) return;
    send_batch();
}

void BamScanner::require_partition() const {
    // Both calls are required, because the failure is otherwise invisible: with no junction table every observed
    // intron reads as unannotated, and the tally still looks well-formed.
    if (!acc_set_) {
        throw std::runtime_error(
            "scan: set_regions was not called. The calibration tally is not optional — install the "
            "index's region partition with set_regions, then its sj with set_sj.");
    }
    if (!partition_.has_sj) {
        throw std::runtime_error(
            "scan: set_regions was called but set_sj was not. The sj table is not "
            "optional — without it every observed intron reads as unannotated, so the splice junctions "
            "and the spliced banks are silently empty. Call set_sj with empty arrays if this "
            "annotation genuinely has no sj.");
    }
}

void BamScanner::merge_workers(std::vector<std::unique_ptr<Worker>>& workers) {
    for (auto& w : workers) {
        stats_.merge_from(w->stats);
        strand_obs_.merge_from(w->strand_obs);
        acc_set_->merge_from(w->accumulators());
    }
}

namespace {

/// The main thread's part: hand each finished chunk to Python, holding the GIL only while calling back. An exception
/// aborts both queues, so the reader and the workers stop, and is returned.
template <class InputQueue>
std::exception_ptr hand_chunks_to_python(BoundedQueue<FragmentChunk>& output, InputQueue& input,
                                         nb::callable& chunk_callback) {
    try {
        FragmentChunk chunk;
        nb::gil_scoped_release release;
        while (output.pop(chunk)) {
            nb::gil_scoped_acquire acquire;
            chunk_callback(chunk.finalize());
        }
    } catch (...) {
        input.abort();
        output.abort();
        return std::current_exception();
    }
    return nullptr;
}

}  // namespace

nb::dict BamScanner::scan(const std::string& bam_path, nb::callable chunk_callback, int64_t chunk_size, int n_workers,
                          int n_decomp_threads, int qname_batch_size) {
    if (n_workers < 1 || chunk_size < 1 || qname_batch_size < 1) {
        throw std::invalid_argument("scan: n_workers, chunk_size and qname_batch_size must each be >= 1; got " +
                                    std::to_string(n_workers) + ", " + std::to_string(chunk_size) + " and " +
                                    std::to_string(qname_batch_size));
    }
    require_partition();

    BoundedQueue<QnameBatch> input(static_cast<size_t>(n_workers * INPUT_BATCHES_PER_WORKER));
    BoundedQueue<FragmentChunk> output(static_cast<size_t>(n_workers * OUTPUT_CHUNKS_PER_WORKER));
    std::vector<std::unique_ptr<Worker>> workers;
    for (int i = 0; i < n_workers; i++)
        workers.push_back(std::make_unique<Worker>(*resolver_, partition_, include_multimap_, chunk_size));

    // Workers take batches and push full chunks; the first worker exception is kept and aborts both queues.
    FirstError worker_error;
    std::vector<std::thread> worker_threads;
    for (auto& worker : workers) {
        worker_threads.emplace_back([&input, &output, &worker_error, w = worker.get()] {
            try {
                w->run(input, output);
            } catch (...) {
                worker_error.keep(std::current_exception());
                input.abort();
                output.abort();
            }
        });
    }
    const auto join_workers = [&worker_threads] {
        for (auto& t : worker_threads)
            if (t.joinable()) t.join();
    };

    // The reader reads the BAM, then waits for the workers and merges what they found.
    std::exception_ptr reader_error;
    std::thread reader([&] {
        try {
            read_bam(bam_path, n_decomp_threads, qname_batch_size, input);
            input.close();
            join_workers();
            merge_workers(workers);
            output.close();
        } catch (...) {
            reader_error = std::current_exception();
            input.abort();
            join_workers();
            output.abort();
        }
    });
    const std::exception_ptr main_error = hand_chunks_to_python(output, input, chunk_callback);
    reader.join();

    // The reader's exception first (an unreadable BAM is the root cause of what followed), then a worker's, then the
    // main thread's. An aborted push or pop returns false rather than throwing, so no stored exception echoes another.
    if (reader_error) std::rethrow_exception(reader_error);
    if (worker_error.get()) std::rethrow_exception(worker_error.get());
    if (main_error) std::rethrow_exception(main_error);
    return build_result();
}

}  // namespace rigel
