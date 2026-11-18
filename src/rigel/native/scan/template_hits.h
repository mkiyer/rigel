// template_hits.h — one template (a read pair's records) judged as both passes judge it: its candidate hits, found by
// geometry; each hit's verdict; the template's outcome, booked on the ledger both passes keep; and a usable hit's
// fragment, resolved.
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include <nanobind/nanobind.h>

#include "constants.h"
#include "fragment_blocks.h"
#include "resolve/fragment_resolver.h"
#include "scan/bam_record.h"

namespace rigel {

/// One candidate hit: an R1 location and an R2 location, either absent for a one-mate hit.
struct Hit {
    ParsedAlignment* r1 = nullptr;
    ParsedAlignment* r2 = nullptr;
};

/// A template's candidate hits; and, per record, the record it is written with in the annotated BAM when it carries no
/// hit of its own, an alias of an alignment another record represents (-1 otherwise).
struct TemplateHits {
    std::vector<Hit> hits;
    std::vector<int32_t> alias_of;
    int32_t n_nm_conflicting = 0;  // locations whose records disagree on NM
};

/// The read pairs and their records as a pass judged them, and the digest of the read pairs (`pair_digest`). The scan
/// exports it, and the annotation pass is held to it count by count.
struct PairLedger {
    int64_t n_read_pairs = 0;
    int64_t n_excluded = 0;  // every record filtered
    std::array<int64_t, N_WITHHELD> withheld{};            // by `WITHHELD_*` outcome
    std::array<int64_t, N_HIT_VERDICTS> unpaired_reason{};  // the unpaired ones by `unpaired_reason`
    int64_t n_excluded_by_option = 0;
    int64_t n_placed_unique = 0;
    int64_t n_placed_multimapping = 0;
    // One status per eligible record: used (its alignment is a location of a usable hit), unused (in a placed read
    // pair, in no usable hit), unplaced (in a read pair with any other outcome)
    int64_t n_records_used = 0;
    int64_t n_records_unused = 0;
    int64_t n_records_unplaced = 0;
    int64_t filtered = 0;  // records dropped as unmapped, QC-failed or a skipped duplicate
    uint64_t template_digest = 0;

    void merge_from(const PairLedger& o);
    /// Every count under its stats-export key.
    void export_to(nanobind::dict& out) const;
};

/// A template's judgment: its hits, its outcome (a `WITHHELD_*` reason, `PAIR_EXCLUDED_BY_OPTION` or a
/// `PAIR_PLACED_*`) and its usable hits, its multiplicity.
struct TemplateJudgment {
    TemplateHits hits;
    int32_t outcome = WITHHELD_UNPAIRED;
    int32_t n_usable = 0;

    bool placed() const { return outcome == PAIR_PLACED_UNIQUE || outcome == PAIR_PLACED_MULTIMAPPING; }
};

/// The working memory of judging one template, beside the resolver's: a hit's junctions while it is judged, each hit's
/// verdict, and which records are used.
struct JudgeScratch {
    std::vector<IntronBlock> junctions;
    std::vector<int32_t> verdicts;
    std::vector<char> used;
};

/// Judge one template and book it on ``ledger``. Its hits are found by geometry (`template_hits`) and each is judged
/// (`hit_verdict`, into ``js.verdicts``): a hit is usable when its two mates can be one molecule and one annotated
/// transcript holds every junction it observes. A template with no usable hit is withheld, by precedence off the
/// index, else chimeric, else an unannotated splice, else unpaired; the multimap option excludes one with more than
/// one; any other is placed, its records' statuses into ``js.used`` (`records_used`). The ledger takes its outcome,
/// its records and its part of the digest. Both passes judge here, so the annotated BAM says what the scan counted.
TemplateJudgment judge_template(std::vector<ParsedAlignment>& records, int64_t frag_id, size_t name_hash,
                                bool include_multimap, const FragmentResolver& resolver, ResolverScratch& scratch,
                                JudgeScratch& js, PairLedger& ledger);

/// One hit's mates as one fragment: their aligned blocks merged per reference and strand, R2's strand flipped, so a
/// block carries the fragment's strand; the hit's junctions (sorted and unique); and its sequence summed: NM, the
/// soft-clipped bases and the aligned ones.
struct AssembledFragment {
    std::vector<ExonBlock> exons;
    std::vector<IntronBlock> introns;
    int32_t nm = 0;
    int32_t clip = 0;
    int32_t aligned = 0;  // M, =, X and I
};

/// One usable hit's fragment and its resolution, as both passes resolve it: the fragment assembled and resolved
/// against the annotation.
void resolve_hit(const Hit& hit, const FragmentResolver& resolver, ResolverScratch& scratch, AssembledFragment& frag,
                 ResolvedFragment& resolved);

/// A template's NH: the first record's above 1, else 1 (a missing NH reads 1).
int32_t template_nh(const std::vector<ParsedAlignment>& records);

}  // namespace rigel
