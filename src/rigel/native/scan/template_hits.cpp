// template_hits.cpp — a template's hits found by geometry, judged, booked, and a usable hit's fragment resolved.
#include "scan/template_hits.h"

#include <algorithm>
#include <string>
#include <tuple>

#include "scan/labels.h"

namespace nb = nanobind;

namespace rigel {

void PairLedger::merge_from(const PairLedger& o) {
    n_read_pairs += o.n_read_pairs;
    n_excluded += o.n_excluded;
    for (int32_t w = 0; w < N_WITHHELD; w++) withheld[w] += o.withheld[w];
    for (int32_t v = 0; v < N_HIT_VERDICTS; v++) unpaired_reason[v] += o.unpaired_reason[v];
    n_excluded_by_option += o.n_excluded_by_option;
    n_placed_unique += o.n_placed_unique;
    n_placed_multimapping += o.n_placed_multimapping;
    n_records_used += o.n_records_used;
    n_records_unused += o.n_records_unused;
    n_records_unplaced += o.n_records_unplaced;
    filtered += o.filtered;
    template_digest += o.template_digest;
}

void PairLedger::export_to(nb::dict& out) const {
    out["n_read_pairs"] = n_read_pairs;
    out["n_excluded"] = n_excluded;
    for (int32_t w = WITHHELD_NONE + 1; w < N_WITHHELD; w++)
        out[("n_withheld_" + std::string(withheld_label(w))).c_str()] = withheld[w];
    for (int32_t v : {HIT_MISSING_MATE, HIT_REFERENCE, HIT_STRAND, HIT_OUTWARD, HIT_LENGTH})
        out[("n_withheld_unpaired_" + std::string(hit_verdict_label(v))).c_str()] = unpaired_reason[v];
    out["n_excluded_by_option"] = n_excluded_by_option;
    out["n_placed_unique"] = n_placed_unique;
    out["n_placed_multimapping"] = n_placed_multimapping;
    out["n_records_used"] = n_records_used;
    out["n_records_unused"] = n_records_unused;
    out["n_records_unplaced"] = n_records_unplaced;
    out["filtered"] = filtered;
    out["template_digest"] = template_digest;
}

namespace {

/// One record's evidence as both passes read it: the index's reference and the BAM's, its leftmost position, its
/// flags, its CIGAR, its junctions with their motif strands, its soft clips and the tags read (NM, NH, SA). Two passes
/// over different alignments disagree here even where every counter agrees.
uint64_t record_digest(const ParsedAlignment& r) {
    size_t h = 0;
    for (const int32_t v : {r.ref_id, r.tid, r.ref_start, static_cast<int32_t>(r.flag), r.nm, r.nh, r.clip_lead,
                            r.clip_trail, static_cast<int32_t>(r.has_sa)})
        hash_combine(h, v);
    for (const uint32_t op : r.cigar) hash_combine(h, static_cast<int32_t>(op));
    for (const SJCigarEntry& sj : r.sjs) {
        hash_combine(h, sj.start);
        hash_combine(h, sj.end);
        hash_combine(h, sj.strand);
    }
    return static_cast<uint64_t>(h);
}

/// One read pair's part of the digest both passes fold their read pairs into: its frag_id, its name (``name_hash``),
/// its records and their evidence, its outcome (with an unpaired pair's reason) and its used records. Summed over the
/// read pairs, and the evidence over the records, so the order in which the workers take them does not matter.
uint64_t pair_digest(int64_t frag_id, size_t name_hash, uint64_t evidence, int64_t n_records, int32_t outcome,
                     int32_t reason, int64_t n_used) {
    size_t h = name_hash;
    hash_combine(h, static_cast<int32_t>(evidence));
    hash_combine(h, static_cast<int32_t>(evidence >> 32));
    hash_combine(h, static_cast<int32_t>(frag_id));
    hash_combine(h, static_cast<int32_t>(frag_id >> 32));
    hash_combine(h, static_cast<int32_t>(n_records));
    hash_combine(h, outcome);
    hash_combine(h, reason);
    hash_combine(h, static_cast<int32_t>(n_used));
    return static_cast<uint64_t>(h);
}

// ================================================================
// A template's hits: every R1 location with every R2 location
// ================================================================

/// Whether two records carry one alignment: one mate, one BAM reference, one strand, one leftmost aligned position and
/// one CIGAR (`canonical_cigar`). STAR repeats a mate's record in every pair alignment that shares it.
bool same_alignment(const ParsedAlignment& x, const ParsedAlignment& y) {
    return x.is_read1() == y.is_read1() && x.tid == y.tid && x.is_reverse() == y.is_reverse() &&
           x.ref_start == y.ref_start && x.cigar == y.cigar;
}

/// Whether a record's read is hard-clipped at its 5' end: the CIGAR's first operation on a forward record, its last on
/// a reverse one. Those bases are not in the record, so the molecule's end is not either.
bool five_prime_hard_clipped(const ParsedAlignment& r) {
    if (r.cigar.empty()) return false;
    return bam_cigar_op(r.is_reverse() ? r.cigar.back() : r.cigar.front()) == BAM_CHARD_CLIP;
}

/// One alignment of one mate: its records, in record order, and the one that represents it.
struct Location {
    std::vector<int32_t> members;
    int32_t rep = -1;
};

/// One mate's records sorted out: its locations, one per alignment (`same_alignment`); its supplementary records; and
/// its records with no aligned base.
struct MateRecords {
    std::vector<Location> locations;
    std::vector<int32_t> supplementary;
    std::vector<int32_t> no_base;
};

std::array<MateRecords, 2> sort_by_mate(const std::vector<ParsedAlignment>& records) {
    std::array<MateRecords, 2> mates;
    for (int32_t k = 0; k < static_cast<int32_t>(records.size()); ++k) {
        const ParsedAlignment& r = records[k];
        MateRecords& mate = mates[r.is_read1() ? 0 : 1];
        if (r.is_supplementary()) {
            mate.supplementary.push_back(k);
        } else if (r.exons.empty()) {
            mate.no_base.push_back(k);
        } else {
            auto it = std::find_if(mate.locations.begin(), mate.locations.end(),
                                   [&](const Location& l) { return same_alignment(records[l.members.front()], r); });
            if (it == mate.locations.end()) mate.locations.push_back({{k}, -1});
            else it->members.push_back(k);
        }
    }
    return mates;
}

/// Each location's representative, its primary record or else its first, which the other members alias. Records of
/// one alignment that disagree on NM are an input defect: counted, and the larger NM kept. Returns the count.
int32_t choose_representatives(std::vector<ParsedAlignment>& records, MateRecords& mate,
                               std::vector<int32_t>& alias_of) {
    int32_t n_nm_conflicting = 0;
    for (Location& l : mate.locations) {
        l.rep = l.members.front();
        for (int32_t k : l.members) {
            if (!records[k].is_secondary()) {
                l.rep = k;
                break;
            }
        }
        for (int32_t k : l.members)
            if (k != l.rep) alias_of[k] = l.rep;
        int32_t lo = records[l.rep].nm, hi = lo;
        for (int32_t k : l.members) {
            lo = std::min(lo, records[k].nm);
            hi = std::max(hi, records[k].nm);
        }
        if (lo != hi) {
            ++n_nm_conflicting;
            records[l.rep].nm = hi;
        }
    }
    return n_nm_conflicting;
}

/// A location is chimeric when a record of it carries SA or is hard-clipped at its 5' end, or when its mate has a
/// supplementary record and none of the mate's records carries SA, since the piece's alignment is then unknown. A
/// supplementary record, or one with no aligned base, is no location: it is written with the mate's primary location,
/// or its first.
void mark_chimeric_and_strays(std::vector<ParsedAlignment>& records, const MateRecords& mate,
                              std::vector<int32_t>& alias_of) {
    bool sa = false;
    for (int32_t k : mate.supplementary) sa |= records[k].has_sa;
    for (const Location& l : mate.locations)
        for (int32_t k : l.members) sa |= records[k].has_sa;
    const bool unnamed_chain = !mate.supplementary.empty() && !sa;
    for (const Location& l : mate.locations) {
        ParsedAlignment& rep = records[l.rep];
        rep.chimeric = unnamed_chain || five_prime_hard_clipped(rep);
        for (int32_t k : l.members) rep.chimeric |= records[k].has_sa;
    }
    if (mate.locations.empty()) return;
    auto primary = std::find_if(mate.locations.begin(), mate.locations.end(),
                                [&](const Location& l) { return !records[l.rep].is_secondary(); });
    const int32_t owner = (primary == mate.locations.end() ? mate.locations.front() : *primary).rep;
    for (int32_t k : mate.supplementary) alias_of[k] = owner;
    for (int32_t k : mate.no_base) alias_of[k] = owner;
}

/// A template's candidate hits, found by geometry alone (the aligner's own pairing, STAR's HI, is not read). Each R1
/// location pairs with each R2 location on the same BAM reference and opposite strand. With one mate unlocated, each
/// location of the other is a one-mate hit (`HIT_MISSING_MATE`); when no pair qualifies, the first location of each is
/// offered, so the verdict counts the template's reason once.
TemplateHits template_hits(std::vector<ParsedAlignment>& records) {
    TemplateHits out;
    out.alias_of.assign(records.size(), -1);
    std::array<MateRecords, 2> mates = sort_by_mate(records);
    for (MateRecords& mate : mates) out.n_nm_conflicting += choose_representatives(records, mate, out.alias_of);
    for (const MateRecords& mate : mates) mark_chimeric_and_strays(records, mate, out.alias_of);

    const std::vector<Location>& r1 = mates[0].locations;
    const std::vector<Location>& r2 = mates[1].locations;
    if (r1.empty() || r2.empty()) {
        for (const Location& l : r1) out.hits.push_back({&records[l.rep], nullptr});
        for (const Location& l : r2) out.hits.push_back({nullptr, &records[l.rep]});
        return out;
    }
    for (const Location& l1 : r1) {
        ParsedAlignment& a = records[l1.rep];
        for (const Location& l2 : r2) {
            ParsedAlignment& b = records[l2.rep];
            if (a.tid == b.tid && a.is_reverse() != b.is_reverse()) out.hits.push_back({&a, &b});
        }
    }
    if (out.hits.empty()) out.hits.push_back({&records[r1.front().rep], &records[r2.front().rep]});
    return out;
}

// ================================================================
// The hit verdict: two mates that can be one molecule, explained by one annotated transcript when spliced
// ================================================================

/// Every junction a hit observes: the kept junctions of both mates' records, sorted and unique. The pair validator and
/// the fragment the resolver reads take this one set.
void observed_junctions(const Hit& hit, std::vector<IntronBlock>& out) {
    out.clear();
    for (const ParsedAlignment* rec : {hit.r1, hit.r2})
        if (rec)
            for (const auto& sj : rec->sjs) out.push_back({rec->ref_id, sj.start, sj.end, sj.strand});
    sort_unique_junctions(out);
}

/// `HIT_USABLE` when a hit's two mates can be one molecule (`FragmentResolver::validate_pair`), else why not: one mate
/// gives no fragment length (`HIT_MISSING_MATE`), a reference outside the index has no partition (`HIT_OFF_INDEX`), and
/// a mate aligned in pieces is chimeric, which Rigel does not quantify (`HIT_CHIMERIC`). References are told apart by
/// the BAM's ids, since every reference outside the index shares index id -1. ``f`` and ``r`` are set to the forward
/// and reverse mates when the rule judged them.
int32_t pair_verdict(const Hit& hit, const std::vector<IntronBlock>& junctions, const FragmentResolver& resolver,
                     ResolverScratch& scratch, const ParsedAlignment*& f, const ParsedAlignment*& r) {
    f = r = nullptr;
    if (!hit.r1 || !hit.r2) return HIT_MISSING_MATE;
    const ParsedAlignment* m1 = hit.r1;
    const ParsedAlignment* m2 = hit.r2;
    if (m1->ref_id < 0 || m2->ref_id < 0) return m1->tid != m2->tid ? HIT_REFERENCE : HIT_OFF_INDEX;
    if (m1->chimeric || m2->chimeric) return HIT_CHIMERIC;
    if (m1->ref_id != m2->ref_id) return HIT_REFERENCE;
    if (m1->is_reverse() == m2->is_reverse()) return HIT_STRAND;
    f = m1->is_reverse() ? m2 : m1;
    r = m1->is_reverse() ? m1 : m2;
    const PairMate pf{&f->exons, f->clip_lead, f->clip_trail};
    const PairMate pr{&r->exons, r->clip_lead, r->clip_trail};
    return resolver.validate_pair(f->ref_id, pf, pr, junctions, scratch);
}

/// Why a hit is not usable, or `HIT_USABLE`: its mates must be able to be one molecule (`pair_verdict`), and a spliced
/// hit needs one annotated transcript holding every junction it observes and overlapping both mates
/// (`FragmentResolver::pair_holders`), the annotated-path rule. The molecule is judged first, so a hit the rule
/// withholds always has a fragment length.
int32_t hit_verdict(const Hit& hit, const FragmentResolver& resolver, ResolverScratch& scratch,
                    std::vector<IntronBlock>& junctions) {
    observed_junctions(hit, junctions);
    const ParsedAlignment* f = nullptr;
    const ParsedAlignment* r = nullptr;
    const int32_t pair = pair_verdict(hit, junctions, resolver, scratch, f, r);
    if (pair != HIT_USABLE || junctions.empty()) return pair;
    return resolver.pair_holders(junctions, f->exons.front().first, f->exons.back().second, r->exons.front().first,
                                 r->exons.back().second, scratch, scratch.pair_candidates);
}

/// Why a template with no usable hit is withheld, by precedence: off the index, else chimeric, else an unannotated
/// splice, else unpaired (as is a template with no hit); WITHHELD_NONE when one is usable.
int32_t withheld_outcome(const std::vector<int32_t>& verdicts) {
    bool usable = false, off_index = false, chimeric = false, unannotated_splice = false;
    for (int32_t v : verdicts) {
        usable |= v == HIT_USABLE;
        off_index |= v == HIT_OFF_INDEX;
        chimeric |= v == HIT_CHIMERIC;
        unannotated_splice |= is_unannotated_splice(v);
    }
    if (usable) return WITHHELD_NONE;
    if (off_index) return WITHHELD_OFF_INDEX;
    if (chimeric) return WITHHELD_CHIMERIC;
    return unannotated_splice ? WITHHELD_UNANNOTATED_SPLICE : WITHHELD_UNPAIRED;
}

/// Why an unpaired template's mates are not one molecule: the furthest check of the pair rule any of its hits reached,
/// in the order the rule asks them, `HIT_MISSING_MATE` < `HIT_REFERENCE` < `HIT_STRAND` < `HIT_OUTWARD` < `HIT_LENGTH`.
/// A template with no hit at all has no mate.
int32_t unpaired_reason(const std::vector<int32_t>& verdicts) {
    static_assert(HIT_REFERENCE < HIT_STRAND && HIT_STRAND < HIT_OUTWARD && HIT_OUTWARD < HIT_LENGTH,
                  "the pair rule's checks are numbered in the order it asks them");
    int32_t reason = HIT_MISSING_MATE;
    for (int32_t v : verdicts)
        if (v != HIT_MISSING_MATE && (reason == HIT_MISSING_MATE || v > reason)) reason = v;
    return reason;
}

/// Which records of a template are used, ``used[k]`` 1: those whose alignment is a location of a usable hit, an alias
/// of the location included. A supplementary record, or one with no aligned base, joins no hit, so it is never used.
void records_used(const std::vector<ParsedAlignment>& records, const TemplateHits& th,
                  const std::vector<int32_t>& verdicts, std::vector<char>& used) {
    used.assign(records.size(), 0);
    for (size_t h = 0; h < th.hits.size(); ++h) {
        if (verdicts[h] != HIT_USABLE) continue;
        for (const ParsedAlignment* r : {th.hits[h].r1, th.hits[h].r2})
            if (r) used[r - records.data()] = 1;
    }
    for (size_t k = 0; k < records.size(); ++k)
        if (th.alias_of[k] >= 0)
            used[k] = records[k].is_supplementary() || records[k].exons.empty() ? 0 : used[th.alias_of[k]];
}

/// The read bases a record aligns: its M operations (`canonical_cigar` merges = and X into them) and its insertions.
int32_t aligned_bases(const std::vector<uint32_t>& cigar) {
    int32_t n = 0;
    for (const uint32_t c : cigar) {
        const int op = bam_cigar_op(c);
        if (op == BAM_CMATCH || op == BAM_CINS) n += static_cast<int32_t>(bam_cigar_oplen(c));
    }
    return n;
}

AssembledFragment build_fragment(const Hit& hit) {
    AssembledFragment frag;
    std::vector<ExonBlock> blocks;
    for (const bool r2 : {false, true}) {
        const ParsedAlignment* rec = r2 ? hit.r2 : hit.r1;
        if (!rec) continue;
        const int32_t strand = (rec->is_reverse() != r2) ? STRAND_NEG : STRAND_POS;
        for (const auto& [start, end] : rec->exons) blocks.push_back({rec->ref_id, start, end, strand});
        frag.nm += rec->nm;
        frag.clip += rec->clip_lead + rec->clip_trail;
        frag.aligned += aligned_bases(rec->cigar);
    }
    // Overlapping and abutting blocks merge within each reference and strand
    std::sort(blocks.begin(), blocks.end(), [](const ExonBlock& a, const ExonBlock& b) {
        return std::tie(a.ref_id, a.strand, a.start, a.end) < std::tie(b.ref_id, b.strand, b.start, b.end);
    });
    for (const ExonBlock& b : blocks) {
        ExonBlock* last = frag.exons.empty() ? nullptr : &frag.exons.back();
        if (last && last->ref_id == b.ref_id && last->strand == b.strand && b.start <= last->end)
            last->end = std::max(last->end, b.end);
        else
            frag.exons.push_back(b);
    }
    std::sort(frag.exons.begin(), frag.exons.end(), exon_block_less);
    observed_junctions(hit, frag.introns);
    return frag;
}

}  // namespace

// ================================================================
// The judgment
// ================================================================

TemplateJudgment judge_template(std::vector<ParsedAlignment>& records, int64_t frag_id, size_t name_hash,
                                bool include_multimap, const FragmentResolver& resolver, ResolverScratch& scratch,
                                JudgeScratch& js, PairLedger& ledger) {
    TemplateJudgment j;
    uint64_t evidence = 0;
    for (const ParsedAlignment& r : records) evidence += record_digest(r);
    j.hits = template_hits(records);
    js.verdicts.clear();
    for (const Hit& hit : j.hits.hits) js.verdicts.push_back(hit_verdict(hit, resolver, scratch, js.junctions));

    const int64_t n_records = static_cast<int64_t>(records.size());
    const int32_t withheld = withheld_outcome(js.verdicts);
    if (withheld != WITHHELD_NONE) {
        const int32_t reason = withheld == WITHHELD_UNPAIRED ? unpaired_reason(js.verdicts) : 0;
        ledger.withheld[withheld]++;
        if (withheld == WITHHELD_UNPAIRED) ledger.unpaired_reason[reason]++;
        ledger.n_records_unplaced += n_records;
        ledger.template_digest += pair_digest(frag_id, name_hash, evidence, n_records, withheld, reason, 0);
        j.outcome = withheld;
        return j;
    }
    // The multimap option excludes a template with more than one usable hit, its multiplicity
    j.n_usable = static_cast<int32_t>(std::count(js.verdicts.begin(), js.verdicts.end(), HIT_USABLE));
    if (!include_multimap && j.n_usable > 1) {
        ledger.n_excluded_by_option++;
        ledger.n_records_unplaced += n_records;
        ledger.template_digest += pair_digest(frag_id, name_hash, evidence, n_records, PAIR_EXCLUDED_BY_OPTION, 0, 0);
        j.outcome = PAIR_EXCLUDED_BY_OPTION;
        return j;
    }
    records_used(records, j.hits, js.verdicts, js.used);
    const int64_t n_used = std::count(js.used.begin(), js.used.end(), 1);
    j.outcome = j.n_usable == 1 ? PAIR_PLACED_UNIQUE : PAIR_PLACED_MULTIMAPPING;
    ++(j.n_usable == 1 ? ledger.n_placed_unique : ledger.n_placed_multimapping);
    ledger.n_records_used += n_used;
    ledger.n_records_unused += n_records - n_used;
    ledger.template_digest += pair_digest(frag_id, name_hash, evidence, n_records, j.outcome, 0, n_used);
    return j;
}

void resolve_hit(const Hit& hit, const FragmentResolver& resolver, ResolverScratch& scratch, AssembledFragment& frag,
                 ResolvedFragment& resolved) {
    frag = build_fragment(hit);
    resolver.resolve(frag.exons, frag.introns, resolved, scratch);
}

int32_t template_nh(const std::vector<ParsedAlignment>& records) {
    for (const auto& r : records)
        if (r.nh > 1) return r.nh;
    return 1;
}

}  // namespace rigel
