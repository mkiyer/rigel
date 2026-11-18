// bam_record.cpp — reading one BAM record: its tags, and its CIGAR's blocks and junctions.
#include "scan/bam_record.h"

#include <limits>

#include "constants.h"

namespace rigel {

BamRecord new_bam_record() {
    BamRecord b(bam_init1());
    if (!b) throw std::runtime_error("Failed to allocate bam1_t");
    return b;
}

std::runtime_error read_error(const std::string& path, int ret, int64_t n_records) {
    return std::runtime_error("Failed to read " + path + " after " + std::to_string(n_records) +
                              " records (htslib error " + std::to_string(ret) +
                              "): a truncated or corrupt BAM, or a CRAM whose reference cannot be found");
}

SJTagMode parse_sj_tag_spec(const std::string& spec) {
    if (spec == "auto")  return SJTagMode::AUTO;
    if (spec == "XS")    return SJTagMode::XS_ONLY;
    if (spec == "ts")    return SJTagMode::TS_ONLY;
    if (spec == "XS,ts") return SJTagMode::XS_TS;
    if (spec == "ts,XS") return SJTagMode::TS_XS;
    if (spec == "none" || spec.empty()) return SJTagMode::NONE;
    throw std::invalid_argument(
        "sj_tag_spec must be 'auto', 'XS', 'ts', 'XS,ts', 'ts,XS' or 'none'; got '" + spec + "'");
}

namespace {

// ================================================================
// The CIGAR
// ================================================================

/// Whether a CIGAR operation other than N advances the reference position.
bool cigar_advances_ref(int op) {
    return op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CEQUAL || op == BAM_CDIFF;
}

/// The soft clips at the CIGAR's first and last operation.
void read_soft_clips(const bam1_t* b, int32_t& lead, int32_t& trail) {
    lead = trail = 0;
    const int32_t n_cigar = b->core.n_cigar;
    if (n_cigar == 0) return;
    const uint32_t* cigar = bam_get_cigar(b);
    if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) lead = bam_cigar_oplen(cigar[0]);
    if (bam_cigar_op(cigar[n_cigar - 1]) == BAM_CSOFT_CLIP) trail = bam_cigar_oplen(cigar[n_cigar - 1]);
}

/// The CIGAR's aligned blocks and junctions. A deletion lies inside a block, since it occupies reference positions;
/// an insertion occupies none.
void parse_cigar(const bam1_t* b, std::vector<std::pair<int32_t, int32_t>>& exons, std::vector<SJCigarEntry>& sjs) {
    exons.clear();
    sjs.clear();
    const int32_t n_cigar = b->core.n_cigar;
    if (n_cigar == 0) return;

    const uint32_t* cigar = bam_get_cigar(b);
    int32_t pos = b->core.pos;  // 0-based genomic position
    int32_t start = pos;        // current block start
    for (int32_t i = 0; i < n_cigar; i++) {
        const int op = bam_cigar_op(cigar[i]);
        const int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CREF_SKIP) {
            if (pos > start) exons.push_back({start, pos});
            sjs.push_back({pos, pos + len, STRAND_NONE});
            start = pos + len;
            pos = start;
        } else if (cigar_advances_ref(op)) {
            pos += len;
        }
    }
    if (pos > start) exons.push_back({start, pos});
}

/// A record's CIGAR as its alignment's identity: M, = and X read as one operation, adjacent runs of it merged (`50=50X`
/// and `100M` are one alignment); every other operation, soft and hard clips included, as written.
void canonical_cigar(const bam1_t* b, std::vector<uint32_t>& out) {
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

// ================================================================
// The tags
// ================================================================

/// A strand from tag ``tag`` (type 'A', or 'Z'/'H' read by its first character), or -1 when the record lacks it or it
/// holds neither '+' nor '-' (counted invalid). minimap2's ``ts`` is alignment-relative, so it is flipped on a reverse
/// read.
int32_t read_strand_tag(const bam1_t* b, const char* tag, bool is_ts, SJTagCounts* counts) {
    uint8_t* aux = bam_aux_get(b, tag);
    if (!aux) return -1;
    char val = 0;
    if (*aux == 'A') {
        val = bam_aux2A(aux);
    } else if (*aux == 'Z' || *aux == 'H') {
        const char* s = bam_aux2Z(aux);
        if (s && s[0]) val = s[0];
    }
    int32_t strand;
    if (val == '+') {
        strand = STRAND_POS;
    } else if (val == '-') {
        strand = STRAND_NEG;
    } else {
        if (counts) ++counts->invalid;
        return -1;
    }
    if (is_ts && (b->core.flag & BAM_FREVERSE)) strand = (strand == STRAND_POS) ? STRAND_NEG : STRAND_POS;
    return strand;
}

/// The junction strand a spliced record's tags give under ``mode``: STRAND_NONE when none does, STRAND_AMBIGUOUS when
/// under AUTO its XS and ts disagree (neither name proves which is right).
int32_t read_sj_strand(const bam1_t* b, SJTagMode mode, SJTagCounts* counts) {
    int32_t result = -1;
    switch (mode) {
        case SJTagMode::NONE:
            break;
        case SJTagMode::XS_ONLY:
            result = read_strand_tag(b, "XS", false, counts);
            break;
        case SJTagMode::TS_ONLY:
            result = read_strand_tag(b, "ts", true, counts);
            break;
        case SJTagMode::XS_TS:
            result = read_strand_tag(b, "XS", false, counts);
            if (result < 0) result = read_strand_tag(b, "ts", true, counts);
            break;
        case SJTagMode::TS_XS:
            result = read_strand_tag(b, "ts", true, counts);
            if (result < 0) result = read_strand_tag(b, "XS", false, counts);
            break;
        case SJTagMode::AUTO: {
            const int32_t xs = read_strand_tag(b, "XS", false, counts);
            const int32_t ts = read_strand_tag(b, "ts", true, counts);
            if (xs >= 0 && ts >= 0 && xs != ts) {
                if (counts) ++counts->conflicting;
                return STRAND_AMBIGUOUS;
            }
            result = (xs >= 0) ? xs : ts;
            break;
        }
    }
    if (result < 0) {
        if (counts) ++counts->untagged;
        return STRAND_NONE;
    }
    if (counts) ++counts->tagged;
    return result;
}

/// A record's integer tag ``tag``: its value, or ``absent`` when the record lacks it or holds it as anything but an
/// integer of 0..2^31-1, each counted in ``counts``.
int32_t read_int_tag(const bam1_t* b, const char* tag, int32_t absent, IntTagCounts* counts) {
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

}  // namespace

// ================================================================
// The parser
// ================================================================

RecordParser::RecordParser(const bam_hdr_t* hdr, const FragmentResolver& resolver, SJTagMode sj_tag_mode)
    : tid_to_ref_id_(static_cast<size_t>(hdr->n_targets)), sj_tag_mode_(sj_tag_mode) {
    for (int32_t i = 0; i < hdr->n_targets; i++) tid_to_ref_id_[i] = resolver.ref_id(hdr->target_name[i]);
}

ParsedAlignment RecordParser::parse(const bam1_t* b, RecordCounts* counts) const {
    ParsedAlignment rec;
    rec.flag = b->core.flag;
    rec.tid = b->core.tid;
    rec.ref_id = (rec.tid >= 0 && rec.tid < static_cast<int32_t>(tid_to_ref_id_.size())) ? tid_to_ref_id_[rec.tid] : -1;
    rec.ref_start = b->core.pos;
    canonical_cigar(b, rec.cigar);
    rec.has_sa = bam_aux_get(b, "SA") != nullptr;

    // A missing NM reads 0, a disclosed approximation (counted, and warned of); a missing NH reads 1
    rec.nm = read_int_tag(b, "NM", 0, counts ? &counts->nm_tag : nullptr);
    rec.nh = read_int_tag(b, "NH", 1, counts ? &counts->nh_tag : nullptr);

    parse_cigar(b, rec.exons, rec.sjs);
    read_soft_clips(b, rec.clip_lead, rec.clip_trail);
    if (!rec.sjs.empty()) {
        const int32_t sj_strand = read_sj_strand(b, sj_tag_mode_, counts ? &counts->sj_tags : nullptr);
        for (auto& sj : rec.sjs) sj.strand = sj_strand;
    }

    if (counts) counts->n_sj_observed += static_cast<int64_t>(rec.sjs.size());
    return rec;
}

}  // namespace rigel
