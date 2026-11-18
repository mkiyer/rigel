// bam_record.h — one BAM record as both passes read it (the scan and the annotated-BAM pass): the htslib handles, the
// junction-strand and integer tags, the CIGAR's blocks and junctions, and where read pairs begin.
#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "resolve/fragment_resolver.h"

namespace rigel {

// ================================================================
// htslib handles, released on every path
// ================================================================

struct HtsFileClose {
    void operator()(htsFile* f) const noexcept { hts_close(f); }
};
struct BamHeaderDestroy {
    void operator()(bam_hdr_t* h) const noexcept { bam_hdr_destroy(h); }
};
struct BamRecordDestroy {
    void operator()(bam1_t* b) const noexcept { bam_destroy1(b); }
};
using HtsFile = std::unique_ptr<htsFile, HtsFileClose>;
using BamHeader = std::unique_ptr<bam_hdr_t, BamHeaderDestroy>;
using BamRecord = std::unique_ptr<bam1_t, BamRecordDestroy>;

/// A fresh record; throws when htslib cannot allocate one.
BamRecord new_bam_record();

/// htslib's sam_read1 returns -1 at the end of the file and less than -1 on an error. A loop that stops on either would
/// quantify a truncated or corrupt BAM, or a CRAM whose reference cannot be found, on the records read before the
/// error, with no sign of it; each read loop keeps the return and throws this instead.
std::runtime_error read_error(const std::string& path, int ret, int64_t n_records);

/// Whether both passes drop a record: unmapped, QC-failed, or a duplicate when duplicates are skipped.
inline bool is_filtered(uint16_t flag, bool skip_duplicates) {
    return (flag & (BAM_FQCFAIL | BAM_FUNMAP)) || ((flag & BAM_FDUP) && skip_duplicates);
}

// ================================================================
// Tags
// ================================================================

enum class SJTagMode : uint8_t {
    NONE = 0,     // No SJ strand tag
    XS_ONLY = 1,  // Check XS tag only
    TS_ONLY = 2,  // Check ts tag only (minimap2)
    XS_TS = 3,    // Check XS first, then ts
    TS_XS = 4,    // Check ts first, then XS
    AUTO = 5,     // Read both on every record; a disagreement is AMBIGUOUS
};

/// The junction-strand tag spec as a mode. Any other spec is refused.
SJTagMode parse_sj_tag_spec(const std::string& spec);

/// What the junction-strand tag said on each spliced record. A tag holding neither '+' nor '-' reads as absent; under
/// AUTO a record whose XS and ts both give a strand and disagree is AMBIGUOUS.
struct SJTagCounts {
    int64_t tagged = 0;       // spliced records given a strand
    int64_t untagged = 0;     // spliced records no tag gave a strand
    int64_t conflicting = 0;  // spliced records whose XS and ts disagree (AUTO only)
    int64_t invalid = 0;      // XS or ts tags holding neither '+' nor '-'

    void merge_from(const SJTagCounts& o) {
        tagged += o.tagged;
        untagged += o.untagged;
        conflicting += o.conflicting;
        invalid += o.invalid;
    }
};

/// An integer tag on each record (NM, NH): absent, or held as anything but an integer of 0..2^31-1 (read as absent).
struct IntTagCounts {
    int64_t missing = 0;
    int64_t invalid = 0;

    void merge_from(const IntTagCounts& o) {
        missing += o.missing;
        invalid += o.invalid;
    }
};

/// What parsing the records found in their tags and junctions. The scan counts it; the annotation pass reads the same
/// records and counts nothing.
struct RecordCounts {
    SJTagCounts sj_tags;   // on each spliced record
    IntTagCounts nm_tag;   // on each record
    IntTagCounts nh_tag;   // on each record
    int64_t n_sj_observed = 0;  // CIGAR junctions, per read, across all alignments

    void merge_from(const RecordCounts& o) {
        sj_tags.merge_from(o.sj_tags);
        nm_tag.merge_from(o.nm_tag);
        nh_tag.merge_from(o.nh_tag);
        n_sj_observed += o.n_sj_observed;
    }
};

// ================================================================
// A parsed record
// ================================================================

/// A CIGAR splice junction: the N's reference span, and the strand the record's junction-strand tag gives it.
struct SJCigarEntry {
    int32_t start;
    int32_t end;
    int32_t strand;
};

/// Every field of one record the passes read, so nothing downstream touches a bam1_t and no raw record crosses a
/// thread.
struct ParsedAlignment {
    int32_t ref_id;    // the index's reference id; -1 for a reference outside the index
    int32_t tid = -1;  // the BAM's reference id
    int32_t ref_start;
    uint16_t flag;
    int32_t nm;  // NM tag (edit distance), 0 if absent
    int32_t nh;  // NH tag, 1 if absent

    std::vector<std::pair<int32_t, int32_t>> exons;  // the CIGAR's aligned blocks, (start, end)
    std::vector<SJCigarEntry> sjs;                   // its junctions

    // The soft clips at the CIGAR's two ends, which the pair validator walks from the read's 5' end
    int32_t clip_lead = 0;
    int32_t clip_trail = 0;

    std::vector<uint32_t> cigar;  // the CIGAR as the alignment's identity (`canonical_cigar`)

    // Whether the record carries SA, naming the other pieces of a chimeric alignment; and, on a location's
    // representative, whether the location is one such piece (`template_hits`)
    bool has_sa = false;
    bool chimeric = false;

    bool is_read1()         const { return flag & BAM_FREAD1; }
    bool is_secondary()     const { return flag & BAM_FSECONDARY; }
    bool is_supplementary() const { return flag & BAM_FSUPPLEMENTARY; }
    bool is_reverse()       const { return flag & BAM_FREVERSE; }
};

/// One BAM's records parsed as both passes read them: on the index's reference ids, with the junction-strand tag read
/// under one mode.
class RecordParser {
public:
    RecordParser(const bam_hdr_t* hdr, const FragmentResolver& resolver, SJTagMode sj_tag_mode);

    /// ``b`` parsed; ``counts`` takes what its tags and junctions held, when given.
    ParsedAlignment parse(const bam1_t* b, RecordCounts* counts) const;

private:
    std::vector<int32_t> tid_to_ref_id_;
    SJTagMode sj_tag_mode_;
};

// ================================================================
// Read pairs
// ================================================================

/// Where read pairs begin. A read pair is a run of records sharing a name, filtered records included, so a name whose
/// every record is filtered is a read pair too, and the annotation pass, which writes a filtered record where it
/// stands, closes a pair where the scan does.
class ReadPairRuns {
public:
    /// Whether ``b`` opens a new read pair; call on every record, in order.
    bool opens(const bam1_t* b) {
        const char* name = bam_get_qname(b);
        if (open_ && name_ == name) return false;
        name_ = name;
        open_ = true;
        return true;
    }
    bool open() const { return open_; }
    /// `std::hash` of the current read pair's name, for the digest both passes fold (`pair_digest`).
    size_t name_hash() const { return std::hash<std::string>{}(name_); }

private:
    std::string name_;
    bool open_ = false;
};

}  // namespace rigel
