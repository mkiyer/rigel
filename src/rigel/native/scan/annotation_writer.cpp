// annotation_writer.cpp — the annotated-BAM pass: read pair by read pair, judged as the scan judged it, each record
// written once with its fragment's tags.
#include "scan/annotation_writer.h"

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "constants.h"
#include "scan/labels.h"
#include "scan/template_hits.h"

namespace nb = nanobind;

namespace rigel {

BamAnnotationWriter::BamAnnotationWriter(const FragmentResolver& resolver, const std::string& sj_tag_spec,
                                         bool skip_duplicates, bool include_multimap, int32_t max_fragment_length,
                                         int32_t intron_overhang_tolerance)
    : resolver_(&resolver),
      sj_tag_mode_(parse_sj_tag_spec(sj_tag_spec)),
      skip_duplicates_(skip_duplicates),
      include_multimap_(include_multimap),
      policy_(max_fragment_length, intron_overhang_tolerance) {}

namespace {

/// A fragment class code as the ZC tag names it.
const char* frag_class_label(int code) {
    switch (code) {
        case FRAG_UNAMBIG:           return "unambig";
        case FRAG_AMBIG_SAME_STRAND: return "ambig_same_strand";
        case FRAG_AMBIG_OPP_STRAND:  return "ambig_opp_strand";
        case FRAG_MULTIMAPPER:       return "multimapper";
        default:                     return "unknown";
    }
}

// Rigel owns the tags it stamps: an aligner's tag of the same name (HISAT2 writes ZS:i, its second-best score) is
// replaced, never kept (htslib will not update a tag in place across types), and a failed write throws.
template <class Update>
void replace_tag(bam1_t* r, const char* tag, Update update) {
    if (uint8_t* old = bam_aux_get(r, tag)) bam_aux_del(r, old);
    if (update() < 0) throw std::runtime_error(std::string("Failed to write the ") + tag + " tag");
}
void set_tag_int(bam1_t* r, const char* tag, int64_t v) {
    replace_tag(r, tag, [&] { return bam_aux_update_int(r, tag, v); });
}
void set_tag_float(bam1_t* r, const char* tag, float v) {
    replace_tag(r, tag, [&] { return bam_aux_update_float(r, tag, v); });
}
void set_tag_str(bam1_t* r, const char* tag, const char* v) {
    replace_tag(r, tag, [&] { return bam_aux_update_str(r, tag, static_cast<int>(strlen(v) + 1), v); });
}

/// The tags every record of one fragment carries: its winning transcript and gene (ZT, ZG, ZR by name; ZI, ZJ by
/// index, -1 unassigned), its outcome flags (ZF), posterior (ZW), class (ZC), candidates (ZN) and locus (ZL, -1 none).
/// The default is the intergenic one: gDNA, ZW 1.
struct FragmentTags {
    const char* transcript = ".";
    const char* gene = ".";
    const char* gene_name = ".";
    int32_t transcript_index = -1;
    int32_t gene_index = -1;
    int32_t flags = AF_GDNA_INTERGENIC;
    float posterior = 1.0f;
    const char* frag_class = ".";
    int32_t n_candidates = 0;
    int32_t locus = -1;
};

/// One write() call: the files, the read-pair walk, the current read pair's records, and what the pass has counted.
class AnnotationPass {
public:
    AnnotationPass(const std::string& bam_path, const std::string& output_path,
                   const BamAnnotationWriter::Annotations& annotations, const std::vector<std::string>& t_ids,
                   const std::vector<std::string>& g_ids, const std::vector<std::string>& g_names,
                   const FragmentResolver& resolver, ResolverPolicy policy, SJTagMode sj_tag_mode,
                   bool skip_duplicates, bool include_multimap);

    nb::dict run();

private:
    void close_pair();
    void write_template();
    void write_unplaced(int32_t flags, const char* frag_class);
    void write_placed(const TemplateJudgment& judgment);
    FragmentTags tags_of(int64_t frag_id, bool& annotated) const;
    void write_records(const std::vector<size_t>& indices, const FragmentTags& tags, int32_t primary,
                       const char* splice_label);

    const std::string& bam_path_;
    const std::string& output_path_;
    const BamAnnotationWriter::Annotations& ann_;
    const std::vector<std::string>& t_ids_;
    const std::vector<std::string>& g_ids_;
    const std::vector<std::string>& g_names_;
    const FragmentResolver& resolver_;
    const bool skip_duplicates_;
    const bool include_multimap_;
    std::unordered_map<int64_t, int64_t> row_of_;  // frag_id -> annotation row

    HtsFile in_;
    BamHeader hdr_;
    HtsFile out_;
    RecordParser parser_;

    // The current read pair's eligible records, raw for writing and parsed for judging
    std::vector<BamRecord> raw_;
    std::vector<ParsedAlignment> records_;
    size_t name_hash_ = 0;
    int64_t frag_id_ = 0;

    ResolverScratch scratch_;
    JudgeScratch judge_;
    PairLedger ledger_;
    int64_t n_records_written_ = 0;
};

HtsFile open_input(const std::string& path) {
    HtsFile fp(hts_open(path.c_str(), "rb"));
    if (!fp) throw std::runtime_error("Failed to open BAM: " + path);
    return fp;
}

BamHeader read_header(htsFile* fp, const std::string& path) {
    BamHeader hdr(sam_hdr_read(fp));
    if (!hdr) throw std::runtime_error("Failed to read header: " + path);
    return hdr;
}

HtsFile open_output(const std::string& path, const bam_hdr_t* hdr) {
    HtsFile out(hts_open(path.c_str(), "wb"));
    if (!out) throw std::runtime_error("Failed to open output: " + path);
    if (sam_hdr_write(out.get(), hdr) < 0) throw std::runtime_error("Failed to write header: " + path);
    return out;
}

AnnotationPass::AnnotationPass(const std::string& bam_path, const std::string& output_path,
                               const BamAnnotationWriter::Annotations& annotations,
                               const std::vector<std::string>& t_ids, const std::vector<std::string>& g_ids,
                               const std::vector<std::string>& g_names, const FragmentResolver& resolver,
                               ResolverPolicy policy, SJTagMode sj_tag_mode, bool skip_duplicates,
                               bool include_multimap)
    : bam_path_(bam_path),
      output_path_(output_path),
      ann_(annotations),
      t_ids_(t_ids),
      g_ids_(g_ids),
      g_names_(g_names),
      resolver_(resolver),
      skip_duplicates_(skip_duplicates),
      include_multimap_(include_multimap),
      in_(open_input(bam_path)),
      hdr_(read_header(in_.get(), bam_path)),
      out_(open_output(output_path, hdr_.get())),
      parser_(hdr_.get(), resolver, sj_tag_mode),
      scratch_(resolver.n_transcripts(), policy) {
    row_of_.reserve(static_cast<size_t>(ann_.size));
    for (int64_t i = 0; i < ann_.size; i++) row_of_[ann_.frag_ids.data()[i]] = i;
}

nb::dict AnnotationPass::run() {
    BamRecord b = new_bam_record();
    ReadPairRuns runs;
    int read_ret = 0;
    while ((read_ret = sam_read1(in_.get(), hdr_.get(), b.get())) >= 0) {
        const bool pair_open = runs.open();
        if (runs.opens(b.get())) {
            if (pair_open) close_pair();
            ledger_.n_read_pairs++;
        }
        // A filtered record is written through as-is, without annotation tags, so the output keeps every input record.
        // The scan skips it too, so it gets no frag_id. (An unpaired record the scan did not filter stops the scan.)
        if (is_filtered(b->core.flag, skip_duplicates_)) {
            if (sam_write1(out_.get(), hdr_.get(), b.get()) < 0)
                throw std::runtime_error("Failed to write filtered BAM record");
            ledger_.filtered++;
            n_records_written_++;
            continue;
        }
        if (records_.empty()) name_hash_ = runs.name_hash();
        raw_.emplace_back(bam_dup1(b.get()));
        // The scan's own parser, so the annotated BAM sees the evidence the scan saw
        records_.push_back(parser_.parse(b.get(), nullptr));
    }
    if (read_ret < -1) throw read_error(bam_path_, read_ret, n_records_written_);
    if (runs.open()) close_pair();
    if (hts_close(out_.release()) < 0) throw std::runtime_error("Failed to close output: " + output_path_);

    nb::dict summary;
    ledger_.export_to(summary);
    summary["n_records_written"] = n_records_written_;
    return summary;
}

/// Close the read pair just read: one with an eligible record is written as a template, given the next frag_id as the
/// scan gives it; one without is excluded.
void AnnotationPass::close_pair() {
    if (records_.empty()) {
        ledger_.n_excluded++;
        return;
    }
    write_template();
    frag_id_++;
    raw_.clear();
    records_.clear();
}

/// The scan's judgment, replayed: a withheld template is written unresolved, its class the outcome's label; one the
/// multimap option excludes is written as dropped; a placed one with its annotation.
void AnnotationPass::write_template() {
    const TemplateJudgment judgment = judge_template(records_, frag_id_, name_hash_, include_multimap_, resolver_,
                                                     scratch_, judge_, ledger_);
    if (judgment.placed()) write_placed(judgment);
    else if (judgment.outcome == PAIR_EXCLUDED_BY_OPTION) write_unplaced(AF_MULTIMAPPER_DROP, ".");
    else write_unplaced(AF_UNRESOLVED, withheld_label(judgment.outcome));
}

void AnnotationPass::write_unplaced(int32_t flags, const char* frag_class) {
    FragmentTags tags;
    tags.flags = flags;
    tags.posterior = 0.0f;
    tags.frag_class = frag_class;
    std::vector<size_t> all(records_.size());
    for (size_t i = 0; i < all.size(); ++i) all[i] = i;
    write_records(all, tags, 0, "unknown");
}

/// A placed read pair: a used record is written with the first usable hit holding its alignment, flagged as the
/// winning placement or not (ZH), with that hit's splice type; every other record is an unused alignment. All carry the
/// fragment's annotation row or, with none, the intergenic default. Each input record is written exactly once; an alias
/// of an alignment is written with the record that represents it (`TemplateHits::alias_of`).
void AnnotationPass::write_placed(const TemplateJudgment& judgment) {
    const std::vector<int32_t>& verdicts = judge_.verdicts;
    const std::vector<char>& used = judge_.used;
    const auto first_usable =
        static_cast<int>(std::find(verdicts.begin(), verdicts.end(), HIT_USABLE) - verdicts.begin());
    bool annotated = false;
    const FragmentTags tags = tags_of(frag_id_, annotated);

    std::vector<bool> written(records_.size(), false);
    std::vector<std::vector<size_t>> aliases(records_.size());
    for (size_t i = 0; i < judgment.hits.alias_of.size(); ++i)
        if (judgment.hits.alias_of[i] >= 0) aliases[judgment.hits.alias_of[i]].push_back(i);
    const auto take = [&](size_t idx, std::vector<size_t>& into) {
        if (written[idx]) return;
        written[idx] = true;
        into.push_back(idx);
    };

    for (int h = 0; h < static_cast<int>(judgment.hits.hits.size()); ++h) {
        if (verdicts[h] != HIT_USABLE) continue;
        const Hit& hit = judgment.hits.hits[h];
        std::vector<size_t> hit_records;
        for (const ParsedAlignment* r : {hit.r1, hit.r2}) {
            if (!r) continue;
            const size_t idx = static_cast<size_t>(r - records_.data());
            take(idx, hit_records);
            for (size_t a : aliases[idx])
                if (used[a]) take(a, hit_records);
        }
        if (hit_records.empty()) continue;

        AssembledFragment frag;
        ResolvedFragment resolved;
        resolve_hit(hit, resolver_, scratch_, frag, resolved);
        bool primary = judgment.n_usable == 1;
        if (!primary && annotated && tags.transcript_index >= 0)
            primary = std::find(resolved.t_inds.begin(), resolved.t_inds.end(), tags.transcript_index) !=
                      resolved.t_inds.end();
        else if (!primary && !annotated)
            primary = h == first_usable;
        write_records(hit_records, tags, primary ? 1 : 0, splice_type_label(resolved.splice_type));
    }

    std::vector<size_t> unused;
    for (size_t i = 0; i < records_.size(); ++i) take(i, unused);
    if (unused.empty()) return;
    FragmentTags unused_tags = tags;
    unused_tags.frag_class = "unused_alignment";
    write_records(unused, unused_tags, 0, "unknown");
}

FragmentTags AnnotationPass::tags_of(int64_t frag_id, bool& annotated) const {
    FragmentTags tags;
    const auto it = row_of_.find(frag_id);
    annotated = it != row_of_.end();
    if (!annotated) return tags;
    const int64_t row = it->second;
    tags.transcript_index = ann_.best_tid.data()[row];
    tags.gene_index = ann_.best_gid.data()[row];
    tags.flags = ann_.tx_flags.data()[row];
    tags.posterior = ann_.posterior.data()[row];
    tags.frag_class = frag_class_label(ann_.frag_class.data()[row]);
    tags.n_candidates = ann_.n_candidates.data()[row];
    tags.locus = ann_.locus_id.data()[row];
    const auto in = [](int32_t i, const std::vector<std::string>& names) {
        return i >= 0 && i < static_cast<int32_t>(names.size());
    };
    if (in(tags.transcript_index, t_ids_)) tags.transcript = t_ids_[tags.transcript_index].c_str();
    if (in(tags.gene_index, g_ids_)) tags.gene = g_ids_[tags.gene_index].c_str();
    if (in(tags.gene_index, g_names_)) tags.gene_name = g_names_[tags.gene_index].c_str();
    return tags;
}

/// Stamp the records at ``indices`` (the fragment's tags, then the hit's ZH and ZS) and write them, in order.
void AnnotationPass::write_records(const std::vector<size_t>& indices, const FragmentTags& tags, int32_t primary,
                                   const char* splice_label) {
    for (const size_t i : indices) {
        bam1_t* r = raw_[i].get();
        set_tag_str(r, "ZT", tags.transcript);
        set_tag_str(r, "ZG", tags.gene);
        set_tag_str(r, "ZR", tags.gene_name);
        set_tag_int(r, "ZI", tags.transcript_index);
        set_tag_int(r, "ZJ", tags.gene_index);
        set_tag_int(r, "ZF", tags.flags);
        set_tag_float(r, "ZW", tags.posterior);
        set_tag_str(r, "ZC", tags.frag_class);
        set_tag_int(r, "ZH", primary);
        set_tag_int(r, "ZN", tags.n_candidates);
        set_tag_str(r, "ZS", splice_label);
        set_tag_int(r, "ZL", tags.locus);
        if (sam_write1(out_.get(), hdr_.get(), r) < 0) throw std::runtime_error("Failed to write BAM record");
    }
    n_records_written_ += static_cast<int64_t>(indices.size());
}

}  // namespace

nb::dict BamAnnotationWriter::write(const std::string& bam_path, const std::string& output_path,
                                    const Annotations& annotations, const std::vector<std::string>& t_ids,
                                    const std::vector<std::string>& g_ids,
                                    const std::vector<std::string>& g_names) const {
    AnnotationPass pass(bam_path, output_path, annotations, t_ids, g_ids, g_names, *resolver_, policy_,
                        sj_tag_mode_, skip_duplicates_, include_multimap_);
    return pass.run();
}

}  // namespace rigel
