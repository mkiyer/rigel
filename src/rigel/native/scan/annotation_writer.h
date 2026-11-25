// annotation_writer.h — the annotated-BAM pass (rigel._bam_impl.BamAnnotationWriter): it re-reads the BAM, judges every
// template as the scan judged it, and writes every record back with its fragment's assignment tags.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "ndarray_util.h"
#include "resolve/fragment_resolver.h"
#include "scan/bam_record.h"

namespace rigel {

class BamAnnotationWriter {
public:
    /// The annotation table: one row per annotated fragment, sliced to its rows; a fragment with no row (intergenic)
    /// gets the intergenic default.
    struct Annotations {
        Array1<int64_t> frag_ids;
        Array1<int32_t> best_tid;
        Array1<int32_t> best_gid;
        Array1<uint8_t> tx_flags;
        Array1<float> posterior;
        Array1<int8_t> frag_class;
        Array1<int32_t> n_candidates;
        Array1<int32_t> locus_id;
        int64_t size;
    };

    /// The scan's own settings (pass 1's), so the pass reads and judges the BAM as the scan did.
    BamAnnotationWriter(const FragmentResolver& resolver, const std::string& sj_tag_spec, bool skip_duplicates,
                        bool include_multimap);

    /// Write ``bam_path``'s records to ``output_path``, each stamped with its fragment's tags (the transcript, gene and
    /// gene name ids from ``t_ids``, ``g_ids`` and ``g_names``). Returns the read pairs and records as this pass judged
    /// them, counted as the scan counts them, the digest of its read pairs, and the records written.
    nanobind::dict write(const std::string& bam_path, const std::string& output_path, const Annotations& annotations,
                         const std::vector<std::string>& t_ids, const std::vector<std::string>& g_ids,
                         const std::vector<std::string>& g_names) const;

private:
    const FragmentResolver* resolver_;
    SJTagMode sj_tag_mode_;
    bool skip_duplicates_;
    bool include_multimap_;
};

}  // namespace rigel
