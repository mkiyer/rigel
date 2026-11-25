// resolved_fragment.h — what the resolver concludes about one fragment (ResolvedFragment), and the columnar chunk the
// scan gathers the rows into for rigel.buffer (FragmentChunk).
#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>

#include "constants.h"
#include "fragment_blocks.h"
#include "ndarray_util.h"

namespace nb = nanobind;

namespace rigel {

/// One fragment resolved against the annotation (`FragmentResolver::resolve`), with the hit's own values the scan
/// adds: the buffer's row, and what the scan offers calibration.
struct ResolvedFragment {
    // The candidate transcripts, sorted and unique; per candidate, the fragment's length on the candidate's own path
    // (-1 when the path has no base) and its aligned bases on the candidate's exons.
    std::vector<int32_t> t_inds;
    std::vector<int32_t> frag_lengths;
    std::vector<int32_t> exon_bp;

    int32_t splice_type = SPLICE_UNSPLICED;
    int32_t align_strand = STRAND_NONE;  // the blocks' strands, OR-ed
    int32_t sj_strand = STRAND_NONE;     // the observed junctions' motif strands, OR-ed
    int32_t ambig_strand = 0;            // 1 when the candidates lie on both strands
    int32_t read_length = 0;             // the blocks' aligned bases
    int32_t genomic_footprint = -1;      // the molecule's extent, the mate gap included: its bases
    int32_t genomic_start = -1;          // and its leftmost base

    int32_t genomic_end() const noexcept { return genomic_start + genomic_footprint; }

    // The per-junction strand table's key: the leftmost observed junction (-1 when none), completed by `sj_strand`. A
    // strand-qualified fragment carries one strand bit (mates whose XS/ts tags disagree OR to STRAND_AMBIGUOUS and are
    // rejected), so all its junctions agree, and crediting one junction per fragment keeps the 2×2 marginal identity.
    // Every observed junction is annotated (the scan's annotated-path rule), so the pick is an annotated junction.
    int32_t sj_key_ref = -1;
    int32_t sj_key_start = -1;
    int32_t sj_key_end = -1;

    // Every explanation of the unsequenced gaps, flat: hypothesis h implies introns
    // [gap_intron_offsets[h], gap_intron_offsets[h+1]) and carries strand gap_sj_strand[h]. The empty hypothesis is the
    // unspliced, genomic one; the synthetic span is this hypothesis, never a separate one.
    std::vector<IntronBlock> gap_introns;
    std::vector<int32_t> gap_intron_offsets;  // n_hypotheses + 1, starting at 0
    std::vector<int32_t> gap_sj_strand;       // n_hypotheses

    // The hit's own, set by the scan: the template's usable hits (`judge_template`'s count, which can pass 16 bits in
    // a repeat), and the hit's NM and soft-clipped bases over its records.
    int32_t num_hits = 1;
    int32_t nm = 0;
    int32_t clip = 0;

    int32_t n_gap_hypotheses() const {
        return gap_intron_offsets.empty() ? 0 : static_cast<int32_t>(gap_intron_offsets.size()) - 1;
    }

    bool is_same_strand() const { return !ambig_strand; }

    /// A fragment the spliced strand model trains on: annotated-spliced, on one candidate strand, with one aligned
    /// strand and one motif strand.
    bool is_strand_qualified() const {
        return splice_type == SPLICE_SPLICED_ANNOT && !ambig_strand &&
               (align_strand == STRAND_POS || align_strand == STRAND_NEG) &&
               (sj_strand == STRAND_POS || sj_strand == STRAND_NEG);
    }
};

/// The rows of one chunk, column by column, as rigel.buffer reads them.
class FragmentChunk {
public:
    FragmentChunk() {
        t_offsets_.push_back(0);
        place_offsets_.push_back(0);
    }

    /// Pre-allocate for ``n_fragments`` rows and ``n_candidates`` CSR entries (transcript candidates over all rows).
    void reserve(int64_t n_fragments, int64_t n_candidates) {
        splice_type_.reserve(n_fragments);
        align_strand_.reserve(n_fragments);
        ambig_strand_.reserve(n_fragments);
        num_hits_.reserve(n_fragments);
        frag_id_.reserve(n_fragments);
        read_length_.reserve(n_fragments);
        genomic_footprint_.reserve(n_fragments);
        genomic_start_.reserve(n_fragments);
        nm_.reserve(n_fragments);
        clip_.reserve(n_fragments);
        t_indices_.reserve(n_candidates);
        t_offsets_.reserve(n_fragments + 1);
        frag_lengths_.reserve(n_candidates);
        exon_bp_.reserve(n_candidates);
    }

    /// One object of the next appended row's place: call before `append`.
    void add_place_mass(int8_t kind, int64_t index, double mass) {
        place_kind_.push_back(kind);
        place_index_.push_back(index);
        place_mass_.push_back(static_cast<float>(mass));
    }

    void append(const ResolvedFragment& r, int64_t frag_id) {
        splice_type_.push_back(static_cast<uint8_t>(r.splice_type));
        align_strand_.push_back(static_cast<uint8_t>(r.align_strand));
        ambig_strand_.push_back(static_cast<uint8_t>(r.ambig_strand));
        num_hits_.push_back(r.num_hits);

        t_indices_.insert(t_indices_.end(), r.t_inds.begin(), r.t_inds.end());
        frag_lengths_.insert(frag_lengths_.end(), r.frag_lengths.begin(), r.frag_lengths.end());
        for (int32_t bp : r.exon_bp) exon_bp_.push_back(checked_u16(bp, "exon_bp"));
        t_offsets_.push_back(static_cast<int32_t>(t_indices_.size()));

        frag_id_.push_back(frag_id);
        read_length_.push_back(checked_u16(r.read_length, "read_length"));
        genomic_footprint_.push_back(r.genomic_footprint);
        genomic_start_.push_back(r.genomic_start);
        nm_.push_back(checked_u16(r.nm, "nm"));
        clip_.push_back(checked_u16(r.clip, "clip"));
        place_offsets_.push_back(static_cast<int64_t>(place_kind_.size()));
        size_++;
    }

    int32_t size() const { return size_; }

    /// Move the columns into capsule-backed numpy arrays, leaving the accumulator as a new one.
    nb::dict finalize() {
        nb::dict result;
        result["splice_type"]       = vec_to_ndarray(std::move(splice_type_));
        result["align_strand"]      = vec_to_ndarray(std::move(align_strand_));
        result["num_hits"]          = vec_to_ndarray(std::move(num_hits_));
        result["t_indices"]         = vec_to_ndarray(std::move(t_indices_));
        result["t_offsets"]         = vec_to_ndarray(std::move(t_offsets_));
        result["frag_lengths"]      = vec_to_ndarray(std::move(frag_lengths_));
        result["exon_bp"]           = vec_to_ndarray(std::move(exon_bp_));
        result["ambig_strand"]      = vec_to_ndarray(std::move(ambig_strand_));
        result["frag_id"]           = vec_to_ndarray(std::move(frag_id_));
        result["read_length"]       = vec_to_ndarray(std::move(read_length_));
        result["genomic_footprint"] = vec_to_ndarray(std::move(genomic_footprint_));
        result["genomic_start"]     = vec_to_ndarray(std::move(genomic_start_));
        result["nm"]                = vec_to_ndarray(std::move(nm_));
        result["clip"]              = vec_to_ndarray(std::move(clip_));
        result["place_offsets"]     = vec_to_ndarray(std::move(place_offsets_));
        result["place_kind"]        = vec_to_ndarray(std::move(place_kind_));
        result["place_index"]       = vec_to_ndarray(std::move(place_index_));
        result["place_mass"]        = vec_to_ndarray(std::move(place_mass_));
        result["size"]              = nb::cast(size_);
        *this = FragmentChunk();  // the CSR offsets' leading 0s back, so the next row starts a valid chunk
        return result;
    }

private:
    /// ``value`` as a uint16 column entry; one that does not fit throws, naming the column, rather than wrapping.
    static uint16_t checked_u16(int32_t value, const char* column) {
        if (value < 0 || value > static_cast<int32_t>(std::numeric_limits<uint16_t>::max())) {
            throw std::runtime_error(std::string("Fragment buffer column '") + column + "' cannot store value " +
                                     std::to_string(value) + " as uint16");
        }
        return static_cast<uint16_t>(value);
    }

    std::vector<uint8_t>  splice_type_;
    std::vector<uint8_t>  align_strand_;
    std::vector<uint8_t>  ambig_strand_;
    std::vector<int32_t>  num_hits_;
    std::vector<int32_t>  t_indices_;
    std::vector<int32_t>  t_offsets_;
    std::vector<int32_t>  frag_lengths_;
    std::vector<uint16_t> exon_bp_;
    std::vector<int64_t>  frag_id_;
    std::vector<uint16_t> read_length_;
    std::vector<int32_t>  genomic_footprint_;
    std::vector<int32_t>  genomic_start_;
    std::vector<uint16_t> nm_;
    std::vector<uint16_t> clip_;
    // Each row's place masses, CSR: the deposit rule's mass at every object (region, boundary or sj, flat ids) a
    // multimapper's genic place deposits on, which its capture weight reads. Empty on every other row.
    std::vector<int64_t>  place_offsets_;
    std::vector<int8_t>   place_kind_;
    std::vector<int64_t>  place_index_;
    std::vector<float>    place_mass_;
    int32_t size_ = 0;
};

}  // namespace rigel
