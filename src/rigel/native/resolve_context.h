/**
 * resolve_context.h — FragmentResolver, ResolvedFragment, FragmentAccumulator.
 *
 * Extracted from resolve.cpp so that bam_scanner.cpp can call
 * _resolve_core() directly in C++ and accumulate results into
 * FragmentAccumulator without crossing the Python boundary.
 *
 * Depends on constants.h and cgranges.h.
 */

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

extern "C" {
#include "cgranges.h"
}

#include "constants.h"

namespace nb = nanobind;

namespace rigel {

// ================================================================
// ResolvedFragment — C++ result object exposed to Python
// ================================================================

class ResolvedFragment {
public:
    std::vector<int32_t> t_inds;
    int32_t ambig_strand = 0;
    int32_t splice_type = 0;
    int32_t exon_strand = 0;
    int32_t sj_strand = 0;
    int32_t genomic_footprint = -1;
    int32_t genomic_start = -1;
    int32_t merge_criteria = MC_EMPTY;
    int32_t read_length = 0;
    int32_t chimera_type = CHIMERA_NONE;
    int32_t chimera_gap = -1;
    int32_t num_hits = 1;
    int32_t nm = 0;

    // Per-transcript parallel arrays (same order as t_inds)
    std::vector<int32_t> frag_lengths;      // -1 if missing
    std::vector<int32_t> exon_bp;
    std::vector<int32_t> intron_bp;
    std::vector<int32_t> unambig_intron_bp;

    // --- Properties for Python model training ---

    bool get_is_chimeric() const {
        return chimera_type != CHIMERA_NONE;
    }

    bool get_is_same_strand() const {
        return !ambig_strand;
    }

    bool get_is_strand_qualified() const {
        return splice_type == SPLICE_SPLICED_ANNOT
            && !ambig_strand
            && (exon_strand == STRAND_POS || exon_strand == STRAND_NEG)
            && (sj_strand == STRAND_POS || sj_strand == STRAND_NEG)
            && chimera_type == CHIMERA_NONE;
    }

    int32_t get_first_t_ind() const {
        return t_inds.empty() ? -1 : t_inds[0];
    }

    bool get_has_frag_lengths() const {
        for (auto fl : frag_lengths)
            if (fl > 0) return true;
        return false;
    }

    int32_t get_unique_frag_length() const {
        int32_t val = -1;
        for (auto fl : frag_lengths) {
            if (fl > 0) {
                if (val == -1) val = fl;
                else if (fl != val) return -1;
            }
        }
        return val;
    }

    // Return t_inds as a Python frozenset for compatibility
    nb::object get_t_inds() const {
        nb::set s;
        for (int32_t t : t_inds) s.add(nb::cast(t));
        return nb::frozenset(s);
    }

    // Return frag_lengths as a Python dict for compatibility
    nb::dict get_frag_lengths() const {
        nb::dict d;
        for (size_t i = 0; i < t_inds.size(); i++) {
            if (frag_lengths[i] > 0)
                d[nb::cast(t_inds[i])] = nb::cast(frag_lengths[i]);
        }
        return d;
    }

    // Return overlap_bp as a Python dict for compatibility
    nb::dict get_overlap_bp() const {
        nb::dict d;
        for (size_t i = 0; i < t_inds.size(); i++) {
            d[nb::cast(t_inds[i])] = nb::make_tuple(
                exon_bp[i], intron_bp[i], unambig_intron_bp[i]);
        }
        return d;
    }

    // Build from RawResolveResult (used by resolve_fragment and bam_scanner)
    static ResolvedFragment from_core(RawResolveResult& cr) {
        ResolvedFragment r;
        r.t_inds = std::move(cr.t_inds);
        r.ambig_strand = cr.ambig_strand;
        r.splice_type = cr.splice_type;
        r.exon_strand = cr.exon_strand;
        r.sj_strand = cr.sj_strand;
        r.genomic_footprint = cr.genomic_footprint;
        r.genomic_start = cr.genomic_start;
        r.merge_criteria = cr.merge_criteria;
        r.read_length = cr.read_length;
        r.chimera_type = cr.chimera_type;
        r.chimera_gap = cr.chimera_gap;
        r.exon_bp = std::move(cr.t_exon_bp);
        r.intron_bp = std::move(cr.t_intron_bp);
        r.unambig_intron_bp = std::move(cr.t_unambig_intron_bp);
        // Convert frag_length_map to parallel array
        r.frag_lengths.reserve(r.t_inds.size());
        for (int32_t t : r.t_inds) {
            auto it = cr.frag_length_map.find(t);
            r.frag_lengths.push_back(
                (it != cr.frag_length_map.end()) ? it->second : -1);
        }
        return r;
    }
};

// ================================================================
// FragmentAccumulator — C++ columnar buffer replacing _AccumulatorChunk
// ================================================================

class FragmentAccumulator {
public:
    std::vector<uint8_t>  splice_type_;
    std::vector<uint8_t>  exon_strand_;
    std::vector<uint8_t>  sj_strand_;
    std::vector<uint16_t> num_hits_;
    std::vector<uint8_t>  merge_criteria_;
    std::vector<uint8_t>  chimera_type_;
    std::vector<int32_t>  t_indices_;
    std::vector<int64_t>  t_offsets_;
    std::vector<int32_t>  frag_lengths_;
    std::vector<int32_t>  exon_bp_;
    std::vector<int32_t>  intron_bp_;
    std::vector<int32_t>  unambig_intron_bp_;
    std::vector<int64_t>  frag_id_;
    std::vector<uint32_t> read_length_;
    std::vector<int32_t>  genomic_footprint_;
    std::vector<int32_t>  genomic_start_;
    std::vector<uint16_t> nm_;
    int32_t size_ = 0;

    FragmentAccumulator() {
        t_offsets_.push_back(0);
    }

    void append(const ResolvedFragment& r, int64_t frag_id) {
        splice_type_.push_back(static_cast<uint8_t>(r.splice_type));
        exon_strand_.push_back(static_cast<uint8_t>(r.exon_strand));
        sj_strand_.push_back(static_cast<uint8_t>(r.sj_strand));
        num_hits_.push_back(static_cast<uint16_t>(r.num_hits));
        merge_criteria_.push_back(static_cast<uint8_t>(r.merge_criteria));
        chimera_type_.push_back(static_cast<uint8_t>(r.chimera_type));

        // CSR data — parallel arrays to t_indices
        t_indices_.insert(t_indices_.end(),
                          r.t_inds.begin(), r.t_inds.end());
        frag_lengths_.insert(frag_lengths_.end(),
                             r.frag_lengths.begin(), r.frag_lengths.end());
        exon_bp_.insert(exon_bp_.end(),
                        r.exon_bp.begin(), r.exon_bp.end());
        intron_bp_.insert(intron_bp_.end(),
                          r.intron_bp.begin(), r.intron_bp.end());
        unambig_intron_bp_.insert(unambig_intron_bp_.end(),
                                  r.unambig_intron_bp.begin(),
                                  r.unambig_intron_bp.end());
        t_offsets_.push_back(static_cast<int64_t>(t_indices_.size()));

        frag_id_.push_back(frag_id);
        read_length_.push_back(static_cast<uint32_t>(r.read_length));
        genomic_footprint_.push_back(r.genomic_footprint);
        genomic_start_.push_back(r.genomic_start);
        nm_.push_back(static_cast<uint16_t>(r.nm));
        size_++;
    }

    int32_t get_size() const { return size_; }

    /// Finalize accumulator → dict of raw bytes (for numpy frombuffer).
    /// Also computes ambig_strand from transcript strand mapping.
    nb::dict finalize(const std::vector<int32_t>& t_strand_arr) {
        int32_t n = size_;

        // Compute ambig_strand per fragment from transcript strand array
        std::vector<uint8_t> ambig_strand(n, 0);
        for (int32_t i = 0; i < n; i++) {
            int64_t start = t_offsets_[i];
            int64_t end = t_offsets_[i + 1];
            if (end - start <= 1) continue;  // 0 or 1 transcript: not mixed
            int32_t first_strand = -999;
            for (int64_t j = start; j < end; j++) {
                int32_t t = t_indices_[j];
                if (t >= 0 && t < static_cast<int32_t>(t_strand_arr.size())) {
                    int32_t s = t_strand_arr[t];
                    if (first_strand == -999) {
                        first_strand = s;
                    } else if (s != first_strand) {
                        ambig_strand[i] = 1;
                        break;
                    }
                }
            }
        }

        auto to_bytes = [](const void* data, size_t nbytes) -> nb::bytes {
            if (nbytes == 0)
                return nb::bytes(static_cast<const char*>(nullptr), 0);
            return nb::bytes(
                reinterpret_cast<const char*>(data), nbytes);
        };

        nb::dict result;
        result["splice_type"] = to_bytes(
            splice_type_.data(), splice_type_.size());
        result["exon_strand"] = to_bytes(
            exon_strand_.data(), exon_strand_.size());
        result["sj_strand"] = to_bytes(
            sj_strand_.data(), sj_strand_.size());
        result["num_hits"] = to_bytes(
            num_hits_.data(), num_hits_.size() * sizeof(uint16_t));
        result["merge_criteria"] = to_bytes(
            merge_criteria_.data(), merge_criteria_.size());
        result["chimera_type"] = to_bytes(
            chimera_type_.data(), chimera_type_.size());
        result["t_indices"] = to_bytes(
            t_indices_.data(), t_indices_.size() * sizeof(int32_t));
        result["t_offsets"] = to_bytes(
            t_offsets_.data(), t_offsets_.size() * sizeof(int64_t));
        result["frag_lengths"] = to_bytes(
            frag_lengths_.data(), frag_lengths_.size() * sizeof(int32_t));
        result["exon_bp"] = to_bytes(
            exon_bp_.data(), exon_bp_.size() * sizeof(int32_t));
        result["intron_bp"] = to_bytes(
            intron_bp_.data(), intron_bp_.size() * sizeof(int32_t));
        result["unambig_intron_bp"] = to_bytes(
            unambig_intron_bp_.data(),
            unambig_intron_bp_.size() * sizeof(int32_t));
        result["ambig_strand"] = to_bytes(ambig_strand.data(), ambig_strand.size());
        result["frag_id"] = to_bytes(
            frag_id_.data(), frag_id_.size() * sizeof(int64_t));
        result["read_length"] = to_bytes(
            read_length_.data(), read_length_.size() * sizeof(uint32_t));
        result["genomic_footprint"] = to_bytes(
            genomic_footprint_.data(),
            genomic_footprint_.size() * sizeof(int32_t));
        result["genomic_start"] = to_bytes(
            genomic_start_.data(),
            genomic_start_.size() * sizeof(int32_t));
        result["nm"] = to_bytes(
            nm_.data(), nm_.size() * sizeof(uint16_t));
        result["size"] = nb::cast(size_);

        return result;
    }
};

// ================================================================
// ResolverScratch — per-thread mutable scratch buffers
// ================================================================
// Extracted from FragmentResolver so that multiple threads can call
// _resolve_core() and compute_frag_lengths() concurrently, each with
// its own scratch. The single-threaded path uses FragmentResolver's
// internal scratch_ member.

struct ResolverScratch {
    // Per-transcript bp counters (sparse-cleaned)
    std::vector<int32_t> t_exon_bp;
    std::vector<int32_t> t_transcript_bp;
    std::vector<int32_t> t_unambig_intron_bp;
    std::vector<uint8_t> t_dirty;
    std::vector<int32_t> dirty_indices;

    // cgranges query buffers (reusable per-call)
    int64_t* buf = nullptr;
    int64_t  buf_cap = 0;
    int64_t* sj_buf = nullptr;
    int64_t  sj_buf_cap = 0;

    ResolverScratch() = default;

    explicit ResolverScratch(int32_t n_transcripts)
        : t_exon_bp(n_transcripts, 0),
          t_transcript_bp(n_transcripts, 0),
          t_unambig_intron_bp(n_transcripts, 0),
          t_dirty(n_transcripts, 0)
    {
        dirty_indices.reserve(512);
    }

    ~ResolverScratch() {
        free(buf);
        free(sj_buf);
    }

    // Non-copyable, movable
    ResolverScratch(const ResolverScratch&) = delete;
    ResolverScratch& operator=(const ResolverScratch&) = delete;

    ResolverScratch(ResolverScratch&& o) noexcept
        : t_exon_bp(std::move(o.t_exon_bp)),
          t_transcript_bp(std::move(o.t_transcript_bp)),
          t_unambig_intron_bp(std::move(o.t_unambig_intron_bp)),
          t_dirty(std::move(o.t_dirty)),
          dirty_indices(std::move(o.dirty_indices)),
          buf(o.buf), buf_cap(o.buf_cap),
          sj_buf(o.sj_buf), sj_buf_cap(o.sj_buf_cap)
    {
        o.buf = nullptr; o.buf_cap = 0;
        o.sj_buf = nullptr; o.sj_buf_cap = 0;
    }

    ResolverScratch& operator=(ResolverScratch&& o) noexcept {
        if (this != &o) {
            free(buf); free(sj_buf);
            t_exon_bp = std::move(o.t_exon_bp);
            t_transcript_bp = std::move(o.t_transcript_bp);
            t_unambig_intron_bp = std::move(o.t_unambig_intron_bp);
            t_dirty = std::move(o.t_dirty);
            dirty_indices = std::move(o.dirty_indices);
            buf = o.buf; buf_cap = o.buf_cap;
            sj_buf = o.sj_buf; sj_buf_cap = o.sj_buf_cap;
            o.buf = nullptr; o.buf_cap = 0;
            o.sj_buf = nullptr; o.sj_buf_cap = 0;
        }
        return *this;
    }

    void mark_dirty(int32_t t_idx) {
        if (!t_dirty[t_idx]) {
            t_dirty[t_idx] = 1;
            dirty_indices.push_back(t_idx);
        }
    }

    void clean() {
        for (int32_t t : dirty_indices) {
            t_exon_bp[t] = 0;
            t_transcript_bp[t] = 0;
            t_unambig_intron_bp[t] = 0;
            t_dirty[t] = 0;
        }
        dirty_indices.clear();
    }
};

// ================================================================
// FragmentResolver — holds all index data for C++ resolution
// ================================================================

class FragmentResolver {
public:
    // --- Overlap index (main cgranges) ---
    cgranges_t* cr_ = nullptr;
    std::vector<int8_t>  iv_type_;
    std::vector<int32_t> t_set_data_;
    std::vector<int32_t> t_set_offsets_;

    // --- SJ exact-match map ---
    using SJMap = std::unordered_map<SJKey, std::pair<int32_t, int32_t>, SJKeyHash>;
    SJMap sj_map_;
    std::vector<int32_t> sj_map_data_;

    // --- SJ gap index (for fragment-length computation) ---
    cgranges_t* sj_cr_ = nullptr;
    std::vector<int32_t> sj_t_index_;
    std::vector<int32_t> sj_strand_;

    // --- Transcript metadata ---
    std::vector<int32_t> t_to_g_arr_;
    int32_t n_transcripts_ = 0;

    // --- Gene strand metadata (for BAM scanner model training) ---
    std::vector<int32_t> g_to_strand_arr_;

    // --- Transcript strand metadata (direct lookup, no gene indirection) ---
    std::vector<int32_t> t_strand_arr_;

    // --- Reference name <-> ID ---
    std::unordered_map<std::string, int32_t> ref_to_id_;
    std::vector<std::string> id_to_ref_;

    // --- Per-call scratch (mutable, NOT thread-safe) ---
    // Single-threaded callers use this. Multi-threaded callers pass
    // their own ResolverScratch to the overloaded methods.
    ResolverScratch scratch_;

    // ----------------------------------------------------------------
    FragmentResolver() = default;

    ~FragmentResolver() {
        if (cr_) cr_destroy(cr_);
        if (sj_cr_) cr_destroy(sj_cr_);
    }

    // Non-copyable
    FragmentResolver(const FragmentResolver&) = delete;
    FragmentResolver& operator=(const FragmentResolver&) = delete;

    // ----------------------------------------------------------------
    // Ref-ID helpers
    // ----------------------------------------------------------------

    int32_t get_or_create_ref_id(const std::string& ref) {
        auto it = ref_to_id_.find(ref);
        if (it != ref_to_id_.end()) return it->second;
        int32_t id = static_cast<int32_t>(id_to_ref_.size());
        ref_to_id_[ref] = id;
        id_to_ref_.push_back(ref);
        return id;
    }

    // ----------------------------------------------------------------
    // Build methods — called once at index load time
    // ----------------------------------------------------------------

    void build_overlap_index(
        const std::vector<std::string>& refs,
        const std::vector<int32_t>& starts,
        const std::vector<int32_t>& ends,
        const std::vector<int32_t>& iv_types,
        const std::vector<int32_t>& tset_data,
        const std::vector<int32_t>& tset_offsets)
    {
        size_t n = refs.size();
        for (const auto& ref : refs) get_or_create_ref_id(ref);

        cr_ = cr_init();
        for (size_t i = 0; i < n; i++)
            cr_add(cr_, refs[i].c_str(), starts[i], ends[i],
                   static_cast<int32_t>(i));
        cr_index(cr_);

        iv_type_.resize(n);
        for (size_t i = 0; i < n; i++)
            iv_type_[i] = static_cast<int8_t>(iv_types[i]);

        t_set_data_ = tset_data;
        t_set_offsets_ = tset_offsets;
    }

    void build_sj_map(
        const std::vector<std::string>& refs,
        const std::vector<int32_t>& starts,
        const std::vector<int32_t>& ends,
        const std::vector<int32_t>& strands,
        const std::vector<int32_t>& tset_data,
        const std::vector<int32_t>& tset_offsets)
    {
        sj_map_data_ = tset_data;
        size_t n = refs.size();
        for (size_t i = 0; i < n; i++) {
            int32_t ref_id = get_or_create_ref_id(refs[i]);
            SJKey key{ref_id, starts[i], ends[i], strands[i]};
            sj_map_[key] = {tset_offsets[i],
                            tset_offsets[i + 1] - tset_offsets[i]};
        }
    }

    void build_sj_gap_index(
        const std::vector<std::string>& refs,
        const std::vector<int32_t>& starts,
        const std::vector<int32_t>& ends,
        const std::vector<int32_t>& t_indices,
        const std::vector<int32_t>& strands)
    {
        size_t n = refs.size();
        for (const auto& ref : refs) get_or_create_ref_id(ref);

        sj_cr_ = cr_init();
        for (size_t i = 0; i < n; i++)
            cr_add(sj_cr_, refs[i].c_str(), starts[i], ends[i],
                   static_cast<int32_t>(i));
        cr_index(sj_cr_);

        sj_t_index_ = t_indices;
        sj_strand_ = strands;
    }

    void set_metadata(const std::vector<int32_t>& t_to_g,
                      int32_t n_transcripts) {
        t_to_g_arr_ = t_to_g;
        n_transcripts_ = n_transcripts;

        // Initialize internal scratch buffers
        scratch_ = ResolverScratch(n_transcripts_);
    }

    /// Set gene strand mapping (for BAM scanner model training)
    void set_gene_strands(const std::vector<int32_t>& g_to_strand) {
        g_to_strand_arr_ = g_to_strand;
    }

    /// Set per-transcript strand array (direct lookup, no gene indirection).
    void set_transcript_strands(const std::vector<int32_t>& t_strand) {
        t_strand_arr_ = t_strand;
    }

    nb::dict get_ref_to_id() const {
        nb::dict d;
        for (const auto& [k, v] : ref_to_id_)
            d[nb::cast(k)] = nb::cast(v);
        return d;
    }

    // ----------------------------------------------------------------
    // Fragment-length computation (ports compute_frag_lengths)
    // ----------------------------------------------------------------

    /// Thread-safe overload: uses caller-supplied scratch for sj_buf.
    std::unordered_map<int32_t, int32_t> compute_frag_lengths(
        const std::vector<ExonBlock>& exons,
        const std::vector<IntronBlock>& introns,
        const std::vector<int32_t>& t_inds,
        ResolverScratch& scratch) const
    {
        std::unordered_map<int32_t, int32_t> result;
        if (exons.empty() || t_inds.empty()) return result;

        // Single exon block: shared fragment length
        if (exons.size() == 1) {
            int32_t fl = exons[0].end - exons[0].start;
            if (fl > 0)
                for (int32_t t : t_inds) result[t] = fl;
            return result;
        }

        int32_t ref_id = exons[0].ref_id;

        // Sort exons by start
        std::vector<std::pair<int32_t, int32_t>> sorted_ex;
        sorted_ex.reserve(exons.size());
        for (const auto& e : exons)
            sorted_ex.push_back({e.start, e.end});
        std::sort(sorted_ex.begin(), sorted_ex.end());

        int32_t footprint = sorted_ex.back().second - sorted_ex.front().first;

        // Observed introns on this reference
        std::set<std::pair<int32_t, int32_t>> obs_introns;
        int32_t obs_total = 0;
        for (const auto& intr : introns) {
            if (intr.ref_id == ref_id) {
                auto p = std::make_pair(intr.start, intr.end);
                if (obs_introns.insert(p).second)
                    obs_total += intr.end - intr.start;
            }
        }
        int32_t upper = footprint - obs_total;

        // Gaps between consecutive blocks not in observed introns
        std::vector<std::pair<int32_t, int32_t>> gaps;
        for (size_t i = 0; i + 1 < sorted_ex.size(); i++) {
            int32_t gs = sorted_ex[i].second;
            int32_t ge = sorted_ex[i + 1].first;
            if (gs < ge && obs_introns.find({gs, ge}) == obs_introns.end())
                gaps.push_back({gs, ge});
        }

        if (gaps.empty()) {
            if (upper > 0)
                for (int32_t t : t_inds) result[t] = upper;
            return result;
        }

        // Per-transcript gap correction using SJ gap index
        std::unordered_set<int32_t> t_set(t_inds.begin(), t_inds.end());
        std::unordered_map<int32_t, int32_t> t_gap_size;

        if (ref_id >= 0 && ref_id < static_cast<int32_t>(id_to_ref_.size())) {
            const char* ref_str = id_to_ref_[ref_id].c_str();
            for (const auto& [gs, ge] : gaps) {
                int64_t n = cr_overlap(sj_cr_, ref_str, gs, ge,
                                       &scratch.sj_buf, &scratch.sj_buf_cap);
                for (int64_t i = 0; i < n; i++) {
                    int64_t idx = scratch.sj_buf[i];
                    int32_t hs = cr_start(sj_cr_, idx);
                    int32_t he = cr_end(sj_cr_, idx);
                    int32_t label = cr_label(sj_cr_, idx);
                    if (hs >= gs && he <= ge) {
                        int32_t ti = sj_t_index_[label];
                        if (t_set.count(ti))
                            t_gap_size[ti] += (he - hs);
                    }
                }
            }
        }

        for (int32_t t : t_inds) {
            auto it = t_gap_size.find(t);
            int32_t corr = (it != t_gap_size.end()) ? it->second : 0;
            int32_t sz = upper - corr;
            if (sz > 0) result[t] = sz;
        }
        return result;
    }

    /// Backward-compatible wrapper using internal scratch.
    std::unordered_map<int32_t, int32_t> compute_frag_lengths(
        const std::vector<ExonBlock>& exons,
        const std::vector<IntronBlock>& introns,
        const std::vector<int32_t>& t_inds)
    {
        return compute_frag_lengths(exons, introns, t_inds, scratch_);
    }

    // ----------------------------------------------------------------
    // SJ map lookup helper
    // ----------------------------------------------------------------

    std::vector<int32_t> sj_lookup(int32_t ref_id, int32_t start,
                                   int32_t end, int32_t strand) const {
        SJKey key{ref_id, start, end, strand};
        auto it = sj_map_.find(key);
        if (it != sj_map_.end()) {
            auto [off, cnt] = it->second;
            std::vector<int32_t> out(sj_map_data_.begin() + off,
                                     sj_map_data_.begin() + off + cnt);
            std::sort(out.begin(), out.end());
            return out;
        }
        // Strand-agnostic fallback
        if (strand != STRAND_POS && strand != STRAND_NEG) {
            SJKey pk{ref_id, start, end, STRAND_POS};
            SJKey nk{ref_id, start, end, STRAND_NEG};
            auto pi = sj_map_.find(pk);
            auto ni = sj_map_.find(nk);
            if (pi != sj_map_.end() || ni != sj_map_.end()) {
                std::unordered_set<int32_t> combined;
                if (pi != sj_map_.end()) {
                    auto [off, cnt] = pi->second;
                    for (int32_t k = 0; k < cnt; k++)
                        combined.insert(sj_map_data_[off + k]);
                }
                if (ni != sj_map_.end()) {
                    auto [off, cnt] = ni->second;
                    for (int32_t k = 0; k < cnt; k++)
                        combined.insert(sj_map_data_[off + k]);
                }
                std::vector<int32_t> out(combined.begin(), combined.end());
                std::sort(out.begin(), out.end());
                return out;
            }
        }
        return {};
    }

    // ----------------------------------------------------------------
    // Scratch buffer management (delegates to scratch_)
    // ----------------------------------------------------------------

    void mark_dirty(int32_t t_idx) {
        scratch_.mark_dirty(t_idx);
    }

    void clean_scratch() {
        scratch_.clean();
    }

    // ================================================================
    // _resolve_core — shared logic for resolve() and resolve_fragment()
    // ================================================================

    /// Thread-safe overload: uses caller-supplied scratch buffers.
    bool _resolve_core(
        const std::vector<ExonBlock>& exons,
        const std::vector<IntronBlock>& introns,
        int32_t genomic_footprint,
        RawResolveResult& cr,
        ResolverScratch& scratch)
    {
        int n_exons = static_cast<int>(exons.size());
        if (n_exons == 0) return false;

        cr.genomic_footprint = genomic_footprint;

        // --- Interchromosomal chimera detection ---
        cr.chimera_type = CHIMERA_NONE;
        cr.chimera_gap = -1;
        std::unordered_set<int32_t> ref_set;
        for (const auto& e : exons) ref_set.insert(e.ref_id);
        bool is_interchromosomal = (ref_set.size() > 1);
        if (is_interchromosomal) cr.chimera_type = CHIMERA_TRANS;

        // --- Query each exon block ---
        int n_introns = static_cast<int>(introns.size());
        std::vector<std::vector<int32_t>> exon_t_sets(n_exons);
        std::vector<std::vector<int32_t>> transcript_t_sets(n_exons);
        cr.exon_strand = STRAND_NONE;
        cr.read_length = 0;

        for (int bi = 0; bi < n_exons; bi++) {
            const auto& eb = exons[bi];
            cr.exon_strand |= eb.strand;
            int32_t bstart = eb.start;
            int32_t bend = eb.end;
            cr.read_length += bend - bstart;

            if (eb.ref_id < 0 ||
                eb.ref_id >= static_cast<int32_t>(id_to_ref_.size())) {
                continue;
            }

            std::unordered_set<int32_t> block_exon_t;
            std::unordered_set<int32_t> block_transcript_t;

            const char* ref_str = id_to_ref_[eb.ref_id].c_str();
            int64_t n = cr_overlap(cr_, ref_str, bstart, bend,
                                   &scratch.buf, &scratch.buf_cap);

            for (int64_t hi = 0; hi < n; hi++) {
                int64_t idx = scratch.buf[hi];
                int32_t h_start = cr_start(cr_, idx);
                int32_t h_end   = cr_end(cr_, idx);
                int32_t label   = cr_label(cr_, idx);
                int8_t  itype   = iv_type_[label];

                int32_t off = t_set_offsets_[label];
                int32_t cnt = t_set_offsets_[label + 1] - off;

                if (itype == ITYPE_EXON) {
                    for (int32_t k = 0; k < cnt; k++)
                        block_exon_t.insert(t_set_data_[off + k]);
                    int32_t clo = std::max(bstart, h_start);
                    int32_t chi = std::min(bend, h_end);
                    if (chi > clo) {
                        int32_t bp = chi - clo;
                        for (int32_t k = 0; k < cnt; k++) {
                            int32_t ti = t_set_data_[off + k];
                            scratch.mark_dirty(ti);
                            scratch.t_exon_bp[ti] += bp;
                        }
                    }
                } else if (itype == ITYPE_TRANSCRIPT) {
                    for (int32_t k = 0; k < cnt; k++)
                        block_transcript_t.insert(t_set_data_[off + k]);
                    int32_t clo = std::max(bstart, h_start);
                    int32_t chi = std::min(bend, h_end);
                    if (chi > clo) {
                        int32_t bp = chi - clo;
                        for (int32_t k = 0; k < cnt; k++) {
                            int32_t ti = t_set_data_[off + k];
                            scratch.mark_dirty(ti);
                            scratch.t_transcript_bp[ti] += bp;
                        }
                    }
                } else if (itype == ITYPE_UNAMBIG_INTRON) {
                    int32_t clo = std::max(bstart, h_start);
                    int32_t chi = std::min(bend, h_end);
                    if (chi > clo) {
                        int32_t bp = chi - clo;
                        for (int32_t k = 0; k < cnt; k++) {
                            int32_t ti = t_set_data_[off + k];
                            scratch.mark_dirty(ti);
                            scratch.t_unambig_intron_bp[ti] += bp;
                        }
                    }
                }
            }

            exon_t_sets[bi].assign(block_exon_t.begin(), block_exon_t.end());
            std::sort(exon_t_sets[bi].begin(), exon_t_sets[bi].end());
            transcript_t_sets[bi].assign(block_transcript_t.begin(),
                                         block_transcript_t.end());
            std::sort(transcript_t_sets[bi].begin(),
                      transcript_t_sets[bi].end());
        }

        // --- Intrachromosomal chimera detection ---
        if (!is_interchromosomal) {
            auto cr_res = detect_chimera(exons, exon_t_sets);
            if (cr_res.type != CHIMERA_NONE) {
                cr.chimera_type = cr_res.type;
                cr.chimera_gap = cr_res.gap;
            }
        }

        // --- SJ matching ---
        std::vector<std::vector<int32_t>> sj_t_sets_vec;
        bool has_annotated_sj = false;
        bool has_unannotated_sj = false;
        cr.sj_strand = STRAND_NONE;

        for (int ii = 0; ii < n_introns; ii++) {
            auto sj_t = sj_lookup(introns[ii].ref_id, introns[ii].start,
                                   introns[ii].end, introns[ii].strand);
            if (!sj_t.empty()) {
                sj_t_sets_vec.push_back(std::move(sj_t));
                has_annotated_sj = true;
                cr.sj_strand |= introns[ii].strand;
            } else {
                has_unannotated_sj = true;
            }
        }

        // --- Resolution: merge and classify ---
        bool any_exon = false, any_transcript = false;
        for (const auto& s : exon_t_sets)
            if (!s.empty()) { any_exon = true; break; }
        for (const auto& s : transcript_t_sets)
            if (!s.empty()) { any_transcript = true; break; }

        if (has_annotated_sj && any_exon) {
            auto exon_merge = merge_sets(exon_t_sets);
            auto sj_merge   = merge_sets(sj_t_sets_vec);

            if (!sj_merge.is_empty()) {
                if (!exon_merge.is_empty()) {
                    auto isect = vec_intersect(exon_merge.t_inds,
                                               sj_merge.t_inds);
                    cr.t_inds = isect.empty() ? std::move(sj_merge.t_inds)
                                              : std::move(isect);
                } else {
                    cr.t_inds = std::move(sj_merge.t_inds);
                }
                cr.merge_criteria = sj_merge.criteria;
            } else {
                std::vector<std::vector<int32_t>> all_sets;
                all_sets.insert(all_sets.end(),
                                exon_t_sets.begin(), exon_t_sets.end());
                all_sets.insert(all_sets.end(),
                                sj_t_sets_vec.begin(), sj_t_sets_vec.end());
                auto mr = merge_sets(all_sets);
                cr.t_inds = std::move(mr.t_inds);
                cr.merge_criteria = mr.criteria;
            }
            cr.splice_type = SPLICE_SPLICED_ANNOT;

        } else if (any_exon || any_transcript) {
            std::vector<std::vector<int32_t>> parts;
            int32_t best_criteria = MC_UNION;

            if (any_exon) {
                auto em = merge_sets(exon_t_sets);
                if (!em.is_empty()) {
                    parts.push_back(std::move(em.t_inds));
                    best_criteria = em.criteria;
                }
            }
            if (any_transcript) {
                auto tm = merge_sets(transcript_t_sets);
                if (!tm.is_empty())
                    parts.push_back(std::move(tm.t_inds));
            }

            if (!parts.empty()) {
                std::unordered_set<int32_t> all;
                for (const auto& p : parts)
                    for (int32_t v : p) all.insert(v);
                cr.t_inds.assign(all.begin(), all.end());
                std::sort(cr.t_inds.begin(), cr.t_inds.end());
                cr.merge_criteria =
                    (parts.size() == 1) ? best_criteria : MC_UNION;
            }

            cr.splice_type = has_unannotated_sj ? SPLICE_SPLICED_UNANNOT
                                                : SPLICE_UNSPLICED;

        } else {
            scratch.clean();
            return false;
        }

        if (cr.t_inds.empty()) {
            scratch.clean();
            return false;
        }

        // --- Overlap profiles (parallel to t_inds) ---
        std::unordered_set<int32_t> all_overlap_t;
        for (const auto& s : exon_t_sets)
            for (int32_t v : s) all_overlap_t.insert(v);
        for (const auto& s : transcript_t_sets)
            for (int32_t v : s) all_overlap_t.insert(v);

        bool any_overlap = !all_overlap_t.empty() && cr.read_length > 0;
        cr.t_exon_bp.reserve(cr.t_inds.size());
        cr.t_intron_bp.reserve(cr.t_inds.size());
        cr.t_unambig_intron_bp.reserve(cr.t_inds.size());
        for (int32_t t : cr.t_inds) {
            if (any_overlap && all_overlap_t.count(t)) {
                cr.t_exon_bp.push_back(scratch.t_exon_bp[t]);
                cr.t_intron_bp.push_back(
                    std::max(scratch.t_transcript_bp[t] - scratch.t_exon_bp[t], 0));
                cr.t_unambig_intron_bp.push_back(scratch.t_unambig_intron_bp[t]);
            } else if (!any_overlap) {
                cr.t_exon_bp.push_back(0);
                cr.t_intron_bp.push_back(cr.read_length);
                cr.t_unambig_intron_bp.push_back(cr.read_length);
            } else {
                cr.t_exon_bp.push_back(cr.read_length);
                cr.t_intron_bp.push_back(0);
                cr.t_unambig_intron_bp.push_back(0);
            }
        }

        // --- ambig_strand ---
        int32_t first_strand = -999;
        bool mixed = false;
        for (int32_t t : cr.t_inds) {
            if (t >= 0 && t < static_cast<int32_t>(t_strand_arr_.size())) {
                int32_t s = t_strand_arr_[t];
                if (first_strand == -999) {
                    first_strand = s;
                } else if (s != first_strand) {
                    mixed = true;
                    break;
                }
            }
        }
        cr.ambig_strand = mixed ? 1 : 0;

        // --- Fragment lengths ---
        if (cr.chimera_type == CHIMERA_NONE) {
            cr.frag_length_map = compute_frag_lengths(exons, introns, cr.t_inds, scratch);
        }

        // --- Genomic start ---
        cr.genomic_start = exons[0].start;
        for (const auto& e : exons)
            cr.genomic_start = std::min(cr.genomic_start, e.start);

        scratch.clean();
        return true;
    }

    /// Backward-compatible wrapper using internal scratch.
    bool _resolve_core(
        const std::vector<ExonBlock>& exons,
        const std::vector<IntronBlock>& introns,
        int32_t genomic_footprint,
        RawResolveResult& cr)
    {
        return _resolve_core(exons, introns, genomic_footprint, cr, scratch_);
    }

    // ================================================================
    // resolve — legacy entry point (returns 13-element tuple)
    // ================================================================

    nb::object resolve(
        const std::vector<int32_t>& exon_ref_ids,
        const std::vector<int32_t>& exon_starts,
        const std::vector<int32_t>& exon_ends,
        const std::vector<int32_t>& exon_strands,
        const std::vector<int32_t>& intron_ref_ids,
        const std::vector<int32_t>& intron_starts,
        const std::vector<int32_t>& intron_ends,
        const std::vector<int32_t>& intron_strands,
        int32_t genomic_footprint)
    {
        int n_exons = static_cast<int>(exon_ref_ids.size());
        if (n_exons == 0) return nb::none();

        std::vector<ExonBlock> exons(n_exons);
        for (int i = 0; i < n_exons; i++) {
            exons[i] = {exon_ref_ids[i], exon_starts[i],
                        exon_ends[i], exon_strands[i]};
        }
        int n_introns = static_cast<int>(intron_ref_ids.size());
        std::vector<IntronBlock> introns(n_introns);
        for (int i = 0; i < n_introns; i++) {
            introns[i] = {intron_ref_ids[i], intron_starts[i],
                          intron_ends[i], intron_strands[i]};
        }

        RawResolveResult cr;
        if (!_resolve_core(exons, introns, genomic_footprint, cr))
            return nb::none();

        // Convert to legacy tuple format
        nb::list t_inds_list;
        for (int32_t t : cr.t_inds) t_inds_list.append(t);

        nb::dict frag_lengths_dict;
        for (const auto& [k, v] : cr.frag_length_map)
            frag_lengths_dict[nb::cast(k)] = nb::cast(v);

        nb::dict overlap_bp_dict;
        for (size_t i = 0; i < cr.t_inds.size(); i++) {
            overlap_bp_dict[nb::cast(cr.t_inds[i])] = nb::make_tuple(
                cr.t_exon_bp[i], cr.t_intron_bp[i],
                cr.t_unambig_intron_bp[i]);
        }

        return nb::make_tuple(
            t_inds_list,
            cr.ambig_strand,
            cr.splice_type,
            cr.exon_strand,
            cr.sj_strand,
            frag_lengths_dict,
            cr.genomic_footprint,
            cr.genomic_start,
            cr.merge_criteria,
            overlap_bp_dict,
            cr.read_length,
            cr.chimera_type,
            cr.chimera_gap
        );
    }

    // ================================================================
    // resolve_fragment — accepts Python Fragment, returns ResolvedFragment
    // ================================================================

    nb::object resolve_fragment(nb::object frag) {
        nb::object exons_obj = frag.attr("exons");
        size_t n_exons = nb::len(exons_obj);
        if (n_exons == 0) return nb::none();

        std::vector<ExonBlock> exons(n_exons);
        for (size_t i = 0; i < n_exons; i++) {
            nb::object e = exons_obj[nb::int_(i)];
            std::string ref = nb::cast<std::string>(e[nb::int_(0)]);
            auto it = ref_to_id_.find(ref);
            exons[i] = {
                (it != ref_to_id_.end()) ? it->second : -1,
                nb::cast<int32_t>(e[nb::int_(1)]),
                nb::cast<int32_t>(e[nb::int_(2)]),
                nb::cast<int32_t>(e[nb::int_(3)])
            };
        }

        nb::object introns_obj = frag.attr("introns");
        size_t n_introns = nb::len(introns_obj);
        std::vector<IntronBlock> introns(n_introns);
        for (size_t i = 0; i < n_introns; i++) {
            nb::object intr = introns_obj[nb::int_(i)];
            std::string ref = nb::cast<std::string>(intr[nb::int_(0)]);
            auto it = ref_to_id_.find(ref);
            introns[i] = {
                (it != ref_to_id_.end()) ? it->second : -1,
                nb::cast<int32_t>(intr[nb::int_(1)]),
                nb::cast<int32_t>(intr[nb::int_(2)]),
                nb::cast<int32_t>(intr[nb::int_(3)])
            };
        }

        int32_t genomic_footprint = nb::cast<int32_t>(
            frag.attr("genomic_footprint"));

        RawResolveResult cr;
        if (!_resolve_core(exons, introns, genomic_footprint, cr))
            return nb::none();

        return nb::cast(ResolvedFragment::from_core(cr));
    }
};

}  // namespace rigel
