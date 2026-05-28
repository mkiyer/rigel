/**
 * accumulator.h - Per-worker fractional calibration observation accumulator.
 *
 * Phase 3 fractional cutover. The accumulator routes per-fragment
 * fractional mass into:
 *   - region_counts[R, 12]    : float32, the 12-channel signal
 *     (compartment x splice x strand) per fine region.
 *   - fl_pool_mass[6, 1024]   : float64, six named FL pools keyed by
 *     (coarse_class_of_region, contained|boundary). Unspliced, strand-
 *     collapsed. gDNA FL aggregate is the sum of the four non-EXON pools.
 *   - channel_mass[12]        : float64, global per-channel mass.
 *   - signature_mass[16]      : float64, global per-signature mass.
 *   - counters                : observed, excluded_*, unobserved, etc.
 *
 * The hot path (``observe()``) is allocation-free after warm-up: it
 * reuses scratch vectors for region hits and per-block overlap. The
 * one-fragment-one-unit invariant holds exactly: total mass added to
 * region_counts for any single observed fragment is 1.0 (modulo
 * float32 rounding) when total_aligned_bp > 0.
 *
 * Compartment routing (per region within a fragment):
 *   - cross_left  = (frag_start < region.start)
 *   - cross_right = (frag_end   > region.end)
 *   - both true   -> fully spans; split w in half between BOUNDARY_LEFT
 *                    and BOUNDARY_RIGHT
 *   - cross_left  -> BOUNDARY_LEFT  (mass = w)
 *   - cross_right -> BOUNDARY_RIGHT (mass = w)
 *   - neither     -> CONTAINED      (mass = w)
 *
 * FL routing (only for unspliced fragments with fl_idx >= 0):
 *   - compartment determines contained vs boundary pool slot
 *   - region signature determines coarse class (INTERGENIC/INTRON/EXON)
 *   - mass is strand-collapsed before adding to fl_pool_mass
 *   - boundary_left and boundary_right both feed the same boundary pool
 *
 * splicing_anchor_tolerance is intentionally absent: in Phase 3 K is
 * resolver-only (implicit-splice slack). The accumulator does not gate
 * boundary mass by overlap thresholds.
 */

#pragma once

#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "constants.h"
#include "calibration/region_index.h"
#include "calibration/region_signature.h"

namespace rigel::calibration {

struct CalibrationPayload {
    static constexpr int32_t kFlBins   = 1024;
    static constexpr int32_t kChannels = region_signature::kNChannels;       // 12
    static constexpr int32_t kSigs     = region_signature::kNSignatures;     // 16
    static constexpr int32_t kFlPools  = region_signature::kNFlPools;        // 6

    // Per-region per-channel fractional mass: shape (n_regions, 12) row-major.
    std::vector<float> region_counts;

    // Global per-channel and per-signature mass (sums over regions and FL).
    std::array<double, kChannels> channel_mass {};
    std::array<double, kSigs>     signature_mass {};

    // Six named FL pools, float64. Unspliced, strand-collapsed.
    // Layout: row pool, column fl_idx. Index = pool * kFlBins + fl_idx.
    std::vector<double>             fl_pool_mass;       // size = 6 * 1024
    std::array<double, kFlPools>    fl_pool_total {};

    // Per-region physical fragment support counts, partitioned by the
    // fragment's splice class. A fragment is counted at most once per
    // region (compartment-agnostic), and only when it contributes
    // positive overlap mass to that region. These are the effective
    // sample sizes consumed by the EB exposure model.
    std::vector<std::uint64_t> region_unspliced_support;  // size = n_regions
    std::vector<std::uint64_t> region_spliced_support;    // size = n_regions

    // Counters
    int64_t n_observed                 = 0;
    int64_t n_excluded_multimap        = 0;
    int64_t n_excluded_chimera         = 0;
    int64_t n_excluded_artifact        = 0;
    int64_t n_excluded_strand_ambig    = 0;
    int64_t n_unobserved               = 0;
    int64_t n_unannotated_ref          = 0;  // observed but no region hits
    int64_t n_fl_unavailable           = 0;  // unspliced observed w/ fl_idx out of range
};

class CalibrationAccumulator {
public:
    explicit CalibrationAccumulator(int64_t n_regions) {
        if (n_regions < 0) {
            throw std::invalid_argument(
                "CalibrationAccumulator: n_regions must be >= 0");
        }
        n_regions_ = n_regions;
        payload_.region_counts.assign(
            static_cast<size_t>(n_regions) * CalibrationPayload::kChannels, 0.0f);
        payload_.fl_pool_mass.assign(
            static_cast<size_t>(CalibrationPayload::kFlPools) *
                CalibrationPayload::kFlBins,
            0.0);
        payload_.region_unspliced_support.assign(
            static_cast<size_t>(n_regions), 0);
        payload_.region_spliced_support.assign(
            static_cast<size_t>(n_regions), 0);
        // Warm reserve to amortize the first ~10 fragments' allocations.
        block_hits_.reserve(16);
        touched_regions_.reserve(16);
    }

    /**
     * Hot-path observation. ``exons`` is the chosen fragment's
     * R1UR2 exon blocks (sorted, single-ref). The caller must have
     * gated on: unique mapper, non-chimeric, non-artifact splice.
     *
     * The accumulator increments ``n_observed`` for every accepted
     * fragment (even if no region overlaps). Fragments whose every
     * block falls outside any region increment ``n_unannotated_ref``.
     * Fragments whose ``fragment_strand`` is neither STRAND_POS nor
     * STRAND_NEG are routed to ``n_excluded_strand_ambig`` and emit no
     * mass.
     */
    void observe(int8_t splice_type,
                 int32_t ref_id,
                 int64_t frag_start,
                 int64_t frag_end,
                 const ExonBlock* exons,
                 int32_t n_exons,
                 int8_t fragment_strand,
                 const RegionIndex& regions);

    inline void note_multimap()        { ++payload_.n_excluded_multimap; }
    inline void note_chimera()         { ++payload_.n_excluded_chimera;  }
    inline void note_artifact()        { ++payload_.n_excluded_artifact; }
    inline void note_unobserved()      { ++payload_.n_unobserved;        }
    inline void note_strand_ambig()    { ++payload_.n_excluded_strand_ambig; }

    /// Worker merge: element-wise additive into *this.
    void merge_from(const CalibrationAccumulator& other);

    const CalibrationPayload& payload() const { return payload_; }
    CalibrationPayload&       payload()       { return payload_; }

    int64_t n_regions() const { return n_regions_; }

private:
    CalibrationPayload payload_;
    int64_t n_regions_ = 0;

    // Per-block scratch (reused across observe() calls).
    std::vector<int32_t> block_hits_;
    // Per-fragment scratch: region IDs that received positive overlap
    // mass from any block of the current fragment. Sorted+uniqued at
    // end of observe() to drive the per-region support increment.
    std::vector<int32_t> touched_regions_;
};

}  // namespace rigel::calibration
