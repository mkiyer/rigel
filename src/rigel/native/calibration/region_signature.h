#pragma once

#include <cstdint>

namespace rigel::calibration::region_signature {

inline constexpr std::uint8_t kBitIntronPos = 0x8;
inline constexpr std::uint8_t kBitIntronNeg = 0x4;
inline constexpr std::uint8_t kBitExonPos = 0x2;
inline constexpr std::uint8_t kBitExonNeg = 0x1;
inline constexpr int kNSignatures = 16;

inline constexpr std::uint8_t kRegionTypeIntergenic = 0;
inline constexpr std::uint8_t kRegionTypeIntron = 1;
inline constexpr std::uint8_t kRegionTypeExon = 2;

inline constexpr std::uint8_t kRegionStrandNone = 0;
inline constexpr std::uint8_t kRegionStrandPos = 1;
inline constexpr std::uint8_t kRegionStrandNeg = 2;
inline constexpr std::uint8_t kRegionStrandAmbig = 3;

inline constexpr int kCompartmentContained = 0;
inline constexpr int kCompartmentBoundaryLeft = 1;
inline constexpr int kCompartmentBoundaryRight = 2;
inline constexpr int kNCompartments = 3;

inline constexpr int kSpliceUnspliced = 0;
inline constexpr int kSpliceSpliced = 1;
inline constexpr int kNSpliceStates = 2;

inline constexpr int kChannelStrandPos = 0;
inline constexpr int kChannelStrandNeg = 1;
inline constexpr int kNChannelStrands = 2;
inline constexpr int kNChannels = kNCompartments * kNSpliceStates * kNChannelStrands;

inline constexpr std::uint8_t pack_signature(bool intron_pos,
                                             bool intron_neg,
                                             bool exon_pos,
                                             bool exon_neg) {
    return static_cast<std::uint8_t>(
        (intron_pos ? kBitIntronPos : 0) |
        (intron_neg ? kBitIntronNeg : 0) |
        (exon_pos ? kBitExonPos : 0) |
        (exon_neg ? kBitExonNeg : 0));
}

inline constexpr bool has_intron_flag(std::uint8_t signature) {
    return (signature & (kBitIntronPos | kBitIntronNeg)) != 0;
}

inline constexpr bool has_exon_flag(std::uint8_t signature) {
    return (signature & (kBitExonPos | kBitExonNeg)) != 0;
}

inline constexpr std::uint8_t coarse_type_from_signature(std::uint8_t signature) {
    if (has_exon_flag(signature)) return kRegionTypeExon;
    if (has_intron_flag(signature)) return kRegionTypeIntron;
    return kRegionTypeIntergenic;
}

inline constexpr std::uint8_t coarse_strand_from_signature(std::uint8_t signature) {
    const bool has_pos = (signature & (kBitIntronPos | kBitExonPos)) != 0;
    const bool has_neg = (signature & (kBitIntronNeg | kBitExonNeg)) != 0;
    if (has_pos && has_neg) return kRegionStrandAmbig;
    if (has_pos) return kRegionStrandPos;
    if (has_neg) return kRegionStrandNeg;
    return kRegionStrandNone;
}

inline constexpr bool is_ambiguous_signature(std::uint8_t signature) {
    return coarse_strand_from_signature(signature) == kRegionStrandAmbig;
}

inline constexpr int channel_index(int compartment, int splice_idx, int strand_idx) {
    return compartment * (kNSpliceStates * kNChannelStrands) + splice_idx * 2 + strand_idx;
}

// FL pool enum used by the Phase 3 fractional accumulator.
// Indexed via fl_pool_index(signature, compartment): boundary_left and
// boundary_right collapse to a single BOUNDARY pool because the
// receiving region's coarse class is what matters for gDNA FL
// extraction; the per-side discrimination has already happened via
// mass routing.
inline constexpr int kFlPoolIntergenicContained = 0;
inline constexpr int kFlPoolIntergenicBoundary  = 1;
inline constexpr int kFlPoolIntronicContained   = 2;
inline constexpr int kFlPoolIntronicBoundary    = 3;
inline constexpr int kFlPoolExonicContained     = 4;
inline constexpr int kFlPoolExonicBoundary      = 5;
inline constexpr int kNFlPools                  = 6;

inline constexpr int fl_pool_index(std::uint8_t signature, int compartment) {
    // coarse: 0=INTERGENIC, 1=INTRON, 2=EXON.  Ambiguous intron 0xC
    // folds into INTRON; any signature with an exon bit (pure exon or
    // mixed exon+intron) folds into EXON.
    const int coarse = static_cast<int>(coarse_type_from_signature(signature));
    const int contained_offset =
        (compartment == kCompartmentContained) ? 0 : 1;
    return coarse * 2 + contained_offset;
}

}  // namespace rigel::calibration::region_signature
