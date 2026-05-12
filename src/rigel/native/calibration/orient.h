/**
 * orient.h — Calibration strand-orientation helper.
 *
 * Numeric values must match src/rigel/calibration/_orient.py.
 */

#pragma once

#include <cstdint>

#include "constants.h"

namespace rigel::calibration::orient {

static constexpr uint8_t SAME  = 0;
static constexpr uint8_t OPP   = 1;
static constexpr uint8_t UNINF = 2;
static constexpr uint8_t N     = 3;

inline uint8_t classify(uint8_t region_strand, int32_t fragment_strand) {
    const bool region_ok =
        region_strand == static_cast<uint8_t>(STRAND_POS) ||
        region_strand == static_cast<uint8_t>(STRAND_NEG);
    const bool fragment_ok =
        fragment_strand == STRAND_POS || fragment_strand == STRAND_NEG;
    if (!region_ok || !fragment_ok) return UNINF;
    return fragment_strand == static_cast<int32_t>(region_strand) ? SAME : OPP;
}

}  // namespace rigel::calibration::orient