// labels.h — the names codes take in the scan's stats export and the annotated BAM's tags.
#pragma once

#include <cstdint>

#include "constants.h"

namespace rigel {

/// A `SpliceType` as the census keys and the ZS tag name it: the member name lower-cased, so `rigel.splice` derives the
/// same key in Python without a second table.
inline const char* splice_type_label(int32_t code) {
    switch (code) {
        case SPLICE_UNSPLICED:       return "unspliced";
        case SPLICE_SPLICED_ANNOT:   return "spliced_annot";
        case SPLICE_IMPLICIT:        return "spliced_implicit";
        case SPLICE_ARTIFACT:        return "splice_artifact";
        default:                     return "unknown";
    }
}

/// Why a hit is not usable (`HIT_*`), as the stats export names it.
inline const char* hit_verdict_label(int32_t verdict) {
    switch (verdict) {
        case HIT_REFERENCE:      return "reference";
        case HIT_STRAND:         return "strand";
        case HIT_OUTWARD:        return "outward";
        case HIT_LENGTH:         return "length";
        case HIT_MISSING_MATE:   return "missing_mate";
        case HIT_CONTRADICTORY:  return "contradictory";
        case HIT_NOVEL_JUNCTION: return "novel_junction";
        case HIT_UNSUPPORTED:    return "unsupported";
        case HIT_OFF_INDEX:      return "off_index";
        case HIT_CHIMERIC:       return "chimeric";
        default:                 return "unknown";
    }
}

/// Why a template is withheld (`WITHHELD_*`): the annotated BAM's ZC and the stats export's name.
inline const char* withheld_label(int32_t outcome) {
    switch (outcome) {
        case WITHHELD_UNPAIRED:           return "unpaired";
        case WITHHELD_UNANNOTATED_SPLICE: return "unannotated_splice";
        case WITHHELD_OFF_INDEX:          return "off_index";
        case WITHHELD_CHIMERIC:           return "chimeric";
        default:                          return "unknown";
    }
}

}  // namespace rigel
