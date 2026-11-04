// fragment_blocks.h — a fragment's geometry as the scan, the resolver and the accumulator share it: its aligned blocks,
// the introns between them, and the junction key.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>

namespace rigel {

/// One aligned block of a fragment, on the fragment's strand.
struct ExonBlock {
    int32_t ref_id;
    int32_t start;
    int32_t end;
    int32_t strand;
};

/// One intron: a junction a read observed (its motif strand) or one a hypothesis implies (its transcript's strand).
struct IntronBlock {
    int32_t ref_id;
    int32_t start;
    int32_t end;
    int32_t strand;
};

/// Blocks ordered by reference, then coordinates, then strand.
inline bool exon_block_less(const ExonBlock& a, const ExonBlock& b) {
    return std::tie(a.ref_id, a.start, a.end, a.strand) < std::tie(b.ref_id, b.start, b.end, b.strand);
}

/// A fragment's extent [start, end): its leftmost block start to its rightmost block end, the mate gap included.
/// ``blocks`` is non-empty.
inline std::pair<int32_t, int32_t> fragment_extent(const std::vector<ExonBlock>& blocks) {
    int32_t start = blocks.front().start, end = blocks.front().end;
    for (const auto& b : blocks) {
        start = std::min(start, b.start);
        end = std::max(end, b.end);
    }
    return {start, end};
}

/// Junctions as the resolver reads them: ascending by reference, coordinates and motif strand, each once.
inline void sort_unique_junctions(std::vector<IntronBlock>& junctions) {
    const auto key = [](const IntronBlock& x) { return std::tie(x.ref_id, x.start, x.end, x.strand); };
    std::sort(junctions.begin(), junctions.end(),
              [&](const IntronBlock& x, const IntronBlock& y) { return key(x) < key(y); });
    junctions.erase(std::unique(junctions.begin(), junctions.end(),
                                [&](const IntronBlock& x, const IntronBlock& y) { return key(x) == key(y); }),
                    junctions.end());
}

/// Fold one int32 into a running hash (the boost hash_combine step). 0x9e3779b9 is 2^32/φ, the golden-ratio increment
/// that spreads consecutive values across the bits.
inline void hash_combine(size_t& h, int32_t v) {
    h ^= std::hash<int32_t>()(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
}

/// A junction by reference, coordinates and strand: the resolver's junction map and the scan's strand table key.
struct SJKey {
    int32_t ref_id;
    int32_t start;
    int32_t end;
    int32_t strand;

    bool operator==(const SJKey& o) const {
        return ref_id == o.ref_id && start == o.start && end == o.end && strand == o.strand;
    }
    bool operator<(const SJKey& o) const {
        return std::tie(ref_id, start, end, strand) < std::tie(o.ref_id, o.start, o.end, o.strand);
    }
};

struct SJKeyHash {
    size_t operator()(const SJKey& k) const {
        size_t h = std::hash<int32_t>()(k.ref_id);
        hash_combine(h, k.start);
        hash_combine(h, k.end);
        hash_combine(h, k.strand);
        return h;
    }
};

}  // namespace rigel
