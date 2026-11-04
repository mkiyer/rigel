// em_scatter.h — the global candidate CSR scattered into the per-locus partitions (`rigel.locus_partition`).
//
// Each locus is a list of global unit indices. The per-candidate arrays are copied segment by segment into one array
// per locus, the per-unit arrays element by element. A list element must be of the exact dtype: its data is read
// after the cast that checked it.
#pragma once

#include <cstdint>
#include <cstring>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "ndarray_util.h"

namespace nb = nanobind;

namespace rigel {

// Locus li's global unit indices.
inline Array1<int32_t> locus_units_of(const nb::list& locus_units, int li) {
    return nb::cast<Array1<int32_t>>(locus_units[li], false);
}

// Per-locus CSR offsets from the global offsets and each locus's unit list: offsets[li][0] = 0 and
// offsets[li][k+1] = offsets[li][k] + (g_offsets[u+1] - g_offsets[u]) for its k-th unit u. One int64 array per locus.
inline nb::list build_partition_offsets(Array1<int64_t> g_offsets, nb::list locus_units, int n_loci) {
    const int64_t* goff = g_offsets.data();
    nb::list out;
    for (int li = 0; li < n_loci; ++li) {
        const auto units = locus_units_of(locus_units, li);
        const int32_t* u = units.data();
        std::vector<int64_t> offsets(units.shape(0) + 1);
        offsets[0] = 0;
        for (size_t k = 0; k < units.shape(0); ++k) offsets[k + 1] = offsets[k] + (goff[u[k] + 1] - goff[u[k]]);
        out.append(vec_to_ndarray(std::move(offsets)));
    }
    return out;
}

// A per-candidate global CSR array scattered into per-locus arrays: the segment of the locus's k-th unit u,
// src[g_off[u], g_off[u+1]), is copied to dst + p_off[k]. One array per locus.
template <typename T>
nb::list scatter_candidates(Array1<T> global_arr, Array1<int64_t> g_offsets, nb::list locus_units,
                            nb::list partition_offsets, int n_loci) {
    const T* src = global_arr.data();
    const int64_t* goff = g_offsets.data();
    nb::list out;
    for (int li = 0; li < n_loci; ++li) {
        const auto units = locus_units_of(locus_units, li);
        const int64_t* p_off = nb::cast<Array1<int64_t>>(partition_offsets[li], false).data();
        const size_t n_units = units.shape(0);
        std::vector<T> dst(static_cast<size_t>(p_off[n_units]));
        for (size_t k = 0; k < n_units; ++k) {
            const int32_t u = units.data()[k];
            const int64_t seg_len = goff[u + 1] - goff[u];
            if (seg_len > 0) std::memcpy(dst.data() + p_off[k], src + goff[u], static_cast<size_t>(seg_len) * sizeof(T));
        }
        out.append(vec_to_ndarray(std::move(dst)));
    }
    return out;
}

// A per-unit global array gathered into per-locus arrays: dst[k] = src[units[k]]. One array per locus.
template <typename T>
nb::list scatter_units(Array1<T> global_arr, nb::list locus_units, int n_loci) {
    const T* src = global_arr.data();
    nb::list out;
    for (int li = 0; li < n_loci; ++li) {
        const auto units = locus_units_of(locus_units, li);
        std::vector<T> dst(units.shape(0));
        for (size_t k = 0; k < dst.size(); ++k) dst[k] = src[units.data()[k]];
        out.append(vec_to_ndarray(std::move(dst)));
    }
    return out;
}

}  // namespace rigel
