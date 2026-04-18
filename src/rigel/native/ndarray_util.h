/**
 * ndarray_util.h — Zero-copy std::vector<T> → numpy ndarray transfer.
 *
 * Moves a vector to heap-allocated storage and returns a capsule-backed
 * 1-D numpy array.  The capsule destructor owns the memory; the vector
 * is consumed (left moved-from) after the call.
 *
 * Two overloads:
 *   vec_to_ndarray(std::vector<T>&& v)  — move from rvalue
 *   vec_to_ndarray(std::vector<T>* v)   — take ownership of heap pointer
 */

#pragma once

#include <cstddef>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

namespace rigel {

/// Move a std::vector<T> to the heap and return a capsule-backed 1-D ndarray.
/// The source vector is consumed (moved-from) after this call.
template <typename T>
nb::object vec_to_ndarray(std::vector<T>&& v) {
    auto* heap = new std::vector<T>(std::move(v));
    size_t n = heap->size();
    nb::capsule del(heap, [](void* p) noexcept {
        delete static_cast<std::vector<T>*>(p);
    });
    return nb::ndarray<nb::numpy, T, nb::ndim<1>>(
        heap->data(), {n}, del).cast();
}

/// Take ownership of a heap-allocated std::vector<T>* and return a
/// capsule-backed 1-D ndarray.  The pointer must not be used after this call.
template <typename T>
nb::object vec_to_ndarray(std::vector<T>* v) {
    size_t n = v->size();
    nb::capsule del(v, [](void* p) noexcept {
        delete static_cast<std::vector<T>*>(p);
    });
    return nb::ndarray<nb::numpy, T, nb::ndim<1>>(
        v->data(), {n}, del).cast();
}

}  // namespace rigel
