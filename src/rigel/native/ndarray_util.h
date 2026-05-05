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
#include <stdexcept>
#include <string>
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

/// Move a std::vector<T> to the heap and return a capsule-backed 2-D
/// ndarray with shape (rows, cols), row-major.  Requires
/// rows * cols == v.size().  The source vector is consumed.
template <typename T>
nb::object vec_to_ndarray2d(std::vector<T>&& v, size_t rows, size_t cols) {
    auto* heap = new std::vector<T>(std::move(v));
    if (heap->size() != rows * cols) {
        size_t actual = heap->size();
        delete heap;
        throw std::runtime_error(
            "vec_to_ndarray2d: size mismatch (expected " +
            std::to_string(rows * cols) + ", got " +
            std::to_string(actual) + ")");
    }
    nb::capsule del(heap, [](void* p) noexcept {
        delete static_cast<std::vector<T>*>(p);
    });
    return nb::ndarray<nb::numpy, T, nb::ndim<2>>(
        heap->data(), {rows, cols}, del).cast();
}

}  // namespace rigel
