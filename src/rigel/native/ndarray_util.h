// ndarray_util.h — a std::vector handed to Python as a numpy array without a copy.
//
// The vector is moved to the heap and owned by the array's capsule; the source is left moved-from.
#pragma once

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

namespace rigel {

// A read-only, C-contiguous 1-D array argument.
template <typename T>
using Array1 = nb::ndarray<const T, nb::ndim<1>, nb::c_contig>;

// A vector as a numpy array of `shape`, whose sizes multiply to the vector's length; `Elem` is the array's element
// type, of T's size.
template <typename Elem, typename T>
nb::object owned_ndarray(std::vector<T>&& v, std::initializer_list<size_t> shape) {
    static_assert(sizeof(Elem) == sizeof(T), "the array's elements are the vector's");
    size_t n = 1;
    for (size_t d : shape) n *= d;
    if (n != v.size()) throw std::logic_error("owned_ndarray: the shape does not hold the vector");
    auto* heap = new std::vector<T>(std::move(v));
    nb::capsule owner(heap, [](void* p) noexcept { delete static_cast<std::vector<T>*>(p); });
    return nb::ndarray<nb::numpy, Elem>(reinterpret_cast<Elem*>(heap->data()), shape, owner).cast();
}

// A vector as a numpy array of `shape`.
template <typename T>
nb::object vec_to_ndarray(std::vector<T>&& v, std::initializer_list<size_t> shape) {
    return owned_ndarray<T>(std::move(v), shape);
}

// A vector as a 1-D numpy array.
template <typename T>
nb::object vec_to_ndarray(std::vector<T>&& v) {
    const size_t n = v.size();
    return owned_ndarray<T>(std::move(v), {n});
}

// A byte vector as a numpy bool array of `shape`.
inline nb::object bool_ndarray(std::vector<uint8_t>&& v, std::initializer_list<size_t> shape) {
    return owned_ndarray<bool>(std::move(v), shape);
}

// A byte vector as a 1-D numpy bool array.
inline nb::object bool_ndarray(std::vector<uint8_t>&& v) {
    const size_t n = v.size();
    return owned_ndarray<bool>(std::move(v), {n});
}

}  // namespace rigel
