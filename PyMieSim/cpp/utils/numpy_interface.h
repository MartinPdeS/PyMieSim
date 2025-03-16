#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

template <typename T>
std::vector<size_t> get_stride(const std::vector<size_t>& dimensions)
{
    if (dimensions.empty())
        return {};

    std::vector<size_t> stride(dimensions.size());
    stride.back() = sizeof(T);  // Start from the innermost dimension

    for (ssize_t i = dimensions.size() - 2; i >= 0; --i)
        stride[i] = stride[i + 1] * dimensions[i + 1];

    return stride;
}

template<typename T>
pybind11::array_t<T> _vector_to_numpy(const std::vector<T> input_vector, std::vector<size_t> shape = {})
{
    // Set shape based on input_vector size if shape is empty
    if (shape.empty()) {
        shape.push_back(input_vector.size());
    }

    // Calculate strides
    std::vector<size_t> strides(shape.size(), sizeof(T));
    for (ssize_t i = shape.size() - 2; i >= 0; --i)
        strides[i] = strides[i + 1] * shape[i + 1];

    // Create numpy array directly from input_vector data, no copying
    pybind11::array_t<T> numpy_array(
        shape, strides, input_vector.data()
    );

    return numpy_array;
}
