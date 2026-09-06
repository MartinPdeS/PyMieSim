#pragma once

#include <cstring>
#include <limits>
#include <sstream>
#include <typeinfo>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;
typedef std::complex<double> complex128;


template <typename T>
std::vector<size_t> get_stride(const std::vector<size_t>& dimensions)
{
    if (dimensions.empty())
        return {};

    std::vector<size_t> stride(dimensions.size());
    stride.back() = sizeof(T);  // Start from the innermost dimension

    for (int64_t i = dimensions.size() - 2; i >= 0; --i)
        stride[i] = stride[i + 1] * dimensions[i + 1];

    return stride;
}

/*
    @brief Creates a zero-copy NumPy array view of a std::vector.
    @tparam Owner The type of the owner object that manages the lifetime of the vector.
    @tparam dtype The data type of the elements in the vector and NumPy array.
    @param owner The owner object that manages the lifetime of the vector.
    @param vector The std::vector to be viewed as a NumPy array.
    @return A NumPy array that is a zero-copy view of the input vector.
    @note The lifetime of the returned NumPy array is tied to the lifetime of the owner object.
*/
template <class Owner, class dtype> inline pybind11::array_t<dtype> vector_as_numpy_view(
    Owner& owner,
    std::vector<dtype>& vector
) {
    const int64_t n = static_cast<int64_t>(vector.size());
    const int64_t stride = static_cast<int64_t>(sizeof(dtype));
    // tie lifetime to the owner
    return pybind11::array_t<dtype>({n}, {stride}, vector.data(), pybind11::cast(&owner));
}


/*
    @brief Assigns data from a NumPy array to a std::vector.
    @tparam dtype The data type of the elements in the vector and NumPy array.
    @param vector The std::vector to which data will be assigned.
    @param numpy_array The NumPy array from which data will be copied.
    @note This function resizes the vector to match the size of the NumPy array and copies the data from the array to the vector.
*/
template <class dtype> inline void vector_assign_from_numpy(
    std::vector<dtype>& vector,
    pybind11::array_t<dtype, pybind11::array::c_style | pybind11::array::forcecast> numpy_array
) {
    const size_t n = static_cast<size_t>(numpy_array.shape(0));
    vector.resize(n);
    std::memcpy(vector.data(), numpy_array.data(), n * sizeof(dtype));
}

/*
    @brief Copies a std::vector into an independently owned NumPy array.
    @tparam T The data type of the elements in the vector and NumPy array.
    @param data The std::vector whose data will be copied to the NumPy array.
    @param shape The desired shape of the resulting NumPy array.
    @return A NumPy array containing the data from the input vector.
    @note The input vector is unchanged; the result owns its memory.
*/
template <class T>
inline pybind11::array_t<T> vector_to_numpy_copy(
    const std::vector<T>& data,
    const std::vector<size_t>& shape
)
{
    // product(shape) in size_t with a simple overflow check
    size_t expect = 1;
    for (size_t d : shape) {
        if (d != 0 && expect > std::numeric_limits<size_t>::max() / d) {
            std::ostringstream os;
            os << "Shape product overflows size_t. shape=[";
            for (size_t i = 0; i < shape.size(); ++i) {
                os << shape[i] << (i + 1 < shape.size() ? ", " : "");
            }
            os << "]";
            throw pybind11::value_error(os.str());
        }
        expect *= d;
    }

    const size_t got = data.size();

    if (expect != got) {
        std::ostringstream os;
        os << "Shape does not match vector size. "
           << "shape=[";
        for (size_t i = 0; i < shape.size(); ++i) {
            os << shape[i] << (i + 1 < shape.size() ? ", " : "");
        }
        os << "] product=" << expect
           << " but vector size=" << got
           << " (dtype=" << typeid(T).name() << ")";
        throw pybind11::value_error(os.str());
    }

    // pybind11 wants int64_t for the constructor, cast once here
    std::vector<pybind11::ssize_t> shape_ss;
    shape_ss.reserve(shape.size());
    for (size_t d : shape) shape_ss.push_back(static_cast<pybind11::ssize_t>(d));

    pybind11::array_t<T> out(shape_ss);
    if (got > 0) {
        std::memcpy(out.mutable_data(),
                    data.data(),
                    got * sizeof(T));
    }
    return out;
}

// Convenience overload for a one-dimensional copy.
template <class T>
inline pybind11::array_t<T> vector_to_numpy_copy(const std::vector<T>& data) {
    return vector_to_numpy_copy(data, {data.size()});
}
