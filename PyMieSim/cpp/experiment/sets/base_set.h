#pragma once

#include "single/scatterer/sphere/sphere.h"
#include "single/scatterer/cylinder/cylinder.h"
#include "single/scatterer/coreshell/coreshell.h"
#include "single/detector/detector.h"
#include "single/source/source.h"




#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <pybind11/complex.h> // For std::complex support
#include <limits>

#include <pint/pint.h>


// namespace py = pybind11;
// typedef std::complex<double> complex128;

// std::vector<double> cast_scalar_or_array_to_vector_double(const py::object& obj) {

//     // If already iterable (NumPy array, list, tuple)
//     if (py::isinstance<py::sequence>(obj) && !py::isinstance<py::str>(obj)) {
//         return obj.cast<std::vector<double>>();
//     }

//     // Otherwise treat as scalar
//     return { obj.cast<double>() };
// }


// std::vector<unsigned> cast_scalar_or_array_to_vector_unsigned(const py::object& obj) {

//     // If already iterable (NumPy array, list, tuple)
//     if (py::isinstance<py::sequence>(obj) && !py::isinstance<py::str>(obj)) {
//         return obj.cast<std::vector<unsigned>>();
//     }

//     // Otherwise treat as scalar
//     return { obj.cast<unsigned>() };
// }

// std::vector<std::string> cast_scalar_or_array_to_vector_string(const py::object& obj) {

//     // If already iterable (NumPy array, list, tuple)
//     if (py::isinstance<py::sequence>(obj) && !py::isinstance<py::str>(obj)) {
//         return obj.cast<std::vector<std::string>>();
//     }

//     // Otherwise treat as scalar
//     return { obj.cast<std::string>() };
// }


// Base class to reduce redundancy
class BaseSet{
    public:
        bool is_sequential;
        std::vector<size_t> shape = {1};
        size_t current_index = 0;
        size_t total_combinations = 1;
        bool is_empty = true;
        virtual void update_shape() {};

        BaseSet() = default;
        virtual ~BaseSet() = default;

        BaseSet(bool is_sequential){this->is_sequential = is_sequential;}

        size_t get_vector_sigma(const std::vector<size_t> &vector)
        {
            size_t sigma = 1;
            for (auto e: vector)
              sigma *= e;

            return sigma;
        }

        // Calculate the multi-dimensional indices for the current index
        std::vector<size_t> calculate_indices(size_t flat_index) const
        {
            std::vector<size_t> indices(shape.size());
            for (size_t i = shape.size(); i-- > 0;)
            {
                indices[i] = flat_index % shape[i];
                flat_index /= shape[i];
            }

            return indices;
        }

        template <typename Container> void check_size(const Container& container, size_t expected_size, const std::string& name) const {
            if (container.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. " + name + " has size " + std::to_string(container.size()) + ", but wavelength size is " + std::to_string(expected_size) + ".");
        }
    };
