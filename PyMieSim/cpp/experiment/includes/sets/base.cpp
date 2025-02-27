#pragma once

#include "single/includes/sphere.cpp"
#include "single/includes/cylinder.cpp"
#include "single/includes/coreshell.cpp"
#include "single/includes/detectors.cpp"
#include "single/includes/sources.cpp"

using complex128 = std::complex<double>;

// Base class to reduce redundancy
class BaseSet{
    public:
        std::vector<size_t> shape;
        size_t current_index = 0;
        size_t total_combinations = 1;
        virtual void update_shape() {};

        BaseSet() = default;
        virtual ~BaseSet() = default;

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

        template <typename Container>
        void check_size(const Container& container, size_t expected_size, const std::string& name) const {
            if (container.size() != expected_size)
                throw std::runtime_error("Error: Vector size mismatch in sequential computation. " + name + " has size " + std::to_string(container.size()) + ", but wavelength size is " + std::to_string(expected_size) + ".");
        }
};
