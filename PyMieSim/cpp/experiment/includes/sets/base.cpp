#pragma once

#include "single/includes/sphere.cpp"
#include "single/includes/cylinder.cpp"
#include "single/includes/coreshell.cpp"
#include "single/includes/detectors.cpp"
#include "single/includes/sources.cpp"

using complex128 = std::complex<double>;

// Helper to get size from variant
template<typename T>
size_t get_variant_size(const std::variant<std::vector<T>, std::vector<std::vector<T>>>& var) {
    if (std::holds_alternative<std::vector<std::vector<T>>>(var))
        return std::get<std::vector<std::vector<T>>>(var).size();
    else
        return std::get<std::vector<T>>(var).size();
}

template<typename T>
T get_variant_value(const std::variant<std::vector<T>, std::vector<std::vector<T>>>& var, size_t idx1, size_t idx2 = 0) {
    if (std::holds_alternative<std::vector<std::vector<T>>>(var))
        return std::get<std::vector<std::vector<T>>>(var)[idx1][idx2];
    else
        return std::get<std::vector<T>>(var)[idx1];
}


// Base class to reduce redundancy
class BaseSet{
    public:
        std::vector<size_t> shape;
        size_t current_index = 0;
        size_t total_combinations = 1;
        void increment_index() { current_index = (current_index + 1) % total_combinations;}
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

        size_t flatten_multi_index(const std::vector<size_t>& multi_index, const std::vector<size_t>& dimensions)
        {
            size_t flatten_index = 0;
            size_t stride = 1;

            const size_t* multi_index_ptr = multi_index.data();
            const size_t* dimensions_ptr = dimensions.data();

            // Iterate from the last dimension to the first
            for (int i = dimensions.size() - 1; i >= 0; --i) {
                flatten_index += multi_index_ptr[i] * stride;
                stride *= dimensions_ptr[i];
            }

            return flatten_index;
        }

        template <typename T>
        T get_vector_sigma(const std::vector<T> &vector)
        {
            T sigma = 1;
            for (auto e: vector)
            sigma *= e;

            return sigma;
        }

};
