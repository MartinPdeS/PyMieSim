#pragma once

#include <vector>

class Material {
    public:
        std::vector<double> refractive_index_spectrum;

        Material(const std::vector<double>& _refractive_index_spectrum)
            : refractive_index_spectrum(_refractive_index_spectrum) {}

        Material() = default;
        ~Material() = default;
};