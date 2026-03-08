#pragma once

#include <vector>
#include <complex>
#include <string>
#include <stdexcept>

typedef std::complex<double> complex128;

enum class PropertyMode {
    Constant,
    Spectral
};

template<typename T>
class Properties {
public:

    PropertyMode mode;

    std::vector<T> constant_values;
    std::vector<std::vector<T>> spectral_values;

    // Only meaningful when mode == Spectral
    std::vector<std::string> material_names;

public:

    Properties() = default;

    Properties(const std::vector<T>& values)
        : mode(PropertyMode::Constant),
          constant_values(values) {}

    Properties(
        const std::vector<std::vector<T>>& values,
        const std::vector<std::string>& names
    )
        : mode(PropertyMode::Spectral),
          spectral_values(values),
          material_names(names)
    {
        if (spectral_values.size() != material_names.size()) {
            throw std::runtime_error(
                "Number of spectral property sets must match number of material names."
            );
        }
    }

    size_t size() const {

        if (mode == PropertyMode::Constant)
            return constant_values.size();

        return spectral_values.size();
    }

    T get(const size_t& index, const size_t& wl_index) const {

        if (mode == PropertyMode::Constant)
            return constant_values[index];

        return spectral_values[index][wl_index];
    }

    const std::string& get_material_name(const size_t& index) const {

        if (mode != PropertyMode::Spectral)
            throw std::runtime_error(
                "Material names are only defined for spectral properties."
            );

        return material_names[index];
    }

    bool is_constant() const {
        return mode == PropertyMode::Constant;
    }

    bool is_spectral() const {
        return mode == PropertyMode::Spectral;
    }
};

using ScattererProperties = Properties<complex128>;
using MediumProperties = Properties<double>;