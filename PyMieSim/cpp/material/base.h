#pragma once

#include <vector>
#include <complex>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cstddef>

using complex128 = std::complex<double>;


template<typename RefractiveIndexType>
class Base {
public:
    virtual ~Base() = default;

    void initialize(const double wavelength) {
        this->refractive_index = this->compute_refractive_index(wavelength);
        this->is_initialized = true;
    }

    RefractiveIndexType get_refractive_index() const {
        if (!this->is_initialized) {
            throw std::runtime_error(
                "Material refractive index was requested before initialization."
            );
        }

        return this->refractive_index;
    }

    RefractiveIndexType get_refractive_index(const double wavelength) const {
        return this->compute_refractive_index(wavelength);
    }

    virtual std::shared_ptr<Base<RefractiveIndexType>> clone() const = 0;

public:
    RefractiveIndexType refractive_index;
    bool is_initialized = false;

    virtual RefractiveIndexType compute_refractive_index(double wavelength) const = 0;
};



using BaseMaterial = Base<complex128>;
using BaseMedium = Base<double>;