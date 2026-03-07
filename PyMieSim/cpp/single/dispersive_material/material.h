#pragma once

#include <vector>
#include <complex>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cstddef>

using complex128 = std::complex<double>;

class BaseMaterial {
public:
    virtual ~BaseMaterial() = default;

    void initialize(const double wavelength) {
        this->refractive_index = this->compute_refractive_index(wavelength);
        this->is_initialized = true;
    }

    complex128 get_refractive_index() const {
        if (!this->is_initialized) {
            throw std::runtime_error(
                "Material refractive index was requested before initialization."
            );
        }

        return this->refractive_index;
    }

    complex128 get_refractive_index(const double wavelength) const {
        return this->compute_refractive_index(wavelength);
    }

public:
    complex128 refractive_index = complex128(0.0, 0.0);
    bool is_initialized = false;

    virtual complex128 compute_refractive_index(double wavelength) const = 0;
};

class ConstantMaterial : public BaseMaterial {
public:
    complex128 constant_refractive_index;

    ConstantMaterial() = default;

    explicit ConstantMaterial(const complex128& refractive_index)
        : constant_refractive_index(refractive_index) {}

protected:
    complex128 compute_refractive_index(const double wavelength) const override {
        (void)wavelength;
        return constant_refractive_index;
    }
};

class DispersiveMaterial : public BaseMaterial {
public:
    std::string name;
    std::vector<double> wavelengths;
    std::vector<complex128> refractive_indices;
    bool allow_extrapolation = false;

    DispersiveMaterial() = default;

    DispersiveMaterial(
        const std::string& name,
        const std::vector<double>& wavelengths,
        const std::vector<complex128>& refractive_indices,
        const bool allow_extrapolation = false
    )
        : name(name),
          wavelengths(wavelengths),
          refractive_indices(refractive_indices),
          allow_extrapolation(allow_extrapolation)
    {
        this->validate_input_data();
    }

protected:
    complex128 compute_refractive_index(const double wavelength) const override {
        this->validate_query_preconditions();

        const std::size_t number_of_points = wavelengths.size();

        if (number_of_points == 1) {
            return refractive_indices.front();
        }

        if (wavelength < wavelengths.front()) {
            if (!allow_extrapolation) {
                throw std::out_of_range(
                    "Wavelength is below the range of the material data."
                );
            }

            return interpolate_between_indices(wavelength, 0, 1);
        }

        if (wavelength > wavelengths.back()) {
            if (!allow_extrapolation) {
                throw std::out_of_range(
                    "Wavelength is above the range of the material data."
                );
            }

            return interpolate_between_indices(
                wavelength,
                number_of_points - 2,
                number_of_points - 1
            );
        }

        const auto lower_iterator = std::lower_bound(
            wavelengths.begin(),
            wavelengths.end(),
            wavelength
        );

        if (lower_iterator == wavelengths.begin()) {
            return refractive_indices.front();
        }

        if (lower_iterator == wavelengths.end()) {
            return refractive_indices.back();
        }

        const std::size_t upper_index = static_cast<std::size_t>(
            std::distance(wavelengths.begin(), lower_iterator)
        );

        if (*lower_iterator == wavelength) {
            return refractive_indices[upper_index];
        }

        const std::size_t lower_index = upper_index - 1;

        return interpolate_between_indices(
            wavelength,
            lower_index,
            upper_index
        );
    }

private:
    void validate_input_data() const {
        if (wavelengths.empty() || refractive_indices.empty()) {
            throw std::runtime_error("Material data is empty.");
        }

        if (wavelengths.size() != refractive_indices.size()) {
            throw std::runtime_error(
                "Wavelength and refractive index arrays must have the same size."
            );
        }

        for (std::size_t index = 1; index < wavelengths.size(); ++index) {
            if (wavelengths[index] <= wavelengths[index - 1]) {
                throw std::runtime_error(
                    "Wavelength data must be strictly increasing."
                );
            }
        }
    }

    void validate_query_preconditions() const {
        if (wavelengths.empty() || refractive_indices.empty()) {
            throw std::runtime_error("Material data is empty.");
        }

        if (wavelengths.size() != refractive_indices.size()) {
            throw std::runtime_error(
                "Wavelength and refractive index arrays must have the same size."
            );
        }
    }

    complex128 interpolate_between_indices(
        const double wavelength,
        const std::size_t lower_index,
        const std::size_t upper_index
    ) const {
        const double lower_wavelength = wavelengths[lower_index];
        const double upper_wavelength = wavelengths[upper_index];

        const complex128 lower_refractive_index = refractive_indices[lower_index];
        const complex128 upper_refractive_index = refractive_indices[upper_index];

        const double interpolation_fraction =
            (wavelength - lower_wavelength) / (upper_wavelength - lower_wavelength);

        return lower_refractive_index * (1.0 - interpolation_fraction)
             + upper_refractive_index * interpolation_fraction;
    }
};