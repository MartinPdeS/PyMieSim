#pragma once

#include "./base.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>


template <typename RefractiveIndexType>
class BaseSellmeier : public Base<RefractiveIndexType> {
public:
    std::string name;
    std::vector<double> coefficients;
    std::vector<double> wavelength_bound;
    bool has_wavelength_bound = false;
    bool allow_extrapolation = false;
    std::size_t formula_type = 0;

    BaseSellmeier() = default;

    BaseSellmeier(
        const std::string& name,
        const std::vector<double>& coefficients,
        const std::size_t formula_type,
        const std::vector<double>& wavelength_bound = {},
        const bool allow_extrapolation = false
    )
        : name(name),
          coefficients(coefficients),
          wavelength_bound(wavelength_bound),
          has_wavelength_bound(!wavelength_bound.empty()),
          allow_extrapolation(allow_extrapolation),
          formula_type(formula_type)
    {
        this->validate_input_data();
    }

    std::shared_ptr<Base<RefractiveIndexType>> clone() const override {
        return std::make_shared<BaseSellmeier<RefractiveIndexType>>(
            this->name,
            this->coefficients,
            this->formula_type,
            this->wavelength_bound,
            this->allow_extrapolation
        );
    }

protected:
    RefractiveIndexType compute_refractive_index(const double wavelength) const override {
        this->validate_input_data();
        this->check_wavelength(wavelength);

        if (wavelength <= 0.0) {
            throw std::runtime_error("Wavelength must be strictly positive.");
        }

        const double wavelength_in_micrometer = wavelength * 1e6;
        const double wavelength_squared = wavelength_in_micrometer * wavelength_in_micrometer;

        switch (this->formula_type) {
            case 1:
                return static_cast<RefractiveIndexType>(
                    this->compute_formula_1(wavelength_squared)
                );

            case 2:
                return static_cast<RefractiveIndexType>(
                    this->compute_formula_2(wavelength_squared)
                );

            case 5:
                return static_cast<RefractiveIndexType>(
                    this->compute_formula_5(wavelength_in_micrometer)
                );

            case 6:
                return static_cast<RefractiveIndexType>(
                    this->compute_formula_6(wavelength_in_micrometer)
                );

            default:
                throw std::runtime_error("Unsupported Sellmeier formula type.");
        }
    }

private:
    void validate_input_data() const {
        if (coefficients.empty()) {
            throw std::runtime_error("Sellmeier coefficient array is empty.");
        }

        if (formula_type != 1 && formula_type != 2 && formula_type != 5 && formula_type != 6) {
            throw std::runtime_error("Unsupported Sellmeier formula type.");
        }

        if (has_wavelength_bound && wavelength_bound.size() != 2) {
            throw std::runtime_error(
                "Wavelength bound must contain exactly two values: lower and upper bound."
            );
        }
    }

    void check_wavelength(const double wavelength) const {
        if (!has_wavelength_bound || allow_extrapolation) {
            return;
        }

        const double wavelength_in_micrometer = wavelength * 1e6;

        if (wavelength_in_micrometer < wavelength_bound[0] || wavelength_in_micrometer > wavelength_bound[1]) {
            throw std::out_of_range("Wavelength is outside the allowed range.");
        }
    }

    double compute_formula_1(const double wavelength_squared) const {
        double refractive_index_squared = 1.0;

        for (std::size_t coefficient_index = 1; coefficient_index + 1 < coefficients.size(); coefficient_index += 2) {
            const double numerator_coefficient = coefficients[coefficient_index];
            const double denominator_coefficient = coefficients[coefficient_index + 1];
            const double denominator = wavelength_squared - denominator_coefficient * denominator_coefficient;

            if (std::abs(denominator) < 1e-30) {
                throw std::runtime_error(
                    "Sellmeier Formula 1 denominator is too close to zero."
                );
            }

            refractive_index_squared +=
                (numerator_coefficient * wavelength_squared) / denominator;
        }

        if (refractive_index_squared < 0.0) {
            throw std::runtime_error(
                "Sellmeier Formula 1 produced a negative refractive_index_squared value."
            );
        }

        return std::sqrt(refractive_index_squared);
    }

    double compute_formula_2(const double wavelength_squared) const {
        double refractive_index_squared = 1.0 + coefficients[0];

        for (std::size_t coefficient_index = 1; coefficient_index + 1 < coefficients.size(); coefficient_index += 2) {
            const double numerator_coefficient = coefficients[coefficient_index];
            const double denominator_coefficient = coefficients[coefficient_index + 1];
            const double denominator = wavelength_squared - denominator_coefficient;

            if (std::abs(denominator) < 1e-30) {
                throw std::runtime_error(
                    "Sellmeier Formula 2 denominator is too close to zero."
                );
            }

            refractive_index_squared +=
                (numerator_coefficient * wavelength_squared) / denominator;
        }

        if (refractive_index_squared < 0.0) {
            throw std::runtime_error(
                "Sellmeier Formula 2 produced a negative refractive_index_squared value."
            );
        }

        return std::sqrt(refractive_index_squared);
    }

    double compute_formula_5(const double wavelength_in_micrometer) const {
        double refractive_index = 1.0 + coefficients[0];

        for (std::size_t coefficient_index = 1; coefficient_index + 1 < coefficients.size(); coefficient_index += 2) {
            const double amplitude_coefficient = coefficients[coefficient_index];
            const double exponent_coefficient = coefficients[coefficient_index + 1];

            refractive_index += amplitude_coefficient * std::pow(wavelength_in_micrometer, exponent_coefficient);
        }

        return refractive_index;
    }

    double compute_formula_6(const double wavelength_in_micrometer) const {
        double refractive_index = 1.0 + coefficients[0];

        for (std::size_t coefficient_index = 1; coefficient_index + 1 < coefficients.size(); coefficient_index += 2) {
            const double numerator_coefficient = coefficients[coefficient_index];
            const double denominator_coefficient = coefficients[coefficient_index + 1];

            const double inverse_wavelength_squared =
                1.0 / (wavelength_in_micrometer * wavelength_in_micrometer);

            const double denominator = denominator_coefficient - inverse_wavelength_squared;

            if (std::abs(denominator) < 1e-30) {
                throw std::runtime_error(
                    "Sellmeier Formula 6 denominator is too close to zero."
                );
            }

            refractive_index = numerator_coefficient / denominator;
        }

        return refractive_index;
    }
};


using SellmeierMedium = BaseSellmeier<double>;
using SellmeierMaterial = BaseSellmeier<complex128>;