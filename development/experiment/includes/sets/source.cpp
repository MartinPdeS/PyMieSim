#pragma once

#include "experiment/includes/sets/base.cpp"
#include "single/includes/sources.cpp"

using complex128 = std::complex<double>;

namespace SOURCE
{
    class Set : public BaseSet
    {
    public:
        std::vector<double> wavelength;
        std::vector<std::vector<complex128>> jones_vector;
        std::vector<double> NA;
        std::vector<double> optical_power;
        std::vector<double> amplitude;
        bool is_gaussian;

        Set() = default;

        // Constructor for Gaussian source
        Set(const std::vector<double>& wavelength,
            const std::vector<std::vector<complex128>>& jones_vector,
            const std::vector<double>& NA,
            const std::vector<double>& optical_power,
            const bool is_sequential)
            : BaseSet(is_sequential), wavelength(wavelength), jones_vector(jones_vector), NA(NA), optical_power(optical_power), is_gaussian(true)
        {
            shape = {wavelength.size(), jones_vector.size(), NA.size(), optical_power.size()};

            total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);

        }

        // Constructor for Planewave source
        Set(const std::vector<double>& wavelength,
            const std::vector<std::vector<complex128>>& jones_vector,
            const std::vector<double>& amplitude,
            const bool is_sequential)
            : BaseSet(is_sequential), wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude), is_gaussian(false)
        {
            this->shape = {wavelength.size(), jones_vector.size(), amplitude.size()};
            total_combinations = is_sequential ? shape[0] : get_vector_sigma(shape);
        }

        BaseSource get_source_by_index(const size_t flat_index) const {
            std::vector<size_t> indices = calculate_indices(flat_index);

            BaseSource source;

            if (is_gaussian)
                source = Gaussian(
                    this->wavelength[indices[0]],
                    this->jones_vector[indices[1]],
                    this->NA[indices[2]],
                    this->optical_power[indices[3]]
                );
            else
                source = Planewave(
                    this->wavelength[indices[0]],
                    this->jones_vector[indices[1]],
                    this->amplitude[indices[2]]
                );

            source.indices = indices;
            source.wavelength_index = indices[0];
            return source;
        }


        void validate_sequential_data(const size_t expected_size) const {
            // Check each vector's size and throw an error with the specific vector name if sizes don't match
            this->check_size(this->wavelength, expected_size, "wavelength");
            this->check_size(this->jones_vector, expected_size, "jones_vector");

            if (is_gaussian) {
                this->check_size(this->NA, expected_size, "NA");
                this->check_size(this->optical_power, expected_size, "optical_power");
            }
            else
                this->check_size(this->amplitude, expected_size, "amplitude");
        }


        BaseSource get_source_by_index_sequential(const size_t index) const {
            if (is_gaussian)
                return  Gaussian(
                    this->wavelength[index],
                    this->jones_vector[index],
                    this->NA[index],
                    this->optical_power[index]
                );
            else
                return Planewave(
                    this->wavelength[index],
                    this->jones_vector[index],
                    this->amplitude[index]
                );

        }
    };
}
