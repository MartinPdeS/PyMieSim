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
            const std::vector<double>& optical_power)
            : wavelength(wavelength), jones_vector(jones_vector), NA(NA), optical_power(optical_power), is_gaussian(true)
        {
            shape = {wavelength.size(), jones_vector.size(), NA.size(), optical_power.size()};
            total_combinations = get_vector_sigma(shape);
        }

        // Constructor for Planewave source
        Set(const std::vector<double>& wavelength,
            const std::vector<std::vector<complex128>>& jones_vector,
            const std::vector<double>& amplitude)
            : wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude), is_gaussian(false)
        {
            this->shape = {wavelength.size(), jones_vector.size(), amplitude.size()};
            total_combinations = get_vector_sigma(shape);
        }

        BaseSource get_source_by_index(size_t flat_index) const {
            std::vector<size_t> indices = calculate_indices(flat_index);

            if (is_gaussian)
                return to_gaussian(indices);
            else
                return to_planewave(indices);
        }

    private:
        Planewave to_planewave(const std::vector<size_t>& indices) const
        {
            Planewave source(
                this->wavelength[indices[0]],
                this->jones_vector[indices[1]],
                this->amplitude[indices[2]]
            );

            source.indices = indices;

            return source;
        }

        Gaussian to_gaussian(const std::vector<size_t>& indices) const
        {
            Gaussian source(
                this->wavelength[indices[0]],
                this->jones_vector[indices[1]],
                this->NA[indices[2]],
                this->optical_power[indices[3]]
            );

            source.indices = indices;

            return source;
        }
    };
}
