#pragma once

#include <vector>
#include <complex>
#include <cmath> // For std::isnan and std::pow
#include "utils/special_function.cpp"
#include "single/includes/fibonacci_mesh.cpp"
#include "utils/utils.cpp"
#include "utils/numpy_interface.cpp"
#include <single/headers/LG_modes.h>
#include <single/headers/HG_modes.h>
#include <single/headers/LP_modes.h>
#include <stdexcept>


namespace DETECTOR {

    using complex128 = std::complex<double>;

    class Detector {
        public:
            std::string mode_number;
            size_t sampling = 0;
            double NA = 0.0;
            double cache_NA = 0.0;
            double phi_offset = 0.0;
            double gamma_offset = 0.0;
            double polarization_filter = 0.0;
            double rotation = 0.0;
            bool coherent = true;
            bool mean_coupling = true;
            double medium_refractive_index;
            double max_angle = 0;
            double min_angle = 0;
            std::vector<complex128> scalar_field;
            FibonacciMesh fibonacci_mesh;
            std::vector<size_t> indices;

            Detector() = default;

            Detector(
                const std::string &mode_number,
                const size_t sampling,
                const double NA,
                const double cache_NA,
                const double phi_offset,
                const double gamma_offset,
                const double polarization_filter,
                const double rotation,
                const bool coherent,
                const bool mean_coupling,
                const double medium_refractive_index = 1.0)
            : mode_number(mode_number), sampling(sampling), NA(NA), cache_NA(cache_NA), phi_offset(phi_offset), gamma_offset(gamma_offset), polarization_filter(polarization_filter),
                rotation(rotation), coherent(coherent), mean_coupling(mean_coupling), medium_refractive_index(medium_refractive_index)
            {
                this->initialize(medium_refractive_index);
            }

            void initialize(const double &medium_refractive_index) {
                this->max_angle = NA2Angle(this->NA / medium_refractive_index);
                this->min_angle = NA2Angle(this->cache_NA / medium_refractive_index);

                if (this->max_angle < this->min_angle)
                    throw std::invalid_argument("Cache NA cannot be larger than detector NA.");


                this->fibonacci_mesh = FibonacciMesh(
                    this->sampling,
                    this->max_angle,
                    this->min_angle,
                    this->phi_offset,
                    this->gamma_offset,
                    this->rotation
                );

                int number_0 = mode_number[2] - '0';
                int number_1 = mode_number[3] - '0';

                std::string mode_prefix = mode_number.substr(0, 2);

                if (mode_prefix == "LG")
                    this->scalar_field = get_LG_mode_field(this->fibonacci_mesh.base_cartesian_coordinates.x, this->fibonacci_mesh.base_cartesian_coordinates.y, number_0, number_1);
                else if (mode_prefix == "HG")
                    this->scalar_field = get_HG_mode_field(this->fibonacci_mesh.base_cartesian_coordinates.x, this->fibonacci_mesh.base_cartesian_coordinates.y, number_0, number_1);
                else if (mode_prefix == "LP")
                    this->scalar_field = get_LP_mode_field(this->fibonacci_mesh.base_cartesian_coordinates.x, this->fibonacci_mesh.base_cartesian_coordinates.y, number_0, number_1);
                else if (mode_prefix == "NC")
                    this->scalar_field = std::vector<complex128>(this->fibonacci_mesh.base_cartesian_coordinates.x.size(), 1.0);
                else
                    throw std::invalid_argument("Invalid mode family name");
            }

            template <typename T>
            double get_coupling(T& scatterer) {
                if (this->coherent)
                    return this->mean_coupling ? get_coupling_mean_coherent(scatterer) : get_coupling_point_coherent(scatterer);
                else
                    return this->mean_coupling ? get_coupling_mean_no_coherent(scatterer) : get_coupling_point_no_coherent(scatterer);
            }

            template <typename T> double get_coupling_point_no_coherent(T& scatterer);
            template <typename T> double get_coupling_mean_no_coherent(T& scatterer);
            template <typename T> double get_coupling_point_coherent(T& scatterer);
            template <typename T> double get_coupling_mean_coherent(T& scatterer);

        private:
            template <typename T> double calculate_coupling(const T& scatterer, bool point, bool coherent);
            std::tuple<std::vector<complex128>, std::vector<complex128>> get_projected_fields(const std::vector<complex128>& theta_field, const std::vector<complex128>& phi_field) const;
            void apply_scalar_field(std::vector<complex128> &field0, std::vector<complex128> &field1) const;
            template <typename T> inline double get_norm1_squared(const std::vector<T>& array) const;
            template <typename T> inline double get_norm2_squared(const std::vector<T>& array) const;
            template <typename T> inline void square_array(std::vector<T>& array);
            template <typename T> inline void apply_polarization_filter(T& coupling_theta, T& coupling_phi, double polarization_filter) const;
    };
} // namespace DETECTOR


