#pragma once

#include <vector>
#include <complex>
#include <cmath> // For std::isnan and std::pow
#include <stdexcept>
#include "utils/special_function.cpp"
#include "fibonacci/fibonacci.h"
#include "utils/numpy_interface.h"
#include <scatterer/base_scatterer/base_scatterer.h>
#include "../mode_field/mode_field.h"
#include <iostream>




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

        ModeID mode_id;
        ModeField mode_field;

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

        double NA2Angle(double NA) const;

        void initialize(const double &medium_refractive_index) {
            this->parse_mode(this->mode_number);
            this->mode_field = ModeField(this->mode_id);

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

            this->scalar_field = this->mode_field.get_unstructured(
                this->fibonacci_mesh.base_cartesian_coordinates.x,
                this->fibonacci_mesh.base_cartesian_coordinates.y
            );

        }

        double get_coupling(const BaseScatterer& scatterer) const {
            if (this->coherent)
                return this->mean_coupling ? get_coupling_mean_coherent(scatterer) : get_coupling_point_coherent(scatterer);
            else
                return this->mean_coupling ? get_coupling_mean_no_coherent(scatterer) : get_coupling_point_no_coherent(scatterer);
        }

        double get_coupling_point_no_coherent(const BaseScatterer& scatterer) const;
        double get_coupling_mean_no_coherent(const BaseScatterer& scatterer) const;
        double get_coupling_point_coherent(const BaseScatterer& scatterer) const;
        double get_coupling_mean_coherent(const BaseScatterer& scatterer) const;

        [[nodiscard]] pybind11::array_t<complex128> get_structured_scalarfield(const size_t sampling) const;

    private:
        void parse_mode(const std::string& mode_number);

        std::tuple<std::vector<complex128>, std::vector<complex128>>
        get_projected_fields(const std::vector<complex128>& theta_field, const std::vector<complex128>& phi_field) const;

        void apply_scalar_field(std::vector<complex128> &field0, std::vector<complex128> &field1) const;
        template <typename T> inline double get_norm1_squared(const std::vector<T>& array) const;
        template <typename T> inline double get_norm2_squared(const std::vector<T>& array) const;
        template <typename T> inline void square_array(std::vector<T>& array);
        template <typename T> inline void apply_polarization_filter(T& coupling_theta, T& coupling_phi, double polarization_filter) const;
};



