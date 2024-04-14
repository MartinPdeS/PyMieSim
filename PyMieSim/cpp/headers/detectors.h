#pragma once

#include <vector>
#include <complex>
#include <cmath> // For std::isnan and std::pow
#include "special_function.cpp"
#include "fibonacci_mesh.cpp"
#include "utils.cpp"
#include "numpy_interface.cpp"

namespace DETECTOR {

    using complex128 = std::complex<double>;

    class Set
    {
        public:
            std::vector<std::vector<complex128>> scalar_fields;
            std::vector<double> NA;
            std::vector<double> phi_offset;
            std::vector<double> gamma_offset;
            std::vector<double> polarization_filter;
            std::vector<double> rotation_angle;

            bool coherent;
            bool point_coupling;

            std::vector<size_t> shape;

            Set() = default;

            Set(const std::vector<std::vector<complex128>> &scalar_fields,
                const std::vector<double> &NA,
                const std::vector<double> &phi_offset,
                const std::vector<double> &gamma_offset,
                const std::vector<double> &polarization_filter,
                const std::vector<double> &rotation_angle,
                const bool &coherent,
                const bool &point_coupling)
            : scalar_fields(scalar_fields), NA(NA), phi_offset(phi_offset), gamma_offset(gamma_offset),
              polarization_filter(polarization_filter), rotation_angle(rotation_angle), coherent(coherent), point_coupling(point_coupling)
              {
                this->shape = {(size_t) this->scalar_fields.size(), this->NA.size(), this->phi_offset.size(), this->gamma_offset.size(), this->polarization_filter.size()};
              }
    };

    class Detector {
        public:
            double NA = 0.0;
            double phi_offset = 0.0;
            double gamma_offset = 0.0;
            double polarization_filter = 0.0;
            double rotation_angle = 0.0;
            bool coherent = true;
            bool point_coupling = true;
            size_t sampling = 0;
            double max_angle = 0;
            std::vector<complex128> scalar_field;
            FibonacciMesh fibonacci_mesh;

            Detector() = default;

            Detector(
                const std::vector<complex128>& scalar_field, double NA, double phi_offset,
                double gamma_offset, double polarization_filter, double rotation_angle, bool coherent, bool point_coupling
            ) : scalar_field(scalar_field), NA(NA), phi_offset(phi_offset), gamma_offset(gamma_offset), polarization_filter(polarization_filter), 
            rotation_angle(rotation_angle), coherent(coherent), point_coupling(point_coupling)
            {
                this->max_angle = NA2Angle(this->NA);
                this->sampling = scalar_field.size();
                this->scalar_field = scalar_field;

                this->fibonacci_mesh = FibonacciMesh(
                    this->sampling,
                    this->max_angle,
                    this->phi_offset,
                    this->gamma_offset,
                    this->rotation_angle
                );
            }

            Detector(
                size_t sampling, double NA, double phi_offset,
                double gamma_offset, double polarization_filter, double rotation_angle, bool coherent, bool point_coupling
            ) : sampling(sampling), NA(NA), phi_offset(phi_offset), gamma_offset(gamma_offset), polarization_filter(polarization_filter), 
            rotation_angle(rotation_angle), coherent(coherent), point_coupling(point_coupling)
            {
                this->max_angle = NA2Angle(this->NA);

                this->fibonacci_mesh = FibonacciMesh(
                    this->sampling,
                    this->max_angle,
                    this->phi_offset,
                    this->gamma_offset,
                    this->rotation_angle
                );
            }

            template <typename T>
            double get_coupling(T& scatterer) {
                if (this->coherent) {
                    return this->point_coupling ? get_coupling_point_coherent(scatterer) : get_coupling_mean_coherent(scatterer);
                } else {
                    return this->point_coupling ? get_coupling_point_no_coherent(scatterer) : get_coupling_mean_no_coherent(scatterer);
                }
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


