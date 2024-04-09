#pragma once

#include "fibonacci_mesh.cpp"
#include "utils.cpp"
#include "detectors.h"

#include <vector>
#include <complex>
#include <cmath> // For std::isnan and std::pow

namespace DETECTOR {

    using complex128 = std::complex<double>;

    template <class T>
    double Detector::get_coupling_point_no_coherent(T &scatterer)
    {
        auto [theta_field, phi_field] = scatterer.compute_unstructured_fields(this->state.fibonacci_mesh);

        double
            coupling_theta = this->get_norm2_squared(theta_field),
            coupling_phi = this->get_norm2_squared(phi_field);

        this->apply_polarization_filter(
            coupling_theta,
            coupling_phi,
            state.polarization_filter
        );

        return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * this->state.fibonacci_mesh.dOmega;
    }



    template <class T>
    double Detector::get_coupling_mean_no_coherent(T &scatterer)
    {
        return get_coupling_point_no_coherent(scatterer);
    }

    template <class T>
    double Detector::get_coupling_point_coherent(T &scatterer)
    {
        auto [theta_field, phi_field] = scatterer.compute_unstructured_fields(this->state.fibonacci_mesh);

        auto [horizontal_projection, vertical_projection] = this->get_projected_fields(theta_field, phi_field);

        this->apply_scalar_field(horizontal_projection, vertical_projection);

        double
            coupling_theta = get_norm1_squared(horizontal_projection),
            coupling_phi = get_norm1_squared(vertical_projection);

        this->apply_polarization_filter(
            coupling_theta,
            coupling_phi,
            state.polarization_filter
        );

        return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * this->state.fibonacci_mesh.dOmega;
    }


    template <class T> double
    Detector::get_coupling_mean_coherent(T &scatterer)
    {
        auto [theta_field, phi_field] = scatterer.compute_unstructured_fields(this->state.fibonacci_mesh);

        auto [horizontal_projection, vertical_projection] = this->get_projected_fields(theta_field, phi_field);

        this->apply_scalar_field(horizontal_projection, vertical_projection);

        double
            coupling_theta = this->get_norm2_squared(horizontal_projection),
            coupling_phi = this->get_norm2_squared(vertical_projection);

        this->apply_polarization_filter(
            coupling_theta,
            coupling_phi,
            state.polarization_filter
        );

        return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * this->state.fibonacci_mesh.dOmega / this->state.fibonacci_mesh.Omega;
    }

    std::tuple<std::vector<complex128>, std::vector<complex128>>
    Detector::get_projected_fields(const std::vector<complex128> &theta_field, const std::vector<complex128> &phi_field) const
    {

        std::vector<complex128>
            horizontal_projection(theta_field.size()),
            vertical_projection(theta_field.size());

        for (size_t i=0; i<theta_field.size(); ++i)
        {
            vertical_projection[i] =
                theta_field[i] * this->state.fibonacci_mesh.vertical_perpendicular_projection[i] +
                phi_field[i] * this->state.fibonacci_mesh.vertical_parallel_projection[i] ;  // new_version

            horizontal_projection[i] =
                theta_field[i] * this->state.fibonacci_mesh.horizontal_perpendicular_projection[i] +
                phi_field[i] * this->state.fibonacci_mesh.horizontal_parallel_projection[i] ; // new_version
        }

        return std::make_tuple(horizontal_projection, vertical_projection);

    }


    void Detector::apply_scalar_field(std::vector<complex128> &field0, std::vector<complex128> &field1) const //Theta = Para
    {
        auto buffer_scalar_field = state.scalar_field.request();

        complex128 *scalar_field_ptr = (complex128 *) buffer_scalar_field.ptr;

        for (size_t i=0; i<field0.size(); i++)
        {
            field0[i] *= scalar_field_ptr[i];
            field1[i] *= scalar_field_ptr[i];
        }
    }

    template <class T> inline
    double Detector::get_norm1_squared(const std::vector<T> &array) const
    {
        T sum  = 0.0;

        for (auto v : array)
            sum += v;

        return pow(abs(sum), 2);
    }

    template <class T> inline
    double Detector::get_norm2_squared(const std::vector<T> &array) const
    {
      T sum  = 0.0;

      for (auto v : array)
        sum += pow( abs(v), 2 );

      return abs(sum);
    }


    template <class T> inline
    void Detector::square_array(std::vector<T> &array)
    {
      for (T &v : array)
         v = pow( abs(v), 2);
    }

    template <typename T> inline
    void Detector::apply_polarization_filter(T &coupling_theta, T &coupling_phi, double polarization_filter) const
    {

        if (isnan(state.polarization_filter))
            return;

        double
            theta_polarization_filtering  = pow( sin(polarization_filter), 2 ),
            phi_polarization_filtering    = pow( cos(polarization_filter), 2 );

        coupling_theta *= theta_polarization_filtering;
        coupling_phi *= phi_polarization_filtering;
    }

} // namespace DETECTOR


