#pragma once

#include "fibonnaci_mesh.cpp"
#include "utils.cpp"

namespace DETECTOR
{
    struct State
    {
        std::vector<complex128>
            scalar_field;

        double
            NA,
            phi_offset,
            gamma_offset,
            polarization_filter,
            rotation_angle;

        bool
            coherent,
            point_coupling;

        State(){}

        State(
            const std::vector<complex128> &scalar_field,
            const double &NA,
            const double &phi_offset,
            const double &gamma_offset,
            const double &polarization_filter,
            const double &rotation_angle,
            const bool &coherent,
            const bool &point_coupling)
            : scalar_field(scalar_field),
                NA(NA),
                phi_offset(phi_offset),
                gamma_offset(gamma_offset),
                polarization_filter(polarization_filter),
                rotation_angle(rotation_angle),
                coherent(coherent),
                point_coupling(point_coupling)
            {

            }
    };


    class Detector
    {
    public:
        FibonacciMesh
            fibonacci_mesh;

        State
            state;

        bool
            coherent;

        Detector(){}


        Detector(
            const std::vector<complex128> &scalar_field,
            const double &NA,
            const double& phi_offset,
            const double &gamma_offset,
            const double &polarization_filter,
            const double &rotation_angle,
            const bool &coherent,
            const bool &point_coupling)
            : state(scalar_field, NA, phi_offset, gamma_offset, polarization_filter, rotation_angle, coherent, point_coupling)
        {
            this->fibonacci_mesh = FibonacciMesh(
                state.scalar_field.size(),
                NA2Angle(state.NA),
                state.phi_offset,
                state.gamma_offset,
                state.rotation_angle
            );
        }

    Detector(State &state) : state(state)
    {
        this->fibonacci_mesh = FibonacciMesh(
            state.scalar_field.size(),
            NA2Angle(state.NA),
            state.phi_offset,
            state.gamma_offset,
            state.rotation_angle
        );
    }


    void rotate_around_axis(double &rotation_angle)
    {
        this->fibonacci_mesh.rotate_around_axis(rotation_angle);
    }

    template <class T>
    double get_coupling(T &scatterer)
    {
        if (state.point_coupling && state.coherent)
            return this->get_coupling_point_coherent(scatterer);

        if (!state.point_coupling && state.coherent)
            return this->get_coupling_mean_coherent(scatterer);

        if (state.point_coupling && !state.coherent)
            return this->get_coupling_point_no_coherent(scatterer);

        if (!state.point_coupling && !state.coherent)
            return this->get_coupling_mean_no_coherent(scatterer);
    }

    template <class T> double get_coupling_point_no_coherent(T &scatterer)
    {
        auto [theta_field, phi_field] = scatterer.compute_unstructured_fields(this->fibonacci_mesh);

        double
            coupling_theta = this->get_norm2_squared(theta_field),
            coupling_phi = this->get_norm2_squared(phi_field);

        this->apply_polarization_filter(
            coupling_theta,
            coupling_phi,
            state.polarization_filter
        );

        return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * this->fibonacci_mesh.dOmega;
    }


    template <class T> double get_coupling_mean_no_coherent(T &scatterer)
    {
        return get_coupling_point_no_coherent(scatterer);
    }


    template <class T> double get_coupling_point_coherent(T &scatterer)
    {
        auto [theta_field, phi_field] = scatterer.compute_unstructured_fields(this->fibonacci_mesh);

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

        return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * this->fibonacci_mesh.dOmega;
    }


    template <class T> double get_coupling_mean_coherent(T &scatterer)
    {
        auto [theta_field, phi_field] = scatterer.compute_unstructured_fields(this->fibonacci_mesh);

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

        return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * this->fibonacci_mesh.dOmega / this->fibonacci_mesh.Omega;
    }


    std::tuple<std::vector<complex128>, std::vector<complex128>>
    get_projected_fields(std::vector<complex128> &theta_field, std::vector<complex128> &phi_field) const
    {

        std::vector<complex128>
            horizontal_projection(theta_field.size()),
            vertical_projection(theta_field.size());

        for (size_t i=0; i<theta_field.size(); ++i)
        {
            vertical_projection[i] =
                theta_field[i] * this->fibonacci_mesh.vertical_perpendicular_projection[i] +
                phi_field[i] * this->fibonacci_mesh.vertical_parallel_projection[i] ;  // new_version

            horizontal_projection[i] =
                theta_field[i] * this->fibonacci_mesh.horizontal_perpendicular_projection[i] +
                phi_field[i] * this->fibonacci_mesh.horizontal_parallel_projection[i] ; // new_version
        }

        return std::make_tuple(horizontal_projection, vertical_projection);

    }

    void apply_scalar_field(std::vector<complex128> &field0, std::vector<complex128> &field1) const //Theta = Para
    {
        for (size_t i=0; i<field0.size(); i++)
        {
            field0[i] *= state.scalar_field[i];
            field1[i] *= state.scalar_field[i];
        }
    }


    template <class T> inline double
    get_norm1_squared(const std::vector<T> &array) const
    {
        T sum  = 0.0;

        for (auto v : array)
            sum += v;

        return pow(abs(sum), 2);
    }

    template <class T> inline double
    get_norm2_squared(const std::vector<T> &array) const
    {
      T sum  = 0.0;

      for (auto v : array)
        sum += pow( abs(v), 2 );

      return abs(sum);
    }


    template <class T> inline void square_array(std::vector<T> &array)
    {
      for (T &v : array)
         v = pow( abs(v), 2);
    }

    template <class T> inline void
    apply_polarization_filter(
        T &coupling_theta,
        T &coupling_phi,
        const double &polarization_filter) const
    {

        if (isnan(state.polarization_filter))
            return;

        double
            theta_polarization_filtering  = pow( sin(polarization_filter), 2 ),
            phi_polarization_filtering    = pow( cos(polarization_filter), 2 );

        coupling_theta *= theta_polarization_filtering;
        coupling_phi *= phi_polarization_filtering;
    }


  };
}
