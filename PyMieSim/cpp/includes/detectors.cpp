#pragma once

#include "fibonnaci_mesh.cpp"
#include "utils.cpp"

namespace DETECTOR
{

    struct State
    {
        std::vector<complex128> scalar_field;

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
            std::vector<complex128> &scalar_field,
            double &NA,
            double &phi_offset,
            double &gamma_offset,
            double &polarization_filter,
            double &rotation_angle,
            bool   &coherent,
            bool   &point_coupling)
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
        FibonacciMesh Mesh;
        State state;
        bool coherent;

        Detector(){}


        Detector(
            std::vector<complex128> &scalar_field,
            double &NA,
            double& phi_offset,
            double &gamma_offset,
            double &polarization_filter,
            double &rotation_angle,
            bool &coherent,
            bool &point_coupling)
            : state(scalar_field, NA, phi_offset, gamma_offset, polarization_filter, rotation_angle, coherent, point_coupling)
        {
            this->Mesh = FibonacciMesh(
                state.scalar_field.size(),
                NA2Angle(state.NA),
                state.phi_offset,
                state.gamma_offset,
                state.rotation_angle
            );
        }

    Detector(State &state) : state(state)
    {
        this->Mesh = FibonacciMesh(
            state.scalar_field.size(),
            NA2Angle(state.NA),
            state.phi_offset,
            state.gamma_offset,
            state.rotation_angle
        );
    }


    void rotate_around_axis(double &rotation_angle)
    {
        this->Mesh.rotate_around_axis(rotation_angle);
    }

    template <class T>
    double Coupling(T &scatterer)
    {
        if (state.point_coupling && state.coherent)
            return get_coupling_point_coherent(scatterer);

        if (!state.point_coupling && state.coherent)
            return get_coupling_mean_coherent(scatterer);

        if (state.point_coupling && !state.coherent)
            return get_coupling_point_no_coherent(scatterer);

        if (!state.point_coupling && !state.coherent)
            return get_coupling_mean_no_coherent(scatterer);
    }

    template <class T> double get_coupling_point_no_coherent(T &scatterer)
    {
        auto [theta_field, phi_field] = scatterer.compute_unstructured_fields(Mesh);

        double coupling_theta = get_norm2_squared(theta_field),
               coupling_phi = get_norm2_squared(phi_field);

        apply_polarization_filter(
            coupling_theta,
            coupling_phi,
            state.polarization_filter
        );

        return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * Mesh.dOmega;
    }


    template <class T> double get_coupling_mean_no_coherent(T &scatterer)
    {
        return get_coupling_point_no_coherent(scatterer);
    }


    template <class T> double get_coupling_point_coherent(T &scatterer)
    {
        auto [theta_field, phi_field] = scatterer.compute_unstructured_fields(Mesh);

        auto [horizontal_projection, vertical_projection] = get_projected_fields(theta_field, phi_field);

        apply_scalar_field(horizontal_projection, vertical_projection);

        double coupling_theta = get_norm1_squared(horizontal_projection),
               coupling_phi = get_norm1_squared(vertical_projection);

        apply_polarization_filter(
            coupling_theta,
            coupling_phi,
            state.polarization_filter
        );

        return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * Mesh.dOmega;
    }


    template <class T> double get_coupling_mean_coherent(T &scatterer)
    {
        auto [theta_field, phi_field] = scatterer.compute_unstructured_fields(Mesh);

        auto [horizontal_projection, vertical_projection] = get_projected_fields(theta_field, phi_field);

        apply_scalar_field(horizontal_projection, vertical_projection);

        double coupling_theta = get_norm2_squared(horizontal_projection),
               coupling_phi = get_norm2_squared(vertical_projection);

        apply_polarization_filter(
            coupling_theta,
            coupling_phi,
            state.polarization_filter
        );

        return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * Mesh.dOmega / Mesh.Omega;
    }


    std::tuple<std::vector<complex128>, std::vector<complex128>>
    get_projected_fields(std::vector<complex128> &theta_field, std::vector<complex128> &phi_field)
    {

        std::vector<complex128> horizontal_projection(theta_field.size()),
                                vertical_projection(theta_field.size());

        for (size_t i=0; i<theta_field.size(); ++i)
        {
            vertical_projection[i] = theta_field[i] * Mesh.VPara[i] + phi_field[i] * Mesh.VPerp[i] ;
            horizontal_projection[i] = theta_field[i] * Mesh.HPara[i] + phi_field[i] * Mesh.HPerp[i] ;
        }

        return std::make_tuple(horizontal_projection, vertical_projection);

    }

    void apply_scalar_field(std::vector<complex128> &Field0, std::vector<complex128> &Field1) //Theta = Para
    {
        for (size_t i=0; i<Field0.size(); i++)
        {
            Field0[i] *= state.scalar_field[i];
            Field1[i] *= state.scalar_field[i];
        }
    }


    template <class T> inline double get_norm1_squared(const std::vector<T> &array)
    {
        T sum  = 0.0;

        for (auto v : array)
            sum += v;

        return pow(abs(sum), 2);
    }

    template <class T> inline double get_norm2_squared(const std::vector<T> &array)
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

    template <class T> inline void apply_polarization_filter(T &coupling_theta,
                                                             T &coupling_phi,
                                                             const double &polarization_filter)
    {

        if (isnan(state.polarization_filter))
            return;

        double theta_polarization_filtering  = pow( sin(polarization_filter), 2 ),
               phi_polarization_filtering    = pow( cos(polarization_filter), 2 );

        coupling_theta *= theta_polarization_filtering;
        coupling_phi *= phi_polarization_filtering;
    }


  };
}
