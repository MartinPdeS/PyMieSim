#pragma once

#include "definitions.cpp"
#include "utils.cpp"
#include "sources.cpp"
#include "fibonacci_mesh.cpp"
#include "VSH.cpp"
#include <cmath>
#include <vector>
#include <complex>

class BaseSphericalScatterer
{
public:
    size_t max_order;
    double size_parameter;
    double area;
    double medium_index;

    std::vector<complex128> an;
    std::vector<complex128> bn;
    std::vector<complex128> cn;
    std::vector<complex128> dn;

    SOURCE::BaseSource source;

    BaseSphericalScatterer() = default;

    BaseSphericalScatterer(const double wavelength, const std::vector<complex128> jones_vector, const double amplitude, const double medium_index)
    : source(wavelength, jones_vector, amplitude), medium_index(medium_index){}

    BaseSphericalScatterer(const SOURCE::BaseSource &source, const double medium_index) : source(source), medium_index(medium_index){}

    std::vector<complex128> get_an() const { return an; };
    std::vector<complex128> get_bn() const { return bn; };
    std::vector<complex128> get_cn() const { return cn; };
    std::vector<complex128> get_dn() const { return dn; };

    pybind11::array_t<complex128> get_an_py() { return vector_to_numpy(an, {max_order}); }
    pybind11::array_t<complex128> get_bn_py() { return vector_to_numpy(bn, {max_order}); }
    pybind11::array_t<complex128> get_cn_py() { return vector_to_numpy(cn, {max_order}); }
    pybind11::array_t<complex128> get_dn_py() { return vector_to_numpy(dn, {max_order}); }

    double get_Qforward() const {return get_Qsca() - get_Qback();};
    double get_Qpr() const {return get_Qext() - get_g() * get_Qsca();};
    double get_Qratio() const {return get_Qback() / get_Qsca();};
    double get_Qabs() const {return get_Qext() - get_Qsca();};
    double get_Csca() const {return get_Qsca() * area;};
    double get_Cext() const {return get_Qext() * area;};
    double get_Cabs() const {return get_Qabs() * area;};
    double get_Cback() const {return get_Qback() * area;};
    double get_Cforward() const {return get_Qforward() * area;};
    double get_Cpr() const {return get_Qpr() * area;};
    double get_Cratio() const {return get_Qratio() * area;};

    size_t get_wiscombe_criterion(const double size_parameter) const {
        return static_cast<size_t>(2 + size_parameter + 4 * std::cbrt(size_parameter)) + 16;
    }

    std::vector<double> get_prefactor() const {
        std::vector<double> output;

        output.reserve(max_order);

        for (size_t m = 0; m < max_order ; ++m)
            output[m] = (double) ( 2 * (m+1) + 1 ) / ( (m+1) * ( (m+1) + 1 ) );

        return output;
    }

    complex128 compute_propagator(const double &radius) const {
        return source.amplitude / (source.k * radius) * exp(-JJ * source.k * radius);
    }

    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_structured_fields(const std::vector<complex128>& S1, const std::vector<complex128>& S2, const std::vector<double>& theta, const double& radius) const {
        std::vector<complex128> phi_field, theta_field;

        size_t full_size = theta.size() * S1.size();

        phi_field.reserve(full_size);
        theta_field.reserve(full_size);

        complex128 propagator = this->compute_propagator(radius);

        for (unsigned int p=0; p < S1.size(); p++ )
            for (unsigned int t=0; t < theta.size(); t++ )
            {
                complex128 phi_point_field = propagator * S1[p] * (source.jones_vector[0] * cos(theta[t]) + source.jones_vector[1] * sin(theta[t]));
                complex128 thetea_point_field = propagator * S2[p] * (source.jones_vector[0] * sin(theta[t]) - source.jones_vector[1] * cos(theta[t]));

                phi_field.push_back(phi_point_field);
                theta_field.push_back(thetea_point_field);
            }

        return std::make_tuple(phi_field, theta_field);
    }


    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_unstructured_fields(const std::vector<double>& phi, const std::vector<double>& theta, const double radius=1.0) const
    {
        auto [S1, S2] = this->compute_s1s2(phi);

        std::vector<complex128> phi_field, theta_field;

        size_t full_size = theta.size();

        phi_field.reserve(full_size);
        theta_field.reserve(full_size);

        complex128 propagator = this->compute_propagator(radius);

        for (unsigned int idx=0; idx < full_size; idx++)
        {
            complex128 phi_field_point = propagator * S1[idx] * (source.jones_vector[0] * cos(theta[idx]) + source.jones_vector[1] * sin(theta[idx])),
            theta_field_point = propagator * S2[idx] * (source.jones_vector[0] * sin(theta[idx]) - source.jones_vector[1] * cos(theta[idx]));

            phi_field.push_back(phi_field_point);
            theta_field.push_back(theta_field_point);
        }

        return std::make_tuple(phi_field, theta_field);
    }

    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_unstructured_fields(const FibonacciMesh& fibonacci_mesh, const double radius=1.0) const
    {
        return this->compute_unstructured_fields(
            fibonacci_mesh.spherical_coordinates.Phi,
            fibonacci_mesh.spherical_coordinates.Theta,
            radius
        );
    }


    std::tuple<std::vector<complex128>, std::vector<complex128>, std::vector<double>, std::vector<double>>
    compute_full_structured_fields(const size_t& sampling, const double& radius) const
    {
        FullSteradian full_mesh = FullSteradian(sampling);

        auto [S1, S2] = this->compute_s1s2(full_mesh.spherical_coordinates.Phi);

        auto [phi_field, theta_field] = this->compute_structured_fields(
            S1,
            S2,
            full_mesh.spherical_coordinates.Theta,
            radius
        );

        return std::make_tuple(
            phi_field,
            theta_field,
            full_mesh.spherical_coordinates.Phi,
            full_mesh.spherical_coordinates.Theta
        );
    }


    std::tuple<std::vector<complex128>, FullSteradian>
    compute_full_structured_spf(const size_t &sampling, const double &radius = 1.0) const
    {
        FullSteradian full_mesh = FullSteradian(sampling);

        auto [phi_field, theta_field] = this->compute_structured_fields(
            full_mesh.spherical_coordinates.Phi,
            full_mesh.spherical_coordinates.Theta,
            radius
        );

        Squared(phi_field);
        Squared(theta_field);

        std::vector<complex128> spf = Add(phi_field, theta_field);

        return std::make_tuple(spf, full_mesh);
    }

    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_structured_fields(const std::vector<double>& phi, const std::vector<double>& theta, const double radius) const
    {
        complex128 propagator = this->compute_propagator(radius);

        auto [S1, S2] = this->compute_s1s2(phi);

        return this->compute_structured_fields(S1, S2, theta, radius);
    }

    double get_g() const {
        double value = 0;

          for(size_t it = 0; it < max_order-1; ++it) {
             double n = (double) it + 1;

              value += ( n * (n + 2.) / (n + 1.) ) * std::real(this->an[it] * std::conj(this->an[it+1]) + this->bn[it] * std::conj(this->bn[it+1]) );
              value += ( (2. * n + 1. ) / ( n * (n + 1.) ) )  * std::real( this->an[it] * std::conj(this->bn[it]) );
          }
          return value * 4. / ( get_Qsca() * pow(size_parameter, 2) );
      }

    double get_Qsca() const {
        // std::cout<<area<< "  " <<this->an[0]<<" dsdas\n";
        double value = 0;

        for(size_t it = 0; it < max_order; ++it){
            double n = (double) it + 1;
            value += (2.* n + 1.) * ( pow( std::abs(this->an[it]), 2) + pow( std::abs(this->bn[it]), 2)  );
        }
        // std::cout<<"dsdas\n";
        return value * 2. / pow( size_parameter, 2.); 
    }

    double get_Qext() const {
        double value = 0;
        for(size_t it = 0; it < max_order; ++it)
        {
            double n = (double) it + 1;
            value += (2.* n + 1.) * std::real( this->an[it] + this->bn[it] );

        }
        return value * 2. / pow( size_parameter, 2.);
    }

    double get_Qback() const {
        complex128 value = 0;

        for(size_t it = 0; it < max_order-1; ++it)
        {
            double n = (double) it + 1;

            value += (2. * n + 1) * pow(-1., n) * ( this->an[it] - this->bn[it] ) ;
        }

        value = pow( std::abs(value), 2. ) / pow( size_parameter, 2. );
        return std::abs(value);
    }

    //--------------------------------------------------------PYTHON-------------------
    std::tuple<py::array_t<complex128>, py::array_t<complex128>>
    get_s1s2_py(const std::vector<double> &phi) const
    {
        auto [S1, S2] = this->compute_s1s2(phi);

        return std::make_tuple(vector_to_numpy(S1), vector_to_numpy(S2));
    }


    std::tuple<py::array_t<complex128>, py::array_t<complex128>, py::array_t<double>, py::array_t<double>>
    get_full_structured_fields_py(size_t &sampling, double& distance) const
    {
        auto [phi_field, theta_field, theta, phi] = this->compute_full_structured_fields(sampling, distance);

        py::array_t<complex128>
            phi_field_py = vector_to_numpy(phi_field, {sampling, sampling}),
            theta_field_py = vector_to_numpy(theta_field, {sampling, sampling});

        py::array_t<double>
            theta_py = vector_to_numpy(theta),
            phi_py = vector_to_numpy(phi);

        phi_field_py = phi_field_py.attr("transpose")();
        theta_field_py = theta_field_py.attr("transpose")();

        return std::make_tuple(
            phi_field_py,
            theta_field_py,
            phi_py,
            theta_py
        );
    }

    std::tuple<py::array_t<complex128>, py::array_t<complex128>>
    get_unstructured_fields_py(const std::vector<double>& phi, const std::vector<double>& theta, const double radius) const
    {
        auto [theta_field, phi_field] = this->compute_unstructured_fields(phi, theta, radius);

        return std::make_tuple(vector_to_numpy(phi_field), vector_to_numpy(theta_field));
    }

    std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const {
        std::vector<complex128>
            S1(phi.size(), 0.0),
            S2(phi.size(), 0.0);

        std::vector<double>
            prefactor = get_prefactor(),
            mu;

        mu.reserve(phi.size());

        for (double phi : phi)
            mu.push_back( cos( phi - PI / 2.0 ) );


        for (unsigned int i = 0; i < phi.size(); i++){
            auto [pin, taun] = VSH::SPHERICAL::MiePiTau(mu[i], max_order);

            for (unsigned int m = 0; m < max_order ; m++){
                S1[i] += prefactor[m] * ( this->an[m] * pin[m] +  this->bn[m] * taun[m] );
                S2[i] += prefactor[m] * ( this->an[m] * taun[m] + this->bn[m] * pin[m]  );
            }
        }

        return std::make_tuple(S1, S2);
    }

};

