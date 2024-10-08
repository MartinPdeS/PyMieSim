#pragma once

#include "utils/utils.cpp"
#include "single/includes/sources.cpp"
#include "single/includes/fibonacci_mesh.cpp"
#include "utils/VSH.cpp"
#include <cmath>
#include <vector>
#include <complex>
#include <single/includes/base_class.cpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

typedef std::complex<double> complex128;

class BaseSphericalScatterer : public BaseScatterer
{
public:
    std::vector<complex128> an;
    std::vector<complex128> bn;
    std::vector<complex128> cn;
    std::vector<complex128> dn;

    BaseSphericalScatterer() = default;
    virtual ~BaseSphericalScatterer() = default;

    BaseSphericalScatterer(const SOURCE::BaseSource &source, const size_t max_order, const double medium_index) :
    BaseScatterer(max_order, source, medium_index){}

    std::vector<complex128> get_an() const { return an; };
    std::vector<complex128> get_bn() const { return bn; };
    std::vector<complex128> get_cn() const { return cn; };
    std::vector<complex128> get_dn() const { return dn; };

    pybind11::array_t<complex128> get_an_py(size_t _max_order) { _max_order == 0 ? _max_order = this->max_order : _max_order = max_order; return vector_to_numpy(an, {_max_order}); }
    pybind11::array_t<complex128> get_bn_py(size_t _max_order) { _max_order == 0 ? _max_order = this->max_order : _max_order = max_order; return vector_to_numpy(bn, {_max_order}); }
    pybind11::array_t<complex128> get_cn_py(size_t _max_order) { _max_order == 0 ? _max_order = this->max_order : _max_order = max_order; return vector_to_numpy(cn, {_max_order}); }
    pybind11::array_t<complex128> get_dn_py(size_t _max_order) { _max_order == 0 ? _max_order = this->max_order : _max_order = max_order; return vector_to_numpy(dn, {_max_order}); }

    double get_Qforward() const {return get_Qsca() - get_Qback();};
    double get_Qratio() const {return get_Qback() / get_Qsca();};
    double get_Cback() const {return get_Qback() * area;};
    double get_Cforward() const {return get_Qforward() * area;};
    double get_Cratio() const {return get_Qratio() * area;};


    complex128 get_a1() const { return this->an[0]; }
    complex128 get_a2() const { return this->an[1]; }
    complex128 get_a3() const { return this->an[2]; }
    complex128 get_b1() const { return this->bn[0]; }
    complex128 get_b2() const { return this->bn[1]; }
    complex128 get_b3() const { return this->bn[2]; }

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
        double value = 0;

        for(size_t it = 0; it < max_order; ++it){
            double n = (double) it + 1;
            value += (2.* n + 1.) * ( pow( std::abs(this->an[it]), 2) + pow( std::abs(this->bn[it]), 2)  );
        }
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

        value = pow( std::abs(value), 2. ) / pow(size_parameter, 2.);
        return std::abs(value);
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

        return std::make_tuple(std::move(S1), std::move(S2));
    }

};

