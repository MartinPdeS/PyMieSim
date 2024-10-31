#pragma once


#define DEFINE_COMPLEX_VECTOR(name) \
    std::vector<complex128> name##n; \
    std::vector<complex128> get_##name##n() const { return name##n; }

#define DEFINE_GETTER(name, index) \
    complex128 get_##name##index() const { return this->name##n[index - 1]; }

#define DEFINE_GETTER_ABS(name, index) \
    double get_##name##index##_abs() const { return abs(this->name##n[index - 1]); }

#define DEFINE_COEFFICIENTS(name) \
    DEFINE_COMPLEX_VECTOR(name) \
    DEFINE_GETTER(name, 1) \
    DEFINE_GETTER(name, 2) \
    DEFINE_GETTER(name, 3) \
    DEFINE_GETTER(name, 4) \
    DEFINE_GETTER_ABS(name, 1) \
    DEFINE_GETTER_ABS(name, 2) \
    DEFINE_GETTER_ABS(name, 3) \
    DEFINE_GETTER_ABS(name, 4) \
    pybind11::array_t<complex128> get_##name##n_py(size_t _max_order) { _max_order = (_max_order == 0 ? this->max_order : _max_order); return _vector_to_numpy(name##n, {_max_order}); }

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
    DEFINE_COEFFICIENTS(a)
    DEFINE_COEFFICIENTS(b)
    DEFINE_COEFFICIENTS(c)
    DEFINE_COEFFICIENTS(d)

    BaseSphericalScatterer() = default;
    virtual ~BaseSphericalScatterer() = default;

    BaseSphericalScatterer(const SOURCE::BaseSource &source, const size_t max_order, const double medium_index) :
    BaseScatterer(max_order, source, medium_index){}

    double get_Qforward() const {return get_Qsca() - get_Qback();};
    double get_Qratio() const {return get_Qback() / get_Qsca();};
    double get_Cback() const {return get_Qback() * area;};
    double get_Cforward() const {return get_Qforward() * area;};
    double get_Cratio() const {return get_Qratio() * area;};

    double get_g() const {
        double value = 0;

          for(size_t it = 0; it < max_order-1; ++it) {
             double n = (double) it + 1;

              value += ( n * (n + 2.) / (n + 1.) ) * std::real(this->an[it] * std::conj(this->an[it+1]) + this->bn[it] * std::conj(this->bn[it+1]) );
              value += ( (2. * n + 1. ) / ( n * (n + 1.) ) )  * std::real( this->an[it] * std::conj(this->bn[it]) );
          }
          return value * 4. / ( get_Qsca() * size_parameter_squared );
      }

    double get_Qsca() const {
        double value = 0;

        for(size_t it = 0; it < max_order; ++it){
            double n = (double) it + 1;
            value += (2. * n + 1.) * ( pow( std::abs(this->an[it]), 2) + pow( std::abs(this->bn[it]), 2)  );
        }
        return value * 2. / size_parameter_squared;
    }

    double get_Qext() const {
        double value = 0;
        for(size_t it = 0; it < max_order; ++it)
        {
            double n = (double) it + 1;
            value += (2.* n + 1.) * std::real( this->an[it] + this->bn[it] );

        }
        return value * 2. / size_parameter_squared;
    }

    double get_Qback() const {
        complex128 value = 0;

        for(size_t it = 0; it < max_order-1; ++it)
        {
            double n = (double) it + 1;

            value += (2. * n + 1) * pow(-1., n) * ( this->an[it] - this->bn[it] ) ;
        }

        value = pow( std::abs(value), 2. ) / size_parameter_squared;
        return std::abs(value);
    }

    std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const {
        std::vector<complex128> S1, S2;

        S1.reserve(phi.size());
        S2.reserve(phi.size());

        std::vector<double> mu, prefactor = get_prefactor();

        mu.reserve(phi.size());

        for (const double phi : phi)
            mu.push_back( cos( phi - PI / 2.0 ) );

        for (size_t i = 0; i < phi.size(); i++){
            auto [pin, taun] = VSH::SPHERICAL::MiePiTau(mu[i], max_order);
            complex128 S1_temp = 0., S2_temp = 0.;

            for (size_t m = 0; m < max_order ; m++){
                S1_temp += prefactor[m] * ( this->an[m] * pin[m] +  this->bn[m] * taun[m] );
                S2_temp += prefactor[m] * ( this->an[m] * taun[m] + this->bn[m] * pin[m]  );
            }
            S1.push_back(S1_temp);
            S2.push_back(S2_temp);
        }

        return std::make_tuple(std::move(S1), std::move(S2));
    }

};

