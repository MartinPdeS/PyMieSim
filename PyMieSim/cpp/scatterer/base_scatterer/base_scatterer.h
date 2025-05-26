#pragma once

#include <vector>
#include "source/source.h"
#include "utils/numpy_interface.h"
#include "fibonacci/fibonacci.h"
#include "full_mesh/full_mesh.h"

typedef std::complex<double> complex128;


#define DEFINE_COEFFICIENTS_GETTERS(name) \
    std::vector<complex128> name##n; \
    pybind11::array_t<complex128> get_##name##n_list_py() {return _vector_to_numpy(this->name##n); } \
    double get_##name##1() const { return abs(this->name##n[0]); } \
    double get_##name##2() const { return abs(this->name##n[1]); } \
    double get_##name##3() const { return abs(this->name##n[2]); } \
    complex128 get_##name##1_complex128() const { return this->name##n[0]; } \
    complex128 get_##name##2_complex128() const { return this->name##n[1]; } \
    complex128 get_##name##3_complex128() const { return this->name##n[2]; }


class BaseScatterer {
public:
    size_t max_order;
    BaseSource source;
    double size_parameter;
    double size_parameter_squared;
    double cross_section;
    double medium_refractive_index;
    std::vector<size_t> indices;

    BaseScatterer() = default;

    // VIRTUAL INTERFACE ----------------------------------------------
    virtual ~BaseScatterer() = default;

    virtual std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const = 0;
    virtual double get_Qsca() const {throw std::logic_error{"Function not implementend!"};};
    virtual double get_Qext() const {throw std::logic_error{"Function not implementend!"};};
    virtual double get_Qback() const {throw std::logic_error{"Function not implementend!"};};
    virtual double get_g() const {throw std::logic_error{"Function not implementend!"};};

    virtual void compute_size_parameter() = 0;
    virtual void compute_cross_section() = 0;


    // SPHERICAL SCATTERER GETTERS --------------------------------------------------------
    DEFINE_COEFFICIENTS_GETTERS(a)
    DEFINE_COEFFICIENTS_GETTERS(b)
    DEFINE_COEFFICIENTS_GETTERS(c)
    DEFINE_COEFFICIENTS_GETTERS(d)

    // CYLINDRICAL SCATTERER GETTERS --------------------------------------------------------
    DEFINE_COEFFICIENTS_GETTERS(a1)
    DEFINE_COEFFICIENTS_GETTERS(b1)
    DEFINE_COEFFICIENTS_GETTERS(a2)
    DEFINE_COEFFICIENTS_GETTERS(b2)


    // COEFFICIENT ----------------------------------------------------
    double get_g_with_fields(size_t sampling) const;
    double get_Qabs() const {return get_Qext() - get_Qsca();};
    double get_Qpr() const {return get_Qext() - get_g() * get_Qsca();};
    double get_Qforward() const {return get_Qsca() - get_Qback();};
    double get_Qratio() const {return get_Qback() / get_Qsca();};
    double get_Cback() const {return get_Qback() * this->cross_section;};
    double get_Cforward() const {return get_Qforward() * this->cross_section;};
    double get_Cratio() const {return get_Qratio() * this->cross_section;};

    double get_Csca() const {return get_Qsca() * this->cross_section;};
    double get_Cext() const {return get_Qext() * this->cross_section;};
    double get_Cabs() const {return get_Qabs() * this->cross_section;};
    double get_Cpr() const {return get_Qpr() * this->cross_section;};


    // GENERAL METHODS ----------------------------------------------------
    BaseScatterer(const size_t max_order, const BaseSource &source, const double medium_refractive_index)
    : max_order(max_order), source(source), medium_refractive_index(medium_refractive_index){}

    size_t get_wiscombe_criterion(const double size_parameter) const;

    std::vector<double> get_prefactor() const;


    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_structured_fields(const std::vector<complex128>& S1, const std::vector<complex128>& S2, const std::vector<double>& theta, const double& radius = 1.) const;

    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_structured_fields(const std::vector<double>& phi, const std::vector<double>& theta, const double radius) const;

    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_unstructured_fields(const std::vector<double>& phi, const std::vector<double>& theta, const double radius) const;

    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_unstructured_fields(const FibonacciMesh& fibonacci_mesh, const double radius) const;

    std::tuple<std::vector<complex128>, std::vector<complex128>, std::vector<double>, std::vector<double>>
    compute_full_structured_fields(const size_t& sampling, const double& radius) const;

    std::tuple<std::vector<complex128>, std::vector<complex128>>
    get_pi_tau(const double& mu, const size_t& max_order) const;

    void
    get_pi_tau(double mu, size_t max_order, complex128 *pin, complex128 *taun) const;

    complex128 get_propagator(const double &radius) const;

    std::vector<complex128>
    compute_dn(double nmx, complex128 z)  const;

    //- PYTHON INTERFACE --------------------------------------------------------------------------
    std::tuple<py::array_t<complex128>, py::array_t<complex128>>
    get_unstructured_fields_py(const std::vector<double>& phi, const std::vector<double>& theta, const double radius) const;

    std::tuple<py::array_t<complex128>, py::array_t<complex128>>
    get_s1s2_py(const std::vector<double> &phi) const;

    std::tuple<py::array_t<complex128>, py::array_t<complex128>, py::array_t<double>, py::array_t<double>>
    get_full_structured_fields_py(size_t &sampling, double& distance) const;

    std::tuple<std::vector<double>, FullSteradian>
    compute_full_structured_spf(const size_t sampling, const double radius = 1.0) const;

    complex128
    get_coefficient_py(const std::string &type, const size_t order);
};