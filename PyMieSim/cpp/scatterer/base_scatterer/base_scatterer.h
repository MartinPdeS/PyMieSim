#pragma once

#include <vector>
#include "source/source.h"
#include "utils/numpy_interface.h"
#include "fibonacci/fibonacci.h"
#include "full_mesh/full_mesh.h"

typedef std::complex<double> complex128;


class BaseScatterer {
public:
    size_t max_order;
    BaseSource source;
    double size_parameter;
    double size_parameter_squared;
    double area;
    double medium_refractive_index;
    std::vector<size_t> indices;
    std::vector<complex128> an, bn, cn, dn;  // For spherical scatterers

    BaseScatterer() = default;

    // VIRTUAL INTERFACE ----------------------------------------------
    virtual ~BaseScatterer() = default;

    virtual std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const = 0;
    virtual double get_Qsca() const {throw std::logic_error{"Function not implementend!"};};
    virtual double get_Qext() const {throw std::logic_error{"Function not implementend!"};};
    virtual double get_Qback() const {throw std::logic_error{"Function not implementend!"};};
    virtual double get_g() const {throw std::logic_error{"Function not implementend!"};};

    virtual void compute_size_parameter() = 0;
    virtual void compute_area() = 0;


    // GETTERS --------------------------------------------------------
    std::vector<complex128> get_an() const { return this->an; }
    std::vector<complex128> get_bn() const { return this->bn; }
    std::vector<complex128> get_cn() const { return this->cn; }
    std::vector<complex128> get_dn() const { return this->dn; }



    complex128 get_a() const { return this->an.back(); }
    complex128 get_b() const { return this->bn.back(); }
    complex128 get_c() const { return this->cn.back(); }
    complex128 get_d() const { return this->dn.back(); }

    double get_a1_abs() const { return abs(this->an[0]); }
    double get_b1_abs() const { return abs(this->bn[0]); }
    double get_c1_abs() const { return abs(this->cn[0]); }
    double get_d1_abs() const { return abs(this->dn[0]); }

    double get_a2_abs() const { return abs(this->an[1]); }
    double get_b2_abs() const { return abs(this->bn[1]); }
    double get_c2_abs() const { return abs(this->cn[1]); }
    double get_d2_abs() const { return abs(this->dn[1]); }

    double get_a3_abs() const { return abs(this->an[2]); }
    double get_b3_abs() const { return abs(this->bn[2]); }
    double get_c3_abs() const { return abs(this->cn[2]); }
    double get_d3_abs() const { return abs(this->dn[2]); }


    double get_an_abs() const { return abs(this->an.back()); }
    double get_bn_abs() const { return abs(this->bn.back()); }
    double get_cn_abs() const { return abs(this->cn.back()); }
    double get_dn_abs() const { return abs(this->dn.back()); }

    pybind11::array_t<complex128> get_an_py(size_t _max_order) { _max_order = (_max_order == 0 ? this->max_order : _max_order); return _vector_to_numpy(this->an, {_max_order}); }
    pybind11::array_t<complex128> get_bn_py(size_t _max_order) { _max_order = (_max_order == 0 ? this->max_order : _max_order); return _vector_to_numpy(this->bn, {_max_order}); }
    pybind11::array_t<complex128> get_cn_py(size_t _max_order) { _max_order = (_max_order == 0 ? this->max_order : _max_order); return _vector_to_numpy(this->cn, {_max_order}); }
    pybind11::array_t<complex128> get_dn_py(size_t _max_order) { _max_order = (_max_order == 0 ? this->max_order : _max_order); return _vector_to_numpy(this->dn, {_max_order}); }



    // COEFFICIENT ----------------------------------------------------
    double get_g_with_fields(size_t sampling) const;
    double get_Qabs() const {return get_Qext() - get_Qsca();};
    double get_Qpr() const {return get_Qext() - get_g() * get_Qsca();};
    double get_Qforward() const {return get_Qsca() - get_Qback();};
    double get_Qratio() const {return get_Qback() / get_Qsca();};
    double get_Cback() const {return get_Qback() * area;};
    double get_Cforward() const {return get_Qforward() * area;};
    double get_Cratio() const {return get_Qratio() * area;};

    double get_Csca() const {return get_Qsca() * area;};
    double get_Cext() const {return get_Qext() * area;};
    double get_Cabs() const {return get_Qabs() * area;};
    double get_Cpr() const {return get_Qpr() * area;};


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