#pragma once

#include <vector>
#include "source/source.h"
#include "utils/numpy_interface.h"
#include "fibonacci/fibonacci.h"
#include "full_mesh/full_mesh.h"

typedef std::complex<double> complex128;
namespace py = pybind11;



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
    virtual double get_Qforward() const {throw std::logic_error{"Function not implementend!"};};
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

    /**
     * @brief Computes the asymmetry factor g using the fields at a given sampling.
     * @param sampling The number of sampling points to use for the computation.
     * @return The asymmetry factor g.
     * @note This method computes the asymmetry factor g based on the fields at a specified sampling.
     * It uses the scattering amplitudes S1 and S2 to calculate the asymmetry factor.
     */
    double get_g_with_fields(size_t sampling) const;


    /**
     * @brief Computes the scattering efficiency Qsca.
     * @return The scattering efficiency Qsca.
     * @note This method computes the scattering efficiency based on the coefficients an and bn.
     */
    double get_Qabs() const {return get_Qext() - get_Qsca();};

    /**
     * @brief Computes the scattering efficiency Qpr.
     * @return The scattering efficiency Qpr.
     * @note This method computes the scattering efficiency based on the coefficients an and bn.
     */
    double get_Qpr() const {return get_Qext() - get_g() * get_Qsca();};

    // /**
    //  * @brief Computes the scattering efficiency Qforward.
    //  * @return The scattering efficiency Qforward.
    //  * @note This method computes the scattering efficiency based on the coefficients an and bn.
    //  */
    // double get_Qforward() const {return get_Qsca() - get_Qback();};

    /**
     * @brief Computes the backscattering efficiency Qback.
     * @return The backscattering efficiency Qback.
     * @note This method computes the backscattering efficiency based on the coefficients an and bn.
     */
    double get_Qratio() const {return get_Qback() / get_Qsca();};

    /**
     * @brief Computes the backscattering coefficient Cback.
     * @return The backscattering coefficient Cback.
     * @note This method computes the backscattering coefficient based on the Qback and the cross section.
     */
    double get_Cback() const {return get_Qback() * this->cross_section;};

    /**
     * @brief Computes the forward scattering coefficient Cforward.
     * @return The forward scattering coefficient Cforward.
     * @note This method computes the forward scattering coefficient based on the Qforward and the cross section.
     */
    double get_Cforward() const {return get_Qforward() * this->cross_section;};

    /**
     * @brief Computes the scattering coefficient Cratio.
     * @return The scattering coefficient Cratio.
     * @note This method computes the scattering coefficient based on the Qratio and the cross section.
     */
    double get_Cratio() const {return get_Qratio() * this->cross_section;};

    /**
     * @brief Computes the scattering coefficient Csca.
     * @return The scattering coefficient Csca.
     * @note This method computes the scattering coefficient based on the Qsca and the cross section.
     */
    double get_Csca() const {return get_Qsca() * this->cross_section;};

    /**
     * @brief Computes the extinction coefficient Cext.
     * @return The extinction coefficient Cext.
     * @note This method computes the extinction coefficient based on the Qext and the cross section.
     */
    double get_Cext() const {return get_Qext() * this->cross_section;};

    /**
     * @brief Computes the absorption coefficient Cabs.
     * @return The absorption coefficient Cabs.
     * @note This method computes the absorption coefficient based on the Qabs and the cross section.
     */
    double get_Cabs() const {return get_Qabs() * this->cross_section;};

    /**
     * @brief Computes the forward scattering coefficient Cpr.
     * @return The forward scattering coefficient Cpr.
     * @note This method computes the forward scattering coefficient based on the Qpr and the cross section.
     */
    double get_Cpr() const {return get_Qpr() * this->cross_section;};

    // GENERAL METHODS ----------------------------------------------------
    BaseScatterer(const size_t max_order, const BaseSource &source, const double medium_refractive_index)
    : max_order(max_order), source(source), medium_refractive_index(medium_refractive_index){}

    /**
     * @brief Computes the Wiscombe criterion.
     * @param size_parameter The size parameter.
     * @return The Wiscombe criterion.
     */
    size_t get_wiscombe_criterion(const double size_parameter) const;

    /**
     * @brief Computes the coefficients cn and dn for a scatterer.
     * @param max_order The maximum order of the coefficients to compute.
     */
    std::vector<double> get_prefactor() const;

    /**
     * @brief Computes the structured fields.
     * @param S1 The first set of scattering coefficients.
     * @param S2 The second set of scattering coefficients.
     * @param theta The scattering angles.
     * @param radius The radius of the scatterer.
     * @return A tuple containing the phi and theta fields.
     */
    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_structured_fields(const std::vector<complex128>& S1, const std::vector<complex128>& S2, const std::vector<double>& theta, const double& radius = 1.) const;

    /**
     * @brief Computes the structured fields based on phi and theta angles.
     * @param phi The azimuthal angles in radians.
     * @param theta The polar angles in radians.
     * @param radius The radius of the scatterer.
     * @return A tuple containing the phi and theta fields.
     */
    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_structured_fields(const std::vector<double>& phi, const std::vector<double>& theta, const double radius) const;

    /**
     * @brief Computes the unstructured fields based on phi and theta angles.
     * @param phi The azimuthal angles in radians.
     * @param theta The polar angles in radians.
     * @param radius The radius of the scatterer.
     * @return A tuple containing the phi and theta fields.
     */
    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_unstructured_fields(const std::vector<double>& phi, const std::vector<double>& theta, const double radius) const;

    /**
     * @brief Computes the unstructured fields for a Fibonacci mesh.
     * @param fibonacci_mesh The Fibonacci mesh object.
     * @param radius The radius of the scatterer.
     * @return A tuple containing the phi and theta fields.
     */
    std::tuple<std::vector<complex128>, std::vector<complex128>>
    compute_unstructured_fields(const FibonacciMesh& fibonacci_mesh, const double radius) const;

    /**
     * @brief Computes the full structured fields for a given sampling and radius.
     * @param sampling The number of sampling points.
     * @param radius The radius of the scatterer.
     * @return A tuple containing the phi field, theta field, phi angles, and theta angles.
     */
    std::tuple<std::vector<complex128>, std::vector<complex128>, std::vector<double>, std::vector<double>>
    compute_full_structured_fields(const size_t& sampling, const double& radius) const;

    /**
     * @brief Computes the pi and tau coefficients for a given mu and maximum order.
     * @param mu The cosine of the scattering angle.
     * @param max_order The maximum order of the coefficients to compute.
     * @return A tuple containing the pi and tau coefficients as vectors.
     */
    std::tuple<std::vector<complex128>, std::vector<complex128>>
    get_pi_tau(const double& mu, const size_t& max_order) const;

    /**
     * @brief Computes the pi and tau coefficients for a given mu and maximum order.
     * @param mu The cosine of the scattering angle.
     * @param max_order The maximum order of the coefficients to compute.
     * @param pin Output vector for pi coefficients.
     * @param taun Output vector for tau coefficients.
     */
    void get_pi_tau(double mu, size_t max_order, complex128 *pin, complex128 *taun) const;

    /**
     * @brief Computes the propagator for a given radius.
     * @param radius The radius of the scatterer.
     * @return The propagator as a complex128 value.
     */
    complex128 get_propagator(const double &radius) const;

    /**
     * @brief Computes the dn coefficients for a given nmx and z.
     * @param nmx The maximum size parameter.
     * @param z The complex refractive index.
     * @return A vector of dn coefficients.
     */
    std::vector<complex128> compute_dn(double nmx, complex128 z)  const;

    //- PYTHON INTERFACE --------------------------------------------------------------------------
    std::tuple<py::array_t<complex128>, py::array_t<complex128>>
    get_unstructured_fields_py(const std::vector<double>& phi, const std::vector<double>& theta, const double radius) const;

    std::tuple<py::array_t<complex128>, py::array_t<complex128>>
    get_s1s2_py(const std::vector<double> &phi) const;

    std::tuple<py::array_t<complex128>, py::array_t<complex128>, py::array_t<double>, py::array_t<double>>
    get_full_structured_fields_py(size_t &sampling, double& distance) const;

    std::tuple<std::vector<double>, FullSteradian>
    compute_full_structured_spf(const size_t sampling, const double radius = 1.0) const;

    complex128 get_coefficient_py(const std::string &type, const size_t order);
};