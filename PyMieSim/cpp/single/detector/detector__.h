#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <string>

#include <fibonacci/fibonacci.h>
#include <single/source/source.h>
#include <single/optical_interface/optical_interface.h>
#include <single/scatterer/base_scatterer/base_scatterer.h>
#include <mode_field/mode_field.h>
#include <utils/math.h>
#include <utils/constants.h>

using complex128 = std::complex<double>;


class BaseDetector {
public:
    size_t sampling = 0;
    double numerical_aperture = 0.0;
    double cache_numerical_aperture = 0.0;
    double phi_offset = 0.0;
    double gamma_offset = 0.0;
    double polarization_filter = 0.0;
    double medium_refractive_index = 1.0;
    std::shared_ptr<BaseSource> source;
    bool mean_coupling = true;


    std::vector<size_t> indices;
    double max_angle = 0.0;
    double min_angle = 0.0;
    std::vector<complex128> scalar_field;
    FibonacciMesh fibonacci_mesh;
    bool is_coherent = true;

    ModeID mode_id;
    ModeField mode_field;

    // Interface helper + cached per-ray Fresnel amplitude transmissions
    SnellFresnelInterface snell_interface;
    std::vector<double> interface_t_s;
    std::vector<double> interface_t_p;

    BaseDetector() = default;
    virtual ~BaseDetector() = default;

    BaseDetector(
        const size_t _sampling,
        const double _numerical_aperture,
        const double _cache_numerical_aperture,
        const double _phi_offset,
        const double _gamma_offset,
        const double _polarization_filter,
        const double _medium_refractive_index,
        std::shared_ptr<BaseSource> _source,
        const bool _mean_coupling = true)
    :   sampling(_sampling),
        numerical_aperture(_numerical_aperture),
        cache_numerical_aperture(_cache_numerical_aperture),
        phi_offset(_phi_offset),
        gamma_offset(_gamma_offset),
        polarization_filter(_polarization_filter),
        medium_refractive_index(_medium_refractive_index),
        source(std::move(_source)),
        mean_coupling(_mean_coupling)
    {}

    // Media refractive indices for interface calculations
    double scatterer_medium_refractive_index = 1.0; // n_s

    // Detector medium refractive index
    double detector_medium_refractive_index  = 1.0; // n_d

    /**
     * @brief Computes the coupling coefficient for a given scatterer.
     * @param scatterer The scatterer for which the coupling coefficient is computed.
     * @return The coupling coefficient.
     * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
     */
    virtual double get_coupling(const BaseScatterer& scatterer) const = 0;

    /**
     * @brief Computes the structured scalar field for the detector.
     * @param sampling The number of samples to use for the structured scalar field.
     * @return A numpy array containing the structured scalar field.
     * @note This function computes the scalar field based on the mode field and the Fibonacci mesh.
     */
    virtual std::vector<complex128> get_structured_scalarfield(const size_t sampling) const = 0;

    /**
     * @brief Squares each element in the input vector.
     * @tparam T Type of the elements in the vector.
     * @param array The input vector to be squared.
     * @note This function modifies the input vector in place, squaring each element.
     */
    template <typename T> inline void square_array(std::vector<T>& array);

    /**
     * @brief Computes the squared norm of a 1-norm vector.
     * @tparam T Type of the elements in the vector.
     * @param array The input vector.
     * @return The squared norm of the vector.
     * @note This function computes the squared 1-norm, which is the square of the absolute value of the sum of the elements.
     */
    template <typename T> inline double get_norm1_squared(const std::vector<T>& array) const;

    /**
     * @brief Applies the interface transmission coefficients to the theta and phi fields.
     * @param theta_field The field in the theta direction.
     * @param phi_field The field in the phi direction.
     * @note This function modifies the input fields in place, applying the Fresnel transmission coefficients.
     */
    void apply_interface_transmission_to_fields(std::vector<complex128>& theta_field, std::vector<complex128>& phi_field) const;

    /**
     * @brief Computes the squared norm of a 2-norm vector.
     * @tparam T Type of the elements in the vector.
     * @param array The input vector.
     * @return The squared norm of the vector.
     * @note This function computes the squared 2-norm, which is the sum of the squares of the absolute values of the elements.
     */
    template <typename T> inline double get_norm2_squared(const std::vector<T>& array) const;

    /**
     * @brief Computes the Poynting vector value for a given scatterer at a specified distance.
     * @param scatterer The scatterer for which the Poynting vector is computed.
     * @param distance The distance at which to compute the Poynting vector (default is 1).
     * @return The Poynting vector value.
     */
    std::vector<double> get_poynting_field(const BaseScatterer& scatterer, double distance = 1) const;

    /**
     * @brief Computes the energy flow for a given scatterer at a specified distance.
     * @param scatterer The scatterer for which the energy flow is computed.
     * @return The total energy flow value.
     */
    double get_energy_flow(const BaseScatterer& scatterer) const;

private:
    static inline double clamp_m1_p1(double x) {
        return std::max(-1.0, std::min(1.0, x));
    }
};


class Photodiode : public BaseDetector {
public:
    std::string mode_number;
    bool is_coherent = true;
    double rotation = 0.0;

public:
    Photodiode() = default;

    Photodiode(
        const std::string &_mode_number,
        const size_t _sampling,
        const double _numerical_aperture,
        const double _cache_numerical_aperture,
        const double _phi_offset,
        const double _gamma_offset,
        const double _polarization_filter,
        const bool _mean_coupling,
        const double _medium_refractive_index = 1.0)
    :   BaseDetector(
            _sampling,
            _numerical_aperture,
            _cache_numerical_aperture,
            _phi_offset,
            _gamma_offset,
            _polarization_filter,
            _medium_refractive_index,
            nullptr,
            _mean_coupling
        ),
        mode_number(_mode_number)
    {
        this->is_coherent = false;
        this->initialize();
    }

    double get_coupling(const BaseScatterer& scatterer) const override;

    [[nodiscard]] std::vector<complex128> get_structured_scalarfield(const size_t sampling) const override;

private:
    /**
     * @brief Computes the coupling coefficient for a given scatterer, considering coherent or non-coherent modes.
     * @param scatterer The scatterer for which the coupling coefficient is computed.
     * @return The coupling coefficient.
     * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
     */
    double get_coupling_point(const BaseScatterer& scatterer) const;

    /**
     * @brief Computes the coupling coefficient for a given scatterer, considering coherent modes.
     * @param scatterer The scatterer for which the coupling coefficient is computed.
     * @return The coupling coefficient.
     * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
     */
    double get_coupling_mean(const BaseScatterer& scatterer) const;

    /**
     * @brief Initializes the detector with the given medium refractive index.
     * @note This function sets up the mode field and Fibonacci mesh based on the provided parameters.
     */
    void initialize();

    /**
     * @brief Parses the mode number string to extract mode family and numbers.
     * @param mode_number The mode number string in the format "LP01", "HG12", etc.
     * @note This function sets the mode_id based on the parsed values.
     * @throws std::invalid_argument if the mode number string is not in the expected format.
     */
    void parse_mode(const std::string& mode_number);

    /**
     * @brief Computes the projected fields based on the theta and phi fields.
     * @param theta_field The field in the theta direction.
     * @param phi_field The field in the phi direction.
     * @return A tuple containing the horizontal and vertical projections of the fields.
     * @note This function computes the projections based on the spherical coordinates of the Fibonacci mesh.
     */
    std::tuple<std::vector<complex128>, std::vector<complex128>>
    get_projected_farfields(const std::vector<complex128>& theta_field, const std::vector<complex128>& phi_field) const;

    /**
     * @brief Applies a scalar field to the horizontal and vertical projections of the coupling coefficients.
     * @param field0 The horizontal projection field.
     * @param field1 The vertical projection field.
     * @note This function modifies the input fields in place, applying the scalar field to both projections.
     */
    void apply_scalar_field(std::vector<complex128> &field0, std::vector<complex128> &field1) const;

    /**
     * @brief Applies a polarization filter to the coupling coefficients.
     * @tparam T Type of the coupling coefficients.
     * @param coupling_theta The coupling coefficient for the theta direction.
     * @param coupling_phi The coupling coefficient for the phi direction.
     * @param polarization_filter The polarization filter value.
     */
    template <typename T> inline void apply_polarization_filter(T& coupling_theta, T& coupling_phi, double polarization_filter) const;

public:
    void print_properties(int precision) const {
        printf("Photodiode Detector Properties:\n");
        printf("  Mode Number: %s\n", this->mode_number.c_str());
        printf("  Sampling: %zu\n", this->sampling);
        printf("  Numerical Aperture: %.*f\n", precision, this->numerical_aperture);
        printf("  Cache Numerical Aperture: %.*f\n", precision, this->cache_numerical_aperture);
        printf("  Phi Offset (radians): %.*f\n", precision, this->phi_offset);
        printf("  Gamma Offset (radians): %.*f\n", precision, this->gamma_offset);
        printf("  Polarization Filter: %.*f\n", precision, this->polarization_filter);
        printf("  Medium Refractive Index: %.*f\n", precision, this->medium_refractive_index);
        printf("  Scatterer Medium Refractive Index: %.*f\n", precision, this->scatterer_medium_refractive_index);
        printf("  Detector Medium Refractive Index: %.*f\n", precision, this->detector_medium_refractive_index);
        printf("  Is Coherent: %s\n", this->is_coherent ? "True" : "False");
        printf("  Mean Coupling: %s\n", this->mean_coupling ? "True" : "False");
    }

};



class CoherentMode : public BaseDetector {
public:
    std::string mode_number;
    double rotation = 0.0;
    bool is_coherent = true;

public:
    CoherentMode() = default;

    CoherentMode(
        const std::string &_mode_number,
        const size_t _sampling,
        const double _numerical_aperture,
        const double _cache_numerical_aperture,
        const double _phi_offset,
        const double _gamma_offset,
        const double _polarization_filter,
        const double _rotation,
        const bool _mean_coupling,
        const double _medium_refractive_index = 1.0)
    :   BaseDetector(
            _sampling,
            _numerical_aperture,
            _cache_numerical_aperture,
            _phi_offset,
            _gamma_offset,
            _polarization_filter,
            _medium_refractive_index,
            nullptr,
            _mean_coupling
        ),
        mode_number(_mode_number),
        rotation(_rotation)
    {
        this->is_coherent = true;
        this->initialize();

    }

    double get_coupling(const BaseScatterer& scatterer) const override;

    [[nodiscard]] std::vector<complex128> get_structured_scalarfield(const size_t sampling) const override;

private:

    /**
     * @brief Computes the coupling coefficient for a given scatterer, considering coherent modes.
     * @param scatterer The scatterer for which the coupling coefficient is computed.
     * @return The coupling coefficient.
     * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
     */
    double get_coupling_point(const BaseScatterer& scatterer) const;

    /**
     * @brief Computes the coupling coefficient for a given scatterer, considering coherent modes.
     * @param scatterer The scatterer for which the coupling coefficient is computed.
     * @return The coupling coefficient.
     * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
     */
    double get_coupling_mean(const BaseScatterer& scatterer) const;

    /**
     * @brief Initializes the detector with the given medium refractive index.
     * @note This function sets up the mode field and Fibonacci mesh based on the provided parameters.
     */
    void initialize();

    /**
     * @brief Parses the mode number string to extract mode family and numbers.
     * @param mode_number The mode number string in the format "LP01", "HG12", etc.
     * @note This function sets the mode_id based on the parsed values.
     * @throws std::invalid_argument if the mode number string is not in the expected format.
     */
    void parse_mode(const std::string& mode_number);

    /**
     * @brief Computes the projected fields based on the theta and phi fields.
     * @param theta_field The field in the theta direction.
     * @param phi_field The field in the phi direction.
     * @return A tuple containing the horizontal and vertical projections of the fields.
     * @note This function computes the projections based on the spherical coordinates of the Fibonacci mesh.
     */
    std::tuple<std::vector<complex128>, std::vector<complex128>>
    get_projected_farfields(const std::vector<complex128>& theta_field, const std::vector<complex128>& phi_field) const;

    /**
     * @brief Applies a scalar field to the horizontal and vertical projections of the coupling coefficients.
     * @param field0 The horizontal projection field.
     * @param field1 The vertical projection field.
     * @note This function modifies the input fields in place, applying the scalar field to both projections.
     */
    void apply_scalar_field(std::vector<complex128> &field0, std::vector<complex128> &field1) const;

    /**
     * @brief Applies a polarization filter to the coupling coefficients.
     * @tparam T Type of the coupling coefficients.
     * @param coupling_theta The coupling coefficient for the theta direction.
     * @param coupling_phi The coupling coefficient for the phi direction.
     * @param polarization_filter The polarization filter value.
     */
    template <typename T> inline void apply_polarization_filter(T& coupling_theta, T& coupling_phi, double polarization_filter) const;

public:
    void print_properties(int precision) const {
        printf("Photodiode Detector Properties:\n");
        printf("  Mode Number: %s\n", this->mode_number.c_str());
        printf("  Sampling: %zu\n", this->sampling);
        printf("  Numerical Aperture: %.*f\n", precision, this->numerical_aperture);
        printf("  Cache Numerical Aperture: %.*f\n", precision, this->cache_numerical_aperture);
        printf("  Phi Offset (radians): %.*f\n", precision, this->phi_offset);
        printf("  Gamma Offset (radians): %.*f\n", precision, this->gamma_offset);
        printf("  Polarization Filter: %.*f\n", precision, this->polarization_filter);
        printf("  Medium Refractive Index: %.*f\n", precision, this->medium_refractive_index);
        printf("  Scatterer Medium Refractive Index: %.*f\n", precision, this->scatterer_medium_refractive_index);
        printf("  Detector Medium Refractive Index: %.*f\n", precision, this->detector_medium_refractive_index);
        printf("  Is Coherent: %s\n", this->is_coherent ? "True" : "False");
        printf("  Mean Coupling: %s\n", this->mean_coupling ? "True" : "False");
    }

};

class IntegratingSphere : public BaseDetector {
public:
    /**
     * @brief Construct an integrating sphere detector that collects over 4π steradians.
     *
     * This detector samples directions over the full sphere (theta in [0, pi]).
     * It is intended as a power collector: it integrates the scattered field intensity
     * over all directions, optionally applying a transmission model at a planar interface
     * between scatterer medium and detector medium (SnellFresnelInterface).
     *
     * Notes:
     * - numerical_aperture and cache_numerical_aperture are ignored for this detector.
     * - phi_offset and gamma_offset can still be used to rotate the sampling pattern,
     *   but they do not restrict the acceptance.
     * - If you want an ideal integrating sphere with no interface physics, disable
     *   Fresnel transmission in SnellFresnelInterface (or set n_s = n_d and disable).
     */
    IntegratingSphere() = default;

    IntegratingSphere(
        const size_t _sampling,
        const double _polarization_filter)
    :   BaseDetector(
            _sampling,
            0.0,  /* numerical_aperture */
            0.0,  /* cache_numerical_aperture */
            0.0,  /* phi_offset */
            0.0,  /* gamma_offset */
            _polarization_filter,
            1.0,
            nullptr,
            false
        )
    {
        this->is_coherent = false;
        this->initialize();
    }

    /**
     * @brief Compute the collected power from the scatterer over 4π.
     *
     * For incoherent detection, the coupling is proportional to the integral
     * of |E_theta|^2 + |E_phi|^2 over the full solid angle.
     *
     * If Fresnel transmission is enabled in snell_interface, a scalar amplitude
     * transmission is applied to the fields before intensity computation.
     */
    double get_coupling(const BaseScatterer& scatterer) const override;

    /**
     * @brief Not used for integrating sphere.
     *
     * This detector does not have a mode field. The function returns an empty vector.
     */
    [[nodiscard]] std::vector<complex128> get_structured_scalarfield(const size_t sampling) const override;

    void print_properties(int precision) const {
        printf("Integrating Sphere Detector Properties:\n");
        printf("  Sampling: %zu\n", this->sampling);
        printf("  Polarization Filter: %.*f\n", precision, this->polarization_filter);
        printf("  Medium Refractive Index: %.*f\n", precision, this->medium_refractive_index);
        printf("  Scatterer Medium Refractive Index: %.*f\n", precision, this->scatterer_medium_refractive_index);
        printf("  Detector Medium Refractive Index: %.*f\n", precision, this->detector_medium_refractive_index);
        printf("  Is Coherent: %s\n", this->is_coherent ? "True" : "False");
        printf("  Mean Coupling: %s\n", this->mean_coupling ? "True" : "False");
    }

private:
    /**
     * @brief Initialize 4π sampling mesh and cache Fresnel transmission coefficients.
     *
     * Sets min_angle = 0 and max_angle = pi to cover 4π.
     */
    void initialize();

};
