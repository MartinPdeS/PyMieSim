#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <string>

#include <utils/math.h>
#include <utils/constants.h>
#include <material/material.h>
#include <mesh/fibonacci.h>

#include <single/source/source.h>
#include <single/optical_interface/optical_interface.h>
#include <single/scatterer/base_scatterer/base_scatterer.h>
#include <single/mode_field/mode_field.h>


using complex128 = std::complex<double>;


class BaseDetector {
public:
    size_t sampling = 0;
    double numerical_aperture = 0.0;
    double cache_numerical_aperture = 0.0;
    double phi_offset = 0.0;
    double gamma_offset = 0.0;
    double polarization_filter = 0.0;
    std::shared_ptr<BaseMedium> medium;
    bool mean_coupling = true;


    std::vector<size_t> indices;
    double max_angle = 0.0;
    double min_angle = 0.0;
    std::vector<complex128> scalar_field;
    FibonacciMesh fibonacci_mesh;

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
        std::shared_ptr<BaseMedium> _medium,
        const bool _mean_coupling)
    :   sampling(_sampling),
        numerical_aperture(_numerical_aperture),
        cache_numerical_aperture(_cache_numerical_aperture),
        phi_offset(_phi_offset),
        gamma_offset(_gamma_offset),
        polarization_filter(_polarization_filter),
        medium(std::move(_medium)),
        mean_coupling(_mean_coupling)
    {}

    /**
     * @brief Initializes the detector's mesh based on the scatterer's properties and the interface conditions.
     * @param scatterer The scatterer for which the mesh is initialized.
     * @param rotation The rotation angle to apply to the mesh.
     * @note This function sets up the Fibonacci mesh for the detector, computes the scalar field based on the mode field, and applies the Fresnel transmission coefficients to the fields.
     */
    virtual void initialize_mesh(const std::shared_ptr<BaseScatterer> scatterer) = 0;


    /**
     * @brief Computes the coupling coefficient for a given scatterer.
     * @param scatterer The scatterer for which the coupling coefficient is computed.
     * @return The coupling coefficient.
     * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
     */
    virtual double get_coupling(std::shared_ptr<BaseScatterer> scatterer, std::shared_ptr<BaseSource> source) = 0;

    /**
     * @brief Parses the mode number string to extract mode family and numbers.
     * @param mode_number The mode number string in the format "LP01", "HG12", etc.
     * @note This function sets the mode_id based on the parsed values.
     * @throws std::invalid_argument if the mode number string is not in the expected format.
     */
    virtual void parse_mode(const std::string& mode_number);

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
    std::vector<double> get_poynting_field(std::shared_ptr<BaseScatterer> scatterer, std::shared_ptr<BaseSource> source, double distance = 1) const;

    /**
     * @brief Computes the energy flow for a given scatterer at a specified distance.
     * @param scatterer The scatterer for which the energy flow is computed.
     * @return The total energy flow value.
     */
    double get_energy_flow(std::shared_ptr<BaseScatterer> scatterer, std::shared_ptr<BaseSource> source) const;

    /**
     * @brief Prints the properties of the detector with a specified precision.
     * @param precision The number of decimal places to display for floating-point values.
     */
    void print_properties(int precision = 4) const {
        printf("Detector Properties:\n");
        printf("  Sampling: %zu\n", this->sampling);
        printf("  Numerical Aperture: %.*f\n", precision, this->numerical_aperture);
        printf("  Cache Numerical Aperture: %.*f\n", precision, this->cache_numerical_aperture);
        printf("  Phi Offset (radians): %.*f\n", precision, this->phi_offset);
        printf("  Gamma Offset (radians): %.*f\n", precision, this->gamma_offset);
        printf("  Polarization Filter: %.*f\n", precision, this->polarization_filter);
    }

private:
    static inline double clamp_m1_p1(double x) {
        return std::max(-1.0, std::min(1.0, x));
    }
};
