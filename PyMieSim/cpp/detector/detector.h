#pragma once

#include <vector>
#include <complex>
#include <cmath> // For std::isnan and std::pow
#include <stdexcept>
#include "utils/special_function.cpp"
#include "fibonacci/fibonacci.h"
#include "utils/numpy_interface.h"
#include <scatterer/base_scatterer/base_scatterer.h>
#include "../mode_field/mode_field.h"

using complex128 = std::complex<double>;
#define EPSILON0 (double)8.854187817620389e-12
#define C (double)299792458.0


class Detector {
    public:
        std::string mode_number;
        size_t sampling = 0;
        double numerical_aperture = 0.0;
        double cache_numerical_aperture = 0.0;
        double phi_offset = 0.0;
        double gamma_offset = 0.0;
        double polarization_filter = 0.0;
        double rotation = 0.0;
        bool coherent = true;
        bool mean_coupling = true;
        double medium_refractive_index;
        double max_angle = 0;
        double min_angle = 0;
        std::vector<complex128> scalar_field;
        FibonacciMesh fibonacci_mesh;
        std::vector<size_t> indices;

        ModeID mode_id;
        ModeField mode_field;

        Detector() = default;

        Detector(
            const std::string &mode_number,
            const size_t sampling,
            const double numerical_aperture,
            const double cache_numerical_aperture,
            const double phi_offset,
            const double gamma_offset,
            const double polarization_filter,
            const double rotation,
            const bool coherent,
            const bool mean_coupling,
            const double medium_refractive_index = 1.0)
        :   mode_number(mode_number),
            sampling(sampling),
            numerical_aperture(numerical_aperture),
            cache_numerical_aperture(cache_numerical_aperture),
            phi_offset(phi_offset),
            gamma_offset(gamma_offset),
            polarization_filter(polarization_filter),
            rotation(rotation),
            coherent(coherent),
            mean_coupling(mean_coupling),
            medium_refractive_index(medium_refractive_index)
        {
            this->initialize(medium_refractive_index);
        }

        /**
         * @brief Computes the coupling coefficient for a given scatterer.
         * @param scatterer The scatterer for which the coupling coefficient is computed.
         * @return The coupling coefficient.
         * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
         */
        double get_coupling(const BaseScatterer& scatterer) const;

        /**
         * @brief Computes the structured scalar field for the detector.
         * @param sampling The number of samples to use for the structured scalar field.
         * @return A numpy array containing the structured scalar field.
         * @note This function computes the scalar field based on the mode field and the Fibonacci mesh.
         */
        [[nodiscard]] pybind11::array_t<complex128> get_structured_scalarfield(const size_t sampling) const;

    private:
        /**
         * @brief Computes the coupling coefficient for a given scatterer, considering coherent or non-coherent modes.
         * @param scatterer The scatterer for which the coupling coefficient is computed.
         * @return The coupling coefficient.
         * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
         */
        double get_coupling_point_no_coherent(const BaseScatterer& scatterer) const;

        /**
         * @brief Computes the coupling coefficient for a given scatterer, considering coherent modes.
         * @param scatterer The scatterer for which the coupling coefficient is computed.
         * @return The coupling coefficient.
         * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
         */
        double get_coupling_mean_no_coherent(const BaseScatterer& scatterer) const;

        /**
         * @brief Computes the coupling coefficient for a given scatterer, considering coherent modes.
         * @param scatterer The scatterer for which the coupling coefficient is computed.
         * @return The coupling coefficient.
         * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
         */
        double get_coupling_point_coherent(const BaseScatterer& scatterer) const;

        /**
         * @brief Computes the coupling coefficient for a given scatterer, considering coherent modes.
         * @param scatterer The scatterer for which the coupling coefficient is computed.
         * @return The coupling coefficient.
         * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
         */
        double get_coupling_mean_coherent(const BaseScatterer& scatterer) const;

            /**
         * @brief Converts numerical aperture to angle in radians.
         * @param NA The numerical aperture to convert.
         * @return The angle in radians corresponding to the numerical aperture.
         * @note This function uses the formula: angle = asin(NA).
         */
        double numercical_aperture_to_angle(double numerical_aperture) const;

        /**
         * @brief Initializes the detector with the given medium refractive index.
         * @param medium_refractive_index The refractive index of the medium in which the detector operates.
         * @note This function sets up the mode field and Fibonacci mesh based on the provided parameters.
         */
        void initialize(const double &medium_refractive_index);

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
        get_projected_fields(const std::vector<complex128>& theta_field, const std::vector<complex128>& phi_field) const;


        /**
         * @brief Applies a scalar field to the horizontal and vertical projections of the coupling coefficients.
         * @param field0 The horizontal projection field.
         * @param field1 The vertical projection field.
         * @note This function modifies the input fields in place, applying the scalar field to both projections.
         */
        void apply_scalar_field(std::vector<complex128> &field0, std::vector<complex128> &field1) const;

        /**
         * @brief Computes the squared norm of a 1-norm vector.
         * @tparam T Type of the elements in the vector.
         * @param array The input vector.
         * @return The squared norm of the vector.
         * @note This function computes the squared 1-norm, which is the square of the absolute value of the sum of the elements.
         */
        template <typename T> inline double get_norm1_squared(const std::vector<T>& array) const;

        /**
         * @brief Computes the squared norm of a 2-norm vector.
         * @tparam T Type of the elements in the vector.
         * @param array The input vector.
         * @return The squared norm of the vector.
         * @note This function computes the squared 2-norm, which is the sum of the squares of the absolute values of the elements.
         */
        template <typename T> inline double get_norm2_squared(const std::vector<T>& array) const;

        /**
         * @brief Squares each element in the input vector.
         * @tparam T Type of the elements in the vector.
         * @param array The input vector to be squared.
         * @note This function modifies the input vector in place, squaring each element.
         */
        template <typename T> inline void square_array(std::vector<T>& array);

        /**
         * @brief Applies a polarization filter to the coupling coefficients.
         * @tparam T Type of the coupling coefficients.
         * @param coupling_theta The coupling coefficient for the theta direction.
         * @param coupling_phi The coupling coefficient for the phi direction.
         * @param polarization_filter The polarization filter value.
         */
        template <typename T> inline void apply_polarization_filter(T& coupling_theta, T& coupling_phi, double polarization_filter) const;
};



