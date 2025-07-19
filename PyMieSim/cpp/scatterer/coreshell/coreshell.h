#pragma once

#include <complex>
#include "scatterer/base_scatterer/base_scatterer.h"
#include "utils/special_function.cpp"

using complex128 = std::complex<double>;


class CoreShell: public BaseScatterer
{
    public:
        double core_diameter;
        double shell_thickness;
        double shell_diameter;
        complex128 core_refractive_index;
        complex128 shell_refractive_index;
        double x_core;
        double x_shell;

        /**
         * @brief Constructs a CoreShell object.
         * @param core_diameter The diameter of the core.
         * @param shell_thickness The thickness of the shell.
         * @param core_refractive_index The refractive index of the core.
         * @param shell_refractive_index The refractive index of the shell.
         * @param medium_refractive_index The refractive index of the medium.
         * @param source The light source.
         * @param max_order The maximum order of the scattering coefficients (default is 0, which means it will be computed).
         * @param compute_c_d Whether to compute the cn and dn coefficients (default is false).
         */
        CoreShell(double core_diameter, double shell_thickness, complex128 core_refractive_index, complex128 shell_refractive_index, double medium_refractive_index, const BaseSource &source, size_t max_order = 0, bool compute_c_d = false);

        /**
         * @brief Computes the size parameter for the sphere.
         * The size parameter is defined as (wavenumber * diameter / 2) * medium_refractive_index.
         */
        void compute_size_parameter() override;

        /**
         * @brief Computes the cross-sectional area of the sphere.
         * This is calculated as π * (diameter / 2)².
         */
        void compute_cross_section() override;

        /**
         * @brief Computes the coefficients an and bn for a sphere.
         * @param _max_order The maximum order of the coefficients to compute.
         */
        void compute_an_bn(const size_t max_order);

        /**
         * @brief Computes the coefficients cn and dn for a sphere.
         * @param _max_order The maximum order of the coefficients to compute.
         */
        void compute_cn_dn(size_t max_order);

        /**
         * @brief Computes the scattering efficiency Qsca for a sphere.
         * @return The scattering efficiency Qsca.
         */
        double get_Qsca() const override;

        /**
         * @brief Computes the extinction efficiency Qext for a sphere.
         * @return The extinction efficiency Qext.
         */
        double get_Qext() const override;

        /**
         * @brief Computes the backscattering efficiency Qback for a sphere.
         * @return The backscattering efficiency Qback.
         */
        double get_Qback() const override;

        /**
         * @brief Computes the forward scattering efficiency Qforward for a sphere.
         * @return The forward scattering efficiency Qforward.
         */
        double get_Qforward() const override;

        /**
         * @brief Computes the asymmetry factor g for a sphere.
         * @return The asymmetry factor g.
         */
        double get_g() const override;

        /**
         * @brief Computes the scattering amplitudes S1 and S2 for a sphere.
         * @param phi A vector of angles in radians at which to compute the scattering amplitudes.
         * @return A tuple containing two vectors: S1 and S2, which are the scattering amplitudes.
         */
        std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const override;

    private:
        /**
         * @brief Applies the medium refractive index to the core and shell refractive indices and diameters.
         */
        void apply_medium();
};


// -
