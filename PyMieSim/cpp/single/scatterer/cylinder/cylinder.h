#pragma once

#include <single/scatterer/base_scatterer/base_scatterer.h>

using complex128 = std::complex<double>;


class InfiniteCylinder: public BaseScatterer
{
    public:
        double diameter;
        complex128 refractive_index;
        std::vector<std::string> property_names = {
            "size_parameter",
            "radius",
            "cross_section",
            "g",
            "Qsca",
            "Qext",
            "Qabs",
            "Csca",
            "Cext",
            "Cabs",
        };

        /**
         * @brief Constructs a InfiniteCylinder object.
         * @param diameter The diameter of the cylinder.
         * @param refractive_index The refractive index of the cylinder.
         * @param medium_refractive_index The refractive index of the medium.
         * @param source The light source.
         * @param max_order The maximum order of the scattering coefficients.
         */
        InfiniteCylinder(double diameter, complex128 refractive_index, double medium_refractive_index, std::shared_ptr<BaseSource> source, size_t max_order = 0);

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
        void compute_an_bn(const size_t max_order = 0) override;

        /**
         * @brief Computes the coefficients cn and dn for a sphere.
         * @param _max_order The maximum order of the coefficients to compute.
         * * @note This function is not implemented in the original code, so it throws an exception.
         */
        void compute_cn_dn(const size_t max_order = 0) override;

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
         * @brief Computes the backscattering efficiency Qback for a cylinder.
         * @return The backscattering efficiency Qback.
         * * @note This function is not implemented in the original code, so it throws an exception.
         */
        double get_Qback() const override {throw std::logic_error{"Function not implemented!"};}

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
         * @brief Processes the polarization of the scattered light.
         * @param value_0 The first polarization component.
         * @param value_1 The second polarization component.
         * @return The processed polarization value.
         */
        double process_polarization(const complex128 value_0, const complex128 value_1) const;

        /**
         * @brief Computes the dn coefficients for the cylinder.
         * This is based on the formula from Bohren and Huffman, page 205.
         * @param nmx The maximum order of the coefficients.
         * @param z The complex argument for the Bessel functions.
         * @return A vector of dn coefficients.
         */
        std::vector<complex128> compute_dn(double nmx, complex128 z) const;  // Page 205 of BH

        /**
         * @brief Computes the near-field electromagnetic fields for a cylinder.
         * @param x A vector of x-coordinates where the fields are computed.
         * @param y A vector of y-coordinates where the fields are computed.
         * @param z A vector of z-coordinates where the fields are computed.
         * @param field_type The type of field to compute (e.g., "E", "H").
         * @return A vector of complex128 values representing the near-field electromagnetic fields.
         */
        std::vector<complex128> compute_total_nearfields(
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<double>& z,
            const std::string& field_type
        ) override {
            throw std::logic_error{"Function not implemented!"};
            return std::vector<complex128>{};
        };

        /**
         * @brief Computes scattered near-field electromagnetic fields using an and bn coefficients.
         * @param x A vector of x-coordinates where the fields are computed.
         * @param y A vector of y-coordinates where the fields are computed.
         * @param z A vector of z-coordinates where the fields are computed.
         * @param field_type The type of field to compute (e.g., "E", "H").
         * @return A vector of complex128 values representing the scattered near-field electromagnetic fields.
         */
        std::vector<complex128> compute_scattered_nearfields(
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<double>& z,
            const std::string& field_type
        ) override {
            throw std::logic_error{"Function not implemented!"};
            return std::vector<complex128>{};
        };

    public:
        void print_properties(int precision) const
        {
            const double radius_m = 0.5 * this->diameter;

            // If you already store cross_section, use it; otherwise compute consistent with your code.
            // Your docstring says compute_cross_section() is pi*(d/2)^2, so:
            const double cross_section_m2 = Constants::PI * radius_m * radius_m;

            std::printf("\n");
            std::printf("%-16s | %.*e | %-6s\n", "property", precision, 0.0, "unit");
            std::printf("-----------------|--------------------|------\n");

            std::printf("%-16s | %.*e | %-6s\n", "size_parameter", precision, this->size_parameter, "");
            std::printf("%-16s | %.*e | %-6s\n", "radius",         precision, radius_m,            "m");
            std::printf("%-16s | %.*e | %-6s\n", "cross_section",  precision, cross_section_m2,    "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "g",              precision, this->get_g(),             "");

            std::printf("%-16s | %.*e | %-6s\n", "Qsca", precision, this->get_Qsca(), "");
            std::printf("%-16s | %.*e | %-6s\n", "Qext", precision, this->get_Qext(), "");
            std::printf("%-16s | %.*e | %-6s\n", "Qabs", precision, this->get_Qabs(), "");

            std::printf("%-16s | %.*e | %-6s\n", "Csca", precision, this->get_Csca(), "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "Cext", precision, this->get_Cext(), "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "Cabs", precision, this->get_Cabs(), "m^2");

            std::printf("\n");
        }

    };
