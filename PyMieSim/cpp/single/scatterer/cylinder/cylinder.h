#pragma once

#include <single/scatterer/base_scatterer/base_scatterer.h>
#include <material/material.h>

using complex128 = std::complex<double>;


class InfiniteCylinder: public BaseScatterer
{
    public:
        double diameter;
        std::shared_ptr<BaseMaterial> material;
        std::vector<complex128> jones_vector;
        inline static const std::vector<std::string> property_names = {
            "size_parameter",
            "cross_section",
            "g_with_farfields",
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
         * @param medium The medium in which the cylinder is embedded.
         * @param max_order The maximum order of the scattering coefficients.
         */
        InfiniteCylinder(
            double _diameter,
            std::shared_ptr<BaseMaterial> _material,
            std::shared_ptr<BaseMedium> _medium,
            size_t _max_order = 0
        ) : BaseScatterer(_max_order, std::move(_medium)),
            diameter(_diameter),
            material(_material)
        {
        }

        /**
         * @brief Constructs a InfiniteCylinder object with constant material.
         * @param diameter The diameter of the cylinder.
         * @param material The refractive index of the cylinder material.
         * @param medium A shared pointer to the surrounding medium.
         * @param max_order The maximum order of the coefficients to compute.
         */
        InfiniteCylinder(
            double _diameter,
            complex128 _material,
            std::shared_ptr<BaseMedium> _medium,
            size_t _max_order = 0) :
            InfiniteCylinder(
                _diameter,
                std::make_shared<ConstantMaterial>(_material),
                std::move(_medium),
                _max_order
            )
        {}

        /**
         * @brief Constructs a InfiniteCylinder object with constant material and medium.
         * @param diameter The diameter of the cylinder.
         * @param material The refractive index of the cylinder material.
         * @param medium The refractive index of the surrounding medium.
         * @param max_order The maximum order of the scattering coefficients.
         */
        InfiniteCylinder(
            double _diameter,
            complex128 _material,
            double _medium,
            size_t _max_order = 0) :
            InfiniteCylinder(
                _diameter,
                std::make_shared<ConstantMaterial>(_material),
                std::make_shared<ConstantMedium>(_medium),
                _max_order
            )
        {}

        /**
         * @brief Constructs a InfiniteCylinder object with constant medium.
         * @param diameter The diameter of the cylinder.
         * @param material A shared pointer to the material of the cylinder.
         * @param medium The refractive index of the surrounding medium.
         * @param max_order The maximum order of the coefficients to compute.
         */
        InfiniteCylinder(
            double _diameter,
            std::shared_ptr<BaseMaterial> _material,
            double _medium,
            size_t _max_order = 0) :
            InfiniteCylinder(
                _diameter,
                std::move(_material),
                std::make_shared<ConstantMedium>(_medium),
                _max_order
            )
        {}

        /**
         * @brief Constructs a InfiniteCylinder object with constant material and medium.
         * @param diameter The diameter of the cylinder.
         * @param material The refractive index of the cylinder material.
         * @param medium The refractive index of the surrounding medium.
         * @param max_order The maximum order of the scattering coefficients.
         */
        void init(const std::shared_ptr<BaseSource>& source, size_t _max_order = 0) override;

        /**
         * @brief Computes the size parameter for the sphere.
         * The size parameter is defined as (wavenumber * diameter / 2) * medium_refractive_index.
         */
        void compute_size_parameter(const std::shared_ptr<BaseSource>& source) override;

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
        double get_g() const override {throw std::logic_error{"Function not implemented!"};}

        /**
         * @brief Computes the scattering amplitudes S1 and S2 for a sphere.
         * @param phi A vector of angles in radians at which to compute the scattering amplitudes.
         * @return A tuple containing two vectors: S1 and S2, which are the scattering amplitudes.
         */
        std::pair<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const override;

    private:
        /**
         * @brief Processes the polarization of the scattered light.
         * @param value_0 The first polarization component.
         * @param value_1 The second polarization component.
         * @return The processed polarization value.
         */
        double process_polarization(const double value_0, const double value_1) const;

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
        std::vector<complex128> get_total_nearfields(
            const std::vector<double>&,
            const std::vector<double>&,
            const std::vector<double>&,
            const std::string&,
            const std::shared_ptr<BaseSource>&
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
        std::vector<complex128> get_scattered_nearfields(
            const std::vector<double>&,
            const std::vector<double>&,
            const std::vector<double>&,
            const std::string&,
            const std::shared_ptr<BaseSource>&
        ) override {
            throw std::logic_error{"Function not implemented!"};
            return std::vector<complex128>{};
        };

    public:
        void print_properties(int precision) const override
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

            std::printf("%-16s | %.*e | %-6s\n", "Qsca", precision, this->get_Qsca(), "");
            std::printf("%-16s | %.*e | %-6s\n", "Qext", precision, this->get_Qext(), "");
            std::printf("%-16s | %.*e | %-6s\n", "Qabs", precision, this->get_Qabs(), "");

            std::printf("%-16s | %.*e | %-6s\n", "Csca", precision, this->get_Csca(), "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "Cext", precision, this->get_Cext(), "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "Cabs", precision, this->get_Cabs(), "m^2");

            std::printf("\n");
        }

    };
