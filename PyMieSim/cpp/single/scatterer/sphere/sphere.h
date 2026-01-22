#pragma once

#include <single/scatterer/base_scatterer/base_scatterer.h>
#include <utils/constants.h>

using complex128 = std::complex<double>;


class Sphere: public BaseScatterer
{
    public:
        double diameter;
        complex128 refractive_index;
        std::vector<std::string> property_names = {
            "size_parameter",
            "radius",
            "volume",
            "cross_section",
            "g",
            "Qsca",
            "Qext",
            "Qabs",
            "Qback",
            "Qratio",
            "Qpr",
            "Csca",
            "Cext",
            "Cabs",
            "Cback",
            "Cratio",
            "Cpr"
        };

        /**
         * @brief Constructs a Sphere object.
         * @param diameter The diameter of the sphere.
         * @param refractive_index The refractive index of the sphere.
         * @param medium_refractive_index The refractive index of the medium.
         * @param source The light source.
         * @param max_order The maximum order of the scattering coefficients (default is 0, which means it will be computed).
         * @param compute_cn_dn Whether to compute the cn and dn coefficients (default is false).
         */
        Sphere(const double diameter, const complex128 refractive_index, const double medium_refractive_index, std::shared_ptr<BaseSource> source, size_t max_order = 0);

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

        /**
         * @brief Computes the near-field electromagnetic fields for a sphere.
         * @param x A vector of x-coordinates where the fields are computed.
         * @param y A vector of y-coordinates where the fields are computed.
         * @param z A vector of z-coordinates where the fields are computed.
         * @param field_type The type of field to compute (e.g., "E", "H").
         * @param radius The radius of the sphere.
         * @return A vector of complex128 values representing the near-field electromagnetic fields.
         */
        std::vector<complex128> compute_total_nearfields(
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<double>& z,
            const std::string& field_type
        ) override;

        std::vector<complex128> compute_scattered_nearfields(
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<double>& z,
            const std::string& field_type
        ) override;

    public:
        void print_properties(int precision) const
        {
            const double radius_m = 0.5 * this->diameter;
            const double volume_m3 = (4.0 / 3.0) * Constants::PI * std::pow(radius_m, 3);
            const double cross_section_m2 = Constants::PI * std::pow(radius_m, 2);

            // Header
            std::printf("\n");
            std::printf("%-16s | %.*e | %-6s\n", "property", precision, 0.0, "unit");
            std::printf("-----------------|--------------------|------\n");

            // Geometry
            std::printf("%-16s | %.*e | %-6s\n", "size_parameter", precision, this->size_parameter, "");
            std::printf("%-16s | %.*e | %-6s\n", "radius",         precision, radius_m,            "m");
            std::printf("%-16s | %.*e | %-6s\n", "volume",         precision, volume_m3,           "m^3");
            std::printf("%-16s | %.*e | %-6s\n", "cross_section",  precision, cross_section_m2,    "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "g",              precision, this->get_g(),             "");

            // Efficiencies (dimensionless)
            std::printf("%-16s | %.*e | %-6s\n", "Qsca",   precision, this->get_Qsca(),   "");
            std::printf("%-16s | %.*e | %-6s\n", "Qext",   precision, this->get_Qext(),   "");
            std::printf("%-16s | %.*e | %-6s\n", "Qabs",   precision, this->get_Qabs(),   "");
            std::printf("%-16s | %.*e | %-6s\n", "Qback",  precision, this->get_Qback(),  "");
            std::printf("%-16s | %.*e | %-6s\n", "Qratio", precision, this->get_Qratio(), "");
            std::printf("%-16s | %.*e | %-6s\n", "Qpr",    precision, this->get_Qpr(),    "");

            // Cross sections (m^2)
            std::printf("%-16s | %.*e | %-6s\n", "Csca",   precision, this->get_Csca(),   "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "Cext",   precision, this->get_Cext(),   "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "Cabs",   precision, this->get_Cabs(),   "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "Cback",  precision, this->get_Cback(),  "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "Cratio", precision, this->get_Cratio(), "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "Cpr",    precision, this->get_Cpr(),    "m^2");

            std::printf("\n");
        }

    };
