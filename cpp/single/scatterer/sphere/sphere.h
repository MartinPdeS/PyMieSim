#pragma once

#include <single/scatterer/base_scatterer/base_scatterer.h>
#include <utils/constants.h>
#include <material/material.h>

using complex128 = std::complex<double>;


class Sphere: public BaseScatterer
{
    public:
        double diameter;
        std::shared_ptr<BaseMaterial> material;
        inline static const std::vector<std::string> property_names = {
            "size_parameter",
            "cross_section",
            "g",
            "g_with_farfields",
            "Qsca",
            "Qext",
            "Qabs",
            "Qback",
            "Qratio",
            "Qforward",
            "Qpr",
            "Csca",
            "Cext",
            "Cabs",
            "Cback",
            "Cratio",
            "Cforward",
            "Cpr",
        };

        /**
         * @brief Constructs a Sphere scatterer with the given parameters.
         * @param diameter The diameter of the sphere.
         * @param material A shared pointer to the material of the sphere.
         * @param medium A shared pointer to the surrounding medium.
         * @param max_order The maximum order of the coefficients to compute.
         */
        Sphere(
            const double diameter,
            std::shared_ptr<BaseMaterial> material,
            std::shared_ptr<BaseMedium> medium,
            const size_t max_order = 0
        )
        :   BaseScatterer(max_order, std::move(medium)),
            diameter(diameter),
            material(std::move(material))
        {}

        /**
         * @brief Constructs a Sphere scatterer with the given parameters, using a constant medium.
         * @param diameter The diameter of the sphere.
         * @param material A shared pointer to the material of the sphere.
         * @param medium The refractive index of the surrounding medium.
         * @param max_order The maximum order of the coefficients to compute.
         */
        Sphere(
            const double diameter,
            std::shared_ptr<BaseMaterial> material,
            const double medium,
            const size_t max_order = 0
        )
        : Sphere(
            diameter,
            std::move(material),
            std::make_shared<ConstantMedium>(medium),
            max_order
        )
        {}

        /**
         * @brief Constructs a Sphere scatterer with the given parameters, using a constant material.
         * @param diameter The diameter of the sphere.
         * @param material The refractive index of the sphere material.
         * @param medium A shared pointer to the surrounding medium.
         * @param max_order The maximum order of the coefficients to compute.
         */
        Sphere(
            const double diameter,
            const complex128 material,
            std::shared_ptr<BaseMedium> medium,
            const size_t max_order = 0
        )
        : Sphere(
            diameter,
            std::make_shared<ConstantMaterial>(material),
            std::move(medium),
            max_order
        )
        {}

        /**
         * @brief Constructs a Sphere scatterer with the given parameters, using constant material and medium.
         * @param diameter The diameter of the sphere.
         * @param material The refractive index of the sphere material.
         * @param medium The refractive index of the surrounding medium.
         * @param max_order The maximum order of the coefficients to compute.
         */
        Sphere(
            const double diameter,
            const complex128 material,
            const double medium,
            const size_t max_order = 0
        )
        : Sphere(
            diameter,
            std::make_shared<ConstantMaterial>(material),
            std::make_shared<ConstantMedium>(medium),
            max_order
        )
        {}

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
        std::pair<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const override;

        /**
         * @brief Computes the near-field electromagnetic fields for a sphere.
         * @param x A vector of x-coordinates where the fields are computed.
         * @param y A vector of y-coordinates where the fields are computed.
         * @param z A vector of z-coordinates where the fields are computed.
         * @param field_type The type of field to compute (e.g., "E", "H").
         * @param radius The radius of the sphere.
         * @return A vector of complex128 values representing the near-field electromagnetic fields.
         */
        std::vector<complex128> get_total_nearfields(
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<double>& z,
            const std::string& field_type,
            const std::shared_ptr<BaseSource>& source
        ) override;

        std::vector<complex128> get_scattered_nearfields(
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<double>& z,
            const std::string& field_type,
            const std::shared_ptr<BaseSource>& source
        ) override;

    public:
        void print_properties(int precision) const override
        {
            const double radius_m = 0.5 * this->diameter;
            const double volume_m3 = (4.0 / 3.0) * Constants::PI * std::pow(radius_m, 3);
            const double cross_section_m2 = Constants::PI * std::pow(radius_m, 2);

            // Header
            std::printf("\n");
            std::printf("Scatterer Type: Sphere\n");
            std::printf("%-16s | %.*e | %-6s\n", "property", precision, 0.0, "unit");
            std::printf("-----------------|--------------------|------\n");

            // Geometry
            std::printf("%-16s | %.*e | %-6s\n", "size_parameter", precision, this->size_parameter, "");
            std::printf("%-16s | %.*e | %-6s\n", "radius",         precision, radius_m,            "m");
            std::printf("%-16s | %.*e | %-6s\n", "RI real",        precision, this->material->get_refractive_index().real(), "");
            std::printf("%-16s | %.*e | %-6s\n", "RI imag",        precision, this->material->get_refractive_index().imag(), "");
            std::printf("%-16s | %.*e | %-6s\n", "medium RI",      precision, this->medium->get_refractive_index(), "");

            std::printf("%-16s | %.*e | %-6s\n", "volume",         precision, volume_m3,           "m^3");
            std::printf("%-16s | %.*e | %-6s\n", "cross_section",  precision, cross_section_m2,    "m^2");
            std::printf("%-16s | %.*e | %-6s\n", "g",              precision, this->get_g(),          "");

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
