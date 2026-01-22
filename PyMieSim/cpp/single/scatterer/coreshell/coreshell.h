#pragma once


#include <single/scatterer/base_scatterer/base_scatterer.h>

using complex128 = std::complex<double>;


class CoreShell: public BaseScatterer
{
    public:
        double core_diameter;
        double shell_thickness;
        double total_diameter;
        complex128 core_refractive_index;
        complex128 shell_refractive_index;
        double x_core;
        double x_shell;

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
        CoreShell(double core_diameter, double shell_thickness, complex128 core_refractive_index, complex128 shell_refractive_index, double medium_refractive_index, std::shared_ptr<BaseSource> source, size_t max_order = 0);

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
         * @return A vector of complex128 values representing the near-field electromagnetic fields.
         * @note This method uses the multipole expansion with vector spherical harmonics to compute the
         * near-field electromagnetic fields inside and near the sphere.
         * The internal fields (r < radius) are computed using cn and dn coefficients, while
         * external fields (r > radius) use an and bn coefficients.
         */
        std::vector<complex128> compute_nearfields(
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
        ) override {
            throw std::logic_error{"Function not implemented!"};
            return std::vector<complex128>{};
        };

    private:
        /**
         * @brief Applies the medium refractive index to the core and shell refractive indices and diameters.
         */
        void apply_medium();

    public:
        void print_properties(int precision) const
        {
            // Table header
            std::printf("\n");
            std::printf("+-------------------+------------------------------+-----------+\n");
            std::printf("| name              | value                        | unit      |\n");
            std::printf("+-------------------+------------------------------+-----------+\n");

            // Geometry and indices (these are CoreShell members)
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "core_diameter",    precision, this->core_diameter,    "m");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "shell_thickness",  precision, this->shell_thickness,  "m");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "total_diameter",   precision, this->total_diameter,   "m");

            std::printf("| %-17s | %14.*e + %14.*ei | %-9s |\n",
                "core_n",
                precision, std::real(this->core_refractive_index),
                precision, std::imag(this->core_refractive_index),
                "RIU"
            );

            std::printf("| %-17s | %14.*e + %14.*ei | %-9s |\n",
                "shell_n",
                precision, std::real(this->shell_refractive_index),
                precision, std::imag(this->shell_refractive_index),
                "RIU"
            );

            double radius = 0.5 * this->total_diameter;
            double volume = (4.0 / 3.0) * Constants::PI * std::pow(radius, 3);
            double cross_section = Constants::PI * std::pow(radius, 2);

            // BaseScatterer computed properties (assumed to exist as members)
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "size_parameter",   precision, this->size_parameter,         "1");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "radius",           precision, radius,                       "m");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "volume",           precision, volume,                     "m^3");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "cross_section",    precision, cross_section,              "m^2");

            std::printf("| %-17s | %14.*e                   | %-9s |\n", "g",                precision, this->get_g(),                "1");

            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Qsca",             precision, this->get_Qsca(),             "1");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Qext",             precision, this->get_Qext(),             "1");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Qabs",             precision, this->get_Qabs(),             "1");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Qback",            precision, this->get_Qback(),            "1");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Qratio",           precision, this->get_Qratio(),           "1");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Qpr",              precision, this->get_Qpr(),              "1");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Csca",             precision, this->get_Csca(),             "m^2");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Cext",             precision, this->get_Cext(),             "m^2");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Cabs",             precision, this->get_Cabs(),             "m^2");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Cback",            precision, this->get_Cback(),            "m^2");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Cratio",           precision, this->get_Cratio(),           "1");
            std::printf("| %-17s | %14.*e                   | %-9s |\n", "Cpr",              precision, this->get_Cpr(),              "m^2");

            std::printf("+-------------------+------------------------------+-----------+\n");
            std::printf("\n");
        }

};


// -
