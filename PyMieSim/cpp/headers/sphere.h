#pragma once

#include "base_spherical_scatterer.cpp"

namespace SPHERE
{
    using complex128 = std::complex<double>;

    class Scatterer: public BaseSphericalScatterer
    {
        public:
            double diameter = 0.0;
            complex128 index = {1.0, 0.0};

            Scatterer() = default;

            Scatterer(double wavelength, double amplitude, double diameter, complex128 index,
                double medium_index, std::vector<complex128> jones_vector, size_t max_order) :
                BaseSphericalScatterer(wavelength, jones_vector, amplitude, medium_index), diameter(diameter), index(index)
            {
                compute_size_parameter();
                this->max_order = max_order;
                compute_area();
                compute_an_bn();
            }

            Scatterer(double wavelength, double amplitude, double diameter, complex128 index,
                double medium_index, std::vector<complex128> jones_vector) :
                BaseSphericalScatterer(wavelength, jones_vector, amplitude, medium_index), diameter(diameter), index(index)
            {
                compute_size_parameter();
                this->max_order = get_wiscombe_criterion(this->size_parameter);
                compute_area();
                compute_an_bn();
            }

            Scatterer(double diameter, complex128 index, double medium_index, SOURCE::BaseSource &source, size_t max_order) :
                BaseSphericalScatterer(source, medium_index), diameter(diameter), index(index)
            {
                compute_size_parameter();
                this->max_order = max_order;
                compute_area();
                compute_an_bn();
            }

            Scatterer(double diameter, complex128 index, double medium_index, SOURCE::BaseSource &source) :
                BaseSphericalScatterer(source, medium_index), diameter(diameter), index(index)
            {
                compute_size_parameter();
                this->max_order = this->get_wiscombe_criterion(this->size_parameter);
                compute_area();
                compute_an_bn();
            }

            void compute_cn_dn();
            void compute_an_bn();
            void compute_size_parameter();
            void compute_area();

    };

    class Set
    {
        public:
            std::vector<double> diameter;
            std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>> scatterer;
            std::variant<std::vector<double>, std::vector<std::vector<double>>> medium;
            std::vector<size_t> shape;

            Set() = default;

            // Unified Constructor for all types of scatterer and medium
            Set(const std::vector<double>& diameter,
                std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>> scatterer_param,
                std::variant<std::vector<double>, std::vector<std::vector<double>>> medium_param)
                : diameter(diameter), scatterer(scatterer_param), medium(medium_param)
            {
                update_shape();
            }

            void update_shape() {
                shape.clear();
                shape.push_back(diameter.size());

                if (std::holds_alternative<std::vector<std::vector<complex128>>>(scatterer)) {

                    shape.push_back(std::get<std::vector<std::vector<complex128>>>(scatterer).size());

                } else {
                    shape.push_back(std::get<std::vector<complex128>>(scatterer).size());
                }

                if (std::holds_alternative<std::vector<std::vector<double>>>(medium)) {
                    shape.push_back(std::get<std::vector<std::vector<double>>>(medium).size());
                } else {
                    shape.push_back(std::get<std::vector<double>>(medium).size());
                }
            }

        Scatterer to_object(size_t d, size_t i, size_t wl, size_t mi, SOURCE::BaseSource& source) const
        {
            return Scatterer(
                diameter[d],
                std::holds_alternative<std::vector<std::vector<complex128>>>(scatterer) ? std::get<std::vector<std::vector<complex128>>>(scatterer)[i][wl] : std::get<std::vector<complex128>>(scatterer)[i],
                std::holds_alternative<std::vector<std::vector<double>>>(medium) ? std::get<std::vector<std::vector<double>>>(medium)[mi][wl] : std::get<std::vector<double>>(medium)[mi],
                source
            );
        }
    };
}

// -
