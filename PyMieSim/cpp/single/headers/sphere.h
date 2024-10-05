#pragma once

#include "single/includes/base_spherical_scatterer.cpp"

using complex128 = std::complex<double>;

namespace SPHERE
{

    class Scatterer: public BaseSphericalScatterer
    {
        public:
            double diameter;
            complex128 index;

            Scatterer(const double diameter, const complex128 index, const double medium_index, const SOURCE::BaseSource &source, size_t max_order = 0) :
                BaseSphericalScatterer(source, max_order, medium_index), diameter(diameter), index(index)
            {
                this->compute_area();
                this->compute_size_parameter();
                this->max_order = (max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : max_order;
                this->compute_an_bn();
            }

            void compute_size_parameter() override {
                this->size_parameter = PI * this->diameter / this->source.wavelength * this->medium_index;
            }

            void compute_area() override {
                this->area = PI * std::pow(this->diameter / 2.0, 2);
            }

            void compute_cn_dn();
            void compute_an_bn();
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

                if (std::holds_alternative<std::vector<std::vector<complex128>>>(scatterer))
                    shape.push_back(std::get<std::vector<std::vector<complex128>>>(scatterer).size());
                else
                    shape.push_back(std::get<std::vector<complex128>>(scatterer).size());

                if (std::holds_alternative<std::vector<std::vector<double>>>(medium))
                    shape.push_back(std::get<std::vector<std::vector<double>>>(medium).size());
                else
                    shape.push_back(std::get<std::vector<double>>(medium).size());
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
