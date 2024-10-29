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
            std::vector<size_t> indices;

            Scatterer(const double diameter, const complex128 index, const double medium_index, const SOURCE::BaseSource &source, size_t max_order = 0) :
                BaseSphericalScatterer(source, max_order, medium_index), diameter(diameter), index(index)
            {
                this->compute_area();
                this->compute_size_parameter();
                this->max_order = (max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : max_order;
                this->compute_an_bn();
            }

            void compute_size_parameter() override {
                this->size_parameter = source.wavenumber * this->diameter / 2 * this->medium_index;
                this->size_parameter_squared = pow(this->size_parameter, 2);
            }

            void compute_area() override {
                this->area = PI * std::pow(this->diameter / 2.0, 2);
            }

            void compute_cn_dn();
            void compute_an_bn();
    };


}

// -
