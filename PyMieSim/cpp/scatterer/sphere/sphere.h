#pragma once

#include <complex>
#include "scatterer/base_spherical_scatterer.cpp"


using complex128 = std::complex<double>;


class Sphere: public BaseSphericalScatterer
{
    public:
        double diameter;
        complex128 refractive_index;

        Sphere(const double diameter, const complex128 refractive_index, const double medium_refractive_index, const BaseSource &source, size_t max_order = 0) :
            BaseSphericalScatterer(source, max_order, medium_refractive_index), diameter(diameter), refractive_index(refractive_index)
        {
            this->compute_area();
            this->compute_size_parameter();
            this->max_order = (max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : max_order;
            this->compute_an_bn();
        }

        void compute_size_parameter() override {
            this->size_parameter = source.wavenumber * this->diameter / 2 * this->medium_refractive_index;
            this->size_parameter_squared = pow(this->size_parameter, 2);
        }

        void compute_area() override {
            this->area = PI * std::pow(this->diameter / 2.0, 2);
        }

        void compute_cn_dn();
        void compute_an_bn();
};

// -
