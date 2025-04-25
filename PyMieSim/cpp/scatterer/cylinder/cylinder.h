#pragma once

#include <complex>
#include "scatterer/base_scatterer/base_scatterer.h"


using complex128 = std::complex<double>;


class Cylinder: public BaseScatterer
{
    public:
        double diameter;
        complex128 refractive_index;

        Cylinder(double diameter, complex128 refractive_index, double medium_refractive_index, const BaseSource &source, size_t max_order = 0) :
        BaseScatterer(max_order, source, medium_refractive_index), diameter(diameter), refractive_index(refractive_index)
        {
            this->compute_area();
            this->compute_size_parameter();
            this->max_order = (max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : max_order;
            this->compute_an_bn(this->max_order);
        }

        void compute_size_parameter() override {
            this->size_parameter = source.wavenumber * this->diameter / 2 * this->medium_refractive_index;
            this->size_parameter_squared = pow(this->size_parameter, 2);
        }

        void compute_area() override {
            this->area = this->diameter;
        }

        void compute_an_bn(const size_t max_order);

        double get_Qsca() const override;
        double get_Qext() const override;
        double get_Qback() const override {throw std::logic_error{"Function not implemented!"};}
        double get_g() const override;

        std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const override;

    private:
        double process_polarization(const complex128 value_0, const complex128 value_1) const;

        std::vector<complex128> compute_dn(double nmx, complex128 z) const;  // Page 205 of BH
};
