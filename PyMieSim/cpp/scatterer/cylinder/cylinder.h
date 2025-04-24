#pragma once

#include "scatterer/base_scatterer/base_scatterer.h"

#define DEFINE_COMPLEX_VECTOR(name) \
    std::vector<complex128> name##n; \
    std::vector<complex128> get_##name##n() const { return name##n; }

#define DEFINE_GETTER(name, index) \
    complex128 get_##name##index() const { return this->name##n[index - 1]; }

#define DEFINE_GETTER_ABS(name, index) \
    double get_##name##index##_abs() const { return abs(this->name##n[index - 1]); }

#define DEFINE_COEFFICIENTS(name) \
    DEFINE_COMPLEX_VECTOR(name) \
    DEFINE_GETTER(name, 1) \
    DEFINE_GETTER(name, 2) \
    DEFINE_GETTER(name, 3) \
    DEFINE_GETTER(name, 4) \
    DEFINE_GETTER_ABS(name, 1) \
    DEFINE_GETTER_ABS(name, 2) \
    DEFINE_GETTER_ABS(name, 3) \
    DEFINE_GETTER_ABS(name, 4) \
    pybind11::array_t<complex128> get_##name##n_list_py(size_t _max_order) { _max_order = (_max_order == 0 ? this->max_order : _max_order); return _vector_to_numpy(name##n, {_max_order}); }


class Cylinder: public BaseScatterer
{
    public:
        double diameter;
        complex128 refractive_index;
        DEFINE_COEFFICIENTS(a1)
        DEFINE_COEFFICIENTS(b1)
        DEFINE_COEFFICIENTS(a2)
        DEFINE_COEFFICIENTS(b2)

        Cylinder(double diameter, complex128 refractive_index, double medium_refractive_index, const BaseSource &source, size_t max_order = 0) :
        BaseScatterer(max_order, source, medium_refractive_index), diameter(diameter), refractive_index(refractive_index)
        {
            this->compute_area();
            this->compute_size_parameter();
            this->max_order = (max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : max_order;
            this->compute_an_bn();
        }

        void compute_size_parameter() override {
            this->size_parameter = PI * this->diameter / source.wavelength;
            this->size_parameter_squared = pow(this->size_parameter, 2);
        }

        void compute_area() override {
            this->area = this->diameter;
        }

        double get_g() const override;
        double get_Qsca() const override;
        double get_Qext() const override;
        double process_polarization(const complex128 value_0, const complex128 value_1) const ;
        void compute_an_bn();

        std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const override;

        std::vector<complex128> compute_dn(double nmx, complex128 z) const { //Page 205 of BH
            std::vector<complex128> Dn(nmx, 0.0);

            for (double n = nmx - 1; n > 0; n--)
                Dn[n-1] = n/z - ( 1. / (Dn[n] + n/z) );

            return Dn;
        }

        double get_Qback() const override {
            throw std::logic_error{"Function not implemented!"};
        }

};
