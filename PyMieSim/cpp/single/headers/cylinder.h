#pragma once

#include "single/includes/base_class.cpp"

namespace CYLINDER
{

    class Scatterer: public BaseScatterer
    {
        public:
            double diameter = 0.0;
            complex128 index = {1.0, 0.0};
            std::vector<size_t> indices;

            std::vector<complex128> a1n;
            std::vector<complex128> b1n;
            std::vector<complex128> a2n;
            std::vector<complex128> b2n;


            pybind11::array_t<complex128> get_a1n_py(size_t _max_order) { _max_order == 0 ? _max_order = this->max_order : _max_order = max_order; return vector_to_numpy(a1n, {_max_order}); }
            pybind11::array_t<complex128> get_b1n_py(size_t _max_order) { _max_order == 0 ? _max_order = this->max_order : _max_order = max_order; return vector_to_numpy(b1n, {_max_order}); }
            pybind11::array_t<complex128> get_a2n_py(size_t _max_order) { _max_order == 0 ? _max_order = this->max_order : _max_order = max_order; return vector_to_numpy(a2n, {_max_order}); }
            pybind11::array_t<complex128> get_b2n_py(size_t _max_order) { _max_order == 0 ? _max_order = this->max_order : _max_order = max_order; return vector_to_numpy(b2n, {_max_order}); }

            std::vector<complex128> get_a1n() const { return a1n; };
            std::vector<complex128> get_b1n() const { return b1n; };
            std::vector<complex128> get_a2n() const { return a2n; };
            std::vector<complex128> get_b2n() const { return b2n; };

            complex128 get_a11() const { return this->a1n[0]; }
            complex128 get_a12() const { return this->a1n[1]; }
            complex128 get_a13() const { return this->a1n[2]; }
            complex128 get_a21() const { return this->a2n[0]; }
            complex128 get_a22() const { return this->a2n[1]; }
            complex128 get_a23() const { return this->a2n[2]; }
            complex128 get_b11() const { return this->b1n[0]; }
            complex128 get_b12() const { return this->b1n[1]; }
            complex128 get_b13() const { return this->b1n[2]; }
            complex128 get_b21() const { return this->b2n[0]; }
            complex128 get_b22() const { return this->b2n[1]; }
            complex128 get_b23() const { return this->b2n[2]; }

            Scatterer(double diameter, complex128 index, double medium_index, const SOURCE::BaseSource &source, size_t max_order = 0) :
            BaseScatterer(max_order, source, medium_index), diameter(diameter), index(index)
            {
                this->compute_area();
                this->compute_size_parameter();
                this->max_order = (max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : max_order;
                this->compute_an_bn();
            }

            void compute_size_parameter() override {
                this->size_parameter = PI * this->diameter / source.wavelength;
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

    };
}

