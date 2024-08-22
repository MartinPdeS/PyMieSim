#pragma once

#include "base_class.cpp"

namespace CYLINDER
{

    class Scatterer: public Base
    {
        public:
            double diameter = 0.0;
            complex128 index = {1.0, 0.0};

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

            Scatterer(double diameter, complex128 index, double medium_index, SOURCE::BaseSource &source, size_t max_order = 0) :
            Base(max_order, source, medium_index), diameter(diameter), index(index)
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

