#pragma once

#include "base_scatterer.cpp"

namespace CYLINDER
{

    class Set
    {
        public:
            std::vector<double> diameter;
            std::vector<double> n_medium;
            std::vector<complex128> index;
            std::vector<std::vector<complex128>> material;
            bool is_material;
            std::vector<size_t> shape;

            Set() = default;
            Set(
                const std::vector<double> &diameter,
                const std::vector<std::vector<complex128>> &material,
                const std::vector<double> &n_medium
            ) : diameter(diameter), material(material), n_medium(n_medium), is_material(true)
            {
                this->shape = {this->diameter.size(), this->material.size(), this->n_medium.size()};
            }

            Set(
                const std::vector<double> &diameter,
                const std::vector<complex128> &index,
                const std::vector<double> &n_medium
            ) : diameter(diameter), index(index), n_medium(n_medium), is_material(false)
            {
                this->shape = {this->diameter.size(), this->index.size(), this->n_medium.size()};
            }
    };

    class Scatterer: public BaseScatterer
    {
        public:
            double diameter = 0.0;
            complex128 index = {1.0, 0.0};

            std::vector<complex128> a1n;
            std::vector<complex128> b1n;
            std::vector<complex128> a2n;
            std::vector<complex128> b2n;


            pybind11::array_t<complex128> get_a1n_py() { return vector_to_numpy(a1n, {max_order}); }
            pybind11::array_t<complex128> get_b1n_py() { return vector_to_numpy(b1n, {max_order}); }
            pybind11::array_t<complex128> get_a2n_py() { return vector_to_numpy(a2n, {max_order}); }
            pybind11::array_t<complex128> get_b2n_py() { return vector_to_numpy(b2n, {max_order}); }

            std::vector<complex128> get_a1n() const { return a1n; };
            std::vector<complex128> get_b1n() const { return b1n; };
            std::vector<complex128> get_a2n() const { return a2n; };
            std::vector<complex128> get_b2n() const { return b2n; };

            Scatterer(double wavelength, double amplitude, double diameter, complex128 index,
            double n_medium, std::vector<complex128> jones_vector, size_t max_order)
            : BaseScatterer(wavelength, jones_vector, amplitude, n_medium), diameter(diameter), index(index)
            {
                this->compute_size_parameter();
                this->max_order = max_order;
                this->compute_area();
                this->compute_an_bn();
            }

            Scatterer(double wavelength, double amplitude, double diameter, complex128 index,
            double n_medium, std::vector<complex128> jones_vector)
            : BaseScatterer(wavelength, jones_vector, amplitude, n_medium), diameter(diameter), index(index)
            {
                this->compute_size_parameter();
                this->max_order = get_wiscombe_criterion(this->size_parameter);
                this->compute_area();
                this->compute_an_bn();
            }

            Scatterer(double diameter, complex128 index, double n_medium, SOURCE::BaseSource &source, size_t max_order) :
                BaseScatterer(source, n_medium), diameter(diameter), index(index)
            {
                this->compute_size_parameter();
                this->max_order = max_order;
                this->compute_area();
                this->compute_an_bn();
            }

            Scatterer(double diameter, complex128 index, double n_medium, SOURCE::BaseSource &source) :
                BaseScatterer(source, n_medium), diameter(diameter), index(index)
            {
                this->compute_size_parameter();
                this->max_order = get_wiscombe_criterion(this->size_parameter);
                this->compute_area();
                this->compute_area();
                this->compute_an_bn();
            }

            void compute_size_parameter();
            void compute_area();
            double get_g() const ;
            double get_Qsca() const ;
            double get_Qext() const ;
            double process_polarization(complex128 &value_0, complex128& value_1) const ;
            void compute_an_bn();

            std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const;

  };

}

