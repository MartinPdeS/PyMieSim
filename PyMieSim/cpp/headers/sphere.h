#pragma once

#include "base_spherical_scatterer.cpp"

namespace SPHERE
{
    using complex128 = std::complex<double>;

    class Set
    {
        public:
            std::vector<double> diameter;
            std::vector<complex128> index;
            std::vector<std::vector<complex128>> material;
            std::vector<double> n_medium;
            bool is_material;
            std::vector<size_t> shape;

            Set() = default;
            Set(const std::vector<double> &diameter, const std::vector<std::vector<complex128>> &material, const std::vector<double> &n_medium) :
            diameter(diameter), material(material), n_medium(n_medium), is_material(true)
            {
                this->shape = {this->diameter.size(), this->material.size(), this->n_medium.size()};
            }

            Set(const std::vector<double> &diameter, const std::vector<complex128> &index, const std::vector<double> &n_medium):
            diameter(diameter), index(index), n_medium(n_medium), is_material(false)
            {
                this->shape = {this->diameter.size(), this->index.size(), this->n_medium.size()};
            }
    };

    class Scatterer: public BaseSphericalScatterer
    {
        public:
            double diameter = 0.0;
            complex128 index = {1.0, 0.0};

            Scatterer() = default;

            Scatterer(double wavelength, double amplitude, double diameter, complex128 index,
                double n_medium, std::vector<complex128> jones_vector, size_t max_order) :
                BaseSphericalScatterer(wavelength, jones_vector, amplitude, n_medium), diameter(diameter), index(index)
            {
                compute_size_parameter();
                this->max_order = max_order;
                compute_area();
                compute_an_bn();
            }

            Scatterer(double wavelength, double amplitude, double diameter, complex128 index,
                double n_medium, std::vector<complex128> jones_vector) :
                BaseSphericalScatterer(wavelength, jones_vector, amplitude, n_medium), diameter(diameter), index(index)
            {
                compute_size_parameter();
                this->max_order = get_wiscombe_criterion(this->size_parameter);
                compute_area();
                compute_an_bn();
            }

            Scatterer(double diameter, complex128 index, double n_medium, SOURCE::BaseSource &source, size_t max_order) :
                BaseSphericalScatterer(source, n_medium), diameter(diameter), index(index)
            {
                compute_size_parameter();
                this->max_order = max_order;
                compute_area();
                compute_an_bn();
            }

            Scatterer(double diameter, complex128 index, double n_medium, SOURCE::BaseSource &source) :
                BaseSphericalScatterer(source, n_medium), diameter(diameter), index(index)
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
}

// -
