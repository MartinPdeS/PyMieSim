#pragma once


#include "base_scatterer.cpp"
#include "numpy_interface.cpp"
#include "VSH.cpp"
#include <cmath>
#include <vector>
#include <complex>


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
            double n_medium = 1.0;
            complex128 index = {1.0, 0.0};

            std::vector<complex128> an;
            std::vector<complex128> bn;
            std::vector<complex128> cn;
            std::vector<complex128> dn;

            Scatterer() = default;

            Scatterer(double wavelength, double amplitude, double diameter, complex128 index,
                double n_medium, std::vector<complex128> jones_vector, size_t max_order = 0) :
                BaseScatterer(wavelength, jones_vector, amplitude), diameter(diameter), index(index), n_medium(n_medium)
            {
                initialize(max_order);
            }

            Scatterer(double diameter, complex128 index, double n_medium, SOURCE::BaseSource &source, size_t max_order = 0) :
                BaseScatterer(source), diameter(diameter), index(index), n_medium(n_medium)
            {
                this->initialize(max_order);
            }

            void initialize(size_t max_order) {
                compute_size_parameter();
                compute_max_order(max_order);
                compute_area();
                compute_an_bn();
            }

            pybind11::array_t<complex128> get_an_py() { return vector_to_numpy(an, {max_order}); }
            pybind11::array_t<complex128> get_bn_py() { return vector_to_numpy(bn, {max_order}); }
            pybind11::array_t<complex128> get_cn_py() { return vector_to_numpy(cn, {max_order}); }
            pybind11::array_t<complex128> get_dn_py() { return vector_to_numpy(dn, {max_order}); }
            std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi) const ;

            std::vector<complex128> get_an() const { return an; };
            std::vector<complex128> get_bn() const { return bn; };
            std::vector<complex128> get_cn() const { return cn; };
            std::vector<complex128> get_dn() const { return dn; };


            double get_g() const ;
            double get_Qsca() const ;
            double get_Qext() const ;
            double get_Qback() const ;
            void compute_cn_dn();
            void compute_an_bn();

            void compute_max_order(size_t max_order);
            void compute_size_parameter();
            void compute_area();

  };
}

// -
