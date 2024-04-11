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

    struct State {
        double diameter = 0.0;
        double n_medium = 1.0;
        complex128 index = {1.0, 0.0};

        State() = default;

        State(double diameter, complex128 index, double n_medium)
            : diameter(diameter), index(index / n_medium), n_medium(n_medium) {}

        void apply_medium() {
            this->index /= this->n_medium;
            this->diameter *= this->n_medium;
        }

    };

    class Set
    {
        public:
            std::vector<double> diameter;
            std::vector<complex128> index;
            std::vector<std::vector<complex128>> material;
            std::vector<double> n_medium;
            bool bounded_index;

            std::vector<State> state_list;

            Set() = default;
            Set(
                const std::vector<double> &diameter,
                const std::vector<std::vector<complex128>> &material,
                const std::vector<double> &n_medium) :
                diameter(diameter), material(material), n_medium(n_medium)
            {
                bounded_index = true;
            }

            Set(const std::vector<double> &diameter, const std::vector<complex128> &index, const std::vector<double> &n_medium) :
            diameter(diameter), index(index), n_medium(n_medium)
            {
                bounded_index = false;
            }

            State operator[](const size_t &idx){return this->state_list[idx];}

            std::vector<size_t> get_array_shape() const
            {
                if (this->bounded_index)
                    return {this->diameter.size(), this->material.size(), this->n_medium.size()};

                if (!this->bounded_index)
                    return {this->diameter.size(), this->index.size(), this->n_medium.size()};

            }

            size_t get_array_size() const
            {
                std::vector<size_t> full_shape = this->get_array_shape();
                size_t full_size = 1;

                for (auto e: full_shape)
                    full_size *= e;

                return full_size;
            }
    };

    class Scatterer: public ScatteringProperties
    {

        public:
            std::vector<complex128> an;
            std::vector<complex128> bn;
            std::vector<complex128> cn;
            std::vector<complex128> dn;

            State state;

            Scatterer() = default;

            Scatterer(
                double wavelength, double amplitude, double diameter, complex128 index,
                double n_medium, std::vector<complex128> jones_vector, size_t max_order = 0) :
                ScatteringProperties(wavelength, jones_vector, amplitude), state(diameter, index, n_medium)
            {
                initialize(max_order);
            }

            Scatterer(State &state, SOURCE::State &source, size_t max_order = 0) :
                ScatteringProperties(source), state(state)
            {
                this->initialize(max_order);
            }

            void initialize(size_t max_order) {
                compute_size_parameter();
                compute_max_order(max_order);
                compute_area();
                compute_an_bn();
            }


            Cndarray get_an_py() { return vector_to_ndarray(an, {max_order}); }
            Cndarray get_bn_py() { return vector_to_ndarray(bn, {max_order}); }
            Cndarray get_cn_py() { return vector_to_ndarray(cn, {max_order}); }
            Cndarray get_dn_py() { return vector_to_ndarray(dn, {max_order}); }
            std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi);

            std::vector<complex128> get_an() { return an; };
            std::vector<complex128> get_bn() { return bn; };
            std::vector<complex128> get_cn() { return cn; };
            std::vector<complex128> get_dn() { return dn; };


            double get_g();
            double get_Qsca();
            double get_Qext();
            double get_Qback();
            void compute_cn_dn();
            void compute_an_bn();

            void compute_max_order(size_t max_order);
            void compute_size_parameter();
            void compute_area();

  };
}

// -
