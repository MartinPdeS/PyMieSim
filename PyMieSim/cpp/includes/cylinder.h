#pragma once

#include "base_scatterer.cpp"

namespace CYLINDER
{

    struct State
    {
        double diameter;
        complex128 index;
        double n_medium;

        State() = default;
        State(const double &diameter, const complex128 &index, const double &n_medium)
        : diameter(diameter), index(index), n_medium(n_medium){}

        void apply_medium() {
            this->index /= this->n_medium;
            this->diameter *= this->n_medium;
        }
    };

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

    class Scatterer: public ScatteringProperties
    {
        public:

            std::vector<complex128> a1n, b1n, a2n, b2n;

            State state;

            pybind11::array_t<complex128> get_a1n_py() { return vector_to_numpy(a1n, {max_order}); }
            pybind11::array_t<complex128> get_b1n_py() { return vector_to_numpy(b1n, {max_order}); }
            pybind11::array_t<complex128> get_a2n_py() { return vector_to_numpy(a2n, {max_order}); }
            pybind11::array_t<complex128> get_b2n_py() { return vector_to_numpy(b2n, {max_order}); }

            std::vector<complex128> get_a1n() { return a1n; };
            std::vector<complex128> get_b1n() { return b1n; };
            std::vector<complex128> get_a2n() { return a2n; };
            std::vector<complex128> get_b2n() { return b2n; };

            double get_diameter(){return state.diameter;}
            complex128 get_index(){return state.index;}
            double get_wavelength(){return source.wavelength;}
            double get_k(){return source.k;}
            std::vector<complex128> get_jones_vector(){return source.jones_vector;}

            Scatterer(double wavelength, double amplitude, double diameter, complex128 index,
            double n_medium, std::vector<complex128> jones_vector, size_t max_order = 0)
            : ScatteringProperties(wavelength, jones_vector, amplitude), state(diameter, index, n_medium){
                initialize(max_order);
            }

            Scatterer(State &state, SOURCE::State &source, size_t max_order = 0)
            : ScatteringProperties(source), state(state){
                initialize(max_order);
            }

            void initialize(size_t &max_order);
            void compute_max_order(size_t &max_order);
            void compute_size_parameter();
            void compute_area();
            double get_g();
            double get_Qsca();
            double get_Qext();

            double process_polarization(complex128 &Value0, complex128& value1);

            void compute_an_bn();

            std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi);

  };

}

